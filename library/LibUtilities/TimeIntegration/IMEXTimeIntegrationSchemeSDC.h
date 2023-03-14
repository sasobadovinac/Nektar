///////////////////////////////////////////////////////////////////////////////
//
// File: IMEXTimeIntegrationSchemeSDC.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Header file of time integration scheme SDC base class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_IMEX_TIME_INTEGRATION_SCHEME_SDC
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_IMEX_TIME_INTEGRATION_SCHEME_SDC

#define LUE LIB_UTILITIES_EXPORT

#include <string>

#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeSDC.h>

namespace Nektar
{
namespace LibUtilities
{
class IMEXTimeIntegrationSchemeSDC : public TimeIntegrationSchemeSDC
{
public:
    IMEXTimeIntegrationSchemeSDC(std::string variant, unsigned int order,
                                 std::vector<NekDouble> freeParams)
        : TimeIntegrationSchemeSDC(variant, order, freeParams)
    {
        ASSERTL0(variant == "GaussRadauLegendre",
                 "only GaussRadauLegendre variant "
                 "(quadrature) type available for IMEXSDC");
        ASSERTL0(0.0 < freeParams[0],
                 "Spectral Deferred Correction Time integration "
                 "scheme bad parameter numbers (> 0.0): " +
                     std::to_string(freeParams[0]));

        std::cerr
            << "WARNING: IMEX Spectral Deferred Correction method has been "
               "implemented but its use is not recommended as the "
               "approach is affected by order-reduction problems."
            << std::endl;

        m_name       = "IMEXSpectralDeferredCorrection";
        m_schemeType = eIMEX;
    }

    static TimeIntegrationSchemeSharedPtr create(
        std::string variant, unsigned int order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<IMEXTimeIntegrationSchemeSDC>::AllocateSharedPtr(
                variant, order, freeParams);

        return p;
    }

    static std::string className;

    TripleArray m_Fexp; /// Array corresponding to the stage Derivatives
    TripleArray m_Fimp; /// Array corresponding to the stage Derivatives
    DoubleArray m_Fexpn;
    DoubleArray m_Fimpn;

protected:
    LUE virtual void v_InitializeScheme(
        const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
        const TimeIntegrationSchemeOperators &op) override;

    LUE virtual ConstDoubleArray &v_TimeIntegrate(
        const int timestep, const NekDouble delta_t,
        const TimeIntegrationSchemeOperators &op) override;

}; // end class IMEXTimeIntegrationSchemeSDC

void IMEXTimeIntegrationSchemeSDC::v_InitializeScheme(
    const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
    const TimeIntegrationSchemeOperators &op)
{
    if (m_initialized)
    {
        return;
    }

    TimeIntegrationSchemeSDC::v_InitializeScheme(deltaT, y_0, time, op);

    // Storage of previous states and associated timesteps.
    m_Fexp = TripleArray(m_nQuadPts);
    m_Fimp = TripleArray(m_nQuadPts);
    for (unsigned int m = 0; m < m_nQuadPts; ++m)
    {
        m_Fexp[m] = DoubleArray(m_nvars);
        m_Fimp[m] = DoubleArray(m_nvars);

        for (unsigned int i = 0; i < m_nvars; ++i)
        {
            m_Fexp[m][i] = SingleArray(m_npoints, 0.0);
            m_Fimp[m][i] = SingleArray(m_npoints, 0.0);
        }
    }

    m_Fexpn = DoubleArray(m_nvars);
    m_Fimpn = DoubleArray(m_nvars);
    for (unsigned int i = 0; i < m_nvars; ++i)
    {
        m_Fexpn[i] = SingleArray(m_npoints, 0.0);
        m_Fimpn[i] = SingleArray(m_npoints, 0.0);
    }
}

ConstDoubleArray &IMEXTimeIntegrationSchemeSDC::v_TimeIntegrate(
    const int timestep, const NekDouble delta_t,
    const TimeIntegrationSchemeOperators &op)
{

    boost::ignore_unused(timestep);

    // 1. Compute initial guess
    for (unsigned int n = 0; n < m_nQuadPts; ++n)
    {
        // 1.1. Compute explicit residual
        op.DoOdeRhs(m_Y[n], m_Fexp[n], m_time + delta_t * m_points[n]);

        // 1.2. Use first-order IMEX as a fist guess
        if (n < m_nQuadPts - 1)
        {
            NekDouble dtn = delta_t * (m_points[n + 1] - m_points[n]);

            // 1.2.1. Update explicit solution
            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                Vmath::Svtvp(m_npoints, dtn, m_Fexp[n][i], 1, m_Y[n][i], 1,
                             m_tmp[i], 1);
            }

            // 1.2.2. Update implicit solution
            op.DoImplicitSolve(m_tmp, m_Y[n + 1],
                               m_time + delta_t * m_points[n + 1], dtn);

            // 1.2.3. Compute implicit flux from updated solution
            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                Vmath::Vsub(m_npoints, m_Y[n + 1][i], 1, m_tmp[i], 1,
                            m_Fimp[n + 1][i], 1);
                Vmath::Smul(m_npoints, 1.0 / dtn, m_Fimp[n + 1][i], 1,
                            m_Fimp[n + 1][i], 1);
            }
        }
    }

    // 2. Return initial guess if m_order = 1
    if (m_order == 1)
    {
        // 2.1. Copy final solution
        for (unsigned int i = 0; i < m_nvars; ++i)
        {
            Vmath::Vcopy(m_npoints, m_Y[m_nQuadPts - 1][i], 1, m_Y[0][i], 1);
        }

        // 2.2. Update time step
        m_time += delta_t;

        // 2.3. Return solution
        return m_Y[m_nQuadPts - 1];
    }

    // 3. Apply SDC correction loop
    // 3.0. Compute total flux
    for (unsigned int n = 0; n < m_nQuadPts - 1; ++n)
    {
        for (unsigned int i = 0; i < m_nvars; ++i)
        {
            Vmath::Vadd(m_npoints, m_Fimp[n + 1][i], 1, m_Fexp[n + 1][i], 1,
                        m_F[n + 1][i], 1);
        }
    }

    for (unsigned int k = 1; k < m_order; ++k)
    {
        // 3.1. Update integrated residual
        UpdateIntegratedFlux(delta_t);

        // 3.2. Loop over quadrature points
        for (unsigned int n = 0; n < m_nQuadPts - 1; ++n)
        {
            NekDouble dtn = delta_t * (m_points[n + 1] - m_points[n]);

            // 3.2.1. Update explicit residual if n > 0
            if (n > 0)
            {
                op.DoOdeRhs(m_Y[n], m_Fexpn, m_time + delta_t * m_points[n]);
            }

            // 3.2.2. Update solution
            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                // Add Fint contribution
                Vmath::Vadd(m_npoints, m_Y[n][i], 1, m_Fint[n][i], 1, m_tmp[i],
                            1);

                if (n > 0)
                {
                    // Add explicit contribution
                    Vmath::Svtvp(m_npoints, m_theta * dtn, m_Fexpn[i], 1,
                                 m_tmp[i], 1, m_tmp[i], 1);
                    Vmath::Svtvp(m_npoints, -m_theta * dtn, m_Fexp[n][i], 1,
                                 m_tmp[i], 1, m_tmp[i], 1);

                    // Copy new explicit residual value to old
                    Vmath::Vcopy(m_npoints, m_Fexpn[i], 1, m_Fexp[n][i], 1);
                }

                // Add implicit contribution
                Vmath::Svtvp(m_npoints, -m_theta * dtn, m_Fimp[n + 1][i], 1,
                             m_tmp[i], 1, m_tmp[i], 1);
            }

            // 3.2.3. Solve implicit system
            op.DoImplicitSolve(m_tmp, m_Y[n + 1],
                               m_time + delta_t * m_points[n + 1],
                               m_theta * dtn);

            // 3.2.4. Compute implicit residual from updated solution
            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                Vmath::Vsub(m_npoints, m_Y[n + 1][i], 1, m_tmp[i], 1,
                            m_Fimp[n + 1][i], 1);
                Vmath::Smul(m_npoints, 1.0 / (m_theta * dtn), m_Fimp[n + 1][i],
                            1, m_Fimp[n + 1][i], 1);
            }
        }
        op.DoOdeRhs(m_Y[m_nQuadPts - 1], m_Fexp[m_nQuadPts - 1],
                    m_time + delta_t * m_points[m_nQuadPts - 1]);

        // 3.3. Compute total residual
        for (unsigned int n = 0; n < m_nQuadPts - 1; ++n)
        {
            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                Vmath::Vadd(m_npoints, m_Fimp[n + 1][i], 1, m_Fexp[n + 1][i], 1,
                            m_F[n + 1][i], 1);
            }
        }
    }

    // 4. Get solution
    // 4.1. Copy final solution
    for (unsigned int i = 0; i < m_nvars; ++i)
    {
        Vmath::Vcopy(m_npoints, m_Y[m_nQuadPts - 1][i], 1, m_Y[0][i], 1);
    }

    // 4.2. Update time step
    m_time += delta_t;

    // 4.3. Return solution
    return m_Y[m_nQuadPts - 1];
}

} // end namespace LibUtilities
} // end namespace Nektar

#endif
