///////////////////////////////////////////////////////////////////////////////
//
// File: ImplicitTimeIntegrationSchemeSDC.h
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

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_IMPLICIT_TIME_INTEGRATION_SCHEME_SDC
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_IMPLICIT_TIME_INTEGRATION_SCHEME_SDC

#define LUE LIB_UTILITIES_EXPORT

#include <string>

#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeSDC.h>

namespace Nektar
{
namespace LibUtilities
{
class ImplicitTimeIntegrationSchemeSDC : public TimeIntegrationSchemeSDC
{
public:
    ImplicitTimeIntegrationSchemeSDC(std::string variant, unsigned int order,
                                     std::vector<NekDouble> freeParams)
        : TimeIntegrationSchemeSDC(variant, order, freeParams)
    {
        ASSERTL0(variant == "GaussRadauLegendre",
                 "only GaussRadauLegendre variant "
                 "(quadrature) type available for ImplicitSDC");
        ASSERTL0(0.0 < freeParams[0],
                 "Spectral Deferred Correction Time integration "
                 "scheme bad parameter numbers (> 0.0): " +
                     std::to_string(freeParams[0]));

        m_name       = "ImplicitSpectralDeferredCorrection";
        m_schemeType = eImplicit;
    }

    static TimeIntegrationSchemeSharedPtr create(
        std::string variant, unsigned int order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<ImplicitTimeIntegrationSchemeSDC>::AllocateSharedPtr(
                variant, order, freeParams);

        return p;
    }

    static std::string className;

protected:
    LUE virtual void v_InitializeScheme(
        const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
        const TimeIntegrationSchemeOperators &op) override;

    LUE virtual ConstDoubleArray &v_TimeIntegrate(
        const int timestep, const NekDouble delta_t,
        const TimeIntegrationSchemeOperators &op) override;

}; // end class ImplicitTimeIntegrationSchemeSDC

void ImplicitTimeIntegrationSchemeSDC::v_InitializeScheme(
    const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
    const TimeIntegrationSchemeOperators &op)
{
    TimeIntegrationSchemeSDC::v_InitializeScheme(deltaT, y_0, time, op);
}

ConstDoubleArray &ImplicitTimeIntegrationSchemeSDC::v_TimeIntegrate(
    const int timestep, const NekDouble delta_t,
    const TimeIntegrationSchemeOperators &op)
{

    boost::ignore_unused(timestep);

    // 1. Compute initial guess
    for (unsigned int n = 0; n < m_nQuadPts - 1; ++n)
    {
        // 1.1. Use implicit Euler as a first guess
        NekDouble dtn = delta_t * (m_points[n + 1] - m_points[n]);
        op.DoImplicitSolve(m_Y[n], m_Y[n + 1],
                           m_time + delta_t * m_points[n + 1], dtn);

        // 1.2. Compute residual from updated solution
        for (unsigned int i = 0; i < m_nvars; ++i)
        {
            Vmath::Vsub(m_npoints, m_Y[n + 1][i], 1, m_Y[n][i], 1,
                        m_F[n + 1][i], 1);
            Vmath::Smul(m_npoints, 1.0 / dtn, m_F[n + 1][i], 1, m_F[n + 1][i],
                        1);
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
    for (unsigned int k = 1; k < m_order; ++k)
    {
        // 3.1. Update integrated residual
        UpdateIntegratedFlux(delta_t);

        // 3.2. Loop over quadrature points
        for (unsigned int n = 0; n < m_nQuadPts - 1; ++n)
        {
            NekDouble dtn = delta_t * (m_points[n + 1] - m_points[n]);

            // Update solution
            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                // Add Fint contribution
                Vmath::Vadd(m_npoints, m_Y[n][i], 1, m_Fint[n][i], 1, m_tmp[i],
                            1);

                // Add implicit contribution
                Vmath::Svtvp(m_npoints, -m_theta * dtn, m_F[n + 1][i], 1,
                             m_tmp[i], 1, m_tmp[i], 1);
            }

            // Solve implicit system
            op.DoImplicitSolve(m_tmp, m_Y[n + 1],
                               m_time + delta_t * m_points[n + 1],
                               m_theta * dtn);

            // Compute residual from updated solution
            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                Vmath::Vsub(m_npoints, m_Y[n + 1][i], 1, m_tmp[i], 1,
                            m_F[n + 1][i], 1);
                Vmath::Smul(m_npoints, 1.0 / (m_theta * dtn), m_F[n + 1][i], 1,
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
