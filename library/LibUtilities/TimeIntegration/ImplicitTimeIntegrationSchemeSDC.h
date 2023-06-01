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

#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeSDC.h>

namespace Nektar
{
namespace LibUtilities
{
class ImplicitTimeIntegrationSchemeSDC : public TimeIntegrationSchemeSDC
{
public:
    ImplicitTimeIntegrationSchemeSDC(std::string variant, size_t order,
                                     std::vector<NekDouble> freeParams)
        : TimeIntegrationSchemeSDC(variant, order, freeParams)
    {
        ASSERTL0(!m_first_quadrature,
                 "Quadrature type that include the left end point (e.g. "
                 "GaussLobattoLegendre) should not be used for ImplicitSDC");

        m_name       = "ImplicitSpectralDeferredCorrection";
        m_schemeType = eImplicit;
    }

    static TimeIntegrationSchemeSharedPtr create(
        std::string variant, size_t order, std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<ImplicitTimeIntegrationSchemeSDC>::AllocateSharedPtr(
                variant, order, freeParams);

        return p;
    }

    static std::string className;

    DoubleArray m_tmp; /// Array for temporary storage

protected:
    LUE virtual void v_InitializeScheme(
        const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
        const TimeIntegrationSchemeOperators &op) override;

    LUE virtual void v_ResidualEval(
        const NekDouble &delta_t, const size_t n,
        const TimeIntegrationSchemeOperators &op) override;

    LUE virtual void v_ResidualEval(
        const NekDouble &delta_t,
        const TimeIntegrationSchemeOperators &op) override;

    LUE virtual void v_ComputeInitialGuess(
        const NekDouble &delta_t,
        const TimeIntegrationSchemeOperators &op) override;

    LUE virtual void v_SDCIterationLoop(
        const NekDouble &delta_t,
        const TimeIntegrationSchemeOperators &op) override;

}; // end class ImplicitTimeIntegrationSchemeSDC

/**
 * @brief Worker method to initialize the integration scheme.
 */
void ImplicitTimeIntegrationSchemeSDC::v_InitializeScheme(
    const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
    const TimeIntegrationSchemeOperators &op)
{
    if (m_initialized)
    {
        for (size_t i = 0; i < m_nvars; ++i)
        {
            // Store the initial values as the first previous state.
            Vmath::Vcopy(m_npoints, y_0[i], 1, m_Y[0][i], 1);
        }
    }
    else
    {
        TimeIntegrationSchemeSDC::v_InitializeScheme(deltaT, y_0, time, op);

        m_tmp = DoubleArray(m_nvars);
        for (size_t i = 0; i < m_nvars; ++i)
        {
            m_tmp[i] = SingleArray(m_npoints, 0.0);
        }
    }
}

/**
 * @brief Worker method to compute the residual.
 */
void ImplicitTimeIntegrationSchemeSDC::v_ResidualEval(
    const NekDouble &delta_t, const size_t n,
    const TimeIntegrationSchemeOperators &op)
{
    if (n == 0)
    {
        // Not implemented, require implicit evaluation for m_F[0].
        // Quadrature type that include the left end point (e.g.
        // GaussLobattoLegendre) should not be used.
    }
    else
    {
        NekDouble dtn = delta_t * (m_tau[n] - m_tau[n - 1]);

        // Update solution
        op.DoImplicitSolve(m_Y[n - 1], m_tmp, m_time + delta_t * m_tau[n],
                           m_theta * dtn);

        // Compute residual from updated solution
        for (size_t i = 0; i < m_nvars; ++i)
        {
            Vmath::Vsub(m_npoints, m_tmp[i], 1, m_Y[n - 1][i], 1, m_F[n][i], 1);
            Vmath::Smul(m_npoints, 1.0 / (m_theta * dtn), m_F[n][i], 1,
                        m_F[n][i], 1);
        }
    }
}

void ImplicitTimeIntegrationSchemeSDC::v_ResidualEval(
    const NekDouble &delta_t, const TimeIntegrationSchemeOperators &op)
{
    for (size_t n = 0; n < m_nQuadPts; ++n)
    {
        v_ResidualEval(delta_t, n, op);
    }
}

/**
 * @brief Worker method to compute the initial SDC guess.
 */
void ImplicitTimeIntegrationSchemeSDC::v_ComputeInitialGuess(
    const NekDouble &delta_t, const TimeIntegrationSchemeOperators &op)
{
    for (size_t n = 0; n < m_nQuadPts; ++n)
    {
        if (n == 0)
        {
            // Not implemented, require implicit evaluation for m_F[0].
            // Quadrature type that include the left end point (e.g.
            // GaussLobattoLegendre) should not be used.
        }
        else
        {
            NekDouble dtn = delta_t * (m_tau[n] - m_tau[n - 1]);

            // Update solution
            op.DoImplicitSolve(m_Y[n - 1], m_Y[n], m_time + delta_t * m_tau[n],
                               dtn);

            // Compute residual from updated solution
            for (size_t i = 0; i < m_nvars; ++i)
            {
                Vmath::Vsub(m_npoints, m_Y[n][i], 1, m_Y[n - 1][i], 1,
                            m_F[n][i], 1);
                Vmath::Smul(m_npoints, 1.0 / dtn, m_F[n][i], 1, m_F[n][i], 1);
            }
        }
    }
}

/**
 * @brief Worker method to compute the SDC iteration.
 */
void ImplicitTimeIntegrationSchemeSDC::v_SDCIterationLoop(
    const NekDouble &delta_t, const TimeIntegrationSchemeOperators &op)
{
    // Update integrated residual
    UpdateIntegratedResidualSFint(delta_t);

    // Add FAS correction to integrated residual
    AddFASCorrectionToSFint();

    // Loop over quadrature points
    for (size_t n = 1; n < m_nQuadPts; ++n)
    {
        NekDouble dtn = delta_t * (m_tau[n] - m_tau[n - 1]);

        // Add rhs terms
        for (size_t i = 0; i < m_nvars; ++i)
        {
            // Add SFint contribution
            Vmath::Vadd(m_npoints, m_Y[n - 1][i], 1, m_SFint[n][i], 1, m_tmp[i],
                        1);

            // Add implicit contribution
            Vmath::Svtvp(m_npoints, -m_theta * dtn, m_F[n][i], 1, m_tmp[i], 1,
                         m_tmp[i], 1);
        }

        // Solve implicit system
        op.DoImplicitSolve(m_tmp, m_Y[n], m_time + delta_t * m_tau[n],
                           m_theta * dtn);

        // Compute residual from updated solution
        for (size_t i = 0; i < m_nvars; ++i)
        {
            Vmath::Vsub(m_npoints, m_Y[n][i], 1, m_tmp[i], 1, m_F[n][i], 1);
            Vmath::Smul(m_npoints, 1.0 / (m_theta * dtn), m_F[n][i], 1,
                        m_F[n][i], 1);
        }
    }
}

} // end namespace LibUtilities
} // end namespace Nektar

#endif
