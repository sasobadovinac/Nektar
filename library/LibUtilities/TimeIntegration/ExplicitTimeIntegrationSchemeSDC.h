///////////////////////////////////////////////////////////////////////////////
//
// File: ExplicitTimeIntegrationSchemeSDC.h
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

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_EXPLICIT_TIME_INTEGRATION_SCHEME_SDC
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_EXPLICIT_TIME_INTEGRATION_SCHEME_SDC

#define LUE LIB_UTILITIES_EXPORT

#include <string>

#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeSDC.h>

namespace Nektar
{
namespace LibUtilities
{
class ExplicitTimeIntegrationSchemeSDC : public TimeIntegrationSchemeSDC
{
public:
    ExplicitTimeIntegrationSchemeSDC(std::string variant, unsigned int order,
                                     std::vector<NekDouble> freeParams)
        : TimeIntegrationSchemeSDC(variant, order, freeParams)
    {
        m_name       = "ExplicitSpectralDeferredCorrection";
        m_schemeType = eExplicit;
    }

    static TimeIntegrationSchemeSharedPtr create(
        std::string variant, unsigned int order,
        std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<ExplicitTimeIntegrationSchemeSDC>::AllocateSharedPtr(
                variant, order, freeParams);

        return p;
    }

    static std::string className;

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

}; // end class ExplicitTimeIntegrationSchemeSDC

/**
 * @brief Worker method to initialize the integration scheme.
 */
void ExplicitTimeIntegrationSchemeSDC::v_InitializeScheme(
    const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
    const TimeIntegrationSchemeOperators &op)
{
    TimeIntegrationSchemeSDC::v_InitializeScheme(deltaT, y_0, time, op);
}

/**
 * @brief Worker method to compute the residual.
 */
void ExplicitTimeIntegrationSchemeSDC::v_ResidualEval(
    const NekDouble &delta_t, const size_t n,
    const TimeIntegrationSchemeOperators &op)
{
    op.DoProjection(m_Y[n], m_Y[n], m_time + delta_t * m_tau[n]);
    op.DoOdeRhs(m_Y[n], m_F[n], m_time + delta_t * m_tau[n]);
}

void ExplicitTimeIntegrationSchemeSDC::v_ResidualEval(
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
void ExplicitTimeIntegrationSchemeSDC::v_ComputeInitialGuess(
    const NekDouble &delta_t, const TimeIntegrationSchemeOperators &op)
{
    for (size_t n = 0; n < m_nQuadPts; ++n)
    {
        // Use explicit Euler as a first guess
        if (n > 0)
        {
            NekDouble dtn = delta_t * (m_tau[n] - m_tau[n - 1]);
            for (size_t i = 0; i < m_nvars; ++i)
            {
                Vmath::Svtvp(m_npoints, dtn, m_F[n - 1][i], 1, m_Y[n - 1][i], 1,
                             m_Y[n][i], 1);
            }
        }

        // Compute residual
        op.DoProjection(m_Y[n], m_Y[n], m_time + delta_t * m_tau[n]);
        op.DoOdeRhs(m_Y[n], m_F[n], m_time + delta_t * m_tau[n]);
    }
}

/**
 * @brief Worker method to compute the SDC iteration.
 */
void ExplicitTimeIntegrationSchemeSDC::v_SDCIterationLoop(
    const NekDouble &delta_t, const TimeIntegrationSchemeOperators &op)
{
    // Update integrated residual
    UpdateIntegratedResidualSFint(delta_t);

    // Add FAS correction to integrated residual
    AddFASCorrectionToSFint();

    // Loop over quadrature points
    for (size_t n = 1; n < m_nQuadPts; ++n)
    {
        // Update solution
        for (size_t i = 0; i < m_nvars; ++i)
        {
            // Add SFint contribution to solution
            Vmath::Vadd(m_npoints, m_Y[n - 1][i], 1, m_SFint[n][i], 1,
                        m_Y[n][i], 1);

            // Add explicit contribution to solution
            if (n > 1)
            {
                NekDouble dtn = delta_t * (m_tau[n] - m_tau[n - 1]);
                Vmath::Svtvp(m_npoints, m_theta * dtn, m_F[n - 1][i], 1,
                             m_Y[n][i], 1, m_Y[n][i], 1);
            }

            // Add explicit contribution to SFint
            if (n < m_nQuadPts - 1)
            {
                NekDouble dtnp = delta_t * (m_tau[n + 1] - m_tau[n]);
                Vmath::Svtvp(m_npoints, -m_theta * dtnp, m_F[n][i], 1,
                             m_SFint[n + 1][i], 1, m_SFint[n + 1][i], 1);
            }
        }

        // Compute residual
        op.DoProjection(m_Y[n], m_Y[n], m_time + delta_t * m_tau[n]);
        op.DoOdeRhs(m_Y[n], m_F[n], m_time + delta_t * m_tau[n]);
    }
}

} // end namespace LibUtilities
} // end namespace Nektar

#endif
