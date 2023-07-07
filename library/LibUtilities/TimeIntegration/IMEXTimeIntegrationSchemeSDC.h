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

#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeSDC.h>

namespace Nektar
{
namespace LibUtilities
{
class IMEXTimeIntegrationSchemeSDC : public TimeIntegrationSchemeSDC
{
public:
    IMEXTimeIntegrationSchemeSDC(std::string variant, size_t order,
                                 std::vector<NekDouble> freeParams)
        : TimeIntegrationSchemeSDC(variant, order, freeParams)
    {

        ASSERTL0(!m_first_quadrature,
                 "Quadrature type that include the left end point (e.g. "
                 "GaussLobattoLegendre) should not be used for IMEXSDC");

        std::cerr
            << "WARNING: IMEX Spectral Deferred Correction method has been "
               "implemented but its use is not recommended as the "
               "approach is affected by order-reduction problems."
            << std::endl;

        m_name       = "IMEXSpectralDeferredCorrection";
        m_schemeType = eIMEX;
    }

    static TimeIntegrationSchemeSharedPtr create(
        std::string variant, size_t order, std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<IMEXTimeIntegrationSchemeSDC>::AllocateSharedPtr(
                variant, order, freeParams);

        return p;
    }

    static std::string className;

    TripleArray m_Fexp; /// Array corresponding to the stage derivatives
    TripleArray m_Fimp; /// Array corresponding to the stage derivatives
    DoubleArray m_tmp;  /// Array for temporary storage

protected:
    LUE virtual void v_InitializeScheme(
        const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
        const TimeIntegrationSchemeOperators &op) override;

    LUE virtual void v_ResidualEval(const NekDouble &delta_t,
                                    const size_t n) override;

    LUE virtual void v_ResidualEval(const NekDouble &delta_t) override;

    LUE virtual void v_ComputeInitialGuess(const NekDouble &delta_t) override;

    LUE virtual void v_SDCIterationLoop(const NekDouble &delta_t) override;

private:
    void ComputeTotalResidual(const size_t n);

}; // end class IMEXTimeIntegrationSchemeSDC

/**
 * @brief Worker method to initialize the integration scheme.
 */
void IMEXTimeIntegrationSchemeSDC::v_InitializeScheme(
    const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
    const TimeIntegrationSchemeOperators &op)
{
    if (m_initialized)
    {
        m_time = time;
        for (size_t i = 0; i < m_nvars; ++i)
        {
            // Store the initial values as the first previous state.
            Vmath::Vcopy(m_npoints, y_0[i], 1, m_Y[0][i], 1);
        }
    }
    else
    {
        TimeIntegrationSchemeSDC::v_InitializeScheme(deltaT, y_0, time, op);

        m_Fexp = TripleArray(m_nQuadPts);
        m_Fimp = TripleArray(m_nQuadPts);
        for (size_t m = 0; m < m_nQuadPts; ++m)
        {
            m_Fexp[m] = DoubleArray(m_nvars);
            m_Fimp[m] = DoubleArray(m_nvars);
            for (size_t i = 0; i < m_nvars; ++i)
            {
                m_Fexp[m][i] = SingleArray(m_npoints, 0.0);
                m_Fimp[m][i] = SingleArray(m_npoints, 0.0);
            }
        }

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
void IMEXTimeIntegrationSchemeSDC::v_ResidualEval(const NekDouble &delta_t,
                                                  const size_t n)
{
    // Compute implicit residual
    if (n == 0)
    {
        // Not implemented, require implicit evaluation for m_Fimp[0].
        // Quadrature type that include the left end point (e.g.
        // GaussLobattoLegendre) should not be used.

        // Apply time-dependent boundary condition
        m_op.DoProjection(m_Y[0], m_Y[0], m_time);
    }
    else
    {
        NekDouble dtn = delta_t * (m_tau[n] - m_tau[n - 1]);

        // Update implicit solution
        m_op.DoImplicitSolve(m_Y[n - 1], m_tmp, m_time + delta_t * m_tau[n],
                             m_theta * dtn);

        // Compute implicit residual from updated solution
        for (size_t i = 0; i < m_nvars; ++i)
        {
            Vmath::Vsub(m_npoints, m_tmp[i], 1, m_Y[n - 1][i], 1, m_Fimp[n][i],
                        1);
            Vmath::Smul(m_npoints, 1.0 / (m_theta * dtn), m_Fimp[n][i], 1,
                        m_Fimp[n][i], 1);
        }
    }

    // Compute explicit residual
    m_op.DoOdeRhs(m_Y[n], m_Fexp[n], m_time + delta_t * m_tau[n]);

    // Compute total residual
    ComputeTotalResidual(n);
}

void IMEXTimeIntegrationSchemeSDC::v_ResidualEval(const NekDouble &delta_t)
{
    for (size_t n = 0; n < m_nQuadPts; ++n)
    {
        v_ResidualEval(delta_t, n);
    }
}

/**
 * @brief Worker method to compute the initial SDC guess.
 */
void IMEXTimeIntegrationSchemeSDC::v_ComputeInitialGuess(
    const NekDouble &delta_t)
{
    for (size_t n = 0; n < m_nQuadPts; ++n)
    {
        if (n == 0)
        {
            // Not implemented, require implicit evaluation for m_Fimp[0].
            // Quadrature type that include the left end point (e.g.
            // GaussLobattoLegendre) should not be used.

            // Apply time-dependent boundary condition
            m_op.DoProjection(m_Y[0], m_Y[0], m_time);
        }
        else
        {
            NekDouble dtn = delta_t * (m_tau[n] - m_tau[n - 1]);

            // Add explicit contribution to rhs
            for (size_t i = 0; i < m_nvars; ++i)
            {
                Vmath::Svtvp(m_npoints, dtn, m_Fexp[n - 1][i], 1, m_Y[n - 1][i],
                             1, m_tmp[i], 1);
            }

            // Solve implicit system from rhs
            m_op.DoImplicitSolve(m_tmp, m_Y[n], m_time + delta_t * m_tau[n],
                                 dtn);

            // Compute implicit flux from updated solution
            for (size_t i = 0; i < m_nvars; ++i)
            {
                Vmath::Vsub(m_npoints, m_Y[n][i], 1, m_tmp[i], 1, m_Fimp[n][i],
                            1);
                Vmath::Smul(m_npoints, 1.0 / dtn, m_Fimp[n][i], 1, m_Fimp[n][i],
                            1);
            }
        }

        // Compute explicit residual
        m_op.DoOdeRhs(m_Y[n], m_Fexp[n], m_time + delta_t * m_tau[n]);

        // Compute total residual
        ComputeTotalResidual(n);
    }
}

/**
 * @brief Worker method to compute the SDC iteration.
 */
void IMEXTimeIntegrationSchemeSDC::v_SDCIterationLoop(const NekDouble &delta_t)
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
            // Add SFint contribution to rhs
            Vmath::Vadd(m_npoints, m_Y[n - 1][i], 1, m_SFint[n][i], 1, m_tmp[i],
                        1);

            // Add explicit contribution to rhs
            if (n > 1)
            {
                Vmath::Svtvp(m_npoints, m_theta * dtn, m_Fexp[n - 1][i], 1,
                             m_tmp[i], 1, m_tmp[i], 1);
            }

            // Add explicit contribution to SFint
            if (n < m_nQuadPts - 1)
            {
                NekDouble dtnp = delta_t * (m_tau[n + 1] - m_tau[n]);
                Vmath::Svtvp(m_npoints, -m_theta * dtnp, m_Fexp[n][i], 1,
                             m_SFint[n + 1][i], 1, m_SFint[n + 1][i], 1);
            }

            // Add implicit contribution to rhs
            Vmath::Svtvp(m_npoints, -m_theta * dtn, m_Fimp[n][i], 1, m_tmp[i],
                         1, m_tmp[i], 1);
        }

        // Solve implicit system from rhs
        m_op.DoImplicitSolve(m_tmp, m_Y[n], m_time + delta_t * m_tau[n],
                             m_theta * dtn);

        // Compute implicit residual from updated solution
        for (size_t i = 0; i < m_nvars; ++i)
        {
            Vmath::Vsub(m_npoints, m_Y[n][i], 1, m_tmp[i], 1, m_Fimp[n][i], 1);
            Vmath::Smul(m_npoints, 1.0 / (m_theta * dtn), m_Fimp[n][i], 1,
                        m_Fimp[n][i], 1);
        }

        // Compute explicit residual
        m_op.DoOdeRhs(m_Y[n], m_Fexp[n], m_time + delta_t * m_tau[n]);

        // Compute total residual
        ComputeTotalResidual(n);
    }
}

/**
 * @brief Worker method to compute the total residual.
 */
void IMEXTimeIntegrationSchemeSDC::ComputeTotalResidual(const size_t n)
{
    for (size_t i = 0; i < m_nvars; ++i)
    {
        Vmath::Vadd(m_npoints, m_Fimp[n][i], 1, m_Fexp[n][i], 1, m_F[n][i], 1);
    }
}

} // end namespace LibUtilities
} // end namespace Nektar

#endif
