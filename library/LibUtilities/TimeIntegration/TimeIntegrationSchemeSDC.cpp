///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationSchemeSDC.cpp
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
// Description: implementation of time integration scheme SDC class
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeSDC.h>

namespace Nektar
{
namespace LibUtilities
{

const TripleArray &TimeIntegrationSchemeSDC::v_GetSolutionVector() const
{
    return m_Y;
}

TripleArray &TimeIntegrationSchemeSDC::v_UpdateSolutionVector()
{
    return m_Y;
}

void TimeIntegrationSchemeSDC::v_SetSolutionVector(const size_t Offset,
                                                   const DoubleArray &y)
{
    m_Y[Offset] = y;
}

std::string TimeIntegrationSchemeSDC::v_GetName() const
{
    return m_name;
}

std::string TimeIntegrationSchemeSDC::v_GetVariant() const
{
    return m_variant;
}

size_t TimeIntegrationSchemeSDC::v_GetOrder() const
{
    return m_order;
}

std::vector<NekDouble> TimeIntegrationSchemeSDC::v_GetFreeParams() const
{
    return m_freeParams;
}

TimeIntegrationSchemeType TimeIntegrationSchemeSDC::v_GetIntegrationSchemeType()
    const
{
    return m_schemeType;
}

NekDouble TimeIntegrationSchemeSDC::v_GetTimeStability() const
{
    return 1.0;
}

size_t TimeIntegrationSchemeSDC::v_GetNumIntegrationPhases() const
{
    return 1;
}

/**
 * @brief Worker method to initialize the integration scheme.
 */
void TimeIntegrationSchemeSDC::v_InitializeScheme(
    const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
    const TimeIntegrationSchemeOperators &op)
{
    boost::ignore_unused(op);
    boost::ignore_unused(deltaT);

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
        m_time    = time;
        m_nvars   = y_0.size();
        m_npoints = y_0[0].size();

        // Compute integration matrix.
        size_t colOffset = m_first_quadrature ? 0 : 1;
        size_t rowOffset = m_first_quadrature ? 0 : m_nQuadPts - 1;
        size_t nCols     = m_nQuadPts - colOffset;
        size_t nRows     = m_nQuadPts;
        m_QMat           = SingleArray(nRows * nCols, 0.0);
        Polylib::Qg(&m_QMat[rowOffset], &m_tau[colOffset], nCols);

        // Compute intepolation coefficient.
        m_interp = SingleArray(m_nQuadPts, 0.0);
        for (size_t i = 0; i < m_nQuadPts; i++)
        {
            m_interp[i] = Polylib::hgj(i, 1.0, &m_tau[0], m_nQuadPts, 0.0, 0.0);
        }

        // Storage of previous states and associated timesteps.
        m_Y     = TripleArray(m_nQuadPts);
        m_F     = TripleArray(m_nQuadPts);
        m_SFint = TripleArray(m_nQuadPts);
        m_QFint = TripleArray(m_nQuadPts);
        for (size_t m = 0; m < m_nQuadPts; ++m)
        {
            m_Y[m]     = DoubleArray(m_nvars);
            m_F[m]     = DoubleArray(m_nvars);
            m_SFint[m] = DoubleArray(m_nvars);
            m_QFint[m] = DoubleArray(m_nvars);
            for (size_t i = 0; i < m_nvars; ++i)
            {
                m_Y[m][i]     = SingleArray(m_npoints, 0.0);
                m_F[m][i]     = SingleArray(m_npoints, 0.0);
                m_SFint[m][i] = SingleArray(m_npoints, 0.0);
                m_QFint[m][i] = SingleArray(m_npoints, 0.0);
                // Store the initial values as the first previous state.
                if (m == 0)
                {
                    Vmath::Vcopy(m_npoints, y_0[i], 1, m_Y[m][i], 1);
                }
            }
        }

        if (!m_last_quadrature)
        {
            m_Y_f = DoubleArray(m_nvars);
            for (size_t i = 0; i < m_nvars; ++i)
            {
                m_Y_f[i] = SingleArray(m_npoints, 0.0);
            }
        }

        if (m_PFASST)
        {
            m_FAScorr = TripleArray(m_nQuadPts);
            for (size_t m = 0; m < m_nQuadPts; ++m)
            {
                m_FAScorr[m] = DoubleArray(m_nvars);
                for (size_t i = 0; i < m_nvars; ++i)
                {
                    m_FAScorr[m][i] = SingleArray(m_npoints, 0.0);
                }
            }
        }

        m_initialized = true;
    }
}

/**
 * @brief Worker method that performs the time integration.
 */
ConstDoubleArray &TimeIntegrationSchemeSDC::v_TimeIntegrate(
    const size_t timestep, const NekDouble delta_t,
    const TimeIntegrationSchemeOperators &op)
{
    boost::ignore_unused(timestep);

    for (size_t k = 0; k < m_order; ++k)
    {
        // Compute initial guess
        if (k == 0)
        {
            ComputeInitialGuess(delta_t, op);
        }
        // Apply SDC correction loop
        else
        {
            SDCIterationLoop(delta_t, op);
        }
    }

    // Update last quadrature
    UpdateLastQuadrature();

    // Update first quadrature
    UpdateFirstQuadrature();

    // Update time step
    m_time += delta_t;

    // Return solution
    return m_Y[0];
}

/**
 * @brief Worker method that update the first quadrature.
 */
void TimeIntegrationSchemeSDC::UpdateFirstQuadrature(void)
{
    DoubleArray Y_f = m_last_quadrature ? m_Y[m_nQuadPts - 1] : m_Y_f;
    for (size_t i = 0; i < m_nvars; ++i)
    {
        Vmath::Vcopy(m_npoints, Y_f[i], 1, m_Y[0][i], 1);
    }
}

/**
 * @brief Worker method that update the last quadrature.
 */
void TimeIntegrationSchemeSDC::UpdateLastQuadrature(void)
{
    if (!m_last_quadrature)
    {
        for (size_t i = 0; i < m_nvars; ++i)
        {
            Vmath::Zero(m_npoints, m_Y_f[i], 1);
            for (size_t n = 0; n < m_nQuadPts; ++n)
            {
                Vmath::Svtvp(m_npoints, m_interp[n], m_Y[n][i], 1, m_Y_f[i], 1,
                             m_Y_f[i], 1);
            }
        }
    }
}

/**
 * @brief Worker method that add the FASCorrection.
 */
void TimeIntegrationSchemeSDC::AddFASCorrectionToSFint(void)
{
    // Add PFASST correction term
    if (m_PFASST)
    {
        for (size_t n = 1; n < m_nQuadPts; ++n)
        {
            for (size_t i = 0; i < m_nvars; ++i)
            {
                Vmath::Vadd(m_npoints, m_SFint[n][i], 1, m_FAScorr[n][i], 1,
                            m_SFint[n][i], 1);
                Vmath::Vsub(m_npoints, m_SFint[n][i], 1, m_FAScorr[n - 1][i], 1,
                            m_SFint[n][i], 1);
            }
        }
    }
}

/**
 * @brief Worker method that compute residual integral.
 */
void TimeIntegrationSchemeSDC::UpdateIntegratedResidualSFint(
    const NekDouble &delta_t)
{
    // Zeroing integrated flux
    for (size_t n = 0; n < m_nQuadPts; ++n)
    {
        for (size_t i = 0; i < m_nvars; ++i)
        {
            Vmath::Zero(m_npoints, m_SFint[n][i], 1);
        }
    }

    // Update integrated flux
    size_t offset = m_first_quadrature ? 0 : 1;
    for (size_t n = 1; n < m_nQuadPts; ++n)
    {
        for (size_t p = 0; p < m_nQuadPts - offset; ++p)
        {
            for (size_t i = 0; i < m_nvars; ++i)
            {
                Vmath::Svtvp(
                    m_npoints,
                    delta_t * (m_QMat[n * (m_nQuadPts - offset) + p] -
                               m_QMat[(n - 1) * (m_nQuadPts - offset) + p]),
                    m_F[p + offset][i], 1, m_SFint[n][i], 1, m_SFint[n][i], 1);
            }
        }
    }
}

void TimeIntegrationSchemeSDC::UpdateIntegratedResidualQFint(
    const NekDouble &delta_t)
{
    // Zeroing integrated flux
    for (size_t n = 0; n < m_nQuadPts; ++n)
    {
        for (size_t i = 0; i < m_nvars; ++i)
        {
            Vmath::Zero(m_npoints, m_QFint[n][i], 1);
        }
    }

    // Update integrated flux
    size_t offset = m_first_quadrature ? 0 : 1;
    for (size_t n = 0; n < m_nQuadPts; ++n)
    {
        for (size_t p = 0; p < m_nQuadPts - offset; ++p)
        {
            for (size_t i = 0; i < m_nvars; ++i)
            {
                Vmath::Svtvp(
                    m_npoints, delta_t * m_QMat[n * (m_nQuadPts - offset) + p],
                    m_F[p + offset][i], 1, m_QFint[n][i], 1, m_QFint[n][i], 1);
            }
        }
    }
}

/**
 * @brief Worker method to print details on the integration scheme
 */
void TimeIntegrationSchemeSDC::v_print(std::ostream &os) const
{
    os << "Time Integration Scheme: " << GetFullName() << std::endl;
}

void TimeIntegrationSchemeSDC::v_printFull(std::ostream &os) const
{
    os << "Time Integration Scheme: " << GetFullName() << std::endl;
}

// Friend Operators
std::ostream &operator<<(std::ostream &os, const TimeIntegrationSchemeSDC &rhs)
{
    rhs.print(os);

    return os;
}

std::ostream &operator<<(std::ostream &os,
                         const TimeIntegrationSchemeSDCSharedPtr &rhs)
{
    os << *rhs.get();

    return os;
}

} // end namespace LibUtilities
} // namespace Nektar
