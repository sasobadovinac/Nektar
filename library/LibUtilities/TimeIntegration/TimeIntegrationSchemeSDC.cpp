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

std::string TimeIntegrationSchemeSDC::v_GetName() const
{
    return m_name;
}

std::string TimeIntegrationSchemeSDC::v_GetVariant() const
{
    return m_variant;
}

unsigned int TimeIntegrationSchemeSDC::v_GetOrder() const
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

unsigned int TimeIntegrationSchemeSDC::v_GetNumIntegrationPhases() const
{
    return 1;
}

/**
 * \brief Gets the solution vector of the ODE
 */
const TripleArray &TimeIntegrationSchemeSDC::v_GetSolutionVector() const
{
    return m_Y;
}

/**
 * \brief Sets the solution vector of the ODE
 */
void TimeIntegrationSchemeSDC::v_SetSolutionVector(const int Offset,
                                                   const DoubleArray &y)
{
    m_Y[Offset] = y;
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
        for (unsigned int i = 0; i < m_nvars; ++i)
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

        // Compute integration matrix
        m_QMat = SingleArray(m_nQuadPts * m_nQuadPts, 0.0);

        // Use first quadrature point in the integral
        if (m_variant == "Equidistant" || m_variant == "GaussLobattoLegendre" ||
            m_variant == "GaussLobattoChebyshev")
        {
            Polylib::Qg(&m_QMat[0], &m_points[0], m_nQuadPts, 0);
        }
        else if (m_variant == "GaussRadauLegendre" ||
                 m_variant == "GaussRadauChebyshev")
        {
            Polylib::Qg(&m_QMat[m_nQuadPts], &m_points[1], m_nQuadPts - 1, 1);
        }

        // Storage of previous states and associated timesteps.
        m_Y = TripleArray(m_nQuadPts);
        m_F = TripleArray(m_nQuadPts);

        for (unsigned int m = 0; m < m_nQuadPts; ++m)
        {
            m_Y[m] = DoubleArray(m_nvars);
            m_F[m] = DoubleArray(m_nvars);

            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                m_Y[m][i] = SingleArray(m_npoints, 0.0);
                m_F[m][i] = SingleArray(m_npoints, 0.0);

                // Store the initial values as the first previous state.
                if (m == 0)
                {
                    Vmath::Vcopy(m_npoints, y_0[i], 1, m_Y[m][i], 1);
                }
            }
        }

        m_Fint = TripleArray(m_nQuadPts - 1);
        for (unsigned int m = 0; m < m_nQuadPts - 1; ++m)
        {
            m_Fint[m] = DoubleArray(m_nvars);

            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                m_Fint[m][i] = SingleArray(m_npoints, 0.0);
            }
        }

        m_Fn  = DoubleArray(m_nvars);
        m_tmp = DoubleArray(m_nvars);
        for (unsigned int i = 0; i < m_nvars; ++i)
        {
            m_Fn[i]  = SingleArray(m_npoints, 0.0);
            m_tmp[i] = SingleArray(m_npoints, 0.0);
        }

        m_initialized = true;
    }
}

/**
 * @brief Worker method that performs the time integration.
 */
ConstDoubleArray &TimeIntegrationSchemeSDC::v_TimeIntegrate(
    const int timestep, const NekDouble delta_t,
    const TimeIntegrationSchemeOperators &op)
{

    boost::ignore_unused(timestep);

    // 1. Compute initial guess
    ComputeInitialGuess(delta_t, op);

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
    v_SDCIterationLoop(delta_t, op);

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

/**
 * @brief Worker method that compute the integrated flux.
 */
void TimeIntegrationSchemeSDC::UpdateIntegratedResidual(
    const NekDouble &delta_t, const int option)
{
    // Zeroing integrated flux
    for (unsigned int n = 0; n < m_nQuadPts - 1; ++n)
    {
        for (unsigned int i = 0; i < m_nvars; ++i)
        {
            Vmath::Zero(m_npoints, m_Fint[n][i], 1);
        }
    }

    // Update integrated flux
    for (unsigned int n = 0; n < m_nQuadPts - 1; ++n)
    {
        for (unsigned int p = 0; p < m_nQuadPts; ++p)
        {
            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                if (option == 0)
                {
                    Vmath::Svtvp(m_npoints,
                                 delta_t * (m_QMat[(n + 1) * m_nQuadPts + p] -
                                            m_QMat[n * m_nQuadPts + p]),
                                 m_F[p][i], 1, m_Fint[n][i], 1, m_Fint[n][i],
                                 1);
                }
                else if (option == 1)
                {
                    Vmath::Svtvp(
                        m_npoints, delta_t * m_QMat[(n + 1) * m_nQuadPts + p],
                        m_F[p][i], 1, m_Fint[n][i], 1, m_Fint[n][i], 1);
                }
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
