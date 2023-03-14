///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationSchemeGEM.cpp
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
// Description: implementation of time integration scheme GEM class
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeGEM.h>

namespace Nektar
{
namespace LibUtilities
{

std::string TimeIntegrationSchemeGEM::v_GetName() const
{
    return m_name;
}

std::string TimeIntegrationSchemeGEM::v_GetVariant() const
{
    return m_variant;
}

unsigned int TimeIntegrationSchemeGEM::v_GetOrder() const
{
    return m_order;
}

std::vector<NekDouble> TimeIntegrationSchemeGEM::v_GetFreeParams() const
{
    return m_freeParams;
}

TimeIntegrationSchemeType TimeIntegrationSchemeGEM::v_GetIntegrationSchemeType()
    const
{
    return m_schemeType;
}

NekDouble TimeIntegrationSchemeGEM::v_GetTimeStability() const
{
    return 1.0;
}

unsigned int TimeIntegrationSchemeGEM::v_GetNumIntegrationPhases() const
{
    return 1;
}

/**
 * \brief Gets the solution vector of the ODE
 */
const TripleArray &TimeIntegrationSchemeGEM::v_GetSolutionVector() const
{
    return m_Y;
}

/**
 * \brief Sets the solution vector of the ODE
 */
void TimeIntegrationSchemeGEM::v_SetSolutionVector(const int Offset,
                                                   const DoubleArray &y)
{
    m_Y[Offset] = y;
}

/**
 * @brief Worker method to initialize the integration scheme.
 */
void TimeIntegrationSchemeGEM::v_InitializeScheme(
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

        int nodes = m_order;
        if (m_variant == "ExplicitMidpoint" || m_variant == "ImplicitMidpoint")
        {
            nodes /= 2;
        }

        // Storage of previous states and associated timesteps.
        m_Y  = TripleArray(m_order + 1);
        m_T  = TripleArray(nodes);
        m_T0 = TripleArray(nodes);

        for (unsigned int m = 0; m <= m_order; ++m)
        {
            m_Y[m] = DoubleArray(m_nvars);

            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                m_Y[m][i] = SingleArray(m_npoints, 0.0);

                // Store the initial values as the first previous state.
                if (m == 0)
                {
                    Vmath::Vcopy(m_npoints, y_0[i], 1, m_Y[m][i], 1);
                }
            }
        }

        if (m_variant == "" || m_variant == "ExplicitEuler" ||
            m_variant == "ExplicitMidpoint")
        {
            op.DoProjection(m_Y[0], m_Y[0], m_time);
        }

        for (unsigned int m = 0; m < nodes; ++m)
        {
            m_T[m]  = DoubleArray(m_nvars);
            m_T0[m] = DoubleArray(m_nvars);

            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                m_T[m][i]  = SingleArray(m_npoints, 0.0);
                m_T0[m][i] = SingleArray(m_npoints, 0.0);
            }
        }

        // Storage for the stage derivative as the data will be re-used to
        // update the solution.
        m_F  = DoubleArray(m_nvars);
        m_F0 = DoubleArray(m_nvars);

        for (unsigned int i = 0; i < m_nvars; ++i)
        {
            m_F[i]  = SingleArray(m_npoints, 0.0);
            m_F0[i] = SingleArray(m_npoints, 0.0);
        }

        m_initialized = true;
    }
}

/**
 * @brief Worker method that performs the time integration.
 */
ConstDoubleArray &TimeIntegrationSchemeGEM::v_TimeIntegrate(
    const int timestep, const NekDouble delta_t,
    const TimeIntegrationSchemeOperators &op)
{

    boost::ignore_unused(timestep);

    // Compute initial residual
    if (m_variant == "" || m_variant == "ExplicitEuler" ||
        m_variant == "ExplicitMidpoint" || m_variant == "IMEXEuler")
    {
        op.DoOdeRhs(m_Y[0], m_F0, m_time);
    }

    // Euler approach
    if (m_variant == "" || m_variant == "ExplicitEuler" ||
        m_variant == "ImplicitEuler" || m_variant == "IMEXEuler")
    {
        // Compute first order approximation
        for (unsigned int m = 1; m <= m_order; ++m)
        {
            for (unsigned int k = 1; k <= m; ++k)
            {
                // Implicit schemes
                if (m_variant == "ImplicitEuler")
                {
                    op.DoImplicitSolve(m_Y[k - 1], m_Y[k],
                                       m_time + k * (delta_t / m), delta_t / m);
                }

                // Explicit schemes
                if (m_variant == "" || m_variant == "ExplicitEuler" ||
                    m_variant == "IMEXEuler")
                {
                    // For the first stage, used pre-computed rhs
                    if (k == 1)
                    {
                        for (unsigned int i = 0; i < m_nvars; ++i)
                        {
                            Vmath::Svtvp(m_npoints, delta_t / m, m_F0[i], 1,
                                         m_Y[k - 1][i], 1, m_Y[k][i], 1);
                        }
                    }
                    // For other stages, compute new rhs
                    else
                    {
                        op.DoOdeRhs(m_Y[k - 1], m_F,
                                    m_time + (k - 1) * (delta_t / m));
                        for (unsigned int i = 0; i < m_nvars; ++i)
                        {
                            Vmath::Svtvp(m_npoints, delta_t / m, m_F[i], 1,
                                         m_Y[k - 1][i], 1, m_Y[k][i], 1);
                        }
                    }
                }
                if (m_variant == "" || m_variant == "ExplicitEuler")
                {
                    op.DoProjection(m_Y[k], m_Y[k], m_time + k * (delta_t / m));
                }

                // IMEX schemes (NOTE: Order reduction problems)
                if (m_variant == "IMEXEuler")
                {
                    op.DoImplicitSolve(m_Y[k], m_Y[k],
                                       m_time + k * (delta_t / m), delta_t / m);
                }
            }

            // Save solution to m_T0
            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                Vmath::Vcopy(m_npoints, m_Y[m][i], 1, m_T0[m - 1][i], 1);
            }
        }

        // No extrapolation required for first-order
        if (m_order == 1)
        {
            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                Vmath::Vcopy(m_npoints, m_Y[1][i], 1, m_Y[0][i], 1);
            }
            m_time += delta_t;
            return m_Y[0];
        }

        // Extrapolate solution
        for (unsigned int m = 1; m < m_order; ++m)
        {
            // Aitken - Neville formula
            for (unsigned int k = m; k < m_order; ++k)
            {
                for (unsigned int i = 0; i < m_nvars; ++i)
                {
                    Vmath::Vsub(m_npoints, m_T0[k][i], 1, m_T0[k - 1][i], 1,
                                m_T[k][i], 1);
                    Vmath::Svtvp(m_npoints,
                                 (k - m + 1.0) / ((k + 1.0) - (k - m + 1.0)),
                                 m_T[k][i], 1, m_T0[k][i], 1, m_T[k][i], 1);
                }
            }

            // Copy new values to old values
            for (unsigned int k = m; k < m_order; ++k)
            {
                for (unsigned int i = 0; i < m_nvars; ++i)
                {
                    Vmath::Vcopy(m_npoints, m_T[k][i], 1, m_T0[k][i], 1);
                }
            }
        }

        // Copy final solution
        for (unsigned int i = 0; i < m_nvars; ++i)
        {
            Vmath::Vcopy(m_npoints, m_T[m_order - 1][i], 1, m_Y[0][i], 1);
        }
    }
    // Midpoint approach
    else if (m_variant == "ExplicitMidpoint" || m_variant == "ImplicitMidpoint")
    {
        // Compute second order approximation
        for (unsigned int m = 1; m <= m_order / 2; ++m)
        {
            // Implicit midpoint
            if (m_variant == "ImplicitMidpoint")
            {
                for (unsigned int k = 1; k <= m; ++k)
                {
                    op.DoImplicitSolve(m_Y[2 * k - 2], m_Y[2 * k - 1],
                                       m_time + (k - 1 + 0.25) * (delta_t / m),
                                       0.25 * delta_t / m);
                    for (unsigned int i = 0; i < m_nvars; ++i)
                    {
                        Vmath::Svtsvtp(m_npoints, 2.0, m_Y[2 * k - 1][i], 1,
                                       -1.0, m_Y[2 * k - 2][i], 1,
                                       m_Y[2 * k][i], 1);
                    }
                    op.DoImplicitSolve(m_Y[2 * k], m_F,
                                       m_time + (k - 0.25) * (delta_t / m),
                                       0.25 * delta_t / m);
                    for (unsigned int i = 0; i < m_nvars; ++i)
                    {
                        Vmath::Vsub(m_npoints, m_F[i], 1, m_Y[2 * k][i], 1,
                                    m_F[i], 1);
                        Vmath::Svtvp(m_npoints, 2.0, m_F[i], 1, m_Y[2 * k][i],
                                     1, m_Y[2 * k][i], 1);
                    }
                }
            }

            // Explicit midpoint
            if (m_variant == "ExplicitMidpoint")
            {
                // Use precomputed rhs for initial Euler stage
                for (unsigned int i = 0; i < m_nvars; ++i)
                {
                    Vmath::Svtvp(m_npoints, delta_t / (2 * m), m_F0[i], 1,
                                 m_Y[0][i], 1, m_Y[1][i], 1);
                }
                op.DoProjection(m_Y[1], m_Y[1], m_time + delta_t / (2 * m));

                // Compute new rhs for midpoint stage
                for (unsigned int k = 2; k <= 2 * m; ++k)
                {
                    op.DoOdeRhs(m_Y[k - 1], m_F,
                                m_time + (k - 1) * (delta_t / (2 * m)));
                    for (unsigned int i = 0; i < m_nvars; ++i)
                    {
                        Vmath::Svtvp(m_npoints, delta_t / m, m_F[i], 1,
                                     m_Y[k - 2][i], 1, m_Y[k][i], 1);
                    }
                    op.DoProjection(m_Y[k], m_Y[k],
                                    m_time + k * (delta_t / (2 * m)));
                }
            }

            // Save solution to m_T0
            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                Vmath::Vcopy(m_npoints, m_Y[2 * m][i], 1, m_T0[m - 1][i], 1);
            }
        }

        // No extrapolation required for second-order
        if (m_order == 2)
        {
            for (unsigned int i = 0; i < m_nvars; ++i)
            {
                Vmath::Vcopy(m_npoints, m_Y[2][i], 1, m_Y[0][i], 1);
            }
            m_time += delta_t;
            return m_Y[0];
        }

        // Extrapolate solution
        for (unsigned int m = 1; m < m_order / 2; ++m)
        {
            // Aitken - Neville formula
            for (unsigned int k = m; k < m_order / 2; ++k)
            {
                for (unsigned int i = 0; i < m_nvars; ++i)
                {
                    Vmath::Vsub(m_npoints, m_T0[k][i], 1, m_T0[k - 1][i], 1,
                                m_T[k][i], 1);
                    Vmath::Svtvp(
                        m_npoints,
                        std::pow(k - m + 1.0, 2) /
                            (std::pow(k + 1.0, 2) - std::pow(k - m + 1.0, 2)),
                        m_T[k][i], 1, m_T0[k][i], 1, m_T[k][i], 1);
                }
            }

            // Copy new values to old values
            for (unsigned int k = m; k < m_order / 2; ++k)
            {
                for (unsigned int i = 0; i < m_nvars; ++i)
                {
                    Vmath::Vcopy(m_npoints, m_T[k][i], 1, m_T0[k][i], 1);
                }
            }
        }

        // Copy final solution
        for (unsigned int i = 0; i < m_nvars; ++i)
        {
            Vmath::Vcopy(m_npoints, m_T[m_order / 2 - 1][i], 1, m_Y[0][i], 1);
        }
    }

    // Return solution
    m_time += delta_t;
    return m_Y[0];
}

/**
 * @brief Worker method to print details on the integration scheme
 */
void TimeIntegrationSchemeGEM::v_print(std::ostream &os) const
{
    os << "Time Integration Scheme: " << GetFullName() << std::endl;
}

void TimeIntegrationSchemeGEM::v_printFull(std::ostream &os) const
{
    os << "Time Integration Scheme: " << GetFullName() << std::endl;
}

// Friend Operators
std::ostream &operator<<(std::ostream &os, const TimeIntegrationSchemeGEM &rhs)
{
    rhs.print(os);

    return os;
}

std::ostream &operator<<(std::ostream &os,
                         const TimeIntegrationSchemeGEMSharedPtr &rhs)
{
    os << *rhs.get();

    return os;
}
} // end namespace LibUtilities
} // namespace Nektar
