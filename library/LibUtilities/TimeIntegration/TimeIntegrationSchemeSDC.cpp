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
        m_coeffs = SingleArray(m_nQuadPts, 0.0);
        m_wMat   = DoubleArray(m_nQuadPts - 1);
        for (unsigned int n = 0; n < m_nQuadPts - 1; ++n)
        {
            m_wMat[n] = SingleArray(m_nQuadPts, 0.0);
        }

        // Use first quadrature point in the integral
        if (m_variant == "Equidistant" || m_variant == "GaussLobattoLegendre" ||
            m_variant == "GaussLobattoChebyshev")
        {
            for (unsigned int m = 0; m < m_nQuadPts; ++m)
            {
                EvaluateCoeffs(&m_coeffs[0], &m_points[0], m, m_nQuadPts);
                for (unsigned int n = 0; n < m_nQuadPts - 1; ++n)
                {
                    m_wMat[n][m] = EvaluateInt(&m_coeffs[0], m_nQuadPts,
                                               m_points[n], m_points[n + 1]);
                }
            }
        }
        // Do not use first quadrature point in the integral
        else if (m_variant == "GaussRadauLegendre" ||
                 m_variant == "GaussRadauChebyshev")
        {
            for (unsigned int m = 0; m < m_nQuadPts - 1; ++m)
            {
                EvaluateCoeffs(&m_coeffs[1], &m_points[1], m, m_nQuadPts - 1);
                for (unsigned int n = 0; n < m_nQuadPts - 1; ++n)
                {
                    m_wMat[n][m + 1] =
                        EvaluateInt(&m_coeffs[1], m_nQuadPts - 1, m_points[n],
                                    m_points[n + 1]);
                }
            }
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

ConstDoubleArray &TimeIntegrationSchemeSDC::v_TimeIntegrate(
    const int timestep, const NekDouble delta_t,
    const TimeIntegrationSchemeOperators &op)
{
    boost::ignore_unused(timestep);
    boost::ignore_unused(delta_t);
    boost::ignore_unused(op);

    ASSERTL0(
        false,
        "Specific version of spectral deferred correction not implemented");

    return m_Y[m_nQuadPts - 1];
}

void TimeIntegrationSchemeSDC::UpdateIntegratedFlux(const NekDouble &delta_t)
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
    if (m_variant == "Equidistant" || m_variant == "GaussLobattoLegendre" ||
        m_variant == "GaussLobattoChebyshev")
    {
        // Use first quadrature point in the integral
        for (unsigned int n = 0; n < m_nQuadPts - 1; ++n)
        {
            for (unsigned int p = 0; p < m_nQuadPts; ++p)
            {
                for (unsigned int i = 0; i < m_nvars; ++i)
                {
                    Vmath::Svtvp(m_npoints, delta_t * m_wMat[n][p], m_F[p][i],
                                 1, m_Fint[n][i], 1, m_Fint[n][i], 1);
                }
            }
        }
    }
    else if (m_variant == "GaussRadauLegendre" ||
             m_variant == "GaussRadauChebyshev")
    {
        // Do not use first quadrature point in the integral
        for (unsigned int n = 0; n < m_nQuadPts - 1; ++n)
        {
            for (unsigned int p = 0; p < m_nQuadPts - 1; ++p)
            {
                for (unsigned int i = 0; i < m_nvars; ++i)
                {
                    Vmath::Svtvp(m_npoints, delta_t * m_wMat[n][p + 1],
                                 m_F[p + 1][i], 1, m_Fint[n][i], 1,
                                 m_Fint[n][i], 1);
                }
            }
        }
    }
}

/**
 * @brief Worker method that evaluates integral.
 */
NekDouble TimeIntegrationSchemeSDC::EvaluateInt(NekDouble *coeffs,
                                                const int npoints,
                                                const NekDouble &t1,
                                                const NekDouble &t2) const
{
    NekDouble val = 0.0;
    for (int n = 0; n < npoints; n++)
    {
        val +=
            coeffs[n] / (n + 1) * (std::pow(t2, n + 1) - std::pow(t1, n + 1));
    }
    return val;
}

/**
 * @brief Worker method that evaluates polynomial coefficients.
 */
void TimeIntegrationSchemeSDC::EvaluateCoeffs(NekDouble *coeffs,
                                              NekDouble *points,
                                              const int index,
                                              const int npoints)
{
    SingleArray tmp1(npoints, 0.0);
    SingleArray tmp2(npoints, 0.0);
    Vmath::Zero(npoints, coeffs, 1);

    // Compute denominator
    NekDouble denom = 1.0;
    for (int n = 0; n < npoints; n++)
    {
        if (index != n)
        {
            denom *= (points[index] - points[n]);
        }
    }

    // Compute coefficient
    coeffs[0] = 1.0;
    int m     = 0;
    for (int n = 0; n < npoints; n++)
    {
        if (index != n)
        {
            m += 1;
            Vmath::Zero(npoints, tmp1, 1);
            Vmath::Zero(npoints, tmp2, 1);
            Vmath::Smul(npoints, -1.0 * points[n], &coeffs[0], 1, &tmp1[0], 1);
            Vmath::Vcopy(m, coeffs, 1, &tmp2[1], 1);
            Vmath::Vadd(npoints, &tmp1[0], 1, &tmp2[0], 1, &coeffs[0], 1);
        }
    }
    Vmath::Smul(npoints, 1.0 / denom, &coeffs[0], 1, &coeffs[0], 1);
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
