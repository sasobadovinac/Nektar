///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationSchemeGLM.cpp
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
// Description: implementation of time integration scheme GLM class
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeGLM.h>

namespace Nektar
{
namespace LibUtilities
{

// Access Methods
std::string TimeIntegrationSchemeGLM::v_GetVariant() const
{
    ASSERTL0(!m_integration_phases.empty(), "No scheme")

    return m_integration_phases.back()->m_variant;
}

size_t TimeIntegrationSchemeGLM::v_GetOrder() const
{
    ASSERTL0(!m_integration_phases.empty(), "No scheme")

    return m_integration_phases.back()->m_order;
}

std::vector<NekDouble> TimeIntegrationSchemeGLM::v_GetFreeParams() const
{
    ASSERTL0(!m_integration_phases.empty(), "No scheme")

    return m_integration_phases.back()->m_freeParams;
}

TimeIntegrationSchemeType TimeIntegrationSchemeGLM::v_GetIntegrationSchemeType()
    const
{
    ASSERTL0(!m_integration_phases.empty(), "No scheme")

    return m_integration_phases.back()->m_schemeType;
}

size_t TimeIntegrationSchemeGLM::v_GetNumIntegrationPhases() const
{
    return m_integration_phases.size();
}

/**
 * @brief Worker method to initialize the integration scheme.
 */
void TimeIntegrationSchemeGLM::v_InitializeScheme(
    const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
    const TimeIntegrationSchemeOperators &op)
{
    for (int i = 0; i < m_integration_phases.size(); i++)
    {
        m_integration_phases[i]->InitializeScheme(op);
    }

    m_solVector =
        m_integration_phases.back()->InitializeData(deltaT, y_0, time);
}

/**
 * @brief Worker method that actually does the time integration.
 */
ConstDoubleArray &TimeIntegrationSchemeGLM::v_TimeIntegrate(
    const size_t timestep, const NekDouble delta_t)
{
    size_t nPhases = m_integration_phases.size();

    TimeIntegrationAlgorithmGLMSharedPtr &algorithm =
        m_integration_phases[std::min(timestep, nPhases - 1)];

    return algorithm->TimeIntegrate(delta_t, m_solVector);
}

void TimeIntegrationSchemeGLM::v_InitializeSecondaryData(
    TimeIntegrationAlgorithmGLM *phase, NekDouble deltaT) const
{
    boost::ignore_unused(phase, deltaT);

    ASSERTL0(false,
             "No InitializeSecondaryData method for scheme " + GetFullName());
}

/**
 * @brief Worker method to print details on the integration scheme
 */
void TimeIntegrationSchemeGLM::v_print(std::ostream &os) const
{
    os << "Time Integration Scheme: " << GetFullName() << std::endl
       << "        Has " << m_integration_phases.size() << " phase(s)"
       << std::endl;

    for (size_t i = 0; i < m_integration_phases.size(); i++)
    {
        os << "            - " << m_integration_phases[i]->m_name << std::endl;
    }
}

void TimeIntegrationSchemeGLM::v_printFull(std::ostream &os) const
{
    os << "Time Integration Scheme: " << GetFullName() << std::endl
       << "        Has " << m_integration_phases.size() << " phase(s)"
       << std::endl;

    for (size_t i = 0; i < m_integration_phases.size(); i++)
    {
        os << "            - " << m_integration_phases[i]->m_name << std::endl;

        os << "            - " << m_integration_phases[i] << std::endl;
    }
}

// Friend Operators
std::ostream &operator<<(std::ostream &os, const TimeIntegrationSchemeGLM &rhs)
{
    rhs.print(os);

    return os;
}

std::ostream &operator<<(std::ostream &os,
                         const TimeIntegrationSchemeGLMSharedPtr &rhs)
{
    os << *rhs.get();

    return os;
}

} // end namespace LibUtilities
} // namespace Nektar
