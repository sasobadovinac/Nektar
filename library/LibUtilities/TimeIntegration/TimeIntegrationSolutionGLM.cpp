///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationSolutionGLM.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeGLM.h>

namespace Nektar
{
namespace LibUtilities
{

TimeIntegrationSolutionGLM::TimeIntegrationSolutionGLM(
    const TimeIntegrationAlgorithmGLM *schemeAlgorithm, const DoubleArray &y,
    const NekDouble time, const NekDouble timestep)
    : m_schemeAlgorithm(schemeAlgorithm),
      m_solVector(m_schemeAlgorithm->m_numsteps),
      m_t(m_schemeAlgorithm->m_numsteps),
      m_setflag(m_schemeAlgorithm->m_numsteps)
{
    m_solVector[0] = y;
    m_t[0]         = time;
    m_setflag[0]   = true;

    size_t nsteps         = m_schemeAlgorithm->m_numsteps;
    size_t nvar           = y.size();
    size_t npoints        = y[0].size();
    size_t nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
    const Array<OneD, const size_t> &timeLevels =
        m_schemeAlgorithm->GetTimeLevelOffset();

    for (size_t i = 1; i < nsteps; i++)
    {
        m_solVector[i] = Array<OneD, Array<OneD, NekDouble>>(nvar);
        for (size_t j = 0; j < nvar; j++)
        {
            m_solVector[i][j] = Array<OneD, NekDouble>(npoints, 0.0);
        }

        if (i < nMultiStepVals)
        {
            m_t[i] = time - i * timestep * timeLevels[i];
        }
        else
        {
            m_t[i] = timestep;
        }

        m_setflag[i] = false;
    }
}

TimeIntegrationSolutionGLM::TimeIntegrationSolutionGLM(
    const TimeIntegrationAlgorithmGLM *schemeAlgorithm, const size_t nvar,
    const size_t npoints)
    : m_schemeAlgorithm(schemeAlgorithm),
      m_solVector(schemeAlgorithm->m_numsteps),
      m_t(schemeAlgorithm->m_numsteps),
      m_setflag(m_schemeAlgorithm->m_numsteps, true)
{
    for (size_t i = 0; i < m_schemeAlgorithm->m_numsteps; i++)
    {
        m_solVector[i] = Array<OneD, Array<OneD, NekDouble>>(nvar);
        for (size_t j = 0; j < nvar; j++)
        {
            m_solVector[i][j] = Array<OneD, NekDouble>(npoints);
        }
    }
}

TimeIntegrationSolutionGLM::TimeIntegrationSolutionGLM(
    const TimeIntegrationAlgorithmGLM *schemeAlgorithm)
    : m_schemeAlgorithm(schemeAlgorithm),
      m_solVector(m_schemeAlgorithm->m_numsteps),
      m_t(m_schemeAlgorithm->m_numsteps),
      m_setflag(m_schemeAlgorithm->m_numsteps, false)
{
}

} // end namespace LibUtilities
} // namespace Nektar
