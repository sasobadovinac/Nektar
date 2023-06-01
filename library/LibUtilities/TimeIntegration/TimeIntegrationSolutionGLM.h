///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationSolutionGLM.h
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
// Description: Header file of time integration solution class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_SOLUTION_GLM
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_SOLUTION_GLM

#define LUE LIB_UTILITIES_EXPORT

#include <LibUtilities/TimeIntegration/TimeIntegrationAlgorithmGLM.h>

namespace Nektar
{
namespace LibUtilities
{

class TimeIntegrationSolutionGLM
{
public:
    LUE TimeIntegrationSolutionGLM(
        const TimeIntegrationAlgorithmGLM *schemeAlgorithm,
        const DoubleArray &y, const NekDouble time, const NekDouble timestep);

    LUE TimeIntegrationSolutionGLM(
        const TimeIntegrationAlgorithmGLM *schemeAlgorithm, const size_t nvar,
        const size_t npoints);

    LUE TimeIntegrationSolutionGLM(
        const TimeIntegrationAlgorithmGLM *schemeAlgorithm);

    inline const TimeIntegrationAlgorithmGLM *GetIntegrationSchemeData() const
    {
        return m_schemeAlgorithm;
    }

    inline const TripleArray &GetSolutionVector() const
    {
        return m_solVector;
    }

    inline TripleArray &UpdateSolutionVector()
    {
        return m_solVector;
    }

    inline const DoubleArray &GetSolution() const
    {
        return m_solVector[0];
    }

    inline DoubleArray &UpdateSolution()
    {
        return m_solVector[0];
    }

    inline void SetSolutionVector(const size_t Offset, const DoubleArray &y)
    {
        m_solVector[Offset] = y;
    }

    inline const Array<OneD, const NekDouble> &GetTimeVector() const
    {
        return m_t;
    }

    inline Array<OneD, NekDouble> &UpdateTimeVector()
    {
        return m_t;
    }

    inline NekDouble GetTime() const
    {
        return m_t[0];
    }

    size_t GetNsteps()
    {
        return m_schemeAlgorithm->m_numsteps;
    }

    inline size_t GetFirstDim() const
    {
        return m_solVector[0].size();
    }

    inline size_t GetSecondDim() const
    {
        return m_solVector[0][0].size();
    }

    // Return the number of entries in the solution vector that correspond to
    // (multi-step) values.
    inline size_t GetNvalues() const
    {
        return m_schemeAlgorithm->GetNmultiStepValues();
    }

    // Return the number of entries in the solution vector that correspond to
    // implicit (multi-step) derivatives.
    inline size_t GetNimplicitderivs() const
    {
        return m_schemeAlgorithm->GetNmultiStepImplicitDerivs();
    }

    // Return the number of entries in the solution vector that correspond to
    // explicit (multi-step) derivatives.
    inline size_t GetNexplicitderivs() const
    {
        return m_schemeAlgorithm->GetNmultiStepExplicitDerivs();
    }

    // Returns an array which indicates to which time-level the entries in the
    // solution vector correspond.
    inline const Array<OneD, const size_t> &GetTimeLevelOffset()
    {
        return m_schemeAlgorithm->GetTimeLevelOffset();
    }

    // Returns the entry in the solution vector which corresponds to the
    // (multi-step) value at the time-level with specified offset
    inline DoubleArray &GetValue(const size_t timeLevelOffset)
    {
        size_t nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
        const Array<OneD, const size_t> &offsetvec =
            m_schemeAlgorithm->GetTimeLevelOffset();

        for (size_t i = 0; i < nMultiStepVals; i++)
        {
            if (timeLevelOffset == offsetvec[i])
            {
                ASSERTL0(m_setflag[i],
                         "Solution vector is not set at this time level");
                return m_solVector[i];
            }
        }
        ASSERTL1(false, "The solution vector of this scheme does not contain a "
                        "value at the requested time-level");
        return m_solVector[0];
    }

    // returns the entry in the solution vector which corresponds to the
    // implicit (multi-step) derivative at the time-level with specified offset
    inline DoubleArray &GetImplicitDerivative(const size_t timeLevelOffset)
    {
        size_t nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
        size_t nMultiStepImplicitDerivs =
            m_schemeAlgorithm->GetNmultiStepImplicitDerivs();
        const Array<OneD, const size_t> &offsetvec =
            m_schemeAlgorithm->GetTimeLevelOffset();

        for (size_t i = nMultiStepVals;
             i < nMultiStepVals + nMultiStepImplicitDerivs; i++)
        {
            if (timeLevelOffset == offsetvec[i])
            {
                ASSERTL0(m_setflag[i],
                         "Implicit derivative solution vector is not set at "
                         "this time level");
                return m_solVector[i];
            }
        }
        ASSERTL1(false, "The solution vector of this scheme does not contain a "
                        "derivative at the requested time-level");
        return m_solVector[0];
    }

    // returns the entry in the solution vector which corresponds to the
    // explicit (multi-step) derivative at the time-level with specified offset
    inline DoubleArray &GetExplicitDerivative(const size_t timeLevelOffset)
    {
        size_t nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
        size_t nMultiStepImplicitDerivs =
            m_schemeAlgorithm->GetNmultiStepImplicitDerivs();
        size_t size = m_schemeAlgorithm->m_numsteps;
        const Array<OneD, const size_t> &offsetvec =
            m_schemeAlgorithm->GetTimeLevelOffset();

        for (size_t i = nMultiStepVals + nMultiStepImplicitDerivs; i < size;
             i++)
        {
            if (timeLevelOffset == offsetvec[i])
            {
                ASSERTL0(m_setflag[i],
                         "Explicit derivative solution vector is not set at "
                         "this time level");
                return m_solVector[i];
            }
        }
        ASSERTL1(false, "The solution vector of this scheme does not contain a "
                        "derivative at the requested time-level");
        return m_solVector[0];
    }

    // returns the time associated with the (multi-step) value at the time-level
    // with the given offset
    inline NekDouble GetValueTime(const size_t timeLevelOffset)
    {
        size_t nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
        const Array<OneD, const size_t> &offsetvec =
            m_schemeAlgorithm->GetTimeLevelOffset();

        for (size_t i = 0; i < nMultiStepVals; i++)
        {
            if (timeLevelOffset == offsetvec[i])
            {
                return m_t[i];
            }
        }
        ASSERTL1(false, "The solution vector of this scheme does not contain a "
                        "value at the requested time-level");
        return m_t[0];
    }

    // sets the (multi-step) value and time in the solution vector which
    // corresponds to the value at the time-level with specified offset
    inline void SetValue(const size_t timeLevelOffset, const DoubleArray &y,
                         const NekDouble t)
    {
        size_t nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
        const Array<OneD, const size_t> &offsetvec =
            m_schemeAlgorithm->GetTimeLevelOffset();

        for (size_t i = 0; i < nMultiStepVals; i++)
        {
            if (timeLevelOffset == offsetvec[i])
            {
                m_solVector[i] = y;
                m_t[i]         = t;
                m_setflag[i]   = true;
                return;
            }
        }
    }

    // sets the (multi-step) derivative and time in the solution vector which
    // corresponds to the implicit derivative at the time-level with specified
    // offset
    inline void SetImplicitDerivative(const size_t timeLevelOffset,
                                      const DoubleArray &y,
                                      const NekDouble timestep)
    {
        size_t nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
        size_t nMultiStepImplicitDerivs =
            m_schemeAlgorithm->GetNmultiStepImplicitDerivs();
        const Array<OneD, const size_t> &offsetvec =
            m_schemeAlgorithm->GetTimeLevelOffset();

        for (size_t i = nMultiStepVals;
             i < nMultiStepVals + nMultiStepImplicitDerivs; i++)
        {
            if (timeLevelOffset == offsetvec[i])
            {
                m_solVector[i] = y;
                m_t[i]         = timestep;
                m_setflag[i]   = true;
                return;
            }
        }
    }

    // sets the (multi-step) derivative and time in the solution vector which
    // corresponds to the explicit derivative at the time-level with specified
    // offset
    inline void SetExplicitDerivative(const size_t timeLevelOffset,
                                      const DoubleArray &y,
                                      const NekDouble timestep)
    {
        size_t nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
        size_t nMultiStepImplicitDerivs =
            m_schemeAlgorithm->GetNmultiStepImplicitDerivs();
        size_t size = m_schemeAlgorithm->m_numsteps;
        const Array<OneD, const size_t> &offsetvec =
            m_schemeAlgorithm->GetTimeLevelOffset();

        for (size_t i = nMultiStepVals + nMultiStepImplicitDerivs; i < size;
             i++)
        {
            if (timeLevelOffset == offsetvec[i])
            {
                m_solVector[i] = y;
                m_t[i]         = timestep;
                m_setflag[i]   = true;
                return;
            }
        }
    }

    // Rotate the solution vector
    // (i.e. updating without calculating/inserting new values)
    inline void RotateSolutionVector()
    {
        size_t nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
        size_t nMultiStepImpDerivs =
            m_schemeAlgorithm->GetNmultiStepImplicitDerivs();
        size_t size = m_schemeAlgorithm->m_numsteps;
        for (size_t i = (nMultiStepVals - 1); i > 0; i--)
        {
            m_solVector[i] = m_solVector[i - 1];
            m_setflag[i]   = m_setflag[i - 1];
        }

        for (size_t i = (nMultiStepVals + nMultiStepImpDerivs - 1);
             i > nMultiStepVals; i--)
        {
            m_solVector[i] = m_solVector[i - 1];
            m_setflag[i]   = m_setflag[i - 1];
        }

        for (size_t i = (size - 1); i > nMultiStepVals + nMultiStepImpDerivs;
             i--)
        {
            m_solVector[i] = m_solVector[i - 1];
            m_setflag[i]   = m_setflag[i - 1];
        }
    }

private:
    const TimeIntegrationAlgorithmGLM *m_schemeAlgorithm;

    TripleArray m_solVector;
    Array<OneD, NekDouble> m_t;
    Array<OneD, bool> m_setflag;

}; // end class TimeIntegrationSolutionGLM

} // end of namespace LibUtilities
} // end of namespace Nektar

#endif
