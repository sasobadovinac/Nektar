///////////////////////////////////////////////////////////////////////////////
//
// File FilterGJPSmoothing.h
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
// Description: Average solution fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERGJPSMOOTHING_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERGJPSMOOTHING_H

#include <SolverUtils/Filters/Filter.h>

namespace Nektar
{
namespace SolverUtils
{

typedef std::tuple<int, int, NekDouble> TraceToCoeffMap;
    
class FilterGJPSmoothing : public Filter
{
public:
    friend class MemoryManager<FilterGJPSmoothing>;

    /// Creates an instance of this class
    static FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem>      &pEquation,
        const std::map<std::string, std::string>   &pParams)
    {
        FilterSharedPtr p = MemoryManager<FilterGJPSmoothing>
                            ::AllocateSharedPtr(pSession, pEquation, pParams);
        return p;
    }

    ///Name of the class
    static std::string className;

    SOLVER_UTILS_EXPORT FilterGJPSmoothing(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem>      &pEquation,
        const ParamMap &pParams);
    SOLVER_UTILS_EXPORT virtual ~FilterGJPSmoothing();

protected:
    SOLVER_UTILS_EXPORT virtual void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);

    SOLVER_UTILS_EXPORT virtual bool v_IsTimeDependent();
private:
    int m_smoothingFrequency; 
    int m_shapeDim; 
    int m_index;
    /// Scale factor for phys values along trace involving the lcoal
    /// normals and tangenet geometric factors on Fwd Trace
    Array<OneD, NekDouble> m_scalFwd;
    /// Scale factor for phys values along trace involving the lcoal
    /// normals and tangenet geometric factors on Bwd Trace
    Array<OneD, NekDouble> m_scalBwd;
    /// Array for every coefficient on trace of a mapping form trace
    /// coefficients to elemental coefficients with a scaling factor
    /// of the derivative of the normal basis along that trace
    Array<OneD, std::pair<NekDouble, std::set<int> > >m_fwdTraceToCoeffMap; 
    /// Array for every coefficient on trace of a mapping form trace
    /// coefficients to elemental coefficients with a scaling factor
    /// of the derivative of the normal basis along that trace
    Array<OneD, std::pair<NekDouble, std::set<int> > >m_bwdTraceToCoeffMap;

    // Link to the trace normals 
    Array<OneD, const Array<OneD, NekDouble> > m_traceNormals; 
};
}
}

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERGJPSMOOTHING_H */
