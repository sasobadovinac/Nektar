///////////////////////////////////////////////////////////////////////////////
//
// File: FilterHistoryPoints.h
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
// Description: Outputs values at specific points during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERHISTORYPOINTS_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERHISTORYPOINTS_H

#include "LibUtilities/BasicConst/NektarUnivTypeDefs.hpp"
#include "LibUtilities/BasicUtils/SharedArray.hpp"
#include <SolverUtils/Filters/Filter.h>

namespace Nektar
{
namespace SolverUtils
{

class FilterHistoryPoints : public Filter
{
public:
    friend class MemoryManager<FilterHistoryPoints>;

    /// Creates an instance of this class
    static FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation,
        const std::map<std::string, std::string> &pParams)
    {
        FilterSharedPtr p =
            MemoryManager<FilterHistoryPoints>::AllocateSharedPtr(
                pSession, pEquation, pParams);
        return p;
    }

    /// Name of the class
    static std::string className;

    SOLVER_UTILS_EXPORT FilterHistoryPoints(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation,
        const ParamMap &pParams);
    SOLVER_UTILS_EXPORT ~FilterHistoryPoints();

protected:
    SOLVER_UTILS_EXPORT virtual void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;
    SOLVER_UTILS_EXPORT virtual void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;
    SOLVER_UTILS_EXPORT virtual void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;
    SOLVER_UTILS_EXPORT virtual bool v_IsTimeDependent() override;
    bool GetPoint(Array<OneD, NekDouble> gloCoord, int I);
    void WriteData(const int &rank, const Array<OneD, NekDouble> &data,
                   const int &numFields, const NekDouble &time);

    Array<OneD, Array<OneD, const NekDouble>> m_historyPoints =
        Array<OneD, Array<OneD, const NekDouble>>(0);
    size_t m_historyPointsSize = 0;
    unsigned int m_index       = 0;
    unsigned int m_outputFrequency;
    /// plane to take history point from if using a homogeneous1D expansion
    int m_outputPlane;
    std::vector<int> m_planeIDs;
    bool m_isHomogeneous1D;
    bool m_waveSpace;
    std::string m_outputFile;
    std::ofstream m_outputStream;
    std::stringstream m_historyPointStream;
    // List of history points that are local to this process
    // Content of tuple:
    // [0] = global coordinates
    // [1] = local coordinates
    // [2] = global index of history point
    // [3] = expansion index of history point
    std::list<std::tuple<Array<OneD, const NekDouble>,
                         Array<OneD, const NekDouble>, int, int>>
        m_historyList;
    std::map<LibUtilities::PtsType, Array<OneD, NekDouble>> m_pointDatMap;
    std::map<LibUtilities::PtsType, Array<OneD, int>> m_pointNumMap;
    unsigned int m_outputIndex = 0;
    bool m_outputOneFile;
    bool m_adaptive;
};

} // namespace SolverUtils
} // namespace Nektar

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERCHECKPOINT_H */
