///////////////////////////////////////////////////////////////////////////////
//
// File: FilterMaxMinFields.h
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
// Description: Maximun/minimum solution fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERMAXMINFIELDS_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERMAXMINFIELDS_H

#include <SolverUtils/Filters/FilterFieldConvert.h>

namespace Nektar
{
namespace SolverUtils
{

enum ProblemType
{
    eCompressible,
    eIncompressible,
    eOthers
};

class FilterMaxMinFields : public FilterFieldConvert
{
public:
    friend class MemoryManager<FilterMaxMinFields>;

    /// Creates an instance of this class
    static FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation,
        const std::map<std::string, std::string> &pParams)
    {
        FilterSharedPtr p =
            MemoryManager<FilterMaxMinFields>::AllocateSharedPtr(
                pSession, pEquation, pParams);
        return p;
    }

    /// Name of the class
    static std::string className;

    SOLVER_UTILS_EXPORT FilterMaxMinFields(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation,
        const ParamMap &pParams);
    SOLVER_UTILS_EXPORT virtual ~FilterMaxMinFields();

protected:
    bool m_isMax;
    ProblemType m_problemType;
    std::vector<Array<OneD, NekDouble>> m_curFieldsPhys;
    std::vector<Array<OneD, NekDouble>> m_outFieldsPhys;

    virtual void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;
    virtual void v_FillVariablesName(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
        override
    {
        FilterFieldConvert::v_FillVariablesName(pFields);
    }
    virtual void v_ProcessSample(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
        const NekDouble &time) override;
    virtual void v_PrepareOutput(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) override;
    virtual NekDouble v_GetScale() override;
    virtual std::string v_GetFileSuffix() override
    {
        if (m_isMax)
        {
            return "_max";
        }
        else
        {
            return "_min";
        }
    }

private:
    bool m_initialized;
};
} // namespace SolverUtils
} // namespace Nektar

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERCHECKPOINT_H */
