///////////////////////////////////////////////////////////////////////////////
//
// File: PreconCfs.h
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
// Description: PreconCfs header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_PRECONCFS
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_PRECONCFS

#include <CompressibleFlowSolver/Preconditioner/PreconCfsOp.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
enum PrecType
{
    eNull, ///< No Solution type specified
    eDiagonal,
    eSparse,
};

//  Forward declaration
class PreconCfs;

/// Declaration of the boundary condition factory
typedef LibUtilities::NekFactory<
    std::string, PreconCfs, const Array<OneD, MultiRegions::ExpListSharedPtr> &,
    const LibUtilities::SessionReaderSharedPtr &,
    const LibUtilities::CommSharedPtr &>
    PreconCfsFactory;

/// Declaration of the boundary condition factory singleton
PreconCfsFactory &GetPreconCfsFactory();

class PreconCfs
{
public:
    PreconCfs(const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
              const LibUtilities::SessionReaderSharedPtr &pSession,
              const LibUtilities::CommSharedPtr &vComm);

    virtual ~PreconCfs()
    {
    }

    inline void InitObject();

    void DoPreconCfs(const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                     const Array<OneD, NekDouble> &pInput,
                     Array<OneD, NekDouble> &pOutput, const bool &flag);

    inline void BuildPreconCfs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, const Array<OneD, NekDouble>> &intmp,
        const NekDouble time, const NekDouble lambda);

    bool UpdatePreconMatCheck(const Array<OneD, const NekDouble> &res,
                              const NekDouble dtLambda);

    inline void SetOperators(const NekPreconCfsOperators &in)
    {
        m_operator = in;
    }

protected:
    LibUtilities::CommSharedPtr m_Comm;
    bool m_verbose;
    int m_spacedim;

    NekPreconCfsOperators m_operator;

    int m_PreconMatFreezNumb;
    int m_PreconTimesCounter;

    NekDouble m_DtLambdaPreconMat = -1.0;
    NekDouble m_BndEvaluateTime;

    bool m_CalcPreconMatFlag = false;

    virtual void v_InitObject() = 0;

    virtual void v_DoPreconCfs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, NekDouble> &pInput, Array<OneD, NekDouble> &pOutput,
        const bool &flag) = 0;

    virtual void v_BuildPreconCfs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, const Array<OneD, NekDouble>> &intmp,
        const NekDouble time, const NekDouble lambda) = 0;

    virtual bool v_UpdatePreconMatCheck(const Array<OneD, const NekDouble> &res,
                                        const NekDouble dtLambda) = 0;
};
typedef std::shared_ptr<PreconCfs> PreconCfsSharedPtr;

/**
 *
 */
inline void PreconCfs::InitObject()
{
    v_InitObject();
}

/**
 *
 */
inline void PreconCfs::DoPreconCfs(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, NekDouble> &pInput, Array<OneD, NekDouble> &pOutput,
    const bool &flag)
{
    ASSERTL0(pInput.size() == pOutput.size(),
             "In and Out not the same size in DoPreconCfs");
    v_DoPreconCfs(pFields, pInput, pOutput, flag);
    m_PreconTimesCounter++;
}

/**
 *
 */
inline void PreconCfs::BuildPreconCfs(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, const Array<OneD, NekDouble>> &intmp,
    const NekDouble time, const NekDouble lambda)
{
    v_BuildPreconCfs(pFields, intmp, time, lambda);
}

/**
 *
 */
inline bool PreconCfs::UpdatePreconMatCheck(
    const Array<OneD, const NekDouble> &res, const NekDouble dtLambda)
{
    return v_UpdatePreconMatCheck(res, dtLambda);
}
} // namespace Nektar

#endif
