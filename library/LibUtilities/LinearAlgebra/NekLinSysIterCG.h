///////////////////////////////////////////////////////////////////////////////
//
// File  NekLinSysIterCG.h
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
// Description: NekLinSysIterCG header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_ITERAT_CG_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_LINSYS_ITERAT_CG_H

#include <LibUtilities/LinearAlgebra/NekLinSysIter.h>
#include <boost/circular_buffer.hpp>
namespace Nektar
{
namespace LibUtilities
{
/// A global linear system.
class NekLinSysIterCG;

typedef std::shared_ptr<NekLinSysIterCG> NekLinSysIterCGSharedPtr;

class NekLinSysIterCG : public NekLinSysIter
{
public:
    /// Support creation through MemoryManager.
    friend class MemoryManager<NekLinSysIterCG>;

    LIB_UTILITIES_EXPORT static NekLinSysIterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm, const int nDimen,
        const NekSysKey &pKey)
    {
        NekLinSysIterCGSharedPtr p =
            MemoryManager<NekLinSysIterCG>::AllocateSharedPtr(pSession, vComm,
                nDimen, pKey);
        p->InitObject();
        return p;
    }
    static std::string className;
    /// Constructor for full direct matrix solve.
    LIB_UTILITIES_EXPORT NekLinSysIterCG(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm, const int nDimen,
        const NekSysKey &pKey);
    LIB_UTILITIES_EXPORT ~NekLinSysIterCG();

protected:
    int m_successiveRHS;
    /// Whether to apply projection technique
    bool m_useProjection = false;

    /// Storage for solutions to previous linear problems
    boost::circular_buffer<Array<OneD, NekDouble>> m_prevLinSol;

    /// Total counter of previous solutions
    int m_numPrevSols = 0;

    virtual void v_InitObject();

    virtual int v_SolveSystem(const int nGlobal,
                              const Array<OneD, const NekDouble> &pInput,
                              Array<OneD, NekDouble> &pOutput, const int nDir,
                              const NekDouble tol, const NekDouble factor);

private:
    /// A-conjugate projection technique
    void DoAconjugateProjection(const int pNumRows,
                                const Array<OneD, const NekDouble> &pInput,
                                Array<OneD, NekDouble> &pOutput,
                                const int pNumDir);

    /// Actual iterative solve
    void DoConjugateGradient(const int pNumRows,
                             const Array<OneD, const NekDouble> &pInput,
                             Array<OneD, NekDouble> &pOutput,
                             const int pNumDir);

    void UpdateKnownSolutions(const int pGlobalBndDofs,
                              const Array<OneD, const NekDouble> &pSolution,
                              const int pNumDirBndDofs);

    NekDouble CalculateAnorm(const int nGlobal,
                             const Array<OneD, const NekDouble> &in,
                             const int nDir);
};
} // namespace LibUtilities
} // namespace Nektar

#endif
