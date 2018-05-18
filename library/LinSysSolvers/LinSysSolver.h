///////////////////////////////////////////////////////////////////////////////
//
// File: LinSysSolver.h
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
// Description: GlobalLinSys header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LINSYSSOLVERS_LINSYSSOLVER
#define NEKTAR_LINSYSSOLVERS_LINSYSSOLVER

#include <functional>

namespace Nektar {
namespace MultiRegions {

class LinSysOperator
{
public:
    LinSysOperator();

    virtual void Apply(const Array<OneD, const NekDouble> &in
                             Array<OneD,       NekDouble> &out) = 0;

    virtual MultiRegions::AssemblyMapSharedPtr GetAssemblyMap()
    {
        ASSERTL0(false, "nope");
    }

    virtual vector<DNekScalMatSharedPtr> const &GetBlocks()
    {
        ASSERTL0(false, "nope");
    }
};

class LinSysOperatorBlock : public LinSysOperator
{
public:
    LinSysOperatorBlock(MultiRegions::AssemblyMapSharedPtr asmMap,
                        std::vector<DNekScalMatSharedPtr> &locBlocks)
        : m_assemblyMap(asmMap), m_locBlocks(locBlocks)
    {
    }

    virtual void Apply(const Array<OneD, const NekDouble> &in
                             Array<OneD,       NekDouble> &out)
    {
        int i, cnt;
        const int nLocal = m_locToGloMap->GetNumLocalBndCoeffs();
        Array<OneD, NekDouble> tmpout(nLocal), tmpin(nLocal);

        m_assemblyMap->GlobalToLocalBnd(in, tmpin);

        for (i = cnt = 0; i < m_locBlocks.size(); ++i)
        {
            const int rows = m_locBlocks[i]->GetNumRows();
            Blas::Dgemv('N', rows, rows, m_locBlocks[i]->GetScale(),
                        m_locBlocks[i]->GetRawPtr(), rows,
                        tmpin.get() + cnt, 1, 0.0, tmpout.get() + cnt, 1);
        }

        m_assemblyMap->AssembleBnd(tmpout, out);
    }

    virtual vector<DNekScalMatSharedPtr> const &GetBlocks()
    {
        return m_locBlocks;
    }

    virtual MultiRegions::AssemblyMapSharedPtr GetAssemblyMap()
    {
        return m_assemblyMap;
    }

protected:
    MultiRegions::AssemblyMapSharedPtr m_assemblyMap;
    std::vector<DNekScalMatSharedPtr> m_locBlocks;
};

class LinSysOperatorFunction : public LinSysOperator
{
    LinSysOperatorFunction(ApplyFunction func) : m_func(func)
    {
    }

    virtual void Apply(
        const Array<OneD, const NekDouble> &in
              Array<OneD,       NekDouble> &out)
    {
        m_func(in, out);
    }

protected:
    ApplyFunction m_func;
};

typedef std::function<void (const Array<OneD, const NekDouble>&,
                            Array<OneD, NekDouble>&)> ApplyFunction;

class LinSysSolver
{
public:
    void Solve(
        const Array<OneD,const NekDouble> &pInput,
              Array<OneD,      NekDouble> &pOutput)
    {
        v_Solve(pInput, pOutput);
    }

protected:
    LinSysSolver(OperatorSharedPtr matOp, OperatorSharedPtr preconOp) :
        m_matOp(matOp), m_preconOp(preconOp)
    {
    }

    virtual void v_Solve(
        const Array<OneD,const NekDouble> &pInput,
              Array<OneD,      NekDouble> &pOutput) = 0;

    OperatorSharedPtr m_matOp;
    OperatorSharedPtr m_preconOp;
};

class LinSysSolverLAPACK : public LinSysSolver
{
public:
    LinSysSolverLAPACK(OperatorSharedPtr matOp, OperatorSharedPtr preconOp)
        : LinSysSolver(matOp, preconOp)
    {
        auto          asmMap     = matOp->GetAssemblyMap();
        const int     totDofs    = asmMap->GetNumGlobalCoeffs();
        const int     nDirDofs   = asmMap->GetNumGlobalDirBndCoeffs();
        unsigned int  rows       = totDofs - nDirDofs;
        MatrixStorage matStorage = ePOSITIVE_DEFINITE_SYMMETRIC; // hardcoded for now
        DNekMatSharedPtr Gmat    = MemoryManager<DNekMat>
            ::AllocateSharedPtr(rows, rows, 0.0, matStorage);

        auto blocks = matOp->GetBlocks();
        auto l2gMap = asmMap->GetLocalToGlobalMap();
        auto l2gSign = asmMap->GetLocalToGlobalSign();

        int cnt = 0;

        // Assemble global matrix from blocks.
        for (auto &locMat : blocks)
        {
            int locRows = mat->GetRows();

            for (i = 0; i < locRows; ++i)
            {
                int gid1 = l2gMap[cnt + i] - nDirDofs, sign1 = l2gSign[cnt + i];

                if (gid1 < 0)
                {
                    continue;
                }

                for (j = 0; j < locRows; ++j)
                {
                    int gid2 = l2gMap[cnt + j] - nDirDofs;
                    int sign2 = l2gSign[cnt + j];

                    if (gid2 < 0)
                    {
                        continue;
                    }

                    // When global matrix is symmetric, only add the value for
                    // the upper triangular part in order to avoid entries to be
                    // entered twice
                    if (matStorage == eFULL || gid2 >= gid1)
                    {
                        Gmat->SetValue(
                            gid1, gid2, Gmat->GetValue(gid1, gid2)
                            + sign1*sign2*(*locMat)(i,j));
                    }
                }
            }

            cnt += locRows;
        }

        m_linSys = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat, eWrapper);
    }

    virtual void v_Solve(
        const Array<OneD,const NekDouble> &pInput,
              Array<OneD,      NekDouble> &pOutput);
    {
        int nHomDofs = matOp->GetAssemblyMap()->Get();
        DNekVec Vin (nHomDofs, pInput), Vout(nHomDofs, tmp, eWrapper);

        m_linSys->Solve(Vin, Vout);
    }

private:
    DNekLinSysSharedPtr m_linSys;
};

}
}

#endif
