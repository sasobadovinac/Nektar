///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysIterativeFull.cpp
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
// Description: GlobalLinSysIterativeFull definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/GlobalLinSysIterativeFull.h>
#include <map>

using namespace std;

namespace Nektar
{
namespace MultiRegions
{
/**
 * @class GlobalLinSysIterativeCG
 *
 * This class implements a conjugate gradient matrix solver.
 * Preconditioning is implemented using a Jacobi (diagonal)
 * preconditioner.
 */

/**
 * Registers the class with the Factory.
 */
string GlobalLinSysIterativeFull::className =
    GetGlobalLinSysFactory().RegisterCreatorFunction(
        "IterativeFull", GlobalLinSysIterativeFull::create,
        "Iterative solver for full matrix system.");

/**
 * Constructor for full direct matrix solve.
 * @param   pKey        Key specifying matrix to solve.
 * @param   pExp        Shared pointer to expansion list for applying
 *                      matrix evaluations.
 * @param   pLocToGloMap Local to global mapping.
 */
GlobalLinSysIterativeFull::GlobalLinSysIterativeFull(
    const GlobalLinSysKey &pKey, const std::weak_ptr<ExpList> &pExp,
    const std::shared_ptr<AssemblyMap> &pLocToGloMap)
    : GlobalLinSys(pKey, pExp, pLocToGloMap),
      GlobalLinSysIterative(pKey, pExp, pLocToGloMap)
{
    ASSERTL1(m_linSysKey.GetGlobalSysSolnType() == eIterativeFull,
             "This routine should only be used when using an Iterative "
             "conjugate gradient matrix solve.");
}

/**
 *
 */
GlobalLinSysIterativeFull::~GlobalLinSysIterativeFull()
{
}

/**
 * Solve a global linear system with Dirichlet forcing using a
 * conjugate gradient method. This routine performs handling of the
 * Dirichlet forcing terms and wraps the underlying iterative solver
 * used for the remaining degrees of freedom.
 *
 * Consider solving for \f$x\f$, the matrix system \f$Ax=b\f$, where
 * \f$b\f$ is known. To enforce the Dirichlet terms we instead solve
 * \f[A(x-x_0) = b - Ax_0 \f]
 * where \f$x_0\f$ is the Dirichlet forcing.
 *
 * @param           pInput      RHS of linear system, \f$b\f$.
 * @param           pOutput     On input, values of dirichlet degrees
 *                              of freedom with initial guess on other values.
 *                              On output, the solution \f$x\f$.
 * @param           pLocToGloMap    Local to global mapping.
 * @param           pDirForcing Precalculated Dirichlet forcing.
 */
void GlobalLinSysIterativeFull::v_Solve(
    const Array<OneD, const NekDouble> &pLocInput,
    Array<OneD, NekDouble> &pLocOutput,
    const AssemblyMapSharedPtr &pLocToGloMap,
    const Array<OneD, const NekDouble> &pDirForcing)
{
    std::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
    bool vCG                                       = false;
    m_locToGloMap                                  = pLocToGloMap;

    if (std::dynamic_pointer_cast<AssemblyMapCG>(pLocToGloMap))
    {
        vCG = true;
    }
    else if (std::dynamic_pointer_cast<AssemblyMapDG>(pLocToGloMap))
    {
        vCG = false;
    }
    else
    {
        NEKERROR(ErrorUtil::efatal, "Unknown map type");
    }

    bool dirForcCalculated = (bool)pDirForcing.size();
    int nDirDofs           = pLocToGloMap->GetNumGlobalDirBndCoeffs();
    int nGlobDofs          = pLocToGloMap->GetNumGlobalCoeffs();
    int nLocDofs           = pLocToGloMap->GetNumLocalCoeffs();

    int nDirTotal = nDirDofs;
    expList->GetComm()->GetRowComm()->AllReduce(nDirTotal,
                                                LibUtilities::ReduceSum);

    Array<OneD, NekDouble> tmp(nLocDofs);
    Array<OneD, NekDouble> tmp1(nLocDofs);

    if (nDirTotal)
    {
        // calculate the Dirichlet forcing
        if (dirForcCalculated)
        {
            Vmath::Vsub(nLocDofs, pLocInput, 1, pDirForcing, 1, tmp1, 1);
        }
        else
        {
            // Calculate the dirichlet forcing B_b (== X_b) and
            // substract it from the rhs
            expList->GeneralMatrixOp(m_linSysKey, pLocOutput, tmp);

            // Iterate over all the elements computing Robin BCs where
            // necessary
            for (auto &r : m_robinBCInfo) // add robin mass matrix
            {
                RobinBCInfoSharedPtr rBC;
                Array<OneD, NekDouble> tmploc;

                int n      = r.first;
                int offset = expList->GetCoeff_Offset(n);

                LocalRegions::ExpansionSharedPtr vExp = expList->GetExp(n);
                // add local matrix contribution
                for (rBC = r.second; rBC; rBC = rBC->next)
                {
                    vExp->AddRobinTraceContribution(
                        rBC->m_robinID, rBC->m_robinPrimitiveCoeffs,
                        pLocOutput + offset, tmploc = tmp + offset);
                }
            }

            Vmath::Vsub(nLocDofs, pLocInput, 1, tmp, 1, tmp1, 1);
        }
        if (vCG)
        {
            // solve for perturbation from initial guess in pOutput
            SolveLinearSystem(nGlobDofs, tmp1, tmp, pLocToGloMap, nDirDofs);

            // Add back initial condition
            Vmath::Vadd(nLocDofs, tmp, 1, pLocOutput, 1, pLocOutput, 1);
        }
        else
        {
            ASSERTL0(false, "Need DG solve if using Dir BCs");
        }
    }
    else
    {
        SolveLinearSystem(nGlobDofs, pLocInput, pLocOutput, pLocToGloMap,
                          nDirDofs);
    }
}

/**
 *
 */
void GlobalLinSysIterativeFull::v_DoMatrixMultiply(
    const Array<OneD, NekDouble> &pInput, Array<OneD, NekDouble> &pOutput)
{
    bool isLocal = m_linsol->IsLocal();

    std::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();

    AssemblyMapSharedPtr asmMap = m_locToGloMap.lock();

    int ncoeffs = expList->GetNcoeffs();
    Array<OneD, NekDouble> InputLoc, OutputLoc;

    if (isLocal)
    {
        InputLoc  = pInput;
        OutputLoc = pOutput;

        // Perform matrix-vector operation A*d_i
        expList->GeneralMatrixOp(m_linSysKey, InputLoc, OutputLoc);
    }
    else
    {
        InputLoc  = Array<OneD, NekDouble>(ncoeffs);
        OutputLoc = Array<OneD, NekDouble>(ncoeffs);

        asmMap->GlobalToLocal(pInput, InputLoc);
        // Perform matrix-vector operation A*d_i
        expList->GeneralMatrixOp(m_linSysKey, InputLoc, OutputLoc);
    }

    // Apply robin boundary conditions to the solution.
    for (auto &r : m_robinBCInfo) // add robin mass matrix
    {
        RobinBCInfoSharedPtr rBC;
        Array<OneD, NekDouble> tmp;

        int n = r.first;

        int offset                            = expList->GetCoeff_Offset(n);
        LocalRegions::ExpansionSharedPtr vExp = expList->GetExp(n);

        // add local matrix contribution
        for (rBC = r.second; rBC; rBC = rBC->next)
        {
            vExp->AddRobinTraceContribution(
                rBC->m_robinID, rBC->m_robinPrimitiveCoeffs, InputLoc + offset,
                tmp = OutputLoc + offset);
        }
    }

    // put back in global coeffs
    if (isLocal == false)
    {
        asmMap->Assemble(OutputLoc, pOutput);
    }
}

/**
 *
 */
void GlobalLinSysIterativeFull::v_UniqueMap()
{
    m_map = m_locToGloMap.lock()->GetGlobalToUniversalMapUnique();
}

/**
 *
 */
void GlobalLinSysIterativeFull::v_SolveLinearSystem(
    const int nGlobal, const Array<OneD, const NekDouble> &pInput,
    Array<OneD, NekDouble> &pOutput, const AssemblyMapSharedPtr &plocToGloMap,
    const int nDir)
{

    if (!m_linsol)
    {
        LibUtilities::CommSharedPtr vComm =
            m_expList.lock()->GetComm()->GetRowComm();
        LibUtilities::SessionReaderSharedPtr pSession =
            m_expList.lock()->GetSession();

        // Check such a module exists for this equation.
        ASSERTL0(LibUtilities::GetNekLinSysIterFactory().ModuleExists(
                     m_linSysIterSolver),
                 "NekLinSysIter '" + m_linSysIterSolver +
                     "' is not defined.\n");
        m_linsol = LibUtilities::GetNekLinSysIterFactory().CreateInstance(
            m_linSysIterSolver, pSession, vComm, nGlobal - nDir,
            LibUtilities::NekSysKey());

        m_linsol->SetSysOperators(m_NekSysOp);
        v_UniqueMap();
        m_linsol->setUniversalUniqueMap(m_map);
    }

    if (!m_precon)
    {
        m_precon = CreatePrecon(plocToGloMap);
        m_precon->BuildPreconditioner();
    }

    m_linsol->setRhsMagnitude(m_rhs_magnitude);

    if (m_useProjection)
    {
        Array<OneD, NekDouble> gloIn(nGlobal);
        Array<OneD, NekDouble> gloOut(nGlobal, 0.0);
        plocToGloMap->Assemble(pInput, gloIn);
        DoProjection(nGlobal, gloIn, gloOut, nDir, m_tolerance, m_isAconjugate);
        plocToGloMap->GlobalToLocal(gloOut, pOutput);
    }
    else
    {
        if (m_linsol->IsLocal())
        {
            int nLocDofs = plocToGloMap->GetNumLocalCoeffs();
            Vmath::Zero(nLocDofs, pOutput, 1);
            m_linsol->SolveSystem(nLocDofs, pInput, pOutput, nDir, m_tolerance);
        }
        else
        {
            Array<OneD, NekDouble> gloIn(nGlobal);
            Array<OneD, NekDouble> gloOut(nGlobal, 0.0);
            plocToGloMap->Assemble(pInput, gloIn);
            m_linsol->SolveSystem(nGlobal, gloIn, gloOut, nDir, m_tolerance);
            plocToGloMap->GlobalToLocal(gloOut, pOutput);
        }
    }
}
} // namespace MultiRegions
} // namespace Nektar
