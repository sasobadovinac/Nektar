///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSys.cpp
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
// Description: GlobalLinSys definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysSaena.h>
#include <MultiRegions/Preconditioner.h>
#include <LibUtilities/Communication/CommMpi.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysSaena
         *
         * Solves a linear system using Saena.
         */
        GlobalLinSysSaena::GlobalLinSysSaena(
            const GlobalLinSysKey                &pKey,
            const std::weak_ptr<ExpList>         &pExp,
            const std::shared_ptr<AssemblyMap>   &pLocToGloMap)
            : GlobalLinSys(pKey, pExp, pLocToGloMap)
        {
        }

        /**
         * @brief Clean up Saena objects.
         *
         * Note that if SessionReader::Finalize is called before the end of the
         * program, Saena may have been finalized already, at which point we
         * cannot deallocate our objects. If that's the case we do nothing and
         * let the kernel clear up after us.
         */
        GlobalLinSysSaena::~GlobalLinSysSaena()
        {
        }

        /**
         * @brief Solve linear system using Saena.
         *
         * The general strategy being a Saena solve is to:
         *
         * - Copy values into the Saena vector #m_b
         * - Solve the system #m_ksp and place result into #m_x.
         * - Scatter results back into #m_locVec using #m_ctx scatter object.
         * - Copy from #m_locVec to output array #pOutput.
         */
        void GlobalLinSysSaena::v_SolveLinearSystem(
            const int                          pNumRows,
            const Array<OneD,const NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const AssemblyMapSharedPtr        &locToGloMap,
            const int                          pNumDir)
        {
            boost::ignore_unused(locToGloMap);

            // @TODO: shouldn't need to but we require a new RHS vector every
            // time this is called.
            saena::vector m_rhs;
            m_rhs.set_comm(m_comm);

            const int nHomDofs = pNumRows - pNumDir;

            m_rhs.set(&m_reorderedMap[0], &pInput[pNumDir], nHomDofs);
            m_rhs.assemble();
            m_amg.set_rhs(m_rhs);

            // Temporary solution storage?
            NekDouble *sol = nullptr;

            // Solve with pCG method
            m_amg.solve_pCG(sol, &m_opts);

            Vmath::Vcopy(nHomDofs, sol, 1, &pOutput[pNumDir], 1);

            if(sol != nullptr)
            {
                free(sol);
                sol = nullptr;
            }
        }

        /**
         * @brief Calculate a reordering of universal IDs for Saena.
         *
         * Saena requires a unique, contiguous index of all global and universal
         * degrees of freedom which represents its position inside the
         * matrix. Presently Gs does not guarantee this, so this routine
         * constructs a new universal mapping.
         *
         * @param glo2uniMap    Global to universal map
         * @param glo2unique    Global to unique map
         * @param pLocToGloMap  Assembly map for this system
         */
        void GlobalLinSysSaena::CalculateReordering(
            const Array<OneD, const int> &glo2uniMap,
            const Array<OneD, const int> &glo2unique,
            const AssemblyMapSharedPtr   &pLocToGloMap)
        {
            LibUtilities::CommSharedPtr vComm
                = m_expList.lock()->GetSession()->GetComm();

            const int nDirDofs = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            const int nHomDofs = glo2uniMap.size() - nDirDofs;
            const int nProc    = vComm->GetSize();
            const int rank     = vComm->GetRank();

            int n, cnt;

            // Count number of unique degrees of freedom on each process.
            m_nLocal = Vmath::Vsum(nHomDofs, glo2unique + nDirDofs, 1);
            m_reorderedMap.resize(nHomDofs);

            // Reduce coefficient counts across all processors.
            Array<OneD, int> localCounts(nProc, 0), localOffset(nProc, 0);
            localCounts[rank] = nHomDofs;
            vComm->AllReduce(localCounts, LibUtilities::ReduceSum);

            for (n = 1; n < nProc; ++n)
            {
                localOffset[n] = localOffset[n-1] + localCounts[n-1];
            }

            int totHomDofs = Vmath::Vsum(nProc, localCounts, 1);
            vector<unsigned int> allUniIds(totHomDofs, 0);

            // Assemble list of universal IDs
            for (n = 0; n < nHomDofs; ++n)
            {
                int gid = n + nDirDofs;
                allUniIds[n + localOffset[rank]] = glo2uniMap[gid];
            }

            // Reduce this across processors so that each process has a list of
            // all universal IDs.
            vComm->AllReduce(allUniIds, LibUtilities::ReduceSum);
            std::sort(allUniIds.begin(), allUniIds.end());
            map<int,int> uniIdReorder;

            // Renumber starting from 0.
            for (cnt = n = 0; n < allUniIds.size(); ++n)
            {
                if (uniIdReorder.count(allUniIds[n]) > 0)
                {
                    continue;
                }

                uniIdReorder[allUniIds[n]] = cnt++;
            }

            // Populate reordering map.
            for (n = 0; n < nHomDofs; ++n)
            {
                int gid = n + nDirDofs;
                int uniId = glo2uniMap[gid];
                ASSERTL0(uniIdReorder.count(uniId) > 0, "Error in ordering");
                m_reorderedMap[n] = uniIdReorder[uniId];
            }

            m_bdydof = nDirDofs;
        }

        /**
         * @brief Construct Saena matrix and vector handles.
         *
         * @todo Preallocation should be done at this point, since presently
         *       matrix allocation takes a significant amount of time.
         *
         * @param nGlobal  Number of global degrees of freedom in the system (on
         *                 this processor)
         * @param nDir     Number of Dirichlet degrees of freedom (on this
         *                 processor).
         */
        void GlobalLinSysSaena::SetUpMatVec()
        {
            LibUtilities::CommSharedPtr comm =
                m_expList.lock()->GetSession()->GetComm();
            auto mpiComm = std::dynamic_pointer_cast<
                LibUtilities::CommMpi>(comm);

            m_comm = mpiComm->GetComm();
            m_matrix.set_comm(m_comm);
            m_matrix.add_duplicates(true);
            m_rhs.set_comm(m_comm);

            int nummodes = m_expList.lock()->GetFieldDefinitions()[0]->m_numModes[0];
            int p_order  = nummodes - 1;
            int prodim   = m_expList.lock()->GetCoordim(0);

            m_matrix.set_p_order(p_order);
            m_matrix.set_prodim(prodim);

            // set p_coarsen levels computation. subtract by a constant.
            vector<int> order_dif;
            for(int i = 0; i < p_order - 1; ++i)
            {
                order_dif.emplace_back(1);
            }

            // set number of multigrid levels
            int max_h_level = 1; // h-multigrid levels
            m_amg.set_multigrid_max_level(
                static_cast<int>(order_dif.size()) + max_h_level);

            m_amg.set_scale(m_scale);
            m_amg.set_matrix(
                &m_matrix, &m_opts, m_l2g, m_reorderedMap, m_bdydof, order_dif);
        }

        /**
         * @brief Set up KSP solver object.
         *
         * This is reasonably generic setup -- most solver types can be changed
         * using the ? file.
         *
         * @param tolerance  Residual tolerance to converge to.
         */
        void GlobalLinSysSaena::SetUpSolver(NekDouble tolerance)
        {
            m_scale = false;
            m_opts.set_relative_tolerance(tolerance);
            // m_opts.set_dynamic_levels(false);
            // m_opts.set_max_lev(5);
            // m_opts.set_vcycle_num(400);
            // m_opts.set_smoother("chebyshev"); // chebyshev, jacobi
            // m_opts.set_preSmooth(3);
            // m_opts.set_postSmooth(3);
        }
    }
}
