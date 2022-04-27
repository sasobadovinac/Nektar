///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysSaenaFull.cpp
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
// Description: GlobalLinSysSaenaFull definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysSaenaFull.h>
#include <MultiRegions/ExpList.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysSaenaFull
         */

        /**
         * Registers the class with the Factory.
         */
        string GlobalLinSysSaenaFull::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "SaenaFull",
                    GlobalLinSysSaenaFull::create,
                    "Saena Full Matrix.");


        /// Constructor for full direct matrix solve.
        GlobalLinSysSaenaFull::GlobalLinSysSaenaFull(
            const GlobalLinSysKey                &pLinSysKey,
            const std::weak_ptr<ExpList>         &pExp,
            const std::shared_ptr<AssemblyMap>   &pLocToGloMap)
            : GlobalLinSys     (pLinSysKey, pExp, pLocToGloMap),
              GlobalLinSysSaena(pLinSysKey, pExp, pLocToGloMap)
        {
            // SET UP VECTORS AND MATRIX
            SetUpMatVec();

            int rank = 0, nprocs = 0;
            MPI_Comm_size(m_comm, &nprocs);
            MPI_Comm_rank(m_comm, &rank);

            auto tbegin = clock();

            const int nDirDofs = pLocToGloMap->GetNumGlobalDirBndCoeffs();

            int i, j, n, cnt, gid1, gid2, loc_lda;
            NekDouble sign1, sign2, value;
            DNekScalMatSharedPtr loc_mat;

            // CALCULATE REORDERING MAPPING
            CalculateReordering(pLocToGloMap->GetGlobalToUniversalMap(),
                                pLocToGloMap->GetGlobalToUniversalMapUnique(),
                                pLocToGloMap);

            // STORE MESH INFO TO BE PASSED TO SAENA
//            int total_elm = this->GetExp()->size();
            auto ExpTmp = m_expList.lock()->GetExp();
            int total_elm = ExpTmp->size();
//            std::cout << total_elm << "\n";

            int counter = 0;
            vector<int> dof_elems;
            for (i = 0; i < total_elm; ++i){
//                std::cout << ExpTmp->at(i)->GetNcoeffs() << std::endl;
                for (j = 0; j < ExpTmp->at(i)->GetNcoeffs(); ++j){
//                    printf("%i\t", pLocToGloMap->GetLocalToGlobalMap()[counter]);
                    dof_elems.emplace_back(pLocToGloMap->GetLocalToGlobalMap()[counter] + 1);
                    ++counter;
                }
                m_l2g.emplace_back(dof_elems);
                dof_elems.clear();
            }

            auto tend = clock();
            auto t = double(tend - tbegin) / CLOCKS_PER_SEC;
            double t_ave = 0.0;
            MPI_Reduce(&t, &t_ave, 1, MPI_DOUBLE, MPI_SUM, 0, m_comm);
            if(!rank) printf("Saena mesh info generation time: %f\n", t_ave / nprocs);

            // CONSTRUCT KSP OBJECT
            SetUpSolver(pLocToGloMap->GetIterativeTolerance());

            tbegin = clock();

            m_matrix.erase_no_shrink_to_fit();

            // POPULATE MATRIX
            for(n = cnt = 0; n < m_expList.lock()->GetNumElmts(); ++n)
            {
                loc_mat = GetBlock(n);
                loc_lda = loc_mat->GetRows();

                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalMap(cnt+i) - nDirDofs;
                    sign1 = pLocToGloMap->GetLocalToGlobalSign(cnt + i);
                    if(gid1 >= 0)
                    {
                        int gid1ro = m_reorderedMap[gid1];
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = pLocToGloMap->GetLocalToGlobalMap(cnt + j)
                                   - nDirDofs;
                            sign2 = pLocToGloMap->GetLocalToGlobalSign(cnt + j);
                            if(gid2 >= 0)
                            {
                                int gid2ro = m_reorderedMap[gid2];
                                value = sign1*sign2*(*loc_mat)(i,j);
                                m_matrix.set(gid1ro, gid2ro, value);
                            }
                        }
                    }
                }
                cnt += loc_lda;
            }

            // timing
            tend = clock();
            t = double(tend - tbegin) / CLOCKS_PER_SEC;
            MPI_Reduce(&t, &t_ave, 1, MPI_DOUBLE, MPI_SUM, 0, m_comm);
            if(!rank) printf("nektar assembly time: %f\n", t_ave / nprocs);

            tbegin = clock();

            // ASSEMBLE MATRIX
//            m_matrix.set_num_threads(1);
            m_matrix.assemble(m_scale);
//            m_matrix.assemble_writeToFile("matrix_folder");

            // timing
            tend = clock();
            t = double(tend - tbegin) / CLOCKS_PER_SEC;
            MPI_Reduce(&t, &t_ave, 1, MPI_DOUBLE, MPI_SUM, 0, m_comm);
            if(!rank) printf("Saena matrix assembly time: %f\n", t_ave / nprocs);
        }


        GlobalLinSysSaenaFull::~GlobalLinSysSaenaFull()
        {

        }


        /**
        * Solve the linear system using a full global matrix system.
        */
        void GlobalLinSysSaenaFull::v_Solve(
            const Array<OneD, const NekDouble>  &pLocInput,
            Array<OneD,       NekDouble>  &pLocOutput,
            const AssemblyMapSharedPtr &pLocToGloMap,
            const Array<OneD, const NekDouble>  &pDirForcing)
        {
            std::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
            bool dirForcCalculated = (bool) pDirForcing.size();
            int nDirDofs  = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobDofs = pLocToGloMap->GetNumGlobalCoeffs();
            int nLocDofs  = pLocToGloMap->GetNumLocalCoeffs();

//            m_locToGloMap = pLocToGloMap; // required for DoMatrixMultiply

            Array<OneD, NekDouble> tmp(nLocDofs);
            Array<OneD, NekDouble> tmp1(nLocDofs);
            Array<OneD, NekDouble> global(nGlobDofs,0.0);

            int nDirTotal = nDirDofs;
            expList->GetComm()->GetRowComm()
                ->AllReduce(nDirTotal, LibUtilities::ReduceSum);

            if(nDirTotal)
            {
                // calculate the dirichlet forcing
                if(dirForcCalculated)
                {
                    // assume pDirForcing is in local space
                    ASSERTL0(pDirForcing.size() >= nLocDofs,
                             "DirForcing is not of sufficient size. Is it in local space?");
                    Vmath::Vsub(nLocDofs, pLocInput, 1,
                                pDirForcing, 1,tmp1, 1);
                }
                else
                {
                    // Calculate the dirichlet forcing and substract it
                    // from the rhs
                    expList->GeneralMatrixOp(
                        m_linSysKey, pLocOutput, tmp);

                    // Apply robin boundary conditions to the solution.
                    for(auto &r : m_robinBCInfo) // add robin mass matrix
                    {
                        RobinBCInfoSharedPtr rBC;
                        Array<OneD, NekDouble> tmploc;

                        int n  = r.first;

                        int offset = expList->GetCoeff_Offset(n);
                        LocalRegions::ExpansionSharedPtr vExp = expList->GetExp(n);

                        // add local matrix contribution
                        for(rBC = r.second;rBC; rBC = rBC->next)
                        {
                            vExp->AddRobinEdgeContribution(rBC->m_robinID,
                                                           rBC->m_robinPrimitiveCoeffs,
                                                           pLocOutput + offset,
                                                           tmploc = tmp + offset);
                        }
                    }

                    Vmath::Vsub(nLocDofs, pLocInput, 1, tmp, 1, tmp1, 1);
                }

                pLocToGloMap->Assemble(tmp1,tmp);

                SolveLinearSystem(nGlobDofs,tmp, global, pLocToGloMap, nDirDofs);

                pLocToGloMap->GlobalToLocal(global,tmp);

                // Add back initial and boundary condition
                Vmath::Vadd(nLocDofs, tmp, 1, pLocOutput, 1, pLocOutput, 1);
            }
            else
            {
                pLocToGloMap->Assemble(pLocInput,tmp);
                SolveLinearSystem(nGlobDofs, tmp,global, pLocToGloMap);
                pLocToGloMap->GlobalToLocal(global,pLocOutput);
            }
        }
    }
}
