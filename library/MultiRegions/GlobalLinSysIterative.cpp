///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysIterative.cpp
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
// Description: GlobalLinSysIterative definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysIterative.h>
#include <MultiRegions/PreconditionerDiagonal.h>
#include <MultiRegions/Preconditioner.h>

#include <MultiRegions/ConjugateGradient_plain.hxx>
#include <MultiRegions/ConjugateGradient_Kokkos.hxx>
#include <MultiRegions/ConjugateGradient_OpenMP.hxx>

#include <mutex>
#include "omp.h"

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysIterative
         *
         * Solves a linear system using iterative methods.
         */

        /// Constructor for full direct matrix solve.
        GlobalLinSysIterative::GlobalLinSysIterative(
                const GlobalLinSysKey &pKey,
                const boost::weak_ptr<ExpList> &pExpList,
                const boost::shared_ptr<AssemblyMap>
                &pLocToGloMap)
                : GlobalLinSys(pKey, pExpList, pLocToGloMap),
                  m_rhs_magnitude(NekConstants::kNekUnsetDouble),
                  m_rhs_mag_sm(0.9),
                  m_precon(NullPreconditionerSharedPtr),
                  m_totalIterations(0),
                  m_useProjection(false),
                  m_numPrevSols(0)
        {
            m_tolerance = pLocToGloMap->GetIterativeTolerance();
            m_maxiter   = pLocToGloMap->GetMaxIterations();

            LibUtilities::CommSharedPtr vComm = m_expList.lock()->GetComm()->GetRowComm();
            m_root    = (vComm->GetRank())? false : true;

            int successiveRHS;
            
            if((successiveRHS = pLocToGloMap->GetSuccessiveRHS()))
            {
                m_prevLinSol.set_capacity(successiveRHS);
                m_useProjection = true;
            }
            else
            {
                m_useProjection = false;
            }
        }

        GlobalLinSysIterative::~GlobalLinSysIterative()
        {
        }


        /**
         * 
         */
        void GlobalLinSysIterative::v_SolveLinearSystem(
                    const int nGlobal,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &plocToGloMap,
                    const int nDir)
        {
            printf("Within GlobalLinSysIterative::v_SolveLinearSystem\n" );        

            if (m_useProjection)
            {
                DoAconjugateProjection(nGlobal, pInput, pOutput, plocToGloMap, nDir);
            }
            else
            {
                int version;
                cout << "What parallelisation type do you wish to choose?" << endl
                     << "(1) Original Version" << endl 
                     << "(2) Plain Version" << endl 
                     << "(3) Full Kokkos Version" << endl 
                     << "(4) OpenMP Version" << endl;
                //cin >> version;
                version = 3;

                if(version == 1)
                {
                    DoConjugateGradient(nGlobal, pInput, pOutput, plocToGloMap, nDir);
                }
                else if(version == 2)
                {
                    DoConjugateGradient_plain(nGlobal, pInput, pOutput, plocToGloMap, nDir);
                }
                else if(version == 3)
                {
                    DoConjugateGradient_Kokkos(nGlobal, pInput, pOutput, plocToGloMap, nDir);
                }
                else if(version == 4)
                {
                    DoConjugateGradient_OpenMP(nGlobal, pInput, pOutput, plocToGloMap, nDir);
                }
                else
                {
                    ASSERTL0(false, "no parallelisation type set");    
                }
            }
        }


        /**
         * This method implements A-conjugate projection technique
         * in order to speed up successive linear solves with
         * right-hand sides arising from time-dependent discretisations.
         * (P.F.Fischer, Comput. Methods Appl. Mech. Engrg. 163, 1998)
         */
        void GlobalLinSysIterative::DoAconjugateProjection(
                    const int nGlobal,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &plocToGloMap,
                    const int nDir)
        {
            printf("Within GlobalLinSysIterative::DoAconjugateProjection\n" );

            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                                = m_expList.lock()->GetComm()->GetRowComm();

            // Get vector sizes
            int nNonDir = nGlobal - nDir;
            Array<OneD, NekDouble> tmp;

            if (0 == m_numPrevSols)
            {
                // no previous solutions found, call CG

                DoConjugateGradient(nGlobal, pInput, pOutput, plocToGloMap, nDir);

                UpdateKnownSolutions(nGlobal, pOutput, nDir);
            }
            else
            {
                // Create NekVector wrappers for linear algebra operations
                NekVector<NekDouble> b     (nNonDir, pInput  + nDir, eWrapper);
                NekVector<NekDouble> x     (nNonDir, tmp = pOutput + nDir, eWrapper);

                // check the input vector (rhs) is not zero

                NekDouble rhsNorm = Vmath::Dot2(nNonDir,
                                                pInput + nDir,
                                                pInput + nDir,
                                                m_map + nDir);

                vComm->AllReduce(rhsNorm, Nektar::LibUtilities::ReduceSum);

                if (rhsNorm < NekConstants::kNekZeroTol)
                {
                    Array<OneD, NekDouble> tmp = pOutput+nDir;
                    Vmath::Zero(nNonDir, tmp, 1);
                    return;
                }

                // Allocate array storage
                Array<OneD, NekDouble> px_s       (nGlobal, 0.0);
                Array<OneD, NekDouble> pb_s       (nGlobal, 0.0);
                Array<OneD, NekDouble> tmpAx_s    (nGlobal, 0.0);
                Array<OneD, NekDouble> tmpx_s     (nGlobal, 0.0);

                NekVector<NekDouble> pb    (nNonDir, tmp = pb_s    + nDir, eWrapper);
                NekVector<NekDouble> px    (nNonDir, tmp = px_s    + nDir, eWrapper);
                NekVector<NekDouble> tmpAx (nNonDir, tmp = tmpAx_s + nDir, eWrapper);
                NekVector<NekDouble> tmpx  (nNonDir, tmp = tmpx_s  + nDir, eWrapper);


                // notation follows the paper cited:
                // \alpha_i = \tilda{x_i}^T b^n
                // projected x, px = \sum \alpha_i \tilda{x_i}

                Array<OneD, NekDouble> alpha     (m_prevLinSol.size(), 0.0);
                for (int i = 0; i < m_prevLinSol.size(); i++)
                {
                    alpha[i] = Vmath::Dot2(nNonDir,
                                           m_prevLinSol[i],
                                           pInput + nDir,
                                           m_map + nDir);
                }
                vComm->AllReduce(alpha, Nektar::LibUtilities::ReduceSum);

                for (int i = 0; i < m_prevLinSol.size(); i++)
                {
                    if (alpha[i] < NekConstants::kNekZeroTol)
                    {
                        continue;
                    }

                    NekVector<NekDouble> xi (nNonDir, m_prevLinSol[i], eWrapper);
                    px += alpha[i] * xi;
                }

                // pb = b^n - A px
                Vmath::Vcopy(nNonDir,
                             pInput.get() + nDir, 1,
                             pb_s.get()   + nDir, 1);

                v_DoMatrixMultiply(px_s, tmpAx_s);

                pb -= tmpAx;


                // solve the system with projected rhs
                DoConjugateGradient(nGlobal, pb_s, tmpx_s, plocToGloMap, nDir);


                // remainder solution + projection of previous solutions
                x = tmpx + px;

                // save the auxiliary solution to prev. known solutions
                UpdateKnownSolutions(nGlobal, tmpx_s, nDir);
            }
        }

        
        /**
         * Calculating A-norm of an input vector,
         * A-norm(x) := sqrt( < x, Ax > )
         */
        NekDouble GlobalLinSysIterative::CalculateAnorm(
                                                        const int nGlobal,
                                                        const Array<OneD,const NekDouble> &in,
                                                        const int nDir)
        {
            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                = m_expList.lock()->GetComm()->GetRowComm();

            // Get vector sizes
            int nNonDir = nGlobal - nDir;

            // Allocate array storage
            Array<OneD, NekDouble> tmpAx_s    (nGlobal, 0.0);

            v_DoMatrixMultiply(in, tmpAx_s);

            NekDouble anorm_sq = Vmath::Dot2(nNonDir,
                                             in      + nDir,
                                             tmpAx_s + nDir,
                                             m_map   + nDir);
            vComm->AllReduce(anorm_sq, Nektar::LibUtilities::ReduceSum);
            return std::sqrt(anorm_sq);
        }

        /**
         * Updates the storage of previously known solutions.
         * Performs normalisation of input vector wrt A-norm.
         */
        void GlobalLinSysIterative::UpdateKnownSolutions(
                                                         const int nGlobal,
                                                         const Array<OneD,const NekDouble> &newX,
                                                         const int nDir)
        {
            // Get vector sizes
            int nNonDir = nGlobal - nDir;

            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                = m_expList.lock()->GetComm()->GetRowComm();

            // Check the solution is non-zero
            NekDouble solNorm = Vmath::Dot2(nNonDir,
                                            newX + nDir,
                                            newX + nDir,
                                            m_map + nDir);
            vComm->AllReduce(solNorm, Nektar::LibUtilities::ReduceSum);

            if (solNorm < NekConstants::kNekZeroTol)
            {
                return;
            }


            // Allocate array storage
            Array<OneD, NekDouble> tmpAx_s    (nGlobal, 0.0);
            Array<OneD, NekDouble> px_s       (nGlobal, 0.0);
            Array<OneD, NekDouble> tmp1, tmp2;

            // Create NekVector wrappers for linear algebra operations
            NekVector<NekDouble> px           (nNonDir, tmp1 = px_s    + nDir, eWrapper);
            NekVector<NekDouble> tmpAx        (nNonDir, tmp2 = tmpAx_s + nDir, eWrapper);


            // calculating \tilda{x} - sum \alpha_i\tilda{x}_i

            Vmath::Vcopy(nNonDir,
                         tmp1 = newX + nDir, 1,
                         tmp2 = px_s + nDir, 1);

            if (m_prevLinSol.size() > 0)
            {
                v_DoMatrixMultiply(newX, tmpAx_s);
            }

            Array<OneD, NekDouble> alpha (m_prevLinSol.size(), 0.0);
            for (int i = 0; i < m_prevLinSol.size(); i++)
            {
                alpha[i] = Vmath::Dot2(nNonDir,
                                       m_prevLinSol[i],
                                       tmpAx_s + nDir,
                                       m_map + nDir);
            }
            vComm->AllReduce(alpha, Nektar::LibUtilities::ReduceSum);

            for (int i = 0; i < m_prevLinSol.size(); i++)
            {
                if (alpha[i] < NekConstants::kNekZeroTol)
                {
                    continue;
                }

                NekVector<NekDouble> xi (nNonDir, m_prevLinSol[i], eWrapper);
                px -= alpha[i] * xi;
            }


            // Some solutions generated by CG are identical zeros, see
            // solutions generated for Test_Tet_equitri.xml (IncNavierStokesSolver).
            // Not going to store identically zero solutions.

            NekDouble anorm = CalculateAnorm(nGlobal, px_s, nDir);
            if (anorm < NekConstants::kNekZeroTol)
            {
                return;
            }

            // normalisation of new solution
            Vmath::Smul(nNonDir, 1.0/anorm, px_s.get() + nDir, 1, px_s.get() + nDir, 1);

            // updating storage with non-Dirichlet-dof part of new solution vector
            m_prevLinSol.push_back(px_s + nDir);
            m_numPrevSols++;
        }



        /**  
         * Solve a global linear system using the conjugate gradient method.  
         * We solve only for the non-Dirichlet modes. The operator is evaluated  
         * using an auxiliary function v_DoMatrixMultiply defined by the  
         * specific solver. Distributed math routines are used to support  
         * parallel execution of the solver.  
         *  
         * The implemented algorithm uses a reduced-communication reordering of  
         * the standard PCG method (Demmel, Heath and Vorst, 1993)  
         *  
         * @param       pInput      Input residual  of all DOFs.  
         * @param       pOutput     Solution vector of all DOFs.  
         */
        void GlobalLinSysIterative::DoConjugateGradient(
            const int                          nGlobal,
            const Array<OneD,const NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const AssemblyMapSharedPtr        &plocToGloMap,
            const int                          nDir)
        {
            if (!m_precon)
            {
                v_UniqueMap();
                m_precon = CreatePrecon(plocToGloMap);
                m_precon->BuildPreconditioner();
            }

            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                = m_expList.lock()->GetComm()->GetRowComm();

            // Get vector sizes
            int nNonDir = nGlobal - nDir;

            // Allocate array storage
            Array<OneD, NekDouble> w_A    (nGlobal, 0.0);
            Array<OneD, NekDouble> s_A    (nGlobal, 0.0);
            Array<OneD, NekDouble> p_A    (nNonDir, 0.0);
            Array<OneD, NekDouble> r_A    (nNonDir, 0.0);
            Array<OneD, NekDouble> q_A    (nNonDir, 0.0);
            Array<OneD, NekDouble> tmp;

            // Create NekVector wrappers for linear algebra operations
            NekVector<NekDouble> in (nNonDir,pInput  + nDir,      eWrapper);
            NekVector<NekDouble> out(nNonDir,tmp = pOutput + nDir,eWrapper);
            NekVector<NekDouble> w  (nNonDir,tmp = w_A + nDir,    eWrapper);
            NekVector<NekDouble> s  (nNonDir,tmp = s_A + nDir,    eWrapper);
            NekVector<NekDouble> p  (nNonDir,p_A,                 eWrapper);
            NekVector<NekDouble> r  (nNonDir,r_A,                 eWrapper);
            NekVector<NekDouble> q  (nNonDir,q_A,                 eWrapper);

            int k;
            NekDouble alpha, beta, rho, rho_new, mu, eps,  min_resid;
            Array<OneD, NekDouble> vExchange(3,0.0);

            // Copy initial residual from input
            r = in;
            // zero homogeneous out array ready for solution updates
            // Should not be earlier in case input vector is same as
            // output and above copy has been peformed
            Vmath::Zero(nNonDir,tmp = pOutput + nDir,1);


            // evaluate initial residual error for exit check
            vExchange[2] = Vmath::Dot2(nNonDir,
                                       r_A,
                                       r_A,
                                       m_map + nDir);

            vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);
            
            eps       = vExchange[2];

            if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
            {
                NekVector<NekDouble> inGlob (nGlobal, pInput, eWrapper);
                Set_Rhs_Magnitude(inGlob);
            }

            m_totalIterations = 0;

            // If input residual is less than tolerance skip solve.
            if (eps < m_tolerance * m_tolerance * m_rhs_magnitude)
            {
                if (m_verbose && m_root)
                {
                    cout << "CG iterations made = " << m_totalIterations 
                         << " using tolerance of "  << m_tolerance 
                         << " (error = " << sqrt(eps/m_rhs_magnitude) 
                         << ", rhs_mag = " << sqrt(m_rhs_magnitude) <<  ")" 
                         << endl;
                }
                return;
            }

            // Timing
            Timer t;
            t.Start();

            m_precon->DoPreconditioner(r_A, tmp = w_A + nDir);

            v_DoMatrixMultiply(w_A, s_A);

            k = 0;

            vExchange[0] = Vmath::Dot2(nNonDir,
                                       r_A,
                                       w_A + nDir,
                                       m_map + nDir);

            vExchange[1] = Vmath::Dot2(nNonDir,
                                       s_A + nDir,
                                       w_A + nDir,
                                       m_map + nDir);

            vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

            rho               = vExchange[0];
            mu                = vExchange[1];
            min_resid         = m_rhs_magnitude;
            beta              = 0.0;
            alpha             = rho/mu;
            m_totalIterations = 1;

            // Continue until convergence
            while (true)
            {
                if(k >= m_maxiter)
                {
                    if (m_root)
                    {
                        cout << "CG iterations made = " << m_totalIterations 
                             << " using tolerance of "  << m_tolerance 
                             << " (error = " << sqrt(eps/m_rhs_magnitude)
                             << ", rhs_mag = " << sqrt(m_rhs_magnitude) <<  ")"
                             << endl;
                    }
                    ROOTONLY_NEKERROR(ErrorUtil::efatal,
                                      "Exceeded maximum number of iterations");
                }

                // Compute new search direction p_k, q_k
                Vmath::Svtvp(nNonDir, beta, &p_A[0], 1, &w_A[nDir], 1, &p_A[0], 1);
                Vmath::Svtvp(nNonDir, beta, &q_A[0], 1, &s_A[nDir], 1, &q_A[0], 1);

                // Update solution x_{k+1}
                Vmath::Svtvp(nNonDir, alpha, &p_A[0], 1, &pOutput[nDir], 1, &pOutput[nDir], 1);

                // Update residual vector r_{k+1}
                Vmath::Svtvp(nNonDir, -alpha, &q_A[0], 1, &r_A[0], 1, &r_A[0], 1);

                // Apply preconditioner
                m_precon->DoPreconditioner(r_A, tmp = w_A + nDir);

                // Perform the method-specific matrix-vector multiply operation.
                 v_DoMatrixMultiply(w_A, s_A);

                // <r_{k+1}, w_{k+1}>
                vExchange[0] = Vmath::Dot2(nNonDir,
                                           r_A,
                                           w_A + nDir,
                                           m_map + nDir);
                // <s_{k+1}, w_{k+1}>
                vExchange[1] = Vmath::Dot2(nNonDir,
                                           s_A + nDir,
                                           w_A + nDir,
                                           m_map + nDir);

                // <r_{k+1}, r_{k+1}>
                vExchange[2] = Vmath::Dot2(nNonDir,
                                           r_A,
                                           r_A,
                                           m_map + nDir);

                // Perform inner-product exchanges
                vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

                rho_new = vExchange[0];
                mu      = vExchange[1];
                eps     = vExchange[2];

                m_totalIterations++;
                // test if norm is within tolerance
                if (eps < m_tolerance * m_tolerance * m_rhs_magnitude)
                {
                    if (m_verbose && m_root)
                    {
                        cout << "CG iterations made = " << m_totalIterations 
                             << " using tolerance of "  << m_tolerance 
                             << " (error = " << sqrt(eps/m_rhs_magnitude)
                             << ", rhs_mag = " << sqrt(m_rhs_magnitude) <<  ")"
                             << endl;
                    }
                    break;
                }
                min_resid = min(min_resid, eps);

                // Compute search direction and solution coefficients
                beta  = rho_new/rho;
                alpha = rho_new/(mu - rho_new*beta/alpha);
                rho   = rho_new;
                k++;
            }

            t.Stop();
            cout << "Time to compute: " << t.TimePerTest(1) << endl;
        }


        void GlobalLinSysIterative::Set_Rhs_Magnitude(
            const NekVector<NekDouble> &pIn)
        {
            Array<OneD, NekDouble> vExchange(1, 0.0);
            if (m_map.num_elements() > 0)
            {
                vExchange[0] = Vmath::Dot2(pIn.GetDimension(),
                                        &pIn[0],&pIn[0],&m_map[0]);
            }

            m_expList.lock()->GetComm()->GetRowComm()->AllReduce(
                vExchange, Nektar::LibUtilities::ReduceSum);

            // To ensure that very different rhs values are not being
            // used in subsequent solvers such as the velocit solve in
            // INC NS. If this works we then need to work out a better
            // way to control this.
            NekDouble new_rhs_mag = (vExchange[0] > 1e-6)? vExchange[0] : 1.0;

            if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
            {
                m_rhs_magnitude = new_rhs_mag;
            }
            else
            {
                m_rhs_magnitude = (m_rhs_mag_sm*(m_rhs_magnitude) + 
                                   (1.0-m_rhs_mag_sm)*new_rhs_mag); 
            }
        }

        std::vector<std::vector<int> > GlobalLinSysIterative::CreateColours(
                Array<OneD, const int> &localToGlobalMap, int ncoeffs, int el_begin, int el_end)
        {
            // ===== create graph of connected elements =====
            std::vector<std::vector<int>> el_connect;
            el_connect.resize(el_end); // has lots of zeros at the beginning, could refactor indexing to save memory
            for (int e_n = el_begin; e_n < el_end; ++e_n)
            {                
                for (int e_m = e_n+1; e_m < el_end; ++e_m)
                {                    
                    for (int n = 0; n < ncoeffs; ++n)
                    {
                        int map_n = localToGlobalMap[n + e_n * ncoeffs];
                        for (int m = 0; m < ncoeffs; ++m)
                        {                            
                            int map_m = localToGlobalMap[m + e_m * ncoeffs];
                            if (map_n == map_m)
                            {
                                el_connect[e_n].push_back(e_m);
                                el_connect[e_m].push_back(e_n);
                                n = ncoeffs; // to break out of second level loop
                                break;
                            }
                        }
                    }
                }
            }
            /*for(int e_n = el_begin; e_n < el_end; ++e_n)
            {
                printf("element %i: ",e_n);
                for (std::vector<int>::iterator it = el_connect[e_n].begin() ; it != el_connect[e_n].end(); ++it)
                {
                    printf("%i ",*it );
                }
                printf("\n");
            }
            printf("completed connecting of partition %i\n", par );*/

            // ===== create element colourgroups based on graph =====
            std::vector<int> remain(el_end - el_begin);
            for (int i = 0; i < remain.size(); ++i) // or use std::iota
            {
                remain[i] = el_begin + i;
            }
            std::vector<std::vector<int> > coloursets;
            
            // loop until all free elements have been sorted
            while (remain.size() > 0)
            {
                std::vector<int> colour;
                std::set<int> locked;
                std::set<int> completed;
                for (int i = 0; i < remain.size(); i++)
                {
                    if (locked.find(remain[i]) == locked.end()) // if the element is not locked
                    {                            
                        colour.push_back(remain[i]);   // insert element into colour
                        completed.insert(remain[i]);    // insert sorted element into "completed" list
                        for (int j = 0; j < el_connect[remain[i]].size(); j++)   // loop over all connected element
                        {
                            locked.insert(el_connect[remain[i]][j]);   // and flag these elements as locked
                        }
                    }                       
                }
                // identify elements which are not sorted, yet and create new "remain" vector
                std::vector<int> tmp = remain;
                remain.clear();
                for (int i = 0; i < tmp.size(); i++)
                {
                    if (completed.find(tmp[i]) == completed.end()) // if the element is not completed
                    {
                        remain.push_back(tmp[i]);
                    }
                }
                // include colour into vector of coloursets
                coloursets.push_back(colour);
            }

            // ===== distribute elements more evenly between groups =====
            int iter = coloursets.size() / 2 + 1;
            for (int it = 0; it < iter; ++it)
            {
                // find smallest and largest colourset
                int set_min = 0;
                int set_min_size = coloursets[0].size();
                int set_max = 0;
                int set_max_size = coloursets[0].size();
                for (int set = 1; set < coloursets.size(); ++set)
                {
                    set_min = set_min_size < coloursets[set].size() ? set_min : set;
                    set_min_size = coloursets[set_min].size();
                    set_max = set_max_size > coloursets[set].size() ? set_max : set;
                    set_max_size = coloursets[set_max].size();
                }

                // move elements from largest set to smallest set
                int move_size = (set_max_size - set_min_size) / 2;
                int move = 0;
                // check connections of all all elements of largest set with all elements of smallest set               
                for (int i = set_max_size-1; i >= 0 ; i--)
                {
                    // take element in largest set
                    int el_max = coloursets[set_max][i];
                    bool connection = false;
                    for (int j = 0; j < set_min_size; ++j)
                    {
                        // take alement in smallest set
                        int el_min = coloursets[set_min][j];
                        // if there is a connection between the two elements
                        if (std::find(el_connect[el_max].begin(), el_connect[el_max].end(), el_min)
                                 != el_connect[el_max].end() )
                        {
                            connection = true;
                            break;
                        }
                    }
                    if (connection == false) // if there is no connection
                    {
                        // move element from largest set to smallest set
                        coloursets[set_min].push_back(el_max);
                        coloursets[set_max].erase(coloursets[set_max].begin() + i);
                        move++;
                    }
                    if (move == move_size) // dont move more elements than necessary
                    {
                        break;
                    }
                }
            }

            return coloursets;
        }

        std::vector<std::vector<int> > GlobalLinSysIterative::CreateColourSets(
                Array<OneD, const int> &localToGlobalMap, int ncoeffs, int total_elmts, int threads)
        {
            // divide the list of elements into partitions
            int min_blocks = 896; // minimum number of blocks to fully occupy a GPU, 896 for GTX 1080
            int num_cs = 8; // typical number of coloursets depending on element type, 8 for triangular elements
            int partitions = std::max(1,total_elmts / (min_blocks*num_cs));
            //int partitions = threads;

            int part_size = total_elmts / partitions;
            std::vector<int> elmts(partitions+1);            
            elmts[0] = 0;
            for (int i = 1; i < partitions; ++i)
            {
                elmts[i] = i * part_size;
            }
            elmts[partitions] = total_elmts;
            
            // colour each partition
            std::vector<std::vector<std::vector<int> > > coloursets;
            coloursets.resize(partitions);

            #pragma omp parallel for
            for (int par = 0; par < partitions; ++par)
            {
                int el_begin = elmts[par];
                int el_end = elmts[par+1];
                coloursets[par] = CreateColours(localToGlobalMap, ncoeffs, el_begin, el_end);
            }

            // combine the coloursets of all partitions
            std::vector<std::vector<int> > all_coloursets;
            for (int i = 0; i < partitions; ++i)
            {
                all_coloursets.insert(
                    all_coloursets.end(), coloursets[i].begin(), coloursets[i].end());
            }

            return all_coloursets;
        }        

    }
}
