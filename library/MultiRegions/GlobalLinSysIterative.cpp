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

#include <Kokkos_Core.hpp>

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
                // applying plain Conjugate Gradient
                //DoConjugateGradient(nGlobal, pInput, pOutput, plocToGloMap, nDir);
                // Conjugate Gradient with plain arrays (Jan)
                DoConjugateGradient_plain(nGlobal, pInput, pOutput, plocToGloMap, nDir);
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



        void GlobalLinSysIterative::DoConjugateGradient_plain(
            const int                          nGlobal,
            const Array<OneD,const NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const AssemblyMapSharedPtr        &plocToGloMap,
            const int                          nDir)
        {
            printf("Within GlobalLinSysIterative::DoConjugateGradient_plain\n" );

            //Initialise Kokkos
            Kokkos::InitArguments args;
            args.num_threads = 1;
            Kokkos::initialize(args);   

            // create preconditioner
            Array<OneD, NekDouble> diagonals;
            if (!m_precon)
            {
                v_UniqueMap();
                m_precon = CreatePrecon(plocToGloMap);

                //m_precon->BuildPreconditioner();
                diagonals = m_precon->Preconditioner::DiagonalPreconditionerSum_plain();
            }

            // Get vector sizes
            int nNonDir = nGlobal - nDir;
            
            // Allocate array storage
            Array<OneD, NekDouble> w_A    (nGlobal, 0.0);
            Array<OneD, NekDouble> s_A    (nGlobal, 0.0);
            Array<OneD, NekDouble> p_A    (nNonDir, 0.0);
            Array<OneD, NekDouble> r_A    (nNonDir, 0.0);
            Array<OneD, NekDouble> q_A    (nNonDir, 0.0);
            
            NekDouble alpha, beta, rho, rho_new, mu, eps,  min_resid;
            int map[nGlobal], maxiter, totalIterations, k;
            NekDouble rhs_magnitude, tolerance;

            // copy member variables            
            for (int i = 0; i < nGlobal; ++i)
            {
                map[i] = m_map[i];
            }
            tolerance = m_tolerance;
            maxiter = m_maxiter;

            // Gather Data for Helmholtz Matrix Multiplication
            NekDouble lambda = m_linSysKey.GetConstFactor(StdRegions::eFactorLambda);

            int nquad0, nquad1, nmodes0, nmodes1, ncoeffs;
            Array<OneD, const NekDouble> base0, base1, dbase0, dbase1;
            DNekMatSharedPtr D0, D1;
            int elmts;
            Array<OneD, int> coeff_offset;

            boost::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
            expList->v_GetStdExpansionMetric(
                nquad0, nquad1, nmodes0, nmodes1, ncoeffs, 
                coeff_offset, elmts,
                base0, base1, dbase0, dbase1,
                D0, D1);
            
            int metricSize = elmts * nquad0 * nquad1;
            Array<OneD, NekDouble> quadMetricGlo (4*metricSize);
            Array<OneD, NekDouble> laplacian00Glo(quadMetricGlo+metricSize);
            Array<OneD, NekDouble> laplacian01Glo(quadMetricGlo+2*metricSize);
            Array<OneD, NekDouble> laplacian11Glo(quadMetricGlo+3*metricSize);

            int numLocalCoeffs, numGlobalCoeffs;
            Array<OneD, const int> localToGlobalMap;
            Array<OneD, const NekDouble> localToGlobalSign; 
            GetMatrixMultiplyMetrics(w_A, s_A, lambda,
                        quadMetricGlo, laplacian00Glo, laplacian01Glo, laplacian11Glo,
                        nquad0, nquad1, nmodes0, nmodes1, ncoeffs, 
                        coeff_offset, elmts,
                        base0, base1, dbase0, dbase1,
                        D0, D1,
                        numLocalCoeffs, numGlobalCoeffs,
                        localToGlobalMap, localToGlobalSign); 
            printf("%s\n", "==== completed data gathering ====");

            // Copy initial residual from input
            for (int i = 0; i < nNonDir; ++i)
            {
                r_A[i] = pInput[i+nDir];
            }
            
            // zero homogeneous out array ready for solution updates
            // Should not be earlier in case input vector is same as
            // output and above copy has been peformed
            for (int i = 0; i < nNonDir; ++i)
            {
                pOutput[i+nDir] = 0.0;
            }            

            // evaluate initial residual error for exit check
            eps = 0.0;
            for (int i = 0; i < nNonDir; ++i)
            {
                eps += r_A[i] * r_A[i] * map[i+nDir];
            }
            
            rhs_magnitude = 0.0;
            for (int i = 0; i < nGlobal; ++i)
            {
                rhs_magnitude += pInput[i] * pInput[i] * map[i];
            }                       

            // If input residual is less than tolerance skip solve.
            totalIterations = 0;
            if (eps < tolerance * tolerance * rhs_magnitude)
            {
                cout << "CG iterations made = " << totalIterations 
                     << " using tolerance of "  << tolerance 
                     << " (error = " << sqrt(eps/rhs_magnitude) 
                     << ", rhs_mag = " << sqrt(rhs_magnitude) <<  ")" 
                     << endl;
                return;
            }

            // Timing
            Timer t;
            t.Start();


            //preconditioner
            //m_precon->DoPreconditioner(r_A, tmp = w_A + nDir);
            //Vmath::Vmul(nNonDir, &r_A[0], 1, &diagonals[0], 1, &w_A[nDir], 1);
            for (int i = 0; i < nNonDir; ++i)
            {
                w_A[i+nDir] = r_A[i] * diagonals[i];
            }

            GeneralMatrixOp_plain(
                    w_A, s_A, lambda,
                    quadMetricGlo, laplacian00Glo, laplacian01Glo, laplacian11Glo,
                    nquad0, nquad1, nmodes0, nmodes1, ncoeffs, 
                    coeff_offset, elmts,
                    base0, base1, dbase0, dbase1,
                    D0, D1,
                    numLocalCoeffs, numGlobalCoeffs,
                    localToGlobalMap, localToGlobalSign);

            rho = 0.0;
            for (int i = 0; i < nNonDir; ++i)
            {
                rho += r_A[i] * w_A[i+nDir] * map[i+nDir];
            }

            mu = 0.0;
            for (int i = 0; i < nNonDir; ++i)
            {
                mu += s_A[i+nDir] * w_A[i+nDir] * map[i+nDir];
            }

            min_resid         = rhs_magnitude;
            beta              = 0.0;
            alpha             = rho/mu;            

            // Continue until convergence
            k = 0;
            totalIterations = 1;
            while (true)
            {
                if(k >= maxiter)
                {
                    cout << "CG iterations made = " << totalIterations 
                         << " using tolerance of "  << tolerance 
                         << " (error = " << sqrt(eps/rhs_magnitude)
                         << ", rhs_mag = " << sqrt(rhs_magnitude) <<  ")"
                         << endl;                    
                    ROOTONLY_NEKERROR(ErrorUtil::efatal,
                                      "Exceeded maximum number of iterations");
                }

                // Compute new search direction p_k, q_k
                //Vmath::Svtvp(nNonDir, beta, &p_A[0], 1, &w_A[nDir], 1, &p_A[0], 1);
                for (int i = 0; i < nNonDir; ++i)
                {
                    p_A[i] = beta * p_A[i] + w_A[i+nDir];                    
                }                
                
                //Vmath::Svtvp(nNonDir, beta, &q_A[0], 1, &s_A[nDir], 1, &q_A[0], 1);
                for (int i = 0; i < nNonDir; ++i)
                {
                    q_A[i] = beta * q_A[i] + s_A[i+nDir];                    
                } 
                // Update solution x_{k+1}
                //Vmath::Svtvp(nNonDir, alpha, &p_A[0], 1, &pOutput[nDir], 1, &pOutput[nDir], 1);
                for (int i = 0; i < nNonDir; ++i)
                {
                    pOutput[i+nDir] = alpha * p_A[i] + pOutput[i+nDir];                    
                } 
                // Update residual vector r_{k+1}
                //Vmath::Svtvp(nNonDir, -alpha, &q_A[0], 1, &r_A[0], 1, &r_A[0], 1);
                for (int i = 0; i < nNonDir; ++i)
                {
                    r_A[i] = -alpha * q_A[i] + r_A[i];                    
                } 


                // Apply preconditioner
                //m_precon->DoPreconditioner(r_A, tmp = w_A + nDir);
                //Vmath::Vmul(nNonDir, &r_A[0], 1, &diagonals[0], 1, &w_A[nDir], 1);
                for (int i = 0; i < nNonDir; ++i)
                {
                    w_A[i+nDir] = r_A[i] * diagonals[i];
                }

                // Perform the method-specific matrix-vector multiply operation.
                printf("CG iteration %i ==================================\n", totalIterations);
                GeneralMatrixOp_plain(
                    w_A, s_A, lambda,
                    quadMetricGlo, laplacian00Glo, laplacian01Glo, laplacian11Glo,
                    nquad0, nquad1, nmodes0, nmodes1, ncoeffs, 
                    coeff_offset, elmts,
                    base0, base1, dbase0, dbase1,
                    D0, D1,
                    numLocalCoeffs, numGlobalCoeffs,
                    localToGlobalMap, localToGlobalSign);
                

                rho_new = 0.0;
                for (int i = 0; i < nNonDir; ++i)
                {
                    rho_new += r_A[i] * w_A[i+nDir] * map[i+nDir];
                }

                mu = 0.0;
                for (int i = 0; i < nNonDir; ++i)
                {
                    mu += s_A[i+nDir] * w_A[i+nDir] * map[i+nDir];
                }

                eps = 0.0;
                for (int i = 0; i < nNonDir; ++i)
                {
                    eps += r_A[i] * r_A[i] * map[i+nDir];
                }
                // Compute search direction and solution coefficients
                beta  = rho_new/rho;
                alpha = rho_new/(mu - rho_new*beta/alpha);
                rho   = rho_new;
                min_resid = (eps < min_resid) ? eps : min_resid;

                k++;
                totalIterations++;
                // test if norm is within tolerance
                if (eps < tolerance * tolerance * rhs_magnitude)
                {
                    cout << "CG iterations made = " << totalIterations 
                         << " using tolerance of "  << tolerance 
                         << " (error = " << sqrt(eps/rhs_magnitude)
                         << ", rhs_mag = " << sqrt(rhs_magnitude) <<  ")"
                         << endl;
                    break;
                }
            }
            t.Stop();
            cout << "Time to compute: " << t.TimePerTest(1) << endl;

            Kokkos::finalize();     
        }




        void GlobalLinSysIterative::GetMatrixMultiplyMetrics(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput,
                const NekDouble &lambda,
                      Array<OneD, NekDouble> &quadMetricGlo,                
                Array<OneD, NekDouble> &laplacian00Glo,
                Array<OneD, NekDouble> &laplacian01Glo,
                Array<OneD, NekDouble> &laplacian11Glo,
                int &nquad0, int &nquad1, int &nmodes0, int &nmodes1, int &ncoeffs, 
                        Array<OneD, const int>  &coeff_offset, int &elmts,
                        Array<OneD, const NekDouble> &base0,
                        Array<OneD, const NekDouble> &base1,
                        Array<OneD, const NekDouble> &dbase0,
                        Array<OneD, const NekDouble> &dbase1,
                        DNekMatSharedPtr &D0, DNekMatSharedPtr &D1,
                        int &numLocalCoeffs, int &numGlobalCoeffs,
            Array<OneD, const int> &localToGlobalMap,
            Array<OneD, const NekDouble> &localToGlobalSign)
        {
            printf("Within GlobalLinSysIterative::GetMatrixMultiplyMetrics\n" );           

            boost::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
            expList->GeneralMatrixOp_plain(pInput, pOutput, lambda,
                        quadMetricGlo, laplacian00Glo, laplacian01Glo, laplacian11Glo,
                        nquad0, nquad1, nmodes0, nmodes1, ncoeffs, 
                        coeff_offset, elmts,
                        base0, base1, dbase0, dbase1,
                        D0, D1,
                        numLocalCoeffs, numGlobalCoeffs,
                        localToGlobalMap, localToGlobalSign);

        }


        void GlobalLinSysIterative::GeneralMatrixOp_plain(
                const Array<OneD,const NekDouble>  &inarray,
                      Array<OneD,      NekDouble>  &outarray,
                      const NekDouble &lambda,
                      Array<OneD, NekDouble> &quadMetricGlo,                
                Array<OneD, NekDouble> &laplacian00Glo,
                Array<OneD, NekDouble> &laplacian01Glo,
                Array<OneD, NekDouble> &laplacian11Glo,
                int &nquad0, int &nquad1, int &nmodes0, int &nmodes1, int &ncoeffs, 
                Array<OneD, const int>  &coeff_offset, int &elmts,
                Array<OneD, const NekDouble> &base0,
                Array<OneD, const NekDouble> &base1,
                Array<OneD, const NekDouble> &dbase0,
                Array<OneD, const NekDouble> &dbase1,
                DNekMatSharedPtr &D0, DNekMatSharedPtr &D1,
                int &numLocalCoeffs, int &numGlobalCoeffs,
                Array<OneD, const int> &localToGlobalMap,
                Array<OneD, const NekDouble> &localToGlobalSign)
        {
            printf("%s\n", "do the global to local mapping");
            Array<OneD,NekDouble> tmp1(2*numLocalCoeffs);
            Array<OneD,NekDouble> tmp2(tmp1+numLocalCoeffs);

            //GlobalToLocal_plain(inarray,tmp1);
            //Vmath::Gathr(numLocalCoeffs, localToGlobalSign.get(),
            //         inarray.get(), localToGlobalMap.get(), tmp1.get());
            for (int i = 0; i < numLocalCoeffs; ++i)
            {
                tmp1[i] = localToGlobalSign[i] * inarray[localToGlobalMap[i]];
            }

            GeneralMatrixOp_IterPerExp_plain(tmp1,tmp2,lambda,
                    quadMetricGlo, laplacian00Glo,laplacian01Glo,laplacian11Glo,                    
                    nquad0, nquad1, nmodes0, nmodes1, ncoeffs,
                    coeff_offset, elmts,
                    base0, base1, dbase0, dbase1, D0, D1);

            printf("%s\n", "do the local to global mapping");
            //Assemble_plain(tmp2,outarray);  
            //Vmath::Zero(numGlobalCoeffs, outarray.get(), 1);
            for (int i = 0; i < numGlobalCoeffs; ++i)
            {
                outarray[i] = 0.0;
            }
            //Vmath::Assmb(numLocalCoeffs, localToGlobalSign.get(), 
            //            tmp2.get(), localToGlobalMap.get(), outarray.get());
            for (int i = 0; i < numLocalCoeffs; ++i)
            {
                outarray[localToGlobalMap[i]] += localToGlobalSign[i] * tmp2[i]; 
            }
            
        }


        void GlobalLinSysIterative::GeneralMatrixOp_IterPerExp_plain(
                    const Array<OneD,const NekDouble> &inarray,
                    Array<OneD,      NekDouble> &outarray,
                    const NekDouble &lambda,
                    Array<OneD, NekDouble> &quadMetricGlo,                
                    Array<OneD, NekDouble> &laplacian00Glo,
                    Array<OneD, NekDouble> &laplacian01Glo,
                    Array<OneD, NekDouble> &laplacian11Glo,
                    int &nquad0, int &nquad1, int &nmodes0, int &nmodes1, int &ncoeffs, 
                    Array<OneD, const int>  &coeff_offset, int &elmts,
                    Array<OneD, const NekDouble> &base0,
                    Array<OneD, const NekDouble> &base1,
                    Array<OneD, const NekDouble> &dbase0,
                    Array<OneD, const NekDouble> &dbase1,
                    DNekMatSharedPtr &D0, DNekMatSharedPtr &D1)
        {
            printf("%s\n", "perform operations by element");            

            //Array<OneD, NekDouble> quadMetric (nquad0*nquad1);
            //Array<OneD, NekDouble> laplacian00(nquad0*nquad1);
            //Array<OneD, NekDouble> laplacian01(nquad0*nquad1);
            //Array<OneD, NekDouble> laplacian11(nquad0*nquad1);
            Array<OneD, NekDouble> tmp_outarray(ncoeffs);
            // Calculating
            for(int el = 0; el < elmts; ++el)
            {
                printf("%i ", el);
                /*for (int i = 0; i < nquad0*nquad1; ++i)
                {
                    quadMetric[i] = quadMetricGlo[el*nquad0*nquad1+i];
                    laplacian00[i] = laplacian00Glo[el*nquad0*nquad1+i];
                    laplacian01[i] = laplacian01Glo[el*nquad0*nquad1+i];
                    laplacian11[i] = laplacian11Glo[el*nquad0*nquad1+i];
                }*/
                HelmholtzMatrixOp_MatFree_plain(
                    inarray + coeff_offset[el],
                    tmp_outarray,
                    lambda,
                    quadMetricGlo + el*nquad0*nquad1 ,
                    laplacian00Glo + el*nquad0*nquad1,
                    laplacian01Glo + el*nquad0*nquad1,
                    laplacian11Glo + el*nquad0*nquad1,
                    nquad0, nquad1, nmodes0, nmodes1, ncoeffs,
                    base0, base1, dbase0, dbase1, D0, D1);

                for (int i = 0; i < ncoeffs; ++i)
                {
                    outarray[coeff_offset[el]+i] = tmp_outarray[i];
                }                
            }
            printf("\n");           
        }


        void GlobalLinSysIterative::HelmholtzMatrixOp_MatFree_plain(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble>  &outarray,
                const NekDouble &lambda,
                const Array<OneD, NekDouble> &quadMetric,
                const Array<OneD, NekDouble> &laplacian00,
                const Array<OneD, NekDouble> &laplacian01,
                const Array<OneD, NekDouble> &laplacian11,
                int &nquad0, int &nquad1, int &nmodes0, int &nmodes1, int &ncoeffs,
                const Array<OneD, const NekDouble> &base0,
                const Array<OneD, const NekDouble> &base1,
                const Array<OneD, const NekDouble> &dbase0,
                const Array<OneD, const NekDouble> &dbase1,
                DNekMatSharedPtr &D0, DNekMatSharedPtr &D1)
        {
            //printf("%s\n", "within GlobalLinSysIterative::HelmholtzMatrixOp_MatFree_plain");
            
            int       nqtot   = nquad0*nquad1;
            //int       wspsize = std::max(std::max(std::max(nqtot,ncoeffs),nquad1*nmodes0), nquad0*nmodes1);
            int max1 = (nqtot >= ncoeffs) ? nqtot : ncoeffs;
            int max2 = (nquad1*nmodes0 >= nquad0*nmodes1) ? nquad1*nmodes0 : nquad0*nmodes1;
            int wspsize = (max1 >= max2) ? max1 : max2;

            // Allocate temporary storage
            Array<OneD,NekDouble> wsp0(5*wspsize);       // size wspsize
            Array<OneD,NekDouble> wsp1(wsp0 + wspsize);  // size wspsize // u_hat
            Array<OneD,NekDouble> wsp2(wsp0 + 2*wspsize);// size 3*wspsize
            
            //BwdTrans_SumFacKernel(base0, base1, inarray, wsp0, wsp2,true,true);            
            BwdTrans_SumFacKernel_plain(base0, base1, inarray, wsp0, wsp2,
                nmodes0, nmodes1, nquad0, nquad1);

            //MultiplyByQuadratureMetric(wsp0, wsp1);
            //Vmath::Vmul(nqtot, quadMetric, 1, wsp0, 1, wsp1, 1);            
            for (int i = 0; i < nqtot; ++i)
            {
                wsp1[i] = quadMetric[i] * wsp0[i];
            }

            //IProductWRTBase_SumFacKernel(base0, base1, wsp1, outarray, wsp2, true, true);
            IProductWRTBase_SumFacKernel_plain(base0, base1, wsp1, outarray,
                                         wsp2, nmodes0, nmodes1, nquad0, nquad1);

            //LaplacianMatrixOp_MatFree_Kernel(wsp0, wsp1, wsp2);
            // Allocate temporary storage
            Array<OneD,NekDouble> wsp0L(wsp2);//
            Array<OneD,NekDouble> wsp1L(wsp2+wspsize);
            Array<OneD,NekDouble> wsp2L(wsp2+2*wspsize);
            
            //StdExpansion2D::PhysTensorDeriv(wsp0,wsp1L,wsp2L);
            PhysTensorDeriv_plain(wsp0,wsp1L,wsp2L, nquad0, nquad1, D0, D1);

            //Vmath::Vvtvvtp(nqtot,&metric00[0],1,&wsp1L[0],1,&metric01[0],1,&wsp2L[0],1,&wsp0L[0],1);
            for (int i = 0; i < nqtot; ++i)
            {
                wsp0L[i] = laplacian00[i] * wsp1L[i] + laplacian01[i] * wsp2L[i];
            }
            //Vmath::Vvtvvtp(nqtot,&metric01[0],1,&wsp1L[0],1,&metric11[0],1,&wsp2L[0],1,&wsp2L[0],1);
            for (int i = 0; i < nqtot; ++i)
            {
                wsp2L[i] = laplacian01[i] * wsp1L[i] + laplacian11[i] * wsp2L[i];
            }

            //IProductWRTBase_SumFacKernel(dbase0, base1,wsp0L,wsp1 ,wsp1L);
            //IProductWRTBase_SumFacKernel( base0,dbase1,wsp2L,wsp1L,wsp0L);
            IProductWRTBase_SumFacKernel_plain(dbase0, base1,wsp0L,wsp1 ,wsp1L,
                         nmodes0, nmodes1, nquad0, nquad1);
            IProductWRTBase_SumFacKernel_plain( base0,dbase1,wsp2L,wsp1L,wsp0L,
                         nmodes0, nmodes1, nquad0, nquad1);

            //Vmath::Vadd(m_ncoeffs,wsp1L.get(),1,wsp1.get(),1,wsp1.get(),1);
            for (int i = 0; i < ncoeffs; ++i)
            {
                wsp1[i] += wsp1L[i];
            }
            //end LaplacianMatrixOp_MatFree_Kernel
           
            // outarray = lambda * outarray + wsp1
            //          = (lambda * M + L ) * u_hat
            //Vmath::Svtvp(m_ncoeffs, lambda, &outarray[0], 1, &wsp1[0], 1, &outarray[0], 1);
            for (int i = 0; i < ncoeffs; ++i)
            {
                outarray[i] = lambda * outarray[i] + wsp1[i];
            }    
            

        }

        void GlobalLinSysIterative::IProductWRTBase_SumFacKernel_plain(
                const Array<OneD, const NekDouble>& base0,
                const Array<OneD, const NekDouble>& base1,
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray,
                Array<OneD, NekDouble> &wsp,
                int &nmodes0, int &nmodes1,
                int &nquad0, int &nquad1)
        {
            Blas::Dgemm('T','N',nquad1,nmodes0,nquad0,1.0,inarray.get(),nquad0,
                        base0.get(),nquad0,0.0,wsp.get(),nquad1);
            int i, mode;
            for (mode=i=0; i < nmodes0; ++i)
            {
                Blas::Dgemv('T',nquad1,nmodes1-i,1.0, base1.get()+mode*nquad1,
                            nquad1,wsp.get()+i*nquad1,1, 0.0,
                            outarray.get() + mode,1);
                mode += nmodes1 - i;
            }
            outarray[1] += Blas::Ddot(nquad1,base1.get()+nquad1,1,
                                          wsp.get()+nquad1,1);
        }

        void GlobalLinSysIterative::BwdTrans_SumFacKernel_plain(
                const Array<OneD, const NekDouble>& base0,
                const Array<OneD, const NekDouble>& base1,
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray,
                Array<OneD, NekDouble> &wsp,
                int &nmodes0, int &nmodes1,
                int &nquad0, int &nquad1)
        {
            int i, mode;
            for (i = mode = 0; i < nmodes0; ++i)
            {
                Blas::Dgemv('N', nquad1,nmodes1-i,1.0,base1.get()+mode*nquad1,
                            nquad1,&inarray[0]+mode,1,0.0,&wsp[0]+i*nquad1,1);
                mode += nmodes1-i;
            }
            Blas::Daxpy(nquad1,inarray[1],base1.get()+nquad1,1,
                            &wsp[0]+nquad1,1);
            Blas::Dgemm('N','T', nquad0,nquad1,nmodes0,1.0, base0.get(),nquad0,
                        &wsp[0], nquad1,0.0, &outarray[0], nquad0);            
        }

        void GlobalLinSysIterative::PhysTensorDeriv_plain(
            const Array<OneD, const NekDouble>& inarray,
                         Array<OneD, NekDouble> &outarray_d0,
                         Array<OneD, NekDouble> &outarray_d1,
                         int &nquad0, int &nquad1,
                         DNekMatSharedPtr &D0, DNekMatSharedPtr &D1)
        {
            Blas::Dgemm('N', 'N', nquad0, nquad1, nquad0, 1.0,
                        &(D0->GetPtr())[0], nquad0, &inarray[0], nquad0, 0.0,
                        &outarray_d0[0], nquad0);
            Blas:: Dgemm('N', 'T', nquad0, nquad1, nquad1, 1.0, &inarray[0], nquad0,
                         &(D1->GetPtr())[0], nquad1, 0.0, &outarray_d1[0], nquad0);
        }




        

    }
}
