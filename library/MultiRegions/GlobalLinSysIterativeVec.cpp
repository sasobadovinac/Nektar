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
// Description: GlobalLinSysIterativeVec definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysIterativeVec.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysIterativeVec
         *
         * Solves a linear system using iterative methods.
         */

        /// Constructor for full direct matrix solve.
        GlobalLinSysIterativeVec::GlobalLinSysIterativeVec(
                const GlobalLinSysKey &pKey,
                const Array<OneD, std::weak_ptr<ExpList> > &pVecExpList,
                const Array<OneD, std::shared_ptr<AssemblyMap> > &pVecLocToGloMap)
                : GlobalLinSys(pKey, pVecExpList[0], pVecLocToGloMap[0]),
                  m_expListVec(pVecExpList),
                  m_rhs_magnitude(NekConstants::kNekUnsetDouble),
                  m_rhs_mag_sm(0.9),
                  m_totalIterations(0)
        {
            m_tolerance = pVecLocToGloMap[0]->GetIterativeTolerance();
            m_maxiter   = pVecLocToGloMap[0]->GetMaxIterations();

            LibUtilities::CommSharedPtr vComm = m_expListVec[0].lock()->GetComm()->GetRowComm();
            m_root    = (vComm->GetRank())? false : true;

            // possibly see about adding sucessive RHS later
            m_useProjection = false;

        }

        GlobalLinSysIterativeVec::~GlobalLinSysIterativeVec()
        {
        }


        /**
         * 
         */
        void GlobalLinSysIterativeVec::v_SolveVecLinearSystem(
                    const Array<OneD, int >& nGlobal,
                    const Array<OneD, const Array<OneD, NekDouble> > &pvecInput,
                    Array <OneD, Array<OneD, NekDouble> > &pvecOutput,
                    const Array<OneD, AssemblyMapSharedPtr > &pvecLocToGloMap,
                    const Array<OneD, int> & nDir)
        {
            // applying plain Conjugate Gradient
            DoConjugateGradient(nGlobal, pvecInput, pvecOutput, pvecLocToGloMap, nDir);
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
        void GlobalLinSysIterativeVec::DoConjugateGradient(
            const Array<OneD, int >                    &nGlobal,
            const Array<OneD, Array<OneD, NekDouble> > &pvecInput,
                  Array<OneD, Array<OneD, NekDouble> > &pvecOutput,
            const Array<OneD, AssemblyMapSharedPtr >   &pvecLocToGloMap,
            const Array<OneD, int >                    &nDir)
        {

            ASSERTL1(nGlobal.num_elements() == pvecInput.num_elements(),
                     "Exect nGlobal to be of same dimension as pvecInput");
            ASSERTL1(pvecOutput.num_elements() == pvecInput.num_elements(),
                     "Exect pvecInput to be of same dimension as pvecInput");
            ASSERTL1(pvecLocToGloMap.num_elements() == pvecInput.num_elements(),
                "Exect pvecLocToGloMap to be of same dimension as pvecInput");
            ASSERTL1(nDir.num_elements() == pvecInput.num_elements(),
                "Exect nDir to be of same dimension as pvecInput");

            int nvec = pvecInput.num_elements();
            
            if (m_preconVec.num_elements() == 0)
            {
                v_UniqueMap();
                m_preconVec = Array<OneD, PreconditionerSharedPtr>(nvec);
                for(int i = 0; i < nvec; ++i)
                {
                    m_preconVec[i] = CreatePrecon(pvecLocToGloMap[i]);
                    m_preconVec[i]->BuildPreconditioner();
                }
            }

            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                = m_expListVec[0].lock()->GetComm()->GetRowComm();

            // Get vector sizes
            int nNonDirTot, nGlobalTot, nDirTot;
            Array<OneD, int> nNonDir(nvec,0);

            nNonDirTot = nGlobalTot = nDirTot = 0; 
            for(int i =0; i < nvec; ++i)
            {
                nGlobalTot += nGlobal[i];
                nDirTot    += nDir[i];
                nNonDir[i]  = nGlobal[i] - nDir[i];
                nNonDirTot += nNonDir[i];
            }

            // Allocate array storage
            Array<OneD, NekDouble> w_A    (nGlobalTot, 0.0);
            Array<OneD, NekDouble> s_A    (nGlobalTot, 0.0);
            Array<OneD, NekDouble> p_A    (nNonDirTot, 0.0);
            Array<OneD, NekDouble> r_A    (nNonDirTot, 0.0);
            Array<OneD, NekDouble> q_A    (nNonDirTot, 0.0);
            Array<OneD, NekDouble> tmp,tmp1;

            int i,k;
            NekDouble alpha, beta, rho, rho_new, mu, eps,  min_resid;
            Array<OneD, NekDouble> vExchange(3,0.0);

            int cnt  = 0;
            int cnt1 = 0; 
            for(i = 0; i < nvec; ++i, cnt += nNonDir[i])
            {
                // Copy initial residual from input
                Vmath::Vcopy(nNonDir[i],pvecInput[i],1,tmp = r_A+cnt,1);

                // zero homogeneous out array ready for solution updates
                // Should not be earlier in case input vector is same as
                // output and above copy has been peformed
                Vmath::Zero(nNonDir[i],tmp = pvecOutput[i] + nDir[i],1);
                
                // evaluate initial residual error for exit check
                vExchange[2] += Vmath::Dot2(nNonDir[i],
                                           tmp = r_A+cnt,
                                           tmp1 = r_A+cnt,
                                           m_mapVec[i] + nDir[i]);
            }
            
            vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);
                
            eps       = vExchange[2];
            
            if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
            {
                Set_Rhs_Magnitude(pvecInput[0].num_elements(), pvecInput);
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

            {
                m_preconVec[i]->DoPreconditioner(r_A+cnt, tmp = w_A + cnt +  nDir[i]);
            }
            
            v_DoMatrixMultiply(nGlobal,w_A, s_A);

            k = 0;


            vExchange[0]  = vExchange[1] = 0.0;
            for(i = 0,cnt=0,cnt1=0; i < nvec;
                        cnt +=nNonDir[i],cnt1 += nGlobal[i],++i)
            {
                vExchange[0] += Vmath::Dot2(nNonDir[i],
                                           tmp  = r_A + cnt,
                                           tmp1 = w_A + cnt1 + nDir[i],
                                           m_mapVec[i] + nDir[i]);
                
                vExchange[1] += Vmath::Dot2(nNonDir[i],
                                            tmp  = s_A + cnt1 + nDir[i],
                                            tmp1 = w_A + cnt1 + nDir[i],
                                            m_mapVec[i] + nDir[i]);
            }
            
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

                
                for(i = 0,cnt=0,cnt1=0; i < nvec;
                               cnt +=nNonDir[i],cnt1 += nGlobal[i],++i)
                {
                    // Compute new search direction p_k, q_k
                    Vmath::Svtvp(nNonDir[i], beta, &p_A[cnt], 1,
                                 &w_A[cnt1+nDir[i]], 1,
                                 &p_A[cnt], 1);
                    Vmath::Svtvp(nNonDir[i], beta, &q_A[cnt], 1,
                                 &s_A[cnt1+nDir[i]], 1,
                                 &q_A[cnt], 1);

                    // Update solution x_{k+1}
                    Vmath::Svtvp(nNonDir[i], alpha, &p_A[cnt], 1,
                                 &pvecOutput[i][nDir[i]], 1,
                                 &pvecOutput[i][nDir[i]], 1);

                    // Update residual vector r_{k+1}
                    Vmath::Svtvp(nNonDir[i], -alpha, &q_A[cnt], 1, &r_A[cnt], 1,
                                 &r_A[cnt], 1);

                    // Apply preconditioner
                    m_preconVec[i]->DoPreconditioner(r_A+cnt,
                                               tmp = w_A + cnt1 + nDir[i]);
                }
                
                // Perform the method-specific matrix-vector multiply operation.
                v_DoMatrixMultiply(nGlobal,w_A, s_A);
                
                vExchange[0] = vExchange[1] = vExchange[2] = 0.0;
                for(i = 0,cnt=0,cnt1=0; i < nvec;
                               cnt +=nNonDir[i],cnt1 += nGlobal[i],++i)
                {
                    // <r_{k+1}, w_{k+1}>
                    vExchange[0] += Vmath::Dot2(nNonDir[i],
                                                tmp  = r_A + cnt,
                                                tmp1 = w_A + cnt1 + nDir[i],
                                                m_mapVec[i] + nDir[i]);
                    // <s_{k+1}, w_{k+1}>
                    vExchange[1] += Vmath::Dot2(nNonDir[i],
                                                tmp  = s_A + cnt1 + nDir[i],
                                                tmp1 = w_A + cnt1 + nDir[i],
                                                m_mapVec[i] + nDir[i]);
                    // <r_{k+1}, r_{k+1}>
                    vExchange[2] += Vmath::Dot2(nNonDir[i],
                                                tmp  = r_A + cnt,
                                                tmp1 = r_A + cnt,
                                                m_mapVec[i] + nDir[i]);
                }
                
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
        }

        void GlobalLinSysIterativeVec::Set_Rhs_Magnitude(
                                       int ntest,
                                       const Array<OneD, Array<OneD, NekDouble> > &pIn)
        {
            Array<OneD, NekDouble> vExchange(1, 0.0);
            if (m_mapVec[0].num_elements() > 0)
            {
                for(int i = 0;  i < pIn.num_elements(); ++i)
                {
                    vExchange[0] += Vmath::Dot2(ntest,
                                               &pIn[i][0],&pIn[i][0],&m_mapVec[i][0]);
                }
            }

            m_expListVec[0].lock()->GetComm()->GetRowComm()->AllReduce(
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

    }
}
