///////////////////////////////////////////////////////////////////////////////
//
// File: ConjugateGradient_plain.hxx
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_CONJUGATEGRADIENT_PLAIN_HXX
#define NEKTAR_LIB_MULTIREGIONS_CONJUGATEGRADIENT_PLAIN_HXX

#include <MultiRegions/GlobalLinSysIterative.h>
#include <MultiRegions/ConjugateGradient_BLAS.hxx>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
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
            cout << "How many threads?" << endl;
            int threads;
            cin >> threads;
            args.num_threads = threads;
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
            int map[nGlobal], maxiter, k;
            NekDouble rhs_magnitude, tolerance;
            int totalIterations = 0;

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

            printf("base0: %i\n",base0.num_elements());
            printf("base1: %i\n",base1.num_elements());
            printf("nquad0: %i\n",nquad0);
            printf("nquad1: %i\n",nquad1);
            printf("nmodes0: %i\n",nmodes0);
            printf("nmodes1: %i\n",nmodes1);

            
            int metricSize = elmts * nquad0 * nquad1;
            Array<OneD, NekDouble> quadMetricGlo (4*metricSize);
            Array<OneD, NekDouble> laplacian00Glo(quadMetricGlo+metricSize);
            Array<OneD, NekDouble> laplacian01Glo(quadMetricGlo+2*metricSize);
            Array<OneD, NekDouble> laplacian11Glo(quadMetricGlo+3*metricSize);

            int numLocalCoeffs, numGlobalCoeffs;
            Array<OneD, const int> localToGlobalMap;
            Array<OneD, const NekDouble> localToGlobalSign; 
            GetMatrixMultiplyMetrics(
                        quadMetricGlo, laplacian00Glo, laplacian01Glo, laplacian11Glo,
                        nquad0, nquad1, elmts,
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
                Kokkos::finalize(); 
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
                    localToGlobalMap, localToGlobalSign, totalIterations);

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
                    Kokkos::finalize(); 
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
                    localToGlobalMap, localToGlobalSign, totalIterations);
                

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
                Array<OneD, NekDouble> &quadMetricGlo,                
                Array<OneD, NekDouble> &laplacian00Glo,
                Array<OneD, NekDouble> &laplacian01Glo,
                Array<OneD, NekDouble> &laplacian11Glo,
                int &nquad0, int &nquad1, int &elmts,
                int &numLocalCoeffs, int &numGlobalCoeffs,
                Array<OneD, const int> &localToGlobalMap,
                Array<OneD, const NekDouble> &localToGlobalSign)
        {
            printf("Within GlobalLinSysIterative::GetMatrixMultiplyMetrics\n" );           

            boost::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
            expList->GeneralMatrixOp_plain(
                        quadMetricGlo, laplacian00Glo, laplacian01Glo, laplacian11Glo,
                        nquad0, nquad1, elmts,
                        numLocalCoeffs, numGlobalCoeffs,
                        localToGlobalMap, localToGlobalSign);

        }


        void GlobalLinSysIterative::GeneralMatrixOp_plain(
                const Array<OneD,const NekDouble> &inarray,
                Array<OneD,      NekDouble> &outarray,
                const NekDouble &lambda,
                const Array<OneD, const NekDouble> &quadMetricGlo,                
                const Array<OneD, const NekDouble> &laplacian00Glo,
                const Array<OneD, const NekDouble> &laplacian01Glo,
                const Array<OneD, const NekDouble> &laplacian11Glo,
                const int &nquad0, const int &nquad1, 
                const int &nmodes0, const int &nmodes1, const int &ncoeffs, 
                const Array<OneD, const int>  &coeff_offset, const int &elmts,
                const Array<OneD, const NekDouble> &base0,
                const Array<OneD, const NekDouble> &base1,
                const Array<OneD, const NekDouble> &dbase0,
                const Array<OneD, const NekDouble> &dbase1,
                const DNekMatSharedPtr &D0, const DNekMatSharedPtr &D1,
                const int &numLocalCoeffs, const int &numGlobalCoeffs,
                const Array<OneD, const int> &localToGlobalMap,
                const Array<OneD, const NekDouble> &localToGlobalSign,
                const int iteration)
        {
            //printf("%s\n", "do the global to local mapping");
            Array<OneD,NekDouble> tmp1(numLocalCoeffs);
            
            //GlobalToLocal_plain(inarray,tmp1);
            //Vmath::Gathr(numLocalCoeffs, localToGlobalSign.get(),
            //         inarray.get(), localToGlobalMap.get(), tmp1.get());
            for (int i = 0; i < numLocalCoeffs; ++i)
            {
                tmp1[i] = localToGlobalSign[i] * inarray[localToGlobalMap[i]];
            }

            Kokkos::View<double*, Kokkos::HostSpace> transfer_out;
            transfer_out = Kokkos::View<double*, Kokkos::HostSpace>("transfer_out", numLocalCoeffs);                        

            GeneralMatrixOp_IterPerExp_plain(tmp1,transfer_out,lambda,
                    quadMetricGlo, laplacian00Glo,laplacian01Glo,laplacian11Glo,                    
                    nquad0, nquad1, nmodes0, nmodes1, ncoeffs,
                    coeff_offset, elmts,
                    base0, base1, dbase0, dbase1, D0, D1);
            
            //printf("%s\n", "do the local to global mapping");
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
                outarray[localToGlobalMap[i]] += localToGlobalSign[i] * transfer_out[i];                
            }
        }


        void GlobalLinSysIterative::GeneralMatrixOp_IterPerExp_plain(
                    const Array<OneD,const NekDouble> &inarray,
                    //Array<OneD,      NekDouble> &outarray,
                    Kokkos::View<double*, Kokkos::HostSpace> transfer_out,
                    const NekDouble &lambda,
                    const Array<OneD, const NekDouble> &quadMetricGlo,                
                    const Array<OneD, const NekDouble> &laplacian00Glo,
                    const Array<OneD, const NekDouble> &laplacian01Glo,
                    const Array<OneD, const NekDouble> &laplacian11Glo,
                    const int &nquad0, const int &nquad1, 
                    const int &nmodes0, const int &nmodes1, const int &ncoeffs, 
                    const Array<OneD, const int>  &coeff_offset, const int &elmts,
                    const Array<OneD, const NekDouble> &base0,
                    const Array<OneD, const NekDouble> &base1,
                    const Array<OneD, const NekDouble> &dbase0,
                    const Array<OneD, const NekDouble> &dbase1,
                    const DNekMatSharedPtr &D0, const DNekMatSharedPtr &D1)
        {
            printf("%s\n", "perform operations by element");            
            Kokkos::parallel_for(range_policy_host(0,elmts),KOKKOS_LAMBDA (const int el)
            {                                    
                printf("%i ", el);
                Array<OneD, NekDouble> tmp_inarray (ncoeffs);
                for (int i = 0; i < ncoeffs; ++i)
                {
                    tmp_inarray[i] = inarray[coeff_offset[el]+i];
                }
                HelmholtzMatrixOp_MatFree_plain(
                    tmp_inarray,
                    transfer_out,
                    el, coeff_offset,
                    lambda,
                     quadMetricGlo,
                    laplacian00Glo,
                    laplacian01Glo,
                    laplacian11Glo,
                    nquad0, nquad1, nmodes0, nmodes1, ncoeffs,
                    base0, base1, dbase0, dbase1, D0, D1);
            });
            printf("\n");             
        }


        void GlobalLinSysIterative::HelmholtzMatrixOp_MatFree_plain(
                const Array<OneD, const NekDouble> &inarray,
                Kokkos::View<double*, Kokkos::HostSpace> transfer_out,
                const int &el, const Array<OneD,
                const int>  &coeff_offset,
                const NekDouble &lambda,
                const Array<OneD, const NekDouble> &quadMetric,
                const Array<OneD, const NekDouble> &laplacian00,
                const Array<OneD, const NekDouble> &laplacian01,
                const Array<OneD, const NekDouble> &laplacian11,
                const int &nquad0, const int &nquad1,
                const int &nmodes0, const int &nmodes1, const int &ncoeffs,
                const Array<OneD, const NekDouble> &base0,
                const Array<OneD, const NekDouble> &base1,
                const Array<OneD, const NekDouble> &dbase0,
                const Array<OneD, const NekDouble> &dbase1,
                const DNekMatSharedPtr &D0, const DNekMatSharedPtr &D1)
        {
            //printf("within GlobalLinSysIterative::HelmholtzMatrixOp_MatFree_plain \n");
            
            Array<OneD, NekDouble> t_outarray(ncoeffs);
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
                wsp1[i] = quadMetric[el*nqtot+i] * wsp0[i];
            }            
            
            //IProductWRTBase_SumFacKernel(base0, base1, wsp1, outarray, wsp2, true, true);
            IProductWRTBase_SumFacKernel_plain(base0, base1, wsp1, t_outarray,
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
                wsp0L[i] = laplacian00[el*nqtot+i] * wsp1L[i] 
                         + laplacian01[el*nqtot+i] * wsp2L[i];
            }
            //Vmath::Vvtvvtp(nqtot,&metric01[0],1,&wsp1L[0],1,&metric11[0],1,&wsp2L[0],1,&wsp2L[0],1);
            for (int i = 0; i < nqtot; ++i)
            {
                wsp2L[i] = laplacian01[el*nqtot+i] * wsp1L[i] 
                         + laplacian11[el*nqtot+i] * wsp2L[i];
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
                //outarray[i] = lambda * t_outarray[i] + wsp1[i];
                transfer_out[coeff_offset[el] + i] = lambda * t_outarray[i] + wsp1[i];
            }    
            

        }

        void GlobalLinSysIterative::IProductWRTBase_SumFacKernel_plain(
                const Array<OneD, const NekDouble>& base0,
                const Array<OneD, const NekDouble>& base1,
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray,
                Array<OneD, NekDouble> &wsp,
                const int &nmodes0, const int &nmodes1,
                const int &nquad0, const int &nquad1)
        {
            //printf("within GlobalLinSysIterative::IProductWRTBase_SumFacKernel_plain \n");
            plainDgemm('T','N',nquad1,nmodes0,nquad0,1.0,inarray.get(),nquad0,
                        base0.get(),nquad0,0.0,wsp.get(),nquad1);
            int i, mode;
            for (mode=i=0; i < nmodes0; ++i)
            {
                plainDgemv('T',nquad1,nmodes1-i,1.0, base1.get()+mode*nquad1,
                            nquad1,wsp.get()+i*nquad1,1, 0.0,
                            outarray.get() + mode,1);
                mode += nmodes1 - i;
            }
            outarray[1] += plainDdot(nquad1,base1.get()+nquad1,1,
                                          wsp.get()+nquad1,1);
        }

        void GlobalLinSysIterative::BwdTrans_SumFacKernel_plain(
                const Array<OneD, const NekDouble>& base0,
                const Array<OneD, const NekDouble>& base1,
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray,
                Array<OneD, NekDouble> &wsp,
                const int &nmodes0, const int &nmodes1,
                const int &nquad0, const int &nquad1)
        {
            //printf("within GlobalLinSysIterative::BwdTrans_SumFacKernel_plain \n");
            int i, mode;
            for (i = mode = 0; i < nmodes0; ++i)
            {
                plainDgemv('N', nquad1,nmodes1-i,1.0,base1.get()+mode*nquad1,
                            nquad1,&inarray[0]+mode,1,0.0,&wsp[0]+i*nquad1,1);
                mode += nmodes1-i;
            }
            plainDaxpy(nquad1,inarray[1],base1.get()+nquad1,1,
                            &wsp[0]+nquad1,1);
            plainDgemm('N','T', nquad0,nquad1,nmodes0,1.0, base0.get(),nquad0,
                        &wsp[0], nquad1,0.0, &outarray[0], nquad0);            
        }

        void GlobalLinSysIterative::PhysTensorDeriv_plain(
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray_d0,
                Array<OneD, NekDouble> &outarray_d1,
                const int &nquad0, const int &nquad1,
                const DNekMatSharedPtr &D0, const DNekMatSharedPtr &D1)
        {
            //printf("within GlobalLinSysIterative::PhysTensorDeriv_plain \n");
            plainDgemm('N', 'N', nquad0, nquad1, nquad0, 1.0,
                        &(D0->GetPtr())[0], nquad0, &inarray[0], nquad0, 0.0,
                        &outarray_d0[0], nquad0);
            plainDgemm('N', 'T', nquad0, nquad1, nquad1, 1.0, &inarray[0], nquad0,
                         &(D1->GetPtr())[0], nquad1, 0.0, &outarray_d1[0], nquad0);
        }
        
    }
}

#endif