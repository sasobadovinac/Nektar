///////////////////////////////////////////////////////////////////////////////
//
// File: ConjugateGradient_Kokkos.hxx
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
#ifndef NEKTAR_LIB_MULTIREGIONS_CONJUGATEGRADIENT_KOKKOS_HXX
#define NEKTAR_LIB_MULTIREGIONS_CONJUGATEGRADIENT_KOKKOS_HXX

#include <MultiRegions/GlobalLinSysIterative.h>
#include <MultiRegions/ConjugateGradient_BLAS.hxx>

//#ifdef __CUDA_ARCH__
#include "cudaProfiler.h"
//#endif

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        void GlobalLinSysIterative::DoConjugateGradient_Kokkos(
            const int                          nGlobal,
            const Array<OneD,const NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const AssemblyMapSharedPtr        &plocToGloMap,
            const int                          nDir)
    	{
    		printf("Within GlobalLinSysIterative::DoConjugateGradient_Kokkos\n" ); 

            //Initialise Kokkos]
            cout << "How many threads?" << endl;
            int threads;
            //cin >> threads;
            threads = 6;
            Kokkos::InitArguments args;
            args.num_threads = threads;
            args.device_id = 0; // 0 - 1080, 1 - 970
            Kokkos::initialize(args);
            printf("%s\n","kokkos initialised");

            // Get vector sizes
            int nNonDir = nGlobal - nDir;

            // create preconditioner
            Array<OneD, NekDouble> NekDiagonals(nNonDir);
            if (!m_precon)
            {
                v_UniqueMap();
                m_precon = CreatePrecon(plocToGloMap);

                //m_precon->BuildPreconditioner();
                NekDiagonals = m_precon->Preconditioner::DiagonalPreconditionerSum_plain();
            }
            // copy preconditioner to device
            Kokkos::View<double*> diagonals ("diagonals", nNonDir); 
            typename Kokkos::View< double*>::HostMirror h_diagonals;
            h_diagonals = Kokkos::create_mirror_view(diagonals);
            Kokkos::parallel_for(range_policy_host(0,nNonDir), KOKKOS_LAMBDA (const int i)
    		{
            	h_diagonals(i) = NekDiagonals[i];
            });
            Kokkos::deep_copy(diagonals,h_diagonals);                        

            // copy mapping            
            Kokkos::View<double*> map ("map", nGlobal);                        
            typename Kokkos::View< double*>::HostMirror h_map;
            h_map = Kokkos::create_mirror_view(map);
            Kokkos::parallel_for(range_policy_host(0,nGlobal), KOKKOS_LAMBDA (const int i)
    		{
                h_map(i) = m_map[i];
            });
            Kokkos::deep_copy(map,h_map);

             // copying lambda
            Kokkos::View<double[1],random_memory> lambda ("lambda");                        
            typename Kokkos::View< double[1]>::HostMirror h_lambda;
            h_lambda = Kokkos::create_mirror_view(lambda);
            h_lambda[0] = m_linSysKey.GetConstFactor(StdRegions::eFactorLambda);
            Kokkos::deep_copy(lambda,h_lambda);

            // -------------------------------------------------------------------------------
            // Gather constant element data for Helmholtz Matrix Multiplication
            int nquad0, nquad1, nmodes0, nmodes1, ncoeffs, elmts;
            Array<OneD, const NekDouble> Nekbase0, Nekbase1, Nekdbase0, Nekdbase1;
            DNekMatSharedPtr NekD0, NekD1; // length nquad0, nquad1
            Array<OneD, int> Nekcoeff_offset; //ncoeffs

            boost::shared_ptr<MultiRegions::ExpList> expList = m_expList.lock();
            expList->v_GetStdExpansionMetric(
                nquad0, nquad1, nmodes0, nmodes1, ncoeffs, 
                Nekcoeff_offset, elmts,
                Nekbase0, Nekbase1, Nekdbase0, Nekdbase1,
                NekD0, NekD1);

            // copying base0 and dbase0
            int base0_len = Nekbase0.num_elements();
            printf("base0_len = %i \n",base0_len );
            
            Kokkos::View<double*,random_memory> base0 ("base0", base0_len);                        
            typename Kokkos::View< double*>::HostMirror h_base0;
            h_base0 = Kokkos::create_mirror_view(base0);
            Kokkos::View<double*,random_memory> dbase0 ("dbase0", base0_len);                        
            typename Kokkos::View< double*>::HostMirror h_dbase0;
            h_dbase0 = Kokkos::create_mirror_view(dbase0);
            Kokkos::parallel_for(range_policy_host(0,base0_len), KOKKOS_LAMBDA (const int i)
    		{
            	h_base0(i) = Nekbase0[i];
            	h_dbase0(i) = Nekdbase0[i];
            });
            Kokkos::deep_copy(base0,h_base0);
            Kokkos::deep_copy(dbase0,h_dbase0);

            // copying base1 and dbase1
            int base1_len = Nekbase1.num_elements(); //nmodes0*(nmodes1-0.5*(nmodes0-1))*nquad1;
            printf("base1_len = %i \n",base1_len );

            Kokkos::View<double*,random_memory> base1 ("base1", base1_len);                        
            typename Kokkos::View< double*>::HostMirror h_base1;
            h_base1 = Kokkos::create_mirror_view(base1);
            Kokkos::View<double*,random_memory> dbase1 ("dbase1", base1_len);                        
            typename Kokkos::View< double*>::HostMirror h_dbase1;
            h_dbase1 = Kokkos::create_mirror_view(dbase1);
            Kokkos::parallel_for(range_policy_host(0,base1_len), KOKKOS_LAMBDA (const int i)
    		{
            	h_base1(i) = Nekbase1[i];
            	h_dbase1(i) = Nekdbase1[i];
            });
            Kokkos::deep_copy(base1,h_base1);
            Kokkos::deep_copy(dbase1,h_dbase1);

            // copy coeff_offset to device 
            Kokkos::View<int*,random_memory> coeff_offset ("coeff_offset", elmts);                        
            typename Kokkos::View< int*>::HostMirror h_coeff_offset;
            h_coeff_offset = Kokkos::create_mirror_view(coeff_offset);
            Kokkos::parallel_for(range_policy_host(0,elmts), KOKKOS_LAMBDA (const int i)
    		{
            	h_coeff_offset(i) = Nekcoeff_offset[i];
            });
            Kokkos::deep_copy(coeff_offset,h_coeff_offset);

            // copying DNekMatSharedPtr NekD0
            Kokkos::View<double*,random_memory> D0 ("D0", nquad0*nquad0);                        
            typename Kokkos::View< double*>::HostMirror h_D0;
            h_D0 = Kokkos::create_mirror_view(D0);
            NekDouble* D0_raw = &(NekD0->GetPtr())[0];
            Kokkos::parallel_for(range_policy_host(0,nquad0*nquad0), KOKKOS_LAMBDA (const int i)
    		{
            	h_D0(i) = *(D0_raw+i);
            });
            Kokkos::deep_copy(D0,h_D0);
            // copying DNekMatSharedPtr NekD1
            Kokkos::View<double*,random_memory> D1 ("D1", nquad1*nquad1);                        
            typename Kokkos::View< double*>::HostMirror h_D1;
            h_D1 = Kokkos::create_mirror_view(D1);
            NekDouble* D1_raw = &(NekD1->GetPtr())[0];
            Kokkos::parallel_for(range_policy_host(0,nquad1*nquad1), KOKKOS_LAMBDA (const int i)
    		{
            	h_D1(i) = *(D1_raw+i);
            });
            Kokkos::deep_copy(D1,h_D1);


            // --------------------------------------------------------------------------------
            // Gathering  quadrature and laplacian metrics for each element
            int metricSize = nquad0 * nquad1;
            Kokkos::View<double**,random_memory> quadMetricGlo  ("quadMetricGlo" , metricSize, elmts);                        
            typename Kokkos::View< double**>::HostMirror h_quadMetricGlo ;
            h_quadMetricGlo  = Kokkos::create_mirror_view(quadMetricGlo );
            Kokkos::View<double**,random_memory> laplacian00Glo ("laplacian00Glo", metricSize, elmts);                        
            typename Kokkos::View< double**>::HostMirror h_laplacian00Glo;
            h_laplacian00Glo = Kokkos::create_mirror_view(laplacian00Glo);
            Kokkos::View<double**,random_memory> laplacian01Glo ("laplacian01Glo", metricSize, elmts);                        
            typename Kokkos::View< double**>::HostMirror h_laplacian01Glo;
            h_laplacian01Glo = Kokkos::create_mirror_view(laplacian01Glo);
            Kokkos::View<double**,random_memory> laplacian11Glo ("laplacian11Glo", metricSize, elmts);                        
            typename Kokkos::View< double**>::HostMirror h_laplacian11Glo;
            h_laplacian11Glo = Kokkos::create_mirror_view(laplacian11Glo);

            int numLocalCoeffs, numGlobalCoeffs;
            Array<OneD, const int> NeklocalToGlobalMap; //size = numLocalCoeffs
            Array<OneD, const NekDouble> NeklocalToGlobalSign; //size = numLocalCoeffs            

            Array<OneD, NekDouble> NekquadMetricGlo (4*metricSize*elmts);
            Array<OneD, NekDouble> Neklaplacian00Glo(NekquadMetricGlo+metricSize*elmts);
            Array<OneD, NekDouble> Neklaplacian01Glo(NekquadMetricGlo+2*metricSize*elmts);
            Array<OneD, NekDouble> Neklaplacian11Glo(NekquadMetricGlo+3*metricSize*elmts);
            
            GetMatrixMultiplyMetrics(
                        NekquadMetricGlo, Neklaplacian00Glo, Neklaplacian01Glo, Neklaplacian11Glo,
                        nquad0, nquad1, elmts,
                        numLocalCoeffs, numGlobalCoeffs,
                        NeklocalToGlobalMap, NeklocalToGlobalSign);

            // copying metrics

            Kokkos::parallel_for(range_policy_host(0,metricSize), KOKKOS_LAMBDA (const int i)
    		{
            	for (int el = 0; el < elmts; ++el)
                {
                    h_quadMetricGlo (i,el) = NekquadMetricGlo [i+el*metricSize];
            	    h_laplacian00Glo(i,el) = Neklaplacian00Glo[i+el*metricSize];
            	    h_laplacian01Glo(i,el) = Neklaplacian01Glo[i+el*metricSize];
            	    h_laplacian11Glo(i,el) = Neklaplacian11Glo[i+el*metricSize];
                }
            });
            Kokkos::deep_copy(quadMetricGlo,h_quadMetricGlo);
            Kokkos::deep_copy(laplacian00Glo,h_laplacian00Glo);
            Kokkos::deep_copy(laplacian01Glo,h_laplacian01Glo);
            Kokkos::deep_copy(laplacian11Glo,h_laplacian11Glo);
            
            Kokkos::View<int*> localToGlobalMap ("localToGlobalMap", numLocalCoeffs);                        
            typename Kokkos::View< int*>::HostMirror h_localToGlobalMap;
            h_localToGlobalMap = Kokkos::create_mirror_view(localToGlobalMap);
            Kokkos::View<double*> localToGlobalSign ("localToGlobalSign", numLocalCoeffs);                        
            typename Kokkos::View< double*>::HostMirror h_localToGlobalSign;
            h_localToGlobalSign = Kokkos::create_mirror_view(localToGlobalSign);
            Kokkos::parallel_for(range_policy_host(0,numLocalCoeffs), KOKKOS_LAMBDA (const int i)
    		{
            	h_localToGlobalMap(i) = NeklocalToGlobalMap[i];
            	h_localToGlobalSign(i) = NeklocalToGlobalSign[i];
            });
            Kokkos::deep_copy(localToGlobalMap,h_localToGlobalMap);
            Kokkos::deep_copy(localToGlobalSign,h_localToGlobalSign);


            // creating and copying coloursets
            std::vector<std::vector<int> > coloursets = CreateColourSets(
                    NeklocalToGlobalMap, ncoeffs, elmts, threads);
            /*for(int cs = 0; cs < coloursets.size(); ++cs)
            {
                printf("colourset %i: ",cs);
                for (std::vector<int>::iterator it = coloursets[cs].begin() ; it != coloursets[cs].end(); ++it)
                {
                    printf("%i ",*it );
                }
                printf("\n");
            }*/
            for(int cs = 0; cs < coloursets.size(); ++cs)
            {
                printf("colourset %i: %i\n",cs, coloursets[cs].size() );
            }

            int ncs = coloursets.size();
            int max_cs = 0;
            int cs_sizes[ncs];
            for (int i = 0; i < ncs; ++i)
            {
                cs_sizes[i] = coloursets[i].size();
                max_cs = (max_cs < cs_sizes[i] ? cs_sizes[i] : max_cs);
            }            

            Kokkos::View<int**> coloursetArray ("coloursetArray", ncs, max_cs);                        
            typename Kokkos::View<int**>::HostMirror h_coloursetArray;
            h_coloursetArray = Kokkos::create_mirror_view(coloursetArray);
            for (int cs = 0; cs < ncs; ++cs)
            {
                for (int el = 0; el < cs_sizes[cs]; ++el)
                {
                    h_coloursetArray(cs,el) = coloursets[cs][el];
                }
            }        
            Kokkos::deep_copy(coloursetArray,h_coloursetArray);

            
            printf("%s\n", "finished data gathering");
            //--------------------------------------------------------------------------

            // CG variables
            NekDouble alpha, beta, rho, rho_new, mu, eps, min_resid, rhs_magnitude;
            int totalIterations = 0;
            int k;

            NekDouble tolerance = m_tolerance;
            int maxiter = m_maxiter;
            
            // Allocate array storage
            Kokkos::View<double*> w_A;
            w_A = Kokkos::View<double*>("w_A", nGlobal);                        
            typename Kokkos::View< double*>::HostMirror h_w_A;
            h_w_A = Kokkos::create_mirror_view(w_A);

            Kokkos::View<double*> s_A ("s_A", nGlobal);                        
            typename Kokkos::View< double*>::HostMirror h_s_A;
            h_s_A = Kokkos::create_mirror_view(s_A);

            Kokkos::View<double*> p_A ("p_A", nNonDir);                        
            typename Kokkos::View< double*>::HostMirror h_p_A;
            h_p_A = Kokkos::create_mirror_view(p_A);

            Kokkos::View<double*> r_A ("r_A", nNonDir);                        
            typename Kokkos::View< double*>::HostMirror h_r_A;
            h_r_A = Kokkos::create_mirror_view(r_A);

            Kokkos::View<double*> q_A ("q_A", nNonDir);                        
            typename Kokkos::View< double*>::HostMirror h_q_A;
            h_q_A = Kokkos::create_mirror_view(q_A);

            // Copy initial residual from input
            Kokkos::parallel_for(range_policy_host(0,nNonDir), KOKKOS_LAMBDA (const int i)
    		{
                h_r_A[i] = pInput[i+nDir];
            });
            Kokkos::deep_copy(r_A,h_r_A);

            // zero homogeneous out array ready for solution updates
            // Should not be earlier in case input vector is same as
            // output and above copy has been peformed
            Kokkos::View<double*> Output ("Output", nGlobal);                        
            typename Kokkos::View< double*>::HostMirror h_Output;
            h_Output = Kokkos::create_mirror_view(Output);
            Kokkos::parallel_for(range_policy_host(0,nNonDir), KOKKOS_LAMBDA (const int i)
    		{
                h_Output[i+nDir] = 0.0;
            });          
            Kokkos::deep_copy(Output,h_Output);            

            printf("%s\n", "finished data copying");
            //--------------------------------------------------------------------------

            // set size of heap memory on device
            /*#ifdef __CUDA_ARCH__
                int nqtot   = nquad0*nquad1;
                int max1 = (nqtot >= ncoeffs) ? nqtot : ncoeffs;
                int max2 = (nquad1*nmodes0 >= nquad0*nmodes1) ? nquad1*nmodes0 : nquad0*nmodes1;
                int wspsize = (max1 >= max2) ? max1 : max2;
                int el_size = (2*ncoeffs + 4* wspsize);
                size_t size;
                cudaDeviceGetLimit(&size, cudaLimitMallocHeapSize);
                printf("Heap size found to be %d\n",(int)size);
                //cudaDeviceSetLimit(cudaLimitMallocHeapSize, 50 * elmts * el_size * sizeof(double));                
                cudaDeviceGetLimit(&size, cudaLimitMallocHeapSize);
                printf("New heap size set to be %d\n",(int)size);
            #endif*/

            // evaluate initial residual error for exit check            
            eps = 0.0;
            Kokkos::parallel_reduce(range_policy(0,nNonDir),KOKKOS_LAMBDA(const int i, NekDouble &ieps)
            {
                ieps += r_A[i] * r_A[i] * map[i+nDir];
            },eps);
            

            rhs_magnitude = 0.0;
            /*Kokkos::parallel_reduce(range_policy_host(0,nGlobal),KOKKOS_LAMBDA(const int i, NekDouble &mag)
            {
                mag += pInput[i] * pInput[i] * h_map[i];
            },rhs_magnitude);*/
            for (int i = 0; i < nGlobal; ++i)
            {
                rhs_magnitude += pInput[i] * pInput[i] * h_map[i];
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
            Kokkos::parallel_for(range_policy(0,nNonDir), KOKKOS_LAMBDA (const int i)
    		{
                w_A[i+nDir] = r_A[i] * diagonals[i];
            });

            GeneralMatrixOp_Kokkos(
                    w_A, s_A, lambda,
                    quadMetricGlo, laplacian00Glo, laplacian01Glo, laplacian11Glo,
                    nquad0, nquad1, nmodes0, nmodes1, ncoeffs, 
                    coeff_offset, elmts,
                    base0, base1, dbase0, dbase1,
                    D0, D1,
                    numLocalCoeffs, numGlobalCoeffs,
                    localToGlobalMap, localToGlobalSign,totalIterations,
                    coloursetArray, cs_sizes, ncs,
                    base0_len, base1_len);

            rho = 0.0;
            Kokkos::parallel_reduce(range_policy(0,nNonDir),KOKKOS_LAMBDA(const int i, NekDouble &irho)
            {
                irho += r_A[i] * w_A[i+nDir] * map[i+nDir];
            },rho);

            mu = 0.0;
            Kokkos::parallel_reduce(range_policy(0,nNonDir),KOKKOS_LAMBDA(const int i, NekDouble &imu)
            {
                imu += s_A[i+nDir] * w_A[i+nDir] * map[i+nDir];
            },mu);

            min_resid         = rhs_magnitude;
            beta              = 0.0;
            alpha             = rho/mu;            

            // Continue until convergence
            k = 0;
            totalIterations = 1;
            while (true)
            {
                // CUDA Profiling
                //#ifdef __CUDA_ARCH__
                if (k == 0)
                {
                    cuProfilerStart();
                }
                else if(k == 3)
                {
                    cudaDeviceSynchronize();
                    cuProfilerStop();
                }
                //#endif

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

                Kokkos::parallel_for(range_policy(0,nNonDir), KOKKOS_LAMBDA (const int i)
    		    {
                    // Compute new search direction p_k, q_k
                    p_A[i] = beta * p_A[i] + w_A[i+nDir]; 
                    q_A[i] = beta * q_A[i] + s_A[i+nDir];
                    // Update solution x_{k+1}
                    Output[i+nDir] = alpha * p_A[i] + Output[i+nDir];
                    // Update residual vector r_{k+1}
                    r_A[i] = -alpha * q_A[i] + r_A[i];
                    // Apply preconditioner  
                    w_A[i+nDir] = r_A[i] * diagonals[i];                                      
                });               
                
                // Perform the method-specific matrix-vector multiply operation.
                printf("CG iteration %i ==================================\n", totalIterations);
                GeneralMatrixOp_Kokkos(
                    w_A, s_A, lambda,
                    quadMetricGlo, laplacian00Glo, laplacian01Glo, laplacian11Glo,
                    nquad0, nquad1, nmodes0, nmodes1, ncoeffs, 
                    coeff_offset, elmts,
                    base0, base1, dbase0, dbase1,
                    D0, D1,
                    numLocalCoeffs, numGlobalCoeffs,
                    localToGlobalMap, localToGlobalSign,totalIterations,
                    coloursetArray, cs_sizes, ncs,
                    base0_len, base1_len);               

                rho_new = 0.0;
                Kokkos::parallel_reduce(range_policy(0,nNonDir),KOKKOS_LAMBDA(const int i, NekDouble &irho_new)
            	{
                    irho_new += r_A[i] * w_A[i+nDir] * map[i+nDir];
                },rho_new);

                mu = 0.0;
                Kokkos::parallel_reduce(range_policy(0,nNonDir),KOKKOS_LAMBDA(const int i, NekDouble &imu)
            	{
                    imu += s_A[i+nDir] * w_A[i+nDir] * map[i+nDir];
                },mu);

                eps = 0.0;
                Kokkos::parallel_reduce(range_policy(0,nNonDir),KOKKOS_LAMBDA(const int i, NekDouble &ieps)
                {
                    ieps += r_A[i] * r_A[i] * map[i+nDir];
                },eps);

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
            
			Kokkos::deep_copy(h_Output,Output);
			for (int i = 0; i < nNonDir; ++i)
            {
                pOutput[i+nDir] = h_Output[i+nDir];
            } 
            Kokkos::finalize();
            printf("%s\n", "end of conjugate gradient Kokkos");
            
    	}


        void GlobalLinSysIterative::GeneralMatrixOp_Kokkos(
                const Kokkos::View<double*> inarray,
                Kokkos::View<double*> outarray,
                const Kokkos::View<double[1]> lambda,                
                const Kokkos::View<double**> quadMetricGlo,
                const Kokkos::View<double**> laplacian00Glo,
                const Kokkos::View<double**> laplacian01Glo,
                const Kokkos::View<double**> laplacian11Glo,
                const int &nquad0, const int &nquad1,
                const int &nmodes0, const int &nmodes1, const int &ncoeffs, 
                const Kokkos::View<int*> coeff_offset,
                const int &elmts,
                const Kokkos::View<double*> base0,
                const Kokkos::View<double*> base1,
                const Kokkos::View<double*> dbase0,
                const Kokkos::View<double*> dbase1,
                const Kokkos::View<double*> D0,
                const Kokkos::View<double*> D1,
                const int &numLocalCoeffs, const int &numGlobalCoeffs,
                const Kokkos::View<int*> localToGlobalMap,
                const Kokkos::View<double*> localToGlobalSign,
                const int iteration,
                const Kokkos::View<int**> coloursetArray, int cs_sizes[], int ncs,
                const int base0_len, const int base1_len)
        {
            //printf("%s %i\n", "do the global to local mapping", numLocalCoeffs);
            Kokkos::View<double*> transfer_in("transfer_in", numLocalCoeffs);
            Kokkos::View<double*> transfer_out("transfer_out",numLocalCoeffs);
             
            Kokkos::parallel_for(range_policy(0,numLocalCoeffs), KOKKOS_LAMBDA (const int i)    		
            {
                transfer_in[i] = localToGlobalSign[i] * inarray[localToGlobalMap[i]];
            });
            
            GeneralMatrixOp_IterPerExp_Kokkos(transfer_in,transfer_out,lambda,
                    quadMetricGlo, laplacian00Glo,laplacian01Glo,laplacian11Glo,                    
                    nquad0, nquad1, nmodes0, nmodes1, ncoeffs,
                    coeff_offset, elmts,
                    base0, base1, dbase0, dbase1, D0, D1,
                    base0_len, base1_len);
            
            //printf("\n%s\n", "do the local to global mapping");
            Kokkos::parallel_for(range_policy(0,numGlobalCoeffs), KOKKOS_LAMBDA (const int i)    		
            {
                outarray[i] = 0.0;
            });
            
            // do mapping on per element basis and using element colourgroups
            // ( numLocalCoeffs = elmts * ncoeffs, elmts = sum(cs_sizes) )
            // need to use element-colouring to prevent race condition
            for (int cs = 0; cs < ncs; ++cs)
            {
                Kokkos::parallel_for( team_policy( cs_sizes[cs] , ncoeffs )
                , KOKKOS_LAMBDA ( const member_type& teamMember)
                {
                    const int i = teamMember.league_rank();
                    {
                        const int el = coloursetArray(cs,i);
                        Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ncoeffs ), [&] ( const int j)                         
                        { 
                            const int n = el * ncoeffs + j;
                            outarray[localToGlobalMap[n]] += localToGlobalSign[n] * transfer_out[n];
                        });
                    }
                });
            }           
            printf("\n");            
        }


        void GlobalLinSysIterative::GeneralMatrixOp_IterPerExp_Kokkos(
                const Kokkos::View<double*> transfer_in,
                Kokkos::View<double*> transfer_out,
                const Kokkos::View<double[1]> lambda,
                const Kokkos::View<double**> quadMetricGlo,
                const Kokkos::View<double**> laplacian00Glo,
                const Kokkos::View<double**> laplacian01Glo,
                const Kokkos::View<double**> laplacian11Glo,
                const int &nquad0, const int &nquad1,
                const int &nmodes0, const int &nmodes1, const int &ncoeffs, 
                const Kokkos::View<int*> coeff_offset,
                const int &elmts,
                const Kokkos::View<double*> base0,
                const Kokkos::View<double*> base1,
                const Kokkos::View<double*> dbase0,
                const Kokkos::View<double*> dbase1,
                const Kokkos::View<double*> D0,
                const Kokkos::View<double*> D1,
                const int base0_len, const int base1_len)
        {
            int nqtot   = nquad0*nquad1;
            int max1 = (nqtot >= ncoeffs) ? nqtot : ncoeffs;
            int max2 = (nquad1*nmodes0 >= nquad0*nmodes1) ? nquad1*nmodes0 : nquad0*nmodes1;
            int wspsize = (max1 >= max2) ? max1 : max2;

            printf("%s %i\n", "perform operations by element, elements in total: ", elmts);
            int max_threads = 32;
            int no_teams = (elmts + max_threads - 1) / max_threads;
            //int scratch_size_thread = ScratchViewType::shmem_size(4*wspsize);
            int scratch_size_team0 = ScratchViewType::shmem_size(
                2 * base0_len + 2 * base1_len + nquad0*nquad0 + nquad1*nquad1);
            int scratch_size_team1 = ScratchViewType::shmem_size((2*ncoeffs+4*wspsize)*max_threads);

            //Kokkos::parallel_for(range_policy(0,elmts),KOKKOS_LAMBDA (const int el)   {    
            Kokkos::parallel_for( team_policy( no_teams , max_threads )
                //.set_scratch_size(1,Kokkos::PerThread(scratch_size_thread))
                .set_scratch_size(0,Kokkos::PerTeam(scratch_size_team0)) // fast shared memory
                .set_scratch_size(1,Kokkos::PerTeam(scratch_size_team1)) // L2 cache?
                , KOKKOS_LAMBDA ( const member_type& teamMember)
            {
                const int el_o = teamMember.league_rank();

                // load per team scratch views
                ScratchViewType s_base0(teamMember.team_scratch(0),base0_len);
                ScratchViewType s_base1(teamMember.team_scratch(0),base1_len);
                ScratchViewType s_dbase0(teamMember.team_scratch(0),base0_len);
                ScratchViewType s_dbase1(teamMember.team_scratch(0),base1_len);
                ScratchViewType s_D0(teamMember.team_scratch(0),nquad0*nquad0);
                ScratchViewType s_D1(teamMember.team_scratch(0),nquad1*nquad1);
                Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , base0_len ), [&] ( const int i)
                {
                    s_base0[i] = base0[i];
                });
                Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , base1_len ), [&] ( const int i)
                {
                    s_base1[i] = base1[i];
                });
                Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , base0_len ), [&] ( const int i)
                {
                    s_dbase0[i] = dbase0[i];
                });
                Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , base1_len ), [&] ( const int i)
                {
                    s_dbase1[i] = dbase1[i];
                });
                Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , nquad0*nquad0 ), [&] ( const int i)
                {
                    s_D0[i] = D0[i];
                });
                Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , nquad1*nquad1 ), [&] ( const int i)
                {
                    s_D1[i] = D1[i];
                });

                int level = 1;
                ScratchViewType32 s_outarray(teamMember.team_scratch(level),ncoeffs);
                ScratchViewType32 s_inarray(teamMember.team_scratch(level),ncoeffs);
                Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , max_threads ), [&] ( const int el_i)
                {
                    const int el = el_o * max_threads + el_i;
                    if (el < elmts)
                    {
                        for (int i = 0; i < ncoeffs; ++i)
                        {
                            s_inarray(i,el_i) = transfer_in[el*ncoeffs+i];
                        }
                    }
                });
                
                ScratchViewType32 s_wsp0(teamMember.team_scratch(level),wspsize);
                ScratchViewType32 s_wsp1(teamMember.team_scratch(level),wspsize);
                ScratchViewType32 s_wsp2(teamMember.team_scratch(level),wspsize);
                ScratchViewType32 s_wsp1L(teamMember.team_scratch(level),wspsize);

                teamMember.team_barrier();
                
                Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , max_threads ), [&] ( const int el_i)
                {
                    const int el = el_o * max_threads + el_i;
                    if (el < elmts)
                    {
                        HelmholtzMatrixOp_MatFree_Kokkos(
                            s_inarray,
                            s_outarray,
                            s_wsp0, s_wsp1, s_wsp2, s_wsp1L,
                            el, el_i, max_threads,
                            coeff_offset,
                            lambda,
                             quadMetricGlo,
                            laplacian00Glo,
                            laplacian01Glo,
                            laplacian11Glo,
                            nquad0, nquad1, nmodes0, nmodes1, ncoeffs,
                            s_base0, s_base1, s_dbase0, s_dbase1, s_D0, s_D1,
                            teamMember, wspsize);
                    }
                });

                Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , max_threads ), [&] ( const int el_i)
                {
                    const int el = el_o * max_threads + el_i;
                    if (el < elmts)
                    {
                        for (int i = 0; i < ncoeffs; ++i)
                        {
                            transfer_out[el*ncoeffs + i] = s_outarray(i,el_i);
                        }
                    }
                });

            });
        }

        KOKKOS_INLINE_FUNCTION
        void GlobalLinSysIterative::HelmholtzMatrixOp_MatFree_Kokkos(
                const ScratchViewType32 s_inarray,
                      ScratchViewType32 s_outarray,
                      ScratchViewType32 s32_wsp0,
                      ScratchViewType32 s32_wsp1,
                      ScratchViewType32 s32_wsp2,
                      ScratchViewType32 s32_wsp1L,
                const int &el, const int el_i, const int max_threads,
                const Kokkos::View<int*>  coeff_offset,
                const Kokkos::View<double[1]> lambda,
				const Kokkos::View<double**>  quadMetricGlo,
                const Kokkos::View<double**> laplacian00Glo,
                const Kokkos::View<double**> laplacian01Glo,
                const Kokkos::View<double**> laplacian11Glo,
                const int &nquad0, const int &nquad1,
                const int &nmodes0, const int &nmodes1, const int &ncoeffs,
                const ScratchViewType s_base0,
                const ScratchViewType s_base1,
                const ScratchViewType s_dbase0,
                const ScratchViewType s_dbase1,
                const ScratchViewType s_D0,
                const ScratchViewType s_D1,
                const member_type &teamMember, const int &wspsize)
        {
            //reordering possible to use less memory?

            int nqtot   = nquad0*nquad1;            
            
            BwdTrans_SumFacKernel_Kokkos_s32(s_base0, s_base1, s_inarray, s32_wsp0, s32_wsp2,
                nmodes0, nmodes1, nquad0, nquad1, el_i, max_threads);
            
            for (int i = 0; i < nqtot; ++i)
            {
                s32_wsp1(i,el_i) = quadMetricGlo(i,el) * s32_wsp0(i,el_i);
            }
            
            IProductWRTBase_SumFacKernel_Kokkos_s32(s_base0, s_base1, s32_wsp1, s_outarray,
                                         s32_wsp2, nmodes0, nmodes1, nquad0, nquad1, el_i, max_threads);                 	
                  	
            // === LaplacianMatrixOp_MatFree_Kernel ===
            //ScratchViewType32 s32_wsp2L = s32_wsp2;

            PhysTensorDeriv_Kokkos_s32(s32_wsp0,s32_wsp1L,s32_wsp2, nquad0, nquad1, s_D0, s_D1, el_i, max_threads);

            //ScratchViewType32 s32_wsp0L = s32_wsp0;          
        
            for (int i = 0; i < nqtot; ++i)
            {
                s32_wsp0(i,el_i) = laplacian00Glo(i,el) * s32_wsp1L(i,el_i) 
                                 + laplacian01Glo(i,el) * s32_wsp2(i,el_i);
            }
            for (int i = 0; i < nqtot; ++i)
            {
                s32_wsp2(i,el_i) = laplacian01Glo(i,el) * s32_wsp1L(i,el_i)
                                 + laplacian11Glo(i,el) * s32_wsp2(i,el_i);
            }
            IProductWRTBase_SumFacKernel_Kokkos_s32(s_dbase0, s_base1,s32_wsp0,s32_wsp1 ,s32_wsp1L,
                         nmodes0, nmodes1, nquad0, nquad1, el_i, max_threads);
            IProductWRTBase_SumFacKernel_Kokkos_s32( s_base0,s_dbase1,s32_wsp2,s32_wsp1L,s32_wsp0,
                         nmodes0, nmodes1, nquad0, nquad1, el_i, max_threads);
            
            for (int i = 0; i < ncoeffs; ++i)
            {
                s32_wsp1(i,el_i) += s32_wsp1L(i,el_i);
            }
            // === end LaplacianMatrixOp_MatFree_Kernel ===

            for (int i = 0; i < ncoeffs; ++i)
            {
                s_outarray(i,el_i) = lambda[0] * s_outarray(i,el_i) + s32_wsp1(i,el_i);
            }
        }

        KOKKOS_INLINE_FUNCTION
        void GlobalLinSysIterative::IProductWRTBase_SumFacKernel_Kokkos(
                const ScratchViewType base0,
                const ScratchViewType base1,
                const ScratchViewType inarray,
                ScratchViewType outarray,
                ScratchViewType wsp,
                const int &nmodes0, const int &nmodes1,
                const int &nquad0, const int &nquad1,
                const int el_i, const int max_threads)
        {               
            const double alpha = 1.0;
            const double beta = 0.0;
            plainerDgemm('T','N',nquad1,nmodes0,nquad0,alpha,inarray.ptr_on_device(),nquad0,
                        base0.ptr_on_device(),nquad0,beta,wsp.ptr_on_device(),nquad1);
            int i, mode;
            for (mode=i=0; i < nmodes0; ++i)
            {
                plainDgemv('T',nquad1,nmodes1-i,alpha, base1.ptr_on_device()+mode*nquad1,
                            nquad1,wsp.ptr_on_device()+i*nquad1,1, beta,
                            outarray.ptr_on_device() + mode,1);
                mode += nmodes1 - i;
            }
            double result = plainDdot(nquad1,base1.ptr_on_device()+nquad1,1,
                                          wsp.ptr_on_device()+nquad1,1);
            outarray[1] += result;            
        }

        KOKKOS_INLINE_FUNCTION
        void GlobalLinSysIterative::IProductWRTBase_SumFacKernel_Kokkos_s32(
                const ScratchViewType base0,
                const ScratchViewType base1,
                const ScratchViewType32 inarray32,
                ScratchViewType32 outarray32,
                ScratchViewType32 wsp32,
                const int &nmodes0, const int &nmodes1,
                const int &nquad0, const int &nquad1,
                const int el_i, const int max_threads)
        {               
            const double alpha = 1.0;
            const double beta = 0.0;
            stridedDgemm('T','N',nquad1,nmodes0,nquad0,
                alpha,
                inarray32.ptr_on_device(),nquad0,max_threads,el_i,
                base0.ptr_on_device(),nquad0, 1,0,
                beta,
                wsp32.ptr_on_device(),nquad1,max_threads,el_i);
            
            int i, mode;
            for (mode=i=0; i < nmodes0; ++i)
            {
                plainDgemv('T',nquad1,nmodes1-i,alpha, base1.ptr_on_device()+mode*nquad1,
                            nquad1,wsp32.ptr_on_device()+i*nquad1*max_threads+el_i,max_threads, beta,
                            outarray32.ptr_on_device() + mode*max_threads+el_i,max_threads);
                mode += nmodes1 - i;
            }
            double result = plainDdot(nquad1,base1.ptr_on_device()+nquad1,1,
                                          wsp32.ptr_on_device()+nquad1*max_threads+el_i,max_threads);
            outarray32(1,el_i) += result;            
        }



        KOKKOS_INLINE_FUNCTION
        void GlobalLinSysIterative::BwdTrans_SumFacKernel_Kokkos(
                const ScratchViewType base0,
                const ScratchViewType base1,
                const ScratchViewType inarray,
                ScratchViewType outarray,
                ScratchViewType wsp,
                const int &nmodes0, const int &nmodes1,
                const int &nquad0, const int &nquad1,
                const int el_i, const int max_threads)
        {
            const double alpha = 1.0;
            const double beta = 0.0;
            int i, mode;
            for (i = mode = 0; i < nmodes0; ++i)
            {
                plainDgemv('N', nquad1,nmodes1-i,alpha,base1.ptr_on_device()+mode*nquad1,
                            nquad1,inarray.ptr_on_device()+mode,1,beta,wsp.ptr_on_device()+i*nquad1,1);
                mode += nmodes1-i;
            }
            plainDaxpy(nquad1,inarray[1],base1.ptr_on_device()+nquad1,1,
                            wsp.ptr_on_device()+nquad1,1);
            
            plainerDgemm('N','T', nquad0,nquad1,nmodes0,alpha, base0.ptr_on_device(),nquad0,
                        wsp.ptr_on_device(), nquad1,beta,outarray.ptr_on_device(), nquad0);
        }

        KOKKOS_INLINE_FUNCTION
        void GlobalLinSysIterative::BwdTrans_SumFacKernel_Kokkos_s32(
                const ScratchViewType base0,
                const ScratchViewType base1,
                const ScratchViewType32 inarray32,
                ScratchViewType32 outarray32,
                ScratchViewType32 wsp32,
                const int &nmodes0, const int &nmodes1,
                const int &nquad0, const int &nquad1,
                const int el_i, const int max_threads)
        {
            const double alpha = 1.0;
            const double beta = 0.0;
            int i, mode;
            for (i = mode = 0; i < nmodes0; ++i)
            {
                plainDgemv('N', nquad1,nmodes1-i,alpha,base1.ptr_on_device()+mode*nquad1,
                            nquad1,inarray32.ptr_on_device()+mode*max_threads+el_i,max_threads,beta,wsp32.ptr_on_device()+i*nquad1*max_threads+el_i,max_threads);
                mode += nmodes1-i;
            }
            plainDaxpy(nquad1,inarray32(1,el_i),base1.ptr_on_device()+nquad1,1,
                            wsp32.ptr_on_device()+nquad1*max_threads+el_i,max_threads);
            
            stridedDgemm('N','T', nquad0,nquad1,nmodes0,
                alpha,
                base0.ptr_on_device(),nquad0,1,0,
                wsp32.ptr_on_device(), nquad1,max_threads,el_i,
                beta,
                outarray32.ptr_on_device(), nquad0,max_threads,el_i);            
        }


        KOKKOS_INLINE_FUNCTION
        void GlobalLinSysIterative::PhysTensorDeriv_Kokkos(
                const ScratchViewType inarray,
                ScratchViewType outarray_d0,
                ScratchViewType outarray_d1,
                const int &nquad0, const int &nquad1,
                const ScratchViewType D0,
                const ScratchViewType D1,
                const int el_i, const int max_threads)
        {
            const double alpha = 1.0;
            const double beta = 0.0;
            plainerDgemm('N', 'N', nquad0, nquad1, nquad0, alpha,
                    (D0.ptr_on_device()), nquad0, inarray.ptr_on_device(), nquad0, beta,
                    outarray_d0.ptr_on_device(), nquad0);
            plainerDgemm('N', 'T', nquad0, nquad1, nquad1, alpha, inarray.ptr_on_device(), nquad0,
                     (D1.ptr_on_device()), nquad1, beta, outarray_d1.ptr_on_device(), nquad0);
        }




        KOKKOS_INLINE_FUNCTION
        void GlobalLinSysIterative::PhysTensorDeriv_Kokkos_s32(
                const ScratchViewType32 inarray,
                ScratchViewType32 outarray_d0,
                ScratchViewType32 outarray_d1,
                const int &nquad0, const int &nquad1,
                const ScratchViewType D0,
                const ScratchViewType D1,
                const int el_i, const int max_threads)
        {
            const double alpha = 1.0;
            const double beta = 0.0;
            
            stridedDgemm('N', 'N', nquad0, nquad1, nquad0,
                alpha,
                D0.ptr_on_device(), nquad0,1,0,
                inarray.ptr_on_device(), nquad0,max_threads,el_i,
                beta,
                outarray_d0.ptr_on_device(), nquad0,max_threads,el_i);

            stridedDgemm('N', 'T', nquad0, nquad1, nquad1,
                alpha,
                inarray.ptr_on_device(), nquad0,max_threads,el_i,
                D1.ptr_on_device(), nquad1,1,0,
                beta,
                outarray_d1.ptr_on_device(), nquad0,max_threads,el_i);
        }

    }
}

#endif