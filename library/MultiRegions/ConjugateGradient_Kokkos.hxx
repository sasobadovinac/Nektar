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

#include <MultiRegions/GlobalLinSysIterative.h>

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
    		
            //Initialise Kokkos
            Kokkos::InitArguments args;
            args.num_threads = 1;
            Kokkos::initialize(args);   
            

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
            Kokkos::View<double[1]> lambda ("lambda");                        
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
            Kokkos::View<double*> base0 ("base0", nquad0);                        
            typename Kokkos::View< double*>::HostMirror h_base0;
            h_base0 = Kokkos::create_mirror_view(base0);
            Kokkos::View<double*> dbase0 ("dbase0", nquad0);                        
            typename Kokkos::View< double*>::HostMirror h_dbase0;
            h_dbase0 = Kokkos::create_mirror_view(dbase0);
            Kokkos::parallel_for(range_policy_host(0,nquad0), KOKKOS_LAMBDA (const int i)
    		{
            	h_base0(i) = Nekbase0[i];
            	h_dbase0(i) = Nekdbase0[i];
            });
            Kokkos::deep_copy(base0,h_base0);
            Kokkos::deep_copy(dbase0,h_dbase0);

            // copying base1 and dbase1
            int base1_len = nmodes0*(nmodes1-0.5*(nmodes0-1))*nquad1;
            Kokkos::View<double*> base1 ("base1", base1_len);                        
            typename Kokkos::View< double*>::HostMirror h_base1;
            h_base1 = Kokkos::create_mirror_view(base1);
            Kokkos::View<double*> dbase1 ("dbase1", base1_len);                        
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
            Kokkos::View<int*> coeff_offset ("coeff_offset", ncoeffs);                        
            typename Kokkos::View< int*>::HostMirror h_coeff_offset;
            h_coeff_offset = Kokkos::create_mirror_view(coeff_offset);
            Kokkos::parallel_for(range_policy_host(0,ncoeffs), KOKKOS_LAMBDA (const int i)
    		{
            	h_coeff_offset(i) = Nekcoeff_offset[i];
            });
            Kokkos::deep_copy(coeff_offset,h_coeff_offset);

            // copying DNekMatSharedPtr NekD0
            Kokkos::View<double*> D0 ("D0", nquad0);                        
            typename Kokkos::View< double*>::HostMirror h_D0;
            h_D0 = Kokkos::create_mirror_view(D0);
            NekDouble* D0_raw = &(NekD0->GetPtr())[0];
            Kokkos::parallel_for(range_policy_host(0,nquad0), KOKKOS_LAMBDA (const int i)
    		{
            	h_D0(i) = *(D0_raw+i);
            });
            Kokkos::deep_copy(D0,h_D0);
            // copying DNekMatSharedPtr NekD1
            Kokkos::View<double*> D1 ("D1", nquad1);                        
            typename Kokkos::View< double*>::HostMirror h_D1;
            h_D1 = Kokkos::create_mirror_view(D1);
            NekDouble* D1_raw = &(NekD1->GetPtr())[0];
            Kokkos::parallel_for(range_policy_host(0,nquad1), KOKKOS_LAMBDA (const int i)
    		{
            	h_D1(i) = *(D1_raw+i);
            });
            Kokkos::deep_copy(D1,h_D1);



            // --------------------------------------------------------------------------------
            // Gathering  quadrature and laplacian metrics for each element
            int metricSize = elmts * nquad0 * nquad1;
            Array<OneD, NekDouble> NekquadMetricGlo (4*metricSize);
            Array<OneD, NekDouble> Neklaplacian00Glo(NekquadMetricGlo+metricSize);
            Array<OneD, NekDouble> Neklaplacian01Glo(NekquadMetricGlo+2*metricSize);
            Array<OneD, NekDouble> Neklaplacian11Glo(NekquadMetricGlo+3*metricSize);

            int numLocalCoeffs, numGlobalCoeffs;
            Array<OneD, const int> NeklocalToGlobalMap; //numLocalCoeffs
            Array<OneD, const NekDouble> NeklocalToGlobalSign; //numLocalCoeffs
            GetMatrixMultiplyMetrics(
                        NekquadMetricGlo, Neklaplacian00Glo, Neklaplacian01Glo, Neklaplacian11Glo,
                        nquad0, nquad1, elmts,
                        numLocalCoeffs, numGlobalCoeffs,
                        NeklocalToGlobalMap, NeklocalToGlobalSign);

            // copying metrics
            Kokkos::View<double*> quadMetricGlo ("quadMetricGlo", metricSize);                        
            typename Kokkos::View< double*>::HostMirror h_quadMetricGlo;
            h_quadMetricGlo = Kokkos::create_mirror_view(quadMetricGlo);
            Kokkos::View<double*> laplacian00Glo ("laplacian00Glo", metricSize);                        
            typename Kokkos::View< double*>::HostMirror h_laplacian00Glo;
            h_laplacian00Glo = Kokkos::create_mirror_view(laplacian00Glo);
            Kokkos::View<double*> laplacian01Glo ("laplacian01Glo", metricSize);                        
            typename Kokkos::View< double*>::HostMirror h_laplacian01Glo;
            h_laplacian01Glo = Kokkos::create_mirror_view(laplacian01Glo);
            Kokkos::View<double*> laplacian11Glo ("laplacian11Glo", metricSize);                        
            typename Kokkos::View< double*>::HostMirror h_laplacian11Glo;
            h_laplacian11Glo = Kokkos::create_mirror_view(laplacian11Glo);            
            Kokkos::parallel_for(range_policy_host(0,metricSize), KOKKOS_LAMBDA (const int i)
    		{
            	h_quadMetricGlo(i) = NekquadMetricGlo[i];
            	h_laplacian00Glo(i) = Neklaplacian00Glo[i];
            	h_laplacian01Glo(i) = Neklaplacian01Glo[i];
            	h_laplacian11Glo(i) = Neklaplacian11Glo[i];
            });
            Kokkos::deep_copy(quadMetricGlo,h_quadMetricGlo);
            Kokkos::deep_copy(laplacian00Glo,h_laplacian00Glo);
            Kokkos::deep_copy(laplacian01Glo,h_laplacian01Glo);
            Kokkos::deep_copy(laplacian11Glo,h_laplacian11Glo);
            
            //Array<OneD, const int> NeklocalToGlobalMap; //numLocalCoeffs
            //Array<OneD, const NekDouble> NeklocalToGlobalSign; //numLocalCoeffs
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

            printf("%s\n", "finished data gathering");
            //--------------------------------------------------------------------------




            // CG variables
            NekDouble alpha, beta, rho, rho_new, mu, eps, min_resid, rhs_magnitude;
            int totalIterations, k;

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


            // evaluate initial residual error for exit check            
            eps = 0.0;
            Kokkos::parallel_reduce(range_policy(0,nNonDir),KOKKOS_LAMBDA(const int i, NekDouble &ieps)
            {
                ieps += r_A[i] * r_A[i] * map[i+nDir];
            },eps);
            
            rhs_magnitude = 0.0;
            Kokkos::parallel_reduce(range_policy(0,nGlobal),KOKKOS_LAMBDA(const int i, NekDouble &mag)
            {
                mag += pInput[i] * pInput[i] * map[i];
            },rhs_magnitude);                      

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
                    localToGlobalMap, localToGlobalSign);

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
                Kokkos::parallel_for(range_policy(0,nNonDir), KOKKOS_LAMBDA (const int i)
    		    {
                    p_A[i] = beta * p_A[i] + w_A[i+nDir];                    
                });               
                
                //Vmath::Svtvp(nNonDir, beta, &q_A[0], 1, &s_A[nDir], 1, &q_A[0], 1);
                Kokkos::parallel_for(range_policy(0,nNonDir), KOKKOS_LAMBDA (const int i)
    			{
                    q_A[i] = beta * q_A[i] + s_A[i+nDir];                    
                });
                // Update solution x_{k+1}
                //Vmath::Svtvp(nNonDir, alpha, &p_A[0], 1, &pOutput[nDir], 1, &pOutput[nDir], 1);
                Kokkos::parallel_for(range_policy(0,nNonDir), KOKKOS_LAMBDA (const int i)    		
                {
                    Output[i+nDir] = alpha * p_A[i] + Output[i+nDir];                    
                }); 
                // Update residual vector r_{k+1}
                //Vmath::Svtvp(nNonDir, -alpha, &q_A[0], 1, &r_A[0], 1, &r_A[0], 1);
                Kokkos::parallel_for(range_policy(0,nNonDir), KOKKOS_LAMBDA (const int i)    		
                {
                    r_A[i] = -alpha * q_A[i] + r_A[i];                    
                });

                // Apply preconditioner
                //m_precon->DoPreconditioner(r_A, tmp = w_A + nDir);
                //Vmath::Vmul(nNonDir, &r_A[0], 1, &diagonals[0], 1, &w_A[nDir], 1);
                Kokkos::parallel_for(range_policy(0,nNonDir), KOKKOS_LAMBDA (const int i)    		
                {
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
                    localToGlobalMap, localToGlobalSign);               

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
                //const Array<OneD,const NekDouble> &inarray,
                //Array<OneD,      NekDouble> &outarray,
                Kokkos::View<double*> inarray,
                Kokkos::View<double*> outarray,
                //    const NekDouble &lambda,
                Kokkos::View<double[1]> lambda,                
                //const Array<OneD, const NekDouble> &quadMetricGlo,                
                //const Array<OneD, const NekDouble> &laplacian00Glo,
                //const Array<OneD, const NekDouble> &laplacian01Glo,
                //const Array<OneD, const NekDouble> &laplacian11Glo,
                Kokkos::View<double*> quadMetricGlo,
                Kokkos::View<double*> laplacian00Glo,
                Kokkos::View<double*> laplacian01Glo,
                Kokkos::View<double*> laplacian11Glo,
                const int &nquad0, const int &nquad1,
                const int &nmodes0, const int &nmodes1, const int &ncoeffs, 
                //const Array<OneD, const int>  &coeff_offset, 
                Kokkos::View<int*> coeff_offset,
                const int &elmts,
                //const Array<OneD, const NekDouble> &base0,
                //const Array<OneD, const NekDouble> &base1,
                //const Array<OneD, const NekDouble> &dbase0,
                //const Array<OneD, const NekDouble> &dbase1,
                Kokkos::View<double*> base0,
                Kokkos::View<double*> base1,
                Kokkos::View<double*> dbase0,
                Kokkos::View<double*> dbase1,
                //const DNekMatSharedPtr &D0, const DNekMatSharedPtr &D1,
                Kokkos::View<double*> D0,
                Kokkos::View<double*> D1,
                const int &numLocalCoeffs, const int &numGlobalCoeffs,
                //const Array<OneD, const int> &localToGlobalMap,
                //const Array<OneD, const NekDouble> &localToGlobalSign,
                Kokkos::View<int*> localToGlobalMap,
                Kokkos::View<double*> localToGlobalSign
                )
        {
            printf("%s\n", "do the global to local mapping");
            Kokkos::View<double*> transfer_in("transfer_in", numLocalCoeffs);
            Kokkos::View<double*> transfer_out("transfer_out", elmts*ncoeffs);

            //GlobalToLocal_plain(inarray,tmp1);
            //Vmath::Gathr(numLocalCoeffs, localToGlobalSign.get(),
            //         inarray.get(), localToGlobalMap.get(), tmp1.get());
            Kokkos::parallel_for(range_policy(0,numLocalCoeffs), KOKKOS_LAMBDA (const int i)    		
            {
                transfer_in[i] = localToGlobalSign[i] * inarray[localToGlobalMap[i]];            
            });
            
            GeneralMatrixOp_IterPerExp_Kokkos(transfer_in,transfer_out,lambda,
                    quadMetricGlo, laplacian00Glo,laplacian01Glo,laplacian11Glo,                    
                    nquad0, nquad1, nmodes0, nmodes1, ncoeffs,
                    coeff_offset, elmts,
                    base0, base1, dbase0, dbase1, D0, D1);
            //double *out_raw = transfer_out.ptr_on_device();
            
            printf("%s\n", "do the local to global mapping");
            //Assemble_plain(tmp2,outarray);  
            //Vmath::Zero(numGlobalCoeffs, outarray.get(), 1);
            Kokkos::parallel_for(range_policy(0,numGlobalCoeffs), KOKKOS_LAMBDA (const int i)    		
            {
                outarray[i] = 0.0;
            });
            //Vmath::Assmb(numLocalCoeffs, localToGlobalSign.get(), 
            //            tmp2.get(), localToGlobalMap.get(), outarray.get());
            Kokkos::parallel_for(range_policy(0,numLocalCoeffs), KOKKOS_LAMBDA (const int i)    		
            {
                outarray[localToGlobalMap[i]] += localToGlobalSign[i] * transfer_out[i]; 
            });
            
        }


        void GlobalLinSysIterative::GeneralMatrixOp_IterPerExp_Kokkos(
                Kokkos::View<double*> transfer_in,
                Kokkos::View<double*> transfer_out,
                Kokkos::View<double[1]> lambda,
                Kokkos::View<double*> quadMetricGlo,
                Kokkos::View<double*> laplacian00Glo,
                Kokkos::View<double*> laplacian01Glo,
                Kokkos::View<double*> laplacian11Glo,
                const int &nquad0, const int &nquad1,
                const int &nmodes0, const int &nmodes1, const int &ncoeffs, 
                Kokkos::View<int*> coeff_offset,
                const int &elmts,
                Kokkos::View<double*> base0,
                Kokkos::View<double*> base1,
                Kokkos::View<double*> dbase0,
                Kokkos::View<double*> dbase1,
                Kokkos::View<double*> D0,
                Kokkos::View<double*> D1)
        {
            printf("%s\n", "perform operations by element");            
            // Calculating
            //Kokkos::parallel_for(range_policy(0,elmts),KOKKOS_LAMBDA (const int el)
            for (int el = 0; el < elmts; ++el)
            {
                                                
                printf("%i ", el);
                /*NekDouble tmp_inarray[ncoeffs];
                for (int i = 0; i < ncoeffs; ++i)
                {
                    tmp_inarray[i] = transfer_in[coeff_offset[el]+i];
                }
                NekDouble quadMetric [nquad0*nquad1];
                NekDouble laplacian00[nquad0*nquad1];
                NekDouble laplacian01[nquad0*nquad1];
                NekDouble laplacian11[nquad0*nquad1];
                for (int i = 0; i < nquad0*nquad1; ++i)
                {
                     quadMetric[i] =  quadMetricGlo[el*nquad0*nquad1+i];
                    laplacian00[i] = laplacian00Glo[el*nquad0*nquad1+i];
                    laplacian01[i] = laplacian01Glo[el*nquad0*nquad1+i];
                    laplacian11[i] = laplacian11Glo[el*nquad0*nquad1+i];
                }*/
                HelmholtzMatrixOp_MatFree_Kokkos(
                    transfer_in,//tmp_inarray,
                    transfer_out,
                    el, coeff_offset,
                    lambda,
                     quadMetricGlo,
                    laplacian00Glo,
                    laplacian01Glo,
                    laplacian11Glo,
                    nquad0, nquad1, nmodes0, nmodes1, ncoeffs,
                    base0, base1, dbase0, dbase1, D0, D1);
            }
            printf("\n completed all elements\n");             
        }

        KOKKOS_INLINE_FUNCTION
        void GlobalLinSysIterative::HelmholtzMatrixOp_MatFree_Kokkos(
                Kokkos::View<double*> inarray,
                Kokkos::View<double*> outarray,
                const int &el,
                const Kokkos::View<int*>  coeff_offset,
                const Kokkos::View<double[1]> lambda,
				const Kokkos::View<double*>  quadMetricGlo,
                const Kokkos::View<double*> laplacian00Glo,
                const Kokkos::View<double*> laplacian01Glo,
                const Kokkos::View<double*> laplacian11Glo,
                const int &nquad0, const int &nquad1,
                const int &nmodes0, const int &nmodes1, const int &ncoeffs,
                const Kokkos::View<double*> base0,
                const Kokkos::View<double*> base1,
                const Kokkos::View<double*> dbase0,
                const Kokkos::View<double*> dbase1,
                const Kokkos::View<double*> D0,
                const Kokkos::View<double*> D1)
        {
            int nqtot   = nquad0*nquad1;
            //int       wspsize = std::max(std::max(std::max(nqtot,ncoeffs),nquad1*nmodes0), nquad0*nmodes1);
            int max1 = (nqtot >= ncoeffs) ? nqtot : ncoeffs;
            int max2 = (nquad1*nmodes0 >= nquad0*nmodes1) ? nquad1*nmodes0 : nquad0*nmodes1;
            int wspsize = (max1 >= max2) ? max1 : max2;

            // Allocate temporary storage
            NekDouble* wsp0; // use mallloc ?
            NekDouble* wsp1;
            NekDouble* wsp2;

            NekDouble* tmp_outarray;
            
            BwdTrans_SumFacKernel_Kokkos(base0, base1, inarray, wsp0, wsp2,
                nmodes0, nmodes1, nquad0, nquad1);

            for (int i = 0; i < nqtot; ++i)
            {
                wsp1[i] = quadMetricGlo[el*nqtot+i] * wsp0[i];
            }

            IProductWRTBase_SumFacKernel_Kokkos(base0, base1, wsp1, tmp_outarray,
                                         wsp2, nmodes0, nmodes1, nquad0, nquad1);

            //LaplacianMatrixOp_MatFree_Kernel
            // Allocate temporary storage
            NekDouble* wsp0L;
            NekDouble* wsp1L;
            NekDouble* wsp2L;
            
            PhysTensorDeriv_Kokkos(wsp0,wsp1L,wsp2L, nquad0, nquad1, D0, D1);

            for (int i = 0; i < nqtot; ++i)
            {
                wsp0L[i] = laplacian00Glo[el*nqtot+i] * wsp1L[i] 
                         + laplacian01Glo[el*nqtot+i] * wsp2L[i];
            }
            for (int i = 0; i < nqtot; ++i)
            {
                wsp2L[i] = laplacian01Glo[el*nqtot+i] * wsp1L[i]
                		 + laplacian11Glo[el*nqtot+i] * wsp2L[i];
            }

            IProductWRTBase_SumFacKernel_Kokkos(dbase0, base1,wsp0L,wsp1 ,wsp1L,
                         nmodes0, nmodes1, nquad0, nquad1);
            IProductWRTBase_SumFacKernel_Kokkos( base0,dbase1,wsp2L,wsp1L,wsp0L,
                         nmodes0, nmodes1, nquad0, nquad1);

			for (int i = 0; i < ncoeffs; ++i)
            {
                wsp1[i] += wsp1L[i];
            }
            //end LaplacianMatrixOp_MatFree_Kernel
           
            // outarray = lambda * outarray + wsp1
            //          = (lambda * M + L ) * u_hat
            for (int i = 0; i < ncoeffs; ++i)
            {
                outarray[coeff_offset[el] + i] = lambda[0] * tmp_outarray[i] + wsp1[i];
            }    
            
        }

        KOKKOS_INLINE_FUNCTION
        void GlobalLinSysIterative::IProductWRTBase_SumFacKernel_Kokkos(
                const Kokkos::View<double*> base0,
                const Kokkos::View<double*> base1,
                NekDouble* inarray,
                NekDouble* outarray,
                NekDouble* wsp,
                const int &nmodes0, const int &nmodes1,
                const int &nquad0, const int &nquad1)
        {
            Blas::Dgemm('T','N',nquad1,nmodes0,nquad0,1.0,&inarray[0],nquad0,
                        base0.ptr_on_device(),nquad0,0.0,&wsp[0],nquad1);
            int i, mode;
            for (mode=i=0; i < nmodes0; ++i)
            {
                Blas::Dgemv('T',nquad1,nmodes1-i,1.0, base1.ptr_on_device()+mode*nquad1,
                            nquad1,&wsp[0]+i*nquad1,1, 0.0,
                            &outarray[0] + mode,1);
                mode += nmodes1 - i;
            }
            outarray[1] += Blas::Ddot(nquad1,base1.ptr_on_device()+nquad1,1,
                                          &wsp[0]+nquad1,1);
        }

        KOKKOS_INLINE_FUNCTION
        void GlobalLinSysIterative::BwdTrans_SumFacKernel_Kokkos(
                const Kokkos::View<double*> base0,
                const Kokkos::View<double*> base1,
                Kokkos::View<double*> inarray,
                NekDouble* outarray,
                NekDouble* wsp,
                const int &nmodes0, const int &nmodes1,
                const int &nquad0, const int &nquad1)
        {
            int i, mode;
            for (i = mode = 0; i < nmodes0; ++i)
            {
                Blas::Dgemv('N', nquad1,nmodes1-i,1.0,base1.ptr_on_device()+mode*nquad1,
                            nquad1,&inarray[0]+mode,1,0.0,&wsp[0]+i*nquad1,1);
                mode += nmodes1-i;
            }
            Blas::Daxpy(nquad1,inarray[1],base1.ptr_on_device()+nquad1,1,
                            &wsp[0]+nquad1,1);
            Blas::Dgemm('N','T', nquad0,nquad1,nmodes0,1.0, base0.ptr_on_device(),nquad0,
                        &wsp[0], nquad1,0.0, &outarray[0], nquad0);          
        }

        KOKKOS_INLINE_FUNCTION
        void GlobalLinSysIterative::PhysTensorDeriv_Kokkos(
                NekDouble* inarray,
                NekDouble* outarray_d0,
                NekDouble* outarray_d1,
                const int &nquad0, const int &nquad1,
                const Kokkos::View<double*> D0,
                const Kokkos::View<double*> D1)
        {
            Blas::Dgemm('N', 'N', nquad0, nquad1, nquad0, 1.0,
                        (D0.ptr_on_device()), nquad0, &inarray[0], nquad0, 0.0,
                        &outarray_d0[0], nquad0);
            Blas:: Dgemm('N', 'T', nquad0, nquad1, nquad1, 1.0, &inarray[0], nquad0,
                         (D1.ptr_on_device()), nquad1, 0.0, &outarray_d1[0], nquad0);
        }

    }
}