////////////////////////////////////////////////////////////////////////////////
//
//  File: GetFunctional_eHypEl.hxx
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_NODEOPTI_GETFUNCTIONAL_EHYPEL
#define UTILITIES_NEKMESH_NODEOPTI_GETFUNCTIONAL_EHYPEL

#include "Hessian.hxx"

namespace Nektar
{
namespace Utilities
{


template<const int DIM, const bool gradient>
struct ProcessVarOpti::GetFunctional<DIM,gradient,eHypEl>
{
    KOKKOS_INLINE_FUNCTION
    NekDouble operator()(const DerivUtilGPU &derivUtilGPU,
             const NodesGPU &nodes, const ElUtilGPU &elUtil, 
             const Grad &grad, int nElmt, int node, int cs,//const int elId, const int localNodeId,
             const double ep, const member_type &teamMember)
{ 
    //printf("%s\n", "in GetFunctional");
    const double nu = 0.45;
    const double mu = 1.0 / 2.0 / (1.0+nu);
    const double K  = 1.0 / 3.0 / (1.0 - 2.0 * nu); 

    const int ptsLowGPU = derivUtilGPU.ptsLow;
    const int ptsHighGPU = derivUtilGPU.ptsHigh;
    
    grad.minJacNew(node) = DBL_MAX;

    grad.integral(node) = 0.0;
    
    
    for (int el = 0; el < nElmt; ++el)
    {
        const int elId = nodes.elIdArray(cs,node,el);
        const int localNodeId = nodes.localNodeIdArray(cs,node,el);

        // put node coordinates in shared memory
        ScratchViewType s_nodesX(teamMember.team_scratch(0),DIM*ptsLowGPU);

        Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsLowGPU ), [&] ( const int n)
        {
            s_nodesX(n) = nodes.X(elId,n);
        });
        Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsLowGPU ), [&] ( const int n)
        {
            s_nodesX(n+ptsLowGPU) = nodes.Y(elId,n);            
        });
        if (DIM == 3)
        {
            Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsLowGPU ), [&] ( const int n)
            {
                s_nodesX(n+2*ptsLowGPU) = nodes.Z(elId,n);
            }); 
        }

        teamMember.team_barrier();


        Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsHighGPU ), [&] ( const int k)
        {               
            // Timings
            /*#ifdef __CUDA_ARCH__
            long long int start = clock64();
            #endif*/        

            double derivGPU[DIM*DIM];
            derivGPU[0] = 0.0;
            derivGPU[1] = 0.0;
            derivGPU[2] = 0.0;
            derivGPU[3] = 0.0;
            if (DIM == 3)
            {
                derivGPU[4] = 0.0;
                derivGPU[5] = 0.0;
                derivGPU[6] = 0.0;
                derivGPU[7] = 0.0;
                derivGPU[8] = 0.0;
            }

            for (int n = 0; n < ptsLowGPU; ++n)
            {
                derivGPU[0] += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+0*ptsLowGPU);
                derivGPU[1] += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+1*ptsLowGPU);
                if (DIM == 3)
                {
                    derivGPU[2] += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+2*ptsLowGPU);
                }
                derivGPU[3] += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+0*ptsLowGPU);
                derivGPU[4] += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+1*ptsLowGPU);

                if (DIM == 3)
                {
                    derivGPU[5] += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+2*ptsLowGPU);
                    
                    derivGPU[6] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+0*ptsLowGPU);
                    derivGPU[7] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+1*ptsLowGPU);
                    derivGPU[8] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+2*ptsLowGPU);
                }
            } 
            /*#ifdef __CUDA_ARCH__
            long long int stop = clock64();
            long long int cycles = stop - start;
            if (k ==0){printf("GPUTimer Node %i, QP %i: %lli\n", node, k, cycles); }

            #endif*/
            
            /*
            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( teamMember , ptsLowGPU ), [&] ( const int n, double &update)
            {
                update += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+0*ptsLowGPU);
                
            }, derivGPU[0]);          */
            

            double quadW = derivUtilGPU.quadW(k);

            double phiM[DIM][DIM];
            for (int l = 0; l < DIM; ++l)
            {
                for (int n = 0; n < DIM; ++n)
                {
                    phiM[n][l] = derivGPU[l*DIM+n];
                    //phiM[n][l] = derivGPU[l*DIM+n][k];
                }
            }
            // begin CalcIdealJac
            double jacIdeal[DIM][DIM];
            for (int m = 0; m < DIM; ++m)
            {
                for (int n = 0; n < DIM; ++n)
                {
                    jacIdeal[n][m] = 0.0;
                    for (int l = 0; l < DIM; ++l)
                    {
                        jacIdeal[n][m] += phiM[n][l] *
                            elUtil.idealMap(elId,k,m * 3 + l);               
                    }
                }
            }
            double jacDet = Determinant<DIM>(jacIdeal);
            // end CalcIdealJac

            Kokkos::atomic_fetch_min(&grad.minJacNew(node), jacDet);

            double absIdealMapDet = fabs(elUtil.idealMap(elId,k,9));

            double sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*ep*ep)); // the regularised Jacobian

            if(sigma < DBL_MIN && !gradient)
            {
                return DBL_MAX;
            }
            ASSERTL0(sigma > DBL_MIN,"dividing by zero");                      
                                
            double lsigma = log(sigma);
            double I1 = FrobeniusNorm<DIM>(jacIdeal);                           
            double inc = quadW * absIdealMapDet *
                        (0.5 * mu * (I1 - 3.0 - 2.0*lsigma) +
                         0.5 * K * lsigma * lsigma);
            Kokkos::atomic_add(&grad.integral(node), inc);
         
            if (gradient)
            {
                // Derivative of basis function in each direction
                NekDouble jacDerivPhi[DIM];
                NekDouble jacDetDeriv[DIM]; 

                NekDouble derivDet = Determinant<DIM>(phiM);
                NekDouble jacInvTrans[DIM][DIM];                                
                InvTrans<DIM>(phiM, jacInvTrans);
                
                NekDouble basisDeriv [DIM];
                basisDeriv[0] = derivUtilGPU.VdmD_0(k,localNodeId);
                basisDeriv[1] = derivUtilGPU.VdmD_1(k,localNodeId);
                if (DIM == 3)
                {
                    basisDeriv[2] = derivUtilGPU.VdmD_2(k,localNodeId);                        
                }
                // jacDeriv is actually a tensor,
                // but can be stored as a vector, as 18 out of 27 entries are zero
                // and the other 9 entries are three triplets
                // this is due to the delta function in jacDeriv
                NekDouble jacDeriv[DIM];
                for (int l = 0; l < DIM; ++l)
                {
                    jacDeriv[l] = basisDeriv[l];
                }

                // jacDerivPhi is actually a tensor,
                // but can be stored as a vector due to the simple form of jacDeriv 
                for (int n = 0; n < DIM; ++n)
                {
                    jacDerivPhi[n] = 0.0;                                   
                    for (int l = 0; l < DIM; ++l)
                    {
                        jacDerivPhi[n] += jacDeriv[l] * elUtil.idealMap(elId,k,l + 3*n);                        
                    }                    
                }

                for (int m = 0; m < DIM; ++m)
                {
                    jacDetDeriv[m] = 0.0;
                    for (int n = 0; n < DIM; ++n)
                    {
                        jacDetDeriv[m] += jacInvTrans[m][n] * basisDeriv[n];
                    }
                    jacDetDeriv[m] *= derivDet / absIdealMapDet;
                }
                // end of common part to all four versions

                for (int m = 0; m < DIM; ++m)
                {
                    // because of the zero entries of the tensor jacDerivPhi,
                    // the Frobenius-product becomes a scalar product
                    double frobProd = ScalarProd<DIM>(jacIdeal[m],jacDerivPhi);
                                        
                    double inc = quadW * absIdealMapDet * (
                        mu * frobProd + (jacDetDeriv[m] / (2.0*sigma - jacDet)
                                            * (K * lsigma - mu)));
                    Kokkos::atomic_add(&grad.G(node,m), inc);
                }
                
                int ct = 0;
                for (int m = 0; m < DIM; ++m)
                {
                    for(int l = m; l < DIM; ++l, ct++)
                    {
                        double frobProdHes = 0.0;
                        // because of the zero entries of the tensor jacDerivPhi,
                        // the matrix frobProdHes has only diagonal entries
                        if (m == l)
                        {
                            // because of the zero entries of the tensor jacDerivPhi,
                            // the Frobenius-product becomes a scalar product
                            frobProdHes = ScalarProd<DIM>(jacDerivPhi,jacDerivPhi);                                
                        }                        
        
                        double inc = quadW * absIdealMapDet 
                            * (mu * frobProdHes
                            + jacDetDeriv[m]*jacDetDeriv[l]/(2.0*sigma-jacDet)/(2.0*sigma-jacDet)
                                * (K- jacDet*(K*lsigma-mu)/(2.0*sigma-jacDet)));
                        Kokkos::atomic_add(&grad.G(node,ct+DIM), inc);
                    }
                }                        
            } //gradient
        
        //}   // loop over k
        });
    
    }   // loop over el
    //});     
       
    return grad.integral(node);

}

}; // struct


} // Utilities

} // Nektar

#endif
