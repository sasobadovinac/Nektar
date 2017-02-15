////////////////////////////////////////////////////////////////////////////////
//
//  File: GetFunctional_eRoca.hxx
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

#ifndef UTILITIES_NEKMESH_NODEOPTI_GETFUNCTIONAL_EROCA
#define UTILITIES_NEKMESH_NODEOPTI_GETFUNCTIONAL_EROCA

#include "Hessian.hxx"

namespace Nektar
{
namespace Utilities
{

/*template<const bool gradient>
struct ProcessVarOpti::GetFunctional<3,gradient,eRoca>
{
    KOKKOS_INLINE_FUNCTION
    NekDouble operator()(const DerivUtilGPU &derivUtilGPU,
             const NodesGPU &nodes, const ElUtilGPU &elUtil, 
             const Grad &grad, int nElmt, int node, int cs,
             const NekDouble ep, const member_type &teamMember)
{   
    const int ptsLowGPU = derivUtilGPU.ptsLow;
    const int ptsHighGPU = derivUtilGPU.ptsHigh;

    grad.minJacNew(node) = DBL_MAX;
    
    grad.integral(node) = 0.0;

    for (int el = 0; el < nElmt; ++el)
    {
        const int elId = nodes.elIdArray(cs,node,el);
        const int localNodeId = nodes.localNodeIdArray(cs,node,el);

        // put node coordinates in shared memory
        ScratchViewType s_nodesX(teamMember.team_scratch(0),3*ptsLowGPU);

        Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsLowGPU ), [&] ( const int n)
        {
            s_nodesX(n) = nodes.X(elId,n);
        });
        Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsLowGPU ), [&] ( const int n)
        {
            s_nodesX(n+ptsLowGPU) = nodes.Y(elId,n);            
        }); 
        Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsLowGPU ), [&] ( const int n)
        {
            s_nodesX(n+2*ptsLowGPU) = nodes.Z(elId,n);
        });  

        teamMember.team_barrier();


        Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsHighGPU ), [&] ( const int k)
        {               
            NekDouble phiM[3][3];
            phiM[0][0] = 0.0;
            phiM[1][0] = 0.0;
            phiM[2][0] = 0.0;
            phiM[0][1] = 0.0;
            phiM[1][1] = 0.0;
            phiM[2][1] = 0.0;
            phiM[0][2] = 0.0;
            phiM[1][2] = 0.0;
            phiM[2][2] = 0.0;

            for (int n = 0; n < ptsLowGPU; ++n)
            {
                phiM[0][0] += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+0*ptsLowGPU);
                phiM[1][0] += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+1*ptsLowGPU);
                phiM[2][0] += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+2*ptsLowGPU);

                phiM[0][1] += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+0*ptsLowGPU);
                phiM[1][1] += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+1*ptsLowGPU);
                phiM[2][1] += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+2*ptsLowGPU);

                phiM[0][2] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+0*ptsLowGPU);
                phiM[1][2] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+1*ptsLowGPU);
                phiM[2][2] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+2*ptsLowGPU);
            }

            NekDouble jacIdeal[3][3];
            jacIdeal[0][0] = phiM[0][0] * elUtil.idealMap(elId,k,0 * 3 + 0)
                           + phiM[0][1] * elUtil.idealMap(elId,k,0 * 3 + 1)
                           + phiM[0][2] * elUtil.idealMap(elId,k,0 * 3 + 2);
            jacIdeal[0][1] = phiM[0][0] * elUtil.idealMap(elId,k,1 * 3 + 0)
                           + phiM[0][1] * elUtil.idealMap(elId,k,1 * 3 + 1)
                           + phiM[0][2] * elUtil.idealMap(elId,k,1 * 3 + 2);
            jacIdeal[0][2] = phiM[0][0] * elUtil.idealMap(elId,k,2 * 3 + 0)
                           + phiM[0][1] * elUtil.idealMap(elId,k,2 * 3 + 1)
                           + phiM[0][2] * elUtil.idealMap(elId,k,2 * 3 + 2);

            jacIdeal[1][0] = phiM[1][0] * elUtil.idealMap(elId,k,0 * 3 + 0)
                           + phiM[1][1] * elUtil.idealMap(elId,k,0 * 3 + 1)
                           + phiM[1][2] * elUtil.idealMap(elId,k,0 * 3 + 2);
            jacIdeal[1][1] = phiM[1][0] * elUtil.idealMap(elId,k,1 * 3 + 0)
                           + phiM[1][1] * elUtil.idealMap(elId,k,1 * 3 + 1)
                           + phiM[1][2] * elUtil.idealMap(elId,k,1 * 3 + 2);
            jacIdeal[1][2] = phiM[1][0] * elUtil.idealMap(elId,k,2 * 3 + 0)
                           + phiM[1][1] * elUtil.idealMap(elId,k,2 * 3 + 1)
                           + phiM[1][2] * elUtil.idealMap(elId,k,2 * 3 + 2);

            jacIdeal[2][0] = phiM[2][0] * elUtil.idealMap(elId,k,0 * 3 + 0)
                           + phiM[2][1] * elUtil.idealMap(elId,k,0 * 3 + 1)
                           + phiM[2][2] * elUtil.idealMap(elId,k,0 * 3 + 2);
            jacIdeal[2][1] = phiM[2][0] * elUtil.idealMap(elId,k,1 * 3 + 0)
                           + phiM[2][1] * elUtil.idealMap(elId,k,1 * 3 + 1)
                           + phiM[2][2] * elUtil.idealMap(elId,k,1 * 3 + 2);
            jacIdeal[2][2] = phiM[2][0] * elUtil.idealMap(elId,k,2 * 3 + 0)
                           + phiM[2][1] * elUtil.idealMap(elId,k,2 * 3 + 1)
                           + phiM[2][2] * elUtil.idealMap(elId,k,2 * 3 + 2);


            NekDouble jacDet = Determinant<3>(jacIdeal);            
            Kokkos::atomic_fetch_min(&grad.minJacNew(node), jacDet);
            NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*ep*ep)); // the regularised Jacobian
            if(sigma < DBL_MIN && !gradient)
            {
                return DBL_MAX;
            }
            ASSERTL0(sigma > DBL_MIN,"dividing by zero");

            NekDouble absIdealMapDet = fabs(elUtil.idealMap(elId,k,9));
            NekDouble quadW = derivUtilGPU.quadW(k);
            NekDouble frob = FrobeniusNorm(jacIdeal);
            //NekDouble W = frob / 3.0 / pow(fabs(sigma), 2.0/3.0);
            NekDouble pow_3 = rcbrt(sigma);
            NekDouble W = frob / 3.0 * pow_3 * pow_3;
            NekDouble inc = quadW * absIdealMapDet * W;
            Kokkos::atomic_add(&grad.integral(node), inc);

            // Derivative of basis function in each direction
            if(gradient)
            {
                NekDouble jacDerivPhi[3];
                NekDouble jacDetDeriv[3];   

                NekDouble derivDet = Determinant<3>(phiM);
                NekDouble jacInvTrans[3][3];                                
                InvTrans<3>(phiM, jacInvTrans, derivDet);
                
                
                // jacDeriv is actually a tensor,
                // but can be stored as a vector, as 18 out of 27 entries are zero
                // and the other 9 entries are three triplets
                // this is due to the delta function in jacDeriv
                NekDouble jacDeriv[3];
                jacDeriv[0] = derivUtilGPU.VdmD_0(k,localNodeId);
                jacDeriv[1] = derivUtilGPU.VdmD_1(k,localNodeId);
                jacDeriv[2] = derivUtilGPU.VdmD_2(k,localNodeId);
                
                // jacDerivPhi is actually a tensor,
                // but can be stored as a vector due to the simple form of jacDeriv                    
                jacDerivPhi[0] = jacDeriv[0] * elUtil.idealMap(elId,k,0 + 3*0)
                               + jacDeriv[1] * elUtil.idealMap(elId,k,1 + 3*0)
                               + jacDeriv[2] * elUtil.idealMap(elId,k,2 + 3*0);
                jacDerivPhi[1] = jacDeriv[0] * elUtil.idealMap(elId,k,0 + 3*1)
                               + jacDeriv[1] * elUtil.idealMap(elId,k,1 + 3*1)
                               + jacDeriv[2] * elUtil.idealMap(elId,k,2 + 3*1); 
                jacDerivPhi[2] = jacDeriv[0] * elUtil.idealMap(elId,k,0 + 3*2)
                               + jacDeriv[1] * elUtil.idealMap(elId,k,1 + 3*2)
                               + jacDeriv[2] * elUtil.idealMap(elId,k,2 + 3*2);  
                
                jacDetDeriv[0] = jacInvTrans[0][0] * jacDeriv[0]
                               + jacInvTrans[0][1] * jacDeriv[1]
                               + jacInvTrans[0][2] * jacDeriv[2];
                jacDetDeriv[0] *= derivDet / absIdealMapDet;
                jacDetDeriv[1] = jacInvTrans[1][0] * jacDeriv[0]
                               + jacInvTrans[1][1] * jacDeriv[1]
                               + jacInvTrans[1][2] * jacDeriv[2];
                jacDetDeriv[1] *= derivDet / absIdealMapDet;
                jacDetDeriv[2] = jacInvTrans[2][0] * jacDeriv[0]
                               + jacInvTrans[2][1] * jacDeriv[1]
                               + jacInvTrans[2][2] * jacDeriv[2];
                jacDetDeriv[2] *= derivDet / absIdealMapDet;
                

                NekDouble frobProd[3];
                NekDouble inc[3];
                // because of the zero entries of the tensor jacDerivPhi,
                // the Frobenius-product becomes a scalar product
                frobProd[0] = ScalarProd<3>(jacIdeal[0],jacDerivPhi);                
                inc[0] = quadW * absIdealMapDet * (
                        2.0*W*(frobProd[0]/frob -
                                jacDetDeriv[0]/3.0/(2.0*sigma-jacDet)));
                Kokkos::atomic_add(&grad.G(node,0), inc[0]);

                frobProd[1] = ScalarProd<3>(jacIdeal[1],jacDerivPhi);                
                inc[1] = quadW * absIdealMapDet * (
                        2.0*W*(frobProd[1]/frob -
                                jacDetDeriv[1]/3.0/(2.0*sigma-jacDet)));
                Kokkos::atomic_add(&grad.G(node,1), inc[1]);

                frobProd[2] = ScalarProd<3>(jacIdeal[2],jacDerivPhi);                
                inc[2] = quadW * absIdealMapDet * (
                        2.0*W*(frobProd[2]/frob -
                                jacDetDeriv[2]/3.0/(2.0*sigma-jacDet)));
                Kokkos::atomic_add(&grad.G(node,2), inc[2]);
                

                NekDouble incHes[6];
                NekDouble frobProdHes = ScalarProd<3>(jacDerivPhi,jacDerivPhi);                                 
                
                incHes[0] = quadW * absIdealMapDet * (
                        inc[0] * inc[0] / W + 2.0*W*(frobProdHes/frob
                        - 2.0 * frobProd[0]*frobProd[0]/frob/frob
                        + jacDetDeriv[0]*jacDetDeriv[0] * jacDet/(2.0*sigma-jacDet)/
                        (2.0*sigma-jacDet)/(2.0*sigma-jacDet)/3.0));
                Kokkos::atomic_add(&grad.G(node,3), incHes[0]);

                incHes[1] = quadW * absIdealMapDet * (
                        inc[0] * inc[1] / W + 2.0*W*(
                        - 2.0 * frobProd[0]*frobProd[1]/frob/frob
                        + jacDetDeriv[0]*jacDetDeriv[1] * jacDet/(2.0*sigma-jacDet)/
                        (2.0*sigma-jacDet)/(2.0*sigma-jacDet)/3.0));
                Kokkos::atomic_add(&grad.G(node,4), incHes[1]);

                incHes[2] = quadW * absIdealMapDet * (
                        inc[0] * inc[2] / W + 2.0*W*(
                        - 2.0 * frobProd[0]*frobProd[2]/frob/frob
                        + jacDetDeriv[0]*jacDetDeriv[2] * jacDet/(2.0*sigma-jacDet)/
                        (2.0*sigma-jacDet)/(2.0*sigma-jacDet)/3.0));
                Kokkos::atomic_add(&grad.G(node,5), incHes[2]);


                incHes[3] = quadW * absIdealMapDet * (
                        inc[1] * inc[1] / W + 2.0*W*(frobProdHes/frob
                        - 2.0 * frobProd[1]*frobProd[1]/frob/frob
                        + jacDetDeriv[1]*jacDetDeriv[1] * jacDet/(2.0*sigma-jacDet)/
                        (2.0*sigma-jacDet)/(2.0*sigma-jacDet)/3.0));
                Kokkos::atomic_add(&grad.G(node,6), incHes[3]);

                incHes[4] = quadW * absIdealMapDet * (
                        inc[1] * inc[2] / W + 2.0*W*(
                        - 2.0 * frobProd[1]*frobProd[2]/frob/frob
                        + jacDetDeriv[1]*jacDetDeriv[2] * jacDet/(2.0*sigma-jacDet)/
                        (2.0*sigma-jacDet)/(2.0*sigma-jacDet)/3.0));
                Kokkos::atomic_add(&grad.G(node,7), incHes[4]);


                incHes[5] = quadW * absIdealMapDet * (
                        inc[2] * inc[2] / W + 2.0*W*(frobProdHes/frob
                        - 2.0 * frobProd[2]*frobProd[2]/frob/frob
                        + jacDetDeriv[2]*jacDetDeriv[2] * jacDet/(2.0*sigma-jacDet)/
                        (2.0*sigma-jacDet)/(2.0*sigma-jacDet)/3.0));
                Kokkos::atomic_add(&grad.G(node,8), incHes[5]);





            }   // gradient

        
        //}   // loop over k
        });
    
    }   // loop over el
    //});     
       
    return grad.integral(node);

}

}; // struct



template<const bool gradient>
struct ProcessVarOpti::GetFunctional<2,gradient,eRoca>
{
    KOKKOS_INLINE_FUNCTION
    NekDouble operator()(const DerivUtilGPU &derivUtilGPU,
             const NodesGPU &nodes, const ElUtilGPU &elUtil, 
             const Grad &grad, int nElmt, int node, int cs,
             const NekDouble ep, const member_type &teamMember)
{   
    const int DIM = 2;
    
    const int ptsLowGPU = derivUtilGPU.ptsLow;
    const int ptsHighGPU = derivUtilGPU.ptsHigh;

    grad.minJacNew(node) = DBL_MAX;
    
    grad.integral(node) = 0.0;

    for (int el = 0; el < nElmt; ++el)
    {
        const int elId = nodes.elIdArray(cs,node,el);
        const int localNodeId = nodes.localNodeIdArray(cs,node,el);

        // put node coordinates in shared memory
        ScratchViewType s_nodesX(teamMember.team_scratch(0),3*ptsLowGPU);

        Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsLowGPU ), [&] ( const int n)
        {
            s_nodesX(n) = nodes.X(elId,n);
        });
        Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsLowGPU ), [&] ( const int n)
        {
            s_nodesX(n+ptsLowGPU) = nodes.Y(elId,n);            
        }); 
        Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsLowGPU ), [&] ( const int n)
        {
            s_nodesX(n+2*ptsLowGPU) = nodes.Z(elId,n);
        });  

        teamMember.team_barrier();


        Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsHighGPU ), [&] ( const int k)
        {               
            NekDouble derivGPU[DIM*DIM];
            derivGPU[0] = 0.0;
            derivGPU[1] = 0.0;
            derivGPU[2] = 0.0;
            derivGPU[3] = 0.0;
            derivGPU[4] = 0.0;
            derivGPU[5] = 0.0;
            derivGPU[6] = 0.0;
            derivGPU[7] = 0.0;
            derivGPU[8] = 0.0;

            for (int n = 0; n < ptsLowGPU; ++n)
            {
                derivGPU[0] += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+0*ptsLowGPU);
                derivGPU[1] += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+1*ptsLowGPU);
                derivGPU[2] += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+2*ptsLowGPU);

                derivGPU[3] += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+0*ptsLowGPU);
                derivGPU[4] += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+1*ptsLowGPU);
                derivGPU[5] += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+2*ptsLowGPU);

                derivGPU[6] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+0*ptsLowGPU);
                derivGPU[7] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+1*ptsLowGPU);
                derivGPU[8] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+2*ptsLowGPU);
            }          
            
            NekDouble quadW = derivUtilGPU.quadW(k);

            NekDouble phiM[DIM][DIM];
            for (int l = 0; l < DIM; ++l)
            {
                for (int n = 0; n < DIM; ++n)
                {
                    phiM[n][l] = derivGPU[l*DIM+n];
                    //phiM[n][l] = derivGPU[l*DIM+n][k];
                }
            }
            // begin CalcIdealJac
            NekDouble jacIdeal[DIM][DIM];
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
            NekDouble jacDet = Determinant<DIM>(jacIdeal);
            // end CalcIdealJac

            Kokkos::atomic_fetch_min(&grad.minJacNew(node), jacDet);

            NekDouble absIdealMapDet = fabs(elUtil.idealMap(elId,k,9));

            NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*ep*ep)); // the regularised Jacobian
            if(sigma < DBL_MIN && !gradient)
            {
                return DBL_MAX;
            }
            ASSERTL0(sigma > DBL_MIN,"dividing by zero");
           
            NekDouble frob = FrobeniusNorm(jacIdeal);
            //NekDouble W = frob / 2.0 / pow(fabs(sigma), 2.0/2.0);
            NekDouble W = frob / 2.0 / sigma;
            NekDouble inc = quadW * absIdealMapDet * W;
            Kokkos::atomic_add(&grad.integral(node), inc);

            // Derivative of basis function in each direction
            if(gradient)
            {
                NekDouble jacDerivPhi[DIM];
                NekDouble jacDetDeriv[DIM];   

                NekDouble derivDet = Determinant<DIM>(phiM);
                NekDouble jacInvTrans[DIM][DIM];                                
                InvTrans<DIM>(phiM, jacInvTrans, derivDet);

                NekDouble basisDeriv [DIM];
                basisDeriv[0] = derivUtilGPU.VdmD_0(k,localNodeId);
                basisDeriv[1] = derivUtilGPU.VdmD_1(k,localNodeId);
                basisDeriv[2] = derivUtilGPU.VdmD_2(k,localNodeId);                        

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

                NekDouble frobProd[DIM];
                NekDouble inc[DIM];
                for (int m = 0; m < DIM; ++m)
                {
                    // because of the zero entries of the tensor jacDerivPhi,
                    // the Frobenius-product becomes a scalar product
                    frobProd[m] = ScalarProd<DIM>(jacIdeal[m],jacDerivPhi);
                
                    inc[m] = quadW * absIdealMapDet * (
                            2.0*W*(frobProd[m]/frob -
                                    jacDetDeriv[m]/2.0/(2.0*sigma-jacDet)));
                    Kokkos::atomic_add(&grad.G(node,m), inc[m]);
                }
                
                int ct = 0;
                for (int m = 0; m < DIM; ++m)
                {
                    for(int l = m; l < DIM; ++l,  ct++)
                    {
                        NekDouble frobProdHes = 0.0;
                        // because of the zero entries of the tensor jacDerivPhi,
                        // the matrix frobProdHes has only diagonal entries
                        if (m == l)
                        {
                            // because of the zero entries of the tensor jacDerivPhi,
                            // the Frobenius-product becomes a scalar product
                            frobProdHes = ScalarProd<DIM>(jacDerivPhi,jacDerivPhi);                                
                        }
                    
                        NekDouble incHes = quadW * absIdealMapDet * (
                                inc[m] * inc[l] / W + 2.0*W*(frobProdHes/frob
                                - 2.0 * frobProd[m]*frobProd[l]/frob/frob
                                + jacDetDeriv[m]*jacDetDeriv[l] * jacDet/(2.0*sigma-jacDet)/
                                (2.0*sigma-jacDet)/(2.0*sigma-jacDet)/2.0));
                        Kokkos::atomic_add(&grad.G(node,ct+DIM), incHes);
                    }
                }   

            }   // gradient

        
        //}   // loop over k
        });
    
    }   // loop over el
    //});     
       
    return grad.integral(node);

}

}; // struct

*/


template<const int DIM, const bool gradient>
struct ProcessVarOpti::GetFunctional<DIM,gradient,eRoca>
{
    KOKKOS_INLINE_FUNCTION
    NekDouble operator()(const DerivUtilGPU &derivUtilGPU,
             const NodesGPU &nodes, const ElUtilGPU &elUtil, 
             const Grad &grad, int nElmt, int node, int cs,
             const NekDouble ep, const member_type &teamMember)
{   
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

            NekDouble derivGPU[DIM*DIM];
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

            if (DIM == 3)
            {
                for (int n = 0; n < ptsLowGPU; ++n)
                {
                    derivGPU[0] += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+0*ptsLowGPU);
                    derivGPU[1] += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+1*ptsLowGPU);
                    derivGPU[2] += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+2*ptsLowGPU);
                    
                    derivGPU[3] += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+0*ptsLowGPU);
                    derivGPU[4] += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+1*ptsLowGPU);                
                    derivGPU[5] += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+2*ptsLowGPU);
                    
                    derivGPU[6] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+0*ptsLowGPU);
                    derivGPU[7] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+1*ptsLowGPU);
                    derivGPU[8] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+2*ptsLowGPU);
                    
                }
            }
            else if (DIM == 2)
            {
                for (int n = 0; n < ptsLowGPU; ++n)
                {
                    derivGPU[0] += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+0*ptsLowGPU);
                    derivGPU[1] += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+1*ptsLowGPU);
                    
                    derivGPU[2] += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+0*ptsLowGPU);
                    derivGPU[3] += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+1*ptsLowGPU);                    
                }
            }    
            
            NekDouble quadW = derivUtilGPU.quadW(k);

            NekDouble phiM[DIM][DIM];
            for (int l = 0; l < DIM; ++l)
            {
                for (int n = 0; n < DIM; ++n)
                {
                    phiM[n][l] = derivGPU[l*DIM+n];
                    //phiM[n][l] = derivGPU[l*DIM+n][k];
                }
            }
            // begin CalcIdealJac
            NekDouble jacIdeal[DIM][DIM];
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
            NekDouble jacDet = Determinant<DIM>(jacIdeal);
            // end CalcIdealJac

            Kokkos::atomic_fetch_min(&grad.minJacNew(node), jacDet);

            NekDouble absIdealMapDet = fabs(elUtil.idealMap(elId,k,9));

            NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*ep*ep)); // the regularised Jacobian
            if(sigma < DBL_MIN && !gradient)
            {
                return DBL_MAX;
            }
            ASSERTL0(sigma > DBL_MIN,"dividing by zero");
           
            NekDouble frob = FrobeniusNorm(jacIdeal);
            NekDouble W;
            if (DIM == 3)
            {
                NekDouble pow_3 = rcbrt(sigma);
                W = frob / 3.0 * pow_3 * pow_3;
            }
            else if (DIM == 2)
            {
                W = frob / 2.0 / sigma;
            }
            
            NekDouble inc = quadW * absIdealMapDet * W;
            Kokkos::atomic_add(&grad.integral(node), inc);

            // Derivative of basis function in each direction
            if(gradient)
            {
                NekDouble jacDerivPhi[DIM];
                NekDouble jacDetDeriv[DIM];   

                NekDouble derivDet = Determinant<DIM>(phiM);
                NekDouble jacInvTrans[DIM][DIM];                                
                InvTrans<DIM>(phiM, jacInvTrans, derivDet);

                NekDouble basisDeriv [DIM];
                basisDeriv[0] = derivUtilGPU.VdmD_0(k,localNodeId);
                basisDeriv[1] = derivUtilGPU.VdmD_1(k,localNodeId);
                basisDeriv[2] = derivUtilGPU.VdmD_2(k,localNodeId);                        

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

                NekDouble frobProd[DIM];
                NekDouble inc[DIM];
                for (int m = 0; m < DIM; ++m)
                {
                    // because of the zero entries of the tensor jacDerivPhi,
                    // the Frobenius-product becomes a scalar product
                    frobProd[m] = ScalarProd<DIM>(jacIdeal[m],jacDerivPhi);
                
                    inc[m] = quadW * absIdealMapDet * (
                            2.0*W*(frobProd[m]/frob -
                                    jacDetDeriv[m]/DIM/(2.0*sigma-jacDet)));
                    Kokkos::atomic_add(&grad.G(node,m), inc[m]);
                }
                
                int ct = 0;
                for (int m = 0; m < DIM; ++m)
                {
                    for(int l = m; l < DIM; ++l,  ct++)
                    {
                        NekDouble frobProdHes = 0.0;
                        // because of the zero entries of the tensor jacDerivPhi,
                        // the matrix frobProdHes has only diagonal entries
                        if (m == l)
                        {
                            // because of the zero entries of the tensor jacDerivPhi,
                            // the Frobenius-product becomes a scalar product
                            frobProdHes = ScalarProd<DIM>(jacDerivPhi,jacDerivPhi);                                
                        }
                    
                        NekDouble incHes = quadW * absIdealMapDet * (
                                inc[m] * inc[l] / W + 2.0*W*(frobProdHes/frob
                                - 2.0 * frobProd[m]*frobProd[l]/frob/frob
                                + jacDetDeriv[m]*jacDetDeriv[l] * jacDet/(2.0*sigma-jacDet)/
                                (2.0*sigma-jacDet)/(2.0*sigma-jacDet)/DIM));
                        Kokkos::atomic_add(&grad.G(node,ct+DIM), incHes);
                    }
                }   

            }   // gradient

        
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
