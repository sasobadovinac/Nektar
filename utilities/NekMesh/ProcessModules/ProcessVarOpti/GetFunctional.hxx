////////////////////////////////////////////////////////////////////////////////
//
//  File: GetFunctional.hxx
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
//  Description: Calculate jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_NODEOPTI_GETFUNCTIONAL
#define UTILITIES_NEKMESH_NODEOPTI_GETFUNCTIONAL

#include "Hessian.hxx"

namespace Nektar
{
namespace Utilities
{



/*template<const int DIM, const bool gradient, const optimiser opti>
struct ProcessVarOpti::GetFunctional
{
    KOKKOS_INLINE_FUNCTION
    NekDouble operator()(const DerivUtilGPU &derivUtilGPU,
             const NodesGPU &nodes, const ElUtilGPU &elUtil, 
             const Grad &grad, int nElmt, int node, int cs,//const int elId, const int localNodeId,
             const double ep, const member_type &teamMember)
    {
    }
};*/

/*template<const int DIM, const bool gradient>
struct ProcessVarOpti::GetFunctional<DIM,gradient,eHypEl>
{
    KOKKOS_INLINE_FUNCTION
    NekDouble operator()(const DerivUtilGPU &derivUtilGPU,
             const NodesGPU &nodes, const ElUtilGPU &elUtil, 
             const Grad &grad, int nElmt, int node, int cs,//const int elId, const int localNodeId,
             const double ep, const member_type &teamMember)
*/

template<const int DIM, const bool gradient>
KOKKOS_INLINE_FUNCTION
NekDouble ProcessVarOpti::GetFunctional(const DerivUtilGPU &derivUtilGPU,
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
            double derivGPU[DIM*DIM];
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
            //}
            //for (int n = 0; n < ptsLowGPU; ++n)
            //{
                derivGPU[6] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+0*ptsLowGPU);
                derivGPU[7] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+1*ptsLowGPU);
                derivGPU[8] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+2*ptsLowGPU);
            }          
            
                        /*
            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( teamMember , ptsLowGPU ), [&] ( const int n, double &update)
            {
                update += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+0*ptsLowGPU);
                
            }, derivGPU[0]);

            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( teamMember , ptsLowGPU ), [&] ( const int n, double &update)
            {
                update += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+1*ptsLowGPU);
                
            }, derivGPU[1]);

            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( teamMember , ptsLowGPU ), [&] ( const int n, double &update)
            {
                update += derivUtilGPU.VdmD_0(k,n) * s_nodesX(n+2*ptsLowGPU);
                
            }, derivGPU[2]);

            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( teamMember , ptsLowGPU ), [&] ( const int n, double &update)
            {
                update += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+0*ptsLowGPU);
                
            }, derivGPU[3]);

            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( teamMember , ptsLowGPU ), [&] ( const int n, double &update)
            {
                update += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+1*ptsLowGPU);
                
            }, derivGPU[4]);

            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( teamMember , ptsLowGPU ), [&] ( const int n, double &update)
            {
                update += derivUtilGPU.VdmD_1(k,n) * s_nodesX(n+2*ptsLowGPU);
                
            }, derivGPU[5]);

            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( teamMember , ptsLowGPU ), [&] ( const int n, double &update)
            {
                update += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+0*ptsLowGPU);
                
            }, derivGPU[6]);

            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( teamMember , ptsLowGPU ), [&] ( const int n, double &update)
            {
                update += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+1*ptsLowGPU);
                
            }, derivGPU[7]);

            Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( teamMember , ptsLowGPU ), [&] ( const int n, double &update)
            {
                update += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+2*ptsLowGPU);
                
            }, derivGPU[8]);
            */
            

            double quadW = derivUtilGPU.quadW(k);

            double phiM[DIM][DIM];
            // begin CalcIdealJac
            double jacIdeal[DIM][DIM];
            for (int m = 0; m < DIM; ++m)
            {
                for (int n = 0; n < DIM; ++n)
                {
                    jacIdeal[n][m] = 0.0;
                    for (int l = 0; l < DIM; ++l)
                    {
                        phiM[n][l] = derivGPU[l*DIM+n];
                        //phiM[n][l] = derivGPU[l*DIM+n][k];

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
            }
        
        //}   // loop over k
        });
    
    }   // loop over el
    //});     
       
    return grad.integral(node);

}

//}; // struct


/*
template<const int DIM, const bool gradient>
struct ProcessVarOpti::GetFunctional<DIM,gradient,eRoca>
{
    KOKKOS_INLINE_FUNCTION
    NekDouble operator()(const DerivUtilGPU &derivUtilGPU,
             const NodesGPU &nodes, const ElUtilGPU &elUtil, 
             const Grad &grad, int nElmt, int node, int cs,//const int elId, const int localNodeId,
             const double ep, const member_type &teamMember, double &minJacNew)
{     
    
    optimiser opti = eHypEl;

    //printf("%s\n", "in GetFunctional");
    const double nu = 0.45;
    const double mu = 1.0 / 2.0 / (1.0+nu);
    const double K  = 1.0 / 3.0 / (1.0 - 2.0 * nu); 

    const int ptsLowGPU = derivUtilGPU.ptsLow;
    const int ptsHighGPU = derivUtilGPU.ptsHigh;
    // Storage for derivatives, ordered by:
    //   - standard coordinate direction
    //   - cartesian coordinate direction, combined
    //   - quadrature points    
    //double derivGPU[DIM*DIM][ptsHighGPU];
    
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
            double derivGPU[DIM*DIM];
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
            //}
            //for (int n = 0; n < ptsLowGPU; ++n)
            //{
                derivGPU[6] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+0*ptsLowGPU);
                derivGPU[7] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+1*ptsLowGPU);
                derivGPU[8] += derivUtilGPU.VdmD_2(k,n) * s_nodesX(n+2*ptsLowGPU);
            }          
            


        //});
        //Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsHighGPU ), [&] ( const int k)
        //{

            double quadW = derivUtilGPU.quadW(k);

            double phiM[DIM][DIM];
            // begin CalcIdealJac
            double jacIdeal[DIM][DIM];
            for (int m = 0; m < DIM; ++m)
            {
                for (int n = 0; n < DIM; ++n)
                {
                    jacIdeal[n][m] = 0.0;
                    for (int l = 0; l < DIM; ++l)
                    {
                        phiM[n][l] = derivGPU[l*DIM+n];
                        //phiM[n][l] = derivGPU[l*DIM+n][k];


                        //jacIdeal[n][m] += derivGPU(l*DIM+n,k) *
                        //    elUtil.idealMap(elId,k,m * 3 + l);
                        jacIdeal[n][m] += phiM[n][l] *
                            elUtil.idealMap(elId,k,m * 3 + l);               
                    }
                }
            }
            double jacDet = Determinant<DIM>(jacIdeal);
            // end CalcIdealJac

            double absIdealMapDet = fabs(elUtil.idealMap(elId,k,9));

            double sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*ep*ep)); // the regularised Jacobian
            if(sigma < DBL_MIN && !gradient)
            {
                return DBL_MAX;
            }
            ASSERTL0(sigma > DBL_MIN,"dividing by zero");
            
            // Derivative of basis function in each direction
            NekDouble jacDerivPhi[DIM];
            NekDouble jacDetDeriv[DIM];                  
            if(gradient)
            {
                NekDouble jacInvTrans[DIM][DIM];                                
                InvTrans<DIM>(phiM, jacInvTrans);
                NekDouble derivDet = Determinant<DIM>(phiM);

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
            }


            switch(opti)
            {               
                case eLinEl:
                {        
                    double lsigma = log(sigma);
                    NekDouble Emat[DIM][DIM];
                    EMatrix<DIM>(jacIdeal,Emat);
                    NekDouble trEtE = FrobeniusNorm<DIM>(Emat);
                    
                    double inc = quadW * absIdealMapDet *
                                (K * 0.5 * lsigma * lsigma + mu * trEtE);
                    Kokkos::atomic_add(&grad.integral(node), inc);

                    // Derivative of basis function in each direction
                    if(gradient)
                    {                       
                        NekDouble dEdxi[DIM][DIM][DIM];
                        // use the delta function in jacDeriv and do some tensor calculus
                        // to come up with this simplified expression for:
                        //CalcLinElGrad<DIM>(jacIdeal,jacDerivPhi,dEdxi);
                        for(int d = 0; d < DIM; d++)
                        {
                            for (int m = 0; m < DIM; ++m)
                            {
                                for (int n = 0; n < DIM; ++n)
                                {
                                    dEdxi[d][m][n] = 0.5*(jacDerivPhi[m] * jacIdeal[d][n]
                                                        + jacIdeal[d][m] * jacDerivPhi[n]);
                                }
                            }
                        }
                    
                        for (int m = 0; m < DIM; ++m)
                        {   
                            double frobProd = FrobProd<DIM>(dEdxi[m],Emat);
                        
                            double inc = quadW * absIdealMapDet * (
                                    2.0 * mu * frobProd + K * lsigma * jacDetDeriv[m] / (2.0*sigma - jacDet));
                            Kokkos::atomic_add(&grad.G(node,m), inc);
                        }                        
                        
                        //NekDouble frobProdHes[DIM][DIM]; //holder for the hessian frobprods
                        int ct = 0;
                        for (int m = 0; m < DIM; ++m)
                        {
                            for(int l = m; l < DIM; ++l,  ct++)
                            {
                                double frobProdHes = FrobProd<DIM>(dEdxi[m],dEdxi[l]);
                                
                                NekDouble d2Edxi[DIM][DIM];
                                // use the delta function in jacDeriv and do some tensor calculus
                                // to come up with this simplified expression for:
                                //CalcLinElSecGrad<DIM>(jacDerivPhi[m],jacDerivPhi[l],d2Edxi);
                                if (m == l)
                                {
                                    for (int p = 0; p < DIM; ++p)
                                    {
                                        for (int q = 0; q < DIM; ++q)
                                        {                                                
                                            d2Edxi[p][q] = jacDerivPhi[p] * jacDerivPhi[q];                                                
                                        }
                                    }    
                                    frobProdHes += FrobProd<DIM>(d2Edxi,Emat);
                                }
                       
                                double inc = quadW * absIdealMapDet * (2.0 * mu * frobProdHes +
                                    jacDetDeriv[m]*jacDetDeriv[l]*K/(2.0*sigma-jacDet)/(2.0*sigma-jacDet)*(
                                        1.0 - jacDet*lsigma/(2.0*sigma-jacDet)));                                        
                                Kokkos::atomic_add(&grad.G(node,ct+DIM), inc);
                            }
                        }                           
                        
                    }                                 
                    break;
                }

                case eHypEl:
                {                    
                    double lsigma = log(sigma);
                    double I1 = FrobeniusNorm<DIM>(jacIdeal);                           
                    double inc = quadW * absIdealMapDet *
                                (0.5 * mu * (I1 - 3.0 - 2.0*lsigma) +
                                 0.5 * K * lsigma * lsigma);
                    Kokkos::atomic_add(&grad.integral(node), inc);
                 
                    if (gradient)
                    {
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
                    }
                    break;                      
                }

                case eRoca:
                {
                    NekDouble frob = FrobeniusNorm(jacIdeal);
                    NekDouble W = frob / DIM / pow(fabs(sigma), 2.0/DIM);
                    double inc = quadW * absIdealMapDet * W;
                    Kokkos::atomic_add(&grad.integral(node), inc);

                    // Derivative of basis function in each direction
                    if(gradient)
                    {
                        NekDouble frobProd[DIM];
                        double inc[DIM];
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
                                double frobProdHes = 0.0;
                                // because of the zero entries of the tensor jacDerivPhi,
                                // the matrix frobProdHes has only diagonal entries
                                if (m == l)
                                {
                                    // because of the zero entries of the tensor jacDerivPhi,
                                    // the Frobenius-product becomes a scalar product
                                    frobProdHes = ScalarProd<DIM>(jacDerivPhi,jacDerivPhi);                                
                                }
                            
                                double incHes = quadW * absIdealMapDet * (
                                    //grad.G(node,m)*grad.G(node,l) / W + 2.0*W*(frobProdHes/frob
                                        inc[m] * inc[l] / W + 2.0*W*(frobProdHes/frob
                                        - 2.0 * frobProd[m]*frobProd[l]/frob/frob
                                        + jacDetDeriv[m]*jacDetDeriv[l] * jacDet/(2.0*sigma-jacDet)/
                                        (2.0*sigma-jacDet)/(2.0*sigma-jacDet)/DIM));
                                Kokkos::atomic_add(&grad.G(node,ct+DIM), incHes);
                            }
                        }                        
                    }
                    break;
                }

                case eWins:
                {
                    NekDouble frob = FrobeniusNorm(jacIdeal);
                    NekDouble W = frob / sigma;
                    double inc = quadW * absIdealMapDet * W;
                    Kokkos::atomic_add(&grad.integral(node), inc);

                    // Derivative of basis function in each direction
                    if(gradient)
                    {
                        NekDouble frobProd[DIM];  
                        double inc[DIM];
                        for (int m = 0; m < DIM; ++m)
                        {
                            // because of the zero entries of the tensor jacDerivPhi,
                            // the Frobenius-product becomes a scalar product
                            frobProd[m] = ScalarProd<DIM>(jacIdeal[m],jacDerivPhi);                        
                        
                            inc[m] = quadW * absIdealMapDet * (
                                    W*(2.0*frobProd[m]/frob -
                                            jacDetDeriv[m]/(2.0*sigma-jacDet)));
                            Kokkos::atomic_add(&grad.G(node,m), inc[m]);
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
                                double incHes = quadW * absIdealMapDet * (
                                    //grad.G(node,m)*grad.G(node,l) / W + 2.0*W*(frobProdHes/frob
                                    inc[m] * inc[l] / W + 2.0*W*(frobProdHes/frob
                                        - 2.0 * frobProd[m]*frobProd[l]/frob/frob
                                        + 0.5*jacDetDeriv[m]*jacDetDeriv[l] * jacDet/(2.0*sigma-jacDet)
                                        /(2.0*sigma-jacDet)/(2.0*sigma-jacDet)));
                                Kokkos::atomic_add(&grad.G(node,ct+DIM), incHes);
                            }
                        }                        
                    }
                    break;
                }

            }   // switch

        
        //}   // loop over k
        });
    
    }   // loop over el
    //});     
       
    return grad.integral(node);

}

}; // struct
*/

} // Utilities

} // Nektar

#endif
