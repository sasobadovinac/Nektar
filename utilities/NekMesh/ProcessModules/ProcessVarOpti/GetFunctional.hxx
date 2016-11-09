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


template<int DIM>
KOKKOS_INLINE_FUNCTION
NekDouble ProcessVarOpti::GetFunctional(const DerivUtilGPU &derivUtilGPU,
         const NodesGPU &nodes, const ElUtilGPU &elUtil, 
         const Grad &grad, int nElmt, int node, int cs,//const int elId, const int localNodeId,
         const double ep, const member_type &teamMember, const int opti,
         bool gradient, bool hessian)

{ 
    //printf("%s\n", "in GetFunctional");
    const double nu = 0.4;
    const double mu = 1.0 / 2.0 / (1.0+nu);
    const double K  = 1.0 / 3.0 / (1.0 - 2.0 * nu); 

    const int ptsLowGPU = derivUtilGPU.ptsLow;
    const int ptsHighGPU = derivUtilGPU.ptsHigh;
    // Storage for derivatives, ordered by:
    //   - standard coordinate direction
    //   - cartesian coordinate direction, combined
    //   - quadrature points    
    //Kokkos::View<double**> derivGPU("derivGPU", DIM*DIM, ptsHighGPU);
       
    
    grad.integral(node) = 0.0;

    //Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , nElmt ), [&] ( const int el)
    //{
    for (int el = 0; el < nElmt; ++el)
    {
        const int elId = nodes.elIdArray(cs,node,el);
        const int localNodeId = nodes.localNodeIdArray(cs,node,el);   
        
        Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsHighGPU ), [&] ( const int k)
        {
        //Kokkos::parallel_for( Kokkos::ThreadVectorRange( teamMember , ptsHighGPU ), [&] ( const int k)        
        //{        
        //for (int k = 0; k < ptsHighGPU; ++k)
        //{
            /*
            for(int i = 0; i < DIM*DIM; i++)
            {
                derivGPU(i,k) = 0.0;
            }
            
            for (int n = 0; n < ptsLowGPU; ++n)
            {
                derivGPU(0,k) += derivUtilGPU.VdmD_0(k,n) * nodes.X(elId,n);
                derivGPU(1,k) += derivUtilGPU.VdmD_0(k,n) * nodes.Y(elId,n);
                derivGPU(2,k) += derivUtilGPU.VdmD_0(k,n) * nodes.Z(elId,n);

                derivGPU(3,k) += derivUtilGPU.VdmD_1(k,n) * nodes.X(elId,n);
                derivGPU(4,k) += derivUtilGPU.VdmD_1(k,n) * nodes.Y(elId,n);
                derivGPU(5,k) += derivUtilGPU.VdmD_1(k,n) * nodes.Z(elId,n);

                derivGPU(6,k) += derivUtilGPU.VdmD_2(k,n) * nodes.X(elId,n);
                derivGPU(7,k) += derivUtilGPU.VdmD_2(k,n) * nodes.Y(elId,n);
                derivGPU(8,k) += derivUtilGPU.VdmD_2(k,n) * nodes.Z(elId,n);
            }*/
            double derivGPU[DIM*DIM];
            for(int i = 0; i < DIM*DIM; i++)
            {
                derivGPU[i] = 0.0;
            }            
            
            for (int n = 0; n < ptsLowGPU; ++n)
            {
                derivGPU[0] += derivUtilGPU.VdmD_0(k,n) * nodes.X(elId,n);
                derivGPU[1] += derivUtilGPU.VdmD_0(k,n) * nodes.Y(elId,n);
                derivGPU[2] += derivUtilGPU.VdmD_0(k,n) * nodes.Z(elId,n);

                derivGPU[3] += derivUtilGPU.VdmD_1(k,n) * nodes.X(elId,n);
                derivGPU[4] += derivUtilGPU.VdmD_1(k,n) * nodes.Y(elId,n);
                derivGPU[5] += derivUtilGPU.VdmD_1(k,n) * nodes.Z(elId,n);

                derivGPU[6] += derivUtilGPU.VdmD_2(k,n) * nodes.X(elId,n);
                derivGPU[7] += derivUtilGPU.VdmD_2(k,n) * nodes.Y(elId,n);
                derivGPU[8] += derivUtilGPU.VdmD_2(k,n) * nodes.Z(elId,n);
            }
        //});
        
        //Kokkos::parallel_for (range_policy(0,ptsHighGPU), KOKKOS_LAMBDA (const int k)
        //{
            
            double absIdealMapDet = fabs(elUtil.idealMap(elId,k,9));
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
                        //jacIdeal[n][m] += derivGPU(l*DIM+n,k) *
                        //    elUtil.idealMap(elId,k,m * 3 + l);
                        jacIdeal[n][m] += phiM[n][l] *
                            elUtil.idealMap(elId,k,m * 3 + l);               
                    }
                }
            }
            double jacDet = Determinant<DIM>(jacIdeal);
            // end CalcIdealJac

            double I1 = FrobeniusNorm<DIM>(jacIdeal);        
            double sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*ep*ep));

            if(sigma < DBL_MIN && !gradient)
            {
                return DBL_MAX;
            }
            ASSERTL0(sigma > DBL_MIN,"dividing by zero");

            double lsigma = log(sigma);        
            double inc = quadW * absIdealMapDet *
                        (0.5 * mu * (I1 - 3.0 - 2.0*lsigma) +
                         0.5 * K * lsigma * lsigma);
            Kokkos::atomic_add(&grad.integral(node), inc);
            
            // Derivative of basis function in each direction
            if(gradient)
            {
                /*for (int n = 0; n < DIM; ++n)
                {
                    for (int l = 0; ml< DIM; ++l)
                    {
                        //phiM[n][l] = derivGPU(l*DIM+n,k);
                        phiM[n][l] = derivGPU[l*DIM+n];
                    }
                }*/

                double jacInvTrans [DIM][DIM];
                InvTrans<DIM>(phiM, jacInvTrans);
                double derivDet = Determinant<DIM>(phiM);

                double basisDeriv [DIM];
                basisDeriv[0] = derivUtilGPU.VdmD_0(k,localNodeId);
                basisDeriv[1] = derivUtilGPU.VdmD_1(k,localNodeId);
                basisDeriv[2] = derivUtilGPU.VdmD_2(k,localNodeId);

                double jacDetDeriv [DIM];
                for (int m = 0; m < DIM; ++m)
                {
                    jacDetDeriv[m] = 0.0;
                    for (int n = 0; n < DIM; ++n)
                    {
                        jacDetDeriv[m] += jacInvTrans[m][n] * basisDeriv[n];
                    }
                    jacDetDeriv[m] *= derivDet / absIdealMapDet;
                }

                //double jacDeriv [DIM][DIM][DIM];
                /*for (int m = 0; m < DIM; ++m)
                {
                    for (int n = 0; n < DIM; ++n)
                    {
                        double delta = (m == n ? 1.0 : 0.0);
                        for (int l = 0; l < DIM; ++l)
                        {
                            jacDeriv[m][n][l] = delta * basisDeriv[l];
                        }
                    }
                }*/   

                double jacDerivPhi[DIM];                
                for (int n = 0; n < DIM; ++n)
                {
                    jacDerivPhi[n] = 0.0;                                   
                    for (int l = 0; l < DIM; ++l)
                    {
                        // want phi_I^{-1} (l,n)
                        jacDerivPhi[n] += basisDeriv[l] * elUtil.idealMap(elId,k,l + 3*n);
                        /*if (cs == 0 && node == 0 && el == 0 && k==0)
                        {
                            Kokkos::single(Kokkos::PerThread(teamMember),[&] ()
                            {
                                printf("jacDerivPhi[%i][%i]%e\n",m,n, jacDerivPhi[m][n]);
                            });
                        }*/
                    }
                    
                }
                
                
                /*double jacDerivPhi [DIM][DIM][DIM];
                for (int p = 0; p < DIM; ++p)
                {
                    for (int m = 0; m < DIM; ++m)
                    {
                        for (int n = 0; n < DIM; ++n)
                        {
                            jacDerivPhi[p][m][n] = 0.0;
                            for (int l = 0; l < DIM; ++l)
                            {
                                // want phi_I^{-1} (l,n)
                                jacDerivPhi[p][m][n] +=
                                    jacDeriv[p][m][l] * elUtil.idealMap(elId,k,l + 3*n);
                            }
                        }
                    }
                }*/
            
                double frobProd [DIM];
                double frobProd1 [DIM];
                for (int m = 0; m < DIM; ++m)
                {
                    double jacIdeal_m[DIM];
                    frobProd[m] = 0.0;
                    for(int n = 0; n<DIM; n++)
                    {
                        jacIdeal_m[n] = jacIdeal[m][n];
                    }
                    frobProd[m] = ScalarProd<DIM>(jacIdeal_m,jacDerivPhi);

                    ////frobProd[m] = FrobProd<DIM>(jacIdeal,jacDerivPhi[m]);
                    //frobProd[m] = 0.0;
                    for(int n = 0; n<DIM; n++)
                    {
                        //frobProd[m] += jacIdeal[m][n] * jacDerivPhi[n];
                        //frobProd[m] += jacIdeal_m[n] * jacDerivPhi[n];
                        /*if (frobProd1[m] != frobProd[m] && cs ==1 && node == 0 && k== 0)
                        {
                            printf("jacIdeal_m[0] = %e, jacIdeal[%i][0] = %e \t,jacIdeal_m[1] = %e, jacIdeal[%i][1] = %e \t,jacIdeal_m[2] = %e, jacIdeal[%i][2] = %e \t, frobProd[%i] = %e, frobProd1[%i] = %e\n",
                                 jacIdeal_m[0], m, jacIdeal[m][0], jacIdeal_m[1], m, jacIdeal[m][1],jacIdeal_m[2], m, jacIdeal[m][2],
                                 m, frobProd[m], m, frobProd1[m]);
                        }*/
                    }

                }

                for (int m = 0; m < DIM; ++m)
                {
                    double inc = quadW * absIdealMapDet * (
                        mu * frobProd[m] + (jacDetDeriv[m] / (2.0*sigma - jacDet)
                                            * (K * lsigma - mu)));
                    Kokkos::atomic_add(&grad.G(node,m), inc);
                }

                if(hessian)
                {
                    //holder for the hessian frobprods
                    double frobProdHes[DIM][DIM];
                    for (int m = 0; m < DIM; ++m)
                    {
                        for(int l = m; l < DIM; ++l)
                        {
                            //frobProdHes[m][l] = FrobProd<DIM>(jacDerivPhi[m],jacDerivPhi[l]);
                            frobProdHes[m][l] = 0.0;
                            if (m == l)
                            {
                                frobProdHes[m][l] = ScalarProd<DIM>(jacDerivPhi,jacDerivPhi);
                                /*for(int n = 0; n<DIM; n++)
                                { 
                                    frobProdHes[m][l] += jacDerivPhi[n] * jacDerivPhi[n] ;
                                }*/
                            }                        
                        }
                    }

                    int ct = 0;
                    for (int m = 0; m < DIM; ++m)
                    {
                        for(int l = m; l < DIM; ++l, ct++)
                        {
                            double inc = quadW * absIdealMapDet 
                                * (mu * frobProdHes[m][l]
                                + jacDetDeriv[m]*jacDetDeriv[l]/(2.0*sigma-jacDet)/(2.0*sigma-jacDet)
                                    * (K- jacDet*(K*lsigma-mu)/(2.0*sigma-jacDet)));
                            Kokkos::atomic_add(&grad.G(node,ct+DIM), inc);
                        }
                    }
                }
            }
        });
    }

    //});     
       
    return grad.integral(node);

}


} // Utilities

} // Nektar

#endif
