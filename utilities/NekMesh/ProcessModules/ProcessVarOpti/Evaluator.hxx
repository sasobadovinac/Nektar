////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.h
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

#ifndef UTILITIES_NEKMESH_NODEOPTI_EVALUATOR
#define UTILITIES_NEKMESH_NODEOPTI_EVALUATOR

namespace Nektar
{
namespace Utilities
{

typedef Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type  member_type;

template<int DIM> inline NekDouble Determinant(NekDouble jac[DIM][DIM])
{
    return 0.0;
}

template<> inline NekDouble Determinant<2>(NekDouble jac[2][2])
{
    return jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
}

template<> inline NekDouble Determinant<3>(NekDouble jac[3][3])
{
    return jac[0][0] * (jac[1][1]*jac[2][2] - jac[2][1]*jac[1][2])
          -jac[0][1] * (jac[1][0]*jac[2][2] - jac[1][2]*jac[2][0])
          +jac[0][2] * (jac[1][0]*jac[2][1] - jac[1][1]*jac[2][0]);
}


template<int DIM> inline NekDouble LinElasTrace(NekDouble jac[DIM][DIM])
{
    return 0.0;
}

template<> inline NekDouble LinElasTrace<2>(NekDouble jac[2][2])
{
    return 0.25 * (
        (jac[0][0]*jac[0][0] + jac[1][0]*jac[1][0] - 1.0) *
        (jac[0][0]*jac[0][0] + jac[1][0]*jac[1][0] - 1.0) +
        (jac[0][1]*jac[0][1] + jac[1][1]*jac[1][1] - 1.0)*
        (jac[0][1]*jac[0][1] + jac[1][1]*jac[1][1] - 1.0))
        + 0.5 * (
            (jac[0][0]*jac[0][1] + jac[1][0]*jac[1][1])*
            (jac[0][0]*jac[0][1] + jac[1][0]*jac[1][1]));
}

template<> inline NekDouble LinElasTrace<3>(NekDouble jac[3][3])
{
    return 0.25 *(
        (jac[0][0]*jac[0][0]+jac[1][0]*jac[1][0]+jac[2][0]*jac[2][0]-1.0)*
        (jac[0][0]*jac[0][0]+jac[1][0]*jac[1][0]+jac[2][0]*jac[2][0]-1.0) +
        (jac[0][1]*jac[0][1]+jac[1][1]*jac[1][1]+jac[2][1]*jac[2][1]-1.0)*
        (jac[0][1]*jac[0][1]+jac[1][1]*jac[1][1]+jac[2][1]*jac[2][1]-1.0) +
        (jac[0][2]*jac[0][2]+jac[1][2]*jac[1][2]+jac[2][2]*jac[2][2]-1.0)*
        (jac[0][2]*jac[0][2]+jac[1][2]*jac[1][2]+jac[2][2]*jac[2][2]-1.0))
        + 0.5 * (
            (jac[0][0]*jac[0][2]+jac[1][0]*jac[1][2]+jac[2][0]*jac[2][2])*
            (jac[0][0]*jac[0][2]+jac[1][0]*jac[1][2]+jac[2][0]*jac[2][2])+
            (jac[0][1]*jac[0][2]+jac[1][1]*jac[1][2]+jac[2][1]*jac[2][2])*
            (jac[0][1]*jac[0][2]+jac[1][1]*jac[1][2]+jac[2][1]*jac[2][2])+
            (jac[0][0]*jac[0][1]+jac[1][0]*jac[1][1]+jac[0][1]*jac[2][1])*
            (jac[0][0]*jac[0][1]+jac[1][0]*jac[1][1]+jac[0][1]*jac[2][1]));
}

template<int DIM>
inline void InvTrans(NekDouble in[DIM][DIM],
                                       NekDouble out[DIM][DIM])
{
}

template<>
inline void InvTrans<2>(NekDouble in[2][2], NekDouble out[2][2])
{
    NekDouble invDet = 1.0 / Determinant(in);
    out[0][0] =  in[1][1] * invDet;
    out[1][0] = -in[0][1] * invDet;
    out[0][1] = -in[1][0] * invDet;
    out[1][1] =  in[0][0] * invDet;
}

template<>
inline void InvTrans<3>(NekDouble in[3][3], NekDouble out[3][3])
{
    NekDouble invdet = 1.0 / Determinant(in);
    out[0][0] =  (in[1][1]*in[2][2]-in[2][1]*in[1][2])*invdet;
    out[1][0] = -(in[0][1]*in[2][2]-in[0][2]*in[2][1])*invdet;
    out[2][0] =  (in[0][1]*in[1][2]-in[0][2]*in[1][1])*invdet;
    out[0][1] = -(in[1][0]*in[2][2]-in[1][2]*in[2][0])*invdet;
    out[1][1] =  (in[0][0]*in[2][2]-in[0][2]*in[2][0])*invdet;
    out[2][1] = -(in[0][0]*in[1][2]-in[1][0]*in[0][2])*invdet;
    out[0][2] =  (in[1][0]*in[2][1]-in[2][0]*in[1][1])*invdet;
    out[1][2] = -(in[0][0]*in[2][1]-in[2][0]*in[0][1])*invdet;
    out[2][2] =  (in[0][0]*in[1][1]-in[1][0]*in[0][1])*invdet;
}


template<int DIM>
inline NekDouble FrobProd(NekDouble in1[DIM][DIM],
                          NekDouble in2[DIM][DIM])
{
    return 0.0;
}

template<>
inline NekDouble FrobProd<2>(NekDouble in1[2][2], NekDouble in2[2][2])
{
    return    in1[0][0] * in2[0][0]
            + in1[0][1] * in2[0][1]
            + in1[1][0] * in2[1][0]
            + in1[1][1] * in2[1][1] ;
}

template<>
inline NekDouble FrobProd<3>(NekDouble in1[3][3], NekDouble in2[3][3])
{
    return    in1[0][0] * in2[0][0]
            + in1[0][1] * in2[0][1]
            + in1[0][2] * in2[0][2]
            + in1[1][0] * in2[1][0]
            + in1[1][1] * in2[1][1]
            + in1[1][2] * in2[1][2]
            + in1[2][0] * in2[2][0]
            + in1[2][1] * in2[2][1]
            + in1[2][2] * in2[2][2] ;
}


template<int DIM>
inline NekDouble FrobeniusNorm(NekDouble inarray[DIM][DIM])
{
    return 0.0;
}

template<>
inline NekDouble FrobeniusNorm<2>(NekDouble inarray[2][2])
{
    return    inarray[0][0] * inarray[0][0]
            + inarray[0][1] * inarray[0][1]
            + inarray[1][0] * inarray[1][0]
            + inarray[1][1] * inarray[1][1] ;
}

template<>
inline NekDouble FrobeniusNorm<3>(NekDouble inarray[3][3])
{
    return    inarray[0][0] * inarray[0][0]
            + inarray[0][1] * inarray[0][1]
            + inarray[0][2] * inarray[0][2]
            + inarray[1][0] * inarray[1][0]
            + inarray[1][1] * inarray[1][1]
            + inarray[1][2] * inarray[1][2]
            + inarray[2][0] * inarray[2][0]
            + inarray[2][1] * inarray[2][1]
            + inarray[2][2] * inarray[2][2] ;
}




template<int DIM>
KOKKOS_INLINE_FUNCTION
NekDouble NodeOpti::GetFunctional(const DerivUtilGPU &derivUtilGPU,
         const NodesGPU &nodes, const ElUtilGPU &elUtil, 
         const Grad &grad, int nElmt, //const int elId, const int localNodeId,
         const double ep, const member_type &teamMember,
         bool gradient, bool hessian)
{ 
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
    
    //grad.h_integral[0] = 0.0;
    //Kokkos::deep_copy(grad.integral,grad.h_integral);
grad.integral[0] = 0.0;
for (int el = 0; el < nElmt; ++el)
{
//Kokkos::parallel_for( team_policy( nElmt, Kokkos::AUTO ), KOKKOS_LAMBDA ( const member_type& teamMember)
//{ 
    int elId = nodes.elIdArray[el];
    int localNodeId = nodes.localNodeIdArray[el];   
    
    Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsHighGPU ), [&] ( const int k)
    {        
    //Kokkos::parallel_for (range_policy(0,ptsHighGPU), KOKKOS_LAMBDA (const int k)
    //{
    //for (int k = 0; k < ptsHighGPU; ++k)
    //{
        /*derivGPU(0,k) = 0.0;
        derivGPU(1,k) = 0.0;
        derivGPU(2,k) = 0.0;

        derivGPU(3,k) = 0.0;
        derivGPU(4,k) = 0.0;
        derivGPU(5,k) = 0.0;

        derivGPU(6,k) = 0.0;
        derivGPU(7,k) = 0.0;
        derivGPU(8,k) = 0.0;
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
        double jacIdeal[3][3];
        for (int m = 0; m < DIM; ++m)
        {
            for (int n = 0; n < DIM; ++n)
            {
                jacIdeal[n][m] = 0.0;
                for (int l = 0; l < DIM; ++l)
                {
                    //jacIdeal[n][m] += derivGPU(l*DIM+n,k) *
                    //    elUtil.idealMap(elId,k,m * 3 + l);
                    jacIdeal[n][m] += derivGPU[l*DIM+n] *
                        elUtil.idealMap(elId,k,m * 3 + l);               
                }
            }
        }
        double jacDet = Determinant<DIM>(jacIdeal);
        // end CalcIdealJac

        double I1 = FrobeniusNorm<DIM>(jacIdeal);        
        double sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*ep*ep));
        double lsigma = log(sigma);        
        double inc = quadW * absIdealMapDet *
                    (0.5 * mu * (I1 - 3.0 - 2.0*lsigma) +
                     0.5 * K * lsigma * lsigma);
        Kokkos::atomic_add(&grad.integral[0], inc);
        
        // Derivative of basis function in each direction
        if(gradient)
        {
            double phiM [DIM][DIM]; 
            for (int m = 0; m < DIM; ++m)
            {
                for (int n = 0; n < DIM; ++n)
                {
                    //phiM[n][m] = derivGPU(m*DIM+n,k);
                    phiM[n][m] = derivGPU[m*DIM+n];
                }
            }

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

            double jacDeriv [DIM][DIM][DIM];
            for (int m = 0; m < DIM; ++m)
            {
                for (int n = 0; n < DIM; ++n)
                {
                    double delta = (m == n ? 1.0 : 0.0);
                    for (int l = 0; l < DIM; ++l)
                    {
                        jacDeriv[m][n][l] = delta * basisDeriv[l];
                    }
                }
            }

            double jacDerivPhi [DIM][DIM][DIM];
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
            }
        
            double frobProd [DIM];
            for (int m = 0; m < DIM; ++m)
            {
                frobProd[m] = FrobProd<DIM>(jacIdeal,jacDerivPhi[m]);
            }

            for (int j = 0; j < DIM; ++j)
            {
                double inc = quadW * absIdealMapDet * (
                    mu * frobProd[j] + (jacDetDeriv[j] / (2.0*sigma - jacDet)
                                        * (K * lsigma - mu)));
                Kokkos::atomic_add(&grad.G(j), inc);
            }

            if(hessian)
            {
                //holder for the hessian frobprods
                double frobProdHes[DIM][DIM];
                for (int m = 0; m < DIM; ++m)
                {
                    for(int l = m; l < DIM; ++l)
                    {
                        frobProdHes[m][l] = FrobProd<DIM>(jacDerivPhi[m],jacDerivPhi[l]);
                    }
                }

                int ct = 0;
                for (int m = 0; m < DIM; ++m)
                {
                    for(int l = m; l < DIM; ++l, ct++)
                    {
                        double inc = quadW * absIdealMapDet * (
                            mu * frobProdHes[m][l] +
                            jacDetDeriv[m]*jacDetDeriv[l]/(2.0*sigma-jacDet)/(2.0*sigma-jacDet)*(
                                K- jacDet*(K*lsigma-mu)/(2.0*sigma-jacDet)));
                        Kokkos::atomic_add(&grad.G(ct+DIM), inc);
                    }
                }
            }
        }
    });
}

    //});
     
    //Kokkos::deep_copy(grad.h_integral, grad.integral);
        
    return grad.integral[0];

}

}
}

#endif
