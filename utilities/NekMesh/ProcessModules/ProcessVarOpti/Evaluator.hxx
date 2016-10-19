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


/*template<int DIM>
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
}*/
template<int DIM> inline NekDouble FrobProd(NekDouble in1[DIM][DIM],
                                            NekDouble in2[DIM][DIM])
{
    NekDouble ret = 0;
    for (int n = 0; n < DIM; ++n)
    {
        for (int l = 0; l < DIM; ++l)
        {
            ret += in1[n][l] * in2[n][l];
        }
    }
    return ret;
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
NekDouble NodeOpti::GetFunctional(DerivUtilGPU &derivUtilGPU,
         NodesGPU &nodes, ElUtilGPU &elUtil, NodeMap &nodeMap, 
         Grad &grad,
         bool gradient, bool hessian)
{
    LibUtilities::ShapeType st = data[0]->GetEl()->GetShapeType(); // only for derivUtil->quadW
    const int nElmt = data.size();    
    const int ptsLowGPU = derivUtilGPU.ptsLow;
    const int ptsHighGPU = derivUtilGPU.ptsHigh;

    int elIdArray[nElmt];
    // using the node ID and the node map find the local coordinates of the node,
                // depending on the considered element
    int globalNodeId = node->m_id;
    int localNodeIdArray[nElmt];
    NodeMap::const_iterator coeffs;
    coeffs = nodeMap.find(globalNodeId); 
    for (int el = 0; el < nElmt; ++el)
    {
        elIdArray[el] = data[el]->GetId();

        int elmt = std::get<0>(coeffs->second);
        if (elmt == elIdArray[el])
        {
            localNodeIdArray[el] = std::get<1>(coeffs->second);
            //printf("localNodeIdArray[%i] = %i\n", el, localNodeIdArray[el]);
        }            
        coeffs++;
    }

    double minJac = CalcMinJac(elUtil, nElmt, elIdArray);
    double ep = minJac < 0.0 ? sqrt(1e-9 + 0.04*minJac*minJac) : sqrt(1e-9);    
        
    const double nu = 0.4;
    const double mu = 1.0 / 2.0 / (1.0+nu);
    const double K  = 1.0 / 3.0 / (1.0 - 2.0 * nu); 

    
    // Storage for derivatives, ordered by:
    //   - standard coordinate direction
    //   - cartesian coordinate direction, combined
    //   - quadrature points    
    Kokkos::View<double**> derivGPU("derivGPU", DIM*DIM, ptsHighGPU);
    typename Kokkos::View< double**>::HostMirror h_derivGPU = Kokkos::create_mirror_view(derivGPU);

    grad.h_integral[0] = 0.0; 


    for (int el = 0; el < nElmt; ++el)
    {   
        int elId = elIdArray[el];
        int localNodeId = localNodeIdArray[el];
        if (node->m_id == 165  || node->m_id == 3905 || node->m_id == 189)
        {
            printf("elId = %i\n", elId);
            printf("localNodeId = %i\n", localNodeId);
        }

        Kokkos::parallel_for (range_policy(0,ptsHighGPU), KOKKOS_LAMBDA (const int j)
        {
            derivGPU(0,j) = 0.0;
            derivGPU(1,j) = 0.0;
            derivGPU(2,j) = 0.0;

            derivGPU(3,j) = 0.0;
            derivGPU(4,j) = 0.0;
            derivGPU(5,j) = 0.0;

            derivGPU(6,j) = 0.0;
            derivGPU(7,j) = 0.0;
            derivGPU(8,j) = 0.0;
            for (int k = 0; k < ptsLowGPU; ++k)
            {
                derivGPU(0,j) += derivUtilGPU.VdmD_0(j,k) * nodes.X(elId,k);
                derivGPU(1,j) += derivUtilGPU.VdmD_0(j,k) * nodes.Y(elId,k);
                derivGPU(2,j) += derivUtilGPU.VdmD_0(j,k) * nodes.Z(elId,k);

                derivGPU(3,j) += derivUtilGPU.VdmD_1(j,k) * nodes.X(elId,k);
                derivGPU(4,j) += derivUtilGPU.VdmD_1(j,k) * nodes.Y(elId,k);
                derivGPU(5,j) += derivUtilGPU.VdmD_1(j,k) * nodes.Z(elId,k);

                derivGPU(6,j) += derivUtilGPU.VdmD_2(j,k) * nodes.X(elId,k);
                derivGPU(7,j) += derivUtilGPU.VdmD_2(j,k) * nodes.Y(elId,k);
                derivGPU(8,j) += derivUtilGPU.VdmD_2(j,k) * nodes.Z(elId,k);
            }
        });

        Kokkos::deep_copy(h_derivGPU,derivGPU);

        Kokkos::parallel_for (range_policy_host(0,ptsHighGPU), KOKKOS_LAMBDA (const int k)
        {
            
            double absIdealMapDet = fabs(elUtil.h_idealMap(elId,k,9));
            double quadW = derivUtilGPU.h_quadW(k);
            double jacIdeal[3][3];
            for (int m = 0; m < DIM; ++m)
            {
                for (int n = 0; n < DIM; ++n)
                {
                    jacIdeal[n][m] = 0.0;
                    for (int l = 0; l < DIM; ++l)
                    {
                        jacIdeal[n][m] += h_derivGPU(l*DIM+n,k) *
                            elUtil.h_idealMap(elId,k,m * 3 + l);               
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
            Kokkos::atomic_add(&grad.h_integral[0], inc);
            //printf("inc = %e\n", inc);
            //printf("grad.h_integral[0] = %e\n", grad.h_integral[0]);

            // Derivative of basis function in each direction
            if(gradient)
            {
                double phiM [DIM][DIM]; 
                for (int m = 0; m < DIM; ++m)
                {
                    for (int n = 0; n < DIM; ++n)
                    {
                        //phiM[n][m] = derivGPU(m,n,k);
                        phiM[n][m] = h_derivGPU(m*DIM+n,k);
                    }
                }

                double jacInvTrans [DIM][DIM];
                InvTrans<DIM>(phiM, jacInvTrans);
                double derivDet = Determinant<DIM>(phiM);

                double basisDeriv [DIM];
                if (localNodeId != nodeIds[el])
                    printf("%s%i\n", "localNodeId, " ,localNodeId);
                basisDeriv[0] = derivUtilGPU.h_VdmD_0(k,localNodeId);
                basisDeriv[1] = derivUtilGPU.h_VdmD_1(k,localNodeId);
                basisDeriv[2] = derivUtilGPU.h_VdmD_2(k,localNodeId);

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
                                    jacDeriv[p][m][l] * elUtil.h_idealMap(elId,k,l + 3*n);
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
                    Kokkos::atomic_add(&grad.h_G(j), inc);
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
                            Kokkos::atomic_add(&grad.h_G(ct+DIM), inc);
                        }
                    }
                }
            }
        });
    } 
    //Kokkos::deep_copy(grad.h_integral, grad.integral);
    //Kokkos::deep_copy(grad.h_G, grad.G);

    //printf("grad.h_integral[0] = %e\n", grad.h_integral[0]);

    return grad.h_integral[0];
}

}
}

#endif
