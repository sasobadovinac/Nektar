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


template<int DIM>
inline double Determinant(Kokkos::View<double[DIM][DIM], host_space> &jac)
{
    return 0.0;
}

template<> inline double Determinant<2>(Kokkos::View<double[2][2], host_space> &jac)
{
    return jac(0,0) * jac(1,1) - jac(0,1) * jac(1,0);
}

template<> inline double Determinant<3>(Kokkos::View<double[3][3], host_space> &jac)
{
    return jac(0,0) * (jac(1,1)*jac(2,2) - jac(2,1)*jac(1,2))
          -jac(0,1) * (jac(1,0)*jac(2,2) - jac(1,2)*jac(2,0))
          +jac(0,2) * (jac(1,0)*jac(2,1) - jac(1,1)*jac(2,0));
}


template<int DIM>
inline double LinElasTrace(Kokkos::View<double[DIM][DIM], host_space> &jac)
{
    return 0.0;
}

template<> inline double LinElasTrace<2>(Kokkos::View<double[2][2], host_space> &jac)
{
    return 0.25 * (
        (jac(0,0)*jac(0,0) + jac(1,0)*jac(1,0) - 1.0)*
        (jac(0,0)*jac(0,0) + jac(1,0)*jac(1,0) - 1.0) +
        (jac(0,1)*jac(0,1) + jac(1,1)*jac(1,1) - 1.0)*
        (jac(0,1)*jac(0,1) + jac(1,1)*jac(1,1) - 1.0) )
        + 0.5 * (
            (jac(0,0)*jac(0,1) + jac(1,0)*jac(1,1))*
            (jac(0,0)*jac(0,1) + jac(1,0)*jac(1,1)) );
}

template<> inline double LinElasTrace<3>(Kokkos::View<double[3][3], host_space> &jac)
{
    return 0.25 *(
        (jac(0,0)*jac(0,0)+jac(1,0)*jac(1,0)+jac(2,0)*jac(2,0)-1.0)*
        (jac(0,0)*jac(0,0)+jac(1,0)*jac(1,0)+jac(2,0)*jac(2,0)-1.0) +
        (jac(0,1)*jac(0,1)+jac(1,1)*jac(1,1)+jac(2,1)*jac(2,1)-1.0)*
        (jac(0,1)*jac(0,1)+jac(1,1)*jac(1,1)+jac(2,1)*jac(2,1)-1.0) +
        (jac(0,2)*jac(0,2)+jac(1,2)*jac(1,2)+jac(2,2)*jac(2,2)-1.0)*
        (jac(0,2)*jac(0,2)+jac(1,2)*jac(1,2)+jac(2,2)*jac(2,2)-1.0) )
        + 0.5 * (
            (jac(0,0)*jac(0,2)+jac(1,0)*jac(1,2)+jac(2,0)*jac(2,2))*
            (jac(0,0)*jac(0,2)+jac(1,0)*jac(1,2)+jac(2,0)*jac(2,2)) +
            (jac(0,1)*jac(0,2)+jac(1,1)*jac(1,2)+jac(2,1)*jac(2,2))*
            (jac(0,1)*jac(0,2)+jac(1,1)*jac(1,2)+jac(2,1)*jac(2,2)) +
            (jac(0,0)*jac(0,1)+jac(1,0)*jac(1,1)+jac(0,1)*jac(2,1))*
            (jac(0,0)*jac(0,1)+jac(1,0)*jac(1,1)+jac(0,1)*jac(2,1)) );
}


template<int DIM>
inline void InvTrans(Kokkos::View<double[DIM][DIM], host_space> &in,
                                       Kokkos::View<double[DIM][DIM], host_space> &out)
{
}

template<>
inline void InvTrans<2>(Kokkos::View<double[2][2], host_space> &in, Kokkos::View<double[2][2], host_space> &out)
{
    double invDet = 1.0 / Determinant(in);
    out(0,0) =  in(1,1) * invDet;
    out(1,0) = -in(0,1) * invDet;
    out(0,1) = -in(1,0) * invDet;
    out(1,1) =  in(0,0) * invDet;
}

template<>
inline void InvTrans<3>(Kokkos::View<double[3][3], host_space> &in, Kokkos::View<double[3][3], host_space> &out)
{
    double invdet = 1.0 / Determinant(in);
    out(0,0) =  (in(1,1)*in(2,2)-in(2,1)*in(1,2))*invdet;
    out(1,0) = -(in(0,1)*in(2,2)-in(0,2)*in(2,1))*invdet;
    out(2,0) =  (in(0,1)*in(1,2)-in(0,2)*in(1,1))*invdet;
    out(0,1) = -(in(1,0)*in(2,2)-in(1,2)*in(2,0))*invdet;
    out(1,1) =  (in(0,0)*in(2,2)-in(0,2)*in(2,0))*invdet;
    out(2,1) = -(in(0,0)*in(1,2)-in(1,0)*in(0,2))*invdet;
    out(0,2) =  (in(1,0)*in(2,1)-in(2,0)*in(1,1))*invdet;
    out(1,2) = -(in(0,0)*in(2,1)-in(2,0)*in(0,1))*invdet;
    out(2,2) =  (in(0,0)*in(1,1)-in(1,0)*in(0,1))*invdet;
}


template<int DIM>
inline double FrobProd(Kokkos::View<double[DIM][DIM], host_space> &in1,
                                            Kokkos::View<double[DIM][DIM], host_space> &in2)
{
    double ret = 0;
    for (int n = 0; n < DIM; ++n)
    {
        for (int l = 0; l < DIM; ++l)
        {
            ret += in1(n,l)* in2(n,l);
        }
    }
    return ret;
}

template<int DIM>
inline double FrobeniusNorm(Kokkos::View<double[DIM][DIM], host_space> &inarray)
{
    double ret = 0.0;
    double *start = &inarray(0,0);
    for (int i = 0; i < DIM*DIM; ++i, ++start)
    {
        ret += (*start) * (*start);
    }
    return ret;
}


template<int DIM>
inline double CalcIdealJac(int iD, int point,
                  typename Kokkos::View< double**>::HostMirror &deriv,                  
                  Kokkos::View<double[DIM][DIM], host_space> &jacIdeal,
                  ElUtilGPU &elUtil)
{
    //int iD = data[elmt]->GetId();
    for (int m = 0; m < DIM; ++m)
    {
        for (int n = 0; n < DIM; ++n)
        {
            jacIdeal(n,m) = 0.0;
            for (int l = 0; l < DIM; ++l)
            {
                jacIdeal(n,m) += deriv(l*DIM+n,point) *
                    elUtil.h_idealMap(iD,point,m * 3 + l);                
            }
        }
    }
    return Determinant(jacIdeal);
}


template<int DIM>
double NodeOpti::GetFunctional(DerivUtilGPU &derivUtilGPU,
         NodesGPU &nodes, ElUtilGPU &elUtil, NodeMap &nodeMap, 
         Kokkos::View<double[DIM == 2 ? 5 : 9], host_space> &G,
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
        }            
        coeffs++;
    }

    double minJac = CalcMinJac(elUtil, nElmt, elIdArray);
    double ep = minJac < 0.0 ? sqrt(1e-9 + 0.04*minJac*minJac) : sqrt(1e-9);    
        
    const double nu = 0.4;
    const double mu = 1.0 / 2.0 / (1.0+nu);
    const double K  = 1.0 / 3.0 / (1.0 - 2.0 * nu); 

    //double integral = 0.0;
    Kokkos::View<double[1]> integral("integral");
    typename Kokkos::View< double[1]>::HostMirror h_integral = Kokkos::create_mirror_view(integral);
    h_integral[0] = 0.0; 

    // Storage for derivatives, ordered by:
    //   - standard coordinate direction
    //   - cartesian coordinate direction, combined
    //   - quadrature points    
    Kokkos::View<double**> derivGPU("derivGPU", DIM*DIM, ptsHighGPU);
    typename Kokkos::View< double**>::HostMirror h_derivGPU = Kokkos::create_mirror_view(derivGPU);

    for (int el = 0; el < nElmt; ++el)
    {   
        int elId = elIdArray[el];
        int localNodeId = localNodeIdArray[el];

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
        
    
        //for(int k = 0; k < ptsHighGPU; ++k)
        Kokkos::parallel_for (range_policy_host(0,ptsHighGPU), KOKKOS_LAMBDA (const int k)
        {
            double absIdealMapDet = fabs(elUtil.h_idealMap(elId,k,9));
            double quadW = derivUtilGPU.h_quadW(k);

            Kokkos::View<double[DIM][DIM], host_space> jacIdeal("jacIdeal");
            //double jacDet = CalcIdealJac(elId, k, h_derivGPU, jacIdeal, elUtil);
            for (int m = 0; m < DIM; ++m)
            {
                for (int n = 0; n < DIM; ++n)
                {
                    jacIdeal(n,m) = 0.0;
                    for (int l = 0; l < DIM; ++l)
                    {
                        jacIdeal(n,m) += h_derivGPU(l*DIM+n,k) *
                            elUtil.h_idealMap(elId,k,m * 3 + l);                
                    }
                }
            }
            double jacDet = Determinant(jacIdeal);
            // end CalcIdealJac
            
            double I1 = FrobeniusNorm(jacIdeal);

            double sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*ep*ep));
            double lsigma = log(sigma);

            double inc = quadW * absIdealMapDet *
                        (0.5 * mu * (I1 - 3.0 - 2.0*lsigma) +
                         0.5 * K * lsigma * lsigma);
            Kokkos::atomic_add(&h_integral[0], inc);
            


            // Derivative of basis function in each direction
            if(gradient)
            {
                Kokkos::View<double[DIM][DIM], host_space> phiM("phiM");
                for (int m = 0; m < DIM; ++m)
                {
                    for (int n = 0; n < DIM; ++n)
                    {
                        //phiM[n][m] = h_derivGPU(m,n,k);
                        phiM(n,m) = h_derivGPU(m*DIM+n,k);
                    }
                }

                Kokkos::View<double[DIM][DIM], host_space> jacInvTrans("jacInvTrans");
                InvTrans<DIM>(phiM, jacInvTrans);
                double derivDet = Determinant<DIM>(phiM);

                Kokkos::View<double[DIM], host_space> basisDeriv("basisDeriv");
                basisDeriv(0) = derivUtilGPU.h_VdmD_0(k,localNodeId);
                basisDeriv(1) = derivUtilGPU.h_VdmD_1(k,localNodeId);
                basisDeriv(2) = derivUtilGPU.h_VdmD_2(k,localNodeId);

                Kokkos::View<double[DIM], host_space> jacDetDeriv("jacDetDeriv");
                for (int m = 0; m < DIM; ++m)
                {
                    jacDetDeriv(m) = 0.0;
                    for (int n = 0; n < DIM; ++n)
                    {
                        jacDetDeriv(m) += jacInvTrans(m,n) * basisDeriv(n);
                    }
                    jacDetDeriv(m) *= derivDet / absIdealMapDet;
                }

                Kokkos::View<double[DIM][DIM][DIM], host_space> jacDeriv("jacDeriv");
                for (int m = 0; m < DIM; ++m)
                {
                    for (int n = 0; n < DIM; ++n)
                    {
                        double delta = (m == n ? 1.0 : 0.0);
                        for (int l = 0; l < DIM; ++l)
                        {
                            jacDeriv(m,n,l) = delta * basisDeriv(l);
                        }
                    }
                }

                Kokkos::View<double[DIM][DIM][DIM], host_space> jacDerivPhi("jacDerivPhi");
                for (int p = 0; p < DIM; ++p)
                {
                    for (int m = 0; m < DIM; ++m)
                    {
                        for (int n = 0; n < DIM; ++n)
                        {
                            jacDerivPhi(p,m,n) = 0.0;
                            for (int l = 0; l < DIM; ++l)
                            {
                                // want phi_I^{-1} (l,n)
                                jacDerivPhi(p,m,n) +=
                                    jacDeriv(p,m,l) * elUtil.h_idealMap(elId,k,l + 3*n);
                            }
                        }
                    }
                }

                Kokkos::View<double[DIM], host_space> frobProd("frobProd");
                for (int m = 0; m < DIM; ++m)
                {
                    Kokkos::View<double[DIM][DIM], host_space> jacDerivPhi_sub 
                            = Kokkos::subview(jacDerivPhi, m, Kokkos::ALL(), Kokkos::ALL());
                    frobProd(m) = FrobProd<DIM>(jacIdeal,jacDerivPhi_sub);
                }

                for (int j = 0; j < DIM; ++j)
                {
                    double inc = quadW * absIdealMapDet * (
                        mu * frobProd(j) + (jacDetDeriv(j) / (2.0*sigma - jacDet)
                                            * (K * lsigma - mu)));
                    Kokkos::atomic_add(&G(j), inc);
                }

                if(hessian)
                {
                    Kokkos::View<double[DIM][DIM], host_space> frobProdHes("frobProdHes"); //holder for the hessian frobprods
                    for (int m = 0; m < DIM; ++m)
                    {
                        for(int l = m; l < DIM; ++l)
                        {
                            Kokkos::View<double[DIM][DIM], host_space> jacDerivPhi_m 
                                = Kokkos::subview(jacDerivPhi, m, Kokkos::ALL(), Kokkos::ALL());
                            Kokkos::View<double[DIM][DIM], host_space> jacDerivPhi_l 
                                = Kokkos::subview(jacDerivPhi, l, Kokkos::ALL(), Kokkos::ALL());
                            frobProdHes(m,l) = FrobProd<DIM>(jacDerivPhi_m,jacDerivPhi_l);
                        }
                    }


                    int ct = 0;
                    for (int m = 0; m < DIM; ++m)
                    {
                        for(int l = m; l < DIM; ++l, ct++)
                        {
                            double inc = quadW * absIdealMapDet * (
                                mu * frobProdHes(m,l) +
                                jacDetDeriv(m)*jacDetDeriv(l)/(2.0*sigma-jacDet)/(2.0*sigma-jacDet)*(
                                    K- jacDet*(K*lsigma-mu)/(2.0*sigma-jacDet)));
                            Kokkos::atomic_add(&G(ct+DIM), inc);
                        }
                    }
                }
            }
        });
    } 

    return h_integral[0];
}

}
}

#endif
