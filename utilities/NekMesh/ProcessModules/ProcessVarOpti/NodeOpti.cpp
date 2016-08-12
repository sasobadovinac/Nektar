////////////////////////////////////////////////////////////////////////////////
//
//  File: NodeOpti.cpp
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

#include <boost/multi_array.hpp>

#include "NodeOpti.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

NodeOptiFactory &GetNodeOptiFactory()
{
    /*
    typedef Loki::SingletonHolder<NodeOptiFactory, Loki::CreateUsingNew,
                                  Loki::NoDestroy, Loki::ClassLevelLockable> Type;
    return Type::Instance();
    */
    //typedef Loki::SingletonHolder<NodeOptiFactory, Loki::CreateUsingNew,
    //                              Loki::NoDestroy, Loki::ClassLevelLockable> Type;
    //return Type::Instance();

    static NodeOptiFactory asd;
    return asd;
}

void NodeOpti::CalcDX()
{
    dx = numeric_limits<double>::max();

    for(int i = 0; i < data.size(); i++)
    {
        dx = min(dx, data[i]->delta);
    }
}

void NodeOpti::CalcMinJac()
{
    minJac = numeric_limits<double>::max();

    for(int i = 0; i < data.size(); i++)
    {
        minJac = min(minJac, data[i]->minJac);
    }
}

int NodeOpti2D2D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    22, NodeOpti2D2D::create, "2D2D");

void NodeOpti2D2D::Optimise()
{
    CalcDX();

    CalcMinJac();

    Array<OneD, NekDouble> G = GetGrad();
    Array<OneD, NekDouble> GA = GetGrad(true);

    //NekDouble tmp[5] = { G[0], G[1], G[2], G[3], G[4] };
    G[0] = GA[0];
    G[1] = GA[1];
    G[2] = GA[2];
    G[3] = GA[4];
    G[4] = GA[3];

    if(G[0]*G[0] + G[1]*G[1] > 1e-12)
    {
        // cout << endl;
        // cout << "approx " << G[0] << " " << G[1] << endl;
        // cout << "approx " << G[2] << " " << G[3] << " " << G[4] << endl;
        // cout << "analytic " << GA[0] << " " << GA[1] << endl;
        // cout << "analytic " << GA[2] << " " << GA[3] << " " << GA[4] << endl;

        //needs to optimise
        NekDouble currentW = GetFunctional<2>();
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;
        NekDouble alpha    = 1.0;
        NekDouble delX = 1.0/(G[2]*G[3]-G[4]*G[4])*(G[3]*G[0] - G[4]*G[1]);
        NekDouble delY = 1.0/(G[2]*G[3]-G[4]*G[4])*(G[2]*G[1] - G[4]*G[0]);

        bool found = false;
        while(alpha > 1e-10)
        {
            node->m_x = xc - alpha * delX;
            node->m_y = yc - alpha * delY;
            if(GetFunctional<2>() < currentW)
            {
                found = true;
                break;
            }

            alpha /= 2.0;
        }

        if(!found)
        {
            //reset the node
            node->m_x = xc;
            node->m_y = yc;
            //cout << "warning: had to reset node " << " " << nodeIds[0] << endl;
            //cout << "    approx grad  : " << tmp[0] << " " << tmp[1] << endl;
            // cout << "    approx hess  : " << tmp[2] << " " << tmp[3] << " " << tmp[4] << endl;
            //cout << "    analytic grad: " << GA[0] << " " << GA[1] << endl;
            // cout << "    analytic hess: " << GA[2] << " " << GA[3] << " " << GA[4] << endl;
        }
        mtx.lock();
        res->val = max(sqrt((node->m_x-xc)*(node->m_x-xc)+(node->m_y-yc)*(node->m_y-yc)),res->val);
        res->nReset++;
        mtx.unlock();
    }
}

Array<OneD, NekDouble> NodeOpti2D2D::GetGrad(bool analytic)
{

    if(analytic)
    {
        GetFunctional<2>(true);
        return grad;
    }

    NekDouble xc = node->m_x;
    NekDouble yc = node->m_y;
    vector<NekDouble> w(9);

    for(int i = 0; i < 7; i++)
    {
        node->m_x = xc + dir[i][0] * dx;
        node->m_y = yc + dir[i][1] * dx;
        w[i] = GetFunctional<2>();
    }
    node->m_x = xc;
    node->m_y = yc;

    Array<OneD, NekDouble> ret(5,0.0);

    //ret[0] d/dx
    //ret[1] d/dy

    //ret[3] d2/dx2
    //ret[4] d2/dxdy
    //ret[5] d2/dy2

    ret[0] = (w[1] - w[4]) / 2.0 / dx;
    ret[1] = (w[3] - w[6]) / 2.0 / dx;

    ret[2] = (w[1] + w[4] - 2.0*w[0]) / dx / dx;
    ret[3] = (w[3] + w[6] - 2.0*w[0]) / dx / dx;
    ret[4] = (w[2] - w[1] - w[3] + 2.0*w[0] - w[4] - w[6] + w[5]) / 2.0 / dx / dx;

    return ret;
}

int NodeOpti3D3D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    33, NodeOpti3D3D::create, "3D3D");

void NodeOpti3D3D::Optimise()
{
    CalcDX();

    CalcMinJac();

    Array<OneD, NekDouble> G = GetGrad();

    if(sqrt(G[0]*G[0] + G[1]*G[1] + G[2]*G[2]) > 1e-10)
    {
        //needs to optimise
        NekDouble currentW = GetFunctional<3>();
        NekDouble functional;
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;
        NekDouble zc       = node->m_z;
        NekDouble alpha    = 1.0;
        NekDouble delX;
        NekDouble delY;
        NekDouble delZ;

        NekDouble det = G[3]*(G[4]*G[5]-G[8]*G[8])
                       -G[6]*(G[6]*G[5]-G[7]*G[8])
                       +G[7]*(G[6]*G[8]-G[7]*G[4]);

        delX = G[0]*(G[4]*G[5]-G[8]*G[8]) +
               G[1]*(G[7]*G[8]-G[6]*G[5]) +
               G[2]*(G[6]*G[8]-G[7]*G[4]);
        delY = G[0]*(G[8]*G[7]-G[6]*G[5]) +
               G[1]*(G[3]*G[5]-G[7]*G[7]) +
               G[2]*(G[6]*G[7]-G[3]*G[8]);
        delZ = G[0]*(G[6]*G[8]-G[4]*G[7]) +
               G[1]*(G[6]*G[7]-G[3]*G[8]) +
               G[2]*(G[3]*G[4]-G[6]*G[6]);

        delX /= det;
        delY /= det;
        delZ /= det;

        bool found = false;
        while(alpha > 1e-10)
        {
            node->m_x = xc - alpha * delX;
            node->m_y = yc - alpha * delY;
            node->m_z = zc - alpha * delZ;
            functional = GetFunctional<3>();
            if(functional < currentW)
            {
                found = true;
                break;
            }

            alpha /= 2.0;
        }

        if(!found)
        {
            //reset the node
            node->m_x = xc;
            node->m_y = yc;
            node->m_z = zc;
            //functional = currentW;
            //cout << "warning: had to reset node" << endl;
        }
        mtx.lock();
        res->val = max(sqrt((node->m_x-xc)*(node->m_x-xc)+(node->m_y-yc)*(node->m_y-yc)+
                            (node->m_z-zc)*(node->m_z-zc)),res->val);
        mtx.unlock();
    }
}

Array<OneD, NekDouble> NodeOpti3D3D::GetGrad(bool analytic)
{
    NekDouble xc = node->m_x;
    NekDouble yc = node->m_y;
    NekDouble zc = node->m_z;

    vector<NekDouble> w;

    for(int i = 0; i < 13; i++)
    {
        node->m_x = xc + dir[i][0] * dx;
        node->m_y = yc + dir[i][1] * dx;
        node->m_z = zc + dir[i][2] * dx;
        w.push_back(GetFunctional<3>());
    }
    node->m_x = xc;
    node->m_y = yc;
    node->m_z = zc;

    w[0] = GetFunctional<3>();

    //cout << "ANALYTIC: " << gradient[0] << " " << gradient[1] << " " << gradient[2] << endl;

    Array<OneD, NekDouble> ret(9,0.0);

    //ret[0] d/dx
    //ret[1] d/dy
    //ret[2] d/dz

    //ret[3] d2/dx2
    //ret[4] d2/dy2
    //ret[5] d2/dz2
    //ret[6] d2/dxdy
    //ret[7] d2/dxdz
    //ret[8] d2/dydz

    //ret[0] = gradient[0];
    //ret[1] = gradient[1];
    //ret[2] = gradient[2];

    ret[0] = (w[1] - w[4]) / 2.0 / dx;
    ret[1] = (w[3] - w[6]) / 2.0 / dx;
    ret[2] = (w[9] - w[8]) / 2.0 / dx;

    //cout << "APPROX: " << ret[0] << " " << ret[1] << " " << ret[2] << endl;
    //cout << endl;

    ret[3] = (w[1] + w[4] - 2.0*w[0]) / dx / dx;
    ret[4] = (w[3] + w[6] - 2.0*w[0]) / dx / dx;
    ret[5] = (w[9] + w[8] - 2.0*w[0]) / dx / dx;

    ret[6] = (w[2] - w[1] - w[3] + 2.0*w[0] - w[4] - w[6] + w[5]) / 2.0 / dx / dx;
    ret[7] = (w[10] - w[1] - w[9] + 2.0*w[0] - w[4] - w[8] + w[7]) / 2.0 / dx / dx;
    ret[8] = (w[11] - w[3] - w[9] + 2.0*w[0] - w[6] - w[8] + w[12]) / 2.0 / dx / dx;

    return ret;
}

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

template<int DIM> inline void InvTrans(NekDouble in[DIM][DIM],
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

// Typedef for derivative storage, we use boost::multi_array so we can pass this
// to functions easily
typedef boost::multi_array<NekDouble, 4> DerivArray;

template<int DIM>
inline NekDouble CalcIdealJac(int elmt,
                              int point,
                              DerivArray &deriv,
                              std::vector<ElUtilSharedPtr> &data,
                              NekDouble jacIdeal[DIM][DIM])
{
    for (int m = 0; m < DIM; ++m)
    {
        for (int n = 0; n < DIM; ++n)
        {
            jacIdeal[n][m] = 0.0;
            for (int l = 0; l < DIM; ++l)
            {
                jacIdeal[n][m] += deriv[l][elmt][n][point] *
                    data[elmt]->maps[point][m * 3 + l];
            }
        }
    }

    return Determinant(jacIdeal);
}

template<int DIM>
inline NekDouble FrobeniusNorm(NekDouble inarray[DIM][DIM])
{
    NekDouble ret = 0.0, *start = &inarray[0][0];
    for (int i = 0; i < DIM*DIM; ++i, ++start)
    {
        ret += (*start) * (*start);
    }
    return ret;
}

template<int DIM>
NekDouble NodeOpti::GetFunctional(bool analytic)
{
    const int nElmt  = data.size();
    const int totpts = derivUtil->ptsLow * nElmt;
    NekDouble X[DIM * totpts];

    // Store x/y components of each element sequentially in memory
    for (int i = 0, cnt = 0; i < nElmt; ++i)
    {
        for (int j = 0; j < derivUtil->ptsLow; ++j)
        {
            for (int d = 0; d < DIM; ++d)
            {
                X[cnt + d*derivUtil->ptsLow + j] = *(data[i]->nodes[j][d]);
            }
        }

        cnt += DIM*derivUtil->ptsLow;
    }

    // Storage for derivatives, ordered by:
    //   - standard coordinate direction
    //   - number of elements
    //   - cartesian coordinate direction
    //   - quadrature points
    DerivArray deriv(boost::extents[DIM][nElmt][DIM][derivUtil->ptsHigh]);

    // Calculate x- and y-gradients
    for (int d = 0; d < DIM; ++d)
    {
        Blas::Dgemm(
            'N', 'N', derivUtil->ptsHigh, DIM * nElmt, derivUtil->ptsLow, 1.0,
            derivUtil->VdmD[d].GetRawPtr(), derivUtil->ptsHigh, X,
            derivUtil->ptsLow, 0.0, &deriv[d][0][0][0], derivUtil->ptsHigh);
    }

    NekDouble integral = 0.0;
    NekDouble gam = numeric_limits<float>::epsilon();
    NekDouble ep = minJac < gam ? sqrt(gam*(gam-minJac)) : 0.0;
    NekDouble jacIdeal[DIM][DIM], jacDet;

    switch(opti)
    {
        case eLinEl:
        {
            const NekDouble nu = 0.45;
            const NekDouble mu = 1.0 / 2.0 / (1.0+nu);
            const NekDouble K  = 1.0 / 3.0 / (1.0 - 2.0 * nu);

            for (int i = 0; i < nElmt; ++i)
            {
                for(int k = 0; k < derivUtil->ptsHigh; ++k)
                {
                    jacDet = CalcIdealJac(i, k, deriv, data, jacIdeal);

                    NekDouble trEtE = LinElasTrace<DIM>(jacIdeal);
                    NekDouble sigma =
                        0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*ep*ep));
                    NekDouble lsigma = log(sigma);
                    integral += derivUtil->quadW[k] *
                                fabs(data[i]->maps[k][9]) *
                                (K * 0.5 * lsigma * lsigma + mu * trEtE);
                }
            }
            break;
        }

        case eHypEl:
        {
            const NekDouble nu = 0.4;
            const NekDouble mu = 1.0 / 2.0 / (1.0+nu);
            const NekDouble K  = 1.0 / 3.0 / (1.0 - 2.0 * nu);

            grad = Array<OneD, NekDouble>(DIM == 2 ? 5 : 9, 0.0);

            for (int i = 0; i < nElmt; ++i)
            {
                for(int k = 0; k < derivUtil->ptsHigh; ++k)
                {
                    jacDet = CalcIdealJac(i, k, deriv, data, jacIdeal);
                    NekDouble I1 = FrobeniusNorm(jacIdeal);

                    NekDouble sigma =
                        0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*ep*ep));
                    NekDouble lsigma = log(sigma);
                    integral += derivUtil->quadW[k]*
                        fabs(data[i]->maps[k][9]) *
                                (0.5 * mu * (I1 - 3.0 - 2.0*lsigma) +
                                 0.5 * K * lsigma * lsigma);

                    // Derivative of basis function in each direction
                    if(analytic)
                    {
                        NekDouble jacInvTrans[DIM][DIM];
                        NekDouble jacDetDeriv[DIM];

                        NekDouble phiM[DIM][DIM];
                        for (int m = 0; m < DIM; ++m)
                        {
                            for (int n = 0; n < DIM; ++n)
                            {
                                phiM[n][m] = deriv[m][i][n][k];
                            }
                        }

                        InvTrans<DIM>(phiM, jacInvTrans);
                        NekDouble derivDet = Determinant<DIM>(phiM);

                        NekDouble basisDeriv[DIM];
                        for (int m = 0; m < DIM; ++m)
                        {
                            basisDeriv[m] = *(derivUtil->VdmD[m])(k,nodeIds[i]);
                        }

                        for (int m = 0; m < DIM; ++m)
                        {
                            jacDetDeriv[m] = 0.0;
                            for (int n = 0; n < DIM; ++n)
                            {
                                jacDetDeriv[m] += jacInvTrans[m][n] * basisDeriv[n];
                            }
                            jacDetDeriv[m] *= derivDet / fabs(data[i]->maps[k][9]);
                        }

                        NekDouble jacDeriv[DIM][DIM][DIM];
                        for (int m = 0; m < DIM; ++m)
                        {
                            for (int n = 0; n < DIM; ++n)
                            {
                                NekDouble delta = m == n ? 1.0 : 0.0;
                                for (int l = 0; l < DIM; ++l)
                                {
                                    jacDeriv[m][n][l] = delta * basisDeriv[l];
                                }
                            }
                        }

                        NekDouble jacDerivPhi[DIM][DIM][DIM];
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
                                            jacDeriv[p][m][l] * data[i]->maps[k][l + 3*n];
                                    }
                                }
                            }
                        }

                        NekDouble frobProd[DIM];
                        for (int m = 0; m < DIM; ++m)
                        {
                            frobProd[m] = FrobProd<DIM>(jacIdeal,jacDerivPhi[m]);
                        }

                        NekDouble frobProdHes[DIM][DIM]; //holder for the hessian frobprods
                        for (int m = 0; m < DIM; ++m)
                        {
                            for(int l = m; l < DIM; ++l)
                            {
                                frobProdHes[m][l] = FrobProd<DIM>(jacDerivPhi[m],jacDerivPhi[l]);
                            }
                        }

                        for (int j = 0; j < DIM; ++j)
                        {
                            grad[j] += derivUtil->quadW[k] * fabs(data[i]->maps[k][9]) * (
                                mu * frobProd[j] + (jacDetDeriv[j] / (2.0*sigma - jacDet)
                                                    * (K * lsigma - mu)));
                        }
                        int ct = 0;
                        for (int m = 0; m < DIM; ++m)
                        {
                            for(int l = m; l < DIM; ++l, ct++)
                            {
                                grad[ct+DIM] += derivUtil->quadW[k] * fabs(data[i]->maps[k][9]) * (
                                    mu * frobProdHes[m][l] + jacDetDeriv[m]*jacDetDeriv[l]*(
                                        K/(2.0*sigma-jacDet)/(2.0*sigma-jacDet) - jacDet*(K*lsigma-mu)));
                            }
                        }
                    }
                }
            }
            break;
        }

        case eRoca:
        {
            for (int i = 0; i < nElmt; ++i)
            {
                for(int k = 0; k < derivUtil->ptsHigh; ++k)
                {
                    jacDet = CalcIdealJac(i, k, deriv, data, jacIdeal);
                    NekDouble frob = FrobeniusNorm(jacIdeal);
                    NekDouble sigma = 0.5*(jacDet +
                                    sqrt(jacDet*jacDet + 4.0*ep*ep));
                    integral += derivUtil->quadW[k] *
                                fabs(data[i]->maps[k][9]) *
                        (frob / DIM / pow(fabs(sigma), 2.0/DIM) -1.0);
                }
            }
            break;
        }

        case eWins:
        {
            for (int i = 0; i < nElmt; ++i)
            {
                for(int k = 0; k < derivUtil->ptsHigh; ++k)
                {
                    jacDet = CalcIdealJac(i, k, deriv, data, jacIdeal);
                    NekDouble frob = FrobeniusNorm(jacIdeal);
                    NekDouble sigma = 0.5*(jacDet +
                                    sqrt(jacDet*jacDet + 4.0*ep*ep));
                    integral += derivUtil->quadW[k]*
                                fabs(data[i]->maps[k][9])*
                                (frob / sigma);
                }
            }
            break;
        }
    }

    return integral;
}

NodeOptiJob* NodeOpti::GetJob()
{
    return new NodeOptiJob(this);
}

}
}
