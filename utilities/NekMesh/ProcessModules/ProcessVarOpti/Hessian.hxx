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

#ifndef UTILITIES_NEKMESH_NODEOPTI_HESSIAN
#define UTILITIES_NEKMESH_NODEOPTI_HESSIAN

#define     PI   3.14159265358979323846


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



//returns zero if hessian is good, 1 if indefinite
template<int DIM> int ProcessVarOpti::IsIndefinite(const Grad &grad, const int node)
{
    ASSERTL0(false,"DIM error");
}

template<> int ProcessVarOpti::IsIndefinite<1>(const Grad &grad, const int node)
{

}

template<> int ProcessVarOpti::IsIndefinite<2>(const Grad &grad, const int node)
{

}

template<> int ProcessVarOpti::IsIndefinite<3>(const Grad &grad, const int node)
{
    double H[3][3];
    H[0][0] = grad.G(node,3);
    H[1][0] = grad.G(node,4);
    H[0][1] = H[1][0];
    H[2][0] = grad.G(node,5);
    H[0][2] = H[2][0];
    H[1][1] = grad.G(node,6);
    H[2][1] = grad.G(node,7);
    H[1][2] = H[2][1];
    H[2][2] = grad.G(node,8);

    double eval[3]; // the eigenvalues

    double p1 = H[0][1] * H[0][1] + H[0][2] * H[0][2] + H[1][2] * H[1][2];
    if (p1 == 0.0) 
    {  
        // H is diagonal
        eval[0] = H[0][0];
        eval[1] = H[1][1];
        eval[2] = H[2][2];
    }
    else
    {
        double q  = (H[0][0] + H[1][1] + H[2][2]) / 3.0;
        double p2 =    (H[1][1] - q)*(H[1][1] - q)
                     + (H[2][2] - q)*(H[2][2] - q)
                     + (H[3][3] - q)*(H[3][3] - q)
                     + 2.0 * p1;
        double p = sqrt(p2 / 6.0);

        double B[3][3];   // B = (1.0 / p) * (H - q * I)   with I being the identity matrix
        B[0][0] = (1.0 / p) * (H[0][0] - q);
        B[1][1] = (1.0 / p) * (H[1][1] - q);
        B[2][2] = (1.0 / p) * (H[2][2] - q);
        B[0][1] = 1.0/p * H[0][1];
        B[1][0] = B[0][1];
        B[0][2] = 1.0/p * H[0][2];
        B[2][0] = B[0][2];
        B[1][2] = 1.0/p * H[1][2];
        B[2][1] = B[1][2];

        double r = Determinant<3>(B) / 2.0;

        // In exact arithmetic for h symmetric matrix  -1 <= r <= 1
        // but computation error can leave it slightly outside this range.
        double phi;
        if (r <= -1)
        { 
            phi = PI / 3.0;
        }
        else if (r >= 1)
        {
            phi = 0.0;
        }
        else
        {
            phi = acos(r) / 3.0;
        }    

        // the eigenvalues satisfy eval[2] <= eval[1] <= eval[0]
        eval[0] = q + 2.0 * p * cos(phi);
        eval[2] = q + 2.0 * p * cos(phi + (2.0*PI/3.0));
        eval[1] = 3.0 * q - eval[0] - eval[2];     // since trace(H) = eval[0] + eval[1] + eval[2]
    }

    NekMatrix<NekDouble> HH(3,3);
    HH(0,0) = H[0][0];
    HH(1,0) = H[1][0];
    HH(0,1) = H[0][1];
    HH(2,0) = H[2][0];
    HH(0,2) = H[0][2];
    HH(1,1) = H[1][1];
    HH(2,1) = H[2][1];
    HH(1,2) = H[1][2];
    HH(2,2) = H[2][2];

    int nVel = 3;
    char jobvl = 'N', jobvr = 'N';
    int worklen = 8*nVel, info;
    NekDouble dum;

    DNekMat evalDia   (nVel, nVel, 0.0, eDIAGONAL);
    Array<OneD, NekDouble> vl  (nVel*nVel);
    Array<OneD, NekDouble> work(worklen);
    Array<OneD, NekDouble> wi  (nVel);

    Lapack::Dgeev(jobvl, jobvr, nVel, HH.GetRawPtr(), nVel,
                  &(evalDia.GetPtr())[0], &wi[0], &vl[0], nVel,
                  &dum, nVel,
                  &work[0], worklen, info);

    if ( abs(evalDia(0,0) - eval[0]) <= 1e-08 &&
         abs(evalDia(1,1) - eval[1]) <= 1e-08 &&
         abs(evalDia(2,2) - eval[2]) <= 1e-08 )
    {}
    else
    {
        printf("%s\n", "wrong eigenvalues");
    }

    if(eval[0] < 0.0 || eval[1] < 0.0 || eval[2] < 0.0)
    {
        if(eval[0] < 0.0 && eval[1] < 0.0 && eval[2])
        {
            return 2;
        }
        else
        {
            return 1;
        }
    }

    return 0;
}

template<int DIM> void NodeOpti::MinEigen(NekDouble &val, NekDouble (&vec)[DIM], Grad &grad)
{
    ASSERTL0(false,"DIM error");
}

template<> void NodeOpti::MinEigen<1>(NekDouble &val, NekDouble (&vec)[1], Grad &grad)
{
}

template<> void NodeOpti::MinEigen<2>(NekDouble &val, NekDouble (&vec)[2], Grad &grad)
{

}

template<> void NodeOpti::MinEigen<3>(NekDouble &val, NekDouble (&vec)[3], Grad &grad)
{
    /*NekMatrix<NekDouble> H(3,3);
    H(0,0) = grad.h_G[3];
    H(1,0) = grad.h_G[4];
    H(0,1) = H(1,0);
    H(2,0) = grad.h_G[5];
    H(0,2) = H(2,0);
    H(1,1) = grad.h_G[6];
    H(2,1) = grad.h_G[7];
    H(1,2) = H(2,1);
    H(2,2) = grad.h_G[8];

    int nVel = 3;
    int worklen = 8*nVel, info;

    DNekMat eval   (nVel, nVel, 0.0, eDIAGONAL);
    DNekMat evec   (nVel, nVel, 0.0, eFULL);
    Array<OneD, NekDouble> vl  (nVel*nVel);
    Array<OneD, NekDouble> work(worklen);
    Array<OneD, NekDouble> wi  (nVel);

    Lapack::Dgeev('N', 'V', nVel, H.GetRawPtr(), nVel,
                  &(eval.GetPtr())[0], &wi[0], &vl[0], nVel,
                  &(evec.GetPtr())[0], nVel,
                  &work[0], worklen, info);

    ASSERTL0(!info,"dgeev failed");

    int minI;
    NekDouble tmp = std::numeric_limits<double>::max();
    for(int i = 0; i < 3; i++)
    {
        if(eval(i,i) < tmp)
        {
            minI = i;
            tmp = eval(i,i);
        }
    }

    val = eval(minI,minI);
    vec[0] = evec(0,minI);
    vec[1] = evec(1,minI);
    vec[2] = evec(2,minI);*/

}

}
}

#endif
