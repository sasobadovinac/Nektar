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

namespace Nektar
{
namespace Utilities
{

//returns zero if hessian is good, 1 if indefinite
template<int DIM> int NodeOpti::IsIndefinite(Grad &grad)
{
    ASSERTL0(false,"DIM error");
}

template<> int NodeOpti::IsIndefinite<1>(Grad &grad)
{

}

template<> int NodeOpti::IsIndefinite<2>(Grad &grad)
{

}

template<> int NodeOpti::IsIndefinite<3>(Grad &grad)
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
    NekDouble dum;

    DNekMat eval   (nVel, nVel, 0.0, eDIAGONAL);
    Array<OneD, NekDouble> vl  (nVel*nVel);
    Array<OneD, NekDouble> work(worklen);
    Array<OneD, NekDouble> wi  (nVel);

    Lapack::Dgeev('N', 'N', nVel, H.GetRawPtr(), nVel,
                  &(eval.GetPtr())[0], &wi[0], &vl[0], nVel,
                  &dum, nVel,
                  &work[0], worklen, info);

    ASSERTL0(!info,"dgeev failed");

    if(eval(0,0) < 0.0 || eval(1,1) < 0.0 || eval(2,2) < 0.0)
    {
        if(eval(0,0) < 0.0 && eval(1,1) < 0.0 && eval(2,2))
        {
            return 2;
        }
        else
        {
            return 1;
        }
    }

    return 0;*/
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
