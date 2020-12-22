///////////////////////////////////////////////////////////////////////////////
//
// File: StdInterpBasis.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Demo for testing functionality of PhysEvaluateBasis
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Timer.h>
#include <iostream>

#include "StdDemoSupport.hpp"

int main(int argc, char *argv[])
{
    DemoSupport demo;
    demo.ParseArguments(argc, argv);
    StdExpansion *E = demo.CreateStdExpansion();

    if (E == nullptr)
    {
        return 1;
    }

    int nCoeffs   = E->GetNcoeffs();
    int dimension = E->GetShapeDimension();

    Array<OneD, Array<OneD, NekDouble>> coords = demo.GetCoords(E);

    int nPts = coords[0].size();
    Array<OneD, NekDouble>
        sol(nPts * nCoeffs),
        hold1(nPts),
        temp(nPts),
        temp1(nPts),
        temp2(nPts),
        out_eval(nPts * nCoeffs),
        out_eval1(nPts * nCoeffs),
        phys(nPts * nCoeffs),
        phys1(nCoeffs * nPts),
        phys2(nPts * nCoeffs),
        sol1(nPts * nCoeffs),
        sol2(nPts * nCoeffs);

    auto storage = E->GetPhysEvaluateStorage();

    Array<OneD, NekDouble> out_eval2 = E->PhysEvaluateBasis(coords, storage, sol, sol1, sol2);

    for (int k = 0; k < nCoeffs; ++k)
    {
        // Fill the 'solution' field with each of the modes using FillMode.
        E->FillMode(k, hold1);
        E->PhysDeriv(hold1, temp, temp1, temp2);

        switch (dimension) // It is meant to fall down
        {
            case 3:
                Vmath::Vcopy(nPts, &temp2[0], 1, &phys2[k * nPts], 1);
            case 2:
                Vmath::Vcopy(nPts, &temp1[0], 1, &phys1[k * nPts], 1);
            case 1:
                Vmath::Vcopy(nPts, &temp[0], 1, &phys[k * nPts], 1);
                Vmath::Vcopy(nPts, &hold1[0], 1, &out_eval1[k * nPts], 1);
            default:
                break;
        }
    }

    // Separate modes 0 to nCoeffs
    Array<OneD, NekDouble> tmp(nPts);
    Array<OneD, NekDouble> tmp2(nPts);
    NekDouble errL2 = 0, errLinf = 0;
    for (int i = 0; i < nCoeffs; i++)
    {
        Vmath::Vcopy(nPts, &out_eval1[0] + i * nPts, 1, &tmp[0], 1);
        Vmath::Vcopy(nPts, &out_eval2[0] + i * nPts, 1, &tmp2[0], 1);

        errL2 += E->L2(tmp, tmp2);
        errLinf += E->Linf(tmp, tmp2);

        Vmath::Vcopy(nPts, &phys[0] + i, nCoeffs, &tmp[0], 1);
        Vmath::Vcopy(nPts, &sol[0] + i, nCoeffs, &tmp2[0], 1);
        errL2 += E->L2(tmp, tmp2);
        errLinf += E->Linf(tmp, tmp2);

        Vmath::Vcopy(nPts, &phys1[0] + i, nCoeffs, &tmp[0], 1);
        Vmath::Vcopy(nPts, &sol1[0] + i, nCoeffs, &tmp2[0], 1);
        errL2 += E->L2(tmp, tmp2);
        errLinf += E->Linf(tmp, tmp2);

        Vmath::Vcopy(nPts, &phys2[0] + i, nCoeffs, &tmp[0], 1);
        Vmath::Vcopy(nPts, &sol2[0] + i, nCoeffs, &tmp2[0], 1);
        errL2 += E->L2(tmp, tmp2);
        errLinf += E->Linf(tmp, tmp2);
    }

    cout << "L infinity error : " << scientific << errLinf << endl;
    cout << "L 2 error        : " << scientific << errL2 << endl;

    return 0;
}
