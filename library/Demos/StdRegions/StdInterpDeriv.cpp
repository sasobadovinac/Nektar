///////////////////////////////////////////////////////////////////////////////
//
// File: StdInterpDeriv.cpp
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
// Description: Demo for testing functionality of PhysEvaluateDeriv
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <LibUtilities/BasicUtils/Timer.h>

#include "StdDemoSupport.hpp"

// polynomial = 2*x^4 + 3*y^2 - z^3
Array<OneD, NekDouble> EvalPoly(const Array<OneD, const Array<OneD, NekDouble>> &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    unsigned dim = pts.size();
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] = 2 * (pts[0][i] * pts[0][i] * pts[0][i] * pts[0][i])
            + (dim >= 2 ? 3 * (pts[1][i] * pts[1][i]) : 0.0)
            - (dim >= 3 ? (pts[2][i] * pts[2][i] * pts[2][i]) : 0.0);
    }

    return ret;
}

// derivative in x = 8x^3
Array<OneD, NekDouble> EvalPolyDerivx(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  8*(pts[0][i] * pts[0][i] * pts[0][i]);
    }

    return ret;
}

// derivative in y = 6y
Array<OneD, NekDouble> EvalPolyDerivy(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  6*pts[1][i];
    }

    return ret;
}

// derivative in z = -3*z^2
Array<OneD, NekDouble> EvalPolyDerivz(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  -3*(pts[2][i] * pts[2][i]);
    }

    return ret;
}

int main(int argc, char *argv[])
{
    DemoSupport demo;
    demo.ParseArguments(argc, argv);
    StdExpansion *E = demo.CreateStdExpansion();

    const auto dimension = (unsigned) E->GetShapeDimension();

    const auto totPoints = (unsigned)E->GetTotPoints();
    Array<OneD, NekDouble> physIn(totPoints, 0.0), physOut(totPoints, 0.0),
        physOut0(totPoints, 0.0), physOut1(totPoints, 0.0),
        physOut2(totPoints, 0.0), sol(totPoints, 0.0), sol0(totPoints, 0.0),
        sol1(totPoints, 0.0), sol2(totPoints, 0.0);

    Array<OneD, Array<OneD, NekDouble> > coordsE = demo.GetCoords(E);
    physIn = EvalPoly(coordsE);

    // Create a new element but with PolyEvenlySpaced points, so that we
    // perform a PhysEvaluateDeriv at a different set of nodal points
    // (i.e. non-collocated interpolation), all tests use default types
    vector<string> &ptypes = demo.GetPointsType();
    for (int i = 0; i < dimension; ++i)
    {
        ptypes[i] = "PolyEvenlySpaced";
    }
    StdExpansion *F = demo.CreateStdExpansion();
    const Array<OneD, const Array<OneD, NekDouble> > coordsF = demo.GetCoords(F);

    for (int i = 0; i < totPoints; ++i)
    {
        // Fill coords array
        Array<OneD, NekDouble> coordIn(dimension, 0.0);
        for (int j = 0; j < dimension; ++j)
        {
            coordIn[j] = coordsF[j][i];
        }

        // Operating on expansion E so using physIn from E, and coordIn from F
        physOut[i] = E->PhysEvaluate(coordIn, physIn, physOut0[i], physOut1[i], physOut2[i]);
    }

    switch(dimension)
    {
        case 3:
            sol2 = EvalPolyDerivz(coordsF);
            /* fall through */
        case 2:
            sol1 = EvalPolyDerivy(coordsF);
            /* fall through */
        case 1:
            sol0 = EvalPolyDerivx(coordsF);
            /* fall through */
        default:
            sol = EvalPoly(coordsF);
            break;
    }

    cout << "\nL infinity error: " << scientific << F->Linf(physOut, sol) +
                                                    F->Linf(physOut0, sol0) +
                                                    F->Linf(physOut1, sol1) +
                                                    F->Linf(physOut2, sol2)
                                                  << endl;
    cout << "L 2 error         : " << scientific << F->L2(physOut, sol) +
                                                    F->L2(physOut0, sol0) +
                                                    F->L2(physOut1, sol1) +
                                                    F->L2(physOut2, sol2)
                                                 << endl;

    return 0;
}
