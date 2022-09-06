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

#include <LibUtilities/BasicUtils/Timer.h>
#include <iostream>

#include "StdDemoSupport.hpp"

// polynomial = 4x^4 + 3y^3 + 2z^2 + 1
Array<OneD, NekDouble> EvalPoly(
    const Array<OneD, const Array<OneD, NekDouble>> &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    unsigned dim = pts.size();
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] = (4 * pts[0][i] * pts[0][i] * pts[0][i] * pts[0][i]) +
                 (dim >= 2 ? (3 * pts[1][i] * pts[1][i] * pts[1][i]) : 0.0) +
                 (dim >= 3 ? (2 * pts[2][i] * pts[2][i]) : 0.0) + 1;
    }

    return ret;
}

// derivative in x = 16x^3
Array<OneD, NekDouble> EvalPolyDerivx(
    const Array<OneD, const Array<OneD, NekDouble>> &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] = 16 * pts[0][i] * pts[0][i] * pts[0][i];
    }

    return ret;
}

// derivative in y = 9y^2
Array<OneD, NekDouble> EvalPolyDerivy(
    const Array<OneD, const Array<OneD, NekDouble>> &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] = 9 * pts[1][i] * pts[1][i];
    }

    return ret;
}

// derivative in z = 4z
Array<OneD, NekDouble> EvalPolyDerivz(
    const Array<OneD, const Array<OneD, NekDouble>> &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] = 4 * pts[2][i];
    }

    return ret;
}

int main(int argc, char *argv[])
{
    DemoSupport demo;
    demo.ParseArguments(argc, argv);
    StdExpansion *E = demo.CreateStdExpansion();

    const auto dimension = (unsigned)E->GetShapeDimension();

    const auto totPoints = (unsigned)E->GetTotPoints();
    Array<OneD, NekDouble> physIn(totPoints, 0.0), physOut(totPoints, 0.0),
        sol(totPoints, 0.0), sol0(totPoints, 0.0), sol1(totPoints, 0.0),
        sol2(totPoints, 0.0);

    Array<OneD, std::array<NekDouble, 3>> physOutDeriv(totPoints);

    Array<OneD, Array<OneD, NekDouble>> coordsE = demo.GetCoords(E);
    physIn                                      = EvalPoly(coordsE);

    // Create a new element but with GaussGaussChebyshev points so that we
    // perform a PhysEvaluateDeriv at a different set of nodal points
    // (i.e. non-collocated interpolation), all tests use default types
    // initially
    vector<string> &ptypes = demo.GetPointsType();
    for (int i = 0; i < dimension; ++i)
    {
        ptypes[i] = "PolyEvenlySpaced";
    }

    StdExpansion *F = demo.CreateStdExpansion();
    const Array<OneD, const Array<OneD, NekDouble>> coordsF = demo.GetCoords(F);

    for (int i = 0; i < totPoints; ++i)
    {
        // Fill coords array
        Array<OneD, NekDouble> coordIn(dimension, 0.0);
        for (int j = 0; j < dimension; ++j)
        {
            coordIn[j] = coordsF[j][i];
        }

        // Operating on expansion E so using physIn from E, and coordIn from F
        physOut[i] = E->PhysEvaluate(coordIn, physIn, physOutDeriv[i]);
    }

    // Extract derivatives in to array for error calculation
    Array<OneD, NekDouble> physOut0(totPoints, 0.0), physOut1(totPoints, 0.0),
        physOut2(totPoints, 0.0);
    switch (dimension)
    {
        case 3:
            for (int j = 0; j < totPoints; ++j)
            {
                physOut2[j] = physOutDeriv[j][2];
            }
            sol2 = EvalPolyDerivz(coordsF);
            /* fall through */
        case 2:
            for (int j = 0; j < totPoints; ++j)
            {
                physOut1[j] = physOutDeriv[j][1];
            }
            sol1 = EvalPolyDerivy(coordsF);
            /* fall through */
        case 1:
            for (int j = 0; j < totPoints; ++j)
            {
                physOut0[j] = physOutDeriv[j][0];
            }
            sol0 = EvalPolyDerivx(coordsF);
            /* fall through */
        default:
            sol = EvalPoly(coordsF);
            break;
    }

    cout << "\nL infinity error: " << scientific
         << E->Linf(physOut, sol) + E->Linf(physOut0, sol0) +
                E->Linf(physOut1, sol1) + E->Linf(physOut2, sol2)
         << endl;
    cout << "L 2 error         : " << scientific
         << E->L2(physOut, sol) + E->L2(physOut0, sol0) +
                E->L2(physOut1, sol1) + E->L2(physOut2, sol2)
         << endl;

    return 0;
}
