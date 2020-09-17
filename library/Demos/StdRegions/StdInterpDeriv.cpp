///////////////////////////////////////////////////////////////////////////////
//
// File: StdInterp.cpp
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

// polynomial = x^2 + y^2 - z^2
Array<OneD, NekDouble> EvalPoly(Array<OneD, Array<OneD, NekDouble>> &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    unsigned dim = pts.size();
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] = pow(pts[0][i],2) 
            + (dim >= 2 ? pow(pts[1][i], 2) : 0.0)
            - (dim >= 3 ? pow(pts[2][i], 2) : 0.0);
    }

    return ret;
}

// derivative in x = 2x
Array<OneD, NekDouble> EvalPolyDerivx(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  2*pts[0][i];
    }

    return ret;
}

// derivative in y = 2y
Array<OneD, NekDouble> EvalPolyDerivy(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  2*pts[1][i];
    }

    return ret;
}

// derivative in z = -2z
Array<OneD, NekDouble> EvalPolyDerivz(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  -2*pts[2][i] ;
    }

    return ret;
}

int main(int argc, char *argv[])
{
    DemoSupport demo;
    demo.ParseArguments(argc, argv);
    StdExpansion *E = demo.CreateStdExpansion();

    const auto dimension = (unsigned) E->GetShapeDimension();

    // Create a new element but with the evenly-spaced points type, so that we
    // perform a PhysEvaluateDeriv at a different set of nodal points
    // (i.e. non-collocated interpolation).
    vector<string> &ptypes = demo.GetPointsType();
    ptypes[0] = "GaussRadauPLegendre";
    ptypes[1] = "GaussRadauMLegendre";
    ptypes[2] = "GaussRadauMLegendre";
    StdExpansion *F = demo.CreateStdExpansion();

    const auto totPoints = (unsigned) E->GetTotPoints();

    Array<OneD, Array<OneD, NekDouble> > coordsE = demo.GetCoords(E);
 
    int totpts = coordsE[0].size();
    Array<OneD, NekDouble> physIn(totPoints), physOut0(totpts), physOut1(totpts), physOut2(totpts);
    Array<OneD, NekDouble>  sol0(totpts), sol1(totpts), sol2(totpts);
    
    physIn = EvalPoly(coordsE);
   
    const Array<OneD, const Array<OneD, NekDouble> > coordsF = demo.GetCoords(F);

    E->PhysEvalGrad(coordsF, physIn, physOut0, physOut1, physOut2);

    switch(dimension)
    {
        case 3:
            sol2 = EvalPolyDerivz(coordsF);
        case 2:
            sol1 = EvalPolyDerivy(coordsF);
        case 1:
            sol0 = EvalPolyDerivx(coordsF);
        default:
            break;
    }

    cout << "\nL infinity error: " << scientific << E->Linf(physOut0, sol0) +
                                                    E->Linf(physOut1, sol1) +
                                                    E->Linf(physOut2, sol2)
                                                 << endl;
    cout << "L 2 error         : " << scientific << E->L2(physOut0, sol0) +
                                                    E->L2(physOut1, sol1) +
                                                    E->L2(physOut2, sol2)
                                                 << endl;

    return 0;
}

