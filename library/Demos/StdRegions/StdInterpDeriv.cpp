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

// Evaluate polynomial for testing and save in ret (size same as pts[0]) if
Array<OneD, NekDouble> EvalPoly(Array<OneD, Array<OneD, NekDouble>> &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    unsigned dim = pts.size();

    // check if pts[0] and pts[1] have same size
    // polynomial = x^2 + y^2 - z^2 
   
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] = pow(pts[0][i],2) 
            + (dim >= 2 ? pow(pts[1][i], 2) : 0.0)
            - (dim >= 3 ? pow(pts[2][i], 2) : 0.0);
    }
    return ret;
}



// Evaluate polynomial for testing and save in ret (size same as pts[0]) if
// tensorp = 0, we need tensorprod else just eval at pts
Array<OneD, NekDouble> EvalPolyDerivx(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    //    unsigned dim = pts.size();

    // check if pts[0] and pts[1] have same size
    // polynomial = x^2 + y^2 - z^2 
    // derivative in x = 2x 
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  2*pts[0][i];
    }
    return ret;
}

Array<OneD, NekDouble> EvalPolyDerivy(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());

    // check if pts[0] and pts[1] have same size
    // polynomial = x^2 + y^2 - z^2
    // derivative in y = 2y

    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  2*pts[1][i];
    }
    return ret;
}

Array<OneD, NekDouble> EvalPolyDerivz(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());

    // check if pts[0] and pts[1] have same size
    // polynomial = x^2 + y^2 - z^2 
    // derivative in z = -2z 
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
   
    // Doubt: should it work for polyevenlyspaced?
    // Only with these ptypes, this demo works for triangles
    ptypes[0] = "GaussRadauPLegendre"; // if ptypes[0] is MLegendre, hell breaks lose!
    
    ptypes[1] = "GaussRadauMLegendre"; // if ptypes[1] is PLegendre, another hell breaks lose!
    
    ptypes[2] = "GaussRadauMLegendre";
    
    
    StdExpansion *F = demo.CreateStdExpansion();

    const auto totPoints = (unsigned) E->GetTotPoints();


    Array<OneD, Array<OneD, NekDouble> > coordsE = demo.GetCoords(E);
 
    int totpts = coordsE[0].size();
    Array<OneD, NekDouble> physIn(totPoints), physOut0(totpts), physOut1(totpts), physOut2(totpts);
    Array<OneD, NekDouble>  sol0(totpts), sol1(totpts), sol2(totpts);
    
    physIn = EvalPoly(coordsE);
   
    const Array<OneD, const Array<OneD, NekDouble> > coordsF = demo.GetCoords(F);
    
    if(dimension>2)
    {

        E->PhysEvalGrad(coordsF, physIn, physOut0, physOut1, physOut2);

    }
    
    else if(dimension>1)
    {
        
        E->PhysEvalGrad(coordsF, physIn, physOut0, physOut1, NullNekDouble1DArray);

    }
    
    else if(dimension>0)
    {
        E->PhysEvalGrad(coordsF, physIn, physOut0, NullNekDouble1DArray, NullNekDouble1DArray);    
    }
 
    if(dimension>2)
        sol2 = EvalPolyDerivz(coordsF);
    if(dimension>1)
        sol1 = EvalPolyDerivy(coordsF);
    if(dimension>0)
        sol0 = EvalPolyDerivx(coordsF);
    

    
    cout << "\nL infinity error: " << scientific << E->Linf(physOut2, sol2) +E->Linf(physOut1, sol1) + E->Linf(physOut0, sol0)  << endl;
    cout << "L 2 error         : " << scientific << E->L2  (physOut2, sol2) + E->L2  (physOut1, sol1) + E->L2  (physOut0, sol0)<< endl;
    

    return 0;
}

