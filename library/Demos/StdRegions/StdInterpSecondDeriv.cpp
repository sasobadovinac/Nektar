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
    // polynomial = x^3 + y^3 - z^2 + xy -3yz + 2y^2z +z 
   
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] = pow(pts[0][i],3) 
            + (dim >= 2 ?  pow(pts[1][i],3) + pts[0][i]*pts[1][i] : 0.0)
            + (dim >= 3 ? -1.0* pow(pts[2][i], 2) - 3*pts[1][i]*pts[2][i]  + 2*pow(pts[1][i],2)*pts[2][i] + pts[2][i]: 0.0);
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
    // polynomial = x^3 + y^3 - z^2 + xy -3yz + 2y^2z +z 
    // derivative in x = 3x^2 + y
    // second der in x = 6x
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  6.0*pts[0][i];//*pow(pts[0][i],);
    }
    return ret;
}

Array<OneD, NekDouble> EvalPolyDerivy(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    unsigned dim = pts.size();
    // polynomial = x^4 + y^3 - z^2 + xy -3yz + 2y^2z +z 
    // derivative in y = 3y^2 + x - 3z + 4yz
    // second der in y = 6y + 4z
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  6*pts[1][i]             + (dim >= 3 ?  4*pts[2][i] : 0.0);

    }
    return ret;
}

Array<OneD, NekDouble> EvalPolyDerivz(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());

    // polynomial = x^4 + y^3 - z^2 + xy -3yz + 2y^2z +z 
    // derivative in z = -2z - 3y +2y^2 + 1
    // second der in z = -2 
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
    
    if(dimension>2)
    {
     
        ASSERTL0(false, "Not implemented for dim > 1 yet");       
    }
    
    else if(dimension>1)
    {
        ASSERTL0(false, "Not implemented for dim > 1 yet");

    }
    
    else if(dimension>0)
    {
        for(int ii = 0; ii < coordsF[0].size(); ii++)
            E->PhysEvalSecondDeriv(coordsF, physIn, physOut0, NullNekDouble1DArray, NullNekDouble1DArray);    
    }

    
    if(dimension>2)
    {
        ASSERTL0(false, "Not implemented for dim > 1 yet");
    }
    if(dimension>1)
    {
        ASSERTL0(false, "Not implemented for dim > 1 yet");
    }
    if(dimension>0)
    {
        sol0 = EvalPolyDerivx(coordsF);
    }// cout<<"\n sol0:\n\n";

    // for(int k = 0; k < sol0.size(); k++)
    // {
    //     cout<<sol0[k]<<" ";
    // }
    // cout<<"\n physout0:\n";

    // for(int k = 0; k < sol0.size(); k++)
    // {
    //     cout<<physOut0[k]<<" ";
    // }
    
    cout << "\nL infinity error: " << scientific << E->Linf(physOut0, sol0)+ E->Linf(physOut1, sol1)+E->Linf(physOut2, sol2)  << endl;
cout << "L 2 error         : " << scientific << E->L2  (physOut0, sol0)+ E->L2  (physOut1, sol1)+E->L2  (physOut2, sol2)<< endl;
    

    return 0;
}

