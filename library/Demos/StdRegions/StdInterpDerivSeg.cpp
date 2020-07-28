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
        ret[i] = pow(pts[0][i],4) 
            + (dim >= 2 ? pow(pts[1][i], 2) : 0.0)
            + (dim >= 3 ? pow(pts[2][i], 2) : 0.0);
    }
    return ret;
}


// Evaluate polynomial for testing and save in ret (size same as pts[0]) if
// tensorp = 0, we need tensorprod else just eval at pts
Array<OneD, NekDouble> EvalPolyDerivx(Array<OneD, Array<OneD, NekDouble>> &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    //    unsigned dim = pts.size();

    // check if pts[0] and pts[1] have same size
    // polynomial = x^2 + y^2 - z^2 
    // derivative in x = 2x 
    unsigned dim = pts.size();

    

    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] = 4.0*pow(pts[0][i],3) 
            + (dim >= 2 ? pow(pts[1][i], 2) : 0.0)
            + (dim >= 3 ? pow(pts[2][i], 2) : 0.0);
    }
    return ret;
}



int main(int argc, char *argv[])
{
    DemoSupport demo;
    demo.ParseArguments(argc, argv);
    StdExpansion *E = demo.CreateStdExpansion();
    
    //    const auto totPoints = (unsigned) E->GetTotPoints();
                //    NekDouble totPoints = 1;
    const auto dimension = (unsigned) E->GetShapeDimension();

    // Create a new element but with the evenly-spaced points type, so that we
    // perform a PhysEvaluateDeriv at a different set of nodal points
    // (i.e. non-collocated interpolation).
    vector<string> &ptypes = demo.GetPointsType();
    for (int i = 0; i < dimension; ++i)
    {
        ptypes[i] = "NodalTriFekete" ;
    }
    
    LibUtilities::BasisType btype;
    btype = LibUtilities::eOrtho_A;
    LibUtilities::PointsType quadPointsType =
	LibUtilities::eGaussLobattoLegendre;
    LibUtilities::PointsKey quadpkCheb(10, quadPointsType);
    LibUtilities::BasisKey bkey0(
                                 btype,
                                 
                                 4,quadpkCheb);
    StdRegions::StdSegExp tmp = StdSegExp(bkey0);
      
    StdRegions::StdSegExp *F = &tmp;
    //StdExpansion *F = demo.CreateStdExpansion();
    Array<OneD, Array<OneD, NekDouble>> coordsE = demo.GetCoords(E);
    Array<OneD, Array<OneD, NekDouble>> coordsF(2);
    coordsF[0] = Array<OneD,NekDouble>(F->GetTotPoints());
    coordsF[1] = Array<OneD,NekDouble>(F->GetTotPoints());
    F->GetCoords(coordsF[0],coordsF[1]);

    int totevalpts = F->GetTotPoints();

    Array<OneD, NekDouble> physIn(totevalpts), physOut0(totevalpts);
    Array<OneD, NekDouble> tmpIn(dimension), sol0(totevalpts);
    Array<OneD, NekDouble> collcoor(2);
     //Make sure that there are no  overlapping elements in coordsF and coordsE 
    for(int i = 0; i<coordsF[0].size(); i++)
    {
        for(int j = 0; j<coordsE[0].size(); j++)
        {
            if(coordsE[0][j] == coordsF[0][i])
            {
                cout<<"\n here\n";
                coordsF[0][i] +=0.001;
            }
        }

    }
    physIn = EvalPoly(coordsE);
    
    Array<OneD, NekDouble> tmp2(2);
    // Evaluate polynomial at the set of elemental solution points.
    //cout<<"\n physout:\n";
    for (int i = 0; i < sol0.size(); ++i)
    {
        for (int d = 0; d < dimension; ++d)
        {
            tmpIn[d] = coordsF[d][i];
            
        }
        collcoor = tmpIn;

        if(dimension>0)
        {
      
            physOut0[i] = E->PhysEvaluatedx(collcoor, physIn);
        }
        
    }    
    sol0 = EvalPolyDerivx(coordsF);
    

    cout<<"\n sol:  ";
    for(int i = 0; i<sol0.size(); i++)
        cout<<sol0[i]<<" ";
    cout<<"\nphysout = ";
    for(int i = 0; i<physOut0.size(); i++)
        cout<<physOut0[i]<<" ";
    cout<<"\n";

    
    cout << "\n\nL infinity error : " << scientific << F->Linf  (physOut0, sol0);
    cout << "L 2 error        : " << scientific << F->L2  (physOut0, sol0);
    

    return 0;
}

