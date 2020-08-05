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
Array<OneD, NekDouble> EvalPolyDerivx(Array<OneD, Array<OneD, NekDouble>> &pts)
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

Array<OneD, NekDouble> EvalPolyDerivy(Array<OneD, Array<OneD, NekDouble>> &pts)
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

Array<OneD, NekDouble> EvalPolyDerivz(Array<OneD, Array<OneD, NekDouble>> &pts)
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
    for (int i = 0; i < dimension; ++i)
    {
        ptypes[i] = "PolyEvenlySpaced";
    }
    StdExpansion *F = demo.CreateStdExpansion();

    const auto totPoints = (unsigned) F->GetTotPoints();


    Array<OneD, Array<OneD, NekDouble>> coordsE = demo.GetCoords(E);
    Array<OneD, Array<OneD, NekDouble>> coordsF = demo.GetCoords(F);
    Array<OneD, NekDouble> physIn(totPoints), physOut0(totPoints), physOut1(totPoints), physOut2(totPoints);
    Array<OneD, NekDouble> tmpIn(dimension), sol0(totPoints), sol1(totPoints), sol2(totPoints);


    Array<OneD, NekDouble> collcoor(coordsF.size());
    
    if(coordsF.size()>1)
    {
        //collapse coordsE:
        
        for(int i = 0; i<coordsE[0].size(); i++)
        {
        
            Array<OneD, NekDouble> tmp(coordsE.size());
            for(int j = 0; j<coordsF.size(); j++)
            {
                collcoor[j] = coordsE[j][i];
            }
            E->LocCoordToLocCollapsed(collcoor,tmp); 
            for(int j = 0; j<coordsE.size(); j++)
            {
                coordsE[j][i] = tmp[j];
            }
        }

        for(int i = 0; i<coordsF[0].size(); i++)
        {
             
            Array<OneD, NekDouble> tmp(coordsF.size());
            for(int j = 0; j<coordsF.size(); j++)
            {
                collcoor[j] = coordsF[j][i];
            }
            F->LocCoordToLocCollapsed(collcoor,tmp); 
            for(int j = 0; j<coordsF.size(); j++)
            {
                coordsF[j][i] = tmp[j];

            }
        }

    }

    if(dimension>2)
        sol2 = EvalPolyDerivz(coordsF);
    if(dimension>1)
        sol1 = EvalPolyDerivy(coordsF);
    if(dimension>0)
        sol0 = EvalPolyDerivx(coordsF);


    physIn = EvalPoly(coordsE);
    //cout<<"\n physIn.size = "<<physIn.size()<<" ";
    //for(int i = 0; i<physIn.size(); i++)
    //    cout<<" "<<physIn[i]<<" ";
    //cout<<"\n";

    // Evaluate polynomial at the set of elemental solution points.
    //cout<<"\n physout:\n";
    for (int i = 0; i < totPoints; ++i)
    {
        /*        for (int d = 0; d < dimension; ++d)
        {
            tmpIn[d] = coordsF[d][i];
            }*/
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
        
    }

        /*for(int i = 0; i<coordsF[0].size(); i++)
        {
             
            Array<OneD, NekDouble> tmp(coordsF.size());
            for(int j = 0; j<coordsF.size(); j++)
            {
                collcoor[j] = coordsF[j][i];
            }
            F->LocCoordToLocCollapsed(collcoor,tmp); 
            for(int j = 0; j<coordsF.size(); j++)
            {
                coordsF[j][i] = tmp[j];
            }
        }*/

    /*    cout<<"\n sol:  ";
    for(int i = 0; i<sol2.size(); i++)
        cout<<sol2[i]<<" ";
    cout<<"\nphysout = ";
    for(int i = 0; i<physOut2.size(); i++)
        cout<<physOut2[i]<<" ";
    cout<<"\n";

    
    cout<<"\n sol:  ";
    for(int i = 0; i<sol1.size(); i++)
        cout<<sol1[i]<<" ";
    cout<<"\nphysout = ";
    for(int i = 0; i<physOut1.size(); i++)
        cout<<physOut1[i]<<" ";

    
    cout<<"\n sol x:  ";
    for(int i = 0; i<sol0.size(); i++)
        cout<<sol0[i]<<" ";
    cout<<"\nphysout x = ";
    for(int i = 0; i<physOut0.size(); i++)
        cout<<physOut0[i]<<" ";


    cout<<"\n coordsE = ";
    for(int i = 0; i<coordsE[0].size(); i++)
        cout<<coordsE[0][i]<<" "<<coordsE[1][i]<<"\t";
    cout<<"\n coordsF = ";
    for(int i = 0; i<coordsF[0].size(); i++)
        cout<<coordsF[0][i]<<" "<<coordsF[1][i]<<"\t";

    */ 
    cout << "\n\nL infinity error: " << scientific << E->Linf(physOut0, sol0)  +E->Linf(physOut1, sol1) +E->Linf(physOut2, sol2) << endl;
  cout << "L 2 error         : " << scientific << E->L2  (physOut0, sol0)+E->L2  (physOut1, sol1)+E->L2  (physOut2, sol2)<< endl;
    

    return 0;
}

