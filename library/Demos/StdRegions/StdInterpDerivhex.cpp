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
// tensorp = 0, we need tensorprod else just eval at pts
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
            
        ret[i] =  2.0*pts[0][i];
    }
    return ret;
}

Array<OneD, NekDouble> EvalPolyDerivy(Array<OneD, Array<OneD, NekDouble>> &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    //unsigned dim = pts.size();

    // check if pts[0] and pts[1] have same size
    // polynomial = x^2 + y^2 - z^2
    // derivative in y = 2y

    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  2.0*pts[1][i];
    }
    return ret;
}

Array<OneD, NekDouble> EvalPolyDerivz(Array<OneD, Array<OneD, NekDouble>> &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    //unsigned dim = pts.size();

    // check if pts[0] and pts[1] have same size
    // polynomial = x^2 + y^2 - z^2
    // derivative in z = -2z9

    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  -2*pts[2][i];
    }
    return ret;
}

int main(int argc, char *argv[])
{
    DemoSupport demo;
    demo.ParseArguments(argc, argv);
    StdExpansion *E = demo.CreateStdExpansion();
    
    const auto totPoints = (unsigned) E->GetTotPoints();
    const auto dimension = (unsigned) E->GetShapeDimension();

    // Create a new element but with the evenly-spaced points type, so that we
    // perform a PhysEvaluate at a different set of nodal points
    // (i.e. non-collocated interpolation).
    vector<string> &ptypes = demo.GetPointsType();
    for (int i = 0; i < dimension; ++i)
    {
        ptypes[i] = "PolyEvenlySpaced";
    }


    LibUtilities::BasisType btype0;
    btype0 = LibUtilities::eOrtho_A;

    LibUtilities::BasisType btype1;
    btype1 = LibUtilities::eOrtho_A;

    LibUtilities::BasisType btype2;
    btype2 = LibUtilities::eOrtho_B;
    //    LibUtilities::BasisType btype2;
    //btype2 = LibUtilities::eOrtho_A;
    LibUtilities::PointsType quadPointsType =
	LibUtilities::eGaussLobattoLegendre;

    

    LibUtilities::PointsKey quadpkCheb(6, quadPointsType);
    LibUtilities::BasisKey bkey0(
                                 btype0,
                                 4,
                                 quadpkCheb);
    LibUtilities::BasisKey bkey1(
                                 btype1,
                                 4,
                                 quadpkCheb);
    LibUtilities::BasisKey bkey2(
                                   btype2,
                                   4,
                                   quadpkCheb);
    StdRegions::StdPrismExp *F = new StdPrismExp(bkey0,bkey1,bkey2);
      

    Array<OneD, Array<OneD, NekDouble>> coordsE = demo.GetCoords(E);
    Array<OneD, Array<OneD, NekDouble>> coordsF(3);
    coordsF[0] = Array<OneD,NekDouble>(F->GetTotPoints());
    coordsF[1] = Array<OneD,NekDouble>(F->GetTotPoints());
    coordsF[2] = Array<OneD,NekDouble>(F->GetTotPoints());

    F->GetCoords(coordsF[0],coordsF[1],coordsF[2]);

    int totevalpts = F->GetTotPoints();

    Array<OneD, NekDouble> physIn(totevalpts), physOut0(totevalpts), physOut1(totevalpts), physOut2(totPoints);
    Array<OneD, NekDouble> tmpIn(dimension), sol0(totevalpts), sol1(totevalpts), sol2(totPoints);

    Array<OneD, NekDouble> collcoor(3);
  
    //collapse coordsE and coordsF:
      for(int i = 0; i<coordsF[0].size(); i++)
    {
        Array<OneD, NekDouble> tmp(3);
        collcoor[0] = coordsF[0][i];
        collcoor[1] = coordsF[1][i];
        collcoor[2] = coordsF[2][i];
        F->LocCoordToLocCollapsed(collcoor,tmp); 
        coordsF[0][i] = tmp[0];
        coordsF[1][i] = tmp[1];
        coordsF[2][i] = tmp[2];

    }
      for(int i = 0; i<coordsE[0].size(); i++)
    {
        Array<OneD, NekDouble> tmp(3);
        collcoor[0] = coordsE[0][i];
        collcoor[1] = coordsE[1][i];
        collcoor[2] = coordsE[2][i];
        E->LocCoordToLocCollapsed(collcoor,tmp); 
        coordsE[0][i] = tmp[0];
        coordsE[1][i] = tmp[1];
        coordsE[2][i] = tmp[2];

        }

  

    //Make sure that there are no  overlapping elements in coordsF and coordsE 
    /*for(int i = 0; i<coordsF[0].size(); i++)
    {
        for(int j = 0; j<coordsE[0].size(); j++)
        {
            if(abs(coordsE[0][j]) == abs(coordsF[0][i]) )
            {
                coordsF[0][i] -=(rand()%100)*0.004;
            }
            if( abs(coordsE[1][j]) == abs(coordsF[1][i]) )
            {
            coordsF[1][i] -=(rand()%100)*0.002;
            }
            if(abs(coordsE[2][j]) == abs(coordsF[2][i]) )
            {
                coordsF[2][i] -=(rand()%100)*0.003;
            }
                                if(abs(abs(coordsE[0][j]) -abs( coordsF[0][i]))<1e-9)
        {
            cout<<"\n impossibru!!! F[0]["<<i<<"]="<<coordsF[0][i];
            cout<<"\n here coordsE[0][j]="<<coordsE[0][j]<<" coordsF[0][i] ="<<coordsF[0][i] <<" coordsE[1][j] = "<<coordsE[1][j] <<" coordsF[1][i]="<<coordsF[1][i]<<"\n";
            
            } 

        }
        
        }*/


    //Make sure that there are no  overlapping elements in coordsF and coordsE 
    for(int i = 0; i<coordsF[0].size(); i++)
    {
        for(int j = 0; j<coordsE[0].size(); j++)
        {
            if(abs(coordsE[0][j]) == abs(coordsF[0][i]) || abs(coordsE[1][j]) == abs(coordsF[1][i]) || abs(coordsE[2][j]) == abs(coordsF[2][i]))
            {
                coordsF[0][i] -=(rand()%100)*0.001;
                coordsF[1][i] -=(rand()%100)*0.002;
             coordsF[2][i] -=(rand()%100)*0.003;
            }
            /*            if(abs(abs(coordsE[0][j]) -abs( coordsF[0][i]))<1e-9)
            {
                cout<<"\n impossibru!!! F[0]["<<i<<"]="<<coordsF[0][i];
                cout<<"\n here coordsE[0][j]="<<coordsE[0][j]<<" coordsF[0][i] ="<<coordsF[0][i] <<" coordsE[1][j] = "<<coordsE[1][j] <<" coordsF[1][i]="<<coordsF[1][i]<<"\n";

                } 
            if(abs(coordsF[0][i])==1 && abs(coordsE[0][j])==1)
            {
                cout<<"\n coordsE[0][j]="<<coordsE[0][j]<<" coordsF[0][i]="<<coordsF[0][i]<<" ";
                }*/

        }

    }


    
 physIn = EvalPoly(coordsE);
    Array<OneD, NekDouble> tmp2(3);
    // Evaluate polynomial at the set of elemental solution points.

    for (int i = 0; i < sol0.size(); ++i)
    {
        for (int d = 0; d < dimension; ++d)
        {
            tmpIn[d] = coordsF[d][i];
            
        }
        collcoor = tmpIn;
        if(dimension>2)
        {
            physOut2[i] = E->PhysEvaluatedz(collcoor, physIn);
        }
        if(dimension>1)
        {
      
            physOut1[i] = E->PhysEvaluatedy(collcoor, physIn);
        }

        if(dimension>0)
        {
      
            physOut0[i] = E->PhysEvaluatedx(collcoor, physIn);
        }
        
    }    

    
    if(dimension>2)
        sol2 = EvalPolyDerivz(coordsF);
    if(dimension>1)
        sol1 = EvalPolyDerivy(coordsF);
    if(dimension>0)
        sol0 = EvalPolyDerivx(coordsF);
    
    /*  cout<<"\n\n\ndghajsdgjhasgdjyadgyagd2222222\n\n\n";
    
     */
      cout<<"\n sol:  ";
    for(int i = 0; i<sol0.size(); i++)
        cout<<sol0[i]<<" ";
    cout<<"\nphysout = ";
    for(int i = 0; i<physOut0.size(); i++)
        cout<<physOut0[i]<<" ";
    cout<<"\n";

      cout<<"\n sol:  ";
    for(int i = 0; i<sol1.size(); i++)
        cout<<sol1[i]<<" ";
    cout<<"\nphysout = ";
    for(int i = 0; i<physOut1.size(); i++)
        cout<<physOut1[i]<<" ";
    cout<<"\n";


    cout<<"\n sol:  ";
    for(int i = 0; i<sol2.size(); i++)
        cout<<sol2[i]<<" ";
    cout<<"\nphysout = ";
    for(int i = 0; i<physOut2.size(); i++)
        cout<<physOut2[i]<<" ";

    cout << "\n\nL infinity error : " << scientific << F->Linf  (physOut1, sol1)+F->Linf  (physOut0, sol0)+ F->Linf  (physOut2, sol2)   << endl;
    cout << "L 2 error        : " << scientific << F->L2  (physOut1, sol1)+F->L2  (physOut0, sol0)    +F->L2  (physOut2, sol2)    << endl;
    F = NULL;;
    E = NULL;
    
    return 0;
}

