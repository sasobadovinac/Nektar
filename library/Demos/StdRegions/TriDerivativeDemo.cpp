///////////////////////////////////////////////////////////////////////////////
//
// File: NodalDemo.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// ermission is hereby granted, free of charge, to any person obtaining a
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
// Description: Demo for testing functionality of StdProject
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "StdDemoSupport.hpp"
#include <fstream>
#include <time.h>
#include <LibUtilities/BasicUtils/Timer.h>

#define LOGNAME_FORMAT "%Y%m%d_%H%M%S"
#define LOGNAME_SIZE 20


namespace po = boost::program_options;

NekDouble Shape_sol(NekDouble x, NekDouble y, NekDouble z, vector<int> order,
                    vector<BasisType> btype, ShapeType stype, bool diff);

// //Modification to deal with exact solution for diff. Return 1 if integer < 0.
// static double pow_loc(const double val, const int i)
// {
//     return (i < 0) ? 1.0 : pow(val, i);
// }

Array<OneD, Array<OneD, NekDouble> > storage2d;

int iterstaken;

string funcdef;
// Colleague matrix
int numedges, numsurfaces;
StdExpansion *E;

//StdExpansion *E3quad;
Array<OneD, NekDouble> qZin; //inner point grid for 2D rootfinding
Array<OneD, NekDouble> qZinmid; //inner mid point grid for 2D rootfinding

// for tets edges root finding

// edge front left (AD) (x = -1) (y = -1)
Array<OneD, NekDouble> Vxm1ym1z   ;
Array<OneD, NekDouble> Vdyxm1ym1z   ;
Array<OneD, NekDouble> Vdxxm1ym1z   ;
Array<OneD, NekDouble> Vdzxm1ym1z   ;
               
//edge front hypt (DB) (y = -1) (z = -x)
Array<OneD, NekDouble> Vym1xmz   ;
Array<OneD, NekDouble> Vdxym1xmz   ;
Array<OneD, NekDouble> Vdyym1xmz   ;
Array<OneD, NekDouble> Vdzym1xmz   ;
               
//edge front bot (AB) (y = -1) (z = -1)
Array<OneD, NekDouble> Vym1xzm1   ;
Array<OneD, NekDouble> Vdyym1xzm1   ;
Array<OneD, NekDouble> Vdxym1xzm1   ;
Array<OneD, NekDouble> Vdzym1xzm1   ;

               
//edge left hypt (DC) ( x = -1) (z = -y)
Array<OneD, NekDouble> Vxm1ymz   ;
Array<OneD, NekDouble> Vdxxm1ymz   ;
Array<OneD, NekDouble> Vdyxm1ymz   ;
Array<OneD, NekDouble> Vdzxm1ymz   ;
               
// edge bot diag (BC) (z = -1) (y = -x)
Array<OneD, NekDouble> Vxmyzm1   ;
Array<OneD, NekDouble> Vdxxmyzm1   ;
Array<OneD, NekDouble> Vdyxmyzm1   ;
Array<OneD, NekDouble> Vdzxmyzm1   ;

//edge CA bot left (x = -1) (z = -1)
Array<OneD, NekDouble> Vxm1yzm1   ;
Array<OneD, NekDouble> Vdxxm1yzm1   ;
Array<OneD, NekDouble> Vdyxm1yzm1   ;
Array<OneD, NekDouble> Vdzxm1yzm1   ;


//surface bot z = -1, x + y = 0 (ABC)
Array<OneD, NekDouble>Vxmyzm1surf ;


//surface left x = -1 y + z = 0 (DAC)
Array<OneD, NekDouble>Vxm1yz ;

//surf front y = -1, x + z = 0 (DAB)
Array<OneD, NekDouble> Vxmzym1;

//surf DCB (x + y + z = -1)
Array<OneD, NekDouble> Vxyzm1;


int dimension ;
            
DemoSupport demo;

int main(int argc, char *argv[])
{
    demo.GetOptions().add_options()("optm,z", "positivity preserving optimizer");

    demo.ParseArguments(argc, argv);
    
    po::variables_map vm = demo.GetVariableMap();

       
    E = demo.CreateStdExpansion();
    if (E == nullptr)
    {
        return 1;
    }
    dimension = E->GetShapeDimension();
    // if(dimension < 3 )
    // {
    //     cout<<"\n dimension should be 3, try using StdProjectPositivityPres1D or StdProjectPositivityPres2D for other dimensions\n\n";
    //     exit(0);
    // }

    std::vector<int> order;
    std::vector<BasisType> btype(3, eNoBasisType);
    //    LibUtilities::PointsType pointsTypetriA = LibUtilities::eGaussLobattoLegendre;    

    for (int i = 0; i < dimension; ++i)
    {
        btype[i] = E->GetBasisType(i);
        order.push_back(E->GetBasisNumModes(i));
    }

    switch(E->DetShapeType())
    {
    case LibUtilities::eSegment:
        
        numedges = 1;
        numsurfaces = 0;
        break;
    

    case LibUtilities::eTriangle:
        
        numedges = 3;
        numsurfaces = 0;
        break;

    case LibUtilities::ePrism:
        
        numedges = 9;
        numsurfaces = 5;
        break;

    case LibUtilities::eTetrahedron:
        
        numedges = 6;
        numsurfaces = 4;
        break;

    
    case LibUtilities::eQuadrilateral:
        
        numedges = 4;
        numsurfaces = 1;
        break;

    case LibUtilities::eHexahedron:
        
        numedges = 12;
        numsurfaces = 6;
        break;
        

    
    default: cout<<"\n unknown shape type\n\n";exit(0);
    
    }
    
    const auto totPoints = (unsigned) E->GetTotPoints();

    Array<OneD, NekDouble> x = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> y = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> z = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> dx = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> dy = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> dz = Array<OneD, NekDouble>(totPoints);
    Array<OneD, NekDouble> sol = Array<OneD, NekDouble>(totPoints);


    switch (dimension)
    {
    case 1:
        {
            E->GetCoords(x);
            break;
        }

    case 2:
        {
            E->GetCoords(x, y);
            break;
        }

    case 3:
        {
            E->GetCoords(x, y, z);
            break;
        }
    default:
        break;
    }
    //get solution array
    // for (int i = 0; i < totPoints; ++i)
    // {
    //     if(dimension ==2) //hex
    //     {
    // 	  sol[i] = pow((x[i]-0.5),2)+ pow((y[i]+0.6),2) +  demo.GetEps();
    // 	  NekDouble eps =  demo.GetEps();
    // 	  string s = to_string(eps);
    // 	  funcdef =" pow((x[i]-0.5),2)+ pow((y[i]+0.6),2) + demo.GetEps();";
    // 	}
    // }

      storage2d = E->GetPhysEvaluateStorage();

      // timer.Start();
      // int dimension = E->GetShapeDimension(); 
      // demo.testcoord2dtqpts = demo.GetCoords(E);
      // demo.testcoord2dtqmidpts = demo.GetQuadratureMidCoords(E);
      // demo.testcoord2dtlattice = demo.GetLatticeCoords(E);
      

      //For 2D:
      Array<OneD, Array<OneD, NekDouble> > coords(2);
      coords[0] = Array<OneD, NekDouble>(2);
      coords[1] = Array<OneD, NekDouble>(2);
      Array<OneD, NekDouble> out_eval;
      coords[0][0] = -1.0;
      coords[1][0] = 1.0;

      coords[0][1] = -0.0641299;
      coords[1][1] = -0.73172;
      Array<OneD, NekDouble> out_d0(coords[0].size()*E->GetNcoeffs());
      Array<OneD, NekDouble> out_d1(coords[0].size()*E->GetNcoeffs());


	out_eval = E->PhysEvaluateBasis(  coords, storage2d, out_d0, out_d1, NullNekDouble1DArray   );

	cout<<"\n out_eval: ";
      for(int k =0; k < out_eval.size(); k++)
	cout<<" "<<out_eval[k]<<" ";
      cout<<"\n";
	cout<<"\n out_d0: ";
      for(int k =0; k < out_d0.size(); k++)
	cout<<" "<<out_d0[k]<<" ";
      cout<<"\n";
	cout<<"\n out_d1: ";
      for(int k =0; k < out_d1.size(); k++)
	cout<<" "<<out_d1[k]<<" ";
      cout<<"\n";
}
