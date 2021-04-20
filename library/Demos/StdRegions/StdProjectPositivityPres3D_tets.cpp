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

//Modification to deal with exact solution for diff. Return 1 if integer < 0.
static double pow_loc(const double val, const int i)
{
  return (i < 0) ? 1.0 : pow(val, i);
}


//declare Do_optimize
void Do_optimize(Array<OneD, NekDouble> &uhats);

// declare caller routine to find_roots
// flag = 0 -> opt_needed calls
// flag = 1 -> sphere_rot calls
Array<OneD, Array< OneD,  NekDouble> >  call_find_roots(Array<OneD, NekDouble> &uhatsall, NekDouble &avgiterGD, int d = 0, Array< OneD, Array<OneD, NekDouble> >&uhatsedges =NullNekDoubleArrayofArray , Array< OneD, Array<OneD, NekDouble> >&surfaceuhats =NullNekDoubleArrayofArray );

//declare Opt_needed
int Opt_needed(Array<OneD, NekDouble> uhats, int flag = 0);


// for 2D elements to get uhats at edges
// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble> uhats,Array<OneD, Array<OneD, NekDouble> >&ret , int d = 0);

void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret, int modes, Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1 = NullNekDouble1DArray , Array<OneD, NekDouble> Vxyd2 = NullNekDouble1DArray);


void surfuhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz );


void deruhatsedges(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble >&ret,  Array<OneD, NekDouble> Vxyz, Array<OneD, NekDouble> Vdxyz, Array<OneD, NekDouble> Vxdyz, Array<OneD, NekDouble> Vxydz);

void project_surfaces( Array<OneD, NekDouble>uhats,  Array<OneD, Array<OneD, NekDouble> >&ret, int d = 0);

vector< vector<NekDouble> > gradient_descent(Array<OneD, NekDouble> uhats, Array<OneD, NekDouble> evalxyz, int sig, int surfflag , int surfid = 0);


FILE *csvfile(int flag);
void WriteCSV( Array<OneD, NekDouble>  uhats, int flag = 0);


Array<OneD, Array<OneD, NekDouble> > storage3d;
Array<OneD, Array<OneD, NekDouble> > storage2dtri;

int iterstaken;

string funcdef;
// Colleague matrix
int numedges, numsurfaces;
StdExpansion *E;
StdExpansion *E3seg;
StdExpansion *Etri;

//StdExpansion *E3quad;
Array<OneD, NekDouble> qZin; //inner point grid for 2D rootfinding
Array<OneD, Array<OneD, NekDouble> > edgeptsin;
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
NekDouble startval, startcoordx, startcoordy, startcoordz,  itersGD1, itersGD2, itersGD3;
            
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
  if(dimension < 3 )
    {
      cout<<"\n dimension should be 3, try using StdProjectPositivityPres1D or StdProjectPositivityPres2D for other dimensions\n\n";
      exit(0);
    }

  std::vector<int> order;
  std::vector<BasisType> btype(3, eNoBasisType);
  LibUtilities::ShapeType stype = E->DetShapeType();
  LibUtilities::PointsType pointsTypeCheb = LibUtilities::eGaussGaussChebyshev;    
  LibUtilities::PointsType pointsTypetriA = LibUtilities::eGaussLobattoLegendre;    

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
  for (int i = 0; i < totPoints; ++i)
    {
      if(dimension ==3) //hex
        {
	  sol[i] = pow((x[i]-1),2)+ pow((y[i]+1),2) + pow((z[i]+1),2) + demo.GetEps();//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i]) +0.75;//pow((x[i]+0.3),2)+ pow((y[i]-0.2),2) + pow((z[i]+0.5),2) + demo.GetEps();//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i]) +0.5;//pow((x[i]+0.5),2)+ pow((y[i]-0.5),2) + pow((z[i]-0.2),2) + demo.GetEps();//sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i]) +0.2  ;///sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(x[i])*sin(y[i])*sin(z[i]);//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(x[i])*sin(y[i])*sin(z[i]);//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(x[i])*sin(y[i])*sin(z[i]);//sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);//sin(x[i])*sin(y[i])*sin(z[i])+0.1;//sin(x[i])*sin(y[i])*sin(z[i]);// +0.1;//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(x[i])*sin(y[i])*sin(z[i])+0.1;//sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);//sin(x[i])*sin(y[i])*sin(z[i]);//sin(x[i])*sin(y[i])*sin(z[i]) +0.1;//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);//sin(x[i])*sin(y[i])*sin(z[i]);//sin(x[i])*sin(y[i])*sin(z[i]) ;//sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);//sin(x[i])*sin(y[i])*sin(z[i]);//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2;////sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1;
	  NekDouble eps =  demo.GetEps();
	  string s = to_string(eps);
	  funcdef =" pow((x[i]+0.5),2)+ pow((y[i]+1),2) + pow((z[i]+0.2),2) + demo.GetEps();";//"pow((x[i]+0.3),2)+ pow((y[i]-0.2),2) + pow((z[i]+0.5),2) +  demo.GetEps()";;//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i]) +0.75";//"pow((x[i]+0.y),2)+ pow((y[i]-n0.2),2) + pow((z[i]+0.5),2) + demo.GetEps();";//sin(M_PI*x[i])*sin(M_PI*y[";;//sin(x[i])*sin(y[i])*sin(z[i]);";//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.5";//sin(x[i])*sin(y[i])*sin(z[i]);";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])";//"sin(x[i])*sin(y[i])*sin(z[i])+0.1";//"sin(x[i])*sin(y[i])*sin(z[i]) ";//+0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i])+0.1";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(x[i])*sin(y[i])*sin(z[i]) +0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])";//"sin(x[i])*sin(y[i])*sin(z[i]) ";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i]);";//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1;";//"sin(x[i])*sin(y[i])*sin(z[i];";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2 ";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.5";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);";// "sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])  -2";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";// -  exp(-1) + 1.37";
	}
    }

  Array<OneD, NekDouble> phys(totPoints);
  Array<OneD, NekDouble> coeffs((unsigned) E->GetNcoeffs());
  LibUtilities::Timer     timer;
  //    NekDouble elapsed       = 0.0;

  //Project onto expansion
  E->FwdTrans(sol, coeffs);

  //Backward transform solution to get projected values
  E->BwdTrans(coeffs, phys);
    
  if (vm.count("optm"))
    {
      storage3d = E->GetPhysEvaluateStorage();

      timer.Start();
      int dimension = E->GetShapeDimension(); 
      vector<LibUtilities::BasisKey> tmpBasisKey = demo.GetBasisKey();       
      vector<LibUtilities::PointsKey> tmpPtsKey = demo.GetPointsKey();       
      LibUtilities::PointsKey pkeycheb(E->GetBasis(0)->GetNumPoints(), pointsTypeCheb);
      LibUtilities::PointsKey pkeytri(E->GetBasis(0)->GetZ().size(), pointsTypetriA);
      LibUtilities::BasisKey bkeycheb(LibUtilities::eChebyshev, 3*(E->GetBasis(0)->GetNumModes()),  pkeycheb);
      LibUtilities::BasisKey bkeytriB(LibUtilities::eOrtho_B,(E->GetBasis(0)->GetNumModes()), pkeytri);
      LibUtilities::BasisKey bkeytriA(LibUtilities::eOrtho_A, (E->GetBasis(0)->GetNumModes()), pkeytri);
      E3seg = new StdSegExp(bkeycheb);
       
      if (E3seg == nullptr)
        {
	  return 1;
        }
      Etri = new StdTriExp(bkeytriA, bkeytriB); 
      if (Etri == nullptr)
        {
	  return 1;
        }

      storage2dtri = Etri->GetPhysEvaluateStorage();
      demo.testcoord2dtqpts = demo.GetCoords(Etri);
      demo.testcoord2dtqmidpts = demo.GetQuadratureMidCoords(Etri, demo.testcoord2dtqpts);
      demo.testcoord2dtlattice = demo.GetLatticeCoords(demo.testcoord2dtqpts,  demo.testcoord2dtqmidpts);
      demo.testcoord3dqpts  =  demo.GetCoords(E); 
      demo.testcoord3dqmidpts  = demo.GetQuadratureMidCoords(E, demo.testcoord3dqpts);
      demo.testcoord3dlattice = demo.GetLatticeCoords(demo.testcoord3dqpts,  demo.testcoord3dlattice);
	

      //For 2D:
        
      demo.interioreval2dtqmidpts = Etri->PhysEvaluateBasis(  demo.testcoord2dtqmidpts, storage2dtri, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray   );

      // For 3D:
      demo.interioreval3dqmidpts = E->PhysEvaluateBasis(demo.testcoord3dqmidpts, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray   );

            
      if(numedges == 6) //tet 
	{
	  Array<OneD, NekDouble> edgexyztemp =   E3seg->GetBasis(0)->GetZ();

	  int totszedges1d =  edgexyztemp.size()*(E->GetNcoeffs());
	  edgeptsin = Array<OneD, Array<OneD, NekDouble> >(dimension);         
	  for(int p = 0; p < dimension; p++)
	    {
	      edgeptsin[p] =   Array<OneD, NekDouble>(edgexyztemp);            
	    }

	  // edge front left (AD) (x = -1) (y = -1)
	  Vxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
	  Vdxxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
	  Vdyxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
	  Vdzxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
	  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
	  edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

	  Vxm1ym1z = E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);  
                
	  //edge front hypt (DB) (y = -1) (z = -x)
	  Vym1xmz = Array<OneD, NekDouble>(totszedges1d);
	  Vdxym1xmz = Array<OneD, NekDouble>(totszedges1d);
	  Vdyym1xmz = Array<OneD, NekDouble>(totszedges1d);
	  Vdzym1xmz = Array<OneD, NekDouble>(totszedges1d);

	  Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[2][0], 1, &edgeptsin[0][0], 1);

	  Vym1xmz = E->PhysEvaluateBasis(edgeptsin,storage3d, Vdxym1xmz,Vdyym1xmz,Vdzym1xmz);  

                
	  //edge front bot (AB) (y = -1) (z = -1)
	  Vym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
	  Vdxym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
	  Vdyym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
	  Vdzym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
	  edgeptsin[0] = edgexyztemp;
	  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);           
            
	  Vym1xzm1 = E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);  

	  //edge left hypt (DC) ( x = -1) (z = -y)
	  Vxm1ymz = Array<OneD, NekDouble>(totszedges1d);
	  Vdxxm1ymz = Array<OneD, NekDouble>(totszedges1d);
	  Vdyxm1ymz = Array<OneD, NekDouble>(totszedges1d);
	  Vdzxm1ymz = Array<OneD, NekDouble>(totszedges1d);
	  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);           
	  edgeptsin[1] = edgexyztemp;
	  Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[1][0], 1, &edgeptsin[2][0], 1);
            
	  Vxm1ymz = E->PhysEvaluateBasis(edgeptsin,storage3d, Vdxxm1ymz, Vdyxm1ymz, Vdzxm1ymz);  

	  // edge bot diag (BC) (z = -1) (y = -x)  
	  Vxmyzm1   = Array<OneD, NekDouble>(totszedges1d);
	  Vdxxmyzm1 = Array<OneD, NekDouble>(totszedges1d)   ;
	  Vdyxmyzm1 = Array<OneD, NekDouble>(totszedges1d)    ;
	  Vdzxmyzm1 = Array<OneD, NekDouble>(totszedges1d)    ;
	  edgeptsin[0] = edgexyztemp;
	  Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[0][0], 1, &edgeptsin[1][0], 1); 
	  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[2].size(), -1.0);    
	  Vxmyzm1 = E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxxmyzm1, Vdyxmyzm1, Vdzxmyzm1);  

	  //edge CA bot left (x = -1) (z = -1)
	  Vxm1yzm1   = Array<OneD, NekDouble>(totszedges1d);
	  Vdxxm1yzm1 = Array<OneD, NekDouble>(totszedges1d)   ;
	  Vdyxm1yzm1 = Array<OneD, NekDouble>(totszedges1d)    ;
	  Vdzxm1yzm1 = Array<OneD, NekDouble>(totszedges1d)    ;
	  edgeptsin[1] = edgexyztemp;
	  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[2].size(), -1.0);    
                
	  Vxm1yzm1 = 		E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);  

	}
      if(numsurfaces == 4) // tet
	{

	  int totsurf2d =  (demo.testcoord2dtqpts[0].size())*E->GetNcoeffs();
	  Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);
	      
	  for(int k = 0; k < dimension-1; k++)
	    {
	      surfptsin[k] = Array<OneD, NekDouble>(demo.testcoord2dtqpts[k]);
	      surfptsintemp[k] = Array<OneD, NekDouble>(demo.testcoord2dtqpts[k]);
	    }

	  surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.testcoord2dtqpts[0].size()); 
	  surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.testcoord2dtqpts[0].size());
	  int totpt = surfptsin[0].size();
	  //surface bot z = -1, (ABC)
	  Vxmyzm1surf = Array<OneD, NekDouble>(totsurf2d);
	  surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0); 
	        
	  surfptsintemp[0] = surfptsin[0]; 
	  surfptsintemp[1] = surfptsin[1]; 
	  Vxmyzm1surf = E->PhysEvaluateBasis(surfptsintemp, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

	  //surface left x = -1  (DAC)
	  Vxm1yz =  Array<OneD, NekDouble>(totsurf2d);  
	  surfptsintemp[0] = Array<OneD, NekDouble> (totpt, -1.0);
	  surfptsintemp[1] = surfptsin[0];
	  surfptsintemp[2] = surfptsin[1];
        
	  Vxm1yz = E->PhysEvaluateBasis(surfptsintemp, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
		
	  //surf front y = -1,  (DAB)
	  Vxmzym1 =  Array<OneD, NekDouble>(totsurf2d);
	  surfptsintemp[1] = Array<OneD, NekDouble> (totpt, -1.0);   
	  surfptsintemp[0] = surfptsin[0];
	  surfptsintemp[2] = surfptsin[1];
	  Vxmzym1 = E->PhysEvaluateBasis(surfptsintemp, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

	  //surf DCB (x + y + z = -1), 
	  Vxyzm1 =  Array<OneD, NekDouble>(totsurf2d);
	  surfptsintemp[1] = surfptsin[1];;
	  Vmath::Vadd(totpt, &surfptsin[0][0], 1, &surfptsin[1][0], 1, &surfptsintemp[2][0], 1);
	  Vmath::Smul(totpt, -1.0, &surfptsintemp[2][0], 1, &surfptsintemp[2][0], 1);
	  Vmath::Sadd(totpt, -1.0, &surfptsintemp[2][0], 1, &surfptsintemp[2][0], 1); 
	  Vxyzm1 = E->PhysEvaluateBasis(surfptsintemp,storage3d , NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
		
	}
      timer.Stop();
      //	elapsed  = ;
      cout<<"\n setup phase took "<<timer.TimePerTest(1)<<"s\n\n";
      timer.Start();
	
      timer.Start();
	
	
      int opt_need = Opt_needed(coeffs);            
	
      timer.Stop();
      NekDouble timeneg = timer.TimePerTest(1);  
      cout<<"\n checking for -vity took "<<timeneg<<"s\n\n";

	
      if(opt_need)
        {

	  cout<<"\n need optimization";
	  //WriteCSV(coeffs);
	  timer.Start();
	  Do_optimize(coeffs);
	  timer.Stop();
	  //	      elapsed += timer.TimePerTest(1);
	  cout<<"\n optimizertook "<<timer.TimePerTest(1)<<"s\n\n";
            
	  cout<<"\n doopt done\n verifying at points on lattice...\n";
	    
	  timer.Start();
	  opt_need = Opt_needed(coeffs, 1 );
	  timer.Stop();
	  //	    timeverify = ;
	  cout<<"\n Verification took "<<timer.TimePerTest(1)<<"s\n\n";
	    
	  if(opt_need)
            {
	      
	      cout<<"func="<<funcdef;
	      
	      cout<<"\n fail\n";
	      //cout<<" "<<" "<<startval<<", ("<<startcoordx<<" "<<startcoordy<<" "<<startcoordz<<"),fail, "<<" "<<itersGD1<<", "<<itersGD2<<", "<<itersGD3<<", "<<timeneg<<", "<<timeverify<<","<<elapref<<" \n";;

	      exit(0);
            }
	  else
            {

	      
	      cout<<"func="<<funcdef;
	      
	      cout<<"\n pass\n\n";
	      //	      cout<<" "<<" "<<startval<<", ("<<startcoordx<<" "<<startcoordy<<" "<<startcoordz<<"), , pass, "<<" "<<itersGD1<<", "<<itersGD2<<", "<<itersGD3<<", "<<iterstaken<<", "<<timeneg<<", "<<timeverify<<", "<<elapref<<"\n";
 	      //WriteCSV(coeffs, 1);
 
	      sol = phys; 
            }
	}
      else{
	cout<<"\n optimizer no need\n\n";
      }
      //Backward transform solution to get projected values
      E->BwdTrans(coeffs, phys);

    }

  //Calculate L_inf & L_2 error
  cout << "L infinity error: \t" << E->Linf(phys, sol) << endl;
  if (stype != ePoint)
    {
      cout << "L 2 error: \t \t \t" << E->L2(phys, sol) << endl;
    }

  if (!vm.count("diff") && stype != ePoint)
    {
      //Evaluate solution at x = y = 0 and print error
      Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
      t[0] = -0.5;
      t[1] = -0.25;
      t[2] = -0.3;
      sol[0] = Shape_sol(t[0], t[1], t[2], order, btype, stype, false);

      NekDouble nsol = E->PhysEvaluate(t, phys);

      cout << "Error at x = (";
      for (int i = 0; i < dimension; ++i)
        {
	  cout << t[i] << ", ";
        }
      cout << "\b\b): " << nsol - sol[0] << endl;
    }

  // Calculate integral of known function to test different quadrature
  // distributions on each element.
  for (int i = 0; i < totPoints; ++i)
    {
      sol[i] = dimension == 1 ? exp(x[i]) : dimension == 2 ?
	exp(x[i]) * sin(y[i]) : exp(x[i] + y[i] + z[i]);
    }

  NekDouble exact = 0.0;
  switch(stype)
    {
    case eSegment:
      exact = M_E - 1.0 / M_E;
      break;
    case eTriangle:
      exact = -0.5 * (sin(1.0) + cos(1.0) + M_E * M_E *
		      (sin(1.0) - cos(1.0))) / M_E;
      break;
    case eQuadrilateral:
      exact = 2.0 * (M_E - 1.0 / M_E) * sin(1.0);
      break;
    case eTetrahedron:
      exact = 1.0 / M_E - 1.0 / M_E / M_E / M_E;
      break;
    case ePrism:
      exact = M_E - 1.0 / M_E / M_E / M_E;
      break;
    case ePyramid:
      exact = - 1.0 / M_E / M_E / M_E - 4.0 / M_E + M_E;
      break;
    case eHexahedron:
      exact = pow((M_E * M_E - 1.0) / M_E, 3.0);
      break;
    default:
      ASSERTL0(false, "Exact solution not known.");
      break;
    }
  std::cout << "Integral error: " << fabs(exact - E->Integral(sol))
	    << std::endl;

  return 0;
}


void project_surfaces( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret, int d)
{
  boost::ignore_unused(d);

  //surface bot z = -1 (ABC)
  surfuhats(uhats, ret[0], Vxmyzm1surf);
    
  //surface left x = -1 (DAC)  
  surfuhats(uhats, ret[1], Vxm1yz);

  //surf front y = -1  (DAB)
  surfuhats(uhats, ret[2], Vxmzym1);
    
			
    
  //surf DCB (x + y + z = -1)
  surfuhats(uhats, ret[3], Vxyzm1);
    
}


// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret , int d)
{
  // assert d = 0 or 1
  if(d == 0)
    {   
      int modes =  E3seg->GetNcoeffs();//(3*(E->GetBasis(0)->GetNumModes()));

      for(int k= 0; k < numedges; k++)
        {
	  ret[k] = Array<OneD, NekDouble>(modes);
        }
    }
  if(numedges == 6) // tet
    {

      if(d == 1)
        {

	  // edge front left (AD) (x = -1) (y = -1)
	  edgederpquhats(uhats, ret[0], E->GetNcoeffs(),  Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);
                        
	  //edge front hypt (DB) (y = -1) (z = -x) 
	  edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vym1xmz, Vdxym1xmz, Vdyym1xmz, Vdzym1xmz  );
            
	  //edge front bot (AB) (y = -1) (z = -1)  
	  edgederpquhats(uhats, ret[2], E->GetNcoeffs(), Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1 );

	  //edge left hypt (DC) ( x = -1) (z = -y)  
	  edgederpquhats(uhats, ret[3], E->GetNcoeffs(), Vxm1ymz, Vdxxm1ymz, Vdyxm1ymz, Vdzxm1ymz);
     
	  // edge bot diag (BC) (z = -1) (y = -x)
	  edgederpquhats(uhats, ret[4], E->GetNcoeffs(),  Vxm1ymz, Vdxxm1ymz, Vdyxm1ymz, Vdzxm1ymz);

	  //edge CA bot left (x = -1) (z = -1)
	  edgederpquhats(uhats, ret[5], E->GetNcoeffs(),Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);
            
     
        }
      else if( d == 0)
        {
	  // edge front left (AD) (x = -1) (y = -1)
	  deruhatsedges(uhats, ret[0],  Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);

	  //edge front hypt (DB) (y = -1) (z = -x)
	  deruhatsedges(uhats, ret[1],Vym1xmz, Vdxym1xmz, Vdyym1xmz, Vdzym1xmz);
            
	  //edge front bot (AB) (y = -1) (z = -1)
	  deruhatsedges(uhats, ret[2],  Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);

	  //edge left hypt (DC) ( x = -1) (z = -y)
	  deruhatsedges(uhats, ret[3], Vxm1ymz, Vdxxm1ymz, Vdyxm1ymz, Vdzxm1ymz);

	  // edge bot diag (BC) (z = -1) (y = -x)
	  deruhatsedges(uhats, ret[4],  Vxm1ymz, Vdxxm1ymz, Vdyxm1ymz, Vdzxm1ymz);
                   
	  //edge CA bot left (x = -1) (z = -1)
	  deruhatsedges(uhats, ret[5], Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);
          

        }
    }
   
}

// FILE *csvfile(int flag)
// {
//     static char name[LOGNAME_SIZE];
//     char fullname[] = "Pri_";//[LOGNAME_SIZE+4];                                                                                                                                                           
//     time_t now = time(0);
//     strftime(name, sizeof(name), LOGNAME_FORMAT, localtime(&now));
//     char * newArray = new char[std::strlen(fullname)+std::strlen(name)+10];
//     std::strcpy(newArray,fullname);
//     std::strcat(newArray,name);
//     char integer_string[20];
//     sprintf(integer_string, "_P%d_Q%d_%d.csv" , E->GetBasis(0)->GetNumModes(), E->GetBasis(0)->GetNumPoints(), flag);
//     std::strcat(newArray,integer_string);
//     cout<<"\n"<< newArray<<"\n\n";
//     return fopen(newArray, "w");

// }

// void WriteCSV(Array<OneD, NekDouble> uhats, int flag)
// {

//     // if flag == 0, before opt
//     // if flag == 1, after opt
//     FILE *file = csvfile(flag);
//     int uhatstot = uhats.size();
//     int totpts = demo.testcoord3dqpts[0].size();
//     Array<OneD, NekDouble> temp(uhatstot);
//     for(int i = 0; i<totpts; i++)
//     {
//         Vmath::Vmul(uhatstot, &demo.interioreval3d[i], totpts, &uhats[0], 1, &temp[0], 1);
//         NekDouble v = Vmath::Vsum(uhatstot, temp, 1);
//         string s = std::to_string(demo.testcoord3dqpts[0][i])+","+std::to_string(demo.testcoord3dqpts[1][i])+","+std::to_string(demo.testcoord3dqpts[2][i])+","+std::to_string(v)+"\n";
//         const char* c = s.c_str();
//         fprintf (file, c); 
//     }
//     fclose(file);
// }


void surfuhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz )
{
  int modes = uhats.size();

  Array<OneD, NekDouble> temp(uhats.size());  

  int totpts = Etri->GetTotPoints();
  Array<OneD, NekDouble> vals(totpts);    
    
  //Vxyz*uhats -> project to -> E3quad
  for(int k = 0; k < totpts; k++)
    {    
      Vmath::Vmul(uhats.size(), &Vxyz[k], totpts, &uhats[0], 1, &temp[0], 1);
      vals[k]  = Vmath::Vsum(modes, temp, 1);
    }
    
  Etri->FwdTrans(vals, ret);
    
}

void deruhatsedges(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz, Array<OneD, NekDouble> Vdxyz, Array<OneD, NekDouble> Vxdyz, Array<OneD, NekDouble> Vxydz)
{
  boost::ignore_unused(Vxyz);
  int totpts =  edgeptsin[0].size();//E3seg->GetTotPoints();
  int modes = uhats.size();
  NekDouble v1;   
  Array<OneD, NekDouble> vals(totpts), hold(totpts);
  Array<OneD, NekDouble> temp(uhats.size());  
    
  for(int k = 0; k < totpts; k++)
    {    
      Vmath::Vmul(uhats.size(), &Vdxyz[k], totpts, &uhats[0], 1, &temp[0], 1);
      v1  = Vmath::Vsum(modes, temp, 1);
      Vmath::Vmul(uhats.size(), &Vxdyz[k], totpts, &uhats[0], 1, &temp[0], 1);
      v1  = v1 + Vmath::Vsum(modes, temp, 1);  
        
      Vmath::Vmul(uhats.size(), &Vxydz[k], totpts, &uhats[0], 1, &temp[0], 1);
      vals[k]  = v1 + Vmath::Vsum(modes, temp, 1);
    }    
    
  E3seg->FwdTrans(vals,ret);
}
    


void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  int modes,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2)
{
  boost::ignore_unused(modes);
  int uhatstot = uhats.size();
  int totpts = edgeptsin[0].size();//E3seg->GetTotPoints();
  Array<OneD, NekDouble> temp(totpts), temp2(uhatstot), temp3(uhatstot);
  Array<OneD, NekDouble> pqeval(totpts);
  NekDouble v1, v2;
  for(int i = 0; i<totpts; i++)
    {
      Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd0[0]+i, totpts, &temp2[0], 1);
      v1  = Vmath::Vsum(uhatstot, temp2, 1); 
      Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);

      v2  = Vmath::Vsum(uhatstot, temp3, 1);  

      v1 = v2*v1;

      // At this point,                 

      // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

      // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
      Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &temp2[0], 1);

 
        
      v2  = Vmath::Vsum(uhatstot, temp2, 1);  
 
      Vmath::Vmul(uhats.size(), &Vxyd0[i], totpts, &uhats[0], 1, &temp2[0], 1);

      v1= v2*Vmath::Vsum(uhats.size(), temp2, 1)- v1;
  
      pqeval[i] = v1;
 

      if(dimension > 1)
        {
	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd1[0]+i, totpts, &temp2[0], 1);
	  v1  = Vmath::Vsum(uhatstot, temp2, 1); 
	  Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);

	  v2  = Vmath::Vsum(uhatstot, temp3, 1);  

	  v1 = v2*v1;

	  // At this point,                 

	  // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

	  // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &temp2[0], 1);

 
        
	  v2  = Vmath::Vsum(uhatstot, temp2, 1);  
 
	  Vmath::Vmul(uhats.size(), &Vxyd1[i], totpts, &uhats[0], 1, &temp2[0], 1);

	  v1=  v2*Vmath::Vsum(uhats.size(), temp2, 1) - v1;
 
	  pqeval[i] += v1;
 
        }
      if(dimension == 3)
        {
	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd2[0]+i, totpts, &temp2[0], 1);
	  v1  = Vmath::Vsum(uhatstot, temp2, 1); 
	  Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);

	  v2  = Vmath::Vsum(uhatstot, temp3, 1);  

	  v1 = v2*v1;

	  // At this point,                 

	  // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

	  // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &temp2[0], 1);

 
        
	  v2  = Vmath::Vsum(uhatstot, temp2, 1);  
 
	  Vmath::Vmul(uhats.size(), &Vxyd2[i], totpts, &uhats[0], 1, &temp2[0], 1);

	  v1= v2*Vmath::Vsum(uhats.size(), temp2, 1) - v1;
 
	  pqeval[i] += v1;
 
        }
    }

  E3seg->FwdTrans(pqeval, ret);
    
}


// sig = 0 -> called by opt_needed
// sig = 1 -> called by sphere_rotation
Array< OneD,  NekDouble>  call_find_roots(Array<OneD,  NekDouble> &uhats , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats,  NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime)
{

  int dimension = E->GetShapeDimension(); 

  NekDouble dummy;
  Array<OneD, NekDouble> temp(uhats.size());
  roots1dtime = 0.0;
  roots2dtime = 0.0;
  roots3dtime = 0.0;
  NekDouble  inf = numeric_limits<double>::infinity();
  minv = inf;
  NekDouble avgiterGDhold;
  Timer t;
  Array<OneD, Array<OneD, NekDouble > >  rethold(dimension);
  vector<vector<NekDouble> > retall(dimension);
  Array<OneD, NekDouble> ret(dimension);
    
  if(numedges == 6) //tet
    {

      for(int ii = 0; ii < numedges; ii++)
        {
	
    	  t.Start();
	  rethold  = (demo.find_roots(uhatsedges[ii], E, NullNekDoubleArrayofArray, dummy,  0, 0, 0, 0)) ;
	  t.Stop();

	  roots1dtime +=  t.TimePerTest(1);
	  for(int p = 0; p < rethold[0].size(); p++)
            {
	      switch(ii)
		{
		case 0:
		  retall[0].push_back(-1);
		  retall[1].push_back(-1);
		  retall[2].push_back(rethold[0][p]);
		  break;
		case 1:  //edge front hypt (DB) (y = -1) (z = -x)
                
		  
		  retall[0].push_back(rethold[0][p]);
		  retall[2].push_back(-rethold[0][p]);
		  retall[1].push_back(-1);
		  
		  retall[0].push_back(-rethold[0][p]);
		  retall[2].push_back(rethold[0][p]);
		  retall[1].push_back(-1);
		  break;
		case 2: //edge front bot (AB) (y = -1) (z = -1)  
		
                
		  retall[0].push_back(rethold[0][p]);
		  retall[1].push_back(-1);
		  retall[2].push_back(-1);
                  
                  
		  break;
		case 3:    //edge left hypt (DC) ( x = -1) (z = -y)    
                
		  retall[0].push_back(-1);
		  retall[1].push_back(rethold[0][p]);
		  retall[2].push_back(-rethold[0][p]);

		  retall[0].push_back(-1);
		  retall[1].push_back(-rethold[0][p]);
		  retall[2].push_back(rethold[0][p]);
		  break;
		case 4:   // edge bot diag (BC) (z = -1) (y = -x)  
                
		  retall[2].push_back(-1);
		  retall[1].push_back(-rethold[0][p]);
		  retall[0].push_back(rethold[0][p]);
		  
		  retall[2].push_back(-1);
		  retall[1].push_back(rethold[0][p]);
		  retall[0].push_back(-rethold[0][p]);
		  
		  break;
		case 5:  //edge CA bot left (x = -1) (z = -1)    
                
		  retall[1].push_back(rethold[0][p]);
		  retall[0].push_back(-1);
		  retall[2].push_back(-1);
		  break;
		default:
		  cout<<"\n wrong edge id\n";
		  exit(0);
	            
		}
            
	    }
            
        }
    }
  if(numsurfaces == 4) //tet
    {
      NekDouble avgiterGDhold;

      // call 2D rootfinder on each surface:
      for(int ii = 0; ii < numsurfaces; ii++)
        {
	  t.Start();   
	  rethold = demo.find_roots(surfaceuhats[ii], Etri, storage2dtri,  avgiterGDhold, 0, 1, 0 , 1); // last arg = 1 coz all 2d ele are triangles

	  t.Stop();
	  roots2dtime +=  t.TimePerTest(1);
	  avgiterGD += avgiterGDhold;

	  switch(ii)
	    {
	    case 0:     //surface bot z = -1, x + y = 0 (ABC)   
              
	      retall[0].push_back( (rethold[0][0]));
	      retall[1].push_back( (rethold[1][0]));
	      retall[2].push_back( (-1.0));
	      break;
	    case 1:
	      //surface left x = -1 y + z = 0 (DAC)  
		  
	      retall[0].push_back(-1);
	      retall[2].push_back(rethold[1][0]);
	      retall[1].push_back(rethold[0][0]);
	      break;
	    case 2:  //surf front y = -1, x + z = 0 (DAB)  
            
	      retall[2].push_back(rethold[1][0]);
	      retall[1].push_back(-1);
	      retall[0].push_back(rethold[0][0]);
	      break;
	    case 3: //surf DCB (x + y + z = -1)  
		  
	      
	      retall[2].push_back(-1.0 - rethold[1][0] - rethold[0][0]);
	      retall[1].push_back(rethold[1][0]);
	      retall[0].push_back(rethold[0][0]);
	      break;
	    default: cout<<"\n invalid surface id\n";
	      exit(0);
	    }
	}
    }
    
    
  //find interior roots:
  t.Start();
  rethold =  (demo.find_roots(uhats, E,storage3d,  avgiterGDhold, 0, 0, 1, 0)) ;
  t.Stop();
  //      cout<<"\n roots3dtime at ctr = "<<t.TimePerTest(1);
  roots3dtime +=  t.TimePerTest(1);  

  //cout<<" rethold = "<<rethold[0][0]<<" "<<rethold[1][0]<<" "<<rethold[2][0];
  avgiterGD += avgiterGDhold;
  Array<OneD, Array<OneD, NekDouble> > tmpcoord(dimension);
	  
  if(rethold[0][0] < inf && rethold[1][0] < inf && rethold[2][0] < inf)
    {
      retall[0].push_back(rethold[0][0]);
      retall[1].push_back(rethold[1][0]);
      retall[2].push_back(rethold[2][0]);
      // cout<<"\n 3d rootfinder ret: "<<rethold[0][0]<<" "<<rethold[1][0]<<" "<<rethold[2][0]<<"\n";
    }
  for(int p = 0; p < dimension; p++)
    {
      tmpcoord[p] = Array<OneD, NekDouble>(retall[0].size());
    }
  for(int p = 0; p < retall[0].size(); p++)
    {
      tmpcoord[0][p] = retall[0][p];
      tmpcoord[1][p] = retall[1][p];
      tmpcoord[2][p] = retall[2][p];
    }
  NekDouble tempmin = inf;
  Array<OneD, NekDouble> evalroots = E->PhysEvaluateBasis(tmpcoord, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
     
      
  Array<OneD, NekDouble> minvandx(dimension+1);
  Array<OneD, NekDouble> tmparr(retall[0].size());
      
      
  demo.pq(uhats, tmpcoord, evalroots, minvandx, tmparr);
  tempmin = minvandx[0];
  ret[0] = minvandx[1];
  ret[1] = minvandx[2];
  ret[2] = minvandx[3];
	
  minv = tempmin;
  return ret;
}


int Opt_needed(Array<OneD, NekDouble> uhats, int flag)
{
  // check for -ve vals on lattice (E->testcoord3dlat)
  // ? should we check on 2dqlattice as well? Find out
      
  int nq = E->GetTotPoints();//demo.testcoord3dlat[0].size();
  Array<OneD, NekDouble> temp1(nq);
  Array<OneD, NekDouble> temp2(uhats.size());
  Array<OneD, Array<OneD, NekDouble> > coordquad = demo.GetCoords(E);
  //interioreval3dlat*uhats (todo: change to blas  call if needed)
  NekDouble minv = numeric_limits<double>::infinity();
  int idxmin = -1;
  for (int i = 0; i<nq; i++)
    {
      Vmath::Vmul(uhats.size(), &storage3d[0][i], nq, &uhats[0], 1, &temp2[0], 1);
      temp1[i] = Vmath::Vsum(uhats.size(), temp2, 1);
      if(temp1[i] < minv)
	{
	  minv = temp1[i];
	  idxmin = i;
	}
    }
      
  //  cout<<demo.testcoord3dlattice[0][idxmin]<<" "<<demo.testcoord3dlattice[1][idxmin]<<" "<<demo.testcoord3dlattice[2][idxmin];
  int ct = nq;
  nq = demo.testcoord3dqmidpts[0].size();
  temp1=Array<OneD, NekDouble>(nq);
  for (int i = 0; i<nq; i++)
    {
      Vmath::Vmul(uhats.size(), &demo.interioreval3dqmidpts[i], nq, &uhats[0], 1, &temp2[0], 1);
      temp1[i] = Vmath::Vsum(uhats.size(), temp2, 1);
      if(temp1[i] < minv)
	{
	  minv = temp1[i];
	  idxmin = i+ct;
	}
    }
  if(minv < 0 && abs(minv) > 1e-10 && idxmin != -1)
    { 
      //      cout<<"\n minv = "<<minv<<" at idxmin = "<<idxmin<< " ";
      //cout<<demo.testcoord3dlattice[0][idxmin]<<" "<<demo.testcoord3dlattice[1][idxmin]<<" "<<demo.testcoord3dlattice[2][idxmin];
      if(flag == 0)
	{
	  startval = minv;
	  startcoordx = demo.testcoord3dlattice[0][idxmin];
	  startcoordy = demo.testcoord3dlattice[1][idxmin];
	  startcoordz = demo.testcoord3dlattice[2][idxmin];
	}
      cout<<"\n minv = "<<minv<<" at "<< demo.testcoord3dlattice[0][idxmin]<<  " "<<demo.testcoord3dlattice[1][idxmin]<<" "<< demo.testcoord3dlattice[2][idxmin];;;
	  
      return 1;   
	  
    }
  cout<<"\n minv = "<<minv<<" at "<< demo.testcoord3dlattice[0][idxmin]<<  " "<<demo.testcoord3dlattice[1][idxmin]<<" "<< demo.testcoord3dlattice[2][idxmin];;;

  return 0;
      
}




void Do_optimize(Array<OneD, NekDouble> &uhats)
{

  int dim = E->GetShapeDimension();
  double inf = numeric_limits<double>::infinity();
    
  //assert(size(constraints, 1) == N+2);
    
  int N1 = uhats.size();
  vector<Array<OneD,NekDouble> > d;
 
  d.push_back(uhats);
  int counter = 0;
      
  //    int NC = ck[1].size();          // number of constraints
  vector<double> tols;         // constraint specific tolerances

  int niter = 1e3;
  NekDouble minv;
  Array<OneD,NekDouble> optima(dim);
     
  //    int NC = 1; //number of constraints, only positivity for now
  tols.push_back(1e-13);
    
  Array<OneD, NekDouble>  pqvalcoords(dim+1), xastarr(dim), utemp(N1), wsp1(1);//xastarrhist(dim), wsp1;   
    

  Array<OneD, Array<OneD, NekDouble> > surfaceuhats  (numsurfaces), tmpcoord(dim);
  Array<OneD, Array<OneD, NekDouble> > Pf(numedges);
  for(int k= 0; k < numedges; k++)
    {
            
      Pf[k] = Array<OneD, NekDouble>(E3seg->GetNcoeffs());//(3*(E->GetBasis(0)->GetNumModes()-1));        
    }
  for(int k = 0; k < numsurfaces; k++)
    {
      surfaceuhats[k] = Array<OneD, NekDouble>(Etri->GetNcoeffs());
    }
        
  Array<OneD, NekDouble> Vtmp(  E->GetNcoeffs());    
  NekDouble pqval, timeprojectsurf = 0.0, timeprojectedges = 0.0;
  NekDouble avgiterGD = 0.0, avgiterhold = 0.0, roots1dtimehold = 0.0, roots2dtimehold = 0.0, roots3dtimehold = 0.0;//, evalrootshold;   
  Timer t;   
  while (counter <= niter)
    {
      NekDouble roots1dtime = 0.0, roots2dtime = 0.0, roots3dtime = 0.0 ;
      pqval = inf;
      utemp = d.back();
      if (counter > 0)
	{

	  t.Start();
	  //pq = @(xx) cfun(xx) .* (constraints(1,jj).*Vc(xx)*uhats);  
	  project_edges(utemp, Pf, 1);
	  t.Stop();
	  timeprojectedges += t.TimePerTest(1);

	  t.Start();
	  project_surfaces(utemp, surfaceuhats);
	  t.Stop();

	  timeprojectsurf += t.TimePerTest(1);
	  optima = call_find_roots(utemp, avgiterhold, Pf, surfaceuhats,  minv, roots1dtime, roots2dtime, roots3dtime );

	  roots1dtimehold += roots1dtime;
	  roots2dtimehold += roots2dtime;

	  roots3dtimehold += roots3dtime;
	  //	  cout<<"\n counter ="<<counter<<" roots3dtimehold="<<roots3dtimehold<<"\n";	
 	
	  //	cout<<"\n avgiterhold="<<avgiterhold<<" ";
	  avgiterGD += avgiterhold;    
    
	  //        wsp1 = Array<OneD, NekDouble>(optima[0].size());
	  // cout<<"\n optima:\n";
	  // for(int pp = 0; pp < optima.size(); pp++)
	  // {
	  //   //    for(int jj = 0; jj < optima[pp].size(); jj++)
	  //   //{
	  //         cout<<optima[pp]<<" ";
	  // 	  // }
	  // 	  // cout<<"\n";
	  // }
	  //        cout<<"\n";
	  // Vtmp is evaluation of basis fn at optima
	  for(int k = 0; k < dim; k++)
	    {
	      tmpcoord[k] = Array<OneD, NekDouble>(1, optima[k]);
	    }
	
	  //	t.Start();

	  Vtmp = E->PhysEvaluateBasis(tmpcoord, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        
	  //	  demo.pq(utemp, tmpcoord, Vtmp, pqvalcoords, wsp1);

	}// end if(counter > 0)
      else //if(counter == 0)
	{
	  for(int k = 0; k < dim; k++)
	    {
	      tmpcoord[k] = Array<OneD, NekDouble>(1);
	    }

	  tmpcoord[0][0] = startcoordx;
	  tmpcoord[1][0] = startcoordy;
	  tmpcoord[2][0] = startcoordz;
	  optima[0] = startcoordx;
	  optima[1] = startcoordy;
	  optima[2] = startcoordz;
	  Vtmp = E->PhysEvaluateBasis(tmpcoord, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        

	}
      demo.pq(utemp, tmpcoord, Vtmp, pqvalcoords, wsp1);
      // cout<<"\n fvals: "<<wsp1.size()<<"\n";;
      // for(int kk = 0 ; kk<wsp1.size(); kk++)
      //     cout<<" "<<wsp1[kk]<<" ";
      // cout<<"\n";

      if (pqvalcoords[0] < pqval)
        {
	  xastarr = optima;
	  pqval = pqvalcoords[0];
        }
        
      cout<<"\n at counter="<<counter<<" min val="<<pqval<<" xast ="<<xastarr[0]<<" "<<xastarr[1]<<" "<<xastarr[2]<<
	" ";
      // If minimum is non-negative, we're done
      if (pqval >= -tols.at(0))
        {
	  break;
        }
      Array<OneD, NekDouble>  Vastsq (N1);
      Array<OneD, NekDouble>  Vast (N1);
      Array<OneD, Array<OneD, NekDouble> > xastarrofarr(3);
      xastarrofarr[0] = Array<OneD, NekDouble>(1, xastarr[0]);
      xastarrofarr[1] = Array<OneD, NekDouble>(1, xastarr[1]);
      xastarrofarr[2] = Array<OneD, NekDouble>(1, xastarr[2]);
      Array<OneD, NekDouble> tmp;
      NekDouble vastsqsum;
      
      for( int i = 0; i < N1; i++)
        {
	  tmp = E->PhysEvaluateBasis(xastarrofarr, storage3d, i);
	  Vast[i] = tmp[0];
	  Vastsq[i] = (Vast[i]*Vast[i]);
	  
        }
      vastsqsum = Vmath::Vsum(N1, &Vastsq[0], 1);

      Array<OneD, NekDouble>  qast(N1);

      for(int i = 0; i<N1; i++)
        {
	  qast[i] = ((1/sqrt(vastsqsum))*(Vast[i]));
        }
      Vmath::Smul(N1, pqval, &qast[0], 1, &qast[0], 1);
        
      Vmath::Vsub(utemp.size(), &utemp[0], 1, &qast[0], 1, &qast[0], 1);
      d.push_back(qast);
      //	t.Stop();
      //evalrootshold+= t.TimePerTest(1);  

      counter = counter + 1;
    
    
    }
  roots1dtimehold = roots1dtimehold/(counter-1);
  roots2dtimehold = roots2dtimehold/(counter-1);
  roots3dtimehold = roots3dtimehold/(counter-1);
  timeprojectedges = timeprojectedges/(counter);
  timeprojectsurf = timeprojectsurf/(counter);    
  itersGD2 = (avgiterGD)/(counter);
  iterstaken = counter;
  cout<<"sphere_rotation took "<<counter<<"iterations, tets";//
  uhats = d.back();
  // cout<<"\n checking at (-1,-1,1):\n";
  // Array<OneD, Array<OneD, NekDouble> > coordtemp(3);
  // coordtemp[0] = Array<OneD, NekDouble>(1,-1);
  // coordtemp[1] = Array<OneD, NekDouble>(1,-1);
  // coordtemp[2] = Array<OneD, NekDouble>(1,1);
  // Array<OneD, NekDouble> valbarr = E->PhysEvaluateBasis(coordtemp, storage3d, NullNekDouble1DArray , NullNekDouble1DArray, NullNekDouble1DArray);
  // Array<OneD, NekDouble> tempval(uhats.size());
  // Vmath::Vmul(uhats.size(), uhats, 1, valbarr, 1, tempval, 1);
  // NekDouble valat = Vmath::Vsum(tempval.size(), tempval, 1);
  // cout<<" val = "<< valat<<"\n";
}


NekDouble Shape_sol(NekDouble x, NekDouble y, NekDouble z, vector<int> order,
                    vector<BasisType> btype, ShapeType stype, bool diff)
{
  map<ShapeType, function<int(int, const vector<int> &)>> shapeConstraint2;
  shapeConstraint2[ePoint] =
    [](int,   const vector<int> &     ) { return 1; };
  shapeConstraint2[eSegment] =
    [](int,   const vector<int> &     ) { return 1; };
  shapeConstraint2[eTriangle] =
    [](int k, const vector<int> &order) { return order[1] - k; };
  shapeConstraint2[eQuadrilateral] =
    [](int,   const vector<int> &order) { return order[1]; };
  shapeConstraint2[eTetrahedron] =
    [](int k, const vector<int> &order) { return order[1] - k; };
  shapeConstraint2[ePyramid] =
    [](int k, const vector<int> &order) { return order[1] - k; };
  shapeConstraint2[ePrism] =
    [](int,   const vector<int> &order) { return order[1]; };
  shapeConstraint2[eHexahedron] =
    [](int,   const vector<int> &order) { return order[1]; };

  map<ShapeType, function<int(int, int, const vector<int> &order)>>
    shapeConstraint3;
  shapeConstraint3[ePoint] =
    [](int,   int,   const vector<int> &     ) { return 1; };
  shapeConstraint3[eSegment] =
    [](int,   int,   const vector<int> &     ) { return 1; };
  shapeConstraint3[eTriangle] =
    [](int,   int,   const vector<int> &     ) { return 1; };
  shapeConstraint3[eQuadrilateral] =
    [](int,   int,   const vector<int> &     ) { return 1; };
  shapeConstraint3[eTetrahedron] =
    [](int k, int l, const vector<int> &order) { return order[2] - k - l; };
  shapeConstraint3[ePyramid] =
    [](int k, int l, const vector<int> &order) { return order[2] - k - l; };
  shapeConstraint3[ePrism] =
    [](int k, int,   const vector<int> &order) { return order[2] - k; };
  shapeConstraint3[eHexahedron] =
    [](int,   int,   const vector<int> &order) { return order[2]; };

  NekDouble sol = 0.0;
  if (!diff)
    {
      if (btype[0] == eFourier && stype == eSegment)
        {
	  for (int k = 0; k < order[0] / 2 - 1; ++k)
            {
	      sol += sin(k * M_PI * x) + cos(k * M_PI * x);
            }
        }
      else if (btype[0] == eFourierSingleMode && stype == eSegment)
        {
	  sol += 0.25 * sin(M_PI * x) + 0.25 * cos(M_PI * x);
        }
      else if (btype[0] == eFourier && stype == eQuadrilateral)
        {
	  if (btype[1] == eFourier)
            {
	      for (int k = 0; k < order[0] / 2; ++k)
                {
		  for (int l = 0; l < order[1] / 2; ++l)
                    {
		      sol += sin(k * M_PI * x) * sin(l * M_PI * y) +
			sin(k * M_PI * x) * cos(l * M_PI * y) +
			cos(k * M_PI * x) * sin(l * M_PI * y) +
			cos(k * M_PI * x) * cos(l * M_PI * y);
                    }
                }
            }
	  else if (btype[1] == eFourierSingleMode)
            {
	      for (int k = 0; k < order[0] / 2; ++k)
                {
		  sol += sin(k * M_PI * x) * sin(M_PI * y) +
		    sin(k * M_PI * x) * cos(M_PI * y) +
		    cos(k * M_PI * x) * sin(M_PI * y) +
		    cos(k * M_PI * x) * cos(M_PI * y);
                }
            }
	  else
            {
	      for (int k = 0; k < order[0] / 2; ++k)
                {
		  for (int l = 0; l < order[1]; ++l)
                    {
		      sol += sin(k * M_PI * x) * pow_loc(y, l) +
			cos(k * M_PI * x) * pow_loc(y, l) ;
                    }
                }
            }
        }
      else if (btype[0] == eFourierSingleMode && stype == eQuadrilateral)
        {
	  if (btype[1] == eFourier)
            {
	      for (int l = 0; l < order[1] / 2; ++l)
                {
		  sol += sin(M_PI * x) * sin(l * M_PI * y) +
		    sin(M_PI * x) * cos(l * M_PI * y) +
		    cos(M_PI * x) * sin(l * M_PI * y) +
		    cos(M_PI * x) * cos(l * M_PI * y);
                }

            }
	  else if (btype[1] == eFourierSingleMode)
            {
	      sol += sin(M_PI * x) * sin(M_PI * y) +
		sin(M_PI * x) * cos(M_PI * y) +
		cos(M_PI * x) * sin(M_PI * y) +
		cos(M_PI * x) * cos(M_PI * y);
            }
	  else
            {
	      for (int l = 0; l < order[1]; ++l)
                {
		  sol += sin(M_PI * x) * pow_loc(y, l) +
		    cos(M_PI * x) * pow_loc(y, l);
                }
            }
        }
      else if (btype[1] == eFourier && stype == eQuadrilateral)
        {
	  for (int k = 0; k < order[0]; ++k)
            {
	      for (int l = 0; l < order[1] / 2; ++l)
                {
		  sol += sin(l * M_PI * y) * pow_loc(x, k) +
		    cos(l * M_PI * y) * pow_loc(x, k);
                }
            }
        }
      else if (btype[1] == eFourierSingleMode && stype == eQuadrilateral)
        {
	  for (int k = 0; k < order[0]; ++k)
            {
	      sol += sin(M_PI * y) * pow_loc(x, k) +
		cos(M_PI * y) * pow_loc(x, k);
            }
        }
      else
        {
	  for (int k = 0;
	       k < order[0]; ++k) //ShapeConstraint 1 is always < order1
            {
	      for (int l = 0; l < shapeConstraint2[stype](k, order); ++l)
                {
		  for (int m = 0;
		       m < shapeConstraint3[stype](k, l, order); ++m)
                    {
		      sol += pow_loc(x, k) * pow_loc(y, l) * pow_loc(z, m);
                    }
                }
            }
        }
    }
  else if (diff)
    {
      if (btype[0] == eFourier && stype == eSegment)
        {
	  for (int k = 0; k < order[0] / 2 - 1; ++k)
            {
	      sol += k * M_PI * (cos(k * M_PI * z) - sin(k * M_PI * z));
            }
        }
      else if (btype[0] != eFourier && btype[1] == eFourier &&
	       stype == eQuadrilateral)
        {
	  for (int k = 0; k < order[0]; ++k)
            {
	      for (int l = 0; l < order[1] / 2; ++l)
                {
		  sol += k * pow_loc(x, k - 1) * sin(M_PI * l * y)
		    + M_PI * l * pow_loc(x, k) * cos(M_PI * l * y) +
		    +k * pow_loc(x, k - 1) * cos(M_PI * l * y)
		    - M_PI * l * pow_loc(x, k) * sin(M_PI * l * y);
                }
            }
        }
      else if (btype[0] == eFourier && btype[1] != eFourier &&
	       stype == eQuadrilateral)
        {
	  for (int k = 0; k < order[0] / 2; ++k)
            {
	      for (int l = 0; l < order[1]; ++l)
                {
		  sol += M_PI * k * cos(M_PI * k * x) * pow_loc(y, l)
		    + l * sin(M_PI * k * x) * pow_loc(y, l - 1) +
		    -M_PI * k * sin(M_PI * k * x) * pow_loc(y, l)
		    + l * sin(M_PI * k * x) * pow_loc(y, l - 1);
                }
            }
        }
      else if (btype[0] == eFourier && btype[1] == eFourier &&
	       stype == eQuadrilateral)
        {
	  for (int k = 0; k < order[0] / 2; ++k)
            {
	      for (int l = 0; l < order[1] / 2; ++l)
                {
		  sol += M_PI * k * cos(M_PI * k * x) * sin(M_PI * l * y)
		    + M_PI * l * sin(M_PI * k * x) * cos(M_PI * l * y)
		    + M_PI * k * cos(M_PI * k * x) * cos(M_PI * l * y)
		    - M_PI * l * sin(M_PI * k * x) * sin(M_PI * l * y)
		    - M_PI * k * sin(M_PI * k * x) * sin(M_PI * l * y)
		    + M_PI * l * cos(M_PI * k * x) * cos(M_PI * l * y)
		    - M_PI * k * sin(M_PI * k * x) * cos(M_PI * l * y)
		    - M_PI * l * cos(M_PI * k * x) * sin(M_PI * l * y);
                }
            }
        }
      else
        {
	  NekDouble a;
	  for (int k = 0;
	       k < order[0]; ++k) //ShapeConstraint 1 is always < order1
            {
	      for (int l = 0; l < shapeConstraint2[stype](k, order); ++l)
                {
		  for (int m = 0;
		       m < shapeConstraint3[stype](k, l, order); ++m)
                    {
		      a = k * pow_loc(x, k - 1) * pow_loc(y, l) *
			pow_loc(z, m);
		      sol += a;
		      a = l * pow_loc(x, k) * pow_loc(y, l - 1) *
			pow_loc(z, m);
		      sol += a;
		      a = m * pow_loc(x, k) * pow_loc(y, l) *
			pow_loc(z, m - 1);
		      sol += a;
                    }
                }
            }
        }
    }

  return sol;
}
