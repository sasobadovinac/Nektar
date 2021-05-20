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
// Description: Demo for testing functionality of StdProject
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "StdDemoSupport.hpp"
#include <LibUtilities/BasicUtils/Timer.h>

namespace po = boost::program_options;

//declare Do_optimize
void Do_optimize(Array<OneD, NekDouble> &uhats);

//surface roots
void surfuhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz, StdExpansion* Etemp );

void call_setup_quad();
void call_setup_hex();
void call_setup_tet();
void call_setup_pyr();
void call_setup_tri();  
// declare caller routine to find_roots
Array<OneD, NekDouble> call_find_roots(Array<OneD, NekDouble> &uhatsall, NekDouble &avgiterGD, Array< OneD, Array<OneD, NekDouble> >&uhatsedges, Array< OneD, Array<OneD, NekDouble> >&surfaceuhats,  NekDouble &minv, NekDouble &roots1dtime,  NekDouble &roots2dtime , NekDouble &roots3dtime );

// for 2D elements to get uhats at edges
void project_edges( Array<OneD, NekDouble> uhats,Array<OneD, Array<OneD, NekDouble> >&ret);

int Opt_needed(Array<OneD, NekDouble> uhats, int flag = 0);
  
// derivative of scaled function stuff
void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret, int modes,   Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1 = NullNekDouble1DArray , Array<OneD, NekDouble> Vxyd2 = NullNekDouble1DArray);

void project_surfaces( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret);

string funcdef;
int iterstaken;
Array<OneD, Array<OneD, NekDouble> > storage3d;
Array<OneD, Array<OneD, NekDouble> > storage2dt;
Array<OneD, Array<OneD, NekDouble> > storage2dq;

// Confederate matrix
int numedges, numsurfaces;
StdExpansion *E, *Eorth;
//StdExpansion *E3seg;
StdExpansion *Equad;
StdExpansion *Etri;

Array<OneD, Array<OneD, NekDouble> > edgeptsin;

// for edges root finding hex:
// edge front left (x = -1) (y = -1)
Array<OneD, NekDouble> Vxm1ym1z   ;
Array<OneD, NekDouble> Vdxxm1ym1z   ;
Array<OneD, NekDouble> Vdyxm1ym1z   ;
Array<OneD, NekDouble> Vdzxm1ym1z   ;
                
//edge front right (x = 1) (y = -1)
Array<OneD, NekDouble> Vx1ym1z   ;
Array<OneD, NekDouble> Vdxx1ym1z   ;
Array<OneD, NekDouble> Vdyx1ym1z   ;
Array<OneD, NekDouble> Vdzx1ym1z   ;
               
//edge front top (y = -1) (z = 1)
Array<OneD, NekDouble> Vym1xz1   ;
Array<OneD, NekDouble> Vdxym1xz1   ;
Array<OneD, NekDouble> Vdyym1xz1   ;
Array<OneD, NekDouble> Vdzym1xz1   ;
               
//edge front bot (y = -1) (z = -1)
Array<OneD, NekDouble> Vym1xzm1   ;
Array<OneD, NekDouble> Vdxym1xzm1   ;
Array<OneD, NekDouble> Vdyym1xzm1   ;
Array<OneD, NekDouble> Vdzym1xzm1   ;

// edge back left (y = 1), (x = -1)
Array<OneD, NekDouble> Vxm1y1z   ;
Array<OneD, NekDouble> Vdxxm1y1z   ;
Array<OneD, NekDouble> Vdyxm1y1z   ;
Array<OneD, NekDouble> Vdzxm1y1z   ;
                
//edge back right (x = 1), (y = 1))
Array<OneD, NekDouble> Vx1y1z   ;
Array<OneD, NekDouble> Vdxx1y1z   ;
Array<OneD, NekDouble> Vdyx1y1z   ;
Array<OneD, NekDouble> Vdzx1y1z   ;
               
//edge back top ( y = 1) (z = 1)
Array<OneD, NekDouble> Vy1xz1   ;
Array<OneD, NekDouble> Vdxy1xz1   ;
Array<OneD, NekDouble> Vdyy1xz1   ;
Array<OneD, NekDouble> Vdzy1xz1   ;
               
//edge back bot (y = 1) (z = -1))
Array<OneD, NekDouble> Vy1xzm1   ;
Array<OneD, NekDouble> Vdxy1xzm1   ;
Array<OneD, NekDouble> Vdyy1xzm1   ;
Array<OneD, NekDouble> Vdzy1xzm1   ;

// edge left bot (z = -1), (x = -1)
Array<OneD, NekDouble> Vxm1yzm1   ;
Array<OneD, NekDouble> Vdxxm1yzm1   ;
Array<OneD, NekDouble> Vdyxm1yzm1   ;
Array<OneD, NekDouble> Vdzxm1yzm1   ;
                
//edge left top (x = -1), (z = 1))
Array<OneD, NekDouble> Vxm1yz1   ;
Array<OneD, NekDouble> Vdxxm1yz1   ;
Array<OneD, NekDouble> Vdyxm1yz1   ;
Array<OneD, NekDouble> Vdzxm1yz1   ;
               
//edge right bot ( z = -1) (x = 1)
Array<OneD, NekDouble> Vx1yzm1   ;
Array<OneD, NekDouble> Vdxx1yzm1   ;
Array<OneD, NekDouble> Vdyx1yzm1   ;
Array<OneD, NekDouble> Vdzx1yzm1   ;
               
//edge right top (z  1) (x  1))
Array<OneD, NekDouble>Vx1yz1 ;
Array<OneD, NekDouble>Vdxx1yz1 ;
Array<OneD, NekDouble>Vdyx1yz1 ;
Array<OneD, NekDouble> Vdzx1yz1;

//hex edges:

//surface bot z = -1
Array<OneD, NekDouble>Vxyzm1 ;

//surface right x = 1
Array<OneD, NekDouble>Vx1yz ;

//surface top z = 1
Array<OneD, NekDouble>Vxyz1 ; 

//surface left x = -1
Array<OneD, NekDouble>Vxm1yz ;

//surface front y = -1
Array<OneD, NekDouble>Vxym1z ;

//surface back y = 1
Array<OneD, NekDouble>Vxy1z ;

// only tets edges:

// edge front left (x = -1) (y = -1)
Array<OneD, NekDouble> Vxm1ym1ztet   ;
Array<OneD, NekDouble> Vdxxm1ym1ztet   ;
Array<OneD, NekDouble> Vdyxm1ym1ztet   ;
Array<OneD, NekDouble> Vdzxm1ym1ztet   ;


//edge front bot (y = -1) (z = -1)
Array<OneD, NekDouble> Vym1xzm1tet   ;
Array<OneD, NekDouble> Vdxym1xzm1tet   ;
Array<OneD, NekDouble> Vdyym1xzm1tet   ;
Array<OneD, NekDouble> Vdzym1xzm1tet   ;


// edge left bot (z = -1), (x = -1)
Array<OneD, NekDouble> Vxm1yzm1tet   ;
Array<OneD, NekDouble> Vdxxm1yzm1tet   ;
Array<OneD, NekDouble> Vdyxm1yzm1tet   ;
Array<OneD, NekDouble> Vdzxm1yzm1tet   ;


//edge front hypt (DB) (y = -1) (z = -x)                                                       
Array<OneD, NekDouble> Vym1xmztet   ;
Array<OneD, NekDouble> Vdxym1xmztet   ;
Array<OneD, NekDouble> Vdyym1xmztet   ;
Array<OneD, NekDouble> Vdzym1xmztet   ;

//edge left hypt (DC) ( x = -1) (z = -y)                                                       
Array<OneD, NekDouble> Vxm1ymztet   ;
Array<OneD, NekDouble> Vdxxm1ymztet   ;
Array<OneD, NekDouble> Vdyxm1ymztet   ;
Array<OneD, NekDouble> Vdzxm1ymztet   ;

// edge bot diag (BC) (z = -1) (y = -x)                                                        
Array<OneD, NekDouble> Vxmyzm1tet   ;
Array<OneD, NekDouble> Vdxxmyzm1tet   ;
Array<OneD, NekDouble> Vdyxmyzm1tet   ;
Array<OneD, NekDouble> Vdzxmyzm1tet   ;

//tet surfaces:
//surface bot z = -1 restriction on GD: x+y = 0 (ABC)
Array<OneD, NekDouble>  Vxyzm1tet ;

//surface left x = -1 restriction on GD: y + z = 0 (DAC)
Array<OneD, NekDouble>Vxm1ypz0tet;

//surf front y = -1   restriction on GD: x + z = 0 (DAB)
Array<OneD, NekDouble> Vxpz0ym1tet;

//surf DCB restriction  on GD: (x + y + z = -1)  
Array<OneD, NekDouble> Vxpypzm1tet;

//pyr dges:
// edge front left (x = -1) (y = -1) EA
Array<OneD, NekDouble> Vxm1ym1zpyr;
Array<OneD, NekDouble> Vdyxm1ym1zpyr;
Array<OneD, NekDouble> Vdxxm1ym1zpyr;
Array<OneD, NekDouble> Vdzxm1ym1zpyr;

//edge front hypt (y = -1) (z = -x) EB
Array<OneD, NekDouble> Vym1xmzpyr   ;
Array<OneD, NekDouble> Vdxym1xmzpyr   ;
Array<OneD, NekDouble> Vdyym1xmzpyr   ;
Array<OneD, NekDouble> Vdzym1xmzpyr   ;


//edge back hypt (z = -x  or x = -z) (x = y) EC
Array<OneD, NekDouble> Vxeyxmzpyr   ;
Array<OneD, NekDouble> Vdxxeyxmzpyr   ;
Array<OneD, NekDouble> Vdyxeyxmzpyr   ;
Array<OneD, NekDouble> Vdzxeyxmzpyr   ;

//edge back left (y = -z) (x = -1) ED
Array<OneD, NekDouble> Vx1ymzpyr   ;
Array<OneD, NekDouble> Vdxx1ymzpyr   ;
Array<OneD, NekDouble> Vdyx1ymzpyr   ;
Array<OneD, NekDouble> Vdzx1ymzpyr   ;


//edge front bot (y = -1) (z = -1) AB
Array<OneD, NekDouble> Vym1xzm1pyr   ;
Array<OneD, NekDouble> Vdyym1xzm1pyr   ;
Array<OneD, NekDouble> Vdxym1xzm1pyr   ;
Array<OneD, NekDouble> Vdzym1xzm1pyr   ;


//edge back bot (y = 1) (z = -1)) DC
Array<OneD, NekDouble> Vy1xzm1pyr   ;
Array<OneD, NekDouble> Vdyy1xzm1pyr   ;
Array<OneD, NekDouble> Vdxy1xzm1pyr   ;
Array<OneD, NekDouble> Vdzy1xzm1pyr   ;

// edge left bot (z = -1), (x = -1) AD
Array<OneD, NekDouble> Vxm1yzm1pyr   ;
Array<OneD, NekDouble> Vdyxm1yzm1pyr   ;
Array<OneD, NekDouble> Vdxxm1yzm1pyr   ;
Array<OneD, NekDouble> Vdzxm1yzm1pyr   ;

//edge right bot ( z = -1) (x = 1) BC
Array<OneD, NekDouble> Vx1yzm1pyr   ;
Array<OneD, NekDouble> Vdyx1yzm1pyr   ;
Array<OneD, NekDouble> Vdxx1yzm1pyr   ;
Array<OneD, NekDouble> Vdzx1yzm1pyr   ;

//pyr surfaces
//surface bot z = -1
Array<OneD, NekDouble> Vxyzm1pyr ;

//surface hypt x+z = 0
Array<OneD, NekDouble> Vxmzypyr;

//surface left x = -1
Array<OneD, NekDouble> Vxm1yzpyr;


//surface front y = -1
Array<OneD, NekDouble> Vxym1zpyr;

//surface back y +z = 0
Array<OneD, NekDouble> Vxymzpyr;



int dimension ;
NekDouble  startval, startcoordx, startcoordy, startcoordz, itersGD1, itersGD2, itersGD3;

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
  storage3d = E->GetPhysEvaluateStorage();
  if(dimension < 3 )
    {
      cout<<"\n dimension should be 3, try using StdProjectPositivityPres1D or StdProjectPositivityPres2D for other dimensions\n\n";
      exit(0);
    }

  std::vector<int> order;
  std::vector<BasisType> btype(3, eNoBasisType);
  LibUtilities::ShapeType stype = E->DetShapeType();

  for (int i = 0; i < dimension; ++i)
    {
      btype[i] = E->GetBasisType(i);
      order.push_back(E->GetBasisNumModes(i));
    }
  if(dimension > 2)
    {
      demo.testcoord3dhqpts  =  demo.GetCoords(E);
      demo.testcoord3dhqmidpts  = demo.GetQuadratureMidCoords(E, demo.testcoord3dhqpts);
      demo.testcoord3dhlattice = demo.GetLatticeCoords(demo.testcoord3dhqpts, demo.testcoord3dhqmidpts);
      demo.interioreval3hqmidpts = E->PhysEvaluateBasis(demo.testcoord3dhqmidpts, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray   );
    }
  switch(E->DetShapeType())
    {
    case LibUtilities::eQuadrilateral:
      break;
    case LibUtilities::eTriangle:
      break;
    
    case LibUtilities::eTetrahedron:
      {
      numedges = 6;
      numsurfaces = 4;

      }
      break;

    case LibUtilities::eHexahedron:
        
      numedges = 12;
      numsurfaces = 6;

      break;

    case LibUtilities::ePyramid:
        
      numedges = 8;
      numsurfaces = 5;
      
      break;

      default: cout<<"\n unknown shape typen\n";exit(0);
    
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
  for (int i = 0; i < x.size(); ++i)
    {
      if(dimension == 3)
        {
	  sol[i] =   (pow(x[i] +0.3, 2) + pow(y[i]+0.290413,2)) + pow(z[i]+0.419174,2) -0.1;//pow((x[i]+0.6),2)+ pow((y[i]+0.6),2) + pow((z[i]-0.3),2) + demo.GetEps();
	  funcdef ="pow((x[i]+0.2),2)+ pow((y[i]-0.3),2) + pow((z[i]-0.5),2) + demo.GetEps()";
        }
    }
    
  Array<OneD, NekDouble> phys(totPoints);
  Array<OneD, NekDouble> coeffs((unsigned) E->GetNcoeffs());
    
  //Project onto expansion
  E->FwdTrans(sol, coeffs);

  //Backward transform solution to get projected values
  E->BwdTrans(coeffs, phys);

  LibUtilities::Timer     timer;

  Array<OneD, NekDouble> physold(phys.size());
  Vmath::Vcopy(phys.size(), phys, 1, physold, 1);

  // check for -ve values and apply opt if necessary
  if (vm.count("optm"))
    {

      timer.Start();

      int opt_need = Opt_needed(coeffs);            
      timer.Stop();
      NekDouble  timeneg = timer.TimePerTest(1);  
      cout<<"\nchecking for -vity took "<<timeneg<<"s\n";

      Eorth = E;      
      if(opt_need)
	{
	  NekDouble setup = 0, timeverify;
	  timer.Start();
	  int orthoflag = (E->GetBasisType(0) != LibUtilities::eOrtho_A);
	  Array<OneD, BasisType> btypearr(dimension);
	  btypearr[0] = LibUtilities::eOrtho_A; 

	  timer.Start();
	  if(orthoflag)
	    {
	      switch(E->DetShapeType())
		{
		case LibUtilities::eHexahedron:
		  {
		    	LibUtilities::BasisKey bkeyorthA ( LibUtilities::eOrtho_A,
							   E->GetBasis(0)->GetNumModes(),
							   E->GetBasis(0)->GetPointsKey());
			LibUtilities::BasisKey bkeyorthB ( LibUtilities::eOrtho_A,
							   E->GetBasis(1)->GetNumModes(),
							   E->GetBasis(1)->GetPointsKey());
			LibUtilities::BasisKey bkeyorthC ( LibUtilities::eOrtho_A,
							   E->GetBasis(2)->GetNumModes(),
							   E->GetBasis(2)->GetPointsKey());
	  
			demo.Eorthele =  new StdHexExp(bkeyorthA, bkeyorthB, bkeyorthC);
			demo.interiorevalorth3hqmidpts = demo.Eorthele->PhysEvaluateBasis(demo.testcoord3dhqmidpts , demo.Eorthele->GetPhysEvaluateStorage(), NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
		      

		    if( demo.Eorthele != nullptr)  
		    {
		      Eorth = demo.Eorthele;
		      btypearr[1] = LibUtilities::eOrtho_A;
		      btypearr[2] = LibUtilities::eOrtho_A;  
		      
		    }
		  }
		    break;
		case LibUtilities::eTetrahedron: 
		  {
		          // check if we need to activate the ortho version of this element
		    LibUtilities::BasisKey bkeyorthA ( LibUtilities::eOrtho_A,
						       E->GetBasis(0)->GetNumModes(),
						       E->GetBasis(0)->GetPointsKey());
		    LibUtilities::BasisKey bkeyorthB ( LibUtilities::eOrtho_B,
						       E->GetBasis(1)->GetNumModes(),
						       E->GetBasis(1)->GetPointsKey());
		    LibUtilities::BasisKey bkeyorthC ( LibUtilities::eOrtho_C,
						       E->GetBasis(2)->GetNumModes(),
						       E->GetBasis(2)->GetPointsKey());
	  
		    demo.Eorthele =  new StdTetExp(bkeyorthA, bkeyorthB, bkeyorthC);
		    demo.interiorevalorth3hqmidpts = demo.Eorthele->PhysEvaluateBasis(demo.testcoord3dhqmidpts , demo.Eorthele->GetPhysEvaluateStorage(), NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
		  

		    if( demo.Eorthele != nullptr)
		      {
			Eorth = demo.Eorthele;   
			btypearr[1] = LibUtilities::eOrtho_B;
			btypearr[2] = LibUtilities::eOrtho_C;
		      }
		  }
		  break;
		case LibUtilities::ePyramid:
		  {
		       // check if we need to activate the ortho version of this element
		    LibUtilities::BasisKey bkeyorthA ( LibUtilities::eOrtho_A,
						       E->GetBasis(0)->GetNumModes(),
						       E->GetBasis(0)->GetPointsKey());
		    LibUtilities::BasisKey bkeyorthB ( LibUtilities::eOrtho_A,
						       E->GetBasis(1)->GetNumModes(),
						       E->GetBasis(1)->GetPointsKey());
		    LibUtilities::BasisKey bkeyorthC ( LibUtilities::eOrthoPyr_C,
						       E->GetBasis(2)->GetNumModes(),
						       E->GetBasis(2)->GetPointsKey());
	  
		    demo.Eorthele =  new StdPyrExp(bkeyorthA, bkeyorthB, bkeyorthC);
		    demo.interiorevalorth3hqmidpts = demo.Eorthele->PhysEvaluateBasis(demo.testcoord3dhqmidpts , demo.Eorthele->GetPhysEvaluateStorage(), NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
		  

		    if( demo.Eorthele != nullptr)
		      {
			Eorth = demo.Eorthele;   
			btypearr[1] = LibUtilities::eOrtho_B;
			btypearr[2] = LibUtilities::eOrtho_C;
		      }
		  }
		  break;
	        default:
		  cout<<"\n invalid element type!\n";
		  exit(0);
		}
	    }
	
	  storage3d = Eorth->GetPhysEvaluateStorage();
	  
	  if(stype != LibUtilities::eTetrahedron)
	    {
	      call_setup_quad();
	    }
	  if(stype != LibUtilities::eHexahedron)
	    {
	      call_setup_tri();
	    }
      

	  if(numedges == 12) //hex
	    {
	      call_setup_hex();
	    }
	  else if(numedges == 6) //tet
	    {
	      call_setup_tet();
	    }
	  else if(numedges == 8) //pyr
	    {
	      call_setup_pyr();
	    }
	  timer.Stop();
	  setup += timer.TimePerTest(1);

	  cout<<"\nsetup phase took "<<setup<<"s";
	  if(orthoflag)
	    {
	      E->BwdTrans(coeffs, phys);
	      Eorth->FwdTrans(phys, coeffs);
	    }
	  
	  NekDouble elap = 0;
	  timer.Start();

	  Do_optimize(coeffs);

	  timer.Stop();
	  elap += timer.TimePerTest(1); 
	  cout<<"\noptimizertook "<<elap<<"s";

	  if(orthoflag)
	    {
	      Eorth->BwdTrans(coeffs, phys);
	      E->FwdTrans(phys, coeffs);
	      storage3d = E->GetPhysEvaluateStorage();
	    }
	  timer.Start();
	  opt_need = Opt_needed(coeffs, 1 );
	  // check val now at the min ppreviously found:
	  Array<OneD, Array< OneD, NekDouble> > checkcoord(dimension);
	  Array<OneD, NekDouble> saveeval, holdv(coeffs.size());
	  checkcoord[0] = Array<OneD, NekDouble>(1, startcoordx);
	  checkcoord[1] = Array<OneD, NekDouble>(1, startcoordy);
	  checkcoord[2] = Array<OneD, NekDouble>(1, startcoordz);
	  saveeval = E->PhysEvaluateBasis(checkcoord, storage3d,NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
	  Vmath::Vmul(coeffs.size(), saveeval, 1, coeffs, 1, holdv, 1);
	  NekDouble valnow =  Vmath::Vsum(coeffs.size(), holdv, 1);
	  cout<<"\n val now at "<<startcoordx<<" "<<startcoordy<<" "<<startcoordz<<" ="<<valnow;
	  if(valnow < 0 && abs(valnow)>1e-9)
	    {
	      opt_need = 1;
	    }
	  timer.Stop();
	  timeverify = timer.TimePerTest(1);
	  cout<<"\nVerification took "<<timeverify<<"s";
	
	  if(opt_need)
	    {
	  
	      //cout<<"func="<<funcdef;
	  
	      cout<<"\n\n fail\n";
	  
	      //cout<<" "<<" "<<startval<<", ("<<startcoordx<<" "<<startcoordy<<" "<<startcoordz<<"), "<<setup<<", fail, "<<" "<<itersGD1<<", "<<itersGD2<<", "<<itersGD3<<", "<<timeneg<<", "<<elap<<", "<<timeverify<<","<<elapref<<" \n";;
	  
	    }
	  else
	    {
	      //cout<<"func="<<funcdef;
	      cout<<"\n\n pass\n\n";
	  
	      //	cout<<" "<<" "<<startval<<", ("<<startcoordx<<" "<<startcoordy<<" "<<startcoordz<<"), "<<setup<<", pass, "<<" "<<itersGD1<<", "<<itersGD2<<", "<<itersGD3<<", "<<iterstaken<<", "<<timeneg<<", "<<elap<<", "<<timeverify<<", "<<elapref<<"\n";;
	      //WriteCSV(coeffs, 1);
	  
	      sol = phys; 
	    }
          
	}
      else{
    
	cout<<"\noptimizer no need\n\n";
      }
    }
  //Backward transform solution to get projected values
  E->BwdTrans(coeffs, phys);
  

  
  //Calculate L_inf & L_2 error
  cout << "L infinity error: \t" << E->Linf(phys, sol) << endl;
  cout << "L 2 error: \t \t \t" << E->L2(phys, sol) << endl;
  cout<<"\n***********\n";
  return 0;
    
}

  void call_setup_tri()
  {
    // LibUtilities::PointsType pointsTypetriA = LibUtilities::eGaussLobattoLegendre;
    // LibUtilities::PointsKey pkeytri(E->GetBasis(0)->GetZ().size(), pointsTypetriA);
    LibUtilities::BasisKey bkeytriB(LibUtilities::eOrtho_B,(E->GetBasis(0)->GetNumModes()), E->GetBasis(0)->GetPointsKey()); 
    LibUtilities::BasisKey bkeytriA(LibUtilities::eOrtho_A, (E->GetBasis(1)->GetNumModes()),E->GetBasis(1)->GetPointsKey());
    Etri = new StdTriExp(bkeytriA, bkeytriB);
    if (Etri== nullptr)
      {
	cout<<"\n error!\n";
	exit(0);
      }
    storage2dt = Etri->GetPhysEvaluateStorage();
    // to-do : move this in demo constructor
    demo.testcoord2dtqpts = demo.GetCoords(Etri);
    demo.testcoord2dtqmidpts = demo.GetQuadratureMidCoords( Etri, demo.testcoord2dtqpts);
    demo.testcoord2dtlattice = demo.GetLatticeCoords(demo.testcoord2dtqpts, demo.testcoord2dtqmidpts);
    demo.interioreval2dtqmidpts = Etri->PhysEvaluateBasis(  demo.testcoord2dtqmidpts, storage2dt, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray   );

  }

  void call_setup_quad()
  {
    LibUtilities::BasisKey bkeyA ( LibUtilities::eOrtho_A, (E->GetBasis(0)->GetNumModes()),  E->GetBasis(0)->GetPointsKey());
    Equad = new StdQuadExp(bkeyA, bkeyA);
    if (Equad == nullptr)
      {
	cout<<"\n error!\n";
	exit(0);
      }
    storage2dq = Equad->GetPhysEvaluateStorage();
    // to-do : move this in demo constructor
    demo.testcoord2dqqpts = demo.GetCoords(Equad);
    demo.testcoord2dqqmidpts = demo.GetQuadratureMidCoords(Equad,  demo.testcoord2dqqpts);
    demo.testcoord2dqlattice = demo.GetLatticeCoords(  demo.testcoord2dqqpts, demo.testcoord2dqqmidpts);
  
    demo.interioreval2dqqmidpts = Equad->PhysEvaluateBasis(  demo.testcoord2dqqmidpts, storage2dq, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray );
  
}

void call_setup_tet()
{
  Array<OneD, NekDouble> edgexyztemp =   demo.E3seg->GetBasis(0)->GetZ();

  int totszedges1d =  edgexyztemp.size()*(E->GetNcoeffs());
  edgeptsin = Array<OneD, Array<OneD, NekDouble> >(dimension);
  for(int p = 0; p < dimension; p++)
    {
      edgeptsin[p] =   Array<OneD, NekDouble>(edgexyztemp);
    }

  // edge front left (AD) (x = -1) (y = -1)                                                
  Vxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

  Vxm1ym1ztet = E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxxm1ym1ztet, Vdyxm1ym1ztet, Vdzxm1ym1ztet);

  //edge front hypt (DB) (y = -1) (z = -x)                                                 
  Vym1xmztet = Array<OneD, NekDouble>(totszedges1d);
  Vdxym1xmztet = Array<OneD, NekDouble>(totszedges1d);
  Vdyym1xmztet = Array<OneD, NekDouble>(totszedges1d);
  Vdzym1xmztet = Array<OneD, NekDouble>(totszedges1d);

  Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[2][0], 1, &edgeptsin[0][0], 1);

  Vym1xmztet = E->PhysEvaluateBasis(edgeptsin,storage3d, Vdxym1xmztet, Vdyym1xmztet, Vdzym1xmztet);

  //edge front bot (AB) (y = -1) (z = -1)                                                
  Vym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
  Vdxym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
  Vdyym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
  Vdzym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = edgexyztemp;
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

  Vym1xzm1tet = E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxym1xzm1tet, Vdyym1xzm1tet, Vdzym1xzm1tet);

  //edge left hypt (DC) ( x = -1) (z = -y)                                               
  Vxm1ymztet = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1ymztet= Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1ymztet = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1ymztet= Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[1] = edgexyztemp;
  Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[1][0], 1, &edgeptsin[2][0], 1);
  
  Vxm1ymztet = E->PhysEvaluateBasis(edgeptsin,storage3d, Vdxxm1ymztet, Vdyxm1ymztet, Vdzxm1ymztet);

  // edge bot diag (BC) (z = -1) (y = -x)                                                
  Vxmyzm1tet = Array<OneD, NekDouble>(totszedges1d);
  Vdxxmyzm1tet = Array<OneD, NekDouble>(totszedges1d)   ;
  Vdyxmyzm1tet = Array<OneD, NekDouble>(totszedges1d)    ;
  Vdzxmyzm1tet = Array<OneD, NekDouble>(totszedges1d)    ;
  edgeptsin[0] = edgexyztemp;
  Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[0][0], 1, &edgeptsin[1][0], 1);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[2].size(), -1.0);
  Vxmyzm1tet = E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxxmyzm1tet, Vdyxmyzm1tet, Vdzxmyzm1tet);
  
  //edge CA bot left (x = -1) (z = -1)                                                   
  Vxm1yzm1tet   = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1yzm1tet = Array<OneD, NekDouble>(totszedges1d)   ;
  Vdyxm1yzm1tet = Array<OneD, NekDouble>(totszedges1d)    ;
  Vdzxm1yzm1tet = Array<OneD, NekDouble>(totszedges1d)    ;
  edgeptsin[1] = edgexyztemp;
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[2].size(), -1.0);
  
  Vxm1yzm1tet = E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxxm1yzm1tet, Vdyxm1yzm1tet, Vdzxm1yzm1tet);
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
  Vxyzm1tet = Array<OneD, NekDouble>(totsurf2d);
  surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[1] = surfptsin[1];
  Vxyzm1tet = E->PhysEvaluateBasis(surfptsintemp, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  //surface left x = -1  (DAC)  
  Vxm1ypz0tet =  Array<OneD, NekDouble>(totsurf2d);
  surfptsintemp[0] = Array<OneD, NekDouble> (totpt, -1.0);
  surfptsintemp[1] = surfptsin[0];
  surfptsintemp[2] = surfptsin[1];
  
  Vxm1ypz0tet = E->PhysEvaluateBasis(surfptsintemp, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  //surf front y = -1,  (DAB) 
  Vxpz0ym1tet =  Array<OneD, NekDouble>(totsurf2d);
  surfptsintemp[1] = Array<OneD, NekDouble> (totpt, -1.0);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[2] = surfptsin[1];
  Vxpz0ym1tet = E->PhysEvaluateBasis(surfptsintemp, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  //surf DCB (x + y + z = -1), 
  Vxpypzm1tet =  Array<OneD, NekDouble>(totsurf2d);
  surfptsintemp[1] = surfptsin[1];;

  Vmath::Vadd(totpt, &surfptsin[0][0], 1, &surfptsin[1][0], 1, &surfptsintemp[2][0], 1);
  Vmath::Smul(totpt, -1.0, &surfptsintemp[2][0], 1, &surfptsintemp[2][0], 1);
  Vmath::Sadd(totpt, -1.0, &surfptsintemp[2][0], 1, &surfptsintemp[2][0], 1);
  Vxpypzm1tet = E->PhysEvaluateBasis(surfptsintemp,storage3d , NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
  
}

void call_setup_hex()
{
  Array<OneD, NekDouble> edgexyztemp =   demo.E3seg->GetBasis(0)->GetZ();

  int totszedges1d = edgexyztemp.size()*(Eorth->GetNcoeffs());
  edgeptsin = Array<OneD, Array<OneD, NekDouble> >(dimension);    
  for(int p = 0; p < dimension; p++)
    {
      edgeptsin[p] = Array<OneD, NekDouble>(edgexyztemp);
    } 

  // edge front left (x = -1) (y = -1)
  Vxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

            
  Vxm1ym1z = Eorth->PhysEvaluateBasis(edgeptsin, storage3d, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);  
                
  //edge front right (x = 1) (y = -1)
  Vx1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1ym1z = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);       
  Vx1ym1z = Eorth->PhysEvaluateBasis(edgeptsin, storage3d, Vdxx1ym1z, Vdyx1ym1z, Vdzx1ym1z);  


  //edge front top (y = -1) (z = 1)
  Vym1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxym1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyym1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzym1xz1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = edgexyztemp; 
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
  Vym1xz1 = Eorth->PhysEvaluateBasis(edgeptsin,storage3d, Vdxym1xz1, Vdyym1xz1, Vdzym1xz1);  

                
  //edge front bot (y = -1) (z = -1)
  Vym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  Vym1xzm1 =             Eorth->PhysEvaluateBasis( edgeptsin, storage3d, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);  
	    
	    
  // edge back left (y = 1), (x = -1)
  Vxm1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1y1z = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);      
  edgeptsin[2] =    edgexyztemp;


  Vxm1y1z = Eorth->PhysEvaluateBasis(edgeptsin,storage3d, Vdxxm1y1z, Vdyxm1y1z, Vdzxm1y1z);  

  //edge back right (x = 1), (y = 1))
  Vx1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1y1z = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);                
         
  Vx1y1z  = Eorth->PhysEvaluateBasis(edgeptsin,storage3d, Vdxx1y1z, Vdyx1y1z, Vdzx1y1z);  

                
  //edge back top ( y = 1) (z = 1)
  Vy1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxy1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyy1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzy1xz1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = edgexyztemp;
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
  Vy1xz1 =             Eorth->PhysEvaluateBasis(edgeptsin, storage3d, Vdxy1xz1, Vdyy1xz1, Vdzy1xz1 );  

      
  //edge back bot (y = 1) (z = -1))
  Vy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
            
  Vy1xzm1 = Eorth->PhysEvaluateBasis(edgeptsin, storage3d, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);  

  // edge left bot (z = -1), (x = -1)
  Vxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[1] = edgexyztemp;
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

  Vxm1yzm1 = Eorth->PhysEvaluateBasis(edgeptsin,storage3d, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);  
                
  //edge left top (x = -1), (z = 1))
  Vxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
                
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);           
  Vxm1yz1 = Eorth->PhysEvaluateBasis(edgeptsin,storage3d, Vdxxm1yz1, Vdyxm1yz1, Vdzxm1yz1);  

  //edge right bot ( z = -1) (x = 1)
  Vx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1yzm1 = Array<OneD, NekDouble>(totszedges1d);

  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);           
  Vx1yzm1 = Eorth->PhysEvaluateBasis(edgeptsin, storage3d, Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);   

  //edge right top (z  1) (x  1))
  Vx1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1yz1 = Array<OneD, NekDouble>(totszedges1d);
                
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);

  Vx1yz1 = Eorth->PhysEvaluateBasis(edgeptsin,storage3d, Vdxx1yz1, Vdyx1yz1, Vdzx1yz1);   

              
  int totszsurf2d = (demo.testcoord2dqqpts[0].size())*Eorth->GetNcoeffs();
            
  Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);
  for(int k = 0; k < dimension-1; k++)
    {
      surfptsin[k] = Array<OneD, NekDouble>(demo.testcoord2dqqpts[k]);
      surfptsintemp[k] = Array<OneD, NekDouble>(demo.testcoord2dqqpts[k]);
    }
            
  surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.testcoord2dqqpts[0].size()); 
  surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.testcoord2dqqpts[0].size()); 
            
  //surface bot z = -1
  Vxyzm1 =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[1] = surfptsin[1];
  Vxyzm1 =             Eorth->PhysEvaluateBasis(surfptsintemp,storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

  //surface right x = 1
  Vx1yz =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[0] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
  surfptsintemp[1] =  surfptsin[0];
  surfptsintemp[2] =  surfptsin[1];
  Vx1yz = Eorth->PhysEvaluateBasis(surfptsintemp,storage3d,NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   


  //surface top z = 1
  Vxyz1 =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[1] = surfptsin[1];
  surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
  Vxyz1 = Eorth->PhysEvaluateBasis(surfptsintemp, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

  //surface left x = -1
  Vxm1yz =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[0] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
  surfptsintemp[1] = surfptsin[0];
  surfptsintemp[2] = surfptsin[1];
  Vxm1yz = Eorth->PhysEvaluateBasis(surfptsintemp,storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

  //surface front y = -1
  Vxym1z =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[1] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[2] = surfptsin[1];
  Vxym1z = Eorth->PhysEvaluateBasis(surfptsintemp,storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   
      
  //surface back y = 1
  Vxy1z =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[1] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[2] = surfptsin[1];
  Vxy1z = Eorth->PhysEvaluateBasis(surfptsintemp,storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
}

void call_setup_pyr()
{

  //we need qzin + qzinmid as input points                                                         
  Array<OneD, NekDouble> edgexyztemp =   demo.E3seg->GetBasis(0)->GetZ();

  int totszedges1d =  edgexyztemp.size()*(E->GetNcoeffs());      
  edgeptsin = Array<OneD, Array<OneD, NekDouble> >(dimension);
  for(int p = 0; p < dimension; p++)
    {
      edgeptsin[p] = Array<OneD, NekDouble>(edgexyztemp);
    }
		
  // edge front left EA (x = -1) (y = -1)                 
  Vxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);

  edgeptsin[0]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[1]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  Vxm1ym1zpyr   = E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxxm1ym1zpyr, Vdyxm1ym1zpyr, Vdzxm1ym1zpyr);

  //edge front hypt EB (y = -1) (z + x = 0)
  Vym1xmzpyr    = Array<OneD, NekDouble>(totszedges1d);
  Vdxym1xmzpyr  = Array<OneD, NekDouble>(totszedges1d);
  Vdyym1xmzpyr  = Array<OneD, NekDouble>(totszedges1d);
  Vdzym1xmzpyr  = Array<OneD, NekDouble>(totszedges1d);
  Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[2][0], 1, &edgeptsin[0][0], 1);

  Vxm1ym1zpyr =  E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxym1xmzpyr, Vdxym1xmzpyr, Vdzym1xmzpyr);
  //edge back hypt (z = -x  or x = -z) (x = y) EC
  Vxeyxmzpyr    = Array<OneD, NekDouble>(totszedges1d);
  Vdxxeyxmzpyr  = Array<OneD, NekDouble>(totszedges1d);
  Vdyxeyxmzpyr  = Array<OneD, NekDouble>(totszedges1d);
  Vdzxeyxmzpyr  = Array<OneD, NekDouble>(totszedges1d);

  edgeptsin[1]  =  Array<OneD, NekDouble>(edgexyztemp);
  Vxeyxmzpyr    =  E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxxeyxmzpyr, Vdyxeyxmzpyr, Vdzxeyxmzpyr);

  //edge back left (y = -z) (x = -1) ED
  Vx1ymzpyr     = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1ymzpyr   = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1ymzpyr   = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1ymzpyr   = Array<OneD, NekDouble>(totszedges1d); ;
  Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[2][0], 1, &edgeptsin[1][0], 1);

  edgeptsin[0]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  Vx1ymzpyr     = E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxx1ymzpyr, Vdyx1ymzpyr, Vdzx1ymzpyr);
  
  //edge front bot (y = -1) (z = -1) AB
  Vym1xzm1pyr   = Array<OneD, NekDouble>(totszedges1d);
  Vdxym1xzm1pyr = Array<OneD, NekDouble>(totszedges1d)   ;
  Vdyym1xzm1pyr = Array<OneD, NekDouble>(totszedges1d)    ;
  Vdzym1xzm1pyr = Array<OneD, NekDouble>(totszedges1d)    ;

  edgeptsin[1]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[2]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[0]  = Array<OneD, NekDouble>(edgexyztemp);
  Vym1xzm1pyr   = E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxym1xzm1pyr, Vdyym1xzm1pyr,Vdzym1xzm1pyr);

  
  //edge back bot (y = 1) (z = -1)) DC
  Vy1xzm1pyr    = Array<OneD, NekDouble>(totszedges1d);
  Vdyy1xzm1pyr  = Array<OneD, NekDouble>(totszedges1d)   ;
  Vdxy1xzm1pyr  = Array<OneD, NekDouble>(totszedges1d);
  Vdzy1xzm1pyr  = Array<OneD, NekDouble>  (totszedges1d) ;

  edgeptsin[0]  = edgexyztemp;
  edgeptsin[1]  = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
  Vy1xzm1pyr   = E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxy1xzm1pyr, Vdyy1xzm1pyr, Vdzy1xzm1pyr);

  
  // edge left bot (z = -1), (x = -1) AD
  Vxm1yzm1pyr   = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1yzm1pyr = Array<OneD, NekDouble>(totszedges1d)   ;
  Vdxxm1yzm1pyr = Array<OneD, NekDouble>(totszedges1d)   ;
  Vdzxm1yzm1pyr = Array<OneD, NekDouble>(totszedges1d)   ;

  edgeptsin[1]  = edgexyztemp;
  edgeptsin[0]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[2]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  Vxm1yzm1pyr   =  E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxxm1yzm1pyr, Vdyxm1yzm1pyr, Vdzxm1yzm1pyr);

  
  //edge right bot ( z = -1) (x = 1) BC
  Vx1yzm1pyr    = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1yzm1pyr  = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1yzm1pyr  = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1yzm1pyr  = Array<OneD, NekDouble>(totszedges1d);

  edgeptsin[1]  = edgexyztemp;
  edgeptsin[0]  = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
  edgeptsin[2]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  Vx1yzm1pyr    =  E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxx1yzm1pyr, Vdyx1yzm1pyr, Vdzx1yzm1pyr);

  int totsurf2d = (demo.testcoord2dqqpts[0].size())*E->GetNcoeffs();
            
  Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);
  for(int k = 0; k < dimension-1; k++)
    {
      surfptsin[k] = Array<OneD, NekDouble>(demo.testcoord2dqqpts[k]);
      surfptsintemp[k] = Array<OneD, NekDouble>(demo.testcoord2dqqpts[k]);
    }
            
  surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.testcoord2dqqpts[0].size()); 
  surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.testcoord2dqqpts[0].size()); 
            
  //surface bot z = -1, (ABC)  
  Vxyzm1pyr = Array<OneD, NekDouble>(totsurf2d);
  surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[1] = surfptsin[1];
  Vxyzm1pyr = E->PhysEvaluateBasis(surfptsintemp, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  
  totsurf2d =  (demo.testcoord2dtqpts[0].size())*E->GetNcoeffs();
  //Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);

  for(int k = 0; k < dimension-1; k++)
    {
      surfptsin[k] = Array<OneD, NekDouble>(demo.testcoord2dtqpts[k]);
      surfptsintemp[k] = Array<OneD, NekDouble>(demo.testcoord2dtqpts[k]);
    }
  surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.testcoord2dtqpts[0].size());
  surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.testcoord2dtqpts[0].size());

  int totpt = surfptsin[0].size();
  
  //surface hypt (tri) x+z = 0 && y + z = 0   
  Vxmzypyr = Array<OneD, NekDouble>(totsurf2d);
  Vmath::Smul(totpt, -1.0, &surfptsintemp[0][0], 1, &surfptsintemp[2][0], 1);
  Vmath::Vcopy(totpt,  &surfptsintemp[0][0], 1,  &surfptsintemp[1][0], 1);
  Vxmzypyr = E->PhysEvaluateBasis(surfptsintemp, storage3d, NullNekDouble1DArray, NullNekDouble1DArray\
, NullNekDouble1DArray);

  //surface left (tri) x = -1   
  Vxm1yzpyr = Array<OneD, NekDouble>(totsurf2d);
  surfptsintemp[1] = surfptsin[1];
  surfptsintemp[2] = surfptsin[2];
  surfptsintemp[0] = Array<OneD, NekDouble>(totpt, -1.0);
  Vxm1yzpyr = E->PhysEvaluateBasis(surfptsintemp, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
  
  //surface front (tri) y = -1     
  Vxym1zpyr = Array<OneD, NekDouble>(totsurf2d);
  surfptsintemp[1] = Array<OneD, NekDouble>(totpt, -1.0);
  surfptsintemp[0] = surfptsin[0];
  Vxym1zpyr = E->PhysEvaluateBasis(surfptsintemp, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

 //surface back (tri) y = -z
  Vxymzpyr = Array<OneD, NekDouble>(totsurf2d);
  Vmath::Smul(totpt, -1.0, &surfptsintemp[2][0], 1, &surfptsintemp[1][0], 1);
  Vxymzpyr =  E->PhysEvaluateBasis(surfptsintemp, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

}



// FILE *csvfile(int flag)
// {
//     static char name[LOGNAME_SIZE];
//     char fullname[] = "Hex_";//[LOGNAME_SIZE+4];                                                                                                                                                           
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

//void WriteCSV(Array<OneD, NekDouble> uhats, int flag)
//{
  
// if flag == 0, before opt
// if flag == 1, after opt
// FILE *file = csvfile(flag);

// int uhatstot = uhats.size();
// int totpts = demo.testcoord3dlattice[0].size();
// Array<OneD, NekDouble> temp(uhatstot);
// for(int i = 0; i<totpts; i++)
// {
//     Vmath::Vmul(uhatstot, &demo.testcoord3dlattice[i], totpts, &uhats[0], 1, &temp[0], 1);
//     NekDouble v = Vmath::Vsum(uhatstot, temp, 1);
//     string s = std::to_string(demo.testcoord3dlattice[0][i])+","+std::to_string(demo.testcoord3dlattice[1][i])+","+std::to_string(demo.testcoord3dlattice[2][i])+","+std::to_string(v)+"\n";
//     const char* c = s.c_str();
//     fprintf (file, c); 
// }
//    fclose(file);
//}


void project_surfaces( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret)
{

  if(numsurfaces == 6) //hex
    {
      //bot surface  z = -1
      surfuhats(uhats, ret[0], Vxyzm1, Equad);

      //right surface x = 1
      surfuhats(uhats, ret[1], Vx1yz, Equad);

      //surface top z = 1
      surfuhats(uhats, ret[2], Vxyz1, Equad);
            
      //left surface x = -1
      surfuhats(uhats, ret[3], Vxm1yz, Equad);
       
      //surface front y = -1
      surfuhats(uhats, ret[4], Vxym1z, Equad);
       
      //surface back y = 1
      surfuhats(uhats, ret[5], Vxy1z, Equad);
    }
  else if(numsurfaces == 4) //tets
    {

      //surface bot z = -1 restriction on GD: x+y = 0 (ABC)
      surfuhats(uhats, ret[0],  Vxyzm1tet, Etri);

      //surface left x = -1 restriction on GD: y + z = 0 (DAC)
      surfuhats(uhats, ret[1], Vxm1ypz0tet, Etri);

      //surf front y = -1   restriction on GD: x + z = 0 (DAB)
      surfuhats(uhats, ret[2], Vxpz0ym1tet, Etri);

      //surf DCB restriction  on GD: (x + y + z = -1)
      surfuhats(uhats, ret[3], Vxpypzm1tet, Etri);

      
    }
  else if(numsurfaces == 5)
    {
      //surface bot z = -1   
      surfuhats(uhats, ret[0], Vxyzm1pyr, Equad);

      //surface hypt x+z = 0
      surfuhats(uhats, ret[1], Vxmzypyr, Etri);

	//surface left x = -1
      surfuhats(uhats, ret[2], Vxm1yzpyr, Etri);

	//surface front y = -1
      surfuhats(uhats, ret[3], Vxym1zpyr, Etri);
	
	//surface back y +z =  0
      surfuhats(uhats, ret[4], Vxymzpyr, Etri);
    }
}

void surfuhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz, StdExpansion *Etemp)
{
  int modes = uhats.size();

  Array<OneD, NekDouble> temp(uhats.size());  
  int totpts = Etemp->GetTotPoints();
  Array<OneD, NekDouble> vals(totpts);    
  //Vxyz*uhats -> project to -> E3tri or Equad
  for(int k = 0; k < totpts; k++)
    {
      
      Vmath::Vmul(uhats.size(), &Vxyz[k], totpts, &uhats[0], 1, &temp[0], 1);
      
      vals[k]  = Vmath::Vsum(modes, temp, 1);
      
    }
  //    vals = demo.blasmatvec(Vxyz, uhats, totpts, modes);
  
  Etemp->FwdTrans(vals, ret);
}



// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret)
{
  if(numedges == 12) // hex
    {

      // edge front left (x = -1) (y = -1)
      edgederpquhats(uhats, ret[0], Eorth->GetNcoeffs(), Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);                     
      //edge front right (x = 1) (y = -1)
      edgederpquhats(uhats, ret[1], Eorth->GetNcoeffs(), Vx1ym1z, Vdxx1ym1z, Vdyx1ym1z, Vdzx1ym1z);
            
      //edge front top (y = -1) (z = 1)
      edgederpquhats(uhats, ret[2], Eorth->GetNcoeffs(),  Vym1xz1, Vdxym1xz1, Vdyym1xz1, Vdzym1xz1);
    
      //edge front bot (y = -1) (z = -1)
      edgederpquhats(uhats, ret[3], Eorth->GetNcoeffs(),  Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);
     
      // edge back left (y = 1), (x = -1)
      edgederpquhats(uhats, ret[4], Eorth->GetNcoeffs(), Vxm1y1z, Vdxxm1y1z, Vdyxm1y1z, Vdzxm1y1z);
     
      //edge back right (x = 1), (y = 1))
      edgederpquhats(uhats, ret[5], Eorth->GetNcoeffs(),  Vx1y1z, Vdxx1y1z, Vdyx1y1z, Vdzx1y1z);
   
      //edge back top ( y = 1) (z = 1)
      edgederpquhats(uhats, ret[6], Eorth->GetNcoeffs(),  Vy1xz1, Vdxy1xz1, Vdyy1xz1, Vdzy1xz1);
   
      //edge back bot (y = 1) (z = -1))
      edgederpquhats(uhats, ret[7], Eorth->GetNcoeffs(),  Vy1xzm1, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);
            
      // edge left bot (z = -1), (x = -1)
      edgederpquhats(uhats, ret[8], Eorth->GetNcoeffs(), Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);
            
      //edge left top (x = -1), (z = 1))
      edgederpquhats(uhats, ret[9], Eorth->GetNcoeffs(),   Vxm1yz1,Vdxxm1yz1, Vdyxm1yz1, Vdzxm1yz1);
            
      //edge right bot ( z = -1) (x = 1)
      edgederpquhats(uhats, ret[10], Eorth->GetNcoeffs(),   Vx1yzm1, Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);
        
      
      //edge right top (z  1) (x  1))
      edgederpquhats(uhats, ret[11], Eorth->GetNcoeffs(),Vx1yz1, Vdxx1yz1, Vdyx1yz1, Vdzx1yz1); 

        
    }
  else if(numedges == 6) //tet
    {
      // edge front left (AD) (x = -1) (y = -1)
      edgederpquhats(uhats, ret[0], Eorth->GetNcoeffs(),  Vxm1ym1ztet, Vdxxm1ym1ztet, Vdyxm1ym1ztet, Vdzxm1ym1ztet);

      //edge front hypt (DB) (y = -1) (z = -x)
      edgederpquhats(uhats, ret[1], Eorth->GetNcoeffs(), Vym1xmztet, Vdxym1xmztet, Vdyym1xmztet, Vdzym1xmztet);

      //edge front bot (AB) (y = -1) (z = -1)
      edgederpquhats(uhats, ret[2], Eorth->GetNcoeffs(), Vym1xzm1tet, Vdxym1xzm1tet, Vdyym1xzm1tet, Vdzym1xzm1tet );

      //edge left hypt (DC) ( x = -1) (z = -y)
      edgederpquhats(uhats, ret[3], Eorth->GetNcoeffs(), Vxm1ymztet, Vdxxm1ymztet, Vdyxm1ymztet, Vdzxm1ymztet);

      // edge bot diag (BC) (z = -1) (y = -x)
      edgederpquhats(uhats, ret[4], Eorth->GetNcoeffs(),  Vxmyzm1tet, Vdxxmyzm1tet, Vdyxmyzm1tet, Vdzxmyzm1tet);

      //edge CA bot left (x = -1) (z = -1)
      edgederpquhats(uhats, ret[5], Eorth->GetNcoeffs(),Vxm1yzm1tet, Vdxxm1yzm1tet,  Vdyxm1yzm1tet,  Vdzxm1yzm1tet);

	
    }
  else if(numedges == 8) //pyr
    {
                                                                     
      // edge front left (x = -1) (y = -1) EA
      edgederpquhats(uhats, ret[0], Eorth->GetNcoeffs(), Vxm1ym1zpyr, Vdxxm1ym1zpyr, Vdyxm1ym1zpyr, Vdzxm1ym1zpyr);

      //edge front hypt (y = -1) (z = -x) EB
      edgederpquhats(uhats, ret[1], Eorth->GetNcoeffs(), Vym1xmzpyr, Vdxym1xmzpyr, Vdyym1xmzpyr, Vdzym1xmzpyr);

      //edge back hypt (z = -x  or x = -z) (x = y) EC
      edgederpquhats(uhats, ret[2], Eorth->GetNcoeffs(), Vxeyxmzpyr, Vdxxeyxmzpyr, Vdyxeyxmzpyr, Vdzxeyxmzpyr);

      //edge back left (y = -z) (x = -1) ED
      edgederpquhats(uhats, ret[3], Eorth->GetNcoeffs(), Vxeyxmzpyr, Vdxxeyxmzpyr, Vdyxeyxmzpyr, Vdzxeyxmzpyr); 

      //edge front bot (y = -1) (z = -1) AB
      edgederpquhats(uhats, ret[4], Eorth->GetNcoeffs(), Vym1xzm1pyr, Vdxym1xzm1pyr, Vdyym1xzm1pyr, Vdzym1xzm1pyr); 

      //edge back bot (y = 1) (z = -1)) DC
      edgederpquhats(uhats, ret[5], Eorth->GetNcoeffs(), Vy1xzm1pyr, Vdxy1xzm1pyr, Vdyy1xzm1pyr, Vdzy1xzm1pyr);

      // edge left bot (z = -1), (x = -1) AD
      edgederpquhats(uhats, ret[6], Eorth->GetNcoeffs(), Vxm1yzm1pyr, Vdxxm1yzm1pyr, Vdyxm1yzm1pyr, Vdzxm1yzm1pyr);

      //edge right bot ( z = -1) (x = 1) BC
      edgederpquhats(uhats, ret[7], Eorth->GetNcoeffs(), Vx1yzm1pyr, Vdxx1yzm1pyr, Vdyx1yzm1pyr, Vdzx1yzm1pyr);
    }
}
    
void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  int modes,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2)
{
  int uhatstot = uhats.size();
  int totpts = edgeptsin[0].size();

  Array<OneD, NekDouble> temp(totpts), temp2(modes), temp3(uhatstot);
  Array<OneD, NekDouble> pqeval(totpts);
  NekDouble v1, v2;
    
  for(int i = 0; i<totpts; i++)
    {
      Vmath::Vmul(modes, &Vxy[0]+i, totpts, &Vxyd0[0]+i, totpts, &temp2[0], 1);
      v1  = Vmath::Vsum(modes, temp2, 1); 
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
	  Vmath::Vmul(modes, &Vxy[0]+i, totpts, &Vxyd1[0]+i, totpts, &temp2[0], 1);
	  v1  = Vmath::Vsum(modes, temp2, 1); 
	  Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);

	  v2  = Vmath::Vsum(uhatstot, temp3, 1);  
	  v1 = v2*v1;

	  // At this point,                 
	  // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

	  // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &temp2[0], 1);

 
        
	  v2  = Vmath::Vsum(uhatstot, temp2, 1);  
	  Vmath::Vmul(uhats.size(), &Vxyd1[i], totpts, &uhats[0], 1, &temp2[0], 1);

	  v1= v2*Vmath::Vsum(uhats.size(), temp2, 1) - v1;
	  pqeval[i] += v1;
 
        }
      if(dimension == 3)
        {
	  Vmath::Vmul(modes, &Vxy[0]+i, totpts, &Vxyd2[0]+i, totpts, &temp2[0], 1);
	  v1  = Vmath::Vsum(modes, temp2, 1); 
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
    
  demo.E3seg->FwdTrans(pqeval, ret);
  
}

Array<OneD, NekDouble>  call_find_roots(Array<OneD,  NekDouble> &uhats , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats, NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime)
{
  int dimension = Eorth->GetShapeDimension(); 

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
  // EDGES
  if(numedges == 6) // tet
    {
      for(int i = 0; i < numedges; i++)
        {

	  t.Start();
	  rethold  = (demo.call_companion_rf(uhatsedges[i]));//, Eorth, NullNekDoubleArrayOfArray, dum,  0, 0, 0, 0));
          t.Stop();

          roots1dtime +=  t.TimePerTest(1);
	  for(int p = 0; p < rethold[0].size(); p++)
            {

	      switch(i)
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
  else if(numedges == 12) //hex
    {
      for(int i = 0; i < numedges; i++)
	{

	  t.Start();
	  rethold  = demo.call_companion_rf(uhatsedges[i]);//, Eorth, NullNekDoubleArrayOfArray, dum,  0, 0, 0, 0));
	  t.Stop();

	  roots1dtime +=  t.TimePerTest(1);
	  
	  for(int p = 0; p < rethold[0].size(); p++)
	    {
	      switch(i)
		{
		case 0: // edge front left (x = -1) (y = -1)
		  {
		    retall[0].push_back(   (-1.0));
		    retall[1].push_back( (-1.0));
		    retall[2].push_back( (rethold[0][p])); 

		  }
		  break;
		case 1: //edge front right (x = 1) (y = -1)
		  {
		    retall[0].push_back(  (1.0));
		    retall[1].push_back( (-1.0));
		    retall[2].push_back( (rethold[0][p]));

		  }
		  break;
		case 2: //edge front top (y = -1) (z = 1)
		  {
		    retall[0].push_back(  (rethold[0][p]));
		    retall[1].push_back( (-1.0));
		    retall[2].push_back( (1.0));
		  }
		  break;
		case 3: //edge front bot (y = -1) (z = -1)
		  {
		    retall[0].push_back( (rethold[0][p])); 
		    retall[1].push_back( (-1.0));
		    retall[2].push_back( (-1.0));
		  }
		  break;
		case 4: //edge back left (y = 1), (x = -1)
		  {
		    retall[0].push_back(  (-1.0));
		    retall[1].push_back( (1.0));
		    retall[2].push_back( (rethold[0][p]));
		  }
		  break;
		case 5: //edge back right (x = 1), (y = 1)
		  {
		    retall[0].push_back(  (1.0));
		    retall[1].push_back( (1.0));
		    retall[2].push_back( (rethold[0][p]));
		  }
		  break;
		case 6: //edge back top ( y = 1) (z = 1)
		  {
		    retall[0].push_back(  (rethold[0][p])); 
		    retall[1].push_back( (1.0));
		    retall[2].push_back( (1.0));

		  }
		  break;
		case 7: //edge back bot (y = 1) (z = -1)
		  {
		    retall[0].push_back( (rethold[0][p])); 
		    retall[1].push_back( (1.0));
		    retall[2].push_back( (-1.0));
		  }
		  break;
		case 8: //edge left bot (z = -1), (x = -1)
		  {
		    retall[0].push_back(  (-1.0));
		    retall[1].push_back( (rethold[0][p])); 
		    retall[2].push_back( (-1.0));
		  }
		  break;
		case 9: //edge left top (x = -1), (z = 1)
		  {
		    retall[0].push_back(  (-1.0));
		    retall[1].push_back( (rethold[0][p])); 
		    retall[2].push_back( (1.0));

		  }
		  break;
		case 10: //edge right bot (z = -1) (x = 1)
		  {
		    retall[0].push_back(  (1.0));
		    retall[1].push_back( (rethold[0][p])); 
		    retall[2].push_back( (-1.0));
		  }
		  break;
		case 11: //edge right top (z  1) (x  1)
		  {
		    retall[0].push_back(  (1.0));
		    retall[1].push_back( (rethold[0][p])); 
		    retall[2].push_back( (1.0));;
		  }
		  break;
		default:
		  {
		    cout<<"\n edge number not valid!\n\n";
		    exit(0);
		  }
		}
	    }
    	    
	}
    }
  else if(numedges == 8) //pyr
    {
       for(int i = 0; i < numedges; i++)
        {
	  t.Start();
          rethold  = demo.call_companion_rf(uhatsedges[i]);//, Eorth, NullNekDoubleArrayOfArray, dum,  0, 0, 0, 0));
          t.Stop();
	  
          roots1dtime +=  t.TimePerTest(1);
	  if(rethold != NullNekDoubleArrayOfArray )
	    {
	      for(int p = 0; p < rethold[0].size(); p++)
		{
		  if(i == 0)
		    { 
		      retall[0].push_back(-1);
		      retall[1].push_back(-1);
		      retall[2].push_back(rethold[0][p]);
		    }
		  else if(i == 1) //edge front hypt (y = -1) (z = -x)                                    
		    {
		      retall[0].push_back(rethold[0][p]);
		      retall[2].push_back(-rethold[0][p]);
		      retall[1].push_back(-1);

		      retall[0].push_back(-rethold[0][p]);
		      retall[2].push_back(rethold[0][p]);
		      retall[1].push_back(-1);
		    }
		  else if(i == 2) //edge back hypt ( y = -z) (z +x = 0)                                  
		    {
		      retall[0].push_back(rethold[0][p]);
		      retall[1].push_back(rethold[0][p]);
		      retall[2].push_back(-rethold[0][p]);

		      retall[0].push_back(-rethold[0][p]);
		      retall[1].push_back(-rethold[0][p]);
		      retall[2].push_back(rethold[0][p]);

		    }
		  else if(i == 3)  //edge back left (y = -z) (x = -1) ED                                 
		    {
		      retall[0].push_back(-1);
		      retall[1].push_back(-rethold[0][p]);
		      retall[2].push_back(rethold[0][p]);

		      retall[0].push_back(-1);
		      retall[1].push_back(rethold[0][p]);
		      retall[2].push_back(-rethold[0][p]);
		    }
		  else if(i == 4)  //edge front bot (y = -1) (z = -1)                                    
		    {
		      retall[1].push_back(-1);
		      retall[0].push_back(rethold[0][p]);
		      retall[2].push_back(-1);
		    }
		  else if(i == 5) //edge back bot (y = 1) (z = -1)                                       
		    {
		      retall[1].push_back(1);
		      retall[0].push_back(rethold[0][p]);
		      retall[2].push_back(-1);
		    }
		  else if(i == 6)//edge left bot (z = -1), (x = -1)                                      
		    {
		      retall[0].push_back(-1);
		      retall[1].push_back(rethold[0][p]);
		      retall[2].push_back(-1);
		    }
		  else if(i == 7) //edge right bot ( z = -1) (x = 1)                                     
		    {
		      retall[2].push_back(-1);
		      retall[1].push_back(rethold[0][p]);
		      retall[0].push_back(1);
		    }
		}
	    }
	}
    }

  //SURFACES:
  if(numsurfaces == 4) //tet
    {

      NekDouble avgiterGDhold;

      // call 2D rootfinder on each surface:
      for(int i = 0; i < numsurfaces; i++)
	{
	  t.Start();

	  demo.steepestgradient_descent2Dtri(surfaceuhats[i], Etri, storage2dt, rethold,  avgiterGDhold);
	  if(rethold != NullNekDoubleArrayOfArray)
	    {
	      //find_roots(surfaceuhats[i]);//, Etri, storage2dt,  avgiterGDhold, 0, 1, 0, 1);
	      // last arg = 1 coz all 2d ele are triangles

	      t.Stop();
	      roots2dtime +=  t.TimePerTest(1);
	      avgiterGD += avgiterGDhold;

	      switch(i)
		{
		case 0:
		  //surface bot z = -1, x + y = 0 (ABC)
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
		case 2:
		  //surf front y = -1, x + z = 0 (DAB)
		  retall[2].push_back(rethold[1][0]);
		  retall[1].push_back(-1);
		  retall[0].push_back(rethold[0][0]);
		  break;
		case 3:
		  //surf DCB (x + y + z = -1)
		  retall[2].push_back(-1.0 - rethold[1][0] - rethold[0][0]);
		  retall[1].push_back(rethold[1][0]);
		  retall[0].push_back(rethold[0][0]);
		  break;
		default: cout<<"\n invalid surface id\n";
		  exit(0);
		}
	    }
	}
	  
    }
  else if(numsurfaces == 6) //hex
    {
      NekDouble avgiterGDhold;
      for(int i = 0; i < numsurfaces; i++)
	{
	  
	  t.Start();   
	  demo.steepestgradient_descent2Dquad(surfaceuhats[i], Equad, storage2dq, rethold,  avgiterGDhold);//find_roots(surfaceuhats[i], Equad, storage2dq,  avgiterGDhold, 0, 1, 0 , 0);
	  t.Stop();

	  roots2dtime +=  t.TimePerTest(1);
	  avgiterGD += avgiterGDhold;
	  if(rethold == NullNekDoubleArrayOfArray)
	    {
	      continue;
	      
	    }
	  if(i == 0 ) // surf bot (z = -1)
	    {

	      retall[0].push_back( (rethold[0][0]));
	      retall[1].push_back( (rethold[1][0]));
	      retall[2].push_back( (-1.0));

	    }
	  else if(i == 1) //surf right (x = 1)
	    {
	      retall[0].push_back  (1.0);
	      retall[1].push_back (rethold[0][0]);
	      retall[2].push_back (rethold[1][0]);
	    
	    }
	  else if(i == 2) //surf top (z = 1)
	    {
	      retall[0].push_back  (rethold[0][0]);
	      retall[1].push_back (rethold[1][0]);
	      retall[2].push_back (1.0);
	    }
	  else if(i == 3) //surf left (x = -1)
	    {
	      retall[0].push_back  (-1.0);
	      retall[1].push_back (rethold[0][0]);
	      retall[2].push_back (rethold[1][0]);
	    }
	  else if(i == 4) //surf front (y = -1)           
	    {
	      retall[0].push_back  (rethold[0][0]);
	      retall[1].push_back (-1.0);
	      retall[2].push_back (rethold[1][0]);
	    }
	  else if(i == 5) //surf back (y = 1)
	    {
	      retall[0].push_back  (rethold[0][0]);
	      retall[1].push_back (1.0);
	      retall[2].push_back (rethold[1][0]);
	    }
	    
	}
      
    }
  else if(numsurfaces == 5) //pyr
    {
      NekDouble avgiterGDhold;
      for(int i = 0; i < numsurfaces; i++)
	{
	  if( i == 0) // bot surface: quad, z = -1
	    {
	      t.Start();   
	      demo.steepestgradient_descent2Dquad(surfaceuhats[i], Equad, storage2dq, rethold,  avgiterGDhold);
	      t.Stop();
	      
	      roots2dtime +=  t.TimePerTest(1);
	      avgiterGD += avgiterGDhold;
	      if(rethold != NullNekDoubleArrayOfArray)
		{
		  retall[2].push_back(-1);
		  retall[1].push_back(rethold[1][0]);
		  retall[0].push_back(rethold[0][0]);
		}
	    }
	  else  // rest surfaces: tri
	    {
	      t.Start();
	      demo.steepestgradient_descent2Dtri(surfaceuhats[i], Etri, storage2dt, rethold,  avgiterGDhold);
	      //rethold = demo.find_roots(surfaceuhats[i], Etri, storage2dt,  avgiterGDhold, 0, 1, 0 , 0);
              t.Stop();
	      roots2dtime +=  t.TimePerTest(1);
              avgiterGD += avgiterGDhold;
	      if(rethold != NullNekDoubleArrayOfArray)
		{
		  switch(i)
		    {
		    case 1: // x+z = 0
		      retall[0].push_back(rethold[0][0]);
		      retall[1].push_back(rethold[1][0]);
		      retall[2].push_back(-rethold[0][0]);
		      break;
		    case 2: // x = -1 
		      retall[2].push_back(rethold[1][0]);
		      retall[1].push_back(rethold[0][0]);
		      retall[0].push_back(-1.0);
		      break;
		    case 3: //y = -1
		      retall[2].push_back(rethold[1][0]);
		      retall[1].push_back(-1.0);
		      retall[0].push_back(rethold[0][0]);
		      break;
		    case 4: // y+z = 0
		      retall[2].push_back(-rethold[1][0]);
		      retall[1].push_back(rethold[1][0]);
		      retall[0].push_back(rethold[0][0]);
		  
		    }
		}
	    }
	}
    }
    
  //find interior roots:
  t.Start();
  demo.steepestgradientdescent3D(uhats, Eorth, storage3d, rethold, avgiterGDhold);
  //(demo.find_roots(uhats, E,storage3d,  avgiterGDhold, 0, 0, 1, 0)) ;
  t.Stop();
  roots3dtime +=  t.TimePerTest(1);  

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
  Array<OneD, NekDouble> evalroots = Eorth->PhysEvaluateBasis(tmpcoord, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
     
  Array<OneD, NekDouble> minvandx(dimension+1);
  Array<OneD, NekDouble> tmparr(retall[0].size());


  demo.pq(uhats, tmpcoord, evalroots, minvandx, tmparr);


  tempmin = minvandx[0];
  ret[0] = minvandx[1];
  ret[1] = minvandx[2];
  ret[2] = minvandx[3];
	
  //}
  minv = tempmin;
  return ret;
    
}    

int Opt_needed(Array<OneD, NekDouble> uhats, int flag)
{
  int nq = E->GetTotPoints();
  Array<OneD, NekDouble> temp1(nq);
  Array<OneD, NekDouble> temp2(uhats.size());
  Array<OneD, Array<OneD, NekDouble> > coordquad = demo.GetCoords(E);
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
  
  int ct = nq;
  nq = demo.testcoord3dhqmidpts[0].size();
  temp1=Array<OneD, NekDouble>(nq);
  cout<<"\n all opt and vals:\n";
  for (int i = 0; i<nq; i++)
    {
      Vmath::Vmul(uhats.size(), &demo.interioreval3hqmidpts[i], nq, &uhats[0], 1, &temp2[0], 1);
      temp1[i] = Vmath::Vsum(uhats.size(), temp2, 1);
      if(temp1[i] < minv)
	{
	  minv = temp1[i];
	  idxmin = i+ct;
	}
    }
  if(minv < 0 && abs(minv) > 1e-10 && idxmin != -1)
    { 
      if(flag == 0)
	{
	  startval = minv;
	  startcoordx = demo.testcoord3dhlattice[0][idxmin];
	  startcoordy = demo.testcoord3dhlattice[1][idxmin];
	  startcoordz = demo.testcoord3dhlattice[2][idxmin];
	}
      cout<<"\n minv = "<<minv<<" at "<< demo.testcoord3dhlattice[0][idxmin]<<  " "<<demo.testcoord3dhlattice[1][idxmin]<<" "<< demo.testcoord3dhlattice[2][idxmin];;;
	  
      return 1;   

    }
  cout<<"\n minv = "<<minv<<" at "<< demo.testcoord3dhlattice[0][idxmin]<<  " "<<demo.testcoord3dhlattice[1][idxmin]<<" "<< demo.testcoord3dhlattice[2][idxmin];;;

  return 0;
  
}

void Do_optimize(Array<OneD, NekDouble> &uhats)
{

  int dim = Eorth->GetShapeDimension();
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

  //    int NC = 1; //number of constraints, only positivity for now
  tols.push_back(1e-12);

  Array<OneD, NekDouble> optima(dim), pqvalcoords(dim+1), xastarr(dim), utemp(N1), wsp1(1); 
  Array<OneD, Array<OneD, NekDouble> > surfaceuhats  (numsurfaces), tmpcoord(dim), Pf(numedges);
  for(int k= 0; k < numedges; k++)
    {
      Pf[k] = Array<OneD, NekDouble>(demo.E3seg->GetNcoeffs());
    }

	  
  for(int k = 0; k < numsurfaces; k++)
    {
      if(Eorth->DetShapeType() == LibUtilities::eHexahedron)
	{
	  surfaceuhats[k] = Array<OneD, NekDouble>(Equad->GetNcoeffs());
	}
      if( Eorth->DetShapeType() == LibUtilities::eTetrahedron)
	{
	  surfaceuhats[k] = Array<OneD, NekDouble>(Etri->GetNcoeffs());
	}
      if(Eorth->DetShapeType() == LibUtilities::ePyramid)
	{
	  if(k == 0 )
	    {
	       surfaceuhats[k] = Array<OneD, NekDouble>(Equad->GetNcoeffs());
	    }
	  else
	    {
	      surfaceuhats[k] = Array<OneD, NekDouble>(Etri->GetNcoeffs());
	    }
	}
	  
    }

  NekDouble pqval, timeprojectsurf = 0.0, timeprojectedges = 0.0;
  NekDouble avgiterGD = 0.0, avgiterhold = 0.0, roots1dtimehold = 0.0, roots2dtimehold = 0.0, roots3dtimehold = 0.0;
  Array<OneD, NekDouble> Vtmp(  Eorth->GetNcoeffs());

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

	  project_edges(utemp, Pf);

	  t.Stop();
	  timeprojectedges += t.TimePerTest(1);
	  t.Start();
	  project_surfaces(utemp, surfaceuhats);
	  t.Stop();
	  timeprojectsurf += t.TimePerTest(1);
	  optima = call_find_roots(utemp, avgiterhold, Pf, surfaceuhats,  minv, roots1dtime, roots2dtime, roots3dtime);
	  roots1dtimehold += roots1dtime;
	  roots2dtimehold += roots2dtime;
	  roots3dtimehold += roots3dtime;
	  avgiterGD += avgiterhold;    
    
	  //	  cout<<"\n opt = ";
	  // Vtmp is evaluation of basis fn at optima
	  for(int k = 0; k < dim; k++)
	    {
	      //  cout<<" "<<optima[k];
	      tmpcoord[k] = Array<OneD, NekDouble>(1, optima[k]);
	    }
	
	  Vtmp = Eorth->PhysEvaluateBasis(tmpcoord, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

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

	  Vtmp = Eorth->PhysEvaluateBasis(tmpcoord, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

	}
      demo.pq(utemp, tmpcoord, Vtmp, pqvalcoords, wsp1);
      //      cout<<"\n pqval = "<<pqval<<" " <<pqvalcoords[0];
      if (pqvalcoords[0] < pqval)
        {
	  xastarr = optima;
	  pqval = pqvalcoords[0];
        }
        
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
	  tmp = Eorth->PhysEvaluateBasis(xastarrofarr, storage3d, i);
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

      counter = counter + 1;
    
    
    }
  roots1dtimehold = roots1dtimehold;///(counter-1);
  roots2dtimehold = roots2dtimehold;///(counter-1);
  roots3dtimehold = roots3dtimehold;///(counter-1);
  timeprojectedges = timeprojectedges;//s/(counter);
  timeprojectsurf = timeprojectsurf;///(counter);    
  itersGD2 = (avgiterGD);///(counter);
  iterstaken = counter;
  cout<<"\nsphere_rotation took "<<counter<<"iterations";//

  //cout<<"\n  startv = "<<startval<<" at startcoords("<<startcoordx<<","<<startcoordy<<" ," <<startcoordz<<")    GD iter stepszie ="<<demo.GetTol()<<" maxiters = "<<demo.GetMaxIter()<<" o = "<<Eorth->GetBasis(0)->GetNumModes()<<" p = "<<Eorth->GetBasis(0)->GetNumPoints()<<" eps = "<<demo.GetEps()<<"  avgiterGD = "<<itersGD2;//<<" total time taken post rootfinding per iter = "<< evalrootshold/(counter)<<" hex \n";;

  cout<<"\nAvg times per iter (1drootfinder,2drootfinder,3drootfinder) = : "  << roots1dtimehold<<", "<<roots2dtimehold<<", "<<roots3dtimehold<<" itergd2 = "<<itersGD2; ;//<< timeprojectedges<<", "<<timeprojectsurf<<",

  uhats = d.back();
}


 
      // // calculating reference time here:
      //      NekDouble elapref = 0;

      // Array<OneD, NekDouble> physder(phys.size()), physder0(phys.size()), uhatsnew(coeffs.size(), 1.0);

      // // demo.GetAvgNum() refers to the user defined parameter to get timing values for demo.GetAvgNum() number of runs 
      // for(int k = 0; k < demo.GetAvgNum(); k++)
      // 	{
      
      // 	  timer.Start();
      // 	  Eorth->BwdTrans(uhatsnew, physder);
      // 	  Eorth->PhysDeriv(physder, physder0);
      // 	  Eorth->FwdTrans(physder0, uhatsnew); 
      // 	  timer.Stop();
      // 	  elapref += timer.TimePerTest(1);
      // 	}

      // elapref = elapref/demo.GetAvgNum();
      // cout<<"\n Ref time (first derivative uhats from f(x) uhats took) "<<(elapref)<<" s\n\n";

      //      int dimension = Eorth->GetShapeDimension(); 
        
