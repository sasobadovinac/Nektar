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
void Do_optimize(StdExpansion *exp, Array<OneD, NekDouble> &uhats);

//surface roots
void surfuhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz, StdExpansion* Etemp );

void call_setup_quad(StdExpansion *E);
void call_setup_hex(StdExpansion *E);
void call_setup_tet(StdExpansion *E);
void call_setup_pyr(StdExpansion *E);
void call_setup_tri(StdExpansion *E);  

// declare caller routine to find_roots
Array<OneD, NekDouble> call_find_roots(Array<OneD, NekDouble> &uhatsall, NekDouble &avgiterGD, Array< OneD, Array<OneD, NekDouble> >&uhatsedges, Array< OneD, Array<OneD, NekDouble> >&surfaceuhats,  NekDouble &minv, NekDouble &roots1dtime,  NekDouble &roots2dtime , NekDouble &roots3dtime );

// for 2D elements to get uhats at edges
void project_edges( Array<OneD, NekDouble> uhats,Array<OneD, Array<OneD, NekDouble> >&ret, StdExpansion *Eseghold);

int Opt_needed(Array<OneD, NekDouble> uhats, Array<OneD, NekDouble> &latticeeval, int flag = 0);
  
// derivative of scaled function stuff
void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret, int modes,   Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1 = NullNekDouble1DArray , Array<OneD, NekDouble> Vxyd2 = NullNekDouble1DArray);

void project_surfaces( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret);

Array<OneD, NekDouble> FindIds(Array<OneD, NekDouble> &uhatslocal, StdExpansion *Eorth, NekDouble tol);

Array<OneD, NekDouble> eval_sol(Array<OneD, NekDouble> x, Array<OneD, NekDouble> y, Array<OneD, NekDouble> z);

string funcdef;
int iterstaken, verboseflag;
Array<OneD, Array<OneD, NekDouble> > storage3dhold;
Array<OneD, Array<OneD, NekDouble> > storage2dt;
Array<OneD, Array<OneD, NekDouble> > storage2dq;

// Confederate matrix
int numedges, numsurfaces;
StdExpansion *E, *Eorth;
StdExpansion *Equad;
StdExpansion *Etri;

Array<OneD, Array<OneD, NekDouble> > edgeptsin;
// for edges root finding quad:
// left (x = -1)
Array<OneD, NekDouble> Vxm1q;
Array<OneD, NekDouble> Vdyxm1q;
Array<OneD, NekDouble> Vdxxm1q;

// bot (y = -1)
Array<OneD, NekDouble> Vym1q;
Array<OneD, NekDouble> Vdxym1q;
Array<OneD, NekDouble> Vdyym1q;
      
// right quad (x = 1)
Array<OneD, NekDouble> Vx1q;
Array<OneD, NekDouble> Vdxx1q;
Array<OneD, NekDouble> Vdyx1q;

// top quad (y = 1)
Array<OneD, NekDouble> Vdxy1q;
Array<OneD, NekDouble> Vdyy1q;
Array<OneD, NekDouble> Vy1q;

// for edges root finding tri:
// left (x = -1)
Array<OneD, NekDouble> Vxm1t;
Array<OneD, NekDouble> Vdyxm1t;
Array<OneD, NekDouble> Vdxxm1t;

// bot (y = -1)
Array<OneD, NekDouble> Vym1t;
Array<OneD, NekDouble> Vdxym1t;
Array<OneD, NekDouble> Vdyym1t;

// hypt tri (y = -x)
Array<OneD, NekDouble> Vxyhypt;
Array<OneD, NekDouble> Vdxxyhypt;
Array<OneD, NekDouble> Vdyxyhypt;
      
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

//surface hypt x+z = 0, y+z = 0
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
Array<OneD, Array<OneD, NekDouble> > allptslattice, coordsE, midcoordsE;
int main(int argc, char *argv[])
{
  demo.GetOptions().add_options()("optm,z", "positivity preserving optimizer");
  demo.ParseArguments(argc, argv);
  po::variables_map vm = demo.GetVariableMap();
   
  verboseflag = demo.GetVerbose();
  E = demo.CreateStdExpansion();

  coordsE = demo.GetCoords(E);
  midcoordsE = demo.GetQuadratureMidCoords(coordsE);
  allptslattice = demo.GetLatticeCoords(coordsE, midcoordsE);

  LibUtilities::PointsKey pkeycheb( 3*(E->GetBasis(0)->GetNumPoints()),LibUtilities::eGaussLobattoChebyshev);
  LibUtilities::BasisKey bkeyorth(LibUtilities::eOrtho_A, 3*(E->GetBasis(0)->GetNumModes()),  pkeycheb);
  demo.E3seg = new StdSegExp(bkeyorth);
  if (demo.E3seg == nullptr)
    {
      cout<<"\n Could not initialize segment object! \n";
      exit(0);
    }
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

  for (int i = 0; i < dimension; ++i)
    {
      btype[i] = E->GetBasisType(i);
      order.push_back(E->GetBasisNumModes(i));
    }

  switch(E->DetShapeType())
    {
    case LibUtilities::eQuadrilateral:
      call_setup_quad(E);

      break;
    case LibUtilities::eTriangle:
      call_setup_tri(E);
      break;
    
    case LibUtilities::eTetrahedron:
      {
	call_setup_tri(E); 
	call_setup_tet(E);
	numedges = 6;
	numsurfaces = 4;

      }
      break;

    case LibUtilities::eHexahedron:
      call_setup_quad(E);
      call_setup_hex(E);
        
      numedges = 12;
      numsurfaces = 6;

      break;

    case LibUtilities::ePyramid:
      call_setup_tri(E);
      call_setup_quad(E);
      call_setup_pyr(E);
      numedges = 8;
      numsurfaces = 5;
      
      break;

    default: cout<<"\n unknown shape typen\n";exit(0);
    
    }
  storage3dhold = Eorth->GetPhysEvaluateStorage();

      
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
    
  sol = eval_sol(x, y, z);
  
  Array<OneD, NekDouble> phys(totPoints), sollattice;
  Array<OneD, NekDouble> coeffs((unsigned) E->GetNcoeffs());
  //NekDouble L2errlattice = 0, Linferrlattice = 0;

  sollattice = eval_sol(allptslattice[0], allptslattice[1], allptslattice[2] );

  
  //Project onto expansion
  E->FwdTrans(sol, coeffs);

  //Backward transform solution to get projected values
  E->BwdTrans(coeffs, phys);
  
  LibUtilities::Timer     timer;

  Array<OneD, NekDouble> retvalslattice, solhold(sollattice.size());
  
  
  
  // check for -ve values and apply opt if necessary
  if (vm.count("optm"))
    {

      int opt_need = 0;
      if(1)
	{
	  E->BwdTrans(coeffs, phys);
	  Eorth->FwdTrans(phys, coeffs);

	}
      cout<<"\n1\n\n";

      timer.Start();
      Array<OneD, NekDouble> chk = FindIds(coeffs, Eorth, 1e-7);
      if(chk.size() != 0)
	{
	  opt_need = 1;
	}

      //int opt_need = Opt_needed(coeffs, retvalslattice, 0);            
      timer.Stop();
      NekDouble  timeneg = timer.TimePerTest(1); 
      if(verboseflag)
	{
	  cout<<"\nchecking for -vity took "<<timeneg<<"s\n";
	}
       
      if(opt_need)
	{
	  NekDouble  timeverify;
	  Array<OneD, BasisType> btorth(3);

	  NekDouble elap = 0;
	  timer.Start();
	  Do_optimize(Eorth, coeffs);
	  timer.Stop();
	  elap += timer.TimePerTest(1); 
	  if(verboseflag)
	    {
	      cout<<"\noptimizertook "<<elap<<"s";
	    }
	  Array<OneD, NekDouble > tmp = FindIds(coeffs, Eorth, 1e-7);
	  
	  if(1)//orthoflag)
	    {
	      // Eorth->BwdTrans(coeffs, phys);
	      // E->FwdTrans(phys, coeffs);
	      // demo.OrthoNormalize(E, coeffs, btorth, 1); 
	      // storage3dhold = E->GetPhysEvaluateStorage();
	      Eorth->BwdTrans(coeffs, phys);
	      E->FwdTrans(phys, coeffs);

	    }
	  timer.Start();
	  if(tmp.size() == 0)//valnow < 0 && abs(valnow)>1e-7)
	    {
	      opt_need = 0;
	    }
	  timer.Stop();
	  timeverify = timer.TimePerTest(1);
	  cout<<"\nVerification took "<<timeverify<<"s";
	  //L2 error between retvalslattice and eval_sol on lattice:
	  // Vmath::Vcopy(sollattice.size(), sollattice, 1, solhold, 1);
	  // Vmath::Vsub(latticesz, retvalslattice, 1, sollattice, 1,
	  // 	      sollattice, 1);
	  // Vmath::Vabs(latticesz, sollattice, 1, sollattice, 1);
	  // Linferrlattice = Vmath::Vmax(latticesz, sollattice, 1);
	  // Vmath::Vmul(latticesz, sollattice, 1, sollattice, 1, sollattice, 1);
	  // L2errlattice = sqrt(Vmath::Vsum(latticesz, sollattice, 1) );
	  
	  if(opt_need)
	    {
	      //cout<<"func="<<funcdef;
	  
	      cout<<"\n\n fail\n";
	      //sabotage result of test here:
	      cout <<"L infinity error: \t" << 1 << endl;
	      cout << "L2 error \t: \t" << 1 << endl;
	      return 0;
	    }
	  else
	    {
	      cout<<"\n\n pass\n\n";

	      if(verboseflag)
	       	{
		  // cout << "constrained L infinity error: \t" << Linferrlattice << endl;
		  // cout << "constrained L2 error \t: \t" << L2errlattice << endl;
      		}  
	      //		    cout<<"\n "<<" "<<startval<<", ("<<startcoordx<<" "<<startcoordy<<" "<<startcoordz<<"), "<<setup<<", pass, "<<" "<<itersGD1<<", "<<itersGD2<<", "<<itersGD3<<", "<<iterstaken<<", "<<timeneg<<", "<<elap<<", "<<timeverify;//<<", "<<elapref<<"\n";;
		  
	      //		}
	      //Calculate L_inf & L_2 error
	      cout << "L infinity error: \t" << 0.0 << endl;
	      if (stype != ePoint)
		{
		  cout << "L 2 error: \t \t \t" << 0.0 << endl;
		}
	      //WriteCSV(coeffs, 1);
	      return 0;
	  
	    }
        }
      else
	{
	  
	  cout<<"\\n Filter found no negative vals!\n\n";
	  //Calculate L_inf & L_2 error
	  cout << "L infinity error: \t" << 0.0 << endl;
	  if (stype != ePoint)
	    {
	      cout << "L 2 error: \t \t \t" << 0.0 << endl;
	    }
	}
      
    }
  else
    {
      cout<<"\n optimizer flag not set. Add -z option\n";
    }

  //Calculate L_inf & L_2 error
  cout << "L infinity error: \t" << E->Linf(phys, sol) << endl;
  if (stype != ePoint)
    {
      cout << "L 2 error: \t \t \t" << E->L2(phys, sol) << endl;
    }
  
  return 0;
}

void call_setup_tri(StdExpansion * exp)
{
  int nmodes0 = exp->GetBasis(0)->GetNumModes();
      
  int nmodes1 = exp->GetBasis(1)->GetNumModes();
  int npts0 = exp->GetBasis(0)->GetNumPoints();
  int npts1 = exp->GetBasis(1)->GetNumPoints();
  PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
  PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
  BasisKey  b0(LibUtilities::eOrtho_A,  nmodes0, p0);
  BasisKey  b1(LibUtilities::eOrtho_B,  nmodes1, p1);
  Etri = new StdTriExp(b0, b1);
      
  demo.storage2dt =  Etri->GetPhysEvaluateStorage(); 
  demo.coordtri = demo.GetCoords(Etri);
  demo.coordmidtri = demo.GetQuadratureMidCoords(demo.coordtri);
  demo.coordlatticetri = demo.GetLatticeCoords(demo.coordtri, demo.coordmidtri);

  demo.midptevaltri = Etri->PhysEvaluateBasis(demo.coordmidtri, demo.storage2dt, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  Array<OneD, Array<OneD, NekDouble> > edgeptsin(2), edgexy(2);
  for(int k = 0; k <2; k++)
    {
      edgexy[k] = Array<OneD, NekDouble>(demo.E3seg->GetBasis(0)->GetNumPoints());
      edgeptsin[k] =   Array<OneD, NekDouble>(edgexy[k]);
    }
  Array<OneD, NekDouble> edgexytemp = demo.E3seg->GetBasis(0)->GetZ();
  int totszedges = edgexytemp.size()*(Etri->GetNcoeffs());
      
  // left (x = -1)
  Vxm1t = Array<OneD, NekDouble>(totszedges);
  Vdyxm1t = Array<OneD, NekDouble>(totszedges);
  Vdxxm1t = Array<OneD, NekDouble>(totszedges);
      
  // left x = -1
  edgexy[1] = edgexytemp;
  edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
  Vxm1t = Etri->PhysEvaluateBasis(edgexy, demo.storage2dt, Vdxxm1t, Vdyxm1t, NullNekDouble1DArray);
      
  // bot (y = -1)
  Vym1t  = Array<OneD, NekDouble>(totszedges);
  Vdxym1t  = Array<OneD, NekDouble>(totszedges);
  Vdyym1t  = Array<OneD, NekDouble>(totszedges);
      
  // bot y = -1
  edgexy[0] = edgexytemp;             
  edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
      
  Vym1t = Etri->PhysEvaluateBasis(edgexy,  demo.storage2dt, Vdxym1t, Vdyym1t, NullNekDouble1DArray);
      
  // hypt tri (y = -x)
  Vxyhypt = Array<OneD, NekDouble>(totszedges);
  Vdxxyhypt = Array<OneD, NekDouble>(totszedges);
  Vdyxyhypt = Array<OneD, NekDouble>(totszedges);
      
  edgexy[0] = edgexytemp;
  Vmath::Smul(edgexy[0].size(), -1.0, edgexy[0], 1, edgexy[1] , 1);
  Vxyhypt = Etri->PhysEvaluateBasis(edgexy, demo.storage2dt,  Vdxxyhypt, Vdyxyhypt, NullNekDouble1DArray);
      
}

    
void call_setup_quad(StdExpansion* exp)
{
  // get # of quad pts in curr exp
  // get order of current exp
  // get ptypes and btypes
  int nmodes0 = exp->GetBasis(0)->GetNumModes();
  int nmodes1 = exp->GetBasis(1)->GetNumModes();
  int npts0 = exp->GetBasis(0)->GetNumPoints();
  int npts1 = exp->GetBasis(1)->GetNumPoints();
      
  PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
  PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
  BasisKey  b0(LibUtilities::eOrtho_A,  nmodes0, p0);
  BasisKey  b1(LibUtilities::eOrtho_A,  nmodes1, p1);
  Equad = new StdQuadExp(b0, b1);
  demo.storage2dq =  Equad->GetPhysEvaluateStorage();
  demo.coordquad = demo.GetCoords(Equad);
  demo.coordmidquad = demo.GetQuadratureMidCoords(demo.coordquad);
  demo.coordlatticequad = demo.GetLatticeCoords(demo.coordquad, demo.coordmidquad);

  demo.midptevalquad = Equad->PhysEvaluateBasis(demo.coordmidquad, demo.storage2dq, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
      

  Array<OneD, Array<OneD, NekDouble> > edgeptsin(2), edgexy(2);
  for(int k = 0; k <2; k++)
    {
      edgexy[k] = Array<OneD, NekDouble>(demo.E3seg->GetBasis(0)->GetNumPoints());
      edgeptsin[k] =   Array<OneD, NekDouble>(edgexy[k]);
    }
  Array<OneD, NekDouble> edgexytemp =   demo.E3seg->GetBasis(0)->GetZ();
  int totszedges = edgexytemp.size()*(Equad->GetNcoeffs());
      
  // left (x = -1)
  Vxm1q = Array<OneD, NekDouble>(totszedges);
  Vdyxm1q = Array<OneD, NekDouble>(totszedges);
  Vdxxm1q = Array<OneD, NekDouble>(totszedges);
      
  // left x = -1
  edgexy[1] = edgexytemp;
  edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
  Vxm1q = Equad->PhysEvaluateBasis(edgexy, demo.storage2dq, Vdxxm1q, Vdyxm1q, NullNekDouble1DArray);
      
  // bot (y = -1)
  Vym1q  = Array<OneD, NekDouble>(totszedges);
  Vdxym1q  = Array<OneD, NekDouble>(totszedges);
  Vdyym1q  = Array<OneD, NekDouble>(totszedges);
      
  // bot y = -1
  edgexy[0] = edgexytemp;             
  edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
      
  Vym1q = Equad->PhysEvaluateBasis(edgexy,  demo.storage2dq, Vdxym1q, Vdyym1q, NullNekDouble1DArray);

      	  
  // right quad (x = 1)
  Vx1q = Array<OneD, NekDouble>(totszedges);
  Vdxx1q = Array<OneD, NekDouble>(totszedges);
  Vdyx1q = Array<OneD, NekDouble>(totszedges);
      
  // top quad (y = 1)
  Vdxy1q = Array<OneD, NekDouble>(totszedges);
  Vdyy1q = Array<OneD, NekDouble>(totszedges);
  Vy1q = Array<OneD, NekDouble>(totszedges);
      
  // right x = 1
  edgexy[1] = edgexytemp; 
  edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), 1.0); 
  Vx1q = Equad->PhysEvaluateBasis(edgexy, demo.storage2dq, Vdxx1q, Vdyx1q, NullNekDouble1DArray);
      
  //top y = 1
  edgexy[0] = edgexytemp; 
  edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), 1.0); 
  Vy1q = Equad->PhysEvaluateBasis(edgexy, demo.storage2dq, Vdxy1q, Vdyy1q, NullNekDouble1DArray);
      
}

void call_setup_tet(StdExpansion *exp)
{
  int orthoflag = (exp->GetBasisType(0) != LibUtilities::eOrtho_A);
  if(!orthoflag)
    {
      int nmodes0 = exp->GetBasis(0)->GetNumModes();
      int nmodes1 = exp->GetBasis(1)->GetNumModes();
      int nmodes2 = exp->GetBasis(2)->GetNumModes();
      int npts0 = exp->GetBasis(0)->GetNumPoints();	
      int npts1 = exp->GetBasis(1)->GetNumPoints();
      int npts2 = exp->GetBasis(2)->GetNumPoints();
      
      PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
      PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
      PointsKey p2(npts2, exp->GetBasis(2)->GetPointsType());
      BasisKey  b0(LibUtilities::eOrtho_A,  nmodes0, p0);
      BasisKey  b1(LibUtilities::eOrtho_B,  nmodes1, p1);
      BasisKey  b2(LibUtilities::eOrtho_C,  nmodes2, p2);
      
      Eorth = new StdTetExp(b0, b1, b2);
    }
  else
    {
      Eorth = exp;
    }
  
  demo.storage3dtet = Eorth->GetPhysEvaluateStorage();
  demo.coordtet = demo.GetCoords(Eorth);
  demo.coordmidtet = demo.GetQuadratureMidCoords(demo.coordtet);
  demo.coordlatticetri = demo.GetLatticeCoords(demo.coordtet, demo.coordmidtet);
  
  demo.midptevaltet = Eorth->PhysEvaluateBasis(demo.coordmidtet, demo.storage3dtet, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  Array<OneD, NekDouble> edgexyz = (demo.E3seg->GetBasis(0)->GetZ());
  int totszedges1d =   edgexyz.size()*(Eorth->GetNcoeffs());
  
  edgeptsin = Array<OneD, Array<OneD, NekDouble> >(dimension);
  for(int p = 0; p < dimension; p++)
    {
      edgeptsin[p] = Array<OneD, NekDouble>(edgexyz.size(), &edgexyz[0]);
    } 

  // edge front left (AD) (x = -1) (y = -1)                                                
  Vxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

  Vxm1ym1ztet = Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dtet, Vdxxm1ym1ztet, Vdyxm1ym1ztet, Vdzxm1ym1ztet);

  //edge front hypt (DB) (y = -1) (z = -x)                                                 
  Vym1xmztet = Array<OneD, NekDouble>(totszedges1d);
  Vdxym1xmztet = Array<OneD, NekDouble>(totszedges1d);
  Vdyym1xmztet = Array<OneD, NekDouble>(totszedges1d);
  Vdzym1xmztet = Array<OneD, NekDouble>(totszedges1d);
  
  edgeptsin[0] =   Array<OneD, NekDouble>(edgexyz);   
  Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[1][0], 1, &edgeptsin[2][0], 1);

  Vym1xmztet = Eorth->PhysEvaluateBasis(edgeptsin,demo.storage3dtet, Vdxym1xmztet, Vdyym1xmztet, Vdzym1xmztet);

  //edge front bot (AB) (y = -1) (z = -1)                                                
  Vym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
  Vdxym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
  Vdyym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
  Vdzym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[0] =   Array<OneD, NekDouble>(edgexyz);   

  Vym1xzm1tet = Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dtet, Vdxym1xzm1tet, Vdyym1xzm1tet, Vdzym1xzm1tet);

  //edge left hypt (DC) ( x = -1) (z = -y)                                               
  Vxm1ymztet = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1ymztet= Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1ymztet = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1ymztet= Array<OneD, NekDouble>(totszedges1d);

  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[1] =  Array<OneD, NekDouble>(edgexyz);
  Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[1][0], 1, &edgeptsin[2][0], 1);
  
  Vxm1ymztet = Eorth->PhysEvaluateBasis(edgeptsin,demo.storage3dtet, Vdxxm1ymztet, Vdyxm1ymztet, Vdzxm1ymztet);

  // edge bot diag (BC) (z = -1) (y = -x)                                                
  Vxmyzm1tet = Array<OneD, NekDouble>(totszedges1d);
  Vdxxmyzm1tet = Array<OneD, NekDouble>(totszedges1d)   ;
  Vdyxmyzm1tet = Array<OneD, NekDouble>(totszedges1d)    ;
  Vdzxmyzm1tet = Array<OneD, NekDouble>(totszedges1d)    ;

  edgeptsin[0] =   Array<OneD, NekDouble>(edgexyz);  
  edgeptsin[2] =   Array<OneD, NekDouble>(edgexyz);  
  Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[0][0], 1, &edgeptsin[1][0], 1);
  Vxmyzm1tet = Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dtet, Vdxxmyzm1tet, Vdyxmyzm1tet, Vdzxmyzm1tet);
  
  //edge CA bot left (x = -1) (z = -1)                                                   
  Vxm1yzm1tet   = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1yzm1tet = Array<OneD, NekDouble>(totszedges1d)   ;
  Vdyxm1yzm1tet = Array<OneD, NekDouble>(totszedges1d)    ;
  Vdzxm1yzm1tet = Array<OneD, NekDouble>(totszedges1d)    ;
  edgeptsin[1] =   Array<OneD, NekDouble>(edgexyz);  ;
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[2].size(), -1.0);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  
  Vxm1yzm1tet = Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dtet, Vdxxm1yzm1tet, Vdyxm1yzm1tet, Vdzxm1yzm1tet);
  int totsurf2d =  (demo.coordtri[0].size())*Eorth->GetNcoeffs();
  Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);
  
  for(int k = 0; k < dimension-1; k++)
    {
      surfptsin[k] = Array<OneD, NekDouble>(demo.coordtri[k]);
      surfptsintemp[k] = Array<OneD, NekDouble>(demo.coordtri[k]);
    }
  surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.coordtri[0].size());
  surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.coordtri[0].size());

  int totpt = surfptsin[0].size();

  //surface bot z = -1, (ABC) 
  Vxyzm1tet = Array<OneD, NekDouble>(totsurf2d);
  surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[1] = surfptsin[1];
  Vxyzm1tet = Eorth->PhysEvaluateBasis(surfptsintemp, demo.storage3dtet, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  //surface left x = -1  (DAC)  
  Vxm1ypz0tet =  Array<OneD, NekDouble>(totsurf2d);
  surfptsintemp[0] = Array<OneD, NekDouble> (totpt, -1.0);
  surfptsintemp[1] = surfptsin[0];
  surfptsintemp[2] = surfptsin[1];
  
  Vxm1ypz0tet = Eorth->PhysEvaluateBasis(surfptsintemp, demo.storage3dtet, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  //surf front y = -1,  (DAB) 
  Vxpz0ym1tet =  Array<OneD, NekDouble>(totsurf2d);
  surfptsintemp[1] = Array<OneD, NekDouble> (totpt, -1.0);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[2] = surfptsin[1];
  Vxpz0ym1tet = Eorth->PhysEvaluateBasis(surfptsintemp, demo.storage3dtet, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  //surf DCB (x + y + z = -1), 
  Vxpypzm1tet =  Array<OneD, NekDouble>(totsurf2d);
  surfptsintemp[1] = surfptsin[1];;

  Vmath::Vadd(totpt, &surfptsin[0][0], 1, &surfptsin[1][0], 1, &surfptsintemp[2][0], 1);
  Vmath::Smul(totpt, -1.0, &surfptsintemp[2][0], 1, &surfptsintemp[2][0], 1);
  Vmath::Sadd(totpt, -1.0, &surfptsintemp[2][0], 1, &surfptsintemp[2][0], 1);
  Vxpypzm1tet = Eorth->PhysEvaluateBasis(surfptsintemp,demo.storage3dtet , NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
  
}

void call_setup_hex(StdExpansion *exp)
{
  
  // get # of quad pts in curr exp
  // get order of current exp
  // get ptypes and btypes
  int orthoflag = (exp->GetBasisType(0) != LibUtilities::eOrtho_A);
  if(!orthoflag)
    {

      int nmodes0 = exp->GetBasis(0)->GetNumModes();
      int nmodes1 = exp->GetBasis(1)->GetNumModes();
      int nmodes2 = exp->GetBasis(2)->GetNumModes();
      
      int npts0 = exp->GetBasis(0)->GetNumPoints();
      int npts1 = exp->GetBasis(1)->GetNumPoints();
      int npts2 = exp->GetBasis(2)->GetNumPoints();
      
      PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
      PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
      PointsKey p2(npts2, exp->GetBasis(2)->GetPointsType());
      BasisKey  b0(LibUtilities::eOrtho_A,  nmodes0, p0);
      BasisKey  b1(LibUtilities::eOrtho_A,  nmodes1, p1);
      BasisKey  b2(LibUtilities::eOrtho_A,  nmodes2, p2);
      Eorth = new StdHexExp(b0, b1, b2);
    }
  else
    {
      Eorth = exp;
    }
  demo.storage3dhex =  Eorth->GetPhysEvaluateStorage();
  
  demo.coordhex = demo.GetCoords(Eorth);
  demo.coordmidhex = demo.GetQuadratureMidCoords(demo.coordhex);
  demo.coordlatticehex = demo.GetLatticeCoords(demo.coordhex, demo.coordmidhex);
  
  demo.midptevalhex = Eorth->PhysEvaluateBasis(demo.coordmidhex, demo.storage3dhex, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  Array<OneD, NekDouble> edgexyz = demo.E3seg->GetBasis(0)->GetZ();

  int totszedges1d = edgexyz.size()*(Eorth->GetNcoeffs());
  edgeptsin = Array<OneD, Array<OneD, NekDouble> >(dimension);    
  for(int p = 0; p < dimension; p++)
    {
      edgeptsin[p] = Array<OneD, NekDouble>(edgexyz.size(), &edgexyz[0]);
    } 

  // edge front left (x = -1) (y = -1)
  Vxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

            
  Vxm1ym1z = Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);  
                
  //edge front right (x = 1) (y = -1)
  Vx1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1ym1z = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);       
  Vx1ym1z = Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxx1ym1z, Vdyx1ym1z, Vdzx1ym1z);  


  //edge front top (y = -1) (z = 1)
  Vym1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxym1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyym1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzym1xz1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = edgexyz; 
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
  Vym1xz1 = Eorth->PhysEvaluateBasis(edgeptsin,demo.storage3dhex, Vdxym1xz1, Vdyym1xz1, Vdzym1xz1);  

                
  //edge front bot (y = -1) (z = -1)
  Vym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  Vym1xzm1 =             Eorth->PhysEvaluateBasis( edgeptsin, demo.storage3dhex, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);  
	    
	    
  // edge back left (y = 1), (x = -1)
  Vxm1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1y1z = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);      
  edgeptsin[2] =    edgexyz;


  Vxm1y1z = Eorth->PhysEvaluateBasis(edgeptsin,demo.storage3dhex, Vdxxm1y1z, Vdyxm1y1z, Vdzxm1y1z);  

  //edge back right (x = 1), (y = 1))
  Vx1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1y1z = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);                
         
  Vx1y1z  = Eorth->PhysEvaluateBasis(edgeptsin,demo.storage3dhex, Vdxx1y1z, Vdyx1y1z, Vdzx1y1z);  

                
  //edge back top ( y = 1) (z = 1)
  Vy1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxy1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyy1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzy1xz1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = edgexyz;
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
  Vy1xz1 =             Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxy1xz1, Vdyy1xz1, Vdzy1xz1 );  

      
  //edge back bot (y = 1) (z = -1))
  Vy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
            
  Vy1xzm1 = Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);  

  // edge left bot (z = -1), (x = -1)
  Vxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[1] = edgexyz;
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

  Vxm1yzm1 = Eorth->PhysEvaluateBasis(edgeptsin,demo.storage3dhex, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);  
                
  //edge left top (x = -1), (z = 1))
  Vxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
                
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);           
  Vxm1yz1 = Eorth->PhysEvaluateBasis(edgeptsin,demo.storage3dhex, Vdxxm1yz1, Vdyxm1yz1, Vdzxm1yz1);  

  //edge right bot ( z = -1) (x = 1)
  Vx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1yzm1 = Array<OneD, NekDouble>(totszedges1d);

  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);           
  Vx1yzm1 = Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);   

  //edge right top (z  1) (x  1))
  Vx1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1yz1 = Array<OneD, NekDouble>(totszedges1d);
                
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);

  Vx1yz1 = Eorth->PhysEvaluateBasis(edgeptsin,demo.storage3dhex, Vdxx1yz1, Vdyx1yz1, Vdzx1yz1);   

              
  int totszsurf2d = (demo.coordquad[0].size())*Eorth->GetNcoeffs();
            
  Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);
  for(int k = 0; k < dimension-1; k++)
    {
      surfptsin[k] = Array<OneD, NekDouble>(demo.coordquad[k]);
      surfptsintemp[k] = Array<OneD, NekDouble>(demo.coordquad[k]);
    }
            
  surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.coordquad[0].size()); 
  surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.coordquad[0].size()); 
            
  //surface bot z = -1
  Vxyzm1 =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[1] = surfptsin[1];
  Vxyzm1 =             Eorth->PhysEvaluateBasis(surfptsintemp,demo.storage3dhex, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

  //surface right x = 1
  Vx1yz =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[0] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
  surfptsintemp[1] =  surfptsin[0];
  surfptsintemp[2] =  surfptsin[1];
  Vx1yz = Eorth->PhysEvaluateBasis(surfptsintemp,demo.storage3dhex,NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   


  //surface top z = 1
  Vxyz1 =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[1] = surfptsin[1];
  surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
  Vxyz1 = Eorth->PhysEvaluateBasis(surfptsintemp, demo.storage3dhex, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

  //surface left x = -1
  Vxm1yz =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[0] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
  surfptsintemp[1] = surfptsin[0];
  surfptsintemp[2] = surfptsin[1];
  Vxm1yz = Eorth->PhysEvaluateBasis(surfptsintemp,demo.storage3dhex, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

  //surface front y = -1
  Vxym1z =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[1] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[2] = surfptsin[1];
  Vxym1z = Eorth->PhysEvaluateBasis(surfptsintemp,demo.storage3dhex, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   
      
  //surface back y = 1
  Vxy1z =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[1] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[2] = surfptsin[1];
  Vxy1z = Eorth->PhysEvaluateBasis(surfptsintemp,demo.storage3dhex, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
}

void call_setup_pyr(StdExpansion *exp)
{
  int orthoflag = (exp->GetBasisType(0) != LibUtilities::eOrtho_A);
  if(!orthoflag)
    {
      int nmodes0 = exp->GetBasis(0)->GetNumModes();
      int nmodes1 = exp->GetBasis(1)->GetNumModes();
      int nmodes2 = exp->GetBasis(2)->GetNumModes();
      int npts0 = exp->GetBasis(0)->GetNumPoints();
      int npts1 = exp->GetBasis(1)->GetNumPoints();
      int npts2 = exp->GetBasis(2)->GetNumPoints();
      PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
      PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
      PointsKey p2(npts2, exp->GetBasis(2)->GetPointsType());
      BasisKey  b0(LibUtilities::eOrtho_A,  nmodes0, p0);
      BasisKey  b1(LibUtilities::eOrtho_A,  nmodes1, p1);
      BasisKey  b2(LibUtilities::eOrthoPyr_C,  nmodes2, p2);
      Eorth = new StdPyrExp(b0, b1, b2);
    }
  else
    {
      Eorth = E;
    }
  demo.storage3dpyr = Eorth->GetPhysEvaluateStorage();
  demo.coordpyr = demo.GetCoords(Eorth);
  demo.coordmidpyr = demo.GetQuadratureMidCoords(demo.coordpyr);
  demo.coordlatticepyr = demo.GetLatticeCoords(demo.coordpyr, demo.coordmidpyr);
  
  demo.midptevalpyr = Eorth->PhysEvaluateBasis(demo.coordmidpyr, demo.storage3dpyr, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
  
  Array<OneD, NekDouble> edgexyz = (demo.E3seg->GetBasis(0)->GetZ());
  int totszedges1d =   edgexyz.size()*(Eorth->GetNcoeffs());
  
  edgeptsin = Array<OneD, Array<OneD, NekDouble> >(dimension);    
  for(int p = 0; p < dimension; p++)
    {
      edgeptsin[p] = Array<OneD, NekDouble>(edgexyz.size(), &edgexyz[0]);
    } 
  
  // edge front left EA (x = -1) (y = -1)                 
  Vxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);

  edgeptsin[0]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[1]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  Vxm1ym1zpyr   = Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxxm1ym1zpyr, Vdyxm1ym1zpyr, Vdzxm1ym1zpyr);

  //edge front hypt EB (y = -1) (z + x = 0)
  Vym1xmzpyr    = Array<OneD, NekDouble>(totszedges1d);
  Vdxym1xmzpyr  = Array<OneD, NekDouble>(totszedges1d);
  Vdyym1xmzpyr  = Array<OneD, NekDouble>(totszedges1d);
  Vdzym1xmzpyr  = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] =  Array<OneD, NekDouble>(edgexyz);;
  Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[0][0], 1, &edgeptsin[2][0], 1);

  Vym1xmzpyr =  Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxym1xmzpyr, Vdyym1xmzpyr, Vdzym1xmzpyr);

  //edge back hypt EC: x+z = 0  and x = y // (z = -x  or x = -z) (x = y) EC
  Vxeyxmzpyr    = Array<OneD, NekDouble>(totszedges1d);
  Vdxxeyxmzpyr  = Array<OneD, NekDouble>(totszedges1d);
  Vdyxeyxmzpyr  = Array<OneD, NekDouble>(totszedges1d);
  Vdzxeyxmzpyr  = Array<OneD, NekDouble>(totszedges1d);

  edgeptsin[1]  =  Array<OneD, NekDouble>(edgexyz);
  Vxeyxmzpyr    =  Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxxeyxmzpyr, Vdyxeyxmzpyr, Vdzxeyxmzpyr);
  //edge back left (y = -z) (x = -1) ED
  Vx1ymzpyr     = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1ymzpyr   = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1ymzpyr   = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1ymzpyr   = Array<OneD, NekDouble>(totszedges1d); ;
  edgeptsin[0]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[1]  = Array<OneD, NekDouble>(edgexyz);
  Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[1][0], 1, &edgeptsin[2][0], 1);
  Vx1ymzpyr     = Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxx1ymzpyr, Vdyx1ymzpyr, Vdzx1ymzpyr);
  
  //edge front bot (y = -1) (z = -1) AB
  Vym1xzm1pyr   = Array<OneD, NekDouble>(totszedges1d);
  Vdxym1xzm1pyr = Array<OneD, NekDouble>(totszedges1d)   ;
  Vdyym1xzm1pyr = Array<OneD, NekDouble>(totszedges1d)    ;
  Vdzym1xzm1pyr = Array<OneD, NekDouble>(totszedges1d)    ;

  edgeptsin[1]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[2]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[0]  = Array<OneD, NekDouble>(edgexyz);
  Vym1xzm1pyr   = Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxym1xzm1pyr, Vdyym1xzm1pyr,Vdzym1xzm1pyr);

  
  //edge back bot (y = 1) (z = -1)) DC
  Vy1xzm1pyr    = Array<OneD, NekDouble>(totszedges1d);
  Vdyy1xzm1pyr  = Array<OneD, NekDouble>(totszedges1d)   ;
  Vdxy1xzm1pyr  = Array<OneD, NekDouble>(totszedges1d);
  Vdzy1xzm1pyr  = Array<OneD, NekDouble>  (totszedges1d) ;

  edgeptsin[1]  = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
  Vy1xzm1pyr   = Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxy1xzm1pyr, Vdyy1xzm1pyr, Vdzy1xzm1pyr);

  
  // edge left bot (z = -1), (x = -1) AD
  Vxm1yzm1pyr   = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1yzm1pyr = Array<OneD, NekDouble>(totszedges1d)   ;
  Vdxxm1yzm1pyr = Array<OneD, NekDouble>(totszedges1d)   ;
  Vdzxm1yzm1pyr = Array<OneD, NekDouble>(totszedges1d)   ;

  edgeptsin[1]  =  Array<OneD, NekDouble>(edgexyz);
  edgeptsin[0]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  Vxm1yzm1pyr   = Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxxm1yzm1pyr, Vdyxm1yzm1pyr, Vdzxm1yzm1pyr);

  
  //edge right bot ( z = -1) (x = 1) BC
  Vx1yzm1pyr    = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1yzm1pyr  = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1yzm1pyr  = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1yzm1pyr  = Array<OneD, NekDouble>(totszedges1d);

  edgeptsin[0]  = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
  Vx1yzm1pyr    = Eorth->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxx1yzm1pyr, Vdyx1yzm1pyr, Vdzx1yzm1pyr);

  int totsurf2d = (demo.coordquad[0].size())*Eorth->GetNcoeffs();
            
  Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);
  for(int k = 0; k < dimension-1; k++)
    {
      surfptsin[k] = Array<OneD, NekDouble>(demo.coordquad[k]);
      surfptsintemp[k] = Array<OneD, NekDouble>(demo.coordquad[k]);
    }
            
  surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.coordquad[0].size()); 
  surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.coordquad[0].size()); 
            
  //surface bot z = -1, (ABC)  
  Vxyzm1pyr = Array<OneD, NekDouble>(totsurf2d);
  surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[1] = surfptsin[1];
  Vxyzm1pyr = Eorth->PhysEvaluateBasis(surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  
  totsurf2d =  (demo.coordtri[0].size())*Eorth->GetNcoeffs();
  //Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);

  for(int k = 0; k < dimension-1; k++)
    {
      surfptsin[k] = Array<OneD, NekDouble>(demo.coordtri[k]);
      surfptsintemp[k] = Array<OneD, NekDouble>(demo.coordtri[k]);
    }
  surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.coordtri[0].size());
  surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.coordtri[0].size());

  int totpt = surfptsin[0].size();
  
  //surface hypt (tri)   z = -x   //old:x+z = 0 && y + z = 0   
  Vxmzypyr = Array<OneD, NekDouble>(totsurf2d);
  //leave x and y alone: they only need to be x+y <= 0 inside GD and z = -1 here
  Vmath::Smul(totpt, -1.0, &surfptsintemp[0][0], 1, &surfptsintemp[2][0], 1);
  //  surfptsintemp[1] = surfptsin[0];
  //Vmath::Vcopy(totpt,  &surfptsintemp[0][0], 1,  &surfptsintemp[1][0], 1);
  Vxmzypyr = Eorth->PhysEvaluateBasis(surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  //surface left (tri) x = -1   
  Vxm1yzpyr = Array<OneD, NekDouble>(totsurf2d);
  surfptsintemp[1] = surfptsin[1];
  surfptsintemp[2] = surfptsin[2];
  surfptsintemp[0] = Array<OneD, NekDouble>(totpt, -1.0);
  Vxm1yzpyr = Eorth->PhysEvaluateBasis(surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
  
  //surface front (tri) y = -1     
  Vxym1zpyr = Array<OneD, NekDouble>(totsurf2d);
  surfptsintemp[1] = Array<OneD, NekDouble>(totpt, -1.0);
  surfptsintemp[0] = surfptsin[0];
  Vxym1zpyr = Eorth->PhysEvaluateBasis(surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  //surface back (tri) y = -z
  Vxymzpyr = Array<OneD, NekDouble>(totsurf2d);
  Vmath::Smul(totpt, -1.0, &surfptsintemp[2][0], 1, &surfptsintemp[1][0], 1);
  Vxymzpyr =  Eorth->PhysEvaluateBasis(surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

}

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
  else if(numsurfaces == 5)// pyr
    {
      //surface bot z = -1   
      surfuhats(uhats, ret[0], Vxyzm1pyr, Equad);

      //surface hypt x+z = 0 and y +z = 0 -> x-y = 0
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



void project_edges( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret, StdExpansion *E)
{
  if(E->DetShapeType() == LibUtilities::eTriangle) // tri
    {
	  
      // bot edge  
      edgederpquhats(uhats, ret[0], E->GetNcoeffs(),  Vym1t, Vdxym1t, Vdyym1t,NullNekDouble1DArray);
      // left edge
      edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vxm1t, Vdxxm1t, Vdyxm1t, NullNekDouble1DArray);
	  
      // hypto
      edgederpquhats(uhats, ret[2], E->GetNcoeffs(),  Vxyhypt, Vdxxyhypt, Vdyxyhypt, NullNekDouble1DArray);
	        
      return;
    }
      
  if(E->DetShapeType() == LibUtilities::eHexahedron) // hex
    {

      // edge front left (x = -1) (y = -1)                                                         
      edgederpquhats(uhats, ret[0], E->GetNcoeffs(), Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);
      //edge front right (x = 1) (y = -1)                                                          
      edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vx1ym1z, Vdxx1ym1z, Vdyx1ym1z, Vdzx1ym1z);

      //edge front top (y = -1) (z = 1)                                                            
      edgederpquhats(uhats, ret[2], E->GetNcoeffs(),  Vym1xz1, Vdxym1xz1, Vdyym1xz1, Vdzym1xz1);

      //edge front bot (y = -1) (z = -1)                                                           
      edgederpquhats(uhats, ret[3], E->GetNcoeffs(), Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);

      // edge back left (y = 1), (x = -1)                                                          
      edgederpquhats(uhats, ret[4], E->GetNcoeffs(), Vxm1y1z, Vdxxm1y1z, Vdyxm1y1z, Vdzxm1y1z);

      //edge back right (x = 1), (y = 1))                                                          
      edgederpquhats(uhats, ret[5], E->GetNcoeffs(),  Vx1y1z, Vdxx1y1z, Vdyx1y1z, Vdzx1y1z);

      //edge back top ( y = 1) (z = 1)                                                             
      edgederpquhats(uhats, ret[6], E->GetNcoeffs(),  Vy1xz1, Vdxy1xz1, Vdyy1xz1, Vdzy1xz1);

      //edge back bot (y = 1) (z = -1))                                                            
      edgederpquhats(uhats, ret[7], E->GetNcoeffs(),  Vy1xzm1, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);

      // edge left bot (z = -1), (x = -1)                                                          
      edgederpquhats(uhats, ret[8], E->GetNcoeffs(), Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);
  
      //edge left top (x = -1), (z = 1))                                                           
      edgederpquhats(uhats, ret[9], E->GetNcoeffs(),   Vxm1yz1,Vdxxm1yz1, Vdyxm1yz1, Vdzxm1yz1);
  
      //edge right bot ( z = -1) (x = 1)                                                           
      edgederpquhats(uhats, ret[10], E->GetNcoeffs(),   Vx1yzm1, Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);


      //edge right top (z  1) (x  1))                                                              
      edgederpquhats(uhats, ret[11], E->GetNcoeffs(),Vx1yz1, Vdxx1yz1, Vdyx1yz1, Vdzx1yz1);

  
    }

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
  Array<OneD, Array<OneD, NekDouble > >  rethold(dimension), storage;
  vector<vector<NekDouble> > retall(dimension);
  Array<OneD, NekDouble> ret(dimension);
  Timer t1;
  // EDGES
  if(numedges == 4)
    {
      dimension = 2;
      retall = vector<vector<NekDouble> >(dimension); 
      int numedges = 4;
      //	  NekDouble dum;
      for(int i = 0; i < numedges; i++)
	{
	  storage = demo.storage2dq;
	  t1.Start();
	  rethold  = demo.call_companion_rf(uhatsedges[i]);
	  t1.Stop();
	  roots1dtime +=  t1.TimePerTest(1);
	  if(rethold != NullNekDoubleArrayOfArray)
	    {
	      for(int p = 0; p < rethold[0].size(); p++)
		{

		  switch(i)
		    {
		    case 0:
		      //if(i == 0) // edge bot (y = -1)                                               
		      retall[1].push_back(   (-1.0));
		      retall[0].push_back( (rethold[0][p]));
		      break;
		    case 1:
		      //if(i == 1) // edge right (x = 1)                                              
		      retall[0].push_back(   (1.0));
		      retall[1].push_back( (rethold[0][p]));
		      break;
		    case 2:
		      //if(i == 2) // edge top (y = 1)                                                
		      retall[1].push_back(   (1.0));
		      retall[0].push_back( (rethold[0][p]));
		      break;
		    case 3:
		      //if(i == 3) // edge left (x = -1)                                              
		      retall[0].push_back(   (-1.0));
		      retall[1].push_back( (rethold[0][p]));
		      break;
		    default:
		      break;
		    }
		}
	    }
	}
      // add 4 corners                                                                             
      retall[0].push_back(-1);
      retall[1].push_back(-1);
	  
      retall[0].push_back(1);
      retall[1].push_back(1);
	  
      retall[0].push_back(-1);
      retall[1].push_back(1);
	  
      retall[0].push_back(1);
      retall[1].push_back(-1);

      rethold = NullNekDoubleArrayOfArray;

      t1.Start();
      demo.steepestgradient_descent2D(uhats, Equad, rethold, avgiterGDhold, 1e-7); 
      t1.Stop();
      roots2dtime +=  t1.TimePerTest(1);
      avgiterGD += avgiterGDhold;
      //	  cout<<"\n GD: quad returns:\n";
      if(rethold != NullNekDoubleArrayOfArray)
	{
	  retall[0].push_back(rethold[0][0]);
	  retall[1].push_back(rethold[1][0]);
	  //  cout<<rethold[0][0]<<" ,"<<rethold[1][0];
	}
		  
    }
  if(numedges == 6) // tet
    {
      for(int i = 0; i < numedges; i++)
        {
	  t.Start();
	  rethold  = (demo.call_companion_rf(uhatsedges[i]));
          t.Stop();

          roots1dtime +=  t.TimePerTest(1);
	  if(rethold != NullNekDoubleArrayOfArray && rethold[0] != NullNekDouble1DArray)
	    {
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
    }
  else if(numedges == 12) //hex
    {
      for(int i = 0; i < numedges; i++)
	{

	  t.Start();
	  rethold  = demo.call_companion_rf(uhatsedges[i]);
	  t.Stop();

	  roots1dtime +=  t.TimePerTest(1);
	  if(rethold != NullNekDoubleArrayOfArray)
	    {
	      
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
    }
  else if(numedges == 8) //pyr
    {
      for(int i = 0; i < numedges; i++)
        {
	  t.Start();
          rethold  = demo.call_companion_rf(uhatsedges[i]);
          t.Stop();
	  
          roots1dtime +=  t.TimePerTest(1);
	  if(rethold != NullNekDoubleArrayOfArray )
	    {
	      for(int p = 0; p < rethold[0].size(); p++)
		{
		  if(i == 0)  // edge front left EA (x = -1) (y = -1) 
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
		  else if(i == 2) //  edge back hypt ( y = -z) (z +x = 0)                                  
		    {
		      retall[0].push_back(rethold[0][p]);
		      retall[1].push_back(rethold[0][p]);
		      retall[2].push_back(-rethold[0][p]);

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
	  demo.steepestgradient_descent2D(surfaceuhats[i], Etri, rethold, avgiterGDhold, 1e-7); 

	  t.Stop();
	  roots2dtime +=  t.TimePerTest(1);
	  avgiterGD += avgiterGDhold;

	  
	  if(rethold != NullNekDoubleArrayOfArray)
	    {

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

      //find interior roots:
      t.Start();
      demo.steepestgradientdescent3D(uhats, Eorth, storage3dhold, demo.coordtet, demo.coordmidtet, demo.midptevaltet, rethold, avgiterGDhold);
      t.Stop();
      roots3dtime +=  t.TimePerTest(1);  
      
      avgiterGD += avgiterGDhold;
    }
  else if(numsurfaces == 6) //hex
    {
      NekDouble avgiterGDhold;
      for(int i = 0; i < numsurfaces; i++)
	{
	  
	  t.Start();   
	  demo.steepestgradient_descent2D(surfaceuhats[i], Equad, rethold,  avgiterGDhold, 1e-7);
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
      
      //find interior roots:
      t.Start();
      demo.steepestgradientdescent3D(uhats, Eorth, storage3dhold, demo.coordhex, demo.coordmidhex, demo.midptevalhex, rethold, avgiterGDhold);
      t.Stop();
      roots3dtime +=  t.TimePerTest(1);  
      
      avgiterGD += avgiterGDhold;
    }
  else if(numsurfaces == 5) //pyr
    {
      NekDouble avgiterGDhold;
      for(int i = 0; i < numsurfaces; i++)
	{
	  if( i == 0) // bot surface: quad, z = -1
	    {
	      t.Start();   
	      demo.steepestgradient_descent2D(surfaceuhats[i], Equad, rethold,  avgiterGDhold, 1e-7);
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
	      demo.steepestgradient_descent2D(surfaceuhats[i], Etri, rethold,  avgiterGDhold, 1e-7);
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
      //find interior roots:
      t.Start();
      demo.steepestgradientdescent3D(uhats, Eorth, storage3dhold, demo.coordpyr, demo.coordmidpyr, demo.midptevalpyr, rethold, avgiterGDhold);
      t.Stop();
      roots3dtime +=  t.TimePerTest(1);  
      
      avgiterGD += avgiterGDhold;
    }
    

  Array<OneD, Array<OneD, NekDouble> > tmpcoord(dimension);
  if(rethold!=NullNekDoubleArrayOfArray)
    {
      if(rethold[0][0] < inf && rethold[1][0] < inf && rethold[2][0] < inf)
	{
	  retall[0].push_back(rethold[0][0]);
	  retall[1].push_back(rethold[1][0]);
	  retall[2].push_back(rethold[2][0]);

	}
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
  Array<OneD, NekDouble> evalroots = Eorth->PhysEvaluateBasis(tmpcoord, storage3dhold, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
     
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

void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  int modes,  Array <OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2)
{
  int uhatstot = uhats.size();
  int totpts =  demo.E3seg->GetBasis(0)->GetZ().size();

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

Array<OneD, NekDouble> FindIds(Array<OneD, NekDouble> &uhatslocal, StdExpansion *exp, NekDouble tol)
{
  boost::ignore_unused(tol);
  Array<OneD, NekDouble> minvalandpthold;
  // // check on mid-quadrature pts of staggered lattice

  switch(exp->DetShapeType())
    {
    case LibUtilities::eQuadrilateral:
      minvalandpthold = demo.FindLatticeEval(uhatslocal, demo.storage2dq, demo.coordquad, demo.midptevalquad, demo.coordmidquad);
      break;
    case LibUtilities::eTriangle:
      minvalandpthold = demo.FindLatticeEval(uhatslocal, demo.storage2dt, demo.coordtri, demo.midptevaltri, demo.coordmidtri);
      break;
    case LibUtilities::eHexahedron:

      minvalandpthold = demo.FindLatticeEval(uhatslocal, demo.storage3dhex, demo.coordhex, demo.midptevalhex, demo.coordmidhex);
      break;
    case LibUtilities::eTetrahedron:
      minvalandpthold = demo.FindLatticeEval(uhatslocal, demo.storage3dtet, demo.coordtet, demo.midptevaltet, demo.coordmidtet);
      break;
    case LibUtilities::ePyramid:
      minvalandpthold = demo.FindLatticeEval(uhatslocal, demo.storage3dpyr, demo.coordpyr, demo.midptevalpyr, demo.coordmidpyr);
      break;
    case LibUtilities::ePrism:
      minvalandpthold = demo.FindLatticeEval(uhatslocal, demo.storage3dpri, demo.coordpri, demo.midptevalpri, demo.coordmidpri);
      break;
    default:
      cout<<"\n not implemented for this shape type yet!";
      exit(0);
	  
    }
  int sz = minvalandpthold.size();
  for(int l = 0; l < minvalandpthold.size(); l++)
    {
      cout<<" "<<minvalandpthold[l];
    }
  cout<<"\n\n";
  
  if(minvalandpthold[sz-1] < 0 && abs(minvalandpthold[sz-1])>tol)
    {
      return minvalandpthold;
    }	    
  return NullNekDouble1DArray;
}



Array<OneD, NekDouble> eval_sol(Array<OneD, NekDouble> x, Array<OneD, NekDouble> y, Array<OneD, NekDouble> z)
{
  boost::ignore_unused(x,y,z);
  int totPoints = x.size();
  Array<OneD, NekDouble> sol(totPoints);
  //get solution array

  for (int i = 0; i < totPoints; ++i)
    {
      if (dimension == 3)
	{
	  sol[i] = floor(0.5*floor(0.5*(floor(0.5*(x[i] >= -0.8) + 0.5*(x[i]<=-0.5))) +  0.5*(floor(0.5*( y[i]>=0) + 0.5*(y[i]<=0.3)))) + 0.5*floor(0.5*(z[i] <= -0.1) + 0.5*(z[i] >=-0.4)))*cos(3*M_PI*(x[i]))*sin(3*M_PI*(y[i]));
	  //(pow(x[i] +0.305082 , 2) + pow(y[i]+0.20068 ,2)) + pow(z[i]+0.59864,2)-0.0001;
	}
      else
	{
	  cout<<"\n This demo only tests 2-D and 3-D root-fnding\n";
	  exit(0);
	}
    }
  return sol;
}


void Do_optimize(StdExpansion *exp, Array<OneD, NekDouble> &uhats)
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
  tols.push_back(1e-9);

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
  Array<OneD, NekDouble> Vtmp( Eorth->GetNcoeffs());

  Timer t;
  while (counter <= niter)
    {
      NekDouble roots1dtime = 0.0, roots2dtime = 0.0, roots3dtime = 0.0 ;
      pqval = inf;
      utemp = d.back();

	  
      t.Start();
      //pq = @(xx) cfun(xx) .* (constraints(1,jj).*Vc(xx)*uhats);  

      project_edges(utemp, Pf, exp);

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

      //if(verboseflag)
      //{
      //    cout<<" \n optima:\n";
      // }
	  
      for(int k = 0; k < dim; k++)
	{
	  //  if(verboseflag)
	  //	{
	  //	  cout<<" "<<optima[k];
	  //	}
	  tmpcoord[k] = Array<OneD, NekDouble>(1, optima[k]);
	}
      // Vtmp is evaluation of basis fn at optima
	  
      Vtmp = Eorth->PhysEvaluateBasis(tmpcoord, storage3dhold, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

      demo.pq(utemp, tmpcoord, Vtmp, pqvalcoords, wsp1);

      // if(verboseflag)
      // 	{
      cout<<"\nat iter = "<<counter<<"  pqval = "<<pqvalcoords[0]<< " at "<<  optima[0]<<" "<< optima[1]<<" "<< optima[2];	  
      //	}
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
	  tmp = Eorth->PhysEvaluateBasis(xastarrofarr, storage3dhold, i);
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
  if(verboseflag)
    {
      cout<<"\nsphere_rotation took "<<counter<<"iterations";//
    }
  if(verboseflag)
    {
      //      cout<<"\n  startv = "<<startval<<" at startcoords("<<startcoordx<<","<<startcoordy<<" ," <<startcoordz<<")   \n GD iter stepszie ="<<demo.GetTol()<<" maxiters = "<<demo.GetMaxIter()<<" o = "<<Eorth->GetBasis(0)->GetNumModes()<<" p = "<<Eorth->GetBasis(0)->GetNumPoints()<<" eps = "<<demo.GetEps()<<"  avgiterGD = "<<itersGD2;//<<" total time taken post rootfinding per iter = "<< evalrootshold/(counter)<<" hex \n";;
    
      cout<<"\nAvg times per iter (1drootfinder,2drootfinder,3drootfinder) = : "  << roots1dtimehold<<", "<<roots2dtimehold<<", "<<roots3dtimehold<<" itergd2 = "<<itersGD2; ;//<< timeprojectedges<<", "<<timeprojectsurf<<",
    }
  uhats = d.back();
}

