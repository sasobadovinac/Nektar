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
void surfuhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz );

void call_setup_hex();

void call_setup_tet();

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
StdExpansion *E;
StdExpansion *E3seg;
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
  if(dimension < 3 )
    {
      cout<<"\n dimension should be 3, try using StdProjectPositivityPres1D or StdProjectPositivityPres2D for other dimensions\n\n";
      exit(0);
    }

  std::vector<int> order;
  std::vector<BasisType> btype(3, eNoBasisType);
  LibUtilities::ShapeType stype = E->DetShapeType();
  LibUtilities::PointsType pointsTypeCheb = LibUtilities::eGaussGaussChebyshev;

  for (int i = 0; i < dimension; ++i)
    {
      btype[i] = E->GetBasisType(i);
      order.push_back(E->GetBasisNumModes(i));
    }

  switch(E->DetShapeType())
    {
    case LibUtilities::eTetrahedron:
        
      numedges = 6;
      numsurfaces = 4;
      break;

    case LibUtilities::eHexahedron:
        
      numedges = 12;
      numsurfaces = 6;
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
	  sol[i] = pow((x[i]+0.6),2)+ pow((y[i]+0.6),2) + pow((z[i]-0.3),2) + demo.GetEps();
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
      storage3d = E->GetPhysEvaluateStorage();
      NekDouble setup, timeneg, timeverify;

      //      NekDouble elapref = 0;
 
      // // calculating reference time here:
      // Array<OneD, NekDouble> physder(phys.size()), physder0(phys.size()), uhatsnew(coeffs.size(), 1.0);

      // // demo.GetAvgNum() refers to the user defined parameter to get timing values for demo.GetAvgNum() number of runs 
      // for(int k = 0; k < demo.GetAvgNum(); k++)
      // 	{
      
      // 	  timer.Start();
      // 	  E->BwdTrans(uhatsnew, physder);
      // 	  E->PhysDeriv(physder, physder0);
      // 	  E->FwdTrans(physder0, uhatsnew); 
      // 	  timer.Stop();
      // 	  elapref += timer.TimePerTest(1);
      // 	}

      // elapref = elapref/demo.GetAvgNum();
      // cout<<"\n Ref time (first derivative uhats from f(x) uhats took) "<<(elapref)<<" s\n\n";

      //      int dimension = E->GetShapeDimension(); 
        
      vector<LibUtilities::BasisKey> tmpBasisKey = demo.GetBasisKey();       
      LibUtilities::PointsKey pkeycheb(E->GetBasis(0)->GetNumPoints(), pointsTypeCheb);
      
      LibUtilities::BasisKey bkeycheb(LibUtilities::eChebyshev, 3*(E->GetBasis(0)->GetNumModes()),  pkeycheb);
      
      E3seg = new StdSegExp(bkeycheb);
      if (E3seg == nullptr)
        {
	  return 1;
        }
   
      if(stype != LibUtilities::eTetrahedron)
	{
	  LibUtilities::BasisKey bkeycheb2 (LibUtilities::eChebyshev, (E->GetBasis(0)->GetNumModes()), pkeycheb);
	  Equad = new StdQuadExp(bkeycheb2, bkeycheb2);
	  if (Equad == nullptr)
	    {
	      return 1;
	    }
	  storage2dq = Equad->GetPhysEvaluateStorage();
	  // to-do : move this in demo constructor
	  demo.testcoord2dqqpts = demo.GetCoords(Equad);
	  demo.testcoord2dqqmidpts = demo.GetQuadratureMidCoords(Equad,  demo.testcoord2dqqpts);
	  demo.testcoord2dqlattice = demo.GetLatticeCoords(  demo.testcoord2dqqpts, demo.testcoord2dqqmidpts);

	  demo.interioreval2dqqmidpts = Equad->PhysEvaluateBasis(  demo.testcoord2dqqmidpts, storage2dq, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray );
	

	}
      if(stype != LibUtilities::eHexahedron)
	{
	  LibUtilities::PointsType pointsTypetriA = LibUtilities::eGaussLobattoLegendre;
	  LibUtilities::PointsKey pkeytri(E->GetBasis(0)->GetZ().size(), pointsTypetriA);
	  LibUtilities::BasisKey bkeytriB(LibUtilities::eOrtho_B,(E->GetBasis(0)->GetNumModes()), pkeytri);
	  LibUtilities::BasisKey bkeytriA(LibUtilities::eOrtho_A, (E->GetBasis(0)->GetNumModes()),pkeytri);
	  Etri = new StdTriExp(bkeytriA, bkeytriB);
	  if (Etri== nullptr)
	    {
	      return 1;
	    }
	  storage2dt = Etri->GetPhysEvaluateStorage();
	  // to-do : move this in demo constructor
	  demo.testcoord2dtqpts = demo.GetCoords(Etri);
	  demo.testcoord2dtqmidpts = demo.GetQuadratureMidCoords( Etri, demo.testcoord2dtqpts);
	  demo.testcoord2dtlattice = demo.GetLatticeCoords(demo.testcoord2dtqpts, demo.testcoord2dtqmidpts);
	  demo.interioreval2dtqmidpts = Etri->PhysEvaluateBasis(  demo.testcoord2dtqmidpts, storage2dt, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray   );
	  
	}
      
      demo.testcoord3dqpts  =  demo.GetCoords(E);
      demo.testcoord3dqmidpts  = demo.GetQuadratureMidCoords(E, demo.testcoord3dqpts);
      demo.testcoord3dlattice = demo.GetLatticeCoords(demo.testcoord3dqpts, demo.testcoord3dqmidpts);
      demo.interioreval3dqmidpts = E->PhysEvaluateBasis(demo.testcoord3dqmidpts, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray   );

      timer.Stop();
      setup = timer.TimePerTest(1);  

      timer.Start();
      
      int opt_need = Opt_needed(coeffs);            
      
      timer.Stop();
      timeneg = timer.TimePerTest(1);  
      cout<<"\nchecking for -vity took "<<timeneg<<"s\n";

      cout<<"need optimization";

      if(opt_need)
	{
	  timer.Start();	  
	  if(numedges == 12) //hex
	    {
	      call_setup_hex();
	    }
	  else if(numedges == 6) //tet
	    {
	      call_setup_tet();
	    }
	  timer.Stop();
	  setup += timer.TimePerTest(1);

	  cout<<"\nsetup phase took "<<setup<<"s";
	
	  NekDouble elap = 0;
	  timer.Start();
	  Do_optimize(coeffs);
	  timer.Stop();
	  elap += timer.TimePerTest(1); 
	  cout<<"\noptimizertook "<<elap<<"s";
	    
	  timer.Start();
	  opt_need = Opt_needed(coeffs, 1 );
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
	cout<<"\nchecking for -vity took "<<timer.TimePerTest(1)<<"s\n\n";
	
	cout<<"\noptimizer no need\n\n";
      }
      //Backward transform solution to get projected values
      E->BwdTrans(coeffs, phys);
      
    }
    
  //Calculate L_inf & L_2 error
  cout << "L infinity error: \t" << E->Linf(phys, sol) << endl;
  cout << "L 2 error: \t \t \t" << E->L2(phys, sol) << endl;
  cout<<"\n***********\n";
  return 0;
}

void call_setup_tet()
{
  Array<OneD, NekDouble> edgexyztemp =   E3seg->GetBasis(0)->GetZ();

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
  Array<OneD, NekDouble> edgexyztemp =   E3seg->GetBasis(0)->GetZ();

  int totszedges1d = edgexyztemp.size()*(E->GetNcoeffs());
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

            
  Vxm1ym1z =             E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);  
                
  //edge front right (x = 1) (y = -1)
  Vx1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1ym1z = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1ym1z = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);       
  Vx1ym1z = E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxx1ym1z, Vdyx1ym1z, Vdzx1ym1z);  


  //edge front top (y = -1) (z = 1)
  Vym1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxym1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyym1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzym1xz1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = edgexyztemp; 
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
  Vym1xz1 = E->PhysEvaluateBasis(edgeptsin,storage3d, Vdxym1xz1, Vdyym1xz1, Vdzym1xz1);  

                
  //edge front bot (y = -1) (z = -1)
  Vym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  Vym1xzm1 =             E->PhysEvaluateBasis( edgeptsin, storage3d, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);  
	    
	    
  // edge back left (y = 1), (x = -1)
  Vxm1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1y1z = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);      
  edgeptsin[2] =    edgexyztemp;


  Vxm1y1z = E->PhysEvaluateBasis(edgeptsin,storage3d, Vdxxm1y1z, Vdyxm1y1z, Vdzxm1y1z);  

  //edge back right (x = 1), (y = 1))
  Vx1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1y1z = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1y1z = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);                
         
  Vx1y1z  = E->PhysEvaluateBasis(edgeptsin,storage3d, Vdxx1y1z, Vdyx1y1z, Vdzx1y1z);  

                
  //edge back top ( y = 1) (z = 1)
  Vy1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxy1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyy1xz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzy1xz1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[0] = edgexyztemp;
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
  Vy1xz1 =             E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxy1xz1, Vdyy1xz1, Vdzy1xz1 );  

      
  //edge back bot (y = 1) (z = -1))
  Vy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
            
  Vy1xzm1 = E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);  

  // edge left bot (z = -1), (x = -1)
  Vxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  edgeptsin[1] = edgexyztemp;
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

  Vxm1yzm1 = E->PhysEvaluateBasis(edgeptsin,storage3d, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);  
                
  //edge left top (x = -1), (z = 1))
  Vxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
                
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);           
  Vxm1yz1 = E->PhysEvaluateBasis(edgeptsin,storage3d, Vdxxm1yz1, Vdyxm1yz1, Vdzxm1yz1);  

  //edge right bot ( z = -1) (x = 1)
  Vx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1yzm1 = Array<OneD, NekDouble>(totszedges1d);

  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
  edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);           
  Vx1yzm1 = E->PhysEvaluateBasis(edgeptsin, storage3d, Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);   

  //edge right top (z  1) (x  1))
  Vx1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdxx1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdyx1yz1 = Array<OneD, NekDouble>(totszedges1d);
  Vdzx1yz1 = Array<OneD, NekDouble>(totszedges1d);
                
  edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);

  Vx1yz1 = E->PhysEvaluateBasis(edgeptsin,storage3d, Vdxx1yz1, Vdyx1yz1, Vdzx1yz1);   

              
  int totszsurf2d = (demo.testcoord2dqqpts[0].size())*E->GetNcoeffs();
            
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
  Vxyzm1 =             E->PhysEvaluateBasis(surfptsintemp,storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

  //surface right x = 1
  Vx1yz =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[0] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
  surfptsintemp[1] =  surfptsin[0];
  surfptsintemp[2] =  surfptsin[1];
  Vx1yz = E->PhysEvaluateBasis(surfptsintemp,storage3d,NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   


  //surface top z = 1
  Vxyz1 =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[1] = surfptsin[1];
  surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
  Vxyz1 = E->PhysEvaluateBasis(surfptsintemp, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

  //surface left x = -1
  Vxm1yz =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[0] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
  surfptsintemp[1] = surfptsin[0];
  surfptsintemp[2] = surfptsin[1];
  Vxm1yz = E->PhysEvaluateBasis(surfptsintemp,storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

  //surface front y = -1
  Vxym1z =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[1] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[2] = surfptsin[1];
  Vxym1z = E->PhysEvaluateBasis(surfptsintemp,storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   
      
  //surface back y = 1
  Vxy1z =  Array<OneD, NekDouble>(totszsurf2d);
  surfptsintemp[1] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
  surfptsintemp[0] = surfptsin[0];
  surfptsintemp[2] = surfptsin[1];
  Vxy1z = E->PhysEvaluateBasis(surfptsintemp,storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
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
      surfuhats(uhats, ret[0], Vxyzm1);

      //right surface x = 1
      surfuhats(uhats, ret[1], Vx1yz);

      //surface top z = 1
      surfuhats(uhats, ret[2], Vxyz1);
            
      //left surface x = -1
      surfuhats(uhats, ret[3], Vxm1yz);
       
      //surface front y = -1
      surfuhats(uhats, ret[4], Vxym1z);
       
      //surface back y = 1
      surfuhats(uhats, ret[5], Vxy1z);
    }
  else if(numsurfaces == 4) //tets
    {

      //surface bot z = -1 restriction on GD: x+y = 0 (ABC)
      surfuhats(uhats, ret[0],  Vxyzm1tet);

      //surface left x = -1 restriction on GD: y + z = 0 (DAC)
      surfuhats(uhats, ret[1], Vxm1ypz0tet);

      //surf front y = -1   restriction on GD: x + z = 0 (DAB)
      surfuhats(uhats, ret[2], Vxpz0ym1tet);

      //surf DCB restriction  on GD: (x + y + z = -1)
      surfuhats(uhats, ret[3], Vxpypzm1tet);

      
    }
}

void surfuhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz)
{
  int modes = uhats.size();

  Array<OneD, NekDouble> temp(uhats.size());  
  int totpts;
  if(E->DetShapeType() == LibUtilities::eHexahedron)
    {
      totpts = Equad->GetTotPoints();
    }
  else if (E->DetShapeType() == LibUtilities::eTetrahedron)
    {
      totpts = Etri->GetTotPoints();

    }
  Array<OneD, NekDouble> vals(totpts);    
  //Vxyz*uhats -> project to -> E3tri
  for(int k = 0; k < totpts; k++)
    {
      
      Vmath::Vmul(uhats.size(), &Vxyz[k], totpts, &uhats[0], 1, &temp[0], 1);
      
      vals[k]  = Vmath::Vsum(modes, temp, 1);
    }
  //    vals = demo.blasmatvec(Vxyz, uhats, totpts, modes);
  if(E->DetShapeType() == LibUtilities::eHexahedron)
    {
      Equad->FwdTrans(vals, ret);
    }
  else if (E->DetShapeType() == LibUtilities::eTetrahedron)
    {
      Etri->FwdTrans(vals, ret);

    }

}



// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret)
{
  if(numedges == 12) // hex
    {

      // edge front left (x = -1) (y = -1)
      edgederpquhats(uhats, ret[0], E->GetNcoeffs(), Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);                     
      //edge front right (x = 1) (y = -1)
      edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vx1ym1z, Vdxx1ym1z, Vdyx1ym1z, Vdzx1ym1z);
            
      //edge front top (y = -1) (z = 1)
      edgederpquhats(uhats, ret[2], E->GetNcoeffs(),  Vym1xz1, Vdxym1xz1, Vdyym1xz1, Vdzym1xz1);
    
      //edge front bot (y = -1) (z = -1)
      edgederpquhats(uhats, ret[3], E->GetNcoeffs(),  Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);
     
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
  else if(numedges == 6) //tet
    {
      // edge front left (AD) (x = -1) (y = -1)
      edgederpquhats(uhats, ret[0], E->GetNcoeffs(),  Vxm1ym1ztet, Vdxxm1ym1ztet, Vdyxm1ym1ztet, Vdzxm1ym1ztet);

      //edge front hypt (DB) (y = -1) (z = -x)
      edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vym1xmztet, Vdxym1xmztet, Vdyym1xmztet, Vdzym1xmztet);

      //edge front bot (AB) (y = -1) (z = -1)
      edgederpquhats(uhats, ret[2], E->GetNcoeffs(), Vym1xzm1tet, Vdxym1xzm1tet, Vdyym1xzm1tet, Vdzym1xzm1tet );

      //edge left hypt (DC) ( x = -1) (z = -y)
      edgederpquhats(uhats, ret[3], E->GetNcoeffs(), Vxm1ymztet, Vdxxm1ymztet, Vdyxm1ymztet, Vdzxm1ymztet);

      // edge bot diag (BC) (z = -1) (y = -x)
      edgederpquhats(uhats, ret[4], E->GetNcoeffs(),  Vxmyzm1tet, Vdxxmyzm1tet, Vdyxmyzm1tet, Vdzxmyzm1tet);

      //edge CA bot left (x = -1) (z = -1)
      edgederpquhats(uhats, ret[5], E->GetNcoeffs(),Vxm1yzm1tet, Vdxxm1yzm1tet,  Vdyxm1yzm1tet,  Vdzxm1yzm1tet);

	
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
    
  E3seg->FwdTrans(pqeval, ret);
  
}

Array<OneD, NekDouble>  call_find_roots(Array<OneD,  NekDouble> &uhats , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats, NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime)
{
  int dimension = E->GetShapeDimension(); 
  NekDouble dum;
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
  if(numedges == 6) // tet
    {
      for(int ii = 0; ii < numedges; ii++)
        {
	  t.Start();
	  rethold  = (demo.find_roots(uhatsedges[ii], E, NullNekDoubleArrayOfArray, dum,  0, 0, 0, 0));
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
  else if(numedges == 12) //hex
    {
      for(int ii = 0; ii < numedges; ii++)
	{
	  
	  t.Start();
	  rethold  = (demo.find_roots(uhatsedges[ii], E, NullNekDoubleArrayOfArray, dum,  0, 0, 0, 0));
	  t.Stop();

	  roots1dtime +=  t.TimePerTest(1);
	  
	  for(int p = 0; p < rethold[0].size(); p++)
	    {
	      switch(ii)
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
  if(numsurfaces == 4)
    {
      NekDouble avgiterGDhold;

      // call 2D rootfinder on each surface:
      for(int ii = 0; ii < numsurfaces; ii++)
	{
	  t.Start();
	  rethold = demo.find_roots(surfaceuhats[ii], Etri, storage2dt,  avgiterGDhold, 0, 1, 0, 1);
	  // last arg = 1 coz all 2d ele are triangles

	  t.Stop();
	  roots2dtime +=  t.TimePerTest(1);
	  avgiterGD += avgiterGDhold;

	  switch(ii)
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
  else if(numsurfaces == 6)
    {
      NekDouble avgiterGDhold;
      for(int ii = 0; ii < numsurfaces; ii++)
	{
	  t.Start();   
	  rethold = demo.find_roots(surfaceuhats[ii], Equad, storage2dq,  avgiterGDhold, 0, 1, 0 , 0);
	  t.Stop();

	  roots2dtime +=  t.TimePerTest(1);
	  avgiterGD += avgiterGDhold;
	      
	  if(ii == 0) // surf bot (z = -1)
	    {
	      retall[0].push_back( (rethold[0][0]));
	      retall[1].push_back( (rethold[1][0]));
	      retall[2].push_back( (-1.0));
	    }
	  else if(ii == 1) //surf right (x = 1)
	    {
	      retall[0].push_back  (1.0);
	      retall[1].push_back (rethold[0][0]);
	      retall[2].push_back (rethold[1][0]);
	    
	    }
	  else if(ii == 2) //surf top (z = 1)
	    {
	      retall[0].push_back  (rethold[0][0]);
	      retall[1].push_back (rethold[1][0]);
	      retall[2].push_back (1.0);
	    }
	  else if(ii == 3) //surf left (x = -1)
	    {
	      retall[0].push_back  (-1.0);
	      retall[1].push_back (rethold[0][0]);
	      retall[2].push_back (rethold[1][0]);
	    }
	  else if(ii == 4) //surf front (y = -1)           
	    {
	      retall[0].push_back  (rethold[0][0]);
	      retall[1].push_back (-1.0);
	      retall[2].push_back (rethold[1][0]);
	    }
	  else if(ii == 5) //surf back (y = 1)
	    {
	      retall[0].push_back  (rethold[0][0]);
	      retall[1].push_back (1.0);
	      retall[2].push_back (rethold[1][0]);
	    }
	    
	}
      
    }                     
    
  //find interior roots:
      
  t.Start();
  rethold =  (demo.find_roots(uhats, E,storage3d,  avgiterGDhold, 0, 0, 1, 0)) ;
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
  Array<OneD, NekDouble> evalroots = E->PhysEvaluateBasis(tmpcoord, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
     
  Array<OneD, NekDouble> minvandx(dimension+1);
  Array<OneD, NekDouble> tmparr(retall[0].size());


  demo.pq(uhats, tmpcoord, evalroots, minvandx, tmparr);

  // cout<<"\n vals at roots:\n";
  // for(int k = 0; k < tmparr.size(); k++)
  //   cout<<" "<<tmparr[k]<<" ";
  // cout<<"\n";

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

  //    int NC = 1; //number of constraints, only positivity for now
  tols.push_back(1e-12);

  Array<OneD, NekDouble> optima(dim), pqvalcoords(dim+1), xastarr(dim), utemp(N1), wsp1(1); 
  Array<OneD, Array<OneD, NekDouble> > surfaceuhats  (numsurfaces), tmpcoord(dim), Pf(numedges);
  for(int k= 0; k < numedges; k++)
    {
      Pf[k] = Array<OneD, NekDouble>(E3seg->GetNcoeffs());
    }

  for(int k = 0; k < numsurfaces; k++)
    {
      if(E->DetShapeType() == LibUtilities::eHexahedron)
	{
	  surfaceuhats[k] = Array<OneD, NekDouble>(Equad->GetNcoeffs());
	}
      else if( E->DetShapeType() == LibUtilities::eTetrahedron)
	{
	  surfaceuhats[k] = Array<OneD, NekDouble>(Etri->GetNcoeffs());
	}
    }

  NekDouble pqval, timeprojectsurf = 0.0, timeprojectedges = 0.0;
  NekDouble avgiterGD = 0.0, avgiterhold = 0.0, roots1dtimehold = 0.0, roots2dtimehold = 0.0, roots3dtimehold = 0.0;
  Array<OneD, NekDouble> Vtmp(  E->GetNcoeffs());


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
    
	  // Vtmp is evaluation of basis fn at optima
	  for(int k = 0; k < dim; k++)
	    {
	      tmpcoord[k] = Array<OneD, NekDouble>(1, optima[k]);
	    }
	
	  Vtmp = E->PhysEvaluateBasis(tmpcoord, storage3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

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

      counter = counter + 1;
    
    
    }
  roots1dtimehold = roots1dtimehold/(counter-1);
  roots2dtimehold = roots2dtimehold/(counter-1);
  roots3dtimehold = roots3dtimehold/(counter-1);
  timeprojectedges = timeprojectedges/(counter);
  timeprojectsurf = timeprojectsurf/(counter);    
  itersGD2 = (avgiterGD)/(counter);
  iterstaken = counter;
  cout<<"\nsphere_rotation took "<<counter<<"iterations";//

  //cout<<"\n  startv = "<<startval<<" at startcoords("<<startcoordx<<","<<startcoordy<<" ," <<startcoordz<<")    GD iter stepszie ="<<demo.GetTol()<<" maxiters = "<<demo.GetMaxIter()<<" o = "<<E->GetBasis(0)->GetNumModes()<<" p = "<<E->GetBasis(0)->GetNumPoints()<<" eps = "<<demo.GetEps()<<"  avgiterGD = "<<itersGD2;//<<" total time taken post rootfinding per iter = "<< evalrootshold/(counter)<<" hex \n";;

  cout<<"\nAvg times per iter (1drootfinder,2drootfinder,3drootfinder) = : "  << roots1dtimehold<<", "<<roots2dtimehold<<", "<<roots3dtimehold ;//<< timeprojectedges<<", "<<timeprojectsurf<<",

  uhats = d.back();
}
