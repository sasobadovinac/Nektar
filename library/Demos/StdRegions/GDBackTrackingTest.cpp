///////////////////////////////////////////////////////////////////////////////
//
// File: RootFinderTest.cpp
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
// Description: Demo to test functionality of 2D rootfinder (bt line search)
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "StdDemoSupport.hpp"
#include <LibUtilities/BasicUtils/Timer.h>

namespace po = boost::program_options;

Array<OneD,NekDouble>  callGD(Array<OneD,  NekDouble> &uhats1 , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats , NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime);

Array<OneD, NekDouble> find_rough_min(Array<OneD, NekDouble> uhats);


//string funcdef;
Array<OneD, Array<OneD, NekDouble> > storage2d, storage3d;

// Colleague matrix
int numedges, numsurfaces;
//Array<OneD, Array<OneD, NekDouble> > C;
StdExpansion *E, *Eorth;
StdExpansion *Equad;
StdExpansion *Etri;

void call_setup_quad(StdExpansion *E);
void call_setup_hex(StdExpansion *E);
void call_setup_tet(StdExpansion *E);
void call_setup_pyr(StdExpansion *E);
void call_setup_tri(StdExpansion *E);  

int dimension ;

DemoSupport demo;
LibUtilities::ShapeType stype;
Array<OneD, Array<OneD, NekDouble> > edgeptsin;

int main(int argc, char *argv[])
{
  demo.GetOptions().add_options()("optm,z", "positivity preserving optimizer");
  demo.ParseArguments(argc, argv);

  po::variables_map vm = demo.GetVariableMap();

  //only for 2D
  E = demo.CreateStdExpansion();

  stype = E->DetShapeType();

  if (E == nullptr)
    {
      return 1;
    }
  dimension = E->GetShapeDimension();
std::vector<int> order;
  std::vector<BasisType> btype(3, eNoBasisType);

  for (int i = 0; i < dimension; ++i)
    {
      btype[i] = E->GetBasisType(i);
      order.push_back(E->GetBasisNumModes(i));
    }

  switch(E->DetShapeType())
    {
    case LibUtilities::eQuadrilateral:
      call_setup_quad(E);
      numsurfaces = 1;
      numedges = 4;
      break;
    case LibUtilities::eTriangle:
      call_setup_tri(E);

      numsurfaces = 1;
      numedges = 3;
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
 if(dimension == 1 )
    {
      cout<<"\n dimension should be 2 or 3 for this demo!";
      exit(0);
    }
  else if(dimension == 2)
    {
      storage2d = E->GetPhysEvaluateStorage(); 
    }
  else
    {
      storage3d = E->GetPhysEvaluateStorage();
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
  demo.funcdef = "-2x + (x+0.6)^3 + y^2 - 0.2";

  //get solution array
  for (int i = 0; i < totPoints; ++i)
    {
      if(dimension ==2) 
        {
	  // f3
	  sol[i] = -2*x[i] + pow((x[i] + 0.6),3) + (pow(y[i],2) - 0.2);
	  
	}
      else if (dimension == 3)
	{
	  sol[i] =-2*x[i] + pow((x[i] + 0.6),3) + (pow(y[i],2) - 0.2);
	  //floor(x[i]<=0 && y[i] <=0)+0.1;
	  //-1*(sin((x[i]-0.1)+0.5*M_PI))*cos((y[i]-0.2))+0.95; 
	}
      else
	{
	  cout<<"\n This demo only tests 2-D and 3-D root-fnding\n";
	  exit(0);
	}
    }
  Array<OneD, NekDouble> phys(totPoints);
  Array<OneD, NekDouble> coeffs((unsigned) E->GetNcoeffs());
 
  //Project onto expansion
  E->FwdTrans(sol, coeffs);

  //Backward transform solution to get projected values
  E->BwdTrans(coeffs, phys);
  NekDouble v1;
  Array<OneD, Array<OneD, NekDouble> > dummy, surf(1);
  surf[0] = Array<OneD, NekDouble>(coeffs.size());
  Vmath::Vcopy(coeffs.size(), &coeffs[0], 1, &surf[0][0], 1); 

  Array<OneD, NekDouble> optima = callGD(coeffs, v1, dummy, surf,  v1, v1, v1, v1);

}


Array<OneD,  NekDouble> callGD(Array<OneD,  NekDouble> &coeff , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats, NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime)
{
  
  boost::ignore_unused(surfaceuhats, roots2dtime, roots3dtime, roots1dtime, minv, uhatsedges);
  int dimension = E->GetShapeDimension();
  Array<OneD, NekDouble> coords(dimension);
  Array<OneD, Array<OneD, NekDouble > >  rethold(dimension), holdstorage;
  // for(NekDouble  g = 0.1; g < 0.9; g+=0.1)
  //   {
  // 	for(NekDouble c = 0.1; c< 0.9; c+=0.1)
  // 	  {
  Array<OneD, NekDouble>  startA = find_rough_min(coeff);

  demo.SetStartArr(startA);
  // demo.Setchold(c);
  // demo.Setgamhold(g);
  // cout<<"\n c = "<<c<<" gam"<<g;
  if(stype == LibUtilities::ePoint)
    {
      cout<<"\n This demo works only for dim = 2 and 3!";
      exit(0);
    }
	    
  if(numsurfaces == 1)
    {
      holdstorage = storage2d;
      if(stype == 4)	  //ifstarting pt = "<<xnew[0]<<","<<xnew[1]<<" N = 5"; quad
	{
      demo.steepestgradient_descent2D(coeff, E, rethold, avgiterGD, 1e-7); 
	  //, demo.testcoord2dqqpts, demo.testcoord2dqqmidpts,demo. interioreval2dqqmidpts, demo.testcoord2dqlattice);
	}
      else if(stype == 3)//else tri
	{
      demo.steepestgradient_descent2D(coeff, E, rethold, avgiterGD, 1e-7); 
	}
    }
  else //3d
    {
      if(stype == LibUtilities::eHexahedron)
	{
	  storage3d = demo.storage3dhex;
	  demo.steepestgradientdescent3D(coeff, E, storage3d, demo.coordhex, demo.coordmidhex, demo.midptevalhex, rethold, avgiterGD);  

	}
      else if( stype == LibUtilities::eTetrahedron)
	{
	  storage3d = demo.storage3dtet;
	 
	  demo.steepestgradientdescent3D(coeff, E, storage3d, demo.coordtet, demo.coordmidtet, demo.midptevaltet, rethold, avgiterGD);  
	}
      else if( stype == LibUtilities::ePyramid)    //if hex
	{
	  storage3d = demo.storage3dpyr;
	  demo.steepestgradientdescent3D(coeff, E, storage3d, demo.coordpyr, demo.coordmidpyr, demo.midptevalpyr, rethold, avgiterGD);  
    
	}
       holdstorage = storage3d;

    }
  cout<<"\n iters taken = "<<avgiterGD;
  cout<<"\n expected ="<< startA[dimension]<<" rethold.size()="<<rethold.size()<<" rethold[0].size()="<<rethold[0][0]<<", "<<rethold[1][0]<<",";//<<rethold[2][0]<<"\n\n";

  coords[0] = rethold[0][0];
  coords[1] = rethold[1][0];
    
  if(dimension == 3)
    coords[2] = rethold[2][0];

  Array<OneD, NekDouble> evalbasis(coeff.size()), tmp(1);
  evalbasis = E->PhysEvaluateBasis(rethold, holdstorage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
  demo.pq(coeff, rethold, evalbasis, NullNekDouble1DArray, tmp);
  Array<OneD, NekDouble> s1(E->GetTotPoints(),0.0), p1(E->GetTotPoints(),0.0);
  cout<<"\n Min value on dense lattice = "<<startA[dimension];
  if(startA[dimension] > tmp[0] )
    {
      startA[dimension] = tmp[0];
    }
	    
  p1[0] = startA[dimension];;
  s1[0] = tmp[0];
	    
  cout<<"\n expected ="<< startA[dimension]<<" got"<<tmp[0]<<" \n";
  NekDouble eLinf, eL2;
	    
  eLinf = E->Linf(s1,p1);
  eL2 = E->L2(s1,p1);
	    
  cout << "L infinity error : " << scientific << eLinf << endl;
  cout << "L 2 error        : " << scientific << eL2 << endl<<" *****\n";
	 
  return coords;	
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
      E = new StdHexExp(b0, b1, b2);
    }
  else
    {
      E = exp;
    }
  demo.storage3dhex =  E->GetPhysEvaluateStorage();
  
  demo.coordhex = demo.GetCoords(E);
  demo.coordmidhex = demo.GetQuadratureMidCoords(demo.coordhex);
  demo.coordlatticehex = demo.GetLatticeCoords(demo.coordhex, demo.coordmidhex);
  
  demo.midptevalhex = E->PhysEvaluateBasis(demo.coordmidhex, demo.storage3dhex, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
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
      E = new StdPyrExp(b0, b1, b2);
    }
  else
    {
      E = E;
    }
  demo.storage3dpyr = E->GetPhysEvaluateStorage();
  demo.coordpyr = demo.GetCoords(E);
  demo.coordmidpyr = demo.GetQuadratureMidCoords(demo.coordpyr);
  demo.coordlatticepyr = demo.GetLatticeCoords(demo.coordpyr, demo.coordmidpyr);
  
  demo.midptevalpyr = E->PhysEvaluateBasis(demo.coordmidpyr, demo.storage3dpyr, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
  

}

Array<OneD, NekDouble> find_rough_min(Array<OneD, NekDouble> uhats)
{
  int dim = E->GetShapeDimension();
  Array<OneD, NekDouble> roughmin(dim+1);
  Array<OneD, Array<OneD, NekDouble> > tempdense (dim), coeffhold(dim);

  for(int k = 0; k < dim; k++)
    {
      if(dim == 2)
	tempdense[k] = Array<OneD, NekDouble>(11*11);
      else
	tempdense[k] = Array<OneD, NekDouble>(11*11*11);
      coeffhold[k] = Array<OneD, NekDouble>(1);

    }
  vector<NekDouble> idxx, idxy, idxz;
  int ct = 0, flg = 1;
  double inf = numeric_limits<double>::infinity();
  Array<OneD, NekDouble> holdpq(1);

  NekDouble allmin = inf;
  Array<OneD, Array<OneD, NekDouble> > storage_t;

  if(dim == 2)
    {
      storage_t = storage2d;

    }
  else if(stype == LibUtilities::eHexahedron || stype == LibUtilities::eTetrahedron)
    {
      storage_t = storage3d;
    }
  else
    {
      cout<<"\n shape type not impl yet!\n";
      exit(0);
    }
    
  if(dim == 3)
    {
    
      for(int y  = 0; y < 11; y++)
	{
	  NekDouble valtmpx = ( 1.0*y)/5 - 1.0;
	    
	  for(int u = 0; u < 11; u++)
	    {
	      NekDouble valtmpy = ( 1.0*u)/5 - 1.0;
	      for(int k = 0; k < 11; k++)
		{
		  tempdense[0][ct] = valtmpx;
		  coeffhold[0][0] = valtmpx;
		  tempdense[1][ct] = valtmpy;
		  coeffhold[1][0] = valtmpy;
		    
		  tempdense[2][ct] =(1.0*k)/5 - 1.0;
		  coeffhold[2][0] =(1.0*k)/5 - 1.0;

		  Array<OneD, NekDouble > tmp  = E->PhysEvaluateBasis(coeffhold, storage_t, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray  );
		  demo.pq(uhats, coeffhold, tmp, NullNekDouble1DArray, holdpq);
		  if(holdpq[0]<=allmin )
		    {
		      if(abs(holdpq[0] - allmin)<1e-10)
			{
			  idxx.push_back(tempdense[0][ct]);
			  idxy.push_back(tempdense[1][ct]);
			  idxz.push_back(tempdense[2][ct]);

			}
		      else
			{
			  idxx.clear();
			  idxy.clear();
			  idxz.clear();
			  allmin = holdpq[0];

			  idxx.push_back(tempdense[0][ct]);
			  idxy.push_back(tempdense[1][ct]);
			  idxz.push_back(tempdense[2][ct]);
			    
			}
		    }
	      
		  ct++;
		
		}
	    
	    }
	}
    }
  else if(dim == 2)
    {
    
      for(int y  = 0; y < 11; y++)
	{
	  NekDouble valtmpx = ( 1.0*y)/5 - 1.0;
	    
	  for(int u = 0; u < 11; u++)
	    {
	      NekDouble valtmpy = ( 1.0*u)/5 - 1.0;
		
	      if(stype == 3)
		{
		  flg = (valtmpx+valtmpy <= 1);
		}
	      else
		{
		  flg = 1;
		}
	      if(flg)
		{
		  tempdense[0][ct] = valtmpx;
		  coeffhold[0][0] = valtmpx;
		    
		  tempdense[1][ct] = valtmpy;
		  coeffhold[1][0] = valtmpy;
		}

	      Array<OneD, NekDouble > tmp  = E->PhysEvaluateBasis(coeffhold, storage_t, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray  );
	      demo.pq(uhats, coeffhold, tmp, NullNekDouble1DArray, holdpq);
		
	      if(holdpq[0]<=allmin )
		{
		  if(abs(holdpq[0] - allmin)<1e-10)
		    {
		      idxx.push_back(tempdense[0][ct]);
		      idxy.push_back(tempdense[1][ct]);
		    }
		  else
		    {
		      idxx.clear();
		      idxy.clear();
		      allmin = holdpq[0];
		      idxx.push_back(tempdense[0][ct]);
		      idxy.push_back(tempdense[1][ct]);
			
		    }
		    
		}
	      ct++;
		
	    }
	}
    }
    
  roughmin[0] = idxx[0];
  roughmin[1] = idxy[0];
  if(dim == 3)
    {
      roughmin[2] = idxz[0];
    }
  roughmin[dim] = allmin;

  return roughmin;

}
