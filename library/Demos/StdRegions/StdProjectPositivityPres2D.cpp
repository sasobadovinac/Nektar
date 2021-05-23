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

//declare filter call
void Do_optimize(Array<OneD, NekDouble> &uhats);

Array<OneD,NekDouble>  call_find_roots(Array<OneD,  NekDouble> &uhats1 , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats , NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime);

// for 2D elements to get uhats at edges
void project_edges( Array<OneD, NekDouble> uhats,Array<OneD, Array<OneD, NekDouble> >&ret );


int Opt_needed(Array<OneD, NekDouble> uhats, Array<OneD, NekDouble> &evallattice, int flag ); 

Array<OneD, NekDouble> eval_sol(Array<OneD, NekDouble> x, Array<OneD, NekDouble> y, Array<OneD, NekDouble> z);
void call_setup_quad();
void call_setup_tri();

string funcdef;
NekDouble  itersGD2 ;

int iterstaken, verboseflag;
NekDouble  startcoordx, startcoordy, startcoordz, startval;
Array<OneD, Array<OneD, NekDouble> > edgexy(2);
Array<OneD, Array<OneD, NekDouble> > storage2dt, storage2dq, storage2d;

int numedges, numsurfaces;
StdExpansion *E;
StdExpansion *Eorth;

Array<OneD, NekDouble> Vxm1;
Array<OneD, NekDouble> Vdxxm1;
Array<OneD, NekDouble> Vdyxm1;
Array<OneD, NekDouble> Vx1 ;
Array<OneD, NekDouble> Vdyx1;
Array<OneD, NekDouble> Vdxx1;
Array<OneD, NekDouble> Vdxy1;
Array<OneD, NekDouble> Vdyy1;
Array<OneD, NekDouble> Vy1;
Array<OneD, NekDouble> Vym1;
Array<OneD, NekDouble> Vdxym1;
Array<OneD, NekDouble> Vdyym1;

// for triangle edges root-finding:
Array<OneD, NekDouble> Vxyhyp;
Array<OneD, NekDouble> Vdxxyhyp;
Array<OneD, NekDouble> Vdyxyhyp;

int dimension ;
            
DemoSupport demo;
LibUtilities::ShapeType stype;
Array<OneD, Array<OneD, NekDouble> > edgeptsin;

int main(int argc, char *argv[])
{
  
  demo.GetOptions().add_options()("optm,z", "positivity preserving optimizer");
  demo.ParseArguments(argc, argv);

  po::variables_map vm = demo.GetVariableMap();
  verboseflag = demo.GetVerbose(); 
  E = demo.CreateStdExpansion();
  
  if (E == nullptr)
    {
      return 1;
    }
  dimension = E->GetShapeDimension();
  storage2d = E->GetPhysEvaluateStorage();
  if(dimension != 2 )
    {
      cout<<"\n dimension should be 2, try using StdProjectPositivityPres1D or StdProjectPositivityPres3DD for other dimensions\n\n";
      exit(0);
    }

  std::vector<int> order;
  std::vector<BasisType> btype(3, eNoBasisType);
  stype = E->DetShapeType();
  
  for (int i = 0; i < dimension; ++i)
    {
      btype[i] = E->GetBasisType(i);
      order.push_back(E->GetBasisNumModes(i));
    }

    
  switch(E->DetShapeType())
    {
    case LibUtilities::eSegment:
        
      cout<<"\n This demo works only for quads and tri,  try using StdProjectPositivityPres1D or StdProjectPositivityPresDD for other shapes\n\n";
      exit(0);
      break;
    

    case LibUtilities::eTriangle:
      numedges = 3;
      numsurfaces = 1;
      edgeptsin = Array<OneD, Array<OneD, NekDouble> >(2);
      
      break;
    
    case LibUtilities::eQuadrilateral:

      numedges = 4;
      numsurfaces = 1;
      edgeptsin = Array<OneD, Array<OneD, NekDouble> >(2);
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

  sol = eval_sol(x, y, z);


  Array<OneD, NekDouble> phys(totPoints),  sollattice;
  Array<OneD, NekDouble> coeffs((unsigned) E->GetNcoeffs());
  NekDouble L2errlattice = 0, Linferrlattice = 0;
  
  
  //Project onto expansion
  E->FwdTrans(sol, coeffs);

  //Backward transform solution to get projected values
  E->BwdTrans(coeffs, phys);
  LibUtilities::Timer     timer;
  NekDouble elapsed       = 0.0;


  Array<OneD, NekDouble> retvalslattice, solhold;
  // check for -ve values and apply opt if necessary
  if (vm.count("optm"))
    {
      if(E->DetShapeType() == LibUtilities::eTriangle)
	{
	  call_setup_tri();
	  sollattice = eval_sol(demo.testcoord2dtlattice[0], demo.testcoord2dtlattice[1], NullNekDouble1DArray );
	  
	}
      if(E->DetShapeType() == LibUtilities::eQuadrilateral)
	{
	  call_setup_quad();
	  sollattice = eval_sol(demo.testcoord2dqlattice[0], demo.testcoord2dqlattice[1], NullNekDouble1DArray );

	}
      Array<OneD, NekDouble> retvalslatticesave;

      timer.Start();
      int found_negb4 = Opt_needed(coeffs, retvalslattice, 0);
      timer.Stop();
      
      elapsed += timer.TimePerTest(1);
      if(verboseflag)
	{
	  cout<<"\nchecking for -vity took "<<timer.TimePerTest(1)<<"s\n";
	}

      
      solhold = Array<OneD, NekDouble>(sollattice);
      
      retvalslatticesave = Array<OneD, NekDouble>(retvalslattice);
      int latticesz = retvalslattice.size();
      Vmath::Vsub(latticesz, retvalslatticesave, 1, solhold, 1, solhold, 1);
      Vmath::Vabs(latticesz, solhold, 1, solhold, 1);
      Linferrlattice = Vmath::Vmax(latticesz, solhold, 1);
      Vmath::Vmul(latticesz, solhold, 1, solhold, 1, solhold, 1);
      L2errlattice = sqrt(Vmath::Vsum(latticesz, solhold, 1) );
      cout << "Unconstrained L infinity error: \t" << Linferrlattice << endl;
      cout << "Unconstrained L2 error \t: \t" << L2errlattice << endl;
      
      if(found_negb4 == 1)
	{
	  int orthoflag = (E->GetBasisType(0) != LibUtilities::eOrtho_A);
	  Array<OneD, BasisType> btypearr(dimension);
	  btypearr[0] = LibUtilities::eOrtho_A;
	  timer.Start();
	  if(orthoflag)
	    {
	      if(numedges == 4 && demo.Eorthquad != nullptr)
		{
		  Eorth = demo.Eorthquad;
		  btypearr[1] = LibUtilities::eOrtho_A;
		  
		}
	      else //if (numedges == 3  && demo.Eorthtri != nullptr)
		{
		  Eorth = demo.Eorthtri;
		  btypearr[1] = LibUtilities::eOrtho_B;     
		}
	    }
	  else
	    {
	      Eorth = E;
	    }

	  storage2d = Eorth->GetPhysEvaluateStorage();

	  for(int k = 0; k <2; k++)
	    {
	      edgexy[k] = Array<OneD, NekDouble>(demo.E3seg->GetBasis(0)->GetNumPoints());
	      edgeptsin[k] =   Array<OneD, NekDouble>(edgexy[k]);
	    }
	  Array<OneD, NekDouble> edgexytemp =   demo.E3seg->GetBasis(0)->GetZ();
	  int totszedges = edgexytemp.size()*(Eorth->GetNcoeffs());

	  // left (x = -1)
	  Vxm1 = Array<OneD, NekDouble>(totszedges);
	  Vdyxm1 = Array<OneD, NekDouble>(totszedges);
	  Vdxxm1 = Array<OneD, NekDouble>(totszedges);

	  // left x = -1
	  edgexy[1] = edgexytemp;
	  edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
	  Vxm1 = Eorth->PhysEvaluateBasis(edgexy, storage2d, Vdxxm1, Vdyxm1, NullNekDouble1DArray);

	  // bot (y = -1)
	  Vym1  = Array<OneD, NekDouble>(totszedges);
	  Vdxym1  = Array<OneD, NekDouble>(totszedges);
	  Vdyym1  = Array<OneD, NekDouble>(totszedges);

	  // bot y = -1
	  edgexy[0] = edgexytemp;             
	  edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);

	  Vym1 = Eorth->PhysEvaluateBasis(edgexy,  storage2d, Vdxym1, Vdyym1, NullNekDouble1DArray);
      
	  if(numedges == 4)
	    {
	  
	      // right quad (x = 1)
	      Vx1 = Array<OneD, NekDouble>(totszedges);
	      Vdxx1 = Array<OneD, NekDouble>(totszedges);
	      Vdyx1 = Array<OneD, NekDouble>(totszedges);

	      // top quad (y = 1)
	      Vdxy1 = Array<OneD, NekDouble>(totszedges);
	      Vdyy1 = Array<OneD, NekDouble>(totszedges);
	      Vy1 = Array<OneD, NekDouble>(totszedges);

	      // right x = 1
	      edgexy[1] = edgexytemp; 
	      edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), 1.0); 
	      Vx1 = Eorth->PhysEvaluateBasis(edgexy, storage2d, Vdxx1, Vdyx1, NullNekDouble1DArray);
	  
	      //top y = 1
	      edgexy[0] = edgexytemp; 
	      edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), 1.0); 
	      Vy1 = Eorth->PhysEvaluateBasis(edgexy, storage2d,Vdxy1,Vdyy1, NullNekDouble1DArray);
	  
	    }
	  else if(numedges == 3)
	    {
	      // hypt tri (y = -x)
	      Vxyhyp = Array<OneD, NekDouble>(totszedges);
	      Vdxxyhyp = Array<OneD, NekDouble>(totszedges);
	      Vdyxyhyp = Array<OneD, NekDouble>(totszedges);
	  
	      edgexy[0] = edgexytemp;
	      Vmath::Smul(edgexy[0].size(), -1.0, edgexy[0], 1, edgexy[1] , 1);
	      Vxyhyp = Eorth->PhysEvaluateBasis(edgexy, storage2d,  Vdxxyhyp, Vdyxyhyp, NullNekDouble1DArray);
	    }

	  timer.Stop();
	  elapsed  = timer.TimePerTest(1);

	  if(verboseflag)
	    {
	      cout<<"\nsetup phase took "<<elapsed<<"s";
	    
	      cout<<"\nneed optimization. Negative found at "<<startcoordx<<","<<startcoordy<<" ="<<startval<<"\n";
	    }
	  
	  timer.Start();
	  if(orthoflag)
	    {
	      E->BwdTrans(coeffs, phys);
	      Eorth->FwdTrans(phys, coeffs);
	    }
	  Do_optimize(coeffs);
	  if(orthoflag)
	    {
	      Eorth->BwdTrans(coeffs, phys);
	      E->FwdTrans(phys, coeffs);
	      storage2d = E->GetPhysEvaluateStorage();
	      //	      demo.OrthoNormalize(E, coeffs, btypearr, 1);
	    }
	  timer.Stop();
	  elapsed += timer.TimePerTest(1);
	  if(verboseflag)
	    {
	      cout<<"\noptimizertook "<<timer.TimePerTest(1)<<"s";
	    }
	  //Backward transform solution to get projected values
	  //E->BwdTrans(coeffs, phys);
	  
	  int found_neg = Opt_needed(coeffs, retvalslattice, 1);

	  //L2 error between retvalslattice and eval_sol on lattice:
	  Vmath::Vcopy(sollattice.size(), sollattice, 1, solhold, 1);
	  Vmath::Vsub(latticesz, retvalslattice, 1, sollattice, 1,
		      sollattice, 1);
	  Vmath::Vabs(latticesz, sollattice, 1, sollattice, 1);
	  Linferrlattice = Vmath::Vmax(latticesz, sollattice, 1);
	  Vmath::Vmul(latticesz, sollattice, 1, sollattice, 1, sollattice, 1);
	  L2errlattice = sqrt(Vmath::Vsum(latticesz, sollattice, 1) );

	  cout<<"\n constrained:\n";
	  cout << "\nconstrained L infinity error: \t" << Linferrlattice << endl;
	  if (stype != ePoint) 
	    {
	      cout << "constrained L 2 error: \t \t " << L2errlattice<< endl;
	    }
	  
	  // check val now at the min ppreviously found:
	  Array<OneD, Array< OneD, NekDouble> > checkcoord(2);
	  Array<OneD, NekDouble> saveeval, holdv(coeffs.size());
	  checkcoord[0] = Array<OneD, NekDouble>(1, startcoordx);
	  checkcoord[1] = Array<OneD, NekDouble>(1, startcoordy);
	  
	  saveeval = E->PhysEvaluateBasis(checkcoord, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
	  
	  Vmath::Vmul(coeffs.size(), saveeval, 1, coeffs, 1, holdv, 1);
	  NekDouble valnow =  Vmath::Vsum(coeffs.size(), holdv, 1);
	  if(verboseflag)
	    {
	      cout<<"\n val now at "<<startcoordx<<" "<<startcoordy<<" "<<valnow;
	    }
	  if(valnow < 0 && abs(valnow)>1e-9)
	    {
	      found_neg = 1;
	    }
	  
	  if(found_neg)
	    {
	      timer.Stop();
	      elapsed += timer.TimePerTest(1);
	      if(verboseflag)
		{
		  cout<<"\nVerification took "<<timer.TimePerTest(1)<<"s";
		}
	      cout<<" \nfail";
	      //fail unit test:
	      
	      cout << "\nL infinity error: \t" << 1 << endl;
	      if (stype != ePoint)
		{
		  cout << "L 2 error: \t \t " << 1<< endl;
		}
	      return 0;
	    }
	  else
	    {
	      timer.Stop();
	      elapsed += timer.TimePerTest(1);
	      if(verboseflag)
		{
		  cout<<"\nVerification took "<<timer.TimePerTest(1)<<"s";
		}
	      cout<<"\npass";
		  
	      //pass  unit test:
	      
	      cout << "\nL infinity error: \t" << 0 << endl;
	      if (stype != ePoint)
		{
		  cout << "L 2 error: \t \t " << 0<< endl;
		}
	      return 0;

	    }
	  
	}
      else
	{
	  cout<<"\n Filter found no negative values\n";
	  // for unit test to pass
	  cout << "\nL infinity error: \t" << 0 << endl;
	  cout << "L 2 error: \t \t " << 0<< endl;
	  cout<<"\n***********\n";
	  return 0;
      
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
      cout << "L 2 error: \t \t " << E->L2(phys, sol) << endl;
    }
  
  return 0;
  
    
}

void project_edges( Array<OneD, NekDouble>uhats1,    Array<OneD, Array<OneD, NekDouble> >&ret)
{
    
  Array<OneD, NekDouble> wsp(uhats1.size());

  if(numedges == 4) // quad
    {
      // bot edge
      demo.edgederpquhats(uhats1, ret[0],  Vym1, Vdxym1, Vdyym1, NullNekDouble1DArray);
      // right edge
      demo.edgederpquhats(uhats1, ret[1], Vx1, Vdxx1, Vdyx1, NullNekDouble1DArray);
      // top edge
      demo.edgederpquhats(uhats1, ret[2],  Vy1, Vdxy1, Vdyy1, NullNekDouble1DArray);
      // left edge
      demo.edgederpquhats(uhats1, ret[3],   Vxm1, Vdxxm1, Vdyxm1, NullNekDouble1DArray);

    }
  if(numedges == 3) // tri
    {
      // bot edge  
      demo.edgederpquhats(uhats1, ret[0],   Vym1, Vdxym1, Vdyym1,NullNekDouble1DArray);

      // left edge
      demo.edgederpquhats(uhats1, ret[1],  Vxm1, Vdxxm1, Vdyxm1, NullNekDouble1DArray);

      // hypto
      demo.edgederpquhats(uhats1, ret[2],  Vxyhyp, Vdxxyhyp, Vdyxyhyp, NullNekDouble1DArray);
    }

}
 
Array<OneD, NekDouble> eval_sol(Array<OneD, NekDouble> x, Array<OneD, NekDouble> y, Array<OneD, NekDouble> z )
{
  boost::ignore_unused(x,y,z);
  int totPoints = x.size();
  Array<OneD, NekDouble> sol(totPoints);
  //get solution array
  for (int i = 0; i < totPoints; ++i)
    {
      if(dimension ==2) 
	{
	  // f3
	  sol[i] = floor(0.5*(floor(0.5*(x[i] >= -0.8) + 0.5*(x[i]<=-0.5))) +  0.5*(floor(0.5*( y[i]>=0) + 0.5*(y[i]<=0.3))))*cos(3*M_PI*(x[i]))*sin(3*M_PI*(y[i]));
	  //-1 to 1 golden fn :
	  //	    f2=@(x,y) floor(0.5*(floor(0.5*(x >= -0.8) + 0.5*(x<=-0.5))) +  0.5*(floor(0.5*( y>=0) + 0.5*(y<=0.3)))).*cos(3*pi*(x)).*sin(3*pi*(y))
	}	  
      else
	{
	  cout<<"\n This demo only tests 2-D and 3-D root-fnding\n";
	  exit(0);
	}
    }
  return sol;
}
 
int Opt_needed(Array<OneD, NekDouble> coeffs, Array<OneD, NekDouble> &valslattice, int flag) 
{

  int nq = E->GetTotPoints();
  Array<OneD, NekDouble> temp2(coeffs.size()), holdmidpteval;
  Array<OneD, Array<OneD, NekDouble> > holdlattice;
  NekDouble minv = numeric_limits<double>::infinity();
  int idxmin = -1;

  if(E->DetShapeType() == LibUtilities::eQuadrilateral)
    {
      holdmidpteval = demo.interioreval2dqqmidpts;
      holdlattice = demo.testcoord2dqlattice;
    }
  else if(E->DetShapeType() == LibUtilities::eTriangle)
    {
      holdmidpteval = demo.interioreval2dtqmidpts;
      holdlattice = demo.testcoord2dtlattice;
    } 
  
  Array<OneD, NekDouble> holdv(holdlattice[0].size());

  for (int i = 0; i<nq; i++)
    {
      Vmath::Vmul(coeffs.size(), &storage2d[0][i], nq, &coeffs[0], 1, &temp2[0], 1);  
      holdv[i] = Vmath::Vsum(coeffs.size(), temp2, 1);    
      if(holdv[i] < minv)
	{
	  minv = holdv[i];
	  idxmin = i;
	}
    }
  
  int ct = nq;
  int midptsize = holdlattice[0].size() - nq;
  
  for(int i = 0,  k = ct; k < holdlattice[0].size(); ++k, i++)
    {
      Vmath::Vmul(coeffs.size(), &holdmidpteval[i], midptsize, &coeffs[0], 1, &temp2[0], 1);   
      holdv[k] = Vmath::Vsum(coeffs.size(), temp2, 1);
      if(holdv[k] < minv)
	{
	  minv = holdv[k];
	  idxmin = k;
	}
    }
  if(verboseflag) 
    cout<<"\n minv = "<<minv<<" at "<< holdlattice[0][idxmin]<< " "<<holdlattice[1][idxmin];
 
  valslattice = Array<OneD, NekDouble>(holdv);
  
  if(minv < 0 && abs(minv) > 1e-10 && idxmin != -1)
    {
      if(flag == 0)
	{
	  startval = minv;
	  startcoordx = holdlattice[0][idxmin];
	  startcoordy = holdlattice[1][idxmin];
	}
      return 1;
    }
  return 0;
}

void call_setup_tri()
{
  storage2dt = E->GetPhysEvaluateStorage(); 
 
  demo.testcoord2dtqpts = demo.GetCoords(E);
  demo.testcoord2dtqmidpts = demo.GetQuadratureMidCoords( E, demo.testcoord2dtqpts);
  demo.testcoord2dtlattice = demo.GetLatticeCoords(demo.testcoord2dtqpts, demo.testcoord2dtqmidpts);
  demo.interioreval2dtqmidpts = E->PhysEvaluateBasis( demo.testcoord2dtqmidpts, storage2dt, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
  if(E->GetBasisType(0) != eOrtho_A )//demo.Eorthquad != nullptr)
    {
      LibUtilities::BasisKey bkeyorthA ( LibUtilities::eOrtho_A,
					 E->GetBasis(0)->GetNumModes(),
					 E->GetBasis(0)->GetPointsKey());
      LibUtilities::BasisKey bkeyorthB ( LibUtilities::eOrtho_B,
					 E->GetBasis(1)->GetNumModes(),
					 E->GetBasis(1)->GetPointsKey());
      demo.Eorthtri =  new StdTriExp(bkeyorthA, bkeyorthB);
	  
	    
      demo.interiorevalorth2dtqmidpts = demo.Eorthtri->PhysEvaluateBasis(demo.testcoord2dtqmidpts , demo.Eorthtri->GetPhysEvaluateStorage(), NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
    }

}

void call_setup_quad()
{
  storage2dq = E->GetPhysEvaluateStorage(); 
  
  demo.testcoord2dqqpts = demo.GetCoords(E);
  demo.testcoord2dqqmidpts = demo.GetQuadratureMidCoords( E, demo.testcoord2dqqpts);
  demo.testcoord2dqlattice = demo.GetLatticeCoords( demo.testcoord2dqqpts, demo.testcoord2dqqmidpts);
  demo.interioreval2dqqmidpts = E->PhysEvaluateBasis( demo.testcoord2dqqmidpts, storage2dq, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
  if(E->GetBasisType(0) != eOrtho_A )
    {
      LibUtilities::BasisKey bkeyorthA ( LibUtilities::eOrtho_A,
					 E->GetBasis(0)->GetNumModes(),
					 E->GetBasis(0)->GetPointsKey());
      LibUtilities::BasisKey bkeyorthB ( LibUtilities::eOrtho_A,
					 E->GetBasis(1)->GetNumModes(),
					 E->GetBasis(1)->GetPointsKey());
      demo.Eorthquad =  new StdQuadExp(bkeyorthA, bkeyorthB);
      
      
      demo.interiorevalorth2dqqmidpts = demo.Eorthquad->PhysEvaluateBasis(demo.testcoord2dqqmidpts , demo.Eorthquad->GetPhysEvaluateStorage(), NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
    }
  
}


//try replacing storage with the global var -> storage2d
Array<OneD,NekDouble>  call_find_roots(Array<OneD,  NekDouble> &uhats , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats, NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime)
{
  
  boost::ignore_unused(surfaceuhats, roots2dtime, roots3dtime, roots1dtime);
  int dimension = Eorth->GetShapeDimension(); 
  Array<OneD, NekDouble> temp(uhats.size());

  NekDouble  inf = numeric_limits<double>::infinity();
  Timer t;
  Array<OneD, Array<OneD, NekDouble > >  rethold(dimension);
  vector<vector<NekDouble> > retall(dimension);
  Array<OneD, Array<OneD, NekDouble> > tmpcoord(dimension);

  Array<OneD, NekDouble> ret(dimension);

  if(numedges == 4) // quad
    {
          
      for(int i = 0; i < numedges; i++)
	{

	  t.Start();
	  rethold  = demo.call_companion_rf(uhatsedges[i]);//, nullptr, NullNekDoubleArrayOfArray, dummy,  0, 0, 0, 0)) ;  
	  t.Stop();

	  roots1dtime +=  t.TimePerTest(1);
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
      // add 4 corners
      retall[0].push_back(-1);
      retall[1].push_back(-1);
      
      retall[0].push_back(1);
      retall[1].push_back(1);
       
      retall[0].push_back(-1);
      retall[1].push_back(1);
        
      retall[0].push_back(1);
      retall[1].push_back(-1);

    }
  else if(numedges == 3) //tri
    {
      
      for(int i = 0; i < numedges; i++)
	{
	  t.Start();

	  rethold  = demo.call_companion_rf(uhatsedges[i]);//(demo.find_roots(uhatsedges[i], nullptr, NullNekDoubleArrayOfArray, dummy,  0, 0, 0, 0)) ;
	  
	  t.Stop();
	  roots1dtime +=  t.TimePerTest(1);
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
		  //if(i == 1) // edge left (x = -1) 
		  retall[0].push_back(   (-1.0));
		  retall[1].push_back( (rethold[0][p]));  
		  
		  break;
		case 2: 
		  //if(i == 2) // edge hypt (y = -x) 
		  retall[1].push_back(   (-rethold[0][p]));
		  retall[0].push_back( (rethold[0][p]));  
		  
		  break;
		default:
		  break;
		}
	    }
	}
      // add 3 corners
      retall[0].push_back(-1);
      retall[1].push_back(-1);
   
      retall[0].push_back(-1);
      retall[1].push_back(1);
        
      retall[0].push_back(1);
      retall[1].push_back(-1);

    }
  if(numsurfaces == 1)
    {
      //	NekDouble avgiterGDhold;
      t.Start();   
      if(stype == 3)	  //if tri
	{
	  demo.steepestgradient_descent2Dtri(uhats, Eorth, storage2d, rethold, avgiterGD);  
	}
      else if(stype == 4)   //else quad
	{
	  demo.steepestgradient_descent2Dquad(uhats, Eorth, storage2d, rethold, avgiterGD);  
	}
      //rethold = demo.find_roots(uhats, Eorth, storage2d,  avgiterGD, 0, 1, 0 , 0);
	  
      t.Stop();
      roots2dtime +=  t.TimePerTest(1);

      if(rethold.size() > 0)
	{
	  if(rethold[0][0] < inf && rethold[1][0] < inf)
	    {
	      retall[0].push_back(rethold[0][0]);
	      retall[1].push_back(rethold[1][0]);
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
	}
      
      NekDouble tempmin = inf;
      Array<OneD, NekDouble> evalroots;
      if(stype == 3)
	evalroots = Eorth->PhysEvaluateBasis(tmpcoord, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
      else if(stype == 4)
	evalroots = Eorth->PhysEvaluateBasis(tmpcoord, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
	    
      Array<OneD, NekDouble> minvandx(dimension+1);
      Array<OneD, NekDouble> tmparr(retall[0].size());
      demo.pq(uhats, tmpcoord, evalroots, minvandx, tmparr);
	  
      tempmin = minvandx[0];
      ret[0] = minvandx[1];
      ret[1] = minvandx[2];
	  
      minv = tempmin;
    }
  return ret;

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
  Array<OneD,NekDouble> optima(dim);
  //    int NC = 1; //number of constraints, only positivity for now
  tols.push_back(1e-12);
    
  Array<OneD, NekDouble>  pqvalcoords(dim+1), xastarr(dim), utemp(N1), wsp1(1);   
    

  Array<OneD, Array<OneD, NekDouble> > surfaceuhats  (numsurfaces), tmpcoord(dim);
  Array<OneD, Array<OneD, NekDouble> > Pf(numedges);

  for(int k= 0; k < numedges; k++)
    {
      Pf[k] = Array<OneD, NekDouble>(demo.E3seg->GetNcoeffs()); 
    }

  for(int k = 0; k < numsurfaces; k++)
    {
      surfaceuhats[k] = Array<OneD, NekDouble>(Eorth->GetNcoeffs());
    }

  NekDouble pqval, timeprojectsurf = 0.0, timeprojectedges = 0.0;
  NekDouble avgiterGD = 0.0, avgiterhold = 0.0, roots1dtimehold = 0.0, roots2dtimehold = 0.0, roots3dtimehold = 0.0;
  Array<OneD, NekDouble> Vtmp(Eorth->GetNcoeffs());

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

	  optima = call_find_roots(utemp, avgiterhold, Pf, surfaceuhats, minv, roots1dtime, roots2dtime, roots3dtime);
	  roots1dtimehold += roots1dtime;
	  roots2dtimehold += roots2dtime;
	  roots3dtimehold += roots3dtime;
	  avgiterGD += avgiterhold;
	  // Vtmp is evaluation of basis fn at optima
	  for(int k = 0; k < dim; k++)
	    {
	      tmpcoord[k] = Array<OneD, NekDouble>(1, optima[k]);
	    }
	
	  Vtmp = Eorth->PhysEvaluateBasis(tmpcoord, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        
	}// end if(counter > 0)
      else 
	{
	  for(int k = 0; k < dim; k++)
	    {
	      tmpcoord[k] = Array<OneD, NekDouble>(1);
	    }

	  tmpcoord[0][0] = startcoordx;
	  tmpcoord[1][0] = startcoordy;

	  optima[0] = startcoordx;
	  optima[1] = startcoordy;
	  Vtmp = Eorth->PhysEvaluateBasis(tmpcoord, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        

	}
      demo.pq(utemp, tmpcoord, Vtmp, pqvalcoords, wsp1);

      if (pqvalcoords[0] < pqval)
	{
	  xastarr = optima;
	  pqval = pqvalcoords[0];
	}
      if(verboseflag)
	{
	  cout<<"\n pqval = "<<pqval<<" at "<<xastarr[0]<<" "<<xastarr[1]<<"\n";        
	}
        
      // If minimum is non-negative, we're done
      if (pqval >= -tols.at(0))
	{
	  break;
	}
        
      Array<OneD, NekDouble>  Vastsq (N1);
      Array<OneD, NekDouble>  Vast (N1);
      Array<OneD, Array<OneD, NekDouble> > xastarrofarr(2);
      xastarrofarr[0] = Array<OneD, NekDouble>(1, xastarr[0]);
      xastarrofarr[1] = Array<OneD, NekDouble>(1, xastarr[1]);

      Array<OneD, NekDouble> tmp;
      NekDouble vastsqsum;
      Vast = Eorth->PhysEvaluateBasis(xastarrofarr, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

      for( int i = 0; i < N1; i++)
	{
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
  if(verboseflag) 
    cout<<"\nfilter took "<<counter<<"iterations. avgiterGD = "<<avgiterGD;
  uhats = d.back();
}


