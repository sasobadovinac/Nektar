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

// for 2D elements to get uhats at edges
void project_edges( Array<OneD, NekDouble> uhats,Array<OneD, Array<OneD, NekDouble> >&ret );


//int Opt_needed(Array<OneD, NekDouble> uhats, Array<OneD, NekDouble> &evallattice, int flag ); 
int Opt_needed(Array<OneD, NekDouble> &uhatslocal, Array <OneD, NekDouble>& t,int flag );
Array<OneD, NekDouble> eval_sol(Array<OneD, NekDouble> x, Array<OneD, NekDouble> y, Array<OneD, NekDouble> z);
void call_setup_quad(StdExpansion *E);
void call_setup_tri(StdExpansion *E);

string funcdef;
NekDouble  itersGD2 ;

int iterstaken, verboseflag;
NekDouble  startcoordx, startcoordy, startcoordz, startval;
Array<OneD, Array<OneD, NekDouble> > edgexy(2);
Array<OneD, Array<OneD, NekDouble> > storage2dt, storage2dq, storage2d;

int numedges, numsurfaces;
StdExpansion *E;
StdExpansion *Eorth;

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
  LibUtilities::PointsKey pkeycheb( 3*(E->GetBasis(0)->GetNumPoints()),LibUtilities::eGaussLobattoChebyshev);
  LibUtilities::BasisKey bkeyorth(LibUtilities::eOrtho_A, 3*(E->GetBasis(0)->GetNumModes()),  pkeycheb);
  demo.E3seg = new StdSegExp(bkeyorth);

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
      call_setup_tri(E);
	  
      numedges = 3;
      numsurfaces = 1;
      edgeptsin = Array<OneD, Array<OneD, NekDouble> >(2);
      
      break;
    
    case LibUtilities::eQuadrilateral:
      call_setup_quad(E);
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
  cout<<"\n unconstrained sol on quadrature points (Linf, L2)::"<< E->Linf(phys, sol)<<", "<< E->L2(phys, sol)   ;
  
  Array<OneD, NekDouble> retvalslattice, solhold;
  // check for -ve values and apply opt if necessary
  if (vm.count("optm"))
    {
      if(E->DetShapeType() == LibUtilities::eTriangle)
	{
	  sollattice = eval_sol(demo.coordlatticetri[0], demo.coordlatticetri[1], NullNekDouble1DArray );
	}
      if(E->DetShapeType() == LibUtilities::eQuadrilateral)
	{
	  sollattice = eval_sol(      demo.coordlatticequad[0],       demo.coordlatticequad[1], NullNekDouble1DArray );

	}
      
      Eorth = E;
      Array<OneD, NekDouble> retvalslatticesave;
      storage2d = Eorth->GetPhysEvaluateStorage();

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
      cout << "\nUnconstrained L infinity error: \t" << Linferrlattice << endl;
      cout << "Unconstrained L2 error \t: \t" << L2errlattice << endl;
      
      if(found_negb4 == 1)
	{
	  int orthoflag = (E->GetBasisType(0) != LibUtilities::eOrtho_A);
	  
	  if(orthoflag)
	    {
	      timer.Start();
	      
	      E->BwdTrans(coeffs, phys);
	      Eorth->FwdTrans(phys, coeffs);
	      timer.Stop();
	      cout<<"\northonorm 1 took "<<timer.TimePerTest(1)<<"  s";
	      
	    }
	  timer.Start();
	  Do_optimize(coeffs);
	  timer.Stop();
	  //if(verboseflag)
	  //{
	  cout<<"\noptimizertook "<<timer.TimePerTest(1)<<"  s";
	  // }
	  if(orthoflag)
	    {
	      timer.Start();
	 
	      Eorth->BwdTrans(coeffs, phys);
	      E->FwdTrans(phys, coeffs);
	      storage2d = E->GetPhysEvaluateStorage();
	      timer.Stop();
	      cout<<"\northonorm 2 took  "<<timer.TimePerTest(1)<<"  s";
	    }
	  else
	    {
	      E->BwdTrans(coeffs, phys);
	    }

	  // cout<<"\ncsv for exact sol on lattice:\n";
	  // for(int p = 0; p < sollattice.size(); p++)
	  //   {
	  //     cout<<" \n"<< demo.coordlatticequad[0][p]<<"," << demo.coordlatticequad[1][p]<<",0,"<<sollattice[p];
	  //   }
	  
 	  cout<<"\n Constrained sol on quad pts (Linf, L2)::"<< E->Linf(phys, sol)<<" ,"<< E->L2(phys, sol) ;
	  //Backward transform solution to get projected values
	  int found_neg  = Opt_needed(coeffs, retvalslattice, 1);

	  //L2 error between retvalslattice and eval_sol on lattice:
	  Vmath::Vcopy(sollattice.size(), sollattice, 1, solhold, 1);
	  Vmath::Vsub(latticesz, retvalslattice, 1, sollattice, 1,
		      sollattice, 1);
	  Vmath::Vabs(latticesz, sollattice, 1, sollattice, 1);
	  Linferrlattice = Vmath::Vmax(latticesz, sollattice, 1);
	  Vmath::Vmul(latticesz, sollattice, 1, sollattice, 1, sollattice, 1);
	  L2errlattice = sqrt(Vmath::Vsum(latticesz, sollattice, 1) );

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
		  cout<<"\nVerification took "<<timer.TimePerTest(1)<<"  s";
		}
	      cout<<" \nfail tests";
	      //	      fail unit test:
	      
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
		  cout<<"\nVerification took "<<timer.TimePerTest(1)<<"  s";
		}
	      cout<<"\npass tests";
		  
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



void project_edges( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret)
{

  if(Eorth->DetShapeType()  == LibUtilities::eQuadrilateral) // quad
    {
      // bot edge
      edgederpquhats(uhats, ret[0], Eorth->GetNcoeffs(),  Vym1q, Vdxym1q, Vdyym1q, NullNekDouble1DArray);
      // right edge
      edgederpquhats(uhats, ret[1], Eorth->GetNcoeffs(), Vx1q, Vdxx1q, Vdyx1q, NullNekDouble1DArray);
      // top edge
      edgederpquhats(uhats, ret[2], Eorth->GetNcoeffs(), Vy1q, Vdxy1q, Vdyy1q, NullNekDouble1DArray);
      // left edge
      edgederpquhats(uhats, ret[3], Eorth->GetNcoeffs(),  Vxm1q, Vdxxm1q, Vdyxm1q, NullNekDouble1DArray);

      return;
    }
  if(Eorth->DetShapeType() == LibUtilities::eTriangle) // tri
    {
	  
      // bot edge  
      edgederpquhats(uhats, ret[0], Eorth->GetNcoeffs(),  Vym1t, Vdxym1t, Vdyym1t,NullNekDouble1DArray);
      // left edge
      edgederpquhats(uhats, ret[1], Eorth->GetNcoeffs(), Vxm1t, Vdxxm1t, Vdyxm1t, NullNekDouble1DArray);
	  
      // hypto
      edgederpquhats(uhats, ret[2], Eorth->GetNcoeffs(),  Vxyhypt, Vdxxyhypt, Vdyxyhypt, NullNekDouble1DArray);
	        
      return;
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
	  sol[i] = floor(0.5*(floor(0.5*(x[i]>= -0.8) + 0.5*(x[i]<=0.2))) +  0.5*(floor(0.5*( y[i]>=-0.2) + 0.5*(y[i]<=0.8))))*sin(M_PI*(0.2-x[i]))*sin(M_PI*(y[i]+0.2));
	}	  
      else
	{
	  cout<<"\n This demo only tests 2-D and 3-D root-fnding\n";
	  exit(0);
	}
    }
  return sol;
}
 
// FLAG == 1 - > verification,
// flag == 0 - > identification
int Opt_needed(Array<OneD, NekDouble> &coeffs, Array<OneD, NekDouble> &valslattice, int flag)
{
  int nq = E->GetTotPoints();
  Array<OneD, NekDouble> temp2(coeffs.size()), holdmidpteval;
  Array<OneD, Array<OneD, NekDouble> > holdlattice, storage2dhold;
  NekDouble minv = numeric_limits<double>::infinity();
  int idxmin = -1;

  if(E->DetShapeType() == LibUtilities::eQuadrilateral)
    {
      holdmidpteval = demo.midptevalquad;
      holdlattice = demo.coordlatticequad;
      storage2dhold = demo.storage2dq;
    }
  else if(E->DetShapeType() == LibUtilities::eTriangle)
    {
      holdmidpteval = demo.midptevaltri;
      holdlattice = demo.coordlatticetri;
      storage2dhold = demo.storage2dt; 
    }

  Array<OneD, NekDouble> holdv(holdlattice[0].size());

  for (int i = 0; i<nq; i++)
    {
      Vmath::Vmul(coeffs.size(), &storage2dhold[0][i], nq, &coeffs[0], 1, &temp2[0], 1);
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


void call_setup_tri(StdExpansion * exp)
{
  int orthoflag = (exp->GetBasisType(0) != LibUtilities::eOrtho_A);
  if(!orthoflag)
    {
      int nmodes0 = exp->GetBasis(0)->GetNumModes();
	  
      int nmodes1 = exp->GetBasis(1)->GetNumModes();
      int npts0 = exp->GetBasis(0)->GetNumPoints();
      int npts1 = exp->GetBasis(1)->GetNumPoints();
      PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
      PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
      BasisKey  b0(LibUtilities::eOrtho_A,  nmodes0, p0);
      BasisKey  b1(LibUtilities::eOrtho_B,  nmodes1, p1);
      E = new StdTriExp(b0, b1);
    }
  else
    {
      E = exp;
    }
  demo.storage2dt =  E->GetPhysEvaluateStorage(); 
  demo.coordtri = demo.GetCoords(E);
  demo.coordmidtri = demo.GetQuadratureMidCoords(demo.coordtri);
  demo.coordlatticetri = demo.GetLatticeCoords(demo.coordtri, demo.coordmidtri);

  demo.midptevaltri = E->PhysEvaluateBasis(demo.coordmidtri, demo.storage2dt, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

  Array<OneD, Array<OneD, NekDouble> > edgeptsin(2), edgexy(2);
  for(int k = 0; k <2; k++)
    {
      edgexy[k] = Array<OneD, NekDouble>(demo.E3seg->GetBasis(0)->GetNumPoints());
      edgeptsin[k] =   Array<OneD, NekDouble>(edgexy[k]);
    }
  Array<OneD, NekDouble> edgexytemp = demo.E3seg->GetBasis(0)->GetZ();
  int totszedges = edgexytemp.size()*(E->GetNcoeffs());
      
  // left (x = -1)
  Vxm1t = Array<OneD, NekDouble>(totszedges);
  Vdyxm1t = Array<OneD, NekDouble>(totszedges);
  Vdxxm1t = Array<OneD, NekDouble>(totszedges);
      
  // left x = -1
  edgexy[1] = edgexytemp;
  edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
  Vxm1t = E->PhysEvaluateBasis(edgexy, demo.storage2dt, Vdxxm1t, Vdyxm1t, NullNekDouble1DArray);
      
  // bot (y = -1)
  Vym1t  = Array<OneD, NekDouble>(totszedges);
  Vdxym1t  = Array<OneD, NekDouble>(totszedges);
  Vdyym1t  = Array<OneD, NekDouble>(totszedges);
      
  // bot y = -1
  edgexy[0] = edgexytemp;             
  edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
      
  Vym1t = E->PhysEvaluateBasis(edgexy,  demo.storage2dt, Vdxym1t, Vdyym1t, NullNekDouble1DArray);
      
  // hypt tri (y = -x)
  Vxyhypt = Array<OneD, NekDouble>(totszedges);
  Vdxxyhypt = Array<OneD, NekDouble>(totszedges);
  Vdyxyhypt = Array<OneD, NekDouble>(totszedges);
      
  edgexy[0] = edgexytemp;
  Vmath::Smul(edgexy[0].size(), -1.0, edgexy[0], 1, edgexy[1] , 1);
  Vxyhypt = E->PhysEvaluateBasis(edgexy, demo.storage2dt,  Vdxxyhypt, Vdyxyhypt, NullNekDouble1DArray);
      
}

    
void call_setup_quad(StdExpansion* exp)
{
  int orthoflag = (exp->GetBasisType(0) != LibUtilities::eOrtho_A);
  if(!orthoflag)
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
      E = new StdQuadExp(b0, b1);
    }
  else
    {
      E = exp;
    }
  demo.storage2dq =  E->GetPhysEvaluateStorage();
  demo.coordquad = demo.GetCoords(E);
  demo.coordmidquad = demo.GetQuadratureMidCoords(demo.coordquad);
  demo.coordlatticequad = demo.GetLatticeCoords(demo.coordquad, demo.coordmidquad);

  demo.midptevalquad = E->PhysEvaluateBasis(demo.coordmidquad, demo.storage2dq, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
      
  Array<OneD, Array<OneD, NekDouble> > edgeptsin(2), edgexy(2);
  for(int k = 0; k <2; k++)
    {
      edgexy[k] = Array<OneD, NekDouble>(demo.E3seg->GetBasis(0)->GetNumPoints());
      edgeptsin[k] =   Array<OneD, NekDouble>(edgexy[k]);
    }
  Array<OneD, NekDouble> edgexytemp =   demo.E3seg->GetBasis(0)->GetZ();
  int totszedges = edgexytemp.size()*(E->GetNcoeffs());
      
  // left (x = -1)
  Vxm1q = Array<OneD, NekDouble>(totszedges);
  Vdyxm1q = Array<OneD, NekDouble>(totszedges);
  Vdxxm1q = Array<OneD, NekDouble>(totszedges);
      
  // left x = -1
  edgexy[1] = edgexytemp;
  edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
  Vxm1q = E->PhysEvaluateBasis(edgexy, demo.storage2dq, Vdxxm1q, Vdyxm1q, NullNekDouble1DArray);
      
  // bot (y = -1)
  Vym1q  = Array<OneD, NekDouble>(totszedges);
  Vdxym1q  = Array<OneD, NekDouble>(totszedges);
  Vdyym1q  = Array<OneD, NekDouble>(totszedges);
      
  // bot y = -1
  edgexy[0] = edgexytemp;             
  edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
      
  Vym1q = E->PhysEvaluateBasis(edgexy,  demo.storage2dq, Vdxym1q, Vdyym1q, NullNekDouble1DArray);

      	  
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
  Vx1q = E->PhysEvaluateBasis(edgexy, demo.storage2dq, Vdxx1q, Vdyx1q, NullNekDouble1DArray);
      
  //top y = 1
  edgexy[0] = edgexytemp; 
  edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), 1.0); 
  Vy1q = E->PhysEvaluateBasis(edgexy, demo.storage2dq, Vdxy1q, Vdyy1q, NullNekDouble1DArray);
      
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
	  rethold  = demo.call_companion_rf(uhatsedges[i]);  
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

	  rethold  = demo.call_companion_rf(uhatsedges[i]);
	  
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
	  demo.steepestgradient_descent2D(uhats, Eorth, rethold, avgiterGD, 1e-7);  
	}
      else if(stype == 4)   //else quad
	{
	  demo.steepestgradient_descent2D(uhats, Eorth, rethold, avgiterGD, 1e-7);  
	}
	  
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

  NekDouble pqval;//, timeprojectsurf = 0.0, timeprojectedges = 0.0;
  NekDouble avgiterGD = 0.0, avgiterhold = 0.0;// roots1dtimehold = 0.0, roots2dtimehold = 0.0, roots3dtimehold = 0.0;
  Array<OneD, NekDouble> Vtmp(Eorth->GetNcoeffs());
  NekDouble roots1dtime = 0.0, roots2dtime = 0.0, roots3dtime = 0.0 ;

  Timer t;
  while (counter <= niter)
    {
      pqval = inf;
      utemp = d.back();
      //t.Start();
      //pq = @(xx) cfun(xx) .* (constraints(1,jj).*Vc(xx)*uhats);  
      project_edges(utemp, Pf);
      //t.Stop();
      //	  timeprojectedges += t.TimePerTest(1);

      optima = call_find_roots(utemp, avgiterhold, Pf, surfaceuhats, minv, roots1dtime, roots2dtime, roots3dtime);
      avgiterGD += avgiterhold;
      // Vtmp is evaluation of basis fn at optima
      for(int k = 0; k < dim; k++)
	{
	  tmpcoord[k] = Array<OneD, NekDouble>(1, optima[k]);
	}
      
      Vtmp = Eorth->PhysEvaluateBasis(tmpcoord, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
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
  // roots1dtimehold = roots1dtimehold/(counter-1);
  // roots2dtimehold = roots2dtimehold/(counter-1);
  // roots3dtimehold = roots3dtimehold/(counter-1);
  // timeprojectedges = timeprojectedges/(counter);
  // timeprojectsurf = timeprojectsurf/(counter);    
  // itersGD2 = (avgiterGD)/(counter);
  // iterstaken = counter;
  //  if(verboseflag) 
  cout<<"\nfilter took "<<counter<<"iterations.\nTotal iterations taken by GD = "<<avgiterGD<<"\nroots 1d time = "<<roots1dtime<<"\nroots2d time = "<< roots2dtime<<"\n";
  uhats = d.back();
}


