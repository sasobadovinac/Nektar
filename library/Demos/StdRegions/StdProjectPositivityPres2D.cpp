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
void Do_optimize(StdExpansion *E, Array<OneD, Array<OneD, NekDouble> > &storage, Array<OneD, NekDouble> &uhats);

Array<OneD,NekDouble>  call_find_roots(Array<OneD,  NekDouble> &uhats1 , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats ,Array<OneD, Array<OneD, NekDouble> >&storage,  StdExpansion *E, NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime);

// for 2D elements to get uhats at edges
void project_edges( Array<OneD, NekDouble> uhats,Array<OneD, Array<OneD, NekDouble> >&ret );

// derivative of scaled fn stuff: 
void edgederpquhats(Array<OneD, NekDouble> &uhats1, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> &wsp,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2);

int Opt_needed(Array<OneD, NekDouble> uhats, int flag); 


string funcdef;
NekDouble  itersGD2 ;

int iterstaken;
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
      
      edgeptsin = Array<OneD, Array<OneD, NekDouble> >(2);
      
      break;
    
    case LibUtilities::eQuadrilateral:

      numedges = 4;
      numsurfaces = 1;
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

  //get solution array
  for (int i = 0; i < totPoints; ++i)
    {
      if(dimension ==2) 
	{
	  // f3
	  sol[i] =  (pow(x[i] -0.1, 2) + pow(y[i]-0.2,2))-0.1;///-2*x[i] + pow((x[i] + 0.6),3) + (pow(y[i],2) - 0.2);
	  
	}
      else if (dimension == 3)
	{
	  sol[i] = sqrt(pow(x[i] -0.1, 2) + pow(y[i]-0.2,2))+0.1;//-2*x[i] + pow((x[i] + 0.6),3) + (pow(y[i],2) - 0.2);
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
  LibUtilities::Timer     timer;
  NekDouble elapsed       = 0.0;

  Array<OneD, NekDouble> physold(phys.size());
  Vmath::Vcopy(phys.size(), phys, 1, physold, 1);
  int orthoflag = (E->GetBasisType(0) != LibUtilities::eOrtho_A);
  Array<OneD, BasisType> btypearr(dimension);
  btypearr[0] = LibUtilities::eOrtho_A;

  // check for -ve values and apply opt if necessary
  if (vm.count("optm"))
    {
      timer.Start();
      int found_negb4 = Opt_needed(coeffs, 0);
      
      elapsed += timer.TimePerTest(1);
      cout<<"\nchecking for -vity took "<<timer.TimePerTest(1)<<"s";

	  

      if(found_negb4 == 1)
	{

	  timer.Stop();

	  timer.Start();

	  if(numedges == 4 && demo.Eorthquad != nullptr)
	    {
	      Eorth = demo.Eorthquad;
	      btypearr[1] = LibUtilities::eOrtho_A;
	  
	    }
	  else if (numedges == 3  && demo.Eorthtri != nullptr)
	    {
	      Eorth = demo.Eorthtri;
	      btypearr[1] = LibUtilities::eOrtho_B;     
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
	  cout<<"\nN = "<<Eorth->GetBasis(0)->GetNumModes()<<" p = "<<Eorth->GetBasis(0)->GetNumPoints();
	  cout<<"\nsetup phase took "<<elapsed<<"s";

	  cout<<"\nneed optimization. Negative found at "<<startcoordx<<","<<startcoordy<<" ="<<startval<<"\n";
	  Array<OneD, NekDouble> startarr(dimension);
	  
	  // startarr[0] = startcoordx;
	  // startarr[1] = startcoordy;
	  // demo.SetStartArr(startarr);
	  
	  timer.Start();
	  if(orthoflag)
	    {
	      E->BwdTrans(coeffs, phys);
	      Eorth->FwdTrans(phys, coeffs);
	      //	      demo.OrthoNormalize(E, coeffs, btypearr, 0);
	    }
	  //maybe don't pass storage2d as it is glo var?
	  Do_optimize(Eorth, storage2d,  coeffs);
	  if(orthoflag)
	    {
	      Eorth->BwdTrans(coeffs, phys);
	      E->FwdTrans(phys, coeffs);
	      storage2d = E->GetPhysEvaluateStorage();
	      //	      demo.OrthoNormalize(E, coeffs, btypearr, 1);
	    }
	  timer.Stop();
	  elapsed += timer.TimePerTest(1);
	  cout<<"\noptimizertook "<<timer.TimePerTest(1)<<"s";

	}
      else{
	elapsed += timer.TimePerTest(1);	  
	cout<<"\noptimizer not needed\n";
      }
      //Backward transform solution to get projected values
      E->BwdTrans(coeffs, phys);

      
      int found_neg = Opt_needed(coeffs, 1);

      
      // check val now at the min ppreviously found:
      Array<OneD, Array< OneD, NekDouble> > checkcoord(2);
      Array<OneD, NekDouble> saveeval, holdv(coeffs.size());
      checkcoord[0] = Array<OneD, NekDouble>(1, startcoordx);
      checkcoord[1] = Array<OneD, NekDouble>(1, startcoordy);
      
      saveeval = E->PhysEvaluateBasis(checkcoord, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

      Vmath::Vmul(coeffs.size(), saveeval, 1, coeffs, 1, holdv, 1);
      NekDouble valnow =  Vmath::Vsum(coeffs.size(), holdv, 1);
      cout<<"\n val now at "<<startcoordx<<" "<<startcoordy<<" "<<valnow;
      if(valnow < 0 && abs(valnow)>1e-9)
	{
	  found_neg = 1;
	}
      
      if(found_neg)
	{
	  timer.Stop();
	  elapsed += timer.TimePerTest(1);
	  cout<<"\nVerification took "<<timer.TimePerTest(1)<<"s";
		  
	  cout<<" \nfail";
	  if(numedges == 3)
	    cout<<"\ntri\n";
	  else
	    cout<<"\nquad";

	  // in order for test to fail, sabotahe L2 and Linf errors here:
	  cout<<"\nConstrained:";
	  NekDouble dummyv = 1;
	  cout << "\nL infinity error: \t" << dummyv << endl;
	  if (stype != ePoint)
	    {
	      cout << "L 2 error: \t \t \t" << 1 << dummyv << endl;
	    }
	  
	}
      else
	{
	  timer.Stop();
	  elapsed += timer.TimePerTest(1);
	  cout<<"\nVerification took "<<timer.TimePerTest(1)<<"s";
	  cout<<"\npass";
	  if(numedges == 3)
	    cout<<"\ntri";
	  else
	    cout<<"\nquad";
	      
	  phys = sol;
	  cout<<"\nConstrained:";
	  NekDouble dummyv = 0;
	  cout << "\nL infinity error: \t" << dummyv << endl;
	  if (stype != ePoint)
	    {
	      cout << "L 2 error: \t \t \t" << dummyv << endl;
	    }
	  
	}
	  
    }
   cout<<"\n*******************\n";
  return 0;
    
}

Array<OneD, Array<OneD, NekDouble> > vectoarr(vector<vector <NekDouble> > vec)
{
  Array<OneD, Array<OneD, NekDouble> > ret(vec.size());
  for (int k = 0; k < vec.size(); k++)
    {
      ret[k ] = Array<OneD, NekDouble>(vec[k].size(), vec[k].data());
    }
  return ret;
}


// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble>uhats1,    Array<OneD, Array<OneD, NekDouble> >&ret)
{
    
  Array<OneD, NekDouble> wsp(uhats1.size());


  //    int demo.E3segncoeff = demo.E3seg->GetNcoeffs();

  // for(int k= 0; k < numedges; k++)
  //   {
  // 	ret[k] = Array<OneD, NekDouble>(demo.E3segncoeff);//(3*(E->GetBasis(0)->GetNumModes()-1));        
      
  //   }
  // cout<<"\nhere2\n\n";

  if(numedges == 4) // quad
    {
      // bot edge
      edgederpquhats(uhats1, ret[0], wsp, Vym1, Vdxym1, Vdyym1, NullNekDouble1DArray);
      // right edge
      edgederpquhats(uhats1, ret[1],wsp, Vx1, Vdxx1, Vdyx1, NullNekDouble1DArray);
      // top edge
      edgederpquhats(uhats1, ret[2], wsp, Vy1, Vdxy1, Vdyy1, NullNekDouble1DArray);
      // left edge
      edgederpquhats(uhats1, ret[3], wsp, Vxm1, Vdxxm1, Vdyxm1, NullNekDouble1DArray);

    }
  if(numedges == 3) // tri
    {
      // bot edge  
      edgederpquhats(uhats1, ret[0], wsp,  Vym1, Vdxym1, Vdyym1,NullNekDouble1DArray);

      // left edge
      edgederpquhats(uhats1, ret[1], wsp, Vxm1, Vdxxm1, Vdyxm1, NullNekDouble1DArray);

      // hypto
      edgederpquhats(uhats1, ret[2], wsp, Vxyhyp, Vdxxyhyp, Vdyxyhyp, NullNekDouble1DArray);
    }

}

int Opt_needed(Array<OneD, NekDouble> coeffs, int flag) 
{

  int nq = E->GetTotPoints();
  Array<OneD, NekDouble> temp2(coeffs.size()), holdmidpteval;
  Array<OneD, Array<OneD, NekDouble> > holdlattice;
  NekDouble minv = numeric_limits<double>::infinity(), holdv;
  int idxmin = -1;
  
  for (int i = 0; i<nq; i++)
    {
      Vmath::Vmul(coeffs.size(), &storage2d[0][i], nq, &coeffs[0], 1, &temp2[0], 1);  
      holdv = Vmath::Vsum(coeffs.size(), temp2, 1);    
      if(holdv < minv)
	{
	  minv = holdv;
	  idxmin = i;
	}
    }
  
  int ct = nq;
 
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
  int midptsize = holdlattice[0].size() - nq;
  for(int i = 0,  k = ct; k < holdlattice[0].size(); ++k, i++)
    {
      Vmath::Vmul(coeffs.size(), &holdmidpteval[i], midptsize, &coeffs[0], 1, &temp2[0], 1);   
      holdv = Vmath::Vsum(coeffs.size(), temp2, 1);
      if(holdv < minv)
	{
	  minv = holdv;
	  idxmin = k;
	}
    }
  cout<<"\n minv = "<<minv<<" at "<< holdlattice[0][idxmin]<< " "<<holdlattice[1][idxmin];

  if(minv < 0 && abs(minv) > 1e-10 && idxmin != -1)
    {
      if(flag == 0)
	{
	  startval = minv;
	  startcoordx = holdlattice[0][idxmin];
	  startcoordy = holdlattice[1][idxmin];
	}
      cout<<"\n minv = "<<minv;
      return 1;
    }
  return 0;
}

void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> &wsp,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2)
{
  int uhatstot = uhats.size();
  int totpts = edgeptsin[0].size();//demo.E3seg->GetTotPoints();

  Array<OneD, NekDouble> pqeval(totpts);

  NekDouble v1, v2;
  for(int i = 0; i<totpts; i++)
    {

      Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd0[0]+i, totpts, &wsp[0], 1);
      v1  = Vmath::Vsum(uhatstot, wsp, 1); 
      Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &wsp[0], 1);

      v2  = Vmath::Vsum(uhatstot, wsp, 1);  

      v1 = v2*v1;

      // At this point,                 

      // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

      // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
      Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &wsp[0], 1);

 
        
      v2  = Vmath::Vsum(uhatstot, wsp, 1);  
 
      Vmath::Vmul(uhats.size(), &Vxyd0[i], totpts, &uhats[0], 1, &wsp[0], 1);

      v1= v2*Vmath::Vsum(uhats.size(), wsp, 1)- v1;
  
      pqeval[i] = v1;
 
      if(dimension > 1)
	{
	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd1[0]+i, totpts, &wsp[0], 1);
	  v1  = Vmath::Vsum(uhatstot, wsp, 1); 
	  Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &wsp[0], 1);

	  v2  = Vmath::Vsum(uhatstot, wsp, 1);  

	  v1 = v2*v1;

	  // At this point,                 

	  // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

	  // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &wsp[0], 1);

 
        
	  v2  = Vmath::Vsum(uhatstot, wsp, 1);  
 
	  Vmath::Vmul(uhats.size(), &Vxyd1[i], totpts, &uhats[0], 1, &wsp[0], 1);

	  v1=  v2*Vmath::Vsum(uhats.size(), wsp, 1) - v1;
 
	  pqeval[i] += v1;
 
	}
      if(dimension == 3)
	{
	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd2[0]+i, totpts, &wsp[0], 1);
	  v1  = Vmath::Vsum(uhatstot, wsp, 1); 
	  Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &wsp[0], 1);

	  v2  = Vmath::Vsum(uhatstot, wsp, 1);  

	  v1 = v2*v1;

	  // At this point,                 

	  // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

	  // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &wsp[0], 1);

 
        
	  v2  = Vmath::Vsum(uhatstot, wsp, 1);  
 
	  Vmath::Vmul(uhats.size(), &Vxyd2[i], totpts, &uhats[0], 1, &wsp[0], 1);

	  v1= v2*Vmath::Vsum(uhats.size(), wsp, 1) - v1;
 
	  pqeval[i] += v1;
 
	}
    }

  demo.E3seg->FwdTrans(pqeval, ret);
    
}

//try replacing storage with the global var -> storage2d
Array<OneD,NekDouble>  call_find_roots(Array<OneD,  NekDouble> &uhats1 , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats ,Array<OneD, Array<OneD, NekDouble> >&storage,  StdExpansion *E, NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime)
{
  
  boost::ignore_unused(surfaceuhats, roots2dtime, roots3dtime, roots1dtime);
  int dimension = E->GetShapeDimension(); 
  NekDouble dummy;
  Array<OneD, NekDouble> temp(uhats1.size());

  NekDouble  inf = numeric_limits<double>::infinity();
  //    NekDouble avgiterGDhold;
  Timer t;
  Array<OneD, Array<OneD, NekDouble > >  rethold(dimension);
  vector<vector<NekDouble> > retall(dimension);
  Array<OneD, Array<OneD, NekDouble> > tmpcoord(dimension);

  Array<OneD, NekDouble> ret(dimension);

  if(numedges == 4) // quad
    {
          
      for(int ii = 0; ii < numedges; ii++)
	{

	  t.Start();
	  rethold  = (demo.find_roots(uhatsedges[ii], nullptr, NullNekDoubleArrayOfArray, dummy,  0, 0, 0, 0)) ;  
	  t.Stop();

	  roots1dtime +=  t.TimePerTest(1);
	  for(int p = 0; p < rethold[0].size(); p++)
	    {
	      switch(ii)
		{
		case 0:
		  //if(ii == 0) // edge bot (y = -1) 
		  retall[1].push_back(   (-1.0));
		  retall[0].push_back( (rethold[0][p]));  
		  break;
		case 1:
		  //if(ii == 1) // edge right (x = 1) 
		  retall[0].push_back(   (1.0));
		  retall[1].push_back( (rethold[0][p]));  
		  break;
		case 2: 
		  //if(ii == 2) // edge top (y = 1) 
		  retall[1].push_back(   (1.0));
		  retall[0].push_back( (rethold[0][p]));  
		  break;
		case 3:
		  //if(ii == 3) // edge left (x = -1) 
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

      // cout<<"\n edge roots:\n";
      // for(int k = 0; k < retall[0].size(); k++)
      // 	cout<<" "<<retall[0][k]<<" "<<retall[1][k]<<" \n";
    }
  else if(numedges == 3) //tri
    {
      
      for(int ii = 0; ii < numedges; ii++)
	{
	  t.Start();
	  rethold  = (demo.find_roots(uhatsedges[ii], nullptr, NullNekDoubleArrayOfArray, dummy,  0, 0, 0, 0)) ;
	  
	  t.Stop();
	  roots1dtime +=  t.TimePerTest(1);
	  for(int p = 0; p < rethold[0].size(); p++)
	    {
	      switch(ii)
		{
		case 0:
		  //if(ii == 0) // edge bot (y = -1) 
		  retall[1].push_back(   (-1.0));
		  retall[0].push_back( (rethold[0][p]));  
		  break;
		case 1:
		  //if(ii == 1) // edge left (x = -1) 
		  retall[0].push_back(   (-1.0));
		  retall[1].push_back( (rethold[0][p]));  
		  
		  break;
		case 2: 
		  //if(ii == 2) // edge hypt (y = -x) 
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
	rethold = demo.find_roots(uhats1, E, storage,  avgiterGD, 0, 1, 0 , 1);
      else if(stype == 4)   //else quad
	rethold = demo.find_roots(uhats1, E, storage,  avgiterGD, 0, 1, 0 , 0);
	  
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
	evalroots = E->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
      else if(stype == 4)
	evalroots = E->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
	    
      Array<OneD, NekDouble> minvandx(dimension+1);
      Array<OneD, NekDouble> tmparr(retall[0].size());
      demo.pq(uhats1, tmpcoord, evalroots, minvandx, tmparr);
      // cout<<"\n vals at roots:\n";
      // for(int k = 0; k < tmparr.size(); k++)
      //   cout<<" "<<tmpcoord[0][k]<<","<<tmpcoord[1][k]<<" = "<<tmparr[k]<<" \n";
      // cout<<"\n";
	  
      tempmin = minvandx[0];
      ret[0] = minvandx[1];
      ret[1] = minvandx[2];
	  
      minv = tempmin;
    }
  return ret;

}

void Do_optimize(StdExpansion *E, Array<OneD, Array<OneD, NekDouble> > &storage, Array<OneD, NekDouble> &uhats)
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
      surfaceuhats[k] = Array<OneD, NekDouble>(E->GetNcoeffs());
    }

  NekDouble pqval, timeprojectsurf = 0.0, timeprojectedges = 0.0;
  NekDouble avgiterGD = 0.0, avgiterhold = 0.0, roots1dtimehold = 0.0, roots2dtimehold = 0.0, roots3dtimehold = 0.0;
  Array<OneD, NekDouble> Vtmp(E->GetNcoeffs());

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

	  optima = call_find_roots(utemp, avgiterhold, Pf, surfaceuhats, storage, E, minv, roots1dtime, roots2dtime, roots3dtime);
	  roots1dtimehold += roots1dtime;
	  roots2dtimehold += roots2dtime;
	  roots3dtimehold += roots3dtime;
	  avgiterGD += avgiterhold;
	  // Vtmp is evaluation of basis fn at optima
	  cout<<"\n opt:\n";
	  for(int k = 0; k < dim; k++)
	    {
	      tmpcoord[k] = Array<OneD, NekDouble>(1, optima[k]);
	      cout<<" "<<optima[k];
	    }
	
	  Vtmp = E->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        
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
	  Vtmp = E->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        

	}
      demo.pq(utemp, tmpcoord, Vtmp, pqvalcoords, wsp1);

      if (pqvalcoords[0] < pqval)
	{
	  xastarr = optima;
	  pqval = pqvalcoords[0];
	}
      
      cout<<"\n pqval = "<<pqval<<" at "<<xastarr[0]<<" "<<xastarr[1]<<"\n";        

        
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
      Vast = E->PhysEvaluateBasis(xastarrofarr, storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

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
  cout<<"\nfilter took "<<counter<<"iterations. avgiterGD = "<<avgiterGD;
  uhats = d.back();
}


