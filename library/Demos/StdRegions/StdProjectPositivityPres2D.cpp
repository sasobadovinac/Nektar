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

//declare find_roots, flag is 0 for confederate matrix approach
Array<OneD, Array<OneD, NekDouble> >  find_roots(Array<OneD, NekDouble> &uhats, int d = 0, int flag = 0,  int sig = 0, int surfflag = 0);

Array<OneD,NekDouble>  call_find_roots(Array<OneD,  NekDouble> &uhats1 , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats , NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime);


//declare Opt_needed
int Opt_needed(Array<OneD, NekDouble> uhats, int flag = 0);


// for 2D elements to get uhats at edges
void project_edges( Array<OneD, NekDouble> uhats,Array<OneD, Array<OneD, NekDouble> >&ret );


// uhatpqd stuff: 
void edgederpquhats(Array<OneD, NekDouble> &uhats1, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> &wsp,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2);

void deruhatsedges(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble >&ret,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2 = NullNekDouble1DArray );

vector< vector<NekDouble> > gradient_descent(Array<OneD, NekDouble> uhats, int sig, int surfid = 0 );

string funcdef;
NekDouble  itersGD1, itersGD2, itersGD3;

int iterstaken;
NekDouble  startval, startcoordx, startcoordy, startcoordz;
Array<OneD, Array<OneD, NekDouble> > edgexy(2);
Array<OneD, Array<OneD, NekDouble> > storage2d;

int numedges, numsurfaces;
StdExpansion *E;
StdExpansion *E3seg;

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
  storage2d = E->GetPhysEvaluateStorage(); 
  
  if (E == nullptr)
    {
      return 1;
    }
  dimension = E->GetShapeDimension();
  if(dimension != 2 )
    {
      cout<<"\n dimension should be 2, try using StdProjectPositivityPres1D or StdProjectPositivityPres3DD for other dimensions\n\n";
      exit(0);
    }

  std::vector<int> order;
  std::vector<BasisType> btype(3, eNoBasisType);
  stype = E->DetShapeType();
  LibUtilities::PointsType pointsTypeCheb = LibUtilities::eGaussGaussChebyshev;

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

      demo.testcoord2dtqpts = demo.GetCoords(E);
      demo.testcoord2dtqmidpts = demo.GetQuadratureMidCoords( E, demo.testcoord2dtqpts);
      demo.testcoord2dtlattice = demo.GetLatticeCoords(demo.testcoord2dtqpts, demo.testcoord2dtqmidpts);
      demo.interioreval2dtqmidpts = E->PhysEvaluateBasis( demo.testcoord2dtqmidpts, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
      edgeptsin = Array<OneD, Array<OneD, NekDouble> >(2);

      break;
    
    case LibUtilities::eQuadrilateral:

      numedges = 4;
      numsurfaces = 1;

      demo.testcoord2dqqpts = demo.GetCoords(E);
      demo.testcoord2dqqmidpts = demo.GetQuadratureMidCoords( E, demo.testcoord2dqqpts);
      demo.testcoord2dqlattice = demo.GetLatticeCoords( demo.testcoord2dqqpts, demo.testcoord2dqqmidpts);
      demo.interioreval2dqqmidpts = E->PhysEvaluateBasis( demo.testcoord2dqqmidpts, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
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
	  
	  sol[i] = floor(0.5*(floor(0.5*(x[i]<=-0.4)+0.5*(x[i]>=-0.8)))+0.5*(floor(0.5*(y[i]>=0.4)+0.5*(y[i]<=0.8))))*sin(2*M_PI*(x[i]-0.15))*cos(2*M_PI*(y[i]-0.6));
	  funcdef = "floor(0.5*(floor(0.5*(x[i]<=-0.4)+0.5*(x[i]>=-0.8)))+0.5*(floor(0.5*(y[i]>=0.4)+0.5*(y[i]<=0.8))))*sin(2*M_PI*(x[i]-0.15))*cos(2*M_PI*(y[i]-0.6))";        
	}
      else
	{
	  cout<<"\n Thiw demo only works for 2D\n";
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


  // check for -ve values and apply opt if necessary
  if (vm.count("optm"))
    {
      timer.Start();
      
      vector<LibUtilities::BasisKey> tmpBasisKey = demo.GetBasisKey();       
      LibUtilities::PointsKey pkeycheb(E->GetBasis(0)->GetNumPoints(), pointsTypeCheb);

      LibUtilities::BasisKey bkeycheb(LibUtilities::eChebyshev, 3*(E->GetBasis(0)->GetNumModes()-1),  pkeycheb);
        
      E3seg = new StdSegExp(bkeycheb);
        
      if (E3seg == nullptr)
        {
	  return 1;
        }
        	  
      for(int k = 0; k <2; k++)
	{
	  edgexy[k] = Array<OneD, NekDouble>(E3seg->GetBasis(0)->GetNumPoints());
	  
	  edgeptsin[k] =   Array<OneD, NekDouble>(edgexy[k]);
	}
      Array<OneD, NekDouble> edgexytemp =   E3seg->GetBasis(0)->GetZ();
      int totszedges = edgexytemp.size()*(E->GetNcoeffs());

      // left (x = -1)
      Vxm1 = Array<OneD, NekDouble>(totszedges);
      Vdyxm1 = Array<OneD, NekDouble>(totszedges);
      Vdxxm1 = Array<OneD, NekDouble>(totszedges);

      // left x = -1
      edgexy[1] = edgexytemp;
      edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
      Vxm1 = E->PhysEvaluateBasis(edgexy, storage2d, Vdxxm1, Vdyxm1, NullNekDouble1DArray);

      // bot (y = -1)
      Vym1  = Array<OneD, NekDouble>(totszedges);
      Vdxym1  = Array<OneD, NekDouble>(totszedges);
      Vdyym1  = Array<OneD, NekDouble>(totszedges);

      // bot y = -1
      edgexy[0] = edgexytemp;             
      edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
      Vym1 = E->PhysEvaluateBasis(edgexy,  storage2d, Vdxym1, Vdyym1, NullNekDouble1DArray);
      
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
	  Vx1 = E->PhysEvaluateBasis(edgexy, storage2d, Vdxx1, Vdyx1, NullNekDouble1DArray);
	  
	  //top y = 1
	  edgexy[0] = edgexytemp; 
	  edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), 1.0); 
	  Vy1 = E->PhysEvaluateBasis(edgexy, storage2d,Vdxy1,Vdyy1, NullNekDouble1DArray);
	  
	}
      else if(numedges == 3)
	{
	  // hypt tri (y = -x)
	  Vxyhyp = Array<OneD, NekDouble>(totszedges);
	  Vdxxyhyp = Array<OneD, NekDouble>(totszedges);
	  Vdyxyhyp = Array<OneD, NekDouble>(totszedges);
	  
	  edgexy[0] = edgexytemp;
	  Vmath::Smul(edgexy[0].size(), -1.0, edgexy[0], 1, edgexy[1] , 1);
	  Vxyhyp = E->PhysEvaluateBasis(edgexy, storage2d,  Vdxxyhyp, Vdyxyhyp, NullNekDouble1DArray);
	}

      timer.Stop();
      elapsed  = timer.TimePerTest(1);
      cout<<"\n setup phase took "<<elapsed<<"s\n\n";

      timer.Start();

      int found_negb4 = 0;
      NekDouble  inf = numeric_limits<double>::infinity();
      NekDouble minall = inf;
      // evaluate fn on quad pts
      Array<OneD, NekDouble> tmp1(coeffs.size());
      NekDouble tempv;
      if(numedges == 4)
	{
	  for(int k = 0; k < phys.size(); k++)
	    {
	      if(phys[k] < 0 && abs(phys[k]) >1e-9)
		{
		  if(minall > phys[k])
		    {
		      minall = phys[k];
		      startcoordx = demo.testcoord2dqqpts[0][k];
		      startcoordy = demo.testcoord2dqqpts[1][k];
		      
		    }
		  // cout<<"\n phys neg = "<<phys[k]<<" at "<<demo.testcoord2dqqpts[0][k]<<" "<<demo.testcoord2dqqpts[1][k]<<"\n";
		  found_negb4 = 1;
		}
	    }

	  //cout<<"\n vals at other pts of lattice:\n";
	  for(int k = 0; k < demo.testcoord2dqqmidpts[0].size(); ++k)
	    {
		  
	      Vmath::Vmul(coeffs.size(), &demo.interioreval2dqqmidpts[k], demo.testcoord2dqqmidpts[0].size(), &coeffs[0], 1, &tmp1[0], 1);
	      tempv = Vmath::Vsum(tmp1.size(), tmp1, 1);
	      //  cout<<" "<<tempv<<" ";
	      if(tempv < 0 &&abs(tempv) > 1e-9)
		{
		  //cout<<" val = "<<tempv<<" at "<<demo.testcoord2dqqmidpts[0][k]<<" ,"<<demo.testcoord2dqqmidpts[1][k];
		  if(minall >  tempv)
		    {
		      minall = tempv;
		      startcoordx = demo.testcoord2dqqmidpts[0][k];//<<" "<<demo.testcoord2dqqmidpts[1][k];
		      startcoordy = demo.testcoord2dqqmidpts[1][k];
		      found_negb4 = 1;
		    }
		}
	    }
	}
      else
	{
	  for(int k = 0; k < phys.size(); k++)
	    {
	      if(phys[k] < 0 && abs(phys[k]) >1e-9)
		{
		  //	  cout<<"\n phys neg = "<<phys[k]<<" at "<<demo.testcoord2dtqpts[0][k]<<" "<<demo.testcoord2dtqpts[1][k]<<"\n";
		  if(minall > phys[k])
                    {
                      minall = phys[k];
                      startcoordx = demo.testcoord2dtqpts[0][k];
                      startcoordy = demo.testcoord2dtqpts[1][k];

                    }

		  found_negb4 = 1;
		}
	    }
	  
	  
	  for(int k = 0; k < demo.testcoord2dtqmidpts[0].size(); ++k)
	    {
	      Vmath::Vmul(coeffs.size(), &demo.interioreval2dtqmidpts[k], demo.testcoord2dtqmidpts[0].size(), &coeffs[0], 1, &tmp1[0], 1);
	      tempv = Vmath::Vsum(tmp1.size(), tmp1, 1);
	      if(tempv < 0 &&abs(tempv) > 1e-9)
		{
		  if(minall > tempv)
                    {
                      minall = tempv;
                      startcoordx = demo.testcoord2dtqmidpts[0][k];
                      startcoordy = demo.testcoord2dtqmidpts[1][k];
		      
                    }

		  found_negb4 = 1;
		}
	    }
	      
	}
      timer.Stop();
      elapsed += timer.TimePerTest(1);
      cout<<"\n checking for -vity took "<<timer.TimePerTest(1)<<"s\n\n";

	  

      if(found_negb4)
	{
	  cout<<"\n need optimization\n\n";
	  timer.Start();
	  Do_optimize(coeffs);
	  timer.Stop();
	  elapsed += timer.TimePerTest(1);
	  cout<<"\n optimizertook "<<timer.TimePerTest(1)<<"s\n\n";

	}
      else{
	elapsed += timer.TimePerTest(1);
	  
	cout<<"\n optimizer no need\n\n";
      }
      //Backward transform solution to get projected values
      E->BwdTrans(coeffs, phys);

      int found_neg = 0;
    
      if(numedges == 4)
	{
	  // evaluate fn on quad pts
	  for(int k = 0; k < phys.size(); k++)
	    {
	      if(phys[k] < 0 && abs(phys[k]) >1e-9)
		{
		  // cout<<"\n phys neg = "<<phys[k]<<" at "<<demo.testcoord2dqqpts[0][k]<<" "<<demo.testcoord2dqqpts[1][k]<<"\n";
		  found_neg = 1;
		}

	      for(int k = 0; k < demo.testcoord2dqqmidpts[0].size(); ++k)
		{
		  Vmath::Vmul(coeffs.size(), &demo.interioreval2dqqmidpts[k], demo.testcoord2dqqmidpts[0].size(), &coeffs[0], 1, &tmp1[0], 1);
		  tempv = Vmath::Vsum(tmp1.size(), tmp1, 1);
		  if(tempv < 0 &&abs(tempv) > 1e-9)
		    {
		      // cout<<"\n tempv = "<<tempv<< " at "<<demo.testcoord2dqqmidpts[0][k]<<","<<demo.testcoord2dqqmidpts[1][k];;
		      found_neg = 1;
		    }
		}
	    }
	}
      else
	{
	  // evaluate fn on quad pts
	  for(int k = 0; k < phys.size(); k++)
	    {
	      if(phys[k] < 0 && abs(phys[k]) >1e-9)
		{
		  // cout<<"\n phys neg = "<<phys[k]<<" at "<<demo.testcoord2dtqpts[0][k]<<" "<<demo.testcoord2dtqpts[1][k]<<"\n";
		  found_neg = 1;
		}
	      
	      for(int k = 0; k < demo.testcoord2dtqmidpts[0].size(); ++k)
		{
		  Vmath::Vmul(coeffs.size(), &demo.interioreval2dtqmidpts[k], demo.testcoord2dtqmidpts[0].size(), &coeffs[0], 1, &tmp1[0], 1);
		  tempv = Vmath::Vsum(tmp1.size(), tmp1, 1);
		  if(tempv < 0 &&abs(tempv) > 1e-9)
		    {
		      // cout<<"\n tempv = "<<tempv<<" at "<<demo.testcoord2dtqmidpts[0][k]<<","<<demo.testcoord2dtqmidpts[1][k];
		      found_neg = 1;
		    }
		}
	    }
	}
      if(found_neg)
	{
	  timer.Stop();
	  elapsed += timer.TimePerTest(1);
	  cout<<"\n Verification took "<<timer.TimePerTest(1)<<"s\n\n";
		  
	  cout<<" fail\n";
	  cout<<"\n func="<<funcdef<<"\n";
	  if(numedges == 3)
	    cout<<"\ntri\n";
	  else
	    cout<<"\n quad\n";
	}
      else
	{
	  timer.Stop();
	  elapsed += timer.TimePerTest(1);
	  cout<<"\n Verification took "<<timer.TimePerTest(1)<<"s\n\n";
	      
	  cout<<"pass\n";
	      
	  cout<<"\n func="<<funcdef<<"\n";;
	  if(numedges == 3)
	    cout<<"\ntri\n";
	  else
	    cout<<"\n quad\n";
	      
	  phys = sol;
	}
	  
    }


  //Calculate L_inf & L_2 error
  cout << "L infinity error: \t" << E->Linf(phys, sol) << endl;
  if (stype != ePoint)
    {
      cout << "L 2 error: \t \t \t" << E->L2(phys, sol) << endl;
    }

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


  //    int E3segncoeff = E3seg->GetNcoeffs();

  // for(int k= 0; k < numedges; k++)
  //   {
  // 	ret[k] = Array<OneD, NekDouble>(E3segncoeff);//(3*(E->GetBasis(0)->GetNumModes()-1));        
      
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
  // if(numedges == 6) // tet
  //   {

  //     if(d == 1)
  //       {

  // 	  // edge front left (AD) (x = -1) (y = -1)
  // 	  edgederpquhats(uhats, ret[0], E->GetNcoeffs(),  Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);
                        
  // 	  //edge front hypt (DB) (y = -1) (z = -x) 
  // 	  edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vym1xmz, Vdxym1xmz, Vdyym1xmz, Vdzym1xmz  );
            
  // 	  //edge front bot (AB) (y = -1) (z = -1)  
  // 	  edgederpquhats(uhats, ret[2], E->GetNcoeffs(), Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1 );

  // 	  //edge left hypt (DC) ( x = -1) (z = -y)  
  // 	  edgederpquhats(uhats, ret[3], E->GetNcoeffs(), Vxm1ymz, Vdxxm1ymz, Vdyxm1ymz, Vdzxm1ymz);
     
  // 	  // edge bot diag (BC) (z = -1) (y = -x)
  // 	  edgederpquhats(uhats, ret[4], E->GetNcoeffs(),  Vxm1ymz, Vdxxm1ymz, Vdyxm1ymz, Vdzxm1ymz);

  // 	  //edge CA bot left (x = -1) (z = -1)
  // 	  edgederpquhats(uhats, ret[5], E->GetNcoeffs(),Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);
            
     
  //       }
  //     else if( d == 0)
  //       {
  // 	  // edge front left (AD) (x = -1) (y = -1)
  // 	  deruhatsedges(uhats, ret[0],  Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);

  // 	  //edge front hypt (DB) (y = -1) (z = -x)
  // 	  deruhatsedges(uhats, ret[1],Vym1xmz, Vdxym1xmz, Vdyym1xmz, Vdzym1xmz);
            
  // 	  //edge front bot (AB) (y = -1) (z = -1)
  // 	  deruhatsedges(uhats, ret[2],  Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);

  // 	  //edge left hypt (DC) ( x = -1) (z = -y)
  // 	  deruhatsedges(uhats, ret[3], Vxm1ymz, Vdxxm1ymz, Vdyxm1ymz, Vdzxm1ymz);

  // 	  // edge bot diag (BC) (z = -1) (y = -x)
  // 	  deruhatsedges(uhats, ret[4],  Vxm1ymz, Vdxxm1ymz, Vdyxm1ymz, Vdzxm1ymz);
                   
  // 	  //edge CA bot left (x = -1) (z = -1)
  // 	  deruhatsedges(uhats, ret[5], Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);
          

  //       }
  //   }
  //return ret;   
}

   


void deruhatsedges(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz, Array<OneD, NekDouble> Vdxyz, Array<OneD, NekDouble> Vxdyz, Array<OneD, NekDouble> Vxydz)
{
  boost::ignore_unused(Vxyz, Vxydz);
  int totpts =   E3seg->GetBasis(0)->GetZ().size();
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
        
    }    
  E3seg->FwdTrans(vals,ret);
}


void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> &wsp,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2)
{
  int uhatstot = uhats.size();
  int totpts = edgeptsin[0].size();//E3seg->GetTotPoints();

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

  E3seg->FwdTrans(pqeval, ret);
    
}

Array<OneD,  NekDouble> call_find_roots(Array<OneD,  NekDouble> &uhats1 , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats, NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime)
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
	  rethold  = (demo.find_roots(uhatsedges[ii], E, storage2d, dummy,  0, 0, 0, 0)) ;  
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
	  rethold  = (demo.find_roots(uhatsedges[ii], E, storage2d, dummy,  0, 0, 0, 0)) ;
	  
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
	rethold = demo.find_roots(uhats1, E, storage2d,  avgiterGD, 0, 1, 0 , 1);
      else if(stype == 4)   //else quad
	rethold = demo.find_roots(uhats1, E, storage2d,  avgiterGD, 0, 1, 0 , 0);
	  
      t.Stop();
      roots2dtime +=  t.TimePerTest(1);
      if(rethold[0][0] < inf && rethold[1][0] < inf)
	{
	  retall[0].push_back(rethold[0][0]);
	  retall[1].push_back(rethold[1][0]);
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
	evalroots = E->PhysEvaluateBasis(tmpcoord, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
      else if(stype == 4)
	evalroots = E->PhysEvaluateBasis(tmpcoord, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
	    
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
      Vmath::Vmul(uhats.size(), &storage2d[0][i], nq, &uhats[0], 1, &temp2[0], 1);
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
      Vmath::Vmul(uhats.size(), &demo.interioreval2dtqmidpts[i], nq, &uhats[0], 1, &temp2[0], 1);
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
	  startcoordx = demo.testcoord2dtlattice[0][idxmin];
	  startcoordy = demo.testcoord2dtlattice[1][idxmin];
	}
      return 1;

    }
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
  tols.push_back(1e-12);
    
  Array<OneD, NekDouble>  pqvalcoords(dim+1), xastarr(dim), utemp(N1), wsp1(1);   
    

  Array<OneD, Array<OneD, NekDouble> > surfaceuhats  (numsurfaces), tmpcoord(dim);
  Array<OneD, Array<OneD, NekDouble> > Pf(numedges);
  for(int k= 0; k < numedges; k++)
    {
      Pf[k] = Array<OneD, NekDouble>(E3seg->GetNcoeffs()); 
    }
  for(int k = 0; k < numsurfaces; k++)
    {
      surfaceuhats[k] = Array<OneD, NekDouble>(E->GetNcoeffs());
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
	
	  Vtmp = E->PhysEvaluateBasis(tmpcoord, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        

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
	  Vtmp = E->PhysEvaluateBasis(tmpcoord, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        

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
      Array<OneD, Array<OneD, NekDouble> > xastarrofarr(2);
      xastarrofarr[0] = Array<OneD, NekDouble>(1, xastarr[0]);
      xastarrofarr[1] = Array<OneD, NekDouble>(1, xastarr[1]);

      Array<OneD, NekDouble> tmp;
      NekDouble vastsqsum;
      Vast = E->PhysEvaluateBasis(xastarrofarr, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

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
  cout<<"sphere_rotation took "<<counter<<"iterations. avgiterGD = "<<avgiterGD;
  uhats = d.back();
}


