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
#include "StdDemoSupportnew.hpp"
#include <LibUtilities/BasicUtils/Timer.h>

namespace po = boost::program_options;

NekDouble Shape_sol(NekDouble x, NekDouble y, NekDouble z, vector<int> order,
                    vector<BasisType> btype, ShapeType stype, bool diff);

//Modification to deal with exact solution for diff. Return 1 if integer < 0.
static double pow_loc(const double val, const int i)
{
  return (i < 0) ? 1.0 : pow(val, i);
}

Array<OneD, Array<OneD, NekDouble> > vectoarr(vector<vector< NekDouble> > vec);


//declare Do_optimize
void Do_optimize(Array<OneD, NekDouble> &uhats);

//declare find_roots, flag is 0 for confederate matrix approach
vector<vector<  NekDouble> > find_roots(Array<OneD, NekDouble> &uhats, int d = 0, int flag = 0,  int sig = 0, int surfflag = 0);

// declare caller routine to find_roots
// flag = 0 -> opt_needed calls
// flag = 1 -> sphere_rot calls
Array<OneD, Array<OneD,  NekDouble> >  call_find_roots(Array<OneD, NekDouble> &uhatsall, NekDouble &avgiterGD, int d = 0, Array< OneD, Array<OneD, NekDouble> >&uhatsedges =NullNekDoubleArrayofArray , int flag = 0 );


//declare Opt_needed
int Opt_needed(Array<OneD, NekDouble> uhats, int flag = 0);


// for 2D elements to get uhats at edges
// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble> uhats,Array<OneD, Array<OneD, NekDouble> >&ret , int d = 0);


// uhatpqd stuff: 
// called by project_edges if d = 1 is passed in project_edges
void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret, int modes,  Array<OneD, NekDouble> V3, Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1 = NullNekDouble1DArray , Array<OneD, NekDouble> Vxyd2 = NullNekDouble1DArray);

void deruhatsedges(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble >&ret,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2 = NullNekDouble1DArray );

vector< vector<NekDouble> > gradient_descent(Array<OneD, NekDouble> uhats, int sig, int surfid = 0 );

string funcdef;
NekDouble  startval, startcoordx, startcoordy, startcoordz;
Array<OneD, Array<OneD, NekDouble> > edgexy(2);

// Colleague matrix
int numedges, numsurfaces;
//Array<OneD, Array<OneD, NekDouble> > C;
StdExpansion *E;
StdExpansion *E3seg;
Array<OneD, NekDouble> qZin; //inner point grid for 2D rootfinding
Array<OneD, NekDouble> qZinmid; //inner mid point grid for 2D rootfinding
Array<OneD, NekDouble> qWin; //inner point grid for 2D rootfinding
Array<OneD, NekDouble> qx;
Array<OneD, NekDouble> qw;

Array<OneD, NekDouble> V;
Array<OneD, NekDouble> Vd;
Array<OneD, NekDouble> V3;
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
Array<OneD, NekDouble > qWinarr;

// for triangle edges root-finding:
Array<OneD, NekDouble> Vxyhyp;
Array<OneD, NekDouble> Vdxxyhyp;
Array<OneD, NekDouble> Vdyxyhyp;

  

int dimension ;
Array<OneD, NekDouble > interioreval2d;
// Array<OneD, NekDouble > interioreval2d0;
// Array<OneD, NekDouble > interioreval2d1;

Array<OneD, Array<OneD, NekDouble> > testcoord2d; 
            
DemoSupportNew demo;

int main(int argc, char *argv[])
{
  demo.GetOptions().add_options()("optm,z", "positivity preserving optimizer");
  //    demo.GetOptions().add_options()("optmtol,t", "tol of opt");
  demo.ParseArguments(argc, argv);

  po::variables_map vm = demo.GetVariableMap();

  //only for 1D
       
  E = demo.CreateStdExpansion();
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
  LibUtilities::ShapeType stype = E->DetShapeType();
  LibUtilities::PointsType pointsTypeCheb = LibUtilities::eGaussGaussChebyshev;

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

      if(dimension ==2) //only quad
        {
	  
	  sol[i] =  sin(M_PI*x[i])*sin(M_PI*y[i]) +0.5;//pow((x[i]-0.5),2)+ pow((y[i]+0.5),2) + pow((z[i]-0.12),2) + demo.GetEps();//sin(2*M_PI*x[i])*cos(2*M_PI*y[i]) - 1;//pow((x[i]-0.1),2)+ pow((y[i]+0.9),2) + pow((z[i]-0.12),2) + demo.GetEps();//    ;//*sin(2*z[i]) ;//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.5;//sin(M_PI*x[i])*sin(M_PI*y[i])-2//sin(M_PI*x[i])*sin(M_PI*y[i]) +0.5;
	  funcdef = "sin(M_PI*x[i])*sin(M_PI*y[i]) +0.5";//pow((x[i]-0.5),2)+ pow((y[i]+0.5),2) + pow((z[i]-0.12),2) + demo.GetEps()";//"sin(2*M_PI*x[i])*cos(2*M_PI*y[i]) - 1";//;//="sin(2*M_PI*x[i])*cos(2*M_PI*y[i]) - 1";//sin(4*x[i])*sin(4*y[i])  -1";//*sin(2*z[i])";//"sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2 ";// "sin(2*M_PI*x[i])*sin(2*M_PI*y[i])-2";// -  exp(-1) + 1.37 ;//- 0.0602135;
	  // x[i] = ((1-0)*(x[i] + 1))/(2) + 0;
	  // y[i] = ((1-0)*(y[i] + 1))/(2) + 0; 
	  // sol[i] = pow(sin(x[i]*2),2)+pow((cos(y[i]*2)),2)+pow((y[i]*x[i]),2)- 0.5;
        
        }
      else{
	// dim = 1
	// for _/\_ function:
	// loop through all elements of the vector
	if (fmod(x[i],1.0) < 0.25 || fmod(x[i],1.0) > 0.75)
	  sol[i] = 0;
	else if( fmod(x[i],1.0) == 0.25 || fmod(x[i],1.0) == 0.75) 
	  sol[i] = 0;
	else if(fmod(x[i],1.0) > 0.25 && fmod(x[i],1.0) < 0.5)
	  sol[i] = (fmod(x[i],1.0)-0.25);
	else
	  sol[i] = 0.75-fmod(x[i],1.0);
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
      
      int dimension = E->GetShapeDimension(); 
      vector<LibUtilities::BasisKey> tmpBasisKey = demo.GetBasisKey();       
      //	cout<<"\n E->GetBasis(0)->GetNumPoints()="<<E->GetBasis(0)->GetNumPoints()<<"\n";
      LibUtilities::PointsKey pkeycheb(E->GetBasis(0)->GetNumPoints(), pointsTypeCheb);

      LibUtilities::BasisKey bkeycheb(LibUtilities::eChebyshev, 3*(E->GetBasis(0)->GetNumModes()-1),  pkeycheb);
        
      E3seg = new StdSegExp(bkeycheb);
        
      if (E3seg == nullptr)
        {
	  return 1;
        }
        
      if(dimension > 1 )
        {
	  testcoord2d = Array<OneD, Array<OneD, NekDouble> > (2); 
	  Array<OneD, Array<OneD, NekDouble> > testcoord2dtemp (2);
	      
	  int ctsz2d = 0;
	      
	  if(numedges == 4)
	    {
	      qZin = E->GetBasis(0)->GetZ();//(LibUtilities::PointsManager()[quadPointsKeyin ])->GetZ();
	      
	    
	      int sz = pow(qZin.size(),2);;
	      int szmid = pow(qZin.size()-1,2);
	      qZinmid = Array<OneD, NekDouble>(qZin.size()-1);
	      for(int i = 0; i <qZinmid.size(); i++)
		{
		  qZinmid[i] = (qZin[i]+qZin[i+1])/2;
		}
	      
	      for(int ii = 0; ii < 2; ++ii)
		{
		  testcoord2dtemp[ii] = Array<OneD, NekDouble>(sz+szmid);
		}
	      for(int ii = 0; ii < qZin.size(); ++ii)
		{
		  
		  for(int jj = 0; jj < qZin.size(); ++jj)
		    {
		      if(numedges == 4)
			{
			  testcoord2dtemp[0][ctsz2d] = qZin[ii];
			  testcoord2dtemp[1][ctsz2d] = qZin[jj];
			  ctsz2d++;
			}
		      else if(numedges == 3)
			{
			  if(qZin[ii]+qZin[jj] <=0)
			    {
			      testcoord2dtemp[0][ctsz2d] = qZin[ii];
			      testcoord2dtemp[1][ctsz2d] = qZin[jj];
			      ctsz2d++;
			    }
			}
		    }
		  
		}            
	      // add midpt grid
	      for(int ii = 0; ii < qZinmid.size(); ++ii)
		{
		  
		  for(int jj = 0; jj < qZinmid.size(); ++jj)
		    {
		      if(numedges == 4)
			{
			  testcoord2dtemp[0][ctsz2d] = qZinmid[ii];
			  testcoord2dtemp[1][ctsz2d] = qZinmid[jj];
			  ctsz2d++;
			}
		      else if(numedges == 3)
			{
			  if(qZinmid[ii]+qZinmid[jj] <=0)
			    {
			      testcoord2dtemp[0][ctsz2d] = qZinmid[ii];
			      testcoord2dtemp[1][ctsz2d] = qZinmid[jj];
			      ctsz2d++;
			    }
			}
		    } 
		  
		}
	      for(int ii = 0; ii < 2; ++ii)
		{
		  testcoord2d[ii] = Array<OneD, NekDouble>(ctsz2d);
		  Vmath::Vcopy(ctsz2d, &testcoord2dtemp[ii][0], 1, &testcoord2d[ii][0], 1);
		  
		}
	      
	    }
	  else
	    {
	      
	      int sz = pow((E->GetBasis(0)->GetZ()).size()-1+ ((E->GetBasis(1)->GetZ()).size()),2);;
	      int szmid = pow((E->GetBasis(0)->GetZ()).size() - 1+ ((E->GetBasis(1)->GetZ()).size()-1),2);

	      for(int ii = 0; ii < 2; ++ii)
		{
		  testcoord2dtemp[ii] = Array<OneD, NekDouble>( sz+szmid);
		}
	      
	      int cts2d = 0;
	      for(int ii = 0; ii < (E->GetBasis(0)->GetZ()).size(); ++ii)
		{
		  for(int jj = 0; jj < (E->GetBasis(1)->GetZ()).size();  ++jj)
		    {
		      if( (E->GetBasis(0)->GetZ())[ii] +  (E->GetBasis(1)->GetZ())[jj] <=0)
			{
			  testcoord2dtemp[0][cts2d] = (E->GetBasis(0)->GetZ())[ii];
			  testcoord2dtemp[1][cts2d] =  (E->GetBasis(1)->GetZ())[jj];
			  cts2d++;
			}
		    }
                
		}   
	      // add midpt grid
	      for(int ii = 0; ii < E->GetBasis(0)->GetZ().size() - 1; ++ii)
		{
		  for(int jj = 0; jj < E->GetBasis(1)->GetZ().size() - 1; ++jj)
		    {
		      NekDouble c1 = ((E->GetBasis(0)->GetZ())[ii] + (E->GetBasis(0)->GetZ())[ii+1])/2;
 
		      NekDouble c2 =  ((E->GetBasis(1)->GetZ())[jj] + (E->GetBasis(1)->GetZ())[jj+1])/2;
		      if(c1+c2 <= 0 )
			{
			  testcoord2dtemp[0][cts2d] = c1;
			  testcoord2dtemp[1][cts2d] = c2;
          
			  cts2d++;
			}
		    }
        
		}  

	      Array<OneD, Array<OneD, NekDouble> > temp(2);

	      temp[0] =  Array<OneD, NekDouble>(cts2d);//testcoord2dtritemp[0].size());
	      temp[1] =  Array<OneD, NekDouble>(cts2d);//testcoord2dtritemp[1].size());
	      int p = 0;

	      for(int k = 0; k < temp[0].size(); k++)
		{
		  if(testcoord2dtemp[0][k]+testcoord2dtemp[1][k]<=0)
		    {
                    
		      temp[0][p] = testcoord2dtemp[0][k];
		      temp[1][p] = testcoord2dtemp[1][k];
                                   
		      p++;
		    }
                
		}

	      for(int ii = 0; ii < 2; ++ii)
		{
		  testcoord2d[ii] = Array<OneD, NekDouble>(p);
		  Vmath::Vcopy(p, &temp[ii][0], 1, &testcoord2d[ii][0], 1);
                
		}
	      cout<<"\n1 \n\n";

	    }
	  //     Array<OneD, Array<OneD, NekDouble> >coord(2);

	  //     coord = demo.GetCoords(E);
	  //     int sz = coord[0].size();
	  //     int midptsz = sz - 1;
	  //     testcoord2dtemp[0] = Array<OneD, NekDouble>(sz+midptsz);
	  //     testcoord2dtemp[1] = Array<OneD, NekDouble>(sz+midptsz);
	      
	  //     Vmath::Vcopy(sz, coord[0], 1, testcoord2dtemp[0], 1);
	  //     Vmath::Vcopy(sz, coord[1], 1, testcoord2dtemp[1], 1);
	  //     int tot = sz;					    
	  //     for(int p = 0; p <  midptsz; p++)
	  // 	{
	  // 	      testcoord2dtemp[0][tot] = (coord[0][p]+coord[0][p+1])/2 ;
	  // 		testcoord2dtemp[1][tot] = (coord[1][p]+coord[1][p+1])/2;
	  // 		tot++;
	  // 	}
	  //     ctsz2d = midptsz+sz;
	  //   }
	  // for(int ii = 0; ii < 2; ++ii)
	  //   {
	  //     testcoord2d[ii] = Array<OneD, NekDouble>(ctsz2d);
	  //     Vmath::Vcopy(ctsz2d, &testcoord2dtemp[ii][0], 1, &testcoord2d[ii][0], 1);
	      
          //   }
	  interioreval2d = Array<OneD, NekDouble >(E->GetNcoeffs()*testcoord2d[0].size());
	  // interioreval2d0 = Array<OneD, NekDouble >(E->GetNcoeffs()*testcoord2d[0].size());
	  // interioreval2d1 = Array<OneD, NekDouble >(E->GetNcoeffs()*testcoord2d[0].size());
	  E->PhysEvalBasisGrad(testcoord2d, interioreval2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray );
	  
	  for(int k = 0; k <2; k++)
	    edgexy[k] = Array<OneD, NekDouble>(E3seg->GetBasis(0)->GetNumPoints());
	  //cout<<"E3seg->GetBasis(0)->GetNumPoints()="<<E3seg->GetBasis(0)->GetNumPoints()<<"\n";
	  Array<OneD, NekDouble> edgexytemp =   E3seg->GetBasis(0)->GetZ();
	  int totszedges = edgexytemp.size()*(E->GetNcoeffs());

	  Vxm1 = Array<OneD, NekDouble>(totszedges);
	  Vdyxm1 = Array<OneD, NekDouble>(totszedges);
	  Vdxxm1 = Array<OneD, NekDouble>(totszedges);
	  Vx1 = Array<OneD, NekDouble>(totszedges);
	  Vdxx1 = Array<OneD, NekDouble>(totszedges);
	  Vdyx1 = Array<OneD, NekDouble>(totszedges);
	  Vdxy1 = Array<OneD, NekDouble>(totszedges);
	  Vdyy1 = Array<OneD, NekDouble>(totszedges);
	  Vy1 = Array<OneD, NekDouble>(totszedges);
	  Vym1  = Array<OneD, NekDouble>(totszedges);
	  Vdxym1  = Array<OneD, NekDouble>(totszedges);
	  Vdyym1  = Array<OneD, NekDouble>(totszedges);
            
	  Vxyhyp = Array<OneD, NekDouble>(totszedges);
	  Vdxxyhyp = Array<OneD, NekDouble>(totszedges);
	  Vdyxyhyp = Array<OneD, NekDouble>(totszedges);

	  // hypotenuse (triangle) y = -x
	  edgexy[0] = edgexytemp;
	  Vmath::Smul(edgexy[0].size(), -1.0, edgexy[0], 1, edgexy[1] , 1);

	  E->PhysEvalBasisGrad(edgexy, Vxyhyp,  Vdxxyhyp, Vdyxyhyp, NullNekDouble1DArray);  
            
            
	  // left x = -1
	  edgexy[1] = edgexytemp; 
            
	  edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
	  E->PhysEvalBasisGrad(edgexy, Vxm1, Vdxxm1, Vdyxm1, NullNekDouble1DArray);

	  // bot y = -1
	  edgexy[0] = edgexytemp; 
            
	  edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
	  E->PhysEvalBasisGrad(edgexy,  Vym1, Vdxym1, Vdyym1, NullNekDouble1DArray);
	  //            cout<<"\n sz = "<<totszedges<<" edgexy[0] sz ="<<edgexy[0].size()<<" Vdxym1 sz = "<<Vdxym1.size()<<" Vym1 sz = "<<Vym1.size()<<"\n\n";
	  // right x = 1
	  edgexy[1] = edgexytemp; 
            
	  edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), 1.0); 
	  E->PhysEvalBasisGrad(edgexy, Vx1, Vdxx1, Vdyx1, NullNekDouble1DArray);
            
	  //top y = 1
	  edgexy[0] = edgexytemp; 
	    
	  //	    cout<<"\n";
	  edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), 1.0); 
	  E->PhysEvalBasisGrad(edgexy, Vy1,Vdxy1,Vdyy1, NullNekDouble1DArray);
	  timer.Stop();
	  elapsed  = timer.TimePerTest(1);
	  cout<<"\n setup phase took "<<elapsed<<"s\n\n";
	  timer.Start();
	  if(Opt_needed(coeffs))
	    {
	      timer.Stop();
	      elapsed += timer.TimePerTest(1);
	      cout<<"\n checking for -vity took "<<timer.TimePerTest(1)<<"s\n\n";

	      cout<<"\n need optimization\n\n";
	      timer.Start();
	      Do_optimize(coeffs);
	      timer.Stop();
	      elapsed += timer.TimePerTest(1);
	      cout<<"\n optimizertook "<<timer.TimePerTest(1)<<"s\n\n";

	      cout<<"\n doopt done\n verifying...\n";//exit(0);
	      timer.Start();
	      if(Opt_needed(coeffs,1))
		{
		  timer.Stop();
		  elapsed += timer.TimePerTest(1);
		  cout<<"\n Verification took "<<timer.TimePerTest(1)<<"s\n\n";
		    
		  cout<<"\n fail\n\n";
		  cout<<"\n func="<<funcdef<<"\n";
		  if(numedges == 3)
		    cout<<"\ntri\n";
		  else
		    cout<<"\n quad\n";
		  exit(0);
		}
	      else
		{
		  timer.Stop();
		  elapsed += timer.TimePerTest(1);
		  cout<<"\n Verification took "<<timer.TimePerTest(1)<<"s\n\n";
		    
		  cout<<"\n pass\n\n";

		  cout<<"\n func="<<funcdef<<"\n";;
		  if(numedges == 3)
		    cout<<"\ntri\n";
		  else
		    cout<<"\n quad\n";
		
		  phys = sol;
		}
	    }
	  else{
	    timer.Stop();
	    elapsed += timer.TimePerTest(1);
	    cout<<"\n checking for -vity took "<<timer.TimePerTest(1)<<"s\n\n";
	  
	    cout<<"\n optimizer no need\n\n";
	  }
	  //Backward transform solution to get projected values
	  E->BwdTrans(coeffs, phys);

	}
    }

  //Calculate L_inf & L_2 error
  cout << "L infinity error: \t" << E->Linf(phys, sol) << endl;
  if (stype != ePoint)
    {
      cout << "L 2 error: \t \t \t" << E->L2(phys, sol) << endl;
    }

  // if (!vm.count("diff") && stype != ePoint)
  // {
  //     //Evaluate solution at x = y = 0 and print error
  //     Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
  //     t[0] = -0.5;
  //     t[1] = -0.25;
  //     t[2] = -0.3;
  //     sol[0] = Shape_sol(t[0], t[1], t[2], order, btype, stype, false);

  //     NekDouble nsol = E->PhysEvaluate(t, phys);

  //     cout << "Error at x = (";
  //     for (int i = 0; i < dimension; ++i)
  //     {
  //         cout << t[i] << ", ";
  //     }
  //     cout << "\b\b): " << nsol - sol[0] << endl;
  // }

  // // Calculate integral of known function to test different quadrature
  // // distributions on each element.
  // for (int i = 0; i < totPoints; ++i)
  // {
  //     sol[i] = dimension == 1 ? exp(x[i]) : dimension == 2 ?
  //         exp(x[i]) * sin(y[i]) : exp(x[i] + y[i] + z[i]);
  // }

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

Array<OneD, Array<OneD, NekDouble> > vectoarr(vector<vector <NekDouble> > vec)
{
  //    cout<<"\n vec size = "<<vec.size()<<" "<<vec[0].size()<<"\n";
  Array<OneD, Array<OneD, NekDouble> > ret(vec.size());
  for (int k = 0; k < vec.size(); k++)
    {
      ret[k ] = Array<OneD, NekDouble>(vec[k].size(), vec[k].data());
      // for(int p = 0; p < vec[k].size(); p++)
      // {
      //     ret[k][p] = vec[k][p];
      // }
    }
  return ret;
}


// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret , int d)
{
  
  int modes = E3seg->GetNcoeffs();
  //    int modes1D = (3*(E->GetBasis(0)->GetNumModes()-1));//3*(E->GetBasis(0)->GetNumModes()-1);
  // In case of tri:
  // hypotenuse = 2, bot = 1, left = 2
    
  // assert d = 0 or 1
     
  if(d == 0)
    {    
      for(int k= 0; k < numedges; k++)
        {
	  ret[k] = Array<OneD, NekDouble>(modes);        
        }
    }

  if(numedges == 4) //quadrilateral
    {
        
      // bot edge = 0, y = -1
      // right edge = 1, x = 1
      // top edge = 2, y = 1
      // left edge = 3, x = -1
        
      //normal projection of der of uhats on quad edges
      if(d == 0) // from do_opt
        {  
	  //bot
	  deruhatsedges(uhats, ret[0], Vym1, Vdxym1, Vdyym1);

	  //right
	  deruhatsedges(uhats, ret[1], Vx1, Vdxx1, Vdyx1);

	  //top
	  deruhatsedges(uhats, ret[2], Vy1, Vdxy1, Vdyy1);
            
	  //left
	  deruhatsedges(uhats, ret[3], Vxm1, Vdxxm1, Vdyxm1);
            
            
        }
        
      else // d = 1 so project der to 1D 3*N space
        {
	  // uhatpqd stuff:
	  //pqeval = -(sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) + (sum(Vxm1.*Vxm1, 2)).*(Vd*c);
        
	  // bot edge
	  edgederpquhats(uhats, ret[0], E->GetNcoeffs(), V3, Vym1, Vdxym1,Vdyym1);
	  // right edge
	  edgederpquhats(uhats, ret[1],E->GetNcoeffs(), V3, Vx1, Vdxx1, Vdyx1);
	  // top edge
	  edgederpquhats(uhats, ret[2], E->GetNcoeffs(),V3, Vy1, Vdxy1, Vdyy1);
	  // left edge
	  edgederpquhats(uhats, ret[3], E->GetNcoeffs(),V3, Vxm1, Vdxxm1, Vdyxm1);
        }
    }
  else if(numedges == 1)
    {
      if(d == 1) // so project der to 3*N space
        {
	  edgederpquhats(uhats, ret[0], E->GetNcoeffs(),V3, V, Vd );
            
        }
      else
        {
	  ret[0] = uhats;
        }
    }
  else if(numedges == 3) //triangle
    {   

      if(d == 1) // so project der to 3*N space
        {
	  // bot edge
	  edgederpquhats(uhats, ret[0], E->GetNcoeffs(), V3, Vym1, Vdxym1,Vdyym1);

	  // left edge
	  edgederpquhats(uhats, ret[1], E->GetNcoeffs(), V3, Vxm1, Vdxxm1, Vdyxm1);
	  // hypotenuse     
	  edgederpquhats(uhats, ret[2], E->GetNcoeffs(), V3, Vxyhyp, Vdxxyhyp, Vdyxyhyp);

        }
      else
        {

	  // bot edge
	  //edgederuhats(uhats, ret[0], Vall, Vdxym1, Vdyym1);
	  deruhatsedges(uhats, ret[0], Vym1, Vdxym1, Vdyym1);

	  // left edge
	  //edgederuhats(uhats, ret[1], Vall, Vxm1, Vdyxm1);
	  deruhatsedges(uhats, ret[1], Vxm1, Vdxxm1, Vdyxm1);
	  // hypotenuse                 
	  //edgederuhats(uhats, ret[2], Vxyhyp,Vdxxyhyp, Vdyxyhyp);
	  deruhatsedges(uhats, ret[2], Vxyhyp,Vdxxyhyp, Vdyxyhyp);
            
        }

    }
   
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
  //    cout<<"\n1\n\n";    
  E3seg->FwdTrans(vals,ret);
  //    cout<<"\n2\n\n";
}



void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  int modes, Array<OneD, NekDouble> V3, Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2)
{
  boost::ignore_unused(Vxyd2, V3);
  int uhatstot = uhats.size();

  int totpts =  E3seg->GetBasis(0)->GetZ().size();
  Array<OneD, NekDouble> temp(totpts), temp2(modes), temp3(uhatstot);
  Array<OneD, NekDouble> pqeval(totpts);
  NekDouble v1, v2;
  for(int i = 0; i<totpts; i++)
    {
      Vmath::Vmul(modes, &Vxy[0]+i, totpts, &Vxyd0[0]+i, totpts, &temp2[0], 1);
      v1  = Vmath::Vsum(modes, temp2, 1); 
      Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);
      //cout<<" v1="<<v1;
      v2  = Vmath::Vsum(uhatstot, temp3, 1);  
      //        cout<<" v2="<<v2;
      v1 = v2*v1;
      //cout<<" v1*v2="<<v1*v2;

      // At this point,                 

      // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

      // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
      Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &temp2[0], 1);
        
      v2  = Vmath::Vsum(uhatstot, temp2, 1);  
 
      Vmath::Vmul(uhats.size(), &Vxyd0[i], totpts, &uhats[0], 1, &temp2[0], 1);

      v1= v2*Vmath::Vsum(uhats.size(), temp2, 1) - v1;
 
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
            
	  v1=  v2*Vmath::Vsum(uhats.size(), temp2, 1)-v1;

            
	  pqeval[i] += v1;
            
        }
    }
  E3seg->FwdTrans(pqeval, ret);
    
}
    

// sig = 0 -> called by opt_needed
// sig = 1 -> called by sphere_rotation
Array<OneD, Array<OneD,  NekDouble> > call_find_roots(Array<OneD,  NekDouble> &uhats , NekDouble &avgiterGD, int d, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, int sig)
{

  boost::ignore_unused(d);
 
    
  int dimension = E->GetShapeDimension(); 
  vector<vector< NekDouble> > ret(dimension);
  Array<OneD, Array<OneD, NekDouble> > retarr;//(dimension);
  NekDouble dummy;
  //find edge roots: (edges = 4 = quad, edges = 3 = tri)

  // to-do: assert 2D elements allowed numedges = 3 or 4
  if( numedges == 4)
    {
      // for 2D, structure of roots vector:
      // row 0 vector = bot edge
      // row 1 vector = right
      // row 2 vector = top
      // row 3 vector = left
      // row 4 vector = {x1,x2} inner root
      for(int ii = 0; ii < numedges; ii++)
	{
	  vector<vector<NekDouble> > tmp;
	  // size of tmp will be dim \times no. of roots
	  tmp  = (demo.find_roots(uhatsedges[ii], E, dummy, 0, sig, 0, 0, 0)) ;

	  for(int p = 0; p < tmp[0].size(); p++)
	    {
	      if(ii == 0) // bot edge
		{
		  ret[0].push_back(tmp[0][p]);
		  ret[1].push_back(-1);
                        
		}else if(ii == 1) // right edge
		{
		  ret[0].push_back(1);
		  ret[1].push_back(tmp[0][p]);
		}else if(ii == 2) // top edge
		{
		  ret[0].push_back(tmp[0][p]);
		  ret[1].push_back(1);
                    
                    
		}else //if(ii == 3) // left edge
		{
		  ret[0].push_back(-1);
		  ret[1].push_back(tmp[0][p]);

		}
                    
	    }
                
	}
      // add 4 corners
      ret[0].push_back(-1);
      ret[1].push_back(-1);
   
      ret[0].push_back(1);
      ret[1].push_back(1);
       
      ret[0].push_back(-1);
      ret[1].push_back(1);
        
      ret[0].push_back(1);
      ret[1].push_back(-1);
    }
  else
    {
      // triangle edge roots code
      // structure of roots vector:
      // row 0 vector = bot edge
      // row 1 vector = left
      // row 2 vector = hypotenuse
      // row 3 vector = {x1,x2} inner root
      for(int ii = 0; ii < numedges; ii++)
	{

	  vector<vector<NekDouble> > tmp;
            
	  // size of tmp will be dim \times no. of roots
                
	  tmp  = (demo.find_roots(uhatsedges[ii], E, dummy, 0, sig, 0, 0, 0)) ;
            
	  for(int p = 0; p < tmp[0].size(); p++)
	    {
	      if(ii == 0) // bot edge
		{
		  ret[0].push_back(tmp[0][p]);
		  ret[1].push_back(-1);
                        
		}else if(ii == 1) // left edge
		{
		  ret[0].push_back(-1);
		  ret[1].push_back(tmp[0][p]);
		}
	      else //hypt
		{

		  ret[0].push_back(tmp[0][p]);
		  ret[1].push_back(-tmp[0][p]);
                        
		  ret[0].push_back(-tmp[0][p]);
		  ret[1].push_back(tmp[0][p]);
                        
		}
	    }
		
	  // ret[0].push_back(0);
	  // ret[1].push_back(0);
	}
      // add 3 corners
      ret[0].push_back(-1);
      ret[1].push_back(-1);
   
      ret[0].push_back(-1);
      ret[1].push_back(1);
        
      ret[0].push_back(1);
      ret[1].push_back(-1);
   
    }

  //find interior roots:
  vector<vector<NekDouble> > tmp;

  if(numedges == 3)
    {

      demo.interioreval2dtri = interioreval2d;
      demo.testcoord2dtri = testcoord2d;
      tmp =  (demo.find_roots(uhats, E, avgiterGD,  0, sig, 1, 0, 1)) ;
    }
  else
    {
      demo.testcoord2dquad = testcoord2d;
      demo.interioreval2dquad = interioreval2d;
      tmp =  (demo.find_roots(uhats, E, avgiterGD,  0, sig, 1, 0, 0)) ;
    }
  if(tmp[0].size() > 0)
    {
      for(int p = 0; p <dimension; p++)
	{
	  ret[p].push_back(tmp[p][0]);
	}
    }
    
  if(ret[0].size() > 0)
    retarr = vectoarr(ret);

  return retarr;
}


int Opt_needed(Array<OneD, NekDouble> uhats, int flag)
{
  int totModes = uhats.size();
  Array<OneD,  Array<OneD, NekDouble> > rootsarr;
  Array<OneD, Array<OneD, NekDouble> > edgeuhats  (numedges);

  NekDouble avgiterGD;
  project_edges(uhats, edgeuhats);
    
  //find roots of der of uhats, that is why d = 0
  rootsarr = call_find_roots(uhats, avgiterGD, 0, edgeuhats);
  // evaluate ortho basis at roots
  // evalBasisRoots is flattened basis eval matrix

  int     evalsz = (rootsarr[0].size())*totModes;
  Array<OneD, NekDouble> evalBasisRoots(evalsz);//, vals(roots[0].size()); 

  E->PhysEvalBasisGrad(rootsarr, evalBasisRoots, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
  
  Array<OneD, NekDouble> tmp(rootsarr.size()),temp(totModes) ;

  NekDouble minv = numeric_limits<double>::infinity();
  int idx = 0;
  for(int k = 0; k < rootsarr[0].size(); ++k)
    {
      NekDouble tempv;
      Vmath::Vmul(totModes, &evalBasisRoots[k], rootsarr[0].size(), &uhats[0], 1, &temp[0], 1);
      tempv = Vmath::Vsum(temp.size(), temp, 1);
      if(tempv < minv)
        {
	  minv = tempv;
	  idx = k;
        }
    }
  cout<<"\n minv = "<<minv<<" at "<<  rootsarr[0][idx]<<" "<<rootsarr[1][idx];
    
  if(minv < 0.0 && abs(minv)>1e-9)
    {
      if(flag == 0)
        {
	  startval = minv;
	  startcoordx = rootsarr[0][idx];
	  startcoordy = rootsarr[1][idx];
     
        }
      return 1;
    }
  return 0;
}


void Do_optimize(Array<OneD, NekDouble> &uhats)
{

  int dim = E->GetShapeDimension();
  double inf = numeric_limits<double>::infinity();
        
  int N1 = uhats.size();
  vector<Array<OneD,NekDouble> > d;
 
  d.push_back(uhats);
  int counter = 0;
      
  //    int NC = ck[1].size();          // number of constraints
  vector<double> tols;         // constraint specific tolerances

  int niter = 1e3;
  //    vector< vector<NekDouble> > optima;
  Array<OneD, Array<OneD,NekDouble> > optima;

  //    int NC = 1; //number of constraints, only positivity for now
  tols.push_back(1e-16);
  Array<OneD, NekDouble>  pqvalcoords(dim+1), xastarr(dim), utemp(N1), wsp1;   
    
  Array<OneD, Array<OneD, NekDouble> > Pf(numedges);

  //Array<OneD, Array<OneD, NekDouble> > surfaceuhats  (numsurfaces);
  for(int k= 0; k < numedges; k++)
    {
      Pf[k] = Array<OneD, NekDouble>(E3seg->GetNcoeffs());        
    }
  NekDouble avgiterGD;
  NekDouble pqval;
  while (counter <= niter)
    {
      NekDouble avgiterhold;
      pqval = inf;
      utemp = d.back();
      //pq = @(xx) cfun(xx) .* (constraints(1,jj).*Vc(xx)*uhats);  
      project_edges(utemp, Pf, 1);
      //cout<<"\n counter:<<"<<counter<<"\n";
      optima = (call_find_roots(utemp, avgiterhold, 0, Pf, 1));
      avgiterGD += avgiterhold;
      //cout<<"\n optima sz = "<<optima[0].size()<<"\n";
      wsp1=Array<OneD, NekDouble>(optima[0].size());
        
      Array<OneD, NekDouble> Vtmp(optima[0].size() * E->GetNcoeffs());
      E->PhysEvalBasisGrad(optima, Vtmp, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

      demo.pq(utemp, optima, Vtmp, pqvalcoords, wsp1);

      // cout<<"\n fvals: ";
      //  for(int kk = 0 ; kk<wsp1.size(); kk++)
      //      cout<<wsp1[kk]<<" ";
      //  cout<<"\n";

      if (pqvalcoords[0] < pqval)
        {
	  for(int k = 0; k  < dimension; k++)
            {
	      xastarr[k] = pqvalcoords[k+1];
            }
	  pqval = pqvalcoords[0];
        }
        
      if(xastarr.size() == 2)
        {
	  cout<<"\n at counter="<<counter<<" min val="<<pqval<<" xast ="<<xastarr[0]<<" "<<xastarr[1];
        }
      else //if(xastarr.size() == 1)
        {
	  cout<<"\n at counter="<<counter<<" min val="<<pqval<<" xast ="<<xastarr[0]<<"";
        
        }
      // If minimum is non-negative, we're done
      if (pqval >= -tols.at(0))
        {
	  break;
        }
        
      vector<NekDouble> Vastsq;
      vector<NekDouble> Vast;
        
      NekDouble vastsqsum;

      for( int ii = 0; ii < N1; ii++)
        {
	  Vast.push_back(E->PhysEvaluateBasis(xastarr, ii));
	  Vastsq.push_back(Vast[ii]*Vast[ii]);
        
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
  cout<<"\n sphere_rotation took "<<counter<<"iterations\n startv = "<<startval<<" at startcoords("<<startcoordx<<","<<startcoordy<<")     GD iter stepszie ="<<demo.GetTol()<<" maxiters = "<<demo.GetMaxIter()<<" o = "<<E->GetBasis(0)->GetNumModes()<<" p = "<<E->GetBasis(0)->GetNumPoints()<< " eps = "<<demo.GetEps()<<" avgiterGD = "<<(avgiterGD)/counter<<" \t";;
  if(numedges == 4)
    cout<<"\n quad\n";
  else
    cout<<"\n tri\n";

  uhats = d.back();
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
