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

// NekDouble Shape_sol(NekDouble x, NekDouble y, NekDouble z, vector<int> order,
//                     vector<BasisType> btype, ShapeType stype, bool diff);

//Modification to deal with exact solution for diff. Return 1 if integer < 0.
// static double pow_loc(const double val, const int i)
// {
//   return (i < 0) ? 1.0 : pow(val, i);
// }

//Array<OneD, Array<OneD, NekDouble> > vectoarr(vector<vector< NekDouble> > vec);


////declare Do_optimize
//void Do_optimize(Array<OneD, NekDouble> &uhats);

//declare find_roots, flag is 0 for confederate matrix approach
// Array<OneD, Array<OneD, NekDouble> >  find_roots(Array<OneD, NekDouble> &uhats, int d = 0, int flag = 0,  int sig = 0, int surfflag = 0);

// declare caller routine to find_roots
// flag = 0 -> opt_needed calls
// flag = 1 -> sphere_rot calls
// Array<OneD, Array<OneD,  NekDouble> >  call_find_roots(Array<OneD, NekDouble> &uhatsall, NekDouble &avgiterGD, int d = 0, Array< OneD, Array<OneD, NekDouble> >&uhatsedges =NullNekDoubleArrayofArray , int flag = 0 );
Array<OneD,NekDouble>  call_find_roots(Array<OneD,  NekDouble> &uhats1 , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats , NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime);


//declare Opt_needed
//int Opt_needed(Array<OneD, NekDouble> uhats, int flag = 0);


// for 2D elements to get uhats at edges
// when d = 1, its the uhatpqd case from do_optimize
//void project_edges( Array<OneD, NekDouble> uhats,Array<OneD, Array<OneD, NekDouble> >&ret );


// uhatpqd stuff: 
// called by project_edges if d = 1 is passed in project_edge
//void edgederpquhats(Array<OneD, NekDouble> &uhats1, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> &wsp,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2);

//void deruhatsedges(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble >&ret,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2 = NullNekDouble1DArray );

//vector< vector<NekDouble> > gradient_descent(Array<OneD, NekDouble> uhats, int sig, int surfid = 0 );

//string funcdef;
Array<OneD, Array<OneD, NekDouble> > storage2d;

// Colleague matrix
int numedges, numsurfaces;
//Array<OneD, Array<OneD, NekDouble> > C;
StdExpansion *E;
// StdExpansion *E3seg;
// Array<OneD, NekDouble> qZin; //inner point grid for 2D rootfinding
// Array<OneD, NekDouble> qZinmid; //inner mid point grid for 2D rootfinding
// Array<OneD, NekDouble> qWin; //inner point grid for 2D rootfinding
// Array<OneD, NekDouble> qx;
// Array<OneD, NekDouble> qw;

// Array<OneD, NekDouble> V;
// Array<OneD, NekDouble> Vd;
// Array<OneD, NekDouble> V3;
// Array<OneD, NekDouble> Vxm1;
// Array<OneD, NekDouble> Vdxxm1;
// Array<OneD, NekDouble> Vdyxm1;
// Array<OneD, NekDouble> Vx1 ;
// Array<OneD, NekDouble> Vdyx1;
// Array<OneD, NekDouble> Vdxx1;
// Array<OneD, NekDouble> Vdxy1;
// Array<OneD, NekDouble> Vdyy1;
// Array<OneD, NekDouble> Vy1;
// Array<OneD, NekDouble> Vym1;
// Array<OneD, NekDouble> Vdxym1;
// Array<OneD, NekDouble> Vdyym1;
// Array<OneD, NekDouble > qWinarr;

// // for triangle edges root-finding:
// Array<OneD, NekDouble> Vxyhyp;
// Array<OneD, NekDouble> Vdxxyhyp;
// Array<OneD, NekDouble> Vdyxyhyp;

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
  storage2d = E->GetPhysEvaluateStorage(); 
  stype = E->DetShapeType();
      
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
 
    
  switch(E->DetShapeType())
    {

    case LibUtilities::eTriangle:

      numedges = 3;
      numsurfaces = 1;

      demo.testcoord2dtqpts = demo.GetCoords(E);
      demo.testcoord2dtqmidpts = demo.GetQuadratureMidCoords( E, demo.testcoord2dtqpts);
      demo.testcoord2dtlattice = demo.GetLatticeCoords(demo.testcoord2dtqpts, demo.testcoord2dtqmidpts);     demo.interioreval2dtqmidpts = E->PhysEvaluateBasis( demo.testcoord2dtqmidpts, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
 
      break;
    
    case LibUtilities::eQuadrilateral:

      demo.testcoord2dqqpts = demo.GetCoords(E);
      demo.testcoord2dqqmidpts = demo.GetQuadratureMidCoords( E, demo.testcoord2dqqpts);
      demo.testcoord2dqlattice = demo.GetLatticeCoords( demo.testcoord2dqqpts, demo.testcoord2dqqmidpts);
      demo.interioreval2dqqmidpts = E->PhysEvaluateBasis( demo.testcoord2dqqmidpts, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
 
      numedges = 4;
      numsurfaces = 1;
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
	  sol[i] = -2*x[i] + pow((x[i] + 0.6),3) + (pow(y[i],2) - 0.2);
	  //f2
	  //floor(x[i]<=0 && y[i] <=0)+0.1;
	  //f1
	  //-1*((sin((x[i]+0.5*M_PI)))*(cos((y[i]))))+0.95;
	  //-2*x[i] + pow((x[i] + 0.6),3) + (pow(y[i],2) - 0.2);
	  // funcdef = tostring(sol);
        }
      // else{
      // 	// dim = 1
      // 	// for _/\_ function:
      // 	// loop through all elements of the vector
      // 	if (fmod(x[i],1.0) < 0.25 || fmod(x[i],1.0) > 0.75)
      // 	  sol[i] = 0;
      // 	else if( fmod(x[i],1.0) == 0.25 || fmod(x[i],1.0) == 0.75) 
      // 	  sol[i] = 0;
      // 	else if(fmod(x[i],1.0) > 0.25 && fmod(x[i],1.0) < 0.5)
      // 	  sol[i] = (fmod(x[i],1.0)-0.25);
      // 	else
      // 	  sol[i] = 0.75-fmod(x[i],1.0);
      // }
    }

  Array<OneD, NekDouble> phys(totPoints);
  Array<OneD, NekDouble> coeffs((unsigned) E->GetNcoeffs());
 
  //Project onto expansion
  E->FwdTrans(sol, coeffs);

  //Backward transform solution to get projected values
  E->BwdTrans(coeffs, phys);

  // LibUtilities::Timer     timer;
  // NekDouble elapsed       = 0.0;


  // check for -ve values and apply opt if necessary
  if (vm.count("optm"))
    {
      // timer.Start();
      
      // // int dimension = E->GetShapeDimension(); 
      // vector<LibUtilities::BasisKey> tmpBasisKey = demo.GetBasisKey();       

      // LibUtilities::PointsKey pkeycheb(E->GetBasis(0)->GetNumPoints(), pointsTypeCheb);

      // LibUtilities::BasisKey bkeycheb(LibUtilities::eChebyshev, 3*(E->GetBasis(0)->GetNumModes()-1),  pkeycheb);
        
      // E3seg = new StdSegExp(bkeycheb);
        
      // if (E3seg == nullptr)
      //   {
      // 	  return 1;
      //   }
        	  
      // for(int k = 0; k <2; k++)
      // 	{
      // 	  edgexy[k] = Array<OneD, NekDouble>(E3seg->GetBasis(0)->GetNumPoints());
	  
      // 	  edgeptsin[k] =   Array<OneD, NekDouble>(edgexy[k]);
      // 	}
      // Array<OneD, NekDouble> edgexytemp =   E3seg->GetBasis(0)->GetZ();
      // int totszedges = edgexytemp.size()*(E->GetNcoeffs());

      // Vxm1 = Array<OneD, NekDouble>(totszedges);
      // Vdyxm1 = Array<OneD, NekDouble>(totszedges);
      // Vdxxm1 = Array<OneD, NekDouble>(totszedges);
      // Vx1 = Array<OneD, NekDouble>(totszedges);
      // Vdxx1 = Array<OneD, NekDouble>(totszedges);
      // Vdyx1 = Array<OneD, NekDouble>(totszedges);
      // Vdxy1 = Array<OneD, NekDouble>(totszedges);
      // Vdyy1 = Array<OneD, NekDouble>(totszedges);
      // Vy1 = Array<OneD, NekDouble>(totszedges);
      // Vym1  = Array<OneD, NekDouble>(totszedges);
      // Vdxym1  = Array<OneD, NekDouble>(totszedges);
      // Vdyym1  = Array<OneD, NekDouble>(totszedges);
            
      // Vxyhyp = Array<OneD, NekDouble>(totszedges);
      // Vdxxyhyp = Array<OneD, NekDouble>(totszedges);
      // Vdyxyhyp = Array<OneD, NekDouble>(totszedges);

      // // hypotenuse (triangle) y = -x
      // edgexy[0] = edgexytemp;
      // Vmath::Smul(edgexy[0].size(), -1.0, edgexy[0], 1, edgexy[1] , 1);

      // Vxyhyp = E->PhysEvaluateBasis(edgexy, storage2d,  Vdxxyhyp, Vdyxyhyp, NullNekDouble1DArray);  
            
            
      // // left x = -1
      // edgexy[1] = edgexytemp; 
            
      // edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
      // Vxm1 = E->PhysEvaluateBasis(edgexy, storage2d, Vdxxm1, Vdyxm1, NullNekDouble1DArray);

      // // bot y = -1
      // edgexy[0] = edgexytemp; 
            
      // edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
      // Vym1 = E->PhysEvaluateBasis(edgexy,  storage2d, Vdxym1, Vdyym1, NullNekDouble1DArray);

      // // right x = 1
      // edgexy[1] = edgexytemp; 
            
      // edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), 1.0); 
      // Vx1 = E->PhysEvaluateBasis(edgexy, storage2d, Vdxx1, Vdyx1, NullNekDouble1DArray);
            
      // //top y = 1
      // edgexy[0] = edgexytemp; 
	    
      // //	    cout<<"\n";
      // edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), 1.0); 
      // Vy1 = E->PhysEvaluateBasis(edgexy, storage2d,Vdxy1,Vdyy1, NullNekDouble1DArray);
      // timer.Stop();
      // elapsed  = timer.TimePerTest(1);
      // cout<<"\n setup phase took "<<elapsed<<"s\n\n";
      // cout<<"\n checking for -vity took "<<timer.TimePerTest(1)<<"s\n\n";
      NekDouble v1;
      Array<OneD, Array<OneD, NekDouble> > dummy, surf(1);
      surf[0] = Array<OneD, NekDouble>(coeffs.size());
      Vmath::Vcopy(coeffs.size(), &coeffs[0], 1, &surf[0][0], 1); 
      Array<OneD, NekDouble> optima = call_find_roots(coeffs, v1, dummy, surf,  v1, v1, v1, v1);
	  
    }


  //Calculate L_inf & L_2 error
  // cout << "L infinity error: \t" << E->Linf(phys, sol) << endl;
  // if (stype != ePoint)
  //   {
  //     cout << "L 2 error: \t \t \t" << E->L2(phys, sol) << endl;
  //   }

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

  // NekDouble exact = 0.0;
  // switch(stype)
  //   {
  //   case eSegment:
  //     exact = M_E - 1.0 / M_E;
  //     break;
  //   case eTriangle:
  //     exact = -0.5 * (sin(1.0) + cos(1.0) + M_E * M_E *
  // 		      (sin(1.0) - cos(1.0))) / M_E;
  //     break;
  //   case eQuadrilateral:
  //     exact = 2.0 * (M_E - 1.0 / M_E) * sin(1.0);
  //     break;
  //   case eTetrahedron:
  //     exact = 1.0 / M_E - 1.0 / M_E / M_E / M_E;
  //     break;
  //   case ePrism:
  //     exact = M_E - 1.0 / M_E / M_E / M_E;
  //     break;
  //   case ePyramid:
  //     exact = - 1.0 / M_E / M_E / M_E - 4.0 / M_E + M_E;
  //     break;
  //   case eHexahedron:
  //     exact = pow((M_E * M_E - 1.0) / M_E, 3.0);
  //     break;
  //   default:
  //     ASSERTL0(false, "Exact solution not known.");
  //     break;
  //   }
  // std::cout << "Integral error: " << fabs(exact - E->Integral(sol))
  // 	    << std::endl;
  // return 0;
    
}

// Array<OneD, Array<OneD, NekDouble> > vectoarr(vector<vector <NekDouble> > vec)
// {
//   //    cout<<"\n vec size = "<<vec.size()<<" "<<vec[0].size()<<"\n";
//   Array<OneD, Array<OneD, NekDouble> > ret(vec.size());
//   for (int k = 0; k < vec.size(); k++)
//     {
//       ret[k ] = Array<OneD, NekDouble>(vec[k].size(), vec[k].data());
//       // for(int p = 0; p < vec[k].size(); p++)
//       // {
//       //     ret[k][p] = vec[k][p];
//       // }
//     }
//   return ret;
// }


// // when d = 1, its the uhatpqd case from do_optimize
// // void project_edges( Array<OneD, NekDouble>uhats1,    Array<OneD, Array<OneD, NekDouble> >&ret)
// // {
    
// //     Array<OneD, NekDouble> wsp(uhats1.size());


// //     //    int E3segncoeff = E3seg->GetNcoeffs();

// //     // for(int k= 0; k < numedges; k++)
// //     //   {
// //     // 	ret[k] = Array<OneD, NekDouble>(E3segncoeff);//(3*(E->GetBasis(0)->GetNumModes()-1));        
      
// //     //   }
// //     // cout<<"\nhere2\n\n";

// //     if(numedges == 4) // quad
// //       {
// // 	// bot edge
// // 	edgederpquhats(uhats1, ret[0], wsp, Vym1, Vdxym1, Vdyym1, NullNekDouble1DArray);
// // 	// right edge
// // 	edgederpquhats(uhats1, ret[1],wsp, Vx1, Vdxx1, Vdyx1, NullNekDouble1DArray);
// // 	// top edge
// // 	edgederpquhats(uhats1, ret[2], wsp, Vy1, Vdxy1, Vdyy1, NullNekDouble1DArray);
// // 	// left edge
// // 	edgederpquhats(uhats1, ret[3], wsp, Vxm1, Vdxxm1, Vdyxm1, NullNekDouble1DArray);

// //       }
// //     if(numedges == 3) // tri
// //       {
// // 	// bot edge  
// // 	edgederpquhats(uhats1, ret[0], wsp,  Vym1, Vdxym1, Vdyym1,NullNekDouble1DArray);

// 	// left edge
// 	edgederpquhats(uhats1, ret[1], wsp, Vxm1, Vdxxm1, Vdyxm1, NullNekDouble1DArray);

// 	// hypto
// 	edgederpquhats(uhats1, ret[2], wsp, Vxyhyp, Vdxxyhyp, Vdyxyhyp, NullNekDouble1DArray);
//       }
//     // if(numedges == 6) // tet
//     //   {

//     //     if(d == 1)
//     //       {

//     // 	  // edge front left (AD) (x = -1) (y = -1)
//     // 	  edgederpquhats(uhats, ret[0], E->GetNcoeffs(),  Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);
                        
//     // 	  //edge front hypt (DB) (y = -1) (z = -x) 
//     // 	  edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vym1xmz, Vdxym1xmz, Vdyym1xmz, Vdzym1xmz  );
            
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
//}

   


// void deruhatsedges(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz, Array<OneD, NekDouble> Vdxyz, Array<OneD, NekDouble> Vxdyz, Array<OneD, NekDouble> Vxydz)
// {
//   boost::ignore_unused(Vxyz, Vxydz);
//   int totpts =   E3seg->GetBasis(0)->GetZ().size();
//   int modes = uhats.size();
//   NekDouble v1;   
//   Array<OneD, NekDouble> vals(totpts), hold(totpts);
//   Array<OneD, NekDouble> temp(uhats.size());  
    
//   for(int k = 0; k < totpts; k++)
//     {
//       Vmath::Vmul(uhats.size(), &Vdxyz[k], totpts, &uhats[0], 1, &temp[0], 1);

//       v1  = Vmath::Vsum(modes, temp, 1);
//       Vmath::Vmul(uhats.size(), &Vxdyz[k], totpts, &uhats[0], 1, &temp[0], 1);
//       v1  = v1 + Vmath::Vsum(modes, temp, 1);  
        
//     }    
//   //    cout<<"\n1\n\n";    
//   E3seg->FwdTrans(vals,ret);
//   //    cout<<"\n2\n\n";
// }



// // void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  int modes, Array<OneD, NekDouble> V3, Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2)
// // {
// //   boost::ignore_unused(Vxyd2, V3);
// //   int uhatstot = uhats.size();

// //   int totpts =  E3seg->GetBasis(0)->GetZ().size();
// //   Array<OneD, NekDouble> temp(totpts), temp2(modes), temp3(uhatstot);
// //   Array<OneD, NekDouble> pqeval(totpts);
// //   NekDouble v1, v2;
// //   for(int i = 0; i<totpts; i++)
// //     {
// //       Vmath::Vmul(modes, &Vxy[0]+i, totpts, &Vxyd0[0]+i, totpts, &temp2[0], 1);
// //       v1  = Vmath::Vsum(modes, temp2, 1); 
// //       Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);
// //       //cout<<" v1="<<v1;
// //       v2  = Vmath::Vsum(uhatstot, temp3, 1);  
// //       //        cout<<" v2="<<v2;
// //       v1 = v2*v1;
// //       //cout<<" v1*v2="<<v1*v2;

//       // At this point,                 

//       // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

//       // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
//       Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &temp2[0], 1);
        
//       v2  = Vmath::Vsum(uhatstot, temp2, 1);  
 
//       Vmath::Vmul(uhats.size(), &Vxyd0[i], totpts, &uhats[0], 1, &temp2[0], 1);

//       v1= v2*Vmath::Vsum(uhats.size(), temp2, 1) - v1;
 
//       pqeval[i] = v1;

//       if(dimension > 1)
//         {
// 	  Vmath::Vmul(modes, &Vxy[0]+i, totpts, &Vxyd1[0]+i, totpts, &temp2[0], 1);
// 	  v1  = Vmath::Vsum(modes, temp2, 1); 
// 	  Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);
            
// 	  v2  = Vmath::Vsum(uhatstot, temp3, 1);  
// 	  v1 = v2*v1;
            
// 	  // At this point,                 
            
// 	  // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point
            
// 	  // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
// 	  Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &temp2[0], 1);
            
// 	  v2  = Vmath::Vsum(uhatstot, temp2, 1);  
            
// 	  Vmath::Vmul(uhats.size(), &Vxyd1[i], totpts, &uhats[0], 1, &temp2[0], 1);
            
// 	  v1=  v2*Vmath::Vsum(uhats.size(), temp2, 1)-v1;

            
// 	  pqeval[i] += v1;
            
//         }
//     }
//   E3seg->FwdTrans(pqeval, ret);
    
// }

// void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> &wsp,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2)
//   {
//     int uhatstot = uhats.size();
//     int totpts = edgeptsin[0].size();//E3seg->GetTotPoints();

//     Array<OneD, NekDouble> pqeval(totpts);

//     NekDouble v1, v2;
//     for(int i = 0; i<totpts; i++)
//       {

// 	Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd0[0]+i, totpts, &wsp[0], 1);
// 	v1  = Vmath::Vsum(uhatstot, wsp, 1); 
// 	Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &wsp[0], 1);

// 	v2  = Vmath::Vsum(uhatstot, wsp, 1);  

// 	v1 = v2*v1;

// 	// At this point,                 

// 	// v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

// 	// calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
// 	Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &wsp[0], 1);

 
        
// 	v2  = Vmath::Vsum(uhatstot, wsp, 1);  
 
// 	Vmath::Vmul(uhats.size(), &Vxyd0[i], totpts, &uhats[0], 1, &wsp[0], 1);

// 	v1= v2*Vmath::Vsum(uhats.size(), wsp, 1)- v1;
  
// 	pqeval[i] = v1;
 
// 	if(dimension > 1)
// 	  {
// 	    Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd1[0]+i, totpts, &wsp[0], 1);
// 	    v1  = Vmath::Vsum(uhatstot, wsp, 1); 
// 	    Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &wsp[0], 1);

// 	    v2  = Vmath::Vsum(uhatstot, wsp, 1);  

// 	    v1 = v2*v1;

// 	    // At this point,                 

// 	    // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

// 	    // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
// 	    Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &wsp[0], 1);

 
        
// 	    v2  = Vmath::Vsum(uhatstot, wsp, 1);  
 
// 	    Vmath::Vmul(uhats.size(), &Vxyd1[i], totpts, &uhats[0], 1, &wsp[0], 1);

// 	    v1=  v2*Vmath::Vsum(uhats.size(), wsp, 1) - v1;
 
// 	    pqeval[i] += v1;
 
// 	  }
// 	if(dimension == 3)
// 	  {
// 	    Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd2[0]+i, totpts, &wsp[0], 1);
// 	    v1  = Vmath::Vsum(uhatstot, wsp, 1); 
// 	    Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &wsp[0], 1);

// 	    v2  = Vmath::Vsum(uhatstot, wsp, 1);  

// 	    v1 = v2*v1;

// 	    // At this point,                 

// 	    // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

// 	    // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
// 	    Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &wsp[0], 1);

 
        
// 	    v2  = Vmath::Vsum(uhatstot, wsp, 1);  
 
// 	    Vmath::Vmul(uhats.size(), &Vxyd2[i], totpts, &uhats[0], 1, &wsp[0], 1);

// 	    v1= v2*Vmath::Vsum(uhats.size(), wsp, 1) - v1;
 
// 	    pqeval[i] += v1;
 
// 	  }
//       }

//     E3seg->FwdTrans(pqeval, ret);
    
//   }




Array<OneD,  NekDouble> call_find_roots(Array<OneD,  NekDouble> &uhats1 , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats, NekDouble &minv, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime)
{
  
  boost::ignore_unused(surfaceuhats, roots2dtime, roots3dtime, roots1dtime, minv, uhatsedges);
  int dimension = E->GetShapeDimension(); 
  Array<OneD, NekDouble> coords(dimension);
  Array<OneD, Array<OneD, NekDouble > >  rethold(dimension);
  if(numsurfaces == 1)
    {
      if(stype == 3)	  //if tri
	rethold = demo.find_roots(uhats1, E, storage2d,  avgiterGD, 0, 1, 0 , 1);
      else if(stype == 4)//else quad
	rethold = demo.find_roots(uhats1, E, storage2d,  avgiterGD, 0, 1, 0 , 0);
    }
  coords[0] = rethold[0][0];
  coords[1] = rethold[1][0];
  return coords;	
}


// int Opt_needed(Array<OneD, NekDouble> uhats, int flag)
// {
//   int totModes = uhats.size();
//   Array<  OneD, Array<OneD, NekDouble> > rootsarr(2);
//   Array<OneD, Array<OneD, NekDouble> > edgeuhats  (numedges);

//   NekDouble avgiterGD, minv, roots1dtime, roots2dtime, roots3dtime;
//   project_edges(uhats, edgeuhats);

//   Array<OneD, NekDouble> rootshold;
//   //find roots of der of uhats, that is why d = 0
//   rootshold = call_find_roots(uhats, avgiterGD, edgeuhats,  NullNekDoubleArrayofArray,minv, roots1dtime, roots2dtime, roots3dtime);
//   // evaluate ortho basis at roots
//   // evalBasisRoots is flattened basis eval matrix
//   rootsarr[0] = Array<OneD, NekDouble>(1,rootshold[0]);
//   rootsarr[1] = Array<OneD, NekDouble>(1,rootshold[1]);
//   int evalsz = totModes;
//   Array<OneD, NekDouble> evalBasisRoots(evalsz);//, vals(roots[0].size()); 
//   evalBasisRoots = E->PhysEvaluateBasis(rootsarr, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
  
//   Array<OneD, NekDouble> tmp(rootsarr.size()),temp(totModes) ;

//    minv = numeric_limits<double>::infinity();
//   int idx = 0;
//   for(int k = 0; k < rootsarr[0].size(); ++k)
//     {
//       NekDouble tempv;
//       Vmath::Vmul(totModes, &evalBasisRoots[k], rootsarr[0].size(), &uhats[0], 1, &temp[0], 1);
//       tempv = Vmath::Vsum(temp.size(), temp, 1);
//       if(tempv < minv)
//         {
// 	  minv = tempv;
// 	  idx = k;
//         }
//     }
//   //cout<<"\n minv = "<<minv<<" at "<<  rootsarr[0][idx]<<" "<<rootsarr[1][idx];
    
//   if(minv < 0.0 && abs(minv)>1e-9)
//     {
//       if(flag == 0)
//         {
// 	  startval = minv;
// 	  startcoordx = rootsarr[0][idx];
// 	  startcoordy = rootsarr[1][idx];
     
//         }
//       return 1;
//     }
//   return 0;
// }



// // void Do_optimize(Array<OneD, NekDouble> &uhats)
// // {

// //     int dim = E->GetShapeDimension();
// //     double inf = numeric_limits<double>::infinity();
    
// //     int N1 = uhats.size();
// //     vector<Array<OneD,NekDouble> > d;
 
// //     d.push_back(uhats);
// //     int counter = 0;
      
    
// //     vector<double> tols;
// //     int niter = 1e3;
// //     NekDouble minv;
// //     Array<OneD,NekDouble> optima(dim);
// //     tols.push_back(1e-12);
    
// //     Array<OneD, NekDouble>  pqvalcoords(dim+1), xastarr(dim), utemp(N1), wsp1(1);   
    

// //     Array<OneD, Array<OneD, NekDouble> > surfaceuhats  (numsurfaces), tmpcoord(dim);
// //     Array<OneD, Array<OneD, NekDouble> > Pf(numedges);
// //     for(int k= 0; k < numedges; k++)
// //     {
            
// //       Pf[k] = Array<OneD, NekDouble>(E3seg->GetNcoeffs());
// //     }
// //     for(int k = 0; k < numsurfaces; k++)
// //       {
// // 	surfaceuhats[k] = Array<OneD, NekDouble>(E->GetNcoeffs());
// //       }
// //     NekDouble pqval, timeprojectsurf = 0.0, timeprojectedges = 0.0;
 //     NekDouble avgiterGD = 0.0, avgiterhold = 0.0, roots1dtimehold = 0.0, roots2dtimehold = 0.0, roots3dtimehold = 0.0;
//     Array<OneD, NekDouble> Vtmp(  E->GetNcoeffs());

//     while (counter <= niter)
//     {
//       NekDouble roots1dtime = 0.0, roots2dtime = 0.0, roots3dtime = 0.0 ;
//       pqval = inf;
//       utemp = d.back();
	
//       if (counter > -1)
// 	{
// 	  project_edges(utemp, Pf);
// 	  optima = call_find_roots(utemp, avgiterhold, Pf, surfaceuhats,  minv, roots1dtime, roots2dtime, roots3dtime);
// 	  roots1dtimehold += roots1dtime;
// 	  roots2dtimehold += roots2dtime;
// 	  roots3dtimehold += roots3dtime;
// 	  avgiterGD += avgiterhold;
// 	  cout<<"\n optima:\n";
// 	  for(int pp = 0; pp < optima.size(); pp++)
// 	  {
// 	    cout<<optima[pp]<<" ";
// 	  }
// 	  cout<<"\n\n";
// 	  exit(0);
// 		 // Vtmp is evaluation of basis fn at optima
// 	  for(int k = 0; k < dim; k++)
// 	    {
// 	      tmpcoord[k] = Array<OneD, NekDouble>(1, optima[k]);
// 	    }
	
// 	  //	t.Start();

// 	  Vtmp = E->PhysEvaluateBasis(tmpcoord, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        
// 	  //	  demo.pq(utemp, tmpcoord, Vtmp, pqvalcoords, wsp1);

// 	}// end if(counter > 0)
//       else //if(counter == 0)
// 	{

// 	  //	  cout<<"\n startcoo="<<startcoordx<<", "<<startcoordy<<"\n";
// 	  for(int k = 0; k < dim; k++)
// 	    {
// 	      tmpcoord[k] = Array<OneD, NekDouble>(1);
// 	    }

// 	  tmpcoord[0][0] = startcoordx;
// 	  tmpcoord[1][0] = startcoordy;
// 	  //	  tmpcoord[2][0] = startcoordz;
// 	  optima[0] = startcoordx;
// 	  optima[1] = startcoordy;
// 	  //optima[2] = startcoordz;
// 	  Vtmp = E->PhysEvaluateBasis(tmpcoord, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        

// 	}
//       demo.pq(utemp, tmpcoord, Vtmp, pqvalcoords, wsp1);
//       // cout<<"\n fvals: "<<wsp1.size()<<"\n";;
//       // for(int kk = 0 ; kk<wsp1.size(); kk++)
//       //     cout<<" "<<wsp1[kk]<<" ";
//       // cout<<"\n";

//       if (pqvalcoords[0] < pqval)
//         {
// 	  xastarr = optima;
// 	  pqval = pqvalcoords[0];
//         }
        
//       //cout<<"\n at counter="<<counter<<" min val="<<pqval<<" xast ="<<xastarr[0]<<" "<<xastarr[1];//<<" "<<xastarr[2]<<
//       //" ";
//       // If minimum is non-negative, we're done
//       if (pqval >= -tols.at(0))
//         {
// 	  break;
//         }
        
//       Array<OneD, NekDouble>  Vastsq (N1);
//       Array<OneD, NekDouble>  Vast (N1);
//       Array<OneD, Array<OneD, NekDouble> > xastarrofarr(2);
//       xastarrofarr[0] = Array<OneD, NekDouble>(1, xastarr[0]);
//       xastarrofarr[1] = Array<OneD, NekDouble>(1, xastarr[1]);
//       //      xastarrofarr[2] = Array<OneD, NekDouble>(1, xastarr[2]);
//       Array<OneD, NekDouble> tmp;
//       NekDouble vastsqsum;
//       Vast = E->PhysEvaluateBasis(xastarrofarr, storage2d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
	   
//       // for( int i = 0; i < N1; i++)
//       //   {
//       // 	  tmp = E->PhysEvaluateBasis(xastarrofarr, storage2d, i);
//       // 	  Vast[i] = tmp[0];
//       // 	  Vastsq[i] = (Vast[i]*Vast[i]);
	  
//       //   }
//       	for( int i = 0; i < N1; i++)
// 	  {
// 	    // if 3d, use storage 3d
// 	    //tmp = E->PhysEvaluateBasis(xastarrofarr, storage2d, i);
// 	    Vastsq[i] = (Vast[i]*Vast[i]);
// 	  }


//         vastsqsum = Vmath::Vsum(N1, &Vastsq[0], 1);

//         Array<OneD, NekDouble>  qast(N1);

//         for(int i = 0; i<N1; i++)
//         {
//             qast[i] = ((1/sqrt(vastsqsum))*(Vast[i]));
//         }
//         Vmath::Smul(N1, pqval, &qast[0], 1, &qast[0], 1);
        
//         Vmath::Vsub(utemp.size(), &utemp[0], 1, &qast[0], 1, &qast[0], 1);
//         d.push_back(qast);
// 	//	t.Stop();
// 	//evalrootshold+= t.TimePerTest(1);  

//         counter = counter + 1;
    
    
//     }
//     roots1dtimehold = roots1dtimehold/(counter-1);
//     roots2dtimehold = roots2dtimehold/(counter-1);
//     roots3dtimehold = roots3dtimehold/(counter-1);
//     timeprojectedges = timeprojectedges/(counter);
//     timeprojectsurf = timeprojectsurf/(counter);    
//     itersGD2 = (avgiterGD)/(counter);
//     iterstaken = counter;
//     cout<<"sphere_rotation took "<<counter<<"iterations. avgiterGD = "<<avgiterGD;//
//     //cout<<"\n  startv = "<<startval<<" at startcoords("<<startcoordx<<","<<startcoordy<<" ," <<startcoordz<<")    GD iter stepszie ="<<demo.GetTol()<<" maxiters = "<<demo.GetMaxIter()<<" o = "<<E->GetBasis(0)->GetNumModes()<<" p = "<<E->GetBasis(0)->GetNumPoints()<<" eps = "<<demo.GetEps()<<"  avgiterGD = "<<itersGD2;//<<" total time taken post rootfinding per iter = "<< evalrootshold/(counter)<<" hex \n";;
//     //    cout<<"\nAvg times per iter (1drootfinder,2drootfinder,3drootfinder) = : "  << roots1dtimehold<<", "<<roots2dtimehold<<", "<<roots3dtimehold ;//<< timeprojectedges<<", "<<timeprojectsurf<<",
//     uhats = d.back();
// }


// NekDouble Shape_sol(NekDouble x, NekDouble y, NekDouble z, vector<int> order,
//                     vector<BasisType> btype, ShapeType stype, bool diff)
// {
//   map<ShapeType, function<int(int, const vector<int> &)>> shapeConstraint2;
//   shapeConstraint2[ePoint] =
//     [](int,   const vector<int> &     ) { return 1; };
//   shapeConstraint2[eSegment] =
//     [](int,   const vector<int> &     ) { return 1; };
//   shapeConstraint2[eTriangle] =
//     [](int k, const vector<int> &order) { return order[1] - k; };
//   shapeConstraint2[eQuadrilateral] =
//     [](int,   const vector<int> &order) { return order[1]; };
//   shapeConstraint2[eTetrahedron] =
//     [](int k, const vector<int> &order) { return order[1] - k; };
//   shapeConstraint2[ePyramid] =
//     [](int k, const vector<int> &order) { return order[1] - k; };
//   shapeConstraint2[ePrism] =
//     [](int,   const vector<int> &order) { return order[1]; };
//   shapeConstraint2[eHexahedron] =
//     [](int,   const vector<int> &order) { return order[1]; };

//   map<ShapeType, function<int(int, int, const vector<int> &order)>>
//     shapeConstraint3;
//   shapeConstraint3[ePoint] =
//     [](int,   int,   const vector<int> &     ) { return 1; };
//   shapeConstraint3[eSegment] =
//     [](int,   int,   const vector<int> &     ) { return 1; };
//   shapeConstraint3[eTriangle] =
//     [](int,   int,   const vector<int> &     ) { return 1; };
//   shapeConstraint3[eQuadrilateral] =
//     [](int,   int,   const vector<int> &     ) { return 1; };
//   shapeConstraint3[eTetrahedron] =
//     [](int k, int l, const vector<int> &order) { return order[2] - k - l; };
//   shapeConstraint3[ePyramid] =
//     [](int k, int l, const vector<int> &order) { return order[2] - k - l; };
//   shapeConstraint3[ePrism] =
//     [](int k, int,   const vector<int> &order) { return order[2] - k; };
//   shapeConstraint3[eHexahedron] =
//     [](int,   int,   const vector<int> &order) { return order[2]; };

//   NekDouble sol = 0.0;
//   if (!diff)
//     {
//       if (btype[0] == eFourier && stype == eSegment)
//         {
// 	  for (int k = 0; k < order[0] / 2 - 1; ++k)
//             {
// 	      sol += sin(k * M_PI * x) + cos(k * M_PI * x);
//             }
//         }
//       else if (btype[0] == eFourierSingleMode && stype == eSegment)
//         {
// 	  sol += 0.25 * sin(M_PI * x) + 0.25 * cos(M_PI * x);
//         }
//       else if (btype[0] == eFourier && stype == eQuadrilateral)
//         {
// 	  if (btype[1] == eFourier)
//             {
// 	      for (int k = 0; k < order[0] / 2; ++k)
//                 {
// 		  for (int l = 0; l < order[1] / 2; ++l)
//                     {
// 		      sol += sin(k * M_PI * x) * sin(l * M_PI * y) +
// 			sin(k * M_PI * x) * cos(l * M_PI * y) +
// 			cos(k * M_PI * x) * sin(l * M_PI * y) +
// 			cos(k * M_PI * x) * cos(l * M_PI * y);
//                     }
//                 }
//             }
// 	  else if (btype[1] == eFourierSingleMode)
//             {
// 	      for (int k = 0; k < order[0] / 2; ++k)
//                 {
// 		  sol += sin(k * M_PI * x) * sin(M_PI * y) +
// 		    sin(k * M_PI * x) * cos(M_PI * y) +
// 		    cos(k * M_PI * x) * sin(M_PI * y) +
// 		    cos(k * M_PI * x) * cos(M_PI * y);
//                 }
//             }
// 	  else
//             {
// 	      for (int k = 0; k < order[0] / 2; ++k)
//                 {
// 		  for (int l = 0; l < order[1]; ++l)
//                     {
// 		      sol += sin(k * M_PI * x) * pow_loc(y, l) +
// 			cos(k * M_PI * x) * pow_loc(y, l) ;
//                     }
//                 }
//             }
//         }
//       else if (btype[0] == eFourierSingleMode && stype == eQuadrilateral)
//         {
// 	  if (btype[1] == eFourier)
//             {
// 	      for (int l = 0; l < order[1] / 2; ++l)
//                 {
// 		  sol += sin(M_PI * x) * sin(l * M_PI * y) +
// 		    sin(M_PI * x) * cos(l * M_PI * y) +
// 		    cos(M_PI * x) * sin(l * M_PI * y) +
// 		    cos(M_PI * x) * cos(l * M_PI * y);
//                 }

//             }
// 	  else if (btype[1] == eFourierSingleMode)
//             {
// 	      sol += sin(M_PI * x) * sin(M_PI * y) +
// 		sin(M_PI * x) * cos(M_PI * y) +
// 		cos(M_PI * x) * sin(M_PI * y) +
// 		cos(M_PI * x) * cos(M_PI * y);
//             }
// 	  else
//             {
// 	      for (int l = 0; l < order[1]; ++l)
//                 {
// 		  sol += sin(M_PI * x) * pow_loc(y, l) +
// 		    cos(M_PI * x) * pow_loc(y, l);
//                 }
//             }
//         }
//       else if (btype[1] == eFourier && stype == eQuadrilateral)
//         {
// 	  for (int k = 0; k < order[0]; ++k)
//             {
// 	      for (int l = 0; l < order[1] / 2; ++l)
//                 {
// 		  sol += sin(l * M_PI * y) * pow_loc(x, k) +
// 		    cos(l * M_PI * y) * pow_loc(x, k);
//                 }
//             }
//         }
//       else if (btype[1] == eFourierSingleMode && stype == eQuadrilateral)
//         {
// 	  for (int k = 0; k < order[0]; ++k)
//             {
// 	      sol += sin(M_PI * y) * pow_loc(x, k) +
// 		cos(M_PI * y) * pow_loc(x, k);
//             }
//         }
//       else
//         {
// 	  for (int k = 0;
// 	       k < order[0]; ++k) //ShapeConstraint 1 is always < order1
//             {
// 	      for (int l = 0; l < shapeConstraint2[stype](k, order); ++l)
//                 {
// 		  for (int m = 0;
// 		       m < shapeConstraint3[stype](k, l, order); ++m)
//                     {
// 		      sol += pow_loc(x, k) * pow_loc(y, l) * pow_loc(z, m);
//                     }
//                 }
//             }
//         }
//     }
//   else if (diff)
//     {
//       if (btype[0] == eFourier && stype == eSegment)
//         {
// 	  for (int k = 0; k < order[0] / 2 - 1; ++k)
//             {
// 	      sol += k * M_PI * (cos(k * M_PI * z) - sin(k * M_PI * z));
//             }
//         }
//       else if (btype[0] != eFourier && btype[1] == eFourier &&
// 	       stype == eQuadrilateral)
//         {
// 	  for (int k = 0; k < order[0]; ++k)
//             {
// 	      for (int l = 0; l < order[1] / 2; ++l)
//                 {
// 		  sol += k * pow_loc(x, k - 1) * sin(M_PI * l * y)
// 		    + M_PI * l * pow_loc(x, k) * cos(M_PI * l * y) +
// 		    +k * pow_loc(x, k - 1) * cos(M_PI * l * y)
// 		    - M_PI * l * pow_loc(x, k) * sin(M_PI * l * y);
//                 }
//             }
//         }
//       else if (btype[0] == eFourier && btype[1] != eFourier &&
// 	       stype == eQuadrilateral)
//         {
// 	  for (int k = 0; k < order[0] / 2; ++k)
//             {
// 	      for (int l = 0; l < order[1]; ++l)
//                 {
// 		  sol += M_PI * k * cos(M_PI * k * x) * pow_loc(y, l)
// 		    + l * sin(M_PI * k * x) * pow_loc(y, l - 1) +
// 		    -M_PI * k * sin(M_PI * k * x) * pow_loc(y, l)
// 		    + l * sin(M_PI * k * x) * pow_loc(y, l - 1);
//                 }
//             }
//         }
//       else if (btype[0] == eFourier && btype[1] == eFourier &&
// 	       stype == eQuadrilateral)
//         {
// 	  for (int k = 0; k < order[0] / 2; ++k)
//             {
// 	      for (int l = 0; l < order[1] / 2; ++l)
//                 {
// 		  sol += M_PI * k * cos(M_PI * k * x) * sin(M_PI * l * y)
// 		    + M_PI * l * sin(M_PI * k * x) * cos(M_PI * l * y)
// 		    + M_PI * k * cos(M_PI * k * x) * cos(M_PI * l * y)
// 		    - M_PI * l * sin(M_PI * k * x) * sin(M_PI * l * y)
// 		    - M_PI * k * sin(M_PI * k * x) * sin(M_PI * l * y)
// 		    + M_PI * l * cos(M_PI * k * x) * cos(M_PI * l * y)
// 		    - M_PI * k * sin(M_PI * k * x) * cos(M_PI * l * y)
// 		    - M_PI * l * cos(M_PI * k * x) * sin(M_PI * l * y);
//                 }
//             }
//         }
//       else
//         {
// 	  NekDouble a;
// 	  for (int k = 0;
// 	       k < order[0]; ++k) //ShapeConstraint 1 is always < order1
//             {
// 	      for (int l = 0; l < shapeConstraint2[stype](k, order); ++l)
//                 {
// 		  for (int m = 0;
// 		       m < shapeConstraint3[stype](k, l, order); ++m)
//                     {
// 		      a = k * pow_loc(x, k - 1) * pow_loc(y, l) *
// 			pow_loc(z, m);
// 		      sol += a;
// 		      a = l * pow_loc(x, k) * pow_loc(y, l - 1) *
// 			pow_loc(z, m);
// 		      sol += a;
// 		      a = m * pow_loc(x, k) * pow_loc(y, l) *
// 			pow_loc(z, m - 1);
// 		      sol += a;
//                     }
//                 }
//             }
//         }
//     }

//   return sol;
// }
