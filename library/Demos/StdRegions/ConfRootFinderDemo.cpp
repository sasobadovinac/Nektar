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
// Description: Demo to test functionality of confederate matrix rootfinder
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "StdDemoSupport.hpp"
#include <LibUtilities/BasicUtils/Timer.h>

namespace po = boost::program_options;

Array<OneD,NekDouble>  callConfRF(Array<OneD,  NekDouble> uhats, Array<OneD,  NekDouble> phys);


//string funcdef;
Array<OneD, Array<OneD, NekDouble> > storage;

// Colleague matrix
int numedges, numsurfaces;
//Array<OneD, Array<OneD, NekDouble> > C;
StdExpansion *E;

int dimension ;

DemoSupport demo;
LibUtilities::ShapeType stype;
Array<OneD, Array<OneD, NekDouble> > edgeptsin;

int main(int argc, char *argv[])
{
  demo.GetOptions().add_options()("optm,z", "positivity preserving optimizer");
  demo.ParseArguments(argc, argv);

  po::variables_map vm = demo.GetVariableMap();

  //only for 1D
  E = demo.CreateStdExpansion();

  stype = E->DetShapeType();

  if (E == nullptr)
    {
      return 1;
    }
  dimension = E->GetShapeDimension();
  if(dimension != 1 )
    {
      cout<<"\n dimension should be 1 for this demo!";
      exit(0);
    }
  else
    {
      storage = E->GetPhysEvaluateStorage();
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

    default:
      cout<<"invalid dimension of problem!\n";
      break;
    }
  demo.funcdef = "x^2 - 0.09";

  //get solution array
  for (int i = 0; i < totPoints; ++i)
    {
      if(dimension == 1)
	{
	  sol[i] = pow(x[i]-0.1,2) - 0.09;//2*pow(x[i],3) -0.5*pow(x[i]-0.1, 2)-x[i];
	  //;//pow((x[i]-0.7),3)+pow((x[i]-0.2),2)-x[i] +0.4;//;//;
	  //pow(x[i]-0.1,2) - 0.09;
	  //floor(x[i]<=0 && y[i] <=0)+0.1;
	  //-1*(sin((x[i]-0.1)+0.5*M_PI))*cos((y[i]-0.2))+0.95; 
	}
      else
	{
	  cout<<"\n This demo only tests 1-D root-fnding\n";
	  exit(0);
	}
    }
  Array<OneD, NekDouble> phys(totPoints);
  Array<OneD, NekDouble> coeffs((unsigned) E->GetNcoeffs());
  Array<OneD, NekDouble> V(coeffs.size()), Vx(coeffs.size());
  //Project onto expansion
  E->FwdTrans(sol, coeffs);

  E->BwdTrans(coeffs, phys);
  
  Array<OneD, Array<OneD, NekDouble> > dummy, surf(1);
  surf[0] = Array<OneD, NekDouble>(coeffs.size());
  Vmath::Vcopy(coeffs.size(), &coeffs[0], 1, &surf[0][0], 1); 

  Array<OneD, NekDouble> optima = callConfRF(coeffs, phys);

}


Array<OneD,  NekDouble> callConfRF(Array<OneD,  NekDouble> coeff, Array<OneD,  NekDouble> phys)
{
  boost::ignore_unused(phys);
  int dimension = E->GetShapeDimension();
  Array<OneD, NekDouble> coords(dimension), tmp, evalbasis(coeff.size());
  Array<OneD, Array<OneD, NekDouble > >  rethold(dimension), holdstorage;
  Array<OneD, NekDouble> copycoeff(coeff.size());
  Vmath:: Vcopy(coeff.size(), coeff, 1, copycoeff, 1);
  Array<OneD, BasisType> btypearr(1);
  btypearr[0] = LibUtilities::eOrtho_A;
  NekDouble dummyv;
  if(stype == LibUtilities::ePoint)
    {
      cout<<"\n This demo works only for dim = 1";

      exit(0);
    }
   if(E->GetBasis(0)->GetBasisType() != LibUtilities::eOrtho_A)
     {
       demo.OrthoNormalize(E, copycoeff, btypearr, 0 ); //non-ortho->ortho
     }
   rethold = demo.find_roots(copycoeff, nullptr, NullNekDoubleArrayOfArray, dummyv , 0, 0, 0, 0);
   
   cout<<"\nconf matrix returns roots =  ";
  
   if(rethold[0].size() > 0)
     {
       for(int k = 0; k < rethold[0].size(); k++)
	 {
	   cout<<" \t"<<rethold[0][k]<<" ";
	 }
     }
   
  //ignore -1 and 1
  tmp = Array<OneD, NekDouble>(rethold[0].size()-2);
  cout<<"\nvals at these roots :\t";
  Array<OneD, Array<OneD, NekDouble> > currx(dimension);
  currx[0] = Array<OneD, NekDouble>(1);
  //StdExpansion *Ehold;
  //  holdstorage = storage;
  // if(E->GetBasis(0)->GetBasisType() != LibUtilities::eOrtho_A)
  //   {
  //     demo.OrthoNormalize(E, copycoeff, 1 );
  //     //  Ehold = Etmp;
  //     //holdstorage = Etmp->GetPhysEvaluateStorage();
  //     //       Etmp->BwdTrans(copycoeff, phys);
  //     //E->FwdTrans(phys, coeff);
  //   }
  // else
  //   {
  //     // Ehold = E;
  //     //holdstorage = storage;
  //   }
  // Vmath::Vcopy(coeff.size(),copycoeff, 1, coeff, 1);
  for(int k = 0; k < tmp.size(); k++)
    {
      currx[0][0] = rethold[0][k];
      evalbasis = E->PhysEvaluateBasis(currx,storage, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
      Vmath::Vmul(evalbasis.size(), evalbasis,  1, coeff, 1, evalbasis, 1 );
      tmp[k] = Vmath::Vsum(evalbasis.size(), evalbasis, 1);
      cout<<"\t "<<tmp[k];

    }
  Array<OneD, NekDouble> s1(E->GetTotPoints(),0.0), p1(E->GetTotPoints(),0.0);	    
  NekDouble eLinf, eL2;
  for(int k = 0; k < tmp.size(); k++)
    {
      p1[k] = tmp[k]; 
    }
  eLinf = E->Linf(s1,p1);
  eL2 = E->L2(s1,p1);
  
  cout << "\nL infinity error : " << scientific << eLinf << endl;
  cout << "L 2 error        : " << scientific << eL2 << endl<<" *****\n";
	 
  return coords;	
}


