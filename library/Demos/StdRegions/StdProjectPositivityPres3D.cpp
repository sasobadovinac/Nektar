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
#include <fstream>
#include <time.h>
#include <LibUtilities/BasicUtils/Timer.h>

#define LOGNAME_FORMAT "%Y%m%d_%H%M%S"
#define LOGNAME_SIZE 20

namespace po = boost::program_options;

NekDouble Shape_sol(NekDouble x, NekDouble y, NekDouble z, vector<int> order,
                    vector<BasisType> btype, ShapeType stype, bool diff);

//Modification to deal with exact solution for diff. Return 1 if integer < 0.
static double pow_loc(const double val, const int i)
{
    return (i < 0) ? 1.0 : pow(val, i);
}



vector< vector<NekDouble> > arrtovec(Array<OneD, Array<OneD, NekDouble> > arr);
Array<OneD, Array<OneD, NekDouble> > vectoarr(vector<vector< NekDouble> > vec);

//declare Do_optimize
void Do_optimize(Array<OneD, NekDouble> &uhats);

void surfuhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz );

// declare caller routine to find_roots
// flag = 0 -> opt_needed calls
// flag = 1 -> sphere_rot calls
Array<OneD, Array<OneD,  NekDouble> > call_find_roots(Array<OneD, NekDouble> &uhatsall, NekDouble &avgiterGD, Array< OneD, Array<OneD, NekDouble> >&uhatsedges, Array< OneD, Array<OneD, NekDouble> >&surfaceuhats, NekDouble &roots1dtime,  NekDouble &roots2dtime , NekDouble &roots3dtime, int flag = 0 );


//declare Opt_needed
int Opt_needed(Array<OneD, NekDouble> uhats, int flag = 0);


// for 2D elements to get uhats at edges
// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble> uhats,Array<OneD, Array<OneD, NekDouble> >&ret , int d = 0);


// uhatpqd stuff: 
// called by project_edges if d = 1 is passed in project_edges
void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret, int modes,   Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1 = NullNekDouble1DArray , Array<OneD, NekDouble> Vxyd2 = NullNekDouble1DArray);

void deruhatsedges(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble >&ret,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0 = NullNekDouble1DArray,  Array<OneD, NekDouble> Vxyd1 = NullNekDouble1DArray, Array<OneD, NekDouble> Vxyd2 = NullNekDouble1DArray);

void project_surfaces( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret);

string funcdef;
int iterstaken;


FILE *csvfile(int flag);
void WriteCSV( Array<OneD, NekDouble>  uhats, int flag = 0);

// Colleague matrix
int numedges, numsurfaces;
StdExpansion *E;
StdExpansion *E3seg;
StdExpansion *Equad;
Array<OneD, NekDouble> qZin;
Array<OneD, Array<OneD, NekDouble> > edgeptsin;

Array<OneD, NekDouble> qZinmid; //inner mid point grid for 2D rootfinding
Array<OneD, Array<OneD, Array<OneD, NekDouble> > > surfcoords;

// for hex edges root finding

// edge front left (x = -1) (y = -1)
Array<OneD, NekDouble> Vxm1ym1z   ;
Array<OneD, NekDouble> Vdyxm1ym1z   ;
Array<OneD, NekDouble> Vdxxm1ym1z   ;
Array<OneD, NekDouble> Vdzxm1ym1z   ;
                
//edge front right (x = 1) (y = -1)
Array<OneD, NekDouble> Vx1ym1z   ;
Array<OneD, NekDouble> Vdyx1ym1z   ;
Array<OneD, NekDouble> Vdxx1ym1z   ;
Array<OneD, NekDouble> Vdzx1ym1z   ;
               
//edge front top (y = -1) (z = 1)
Array<OneD, NekDouble> Vym1xz1   ;
Array<OneD, NekDouble> Vdyym1xz1   ;
Array<OneD, NekDouble> Vdxym1xz1   ;
Array<OneD, NekDouble> Vdzym1xz1   ;
               
//edge front bot (y = -1) (z = -1)
Array<OneD, NekDouble> Vym1xzm1   ;
Array<OneD, NekDouble> Vdyym1xzm1   ;
Array<OneD, NekDouble> Vdxym1xzm1   ;
Array<OneD, NekDouble> Vdzym1xzm1   ;


// edge back left (y = 1), (x = -1)
Array<OneD, NekDouble> Vxm1y1z   ;
Array<OneD, NekDouble> Vdyxm1y1z   ;
Array<OneD, NekDouble> Vdxxm1y1z   ;
Array<OneD, NekDouble> Vdzxm1y1z   ;
                
//edge back right (x = 1), (y = 1))
Array<OneD, NekDouble> Vx1y1z   ;
Array<OneD, NekDouble> Vdyx1y1z   ;
Array<OneD, NekDouble> Vdxx1y1z   ;
Array<OneD, NekDouble> Vdzx1y1z   ;
               
//edge back top ( y = 1) (z = 1)
Array<OneD, NekDouble> Vy1xz1   ;
Array<OneD, NekDouble> Vdyy1xz1   ;
Array<OneD, NekDouble> Vdxy1xz1   ;
Array<OneD, NekDouble> Vdzy1xz1   ;
               
//edge back bot (y = 1) (z = -1))
Array<OneD, NekDouble> Vy1xzm1   ;
Array<OneD, NekDouble> Vdyy1xzm1   ;
Array<OneD, NekDouble> Vdxy1xzm1   ;
Array<OneD, NekDouble> Vdzy1xzm1   ;

// edge left bot (z = -1), (x = -1)
Array<OneD, NekDouble> Vxm1yzm1   ;
Array<OneD, NekDouble> Vdyxm1yzm1   ;
Array<OneD, NekDouble> Vdxxm1yzm1   ;
Array<OneD, NekDouble> Vdzxm1yzm1   ;
                
//edge left top (x = -1), (z = 1))
Array<OneD, NekDouble> Vxm1yz1   ;
Array<OneD, NekDouble> Vdyxm1yz1   ;
Array<OneD, NekDouble> Vdxxm1yz1   ;
Array<OneD, NekDouble> Vdzxm1yz1   ;
               
//edge right bot ( z = -1) (x = 1)
Array<OneD, NekDouble> Vx1yzm1   ;
Array<OneD, NekDouble> Vdyx1yzm1   ;
Array<OneD, NekDouble> Vdxx1yzm1   ;
Array<OneD, NekDouble> Vdzx1yzm1   ;
               
//edge right top (z  1) (x  1))
Array<OneD, NekDouble>Vx1yz1 ;
Array<OneD, NekDouble>Vdxx1yz1 ;
Array<OneD, NekDouble>Vdyx1yz1 ;
Array<OneD, NekDouble> Vdzx1yz1;

//surface bot z = -1
Array<OneD, NekDouble>Vxyzm1 ;
// Array<OneD, NekDouble>Vdyxyzm1 ;
// Array<OneD, NekDouble>Vdxxyzm1 ;
// Array<OneD, NekDouble> Vdzxyzm1 ;

//surface right x = 1
Array<OneD, NekDouble>Vx1yz ;
// Array<OneD, NekDouble>Vdyx1yz ;
// Array<OneD, NekDouble>Vdxx1yz ;
// Array<OneD, NekDouble> Vdzx1yz ;


//surface top z = 1
Array<OneD, NekDouble>Vxyz1 ;
// Array<OneD, NekDouble>Vdyxyz1 ;
// Array<OneD, NekDouble>Vdxxyz1 ;
// Array<OneD, NekDouble> Vdzxyz1 ;
 

//surface left x = -1
Array<OneD, NekDouble>Vxm1yz ;
// Array<OneD, NekDouble>Vdyxm1yz ;
// Array<OneD, NekDouble>Vdxxm1yz ;
// Array<OneD, NekDouble> Vdzxm1yz ;


//surface front y = -1
Array<OneD, NekDouble>Vxym1z ;
// Array<OneD, NekDouble>Vdxxym1z ;
// Array<OneD, NekDouble>Vdyxym1z ;
// Array<OneD, NekDouble> Vdzxym1z ;
  

//surface back y = 1
Array<OneD, NekDouble>Vxy1z ;
// Array<OneD, NekDouble>Vdxxy1z ;
// Array<OneD, NekDouble> Vdyxy1z;
// Array<OneD, NekDouble> Vdzxy1z ;
  

int dimension ;
NekDouble  startval, startcoordx, startcoordy, startcoordz, itersGD1, itersGD2, itersGD3;

DemoSupportNew demo;

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
    for (int i = 0; i < totPoints; ++i)
    {
        if(dimension ==3) //hex
        {
	  sol[i] = pow((x[i]-0.5),2)+ pow((y[i]+0.5),2) + pow((z[i]-0.2),2) + demo.GetEps();//; //sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-1;// ;//sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2;//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9;//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2 ;//sin(x[i])*sin(y[i])*sin(z[i]);//+0.1;
	  funcdef ="pow((x[i]-0.5),2)+ pow((y[i]+0.5),2) + pow(z[i]-0.2),2) + demo.eps";//sin(M_PI*x[i])*cos(M_PI*y[i])*sin(M_PI*z[i]) +0.5; //sin(2*M_PI*x[i])*sin(2*//" sin(M_PI*x[i])*cos(M_PI*y[i])*sin(M_PI*z[i]) +0.5";//sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-1";////sin(x[i])(y[i])*sin(z[i]) -1";//sin(x[i])*sin(y[i])*sin(z[i]);";// sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i]) +0.4";//"sin(x[i])*sin(y[i])*sin(z[i])";//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.5";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2 ";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"=sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);;";//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1;";//"sin(x[i])*sin(y[i])*sin(z[i];";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2 ";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.5";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);";// "sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])  -2";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";// -  exp(-1) + 1.37";

        }
    }
    Array<OneD, NekDouble> phys(totPoints);
    Array<OneD, NekDouble> coeffs((unsigned) E->GetNcoeffs());
    
    //Project onto expansion
    E->FwdTrans(sol, coeffs);

    //Backward transform solution to get projected values
    E->BwdTrans(coeffs, phys);
    //cout<<"\n sol:";
    // for(int k = 0; k < sol.size(); k++)
    //   cout<<sol[k]<< " ";
    // cout<<"\n phys:";
    // for(int k = 0; k < phys.size(); k++)
    //   cout<<phys[k]<< " ";

    LibUtilities::Timer     timer;
    //    NekDouble elapsed       = 0.0;

    // check for -ve values and apply opt if necessary
    if (vm.count("optm"))
    {
      NekDouble setup, timeneg, timeverify;
      NekDouble elapref = 0;
      Array<OneD, NekDouble> physder(phys.size()), physder0(phys.size()), uhatsnew(coeffs.size(), 1.0);
      for(int k = 0; k < demo.GetAvgNum(); k++)
	{
      
	  timer.Start();
	  //Backward transform solution to get projected values
	  E->BwdTrans(uhatsnew, physder);
	  
	  E->PhysDeriv(physder, physder0);
	  //Project onto expansion
	  E->FwdTrans(physder0, uhatsnew); 
	  
	  timer.Stop();
	  elapref += timer.TimePerTest(1);
	}

      elapref = elapref/demo.GetAvgNum();
      cout<<"\n first derivative uhats from f(x) uhats took "<<(elapref)<<" s\n\n";
	  

      timer.Start();
        int dimension = E->GetShapeDimension(); 
        
        vector<LibUtilities::BasisKey> tmpBasisKey = demo.GetBasisKey();       
        LibUtilities::PointsKey pkeycheb(E->GetBasis(0)->GetNumPoints(), pointsTypeCheb);
        
        LibUtilities::BasisKey bkeycheb(LibUtilities::eChebyshev, 3*(E->GetBasis(0)->GetNumModes()-1),  pkeycheb);
        
        E3seg = new StdSegExp(bkeycheb);
        if (E3seg == nullptr)
        {
            return 1;
        }
        LibUtilities::BasisKey bkeycheb2 (LibUtilities::eChebyshev, (E->GetBasis(0)->GetNumModes()), pkeycheb);

        Equad = new StdQuadExp(bkeycheb2, bkeycheb2);
        if (Equad == nullptr)
        {
            return 1;
        }
        // populate edge root finding matrices for quadrilateral
	//        LibUtilities::PointsKey  quadPointsKeyin  (order[0], pointsTypeCheb);
        // for gradient descent
        qZin = Equad->GetBasis(0)->GetZ();
        int sz = pow(qZin.size(), 2);;
        int szmid = pow(qZin.size()-1,2);
        qZinmid = Array<OneD,NekDouble>(qZin.size() - 1);
	demo.testcoord2dquad = Array<OneD, Array<OneD, NekDouble> > (2); 

        int ctsz2d = 0;
        
        for(int ii = 0; ii < 2; ++ii)
        {
            demo.testcoord2dquad[ii] = Array<OneD, NekDouble>((sz+szmid));
        
        }

        for(int ii = 0; ii < qZin.size(); ++ii)
        {
                
            for(int jj = 0; jj < qZin.size(); ++jj)
            {
                demo.testcoord2dquad[0][ctsz2d] = qZin[ii];
                demo.testcoord2dquad[1][ctsz2d] = qZin[jj];
                ctsz2d++;
            }
                
        }        


        // add midpt grid
        for(int ii = 0; ii < qZinmid.size(); ++ii)
        {
            for(int jj = 0; jj < qZinmid.size(); ++jj)
            {
	      demo.testcoord2dquad[0][ctsz2d] =  (qZin[ii] +  qZin[ii+1])/2.0;//qZinmid[ii];
	      demo.testcoord2dquad[1][ctsz2d] = (qZin[jj] +  qZin[jj+1])/2.0;//qZinmid[jj];

                ctsz2d++;
                    
            }
        
        }  
            
        sz = pow((qZin.size()),3)+pow(qZinmid.size(),3);;

        Array<OneD, Array<OneD, NekDouble> >testcoord3dtemp  (dimension); 
	Array<OneD, Array<OneD, NekDouble> > midcoords3d(dimension);
	
	demo.testcoord3d = Array<OneD, Array<OneD, NekDouble> > (dimension); 
	int szmid3d = pow(0+qZinmid.size(),3);
	
        for(int ii = 0; ii < 3; ++ii)
        {
            testcoord3dtemp[ii] = Array<OneD, NekDouble>(sz);

            midcoords3d[ii] = Array<OneD, NekDouble>(szmid3d,3);
        }
        int ctsz = 0;
        for(int ii = 0; ii < qZin.size(); ++ii)
        {
                    
            for(int jj = 0; jj < qZin.size(); ++jj)
            {

                for(int kk = 0; kk < qZin.size(); ++kk)
                {
                    testcoord3dtemp[0][ctsz] = qZin[ii];
                    testcoord3dtemp[1][ctsz] = qZin[jj];
                    testcoord3dtemp[2][ctsz] = qZin[kk];
                    ctsz++;
                }
            }
        }
	int ct;
        for(int ii = 0; ii < qZinmid.size(); ++ii)
        {
                    
            for(int jj = 0; jj < qZinmid.size(); ++jj)
            {

                for(int kk = 0; kk < qZinmid.size(); ++kk)
                {

		  midcoords3d[0][ct] = (qZin[ii] +  qZin[ii+1])/2.0;//qZinmid[ii];
		  midcoords3d[1][ct] = (qZin[jj] +  qZin[jj+1])/2.0;//qZinmid[jj];
		  midcoords3d[2][ct] = (qZin[kk] +  qZin[kk+1])/2.0;//qZinmid[kk];
		  ct++;
		  
		  testcoord3dtemp[0][ctsz] = qZinmid[ii];
                    testcoord3dtemp[1][ctsz] = qZinmid[jj];
                    testcoord3dtemp[2][ctsz] = qZinmid[kk];
                    ctsz++;
                }
            }
        }

        for(int ii = 0; ii < dimension; ++ii)
        {
            demo.testcoord3d[ii] = Array<OneD, NekDouble>(ctsz);
            Vmath::Vcopy(ctsz, &testcoord3dtemp[ii][0], 1, &demo.testcoord3d[ii][0], 1);

        }
        demo.interioreval2dquad = Array<OneD, NekDouble >(Equad->GetNcoeffs()*demo.testcoord2dquad[0].size());
        Equad->PhysEvalBasisGrad(demo.testcoord2dquad, demo.interioreval2dquad , NullNekDouble1DArray , NullNekDouble1DArray, NullNekDouble1DArray);
        demo.interioreval3d = Array<OneD, NekDouble >(E->GetNcoeffs()*demo.testcoord3d[0].size());
	//cout<<"\n eval poly lattice now E totpts="<<E->GetTotPoints()<<" demo.testcoord3d[0] sz = "<<demo.testcoord3d[0].size()<<"demo. interioreval3d sz = "<<E->GetNcoeffs()*demo.testcoord3d[0].size()<<":\n";
	Timer t1;
	
	int nqmid = ct;
	
	int nq = E->GetTotPoints();
	Array<OneD, NekDouble> tmp(nq);
	Array<OneD, NekDouble> tempevalmid( ct* E->GetNcoeffs());
	NekDouble polyavgevaltime = 0.0;
	for(int i = 0; i < E->GetNcoeffs(); i++)
	  {
	    Vmath::Vcopy(nq, &E->m_physevalall[0][i*nq], 1, & demo.interioreval3d[i*ctsz], 1);       
	  }
	cout<<"\n midcoords3d sz="<<midcoords3d.size()<< ","<< midcoords3d[0].size();
	for(int k = 0; k < 1; k++)
	  {
	    t1.Start();
	    E->PhysEvalBasisGrad(midcoords3d, tempevalmid, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
	  
	    t1.Stop();
	    polyavgevaltime += 	  t1.TimePerTest(1);
	  }
	for(int i = 0; i < E->GetNcoeffs(); i++)
	  {
	    Vmath::Vcopy(nqmid, &tempevalmid[i*(ct)], 1, & demo.interioreval3d[i*(ctsz)+nq], 1);
	  }

	cout<<"\n time eval poly lattice: "<<polyavgevaltime/1.0<<" \n";

	//exit(0);

        if(numedges == 12) //hex
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

            
            E->PhysEvalBasisGrad(edgeptsin, Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);  
                
            //edge front right (x = 1) (y = -1)
            Vx1ym1z = Array<OneD, NekDouble>(totszedges1d);
            Vdxx1ym1z = Array<OneD, NekDouble>(totszedges1d);
            Vdyx1ym1z = Array<OneD, NekDouble>(totszedges1d);
            Vdzx1ym1z = Array<OneD, NekDouble>(totszedges1d);
            edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);       
            E->PhysEvalBasisGrad(edgeptsin, Vx1ym1z, Vdxx1ym1z, Vdyx1ym1z, Vdzx1ym1z);  


            //edge front top (y = -1) (z = 1)
            Vym1xz1 = Array<OneD, NekDouble>(totszedges1d);
            Vdxym1xz1 = Array<OneD, NekDouble>(totszedges1d);
            Vdyym1xz1 = Array<OneD, NekDouble>(totszedges1d);
            Vdzym1xz1 = Array<OneD, NekDouble>(totszedges1d);
            edgeptsin[0] = edgexyztemp; 
            edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
            E->PhysEvalBasisGrad(edgeptsin,Vym1xz1, Vdxym1xz1, Vdyym1xz1, Vdzym1xz1);  

                
            //edge front bot (y = -1) (z = -1)
            Vym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
            Vdxym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
            Vdyym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
            Vdzym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
            edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
            E->PhysEvalBasisGrad( edgeptsin, Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);  

            // edge back left (y = 1), (x = -1)
            Vxm1y1z = Array<OneD, NekDouble>(totszedges1d);
            Vdxxm1y1z = Array<OneD, NekDouble>(totszedges1d);
            Vdyxm1y1z = Array<OneD, NekDouble>(totszedges1d);
            Vdzxm1y1z = Array<OneD, NekDouble>(totszedges1d);
            edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
            edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);      
            edgeptsin[2] =    edgexyztemp;


            E->PhysEvalBasisGrad(edgeptsin,Vxm1y1z, Vdxxm1y1z, Vdyxm1y1z, Vdzxm1y1z);  

            //edge back right (x = 1), (y = 1))
            Vx1y1z = Array<OneD, NekDouble>(totszedges1d);
            Vdxx1y1z = Array<OneD, NekDouble>(totszedges1d);
            Vdyx1y1z = Array<OneD, NekDouble>(totszedges1d);
            Vdzx1y1z = Array<OneD, NekDouble>(totszedges1d);
            edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);                
         
            E->PhysEvalBasisGrad(edgeptsin,Vx1y1z, Vdxx1y1z, Vdyx1y1z, Vdzx1y1z);  

                
            //edge back top ( y = 1) (z = 1)
            Vy1xz1 = Array<OneD, NekDouble>(totszedges1d);
            Vdxy1xz1 = Array<OneD, NekDouble>(totszedges1d);
            Vdyy1xz1 = Array<OneD, NekDouble>(totszedges1d);
            Vdzy1xz1 = Array<OneD, NekDouble>(totszedges1d);
            edgeptsin[0] = edgexyztemp;
            edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
            E->PhysEvalBasisGrad(edgeptsin, Vy1xz1, Vdxy1xz1, Vdyy1xz1, Vdzy1xz1 );  

      
            //edge back bot (y = 1) (z = -1))
            Vy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
            Vdxy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
            Vdyy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
            Vdzy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
            edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
            
            E->PhysEvalBasisGrad(edgeptsin, Vy1xzm1, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);  

            // edge left bot (z = -1), (x = -1)
            Vxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
            Vdxxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
            Vdyxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
            Vdzxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
            edgeptsin[1] = edgexyztemp;
            edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

            E->PhysEvalBasisGrad(edgeptsin,Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);  
                
            //edge left top (x = -1), (z = 1))
            Vxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
            Vdxxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
            Vdyxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
            Vdzxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
                
            edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);           
            E->PhysEvalBasisGrad(edgeptsin,Vxm1yz1, Vdxxm1yz1, Vdyxm1yz1, Vdzxm1yz1);  

            //edge right bot ( z = -1) (x = 1)
            Vx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
            Vdxx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
            Vdyx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
            Vdzx1yzm1 = Array<OneD, NekDouble>(totszedges1d);

            edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
            edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);           
            E->PhysEvalBasisGrad(edgeptsin, Vx1yzm1, Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);   

            //edge right top (z  1) (x  1))
            Vx1yz1 = Array<OneD, NekDouble>(totszedges1d);
            Vdxx1yz1 = Array<OneD, NekDouble>(totszedges1d);
            Vdyx1yz1 = Array<OneD, NekDouble>(totszedges1d);
            Vdzx1yz1 = Array<OneD, NekDouble>(totszedges1d);
                
            edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);

            E->PhysEvalBasisGrad(edgeptsin,Vx1yz1, Vdxx1yz1, Vdyx1yz1, Vdzx1yz1);   
            surfcoords = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(numsurfaces);

              
            int totszsurf2d = (demo.testcoord2dquad[0].size())*E->GetNcoeffs();
            
            Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);
            for(int k = 0; k < dimension-1; k++)
            {
                surfptsin[k] = Array<OneD, NekDouble>(demo.testcoord2dquad[k]);
                surfptsintemp[k] = Array<OneD, NekDouble>(demo.testcoord2dquad[k]);
            }
            
            surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.testcoord2dquad[0].size()); 
            surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.testcoord2dquad[0].size()); 
            
            //surface bot z = -1
            Vxyzm1 =  Array<OneD, NekDouble>(totszsurf2d);
            surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
            surfptsintemp[0] = surfptsin[0];
            surfptsintemp[1] = surfptsin[1];
            surfcoords[0] = surfptsintemp;                    
            E->PhysEvalBasisGrad(surfptsintemp,Vxyzm1, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

            //surface right x = 1
            Vx1yz =  Array<OneD, NekDouble>(totszsurf2d);
            surfptsintemp[0] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
            surfptsintemp[1] =  surfptsin[0];
            surfptsintemp[2] =  surfptsin[1];
            surfcoords[1] = surfptsintemp;                        
                
            E->PhysEvalBasisGrad(surfptsintemp,Vx1yz,NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   


            //surface top z = 1
            Vxyz1 =  Array<OneD, NekDouble>(totszsurf2d);
            surfptsintemp[0] = surfptsin[0];
            surfptsintemp[1] = surfptsin[1];
            surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
            surfcoords[2] = surfptsintemp;                        
            E->PhysEvalBasisGrad(surfptsintemp, Vxyz1, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

            //surface left x = -1
            Vxm1yz =  Array<OneD, NekDouble>(totszsurf2d);
            surfptsintemp[0] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
            surfptsintemp[1] = surfptsin[0];
            surfptsintemp[2] = surfptsin[1];
            surfcoords[3] = surfptsintemp; 
            E->PhysEvalBasisGrad(surfptsintemp,Vxm1yz, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

            //surface front y = -1
            Vxym1z =  Array<OneD, NekDouble>(totszsurf2d);

            surfptsintemp[1] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
            surfptsintemp[0] = surfptsin[0];
            surfptsintemp[2] = surfptsin[1];
            surfcoords[4] = surfptsintemp; 
                
            E->PhysEvalBasisGrad(surfptsintemp,Vxym1z, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   
      
            //surface back y = 1
            Vxy1z =  Array<OneD, NekDouble>(totszsurf2d);

            surfptsintemp[1] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
            surfptsintemp[0] = surfptsin[0];
            surfptsintemp[2] = surfptsin[1];
            surfcoords[5] = surfptsintemp; 
                
            E->PhysEvalBasisGrad(surfptsintemp,Vxy1z, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        } 

	  timer.Stop();
	  setup = timer.TimePerTest(1);

	  cout<<"\n setup phase took "<<setup<<"s\n\n";

	  timer.Start();
	  Array<OneD, NekDouble> oldcoeffs;//(coeffs);
	  int opt_need = Opt_needed(coeffs);            
	  timer.Stop();
	  timeneg = timer.TimePerTest(1);  
	  cout<<"\n checking for -vity took "<<timeneg<<"s\n\n";

	  if(opt_need)
	  {
	    //elapsed += timer.TimePerTest(1);

            cout<<"\n need optimization";
            NekDouble elap = 0;
	    //            WriteCSV(coeffs);
	    //for(int k = 0; k < demo.GetAvgNum(); k++)
	    // {
		oldcoeffs = coeffs;
		timer.Start();
		Do_optimize(oldcoeffs);
		timer.Stop();
		//		elapsed += timer.TimePerTest(1);
		elap += timer.TimePerTest(1); 
		//}
	    coeffs = oldcoeffs;
	    //	    elap = elap/demo.GetAvgNum();
	    cout<<"\n optimizertook "<<elap<<"s\n\n";
	      
	    cout<<"\n doopt done\n verifying...\n";
	    
	    timer.Start();
	    opt_need = Opt_needed(coeffs, 1 );
	    timer.Stop();
	    timeverify = timer.TimePerTest(1);
	    cout<<"\n Verification took "<<timeverify<<"s\n\n";
	    
            if(opt_need)
            {
	      
	      cout<<"func="<<funcdef;
	      
	      cout<<"\n fail\n";
	      cout<<" "<<" "<<startval<<", ("<<startcoordx<<" "<<startcoordy<<" "<<startcoordz<<"), "<<setup<<", fail, "<<" "<<itersGD1<<", "<<itersGD2<<", "<<itersGD3<<", "<<timeneg<<", "<<elap<<", "<<timeverify<<","<<elapref<<" \n";;
 
	      exit(0);
            }
            else
            {

	      //elapsed += timer.TimePerTest(1);
	      //	      cout<<"\n Verification took "<<timer.TimePerTest(1)<<"s\n";
	      
	      cout<<"func="<<funcdef;
	      
	      cout<<"\n pass\n\n";
	      cout<<" "<<" "<<startval<<", ("<<startcoordx<<" "<<startcoordy<<" "<<startcoordz<<"), "<<setup<<", pass, "<<" "<<itersGD1<<", "<<itersGD2<<", "<<itersGD3<<", "<<iterstaken<<", "<<timeneg<<", "<<elap<<", "<<timeverify<<", "<<elapref<<"\n";;
 	      //WriteCSV(coeffs, 1);
	      
	      sol = phys; //(? or sol = phys)
            }
        }
        else{
	  timer.Stop();
	  //	  elapsed += timer.TimePerTest(1);
	  cout<<"\n checking for -vity took "<<timer.TimePerTest(1)<<"s\n\n";
	  
	  cout<<"\n optimizer no need\n\n";
        }
        //Backward transform solution to get projected values
        E->BwdTrans(coeffs, phys);

    }

    //Calculate L_inf & L_2 error
    cout << "L infinity error: \t" << E->Linf(phys, sol) << endl;
    if (stype != ePoint)
    {
        cout << "L 2 error: \t \t \t" << E->L2(phys, sol) << endl;
    }

    if (!vm.count("diff") && stype != ePoint)
    {
        //Evaluate solution at x = y = 0 and print error
        Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
        t[0] = -0.5;
        t[1] = -0.25;
        t[2] = -0.3;
        sol[0] = Shape_sol(t[0], t[1], t[2], order, btype, stype, false);

        NekDouble nsol = E->PhysEvaluate(t, phys);

        cout << "Error at x = (";
        for (int i = 0; i < dimension; ++i)
        {
            cout << t[i] << ", ";
        }
        cout << "\b\b): " << nsol - sol[0] << endl;
    }

    // Calculate integral of known function to test different quadrature
    // distributions on each element.
    for (int i = 0; i < totPoints; ++i)
    {
        sol[i] = dimension == 1 ? exp(x[i]) : dimension == 2 ?
            exp(x[i]) * sin(y[i]) : exp(x[i] + y[i] + z[i]);
    }

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

FILE *csvfile(int flag)
{
    static char name[LOGNAME_SIZE];
    char fullname[] = "Hex_";//[LOGNAME_SIZE+4];                                                                                                                                                           
    time_t now = time(0);
    strftime(name, sizeof(name), LOGNAME_FORMAT, localtime(&now));
    char * newArray = new char[std::strlen(fullname)+std::strlen(name)+10];
    std::strcpy(newArray,fullname);
    std::strcat(newArray,name);
    char integer_string[20];
    sprintf(integer_string, "_P%d_Q%d_%d.csv" , E->GetBasis(0)->GetNumModes(), E->GetBasis(0)->GetNumPoints(), flag);
    std::strcat(newArray,integer_string);
    cout<<"\n"<< newArray<<"\n\n";
    return fopen(newArray, "w");
}

void WriteCSV(Array<OneD, NekDouble> uhats, int flag)
{
    // if flag == 0, before opt
    // if flag == 1, after opt
    FILE *file = csvfile(flag);
    cout<<"\n1\n\n";

    int uhatstot = uhats.size();
    int totpts = demo.testcoord3d[0].size();
    Array<OneD, NekDouble> temp(uhatstot);
    for(int i = 0; i<totpts; i++)
    {
        Vmath::Vmul(uhatstot, &demo.interioreval3d[i], totpts, &uhats[0], 1, &temp[0], 1);
        NekDouble v = Vmath::Vsum(uhatstot, temp, 1);
        string s = std::to_string(demo.testcoord3d[0][i])+","+std::to_string(demo.testcoord3d[1][i])+","+std::to_string(demo.testcoord3d[2][i])+","+std::to_string(v)+"\n";
        const char* c = s.c_str();
        fprintf (file, c); 
    }
    fclose(file);
}

vector< vector<NekDouble> > arrtovec(Array<OneD, Array<OneD, NekDouble> > arr)
{
    vector<vector< NekDouble> > ret;

    for(int i = 0; i < arr.size(); i++)
    {
        vector<NekDouble> row;
        for(int k = 0; k < arr[0].size(); k++)
        {
            row.push_back(arr[i][k]);
        }
        ret.push_back(row);
    }
    return ret;
}

void project_surfaces( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret)
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

void surfuhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz)
{
    int modes = uhats.size();

    Array<OneD, NekDouble> temp(uhats.size());  

    int totpts = Equad->GetTotPoints();
    Array<OneD, NekDouble> vals(totpts);    
  
    //Vxyz*uhats -> project to -> E3tri
    for(int k = 0; k < totpts; k++)
    {    
        Vmath::Vmul(uhats.size(), &Vxyz[k], totpts, &uhats[0], 1, &temp[0], 1);
        vals[k]  = Vmath::Vsum(modes, temp, 1);
    }

    Equad->FwdTrans(vals, ret);
    

}



// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret , int d)
{
    int modes = (3*(E->GetBasis(0)->GetNumModes()-1));

    // assert d = 0 or 1
     
    if(d == 0)
    {    
        for(int k= 0; k < numedges; k++)
        {
            ret[k] = Array<OneD, NekDouble>(modes);        
        }
    }
    if(numedges == 12) // hex
    {

        if(d == 1)
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
        else if( d == 0)
        {
            // edge front left (x = -1) (y = -1)
            deruhatsedges(uhats, ret[0], Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z); 
            //edge front right (x = 1) (y = -1)
            deruhatsedges(uhats, ret[1],Vx1ym1z, Vdxx1ym1z, Vdyx1ym1z, Vdzx1ym1z);
            
            //edge front top (y = -1) (z = 1)
            deruhatsedges(uhats, ret[2], Vym1xz1, Vdxym1xz1, Vdyym1xz1, Vdzym1xz1);
    
            //edge front bot (y = -1) (z = -1)
            deruhatsedges(uhats, ret[3],  Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);
            
            // edge back left (y = 1), (x = -1)
            deruhatsedges(uhats, ret[4],Vxm1y1z, Vdxxm1y1z, Vdyxm1y1z, Vdzxm1y1z);
            //edge back right (x = 1), (y = 1))
            deruhatsedges(uhats, ret[5], Vx1y1z, Vdxx1y1z, Vdyx1y1z, Vdzx1y1z);
            //edge back top ( y = 1) (z = 1),
            deruhatsedges(uhats, ret[6], Vy1xz1, Vdxy1xz1, Vdyy1xz1, Vdzy1xz1);
   
            //edge back bot (y = 1) (z = -1))
            deruhatsedges(uhats, ret[7], Vy1xzm1, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);
            
                        
            // edge left bot (z = -1), (x = -1)
            deruhatsedges(uhats, ret[8],Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);
            
            //edge left top (x = -1), (z = 1))
            deruhatsedges(uhats, ret[9],  Vxm1yz1,Vdxxm1yz1, Vdyxm1yz1, Vdzxm1yz1);
            //edge right bot ( z = -1) (x = 1)
            deruhatsedges(uhats, ret[10],  Vx1yzm1, Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);
              
            //edge right top (z  1) (x  1))
            deruhatsedges(uhats, ret[11],Vx1yz1, Vdxx1yz1, Vdyx1yz1, Vdzx1yz1); 


        }
    }
   
}


void deruhatsedges(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz, Array<OneD, NekDouble> Vdxyz, Array<OneD, NekDouble> Vxdyz, Array<OneD, NekDouble> Vxydz)
{

    boost::ignore_unused(Vxyz);
    int totpts =  edgeptsin[0].size();//E3seg->GetTotPoints();
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
        
        Vmath::Vmul(uhats.size(), &Vxydz[k], totpts, &uhats[0], 1, &temp[0], 1);
        vals[k]  = v1 + Vmath::Vsum(modes, temp, 1);
    }    
    
    E3seg->FwdTrans(vals,ret);
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
 
        
        //    }
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


// sig = 0 -> called by opt_needed
// sig = 1 -> called by sphere_rotation
Array<OneD, Array<OneD,  NekDouble> > call_find_roots(Array<OneD,  NekDouble> &uhats , NekDouble &avgiterGD, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats, NekDouble &roots1dtime, NekDouble &roots2dtime, NekDouble &roots3dtime,int sig )
{

    boost::ignore_unused( surfaceuhats);
    int dimension = E->GetShapeDimension(); 
    vector<vector< NekDouble> > ret(dimension);
    NekDouble dummy;
    Array<OneD, Array<OneD, NekDouble> > retarr;//(dimension);
    roots1dtime = 0.0;
    roots2dtime = 0.0;
    roots3dtime = 0;
    Timer t;
        if(numedges == 12) //hex
        {
            for(int ii = 0; ii < numedges; ii++)
            {

                vector<vector<NekDouble> > tmp(dimension);
		t.Start();
                tmp  = (demo.find_roots(uhatsedges[ii], E, dummy,  0, sig, 0, 0, 0)) ;
		t.Stop();
		roots1dtime +=  t.TimePerTest(1);
                for(int p = 0; p < tmp[0].size(); p++)
                {
                    if(ii == 0) // edge front left (x = -1) (y = -1)
                    {
                        ret[2].push_back(tmp[0][p]);
                        ret[1].push_back(-1);
                        ret[0].push_back(-1);
                        
                    }
                    else if(ii == 1)//edge front right (x = 1) (y = -1)
                    {
                        ret[0].push_back(1);
                        ret[1].push_back(-1);
                        ret[2].push_back(tmp[0][p]);
                    }
                    else if(ii == 2)//edge front top (y = -1) (z = 1)
                    {
                        ret[0].push_back(tmp[0][p]);
                        ret[1].push_back(-1);
                        ret[2].push_back(1);
                    
                    }
                    else if(ii == 3) //edge front bot (y = -1) (z = -1)
                    {
                        ret[1].push_back(-1);
                        ret[0].push_back(tmp[0][p]);
                        ret[2].push_back(-1);
                    }
                    else if(ii == 4) //edge back left (y = 1), (x = -1)
                    {
                        ret[0].push_back(-1);
                        ret[1].push_back(1);
                        ret[2].push_back(tmp[0][p]);
                      
      
		    }
                    else if(ii == 5) //edge back right (x = 1), (y = 1)
                    {
                        ret[0].push_back(1);
                        ret[1].push_back(1);
                        ret[2].push_back(tmp[0][p]);
                    }
                    else if(ii == 6) //edge back top ( y = 1) (z = 1)
                    {
                        ret[1].push_back(1);
                        ret[0].push_back(tmp[0][p]);
                        ret[2].push_back(1);
                    }
                    
                    else if(ii == 7) //edge back bot (y = 1) (z = -1)
                    {
                        ret[1].push_back(1);
                        ret[0].push_back(tmp[0][p]);
                        ret[2].push_back(-1);
                    }
                    else if(ii == 8) //edge left bot (z = -1), (x = -1)
                    {
                        ret[1].push_back(tmp[0][p]);
                        ret[2].push_back(-1);
                        ret[0].push_back(-1);
                    }
                    else if(ii == 9) //edge left top (x = -1), (z = 1)
                    {
                        ret[2].push_back(1);
                        ret[1].push_back(tmp[0][p]);
                        ret[0].push_back(-1);
        
 
                    }
                    else if(ii == 10) //edge right bot (z = -1) (x = 1)
                    {
                        ret[1].push_back(tmp[0][p]);
                        ret[2].push_back(-1);
                        ret[0].push_back(1);
                    }
                    else if(ii == 11) //edge right top (z  1) (x  1)
                    {
                        ret[0].push_back(1);
                        ret[1].push_back(tmp[0][p]);
                        ret[2].push_back(1);
                    }
                    
                }
            

            
            }
        }
        
        if(numsurfaces == 6)
        {   

            // call 2D rootfinder on each surface:
            vector<vector<NekDouble> > tmp(dimension);
	    NekDouble avgiterGDhold;
            for(int ii = 0; ii < numsurfaces; ii++)
            {
	      t.Start();
	      tmp = demo.find_roots(surfaceuhats[ii], Equad, avgiterGDhold, 0, sig, 1, 0 , 0); 
	      t.Stop();
	      roots2dtime +=  t.TimePerTest(1);  	      
	      avgiterGD += avgiterGDhold;
	      if(ii == 0) //surf bot (z=-1)
                {
                    for(int p = 0; p < tmp[0].size(); p++)
                    {
                        ret[2].push_back(-1);
                        ret[1].push_back(tmp[1][p]);
                        ret[0].push_back(tmp[0][p]);
                    }
                }
                else if(ii == 1) //surf right (x = 1)
                {
                    for(int p = 0; p < tmp[0].size(); p++)
                    {
                        ret[0].push_back(1);
                        ret[2].push_back(tmp[1][p]);
                        ret[1].push_back(tmp[0][p]);
                    }
                }
                else if(ii == 2) //surf top (z = 1)
                {
		  for(int p = 0; p < tmp[0].size(); p++)
                    {
                        ret[1].push_back(tmp[1][p]);
                        ret[2].push_back(1);
                        ret[0].push_back(tmp[0][p]);
                    }
                }
                else if(ii == 3) //surf left (x = -1)
                {
                    for(int p = 0; p < tmp[0].size(); p++)
                    {
                        ret[0].push_back(-1);
                        ret[2].push_back(tmp[1][p]);
                        ret[1].push_back(tmp[0][p]);
                    } 
                }
                else if(ii == 4) //surf front (y = -1) 
                {
                    for(int p = 0; p < tmp[0].size(); p++)
                    {
                        ret[0].push_back(tmp[0][p]);
                        ret[2].push_back(tmp[1][p]);
                        ret[1].push_back(-1);
                    }
                }
                else if(ii == 5) //surf back (y = 1)
                {
                    for(int p = 0; p < tmp[0].size(); p++)
                    {

                        ret[0].push_back(tmp[0][p]);
                        ret[1].push_back(1);
                        ret[2].push_back(tmp[1][p]);
                    }
                }                     
                
            }
        }
            
        //find interior roots:
        vector<vector<NekDouble> > tmp;
	NekDouble avgiterGDhold;
	t.Start();
        tmp =  (demo.find_roots(uhats, E, avgiterGDhold, 0, sig, 0, 1, 0)) ;
	t.Stop();
	roots3dtime +=  t.TimePerTest(1);  
	avgiterGD += avgiterGDhold;
        if(tmp[0].size() > 0)
        {
            for(int p = 0; p <dimension; p++)
            {
                ret[p].push_back(tmp[p][0]);
            }
        }
    //add vertices:
	//vertex E
	ret[0].push_back(-1.0);
	ret[1].push_back(-1.0);
	ret[2].push_back(1.0);

	//vertex A
	ret[0].push_back(-1.0);
	ret[1].push_back(-1.0);
	ret[2].push_back(-1.0);

	//veryex B
	ret[0].push_back(1.0);
	ret[1].push_back(-1.0);
	ret[2].push_back(-1.0);

	//vertex C
	ret[0].push_back(1.0);
	ret[1].push_back(1.0);
	ret[2].push_back(-1.0);

	//vertex D
	
	ret[0].push_back(-1.0);
	ret[1].push_back(1.0);
	ret[2].push_back(-1.0);

	//vertex H
	
	ret[0].push_back(-1.0);
	ret[1].push_back(1.0);
	ret[2].push_back(1.0);

	//vertex G
	ret[0].push_back(1.0);
	ret[1].push_back(1.0);
	ret[2].push_back(1.0);
	
	//vertex F
	ret[0].push_back(1.0);
	ret[1].push_back(-1.0);
	ret[2].push_back(1.0);

    
    retarr = vectoarr(ret);
    return retarr;

}
    

int Opt_needed(Array<OneD, NekDouble> uhats, int flag)
{
    
    int totModes = uhats.size();
    //vector<vector<  NekDouble> > roots;
    Array<OneD,  Array<OneD, NekDouble> > rootsarr;
    Array<OneD, Array<OneD, NekDouble> > edgeuhats  (numedges);

    Array<OneD, Array<OneD, NekDouble> > surfacederuhats  (numsurfaces);
    for(int t = 0; t <numsurfaces; t++)
    {
        surfacederuhats[t] = Array<OneD, NekDouble>(Equad->GetNcoeffs());
        
    }
    project_surfaces(uhats, surfacederuhats);
    
    //find edge roots: (edges = 12 is hex)

    project_edges(uhats, edgeuhats);

    //find roots of der of uhats, that is why d = 0
    NekDouble dummy, dummy1,dummy2,dummy3;
    rootsarr = call_find_roots(uhats, dummy, edgeuhats, surfacederuhats, dummy1,dummy2,dummy3);
    cout<<"\n to chck -vity or verify, iters of GD = "<<dummy<<" ";
    if(flag == 0)
      itersGD1 = dummy;
    else
      itersGD3 = dummy;
    // evaluate ortho basis at roots
    // evalBasisRoots is flattened basis eval matrix
    int     evalsz = (rootsarr[0].size())*totModes;
    Array<OneD, NekDouble> evalBasisRoots(evalsz); 

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
         
    cout<<"\n minv = "<<minv<<" at ";
    for(int k = 0; k < rootsarr.size(); k++)
    {
        cout<<rootsarr[k][idx]<<" ";

    }
    if(minv < 0.0 && abs(minv)>1e-10)
    {
        if(flag == 0)
        {
            startval = minv;
            startcoordx = rootsarr[0][idx];
            startcoordy = rootsarr[1][idx];
            startcoordz = rootsarr[2][idx];
            
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
    Array<OneD, Array<OneD,NekDouble> > optima;
    //    int NC = 1; //number of constraints, only positivity for now
    tols.push_back(1e-12);
    
    Array<OneD, NekDouble>  pqvalcoords(dim+1), xastarr(dim), utemp(N1), wsp1;   
    

    Array<OneD, Array<OneD, NekDouble> > surfaceuhats  (numsurfaces);
    Array<OneD, Array<OneD, NekDouble> > Pf(numedges);
    for(int k= 0; k < numedges; k++)
    {
            
        Pf[k] = Array<OneD, NekDouble>(E3seg->GetNcoeffs());//(3*(E->GetBasis(0)->GetNumModes()-1));        
    }
    for(int k = 0; k < numsurfaces; k++)
        {
            surfaceuhats[k] = Array<OneD, NekDouble>(Equad->GetNcoeffs());
        }
    NekDouble pqval, timeprojectsurf, timeprojectedges;
    NekDouble avgiterGD = 0.0, avgiterhold = 0.0, roots1dtimehold, roots2dtimehold, roots3dtimehold, evalrootshold;   
    //    cout<<"\n avgiterGD ="<<avgiterGD<<" avgiterhold="<<avgiterhold;
    Timer t;
    while (counter <= niter)
    {
      NekDouble roots1dtime, roots2dtime, roots3dtime ;
      pqval = inf;
      utemp = d.back();
      
      t.Start();
      //pq = @(xx) cfun(xx) .* (constraints(1,jj).*Vc(xx)*uhats);  
      project_edges(utemp, Pf, 1);
      t.Stop();
      timeprojectedges += t.TimePerTest(1);

	t.Start();
        project_surfaces(utemp, surfaceuhats);
	t.Stop();
	timeprojectsurf += t.TimePerTest(1);

	optima = (call_find_roots(utemp, avgiterhold, Pf, surfaceuhats,  roots1dtime, roots2dtime, roots3dtime,1));
	roots1dtimehold += roots1dtime;
	roots2dtimehold += roots2dtime;

	roots3dtimehold += roots3dtime;
	
	
	//	cout<<"\n avgiterhold="<<avgiterhold<<" ";
	avgiterGD += avgiterhold;    
        wsp1 = Array<OneD, NekDouble>(optima[0].size());
        // cout<<"\n optima:\n";
        // for(int pp = 0; pp < optima.size(); pp++)
        // {
        //     for(int jj = 0; jj < optima[pp].size(); jj++)
        //     {
        //         cout<<optima[pp][jj]<<" ";
        //     }
        //     cout<<"\n";
        // }
        //        cout<<"\n";
        // Vtmp is evaluation of basis fn at optima
	t.Start();
        Array<OneD, NekDouble> Vtmp(optima[0].size() * E->GetNcoeffs());
        E->PhysEvalBasisGrad(optima, Vtmp, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        
	demo.pq(utemp, optima, Vtmp, pqvalcoords, wsp1);
        // cout<<"\n fvals: "<<wsp1.size()<<"\n";;
        // for(int kk = 0 ; kk<wsp1.size(); kk++)
        //     cout<<" "<<wsp1[kk]<<" ";
        // cout<<"\n";

        if (pqvalcoords[0] < pqval)
        {
            for(int k = 0; k  < dimension; k++)
            {
                xastarr[k] = pqvalcoords[k+1];
            }
            pqval = pqvalcoords[0];
        }
        
        cout<<"\n at counter="<<counter<<" min val="<<pqval<<" xast ="<<xastarr[0]<<" "<<xastarr[1]<<" "<<xastarr[2]<<
            " ";
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
	t.Stop();
	evalrootshold+= t.TimePerTest(1);  

        counter = counter + 1;
    
    
    }
    roots1dtimehold = roots1dtimehold/(counter+1);
    roots2dtimehold = roots2dtimehold/(counter+1);
    roots3dtimehold = roots3dtimehold/(counter+1);
    timeprojectedges = timeprojectedges/(counter+1);
    timeprojectsurf = timeprojectsurf/(counter+1);    
    itersGD2 = (avgiterGD)/(counter+1);
    iterstaken = counter;
    cout<<"sphere_rotation took "<<counter<<"startv = "<<startval<<" at startcoords("<<startcoordx<<","<<startcoordy<<" ," <<startcoordz<<")   iterations\n  GD iter stepszie ="<<demo.GetTol()<<" maxiters = "<<demo.GetMaxIter()<<" o = "<<E->GetBasis(0)->GetNumModes()<<" p = "<<E->GetBasis(0)->GetNumPoints()<<" eps = "<<demo.GetEps()<<"  avgiterGD = "<<itersGD2<<" total time taken post rootfinding per iter = "<< evalrootshold/(counter)<<" hex \n";;
    cout<<"\nAvg times per iter (edgeprojection,surfprojection,1drootfinder,2drootfinder,3drootfinder) = : "<< timeprojectedges<<", "<<timeprojectsurf<<", " << roots1dtimehold<<", "<<roots2dtimehold<<", "<<roots3dtimehold ;
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
