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
#include <fstream>
#include <time.h>
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


void  pq(
         Array<OneD,NekDouble> uhats,
         Array<OneD, Array<OneD ,NekDouble > >roots,
         Array<OneD,NekDouble> V1,
         Array<OneD,NekDouble> &pqevalxast = NullNekDouble1DArray,
         Array<OneD,NekDouble>&fvals = NullNekDouble1DArray);


//declare Do_optimize
void Do_optimize(Array<OneD, NekDouble> &uhats);

Array<OneD, Array<OneD, NekDouble> > vectoarr(vector<vector< NekDouble> > vec);


//declare find_roots, flag is 0 for confederate matrix approach
vector<vector<  NekDouble> > find_roots(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> evalxyz = NullNekDouble1DArray,  int d = 0, int flag = 0,  int sig = 0, int surfflag = 0, int surfid = 0);

// declare caller routine to find_roots
// flag = 0 -> opt_needed calls
// flag = 1 -> sphere_rot calls
Array<OneD, Array< OneD,  NekDouble> >  call_find_roots(Array<OneD, NekDouble> &uhatsall, int d = 0, Array< OneD, Array<OneD, NekDouble> >&uhatsedges =NullNekDoubleArrayofArray , Array< OneD, Array<OneD, NekDouble> >&surfaceuhats =NullNekDoubleArrayofArray ,  int flag = 0 );


//declare Opt_needed
int Opt_needed(Array<OneD, NekDouble> uhats, int flag = 0);


// for 2D elements to get uhats at edges
// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble> uhats,Array<OneD, Array<OneD, NekDouble> >&ret , int d = 0);

void derpq(Array<OneD, NekDouble> &uhats,  NekDouble &ret,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd);

// uhatpqd stuff: 
// called by project_edges if d = 1 is passed in project_edges
void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret, int modes, Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1 = NullNekDouble1DArray , Array<OneD, NekDouble> Vxyd2 = NullNekDouble1DArray);

void surfacepquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret, Array<OneD, NekDouble> Vxy, int surfid = 0);

void surfuhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz, int surfid = 0 );


void deruhatsedges(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble >&ret,  Array<OneD, NekDouble> Vxyz, Array<OneD, NekDouble> Vdxyz, Array<OneD, NekDouble> Vxdyz, Array<OneD, NekDouble> Vxydz);

void project_surfaces( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&retquad, Array<OneD, Array<OneD, NekDouble> >&rettri, int d = 0);

vector< vector<NekDouble> > gradient_descent(Array<OneD, NekDouble> uhats, Array<OneD, NekDouble> evalxyz, int sig, int surfflag , int surfid = 0);


FILE *csvfile(int flag);
void WriteCSV( Array<OneD, NekDouble>  uhats, int flag = 0);


string funcdef;
// Colleague matrix
int numedges, numsurfaces;
Array<OneD, Array<OneD, NekDouble> > C;
StdExpansion *E;
StdExpansion *E3seg;
StdExpansion *Eseg;
StdExpansion *Equad;
StdExpansion *Etri;

//StdExpansion *E3quad;
Array<OneD, NekDouble> qZin; //inner point grid for 2D rootfinding
Array<OneD, Array<OneD, NekDouble> > edgeptsin;
Array<OneD, NekDouble> qZinmid; //inner mid point grid for 2D rootfinding
Array<OneD, NekDouble> qWin; //inner point grid for 2D rootfinding
Array<OneD, NekDouble> qx;


Array<OneD, Array<OneD, NekDouble> >qxarrquad;
Array<OneD, Array<OneD, NekDouble> > qxarrtri;
Array<OneD, Array<OneD, Array<OneD, NekDouble> > > surfcoords;
Array<OneD, NekDouble> qw;

Array<OneD, NekDouble> V;
Array<OneD, NekDouble> Vd;
//Array<OneD, NekDouble> V3;
Array<OneD, NekDouble > qWinarr;
Array<OneD, NekDouble> Vquad;

// for prism edges root finding

// edge front left (x = -1) (y = -1)
Array<OneD, NekDouble> Vxm1ym1z   ;
Array<OneD, NekDouble> Vdyxm1ym1z   ;
Array<OneD, NekDouble> Vdxxm1ym1z   ;
Array<OneD, NekDouble> Vdzxm1ym1z   ;
               
//edge front hypt (y = -1) (z = -x)
Array<OneD, NekDouble> Vym1xmz   ;
Array<OneD, NekDouble> Vdxym1xmz   ;
Array<OneD, NekDouble> Vdyym1xmz   ;
Array<OneD, NekDouble> Vdzym1xmz   ;
               
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
               
//edge back hypt ( y = 1) (z = -x)
Array<OneD, NekDouble> Vy1xmz   ;
Array<OneD, NekDouble> Vdxy1xmz   ;
Array<OneD, NekDouble> Vdyy1xmz   ;
Array<OneD, NekDouble> Vdzy1xmz   ;
               
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
 
//surface bot z = -1
Array<OneD, NekDouble>Vxyzm1 ;


//surface left x = -1
Array<OneD, NekDouble>Vxm1yz ;


//surface hypt x+z = 0
Array<OneD, NekDouble>Vxzmy ;


//surface front y = -1
Array<OneD, NekDouble>Vxym1z ;

//surface back y = 1
Array<OneD, NekDouble>Vxy1z ;

Array<OneD, NekDouble> Vxy1zproject;  
Array<OneD, NekDouble> Vxym1zproject;  

int dimension ;
Array<OneD, NekDouble > interioreval3d;
Array<OneD, NekDouble > interioreval2dquad;
Array<OneD, NekDouble > interioreval2dtri;
Array<OneD, Array<OneD, NekDouble> > testcoord3d; 
Array<OneD, Array<OneD, NekDouble> > testcoord2dquad; 
Array<OneD, Array<OneD, NekDouble> > testcoord2dtri; 
NekDouble iterGD, secarg, startval, startcoordx, startcoordy, startcoordz;
            
DemoSupport demo;

int main(int argc, char *argv[])
{
    demo.GetOptions().add_options()("optm,z", "positivity preserving optimizer");
    demo.ParseArguments(argc, argv);
    
    po::variables_map vm = demo.GetVariableMap();

    //only for 1D
       
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
    LibUtilities::PointsType pointsTypetriB = LibUtilities:: ePolyEvenlySpaced;    
    LibUtilities::PointsType pointsTypetriA = LibUtilities::ePolyEvenlySpaced;    

    for (int i = 0; i < dimension; ++i)
    {
        btype[i] = E->GetBasisType(i);
        order.push_back(E->GetBasisNumModes(i));
    }

    C = demo.formConf((dimension)*pow(3*order[0]+1,2));
    
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

    case LibUtilities::ePrism:
        
        numedges = 9;
        numsurfaces = 5;
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
        if(dimension ==3) //hex
        {
            sol[i] =sin(2*x[i])*sin(2*y[i])*sin(2*z[i]) ;///sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(x[i])*sin(y[i])*sin(z[i]);//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(x[i])*sin(y[i])*sin(z[i]);//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(x[i])*sin(y[i])*sin(z[i]);//sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);//sin(x[i])*sin(y[i])*sin(z[i])+0.1;//sin(x[i])*sin(y[i])*sin(z[i]);// +0.1;//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(x[i])*sin(y[i])*sin(z[i])+0.1;//sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);//sin(x[i])*sin(y[i])*sin(z[i]);//sin(x[i])*sin(y[i])*sin(z[i]) +0.1;//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);//sin(x[i])*sin(y[i])*sin(z[i]);//sin(x[i])*sin(y[i])*sin(z[i]) ;//sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);//sin(x[i])*sin(y[i])*sin(z[i]);//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2;////sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1;
            funcdef ="sin(2*x[i])*sin(2*y[i])*sin(2*z[i])";;//sin(x[i])*sin(y[i])*sin(z[i]);";//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.5";//sin(x[i])*sin(y[i])*sin(z[i]);";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])";//"sin(x[i])*sin(y[i])*sin(z[i])+0.1";//"sin(x[i])*sin(y[i])*sin(z[i]) ";//+0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i])+0.1";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(x[i])*sin(y[i])*sin(z[i]) +0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])";//"sin(x[i])*sin(y[i])*sin(z[i]) ";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i]);";//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1;";//"sin(x[i])*sin(y[i])*sin(z[i];";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2 ";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.5";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);";// "sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])  -2";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";// -  exp(-1) + 1.37";
	}
    }

    Array<OneD, NekDouble> phys(totPoints);
    Array<OneD, NekDouble> coeffs((unsigned) E->GetNcoeffs());
    
    //Project onto expansion
    E->FwdTrans(sol, coeffs);

    //Backward transform solution to get projected values
    E->BwdTrans(coeffs, phys);
    
    if (vm.count("optm"))
    {
            

        int dimension = E->GetShapeDimension(); 
        //qx = E->GetBasis(0)->GetZ();
        //qw = E->GetBasis(0)->GetW();            

        vector<LibUtilities::BasisKey> tmpBasisKey = demo.GetBasisKey();       
        vector<LibUtilities::PointsKey> tmpPtsKey = demo.GetPointsKey();       
        cout<<"\n E->GetBasis(0)->GetZ().size()="<<E->GetBasis(0)->GetZ().size()<<" E->GetBasis(1)->GetZ() sz = "<<E->GetBasis(1)->GetZ().size()<<" E->GetBasis(2)->GetZ() sz = "<<E->GetBasis(2)->GetZ().size()<<"\n";
        LibUtilities::PointsKey pkeycheb(E->GetBasis(0)->GetNumPoints(), pointsTypeCheb);
        LibUtilities::PointsKey pkeytriA(3*(E->GetBasis(1)->GetZ()).size(), pointsTypetriA);
        LibUtilities::PointsKey pkeytriB(3*(E->GetBasis(2)->GetZ()).size(), pointsTypetriB);

        LibUtilities::BasisKey bkeycheb(LibUtilities::eChebyshev, 3*(E->GetBasis(0)->GetNumModes()),  pkeycheb);
        LibUtilities::BasisKey bkeytriB(LibUtilities::eOrtho_B,(E->GetBasis(1)->GetTotNumModes()), pkeytriB);
        LibUtilities::BasisKey bkeytriA(LibUtilities::eOrtho_A, (E->GetBasis(0)->GetTotNumModes()), pkeytriA);
        E3seg = new StdSegExp(bkeycheb);
       
        if (E3seg == nullptr)
        {
            return 1;
        }
        LibUtilities::BasisKey bkeycheb3 (LibUtilities::eChebyshev, (E->GetBasis(0)->GetNumModes()),  pkeycheb);

        Eseg = new StdSegExp(bkeycheb3);
       
        if (Eseg == nullptr)
        {
            return 1;
        }
        LibUtilities::BasisKey bkeycheb2 (LibUtilities::eChebyshev, (E->GetBasis(0)->GetNumModes()),  pkeycheb);
        Equad = new StdQuadExp(bkeycheb2, bkeycheb2);//tmpBasisKey[0],tmpBasisKey[0]);
        if (Equad == nullptr)
        {
            return 1;
        }
        Etri = new StdTriExp(bkeytriA, bkeytriB);//tmpBasisKey[0], tmpBasisKey[2]);        

        if (Etri == nullptr)
        {
            return 1;
        }


        if(dimension > 1 )
        {

        
            // populate edge root finding matrices for quadrilateral
            LibUtilities::PointsKey  quadPointsKeyin  (order[0], pointsTypeCheb);
            // for gradient descent
            qZin = (LibUtilities::PointsManager()[quadPointsKeyin ])->GetZ();
            int sz = pow(qZin.size(),2);;
            int szmid = pow(qZin.size()-1,2);
            qZinmid = Array<OneD,NekDouble>(qZin.size()-1);
            
            for(int i = 0; i <qZinmid.size(); i++)
            {
                qZinmid[i] = (qZin[i]+qZin[i+1])/2;
            }


            testcoord2dquad = Array<OneD, Array<OneD, NekDouble> > (2); 
            testcoord2dtri = Array<OneD, Array<OneD, NekDouble> > (2); 
            Array<OneD, Array<OneD, NekDouble> > testcoord2dtritemp(2);
               
            int ctsz2dquad = 0;
            int sztri = (Etri->GetBasis(0)->GetZ()).size()*(Etri->GetBasis(1)->GetZ().size())  ;
            int szmidtri = ((Etri->GetBasis(0)->GetZ()).size()-1)*((Etri->GetBasis(1)->GetZ()).size() - 1)  ;
            
            for(int ii = 0; ii < 2; ++ii)
            {
                testcoord2dquad[ii] = Array<OneD, NekDouble>(sz+szmid);
                testcoord2dtritemp[ii] = Array<OneD, NekDouble>(sztri+szmidtri);
            }

            for(int ii = 0; ii < qZin.size(); ++ii)
            {
                
                for(int jj = 0; jj < qZin.size(); ++jj)
                {
                    testcoord2dquad[0][ctsz2dquad] = qZin[ii];
                    testcoord2dquad[1][ctsz2dquad] = qZin[jj];
                        ctsz2dquad++;
                }
                
            }     
            
            // add midpt grid
            for(int ii = 0; ii < qZinmid.size(); ++ii)
            {
                for(int jj = 0; jj < qZinmid.size(); ++jj)
                {
                
                    testcoord2dquad[0][ctsz2dquad] = qZinmid[ii];
                    testcoord2dquad[1][ctsz2dquad] = qZinmid[jj];

                    ctsz2dquad++;
                    
                }
                
            }
       

            int cts2d = 0;
            for(int ii = 0; ii < (Etri->GetBasis(0)->GetZ()).size(); ++ii)
            {
                for(int jj = 0; jj < (Etri->GetBasis(1)->GetZ()).size();  ++jj)
                {
                    if( (Etri->GetBasis(0)->GetZ())[ii] +  (Etri->GetBasis(1)->GetZ())[jj] <=0)
                    {
                        testcoord2dtritemp[0][cts2d] = (Etri->GetBasis(0)->GetZ())[ii];
                        testcoord2dtritemp[1][cts2d] =  (Etri->GetBasis(1)->GetZ())[jj];
                        cts2d++;
                    }
                }
                
            }   

            // add midpt grid
            for(int ii = 0; ii < Etri->GetBasis(0)->GetZ().size() - 1; ++ii)
            {
                for(int jj = 0; jj < Etri->GetBasis(1)->GetZ().size() - 1; ++jj)
                {
                    NekDouble c1 = ((Etri->GetBasis(0)->GetZ())[ii] + (Etri->GetBasis(0)->GetZ())[ii+1])/2;
 
                    NekDouble c2 =  ((Etri->GetBasis(1)->GetZ())[jj] + (Etri->GetBasis(1)->GetZ())[jj+1])/2;
                    if(c1+c2 <= 0 )
                    {
                        testcoord2dtritemp[0][cts2d] = c1;
                        testcoord2dtritemp[1][cts2d] = c2;
          
                        cts2d++;
                    }
                }
        
            }  

            Array<OneD, Array<OneD, NekDouble> > temptri(2);

            temptri[0] =  Array<OneD, NekDouble>(cts2d);//testcoord2dtritemp[0].size());
            temptri[1] =  Array<OneD, NekDouble>(cts2d);//testcoord2dtritemp[1].size());
            int p = 0;

            for(int k = 0; k < temptri[0].size(); k++)
            {
                if(testcoord2dtritemp[0][k]+testcoord2dtritemp[1][k]<=0)
                {
                    
                    temptri[0][p] = testcoord2dtritemp[0][k];
                    temptri[1][p] = testcoord2dtritemp[1][k];
                                   
                    p++;
                }
                
            }

            for(int ii = 0; ii < 2; ++ii)
            {
                testcoord2dtri[ii] = Array<OneD, NekDouble>(p);
                Vmath::Vcopy(p, &temptri[ii][0], 1, &testcoord2dtri[ii][0], 1);
                
            }

            LibUtilities::BasisKey bkeychebtmp(LibUtilities::eChebyshev, tmpBasisKey[0].GetNumModes(),  pkeycheb);
	   

            sz = pow((qZin.size()+qZinmid.size()),3);;

            Array<OneD, Array<OneD, NekDouble> >testcoord3dtemp  (dimension); 
            testcoord3d = Array<OneD, Array<OneD, NekDouble> > (dimension); 
            for(int ii = 0; ii < 3; ++ii)
            {
                testcoord3dtemp[ii] = Array<OneD, NekDouble>(sz);
            }
            int ctsz = 0;
            for(int ii = 0; ii < qZin.size(); ++ii)
            {
                    
                for(int jj = 0; jj < qZin.size(); ++jj)
                {

                    for(int kk = 0; kk < qZin.size(); ++kk)
                    {
                        if(qZin[ii] + qZin[kk] <=0) 
                        {
                            testcoord3dtemp[0][ctsz] = qZin[ii];
                            testcoord3dtemp[1][ctsz] = qZin[jj];
                            testcoord3dtemp[2][ctsz] = qZin[kk];
                            ctsz++;
                        }
                    }
                }
            }
            for(int ii = 0; ii < qZinmid.size(); ++ii)
            {
                    
                for(int jj = 0; jj < qZinmid.size(); ++jj)
                {

                    for(int kk = 0; kk < qZinmid.size(); ++kk)
                    {
                        if(qZinmid[ii] + qZinmid[kk] <=0) 
                        {
                            testcoord3dtemp[0][ctsz] = qZinmid[ii];
                            testcoord3dtemp[1][ctsz] = qZinmid[jj];
                            testcoord3dtemp[2][ctsz] = qZinmid[kk];
                            ctsz++;
                        }
                    }
                }
            }

            for(int ii = 0; ii < dimension; ++ii)
            {
                testcoord3d[ii] = Array<OneD, NekDouble>(ctsz);
                Vmath::Vcopy(ctsz, &testcoord3dtemp[ii][0], 1, &testcoord3d[ii][0], 1);

            }

            interioreval2dquad = Array<OneD, NekDouble >(Equad->GetNcoeffs()*testcoord2dquad[0].size());
            Equad->PhysEvalBasisGrad(testcoord2dquad, interioreval2dquad , NullNekDouble1DArray , NullNekDouble1DArray, NullNekDouble1DArray);
                        
            interioreval2dtri = Array<OneD, NekDouble >(Etri->GetNcoeffs()*testcoord2dtri[0].size());
            // int tmpctr = 0;
            // cout<<"\n 1 Etri interioreval2dtri sz ="<<interioreval2dtri.size()<<" testcoord2dtri sz = "<<testcoord2dtri[0].size()<<"Etri->GetBasis(0)->GetZ() sz =="<<(Etri->GetBasis(0)->GetZ()).size()<<" (Etri->GetBasis(1)->GetZ()).size()= "<<(Etri->GetBasis(1)->GetZ()).size()<<" tmpctr="<<tmpctr<<"\n\n";
            // tmpctr++;         
            // for(int p = 0; p < Etri->GetNcoeffs(); p++)
            // {
            //     cout<<"\n p = "<<p<<"\n";
            //     for(int k = 0; k < testcoord2dtri[0].size(); k++)
            //     {
            //         cout<<"\n k = "<<k<<"\n\n";
            //         Array<OneD, NekDouble> pts(2);
            //         pts[0] = testcoord2dtri[0][k];
            //         pts[1] = testcoord2dtri[1][k];
            //         interioreval2dtri[k*p] = Etri->PhysEvaluateBasis(pts, p);
            //     }
            // }
            Etri->PhysEvalBasisGrad(testcoord2dtri, interioreval2dtri , NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray); 

            interioreval3d = Array<OneD, NekDouble >(E->GetNcoeffs()*testcoord3d[0].size());
                
            E->PhysEvalBasisGrad(testcoord3d, interioreval3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
            
            
            if(numedges == 9) //prism
            {
                //we need qzin + qzinmid as input points
                Array<OneD, NekDouble> edgexyztemp =   E3seg->GetBasis(0)->GetZ();

                int totszedges1d =  edgexyztemp.size()*(E->GetNcoeffs());//E->GetNcoeffs()*(E->GetBasis(0)->GetNumPoints());//(qZin.size()+qZinmid.size())*E->GetNcoeffs(); 
                //                Array<OneD, Array<OneD, NekDouble> >  edgeptsintemp(dimension);
               
                edgeptsin = Array<OneD, Array<OneD, NekDouble> >(dimension);         
                for(int p = 0; p < dimension; p++)
                {
                    edgeptsin[p] =   Array<OneD, NekDouble>(edgexyztemp);            
                }

                // edge front left (x = -1) (y = -1)
                Vxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
                Vdxxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
                Vdyxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
                Vdzxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
                edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
                edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

                E->PhysEvalBasisGrad(edgeptsin, Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);  
                
                //edge front hypt (y = -1) (z + x = 0)
                Vym1xmz = Array<OneD, NekDouble>(totszedges1d);
                Vdxym1xmz = Array<OneD, NekDouble>(totszedges1d);
                Vdyym1xmz = Array<OneD, NekDouble>(totszedges1d);
                Vdzym1xmz = Array<OneD, NekDouble>(totszedges1d);
                Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[2][0], 1, &edgeptsin[0][0], 1);

                E->PhysEvalBasisGrad(edgeptsin,Vym1xmz, Vdxym1xmz,Vdyym1xmz,Vdzym1xmz);  

                
                //edge front bot (y = -1) (z = -1)
                Vym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
                Vdxym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
                Vdyym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
                Vdzym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
                edgeptsin[0] = edgexyztemp;
                edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);           
            
                E->PhysEvalBasisGrad(edgeptsin, Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);  

                // edge back left (y = 1), (x = -1)
                Vxm1y1z = Array<OneD, NekDouble>(totszedges1d);
                Vdxxm1y1z = Array<OneD, NekDouble>(totszedges1d);
                Vdyxm1y1z = Array<OneD, NekDouble>(totszedges1d);
                Vdzxm1y1z = Array<OneD, NekDouble>(totszedges1d);

                edgeptsin[2] = edgexyztemp;
                edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
                edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);           
                E->PhysEvalBasisGrad(edgeptsin,Vxm1y1z, Vdxxm1y1z, Vdyxm1y1z, Vdzxm1y1z);  

                //edge back hypt ( y = 1) (z +x = 0)
                Vy1xmz = Array<OneD, NekDouble>(totszedges1d);
                Vdxy1xmz = Array<OneD, NekDouble>(totszedges1d)   ;
                Vdyy1xmz = Array<OneD, NekDouble>(totszedges1d)    ;
                Vdzy1xmz =Array<OneD, NekDouble>(totszedges1d)    ;

                Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[2][0], 1, &edgeptsin[0][0], 1);
                E->PhysEvalBasisGrad(edgeptsin,Vy1xmz, Vdxy1xmz, Vdyy1xmz, Vdzy1xmz);  
                
                //edge back bot (y = 1) (z = -1))
                Vy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
                Vdyy1xzm1 = Array<OneD, NekDouble>(totszedges1d)   ;
                Vdxy1xzm1   = Array<OneD, NekDouble>(totszedges1d);
                Vdzy1xzm1 = Array<OneD, NekDouble>  (totszedges1d) ;
                edgeptsin[0] = edgexyztemp;
                edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
                E->PhysEvalBasisGrad(edgeptsin,Vy1xzm1, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1 );  

                // edge left bot (z = -1), (x = -1)
                Vxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
                Vdyxm1yzm1 = Array<OneD, NekDouble>(totszedges1d)   ;
                Vdxxm1yzm1 = Array<OneD, NekDouble>(totszedges1d)   ;
                Vdzxm1yzm1 = Array<OneD, NekDouble>(totszedges1d)   ;
                
                edgeptsin[1] = edgexyztemp;
                edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
                edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
                E->PhysEvalBasisGrad(edgeptsin,Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);  
                
                //edge left top (x = -1), (z = 1))
                Vxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
                Vdyxm1yz1 = Array<OneD, NekDouble>  (totszedges1d) ;
                Vdxxm1yz1 = Array<OneD, NekDouble >(totszedges1d)  ;
                Vdzxm1yz1 = Array<OneD, NekDouble>(totszedges1d)   ;

                edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);           
                E->PhysEvalBasisGrad(edgeptsin,Vxm1yz1,Vdxxm1yz1, Vdyxm1yz1, Vdzxm1yz1);  
                
                //edge right bot ( z = -1) (x = 1)
                Vx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
                Vdyx1yzm1 =Array<OneD, NekDouble>(totszedges1d);
                
                Vdxx1yzm1= Array<OneD, NekDouble>(totszedges1d);
                
                Vdzx1yzm1  = Array<OneD, NekDouble>(totszedges1d);
                

                edgeptsin[1] = edgexyztemp;
                edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);           
                edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);           
                E->PhysEvalBasisGrad(edgeptsin,Vx1yzm1, Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);   
            }
            if(numsurfaces == 5) // prism
            {
                surfcoords = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(numsurfaces);
                //we need convolution of qzin + qzinmid as input points
                int totszsurf2d = (testcoord2dquad[0].size())*E->GetNcoeffs();
              
                Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);
                for(int k = 0; k < dimension-1; k++)
                {
                    surfptsin[k] = Array<OneD, NekDouble>(testcoord2dquad[k]);
                    surfptsintemp[k] = Array<OneD, NekDouble>(testcoord2dquad[k]);
                }

                surfptsin[dimension-1] = Array<OneD, NekDouble>(testcoord2dquad[0].size()); 
                surfptsintemp[dimension-1] = Array<OneD, NekDouble>(testcoord2dquad[0].size()); 

                //surface bot z = -1
                Vxyzm1 =  Array<OneD, NekDouble>(totszsurf2d);
                surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
                surfcoords[0] = surfptsintemp;

                E->PhysEvalBasisGrad(surfptsintemp,Vxyzm1, NullNekDouble1DArray,NullNekDouble1DArray,NullNekDouble1DArray);   


                //surface hypt x+z = 0
                Vxzmy =  Array<OneD, NekDouble>(totszsurf2d);
                Vmath::Smul(surfptsin[0].size(), -1.0, &surfptsin[0][0], 1, &surfptsintemp[2][0], 1);
                surfptsintemp[0] = surfptsin[0];
                surfptsintemp[1] = surfptsin[1];
                surfcoords[1] = surfptsintemp;

                E->PhysEvalBasisGrad(surfptsintemp,Vxzmy,NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   
                
                //surface left x = -1
                Vxm1yz =  Array<OneD, NekDouble>(totszsurf2d);
                surfptsintemp[0] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
                surfptsintemp[1] = surfptsin[0];
                surfptsintemp[2] = surfptsin[1];
                surfcoords[2] = surfptsintemp;
                
                E->PhysEvalBasisGrad(surfptsintemp,Vxm1yz, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   
                totszsurf2d = (testcoord2dtri[0].size())*E->GetNcoeffs();
              
                for(int k = 0; k < dimension-1; k++)
                {
                    surfptsin[k] = Array<OneD, NekDouble>(testcoord2dtri[k]);
                    surfptsintemp[k] = Array<OneD, NekDouble>(testcoord2dtri[k]);
                }

                surfptsin[dimension-1] = Array<OneD, NekDouble>(testcoord2dtri[0].size()); 
                surfptsintemp[dimension-1] = Array<OneD, NekDouble>(testcoord2dtri[0].size()); 

                //surface front y = -1
                Vxym1zproject =  Array<OneD, NekDouble>(Etri->GetBasis(0)->GetZ().size()*Etri->GetBasis(1)->GetZ().size()*E->GetNcoeffs());
                Vxym1z =  Array<OneD, NekDouble>(totszsurf2d);
                
                Array<OneD, Array<OneD, NekDouble> >surfproject(3);
                surfproject[0] = Array<OneD, NekDouble> (Etri->GetBasis(0)->GetZ().size()*Etri->GetBasis(1)->GetZ().size());//Etri->GetBasis(0)->GetZ();
                surfproject[1] = Array<OneD, NekDouble> (surfproject[0].size(), -1.0);
                surfproject[2] = Array<OneD, NekDouble> (surfproject[0].size());
                int ctrn = 0;
                for(int k = 0; k < Etri->GetBasis(0)->GetZ().size(); k++)
                {
                    for(int p = 0; p < Etri->GetBasis(1)->GetZ().size(); p++)
                    {
                        surfproject[0][ctrn] = (Etri->GetBasis(0)->GetZ())[k];
                        surfproject[2][ctrn] = (Etri->GetBasis(1)->GetZ())[p];
                        ctrn++;
                    }
                }
                E->PhysEvalBasisGrad(surfproject, Vxym1zproject, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

                surfptsintemp[1] = surfptsintemp[0];
                surfptsintemp[0] = surfptsin[0];
                surfptsintemp[2] = surfptsin[1];
                surfcoords[3] = surfptsintemp;
                
                E->PhysEvalBasisGrad(surfptsintemp,Vxym1z, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   
                
                //surface back y = 1
                Vxy1zproject =  Array<OneD, NekDouble>(Etri->GetBasis(0)->GetZ().size()*Etri->GetBasis(1)->GetZ().size()*E->GetNcoeffs());

                Vxy1z =  Array<OneD, NekDouble>(totszsurf2d);

                surfproject[1] = Array<OneD, NekDouble> (surfproject[0].size(), 1.0);
                E->PhysEvalBasisGrad(surfproject, Vxy1zproject, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   

                
                surfptsintemp[1] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
                surfcoords[4] = surfptsintemp;

                E->PhysEvalBasisGrad(surfptsintemp,Vxy1z, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
            }
        }

        if(Opt_needed(coeffs))
        {

            cout<<"\n need optimization";
            WriteCSV(coeffs);

            Do_optimize(coeffs);
            
            cout<<"\n doopt done\n verifying...\n";
            if(Opt_needed(coeffs,1 ))
            {
                cout<<"func="<<funcdef;

                cout<<"\n fail\n\n";
                exit(0);
            }
            else
            {
                cout<<"func="<<funcdef;
                cout<<"\n pass\n\n";
                WriteCSV(coeffs, 1);

                sol = phys;
            }
        }
        else{
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


void project_surfaces( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&retquad, Array<OneD, Array<OneD, NekDouble> >&rettri, int d)
{
    boost::ignore_unused(d);
    //    int modes = E->GetNcoeffs();
    //    int modesquad = Equad->GetNcoeffs();
    //int modestri = Etri->GetNcoeffs();
    //    if( d == 0)
    //{

        //bot surface  z = -1
        surfuhats(uhats, retquad[0], Vxyzm1, 0);

        //left surface x = -1
        surfuhats(uhats, retquad[1], Vxm1yz, 0);

        //surface hypt x+z = 0
        surfuhats(uhats, retquad[2], Vxzmy,  0);

        //front surface y = -1  
        surfuhats(uhats, rettri[0], Vxym1zproject,  1);
 
        //surface back y = 1  
        surfuhats(uhats, rettri[1], Vxy1zproject, 1);

    // }
    // else
    // {

    //     //bot surface  z = -1
    //     surfacepquhats(uhats, retquad[0], Vxyzm1,0);
        
    //     //left surface x = -1
    //     surfacepquhats(uhats, retquad[1], Vxm1yz, 0);
        
    //     //surface hypt x+z = 0
    //     surfacepquhats(uhats, retquad[2], Vxzmy, 0);
        
    //     //front surface y = -1  
    //     surfacepquhats(uhats, rettri[0], Vxym1z,  1);
        
    //     //surface back y = 1  
    //     surfacepquhats(uhats, rettri[1], Vxy1z, 1);

    // }
}


// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret , int d)
{
    //    int modes = E->GetNcoeffs();
    //    int modes1D = (3*(E->GetBasis(0)->GetNumModes()-1));//3*(E->GetBasis(0)->GetNumModes()-1);
    // In case of tri:
    // hypotenuse = 2, bot = 1, left = 2
    
    // assert d = 0 or 1
    if(d == 0)
    {   
        int modes =  E3seg->GetNcoeffs();//(3*(E->GetBasis(0)->GetNumModes()));

        for(int k= 0; k < numedges; k++)
        {
            ret[k] = Array<OneD, NekDouble>(modes);
        }
    }
    if(numedges == 9) // prism
    {

        if(d == 1)
        {

            // edge front left (x = -1) (y = -1)
            edgederpquhats(uhats, ret[0], E->GetNcoeffs(),  Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);
                        
            //edge front hypt (y = -1) (z + x = 0)
            edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vym1xmz, Vdxym1xmz, Vdyym1xmz, Vdzym1xmz);
            
            //edge front bot (y = -1) (z = -1)
            edgederpquhats(uhats, ret[2], E->GetNcoeffs(),  Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);
     
            // edge back left (y = 1), (x = -1)
            edgederpquhats(uhats, ret[3], E->GetNcoeffs(), Vxm1y1z, Vdxxm1y1z, Vdyxm1y1z, Vdzxm1y1z);
     
            //edge back hypt ( y = 1) (z +x = 0)
            edgederpquhats(uhats, ret[4], E->GetNcoeffs(),  Vy1xmz, Vdxy1xmz, Vdyy1xmz, Vdzy1xmz);
   
            //edge back bot (y = 1) (z = -1))
            edgederpquhats(uhats, ret[5], E->GetNcoeffs(),  Vy1xzm1, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);
            
            // edge left bot (z = -1), (x = -1)
            edgederpquhats(uhats, ret[6], E->GetNcoeffs(), Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);
            
            //edge left top (x = -1), (z = 1))
            edgederpquhats(uhats, ret[7], E->GetNcoeffs(),   Vxm1yz1,Vdxxm1yz1, Vdyxm1yz1, Vdzxm1yz1);
            
            //edge right bot ( z = -1) (x = 1)
            edgederpquhats(uhats, ret[8], E->GetNcoeffs(),   Vx1yzm1, Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);
        
     
        }
        else if( d == 0)
        {
            // edge front left (x = -1) (y = -1)
            deruhatsedges(uhats, ret[0], Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);
        
            //edge front hypt (y = -1) (z + x = 0)
            deruhatsedges(uhats, ret[1],Vym1xmz, Vdxym1xmz, Vdyym1xmz, Vdzym1xmz);
            
            //edge front bot (y = -1) (z = -1)
            deruhatsedges(uhats, ret[2],  Vym1xzm1,  Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);
            
            // edge back left (y = 1), (x = -1)
            deruhatsedges(uhats, ret[3],Vxm1y1z, Vdxxm1y1z, Vdyxm1y1z, Vdzxm1y1z);
            //edge back hypt ( y = 1) (z +x = 0)
            deruhatsedges(uhats, ret[4], Vy1xmz, Vdxy1xmz, Vdyy1xmz, Vdzy1xmz);
                   
            //edge back bot (y = 1) (z = -1))
            deruhatsedges(uhats, ret[5], Vy1xzm1, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);
            
            // edge left bot (z = -1), (x = -1)
            deruhatsedges(uhats, ret[6],Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);
            
            //edge left top (x = -1), (z = 1))
            deruhatsedges(uhats, ret[7],  Vxm1yz1, Vdxxm1yz1, Vdyxm1yz1, Vdzxm1yz1);

            //edge right bot ( z = -1) (x = 1)
            deruhatsedges(uhats, ret[8],  Vx1yzm1,  Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);


        }
    }
   
}

FILE *csvfile(int flag)
{
    static char name[LOGNAME_SIZE];
    char fullname[] = "Pri_";//[LOGNAME_SIZE+4];
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
    int uhatstot = uhats.size();
    int totpts = testcoord3d[0].size();
    Array<OneD, NekDouble> temp(uhatstot);
    for(int i = 0; i<totpts; i++)
    {
        Vmath::Vmul(uhatstot, &interioreval3d[i], totpts, &uhats[0], 1, &temp[0], 1);
        NekDouble v = Vmath::Vsum(uhatstot, temp, 1);
        string s = std::to_string(testcoord3d[0][i])+","+std::to_string(testcoord3d[1][i])+","+std::to_string(testcoord3d[2][i])+","+std::to_string(v)+"\n";
        const char* c = s.c_str();
        fprintf (file, c); 
    }
    fclose(file);
}


void surfuhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz, int surfid )
{
    int modes = uhats.size();

    Array<OneD, NekDouble> temp(uhats.size());  

    if(surfid == 1)
    {
        int totpts = Etri->GetTotPoints();
        Array<OneD, NekDouble> vals(totpts);    
        
        //Vxyz*uhats -> project to -> E3quad
        for(int k = 0; k < totpts; k++)
        {    
            Vmath::Vmul(uhats.size(), &Vxyz[k], totpts, &uhats[0], 1, &temp[0], 1);
            vals[k]  = Vmath::Vsum(modes, temp, 1);
        }
        
        Etri->FwdTrans(vals, ret);

    }
    else
    {
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
    boost::ignore_unused(modes);
    int uhatstot = uhats.size();
    int totpts = edgeptsin[0].size();//E3seg->GetTotPoints();
    Array<OneD, NekDouble> temp(totpts), temp2(uhatstot), temp3(uhatstot);
    Array<OneD, NekDouble> pqeval(totpts);
    NekDouble v1, v2;
    //for(int d = 0; d < dimension; d++)
    //{
    for(int i = 0; i<totpts; i++)
    {
        Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd0[0]+i, totpts, &temp2[0], 1);
        v1  = Vmath::Vsum(uhatstot, temp2, 1); 
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
            Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd1[0]+i, totpts, &temp2[0], 1);
            v1  = Vmath::Vsum(uhatstot, temp2, 1); 
            Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);

            v2  = Vmath::Vsum(uhatstot, temp3, 1);  

            v1 = v2*v1;

            // At this point,                 

            // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

            // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
            Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &temp2[0], 1);

 
        
            v2  = Vmath::Vsum(uhatstot, temp2, 1);  
 
            Vmath::Vmul(uhats.size(), &Vxyd1[i], totpts, &uhats[0], 1, &temp2[0], 1);

            v1=  v2*Vmath::Vsum(uhats.size(), temp2, 1) - v1;
 
            pqeval[i] += v1;
 
        }
        if(dimension == 3)
        {
            Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxyd2[0]+i, totpts, &temp2[0], 1);
            v1  = Vmath::Vsum(uhatstot, temp2, 1); 
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

void surfacepquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret, Array<OneD, NekDouble> Vxy, int surfid)
{
    if(surfid == 0 ) 
    {
        Array<OneD, NekDouble> vals(Equad->GetTotPoints());
        // calculate pq() on surface
        pq(uhats, demo.GetCoords(Equad), Vxy, NullNekDouble1DArray, vals);
        Equad->FwdTrans(vals, ret);
    }
    else
    {
        Array<OneD, NekDouble> vals(Etri->GetTotPoints());
        // calculate pq() on surface
        pq(uhats, demo.GetCoords(Etri), Vxy, NullNekDouble1DArray, vals);
        Equad->FwdTrans(vals, ret);
        
    }
}



// sig = 0 -> called by opt_needed
// sig = 1 -> called by sphere_rotation
Array<OneD, Array< OneD,  NekDouble> > call_find_roots(Array<OneD,  NekDouble> &uhats , int d, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhatsquad,  Array<OneD, Array<OneD, NekDouble> >&surfaceuhatstri,  int sig)
{

    boost::ignore_unused(d);
    int dimension = E->GetShapeDimension(); 
    vector<vector< NekDouble> > ret(dimension);

    Array<OneD, Array<OneD, NekDouble> > retarr;//(dimension);

    if(numedges == 9) //prism
    {

        for(int ii = 0; ii < numedges; ii++)
        {
    
            vector<vector<NekDouble> > tmp(dimension);
 
            tmp  = (find_roots(uhatsedges[ii],  NullNekDouble1DArray,  0, 0, sig)) ;
            for(int p = 0; p < tmp[0].size(); p++)
            {
                if(ii == 0) // edge front left (x = -1) (y = -1)
                {
                    ret[0].push_back(-1);
                    ret[1].push_back(-1);
                    ret[2].push_back(tmp[0][p]);
                        
                }
                else if(ii == 1) //edge front hypt (y = -1) (z = -x)
                {
                    ret[0].push_back(tmp[0][p]);
                    ret[2].push_back(-tmp[0][p]);
                    ret[1].push_back(-1);
			
                    ret[0].push_back(-tmp[0][p]);
                    ret[2].push_back(tmp[0][p]);
                    ret[1].push_back(-1);
                }
                else if(ii == 2) //edge front bot (y = -1) (z = -1)

                {
                    ret[0].push_back(tmp[0][p]);
                    ret[1].push_back(-1);
                    ret[2].push_back(-1);
                    
                    
                }
                else if(ii == 3) //edge back left (y = 1), (x = -1)
                {
                    ret[0].push_back(-1);
                    ret[1].push_back(1);
                    ret[2].push_back(tmp[0][p]);
                }
                else if(ii == 4) //edge backhypt ( y = 1) (z = -x)
                {
                    ret[1].push_back(1);
                    ret[0].push_back(tmp[0][p]);
                    ret[2].push_back(-tmp[0][p]);

                    ret[1].push_back(1);
                    ret[0].push_back(-tmp[0][p]);
                    ret[2].push_back(tmp[0][p]);
                }
                    
                else if(ii == 5) //edge back bot (y = 1) (z = -1)
                {
                    ret[1].push_back(1);
                    ret[0].push_back(tmp[0][p]);
                    ret[2].push_back(-1);
                }
                else if(ii == 6)//edge left bot (z = -1), (x = -1)  
                {
                    ret[0].push_back(-1);
                    ret[1].push_back(tmp[0][p]);
                    ret[2].push_back(-1);
                }
                else if(ii == 7) // edge left top (x = -1), (z = 1))
                {
                    ret[2].push_back(1);
                    ret[1].push_back(tmp[0][p]);
                    ret[0].push_back(-1);
                }
                else if(ii == 8) //edge right bot ( z = -1) (x = 1)
                {
                    ret[2].push_back(-1);
                    ret[1].push_back(tmp[0][p]);
                    ret[0].push_back(1);
                }
                    
            }
            

            
        }
        // ret[0].push_back(1);
        // ret[1].push_back(-7.83575e-06);
        // ret[2].push_back(-1); 
    }
    if(numsurfaces == 5)
    {   
        // call 2D rootfinder on each surface:
        vector<vector<NekDouble> > tmp(dimension);
        for(int ii = 0; ii < numsurfaces; ii++)
        {

                
            if(ii == 0)     //bot surface  z = -1 
            {
                // surfid = 0 for quad;
                tmp  = (find_roots(surfaceuhatsquad[ii], Vxyzm1, 0, 1, sig, 1, 0)) ;  
            

                for(int p = 0; p < tmp[0].size(); p++)
                {

                    ret[2].push_back(-1);
                    ret[1].push_back(tmp[1][p]);
                    ret[0].push_back(tmp[0][p]);
                }
            }
                    
            else if(ii == 1) //surf left x = -1
            {

                // surfid = 0 for quad;
                tmp  = (find_roots(surfaceuhatsquad[ii], Vxzmy,   0, 1, sig, 1, 0)) ;  
                for(int p = 0; p < tmp[0].size(); p++)
                {
                    
                    ret[0].push_back(-1);
                    ret[2].push_back(tmp[1][p]);
                    ret[1].push_back(tmp[0][p]);
                }
            }
            else if(ii == 2) //surface hypt x+z = 0
            {

                // surfid = 0 for quad;
                tmp  = (find_roots(surfaceuhatsquad[ii], Vxm1yz,   0, 1, sig, 1, 0)) ;

                for(int p = 0; p < tmp[0].size(); p++)
                {  

                    ret[2].push_back(-tmp[0][p]);
                    ret[1].push_back(tmp[1][p]);
                    ret[0].push_back(tmp[0][p]);
                    ret[2].push_back(tmp[0][p]);
                    ret[1].push_back(tmp[1][p]);
                    ret[0].push_back(-tmp[0][p]);

                }
            }
            else if(ii == 3) //surf front y = -1
            {
                    
                // surfid = 1 for tri;
                tmp  = (find_roots(surfaceuhatstri[0], Vxym1zproject, 0, 1, sig, 1, 1)) ;  

                for(int p = 0; p < tmp[0].size(); p++)
                {
                    ret[1].push_back(-1);
                    ret[2].push_back(tmp[1][p]);
                    ret[0].push_back(tmp[0][p]);
                
                }
            }
            else if(ii == 4) //surf back y = 1
            {

                // surfid = 0 fortri;
                tmp  = (find_roots(surfaceuhatstri[1], Vxy1zproject, 0, 1, sig, 1, 1)) ;  
                for(int p = 0; p < tmp[0].size(); p++)
                {
                    ret[1].push_back(1);
                    ret[2].push_back(tmp[1][p]);
                    ret[0].push_back(tmp[0][p]);
                }
                    
            }
        }
    }
            
    //find interior roots:
    vector<vector<NekDouble> > tmp;

    tmp =  (find_roots(uhats, interioreval3d, 0, 1, sig)) ;

    if(tmp[0].size() > 0)
    {
        for(int p = 0; p <dimension; p++)
        {
            ret[p].push_back(tmp[p][0]);
        }
    }


    retarr = vectoarr(ret);

    return retarr;
}


// Called by call_find_roots depending on edges or interior region 
// rootfinding
// d=1 :: derivatives
// flag = 0, edges (1D rootfinding using confederate matrix)
// flag = 1, interior region rootfinding using gradient descent
// sig = 0 -> called by opt_needed
// sig = 1 -> called by sphere_rotation
// surfflag = 0 -> only minima in 3D volume in gradient descent
// surffflag = 1 -> minima in the area of 2D surfaces using grad descent
// surfid = 1 is triangle
// surfid = 0 is quad
vector<vector<NekDouble> > find_roots(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> evalxyz, int d, int flag ,int sig, int surfflag, int surfid)
{

    vector<vector< NekDouble> > ret(dimension);
    //Confederate matrix approach
    if(dimension == 1 || flag == 0)
    {
        int N = uhats.size();

        vector<NekDouble> uhatsmon;
        while(abs(uhats[N-1])<1e-9)
        {
            N = N-1;
            
            if(N == 0)
            {
        
                ret[0].push_back(-1.0);
                ret[1].push_back(1.0);
                
                return ret;
            }

        }

        vector<NekDouble> temp(N);
        // convert uhats to monomial, find roots of uhats or der of uhats
        for(int k = 0; k < N; k++)
        {
            for(int jj= 0 ; jj < N; jj++)
            {
                temp[jj ] = C[jj][k];
            }
            Vmath::Vmul(N, &temp[0], 1,  &uhats[0], 1, &temp[0], 1);
            NekDouble temp2 = Vmath::Vsum(N, &temp[0], 1);
            uhatsmon.push_back(temp2);
        }

        if(abs(Vmath::Vmax(uhatsmon.size(), &uhatsmon[0], 1))<1e-10)
        {
            ret[0].push_back(-1.0);
            ret[1].push_back(1.0);
            
            return ret;
         
        }

        // truncate trailing zeros
        while(abs(uhatsmon[N-1])<1e-8)
        {
            N = N-1;
            if(N == 0)
            {
        
                ret[0].push_back(-1.0);
                ret[1].push_back(1.0);
                
                return ret;
            }

        }
        //N = uhatsmon.size();

        if(N == 0 || Vmath::Vmin(N,&uhatsmon[0],1)<-1e10)
        {
            ret[0].push_back(-1.0);
            ret[1].push_back(1.0);
                
            return ret;
        }

      	
        // now size of uhatsmon = N;
        vector<NekDouble> uhatsdiff;
        // if d == 1,

        if(d == 1)
        {
            for(int k = 1; k < N; k++)
            {
                uhatsdiff.push_back(k*uhatsmon[k]);
            }
            N = N-1;
            if(N == 0)
            {
                ret[0].push_back(-1.0);
                ret[1].push_back(1.0);
                
                return ret;
            }
        
        }
        else //d == 0
        {
            for(int k = 0; k<N; k++)
                uhatsdiff.push_back(uhatsmon[k]);
        }
        
        if(N == 1)
        {
            ret[0].push_back(-1.0);
            ret[1].push_back(1.0);
                
            return ret;
        }
        
        Vmath::Smul(N, 1.0/uhatsdiff[N-1], &uhatsdiff[0], 1, &uhatsdiff[0], 1);


        vector<NekDouble> EIG_R = demo.FindEigenval(uhatsdiff, N);


        for(int kk = 0; kk <EIG_R.size(); kk++)
        {   
            ret[0].push_back( EIG_R[kk] );
            
        }


    }
    else if(dimension > 1 && flag == 1 )
    {    
        ret = gradient_descent(uhats, evalxyz,  sig, surfflag, surfid);
       
    }
    


    return ret;    
    
}

// sig = 1 -> from do_opt
// sig = 0 -> from opt_needed
// surfflag = 0 -> volume
// surfflag = 1 -> surfaces
// surfid = 0 -> quad
// surfid = 1 -> tri
vector< vector< NekDouble> >  gradient_descent(Array<OneD, NekDouble> uhats, Array<OneD, NekDouble> evalxyz,  int sig, int surfflag, int surfid  )
{        
    boost::ignore_unused(surfid, evalxyz);
    // Gradient descent
    int temp_dim = dimension;
    if(surfflag == 1)
        temp_dim = dimension - 1;
    // if(surfflag == 1)
    //     temp_dim = dimension - 1;
    // Assert that d = 0
       
    Array<OneD, NekDouble> g(temp_dim);
    double inf = numeric_limits<double>::infinity();
            

    Array<OneD,NekDouble> xnew(temp_dim) ,x(temp_dim);
    NekDouble gprev = inf; // use inf
    int idxgprev;
    Array<OneD, NekDouble> tmp(temp_dim);
    Array<OneD, Array<OneD, NekDouble> > xarr(temp_dim);
    vector< vector< NekDouble> > ret(temp_dim);
    Array<OneD, Array<OneD, NekDouble> > tc(temp_dim);
    StdExpansion *tempE;
    Array<OneD, NekDouble> eval;
    int sz;
    tempE = E;
    Array<OneD, Array<OneD, NekDouble > > tempeval(4);

    if( (dimension == 3 && surfflag == 1) )
    {
        if(surfid == 0)
        {
            tc = testcoord2dquad;
            sz = tc[0].size();
            eval = interioreval2dquad;
            tempE = Equad;
            

        }
        else
        {
            tc = testcoord2dtri;
            sz = tc[0].size();
            eval = interioreval2dtri;
            tempE = Etri;

            // for(int m = 0; m < testcoord2dtri[0].size(); m++)
            // {
            //     if(testcoord2dtri[0][m] + testcoord2dtri[1][m] > 0)
            //     {
            //         exit(0);
            //     }
            // }
            
        }
        // //assert surfid is in 0-4
        // tc = surfcoords[surfid];//testcoord2dquad;
        // sz = tc[0].size();
        // eval = evalxyz;
        
    }
    else //dim = 3 and surfflag is off
    {
        tc =  testcoord3d;
        sz = tc[0].size();

        eval = interioreval3d;// testcoord2d;
    
    }
    for(int i = 0; i < 4; i++)
    {
        tempeval[i] = Array<OneD, NekDouble>(tempE->GetNcoeffs());  
    }
    
    Array<OneD, NekDouble> nullarr(0);
    Array<OneD, NekDouble> temp2(sz);

    if(sig==0)
    {

        Array<OneD, NekDouble> temp(uhats.size());
 
        for(int k = 0; k < sz; k++)
        {

            Vmath::Vmul(uhats.size(), &eval[k], sz, &uhats[0], 1, &temp[0], 1);          
            temp2[k] = Vmath::Vsum(temp.size(), temp, 1);
        }

        gprev = Vmath::Vmin(temp2.size(), temp2, 1);	
        idxgprev = Vmath::Imin(sz, temp2, 1);        
            
    }
    else
    {
        pq(uhats, tc, eval, nullarr, temp2 );
        gprev = Vmath::Vmin(sz, temp2, 1);
        idxgprev = Vmath::Imin(sz, temp2, 1);
    }

    xnew[0] = tc[0][idxgprev];
    xnew[1] = tc[1][idxgprev];
    
    if(surfflag == 0)
    {
        xnew[2] = tc[2][idxgprev]; 
    }

    if(gprev < 0 && abs(gprev)>1e-11)
    {
        int n = 0;
        NekDouble epsl = 1.0;

        iterGD = 1e-3;
        secarg = 1e3;
        NekDouble stepsize = iterGD, iter = secarg;
        Array<OneD, Array<OneD, NekDouble> > xastaa(temp_dim);

        Array<OneD, NekDouble> dereval(temp_dim);
        Array<OneD, NekDouble> xsave(xnew.size());
        for(int k = 0; k < temp_dim; k++)
            xsave[k] = xnew[k];
        // loop for gradient descent 2D:
        
        bool res, res2, res3 = true;
        
        Array<OneD, NekDouble> temp1 (uhats.size());
        while(n < iter)
        {
            epsl = 0.0;
            for(int p = 0; p < temp_dim; p++)
            {

                xastaa[p] = Array<OneD, NekDouble>(1,xnew[p]);
            }   
           
            tempE->PhysEvalBasisGrad(xastaa, tempeval[0], tempeval[1], tempeval[2], tempeval[3]);
            // multiply tempeval by uhats to get g[p]
           
            if(sig == 0) // from opt_needed
            {
                
                for(int p = 0; p < temp_dim; p++)
                {
                 
                    Vmath::Vmul(uhats.size(), &tempeval[p+1][0], 1, &uhats[0], 1, &temp1[0], 1 );
                    g[p] = Vmath::Vsum(temp1.size(), temp1, 1);
    
                    epsl += abs(g[p]);
                }
            }
            else
            {
                for(int p = 0; p < temp_dim; p++)
                {
                    if(surfid == 1)

                    derpq(uhats, dereval[p], tempeval[0], tempeval[p+1]);
                 
                    g[p] = dereval[p];
                    epsl += abs(g[p]);
                }
            }
            x[0] = xnew[0] - stepsize*(g[0]);
            x[1] = xnew[1] - stepsize*(g[1]);
            if(surfflag == 0)
            {
                x[2] = xnew[2] - stepsize*(g[2]);
            }
            res2 = true;
            for(int k = 0; k < temp_dim; k++)
            {   
                               if(abs(x[k])-abs(xnew[k]) > 1e-10)
                {
                    res2 = res2 && false  ;
                    break;
                }
            }
            res3 = true;
            res = true;
            if(surfflag == 1 &&surfid == 1) // or if surf is hypt
            {
                res = (x[0]+x[1])<=0;
                //cout<<"x[0]+x[1]="<<x[0]+x[1]<<"\n";
            }
            if(surfflag == 0)
            {
                //in volume case, check if x+z <= 0
                res3 = (x[0]+x[2])<=0 ;   
            }            
            if( ((std::isnan(epsl)) || (epsl <1e-10) || (epsl>1e3) || (res2)) && (res) && (res3) && abs(x[0]) <=1.0 && abs(x[1])<=1.0  )
	    {
                if(surfflag == 0)
                {
                    if(abs(x[2]) <=1.0) 
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }
            for(int p = 0; p < temp_dim; p++)
	    {
                xnew[p] = x[p];
            }            
            
            n = n+1;
            
            
        }
        res = true;
        if(surfflag == 1 && surfid == 1)
        {
            res = (x[0]+x[1])<=0;
        }
        res3 = true;

        if(surfflag == 0)
        {
            //check if x+z <= 0
            res3 = (x[0]+x[2])<=0;   
        }            
        for(int k = 0; k < temp_dim; k++)
        {
            if(abs(x[k]) >1.0)
            {
                res2 = false;
                break;
            }
        }   
        //check
        if((res2 && res && res3 && n < iter &&(!std::isnan(epsl))))
        {
            for(int p = 0; p < temp_dim; p++)
            {
                ret[p].push_back( (x[p]));
            }
        }
    
        else
        {
            for(int p = 0; p < temp_dim; p++)
            {
                
                ret[p].push_back(xsave[p]);
            }
            
        }

    
    }

    return ret;
}
    
    
    
void derpq(Array<OneD, NekDouble> &uhats, NekDouble &ret,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd)
{
    int N = uhats.size();
    Array<OneD, NekDouble> temp(N);

    //\sum(phi^2)

    Vmath::Vmul(N, Vxy, 1, Vxy, 1, temp, 1); 
    NekDouble tmp = Vmath::Vsum(N, temp, 1);

    //\sum(u*phi')
    Vmath::Vmul(N, uhats, 1,Vxyd, 1, temp, 1);
    NekDouble v1 = Vmath::Vsum(N, temp, 1);
  
    v1 = v1*( Vmath::Vsum(N, temp, 1))*pow(tmp, -0.5);
    //  NekDouble v2 = pow(tmp , -1.5);
    //\sum(phi*phi')
    Vmath::Vmul(N, Vxy, 1,Vxyd, 1, temp, 1);
    NekDouble v2  = Vmath::Vsum(N, temp, 1);

    //\sum(u*phi)
    Vmath::Vmul(N, uhats, 1, Vxy, 1, temp, 1); 
    v2 = v2*Vmath::Vsum(N, temp, 1);
    //\sum(phi*phi)^(-3/2)
    //  Vmath::Vmul(N, Vxy, 1,Vxyd1, 1, temp, 1);
    v2 = v2*( pow(tmp,-1.5));

    ret = v1 - v2;
} 

int Opt_needed(Array<OneD, NekDouble> uhats, int flag )
{
    int totModes = uhats.size();
    Array<OneD,  Array<OneD, NekDouble> > rootsarr;
    Array<OneD, Array<OneD, NekDouble> > edgeuhats  (numedges);
    Array<OneD, Array<OneD, NekDouble> > surfacederuhatsquad  (3);
    Array<OneD, Array<OneD, NekDouble> > surfacederuhatstri  (2);

    for(int t = 0; t < 3; t++)
    {
        surfacederuhatsquad[t] = Array<OneD, NekDouble>(Equad->GetNcoeffs());
    } 
    for(int t = 0; t < 2; t++)
    {
        surfacederuhatstri[t] = Array<OneD, NekDouble>(Etri->GetNcoeffs());
    } 

    project_surfaces(uhats, surfacederuhatsquad, surfacederuhatstri);
            
    project_edges(uhats, edgeuhats);

    //find roots of der of uhats, that is why d = 0


    rootsarr = call_find_roots(uhats, 0, edgeuhats, surfacederuhatsquad, surfacederuhatstri, 0);

    // evaluate ortho basis at roots
    // evalBasisRoots is flattened basis eval matri
        
      

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
        if(tempv < minv )
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

void pq(
        Array<OneD,NekDouble> uhats,
        Array<OneD, Array<OneD,  NekDouble> > roots,
        Array<OneD, NekDouble> V1,
        Array<OneD,NekDouble> &pqevalxast,
        Array<OneD,NekDouble> &fvals
        )
{
    //    boost::ignore_unused(tempE);
    int N = uhats.size();
    Array<OneD,NekDouble> w1(N);
    // Array<OneD, NekDouble> V1(roots[0].size()*N);
  
    vector<NekDouble> Vsumsq;
    
    // tempE->PhysEvalBasisGrad(roots, V1, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
    
    // V1 is flattened basis eval matrix
    
    // V2 = mat_mat_dot(V1,V1);
    // and
    // vector<double> Vsum = mat_sum_dim(V2,2);
 
    for(int i = 0; i < roots[0].size(); i++)
    {

        Vmath::Vmul(N, &V1[i], roots[0].size(), &uhats[0], 1, &w1[0], 1);
        fvals[i] = ( Vmath::Vsum(N, &w1[0], 1));
    
    }
    for( int i = 0; i < roots[0].size(); i++)
    {
        Vmath::Vmul(N, &V1[i], roots[0].size(), &V1[i], roots[0].size(), &w1[0], 1);

        Vsumsq.push_back(pow(Vmath::Vsum(N, &w1[0], 1),-0.5));
    }


    //    cts = vec_vec_dot(Vsumsq, ret);
    Vmath::Vmul(roots[0].size(), &Vsumsq[0], 1, &fvals[0], 1, &fvals[0], 1);
    // cout<<"\nfvals:  ";
    // for(int p = 0; p < fvals.size(); p++)
    //   cout<<fvals[p]<<" ";
    
    if(pqevalxast.size() > 0)
    {
                   
        int minidx = Vmath::Imin(fvals.size(), &fvals[0], 1);
        pqevalxast[0] = fvals[minidx];
        for(int k = 0; k < pqevalxast.size()-1; k++)
        {
            pqevalxast[k+1] = roots[k][minidx];
        }
    }

                   
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
    
    Array<OneD, NekDouble>  pqvalcoords(dim+1), xastarr(dim), utemp(N1), xastarrhist(dim), wsp1;   
    

    Array<OneD, Array<OneD, NekDouble> > surfaceuhatsquad  (3);
    Array<OneD, Array<OneD, NekDouble> > surfaceuhatstri  (2);
    Array<OneD, Array<OneD, NekDouble> > Pf(numedges);
    for(int k= 0; k < numedges; k++)
    {
            
        Pf[k] = Array<OneD, NekDouble>(E3seg->GetNcoeffs());//(3*(E->GetBasis(0)->GetNumModes()-1));        
    }
    
    NekDouble pqval;
    while (counter <= niter)
    {
        pqval = inf;
        utemp = d.back();

        //pq = @(xx) cfun(xx) .* (constraints(1,jj).*Vc(xx)*uhats);  
        project_edges(utemp, Pf, 1);
        
        // if(dimension > 2)
        // {
        //     //find surface roots: (surfaces = 12 is hex)
        for(int k = 0; k < 3; k++)
        {
            surfaceuhatsquad[k] = Array<OneD, NekDouble>(Equad->GetNcoeffs());
        }
        for(int k = 0; k < 2; k++)
        {
            surfaceuhatstri[k] = Array<OneD, NekDouble>(Etri->GetNcoeffs());
        }
        
        project_surfaces(utemp, surfaceuhatsquad, surfaceuhatstri, 1);

        // }
        optima = (call_find_roots(utemp, 0, Pf, surfaceuhatsquad, surfaceuhatstri,  1));
	// cout<<" roots:  \n";
	// for(int k = 0; k < optima.size(); k++)
	//   {
	//     for(int p = 0; p <optima[0].size(); p++)
	//       {
	// 	cout<<optima[k][p]<<" ";
	//       }
	//     cout<<"\n";
	    
	//   }
	// cout<<"\n";
        wsp1 = Array<OneD, NekDouble>(optima[0].size());
        
        // Vtmp is evaluation of basis fn at optima
        Array<OneD, NekDouble> Vtmp(optima[0].size() * E->GetNcoeffs());
        E->PhysEvalBasisGrad(optima, Vtmp, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        pq(utemp, optima, Vtmp,  pqvalcoords, wsp1);


        if (pqvalcoords[0] < pqval)
        {
            for(int k = 0; k  < dimension; k++)
            {
                xastarr[k] = pqvalcoords[k+1];
            }
            pqval = pqvalcoords[0];
        }
        
        cout<<"\n at counter="<<counter<<" min val="<<pqval<<" xast ="<<xastarr[0]<<" "<<xastarr[1]<<" "<<xastarr[2];//<< "diff ="<<xastarrhist[0] - xastarr[0]<<" "<<xastarrhist[1] - xastarr[1] <<" "<<xastarrhist[2] - xastarr[2];
        xastarrhist[0] = xastarr[0];   
        xastarrhist[1] = xastarr[1];        
        xastarrhist[2] = xastarr[2];
        // If minimum is non-negative, we're done
        if (pqval >= -tols.at(0))
        {
            break;
        }
        
        //vector<NekDouble> Vastsq;
        //        vector<NekDouble> Vast;
        Array<OneD, NekDouble> Vast(uhats.size()), Vastsq(uhats.size());
        Array<OneD, Array<OneD, NekDouble> > xastarr2(3);
        xastarr2[0] = Array<OneD, NekDouble>(1, xastarr[0]);
        xastarr2[1] = Array<OneD, NekDouble>(1, xastarr[1]);
        xastarr2[2] = Array<OneD, NekDouble>(1, xastarr[2]);
        NekDouble vastsqsum;

        // for( int ii = 0; ii < N1; ii++)
        // {
        //     Vast.push_back(E->PhysEvaluateBasis(xastarr, ii));
        //     Vastsq.push_back(Vast[ii]*Vast[ii]);
        
        // }
        E->PhysEvalBasisGrad(xastarr2, Vast, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        Vmath::Vmul(uhats.size(), Vast, 1, Vast, 1, Vastsq, 1);
        vastsqsum = Vmath::Vsum(N1, &Vastsq[0], 1);

        Array<OneD, NekDouble>  qast(N1);
        //cout<<"\n correction:";
        for(int i = 0; i<N1; i++)
        {
            qast[i] = ((1/sqrt(vastsqsum))*(Vast[i]));
            //  cout<<qast[i]<< " ";
        }
        Vmath::Smul(N1, pqval, &qast[0], 1, &qast[0], 1);
        
        Vmath::Vsub(utemp.size(), &utemp[0], 1, &qast[0], 1, &qast[0], 1);
        // cout<<"\n        uhats chng:";
        // for(int i = 0; i<utemp.size(); i++)
        // {
        //     cout<<utemp[i] - qast[i]<<" ";
        // }

        d.push_back(qast);

        counter = counter + 1;
    
    
    }
    cout<<"sphere_rotation took "<<counter<<"startv = "<<startval<<" at startcoords("<<startcoordx<<","<<startcoordy<<" ," <<startcoordz<<")   iterations\n  GD iter stepszie ="<<iterGD<<" maxiters = "<<secarg<<" o = "<<E->GetBasis(0)->GetNumModes()<<" p = "<<E->GetBasis(0)->GetNumPoints()<<" prism \n";;
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
