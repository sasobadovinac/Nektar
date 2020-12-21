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
// ermission is hereby granted, free of charge, to any person obtaining a
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


//declare Do_optimize
void Do_optimize(Array<OneD, NekDouble> &uhats);

Array<OneD, Array<OneD, NekDouble> > vectoarr(vector<vector< NekDouble> > vec);

// declare caller routine to find_roots
// flag = 0 -> opt_needed calls
// flag = 1 -> sphere_rot calls
Array<OneD, Array< OneD,  NekDouble> >  call_find_roots(Array<OneD, NekDouble> &uhatsall, int d = 0, Array< OneD, Array<OneD, NekDouble> >&uhatsedges =NullNekDoubleArrayofArray , Array< OneD, Array<OneD, NekDouble> >&surfaceuhats =NullNekDoubleArrayofArray ,  int flag = 0 );


//declare Opt_needed
int Opt_needed(Array<OneD, NekDouble> uhats, int flag = 0);


// for 2D elements to get uhats at edges
// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble> uhats,Array<OneD, Array<OneD, NekDouble> >&ret , int d = 0);

void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret, int modes, Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1 = NullNekDouble1DArray , Array<OneD, NekDouble> Vxyd2 = NullNekDouble1DArray);


void surfuhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz, int surfid = 0 );


void deruhatsedges(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble >&ret,  Array<OneD, NekDouble> Vxyz, Array<OneD, NekDouble> Vdxyz, Array<OneD, NekDouble> Vxdyz, Array<OneD, NekDouble> Vxydz);

void project_surfaces( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&retquad, Array<OneD, Array<OneD, NekDouble> >&rettri, int d = 0);

vector< vector<NekDouble> > gradient_descent(Array<OneD, NekDouble> uhats, Array<OneD, NekDouble> evalxyz, int sig, int surfflag , int surfid = 0);


FILE *csvfile(int flag);
void WriteCSV( Array<OneD, NekDouble>  uhats, int flag = 0);


string funcdef;
// Colleague matrix
int numedges, numsurfaces;
StdExpansion *E;
StdExpansion *E3seg;
StdExpansion *Equad;
StdExpansion *Etri;

//StdExpansion *E3quad;
Array<OneD, NekDouble> qZin; //inner point grid for 2D rootfinding
Array<OneD, Array<OneD, NekDouble> > edgeptsin;
Array<OneD, NekDouble> qZinmid; //inner mid point grid for 2D rootfinding

// for pyr edges root finding

// edge front left (x = -1) (y = -1) EA
Array<OneD, NekDouble> Vxm1ym1z   ;
Array<OneD, NekDouble> Vdyxm1ym1z   ;
Array<OneD, NekDouble> Vdxxm1ym1z   ;
Array<OneD, NekDouble> Vdzxm1ym1z   ;
               
//edge front hypt (y = -1) (z = -x) EB
Array<OneD, NekDouble> Vym1xmz   ;
Array<OneD, NekDouble> Vdxym1xmz   ;
Array<OneD, NekDouble> Vdyym1xmz   ;
Array<OneD, NekDouble> Vdzym1xmz   ;
               

//edge back hypt (z = -x  or x = -z) (x = y) EC
Array<OneD, NekDouble> Vxeyxmz   ;
Array<OneD, NekDouble> Vdxxeyxmz   ;
Array<OneD, NekDouble> Vdyxeyxmz   ;
Array<OneD, NekDouble> Vdzxeyxmz   ;

//edge back left (y = -z) (x = -1) ED
Array<OneD, NekDouble> Vx1ymz   ;
Array<OneD, NekDouble> Vdxx1ymz   ;
Array<OneD, NekDouble> Vdyx1ymz   ;
Array<OneD, NekDouble> Vdzx1ymz   ;


//edge front bot (y = -1) (z = -1) AB
Array<OneD, NekDouble> Vym1xzm1   ;
Array<OneD, NekDouble> Vdyym1xzm1   ;
Array<OneD, NekDouble> Vdxym1xzm1   ;
Array<OneD, NekDouble> Vdzym1xzm1   ;
               
               
//edge back bot (y = 1) (z = -1)) DC
Array<OneD, NekDouble> Vy1xzm1   ;
Array<OneD, NekDouble> Vdyy1xzm1   ;
Array<OneD, NekDouble> Vdxy1xzm1   ;
Array<OneD, NekDouble> Vdzy1xzm1   ;

// edge left bot (z = -1), (x = -1) AD
Array<OneD, NekDouble> Vxm1yzm1   ; 
Array<OneD, NekDouble> Vdyxm1yzm1   ;
Array<OneD, NekDouble> Vdxxm1yzm1   ;
Array<OneD, NekDouble> Vdzxm1yzm1   ;
                
//edge right bot ( z = -1) (x = 1) BC
Array<OneD, NekDouble> Vx1yzm1   ;
Array<OneD, NekDouble> Vdyx1yzm1   ;
Array<OneD, NekDouble> Vdxx1yzm1   ;
Array<OneD, NekDouble> Vdzx1yzm1   ;
 
//surface bot z = -1
Array<OneD, NekDouble>Vxyzm1 ;


//surface hypt x+z = 0
Array<OneD, NekDouble> Vxmzyproject;  


//surface left x = -1
Array<OneD, NekDouble> Vxm1yzproject;  


//surface front y = -1
Array<OneD, NekDouble> Vxym1zproject;  

//surface back y +z = 0 
Array<OneD, NekDouble> Vxymzproject;  


int dimension ;
NekDouble startval, startcoordx, startcoordy, startcoordz;
            
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
    LibUtilities::PointsType pointsTypetriB = LibUtilities:: ePolyEvenlySpaced;    
    LibUtilities::PointsType pointsTypetriA = LibUtilities::ePolyEvenlySpaced;    

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

    case LibUtilities::ePrism:
        
        numedges = 9;
        numsurfaces = 5;
        break;

    case LibUtilities::ePyramid:
        
        numedges = 8;
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
	  sol[i] =sin(x[i])*sin(y[i])*sin(z[i])-1;//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.5   ;///sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(x[i])*sin(y[i])*sin(z[i]);//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(x[i])*sin(y[i])*sin(z[i]);//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(x[i])*sin(y[i])*sin(z[i]);//sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);//sin(x[i])*sin(y[i])*sin(z[i])+0.1;//sin(x[i])*sin(y[i])*sin(z[i]);// +0.1;//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(x[i])*sin(y[i])*sin(z[i])+0.1;//sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);//sin(x[i])*sin(y[i])*sin(z[i]);//sin(x[i])*sin(y[i])*sin(z[i]) +0.1;//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);//sin(x[i])*sin(y[i])*sin(z[i]);//sin(x[i])*sin(y[i])*sin(z[i]) ;//sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);//sin(x[i])*sin(y[i])*sin(z[i]);//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2;//sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2;////sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1;
	  funcdef ="sin(x[i])*sin(y[i])*sin(z[i]) -1;";////sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.5";//sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])+0.9 ";;//sin(x[i])*sin(y[i])*sin(z[i]);";//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.5";//sin(x[i])*sin(y[i])*sin(z[i]);";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])";//"sin(x[i])*sin(y[i])*sin(z[i])+0.1";//"sin(x[i])*sin(y[i])*sin(z[i]) ";//+0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i])+0.1";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(x[i])*sin(y[i])*sin(z[i]) +0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])";//"sin(x[i])*sin(y[i])*sin(z[i]) ";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(x[i])*sin(y[i])*sin(z[i]);";//sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1;";//"sin(x[i])*sin(y[i])*sin(z[i];";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2 ";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.9";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])+0.5";//"sin(x[i])*sin(y[i])*sin(z[i])";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";//"sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i])-2";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i]);";// "sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])  -2";//"sin(3*x[i])*sin(3*y[i])*sin(3*z[i])+0.1";//"sin(M_PI*x[i])*sin(M_PI*y[i])*sin(M_PI*z[i])-2";// -  exp(-1) + 1.37";
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
        vector<LibUtilities::BasisKey> tmpBasisKey = demo.GetBasisKey();       
        vector<LibUtilities::PointsKey> tmpPtsKey = demo.GetPointsKey();       
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


            demo.testcoord2dquad = Array<OneD, Array<OneD, NekDouble> > (2); 
            demo.testcoord2dtri = Array<OneD, Array<OneD, NekDouble> > (2); 
            Array<OneD, Array<OneD, NekDouble> > testcoord2dtritemp(2);
               
            int ctsz2dquad = 0;
            int sztri = (Etri->GetBasis(0)->GetZ()).size()*(Etri->GetBasis(1)->GetZ().size())  ;
            int szmidtri = ((Etri->GetBasis(0)->GetZ()).size()-1)*((Etri->GetBasis(1)->GetZ()).size() - 1)  ;
            
            for(int ii = 0; ii < 2; ++ii)
            {
                demo.testcoord2dquad[ii] = Array<OneD, NekDouble>(sz+szmid);
                testcoord2dtritemp[ii] = Array<OneD, NekDouble>(sztri+szmidtri);
            }

            for(int ii = 0; ii < qZin.size(); ++ii)
            {
                
                for(int jj = 0; jj < qZin.size(); ++jj)
                {
                    demo.testcoord2dquad[0][ctsz2dquad] = qZin[ii];
                    demo.testcoord2dquad[1][ctsz2dquad] = qZin[jj];
                        ctsz2dquad++;
                }
                
            }     
            
            // add midpt grid
            for(int ii = 0; ii < qZinmid.size(); ++ii)
            {
                for(int jj = 0; jj < qZinmid.size(); ++jj)
                {
                
                    demo.testcoord2dquad[0][ctsz2dquad] = qZinmid[ii];
                    demo.testcoord2dquad[1][ctsz2dquad] = qZinmid[jj];

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
                demo.testcoord2dtri[ii] = Array<OneD, NekDouble>(p);
                Vmath::Vcopy(p, &temptri[ii][0], 1, &demo.testcoord2dtri[ii][0], 1);
                
            }

            LibUtilities::BasisKey bkeychebtmp(LibUtilities::eChebyshev, tmpBasisKey[0].GetNumModes(),  pkeycheb);
	   

            sz = pow((qZin.size()+qZinmid.size()),3);;

            Array<OneD, Array<OneD, NekDouble> >testcoord3dtemp  (dimension); 
            demo.testcoord3d = Array<OneD, Array<OneD, NekDouble> > (dimension); 
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
                        if(qZin[ii] + qZin[kk] <=0 &&qZin[jj] + qZin[kk] <=0) 
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
                        if(qZinmid[ii] + qZinmid[kk] <=0 && qZinmid[jj] + qZinmid[kk] <=0 ) 
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
                demo.testcoord3d[ii] = Array<OneD, NekDouble>(ctsz);
                Vmath::Vcopy(ctsz, &testcoord3dtemp[ii][0], 1, &demo.testcoord3d[ii][0], 1);

            }

            demo.interioreval2dquad = Array<OneD, NekDouble >(Equad->GetNcoeffs()*demo.testcoord2dquad[0].size());
            Equad->PhysEvalBasisGrad(demo.testcoord2dquad, demo.interioreval2dquad , NullNekDouble1DArray , NullNekDouble1DArray, NullNekDouble1DArray);
                        
            demo.interioreval2dtri = Array<OneD, NekDouble >(Etri->GetNcoeffs()*demo.testcoord2dtri[0].size());

            Etri->PhysEvalBasisGrad(demo.testcoord2dtri, demo.interioreval2dtri , NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray); 

            demo.interioreval3d = Array<OneD, NekDouble >(E->GetNcoeffs()*demo.testcoord3d[0].size());
                
            E->PhysEvalBasisGrad(demo.testcoord3d, demo.interioreval3d, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
            
            
            if(numedges == 8) //pyr
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

                // edge front left EA (x = -1) (y = -1)
                Vxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
                Vdxxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
                Vdyxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
                Vdzxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
                edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
                edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

                E->PhysEvalBasisGrad(edgeptsin, Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);  
                
                //edge front hypt EB (y = -1) (z + x = 0)
                Vym1xmz = Array<OneD, NekDouble>(totszedges1d);
                Vdxym1xmz = Array<OneD, NekDouble>(totszedges1d);
                Vdyym1xmz = Array<OneD, NekDouble>(totszedges1d);
                Vdzym1xmz = Array<OneD, NekDouble>(totszedges1d);
                Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[2][0], 1, &edgeptsin[0][0], 1);

                E->PhysEvalBasisGrad(edgeptsin,Vym1xmz, Vdxym1xmz,Vdyym1xmz,Vdzym1xmz);  

		//edge back hypt (z = -x  or x = -z) (x = y) EC
		Vxeyxmz     = Array<OneD, NekDouble>(totszedges1d);
		Vdxxeyxmz   = Array<OneD, NekDouble>(totszedges1d);
		Vdyxeyxmz   = Array<OneD, NekDouble>(totszedges1d);
		Vdzxeyxmz   = Array<OneD, NekDouble>(totszedges1d);
		edgeptsin[1] =  Array<OneD, NekDouble>(edgexyztemp);

		E->PhysEvalBasisGrad(edgeptsin,Vxeyxmz, Vdxxeyxmz, Vdyxeyxmz, Vdzxeyxmz);  

		//edge back left (y = -z) (x = -1) ED
	        Vx1ymz     = Array<OneD, NekDouble>(totszedges1d);;
		Vdxx1ymz   = Array<OneD, NekDouble>(totszedges1d);;
		Vdyx1ymz   = Array<OneD, NekDouble>(totszedges1d);
		Vdzx1ymz   = Array<OneD, NekDouble>(totszedges1d); ;
                Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[2][0], 1, &edgeptsin[1][0], 1);
                edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

		E->PhysEvalBasisGrad(edgeptsin,Vx1ymz, Vdxx1ymz, Vdyx1ymz, Vdzx1ymz);  
		
		
		//edge front bot (y = -1) (z = -1) AB
                Vym1xzm1 = Array<OneD, NekDouble>(totszedges1d);
                Vdxym1xzm1 = Array<OneD, NekDouble>(totszedges1d)   ;
                Vdyym1xzm1 = Array<OneD, NekDouble>(totszedges1d)    ;
                Vdzym1xzm1 =Array<OneD, NekDouble>(totszedges1d)    ;
                edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
                edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
                edgeptsin[0] = Array<OneD, NekDouble>(edgexyztemp);

                E->PhysEvalBasisGrad(edgeptsin,Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);  
                
		//edge back bot (y = 1) (z = -1)) DC
                Vy1xzm1 = Array<OneD, NekDouble>(totszedges1d);
                Vdyy1xzm1 = Array<OneD, NekDouble>(totszedges1d)   ;
                Vdxy1xzm1   = Array<OneD, NekDouble>(totszedges1d);
                Vdzy1xzm1 = Array<OneD, NekDouble>  (totszedges1d) ;
                edgeptsin[0] = edgexyztemp;
                edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
                E->PhysEvalBasisGrad(edgeptsin,Vy1xzm1, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1 );  

                // edge left bot (z = -1), (x = -1) AD
                Vxm1yzm1 = Array<OneD, NekDouble>(totszedges1d);
                Vdyxm1yzm1 = Array<OneD, NekDouble>(totszedges1d)   ;
                Vdxxm1yzm1 = Array<OneD, NekDouble>(totszedges1d)   ;
                Vdzxm1yzm1 = Array<OneD, NekDouble>(totszedges1d)   ;
                
                edgeptsin[1] = edgexyztemp;
                edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
                edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
                E->PhysEvalBasisGrad(edgeptsin,Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);  
                
                //edge right bot ( z = -1) (x = 1) BC
                Vx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
                Vdyx1yzm1 =Array<OneD, NekDouble>(totszedges1d);
                Vdxx1yzm1= Array<OneD, NekDouble>(totszedges1d);
                Vdzx1yzm1  = Array<OneD, NekDouble>(totszedges1d);
                

                edgeptsin[1] = edgexyztemp;
                edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);           
                edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);           
                E->PhysEvalBasisGrad(edgeptsin,Vx1yzm1, Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);   

            }
            if(numsurfaces == 5) // pyr
            {
                //we need convolution of qzin + qzinmid as input points
                int totszsurf2d = (demo.testcoord2dquad[0].size())*E->GetNcoeffs();
              
                Array<OneD, Array<OneD, NekDouble> > surfptsin(dimension), surfptsintemp(dimension);
                for(int k = 0; k < dimension-1; k++)
                {
                    surfptsin[k] = Array<OneD, NekDouble>(demo.testcoord2dquad[k]);
                    surfptsintemp[k] = Array<OneD, NekDouble>(demo.testcoord2dquad[k]);
                }

                surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.testcoord2dquad[0].size()); 
                surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.testcoord2dquad[0].size()); 

                //surface bot (quad)z = -1
                Vxyzm1 =  Array<OneD, NekDouble>(totszsurf2d);
                surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);

                E->PhysEvalBasisGrad(surfptsintemp,Vxyzm1, NullNekDouble1DArray,NullNekDouble1DArray,NullNekDouble1DArray);   


                totszsurf2d = (demo.testcoord2dtri[0].size())*E->GetNcoeffs();
              
                for(int k = 0; k < dimension-1; k++)
                {
                    surfptsin[k] = Array<OneD, NekDouble>(demo.testcoord2dtri[k]);
                    surfptsintemp[k] = Array<OneD, NekDouble>(demo.testcoord2dtri[k]);
                }

                surfptsin[dimension-1] = Array<OneD, NekDouble>(demo.testcoord2dtri[0].size()); 
                surfptsintemp[dimension-1] = Array<OneD, NekDouble>(demo.testcoord2dtri[0].size()); 

                //surface hypt (tri) x+z = 0 && y + z = 0
		Array<OneD, Array<OneD, NekDouble> >surfprojecttmp(3);
		Array<OneD, Array<OneD, NekDouble> >surfproject(3);
		surfprojecttmp[0] = Array<OneD, NekDouble> (Etri->GetBasis(0)->GetZ().size()*Etri->GetBasis(1)->GetZ().size());//Etri->GetBasis(0)->GetZ();
		surfprojecttmp[1] = Array<OneD, NekDouble> (surfprojecttmp[0].size());
		surfprojecttmp[2] = Array<OneD, NekDouble> (surfprojecttmp[0].size());
		int ctrn = 0;
		for(int k = 0; k < Etri->GetBasis(0)->GetZ().size(); k++)
		  {
		    for(int p = 0; p < Etri->GetBasis(1)->GetZ().size(); p++)
		      {
			if( (Etri->GetBasis(0)->GetZ())[k] + (Etri->GetBasis(1)->GetZ())[p] <=0)
			  {
			    surfprojecttmp[0][ctrn] = (Etri->GetBasis(0)->GetZ())[k];
			    surfprojecttmp[1][ctrn] = (Etri->GetBasis(1)->GetZ())[p];
			    ctrn++;
			  }
		      }
		  }
		surfproject[0] = Array<OneD, NekDouble> (ctrn);
		surfproject[1] = Array<OneD, NekDouble> (ctrn);
		surfproject[2] = Array<OneD, NekDouble> (ctrn);
		Vmath::Vcopy(ctrn, &surfprojecttmp[0][0], 1, &surfproject[0][0], 1);
		
		Vmath::Smul(surfproject[0].size(), -1.0, & surfproject[0][0], 1,  &surfproject[2][0], 1);
		Vmath::Vcopy(surfproject[0].size(), &surfproject[0][0], 1,  &surfproject[1][0], 1);

		Vxmzyproject  =   Array<OneD, NekDouble>(Etri->GetBasis(0)->GetZ().size()*Etri->GetBasis(1)->GetZ().size()*E->GetNcoeffs());
		E->PhysEvalBasisGrad(surfproject, Vxmzyproject, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

		
                //surface left (tri) x = -1
                Vxm1yzproject =  Array<OneD, NekDouble>(Etri->GetBasis(0)->GetZ().size()*Etri->GetBasis(1)->GetZ().size()*E->GetNcoeffs());
                surfptsintemp[0] = Array<OneD, NekDouble>(surfproject[0].size(), -1.0);
                surfptsintemp[1] = Array<OneD, NekDouble>(surfproject[0]);

                surfptsintemp[2] = Array<OneD, NekDouble>(surfproject[1]);
                
                E->PhysEvalBasisGrad(surfptsintemp, Vxm1yzproject, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);   


		//surface front (tri) y = -1
                Vxym1zproject =  Array<OneD, NekDouble>(Etri->GetBasis(0)->GetZ().size()*Etri->GetBasis(1)->GetZ().size()*E->GetNcoeffs());
                surfptsintemp[1] = Array<OneD, NekDouble>(surfproject[0].size(), -1.0);
		surfptsintemp[0] = Array<OneD, NekDouble>(surfproject[0]);
		surfptsintemp[2] =  Array<OneD, NekDouble>(surfproject[1]);

                E->PhysEvalBasisGrad(surfptsintemp, Vxym1zproject, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

                
		//surface back (tri) y = -z
		Vxymzproject =  Array<OneD, NekDouble>(Etri->GetBasis(0)->GetZ().size()*Etri->GetBasis(1)->GetZ().size()*E->GetNcoeffs());
		surfptsintemp[1] = Array<OneD, NekDouble>(surfproject[1]);
		Vmath::Smul(surfproject[0].size(), -1.0, & surfproject[1][0], 1,  &surfptsintemp[2][0], 1);
		E->PhysEvalBasisGrad(surfptsintemp,Vxymzproject, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);

                
            }
        }
        
        if(Opt_needed(coeffs))
        {

            cout<<"\n need optimization";
	    //            WriteCSV(coeffs);

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
                //WriteCSV(coeffs, 1);

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
  //  cout<<"\n begin: "<<vec[0].size()<<"\n";;
  //    cout<<"\n vec size = "<<vec.size()<<" "<<vec[0].size()<<"\n";
    Array<OneD, Array<OneD, NekDouble> > ret(vec.size());
    for (int k = 0; k < vec.size(); k++)
      ret[k] = Array<OneD, NekDouble>(vec[0].size());
    int ct = 0;
    for (int k = 0; k < vec[0].size(); k++)
    {
      //   ret[k ] = Array<OneD, NekDouble>(vec[k].size(), vec[k].data());
         

      if(vec[0][k] + vec[2][k] <= 0 && vec[2][k] + vec[1][k] <= 0 )
	{

	  ret[0][ct] = vec[0][k];
	  ret[1][ct] = vec[1][k];
	  ret[2][ct] = vec[2][k];
	  ct++;
	}
        
    }
    //cout<<"\n ct = "<<ct<<" ";
    Array<OneD, Array<OneD, NekDouble> > ret2(vec.size());
    ret2[0] = Array<OneD, NekDouble>(ct);
    ret2[1] = Array<OneD, NekDouble>(ct);
    ret2[2] = Array<OneD, NekDouble>(ct);
    Vmath::Vcopy(ct, ret[0], 1, ret2[0], 1);
    
    Vmath::Vcopy(ct, ret[1], 1, ret2[1], 1);
    
    Vmath::Vcopy(ct, ret[2], 1, ret2[2], 1);
      
    return ret2;
}


void project_surfaces( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&retquad, Array<OneD, Array<OneD, NekDouble> >&rettri, int d)
{
    boost::ignore_unused(d);

    //bot surface  z = -1
    surfuhats(uhats, retquad[0], Vxyzm1, 0);

    //left surface x = -1
    surfuhats(uhats, rettri[0], Vxmzyproject, 1);

    //surface hypt x+z = 0
    surfuhats(uhats, rettri[1], Vxm1yzproject,  1);

    //front surface y = -1  
    surfuhats(uhats, rettri[2], Vxym1zproject,  1);
 
    //surface back y = 1  
    surfuhats(uhats, rettri[3], Vxymzproject, 1);

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
    if(numedges == 8) // pyr
    {

        if(d == 1)
        {

            // edge front left (x = -1) (y = -1)
            edgederpquhats(uhats, ret[0], E->GetNcoeffs(),  Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);
                        
            //edge front hypt (y = -1) (z + x = 0)
            edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vym1xmz, Vdxym1xmz, Vdyym1xmz, Vdzym1xmz);

	    //edge back hypt ( y = -z) (z +x = 0)
            edgederpquhats(uhats, ret[2], E->GetNcoeffs(),  Vxeyxmz, Vdxxeyxmz, Vdyxeyxmz, Vdzxeyxmz);
   	    
	    //edge back left (y = -z) (x = 1) ED
            edgederpquhats(uhats, ret[3], E->GetNcoeffs(), Vx1ymz, Vdxx1ymz, Vdyx1ymz, Vdzx1ymz);
	        
            //edge front bot (y = -1) (z = -1)
            edgederpquhats(uhats, ret[4], E->GetNcoeffs(), Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);
     
            //edge back bot (y = 1) (z = -1))
            edgederpquhats(uhats, ret[5], E->GetNcoeffs(), Vy1xzm1,  Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);
            
            // edge left bot (z = -1), (x = -1)
            edgederpquhats(uhats, ret[6], E->GetNcoeffs(), Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);
            
            //edge right bot ( z = -1) (x = 1)
            edgederpquhats(uhats, ret[7], E->GetNcoeffs(),   Vx1yzm1, Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);
        
     
        }
        else if( d == 0)
        {
            // edge front left (x = -1) (y = -1)
            deruhatsedges(uhats, ret[0], Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);
        
            //edge front hypt (y = -1) (z + x = 0)
            deruhatsedges(uhats, ret[1],Vym1xmz, Vdxym1xmz, Vdyym1xmz, Vdzym1xmz);

	    //edge back hypt ( y = 1) (z +x = 0)
	    deruhatsedges(uhats, ret[2],  Vxeyxmz, Vdxxeyxmz, Vdyxeyxmz, Vdzxeyxmz);

	    //edge back left (y = -z) (x = -1) ED
	    deruhatsedges(uhats, ret[3], Vx1ymz, Vdxx1ymz, Vdyx1ymz, Vdzx1ymz); 

            //edge front bot (y = -1) (z = -1)
            deruhatsedges(uhats, ret[4], Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);

	    //edge back bot (y = 1) (z = -1)) 
	    deruhatsedges(uhats, ret[5],Vy1xzm1,  Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);

            // edge left bot (z = -1), (x = -1)
	    deruhatsedges(uhats, ret[6], Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);
            
            //edge right bot ( z = -1) (x = 1)
            deruhatsedges(uhats, ret[7],  Vx1yzm1,  Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);


        }
    }
   
}

FILE *csvfile(int flag)
{
    static char name[LOGNAME_SIZE];
    char fullname[] = "Pyr_";//[LOGNAME_SIZE+4];                                                                                                                                                           
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


// sig = 0 -> called by opt_needed
// sig = 1 -> called by sphere_rotation
Array<OneD, Array< OneD,  NekDouble> > call_find_roots(Array<OneD,  NekDouble> &uhats , int d, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhatsquad,  Array<OneD, Array<OneD, NekDouble> >&surfaceuhatstri,  int sig)
{

    boost::ignore_unused(d);
    int dimension = E->GetShapeDimension(); 
    vector<vector< NekDouble> > ret(dimension);

    Array<OneD, Array<OneD, NekDouble> > retarr;//(dimension);

    if(numedges == 8) //pyr
    {
        for(int ii = 0; ii < numedges; ii++)
        {
	
            vector<vector<NekDouble> > tmp(dimension);
 
            tmp  = (demo.find_roots(uhatsedges[ii], E, 0, sig, 0, 0, 0)) ;
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
                else if(ii == 2) //edge back hypt ( y = -z) (z +x = 0)    
                {
                    ret[0].push_back(tmp[0][p]);
                    ret[1].push_back(tmp[0][p]);
                    ret[2].push_back(-tmp[0][p]);

		    ret[0].push_back(-tmp[0][p]);
                    ret[1].push_back(-tmp[0][p]);
                    ret[2].push_back(tmp[0][p]);
                    
                    
                }
                else if(ii == 3)  //edge back left (y = -z) (x = -1) ED
                {
                    ret[0].push_back(-1);
                    ret[1].push_back(-tmp[0][0]);
                    ret[2].push_back(tmp[0][p]);

		    ret[0].push_back(-1);
                    ret[1].push_back(tmp[0][0]);
                    ret[2].push_back(-tmp[0][p]);
		}
                else if(ii == 4)  //edge front bot (y = -1) (z = -1)   
                {
                    ret[1].push_back(-1);
                    ret[0].push_back(tmp[0][p]);
                    ret[2].push_back(-1);
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
                else if(ii == 7) //edge right bot ( z = -1) (x = 1)
                {
                    ret[2].push_back(-1);
                    ret[1].push_back(tmp[0][p]);
                    ret[0].push_back(1);
                }
                    
            }
            

            
        }
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
	      tmp  = (demo.find_roots(surfaceuhatsquad[ii], Equad, 0, sig, 1, 0, 0)) ;  
            

                for(int p = 0; p < tmp[0].size(); p++)
                {

                    ret[2].push_back(-1);
                    ret[1].push_back(tmp[1][p]);
                    ret[0].push_back(tmp[0][p]);
                }
            }
                    
            else if(ii == 1) //surf left x = -1 (tri)
            {

                // surfid = 0 for quad;
	      tmp  = (demo.find_roots(surfaceuhatstri[0], Etri, 0, sig, 1, 0, 1)) ;  
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
	      tmp  = (demo.find_roots(surfaceuhatstri[1], Etri,  0, sig, 1, 0, 1)) ;

                for(int p = 0; p < tmp[0].size(); p++)
                {  

                    ret[2].push_back(-tmp[0][p]);
                    ret[1].push_back(tmp[0][p]);
                    ret[0].push_back(tmp[0][p]);
                    ret[2].push_back(tmp[0][p]);
		    ret[1].push_back(tmp[0][p]);
		    ret[0].push_back(-tmp[0][p]);

                }
            }
            else if(ii == 3) //surf front y = -1
            {
                    
                // surfid = 1 for tri;
	      tmp  = (demo.find_roots(surfaceuhatstri[2], Etri, 0, sig, 1, 0, 1)) ;  

                for(int p = 0; p < tmp[0].size(); p++)
                {
                    ret[1].push_back(-1);
                    ret[2].push_back(tmp[1][p]);
                    ret[0].push_back(tmp[0][p]);
                
                }
            }
            else if(ii == 4) //surf back y = -z, x = y
            {

                // surfid = 0 fortri;
	      tmp  = (demo.find_roots(surfaceuhatstri[3], Etri, 0, sig, 1, 0, 1)) ;  
                for(int p = 0; p < tmp[0].size(); p++)
                {
                    ret[1].push_back(tmp[0][p]);
                    ret[2].push_back(-tmp[0][p]);
                    ret[0].push_back(tmp[0][p]);
                }
                    
            }
        }
    }
            
    //find interior roots:
    vector<vector<NekDouble> > tmp;
    tmp =  (demo.find_roots(uhats, E, 0, sig, 0, 1, 0)) ;
    //    cout<<"\n 3d GD root:\n";
    if(tmp[0].size() > 0)
    {
        for(int p = 0; p <dimension; p++)
        {
            ret[p].push_back(tmp[p][0]);
	    //	    cout<< " "<<tmp[p][0]<<" ";
	}
    }
    //cout<<"\n";

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
      

    retarr = vectoarr(ret);

    return retarr;
}


int Opt_needed(Array<OneD, NekDouble> uhats, int flag )
{
    int totModes = uhats.size();
    Array<OneD,  Array<OneD, NekDouble> > rootsarr;
    Array<OneD, Array<OneD, NekDouble> > edgeuhats  (numedges);
    Array<OneD, Array<OneD, NekDouble> > surfacederuhatsquad  (1);
    Array<OneD, Array<OneD, NekDouble> > surfacederuhatstri  (4);

    for(int t = 0; t < 1; t++)
    {
        surfacederuhatsquad[t] = Array<OneD, NekDouble>(Equad->GetNcoeffs());
    } 
    for(int t = 0; t < 4; t++)
    {
        surfacederuhatstri[t] = Array<OneD, NekDouble>(Etri->GetNcoeffs());
    } 

    project_surfaces(uhats, surfacederuhatsquad, surfacederuhatstri);
            
    project_edges(uhats, edgeuhats);

    //find roots of der of uhats, that is why d = 0


    rootsarr = call_find_roots(uhats, 0, edgeuhats, surfacederuhatsquad, surfacederuhatstri, 0);
    // cout<<"\n roots:\n";
    // for(int p = 0; p < rootsarr[0].size(); p++)
    //   {
    // 	cout<<" "<<rootsarr[0][p]<<" "<<rootsarr[1][p]<<" "<<rootsarr[2][p]<<"\n";
    // 	if(rootsarr[0][p]+rootsarr[2][p] > 0 && abs(rootsarr[0][p]+rootsarr[2][p])>1e-10)
    // 	  {
    // 	    cout<<"\n panik!x = "<<rootsarr[0][p]<<" y = "<<rootsarr[1][p]<<" z ="<<rootsarr[2][p]<<"\n \n";
    // 	    //		exit(0);
    // 	      }
    // 	if(rootsarr[1][p]+rootsarr[2][p] > 0 &&abs(rootsarr[1][p]+rootsarr[2][p] )>1e-10)
    // 	      {
    // 		cout<<"\n panik! x = "<<rootsarr[0][p]<<" y="<<rootsarr[1][p]<<" z ="<<rootsarr[2][p]<<"\n \n";
    // 		//	exit(0);
    // 	      }
    //   }
    // exit(0);
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
    

    Array<OneD, Array<OneD, NekDouble> > surfaceuhatsquad  (1);
    Array<OneD, Array<OneD, NekDouble> > surfaceuhatstri  (4);
    Array<OneD, Array<OneD, NekDouble> > Pf(numedges);
    for(int k= 0; k < numedges; k++)
    {
            
        Pf[k] = Array<OneD, NekDouble>(E3seg->GetNcoeffs());//(3*(E->GetBasis(0)->GetNumModes()-1));        
    }
    for(int k = 0; k < 1; k++)
      {
	surfaceuhatsquad[k] = Array<OneD, NekDouble>(Equad->GetNcoeffs());
      }
    for(int k = 0; k < 4; k++)
      {
	surfaceuhatstri[k] = Array<OneD, NekDouble>(Etri->GetNcoeffs());
      }
    
    NekDouble pqval;
    while (counter <= niter)
    {
        pqval = inf;
        utemp = d.back();

        //pq = @(xx) cfun(xx) .* (constraints(1,jj).*Vc(xx)*uhats);  
        project_edges(utemp, Pf, 1);
        
        
        project_surfaces(utemp, surfaceuhatsquad, surfaceuhatstri, 1);
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
	//	cout<<"\n";
	// int flgpnk = 0;
	//	cout<<"\n optima sz = "<<optima[0].size()<<" ";
	for(int p = 0; p < optima[0].size(); p++)
	  {
	    
	    //		cout<<"\n x = "<<optima[0][p]<<"y = "<<optima[1][p]<<"  z ="<<optima[2][p]<<"\n \n";

	    if(optima[0][p]+optima[2][p] > 0)
	      {
		cout<<"\n panik!x = "<<optima[0][p]<<"y = "<<optima[1][p]<<"  z ="<<optima[2][p]<<"\n \n";
		//flgpnk = 1;
				exit(0);
	      }
	    if(optima[1][p]+optima[2][p] > 0)
	      {
		cout<<"\n panik! x = "<<optima[0][p]<<"  y="<<optima[1][p]<<" z ="<<optima[2][p]<<"\n \n";
		//flgpnk = 1;
		exit(0);
	      }
	   
	  }
	// if(flgpnk == 1)
	//   exit(0);
	wsp1 = Array<OneD, NekDouble>(optima[0].size());
        
        // Vtmp is evaluation of basis fn at optima
        Array<OneD, NekDouble> Vtmp(optima[0].size() * E->GetNcoeffs());
        E->PhysEvalBasisGrad(optima, Vtmp, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
        demo.pq(utemp, optima, Vtmp,  pqvalcoords, wsp1);


        if (pqvalcoords[0] < pqval)
        {
            for(int k = 0; k  < dimension; k++)
            {
                xastarr[k] = pqvalcoords[k+1];
            }
            pqval = pqvalcoords[0];
        }
        
        cout<<"\n at counter="<<counter<<" min val="<<pqval<<" xast ="<<xastarr[0]<<" "<<xastarr[1]<<" "<<xastarr[2]<< " o = "<<E->GetBasis(0)->GetNumModes()<<" p = "<<E->GetBasis(0)->GetNumPoints();
        xastarrhist[0] = xastarr[0];   
        xastarrhist[1] = xastarr[1];        
        xastarrhist[2] = xastarr[2];
        // If minimum is non-negative, we're done
        if (pqval >= -tols.at(0))
        {
            break;
        }
        
        Array<OneD, NekDouble> Vast(uhats.size()), Vastsq(uhats.size());
        Array<OneD, Array<OneD, NekDouble> > xastarr2(3);
        xastarr2[0] = Array<OneD, NekDouble>(1, xastarr[0]);
        xastarr2[1] = Array<OneD, NekDouble>(1, xastarr[1]);
        xastarr2[2] = Array<OneD, NekDouble>(1, xastarr[2]);
        NekDouble vastsqsum;

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

        d.push_back(qast);

        counter = counter + 1;
    
    
    }
    cout<<"sphere_rotation took "<<counter<<"startv = "<<startval<<" at startcoords("<<startcoordx<<","<<startcoordy<<" ," <<startcoordz<<")   iterations\n  GD iter stepszie ="<<demo.GetTol()<<" maxiters = "<<demo.GetMaxIter()<<" o = "<<E->GetBasis(0)->GetNumModes()<<" p = "<<E->GetBasis(0)->GetNumPoints()<<" pyr \n";;
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
