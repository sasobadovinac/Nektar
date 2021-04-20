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
         vector<vector<  NekDouble> >roots,
         StdExpansion *tempE,
         Array<OneD,NekDouble> &pqevalxast = NullNekDouble1DArray,
         Array<OneD,NekDouble>&fvals = NullNekDouble1DArray);


//declare Do_optimize
void Do_optimize(Array<OneD, NekDouble> &uhats);

//declare find_roots, flag is 0 for confederate matrix approach
vector<vector<  NekDouble> > find_roots(Array<OneD, NekDouble> &uhats, int d = 0, int flag = 0,  int sig = 0, int surfflag = 0);

// declare caller routine to find_roots
// flag = 0 -> opt_needed calls
// flag = 1 -> sphere_rot calls
vector<vector<  NekDouble> > call_find_roots(Array<OneD, NekDouble> &uhatsall, int d = 0, Array< OneD, Array<OneD, NekDouble> >&uhatsedges =NullNekDoubleArrayofArray , Array< OneD, Array<OneD, NekDouble> >&surfaceuhats =NullNekDoubleArrayofArray ,  int flag = 0 );


//declare Opt_needed
int Opt_needed(Array<OneD, NekDouble> uhats);


// for 2D elements to get uhats at edges
// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble> uhats,Array<OneD, Array<OneD, NekDouble> >&ret , int d = 0);


void derpq(Array<OneD, NekDouble> &uhats,  NekDouble &ret,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd);

// uhatpqd stuff: 
// called by project_edges if d = 1 is passed in project_edges
void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret, int modes,  Array<OneD, NekDouble> V3, Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1 = NullNekDouble1DArray , Array<OneD, NekDouble> Vxyd2 = NullNekDouble1DArray);

void deruhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble >&ret,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2 = NullNekDouble1DArray, int surfflag = 0 );

void project_surfaces( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret, int sig = 0 );

vector< vector<NekDouble> > gradient_descent(Array<OneD, NekDouble> uhats, int sig, int surfflag );


// Colleague matrix
int numedges, numsurfaces;
Array<OneD, Array<OneD, NekDouble> > C;
StdExpansion *E;
StdExpansion *E3seg;
StdExpansion *E3quad;
Array<OneD, NekDouble> qZin; //inner point grid for 2D rootfinding
Array<OneD, NekDouble> qWin; //inner point grid for 2D rootfinding
Array<OneD, NekDouble> qx;
Array<OneD, NekDouble> qw;

Array<OneD, NekDouble> V;
Array<OneD, NekDouble> Vd;
// Array<OneD, NekDouble >Vall;
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
Array<OneD, Array<OneD, NekDouble> > qZinarr;
Array<OneD, NekDouble > qWinarr;

// for triangle edges root-finding:
Array<OneD, NekDouble> Vxyhyp;
Array<OneD, NekDouble> Vdxxyhyp;
Array<OneD, NekDouble> Vdyxyhyp;

// for hex edges root finding

// edge front left (x = -1) (z = 1)
Array<OneD, NekDouble> Vxm1yz1   ;
Array<OneD, NekDouble> Vdyxm1yz1   ;
Array<OneD, NekDouble> Vdxxm1yz1   ;
Array<OneD, NekDouble> Vdzxm1yz1   ;
                
//edge front right (x = 1) (z = 1)
Array<OneD, NekDouble> Vx1yz1   ;
Array<OneD, NekDouble> Vdyx1yz1   ;
Array<OneD, NekDouble> Vdxx1yz1   ;
Array<OneD, NekDouble> Vdzx1yz1   ;
               
//edge front top (y = 1) (z = 1)
Array<OneD, NekDouble> Vy1xz1   ;
Array<OneD, NekDouble> Vdyy1xz1   ;
Array<OneD, NekDouble> Vdxy1xz1   ;
Array<OneD, NekDouble> Vdzy1xz1   ;
               
//edge front bot (y = -1) (z = 1)
Array<OneD, NekDouble> Vym1xz1   ;
Array<OneD, NekDouble> Vdyym1xz1   ;
Array<OneD, NekDouble> Vdxym1xz1   ;
Array<OneD, NekDouble> Vdzym1xz1   ;


// edge back left (z = -1), (x = -1)
Array<OneD, NekDouble> Vxm1yzm1   ;
Array<OneD, NekDouble> Vdyxm1yzm1   ;
Array<OneD, NekDouble> Vdxxm1yzm1   ;
Array<OneD, NekDouble> Vdzxm1yzm1   ;
                
//edge back right (x = 1), (z = -1))
Array<OneD, NekDouble> Vx1yzm1   ;
Array<OneD, NekDouble> Vdyx1yzm1   ;
Array<OneD, NekDouble> Vdxx1yzm1   ;
Array<OneD, NekDouble> Vdzx1yzm1   ;
               
//edge back top ( y = 1) (z = -1)
Array<OneD, NekDouble> Vy1xzm1   ;
Array<OneD, NekDouble> Vdyy1xzm1   ;
Array<OneD, NekDouble> Vdxy1xzm1   ;
Array<OneD, NekDouble> Vdzy1xzm1   ;
               
//edge back bot (y = -1) (z = -1))
Array<OneD, NekDouble> Vym1xzm1   ;
Array<OneD, NekDouble> Vdyym1xzm1   ;
Array<OneD, NekDouble> Vdxym1xzm1   ;
Array<OneD, NekDouble> Vdzym1xzm1   ;

                

// edge left bot (y = -1), (x = -1)
Array<OneD, NekDouble> Vxm1ym1z   ;
Array<OneD, NekDouble> Vdyxm1ym1z   ;
Array<OneD, NekDouble> Vdxxm1ym1z   ;
Array<OneD, NekDouble> Vdzxm1ym1z   ;
                
//edge left top (x = -1), (y = 1))
Array<OneD, NekDouble> Vxm1y1z   ;
Array<OneD, NekDouble> Vdyxm1y1z   ;
Array<OneD, NekDouble> Vdxxm1y1z   ;
Array<OneD, NekDouble> Vdzxm1y1z   ;
               
//edge right bot ( y = -1) (x = 1)
Array<OneD, NekDouble> Vx1ym1z   ;
Array<OneD, NekDouble> Vdyx1ym1z   ;
Array<OneD, NekDouble> Vdxx1ym1z   ;
Array<OneD, NekDouble> Vdzx1ym1z   ;
               
//edge right top (y  1) (x  1))
Array<OneD, NekDouble>Vy1x1z ;
Array<OneD, NekDouble>Vdyy1x1z ;
Array<OneD, NekDouble>Vdxy1x1z ;
Array<OneD, NekDouble> Vdzy1x1z ;

//surface bot y = -1
Array<OneD, NekDouble>Vym1xz ;
Array<OneD, NekDouble>Vdyym1xz ;
Array<OneD, NekDouble>Vdxym1xz ;
Array<OneD, NekDouble> Vdzym1xz ;


//surface right x = 1
Array<OneD, NekDouble>Vx1yz ;
Array<OneD, NekDouble>Vdyx1yz ;
Array<OneD, NekDouble>Vdxx1yz ;
Array<OneD, NekDouble> Vdzx1yz ;


//surface top y = 1
Array<OneD, NekDouble>Vy1xz ;
Array<OneD, NekDouble>Vdyy1xz ;
Array<OneD, NekDouble>Vdxy1xz ;
Array<OneD, NekDouble> Vdzy1xz ;
 

//surface left x = -1
Array<OneD, NekDouble>Vxm1yz ;
Array<OneD, NekDouble>Vdyxm1yz ;
Array<OneD, NekDouble>Vdxxm1yz ;
Array<OneD, NekDouble> Vdzxm1yz ;


//surface front z = 1
Array<OneD, NekDouble>Vz1xy ;
Array<OneD, NekDouble>Vdyz1xy ;
Array<OneD, NekDouble>Vdxz1xy ;
Array<OneD, NekDouble> Vdzz1xy ;
  

//surface back z = -1
Array<OneD, NekDouble>Vzm1xy ;
Array<OneD, NekDouble>Vdyzm1xy ;
Array<OneD, NekDouble>Vdxzm1xy ;
Array<OneD, NekDouble> Vdzzm1xy ;
  

int dimension ;
Array<OneD, NekDouble > interioreval2d;
Array<OneD, NekDouble > interioreval2d0;
Array<OneD, NekDouble > interioreval2d1;
Array<OneD, NekDouble > interioreval3d;
Array<OneD, NekDouble > interioreval3d0;
Array<OneD, NekDouble > interioreval3d1;
Array<OneD, NekDouble > interioreval3d2;

Array<OneD, Array<OneD, NekDouble> > testcoord3d; 
Array<OneD, Array<OneD, NekDouble> > testcoord2d; 
            
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
    std::vector<int> order;
    std::vector<BasisType> btype(3, eNoBasisType);
    LibUtilities::ShapeType stype = E->DetShapeType();
    LibUtilities::PointsType pointsTypeCheb = LibUtilities::eGaussGaussChebyshev;

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
            //            x[i] = ((1-0)*(x[i] + 1))/(2) + 0;
            //y[i] = ((1-0)*(y[i] + 1))/(2) + 0; 
            //z[i] = ((1-0)*(z[i] + 1))/(2) + 0; 
            sol[i] =(sin(x[i]*2*M_PI))*(sin(y[i]*2*M_PI))*sin(M_PI*2*z[i]) +0.998459440387815;
        
        }

        else if(dimension ==2) //only quad
        {
            sol[i] = sin(2*M_PI*x[i])*sin(2*M_PI*y[i]) -  exp(-1) + 1.37 ;//- 0.0602135;
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
    
    // check for -ve values and apply opt if necessary
    if (vm.count("optm"))
    {
        int dimension = E->GetShapeDimension(); 
        Array<OneD, Array<OneD, NekDouble> > qxarr(dimension);
        qx = E->GetBasis(0)->GetZ();
        qw = E->GetBasis(0)->GetW();            
        vector<LibUtilities::BasisKey> tmpBasisKey = demo.GetBasisKey();       
        LibUtilities::PointsKey pkeycheb(E->GetTotPoints(), pointsTypeCheb);

        LibUtilities::BasisKey bkeycheb(LibUtilities::eChebyshev, 3*(E->GetBasis(0)->GetNumModes()-1),  pkeycheb);
        
        E3seg = new StdSegExp(bkeycheb);
        
        if (E3seg == nullptr)
        {
            return 1;
        }
        // V3 is required by edge root finding (even for 1D)
        V3 = Array<OneD, NekDouble>(E3seg->GetTotPoints()*E3seg->GetNcoeffs());

            qxarr[0] = E->GetBasis(0)->GetZ();
            cout<<"\n qxarr sz = "<<qxarr.size()<<" "<<qxarr[0].size()<<"\n";
            E3seg->PhysEvalBasisGrad(qxarr, V3, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray );
        
        
        
        // Array<OneD, Array<OneD, NekDouble> >temp2(1);
        // temp2[0] = qx;
        
        // for(int pp = 0;  pp<(3*(E->GetNcoeffs()-1)); pp++)
        // {
        //     for(int yy =0; yy < qx.size(); yy++)
        //     {
        //         V3[qx.size()*pp+yy] = E3seg->PhysEvaluateBasis(Array<OneD, NekDouble>(1,qx[yy]), pp);
        //     }   
        // }
      
        if(dimension > 1 )
        {

            int totszedges = (E->GetTotPoints())*(E->GetNcoeffs());

            // populate edge root finding matrices for quadrilateral
            LibUtilities::PointsKey  quadPointsKeyin  (order[0]+5, pointsTypeCheb);
            // for gradient descent
            qZin = (LibUtilities::PointsManager()[quadPointsKeyin ])->GetZ();
            int sz = pow(qZin.size(),2);;
            
            testcoord2d = Array<OneD, Array<OneD, NekDouble> > (2); 
            qxarr = demo.GetCoords(E);
            
            int ctsz2d = 0;
             
            for(int ii = 0; ii < 2; ++ii)
            {
                testcoord2d[ii] = Array<OneD, NekDouble>(sz);
            }
            for(int ii = 0; ii < qZin.size(); ++ii)
            {
                
                for(int jj = 0; jj < qZin.size(); ++jj)
                {
                    if(numedges == 4)
                    {
                        testcoord2d[0][ctsz2d] = qZin[ii];
                        testcoord2d[1][ctsz2d] = qZin[jj];
                        ctsz2d++;
                    }else if(numedges == 3)
                    {
                        if(qZin[ii]+qZin[jj] <=0)
                        {
                            testcoord2d[0][ctsz2d] = qZin[ii];
                            testcoord2d[1][ctsz2d] = qZin[jj];
                            ctsz2d++;
                        }
                    }
                }

            }            

            if(dimension == 2)
            {
                interioreval2d = Array<OneD, NekDouble >(E->GetNcoeffs()*testcoord2d[0].size());
                interioreval2d0 = Array<OneD, NekDouble >(E->GetNcoeffs()*testcoord2d[0].size());
                interioreval2d1 = Array<OneD, NekDouble >(E->GetNcoeffs()*testcoord2d[0].size());
                E->PhysEvalBasisGrad(testcoord2d, interioreval2d, interioreval2d0, interioreval2d1, NullNekDouble1DArray );
             
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
                Vmath::Smul(qxarr[0].size(), -1.0, qxarr[0], 1, qxarr[1] , 1);

            
                //    cout<<"\n11 qxarr sz = "<<qxarr.size()<<" "<<qxarr[0].size()<<"E3quad modes ="<<E3quad->GetNcoeffs()<<" Vxyhyp sz ="<<Vxyhyp.size()<<"\n\n";
                E->PhysEvalBasisGrad(qxarr, Vxyhyp,  Vdxxyhyp, Vdyxyhyp, NullNekDouble1DArray);  

            
            
                // left x = -1
                qxarr[1] = qxarr[0]; 
            
                qxarr[0] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
                E->PhysEvalBasisGrad(qxarr, Vxm1, Vdxxm1, Vdyxm1, NullNekDouble1DArray);
                // bot y = -1
                qxarr[0] = qxarr[1]; 
            
                qxarr[1] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
                E->PhysEvalBasisGrad(qxarr,  Vym1, Vdxym1, Vdyym1, NullNekDouble1DArray);
            
                // right x = 1
                qxarr[1] = qxarr[0]; 
            
                qxarr[0] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0); 
                E->PhysEvalBasisGrad(qxarr, Vx1, Vdxx1, Vdyx1, NullNekDouble1DArray);
            
                //top y = 1
                qxarr[0] = qxarr[1]; 
                qxarr[1] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0); 
                E->PhysEvalBasisGrad(qxarr, Vy1,Vdxy1,Vdyy1, NullNekDouble1DArray);
            
            
            }

            if(dimension == 3)
            {
                
              
                LibUtilities::BasisKey bkeychebtmp(LibUtilities::eChebyshev, tmpBasisKey[0].GetNumModes(),  pkeycheb);
                E3quad = new StdQuadExp(bkeychebtmp,bkeychebtmp);
            
                
                sz = pow(qZin.size(),3);;

                testcoord3d = Array<OneD, Array<OneD, NekDouble> > (dimension); 
                for(int ii = 0; ii < 3; ++ii)
                {
                    testcoord3d[ii] = Array<OneD, NekDouble>(sz);
                }
                int ctsz = 0;
                for(int ii = 0; ii < qZin.size(); ++ii)
                {
                    
                    for(int jj = 0; jj < qZin.size(); ++jj)
                    {

                        for(int kk = 0; kk < qZin.size(); ++kk)
                        {
                            testcoord3d[0][ctsz] = qZin[ii];
                            testcoord3d[1][ctsz] = qZin[jj];
                            testcoord3d[2][ctsz] = qZin[kk];
                            ctsz++;
                        }
                    }
                }
                qxarr = demo.GetCoords(E);

                interioreval2d = Array<OneD, NekDouble >(E3quad->GetNcoeffs()*testcoord2d[0].size());

                interioreval2d0 = Array<OneD, NekDouble >(E3quad->GetNcoeffs()*testcoord2d[0].size());
                interioreval2d1 = Array<OneD, NekDouble >(E3quad->GetNcoeffs()*testcoord2d[0].size());
                
                E3quad->PhysEvalBasisGrad(testcoord2d, interioreval2d, interioreval2d0, interioreval2d1, NullNekDouble1DArray);
                
                interioreval3d = Array<OneD, NekDouble >(E->GetNcoeffs()*testcoord3d[0].size());
                interioreval3d0 = Array<OneD, NekDouble >(E->GetNcoeffs()*testcoord3d[0].size());
                interioreval3d1 = Array<OneD, NekDouble >(E->GetNcoeffs()*testcoord3d[0].size());
                interioreval3d2 = Array<OneD, NekDouble >(E->GetNcoeffs()*testcoord3d[0].size());

                
                E->PhysEvalBasisGrad(testcoord3d, interioreval3d, interioreval3d0, interioreval3d1, interioreval3d2);

            
                if(numedges == 12) //hex
                {
                    int totszedges3d = E->GetTotPoints()*E->GetNcoeffs(); 

                    // edge front left (x = -1) (z = 1)
                    Vxm1yz1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdyxm1yz1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdxxm1yz1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdzxm1yz1 = Array<OneD, NekDouble>(totszedges3d);
                    qxarr[1] = qxarr[0];
                    qxarr[0] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
                    qxarr[2] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0);
            
                    E->PhysEvalBasisGrad(qxarr, Vxm1yz1, Vdyxm1yz1, Vdxxm1yz1, Vdzxm1yz1);  
                
                    //edge front right (x = 1) (z = 1)
                    Vx1yz1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdyx1yz1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdxx1yz1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdzx1yz1 = Array<OneD, NekDouble>(totszedges3d);
                    qxarr[0] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0);
                    qxarr[2] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0);
            
                    E->PhysEvalBasisGrad(qxarr, Vx1yz1, Vdxx1yz1, Vdyx1yz1,  Vdzx1yz1);  



                
                    //edge front top (y = 1) (z = 1)
                    Vy1xz1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdyy1xz1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdxy1xz1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdzy1xz1 = Array<OneD, NekDouble>(totszedges3d);

                    qxarr[0] = qxarr[1];
                    qxarr[2] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0);           
                    qxarr[1] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0);
            
                    E->PhysEvalBasisGrad(qxarr,Vy1xz1, Vdxy1xz1, Vdyy1xz1 , Vdzy1xz1);  

                
                    //edge front bot (y = -1) (z = 1)
                    Vym1xz1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdyym1xz1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdxym1xz1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdzym1xz1 = Array<OneD, NekDouble>(totszedges3d);
                    qxarr[2] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);           
                    qxarr[1] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0);
            
                    E->PhysEvalBasisGrad(qxarr, Vym1xz1, Vdxym1xz1,  Vdyym1xz1, Vdzym1xz1);  

 
                    // edge back left (z = -1), (x = -1)
                    Vxm1yzm1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdyxm1yzm1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdxxm1yzm1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdzxm1yzm1 = Array<OneD, NekDouble>(totszedges3d);

                    qxarr[1] = qxarr[0];
                    qxarr[0] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
                    qxarr[2] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
            
                    E->PhysEvalBasisGrad(qxarr,Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);  


                
                    //edge back right (x = 1), (z = -1))
                    Vx1yzm1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdyx1yzm1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdxx1yzm1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdzx1yzm1 = Array<OneD, NekDouble>(totszedges3d);
                
                    qxarr[0] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0);
                    qxarr[2] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
            
                    E->PhysEvalBasisGrad(qxarr,Vx1yzm1, Vdxx1yzm1, Vdyx1yzm1 , Vdxx1yzm1);  
                
                    //edge back top ( y = 1) (z = -1)
                    Vy1xzm1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdyy1xzm1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdxy1xzm1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdzy1xzm1 = Array<OneD, NekDouble>(totszedges3d);
                    qxarr[0] = qxarr[1];
                    qxarr[1] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0);
                    qxarr[2] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
            
                    E->PhysEvalBasisGrad(qxarr, Vy1xzm1, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1 );  
               
                    //edge back bot (y = -1) (z = -1))
                    Vym1xzm1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdyym1xzm1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdxym1xzm1 = Array<OneD, NekDouble>(totszedges3d);
                    Vdzym1xzm1 = Array<OneD, NekDouble>(totszedges3d);
                    qxarr[1] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
                    qxarr[2] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
            
                    E->PhysEvalBasisGrad(qxarr,Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);  


                    // edge left bot (y = -1), (x = -1)
                    Vxm1ym1z = Array<OneD, NekDouble>(totszedges3d);
                    Vdyxm1ym1z = Array<OneD, NekDouble>(totszedges3d);
                    Vdxxm1ym1z = Array<OneD, NekDouble>(totszedges3d);
                    Vdzxm1ym1z = Array<OneD, NekDouble>(totszedges3d);
                    qxarr[2] = qxarr[0];
                    qxarr[0] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
                    qxarr[1] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
            
                    E->PhysEvalBasisGrad(qxarr,Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);  
                
                    //edge left top (x = -1), (y = 1))
                    Vxm1y1z = Array<OneD, NekDouble>(totszedges3d);
                    Vdyxm1y1z = Array<OneD, NekDouble>(totszedges3d);
                    Vdxxm1y1z = Array<OneD, NekDouble>(totszedges3d);
                    Vdzxm1y1z = Array<OneD, NekDouble>(totszedges3d);
                    qxarr[0] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
                    qxarr[1] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0);
            
                    E->PhysEvalBasisGrad(qxarr,Vxm1y1z, Vdxxm1y1z, Vdyxm1y1z, Vdzxm1y1z);  

                    //edge right bot ( y = -1) (x = 1)
                    Vx1ym1z = Array<OneD, NekDouble>(totszedges3d);
                    Vdyx1ym1z = Array<OneD, NekDouble>(totszedges3d);
                    Vdxx1ym1z = Array<OneD, NekDouble>(totszedges3d);
                    Vdzx1ym1z = Array<OneD, NekDouble>(totszedges3d);
                    qxarr[0] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0);
                    qxarr[1] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
            
                    E->PhysEvalBasisGrad(qxarr,Vx1ym1z,  Vdxx1ym1z, Vdyx1ym1z, Vdzx1ym1z);   

                    //edge right top (y = 1) (x = 1))
                    Vy1x1z = Array<OneD, NekDouble>(totszedges3d);
                    Vdyy1x1z = Array<OneD, NekDouble>(totszedges3d);
                    Vdxy1x1z = Array<OneD, NekDouble>(totszedges3d);
                    Vdzy1x1z = Array<OneD, NekDouble>(totszedges3d);
                    qxarr[0] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0);
                    qxarr[1] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0);
            
                    E->PhysEvalBasisGrad(qxarr,Vy1x1z, Vdxy1x1z, Vdyy1x1z, Vdzy1x1z);   

                    //surface bot y = -1
                    Vym1xz =  Array<OneD, NekDouble>(totszedges3d);
                    Vdxym1xz =  Array<OneD, NekDouble>(totszedges3d);
                    Vdyym1xz =  Array<OneD, NekDouble>(totszedges3d);
                    Vdzym1xz =  Array<OneD, NekDouble>(totszedges3d);
                    qxarr[0] = qxarr[2];
                    qxarr[1] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
                    // qxarr[0] = qxarr[2] = qxarr[0] = qxarr[2];
            
                    E->PhysEvalBasisGrad(qxarr,Vym1xz,Vdxym1xz, Vdyym1xz, Vdzym1xz);   

                    //surface right x = 1
                    Vx1yz =  Array<OneD, NekDouble>(totszedges3d);
                    Vdxx1yz =  Array<OneD, NekDouble>(totszedges3d);
                    Vdyx1yz =  Array<OneD, NekDouble>(totszedges3d);
                    Vdzx1yz =  Array<OneD, NekDouble>(totszedges3d);
                    qxarr[0] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0);
                    qxarr[1] = qxarr[2];
      
                    E->PhysEvalBasisGrad(qxarr,Vx1yz, Vdxx1yz, Vdyx1yz, Vdzx1yz);   

                    //surface top y = 1
                    Vy1xz =  Array<OneD, NekDouble>(totszedges3d);
                    Vdxy1xz =  Array<OneD, NekDouble>(totszedges3d);
                    Vdyy1xz =  Array<OneD, NekDouble>(totszedges3d);
                    Vdzy1xz =  Array<OneD, NekDouble>(totszedges3d);
                    qxarr[1] = qxarr[0];

                    qxarr[0] = qxarr[2];

      
                    E->PhysEvalBasisGrad(qxarr,Vy1xz, Vdxy1xz, Vdyy1xz, Vdzy1xz);   

                    //surface left x = -1
                    Vxm1yz =  Array<OneD, NekDouble>(totszedges3d);
                    Vdxxm1yz =  Array<OneD, NekDouble>(totszedges3d);
                    Vdyxm1yz =  Array<OneD, NekDouble>(totszedges3d);
                    Vdzxm1yz =  Array<OneD, NekDouble>(totszedges3d);
                    qxarr[0] = Array<OneD, NekDouble>(qxarr[0].size(), 1.0);
                    E->PhysEvalBasisGrad(qxarr,Vxm1yz, Vdxxm1yz, Vdyxm1yz, Vdzxm1yz);   

                    //surface front z = 1
                    Vz1xy =  Array<OneD, NekDouble>(totszedges3d);
                    Vdxz1xy =  Array<OneD, NekDouble>(totszedges3d);
                    Vdyz1xy =  Array<OneD, NekDouble>(totszedges3d);
                    Vdzz1xy =  Array<OneD, NekDouble>(totszedges3d);
                    qxarr[2] = qxarr[0];
                    qxarr[0] = qxarr[1];
                    E->PhysEvalBasisGrad(qxarr,Vz1xy, Vdxz1xy, Vdyz1xy, Vdzz1xy);   
      
                    //surface back z = -1
                    Vzm1xy =  Array<OneD, NekDouble>(totszedges3d);
                    Vdxzm1xy =  Array<OneD, NekDouble>(totszedges3d);
                    Vdyzm1xy =  Array<OneD, NekDouble>(totszedges3d);
                    Vdzzm1xy =  Array<OneD, NekDouble>(totszedges3d);
                    qxarr[2] = Array<OneD, NekDouble>(qxarr[0].size(), -1.0);
                    E->PhysEvalBasisGrad(qxarr,Vzm1xy, Vdxzm1xy, Vdyzm1xy, Vdzzm1xy);   
                }

            }
            //            qWin = (LibUtilities::PointsManager()[quadPointsKeyin ])->GetW();
            //qZinarr = Array<OneD, Array<OneD, NekDouble> >(2);

            //qZinarr = demo.GetCoords(E);
            

            
            
        }
        else //(dimension == 1)
        {
            qZinarr = Array<OneD, Array<OneD, NekDouble> >(dimension);
            qZinarr[0] = E->GetBasis(0)->GetZ();
            Vd = Array<OneD, NekDouble>(qx.size()*coeffs.size());
            V = Array<OneD, NekDouble>(qx.size()*coeffs.size());
            
            E->PhysEvalBasisGrad(qZinarr, V, Vd, NullNekDouble1DArray,  NullNekDouble1DArray);
           
        }
        if(Opt_needed(coeffs))
        {

            cout<<"\n need optimization\n\n";
            
            Do_optimize(coeffs);
            
            cout<<"\n doopt done\n verifying...\n";//exit(0);
            cout<<"\n do_opt returning uhats\n";
            for(int ii = 0; ii < coeffs.size(); ii++)
                cout<<coeffs[ii]<<" ";
            cout<<"\n";
            if(Opt_needed(coeffs))
            {
                cout<<"\n fail\n\n";
                exit(0);
            }
            else
            {
                cout<<"\n pass\n\n";exit(0);
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

void project_surfaces( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret, int sig)
{
    boost::ignore_unused(sig);
    int modes = E->GetNcoeffs();
    int modes1D = E3quad->GetNcoeffs();
    
    if(sig == 0)
    {
        for(int k= 0; k < numsurfaces; k++)
        {
            ret[k] = Array<OneD, NekDouble>(modes);
            
        } 
    }
    else if(sig == 1)
    {
        for(int k= 0; k < numsurfaces; k++)
        {
            ret[k] = Array<OneD, NekDouble>(modes1D);
            
        } 
        
    }
    
    if(numsurfaces == 6)
    {
        //bot surface
        deruhats(uhats, ret[0], Vym1xz, Vdxym1xz, Vdyym1xz, Vdzym1xz,1 );
        //right surface
        deruhats(uhats, ret[1], Vx1yz, Vdxx1yz, Vdyx1yz, Vdzx1yz,1 );
        
        //top surface
        deruhats(uhats, ret[2], Vy1xz, Vdxy1xz ,Vdyy1xz, Vdzy1xz,1 );
        
        //left surface
        deruhats(uhats, ret[3], Vxm1yz, Vdxxm1yz ,Vdyxm1yz, Vdzxm1yz, 1);
        
        //front surface
        deruhats(uhats, ret[4], Vz1xy, Vdxz1xy ,Vdyz1xy, Vdzz1xy, 1);
        
        //back surface
        deruhats(uhats, ret[5], Vzm1xy, Vdxzm1xy ,Vdyzm1xy, Vdzzm1xy, 1);
        
    }
}


// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble>uhats,    Array<OneD, Array<OneD, NekDouble> >&ret , int d)
{
    int modes = E->GetNcoeffs();
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
            deruhats(uhats, ret[0], Vym1, Vdxym1, Vdyym1);

            //top
            deruhats(uhats, ret[1], Vx1, Vdxx1, Vdyx1);

            //top
            deruhats(uhats, ret[2], Vy1, Vdxy1, Vdyy1);
            
            //left
            deruhats(uhats, ret[3], Vxm1, Vdxxm1, Vdyxm1);
            
            
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
            deruhats(uhats, ret[0], Vym1, Vdxym1, Vdyym1);

            // left edge
            //edgederuhats(uhats, ret[1], Vall, Vxm1, Vdyxm1);
            deruhats(uhats, ret[1], Vxm1, Vdxxm1, Vdyxm1);
            // hypotenuse                 
            //edgederuhats(uhats, ret[2], Vxyhyp,Vdxxyhyp, Vdyxyhyp);
            deruhats(uhats, ret[2], Vxyhyp,Vdxxyhyp, Vdyxyhyp);
            
        }

    }
    else if(numedges == 12) // hex
    {
            // cout<<"\n d = "<<d<<" \n\n";
            // cout<<"ret sz = "<<ret.size()<<" uhats sz = "<<uhats.size()<<" Vxm1yz1.size() = "<<Vxm1yz1.size()<<"ret[0].size(0 = "<<ret[0].size()<<"\n";

        if(d == 1)
        {
            // edge front left (x = -1) (z = 1)

            edgederpquhats(uhats, ret[0], E->GetNcoeffs(), V3, Vxm1yz1, Vdxxm1yz1, Vdyxm1yz1, Vdzxm1yz1);
            //cout<<"\n out ret[1] sz = "<<Vx1yz1.size()<<"\n\n";
                        
            //edge front right (x = 1) (z = 1)  
            edgederpquhats(uhats, ret[1], E->GetNcoeffs(),V3, Vx1yz1, Vdxx1yz1, Vdyx1yz1, Vdzx1yz1);
            //edge front top (y = 1) (z = 1)      
            edgederpquhats(uhats, ret[2], E->GetNcoeffs(), V3, Vy1xz1, Vdxy1xz1, Vdyy1xz1,Vdzy1xz1);
            //edge front bot (y = -1) (z = 1)
            edgederpquhats(uhats, ret[3], E->GetNcoeffs(), V3, Vym1xz1, Vdxym1xz1, Vdyym1xz1, Vdzym1xz1);
            // edge back left (z = -1), (x = -1)
            edgederpquhats(uhats, ret[4], E->GetNcoeffs(), V3,  Vxm1yzm1, Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1);
            cout<<"\n2\n\n";     

            //edge back right (x = 1), (z = -1))
            edgederpquhats(uhats, ret[5], E->GetNcoeffs(), V3, Vx1yzm1, Vdxx1yzm1, Vdyx1yzm1, Vdzx1yzm1);
            //edge back top ( y = 1) (z = -1)
            edgederpquhats(uhats, ret[6], E->GetNcoeffs(), V3, Vy1xzm1,Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);
            cout<<"\n 3\n\n";
            //edge back bot (y = -1) (z = -1))
            edgederpquhats(uhats, ret[7], E->GetNcoeffs(), V3, Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);
            
            // edge left bot (y = -1), (x = -1)
            edgederpquhats(uhats, ret[8], E->GetNcoeffs(),V3,  Vxm1ym1z, Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);
            
            //edge left top (x = -1), (y = 1))
            edgederpquhats(uhats, ret[9], E->GetNcoeffs(), V3,  Vxm1y1z, Vdxxm1y1z,  Vdyxm1y1z, Vdzxm1y1z);
            
            //edge right bot ( y = -1) (x = 1)
            edgederpquhats(uhats, ret[10], E->GetNcoeffs(), V3,  Vx1ym1z, Vdxx1ym1z, Vdyx1ym1z, Vdzx1ym1z);
        
            //edge right top (y  1) (x  1))
            edgederpquhats(uhats, ret[11], E->GetNcoeffs(),V3, Vy1x1z, Vdxy1x1z, Vdyy1x1z, Vdzy1x1z); 

        }
        else if( d == 0)
        {

            // edge front left (x = -1) (z = 1)
            deruhats(uhats, ret[0], Vxm1yz1, Vdxxm1yz1, Vdyxm1yz1, Vdzxm1yz1);

            
            //edge front right (x = 1) (z = 1)
            deruhats(uhats, ret[1], Vx1yz1,Vdxx1yz1, Vdyx1yz1,  Vdzx1yz1 );

            
            //edge front top (y = 1) (z = 1)
            deruhats(uhats, ret[2],  Vy1xz1, Vdxy1xz1, Vdyy1xz1 , Vdzy1xz1);       

            //edge front bot (y = -1) (z = 1)
            deruhats(uhats, ret[3], Vym1xz1, Vdxym1xz1,  Vdyym1xz1, Vdzym1xz1);   

            // edge back left (z = -1), (x = -1)
            deruhats(uhats, ret[4], Vxm1yzm1,  Vdxxm1yzm1, Vdyxm1yzm1, Vdzxm1yzm1 );

            //edge back right (x = 1), (z = -1))   
            deruhats(uhats, ret[5], Vx1yzm1, Vdxx1yzm1, Vdyx1yzm1 , Vdxx1yzm1);
            //edge back top ( y = 1) (z = -1)
            deruhats(uhats, ret[6], Vy1xzm1, Vdxy1xzm1, Vdyy1xzm1, Vdzy1xzm1);
            //edge back bot (y = -1) (z = -1))
            deruhats(uhats, ret[7], Vym1xzm1, Vdxym1xzm1, Vdyym1xzm1, Vdzym1xzm1);
         
            // edge left bot (y = -1), (x = -1)
            deruhats(uhats, ret[8],Vxm1ym1z , Vdxxm1ym1z, Vdyxm1ym1z, Vdzxm1ym1z);

            //edge left top (x = -1), (y = 1))
            deruhats(uhats, ret[9], Vxm1y1z, Vdxxm1y1z, Vdyxm1y1z, Vdzxm1y1z);


            //edge right bot ( y = -1) (x = 1)
            deruhats(uhats, ret[10], Vx1ym1z, Vdxx1ym1z, Vdyx1ym1z, Vdzx1ym1z);   

            //edge right top (y = 1) (x = 1))
            deruhats(uhats, ret[11], Vy1x1z, Vdxy1x1z, Vdyy1x1z, Vdzy1x1z);   

        }
    }
   
}


void deruhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz, Array<OneD, NekDouble> Vdxyz, Array<OneD, NekDouble> Vxdyz, Array<OneD, NekDouble> Vxydz, int surfflag)
{
    boost::ignore_unused(Vxydz, Vxyz);

    int totpts = E->GetTotPoints();//ret.size();
    Array<OneD, Array<OneD,  NekDouble> >temp(dimension);;
    int modes = uhats.size();
    Array<OneD, NekDouble> vals(totpts), hold(totpts);    
    for(int k = 0; k < dimension; k++)
    {
        temp[k ] = Array<OneD, NekDouble> (modes); 
    }
    // if(surfflag == 0 )//edges
    // {
    //     for(int p = 0; p < totpts; p++)
    //     {
    //         Vmath::Vmul(modes, &Vxyz[p], totpts,  &uhats[0], 1, &temp[0][0] , 1 );
    //         hold[p] = Vmath::Vsum(modes, temp[0],1);
            
    //     }
    //     // cout<<"\nhodld.size()="<<hold.size()<<"\n";;

    //     // for(int k = 0; k < hold.size(); k++)
    //     //     cout<<hold[k]<<" ";
    //     // cout<<"\n\n";

    //     E->FwdTrans(hold, ret);
    //     // cout<<"\n edgesuhats:\n";                
    //     // cout<<"\n  Vxyz sz ="<<Vxyz.size()<<"\n";
    //     // for(int k = 0; k < Vxyz.size(); k++)
    //     //     cout<<Vxyz[k]<<" ";
    //     // cout<<"\n\n";

 
    //     // cout<<"\n uhats.size()="<<uhats.size()<<"\n";;

    //     // for(int k = 0; k < uhats.size(); k++)
    //     //     cout<<uhats[k]<<" ";
    //     // cout<<"\n\n";
        
    //     return;
    // }
    NekDouble v1;
    for(int k = 0; k < totpts; k++)
    {    
        Vmath::Vmul(uhats.size(), &Vdxyz[k], totpts, &uhats[0], 1, &temp[0][0], 1);
        v1  = Vmath::Vsum(modes, temp[0], 1);
        Vmath::Vmul(uhats.size(), &Vxdyz[k], totpts, &uhats[0], 1, &temp[1][0], 1);
        v1  = v1 + Vmath::Vsum(modes, temp[1], 1);  
        
        
        if(dimension > 2)
        {
            Vmath::Vmul(uhats.size(), &Vxydz[k], totpts, &uhats[0], 1, &temp[2][0], 1);
            v1  = v1 + Vmath::Vsum(modes, temp[2], 1);
        }
        vals[k] = v1;
    }    

    E->MultiplyByStdQuadratureMetric(vals,
                                     vals);


    for(int i = 0; i<ret.size(); i++)
    {

        Vmath::Vmul(totpts, &Vxyz[i], ret.size(), &vals[0], 1, &hold[0], 1);
        ret[i] = Vmath::Vsum(vals.size(), hold, 1);
    }
    //E->FwdTrans(vals,ret);

    if(surfflag)
    {
        // cout<<"\n surfuhats:\n";                
        // cout<<"\n ret.size()="<<ret.size()<<"\n";;

        // for(int k = 0; k < ret.size(); k++)
        //     cout<<ret[k]<<" ";
        // cout<<"\n\n";

    }
}


void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  int modes, Array<OneD, NekDouble> V3, Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2)
{
    int uhatstot = uhats.size();
    int totpts = E3seg->GetTotPoints();
    Array<OneD, NekDouble> temp(totpts), temp2(modes), temp3(uhatstot);
    Array<OneD, NekDouble> pqeval(totpts);
    NekDouble v1, v2;
    //for(int d = 0; d < dimension; d++)
    //{
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

        v1= v1 - v2*Vmath::Vsum(uhats.size(), temp2, 1);
 
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

            v1= v1 - v2*Vmath::Vsum(uhats.size(), temp2, 1);
 
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

            v1= v1 - v2*Vmath::Vsum(uhats.size(), temp2, 1);
 
            pqeval[i] += v1;
 
        }
    }
    E->MultiplyByStdQuadratureMetric(pqeval,
                                     pqeval);
    for(int i = 0; i<ret.size(); i++)
    {
        
        Vmath::Vmul(totpts, &V3[i], ret.size(), &pqeval[0], 1, &temp[0], 1);
        
        ret[i] = Vmath::Vsum(totpts, temp, 1);  
        
    }
    
}


// sig = 0 -> called by opt_needed
// sig = 1 -> called by sphere_rotation
vector<vector<  NekDouble> > call_find_roots(Array<OneD,  NekDouble> &uhats , int d, Array<OneD, Array<OneD, NekDouble> >&uhatsedges, Array<OneD, Array<OneD, NekDouble> >&surfaceuhats, int sig)
{

    boost::ignore_unused(d);
    //    cout<<"\nuhatsedges[0] size = "<<uhatsedges[0].size()<<"\n";
    
    int dimension = E->GetShapeDimension(); 
    vector<vector< NekDouble> > ret(dimension);

    if(dimension == 1 )
    {

        ret = find_roots(uhats, 1);

    }
    else if(dimension == 2)
    {
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
                tmp  = (find_roots(uhatsedges[ii], 0, 0, sig)) ;
                
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
                // if(abs(Vmath::Vmin(uhatsedges[ii].size(),uhatsedges[ii],1)) > 1e-9)
                //    {    
                    tmp  = (find_roots(uhatsedges[ii], 0, 0, sig)) ;
                // }else
                // {
                //     break;
                // }
            
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
                ret[0].push_back(0);
                ret[1].push_back(0);
            }
        }

        //find interior roots:
        vector<vector<NekDouble> > tmp;

        tmp =  (find_roots(uhats, 1, 1, sig)) ;
        if(tmp[0].size() > 0)
        {
            for(int p = 0; p <dimension; p++)
            {
                ret[p].push_back(tmp[p][0]);
            }
        }
        
    
    }
    else //3D
    {
        if(numedges == 12) //hex
        {

            for(int ii = 0; ii < numedges; ii++)
            {
    
                vector<vector<NekDouble> > tmp(dimension);
         
                // size of tmp will be dim \times no. of roots
                // if(abs(Vmath::Vmax(uhatsedges[ii].size(),uhatsedges[ii],1)) > 1e-9)
                // {
       
                tmp  = (find_roots(uhatsedges[ii], 0, 0, sig)) ;
                // for(int k = 0; k < tmp.size()-1; k++)
                // {
                //     for(int p = 0; p < tmp[k].size(); p++)
                //     {
                //         cout<<" "<<tmp[k][p]<<" ";
                //     }
                //     cout<<"\n";
                // }
                // cout<<"\n\n";
                // }else
                // {
                //     cout<<"\n in else uhats too small. ";
                //     cout<<"uhatsedgesmin="<<abs(Vmath::Vmin(uhatsedges[ii].size(),uhatsedges[ii],1))<<" d="<<d<<" \n";
                // }
                for(int p = 0; p < tmp[0].size(); p++)
                {
                    if(ii == 0) // edge front left (x = -1) (z = 1)

                    {
                        ret[1].push_back(tmp[0][p]);
                        ret[2].push_back(-1);
                        ret[0].push_back(-1);
                        
                    }
                    else if(ii == 1) // edge front right (x = 1) (z = 1)

                    {
                        ret[0].push_back(1);
                        ret[2].push_back(1);
                        ret[1].push_back(tmp[0][p]);
                    }
                    else if(ii == 2) //edge front top (y = 1) (z = 1)

                    {
                        ret[0].push_back(tmp[0][p]);
                        ret[1].push_back(1);
                        ret[2].push_back(1);
                    
                    
                    }
                    else if(ii == 3) //edge front bot (y = -1) (z = 1)

                    {
                        ret[1].push_back(-1);
                        ret[0].push_back(tmp[0][p]);
                        ret[2].push_back(1);
                    }
                    else if(ii == 4) // edge back left (z = -1), (x = -1)


                    {
                        ret[0].push_back(-1);
                        ret[1].push_back(tmp[0][p]);
                        ret[2].push_back(-1);
                    }
                    else if(ii == 5) //edge back right (x = 1), (z = -1))  
                    {
                        ret[0].push_back(1);
                        ret[1].push_back(tmp[0][p]);
                        ret[2].push_back(-1);
                    }
                    else if(ii == 6) //edge back top ( y = 1) (z = -1)
                    {
                        ret[1].push_back(1);
                        ret[0].push_back(tmp[0][p]);
                        ret[2].push_back(-1);
                    }
                    
                    else if(ii == 7) //edge back bot (y = -1) (z = -1)
                    {
                        ret[1].push_back(-1);
                        ret[0].push_back(tmp[0][p]);
                        ret[2].push_back(-1);
                    }
                    else if(ii == 8) // edge left bot (y = -1), (x = -1)
                    {
                        ret[1].push_back(-1);
                        ret[2].push_back(tmp[0][p]);
                        ret[0].push_back(-1);
                    }
                    else if(ii == 9) // edge left top (x = -1), (y = 1))
                    {
                        ret[1].push_back(1);
                        ret[2].push_back(tmp[0][p]);
                        ret[0].push_back(-1);
                    }
                    else if(ii == 10) //edge right bot ( y = -1) (x = 1)
                    {
                        ret[1].push_back(-1);
                        ret[2].push_back(tmp[0][p]);
                        ret[0].push_back(1);
                    }
                    else if(ii == 11) //edge right top (y  1) (x  1))
                    {
                        ret[1].push_back(1);
                        ret[2].push_back(tmp[0][p]);
                        ret[0].push_back(1);
                    }
                    
                }
            

            
            }
        }
        if(numsurfaces == 6)
        {   
            // call 2D rootfinder on each surface:
            vector<vector<NekDouble> > tmp(dimension);
            
            for(int ii = 0; ii < numsurfaces; ii++)
            {

                // if(abs(Vmath::Vmax(surfaceuhats[ii].size(),surfaceuhats[ii],1)) > 1e-9)
                // {
                  
                    tmp  = (find_roots(surfaceuhats[ii], 0, 1, sig, 1)) ;                  
                // }else
                // {
                //     cout<<"\n in else uhats too small. ";
                //     cout<<"uhatssurfmin="<<abs(Vmath::Vmin(surfaceuhats[ii].size(),surfaceuhats[ii],1));
                // }
                //cout<<"\n here?\n\n";
                for(int p = 0; p < tmp[0].size(); p++)
                {
                     if(ii == 0) //surf bot
                     {
                         ret[1].push_back(-1);
                         ret[2].push_back(tmp[1][p]);
                         ret[0].push_back(tmp[0][p]);
                     }
                     else if(ii == 1) //surf right
                     {
                         ret[0].push_back(1);
                         ret[2].push_back(tmp[1][p]);
                         ret[1].push_back(tmp[0][p]);
                     }
                     else if(ii == 2) //surf top
                     {
                         ret[1].push_back(1);
                         ret[2].push_back(tmp[1][p]);
                         ret[0].push_back(tmp[0][p]);
                     }
                     else if(ii == 3) //surf left
                     {
                         ret[0].push_back(-1);
                         ret[2].push_back(tmp[1][p]);
                         ret[1].push_back(tmp[0][p]);
                     }
                     else if(ii == 4) //surf front
                     {
                         ret[0].push_back(tmp[0][p]);
                         ret[2].push_back(1);
                         ret[1].push_back(tmp[1][p]);
                     }
                     else if(ii == 4) //surf back
                     {
                         ret[0].push_back(tmp[0][p]);
                         ret[2].push_back(-1);
                         ret[1].push_back(tmp[1][p]);
                     }
                     
                }
            }
        }
            
        //find interior roots:
        vector<vector<NekDouble> > tmp;

        tmp =  (find_roots(uhats, 0, 1, sig)) ;

        if(tmp[0].size() > 0)
        {
            for(int p = 0; p <dimension; p++)
            {
                ret[p].push_back(tmp[p][0]);
            }
        }

    }
    //    cout<<"\n ret size ="<<ret.size()<<" "<<ret[0].size()<<"\n\n";
    return ret;
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
vector<vector<NekDouble> > find_roots(Array<OneD, NekDouble> &uhats, int d, int flag ,int sig, int surfflag)
{

    vector<vector< NekDouble> > ret(dimension);
    //Confederate matrix approach
    if(dimension == 1 || flag == 0)
    {
        int N = uhats.size();

        vector<NekDouble> uhatsmon;

        while(abs(uhats[N-1])<1e-8)
        {
            N = N-1;
            
            if(N == 0)
            {
                //                cout<<"\n 1 N = "<<N<<"\n\n  ";
                // cout<<"\n C sz ="<<C.size()<<" "<<C[0].size()<<"\n\n";
        
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
            uhatsmon.push_back(Vmath::Vsum(N, &temp[0], 1));
        }
        
        N = uhatsmon.size();

        // truncate trailing zeros
        
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
            
        }
        else //d == 0
        {
            uhatsdiff = uhatsmon;
            
        }

        Vmath::Smul(N, 1.0/uhatsdiff[N-1], &uhatsdiff[0], 1, &uhatsdiff[0], 1);

        vector<NekDouble> EIG_R = demo.FindEigenval(uhatsdiff, N);
        //        cout<<"\n EIG_R[kk]: "<<EIG_R.size()<<"\n";
        for(int kk = 0; kk <EIG_R.size(); kk++)
        {//cout<<EIG_R[kk]<<" ";
            ret[0].push_back( EIG_R[kk] );
        }

    }
    else if(dimension > 1 && flag == 1)
    {    

        ret = gradient_descent(uhats, sig, surfflag);
       
    }
    //    cout<<"\n dim = "<<dimension<<" surfflag = "<<surfflag<<" flag = "<<flag<<"\n";

    return ret;    
    
}

// sig = 1 -> from do_opt
// sig = 0 -> from opt_needed
// surfflag = 0 -> volume
// surfflag = 1 -> surfaces
vector< vector< NekDouble> >  gradient_descent(Array<OneD, NekDouble> uhats, int sig, int surfflag )
{        
    // Gradient descent
    int temp_dim = dimension;

    if(surfflag == 1)
        temp_dim = dimension - 1;

    // Assert that d = 0
        
    Array<OneD, NekDouble> g(temp_dim);
    double inf = numeric_limits<double>::infinity();

    Array<OneD,NekDouble> xnew(temp_dim) ,x(temp_dim);
    NekDouble gprev = inf; // use inf
    Array<OneD, NekDouble> tmp(temp_dim);

    Array<OneD, Array<OneD, NekDouble> > xarr(temp_dim);
    vector< vector< NekDouble> > ret(temp_dim);
    Array<OneD, Array<OneD, NekDouble> > tc(temp_dim);
    StdExpansion *tempE;
    Array<OneD, Array<OneD, NekDouble> > eval(4);
    int sz;
    if(dimension == 2)
    {
        tempE= E;
        sz = testcoord2d[0].size();//pow(qZin.size(),2);;

        tc = testcoord2d;
        eval[0] = interioreval2d;// testcoord2d;
        eval[1] = interioreval2d0;// testcoord2d;
        eval[2] = interioreval2d1;// testcoord2d;
    }
    else if( (dimension == 3 && surfflag == 1))
    {
        sz = pow(qZin.size(),2);;

        tempE = E3quad;
        tc = testcoord2d;
        eval[0] = interioreval2d;// testcoord2d;
        eval[1] = interioreval2d0;// testcoord2d;
        eval[2] = interioreval2d1;// testcoord2d;
        
    }
    else //dim = 3 and surfflag is off
    {
        tc =  testcoord3d;
        sz = pow(qZin.size(),3);;
        eval[0] = interioreval3d;// testcoord2d;
        eval[1] = interioreval3d0;// testcoord2d;
        eval[2] = interioreval3d1;// testcoord2d;
        eval[3] = interioreval3d2;// testcoord2d;
        tempE = E;
    
    }

    //cout<<"\n sig="<<sig<<" dimension ="<<temp_dim<<"E3quad->GetNcoeffs()"<<E3quad->GetNcoeffs()<<"tc sz ="<<tc.size()<<"\n\n";

    // for(int k = 0; k <temp_dim+1; k++)
    // {
    //     eval[k] = Array<OneD, NekDouble>(tempE->GetNcoeffs());
    
    // }
    // for(int jj = 0; jj < tc[0].size(); ++jj)
    // {       
    //     for(int k = 0; k < temp_dim; k++)
    //     {
    //         x[k] = tc[k][jj];
    //     }
                
        Array<OneD, NekDouble> nullarr(0);
        Array<OneD, NekDouble> temp2(sz);
  

        if(sig==0)
        {

            Array<OneD, NekDouble> temp(uhats.size());
 
            for(int k = 0; k < sz; k++)
            {

                Vmath::Vmul(uhats.size(), &eval[0][k], sz, &uhats[0], 1, &temp[0], 1);          
                temp2[k] = Vmath::Vsum(temp.size(), temp, 1);
            }
            gprev = Vmath::Vmin(temp2.size(), temp2, 1);	
            
            
        }
        else
        {
            // xarr should ve vec<vec> >
            vector<vector<NekDouble> > xarr;
            Array<OneD, NekDouble> temp(uhats.size());
            
            for(int k = 0; k < temp_dim; k++)
            {   
                vector<NekDouble> row;
                 
                for(int p = 0; p < sz; p++)
                {
                    row.push_back(tc[k][p]);
                }
                xarr.push_back(row);
                    
            }
            
            
            pq(uhats,xarr, tempE,nullarr, temp2 );
            
            gprev = Vmath::Vmin(sz, temp2, 1);
        }

        // if(tmp[0]<gprev)
        // {
        //     for(int k =  0; k < temp_dim; k++)
        //     {
        //         xnew[k] = x[k];
                
        //     }
        //     gprev = tmp[0];
        // }
        //    }
    cout<<"\n gprev = "<<gprev<<"\n";    

    if(gprev < 0 && abs(gprev)>1e-12)
    {
        int n = 0;
        NekDouble epsl = 1.0;
       
        cout<<"\n GD: starting... gprev ="<<gprev<<" x = "<<xnew[0]<<" y = "<<xnew[1]<<" ";//<<" z = "<<xnew[2];
        if(temp_dim ==3)
        {
            cout<<" z = "<<xnew[2];
        }        
        NekDouble stepsize = 1e-3, iter = 1e5;
        Array<OneD, Array<OneD, NekDouble> > xastaa(temp_dim);

        Array<OneD, NekDouble> dereval(temp_dim);
        Array<OneD, NekDouble> xsave(xnew.size());
        for(int k = 0; k < temp_dim; k++)
            xsave[k] = xnew[k];
        //        NekDouble derevalx, derevaly;
        // loop for gradient descent 2D:
        
        bool res, res2 = true;
        
        
        while(n < iter)
        {
            epsl = 0.0;

            for(int p = 0; p < temp_dim; p++)
            {

                if(abs(xnew[p])<1e-12)
                    xnew[p] = 0.0;
   
                xastaa[p] = Array<OneD, NekDouble>(1,xnew[p]);
            }         

            tempE->PhysEvalBasisGrad(xastaa, eval[0], eval[1], eval[2], eval[3]);

            if(sig == 0) // from opt_needed
            {

                for(int p = 0; p < temp_dim; p++)
                {
                    
                    g[p] = eval[p+1][0];

                    epsl += abs(g[p]);
                }
            }
            else
            {
                for(int p = 0; p < temp_dim; p++)
                {
                    // cout<<"\n eval sz ="<<eval[0].size()<<"tempE->GetNcoeffs()="<<tempE->GetNcoeffs()<<"\n\n";
                     
                    derpq(uhats, dereval[p], eval[0], eval[p+1]);
                    // cout<<"\n eval sz ="<<eval[3].size()<<"tempE->GetNcoeffs()="<<tempE->GetNcoeffs()<<"\n\n";
      
                    g[p] = dereval[p];
                    epsl += abs(g[p]);
                }
            }
            for(int p = 0; p < temp_dim; p++)
            {
                
                x[p] = xnew[p] - stepsize*(g[p]);

            }
            res2 = true;
            //            cout<<"\n dim = "<<temp_dim<<"\n";
            for(int k = 0; k < temp_dim; k++)
            {   
                if(abs(x[k])-abs(xnew[k]) > 1e-10)
                {
                    res2 = res2 && false;
                    break;
                }
            }

            if( (std::isnan(epsl)) || (epsl <1e-10) || (epsl>1e3)|| (res2))
            {
                cout<<"\n break! n = "<<n<<"surfflag  ="<<surfflag<<"\n"<<" epsl = "<<epsl;//<<" x[0]="<<x[0]<<" x[1]="<<x[1]<<" xnew[0] = "<<xnew[0]<<"xnew[1] = "<<xnew[1]<<" x[0]-xnew[0] = "<<x[0]-xnew[0]<<" x[1]-xnew[1] = "<<x[1]-xnew[1];                ;                
                break;
            }
            for(int p = 0; p < temp_dim; p++)
            {
                xnew[p] = x[p];
            }            
            
            n = n+1;
            
            
        }

        if(numedges == 12)
            res = true;
        if(numedges == 4)
            res = true;
        else if(numedges == 3)
        {
            res = abs(x[0]+x[1])<0;//1e-12;
            cout<<"\n x[0]+x[1]="<<x[0]+x[1];
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
        if((res2 && res && n < iter &&(!std::isnan(epsl))))
        {
            for(int p = 0; p < temp_dim; p++)
            {
                ret[p].push_back( (x[p]));
            }
            cout<<"\n GD ends 1: n = "<<n<<" epsl = "<<epsl<<" x[0] ="<<x[0]<<" x[1]="<<x[1];
            if(temp_dim ==3)
            {
                cout<<"x[2]="<<x[2];
            }        

        }
    
        else
        {
            for(int p = 0; p < temp_dim; p++)
            {
                
                ret[p].push_back(xsave[p]);
            }
            cout<<"\n GD ends 2: res2 = "<<res2<<" res="<<res<<" n = "<<n<<" epsl = "<<epsl<<" xsave[0] ="<<xsave[0]<<" xsave[1]="<<xsave[1]<<"x[0] ="<<x[0]<<" x[1]="<<x[1];//<<"x[2]="<<xsave[2];
            if(temp_dim ==3)
            {
                cout<<"x[2]="<<xsave[2];
            }        
            
        }
    
    }

    return ret;
}
    
    
    
void derpq(Array<OneD, NekDouble> &uhats, NekDouble &ret,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd)
{
  int N = uhats.size();
  Array<OneD, NekDouble> temp(N);
  //    cout<<"\n  N = "<<N<<"Vxy sz = "<<Vxy.size()<<"\n";
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

  ret = v1-v2;

} 

// convert fn val at quad points -> uhats of der of function
// void derpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, Array<OneD,  NekDouble > >&ret,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> V3, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1)
// {
//     int totpts = E3->GetTotPoints();
//     int modes = E3->GetNcoeffs();
//     int uhatstot = uhats.size();
//     ret[0] = Array<OneD, NekDouble>(modes);
//     ret[1] = Array<OneD, NekDouble>(modes);
//     Array<OneD, NekDouble> temp(totpts), temp2(uhats.size()), temp3(uhats.size());

//     Array<OneD, NekDouble> pqevalx(totpts), pqevaly(totpts);
    
//     NekDouble v1, v2;
      

//     for(int i = 0; i<totpts; i++)
//     {

//         Vmath::Vmul(uhats.size(), &Vxyd0[0]+i, totpts, &uhats[0], 1, &temp2[0], 1);
//         v1  = Vmath::Vsum(uhats.size(), temp2, 1);  

//         Vmath::Vmul(uhats.size(), &Vxy[i], totpts, &Vxy[i], totpts, &temp2[0], 1);  
//         v2  = Vmath::Vsum(uhatstot, temp2, 1);  
//         v1 = v1*pow(v2, -1.0/2);

//         v2 = pow(v2, -3.0/2);


//         Vmath::Vmul(uhats.size(), &Vxy[0]+i, totpts, &Vxyd0[i], totpts, &temp2[0], 1);
//         v2 =  v2*Vmath::Vsum(uhatstot, temp2, 1);  

//         Vmath::Vmul(uhats.size(), &Vxy[0]+i, totpts, &uhats[0], 1, &temp2[0], 1);
   
//         v2 = v2* Vmath::Vsum(uhatstot, temp2, 1);  

        
//         pqevalx[i] = v1 - v2;
        

//         //******
//         Vmath::Vmul(uhats.size(), &Vxyd1[0]+i, totpts, &uhats[0], 1, &temp2[0], 1);
//         v1  = Vmath::Vsum(uhats.size(), temp2, 1);  

//         Vmath::Vmul(uhats.size(), &Vxy[i], totpts, &Vxy[i], totpts, &temp2[0], 1);  
//         v2  = Vmath::Vsum(uhatstot, temp2, 1);  
//         v1 = v1*pow(v2, -1.0/2);

//         v2 = pow(v2, -3.0/2);


//         Vmath::Vmul(uhats.size(), &Vxy[0]+i, totpts, &Vxyd1[i], totpts, &temp2[0], 1);
//         v2 =  v2*Vmath::Vsum(uhatstot, temp2, 1);  

//         Vmath::Vmul(uhats.size(), &Vxy[0]+i, totpts, &uhats[0], 1, &temp2[0], 1);
   
//         v2 = v2* Vmath::Vsum(uhatstot, temp2, 1);  

        
//         pqevaly[i] = v1 - v2;

//     } 
    

//     Array<OneD, NekDouble> holdx(pqevalx.size());
//     Array<OneD, NekDouble> holdy(pqevalx.size());
    
//     E3->MultiplyByStdQuadratureMetric(pqevalx, holdx);
//     E3->MultiplyByStdQuadratureMetric(pqevaly, holdy);

//     for(int i = 0; i<modes; i++)     
//     {
       
//         Vmath::Vmul(totpts, &V3[i], modes, &holdx[0], 1, &temp[0], 1);
//         ret[0][i] = Vmath::Vsum(totpts, temp, 1);  
//         Vmath::Vmul(totpts, &V3[i], modes, &holdy[0], 1, &temp[0], 1);
//         ret[1][i] = Vmath::Vsum(totpts, temp, 1);  
        
//     }
    
// }


int Opt_needed(Array<OneD, NekDouble> uhats)
{

    int dimension = E->GetShapeDimension();
    int totModes = uhats.size();
    vector<vector<  NekDouble> > roots;
    Array<OneD, Array<OneD, NekDouble> > edgeuhats  (numedges);
    Array<OneD, Array<OneD, NekDouble> > surfacederuhats  (numsurfaces);
    if( dimension > 1)
    {
        //find edge roots: (edges = 12 is hex)
  
        project_edges(uhats, edgeuhats);

        /*        for(int k = 0; k < edgeuhats.size(); k++)
        {
            for(int t = 0; t < edgeuhats[0].size(); t++)
            {
                cout<<" "<<edgeuhats[k][t]<<" ";
            }
            cout<<"\n";
        }*/
        
        if(dimension > 2)
        {
            //find surface roots: (surfaces = 6 is hex)
            project_surfaces(uhats, surfacederuhats, 0);
            // cout<<"\n surf uhats:\n";
            // for(int y = 0; y < surfacederuhats.size(); y++)
            // {
            //     for(int t = 0; t < surfacederuhats[0].size(); t++)
            //         {
            //             cout<<surfacederuhats[y][t]<< " ";
            //         }
            //     cout<<"\n";
            // }
            // exit(0);
        }

        //find roots of der of uhats, that is why d = 0
        roots = call_find_roots(uhats, 0, edgeuhats, surfacederuhats, 0);

    }
    else if(dimension == 1)
    {
        roots = call_find_roots(uhats, 1);
    }   

    /*    cout<<"\n optima in opt_needed: "<<roots.size()<<"\n";
          for(int pp = 0; pp < roots.size(); pp++)
          {
          for(int jj = 0; jj < roots[0].size(); jj++)
          {
          cout<<roots[pp][jj]<<" ";
          }
          cout<<"\n";
          }
          cout<<"\n";
    */  
    Array<OneD, Array<OneD,  NekDouble> > rootsarr(roots.size());
        
        
    for(int ii = 0; ii < roots.size(); ii++)
    {
        rootsarr[ii] = Array<OneD, NekDouble>(roots[ii].size());
        for(int jj = 0; jj < roots[0].size(); jj++)
        {
            rootsarr[ii][jj] = roots[ii][jj];
        }
                
    } 


    // evaluate ortho basis at roots
    // evalBasisRoots is flattened basis eval matrix

    int     evalsz = (roots[0].size())*totModes;
    Array<OneD, NekDouble> evalBasisRoots(evalsz);//, vals(roots[0].size()); 

    E->PhysEvalBasisGrad(rootsarr, evalBasisRoots, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
  
    Array<OneD, NekDouble> tmp(roots.size()),temp(totModes) ;

    NekDouble minv = numeric_limits<double>::infinity();
    int idx = 0;
    for(int k = 0; k < roots[0].size(); ++k)
    {
        NekDouble tempv;//, tempidx;
        Vmath::Vmul(totModes, &evalBasisRoots[k], roots[0].size(), &uhats[0], 1, &temp[0], 1);
        tempv = Vmath::Vsum(temp.size(), temp, 1);
        if(tempv < minv)
        {
            minv = tempv;
            idx = k;
        }
        cout<<"vals="<<tempv;
    }
         
 
        //    NekDouble minv = Vmath::Vmin(vals.size(), vals, 1);	
        //int idx = Vmath::Imin(vals.size(), vals, 1);

    cout<<"\n minv = "<<minv<<" at ";
    for(int k = 0; k < roots.size(); k++)
    {
        cout<<roots[k][idx]<<" ";
    }
    //    exit(0);
    if(minv < 0.0 && abs(minv)>1e-10)
    {
        return 1;
    }
    return 0;
}

void pq(
        Array<OneD,NekDouble> uhats,
        vector<vector< NekDouble> > roots,
        StdExpansion *tempE,
        Array<OneD,NekDouble> &pqevalxast,
        Array<OneD,NekDouble> &fvals
        )
{

    int N = uhats.size();
    int temp_dim = roots.size();
    
    Array<OneD, Array<OneD, NekDouble> > rootsarr(temp_dim); //dim
    
    for(int d = 0; d <temp_dim; d++)
    {
        rootsarr[d] = Array<OneD, NekDouble>(roots[d].size()); 
        for(int kk = 0; kk < roots[0].size(); kk++)
        {
            rootsarr[d][kk] = roots[d][kk];
            
        }

    }

  
    vector<NekDouble> Vsumsq;
    //cout<<"\n N = "<<N<<" roots[0] sz = "<<roots[0].size()<<" roots[0].size()*N="<<roots[0].size()*N<<"\n";
    // dies after this when o = 4 and p = 6

      
    Array<OneD, NekDouble> V1(roots[0].size()*N);
    //    cout<<"\n N = "<<N<<" roots[0] sz = "<<roots[0].size()<<" roots[0].size()*N="<<roots[0].size()*N<<"\n";
    
    tempE->PhysEvalBasisGrad(rootsarr, V1, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
    
    // V1 is flattened basis eval matrix
    
    
    // V2 = mat_mat_dot(V1,V1);
    // and
    // vector<double> Vsum = mat_sum_dim(V2,2);


    for(int i = 0; i < roots[0].size(); i++)
    {
        Array<OneD,NekDouble> w1(N);
        //cout<<"\n1 i ="<<i<<"\n\n";

        Vmath::Vmul(N, &V1[i], roots[0].size(), &uhats[0], 1, &w1[0], 1);
        fvals[i] = ( Vmath::Vsum(N, &w1[0], 1));
    
    }
    
    for( int i = 0; i < roots[0].size(); i++)
    {
        //        cout<<"\n1 i ="<<i<<"\n\n";
        Array<OneD,NekDouble> w1(N);
        //        Vmath::Smul(N, 0.0, w1, 1, w1, 1); 

        Vmath::Vmul(N, &V1[i], roots[0].size(), &V1[i], roots[0].size(), &w1[0], 1);

        Vsumsq.push_back(pow(Vmath::Vsum(N, &w1[0], 1),-0.5));
    }


    //    cts = vec_vec_dot(Vsumsq, ret);
    Vmath::Vmul(roots[0].size(), &Vsumsq[0], 1, &fvals[0], 1, &fvals[0], 1);
     
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
    vector< vector<NekDouble> > optima;
    
    //    int NC = 1; //number of constraints, only positivity for now
    tols.push_back(1e-11);
    
    Array<OneD, NekDouble>  pqvalcoords(dim+1), xastarr(dim), utemp(N1), wsp1;   
    

    Array<OneD, Array<OneD, NekDouble> > surfaceuhats  (numsurfaces);
    
    NekDouble pqval;
    cout<<"\n1\n\n";
    while (counter <= niter)
    {
        pqval = inf;
        utemp = d.back();
        Array<OneD, Array<OneD, NekDouble> > Pf(numedges);
        for(int k= 0; k < numedges; k++)
        {
            
            Pf[k] = Array<OneD, NekDouble>(3*(E->GetBasis(0)->GetNumModes()-1));        
        }

        //pq = @(xx) cfun(xx) .* (constraints(1,jj).*Vc(xx)*uhats);  
        project_edges(utemp, Pf, 1);
        
        cout<<"\n1 counter="<<counter<<"\n\n";
        if(dimension > 2)
        {
            //find surface roots: (surfaces = 12 is hex)
            cout<<"\n before surfcounter = "<<counter<<"Pf sz ="<<Pf.size()<<" "<<Pf[0].size()<<" ";
        cout<<"\n1 counter="<<counter<<"\n\n";
        
            project_surfaces(utemp, surfaceuhats, 1);

            cout<<"\n after surf counter = "<<counter<<"Pf sz ="<<Pf.size()<<" "<<Pf[0].size()<<" ";
        
        }
        // cout<<"Pf sz in do_opt:"<<Pf[0].size()<<"\n";
        
        // for(int p = 0; p<Pf.size(); p++)
        // {
        //     cout<<" "<<Pf[p].size()<<" ";
        // }
      
        optima = (call_find_roots(utemp, 0, Pf, surfaceuhats, 1));
        wsp1=Array<OneD, NekDouble>(optima[0].size());
        
        cout<<"\n optima in do_opt: "<<optima.size()<<"\n";
        for(int pp = 0; pp < optima.size(); pp++)
        {
            for(int jj = 0; jj < optima[0].size(); jj++)
            {
                cout<<optima[pp][jj]<<" ";
            }
            cout<<"\n";
        }
        cout<<"\n";

        pq(utemp, optima, E, pqvalcoords, wsp1);

        cout<<"\n fvals: ";
         for(int kk = 0 ; kk<wsp1.size(); kk++)
             cout<<wsp1[kk]<<" ";
         cout<<"\n";//exit(0);

        if (pqvalcoords[0] < pqval)
        {
            for(int k = 0; k  < dimension; k++)
            {
                xastarr[k] = pqvalcoords[k+1];
            }
            pqval = pqvalcoords[0];
        }
        
        //}
        if(xastarr.size() == 3)
        {
            cout<<"\n at counter="<<counter<<" min val="<<pqval<<" xast ="<<xastarr[0]<<" "<<xastarr[1]<<" "<<xastarr[2]<<
                "\n\n";
        }
        else if(xastarr.size() == 2)
        {
            cout<<"\n at counter="<<counter<<" min val="<<pqval<<" xast ="<<xastarr[0]<<" "<<xastarr[1]<<
                "\n\n";
        }
        else if(xastarr.size() == 1)
        {
            cout<<"\n at counter="<<counter<<" min val="<<pqval<<" xast ="<<xastarr[0]<<"\n\n";
        
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
    cout<<"sphere_rotation took "<<counter<<"  iterations\n ";
    uhats = d.back();
    /*    cout<<"\n sphere_rot returns:\n";
    for(int k = 0; k < uhats.size(); k++)
        cout<<" "<<uhats[k];
        cout<<"\n";*/
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
