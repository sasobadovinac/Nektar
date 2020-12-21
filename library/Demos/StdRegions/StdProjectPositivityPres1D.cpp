///////////////////////////////////////////////////////////////////////////////
//
// File: StdProjectPositivityPres2D.cpp
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
vector<vector<  NekDouble> > find_roots(Array<OneD, NekDouble> &uhats, int d = 0);


//declare Opt_needed
int Opt_needed(Array<OneD, NekDouble> uhats);

// uhatpqd stuff: 
// called by project_edges if d = 1 is passed in project_edges
void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret, int modes,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd, Array<OneD, NekDouble> V3, Array<OneD, NekDouble> qw);

// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble> uhats,Array<OneD, NekDouble> &ret , int d = 0);


// void deruhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble >&ret,  Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0, Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2 = NullNekDouble1DArray, int surfflag = 0 );


// Colleague matrix

Array<OneD, Array<OneD, NekDouble> > C;
StdExpansion *E;
StdExpansion *E3seg;
Array<OneD, NekDouble> qx;
Array<OneD, NekDouble> qw;
Array<OneD, Array<OneD, NekDouble> > qZinarr;

Array<OneD, NekDouble> V;
Array<OneD, NekDouble> Vd;
Array<OneD, NekDouble> V3;

int dimension;
            
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
    if(dimension > 1)
    {
        cout<<"\n dimension should be 1, try using StdProjectPositivityPres2D or StdProjectPositivityPres3D for other dimensions\n\n";
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

    C = demo.formConf((dimension)*pow(3*order[0]+1,2));
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
            x[i] = ((1-0)*(x[i] + 1))/(2) + 0;
            y[i] = ((1-0)*(y[i] + 1))/(2) + 0; 
            z[i] = ((1-0)*(z[i] + 1))/(2) + 0; 
            sol[i] = pow(sin(x[i]*2),2)+pow((cos(y[i]*2)),2)+pow(sin(z[i]*2),2)+pow((y[i]*x[i]*z[i]),2)- 0.2;
        
        }

        else if(dimension ==2) //only quad
        {
            x[i] = ((1-0)*(x[i] + 1))/(2) + 0;
            y[i] = ((1-0)*(y[i] + 1))/(2) + 0; 
            sol[i] = pow(sin(x[i]*2),2)+pow((cos(y[i]*2)),2)+pow((y[i]*x[i]),2)- 0.5;
        
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
        Array<OneD, Array<OneD, NekDouble> > qxarr(dimension);
        qx = E->GetBasis(0)->GetZ();
        qw = E->GetBasis(0)->GetW();            
        vector<LibUtilities::BasisKey> tmpBasisKey = demo.GetBasisKey();       
        LibUtilities::PointsKey pkeycheb(E->GetBasis(0)->GetZ().size(), pointsTypeCheb);

        LibUtilities::BasisKey bkeycheb(LibUtilities::eChebyshev, 3*(E->GetBasis(0)->GetNumModes()-1),  pkeycheb);
        
        E3seg = new StdSegExp(bkeycheb);
        
        if (E3seg == nullptr)
        {
            return 1;
        }
        // V3 is required by edge root finding (even for 1D)
        V3 = Array<OneD, NekDouble>(E3seg->GetTotPoints()*E3seg->GetNcoeffs());

        qxarr[0] = E->GetBasis(0)->GetZ();
  
        E3seg->PhysEvalBasisGrad(qxarr, V3, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray );
        
            qZinarr = Array<OneD, Array<OneD, NekDouble> >(dimension);
            qZinarr[0] = E->GetBasis(0)->GetZ();
            Vd = Array<OneD, NekDouble>(qx.size()*coeffs.size());
            V = Array<OneD, NekDouble>(qx.size()*coeffs.size());
            
            E->PhysEvalBasisGrad(qZinarr, V, Vd, NullNekDouble1DArray,  NullNekDouble1DArray);
           
        
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
                cout<<"\n pass\n\n";
                phys=sol;
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



// when d = 1, its the uhatpqd case from do_optimize
void project_edges( Array<OneD, NekDouble>uhats,    Array<OneD, NekDouble> &ret , int d)
{
    int modes = E->GetNcoeffs();
    if(d == 1) // so project der to 3*N space
    {
        edgederpquhats(uhats, ret, modes,V, Vd, V3, qw);
        
    }
    else
    {
        ret = Array<OneD, NekDouble>(modes);        

        ret = uhats;
    }

    
   
}


// void deruhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  Array<OneD, NekDouble> Vxyz, Array<OneD, NekDouble> Vdxyz, Array<OneD, NekDouble> Vxdyz, Array<OneD, NekDouble> Vxydz, int surfflag)
// {
//     boost::ignore_unused(Vxydz, Vxyz, surfflag);

//     int totpts = E->GetTotPoints();
//     Array<OneD, Array<OneD,  NekDouble> >temp(dimension);;
//     int modes = uhats.size();
//     Array<OneD, NekDouble> vals(totpts), hold(totpts);    
//     for(int k = 0; k < dimension; k++)
//     {
//         temp[k ] = Array<OneD, NekDouble> (modes); 
//     }
//     NekDouble v1;
//     for(int k = 0; k < totpts; k++)
//     {    
//         Vmath::Vmul(uhats.size(), &Vdxyz[k], totpts, &uhats[0], 1, &temp[0][0], 1);
//         v1  = Vmath::Vsum(modes, temp[0], 1);
//         Vmath::Vmul(uhats.size(), &Vxdyz[k], totpts, &uhats[0], 1, &temp[1][0], 1);
//         v1  = v1 + Vmath::Vsum(modes, temp[1], 1);  
        
        
//     }    
    
//     E->MultiplyByQuadratureMetric(vals,
//                                   temp[0]);

//     for(int i = 0; i<ret.size(); i++)
//     {

//         Vmath::Vmul(totpts, &Vxyz[i], ret.size(), &temp[0][0], 1, &hold[0], 1);
//         ret[i] = Vmath::Vsum(vals.size(), hold, 1);
//     }

// }



void edgederpquhats(Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret,  int modes, Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd, Array<OneD, NekDouble> V3, Array<OneD, NekDouble> qw)
{
    int uhatstot = uhats.size();
    int totpts = E3seg->GetTotPoints();
        
    Array<OneD, NekDouble> temp(totpts), temp2(modes), temp3(uhatstot);
    Array<OneD, NekDouble> pqeval(totpts);
    NekDouble v1, v2;

    for(int i = 0; i<totpts; i++)
    {
        Vmath::Vmul(modes, &Vxy[0]+i, totpts, &Vxyd[0]+i, totpts, &temp2[0], 1);
        v1  = Vmath::Vsum(modes, temp2, 1); 
        Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);

        v2  = Vmath::Vsum(uhatstot, temp3, 1);  

        v1 = v2*v1;

        // At this point,                 

        // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

        // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
        Vmath::Vmul(uhatstot, &Vxy[0]+i, totpts, &Vxy[0]+i, totpts, &temp2[0], 1);

 
        
        v2  = Vmath::Vsum(uhatstot, temp2, 1);  
 
        Vmath::Vmul(uhats.size(), &Vxyd[i], totpts, &uhats[0], 1, &temp2[0], 1);

        v1= v1 - v2*Vmath::Vsum(uhats.size(), temp2, 1);
 
        pqeval[i] = v1;
 
    }

    Vmath::Vmul(totpts, qw, 1, pqeval, 1, pqeval, 1);
    

    for(int i = 0; i<ret.size(); i++)
    {
    
        Vmath::Vmul(totpts, &V3[i], ret.size(), &pqeval[0], 1, &temp[0], 1);
        ret[i] = Vmath::Vsum(totpts, temp, 1);  
        
    }
    
}



// d=1 :: derivatives
// sig = 0 -> called by opt_needed
// sig = 1 -> called by sphere_rotation

vector<vector<NekDouble> > find_roots(Array<OneD, NekDouble> &uhats, int d)
{
    vector<vector< NekDouble> > ret(dimension);
    //Confederate matrix approach
    int N = uhats.size();

    vector<NekDouble> uhatsmon;

    while(abs(uhats[N-1])<1e-8)
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
        
        for(int kk = 0; kk <EIG_R.size(); kk++)
        {
            ret[0].push_back( EIG_R[kk] );
        }

    return ret;    
    
}
    
    

int Opt_needed(Array<OneD, NekDouble> uhats)
{

    int totModes = uhats.size();
    vector<vector<  NekDouble> > roots;
     Array<OneD, NekDouble>  edgeuhats;
    project_edges(uhats, edgeuhats);

    roots = find_roots(edgeuhats,1);

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
    Array<OneD, NekDouble> evalBasisRoots(evalsz); 

    E->PhysEvalBasisGrad(rootsarr, evalBasisRoots, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
  
    Array<OneD, NekDouble> tmp(roots.size());
    
    Array<OneD, NekDouble> vals(roots[0].size()), temp(totModes);
    for(int k = 0; k < vals.size(); ++k)
    {
        Vmath::Vmul(totModes, &evalBasisRoots[k], vals.size(), &uhats[0], 1, &temp[0], 1);
        vals[k] = Vmath::Vsum(temp.size(), temp, 1);
    }
 
    NekDouble minv = Vmath::Vmin(vals.size(), vals, 1);	
    int idx = Vmath::Imin(vals.size(), vals, 1);

    cout<<"\n minv = "<<minv<<" at ";
    for(int k = 0; k < roots.size(); k++)
    {
        cout<<roots[k][idx]<<" ";
    }
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
    Array<OneD,NekDouble> wsp1(N);
    
    Array<OneD, NekDouble> V1(roots[0].size()*N);
    
    tempE->PhysEvalBasisGrad(rootsarr, V1, NullNekDouble1DArray, NullNekDouble1DArray, NullNekDouble1DArray);
    
    // V1 is flattened basis eval matrix
    
    
    // V2 = mat_mat_dot(V1,V1);
    // and
    // vector<double> Vsum = mat_sum_dim(V2,2);

    for( int i = 0; i < roots[0].size(); i++)
    {
        Vmath::Vmul(N, &V1[i], roots[0].size(), &V1[i], roots[0].size(), &wsp1[0], 1);

        Vsumsq.push_back(pow(Vmath::Vsum(N, &wsp1[0], 1),-0.5));
    }

    for(int i = 0; i < roots[0].size(); i++)
    {
        Vmath::Vmul(N, &V1[i], roots[0].size(), &uhats[0], 1, &wsp1[0], 1);
        fvals[i] = ( Vmath::Vsum(N, &wsp1[0], 1));
    
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
    
    NekDouble pqval;
  
    while (counter <= niter)
    {
        pqval = inf;
        utemp = d.back();
        Array<OneD, NekDouble>  Pf(3*(E->GetBasis(0)->GetNumModes()-1));        
        
        //pq = @(xx) cfun(xx) .* (constraints(1,jj).*Vc(xx)*uhats);  
        project_edges(utemp, Pf, 1);
        optima = (find_roots(Pf));
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
        if (pqvalcoords[0] < pqval)
        {
            for(int k = 0; k  < dimension; k++)
            {
                xastarr[k] = pqvalcoords[k+1];
            }
            pqval = pqvalcoords[0];
        }
        
        cout<<"\n at counter="<<counter<<" min val="<<pqval<<" xast ="<<xastarr[0]<<"\n\n";
        
        
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
