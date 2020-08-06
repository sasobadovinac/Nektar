///////////////////////////////////////////////////////////////////////////////
//
// File: StdInterpBasis.cpp
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
// Description: Demo for testing functionality of PhysEvaluateBasis
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <LibUtilities/BasicUtils/Timer.h>

#include "StdDemoSupport.hpp"

int main(int argc, char *argv[])
{
    DemoSupport demo;
    demo.ParseArguments(argc, argv);
    StdExpansion *E = demo.CreateStdExpansion();

    if (E == nullptr)
    {
        return 1;
    }

    int nCoeffs = E->GetNcoeffs(), nPts = E->GetTotPoints();
    int dimension = E->GetShapeDimension();
    vector<string> &ptypes = demo.GetPointsType();
    for (int i = 0; i < dimension; ++i)
    {
        ptypes[i] = "GaussGaussChebyshev";
    }

    cout<<"\nPts = "<<nPts;
    cout<<"\n nCoeffs="<<nCoeffs;
    Array<OneD, Array<OneD, NekDouble>> coords = demo.GetCoords(E);

    Array<OneD, NekDouble> sol(nPts*nCoeffs), soll(nPts), out_eval(nPts*nCoeffs), phys(nPts*nCoeffs);

    Array<OneD, NekDouble> sol1(nPts*nCoeffs), phys1(nCoeffs*nPts);
    Array<OneD, NekDouble> sol2(nPts*nCoeffs), phys2(nPts*nCoeffs);
    NekDouble errL2 = 0, errLinf = 0;
    
    // For each mode, we follow two approaches:
    //
    // 1) Evaluate the der of basis function at each quadrature point using the
    //    StdExpansion::PhysEvaluateBasisdx function.
    // 2) Evaluate the der of basis function at all quadrature points using bary interp of FillModedx and FillModedy.
    //
    // These are then compared to ensure they give the same result.
    
    //for (int k = 0; k < nCoeffs; ++k)
    //{
    // Evaluate each mode at the quadrature points.
    //for (int i = 0; i < nPts; ++i)
    //{
    /*          for (int d = 0; d < dimension; ++d)
                {
                tmpIn[d] = coords[d][i];
                }
    */            
    if(dimension>2)
    {
        //phys2[i] = E->PhysEvaluatedzBasis(tmpIn, k);
        //sol2[i] = E->PhysEvaluatedzBasisBary(tmpIn, k);
    }
    else if(dimension>1)

    {
        E->PhysEvalBasisGradFast(coords,  out_eval, phys, phys1,  NullNekDouble1DArray);
                
        for (int k = 0; k < nCoeffs; ++k)
        {
            // Fill the 'solution' field with each of the modes using FillMode.
            Array<OneD, NekDouble> hold1(nPts);
            Array<OneD, NekDouble> hold2(nPts);
            E->FillMode(k, soll);
                    
            E->PhysDeriv(soll, hold1, hold2, NullNekDouble1DArray);       
            Vmath::Vcopy(nPts, &hold1[0], 1, &sol[k*nPts], 1);
            Vmath::Vcopy(nPts, &hold2[0], 1, &sol1[k*nPts], 1);
        }
                    
    }
    else if(dimension>0)
    {
                E->PhysEvalBasisGradFast(coords, out_eval, phys, NullNekDouble1DArray,  NullNekDouble1DArray);

        for (int k = 0; k < nCoeffs; ++k)
        {
            // Fill the 'solution' field with each of the modes using FillMode.

            Array<OneD, NekDouble> hold1(nPts);
            E->FillMode(k, soll);
            E->PhysDeriv(soll,  hold1, NullNekDouble1DArray, NullNekDouble1DArray);       
            // check if soll = out_eval
                    
            Vmath::Vcopy(nPts, &hold1[0], 1, &sol[k*nPts], 1);  
        }
    }
        
        
    cout<<"\n phys = \n";
    for(int i =  0; i<phys.size(); i++)
        cout<<phys[i]<<" ";
    cout<<" \n\n";


    cout<<"\n sol = \n";
    for(int i =  0; i<sol.size(); i++)
        cout<<sol[i]<<" ";
    cout<<" \n\n";

    cout<<"\n phys1 = \n";
    for(int i =  0; i<phys1.size(); i++)
        cout<<phys1[i]<<" ";
    cout<<" \n\n";
            


    cout<<"\n sol1 = \n";
    for(int i =  0; i<sol1.size(); i++)
        cout<<sol1[i]<<" ";
    cout<<" \n\n";
            
        
    errL2 += E->L2(phys1, sol1);
    errLinf += E->Linf(phys, sol)+ E->Linf(phys1, sol1);
        
        
    cout << "L infinity error : " << scientific << errL2 << endl;
    cout << "L 2 error        : " << scientific << errLinf << endl;
    

    
    return 0;
}
