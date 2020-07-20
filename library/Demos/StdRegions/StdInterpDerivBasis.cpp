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

    Array<OneD, Array<OneD, NekDouble>> coords = demo.GetCoords(E);
    Array<OneD, NekDouble> sol(nPts), phys(nPts), tmpIn(dimension);
    Array<OneD, NekDouble> sol1(nPts), phys1(nPts), tmpIn1(dimension);
    Array<OneD, NekDouble> sol2(nPts), phys2(nPts), tmpIn2(dimension);
    NekDouble errL2 = 0, errLinf = 0;
    
    // For each mode, we follow two approaches:
    //
    // 1) Evaluate the der of basis function at each quadrature point using the
    //    StdExpansion::PhysEvaluateBasisdx function.
    // 2) Evaluate the der of basis function at all quadrature points using bary interp of FillModedx and FillModedy.
    //
    // These are then compared to ensure they give the same result.
    
    for (int k = 0; k < nCoeffs; ++k)
    {
        // Evaluate each mode at the quadrature points.
        for (int i = 0; i < nPts; ++i)
        {
            for (int d = 0; d < dimension; ++d)
            {
                tmpIn[d] = coords[d][i];
            }
            
            if(dimension>2)
            {
                phys2[i] = E->PhysEvaluatedzBasis(tmpIn, k);
                sol2[i] = E->PhysEvaluatedzBasisBary(tmpIn, k);
            }
            if(dimension>1)
            {
                phys1[i] = E->PhysEvaluatedyBasis(tmpIn, k);
                sol1[i] = E->PhysEvaluatedyBasisBary(tmpIn, k);
                
            }
            if(dimension>0)
            {
                phys[i] = E->PhysEvaluatedxBasis(tmpIn, k);       
                sol[i] = E->PhysEvaluatedxBasisBary(tmpIn, k);       

         
            }
        }
        /*
        cout<<"\n phys = \n";
        for(int i =  0; i<phys.size(); i++)
            cout<<phys[i]<<" ";
        cout<<" \n\n";

        cout<<"\n phys1 = \n";
        for(int i =  0; i<phys1.size(); i++)
            cout<<phys1[i]<<" ";
        cout<<" \n\n";

        cout<<"\n sol = \n";
        for(int i =  0; i<sol.size(); i++)
            cout<<sol[i]<<" ";
        cout<<" \n\n";

        cout<<"\n sol1 = \n";
        for(int i =  0; i<sol1.size(); i++)
            cout<<sol1[i]<<" ";
        cout<<" \n\n";

        exit(0);*/
            

        errL2 += E->L2(phys, sol)+ E->L2(phys1, sol1)+ E->L2(phys2, sol2);
        errLinf += E->Linf(phys, sol)+ E->Linf(phys1, sol1)+  E->Linf(phys2, sol2);
    }
        
    cout << "L infinity error : " << scientific << errL2 << endl;
    cout << "L 2 error        : " << scientific << errLinf << endl;
    

    
    return 0;
}
