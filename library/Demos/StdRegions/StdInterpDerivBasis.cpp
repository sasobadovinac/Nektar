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
    int nTot = nCoeffs * nPts, dimension = E->GetShapeDimension();

    Array<OneD, Array<OneD, NekDouble>> coords = demo.GetCoords(E);
    Array<OneD, NekDouble> sol(nTot), phys(nTot), tmpIn(dimension);
    Array<OneD, NekDouble> sol1(nTot), phys1(nTot), tmpIn1(dimension);
    Array<OneD, NekDouble> sol2(nTot), phys2(nTot), tmpIn2(dimension);

    // For each mode, we follow two approaches:
    //
    // 1) Evaluate the basis function at each quadrature point using the
    //    StdExpansion::PhysEvaluateBasis function.
    // 2) Evaluate the basis function at all quadrature points using FillMode.
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
                phys2[k * nPts + i] = E->PhysEvaluatedzBasis(tmpIn, k);
            }
            if(dimension>1)
            {
                phys1[k * nPts + i] = E->PhysEvaluatedyBasis(tmpIn, k);
            }
            if(dimension>0)
            {
                phys[k * nPts + i] = E->PhysEvaluatedxBasis(tmpIn, k);                
            }
        }

        // Fill the 'solution' field with each of the partial der of modes wrt x using FillModedx.
        if(dimension>2)
        {
            Array<OneD, NekDouble> tmp = sol2 + k * nPts;
            E->FillModedz(k, tmp);
        }
        if(dimension>1)
        {
            Array<OneD, NekDouble> tmp = sol1 + k * nPts;
            E->FillModedy(k, tmp);
        }
        // Fill the 'solution' field with each of the partial der of modes wrt y using FillModedx.
        if(dimension>0)
        {
            Array<OneD, NekDouble> tmp = sol + k * nPts;
            E->FillModedx(k, tmp);
            
        }
    }

    cout<<"\n true sol:";
    for (int i = 0; i<sol.size(); i++)
        cout<<" "<<sol[i]<<" ";
    cout<<"\n";
    cout<<"\n phys sol:";
    for (int i = 0; i<phys.size(); i++)
        cout<<" "<<phys[i]<<" ";
    cout<<"\n";


    //cout << "L infinity error : " << scientific << E->Linf(phys, sol) + E->Linf(phys1, sol1)+ E->Linf(phys2, sol2) << endl;
    //cout << "L 2 error        : " << scientific << E->L2  (phys, sol) +  E->L2  (phys1, sol1) +  E->L2  (phys2, sol2) << endl;
        cout << "L infinity error : " << scientific <<  E->Linf(phys1, sol1)<< endl;
    cout << "L 2 error        : " << scientific << E->L2  (phys1, sol1) << endl;

    
    return 0;
}
