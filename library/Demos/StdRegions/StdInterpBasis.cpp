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
    Array<OneD, NekDouble> sol(nTot), phys(nTot), physOut(nTot), tmpIn(dimension);

    // For each mode, we follow two approaches:
    //
    // 1) Evaluate the basis function at each quadrature point using the
    //    StdExpansion::PhysEvaluateBasis function.
    // 2) Evaluate the basis function at all quadrature points using FillMode.
    //
    // These are then compared to ensure they give the same result.
    
    // Evaluate each mode at the quadrature points.

    auto storage = E->GetPhysEvaluateStorage();

    Array<OneD, NekDouble> tmp;
    for (int i = 0; i < nCoeffs; ++i)
    {
        tmp = E->PhysEvaluateBasis(coords, storage, i);
        Vmath::Vcopy(nPts, &tmp[0], 1, &phys[nPts*i], 1);
    }

    // Another approach: Use Nektar++'s approach treating 
    // the whole FillMode on all quad points as a function 
    // evaluation on domain. Do not leverage the multiplicative 
    // separability of basis definitions in each individual direction:
    Array<OneD, NekDouble> ar1(nPts);
    for (int i = 0; i < nCoeffs; ++i)
    {
        Vmath::Vcopy(nPts, &storage[0][i*nPts],1, &ar1[0], 1);
            
        for (int k = 0; k < nPts; ++k)
        {
            for (int d = 0; d < dimension; ++d)
            {
                tmpIn[d] = coords[d][k];
            }

            sol[i*nPts + k] = E->PhysEvaluate(tmpIn, ar1);
        }
    }
    
    Array<OneD, NekDouble> physpts (nPts);
    Array<OneD, NekDouble> solpts (nPts);
    NekDouble errL2 = 0, errLinf = 0;
    // Separate modes 0 to nCoeffs
    for( int i = 0 ; i < nCoeffs; i++)
    {
        Vmath::Vcopy(nPts, &sol[nPts*i], 1, &solpts[0], 1 );
        Vmath::Vcopy(nPts, &phys[nPts*i], 1, &physpts[0], 1 );

        errL2 += E->L2(solpts, physpts);
        errLinf += E->Linf(solpts, physpts);
    }

    cout << "L infinity error : " << scientific << errLinf << endl;
    cout << "L 2 error        : " << scientific << errL2 << endl;

    return 0;
}
