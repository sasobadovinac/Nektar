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
    //for (int i = 0; i < dimension; ++i)
        //{
        ptypes[0] = "Modified_A";
        ptypes[1] = "Modified_A";
        ptypes[2] = "Modified_B";
    //}

    Array<OneD, Array<OneD, NekDouble>> coords = demo.GetCoords(E);

    Array<OneD, NekDouble> sol(nPts*nCoeffs), soll(nPts), out_eval(nPts*nCoeffs), phys(nPts*nCoeffs);

    Array<OneD, NekDouble> sol1(nPts*nCoeffs), phys1(nCoeffs*nPts);
    Array<OneD, NekDouble> sol2(nPts*nCoeffs), phys2(nPts*nCoeffs);
    NekDouble errL2 = 0, errLinf = 0;
    
    if(dimension>2)
    {
     
        E->PhysEvalBasisGradFast(coords, out_eval, phys, phys1, phys2);
        E->PhysEvalBasisGrad(coords, out_eval, sol, sol1, sol2);

    }
    else if(dimension>1)

    {
       E->PhysEvalBasisGradFast(coords,  out_eval, phys, phys1, NullNekDouble1DArray); 
       E->PhysEvalBasisGrad(coords,  out_eval, sol,sol1,  NullNekDouble1DArray);
                    
    }
    else if(dimension>0)
    {
        E->PhysEvalBasisGradFast(coords, out_eval, phys, NullNekDouble1DArray,  NullNekDouble1DArray);
        E->PhysEvalBasisGrad(coords,  out_eval, sol,  NullNekDouble1DArray,  NullNekDouble1DArray);

        // for (int k = 0; k < nCoeffs; ++k)
        // {
        //     // Fill the 'solution' field with each of the modes using FillMode.

        //     Array<OneD, NekDouble> hold1(nPts);
        //     E->FillMode(k, soll);
        //     E->PhysDeriv(soll,  hold1, NullNekDouble1DArray, NullNekDouble1DArray);       
        //     // check if soll = out_eval
                    
        //     Vmath::Vcopy(nPts, &hold1[0], 1, &sol[k*nPts], 1);  
        // }
    }
        
        
    errL2 += E->L2(phys1, sol1)+ E->L2(phys2, sol2)+E->L2(phys, sol);
    errLinf += E->Linf(phys, sol)+ E->Linf(phys1, sol1)+E->Linf(phys2,sol2);
        
        
    cout << "L infinity error : " << scientific << errL2 << endl;
    cout << "L 2 error        : " << scientific << errLinf << endl;
    

    return 0;
}
