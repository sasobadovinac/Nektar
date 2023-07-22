///////////////////////////////////////////////////////////////////////////////
//
// File: ProjectCont1D.cpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ContField.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

// compile using Builds/Demos/StdRegions -> make DEBUG=1 ProjectCont1D

// This routine projects a which has energy in all mdoes of the
// expansions and report an error.

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(argc, argv);

    MultiRegions::ContFieldSharedPtr Exp, Sol;

    int i, j;
    int order, nq;
    int coordim;
    Array<OneD, NekDouble> sol;
    Array<OneD, NekDouble> xc0, xc1, xc2;

    if (argc != 2)
    {
        fprintf(stderr, "Usage: ProjectCont1D mesh \n");
        exit(1);
    }

    //----------------------------------------------
    // read the problem parameters from input file
    SpatialDomains::MeshGraphSharedPtr graph1D =
        SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    const SpatialDomains::ExpansionInfoMap &expansions =
        graph1D->GetExpansionInfos();
    LibUtilities::BasisKey bkey0 =
        expansions.begin()->second->m_basisKeyVector[0];
    int nmodes = bkey0.GetNumModes();
    cout << "Solving 1D Continuous Projection" << endl;
    cout << "    Expansion  : ("
         << LibUtilities::BasisTypeMap[bkey0.GetBasisType()] << ")" << endl;
    cout << "    No. modes  : " << nmodes << endl;
    cout << endl;
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    Exp = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
        vSession, graph1D, vSession->GetVariable(0));
    //----------------------------------------------

    //----------------------------------------------
    // Define solution to be projected
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetTotPoints();
    order   = Exp->GetExp(0)->GetNcoeffs();

    // define coordinates and solution
    sol = Array<OneD, NekDouble>(nq);

    xc0 = Array<OneD, NekDouble>(nq);
    xc1 = Array<OneD, NekDouble>(nq);
    xc2 = Array<OneD, NekDouble>(nq);

    switch (coordim)
    {
        case 1:
            Exp->GetCoords(xc0);
            Vmath::Zero(nq, &xc1[0], 1);
            Vmath::Zero(nq, &xc2[0], 1);
            break;
        case 2:
            Exp->GetCoords(xc0, xc1);
            Vmath::Zero(nq, &xc2[0], 1);
            break;
        case 3:
            Exp->GetCoords(xc0, xc1, xc2);
            break;
    }

    for (i = 0; i < nq; ++i)
    {
        sol[i] = 0.0;
        for (j = 0; j < order; ++j)
        {
            sol[i] += pow(xc0[i], j);
            sol[i] += pow(xc1[i], j);
            sol[i] += pow(xc2[i], j);
        }
    }
    //----------------------------------------------

    //----------------------------------------------
    // Setup Temporary expansion and put in solution
    Sol = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(*Exp);
    Sol->SetPhys(sol);
    //----------------------------------------------

    //---------------------------------------------
    // Project onto Expansion
    Exp->FwdTrans(Sol->GetPhys(), Exp->UpdateCoeffs());
    //---------------------------------------------

    //-------------------------------------------
    // Backward Transform Solution to get projected values
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //-------------------------------------------

    //--------------------------------------------
    // Calculate L_inf error
    cout << "L infinity error: " << Exp->Linf(Sol->GetPhys()) << endl;
    cout << "L 2 error:        " << Exp->L2(Sol->GetPhys()) << endl;
    //--------------------------------------------

    vSession->Finalise();

    return 0;
}
