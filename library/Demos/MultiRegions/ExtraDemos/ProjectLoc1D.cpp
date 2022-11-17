///////////////////////////////////////////////////////////////////////////////
//
// File: ProjectLoc1D.cpp
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

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ExpList.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

// compile using Builds/Demos/StdRegions -> make DEBUG=1 ProjectLoc1D

// This routine projects a polynomial which has energy in all mdoes of
// the expansions and report an error.

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(argc, argv);

    MultiRegions::ExpListSharedPtr Exp, Sol;
    int i, j;
    int nq;
    int coordim;
    Array<OneD, NekDouble> sol;
    Array<OneD, NekDouble> xc0, xc1, xc2;

    // read in mesh
    SpatialDomains::MeshGraphSharedPtr graph1D =
        SpatialDomains::MeshGraph::Read(vSession);

    // Define Expansion
    const SpatialDomains::ExpansionInfoMap &expansions =
        graph1D->GetExpansionInfos();
    LibUtilities::BasisKey bkey0 =
        expansions.begin()->second->m_basisKeyVector[0];
    int nmodes = bkey0.GetNumModes();

    Exp = MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(vSession,
                                                                  graph1D);

    //----------------------------------------------
    // Define solution to be projected
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetTotPoints();

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
        for (j = 0; j < nmodes; ++j)
        {
            sol[i] += pow(xc0[i], j);
            sol[i] += pow(xc1[i], j);
            sol[i] += pow(xc2[i], j);
        }
    }

    //---------------------------------------------
    // Set up ExpList containing the solution
    Sol = MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(*Exp);
    Sol->SetPhys(sol);
    //---------------------------------------------

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
    if (vSession->GetComm()->GetRank() == 0)
    {
        cout << "L infinity error: " << Exp->Linf(Sol->GetPhys()) << endl;
        cout << "L 2 error:        " << Exp->L2(Sol->GetPhys()) << endl;
    }
    //--------------------------------------------

    vSession->Finalise();

    return 0;
}
