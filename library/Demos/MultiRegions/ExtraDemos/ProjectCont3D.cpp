///////////////////////////////////////////////////////////////////////////////
//
// File: ProjectCont3D.cpp
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
#include <MultiRegions/MultiRegions.hpp>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

// This routine projects a polynomial which has energy in all mdoes of
// the expansions and report an error.

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(argc, argv);

    MultiRegions::ContFieldSharedPtr Exp, Fce;
    int i, j, nq, coordim;
    Array<OneD, NekDouble> fce;
    Array<OneD, NekDouble> xc0, xc1, xc2;

    if (argc != 2)
    {
        fprintf(stderr, "Usage: ProjectCont3D  meshfile \n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph3D =
        SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    const SpatialDomains::ExpansionInfoMap &expansions =
        graph3D->GetExpansionInfos();
    LibUtilities::BasisKey bkey =
        expansions.begin()->second->m_basisKeyVector[0];
    int nmodes = bkey.GetNumModes();
    cout << "Solving 3D C0 continuous Projection" << endl;
    cout << "    No. modes  : " << nmodes << endl;
    cout << endl;
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    Exp = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
        vSession, graph3D, vSession->GetVariable(0));
    //----------------------------------------------

    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetTotPoints();

    xc0 = Array<OneD, NekDouble>(nq, 0.0);
    xc1 = Array<OneD, NekDouble>(nq, 0.0);
    xc2 = Array<OneD, NekDouble>(nq, 0.0);

    switch (coordim)
    {
        case 1:
            Exp->GetCoords(xc0);
            break;
        case 2:
            Exp->GetCoords(xc0, xc1);
            break;
        case 3:
            Exp->GetCoords(xc0, xc1, xc2);
            break;
    }
    //----------------------------------------------

    //----------------------------------------------
    // Define forcing function
    fce = Array<OneD, NekDouble>(nq);

    for (i = 0; i < nq; ++i)
    {
        fce[i] = 0.0;
        for (j = 0; j < nmodes; ++j)
        {
            fce[i] += pow(xc0[i], j);
            fce[i] += pow(xc1[i], j);
            fce[i] += pow(xc2[i], j);
        }
    }

    //---------------------------------------------
    // Set up ExpList1D containing the solution
    Fce = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //---------------------------------------------

    //----------------------------------------------
    // Write solution
    // ofstream outfile("ProjectContFileOrig3D.dat");
    // Fce->WriteToFile(outfile,eGnuplot);
    // outfile.close();
    //----------------------------------------------

    //---------------------------------------------
    // Project onto Expansion
    Exp->FwdTrans(Fce->GetPhys(), Exp->UpdateCoeffs());
    //---------------------------------------------

    //-------------------------------------------
    // Backward Transform Solution to get projected values
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //-------------------------------------------

    //----------------------------------------------
    // Write solution
    // ofstream outfile2("ProjectContFile3D.dat");
    // Exp->WriteToFile(outfile2,eGnuplot);
    // outfile2.close();
    //----------------------------------------------

    //--------------------------------------------
    // Calculate L_inf error
    cout << "L infinity error: " << Exp->Linf(Fce->GetPhys()) << endl;
    cout << "L 2 error:        " << Exp->L2(Fce->GetPhys()) << endl;
    //--------------------------------------------

    vSession->Finalise();

    return 0;
}
