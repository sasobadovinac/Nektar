///////////////////////////////////////////////////////////////////////////////
//
// File: Helmholtz3DHomo2D.cpp
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
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();

    MultiRegions::ContField3DHomogeneous2DSharedPtr Exp, Fce;
    int nq;
    Array<OneD, NekDouble> fce;
    Array<OneD, NekDouble> xc0, xc1, xc2;
    StdRegions::ConstFactorMap factors;
    FlagList flags;

    if ((argc != 2) && (argc != 3))
    {
        fprintf(stderr, "Usage: Helmholtz3DHomo2D meshfile [SysSolnType]   \n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph1D =
        SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int nypoints;
    int nzpoints;
    NekDouble ly;
    NekDouble lz;
    int FFT;

    vSession->LoadParameter("HomModesY", nypoints);
    vSession->LoadParameter("HomModesZ", nzpoints);
    vSession->LoadParameter("LY", ly);
    vSession->LoadParameter("LZ", lz);
    vSession->LoadParameter("USEFFT", FFT);

    bool useFFT = false;
    bool deal   = false;
    if (FFT == 1)
    {
        useFFT = true;
    }

    const LibUtilities::PointsKey PkeyY(nypoints,
                                        LibUtilities::eFourierEvenlySpaced);
    const LibUtilities::BasisKey BkeyY(LibUtilities::eFourier, nypoints, PkeyY);

    const LibUtilities::PointsKey PkeyZ(nzpoints,
                                        LibUtilities::eFourierEvenlySpaced);
    const LibUtilities::BasisKey BkeyZ(LibUtilities::eFourier, nzpoints, PkeyZ);

    Exp = MemoryManager<MultiRegions::ContField3DHomogeneous2D>::
        AllocateSharedPtr(vSession, BkeyY, BkeyZ, ly, lz, useFFT, deal, graph1D,
                          vSession->GetVariable(0));
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");

    const SpatialDomains::ExpansionInfoMap &expansions =
        graph1D->GetExpansionInfo();

    LibUtilities::BasisKey bkey0 =
        expansions.begin()->second->m_basisKeyVector[0];

    cout << "Solving 3D Helmholtz (Homogeneous in yz-plane):" << endl;
    cout << "         Lambda          : " << factors[StdRegions::eFactorLambda]
         << endl;
    cout << "         Ly              : " << ly << endl;
    cout << "         Lz              : " << lz << endl;
    cout << "         N.modes         : " << bkey0.GetNumModes() << endl;
    cout << "         N.Y homo modes  : " << BkeyY.GetNumModes() << endl;
    cout << "         N.Z homo modes  : " << BkeyZ.GetNumModes() << endl;
    cout << endl;
    //----------------------------------------------

    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    nq  = Exp->GetTotPoints();
    xc0 = Array<OneD, NekDouble>(nq, 0.0);
    xc1 = Array<OneD, NekDouble>(nq, 0.0);
    xc2 = Array<OneD, NekDouble>(nq, 0.0);

    Exp->GetCoords(xc0, xc1, xc2);
    //----------------------------------------------

    //----------------------------------------------
    // Define forcing function for first variable defined in file
    fce = Array<OneD, NekDouble>(nq);

    LibUtilities::EquationSharedPtr ffunc = vSession->GetFunction("Forcing", 0);

    ffunc->Evaluate(xc0, xc1, xc2, fce);

    //----------------------------------------------

    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<
        MultiRegions::ContField3DHomogeneous2D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------

    //----------------------------------------------
    // Helmholtz solution taking physical forcing
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), factors);
    //----------------------------------------------

    //----------------------------------------------
    // Backward Transform Solution to get solved values at
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //----------------------------------------------

    //----------------------------------------------
    // See if there is an exact solution, if so
    // evaluate and plot errors
    LibUtilities::EquationSharedPtr ex_sol =
        vSession->GetFunction("ExactSolution", 0);

    if (ex_sol)
    {
        //----------------------------------------------
        // evaluate exact solution

        ex_sol->Evaluate(xc0, xc1, xc2, fce);

        //----------------------------------------------

        //--------------------------------------------
        // Calculate L_inf error
        Fce->SetPhys(fce);
        Fce->SetPhysState(true);

        cout << "L infinity error: "
             << Exp->Linf(Exp->GetPhys(), Fce->GetPhys()) << endl;
        cout << "L 2 error:        " << Exp->L2(Exp->GetPhys(), Fce->GetPhys())
             << endl;
        //--------------------------------------------
    }
    //----------------------------------------------

    return 0;
}
