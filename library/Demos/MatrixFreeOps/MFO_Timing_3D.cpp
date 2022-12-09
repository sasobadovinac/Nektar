///////////////////////////////////////////////////////////////////////////////
//
// File: MFO_Timing_3D.cpp
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
// Description: Small demo to run timings of various operators.
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <iomanip>

#include <Collections/Collection.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ContField.h>
#include <MultiRegions/ExpList.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

MultiRegions::ExpListSharedPtr SetupExpList(
    int N, LibUtilities::SessionReaderSharedPtr session,
    SpatialDomains::MeshGraphSharedPtr graph,
    Collections::ImplementationType impType)
{
    graph->SetExpansionInfoToNumModes(N);

    MultiRegions::ExpListSharedPtr expList =
        MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(session, graph);

    expList->CreateCollections(impType);

    return expList;
}

void printOutput(string testname, NekDouble total_sec)
{

    cout << setw(18) << testname << setw(18) << total_sec << endl;
}

int main(int argc, char *argv[])
{
    LibUtilities::Timer timer;

    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(argc, argv);

    // MultiRegions::ContFieldSharedPtr Exp,Fce;
    int nq, coordim;
    Array<OneD, NekDouble> fce;
    Array<OneD, NekDouble> xc0, xc1, xc2;
    StdRegions::ConstFactorMap factors;
    StdRegions::VarCoeffMap varcoeffs;

    try
    {
        LibUtilities::FieldIOSharedPtr fld =
            LibUtilities::FieldIO::CreateDefault(vSession);

        //----------------------------------------------
        // Read in mesh from input file
        SpatialDomains::MeshGraphSharedPtr graph2D =
            SpatialDomains::MeshGraph::Read(vSession);
        //----------------------------------------------

        //----------------------------------------------
        // Define number of times to loop each operation
        int iter;

        if (vSession->DefinesParameter("TimingIterations"))
        {
            iter = vSession->GetParameter("TimingIterations");
        }
        else
        {
            iter = 25;
        }

        //----------------------------------------------

        //----------------------------------------------
        // Print summary of solution details
        factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
        const SpatialDomains::ExpansionInfoMap &expansions =
            graph2D->GetExpansionInfo();
        LibUtilities::BasisKey bkey0 =
            expansions.begin()->second->m_basisKeyVector[0];

        if (vSession->GetComm()->GetRank() == 0)
        {
            cout << "======================================================"
                 << endl;
            cout << "  Timing Helmholtz Solve: " << endl;
            cout << "           Communication: "
                 << vSession->GetComm()->GetType() << endl;
            cout << "           Solver type  : "
                 << vSession->GetSolverInfo("GlobalSysSoln") << endl;
            cout << "           No. modes    : " << bkey0.GetNumModes() << endl;
            cout << "           No. of iters : " << iter << endl;
            cout << "======================================================"
                 << endl;
            cout << endl;
        }
        //----------------------------------------------

        //----------------------------------------------
        // Define Expansion
        int impList[3] = {5, 5, 5}; // 2=IterPerExp, 5=matrixfreeops

        for (int expIter = 0; expIter < 3; expIter++)
        {

            //----------------------------------------------
            // Iterate and run each expansion

            int imp = impList[expIter];

            Collections::ImplementationType impType =
                (Collections::ImplementationType)imp;

            MultiRegions::ContFieldSharedPtr Exp =
                MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
                    vSession, graph2D, vSession->GetVariable(0));
            Exp->CreateCollections(impType);

            std::cout << "Using " << Collections::ImplementationTypeMap[imp]
                      << " Collection Implementation:" << std::endl;

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
                case 2:
                    Exp->GetCoords(xc0, xc1);
                    break;
                case 3:
                    Exp->GetCoords(xc0, xc1, xc2);
                    break;
                default:
                    ASSERTL0(false, "Coordim not valid");
                    break;
            }
            //----------------------------------------------

            //----------------------------------------------
            // Set up variable coefficients if defined
            if (vSession->DefinesFunction("d00"))
            {
                Array<OneD, NekDouble> d00(nq, 0.0);
                LibUtilities::EquationSharedPtr d00func =
                    vSession->GetFunction("d00", 0);
                d00func->Evaluate(xc0, xc1, xc2, d00);
                varcoeffs[StdRegions::eVarCoeffD00] = d00;
            }

            if (vSession->DefinesFunction("d01"))
            {
                Array<OneD, NekDouble> d01(nq, 0.0);
                LibUtilities::EquationSharedPtr d01func =
                    vSession->GetFunction("d01", 0);
                d01func->Evaluate(xc0, xc1, xc2, d01);
                varcoeffs[StdRegions::eVarCoeffD01] = d01;
            }

            if (vSession->DefinesFunction("d11"))
            {
                Array<OneD, NekDouble> d11(nq, 0.0);
                LibUtilities::EquationSharedPtr d11func =
                    vSession->GetFunction("d11", 0);
                d11func->Evaluate(xc0, xc1, xc2, d11);
                varcoeffs[StdRegions::eVarCoeffD11] = d11;
            }

            if (vSession->DefinesFunction("d02"))
            {
                Array<OneD, NekDouble> d02(nq, 0.0);
                LibUtilities::EquationSharedPtr d02func =
                    vSession->GetFunction("d02", 0);
                d02func->Evaluate(xc0, xc1, xc2, d02);
                varcoeffs[StdRegions::eVarCoeffD02] = d02;
            }

            if (vSession->DefinesFunction("d12"))
            {
                Array<OneD, NekDouble> d12(nq, 0.0);
                LibUtilities::EquationSharedPtr d12func =
                    vSession->GetFunction("d12", 0);
                d12func->Evaluate(xc0, xc1, xc2, d12);
                varcoeffs[StdRegions::eVarCoeffD12] = d12;
            }

            if (vSession->DefinesFunction("d22"))
            {
                Array<OneD, NekDouble> d22(nq, 0.0);
                LibUtilities::EquationSharedPtr d22func =
                    vSession->GetFunction("d22", 0);
                d22func->Evaluate(xc0, xc1, xc2, d22);
                varcoeffs[StdRegions::eVarCoeffD22] = d22;
            }
            //----------------------------------------------

            //----------------------------------------------
            // Set up const diffusion coefficients if defined
            if (vSession->DefinesParameter("d00"))
            {
                NekDouble d00;
                vSession->LoadParameter("d00", d00, 1.0);
                factors[StdRegions::eFactorCoeffD00] = d00;
            }

            if (vSession->DefinesParameter("d01"))
            {
                NekDouble d01;
                vSession->LoadParameter("d01", d01, 1.0);
                factors[StdRegions::eFactorCoeffD01] = d01;
            }

            if (vSession->DefinesParameter("d11"))
            {
                NekDouble d11;
                vSession->LoadParameter("d11", d11, 1.0);
                factors[StdRegions::eFactorCoeffD11] = d11;
            }

            if (vSession->DefinesParameter("d02"))
            {
                NekDouble d02;
                vSession->LoadParameter("d02", d02, 1.0);
                factors[StdRegions::eFactorCoeffD02] = d02;
            }

            if (vSession->DefinesParameter("d12"))
            {
                NekDouble d12;
                vSession->LoadParameter("d12", d12, 1.0);
                factors[StdRegions::eFactorCoeffD12] = d12;
            }

            if (vSession->DefinesParameter("d22"))
            {
                NekDouble d22;
                vSession->LoadParameter("d22", d22, 1.0);
                factors[StdRegions::eFactorCoeffD22] = d22;
            }

            if (vSession->DefinesParameter("fn_vardiff"))
            {
                NekDouble tau;
                vSession->LoadParameter("fn_vardiff", tau, 1.0);
                if (tau > 0)
                {
                    factors[StdRegions::eFactorTau] = 1.0;
                }
            }
            //----------------------------------------------

            //----------------------------------------------
            // Define forcing function for first variable defined in file
            fce = Array<OneD, NekDouble>(nq);
            LibUtilities::EquationSharedPtr ffunc =
                vSession->GetFunction("Forcing", 0);
            ffunc->Evaluate(xc0, xc1, xc2, fce);

            //----------------------------------------------

            //----------------------------------------------
            // Setup expansion containing the  forcing function
            MultiRegions::ContFieldSharedPtr Fce =
                MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(*Exp);
            Fce->SetPhys(fce);
            //----------------------------------------------

            //----------------------------------------------
            // Time backwards transform
            NekDouble total_sec;

            Array<OneD, NekDouble> bwdTout(Exp->GetNpoints());

            timer.Start();
            for (int j = 0; j < iter; j++)
            {
                Exp->BwdTrans(Exp->GetCoeffs(), bwdTout);
            }
            timer.Stop();

            total_sec = timer.TimePerTest(iter);
            printOutput("BwdTrans", total_sec);
            //----------------------------------------------

            //----------------------------------------------
            // Time IProd WRT Base

            Array<OneD, NekDouble> ipwrtbout(Exp->GetNpoints());

            timer.Start();
            for (int j = 0; j < iter; j++)
            {
                Exp->IProductWRTBase(bwdTout, ipwrtbout);
            }
            timer.Stop();

            total_sec = timer.TimePerTest(iter);
            printOutput("IProdWRTBase", total_sec);
            //----------------------------------------------

            //----------------------------------------------
            // Time PhysDeriv
            Array<OneD, NekDouble> pdout0(Exp->GetNpoints());
            Array<OneD, NekDouble> pdout1(Exp->GetNpoints());
            Array<OneD, NekDouble> pdout2(Exp->GetNpoints());

            timer.Start();
            for (int j = 0; j < iter; j++)
            {
                Exp->PhysDeriv(bwdTout, pdout0, pdout1, pdout2);
            }
            timer.Stop();

            total_sec = timer.TimePerTest(iter);
            printOutput("PhysDeriv", total_sec);
            //----------------------------------------------

            //----------------------------------------------
            // Time IProd WRT DerivBase
            Array<OneD, Array<OneD, NekDouble>> inputD(Exp->GetCoordim(0));

            inputD[0] = pdout0;
            inputD[1] = pdout1;
            inputD[2] = pdout2;

            timer.Start();
            for (int j = 0; j < iter; j++)
            {
                Exp->IProductWRTDerivBase(inputD, ipwrtbout);
            }
            timer.Stop();

            total_sec = timer.TimePerTest(iter);
            printOutput("IProdWRTDerivB", total_sec);
            //----------------------------------------------

            //----------------------------------------------
            // Helmholtz solution taking physical forcing after setting
            // initial condition to zero

            Vmath::Zero(Exp->GetNcoeffs(), Exp->UpdateCoeffs(), 1);
            timer.Start(); // time Helmholtz
            for (int j = 0; j < iter; j++)
            {
                Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), factors,
                               varcoeffs);
            }
            timer.Stop();

            total_sec = timer.TimePerTest(iter);
            printOutput("HelmSolve", total_sec);

            // redo helmholtz for errors
            Vmath::Zero(Exp->GetNcoeffs(), Exp->UpdateCoeffs(), 1);
            Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), factors,
                           varcoeffs);
            //----------------------------------------------

            //----------------------------------------------
            // Backward Transform Solution to get solved values
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

                Fce->SetPhys(fce);
                Fce->SetPhysState(true);
                //--------------------------------------------

                //--------------------------------------------
                // Calculate errors
                NekDouble vLinfError =
                    Exp->Linf(Exp->GetPhys(), Fce->GetPhys());
                NekDouble vL2Error = Exp->L2(Exp->GetPhys(), Fce->GetPhys());
                NekDouble vH1Error = Exp->H1(Exp->GetPhys(), Fce->GetPhys());
                if (vSession->GetComm()->GetRank() == 0)
                {
                    cout << endl << "Errors : " << endl;
                    cout << "L infinity error: " << vLinfError << endl;
                    cout << "L 2 error:        " << vL2Error << endl;
                    cout << "H 1 error:        " << vH1Error << endl;
                    cout << "=================================================="
                            "===="
                         << endl;
                    cout << endl;
                }
                //--------------------------------------------
            }
            //----------------------------------------------
        }
    }
    catch (const std::runtime_error &)
    {
        cout << "Caught an error" << endl;
        return 1;
    }

    vSession->Finalise();

    return 0;
}
