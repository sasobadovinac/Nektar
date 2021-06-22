///////////////////////////////////////////////////////////////////////////////
//
// File DriverALE.cpp
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
// Description: Incompressible Navier Stokes solver
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <SolverUtils/DriverALE.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

string DriverALE::className = GetDriverFactory().RegisterCreatorFunction("ALE", DriverALE::create);
string DriverALE::driverLookupId = LibUtilities::SessionReader::RegisterEnumValue("Driver","ALE",0);

/**
 *
 */
DriverALE::DriverALE(
    const LibUtilities::SessionReaderSharedPtr pSession,
    const SpatialDomains::MeshGraphSharedPtr pGraph)
    : Driver(pSession, pGraph)
{
}


/**
 *
 */
DriverALE:: ~DriverALE()
{
}


/**
 *
 */
void DriverALE::v_InitObject(ostream &out)
{
    Driver::v_InitObject(out);
}


void DriverALE::v_Execute(ostream &out)
{
    clock_t starttime, endtime;
    NekDouble CPUtime;

    // Store all points
    std::map<int, SpatialDomains::PointGeom> pts;
    auto &ptsMap = m_graph->GetAllPointGeoms();
    for (auto &pt : ptsMap)
    {
        pts[pt.first] = *(pt.second);
    }

    m_equ[0]->PrintSummary(out);

    starttime = clock();

    // Start with a single timestep.
    m_equ[0]->SetPts(pts);
    m_equ[0]->DoInitialise();
    m_equ[0]->DoSolve();

    int nRuns, nVariables = m_equ[0]->GetNvariables();
    m_session->LoadParameter("NumRuns", nRuns, 1);

    NekDouble dt = m_session->GetParameter("TimeStep"), time = dt;

    for (int i = 0; i < nRuns; ++i)
    {
        Array<OneD, Array<OneD, NekDouble> > fielddata(nVariables);
        for (int n = 0; n < nVariables; ++n)
        {
            fielddata[n] = m_equ[0]->UpdateFields()[n]->UpdateCoeffs();
        }

        if (LibUtilities::NekManager<MultiRegions::GlobalLinSysKey,
            MultiRegions::GlobalLinSys>::
            PoolCreated(std::string("GlobalLinSys")))
        {
            std::cout << "clearing globallinsys pool" << std::endl;
            LibUtilities::NekManager<MultiRegions::GlobalLinSysKey,
                                     MultiRegions::GlobalLinSys>::
                ClearManager(std::string("GlobalLinSys"));
        }

        int chkNumber = m_equ[0]->GetCheckpointNumber();
        int chkSteps  = m_equ[0]->GetCheckpointSteps();

        // Initialise driver again
        Driver::v_InitObject(out);

        // Set chkSteps to zero to avoid writing initial condition
        m_equ[0]->SetCheckpointSteps(0);

        // Initialise equation
        m_equ[0]->SetPts(pts);
        m_equ[0]->SetTime(time);
        m_equ[0]->DoInitialise();
        m_equ[0]->SetInitialStep(i);
        m_equ[0]->SetSteps(1);
        m_equ[0]->SetTime(time);
        m_equ[0]->SetBoundaryConditions(time);
        m_equ[0]->SetCheckpointNumber(chkNumber);
        m_equ[0]->SetCheckpointSteps(chkSteps);

        // Project solution to new expansion
        for (int n = 0; n < nVariables; n++)
        {
            m_equ[0]->UpdateFields()[n]->UpdateCoeffs() = fielddata[n];
            m_equ[0]->UpdateFields()[n]->BwdTrans_IterPerExp(
                m_equ[0]->UpdateFields()[n]->GetCoeffs(),
                m_equ[0]->UpdateFields()[n]->UpdatePhys());
        }

        // Solve equation
        m_equ[0]->DoSolve();

        time += dt;
    }

    endtime = clock();

    m_equ[0]->Output();

    if (m_comm->GetRank() == 0)
    {
        CPUtime = (endtime - starttime) / NekDouble(CLOCKS_PER_SEC);
        cout << "-------------------------------------------" << endl;
        cout << "Total Computation Time = " << CPUtime << "s" << endl;
        cout << "-------------------------------------------" << endl;
    }

    // Evaluate and output computation time and solution accuracy.
    // The specific format of the error output is essential for the
    // regression tests to work.
    // Evaluate L2 Error
    for(int i = 0; i < m_equ[0]->GetNvariables(); ++i)
    {
        Array<OneD, NekDouble> exactsoln(m_equ[0]->GetTotPoints(), 0.0);

        // Evaluate "ExactSolution" function, or zero array
        m_equ[0]->EvaluateExactSolution(i, exactsoln,
                                        m_equ[0]->GetFinalTime());

        NekDouble vL2Error   = m_equ[0]->L2Error  (i, exactsoln);
        NekDouble vLinfError = m_equ[0]->LinfError(i, exactsoln);

        if (m_comm->GetRank() == 0)
        {
            out << "L 2 error (variable " << m_equ[0]->GetVariable(i)
                << ") : " << vL2Error << endl;
            out << "L inf error (variable " << m_equ[0]->GetVariable(i)
                << ") : " << vLinfError << endl;
        }
    }
}
}
}

