///////////////////////////////////////////////////////////////////////////////
//
// File DriverParareal.cpp
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

#include <SolverUtils/DriverParareal.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
string DriverParareal::className = GetDriverFactory().RegisterCreatorFunction("Parareal", DriverParareal::create);
string DriverParareal::driverLookupId = LibUtilities::SessionReader::RegisterEnumValue("Driver","Parareal",0);

/**
 *
 */
DriverParareal::DriverParareal(
    const LibUtilities::SessionReaderSharedPtr pSession,
    const SpatialDomains::MeshGraphSharedPtr pGraph)
    : Driver(pSession, pGraph)
{
}

/**
 *
 */
DriverParareal:: ~DriverParareal()
{
}

/**
 *
 */
void DriverParareal::v_InitObject(ostream &out)
{
    Driver::v_InitObject(out);
}

void DriverParareal::RunCoarseSolve(const NekDouble                     time,
                                    const Array<OneD, const NekDouble> &input,
                                          Array<OneD,       NekDouble> &output)
{
    // Set to coarse timestep.
    NekDouble coarseDt = m_timestep * m_coarseSolveFactor;
    int nCoarseSteps = std::max((int)round(m_chunkTime / coarseDt), 1);

    std::cout << "RUNNING COARSE SOLVE: dt = " << coarseDt << " nsteps = " << nCoarseSteps << std::endl;

    ASSERTL0(m_steps % nCoarseSteps == 0,
             "number of coarse steps should divide number of total steps");

    m_equ[0]->SetTime(time);
    m_equ[0]->SetTimeStep(coarseDt);
    m_equ[0]->SetSteps(nCoarseSteps);

    // Copy initial condition from input.
    m_equ[0]->CopyToPhysField(0, input);

    //m_equ[0]->DoInitialise();
    m_equ[0]->DoSolve();

    m_equ[0]->CopyFromPhysField(0, output);
}

void DriverParareal::RunFineSolve(const NekDouble                     time,
                                  const Array<OneD, const NekDouble> &input,
                                        Array<OneD,       NekDouble> &output)
{
    std::cout << "RUNNING FINE SOLVE: dt = " << m_timestep << " nsteps = " << m_steps / m_numChunks << std::endl;

    // Set to fine timestep.
    m_equ[0]->SetTime(time);
    m_equ[0]->SetTimeStep(m_timestep);
    m_equ[0]->SetSteps(m_steps / m_numChunks);

    // Copy initial condition from input.
    m_equ[0]->CopyToPhysField(0, input);

    //m_equ[0]->DoInitialise();
    m_equ[0]->DoSolve();

    m_equ[0]->CopyFromPhysField(0, output);
}

void DriverParareal::v_Execute(ostream &out)
{
    time_t starttime, endtime;
    NekDouble CPUtime;

    m_numChunks = m_session->GetComm()->GetTimeComm()->GetSize();
    m_chunkRank = m_session->GetComm()->GetTimeComm()->GetRank();

    // Set parameters from original session file.
    m_timestep  = m_equ[0]->GetTimeStep();
    m_steps     = m_equ[0]->GetSteps();
    m_totalTime = m_timestep * m_steps;
    m_chunkTime = m_totalTime / m_numChunks;

    ASSERTL0(m_steps % m_numChunks == 0,
             "Total step size should be divisible by number of chunks.");

    m_equ[0]->PrintSummary(out);

    time(&starttime);

    // Allocate storage for coarse solver..
    int nPts = m_equ[0]->GetNpoints();
    Array<OneD, NekDouble> solution;
    Array<OneD, NekDouble> solutionFine;
    Array<OneD, NekDouble> solutionCoarse1;
    Array<OneD, NekDouble> solutionCoarse2;

    // Grab initial condition
    m_equ[0]->DoInitialise();

    if (m_chunkRank == 0) {
        m_equ[0]->CopyFromPhysField(0, solution);
    }

    // Evaluate "ExactSolution" function, or zero array
    Array<OneD, NekDouble> exactsoln(m_equ[0]->GetTotPoints(), 0.0);
    m_equ[0]->EvaluateExactSolution(0, exactsoln,
                                    m_equ[0]->GetFinalTime());

    // Run coarse solution, G(y_j^k) to get initial conditions
    std::cout << "** INITIAL CONDITION **" << std::endl;
    LibUtilities::CommSharedPtr tc = m_session->GetComm()->GetTimeComm();
    int recvProc = m_chunkRank - 1;
    int sendProc = m_chunkRank + 1;

    // Calculate the initial coarse solve approximation
    // This provides each time-slice with its initial condition.
    if (m_chunkRank > 0) {
        tc->Recv(recvProc, solution);
    }
    RunCoarseSolve(m_chunkRank, solution, solution);
    if (m_chunkRank < m_numChunks - 1)
    {
        tc->Send(sendProc, solution);
    }

    // On the last time-slice, we calculate the L2 error of our coarse solution
    if (m_chunkRank == m_numChunks - 1) {
        // Copy the calculated coarse solution back
        m_equ[0]->CopyToPhysField(0, solution);

        NekDouble vL2Error   = m_equ[0]->L2Error(0, exactsoln);

        std::cout << "COARSE SOLVE L2 error = " << vL2Error << std::endl;
    }

    // Iterate to improve on the approximation
    // For the moment, we use the maximum number of iterations
    // We can add a tolerance threshold later.
    // We start the iteration with the current approximation stored in
    // the 'solution' array.
    for (int k = 0; k < m_numChunks; ++k)
    {
        std::cout << "** ITERATION " << k << " **" << std::endl;

        // Calculate the coarse approximation G(y_j^k)
        // For the first iteration, this is already known, so could skip
        // Since this is based on previous solution, no comm necessary
        RunCoarseSolve(m_chunkRank * m_chunkTime, solution, solutionCoarse1);

        // Calculate fine solution, F(y_j^k)
        // Again no communication necessary
        RunFineSolve(m_chunkRank * m_chunkTime, solution, solutionFine);

        // Calculate coarse solve correction G(y_j^{k+1})
        // These are dependent on the previous time slice, so need to compute serially.
        if (m_chunkRank > 0) {
            // All time slices, apart from the first, receive their initial state from
            // the previous time slice.
            tc->Recv(recvProc, solution);
        }
        // Run the coarse solver
        RunCoarseSolve(m_chunkRank * m_chunkTime, solution, solutionCoarse2);
        // Calculate the new approximation y_{j+1}^{k+1}
        // This is calculated point-wise.
        for (int q = 0; q < nPts; ++q)
        {
            solution[q] = solutionCoarse2[q] + solutionFine[q] - solutionCoarse1[q];
        }
        // All but the last time slice should communicate the solution to the next time slice.
        // This will become the initial condition for the next slice.
        if (m_chunkRank < m_numChunks - 1) {
            tc->Send(sendProc, solution);
        }
        
        // On the last time-slice, we calculate the L2 error of our latest approximation
        if (m_chunkRank == m_numChunks - 1) {
            // Copy the calculated coarse solution back
            m_equ[0]->CopyToPhysField(0, solution);

            NekDouble vL2Error   = m_equ[0]->L2Error(0, exactsoln);

            std::cout << "COARSE SOLVE L2 error = " << vL2Error << std::endl;
        }
    }

    time(&endtime);

    m_equ[0]->Output();

    if (m_comm->GetRank() == 0)
    {
        CPUtime = difftime(endtime, starttime);
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

