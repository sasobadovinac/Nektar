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
                                    const Array<OneD, const Array<OneD, NekDouble>> &input,
                                          Array<OneD,       Array<OneD, NekDouble>> &output)
{
    // Compute coarse timestep.
    NekDouble coarseDt = m_timestep * m_coarseSolveFactor;
    int nCoarseSteps = std::max((int)round(m_chunkTime / coarseDt), 1);

    // Output coarse timestep.
    if (m_comm->GetRank() == 0)
    {
        std::cout << "RUNNING COARSE SOLVE: dt = " << coarseDt << " nsteps = " 
		  << nCoarseSteps << std::endl << std::flush;
    }

    ASSERTL0(m_steps % nCoarseSteps == 0,
             "number of coarse steps should divide number of total steps");

    // Set to coarse timestep.
    m_equ[m_nequ-1]->SetTime(time);
    m_equ[m_nequ-1]->SetTimeStep(coarseDt);
    m_equ[m_nequ-1]->SetSteps(nCoarseSteps);

    // Copy initial condition from input.
    for(int i = 0; i < m_equ[m_nequ-1]->GetNvariables(); ++i)
    {
        m_equ[m_nequ-1]->CopyToPhysField(i, input[i]);
    }

    // Solve equations.
    m_equ[m_nequ-1]->DoSolve();

    // Copy solution to output.
    for(int i = 0; i < m_equ[m_nequ-1]->GetNvariables(); ++i)
    {
        m_equ[m_nequ-1]->CopyFromPhysField(i, output[i]);
    }
}

void DriverParareal::RunFineSolve(const NekDouble                     time,
                                  const Array<OneD, const Array<OneD, NekDouble>> &input,
                                        Array<OneD,       Array<OneD, NekDouble>> &output)
{

    // Output fine timestep.
    if (m_comm->GetRank() == 0)
    {
        std::cout << "RUNNING FINE SOLVE: dt = " << m_timestep << " nsteps = " 
		  << m_steps / m_numChunks << std::endl << std::flush;
    }

    // Set to fine timestep.
    m_equ[0]->SetTime(time);
    m_equ[0]->SetTimeStep(m_timestep);
    m_equ[0]->SetSteps(m_steps / m_numChunks);

    // Copy initial condition from input.
    for(int i = 0; i < m_equ[0]->GetNvariables(); ++i)
    {
        m_equ[0]->CopyToPhysField(i, input[i]);
    }

    // Solve equations.
    m_equ[0]->DoSolve();

    // Copy solution to output.
    for(int i = 0; i < m_equ[0]->GetNvariables(); ++i)
    {
        m_equ[0]->CopyFromPhysField(i, output[i]);
    }
}

void DriverParareal::v_Execute(ostream &out)
{
    time_t starttime, endtime;
    NekDouble CPUtime;

    m_numChunks = m_session->GetComm()->GetTimeComm()->GetSize();
    m_chunkRank = m_session->GetComm()->GetTimeComm()->GetRank();
    // Maximum number of parareal iteration
    if (m_session->DefinesParameter("PararealIterMax"))
    {
        m_pararealIterMax = m_session->GetParameter("PararealIterMax");
    }
    else
    {
        m_pararealIterMax = m_numChunks;
    }
    // Coarse solver time factor
    if (m_session->DefinesParameter("CoarseSolveFactor"))
    {
        m_coarseSolveFactor = m_session->GetParameter("CoarseSolveFactor");
    }

    // Set parameters from session file.
    m_timestep  = m_equ[0]->GetTimeStep();
    m_steps     = m_equ[0]->GetSteps();
    m_totalTime = m_timestep * m_steps;
    m_chunkTime = m_totalTime / m_numChunks;

    ASSERTL0(m_steps % m_numChunks == 0,
             "Total step size should be divisible by number of chunks.");

    // Fine solver summary
    if (m_comm->GetRank() == 0)
    {
        std::cout << "=======================================================================" << std::endl << std::flush;
        std::cout << "========================= FINE PROPAGATOR INFO ========================" << std::endl << std::flush;
    }
    m_equ[0]->PrintSummary(out);

    std::cout << std::endl << std::flush;

    // Coarse solver summary
    if (m_nequ-1)
    {
        if (m_comm->GetRank() == 0)
        {
            std::cout << "=======================================================================" << std::endl << std::flush;
            std::cout << "======================== COARSE PROPAGATOR INFO =======================" << std::endl << std::flush;
        }
        m_equ[m_nequ-1]->PrintSummary(out);
    }

    time(&starttime);

    // Allocate storage for coarse solver..
    int nPts = m_equ[0]->GetNpoints();
    int nVar = m_equ[0]->GetNvariables();
    Array<OneD, Array<OneD, NekDouble>> solution(nVar);
    Array<OneD, Array<OneD, NekDouble>> solutionFine(nVar);
    Array<OneD, Array<OneD, NekDouble>> solutionCoarse1(nVar);
    Array<OneD, Array<OneD, NekDouble>> solutionCoarse2(nVar);
    Array<OneD, Array<OneD, NekDouble>> ic(nVar);
    Array<OneD, Array<OneD, NekDouble>> exactsoln(nVar);
    for(int i = 0; i < nVar; ++i)
    {
	    solution[i]        = Array<OneD, NekDouble>(nPts, 0.0);
        solutionFine[i]    = Array<OneD, NekDouble>(nPts, 0.0);
	    solutionCoarse1[i] = Array<OneD, NekDouble>(nPts, 0.0);
	    solutionCoarse2[i] = Array<OneD, NekDouble>(nPts, 0.0);
	    ic[i]              = Array<OneD, NekDouble>(nPts, 0.0);
	    exactsoln[i]       = Array<OneD, NekDouble>(nPts, 0.0);
    }

    // Initialize fine solver
    m_equ[0]->DoInitialise();

    // Initialize coarse solver
    if (m_nequ-1)
    {
        m_equ[m_nequ-1]->DoInitialise();
    }

    // Get initial conditions
    if (m_chunkRank == 0) 
    {
        for(int i = 0; i < nVar; ++i)
        {
            m_equ[0]->CopyFromPhysField(i, ic[i]);
        }
    }

    // Evaluate "ExactSolution" function, or zero array
    if (m_chunkRank == m_numChunks - 1)
    { 
        for(int i = 0; i < nVar; ++i)
        {
            m_equ[0]->EvaluateExactSolution(i, exactsoln[i], m_totalTime);
        }
    }

    // Run coarse solution, G(y_j^k) to get initial conditions
    if (m_chunkRank == 0 && m_comm->GetRank() == 0)
    {
        std::cout << "** INITIAL CONDITION **" << std::endl << std::flush;
    }
    LibUtilities::CommSharedPtr tc = m_session->GetComm()->GetTimeComm();
    int recvProc = m_chunkRank - 1;
    int sendProc = m_chunkRank + 1;

    // Calculate the initial coarse solve approximation
    // This provides each time-slice with its initial condition.
    if (m_chunkRank > 0)
    {
        for(int i = 0; i < nVar; ++i)
        {
            tc->Recv(recvProc, ic[i]);
        }
    }
    RunCoarseSolve(m_chunkRank * m_chunkTime, ic, solution);
    if (m_chunkRank < m_numChunks - 1)
    {
        for(int i = 0; i < nVar; ++i)
        {
            tc->Send(sendProc, solution[i]);
        }
    }

    // On the last time-slice, we calculate the L2 error of our coarse solution
    if (m_chunkRank == m_numChunks - 1) 
    {
        for(int i = 0; i < nVar; ++i)
        {
            // Copy the calculated coarse solution back
            m_equ[0]->CopyToPhysField(i, solution[i]);

            NekDouble vL2Error   = m_equ[0]->L2Error(i, exactsoln[i]);
            NekDouble vLinfError = m_equ[0]->LinfError(i, exactsoln[i]);

            std::cout << "L 2 error (variable " << m_equ[0]->GetVariable(i)
                      << ") : " << vL2Error << std::endl << std::flush;
            std::cout << "L inf error (variable " << m_equ[0]->GetVariable(i)
                      << ") : " << vLinfError << std::endl << std::flush;
    	}
    }

    // Iterate to improve on the approximation
    // For the moment, we use the maximum number of iterations
    // We can add a tolerance threshold later.
    // We start the iteration with the current approximation stored in
    // the 'solution' array.
    for (int k = 0; k < m_pararealIterMax; ++k)
    {
        if (m_chunkRank == 0 && m_comm->GetRank() == 0)
        {
            std::cout << "** ITERATION " << k << " **" << std::endl << std::flush;
        }

        // Calculate the coarse approximation G(y_j^k)
        // For the first iteration, this is already known, so could skip
        // Since this is based on previous solution, no comm necessary
        RunCoarseSolve(m_chunkRank * m_chunkTime, ic, solutionCoarse1);

        // Calculate fine solution, F(y_j^k)
        // Again no communication necessary
        RunFineSolve(m_chunkRank * m_chunkTime, ic, solutionFine);

        // Calculate coarse solve correction G(y_j^{k+1})
        // These are dependent on the previous time slice, so need to compute serially.
        if (m_chunkRank > 0)
        {
            // All time slices, apart from the first, receive their initial state from
            // the previous time slice.
            for(int i = 0; i < nVar; ++i)
            {
                tc->Recv(recvProc, ic[i]);
            }
        }

        // Run the coarse solver
        RunCoarseSolve(m_chunkRank * m_chunkTime, ic, solutionCoarse2);

        // Calculate the new approximation y_{j+1}^{k+1}
        // This is calculated point-wise.
        for(int i = 0; i < nVar; ++i)
        {
            for (int q = 0; q < nPts; ++q)
            {
                solution[i][q] = solutionCoarse2[i][q] + solutionFine[i][q] - solutionCoarse1[i][q];
            }
        }

        // All but the last time slice should communicate the solution to the next time slice.
        // This will become the initial condition for the next slice.
        if (m_chunkRank < m_numChunks - 1)
        {
            for(int i = 0; i < nVar; ++i)
            {
                tc->Send(sendProc, solution[i]);
            }
        }

        // On the last time-slice, we calculate the L2 error of our latest approximation
        if (m_chunkRank == m_numChunks - 1)
        {
            for(int i = 0; i < nVar; ++i)
            {
                // Copy the new approximation back
                m_equ[0]->CopyToPhysField(i, solution[i]);

    	        NekDouble vL2Error   = m_equ[0]->L2Error(i, exactsoln[i]);
                NekDouble vLinfError = m_equ[0]->LinfError(i, exactsoln[i]);

	            std::cout << "L 2 error (variable " << m_equ[0]->GetVariable(i)
                          << ") : " << vL2Error << std::endl << std::flush;
                std::cout << "L inf error (variable " << m_equ[0]->GetVariable(i)
                          << ") : " << vLinfError << std::endl << std::flush;
            }
        }
    }

    time(&endtime);

    m_equ[0]->Output();

    if (m_comm->GetRank() == 0)
    {
        CPUtime = difftime(endtime, starttime);
        std::cout << "-------------------------------------------" << std::endl << std::flush;
        std::cout << "Total Computation Time = " << CPUtime << "s" << std::endl << std::flush;
        std::cout << "-------------------------------------------" << std::endl << std::flush;
    }

}
}
}

