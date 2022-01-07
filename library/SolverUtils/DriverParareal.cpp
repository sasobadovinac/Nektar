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

void DriverParareal::RunCoarseSolve(const Array<OneD, const NekDouble> &input,
                                          Array<OneD,       NekDouble> &output)
{
    // Set to coarse timestep.
    NekDouble coarseDt = m_timestep * m_coarseSolveFactor;
    int nCoarseSteps = std::max(coarseDt / m_numChunks, 1);

    ASSERTL0(m_steps % nCoarseSteps == 0,
             "number of coarse steps should divide number of total steps");

    m_equ[0]->SetTimeStep(coarseDt);
    m_equ[0]->SetNumSteps(nCoarseSteps);

    // Copy initial condition from input.
    m_equ[0]->CopyToPhysField(0, input);

    m_equ[0]->DoInitialise();
    m_equ[0]->DoSolve();

    m_equ[0]->CopyFromPhysField(0, output);
}

void DriverParareal::RunFineSolve(const Array<OneD, const NekDouble> &input,
                                        Array<OneD,       NekDouble> &output)
{
    // Set to fine timestep.
    m_equ[0]->SetTimeStep(m_timestep);
    m_equ[0]->SetNumSteps(m_steps / m_numChunks);

    // Copy initial condition from input.
    m_equ[0]->CopyToPhysField(0, input);

    m_equ[0]->DoInitialise();
    m_equ[0]->DoSolve();

    m_equ[0]->CopyFromPhysField(0, output);
}

void DriverParareal::v_Execute(ostream &out)
{
    time_t starttime, endtime;
    NekDouble CPUtime;

    // Set parameters from original session file.
    m_timestep  = m_equ[0]->GetTimeStep();
    m_steps     = m_equ[0]->GetSteps();
    m_totalTime = m_timestep * m_steps;

    ASSERTL0(m_steps % m_numChunks == 0,
             "Total step size should be divisible by number of chunks.");

    m_equ[0]->PrintSummary(out);

    time(&starttime);

    // Allocate storage for coarse solver..
    int nPts = m_equ[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble>> solution(m_numChunks + 1);
    Array<OneD, Array<OneD, NekDouble>> solutionFine(m_numChunks + 1);
    Array<OneD, Array<OneD, NekDouble>> solutionCoarse1(m_numChunks + 1);
    Array<OneD, Array<OneD, NekDouble>> solutionCoarse2(m_numChunks + 1);

    for (int i = 0; i < solution.size(); ++i)
    {
        solution[i] = Array<OneD, NekDouble>(nPts);
        solutionFine[i] = Array<OneD, NekDouble>(nPts);
        solutionCoarse1[i] = Array<OneD, NekDouble>(nPts);
        solutionCoarse2[i] = Array<OneD, NekDouble>(nPts);
    }

    // Grab initial condition
    m_equ[0]->DoInitialise();
    m_equ[0]->CopyFromPhysField(0, solution[0]);

    // Run coarse solution, G(y_j^k) to get initial conditions
    for (int i = 0; i < m_numChunks; ++i)
    {
        RunCoarseSolve(solution[i], solution[i+1]);
    }

    for (int k = 0; k < m_numChunks; ++k)
    {
        if (k == 0)
        {
            for (int i = 0; i < m_numChunks + 1; ++i)
            {
                solutionCoarse1[i] = solution[i];
            }
        }
        else
        {
            // Run coarse solution, G(y_j^k)
            for (int i = 0; i < m_numChunks; ++i)
            {
                RunCoarseSolve(solution[i], solutionCoarse1[i+1]);
            }
        }

        // Run fine solution, F(y_j^k)
        solutionFine[0] = solution[0];
        for (int i = 0; i < m_numChunks; ++i)
        {
            RunFineSolve(solution[i], solutionFine[i+1]);
        }

        solutionCoarse2[0] = solution[0];
        for (int i = 0; i < m_numChunks; ++i)
        {
            RunCoarseSolve(solution[i], solutionCoarse2[i+1]);
            for (int q = 0; q < nPts; ++q)
            {
                solution[i+1][q] = solutionCoarse2[i+1][q] + solutionFine[i+1][q]
                    - solutionCoarse1[i+1][q];
            }
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

