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
// Description: Driver class for the parareal solver
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
    try
    {
        // Retrieve the equation system to solve.
        ASSERTL0(m_session->DefinesSolverInfo("EqType"),
                 "EqType SolverInfo tag must be defined.");
        std::string vEquation = m_session->GetSolverInfo("EqType");
        if (m_session->DefinesSolverInfo("SolverType"))
        {
            vEquation = m_session->GetSolverInfo("SolverType");
        }

        // Check such a module exists for this equation.
        ASSERTL0(GetEquationSystemFactory().ModuleExists(vEquation),
                 "EquationSystem '" + vEquation + "' is not defined.\n"
                 "Ensure equation name is correct and module is compiled.\n");

        // Retrieve the type of evolution operator to use
        /// @todo At the moment this is Navier-Stokes specific - generalise?
        m_EvolutionOperator =
            m_session->GetSolverInfoAsEnum<EvolutionOperatorType>(
                    "EvolutionOperator");

        m_nequ = 2; 

        m_equ = Array<OneD, EquationSystemSharedPtr>(m_nequ);

        // Set the AdvectiveType tag and create EquationSystem objects.
        switch (m_EvolutionOperator)
        {
            case eNonlinear:
                // Set coarse parareal session file
                SetPararealSessionFile();

                // Set fine parareal solver
                m_session->SetTag("AdvectiveType","Convective");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);

                // Set coarse parareal solver
                m_sessionCoarse->SetTag("AdvectiveType", "Convective");
                m_equ[1] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_sessionCoarse, m_graphCoarse);
                break;
            case eDirect:
                // Set coarse parareal session file
                SetPararealSessionFile();

                // Set fine parareal solver
                m_session->SetTag("AdvectiveType","Linearised");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);

                // Set coarse parareal solver
                m_sessionCoarse->SetTag("AdvectiveType", "Linearised");
                m_equ[1] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_sessionCoarse, m_graphCoarse);
                break;
            case eAdjoint:
                // Set coarse parareal session file
                SetPararealSessionFile();

                // Set fine parareal solver
                m_session->SetTag("AdvectiveType","Adjoint");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);

                // Set coarse parareal solver
                m_sessionCoarse->SetTag("AdvectiveType", "Adjoint");
                m_equ[1] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_sessionCoarse, m_graphCoarse);
                break;
            case eSkewSymmetric:
                // Set coarse parareal session file
                SetPararealSessionFile();

                // Set fine parareal solver
                m_session->SetTag("AdvectiveType","SkewSymmetric");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);

                // Set coarse parareal solver
                m_sessionCoarse->SetTag("AdvectiveType", "SkewSymmetric");
                m_equ[1] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_sessionCoarse, m_graphCoarse);
                break;
            default:
                ASSERTL0(false, "Unrecognised evolution operator.");

        }
    }
    catch (int e)
    {
        ASSERTL0(e == -1, "No such class class defined.");
        out << "An error occurred during driver initialisation." << endl;
    }
}

/// Set the Parareal (coarse solver) session file
void DriverParareal::SetPararealSessionFile(void)
{
    // Get the coarse solver session file
    string meshFile;
    string coarseSolverFile;
    vector<string> coarseSolverFilenames;
    meshFile         = m_session->GetFilenames()[0];
    coarseSolverFile = m_session->GetFilenames().size()>1 ? 
          m_session->GetFilenames()[1] : 
          m_session->GetFilenames()[0];
    coarseSolverFile = coarseSolverFile.substr(0, coarseSolverFile.size() - 4);
    coarseSolverFile += "_coarseSolver.xml";
    std::ifstream f(coarseSolverFile);
    if (f.good())
    {
        // if _coarseSolver.xml exit, read session file
        if (m_session->GetFilenames().size()>1)
        {    
            coarseSolverFilenames.push_back(meshFile);
        }
        coarseSolverFilenames.push_back(coarseSolverFile);
    }
    else
    {
        // if _coarseSolver.xml does not exit, use original session file
        coarseSolverFilenames.push_back(m_session->GetFilenames()[0]);
        if (m_session->GetFilenames().size()>1)
        {    
            coarseSolverFilenames.push_back(m_session->GetFilenames()[1]);
        }
    }

    // Define argument for the coarse parareal solver
    int npx = m_session->DefinesCmdLineArgument("npx")
                  ? m_session->GetCmdLineArgument<int>("npx")
                  : 1;
    int npy = m_session->DefinesCmdLineArgument("npy")
                  ? m_session->GetCmdLineArgument<int>("npy")
                  : 1;
    int npz = m_session->DefinesCmdLineArgument("npz")
                  ? m_session->GetCmdLineArgument<int>("npz")
                  : 1;
    int nsz = m_session->DefinesCmdLineArgument("nsz")
                  ? m_session->GetCmdLineArgument<int>("nsz")
                  : 1;
    int npt = m_session->DefinesCmdLineArgument("npt")
                  ? m_session->GetCmdLineArgument<int>("npt")
                  : 1;

    char *argv[] = {const_cast<char *>("Solver"), // this is just a place holder
                    const_cast<char *>("--npx"),
                    const_cast<char *>(std::to_string(npx).c_str()),
                    const_cast<char *>("--npy"),
                    const_cast<char *>(std::to_string(npy).c_str()),
                    const_cast<char *>("--npz"),
                    const_cast<char *>(std::to_string(npz).c_str()),
                    const_cast<char *>("--nsz"),
                    const_cast<char *>(std::to_string(nsz).c_str()),
                    const_cast<char *>("--npt"),
                    const_cast<char *>(std::to_string(npt).c_str()),
                    nullptr};

    // Set session for coarse solver
    m_sessionCoarse = LibUtilities::SessionReader::CreateInstance(
        11, argv, coarseSolverFilenames, m_session->GetComm());

    // Set graph for coarse solver
    m_graphCoarse = SpatialDomains::MeshGraph::Read(m_sessionCoarse);

    // If a coarse solver session file is not specified, use
    // m_coarseSolveFactor to determine the timestep of the coarse solver
    // (default value is 100.0)
    if (!f.good())
    {
        double TimeStep =
            m_session->GetParameter("TimeStep") * m_coarseSolveFactor;
        int NumSteps =
            m_session->GetParameter("NumSteps") / m_coarseSolveFactor;
        m_sessionCoarse->SetParameter("TimeStep", TimeStep);
        m_sessionCoarse->SetParameter("NumSteps", NumSteps);
    }
}

void DriverParareal::RunCoarseSolve(
    const NekDouble time,
    const Array<OneD, const Array<OneD, NekDouble>> &input,
    Array<OneD, Array<OneD, NekDouble>> &output)
{
    // Output coarse timestep.
    if (m_comm->GetRank() == 0)
    {
        std::cout << "RUNNING COARSE SOLVE: dt = " << m_coarseTimeStep
                  << " nsteps = " << m_coarseSteps / m_numChunks << std::endl
                  << std::flush;
    }

    // Set to coarse timestep.
    m_equ[1]->SetTime(time);
    m_equ[1]->SetSteps(m_coarseSteps / m_numChunks);

    // Copy initial condition from input.
    if (m_equ[0]->GetNpoints() == m_equ[1]->GetNpoints())
    {
        // Interpolation not necessary, directly copy data
        for (int i = 0; i < m_equ[1]->GetNvariables(); ++i)
        {
            m_equ[1]->CopyToPhysField(i, input[i]);
        }
    }
    else
    {
        // Copy data to fine solver and interpolate from fine to coarse (WIP)
        for (int i = 0; i < m_equ[0]->GetNvariables(); ++i)
        {
            m_equ[0]->CopyToPhysField(i, input[i]);
        }
        m_interp.Interpolate(m_equ[0]->UpdateFields(),
                             m_equ[1]->UpdateFields());
    }

    // Solve equations.
    m_equ[1]->DoSolve();

    // Copy solution to output.
    if (m_equ[0]->GetNpoints() == m_equ[1]->GetNpoints())
    {
        // Interpolation not necessary, directly copy data
        for (int i = 0; i < m_equ[1]->GetNvariables(); ++i)
        {
            m_equ[1]->CopyFromPhysField(i, output[i]);
        }
    }
    else
    {
        // Copy data from coarse solver and interpolate coarse to fine (WIP)
        m_interp.Interpolate(m_equ[1]->UpdateFields(),
                             m_equ[0]->UpdateFields());
        for (int i = 0; i < m_equ[1]->GetNvariables(); ++i)
        {
            m_equ[0]->CopyFromPhysField(i, output[i]);
        }
    }
}

void DriverParareal::RunFineSolve(
    const NekDouble time,
    const Array<OneD, const Array<OneD, NekDouble>> &input,
    Array<OneD, Array<OneD, NekDouble>> &output)
{
    // Output fine timestep.
    if (m_comm->GetRank() == 0)
    {
        std::cout << "RUNNING FINE SOLVE: dt = " << m_fineTimeStep
                  << " nsteps = " << m_fineSteps / m_numChunks << std::endl
                  << std::flush;
    }

    // Set to fine timestep.
    m_equ[0]->SetTime(time);
    m_equ[0]->SetSteps(m_fineSteps / m_numChunks);

    // Copy initial condition from input.
    for (int i = 0; i < m_equ[0]->GetNvariables(); ++i)
    {
        m_equ[0]->CopyToPhysField(i, input[i]);
    }

    // Solve equations.
    m_equ[0]->DoSolve();

    // Copy solution to output.
    for (int i = 0; i < m_equ[0]->GetNvariables(); ++i)
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

    // Set parameters from session file.
    m_pararealIterMax= m_session->DefinesParameter("PararealIterMax")
                            ? m_session->GetParameter("PararealIterMax")
                            : m_numChunks;
    m_fineTimeStep   = m_session->GetParameter("TimeStep");
    m_coarseTimeStep = m_sessionCoarse->GetParameter("TimeStep");
    m_fineSteps      = m_session->GetParameter("NumSteps");
    m_coarseSteps    = m_sessionCoarse->GetParameter("NumSteps");
    m_totalTime      = m_fineTimeStep * m_fineSteps;
    m_chunkTime      = m_totalTime / m_numChunks;

    ASSERTL0(m_fineSteps % m_numChunks == 0,
             "Total step size should be divisible by number of chunks.");

    ASSERTL0(m_fineSteps % m_coarseSteps == 0,
             "number of coarse steps should divide number of total steps");

    // Fine solver summary
    if (m_comm->GetRank() == 0)
    {
        std::cout << "========================================================="
                     "=============="
                  << std::endl
                  << std::flush;
        std::cout << "========================= FINE PROPAGATOR INFO "
                     "========================"
                  << std::endl
                  << std::flush;
        m_equ[0]->PrintSummary(out);
    }

    std::cout << std::endl << std::flush;

    // Coarse solver summary
    if (m_comm->GetRank() == 0)
    {
        std::cout << "========================================================="
                     "=============="
                  << std::endl
                  << std::flush;
        std::cout << "======================== COARSE PROPAGATOR INFO "
                     "======================="
                  << std::endl
                  << std::flush;
        m_equ[1]->PrintSummary(out);
    }

    time(&starttime);

    // Allocate storage for parareal solver
    int nPts = m_equ[0]->GetNpoints();
    int nVar = m_equ[0]->GetNvariables();
    Array<OneD, Array<OneD, NekDouble>> solution(nVar);
    Array<OneD, Array<OneD, NekDouble>> solutionFine(nVar);
    Array<OneD, Array<OneD, NekDouble>> solutionCoarse1(nVar);
    Array<OneD, Array<OneD, NekDouble>> solutionCoarse2(nVar);
    Array<OneD, Array<OneD, NekDouble>> ic(nVar);
    Array<OneD, Array<OneD, NekDouble>> exactsoln(nVar);
    for (int i = 0; i < nVar; ++i)
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
    m_equ[1]->DoInitialise();

    // Get initial conditions
    if (m_chunkRank == 0)
    {
        for (int i = 0; i < nVar; ++i)
        {
            m_equ[0]->CopyFromPhysField(i, ic[i]);
        }
    }

    // Evaluate "ExactSolution" function, or zero array
    if (m_chunkRank == m_numChunks - 1)
    {
        for (int i = 0; i < nVar; ++i)
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
        for (int i = 0; i < nVar; ++i)
        {
            tc->Recv(recvProc, ic[i]);
        }
    }

    if (m_chunkRank == 0 && m_comm->GetRank() == 0)
    {
        std::cout << "** ITERATION " << 0 << " **" << std::endl
                  << std::flush;
    }

    RunCoarseSolve(m_chunkRank * m_chunkTime, ic, solution);
    if (m_chunkRank < m_numChunks - 1)
    {
        for (int i = 0; i < nVar; ++i)
        {
            tc->Send(sendProc, solution[i]);
        }
    }

    // On the last time-slice, we calculate the L2 error of our coarse solution
    if (m_chunkRank == m_numChunks - 1)
    {
        for (int i = 0; i < nVar; ++i)
        {
            // Copy the calculated coarse solution back
            m_equ[0]->CopyToPhysField(i, solution[i]);

            NekDouble vL2Error   = m_equ[0]->L2Error(i, exactsoln[i]);
            NekDouble vLinfError = m_equ[0]->LinfError(i, exactsoln[i]);

            std::cout << "L 2 error (variable " << m_equ[0]->GetVariable(i)
                      << ") : " << vL2Error << std::endl
                      << std::flush;
            std::cout << "L inf error (variable " << m_equ[0]->GetVariable(i)
                      << ") : " << vLinfError << std::endl
                      << std::flush;
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
            std::cout << "** ITERATION " << k + 1 << " **" << std::endl
                      << std::flush;
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
            for (int i = 0; i < nVar; ++i)
            {
                tc->Recv(recvProc, ic[i]);
            }
        }

        // Run the coarse solver
        RunCoarseSolve(m_chunkRank * m_chunkTime, ic, solutionCoarse2);

        // Calculate the new approximation y_{j+1}^{k+1}
        // This is calculated point-wise.
        for (int i = 0; i < nVar; ++i)
        {
            for (int q = 0; q < nPts; ++q)
            {
                solution[i][q] = solutionCoarse2[i][q] + solutionFine[i][q] -
                                 solutionCoarse1[i][q];
            }
        }

        // All but the last time slice should communicate the solution to the next time slice.
        // This will become the initial condition for the next slice.
        if (m_chunkRank < m_numChunks - 1)
        {
            for (int i = 0; i < nVar; ++i)
            {
                tc->Send(sendProc, solution[i]);
            }
        }

        // On the last time-slice, we calculate the L2 error of our latest approximation
        if (m_chunkRank == m_numChunks - 1)
        {
            for (int i = 0; i < nVar; ++i)
            {
                // Copy the new approximation back
                m_equ[0]->CopyToPhysField(i, solution[i]);

                NekDouble vL2Error   = m_equ[0]->L2Error(i, exactsoln[i]);
                NekDouble vLinfError = m_equ[0]->LinfError(i, exactsoln[i]);

                std::cout << "L 2 error (variable " << m_equ[0]->GetVariable(i)
                          << ") : " << vL2Error << std::endl
                          << std::flush;
                std::cout << "L inf error (variable "
                          << m_equ[0]->GetVariable(i) << ") : " << vLinfError
                          << std::endl
                          << std::flush;
            }
        }
    }

    time(&endtime);

    m_equ[0]->Output();

    if (m_comm->GetRank() == 0)
    {
        CPUtime = difftime(endtime, starttime);
        std::cout << "-------------------------------------------" << std::endl
                  << std::flush;
        std::cout << "Total Computation Time = " << CPUtime << "s" << std::endl
                  << std::flush;
        std::cout << "-------------------------------------------" << std::endl
                  << std::flush;
    }

}
}
}

