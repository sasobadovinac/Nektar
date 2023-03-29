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

#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/DriverParareal.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
string DriverParareal::className = GetDriverFactory().RegisterCreatorFunction(
    "Parareal", DriverParareal::create);
string DriverParareal::driverLookupId =
    LibUtilities::SessionReader::RegisterEnumValue("Driver", "Parareal", 0);

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
DriverParareal::~DriverParareal()
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
        ASSERTL0(
            GetEquationSystemFactory().ModuleExists(vEquation),
            "EquationSystem '" + vEquation +
                "' is not defined.\n"
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
                // Set coarse parareal session file.
                SetPararealSessionFile();

                // Set fine parareal solver.
                m_session->SetTag("AdvectiveType", "Convective");
                m_session->SetTag("PararealSolver", "FineSolver");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);

                // Set coarse parareal solver.
                m_sessionCoarse->SetTag("AdvectiveType", "Convective");
                m_sessionCoarse->SetTag("PararealSolver", "CoarseSolver");
                m_equ[1] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_sessionCoarse, m_graphCoarse);
                break;
            case eDirect:
                // Set coarse parareal session file.
                SetPararealSessionFile();

                // Set fine parareal solver.
                m_session->SetTag("AdvectiveType", "Linearised");
                m_session->SetTag("PararealSolver", "FineSolver");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);

                // Set coarse parareal solver.
                m_sessionCoarse->SetTag("AdvectiveType", "Linearised");
                m_sessionCoarse->SetTag("PararealSolver", "CoarseSolver");
                m_equ[1] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_sessionCoarse, m_graphCoarse);
                break;
            case eAdjoint:
                // Set coarse parareal session file.
                SetPararealSessionFile();

                // Set fine parareal solver.
                m_session->SetTag("AdvectiveType", "Adjoint");
                m_session->SetTag("PararealSolver", "FineSolver");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);

                // Set coarse parareal solver.
                m_sessionCoarse->SetTag("AdvectiveType", "Adjoint");
                m_sessionCoarse->SetTag("PararealSolver", "CoarseSolver");
                m_equ[1] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_sessionCoarse, m_graphCoarse);
                break;
            case eSkewSymmetric:
                // Set coarse parareal session file.
                SetPararealSessionFile();

                // Set fine parareal solver.
                m_session->SetTag("AdvectiveType", "SkewSymmetric");
                m_session->SetTag("PararealSolver", "FineSolver");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);

                // Set coarse parareal solver.
                m_sessionCoarse->SetTag("AdvectiveType", "SkewSymmetric");
                m_sessionCoarse->SetTag("PararealSolver", "CoarseSolver");
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
    // Get the coarse solver session file.
    string meshFile;
    string coarseSolverFile;
    vector<string> coarseSolverFilenames;
    bool opt         = (m_session->GetFilenames()[0].substr(
                    m_session->GetFilenames()[0].size() - 3) == "opt");
    meshFile         = m_session->GetFilenames()[0 + opt];
    coarseSolverFile = m_session->GetFilenames().size() > 1 + opt
                           ? m_session->GetFilenames()[1 + opt]
                           : m_session->GetFilenames()[0 + opt];
    coarseSolverFile = coarseSolverFile.substr(0, coarseSolverFile.size() - 4);
    coarseSolverFile += "_coarseSolver.xml";
    std::ifstream f(coarseSolverFile);

    if (f.good())
    {
        // if _coarseSolver.xml exit, read session file
        if (m_session->GetFilenames().size() > 1 + opt)
        {
            coarseSolverFilenames.push_back(meshFile);
        }
        coarseSolverFilenames.push_back(coarseSolverFile);
    }
    else
    {
        // if _coarseSolver.xml does not exit, use original session file
        coarseSolverFilenames.push_back(m_session->GetFilenames()[0 + opt]);
        if (m_session->GetFilenames().size() > 1 + opt)
        {
            coarseSolverFilenames.push_back(m_session->GetFilenames()[1 + opt]);
        }
    }

    // Define argument for the coarse parareal solver.
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

    // Convert into string.
    std::string npx_string = std::to_string(npx);
    std::string npy_string = std::to_string(npy);
    std::string npz_string = std::to_string(npz);
    std::string nsz_string = std::to_string(nsz);
    std::string npt_string = std::to_string(npt);

    char *argv[] = {const_cast<char *>("Solver"), // this is just a place holder
                    const_cast<char *>("--npx"),
                    const_cast<char *>(npx_string.c_str()),
                    const_cast<char *>("--npy"),
                    const_cast<char *>(npy_string.c_str()),
                    const_cast<char *>("--npz"),
                    const_cast<char *>(npz_string.c_str()),
                    const_cast<char *>("--nsz"),
                    const_cast<char *>(nsz_string.c_str()),
                    const_cast<char *>("--npt"),
                    const_cast<char *>(npt_string.c_str()),
                    nullptr};

    // Set session for coarse solver.
    m_sessionCoarse = LibUtilities::SessionReader::CreateInstance(
        11, argv, coarseSolverFilenames, m_session->GetComm());

    // Set graph for coarse solver.
    m_graphCoarse = SpatialDomains::MeshGraph::Read(m_sessionCoarse);

    // Set BndRegionOrdering (necessary for DG with periodic BC) FIXME
    m_graphCoarse->SetBndRegionOrdering(m_graph->GetBndRegionOrdering());

    // Set CompositeOrdering (necessary for DG with periodic BC) FIXME
    m_graphCoarse->SetCompositeOrdering(m_graph->GetCompositeOrdering());

    // If a coarse solver session file is not specified, use
    // m_coarseSolveFactor to determine the timestep of the coarse solver
    // (default value is 100.0).
    if (!f.good())
    {
        double timeStep =
            m_session->GetParameter("TimeStep") * m_coarseSolveFactor;
        int numSteps =
            m_session->GetParameter("NumSteps") / m_coarseSolveFactor;
        m_sessionCoarse->SetParameter("TimeStep", timeStep);
        m_sessionCoarse->SetParameter("NumSteps", numSteps);
    }
}

void DriverParareal::RunCoarseSolve(
    const NekDouble time, const int nstep, const int iter,
    const Array<OneD, const Array<OneD, NekDouble>> &input,
    Array<OneD, Array<OneD, NekDouble>> &output)
{
    // Output coarse timestep.
    if (m_chunkRank == iter && m_comm->GetSpaceComm()->GetRank() == 0)
    {
        std::cout << "RUNNING COARSE SOLVE: dt = " << m_coarseTimeStep
                  << " nsteps = " << m_coarseSteps / m_numChunks << std::endl
                  << std::flush;
    }

    // Set to coarse timestep.
    m_equ[1]->SetTime(time);
    m_equ[1]->SetSteps(nstep);

    // Copy initial condition from input.
    if (m_equ[0]->GetNpoints() == m_equ[1]->GetNpoints())
    {
        // Interpolation not necessary, directly copy data.
        for (int i = 0; i < m_equ[1]->GetNvariables(); ++i)
        {
            m_equ[1]->CopyToPhysField(i, input[i]);
        }
    }
    else
    {
        // Copy data to fine solver and interpolate from fine to coarse (WIP).
        for (int i = 0; i < m_equ[0]->GetNvariables(); ++i)
        {
            m_equ[0]->CopyToPhysField(i, input[i]);
        }
        m_interp.InterpExp1ToExp2(m_equ[0]->UpdateFields(),
                                  m_equ[1]->UpdateFields());
    }

    // Solve equations.
    m_equ[1]->DoSolve();

    // Copy solution to output.
    if (m_equ[0]->GetNpoints() == m_equ[1]->GetNpoints())
    {
        // Interpolation not necessary, directly copy data.
        for (int i = 0; i < m_equ[1]->GetNvariables(); ++i)
        {
            m_equ[1]->CopyFromPhysField(i, output[i]);
        }
    }
    else
    {
        // Copy data from coarse solver and interpolate coarse to fine (WIP).
        m_interp.InterpExp1ToExp2(m_equ[1]->UpdateFields(),
                                  m_equ[0]->UpdateFields());
        for (int i = 0; i < m_equ[1]->GetNvariables(); ++i)
        {
            m_equ[0]->CopyFromPhysField(i, output[i]);
        }
    }
}

void DriverParareal::RunFineSolve(
    const NekDouble time, const int nstep, const int iter,
    const Array<OneD, const Array<OneD, NekDouble>> &input,
    Array<OneD, Array<OneD, NekDouble>> &output)
{
    // Output fine timestep.
    if (m_chunkRank == iter && m_comm->GetSpaceComm()->GetRank() == 0)
    {
        std::cout << "RUNNING FINE SOLVE: dt = " << m_fineTimeStep
                  << " nsteps = " << m_fineSteps / m_numChunks << std::endl
                  << std::flush;
    }

    // Number of checkpoint by chunk.
    int nChkPts =
        m_session->GetParameter("IO_CheckSteps")
            ? m_fineSteps /
                  int(m_session->GetParameter("IO_CheckSteps") * m_numChunks)
            : 1;

    // Parareal iteration number.
    int nIter = m_equ[0]->GetPararealIterationNumber();

    // Set to fine timestep.
    m_equ[0]->SetTime(time);
    m_equ[0]->SetSteps(nstep);

    // Reinitialize check point number for each parareal iteration.
    m_equ[0]->SetCheckpointNumber(m_chunkRank * nChkPts + 1);

    // Update parareal iteration number.
    m_equ[0]->SetPararealIterationNumber(++nIter);

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
    Nektar::LibUtilities::Timer timer;
    NekDouble CPUtime = 0.0;

    m_numChunks = m_session->GetComm()->GetTimeComm()->GetSize();
    m_chunkRank = m_session->GetComm()->GetTimeComm()->GetRank();

    // Set parameters from session file.
    m_pararealToler   = m_session->DefinesParameter("PararealToler")
                            ? m_session->GetParameter("PararealToler")
                            : m_pararealToler;
    m_pararealIterMax = m_session->DefinesParameter("PararealIterMax")
                            ? m_session->GetParameter("PararealIterMax")
                            : m_numChunks;
    m_fineTimeStep    = m_session->GetParameter("TimeStep");
    m_coarseTimeStep  = m_sessionCoarse->GetParameter("TimeStep");
    m_fineSteps       = m_session->GetParameter("NumSteps");
    m_coarseSteps     = m_sessionCoarse->GetParameter("NumSteps");
    m_totalTime       = m_fineTimeStep * m_fineSteps;
    m_chunkTime       = m_totalTime / m_numChunks;
    m_exactSolution   = m_session->DefinesParameter("ExactSolution")
                            ? m_session->GetParameter("ExactSolution")
                            : m_exactSolution;

    // Turnoff I/O for coarse solver.
    m_equ[1]->SetInfoSteps(0);
    m_equ[1]->SetCheckpointSteps(0);

    // Check time step inputs.
    ASSERTL0(
        m_fineSteps % m_numChunks == 0,
        "Total number of fine step should be divisible by number of chunks.");

    ASSERTL0(
        m_coarseSteps % m_numChunks == 0,
        "Total number of coarse step should be divisible by number of chunks.");

    ASSERTL0(fabs(m_coarseTimeStep * m_coarseSteps -
                  m_fineTimeStep * m_fineSteps) < 1e-12,
             "Fine and coarse total computational times do not match");

    // Check I/O inputs.
    if (m_session->GetParameter("IO_InfoSteps"))
    {
        ASSERTL0(m_fineSteps % int(m_session->GetParameter("IO_InfoSteps") *
                                   m_numChunks) ==
                     0,
                 "number of IO_InfoSteps should divide number of fine steps "
                 "per time chunk");
    }
    if (m_session->GetParameter("IO_CheckSteps"))
    {
        ASSERTL0(m_fineSteps % int(m_session->GetParameter("IO_CheckSteps") *
                                   m_numChunks) ==
                     0,
                 "number of IO_CheckSteps should divide number of fine steps "
                 "per time chunk");
    }

    // Fine solver summary.
    if (m_chunkRank == 0 && m_comm->GetSpaceComm()->GetRank() == 0)
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

        std::cout << std::endl << std::flush;
    }

    // Coarse solver summary.
    if (m_chunkRank == 0 && m_comm->GetSpaceComm()->GetRank() == 0)
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

    timer.Start();

    // Allocate storage for parareal solver.
    int nPts = m_equ[0]->GetNpoints();
    int nVar = m_equ[0]->GetNvariables();
    Array<OneD, Array<OneD, NekDouble>> solution(nVar);
    Array<OneD, Array<OneD, NekDouble>> solutionCoarsePrev(nVar);
    Array<OneD, Array<OneD, NekDouble>> solutionCoarseCurr(nVar);
    Array<OneD, Array<OneD, NekDouble>> ic(nVar);
    Array<OneD, Array<OneD, NekDouble>> exactsoln(nVar);
    for (int i = 0; i < nVar; ++i)
    {
        solution[i]           = Array<OneD, NekDouble>(nPts, 0.0);
        solutionCoarsePrev[i] = Array<OneD, NekDouble>(nPts, 0.0);
        solutionCoarseCurr[i] = Array<OneD, NekDouble>(nPts, 0.0);
        ic[i]                 = Array<OneD, NekDouble>(nPts, 0.0);
        exactsoln[i]          = Array<OneD, NekDouble>(nPts, 0.0);
    }

    // Initialize fine solver.
    m_equ[0]->DoInitialise();

    // Initialize coarse solver.
    m_equ[1]->DoInitialise(false);

    // Get initial conditions.
    for (int i = 0; i < nVar; ++i)
    {
        m_equ[0]->CopyFromPhysField(i, ic[i]);
    }

    // Run coarse solution, G(y_j^k) to get initial conditions.
    if (m_chunkRank == 0 && m_comm->GetSpaceComm()->GetRank() == 0)
    {
        std::cout << "** INITIAL CONDITION **" << std::endl << std::flush;
    }
    LibUtilities::CommSharedPtr tComm = m_session->GetComm()->GetTimeComm();
    int recvProc                      = m_chunkRank - 1;
    int sendProc                      = m_chunkRank + 1;

    // Iterate to improve on the approximation
    // For the moment, we use the maximum number of iterations
    // We can add a tolerance threshold later.
    // We start the iteration with the current approximation stored in
    // the 'solution' array.
    int k               = 0;
    int kmax            = 0;
    int convergenceCurr = false;
    int convergencePrev = (m_chunkRank == 0);
    NekDouble vL2Error  = 0.0;
    Array<OneD, NekDouble> vL2Errors(nVar, 0.0);
    Array<OneD, NekDouble> vLinfErrors(nVar, 0.0);

    // To compute the initial conditions, the coarse solver is run serially
    // for each time chunk to avoid communication time between processors.
    if (m_chunkRank > 0)
    {
        RunCoarseSolve(0.0, m_chunkRank * m_coarseSteps / m_numChunks, -1, ic,
                       ic);
    }

    // Compute initial coarse solution.
    if (m_chunkRank == 0 && m_comm->GetSpaceComm()->GetRank() == 0)
    {
        std::cout << "** ITERATION " << 0 << " **" << std::endl << std::flush;
    }

    RunCoarseSolve(m_chunkRank * m_chunkTime, m_coarseSteps / m_numChunks, 0,
                   ic, solutionCoarsePrev);

    for (int i = 0; i < nVar; ++i)
    {
        for (int q = 0; q < nPts; ++q)
        {
            solution[i][q] = solutionCoarsePrev[i][q];
        }
    }

    if (m_exactSolution)
    {
        // Evaluate exact solution.
        for (int i = 0; i < nVar; ++i)
        {
            m_equ[0]->EvaluateExactSolution(i, exactsoln[i],
                                            (m_chunkRank + 1) * m_chunkTime);
        }
    }

    timer.Stop();
    CPUtime += timer.Elapsed().count();
    if (m_chunkRank == m_numChunks - 1 &&
        m_comm->GetSpaceComm()->GetRank() == 0)
    {
        std::cout << "-------------------------------------------" << std::endl
                  << std::flush;
        std::cout << "Total Computation Time = " << CPUtime << "s" << std::endl
                  << std::flush;
        std::cout << "-------------------------------------------" << std::endl
                  << std::flush;
    }
    for (int i = 0; i < nVar; ++i)
    {
        vL2Errors[i]   = m_equ[1]->L2Error(i, exactsoln[i], 1);
        vLinfErrors[i] = m_equ[1]->LinfError(i, exactsoln[i]);
        if (m_chunkRank == m_numChunks - 1 &&
            m_comm->GetSpaceComm()->GetRank() == 0)
        {
            std::cout << "L2 error (variable " << m_equ[1]->GetVariable(i)
                      << ") : " << vL2Errors[i] << std::endl
                      << std::flush;
            std::cout << "Linf error (variable " << m_equ[1]->GetVariable(i)
                      << ") : " << vLinfErrors[i] << std::endl
                      << std::flush;
        }
    }
    timer.Start();

    // Start Parareal iteration loop.
    while (k < m_pararealIterMax && !convergenceCurr)
    {
        if (m_chunkRank == min(k, m_numChunks - 1) &&
            m_comm->GetSpaceComm()->GetRank() == 0)
        {
            std::cout << "** ITERATION " << k + 1 << " **" << std::endl
                      << std::flush;
        }

        // Use previous parareal solution as "exact solution".
        if (!m_exactSolution)
        {
            for (int i = 0; i < nVar; ++i)
            {
                for (int q = 0; q < nPts; ++q)
                {
                    exactsoln[i][q] = solution[i][q];
                }
            }
        }

        // Calculate fine solution, F(y_j^k)
        // Again no communication necessary.
        RunFineSolve(m_chunkRank * m_chunkTime, m_fineSteps / m_numChunks, k,
                     ic, solution);

        // Calculate coarse solve correction G(y_j^{k+1})
        // These are dependent on the previous time slice, so need to
        // compute serially.
        if (m_chunkRank > 0 && !convergencePrev)
        {
            // All time slices, apart from the first, receive their initial
            // state from the previous time slice.
            tComm->Recv(recvProc, convergencePrev);
            for (int i = 0; i < nVar; ++i)
            {
                tComm->Recv(recvProc, ic[i]);
            }
        }

        // Run the coarse solver
        RunCoarseSolve(m_chunkRank * m_chunkTime, m_coarseSteps / m_numChunks,
                       k, ic, solutionCoarseCurr);

        // Calculate the new approximation y_{j+1}^{k+1}
        // This is calculated point-wise.
        for (int i = 0; i < nVar; ++i)
        {
            for (int q = 0; q < nPts; ++q)
            {
                solution[i][q] +=
                    solutionCoarseCurr[i][q] - solutionCoarsePrev[i][q];

                solutionCoarsePrev[i][q] = solutionCoarseCurr[i][q];
            }
        }

        // Compute L2 error for each time chunk.
        vL2Error = 0.0;
        for (int i = 0; i < nVar; ++i)
        {
            // Copy the new approximation back.
            m_equ[0]->CopyToPhysField(i, solution[i]);

            vL2Error = max(vL2Error, m_equ[0]->L2Error(i, exactsoln[i], 1));
        }

        // Check convergence of L2 error for each time chunk.
        if ((vL2Error < m_pararealToler && convergencePrev) || m_chunkRank == k)
        {
            convergenceCurr = true;
        }

        // All but the last time slice should communicate the solution to
        // the next time slice. This will become the initial condition for
        // the next slice.
        if (m_chunkRank < m_numChunks - 1)
        {
            tComm->Send(sendProc, convergenceCurr);
            for (int i = 0; i < nVar; ++i)
            {
                tComm->Send(sendProc, solution[i]);
            }
        }

        if (m_exactSolution)
        {
            // Evaluate exact solution.
            for (int i = 0; i < nVar; ++i)
            {
                m_equ[0]->EvaluateExactSolution(
                    i, exactsoln[i], (m_chunkRank + 1) * m_chunkTime);
            }
        }

        timer.Stop();
        CPUtime += timer.Elapsed().count();
        if (m_chunkRank == m_numChunks - 1 &&
            m_comm->GetSpaceComm()->GetRank() == 0)
        {
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
            std::cout << "Total Computation Time = " << CPUtime << "s"
                      << std::endl
                      << std::flush;
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
        }
        for (int i = 0; i < nVar; ++i)
        {
            vL2Errors[i]   = m_equ[0]->L2Error(i, exactsoln[i], 1);
            vLinfErrors[i] = m_equ[0]->LinfError(i, exactsoln[i]);
            if (m_chunkRank == m_numChunks - 1 &&
                m_comm->GetSpaceComm()->GetRank() == 0)
            {
                std::cout << "L2 error (variable " << m_equ[0]->GetVariable(i)
                          << ") : " << vL2Errors[i] << std::endl
                          << std::flush;
                std::cout << "Linf error (variable " << m_equ[0]->GetVariable(i)
                          << ") : " << vLinfErrors[i] << std::endl
                          << std::flush;
            }
        }
        timer.Start();

        // Increment counter.
        k++;
        kmax = k;
    }
    timer.Stop();

    // If already converged, simply copy previous computed solution.
    m_comm->GetTimeComm()->AllReduce(kmax, Nektar::LibUtilities::ReduceMax);
    for (; k < kmax; k++)
    {
        if (m_comm->GetSpaceComm()->GetRank() == 0 &&
            m_session->GetParameter("IO_CheckSteps"))
        {
            std::string olddir = m_equ[0]->GetSessionName() + "_" +
                                 boost::lexical_cast<std::string>(k) + ".pit";
            std::string outdir = m_equ[0]->GetSessionName() + "_" +
                                 boost::lexical_cast<std::string>(k + 1) +
                                 ".pit";

            int nChkPts =
                m_fineSteps /
                int(m_session->GetParameter("IO_CheckSteps") * m_numChunks);
            for (int i = 0; i < nChkPts; i++)
            {
                // Old file name
                std::string oldname =
                    olddir + "/" + m_equ[0]->GetSessionName() + "_" +
                    boost::lexical_cast<std::string>(m_chunkRank * nChkPts + i +
                                                     1) +
                    ".chk";

                // New file name
                std::string outname =
                    outdir + "/" + m_equ[0]->GetSessionName() + "_" +
                    boost::lexical_cast<std::string>(m_chunkRank * nChkPts + i +
                                                     1) +
                    ".chk";

                // Remove file if already existing
                fs::remove_all(outname);

                // Copy from previous converged solution
                fs::copy(oldname, outname);
            }
        }
    }

    // Update solution before printing restart solution.
    for (int i = 0; i < nVar; ++i)
    {
        m_equ[0]->CopyToPhysField(i, solution[i]);
        m_equ[0]->UpdateFields()[i]->FwdTransLocalElmt(
            m_equ[0]->UpdateFields()[i]->GetPhys(),
            m_equ[0]->UpdateFields()[i]->UpdateCoeffs());
    }

    // Print restart solution.
    m_equ[0]->Output();

    // Wait for all processors to finish their writing activities
    m_comm->Block();

    // Print total computational time.
    if (m_chunkRank == m_numChunks - 1)
    {
        if (m_comm->GetSpaceComm()->GetRank() == 0)
        {
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
            std::cout << "Total Computation Time = " << CPUtime << "s"
                      << std::endl
                      << std::flush;
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
        }

        // Print solution errors.
        for (int i = 0; i < nVar; ++i)
        {
            // Copy the new approximation back.
            m_equ[0]->CopyToPhysField(i, solution[i]);

            // Evaluate exact solution.
            m_equ[0]->EvaluateExactSolution(i, exactsoln[i], m_totalTime);

            // Evaluate error norms.
            NekDouble vL2Error   = m_equ[0]->L2Error(i, exactsoln[i]);
            NekDouble vLinfError = m_equ[0]->LinfError(i, exactsoln[i]);

            if (m_comm->GetSpaceComm()->GetRank() == 0)
            {
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

    // Speed-up analysis
    if (true)
    {
        int count = 20;

        if (m_chunkRank == m_numChunks - 1 &&
            m_comm->GetSpaceComm()->GetRank() == 0)
        {
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
            std::cout << "PARAREAL SPEED-UP ANALYSIS" << std::endl
                      << std::flush;
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
        }

        // Mean communication time.
        NekDouble commTime = 0.0;
        timer.Start();
        for (int n = 0; n < count; n++)
        {
            if (m_chunkRank > 0)
            {
                for (int i = 0; i < nVar; ++i)
                {
                    tComm->Recv(recvProc, ic[i]);
                }
            }

            if (m_chunkRank < m_numChunks - 1)
            {
                for (int i = 0; i < nVar; ++i)
                {
                    tComm->Send(sendProc, solution[i]);
                }
            }
        }
        timer.Stop();
        commTime = timer.Elapsed().count() / count;

        // Print communication time.
        if (m_chunkRank == m_numChunks - 1 &&
            m_comm->GetSpaceComm()->GetRank() == 0)
        {
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
            std::cout << "Mean Communication Time = " << commTime << "s"
                      << std::endl
                      << std::flush;
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
        }

        // Mean restriction time.
        NekDouble restTime = 0.0;
        if (m_equ[0]->GetNpoints() != m_equ[1]->GetNpoints())
        {
            timer.Start();
            for (int n = 0; n < count; n++)
            {
                m_interp.InterpExp1ToExp2(m_equ[0]->UpdateFields(),
                                          m_equ[1]->UpdateFields());
            }
            timer.Stop();
            restTime = timer.Elapsed().count() / count;

            // Print restriction time.
            if (m_chunkRank == m_numChunks - 1 &&
                m_comm->GetSpaceComm()->GetRank() == 0)
            {
                std::cout << "-------------------------------------------"
                          << std::endl
                          << std::flush;
                std::cout << "Mean Restriction Time = " << restTime << "s"
                          << std::endl
                          << std::flush;
                std::cout << "-------------------------------------------"
                          << std::endl
                          << std::flush;
            }
        }

        // Mean interpolation time.
        NekDouble interTime = 0.0;
        if (m_equ[0]->GetNpoints() != m_equ[1]->GetNpoints())
        {
            timer.Start();
            for (int n = 0; n < count; n++)
            {
                m_interp.InterpExp1ToExp2(m_equ[1]->UpdateFields(),
                                          m_equ[0]->UpdateFields());
            }
            timer.Stop();
            interTime = timer.Elapsed().count() / count;

            // Print restriction time.
            if (m_chunkRank == m_numChunks - 1 &&
                m_comm->GetSpaceComm()->GetRank() == 0)
            {
                std::cout << "-------------------------------------------"
                          << std::endl
                          << std::flush;
                std::cout << "Mean Interpolation Time = " << interTime << "s"
                          << std::endl
                          << std::flush;
                std::cout << "-------------------------------------------"
                          << std::endl
                          << std::flush;
            }
        }

        // Mean coarse solver time.
        NekDouble coarseSolveTime = 0.0;
        timer.Start();
        for (int n = 0; n < count; n++)
        {
            RunCoarseSolve(0.0, 100, -1, ic, solution);
        }
        timer.Stop();
        coarseSolveTime = 0.01 * timer.Elapsed().count() / count *
                              (m_coarseSteps / m_numChunks) -
                          restTime - interTime;

        // Print restriction time.
        if (m_chunkRank == m_numChunks - 1 &&
            m_comm->GetSpaceComm()->GetRank() == 0)
        {
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
            std::cout << "Mean Coarse Solve Time = " << coarseSolveTime << "s"
                      << std::endl
                      << std::flush;
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
        }

        // Fine solver time.
        NekDouble fineSolveTime = 0.0;
        timer.Start();
        // Turnoff I/O for fine solver.
        m_equ[0]->SetInfoSteps(0);
        m_equ[0]->SetCheckpointSteps(0);
        for (int n = 0; n < count; n++)
        {
            RunFineSolve(0.0, 100, -1, ic, solution);
        }
        timer.Stop();
        fineSolveTime = 0.01 * timer.Elapsed().count() / count *
                        (m_fineSteps / m_numChunks);

        // Print fine solve time.
        if (m_chunkRank == m_numChunks - 1 &&
            m_comm->GetSpaceComm()->GetRank() == 0)
        {
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
            std::cout << "Mean Fine Solve Time = " << fineSolveTime << "s"
                      << std::endl
                      << std::flush;
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
        }

        // Print speedup time.
        if (m_chunkRank == m_numChunks - 1 &&
            m_comm->GetSpaceComm()->GetRank() == 0)
        {
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
            std::cout << "Maximum Speed-up" << std::endl << std::flush;
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
            for (int k = 1; k <= m_numChunks; k++)
            {
                NekDouble ratio      = double(k) / m_numChunks;
                NekDouble ratioSolve = coarseSolveTime / fineSolveTime;
                NekDouble speedup = 1.0 / ((1.0 + ratio) * ratioSolve + ratio);
                std::cout << "Speed-up (" << k << ") = " << speedup << std::endl
                          << std::flush;
            }
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
            std::cout << "Speed-up with comm." << std::endl << std::flush;
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
            for (int k = 1; k <= m_numChunks; k++)
            {
                NekDouble ratio      = double(k) / m_numChunks;
                NekDouble ratioComm  = commTime / fineSolveTime;
                NekDouble ratioSolve = coarseSolveTime / fineSolveTime;
                NekDouble speedup = 1.0 / ((1.0 + ratio) * ratioSolve + ratio +
                                           ratioComm / m_numChunks);
                std::cout << "Speed-up (" << k << ") = " << speedup << std::endl
                          << std::flush;
            }
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
            std::cout << "Speed-up with comm. and interp." << std::endl
                      << std::flush;
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
            for (int k = 1; k <= m_numChunks; k++)
            {
                NekDouble ratio     = double(k) / m_numChunks;
                NekDouble ratioComm = commTime / fineSolveTime;
                NekDouble ratioSolve =
                    (coarseSolveTime + restTime + interTime) / fineSolveTime;
                NekDouble speedup =
                    1.0 / ((1.0 + ratio) * ratioSolve + ratio +
                           ratioComm * k * (2 * m_numChunks - k - 1) / 2.0 /
                               m_numChunks);
                std::cout << "Speed-up (" << k << ") = " << speedup << std::endl
                          << std::flush;
            }
            std::cout << "-------------------------------------------"
                      << std::endl
                      << std::flush;
        }
    }
}
} // namespace SolverUtils
} // namespace Nektar
