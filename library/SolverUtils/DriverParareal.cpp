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
#include <boost/format.hpp>

namespace Nektar
{
namespace SolverUtils
{
std::string DriverParareal::className =
    GetDriverFactory().RegisterCreatorFunction("Parareal",
                                               DriverParareal::create);
std::string DriverParareal::driverLookupId =
    LibUtilities::SessionReader::RegisterEnumValue("Driver", "Parareal", 0);

/**
 *
 */
DriverParareal::DriverParareal(
    const LibUtilities::SessionReaderSharedPtr pSession,
    const SpatialDomains::MeshGraphSharedPtr pGraph)
    : DriverParallelInTime(pSession, pGraph)
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
void DriverParareal::v_Execute(std::ostream &out)
{
    // Set timer.
    Nektar::LibUtilities::Timer timer;
    NekDouble CPUtime = 0.0;

    // Set time communication parameters.
    LibUtilities::CommSharedPtr tComm = m_session->GetComm()->GetTimeComm();
    m_numChunks                       = tComm->GetSize();
    m_chunkRank                       = tComm->GetRank();

    // Get and assert parameters from session file.
    GetParametersFromSession();
    AssertParameters();

    // Print solver info.
    PrintFineSolverInfo(out);
    PrintCoarseSolverInfo(out);

    // Initialization.
    InitialiseEqSystem(false);
    AllocateMemory();
    InitialiseInterpolationField();

    // Allocate storage for Parareal solver.
    Array<OneD, Array<OneD, NekDouble>> solutionCoarseCurr(m_nVar);
    Array<OneD, Array<OneD, NekDouble>> solutionCoarsePrev(m_nVar);
    Array<OneD, Array<OneD, NekDouble>> ic(m_nVar);
    Array<OneD, Array<OneD, NekDouble>> solution(m_nVar);
    for (size_t i = 0; i < m_nVar; ++i)
    {
        solutionCoarseCurr[i] = Array<OneD, NekDouble>(m_fineNpts, 0.0);
        solutionCoarsePrev[i] = Array<OneD, NekDouble>(m_fineNpts, 0.0);
        ic[i]                 = Array<OneD, NekDouble>(m_fineNpts, 0.0);
        solution[i]           = Array<OneD, NekDouble>(m_fineNpts, 0.0);
    }

    // Get initial condition from fields.
    CopyFromFinePhysField(ic);

    // Start iteration windows.
    tComm->Block();
    m_totalTime = m_fineTimeStep * m_fineSteps;
    m_chunkTime = m_totalTime / m_numChunks / m_numWindowsPIT;
    m_fineSteps /= m_numChunks * m_numWindowsPIT;
    m_coarseSteps /= m_numChunks * m_numWindowsPIT;
    for (size_t w = 0; w < m_numWindowsPIT; w++)
    {
        timer.Start();

        // Initialize time for the current window.
        m_time = (w * m_numChunks) * m_chunkTime;

        // Print window number.
        PrintHeaderTitle1((boost::format("WINDOWS #%1%") % (w + 1)).str());

        // Run coarse solver.
        if (m_chunkRank > 0)
        {
            RunCoarseSolve(m_time, m_chunkRank * m_coarseSteps, ic, ic);
            m_time += m_chunkRank * m_chunkTime;
        }
        RunCoarseSolve(m_time, m_coarseSteps, ic, solutionCoarsePrev);
        CopySolutionVector(solutionCoarsePrev, solution);

        // Update fields with solution.
        CopyToFinePhysField(solution);

        // Compute exact solution, if necessary.
        if (m_exactSolution)
        {
            EvaluateExactSolution(m_time + m_chunkTime);
        }

        // Solution convergence monitoring.
        timer.Stop();
        CPUtime += timer.Elapsed().count();
        PrintHeaderTitle2((boost::format("ITERATION %1%") % 0).str());
        SolutionConvergenceMonitoring(CPUtime);
        timer.Start();

        // Start Parareal iteration.
        size_t k            = 0;
        size_t kmax         = 0;
        int convergenceCurr = false;
        int convergencePrev = (m_chunkRank == 0);
        while (k < m_iterMaxPIT && !convergenceCurr)
        {
            // Use previous parareal solution as "exact solution", if necessary.
            if (!m_exactSolution)
            {
                CopySolutionVector(solution, m_exactsoln);
            }

            // Calculate fine solution (parallel-in-time).
            RunFineSolve(m_time, m_fineSteps, k, w, ic, solution);

            // Receive coarse solution from previous processor.
            RecvInitialConditionFromPreviousProc(ic, convergencePrev);

            // Calculate coarse solution (serial-in-time).
            RunCoarseSolve(m_time, m_coarseSteps, ic, solutionCoarseCurr);

            // Calculate corrected solution.
            PararealCorrection(solutionCoarseCurr, solutionCoarsePrev,
                               solution);

            // Save current coarse solution.
            CopySolutionVector(solutionCoarseCurr, solutionCoarsePrev);

            // Update fields with corrected solution.
            CopyToFinePhysField(solution);

            // Solution convergence monitoring.
            timer.Stop();
            CPUtime += timer.Elapsed().count();
            PrintHeaderTitle2((boost::format("ITERATION %1%") % (k + 1)).str());
            SolutionConvergenceMonitoring(CPUtime);
            timer.Start();

            // Check convergence of L2 error for each time chunk.
            if ((vL2ErrorMax() < m_tolerPIT && convergencePrev) ||
                m_chunkRank == k)
            {
                convergenceCurr = true;
            }

            // Send coarse solution to next processor.
            SendSolutionToNextProc(solution, convergenceCurr);

            // Increment index.
            k++;
            kmax = k;
        }
        timer.Stop();

        // Copy converged check points.
        CopyConvergedCheckPoints(w, k, kmax);

        // Update solution before printing solution.
        UpdateSolution(solution);

        // Print solution files.
        PrintSolutionFile();

        // Windowing.
        if (w != m_numWindowsPIT - 1)
        {
            ApplyWindowing(solution, ic);
        }
    }

    // Post-processing.
    tComm->Block();
    PrintHeaderTitle1("SUMMARY");
    EvaluateExactSolution(m_time + m_chunkTime);
    SolutionConvergenceSummary(CPUtime);
    SpeedUpAnalysis();
}

/**
 *
 */
NekDouble DriverParareal::v_EstimateCommunicationTime(void)
{
    // Allocate memory.
    Array<OneD, Array<OneD, NekDouble>> buffer1(m_nVar);
    Array<OneD, Array<OneD, NekDouble>> buffer2(m_nVar);
    for (size_t i = 0; i < m_nVar; ++i)
    {
        buffer1[i] = Array<OneD, NekDouble>(m_fineNpts, 0.0);
        buffer2[i] = Array<OneD, NekDouble>(m_fineNpts, 0.0);
    }

    // Estimate communication time.
    return EstimateCommunicationTime(buffer1, buffer2);
}

/**
 *
 */
NekDouble DriverParareal::v_EstimateRestrictionTime(void)
{
    if (m_fineNpts == m_coarseNpts)
    {
        return 0.0; // No restriction necessary
    }
    else
    {
        // Average restriction time over niter iteration.
        size_t niter = 20;
        Nektar::LibUtilities::Timer timer;
        timer.Start();
        for (size_t n = 0; n < niter; n++)
        {
            Interpolator(m_tmpfine, m_tmpcoarse);
        }
        timer.Stop();
        return timer.Elapsed().count() / niter;
    }
}

/**
 *
 */
NekDouble DriverParareal::v_EstimateInterpolationTime(void)
{
    if (m_fineNpts == m_coarseNpts)
    {
        return 0.0; // No interpolation
    }
    else
    {
        // Average interpolation time over niter iteration.
        size_t niter = 20;
        Nektar::LibUtilities::Timer timer;
        timer.Start();
        for (size_t n = 0; n < niter; n++)
        {
            Interpolator(m_tmpcoarse, m_tmpfine);
        }
        timer.Stop();
        return timer.Elapsed().count() / niter;
    }
}

/**
 *
 */
NekDouble DriverParareal::v_EstimateCoarseSolverTime(void)
{
    // Allocate memory.
    Array<OneD, Array<OneD, NekDouble>> sol(m_nVar);
    for (size_t i = 0; i < m_nVar; ++i)
    {
        sol[i] = Array<OneD, NekDouble>(m_coarseNpts);
    }

    // Turnoff I/O.
    m_coarseEqSys->SetInfoSteps(0);
    m_coarseEqSys->SetCheckpointSteps(0);

    // Get initial condition.
    CopyFromFinePhysField(sol);

    // Estimate coarse solver time.
    Nektar::LibUtilities::Timer timer;
    timer.Start();

    // Set to coarse timestep.
    m_coarseEqSys->SetTime(m_time + m_chunkTime);
    m_coarseEqSys->SetSteps(10);

    // Copy initial condition.
    for (size_t i = 0; i < m_fineEqSys->GetNvariables(); ++i)
    {
        m_coarseEqSys->CopyToPhysField(i, sol[i]);
    }

    // Solve equations.
    m_coarseEqSys->DoSolve();

    // Copy solution.
    for (size_t i = 0; i < m_coarseEqSys->GetNvariables(); ++i)
    {
        m_coarseEqSys->CopyFromPhysField(i, sol[i]);
    }

    timer.Stop();
    return 0.1 * timer.Elapsed().count() * m_coarseSteps;
}

/**
 *
 */
NekDouble DriverParareal::v_EstimateFineSolverTime(void)
{
    // Allocate memory.
    Array<OneD, Array<OneD, NekDouble>> sol(m_nVar);
    for (size_t i = 0; i < m_nVar; ++i)
    {
        sol[i] = Array<OneD, NekDouble>(m_fineNpts);
    }

    // Turnoff I/O.
    m_fineEqSys->SetInfoSteps(0);
    m_fineEqSys->SetCheckpointSteps(0);

    // Get initial condition.
    CopyFromFinePhysField(sol);

    // Estimate fine solver time.
    Nektar::LibUtilities::Timer timer;
    timer.Start();
    RunFineSolve(m_time + m_chunkTime, 100, 0, 0, sol, sol);
    timer.Stop();
    return 0.01 * timer.Elapsed().count() * m_fineSteps;
}

/**
 *
 */
NekDouble DriverParareal::v_EstimatePredictorTime(void)
{
    return v_EstimateCoarseSolverTime();
}

/**
 *
 */
NekDouble DriverParareal::v_EstimateOverheadTime(void)
{
    return 0.0;
}

/**
 *
 */
NekDouble DriverParareal::v_ComputeSpeedUp(
    const size_t iter, NekDouble fineSolveTime, NekDouble coarseSolveTime,
    NekDouble restTime, NekDouble interTime, NekDouble commTime,
    NekDouble predictorTime, NekDouble overheadTime)
{
    // The speed-up estimate is based on "Lunet, T., Bodart, J., Gratton, S., &
    // Vasseur, X. (2018). Time-parallel simulation of the decay of homogeneous
    // turbulence using parareal with spatial coarsening. Computing and
    // Visualization in Science, 19, 31-44".

    size_t nComm             = (iter * (2 * m_numChunks - iter - 1)) / 2;
    NekDouble ratio          = double(iter) / m_numChunks;
    NekDouble ratioPredictor = predictorTime / fineSolveTime;
    NekDouble ratioSolve     = coarseSolveTime / fineSolveTime;
    NekDouble ratioInterp    = (restTime + interTime) / fineSolveTime;
    NekDouble ratioComm      = commTime / fineSolveTime;
    NekDouble ratioOverhead  = overheadTime / fineSolveTime;

    return 1.0 / (ratioPredictor + ratio * (1.0 + ratioSolve + ratioInterp) +
                  (ratioComm * nComm + ratioOverhead) / m_numChunks);
}

/**
 *
 */
void DriverParareal::AssertParameters(void)
{
    // Assert time-stepping parameters
    ASSERTL0(
        m_fineSteps % m_numChunks == 0,
        "Total number of fine step should be divisible by number of chunks.");

    ASSERTL0(
        m_coarseSteps % m_numChunks == 0,
        "Total number of coarse step should be divisible by number of chunks.");

    ASSERTL0(m_fineSteps % (m_numChunks * m_numWindowsPIT) == 0,
             "Total number of fine step should be divisible by number of "
             "windows times number of chunks.");

    ASSERTL0(m_coarseSteps % (m_numChunks * m_numWindowsPIT) == 0,
             "Total number of coarse step should be divisible by number of "
             "windows times number of chunks.");

    ASSERTL0(fabs(m_coarseTimeStep * m_coarseSteps -
                  m_fineTimeStep * m_fineSteps) < 1e-12,
             "Fine and coarse total computational times do not match");

    // Assert I/O parameters
    if (m_infoSteps)
    {
        ASSERTL0(m_fineSteps % (m_infoSteps * m_numChunks * m_numWindowsPIT) ==
                     0,
                 "number of IO_InfoSteps should divide number of fine steps "
                 "per time chunk");
    }

    if (m_checkSteps)
    {
        ASSERTL0(m_fineSteps % (m_checkSteps * m_numChunks * m_numWindowsPIT) ==
                     0,
                 "number of IO_CheckSteps should divide number of fine steps "
                 "per time chunk");
    }
}

/**
 *
 */
void DriverParareal::RunCoarseSolve(
    const NekDouble time, const size_t nstep,
    const Array<OneD, const Array<OneD, NekDouble>> &input,
    Array<OneD, Array<OneD, NekDouble>> &output)
{
    // Set to coarse timestep.
    m_coarseEqSys->SetTime(time);
    m_coarseEqSys->SetSteps(nstep);

    // Restrict initial condition from input.
    Interpolator(input, m_tmpcoarse);

    // Copy initial condition.
    for (size_t i = 0; i < m_fineEqSys->GetNvariables(); ++i)
    {
        m_coarseEqSys->CopyToPhysField(i, m_tmpcoarse[i]);
    }

    // Solve equations.
    m_coarseEqSys->DoSolve();

    // Copy solution.
    for (size_t i = 0; i < m_coarseEqSys->GetNvariables(); ++i)
    {
        m_coarseEqSys->CopyFromPhysField(i, m_tmpcoarse[i]);
    }

    // Interpolate solution to output.
    Interpolator(m_tmpcoarse, output);
}

/**
 *
 */
void DriverParareal::RunFineSolve(
    const NekDouble time, const size_t nstep, const size_t iter,
    const size_t wd, const Array<OneD, const Array<OneD, NekDouble>> &input,
    Array<OneD, Array<OneD, NekDouble>> &output)
{
    // Number of checkpoint by chunk.
    size_t nChkPts = m_checkSteps ? m_fineSteps / m_checkSteps : 1;

    // Checkpoint index.
    size_t iChkPts = (m_chunkRank + wd * m_numChunks) * nChkPts + 1;

    // Set to fine timestep.
    m_fineEqSys->SetTime(time);
    m_fineEqSys->SetSteps(nstep);

    // Reinitialize check point number for each parallel-in-time iteration.
    m_fineEqSys->SetCheckpointNumber(iChkPts);

    // Update parallel-in-time iteration number.
    m_fineEqSys->SetIterationNumberPIT(iter + 1);

    // Update parallel-in-time window number.
    m_fineEqSys->SetWindowNumberPIT(wd);

    // Copy initial condition from input.
    for (size_t i = 0; i < m_fineEqSys->GetNvariables(); ++i)
    {
        m_fineEqSys->CopyToPhysField(i, input[i]);
    }

    // Solve equations.
    m_fineEqSys->DoSolve();

    // Copy solution to output.
    for (size_t i = 0; i < m_fineEqSys->GetNvariables(); ++i)
    {
        m_fineEqSys->CopyFromPhysField(i, output[i]);
    }
}

/**
 *
 */
void DriverParareal::PararealCorrection(
    const Array<OneD, const Array<OneD, NekDouble>> &coarse_new,
    const Array<OneD, const Array<OneD, NekDouble>> &coarse_old,
    Array<OneD, Array<OneD, NekDouble>> &fine)
{
    for (size_t i = 0; i < fine.size(); ++i)
    {
        Vmath::Vadd(fine[i].size(), fine[i], 1, coarse_new[i], 1, fine[i], 1);
        Vmath::Vsub(fine[i].size(), fine[i], 1, coarse_old[i], 1, fine[i], 1);
    }
}

/**
 *
 */
void DriverParareal::PrintSolutionFile(void)
{
    PrintHeaderTitle2("PRINT SOLUTION FILES");
    m_fineEqSys->Output();
}

/**
 *
 */
void DriverParareal::ApplyWindowing(
    const Array<OneD, const Array<OneD, NekDouble>> &in,
    Array<OneD, Array<OneD, NekDouble>> &out)
{
    LibUtilities::CommSharedPtr tComm = m_session->GetComm()->GetTimeComm();

    // Use last chunk solution as initial condition for the next
    // window.
    if (m_chunkRank == m_numChunks - 1)
    {
        for (size_t i = 0; i < in.size(); ++i)
        {
            Vmath::Vcopy(in[i].size(), in[i], 1, out[i], 1);
        }
    }

    // Broadcast I.C. for windowing.
    for (size_t i = 0; i < out.size(); ++i)
    {
        tComm->Bcast(out[i], m_numChunks - 1);
    }
}

/**
 *
 */
void DriverParareal::CopyConvergedCheckPoints(const size_t w, const size_t k,
                                              size_t kmax)
{
    // Determine max number of iteration.
    LibUtilities::CommSharedPtr tComm = m_session->GetComm()->GetTimeComm();
    tComm->AllReduce(kmax, Nektar::LibUtilities::ReduceMax);

    if (m_comm->GetSpaceComm()->GetRank() == 0 && m_checkSteps)
    {
        for (size_t j = k; j < kmax; j++)
        {
            // Copy converged solution files from directory corresponding to
            // iteration j to the directory corresponding to iteration j + 1.

            // Input directory name.
            std::string indir = m_fineEqSys->GetSessionName() + "_" +
                                boost::lexical_cast<std::string>(j) + ".pit";

            /// Output directory name.
            std::string outdir = m_fineEqSys->GetSessionName() + "_" +
                                 boost::lexical_cast<std::string>(j + 1) +
                                 ".pit";

            // Number of checkpoint by chunk.
            size_t nChkPts = m_fineSteps / m_checkSteps;

            // Checkpoint index.
            size_t iChkPts = (m_chunkRank + w * m_numChunks) * nChkPts + 1;

            for (size_t i = 0; i < nChkPts; i++)
            {
                // Filename corresponding to checkpoint iChkPts.
                std::string filename =
                    m_fineEqSys->GetSessionName() + "_" +
                    boost::lexical_cast<std::string>(iChkPts) + ".chk";

                // Intput full file name.
                std::string infullname = indir + "/" + filename;

                // Output full file name.
                std::string outfullname = outdir + "/" + filename;

                // Remove output file if already existing.
                fs::remove_all(outfullname);

                // Copy converged solution files.
                fs::copy(infullname, outfullname);
            }
        }
    }
}

} // namespace SolverUtils
} // namespace Nektar
