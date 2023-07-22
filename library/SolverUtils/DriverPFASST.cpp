///////////////////////////////////////////////////////////////////////////////
//
// File DriverPFASST.cpp
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
// Description: Driver class for the PFASST solver
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <SolverUtils/DriverPFASST.h>
#include <boost/format.hpp>

namespace Nektar
{
namespace SolverUtils
{
std::string DriverPFASST::className =
    GetDriverFactory().RegisterCreatorFunction("PFASST", DriverPFASST::create);
std::string DriverPFASST::driverLookupId =
    LibUtilities::SessionReader::RegisterEnumValue("Driver", "PFASST", 0);

/**
 *
 */
DriverPFASST::DriverPFASST(const LibUtilities::SessionReaderSharedPtr pSession,
                           const SpatialDomains::MeshGraphSharedPtr pGraph)
    : DriverParallelInTime(pSession, pGraph)
{
}

/**
 *
 */
DriverPFASST::~DriverPFASST()
{
}

/**
 *
 */
void DriverPFASST::v_InitObject(std::ostream &out)
{
    DriverParallelInTime::v_InitObject(out);
}

/**
 *
 */
void DriverPFASST::v_Execute(std::ostream &out)
{
    // Set timer.
    Nektar::LibUtilities::Timer timer;
    NekDouble CPUtime = 0.0;

    // Set time communication parameters.
    LibUtilities::CommSharedPtr tComm = m_session->GetComm()->GetTimeComm();
    m_numChunks                       = tComm->GetSize();
    m_chunkRank                       = tComm->GetRank();

    // Set parameters from session file.
    GetParametersFromSession();
    AssertParameters();

    // Print solver summay.
    PrintFineSolverInfo(out);
    PrintCoarseSolverInfo(out);

    // Initialization.
    InitialiseEqSystem(true);
    AllocateMemory();
    InitialiseSDCScheme(true);
    InitialiseInterpolationField();
    SetTimeInterpolator();

    // Start iteration windows.
    tComm->Block();
    size_t chkPts   = m_chunkRank;
    m_totalTime     = m_fineTimeStep * m_fineSteps;
    m_numWindowsPIT = m_fineSteps / m_numChunks;
    m_chunkTime     = m_fineTimeStep;
    for (size_t w = 0; w < m_numWindowsPIT; w++)
    {
        timer.Start();
        PrintHeaderTitle1((boost::format("WINDOWS #%1%") % (w + 1)).str());

        // Compute initial guess for coarse solver.
        m_time = (w * m_numChunks) * m_chunkTime;
        SetTime(m_time);
        CoarseResidualEval(0);
        CopyQuadratureSolutionAndResidual(m_coarseSDCSolver, 0);
        for (size_t k = 0; k < m_chunkRank; k++)
        {
            RunCoarseSweep();
            UpdateCoarseFirstQuadrature();
            CopyQuadratureSolutionAndResidual(m_coarseSDCSolver, 0);
            m_time += m_chunkTime;
            SetTime(m_time);
        }
        RunCoarseSweep();

        // Interpolate coarse solution and residual to fine.
        InterpolateCoarseSolution();
        InterpolateCoarseResidual();

        // Start PFASST iteration.
        size_t k            = 0;
        int convergenceCurr = 0;
        int convergencePrev = (m_chunkRank == 0);
        int convergence     = convergencePrev;
        while (k < m_iterMaxPIT && !convergenceCurr)
        {
            // The PFASST implementation follow "Bolten, M., Moser, D., & Speck,
            // R. (2017). A multigrid perspective on the parallel full
            // approximation scheme in space and time. Numerical Linear Algebra
            // with Applications, 24(6)".

            if (m_chunkRank == m_numChunks - 1 &&
                m_comm->GetSpaceComm()->GetRank() == 0)
            {
                std::cout << "Iteration " << k + 1 << std::endl << std::flush;
            }

            // Performe fine sweep (parallel-in-time).
            FineResidualEval(0);
            RunFineSweep();

            // Compute FAS correction.
            RestrictFineSolution();
            SaveCoarseSolution();
            RestrictFineResidual();
            SaveCoarseResidual();
            ComputeFASCorrection();
            EvaluateSDCResidualNorm();

            // Perform coarse sweep (serial-in-time).
            RecvInitialConditionFromPreviousProc(
                m_coarseSDCSolver->UpdateFirstQuadratureSolutionVector(),
                convergence);
            CoarseResidualEval(0);
            RunCoarseSweep();
            convergenceCurr = (vL2ErrorMax() < m_tolerPIT && convergence);
            SendSolutionToNextProc(
                m_coarseSDCSolver->UpdateLastQuadratureSolutionVector(),
                convergenceCurr);

            // Display L2norm.
            PrintErrorNorm(true);

            // Correct fine solution and residual.
            SendSolutionToNextProc(
                m_fineSDCSolver->UpdateLastQuadratureSolutionVector(),
                convergenceCurr);
            CorrectFineSolution();
            CorrectFineResidual();
            RecvInitialConditionFromPreviousProc(
                m_fineSDCSolver->UpdateFirstQuadratureSolutionVector(),
                convergencePrev);
            CorrectInitialFineSolution();

            k++;
        }

        // Apply windowing.
        if (w < m_numWindowsPIT - 1)
        {
            ApplyWindowing();
        }
        timer.Stop();

        // Update field and write output.
        WriteOutput(chkPts);
        chkPts += m_numChunks;
        CPUtime += timer.Elapsed().count();
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
NekDouble DriverPFASST::v_EstimateCommunicationTime(void)
{
    // Allocate memory.
    Array<OneD, Array<OneD, NekDouble>> buffer1(m_nVar);
    Array<OneD, Array<OneD, NekDouble>> buffer2(m_nVar);
    for (size_t i = 0; i < m_nVar; ++i)
    {
        buffer1[i] = Array<OneD, NekDouble>(m_coarseNpts, 0.0);
        buffer2[i] = Array<OneD, NekDouble>(m_coarseNpts, 0.0);
    }
    // Estimate coarse communication time.
    return EstimateCommunicationTime(buffer1, buffer2);
}

/**
 *
 */
NekDouble DriverPFASST::v_EstimateRestrictionTime(void)
{
    // Average restriction time over niter iteration.
    size_t niter = 20;
    Nektar::LibUtilities::Timer timer;
    timer.Start();
    for (size_t n = 0; n < niter; n++)
    {
        RestrictFineSolution();
        SaveCoarseSolution();
        RestrictFineResidual();
        SaveCoarseResidual();
        ComputeFASCorrection();
        EvaluateSDCResidualNorm();
    }
    timer.Stop();
    return timer.Elapsed().count() / niter;
}

/**
 *
 */
NekDouble DriverPFASST::v_EstimateInterpolationTime(void)
{
    // Average interpolation time over niter iteration.
    size_t niter = 20;
    Nektar::LibUtilities::Timer timer;
    timer.Start();
    for (size_t n = 0; n < niter; n++)
    {
        CorrectFineSolution();
        CorrectFineResidual();
        CorrectInitialFineSolution();
    }
    timer.Stop();
    return timer.Elapsed().count() / niter;
}

/**
 *
 */
NekDouble DriverPFASST::v_EstimateCoarseSolverTime(void)
{
    // Estimate coarse solver time.
    size_t niter = 20;
    Nektar::LibUtilities::Timer timer;
    timer.Start();
    for (int i = 0; i < niter; i++)
    {
        CoarseResidualEval(0);
        RunCoarseSweep();
    }
    timer.Stop();
    return timer.Elapsed().count() / niter;
}

/**
 *
 */
NekDouble DriverPFASST::v_EstimateFineSolverTime(void)
{
    // Estimate coarse solver time.
    size_t niter = 20;
    Nektar::LibUtilities::Timer timer;
    timer.Start();
    for (int i = 0; i < niter; i++)
    {
        FineResidualEval(0);
        RunFineSweep();
    }
    timer.Stop();
    return timer.Elapsed().count() / niter;
}

/**
 *
 */
NekDouble DriverPFASST::v_EstimatePredictorTime(void)
{
    // Estimate coarse overhead time.
    size_t niter = 20;
    Nektar::LibUtilities::Timer timer;
    timer.Start();
    for (int i = 0; i < niter; i++)
    {
        RunCoarseSweep();
        UpdateCoarseFirstQuadrature();
        CopyQuadratureSolutionAndResidual(m_coarseSDCSolver, 0);
    }
    timer.Stop();
    return timer.Elapsed().count() / niter;
}

/**
 *
 */
NekDouble DriverPFASST::v_EstimateOverheadTime(void)
{
    // Estimate overhead time.
    size_t niter = 20;
    Nektar::LibUtilities::Timer timer;
    timer.Start();
    for (int i = 0; i < niter; i++)
    {
        InterpolateCoarseSolution();
        InterpolateCoarseResidual();
        ApplyWindowing();
    }
    timer.Stop();
    return timer.Elapsed().count() / niter;
}

/**
 *
 */
NekDouble DriverPFASST::v_ComputeSpeedUp(
    const size_t iter, NekDouble fineSolveTime, NekDouble coarseSolveTime,
    NekDouble restTime, NekDouble interTime, NekDouble commTime,
    NekDouble predictorTime, NekDouble overheadTime)
{
    // The speed-up estimate based on "Emmett, M., & Minion, M. (2012). Toward
    // an efficient parallel in time method for partial differential equations.
    // Communications in Applied Mathematics and Computational Science, 7(1),
    // 105-132" and on "Lunet, T., Bodart, J., Gratton, S., &
    // Vasseur, X. (2018). Time-parallel simulation of the decay of homogeneous
    // turbulence using parareal with spatial coarsening. Computing and
    // Visualization in Science, 19, 31-44".

    size_t Kiter             = m_fineSDCSolver->GetMaxOrder();
    size_t nComm             = (iter * (2 * m_numChunks - iter - 1)) / 2;
    NekDouble ratio          = double(iter) / m_numChunks;
    NekDouble ratioPredictor = predictorTime / fineSolveTime;
    NekDouble ratioSolve     = coarseSolveTime / fineSolveTime;
    NekDouble ratioFAS       = (restTime + interTime) / fineSolveTime;
    NekDouble ratioComm      = commTime / fineSolveTime;
    NekDouble ratioOverhead  = overheadTime / fineSolveTime;

    // Speed-up relative to SDC.
    return Kiter / (ratioPredictor + ratio * (1.0 + ratioSolve + ratioFAS) +
                    (ratioComm * nComm + ratioOverhead) / m_numChunks);
    // Speed-up relative to MLSDC.
    /*return Kiter * (1.0 + ratioSolve + ratioFAS + ratioOverhead / m_numChunks)
       / (ratioPredictor + ratio * (1.0 + ratioSolve + ratioFAS) +
            (ratioComm * nComm + ratioOverhead) / m_numChunks);*/
}

/**
 *
 */
void DriverPFASST::AssertParameters(void)
{
    // Assert time-stepping parameters.
    ASSERTL0(
        m_fineSteps % m_numChunks == 0,
        "Total number of fine step should be divisible by number of chunks.");

    ASSERTL0(
        m_coarseSteps % m_numChunks == 0,
        "Total number of coarse step should be divisible by number of chunks.");

    ASSERTL0(m_coarseTimeStep == m_fineTimeStep,
             "Please use same timestep for both coarse and fine solver");

    ASSERTL0(m_coarseSteps == m_fineSteps,
             "Please use same timestep for both coarse and fine solver");

    // Assert I/O parameters.
    if (m_checkSteps)
    {
        ASSERTL0(m_fineSteps % m_checkSteps == 0,
                 "number of IO_CheckSteps should divide number of fine steps "
                 "per time chunk");
    }

    if (m_infoSteps)
    {
        ASSERTL0(m_fineSteps % m_infoSteps == 0,
                 "number of IO_InfoSteps should divide number of fine steps "
                 "per time chunk");
    }
}

/**
 *
 */
void DriverPFASST::InitialiseSDCScheme(bool pfasst)
{
    // Cast pointer for TimeIntegrationSchemeSDC.
    m_fineSDCSolver =
        std::dynamic_pointer_cast<LibUtilities::TimeIntegrationSchemeSDC>(
            m_fineEqSys->GetTimeIntegrationScheme());
    m_coarseSDCSolver =
        std::dynamic_pointer_cast<LibUtilities::TimeIntegrationSchemeSDC>(
            m_coarseEqSys->GetTimeIntegrationScheme());

    // Assert if a SDC time-integration is used.
    ASSERTL0(m_fineSDCSolver != nullptr,
             "Should only be run with a SDC method");

    ASSERTL0(m_coarseSDCSolver != nullptr,
             "Should only be run with a SDC method");

    // Initialize fine SDC scheme.
    CopyFromFinePhysField(m_tmpfine);
    m_fineSDCSolver->SetPFASST(pfasst);
    m_fineSDCSolver->InitializeScheme(
        m_fineTimeStep, m_tmpfine, 0.0,
        m_fineEqSys->GetTimeIntegrationSchemeOperators());

    // Initialize coarse SDC scheme.
    CopyFromCoarsePhysField(m_tmpcoarse);
    m_coarseSDCSolver->SetPFASST(pfasst);
    m_coarseSDCSolver->InitializeScheme(
        m_coarseTimeStep, m_tmpcoarse, 0.0,
        m_coarseEqSys->GetTimeIntegrationSchemeOperators());

    // Set some member variables.
    m_fineQuadPts   = m_fineSDCSolver->GetQuadPtsNumber();
    m_coarseQuadPts = m_coarseSDCSolver->GetQuadPtsNumber();

    // Alocate memory.
    m_solutionRest =
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_coarseQuadPts);
    m_residualRest =
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_coarseQuadPts);
    m_tmpfine_arr =
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_coarseQuadPts);
    m_tmpcoarse_arr =
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_coarseQuadPts);
    for (size_t n = 0; n < m_coarseQuadPts; ++n)
    {
        m_solutionRest[n]  = Array<OneD, Array<OneD, NekDouble>>(m_nVar);
        m_residualRest[n]  = Array<OneD, Array<OneD, NekDouble>>(m_nVar);
        m_tmpfine_arr[n]   = Array<OneD, Array<OneD, NekDouble>>(m_nVar);
        m_tmpcoarse_arr[n] = Array<OneD, Array<OneD, NekDouble>>(m_nVar);
        for (size_t i = 0; i < m_nVar; ++i)
        {
            m_solutionRest[n][i]  = Array<OneD, NekDouble>(m_coarseNpts, 0.0);
            m_residualRest[n][i]  = Array<OneD, NekDouble>(m_coarseNpts, 0.0);
            m_tmpfine_arr[n][i]   = Array<OneD, NekDouble>(m_fineNpts, 0.0);
            m_tmpcoarse_arr[n][i] = Array<OneD, NekDouble>(m_coarseNpts, 0.0);
        }
    }
}

/**
 *
 */
void DriverPFASST::SetTimeInterpolator(void)
{
    // Initialize time interpolator.
    LibUtilities::PointsKey fpoints = m_fineSDCSolver->GetPointsKey();
    LibUtilities::PointsKey cpoints = m_coarseSDCSolver->GetPointsKey();
    DNekMatSharedPtr ImatFtoC =
        LibUtilities::PointsManager()[fpoints]->GetI(cpoints);
    DNekMatSharedPtr ImatCtoF =
        LibUtilities::PointsManager()[cpoints]->GetI(fpoints);
    m_ImatFtoC = Array<OneD, NekDouble>(m_fineQuadPts * m_coarseQuadPts, 0.0);
    m_ImatCtoF = Array<OneD, NekDouble>(m_fineQuadPts * m_coarseQuadPts, 0.0);

    // Determine if Radau quadrature are used.
    int i0 = m_fineSDCSolver->HasFirstQuadrature() ? 0 : 1;
    int j0 = m_coarseSDCSolver->HasFirstQuadrature() ? 0 : 1;

    // Adapt fine to coarse time interpolator.
    for (size_t i = i0; i < m_fineQuadPts; ++i)
    {
        for (size_t j = j0; j < m_coarseQuadPts; ++j)
        {
            m_ImatFtoC[i * m_coarseQuadPts + j] =
                (ImatFtoC
                     ->GetPtr())[(i - i0) * (m_coarseQuadPts - j0) + (j - j0)];
        }
    }
    if (j0 == 1)
    {
        m_ImatFtoC[0] = 1.0;
    }

    // Adapt coarse to fine time interpolator.
    for (size_t j = j0; j < m_coarseQuadPts; ++j)
    {
        for (size_t i = i0; i < m_fineQuadPts; ++i)
        {
            m_ImatCtoF[j * m_fineQuadPts + i] =
                (ImatCtoF
                     ->GetPtr())[(j - j0) * (m_fineQuadPts - i0) + (i - i0)];
        }
    }
    if (i0 == 1)
    {
        m_ImatCtoF[0] = 1.0;
    }
}

/**
 *
 */
bool DriverPFASST::IsNotInitialCondition(const size_t n)
{
    return !(n == 0 && m_chunkRank == 0);
}

/**
 *
 */
void DriverPFASST::SetTime(const NekDouble &time)
{
    m_coarseEqSys->SetTime(time);
    m_coarseSDCSolver->SetTime(time);
    m_fineEqSys->SetTime(time);
    m_fineSDCSolver->SetTime(time);
}

/**
 *
 */
void DriverPFASST::CopyQuadratureSolutionAndResidual(
    std::shared_ptr<LibUtilities::TimeIntegrationSchemeSDC> SDCsolver,
    const int index)
{
    int nQuad = SDCsolver->GetQuadPtsNumber();
    int nVar  = SDCsolver->GetNvars();

    for (int n = 0; n < nQuad; ++n)
    {
        if (n != index)
        {
            for (int i = 0; i < nVar; ++i)
            {
                Vmath::Vcopy(SDCsolver->GetNpoints(),
                             SDCsolver->GetSolutionVector()[index][i], 1,
                             SDCsolver->UpdateSolutionVector()[n][i], 1);
                Vmath::Vcopy(SDCsolver->GetNpoints(),
                             SDCsolver->GetResidualVector()[index][i], 1,
                             SDCsolver->UpdateResidualVector()[n][i], 1);
            }
        }
    }
}

/**
 *
 */
void DriverPFASST::UpdateCoarseFirstQuadrature(void)
{
    m_coarseSDCSolver->UpdateFirstQuadrature();
}

/**
 *
 */
void DriverPFASST::UpdateFineFirstQuadrature(void)
{
    m_fineSDCSolver->UpdateFirstQuadrature();
}

/**
 *
 */
void DriverPFASST::ComputeCoarseInitialGuess(void)
{
    m_coarseSDCSolver->ComputeInitialGuess(m_chunkTime);
    m_coarseSDCSolver->UpdateLastQuadrature();
}

/**
 *
 */
void DriverPFASST::ComputeFineInitialGuess(void)
{
    m_fineSDCSolver->ComputeInitialGuess(m_chunkTime);
    m_fineSDCSolver->UpdateLastQuadrature();
}

/**
 *
 */
void DriverPFASST::RunCoarseSweep(void)
{
    // Start SDC iteration loop.
    for (size_t k = 0; k < m_coarseSDCSolver->GetOrder(); k++)
    {
        m_coarseSDCSolver->SDCIterationLoop(m_chunkTime);
    }

    // Update last quadrature point.
    m_coarseSDCSolver->UpdateLastQuadrature();
}

/**
 *
 */
void DriverPFASST::RunFineSweep(void)
{
    // Start SDC iteration loop.
    for (size_t k = 0; k < m_fineSDCSolver->GetOrder(); k++)
    {
        m_fineSDCSolver->SDCIterationLoop(m_chunkTime);
    }

    // Update last quadrature point.
    m_fineSDCSolver->UpdateLastQuadrature();
}

/**
 *
 */
void DriverPFASST::CoarseResidualEval(const size_t n)
{
    m_coarseSDCSolver->ResidualEval(m_chunkTime, n);
}

/**
 *
 */
void DriverPFASST::CoarseResidualEval(void)
{
    m_coarseSDCSolver->ResidualEval(m_chunkTime);
}

/**
 *
 */
void DriverPFASST::FineResidualEval(const size_t n)
{
    m_fineSDCSolver->ResidualEval(m_chunkTime, n);
}

/**
 *
 */
void DriverPFASST::FineResidualEval(void)
{
    m_fineSDCSolver->ResidualEval(m_chunkTime);
}

/**
 *
 */
void DriverPFASST::CoarseIntegratedResidualEval(void)
{
    m_coarseSDCSolver->UpdateIntegratedResidualQFint(m_chunkTime);
}

/**
 *
 */
void DriverPFASST::FineIntegratedResidualEval(void)
{
    m_fineSDCSolver->UpdateIntegratedResidualQFint(m_chunkTime);
}

/**
 *
 */
void DriverPFASST::Interpolate(
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &coarse,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fine, bool forced)
{

    // Interpolate coarse solution in space.
    for (size_t n = 0; n < m_coarseQuadPts; ++n)
    {
        Interpolator(coarse[n], m_tmpfine_arr[n]);
    }

    // Interpolate coarse solution in time.
    for (size_t n = 0; n < m_fineQuadPts; ++n)
    {
        if (forced || IsNotInitialCondition(n))
        {
            for (size_t i = 0; i < m_nVar; ++i)
            {
                Vmath::Zero(m_fineNpts, fine[n][i], 1);
                for (size_t k = 0; k < m_coarseQuadPts; ++k)
                {
                    size_t index = k * m_fineQuadPts + n;
                    Vmath::Svtvp(m_fineNpts, m_ImatCtoF[index],
                                 m_tmpfine_arr[k][i], 1, fine[n][i], 1,
                                 fine[n][i], 1);
                }
            }
        }
    }
}

/**
 *
 */
void DriverPFASST::InterpolateCoarseSolution(void)
{
    Interpolate(m_coarseSDCSolver->GetSolutionVector(),
                m_fineSDCSolver->UpdateSolutionVector(), false);
}

/**
 *
 */
void DriverPFASST::InterpolateCoarseResidual(void)
{
    Interpolate(m_coarseSDCSolver->GetResidualVector(),
                m_fineSDCSolver->UpdateResidualVector(), true);
}

/**
 *
 */
void DriverPFASST::Restrict(
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fine,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &coarse)
{
    // Restrict fine solution in time.
    for (size_t n = 0; n < m_coarseQuadPts; ++n)
    {
        for (size_t i = 0; i < m_nVar; ++i)
        {
            Vmath::Zero(m_fineNpts, m_tmpfine_arr[n][i], 1);
            for (size_t k = 0; k < m_fineQuadPts; ++k)
            {
                size_t index = k * m_coarseQuadPts + n;
                Vmath::Svtvp(m_fineNpts, m_ImatFtoC[index], fine[k][i], 1,
                             m_tmpfine_arr[n][i], 1, m_tmpfine_arr[n][i], 1);
            }
        }
    }

    // Restrict fine solution in space.
    for (size_t n = 0; n < m_coarseQuadPts; ++n)
    {
        Interpolator(m_tmpfine_arr[n], coarse[n]);
    }
}

/**
 *
 */
void DriverPFASST::RestrictFineSolution(void)
{
    Restrict(m_fineSDCSolver->GetSolutionVector(),
             m_coarseSDCSolver->UpdateSolutionVector());
}

/**
 *
 */
void DriverPFASST::RestrictFineResidual(void)
{
    CoarseResidualEval();
}

/**
 *
 */
void DriverPFASST::SaveCoarseSolution(void)
{
    for (size_t n = 0; n < m_coarseQuadPts; ++n)
    {
        CopySolutionVector(m_coarseSDCSolver->GetSolutionVector()[n],
                           m_solutionRest[n]);
    }
}

/**
 *
 */
void DriverPFASST::SaveCoarseResidual(void)
{
    if (!m_updateFineResidual)
    {
        for (size_t n = 0; n < m_coarseQuadPts; ++n)
        {
            CopySolutionVector(m_coarseSDCSolver->GetResidualVector()[n],
                               m_residualRest[n]);
        }
    }
}

/**
 *
 */
void DriverPFASST::ComputeFASCorrection()
{
    // Compute fine integrated residual.
    FineIntegratedResidualEval();

    // Compute coarse integrated residual.
    CoarseIntegratedResidualEval();

    // Restrict fine integratued residual.
    Restrict(m_fineSDCSolver->GetIntegratedResidualQFintVector(),
             m_coarseSDCSolver->UpdateFAScorrectionVector());

    // Compute coarse FAS correction terms.
    for (size_t n = 0; n < m_coarseQuadPts; ++n)
    {
        for (size_t i = 0; i < m_nVar; ++i)
        {
            Vmath::Vsub(
                m_coarseNpts, m_coarseSDCSolver->GetFAScorrectionVector()[n][i],
                1, m_coarseSDCSolver->GetIntegratedResidualQFintVector()[n][i],
                1, m_coarseSDCSolver->UpdateFAScorrectionVector()[n][i], 1);
        }
    }
}

/**
 *
 */
void DriverPFASST::Correct(const Array<OneD, Array<OneD, NekDouble>> &coarse,
                           Array<OneD, Array<OneD, NekDouble>> &fine,
                           bool forced)
{
    if (forced || IsNotInitialCondition(0))
    {
        // Compute difference between coarse solution and restricted
        // solution.
        Interpolator(fine, m_tmpcoarse);
        for (size_t i = 0; i < m_nVar; ++i)
        {
            Vmath::Vsub(m_coarseNpts, coarse[i], 1, m_tmpcoarse[i], 1,
                        m_tmpcoarse[i], 1);
        }

        // Add correction to fine solution.
        Interpolator(m_tmpcoarse, m_tmpfine);
        for (size_t i = 0; i < m_nVar; ++i)
        {
            Vmath::Vadd(m_fineNpts, m_tmpfine[i], 1, fine[i], 1, fine[i], 1);
        }
    }
}

/**
 *
 */
void DriverPFASST::CorrectInitialFineSolution(void)
{
    Correct(m_coarseSDCSolver->GetSolutionVector()[0],
            m_fineSDCSolver->UpdateSolutionVector()[0], false);
}

/**
 *
 */
void DriverPFASST::CorrectInitialFineResidual(void)
{
    Correct(m_coarseSDCSolver->GetResidualVector()[0],
            m_fineSDCSolver->UpdateResidualVector()[0], false);
}

/**
 *
 */
void DriverPFASST::Correct(
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &rest,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &coarse,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fine, bool forced)
{
    // Compute difference between coarse solution and restricted
    // solution.
    for (size_t n = 0; n < m_coarseQuadPts; ++n)
    {
        if (forced || IsNotInitialCondition(n))
        {
            for (size_t i = 0; i < m_nVar; ++i)
            {
                Vmath::Vsub(m_coarseNpts, coarse[n][i], 1, rest[n][i], 1,
                            m_tmpcoarse_arr[n][i], 1);
            }
        }
        else
        {
            for (size_t i = 0; i < m_nVar; ++i)
            {
                Vmath::Zero(m_coarseNpts, m_tmpcoarse_arr[n][i], 1);
            }
        }
    }

    // Interpolate coarse solution delta in space.
    for (size_t n = 0; n < m_coarseQuadPts; ++n)
    {
        Interpolator(m_tmpcoarse_arr[n], m_tmpfine_arr[n]);
    }

    // Interpolate coarse solution delta in time and correct fine solution.
    for (size_t n = 0; n < m_fineQuadPts; ++n)
    {
        if (forced || IsNotInitialCondition(n))
        {
            for (size_t i = 0; i < m_nVar; ++i)
            {
                for (size_t k = 0; k < m_coarseQuadPts; ++k)
                {
                    size_t index = k * m_fineQuadPts + n;
                    Vmath::Svtvp(m_fineNpts, m_ImatCtoF[index],
                                 m_tmpfine_arr[k][i], 1, fine[n][i], 1,
                                 fine[n][i], 1);
                }
            }
        }
    }
}

/**
 *
 */
void DriverPFASST::CorrectFineSolution(void)
{
    Correct(m_solutionRest, m_coarseSDCSolver->GetSolutionVector(),
            m_fineSDCSolver->UpdateSolutionVector(), false);
}

/**
 *
 */
void DriverPFASST::CorrectFineResidual(void)
{
    // Evaluate fine residual.
    if (m_updateFineResidual)
    {
        for (size_t n = 1; n < m_fineQuadPts; ++n)
        {
            FineResidualEval(n);
        }
    }
    // Correct fine residual.
    else
    {
        Correct(m_residualRest, m_coarseSDCSolver->GetResidualVector(),
                m_fineSDCSolver->UpdateResidualVector(), false);
    }
}

/**
 *
 */
void DriverPFASST::ApplyWindowing(void)
{
    LibUtilities::CommSharedPtr tComm = m_session->GetComm()->GetTimeComm();

    // Use last chunk solution as initial condition for the next window.
    if (m_chunkRank == m_numChunks - 1)
    {
        UpdateFineFirstQuadrature();
        Interpolator(m_fineSDCSolver->GetSolutionVector()[0],
                     m_coarseSDCSolver->UpdateSolutionVector()[0]);
    }

    // Broadcast I.C. for windowing.
    for (size_t i = 0; i < m_nVar; ++i)
    {
        tComm->Bcast(m_coarseSDCSolver->UpdateSolutionVector()[0][i],
                     m_numChunks - 1);
        tComm->Bcast(m_fineSDCSolver->UpdateSolutionVector()[0][i],
                     m_numChunks - 1);
    }
}

/**
 *
 */
void DriverPFASST::EvaluateSDCResidualNorm(void)
{
    // Update integrated residual
    FineIntegratedResidualEval();

    // Compute SDC residual norm
    for (size_t i = 0; i < m_nVar; ++i)
    {
        Vmath::Vadd(
            m_fineNpts, m_fineSDCSolver->GetSolutionVector()[0][i], 1,
            m_fineSDCSolver
                ->GetIntegratedResidualQFintVector()[m_fineQuadPts - 1][i],
            1, m_exactsoln[i], 1);
        m_fineEqSys->CopyToPhysField(
            i, m_fineSDCSolver->GetSolutionVector()[m_fineQuadPts - 1][i]);
        m_vL2Errors[i]   = m_fineEqSys->L2Error(i, m_exactsoln[i], 1);
        m_vLinfErrors[i] = m_fineEqSys->LinfError(i, m_exactsoln[i]);
    }
}

/**
 *
 */
void DriverPFASST::WriteOutput(size_t chkPts)
{
    static size_t IOChkStep   = m_checkSteps ? m_checkSteps : m_fineSteps;
    static std::string newdir = m_session->GetSessionName() + ".pit";

    if ((chkPts + 1) % IOChkStep == 0)
    {
        size_t index         = (chkPts + 1) / IOChkStep;
        std::string filename = newdir + "/" + m_session->GetSessionName();

        if (!fs::is_directory(newdir))
        {
            fs::create_directory(newdir);
        }

        UpdateSolution(m_fineSDCSolver->UpdateLastQuadratureSolutionVector());
        m_fineEqSys->WriteFld(filename + "_" +
                              boost::lexical_cast<std::string>(index) + ".fld");
    }
}

} // namespace SolverUtils
} // namespace Nektar
