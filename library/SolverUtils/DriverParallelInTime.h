///////////////////////////////////////////////////////////////////////////////
//
// File DriverParallelInTime.h
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
// Description: Driver class for the parallel-in-time solver
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DRIVERPARALLELINTIME_H
#define NEKTAR_SOLVERUTILS_DRIVERPARALLELINTIME_H

#include <SolverUtils/Driver.h>
#include <SolverUtils/UnsteadySystem.h>

namespace Nektar
{
namespace SolverUtils
{

/// Base class for the development of parallel-in-time solvers.
class DriverParallelInTime : public Driver
{
public:
protected:
    /// Constructor
    SOLVER_UTILS_EXPORT DriverParallelInTime(
        const LibUtilities::SessionReaderSharedPtr pSession,
        const SpatialDomains::MeshGraphSharedPtr pGraph);

    /// Destructor
    SOLVER_UTILS_EXPORT virtual ~DriverParallelInTime();

    /// Virtual function for initialisation implementation.
    SOLVER_UTILS_EXPORT virtual void v_InitObject(
        std::ostream &out = std::cout) override;

    /// Virtual function for solve implementation.
    SOLVER_UTILS_EXPORT virtual void v_Execute(
        std::ostream &out = std::cout) override;

    virtual NekDouble v_EstimateCommunicationTime(void);

    virtual NekDouble v_EstimateRestrictionTime(void);

    virtual NekDouble v_EstimateInterpolationTime(void);

    virtual NekDouble v_EstimateCoarseSolverTime(void);

    virtual NekDouble v_EstimateFineSolverTime(void);

    virtual NekDouble v_EstimatePredictorTime(void);

    virtual NekDouble v_EstimateOverheadTime(void);

    virtual NekDouble v_ComputeSpeedUp(
        const size_t iter, NekDouble fineSolveTime, NekDouble coarseSolveTime,
        NekDouble restTime, NekDouble interTime, NekDouble commTime,
        NekDouble predictorTime, NekDouble overheadTime);

    void SetParallelInTimeSessionFile(void);

    void GetParametersFromSession(void);

    void InitialiseEqSystem(bool turnoff_output);

    void AllocateMemory(void);

    void InitialiseInterpolationField(void);

    void PrintCoarseSolverInfo(std::ostream &out = std::cout);

    void PrintFineSolverInfo(std::ostream &out = std::cout);

    void PrintHeaderTitle1(const std::string &title);

    void PrintHeaderTitle2(const std::string &title);

    void PrintComputationalTime(const NekDouble time);

    void RecvInitialConditionFromPreviousProc(
        Array<OneD, Array<OneD, NekDouble>> &array, int &convergence);

    void RecvInitialConditionFromPreviousProc(
        Array<OneD, Array<OneD, NekDouble>> &array);

    void SendSolutionToNextProc(Array<OneD, Array<OneD, NekDouble>> &array,
                                int &convergence);

    void SendSolutionToNextProc(Array<OneD, Array<OneD, NekDouble>> &array);

    void CopySolutionVector(const Array<OneD, const Array<OneD, NekDouble>> &in,
                            Array<OneD, Array<OneD, NekDouble>> &out);

    void CopyFromFinePhysField(Array<OneD, Array<OneD, NekDouble>> &out);

    void CopyFromCoarsePhysField(Array<OneD, Array<OneD, NekDouble>> &out);

    void CopyToFinePhysField(
        const Array<OneD, const Array<OneD, NekDouble>> &in);

    void CopyToCoarsePhysField(
        const Array<OneD, const Array<OneD, NekDouble>> &in);

    void UpdateSolution(const Array<OneD, const Array<OneD, NekDouble>> &in);

    void EvaluateExactSolution(const NekDouble &time);

    void SolutionConvergenceMonitoring(const NekDouble &time);

    void SolutionConvergenceSummary(const NekDouble &time);

    void UpdateErrorNorm(const bool normalized);

    void PrintErrorNorm(const bool normalized);

    void Interpolator(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                      Array<OneD, Array<OneD, NekDouble>> &outarray);

    NekDouble vL2ErrorMax(void);

    void SpeedUpAnalysis();

    void PrintSpeedUp(NekDouble fineSolveTime, NekDouble coarseSolveTime,
                      NekDouble restTime, NekDouble interTime,
                      NekDouble commTime, NekDouble predictorOverheadTime,
                      NekDouble overheadTime);

    NekDouble EstimateCommunicationTime(
        Array<OneD, Array<OneD, NekDouble>> &buffer1,
        Array<OneD, Array<OneD, NekDouble>> &buffer2);

    /// ParallelInTime (coarse solver) session reader object
    LibUtilities::SessionReaderSharedPtr m_sessionCoarse;

    /// ParallelInTime (coarse solver) MeshGraph object
    SpatialDomains::MeshGraphSharedPtr m_graphCoarse;

    /// Timestep for fine solver.
    NekDouble m_fineTimeStep;

    /// Timestep for coarse solver.
    NekDouble m_coarseTimeStep;

    /// Total time integration interval.
    NekDouble m_totalTime;

    /// Time for chunks
    NekDouble m_chunkTime;

    /// Local time
    NekDouble m_time;

    /// Number of steps for info I/O
    size_t m_infoSteps = 0;

    /// Number of steps for checkpoint
    size_t m_checkSteps = 0;

    /// Number of steps for the fine solver
    size_t m_fineSteps = 1;

    /// Number of steps for the coarse solver
    size_t m_coarseSteps = 1;

    /// Number of time chunks
    size_t m_numChunks = 1;

    /// Rank in time
    size_t m_chunkRank = 0;

    /// Maximum number of parallel-in-time iteration
    size_t m_iterMaxPIT = 0;

    // Number of windows for parallel-in-time time iteration
    size_t m_numWindowsPIT = 1;

    /// Using exact solution to compute error norms
    bool m_exactSolution = 0;

    /// ParallelInTime tolerance
    NekDouble m_tolerPIT = 1e-15;

    // Temporary storage for parallel-in-time iteration
    size_t m_fineQuadPts;
    size_t m_coarseQuadPts;
    size_t m_fineNpts;
    size_t m_coarseNpts;
    size_t m_nVar;
    std::shared_ptr<SolverUtils::UnsteadySystem> m_fineEqSys;
    std::shared_ptr<SolverUtils::UnsteadySystem> m_coarseEqSys;
    Array<OneD, NekDouble> m_vL2Errors;
    Array<OneD, NekDouble> m_vLinfErrors;
    Array<OneD, Array<OneD, NekDouble>> m_exactsoln;
    Array<OneD, Array<OneD, NekDouble>> m_tmpfine;
    Array<OneD, Array<OneD, NekDouble>> m_tmpcoarse;
    Array<OneD, MultiRegions::ExpListSharedPtr> m_fineFields;
    Array<OneD, MultiRegions::ExpListSharedPtr> m_coarseFields;
};

/// Interpolate from an expansion to an expansion
void InterpExp1ToExp2(const Array<OneD, MultiRegions::ExpListSharedPtr> exp1,
                      Array<OneD, MultiRegions::ExpListSharedPtr> &exp2);

} // namespace SolverUtils
} // namespace Nektar

#endif // NEKTAR_SOLVERUTILS_DRIVERPARALLELINTIME_H
