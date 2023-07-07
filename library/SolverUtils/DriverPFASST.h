///////////////////////////////////////////////////////////////////////////////
//
// File DriverPFASST.h
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

#ifndef NEKTAR_SOLVERUTILS_DRIVERPFASST_H
#define NEKTAR_SOLVERUTILS_DRIVERPFASST_H

#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeSDC.h>
#include <SolverUtils/DriverParallelInTime.h>
#include <SolverUtils/UnsteadySystem.h>

namespace Nektar
{
namespace SolverUtils
{

/// Base class for the development of solvers.
class DriverPFASST : public DriverParallelInTime
{
public:
    friend class MemoryManager<DriverPFASST>;

    /// Creates an instance of this class
    static DriverSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        DriverSharedPtr p =
            MemoryManager<DriverPFASST>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of the class
    static std::string className;

protected:
    /// Constructor
    SOLVER_UTILS_EXPORT DriverPFASST(
        const LibUtilities::SessionReaderSharedPtr pSession,
        const SpatialDomains::MeshGraphSharedPtr pGraph);

    /// Destructor
    SOLVER_UTILS_EXPORT virtual ~DriverPFASST();

    /// Virtual function for initialisation implementation.
    SOLVER_UTILS_EXPORT virtual void v_InitObject(
        std::ostream &out = std::cout) override;

    /// Virtual function for solve implementation.
    SOLVER_UTILS_EXPORT virtual void v_Execute(
        std::ostream &out = std::cout) override;

    virtual NekDouble v_EstimateCommunicationTime(void) override;

    virtual NekDouble v_EstimateRestrictionTime(void) override;

    virtual NekDouble v_EstimateInterpolationTime(void) override;

    virtual NekDouble v_EstimateCoarseSolverTime(void) override;

    virtual NekDouble v_EstimateFineSolverTime(void) override;

    virtual NekDouble v_EstimatePredictorTime(void) override;

    virtual NekDouble v_EstimateOverheadTime(void) override;

    virtual NekDouble v_ComputeSpeedUp(
        const size_t iter, NekDouble fineSolveTime, NekDouble coarseSolveTime,
        NekDouble restTime, NekDouble interTime, NekDouble commTime,
        NekDouble predictorOverheadTime, NekDouble overheadTime) override;

    static std::string driverLookupId;

private:
    void AssertParameters(void);

    void InitialiseSDCScheme(bool pfasst);

    void SetTimeInterpolator(void);

    bool IsNotInitialCondition(const size_t n);

    void SetTime(const NekDouble &time);

    void CopyQuadratureSolutionAndResidual(
        std::shared_ptr<LibUtilities::TimeIntegrationSchemeSDC> SDCsolver,
        const int index);

    void UpdateCoarseFirstQuadrature(void);

    void UpdateFineFirstQuadrature(void);

    void ComputeCoarseInitialGuess(void);

    void ComputeFineInitialGuess(void);

    void RunCoarseSweep(void);

    void RunFineSweep(void);

    void CoarseResidualEval(const size_t n);

    void CoarseResidualEval(void);

    void FineResidualEval(const size_t n);

    void FineResidualEval(void);

    void CoarseIntegratedResidualEval(void);

    void FineIntegratedResidualEval(void);

    void Interpolate(
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &coarse,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fine, bool forced);

    void InterpolateCoarseSolution(void);

    void InterpolateCoarseResidual(void);

    void Restrict(const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fine,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &coarse);

    void RestrictFineSolution(void);

    void RestrictFineResidual(void);

    void SaveCoarseSolution(void);

    void SaveCoarseResidual(void);

    void ComputeFASCorrection(void);

    void Correct(const Array<OneD, Array<OneD, NekDouble>> &coarse,
                 Array<OneD, Array<OneD, NekDouble>> &fine, bool forced);

    void CorrectInitialFineSolution(void);

    void CorrectInitialFineResidual(void);

    void Correct(const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &rest,
                 const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &coarse,
                 Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fine,
                 bool forced);

    void CorrectFineSolution(void);

    void CorrectFineResidual(void);

    void ApplyWindowing(void);

    void EvaluateSDCResidualNorm(void);

    void WriteOutput(size_t chkPts);

    // Storage of PFASST
    Array<OneD, NekDouble> m_ImatFtoC;
    Array<OneD, NekDouble> m_ImatCtoF;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_solutionRest;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_residualRest;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_tmpfine_arr;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_tmpcoarse_arr;
    std::shared_ptr<LibUtilities::TimeIntegrationSchemeSDC> m_fineSDCSolver;
    std::shared_ptr<LibUtilities::TimeIntegrationSchemeSDC> m_coarseSDCSolver;

    bool m_updateFineResidual = false;
};

} // namespace SolverUtils
} // namespace Nektar

#endif // NEKTAR_SOLVERUTILS_DRIVERPFASST_H
