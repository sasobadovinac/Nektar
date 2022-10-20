///////////////////////////////////////////////////////////////////////////////
//
// File DriverParareal.h
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

#ifndef NEKTAR_SOLVERUTILS_DRIVERPARAREAL_H
#define NEKTAR_SOLVERUTILS_DRIVERPARAREAL_H

#include <FieldUtils/Interpolator.h>
#include <SolverUtils/Driver.h>

namespace Nektar
{
namespace SolverUtils
{

/// Base class for the development of solvers.
class DriverParareal : public Driver
{
public:
    friend class MemoryManager<DriverParareal>;

    /// Creates an instance of this class
    static DriverSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        DriverSharedPtr p =
            MemoryManager<DriverParareal>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of the class
    static std::string className;

protected:
    // Interpolator
    FieldUtils::Interpolator<Array<OneD, MultiRegions::ExpListSharedPtr>>
        m_interp;

    /// Parareal (coarse solver) session reader object
    LibUtilities::SessionReaderSharedPtr m_sessionCoarse;

    /// Parareal (coarse solver) MeshGraph object
    SpatialDomains::MeshGraphSharedPtr m_graphCoarse;

    /// Timestep for fine solver.
    NekDouble m_fineTimeStep;

    /// Timestep for coarse solver.
    NekDouble m_coarseTimeStep;

    /// Total time integration interval.
    NekDouble m_totalTime;

    /// Time for chunks
    NekDouble m_chunkTime;

    /// Coarse solver time factor
    NekDouble m_coarseSolveFactor = 100.0;

    /// Number of steps for the fine solver
    int m_fineSteps = 1;

    /// Number of steps for the coarse solver
    int m_coarseSteps = 1;

    /// Number of time chunks
    int m_numChunks = 1;

    /// Rank in time
    int m_chunkRank = 0;

    /// Maximum number of parareal iteration
    int m_pararealIterMax = 0;

    /// Constructor
    SOLVER_UTILS_EXPORT DriverParareal(
        const LibUtilities::SessionReaderSharedPtr pSession,
        const SpatialDomains::MeshGraphSharedPtr pGraph);

    /// Destructor
    SOLVER_UTILS_EXPORT virtual ~DriverParareal();

    /// Second-stage initialisation
    SOLVER_UTILS_EXPORT virtual void v_InitObject(
        std::ostream &out = std::cout);

    /// Virtual function for solve implementation.
    SOLVER_UTILS_EXPORT virtual void v_Execute(std::ostream &out = std::cout);

    /// Set the Parareal (coarse solver) session file
    void SetPararealSessionFile(void);

    void RunCoarseSolve(const NekDouble time,
                        const Array<OneD, const Array<OneD, NekDouble>> &input,
                        Array<OneD, Array<OneD, NekDouble>> &output);

    void RunFineSolve(const NekDouble time,
                      const Array<OneD, const Array<OneD, NekDouble>> &input,
                      Array<OneD, Array<OneD, NekDouble>> &output);

    static std::string driverLookupId;
};

} // namespace SolverUtils
} // namespace Nektar

#endif // NEKTAR_SOLVERUTILS_DRIVERSTANDARD_H
