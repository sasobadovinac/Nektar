///////////////////////////////////////////////////////////////////////////////
//
// File: IsentropicVortex.h
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
// Description: Euler equations for isentropic vortex
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_ISENVORTEX_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_ISENVORTEX_H

#include <CompressibleFlowSolver/EquationSystems/EulerCFE.h>

namespace Nektar
{

class IsentropicVortex : public EulerCFE
{
public:
    friend class MemoryManager<IsentropicVortex>;

    /// Creates an instance of this class.
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<IsentropicVortex>::AllocateSharedPtr(pSession,
                                                               pGraph);
        p->InitObject();
        return p;
    }
    /// Name of class.
    static std::string className;

    virtual ~IsentropicVortex();

protected:
    IsentropicVortex(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const SpatialDomains::MeshGraphSharedPtr &pGraph);

    /// Print a summary of time stepping parameters.
    virtual void v_GenerateSummary(SolverUtils::SummaryList &s) override;

    virtual void v_SetInitialConditions(NekDouble initialtime      = 0.0,
                                        bool dumpInitialConditions = true,
                                        const int domain = 0) override;

    virtual void v_EvaluateExactSolution(unsigned int field,
                                         Array<OneD, NekDouble> &outfield,
                                         const NekDouble time = 0.0) override;

private:
    /// Isentropic Vortex Test Case.
    void EvaluateIsentropicVortex(const Array<OneD, NekDouble> &x,
                                  const Array<OneD, NekDouble> &y,
                                  const Array<OneD, NekDouble> &z,
                                  Array<OneD, Array<OneD, NekDouble>> &u,
                                  NekDouble time, const int o = 0);

    /// Maximum strength of the perturbation
    NekDouble m_beta;
    /// Velocity in x-direction
    NekDouble m_u0;
    /// Velocity in y-direction
    NekDouble m_v0;
    /// Origin in x-direction
    NekDouble m_x0;
    /// Origin in y-direction
    NekDouble m_y0;
};
} // namespace Nektar
#endif
