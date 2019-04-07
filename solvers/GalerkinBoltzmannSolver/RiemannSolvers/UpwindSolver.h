///////////////////////////////////////////////////////////////////////////////
//
// File: UpwindSolver.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Upwind Riemann solver for the Galkerin Boltzmann Equations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SHALLOWWATERSOLVER_RIEMANNSOLVER_UPWINDSOLVER
#define NEKTAR_SOLVERS_SHALLOWWATERSOLVER_RIEMANNSOLVER_UPWINDSOLVER

#include <GalerkinBoltzmannSolver/RiemannSolvers/BoltzmannSolver.h>

namespace Nektar
{
    class UpwindSolver : public BoltzmannSolver
    {
    public:
        static RiemannSolverSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        {
            return RiemannSolverSharedPtr(
                new UpwindSolver(pSession));
        }
        
        static std::string solverName;
        
    protected:
        UpwindSolver(const LibUtilities::SessionReaderSharedPtr& pSession);
        
        virtual void v_PointSolve(
                   NekDouble  a0L, NekDouble  a1L, NekDouble a2L,
                   NekDouble  a3L, NekDouble  d4L, NekDouble a5L,
                   NekDouble  a0R, NekDouble  a1R, NekDouble a2R,
                   NekDouble  a3R, NekDouble  a4R, NekDouble a5R,
                   NekDouble &f0,  NekDouble &f1, NekDouble &f2,
                   NekDouble &f3,  NekDouble &f4, NekDouble &f5);
    };
}

#endif
