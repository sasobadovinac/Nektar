///////////////////////////////////////////////////////////////////////////////
//
// File: UpwindSolver.cpp
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
// Description: HLL Riemann solver for the Linear Shallow Water Equations.
//              Only valid for constant depth
//
///////////////////////////////////////////////////////////////////////////////

#include <GalerkinBoltzmannSolver/RiemannSolvers/UpwindSolver.h>

namespace Nektar
{
    std::string UpwindSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "BoltzmannUpwind",
            UpwindSolver::create,
            "Upwind Riemann solver");

    UpwindSolver::UpwindSolver(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : BoltzmannSolver(pSession)
    {

    }

    /**
     * @brief Upwind Riemann solver for the  Boltzmann Equations
     *
     * @param etaL   Free surface elevation left state.
     * @param etaR   Free surface elevation right state.  
     * @param uL     x-velocity  left state.  
     * @param uR     x-velocity  right state.  
     * @param vL     y-velocity  left state.  
     * @param vR     y-velocity  right state. 
     * @param dL     still water depth component left state.  
     * @param dR     still water depth component right state. 
     * @param uf     Computed Riemann flux for x-momentum component 
     * @param vf     Computed Riemann flux for y-momentum component 
     */
    void UpwindSolver::v_PointSolve(
         NekDouble  a0L, NekDouble  a1L, NekDouble a2L,
         NekDouble  a3L, NekDouble  d4L, NekDouble a5L,
         NekDouble  a0R, NekDouble  a1R, NekDouble a2R,
         NekDouble  a3R, NekDouble  a4R, NekDouble a5R,
         NekDouble &f0,  NekDouble &f1, NekDouble &f2,
         NekDouble &f3,  NekDouble &f4, NekDouble &f5)
    {        
        // 1) Set up rotation of variables

        // 2) Upwinding of point

        // 3) Evaluate flux in unrotated frame. 
    }
}
