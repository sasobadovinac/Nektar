///////////////////////////////////////////////////////////////////////////////
//
// File: FORCESolver.cpp
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
// Description: FORCE Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/FORCESolver.h>

namespace Nektar
{
    std::string FORCESolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "FORCE", FORCESolver::create, "FORCE Riemann solver");
    
    FORCESolver::FORCESolver() : CompressibleSolver()
    {
    }
    
    /**
     * @brief FORCE Riemann solver
     *
     * @param rhoL      Density left state.
     * @param rhoR      Density right state.
     * @param rhouL     x-momentum component left state.
     * @param rhouR     x-momentum component right state.
     * @param rhovL     y-momentum component left state.
     * @param rhovR     y-momentum component right state.
     * @param rhowL     z-momentum component left state.
     * @param rhowR     z-momentum component right state.
     * @param EL        Energy left state.
     * @param ER        Energy right state.
     * @param rhof      Computed Riemann flux for density.
     * @param rhouf     Computed Riemann flux for x-momentum component
     * @param rhovf     Computed Riemann flux for y-momentum component
     * @param rhowf     Computed Riemann flux for z-momentum component
     * @param Ef        Computed Riemann flux for energy.
     */
    void FORCESolver::v_PointSolve(
        NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL,
        NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER,
        NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef,
        NekDouble dxL, NekDouble dxR)
    {
        static NekDouble dt      = m_params["dt"]();
        static NekDouble dxForce = m_params["dxForce"]();
        static NekDouble alpha   = m_params["alpha"]();
        static NekDouble gamma   = m_params["gamma"]();

        // Left and Right velocities
        NekDouble uL = rhouL / rhoL;
        NekDouble vL = rhovL / rhoL;
        NekDouble wL = rhowL / rhoL;
        NekDouble uR = rhouR / rhoR;
        NekDouble vR = rhovR / rhoR;
        NekDouble wR = rhowR / rhoR;
        
        // Left and Right pressure
        NekDouble pL = (gamma - 1.0) *
            (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
        NekDouble pR = (gamma - 1.0) *
            (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));
        
        // Left and Right fluxes
        NekDouble rhofL = rhouL;
        NekDouble rhoufL = rhouL * uL + pL;
        NekDouble rhovfL = rhouL * vL;
        NekDouble rhowfL = rhouL * wL;
        NekDouble EfL = uL * (EL + pL);
        NekDouble rhofR = rhouR;
        NekDouble rhoufR = rhouR * uR + pR;
        NekDouble rhovfR = rhouR * vR;
        NekDouble rhowfR = rhouR * wR;
        NekDouble EfR = uR * (ER + pR);

        // dx Force
        NekDouble dxF = std::min( dxL, dxR );
        // NekDouble dxF = dxForce;
        
        // std::cout << "dt     = " << dt    << std::endl;
        // std::cout << "dxL    = " << dxL   << std::endl;
        // std::cout << "dxR    = " << dxR   << std::endl;
        // std::cout << "alpha  = " << alpha << std::endl;
        //std::cout << "dxForce = " << dxForce << std::endl;

        
        // Lax-Wendroff alpha numerical cons vars
        NekDouble rhoLW = 0.5 * (rhoL + rhoR)
                        - 0.5 * alpha * dt / dxF * (rhofR - rhofL);
        NekDouble rhouLW = 0.5 * (rhouL + rhouR)
                         - 0.5 * alpha * dt / dxF * (rhoufR - rhoufL);
        NekDouble rhovLW = 0.5 * ( rhovL + rhovR )
                         - 0.5 * alpha * dt / dxF * (rhovfR - rhovfL);
        NekDouble rhowLW = 0.5 * ( rhowL + rhowR )
                         - 0.5 * alpha * dt / dxF * (rhowfR - rhowfL);
        NekDouble ELW = 0.5 * (EL + ER) - 0.5 * alpha * dt / dxF * (EfR - EfL);

        // Lax-Wendroff alpha velocities
        NekDouble uLW = rhouLW / rhoLW;
        NekDouble vLW = rhovLW / rhoLW;
        NekDouble wLW = rhowLW / rhoLW;
        
        // Lax-Wendroff alpha pressure
        NekDouble pLW = (gamma - 1.0) *
            (ELW - 0.5 * (rhouLW * uLW + rhovLW * vLW + rhowLW * wLW));

        // Lax-Wendroff alpha numerical fluxes
        NekDouble rhofLW = rhouLW;
        NekDouble rhoufLW = rhouLW * uLW + pLW;
        NekDouble rhovfLW = rhouLW * vLW;
        NekDouble rhowfLW = rhouLW * wLW;
        NekDouble EfLW = uLW * (ELW + pLW);

        // Lax-Friedrics alpha numerical fluxes
        NekDouble rhofLF = 0.5 * ( rhofL + rhofR )
            - 0.5 / alpha * dxF / dt * ( rhoR - rhoL );
        NekDouble rhoufLF = 0.5 * ( rhoufL + rhoufR )
            - 0.5 / alpha * dxF / dt * ( rhouR - rhouL );
        NekDouble rhovfLF = 0.5 * ( rhovfL + rhovfR )
            - 0.5 / alpha * dxF / dt * ( rhovR - rhovL );
        NekDouble rhowfLF = 0.5 * ( rhowfL + rhowfR )
            - 0.5 / alpha * dxF / dt * ( rhowR - rhowL );
        NekDouble EfLF = 0.5 * ( EfL + EfR )
            - 0.5 / alpha * dxF / dt * ( ER - EL );

        // FORCE numerical fluxes
        rhof  = 0.5 * ( rhofLW + rhofLF );
        rhouf = 0.5 * ( rhoufLW + rhoufLF );
        rhovf = 0.5 * ( rhovfLW + rhovfLF );
        rhowf = 0.5 * ( rhowfLW + rhowfLF );
        Ef    = 0.5 * ( EfLW + EfLF );
    }
}
