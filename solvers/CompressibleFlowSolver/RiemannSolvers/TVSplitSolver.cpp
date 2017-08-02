///////////////////////////////////////////////////////////////////////////////
//
// File: TVSplitSolver.cpp
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
// Description: Toro-Vazquez flux splitting method.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/TVSplitSolver.h>

namespace Nektar
{
    std::string TVSplitSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "TVSplit", TVSplitSolver::create, "TVSplit Riemann solver");
    
    TVSplitSolver::TVSplitSolver() : CompressibleSolver()
    {
    }
    
    /**
     * @brief TVSplit Riemann solver
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
    void TVSplitSolver::v_PointSolve(
        NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL,
        NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER,
        NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef,
        NekDouble dxL, NekDouble dxR)
    {
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

        // Left and Right sound speed squared
        NekDouble scL = gamma * pL / rhoL;
        NekDouble scR = gamma * pR / rhoR;
        
        // Constants for eigenvalues (FOR IDEAL EOS!)
        NekDouble AL = sqrt(uL * uL + 4. * scL);
        NekDouble AR = sqrt(uR * uR + 4. * scR);
        
        // Constants for Riemann invariants
        NekDouble CL = rhoL * (uL - AL);
        NekDouble CR = rhoR * (uR + AR);
        
        // Intermediate pressure and velocity
        // obtained from Riemann solver for pressure system
        NekDouble u = (CR * uR - CL * uL) / (CR - CL) - 2. * (pR - pL) / (CR - CL);
        NekDouble p = (CR * pL - CL * pR) / (CR - CL) + 0.5 * CR * CL * (uR - uL) / (CR - CL);
        
        // Intermediate density
        // obtained from Riemann solver for pressure system
        NekDouble rho;
        if (u >= 0.0) {
            rho = rhoL;
        }
        else {
            rho = rhoR;
        }
        
        // Pressure flux
        NekDouble rhoeps = p / (gamma-1.);
        NekDouble rhofp  = 0.0;
        NekDouble rhoufp = p;
        NekDouble rhovfp = 0.0;
        NekDouble rhowfp = 0.0;
        NekDouble Efp    = u * (rhoeps + p);
        
        // Left and Right kinetic fluxes
        NekDouble rhofLk  = rhoL;
        NekDouble rhoufLk = rhouL;
        NekDouble rhovfLk = rhovL;
        NekDouble rhowfLk = rhowL;
        NekDouble EfLk    = 0.5 * rhoL * (uL * uL + vL * vL + wL * wL);
        NekDouble rhofRk  = rhoR;
        NekDouble rhoufRk = rhouR;
        NekDouble rhovfRk = rhovR;
        NekDouble rhowfRk = rhowR;
        NekDouble EfRk    = 0.5 * rhoR * (uR * uR + vR * vR + wR * wR);
        
        // Advection flux
        NekDouble rhofa, rhoufa, rhovfa, rhowfa, Efa;
        if (u >= 0.0) {
            rhofa  = u*rhofLk;
            rhoufa = u*rhoufLk;
            rhovfa = u*rhovfLk;
            rhowfa = u*rhowfLk;
            Efa    = u*EfLk;
        }
        else {
            rhofa  = u*rhofRk;
            rhoufa = u*rhoufRk;
            rhovfa = u*rhovfRk;
            rhowfa = u*rhowfRk;
            Efa    = u*EfRk;
        }
       
        // TV numerical flux = pressure flux + advection flux
        rhof  = rhofp  + rhofa;
        rhouf = rhoufp + rhoufa;
        rhovf = rhovfp + rhovfa;
        rhowf = rhowfp + rhowfa;
        Ef    = Efp    + Efa;
    }
}
