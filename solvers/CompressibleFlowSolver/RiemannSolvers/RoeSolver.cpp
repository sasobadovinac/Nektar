///////////////////////////////////////////////////////////////////////////////
//
// File: RoeSolver.cpp
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
// Description: Roe Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/RoeSolver.h>

namespace Nektar
{
    std::string RoeSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "Roe",
            RoeSolver::create,
            "Roe Riemann solver");

    RoeSolver::RoeSolver(const LibUtilities::SessionReaderSharedPtr& pSession)
        : CompressibleSolver(pSession)
    {

    }

    /**
     * @brief Roe Riemann solver.
     *
     * Stated equations numbers are from:
     *
     *   "Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical
     *   Introduction", E. F. Toro (3rd edition, 2009).
     *
     * We follow the algorithm prescribed following equation 11.70.
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
    void RoeSolver::v_PointSolve(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
    {        
        static NekDouble gamma = m_params["gamma"]();

        NekDouble  flux0 = 0.0;
        NekDouble  flux1 = 0.0;
        NekDouble  flux2 = 0.0;
        NekDouble  flux3 = 0.0;
        NekDouble  flux4 = 0.0;
        EulerFlux(rhoR, rhouR, rhovR, rhowR, ER,
                  flux0, flux1, flux2, flux3, flux4);

        EulerFlux(rhoL, rhouL, rhovL, rhowL, EL,
                  rhof, rhouf, rhovf, rhowf, Ef);

        rhof  = (rhof  + flux0) * 0.5;
        rhouf = (rhouf + flux1) * 0.5;
        rhovf = (rhovf + flux2) * 0.5;
        rhowf = (rhowf + flux3) * 0.5;
        Ef    = (Ef    + flux4) * 0.5;

        NekDouble uRoe;
        NekDouble vRoe;
        NekDouble wRoe;
        NekDouble hRoe;
        NekDouble URoe;
        NekDouble cRoe;
        NekDouble ocRoe;
        
        RoeAverage(rhoL, rhouL, rhovL, rhowL, EL,
                   rhoR, rhouR, rhovR, rhowR, ER,
                   uRoe, vRoe, wRoe, hRoe, URoe, cRoe, ocRoe);
        
        // Compute eigenvectors (equation 11.59).
        NekDouble k[5][5] = {
            {1, uRoe - cRoe, vRoe, wRoe, hRoe - uRoe * cRoe},
            {1, uRoe,        vRoe, wRoe, 0.5 * URoe},
            {0, 0,           1,    0,    vRoe},
            {0, 0,           0,    1,    wRoe},
            {1, uRoe+cRoe,  vRoe,  wRoe, hRoe + uRoe*cRoe}
        };

        // Calculate jumps \Delta u_i (defined preceding equation 11.67).
        NekDouble jump[5] = {
            rhoR  - rhoL,
            rhouR - rhouL,
            rhovR - rhovL,
            rhowR - rhowL,
            ER    - EL
        };

        // Define \Delta u_5 (equation 11.70).
        NekDouble jumpbar = jump[4] - (jump[2]-vRoe*jump[0])*vRoe -
            (jump[3]-wRoe*jump[0])*wRoe;
        
        // Compute wave amplitudes (equations 11.68, 11.69).
        NekDouble alpha[5];
        alpha[1] = (gamma-1.0)*(jump[0]*(hRoe - uRoe*uRoe) + uRoe*jump[1] -
                                jumpbar)*ocRoe*ocRoe;
        alpha[0] = (jump[0]*(uRoe + cRoe) - jump[1] - cRoe*alpha[1])*0.5*ocRoe;
        alpha[4] = jump[0] - (alpha[0] + alpha[1]);
        alpha[2] = jump[2] - vRoe * jump[0];
        alpha[3] = jump[3] - wRoe * jump[0];

        // Compute eigenvalues \lambda_i (equation 11.58).
        NekDouble uRoeAbs = fabs(uRoe);
        NekDouble lambda[5] = { 
            fabs(uRoe - cRoe),
            uRoeAbs,
            uRoeAbs,
            uRoeAbs,
            fabs(uRoe + cRoe)
        };

        // Finally perform summation (11.29).
        for (int i = 0; i < 5; ++i)
        {
            uRoeAbs = 0.5*alpha[i]*lambda[i];
			
            rhof  -= uRoeAbs*k[i][0];
            rhouf -= uRoeAbs*k[i][1];
            rhovf -= uRoeAbs*k[i][2];
            rhowf -= uRoeAbs*k[i][3];
            Ef    -= uRoeAbs*k[i][4];
        }
    }

    void RoeSolver::v_PointSolve(
        const Array<OneD, NekDouble> &Fwd, 
        const Array<OneD, NekDouble> &Bwd, 
        Array<OneD, NekDouble>       &flux, 
        const Array<OneD, NekDouble> &FwdJ, 
        const Array<OneD, NekDouble> &BwdJ)
    {        
        static NekDouble gamma = m_params["gamma"]();

        size_t nvars = flux.num_elements();

        Array<OneD, NekDouble>  ftmp {nvars , 0.0};

        EulerFlux(FwdJ[0], FwdJ[1], FwdJ[2], FwdJ[3], FwdJ[4],
                  Fwd[0],  Fwd[1],  Fwd[2],  Fwd[3],  Fwd[4],
                  ftmp[0], ftmp[1], ftmp[2], ftmp[3], ftmp[4]);

        EulerFlux(BwdJ[0], BwdJ[1], BwdJ[2], BwdJ[3], BwdJ[4],
                  Bwd[0],  Bwd[1],  Bwd[2],  Bwd[3],  Bwd[4],
                  flux[0], flux[1], flux[2], flux[3], flux[4]);
        for(size_t i = 0; i < nvars; ++i)
        {
            flux[i] = (flux[i] + ftmp[i])*0.5;
        }

        NekDouble uRoe;
        NekDouble vRoe;
        NekDouble wRoe;
        NekDouble hRoe;
        NekDouble URoe;
        NekDouble cRoe;
        NekDouble ocRoe;
        
        RoeAverage(FwdJ[0], FwdJ[1], FwdJ[2], FwdJ[3], FwdJ[4],
                   BwdJ[0], BwdJ[1], BwdJ[2], BwdJ[3], BwdJ[4],
                   uRoe, vRoe, wRoe, hRoe, URoe, cRoe, ocRoe);
        
        // Compute eigenvectors (equation 11.59).
        NekDouble k[5][5] = {
            {1, uRoe - cRoe, vRoe, wRoe, hRoe - uRoe * cRoe},
            {1, uRoe,        vRoe, wRoe, 0.5 * URoe},
            {0, 0,           1,    0,    vRoe},
            {0, 0,           0,    1,    wRoe},
            {1, uRoe+cRoe,  vRoe,  wRoe, hRoe + uRoe*cRoe}
        };

        // Calculate jumps \Delta u_i (defined preceding equation 11.67).
        NekDouble jump[5] = 
        {
            Bwd[0] - Fwd[0],
            Bwd[1] - Fwd[1],
            Bwd[2] - Fwd[2],
            Bwd[3] - Fwd[3],
            Bwd[4] - Fwd[4]      
        };

        // Define \Delta u_5 (equation 11.70).
        NekDouble jumpbar = jump[4] - (jump[2]-vRoe*jump[0])*vRoe -
            (jump[3]-wRoe*jump[0])*wRoe;
        
        // Compute wave amplitudes (equations 11.68, 11.69).
        NekDouble alpha[5];
        alpha[1] = (gamma-1.0)*(jump[0]*(hRoe - uRoe*uRoe) + uRoe*jump[1] -
                                jumpbar)*ocRoe*ocRoe;
        alpha[0] = (jump[0]*(uRoe + cRoe) - jump[1] - cRoe*alpha[1])*0.5*ocRoe;
        alpha[4] = jump[0] - (alpha[0] + alpha[1]);
        alpha[2] = jump[2] - vRoe * jump[0];
        alpha[3] = jump[3] - wRoe * jump[0];

        // Compute eigenvalues \lambda_i (equation 11.58).
        NekDouble uRoeAbs = fabs(uRoe);
        NekDouble lambda[5] = { 
            fabs(uRoe - cRoe),
            uRoeAbs,
            uRoeAbs,
            uRoeAbs,
            fabs(uRoe + cRoe)
        };

        // Finally perform summation (11.29).
        for (int i = 0; i < 5; ++i)
        {
            uRoeAbs = 0.5*alpha[i]*lambda[i];
			
            flux[0] -= uRoeAbs*k[i][0];
            flux[1] -= uRoeAbs*k[i][1];
            flux[2] -= uRoeAbs*k[i][2];
            flux[3] -= uRoeAbs*k[i][3];
            flux[4] -= uRoeAbs*k[i][4];
        }
    }

    void RoeSolver::RoeAverage(
        double rhoL, double rhouL, double rhovL, double rhowL, double EL,
        double rhoR, double rhouR, double rhovR, double rhowR, double ER,
        double &uRoe, double &vRoe, double &wRoe, double &hRoe, double &URoe,
        double &cRoe, double &ocRoe)
    {        
        static NekDouble gamma = m_params["gamma"]();
        
        // Left and right velocities
        NekDouble orhoL = 1.0/rhoL;
        NekDouble orhoR = 1.0/rhoR;
        NekDouble uL = rhouL * orhoL;
        NekDouble vL = rhovL * orhoL;
        NekDouble wL = rhowL * orhoL;
        NekDouble uR = rhouR * orhoR;
        NekDouble vR = rhovR * orhoR;
        NekDouble wR = rhowR * orhoR;

        // Left and right pressures
        NekDouble pL = (gamma - 1.0) *
            (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
        NekDouble pR = (gamma - 1.0) *
            (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));
        
        // Left and right enthalpy
        NekDouble hL = (EL + pL) * orhoL;
        NekDouble hR = (ER + pR) * orhoR;

        // Square root of rhoL and rhoR.
        NekDouble srL  = sqrt(rhoL);
        NekDouble srR  = sqrt(rhoR);
        NekDouble srLR = srL + srR;
        NekDouble osrLR = 1.0/srLR;
        
        // Velocity, enthalpy and sound speed Roe averages (equation 11.60).
        uRoe   = (srL * uL + srR * uR) * osrLR;
        vRoe   = (srL * vL + srR * vR) * osrLR;
        wRoe   = (srL * wL + srR * wR) * osrLR;
        hRoe   = (srL * hL + srR * hR) * osrLR;
        URoe   = (uRoe * uRoe + vRoe * vRoe + wRoe * wRoe);
        cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 * URoe));
        ocRoe  = 1.0/cRoe;
    }

}
