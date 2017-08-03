///////////////////////////////////////////////////////////////////////////////
//
// File: HLLEMSolver.cpp
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
// Description: HLLEM Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/HLLEMSolver.h>

namespace Nektar
{
    std::string HLLEMSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "HLLEM",
            HLLEMSolver::create,
            "HLLEM Riemann solver");

    HLLEMSolver::HLLEMSolver() : CompressibleSolver()
    {

    }

    /**
     * @brief HLLEM Riemann solver
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
    void HLLEMSolver::v_PointSolve(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef,
        NekDouble dxL, NekDouble dxR)
    {        
        static NekDouble gamma = m_params["gamma"]();
        
        // Left and Right velocities
        NekDouble uL = rhouL / rhoL;
        NekDouble vL = rhovL / rhoL;
        NekDouble wL = rhowL / rhoL;
        NekDouble uR = rhouR / rhoR;
        NekDouble vR = rhovR / rhoR;
        NekDouble wR = rhowR / rhoR;

        // Left and right pressures
        NekDouble pL = (gamma - 1.0) *
            (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
        NekDouble pR = (gamma - 1.0) *
            (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));
        NekDouble cL = sqrt(gamma * pL / rhoL);
        NekDouble cR = sqrt(gamma * pR / rhoR);
        NekDouble hL = (EL + pL) / rhoL;
        NekDouble hR = (ER + pR) / rhoR;
        
        // Square root of rhoL and rhoR.
        NekDouble srL  = sqrt(rhoL);
        NekDouble srR  = sqrt(rhoR);
        NekDouble srLR = srL + srR;
        
        // Velocity Roe averages
        NekDouble uRoe   = (srL * uL + srR * uR) / srLR;
        NekDouble vRoe   = (srL * vL + srR * vR) / srLR;
        NekDouble wRoe   = (srL * wL + srR * wR) / srLR;
        NekDouble hRoe   = (srL * hL + srR * hR) / srLR;
        NekDouble URoe   = (uRoe * uRoe + vRoe * vRoe + wRoe * wRoe);
        NekDouble cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 * URoe));
        NekDouble cRoe2  = (gamma - 1.0)*(hRoe - 0.5 * URoe);
        
        // Maximum wave speeds
        NekDouble SL = std::min(uL-cL, uRoe-cRoe);
        NekDouble SR = std::max(uR+cR, uRoe+cRoe);
        
        // Compute right eigenvectors (Toro book, p. 108)
        NekDouble R[5][5] = {
            {1,                  1,        0,    0,    1},
            {uRoe - cRoe,        uRoe,     0,    0,    uRoe+cRoe},
            {vRoe,               vRoe,     1,    0,    vRoe},
            {wRoe,               wRoe,     0,    1,    wRoe},
            {hRoe - uRoe * cRoe, 0.5*URoe, vRoe, wRoe, hRoe + uRoe*cRoe },
        };
        
        // Compute left eigenvectors (Toro book, p. 108)
        NekDouble gamma1 = gamma-1.;
        NekDouble L[5][5] = {
            {hRoe + cRoe/gamma1 * (uRoe - cRoe), -uRoe - cRoe/gamma1, -vRoe, -wRoe, 1},
            {-2.*hRoe + 4./gamma1 * cRoe2, 2.*uRoe, 2.*vRoe, 2.*wRoe, -2.},
            {-2.*vRoe * cRoe2/gamma1, 0, 2.*cRoe2/gamma1, 0, 0},
            {-2.*wRoe * cRoe2/gamma1, 0, 0, 2.*cRoe2/gamma1, 0},
            {hRoe - cRoe/gamma1 * (uRoe + cRoe), -uRoe + cRoe/gamma1, -vRoe, -wRoe, 1},
        };
        
        NekDouble coeff = 0.5*gamma1/cRoe2;
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                L[i][j] *= coeff;
            }
        }
        
        // Eigenvalues
        NekDouble EigVal[5] = {uRoe - cRoe, uRoe, uRoe, uRoe, uRoe + cRoe};
        
        // Positive and negative eigenvalues
        NekDouble EigValm[5], EigValp[5];
        for (int i = 0; i < 5; ++i) {
            EigValm[i] = 0.5 * (EigVal[i] - std::abs(EigVal[i]));
            EigValp[i] = 0.5 * (EigVal[i] + std::abs(EigVal[i]));
        }
        
        // Antidiffusion coefficients
        NekDouble delta[5];
        for (int i = 0; i < 5; ++i) {
            delta[i] = 1.0 - EigValm[i]/SL - EigValp[i]/SR;
        }
        
        // Compute HLLEM matrix
        int kmin = 1;   // first intermediate field
        int kmax = 3;   // last intermediate field
        NekDouble coeffHLLEM[5][5];
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                coeffHLLEM[i][j] = 0.0;
                for (int k = kmin; k <= kmax; ++k) {
                    coeffHLLEM[i][j] += R[i][k]*delta[k]*L[k][j];
                }
            }
        }
        
        // Define the difference of conservative variables
        NekDouble dU[5];
        dU[0] = rhoR - rhoL;
        dU[1] = rhouR - rhouL;
        dU[2] = rhovR - rhovL;
        dU[3] = rhowR - rhowL;
        dU[4] = ER - EL;
        
        // Compute HLLEM contribution to the flux taking into account only intermediate fields
        NekDouble fHLLEM[5];
        for (int i = 0; i < 5; ++i) {
            fHLLEM[i] = 0.0;
            for (int k = kmin; k <= kmax; ++k) {
                fHLLEM[i] += coeffHLLEM[i][k]*dU[k];
            }
            fHLLEM[i] = SR*SL/(SR-SL) * fHLLEM[i];
        }

        // HLLEM Riemann fluxes (positive case)
        if (SL >= 0)
        {
            rhof  = rhouL;
            rhouf = rhouL * uL + pL;
            rhovf = rhouL * vL;
            rhowf = rhouL * wL;
            Ef    = uL * (EL + pL);
        }
        // HLLEM Riemann fluxes (negative case)
        else if (SR <= 0)
        {
            rhof  = rhouR;
            rhouf = rhouR * uR + pR;
            rhovf = rhouR * vR;
            rhowf = rhouR * wR;
            Ef    = uR * (ER + pR);
        }
        // HLLEM Riemann fluxes (general case (SL < 0 | SR > 0)
        else
        {
            // HLL part
            NekDouble tmp1 = 1.0 / (SR - SL);
            NekDouble tmp2 = SR * SL;
            NekDouble rhofHLL, rhoufHLL, rhovfHLL, rhowfHLL, EfHLL;
            rhofHLL  = (SR * rhouL - SL * rhouR + tmp2 * (rhoR - rhoL)) * tmp1;
            rhoufHLL = (SR * (rhouL * uL + pL) - SL * (rhouR * uR + pR) +
                     tmp2 * (rhouR - rhouL)) * tmp1;
            rhovfHLL = (SR * rhouL * vL - SL * rhouR * vR +
                     tmp2 * (rhovR - rhovL)) * tmp1;
            rhowfHLL = (SR * rhouL * wL - SL * rhouR * wR +
                     tmp2 * (rhowR - rhowL)) * tmp1;
            EfHLL    = (SR * uL * (EL + pL) - SL * uR * (ER + pR) +
                     tmp2 * (ER - EL)) * tmp1;
            
            // HLLEM part
            NekDouble phi = 0.5;  // flattener
            rhof  = rhofHLL  - phi*fHLLEM[0];
            rhouf = rhoufHLL - phi*fHLLEM[1];
            rhovf = rhovfHLL - phi*fHLLEM[2];
            rhowf = rhowfHLL - phi*fHLLEM[3];
            Ef    = EfHLL    - phi*fHLLEM[4];
            
        }
    }
}
