///////////////////////////////////////////////////////////////////////////////
//
// File: LaxFriedrichsSolver.cpp
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
// Description: LaxFriedrichs Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/LaxFriedrichsSolver.h>

namespace Nektar
{
    std::string LaxFriedrichsSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "LaxFriedrichs",
			LaxFriedrichsSolver::create,
            "Lax-Friedrichs Riemann solver");

    LaxFriedrichsSolver::LaxFriedrichsSolver(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : CompressibleSolver(pSession)
    {
        
    }
    
    /**
     * @brief Lax-Friedrichs Riemann solver
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
    void LaxFriedrichsSolver::v_PointSolve(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
    {
        // Left and right velocities
        NekDouble uL = rhouL / rhoL;
        NekDouble vL = rhovL / rhoL;
        NekDouble wL = rhowL / rhoL;
        NekDouble uR = rhouR / rhoR;
        NekDouble vR = rhovR / rhoR;
        NekDouble wR = rhowR / rhoR;

        // Internal energy (per unit mass)
        NekDouble eL =
                (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL)) / rhoL;
        NekDouble eR =
                (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR)) / rhoR;
        // Pressure
        NekDouble pL = m_eos->GetPressure(rhoL, eL);
        NekDouble pR = m_eos->GetPressure(rhoR, eR);

        // Left and right total enthalpy
        NekDouble HL = (EL + pL) / rhoL;
        NekDouble HR = (ER + pR) / rhoR;

        // Square root of rhoL and rhoR.
        NekDouble srL  = sqrt(rhoL);
        NekDouble srR  = sqrt(rhoR);
        NekDouble srLR = srL + srR;
        
        // Roe average state
        NekDouble uRoe   = (srL * uL + srR * uR) / srLR;
        NekDouble vRoe   = (srL * vL + srR * vR) / srLR;
        NekDouble wRoe   = (srL * wL + srR * wR) / srLR;
        NekDouble URoe2  = uRoe*uRoe + vRoe*vRoe + wRoe*wRoe;
        NekDouble HRoe   = (srL * HL + srR * HR) / srLR;
        NekDouble cRoe   = GetRoeSoundSpeed(
                                rhoL, pL, eL, HL, srL,
                                rhoR, pR, eR, HR, srR,
                                HRoe, URoe2, srLR);

		// Maximum eigenvalue
		NekDouble URoe = fabs(uRoe) + cRoe;
		
		// Lax-Friedrichs flux formula
        rhof  = 0.5*(rhouL + rhouR - URoe*(rhoR - rhoL));
        rhouf = 0.5*(pL + rhouL*uL + pR + rhouR*uR - URoe*(rhouR - rhouL));
        rhovf = 0.5*(rhouL*vL + rhouR*vR - URoe*(rhovR - rhovL));
        rhowf = 0.5*(rhouL*wL + rhouR*wR - URoe*(rhowR - rhowL));
        Ef    = 0.5*(uL*(EL + pL) + uR*(ER + pR) - URoe*(ER - EL));
    }


    void LaxFriedrichsSolver::v_ArraySolve(
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
              Array<OneD,       Array<OneD, NekDouble> > &flux,
        const int nDim)
    {
        for (int i = 0; i < Fwd[0].num_elements(); i++)
        {
            NekDouble  rhoL = 0; 
            NekDouble  rhouL = 0; 
            NekDouble  rhovL = 0; 
            NekDouble  rhowL = 0; 
            NekDouble  EL = 0;

                    
            NekDouble  rhoR = 0; 
            NekDouble  rhouR = 0; 
            NekDouble  rhovR = 0; 
            NekDouble  rhowR = 0; 
            NekDouble  ER = 0;

            if(nDim == 1)
            {
                rhoL = Fwd[0][i]; 
                rhouL = Fwd[1][i]; 
                EL = Fwd[2][i];

                        
                rhoR = Bwd[0][i]; 
                rhouR = Bwd[1][i]; 
                ER = Bwd[2][i];
            }
            else if(nDim == 2)
            {
                rhoL = Fwd[0][i]; 
                rhouL = Fwd[1][i]; 
                rhovL = Fwd[2][i]; 
                EL = Fwd[3][i];

                        
                rhoR = Bwd[0][i]; 
                rhouR = Bwd[1][i]; 
                rhovR = Bwd[2][i]; 
                ER = Bwd[3][i];
            }
            else if(nDim == 3)
            {
                rhoL = Fwd[0][i]; 
                rhouL = Fwd[1][i]; 
                rhovL = Fwd[2][i]; 
                rhowL = Fwd[3][i]; 
                EL = Fwd[4][i];

                        
                rhoR = Bwd[0][i]; 
                rhouR = Bwd[1][i]; 
                rhovR = Bwd[2][i]; 
                rhowR = Bwd[3][i]; 
                ER = Bwd[4][i];
            }

            // Left and right velocities
            NekDouble uL = rhouL / rhoL;
            NekDouble vL = rhovL / rhoL;
            NekDouble wL = rhowL / rhoL;
            NekDouble uR = rhouR / rhoR;
            NekDouble vR = rhovR / rhoR;
            NekDouble wR = rhowR / rhoR;

            // Internal energy (per unit mass)
            NekDouble eL =
                    (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL)) / rhoL;
            NekDouble eR =
                    (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR)) / rhoR;
            // Pressure
            NekDouble pL = m_eos->GetPressure(rhoL, eL);
            NekDouble pR = m_eos->GetPressure(rhoR, eR);

            // Left and right total enthalpy
            NekDouble HL = (EL + pL) / rhoL;
            NekDouble HR = (ER + pR) / rhoR;

            // Square root of rhoL and rhoR.
            NekDouble srL  = sqrt(rhoL);
            NekDouble srR  = sqrt(rhoR);
            NekDouble srLR = srL + srR;
            
            // Roe average state
            NekDouble uRoe   = (srL * uL + srR * uR) / srLR;
            NekDouble vRoe   = (srL * vL + srR * vR) / srLR;
            NekDouble wRoe   = (srL * wL + srR * wR) / srLR;
            NekDouble URoe2  = uRoe*uRoe + vRoe*vRoe + wRoe*wRoe;
            NekDouble HRoe   = (srL * HL + srR * HR) / srLR;
            NekDouble cRoe   = GetRoeSoundSpeed(
                                    rhoL, pL, eL, HL, srL,
                                    rhoR, pR, eR, HR, srR,
                                    HRoe, URoe2, srLR);

		    // Maximum eigenvalue
		    NekDouble URoe = fabs(uRoe) + cRoe;

		    // Lax-Friedrichs flux formula
            for(int k=0; k<Fwd.num_elements(); k++)
            {
                flux[k][i] = 0.5*(Fwd[k][i]*uL + Bwd[k][i]*uR - URoe*(Bwd[k][i] - Fwd[k][i]));
            }

            flux[1][i] = 0.5*(pL + Fwd[1][i]*uL + pR + Bwd[1][i]*uR - URoe*(Bwd[1][i] - Fwd[1][i]));
            flux[nDim+1][i] = 0.5*((pL + Fwd[nDim+1][i])*uL + (pR + Bwd[nDim+1][i])*uR - URoe*(Bwd[nDim+1][i] - Fwd[nDim+1][i]));
        }
    }
}
