///////////////////////////////////////////////////////////////////////////////
//
// File: HLLCSolver.cpp
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
// Description: HLLC Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/HLLCSolver.h>

namespace Nektar
{
    std::string HLLCSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "HLLC", HLLCSolver::create, "HLLC Riemann solver");
    
    HLLCSolver::HLLCSolver(const LibUtilities::SessionReaderSharedPtr& pSession)
        : CompressibleSolver(pSession)
    {

    }
    
    /**
     * @brief HLLC Riemann solver
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
    void HLLCSolver::v_PointSolve(
        NekDouble  rhoL, NekDouble  rhouL, NekDouble  rhovL, NekDouble  rhowL, NekDouble  EL,
        NekDouble  rhoR, NekDouble  rhouR, NekDouble  rhovR, NekDouble  rhowR, NekDouble  ER,
        NekDouble &rhof, NekDouble &rhouf, NekDouble &rhovf, NekDouble &rhowf, NekDouble &Ef)
    {
        // Left and Right velocities
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
        // Speed of sound
        NekDouble cL = m_eos->GetSoundSpeed(rhoL, eL);
        NekDouble cR = m_eos->GetSoundSpeed(rhoR, eR);

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

        // Maximum wave speeds
        NekDouble SL = std::min(uL-cL, uRoe-cRoe);
        NekDouble SR = std::max(uR+cR, uRoe+cRoe);
        
        // HLLC Riemann fluxes (positive case)
        if (SL >= 0)
        {
            rhof  = rhouL;
            rhouf = rhouL * uL + pL;
            rhovf = rhouL * vL;
            rhowf = rhouL * wL;
            Ef    = uL * (EL + pL);
        }
        // HLLC Riemann fluxes (negative case)
        else if (SR <= 0)
        {
            rhof  = rhouR;
            rhouf = rhouR * uR + pR;
            rhovf = rhouR * vR;
            rhowf = rhouR * wR;
            Ef    = uR * (ER + pR);
        }
        // HLLC Riemann fluxes (general case (SL < 0 | SR > 0)
        else
        {
            NekDouble SM = (pR - pL + rhouL * (SL - uL) - rhouR * (SR - uR)) /
            (rhoL * (SL - uL) - rhoR * (SR - uR));
            NekDouble rhoML  = rhoL * (SL - uL) / (SL - SM);
            NekDouble rhouML = rhoML * SM;
            NekDouble rhovML = rhoML * vL;
            NekDouble rhowML = rhoML * wL;
            NekDouble EML    = rhoML * (EL / rhoL +
                                        (SM - uL) * (SM + pL / (rhoL * (SL - uL))));
            
            NekDouble rhoMR  = rhoR * (SR - uR) / (SR - SM);
            NekDouble rhouMR = rhoMR * SM;
            NekDouble rhovMR = rhoMR * vR;
            NekDouble rhowMR = rhoMR * wR;
            NekDouble EMR    = rhoMR * (ER / rhoR +
                                        (SM - uR) * (SM + pR / (rhoR * (SR - uR))));
            
            if (SL < 0.0 && SM >= 0.0)
            {
                rhof  = rhouL + SL * (rhoML - rhoL);
                rhouf = rhouL * uL + pL + SL * (rhouML - rhouL);
                rhovf = rhouL * vL + SL * (rhovML - rhovL);
                rhowf = rhouL * wL + SL * (rhowML - rhowL);
                Ef    = uL * (EL + pL) + SL * (EML - EL);
            }
            else if(SM < 0.0 && SR > 0.0)
            {
                rhof  = rhouR + SR * (rhoMR - rhoR);
                rhouf = rhouR * uR + pR + SR * (rhouMR - rhouR);
                rhovf = rhouR * vR + SR * (rhovMR - rhovR);
                rhowf = rhouR * wR + SR * (rhowMR - rhowR);
                Ef    = uR * (ER + pR) + SR * (EMR - ER);
            }
        }
    }

    void HLLCSolver::v_ArraySolve(
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
              Array<OneD,       Array<OneD, NekDouble> > &flux,
        const int nDim)
    {
        for (int i = 0; i < Fwd[0].size(); i++)
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

            // Left and Right velocities
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
            // Speed of sound
            NekDouble cL = m_eos->GetSoundSpeed(rhoL, eL);
            NekDouble cR = m_eos->GetSoundSpeed(rhoR, eR);

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
            
            // Maximum wave speeds
            NekDouble SL = std::min(uL-cL, uRoe-cRoe);
            NekDouble SR = std::max(uR+cR, uRoe+cRoe);

            // HLLC Riemann fluxes (positive case)
            if (SL >= 0)
            {
                for(int k=0; k<Fwd.size(); k++)
                {
                    flux[k][i] = Fwd[k][i]*uL;
                }

                flux[1][i] = Fwd[1][i]*uL+pL;
                flux[nDim+1][i] = (Fwd[nDim+1][i]+pL)*uL;
            }
            // HLLC Riemann fluxes (negative case)
            else if (SR <= 0)
            {
                for(int k=0; k<Bwd.size(); k++)
                {
                    flux[k][i] = Bwd[k][i]*uR;
                }

                flux[1][i] = Bwd[1][i]*uR+pR;
                flux[nDim+1][i] = (Bwd[nDim+1][i]+pL)*uR;
            }
            // HLLC Riemann fluxes (general case (SL < 0 | SR > 0)
            else
            {
                NekDouble SM = (pR - pL + rhouL * (SL - uL) - rhouR * (SR - uR)) /
                (rhoL * (SL - uL) - rhoR * (SR - uR));
                NekDouble rhoML  = rhoL * (SL - uL) / (SL - SM);
                NekDouble rhouML = rhoML * SM;
                NekDouble rhovML = rhoML * vL;
                NekDouble rhowML = rhoML * wL;
                NekDouble EML    = rhoML * (EL / rhoL +
                                            (SM - uL) * (SM + pL / (rhoL * (SL - uL))));
                
                NekDouble rhoMR  = rhoR * (SR - uR) / (SR - SM);
                NekDouble rhouMR = rhoMR * SM;
                NekDouble rhovMR = rhoMR * vR;
                NekDouble rhowMR = rhoMR * wR;
                NekDouble EMR    = rhoMR * (ER / rhoR +
                                            (SM - uR) * (SM + pR / (rhoR * (SR - uR))));

                NekDouble rhof = 0.0;
                NekDouble rhouf = 0.0;
                NekDouble rhovf = 0.0;
                NekDouble rhowf = 0.0;
                NekDouble Ef = 0.0;

                if (SL < 0.0 && SM >= 0.0)
                {
                    rhof  = rhouL + SL * (rhoML - rhoL);
                    rhouf = rhouL * uL + pL + SL * (rhouML - rhouL);
                    rhovf = rhouL * vL + SL * (rhovML - rhovL);
                    rhowf = rhouL * wL + SL * (rhowML - rhowL);
                    Ef    = uL * (EL + pL) + SL * (EML - EL);

                    for(int k=nDim+2; k<Fwd.size(); k++)
                    {
                        flux[k][i] = Fwd[k][i]*uL + Fwd[k][i]*SL*((SL - uL) / (SL - SM) - 1.0);
                    }
                }
                else if(SM < 0.0 && SR > 0.0)
                {
                    rhof  = rhouR + SR * (rhoMR - rhoR);
                    rhouf = rhouR * uR + pR + SR * (rhouMR - rhouR);
                    rhovf = rhouR * vR + SR * (rhovMR - rhovR);
                    rhowf = rhouR * wR + SR * (rhowMR - rhowR);
                    Ef    = uR * (ER + pR) + SR * (EMR - ER);

                    for(int k=nDim+2; k<Bwd.size(); k++)
                    {
                        flux[k][i] = Bwd[k][i]*uR + Bwd[k][i]*SR*((SL - uR) / (SL - SM) - 1.0);
                    }
                }

                if(nDim == 1)
                {
                    flux[0][i] = rhof; 
                    flux[1][i] = rhouf; 
                    flux[2][i] = Ef; 
                 }
                else if(nDim == 2)
                {
                    flux[0][i] = rhof; 
                    flux[1][i] = rhouf; 
                    flux[2][i] = rhovf; 
                    flux[3][i] = Ef; 
                }
                else if(nDim == 3)
                {
                    flux[0][i] = rhof; 
                    flux[1][i] = rhouf; 
                    flux[2][i] = rhovf; 
                    flux[3][i] = rhowf; 
                    flux[4][i] = Ef; 
                }
            }
        }
    }
}
