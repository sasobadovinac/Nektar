///////////////////////////////////////////////////////////////////////////////
//
// File: CompressibleSolver.cpp
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
// Description: Compressible Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include "CompressibleSolver.h"
#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace Nektar
{
    CompressibleSolver::CompressibleSolver(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : RiemannSolver(pSession), m_pointSolve(true)
    {
        m_requiresRotation = true;

        // Create equation of state object
        std::string eosType;
        pSession->LoadSolverInfo("EquationOfState",
                                  eosType, "IdealGas");
        m_eos = GetEquationOfStateFactory()
                                .CreateInstance(eosType, pSession);
        // Check if using ideal gas
        m_idealGas = boost::iequals(eosType, "IdealGas");
    }

    CompressibleSolver::CompressibleSolver()
        : RiemannSolver(), m_idealGas(true), m_pointSolve(true)
        {
            m_requiresRotation = true;
        }

    void CompressibleSolver::v_Solve(
        const int                                         nDim,
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
              Array<OneD,       Array<OneD, NekDouble> > &flux)
    {
        auto gridVelTrace = m_vectors["vgt"]();

        // Create fwd/bwd variables for grid velocity to be applied to
        Array<OneD, Array<OneD, NekDouble>> fwdTmp;
        CopyArray(Fwd, fwdTmp);

        Array<OneD, Array<OneD, NekDouble>> bwdTmp;
        CopyArray(Bwd, bwdTmp);


        if (m_pointSolve)
        {
            int expDim      = nDim;
            int nvariables  = Fwd.size();

            NekDouble rhouf{}, rhovf{};

            // @TODO: Idea here is to multiply the grid velocity by rhoL/rhoR and
            //  @TODO: then subtract from rho(u/v/w)/(L/R), instead of doing it in the HLLCsolver
            // Loop skips first entry which is rho
            for (int i = 1; i < gridVelTrace.size() + 1; ++i)
            {
                for (int j = 0; j < Fwd[0].size(); ++j)
                {
                    fwdTmp[i][j] -= Fwd[0][j] * gridVelTrace[i - 1][j];
                    bwdTmp[i][j] -= Bwd[0][j] * gridVelTrace[i - 1][j];
                }
            }

            if (expDim == 1)
            {
                for (int i = 0; i < Fwd[0].size(); ++i)
                {
                    v_PointSolve(
                        fwdTmp [0][i], fwdTmp [1][i], 0.0,   0.0,   fwdTmp [2][i],
                        bwdTmp [0][i], bwdTmp [1][i], 0.0,   0.0,   bwdTmp [2][i],
                        flux[0][i], flux[1][i], rhouf, rhovf, flux[2][i]);
                }
            }
            else if (expDim == 2)
            {
                for (int i = 0; i < Fwd[0].size(); ++i)
                {
                    v_PointSolve(
                        fwdTmp [0][i], fwdTmp [1][i], fwdTmp [2][i], 0.0,   fwdTmp [3][i],
                        bwdTmp [0][i], bwdTmp [1][i], bwdTmp [2][i], 0.0,   bwdTmp [3][i],
                        flux[0][i], flux[1][i], flux[2][i], rhovf, flux[3][i]);
                }
            }
            else if (expDim == 3)
            {
                for (int i = 0; i < Fwd[0].size(); ++i)
                {
                    v_PointSolve(
                        fwdTmp [0][i], fwdTmp [1][i], fwdTmp [2][i], fwdTmp [3][i], fwdTmp [4][i],
                        bwdTmp [0][i], bwdTmp [1][i], bwdTmp [2][i], bwdTmp [3][i], bwdTmp [4][i],
                        flux[0][i], flux[1][i], flux[2][i], flux[3][i], flux[4][i]);
                }
            }
        }
        else
        {
            v_ArraySolve(Fwd, Bwd, flux);
        }
    }

    NekDouble CompressibleSolver::GetRoeSoundSpeed(
        NekDouble rhoL, NekDouble pL, NekDouble eL, NekDouble HL, NekDouble srL,
        NekDouble rhoR, NekDouble pR, NekDouble eR, NekDouble HR, NekDouble srR,
        NekDouble HRoe, NekDouble URoe2, NekDouble srLR)
    {
        boost::ignore_unused(HL, srL, HR, srR, srLR);

        static NekDouble gamma = m_params["gamma"]();
        NekDouble cRoe;
        if(m_idealGas)
        {
            cRoe = sqrt((gamma - 1.0)*(HRoe - 0.5 * URoe2));
        }
        else
        {
            // Calculate static enthalpy of left and right states
            NekDouble hL = eL + pL/rhoL;
            NekDouble hR = eR + pR/rhoR;

            // Get partial derivatives of P(rho,e)
            NekDouble dpdeL   = m_eos->GetDPDe_rho(rhoL,eL);
            NekDouble dpdeR   = m_eos->GetDPDe_rho(rhoR,eR);
            NekDouble dpdrhoL = m_eos->GetDPDrho_e(rhoL,eL);
            NekDouble dpdrhoR = m_eos->GetDPDrho_e(rhoR,eR);

            // Evaluate chi and kappa parameters
            NekDouble chiL    = dpdrhoL - eL / rhoL * dpdeL;
            NekDouble kappaL  = dpdeL / rhoL;
            NekDouble chiR    = dpdrhoR - eR / rhoR * dpdeR;
            NekDouble kappaR  = dpdeR / rhoR;

            //
            // Calculate interface speed of sound using procedure from
            //    Vinokur, M.; Montagné, J.-L. "Generalized Flux-Vector
            //    Splitting and Roe Average for an Equilibrium Real Gas",
            //    JCP (1990).
            //

            // Calculate averages
            NekDouble avgChi    = 0.5 * (chiL      + chiR);
            NekDouble avgKappa  = 0.5 * (kappaL    + kappaR);
            NekDouble avgKappaH = 0.5 * (kappaL*hL + kappaR*hR);

            // Calculate jumps
            NekDouble deltaP    = pR      - pL;
            NekDouble deltaRho  = rhoR    - rhoL;
            NekDouble deltaRhoe = rhoR*eR - rhoL*eL;

            // Evaluate dP: equation (64) from Vinokur-Montagné
            NekDouble dP = deltaP - avgChi * deltaRho - avgKappa * deltaRhoe;
            // s (eq 66)
            NekDouble s  = avgChi + avgKappaH;
            // D (eq 65)
            NekDouble D  = (s*deltaRho)*(s*deltaRho) + deltaP*deltaP;
            // chiRoe and kappaRoe (eq 66)
            NekDouble chiRoe, kappaRoe;
            NekDouble fac = D - deltaP*deltaRho;
            if( std::abs(fac) > NekConstants::kNekZeroTol)
            {
                chiRoe   = (D*avgChi + s*s*deltaRho*dP) / fac;
                kappaRoe = D*avgKappa / fac;
            }
            else
            {
                chiRoe = avgChi;
                kappaRoe = avgKappa;
            }
            // Speed of sound (eq 53)
            cRoe = sqrt( chiRoe + kappaRoe*(HRoe - 0.5 * URoe2));
        }
        return cRoe;
    }

}
