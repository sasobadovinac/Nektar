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
// Description: Compressible Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/CompressibleSolver.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>

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
        m_idealGas = boost::iequals(eosType,"IdealGas");
    }
    
    void CompressibleSolver::v_Solve(
        const int                                         nDim,
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
              Array<OneD,       Array<OneD, NekDouble> > &flux)
    {
        if (m_pointSolve)
        {
            int expDim      = nDim;
            int nvariables  = Fwd.num_elements();
            
            NekDouble rhouf, rhovf;
            
            // Check if PDE-based SC is used
            if (expDim == 1)
            {
                for (int i = 0; i < Fwd[0].num_elements(); ++i)
                {
                    v_PointSolve(
                        Fwd [0][i], Fwd [1][i], 0.0,   0.0,   Fwd [2][i],
                        Bwd [0][i], Bwd [1][i], 0.0,   0.0,   Bwd [2][i],
                        flux[0][i], flux[1][i], rhouf, rhovf, flux[2][i]);
                }
            }
            else if (expDim == 2)
            {
                if (nvariables == expDim+2)
                {
                    for (int i = 0; i < Fwd[0].num_elements(); ++i)
                    {
                        v_PointSolve(
                            Fwd [0][i], Fwd [1][i], Fwd [2][i], 0.0,   Fwd [3][i],
                            Bwd [0][i], Bwd [1][i], Bwd [2][i], 0.0,   Bwd [3][i],
                            flux[0][i], flux[1][i], flux[2][i], rhovf, flux[3][i]);
                    }
                }
                
                if (nvariables > expDim+2)
                {
                    for (int i = 0; i < Fwd[0].num_elements(); ++i)
                    {
                        v_PointSolveVisc(
                            Fwd [0][i], Fwd [1][i], Fwd [2][i], 0.0, Fwd [3][i], Fwd [4][i],
                            Bwd [0][i], Bwd [1][i], Bwd [2][i], 0.0, Bwd [3][i], Bwd [4][i],
                            flux[0][i], flux[1][i], flux[2][i], rhovf, flux[3][i], flux[4][i]);
                    }
                }
                
            }
            else if (expDim == 3)
            {
                for (int i = 0; i < Fwd[0].num_elements(); ++i)
                {
                    v_PointSolve(
                        Fwd [0][i], Fwd [1][i], Fwd [2][i], Fwd [3][i], Fwd [4][i],
                        Bwd [0][i], Bwd [1][i], Bwd [2][i], Bwd [3][i], Bwd [4][i],
                        flux[0][i], flux[1][i], flux[2][i], flux[3][i], flux[4][i]);
                }
                if (nvariables > expDim+2)
                {
                    for (int i = 0; i < Fwd[0].num_elements(); ++i)
                    {
                        v_PointSolveVisc(
                            Fwd [0][i], Fwd [1][i], Fwd [2][i], Fwd [3][i], Fwd [4][i], Fwd [5][i],
                            Bwd [0][i], Bwd [1][i], Bwd [2][i], Bwd [3][i], Bwd [4][i], Bwd [5][i],
                            flux[0][i], flux[1][i], flux[2][i], flux[3][i], flux[4][i], flux[5][i]);
                    }
                }
            }
        }
        else
        {
            v_ArraySolve(Fwd, Bwd, flux);
        }
    }

    void CompressibleSolver::v_Solve(
        const int                                         nDim,
        const Array<OneD, const Array<OneD, NekDouble> > &FwdJ,
        const Array<OneD, const Array<OneD, NekDouble> > &BwdJ,
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd ,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd ,
              Array<OneD,       Array<OneD, NekDouble> > &flux)
    {
        size_t nvars3D      = 5;
        size_t nvars        = FwdJ.num_elements();
        size_t nvars3DM1    = nvars3D-1;
        size_t nvarsM1      = nvars-1;
        Array<OneD, NekDouble> locFwdJ {nvars3D , 0.0};
        Array<OneD, NekDouble> locBwdJ {nvars3D , 0.0};
        Array<OneD, NekDouble> locFwd  {nvars3D , 0.0};
        Array<OneD, NekDouble> locBwd  {nvars3D , 0.0};
        Array<OneD, NekDouble> locFlux {nvars3D , 0.0};
        
        if (m_pointSolve)
        {
            for (int i = 0; i < Fwd[0].num_elements(); ++i)
            {
                locFwdJ[nvars3DM1] = FwdJ[nvarsM1][i];
                locBwdJ[nvars3DM1] = BwdJ[nvarsM1][i];
                locFwd [nvars3DM1] = Fwd [nvarsM1][i];
                locBwd [nvars3DM1] = Bwd [nvarsM1][i];
                for (size_t n = 0; n < nvarsM1; ++n)
                {
                    locFwdJ[n] = FwdJ[n][i];
                    locBwdJ[n] = BwdJ[n][i];
                    locFwd [n] = Fwd [n][i];
                    locBwd [n] = Bwd [n][i];
                }
                v_PointSolve(locFwd, locBwd, locFlux, locFwdJ, locBwdJ);
                
                flux [nvarsM1][i] = locFlux [nvars3DM1];
                for (size_t n = 0; n < nvarsM1; ++n)
                {
                    flux [n][i] = locFlux [n];
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
            if( abs(fac) > NekConstants::kNekZeroTol)
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

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
    void CompressibleSolver::v_CalcFluxJacobian(
        const int                                         nDim,
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
        const Array<OneD, const Array<OneD, NekDouble> > &normals,
              DNekBlkMatSharedPtr                        &FJac,
              DNekBlkMatSharedPtr                        &BJac)
    {
        int expDim      = nDim;
        int nvariables  = Fwd.num_elements();
        int nnomals     = normals.num_elements();
        // int nvariables3D = nvariables+2;
        int nvariables3D = 5;

        if (nvariables > expDim+2)
        {
            ASSERTL0(false,"nvariables > expDim+2 case not coded")
        }

        DNekMatSharedPtr PointFJac3D = MemoryManager<DNekMat>
            ::AllocateSharedPtr(nvariables3D, nvariables3D);
        DNekMatSharedPtr PointBJac3D = MemoryManager<DNekMat>
            ::AllocateSharedPtr(nvariables3D, nvariables3D);
        
        Array<OneD, NekDouble> PointFwd(nvariables3D,0.0),PointBwd(nvariables3D,0.0);
        Array<OneD, NekDouble> PointNormal(3,0.0);


        Array<OneD, unsigned int> index(nvariables);

        index[nvariables-1] = 4;
        for(int i=0;i<nvariables-1;i++)
        {
            index[i] = i;
        }

        int nj=0;
        int nk=0;
        for (int i = 0; i < Fwd[0].num_elements(); ++i)
        {
            for(int j=0; j< nnomals; j++)
            {
                PointNormal[j] = normals [j][i];
            }

            for(int j=0; j< nvariables; j++)
            {
                nj = index[j];
                PointFwd[nj] = Fwd [j][i];
                PointBwd[nj] = Bwd [j][i];
                // std::cout << "PointFwd["<<i<<"]"<< PointFwd[i]<<std::endl;
            }
            
            v_PointFluxJacobian(PointFwd,PointBwd,PointNormal,PointFJac3D,PointBJac3D);
            
            DNekMatSharedPtr PointFJac = MemoryManager<DNekMat>
                ::AllocateSharedPtr(nvariables, nvariables);
            DNekMatSharedPtr PointBJac = MemoryManager<DNekMat>
                ::AllocateSharedPtr(nvariables, nvariables);

            for(int j=0; j< nvariables; j++)
            {
                nj = index[j];
                for(int k=0; k< nvariables; k++)
                {
                    nk = index[k];
                    (*PointFJac)(j,k) = (*PointFJac3D)(nj,nk); 
                    (*PointBJac)(j,k) = (*PointBJac3D)(nj,nk); 
                }
            }

            FJac->SetBlock(i, i, PointFJac);
            BJac->SetBlock(i, i, PointBJac);
        }
    }

    /// Currently duplacate in compressibleFlowSystem
    /// if fsw=+-1 calculate Steger-Warming flux vector splitting flux Jacobian
    /// if fsw=0   calculate the Jacobian of the exact Euler flux 
    /// efix is the numerical flux entropy fix parameter
    void CompressibleSolver::PointFluxJacobian_pn(
            const Array<OneD, NekDouble>    &Fwd,
            const Array<OneD, NekDouble>    &normals,
            DNekMatSharedPtr                &FJac,
            const NekDouble                 efix,   
            const NekDouble                 fsw)
    {
        int Nvar = 5;
        NekDouble ro,vx,vy,vz,ps,gama,ae ;
        NekDouble a,a2,h,h0,v2,vn,eps,eps2;
        NekDouble nx,ny,nz;
        NekDouble sn,osn,nxa,nya,nza,vna;
        NekDouble l1,l4,l5,al1,al4,al5,x1,x2,x3,y1;
        NekDouble c1,d1,c2,d2,c3,d3,c4,d4,c5,d5;
        NekDouble sml_ssf= 1.0E-12;

        NekDouble fExactorSplt = 2.0-abs(fsw); 

        NekDouble   rhoL  = Fwd[0];
        NekDouble   rhouL = Fwd[1];
        NekDouble   rhovL = Fwd[2];
        NekDouble   rhowL = Fwd[3];
        NekDouble   EL    = Fwd[4];

        Array<OneD, NekDouble> JacData;
        JacData = FJac->GetPtr();

        NekDouble oro;
        ro  = rhoL;
        oro = 1.0/ro;
        vx = rhouL*oro;
        vy = rhovL*oro;
        vz = rhowL*oro;

        // Internal energy (per unit mass)
        NekDouble eL =
                (EL - 0.5 * (rhouL * vx + rhovL * vy + rhowL * vz)) *oro;

        ps      = m_eos->GetPressure(rhoL, eL);
        gama    = m_eos->GetGamma();

        ae = gama - 1.0;
        v2 = vx*vx + vy*vy + vz*vz;
        a2 = gama*ps/ro;
        h = a2/ae;

        h0 = h + 0.5*v2;
        a = sqrt(a2);

        nx = normals[0];
        ny = normals[1];
        nz = normals[2];
        vn = nx*vx + ny*vy + nz*vz;
        sn = std::max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf);
        osn = 1.0/sn;

        nxa = nx * osn;
        nya = ny * osn;
        nza = nz * osn;
        vna = vn * osn;
        l1 = vn;
        l4 = vn + sn*a;
        l5 = vn - sn*a;

        eps = efix*sn;
        eps2 = eps*eps;

        al1 = sqrt(l1*l1 + eps2);
        al4 = sqrt(l4*l4 + eps2);
        al5 = sqrt(l5*l5 + eps2);

        l1 = 0.5*(fExactorSplt*l1 + fsw*al1);
        l4 = 0.5*(fExactorSplt*l4 + fsw*al4);
        l5 = 0.5*(fExactorSplt*l5 + fsw*al5);

        x1 = 0.5*(l4 + l5);
        x2 = 0.5*(l4 - l5);
        x3 = x1 - l1;
        y1 = 0.5*v2;
        c1 = ae*x3/a2;
        d1 = x2/a;

        int ncl0    = 0;
        int ncl1    = Nvar;
        int ncl2    = 2*Nvar;
        int ncl3    = 3*Nvar;
        int ncl4    = 4*Nvar;

        JacData[     ncl0] = c1*y1 - d1*vna + l1;
        JacData[     ncl1] = -c1*vx + d1*nxa;
        JacData[     ncl2] = -c1*vy + d1*nya;
        JacData[     ncl3] = -c1*vz + d1*nza;
        JacData[     ncl4] = c1;

        c2 = c1*vx + d1*nxa*ae;
        d2 = x3*nxa + d1*vx;
        JacData[ 1 + ncl0] = c2*y1 - d2*vna;
        JacData[ 1 + ncl1] = -c2*vx + d2*nxa + l1;
        JacData[ 1 + ncl2] = -c2*vy + d2*nya;
        JacData[ 1 + ncl3] = -c2*vz + d2*nza;
        JacData[ 1 + ncl4] = c2;

        c3 = c1*vy + d1*nya*ae;
        d3 = x3*nya + d1*vy;
        JacData[ 2 + ncl0] = c3*y1 - d3*vna;
        JacData[ 2 + ncl1] = -c3*vx + d3*nxa;
        JacData[ 2 + ncl2] = -c3*vy + d3*nya + l1;
        JacData[ 2 + ncl3] = -c3*vz + d3*nza;
        JacData[ 2 + ncl4] = c3;

        c4 = c1*vz + d1*nza*ae;
        d4 = x3*nza + d1*vz;
        JacData[ 3 + ncl0] = c4*y1 - d4*vna;
        JacData[ 3 + ncl1] = -c4*vx + d4*nxa;
        JacData[ 3 + ncl2] = -c4*vy + d4*nya;
        JacData[ 3 + ncl3] = -c4*vz + d4*nza + l1;
        JacData[ 3 + ncl4] = c4;

        c5 = c1*h0 + d1*vna*ae;
        d5 = x3*vna + d1*h0;
        JacData[ 4 + ncl0] = c5*y1 - d5*vna;
        JacData[ 4 + ncl1] = -c5*vx + d5*nxa;
        JacData[ 4 + ncl2] = -c5*vy + d5*nya;
        JacData[ 4 + ncl3] = -c5*vz + d5*nza;
        JacData[ 4 + ncl4] = c5 + l1;
    }

    void CompressibleSolver::EulerFlux(
        double rhoL, double rhouL, double rhovL, double rhowL, double EL,
        double &f0, double &f1, double &f2, double &f3, double &f4)
    {
        // NekDouble ro  = rhoL;
        // NekDouble oro = 1.0/ro;
        // NekDouble uL = rhouL*oro;
        // NekDouble vL = rhovL*oro;
        // NekDouble wL = rhowL*oro;

        // // Internal energy (per unit mass)
        // NekDouble eL =
        //         (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL)) * oro;


        // NekDouble pL = m_eos->GetPressure(rhoL, eL);

        NekDouble orhoL = 1.0/rhoL;
        NekDouble uL = rhouL * orhoL;
        NekDouble vL = rhovL * orhoL;
        NekDouble wL = rhowL * orhoL;

        static NekDouble gamma = m_params["gamma"]();
        // Left and right pressures
        NekDouble pL = (gamma - 1.0) *
            (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
        f0 = rhouL;
        f1 = rhouL*uL + pL;
        f2 = rhouL*vL;
        f3 = rhouL*wL;
        f4 = uL*(EL + pL); 
    }

    void CompressibleSolver::EulerFlux(
        double rhoJ, double rhouJ, double rhovJ, double rhowJ, double EJ,
        double rhoL, double rhouL, double rhovL, double rhowL, double EL,
        double &f0, double &f1, double &f2, double &f3, double &f4)
    {
        int nvars = 5;
        int nspce = 3;
        NekDouble efix = 0.0;
        NekDouble fsw = 0.0;
        
        Array<OneD, NekDouble> normals{nspce,0.0};
        normals[0] = 1.0;

        Array<OneD, NekDouble> VarJ{nvars};
        VarJ[0] = rhoJ;
        VarJ[1] = rhouJ;
        VarJ[2] = rhovJ;
        VarJ[3] = rhowJ;
        VarJ[4] = EJ;
        DNekMatSharedPtr pointJac = MemoryManager<DNekMat>::AllocateSharedPtr
                        (nvars,nvars,0.0);

        PointFluxJacobian_pn(VarJ,normals,pointJac,efix,fsw);

        NekVector<NekDouble> VarVect{nvars};
        VarVect[0] = rhoL;
        VarVect[1] = rhouL;
        VarVect[2] = rhovL;
        VarVect[3] = rhowL;
        VarVect[4] = EL;

        NekVector<NekDouble> fluxVect (nvars, VarJ, eWrapper);

        fluxVect = (*pointJac)*VarVect;

        f0 = fluxVect[0];
        f1 = fluxVect[1];
        f2 = fluxVect[2];
        f3 = fluxVect[3];
        f4 = fluxVect[4];
    }
#endif
}
