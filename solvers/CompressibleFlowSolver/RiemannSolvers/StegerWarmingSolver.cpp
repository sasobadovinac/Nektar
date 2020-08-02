///////////////////////////////////////////////////////////////////////////////
//
// File: StegerWarmingSolver.cpp
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
// Description: StegerWarming Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/StegerWarmingSolver.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
//#include <>

namespace Nektar
{
    std::string StegerWarmingSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "StegerWarming",
			StegerWarmingSolver::create,
            "StegerWarming Riemann solver");

    StegerWarmingSolver::StegerWarmingSolver(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : CompressibleSolver(pSession)
    {
    }

    /**
     * @brief StegerWarming Riemann solver
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
    void StegerWarmingSolver::v_PointSolve(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef)
    {
        // NekDouble nx,ny,nz;
        // NekDouble f1,f2,f3,f4,f5;
        // nx = 1.0;
        // ny = 0.0;
        // nz = 0.0;

        // NekDouble fsw,efix_StegerWarming;
        // efix_StegerWarming = 0.0;
        // fsw =  1.0;
        // flux_sw_pn(
        //     rhoL,   rhouL,  rhovL,  rhowL,  EL,
        //     nx  ,   ny   ,  nz,
        //     f1  ,   f2   ,  f3   ,  f4   ,  f5,
        //     efix_StegerWarming,   fsw);

        // rhof    = f1;
        // rhouf   = f2;
        // rhovf   = f3;
        // rhowf   = f4;
        // Ef      = f5;

        // fsw =-1.0;
        // flux_sw_pn(
        //     rhoR,   rhouR,  rhovR,  rhowR,  ER,
        //     nx  ,   ny   ,  nz,
        //     f1  ,   f2   ,  f3   ,  f4   ,  f5,
        //     efix_StegerWarming,   fsw);

        // rhof    += f1;
        // rhouf   += f2;
        // rhovf   += f3;
        // rhowf   += f4;
        // Ef      += f5;

        // if(false)
        // {
        //     // int nvariables=5;
        //     // Array<OneD, NekDouble> PointFwd(nvariables,0.0),PointBwd(nvariables,0.0);
        //     // Array<OneD, NekDouble> PointFlux(nvariables,0.0);
        //     // Array<OneD, NekDouble> PointNormal(3,0.0);

        //     // DNekMatSharedPtr PointFJac = MemoryManager<DNekMat>
        //     //     ::AllocateSharedPtr(nvariables, nvariables);
        //     // DNekMatSharedPtr PointBJac = MemoryManager<DNekMat>
        //     //     ::AllocateSharedPtr(nvariables, nvariables);

        //     // PointNormal[0] = nx;
        //     // PointNormal[1] = ny;
        //     // PointNormal[2] = nz;


        //     // PointFwd[0] = rhoL;
        //     // PointFwd[1] = rhouL;
        //     // PointFwd[2] = rhovL;
        //     // PointFwd[3] = rhowL;
        //     // PointFwd[4] = EL;

        //     // PointBwd[0] = rhoR;
        //     // PointBwd[1] = rhouR;
        //     // PointBwd[2] = rhovR;
        //     // PointBwd[3] = rhowR;
        //     // PointBwd[4] = ER;
        //     // v_PointFluxJacobian(nvariables,PointFwd,PointBwd,PointNormal,PointFJac,PointBJac);


        //     // DNekMat &M = (*PointBJac);
        //     // NekVector<NekDouble> VectPrim(nvariables,PointBwd,eWrapper);
        //     // NekVector<NekDouble> VectFlux(nvariables,PointFlux,eWrapper);
        //     // VectFlux = M * VectPrim;
        //     // // for(int i =0;i<nvariables;i++)
        //     // // {
        //     // //     PointFlux[i] = 0.0;
        //     // //     for(int j =0;j<nvariables;j++)
        //     // //     {
        //     // //         PointFlux[i] += M(j,i)*VectPrim[j];
        //     // //     }
        //     // // }


        //     // std::cout   <<std::scientific<<std::setw(12)<<std::setprecision(5)
        //     //             <<"abs(PointFlux[0]-f1)   =   "<<abs(PointFlux[0]-f1)<<"    "<<PointFlux[0]<<"    "<<f1<<std::endl
        //     //             <<"abs(PointFlux[1]-f2)   =   "<<abs(PointFlux[1]-f2)<<"    "<<PointFlux[1]<<"    "<<f2<<std::endl
        //     //             <<"abs(PointFlux[2]-f3)   =   "<<abs(PointFlux[2]-f3)<<"    "<<PointFlux[2]<<"    "<<f3<<std::endl
        //     //             <<"abs(PointFlux[3]-f4)   =   "<<abs(PointFlux[3]-f4)<<"    "<<PointFlux[3]<<"    "<<f4<<std::endl
        //     //             <<"abs(PointFlux[4]-f5)   =   "<<abs(PointFlux[4]-f5)<<"    "<<PointFlux[4]<<"    "<<f5<<std::endl;

        // }

        int nvars = 5;
        Array<OneD, NekDouble> Fwd{nvars};
        Array<OneD, NekDouble> Bwd{nvars};
        Array<OneD, NekDouble> flux{nvars, 0.0};

        Fwd[0] = rhoL  ;
        Fwd[1] = rhouL ;
        Fwd[2] = rhovL ;
        Fwd[3] = rhowL ;
        Fwd[4] = EL    ;

        Bwd[0] = rhoR  ;
        Bwd[1] = rhouR ;
        Bwd[2] = rhovR ;
        Bwd[3] = rhowR ;
        Bwd[4] = ER    ;

        v_PointSolve(
        Fwd,
        Bwd,
        flux,
        Fwd,
        Bwd);

        rhof    = flux[0];
        rhouf   = flux[1];
        rhovf   = flux[2];
        rhowf   = flux[3];
        Ef      = flux[4];
    }


    void StegerWarmingSolver::flux_sw_pn(
        double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
        double  nx  , double  ny   , double  nz   ,
        double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef,
        NekDouble efix, NekDouble fsw)
    {
        NekDouble   ro = rhoL;
        NekDouble   vx = rhouL / rhoL;
        NekDouble   vy = rhovL / rhoL;
        NekDouble   vz = rhowL / rhoL;

        // Internal energy (per unit mass)
        NekDouble eL =
                (EL - 0.5 * (rhouL * vx + rhovL * vy + rhowL * vz)) / rhoL;


        NekDouble   ps = m_eos->GetPressure(rhoL, eL);
        NekDouble   c  = m_eos->GetSoundSpeed(rhoL,eL);
        NekDouble   T  = m_eos->GetTemperature(rhoL, eL);
        NekDouble   h  = m_eos->GetEnthalpy(T);
        NekDouble   c2 = c*c;
        NekDouble   v2 = vx*vx + vy*vy + vz*vz;
        NekDouble   h0 = h + 0.5*v2;
        NekDouble   e0 = eL + 0.5*v2;

        NekDouble sml_ssf= 1.0E-12;


        NekDouble   vn = nx*vx + ny*vy + nz*vz;
        // NekDouble   sn = std::numeric_limits::max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf);
        NekDouble   sn = std::max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf);
        NekDouble   osn = 1.0/sn;
        NekDouble   nxa = nx * osn;
        NekDouble   nya = ny * osn;
        NekDouble   nza = nz * osn;
        NekDouble   vna = vn * osn;
        NekDouble   l1 = vn;
        NekDouble   l4 = vn + sn*c;
        NekDouble   l5 = vn - sn*c;

        NekDouble   eps = efix*sn;
        NekDouble   eps2 = eps*eps;
        NekDouble   al1 = sqrt(l1*l1 + eps2);
        NekDouble   al4 = sqrt(l4*l4 + eps2);
        NekDouble   al5 = sqrt(l5*l5 + eps2);

        l1 = 0.5*(l1 + fsw*al1);
        l4 = 0.5*(l4 + fsw*al4);
        l5 = 0.5*(l5 + fsw*al5);


        NekDouble   c2r = c2 / m_eos->GetGamma();
        NekDouble   x1 = c2r * ( 2.0*l1 - l4 - l5 )/( 2.0 * c2 );
        NekDouble   x2 = c2r * ( l4 - l5 )/( 2.0 * c );

        rhof    = (l1 - x1 ) * ro;
        rhouf   = (l1*vx - x1*vx + nxa*x2 ) * ro;
        rhovf   = (l1*vy - x1*vy + nya*x2 ) * ro;
        rhowf   = (l1*vz - x1*vz + nza*x2 ) * ro;
        Ef      = (l1*e0 - x1*h0 + vna*x2 ) * ro;
    }

    void StegerWarmingSolver::v_PointSolve(
        const Array<OneD, NekDouble> &Fwd,
        const Array<OneD, NekDouble> &Bwd,
        Array<OneD, NekDouble>       &flux,
        const Array<OneD, NekDouble> &FwdJ,
        const Array<OneD, NekDouble> &BwdJ)
    {
        int nvars = 5;
        int nspce = 3;

        Array<OneD, NekDouble> normals{nspce,0.0};
        normals[0] = 1.0;

        DNekMatSharedPtr FJac = MemoryManager<DNekMat>::AllocateSharedPtr
                        (nvars,nvars,0.0);
        DNekMatSharedPtr BJac = MemoryManager<DNekMat>::AllocateSharedPtr
                        (nvars,nvars,0.0);

        v_PointFluxJacobian(FwdJ, BwdJ,normals,FJac,BJac);

        Array<OneD, NekDouble> BwdFlux{nvars,0.0};

        NekVector<NekDouble> FwdVect (nvars, Fwd, eWrapper);
        NekVector<NekDouble> BwdVect (nvars, Bwd, eWrapper);
        NekVector<NekDouble> FwdFluxVect (nvars, flux, eWrapper);
        NekVector<NekDouble> BwdFluxVect (nvars, BwdFlux, eWrapper);

        FwdFluxVect = (*FJac)*FwdVect;
        BwdFluxVect = (*BJac)*BwdVect;

        Vmath::Vadd(nvars, BwdFlux, 1, flux, 1, flux, 1);

        return;
    }

    void StegerWarmingSolver::v_PointFluxJacobian(
            const Array<OneD, NekDouble>    &Fwd,
            const Array<OneD, NekDouble>    &Bwd,
            const Array<OneD, NekDouble>    &normals,
            DNekMatSharedPtr                &FJac,
            DNekMatSharedPtr                &BJac)
    {
        NekDouble fsw,efix_StegerWarming;
        efix_StegerWarming = 0.0;

        fsw = 1.0;
        PointFluxJacobian_pn(Fwd,normals,FJac,efix_StegerWarming,fsw);

        fsw = -1.0;
        PointFluxJacobian_pn(Bwd,normals,BJac,efix_StegerWarming,fsw);
        return;
    }
}
