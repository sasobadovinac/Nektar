/////////////////////////////////////////////////////////////////////////////
//
// File GalerkinBoltzmann.cpp
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
// Description: Unsteady linear advection solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <GalerkinBoltzmannSolver/EquationSystems/GalerkinBoltzmann.h>

using namespace std;

namespace Nektar
{
    string GalerkinBoltzmann::className =
        SolverUtils::GetEquationSystemFactory().
        RegisterCreatorFunction("GalerkinBoltzmann",
                                GalerkinBoltzmann::create,
                                "Galerkin Boltzmann equation.");

    GalerkinBoltzmann::GalerkinBoltzmann(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          AdvectionSystem(pSession, pGraph)
    {
        m_planeNumber = 0;
    }

    /**
     * @brief Initialisation object for the Galerkin Boltzmann equation.
     */
    void GalerkinBoltzmann::v_InitObject()
    {
        // Call to the initialisation object of UnsteadySystem
        AdvectionSystem::v_InitObject();

        ASSERTL0(m_fields.num_elements() >= 5,
                 "Expected there to be five variables in Galerkin "
                 "Boltzmann solver (in 2-dimensions)");
        
        ASSERTL0(m_session->DefinesParameter("RT"),"Must define parameter RT");
        m_sqrtRT = sqrt(m_session->GetParameter("RT"));

        ASSERTL0(m_session->DefinesParameter("tau"),"Must define parameter tau");
        m_tau    = m_session->GetParameter("tau");
        
        // Type of advection class to be used
        switch(m_projectionType)
        {
            // Discontinuous field
            case MultiRegions::eDiscontinuous:
            {
                // Do not forwards transform initial condition
                m_homoInitialFwd = false;
                break;
            }
            default:
            {
                ASSERTL0(false, "Unsupported projection type.");
                break;
            }
        }

        // If explicit it computes RHS and PROJECTION for the time integration
        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs     (&GalerkinBoltzmann::DoOdeRhs,        this);
            m_ode.DefineProjection (&GalerkinBoltzmann::DoOdeProjection, this);
        }
        // Otherwise it gives an error (no implicit integration)
        else
        {
            ASSERTL0(false, "Implicit unsteady Advection not set up.");
        }
    }

    /**
     * @brief Galerkin Boltzmann destructor.
     */
    GalerkinBoltzmann::~GalerkinBoltzmann()
    {
    }

    /**
     * @brief Compute the right-hand side for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void GalerkinBoltzmann::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
              Array<OneD,        Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        int ndim       = m_spacedim;

        ASSERTL0(ndim == 2,"Method is not set up for 3D (or 1D) ");

        int nFields = inarray.num_elements();
        int ncoeffs = m_fields[0]->GetNcoeffs();

        // Calculate inner product with respect to derivative of basis of volume flux

        Array<OneD, Array<OneD, NekDouble> > volflux(nFields);
        for(int i = 0; i < nFields; ++i)
        {
            volflux[i] = Array<OneD, NekDouble> (ncoeffs);
        }
        
        EvaluateIProductWRTDerivBaseVolFlux(inarray,volflux);
        
        // Calculate the normal numerical flux on trace boudaries
        // Add upwind normal flux component
        int nTracePointsTot = m_fields[0]->GetTrace()->GetTotPoints();

        // Store forwards/backwards space along trace space
        Array<OneD, Array<OneD, NekDouble> > Fwd    (nFields);
        Array<OneD, Array<OneD, NekDouble> > Bwd    (nFields);
        Array<OneD, Array<OneD, NekDouble> > numflux(nFields);

        for(int i = 0; i < nFields; ++i)
        {
            Fwd[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
            Bwd[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
            numflux[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
            m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
        }

        EvaluateNormalNumFlux(Fwd,Bwd,numflux);

        // Project back to physical space
        for(int i = 0; i < nFields; ++i)
        {
            Vmath::Neg                        (ncoeffs,    volflux[i], 1);
            m_fields[i]->AddTraceIntegral     (numflux[i], volflux[i]);
            m_fields[i]->MultiplyByElmtInvMass(volflux[i], volflux[i]);
            m_fields[i]->BwdTrans             (volflux[i], outarray[i]);
        }

        //-------------------------------------------------------
        // Negate the RHS 
        int nq = m_fields[0]->GetTotPoints();
        for (int i = 0; i < nFields; ++i)
        {
            Vmath::Neg(nq, outarray[i], 1);
        }

        AddSourceTerms(inarray,outarray);
       
    }

    /**
     * @brief Compute the projection for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void GalerkinBoltzmann::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD,       Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        // Counter variable
        int i;

        // Number of fields (variables of the problem)
        int nVariables = inarray.num_elements();


        // Number of quadrature points
        int nQuadraturePts = GetNpoints();
        
        // Just copy over array
        for(i = 0; i < nVariables; ++i)
        {
            Vmath::Vcopy(nQuadraturePts, inarray[i], 1, outarray[i], 1);
        } 

        // Set the boundary conditions
        SetBoundaryConditions(time);
   }

    /**
     * @brief Return the flux vector for the linear advection equation.
     *
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void GalerkinBoltzmann::EvaluateIProductWRTDerivBaseVolFlux(
        const Array<OneD, const  Array<OneD, NekDouble> >&physfield,
        Array<OneD,        Array<OneD, NekDouble> >&volflux)

    {

        NekDouble sqrt_2RT = sqrt(2)*m_sqrtRT;
        int nq = physfield[0].num_elements();

        Array<OneD, Array<OneD, NekDouble> > flux(2);
        flux[0] = Array<OneD, NekDouble>(nq);
        flux[1] = Array<OneD, NekDouble>(nq);

        //  Volume Flux for first component = sqrt(RT) [a1, a2]
        Vmath::Smul(nq, m_sqrtRT, physfield[1], 1, flux[0], 1);
        Vmath::Smul(nq, m_sqrtRT, physfield[2], 1, flux[1], 1);
        m_fields[0]->IProductWRTDerivBase(flux,volflux[0]);
            
        //  Volume Flux for second component =  sqrt(RT)[a0 + sqrt(2) a4, a3]
        Vmath::Smul(nq, m_sqrtRT, physfield[0], 1, flux[0], 1);
        Vmath::Svtvp(nq,sqrt_2RT,physfield[4],1,flux[0],1,flux[0],1);
        Vmath::Smul(nq, m_sqrtRT, physfield[3], 1, flux[1], 1);
        m_fields[1]->IProductWRTDerivBase(flux,volflux[1]);

        //  Volume Flux for third component = sqrt(RT)[a3, a0 + sqrt(2) a5]
        Vmath::Smul(nq, m_sqrtRT, physfield[3], 1, flux[0], 1);
        Vmath::Smul(nq, m_sqrtRT, physfield[0], 1, flux[1], 1);
        Vmath::Svtvp(nq,sqrt_2RT,physfield[5],1,flux[1],1,flux[1],1);
        m_fields[2]->IProductWRTDerivBase(flux,volflux[2]);

        //  Volume Flux for fourth component = sqrt(RT)[a2, a1]
        Vmath::Smul(nq, m_sqrtRT, physfield[2], 1, flux[0], 1);
        Vmath::Smul(nq, m_sqrtRT, physfield[1], 1, flux[1], 1);
        m_fields[3]->IProductWRTDerivBase(flux,volflux[3]);

        //  Volume Flux for fifth component = sqrt(RT)[sqrt(2) a1, 0]
        Vmath::Smul(nq, sqrt_2RT, physfield[1], 1, flux[0], 1);
        Vmath::Zero(nq,  flux[1], 1);
        m_fields[4]->IProductWRTDerivBase(flux,volflux[4]);

        //  Volume Flux for sixth component = sqrt(RT)[0, sqrt(2) a2]
        Vmath::Zero(nq,  flux[0], 1);
        Vmath::Smul(nq, sqrt_2RT, physfield[2], 1, flux[1], 1);
        m_fields[5]->IProductWRTDerivBase(flux,volflux[5]);

    }


    void GalerkinBoltzmann::EvaluateNormalNumFlux(const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
                                                  const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                                                  Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int nTracePointsTot = m_fields[0]->GetTrace()->GetTotPoints();

        // get hold of normals

        NekDouble sqrt2 = sqrt(2);
        NekDouble sqrt3 = sqrt(3);
        NekDouble sqrt6 = sqrt(6);

        for(int i = 0 ; i < nTracePointsTot; ++i)
        {

            NekDouble nx = m_traceNormals[i][0];
            
            if(nx > 0)
            {
                numflux[0][i] = (sqrt3*Fwd[0][i] + 3*Fwd[1][i]
                                 + sqrt6*Fwd[4][i])*nx/6.0 ;

                numflux[1][i] = 0.5*(Fwd[0][i] + sqrt3*Fwd[1][i]
                                     + sqrt2*Fwd[4][i])*nx;
                
                numflux[2][i] = 0.5*(Fwd[2][i] + Fwd[3][i])*nx;

                numflux[3][i] = 0.5*(Fwd[2][i] + Fwd[3][i])*nx;
                
                numflux[4][i] = (sqrt6*Fwd[0][i] + 3*sqrt2*Fwd[1][i]
                                 + 2*sqrt3*Fwd[4][i])*nx/6.0 ;

                numflux[5][i] = 0.0;
            }
            else
            {
                numflux[0][i] = (-sqrt3*Bwd[0][i] + 3*Bwd[1][i]
                                 - sqrt6*Bwd[4][i])*nx/6.0 ;

                numflux[1][i] = 0.5*(Bwd[0][i] - sqrt3*Bwd[1][i]
                                     + sqrt2*Bwd[4][i])*nx;
                
                numflux[2][i] = 0.5*(-Bwd[2][i] + Bwd[3][i])*nx;

                numflux[3][i] = 0.5*(Bwd[2][i] - Bwd[3][i])*nx;
                
                numflux[4][i] = (-sqrt6*Bwd[0][i] + 3*sqrt2*Bwd[1][i]
                                 - 2*sqrt3*Bwd[4][i])*nx/6.0 ;

                numflux[5][i] = 0.0;
                
            }

            NekDouble ny = m_traceNormals[i][1];
            if(ny > 0)
            {
                numflux[0][i] += (sqrt3*Fwd[0][i] + 3*Fwd[2][i]
                                 + sqrt6*Fwd[5][i])*ny/6.0 ;

                numflux[1][i] += 0.5*(Fwd[1][i] + Fwd[3][i])*ny;

                numflux[2][i] += 0.5*(Fwd[0][i] + sqrt3*Fwd[2][i]
                                     + sqrt2*Fwd[5][i])*ny;

                numflux[3][i] += 0.5*(Fwd[1][i] + Fwd[3][i])*ny;
                
                numflux[5][i] += (sqrt6*Fwd[0][i] + 3*sqrt2*Fwd[2][i]
                                 + 2*sqrt3*Fwd[5][i])*ny/6.0 ;

            }
            else
            {
                numflux[0][i] += (-sqrt3*Bwd[0][i] + 3*Bwd[2][i]
                                 - sqrt6*Bwd[5][i])*ny/6.0 ;

                numflux[1][i] += 0.5*(-Bwd[1][i] + Bwd[3][i])*ny;

                numflux[2][i] += 0.5*(Bwd[0][i] - sqrt3*Bwd[2][i]
                                     + sqrt2*Bwd[5][i])*ny;

                numflux[3][i] += 0.5*(Bwd[1][i] - Bwd[3][i])*ny;
                
                numflux[5][i] += (-sqrt6*Bwd[0][i] + 3*sqrt2*Bwd[2][i]
                                 - 2*sqrt3*Bwd[5][i])*ny/6.0 ;

            }
            for(int j = 0; j < 6; ++j)
            {
                numflux[j][i] *= m_sqrtRT;
            }
        }
    }

    void GalerkinBoltzmann::AddSourceTerms( const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                             Array<OneD,        Array<OneD, NekDouble>        >&outarray)
    {

        int nq = m_fields[0]->GetTotPoints();
        
        //-------------------------------------------------
        // Add "source terms"
        // Input and output in physical space so using a collocation projection 
        //-------------------------------------------------
        Array<OneD, NekDouble> force(nq);
        
        NekDouble invsqrt_2 = 1.0/sqrt(2);
        NekDouble invtau = -1.0/m_tau;
        // outarray[3] = -1/tau (a[3] - a[1] a[2]/a[0])
        Vmath::Vmul(nq,inarray[1],1,inarray[2],1,force,1);
        Vmath::Vdiv(nq,force,1,inarray[0],1,force,1);
        Vmath::Vsub(nq,inarray[3],1,force,1,force,1);
        Vmath::Smul(nq,invtau,force,1,force,1);
        Vmath::Vadd(nq,outarray[3],1,force,1,outarray[3],1);

        // outarray[4] = -1/tau (a[4] - 1/sqrt(2) x a[1] a[1]/a[0])
        Vmath::Vmul(nq,inarray[1],1,inarray[1],1,force,1);
        Vmath::Vdiv(nq,force,1,inarray[0],1,force,1);
        Vmath::Smul(nq,invsqrt_2,force,1,force,1);
        Vmath::Vsub(nq,inarray[4],1,force,1,force,1);
        Vmath::Smul(nq,invtau,force,1,force,1);
        Vmath::Vadd(nq,outarray[4],1,force,1,outarray[4],1);

        // outarray[5] = -1/tau (a[5] - 1/sqrt(2) x a[2] a[2]/a[0])
        Vmath::Vmul(nq,inarray[2],1,inarray[2],1,force,1);
        Vmath::Vdiv(nq,force,1,inarray[0],1,force,1);
        Vmath::Smul(nq,invsqrt_2,force,1,force,1);
        Vmath::Vsub(nq,inarray[5],1,force,1,force,1);
        Vmath::Smul(nq,invtau,force,1,force,1);
        Vmath::Vadd(nq,outarray[5],1,force,1,outarray[5],1);
    }
    
    void GalerkinBoltzmann::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        AdvectionSystem::v_GenerateSummary(s);
    }
}
