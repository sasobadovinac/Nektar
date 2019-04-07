/////////////////////////////////////////////////////////////////////////////
//
// File GalkerinBoltzmann.cpp
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
    string GalkerinBoltzmann::className =
        SolverUtils::GetEquationSystemFactory().
        RegisterCreatorFunction("GalkerinBoltzmann",
                                GalkerinBoltzmann::create,
                                "Galerkin Boltzmann equation.");

    GalkerinBoltzmann::GalkerinBoltzmann(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          AdvectionSystem(pSession, pGraph)
    {
        m_planeNumber = 0;
    }

    /**
     * @brief Initialisation object for the Galkerin Boltzmann equation.
     */
    void GalkerinBoltzmann::v_InitObject()
    {
        // Call to the initialisation object of UnsteadySystem
        AdvectionSystem::v_InitObject();

        ASSERTL0(m_fields.num_elements() >= 5,
                 "Expected there to be five variables in Galerkin "
                 "Boltzmann solver (in 2-dimensions)");
        
        ASSERTL0(m_session->DefinesParameter("RT"),"Must define parameter RT");
        m_sqrtRT = sqrt(m_session->GetParameter("RT"));
        
        // Type of advection class to be used
        switch(m_projectionType)
        {
            // Discontinuous field
            case MultiRegions::eDiscontinuous:
            {
                // Do not forwards transform initial condition
                m_homoInitialFwd = false;
                

                string advName;
                string riemName;
                m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
                m_advObject = SolverUtils::
                    GetAdvectionFactory().CreateInstance(advName, advName);

                if (m_specHP_dealiasing)
                {
                    m_advObject->SetFluxVector(
                        &GalkerinBoltzmann::GetFluxVectorDeAlias, this);
                }
                else
                {
                    m_advObject->SetFluxVector(
                        &GalkerinBoltzmann::GetFluxVector, this);
                }

                m_session->LoadSolverInfo(
                    "UpwindType", riemName, "Upwind");
                m_riemannSolver = SolverUtils::
                    GetRiemannSolverFactory().CreateInstance(
                        riemName, m_session);
                m_riemannSolver->SetScalar(
                    "Vn", &GalkerinBoltzmann::GetNormalVelocity, this);

                m_advObject->SetRiemannSolver(m_riemannSolver);
                m_advObject->InitObject(m_session, m_fields);
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
            m_ode.DefineOdeRhs     (&GalkerinBoltzmann::DoOdeRhs,        this);
            m_ode.DefineProjection (&GalkerinBoltzmann::DoOdeProjection, this);
        }
        // Otherwise it gives an error (no implicit integration)
        else
        {
            ASSERTL0(false, "Implicit unsteady Advection not set up.");
        }
    }

    /**
     * @brief Galkerin Boltzmann destructor.
     */
    GalkerinBoltzmann::~GalkerinBoltzmann()
    {
    }

    /**
     * @brief Get the normal velocity 
     */
    Array<OneD, NekDouble> &GalkerinBoltzmann::GetNormalVelocity()
    {
        // Number of trace (interface) points
        int nTracePts = GetTraceNpoints();

        // Auxiliary variable to compute the normal velocity
        Array<OneD, NekDouble> tmp(nTracePts);

    }

    /**
     * @brief Compute the right-hand side for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void GalkerinBoltzmann::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
              Array<OneD,        Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        int ndim       = m_spacedim;
        int nvariables = inarray.num_elements();
        int nq         = GetTotPoints();

        //-------------------------------------------------------
        // Compute the DG advection including the numerical flux
        // by using SolverUtils/Advection 
        // Input and output in physical space
        Array<OneD, Array<OneD, NekDouble> > advVel;
	
        m_advObject->Advect(nvariables, m_fields, advVel, inarray,
                            outarray, time);
        //-------------------------------------------------------

        // Negate the RHS
        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(nq, outarray[i], 1);
        }
        
        //-------------------------------------------------
        // Add "source terms"
        // Input and output in physical space
        //-------------------------------------------------
    }

    /**
     * @brief Compute the projection for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void GalkerinBoltzmann::DoOdeProjection(
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
    void GalkerinBoltzmann::GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        ASSERTL1(flux[0].num_elements() == m_velocity.num_elements(),
                 "Dimension of flux array and velocity array do not match");

        int nq = physfield[0].num_elements();
        NekDouble sqrt_2RT = sqrt(2)*m_sqrtRT;
        
        //  Flux[0] = [a1, a2]
        Vmath::Smul(nq, m_sqrtRT, physfield[1], 1, flux[0][0], 1);
        Vmath::Smul(nq, m_sqrtRT, physfield[2], 1, flux[0][1], 1);

        //  Flux[1] = [a0 + sqrt(2) a4, a3]
        Vmath::Smul(nq, m_sqrtRT, physfield[0], 1, flux[1][0], 1);
        Vmath::Svtvp(nq,sqrt_2RT,physfield[4],1,flux[1][0],1,flux[1][0],1);
        Vmath::Smul(nq, m_sqrtRT, physfield[3], 1, flux[1][1], 1);

        //  Flux[2] = [a3, a0 + sqrt(2) a5]
        Vmath::Smul(nq, m_sqrtRT, physfield[3], 1, flux[2][0], 1);
        Vmath::Smul(nq, m_sqrtRT, physfield[0], 1, flux[2][1], 1);
        Vmath::Svtvp(nq,sqrt_2RT,physfield[5],1,flux[2][1],1,flux[2][1],1);

        //  Flux[3] = [a2, a1]
        Vmath::Smul(nq, m_sqrtRT, physfield[2], 1, flux[3][0], 1);
        Vmath::Smul(nq, m_sqrtRT, physfield[1], 1, flux[3][1], 1);

        //  Flux[4] = [sqrt(2) a1, 0]
        Vmath::Smul(nq, sqrt_2RT, physfield[1], 1, flux[4][0], 1);
        Vmath::Zero(nq,  flux[4][1], 1);

        //  Flux[5] = [0, sqrt(2) a2]
        Vmath::Zero(nq,  flux[5][0], 1);
        Vmath::Smul(nq, sqrt_2RT, physfield[2], 1, flux[5][1], 1);
    }

    /**
     * @brief Return the flux vector for the linear advection equation using
     * the dealiasing technique.
     *
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void GalkerinBoltzmann::GetFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        ASSERTL0(false,"Methods needs setting up");
    }

    void GalkerinBoltzmann::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        AdvectionSystem::v_GenerateSummary(s);
    }
}
