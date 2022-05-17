/////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyAdvection.cpp
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
// Description: Unsteady linear advection solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <ADRSolver/EquationSystems/UnsteadyAdvection.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/ContField.h>

using namespace std;

namespace Nektar
{
    string UnsteadyAdvection::className = SolverUtils::GetEquationSystemFactory().
        RegisterCreatorFunction("UnsteadyAdvection",
                                UnsteadyAdvection::create,
                                "Unsteady Advection equation.");

    UnsteadyAdvection::UnsteadyAdvection(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          AdvectionSystem(pSession, pGraph), m_planeNumber(0)
    {
    }

    /**
     * @brief Initialisation object for the unsteady linear advection equation.
     */
    void UnsteadyAdvection::v_InitObject()
    {
        // Call to the initialisation object of UnsteadySystem
        AdvectionSystem::v_InitObject();

        m_session->LoadParameter("wavefreq",   m_waveFreq, 0.0);
        // Read the advection velocities from session file

        // check to see if it is explicity turned off
        m_session->MatchSolverInfo("GJPStabilisation", "False",
                                   m_useGJPStabilisation, true);
        
        // if GJPStabilisation set to False bool will be true and
        // if not false so negate/revese bool 
        m_useGJPStabilisation = !m_useGJPStabilisation; 

        m_session->LoadParameter("GJPJumpScale", m_GJPJumpScale, 1.0);

        std::vector<std::string> vel;
        vel.emplace_back("Vx");
        vel.emplace_back("Vy");
        vel.emplace_back("Vz");

        // Resize the advection velocities vector to dimension of the problem
        vel.resize(m_spacedim);

        // Store in the global variable m_velocity the advection velocities
        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        GetFunction( "AdvectionVelocity")->Evaluate(vel,  m_velocity);

        // Type of advection class to be used
        switch(m_projectionType)
        {
            // Continuous field
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                string advName;
                m_session->LoadSolverInfo(
                    "AdvectionType", advName, "NonConservative");
                m_advObject = SolverUtils::
                    GetAdvectionFactory().CreateInstance(advName, advName);
                if (m_specHP_dealiasing)
                {
                    m_advObject->SetFluxVector(
                        &UnsteadyAdvection::GetFluxVectorDeAlias, this);
                }
                else
                {
                    m_advObject->SetFluxVector(
                        &UnsteadyAdvection::GetFluxVector, this);
                }
                break;
            }
            // Discontinuous field
            case MultiRegions::eDiscontinuous:
            {
                // Do not forwards transform initial condition
                m_homoInitialFwd = false;

                // Define the normal velocity fields
                if (m_fields[0]->GetTrace())
                {
                    m_traceVn = Array<OneD, NekDouble>(GetTraceNpoints());
                }

                string advName;
                string riemName;
                m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
                m_advObject = SolverUtils::
                    GetAdvectionFactory().CreateInstance(advName, advName);
                if (m_specHP_dealiasing)
                {
                    m_advObject->SetFluxVector(
                        &UnsteadyAdvection::GetFluxVectorDeAlias, this);
                }
                else
                {
                    m_advObject->SetFluxVector(
                        &UnsteadyAdvection::GetFluxVector, this);
                }
                m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
                m_riemannSolver = SolverUtils::
                    GetRiemannSolverFactory().CreateInstance(
                        riemName, m_session);
                m_riemannSolver->SetScalar(
                    "Vn", &UnsteadyAdvection::GetNormalVelocity, this);
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

        // Forcing terms
        m_forcing = SolverUtils::Forcing::Load(m_session, shared_from_this(),
                                    m_fields, m_fields.size());

        // If explicit it computes RHS and PROJECTION for the time integration
        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs     (&UnsteadyAdvection::DoOdeRhs,        this);
            m_ode.DefineProjection (&UnsteadyAdvection::DoOdeProjection, this);
        }
        // Otherwise it gives an error (no implicit integration)
        else
        {
            ASSERTL0(false, "Implicit unsteady Advection not set up.");
        }
    }

    /**
     * @brief Get the normal velocity for the linear advection equation.
     */
    Array<OneD, NekDouble> &UnsteadyAdvection::GetNormalVelocity()
    {
        // Number of trace (interface) points
        int i;
        int nTracePts = GetTraceNpoints();
        int nPts      = m_velocity[0].size();

        // Auxiliary variable to compute the normal velocity
        Array<OneD, NekDouble> tmp(nPts), tmp2(nTracePts);

        // Reset the normal velocity
        Vmath::Zero(nTracePts, m_traceVn, 1);

        for (i = 0; i < m_velocity.size(); ++i)
        {
            // velocity - grid velocity for ALE before getting trace velocity
            Vmath::Vsub(nPts, m_velocity[i], 1, m_gridVelocity[i], 1, tmp, 1);

            m_fields[0]->ExtractTracePhys(tmp, tmp2);

            Vmath::Vvtvp(nTracePts,
                         m_traceNormals[i], 1,
                         tmp2,              1,
                         m_traceVn,         1,
                         m_traceVn,         1);
        }

        return m_traceVn;
    }

    /**
     * @brief Compute the right-hand side for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyAdvection::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
              Array<OneD,        Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        // Number of fields (variables of the problem)
        int nVariables = inarray.size();

        LibUtilities::Timer timer;
        if (m_ALESolver)
        {
            timer.Start();
            Array<OneD, Array<OneD, NekDouble>> tmpIn(nVariables);
            // If ALE we must take Mu coefficient space to u physical space
            ALEHelper::ALEDoElmtInvMassBwdTrans(inarray, tmpIn);
            m_advObject->AdvectCoeffs(nVariables, m_fields, m_velocity, tmpIn,
                                      outarray, time);
            timer.Stop();
        }
        else
        {
            timer.Start();
            m_advObject->Advect(nVariables, m_fields, m_velocity, inarray, outarray,
                                time);
            timer.Stop();
        }

        // Elapsed time
        timer.AccumulateRegion("Advect");
        
        // Negate the RHS
        for (int i = 0; i < nVariables; ++i)
        {
            Vmath::Neg(outarray[i].size(), outarray[i], 1);
        }

        for (auto &x : m_forcing)
        {
            // set up non-linear terms
            x->Apply(m_fields, inarray, outarray, time);
        }
    }

    /**
     * @brief Compute the projection for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyAdvection::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD,       Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        // Counter variable
        int i;

        // Number of fields (variables of the problem)
        int nVariables = inarray.size();

        // Perform ALE movement
        if (m_ALESolver)
        {
            MoveMesh(time, m_traceNormals);
        }

        // Set the boundary conditions
        SetBoundaryConditions(time);

        // Switch on the projection type (Discontinuous or Continuous)
        switch(m_projectionType)
        {
            // Discontinuous projection
            case MultiRegions::eDiscontinuous:
            {
                // Just copy over array
                for(i = 0; i < nVariables; ++i)
                {
                    Vmath::Vcopy(inarray[i].size(), inarray[i], 1, outarray[i], 1);
                }
                break;
            }

            // Continuous projection
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                int ncoeffs = m_fields[0]->GetNcoeffs();
                Array<OneD, NekDouble> coeffs(ncoeffs,0.0);
                
#if 0
                for(i = 0; i < nVariables; ++i)
                {
                    m_fields[i]->FwdTrans(inarray[i], coeffs);
                    m_fields[i]->BwdTrans_IterPerExp(coeffs, outarray[i]);
                }
#else
                StdRegions::ConstFactorMap factors;
                StdRegions::MatrixType mtype = StdRegions::eMass;

                Array<OneD,NekDouble> wsp(ncoeffs);
                
                for(i = 0; i < nVariables; ++i)
                {
                    MultiRegions::ContFieldSharedPtr cfield =
                        std::dynamic_pointer_cast<
                            MultiRegions::ContField>(m_fields[i]);
                
                    // copy inarray 
                    Array<OneD, NekDouble> in = inarray[i]; 
                    
                    m_fields[i]->IProductWRTBase(in, wsp);
                    
                    if(m_useGJPStabilisation)
                    {
                        const MultiRegions::GJPStabilisationSharedPtr GJPData =
                            cfield->GetGJPForcing(); 
                        
                        factors[StdRegions::eFactorGJP] =
                            m_GJPJumpScale*m_timestep;
                        
                        if(GJPData->IsSemiImplicit())
                        {
                            mtype = StdRegions::eMassGJP;
                        }
                        
                        // to set up forcing need initial guess in
                        // physical space
                        NekDouble scale = -1.0*factors[StdRegions::eFactorGJP];
                        
                        GJPData->Apply(inarray[i],wsp,
                                       NullNekDouble1DArray, scale);
                    }
                
                    // Solve the system
                    MultiRegions::GlobalLinSysKey
                        key(mtype, cfield->GetLocalToGlobalMap(),factors);
                
                    cfield->GlobalSolve(key,wsp,coeffs,
                             NullNekDouble1DArray);
                    
                    m_fields[i]->BwdTrans_IterPerExp(coeffs, outarray[i]);
                }
#endif
                break;
            }
        default:
                ASSERTL0(false,"Unknown projection scheme");
                break;
        }
    }

    /**
     * @brief Return the flux vector for the linear advection equation.
     *
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void UnsteadyAdvection::GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        ASSERTL1(flux[0].size() == m_velocity.size(),
                 "Dimension of flux array and velocity array do not match");

        int i , j;
        int nq = physfield[0].size();

        for (i = 0; i < flux.size(); ++i)
        {
            for (j = 0; j < flux[0].size(); ++j)
            {
                for (int k = 0; k < nq; ++k)
                {
                    // If ALE we need to take off the grid velocity
                    flux[i][j][k] =
                        physfield[i][k] * (m_velocity[j][k] - m_gridVelocity[j][k]);
                }
            }
        }
    }

    /**
     * @brief Return the flux vector for the linear advection equation using
     * the dealiasing technique.
     *
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void UnsteadyAdvection::GetFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        ASSERTL1(flux[0].size() == m_velocity.size(),
                 "Dimension of flux array and velocity array do not match");

        int i, j;
        int nq = physfield[0].size();
        int nVariables = physfield.size();

        // Factor to rescale 1d points in dealiasing
        NekDouble OneDptscale = 2;

        Array<OneD, Array<OneD, NekDouble> >
            advVel_plane(m_velocity.size());

        // Get number of points to dealias a cubic non-linearity
        nq = m_fields[0]->Get1DScaledTotPoints(OneDptscale);

        // Initialisation of higher-space variables
        Array<OneD, Array<OneD, NekDouble> >physfieldInterp(nVariables);
        Array<OneD, Array<OneD, NekDouble> >velocityInterp(m_expdim);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >fluxInterp(nVariables);

        // Interpolation to higher space of physfield
        for (i = 0; i < nVariables; ++i)
        {
            physfieldInterp[i] = Array<OneD, NekDouble>(nq);
            fluxInterp[i] = Array<OneD, Array<OneD, NekDouble> >(m_expdim);
            for (j = 0; j < m_expdim; ++j)
            {
                fluxInterp[i][j] = Array<OneD, NekDouble>(nq);
            }

            m_fields[0]->PhysInterp1DScaled(
                OneDptscale, physfield[i], physfieldInterp[i]);
        }

        // Interpolation to higher space of velocity
        for (j = 0; j < m_expdim; ++j)
        {
            velocityInterp[j] = Array<OneD, NekDouble>(nq);

            m_fields[0]->PhysInterp1DScaled(
                OneDptscale, m_velocity[j], velocityInterp[j]);
        }

        // Evaluation of flux vector in the higher space
        for (i = 0; i < flux.size(); ++i)
        {
            for (j = 0; j < flux[0].size(); ++j)
            {
                Vmath::Vmul(nq, physfieldInterp[i], 1, velocityInterp[j], 1,
                            fluxInterp[i][j], 1);
            }
        }

        // Galerkin project solution back to original space
        for (i = 0; i < nVariables; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                m_fields[0]->PhysGalerkinProjection1DScaled(
                    OneDptscale, fluxInterp[i][j], flux[i][j]);
            }
        }
    }

    void UnsteadyAdvection::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        AdvectionSystem::v_GenerateSummary(s);
        if(m_useGJPStabilisation)
        {
            SolverUtils::AddSummaryItem(s,"GJP Stab. Impl.    ",
                            m_session->GetSolverInfo("GJPStabilisation"));
            SolverUtils::AddSummaryItem(s,"GJP Stab. JumpScale",
                                        m_GJPJumpScale);
        }
    }

    bool UnsteadyAdvection::v_PreIntegrate(int step)
    {
        boost::ignore_unused(step);
        return false;
    }

    void UnsteadyAdvection::v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables)
    {
        bool extraFields;
        m_session->MatchSolverInfo("OutputExtraFields","True",
                                   extraFields, true);

        if (extraFields && m_ALESolver)
        {
            int nCoeffs = m_fields[0]->GetNcoeffs();

            // Adds extra output variables for grid velocity
            std::string gridVarName[3] = {"gridVx", "gridVy", "gridVz"};
            for (int i = 0; i < m_spacedim; ++i)
            {
                Array<OneD, NekDouble> gridVel(nCoeffs, 0.0);
                m_fields[0]->FwdTrans_IterPerExp(m_gridVelocity[i], gridVel);
                fieldcoeffs.emplace_back(gridVel);
                variables.emplace_back(gridVarName[i]);
            }
        }
    }
}
