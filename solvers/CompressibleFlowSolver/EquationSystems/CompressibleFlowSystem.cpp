///////////////////////////////////////////////////////////////////////////////
//
// File CompressibleFlowSystem.cpp
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
// Description: Compressible flow system base class with auxiliary functions
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>

#include <LibUtilities/BasicUtils/Timer.h>

using namespace std;

namespace Nektar
{
    CompressibleFlowSystem::CompressibleFlowSystem(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          AdvectionSystem(pSession, pGraph)
    {
    }

    /**
     * @brief Initialization object for CompressibleFlowSystem class.
     */
    void CompressibleFlowSystem::v_InitObject()
    {
        AdvectionSystem::v_InitObject();

        m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(
                    m_session, m_spacedim);

        ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPE"),
                 "No UPWINDTYPE defined in session.");

        // Do not forwards transform initial condition
        m_homoInitialFwd = false;

        // Set up locations of velocity vector.
        m_vecLocs = Array<OneD, Array<OneD, NekDouble> >(1);
        m_vecLocs[0] = Array<OneD, NekDouble>(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_vecLocs[0][i] = 1 + i;
        }

        // Loading parameters from session file
        InitialiseParameters();

        // Setting up advection and diffusion operators
        InitAdvection();

        // Create artificial diffusion
        if (m_shockCaptureType != "Off")
        {
            m_artificialDiffusion = GetArtificialDiffusionFactory()
                                    .CreateInstance(m_shockCaptureType,
                                                    m_session,
                                                    m_fields,
                                                    m_spacedim);
        }

        // Forcing terms for the sponge region
        m_forcing = SolverUtils::Forcing::Load(m_session, shared_from_this(),
                                        m_fields, m_fields.size());

        // User-defined boundary conditions
        int cnt = 0;
        for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
        {
            std::string type =
                m_fields[0]->GetBndConditions()[n]->GetUserDefined();

            if (m_fields[0]->GetBndConditions()[n]->GetBoundaryConditionType()
                == SpatialDomains::ePeriodic)
            {
                continue;
            }

            if (!type.empty())
            {
                m_bndConds.push_back(GetCFSBndCondFactory().CreateInstance(
                        type,
                        m_session,
                        m_fields,
                        m_traceNormals,
                        m_spacedim,
                        n,
                        cnt));
            }
            cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
        }

        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs    (&CompressibleFlowSystem::DoOdeRhs, this);
            m_ode.DefineProjection(&CompressibleFlowSystem::DoOdeProjection, this);
        }
        else
        {   
            m_ode.DefineOdeRhs    (&CompressibleFlowSystem::DoOdeRhs, this);
            m_ode.DefineProjection(&CompressibleFlowSystem::
                                   DoOdeProjection, this);
            m_ode.DefineImplicitSolve(&CompressibleFlowSystem::
                                      DoImplicitSolvePhysToCoeff, this);

            InitialiseNonlinSysSolver();
        }

        SetBoundaryConditionsBwdWeight();

        string advName;
        m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
    }

    void CompressibleFlowSystem::InitialiseNonlinSysSolver()
    {
        std::string SovlerType = "Newton";
        if (m_session->DefinesSolverInfo("NonlinSysIterSovler"))
        {
            SovlerType = m_session->GetSolverInfo("NonlinSysIterSovler");
        }
        ASSERTL0(LibUtilities::GetNekNonlinSysFactory().
                 ModuleExists(SovlerType),
                "NekNonlinSys '" + SovlerType + "' is not defined.\n");
        LibUtilities::CommSharedPtr v_Comm  = 
                                    m_fields[0]->GetComm()->GetRowComm();
        int ntotal = m_fields[0]->GetNcoeffs() * m_fields.size();
        m_nonlinsol = LibUtilities::GetNekNonlinSysFactory().CreateInstance(
                            SovlerType, m_session, v_Comm,ntotal);

        m_NekSysOp.DefineNekSysRhsEval(&CompressibleFlowSystem
                                                  ::NonlinSysEvaluator1D, this);
        m_NekSysOp.DefineNekSysLhsEval(&CompressibleFlowSystem::
                                         MatrixMultiply_MatrixFreeCoeff, this);
        m_nonlinsol->setSysOperators(m_NekSysOp);
    }

    /**
     * @brief Destructor for CompressibleFlowSystem class.
     */
    CompressibleFlowSystem::~CompressibleFlowSystem()
    {

    }

    /**
     * @brief Load CFS parameters from the session file
     */
    void CompressibleFlowSystem::InitialiseParameters()
    {
        // Get gamma parameter from session file.
        m_session->LoadParameter("Gamma", m_gamma, 1.4);

        // Shock capture
        m_session->LoadSolverInfo("ShockCaptureType",
                                  m_shockCaptureType, "Off");

        // Load parameters for exponential filtering
        m_session->MatchSolverInfo("ExponentialFiltering","True",
                                   m_useFiltering, false);
        if (m_useFiltering)
        {
            m_session->LoadParameter ("FilterAlpha", m_filterAlpha, 36);
            m_session->LoadParameter ("FilterExponent", m_filterExponent, 16);
            m_session->LoadParameter ("FilterCutoff", m_filterCutoff, 0);
        }

        // Load CFL for local time-stepping (for steady state)
        m_session->MatchSolverInfo("LocalTimeStep","True",
                                   m_useLocalTimeStep, false);
        if (m_useLocalTimeStep)
        {
            ASSERTL0(m_cflSafetyFactor != 0,
                    "Local time stepping requires CFL parameter.");
        }
        m_session->LoadParameter ("JacobiFreeEps", m_JacobiFreeEps, 5.0E-8);
        m_session->LoadParameter("NewtonAbsoluteIteTol", 
                                 m_NewtonAbsoluteIteTol, 1.0E-12);
    }

    /**
     * @brief Create advection and diffusion objects for CFS
     */
    void CompressibleFlowSystem::InitAdvection()
    {
        // Check if projection type is correct
        ASSERTL0(m_projectionType == MultiRegions::eDiscontinuous,
                "Unsupported projection type.");

        string advName, riemName;
        m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");

        m_advObject = SolverUtils::GetAdvectionFactory()
                                    .CreateInstance(advName, advName);

        if (m_specHP_dealiasing)
        {
            m_advObject->SetFluxVector(&CompressibleFlowSystem::
                                       GetFluxVectorDeAlias, this);
        }
        else
        {
            m_advObject->SetFluxVector  (&CompressibleFlowSystem::
                                          GetFluxVector, this);
        }

        // Setting up Riemann solver for advection operator
        m_session->LoadSolverInfo("UpwindType", riemName, "Average");

        SolverUtils::RiemannSolverSharedPtr riemannSolver;
        riemannSolver = SolverUtils::GetRiemannSolverFactory()
                                    .CreateInstance(riemName, m_session);

        // Setting up parameters for advection operator Riemann solver
        riemannSolver->SetParam (
            "gamma",   &CompressibleFlowSystem::GetGamma,   this);
        riemannSolver->SetAuxVec(
            "vecLocs", &CompressibleFlowSystem::GetVecLocs, this);
        riemannSolver->SetVector(
            "N",       &CompressibleFlowSystem::GetNormals, this);
        riemannSolver->SetVector(
            "vgt",     &ALEHelper::GetGridVelocityTrace, this);

        // Concluding initialisation of advection / diffusion operators
        m_advObject->SetRiemannSolver   (riemannSolver);
        m_advObject->InitObject         (m_session, m_fields);
    }

    /**
     * @brief Compute the right-hand side.
     */
    void CompressibleFlowSystem::DoOdeRhs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int nvariables = inarray.size();
        //int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

        m_BndEvaluateTime   = time;

        // Store forwards/backwards space along trace space
        Array<OneD, Array<OneD, NekDouble> > Fwd    (nvariables);
        Array<OneD, Array<OneD, NekDouble> > Bwd    (nvariables);

        if (m_HomogeneousType == eHomogeneous1D)
        {
            Fwd = NullNekDoubleArrayofArray;
            Bwd = NullNekDoubleArrayofArray;
        }
        else
        {
            for (int i = 0; i < nvariables; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
        }

        // Calculate advection
LibUtilities::Timer timer;
timer.Start();
        DoAdvection(inarray, outarray, time, Fwd, Bwd);
timer.Stop();
timer.AccumulateRegion("DoAdvection");

        // Negate results
        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(outarray[i].size(), outarray[i], 1);
        }

        // Add diffusion terms
timer.Start();
        DoDiffusion(inarray, outarray, Fwd, Bwd);
timer.Stop();
timer.AccumulateRegion("DoDiffusion");

        // Add forcing terms
        for (auto &x : m_forcing)
        {
            x->Apply(m_fields, inarray, outarray, time);
        }

        if (m_useLocalTimeStep)
        {
            int nElements = m_fields[0]->GetExpSize();
            int nq, offset;
            NekDouble fac;
            Array<OneD, NekDouble> tmp;

            Array<OneD, NekDouble> tstep (nElements, 0.0);
            GetElmtTimeStep(inarray, tstep);

            // Loop over elements
            for (int n = 0; n < nElements; ++n)
            {
                nq     = m_fields[0]->GetExp(n)->GetTotPoints();
                offset = m_fields[0]->GetPhys_Offset(n);
                fac    = tstep[n] / m_timestep;
                for (int i = 0; i < nvariables; ++i)
                {
                    Vmath::Smul(nq, fac, outarray[i] + offset, 1,
                                         tmp = outarray[i] + offset, 1);
                }
            }
        }
    }

    /**
     * @brief Compute the projection and call the method for imposing the
     * boundary conditions in case of discontinuous projection.
     */
    void CompressibleFlowSystem::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int nvariables = inarray.size();

        // Perform ALE movement
        if (m_ALESolver)
        {
            MoveMesh(time, m_traceNormals);
        }

        switch(m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                // Just copy over array
                //int npoints = GetNpoints();

                for (int i = 0; i < nvariables; ++i)
                {
                    //std::cout << "POINTS??? : " << GetNpoints() << " " << inarray[i].size() << " " << outarray[i].size() << std::endl;
                    Vmath::Vcopy(inarray[i].size(), inarray[i], 1, outarray[i], 1); // @TODO: This used to be npoints not inarray[i].size(), not sure why, something to do with ALEHelper::ALEPreMultiplyMass(fields) changing field size

                    if (m_useFiltering)
                    {
                        m_fields[i]->ExponentialFilter(outarray[i],
                            m_filterAlpha, m_filterExponent, m_filterCutoff);
                    }
                }
                SetBoundaryConditions(outarray, time);
                break;
            }
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                ASSERTL0(false, "No Continuous Galerkin for full compressible "
                                "Navier-Stokes equations");
                break;
            }
            default:
                ASSERTL0(false, "Unknown projection scheme");
                break;
        }
    }

    /**
     * @brief Compute the advection terms for the right-hand side
     */
    void CompressibleFlowSystem::DoAdvection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time,
        const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >       &pBwd)
    {
        int nvariables = inarray.size();
        Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);

        if(m_ALESolver)
        {
            ALEHelper::ALEDoOdeRhs(inarray, outarray, time, m_advObject,
                                   advVel);
        }
        else
        {
            m_advObject->Advect(nvariables, m_fields, advVel, inarray, outarray,
                                time, pFwd, pBwd);
        }
    }

    /**
     * @brief Add the diffusions terms to the right-hand side
     */
    void CompressibleFlowSystem::DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd)
    {
        v_DoDiffusion(inarray, outarray, pFwd, pBwd);

    }

    void CompressibleFlowSystem::NonlinSysEvaluator1D(
        const TensorOfArray1D<NekDouble>    &inarray,
        TensorOfArray1D<NekDouble>          &out,
        const bool                          &flag)
    {
        boost::ignore_unused(flag);
        unsigned int nvariables     = m_fields.size();
        unsigned int npoints        = m_fields[0]->GetNcoeffs();
        TensorOfArray2D<NekDouble> in2D(nvariables);
        TensorOfArray2D<NekDouble> out2D(nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            int offset = i * npoints;
            in2D[i]    = inarray + offset;
            out2D[i]   = out + offset;
        }
        NonlinSysEvaluatorCoeff(in2D, out2D);
    }

    void CompressibleFlowSystem::NonlinSysEvaluatorCoeff(
        TensorOfArray2D<NekDouble>  &inarray,
        TensorOfArray2D<NekDouble>  &out)
    {
        const Array<OneD, const NekDouble> refsol 
                                = m_nonlinsol->GetRefSolution();
        unsigned int nvariable  = inarray.size();
        unsigned int ncoeffs    = inarray[nvariable - 1].size();
        unsigned int npoints    = m_fields[0]->GetNpoints();

        TensorOfArray2D<NekDouble> inpnts(nvariable);

        for (int i = 0; i < nvariable; ++i)
        {
            inpnts[i] = TensorOfArray1D<NekDouble>(npoints, 0.0);
            m_fields[i]->BwdTrans(inarray[i], inpnts[i]);
        }
        
        DoOdeProjection(inpnts, inpnts, m_BndEvaluateTime);
        DoOdeRhsCoeff(inpnts, out, m_BndEvaluateTime);

        for (int i = 0; i < nvariable; ++i)
        {
            int noffset = i * ncoeffs;
            Vmath::Svtvp(ncoeffs, m_TimeIntegLambda, out[i], 1,
                         refsol + noffset, 1, out[i], 1);
            Vmath::Vsub(ncoeffs, inarray[i], 1,
                        out[i], 1, out[i], 1);
        }
        return;
    }

    /**
     * @brief Compute the right-hand side.
     */
    void CompressibleFlowSystem::DoOdeRhsCoeff(
        const TensorOfArray2D<NekDouble>    &inarray,
        TensorOfArray2D<NekDouble>          &outarray,
        const NekDouble                     time)
    {
        int nvariables = inarray.size();
        int nTracePts  = GetTraceTotPoints();
        int ncoeffs    = GetNcoeffs();
        // Store forwards/backwards space along trace space
        TensorOfArray2D<NekDouble> Fwd    (nvariables);
        TensorOfArray2D<NekDouble> Bwd    (nvariables);

        if (m_HomogeneousType == eHomogeneous1D)
        {
            Fwd = NullNekDoubleArrayofArray;
            Bwd = NullNekDoubleArrayofArray;
        }
        else
        {
            for (int i = 0; i < nvariables; ++i)
            {
                Fwd[i]     = TensorOfArray1D<NekDouble> (nTracePts, 0.0);
                Bwd[i]     = TensorOfArray1D<NekDouble> (nTracePts, 0.0);
                m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
        }
 
         // Calculate advection
        DoAdvectionCoeff(inarray, outarray, time, Fwd, Bwd);

        // Negate results
        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(ncoeffs, outarray[i], 1);
        }

        DoDiffusionCoeff(inarray, outarray, Fwd, Bwd);

        // Add forcing terms
        for (auto &x : m_forcing)
        {
            x->ApplyCoeff(m_fields, inarray, outarray, time);
        }

        if (m_useLocalTimeStep)
        {
            int nElements = m_fields[0]->GetExpSize();
            int nq, offset;
            NekDouble fac;
            TensorOfArray1D<NekDouble> tmp;

            TensorOfArray1D<NekDouble> tstep(nElements, 0.0);
            GetElmtTimeStep(inarray, tstep);

            // Loop over elements
            for (int n = 0; n < nElements; ++n)
            {
                nq     = m_fields[0]->GetExp(n)->GetNcoeffs();
                offset = m_fields[0]->GetCoeff_Offset(n);
                fac    = tstep[n] / m_timestep;
                for (int i = 0; i < nvariables; ++i)
                {
                    Vmath::Smul(nq, fac, outarray[i] + offset, 1,
                                tmp = outarray[i] + offset, 1);
                }
            }
        }
    }

    /**
     * @brief Compute the advection terms for the right-hand side
     */
    void CompressibleFlowSystem::DoAdvectionCoeff(
        const TensorOfArray2D<NekDouble>    &inarray,
        TensorOfArray2D<NekDouble>          &outarray,
        const NekDouble                     time,
        const TensorOfArray2D<NekDouble>    &pFwd,
        const TensorOfArray2D<NekDouble>    &pBwd)
    {
        int nvariables = inarray.size();
        TensorOfArray2D<NekDouble> advVel(m_spacedim);

        m_advObject->AdvectCoeffs(nvariables, m_fields, advVel, inarray,
                                  outarray, time, pFwd, pBwd);

        // Multiply by inverse mass matrix
        LibUtilities::Timer timer;
        for (int i = 0; i < nvariables; ++i)
        {
            timer.Start();
            m_fields[i]->MultiplyByElmtInvMass(outarray[i], outarray[i]);
            timer.Stop();
            timer.AccumulateRegion("MultiplyByElmtInvMass");
        }
    }

    /**
     * @brief Add the diffusions terms to the right-hand side
     * Similar to DoDiffusion() but with outarray in coefficient space
     */
    void CompressibleFlowSystem::DoDiffusionCoeff(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd)
    {
        v_DoDiffusionCoeff(inarray, outarray, pFwd, pBwd);
    }

    void CompressibleFlowSystem::DoImplicitSolvePhysToCoeff(
        const TensorOfArray2D<NekDouble>    &inpnts,
        TensorOfArray2D<NekDouble>          &outpnt,
        const NekDouble                     time,
        const NekDouble                     lambda)
    {
        unsigned int nvariables  = inpnts.size();
        unsigned int ncoeffs     = m_fields[0]->GetNcoeffs();
        unsigned int ntotal      = nvariables * ncoeffs;

        TensorOfArray1D<NekDouble>  inarray(ntotal);
        TensorOfArray1D<NekDouble>  out(ntotal);
        TensorOfArray1D<NekDouble>  tmpArray;

        for (int i = 0; i < nvariables; ++i)
        {
            int noffset = i * ncoeffs;
            tmpArray = inarray + noffset;
            m_fields[i]->FwdTrans(inpnts[i], tmpArray);
        }

        DoImplicitSolveCoeff(inpnts, inarray, out, time, lambda);

        for (int i = 0; i < nvariables; ++i)
        {
            int noffset = i * ncoeffs;
            tmpArray = out + noffset;
            m_fields[i]->BwdTrans(tmpArray, outpnt[i]);
        }
    }

    void CompressibleFlowSystem::DoImplicitSolveCoeff(
            const TensorOfArray2D<NekDouble>    &inpnts,
            const TensorOfArray1D<NekDouble>    &inarray,
            TensorOfArray1D<NekDouble>          &out,
            const NekDouble                     time,
            const NekDouble                     lambda)
    {
        boost::ignore_unused(inpnts);
        m_TimeIntegLambda               = lambda;
        m_BndEvaluateTime               = time;
        unsigned int ntotal             = inarray.size();

        if (m_inArrayNorm < 0.0)
        {
            CalcRefValues(inarray);
        }
        
        NekDouble tol2 = m_inArrayNorm
                        *m_NewtonAbsoluteIteTol * m_NewtonAbsoluteIteTol;
        m_TotNewtonIts +=  m_nonlinsol->SolveSystem(ntotal,inarray,
                                                    out, 0, tol2);
        
        m_TotImpStages++;
        m_StagesPerStep++;
    }

    void CompressibleFlowSystem::CalcRefValues(
            const TensorOfArray1D<NekDouble>    &inarray)
    {
        unsigned int nvariables         = m_fields.size();
        unsigned int ntotal             = inarray.size();
        unsigned int npoints            = ntotal/nvariables;

        unsigned int ntotalGlobal       = ntotal;
        m_comm->AllReduce(ntotalGlobal, Nektar::LibUtilities::ReduceSum);
        unsigned int ntotalDOF          = ntotalGlobal / nvariables;
        NekDouble ototalDOF             = 1.0 / ntotalDOF;

        m_inArrayNorm = 0.0;
        m_magnitdEstimat = Array<OneD, NekDouble>  (nvariables, 0.0);

        for (int i = 0; i < nvariables; ++i)
        {
            int offset = i * npoints;
            m_magnitdEstimat[i] = Vmath::Dot(npoints, inarray + offset,
                                            inarray + offset);
        }
        m_comm->AllReduce(m_magnitdEstimat, Nektar::LibUtilities::ReduceSum);

        for (int i = 0; i < nvariables; ++i)
        {
            m_inArrayNorm += m_magnitdEstimat[i];
        }

        for (int i = 2; i < nvariables - 1; ++i)
        {
            m_magnitdEstimat[1]   +=   m_magnitdEstimat[i] ;
        }
        for (int i = 2; i < nvariables - 1; ++i)
        {
            m_magnitdEstimat[i]   =   m_magnitdEstimat[1] ;
        }

        for (int i = 0; i < nvariables; ++i)
        {
            m_magnitdEstimat[i] = sqrt(m_magnitdEstimat[i] * ototalDOF);
        }
        bool l_verbose      = m_session->DefinesCmdLineArgument("verbose");
        if (0 == m_session->GetComm()->GetRank() && l_verbose)
        {
            for (int i = 0; i < nvariables; ++i)
            {
                cout << "m_magnitdEstimat[" << i << "]    = "
                     << m_magnitdEstimat[i] << endl;
            }
            cout << "m_inArrayNorm    = " << m_inArrayNorm << endl;
        }
    }

    void CompressibleFlowSystem::MatrixMultiply_MatrixFreeCoeff(
            const TensorOfArray1D<NekDouble>    &inarray,
            TensorOfArray1D<NekDouble>          &out,
            const bool                          &flag)
    {
        boost::ignore_unused(flag);
        const Array<OneD, const NekDouble> refsol 
                                         = m_nonlinsol->GetRefSolution();
        const Array<OneD, const NekDouble> refres 
                                         = m_nonlinsol->GetRefResidual();
        NekDouble eps = m_JacobiFreeEps;
        NekDouble magnitdEstimatMax = 0.0;
        for (int i = 0; i < m_magnitdEstimat.size(); ++i)
        {
            magnitdEstimatMax = max(magnitdEstimatMax, m_magnitdEstimat[i]);
        }
        eps *= magnitdEstimatMax;
        NekDouble oeps = 1.0 / eps;
        unsigned int nvariables = m_fields.size();
        unsigned int ntotal     = inarray.size();
        unsigned int npoints    = ntotal / nvariables;
        Array<OneD, NekDouble > tmp;
        TensorOfArray2D<NekDouble> solplus(nvariables);
        TensorOfArray2D<NekDouble> resplus(nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            solplus[i] = TensorOfArray1D<NekDouble> (npoints, 0.0);
            resplus[i] = TensorOfArray1D<NekDouble> (npoints, 0.0);
        }

        for (int i = 0; i < nvariables; ++i)
        {
            int noffset = i*npoints;
            tmp = inarray + noffset;
            Vmath::Svtvp(npoints, eps, tmp, 1, 
                         refsol + noffset,1, solplus[i], 1);
        }
        
        NonlinSysEvaluatorCoeff(solplus,resplus);

        for (int i = 0; i < nvariables; ++i)
        {
            int noffset = i * npoints;
            tmp = out + noffset;
            Vmath::Vsub(npoints, &resplus[i][0], 1, 
                       &refres[0] + noffset, 1, &tmp[0],1);
            Vmath::Smul(npoints, oeps, &tmp[0], 1, &tmp[0], 1);
        }
       
        return;
    }

    void CompressibleFlowSystem::SetBoundaryConditions(
            Array<OneD, Array<OneD, NekDouble> >             &physarray,
            NekDouble                                         time)
    {

        return; // @TODO: This is only for no boundary case - REMOVE IN FUTURE!!!
        int nTracePts  = GetTraceTotPoints();
        int nvariables = physarray.size();

        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }

        if (m_bndConds.size())
        {
            // Loop over user-defined boundary conditions
            for (auto &x : m_bndConds)
            {
                x->Apply(Fwd, physarray, time);
            }
        }
    }
    /**
     * @brief Set up a weight on physical boundaries for boundary condition
     * applications
     */
    void CompressibleFlowSystem::SetBoundaryConditionsBwdWeight()
    {
        if (m_bndConds.size())
        {
            // Loop over user-defined boundary conditions
            for (auto &x : m_bndConds)
            {
                x->ApplyBwdWeight();
            }
        }
    }

    /**
     * @brief Return the flux vector for the compressible Euler equations.
     *
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >  &physfield,
        TensorOfArray3D<NekDouble>                  &flux)
    {
        int i, j;
        int nq = physfield[0].size();
        int nVariables = m_fields.size();

        Array<OneD, NekDouble> pressure(nq);
        Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);

        // Flux vector for the rho equation
        for (i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nq);
            Vmath::Vcopy(nq, physfield[i+1], 1, flux[0][i], 1);
        }

        m_varConv->GetVelocityVector(physfield, velocity);
        m_varConv->GetPressure(physfield, pressure);

        // Flux vector for the velocity fields
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vmul(nq, velocity[j], 1, physfield[i+1], 1,
                            flux[i+1][j], 1);
            }

            // Add pressure to appropriate field
            Vmath::Vadd(nq, flux[i+1][i], 1, pressure, 1, flux[i+1][i], 1);
        }

        // Flux vector for energy.
        Vmath::Vadd(nq, physfield[m_spacedim+1], 1, pressure, 1,
                    pressure, 1);

        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nq, velocity[j], 1, pressure, 1,
                        flux[m_spacedim+1][j], 1);
        }

        // @TODO : for each row (3 columns) negative grid velocity component (d * d + 2 (4 rows)) rho, rhou, rhov, rhow, E,
        // @TODO : top row is flux for rho etc... each row subtract v_g * conserved variable for that row...
        // For grid velocity subtract v_g * conserved variable
        for (i = 0; i < m_spacedim + 2; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                for (int k = 0; k < nq; ++k)
                {
                    flux[i][j][k] -= physfield[i][k] * m_gridVelocity[j][k];
                }
            }
        }

        // For the smooth viscosity model
        if (nVariables == m_spacedim+3)
        {
            // Add a zero row for the advective fluxes
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Zero(nq, flux[m_spacedim+2][j], 1);
            }
        }
    }

    /**
     * @brief Return the flux vector for the compressible Euler equations
     * by using the de-aliasing technique.
     *
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >      &physfield,
        TensorOfArray3D<NekDouble>                      &flux)
    {
        int i, j;
        int nq = physfield[0].size();
        int nVariables = m_fields.size();

        // Factor to rescale 1d points in dealiasing
        NekDouble OneDptscale = 2;
        nq = m_fields[0]->Get1DScaledTotPoints(OneDptscale);

        Array<OneD, NekDouble> pressure(nq);
        Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);

        Array<OneD, Array<OneD, NekDouble> > physfield_interp(nVariables);
        TensorOfArray3D<NekDouble> flux_interp(nVariables);

        for (i = 0; i < nVariables; ++ i)
        {
            physfield_interp[i] = Array<OneD, NekDouble>(nq);
            flux_interp[i] = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
            m_fields[0]->PhysInterp1DScaled(
                OneDptscale, physfield[i], physfield_interp[i]);

            for (j = 0; j < m_spacedim; ++j)
            {
                flux_interp[i][j] = Array<OneD, NekDouble>(nq);
            }
        }

        // Flux vector for the rho equation
        for (i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nq);

            // Galerkin project solution back to original space
            m_fields[0]->PhysGalerkinProjection1DScaled(
                OneDptscale, physfield_interp[i+1], flux[0][i]);
        }

        m_varConv->GetVelocityVector(physfield_interp, velocity);
        m_varConv->GetPressure      (physfield_interp, pressure);

        // Evaluation of flux vector for the velocity fields
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vmul(nq, velocity[j], 1, physfield_interp[i+1], 1,
                            flux_interp[i+1][j], 1);
            }

            // Add pressure to appropriate field
            Vmath::Vadd(nq, flux_interp[i+1][i], 1, pressure,1,
                        flux_interp[i+1][i], 1);
        }

        // Galerkin project solution back to original space
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                m_fields[0]->PhysGalerkinProjection1DScaled(
                    OneDptscale, flux_interp[i+1][j], flux[i+1][j]);
            }
        }

        // Evaluation of flux vector for energy
        Vmath::Vadd(nq, physfield_interp[m_spacedim+1], 1, pressure, 1,
                    pressure, 1);

        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nq, velocity[j], 1, pressure, 1,
                        flux_interp[m_spacedim+1][j], 1);

            // Galerkin project solution back to original space
            m_fields[0]->PhysGalerkinProjection1DScaled(
                OneDptscale,
                flux_interp[m_spacedim+1][j],
                flux[m_spacedim+1][j]);
        }
    }

    /**
     * @brief Calculate the maximum timestep on each element
     *        subject to CFL restrictions.
     */
    void CompressibleFlowSystem::GetElmtTimeStep(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD, NekDouble> &tstep)
    {
        boost::ignore_unused(inarray);

        int nElements = m_fields[0]->GetExpSize();

        // Change value of m_timestep (in case it is set to zero)
        NekDouble tmp = m_timestep;
        m_timestep    = 1.0;

        Array<OneD, NekDouble> cfl(nElements);
        cfl = GetElmtCFLVals();

        // Factors to compute the time-step limit
        NekDouble alpha     = MaxTimeStepEstimator();

        // Loop over elements to compute the time-step limit for each element
        for (int n = 0; n < nElements; ++n)
        {
            tstep[n] = m_cflSafetyFactor * alpha / cfl[n];
        }

        // Restore value of m_timestep
        m_timestep = tmp;
    }

    /**
     * @brief Calculate the maximum timestep subject to CFL restrictions.
     */
    NekDouble CompressibleFlowSystem::v_GetTimeStep(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray)
    {
        int nElements = m_fields[0]->GetExpSize();
        Array<OneD, NekDouble> tstep (nElements, 0.0);

        GetElmtTimeStep(inarray, tstep);

        // Get the minimum time-step limit and return the time-step
        NekDouble TimeStep = Vmath::Vmin(nElements, tstep, 1);
        m_comm->AllReduce(TimeStep, LibUtilities::ReduceMin);
        return TimeStep;
    }

    /**
     * @brief Set up logic for residual calculation.
     */
    void CompressibleFlowSystem::v_SetInitialConditions(
        NekDouble initialtime,
        bool      dumpInitialConditions,
        const int domain)
    {
        boost::ignore_unused(domain);

        EquationSystem::v_SetInitialConditions(initialtime, false);

        // insert white noise in initial condition
        NekDouble Noise;
        int phystot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> noise(phystot);

        m_session->LoadParameter("Noise", Noise,0.0);
        int m_nConvectiveFields =  m_fields.size();

        if (Noise > 0.0)
        {
            int seed = - m_comm->GetRank()*m_nConvectiveFields;
            for (int i = 0; i < m_nConvectiveFields; i++)
            {
                Vmath::FillWhiteNoise(phystot, Noise, noise, 1,
                                      seed);
                --seed;
                Vmath::Vadd(phystot, m_fields[i]->GetPhys(), 1,
                            noise, 1, m_fields[i]->UpdatePhys(), 1);
                m_fields[i]->FwdTrans_IterPerExp(m_fields[i]->GetPhys(),
                                                 m_fields[i]->UpdateCoeffs());
            }
        }

        if (dumpInitialConditions && m_checksteps)
        {
            Checkpoint_Output(m_nchk);
            m_nchk++;
        }
    }

    /**
     * @brief Compute the advection velocity in the standard space
     * for each element of the expansion.
     */
    Array<OneD, NekDouble> CompressibleFlowSystem::v_GetMaxStdVelocity()
    {
        int nTotQuadPoints = GetTotPoints();
        int n_element      = m_fields[0]->GetExpSize();
        int expdim         = m_fields[0]->GetGraph()->GetMeshDimension();
        int nfields        = m_fields.size();
        int offset;
        Array<OneD, NekDouble> tmp;

        Array<OneD, Array<OneD, NekDouble> > physfields(nfields);
        for (int i = 0; i < nfields; ++i)
        {
            physfields[i] = m_fields[i]->GetPhys();
        }

        Array<OneD, NekDouble> stdV(n_element, 0.0);

        // Getting the velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > velocity   (m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > stdSoundSpeed(m_spacedim);
        Array<OneD, NekDouble>               soundspeed (nTotQuadPoints);
        LibUtilities::PointsKeyVector        ptsKeys;

        for (int i = 0; i < m_spacedim; ++i)
        {
            velocity   [i]   = Array<OneD, NekDouble>(nTotQuadPoints);
            stdVelocity[i]   = Array<OneD, NekDouble>(nTotQuadPoints, 0.0);
            stdSoundSpeed[i] = Array<OneD, NekDouble>(nTotQuadPoints, 0.0);
        }

        m_varConv->GetVelocityVector(physfields, velocity);
        m_varConv->GetSoundSpeed    (physfields, soundspeed);

        for (int el = 0; el < n_element; ++el)
        {
            ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();
            offset  = m_fields[0]->GetPhys_Offset(el);
            int nq = m_fields[0]->GetExp(el)->GetTotPoints();

            const SpatialDomains::GeomFactorsSharedPtr metricInfo =
                m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo();
            const Array<TwoD, const NekDouble> &gmat =
                m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()
                                                  ->GetDerivFactors(ptsKeys);

            // Convert to standard element
            //    consider soundspeed in all directions
            //    (this might overestimate the cfl)
            if (metricInfo->GetGtype() == SpatialDomains::eDeformed)
            {
                // d xi/ dx = gmat = 1/J * d x/d xi
                for (int i = 0; i < expdim; ++i)
                {
                    Vmath::Vmul(nq, gmat[i], 1,
                                    velocity[0] + offset, 1,
                                    tmp = stdVelocity[i] + offset, 1);
                    Vmath::Vmul(nq, gmat[i], 1,
                                    soundspeed + offset, 1,
                                    tmp = stdSoundSpeed[i] + offset, 1);
                    for (int j = 1; j < expdim; ++j)
                    {
                        Vmath::Vvtvp(nq, gmat[expdim*j+i], 1,
                                         velocity[j] + offset, 1,
                                         stdVelocity[i] + offset, 1,
                                         tmp = stdVelocity[i] + offset, 1);
                        Vmath::Vvtvp(nq, gmat[expdim*j+i], 1,
                                         soundspeed + offset, 1,
                                         stdSoundSpeed[i] + offset, 1,
                                         tmp = stdSoundSpeed[i] + offset, 1);
                    }
                }
            }
            else
            {
                for (int i = 0; i < expdim; ++i)
                {
                    Vmath::Smul(nq, gmat[i][0],
                                    velocity[0] + offset, 1,
                                    tmp = stdVelocity[i] + offset, 1);
                    Vmath::Smul(nq, gmat[i][0],
                                    soundspeed + offset, 1,
                                    tmp = stdSoundSpeed[i] + offset, 1);
                    for (int j = 1; j < expdim; ++j)
                    {
                        Vmath::Svtvp(nq, gmat[expdim*j+i][0],
                                         velocity[j] + offset, 1,
                                         stdVelocity[i] + offset, 1,
                                         tmp = stdVelocity[i] + offset, 1);
                        Vmath::Svtvp(nq, gmat[expdim*j+i][0],
                                         soundspeed + offset, 1,
                                         stdSoundSpeed[i] + offset, 1,
                                         tmp = stdSoundSpeed[i] + offset, 1);
                    }
                }
            }

            NekDouble vel;
            for (int i = 0; i < nq; ++i)
            {
                NekDouble pntVelocity = 0.0;
                for (int j = 0; j < expdim; ++j)
                {
                    // Add sound speed
                    vel = std::abs(stdVelocity[j][offset + i]) +
                          std::abs(stdSoundSpeed[j][offset + i]);
                    pntVelocity += vel * vel;
                }
                pntVelocity = sqrt(pntVelocity);
                if (pntVelocity > stdV[el])
                {
                    stdV[el] = pntVelocity;
                }
            }
        }

        return stdV;
    }

    /**
     * @brief Set the denominator to compute the time step when a cfl
     * control is employed. This function is no longer used but is still
     * here for being utilised in the future.
     *
     * @param n   Order of expansion element by element.
     */
    NekDouble CompressibleFlowSystem::GetStabilityLimit(int n)
    {
        ASSERTL0(n <= 20, "Illegal modes dimension for CFL calculation "
                          "(P has to be less then 20)");

        NekDouble CFLDG[21] = {  2.0000,   6.0000,  11.8424,  19.1569,
                                27.8419,  37.8247,  49.0518,  61.4815,
                                75.0797,  89.8181, 105.6700, 122.6200,
                               140.6400, 159.7300, 179.8500, 201.0100,
                               223.1800, 246.3600, 270.5300, 295.6900,
                               321.8300}; //CFLDG 1D [0-20]
        NekDouble CFL = 0.0;

        if (m_projectionType == MultiRegions::eDiscontinuous)
        {
            CFL = CFLDG[n];
        }
        else
        {
            ASSERTL0(false, "Continuous Galerkin stability coefficients "
                            "not introduced yet.");
        }

        return CFL;
    }

    /**
     * @brief Compute the vector of denominators to compute the time step
     * when a cfl control is employed. This function is no longer used but
     * is still here for being utilised in the future.
     *
     * @param ExpOrder   Order of expansion element by element.
     */
    Array<OneD, NekDouble> CompressibleFlowSystem::GetStabilityLimitVector(
        const Array<OneD,int> &ExpOrder)
    {
        int i;
        Array<OneD,NekDouble> returnval(m_fields[0]->GetExpSize(), 0.0);
        for (i =0; i<m_fields[0]->GetExpSize(); i++)
        {
            returnval[i] = GetStabilityLimit(ExpOrder[i]);
        }
        return returnval;
    }

    void CompressibleFlowSystem::v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables)
    {
        bool extraFields;
        m_session->MatchSolverInfo("OutputExtraFields","True",
                                   extraFields, true);
        if (extraFields)
        {
            const int nPhys   = m_fields[0]->GetNpoints();
            const int nCoeffs = m_fields[0]->GetNcoeffs();
            Array<OneD, Array<OneD, NekDouble> > tmp(m_fields.size());

            for (int i = 0; i < m_fields.size(); ++i)
            {
                tmp[i] = m_fields[i]->GetPhys();
            }

            Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);
            Array<OneD, Array<OneD, NekDouble> > velFwd  (m_spacedim);
            for (int i = 0; i < m_spacedim; ++i)
            {
                velocity[i] = Array<OneD, NekDouble> (nPhys);
                velFwd[i]   = Array<OneD, NekDouble> (nCoeffs);
            }

            Array<OneD, NekDouble> pressure(nPhys), temperature(nPhys);
            Array<OneD, NekDouble> entropy(nPhys);
            Array<OneD, NekDouble> soundspeed(nPhys), mach(nPhys);
            Array<OneD, NekDouble> sensor(nPhys), SensorKappa(nPhys);

            m_varConv->GetVelocityVector(tmp, velocity);
            m_varConv->GetPressure  (tmp, pressure);
            m_varConv->GetTemperature(tmp, temperature);
            m_varConv->GetEntropy   (tmp, entropy);
            m_varConv->GetSoundSpeed(tmp, soundspeed);
            m_varConv->GetMach      (tmp, soundspeed, mach);

            int sensorOffset;
            m_session->LoadParameter ("SensorOffset", sensorOffset, 1);
            m_varConv->GetSensor (m_fields[0], tmp, sensor, SensorKappa,
                                    sensorOffset);

            Array<OneD, NekDouble> pFwd(nCoeffs), TFwd(nCoeffs);
            Array<OneD, NekDouble> sFwd(nCoeffs);
            Array<OneD, NekDouble> aFwd(nCoeffs), mFwd(nCoeffs);
            Array<OneD, NekDouble> sensFwd(nCoeffs);

            string velNames[3] = {"u", "v", "w"};
            for (int i = 0; i < m_spacedim; ++i)
            {
                m_fields[0]->FwdTrans_IterPerExp(velocity[i], velFwd[i]);
                variables.push_back(velNames[i]);
                fieldcoeffs.push_back(velFwd[i]);
            }

            m_fields[0]->FwdTrans_IterPerExp(pressure,   pFwd);
            m_fields[0]->FwdTrans_IterPerExp(temperature,TFwd);
            m_fields[0]->FwdTrans_IterPerExp(entropy,    sFwd);
            m_fields[0]->FwdTrans_IterPerExp(soundspeed, aFwd);
            m_fields[0]->FwdTrans_IterPerExp(mach,       mFwd);
            m_fields[0]->FwdTrans_IterPerExp(sensor,     sensFwd);

            variables.push_back  ("p");
            variables.push_back  ("T");
            variables.push_back  ("s");
            variables.push_back  ("a");
            variables.push_back  ("Mach");
            variables.push_back  ("Sensor");
            fieldcoeffs.push_back(pFwd);
            fieldcoeffs.push_back(TFwd);
            fieldcoeffs.push_back(sFwd);
            fieldcoeffs.push_back(aFwd);
            fieldcoeffs.push_back(mFwd);
            fieldcoeffs.push_back(sensFwd);

            if (m_artificialDiffusion)
            {
                // Get min h/p
                m_artificialDiffusion->SetElmtHP(GetElmtMinHP());
                // reuse pressure
                Array<OneD, NekDouble> sensorFwd(nCoeffs);
                m_artificialDiffusion->GetArtificialViscosity(tmp, pressure);
                m_fields[0]->FwdTrans_IterPerExp(pressure,   sensorFwd);

                variables.push_back  ("ArtificialVisc");
                fieldcoeffs.push_back(sensorFwd);
            }
        }
    }

    /**
     *
     */
    void CompressibleFlowSystem::GetPressure(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD, NekDouble>                     &pressure)
    {
        m_varConv->GetPressure(physfield, pressure);
    }

    /**
     *
     */
    void CompressibleFlowSystem::GetDensity(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD, NekDouble>                     &density)
    {
        density = physfield[0];
    }

    /**
     *
     */
    void CompressibleFlowSystem::GetVelocity(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> >       &velocity)
    {
        m_varConv->GetVelocityVector(physfield, velocity);
    }

    void CompressibleFlowSystem::v_DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >       &pBwd)
    {
        boost::ignore_unused(inarray, outarray, pFwd, pBwd);
        if (m_shockCaptureType != "Off")
        {
            m_artificialDiffusion->DoArtificialDiffusion(inarray, outarray);
        }
    }

    void CompressibleFlowSystem::v_DoDiffusionCoeff(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >       &pBwd)
    {
        boost::ignore_unused(inarray, outarray, pFwd, pBwd);
        // Do nothing by default
    }

/**
 * @brief Compute an estimate of minimum h/p
 * for each element of the expansion.
 */
Array<OneD, NekDouble>  CompressibleFlowSystem::GetElmtMinHP(void)
{
    int nElements               = m_fields[0]->GetExpSize();
    Array<OneD, NekDouble> hOverP(nElements, 1.0);

    // Determine h/p scaling
    Array<OneD, int> pOrderElmt = m_fields[0]->EvalBasisNumModesMaxPerExp();
    for (int e = 0; e < nElements; e++)
    {
        NekDouble h = 1.0e+10;
        switch(m_expdim)
        {
            case 3:
            {
                LocalRegions::Expansion3DSharedPtr exp3D;
                exp3D = m_fields[0]->GetExp(e)->as<LocalRegions::Expansion3D>();
                for (int i = 0; i < exp3D->GetNtraces(); ++i)
                {
                    h = min(h, exp3D->GetGeom3D()->GetEdge(i)->GetVertex(0)->
                        dist(*(exp3D->GetGeom3D()->GetEdge(i)->GetVertex(1))));
                }
            break;
            }

            case 2:
            {
                LocalRegions::Expansion2DSharedPtr exp2D;
                exp2D = m_fields[0]->GetExp(e)->as<LocalRegions::Expansion2D>();
                for (int i = 0; i < exp2D->GetNtraces(); ++i)
                {
                    h = min(h, exp2D->GetGeom2D()->GetEdge(i)->GetVertex(0)->
                        dist(*(exp2D->GetGeom2D()->GetEdge(i)->GetVertex(1))));
                }
            break;
            }
            case 1:
            {
                LocalRegions::Expansion1DSharedPtr exp1D;
                exp1D = m_fields[0]->GetExp(e)->as<LocalRegions::Expansion1D>();

                h = min(h, exp1D->GetGeom1D()->GetVertex(0)->
                    dist(*(exp1D->GetGeom1D()->GetVertex(1))));

            break;
            }
            default:
            {
                ASSERTL0(false,"Dimension out of bound.")
            }
        }

        // Determine h/p scaling
        hOverP[e] = h/max(pOrderElmt[e]-1,1);

    }
    return hOverP;
}
}
