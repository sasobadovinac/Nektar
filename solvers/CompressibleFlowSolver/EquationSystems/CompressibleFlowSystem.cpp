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
// Description: Compressible flow system base class with auxiliary functions
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>

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

        for (int i = 0; i < m_fields.num_elements(); i++)
        {
            // Use BwdTrans to make sure initial condition is in solution space
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                  m_fields[i]->UpdatePhys());
        }

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

        m_session->LoadParameter("DebugVolTraceSwitch",    
            m_DebugVolTraceSwitch     ,    0);
        m_CalcTracePartFlag = true;
        m_CalcVolumPartFlag = true;
        switch (m_DebugVolTraceSwitch)
        {
        case 1:
            {
                m_CalcTracePartFlag = false;
            }
            break;
        case 2:
            {
                m_CalcVolumPartFlag = false;
            }
            break;
        }

        m_advObject->SetCalcTracePartFlag( m_CalcTracePartFlag);
        m_advObject->SetCalcVolumPartFlag( m_CalcVolumPartFlag);

        // Create artificial diffusion
        if (m_shockCaptureType != "Off")
        {
            if (m_shockCaptureType == "Physical")
            {
                int nPts = m_fields[0]->GetTotPoints();
                m_muav = Array<OneD, NekDouble>(nPts, 0.0);

                int nTracePts = m_fields[0]->GetTrace()->GetTotPoints();
                m_muavTrace = Array<OneD, NekDouble> (nTracePts,0.0);
            }
            else
            {
                m_artificialDiffusion = GetArtificialDiffusionFactory()
                                        .CreateInstance(m_shockCaptureType,
                                                        m_session,
                                                        m_fields,
                                                        m_spacedim);
            }
        }

        // Forcing terms for the sponge region
        m_forcing = SolverUtils::Forcing::Load(m_session, m_fields,
                                               m_fields.num_elements());

        // User-defined boundary conditions
        int cnt = 0;
        for (int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            std::string type =
                m_fields[0]->GetBndConditions()[n]->GetUserDefined();

            if (m_fields[0]->GetBndConditions()[n]->GetBoundaryConditionType()
                == SpatialDomains::ePeriodic)
            {
                continue;
            }

            if(!type.empty())
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
            m_ode.DefineOdeRhs(
                &CompressibleFlowSystem::DoOdeRhs, this);
            m_ode.DefineProjection(
                &CompressibleFlowSystem::DoOdeProjection, this);
        }
        else
        {
#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF

            m_ode.DefineOdeRhs(
                &CompressibleFlowSystem::DoOdeRhs, this);
            m_ode.DefineProjection(
                &CompressibleFlowSystem::DoOdeProjection, this);
            m_ode.DefineImplicitSolve(
                &CompressibleFlowSystem::DoImplicitSolve_phy2coeff, this);

            //TODO: NekLinSysIterative as a member to avoid repeted initialization
            LibUtilities::CommSharedPtr v_Comm  = 
                m_fields[0]->GetComm()->GetRowComm();
            m_linsol = MemoryManager<NekLinSysIterative>::
                AllocateSharedPtr(m_session,v_Comm); 

            m_maxLinItePerNewton = m_linsol->GetMaxLinIte()*m_MaxNonlinIte + 
                m_MaxNonlinIte;
            m_LinSysOprtors.DefineMatrixMultiply(&CompressibleFlowSystem::
                MatrixMultiply_JacobianFree_coeff_dualtimestep, this);
            m_LinSysOprtors.DefinePrecond(&CompressibleFlowSystem::
                preconditioner_BlkSOR_coeff, this);

            if (boost::iequals(m_session->GetSolverInfo("PRECONDITIONER"),
                               "IncompleteLU"))
            {
                // m_LinSysOprtors.DefinePrecond(&CompressibleFlowSystem::preconditioner_BlkILU_coeff, this);
                m_PrecMatStorage    =   eSparse;

                ASSERTL0(false,"IncompleteLU preconditioner not finished yet");

                // DNekSmvBsrMat::SparseStorageSharedPtr sparseStorage =
                //             MemoryManager<DNekSmvBsrMat::StorageType>::
                //                     AllocateSharedPtr(
                //                         brows, bcols, block_size, bcoMat, matStorage );

                // // Create sparse matrix
                // m_smvbsrmatrix = MemoryManager<DNekSmvBsrMat>::
                //                         AllocateSharedPtr( sparseStorage );

                // matBytes = m_smvbsrmatrix->GetMemoryFootprint();

            }
            else 
            {
                int nvariables  =   m_fields.num_elements();
                // if (boost::iequals(m_session->GetSolverInfo("PRECONDITIONERCFS"),
                //                "MultilvlJacobiIte"))
                // {
                    m_LinSysOprtors.DefinePrecond(&CompressibleFlowSystem::
                        preconditioner_MultiLevel_coeff, this);
                    
                // }

                m_PrecMatStorage    =   eDiagonal;
                m_session->LoadParameter("nPadding", m_nPadding, 4);
                
                int ntmp=0;
                m_session->LoadParameter("PrecondMatDataSingle", ntmp, 1);
                m_flagPrecMatVarsSingle             = true;
                if(0==ntmp)
                {
                    m_flagPrecMatVarsSingle = false;
                }
                if(m_DebugNumJacBSOR)
                {
                    m_flagPrecMatVarsSingle = false;
                }

                m_session->LoadParameter("flagPrecondCacheOptmis", ntmp, 1);
                m_flagPrecondCacheOptmis             = true;
                if(0==ntmp)
                {
                    m_flagPrecondCacheOptmis = false;
                }

                // cout << " flagPrecondCacheOptmis= "<<m_flagPrecondCacheOptmis<<endl;

                if (m_flagPrecMatFree)
                {
                    m_linsolBRJ = MemoryManager<NekLinSysIterative>::
                        AllocateSharedPtr(m_session,v_Comm); 
                    LinSysOperators  tmpLinSysOprtors;

                    tmpLinSysOprtors.DefineMatrixMultiply(
                        &CompressibleFlowSystem::PrecBJRDiagMult, this);
                    tmpLinSysOprtors.DefinePrecond(
                        &CompressibleFlowSystem::PrecBJRDiagPrec, this);
                    m_linsolBRJ->setLinSysOperators(tmpLinSysOprtors);
                }
                
                if(m_flagPrecMatVarsSingle)
                {
                    m_PrecMatVarsSingle = 
                        TensorOfArray2D<SNekBlkMatSharedPtr>(nvariables);
                    for(int i = 0; i < nvariables; i++)
                    {
                        m_PrecMatVarsSingle[i] =  
                            Array<OneD, SNekBlkMatSharedPtr> (nvariables);
                    }
                    AllocatePrecondBlkDiag_coeff(m_PrecMatVarsSingle);

                    if(m_flagPrecondCacheOptmis)
                    {
                        int nelmts  = m_fields[0]->GetNumElmts();
                        int nelmtcoef;
                        Array<OneD, unsigned int > nelmtmatdim(nelmts);
                        for(int i = 0; i < nelmts; i++)
                        {
                            nelmtcoef = m_fields[0]->GetExp(i)->GetNcoeffs();
                            nelmtmatdim[i] = nelmtcoef*nvariables;
                        }
                        AllocateNekBlkMatDig(m_PrecMatSingle,
                            nelmtmatdim,nelmtmatdim);
                    }
                }
                else
                {
                    m_PrecMatVars = 
                        TensorOfArray2D<DNekBlkMatSharedPtr>(nvariables);
                    for(int i = 0; i < nvariables; i++)
                    {
                        m_PrecMatVars[i] = 
                            Array<OneD, DNekBlkMatSharedPtr> (nvariables);
                    }
                    AllocatePrecondBlkDiag_coeff(m_PrecMatVars);
                    if(m_flagPrecondCacheOptmis)
                    {
                        int nelmts  = m_fields[0]->GetNumElmts();
                        int nelmtcoef;
                        Array<OneD, unsigned int > nelmtmatdim(nelmts);
                        for(int i = 0; i < nelmts; i++)
                        {
                            nelmtcoef = m_fields[0]->GetExp(i)->GetNcoeffs();
                            nelmtmatdim[i]  =   nelmtcoef*nvariables;
                        }
                        AllocateNekBlkMatDig(m_PrecMat,nelmtmatdim,nelmtmatdim);
                    }
#ifdef CFS_DEBUGMODE
                    if(m_DebugNumJacBSOR)
                    {
                        int nTotCoeff = GetNcoeffs();
                        m_PrecMatVarsOffDiag = 
                            Array<OneD, Array<OneD, NekDouble> >
                            (nvariables*nTotCoeff);
                        for(int m = 0; m < nvariables*nTotCoeff; m++)
                        {
                            m_PrecMatVarsOffDiag[m] = Array<OneD, NekDouble> 
                                (nvariables*nTotCoeff,0.0);
                        }
                    }
#endif
                }
            }
            m_linsol->setLinSysOperators(m_LinSysOprtors);

            int nvariables  =   m_fields.num_elements();
            Array<OneD, Array<OneD, Array<OneD, int > > >   map;
            bool flag;
            const MultiRegions::LocTraceToTraceMapSharedPtr locTraceToTraceMap =
                m_fields[0]->GetlocTraceToTraceMap();
            m_fields[0]->CalcuTracephysToLeftRightExpphysMap(flag,map);
            locTraceToTraceMap->SetTracephysToLeftRightExpphysMap(map);
            locTraceToTraceMap->SetflagTracephysToLeftRightExpphysMap(flag);

            locTraceToTraceMap->CalcuLocTracephysToTraceIDMap(
                m_fields[0]->GetTrace(),m_spacedim);
            for(int i=1;i<nvariables;i++)
            {
                m_fields[i]->GetlocTraceToTraceMap()->
                    SetTracephysToLeftRightExpphysMap(map);
                m_fields[i]->GetlocTraceToTraceMap()->
                    SetflagTracephysToLeftRightExpphysMap(flag);
                m_fields[i]->GetlocTraceToTraceMap()->
                    SetLocTracephysToTraceIDMap(
                    locTraceToTraceMap->GetLocTracephysToTraceIDMap());
            }

            m_session->LoadParameter("LiniearizationMethod", 
                m_LiniearizationMethod,    0); //0: Fwd Jacobian-Free; 1: Central Jacobian-free; 2: Matrix-free.
            
            if(2==m_LiniearizationMethod||3==m_LiniearizationMethod)
            {
                int ntotElmt = (*m_fields[0]->GetExp()).size();
                m_ElmtFluxJacArray = 
                    TensorOfArray2D<DNekBlkMatSharedPtr> (ntotElmt);
        
                for(int ne=0; ne<ntotElmt;ne++)
                {
                    int nElmtPnt    = GetTotPoints(ne);
                    Array<OneD, unsigned int> n_blks1(nElmtPnt,nvariables);

                    m_ElmtFluxJacArray[ne]  =   
                        Array<OneD, DNekBlkMatSharedPtr >(m_spacedim);
                    for(int ndir=0; ndir<m_spacedim;ndir++)
                    {
                        m_ElmtFluxJacArray[ne][ndir] = MemoryManager<DNekBlkMat>
                            ::AllocateSharedPtr(n_blks1, n_blks1, eDIAGONAL);
                        for(int i=0;i<nElmtPnt;i++)
                        {
                            DNekMatSharedPtr tmpMat = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nvariables, nvariables);
                            m_ElmtFluxJacArray[ne][ndir]->SetBlock(i,i,tmpMat);
                        }
                    }
                }
            }
            
            m_qfield  =   TensorOfArray3D<NekDouble>(m_spacedim);
            int npoints = GetTotPoints();
            for(int i = 0; i< m_spacedim; i++)
            {
                m_qfield[i]   =   
                    Array<OneD, Array<OneD, NekDouble> >(nvariables);
                for(int j = 0; j< nvariables; j++)
                {
                    m_qfield[i][j] = Array<OneD, NekDouble>(npoints,0.0);
                }
            }

            int nTracePts  = GetTraceTotPoints();
            m_MatrixFreeRefFields  =   
                Array<OneD, Array<OneD, NekDouble> > (nvariables);
            m_MatrixFreeRefFwd     =   
                Array<OneD, Array<OneD, NekDouble> > (nvariables);
            m_MatrixFreeRefBwd     =  
                Array<OneD, Array<OneD, NekDouble> > (nvariables);
            for(int j = 0; j< nvariables; j++)
            {
                m_MatrixFreeRefFields[j] = Array<OneD, NekDouble>(npoints);
                m_MatrixFreeRefFwd[j]    = Array<OneD, NekDouble>(nTracePts, 0.0);
                m_MatrixFreeRefBwd[j]    = Array<OneD, NekDouble>(nTracePts, 0.0);
            }
#else
            ASSERTL0(false, "Implicit CFS not set up.");
#endif
        }

        SetBoundaryConditionsBwdWeight();

        string advName;
        m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
        // m_session->LoadSolverInfo("useUnifiedWeakIntegration", m_useUnifiedWeakIntegration, false);
        // m_session->LoadParameter("useUnifiedWeakIntegration", m_useUnifiedWeakIntegration, false);
        m_session->MatchSolverInfo("useUnifiedWeakIntegration", "True", 
            m_useUnifiedWeakIntegration, false);
        if(m_useUnifiedWeakIntegration)
        {
            if(advName=="WeakDG" && m_shockCaptureType!="Smooth"
                &&(eNotHomogeneous == m_HomogeneousType))
            {

            }
            else
            {
                m_useUnifiedWeakIntegration=false;
                if(m_session->DefinesCmdLineArgument("verbose"))
                {
                    WARNINGL0(false, "useUnifiedWeakIntegration not coded");
                }
            }
        }
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
        if(m_useFiltering)
        {
            m_session->LoadParameter ("FilterAlpha", m_filterAlpha, 36);
            m_session->LoadParameter ("FilterExponent", m_filterExponent, 16);
            m_session->LoadParameter ("FilterCutoff", m_filterCutoff, 0);
        }

        // Load CFL for local time-stepping (for steady state)
        m_session->MatchSolverInfo("LocalTimeStep","True",
                                   m_useLocalTimeStep, false);
        if(m_useLocalTimeStep)
        {
            ASSERTL0(m_cflSafetyFactor != 0,
                    "Local time stepping requires CFL parameter.");
        }

        m_session->LoadParameter ("JFEps", m_JFEps, 5.0E-8);

        int ntmp;
        m_session->LoadParameter("DEBUG_ADVECTION_JAC_MAT", ntmp, 1);
        m_DEBUG_ADVECTION_JAC_MAT             = true;
        if(0==ntmp)
        {
            m_DEBUG_ADVECTION_JAC_MAT = false;
        }
        m_session->LoadParameter("DEBUG_VISCOUS_JAC_MAT", ntmp, 1);
        m_DEBUG_VISCOUS_JAC_MAT             = true;
        if(0==ntmp)
        {
            m_DEBUG_VISCOUS_JAC_MAT = false;
        }
        m_session->LoadParameter("DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT", ntmp, 0);
        m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT = false;
        if(1==ntmp)
        {
            m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT = true;
        }

#ifdef CFS_DEBUGMODE
        m_session->LoadParameter("DebugAdvDiffSwitch",     
            m_DebugAdvDiffSwitch      ,    0);
        m_session->LoadParameter("DebugConsDerivSwitch",   
            m_DebugConsDerivSwitch    ,    0);
        m_session->LoadParameter("DebugNumJacMatSwitch",   
            m_DebugNumJacMatSwitch    ,    0);
        m_session->LoadParameter("DebugOutputJacMatSwitch",
            m_DebugOutputJacMatSwitch ,    0);
        m_session->LoadParameter("DebugInvMassSwitch",     
            m_DebugInvMassSwitch      ,    1);
        m_session->LoadParameter("DebugPlusSourceSwitch",  
            m_DebugPlusSourceSwitch   ,    1);

        m_session->LoadParameter("DebugIPSymmFluxJacSwitch",
            m_DebugIPSymmFluxJacSwitch,    0);

        m_session->LoadParameter("DebugNumJacBSOR",
            m_DebugNumJacBSOR,    0);

        m_session->LoadParameter("flagImplItsStatistcs", ntmp, 0);
        m_flagImplItsStatistcs = false;
        if(1==ntmp)
        {
            m_flagImplItsStatistcs = true;
        }
        m_session->LoadParameter("centralDiffTracJac", ntmp , 0);
        m_centralDiffTracJac = false;
        if(1==ntmp)
        {
            m_centralDiffTracJac = true;
        }
        m_session->LoadParameter("flagMultiLvlJacMultiplyMatFree", ntmp , 0);
        m_flagMultiLvlJacMultiplyMatFree = false;
        if(1==ntmp)
        {
            m_flagMultiLvlJacMultiplyMatFree = true;
        }
#endif
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

        m_advObject->SetFluxVectorMF  (&CompressibleFlowSystem::
                                          GetFluxVectorMF, this);
        m_advObject->SetFluxVectorTraceMF (&CompressibleFlowSystem::
                                          GetFluxVectorTraceMF, this);

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

        // Concluding initialisation of advection / diffusion operators
        m_advObject->SetRiemannSolver   (riemannSolver);
        m_advObject->InitObject         (m_session, m_fields);
    }

    /**
     * @brief Compute the right-hand side.
     */
    void CompressibleFlowSystem::DoOdeRhs(
        const TensorOfArray2D<NekDouble>        &inarray,
        Array<OneD, Array<OneD, NekDouble> >    &outarray,
        const NekDouble                         time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();
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
            for(int i = 0; i < nvariables; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
        }

        //Only test solver use the reduced coddes
        if(m_useUnifiedWeakIntegration)
        {
            int nDim        = m_spacedim;
            int nCoeffs     = GetNcoeffs();
            TensorOfArray3D<NekDouble> VolumeFlux1(nvariables);
            TensorOfArray3D<NekDouble> VolumeFlux2(nDim);
            Array<OneD, Array<OneD, NekDouble>> TraceFlux1(nvariables);
            Array<OneD, Array<OneD, NekDouble>> TraceFlux2(nvariables);

            for (int i = 0; i < nvariables; ++i)
            {
                VolumeFlux1[i] = Array<OneD, Array<OneD, NekDouble>>(nDim);
                for (int j= 0; j < nDim; ++j)
                {
                    VolumeFlux1[i][j] = Array<OneD, NekDouble>(npoints,0.0);
                }
            }
            for (int j = 0; j < nDim; ++j)
            {
                VolumeFlux2[j] = 
                    Array<OneD, Array<OneD, NekDouble>>(nvariables);
                for (int i = 0; i < nvariables; ++i)
                {
                    VolumeFlux2[j][i] = Array<OneD, NekDouble>(npoints,0.0);
                }
            }
            for (int i = 0; i < nvariables; ++i)
            {
                TraceFlux1[i]  = Array<OneD, NekDouble>(nTracePts, 0.0);
                TraceFlux2[i]  = Array<OneD, NekDouble>(nTracePts, 0.0);
            }

            Array<OneD, Array<OneD, NekDouble>> advVel(m_spacedim);
            m_advObject->AdvectVolumeFlux(nvariables, m_fields, advVel, 
                inarray, VolumeFlux1, time);
            m_advObject->AdvectTraceFlux(nvariables, m_fields, advVel, 
                inarray,TraceFlux1, time, Fwd, Bwd);
            v_DoDiffusionFlux(inarray, VolumeFlux2, TraceFlux2, Fwd, Bwd);

            // Add Trace and Volume Integral together
            Array<OneD, NekDouble> tmp(nCoeffs, 0.0);
            for (int i = 0; i < nvariables; ++i)
            {
                // Add Volume integral
                for (int j = 0; j < nDim; ++j)
                {
                    // Advection term needs to be negative
                    // Add Advection and Diffusion part
                    Vmath::Vsub(npoints, &VolumeFlux2[j][i][0], 1,
                        &VolumeFlux1[i][j][0], 1, &VolumeFlux1[i][j][0], 1);
                }
                m_fields[i]->IProductWRTDerivBase(VolumeFlux1[i], tmp);

                Vmath::Neg(nCoeffs, &tmp[0], 1);
                // Add Trace integral
                // Advection term needs to be negative
                // Add Advection and Diffusion part
                Vmath::Vsub(nTracePts, &TraceFlux2[i][0], 1, 
                            &TraceFlux1[i][0], 1,
                            &TraceFlux1[i][0], 1);
                m_fields[i]->AddTraceIntegral(TraceFlux1[i], tmp);
                m_fields[i]->MultiplyByElmtInvMass(tmp, tmp);
        
                m_fields[i]->BwdTrans(tmp, outarray[i]);
            }

            // AddDiffusionSymmFluxToPhys(inarray, VolumeFlux2, outarray, Fwd, Bwd);
        }
        else
        {
            //Oringinal CompressibleFlowSolver
            // Calculate advection
            DoAdvection(inarray, outarray, time, Fwd, Bwd);

            // Negate results
            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Neg(npoints, outarray[i], 1);
            }

            // Add diffusion terms
            DoDiffusion(inarray, outarray, Fwd, Bwd);

        }

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
            for(int n = 0; n < nElements; ++n)
            {
                nq     = m_fields[0]->GetExp(n)->GetTotPoints();
                offset = m_fields[0]->GetPhys_Offset(n);
                fac    = tstep[n] / m_timestep;
                for(int i = 0; i < nvariables; ++i)
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
        const TensorOfArray2D<NekDouble>        &inarray,
        Array<OneD, Array<OneD, NekDouble> >    &outarray,
        const NekDouble                         time)
    {
        int i;
        int nvariables = inarray.num_elements();

        switch(m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                // Just copy over array
                int npoints = GetNpoints();

                for(i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
                    if(m_useFiltering)
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

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF

    void CompressibleFlowSystem::preconditioner(
        const Array<OneD, NekDouble>    &inarray,
        Array<OneD, NekDouble >         &out)
    {
        int ntotal     = inarray.num_elements();
        Vmath::Vcopy(ntotal,inarray,1,out,1);
        return;
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::preconditioner_BlkDiag(
        const Array<OneD, NekDouble>                  &inarray,
        Array<OneD, NekDouble >                       &outarray,
        const TensorOfArray2D<TypeNekBlkMatSharedPtr> &PrecMatVars,
        const DataType                                &tmpDataType)
    {
        unsigned int nvariables = m_fields.num_elements();
        unsigned int npoints    = GetNcoeffs();
        Array<OneD, Array<OneD, DataType > >Sinarray(nvariables);
        Array<OneD, NekVector<DataType> >tmpVect(nvariables);
        Array<OneD, DataType > Soutarray(npoints);
        NekVector<DataType> outVect(npoints,Soutarray,eWrapper);

        for(int m = 0; m < nvariables; m++)
        {
            int moffset = m*npoints;
            Sinarray[m] = Array<OneD, DataType > (npoints);
            for(int i=0;i<npoints;i++)
            {
                Sinarray[m][i]  =  DataType(inarray[moffset+i]);
            }
            tmpVect[m] =  NekVector<DataType> (npoints,Sinarray[m],eWrapper);
        }

        Vmath::Fill(outarray.num_elements(),0.0,outarray,1);

        for(int m = 0; m < nvariables; m++)
        {
            Vmath::Zero(npoints,Soutarray,1);
            int moffset = m*npoints;
            for(int n = 0; n < nvariables; n++)
            {
                outVect += (*PrecMatVars[m][n])*tmpVect[n];
            }

            for(int i=0;i<npoints;i++)
            {
                outarray[moffset+i]  =  NekDouble(Soutarray[i]);
            }
        }
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::preconditioner_BlkDiag(
        const Array<OneD, NekDouble>                        &inarray,
        Array<OneD, NekDouble >                             &outarray,
        const TypeNekBlkMatSharedPtr                        &PrecMatVars,
        const DataType                                      &tmpDataType)
    {
        unsigned int nvariables = m_fields.num_elements();
        unsigned int npoints    = GetNcoeffs();
        unsigned int npointsVar = nvariables*npoints;
        Array<OneD, DataType >Sinarray(npointsVar);
        Array<OneD, DataType > Soutarray(npointsVar);
        NekVector<DataType> tmpVect(npointsVar,Sinarray,eWrapper);
        NekVector<DataType> outVect(npointsVar,Soutarray,eWrapper);

        std::shared_ptr<LocalRegions::ExpansionVector> expvect =    
            m_fields[0]->GetExp();
        int ntotElmt            = (*expvect).size();

        for(int m = 0; m < nvariables; m++)
        {
            int nVarOffset = m*npoints;
            for(int ne=0;ne<ntotElmt;ne++)
            {
                int nCoefOffset = GetCoeff_Offset(ne);
                int nElmtCoef   = GetNcoeffs(ne);
                int inOffset    = nVarOffset+nCoefOffset;
                int outOffset   = nCoefOffset*nvariables+m*nElmtCoef;
                for(int i=0;i<nElmtCoef;i++)
                {
                    Sinarray[outOffset+i]  =  DataType(inarray[inOffset+i]);
                }
            }
        }

        outVect = (*PrecMatVars)*tmpVect;

        for(int m = 0; m < nvariables; m++)
        {
            int nVarOffset = m*npoints;
            for(int ne=0;ne<ntotElmt;ne++)
            {
                int nCoefOffset = GetCoeff_Offset(ne);
                int nElmtCoef   = GetNcoeffs(ne);
                int inOffset    = nVarOffset+nCoefOffset;
                int outOffset   = nCoefOffset*nvariables+m*nElmtCoef;
                for(int i=0;i<nElmtCoef;i++)
                {
                    outarray[inOffset+i]  =  NekDouble(Soutarray[outOffset+i]);
                }
            }
        }
    }

    template<typename DataType>
    void CompressibleFlowSystem::preconditioner_BlkDiagMatFree(
        const TensorOfArray1D<NekDouble>     &inarray,
        TensorOfArray1D<NekDouble>           &outarray,
        const DataType                       &tmpDataType)
    {
        unsigned int ntotal = inarray.num_elements();
        LibUtilities::CommSharedPtr v_Comm = m_fields[0]->
                                            GetComm()->GetRowComm();

        NekDouble resnorm = Vmath::Dot(ntotal,inarray,inarray);
        v_Comm->AllReduce(resnorm, Nektar::LibUtilities::ReduceSum);
        NekDouble BRJDiagInvTol = 1.0E-5;
        NekDouble LinSysTol = BRJDiagInvTol*sqrt(resnorm);

        int ntmpGMRESIts =  m_linsolBRJ->SolveLinearSystem(ntotal, 
            inarray,outarray,0,LinSysTol);
    }

    void CompressibleFlowSystem::PrecBJRDiagMult(
        const  TensorOfArray1D<NekDouble>   &inarray,
        TensorOfArray1D<NekDouble>          &outarray,
        const  bool                         &controlFlag)
    {
        size_t nvariables = m_fields.num_elements();
        size_t npoints    = GetNpoints();
        size_t ncoeffs    = GetNcoeffs();
        
        TensorOfArray2D<NekDouble> inpnts(nvariables);
        TensorOfArray2D<NekDouble> out_2d(nvariables);
        for (int m = 0; m < nvariables; ++m)
        {
            int moffset = m*ncoeffs;
            out_2d[m]   = outarray  + moffset;
            inpnts[m]   = TensorOfArray1D<NekDouble> {npoints};

            TensorOfArray1D<NekDouble> tmp = inarray + moffset;
            m_fields[m]->BwdTrans(tmp, inpnts[m]);
        }

        DoOdeRhs_coeffVol(inpnts, out_2d, true);
        // AddTraceIntglToRhs(inarray, outarray);
    }

    void CompressibleFlowSystem::PrecBJRDiagPrec(
        const TensorOfArray1D<NekDouble>    &inarray,
        TensorOfArray1D<NekDouble>          &outarray,
        const bool                          &flag)
    {
        Vmath::Vcopy(inarray.num_elements(), inarray, 1, outarray, 1);
    }

    void CompressibleFlowSystem::preconditioner_NumJac(
        const Array<OneD, NekDouble>                 &inarray,
        Array<OneD, NekDouble >                      &outarray,
        const TensorOfArray2D<DNekBlkMatSharedPtr>   &PrecMatVars,
        const Array<OneD, Array<OneD, NekDouble > >  &PrecMatVarsOffDiag)
    {
        const NekDouble SORParam        =   m_SORRelaxParam;
        const NekDouble OmSORParam      =   1.0-SORParam;

        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int npoints    = m_TimeIntegtSol_n[0].num_elements();
        unsigned int ntotpnt    = inarray.num_elements();
        
        ASSERTL0(nvariables*npoints==ntotpnt,
            "nvariables*npoints!=ntotpnt in preconditioner_BlkSOR");

        Array<OneD, NekDouble> rhs(ntotpnt);

        Array<OneD, NekDouble>  outN(ntotpnt);
        Array<OneD, NekDouble>  outTmp(ntotpnt);
        Vmath::Vcopy(ntotpnt,&inarray[0],1,&rhs[0],1);

        NekDouble tmpDouble = 0.0;
        preconditioner_BlkDiag(rhs,outarray,PrecMatVars,tmpDouble);

        int nSORTot   =   m_JFNKPrecondStep;
        for(int nsor = 0; nsor < nSORTot-1; nsor++)
        {
            Vmath::Smul(ntotpnt,OmSORParam,outarray,1,outN,1);
            
            MinusOffDiag2RhsNumJac(nvariables,npoints,rhs,outarray,
                PrecMatVarsOffDiag);

            preconditioner_BlkDiag(outarray,outTmp,PrecMatVars,tmpDouble);
            Vmath::Svtvp(ntotpnt,SORParam,outTmp,1,outN,1,outarray,1);
        }
    }

    void CompressibleFlowSystem::MinusOffDiag2RhsNumJac(
        const int                                   nvariables,
        const int                                   nCoeffs,
        const Array<OneD, NekDouble>                &rhs,
        Array<OneD, NekDouble>                      &outarray,
        const Array<OneD, Array<OneD, NekDouble > > &PrecMatVarsOffDiag)
    {
        int nTotCoef = nvariables*nCoeffs;
        Array<OneD, NekDouble> tmp (nTotCoef,0.0);
        for(int i=0;i<nTotCoef;i++)
        {
            Vmath::Svtvp(nTotCoef,outarray[i],&PrecMatVarsOffDiag[i][0],1,
                &tmp[0],1,&tmp[0],1);
        }

        // for(int i=0;i<nTotCoef;i++)
        // {
        //     for(int j=0;j<nTotCoef;j++)
        //     {
        //         tmp[j] += outarray[i]*PrecMatVarsOffDiag[j][i];
        //     }
        // }

        Vmath::Vsub(nTotCoef,&rhs[0],1,&tmp[0],1,&outarray[0],1);
    }

    void CompressibleFlowSystem::preconditioner_BlkSOR_coeff(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&outarray,
            const bool                   &Continueflag)
    {
        int nSORTot   =   m_JFNKPrecondStep;
#ifdef CFS_DEBUGMODE
        if(m_DebugNumJacBSOR)
        {

            preconditioner_NumJac(inarray,outarray,m_PrecMatVars,
                m_PrecMatVarsOffDiag);
        }
        else
        {
#endif
            if (0==nSORTot)
            {
                preconditioner(inarray,outarray);
            }
            else
            {
                const NekDouble SORParam        =   m_SORRelaxParam;
                const NekDouble OmSORParam      =   1.0-SORParam;

                unsigned int nvariables = m_fields.num_elements();
                unsigned int npoints    = GetNcoeffs();
                unsigned int ntotpnt    = inarray.num_elements();
                
                ASSERTL0(nvariables*npoints==ntotpnt,
                    "nvariables*npoints!=ntotpnt in preconditioner_BlkSOR");


                Array<OneD, NekDouble> rhs(ntotpnt);
                // Vmath::Vcopy(ntotpnt,inarray,1,rhs,1);

                PointerWrapper pwrapp = eWrapper;

                Array<OneD, NekDouble>  outN(ntotpnt);
                Array<OneD, NekDouble>  outTmp(ntotpnt);
                Array<OneD, Array<OneD, NekDouble> >rhs2d(nvariables);
                Array<OneD, Array<OneD, NekDouble> >out_2d(nvariables);
                Array<OneD, Array<OneD, NekDouble> >outTmp_2d(nvariables);
                for(int m = 0; m < nvariables; m++)
                {
                    int moffset     = m*npoints;
                    rhs2d[m]        = rhs       + moffset;
                    out_2d[m]       = outarray  + moffset;
                    outTmp_2d[m]    = outTmp    + moffset;
                    m_fields[m]->MultiplyByMassMatrix(inarray+moffset,rhs2d[m]);
                }

                int nphysic    = GetNpoints();
                int nTracePts  = GetTraceTotPoints();
                TensorOfArray3D<NekDouble> qfield(m_spacedim);
                for(int i = 0; i< m_spacedim; i++)
                {
                    qfield[i]   =   
                        Array<OneD, Array<OneD, NekDouble> >(nvariables);
                    for(int j = 0; j< nvariables; j++)
                    {
                        qfield[i][j]   =   Array<OneD, NekDouble>(nphysic);
                    }
                }
                int ntmpTrace = 4+2*m_spacedim;
                TensorOfArray3D<NekDouble> tmpTrace(ntmpTrace);
                for(int i = 0; i< ntmpTrace; i++)
                {
                    tmpTrace[i]   =   
                        Array<OneD, Array<OneD, NekDouble> >(nvariables);
                    for(int j = 0; j< nvariables; j++)
                    {
                        tmpTrace[i][j]   =   Array<OneD, NekDouble>(nTracePts);
                    }
                }
                Array<OneD, Array<OneD, NekDouble> > FwdFluxDeriv(nvariables);
                Array<OneD, Array<OneD, NekDouble> > BwdFluxDeriv(nvariables);
                for(int j = 0; j< nvariables; j++)
                {
                    FwdFluxDeriv[j]   =   Array<OneD, NekDouble>(nTracePts);
                    BwdFluxDeriv[j]   =   Array<OneD, NekDouble>(nTracePts);
                }

                bool flagUpdateDervFlux = false;
                if(m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT)
                {
                    flagUpdateDervFlux = true;
                }

                const int nwspTraceDataType = nvariables+1;
                if(m_flagPrecMatVarsSingle)
                {
                    NekSingle tmpSingle;
                    if(m_flagPrecondCacheOptmis)
                    {
                        Array<OneD, Array<OneD, NekSingle> > 
                            wspTraceDataType(nwspTraceDataType);
                        for(int m=0;m<nwspTraceDataType;m++)
                        {
                            wspTraceDataType[m] =
                                Array<OneD, NekSingle>(nTracePts);
                        }

                        if (m_flagPrecMatFree)
                        {
                            preconditioner_BlkDiagMatFree(rhs,outarray,
                                tmpSingle);

                            for(int nsor = 0; nsor < nSORTot-1; nsor++)
                            {
                                Vmath::Smul(ntotpnt,OmSORParam,
                                    outarray,1,outN,1);
                                
                                AddOrMinusPrecTraceFlux(
                                    ePrecDiagOffDiagTagOffdiag, -1, 
                                    rhs2d, out_2d, flagUpdateDervFlux, 
                                    FwdFluxDeriv, BwdFluxDeriv, qfield, 
                                    tmpTrace, wspTraceDataType, 
                                    m_TraceJacArraySingle, 
                                    m_TraceJacDerivArraySingle, 
                                    m_TraceJacDerivSignSingle,
                                    m_TraceIPSymJacArraySingle);

                                // Vmath::Zero(ntotpnt,outTmp,1);
                                preconditioner_BlkDiagMatFree(outarray,outTmp,
                                    tmpSingle);
                                Vmath::Svtvp(ntotpnt,SORParam,outTmp,1,outN,1,
                                    outarray,1);
                                if(m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT)
                                {
                                    flagUpdateDervFlux = true;
                                }
                            }
                        }
                        else
                        {
                            if(Continueflag)
                            {
                                nSORTot=nSORTot+1;
                            }
                            else
                            {
                                preconditioner_BlkDiag(rhs,outarray,
                                    m_PrecMatSingle,tmpSingle);
                            }
                        
                            for(int nsor = 0; nsor < nSORTot-1; nsor++)
                            {
                                Vmath::Smul(ntotpnt,OmSORParam,outarray,1,outN,1);
                                
                                AddOrMinusPrecTraceFlux(
                                    ePrecDiagOffDiagTagOffdiag, -1, 
                                    rhs2d, out_2d, flagUpdateDervFlux, 
                                    FwdFluxDeriv, BwdFluxDeriv, qfield, 
                                    tmpTrace, wspTraceDataType, 
                                    m_TraceJacArraySingle, 
                                    m_TraceJacDerivArraySingle, 
                                    m_TraceJacDerivSignSingle,
                                    m_TraceIPSymJacArraySingle);

                                // Vmath::Zero(ntotpnt,outTmp,1);
                                preconditioner_BlkDiag(outarray,outTmp,
                                    m_PrecMatSingle,tmpSingle);
                                Vmath::Svtvp(ntotpnt,SORParam,outTmp,1,outN,1,
                                    outarray,1);
                                if(m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT)
                                {
                                    flagUpdateDervFlux = true;
                                }
                            }
                        }
                    }
                    else
                    {
                        if(Continueflag)
                        {
                            nSORTot=nSORTot+1;
                        }
                        else
                        {
                            preconditioner_BlkDiag(rhs,outarray,
                                m_PrecMatVarsSingle,tmpSingle);
                        }

                        for(int nsor = 0; nsor < nSORTot-1; nsor++)
                        {
                            Vmath::Smul(ntotpnt,OmSORParam,outarray,1,outN,1);
                            
                            MinusOffDiag2Rhs(nvariables,npoints,rhs2d,out_2d,
                                flagUpdateDervFlux,FwdFluxDeriv,BwdFluxDeriv,
                                qfield,tmpTrace,m_TraceJacSingle, 
                                m_TraceJacDerivSingle, 
                                m_TraceJacDerivSignSingle);

                            // Vmath::Zero(ntotpnt,outTmp,1);
                            preconditioner_BlkDiag(outarray,outTmp,
                                m_PrecMatVarsSingle,tmpSingle);
                            Vmath::Svtvp(ntotpnt,SORParam,outTmp,1,outN,1,
                                outarray,1);
                            if(m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT)
                            {
                                flagUpdateDervFlux = true;
                            }
                        }
                    }
                }
                else
                {
                    NekDouble tmpDouble;
                    if(m_flagPrecondCacheOptmis)
                    {
                        Array<OneD, Array<OneD, NekDouble> > 
                            wspTraceDataType(nwspTraceDataType);
                        for(int m=0;m<nwspTraceDataType;m++)
                        {
                            wspTraceDataType[m] =   
                                Array<OneD, NekDouble>(nTracePts);
                        }
                        
                        if(Continueflag)
                        {
                            nSORTot=nSORTot+1;
                        }
                        else
                        {
                             preconditioner_BlkDiag(rhs,outarray,m_PrecMat,
                                tmpDouble);
                        }

                        for(int nsor = 0; nsor < nSORTot-1; nsor++)
                        {
                            Vmath::Smul(ntotpnt,OmSORParam,outarray,1,outN,1);
                            
                            AddOrMinusPrecTraceFlux(
                                ePrecDiagOffDiagTagOffdiag, -1, 
                                rhs2d, out_2d, flagUpdateDervFlux, 
                                FwdFluxDeriv, BwdFluxDeriv, qfield, 
                                tmpTrace, wspTraceDataType, 
                                m_TraceJacArray, m_TraceJacDerivArray, 
                                m_TraceJacDerivSign,m_TraceIPSymJacArray);

                            // Vmath::Zero(ntotpnt,outTmp,1);
                            preconditioner_BlkDiag(outarray,outTmp,m_PrecMat,
                                tmpDouble);
                            Vmath::Svtvp(ntotpnt,SORParam,outTmp,1,outN,1,
                                outarray,1);
                            if(m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT)
                            {
                                flagUpdateDervFlux = true;
                            }
                        }
                    }
                    else
                    {
                        if(Continueflag)
                        {
                            nSORTot=nSORTot+1;
                        }
                        else
                        {
                            preconditioner_BlkDiag(rhs,outarray,m_PrecMatVars,
                                tmpDouble);
                        }

                        for(int nsor = 0; nsor < nSORTot-1; nsor++)
                        {
                            Vmath::Smul(ntotpnt,OmSORParam,outarray,1,outN,1);
                            
                            MinusOffDiag2Rhs(nvariables,npoints,rhs2d,out_2d,
                                flagUpdateDervFlux,FwdFluxDeriv,BwdFluxDeriv,
                                qfield,tmpTrace,m_TraceJac, m_TraceJacDeriv, 
                                m_TraceJacDerivSign);

                            // Vmath::Zero(ntotpnt,outTmp,1);
                            preconditioner_BlkDiag(outarray,outTmp,
                                m_PrecMatVars,tmpDouble);
                            Vmath::Svtvp(ntotpnt,SORParam,outTmp,1,outN,1,
                                outarray,1);
                            if(m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT)
                            {
                                flagUpdateDervFlux = true;
                            }
                        }
                    }
                }
            }
#ifdef CFS_DEBUGMODE
        }
#endif
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::MinusOffDiag2Rhs(
            const int                                  nvariables,
            const int                                  nCoeffs,
            const TensorOfArray2D<NekDouble>           &inarray,
            Array<OneD, Array<OneD, NekDouble> >       &outarray,
            bool                                       flagUpdateDervFlux,
            Array<OneD, Array<OneD, NekDouble> >       &FwdFluxDeriv,
            Array<OneD, Array<OneD, NekDouble> >       &BwdFluxDeriv,
            TensorOfArray3D<NekDouble>                 &qfield,
            TensorOfArray3D<NekDouble>                 &tmpTrace,
            const Array<OneD, TypeNekBlkMatSharedPtr > &TraceJac,
            const Array<OneD, TypeNekBlkMatSharedPtr > &TraceJacDeriv,
            const Array<OneD, Array<OneD, DataType> >  &TraceJacDerivSign)
    {
        int nTracePts  = GetTraceTotPoints();
        int npoints    = GetNpoints();
        int nDim       = m_spacedim;
        int nConvectiveFields = nvariables;

        Array<OneD, Array<OneD, NekDouble> > outpnts(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            outpnts[i]  =  Array<OneD, NekDouble> (npoints,0.0);
            m_fields[i]->BwdTrans(outarray[i],outpnts[i]);
        }

        for(int i = 0; i< m_spacedim; i++)
        {
            for(int j = 0; j< nConvectiveFields; j++)
            {
                Vmath::Zero(npoints,qfield[i][j],1);
            }
        }
        for(int i = 0; i< tmpTrace.num_elements(); i++)
        {
            for(int j = 0; j< nConvectiveFields; j++)
            {
                Vmath::Zero(nTracePts,tmpTrace[i][j],1);
            }
        }

        // Store forwards/backwards space along trace space
        Array<OneD, Array<OneD, NekDouble> > Fwd    ;
        Array<OneD, Array<OneD, NekDouble> > Bwd    ;
        Array<OneD, Array<OneD, NekDouble> > FwdFlux;
        Array<OneD, Array<OneD, NekDouble> > BwdFlux;
        int indextmpTrace = 0;
        Fwd     =   tmpTrace[indextmpTrace], indextmpTrace++;
        Bwd     =   tmpTrace[indextmpTrace], indextmpTrace++;
        FwdFlux =   tmpTrace[indextmpTrace], indextmpTrace++;
        BwdFlux =   tmpTrace[indextmpTrace], indextmpTrace++;

        for(int i = 0; i < nvariables; ++i)
        {
            m_fields[i]->GetFwdBwdTracePhys(outpnts[i], Fwd[i], Bwd[i]);
        }

#ifdef CFS_DEBUGMODE
        if(1!=m_DebugVolTraceSwitch)
        {
#endif
        int nTracePtsVars = nTracePts*nvariables;
        Array<OneD, DataType> Fwdarray (nTracePtsVars);
        Array<OneD, DataType> Fwdreslt (nTracePtsVars);
        NekVector<DataType> VFwdarray(nTracePtsVars,Fwdarray,eWrapper);
        NekVector<DataType> VFwdreslt(nTracePtsVars,Fwdreslt,eWrapper);

        for(int n = 0; n < nTracePts; ++n)
        {
            int noffset = n*nvariables; 
            for(int i = 0; i < nvariables; ++i)
            {
                Fwdarray[noffset+i] =  DataType( Fwd[i][n] );
            }
        }
        VFwdreslt = (*TraceJac[0])*VFwdarray;
        for(int n = 0; n < nTracePts; ++n)
        {
            int noffset = n*nvariables; 
            for(int i = 0; i < nvariables; ++i)
            {
                FwdFlux[i][n] = NekDouble( Fwdreslt[noffset+i] );
            }
        }

        for(int n = 0; n < nTracePts; ++n)
        {
            int noffset = n*nvariables; 
            for(int i = 0; i < nvariables; ++i)
            {
                Fwdarray[noffset+i] =  DataType( Bwd[i][n] );
            }
        }
        VFwdreslt = (*TraceJac[1])*VFwdarray;
        for(int n = 0; n < nTracePts; ++n)
        {
            int noffset = n*nvariables; 
            for(int i = 0; i < nvariables; ++i)
            {
                BwdFlux[i][n] = NekDouble( Fwdreslt[noffset+i] );
            }
        }

        if(m_DEBUG_VISCOUS_JAC_MAT&&m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT)
        {
            const MultiRegions::AssemblyMapDGSharedPtr TraceMap =
                m_fields[0]->GetTraceMap();
            
            if(flagUpdateDervFlux)
            {
                CalphysDeriv(outpnts,qfield);

                TensorOfArray3D<NekDouble>    numDerivBwd(nDim);
                TensorOfArray3D<NekDouble>    numDerivFwd(nDim);
                for (int nd = 0; nd < nDim; ++nd)
                {
                    numDerivBwd[nd] = tmpTrace[indextmpTrace], indextmpTrace++;
                    numDerivFwd[nd] = tmpTrace[indextmpTrace], indextmpTrace++;
                }

                for (int nd = 0; nd < nDim; ++nd)
                {
                    for (int i = 0; i < nConvectiveFields; ++i)
                    {
                        Vmath::Zero(nTracePts, Bwd[i],1);
                        Vmath::Zero(nTracePts, Fwd[i],1);
                        m_fields[i]->GetFwdBwdTracePhysNoBndFill(
                            qfield[nd][i], numDerivFwd[nd][i], 
                            numDerivBwd[nd][i]);
                        TraceMap->UniversalTraceAssemble(numDerivBwd[nd][i]);
                        TraceMap->UniversalTraceAssemble(numDerivFwd[nd][i]);
                    }
                }

                NekVector<DataType> qFwd(nvariables*nDim),qBwd(nvariables*nDim);

                for(int i = 0; i < nvariables; ++i)
                {
                    Vmath::Zero(nTracePts,FwdFluxDeriv[i],1);
                    Vmath::Zero(nTracePts,BwdFluxDeriv[i],1);
                }  

                NekVector<DataType> VFlux(nvariables);
                for(int n = 0; n < nTracePts; ++n)
                {
                    
                    for (int nd = 0; nd < nDim; ++nd)
                    {
                        for(int j = 0; j< nvariables;j++)
                        {
                            qFwd[j*nDim+nd] = DataType(numDerivFwd[nd][j][n]);
                            qBwd[j*nDim+nd] = DataType(numDerivBwd[nd][j][n]);
                        }
                    }
                
                    VFlux   =   (*TraceJacDeriv[0]->GetBlock(n,n))*qFwd;
                    for(int i = 0; i < nvariables; ++i)
                    {
                        FwdFluxDeriv[i][n] +=  NekDouble( 
                            TraceJacDerivSign[0][n]*VFlux[i]);
                    }  

                    VFlux   =   (*TraceJacDeriv[1]->GetBlock(n,n))*qBwd;
                    for(int i = 0; i < nvariables; ++i)
                    {
                        BwdFluxDeriv[i][n] +=  NekDouble( 
                            TraceJacDerivSign[1][n]*VFlux[i]);
                    }  
                }
            }

            for(int i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(nTracePts,FwdFluxDeriv[i],1,
                    FwdFlux[i],1,FwdFlux[i],1);
                Vmath::Vadd(nTracePts,BwdFluxDeriv[i],1,
                    BwdFlux[i],1,BwdFlux[i],1);
            }  
        }
#ifdef CFS_DEBUGMODE
        }
#endif
        // Evaulate <\phi, \hat{F}\cdot n> - OutField[i]
        for(int i = 0; i < nvariables; ++i)
        {
            Vmath::Fill(nCoeffs,0.0,outarray[i],1);
            // Vmath::Neg                                  (nCoeffs, outarray[i], 1);
            m_fields[i]->AddTraceIntegralToOffDiag(FwdFlux[i],BwdFlux[i], 
                outarray[i]);
            // m_fields[i]->MultiplyByElmtInvMass          (outarray[i], outarray[i]);
        }
        for(int i = 0; i < nvariables; ++i)
        {
            Vmath::Svtvp(nCoeffs,-m_TimeIntegLambda,outarray[i],1,inarray[i],1,
                outarray[i],1);
        }
    }

    template<typename DataType>
    void CompressibleFlowSystem::AddOrMinusPrecTraceFlux(
        const int                            &nswitchDiagOffdiag,
        const NekDouble                      AddOrMinusSign,
        const TensorOfArray2D<NekDouble>     &inarray,
        TensorOfArray2D<NekDouble>           &outarray,
        bool                                 flagUpdateDervFlux,
        TensorOfArray2D<NekDouble>           &FwdFluxDeriv,
        TensorOfArray2D<NekDouble>           &BwdFluxDeriv,
        TensorOfArray3D<NekDouble>           &qfield,
        TensorOfArray3D<NekDouble>           &wspTrace,
        TensorOfArray2D<DataType>            &wspTraceDataType,
        const TensorOfArray4D<DataType>      &TraceJacArray,
        const TensorOfArray4D<DataType>      &TraceJacDerivArray,
        const TensorOfArray2D<DataType>      &TraceJacDerivSign,
        const TensorOfArray5D<DataType>      &TraceIPSymJacArray)
    {
#ifdef CFS_DEBUGMODE
        if(1!=m_DebugVolTraceSwitch)
        {
#endif
        size_t nvariables = m_fields.num_elements();
        size_t nCoeffs    = GetNcoeffs();
        size_t npoints    = GetNpoints();
        size_t nDim       = m_spacedim;
        size_t nTracePts  = GetTraceTotPoints();

        TensorOfArray2D<NekDouble> FwdFlux;
        TensorOfArray2D<NekDouble> BwdFlux;
        TensorOfArray2D<NekDouble> Fwd;
        TensorOfArray2D<NekDouble> Bwd;

        int indexwspTrace = 0;
        FwdFlux = wspTrace[indexwspTrace];
        indexwspTrace++;
        BwdFlux = wspTrace[indexwspTrace];
        indexwspTrace++;
        Fwd = wspTrace[indexwspTrace];
        indexwspTrace++;
        Bwd = wspTrace[indexwspTrace];
        indexwspTrace++;

        TensorOfArray3D<NekDouble> numDerivBwd(nDim);
        TensorOfArray3D<NekDouble> numDerivFwd(nDim);
       
        TensorOfArray2D<NekDouble> outpnts(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            outpnts[i] = TensorOfArray1D<NekDouble> (npoints,0.0);
            m_fields[i]->BwdTrans(outarray[i],outpnts[i]);
            m_fields[i]->GetFwdBwdTracePhys(outpnts[i], Fwd[i], Bwd[i]);
        }

        if(m_DEBUG_VISCOUS_JAC_MAT&&m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT)
        {
            if(flagUpdateDervFlux)
            {
                CalphysDeriv(outpnts,qfield);
            }
        }

        for(int i = 0; i < nvariables; i++)
        {
            outpnts[i] = NullNekDouble1DArray;
        }

        CalTraceFwdBwdFlux(Fwd, Bwd, FwdFlux, BwdFlux, wspTraceDataType,
            TraceJacArray);

        if (m_DEBUG_VISCOUS_JAC_MAT && 
            m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT)
        {
            if (flagUpdateDervFlux || m_DebugIPSymmFluxJacSwitch)
            {
                
                for (int nd = 0; nd < nDim; ++nd)
                {
                    numDerivBwd[nd] =  wspTrace[indexwspTrace];
                    indexwspTrace++;
                    numDerivFwd[nd] =  wspTrace[indexwspTrace];
                    indexwspTrace++;
                }
            }

            if (flagUpdateDervFlux )
            {
                AddTraceFwdBwdDerivFlux(FwdFlux, BwdFlux, numDerivFwd, 
                    numDerivBwd, qfield, wspTraceDataType, TraceJacDerivArray, 
                    TraceJacDerivSign);
            }
        }   

        switch (nswitchDiagOffdiag)
        {
        case ePrecDiagOffDiagTagAll:
            for(int i = 0; i < nvariables; ++i)
            {
                Vmath::Fill(nCoeffs,0.0,outarray[i],1);
                Vmath::Vadd(nTracePts, 
                        FwdFlux[i], 1, 
                        BwdFlux[i], 1, 
                        FwdFlux[i], 1);
                m_fields[i]->AddTraceIntegral(FwdFlux[i], outarray[i]);
            }
            break;
        case ePrecDiagOffDiagTagDiag:
            for(int i = 0; i < nvariables; ++i)
            {
                Vmath::Fill(nCoeffs,0.0,outarray[i],1);
                m_fields[i]->AddTraceIntegralToDiag(FwdFlux[i],BwdFlux[i], 
                    outarray[i]);
            }
            break;
        case ePrecDiagOffDiagTagOffdiag:
            for(int i = 0; i < nvariables; ++i)
            {
                Vmath::Fill(nCoeffs,0.0,outarray[i],1);
                m_fields[i]->AddTraceIntegralToOffDiag(FwdFlux[i],BwdFlux[i], 
                    outarray[i]);
            }
            break;
        default:
            for(int i = 0; i < nvariables; ++i)
            {
                Vmath::Fill(nCoeffs,0.0,outarray[i],1);
                Vmath::Vadd(nTracePts, 
                        FwdFlux[i], 1, 
                        BwdFlux[i], 1, 
                        FwdFlux[i], 1);
                m_fields[i]->AddTraceIntegral(FwdFlux[i], outarray[i]);
            }
            break;
        }
        
#ifdef CFS_DEBUGMODE
        if (m_DEBUG_VISCOUS_JAC_MAT &&
            m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT &&
            m_DebugIPSymmFluxJacSwitch)
        {
            CalTraceFwdBwdSymFlux(Fwd, Bwd, numDerivFwd, numDerivBwd,
                wspTraceDataType, TraceJacDerivSign, TraceIPSymJacArray);

            Array<OneD, int > nonZeroIndex;
            size_t nTracePts  = GetTraceTotPoints();
            switch (nswitchDiagOffdiag)
            {
            case ePrecDiagOffDiagTagAll:
                for (int nd = 0; nd < nDim; ++nd)
                {
                    for(int i = 0; i < nvariables; ++i)
                    {
                        Vmath::Vadd(nTracePts, 
                            numDerivFwd[nd][i], 1, 
                            numDerivBwd[nd][i], 1, 
                            numDerivFwd[nd][i], 1);
                    }
                }
                m_diffusion->v_AddSymmFluxIntegralToCoeff(nvariables,
                    nDim,npoints,nTracePts,m_fields,nonZeroIndex,numDerivFwd,
                    outarray);
                break;
            case ePrecDiagOffDiagTagDiag:
                // for(int i = 0; i < nvariables; ++i)
                // {
                //     m_diffusion->AddSymmFluxIntegralToDiag(nvariables,
                //         nDim,npoints,nTracePts,m_fields,nonZeroIndex,numDerivFwd,
                //         numDerivBwd,outarray);
                // }
                ASSERTL0(false, "AddSymmFluxIntegralToDiag not coded");
                break;
            case ePrecDiagOffDiagTagOffdiag:
                m_diffusion->AddSymmFluxIntegralToOffDiag(nvariables,
                    nDim,npoints,nTracePts,m_fields,nonZeroIndex,numDerivFwd,
                    numDerivBwd,outarray);
                break;
            default:
                for (int nd = 0; nd < nDim; ++nd)
                {
                    for(int i = 0; i < nvariables; ++i)
                    {
                        Vmath::Vadd(nTracePts, 
                            numDerivFwd[nd][i], 1, 
                            numDerivBwd[nd][i], 1, 
                            numDerivFwd[nd][i], 1);
                    }
                }
                m_diffusion->v_AddSymmFluxIntegralToCoeff(nvariables,
                    nDim,npoints,nTracePts,m_fields,nonZeroIndex,numDerivFwd,
                    outarray);
                break;
            }
        }
#endif
                for(int i = 0; i < nvariables; ++i)
        {
            Vmath::Svtvp(nCoeffs,AddOrMinusSign*m_TimeIntegLambda,
                        outarray[i],1,
                        inarray[i],1,
                        outarray[i],1);
        }
#ifdef CFS_DEBUGMODE
        }
#endif
    }

    template<typename DataType>
    void CompressibleFlowSystem::CalTraceFwdBwdFlux(
        const TensorOfArray2D<NekDouble>         &Fwd,
        const TensorOfArray2D<NekDouble>         &Bwd,
        TensorOfArray2D<NekDouble>               &FwdFlux,
        TensorOfArray2D<NekDouble>               &BwdFlux,
        TensorOfArray2D<DataType>                &wspTraceDataType,
        const TensorOfArray4D<DataType>          &TraceJacArray)
    {
        size_t nvariables = m_fields.num_elements();
        size_t nTracePts  = GetTraceTotPoints();

        int indexwspTraceDataType = 0;
        TensorOfArray2D<DataType> Fwdarray (nvariables);
        for(int m = 0; m < nvariables; ++m)
        {
            Fwdarray[m] = wspTraceDataType[indexwspTraceDataType];
            indexwspTraceDataType++;
        }
        TensorOfArray1D<DataType> Fwdreslt;
        Fwdreslt = wspTraceDataType[indexwspTraceDataType];
        indexwspTraceDataType++;

        for(int m = 0; m < nvariables; ++m)
        {
            for(int i = 0; i < nTracePts; ++i)
            {
                Fwdarray[m][i] = DataType(Fwd[m][i]);
            }
        }
        for(int m = 0; m < nvariables; ++m)
        {
            Vmath::Zero(nTracePts, &Fwdreslt[0],1);
            for(int n = 0; n < nvariables; ++n)
            {
                Vmath::Vvtvp(nTracePts,&TraceJacArray[0][m][n][0],1,
                    &Fwdarray[n][0],1,&Fwdreslt[0],1,&Fwdreslt[0],1);
            }

            for(int i = 0; i < nTracePts; ++i)
            {
                FwdFlux[m][i] =  NekDouble( Fwdreslt[i] );
            }
        }

        for(int m = 0; m < nvariables; ++m)
        {
            for(int i = 0; i < nTracePts; ++i)
            {
                Fwdarray[m][i] =  DataType( Bwd[m][i] );
            }
        }
        for(int m = 0; m < nvariables; ++m)
        {
            Vmath::Zero(nTracePts, &Fwdreslt[0],1);
            for(int n = 0; n < nvariables; ++n)
            {
                Vmath::Vvtvp(nTracePts,&TraceJacArray[1][m][n][0],1,
                    &Fwdarray[n][0],1,&Fwdreslt[0],1,&Fwdreslt[0],1);
            }
            for(int i = 0; i < nTracePts; ++i)
            {
                BwdFlux[m][i] =  NekDouble( Fwdreslt[i] );
            }
        }
    }

    template<typename DataType>
    void CompressibleFlowSystem::AddTraceFwdBwdDerivFlux(
        TensorOfArray2D<NekDouble>               &FwdFlux,
        TensorOfArray2D<NekDouble>               &BwdFlux,
        TensorOfArray3D<NekDouble>               &numDerivFwd,
        TensorOfArray3D<NekDouble>               &numDerivBwd,
        TensorOfArray3D<NekDouble>               &qfield,
        TensorOfArray2D<DataType>                &wspTraceDataType,
        const TensorOfArray4D<DataType>          &TraceJacDerivArray,
        const TensorOfArray2D<DataType>          &TraceJacDerivSign)
    {
        size_t nTracePts  = GetTraceTotPoints();
        size_t nDim       = m_spacedim;
        size_t nvariables = m_fields.num_elements();

        int indexwspTraceDataType = 0;
        TensorOfArray2D<DataType> Fwdarray (nvariables);
        for(int m = 0; m < nvariables; ++m)
        {
            Fwdarray[m] = wspTraceDataType[indexwspTraceDataType];
            indexwspTraceDataType++;
        }
        TensorOfArray1D<DataType> Fwdreslt;
        Fwdreslt = wspTraceDataType[indexwspTraceDataType];
        indexwspTraceDataType++;

        const MultiRegions::AssemblyMapDGSharedPtr TraceMap = 
            m_fields[0]->GetTraceMap();
        for (int nd = 0; nd < nDim; ++nd)
        {
            for (int i = 0; i < nvariables; ++i)
            {
                m_fields[i]->GetFwdBwdTracePhysNoBndFill(qfield[nd][i], 
                    numDerivFwd[nd][i], numDerivBwd[nd][i]);
                TraceMap->UniversalTraceAssemble(numDerivBwd[nd][i]);
                TraceMap->UniversalTraceAssemble(numDerivFwd[nd][i]);
            }
        }

        for(int m = 0; m < nvariables; ++m)
        {
            Vmath::Zero(nTracePts, &Fwdreslt[0],1);
            for (int nd = 0; nd < nDim; ++nd)
            {
                for(int n = 0; n < nvariables; ++n)
                {
                    for(int i = 0; i < nTracePts; ++i)
                    {
                        Fwdarray[n][i] =  
                            DataType( numDerivFwd[nd][n][i] );
                    }
                    Vmath::Vvtvp(nTracePts,
                        &TraceJacDerivArray[0][m][n*nDim+nd][0],1,
                        &Fwdarray[n][0],1,
                        &Fwdreslt[0],1,
                        &Fwdreslt[0],1);
                }
            }
            for(int i = 0; i < nTracePts; ++i)
            {
                FwdFlux[m][i] += 
                    NekDouble(TraceJacDerivSign[0][i]*Fwdreslt[i] );
            }
        }

        for(int m = 0; m < nvariables; ++m)
        {
            Vmath::Zero(nTracePts, &Fwdreslt[0],1);
            for (int nd = 0; nd < nDim; ++nd)
            {
                for(int n = 0; n < nvariables; ++n)
                {
                    for(int i = 0; i < nTracePts; ++i)
                    {
                        Fwdarray[n][i] = 
                            DataType(numDerivBwd[nd][n][i]);
                    }
                    Vmath::Vvtvp(nTracePts,
                        &TraceJacDerivArray[1][m][n*nDim+nd][0],1,
                        &Fwdarray[n][0],1,
                        &Fwdreslt[0],1,
                        &Fwdreslt[0],1);
                }
            }
            for(int i = 0; i < nTracePts; ++i)
            {
                BwdFlux[m][i] +=  
                    NekDouble( TraceJacDerivSign[1][i]*Fwdreslt[i] );
            }
        }
    }

    template<typename DataType>
    void CompressibleFlowSystem::CalTraceFwdBwdSymFlux(
        const TensorOfArray2D<NekDouble>         &Fwd,
        const TensorOfArray2D<NekDouble>         &Bwd,
        TensorOfArray3D<NekDouble>               &numDerivFwd,
        TensorOfArray3D<NekDouble>               &numDerivBwd,
        TensorOfArray2D<DataType>                &wspTraceDataType,
        const TensorOfArray2D<DataType>          &TraceJacDerivSign,
        const TensorOfArray5D<DataType>          &TraceIPSymJacArray)
    {
        size_t nTracePts  = GetTraceTotPoints();
        size_t nvariables = m_fields.num_elements();

        int indexwspTraceDataType = 0;
        TensorOfArray2D<DataType> Fwdarray (nvariables);
        for(int m = 0; m < nvariables; ++m)
        {
            Fwdarray[m] = wspTraceDataType[indexwspTraceDataType];
            indexwspTraceDataType++;
        }
        TensorOfArray1D<DataType> Fwdreslt;
        Fwdreslt = wspTraceDataType[indexwspTraceDataType];
        indexwspTraceDataType++;

        for(int m = 0; m < nvariables; ++m)
        {
            for(int i = 0; i < nTracePts; ++i)
            {
                Fwdarray[m][i] =  DataType( Fwd[m][i] );
            }
        }
        for(int nd = 0; nd < m_spacedim; ++nd)
        {
            for(int m = 0; m < nvariables; ++m)
            {
                Vmath::Zero(nTracePts, &Fwdreslt[0],1);
                for(int n = 0; n < nvariables; ++n)
                {
                    Vmath::Vvtvp(nTracePts,
                        &TraceIPSymJacArray[0][nd][m][n][0],1,
                        &Fwdarray[n][0],1,
                        &Fwdreslt[0],1,
                        &Fwdreslt[0],1);
                }

                for(int i = 0; i < nTracePts; ++i)
                {
                    numDerivFwd[nd][m][i] =  NekDouble( 
                        TraceJacDerivSign[0][i]*Fwdreslt[i] );
                }
            }
        }

        for(int m = 0; m < nvariables; ++m)
        {
            for(int i = 0; i < nTracePts; ++i)
            {
                Fwdarray[m][i] =  DataType( Bwd[m][i] );
            }
        }
        for(int nd = 0; nd < m_spacedim; ++nd)
        {
            for(int m = 0; m < nvariables; ++m)
            {
                Vmath::Zero(nTracePts, &Fwdreslt[0],1);
                for(int n = 0; n < nvariables; ++n)
                {
                    Vmath::Vvtvp(nTracePts,
                        &TraceIPSymJacArray[1][nd][m][n][0],1,
                        &Fwdarray[n][0],1,
                        &Fwdreslt[0],1,
                        &Fwdreslt[0],1);
                }

                for(int i = 0; i < nTracePts; ++i)
                {
                    numDerivBwd[nd][m][i] = NekDouble(
                        TraceJacDerivSign[1][i]*Fwdreslt[i] );
                }
            }
        }
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::ElmtVarInvMtrx(
        TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
        TypeNekBlkMatSharedPtr                  &gmtVar,
        const DataType                          &tmpDatatype)
    {
        int n1d = gmtxarray.num_elements();
        int n2d = gmtxarray[0].num_elements();
        int nConvectiveFields = n1d;

        ASSERTL0(n1d==n2d,"ElmtVarInvMtrx requires n1d==n2d");

        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;

        gmtxarray[0][0]->GetBlockSizes(rowSizes,colSizes);
        int ntotElmt  = rowSizes.num_elements();
        int nElmtCoef   =    rowSizes[0]-1;
        int nElmtCoef0  =    -1;
        int blocksize = -1;

        Array<OneD, unsigned int> tmprow(1);
        TypeNekBlkMatSharedPtr tmpGmtx;

        Array<OneD, DataType>    GMatData,ElmtMatData;

        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            int nrows = gmtxarray[0][0]->GetBlock(nelmt,nelmt)->GetRows();
            int ncols = gmtxarray[0][0]->GetBlock(nelmt,nelmt)->GetColumns();
            ASSERTL0(nrows==ncols,"ElmtVarInvMtrx requires nrows==ncols");

            nElmtCoef            = nrows;

            if (nElmtCoef0!=nElmtCoef)
            {
                nElmtCoef0 = nElmtCoef;
                int nElmtCoefVAr = nElmtCoef0*nConvectiveFields;
                blocksize = nElmtCoefVAr*nElmtCoefVAr;
                tmprow[0] = nElmtCoefVAr;
                AllocateNekBlkMatDig(tmpGmtx,tmprow,tmprow);
                // tmpGmtx = MemoryManager<SNekMat>
                //     ::AllocateSharedPtr(nElmtCoefVAr, nElmtCoefVAr,0.0);
                GMatData = tmpGmtx->GetBlock(0,0)->GetPtr();

            }

            for(int n = 0; n < nConvectiveFields; n++)
            {
                for(int m = 0; m < nConvectiveFields; m++)
                {
                    ElmtMatData = 
                        gmtxarray[m][n]->GetBlock(nelmt,nelmt)->GetPtr();

                    for(int ncl = 0; ncl < nElmtCoef; ncl++)
                    {
                        int Goffset = (n*nElmtCoef+ncl)*nConvectiveFields*
                            nElmtCoef+m*nElmtCoef;
                        int Eoffset = ncl*nElmtCoef;

                        Vmath::Vcopy(nElmtCoef,&ElmtMatData[0]+Eoffset,1, 
                            &GMatData[0]+Goffset,1);
                    }
                }
            }

            tmpGmtx->GetBlock(0,0)->Invert();

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    ElmtMatData = 
                        gmtxarray[m][n]->GetBlock(nelmt,nelmt)->GetPtr();

                    for(int ncl = 0; ncl < nElmtCoef; ncl++)
                    {
                        int Goffset = (n*nElmtCoef+ncl)*nConvectiveFields*
                            nElmtCoef+m*nElmtCoef;
                        int Eoffset = ncl*nElmtCoef;

                        Vmath::Vcopy(nElmtCoef, &GMatData[0]+Goffset,1,
                            &ElmtMatData[0]+Eoffset,1);
                    }
                }
            }
            ElmtMatData = gmtVar->GetBlock(nelmt,nelmt)->GetPtr();
            Vmath::Vcopy(blocksize, &GMatData[0],1,&ElmtMatData[0],1);
        }
        return;
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::ElmtVarInvMtrx(
        TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
        const DataType                                    &tmpDatatype)
    {
        int n1d = gmtxarray.num_elements();
        int n2d = gmtxarray[0].num_elements();
        int nConvectiveFields = n1d;

        ASSERTL0(n1d==n2d,"ElmtVarInvMtrx requires n1d==n2d");

        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;

        gmtxarray[0][0]->GetBlockSizes(rowSizes,colSizes);
        int ntotElmt  = rowSizes.num_elements();
        int nElmtCoef   =    rowSizes[0]-1;
        int nElmtCoef0  =    -1;
        int blocksize = -1;

        Array<OneD, unsigned int> tmprow(1);
        TypeNekBlkMatSharedPtr tmpGmtx;

        Array<OneD, DataType>    GMatData,ElmtMatData;

        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            int nrows = gmtxarray[0][0]->GetBlock(nelmt,nelmt)->GetRows();
            int ncols = gmtxarray[0][0]->GetBlock(nelmt,nelmt)->GetColumns();
            ASSERTL0(nrows==ncols,"ElmtVarInvMtrx requires nrows==ncols");

            nElmtCoef            = nrows;
            
            if (nElmtCoef0!=nElmtCoef) 
            {
                nElmtCoef0 = nElmtCoef;
                int nElmtCoefVAr = nElmtCoef0*nConvectiveFields;
                blocksize = nElmtCoefVAr*nElmtCoefVAr;
                tmprow[0] = nElmtCoefVAr;
                AllocateNekBlkMatDig(tmpGmtx,tmprow,tmprow);
                // tmpGmtx = MemoryManager<SNekMat>
                //     ::AllocateSharedPtr(nElmtCoefVAr, nElmtCoefVAr,0.0);
                GMatData = tmpGmtx->GetBlock(0,0)->GetPtr();

            }

            for(int n = 0; n < nConvectiveFields; n++)
            {
                for(int m = 0; m < nConvectiveFields; m++)
                {
                    ElmtMatData = 
                        gmtxarray[m][n]->GetBlock(nelmt,nelmt)->GetPtr();

                    for(int ncl = 0; ncl < nElmtCoef; ncl++)
                    {
                        int Goffset = (n*nElmtCoef+ncl)*nConvectiveFields*
                            nElmtCoef+m*nElmtCoef;
                        int Eoffset = ncl*nElmtCoef;

                        Vmath::Vcopy(nElmtCoef,&ElmtMatData[0]+Eoffset,1, 
                            &GMatData[0]+Goffset,1);
                    }
                }
            }

            tmpGmtx->GetBlock(0,0)->Invert();

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    ElmtMatData = 
                        gmtxarray[m][n]->GetBlock(nelmt,nelmt)->GetPtr();

                    for(int ncl = 0; ncl < nElmtCoef; ncl++)
                    {
                        int Goffset = (n*nElmtCoef+ncl)*nConvectiveFields*
                            nElmtCoef+m*nElmtCoef;
                        int Eoffset = ncl*nElmtCoef;

                        Vmath::Vcopy(nElmtCoef, &GMatData[0]+Goffset,1,
                            &ElmtMatData[0]+Eoffset,1);
                    }
                }
            }
        }
        return;
    }

    void CompressibleFlowSystem::CalcFluxJacVolBnd(
        const TensorOfArray2D<NekDouble>  &inarray,
        const NekDouble                   time)
    {
        size_t nSpaceDim = m_graph->GetSpaceDimension();
        size_t nvariables = m_fields.num_elements();
        size_t npoints = GetTotPoints();

        Array<OneD, Array<OneD, NekDouble> > pnts(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            pnts[i]   =   Array<OneD, NekDouble>(npoints,0.0);
            m_fields[i]->BwdTrans(inarray[i], pnts[i]);
        }

        for(int j = 0; j< m_fields.num_elements(); j++)
        {
            Vmath::Vcopy(npoints, pnts[j],1,m_MatrixFreeRefFields[j],1);
            m_fields[j]->GetFwdBwdTracePhys(m_MatrixFreeRefFields[j], 
                                            m_MatrixFreeRefFwd[j], 
                                            m_MatrixFreeRefBwd[j]);
        }

        if(m_DEBUG_VISCOUS_JAC_MAT)
        {
            m_diffusion->SetRefFields(m_MatrixFreeRefFields);
        }

        DoOdeProjection(m_MatrixFreeRefFields,m_MatrixFreeRefFields,
            time);

        CalphysDeriv(pnts,m_qfield);

        for(int nfluxDir = 0; nfluxDir < nSpaceDim; nfluxDir++)
        {
            if(m_DEBUG_ADVECTION_JAC_MAT)
            {
                GetFluxVectorJacDirctnMat(nfluxDir,pnts, m_ElmtFluxJacArray);
            }

            if(m_DEBUG_VISCOUS_JAC_MAT)
            {
                MinusDiffusionFluxJacDirctnMat(nfluxDir,pnts,m_qfield, 
                    m_ElmtFluxJacArray);
            }
        }

        GetTraceJac(pnts,m_qfield,m_TraceJac,m_TraceJacDeriv,
                    m_TraceJacDerivSign, m_TraceIPSymJacArray);
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::AddMatNSBlkDiag_volume(
        const TensorOfArray2D<NekDouble>         &inarray,
        const TensorOfArray3D<NekDouble>         &qfield,
        TensorOfArray2D<TypeNekBlkMatSharedPtr>  &gmtxarray,
        TensorOfArray4D<DataType>                &StdMatDataDBB,
        TensorOfArray5D<DataType>                &StdMatDataDBDB)
    {
        if(StdMatDataDBB.num_elements()==0)
        {
            CalcVolJacStdMat(StdMatDataDBB,StdMatDataDBDB);
        }
        
        int nSpaceDim = m_graph->GetSpaceDimension();
        int nvariable = inarray.num_elements();
        int npoints   = m_fields[0]->GetTotPoints();
        int nVar2     = nvariable*nvariable;
        std::shared_ptr<LocalRegions::ExpansionVector> expvect =    
            m_fields[0]->GetExp();
        int ntotElmt            = (*expvect).size();

        Array<OneD, NekDouble > mu                 (npoints, 0.0);
        Array<OneD, NekDouble > DmuDT              (npoints, 0.0);
        if(m_DEBUG_VISCOUS_JAC_MAT)
        {
            CalcMuDmuDT(inarray,mu,DmuDT);
        }

        Array<OneD, NekDouble> normals;
        Array<OneD, Array<OneD, NekDouble> > normal3D(3);
        for(int i = 0; i < 3; i++)
        {
            normal3D[i] = Array<OneD, NekDouble>(3,0.0);
        }
        normal3D[0][0] = 1.0;
        normal3D[1][1] = 1.0;
        normal3D[2][2] = 1.0;
        Array<OneD, Array<OneD, NekDouble> > normalPnt(3);
        
        DNekMatSharedPtr wspMat     = MemoryManager<DNekMat>::
            AllocateSharedPtr(nvariable,nvariable,0.0);
        DNekMatSharedPtr wspMatDrv  = MemoryManager<DNekMat>::
            AllocateSharedPtr(nvariable-1,nvariable,0.0);

        Array<OneD, DataType> GmatxData;
        Array<OneD, DataType> MatData;

        Array<OneD, NekDouble> tmppnts;
        TensorOfArray3D<NekDouble> PntJacCons(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        TensorOfArray3D<DataType>  PntJacConsStd(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        Array<OneD, Array<OneD, NekDouble> >    ConsStdd(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> >    ConsCurv(m_spacedim);
        TensorOfArray4D<NekDouble> PntJacDerv(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        TensorOfArray4D<DataType>  PntJacDervStd(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        TensorOfArray3D<NekDouble> DervStdd(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        TensorOfArray3D<NekDouble> DervCurv(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        for(int ndir=0; ndir<m_spacedim;ndir++)
        {
            PntJacDerv[ndir] = TensorOfArray3D<NekDouble>(m_spacedim);
            PntJacDervStd[ndir] = TensorOfArray3D<DataType>(m_spacedim);
            DervStdd[ndir] = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
            DervCurv[ndir] = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        }

        Array<OneD, NekDouble> locmu;
        Array<OneD, NekDouble> locDmuDT;
        Array<OneD, Array<OneD, NekDouble> > locVars(nvariable);
        TensorOfArray3D<NekDouble> locDerv(m_spacedim);
        for(int ndir=0; ndir<m_spacedim;ndir++)
        {
            locDerv[ndir]   =   Array<OneD, Array<OneD, NekDouble> >(nvariable);
        }

        int nElmtCoefOld = -1;
        for(int ne=0; ne<ntotElmt;ne++)
        {
            int nElmtCoef           = (*expvect)[ne]->GetNcoeffs();
            int nElmtCoef2          = nElmtCoef*nElmtCoef;
            int nElmtPnt            = (*expvect)[ne]->GetTotPoints();

            int nQuot = nElmtCoef2/m_nPadding;
            int nRemd = nElmtCoef2- nQuot*m_nPadding;
            int nQuotPlus=nQuot;
            if(nRemd>0)
            {
                nQuotPlus++;
            }
            int nElmtCoef2Paded = nQuotPlus*m_nPadding;

            if(nElmtPnt>PntJacCons[0].num_elements()||nElmtCoef>nElmtCoefOld)
            {
                nElmtCoefOld = nElmtCoef;
                for(int ndir=0; ndir<3;ndir++)
                {
                    normalPnt[ndir] = Array<OneD, NekDouble>(npoints,0.0);
                }
                tmppnts = Array<OneD, NekDouble>  (nElmtPnt);
                MatData = Array<OneD, DataType>  (nElmtCoef2Paded*nVar2);
                for(int ndir=0; ndir<m_spacedim;ndir++)
                {
                    ConsCurv[ndir]      = Array<OneD, NekDouble> (nElmtPnt);
                    ConsStdd[ndir]      = Array<OneD, NekDouble> (nElmtPnt);
                    PntJacCons[ndir]    = 
                        Array<OneD, Array<OneD, NekDouble> > (nElmtPnt);
                    PntJacConsStd[ndir] = 
                        Array<OneD, Array<OneD, DataType> > (nElmtPnt);
                    for(int i=0; i<nElmtPnt;i++)
                    {
                        PntJacCons[ndir][i] = Array<OneD, NekDouble>(nVar2);
                        PntJacConsStd[ndir][i] = Array<OneD, DataType>(nVar2);
                    }
                    
                    for(int ndir1=0; ndir1<m_spacedim;ndir1++)
                    {
                        PntJacDerv[ndir][ndir1] = 
                            Array<OneD, Array<OneD, NekDouble> > (nElmtPnt);
                        PntJacDervStd[ndir][ndir1] = 
                            Array<OneD, Array<OneD, DataType> > (nElmtPnt);
                        DervStdd[ndir][ndir1]       =   
                            Array<OneD, NekDouble> (nElmtPnt);
                        DervCurv[ndir][ndir1]       =   
                            Array<OneD, NekDouble> (nElmtPnt);
                        for(int i=0; i<nElmtPnt;i++)
                        {
                            PntJacDerv[ndir][ndir1][i] = 
                                Array<OneD, NekDouble>(nVar2);
                            PntJacDervStd[ndir][ndir1][i] =
                                Array<OneD, DataType>(nVar2);
                        }
                    }
                }
            }
    
            int noffset = GetPhys_Offset(ne);
            for(int j = 0; j < nvariable; j++)
            {   
                locVars[j] = inarray[j]+noffset;
            }
            
            if(m_DEBUG_ADVECTION_JAC_MAT)
            {
                for(int nfluxDir = 0; nfluxDir < nSpaceDim; nfluxDir++)
                {
                    normals =   normal3D[nfluxDir];
                    GetFluxVectorJacDirctnElmt(nvariable,nElmtPnt,locVars,
                        normals,wspMat,PntJacCons[nfluxDir]);
                }
            }

            if(m_DEBUG_VISCOUS_JAC_MAT)
            {
                for(int j = 0; j < nSpaceDim; j++)
                {   
                    for(int k = 0; k < nvariable; k++)
                    {
                        locDerv[j][k] = qfield[j][k]+noffset;
                    }
                }
                locmu       =   mu      + noffset;
                locDmuDT    =   DmuDT   + noffset;
                for(int nfluxDir = 0; nfluxDir < nSpaceDim; nfluxDir++)
                {
                    normals =   normal3D[nfluxDir];
                    MinusDiffusionFluxJacDirctnElmt(nvariable,nElmtPnt,locVars,
                        locDerv,locmu,locDmuDT,normals,wspMatDrv,
                        PntJacCons[nfluxDir]);
                }
            }

            if(m_DEBUG_VISCOUS_JAC_MAT)
            {
                locmu       =   mu      + noffset;
                for(int nfluxDir = 0; nfluxDir < nSpaceDim; nfluxDir++)
                {
                    Vmath::Fill(npoints,1.0,normalPnt[nfluxDir],1);
                    for(int nDervDir = 0; nDervDir < nSpaceDim; nDervDir++)
                    {
                        GetFluxDerivJacDirctnElmt(nvariable,nElmtPnt,nDervDir,
                            locVars,locmu,normalPnt,wspMatDrv,
                            PntJacDerv[nfluxDir][nDervDir]);
                    }
                    Vmath::Fill(npoints,0.0,normalPnt[nfluxDir],1);
                }
            }

            for(int n=0; n<nvariable;n++)
            {
                for(int m=0; m<nvariable;m++)
                {
                    int nvarOffset = m+n*nvariable;
                    GmatxData = gmtxarray[m][n]->GetBlock(ne,ne)->GetPtr();

                    for(int ndStd0 =0;ndStd0<m_spacedim;ndStd0++)
                    {
                        Vmath::Zero(nElmtPnt,ConsStdd[ndStd0],1);
                    }
                    for(int ndir =0;ndir<m_spacedim;ndir++)
                    {
                        for(int i=0; i<nElmtPnt;i++)
                        {
                            tmppnts[i] =  PntJacCons[ndir][i][nvarOffset];
                        }
                        (*expvect)[ne]->ProjectVectorintoStandardExp(ndir,
                            tmppnts,ConsCurv);
                        for(int nd =0;nd<m_spacedim;nd++)
                        {
                            Vmath::Vadd(nElmtPnt,ConsCurv[nd],1,ConsStdd[nd],1,
                                ConsStdd[nd],1);
                        }
                    }

                    for(int ndir =0;ndir<m_spacedim;ndir++)
                    {
                        (*expvect)[ne]->MultiplyByQuadratureMetric(
                                ConsStdd[ndir],ConsStdd[ndir]); // weight with metric
                        for(int i=0; i<nElmtPnt;i++)
                        {
                            PntJacConsStd[ndir][i][nvarOffset] = 
                                DataType(ConsStdd[ndir][i]);
                        }
                    }
                }
            }

            if(m_DEBUG_VISCOUS_JAC_MAT)
            {
                for(int m=0; m<nvariable;m++)
                {
                    for(int n=0; n<nvariable;n++)
                    {
                        int nvarOffset = m+n*nvariable;
                        for(int ndStd0 =0;ndStd0<m_spacedim;ndStd0++)
                        {
                            for(int ndStd1 =0;ndStd1<m_spacedim;ndStd1++)
                            {
                                Vmath::Zero(nElmtPnt,
                                    DervStdd[ndStd0][ndStd1],1);
                            }
                        }
                        for(int nd0 =0;nd0<m_spacedim;nd0++)
                        {
                            for(int nd1 =0;nd1<m_spacedim;nd1++)
                            {
                                for(int i=0; i<nElmtPnt;i++)
                                {
                                    tmppnts[i] =  
                                        PntJacDerv[nd0][nd1][i][nvarOffset];
                                }

                                (*expvect)[ne]->ProjectVectorintoStandardExp(
                                        nd0,tmppnts,ConsCurv);
                                for(int nd =0;nd<m_spacedim;nd++)
                                {
                                    (*expvect)[ne]->
                                        ProjectVectorintoStandardExp(nd1,
                                        ConsCurv[nd],DervCurv[nd]);
                                }

                                for(int ndStd0 =0;ndStd0<m_spacedim;ndStd0++)
                                {
                                    for(int ndStd1 =0;ndStd1<m_spacedim;ndStd1++)
                                    {
                                        Vmath::Vadd(nElmtPnt,
                                            DervCurv[ndStd0][ndStd1],1,
                                            DervStdd[ndStd0][ndStd1],1,
                                            DervStdd[ndStd0][ndStd1],1);
                                    }
                                }
                            }
                        }
                        for(int nd0 =0;nd0<m_spacedim;nd0++)
                        {
                            for(int nd1 =0;nd1<m_spacedim;nd1++)
                            {
                                (*expvect)[ne]->MultiplyByQuadratureMetric(
                                    DervStdd[nd0][nd1],DervStdd[nd0][nd1]); // weight with metric
                                for(int i=0; i<nElmtPnt;i++)
                                {
                                    PntJacDervStd[nd0][nd1][i][nvarOffset] =
                                        -DataType(DervStdd[nd0][nd1][i]);
                                }
                            }
                        }
                    }
                }
            }
            
            Vmath::Zero(nElmtCoef2Paded*nVar2,&MatData[0],1);
            DataType one = 1.0;
            for(int ndir =0;ndir<m_spacedim;ndir++)
            {
                for(int i=0;i<nElmtPnt;i++)
                {
                    Blas::DoSger (nElmtCoef2Paded,nVar2,one,
                                &StdMatDataDBB[ne][ndir][i][0],1,
                                &PntJacConsStd[ndir][i][0],1,
                                &MatData[0],nElmtCoef2Paded);
                }
            }

            if(m_DEBUG_VISCOUS_JAC_MAT)
            {
                for(int nd0 =0;nd0<m_spacedim;nd0++)
                {
                    for(int nd1 =0;nd1<m_spacedim;nd1++)
                    {
                        for(int i=0;i<nElmtPnt;i++)
                        {
                            Blas::DoSger (nElmtCoef2Paded,nVar2,one,
                                        &StdMatDataDBDB[ne][nd0][nd1][i][0],1,
                                        &PntJacDervStd[nd0][nd1][i][0],1,
                                        &MatData[0],nElmtCoef2Paded);
                        }
                    }
                }
            }

            for(int n=0; n<nvariable;n++)
            {
                for(int m=0; m<nvariable;m++)
                {
                    int nvarOffset = m+n*nvariable;
                    GmatxData = gmtxarray[m][n]->GetBlock(ne,ne)->GetPtr();
                    Vmath::Vcopy(nElmtCoef2,&MatData[nvarOffset*nElmtCoef2Paded]
                        ,1,&GmatxData[0],1);
                }
            }
        }
    }

    template<typename DataType>
    void CompressibleFlowSystem::CalcVolJacStdMat(
        TensorOfArray4D<DataType> &StdMatDataDBB,
        TensorOfArray5D<DataType> &StdMatDataDBDB)
    {
        int nSpaceDim = m_graph->GetSpaceDimension();
        int nvariable = m_fields.num_elements();
        int npoints   = m_fields[0]->GetTotPoints();
        int nVar2     = nvariable*nvariable;
        std::shared_ptr<LocalRegions::ExpansionVector> expvect =
             m_fields[0]->GetExp();
        int ntotElmt            = (*expvect).size();

        StdMatDataDBB   = TensorOfArray4D<DataType> (ntotElmt);
        StdMatDataDBDB  = TensorOfArray5D<DataType> (ntotElmt);

        vector<DNekMatSharedPtr> VectStdDerivBase0;
        vector< TensorOfArray3D<DataType> > VectStdDerivBase_Base;
        vector< TensorOfArray4D<DataType> > VectStdDervBase_DervBase;
        DNekMatSharedPtr MatStdDerivBase0;
        Array<OneD, DNekMatSharedPtr>           ArrayStdMat(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> >    ArrayStdMatData(m_spacedim);
        for(int ne=0; ne<ntotElmt;ne++)
        {
            StdRegions::StdExpansionSharedPtr stdExp;
            stdExp = (*expvect)[ne]->GetStdExp();
            StdRegions::StdMatrixKey  matkey(StdRegions::eDerivBase0,
                                stdExp->DetShapeType(), *stdExp);
            MatStdDerivBase0      =   stdExp->GetStdMatrix(matkey);

            int ntotStdExp = VectStdDerivBase0.size();
            int nfoundStdExp = -1;
            for(int i=0;i<ntotStdExp;i++) 
            {
                if((*VectStdDerivBase0[i])==(*MatStdDerivBase0))
                {
                    nfoundStdExp = i;
                }
            }
            if(nfoundStdExp>=0)
            {
                StdMatDataDBB[ne] = VectStdDerivBase_Base[nfoundStdExp];
                StdMatDataDBDB[ne] = VectStdDervBase_DervBase[nfoundStdExp];
            }
            else
            {
                int nElmtCoef           = (*expvect)[ne]->GetNcoeffs();
                int nElmtCoef2          = nElmtCoef*nElmtCoef;
                int nElmtPnt            = (*expvect)[ne]->GetTotPoints();

                int nQuot = nElmtCoef2/m_nPadding;
                int nRemd = nElmtCoef2- nQuot*m_nPadding;
                int nQuotPlus=nQuot;
                if(nRemd>0)
                {
                    nQuotPlus++;
                }
                int nPaded = nQuotPlus*m_nPadding;

                ArrayStdMat[0] = MatStdDerivBase0;
                if(m_spacedim>1)
                {
                    StdRegions::StdMatrixKey  matkey(StdRegions::eDerivBase1,
                                    stdExp->DetShapeType(), *stdExp);
                    ArrayStdMat[1]  =   stdExp->GetStdMatrix(matkey);
                    
                    if(m_spacedim>2)
                    {
                        StdRegions::StdMatrixKey matkey(StdRegions::eDerivBase2,
                                            stdExp->DetShapeType(), *stdExp);
                        ArrayStdMat[2]  =   stdExp->GetStdMatrix(matkey);
                    }
                }
                for(int nd0=0;nd0<m_spacedim;nd0++)
                {
                    ArrayStdMatData[nd0] =  ArrayStdMat[nd0]->GetPtr();
                }

                StdRegions::StdMatrixKey  matkey(StdRegions::eBwdMat,
                                        stdExp->DetShapeType(), *stdExp);
                DNekMatSharedPtr BwdMat =  stdExp->GetStdMatrix(matkey);
                Array<OneD, NekDouble> BwdMatData = BwdMat->GetPtr();

                TensorOfArray3D<DataType> tmpStdDBB (m_spacedim);
                TensorOfArray4D<DataType> tmpStdDBDB(m_spacedim);

                for(int nd0=0;nd0<m_spacedim;nd0++)
                {
                    tmpStdDBB[nd0] = 
                        Array<OneD, Array<OneD, DataType> > (nElmtPnt);
                    for(int i=0;i<nElmtPnt;i++)
                    {
                        tmpStdDBB[nd0][i] = Array<OneD, DataType> (nPaded,0.0);
                        for(int nc1=0;nc1<nElmtCoef;nc1++)
                        {
                            int noffset = nc1*nElmtCoef;
                            for(int nc0=0;nc0<nElmtCoef;nc0++)
                            {
                                tmpStdDBB[nd0][i][nc0+noffset] = DataType(
                                    ArrayStdMatData[nd0][i*nElmtCoef+nc0]*
                                    BwdMatData[i*nElmtCoef+nc1]);
                            }
                        }
                    }

                    tmpStdDBDB[nd0] = TensorOfArray3D<DataType> (m_spacedim);
                    for(int nd1=0;nd1<m_spacedim;nd1++)
                    {
                        tmpStdDBDB[nd0][nd1] = 
                            Array<OneD, Array<OneD, DataType> > (nElmtPnt);
                        for(int i=0;i<nElmtPnt;i++)
                        {
                            tmpStdDBDB[nd0][nd1][i] = 
                                Array<OneD, DataType> (nPaded,0.0);
                            for(int nc1=0;nc1<nElmtCoef;nc1++)
                            {
                                int noffset = nc1*nElmtCoef;
                                for(int nc0=0;nc0<nElmtCoef;nc0++)
                                {
                                    tmpStdDBDB[nd0][nd1][i][nc0+noffset] = 
                                        DataType(
                                        ArrayStdMatData[nd0][i*nElmtCoef+nc0]*
                                        ArrayStdMatData[nd1][i*nElmtCoef+nc1]);
                                }
                            }
                        }
                    }
                }
                VectStdDerivBase0.push_back(MatStdDerivBase0);
                VectStdDerivBase_Base.push_back(tmpStdDBB);
                VectStdDervBase_DervBase.push_back(tmpStdDBDB);

                StdMatDataDBB[ne]  = tmpStdDBB;
                StdMatDataDBDB[ne] = tmpStdDBDB;
            }
        }
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::AddMatNSBlkDiag_boundary(
        const TensorOfArray2D<NekDouble>        &inarray,
        TensorOfArray3D<NekDouble>              &qfield,
        TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
        Array<OneD, TypeNekBlkMatSharedPtr >    &TraceJac,
        Array<OneD, TypeNekBlkMatSharedPtr >    &TraceJacDeriv,
        Array<OneD, Array<OneD, DataType> >     &TraceJacDerivSign,
        TensorOfArray5D<DataType>               &TraceIPSymJacArray)
    {
        int nvariables = inarray.num_elements();
        GetTraceJac(inarray,qfield,TraceJac,TraceJacDeriv,TraceJacDerivSign,
            TraceIPSymJacArray);

#ifdef CFS_DEBUGMODE
        if(2==m_DebugAdvDiffSwitch&&2==m_DebugConsDerivSwitch)
        {
            DataType tmp=0.0;
            Fill1DArrayOfBlkDiagonalMat(TraceJac,tmp);
        }
#endif        
    
        if(m_DEBUG_VISCOUS_JAC_MAT&&m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT)
        // if(m_DEBUG_VISCOUS_JAC_MAT)
        {
                for(int i = 0; i< m_spacedim; i++)
                {
                    qfield[i]   =   NullNekDoubleArrayofArray;
                }
                m_advObject->AddTraceJacToMat(nvariables,m_spacedim,m_fields, 
                    TraceJac,gmtxarray,TraceJacDeriv,TraceJacDerivSign);
        }
        else
        {
            Array<OneD, TypeNekBlkMatSharedPtr > tmpJac;
            Array<OneD, Array<OneD, DataType> >  tmpSign;

            m_advObject->AddTraceJacToMat(nvariables,m_spacedim,m_fields, 
                TraceJac,gmtxarray,tmpJac,tmpSign);
        }
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::GetTraceJac(
        const TensorOfArray2D<NekDouble>     &inarray,
        const TensorOfArray3D<NekDouble>     &qfield,
        Array<OneD, TypeNekBlkMatSharedPtr > &TraceJac,
        Array<OneD, TypeNekBlkMatSharedPtr > &TraceJacDeriv,
        Array<OneD, Array<OneD, DataType> >  &TraceJacDerivSign,
        TensorOfArray5D<DataType>            &TraceIPSymJacArray)
    {
        int nvariables = inarray.num_elements();
        // int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

        // Store forwards/backwards space along trace space
        Array<OneD, Array<OneD, NekDouble> > Fwd    (nvariables);
        Array<OneD, Array<OneD, NekDouble> > Bwd    (nvariables);

        TypeNekBlkMatSharedPtr FJac,BJac;
        Array<OneD, unsigned int> n_blks1(nTracePts, nvariables);

        if(TraceJac.num_elements()>0)
        {
            FJac = TraceJac[0];
            BJac = TraceJac[1];
        }
        else
        {
            TraceJac    =   Array<OneD, TypeNekBlkMatSharedPtr>(2);

            AllocateNekBlkMatDig(FJac,n_blks1, n_blks1);
            AllocateNekBlkMatDig(BJac,n_blks1, n_blks1);
        }
        

                if (m_HomogeneousType == eHomogeneous1D)
        {
            Fwd = NullNekDoubleArrayofArray;
            Bwd = NullNekDoubleArrayofArray;
        }
        else
        {
            for(int i = 0; i < nvariables; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
        }

        Array<OneD, Array<OneD, NekDouble> > AdvVel(m_spacedim);
        DataType tmp;
#ifdef CFS_DEBUGMODE
        if(m_DebugIPSymmFluxJacSwitch&&(0==TraceIPSymJacArray.num_elements()))
        {
            int nTracePts   = GetTraceTotPoints();
            int nlrTot = 2;
            TraceIPSymJacArray = TensorOfArray5D<DataType> (nlrTot);
            for(int nlr = 0; nlr < nlrTot; ++nlr)
            {
                TraceIPSymJacArray[nlr] = TensorOfArray4D<DataType>(m_spacedim);
                for(int nd = 0; nd < m_spacedim; ++nd)
                {
                    TraceIPSymJacArray[nlr][nd] = 
                        TensorOfArray3D<DataType>(nvariables);
                    for(int i = 0; i < nvariables; ++i)
                    {
                        TraceIPSymJacArray[nlr][nd][i] = 
                            Array<OneD, Array<OneD, DataType> >(nvariables);
                        for(int j = 0; j < nvariables; ++j)
                        {
                            TraceIPSymJacArray[nlr][nd][i][j] = 
                                Array<OneD, DataType>(nTracePts);
                        }
                    }
                }
            }
        }
#endif

        NumCalRiemFluxJac(nvariables, m_fields, AdvVel, inarray,qfield,
                            m_BndEvaluateTime, Fwd, Bwd, FJac, BJac,
                            TraceIPSymJacArray);

        TraceJac[0] = FJac;
        TraceJac[1] = BJac;
        DataType tmpDataType;
        if(m_DEBUG_VISCOUS_JAC_MAT)
        {
            if(m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT)
            {
                int nDeriv = m_spacedim *nvariables; 
                Array<OneD, unsigned int> n_blks2(nTracePts,nDeriv);

                TypeNekBlkMatSharedPtr DerivJac;
                if(TraceJacDeriv.num_elements()>0)
                {
                    DerivJac = TraceJacDeriv[0];
                }
                else
                {
                    AllocateNekBlkMatDig(DerivJac,n_blks1, n_blks2);
                    TraceJacDeriv    = Array<OneD, TypeNekBlkMatSharedPtr>(2);
                    TraceJacDeriv[0] = DerivJac;
                    TraceJacDeriv[1] = DerivJac;
                    // usaully FJac = -BJac for viscous Jacobian of auxilary variables. 
                    // to save memory they share the array element share the same matrix 
                    // using the following sign to give their sign.
                    TraceJacDerivSign = Array<OneD, Array<OneD, DataType> >(2);
                    TraceJacDerivSign[0] = 
                        Array<OneD, DataType> (nTracePts,-1.0); 
                    TraceJacDerivSign[1] = 
                        Array<OneD, DataType> (nTracePts,-1.0); 
                }

                // Riemann flux Jacobian considering implicit boundary condition
                CalVisFluxDerivJac(nvariables, inarray,Fwd, Bwd, DerivJac,
                    tmpDataType);
            }
        }
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::CalVisFluxDerivJac(
        const int                                  nConvectiveFields,
        const Array<OneD, Array<OneD, NekDouble> > &inarray,
        const Array<OneD, Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, Array<OneD, NekDouble> > &Bwd,
        TypeNekBlkMatSharedPtr                     &BJac,
        DataType                                   &tmpDataType)
    {
        MultiRegions::ExpListSharedPtr tracelist = m_fields[0]->GetTrace();
        std::shared_ptr<LocalRegions::ExpansionVector> traceExp = 
            tracelist->GetExp();
        int ntotTrac            = (*traceExp).size();
        int nTracPnt,nTracCoef,noffset,pntoffset;

        int nTracePts  = tracelist->GetTotPoints();

        Array<OneD, NekDouble> tracBwdWeightAver(nTracePts,0.0);
        Array<OneD, NekDouble> tracBwdWeightJump(nTracePts,0.0);
        m_fields[0]->GetBwdWeight(tracBwdWeightAver,tracBwdWeightJump);
        tracBwdWeightAver = NullNekDouble1DArray;
        Vmath::Ssub(nTracePts,2.0,tracBwdWeightJump,1,tracBwdWeightJump,1);
        Vmath::Smul(nTracePts,0.5,tracBwdWeightJump,1,tracBwdWeightJump,1);
        
        Array<OneD, Array<OneD, NekDouble> > solution_jump(nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble> > solution_Aver(nConvectiveFields);
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            solution_jump[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
            solution_Aver[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
        }

        m_diffusion->ConsVarAveJump(nConvectiveFields,nTracePts,Fwd,Bwd,
            solution_Aver,solution_jump);

        const TensorOfArray2D<NekDouble> tracenormals = 
            m_diffusion->GetTraceNormal();

        TensorOfArray2D<DNekMatSharedPtr> ElmtJac;
        ElmtJac =   TensorOfArray2D<DNekMatSharedPtr> (ntotTrac);
        for(int  nelmt = 0; nelmt < ntotTrac; nelmt++)
        {
            int nTracPnt            = (*traceExp)[nelmt]->GetTotPoints();
            ElmtJac[nelmt] =   Array<OneD, DNekMatSharedPtr>(nTracPnt);
            for(int npnt = 0; npnt < nTracPnt; npnt++)
            {
                ElmtJac[nelmt][npnt] = MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nConvectiveFields, nConvectiveFields);
            }
        }
        
        DNekMatSharedPtr tmpMatinn;
        Array<OneD, NekDouble > tmpMatinnData;
        Array<OneD, DataType > tmpMatoutData;
        for(int nDervDir = 0; nDervDir < m_spacedim; nDervDir++)
        {
            GetFluxDerivJacDirctn(tracelist,tracenormals,nDervDir,solution_Aver,
                    ElmtJac);

            for(int  ntrace = 0; ntrace < ntotTrac; ntrace++)
            {
                int nTracPnt    = (*traceExp)[ntrace]->GetTotPoints();
                noffset         = tracelist->GetPhys_Offset(ntrace);
                for(int npnt = 0; npnt < nTracPnt; npnt++)
                {
                    pntoffset = noffset+npnt;
                    NekDouble weight = tracBwdWeightJump[pntoffset];

                    tmpMatinn       = ElmtJac[ntrace][npnt];
                    tmpMatinnData   = tmpMatinn->GetPtr(); 
                    // tmpMatout       = BJac->GetBlock(pntoffset,pntoffset);
                    tmpMatoutData = 
                        BJac->GetBlock(pntoffset,pntoffset)->GetPtr(); 
                    Vmath::Smul(nConvectiveFields*nConvectiveFields,weight,
                            tmpMatinnData,1,tmpMatinnData,1);
                    for(int j = 0; j< nConvectiveFields;j++)
                    {
                        for(int i=0;i<nConvectiveFields;i++)
                        {
                            tmpMatoutData[(j*m_spacedim+nDervDir)
                                *nConvectiveFields+i] =
                                DataType(tmpMatinnData[j*nConvectiveFields+i]);
                        }
                        // Vmath::Vcopy(nConvectiveFields,&tmpMatinnData[j*nConvectiveFields],1,
                        //             &tmpMatoutData[(j*m_spacedim+nDervDir)*nConvectiveFields],1);
                    }

                }
            }
        }
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::NumCalRiemFluxJac(
        const int                                         nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &AdvVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        const TensorOfArray3D<NekDouble>                  &qfield,
        const NekDouble                                   &time,
        const Array<OneD, Array<OneD, NekDouble> >        &Fwd,
        const Array<OneD, Array<OneD, NekDouble> >        &Bwd,
        TypeNekBlkMatSharedPtr                            &FJac,
        TypeNekBlkMatSharedPtr                            &BJac,
        TensorOfArray5D<DataType>                         &TraceIPSymJacArray)
    {
        const NekDouble     PenaltyFactor2  =   0.0;
        int nvariables  = nConvectiveFields;
        int npoints     = GetNpoints();
        // int nPts        = npoints;
        int nTracePts   = GetTraceTotPoints();
        int nDim        = m_spacedim;

        Array<OneD, int > nonZeroIndex;

        Array<OneD,       Array<OneD, NekDouble> >  tmpinarry(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            tmpinarry[i]    =    Array<OneD, NekDouble>(npoints,0.0);
            Vmath::Vcopy(npoints, inarray[i],1,tmpinarry[i],1);
        }

        // DmuDT of artificial diffusion is neglected
        // TODO: to consider the Jacobian of AV seperately
        Array<OneD, NekDouble> muvar        =   NullNekDouble1DArray;
        Array<OneD, NekDouble> MuVarTrace   =   NullNekDouble1DArray;
        if (m_shockCaptureType != "Off" && m_shockCaptureType != "Physical")
        {
            MuVarTrace  =   Array<OneD, NekDouble>(nTracePts, 0.0);
            muvar       =   Array<OneD, NekDouble>(npoints, 0.0);
            m_diffusion->GetAVmu(fields,inarray,muvar,MuVarTrace);
            muvar       =   NullNekDouble1DArray;
        }

        Array<OneD, Array<OneD, NekDouble> > numflux(nvariables);
        for(int i = 0; i < nvariables; ++i)
        {
            numflux[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        }

        const MultiRegions::AssemblyMapDGSharedPtr TraceMap=
            fields[0]->GetTraceMap();
        TensorOfArray3D<NekDouble>    qBwd(nDim);
        TensorOfArray3D<NekDouble>    qFwd(nDim);
        if(m_DEBUG_VISCOUS_JAC_MAT)
        {
            for (int nd = 0; nd < nDim; ++nd)
            {
                qBwd[nd] = 
                    Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                qFwd[nd] = 
                    Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    qBwd[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);
                    qFwd[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);

                    fields[i]->GetFwdBwdTracePhysDeriv_serial(nd,qfield[nd][i], 
                        qFwd[nd][i], qBwd[nd][i]);
                    TraceMap->UniversalTraceAssemble(qBwd[nd][i]);
                    TraceMap->UniversalTraceAssemble(qFwd[nd][i]);
                }
            }
        }

        CalTraceNumericalFlux(nConvectiveFields,nDim,npoints,nTracePts,
            PenaltyFactor2,fields,AdvVel,inarray,time,qfield,Fwd,Bwd,qFwd,qBwd,
            MuVarTrace,nonZeroIndex,numflux);
#ifdef CFS_DEBUGMODE
        TensorOfArray3D<NekDouble> IPSymmFlux(m_spacedim);
        TensorOfArray3D<NekDouble> IPFluxPlus(m_spacedim);
        if(m_DebugIPSymmFluxJacSwitch)
        {
            for(int nd = 0; nd < m_spacedim; ++nd)
            {
                IPSymmFlux[nd] = 
                    Array<OneD, Array<OneD, NekDouble> >(nvariables);
                IPFluxPlus[nd] = 
                    Array<OneD, Array<OneD, NekDouble> >(nvariables);
                for(int i = 0; i < nvariables; ++i)
                {
                    IPSymmFlux[nd][i] = Array<OneD, NekDouble>(nTracePts, 0.0);
                    IPFluxPlus[nd][i] = Array<OneD, NekDouble>(nTracePts, 0.0);
                }

            }
            
            CalTraceIPSymmFlux(nConvectiveFields,nTracePts,fields,inarray,time,
                qfield,Fwd,Bwd,MuVarTrace,nonZeroIndex,IPSymmFlux);
        }
#endif

        int nFields = nvariables;
        // int nPts    = nTracePts;
        Array<OneD, Array<OneD, NekDouble> >  plusFwd(nFields),plusBwd(nFields);
        Array<OneD, Array<OneD, NekDouble> >  Jacvect(nFields);
        Array<OneD, Array<OneD, NekDouble> >  FwdBnd(nFields);
        Array<OneD, Array<OneD, NekDouble> >  plusflux(nFields);
        Array<OneD, Array<OneD, NekDouble> >  minsflux(nFields);
        for(int i = 0; i < nFields; i++)
        {
            Jacvect[i]      =    Array<OneD, NekDouble>(nTracePts,0.0);
            plusFwd[i]      =    Array<OneD, NekDouble>(nTracePts);
            plusBwd[i]      =    Array<OneD, NekDouble>(nTracePts);
            plusflux[i]     =    Array<OneD, NekDouble>(nTracePts,0.0);
            minsflux[i]     =    Array<OneD, NekDouble>(nTracePts,0.0);
            FwdBnd[i]       =    Array<OneD, NekDouble>(nTracePts,0.0);
            Vmath::Vcopy(nTracePts, Fwd[i],1,plusFwd[i],1);
            Vmath::Vcopy(nTracePts, Bwd[i],1,plusBwd[i],1);
        }

        Array<OneD,       Array<OneD, NekDouble> >  minsFwd,minsBwd;
        if(m_centralDiffTracJac)
        {
            minsFwd  = Array<OneD,       Array<OneD, NekDouble> > (nFields);
            minsBwd  = Array<OneD,       Array<OneD, NekDouble> > (nFields);
            for(int i = 0; i < nFields; i++)
            {
                minsFwd[i]      =    Array<OneD, NekDouble>(nTracePts,0.0);
                minsBwd[i]      =    Array<OneD, NekDouble>(nTracePts,0.0);
                Vmath::Vcopy(nTracePts, Fwd[i],1,minsFwd[i],1);
                Vmath::Vcopy(nTracePts, Bwd[i],1,minsBwd[i],1);
            }
        }
        else
        {
            CalTraceNumericalFlux(nConvectiveFields,nDim,npoints,nTracePts,
                PenaltyFactor2,fields,AdvVel,inarray,time,qfield,Fwd,Bwd,qFwd,
                qBwd,MuVarTrace,nonZeroIndex,minsflux);
        }

        // NekDouble eps   =   1.0E-6;
        NekDouble eps   =   1.0E-6;
        
        Array<OneD, DataType> tmpMatData;
        // Fwd Jacobian
        for(int i = 0; i < nFields; i++)
        {
            NekDouble epsvar = eps*m_magnitdEstimat[i];
            NekDouble oepsvar   =   1.0/epsvar;
            Vmath::Sadd(nTracePts,epsvar,Fwd[i],1,plusFwd[i],1);

            if (m_bndConds.size())
            {
                for(int i = 0; i < nFields; i++)
                {
                    Vmath::Vcopy(nTracePts, plusFwd[i],1,FwdBnd[i],1);
                }
                // Loop over user-defined boundary conditions
                for (auto &x : m_bndConds)
                {
                    x->Apply(FwdBnd, tmpinarry, time);
                }
            }

            for(int j = 0; j < nFields; j++)
            {
                m_fields[j]->FillBwdWITHBound(plusFwd[j], plusBwd[j]);
            }

            CalTraceNumericalFlux(nConvectiveFields,nDim,npoints,nTracePts,
                PenaltyFactor2,fields,AdvVel,inarray,time,qfield,plusFwd,
                plusBwd,qFwd,qBwd,MuVarTrace,nonZeroIndex,plusflux);
            
            if(m_centralDiffTracJac)
            {
                Vmath::Sadd(nTracePts, -epsvar,Fwd[i],1,minsFwd[i],1);
                
                if (m_bndConds.size())
                {
                    for(int i = 0; i < nFields; i++)
                    {
                        Vmath::Vcopy(nTracePts, minsFwd[i],1,FwdBnd[i],1);
                    }
                    // Loop over user-defined boundary conditions
                    for (auto &x : m_bndConds)
                    {
                        x->Apply(FwdBnd, tmpinarry, time);
                    }
                }

                for(int j = 0; j < nFields; j++)
                {
                    m_fields[j]->FillBwdWITHBound(minsFwd[j], minsBwd[j]);
                }

                CalTraceNumericalFlux(nConvectiveFields,nDim,npoints,nTracePts,
                    PenaltyFactor2,fields,AdvVel,inarray,time,qfield,minsFwd,
                    minsBwd,qFwd,qBwd,MuVarTrace,nonZeroIndex,minsflux);

                oepsvar   *=   0.5;
            }
            
            for (int n = 0; n < nFields; n++)
            {
                Vmath::Vsub(nTracePts,&plusflux[n][0],1,&minsflux[n][0],1,
                        &Jacvect[n][0],1);
                Vmath::Smul(nTracePts, oepsvar ,&Jacvect[n][0],1,
                        &Jacvect[n][0],1);
            }
            for(int j = 0; j < nTracePts; j++)
            {
                tmpMatData  =   FJac->GetBlock(j,j)->GetPtr();
                for (int n = 0; n < nFields; n++)
                {
                    tmpMatData[n+i*nFields] = DataType(Jacvect[n][j]);
                }
            }
#ifdef CFS_DEBUGMODE
            if(m_DebugIPSymmFluxJacSwitch)
            {
                CalTraceIPSymmFlux(nConvectiveFields,nTracePts,fields,inarray,
                    time,qfield,plusFwd,plusBwd,MuVarTrace,nonZeroIndex,
                    IPFluxPlus);
                for(int nd = 0; nd < m_spacedim; ++nd)
                {
                    for (int n = 0; n < nFields; n++)
                    {
                        Vmath::Vsub(nTracePts,&IPFluxPlus[nd][n][0],1,
                                &IPSymmFlux[nd][n][0],1,&Jacvect[n][0],1);
                        Vmath::Smul(nTracePts, oepsvar ,&Jacvect[n][0],1,
                                &Jacvect[n][0],1);
                        for(int np=0;np<nTracePts;np++)
                        {
                            TraceIPSymJacArray[0][nd][n][i][np] = 
                                DataType( Jacvect[n][np] );
                        }
                    }
                }
            }
#endif

            Vmath::Vcopy(nTracePts, Fwd[i],1,plusFwd[i],1);
            if(m_centralDiffTracJac)
            {
                Vmath::Vcopy(nTracePts, Fwd[i],1,minsFwd[i],1);
            }
        }

        // Reset the boundary conditions
        if (m_bndConds.size())
        {
            for(int i = 0; i < nFields; i++)
            {
                Vmath::Vcopy(nTracePts, Fwd[i],1,FwdBnd[i],1);
            }
            // Loop over user-defined boundary conditions
            for (auto &x : m_bndConds)
            {
                x->Apply(FwdBnd, tmpinarry, time);
            }
        }

        for(int i = 0; i < nFields; i++)
        {
            Vmath::Vcopy(nTracePts, Bwd[i],1,plusBwd[i],1);
        }
        if(m_centralDiffTracJac)
        {
            for(int i = 0; i < nFields; i++)
            {
                Vmath::Vcopy(nTracePts, Bwd[i],1,minsBwd[i],1);
            }
        }

        for(int i = 0; i < nFields; i++)
        {
            NekDouble epsvar    = eps*m_magnitdEstimat[i];
            NekDouble oepsvar   =   1.0/epsvar;

            Vmath::Sadd(nTracePts,epsvar,Bwd[i],1,plusBwd[i],1);

            for(int j = 0; j < nFields; j++)
            {
                m_fields[j]->FillBwdWITHBound(Fwd[j], plusBwd[j]);
            }

            CalTraceNumericalFlux(nConvectiveFields,nDim,npoints,nTracePts,
                PenaltyFactor2,fields,AdvVel,inarray,time,qfield,Fwd,plusBwd,
                qFwd,qBwd,MuVarTrace,nonZeroIndex,plusflux);
            
            if(m_centralDiffTracJac)
            {
                Vmath::Sadd(nTracePts,-epsvar,Bwd[i],1,minsBwd[i],1);

                for(int j = 0; j < nFields; j++)
                {
                    m_fields[j]->FillBwdWITHBound(Fwd[j], minsBwd[j]);
                }

                CalTraceNumericalFlux(nConvectiveFields,nDim,npoints,nTracePts,
                    PenaltyFactor2,fields,AdvVel,inarray,time,qfield,Fwd,
                    minsBwd,qFwd,qBwd,MuVarTrace,nonZeroIndex,minsflux);

                oepsvar   *=   0.5;
            }

            for (int n = 0; n < nFields; n++)
            {
                Vmath::Vsub(nTracePts,&plusflux[n][0],1,&minsflux[n][0],1,
                        &Jacvect[n][0],1);
                Vmath::Smul(nTracePts, oepsvar ,&Jacvect[n][0],1,
                        &Jacvect[n][0],1);
            }
            for(int j = 0; j < nTracePts; j++)
            {
                tmpMatData  =   BJac->GetBlock(j,j)->GetPtr();
                for (int n = 0; n < nFields; n++)
                {
                    tmpMatData[n+i*nFields] = DataType(Jacvect[n][j]);
                }
            }
#ifdef CFS_DEBUGMODE
            if(m_DebugIPSymmFluxJacSwitch)
            {
                CalTraceIPSymmFlux(nConvectiveFields,nTracePts,fields,inarray,
                    time,qfield,Fwd,plusBwd,MuVarTrace,nonZeroIndex,IPFluxPlus);
                for(int nd = 0; nd < m_spacedim; ++nd)
                {
                    for (int n = 0; n < nFields; n++)
                    {
                        Vmath::Vsub(nTracePts,&IPFluxPlus[nd][n][0],1,
                                    &IPSymmFlux[nd][n][0],1,&Jacvect[n][0],1);
                        Vmath::Smul(nTracePts, oepsvar ,&Jacvect[n][0],1,
                                    &Jacvect[n][0],1);
                        for(int np=0;np<nTracePts;np++)
                        {
                            TraceIPSymJacArray[1][nd][n][i][np] = 
                                DataType( Jacvect[n][np] );
                        }
                    }
                }
            }
#endif
            
            Vmath::Vcopy(nTracePts, Bwd[i],1,plusBwd[i],1);
            if(m_centralDiffTracJac)
            {
                Vmath::Vcopy(nTracePts, Bwd[i],1,minsBwd[i],1);
            }
        }
    }

    void CompressibleFlowSystem::CalTraceNumericalFlux(
            const int                                         nConvectiveFields,
            const int                                         nDim,
            const int                                         nPts,
            const int                                         nTracePts,
            const NekDouble                                   PenaltyFactor2,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &AdvVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            const NekDouble                                   time,
            const TensorOfArray3D<NekDouble>                  &qfield,
            const Array<OneD, Array<OneD, NekDouble> >        &vFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &vBwd,
            const TensorOfArray3D<NekDouble>                  &qFwd,
            const TensorOfArray3D<NekDouble>                  &qBwd,
            const Array<OneD, NekDouble >                     &MuVarTrace,
            Array<OneD, int >                                 &nonZeroIndex,
            Array<OneD, Array<OneD, NekDouble> >              &traceflux)
    {
        if (m_DEBUG_ADVECTION_JAC_MAT)
        {
            m_advObject->AdvectTraceFlux(nConvectiveFields, m_fields, AdvVel,
                inarray, traceflux, m_BndEvaluateTime,vFwd, vBwd);
        }
        else
        {
            for (int i = 0; i < nConvectiveFields; i++)
            {
                traceflux[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
            }
        }
        
        if(m_DEBUG_VISCOUS_JAC_MAT)
        {
            Array<OneD, Array<OneD, NekDouble > > visflux(nConvectiveFields);
            for(int i = 0; i < nConvectiveFields; i++)
            {
                visflux[i]  =    Array<OneD, NekDouble>(nTracePts,0.0);
            }

            m_diffusion->DiffuseTraceFlux(nConvectiveFields, fields, inarray,
                qfield,visflux,vFwd, vBwd,qFwd,qBwd,MuVarTrace,nonZeroIndex);
            for(int i = 0; i < nConvectiveFields; i++)
            {
                Vmath::Vsub(nTracePts,traceflux[i],1,visflux[i],1,
                            traceflux[i],1);
            }
        }
    }

    void CompressibleFlowSystem::CalTraceIPSymmFlux(
            const int                                         nConvectiveFields,
            const int                                         nTracePts,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            const NekDouble                                   time,
            const TensorOfArray3D<NekDouble>                  &qfield,
            const Array<OneD, Array<OneD, NekDouble> >        &vFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &vBwd,
            const Array<OneD, NekDouble >                     &MuVarTrace,
            Array<OneD, int >                                 &nonZeroIndex,
            TensorOfArray3D<NekDouble>                        &traceflux)
    {
        if(m_DEBUG_VISCOUS_JAC_MAT)
        {
            TensorOfArray3D<NekDouble>  VolumeFlux;
            m_diffusion->DiffuseTraceSymmFlux(nConvectiveFields,fields,inarray,
                    qfield,VolumeFlux, traceflux,vFwd,vBwd,nonZeroIndex);
            for(int nd = 0; nd < m_spacedim; nd++)
            {
                for(int i = 0; i < nConvectiveFields; i++)
                {
                    Vmath::Ssub(nTracePts,0.0,traceflux[nd][i],1,
                            traceflux[nd][i],1);
                }
            }
            
        }
    }

    void CompressibleFlowSystem::DoImplicitSolve_phy2coeff(
        const TensorOfArray2D<NekDouble>        &inpnts,
        Array<OneD, Array<OneD, NekDouble> >    &outpnt,
        const NekDouble                         time,
        const NekDouble                         lambda)
    {

        unsigned int nvariables  = inpnts.num_elements();
        unsigned int ncoeffs     = m_fields[0]->GetNcoeffs();

        Array<OneD, Array<OneD, NekDouble> > inarray(nvariables);
        Array<OneD, Array<OneD, NekDouble> > out(nvariables);

        for(int i = 0; i < nvariables; i++)
        {
            inarray[i]  =  Array<OneD, NekDouble> (ncoeffs,0.0);
            out[i]      =  Array<OneD, NekDouble> (ncoeffs,0.0);
            m_fields[i]->FwdTrans(inpnts[i],inarray[i]);
        }

        DoImplicitSolve_coeff(inpnts,inarray,out,time,lambda);

        for(int i = 0; i < nvariables; i++)
        {
            m_fields[i]->BwdTrans(out[i],outpnt[i]);
        }
    }

    void CompressibleFlowSystem::DoImplicitSolve_coeff(
        const TensorOfArray2D<NekDouble>        &inpnts,
        const TensorOfArray2D<NekDouble>        &inarray,
        Array<OneD, Array<OneD, NekDouble> >    &out,
        const NekDouble                         time,
        const NekDouble                         lambda)
    {
        m_TimeIntegtSol_n   = inarray;
        m_TimeIntegtSol_k   = out;
        m_BndEvaluateTime   = time;
        m_TimeIntegLambda   = lambda;
        const unsigned int MaxNonlinIte =   m_MaxNonlinIte;
        unsigned int nvariables     = inarray.num_elements();
        unsigned int npoints        = inarray[0].num_elements();
        unsigned int ntotal         = nvariables*npoints;
        unsigned int nElements      = m_fields[0]->GetExpSize();

        LibUtilities::CommSharedPtr v_Comm = 
            m_fields[0]->GetComm()->GetRowComm();

        bool l_verbose      = m_session->DefinesCmdLineArgument("verbose");
        bool l_root=false;
        if(0==v_Comm->GetRank())
        {
            l_root =true;
        }
        unsigned int ntotalGlobal     = ntotal;
        v_Comm->AllReduce(ntotalGlobal, Nektar::LibUtilities::ReduceSum);
        unsigned int ntotalDOF        = ntotalGlobal/nvariables;
        NekDouble ototalGlobal  = 1.0/ntotalGlobal;
        NekDouble ototalDOF     = 1.0/ntotalDOF;

        UpdateSoltnRefNorms(inarray, ototalDOF);

        Array<OneD, NekDouble> NonlinSysRes_1D(ntotal,0.0),
            sol_k_1D(ntotal,0.0),dsol_1D(ntotal,0.0);
        Array<OneD, Array<OneD, NekDouble> > 
            NonlinSysRes(nvariables),dsol(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            int offset = i*npoints;
            NonlinSysRes[i] =  NonlinSysRes_1D + offset;
            dsol[i]         =  dsol_1D + offset;
        }

        m_SysEquatResid_k = NonlinSysRes;
        //sol_k = out;
        for(int i = 0; i < nvariables; i++)
        {
            Vmath::Vcopy(npoints,inarray[i],1,m_TimeIntegtSol_k[i],1);
        }

        if(m_cflLocTimestep>0.0)
        {
            Array<OneD, NekDouble> tstep (nElements, 0.0);
            NekDouble tmpcfl    = m_cflSafetyFactor;
            m_cflSafetyFactor   = m_cflLocTimestep;
            GetElmtTimeStep(inpnts, tstep);
            m_cflSafetyFactor   = tmpcfl;

            m_locTimeStep = tstep;
        }

//Debug
//         NonlinSysEvaluator_coeff(m_TimeIntegtSol_k,m_SysEquatResid_k);
// Fill2DArrayOfBlkDiagonalMat(m_PrecMatVars,0.0);
//         DebugNumCalJac_coeff(m_PrecMatVars);

// Cout2DArrayBlkMat(m_PrecMatVars);

// ASSERTL0(false, "debugstop");
        // ElmtVarInvMtrx_coeff(m_PrecMatVars);
        // TODO: auto precondition recompute. use random vector&L1 norm relative error.
        if(m_CalcuPrecMatFlag||(m_TimeIntegLambdaPrcMat!=m_TimeIntegLambda))
        {
            if(0 < m_JFNKPrecondStep)
            {
                CalPrecMat(inpnts, m_BndEvaluateTime, m_TimeIntegLambda);
            }
        }

        bool converged       = false;
        bool steady_state    = false;
        int nwidthcolm = 10;
        int NtotDoOdeRHS = 0;
        int NtotGMRESIts = 0;
        NekDouble resnorm0 = 0.0;
        NekDouble resnorm  = 0.0;
        NekDouble resmaxm  = 0.0;
        NekDouble resratio = 1.0;
        NekDouble LinSysTol = 0.0;
        int NttlNonlinIte    = 0;
        for (int k = 0; k < MaxNonlinIte; k++)
        {
            NekDouble GMRESRelativeIteTol;
            NonlinSysEvaluator_coeff(m_TimeIntegtSol_k,m_SysEquatResid_k);
        // for(int i = 0; i < nvariables; i++)
        // {

        //     int nwidthcolm = 20;
        //     int nphspnt = m_SysEquatResid_k[nvariables-1].num_elements();
        //     for(int j = 0; j < nphspnt; j++)
        //     {
        //         cout <<right<<scientific<<setw(nwidthcolm)<<setprecision(nwidthcolm-6)
        //             <<"RHS: ("<<i<<" "<<j<<") "<<m_SysEquatResid_k[i][j]
        //             << endl;
        //     }
        // }
            if(2==m_LiniearizationMethod||3==m_LiniearizationMethod)
            {
                m_updateMatFreeJacFlag = true;
            }
        // int ncoefftmp = m_TimeIntegtSol_k[nvariables-1].num_elements();
        // for(int i = 0; i < nvariables; i++)
        // {
        //     int noffset = ncoefftmp*i;
        //     Vmath::Vcopy(ncoefftmp,&m_TimeIntegtSol_k[i][0],1,&dsol_1D[0]+noffset,1);
        // }

        // MatrixMultiply_JacobianFree_coeff(dsol_1D,sol_k_1D);
        // ofstream outfileJF;
        // outfileJF.open("MatrixMultiply_JacobianFree_coeff.dat");
        // for(int i = 0; i < nvariables; i++)
        // {

        //     int nwidthcolm = 20;
        //     int nphspnt = m_SysEquatResid_k[nvariables-1].num_elements();
        //     for(int j = 0; j < nphspnt; j++)
        //     {
        //         outfileJF <<right<<scientific<<setw(nwidthcolm)<<setprecision(nwidthcolm-6)
        //             <<"MF :("<<i<<" "<<j<<") "<<sol_k_1D[i*ncoefftmp+j]
        //             << endl;
        //     }
        // }
        // outfileJF.close();

        // for(int i = 0; i < nvariables; i++)
        // {
        //     int noffset = ncoefftmp*i;
        //     Vmath::Vcopy(ncoefftmp,&m_TimeIntegtSol_k[i][0],1,&dsol_1D[0]+noffset,1);
        // }
        // MatrixMultiply_MatrixFree_coeff(dsol_1D,sol_k_1D);
        // ofstream outfileMF;
        // outfileMF.open("MatrixMultiply_MatrixFree_coeff.dat");
        // for(int i = 0; i < nvariables; i++)
        // {

        //     int nwidthcolm = 20;
        //     int nphspnt = m_SysEquatResid_k[nvariables-1].num_elements();
        //     for(int j = 0; j < nphspnt; j++)
        //     {
        //         outfileMF <<right<<scientific<<setw(nwidthcolm)<<setprecision(nwidthcolm-6)
        //             <<"MF :("<<i<<" "<<j<<") "<<sol_k_1D[i*ncoefftmp+j]
        //             << endl;
              //     }
        // }
        // outfileMF.close();
        // ASSERTL0(false,"debug stop");
            // if(2==m_LiniearizationMethod)
            // {
            //     CalphysDeriv(m_TimeIntegtSol_k,m_qfield);

            //     Array<OneD, Array<OneD, NekDouble> > pnts(nvariables);
            //     int nphspnt = GetTotPoints();
            //     for(int i = 0; i < nvariables; i++)
            //     {
            //         pnts[i]   =   Array<OneD, NekDouble>(nphspnt,0.0);
            //         m_fields[i]->BwdTrans(m_TimeIntegtSol_k[i], pnts[i]);
            //     }
            //     GetTraceJac(pnts,m_qfield,m_TraceJac,m_TraceJacDeriv,m_TraceJacDerivSign);
        
            //     for(int j = 0; j< nvariables; j++)
            //     {
            //         Vmath::Vcopy(nphspnt, inarray[j],1,m_MatrixFreeRefFields[j],1);
            //     }
            // }

          
            // NonlinSysEvaluator_coeff_out(m_TimeIntegtSol_k,m_SysEquatResid_k);
            NtotDoOdeRHS++;
            
            // NonlinSysRes_1D and m_SysEquatResid_k share the same storage
            resnorm = Vmath::Dot(ntotal,NonlinSysRes_1D,NonlinSysRes_1D);
            v_Comm->AllReduce(resnorm, Nektar::LibUtilities::ReduceSum);
            NekDouble resnormOld = 0.0;
            
            converged = NewtonStopCriteria(NonlinSysRes_1D, k, ntotalDOF, 
                        resnorm, resnorm0, resmaxm, resratio);
            
            CalcInexactNewtonForcing(k, resnormOld, resnorm, 
                GMRESRelativeIteTol);
            resnormOld = resnorm;
            
            if (converged)
            {
                break;
            }
            LinSysTol = GMRESRelativeIteTol*sqrt(resnorm);
            int ntmpGMRESIts =  m_linsol->SolveLinearSystem(ntotal,
                NonlinSysRes_1D,dsol_1D,0,LinSysTol);
            NtotDoOdeRHS  +=  ntmpGMRESIts;
            NtotGMRESIts  +=  ntmpGMRESIts;

            for(int i = 0; i < nvariables; i++)
            {
                Vmath::Vsub(npoints,m_TimeIntegtSol_k[i],1,dsol[i],1,
                            m_TimeIntegtSol_k[i],1);
            }
            NttlNonlinIte++;
        }
        if (m_PrcdMatFreezNumb>0)
        {
            m_CalcuPrecMatFlag = false;
        }

        if((l_verbose||(!converged))&&l_root)
        {
            WARNINGL0(converged,
             "     # Nonlinear system solver not converge in DoImplicitSolve ");
            cout <<right<<scientific<<setw(nwidthcolm)
                <<setprecision(nwidthcolm-6)
                <<"     * Newton-Its converged (RES="
                << sqrt(resnorm)<<" Res/Q="<< sqrt(resnorm/m_inArrayNorm)
                <<" ResMax/QPerDOF="<< resmaxm*sqrt(ntotalDOF/m_inArrayNorm)
                <<" Res/(DtRHS): "<<sqrt(resratio)
                <<" with "<<setw(3)<<NttlNonlinIte<<" Non-Its"<<" and "<<setw(4)
                <<NtotDoOdeRHS<<" Lin-Its)"<<endl;
        }

        m_TotNewtonIts  +=  NttlNonlinIte;
        m_TotGMRESIts   +=  NtotGMRESIts;
        m_TotOdeRHS     +=  NtotDoOdeRHS;
        m_TotImpStages++;

        m_TotLinItePerStep += NtotDoOdeRHS;
        m_StagesPerStep++;
        if(!converged)
        {
            m_TotLinItePerStep +=   m_maxLinItePerNewton*100;
        }
    }

    void CompressibleFlowSystem::UpdateSoltnRefNorms(
        const Array<OneD, Array<OneD, NekDouble>>    &inarray,
        const NekDouble                 &ototalDOF)
    {
        int nvariables = m_fields.num_elements();
        int npoints = GetNcoeffs();
        m_inArrayNorm = 0.0;

        LibUtilities::CommSharedPtr v_Comm = 
            m_fields[0]->GetComm()->GetRowComm();
        // estimate the magnitude of each flow variables
        m_magnitdEstimat = Array<OneD, NekDouble>  (nvariables,0.0);

        for(int i = 0; i < nvariables; i++)
        {
            m_magnitdEstimat[i] = Vmath::Dot(npoints,inarray[i],inarray[i]);
        }
        v_Comm->AllReduce(m_magnitdEstimat, Nektar::LibUtilities::ReduceSum);

        for(int i = 0; i < nvariables; i++)
        {
            m_inArrayNorm += m_magnitdEstimat[i];
        }

        for(int i = 2; i < nvariables-1; i++)
        {
            m_magnitdEstimat[1]   +=   m_magnitdEstimat[i] ;
        }
        for(int i = 2; i < nvariables-1; i++)
        {
            m_magnitdEstimat[i]   =   m_magnitdEstimat[1] ;
        }

        for(int i = 0; i < nvariables; i++)
        {
            m_magnitdEstimat[i] = sqrt(m_magnitdEstimat[i]*ototalDOF);
        }
        
    }

    void CompressibleFlowSystem::CalPrecMat(
        const TensorOfArray2D<NekDouble>    &inpnts,
        const NekDouble                     time,
        const NekDouble                     lambda)
    {
        int nvariables = inpnts.num_elements();
#ifdef CFS_DEBUGMODE
        if(m_DebugNumJacBSOR)
        {
            NonlinSysEvaluator_coeff(m_TimeIntegtSol_k,m_SysEquatResid_k);
            Fill2DArrayOfBlkDiagonalMat(m_PrecMatVars,0.0);
            DebugNumCalJac_coeff(m_PrecMatVars,m_PrecMatVarsOffDiag);

            NekDouble zero=0.0;
            ElmtVarInvMtrx(m_PrecMatVars,zero);
        }
        else
        {
#endif
            int nphspnt = inpnts[0].num_elements();
            Array<OneD, Array<OneD, NekDouble> > intmp(nvariables);
            for(int i = 0; i < nvariables; i++)
            {
                intmp[i]    =   Array<OneD, NekDouble>(nphspnt,0.0);
            }

            DoOdeProjection(inpnts,intmp,time);
            if(m_flagPrecMatVarsSingle)
            {
                if(m_flagPrecondCacheOptmis)
                {
                    GetpreconditionerNSBlkDiag_coeff(intmp,m_PrecMatVarsSingle,
                        m_PrecMatSingle, m_TraceJacSingle,
                        m_TraceJacDerivSingle,m_TraceJacDerivSignSingle,
                        m_TraceJacArraySingle,m_TraceJacDerivArraySingle,
                        m_TraceIPSymJacArraySingle,
                        m_StdSMatDataDBB,m_StdSMatDataDBDB);
                }
                else
                {
                    GetpreconditionerNSBlkDiag_coeff(intmp,m_PrecMatVarsSingle,
                        m_TraceJacSingle,m_TraceJacDerivSingle,
                        m_TraceJacDerivSignSingle, m_TraceIPSymJacArraySingle,
                        m_StdSMatDataDBB, m_StdSMatDataDBDB);
                }
            }
            else
            {
                if(m_flagPrecondCacheOptmis)
                {
                    GetpreconditionerNSBlkDiag_coeff(intmp,m_PrecMatVars,
                        m_PrecMat,m_TraceJac,m_TraceJacDeriv,
                        m_TraceJacDerivSign, m_TraceJacArray,
                        m_TraceJacDerivArray,m_TraceIPSymJacArray,
                        m_StdDMatDataDBB,m_StdDMatDataDBDB);
                }
                else
                {
                    GetpreconditionerNSBlkDiag_coeff(intmp,m_PrecMatVars,
                        m_TraceJac,m_TraceJacDeriv,m_TraceJacDerivSign,
                        m_TraceIPSymJacArray, m_StdDMatDataDBB,
                        m_StdDMatDataDBDB);
                }
            }
        
            m_TimeIntegLambdaPrcMat = lambda;

            // to free the storage
            for(int i = 0; i < nvariables; i++)
            {
                intmp[i]    =   NullNekDouble1DArray;
            }

#ifdef CFS_DEBUGMODE
            if(m_DebugNumJacMatSwitch)
            {
                cout <<" m_DebugNumJacMatSwitch "<<endl;
                NonlinSysEvaluator_coeff(m_TimeIntegtSol_k,m_SysEquatResid_k);
                Fill2DArrayOfBlkDiagonalMat(m_PrecMatVars,0.0);
                DebugNumCalJac_coeff(m_PrecMatVars);

                NekDouble zero=0.0;
                ElmtVarInvMtrx(m_PrecMatVars,zero);

                for(int m = 0; m < nvariables; m++)
                {
                    for(int n = 0; n < nvariables; n++)
                    {
                        m_fields[0]->RightMultiplyByElmtInvMassOnDiag(
                            m_PrecMatVars[m][n], m_PrecMatVars[m][n]);
                    }
                }
            }
            if(m_DebugOutputJacMatSwitch)
            {
                Cout2DArrayBlkMat(m_PrecMatVars);
            }
        }
#endif
    }

    bool CompressibleFlowSystem::NewtonStopCriteria(
        const Array<OneD, NekDouble>    &NonlinSysRes_1D,
        const int                       &k,
        const int                       &totalDOF,
        const NekDouble                 &resnorm,
        NekDouble                       &resnorm0,
        NekDouble                       &resmaxm,
        NekDouble                       &resratio)
    {
        bool converged = false;
        LibUtilities::CommSharedPtr v_Comm  = m_fields[0]->GetComm()->GetRowComm();
        bool l_verbose      = m_session->DefinesCmdLineArgument("verbose");
        bool l_root=false;
        if(0==v_Comm->GetRank())
        {
            l_root =true;
        }

        NekDouble tolrnc    = m_NewtonAbsoluteIteTol;
        NekDouble tol2      = m_inArrayNorm*tolrnc*tolrnc;

        NekDouble ototalDOF = 1.0/totalDOF;

        NekDouble ratioTol = m_NewtonRelativeIteTol;
        NekDouble tol2Ratio = ratioTol*ratioTol;
        NekDouble tol2Max   = sqrt(m_inArrayNorm*ototalDOF)*ratioTol;

        int nwidthcolm = 8;
        NekDouble ratioSteps = 1.0;
        if(0==k)
        {
            resnorm0 = resnorm;
            resratio = 1.0;

            if(m_Res0PreviousStep<0.0)
            {
                m_Res0PreviousStep = resnorm0;
                ratioSteps         = 1.0;
            }
            else
            {
                ratioSteps         = m_Res0PreviousStep/resnorm0;
                m_Res0PreviousStep = resnorm0;
            }

        }
        else
        {
            resratio = resnorm/resnorm0;
        }
        if (resratio<tol2Ratio||resnorm<tol2)
        {
            resmaxm = 0.0;
            for(int i=0;i<NonlinSysRes_1D.num_elements();i++)
            {
                resmaxm = max(resmaxm,abs(NonlinSysRes_1D[i]));
            }
            v_Comm->AllReduce(resmaxm, Nektar::LibUtilities::ReduceMax);
            if((resmaxm<tol2Max)&&k>0)
            {
            // if(k>0)
            // {
                /// TODO: m_root
                // if(m_cflLocTimestep>0.0)
                // {
                //     m_cflLocTimestep *= pow(ratioSteps,0.35);
                // }
                // m_cflSafetyFactor *= pow(ratioSteps,0.35);

                converged = true;
                if(resratio>tol2Ratio&&l_root)
                {
                    WARNINGL0(true,
                        "     # resratio>tol2Ratio in DoImplicitSolve ");
                    cout <<right<<scientific<<setw(nwidthcolm)
                         <<setprecision(nwidthcolm-6) <<" resratio= "
                         <<resratio<<" tol2Ratio= "<<tol2Ratio<<endl;
                }
            }
        }
            
        return converged;
    }
    void CompressibleFlowSystem::CalcInexactNewtonForcing(
        const int       &k,
        const NekDouble &resnormOld,
        const NekDouble &resnorm,
        NekDouble       &forcing)
    {
        bool l_verbose      = m_session->DefinesCmdLineArgument("verbose");

        LibUtilities::CommSharedPtr v_Comm = 
            m_fields[0]->GetComm()->GetRowComm();
        bool l_root=false;
        if(0==v_Comm->GetRank())
        {
            l_root =true;
        }
        if (0 == k)
        {
            forcing = m_GMRESRelativeIteTol;
        }
        else
        {
            switch(m_adapGMRESTol)
            {
            case 0:
                {
                    forcing = m_GMRESRelativeIteTol;
                    break;
                }
            case 1: 
                {
                    NekDouble tmpForc = m_ForcingGama * 
                        pow((resnorm / resnormOld), m_ForcingAlpha); 
                    NekDouble tmp = m_ForcingGama * 
                        pow(forcing, m_ForcingAlpha);
                    if (tmp > 0.1)
                    {
                        forcing = min(m_GMRESRelativeIteTol, max(tmp, tmpForc));
                    }
                    else
                    {
                        forcing = min(m_GMRESRelativeIteTol, tmpForc);
                    }
                    if (l_verbose&&l_root)
                    {
                        cout << " resnormRatio= " 
                            << resnorm / resnormOld 
                            << " tmpForc= "
                            << tmpForc;
                    }
                    break;
                }
            }
        }

        // if (l_verbose&&l_root)
        // {
        //     cout << " GMRESRelativeIteTol= " << forcing << endl;
        // }
    }

    template<typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::AllocatePrecondBlkDiag_coeff(
        TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
        const int                               &nscale )
    {

        int nvars = m_fields.num_elements();
        int nelmts  = m_fields[0]->GetNumElmts();
        int nrowsVars,ncolsVars;
        int nelmtcoef;
        Array<OneD, unsigned int > nelmtmatdim(nelmts);
        for(int i = 0; i < nelmts; i++)
        {
            nelmtcoef   =   m_fields[0]->GetExp(i)->GetNcoeffs();
            nelmtmatdim[i]  =   nelmtcoef*nscale;
        }

        for(int i = 0; i < nvars; i++)
        {
            for(int j = 0; j < nvars; j++)
            {
                AllocateNekBlkMatDig(gmtxarray[i][j],nelmtmatdim,nelmtmatdim);
            }
        }
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::GetpreconditionerNSBlkDiag_coeff(
        const TensorOfArray2D<NekDouble>        &inarray,
        TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
        TypeNekBlkMatSharedPtr                  &gmtVar,
        Array<OneD, TypeNekBlkMatSharedPtr >    &TraceJac,
        Array<OneD, TypeNekBlkMatSharedPtr >    &TraceJacDeriv,
        Array<OneD, Array<OneD, DataType> >     &TraceJacDerivSign,
        TensorOfArray4D<DataType>               &TraceJacArray,
        TensorOfArray4D<DataType>               &TraceJacDerivArray,
        TensorOfArray5D<DataType>               &TraceIPSymJacArray,
        TensorOfArray4D<DataType>               &StdMatDataDBB,
        TensorOfArray5D<DataType>               &StdMatDataDBDB)
    {
        TensorOfArray3D<NekDouble> qfield;

        if(m_DEBUG_VISCOUS_JAC_MAT)
        {
            CalphysDeriv(inarray,qfield);
        }

        DataType zero =0.0;
        Fill2DArrayOfBlkDiagonalMat(gmtxarray,zero);
#ifdef CFS_DEBUGMODE
        if(2!=m_DebugVolTraceSwitch)
        {
#endif
            AddMatNSBlkDiag_volume(inarray,qfield,gmtxarray,StdMatDataDBB,
                StdMatDataDBDB);
#ifdef CFS_DEBUGMODE
        }
        if(1!=m_DebugVolTraceSwitch)
        {
#endif
            AddMatNSBlkDiag_boundary(inarray,qfield,gmtxarray,TraceJac,
                    TraceJacDeriv,TraceJacDerivSign,TraceIPSymJacArray);
#ifdef CFS_DEBUGMODE
        }
#endif
        MultiplyElmtInvMass_PlusSource(gmtxarray,m_TimeIntegLambda,zero);

        ElmtVarInvMtrx(gmtxarray,gmtVar,zero);
        
        TransTraceJacMatToArray(TraceJac,TraceJacDeriv,TraceJacArray, 
            TraceJacDerivArray);
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::GetpreconditionerNSBlkDiag_coeff(
        const TensorOfArray2D<NekDouble>        &inarray,
        TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
        Array<OneD, TypeNekBlkMatSharedPtr >    &TraceJac,
        Array<OneD, TypeNekBlkMatSharedPtr >    &TraceJacDeriv,
        Array<OneD, Array<OneD, DataType> >     &TraceJacDerivSign,
        TensorOfArray5D<DataType>               &TraceIPSymJacArray,
        TensorOfArray4D<DataType>               &StdMatDataDBB,
        TensorOfArray5D<DataType>               &StdMatDataDBDB)
    {
        if(m_DEBUG_VISCOUS_JAC_MAT)
        {
            CalphysDeriv(inarray,m_qfield);
        }

        DataType zero =0.0;
        Fill2DArrayOfBlkDiagonalMat(gmtxarray,zero);
#ifdef CFS_DEBUGMODE
        if(2!=m_DebugVolTraceSwitch)
        {
#endif
            AddMatNSBlkDiag_volume(inarray,m_qfield,gmtxarray,StdMatDataDBB,
                                    StdMatDataDBDB);
#ifdef CFS_DEBUGMODE
        }
        if(1!=m_DebugVolTraceSwitch)
        {
#endif
            AddMatNSBlkDiag_boundary(inarray,m_qfield,gmtxarray,TraceJac,
                    TraceJacDeriv,TraceJacDerivSign,TraceIPSymJacArray);
#ifdef CFS_DEBUGMODE
        }
#endif
        MultiplyElmtInvMass_PlusSource(gmtxarray,m_TimeIntegLambda,zero);

        ElmtVarInvMtrx(gmtxarray,zero);
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::MultiplyElmtInvMass_PlusSource(
        TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
        const NekDouble                         dtlamda,
        const DataType                          tmpDataType)
    {
        MultiRegions::ExpListSharedPtr explist = m_fields[0];
        std::shared_ptr<LocalRegions::ExpansionVector> pexp = 
            explist->GetExp();
        int ntotElmt            = (*pexp).size();
        int nElmtCoef;
        int nConvectiveFields = m_fields.num_elements();

        NekDouble Negdtlamda    =   -dtlamda;

        Array<OneD, NekDouble> pseudotimefactor(ntotElmt,0.0);
        if(m_cflLocTimestep>0.0)
        {
            NekDouble Neglamda      =   Negdtlamda/m_timestep;
            Vmath::Svtvp(ntotElmt,Neglamda,m_locTimeStep,1,pseudotimefactor,1,
                        pseudotimefactor,1);
        }
        else
        {
            Vmath::Fill(ntotElmt,Negdtlamda,pseudotimefactor,1);
        }

        Array<OneD, DataType>    GMatData;
        for(int m = 0; m < nConvectiveFields; m++)
        {
            for(int n = 0; n < nConvectiveFields; n++)
            {
                for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                {
                    GMatData = gmtxarray[m][n]->GetBlock(nelmt,nelmt)->GetPtr();
                    DataType factor = DataType(pseudotimefactor[nelmt]);

                    Vmath::Smul(GMatData.num_elements(),factor,GMatData,1,
                                GMatData,1);
                }
            }
        }

        DNekMatSharedPtr MassMat;
        Array<OneD,NekDouble> BwdMatData,MassMatData,tmp;
        Array<OneD,DataType> MassMatDataDataType;
        LibUtilities::ShapeType ElmtTypePrevious = LibUtilities::eNoShapeType;

        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            int nelmtcoef  = GetNcoeffs(nelmt);
            int nelmtpnts  = GetTotPoints(nelmt);
            LibUtilities::ShapeType ElmtTypeNow = 
                    explist->GetExp(nelmt)->DetShapeType();

            if (tmp.num_elements()!=nelmtcoef||(ElmtTypeNow!=ElmtTypePrevious)) 
            {
                StdRegions::StdExpansionSharedPtr stdExp;
                stdExp = explist->GetExp(nelmt)->GetStdExp();
                StdRegions::StdMatrixKey  matkey(StdRegions::eBwdTrans,
                                    stdExp->DetShapeType(), *stdExp);
                    
                DNekMatSharedPtr BwdMat =  stdExp->GetStdMatrix(matkey);
                BwdMatData = BwdMat->GetPtr();

                if(nelmtcoef!=tmp.num_elements())
                {
                    tmp = Array<OneD,NekDouble> (nelmtcoef,0.0);
                    MassMat  =  MemoryManager<DNekMat>
                        ::AllocateSharedPtr(nelmtcoef, nelmtcoef, 0.0);
                    MassMatData = MassMat->GetPtr();
                    MassMatDataDataType = 
                        Array<OneD, DataType> (nelmtcoef*nelmtcoef);
                }

                ElmtTypePrevious    = ElmtTypeNow;
            }
            
            for(int np=0; np<nelmtcoef;np++)
            {
                explist->GetExp(nelmt)->IProductWRTBase
                    (BwdMatData+np*nelmtpnts,tmp);
                Vmath::Vcopy(nelmtcoef,&tmp[0],1,
                    &MassMatData[0]+np*nelmtcoef,1);
            }
            for(int i=0;i<MassMatData.num_elements();i++)
            {
                MassMatDataDataType[i]    =   DataType(MassMatData[i]);
            }

            for(int m = 0; m < nConvectiveFields; m++)
            {
                GMatData = gmtxarray[m][m]->GetBlock(nelmt,nelmt)->GetPtr();
                Vmath::Vadd(MassMatData.num_elements(),MassMatDataDataType,1,
                            GMatData,1,GMatData,1);
            }
        }
        return;
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::TransTraceJacMatToArray(
        const Array<OneD, TypeNekBlkMatSharedPtr > &TraceJac,
        const Array<OneD, TypeNekBlkMatSharedPtr > &TraceJacDeriv,
        TensorOfArray4D<DataType>                  &TraceJacArray,
        TensorOfArray4D<DataType>                  &TraceJacDerivArray)
    {
        int nFwdBwd,nDiagBlks,nvar0Jac,nvar1Jac,nDiagBlksVar1Jac;
        int nvar0Drv,nvar1Drv,nDiagBlksVar1Drv;
        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;
        nFwdBwd = TraceJac.num_elements();
        TraceJac[0]->GetBlockSizes(rowSizes,colSizes);
        nDiagBlks   = rowSizes.num_elements();
        nvar0Jac    = rowSizes[1] - rowSizes[0];
        nvar1Jac    = colSizes[1] - colSizes[0];
        nDiagBlksVar1Jac = nDiagBlks*nvar1Jac;

        if(0==TraceJacArray.num_elements())
        {
            TraceJacArray = TensorOfArray4D<DataType> (nFwdBwd);
            for(int nlr=0;nlr<nFwdBwd;nlr++)
            {
                TraceJacArray[nlr] = TensorOfArray3D<DataType> (nvar0Jac);
                for(int m=0;m<nvar0Jac;m++)
                {
                    TraceJacArray[nlr][m] = 
                        Array<OneD,Array<OneD,DataType > >(nvar1Jac);
                    for(int n=0;n<nvar1Jac;n++)
                    {
                        TraceJacArray[nlr][m][n] = 
                            Array<OneD,DataType >(nDiagBlks);
                    }
                }
            }
        }

        for(int nlr=0;nlr<nFwdBwd;nlr++)
        {
            const TypeNekBlkMatSharedPtr tmpMat = TraceJac[nlr];
            TensorOfArray3D<DataType> tmpaa = TraceJacArray[nlr];
            TranSamesizeBlkDiagMatIntoArray(tmpMat,tmpaa);
        }

        if(m_DEBUG_VISCOUS_JAC_MAT&&m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT)
        {
            TraceJacDeriv[0]->GetBlockSizes(rowSizes,colSizes);
            nvar0Drv    = rowSizes[1] - rowSizes[0];
            nvar1Drv    = colSizes[1] - colSizes[0];
            nDiagBlksVar1Drv = nDiagBlks*nvar1Drv;
            if(0==TraceJacDerivArray.num_elements())
            {
                TraceJacDerivArray = TensorOfArray4D<DataType> (nFwdBwd);
                TraceJacDerivArray[0] = TensorOfArray3D<DataType> (nvar0Jac);
                int nlr =0;
                for(int m=0;m<nvar0Drv;m++)
                {
                    TraceJacDerivArray[nlr][m] = 
                        Array<OneD,Array<OneD,DataType > >(nvar1Drv);
                    for(int n=0;n<nvar1Drv;n++)
                    {
                        TraceJacDerivArray[nlr][m][n] = 
                            Array<OneD,DataType >(nDiagBlks);
                    }
                }
                for(int nlr=1;nlr<nFwdBwd;nlr++)
                {
                    TraceJacDerivArray[nlr]  =   TraceJacDerivArray[0];
                }
            }
            int nlr =0;
            TranSamesizeBlkDiagMatIntoArray(TraceJacDeriv[nlr],
                    TraceJacDerivArray[nlr]);
        }
        return;
    }

    void CompressibleFlowSystem::MatrixMultiply_JacobianFree_coeff(
        const  Array<OneD, NekDouble>   &inarray,
        Array<OneD, NekDouble >         &out)
    {
        NekDouble eps = m_JFEps;
        unsigned int ntotalGlobal     = inarray.num_elements();
        NekDouble magninarray = Vmath::Dot(ntotalGlobal,inarray,inarray);
        LibUtilities::CommSharedPtr v_Comm = 
                m_fields[0]->GetComm()->GetRowComm();
        v_Comm->AllReduce(magninarray, Nektar::LibUtilities::ReduceSum);
        // eps *= magnitdEstimatMax;
        eps *= sqrt( (sqrt(m_inArrayNorm) + 1.0)/magninarray);
        NekDouble oeps = 1.0/eps;
        unsigned int nvariables = m_fields.num_elements();
        unsigned int npoints    = GetNcoeffs();
        Array<OneD, NekDouble > tmp;
        Array<OneD,       Array<OneD, NekDouble> > solplus(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > resplus(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            solplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
            resplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
        }

        for (int i = 0; i < nvariables; i++)
        {
            tmp = inarray + i*npoints;
            Vmath::Svtvp(npoints,eps,tmp,1,m_TimeIntegtSol_k[i],1,solplus[i],1);
        }
        
        NonlinSysEvaluator_coeff(solplus,resplus);

        for (int i = 0; i < nvariables; i++)
        {
            tmp = out + i*npoints;
            Vmath::Vsub(npoints,&resplus[i][0],1,&m_SysEquatResid_k[i][0],1,
                        &tmp[0],1);
            Vmath::Smul(npoints, oeps ,&tmp[0],1,&tmp[0],1);
        }
       
        return;
    }

    void CompressibleFlowSystem::MatrixMultiply_JacobianFree_coeff_noSource(
        const  Array<OneD, NekDouble>   &inarray,
        Array<OneD, NekDouble >         &out)
    {
        NekDouble eps = m_JFEps;
        unsigned int ntotalGlobal     = inarray.num_elements();
        NekDouble magninarray = Vmath::Dot(ntotalGlobal,inarray,inarray);
        LibUtilities::CommSharedPtr v_Comm  = m_fields[0]->GetComm()->GetRowComm();
        v_Comm->AllReduce(magninarray, Nektar::LibUtilities::ReduceSum);
        // eps *= magnitdEstimatMax;
        eps *= sqrt( (sqrt(m_inArrayNorm) + 1.0)/magninarray);
        NekDouble oeps = 1.0/eps;
        unsigned int nvariables = m_fields.num_elements();
        unsigned int npoints    = GetNcoeffs();
        Array<OneD, NekDouble > tmp;
        Array<OneD,       Array<OneD, NekDouble> > solplus(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > resplus(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            solplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
            resplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
        }

        for (int i = 0; i < nvariables; i++)
        {
            tmp = inarray + i*npoints;
            Vmath::Svtvp(npoints,eps,tmp,1,m_TimeIntegtSol_k[i],1,solplus[i],1);
        }
        
        NonlinSysEvaluator_coeff_noSource(solplus,resplus);

        for (int i = 0; i < nvariables; i++)
        {
            tmp = out + i*npoints;
            Vmath::Vsub(npoints,&resplus[i][0],1,&m_LowLevelSysEquatResid[i][0],1,&tmp[0],1);
            Vmath::Smul(npoints, oeps ,&tmp[0],1,&tmp[0],1);
        }
       
        return;
    }

    void CompressibleFlowSystem::MatrixMultiply_JacobianFree_coeff_central(
        const  Array<OneD, NekDouble>   &inarray,
        Array<OneD, NekDouble >         &out)
    {
        NekDouble eps = m_JFEps;
        unsigned int ntotalGlobal     = inarray.num_elements();
        NekDouble magninarray = Vmath::Dot(ntotalGlobal,inarray,inarray);
        LibUtilities::CommSharedPtr v_Comm = 
            m_fields[0]->GetComm()->GetRowComm();
        v_Comm->AllReduce(magninarray, Nektar::LibUtilities::ReduceSum);
        // eps *= magnitdEstimatMax;
        eps *= sqrt( (sqrt(m_inArrayNorm) + 1.0)/magninarray);
        NekDouble oeps = 1.0/eps;
        unsigned int nvariables = m_fields.num_elements();
        unsigned int npoints    = GetNcoeffs();
        Array<OneD, NekDouble > tmp;
        Array<OneD,       Array<OneD, NekDouble> > solplus(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > resplus(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > resminus(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            solplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
            resplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
            resminus[i] =  Array<OneD, NekDouble>(npoints,0.0);
        }

        for (int i = 0; i < nvariables; i++)
        {
            tmp = inarray + i*npoints;
            Vmath::Svtvp(npoints,eps,tmp,1,m_TimeIntegtSol_k[i],1,solplus[i],1);
        }
        NonlinSysEvaluator_coeff_noSource(solplus,resplus);

        for (int i = 0; i < nvariables; i++)
        {
            tmp = inarray + i*npoints;
            Vmath::Svtvp(npoints,-eps,tmp,1,m_TimeIntegtSol_k[i],1,
                        solplus[i],1);
        }
        NonlinSysEvaluator_coeff_noSource(solplus,resminus);

        for (int i = 0; i < nvariables; i++)
        {
            tmp = out + i*npoints;
            Vmath::Vsub(npoints,&resplus[i][0],1,&resminus[i][0],1,&tmp[0],1);
            Vmath::Smul(npoints, 0.5*oeps ,&tmp[0],1,&tmp[0],1);
        }
       
        return;
    }

    void CompressibleFlowSystem::MatrixMultiply_JacobianFree_coeff_FourthCentral(
        const  Array<OneD, NekDouble>   &inarray,
        Array<OneD, NekDouble >         &out)
    {
        NekDouble eps = m_JFEps;
        NekDouble magnitdEstimatMax =0.0;
        for(int i = 0; i < m_magnitdEstimat.num_elements(); i++)
        {
            magnitdEstimatMax = max(magnitdEstimatMax,m_magnitdEstimat[i]);
        }
        eps *= magnitdEstimatMax;
        NekDouble oeps = 1.0/eps;
        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int npoints    = m_TimeIntegtSol_n[0].num_elements();
        Array<OneD, NekDouble > tmp;
        Array<OneD,       Array<OneD, NekDouble> > solplus(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > resplus(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > resPPls(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > resminus(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > resMMmns(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            solplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
            resplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
            resPPls[i] =  Array<OneD, NekDouble>(npoints,0.0);
            resminus[i] =  Array<OneD, NekDouble>(npoints,0.0);
            resMMmns[i] =  Array<OneD, NekDouble>(npoints,0.0);
        }

        for (int i = 0; i < nvariables; i++)
        {
            tmp = inarray + i*npoints;
            Vmath::Svtvp(npoints,eps,tmp,1,m_TimeIntegtSol_k[i],1,
                        solplus[i],1);
        }
        NonlinSysEvaluator_coeff(solplus,resplus);

        for (int i = 0; i < nvariables; i++)
        {
            tmp = inarray + i*npoints;
            Vmath::Svtvp(npoints,-eps,tmp,1,m_TimeIntegtSol_k[i],1,
                        solplus[i],1);
        }
        NonlinSysEvaluator_coeff(solplus,resminus);

        for (int i = 0; i < nvariables; i++)
        {
            tmp = inarray + i*npoints;
            Vmath::Svtvp(npoints,2.0*eps,tmp,1,m_TimeIntegtSol_k[i],1,
                        solplus[i],1);
        }
        NonlinSysEvaluator_coeff(solplus,resPPls);

        for (int i = 0; i < nvariables; i++)
        {
            tmp = inarray + i*npoints;
            Vmath::Svtvp(npoints,-2.0*eps,tmp,1,m_TimeIntegtSol_k[i],1,
                        solplus[i],1);
        }
        NonlinSysEvaluator_coeff(solplus,resMMmns);

        Vmath::Zero(out.num_elements(),out,1);
        for (int i = 0; i < nvariables; i++)
        {
            tmp = out + i*npoints;
            Vmath::Svtvp(npoints, 1.0/12.0,&resMMmns[i][0],1,
                        &tmp[0],1,&tmp[0],1);
            Vmath::Svtvp(npoints,-8.0/12.0,&resminus[i][0],1,
                        &tmp[0],1,&tmp[0],1);
            Vmath::Svtvp(npoints, 8.0/12.0, &resplus[i][0],1,
                        &tmp[0],1,&tmp[0],1);
            Vmath::Svtvp(npoints,-1.0/12.0, &resPPls[i][0],1,
                        &tmp[0],1,&tmp[0],1);
            Vmath::Smul(npoints, oeps ,&tmp[0],1,&tmp[0],1);
        }
       
        return;
    }

    void CompressibleFlowSystem::MatrixMultiply_MatrixFree_coeff(
        const TensorOfArray1D<NekDouble>  &inarray, 
        TensorOfArray1D<NekDouble>        &outarray)
    {
        if (m_flagMultiLvlJacMultiplyMatFree)
        {   
            int nlvl = 1;
            m_EqdriverOperator.MultiLvlJacMultiplyMatFree(
                nlvl,
                inarray,
                outarray,
                m_BndEvaluateTime,
                m_TimeIntegLambda,
                m_TimeIntegtSol_k,
                m_updateMatFreeJacFlag);
        }
        else
        {
            v_MatrixMultiply_MatrixFree_coeff(
                inarray,
                outarray,
                m_BndEvaluateTime,
                m_TimeIntegLambda,
                m_TimeIntegtSol_k,
                m_updateMatFreeJacFlag);
        }
        m_updateMatFreeJacFlag = false;
    }


    void CompressibleFlowSystem::v_MatrixMultiply_MatrixFree_coeff(
        const TensorOfArray1D<NekDouble>  &inarray, 
        TensorOfArray1D<NekDouble>        &out, 
        const NekDouble                   time, 
        const NekDouble                   dtlamda, 
        const TensorOfArray2D<NekDouble>  &refFields, 
        const bool                        flagUpdateJac)
    {
        unsigned int nvariable  = m_fields.num_elements();
        unsigned int ncoeffs    = m_fields[0]->GetNcoeffs();
        unsigned int npoints    = m_fields[0]->GetNpoints();

        if (flagUpdateJac)
        {
            unsigned int ntotalGlobal     = nvariable * ncoeffs;
            m_comm->AllReduce(ntotalGlobal, Nektar::LibUtilities::ReduceSum);
            unsigned int ntotalDOF        = ntotalGlobal/nvariable;
            NekDouble ototalDOF     = 1.0/ntotalDOF;

            UpdateSoltnRefNorms(refFields, ototalDOF);

            CalcFluxJacVolBnd(refFields, time);
        }

        

        Array<OneD, Array<OneD, NekDouble> > inpnts(nvariable);
        Array<OneD, Array<OneD, NekDouble> > in2D(nvariable);
        Array<OneD, Array<OneD, NekDouble> > out2D(nvariable);
        
        for(int i = 0; i < nvariable; i++)
        {
            in2D[i]     = inarray + i*ncoeffs;
            out2D[i]    = out + i*ncoeffs;
            inpnts[i]   =   Array<OneD, NekDouble>(npoints,0.0);
            m_fields[i]->BwdTrans(in2D[i], inpnts[i]);
        }

        bool flag = true;
        DoOdeRhs_coeff(inpnts,out2D,time,flag);

        Vmath::Smul(nvariable*ncoeffs,dtlamda,out,1,out,1);
        Vmath::Vsub(nvariable*ncoeffs,inarray,1,out,1,out,1);
        return;
    }

    void CompressibleFlowSystem::MatrixMultiply_JacobianFree_coeff_dualtimestep(
        const  Array<OneD, NekDouble>   &inarray,
        Array<OneD, NekDouble >         &out,
        const  bool                     &controlFlag)
    {
        switch (m_LiniearizationMethod)
        {
        case 0:
            MatrixMultiply_JacobianFree_coeff(inarray,out);
            break;
        case 1:
            MatrixMultiply_JacobianFree_coeff_central (inarray,out);
            break;
        case 2:
        case 3:
            MatrixMultiply_MatrixFree_coeff (inarray,out);
            break;
        default:
            MatrixMultiply_JacobianFree_coeff(inarray,out);
            break;
        } 
// int nwidthcolm = 16;
// for(int i=0;i<out.num_elements();i++)
// {
//     cout <<std::scientific<<std::setw(nwidthcolm)<<std::setprecision(nwidthcolm-8) 
//          << " out["<<i<<"]= "<<out[i]<<endl;
// }
// ASSERTL0(false, "DebugStop");

        if(m_cflLocTimestep>0.0)
        {
            int nElements = m_fields[0]->GetExpSize();
            int nelmtcoeffs,offset,varoffset;
            NekDouble fac;
            Array<OneD, NekDouble> tmp;
            Array<OneD, NekDouble> tmpout;
            unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
            unsigned int ntotcoeffs = m_TimeIntegtSol_n[0].num_elements();

            Array<OneD, NekDouble> pseudotimefactor(nElements,0.0);
            NekDouble   otimestep   =   1.0/m_timestep;
            Vmath::Smul(nElements,otimestep,m_locTimeStep,1,pseudotimefactor,1);

            Vmath::Vsub(inarray.num_elements(),out,1,inarray,1,out,1);
            for(int i = 0; i < nvariables; ++i)
            {
                varoffset = i*ntotcoeffs;
                tmpout = out + varoffset;
                // Loop over elements
                for(int n = 0; n < nElements; ++n)
                {
                    nelmtcoeffs = m_fields[0]->GetExp(n)->GetNcoeffs();
                    offset      = m_fields[0]->GetCoeff_Offset(n);
                    fac         = pseudotimefactor[n];
                    tmp         = tmpout + offset;
                    Vmath::Smul(nelmtcoeffs,fac,tmp,1,tmp,1);
                }
            }
            Vmath::Vadd(inarray.num_elements(),out,1,inarray,1,out,1);
        }
        return;
    }

    void CompressibleFlowSystem::NonlinSysEvaluator_coeff(
        Array<OneD, Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &out)
    {
        Array<OneD, Array<OneD, NekDouble> > sol_n;
        sol_n                  = m_TimeIntegtSol_n;
        //inforc = m_TimeIntegForce;
        unsigned int nvariable  = m_fields.num_elements();
        unsigned int ncoeffs    = GetNcoeffs();

        NonlinSysEvaluator_coeff_noSource(inarray,out);

        for (int i = 0; i < nvariable; ++i)
        {
            Vmath::Vsub(ncoeffs,out[i],1,sol_n[i],1,out[i],1);
        }
        return;
    }

     void CompressibleFlowSystem::NonlinSysEvaluator_coeff_noSource(
                                                       Array<OneD, Array<OneD, NekDouble> > &inarray,
                                                       Array<OneD, Array<OneD, NekDouble> > &out)
    {
        //inforc = m_TimeIntegForce;
        unsigned int nvariable  = m_fields.num_elements();
        unsigned int ncoeffs    = GetNcoeffs();
        unsigned int npoints    = m_fields[0]->GetNpoints();

        Array<OneD, Array<OneD, NekDouble> > inpnts(nvariable);

        for(int i = 0; i < nvariable; i++)
        {
            inpnts[i]   =   Array<OneD, NekDouble>(npoints,0.0);
            m_fields[i]->BwdTrans(inarray[i], inpnts[i]);
        }

        DoOdeProjection(inpnts,inpnts,m_BndEvaluateTime);
        DoOdeRhs_coeff(inpnts,out,m_BndEvaluateTime);

        for (int i = 0; i < nvariable; i++)
        {
            Vmath::Smul(ncoeffs,m_TimeIntegLambda,out[i],1,out[i],1);
            Vmath::Vsub(ncoeffs,inarray[i],1,out[i],1,out[i],1);
        }
        return;
    }

    /**
     * @brief Compute the right-hand side.
     */
    void CompressibleFlowSystem::DoOdeRhs_coeff(
        const TensorOfArray2D<NekDouble>        &inarray,
        Array<OneD, Array<OneD, NekDouble> >    &outarray,
        const NekDouble                         time,
        const bool                              flagFreezeJac)
    {
        int nvariables = inarray.num_elements();
        int nTracePts  = GetTraceTotPoints();
        int ncoeffs    = GetNcoeffs();
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
            if(flagFreezeJac)
            {
                for(int i = 0; i < nvariables; ++i)
                {
                    Fwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                    Bwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                    m_fields[i]->GetFwdBwdTracePhysNoBndFill(inarray[i], 
                        Fwd[i], Bwd[i]);
                }
            }
            else
            {
                for(int i = 0; i < nvariables; ++i)
                {
                    Fwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                    Bwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                    m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
                }
            }
        }
 
         // Calculate advection
        DoAdvection_coeff(inarray, outarray, time, Fwd, Bwd,flagFreezeJac);
        // Negate results
        
        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(ncoeffs, outarray[i], 1);
        }
#ifdef CFS_DEBUGMODE
        if(2==m_DebugAdvDiffSwitch)
        {
            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Zero(ncoeffs, outarray[i], 1);
            }
        }
#endif

#ifdef CFS_DEBUGMODE
        if(1!=m_DebugAdvDiffSwitch)
        {
#endif
        DoDiffusion_coeff(inarray, outarray, Fwd, Bwd,flagFreezeJac);
#ifdef CFS_DEBUGMODE
        }
#endif
        // Add forcing terms
        for (auto &x : m_forcing)
        {
            if(true==flagFreezeJac)
            {
                ASSERTL0(false,"forcing not coded for DoOdeRhs_coeffMF");
            }
            // x->Apply(m_fields, inarray, outarray, time);
            x->Apply_coeff(m_fields, inarray, outarray, time);
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
            for(int n = 0; n < nElements; ++n)
            {
                nq     = m_fields[0]->GetExp(n)->GetTotPoints();
                offset = m_fields[0]->GetPhys_Offset(n);
                fac    = tstep[n] / m_timestep;
                for(int i = 0; i < nvariables; ++i)
                {
                    Vmath::Smul(nq, fac, outarray[i] + offset, 1,
                                         tmp = outarray[i] + offset, 1);
                }
            }
        }
    }

    void CompressibleFlowSystem::DoOdeRhs_coeffVol(
        const TensorOfArray2D<NekDouble>        &inarray,
        TensorOfArray2D<NekDouble>              &outarray,
        const NekDouble                         time,
        const bool                              flagFreezeJac)
    {
        int nvariables = inarray.num_elements();
        int nTracePts  = GetTraceTotPoints();
        int ncoeffs    = GetNcoeffs();
        
         // Calculate advection
        DoAdvection_coeffVol(inarray, outarray, time, flagFreezeJac);
        // Negate results
        
        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(ncoeffs, outarray[i], 1);
        }
#ifdef CFS_DEBUGMODE
        if(2==m_DebugAdvDiffSwitch)
        {
            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Zero(ncoeffs, outarray[i], 1);
            }
        }
#endif

#ifdef CFS_DEBUGMODE
        if(1!=m_DebugAdvDiffSwitch)
        {
#endif
        DoDiffusion_coeffVol(inarray, outarray, flagFreezeJac);
#ifdef CFS_DEBUGMODE
        }
#endif
        // Add forcing terms
        for (auto &x : m_forcing)
        {
            if(true==flagFreezeJac)
            {
                ASSERTL0(false,"forcing not coded for DoOdeRhs_coeffMF");
            }
            // x->Apply(m_fields, inarray, outarray, time);
            x->Apply_coeff(m_fields, inarray, outarray, time);
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
            for(int n = 0; n < nElements; ++n)
            {
                nq     = m_fields[0]->GetExp(n)->GetTotPoints();
                offset = m_fields[0]->GetPhys_Offset(n);
                fac    = tstep[n] / m_timestep;
                for(int i = 0; i < nvariables; ++i)
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
    void CompressibleFlowSystem::DoAdvection_coeff(
        const TensorOfArray2D<NekDouble>            &inarray,
        Array<OneD, Array<OneD, NekDouble> >        &outarray,
        const NekDouble                             time,
        const Array<OneD, Array<OneD, NekDouble> >  &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >  &pBwd,
        const bool                                  flagFreezeJac)
    {
        int nvariables = inarray.num_elements();
        Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);

        m_advObject->AdvectCoeff(nvariables, m_fields, advVel, inarray,
                            outarray, time, pFwd, pBwd,flagFreezeJac);
    }

    void CompressibleFlowSystem::DoAdvection_coeffVol(
        const TensorOfArray2D<NekDouble>        &inarray,
        Array<OneD, Array<OneD, NekDouble> >    &outarray,
        const NekDouble                         time,
        const bool                              flagFreezeJac)
    {
        int nvariables = inarray.num_elements();
        Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);

        m_advObject->AdvectCoeffVol(nvariables, m_fields, advVel, inarray,
                            outarray, time, flagFreezeJac);
    }

#ifdef CFS_DEBUGMODE

    void CompressibleFlowSystem::DebugNumCalJac_coeff(
        TensorOfArray2D<DNekBlkMatSharedPtr>    &gmtxarray,
        Array<OneD, Array<OneD, NekDouble > >   &JacOffDiagArray)
    {
        MultiRegions::ExpListSharedPtr explist = m_fields[0];
        std::shared_ptr<LocalRegions::ExpansionVector> pexp = 
            explist->GetExp();
        int nElmts    = (*pexp).size();
        int nvariables= gmtxarray.num_elements();

        TensorOfArray2D<DNekMatSharedPtr>  ElmtPrecMatVars(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            ElmtPrecMatVars[i]  =   Array<OneD, DNekMatSharedPtr>(nvariables);
        }
        if(JacOffDiagArray.num_elements()>0)
        {
            for(int i = 0; i < nElmts; i++)
            {
                DebugNumCalElmtJac_coeff(ElmtPrecMatVars,i,JacOffDiagArray);
                for(int m = 0; m < nvariables; m++)
                {
                    for(int n = 0; n < nvariables; n++)
                    {
                        gmtxarray[m][n]->SetBlock(i,i,ElmtPrecMatVars[m][n]);
                    }
                }
            }
        }
        else
        {
            for(int i = 0; i < nElmts; i++)
            {
                DebugNumCalElmtJac_coeff(ElmtPrecMatVars,i);
                for(int m = 0; m < nvariables; m++)
                {
                    for(int n = 0; n < nvariables; n++)
                    {
                        gmtxarray[m][n]->SetBlock(i,i,ElmtPrecMatVars[m][n]);
                    }
                }
            }
        }
        
        
    }

    void CompressibleFlowSystem::DebugNumCalElmtJac_coeff(
        TensorOfArray2D<DNekMatSharedPtr>       &ElmtPrecMatVars ,
        const int                               nelmt,
        Array<OneD, Array<OneD, NekDouble > >   &JacOffDiagArray)
    {
        MultiRegions::ExpListSharedPtr explist = m_fields[0];
        std::shared_ptr<LocalRegions::ExpansionVector> pexp = 
            explist->GetExp();
        int nElmtcoef   = (*pexp)[nelmt]->GetNcoeffs();
        int nElmtOffset = explist->GetCoeff_Offset(nelmt);
        int nTotCoef    = GetNcoeffs();
        
        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int npoints    = m_TimeIntegtSol_n[0].num_elements();
        unsigned int ntotpnt    = nvariables*npoints;
        Array<OneD, NekDouble > tmpinn_1d(ntotpnt,0.0);
        Array<OneD, NekDouble > tmpout_1d(ntotpnt,0.0);
        Array<OneD,       Array<OneD, NekDouble> > tmpinn(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > tmpout(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            int noffset = i*npoints;
            tmpinn[i] = tmpinn_1d+noffset;
            tmpout[i] = tmpout_1d+noffset;
        }


        for(int i = 0; i < nvariables; i++)
        {
            for(int j = 0; j < nvariables; j++)
            {
                ElmtPrecMatVars[i][j] =  MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nElmtcoef, nElmtcoef, 0.0);
            }
        }

        DNekMatSharedPtr    tmpStdMat;
        for (int i = 0; i < nvariables; i++)
        {
            for (int npnt = 0; npnt < nElmtcoef; npnt++)
            {
                tmpinn[i][nElmtOffset+npnt] = 1.0;
                // MatrixMultiply_MatrixFree_coeff_central(tmpinn_1d,tmpout_1d);
                MatrixMultiply_JacobianFree_coeff_FourthCentral(tmpinn_1d,
                                                              tmpout_1d);

                Vmath::Vcopy(ntotpnt,&tmpout_1d[0],1,
                    &JacOffDiagArray[i*nTotCoef+nElmtOffset+npnt][0],1);
                
                for (int j = 0; j < nvariables; j++)
                {
                    for (int npntf = 0; npntf < nElmtcoef; npntf++)
                    {
                        tmpStdMat = ElmtPrecMatVars[j][i];
                        // NekDouble tmp = tmpout[j][nElmtOffset+npntf];
                        (*tmpStdMat)(npntf,npnt) = tmpout[j][nElmtOffset+npntf];

                        JacOffDiagArray[i*nTotCoef+nElmtOffset+npnt]
                            [j*nTotCoef+nElmtOffset+npntf] = 0.0;
                    }
                }

                tmpinn[i][nElmtOffset+npnt] = 0.0;
            }
        }
    }

    void CompressibleFlowSystem::DebugNumCalElmtJac_coeff(
        TensorOfArray2D<DNekMatSharedPtr> &ElmtPrecMatVars ,
        const int                         nelmt)
    {
        MultiRegions::ExpListSharedPtr explist = m_fields[0];
        std::shared_ptr<LocalRegions::ExpansionVector> pexp = 
            explist->GetExp();
        int nElmtcoef   = (*pexp)[nelmt]->GetNcoeffs();
        int nElmtOffset = explist->GetCoeff_Offset(nelmt);
        int nTotCoef    = GetNcoeffs();
        
        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int npoints    = m_TimeIntegtSol_n[0].num_elements();
        unsigned int ntotpnt    = nvariables*npoints;
        Array<OneD, NekDouble > tmpinn_1d(ntotpnt,0.0);
        Array<OneD, NekDouble > tmpout_1d(ntotpnt,0.0);
        Array<OneD,       Array<OneD, NekDouble> > tmpinn(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > tmpout(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            int noffset = i*npoints;
            tmpinn[i] = tmpinn_1d+noffset;
            tmpout[i] = tmpout_1d+noffset;
        }


        for(int i = 0; i < nvariables; i++)
        {
            for(int j = 0; j < nvariables; j++)
            {
                ElmtPrecMatVars[i][j] =  MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nElmtcoef, nElmtcoef, 0.0);
            }
        }

        DNekMatSharedPtr    tmpStdMat;
        for (int i = 0; i < nvariables; i++)
        {
            for (int npnt = 0; npnt < nElmtcoef; npnt++)
            {
                tmpinn[i][nElmtOffset+npnt] = 1.0;
                MatrixMultiply_JacobianFree_coeff_central(tmpinn_1d,tmpout_1d);

                for (int j = 0; j < nvariables; j++)
                {
                    for (int npntf = 0; npntf < nElmtcoef; npntf++)
                    {
                        tmpStdMat = ElmtPrecMatVars[j][i];
                        // NekDouble tmp = tmpout[j][nElmtOffset+npntf];
                        (*tmpStdMat)(npntf,npnt) = tmpout[j][nElmtOffset+npntf];
                    }
                }

                tmpinn[i][nElmtOffset+npnt] = 0.0;
            }
        }
    }
    
    void CompressibleFlowSystem::DebugCheckJac(
        const int                 RowElementID,
        const int                 ColElementID)
    {
        int nelmts  = m_fields[0]->GetNumElmts();
        int nrowsVars   = m_fields[0]->GetExp(RowElementID)->GetNcoeffs();
        int ncolsVars   = m_fields[0]->GetExp(ColElementID)->GetNcoeffs();
        Array<OneD, unsigned int > nelmtmatdim(nelmts);
        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();

        DNekMatSharedPtr  NumericalJacobianMatrix = MemoryManager<DNekMat>::
            AllocateSharedPtr(nrowsVars*nvariables,ncolsVars*nvariables,0.0);
        
        NumJacElemental(NumericalJacobianMatrix, RowElementID, ColElementID);

        cout <<" NumericalJacobianMatrix "<<endl;
        CoutStandardMat(NumericalJacobianMatrix);

        if(RowElementID==ColElementID)
        {
            int intexrw= 0;
            for(int i = 0; i<nvariables;i++)
            {
                intexrw= i;
                cout <<" m_PrecMatVars["<<intexrw<<"][0] "<<endl;
                CoutStandardMat(m_PrecMatVars[intexrw][0]->GetBlock(0,0));
                cout <<" m_PrecMatVars["<<intexrw<<"][1] "<<endl;
                CoutStandardMat(m_PrecMatVars[intexrw][1]->GetBlock(0,0));
                cout <<" m_PrecMatVars["<<intexrw<<"][2] "<<endl;
                CoutStandardMat(m_PrecMatVars[intexrw][2]->GetBlock(0,0));
                cout <<" m_PrecMatVars["<<intexrw<<"][3] "<<endl;
                CoutStandardMat(m_PrecMatVars[intexrw][3]->GetBlock(0,0));
            }
        }
        else
        {
            DNekMatSharedPtr  MinusoffJacobianMatrix = MemoryManager<DNekMat>::
                AllocateSharedPtr(nrowsVars*nvariables,ncolsVars*nvariables,0.0);
            CalOffDiagJacByMinusOffDiagElemental(MinusoffJacobianMatrix,
                RowElementID,ColElementID);

            cout <<" MinusoffJacobianMatrix "<<endl;
            CoutStandardMat(MinusoffJacobianMatrix);
        }
    }

    void CompressibleFlowSystem::CalOffDiagJacByMinusOffDiagElemental(
        DNekMatSharedPtr    &MinusoffJacobianMatrix,
        const int           RowElementID,
        const int           ColElementID)
    {
        unsigned int nvariables     = m_fields.num_elements();
        unsigned int nphspnt        = m_fields[0]->GetNpoints();
        unsigned int   ncoeffs      = m_fields[0]->GetNcoeffs();

        int nelmts  = m_fields[0]->GetNumElmts();
        int nrowsVars   = m_fields[0]->GetExp(RowElementID)->GetNcoeffs();
        int ncolsVars   = m_fields[0]->GetExp(ColElementID)->GetNcoeffs();

        int nTracePts  = GetTraceTotPoints();
        TensorOfArray3D<NekDouble> qfield(m_spacedim);
        for(int i = 0; i< m_spacedim; i++)
        {
            qfield[i]   =   Array<OneD, Array<OneD, NekDouble> >(nvariables);
            for(int j = 0; j< nvariables; j++)
            {
                qfield[i][j]   =   Array<OneD, NekDouble>(nphspnt,0.0);
            }
        }
        int ntmpTrace = 4+2*m_spacedim;
        TensorOfArray3D<NekDouble> tmpTrace(ntmpTrace);
        for(int i = 0; i< ntmpTrace; i++)
        {
            tmpTrace[i]   =   Array<OneD, Array<OneD, NekDouble> >(nvariables);
            for(int j = 0; j< nvariables; j++)
            {
                tmpTrace[i][j]   =   Array<OneD, NekDouble>(nTracePts,0.0);
            }
        }
        Array<OneD, Array<OneD, NekDouble> > FwdFluxDeriv(nvariables);
        Array<OneD, Array<OneD, NekDouble> > BwdFluxDeriv(nvariables);
        for(int j = 0; j< nvariables; j++)
        {
            FwdFluxDeriv[j]   =   Array<OneD, NekDouble>(nTracePts,0.0);
            BwdFluxDeriv[j]   =   Array<OneD, NekDouble>(nTracePts,0.0);
        }

        Array<OneD, Array<OneD, NekDouble> > coeftmp(nvariables);
        Array<OneD, Array<OneD, NekDouble> > rhstmp(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            coeftmp[i]    =   Array<OneD, NekDouble>(ncoeffs,0.0);
            rhstmp[i]    =   Array<OneD, NekDouble>(ncoeffs,0.0);
        }

        bool flagUpdateDervFlux = true;

        int nColElmtOffset  = m_fields[0]->GetCoeff_Offset(ColElementID);
        int nRowElmtOffset  = m_fields[0]->GetCoeff_Offset(RowElementID);
        
        for(int nvcl = 0; nvcl < nvariables; nvcl++)
        {
            int Cloffset = nvcl * ncolsVars;
            for(int nccl = 0; nccl < ncolsVars; nccl++)
            {
                coeftmp[nvcl][nColElmtOffset + nccl] = 1.0;

                MinusOffDiag2Rhs(nvariables,ncoeffs,rhstmp,coeftmp,
                    flagUpdateDervFlux,FwdFluxDeriv,BwdFluxDeriv,qfield,
                    tmpTrace, m_TraceJac, m_TraceJacDeriv, m_TraceJacDerivSign);

                for(int nvrw = 0; nvrw < nvariables; nvrw++)
                {
                    int Rwoffset = nvrw *  nrowsVars;
                    for(int ncrw = 0; ncrw < nrowsVars; ncrw++)
                    {
                        (*MinusoffJacobianMatrix)(Rwoffset + ncrw, 
                            Cloffset + nccl) =
                            coeftmp[nvrw][nRowElmtOffset + ncrw];
                    }

                }
                for(int i = 0; i < nvariables; i++)
                {
                    Vmath::Zero(ncoeffs,coeftmp[i],1);
                }
            }
        }
    }

    void CompressibleFlowSystem::NumJacElemental(
        DNekMatSharedPtr    &NumericalJacobianMatrix,
        const int           RowElementID,
        const int           ColElementID)
    {
        MultiRegions::ExpListSharedPtr explist              = m_fields[0];
        std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
        int nElmts                                          = (*pexp).size();
        int nvariables      = m_fields.num_elements();
        int ntotalcoef      = m_TimeIntegtSol_n[0].num_elements();
        int ntotaltotalcoef = nvariables * ntotalcoef;

        Array<OneD, NekDouble> tmpinn_1d(ntotaltotalcoef, 0.0);
        Array<OneD, NekDouble> tmpout_1d(ntotaltotalcoef, 0.0);
        Array<OneD, Array<OneD, NekDouble>> tmpinn(nvariables);
        Array<OneD, Array<OneD, NekDouble>> tmpout(nvariables);
        for (int i = 0; i < nvariables; i++)
        {
            int noffset = i * ntotalcoef;
            tmpinn[i]   = tmpinn_1d + noffset;
            tmpout[i]   = tmpout_1d + noffset;
        }

        NekDouble magnitdEstimatMax =0.0;
        for(int i = 0; i < m_magnitdEstimat.num_elements(); i++)
        {
            magnitdEstimatMax = max(magnitdEstimatMax,m_magnitdEstimat[i]);
        }
        magnitdEstimatMax = 1.0/magnitdEstimatMax;

        Array<OneD, NekDouble> 
            ScaledmagnitdEstimat(m_magnitdEstimat.num_elements());
        Array<OneD, NekDouble> 
            oScaledmagnitdEstimat(m_magnitdEstimat.num_elements());
        for(int i = 0; i < m_magnitdEstimat.num_elements(); i++)
        {
            ScaledmagnitdEstimat[i] = m_magnitdEstimat[i]*magnitdEstimatMax;
            oScaledmagnitdEstimat[i] = 1.0/ScaledmagnitdEstimat[i];
        }

        int nColumnElmtcoef   = (*pexp)[ColElementID]->GetNcoeffs();
        int nColumnElmtOffset = explist->GetCoeff_Offset(ColElementID);
        int nRowElmtcoef      = (*pexp)[RowElementID]->GetNcoeffs();
        int nRowElmtOffset    = explist->GetCoeff_Offset(RowElementID);

        for (int i = 0; i < nvariables; i++)
        {
            int ioffset = i * nColumnElmtcoef;
            for (int npts = 0; npts < nColumnElmtcoef; npts++)
            {
                tmpinn[i][nColumnElmtOffset + npts] = ScaledmagnitdEstimat[i];
                
                MatrixMultiply_JacobianFree_coeff_central(tmpinn_1d, tmpout_1d);
                
                for (int j = 0; j < nvariables; j++)
                {
                    int joffset = j *  nRowElmtcoef;
                    for (int mpts = 0; mpts < nRowElmtcoef; mpts++)
                    {
                        (*NumericalJacobianMatrix)(joffset+ mpts,ioffset  +npts) 
                            = tmpout[j][nRowElmtOffset + mpts]
                                *oScaledmagnitdEstimat[i];
                    }

                }
                tmpinn[i][nColumnElmtOffset + npts] = 0.0;
            }
        }
    }

    void CompressibleFlowSystem::NonlinSysEvaluator_coeff_out(
        Array<OneD, Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &out)
    {
        Array<OneD, Array<OneD, NekDouble> > sol_n;
        sol_n                  = m_TimeIntegtSol_n;
        //inforc = m_TimeIntegForce;
        unsigned int nvariable  = inarray.num_elements();
        unsigned int ncoeffs    = inarray[0].num_elements();
        unsigned int npoints    = m_fields[0]->GetNpoints();

        Array<OneD, Array<OneD, NekDouble> > inpnts(nvariable);

        for(int i = 0; i < nvariable; i++)
        {
            inpnts[i]   =   Array<OneD, NekDouble>(npoints,0.0);
            m_fields[i]->BwdTrans(inarray[i], inpnts[i]);
        }

        DoOdeProjection(inpnts,inpnts,m_BndEvaluateTime);
        DoOdeRhs_coeff(inpnts,out,m_BndEvaluateTime);


        LibUtilities::CommSharedPtr v_Comm = 
                m_fields[0]->GetComm()->GetRowComm();
        NekDouble resnorm = 0.0;
        for(int i = 0; i < nvariable; i++)
        {
            resnorm += Vmath::Dot(ncoeffs,out[i],out[i]);
        }
        v_Comm->AllReduce(resnorm, Nektar::LibUtilities::ReduceSum);
        if(0==m_session->GetComm()->GetRank())
        {
            cout << "        m_TimeIntegLambda= "<<m_TimeIntegLambda
                 <<" DtRhs^2= "<<m_TimeIntegLambda*m_TimeIntegLambda*resnorm
                 <<endl;
        }

        for (int i = 0; i < nvariable; i++)
        {
            Vmath::Svtvp(ncoeffs,m_TimeIntegLambda,out[i],1,
                         sol_n[i],1,out[i],1);
            Vmath::Vsub(ncoeffs,inarray[i],1,out[i],1,out[i],1);
        }
        return;
    }

    template<typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::Cout1DArrayBlkMat(
        Array<OneD, TypeNekBlkMatSharedPtr> &gmtxarray,
        const unsigned int                  nwidthcolm)
    {
        int nvar1 = gmtxarray.num_elements();


        for(int i = 0; i < nvar1; i++)
        {
            cout << endl 
            << "£$£$£$£$£$£$££$£$£$$$£$££$$£$££$£$$££££$$£$£$£$£$£$£$££$£$$"
            << endl << "Cout1DArrayBlkMat i= " << i << endl;
            CoutBlkMat(gmtxarray[i],nwidthcolm);
        }
    }
    
    template<typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::Cout2DArrayBlkMat(
        TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
        const unsigned int                      nwidthcolm)
    {
        int nvar1 = gmtxarray.num_elements();
        int nvar2 = gmtxarray[0].num_elements();


        for(int i = 0; i < nvar1; i++)
        {
            for(int j = 0; j < nvar2; j++)
            {
                cout << endl 
                << "£$£$£$£$£$£$££$£$£$$$£$££$$£$££$£$$££££$$£$£$£$£$£$£$££$£$$"
                << endl << "Cout2DArrayBlkMat i= " 
                << i << " j=  " << j << endl;
                CoutBlkMat(gmtxarray[i][j],nwidthcolm);
            }
        }
    }

    template<typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::CoutBlkMat(TypeNekBlkMatSharedPtr &gmtx,const unsigned int nwidthcolm)
    {

        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;
        gmtx->GetBlockSizes(rowSizes,colSizes);

        int nelmts  = rowSizes.num_elements();
        
        // int noffset = 0;
        for(int i = 0; i < nelmts; ++i)
        {
            std::cout <<std::endl <<"*********************************"
                      <<std::endl<<"element :   "<<i<<std::endl;
            CoutStandardMat(gmtx->GetBlock(i,i),nwidthcolm);
        }
        return;
    }
   
    template<typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::Cout2DArrayStdMat(
        TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
        const unsigned int                      nwidthcolm)
    {
        int nvar1 = gmtxarray.num_elements();
        int nvar2 = gmtxarray[0].num_elements();


        for(int i = 0; i < nvar1; i++)
        {
            for(int j = 0; j < nvar2; j++)
            {
                cout << endl 
                << "£$£$£$£$£$£$££$£$£$$$£$££$$£$££$£$$££££$$£$£$£$£$£$£$££$£$$"
                << endl << "Cout2DArrayBlkMat i= " 
                << i << " j=  " << j << endl;
                CoutStandardMat(gmtxarray[i][j],nwidthcolm);
            }
        }
    }

    template<typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::CoutStandardMat(
        TypeNekBlkMatSharedPtr  &loc_matNvar,
        const unsigned int      nwidthcolm)
    {
        int nrows = loc_matNvar->GetRows();
        int ncols = loc_matNvar->GetColumns();
        std::cout   <<"ROW="<<std::setw(3)<<-1<<" ";
        for(int k = 0; k < ncols; k++)
        {
            std::cout   <<"   COL="<<std::setw(nwidthcolm-7)<<k;
        }
        std::cout   << endl;

        for(int j = 0; j < nrows; j++)
        {
            std::cout   <<"ROW="<<std::setw(3)<<j<<" ";
            for(int k = 0; k < ncols; k++)
            {
                std::cout <<std::scientific<<std::setw(nwidthcolm)
                          <<std::setprecision(nwidthcolm-8)
                          <<(*loc_matNvar)(j,k);
            }
            std::cout   << endl;
        }
    }
#endif
    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::TranSamesizeBlkDiagMatIntoArray(
        const TypeNekBlkMatSharedPtr &BlkMat,
        TensorOfArray3D<DataType>    &MatArray)
    {
        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;
        BlkMat->GetBlockSizes(rowSizes,colSizes);
        int nDiagBlks   = rowSizes.num_elements();
        int nvar0       = rowSizes[1] - rowSizes[0];
        int nvar1       = colSizes[1] - colSizes[0];

        Array<OneD, DataType>    ElmtMatData;

        for(int i=0;i<nDiagBlks;i++)
        {
            ElmtMatData = BlkMat->GetBlock(i,i)->GetPtr();
            for(int n=0;n<nvar1;n++)
            {
                int noffset = n*nvar0;
                for(int m=0;m<nvar0;m++)
                {
                    MatArray[m][n][i]   =   ElmtMatData[m+noffset];
                }
            }
        }
    }
    
    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::Fill2DArrayOfBlkDiagonalMat(
        TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
        const DataType                          valu)
    {

        int n1d = gmtxarray.num_elements();

        for(int n1 = 0; n1 < n1d; ++n1)
        {
            Fill1DArrayOfBlkDiagonalMat(gmtxarray[n1],valu);
        }
    }

    
    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::Fill1DArrayOfBlkDiagonalMat( 
        Array<OneD, TypeNekBlkMatSharedPtr >    &gmtxarray,
        const DataType                          valu)
    {
        int n1d = gmtxarray.num_elements();

        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;

        Array<OneD, DataType > loc_mat_arr;


        for(int n1 = 0; n1 < n1d; ++n1)
        {
            gmtxarray[n1]->GetBlockSizes(rowSizes,colSizes);
            int nelmts  = rowSizes.num_elements();

            for(int i = 0; i < nelmts; ++i)
            {
                loc_mat_arr = gmtxarray[n1]->GetBlock(i,i)->GetPtr();

                int nrows = gmtxarray[n1]->GetBlock(i,i)->GetRows();
                int ncols = gmtxarray[n1]->GetBlock(i,i)->GetColumns();

                Vmath::Fill(nrows*ncols,valu,loc_mat_arr,1);
            }
        }

    }

    void CompressibleFlowSystem::GetFluxVectorJacDirctn(
        const int                        nDirctn,
        const TensorOfArray2D<NekDouble> &inarray,
        TensorOfArray5D<NekDouble>       &ElmtJacArray)
    {
        int nConvectiveFields   = inarray.num_elements();
        std::shared_ptr<LocalRegions::ExpansionVector> expvect = 
                m_fields[0]->GetExp();
        int ntotElmt            = (*expvect).size();
        LocalRegions::ExpansionSharedPtr pexp = (*expvect)[0];
        // int nElmtPnt            = pexp->GetTotPoints();
        // int nElmtCoef           = pexp->GetNcoeffs();
        Array<OneD, NekDouble> normals;
        Array<OneD, Array<OneD, NekDouble> > normal3D(3);
        for(int i = 0; i < 3; i++)
        {
            normal3D[i] = Array<OneD, NekDouble>(3,0.0);
        }
        normal3D[0][0] = 1.0;
        normal3D[1][1] = 1.0;
        normal3D[2][2] = 1.0;
        normals =   normal3D[nDirctn];

        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);
        Array<OneD, Array<OneD, NekDouble> > locVars(nConvectiveFields);
        DNekMatSharedPtr pointJac = MemoryManager<DNekMat>::
                AllocateSharedPtr(nConvectiveFields,nConvectiveFields,0.0);

        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            // nElmtCoef           = (*expvect)[nelmt]->GetNcoeffs();
            int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();

            for(int j = 0; j < nConvectiveFields; j++)
            {
                locVars[j] = inarray[j]+GetPhys_Offset(nelmt);
            }

            for(int npnt = 0; npnt < nElmtPnt; npnt++)
            {
                for(int j = 0; j < nConvectiveFields; j++)
                {
                    pointVar[j] = locVars[j][npnt];
                }

                GetFluxVectorJacPoint(nConvectiveFields,pointVar,
                                        normals,pointJac);

                for(int m=0;m<nConvectiveFields;m++)
                {
                    for(int n=0;n<nConvectiveFields;n++)
                    {
                        ElmtJacArray[m][n][nDirctn][nelmt][npnt] = 
                                (*pointJac)(m,n);
                    }
                }
            }
        }
        return ;
    }

    void CompressibleFlowSystem::GetFluxVectorJacDirctnElmt(
        const int                                  nConvectiveFields,
        const int                                  nElmtPnt,
        const Array<OneD, Array<OneD, NekDouble> > &locVars,
        const Array<OneD, NekDouble>               &normals,
        DNekMatSharedPtr                           &wspMat,
        Array<OneD, Array<OneD, NekDouble> >       &PntJacArray)
    {
        Array<OneD, NekDouble> wspMatData = wspMat->GetPtr();

        int matsize = nConvectiveFields*nConvectiveFields;

        Array<OneD, NekDouble> pointVar(nConvectiveFields);

        for(int npnt = 0; npnt < nElmtPnt; npnt++)
        {
            for(int j = 0; j < nConvectiveFields; j++)
            {
                pointVar[j] = locVars[j][npnt];
            }

            GetFluxVectorJacPoint(nConvectiveFields,pointVar,normals,wspMat);

            Vmath::Vcopy(matsize, &wspMatData[0],1,&PntJacArray[npnt][0],1);
        }
        return ;
    }

    void CompressibleFlowSystem::GetFluxVectorJacDirctnMat(
        const int                            nDirctn,
        const TensorOfArray2D<NekDouble>     &inarray,
        TensorOfArray2D<DNekBlkMatSharedPtr> &ElmtFluxJacArray)
    {
        int nConvectiveFields   = inarray.num_elements();
        std::shared_ptr<LocalRegions::ExpansionVector> expvect = 
                m_fields[0]->GetExp();
        int ntotElmt            = (*expvect).size();
        LocalRegions::ExpansionSharedPtr pexp = (*expvect)[0];
        // int nElmtPnt            = pexp->GetTotPoints();
        // int nElmtCoef           = pexp->GetNcoeffs();
        Array<OneD, NekDouble> normals;
        Array<OneD, Array<OneD, NekDouble> > normal3D(3);
        for(int i = 0; i < 3; i++)
        {
            normal3D[i] = Array<OneD, NekDouble>(3,0.0);
        }
        normal3D[0][0] = 1.0;
        normal3D[1][1] = 1.0;
        normal3D[2][2] = 1.0;
        normals =   normal3D[nDirctn];

        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);
        Array<OneD, Array<OneD, NekDouble> > locVars(nConvectiveFields);
        DNekMatSharedPtr pointJac = MemoryManager<DNekMat>::
                    AllocateSharedPtr(nConvectiveFields,nConvectiveFields,0.0);
        Array<OneD, NekDouble> pointJacData = pointJac->GetPtr();
        Array<OneD, NekDouble> ElmtFluxJacData;

        int nMatMemoSize = nConvectiveFields*nConvectiveFields;
        NekDouble fsw,efix_StegerWarming;
        efix_StegerWarming = 0.0;
        fsw = 0.0;

        switch (m_spacedim)
        {
        case 3:
            for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
            {
                int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
        
                for(int j = 0; j < nConvectiveFields; j++)
                {   
                    locVars[j] = inarray[j]+GetPhys_Offset(nelmt);
                }

                for(int npnt = 0; npnt < nElmtPnt; npnt++)
                {
                    for(int j = 0; j < nConvectiveFields; j++)
                    {
                        pointVar[j] = locVars[j][npnt];
                    }

                    PointFluxJacobian_pn(pointVar,normals,pointJac,
                                        efix_StegerWarming,fsw);

                    ElmtFluxJacData = ElmtFluxJacArray[nelmt][nDirctn]->
                                        GetBlock(npnt,npnt)->GetPtr();
                    Vmath::Vcopy(nMatMemoSize,pointJacData,1,ElmtFluxJacData,1);
                }
            }
            break;
        case 2:
            for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
            {
                int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
        
                for(int j = 0; j < nConvectiveFields; j++)
                {   
                    locVars[j] = inarray[j]+GetPhys_Offset(nelmt);
                }

                for(int npnt = 0; npnt < nElmtPnt; npnt++)
                {
                    for(int j = 0; j < nConvectiveFields; j++)
                    {
                        pointVar[j] = locVars[j][npnt];
                    }

                    PointFluxJacobian_pn2D(pointVar,normals,pointJac,
                                        efix_StegerWarming,fsw);

                    ElmtFluxJacData = ElmtFluxJacArray[nelmt][nDirctn]->
                                        GetBlock(npnt,npnt)->GetPtr();
                    Vmath::Vcopy(nMatMemoSize,pointJacData,1,ElmtFluxJacData,1);
                }
            }
            break;
        default:
            ASSERTL0(false, "GetFluxVectorJacDirctnMat not coded");
            break;
        }
        return ;
    }
    
    void CompressibleFlowSystem::GetFluxVectorJacPoint(
        const int                                   nConvectiveFields,
        const Array<OneD, NekDouble>                &conservVar, 
        const Array<OneD, NekDouble>                &normals, 
        DNekMatSharedPtr                            &fluxJac)
    {
        int nvariables      = conservVar.num_elements();
        const int nvariables3D    = 5;
        int expDim          = m_spacedim;

        NekDouble fsw,efix_StegerWarming;
        efix_StegerWarming = 0.0;
        fsw = 0.0; // exact flux Jacobian if fsw=0.0
        if (nvariables > expDim+2)
        {
            ASSERTL0(false,"nvariables > expDim+2 case not coded")
        }

        Array<OneD, NekDouble> fluxJacData;
        ;
        fluxJacData = fluxJac->GetPtr();

        if(nConvectiveFields==nvariables3D)
        {
            PointFluxJacobian_pn(conservVar,normals,fluxJac,
                                efix_StegerWarming,fsw);
        }
        else
        {
            DNekMatSharedPtr PointFJac3D = MemoryManager<DNekMat>
                ::AllocateSharedPtr(nvariables3D, nvariables3D,0.0);
            
            Array<OneD, NekDouble> PointFJac3DData;
            PointFJac3DData = PointFJac3D->GetPtr();

            Array<OneD, NekDouble> PointFwd(nvariables3D,0.0);

            Array<OneD, unsigned int> index(nvariables);

            index[nvariables-1] = 4;
            for(int i=0;i<nvariables-1;i++)
            {
                index[i] = i;
            }

            int nj=0;
            int nk=0;
            for(int j=0; j< nvariables; j++)
            {
                nj = index[j];
                PointFwd[nj] = conservVar[j];
            }
            
            PointFluxJacobian_pn(PointFwd,normals,PointFJac3D,
                                efix_StegerWarming,fsw);

            for(int j=0; j< nvariables; j++)
            {
                nj = index[j];
                for(int k=0; k< nvariables; k++)
                {
                    nk = index[k];
                    fluxJacData[j+k*nConvectiveFields] = 
                            PointFJac3DData[nj+nk*nvariables3D]; 
                }
            }
        }
    }

    // Currently duplacate in compressibleFlowSys
    // if fsw=+-1 calculate the steger-Warming flux vector splitting flux Jacobian
    // if fsw=0   calculate the Jacobian of the exact flux
    // efix is the numerical flux entropy fix parameter
    void CompressibleFlowSystem::PointFluxJacobian_pn(
            const Array<OneD, NekDouble> &Fwd,
            const Array<OneD, NekDouble> &normals,
                  DNekMatSharedPtr       &FJac,
            const NekDouble efix,   const NekDouble fsw)
    {
        int Nvar = 5;
        NekDouble ro,vx,vy,vz,ps,gama,ae ;
        NekDouble a,a2,h,h0,v2,vn,eps,eps2;
        NekDouble nx,ny,nz;
        NekDouble sn,osn,nxa,nya,nza,vna;
        NekDouble l1,l4,l5,al1,al4,al5,x1,x2,x3,y1;
        NekDouble c1,d1,c2,d2,c3,d3,c4,d4,c5,d5;
        NekDouble sml_ssf= 1.0E-12;

        NekDouble fExactorSplt = 2.0-abs(fsw); // if fsw=+-1 calculate 

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
        // TODO:
        // ps = m_eos->GetPressure(rhoL, eL);
        // gama = m_eos->GetGamma();
        ps      = m_varConv->Geteos()->GetPressure(rhoL, eL);
        gama    = m_varConv->Geteos()->GetGamma();

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

    void CompressibleFlowSystem::PointFluxJacobian_pn2D(
            const Array<OneD, NekDouble> &Fwd,
            const Array<OneD, NekDouble> &normals,
                  DNekMatSharedPtr       &FJac,
            const NekDouble efix,   const NekDouble fsw)
    {
        int Nvar = 4;
        NekDouble ro,vx,vy,vz,ps,gama,ae ;
        NekDouble a,a2,h,h0,v2,vn,eps,eps2;
        NekDouble nx,ny,nz;
        NekDouble sn,osn,nxa,nya,nza,vna;
        NekDouble l1,l4,l5,al1,al4,al5,x1,x2,x3,y1;
        NekDouble c1,d1,c2,d2,c3,d3,c4,d4,c5,d5;
        NekDouble sml_ssf= 1.0E-12;

        NekDouble fExactorSplt = 2.0-abs(fsw); // if fsw=+-1 calculate 
        
        NekDouble   rhoL  = Fwd[0];
        NekDouble   rhouL = Fwd[1];
        NekDouble   rhovL = Fwd[2];
        NekDouble   rhowL = 0.0;
        NekDouble   EL    = Fwd[Nvar-1];

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
        // TODO:
        // ps = m_eos->GetPressure(rhoL, eL);
        // gama = m_eos->GetGamma();
        ps      = m_varConv->Geteos()->GetPressure(rhoL, eL);
        gama    = m_varConv->Geteos()->GetGamma();

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
        int ncl4    = (Nvar-1)*Nvar;

        JacData[     ncl0] = c1*y1 - d1*vna + l1;
        JacData[     ncl1] = -c1*vx + d1*nxa;
        JacData[     ncl2] = -c1*vy + d1*nya;
        JacData[     ncl4] = c1;

        c2 = c1*vx + d1*nxa*ae;
        d2 = x3*nxa + d1*vx;
        JacData[ 1 + ncl0] = c2*y1 - d2*vna;
        JacData[ 1 + ncl1] = -c2*vx + d2*nxa + l1;
        JacData[ 1 + ncl2] = -c2*vy + d2*nya;
        JacData[ 1 + ncl4] = c2;

        c3 = c1*vy + d1*nya*ae;
        d3 = x3*nya + d1*vy;
        JacData[ 2 + ncl0] = c3*y1 - d3*vna;
        JacData[ 2 + ncl1] = -c3*vx + d3*nxa;
        JacData[ 2 + ncl2] = -c3*vy + d3*nya + l1;
        JacData[ 2 + ncl4] = c3;

        // c4 = c1*vz + d1*nza*ae;
        // d4 = x3*nza + d1*vz;
        // JacData[ 3 + ncl0] = c4*y1 - d4*vna;
        // JacData[ 3 + ncl1] = -c4*vx + d4*nxa;
        // JacData[ 3 + ncl2] = -c4*vy + d4*nya;
        // JacData[ 3 + ncl3] = -c4*vz + d4*nza + l1;
        // JacData[ 3 + ncl4] = c4;

        c5 = c1*h0 + d1*vna*ae;
        d5 = x3*vna + d1*h0;
        JacData[ Nvar-1 + ncl0] = c5*y1 - d5*vna;
        JacData[ Nvar-1 + ncl1] = -c5*vx + d5*nxa;
        JacData[ Nvar-1 + ncl2] = -c5*vy + d5*nya;
        JacData[ Nvar-1 + ncl4] = c5 + l1;
    }

#endif
    /**
     * @brief Compute the advection terms for the right-hand side
     */
    void CompressibleFlowSystem::DoAdvection(
        const TensorOfArray2D<NekDouble>            &inarray,
        Array<OneD, Array<OneD, NekDouble> >        &outarray,
        const NekDouble                             time,
        const Array<OneD, Array<OneD, NekDouble> >  &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >  &pBwd)
    {
        int nvariables = inarray.num_elements();
        Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);

        m_advObject->Advect(nvariables, m_fields, advVel, inarray,
                            outarray, time, pFwd, pBwd);
    }

    /**
     * @brief Add the diffusions terms to the right-hand side
     */
    void CompressibleFlowSystem::DoDiffusion(
        const TensorOfArray2D<NekDouble>            &inarray,
        Array<OneD, Array<OneD, NekDouble> >        &outarray,
        const Array<OneD, Array<OneD, NekDouble> >  &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >  &pBwd)
    {
        v_DoDiffusion(inarray, outarray, pFwd, pBwd);

        if (m_shockCaptureType != "Off" && m_shockCaptureType != "Physical")
        {
            m_artificialDiffusion->DoArtificialDiffusion(inarray, outarray);
        }
    }

    /**
     * @brief Add the diffusions terms to the right-hand side
     */
    void CompressibleFlowSystem::DoDiffusion_coeff(
        const TensorOfArray2D<NekDouble>              &inarray,
        Array<OneD, Array<OneD, NekDouble> >          &outarray,
        const Array<OneD, Array<OneD, NekDouble> >    &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >    &pBwd,
        const bool                                    flagFreezeJac)
    {
        v_DoDiffusion_coeff(inarray, outarray, pFwd, pBwd,flagFreezeJac);
    }

    void CompressibleFlowSystem::SetBoundaryConditions(
            Array<OneD, Array<OneD, NekDouble> >             &physarray,
            NekDouble                                         time)
    {
        int nTracePts  = GetTraceTotPoints();
        int nvariables = physarray.num_elements();

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

    void CompressibleFlowSystem::SetBoundaryConditionsDeriv(
        const TensorOfArray2D<NekDouble>      &physarray,
        const TensorOfArray3D<NekDouble>      &dervarray,
        NekDouble                             time,
        const TensorOfArray2D<NekDouble>      &pFwd,
        const TensorOfArray3D<NekDouble>      &pDervFwd)
    {
        int nTracePts  = GetTraceTotPoints();
        int nvariables = physarray.num_elements();

        Array<OneD, Array<OneD, NekDouble> > Fwd;
        if(pFwd.num_elements())
        {
            Fwd = pFwd;
        }
        else
        {
            Fwd = Array<OneD, Array<OneD, NekDouble> >(nvariables);
            for (int i = 0; i < nvariables; ++i)
            {
                Fwd[i] = Array<OneD, NekDouble>(nTracePts);
                m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
            }
        }

        TensorOfArray3D<NekDouble> DervFwd;
        if(pDervFwd.num_elements())
        {
            DervFwd = pDervFwd;
        }
        else
        {
            int nDim      = m_fields[0]->GetCoordim(0);
            DervFwd = TensorOfArray3D<NekDouble>(nDim);
            for (int nd = 0; nd < nDim; ++nd)
            {
                DervFwd[nd] = Array<OneD, Array<OneD, NekDouble> > (nvariables);
                for (int i = 0; i < nvariables; ++i)
                {
                    DervFwd[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);
                    m_fields[i]->ExtractTracePhys(dervarray[nd][i], 
                                                  DervFwd[nd][i]);
                }
            }
        }

        if (m_bndConds.size())
        {
            // Loop over user-defined boundary conditions
            for (auto &x : m_bndConds)
            {
                x->ApplyDeriv(Fwd, physarray, DervFwd, dervarray, time);
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
        int nq = physfield[0].num_elements();
        int nVariables = m_fields.num_elements();

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

    void CompressibleFlowSystem::GetFluxVectorMF(
        const Array<OneD, Array<OneD, NekDouble> >  &physfield,
        TensorOfArray3D<NekDouble>                  &flux)
    {
        int nConvectiveFields   = m_fields.num_elements();
        std::shared_ptr<LocalRegions::ExpansionVector> expvect = 
                m_fields[0]->GetExp();
        int ntotElmt            = (*expvect).size();
        LocalRegions::ExpansionSharedPtr pexp = (*expvect)[0];

        Array<OneD, NekDouble> elmtalArray;
        Array<OneD, NekDouble> elmtalFlux;

        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            int nElmtPnt = (*expvect)[nelmt]->GetTotPoints();
            int noffset  = GetPhys_Offset(nelmt);
            int nPntVars = nElmtPnt*nConvectiveFields;

            if(elmtalArray.num_elements() < nPntVars)
            {
                elmtalArray = Array<OneD, NekDouble> (nPntVars, 0.0);
                elmtalFlux  = Array<OneD, NekDouble> (nPntVars, 0.0);
            }
            NekVector<NekDouble> elmtalvector (nPntVars,elmtalArray, eWrapper);
            NekVector<NekDouble> Fluxvector (nPntVars,elmtalFlux, eWrapper);

    
            for(int j = 0; j < nConvectiveFields; j++)
            {   
                Vmath::Vcopy(nElmtPnt,&physfield[j][0]+noffset,1,
                            &elmtalArray[0]+j,nConvectiveFields);
            }

            for(int ndir=0;ndir<m_spacedim;ndir++)
            {
                Fluxvector  = (*m_ElmtFluxJacArray[nelmt][ndir])*elmtalvector;

                for(int j = 0; j < nConvectiveFields; j++)
                {   
                    Vmath::Vcopy(nElmtPnt,&elmtalFlux[0]+j, nConvectiveFields, 
                                &flux[j][ndir][0]+noffset,1);
                }
            }
        }
    }

    void CompressibleFlowSystem::GetFluxVectorTraceMF(
        const Array<OneD, Array<OneD, NekDouble> >      &Fwd,
        const Array<OneD, Array<OneD, NekDouble> >      &Bwd,
        Array<OneD, Array<OneD, NekDouble> >            &flux)
    {
        switch (m_LiniearizationMethod)
        {
        case 2:
            GetFluxVectorTraceMFNum(Fwd, Bwd, flux);
            break;
        case 3:
            GetFluxVectorTraceMFAna(Fwd, Bwd, flux);
            break;
        default:
            GetFluxVectorTraceMFNum(Fwd, Bwd, flux);
            break;
        }
    }

    void CompressibleFlowSystem::GetFluxVectorTraceMFNum(
        const Array<OneD, Array<OneD, NekDouble> >      &Fwd,
        const Array<OneD, Array<OneD, NekDouble> >      &Bwd,
        Array<OneD, Array<OneD, NekDouble> >            &flux)
    {
        int nConvectiveFields   = m_fields.num_elements();
        MultiRegions::ExpListSharedPtr tracelist = m_fields[0]->GetTrace();
        std::shared_ptr<LocalRegions::ExpansionVector> traceExp= 
                tracelist->GetExp();
        int ntotTrac            = (*traceExp).size();
        int nTracePnts          = tracelist->GetTotPoints();

        int nTraceVars          = nTracePnts*nConvectiveFields;
        Array<OneD, NekDouble> elmtalArray(nTraceVars);
        Array<OneD, NekDouble> elmtalFwd(nTraceVars);
        Array<OneD, NekDouble> elmtalBwd(nTraceVars);
        for(int i = 0; i < nConvectiveFields; ++i)
        {
            Vmath::Vcopy(nTracePnts, &Fwd[i][0],1, &elmtalArray[0]+i,
                        nConvectiveFields);
        }   
        NekVector<NekDouble> Varray(nTraceVars,elmtalArray, eWrapper);
        NekVector<NekDouble> VfluxFwd (nTraceVars,elmtalFwd, eWrapper);
        VfluxFwd = (*m_TraceJac[0])*Varray;

        for(int i = 0; i < nConvectiveFields; ++i)
        {
            Vmath::Vcopy(nTracePnts, &Bwd[i][0],1, &elmtalArray[0]+i,
                        nConvectiveFields);
        }   
        NekVector<NekDouble> VfluxBwd (nTraceVars,elmtalBwd, eWrapper);
        VfluxBwd = (*m_TraceJac[1])*Varray;

        Vmath::Vadd(nTraceVars,&elmtalFwd[0],1,&elmtalBwd[0],1,&elmtalFwd[0],1);

        for(int i = 0; i < nConvectiveFields; ++i)
        {
            Vmath::Vcopy(nTracePnts, &elmtalFwd[0]+i,nConvectiveFields, 
                        &flux[i][0],1);
        }  
    }

    void CompressibleFlowSystem::GetFluxVectorTraceMFAna(
        const Array<OneD, Array<OneD, NekDouble> >      &Fwd,
        const Array<OneD, Array<OneD, NekDouble> >      &Bwd,
        Array<OneD, Array<OneD, NekDouble> >            &flux)
    {
        m_advObject->GetRiemann()->Solve(m_spacedim, 
                                         m_MatrixFreeRefFwd, 
                                         m_MatrixFreeRefBwd, 
                                         Fwd, 
                                         Bwd, 
                                         flux);
    }

    /**
     * @brief Return the flux vector for the compressible Euler equations
     * by using the de-aliasing technique.
     *
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >  &physfield,
        TensorOfArray3D<NekDouble>                  &flux)
    {
        int i, j;
        int nq = physfield[0].num_elements();
        int nVariables = m_fields.num_elements();

        // Factor to rescale 1d points in dealiasing
        NekDouble OneDptscale = 2;
        nq = m_fields[0]->Get1DScaledTotPoints(OneDptscale);

        Array<OneD, NekDouble> pressure(nq);
        Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);

        Array<OneD, Array<OneD, NekDouble> > physfield_interp(nVariables);
        TensorOfArray3D<NekDouble> flux_interp(
                                                            nVariables);

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

        // Galerkin project solution back to origianl space
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

            // Galerkin project solution back to origianl space
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
        const TensorOfArray2D<NekDouble>    &inarray,
        Array<OneD, NekDouble>              &tstep)
    {
        int n;
        int nElements = m_fields[0]->GetExpSize();

        // Change value of m_timestep (in case it is set to zero)
        NekDouble tmp = m_timestep;
        m_timestep    = 1.0;

        Array<OneD, NekDouble> cfl(nElements);
        cfl = GetElmtCFLVals();

        // Factors to compute the time-step limit
        NekDouble alpha     = MaxTimeStepEstimator();

        // Loop over elements to compute the time-step limit for each element
        for(n = 0; n < nElements; ++n)
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
        const TensorOfArray2D<NekDouble> &inarray)
    {
        int nElements = m_fields[0]->GetExpSize();
        Array<OneD, NekDouble> tstep (nElements, 0.0);

        GetElmtTimeStep(inarray, tstep);

        // Get the minimum time-step limit and return the time-step
        NekDouble TimeStep = Vmath::Vmin(nElements, tstep, 1);
        m_comm->AllReduce(TimeStep, LibUtilities::ReduceMin);

        NekDouble tmp = m_timestep;
        m_timestep    = TimeStep;

        Array<OneD, NekDouble> cflNonAcoustic(nElements,0.0);
        cflNonAcoustic = GetElmtCFLVals(false);

        // Get the minimum time-step limit and return the time-step
        NekDouble MaxcflNonAcoustic = Vmath::Vmax(nElements, cflNonAcoustic, 1);
        m_comm->AllReduce(MaxcflNonAcoustic, LibUtilities::ReduceMax);

        m_cflNonAcoustic = MaxcflNonAcoustic;
        m_timestep = tmp;

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
        EquationSystem::v_SetInitialConditions(initialtime, false);

        // insert white noise in initial condition
        NekDouble Noise;
        int phystot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> noise(phystot);

        m_session->LoadParameter("Noise", Noise,0.0);
        int m_nConvectiveFields =  m_fields.num_elements();

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
    Array<OneD, NekDouble> CompressibleFlowSystem::v_GetMaxStdVelocity(
        const NekDouble SpeedSoundFactor)
    {
        int nTotQuadPoints = GetTotPoints();
        int n_element      = m_fields[0]->GetExpSize();
        int expdim         = m_fields[0]->GetGraph()->GetMeshDimension();
        int nfields        = m_fields.num_elements();
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

        for(int el = 0; el < n_element; ++el)
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
            if(metricInfo->GetGtype() == SpatialDomains::eDeformed)
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
                          SpeedSoundFactor * std::abs(stdSoundSpeed[j][offset + i]);
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
            Array<OneD, Array<OneD, NekDouble> > tmp(m_fields.num_elements());

            for (int i = 0; i < m_fields.num_elements(); ++i)
            {
                tmp[i] = m_fields[i]->GetPhys();
            }

            Array<OneD, NekDouble> pressure(nPhys), temperature(nPhys);
            Array<OneD, NekDouble> entropy(nPhys);
            Array<OneD, NekDouble> soundspeed(nPhys), mach(nPhys);
            Array<OneD, NekDouble> sensor(nPhys), SensorKappa(nPhys);

            m_varConv->GetPressure  (tmp, pressure);
            m_varConv->GetTemperature(tmp, temperature);
            m_varConv->GetEntropy   (tmp, entropy);
            m_varConv->GetSoundSpeed(tmp, soundspeed);
            m_varConv->GetMach      (tmp, soundspeed, mach);

            Array<OneD, Array<OneD, NekDouble> > velocities(m_spacedim);
            for (int i=0;i<m_spacedim;i++)
            {
                velocities[i] = Array<OneD, NekDouble> (nPhys);
            }
            m_varConv->GetVelocityVector(tmp,velocities);

            int sensorOffset;
            m_session->LoadParameter ("SensorOffset", sensorOffset, 1);
            m_varConv->GetSensor (m_fields[0], tmp, sensor, SensorKappa,
                                    sensorOffset);

            Array<OneD, NekDouble> pFwd(nCoeffs), TFwd(nCoeffs);
            Array<OneD, NekDouble> sFwd(nCoeffs);
            Array<OneD, NekDouble> aFwd(nCoeffs), mFwd(nCoeffs);
            Array<OneD, NekDouble> sensFwd(nCoeffs);

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

            Array<OneD, NekDouble> uFwd(nCoeffs);
            m_fields[0]->FwdTrans_IterPerExp(velocities[0],uFwd);
            variables.push_back  ("u");
            fieldcoeffs.push_back(uFwd);

            if(m_spacedim>1)
            {
                Array<OneD, NekDouble> vFwd(nCoeffs);
                variables.push_back  ("v");
                m_fields[0]->FwdTrans_IterPerExp(velocities[1],vFwd);
                fieldcoeffs.push_back(vFwd);
            }
            if(m_spacedim>2)
            {
                Array<OneD, NekDouble> wFwd(nCoeffs);
                variables.push_back  ("w");
                m_fields[0]->FwdTrans_IterPerExp(velocities[2],wFwd);
                fieldcoeffs.push_back(wFwd);
            }

            if (m_artificialDiffusion)
            {
                Array<OneD, NekDouble> sensorFwd(nCoeffs);
                // reuse pressure
                m_artificialDiffusion->GetArtificialViscosity(tmp, pressure);
                m_fields[0]->FwdTrans_IterPerExp(pressure,   sensorFwd);

                variables.push_back  ("ArtificialVisc");
                fieldcoeffs.push_back(sensorFwd);
            }
        }
    }

    void CompressibleFlowSystem::v_GetViscousSymmtrFluxConservVar(
            const int                                   nConvectiveFields,
            const int                                   nSpaceDim,
            const Array<OneD, Array<OneD, NekDouble> >  &inaverg,
            const Array<OneD, Array<OneD, NekDouble > > &inarray,
            TensorOfArray3D<NekDouble>                  &outarray,
            Array< OneD, int >                          &nonZeroIndex,
            const Array<OneD, Array<OneD, NekDouble> >  &normals)
    {
        ASSERTL0(false, "v_GetViscousSymmtrFluxConservVar not coded");
    }
    
    void CompressibleFlowSystem::v_SteadyStateResidual(
                int                         step, 
                Array<OneD, NekDouble>      &L2)
    {
        const int nPoints = GetTotPoints();
        const int nFields = m_fields.num_elements();
        Array<OneD, Array<OneD, NekDouble> > rhs (nFields);
        Array<OneD, Array<OneD, NekDouble> > inarray (nFields);
        for (int i = 0; i < nFields; ++i)
        {
            rhs[i] =   Array<OneD, NekDouble> (nPoints,0.0);
            inarray[i] =   m_fields[i]->UpdatePhys();
        }

        DoOdeRhs(inarray,rhs,m_time);

        // Holds L2 errors.
        Array<OneD, NekDouble> tmp;
        Array<OneD, NekDouble> RHSL2    (nFields);
        Array<OneD, NekDouble> residual(nFields);

        for (int i = 0; i < nFields; ++i)
        {
            tmp = rhs[i];

            Vmath::Vmul(nPoints, tmp, 1, tmp, 1, tmp, 1);
            residual[i] = Vmath::Vsum(nPoints, tmp, 1);
        }

        m_comm->AllReduce(residual , LibUtilities::ReduceSum);

        NekDouble onPoints = 1.0/NekDouble(nPoints);
        for (int i = 0; i < nFields; ++i)
        {
            L2[i] = sqrt(residual[i]*onPoints);
        }
    }

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
    void CompressibleFlowSystem::v_GetFluxDerivJacDirctn(
            const MultiRegions::ExpListSharedPtr &explist,
            const TensorOfArray2D<NekDouble>     &normals,
            const int                            nDervDir,
            const TensorOfArray2D<NekDouble>     &inarray,
            TensorOfArray5D<NekDouble>           &ElmtJacArray,
            const int                            nfluxDir)
    {
        ASSERTL0(false, "v_GetFluxDerivJacDirctn not coded");
    }

    void CompressibleFlowSystem::v_GetFluxDerivJacDirctnElmt(
            const int                                  nConvectiveFields,
            const int                                  nElmtPnt,
            const int                                  nDervDir,
            const Array<OneD, Array<OneD, NekDouble> > &locVars,
            const Array<OneD, NekDouble>               &locmu,
            const Array<OneD, Array<OneD, NekDouble> > &locnormal,
            DNekMatSharedPtr                           &wspMat,
            Array<OneD, Array<OneD, NekDouble> >       &PntJacArray)
    {
        ASSERTL0(false, "v_GetFluxDerivJacDirctn not coded");
    }
    
    void CompressibleFlowSystem::v_GetFluxDerivJacDirctn(
            const MultiRegions::ExpListSharedPtr &explist,
            const TensorOfArray2D<NekDouble>     &normals,
            const int                            nDervDir,
            const TensorOfArray2D<NekDouble>     &inarray,
            TensorOfArray2D<DNekMatSharedPtr>    &ElmtJac)
    {
        ASSERTL0(false, "v_GetFluxDerivJacDirctn not coded");
    }

    // void CompressibleFlowSystem::v_GetFluxDerivJacDirctn(
    //         const MultiRegions::ExpListSharedPtr                            &explist,
    //         const int                                                       nFluxDir,
    //         const int                                                       nDervDir,
    //         const TensorOfArray2D<NekDouble>                &inarray,
    //               TensorOfArray2D<DNekMatSharedPtr>               &ElmtJac)
    // {
    //     ASSERTL0(false, "v_GetFluxDerivJacDirctn not coded");
    // }

    void CompressibleFlowSystem::v_GetDiffusionFluxJacPoint(
            const Array<OneD, NekDouble>     &conservVar, 
            const TensorOfArray2D<NekDouble> &conseDeriv, 
            const NekDouble                  mu,
            const NekDouble                  DmuDT,
            const Array<OneD, NekDouble>     &normals,
            DNekMatSharedPtr                 &fluxJac)
    {
        ASSERTL0(false, "v_GetDiffusionFluxJacPoint not coded");
    }

    void CompressibleFlowSystem::v_MinusDiffusionFluxJacDirctn(
            const int                        nDirctn,
            const TensorOfArray2D<NekDouble> &inarray,
            const TensorOfArray3D<NekDouble> &qfields,
            TensorOfArray5D<NekDouble>       &ElmtJacArray)
    {
        ASSERTL0(false, "v_MinusDiffusionFluxJacDirctn not coded");
    }

    void CompressibleFlowSystem::v_MinusDiffusionFluxJacDirctnMat(
            const int                               nDirctn,
            const TensorOfArray2D<NekDouble>        &inarray,
            const TensorOfArray3D<NekDouble>        &qfields,
            TensorOfArray2D<DNekBlkMatSharedPtr>    &ElmtFluxJacArray)
    {
        ASSERTL0(false, "v_MinusDiffusionFluxJacDirctnMat not coded");
    }

    void CompressibleFlowSystem::v_MinusDiffusionFluxJacDirctnElmt(
            const int                                  nConvectiveFields,
            const int                                  nElmtPnt,
            const Array<OneD, Array<OneD, NekDouble> > &locVars,
            const TensorOfArray3D<NekDouble>           &locDerv,
            const Array<OneD, NekDouble>               &locmu,
            const Array<OneD, NekDouble>               &locDmuDT,
            const Array<OneD, NekDouble>               &normals,
            DNekMatSharedPtr                           &wspMat,
            Array<OneD, Array<OneD, NekDouble> >       &PntJacArray)
    {
        ASSERTL0(false, "not coded");
    }

#endif

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
                for(int i = 0; i < exp3D->GetNedges(); ++i)
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
                for(int i = 0; i < exp2D->GetNedges(); ++i)
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
    
    void CompressibleFlowSystem::v_DoOdeProjection1(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
              const NekDouble                             time)
    {
        DoOdeProjection(inarray,outarray,time);
    } 

    void CompressibleFlowSystem::v_DoOdeRhs1(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
              const NekDouble                             time)
    {
        DoOdeRhs(inarray,outarray,time);
    }

    void CompressibleFlowSystem::preconditioner_MultiLevel_coeff(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&outarray,
                  const bool             &UnusedFlag)
    { 
        int Level=1;
        m_EqdriverOperator.MultiLevel(inarray, outarray, m_UpDateOperatorflag, 
            Level, m_TimeIntegtSol_k, m_BndEvaluateTime, m_TimeIntegLambda);
        m_UpDateOperatorflag = false;
    }

    void CompressibleFlowSystem::v_MutilLvlResidual(
        const TensorOfArray1D<NekDouble>    &inarray,
        TensorOfArray1D<NekDouble>          &outarray,
        TensorOfArray1D<NekDouble>          &outRes,
        const int                           level)
    {
        //Ax: cannot use MatrixMultiply_MatrixFree_coeff in low level.
        //Because it reuse TimeIntegration_Res_n, which does not calculate
        if(level>0)
        {
            MatrixMultiply_JacobianFree_coeff_noSource(outarray,outRes);
        }
        else
        {
            //Avoid repeating calculating Residual_k in highest level
            MatrixMultiply_JacobianFree_coeff(outarray,outRes);
        }

        //b-Ax:Ok
        Vmath::Vsub(inarray.num_elements(),inarray,1,outRes,1,outRes,1);
    }

    void CompressibleFlowSystem::v_MutilLvlPresmooth(
        const TensorOfArray1D<NekDouble>    &inarray,
        TensorOfArray1D<NekDouble>          &outarray)
    {
        preconditioner_BlkSOR_coeff(inarray,outarray,false);
    }

    void CompressibleFlowSystem::v_MutilLvlPostsmooth(
        const TensorOfArray1D<NekDouble>    &inarray,
        TensorOfArray1D<NekDouble>          &outarray)
    {
        preconditioner_BlkSOR_coeff(inarray,outarray,true);
    }

    void CompressibleFlowSystem::v_MutilLvlLowestLvlSolver(
        const TensorOfArray1D<NekDouble>    &inarray,
        TensorOfArray1D<NekDouble>          &outarray)
    {
        preconditioner_BlkSOR_coeff(inarray,outarray,false);
    }

    void CompressibleFlowSystem::v_MutilLvlRestriction1D(
        const TensorOfArray1D<NekDouble>    &inarray,
        TensorOfArray1D<NekDouble>          &outarray)
    {
        RestrictResidual(m_RestrictionMatrix,inarray,outarray);
    }

    void CompressibleFlowSystem::v_MutilLvlRestriction2D(
        const TensorOfArray2D<NekDouble>    &inarray,
        TensorOfArray2D<NekDouble>          &outarray)
    {
        RestrictSolution(m_RestrictionMatrix,inarray,outarray);
    }

    void CompressibleFlowSystem::v_MutilLvlProlong(
        const TensorOfArray1D<NekDouble>    &inarray,
        TensorOfArray1D<NekDouble>          &outarray)
    {
        ProlongateSolution(m_ProlongationMatrix,inarray,outarray);
    }

    void CompressibleFlowSystem::v_MutilLvlUpdateMultilvlMatrix(
        const TensorOfArray2D<NekDouble>    &lowLevelSolution,
        const NekDouble                     time,
        const NekDouble                     dtlamda)
    {
        v_CalculateNextLevelPreconditioner(
            lowLevelSolution, time, dtlamda);
    }

    void  CompressibleFlowSystem::v_CalculateNextLevelPreconditioner(
        const TensorOfArray2D<NekDouble>    &inarrayCoeff,
        const NekDouble                     time,
        const NekDouble                     lambda)
    {
        int nVariables=m_fields.num_elements();
        m_BndEvaluateTime=time;
        m_TimeIntegLambda=lambda;
        Array<OneD,Array<OneD,NekDouble>> inarray(nVariables);
        m_TimeIntegtSol_k=Array<OneD,Array<OneD,NekDouble>> (nVariables);
        m_LowLevelSysEquatResid=Array<OneD,Array<OneD,NekDouble>> (nVariables);
        int npoints=GetNpoints();
        int ncoeffs=GetNcoeffs();
        for(int i=0;i<nVariables;i++)
        {
            inarray[i]=Array<OneD,NekDouble> (npoints,0.0);
            m_LowLevelSysEquatResid[i]=Array<OneD,NekDouble> (ncoeffs,0.0);
            m_TimeIntegtSol_k[i]=inarrayCoeff[i];
            m_fields[i]->BwdTrans(inarrayCoeff[i],inarray[i]);
        }
        NonlinSysEvaluator_coeff_noSource(m_TimeIntegtSol_k,m_LowLevelSysEquatResid);
        
        //This part of codes coping from DoImplicit_Coeff so that get m_magnitude
        int ntotal         = nVariables*GetNcoeffs();
        int nElements      = m_fields[0]->GetExpSize();

        LibUtilities::CommSharedPtr v_Comm  = m_fields[0]->GetComm()->GetRowComm();

        bool l_verbose      = m_session->DefinesCmdLineArgument("verbose");
        bool l_root=false;
        if(0==v_Comm->GetRank())
        {
            l_root =true;
        }
        unsigned int ntotalGlobal     = ntotal;
        v_Comm->AllReduce(ntotalGlobal, Nektar::LibUtilities::ReduceSum);
        unsigned int ntotalDOF        = ntotalGlobal/nVariables;
        NekDouble ototalGlobal  = 1.0/ntotalGlobal;
        NekDouble ototalDOF     = 1.0/ntotalDOF;
        UpdateSoltnRefNorms(inarrayCoeff, ototalDOF);
        CalPrecMat(inarray,m_BndEvaluateTime,m_TimeIntegLambda);
    }
   
   //To Do: check the Restict Residual and Restrict Solution. they seem different
    void CompressibleFlowSystem::RestrictResidual(
        const Array<OneD,DNekMatSharedPtr>         &RestrictionMatrix,
        const Array<OneD, NekDouble>               &inarray,
              Array<OneD, NekDouble>               &outarray)
    {
        std::shared_ptr<LocalRegions::ExpansionVector> pexp = m_fields[0]->GetExp();
        int nVariables=m_fields.num_elements();
        int nElmts=m_fields[0]->GetExpSize();
        
        int RestrictedArrayOffset=0;
        for(int i=0;i<nVariables;i++)
        {
            int VariableOffset=i*GetNcoeffs();
            for(int j=0;j<nElmts;j++)
            {
                int nElmtCoeffs=(*pexp)[j]->GetNcoeffs();
                int nElmtOffset=m_fields[0]->GetCoeff_Offset(j);
                int nRestrictedElmtCoeffs=RestrictionMatrix[j]->GetRows();
                Array<OneD, NekDouble> Residualtmp(nElmtCoeffs,&inarray[VariableOffset+nElmtOffset]);  
                DNekVec Vin(nElmtCoeffs,Residualtmp,eCopy);
                Array<OneD,NekDouble> RestrictedResidualtmp(nRestrictedElmtCoeffs,&outarray[RestrictedArrayOffset]);
                DNekVec Vout(RestrictedResidualtmp,eWrapper);
                Vout=(*RestrictionMatrix[j])*Vin;
                Vmath::Vcopy(nRestrictedElmtCoeffs,&RestrictedResidualtmp[0],1,&outarray[RestrictedArrayOffset],1);
                RestrictedArrayOffset=RestrictedArrayOffset+nRestrictedElmtCoeffs;
            }
        }
    }

    //To Do: check the Restict Residual and Restrict Solution. they seem different
    void CompressibleFlowSystem::RestrictSolution(
        const Array<OneD,DNekMatSharedPtr>         &RestrictionMatrix,
        const Array<OneD, const Array<OneD, NekDouble>>  &inarray,
              Array<OneD, Array<OneD, NekDouble>>  &outarray)
    {
        std::shared_ptr<LocalRegions::ExpansionVector> pexp = m_fields[0]->GetExp();
        int nVariables=m_fields.num_elements();
        int nElmts=m_fields[0]->GetExpSize();

        for(int i=0;i<nVariables;i++)
        {
            int RestrictedArrayOffset=0;
            int VariableOffset=i*GetNcoeffs();
            for(int j=0;j<nElmts;j++)
            {
                int nElmtCoeffs=(*pexp)[j]->GetNcoeffs();
                int nElmtOffset=m_fields[0]->GetCoeff_Offset(j);
                int nRestrictedElmtCoeffs=RestrictionMatrix[j]->GetRows();
                Array<OneD, NekDouble> Solutiontmp(nElmtCoeffs,&inarray[i][nElmtOffset]);  
                DNekVec Vin(nElmtCoeffs,Solutiontmp,eCopy);
                Array<OneD,NekDouble> RestrictedSolutiontmp(nRestrictedElmtCoeffs,&outarray[i][RestrictedArrayOffset]);
                DNekVec Vout(RestrictedSolutiontmp,eWrapper);
                Vout=(*RestrictionMatrix[j])*Vin;
                Vmath::Vcopy(nRestrictedElmtCoeffs,&RestrictedSolutiontmp[0],1,&outarray[i][RestrictedArrayOffset],1);
                RestrictedArrayOffset=RestrictedArrayOffset+nRestrictedElmtCoeffs;
            }
        }
    }

    void CompressibleFlowSystem::ProlongateSolution(
        const Array<OneD, DNekMatSharedPtr>        &ProlongationMatrix,
        const Array<OneD, NekDouble>               &inarray,
              Array<OneD, NekDouble>               &outarray)
    {
        std::shared_ptr<LocalRegions::ExpansionVector> pexp = m_fields[0]->GetExp();
        int nVariables=m_fields.num_elements();
        int nElmts=m_fields[0]->GetExpSize();

        int ProlongatedArrayOffset=0;
        for(int i=0;i<nVariables;i++)
        {
            int VariableOffset=i*GetNcoeffs();
            for(int j=0;j<nElmts;j++)
            {
                int nElmtCoeffs=(*pexp)[j]->GetNcoeffs();
                int nElmtOffset=m_fields[0]->GetCoeff_Offset(j);
                int nProlongatedElmtCoeffs=ProlongationMatrix[j]->GetColumns();
                Array<OneD, NekDouble> Solutiontmp(nProlongatedElmtCoeffs,&inarray[ProlongatedArrayOffset]);  
                DNekVec Vin(nProlongatedElmtCoeffs,Solutiontmp,eCopy);
                Array<OneD,NekDouble> ProlongatedSolutiontmp(nElmtCoeffs,&outarray[VariableOffset+nElmtOffset]);
                DNekVec Vout(nElmtCoeffs,ProlongatedSolutiontmp,eWrapper);
                Vout=(*ProlongationMatrix[j])*Vin;
                Vmath::Vcopy(nElmtCoeffs,&ProlongatedSolutiontmp[0],1,&outarray[VariableOffset+nElmtOffset],1);  
                ProlongatedArrayOffset=ProlongatedArrayOffset+nProlongatedElmtCoeffs;
            }
        }
    }

    void CompressibleFlowSystem::PrintArray(Array<OneD, NekDouble> &Array)
    {
        int nrows                = Array.num_elements();

        for (int i = 0; i < nrows; i++)
        {

            cout << setprecision(1)  << i << "     " << setprecision(16) << Array[i] << endl;
        }
    }

    void CompressibleFlowSystem::OutputArray(Array<OneD, NekDouble> &Array)
    {
        int nrows = Array.num_elements();

        ofstream outfile1;
        outfile1.open("./Array.txt");
        for (int i = 0; i < nrows; i++)
        {

            outfile1 << setprecision(1)  << i+1 << "     " << setprecision(16) << Array[i] << endl;
        }
        outfile1.close();
    }
    
    void CompressibleFlowSystem::OutputConstArray(const Array<OneD, NekDouble> &Array)
    {
        int nrows = Array.num_elements();

        ofstream outfile1;
        outfile1.open("./ConstArray.txt");
        for (int i = 0; i < 10; i++)
        {

            outfile1 << setprecision(1)  << i+1 << "     " << setprecision(16) << Array[i] << endl;
        }
        outfile1.close();
    }

    void CompressibleFlowSystem::PrintMatrix(DNekMatSharedPtr &Matrix)
    {
        int nrows                = Matrix->GetRows();
        int ncols                = Matrix->GetColumns();
        MatrixStorage matStorage = Matrix->GetStorageType();
        for (int i = 0; i < nrows; i++)
        {

            for (int j = 0; j < ncols; j++)
            {
                cout << setprecision(1)  << i  << "    " << j 
                        << "     " << setprecision(16) << (*Matrix)(i, j) << endl;
            }
        
        }
    }

    void CompressibleFlowSystem::OutputMatrix(DNekMatSharedPtr &Matrix)
    {
        int rows = Matrix->GetRows();
        int cols = Matrix->GetColumns();
        ofstream outfile1;
        outfile1.open("./Matrix.txt");
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                outfile1 << i+1<< "     " << j+1  << "    "
                        << std::setprecision(16) << (*Matrix)(i, j) << endl;
            }
        }
        outfile1.close();
    }

    void CompressibleFlowSystem::OutputConstMatrix(const DNekMatSharedPtr &Matrix)
    {
        int rows = Matrix->GetRows();
        int cols = Matrix->GetColumns();
        ofstream outfile1;
        outfile1.open("./ConstMatrix.txt");
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                outfile1 << i+1<< "     " << j+1  << "    "
                        << std::setprecision(16) << (*Matrix)(i, j) << endl;
            }
        }
        outfile1.close();
    }
}
