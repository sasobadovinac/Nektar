///////////////////////////////////////////////////////////////////////////////
//
// File VelocityCorrection.cpp
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
// Description: Velocity Correction Scheme for the Incompressible
// Navier Stokes equations
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/Core/Misc.h>

#include <MultiRegions/GlobalLinSysStaticCondVec.h>
#include <MultiRegions/GlobalLinSysIterativeStaticCondVec.h>
#include <boost/algorithm/string.hpp>

using namespace std;

namespace Nektar
{
    using namespace MultiRegions;

    string VelocityCorrectionScheme::className = 
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "VelocityCorrectionScheme", 
            VelocityCorrectionScheme::create);

    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    VelocityCorrectionScheme::VelocityCorrectionScheme(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession),
          IncNavierStokes(pSession),
          m_varCoeffLap(StdRegions::NullVarCoeffMap)
    {
        
    }

    void VelocityCorrectionScheme::v_InitObject()
    {
        int n;
        
        IncNavierStokes::v_InitObject();
        m_explicitDiffusion = false;

        // Set m_pressure to point to last field of m_fields;
        if (boost::iequals(m_session->GetVariable(m_fields.num_elements()-1), "p"))
        {
            m_nConvectiveFields = m_fields.num_elements()-1;
            m_pressure = m_fields[m_nConvectiveFields];
        }
        else
        {
            ASSERTL0(false,"Need to set up pressure field definition");
        }

        // Determine diffusion coefficients for each field
        m_diffCoeff = Array<OneD, NekDouble> (m_nConvectiveFields, m_kinvis);
        for (n = 0; n < m_nConvectiveFields; ++n)
        {
            std::string varName = m_session->GetVariable(n);
            if ( m_session->DefinesFunction("DiffusionCoefficient", varName))
            {
                LibUtilities::EquationSharedPtr ffunc
                    = m_session->GetFunction("DiffusionCoefficient", varName);
                m_diffCoeff[n] = ffunc->Evaluate();
            }
        }

        // creation of the extrapolation object
        if(m_equationType == eUnsteadyNavierStokes)
        {
            std::string vExtrapolation = v_GetExtrapolateStr();

            if (m_session->DefinesSolverInfo("Extrapolation"))
            {
                vExtrapolation = v_GetSubSteppingExtrapolateStr(
                                 m_session->GetSolverInfo("Extrapolation"));
            }

            m_extrapolation = GetExtrapolateFactory().CreateInstance(
                vExtrapolation,
                m_session,
                m_fields,
                m_pressure,
                m_velocity,
                m_advObject);
        }

        // Integrate only the convective fields
        for (n = 0; n < m_nConvectiveFields; ++n)
        {
            m_intVariables.push_back(n);
        }
        
        m_session->MatchSolverInfo("SpectralVanishingViscosity",
                                   "PowerKernel", m_useSpecVanVisc, false);

        if(m_useSpecVanVisc)
        {
            m_useHomo1DSpecVanVisc = true;
        }
        else
        {
            m_session->MatchSolverInfo("SpectralVanishingViscositySpectralHP",
                                       "PowerKernel", m_useSpecVanVisc, false);
        }

        if(m_useSpecVanVisc)
        {
            m_IsSVVPowerKernel = true;
        }
        else
        {
            m_session->MatchSolverInfo("SpectralVanishingViscosity","DGKernel",
                                       m_useSpecVanVisc, false);
            if(m_useSpecVanVisc)
            {
                m_useHomo1DSpecVanVisc = true;
            }
            else
            {
                m_session->MatchSolverInfo("SpectralVanishingViscositySpectralHP",
                                           "DGKernel", m_useSpecVanVisc, false);
            }
            
            if(m_useSpecVanVisc)
            {
                m_IsSVVPowerKernel = false;
            }
        }

        //set up varcoeff kernel if PowerKernel or DG is specified
        if(m_useSpecVanVisc)
        {
            Array<OneD, Array<OneD, NekDouble> > SVVVelFields = NullNekDoubleArrayofArray;
            if(m_session->DefinesFunction("SVVVelocityMagnitude"))
            {
                if (m_comm->GetRank() == 0)
                {
                    cout << "Seting up SVV velocity from "
                        "SVVVelocityMagnitude section in session file" << endl;
                }
                int nvel = m_velocity.num_elements();
                int phystot = m_fields[0]->GetTotPoints();
                SVVVelFields = Array<OneD, Array<OneD, NekDouble> >(nvel);
                vector<string> vars;
                for(int i = 0; i < nvel; ++i)
                {
                    SVVVelFields[i] = Array<OneD, NekDouble>(phystot);
                    vars.push_back(m_session->GetVariable(m_velocity[i]));
                }
                    
                // Load up files into  m_fields;
                GetFunction("SVVVelocityMagnitude")
                    ->Evaluate(vars,SVVVelFields);
            }

            m_svvVarDiffCoeff = Array<OneD, NekDouble>(m_fields[0]->GetNumElmts());
            SVVVarDiffCoeff(1.0,m_svvVarDiffCoeff,SVVVelFields);
            m_session->LoadParameter("SVVDiffCoeff",  m_sVVDiffCoeff,  1.0);
        }
        else
        {
            m_svvVarDiffCoeff = NullNekDouble1DArray;
            m_session->LoadParameter("SVVDiffCoeff",  m_sVVDiffCoeff,  0.1);
        }

        // Load parameters for Spectral Vanishing Viscosity
        if(m_useSpecVanVisc == false)
        {
            m_session->MatchSolverInfo("SpectralVanishingViscosity","True",
                                       m_useSpecVanVisc, false);
            if(m_useSpecVanVisc == false)
            {
                m_session->MatchSolverInfo("SpectralVanishingViscosity","ExpKernel",
                                           m_useSpecVanVisc, false);
            }
            m_useHomo1DSpecVanVisc = m_useSpecVanVisc;

            if(m_useSpecVanVisc == false)
            {
                m_session->MatchSolverInfo("SpectralVanishingViscositySpectralHP","True",
                                           m_useSpecVanVisc, false);
                if(m_useSpecVanVisc == false)
                {
                    m_session->MatchSolverInfo("SpectralVanishingViscositySpectralHP","ExpKernel",
                                               m_useSpecVanVisc, false);
                }
            }
        }


        // Case of only Homo1D kernel
        if(m_useSpecVanVisc == false)
        {
            m_session->MatchSolverInfo("SpectralVanishingViscosityHomo1D",
                                "True", m_useHomo1DSpecVanVisc, false);
            if(m_useHomo1DSpecVanVisc == false)
            {
                m_session->MatchSolverInfo("SpectralVanishingViscosityHomo1D",
                                       "ExpKernel", m_useHomo1DSpecVanVisc, false);
            }
        }
        else
        {
            bool testForFalse;
            // Case where Homo1D is turned off but has been turned on
            // impliictly by SpectralVanishingViscosity solver info
            m_session->MatchSolverInfo("SpectralVanishingViscosityHomo1D",
                                "False", testForFalse, false);
            if(testForFalse)
            {
                m_useHomo1DSpecVanVisc = false;
            }            
        }

        m_session->LoadParameter("SVVCutoffRatio",m_sVVCutoffRatio,0.75);
        m_session->LoadParameter("SVVCutoffRatioHomo1D",m_sVVCutoffRatioHomo1D,m_sVVCutoffRatio);
        m_session->LoadParameter("SVVDiffCoeffHomo1D",  m_sVVDiffCoeffHomo1D,  m_sVVDiffCoeff);

        if(m_HomogeneousType == eHomogeneous1D)
        {
            ASSERTL0(m_nConvectiveFields > 2,
                "Expect to have three velocity fields with homogenous expansion");

            if(m_useHomo1DSpecVanVisc)
            {
                Array<OneD, unsigned int> planes;
                planes = m_fields[0]->GetZIDs();

                int num_planes = planes.num_elements();
                Array<OneD, NekDouble> SVV(num_planes,0.0);
                NekDouble fac;
                int kmodes = m_fields[0]->GetHomogeneousBasis()->GetNumModes();
                int pstart;

                pstart = m_sVVCutoffRatioHomo1D*kmodes;
                
                for(n = 0; n < num_planes; ++n)
                {
                    if(planes[n] > pstart)
                    {
                        fac = (NekDouble)((planes[n] - kmodes)*(planes[n] - kmodes))/
                            ((NekDouble)((planes[n] - pstart)*(planes[n] - pstart)));
                        SVV[n] = m_sVVDiffCoeffHomo1D*exp(-fac)/m_kinvis;
                    }
                }

                for(int i = 0; i < m_velocity.num_elements(); ++i)
                {
                    m_fields[m_velocity[i]]->SetHomo1DSpecVanVisc(SVV);
                }
            }
        }

        m_session->MatchSolverInfo("SmoothAdvection", "True",
                                    m_SmoothAdvection, false);

        // set explicit time-intregration class operators
        m_ode.DefineOdeRhs(
            &VelocityCorrectionScheme::EvaluateAdvection_SetPressureBCs, this);

        m_extrapolation->SubSteppingTimeIntegration(
            m_intScheme->GetIntegrationMethod(), m_intScheme);
        m_extrapolation->GenerateHOPBCMap(m_session);
        
        // set implicit time-intregration class operators
        m_ode.DefineImplicitSolve(
            &VelocityCorrectionScheme::SolveUnsteadyStokesSystem,this);
    }
    
    /**
     * Destructor
     */
    VelocityCorrectionScheme::~VelocityCorrectionScheme(void)
    {        
    }
    
    /**
     * 
     */
    void VelocityCorrectionScheme::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        AdvectionSystem::v_GenerateSummary(s);
        SolverUtils::AddSummaryItem(s,
                "Splitting Scheme", "Velocity correction (strong press. form)");

        if (m_extrapolation->GetSubStepIntegrationMethod() !=
            LibUtilities::eNoTimeIntegrationMethod)
        {
            SolverUtils::AddSummaryItem(s, "Substepping", 
                             LibUtilities::TimeIntegrationMethodMap[
                              m_extrapolation->GetSubStepIntegrationMethod()]);
        }

        string dealias = m_homogen_dealiasing ? "Homogeneous1D" : "";
        if (m_specHP_dealiasing)
        {
            dealias += (dealias == "" ? "" : " + ") + string("spectral/hp");
        }
        if (dealias != "")
        {
            SolverUtils::AddSummaryItem(s, "Dealiasing", dealias);
        }


        string smoothing = m_useSpecVanVisc ? "spectral/hp" : "";
        if (smoothing != "")
        {
            if(m_svvVarDiffCoeff == NullNekDouble1DArray)
            {
                SolverUtils::AddSummaryItem(
                   s, "Smoothing-SpecHP", "SVV (" + smoothing +
                   " Exp Kernel(cut-off = "
                   + boost::lexical_cast<string>(m_sVVCutoffRatio)
                   + ", diff coeff = "
                   + boost::lexical_cast<string>(m_sVVDiffCoeff)+"))");
            }
            else
            {
                if(m_IsSVVPowerKernel)
                {
                    SolverUtils::AddSummaryItem(
                       s, "Smoothing-SpecHP", "SVV (" + smoothing +
                       " Power Kernel (Power ratio ="
                       + boost::lexical_cast<string>(m_sVVCutoffRatio)
                       + ", diff coeff = "
                       + boost::lexical_cast<string>(m_sVVDiffCoeff)+"*Uh/p))");
                }
                else
                {
                    SolverUtils::AddSummaryItem(
                       s, "Smoothing-SpecHP", "SVV (" + smoothing +
                       " DG Kernel (diff coeff = "
                       + boost::lexical_cast<string>(m_sVVDiffCoeff)+"*Uh/p))");

                }
            }
            
        }

        if (m_useHomo1DSpecVanVisc && (m_HomogeneousType == eHomogeneous1D))
        {
            SolverUtils::AddSummaryItem(
                  s, "Smoothing-Homo1D", "SVV (Homogeneous1D - Exp Kernel(cut-off = "
                  + boost::lexical_cast<string>(m_sVVCutoffRatioHomo1D)
                  + ", diff coeff = "
                  + boost::lexical_cast<string>(m_sVVDiffCoeffHomo1D)+"))");
        }

    }

    /**
     * 
     */
    void VelocityCorrectionScheme::v_DoInitialise(void)
    {
        AdvectionSystem::v_DoInitialise();

        // Set up Field Meta Data for output files
        m_fieldMetaDataMap["Kinvis"]   =
                boost::lexical_cast<std::string>(m_kinvis);
        m_fieldMetaDataMap["TimeStep"] =
                boost::lexical_cast<std::string>(m_timestep);

        // set boundary conditions here so that any normal component
        // correction are imposed before they are imposed on initial
        // field below
        SetBoundaryConditions(m_time);

        m_F = Array<OneD, Array< OneD, NekDouble> > (m_nConvectiveFields);
        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
#if 0  // disabled since we would need to put in a rotation for non-planar periodic BCs 
            m_fields[i]->LocalToGlobal();
            m_fields[i]->ImposeDirichletConditions(m_fields[i]->UpdateCoeffs());
            m_fields[i]->GlobalToLocal();
#endif
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                  m_fields[i]->UpdatePhys());
            m_F[i] = Array< OneD, NekDouble> (m_fields[0]->GetTotPoints(), 0.0);
        }
    }
    

    /**
     * 
     */
    void VelocityCorrectionScheme:: v_TransCoeffToPhys(void)
    {
        int nfields = m_fields.num_elements() - 1;
        for (int k=0 ; k < nfields; ++k)
        {
            //Backward Transformation in physical space for time evolution
            m_fields[k]->BwdTrans_IterPerExp(m_fields[k]->GetCoeffs(),
                                             m_fields[k]->UpdatePhys());
        }
    }
    
    /**
     * 
     */
    void VelocityCorrectionScheme:: v_TransPhysToCoeff(void)
    {
        
        int nfields = m_fields.num_elements() - 1;
        for (int k=0 ; k < nfields; ++k)
        {
            //Forward Transformation in physical space for time evolution
            m_fields[k]->FwdTrans_IterPerExp(m_fields[k]->GetPhys(),
                                             m_fields[k]->UpdateCoeffs());
        }
    }
    
    /**
     * 
     */
    Array<OneD, bool> VelocityCorrectionScheme::v_GetSystemSingularChecks()
    {
        int vVar = m_session->GetVariables().size();
        Array<OneD, bool> vChecks(vVar, false);
        vChecks[vVar-1] = true;
        return vChecks;
    }
    
    /**
     * 
     */
    int VelocityCorrectionScheme::v_GetForceDimension()
    {
        return m_session->GetVariables().size() - 1;
    }

    /**
     * Explicit part of the method - Advection, Forcing + HOPBCs
     */
    void VelocityCorrectionScheme::v_EvaluateAdvection_SetPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        Array<OneD, Array<OneD, NekDouble> > &outarray, 
        const NekDouble time)
    {
        EvaluateAdvectionTerms(inarray, outarray);

        // Smooth advection
        if(m_SmoothAdvection)
        {
            for(int i = 0; i < m_nConvectiveFields; ++i)
            {
                m_pressure->SmoothField(outarray[i]);
            }
        }

        // Add forcing terms
        for (auto &x : m_forcing)
        {
            x->Apply(m_fields, inarray, outarray, time);
        }

        // Calculate High-Order pressure boundary conditions
        m_extrapolation->EvaluatePressureBCs(inarray,outarray,m_kinvis);
    }
    
    /**
     * Implicit part of the method - Poisson + nConv*Helmholtz
     */
    void VelocityCorrectionScheme::SolveUnsteadyStokesSystem(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        Array<OneD, Array<OneD, NekDouble> > &outarray, 
        const NekDouble time, 
        const NekDouble aii_Dt)
    {
        // Substep the pressure boundary condition if using substepping
        m_extrapolation->SubStepSetPressureBCs(inarray,aii_Dt,m_kinvis);

        // Set up forcing term for pressure Poisson equation
        SetUpPressureForcing(inarray, m_F, aii_Dt);

        // Solve Pressure System
        SolvePressure (m_F[0]);

        // Set up forcing term for Helmholtz problems
        SetUpViscousForcing(inarray, m_F, aii_Dt);

        // Solve velocity system
        SolveViscous( m_F, outarray, aii_Dt);
    }
        
    /**
     * Forcing term for Poisson solver solver
     */ 
    void   VelocityCorrectionScheme::v_SetUpPressureForcing(
        const Array<OneD, const Array<OneD, NekDouble> > &fields, 
        Array<OneD, Array<OneD, NekDouble> > &Forcing, 
        const NekDouble aii_Dt)
    {                
        int i;
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_velocity.num_elements();

        m_fields[0]->PhysDeriv(eX,fields[0], Forcing[0]);
        
        for(i = 1; i < nvel; ++i)
        {
            // Use Forcing[1] as storage since it is not needed for the pressure
            m_fields[i]->PhysDeriv(DirCartesianMap[i],fields[i],Forcing[1]);
            Vmath::Vadd(physTot,Forcing[1],1,Forcing[0],1,Forcing[0],1);
        }

        Vmath::Smul(physTot,1.0/aii_Dt,Forcing[0],1,Forcing[0],1);
    }
    
    /**
     * Forcing term for Helmholtz solver
     */
    void   VelocityCorrectionScheme::v_SetUpViscousForcing(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        Array<OneD, Array<OneD, NekDouble> > &Forcing, 
        const NekDouble aii_Dt)
    {
        NekDouble aii_dtinv = 1.0/aii_Dt;
        int phystot = m_fields[0]->GetTotPoints();

        // Grad p
        m_pressure->BwdTrans(m_pressure->GetCoeffs(),m_pressure->UpdatePhys());

        int nvel = m_velocity.num_elements();
        if(nvel == 2)
        {
            m_pressure->PhysDeriv(m_pressure->GetPhys(),
                                  Forcing[m_velocity[0]],
                                  Forcing[m_velocity[1]]);
        }
        else
        {
            m_pressure->PhysDeriv(m_pressure->GetPhys(),
                                  Forcing[m_velocity[0]], 
                                  Forcing[m_velocity[1]],
                                  Forcing[m_velocity[2]]);
        }

        // zero convective fields. 
        for(int i = nvel; i < m_nConvectiveFields; ++i)
        {
            Vmath::Zero(phystot,Forcing[i],1);
        }
        
        // Subtract inarray/(aii_dt) and divide by kinvis. Kinvis will
        // need to be updated for the convected fields.
        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            Blas::Daxpy(phystot,-aii_dtinv,inarray[i],1,Forcing[i],1);
            Blas::Dscal(phystot,1.0/m_diffCoeff[i],&(Forcing[i])[0],1);
        }
    }

    
    /**
     * Solve pressure system
     */
    void   VelocityCorrectionScheme::v_SolvePressure(
        const Array<OneD, NekDouble>  &Forcing)
    {
        StdRegions::ConstFactorMap factors;
        // Setup coefficient for equation
        factors[StdRegions::eFactorLambda] = 0.0;

        // Solver Pressure Poisson Equation
        m_pressure->HelmSolve(Forcing, m_pressure->UpdateCoeffs(),
                              NullFlagList, factors);

        // Add presure to outflow bc if using convective like BCs
        m_extrapolation->AddPressureToOutflowBCs(m_kinvis);
    }
    
    /**
     * Solve velocity system
     */
    void   VelocityCorrectionScheme::v_SolveViscous(
        const Array<OneD, const Array<OneD, NekDouble> > &Forcing,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble aii_Dt)
    {
        StdRegions::ConstFactorMap factors;
        StdRegions::VarCoeffMap varCoeffMap = StdRegions::NullVarCoeffMap;
        MultiRegions::VarFactorsMap varFactorsMap =
            MultiRegions::NullVarFactorsMap;

        if(m_useSpecVanVisc)
        {
            factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
            factors[StdRegions::eFactorSVVDiffCoeff]   = m_sVVDiffCoeff/m_kinvis;
            if(m_svvVarDiffCoeff != NullNekDouble1DArray)
            {
                if(m_IsSVVPowerKernel)
                {
                    varFactorsMap[StdRegions::eFactorSVVPowerKerDiffCoeff] =
                        m_svvVarDiffCoeff;
                }
                else
                {
                    varFactorsMap[StdRegions::eFactorSVVDGKerDiffCoeff] =
                        m_svvVarDiffCoeff;
                }
            }
        }
        
#if 1
        {
            static GlobalLinSysSharedPtr linSys;
            static bool setup_linsys = true;
            static NekDouble save_aii_Dt = 0;


            if(fabs(save_aii_Dt - aii_Dt) > 1e-8)
            {
                setup_linsys = true; 
            }
            
            if(setup_linsys == true)
            {
                setup_linsys = false;
                save_aii_Dt = aii_Dt;

                m_locToGloMapVec = Array<OneD, AssemblyMapSharedPtr>
                    (m_nConvectiveFields);

                Array<OneD, std::weak_ptr<ExpList> > fields(m_nConvectiveFields);
                
                for(int n  = 0; n < m_nConvectiveFields; ++n)
                {
                    factors[StdRegions::eFactorLambda] = 1.0/aii_Dt/m_diffCoeff[n];
                    m_locToGloMapVec[n] = m_fields[n]->GetLocToGloMap();
                    fields[n] = m_fields[n];
                }

                GlobalLinSysKey key(
                                    StdRegions::eHelmholtz, m_locToGloMapVec[0], factors, varCoeffMap, varFactorsMap);

                // note this assumes m_fields hold u,v,w consecutively!
                linSys = MemoryManager<
                    GlobalLinSysIterativeStaticCondVec>::
                    AllocateSharedPtr(key,fields,m_locToGloMapVec);

                linSys->InitObject(); 
            }

            Array<OneD, Array<OneD, NekDouble> >F(m_nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> >coeffs(m_nConvectiveFields);

            // Inner product of forcing
            int contNcoeffs = m_locToGloMapVec[0]->GetNumGlobalCoeffs();
            int Ncoeffs = m_fields[0]->GetNcoeffs();

            // Solve Helmholtz system and put in Physical space
            for(int n = 0; n < m_nConvectiveFields; ++n)
            {
                F[n] = Array<OneD, NekDouble>(Ncoeffs);
                
                m_fields[n]->IProductWRTBase(Forcing[n],F[n]);
                
                // Note -1.0 term necessary to invert forcing function to
                // be consistent with matrix definition
                Vmath::Neg(Ncoeffs, F[n], 1);
                
                // Forcing function with weak boundary conditions
                int i,j;
                int bndcnt = 0;
                NekDouble sign;
                Array<OneD, NekDouble> gamma(Ncoeffs, 0.0);
                
                Array<OneD, const MultiRegions::ExpListSharedPtr>
                    bndCondExpansions = m_fields[n]->GetBndCondExpansions();
                
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> 
                    bndConditions = m_fields[n]->GetBndConditions();
                
                for(i = 0; i < bndCondExpansions.num_elements(); ++i)
                {
                    if(bndConditions[i]->GetBoundaryConditionType()
                       != SpatialDomains::eDirichlet)
                    {
                        Array<OneD, NekDouble> sign = m_locToGloMapVec[n]->GetBndCondCoeffsToLocalCoeffsSign();
                        const Array<OneD, const int> map= m_locToGloMapVec[n]->GetBndCondCoeffsToLocalCoeffsMap();
                        const Array<OneD, NekDouble> bndcoeff = (bndCondExpansions[i])->GetCoeffs(); 

                        if(m_locToGloMapVec[n]->GetSignChange())
                        {
                            for(j = 0; j < (bndCondExpansions[i])->GetNcoeffs(); j++)
                            {
                                gamma[map[j]] += sign[j] * bndcoeff[j]; 
                            }
                        }
                        else
                        {
                            for(j = 0; j < (bndCondExpansions[i])->GetNcoeffs(); j++)
                            {
                                gamma[map[j]] += bndcoeff[j]; 
                            }
                        }
                    }
                    else
                    {
                        bndcnt += bndCondExpansions[i]->GetNcoeffs();
                    }
                    
                    // Add weak boundary conditions to forcing
                    Vmath::Vadd(Ncoeffs, F[n], 1, gamma, 1, F[n], 1);
                }
            }
                    
            for(int n = 0; n < m_nConvectiveFields; ++n)
            {
                coeffs[n] = m_fields[n]->UpdateCoeffs();
            }
                
            // Solve the system
            linSys->SolveVec(F,coeffs);
            
            for(int n = 0; n < m_nConvectiveFields; ++n)
            {
                m_fields[n]->BwdTrans(m_fields[n]->GetCoeffs(),outarray[n]);
            }
            
        }
#else
        // Solve Helmholtz system and put in Physical space
        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            // Setup coefficients for equation
            factors[StdRegions::eFactorLambda] = 1.0/aii_Dt/m_diffCoeff[i];
            m_fields[i]->HelmSolve(Forcing[i], m_fields[i]->UpdateCoeffs(),
                                   NullFlagList,  factors, varCoeffMap,                                            varFactorsMap);
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),outarray[i]);
        }
#endif
    }
    
    void VelocityCorrectionScheme::SVVVarDiffCoeff(
                     const NekDouble velmag, 
                     Array<OneD, NekDouble> &diffcoeff,
                     const  Array<OneD, Array<OneD, NekDouble> >  &vel)
    {
        int phystot = m_fields[0]->GetTotPoints();
        int nel = m_fields[0]->GetNumElmts();
        int nvel,cnt; 
        
        Array<OneD, NekDouble> tmp;
        
        Vmath::Fill(nel,velmag,diffcoeff,1);
        
        if(vel != NullNekDoubleArrayofArray)
        {
            Array<OneD, NekDouble> Velmag(phystot);
            nvel = vel.num_elements();
            // calculate magnitude of v
            Vmath::Vmul(phystot,vel[0],1,vel[0],1,Velmag,1);
            for(int n = 1; n < nvel; ++n)
            {
                Vmath::Vvtvp(phystot,vel[n],1,vel[n],1,Velmag,1,
                             Velmag,1);
            }
            Vmath::Vsqrt(phystot,Velmag,1,Velmag,1);
                

            cnt = 0;
            Array<OneD, NekDouble> tmp;
            // calculate mean value of vel mag. 
            for(int i = 0; i < nel; ++i)
            {
                int nq = m_fields[0]->GetExp(i)->GetTotPoints();
                tmp = Velmag + cnt;
                diffcoeff[i] = m_fields[0]->GetExp(i)->Integral(tmp);
                Vmath::Fill(nq,1.0,tmp,1);
                NekDouble area = m_fields[0]->GetExp(i)->Integral(tmp);
                diffcoeff[i] = diffcoeff[i]/area;
                cnt += nq;
            }
        }
        else
        {
            nvel = m_expdim;
        }
        
        if(m_expdim == 3)
        {
            LocalRegions::Expansion3DSharedPtr exp3D;
            for (int e = 0; e < nel; e++)
            {
                exp3D = m_fields[0]->GetExp(e)->as<LocalRegions::Expansion3D>();
                NekDouble h = 0; 
                for(int i = 0; i < exp3D->GetNedges(); ++i)
                {
                    
                    h = max(h, exp3D->GetGeom3D()->GetEdge(i)->GetVertex(0)->dist(
                             *(exp3D->GetGeom3D()->GetEdge(i)->GetVertex(1))));
                }

                NekDouble p;
                for(int i = 0; i < 3; ++i)
                {
                    p = max(p,exp3D->GetBasisNumModes(i)-1.0);
                }
                
                diffcoeff[e] *= h/p; 
            }
        }
        else
        {
            LocalRegions::Expansion2DSharedPtr exp2D;
            for (int e = 0; e < nel; e++)
            {
                exp2D = m_fields[0]->GetExp(e)->as<LocalRegions::Expansion2D>();
                NekDouble h = 0;
                for(int i = 0; i < exp2D->GetNedges(); ++i)
                {
                    
                   h = max(h, exp2D->GetGeom2D()->GetEdge(i)->GetVertex(0)->dist(
                             *(exp2D->GetGeom2D()->GetEdge(i)->GetVertex(1))));
                }

                NekDouble p;
                for(int i = 0; i < 2; ++i)
                {
                    p = max(p,exp2D->GetBasisNumModes(i)-1.0);
                }
                
                diffcoeff[e] *= h/p; 
            }
        }
    }
} //end of namespace
