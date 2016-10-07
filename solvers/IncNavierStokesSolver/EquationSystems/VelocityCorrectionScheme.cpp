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
#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/Core/Misc.h>
#include <boost/algorithm/string.hpp>

using namespace std;

namespace Nektar
{
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
            std::string vExtrapolation = "Standard";

            if (m_session->DefinesSolverInfo("Extrapolation"))
            {
                vExtrapolation = m_session->GetSolverInfo("Extrapolation");
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
        
        m_saved_aii_Dt = Array<OneD, NekDouble>(m_nConvectiveFields,
                                                NekConstants::kNekUnsetDouble);

        // Load parameters for Spectral Vanishing Viscosity
        m_session->MatchSolverInfo("SpectralVanishingViscosity","True",
                                   m_useSpecVanVisc, false);
        m_useHomo1DSpecVanVisc = m_useSpecVanVisc;
        if(m_useSpecVanVisc == false)
        {
            m_session->MatchSolverInfo("SpectralVanishingViscositySpectralHP",
                                "True", m_useSpecVanVisc, false);
            m_session->MatchSolverInfo("SpectralVanishingViscosityHomo1D",
                                "True", m_useHomo1DSpecVanVisc, false);
        }
        m_session->LoadParameter("SVVCutoffRatio",m_sVVCutoffRatio,0.75);
        m_session->LoadParameter("SVVDiffCoeff",  m_sVVDiffCoeff,  0.1);

        m_session->MatchSolverInfo("SPECTRALHPDEALIASING","True",
                                   m_specHP_dealiasing,false);

        if(m_HomogeneousType == eHomogeneous1D)
        {
            ASSERTL0(m_nConvectiveFields > 2,"Expect to have three velocity fields with homogenous expansion");

            if(m_useHomo1DSpecVanVisc)
            {
                Array<OneD, unsigned int> planes;
                planes = m_fields[0]->GetZIDs();

                int num_planes = planes.num_elements();
                Array<OneD, NekDouble> SVV(num_planes,0.0);
                NekDouble fac;
                int kmodes = m_fields[0]->GetHomogeneousBasis()->GetNumModes();
                int pstart;

                pstart = m_sVVCutoffRatio*kmodes;
                
                for(n = 0; n < num_planes; ++n)
                {
                    if(planes[n] > pstart)
                    {
                        fac = (NekDouble)((planes[n] - kmodes)*(planes[n] - kmodes))/
                            ((NekDouble)((planes[n] - pstart)*(planes[n] - pstart)));
                        SVV[n] = m_sVVDiffCoeff*exp(-fac)/m_kinvis;
                    }
                }

                for(int i = 0; i < m_velocity.num_elements(); ++i)
                {
                    m_fields[m_velocity[i]]->SetHomo1DSpecVanVisc(SVV);
                }
            }
            
        }

        m_session->MatchSolverInfo("SmoothAdvection", "True", m_SmoothAdvection, false);

        if(m_session->DefinesSolverInfo("DynamicViscosity"))
        {

            m_dynamicVisc = MemoryManager<DynamicViscData>::AllocateSharedPtr();
            m_dynamicVisc->m_type = m_session->GetSolverInfo("DynamicViscosity");

            m_dynamicVisc->m_sensorField  = Array<OneD, NekDouble>
                (m_fields[0]->GetTotPoints());

            m_dynamicVisc->m_numSteps     = 0;

            if(boost::iequals(m_dynamicVisc->m_type,"SemiImplicitVariableDiff"))
            {
                ForcingDynamicViscSharedPtr forcing = MemoryManager<ForcingDynamicVisc>::AllocateSharedPtr(m_session);
                forcing->InitObject(m_fields, m_velocity.num_elements(), (TiXmlElement*)NULL);
                // initialisae kinvis as zero
                Array<OneD, NekDouble> Zero(m_fields[0]->GetTotPoints(),0.0);
                forcing->SetKinvis(Zero);
                m_dynamicVisc->m_forcing = forcing;
                m_forcing.push_back(forcing);

                m_session->LoadParameter("DynamicViscEvalSteps",
                                         m_dynamicVisc->m_numStepsAvg,5);
                m_session->LoadParameter("DynamicViscC0",
                                         m_dynamicVisc->m_kinvisC0,0.1);
                m_session->LoadParameter("DynamicViscDefS0",
                                         m_dynamicVisc->m_defS0,3);

                // store original kinvis for later use
                m_dynamicVisc->m_origKinvis = m_kinvis;

                // initialise with zero
                m_dynamicVisc->m_fixedKinvis = 0.0;
                
            }
            else
            {
                m_session->LoadParameter("DynamicViscEvalSteps",
                                         m_dynamicVisc->m_numStepsAvg,100);
            }

        }

        // set explicit time-intregration class operators
        m_ode.DefineOdeRhs(&VelocityCorrectionScheme::EvaluateAdvection_SetPressureBCs, this);

        m_extrapolation->SubSteppingTimeIntegration(m_intScheme->GetIntegrationMethod(), m_intScheme);
        m_extrapolation->GenerateHOPBCMap();
        
        // set implicit time-intregration class operators
        m_ode.DefineImplicitSolve(&VelocityCorrectionScheme::SolveUnsteadyStokesSystem,this);
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
        UnsteadySystem::v_GenerateSummary(s);

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
        if (m_useHomo1DSpecVanVisc && (m_HomogeneousType == eHomogeneous1D))
        {
            smoothing += (smoothing == "" ? "" : " + ") + string("Homogeneous1D");
        }
        if (smoothing != "")
        {
            SolverUtils::AddSummaryItem(
                s, "Smoothing", "SVV (" + smoothing + " SVV (cut-off = "
                + boost::lexical_cast<string>(m_sVVCutoffRatio)
                + ", diff coeff = "
                + boost::lexical_cast<string>(m_sVVDiffCoeff)+")");
        }
    }

    /**
     * 
     */
    void VelocityCorrectionScheme::v_DoInitialise(void)
    {

        UnsteadySystem::v_DoInitialise();

        // Set up Field Meta Data for output files
        m_fieldMetaDataMap["Kinvis"]   = boost::lexical_cast<std::string>(m_kinvis);
        m_fieldMetaDataMap["TimeStep"] = boost::lexical_cast<std::string>(m_timestep);

        // set boundary conditions here so that any normal component
        // correction are imposed before they are imposed on initial
        // field below
        SetBoundaryConditions(m_time);

        m_F = Array<OneD, Array< OneD, NekDouble> > (m_nConvectiveFields);
        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            m_fields[i]->LocalToGlobal();
            m_fields[i]->ImposeDirichletConditions(m_fields[i]->UpdateCoeffs());
            m_fields[i]->GlobalToLocal();
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
            m_fields[k]->FwdTrans_IterPerExp(m_fields[k]->GetPhys(),m_fields[k]->UpdateCoeffs());
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
        std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
        for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
        {
            (*x)->Apply(m_fields, inarray, outarray, time);
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
        // Enforcing boundary conditions on all fields
        SetBoundaryConditions(time);

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
     * Perform the extrapolation.
     */
    bool   VelocityCorrectionScheme::v_PreIntegrate(int step)
    {
        if(m_dynamicVisc)
        {
            DynamicVisc();
        }
        
        return false;
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

        m_fields[0]->PhysDeriv(0,fields[0],
                               Forcing[0]);
        for(i = 1; i < nvel; ++i)
        {
            // Use Forcing[1] as storage since it is not needed for the pressure
            m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[i],fields[i],Forcing[1]);
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
            m_pressure->PhysDeriv(m_pressure->GetPhys(), Forcing[0], Forcing[1]);
        }
        else
        {
            m_pressure->PhysDeriv(m_pressure->GetPhys(), Forcing[0], 
                                  Forcing[1], Forcing[2]);
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
        m_pressure->HelmSolve(Forcing, m_pressure->UpdateCoeffs(), NullFlagList,
                              factors);
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

        if(m_useSpecVanVisc)
        {
            factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
            factors[StdRegions::eFactorSVVDiffCoeff]   = m_sVVDiffCoeff/m_kinvis;            
        }

        if(m_dynamicVisc)
        {
            varCoeffMap = m_dynamicVisc->m_savVarCoeffMap;
        }

        // Solve Helmholtz system and put in Physical space
        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            // Setup coefficients for equation
            factors[StdRegions::eFactorLambda] = 1.0/aii_Dt/m_diffCoeff[i];
            m_fields[i]->HelmSolve(Forcing[i], m_fields[i]->UpdateCoeffs(),
                                   NullFlagList,  factors, varCoeffMap);
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),outarray[i]);
        }
    }

    void VelocityCorrectionScheme::v_ExtraFldOutput(
            std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
            std::vector<std::string>             &variables)
    {
        if(m_dynamicVisc&&m_steps)
        {
            if(m_dynamicVisc->m_savVarCoeffMap.count(StdRegions::eVarCoeffD00) != 0)
            {
                int ncoeffs = m_fields[0]->GetNcoeffs();

                variables.push_back("DynamicViscosity");
                Array<OneD, NekDouble> newfld(ncoeffs); 
                fieldcoeffs.push_back(newfld);  

                m_fields[0]->FwdTrans_IterPerExp(
                       m_dynamicVisc->m_savVarCoeffMap[StdRegions::eVarCoeffD00],
                       fieldcoeffs[fieldcoeffs.size()-1]);
                // rescale with kinvis since we dividing by kinvis earlier
                Vmath::Smul(ncoeffs,m_kinvis,
                            fieldcoeffs[fieldcoeffs.size()-1],1,
                            fieldcoeffs[fieldcoeffs.size()-1],1);
            }
            
            else if(boost::iequals(m_dynamicVisc->m_type,
                                   "SemiImplicitVariableDiff"))
            {
                int ncoeffs = m_fields[0]->GetNcoeffs();
                int nquad   = m_fields[0]->GetTotPoints();
                
                variables.push_back("DynamicViscosity");
                Array<OneD, NekDouble> newfld(ncoeffs);
                Array<OneD, NekDouble> kinvis(nquad);
                fieldcoeffs.push_back(newfld);  

                Vmath::Sadd(nquad,m_kinvis,m_dynamicVisc->m_forcing->GetKinvis(),
                            1,kinvis,1);

                m_fields[0]->FwdTrans_IterPerExp(kinvis,
                       fieldcoeffs[fieldcoeffs.size()-1]);
            }
            
        }
    }
    
    void VelocityCorrectionScheme::DynamicVisc()
    {
        // Update sensor field with u^2
        int nquad = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> energy(nquad,0.0);
        Vmath::Vmul(nquad,m_fields[m_velocity[0]]->GetPhys(),1,
                    m_fields[m_velocity[0]]->GetPhys(),1,energy,1);   
        
        for(int i = 1; i < m_velocity.num_elements(); ++i)
        {
            Vmath::Vvtvp(nquad,m_fields[m_velocity[i]]->GetPhys(),1,
                         m_fields[m_velocity[i]]->GetPhys(),1,energy,1,
                         energy,1); 
        }

        Vmath::Vsqrt(nquad,energy,1,energy,1);
        Vmath::Vadd(nquad,energy,1,m_dynamicVisc->m_sensorField,1,
                    m_dynamicVisc->m_sensorField,1);
        
            
        // Calculate new viscosity if required
        if(boost::iequals(m_dynamicVisc->m_type,"VariableDiff"))
        {
            if(m_dynamicVisc->m_numSteps%m_dynamicVisc->m_numStepsAvg == 0)
            {
                
                if (m_comm->GetRank() == 0)
                {
                    cout << "Reseting up variable viscosity in step:" << m_dynamicVisc->m_numSteps << endl;
                }
                
                
                Vmath::Smul(nquad,1.0/(m_dynamicVisc->m_numSteps +1.0),
                            m_dynamicVisc->m_sensorField,1,energy,1);
                
                Array<OneD, NekDouble> kinvis(nquad);
                
                GetStabiliseKinvis(energy,m_dynamicVisc->m_kinvisC0,
                                   m_dynamicVisc->m_defS0,kinvis);
                
                // Release  current linear sys
                for(int i =0 ; i < m_velocity.num_elements(); ++i)
                {
                    m_fields[m_velocity[i]]->ClearGlobalLinSysManager();
                    m_fields[m_velocity[i]]->ResetLocalManagers();
                }
                
                // add in dynamics viscosity and then divide by
                // dyanmics viscosity for consistency with standard
                // implementation.
                Vmath::Sadd(nquad,m_kinvis,kinvis,1,kinvis,1);
                Vmath::Smul(nquad,1.0/m_kinvis,kinvis,1,kinvis,1);
                m_dynamicVisc->m_savVarCoeffMap[StdRegions::eVarCoeffD00]
                    = kinvis; 
                m_dynamicVisc->m_savVarCoeffMap[StdRegions::eVarCoeffD11]
                    = kinvis;
                if(m_expdim == 3)
                {
                    m_dynamicVisc->m_savVarCoeffMap[StdRegions::eVarCoeffD22]
                        = kinvis;
                }
            }
            m_dynamicVisc->m_numSteps = 0; 
            // reset average field wiht current average
            Vmath::Vcopy(nquad,energy,1,m_dynamicVisc->m_sensorField,1);
            
            m_dynamicVisc->m_numSteps++;
        }
        else if(boost::iequals(m_dynamicVisc->m_type,
                               "SemiImplicitVariableDiff"))
        {
            
            Vmath::Smul(nquad,1.0/(m_dynamicVisc->m_numSteps +1.0),
                        m_dynamicVisc->m_sensorField,1,energy,1);
            
            Array<OneD, NekDouble> kinvis(nquad);
            
            GetStabiliseKinvis(energy,m_dynamicVisc->m_kinvisC0,
                               m_dynamicVisc->m_defS0,kinvis);

            // add in dynamics viscosity and then divide by
            // dyanmics viscosity for consistency with standard
            // implementation.
            Vmath::Sadd(nquad,m_dynamicVisc->m_origKinvis,kinvis,1,kinvis,1);
            
            NekDouble maxkinvis = Vmath::Vmax(nquad,kinvis,1);

            // reset max kinvis to 60% of current and 40% of previous
            // to try and damp large variations. 
            if(m_dynamicVisc->m_numSteps)
            {
                maxkinvis = 0.3*maxkinvis +  0.7*m_dynamicVisc->m_fixedKinvis;
            }
            
            if((maxkinvis > 2*m_dynamicVisc->m_fixedKinvis)||
               (maxkinvis < 0.5*m_dynamicVisc->m_fixedKinvis))
            {
                m_dynamicVisc->m_fixedKinvis = maxkinvis;

                if (m_comm->GetRank() == 0)
                {
                    cout << "Updating fixed viscosity: " <<
                        m_dynamicVisc->m_numSteps << " (value: "
                         << maxkinvis << ")" << endl;
                }
                
                
                // Release  current linear sys
                for(int i =0 ; i < m_velocity.num_elements(); ++i)
                {
                    m_fields[m_velocity[i]]->ClearGlobalLinSysManager();
                    m_fields[m_velocity[i]]->ResetLocalManagers();
                }
            
                for(int i =0 ; i < m_velocity.num_elements(); ++i)
                {
                    m_diffCoeff[m_velocity[i]] = maxkinvis;
                }
                m_kinvis = maxkinvis;
            }
            
            // calculate variable viscosity 
            Vmath::Sadd(nquad,-1*m_dynamicVisc->m_fixedKinvis,
                        kinvis,1,kinvis,1);
            m_dynamicVisc->m_forcing->SetKinvis(kinvis);


            if(m_dynamicVisc->m_numSteps%m_dynamicVisc->m_numStepsAvg == 0)
            {
                m_dynamicVisc->m_numSteps = 0; 
                // reset average field wiht current average
                Vmath::Vcopy(nquad,energy,1,m_dynamicVisc->m_sensorField,1);
            }

            m_dynamicVisc->m_numSteps++;
        }
        else
        {
            ASSERTL0(false,"Unknown definition for SolverInfo DynamicsViscosity");
        }
    }
    void VelocityCorrectionScheme::GetStabiliseKinvis(
                               const Array<OneD, const NekDouble > &physarray,
                               NekDouble C,
                               NekDouble S0_def,
                               Array<OneD, NekDouble>   &StabKinvis)
    {
        int e, NumModesElement, nQuadPointsElement;
        int nTotQuadPoints  = GetTotPoints();
        int nElements       = m_fields[0]->GetExpSize();

        // Find solution (SolP) at p = P;
        // The input array (physarray) is the solution at p = P;
        Array<OneD,int> ExpOrderElement = GetNumExpModesPerExp();
        Array<OneD, NekDouble> tmp;

        int nummodes,PhysCount = 0;
        NekDouble e0,h,sensorVal;
        NekDouble SolPmeanNumerator, SolPmeanDenumerator; 

        NekDouble kappa = 0.5;
        NekDouble s0, kinvis;
        
        int nCoeffsElement, numCutOff;        

        for (e = 0; e < nElements; e++)
        {
            PhysCount = m_fields[0]->GetPhys_Offset(e);
            NumModesElement    = ExpOrderElement[e];
            nQuadPointsElement = m_fields[0]->GetExp(e)->GetTotPoints();
            nCoeffsElement     = m_fields[0]->GetExp(e)->GetNcoeffs();
            numCutOff          = NumModesElement - 1;

            // Set-up of the Orthogonal basis for an element which
            // is needed to obtain thesolution at P =  p - 1;

            Array<OneD, NekDouble> SolPElementPhys  (nQuadPointsElement, 0.0);
            Array<OneD, NekDouble> SolPElementCoeffs(nCoeffsElement,     0.0);
            Array<OneD, NekDouble> SolNorm          (nQuadPointsElement, 0.0);

            Array<OneD, NekDouble> SolPmOneElementPhys(nQuadPointsElement, 0.0);
            Array<OneD, NekDouble> SolPmOneElementCoeffs(nCoeffsElement, 0.0);

            // create vector the save the solution points per element at P = p;
            Vmath::Vcopy(nQuadPointsElement,physarray+PhysCount,1,SolPElementPhys,1);

            m_fields[0]->GetExp(e)->FwdTrans(SolPElementPhys, SolPElementCoeffs);

            // ReduceOrderCoeffs reduces the polynomial order of the solution
            // that is represented by the coeffs given as an inarray. This is
            // done by projecting the higher order solution onto the orthogonal
            // basis and padding the higher order coefficients with zeros.
            m_fields[0]->GetExp(e)->ReduceOrderCoeffs(numCutOff,
                                                      SolPElementCoeffs,
                                                      SolPmOneElementCoeffs);

            m_fields[0]->GetExp(e)->BwdTrans(SolPmOneElementCoeffs,
                                             SolPmOneElementPhys);

            // Determining the norm of the numerator of the Sensor
            Vmath::Vsub(nQuadPointsElement, SolPElementPhys, 1,
                        SolPmOneElementPhys, 1, SolNorm, 1);

            Vmath::Vmul(nQuadPointsElement, SolNorm, 1, SolNorm, 1,
                        SolNorm, 1);

            Vmath::Vmul(nQuadPointsElement, SolPElementPhys, 1,
                        SolPElementPhys, 1, SolPElementPhys, 1);
            
            SolPmeanNumerator   = m_fields[0]->GetExp(e)->Integral(SolNorm);
            SolPmeanDenumerator = m_fields[0]->GetExp(e)->Integral(SolPElementPhys);

            sensorVal = SolPmeanNumerator / SolPmeanDenumerator;
            sensorVal = log10(sensorVal);
            
            nummodes = m_fields[0]->GetExp(e)->GetBasisNumModes(0)-1;
            
            Array<OneD, NekDouble> One(nQuadPointsElement,1.0);
            h = m_fields[0]->GetExp(e)->Integral(One);
            h = pow(h,1.0/(NekDouble) m_spacedim);
            e0 = C*h/nummodes;
            s0 = -(S0_def + 4.0*log10(nummodes));
            
            if(sensorVal < s0 - kappa)
            {
                kinvis = 0.0;
            }
            else if (sensorVal > s0 + kappa)
            {
                kinvis = e0; 
            }
            else
            {
                kinvis = 0.5*e0*(1 + sin(M_PI*(sensorVal-s0)/(2*kappa)));
            }
            
            Vmath::Fill(nQuadPointsElement,kinvis,
                        tmp = StabKinvis + PhysCount,1);
        }
    }

} //end of namespace
