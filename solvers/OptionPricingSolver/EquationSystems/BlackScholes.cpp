///////////////////////////////////////////////////////////////////////////////
//
// File BlackScholes.cpp
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
// Description: Black Scholes solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <OptionPricingSolver/EquationSystems/BlackScholes.h>

using namespace std;

namespace Nektar
{

    string BlackScholes::className = SolverUtils::GetEquationSystemFactory().
        RegisterCreatorFunction("BlackScholes", BlackScholes::create,
                                "Black Scholes equation.");
    
    BlackScholes::BlackScholes(const LibUtilities::SessionReaderSharedPtr& pSession)
            : UnsteadySystem(pSession), AdvectionSystem(pSession)
    {
    }
    
    /**
     * @brief Initialisation object for the unsteady linear advection 
     * diffusion equation.
     */
    void BlackScholes::v_InitObject()
    {
        AdvectionSystem::v_InitObject();
        
        // Risk-free interest rate and volatility
        m_session->LoadParameter("Interest"  , m_interest  , 0.0);
        m_session->LoadParameter("Volatility", m_volatility, 0.0);
        
        // turn on substepping
        m_session->MatchSolverInfo("Extrapolation", "SubStepping",
                                   m_subSteppingScheme, false);
        
        // Define Velocity fields
        const int nq = m_fields[0]->GetNpoints();
        m_stockPrice = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
#if 0 
        std::vector<std::string> stock;
        stock.push_back("Vx");
        stock.push_back("Vy");
        stock.push_back("Vz");
        stock.resize(m_spacedim);
        EvaluateFunction(stock, m_stockPrice, "AdvectionFactor");
#else
        // Load local coords into stockprice
        Array<OneD, Array<OneD, NekDouble> > coords(m_spacedim);
        m_stockPrice[0] = Array<OneD, NekDouble> (nq);
        coords[0]       = Array<OneD, NekDouble> (nq);
        
        switch (m_spacedim)
        {
            case 1:
                m_fields[0]->GetCoords(coords[0]);
                break;
            case 2:
                m_stockPrice[1] = Array<OneD, NekDouble> (nq);
                coords[1]       = Array<OneD, NekDouble> (nq);
                m_fields[0]->GetCoords(coords[0], coords[1]);
                break;
            case 3:
                m_stockPrice[2] = Array<OneD, NekDouble> (nq);
                coords[2]       = Array<OneD, NekDouble> (nq);
                m_fields[0]->GetCoords(coords[0], coords[1], coords[3]);
            default:
                break;
        }
#endif
        
        // Multiply the underlying stock price at time \tau = T - t
        // by the free-risk interest-rate, that is: r * S
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Smul(m_stockPrice[0].num_elements(), m_interest,
                        coords[i], 1, m_stockPrice[i], 1);
        }
        
        // Define diffusion variable containers
        StdRegions::VarCoeffType varCoeffEnum[3] = {
                StdRegions::eVarCoeffD00,
                StdRegions::eVarCoeffD01,
                StdRegions::eVarCoeffD11
        };

        std::string varCoeffString[3] = {"xx","xy","yy"};

        const int nVarDiffCmpts = m_spacedim * (m_spacedim + 1) / 2;

        // Allocate storage for variable coeffs and initialize to 1
        for (int i = 0, k = 0; i < m_spacedim; ++i)
        {
            for (int j = 0; j < i+1; ++j)
            {
                if (i == j)
                {
                    m_vardiff[varCoeffEnum[k]] = Array<OneD,NekDouble>(nq, 1.0);
                }
                else
                {
                    m_vardiff[varCoeffEnum[k]] = Array<OneD,NekDouble>(nq, 0.0);
                }
                ++k;
            }
        }

        // Forcing terms
        m_forcing = SolverUtils::Forcing::Load(m_session, m_fields,
                                               m_fields.num_elements());

        switch (m_spacedim)
        {
            case 1:
                // Set up variable diffusivity as 0.5 * volatility^2 * S
                Vmath::Smul(m_stockPrice[0].num_elements(),
                            0.5*m_volatility*m_volatility,
                            coords[0],1,m_vardiff[varCoeffEnum[0]],1);
                break;
            case 2:
            default:
                ASSERTL0(false,"Set up variable diffusion");
                break;
        }

        
        m_session->MatchSolverInfo(
            "SpectralVanishingViscosity", "True", m_useSpecVanVisc, false);
        
        if (m_useSpecVanVisc)
        {
            m_session->LoadParameter("SVVCutoffRatio", m_sVVCutoffRatio, 0.75);
            m_session->LoadParameter("SVVDiffCoeff"  , m_sVVDiffCoeff  , 0.10);
        }        

        // Type of advection and diffusion classes to be used
        switch(m_projectionType)
        {
            // Discontinuous SEM
            case MultiRegions::eDiscontinuous:
            {
                // Numerical flux and explicit advection term definitions
                string advName, riemName;
                
                // Set numerical flux approach
                m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
                m_riemannSolver = SolverUtils::GetRiemannSolverFactory().
                    CreateInstance(riemName);
                
                // Set additional parameters for numerical flux
                m_riemannSolver->SetScalar("Vn", &BlackScholes::GetNormalStockPrice, this);
                
                // Set explicit advection strategy
                m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
                m_advObject = SolverUtils::GetAdvectionFactory().
                    CreateInstance(advName, advName);
                
                // Set additional parameters for explicit advection
                m_advObject->SetFluxVector   (&BlackScholes::GetFluxVectorAdv, this);
                m_advObject->SetRiemannSolver(m_riemannSolver);
                m_advObject->InitObject      (m_session, m_fields);
                
                // Check if explicit diffusion strategy is ON
                if (m_explicitDiffusion)
                {
                    // Explicit advection term definition
                    std::string diffName;
                    
                    // Set explicit diffusion strategy
                    m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
                    m_diffusion = SolverUtils::GetDiffusionFactory().
                        CreateInstance(diffName, diffName);
                    
                    // Set additional parameters for explicit diffusion
                    m_diffusion->SetFluxVector(&BlackScholes::GetFluxVectorDiff, this);
                    m_diffusion->InitObject(m_session, m_fields);
                }

                // Sub-stepping for DG projection not implemented
                ASSERTL0(m_subSteppingScheme == false,
                         "SubSteppingScheme is not set up for DG projection");
                break;
            }
            // Continuous SEM
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                // Advection term
                std::string advName;
                m_session->LoadSolverInfo("AdvectionType",
                                          advName, "NonConservative");
                m_advObject = SolverUtils::GetAdvectionFactory().
                    CreateInstance(advName, advName);
                m_advObject->SetFluxVector(&BlackScholes::GetFluxVectorAdv, this);
                
                // Add sigma^2/2 to advection term to handle symmetric
                // diffusion oeprator
                for (int i = 0; i < m_spacedim; ++i)
                {
                    Vmath::Svtvp(m_stockPrice[0].num_elements(),
                                -1.0*m_volatility*m_volatility,
                                coords[i],1,m_stockPrice[i], 1, m_stockPrice[i], 1);
                }


                if (advName.compare("WeakDG") == 0)
                {
                    string riemName;
                    m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
                    m_riemannSolver = SolverUtils::GetRiemannSolverFactory().
                        CreateInstance(riemName);
                    m_riemannSolver->SetScalar(
                        "Vn", &BlackScholes::GetNormalStockPrice, this);
                    m_advObject->SetRiemannSolver(m_riemannSolver);
                    m_advObject->InitObject      (m_session, m_fields);
                }

                // In case of Galerkin explicit diffusion gives an error
                if (m_explicitDiffusion)
                {
                    ASSERTL0(false, "Explicit Galerkin diffusion not set up.");
                }
                // In case of Galerkin implicit diffusion: do nothing
                break;
            }
            default:
            {
                ASSERTL0(false, "Unsupported projection type.");
                break;
            }
        }

        m_ode.DefineImplicitSolve(&BlackScholes::DoImplicitSolve, this);
        m_ode.DefineProjection   (&BlackScholes::DoOdeProjection, this);
        m_ode.DefineOdeRhs       (&BlackScholes::DoOdeRhs       , this);

        // Substepping
        if (m_subSteppingScheme)
        {
            ASSERTL0(m_projectionType == MultiRegions::eMixed_CG_Discontinuous,
            "Projection must be set to Mixed_CG_Discontinuous for substepping");
            SetUpSubSteppingTimeIntegration(m_intScheme->GetIntegrationMethod(),
                                            m_intScheme);
        }
    }

    /**
     * @brief Unsteady linear advection diffusion equation destructor.
     */
    BlackScholes::~BlackScholes()
    {
    }
    
    /**
     * @brief Get the normal stock price for the B-S equation.
     */
    Array<OneD, NekDouble> &BlackScholes::GetNormalStockPrice()
    {
        return GetNormalStock(m_stockPrice);
    }


    /**
     * @brief Get the normal stock price for the B-S equation.
     */
    Array<OneD, NekDouble> &BlackScholes::GetNormalStock(
        const Array<OneD, const Array<OneD, NekDouble> > &stockField)
    {
        // Auxiliary variable to compute the normal velocity
        Array<OneD, NekDouble> tmp(GetTraceNpoints());
        m_traceSn = Array<OneD, NekDouble>(GetTraceNpoints(), 0.0);

        // Reset the normal velocity
        Vmath::Zero(GetTraceNpoints(), m_traceSn, 1);

        for (int i = 0; i < stockField.num_elements(); ++i)
        {
            m_fields[0]->ExtractTracePhys(stockField[i], tmp);
            Vmath::Vvtvp(GetTraceNpoints(), m_traceNormals[i], 1,
                         tmp, 1, m_traceSn, 1, m_traceSn, 1);
        }
        
        return m_traceSn;
    }
    
    /**
     * @brief Compute the right-hand side for the B-S problem.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void BlackScholes::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> > &inarray,
              Array<OneD,        Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        // Calculate explicit advection term a(x) * (\nabla u)
        m_advObject->Advect(inarray.num_elements(), m_fields,
                            m_stockPrice, inarray, outarray, time);
        
        // Add forcing terms
        std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
        for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
        {
            (*x)->Apply(m_fields, inarray, outarray, time);
        }

        if (m_explicitDiffusion)
        {
            Array<OneD, Array<OneD, NekDouble> > outD(inarray.num_elements());
            for (int i = 0; i < inarray.num_elements(); ++i)
            {
                outD[i] = Array<OneD, NekDouble>(GetNpoints(), 0.0);
            }

            // Calculate explicit diffusion term eps(x) * (\Delta u)
            m_diffusion->Diffuse(inarray.num_elements(),
                                 m_fields, inarray, outD);

            for (int i = 0; i < inarray.num_elements(); ++i)
            {
                Vmath::Vadd(GetNpoints(),   &outarray[i][0], 1,
                            &outD[i][0], 1, &outarray[i][0], 1);
            }
        }
    }
    
    /**
     * @brief Compute the projection for the unsteady advection 
     * diffusion problem.
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void BlackScholes::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                  time)
    {
        SetBoundaryConditions(time);
        switch(m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                // Just copy over array
                int npoints = GetNpoints();

                for(int i = 0; i < inarray.num_elements(); ++i)
                {
                    Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
                }
                break;
            }
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

                for(int i = 0; i < inarray.num_elements(); ++i)
                {
                    m_fields[i]->FwdTrans(inarray[i], coeffs);
                    m_fields[i]->BwdTrans_IterPerExp(coeffs, outarray[i]);
                }
                break;
            }
            default:
            {
                ASSERTL0(false, "Unknown projection scheme");
                break;
            }
        }
    }
    
    /* @brief Compute the diffusion term implicitly. 
     * 
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     * @param lambda     Diffusion coefficient.
     */
    void BlackScholes::DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD,       Array<OneD, NekDouble> >&outarray,
        const NekDouble time,
        const NekDouble aii_Dt)
    {
        int nvariables = inarray.num_elements();
        int nq          = m_fields[0]->GetNpoints();
        StdRegions::ConstFactorMap factors;

        // time integration terms and source term
        factors[StdRegions::eFactorLambda] = 1.0/aii_Dt + m_interest;
        
        if(m_useSpecVanVisc)
        {
            factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
            factors[StdRegions::eFactorSVVDiffCoeff]   = m_sVVDiffCoeff; 
        }

        if(m_projectionType == MultiRegions::eDiscontinuous)
        {
            factors[StdRegions::eFactorTau] = 1.0;
        }

        //Setting boundary conditions
        SetBoundaryConditions(time);
        
        for (int i = 0; i < nvariables; ++i)
        {
            // input is aii_Dt * F so need to divide by aii_Dt
            Vmath::Smul(nq, -1.0/aii_Dt, inarray[i], 1,
                        m_fields[i]->UpdatePhys(), 1);
            
            // Solve a system of equations with Helmholtz solver
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                   m_fields[i]->UpdateCoeffs(), 
                                   NullFlagList, factors, m_vardiff);
            
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
        }
    }
    
    /**
     * @brief Return the flux vector for the advection part.
     * 
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void BlackScholes::GetFluxVectorAdv(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        ASSERTL1(flux[0].num_elements() == m_stockPrice.num_elements(),
                 "Dimension of flux array and velocity array do not match");

        for (int i = 0; i < flux.num_elements(); ++i)
        {
            for (int j = 0; j < flux[0].num_elements(); ++j)
            {
                Vmath::Vmul(m_fields[0]->GetNpoints(), physfield[i], 1,
                            m_stockPrice[j], 1, flux[i][j], 1);
            }
        }
    }

    /**
     * @brief Return the flux vector for the diffusion part.
     *      
     * @param i           Equation number.
     * @param j           Spatial direction.
     * @param physfield   Fields.
     * @param derivatives First order derivatives.
     * @param flux        Resulting flux.
     */
    void BlackScholes::GetFluxVectorDiff(
        const int i,
        const int j,
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> > &derivatives,
              Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        for (int k = 0; k < flux.num_elements(); ++k)
        {
            Vmath::Zero(GetNpoints(), flux[k], 1);
        }
        Vmath::Vcopy(GetNpoints(), physfield[i], 1, flux[j], 1);
    }
    
    /**
     * Generate summary with information regarding simulation.
     */
    void BlackScholes::v_GenerateSummary( SolverUtils::SummaryList& s)
    {
        AdvectionSystem::v_GenerateSummary(s);

    }

        
    // =========================================================================
    // BELOW ALL STUFF RELATED TO  SUBSTEPPING - NOT NEEDED AT THE MOMENT
    // =========================================================================
    bool BlackScholes::v_PreIntegrate(int step)
    {
        if(m_subSteppingScheme)
        {
            SubStepAdvance(m_intSoln,step,m_time);
        }

        return false;
    }

    /**
     * Advance substepping time-integration.
     */
    void BlackScholes::SubStepAdvance(
        const LibUtilities::TimeIntegrationSolutionSharedPtr &integrationSoln,
        int                                                  nstep,
        NekDouble                                            time)
    {
        int nsubsteps;
        NekDouble dt; 
        
        Array<OneD, Array<OneD, NekDouble> > fields, velfields;
        
        static int ncalls = 1;
        int  nint         = min(ncalls++, m_intSteps);
        
        Array<OneD, NekDouble> CFL(m_fields[0]->GetExpSize(), 
                                   m_cflSafetyFactor);
        
        LibUtilities::CommSharedPtr comm = m_session->GetComm();

        // Get the proper time step with CFL control
        dt = GetSubstepTimeStep();

        nsubsteps = (m_timestep > dt)? ((int)(m_timestep/dt)+1):1; 
        nsubsteps = max(m_minsubsteps, nsubsteps);

        dt = m_timestep/nsubsteps;
        
        if (m_infosteps && !((nstep+1)%m_infosteps) && comm->GetRank() == 0)
        {
            cout << "Sub-integrating using "<< nsubsteps 
                 << " steps over Dt = "     << m_timestep 
                 << " (SubStep CFL="        << m_cflSafetyFactor << ")"<< endl;
        }

        for (int m = 0; m < nint; ++m)
        {
            // We need to update the fields held by the m_integrationSoln
            fields = integrationSoln->UpdateSolutionVector()[m];
            
            // Initialise NS solver which is set up to use a GLM method
            // with calls to EvaluateAdvection_SetPressureBCs and
            // SolveUnsteadyStokesSystem
            LibUtilities::TimeIntegrationSolutionSharedPtr 
                SubIntegrationSoln = m_subStepIntegrationScheme->
                InitializeScheme(dt, fields, time, m_subStepIntegrationOps);
            
            for(int n = 0; n < nsubsteps; ++n)
            {
                fields = m_subStepIntegrationScheme->TimeIntegrate(
                            n, dt,
                            SubIntegrationSoln,
                            m_subStepIntegrationOps);
            }
            
            // Reset time integrated solution in m_integrationSoln 
            integrationSoln->SetSolVector(m,fields);
        }
    }
    
    /**
     * Get maximum substep allowed.
     */
    NekDouble BlackScholes::GetSubstepTimeStep()
    { 
        int n_element = m_fields[0]->GetExpSize();

        const Array<OneD, int> ExpOrder=m_fields[0]->
                                EvalBasisNumModesMaxPerExp();
        Array<OneD, int> ExpOrderList(n_element, ExpOrder);
        
        const NekDouble cLambda = 0.2;
        
        Array<OneD, NekDouble> tstep      (n_element, 0.0);
        Array<OneD, NekDouble> stdVelocity(n_element, 0.0);

        stdVelocity = GetMaxStdVelocity(m_stockPrice);
        
        for (int el = 0; el < n_element; ++el)
        {
            tstep[el] = m_cflSafetyFactor / 
                (stdVelocity[el] * cLambda * 
                 (ExpOrder[el]-1) * (ExpOrder[el]-1));
        }
        
        NekDouble TimeStep = Vmath::Vmin(n_element, tstep, 1);
        m_session->GetComm()->AllReduce(TimeStep,LibUtilities::ReduceMin);        
        
        return TimeStep;
    }

    /**
     * Setup substepping time-integration.
     */
    void BlackScholes::SetUpSubSteppingTimeIntegration(
        int                                                 intMethod,
        const LibUtilities::TimeIntegrationWrapperSharedPtr &IntegrationScheme)
    {
        // Set to 1 for first step and it will then be increased in
        // time advance routines
        switch(intMethod)
        {
        case LibUtilities::eBackwardEuler:
        case LibUtilities::eBDFImplicitOrder1: 
            {
                m_subStepIntegrationScheme = LibUtilities::
                    GetTimeIntegrationWrapperFactory().
                        CreateInstance("ForwardEuler");
                
            }
            break;
        case LibUtilities::eBDFImplicitOrder2:
            {
                m_subStepIntegrationScheme = LibUtilities::
                    GetTimeIntegrationWrapperFactory().
                        CreateInstance("RungeKutta2_ImprovedEuler");
            }
            break;
        default:
            ASSERTL0(0, "Integration method not suitable: "
                        "Options include BackwardEuler or BDFImplicitOrder1");
            break;
        }
        m_intSteps = IntegrationScheme->GetIntegrationSteps();
	
        // set explicit time-integration class operators
        m_subStepIntegrationOps.DefineOdeRhs(&BlackScholes::SubStepAdvection, this);
        m_subStepIntegrationOps.DefineProjection(&BlackScholes::SubStepProjection, this);
    }
    
    /** 
     * Explicit Advection terms used by SubStepAdvance time integration.
     */
    void BlackScholes::SubStepAdvection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,  
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                  time)
    {
        int nVariables = inarray.num_elements();
        
        /// Get the number of coefficients
        int ncoeffs = m_fields[0]->GetNcoeffs(); 
        
        /// Define an auxiliary variable to compute the RHS 
        Array<OneD, Array<OneD, NekDouble> > WeakAdv(nVariables);
        WeakAdv[0] = Array<OneD, NekDouble> (ncoeffs*nVariables);
        for (int i = 1; i < nVariables; ++i)
        {
            WeakAdv[i] = WeakAdv[i-1] + ncoeffs;
        }
        
        // Currently assume velocity field is time independent
        // and does not therefore need extrapolating.
        // RHS computation using the advection base class
        m_advObject->Advect(nVariables, m_fields, m_stockPrice,
                            inarray, outarray, time);

        for (int i = 0; i < nVariables; ++i)
        {
            m_fields[i]->IProductWRTBase(outarray[i],WeakAdv[i]);
            
            // Negation required due to sign of DoAdvection term
            Vmath::Neg(ncoeffs, WeakAdv[i], 1);
        }
        
        AddAdvectionPenaltyFlux(m_stockPrice, inarray, WeakAdv);

        
        /// Operations to compute the RHS
        for (int i = 0; i < nVariables; ++i)
        {
            // Negate the RHS
            Vmath::Neg(ncoeffs, WeakAdv[i], 1);

            /// Multiply the flux by the inverse of the mass matrix
            m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i], WeakAdv[i]);
            
            /// Store in outarray the physical values of the RHS
            m_fields[i]->BwdTrans(WeakAdv[i], outarray[i]);
        }
    }
        
    /** 
     * Projection used by SubStepAdvance time integration.
     */
    void BlackScholes::SubStepProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,  
        Array<OneD, Array<OneD, NekDouble> > &outarray, 
        const NekDouble time)
    {
        ASSERTL1(inarray.num_elements() == outarray.num_elements(),
                 "Inarray and outarray of different sizes ");

        for (int i = 0; i < inarray.num_elements(); ++i)
        {
            Vmath::Vcopy(inarray[i].num_elements(),
                         inarray[i], 1, outarray[i], 1);
        }
    }

    /**
     * Penalty flux for the LDG diffusion term.
     */
    void BlackScholes::AddAdvectionPenaltyFlux(
                                const Array<OneD, const Array<OneD, NekDouble> > &velfield, 
                                const Array<OneD, const Array<OneD, NekDouble> > &physfield, 
                                Array<OneD, Array<OneD, NekDouble> > &Outarray)
    {
        ASSERTL1(physfield.num_elements() == Outarray.num_elements(),
                 "Physfield and outarray are of different dimensions");
        
        /// Number of trace points
        int nTracePts   = m_fields[0]->GetTrace()->GetNpoints();

        /// Forward state array
        Array<OneD, NekDouble> Fwd(3*nTracePts);
        
        /// Backward state array
        Array<OneD, NekDouble> Bwd = Fwd + nTracePts;

        /// upwind numerical flux state array
        Array<OneD, NekDouble> numflux = Bwd + nTracePts;
        
        /// Normal velocity array
        Array<OneD, NekDouble> Sn  = GetNormalStock(velfield);
        
        for (int i = 0; i < physfield.num_elements(); ++i)
        {
            /// Extract forwards/backwards trace spaces
            /// Note: Needs to have correct i value to get boundary conditions
            m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);
            
            /// Upwind between elements
            m_fields[0]->GetTrace()->Upwind(Sn, Fwd, Bwd, numflux);

            /// Construct difference between numflux and Fwd,Bwd
            Vmath::Vsub(nTracePts, numflux, 1, Fwd, 1, Fwd, 1);
            Vmath::Vsub(nTracePts, numflux, 1, Bwd, 1, Bwd, 1);

            /// Calculate the numerical fluxes multipling Fwd, Bwd and
            /// numflux by the normal advection velocity
            Vmath::Vmul(nTracePts, Fwd, 1, Sn, 1, Fwd, 1);
            Vmath::Vmul(nTracePts, Bwd, 1, Sn, 1, Bwd, 1);

            m_fields[0]->AddFwdBwdTraceIntegral(Fwd,Bwd,Outarray[i]);
        }
    }

    /**
     * Get max standard velocity for each element to calculate max dt.
     */
    Array<OneD, NekDouble> BlackScholes::GetMaxStdVelocity(
        const Array<OneD, Array<OneD,NekDouble> > inarray)
    {
        
        int n_points_0 = m_fields[0]->GetExp(0)->GetTotPoints();
        int n_element  = m_fields[0]->GetExpSize();
        int nvel       = inarray.num_elements();
        int cnt; 

        ASSERTL0(nvel >= 2, "Method not implemented for 1D");
        
        NekDouble pntVelocity;
        
        // Getting the standard velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(nvel);
        Array<OneD, NekDouble>               maxV(n_element, 0.0);
        LibUtilities::PointsKeyVector        ptsKeys;
        
        for (int i = 0; i < nvel; ++i)
        {
            stdVelocity[i] = Array<OneD, NekDouble>(n_points_0);
        }
        
        if (nvel == 2)
        {
            cnt = 0.0;
            for (int el = 0; el < n_element; ++el)
            { 
                int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();
                
                // reset local space if necessary
                if(n_points != n_points_0)
                {
                    for (int j = 0; j < nvel; ++j)
                    {
                        stdVelocity[j] = Array<OneD, NekDouble>(n_points);
                    }
                    n_points_0 = n_points;
                }		
                
                Array<TwoD, const NekDouble> gmat = 
                    m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->
                        GetDerivFactors(ptsKeys);
                
                if (m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->
                        GetGtype() == SpatialDomains::eDeformed)
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][i]*inarray[0][i+cnt] 
                            + gmat[2][i]*inarray[1][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][i]*inarray[0][i+cnt] 
                            + gmat[3][i]*inarray[1][i+cnt];
                    }
                }
                else
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][0]*inarray[0][i+cnt] 
                            + gmat[2][0]*inarray[1][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][0]*inarray[0][i+cnt] 
                            + gmat[3][0]*inarray[1][i+cnt];
                    }
                }
                
                cnt += n_points;
                
                
                for (int i = 0; i < n_points; i++)
                {
                    pntVelocity = stdVelocity[0][i]*stdVelocity[0][i] 
                        + stdVelocity[1][i]*stdVelocity[1][i];
                    
                    if (pntVelocity>maxV[el])
                    {
                        maxV[el] = pntVelocity;
                    }
                }
                maxV[el] = sqrt(maxV[el]);
            }
        }
        else
        {
            cnt = 0;
            for (int el = 0; el < n_element; ++el)
            {
                int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();
                
                // reset local space if necessary
                if (n_points != n_points_0)
                {
                    for (int j = 0; j < nvel; ++j)
                    {
                        stdVelocity[j] = Array<OneD, NekDouble>(n_points);
                    }
                    n_points_0 = n_points;
                }		
                
                Array<TwoD, const NekDouble> gmat =
                    m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->
                        GetDerivFactors(ptsKeys);
                
                if (m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->
                    GetGtype() == SpatialDomains::eDeformed)
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][i]*inarray[0][i+cnt] 
                            + gmat[3][i]*inarray[1][i+cnt] 
                            + gmat[6][i]*inarray[2][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][i]*inarray[0][i+cnt] 
                            + gmat[4][i]*inarray[1][i+cnt] 
                            + gmat[7][i]*inarray[2][i+cnt];
                        
                        stdVelocity[2][i] = gmat[2][i]*inarray[0][i+cnt] 
                            + gmat[5][i]*inarray[1][i+cnt] 
                            + gmat[8][i]*inarray[2][i+cnt];
                    }
                }
                else
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][0]*inarray[0][i+cnt] 
                            + gmat[3][0]*inarray[1][i+cnt] 
                            + gmat[6][0]*inarray[2][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][0]*inarray[0][i+cnt] 
                            + gmat[4][0]*inarray[1][i+cnt] 
                            + gmat[7][0]*inarray[2][i+cnt];
                        
                        stdVelocity[2][i] = gmat[2][0]*inarray[0][i+cnt] 
                            + gmat[5][0]*inarray[1][i+cnt] 
                            + gmat[8][0]*inarray[2][i+cnt];
                    }
                }
                cnt += n_points;
                
                for (int i = 0; i < n_points; i++)
                {
                    pntVelocity = stdVelocity[0][i]*stdVelocity[0][i] 
                        + stdVelocity[1][i]*stdVelocity[1][i] 
                        + stdVelocity[2][i]*stdVelocity[2][i];
                    
                    if (pntVelocity > maxV[el])
                    {
                        maxV[el] = pntVelocity;
                    }
                }
                maxV[el] = sqrt(maxV[el]);
            }
        }
        return maxV;
    }
}
