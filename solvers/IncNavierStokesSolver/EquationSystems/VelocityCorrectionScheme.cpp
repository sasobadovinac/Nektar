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
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc

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
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
        : UnsteadySystem(pSession, pGraph),
          IncNavierStokes(pSession, pGraph),
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

        // Integrate only the convective fields
        for (n = 0; n < m_nConvectiveFields; ++n)
        {
            m_intVariables.push_back(n);
        }

        SetUpExtrapolation();
        SetUpSVV();

        SetUpEVM();

        m_session->MatchSolverInfo("SmoothAdvection", "True",
                                    m_SmoothAdvection, false);

        // set explicit time-intregration class operators
        m_ode.DefineOdeRhs(
            &VelocityCorrectionScheme::EvaluateAdvection_SetPressureBCs, this);

        // set implicit time-intregration class operators
        m_ode.DefineImplicitSolve(
            &VelocityCorrectionScheme::SolveUnsteadyStokesSystem, this);

        // Set up bits for flowrate.
        m_session->LoadParameter("Flowrate", m_flowrate, 0.0);
        m_session->LoadParameter("IO_FlowSteps", m_flowrateSteps, 0);
    }

    void VelocityCorrectionScheme::SetUpExtrapolation()
    {
        // creation of the extrapolation object
        if (m_equationType == eUnsteadyNavierStokes ||
            m_equationType == eUnsteadyStokes)
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

            m_extrapolation->SubSteppingTimeIntegration(m_intScheme);
            m_extrapolation->GenerateHOPBCMap(m_session);
        }
    }

    /**
     * @brief Set up the Stokes solution used to impose constant flowrate
     * through a boundary.
     *
     * This routine solves a Stokes equation using a unit forcing direction,
     * specified by the user to be in the desired flow direction. This field can
     * then be used to correct the end of each timestep to impose a constant
     * volumetric flow rate through a user-defined boundary.
     *
     * There are three modes of operation:
     *
     * - Standard two-dimensional or three-dimensional simulations (e.g. pipes
     *   or channels)
     * - 3DH1D simulations where the forcing is not in the homogeneous
     *   direction (e.g. channel flow, where the y-direction of the 2D mesh
     *   is perpendicular to the wall);
     * - 3DH1D simulations where the forcing is in the homogeneous direction
     *   (e.g. pipe flow in the z-direction).
     *
     * In the first two cases, the user should define:
     * - the `Flowrate` parameter, which dictates the volumetric flux through
     *   the reference area
     * - tag a boundary region with the `Flowrate` user-defined type to define
     *   the reference area
     * - define a `FlowrateForce` function with components `ForceX`, `ForceY`
     *   and `ForceZ` that defines a unit forcing in the appropriate direction.
     *
     * In the latter case, the user should define only the `Flowrate`; the
     * reference area is taken to be the homogeneous plane and the force is
     * assumed to be the unit z-vector \f$ \hat{e}_z \f$.
     *
     * This routine solves a single timestep of the Stokes problem
     * (premultiplied by the backwards difference coefficient):
     *
     * \f[ \frac{\partial\mathbf{u}}{\partial t} = -\nabla p +
     * \nu\nabla^2\mathbf{u} + \mathbf{f} \f]
     *
     * with a zero initial condition to obtain a field \f$ \mathbf{u}_s \f$. The
     * flowrate is then corrected at each timestep \f$ n \f$ by adding the
     * correction \f$ \alpha\mathbf{u}_s \f$ where
     *
     * \f[ \alpha = \frac{\overline{Q} - Q(\mathbf{u^n})}{Q(\mathbf{u}_s)} \f]
     *
     * where \f$ Q(\cdot)\f$ is the volumetric flux through the appropriate
     * surface or line, which is implemented in
     * VelocityCorrectionScheme::MeasureFlowrate. For more details, see chapter
     * 3.2 of the thesis of D. Moxey (University of Warwick, 2011).
     */
    void VelocityCorrectionScheme::SetupFlowrate(NekDouble aii_dt)
    {
        m_flowrateBndID = -1;
        m_flowrateArea = 0.0;

        const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bcs =
            m_fields[0]->GetBndConditions();

        std::string forces[] = { "X", "Y", "Z" };
        Array<OneD, NekDouble> flowrateForce(m_spacedim, 0.0);

        // Set up flowrate forces.
        bool defined = true;
        for (int i = 0; i < m_spacedim; ++i)
        {
            std::string varName = std::string("Force") + forces[i];
            defined = m_session->DefinesFunction("FlowrateForce", varName);

            if (!defined && m_HomogeneousType == eHomogeneous1D)
            {
                break;
            }

            ASSERTL0(defined,
                     "A 'FlowrateForce' function must defined with components "
                     "[ForceX, ...] to define direction of flowrate forcing");

            LibUtilities::EquationSharedPtr ffunc
                = m_session->GetFunction("FlowrateForce", varName);
            flowrateForce[i] = ffunc->Evaluate();
        }

        // Define flag for case with homogeneous expansion and forcing not in the
        // z-direction
        m_homd1DFlowinPlane = false;
        if (defined && m_HomogeneousType == eHomogeneous1D)
        {
            m_homd1DFlowinPlane = true;
        }

        // For 3DH1D simulations, if force isn't defined then assume in
        // z-direction.
        if (!defined)
        {
            flowrateForce[2] = 1.0;
        }

        // Find the boundary condition that is tagged as the flowrate boundary.
        for (int i = 0; i < bcs.num_elements(); ++i)
        {
            if (boost::iequals(bcs[i]->GetUserDefined(), "Flowrate"))
            {
                m_flowrateBndID = i;
                break;
            }
        }

        int tmpBr = m_flowrateBndID;
        m_comm->AllReduce(tmpBr, LibUtilities::ReduceMax);
        ASSERTL0(tmpBr >= 0 || m_HomogeneousType == eHomogeneous1D,
                 "One boundary region must be marked using the 'Flowrate' "
                 "user-defined type to monitor the volumetric flowrate.");

        // Extract an appropriate expansion list to represents the boundary.
        if (m_flowrateBndID >= 0)
        {
            // For a boundary, extract the boundary itself.
            m_flowrateBnd = m_fields[0]->GetBndCondExpansions()[m_flowrateBndID];
        }
        else if (m_HomogeneousType == eHomogeneous1D && !m_homd1DFlowinPlane)
        {
            // For 3DH1D simulations with no force specified, find the mean
            // (0th) plane.
            Array<OneD, unsigned int> zIDs = m_fields[0]->GetZIDs();
            int tmpId = -1;

            for (int i = 0; i < zIDs.num_elements(); ++i)
            {
                if (zIDs[i] == 0)
                {
                    tmpId = i;
                    break;
                }
            }

            ASSERTL1(tmpId <= 0, "Should be either at location 0 or -1 if not "
                                 "found");

            if (tmpId != -1)
            {
                m_flowrateBnd = m_fields[0]->GetPlane(tmpId);
            }
        }

        // At this point, some processors may not have m_flowrateBnd set if they
        // don't contain the appropriate boundary. To calculate the area, we
        // integrate 1.0 over the boundary (which has been set up with the
        // appropriate subcommunicator to avoid deadlock), and then communicate
        // this to the other processors with an AllReduce.
        if (m_flowrateBnd)
        {
            Array<OneD, NekDouble> inArea(m_flowrateBnd->GetNpoints(), 1.0);
            m_flowrateArea = m_flowrateBnd->Integral(inArea);
        }
        m_comm->AllReduce(m_flowrateArea, LibUtilities::ReduceMax);

        // In homogeneous case with forcing not aligned to the z-direction,
        // redefine m_flowrateBnd so it is a 1D expansion
        if (m_HomogeneousType == eHomogeneous1D && m_homd1DFlowinPlane &&
            m_flowrateBnd)
        {
            // For 3DH1D simulations with no force specified, find the mean
            // (0th) plane.
            Array<OneD, unsigned int> zIDs = m_fields[0]->GetZIDs();
            m_planeID = -1;

            for (int i = 0; i < zIDs.num_elements(); ++i)
            {
                if (zIDs[i] == 0)
                {
                    m_planeID = i;
                    break;
                }
            }

            ASSERTL1(m_planeID <= 0, "Should be either at location 0 or -1 if not "
                                 "found");

            if (m_planeID != -1)
            {
                m_flowrateBnd = m_fields[0]
                                    ->GetBndCondExpansions()[m_flowrateBndID]
                                    ->GetPlane(m_planeID);
            }
        }

        // Set up some storage for the Stokes solution (to be stored in
        // m_flowrateStokes) and its initial condition (inTmp), which holds the
        // unit forcing.
        int nqTot = m_fields[0]->GetNpoints();
        Array<OneD, Array<OneD, NekDouble> > inTmp(m_spacedim);
        m_flowrateStokes = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);

        for (int i = 0; i < m_spacedim; ++i)
        {
            inTmp[i] = Array<OneD, NekDouble>(
                nqTot, flowrateForce[i] * aii_dt);
            m_flowrateStokes[i] = Array<OneD, NekDouble>(nqTot, 0.0);

            if (m_HomogeneousType == eHomogeneous1D)
            {
                Array<OneD, NekDouble> inTmp2(nqTot);
                m_fields[i]->HomogeneousFwdTrans(inTmp[i], inTmp2);
                m_fields[i]->SetWaveSpace(true);
                inTmp[i] = inTmp2;
            }

            Vmath::Zero(
                m_fields[i]->GetNcoeffs(), m_fields[i]->UpdateCoeffs(), 1);
        }

        // Create temporary extrapolation object to avoid issues with
        // m_extrapolation for HOPBCs using higher order timestepping schemes.
        ExtrapolateSharedPtr tmpExtrap = m_extrapolation;
        m_extrapolation = GetExtrapolateFactory().CreateInstance(
            "Standard", m_session, m_fields, m_pressure, m_velocity,
            m_advObject);

        // Finally, calculate the solution and the flux of the Stokes
        // solution. We set m_greenFlux to maximum numeric limit, which signals
        // to SolveUnsteadyStokesSystem that we don't need to apply a flowrate
        // force.
        m_greenFlux = numeric_limits<NekDouble>::max();
        m_flowrateAiidt = aii_dt;
        SolveUnsteadyStokesSystem(inTmp, m_flowrateStokes, 0.0, aii_dt);
        m_greenFlux = MeasureFlowrate(m_flowrateStokes);

        // If the user specified IO_FlowSteps, open a handle to store output.
        if (m_comm->GetRank() == 0 && m_flowrateSteps &&
            !m_flowrateStream.is_open())
        {
            std::string filename = m_session->GetSessionName();
            filename += ".prs";
            m_flowrateStream.open(filename.c_str());
            m_flowrateStream.setf(ios::scientific, ios::floatfield);
            m_flowrateStream << "# step      time            dP" << endl
                             << "# -------------------------------------------"
                             << endl;
        }

        m_extrapolation = tmpExtrap;
    }

    /**
     * @brief Measure the volumetric flow rate through the volumetric flow rate
     * reference surface.
     *
     * This routine computes the volumetric flow rate
     *
     * \f[
     * Q(\mathbf{u}) = \frac{1}{\mu(R)} \int_R \mathbf{u} \cdot d\mathbf{s}
     * \f]
     *
     * through the boundary region \f$ R \f$.
     */
    NekDouble VelocityCorrectionScheme::MeasureFlowrate(
        const Array<OneD, Array<OneD, NekDouble> > &inarray)
    {
        NekDouble flowrate = 0.0;

        if (m_flowrateBnd && m_flowrateBndID >= 0)
        {
            // If we're an actual boundary, calculate the vector flux through
            // the boundary.
            Array<OneD, Array<OneD, NekDouble> > boundary(m_spacedim);

            if(!m_homd1DFlowinPlane)
            {
                // General case
                for (int i = 0; i < m_spacedim; ++i)
                {
                    m_fields[i]->ExtractPhysToBnd(m_flowrateBndID, inarray[i],
                                                  boundary[i]);
                }
                flowrate = m_flowrateBnd->VectorFlux(boundary);
            }
            else if(m_planeID == 0)
            {
                //Homogeneous with forcing in plane. Calculate flux only on
                // the meanmode - calculateFlux necessary for hybrid
                // parallelisation.
                for (int i = 0; i < m_spacedim; ++i)
                {
                    m_fields[i]->GetPlane(m_planeID)->ExtractPhysToBnd(
                        m_flowrateBndID, inarray[i], boundary[i]);
                }

                // the flowrate is calculated on the mean mode so it needs to be
                // multiplied by LZ to be consistent with the general case.
                flowrate = m_flowrateBnd->VectorFlux(boundary) *
                           m_session->GetParameter("LZ");
            }
        }
        else if (m_flowrateBnd && !m_homd1DFlowinPlane)
        {
            // 3DH1D case with no Flowrate boundary defined: compute flux
            // through the zero-th (mean) plane.
            flowrate = m_flowrateBnd->Integral(inarray[2]);
        }

        // Communication to obtain the total flowrate
        if(!m_homd1DFlowinPlane && m_HomogeneousType == eHomogeneous1D)
        {
            m_comm->GetColumnComm()->AllReduce(flowrate, LibUtilities::ReduceSum);
        }
        else
        {
            m_comm->AllReduce(flowrate, LibUtilities::ReduceSum); 
        }
        return flowrate / m_flowrateArea;
    }

    bool VelocityCorrectionScheme::v_PostIntegrate(int step)
    {
        if (m_flowrateSteps > 0)
        {
            if (m_comm->GetRank() == 0 && (step + 1) % m_flowrateSteps == 0)
            {
                m_flowrateStream << setw(8) << step << setw(16) << m_time
                                 << setw(16) << m_alpha << endl;
            }
        }

        return IncNavierStokes::v_PostIntegrate(step);
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

        if (m_extrapolation->GetSubStepName().size() )
        {
            SolverUtils::AddSummaryItem(s, "Substepping",
                              m_extrapolation->GetSubStepName() );
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
        m_F = Array<OneD, Array<OneD, NekDouble> > (m_nConvectiveFields);

        for (int i = 0; i < m_nConvectiveFields; ++i)
        {
            m_F[i] = Array< OneD, NekDouble> (m_fields[0]->GetTotPoints(), 0.0);
        }

        m_flowrateAiidt = 0.0;

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

        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            m_fields[i]->LocalToGlobal();
            m_fields[i]->ImposeDirichletConditions(m_fields[i]->UpdateCoeffs());
            m_fields[i]->GlobalToLocal();
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                  m_fields[i]->UpdatePhys());
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
        // Set up flowrate if we're starting for the first time or the value of
        // aii_Dt has changed.
        if (m_flowrate > 0.0 && (aii_Dt != m_flowrateAiidt))
        {
            SetupFlowrate(aii_Dt);
        }

        if(m_useEntropyViscosity)
           BackUpSolution(inarray);

        int physTot = m_fields[0]->GetTotPoints();

        // Substep the pressure boundary condition if using substepping
        m_extrapolation->SubStepSetPressureBCs(inarray,aii_Dt,m_kinvis);

        // Set up forcing term for pressure Poisson equation
        SetUpPressureForcing(inarray, m_F, aii_Dt);

        // Solve Pressure System
        SolvePressure (m_F[0]);

        // Set up forcing term for Helmholtz problems
        SetUpViscousForcing(inarray, m_F, aii_Dt);


        if(m_useEntropyViscosity)
        {
//           BackUpSolution();
          if(m_step_counter>m_evm_start_step)
           {
             SetUpEntropyViscosity(inarray, m_F);
           }
           m_step_counter ++;
        }

        // Solve velocity system
        SolveViscous( m_F, outarray, aii_Dt);

        // Apply flowrate correction
        if (m_flowrate > 0.0 && m_greenFlux != numeric_limits<NekDouble>::max())
        {
            NekDouble currentFlux = MeasureFlowrate(outarray);
            m_alpha = (m_flowrate - currentFlux) / m_greenFlux;

            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Svtvp(physTot, m_alpha, m_flowrateStokes[i], 1,
                             outarray[i], 1, outarray[i], 1);
            }
        }
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

        AppendSVVFactors(factors,varFactorsMap);

        // Solve Helmholtz system and put in Physical space
        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            // Setup coefficients for equation
            factors[StdRegions::eFactorLambda] = 1.0/aii_Dt/m_diffCoeff[i];
            m_fields[i]->HelmSolve(Forcing[i], m_fields[i]->UpdateCoeffs(),
                                   NullFlagList,  factors, varCoeffMap,
                                   varFactorsMap);
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),outarray[i]);
        }
    }

    void  VelocityCorrectionScheme::SetUpSVV(void)
    {

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
        if(m_session->DefinesSolverInfo("SpectralVanishingViscosityHomo1D"))
        {
            m_session->MatchSolverInfo("SpectralVanishingViscosityHomo1D",
                                "True", m_useHomo1DSpecVanVisc, false);
            if(m_useHomo1DSpecVanVisc == false)
            {
                m_session->MatchSolverInfo("SpectralVanishingViscosityHomo1D",
                                       "ExpKernel", m_useHomo1DSpecVanVisc, false);
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

                for(int n = 0; n < num_planes; ++n)
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

                int p = 0;
                for(int i = 0; i < 3; ++i)
                {
                    p = max(p,exp3D->GetBasisNumModes(i)-1);
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

                int p = 0;
                for(int i = 0; i < 2; ++i)
                {
                    p = max(p,exp2D->GetBasisNumModes(i)-1);
                }

                diffcoeff[e] *= h/p;
            }
        }
    }

    void VelocityCorrectionScheme::AppendSVVFactors(
                                 StdRegions::ConstFactorMap &factors,
                                 MultiRegions::VarFactorsMap &varFactorsMap)
    {

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

    }

    void  VelocityCorrectionScheme::SetUpEVM(void)
    {

        m_session->MatchSolverInfo("EntropyViscosityMethod",
                                   "True", m_useEntropyViscosity, false);
        m_session->MatchSolverInfo("DynamicEVM",
                                   "True", m_use_dynamic_alpha, false);
      
       // Get EVM paramter
       m_session->LoadParameter("evm_alpha", m_evm_alpha, 0.1);
       m_session->LoadParameter("evm_beta", m_evm_beta, 0.1);
       m_session->LoadParameter("evm_start_step", m_evm_start_step, 2);
       m_session->LoadParameter("evm_reduced_order", m_evm_reduced_order, 1);

       m_step_counter = 0;

       m_use_evm_mapping = false;
       
       int nfields = m_fields.num_elements();
       int phystot = m_fields[0]->GetTotPoints();

        m_solution     = Array<OneD, Array<OneD, NekDouble> >(nfields);
        m_old_solution = Array<OneD, Array<OneD, NekDouble> >(nfields);
   
        if(m_useEntropyViscosity)
         {
          for(int i=0; i<nfields; ++i)
           {
            m_solution[i] = Array<OneD, NekDouble>(phystot);
            m_old_solution[i] = Array<OneD, NekDouble>(phystot);

            Vmath::Zero(phystot, m_solution[i], 1);
            Vmath::Zero(phystot, m_old_solution[i], 1);
           }
            m_evm_visc = Array<OneD, NekDouble>(phystot);
            Vmath::Zero(phystot, m_evm_visc, 1);
          
           if(m_use_dynamic_alpha)
            {
              m_evm_sensor = Array<OneD, NekDouble>(phystot);
              Vmath::Zero(phystot, m_evm_sensor, 1);
           }
        }

   }

   void VelocityCorrectionScheme:: v_BackUpSolution( 
            const Array<OneD, const Array<OneD, NekDouble> > &inarray) 
    {
       int nfields = m_fields.num_elements();
       int phystot = m_fields[0]->GetTotPoints();
       int VelDim  = inarray.num_elements();

       for(int i=0; i<nfields; ++i)
         Vmath::Vcopy(phystot, m_solution[i], 1, m_old_solution[i], 1);
         
//       for(int i=0; i<nfields; ++i)
//         Vmath::Vcopy(phystot, m_fields[i]->GetPhys(), 1, m_solution[i], 1);

        for(int i=0; i<VelDim; ++i)
         Vmath::Vcopy(phystot, inarray[i], 1, m_solution[i], 1);
         
        m_pressure->BwdTrans(m_pressure->GetCoeffs(), m_solution[VelDim]);
//         Vmath::Vcopy(phystot, m_pressure, 1, m_solution[VelDim], 1);
     }

    std::pair<NekDouble, NekDouble >
    VelocityCorrectionScheme::
    v_EvaluateVelocityRange(const Array<OneD, const Array<OneD, NekDouble> > &velocity)
    {
        int nqtot      = m_fields[0]->GetTotPoints();
        int VelDim     = m_velocity.num_elements();

        NekDouble max_velocity,min_velocity;
        
        max_velocity = -numeric_limits<double>::max();
        min_velocity =  numeric_limits<double>::max();

        for(unsigned int j = 0; j < nqtot; ++j)
          {
             NekDouble velocity_mag = 0.;
            for(unsigned int i = 0; i < VelDim; ++i)
                velocity_mag += velocity[i][j]*velocity[i][j];
            velocity_mag = sqrt(velocity_mag);

            max_velocity = max(velocity_mag,max_velocity);
            min_velocity = min(velocity_mag,min_velocity);
          }

         min_velocity *=-1.;

         m_session->GetComm()->AllReduce(max_velocity,LibUtilities::ReduceMax);
         m_session->GetComm()->AllReduce(min_velocity,LibUtilities::ReduceMax);
         
         min_velocity *=-1.;

        return std::make_pair(min_velocity,max_velocity);

    }

    NekDouble VelocityCorrectionScheme::
    v_EvaluateVelocityEntropyVariation(const Array<OneD, const Array<OneD, NekDouble> > &velocity)
    {
        int VelDim     = m_velocity.num_elements();
        StdRegions::StdExpansionSharedPtr elmt;
        
        NekDouble mean_velocity;

        const std::pair<NekDouble, NekDouble> & velocity_range = EvaluateVelocityRange(velocity);

        mean_velocity = 0.5*(velocity_range.first+velocity_range.second);

        double entropy_integrate = 0.,
               element_volume = 0.;
        double max_entropy = -numeric_limits<double>::max();
        double min_entropy =  numeric_limits<double>::max();

        switch(VelDim)
        {
         case 1:
            ASSERTL0(false,"EVM not implemented for one dimensional problem");
           break;
         case 2:
          {
              int nel =  m_fields[0]->GetExpSize();
              int nq  =  m_fields[0]->GetExp(0)->GetTotPoints();

              int cnt = 0;
              for(int i=0; i<nel; ++i,++cnt)
                {
                  elmt   = m_fields[0]->GetExp(i);
                  int offset = m_fields[0]->GetPhys_Offset(cnt);

                  Array<OneD, NekDouble>  entropy_list(nq,0.0);
                  Array<OneD, NekDouble>  one_list(nq,1.);

                  for(int q=0; q<nq; ++q)
                    {
                         NekDouble velocity_mag = 0.;
                         for(unsigned int d = 0; d < VelDim; ++d)
                             velocity_mag += velocity[d][offset+q]*velocity[d][offset+q];
                         velocity_mag = sqrt(velocity_mag);

                      entropy_list[q] = (velocity_mag-mean_velocity)*(velocity_mag-mean_velocity);

                      max_entropy = max(max_entropy,entropy_list[q]);
                      min_entropy = min(min_entropy,entropy_list[q]);
                    }
                  
                  entropy_integrate += elmt->Integral(entropy_list);
                  element_volume  += elmt->Integral(one_list);

                }
           break;
          }
         case 3:
          if( m_fields[0]->GetWaveSpace() == false )
           {
              int nel =  m_fields[0]->GetExpSize();
              int nq  =  m_fields[0]->GetExp(0)->GetTotPoints();

              int cnt = 0;
              for(int i=0; i<nel; ++i,++cnt)
                {
                  elmt   = m_fields[0]->GetExp(i);
                  int offset = m_fields[0]->GetPhys_Offset(cnt);

                  Array<OneD, NekDouble>  entropy_list(nq,0.0);
                  Array<OneD, NekDouble>  one_list(nq,1.);

                  for(int q=0; q<nq; ++q)
                    {
                         NekDouble velocity_mag = 0.;
                         for(unsigned int d = 0; d < VelDim; ++d)
                             velocity_mag += velocity[d][offset+q]*velocity[d][offset+q];
                         velocity_mag = sqrt(velocity_mag);

                      entropy_list[q] = (velocity_mag-mean_velocity)*(velocity_mag-mean_velocity);

                      max_entropy = max(max_entropy,entropy_list[q]);
                      min_entropy = min(min_entropy,entropy_list[q]);
                    }
                  
                  entropy_integrate += elmt->Integral(entropy_list);
                  element_volume  += elmt->Integral(one_list);
                }
           }
          else
           {
            if( m_fields[0]->GetExpType() == MultiRegions::e3DH1D )
             {
              Array<OneD, unsigned int> planes;
              planes = m_fields[0]->GetZIDs();

              int num_planes = planes.num_elements();

              //int nt  =  m_fields[0]->GetPlane(0)->GetNpoints();
              int nel =  m_fields[0]->GetPlane(0)->GetExpSize();
              int nq  =  m_fields[0]->GetPlane(0)->GetExp(0)->GetTotPoints();

              int cnt = 0;
              for(unsigned int n=0; n<num_planes; ++n)
               {
                for(int i=0; i<nel; ++i,++cnt)
                 {
                  elmt   = m_fields[0]->GetPlane(n)->GetExp(i);
                  int offset = m_fields[0]->GetPhys_Offset(cnt);

                  Array<OneD, NekDouble>  entropy_list(nq,0.0);
                  Array<OneD, NekDouble>  one_list(nq,1.);

                  for(int q=0; q<nq; ++q)
                    {
                         NekDouble velocity_mag = 0.;
                         for(unsigned int d = 0; d < VelDim; ++d)
                             velocity_mag += velocity[d][offset+q]*velocity[d][offset+q];
                         velocity_mag = sqrt(velocity_mag);

                      entropy_list[q] = (velocity_mag-mean_velocity)*(velocity_mag-mean_velocity);

                      max_entropy = max(max_entropy,entropy_list[q]);
                      min_entropy = min(min_entropy,entropy_list[q]);
                    }
                  
                  entropy_integrate += elmt->Integral(entropy_list);
                  element_volume  += elmt->Integral(one_list);

                 }
               }
             }
            else
             {
                 ASSERTL0(false,"only fully 3D and 3DH1D are implemented ");
             }
           }

           break;
        default:
          {
              ASSERTL1(false,"Not Implemented !");
          }
        }


        
         m_session->GetComm()->AllReduce(entropy_integrate,LibUtilities::ReduceSum);
         m_session->GetComm()->AllReduce(element_volume,LibUtilities::ReduceSum);

         double mean_entropy = entropy_integrate/element_volume;
         m_session->GetComm()->AllReduce(max_entropy,LibUtilities::ReduceMax);

         min_entropy = -min_entropy;
         m_session->GetComm()->AllReduce(min_entropy,LibUtilities::ReduceMax);
    //     min_entropy =-min_entropy;

         return  max(max_entropy-mean_entropy, mean_entropy-(-min_entropy)); 


    }
    
    void
    VelocityCorrectionScheme::
    v_SetUpEntropyViscosity(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                                  Array<OneD, Array<OneD, NekDouble> > &outarray)
    {     
//          int num_processes = m_session->GetComm()->GetSize();
          int VelDim  = m_velocity.num_elements();
          int ndim    = m_velocity.num_elements();

          NekDouble BBeta[3][3] = {{ 1.0,  0.0, 0.0},{ 2.0, -1.0, 0.0},{ 3.0, -3.0, 1.0}};
          
         // int type = (int) flow_type;
         // Timer     timer;
          NekDouble A = 26.;
          NekDouble damping_length = 50;
          NekDouble damping_factor = 1.;
          NekDouble y_plus = 0.;
          NekDouble max_viscosity = 0.0, max_alpha = 0.0;
          NekDouble max_velocity  = 0.0;

          int nqtots     =  m_fields[0]->GetNpoints();
          int nctots     =  m_fields[0]->GetNcoeffs();
          int nElmts     =  m_fields[0]->GetExpSize();
          int nqElmts    =  m_fields[0]->GetExp(0)->GetTotPoints();

          int nqPtots    =  m_fields[0]->GetNpoints();
          int ncPtots    =  m_fields[0]->GetNcoeffs();
           
          Array<OneD, int> expOrderElement; 
         if(m_fields[0]->GetWaveSpace())
            expOrderElement = m_fields[0]->GetPlane(0)->EvalBasisNumModesMaxPerExp();
         else
            expOrderElement = m_fields[0]->EvalBasisNumModesMaxPerExp();

          int evm_reduced_order = m_evm_reduced_order;

          NekDouble m_homoLen ;
          LibUtilities::TranspositionSharedPtr m_trans;
          Array<OneD, unsigned int> planes;
          int NZ;
         if(m_fields[0]->GetWaveSpace())
           {
            planes     =  m_fields[0]->GetZIDs();
            NZ         =  planes.num_elements();

            nElmts     =  m_fields[0]->GetPlane(0)->GetExpSize();
            nqElmts    =  m_fields[0]->GetPlane(0)->GetExp(0)->GetTotPoints();
            nqPtots    =  m_fields[0]->GetPlane(0)->GetNpoints();
            ncPtots    =  m_fields[0]->GetPlane(0)->GetNcoeffs();
             
            m_trans    =  m_fields[0]->GetTransposition();
            m_homoLen  =  m_fields[0]->GetHomoLen();
           }



         // bool use_van_driest_damping = false;
         // m_session->MatchSolverInfo("VanDriestDamping","True",use_van_driest_damping,false);
         // m_session->LoadParameter("DampingLength", damping_length, 50.);
         // m_session->LoadParameter("HomModesZ", num_points_z, 32);

          Array<OneD, Array<OneD, NekDouble> > velocity_solution(VelDim),old_velocity_solution(VelDim) ;
          Array<OneD, Array<OneD, NekDouble> > velocity(VelDim),old_velocity(VelDim);
          Array<OneD, NekDouble> pressure_solution, old_pressure_solution;
          Array<OneD, NekDouble> pressure, old_pressure;

        

        
         for(int d = 0; d < VelDim; ++d)
          {
            velocity_solution[d] = Array<OneD, NekDouble>(nqtots,0.0);
            old_velocity_solution[d] = Array<OneD, NekDouble>(nqtots,0.0);
            
            Vmath::Vcopy(nqtots, m_solution[d], 1, velocity_solution[d], 1);
            Vmath::Vcopy(nqtots, m_old_solution[d], 1, old_velocity_solution[d], 1);

  //           m_fields[d]->BwdTrans(m_intSoln->GetValue(0)[d],velocity_solution[d]);
  //           m_fields[d]->BwdTrans(m_intSoln->GetValue(1)[d],old_velocity_solution[d]);
          }

          {
             pressure_solution = Array<OneD, NekDouble>(nqtots,0.0);
             old_pressure_solution = Array<OneD, NekDouble>(nqtots,0.0);

            Vmath::Vcopy(nqtots, m_solution[VelDim], 1, pressure_solution, 1);
            Vmath::Vcopy(nqtots, m_old_solution[VelDim], 1, old_pressure_solution, 1);
//             m_fields[VelDim]->BwdTrans(m_intSoln->GetValue(0)[VelDim],pressure_solution);
//             m_fields[VelDim]->BwdTrans(m_intSoln->GetValue(1)[VelDim],old_pressure_solution);
          }


         if(m_fields[0]->GetWaveSpace())
          {
            for(int d = 0; d < VelDim; ++d)
             {
              velocity[d] = Array<OneD, NekDouble>(nqtots,0.0);
              old_velocity[d] = Array<OneD, NekDouble>(nqtots,0.0);

              m_fields[d]->HomogeneousBwdTrans(velocity_solution[d],velocity[d]);
              m_fields[d]->HomogeneousBwdTrans(old_velocity_solution[d],old_velocity[d]);
             }

              pressure = Array<OneD, NekDouble>(nqtots,0.0);
              old_pressure = Array<OneD, NekDouble>(nqtots,0.0);

              m_fields[0]->HomogeneousBwdTrans(pressure_solution,pressure);
              m_fields[0]->HomogeneousBwdTrans(old_pressure_solution,old_pressure);
             
          }
         else
          {
             for(int d = 0; d < VelDim; ++d)
              {
               velocity[d] = velocity_solution[d];
               old_velocity[d] = old_velocity_solution[d];
              }

              pressure = pressure_solution;
              old_pressure = old_pressure_solution;
          }
       
       // NOTE upon here 
       // velocity[], old_velocity[], pressure and old_pressure store velocity and pressure in (Q,P)-space 
       //

       //additional working space
        Array<OneD, NekDouble> wkSp0,wkSp1,wkSp2;
        wkSp0 = Array<OneD, NekDouble>(nqtots,0.0);
        wkSp1 = Array<OneD, NekDouble>(nqtots,0.0);
        wkSp2 = Array<OneD, NekDouble>(nqtots,0.0);

        Array<OneD, NekDouble> Grad0,Grad1,Grad2;
        Array<OneD, NekDouble> wkCoef0,wkCoef1,wkCoef2;

        Array<OneD, NekDouble> entropy_residual;
        entropy_residual = Array<OneD, NekDouble>(nqtots,0.0);
        Vmath::Zero(nqtots, entropy_residual, 1);

  
        switch(ndim)
        {
         case 1:
            ASSERTL0(false,"EVM not implemented for one dimensional problem");
            //U = 0.5*(u*u)   

            Vmath::Vmul(nqtots,velocity[0],1,velocity[0],1, wkSp0,1);
            Vmath::Smul(nqtots, 0.5, wkSp0, 1, wkSp0, 1); // U

            Vmath::Vmul(nqtots, old_velocity[0],1, old_velocity[0],1, wkSp1,1);
            Vmath::Smul(nqtots, 0.5, wkSp1, 1, wkSp1,1); // U at previous time step

            Vmath::Vsub(nqtots, wkSp0,1, wkSp1,1,entropy_residual,1);
            Vmath::Smul(nqtots,1./m_timestep,entropy_residual,1,entropy_residual,1); // dU/dt

            Vmath::Vadd(nqtots, wkSp0, 1, wkSp1, 1, wkSp2,1);
            Vmath::Smul(nqtots, 0.5, wkSp2, 1, wkSp2,1); // U at (n+1/2)
            
            Vmath::Vadd(nqtots, pressure, 1, old_pressure, 1, wkSp0,1); //p at (n+1/2)
            Vmath::Vadd(nqtots, wkSp0, 1, wkSp2, 1, wkSp0,1); //(p+U)

            //wkSp1 is free to use
            Vmath::Vadd(nqtots, velocity[0], 1, old_velocity[0], 1, wkSp1,1);
            Vmath::Smul(nqtots, 0.5, wkSp1, 1, wkSp1,1); // u at (n+1/2)

            Vmath::Vmul(nqtots, wkSp0,1, wkSp1,1, wkSp0,1); //(U+p)*u

            break;
        case 2:
          {
            Grad0 = Array<OneD, NekDouble>(nqtots,0.0);
            Grad1 = Array<OneD, NekDouble>(nqtots,0.0);
            //U = 0.5*(u*u)   
            Vmath::Vmul(nqtots,velocity[0], 1, velocity[0], 1, wkSp0, 1);
            Vmath::Vvtvp(nqtots,velocity[1], 1, velocity[1], 1, wkSp0, 1, wkSp0, 1); //U
            Vmath::Smul(nqtots, 0.5, wkSp0, 1, wkSp0, 1); // U at (n+1)

            Vmath::Vmul(nqtots, old_velocity[0],1, old_velocity[0],1, wkSp1,1);
            Vmath::Vvtvp(nqtots, old_velocity[1], 1, old_velocity[1], 1, wkSp1, 1, wkSp1, 1); //U
            Vmath::Smul(nqtots, 0.5, wkSp1, 1, wkSp1,1); // U at (n)

            Vmath::Vsub(nqtots, wkSp0,1, wkSp1,1,entropy_residual,1);
            Vmath::Smul(nqtots,1./m_timestep,entropy_residual,1,entropy_residual,1); // dU/dt
            
            Vmath::Vadd(nqtots, wkSp0, 1, wkSp1, 1, wkSp2,1);
            Vmath::Smul(nqtots, 0.5, wkSp2, 1, wkSp2,1); // U at (n+1/2)

            Vmath::Vadd(nqtots, pressure, 1, old_pressure, 1, wkSp0,1); 
            Vmath::Smul(nqtots, 0.5, wkSp0, 1, wkSp0,1); // P at (n+1/2)
            Vmath::Vadd(nqtots, wkSp0, 1, wkSp2, 1, wkSp0,1); //(p+U) at (n+1/2)

            //wkSp1 is free to use
            Vmath::Vadd(nqtots, velocity[0], 1, old_velocity[0], 1, wkSp1,1);
            Vmath::Smul(nqtots, 0.5, wkSp1, 1, wkSp1,1); // u at (n+1/2)
            Vmath::Vmul(nqtots, wkSp0,1, wkSp1,1, wkSp1, 1); //(U+p)*u
             
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], wkSp1, Grad0);

            Vmath::Vadd(nqtots, velocity[1], 1, old_velocity[1], 1, wkSp1,1);
            Vmath::Smul(nqtots, 0.5, wkSp1, 1, wkSp1,1); // v at (n+1/2)
            Vmath::Vmul(nqtots, wkSp0,1, wkSp1,1, wkSp1, 1); //(U+p)*v
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], wkSp1, Grad1);

            Vmath::Vadd(nqtots, entropy_residual,1,  Grad0,1, entropy_residual, 1);
            Vmath::Vadd(nqtots, entropy_residual,1,  Grad1,1, entropy_residual, 1);

            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], wkSp2, wkSp1);
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], wkSp1, Grad0);

            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], wkSp2, wkSp1);
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], wkSp1, Grad1);

            Vmath::Vadd(nqtots, Grad0, 1,  Grad1, 1, Grad0, 1);
            Vmath::Smul(nqtots,m_kinvis, Grad0, 1, Grad0,1); // laplacian term
            Vmath::Vsub(nqtots, entropy_residual,1,  Grad0, 1, entropy_residual, 1);

            //
            Vmath::Vadd(nqtots, velocity[0], 1, old_velocity[0], 1, wkSp0,1);
            Vmath::Smul(nqtots, 0.5, wkSp0, 1, wkSp0,1); // u at (n+1/2)
            m_fields[0]->PhysDeriv(wkSp0, Grad0, Grad1);
            Vmath::Vmul(nqtots, Grad0,1,  Grad0, 1, wkSp0, 1);
            Vmath::Vvtvp(nqtots, Grad1, 1, Grad1, 1, wkSp0, 1, wkSp0, 1); //
            Vmath::Smul(nqtots, m_kinvis, wkSp0, 1, wkSp0,1); // 
            Vmath::Vadd(nqtots, entropy_residual,1,  wkSp0, 1, entropy_residual, 1);

            Vmath::Vadd(nqtots, velocity[1], 1, old_velocity[1], 1, wkSp1,1);
            Vmath::Smul(nqtots, 0.5, wkSp1, 1, wkSp1,1); // v at (n+1/2)
            m_fields[0]->PhysDeriv(wkSp1, Grad0, Grad1);
            Vmath::Vmul(nqtots, Grad0,1,  Grad0, 1, wkSp0, 1);
            Vmath::Vvtvp(nqtots, Grad1, 1, Grad1, 1, wkSp0, 1, wkSp0, 1); //
            Vmath::Smul(nqtots, m_kinvis, wkSp0, 1, wkSp0,1); // 
            Vmath::Vadd(nqtots, entropy_residual,1,  wkSp0, 1, entropy_residual, 1);
              
            
            //body force
            //wkSp0 wkSp1 wkSp2 are free to use
            Vmath::Vadd(nqtots, velocity[0], 1, old_velocity[0], 1, wkSp0,1);
            Vmath::Smul(nqtots, 0.5, wkSp0, 1, wkSp0,1); // u at (n+1/2)
            Vmath::Vadd(nqtots, velocity[1], 1, old_velocity[1], 1, wkSp1,1);
            Vmath::Smul(nqtots, 0.5, wkSp1, 1, wkSp1,1); // v at (n+1/2)

            Array<OneD, Array<OneD, NekDouble> > velocity_tmp(VelDim);
            velocity_tmp[0] = Grad0;
            velocity_tmp[1] = Grad1;
           
            /*
            std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
            for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
              {
                 Vmath::Zero(nqtots, Grad0, 1);
                 Vmath::Zero(nqtots, Grad1, 1);

                 (*x)->Apply(m_fields, velocity_solution, velocity_tmp, m_time);
                 
                 Vmath::Vmul(nqtots,  wkSp0,1, velocity_tmp[0],1,wkSp2,1);
                 Vmath::Vvtvp(nqtots, wkSp1,1, velocity_tmp[1],1,wkSp2,1, wkSp2,1);

                 Vmath::Vsub(nqtots, entropy_residual,1, wkSp2,1, entropy_residual,1); // subract the forcing 
              }
            */

            break;
          }
        case 3:
          {
            Grad0 = Array<OneD, NekDouble>(nqtots,0.0);
            Grad1 = Array<OneD, NekDouble>(nqtots,0.0);
            Grad2 = Array<OneD, NekDouble>(nqtots,0.0);

            
            //U = 0.5*(u*u)   
            Vmath::Vmul(nqtots,velocity[0], 1, velocity[0], 1, wkSp0, 1);
            Vmath::Vvtvp(nqtots,velocity[1], 1, velocity[1], 1, wkSp0, 1, wkSp0, 1); 
            Vmath::Vvtvp(nqtots,velocity[2], 1, velocity[2], 1, wkSp0, 1, wkSp0, 1); //U
            Vmath::Smul(nqtots, 0.5, wkSp0, 1, wkSp0, 1); // U at (n+1)

            Vmath::Vmul(nqtots, old_velocity[0],1, old_velocity[0],1, wkSp1,1);
            Vmath::Vvtvp(nqtots, old_velocity[1], 1, old_velocity[1], 1, wkSp1, 1, wkSp1, 1); //U
            Vmath::Vvtvp(nqtots, old_velocity[2], 1, old_velocity[2], 1, wkSp1, 1, wkSp1, 1); //U
            Vmath::Smul(nqtots, 0.5, wkSp1, 1, wkSp1,1); // U at (n)

            Vmath::Vsub(nqtots, wkSp0,1, wkSp1,1,entropy_residual,1);
            Vmath::Smul(nqtots,1./m_timestep,entropy_residual,1,entropy_residual,1); // dU/dt
            
            Vmath::Vadd(nqtots, wkSp0, 1, wkSp1, 1, wkSp2,1);
            Vmath::Smul(nqtots, 0.5, wkSp2, 1, wkSp2,1); // U at (n+1/2)

            Vmath::Vadd(nqtots, pressure, 1, old_pressure, 1, wkSp0,1); 
            Vmath::Smul(nqtots, 0.5, wkSp0, 1, wkSp0,1); // P at (n+1/2)
            Vmath::Vadd(nqtots, wkSp0, 1, wkSp2, 1, wkSp0,1); //(p+U)

            //wkSp1 is free to use
            Vmath::Vadd(nqtots, velocity[0], 1, old_velocity[0], 1, wkSp1,1);
            Vmath::Smul(nqtots, 0.5, wkSp1, 1, wkSp1,1); // u at (n+1/2)
            Vmath::Vmul(nqtots, wkSp0,1, wkSp1,1, wkSp1, 1); //(U+p)*u
             
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], wkSp1, Grad0);

            Vmath::Vadd(nqtots, velocity[1], 1, old_velocity[1], 1, wkSp1,1);
            Vmath::Smul(nqtots, 0.5, wkSp1, 1, wkSp1,1); // v at (n+1/2)
            Vmath::Vmul(nqtots, wkSp0,1, wkSp1,1, wkSp1, 1); //(U+p)*v
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], wkSp1, Grad1);

            Vmath::Vadd(nqtots, entropy_residual, 1,  Grad0,1, entropy_residual, 1);
            Vmath::Vadd(nqtots, entropy_residual, 1,  Grad1,1, entropy_residual, 1);

            Vmath::Vadd(nqtots, velocity[2], 1, old_velocity[2], 1, wkSp1,1);
            Vmath::Smul(nqtots, 0.5, wkSp1, 1, wkSp1,1); // w at (n+1/2)
            Vmath::Vmul(nqtots, wkSp0,1, wkSp1,1, wkSp1, 1); //(U+p)*w
            //wkSp0 is free to use
           if ( m_fields[0]->GetWaveSpace() == false )
               m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], wkSp1,Grad2);
           else 
            {
             if( m_fields[0]->GetExpType() == MultiRegions::e3DH1D )
              {
               m_fields[0]->HomogeneousFwdTrans(wkSp1, wkSp0);
               m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], wkSp0, wkSp1);
               m_fields[0]->HomogeneousBwdTrans(wkSp1, Grad2);
              }
             else
              {
                 ASSERTL0(false,"only fully 3D and 3DH1D are implemented ");
              }
            } 

            Vmath::Vadd(nqtots, entropy_residual, 1,  Grad2,1, entropy_residual, 1);

            //wkSp2 still has U at (n+1/2)
            //wkSp0 and wkSp1 are free to use
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], wkSp2, wkSp1);
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], wkSp1, Grad0);

            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], wkSp2, wkSp1);
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], wkSp1, Grad1);

            if( m_fields[0]->GetWaveSpace() == false )
             {
               m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], wkSp2, wkSp1);
               m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], wkSp1, Grad2);
             }
            else 
            {
              if( m_fields[0]->GetExpType() == MultiRegions::e3DH1D )
               {
                m_fields[0]->HomogeneousFwdTrans(wkSp2, wkSp0);
                m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], wkSp0,wkSp1);
                m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], wkSp1,wkSp0);
                m_fields[0]->HomogeneousBwdTrans(wkSp0, Grad2);
               }
             else
               {
                 ASSERTL0(false,"only fully 3D and 3DH1D are implemented ");
               }
            } 

            Vmath::Vadd(nqtots, Grad0, 1,  Grad1, 1, Grad0, 1);
            Vmath::Vadd(nqtots, Grad2, 1,  Grad0, 1, Grad0, 1);
            Vmath::Smul(nqtots,m_kinvis, Grad0, 1, Grad0,1); // laplacian term
            Vmath::Vsub(nqtots, entropy_residual,1,  Grad0, 1, entropy_residual, 1);

            //
            Vmath::Vadd(nqtots, velocity[0], 1, old_velocity[0], 1, wkSp0,1);
            Vmath::Smul(nqtots, 0.5, wkSp0, 1, wkSp0,1); // u at (n+1/2)
            m_fields[0]->PhysDeriv(wkSp0, Grad0, Grad1);
            Vmath::Vmul(nqtots, Grad0,1,  Grad0, 1, wkSp2, 1);
            Vmath::Vvtvp(nqtots, Grad1, 1, Grad1, 1, wkSp2, 1, wkSp2, 1); //
            Vmath::Smul(nqtots, m_kinvis, wkSp2,  1, Grad0,1); // 
            Vmath::Vadd(nqtots, entropy_residual, 1,  Grad0, 1, entropy_residual, 1);

            Vmath::Vadd(nqtots, velocity[1], 1, old_velocity[1], 1, wkSp1,1);
            Vmath::Smul(nqtots, 0.5, wkSp1, 1, wkSp1,1); // v at (n+1/2)
            m_fields[0]->PhysDeriv(wkSp1, Grad0, Grad1);
            Vmath::Vmul(nqtots, Grad0,1,  Grad0, 1, wkSp2, 1);
            Vmath::Vvtvp(nqtots, Grad1, 1, Grad1, 1, wkSp2, 1, wkSp2, 1); //
            Vmath::Smul(nqtots, m_kinvis, wkSp2, 1, Grad0, 1); // 
            Vmath::Vadd(nqtots, entropy_residual, 1,  Grad0, 1, entropy_residual, 1);

            if( m_fields[0]->GetWaveSpace() == false )
             {
              m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], wkSp0, Grad2);
              Vmath::Vmul(nqtots, Grad2, 1, Grad2, 1, wkSp2, 1);
              Vmath::Smul(nqtots, m_kinvis, wkSp2, 1, Grad0, 1); // 
              Vmath::Vadd(nqtots, entropy_residual, 1,  Grad0, 1, entropy_residual, 1);

              m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], wkSp1, Grad2);
              Vmath::Vmul(nqtots, Grad2, 1, Grad2, 1, wkSp2, 1);
              Vmath::Smul(nqtots, m_kinvis, wkSp2, 1, Grad0,1); // 
              Vmath::Vadd(nqtots, entropy_residual,1,  Grad0, 1, entropy_residual, 1);

              Vmath::Vadd(nqtots, velocity[2], 1, old_velocity[2], 1, wkSp2,1);
              Vmath::Smul(nqtots, 0.5, wkSp2, 1, wkSp2,1); // w at (n+1/2)
              m_fields[0]->PhysDeriv(wkSp2, Grad0, Grad1);
              Vmath::Vmul(nqtots, Grad0, 1,  Grad0, 1, wkSp0, 1);
              Vmath::Vvtvp(nqtots, Grad1, 1, Grad1, 1, wkSp0, 1, wkSp0, 1); //
              Vmath::Smul(nqtots, m_kinvis, wkSp0,  1, Grad0,1); // 
              Vmath::Vadd(nqtots, entropy_residual, 1,  Grad0, 1, entropy_residual, 1);
              
              m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], wkSp2, Grad2);
              Vmath::Vmul(nqtots, Grad2, 1, Grad2,  1, wkSp0, 1);
              Vmath::Smul(nqtots, m_kinvis, wkSp0,  1, Grad0,1); // 
              Vmath::Vadd(nqtots, entropy_residual, 1,  Grad0, 1, entropy_residual, 1);
             }
            else
             {
              if( m_fields[0]->GetExpType() == MultiRegions::e3DH1D )
               {
                
                 //u transforms to Fourier space
                 m_fields[0]->HomogeneousFwdTrans(wkSp0, wkSp2);
                 m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], wkSp2, wkSp0);
                 m_fields[0]->HomogeneousBwdTrans(wkSp0, Grad2);
                 Vmath::Vmul(nqtots, Grad2, 1, Grad2, 1, wkSp2, 1);
                 Vmath::Smul(nqtots, m_kinvis, wkSp2, 1, Grad0,1); // 
                 Vmath::Vadd(nqtots, entropy_residual, 1,  Grad0, 1, entropy_residual, 1);

                 //v transforms to Fourier space
                 m_fields[0]->HomogeneousFwdTrans(wkSp1, wkSp2);
                 m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], wkSp2, wkSp0);
                 m_fields[0]->HomogeneousBwdTrans(wkSp0, Grad2);
                 Vmath::Vmul(nqtots, Grad2, 1, Grad2, 1, wkSp2, 1);
                 Vmath::Smul(nqtots, m_kinvis, wkSp2, 1, Grad0,1); // 
                 Vmath::Vadd(nqtots, entropy_residual, 1,  Grad0, 1, entropy_residual, 1);

                 //bellow for w component
                 Vmath::Vadd(nqtots, velocity[2], 1, old_velocity[2], 1, wkSp1,1);
                 Vmath::Smul(nqtots, 0.5, wkSp1, 1, wkSp1,1); // w at (n+1/2)
                 m_fields[0]->PhysDeriv(wkSp1, Grad0, Grad1);

                 Vmath::Vmul(nqtots, Grad0, 1,  Grad0, 1, wkSp0, 1);
                 Vmath::Vvtvp(nqtots, Grad1, 1, Grad1, 1, wkSp0, 1, wkSp0, 1); //
                 Vmath::Smul(nqtots, m_kinvis, wkSp0,  1, Grad0,1); // 
                 Vmath::Vadd(nqtots, entropy_residual, 1,  Grad0, 1, entropy_residual, 1);

                 m_fields[0]->HomogeneousFwdTrans(wkSp1, wkSp2);
                 m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], wkSp2, wkSp0);
                 m_fields[0]->HomogeneousBwdTrans(wkSp0, Grad2);
                 //d/dz transformed into (Q,P)-space
                 Vmath::Vmul(nqtots, Grad2, 1, Grad2,  1, wkSp0, 1);
                 Vmath::Smul(nqtots, m_kinvis, wkSp0,  1, Grad0,1); // 
                 Vmath::Vadd(nqtots, entropy_residual, 1,  Grad0, 1, entropy_residual, 1);
            
               }
              else
               {
                 ASSERTL0(false,"only fully 3D and 3DH1D are implemented ");
               }
             }
              
              //body force
                //wkSp0 wkSp1 wkSp2 are free to use
                 Vmath::Vadd(nqtots, velocity[0], 1, old_velocity[0], 1, wkSp0,1);
                 Vmath::Smul(nqtots, 0.5, wkSp0, 1, wkSp0,1); // u at (n+1/2)

                 Vmath::Vadd(nqtots, velocity[1], 1, old_velocity[1], 1, wkSp1,1);
                 Vmath::Smul(nqtots, 0.5, wkSp1, 1, wkSp1,1); // v at (n+1/2)

                 Vmath::Vadd(nqtots, velocity[2], 1, old_velocity[2], 1, wkSp2,1);
                 Vmath::Smul(nqtots, 0.5, wkSp2, 1, wkSp2,1); // w at (n+1/2)

                 Array<OneD, Array<OneD, NekDouble> > velocity_tmp(VelDim);
                 velocity_tmp[0] = Grad0;
                 velocity_tmp[1] = Grad1;
                 velocity_tmp[2] = Grad2;
               
                 /*
              
                 std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
                 for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
                  {
                     for(int d=0; d<VelDim; ++d)
                         Vmath::Zero(nqtots, velocity_tmp[d], 1);

                        (*x)->Apply(m_fields, velocity, velocity_tmp, m_time);
                 
                          Vmath::Vmul(nqtots,  Grad0,1, wkSp0,1,Grad0,1);
                          Vmath::Vvtvp(nqtots, Grad1,1, wkSp1,1,Grad0,1,Grad0,1);
                          Vmath::Vvtvp(nqtots, Grad2,1, wkSp2,1,Grad0,1,Grad0,1);

                    Vmath::Vsub(nqtots, entropy_residual,1, Grad0,1, entropy_residual,1); // subract the forcing 
                  }
                */
                
                 if(m_use_evm_mapping)
                  {
                      Vmath::Vmul(nqtots,  m_evm_mapping_force[0],1, wkSp0,1,Grad0,1);
                      Vmath::Vvtvp(nqtots, m_evm_mapping_force[1],1, wkSp1,1,Grad0,1,Grad0,1);
                      Vmath::Vvtvp(nqtots, m_evm_mapping_force[2],1, wkSp2,1,Grad0,1,Grad0,1);
                      
                      Vmath::Vsub(nqtots, entropy_residual,1, Grad0,1, entropy_residual,1); // subract the forcing 
                  }
            break;
          }
         default:
            ASSERTL0(false,"dimension unknown");
        }

              

        //const std::pair<NekDouble, NekDouble> & velocity_range = EvaluateVelocityRange(velocity);
        const NekDouble VelocityEntropyVariation = EvaluateVelocityEntropyVariation(velocity);

//        if ( m_session->GetComm()->GetRank() == 0 )
//           cout<<"entropy = "<<VelocityEntropyVariation<<endl;

  //      if ( (m_session->GetComm()->GetRank() == 0 )
  //         && ( ! ((m_step_counter + 1)%m_infosteps )) )
  //           cout<<"Maximal Velocity    : "<<setw(8)<<left<<velocity_range.second<<endl;

            LibUtilities::PointsType stdPointsType
              = m_fields[0]->GetExp(0)->GetPointsType(0);
            int nQuadPoints = m_fields[0]->GetExp(0)->GetNumPoints(0);
            const LibUtilities::PointsKey quadPointsKey(nQuadPoints, stdPointsType);
                  Array<OneD, NekDouble> quadZeros(nQuadPoints);
            quadZeros   = (LibUtilities::PointsManager()[quadPointsKey])->GetZ();


        Array<OneD, NekDouble> hOverP;
        hOverP = Array<OneD, NekDouble>(nElmts,0.0);

        switch(ndim)
        {
         case 1:
            ASSERTL0(false,"EVM not implemented for one dimensional problem");
            break;
         case 2:
           {
              wkCoef0 = Array<OneD, NekDouble>(nctots,0.0);
              wkCoef1 = Array<OneD, NekDouble>(nctots,0.0);
               
              Vmath::Smul(nqtots, BBeta[1][0], velocity[0], 1, wkSp0,1);
              Vmath::Svtvp(nqtots, BBeta[1][1], old_velocity[0], 1, wkSp0,1, wkSp0, 1);
              
              Vmath::Smul(nqtots, BBeta[1][0], velocity[1], 1, wkSp1,1);
              Vmath::Svtvp(nqtots, BBeta[1][1], old_velocity[1], 1, wkSp1,1, wkSp1, 1);
              
              Vmath::Vmul(nqtots,wkSp0,1, wkSp0,1, wkSp2,1);
              Vmath::Vvtvp(nqtots, wkSp1, 1, wkSp1, 1, wkSp2, 1, wkSp2, 1); 
              Vmath::Smul(nqtots, 0.5, wkSp2, 1, wkSp2,1); //0.5*(u^2+v^2)
                 
              max_viscosity = 0.0;
              max_alpha = 0.0;
              max_velocity = 0.0;
             for(int e=0; e<nElmts; ++e)
              {
                 int offset = m_fields[0]->GetPhys_Offset(e);
                 Array<OneD, NekDouble> tmp;
                 NekDouble elmtSensor = m_evm_alpha;
                if(m_use_dynamic_alpha) 
                 {
                  int numModesElement = expOrderElement[e];
                  int nElmtCoeffs     = m_fields[0]->GetExp(e)->GetNcoeffs();
                  int numCutOff       = numModesElement - evm_reduced_order;
                
                  Array<OneD, NekDouble> elmtPhys(nqElmts);
                  Vmath::Vcopy(nqElmts, tmp = wkSp2+offset, 1, elmtPhys, 1);

                  // Compute coefficients
                  Array<OneD, NekDouble> elmtCoeffs(nElmtCoeffs, 0.0);
                  m_fields[0]->GetExp(e)->FwdTrans(elmtPhys, elmtCoeffs);

                  Array<OneD, NekDouble> reducedElmtCoeffs(nElmtCoeffs, 0.0);
                  m_fields[0]->GetExp(e)->ReduceOrderCoeffs(numCutOff, elmtCoeffs,
                                                         reducedElmtCoeffs);

                  Array<OneD, NekDouble> reducedElmtPhys(nqElmts, 0.0);
                  m_fields[0]->GetExp(e)->BwdTrans(reducedElmtCoeffs, reducedElmtPhys);
               

                  NekDouble numerator   = 0.0;
                  NekDouble denominator = 0.0;

                // Determining the norm of the numerator of the Sensor
                  Array<OneD, NekDouble> difference(nqElmts, 0.0);
                  Vmath::Vsub(nqElmts, elmtPhys, 1, reducedElmtPhys, 1, difference, 1);

                  numerator = Vmath::Dot(nqElmts, difference, difference);
                  denominator = Vmath::Dot(nqElmts, elmtPhys, elmtPhys);

                  elmtSensor = sqrt(numerator / denominator);
   //              elmtSensor = log10(max(elmtSensor, NekConstants::kNekSqrtTol));
                  max_alpha = max(max_alpha, elmtSensor);
                 
                  Vmath::Fill(nqElmts,elmtSensor, tmp = m_evm_sensor+offset,1);
                }

                NekDouble h = 1.0e+10;
                LocalRegions::Expansion2DSharedPtr exp2D;
                exp2D = m_fields[0]->GetExp(e)->as<LocalRegions::Expansion2D>();
                for(int i = 0; i < exp2D->GetNedges(); ++i)
                {
                    h = min(h, exp2D->GetGeom2D()->GetEdge(i)->GetVertex(0)->
                        dist(*(exp2D->GetGeom2D()->GetEdge(i)->GetVertex(1))));
                }

                h *= (0.5*(quadZeros[1]-quadZeros[0]));
              
                
                NekDouble first_order_viscosity = 0.0;
                for(int q=0; q<nqElmts; ++q)
                  {
                    NekDouble velocity_magnitude = 0.;
                   for(int d=0; d<VelDim; ++d)
                      velocity_magnitude += velocity[d][offset+q]*velocity[d][offset+q];

                     first_order_viscosity = max(first_order_viscosity, m_evm_beta*h*sqrt(velocity_magnitude));

                     max_velocity = max(max_velocity,sqrt(velocity_magnitude));
                  }

                NekDouble evm_alpha = m_evm_alpha;
                if(m_use_dynamic_alpha)
                  evm_alpha = elmtSensor;

                NekDouble elmt_entropy_viscosity = 0.0;
                for(int q=0; q<nqElmts; ++q)
                 {
                   NekDouble absResidual = fabs(entropy_residual[offset+q]);
                   
                   elmt_entropy_viscosity = max(elmt_entropy_viscosity,
//                                  m_evm_alpha*h*h*absResidual/VelocityEntropyVariation);
                                    evm_alpha*h*h*absResidual/VelocityEntropyVariation);
                   
                 //  if(elmt_entropy_viscosity>0.1)
                 //  cout<<"elmt_entropy_viscosity = "<<elmt_entropy_viscosity<<" absResidual = "<<absResidual<<" entropy variation = "
                 //     <<VelocityEntropyVariation<<" first order viscosity = "<<first_order_viscosity<<" h = "<<h<<endl;
                 }
                  
                   elmt_entropy_viscosity = min(first_order_viscosity, elmt_entropy_viscosity);

                   Vmath::Fill(nqElmts,elmt_entropy_viscosity, tmp = m_evm_visc+offset,1);

                   max_viscosity = max(max_viscosity, elmt_entropy_viscosity);

                   elmt_entropy_viscosity /= m_kinvis;

                   Array<OneD, NekDouble> tmp0,tmp1;
                   Vmath::Smul(nqElmts, elmt_entropy_viscosity,  tmp0 = wkSp0 + offset, 1, tmp1 = wkSp0 + offset, 1);
                   Vmath::Smul(nqElmts, elmt_entropy_viscosity,  tmp0 = wkSp1 + offset, 1, tmp1 = wkSp1 + offset, 1);
               }
                   
                   
             //       bool waveSpace = m_fields[0]->GetWaveSpace();
             //       m_fields[0]->SetWaveSpace(true);

                    m_fields[0]->PhysDeriv(wkSp0, Grad0, Grad1);

                    m_fields[0]->IProductWRTDerivBase(0, Grad0, wkCoef0);
                    m_fields[0]->IProductWRTDerivBase(1, Grad1, wkCoef1);
                    Vmath::Vadd(nctots, wkCoef0, 1,  wkCoef1, 1, wkCoef0, 1);

                    m_fields[0]->SetPhysState (false);
                    m_fields[0]->MultiplyByElmtInvMass(wkCoef0, wkCoef0);
                    m_fields[0]->BwdTrans(wkCoef0, Grad0);
                    Vmath::Vadd(nqtots, outarray[0], 1,  Grad0, 1, outarray[0], 1);
                    
                   
//                    m_fields[0]->SetWaveSpace(true);
                    m_fields[0]->PhysDeriv(wkSp1, Grad0, Grad1);
                    m_fields[0]->IProductWRTDerivBase(0, Grad0, wkCoef0);
                    m_fields[0]->IProductWRTDerivBase(1, Grad1, wkCoef1);
                    Vmath::Vadd(nctots, wkCoef0, 1,  wkCoef1, 1, wkCoef0, 1);

                    m_fields[0]->SetPhysState (false);
                    m_fields[0]->MultiplyByElmtInvMass(wkCoef0, wkCoef1);
                    m_fields[0]->BwdTrans(wkCoef1, Grad1);
                    Vmath::Vadd(nqtots, outarray[1], 1,  Grad1, 1, outarray[1], 1);

           //         m_fields[0]->SetWaveSpace(waveSpace);

             break;
          }
         case 3:
          {
              wkCoef0 = Array<OneD, NekDouble>(nctots,0.0);
              wkCoef1 = Array<OneD, NekDouble>(nctots,0.0);
              wkCoef2 = Array<OneD, NekDouble>(nctots,0.0);

              Vmath::Smul(nqtots, BBeta[1][0], velocity[0], 1, wkSp0,1);
              Vmath::Svtvp(nqtots, BBeta[1][1], old_velocity[0], 1, wkSp0,1, wkSp0, 1);
              
              Vmath::Smul(nqtots, BBeta[1][0], velocity[1], 1, wkSp1,1);
              Vmath::Svtvp(nqtots, BBeta[1][1], old_velocity[1], 1, wkSp1,1, wkSp1, 1);

              Vmath::Smul(nqtots, BBeta[1][0], velocity[2], 1, wkSp2,1);
              Vmath::Svtvp(nqtots, BBeta[1][1], old_velocity[2], 1, wkSp2,1, wkSp2, 1);
            
              Vmath::Vmul(nqtots,wkSp0,1, wkSp0,1, Grad0,1);
              Vmath::Vvtvp(nqtots, wkSp1, 1, wkSp1, 1, Grad0, 1, Grad0, 1); 
              Vmath::Vvtvp(nqtots, wkSp2, 1, wkSp2, 1, Grad0, 1, Grad0, 1); 
              Vmath::Smul(nqtots, 0.5, Grad0, 1, Grad0,1); //0.5*(u^2+v^2)

            if( m_fields[0]->GetWaveSpace() == false )
             {
               max_viscosity = 0.0;
               max_alpha = 0.0;
               max_velocity = 0.0;
              for(int e=0; e<nElmts; ++e)
               {

                Array<OneD, NekDouble> tmp;
                int offset = m_fields[0]->GetPhys_Offset(e);
                NekDouble elmtSensor = m_evm_alpha;
                if(m_use_dynamic_alpha) 
                 {
                  int numModesElement = expOrderElement[e];
                  int nElmtCoeffs     = m_fields[0]->GetExp(e)->GetNcoeffs();
                  int numCutOff       = numModesElement - evm_reduced_order;
                
                  Array<OneD, NekDouble> elmtPhys(nqElmts);
                  Vmath::Vcopy(nqElmts, tmp = Grad0+offset, 1, elmtPhys, 1);

                  // Compute coefficients
                  Array<OneD, NekDouble> elmtCoeffs(nElmtCoeffs, 0.0);
                  m_fields[0]->GetExp(e)->FwdTrans(elmtPhys, elmtCoeffs);

                  Array<OneD, NekDouble> reducedElmtCoeffs(nElmtCoeffs, 0.0);
                  m_fields[0]->GetExp(e)->ReduceOrderCoeffs(numCutOff, elmtCoeffs,
                                                         reducedElmtCoeffs);

                  Array<OneD, NekDouble> reducedElmtPhys(nqElmts, 0.0);
                  m_fields[0]->GetExp(e)->BwdTrans(reducedElmtCoeffs, reducedElmtPhys);
               

                  NekDouble numerator   = 0.0;
                  NekDouble denominator = 0.0;

                  // Determining the norm of the numerator of the Sensor
                  Array<OneD, NekDouble> difference(nqElmts, 0.0);
                  Vmath::Vsub(nqElmts, elmtPhys, 1, reducedElmtPhys, 1, difference, 1);

                  numerator = Vmath::Dot(nqElmts, difference, difference);
                  denominator = Vmath::Dot(nqElmts, elmtPhys, elmtPhys);

                  NekDouble elmtSensor = sqrt(numerator / denominator);
   //              elmtSensor = log10(max(elmtSensor, NekConstants::kNekSqrtTol));
                  max_alpha = max(max_alpha, elmtSensor);
                 
                  Vmath::Fill(nqElmts,elmtSensor, tmp = m_evm_sensor+offset,1);
                }

                  NekDouble h = 1.0e+10;
                  LocalRegions::Expansion3DSharedPtr exp3D;
                  exp3D = m_fields[0]->GetExp(e)->as<LocalRegions::Expansion3D>();
                for(int i = 0; i < exp3D->GetNedges(); ++i)
                {
                    h = min(h, exp3D->GetGeom3D()->GetEdge(i)->GetVertex(0)->
                        dist(*(exp3D->GetGeom3D()->GetEdge(i)->GetVertex(1))));
                }

                h *= (0.5*(quadZeros[1]-quadZeros[0]));
                
//                int offset = m_fields[0]->GetPhys_Offset(e);
                NekDouble first_order_viscosity = 0.0;
                for(int q=0; q<nqElmts; ++q)
                  {
                    NekDouble velocity_magnitude = 0.;
                   for(int d=0; d<VelDim; ++d)
                      velocity_magnitude += velocity[d][offset+q]*velocity[d][offset+q];

                     first_order_viscosity = max(first_order_viscosity, m_evm_beta*h*sqrt(velocity_magnitude));
                     
                     max_velocity = max(max_velocity,sqrt(velocity_magnitude));
                  }

                NekDouble evm_alpha = m_evm_alpha;
                if(m_use_dynamic_alpha)
                  evm_alpha = elmtSensor;

                NekDouble elmt_entropy_viscosity = 0.0;
                for(int q=0; q<nqElmts; ++q)
                 {
                   NekDouble absResidual = fabs(entropy_residual[offset+q]);
                   
                   elmt_entropy_viscosity = max(elmt_entropy_viscosity,
                                  evm_alpha*h*h*absResidual/VelocityEntropyVariation);
                 }
                  
                   elmt_entropy_viscosity = min(first_order_viscosity, elmt_entropy_viscosity);
                   
                   Vmath::Fill(nqElmts,elmt_entropy_viscosity, tmp = m_evm_visc+offset,1);

                   max_viscosity = max(max_viscosity, elmt_entropy_viscosity);

                   elmt_entropy_viscosity /= m_kinvis;

                   Array<OneD, NekDouble> tmp0,tmp1;
                   Vmath::Smul(nqElmts, elmt_entropy_viscosity,  tmp0 = wkSp0 + offset, 1, tmp1 = wkSp0 + offset, 1);
                   Vmath::Smul(nqElmts, elmt_entropy_viscosity,  tmp0 = wkSp1 + offset, 1, tmp1 = wkSp1 + offset, 1);
                   Vmath::Smul(nqElmts, elmt_entropy_viscosity,  tmp0 = wkSp2 + offset, 1, tmp1 = wkSp2 + offset, 1);
               }
                    //Grad0,Grad1,Grad2 are safe to use
                    m_fields[0]->PhysDeriv(wkSp0, Grad0, Grad1, Grad2);

                    m_fields[0]->IProductWRTDerivBase(0, Grad0, wkCoef0);
                    m_fields[0]->IProductWRTDerivBase(1, Grad1, wkCoef1);
                    m_fields[0]->IProductWRTDerivBase(2, Grad2, wkCoef2);
                    Vmath::Vadd(nctots, wkCoef0, 1,  wkCoef1, 1, wkCoef0, 1);
                    Vmath::Vadd(nctots, wkCoef0, 1,  wkCoef2, 1, wkCoef0, 1);

                    m_fields[0]->SetPhysState (false);
                    m_fields[0]->MultiplyByElmtInvMass(wkCoef0, wkCoef0);
                    m_fields[0]->BwdTrans(wkCoef0, Grad0);
                    Vmath::Vadd(nqtots, outarray[0], 1,  Grad0, 1, outarray[0], 1);
                    
                   
                    m_fields[0]->PhysDeriv(wkSp1, Grad0, Grad1, Grad2);
                    m_fields[0]->IProductWRTDerivBase(0, Grad0, wkCoef0);
                    m_fields[0]->IProductWRTDerivBase(1, Grad1, wkCoef1);
                    m_fields[0]->IProductWRTDerivBase(2, Grad2, wkCoef2);
                    Vmath::Vadd(nctots, wkCoef0, 1,  wkCoef1, 1, wkCoef0, 1);
                    Vmath::Vadd(nctots, wkCoef0, 1,  wkCoef2, 1, wkCoef0, 1);

                    m_fields[0]->SetPhysState (false);
                    m_fields[0]->MultiplyByElmtInvMass(wkCoef0, wkCoef1);
                    m_fields[0]->BwdTrans(wkCoef1, Grad1);
                    Vmath::Vadd(nqtots, outarray[1], 1,  Grad1, 1, outarray[1], 1);

                   // m_fields[0]->SetWaveSpace(true);
                    m_fields[0]->PhysDeriv(wkSp2, Grad0, Grad1, Grad2);
                    m_fields[0]->IProductWRTDerivBase(0, Grad0, wkCoef0);
                    m_fields[0]->IProductWRTDerivBase(1, Grad1, wkCoef1);
                    m_fields[0]->IProductWRTDerivBase(2, Grad2, wkCoef2);
                    Vmath::Vadd(nctots, wkCoef0, 1,  wkCoef1, 1, wkCoef0, 1);
                    Vmath::Vadd(nctots, wkCoef0, 1,  wkCoef2, 1, wkCoef0, 1);

                    m_fields[0]->SetPhysState (false);
                    m_fields[0]->MultiplyByElmtInvMass(wkCoef0, wkCoef1);
                    m_fields[0]->BwdTrans(wkCoef1, Grad1);
                    Vmath::Vadd(nqtots, outarray[2], 1,  Grad1, 1, outarray[2], 1);


                }//fully 3D
            else
            {
             
              if( m_fields[0]->GetExpType() == MultiRegions::e3DH1D )
               {
              

                Array<OneD, NekDouble> visc_tmp, visc_ave;

                visc_tmp = Array<OneD, NekDouble>(nElmts*NZ,0.0);
                visc_ave = Array<OneD, NekDouble>(nElmts,0.0);

                int cnt = 0;
                max_velocity = 0.0;
               for(int k=0; k<NZ; ++k)
                {

                 for(int e=0; e<nElmts; ++e, ++cnt)
                  {
                    Array<OneD, NekDouble> tmp;
                    int offset = m_fields[0]->GetPhys_Offset(cnt);
                    NekDouble elmtSensor = m_evm_alpha;
                   if(m_use_dynamic_alpha) 
                    {
                      int numModesElement = expOrderElement[e];
                      int nElmtCoeffs     = m_fields[0]->GetPlane(0)->GetExp(e)->GetNcoeffs();
                      int numCutOff       = numModesElement - evm_reduced_order;

                      Array<OneD, NekDouble> elmtPhys(nqElmts);
                      Vmath::Vcopy(nqElmts, tmp = Grad0+offset, 1, elmtPhys, 1);

                     // Compute coefficients
                      Array<OneD, NekDouble> elmtCoeffs(nElmtCoeffs, 0.0);
                      m_fields[0]->GetPlane(0)->GetExp(e)->FwdTrans(elmtPhys, elmtCoeffs);

                      Array<OneD, NekDouble> reducedElmtCoeffs(nElmtCoeffs, 0.0);
                      m_fields[0]->GetPlane(0)->GetExp(e)->ReduceOrderCoeffs(numCutOff, elmtCoeffs,
                                                         reducedElmtCoeffs);

                      Array<OneD, NekDouble> reducedElmtPhys(nqElmts, 0.0);
                      m_fields[0]->GetPlane(0)->GetExp(e)->BwdTrans(reducedElmtCoeffs, reducedElmtPhys);

                      NekDouble numerator   = 0.0;
                      NekDouble denominator = 0.0;

                     // Determining the norm of the numerator of the Sensor
                      Array<OneD, NekDouble> difference(nqElmts, 0.0);
                      Vmath::Vsub(nqElmts, elmtPhys, 1, reducedElmtPhys, 1, difference, 1);

                      numerator = Vmath::Dot(nqElmts, difference, difference);
                      denominator = Vmath::Dot(nqElmts, elmtPhys, elmtPhys);

                      NekDouble elmtSensor = sqrt(numerator / denominator);
                      max_alpha = max(max_alpha, elmtSensor);
                 
                      Vmath::Fill(nqElmts,elmtSensor, tmp = m_evm_sensor+offset,1);
                   }

                    NekDouble h = 1.0e+10;
                    LocalRegions::Expansion2DSharedPtr exp2D;
                    exp2D = m_fields[0]->GetPlane(0)->GetExp(e)->as<LocalRegions::Expansion2D>();
                 
                   for(int i = 0; i < exp2D->GetNedges(); ++i)
                    {
                       h = min(h, exp2D->GetGeom2D()->GetEdge(i)->GetVertex(0)->
                           dist(*(exp2D->GetGeom2D()->GetEdge(i)->GetVertex(1))));
                    }

                    h *= (0.5*(quadZeros[1]-quadZeros[0]));
                
                    NekDouble first_order_viscosity = 0.0;
                   for(int q=0; q<nqElmts; ++q)
                     {
                         NekDouble velocity_magnitude = 0.;
                       for(int d=0; d<VelDim; ++d)
                         velocity_magnitude += velocity[d][offset+q]*velocity[d][offset+q];

                          first_order_viscosity = max(first_order_viscosity, m_evm_beta*h*sqrt(velocity_magnitude));
                          
                          max_velocity = max(max_velocity, sqrt(velocity_magnitude));
                     }

                    NekDouble evm_alpha = m_evm_alpha;
                    if(m_use_dynamic_alpha)
                       evm_alpha = elmtSensor;

                   NekDouble elmt_entropy_viscosity = 0.0;
                  for(int q=0; q<nqElmts; ++q)
                    {
                      NekDouble absResidual = fabs(entropy_residual[offset+q]);
                   
                      elmt_entropy_viscosity = max(elmt_entropy_viscosity,
                                  evm_alpha*h*h*absResidual/VelocityEntropyVariation);
//                                  elmtSensor*h*h*absResidual/VelocityEntropyVariation);
//                       if(elmt_entropy_viscosity>0.1)
//                      cout<<" elmt_entropy_viscosity = "<<elmt_entropy_viscosity<<" residual = "<<absResidual<<" entropy_variation = "
//                        <<VelocityEntropyVariation<<" h = "<<h<<" sensor = "<<elmtSensor<<endl;
                    }
                  
                      elmt_entropy_viscosity = min(first_order_viscosity, elmt_entropy_viscosity);

                      Vmath::Fill(nqElmts,elmt_entropy_viscosity, tmp = m_evm_visc+offset,1);

                      elmt_entropy_viscosity /= m_kinvis;
                      
                      visc_tmp[k*nElmts+e] = elmt_entropy_viscosity;
                  }
               }//end of NZ planes
                   

                     for(int e=0; e<nElmts; ++e)
                       visc_ave[e] = visc_tmp[e];
                       
                     for(int k=1; k<NZ; ++k)
                      for(int e=0; e<nElmts; ++e)
                        visc_ave[e] = max(visc_ave[e], visc_tmp[k*nElmts+e]);
//

                     for(int e=0; e<nElmts; ++e)
                        m_comm->GetColumnComm()->AllReduce(visc_ave[e], LibUtilities::ReduceMax);

                      max_viscosity = 0.0;
                      for(int e=0; e<nElmts; ++e) max_viscosity = max(max_viscosity, visc_ave[e]);

                      max_viscosity *=m_kinvis;

                  // if (m_session->GetComm()->GetRank() == 0 )
                  //   for(int e=0; e<nElmts; ++e)
                  //   cout<<"visc["<<e<<"] = "<<visc_ave[e]*m_kinvis<<endl;

                     cnt = 0;
                    for(int k=0; k<NZ; ++k)
                    for(int e=0; e<nElmts; ++e, ++cnt)
                     {
                       int offset = m_fields[0]->GetPhys_Offset(cnt);

                       Array<OneD, NekDouble> tmp0,tmp1;
                       Vmath::Smul(nqElmts, visc_ave[e],  tmp0 = wkSp0 + offset, 1, tmp1 = wkSp0 + offset, 1);
                       Vmath::Smul(nqElmts, visc_ave[e],  tmp0 = wkSp1 + offset, 1, tmp1 = wkSp1 + offset, 1);
                       Vmath::Smul(nqElmts, visc_ave[e],  tmp0 = wkSp2 + offset, 1, tmp1 = wkSp2 + offset, 1);
                     }


                    m_fields[0]->HomogeneousFwdTrans(wkSp0, Grad2);
                    m_fields[0]->PhysDeriv(Grad2, Grad0, Grad1);

                    //u transforms to Fourier space
                    for(int k=0; k<NZ; ++k)
//                    if(k != 1 || m_trans->GetK(k) != 0)
                     {
                      Array<OneD, NekDouble> tmp0,tmp1;
                      m_fields[0]->GetPlane(k)->IProductWRTDerivBase(0, tmp0 = Grad0+k*nqPtots, tmp1 = wkCoef0+k*ncPtots);
                      m_fields[0]->GetPlane(k)->IProductWRTDerivBase(1, tmp0 = Grad1+k*nqPtots, tmp1 = wkCoef1+k*ncPtots);
                     }

                    Vmath::Vadd(nctots, wkCoef0, 1,  wkCoef1, 1, wkCoef0, 1);

                    for(int k=0; k<NZ; ++k)
//                    if(k != 1 || m_trans->GetK(k) != 0)
                     {
                      Array<OneD, NekDouble> tmp0,tmp1;
                      bool physState = m_fields[0]->GetPhysState();
                      m_fields[0]->GetPlane(k)->SetPhysState (false);
                      m_fields[0]->GetPlane(k)->MultiplyByElmtInvMass(tmp0 = wkCoef0+k*ncPtots, tmp1 =wkCoef0+k*ncPtots);
                      m_fields[0]->GetPlane(k)->SetPhysState (physState);
                     }

                     m_fields[0]->BwdTrans(wkCoef0, Grad0);

                     Vmath::Vadd(nqtots, outarray[0], 1,  Grad0, 1, outarray[0], 1);
                    
                     for (int k = 0; k <NZ; ++k)
//                    if(k != 1 || m_trans->GetK(k) != 0)
                       {
                          Array<OneD, NekDouble> tmp1,tmp2; 
                          double beta  = 2*M_PI*m_trans->GetK(k)/m_homoLen;
                          beta *= beta;
                          Vmath::Smul(nqPtots, beta, tmp1 = Grad2 + k*nqPtots, 1, tmp2 = wkSp0 + k*nqPtots, 1);
                       }
                       Vmath::Vadd(nqtots, outarray[0], 1, wkSp0, 1, outarray[0], 1);
                      
                   //
                   //
                   //
                    //v transforms to Fourier space
                    m_fields[0]->HomogeneousFwdTrans(wkSp1, Grad2);
                    m_fields[0]->PhysDeriv(Grad2, Grad0, Grad1);

                    for(int k=0; k<NZ; ++k)
//                    if(k != 1 || m_trans->GetK(k) != 0)
                     {
                      Array<OneD, NekDouble> tmp0,tmp1;
                      m_fields[0]->GetPlane(k)->IProductWRTDerivBase(0, tmp0 = Grad0+k*nqPtots, tmp1 = wkCoef0+k*ncPtots);
                      m_fields[0]->GetPlane(k)->IProductWRTDerivBase(1, tmp0 = Grad1+k*nqPtots, tmp1 = wkCoef1+k*ncPtots);
                     }

                    Vmath::Vadd(nctots, wkCoef0, 1,  wkCoef1, 1, wkCoef0, 1);

          //          m_fields[0]->SetPhysState (false);
                    for(int k=0; k<NZ; ++k)
//                    if(k != 1 || m_trans->GetK(k) != 0)
                     {
                      Array<OneD, NekDouble> tmp0,tmp1;
                      bool physState = m_fields[0]->GetPhysState();
                      m_fields[0]->GetPlane(k)->SetPhysState (false);
                      m_fields[0]->GetPlane(k)->MultiplyByElmtInvMass(tmp0 = wkCoef0+k*ncPtots, tmp1 =wkCoef0+k*ncPtots);
                      m_fields[0]->GetPlane(k)->SetPhysState (physState);
                     }

                     m_fields[0]->BwdTrans(wkCoef0, Grad0);
                     
                     Vmath::Vadd(nqtots, outarray[1], 1,  Grad0, 1, outarray[1], 1);
                   //
                     for (int k = 0; k <NZ; ++k)
//                     if(k != 1 || m_trans->GetK(k) != 0)
                       {
                          Array<OneD, NekDouble> tmp1,tmp2; 
                          double beta  = 2*M_PI*m_trans->GetK(k)/m_homoLen;
                          beta *= beta;
                          Vmath::Smul(nqPtots, beta, tmp1 = Grad2 + k*nqPtots, 1, tmp2 = wkSp0 + k*nqPtots, 1);
                       }
                        Vmath::Vadd(nqtots, outarray[1], 1, wkSp0, 1, outarray[1], 1);
                   //
                   //
                   //
                     //w transforms to Fourier space
                    m_fields[0]->HomogeneousFwdTrans(wkSp2, Grad2);
                    m_fields[0]->PhysDeriv(Grad2, Grad0, Grad1);

                    for(int k=0; k<NZ; ++k)
//                    if(k != 1 || m_trans->GetK(k) != 0)
                     {
                      Array<OneD, NekDouble> tmp0,tmp1;
                      m_fields[0]->GetPlane(k)->IProductWRTDerivBase(0, tmp0 = Grad0+k*nqPtots, tmp1 = wkCoef0+k*ncPtots);
                      m_fields[0]->GetPlane(k)->IProductWRTDerivBase(1, tmp0 = Grad1+k*nqPtots, tmp1 = wkCoef1+k*ncPtots);
                     }

                    Vmath::Vadd(nctots, wkCoef0, 1,  wkCoef1, 1, wkCoef0, 1);

                   // m_fields[0]->SetPhysState (false);
                    for(int k=0; k<NZ; ++k)
//                    if(k != 1 || m_trans->GetK(k) != 0)
                     {
                      Array<OneD, NekDouble> tmp0,tmp1;
                      bool physState = m_fields[0]->GetPhysState();
                      m_fields[0]->GetPlane(k)->SetPhysState (false);
                      m_fields[0]->GetPlane(k)->MultiplyByElmtInvMass(tmp0 = wkCoef0+k*ncPtots, tmp1 =wkCoef0+k*ncPtots);
                      m_fields[0]->GetPlane(k)->SetPhysState (physState);
                     }

                     m_fields[0]->BwdTrans(wkCoef0, Grad0);
                     
                     Vmath::Vadd(nqtots, outarray[2], 1,  Grad0, 1, outarray[2], 1);
                   
                     for (int k = 0; k <NZ; ++k)
 //                      if(k != 1 || m_trans->GetK(k) != 0)
                       {
                          Array<OneD, NekDouble> tmp1,tmp2; 
                          double beta  = 2*M_PI*m_trans->GetK(k)/m_homoLen;
                          beta *= beta;
                          Vmath::Smul(nqPtots, beta, tmp1 = Grad2 + k*nqPtots, 1, tmp2 = wkSp0 + k*nqPtots, 1);
                       }
                      Vmath::Vadd(nqtots, outarray[2], 1, wkSp0, 1, outarray[2], 1);


               }
              else
               {
                 ASSERTL0(false,"only fully 3D and 3DH1D are implemented ");
               }
             }

            break;
          }
         default:
            ASSERTL0(false,"dimension unknown");
        }

             
               if( m_fields[0]->GetExpType() == MultiRegions::e3DH1D )
                {
                   m_comm->GetColumnComm()->AllReduce(max_velocity, LibUtilities::ReduceMax);
                   m_comm->GetColumnComm()->AllReduce(max_viscosity, LibUtilities::ReduceMax);
                   m_comm->GetColumnComm()->AllReduce(max_alpha, LibUtilities::ReduceMax);
                   m_comm->GetRowComm()->AllReduce(max_velocity, LibUtilities::ReduceMax);
                   m_comm->GetRowComm()->AllReduce(max_viscosity, LibUtilities::ReduceMax);
                   m_comm->GetRowComm()->AllReduce(max_alpha, LibUtilities::ReduceMax);
                }
               else
                {
                  m_comm->AllReduce(max_viscosity,LibUtilities::ReduceMax);
                  m_comm->AllReduce(max_alpha,LibUtilities::ReduceMax);
                  m_comm->AllReduce(max_velocity,LibUtilities::ReduceMax);
                }
              
              
             if(!m_use_dynamic_alpha) max_alpha = m_evm_alpha;

             if ( (m_session->GetComm()->GetRank() == 0 )
                 && ( ! (m_step_counter%m_infosteps )) )
                 cout<<" Max Velocity  : "<<setw(12)<<left<<max_velocity<<"Max EVM Viscosity   : "<<setw(12)<<left
                       <<max_viscosity<<" Max EVM alpha  :"<<setw(12)<<left<<max_alpha<<endl;

    }


    void VelocityCorrectionScheme::
      v_ExtraFldOutput(
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
          if(m_useEntropyViscosity)
           {

            Array<OneD, NekDouble> evm_visc_Fwd(nCoeffs);

            m_fields[0]->FwdTrans_IterPerExp(m_evm_visc,   evm_visc_Fwd);
            variables.push_back  ("visco");
            fieldcoeffs.push_back(evm_visc_Fwd);

            if(m_use_dynamic_alpha)
             {
              Array<OneD, NekDouble> evm_sensor_Fwd(nCoeffs);
              m_fields[0]->FwdTrans_IterPerExp(m_evm_sensor,  evm_sensor_Fwd);
              variables.push_back  ("alpha");
              fieldcoeffs.push_back(evm_sensor_Fwd);
             }
          }

        }
    }
} //end of namespace
