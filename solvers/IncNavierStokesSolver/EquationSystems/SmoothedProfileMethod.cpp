///////////////////////////////////////////////////////////////////////////////
//
// File SmoothedProfileMethod.cpp
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
// Description: Smoothed Profile Method for the Incompressible
// Navier Stokes equations
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/SmoothedProfileMethod.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>

using namespace std;

namespace Nektar
{
    using namespace MultiRegions;

    string SmoothedProfileMethod::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "SmoothedProfileMethod",
            SmoothedProfileMethod::create);

    /**
     * @brief Construct a new Smoothed Profile Method object
     *
     * @param pSession
     * @param pGraph
     */
    SmoothedProfileMethod::SmoothedProfileMethod(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
        : UnsteadySystem(pSession, pGraph),
          VelocityCorrectionScheme(pSession, pGraph)
    {

    }

    /**
     * @brief Destroy the Smoothed Profile Method object
     *
     */
    SmoothedProfileMethod::~SmoothedProfileMethod(void)
    {

    }

    void SmoothedProfileMethod::v_InitObject()
    {
        VelocityCorrectionScheme::v_InitObject();

        // Update implicit time-intregration class operators
        m_ode.DefineImplicitSolve(
            &SmoothedProfileMethod::SolveUnsteadyStokesSystem, this);

        // Number of dims as number of velocity vectors
        int nvel = m_velocity.num_elements();

        // Initialization of correction pressure and shape function
        switch (nvel)
        {
            case 1:
                if (m_projectionType == eGalerkin)
                {
                    m_pressureP = MemoryManager<ContField1D>::
                        AllocateSharedPtr(
                            *dynamic_pointer_cast<ContField1D>(m_pressure));
                    m_phi = MemoryManager<ContField1D>::
                        AllocateSharedPtr(
                            *dynamic_pointer_cast<ContField1D>(m_pressure));
                }
                else if (m_projectionType == eDiscontinuous)
                {
                    m_pressureP = MemoryManager<DisContField1D>::
                        AllocateSharedPtr(
                            *dynamic_pointer_cast<DisContField1D>(m_pressure));
                    m_phi = MemoryManager<DisContField1D>::
                        AllocateSharedPtr(
                            *dynamic_pointer_cast<DisContField1D>(m_pressure));
                }
                break;

            case 2:
                if (m_projectionType == eGalerkin)
                {
                    m_pressureP = MemoryManager<ContField2D>::
                        AllocateSharedPtr(
                            *dynamic_pointer_cast<ContField2D>(m_pressure));
                    m_phi = MemoryManager<ContField2D>::
                        AllocateSharedPtr(
                            *dynamic_pointer_cast<ContField2D>(m_pressure));
                }
                else if (m_projectionType == eDiscontinuous)
                {
                    m_pressureP = MemoryManager<DisContField2D>::
                        AllocateSharedPtr(
                            *dynamic_pointer_cast<DisContField2D>(m_pressure));
                    m_phi = MemoryManager<DisContField2D>::
                        AllocateSharedPtr(
                            *dynamic_pointer_cast<DisContField2D>(m_pressure));
                }
                break;

            case 3:
                if (m_projectionType == eGalerkin)
                {
                    if (m_HomogeneousType == EquationSystem::eNotHomogeneous)
                    {
                        m_pressureP = MemoryManager<ContField3D>::
                            AllocateSharedPtr(
                                *dynamic_pointer_cast<ContField3D>(
                                    m_pressure));
                        m_phi = MemoryManager<ContField3D>::
                            AllocateSharedPtr(
                                *dynamic_pointer_cast<ContField3D>(
                                    m_pressure));
                    }
                    else if (m_HomogeneousType == EquationSystem::eHomogeneous1D)
                    {
                        m_pressureP = MemoryManager<ContField3DHomogeneous1D>::
                            AllocateSharedPtr(
                                *dynamic_pointer_cast<ContField3DHomogeneous1D>(
                                    m_pressure));
                        m_phi = MemoryManager<ContField3DHomogeneous1D>::
                            AllocateSharedPtr(
                                *dynamic_pointer_cast<ContField3DHomogeneous1D>(
                                    m_pressure));
                    }
                    else if (m_HomogeneousType == EquationSystem::eHomogeneous2D
                          || m_HomogeneousType == EquationSystem::eHomogeneous3D)
                    {
                        m_pressureP = MemoryManager<ContField3DHomogeneous2D>::
                            AllocateSharedPtr(
                                *dynamic_pointer_cast<ContField3DHomogeneous2D>(
                                    m_pressure));
                        m_phi = MemoryManager<ContField3DHomogeneous2D>::
                            AllocateSharedPtr(
                                *dynamic_pointer_cast<ContField3DHomogeneous2D>(
                                    m_pressure));
                    }
                }
                else if (m_projectionType == eDiscontinuous)
                {
                    if (m_HomogeneousType == EquationSystem::eNotHomogeneous)
                    {
                        m_pressureP = MemoryManager<DisContField3D>::
                            AllocateSharedPtr(
                                *dynamic_pointer_cast<DisContField3D>(
                                    m_pressure));
                        m_phi = MemoryManager<DisContField3D>::
                            AllocateSharedPtr(
                                *dynamic_pointer_cast<DisContField3D>(
                                    m_pressure));
                    }
                    else if (m_HomogeneousType == EquationSystem::eHomogeneous1D)
                    {
                        m_pressureP = MemoryManager<
                            DisContField3DHomogeneous1D>::AllocateSharedPtr(
                                *dynamic_pointer_cast<
                                    DisContField3DHomogeneous1D>(m_pressure));
                        m_phi = MemoryManager<DisContField3DHomogeneous1D>::
                            AllocateSharedPtr(
                                *dynamic_pointer_cast<
                                    DisContField3DHomogeneous1D>(m_pressure));
                    }
                    else if (m_HomogeneousType == EquationSystem::eHomogeneous2D
                          || m_HomogeneousType == EquationSystem::eHomogeneous3D)
                    {
                        m_pressureP = MemoryManager<
                            DisContField3DHomogeneous2D>::AllocateSharedPtr(
                                *dynamic_pointer_cast<
                                    DisContField3DHomogeneous2D>(m_pressure));
                        m_phi = MemoryManager<DisContField3DHomogeneous2D>::
                            AllocateSharedPtr(
                                *dynamic_pointer_cast<
                                    DisContField3DHomogeneous2D>(m_pressure));
                    }
                }
                break;
        }

        // Get function evaluator for 'm_phi' from the session file
        if (m_session->DefinesFunction("ShapeFunction"))
        {
            // Get the evaluator
            m_phiEvaluator = GetFunction("ShapeFunction");
            m_timeDependentPhi = GetFunctionTimeDependence("ShapeFunction",
                                                           "E");
        }
        else
        {
            ASSERTL0(true, "ShapeFunction must be defined in "
                           "the session file.");
        }

        // Read 'u_p' from session file
        int physTot = m_pressureP->GetTotPoints();
        m_velName.push_back("Up");
        if (nvel > 1)
        {
            m_velName.push_back("Vp");
        }
        if (nvel > 2)
        {
            m_velName.push_back("Wp");
        }

        m_up = Array<OneD, Array<OneD, NekDouble> >(nvel);
        for (int i = 0; i < nvel; ++i)
        {
            m_up[i] = Array<OneD, NekDouble>(physTot);
        }
        if (m_session->DefinesFunction("ParticleVelocity"))
        {
            // Get the evaluator
            m_upEvaluator = GetFunction("ParticleVelocity");
            m_timeDependentUp = GetFunctionTimeDependence("ParticleVelocity",
                                                          "E");
        }
        else
        {
            ASSERTL0(true, "ParticleVelocity must be defined in "
                           "the session file.");
        }

        // Make sure that m_phi and m_up are defined
        UpdatePhiUp(0.0);

        // Select m_gamma0 depending on IMEX order
        string type = m_session->GetSolverInfo("TimeIntegrationMethod");
        switch (type.back()-'0')
        {
        case 1:
            m_gamma0 = 1.0;
            break;

        case 2:
            m_gamma0 = 3.0/2.0;
            break;

        case 3:
            m_gamma0 = 11.0/6.0;
            break;
        }
    }

    /**
     * @brief Generates the summary of the current simulation
     *
     * @param s
     */
    void SmoothedProfileMethod::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        VelocityCorrectionScheme::v_GenerateSummary(s);
        SolverUtils::AddSummaryItem(s, "IB formulation",
                                       "Smoothed Profile Method (SPM)");
    }

    /**
     * @brief Linear terms due to pressure and visosity are calculated here.
     * After solving the velocity filed without taking into account the
     * immersed boundaries, a new correction is applied through the force
     * \f[f_s\f]:
     *
     * \f[ \mathbf{f_s} = \frac{\Phi^{n+1}(\mathbf{u_p}-\mathbf{u^*})}
     * {\Delta t} \f]
     *
     * @param inarray
     * @param outarray
     * @param time
     * @param a_iixDt
     */
    void SmoothedProfileMethod::v_SolveUnsteadyStokesSystem(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble time,
                    const NekDouble a_iixDt)
    {
        VelocityCorrectionScheme::SolveUnsteadyStokesSystem(inarray,
                                                            outarray,
                                                            time,
                                                            a_iixDt);

        int physTot = m_pressureP->GetTotPoints();

        /* SPM correction of velocity */
        // Update 'm_phi' and 'm_up' if needed (evaluated at next time step)
        UpdatePhiUp(time + a_iixDt);
        // DEBUG: Test EstimateForces function
        EstimateForces(outarray, a_iixDt);
        // Set BC conditions for pressure p_p
        SetUpCorrectionPressure(outarray, m_F, time, a_iixDt);
        // Solve Poisson equation for pressure p_p
        SolveCorrectionPressure(m_F[0]);
        // Solve velocity in the next step with IB
        SolveCorrectedVelocity(m_F, outarray, time, a_iixDt);

        // Add pressures to get final value
        Vmath::Vadd(physTot, m_pressure->GetPhys(), 1,
                             m_pressureP->GetPhys(), 1,
                             m_pressure->UpdatePhys(), 1);
    }

    /**
     * @brief Sets the forcing term of the equation for the correction pressure
     * \f[p_p\f]:
     *
     * \f[ \nabla\cdot\mathbf{f_s} \f]
     *
     * @param fields
     * @param Forcing
     * @param aii_Dt
     */
    void SmoothedProfileMethod::SetUpCorrectionPressure(
                    const Array<OneD, const Array<OneD, NekDouble> > &fields,
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    NekDouble time,
                    NekDouble aii_Dt)
    {
        int physTot = m_pressureP->GetTotPoints();
        int nvel    = m_velocity.num_elements();

        // DEBUG: Set boundary conditions
        SetCorrectionPressureBCs(time, aii_Dt);

        // Virtual force 'fs'
        Array<OneD, Array<OneD, NekDouble> > f_s;
        IBForcing(fields, time, aii_Dt, f_s);
        m_fields[m_velocity[0]]->PhysDeriv(eX, f_s[0], Forcing[0]);

        // Using 'Forcing[1]' as storage
        for (int i = 1; i < nvel; ++i)
        {
            int ind = m_velocity[i];
            m_fields[ind]->PhysDeriv(DirCartesianMap[i], f_s[i], Forcing[1]);
            Vmath::Vadd(physTot, Forcing[1], 1, Forcing[0], 1, Forcing[0], 1);
        }
    }

    /**
     * @brief Solves the Poisson equation for the correction pressure
     * \f[p_p\f]:
     *
     * \f[ \nabla^2 p_p = \nabla\cdot\mathbf{f_s} \f]
     *
     * @param Forcing
     */
    void SmoothedProfileMethod::SolveCorrectionPressure(
                    const Array<OneD, NekDouble> &Forcing)
    {
        StdRegions::ConstFactorMap factors;
        // Factor 'lambda=0' in Helmholtz equation to get the Poisson form
        factors[StdRegions::eFactorLambda] = 0.0;

        // Solve the Poisson equation
        m_pressureP->HelmSolve(Forcing, m_pressureP->UpdateCoeffs(),
                               NullFlagList, factors);
        // Update node values from coefficients
        m_pressureP->BwdTrans(m_pressureP->GetCoeffs(),
                              m_pressureP->UpdatePhys());

        // DEBUG: AddPressureToOutflowBCs?
    }

    /**
     * @brief
     *
     * @param Forcing
     * @param fields
     * @param dt
     */
    void SmoothedProfileMethod::SolveCorrectedVelocity(
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    Array<OneD, Array<OneD, NekDouble> > &fields,
                    NekDouble time,
                    NekDouble dt)
    {
        int physTot = m_fields[0]->GetTotPoints();

        // Gradient of p_p
        int nvel = m_velocity.num_elements();
        if (nvel == 2)
        {
            m_pressureP->PhysDeriv(m_pressureP->GetPhys(),
                                   Forcing[0],
                                   Forcing[1]);
        }
        else
        {
            m_pressureP->PhysDeriv(m_pressureP->GetPhys(),
                                   Forcing[0],
                                   Forcing[1],
                                   Forcing[2]);
        }

        // Virtual force 'fs'
        Array<OneD, Array<OneD, NekDouble> > f_s;
        IBForcing(fields, time, dt, f_s);

        // Velocity correction
        for (int i = 0; i < nvel; ++i)
        {
            int ind = m_velocity[i];

            // DEBUG: Try adding -(1-m_phi)*grad(p_p) instead of -grad(p_p)
            // Vmath::Vvtvm(physTot, m_phi->GetPhys(), 1, Forcing[i], 1,
            //                                            Forcing[i], 1,
            //                                            Forcing[i], 1);
            // Vmath::Vadd(physTot, f_s[i], 1, Forcing[i], 1, Forcing[i], 1);
            Vmath::Vsub(physTot, f_s[i], 1, Forcing[i], 1, Forcing[i], 1);
            Blas::Daxpy(physTot, dt/m_gamma0, Forcing[i], 1, fields[ind], 1);
        }
    }

    /**
     * @brief DEBUG: Updates the BCs for boundaries with Dirichlet BCs in the
     * velocity:
     *
     * \f[ \frac{\partial p_p}{\partial\mathbf{n}} =
     *     \mathbf{f_s}\cdot\mathbf{n} \f]
     *
     * @param dt
     */
    void SmoothedProfileMethod::SetCorrectionPressureBCs(NekDouble time,
                                                           NekDouble dt)
    {
        Array<OneD, ExpListSharedPtr> BndExp;
        Array<OneD, SpatialDomains::BoundaryConditionShPtr> BndCond;

        // Get the BC expansions
        BndExp  = m_pressureP->GetBndCondExpansions();
        BndCond = m_pressureP->GetBndConditions();

        // For each boundary...
        for (int b = 0; b < BndExp.num_elements(); ++b)
        {
            // Skip this step for non Neumann BCs
            if (BndCond[b]->GetBoundaryConditionType() !=
                SpatialDomains::eNeumann)
            {
                continue;
            }

            // Calculate f_s values
            Array<OneD, Array<OneD, NekDouble> > f_s;
            IBForcingBC(b, BndExp[b], time, dt, f_s);

            // BC is f_s * n
            BndExp[b]->NormVectorIProductWRTBase(f_s, BndExp[b]->UpdatePhys());
        }
    }

    /**
     * @brief Calculates the values of the shape function
     *
     * @param expansion
     * @param t
     * @param phi
     */
    void SmoothedProfileMethod::UpdatePhiUp(NekDouble t)
    {
        // Calculate at least once (more times if m_phi or m_up depend on time)
        if (m_timeDependentPhi || t <= 0.0)
        {
            m_phiEvaluator->Evaluate("Phi", m_phi->UpdatePhys(), t);
        }
        if (m_timeDependentUp || t <= 0.0)
        {
            if (t <= 0.0)
            {
                // Initialize both variables for the first step
                m_upEvaluator->Evaluate(m_velName, m_up, t);
                m_upPrev = m_up;
            }
            else
            {
                // Store previous value of u_p during simulation
                m_upPrev = m_up;
                m_upEvaluator->Evaluate(m_velName, m_up, t);
            }
        }
    }

    /**
     * @brief For a body with a constant velocity \f[\mathbf{u_p}\f], the force
     * \f[\mathbf{f_s}\f] applied to the fluid ensures that the IBC are met:
     *
     * \f[ \mathbf{f_s} = \frac{\Phi^{n+1}\left(\mathbf{u_p} -
     * \mathbf{u^*}\right)}{\Delta t} \f]
     *
     * @param fields
     * @param dt
     * @param f_s
     */
    void SmoothedProfileMethod::IBForcing(
                    const Array<OneD, const Array<OneD, NekDouble> > &fields,
                    NekDouble time,
                    NekDouble dt,
                    Array<OneD, Array<OneD, NekDouble> > &f_s)
    {
        int nvel = m_velocity.num_elements();
        int nq   = m_pressureP->GetTotPoints();

        // Vector f_s
        f_s = Array<OneD, Array<OneD, NekDouble> >(nvel);
        for (int i = 0; i < nvel; ++i)
        {
            f_s[i] = Array<OneD, NekDouble>(nq);
        }

        for (int i = 0; i < nvel; ++i)
        {
            Vmath::Vsub(nq, m_up[i], 1, fields[m_velocity[i]], 1, f_s[i], 1);
            Vmath::Vmul(nq, m_phi->GetPhys(), 1, f_s[i], 1, f_s[i], 1);
            Vmath::Smul(nq, m_gamma0/dt, f_s[i], 1, f_s[i], 1);
        }
    }

    /**
     * @brief For a body with a constant velocity \f[\mathbf{u_p}\f], the force
     * \f[\mathbf{f_s}\f] applied to the fluid ensures that the IBC are met.
     * Calculated in the boundary defined by 'BndExp'
     *
     * @param bndInd
     * @param BndExp
     * @param dt
     * @param f_s
     */
    void SmoothedProfileMethod::IBForcingBC(int bndInd,
                                const ExpListSharedPtr &BndExp,
                                NekDouble time,
                                NekDouble dt,
                                Array<OneD, Array<OneD, NekDouble> > &f_s)
    {
        int nvel = m_velocity.num_elements();
        int nq   = BndExp->GetTotPoints();

        // Vector f_s
        f_s = Array<OneD, Array<OneD, NekDouble> >(nvel);
        for (int i = 0; i < nvel; ++i)
        {
            f_s[i] = Array<OneD, NekDouble>(nq);
        }

        for (int i = 0; i < nvel; ++i)
        {
            ExpListSharedPtr velExp =
                (m_fields[m_velocity[i]]->GetBndCondExpansions())[bndInd];
            Vmath::Vsub(nq, m_up[i], 1, velExp->GetPhys(), 1, f_s[i], 1);
            Vmath::Vmul(nq, m_phi->GetBndCondExpansions()[bndInd]->GetPhys(),
                        1, f_s[i],
                        1, f_s[i], 1);
            Vmath::Smul(nq, m_gamma0/dt, f_s[i], 1, f_s[i], 1);
        }
    }

    /**
     * @brief True if the function is timedependent, false otherwise
     *
     * @param name
     * @param type
     * @param attribute
     * @return string
     */
    bool SmoothedProfileMethod::GetFunctionTimeDependence(string name,
                                                          string type)
    {
        // Get the handler of first function block
        TiXmlElement *conds = m_session->GetElement("Nektar/Conditions");
        TiXmlElement *function = conds->FirstChildElement("FUNCTION");

        // Loop over functions until 'name' block
        string functionType = function->Attribute("NAME");
        while (function && !boost::iequals(functionType, name))
        {
            function = function->NextSiblingElement("FUNCTION");
            functionType = function->Attribute("NAME");
        }

        // Go to the first element
        TiXmlElement *functionDef = function->FirstChildElement(type.c_str());
        ASSERTL0(functionDef,"At least one element must be defined in" + name);

        // And return the value of TIMEDEPENDENT
        bool output;
        int err = functionDef->QueryBoolAttribute("TIMEDEPENDENT", &output);
        ASSERTL0(err != TIXML_WRONG_TYPE, "TIMEDEPENDENT must be 1/0, "
                                          "true/false or yes/no.")
        if (err == TIXML_NO_ATTRIBUTE)
        {
            // Set output to 'false' if no time-dependence specified
            output = false;
        }

        return output;
    }

    /**
     * @brief Determine the total force on the body defined by \f[\Phi\f]
     * (note that if the shape function represents more than one
     * body, this function calculates the value of the final force after adding
     * up the values for each body). This value must be scaled with the
     * density to get the real force vector.
     * 
     * SHOULD BE IMPLEMENTED AS A FILTER!!
     *
     * @param inarray
     * @param dt
     */
    void SmoothedProfileMethod::EstimateForces(
                    const Array<OneD, const Array<OneD, NekDouble> > &velocity,
                    NekDouble dt)
    {
        // DEBUG: Example of integration in FilterAeroForces.cpp
        // DEBUG: Only when forcing terms are added via Filters
        int nq   = m_pressureP->GetTotPoints();
        int nvel = m_velocity.num_elements();
        Array<OneD, NekDouble> tmp(nq);
        Array<OneD, NekDouble> forceTmp(nvel, 0.0);

        for (int i = 0; i < nvel; ++i)
        {
            // "Scalar" field to be integrated
            Vmath::Vsub(nq, velocity[m_velocity[i]], 1, m_upPrev[i], 1,
                        tmp, 1);
            Vmath::Vmul(nq, m_phi->GetPhys(), 1, tmp, 1, tmp, 1);

            // Integration of force throughout the domain
            forceTmp[i] = m_pressureP->Integral(tmp);
            forceTmp[i] /= dt;
        }
        m_Forces.push_back(forceTmp);
    }

} // end of namespace
