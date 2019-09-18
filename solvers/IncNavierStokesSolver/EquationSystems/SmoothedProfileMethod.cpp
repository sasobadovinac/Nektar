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
#include <IncNavierStokesSolver/Filters/FilterAeroForcesSPM.h>
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

        // Read 'm_phi' and its velocity
        if (m_session->DefinesFunction("ShapeFunction"))
        {
            ReadPhi();
        }
        else
        {
            ASSERTL0(false, "ShapeFunction must be defined in "
                            "the session file.")
        }
        // Allocate the vector 'm_up'
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

        // Check if the aeroforces filter is active, negative if inactive
        m_forcesFilter = -1;
        for (int i = 0; i < m_session->GetFilters().size(); ++i)
        {
            if (m_session->GetFilters()[i].first == "AeroForcesSPM")
            {
                m_forcesFilter = i;
            }
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
        // Estimate forces only if requested
        if (m_forcesFilter >= 0)
        {
            dynamic_pointer_cast<FilterAeroForcesSPM>(
                m_filters[m_forcesFilter])->CalculateForces(outarray, m_upPrev,
                                            m_phi, time, a_iixDt);
        }
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
        m_pressure->FwdTrans(m_pressure->GetPhys(),
                             m_pressure->UpdateCoeffs());
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
     * @param time
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
     * @param time
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
     * @param t
     */
    void SmoothedProfileMethod::UpdatePhiUp(NekDouble t)
    {
        // Initialise 'm_up' and 'm_phi' during first step
        if (t <= 0.0)
        {
            // Update 'm_phi' only if it was provided as an analytical function
            if (!m_filePhi)
            {
                m_phiEvaluator->Evaluate("Phi", m_phi->UpdatePhys(), t);
            }
            // Initialize both variables for the first step
            m_phiEvaluator->Evaluate(m_velName, m_up, t);
            m_upPrev = m_up;
        }
        // If timedependent 'm_phi'...
        else if (m_timeDependentPhi)
        {
            // ...and not loaded from an external file
            if (!m_filePhi)
            {
                m_phiEvaluator->Evaluate("Phi", m_phi->UpdatePhys(), t);
                // And if velocities are timedependent as well
                if (m_timeDependentUp)
                {
                    // Store previous value of u_p during simulation
                    m_upPrev = m_up;
                    m_phiEvaluator->Evaluate(m_velName, m_up, t);
                }
            }
            // DEBUG ...and loaded from a file
            else
            {
                // With constant velocity, the new position of the IBs is just
                // x = v*dt
                if (!m_timeDependentUp)
                {
                    
                }
                // Otherwise, the values must be integrated in time by means of
                // an algorithm such as RK...
                else
                {
                    
                }
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
    bool SmoothedProfileMethod::GetVarTimeDependence(string funcName,
                                                     string elemName)
    {
        // Get the handler of the function
        TiXmlElement *function = GetFunctionHdl(funcName);

        // Go to the first element
        TiXmlElement *functionDef = function->FirstChildElement();
        ASSERTL0(functionDef, "At least one element must be defined in " +
                              funcName)

        // And search the element with name 'elemName'
        string varName = functionDef->Attribute("VAR");
        while(functionDef && !boost::iequals(varName, elemName))
        {
            functionDef = functionDef->NextSiblingElement();
            varName = functionDef->Attribute("VAR");
        }

        ASSERTL0(functionDef, "Variable " + elemName + " must be defined in " +
                              funcName + ".");

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
     * @brief Returns a handle to the requested function. Returns NULL if it
     * does not exist
     * 
     * @param functionName 
     * @return TiXmlElement* 
     */
    TiXmlElement* SmoothedProfileMethod::GetFunctionHdl(string functionName)
    {
        // Get the handler of first function block
        TiXmlElement *conds = m_session->GetElement("Nektar/Conditions");
        TiXmlElement *function = conds->FirstChildElement("FUNCTION");

        // Loop over functions until the block 'name' is found
        string functionType = function->Attribute("NAME");
        while (function && !boost::iequals(functionType, functionName))
        {
            function = function->NextSiblingElement("FUNCTION");
            functionType = function->Attribute("NAME");
        }
        
        return function;
    }

    void SmoothedProfileMethod::ReadPhi()
    {
        // Function evaluator for Phi and Up
        m_phiEvaluator = GetFunction("ShapeFunction");

        TiXmlElement *function = GetFunctionHdl("ShapeFunction");
        TiXmlElement *child    = function->FirstChildElement();
        m_filePhi = false;

        // If defined by using a file
        if (boost::iequals(child->ValueStr(), "F"))
        {
            // Get name of STL file
            string fileName;
            int status = child->QueryStringAttribute("FILE", &fileName);
            ASSERTL0(status == TIXML_SUCCESS, "An STL file with the geometry "
                     "of the immersed bodies has to be supplied.")
            ASSERTL0(boost::iequals(fileName.substr(fileName.length()-4),
                     ".stl"), "A valid STL file must be supplied in the "
                     "'ShapeFunction' field.")
            
            // Get phi values from file
            STLfile file = ReadSTL(fileName);
            GetPhifromSTL(file);
            m_filePhi = true;
        }

        // Check if Phi is timedependent
        m_timeDependentPhi = GetVarTimeDependence("ShapeFunction", "Phi");

        // If so, check if its velocity changes as well
        m_timeDependentUp = GetVarTimeDependence("ShapeFunction", "Up");
    }

    /**
     * @brief Read one 3D vector from a STL file, starting from the next line
     * of the input 'ifstream'. Numbers in ifstream are defined as 'float'
     * 
     * @param in 
     * @return Array<OneD, NekDouble> 
     */
    Array<OneD, NekDouble> SmoothedProfileMethod::ReadVector(ifstream &in)
    {
        Array<OneD, NekDouble> out(3);
        char buf[4];

        in.read(buf, 4);
        out[0] = *((float*) buf);
        in.read(buf, 4);
        out[1] = *((float*) buf);
        in.read(buf, 4);
        out[2] = *((float*) buf);

        return out;
    }

    /**
     * @brief Read an STL binary file and returns a struct of type 'STLfile'
     * containing the parsed data
     * 
     * @param filename 
     * @return SmoothedProfileMethod::STLfile 
     */
    SmoothedProfileMethod::STLfile SmoothedProfileMethod::ReadSTL(string filename)
    {
        STLfile out;

        // Open file
        ifstream fileStl(filename.c_str(), ios::binary);
        ASSERTL0(fileStl, "An error occurred while trying to open the STL file.")

        // Buffers
        char headerBuf[80];
        char numTriBuf[4];
        char dumpBuf[2];

        // Read header and num of triangles
        fileStl.read(headerBuf, 80);
        fileStl.read(numTriBuf, 4);
        unsigned int numTri = *((unsigned int*) numTriBuf);

        out.header = headerBuf;
        out.numTri = numTri;

        // Read triangle data
        out.triangles = Array<OneD, triangle>(numTri);
        for (uint i = 0; i < numTri; ++i)
        {
            // Read normal vector
            triangle tmpTri;
            tmpTri.normal = ReadVector(fileStl);

            // Read three vertices
            tmpTri.v0 = ReadVector(fileStl);
            tmpTri.v1 = ReadVector(fileStl);
            tmpTri.v2 = ReadVector(fileStl);
            out.triangles[i] = tmpTri;

            // Dump triangle type
            fileStl.read(dumpBuf, 2);
        }

        // Close the file
        fileStl.close();

        return out;
    }

    /**
     * @brief Smoothing function for the SPM method given a distance value
     * and a scaling coefficient
     * 
     * @param dist 
     * @param coeff 
     * @return double 
     */
    double SmoothedProfileMethod::PhiFunction(double dist, double coeff)
    {
        return -0.5*(std::tanh(dist/coeff)-1.0);
    }

    /**
     * @brief Assigns to 'm_phi' the corresponding values of Phi
     * 
     * @param file 
     * @param phi 
     */
    void SmoothedProfileMethod::GetPhifromSTL(
                                    const SmoothedProfileMethod::STLfile &file)
    {
        // Get size of domain points
        size_t n = m_phi->GetTotPoints();
        int nvel = m_velocity.num_elements();

        // Get coordinates of nodes
        Array<OneD, NekDouble> x0(n);
        Array<OneD, NekDouble> x1(n);
        Array<OneD, NekDouble> x2(n);
        m_phi->GetCoords(x0, x1, x2);
        
        // Parallelisation is highly recommended here
        for (size_t i = 0; i < n; ++i)
        {
            // Cet current coords of the point
            Array<OneD, NekDouble> coords(3);
            coords[0] = x0[i]; coords[1] = 0.0; coords[2] = 0.0;
            if (nvel > 1)
            {
                coords[1] = x1[i];
            }
            if (nvel > 2)
            {
                coords[2] = x2[i];
            }

            // Find the shortest distance to the body(ies)
            double dist;
            bool inside = IsInterior(file, coords);
            FindShortestDist(file, coords, dist);
            if (inside)
            {
                dist = -dist;
            }

            // Get corresponding value of Phi
            m_phi->UpdatePhys()[i] = PhiFunction(dist, 0.01);
        }
    }

    /**
     * @brief Checks if a ray traced from 'Origin' with direction 'Dvec' hits
     * the triangle defined by 'tri'. Returns the distance to the plane
     * defined by 'tri' in any case. A negative distance means that the hit
     * happend in the direction oposite that of the ray
     * 
     * @param tri 
     * @param Origin 
     * @param Dvec 
     * @param distance 
     * @param u 
     * @param v 
     * @return true 
     * @return false 
     */
    bool SmoothedProfileMethod::CheckHit(
                                    const SmoothedProfileMethod::triangle &tri,
                                    const Array<OneD, NekDouble> &Origin,
                                    const Array<OneD, NekDouble> &Dvec,
                                    double &distance, double &u, double &v)
    {
        // Edge vectors
        Array<OneD, NekDouble> E1(3);
        Array<OneD, NekDouble> E2(3);
        for (int i = 0; i < 3; ++i)
        {
            E1[i] = tri.v1[i]-tri.v0[i];
            E2[i] = tri.v2[i]-tri.v0[i];
        }

        // If det == 0, ray parallel to triangle
        Array<OneD, NekDouble> Pvec = Cross(Dvec, E2);
        double det = Vmath::Dot(3, Pvec, E1);
        double inv_det = 1.0 / det;
        if (IsZero(det))
        {
            distance = std::numeric_limits<double>::infinity();
            u        = std::numeric_limits<double>::infinity();
            v        = std::numeric_limits<double>::infinity();
            return false;
        }

        // Vector T and parameter u = (0.0, 1.0)
        Array<OneD, NekDouble> Tvec(3);
        for (int i = 0; i < 3; ++i)
        {
            Tvec[i] = Origin[i]-tri.v0[i];
        }
        u = Vmath::Dot(3, Pvec, Tvec) * inv_det;

        // Vector Q and parameter v = (0.0, 1.0)
        Array<OneD, NekDouble> Qvec = Cross(Tvec, E1);
        v = Vmath::Dot(3, Qvec, Dvec) * inv_det;

        // There is a hit if (u,v) coordinates are bounded
        distance = Vmath::Dot(3, Qvec, E2) * inv_det;
        if ((u < 0.0 || u > 1.0) || (v < 0.0 || u+v > 1.0))
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    /**
     * @brief Returns true if a point is inside the 3D object defined in the
     * STL file. It is based in the idea that a ray traced from inside a closed
     * surface will go through an odd number of surfaces, no matter how complex
     * the geometry is
     * 
     * @param file 
     * @param x 
     * @return true 
     * @return false 
     */
    bool SmoothedProfileMethod::IsInterior(const STLfile &file,
                                           const Array<OneD, NekDouble> &x)
    {
        // Direction does not matter
        Array<OneD, NekDouble> dir(3);
        dir[0]   = 1.0;
        dir[1]   = 0.0;
        dir[2]   = 0.0;
        int hits = 0;

        // Stores the distances of the hits with each surface
        // It has to be a dynamic container, it will always be small
        vector<NekDouble> distVec;

        // Check hits with all the triangles
        for (triangle tri : file.triangles)
        {
            double dist;
            double u, v;
            bool hit = CheckHit(tri, x, dir, dist, u, v);

            if (hit && dist > 0.0 &&
                std::find_if(distVec.begin(), distVec.end(),
                    [&](double x){ return IsZero(x-dist); }) == distVec.end())
            {
                distVec.push_back(dist);
                hits++;
            }
        }

        // Odd number of hits -> the point lies INSIDE
        if (hits % 2)
        {
            return true;
        }
        // Otherwise, it falls OUTSIDE
        else
        {
            return false;
        }
    }

    void SmoothedProfileMethod::FindShortestDist(
                                    const SmoothedProfileMethod::STLfile &file,
                                    const Array<OneD, NekDouble> &x,
                                    double &dist)
    {
        // Set 'dist' to an unreal value
        dist = numeric_limits<double>::infinity();

        for (triangle tri : file.triangles)
        {
            double tmpDist;
            double u, v;
            bool hit = CheckHit(tri, x, tri.normal, tmpDist, u, v);
            tmpDist  = abs(tmpDist);

            if (!hit)
            {
                // The minimum has to be in one of the edges
                if (v < 0)   // Edge V0-V1
                {
                    tmpDist = Distance2edge(x, tri.v0, tri.v1);
                }
                else if (u < 0)   // Edge V0-V2
                {
                    tmpDist = Distance2edge(x, tri.v0, tri.v2);
                }
                else   // Edge V1-V2
                {
                    tmpDist = Distance2edge(x, tri.v1, tri.v2);
                }
            }

            // Update 'dist'
            if (tmpDist < dist)
            {
                dist = tmpDist;
            }
        }
    }

    /**
     * @brief Returns true if the argument is CLOSE to zero. Tuned for
     * the STL parsing module, do not use for other purposes
     * 
     * @param x 
     * @return true 
     * @return false 
     */
    bool SmoothedProfileMethod::IsZero(double x)
    {
        double EPS = 0.000001;
        return (x > -EPS && x < EPS);
    }

    Array<OneD, NekDouble> SmoothedProfileMethod::Cross(
                                    const Array<OneD, NekDouble> &v0,
                                    const Array<OneD, NekDouble> &v1)
    {
        Array<OneD, NekDouble> out(3);

        out[0] = v0[1]*v1[2] - v0[2]*v1[1];
        out[1] = v0[2]*v1[0] - v0[0]*v1[2];
        out[2] = v0[0]*v1[1] - v0[1]*v1[0];

        return out;
    }

    /**
     * @brief Calculates the distance between two n-dimensional points
     * 
     * @param v0 
     * @param v1 
     * @return double 
     */
    double SmoothedProfileMethod::Distance2point(
                                    const Array<OneD, NekDouble> &v0,
                                    const Array<OneD, NekDouble> &v1)
    {
        size_t n   = v0.num_elements();
        double out = 0.0;

        for (size_t i = 0; i < n; ++i)
        {
            out += (v1[i]-v0[i]) * (v1[i]-v0[i]);
        }
        
        return sqrt(out);
    }

    /**
     * @brief Determines the shortest distance from a point 'x' to the segment
     * defined by the points 'e1' and 'e2'. Note that this distance may be
     * equal to that to one of the end points
     * 
     * @param x 
     * @param e1 
     * @param e2 
     * @return double 
     */
    double SmoothedProfileMethod::Distance2edge(
                                    const Array<OneD, NekDouble> &x,
                                    const Array<OneD, NekDouble> &e1,
                                    const Array<OneD, NekDouble> &e2)
    {
        size_t n = x.num_elements();
        Array<OneD, NekDouble> e1x(n);
        Array<OneD, NekDouble> e1e2(n);
        for (size_t i = 0; i < n; ++i)
        {
            e1x[i]  = x[i]-e1[i];
            e1e2[i] = e2[i]-e1[i];
        }
        double norm = sqrt(Vmath::Dot(n, e1e2, e1e2));
        for (size_t i = 0; i < n; ++i)
        {
            e1e2[i] /= norm;
        }

        double proj = Vmath::Dot(n, e1x, e1e2);
        if (proj < 0.0)
        {
            return Distance2point(x, e1);
        }
        else if (proj > norm)
        {
            return Distance2point(x, e2);
        }

        Array<OneD, NekDouble> distVec(n);
        for (size_t i = 0; i < n; ++i)
        {
            distVec[i] = e1x[i]-proj*e1e2[i];
        }

        return sqrt(Vmath::Dot(n, distVec, distVec));
    }

} // end of namespace
