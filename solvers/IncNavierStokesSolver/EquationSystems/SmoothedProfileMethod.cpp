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
        ASSERTL0(m_session->DefinesFunction("ShapeFunction"),
                 "ShapeFunction must be defined in the session file.")
        ReadPhi();

        // Allocate the vector 'm_up'
        int physTot = m_pressureP->GetTotPoints();
        m_velName.push_back("Up");
        switch (nvel)
        {
            case (3):
                m_velName.push_back("Wp");
            case (2):
                m_velName.push_back("Vp");
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
    void SmoothedProfileMethod::v_GenerateSummary(SolverUtils::SummaryList &s)
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
        SetUpCorrectionPressure(outarray, m_F, a_iixDt);
        // Solve Poisson equation for pressure p_p
        SolveCorrectionPressure(m_F[0]);
        // Solve velocity in the next step with IB
        SolveCorrectedVelocity(m_F, outarray, a_iixDt);

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
                    NekDouble aii_Dt)
    {
        int physTot = m_pressureP->GetTotPoints();
        int nvel    = m_velocity.num_elements();

        // DEBUG: Set boundary conditions
        SetCorrectionPressureBCs(aii_Dt);

        // Virtual force 'fs'
        Array<OneD, Array<OneD, NekDouble> > f_s;
        IBForcing(fields, aii_Dt, f_s);
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

        // Update node values from coefficients only if not in homogeneous case
        if (m_HomogeneousType == EquationSystem::eNotHomogeneous)
        {
            m_pressureP->BwdTrans(m_pressureP->GetCoeffs(),
                                  m_pressureP->UpdatePhys());
        }
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
        IBForcing(fields, dt, f_s);

        // Velocity correction
        for (int i = 0; i < nvel; ++i)
        {
            int ind = m_velocity[i];

            // Adding -(1-m_phi)*grad(p_p) instead of -grad(p_p) reduces the
            // flux through the walls, but the flow is not incompressible
            if (m_session->DefinesSolverInfo("ForceBoundary") &&
                boost::iequals(m_session->GetSolverInfo("ForceBoundary"), "1"))
            {
                Vmath::Vvtvm(physTot, m_phi->GetPhys(), 1, Forcing[i], 1,
                                                        Forcing[i], 1,
                                                        Forcing[i], 1);
                Vmath::Vadd(physTot, f_s[i], 1, Forcing[i], 1, Forcing[i], 1);
            }
            else
            {
                Vmath::Vsub(physTot, f_s[i], 1, Forcing[i], 1, Forcing[i], 1);
                Blas::Daxpy(physTot, dt/m_gamma0, Forcing[i], 1, fields[ind], 1);
            }
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
    void SmoothedProfileMethod::SetCorrectionPressureBCs(NekDouble dt)
    {
        Array<OneD, ExpListSharedPtr> BndExp;
        Array<OneD, SpatialDomains::BoundaryConditionShPtr> BndCond;

        // Get the BC expansions
        BndExp  = m_pressureP->GetBndCondExpansions();
        BndCond = m_pressureP->GetBndConditions();

        // For each boundary...
        for (int b = 0; b < BndExp.num_elements(); ++b)
        {
            // Only for BCs based on the derivative
            if (BndCond[b]->GetBoundaryConditionType() ==
                SpatialDomains::eNeumann ||
                BndCond[b]->GetBoundaryConditionType() ==
                SpatialDomains::ePeriodic)
            {
                // Calculate f_s values
                Array<OneD, Array<OneD, NekDouble> > f_s;
                IBForcingBC(b, BndExp[b], dt, f_s);

                // BC is f_s * n
                BndExp[b]->NormVectorIProductWRTBase(f_s, BndExp[b]->UpdatePhys());
            }
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

        // And return the value of USERDEFINEDTYPE
        string attr;
        int err = functionDef->QueryStringAttribute("USERDEFINEDTYPE", &attr);
        bool output = boost::iequals(attr, "TimeDependent");

        ASSERTL0((err == TIXML_NO_ATTRIBUTE) ||
                 (err == TIXML_SUCCESS && output), "USERDEFINEDTYPE in " +
                 elemName + " must be TimeDependent if defined");

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
            
            // Get phi values from XML file (after "FieldConvert" the STL file)
            // First, load the data
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> fieldDef;
            std::vector<std::vector<NekDouble> > fieldData;
            LibUtilities::FieldMetaDataMap fieldMetaData;
            LibUtilities::FieldIOSharedPtr phiFile =
                LibUtilities::FieldIO::CreateForFile(m_session, fileName);
            phiFile->Import(fileName, fieldDef, fieldData, fieldMetaData);

            // Extract Phi field to output
            string tmp("Phi");
            for (int i = 0; i < fieldData.size(); ++i)
            {
                m_phi->ExtractDataToCoeffs(
                    fieldDef[i],
                    fieldData[i],
                    tmp,
                    m_phi->UpdateCoeffs());
            }
            m_filePhi = true;
        }

        // Check if Phi is timedependent
        m_timeDependentPhi = GetVarTimeDependence("ShapeFunction", "Phi");

        // If so, check if its velocity changes as well
        m_timeDependentUp = GetVarTimeDependence("ShapeFunction", "Up");
        switch (m_velocity.num_elements())
        {
            case 3:
                m_timeDependentUp |= GetVarTimeDependence("ShapeFunction", "Wp");
            case 2:
                m_timeDependentUp |= GetVarTimeDependence("ShapeFunction", "Vp");
        }
    }

} // end of namespace
