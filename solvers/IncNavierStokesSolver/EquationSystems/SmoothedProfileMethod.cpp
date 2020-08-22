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

#include <iostream>
#include <fstream>

using namespace std;
# define my_sizeof(type) ((char *)(&type+1)-(char*)(&type))
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
        int nvel = m_velocity.size();

        // Initialization of correction pressure and shape function
        switch (nvel)
        {
            case 1:
                if (m_projectionType == eGalerkin)
                {
                    SetUpExpansions<ContField1D>(nvel);
                }
                else if (m_projectionType == eDiscontinuous)
                {
                    SetUpExpansions<DisContField1D>(nvel);
                }
                break;

            case 2:
                if (m_projectionType == eGalerkin)
                {
                    SetUpExpansions<ContField2D>(nvel);
                }
                else if (m_projectionType == eDiscontinuous)
                {
                    SetUpExpansions<DisContField2D>(nvel);
                }
                break;

            case 3:
                if (m_projectionType == eGalerkin)
                {
                    if (m_HomogeneousType == EquationSystem::eNotHomogeneous)
                    {
                        SetUpExpansions<ContField3D>(nvel);
                    }
                    else if (m_HomogeneousType ==
                             EquationSystem::eHomogeneous1D)
                    {
                        SetUpExpansions<ContField3DHomogeneous1D>(nvel);
                    }
                    else if (m_HomogeneousType ==
                             EquationSystem::eHomogeneous2D ||
                             m_HomogeneousType ==
                             EquationSystem::eHomogeneous3D)
                    {
                        SetUpExpansions<ContField3DHomogeneous2D>(nvel);
                    }
                }
                else if (m_projectionType == eDiscontinuous)
                {
                    if (m_HomogeneousType == EquationSystem::eNotHomogeneous)
                    {
                        SetUpExpansions<DisContField3D>(nvel);
                    }
                    else if (m_HomogeneousType ==
                             EquationSystem::eHomogeneous1D)
                    {
                        SetUpExpansions<DisContField3DHomogeneous1D>(nvel);
                    }
                    else if (m_HomogeneousType ==
                             EquationSystem::eHomogeneous2D ||
                             m_HomogeneousType ==
                             EquationSystem::eHomogeneous3D)
                    {
                        SetUpExpansions<DisContField3DHomogeneous2D>(nvel);
                    }
                }
                break;
        }

        // Read 'm_phi' and its velocity
        ASSERTL0(m_session->DefinesFunction("ShapeFunction"),
                 "ShapeFunction must be defined in the session file.")
        ReadPhi();

        // Allocate the vector 'm_up'
        int physTot = m_phi->GetTotPoints();
        m_velName.push_back("Up");
        if (nvel > 1)
        {
            m_velName.push_back("Vp");
        }
        if (nvel == 3)
        {
            m_velName.push_back("Wp");
        }

        m_up = Array<OneD, Array<OneD, NekDouble> >(nvel);
        for (int i = 0; i < nvel; ++i)
        {
            m_up[i] = Array<OneD, NekDouble>(physTot, 0.0);
        }

        // Make sure that m_phi and m_up are defined
        UpdatePhiUp(0.0);

        // Select 'm_gamma0' depending on IMEX order
        string intType = m_session->GetSolverInfo("TimeIntegrationMethod");
        ASSERTL0(boost::iequals(intType.substr(0, 9), "IMEXOrder"),
                 "TimeIntegrationMethod must be 'IMEXOrder1' to '4'.")
        switch (intType.back()-'0')
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

        case 4:
            m_gamma0 = 25.0/12.0;
            break;
        }

        // Check if the aeroforces filter is active, negative if inactive
        m_forcesFilter = -1;
        for (int i = 0; i < m_session->GetFilters().size(); ++i)
        {
            if (boost::iequals(m_session->GetFilters()[i].first,
                               "AeroForcesSPM"))
            {
                m_forcesFilter = i;
                break;
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
     * \f$f_s\f$:
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

        int physTot = m_pressureP->GetNpoints();

        /* SPM correction of velocity */
        // Update 'm_phi' and 'm_up' if needed (evaluated at next time step)
        UpdatePhiUp(time + a_iixDt);
        // Update calculation of IB forcing 'm_fs'
        UpdateForcing(outarray, a_iixDt);
        // Estimate forces only if requested
        if (m_forcesFilter >= 0)
        {
            static_pointer_cast<FilterAeroForcesSPM>(
                m_filters[m_forcesFilter].second)->CalculateForces(
                    outarray, m_upPrev, m_phi, time, a_iixDt);
        }
        // Set BC conditions for pressure p_p
        SetUpCorrectionPressure(outarray, m_F);
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

        // Add presure to outflow bc if using convective like BCs
        m_extrapolation->AddPressureToOutflowBCs(m_kinvis);

    }

    /**
     * @brief Sets the forcing term of the equation for the correction pressure
     * \f$p_p\f$:
     *
     * \f[ \nabla\cdot\mathbf{f_s} \f]
     *
     * @param fields
     * @param Forcing
     */
    void SmoothedProfileMethod::SetUpCorrectionPressure(
                    const Array<OneD, const Array<OneD, NekDouble> > &fields,
                    Array<OneD, Array<OneD, NekDouble> > &Forcing)
    {
        int physTot = m_fs[0]->GetNpoints();
        int nvel    = m_velocity.size();

        // Set boundary conditions
        SetCorrectionPressureBCs();

        // Divergence of 'fs'
        m_fields[m_velocity[0]]->PhysDeriv(eX, m_fs[0]->GetPhys(), Forcing[0]);

        // Using 'Forcing[1]' as storage
        for (int i = 1; i < nvel; ++i)
        {
            int ind = m_velocity[i];
            m_fields[ind]->PhysDeriv(DirCartesianMap[i],
                                     m_fs[i]->GetPhys(), Forcing[1]);
            Vmath::Vadd(physTot, Forcing[1], 1, Forcing[0], 1, Forcing[0], 1);
        }
    }

    /**
     * @brief Solves the Poisson equation for the correction pressure
     * \f$p_p\f$:
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
        m_pressureP->HelmSolve(Forcing, m_pressureP->UpdateCoeffs(), factors);

        // Update node values from coefficients
        m_pressureP->BwdTrans(m_pressureP->GetCoeffs(),
                              m_pressureP->UpdatePhys());
    }

    /**
     * @brief Corrects the velocity field so that the IBs are taken into
     * account. Solves the explicit equation:
     *
     * \f[ \frac{\gamma_0(\mathbf{u_p}^{n+1} - \mathbf{u}^*)}{\Delta t} =
     *     \mathbf{f_s} - \nabla p_p \f]
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
        int physTot = m_phi->GetNpoints();

        // Gradient of p_p
        int nvel = m_velocity.size();
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

        // Velocity correction
        for (int i = 0; i < nvel; ++i)
        {
            // Adding -(1-m_phi)*grad(p_p) instead of -grad(p_p) reduces the
            // flux through the walls, but the flow is not incompressible
            if (m_session->DefinesSolverInfo("ForceBoundary") &&
                boost::iequals(m_session->GetSolverInfo("ForceBoundary"),
                               "True"))
            {
                Vmath::Vvtvm(physTot, m_phi->GetPhys(), 1, Forcing[i], 1,
                             Forcing[i], 1, Forcing[i], 1);
                Vmath::Vadd(physTot, m_fs[i]->GetPhys(), 1, Forcing[i], 1,
                            Forcing[i], 1);
            }
            else
            {
                Vmath::Vsub(physTot, m_fs[i]->GetPhys(), 1, Forcing[i], 1,
                            Forcing[i], 1);
            }
            Blas::Daxpy(physTot, dt/m_gamma0, Forcing[i], 1, fields[i], 1);
        }
    }

    /**
     * @brief Updates the BCs for boundaries with Neumann or periodic BCs in
     * the pressure:
     *
     * \f[ \frac{\partial p_p}{\partial\mathbf{n}} =
     *     \mathbf{f_s}\cdot\mathbf{n} \f]
     */
    void SmoothedProfileMethod::SetCorrectionPressureBCs()
    {
        int nvel = m_velocity.size();
        Array<OneD, ExpListSharedPtr> BndExp;
        Array<OneD, SpatialDomains::BoundaryConditionShPtr> BndCond;

        // Get the BC expansions
        BndExp  = m_pressureP->GetBndCondExpansions();
        BndCond = m_pressureP->GetBndConditions();

        // For each boundary...
        for (int b = 0; b < BndExp.size(); ++b)
        {
            // Only for BCs based on the derivative
            if (BndCond[b]->GetBoundaryConditionType() ==
                SpatialDomains::eNeumann ||
                BndCond[b]->GetBoundaryConditionType() ==
                SpatialDomains::ePeriodic)
            {
                // Calculate f_s values
                Array<OneD, Array<OneD, NekDouble> > f_s(nvel);
                for (int i = 0; i < nvel; ++i)
                {
                    f_s[i] = m_fs[0]->GetBndCondExpansions()[b]->GetPhys();
                }

                // BC is f_s * n
                BndExp[b]->NormVectorIProductWRTBase(f_s,
                    BndExp[b]->UpdatePhys());
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
            if (!m_filePhi)
            {
                // Update 'm_phi' only if it was provided as a function
                m_phiEvaluator->Evaluate("Phi", m_phi->UpdatePhys(), t);
            }

            // Initialize both variables for the first step
            m_phiEvaluator->Evaluate(m_velName, m_up, t);

            // Initialise 'm_upPrev' in all cases
            m_upPrev = m_up;
        }
        // If timedependent 'm_phi'
        // Phi functions from files are not timedependent
        else if (m_timeDependentPhi)
        {
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
            
            else
            {
              //Fourier Interpolation
              npoints = m_phi -> GetTotPoints(); 
              //nCoeffs = m_phi -> GetNcoeffs();

              Array< OneD, NekDouble > temp;
              temp = Array<OneD, NekDouble>(npoints, 0.0);
              //Array< OneD, NekDouble > tempCoeffs;
              //tempCoeffs= Array<OneD, NekDouble>(nCoeffs, 0.0);

              NekDouble m_period = M_PI;
              NekDouble BetaT = 2*M_PI*fmod (t, m_period) / m_period;
              NekDouble phase;

              
              Vmath::Vcopy(npoints,&m_phiInterp[0],1,&temp[0],1);
              Vmath::Svtvp(npoints,cos(0.5*nSamples*BetaT),&m_phiInterp[npoints],1,&temp[0],1,&temp[0],1);
              Vmath::Svtvp(npoints,-sin(BetaT),&m_phiInterp[npoints],1,&temp[0],1,&temp[0],1);
              //Vmath::Vcopy(nCoeffs,&m_phiInterpCoeffs[0],1,&tempCoeffs[0],1);             
              //Vmath::Svtvp(nCoeffs,cos(0.5*nSamples*BetaT),&m_phiInterpCoeffs[nCoeffs],1,&tempCoeffs[0],1,&tempCoeffs[0],1);
              cout << nSamples << endl;
              for (int i = 2; i < nSamples; i += 2)
              {
                  phase = (i>>1) * BetaT;
                  Vmath::Svtvp(npoints, cos(phase),&m_phiInterp[i*npoints],1,&temp[0],1,&temp[0],1);
                  Vmath::Svtvp(npoints, -sin(phase),&m_phiInterp[(i+1)*npoints], 1, &temp[0], 1,&temp[0],1);
                  //Vmath::Svtvp(nCoeffs, cos(phase),&m_phiInterpCoeffs[i*nCoeffs],1,&tempCoeffs[0],1,&tempCoeffs[0],1);
                  //Vmath::Svtvp(nCoeffs, -sin(phase),&m_phiInterpCoeffs[(i+1)*nCoeffs], 1, &tempCoeffs[0], 1,&tempCoeffs[0],1);
              }

              /*for (int i = 1; i < nSamples; i += 1)
              {
                  phase = i * BetaT;
                  //Vmath::Svtvp(npoints, cos(phase),&m_phiInterp[i*npoints],1,&temp[0],1,&temp[0],1);
                  //Vmath::Svtvp(npoints, -sin(phase),&m_phiInterp[(i+1)*npoints], 1, &temp[0], 1,&temp[0],1);
                  Vmath::Svtvp(nCoeffs, cos(phase),&m_phiInterpCoeffs[i*nCoeffs],1,&tempCoeffs[0],1,&tempCoeffs[0],1);
                  //Vmath::Svtvp(nCoeffs, sin(phase),&m_phiInterpCoeffs[(i+1)*nCoeffs], 1, &tempCoeffs[0], 1,&tempCoeffs[0],1);
              }*/
              
              //Copy the coefficients to m_phi
              for (int i = 0; i < npoints; i+=1)
              {
                m_phi -> SetPhys(i,temp[i]);
              }

              /*for (int i = 0; i < nCoeffs; i+=1)
              {
                m_phi -> SetCoeffs(i,tempCoeffs[i]);
              }*/

              //Output Phi field
              m_phi->FwdTrans_IterPerExp(m_phi->GetPhys(), m_phi->UpdateCoeffs());
              m_phi->BwdTrans_IterPerExp(m_phi->GetCoeffs(), m_phi->UpdatePhys());
              //m_phi->BwdTrans_IterPerExp(m_phi->GetCoeffs(), m_phi->UpdatePhys());
              //m_phi->FwdTrans_IterPerExp(m_phi->GetPhys(), m_phi->UpdateCoeffs());
              std::vector<Array<OneD, NekDouble> > phiOutputVectorOfArray;
              phiOutputVectorOfArray.push_back(m_phi -> UpdateCoeffs()); 
              std::vector<std::string> variableName;
              variableName.push_back("phi");
              string tstring = to_string(t);
              EquationSystem::WriteFld("phi_"+ tstring +".fld",m_phi,phiOutputVectorOfArray, variableName); 
              //cout <<"After output tranformation"<< endl;

              //cout << "(m_phi -> UpdateCoeffs())[10] = " << (m_phi -> UpdateCoeffs())[10] << endl;
              //cout << "tempCoeffs[10] = " << tempCoeffs[10] << endl;

              //cout << "(m_phi -> UpdatePhys())[10] = " <<(m_phi -> UpdatePhys())[10] << endl;
              //cout << "temp[10] = " << temp[10] << endl;
            }            
        
        }
    }

    /**
     * @brief For a body with a velocity \f$\mathbf{u_p}\f$, the force
     * \f$\mathbf{f_s}\f$ applied to the fluid ensures that the IBC are met:
     *
     * \f[ \mathbf{f_s} = \frac{\Phi^{n+1}\left(\mathbf{u_p}^{n+1} -
     * \mathbf{u^*}\right)}{\Delta t} \f]
     *
     * @param fields
     * @param dt
     * @param f_s
     */
    void SmoothedProfileMethod::UpdateForcing(
                    const Array<OneD, const Array<OneD, NekDouble> > &fields,
                    NekDouble dt)
    {
        int nvel = m_velocity.size();
        int nq   = m_phi->GetNpoints();

        for (int i = 0; i < nvel; ++i)
        {
            // In homogeneous cases, switch out of wave space
            Array<OneD, NekDouble> tmpField(nq);
            int ind = m_velocity[i];

            if (m_HomogeneousType != EquationSystem::eNotHomogeneous &&
                m_fields[ind]->GetWaveSpace())
            {
                m_fields[ind]->HomogeneousBwdTrans(fields[i], tmpField);
                m_fs[i]->HomogeneousBwdTrans(m_fs[i]->GetPhys(),
                                             m_fs[i]->UpdatePhys());
            }
            else
            {
                tmpField = fields[i];
            }

            Vmath::Vsub(nq, m_up[i], 1, tmpField, 1, m_fs[i]->UpdatePhys(), 1);
            Vmath::Vmul(nq, m_phi->GetPhys(), 1, m_fs[i]->GetPhys(), 1,
                        m_fs[i]->UpdatePhys(), 1);
            Vmath::Smul(nq, m_gamma0/dt, m_fs[i]->GetPhys(), 1,
                        m_fs[i]->UpdatePhys(), 1);

            // And go back to wave space if the 'if' was executed
            if (m_HomogeneousType != EquationSystem::eNotHomogeneous &&
                m_fields[ind]->GetWaveSpace())
            {
                m_fs[i]->HomogeneousFwdTrans(m_fs[i]->GetPhys(),
                                             m_fs[i]->UpdatePhys());
            }
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
        SPMTest();
        // Function evaluator for Phi and Up
        m_phiEvaluator = GetFunction("ShapeFunction");

        TiXmlElement *function = GetFunctionHdl("ShapeFunction");
        TiXmlElement *child    = function->FirstChildElement();
        m_filePhi = false;
        m_timeDependentPhi = false;

        // If defined by using a file
        if (boost::iequals(child->ValueStr(), "F"))
        {
            // Get name of STL file
            string fileName;
            int status = child->QueryStringAttribute("FILE", &fileName);
            ASSERTL0(status == TIXML_SUCCESS, "A folder or or an FLD file with the values "
                     "of the phi function has to be supplied.")
          
            // Check if Phi is timedependent
            m_timeDependentPhi = GetVarTimeDependence("ShapeFunction", "Phi");

            if(m_timeDependentPhi)
            {             
                //Read number of samples
                TiXmlAttribute *childAttr = child->LastAttribute();
                std::string attrName(childAttr->Name());
                ASSERTL0(attrName == "NSAMPLES", "Unable to read attribute number of samples.");

                status = childAttr->QueryIntValue(&nSamples);
                ASSERTL0(status == TIXML_SUCCESS, "The number of samples "
                     "has to be specified.")

                // Import the STL samples into auxiliary vector

                // The STL samples should be stored in the form "%d.stl" where %d is
                // the rotation angle from the initial position

                // A subdirectory can also be included, such as "dir/%d.stl"                
                int angle = 0;
                int dAngle = 360/nSamples;
                ASSERTL0(typeid(dAngle) == typeid(int),"Angular postion of samples must be integers");
                
                string sampleFileName;
                string angleString;

                //Extract the value of npoints
                int npoints;
                //int nCoeffs;
                std::vector<LibUtilities::FieldDefinitionsSharedPtr> fieldDef;
                std::vector<std::vector<NekDouble> > fieldData;
                LibUtilities::FieldMetaDataMap fieldMetaData;
                sampleFileName = "./" + fileName + "/" + "0.fld";
                LibUtilities::FieldIOSharedPtr phiFile =
                    LibUtilities::FieldIO::CreateForFile(m_session, sampleFileName);
                phiFile->Import(sampleFileName, fieldDef, fieldData, fieldMetaData);

                // Extract Phi field to output
                string tmp("phi");
                m_phi->ExtractDataToCoeffs(fieldDef[0], fieldData[0],
                                            tmp, m_phi->UpdateCoeffs());
                m_phi->BwdTrans(m_phi->GetCoeffs(), m_phi->UpdatePhys());
                    
                npoints = m_phi -> GetTotPoints();
                //nCoeffs = m_phi -> GetNcoeffs();

                // Table containing phi fields of all samples 
                m_phiInterp = Array<OneD, NekDouble> (npoints*nSamples, 0.0);
                //m_phiInterpCoeffs = Array<OneD, NekDouble> (nCoeffs*nSamples, 0.0);

                for (int i = 0; i < nSamples; ++i)
                {
                    // Get phi values from XML file (after "FieldConvert" the STL file)
                    // First, load the data
                    std::vector<LibUtilities::FieldDefinitionsSharedPtr> fieldDef;
                    std::vector<std::vector<NekDouble> > fieldData;
                    LibUtilities::FieldMetaDataMap fieldMetaData;
                    angleString = to_string(angle);
                    
                    cout << "Imported Phi field for sample angle " + angleString + " degrees" << endl;
                    
                    sampleFileName = "./" + fileName + "/" + angleString + ".fld";
                    LibUtilities::FieldIOSharedPtr phiFile =
                        LibUtilities::FieldIO::CreateForFile(m_session, sampleFileName);
                    phiFile->Import(sampleFileName, fieldDef, fieldData, fieldMetaData);

                    // Only Phi field should be defined in the file
                    ASSERTL0(fieldData.size() == 1, "Only one field (phi) must be "
                                                    "defined in the FLD file.")

                    // Extract Phi field to output
                    string tmp("phi");
                    m_phi->ExtractDataToCoeffs(fieldDef[0], fieldData[0],
                                               tmp, m_phi->UpdateCoeffs());
                    m_phi->BwdTrans(m_phi->GetCoeffs(), m_phi->UpdatePhys());
                    
                    // Store Phi in table of Phi fields
                    Vmath::Vcopy(npoints,&((m_phi -> UpdatePhys())[0]),1,&m_phiInterp[i*npoints],1);
                    //Vmath::Vcopy(nCoeffs,&((m_phi -> UpdateCoeffs())[0]),1,&m_phiInterpCoeffs[i*nCoeffs],1);

                    //Update angle name
                    /*if(angle == 0)
                    {
                        angle = angle+dAngle;
                    }
                    else if(angle < 0)
                    {
                        angle = -angle + 10;
                    }
                    else
                    {
                        angle = -angle;
                    }*/

                    angle = angle + dAngle;
                }

                sampleFileName = "./" + fileName + "/0.fld";
                phiFile = LibUtilities::FieldIO::CreateForFile(m_session, sampleFileName);
                phiFile->Import(sampleFileName, fieldDef, fieldData, fieldMetaData);

                // Extract Phi field to output
                m_phi->ExtractDataToCoeffs(fieldDef[0], fieldData[0],
                                               tmp, m_phi->UpdateCoeffs());
                m_phi->BwdTrans(m_phi->GetCoeffs(), m_phi->UpdatePhys());
                
                //ASSERTL0(m_phiExp.size() == nSamples, "The number of samples provided does not match"
                //                                        " NSAMPLES.")

//#ifdef NEKTAR_USING_FFTW              

                //Discrete Fourier Transform using FFTW
                Array<OneD, NekDouble> fft_in(npoints*nSamples);
                Array<OneD, NekDouble> fft_out(npoints*nSamples);

                //Array<OneD, NekDouble> fft_inCoeffs(nCoeffs*nSamples);
                //Array<OneD, NekDouble> fft_outCoeffs(nCoeffs*nSamples);

                Array<OneD, NekDouble> m_tmpIN(nSamples);
                Array<OneD, NekDouble> m_tmpOUT(nSamples);

                //Array<OneD, NekDouble> m_tmpINCoeffs(nCoeffs);
                //Array<OneD, NekDouble> m_tmpOUTCoeffs(nCoeffs);

                //Shuffle the data
                for(int j= 0; j < nSamples; ++j)
                {
                    Vmath::Vcopy(npoints,&m_phiInterp[j*npoints],1,&(fft_in[j]),nSamples);
                    //Vmath::Vcopy(nCoeffs,&m_phiInterpCoeffs[j*nCoeffs],1,&(fft_inCoeffs[j*nCoeffs]),1);
                }

                m_FFT = LibUtilities::GetNektarFFTFactory().CreateInstance("NekFFTW", nSamples);
                //m_FFTCoeffs = LibUtilities::GetNektarFFTFactory().CreateInstance("NekFFTW", nSamples);

                //FFT Transform
                for(int i=0; i<npoints; i++)
                {
                    m_FFT->FFTFwdTrans(m_tmpIN =fft_in + i*nSamples, m_tmpOUT =fft_out + i*nSamples);

                }

                //for(int i=0; i<nCoeffs; i++)
                //{
                //    m_FFTCoeffs->FFTFwdTrans(m_tmpINCoeffs =fft_inCoeffs + i*nSamples, m_tmpOUTCoeffs =fft_outCoeffs + i*nSamples);
                //}

                //Reshuffle data
                for(int s = 0; s < nSamples; ++s)
                {
                    Vmath::Vcopy(npoints,&fft_out[s],nSamples,&m_phiInterp[s*npoints],1);
                    //Vmath::Vcopy(nCoeffs,&fft_outCoeffs[s*nCoeffs],1,&m_phiInterpCoeffs[s*nCoeffs],1);
                }

                Vmath::Zero(fft_in.size(),&fft_in[0],1);
                Vmath::Zero(fft_out.size(),&fft_out[0],1);

                //Vmath::Zero(fft_inCoeffs.size(),&fft_inCoeffs[0],1);
                //Vmath::Zero(fft_outCoeffs.size(),&fft_outCoeffs[0],1);

                //scaling of the Fourier coefficients
                NekDouble j=-1;
                for (int i = 2; i < nSamples; i += 2)
                {
                    Vmath::Smul(2*npoints,j,&m_phiInterp[i*npoints],1,&m_phiInterp[i*npoints],1);
                    //Vmath::Smul(2*nCoeffs,j,&m_phiInterpCoeffs[i*nCoeffs],1,&m_phiInterpCoeffs[i*nCoeffs],1);
                    j=-j;

                }
                
                m_timeDependentPhi = true;
                m_timeDependentUp  = false;

              //Output Phi field
              m_phi->FwdTrans_IterPerExp(m_phi->GetPhys(), m_phi->UpdateCoeffs());
              m_phi->BwdTrans_IterPerExp(m_phi->GetCoeffs(), m_phi->UpdatePhys());
              
              //m_phi->BwdTrans_IterPerExp(m_phi->GetCoeffs(), m_phi->UpdatePhys());
              //m_phi->FwdTrans_IterPerExp(m_phi->GetPhys(), m_phi->UpdateCoeffs());
              
              std::vector<Array<OneD, NekDouble> > phiOutputVectorOfArray;
              phiOutputVectorOfArray.push_back(m_phi -> UpdateCoeffs()); 
              std::vector<std::string> variableName;
              variableName.push_back("phi");
              EquationSystem::WriteFld("phi_0.fld",m_phi,phiOutputVectorOfArray, variableName); 

            }
            else
            {
                ASSERTL0(boost::iequals(fileName.substr(fileName.length()-4),
                         ".fld"), "A valid FLD file must be supplied in the "
                         "'ShapeFunction' field.")

                // Get phi values from XML file (after "FieldConvert" the STL file)
                // First, load the data
                std::vector<LibUtilities::FieldDefinitionsSharedPtr> fieldDef;
                std::vector<std::vector<NekDouble> > fieldData;
                LibUtilities::FieldMetaDataMap fieldMetaData;
                LibUtilities::FieldIOSharedPtr phiFile =
                    LibUtilities::FieldIO::CreateForFile(m_session, fileName);
                phiFile->Import(fileName, fieldDef, fieldData, fieldMetaData);

                // Only Phi field should be defined in the file
                ASSERTL0(fieldData.size() == 1, "Only one field (phi) must be "
                                                "defined in the FLD file.")

                // Extract Phi field to output
                string tmp("phi");
                m_phi->ExtractDataToCoeffs(fieldDef[0], fieldData[0],
                                           tmp, m_phi->UpdateCoeffs());
                m_phi->BwdTrans(m_phi->GetCoeffs(), m_phi->UpdatePhys());
                m_timeDependentPhi = false;
                m_timeDependentUp  = false;

                //Output Phi field
                m_phi->FwdTrans_IterPerExp(m_phi->GetPhys(), m_phi->UpdateCoeffs());
                std::vector<Array<OneD, NekDouble> > phiOutputVectorOfArray;
                phiOutputVectorOfArray.push_back(m_phi -> UpdateCoeffs()); 
                std::vector<std::string> variableName;
                variableName.push_back("phi");
                EquationSystem::WriteFld("phi.fld",m_phi,phiOutputVectorOfArray, variableName); 
                m_phi->BwdTrans_IterPerExp(m_phi->GetCoeffs(), m_phi->UpdatePhys());
                
            }
            m_filePhi = true;
        }
        else
        {
            // Check if Phi is timedependent
            m_timeDependentPhi = GetVarTimeDependence("ShapeFunction", "Phi");

            // If so, check if its velocity changes as well
            m_timeDependentUp = GetVarTimeDependence("ShapeFunction", "Up");
            switch (m_velocity.size())
            {
                case 3:
                    m_timeDependentUp |=
                        GetVarTimeDependence("ShapeFunction", "Wp");
                case 2:
                    m_timeDependentUp |=
                        GetVarTimeDependence("ShapeFunction", "Vp");
            }
        }
    }

    void SmoothedProfileMethod::SPMTest()
    {
        int nX = 10;
        int nSamplesD = 10;

        Array<OneD, NekDouble> Xspace;
        Xspace = Array<OneD, NekDouble>(nX);
        /*Array<OneD, NekDouble> func;
        func = Array<OneD, NekDouble>(nX, 0.0);*/
        const Array<OneD, NekDouble> func;
        //func = Array<OneD, NekDouble>(nX, 0.0);
        Array<OneD, NekDouble> funcInterp;
        funcInterp = Array<OneD, NekDouble>(nX*nSamplesD);

        for(int i=0; i<nX; ++i)
        {
            Xspace[i] = double(i)/double(nX);
        }

        ofstream myfile;
        myfile.open ("output_base.csv");

        NekDouble mu;
        NekDouble sigma = 0.1; 
        for(int i=0; i<nSamplesD; ++i)
        {   
            mu = double(i)/double(nSamplesD);
            for(int j=0; j<nX; ++j)
            {
                funcInterp[i*nX+j] = (sigma*2.51)*exp(-(Xspace[j]-mu)*(Xspace[j]-mu)/(2*sigma*sigma));
                if(i==0)
                {
                    myfile << Xspace[j] <<","<< funcInterp[i*nX+j] << "\n"; 
                }
                
            }        
        }

        myfile.close();

        Array<OneD, NekDouble> fft_inD(nX*nSamplesD);
        Array<OneD, NekDouble> fft_outD(nX*nSamplesD);

        Array<OneD, NekDouble> m_tmpIND(nSamplesD);
        Array<OneD, NekDouble> m_tmpOUTD(nSamplesD);

        /*for(int j= 0; j < nSamplesD; ++j)
        {
            Vmath::Vcopy(nX,&funcInterp[j*nX],1,&(fft_inD[j*nX]),1);
        }

        m_FFTD = LibUtilities::GetNektarFFTFactory().CreateInstance("NekFFTW", nSamplesD);

        //FFT Transform
        for(int i=0; i<nX; i++)
        {
            m_FFTD->FFTFwdTrans(m_tmpIND =fft_inD + i*nSamplesD, m_tmpOUTD =fft_outD + i*nSamplesD);
        }

        //scaling of the Fourier coefficients
        NekDouble j=-1;
        for (int i = 2; i < nSamplesD; i += 2)
        {
            Vmath::Smul(2*nX,j,&funcInterp[i*nX],1,&funcInterp[i*nX],1);
            j=-j;
        }


        //Reshuffle data
        for(int s = 0; s < nSamplesD; ++s)
        {
            Vmath::Vcopy(nSamplesD,&fft_outD[s*nX],1,&funcInterp[s*nX],1);
        }

        Vmath::Zero(fft_inD.size(),&fft_inD[0],1);
        Vmath::Zero(fft_outD.size(),&fft_outD[0],1);*/

        //Shuffle the data
        for(int j= 0; j < nSamplesD; ++j)
        {
            Vmath::Vcopy(nX,&funcInterp[j*nX],1,&(fft_inD[j]),nSamplesD);
        }

        m_FFTD = LibUtilities::GetNektarFFTFactory().CreateInstance("NekFFTW", nSamplesD);

        //FFT Transform
        for(int i=0; i<nX; i++)
        {
            m_FFTD->FFTFwdTrans(m_tmpIND =fft_inD + i*nSamplesD, m_tmpOUTD =fft_outD + i*nSamplesD);

        }

        //Reshuffle data
        for(int s = 0; s < nSamplesD; ++s)
        {
            Vmath::Vcopy(nX,&fft_outD[s],nSamplesD,&funcInterp[s*nX],1);

        }

        Vmath::Zero(fft_inD.size(),&fft_inD[0],1);
        Vmath::Zero(fft_outD.size(),&fft_outD[0],1);

        //scaling of the Fourier coefficients
        /*NekDouble j=-1;
        for (int i = 2; i < nSamplesD; i += 2)
        {
            Vmath::Smul(2*nX,j,&funcInterp[i*nX],1,&funcInterp[i*nX],1);
            j=-j;
        }*/

        Array< OneD, NekDouble > tempD;
        tempD = Array<OneD, NekDouble>(nX, 0.0);

        NekDouble m_period = 1.0;
        NekDouble phase;
        NekDouble BetaT;
        NekDouble t;

        for(int i=0; i<nX*nSamplesD; ++i)
        {
            cout << "funcInterp[" << i << "] = " << funcInterp[i] << endl;
        }
        
        for (int j=0; j<21; ++j)
        {   
            t = 0.1*j;
            BetaT = 2*M_PI*fmod (t, m_period) / m_period;
            
            Vmath::Vcopy(nX,&funcInterp[0],1,&tempD[0],1);
            Vmath::Svtvp(nX,cos(0.5*nSamplesD*BetaT),&funcInterp[nX],1,&tempD[0],1,&tempD[0],1);
            //Vmath::Svtvp(nX,-sin(BetaT),&funcInterp[nX],1,&tempD[0],1,&tempD[0],1);
            for (int i = 2; i < nSamplesD; i += 2)
            {
                phase = (i>>1) * BetaT;
                Vmath::Svtvp(nX, cos(phase),&funcInterp[i*nX],1,&tempD[0],1,&tempD[0],1);
                Vmath::Svtvp(nX, -sin(phase),&funcInterp[(i+1)*nX], 1, &tempD[0], 1,&tempD[0],1);
            }
            
            string tstring = to_string(t);
            myfile.open ("output_t" + tstring + ".csv");
            for(int i=0; i<nX; ++i)
            {
                myfile << Xspace[i] <<","<< tempD[i] << "\n";
            }
            
            myfile.close();
        }
            

        cout << "1D DONE" << endl;
    }


} // end of namespace
