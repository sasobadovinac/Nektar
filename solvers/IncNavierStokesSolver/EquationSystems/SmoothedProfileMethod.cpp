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
        
        // Call it at the begining to make sure 'm_phi' is defined
        CalcPhi(m_pressureP, 0.0, false);

        // DEBUG: 'm_up' equal to 0 at all times
        m_up = Array<OneD, NekDouble>(nvel, 0.0);

        // DEGUB: select m_gamma0 depending on IMEX order
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

        // SPM correction of velocity
        v_SetUpCorrectionPressure(outarray, m_F, time, a_iixDt);
        v_SolveCorrectionPressure(m_F[0]);
        v_SolveCorrectedVelocity(m_F, outarray, time, a_iixDt);

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
    void SmoothedProfileMethod::v_SetUpCorrectionPressure(
                    const Array<OneD, const Array<OneD, NekDouble> > &fields,
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    NekDouble time,
                    NekDouble aii_Dt)
    {
        int physTot = m_pressureP->GetTotPoints();
        int nvel    = m_velocity.num_elements();

        // DEBUG: Set boundary conditions
        v_SetCorrectionPressureBCs(time, aii_Dt);

        // Virtual force 'fs'
        Array<OneD, Array<OneD, NekDouble> > f_s;
        IBForcing(fields, time, aii_Dt, f_s);
        m_fields[0]->PhysDeriv(eX, f_s[0], Forcing[0]);

        // Using 'Forcing[1]' as storage
        for (int i = 1; i < nvel; ++i)
        {
            m_fields[i]->PhysDeriv(DirCartesianMap[i], f_s[i], Forcing[1]);
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
    void SmoothedProfileMethod::v_SolveCorrectionPressure(
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
    void SmoothedProfileMethod::v_SolveCorrectedVelocity(
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
                                  Forcing[m_velocity[0]],
                                  Forcing[m_velocity[1]]);
        }
        else
        {
            m_pressureP->PhysDeriv(m_pressureP->GetPhys(),
                                  Forcing[m_velocity[0]],
                                  Forcing[m_velocity[1]],
                                  Forcing[m_velocity[2]]);
        }

        // Virtual force 'fs'
        Array<OneD, Array<OneD, NekDouble> > f_s;
        IBForcing(fields, time, dt, f_s);

        // Velocity correction
        for (int i = 0; i < nvel; ++i)
        {
            int ind = m_velocity[i];
            Vmath::Vsub(physTot, f_s[ind], 1, Forcing[ind], 1, Forcing[ind], 1);
            Blas::Daxpy(physTot, dt/m_gamma0, Forcing[ind], 1, fields[i], 1);
        }

        // DEBUG: What about other convective fields?
    }

    /**
     * @brief DEBUG: Updates the BCs for boundaries with Dirichlet BCs in the
     * velocity:
     *
     * \f[ \frac{\partial p_p}{\partial\mathbf{n}} = \mathbf{f_s}\cdot\mathbf{n} \f]
     *
     * @param dt
     */
    void SmoothedProfileMethod::v_SetCorrectionPressureBCs(NekDouble time,
                                                           NekDouble dt)
    {
        Array<OneD, ExpListSharedPtr> BndExp;

        // Get the BC expansions
        BndExp = m_pressureP->GetBndCondExpansions();

        // For each boundary...
        for (int b = 0; b < BndExp.num_elements(); ++b)
        {
            // Calculate f_s values
            Array<OneD, Array<OneD, NekDouble> > f_s;
            IBForcingBC(b, BndExp[b], time, dt, f_s);

            // BC is f_s * n
            BndExp[b]->NormVectorIProductWRTBase(f_s, BndExp[b]->UpdatePhys());
        }
    }

    /**
     * @brief DEBUG: Calculates the values of the shape function
     *
     * @param expansion
     * @param t
     * @param phi
     */
    void SmoothedProfileMethod::CalcPhi(const ExpListSharedPtr &expansion,
                                        NekDouble t,
                                        bool timeDependent)
    {
        // Calculate only once if not time-dependent
        if (timeDependent || t <= 0.0)
        {
            int physTot = expansion->GetTotPoints();
            Array<OneD, NekDouble> coord_0(physTot);
            Array<OneD, NekDouble> coord_1(physTot);
            Array<OneD, NekDouble> coord_2(physTot);

            expansion->GetCoords(coord_0, coord_1, coord_2);
            for (int i = 0; i < physTot; ++i)
            {
                m_phi->UpdatePhys()[i] = -0.5 * (tanh(
                        ((coord_0[i]-5.0)*(coord_0[i]-5.0) +
                        (coord_1[i]-5.0)*(coord_1[i]-5.0) - 0.25)/0.047386)
                        - 1.0);
            }
        }
    }

    /**
     * @brief For a body with a constant velocity \f[\mathbf{u_p}\f], the force
     * \f[\mathbf{f_s}\f] applied to the fluid ensures that the IBC are met
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

        // Vector phi
        CalcPhi(m_pressureP, time, false);

        // Vector f_s
        f_s = Array<OneD, Array<OneD, NekDouble> >(nvel);
        for (int i = 0; i < nvel; ++i)
        {
            f_s[i] = Array<OneD, NekDouble>(nq);
        }

        for (int i = 0; i < nvel; ++i)
        {
            Vmath::Sadd(nq, -m_up[m_velocity[i]],
                        fields[m_velocity[i]], 1,
                        f_s[i], 1);
            Vmath::Vmul(nq, m_phi->GetPhys(), 1, f_s[i], 1, f_s[i], 1);
            Vmath::Smul(nq, -m_gamma0/dt, f_s[i], 1, f_s[i], 1);
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

        // Vector phi
        CalcPhi(BndExp, time, false);

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
            Vmath::Sadd(nq, -m_up[m_velocity[i]],
                        velExp->GetPhys(), 1,
                        f_s[i], 1);
            Vmath::Vmul(nq, m_phi->GetBndCondExpansions()[bndInd]->GetPhys(), 1, f_s[i], 1, f_s[i], 1);
            Vmath::Smul(nq, -m_gamma0/dt, f_s[i], 1, f_s[i], 1);
        }
    }

} // end of namespace
