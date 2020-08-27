#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Advection/Advection.h>
#include <iomanip>

using namespace std;

namespace Nektar
{

class ALEDemo2 : public SolverUtils::EquationSystem
{
public:
    /// Class may only be instantiated through the MemoryManager.
    friend class MemoryManager<ALEDemo2>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<ALEDemo2>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

    virtual ~ALEDemo2() = default;

protected:
    NekDouble m_prevStageTime = 0.0;
    SolverUtils::RiemannSolverSharedPtr m_riemannSolver;
    Array<OneD, Array<OneD, NekDouble>> m_velocity;
    Array<OneD, Array<OneD, NekDouble>> m_gridVelocity;
    Array<OneD, Array<OneD, NekDouble>> m_traceGridVelocity;
    Array<OneD, NekDouble> m_traceVn;
    SolverUtils::AdvectionSharedPtr m_advObject;
    int m_infosteps;

    /// Wrapper to the time integration scheme
    LibUtilities::TimeIntegrationSchemeSharedPtr    m_intScheme;
    /// The time integration scheme operators to use.
    LibUtilities::TimeIntegrationSchemeOperators    m_ode;
    ///
    LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr  m_intSoln;

    ALEDemo2(const LibUtilities::SessionReaderSharedPtr &pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph)
        : EquationSystem(pSession, pGraph)
    {
    }

    void v_InitObject()
    {
        EquationSystem::v_InitObject();

        std::vector<std::string> vel;
        vel.push_back("Vx");
        vel.push_back("Vy");
        vel.push_back("Vz");

        // Resize the advection velocities vector to dimension of the problem
        vel.resize(m_spacedim);

        // Store in the global variable m_velocity the advection velocities
        m_velocity = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        GetFunction("AdvectionVelocity")->Evaluate(vel, m_velocity);

        // Define the normal velocity fields
        if (m_fields[0]->GetTrace())
        {
            m_traceVn = Array<OneD, NekDouble>(GetTraceNpoints());
        }

        // Construct advection object
        m_advObject = SolverUtils::GetAdvectionFactory().CreateInstance(
            "WeakDG", "WeakDG");
        m_advObject->SetFluxVector(&ALEDemo2::GetFluxVector, this);

        // Construct Riemann solver
        m_riemannSolver =
            SolverUtils::GetRiemannSolverFactory().CreateInstance(
                "Upwind", m_session);
        m_riemannSolver->SetScalar("Vn", &ALEDemo2::GetNormalVelocity, this);
        m_advObject->SetRiemannSolver(m_riemannSolver);
        m_advObject->InitObject(m_session, m_fields);

        // Create time integration scheme
        std::string methodName =
            m_session->GetSolverInfo("TIMEINTEGRATIONMETHOD");
        m_intScheme = LibUtilities::GetTimeIntegrationSchemeFactory()
            .CreateInstance(methodName, "", 0, std::vector<NekDouble>());

        // Define RHS and projection methods
        m_ode.DefineOdeRhs     (&ALEDemo2::DoOdeRhs,        this);
        m_ode.DefineProjection (&ALEDemo2::DoOdeProjection, this);

        m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
    }

    virtual void v_DoInitialise()
    {
        SetBoundaryConditions(m_time);
        SetInitialConditions(m_time);
    }

    /*
     * A very simplified version of UnsteadySystem::v_DoSolve, which instead of
     * storing things in physical space stores them in modal space instead. It
     * also explicitly breaks the link between m_fields (i.e. the ExpList) and
     * the fields class that is used to store the solution at each timestep.
     *
     * It also pre-multiplies the initial conditions by the mass matrix, as this
     * class is designed to timestep the quantity (Mu) rather than u.
     */
    virtual void v_DoSolve()
    {
        const int nm = GetNcoeffs();
        const int nvar = m_fields.size();

        // Set initial condition
        MultiRegions::GlobalMatrixKey mkey(StdRegions::eMass);
        Array<OneD, Array<OneD, NekDouble>> fields(nvar);

        // Premultiply each field by the mass matrix
        for (int i = 0; i < nvar; ++i)
        {
            fields[i] = Array<OneD, NekDouble>(nm);
            m_fields[i]->GeneralMatrixOp_IterPerExp(
                mkey, m_fields[i]->GetCoeffs(), fields[i]);
        }

        m_intSoln = m_intScheme->InitializeScheme(
            m_timestep, fields, m_time, m_ode);

        int step = 0;

        for (int i = 0; i < m_steps; ++i)
        {
            fields = m_intScheme->TimeIntegrate(
                step, m_timestep, m_intSoln, m_ode);

            ++step;
            m_time += m_timestep;

            if (m_session->GetComm()->GetRank() == 0 &&
                !((step+1) % m_infosteps))
            {
                cout << "Steps: " << setw(8)  << left << step << " "
                     << "Time: "  << setw(12) << left << m_time << endl;
            }

            // Update m_fields with u^n by multiplying by inverse mass
            // matrix. That's then used in e.g. checkpoint output and L^2 error
            // calculation.
            for (int j = 0; j < nvar; ++j)
            {
                m_fields[j]->MultiplyByElmtInvMass(
                    fields[j], m_fields[j]->UpdateCoeffs());
                m_fields[j]->BwdTrans(
                    m_fields[j]->GetCoeffs(), m_fields[j]->UpdatePhys());
            }

            if ((m_checksteps && !((step + 1) % m_checksteps)))
            {
                Checkpoint_Output(m_nchk);
                m_nchk++;
            }
        }
    }

    /**
     * @brief Compute the right-hand side for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
              Array<OneD,        Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        boost::ignore_unused(time);
        const int nVariables = inarray.size();
        const int nc = GetNcoeffs();
        int i;

        // General idea is that we are time-integrating the quantity (Mu), so we
        // need to multiply input by inverse mass matrix to get coefficients u,
        // and then backwards transform so we can apply the DG operator.
        Array<OneD, NekDouble> tmp(nc);
        Array<OneD, Array<OneD, NekDouble>> tmpin(nVariables);

        for (i = 0; i < nVariables; ++i)
        {
            tmpin[i] = Array<OneD, NekDouble>(GetNpoints());
            m_fields[i]->MultiplyByElmtInvMass(inarray[i], tmp);
            m_fields[i]->BwdTrans(tmp, tmpin[i]);
        }

        // RHS computation using the new advection base class
        m_advObject->AdvectCoeffs(nVariables, m_fields, m_velocity,
                                  tmpin, outarray, m_time);

        // Negate the RHS. Note the output here is _not_ in physical space, but
        // in coefficient space.
        for (i = 0; i < nVariables; ++i)
        {
            Vmath::Neg(outarray[i].size(), outarray[i], 1);
        }
    }

    /**
     * @brief Compute the projection for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD,       Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        // Counter variable
        int i;

        // Number of fields (variables of the problem)
        int nVariables = inarray.size();

        // Get rid of multiplication by mass matrix
        Array<OneD, NekDouble> tmp(GetNcoeffs());

        // We need to update the geometry for the next stage, if necessary.
        if (time != m_prevStageTime)
        {
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->GetInterfaces()->PerformMovement(time - m_prevStageTime);
                m_fields[i]->Reset();
                m_fields[i]->GetInterfaces()->GenGeomFactors();
            }

            m_prevStageTime = time;
        }

        // For DG the projection is just a straightforward copy.
        for (i = 0; i < nVariables; ++i)
        {
            Vmath::Vcopy(inarray[i].size(), inarray[i], 1, outarray[i], 1);
        }

        // Set the boundary conditions
        SetBoundaryConditions(time);
    }

    /**
     * @brief Return the flux vector for the linear advection equation.
     *
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
    {
        ASSERTL1(flux[0].size() == m_velocity.size(),
                 "Dimension of flux array and velocity array do not match");

        int i, j;
        int nq = physfield[0].size();

        GetGridVelocity();

        for (i = 0; i < flux.size(); ++i)
        {
            for (j = 0; j < flux[0].size(); ++j)
            {
                // This is u * vel - u * gridvel
                Vmath::Vvtvvtm(nq, physfield[i], 1, m_velocity[j], 1,
                               physfield[i], 1, m_gridVelocity[j], 1,
                               flux[i][j], 1);
                //Vmath::Vmul(nq, physfield[i], 1, m_velocity[j], 1, flux[i][j], 1);
            }
        }

        // Here you're going to need something to add onto the flux the
        // contribution -vu
        //
        // v = grid velocity
        // u = physfield
    }

    virtual bool v_PostIntegrate(int step)
    {
        boost::ignore_unused(step);

        /*
        for (int i = 0; i < m_fields.size(); ++i)
        {
            m_fields[i]->MultiplyByElmtInvMass(m_fields[i]->GetCoeffs(),
                                               m_fields[i]->UpdateCoeffs());
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                  m_fields[i]->UpdatePhys());
        }
        */

        return false;
    }

    /**
     * @brief Get the normal velocity for the linear advection equation.
     */
    Array<OneD, NekDouble> &GetNormalVelocity()
    {
        // Number of trace (interface) points
        int i;
        int nTracePts = GetTraceNpoints();

        // Auxiliary variable to compute the normal velocity
        Array<OneD, NekDouble> tmp(nTracePts);

        // Reset the normal velocity
        Vmath::Zero(nTracePts, m_traceVn, 1);

        GetTraceGridVelocity();

        for (i = 0; i < m_velocity.size(); ++i)
        {
            m_fields[0]->ExtractTracePhys(m_velocity[i], tmp);

            // Subtract grid velocity here from trace velocity
            // tmp - grid velocity
            Vmath::Vsub(nTracePts, tmp, 1, m_traceGridVelocity[i], 1, tmp, 1);

            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, tmp, 1, m_traceVn,
                         1, m_traceVn, 1);
        }


        return m_traceVn;
    }

    const Array<OneD, const Array<OneD, NekDouble>> &GetGridVelocity()
    {
        // Initialise grid velocity as 0s
        m_gridVelocity = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_gridVelocity[i] = Array<OneD, NekDouble>(m_fields[0]->GetTotPoints(), 0.0);
        }

        // Do I need to do this for each field? I think it is consistent?
        auto exp = m_fields[0]->GetExp();

        // Create map from element ID to expansion ID
        std::map<int, int> elmtToExpId;
        for (int i = (*exp).size() - 1; i >= 0; --i)
        {
            elmtToExpId[(*exp)[i]->GetGeom()->GetGlobalID()] = i;
        }

        auto intVec = m_fields[0]->GetInterfaces()->GetInterfaceVector();
        for (auto &interface : intVec)
        {
            // If the interface domain is fixed then grid velocity is left at 0
            if (interface->GetInterfaceType() == SpatialDomains::eFixed)
            {
                continue;
            }
            else if (interface->GetInterfaceType() == SpatialDomains::eRotating)
            {
                SpatialDomains::RotatingInterfaceShPtr interfaceRotate =
                    std::static_pointer_cast<
                        SpatialDomains::RotatingInterface>(interface);
                NekDouble angVel = interfaceRotate->GetAngularVel();

                auto ids = interface->GetElementIds();
                for (auto id : ids)
                {
                    int indx       = elmtToExpId[id];
                    int offset     = m_fields[0]->GetPhys_Offset(indx);
                    auto expansion = (*exp)[indx];

                    int nq = expansion->GetTotPoints();
                    Array<OneD, NekDouble> xc(nq, 0.0), yc(nq, 0.0),
                        zc(nq, 0.0);
                    expansion->GetCoords(xc, yc, zc);

                    for (int i = 0; i < nq; ++i)
                    {
                        //std::cout << "Coordinate: (" << xc[i] << ", " << yc[i] << ") has velocity = " << -yc[i] << ", " << xc[i] << std::endl;
                        m_gridVelocity[0][offset + i] = -angVel * yc[i];
                        m_gridVelocity[1][offset + i] = angVel * xc[i];
                    }
                }
            }
            else if (interface->GetInterfaceType() == SpatialDomains::eSliding)
            {
                SpatialDomains::SlidingInterfaceShPtr interfaceSlide =
                    std::static_pointer_cast<
                        SpatialDomains::SlidingInterface>(interface);
                std::vector<NekDouble> velocity = interfaceSlide->GetVel();

                auto ids = interface->GetElementIds();
                for (auto id : ids)
                {
                    int indx       = elmtToExpId[id];
                    int offset     = m_fields[0]->GetPhys_Offset(indx);
                    auto expansion = (*exp)[indx];

                    int nq = expansion->GetTotPoints();
                    for (int i = 0; i < nq; ++i)
                    {
                        m_gridVelocity[0][offset + i] = velocity[0];
                        m_gridVelocity[1][offset + i] = velocity[1];
                    }
                }
            }
        }

        return m_gridVelocity;
    }

    const Array<OneD, const Array<OneD, NekDouble>> &GetTraceGridVelocity()
    {
        // Initialise grid velocity as 0s
        m_traceGridVelocity = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int i =0; i < m_spacedim; ++i)
        {
            m_traceGridVelocity[i] = Array<OneD, NekDouble>(m_fields[0]->GetTrace()->GetTotPoints(), 0.0);
        }

        // Do I need to do this for each field? I think it is consistent?
        auto exp = m_fields[0]->GetTrace()->GetExp();

        // Create map from element ID to expansion ID
        std::map<int, int> elmtToTraceId;
        for (int i = (*exp).size() - 1; i >= 0; --i)
        {
            elmtToTraceId[(*exp)[i]->GetGeom()->GetGlobalID()] = i;
        }

        auto intVec = m_fields[0]->GetInterfaces()->GetInterfaceVector();
        for (auto &interface : intVec)
        {
            // If the interface domain is fixed then grid velocity is left at 0
            if (interface->GetInterfaceType() == SpatialDomains::eFixed)
            {
                continue;
            }
            else if (interface->GetInterfaceType() == SpatialDomains::eRotating)
            {

                SpatialDomains::RotatingInterfaceShPtr interfaceRotate =
                    std::static_pointer_cast<
                        SpatialDomains::RotatingInterface>(interface);
                NekDouble angVel = interfaceRotate->GetAngularVel();

                auto ids = interface->GetElementIds();
                for (auto id : ids)
                {
                    int ne = m_graph->GetGeometry2D(id)->GetNumEdges();

                    for (int i = 0; i < ne; ++i)
                    {
                        auto edge = m_graph->GetGeometry2D(id)->GetEdge(i);
                        int indx  = elmtToTraceId[edge->GetGlobalID()];
                        int offset = m_fields[0]->GetTrace()->GetPhys_Offset(indx);
                        auto expansion = (*exp)[indx];

                        int nq = expansion->GetTotPoints();
                        Array<OneD, NekDouble> xc(nq, 0.0), yc(nq, 0.0), zc(nq, 0.0);
                        expansion->GetCoords(xc, yc, zc);

                        for (int j = 0; j < nq; ++j)
                        {
                            m_traceGridVelocity[0][offset + j] = -angVel * yc[j];
                            m_traceGridVelocity[1][offset + j] = angVel * xc[j];
                        }
                    }
                }
            }
            else if (interface->GetInterfaceType() == SpatialDomains::eSliding)
            {
                SpatialDomains::SlidingInterfaceShPtr interfaceSlide =
                    std::static_pointer_cast<
                        SpatialDomains::SlidingInterface>(interface);
                std::vector<NekDouble> velocity = interfaceSlide->GetVel();

                auto ids = interface->GetElementIds();
                for (auto id : ids)
                {
                    int ne = m_graph->GetGeometry2D(id)->GetNumEdges();
                    for (int i = 0; i < ne; ++i)
                    {
                        auto edge = m_graph->GetGeometry2D(id)->GetEdge(i);
                        int indx       = elmtToTraceId[edge->GetGlobalID()];
                        int offset = m_fields[0]->GetTrace()->GetPhys_Offset(indx);
                        auto expansion = (*exp)[indx];

                        int nq = expansion->GetTotPoints();
                        for (int j = 0; j < nq; ++j)
                        {
                            m_traceGridVelocity[0][offset + j] = velocity[0];
                            m_traceGridVelocity[1][offset + j] = velocity[1];
                        }
                    }
                }
            }
        }

        return m_traceGridVelocity;
    }
};

std::string ALEDemo2::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "ALEDemo2", ALEDemo2::create,
        "Arbitrary-Lagrangian-Eulerian advection system");
}
