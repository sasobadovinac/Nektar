#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Advection/Advection.h>
#include <SolverUtils/EquationSystem.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeOperators.h>
#include <iomanip>
#include <SolverUtils/ALEHelper.h>

#include <StdRegions/StdQuadExp.h>

using namespace std;

namespace Nektar
{

class ALEDemo2 : public SolverUtils::EquationSystem, public SolverUtils::ALEHelper
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
    Array<OneD, NekDouble> m_traceVn, m_xc, m_yc;
    SolverUtils::AdvectionSharedPtr m_advObject;
    int m_infosteps;
    std::map<int, SpatialDomains::PointGeom> m_pts;
    SpatialDomains::CurveMap m_curveEdge, m_curveFace;

    /// Wrapper to the time integration scheme
    LibUtilities::TimeIntegrationSchemeSharedPtr    m_intScheme;
    /// The time integration scheme operators to use.
    LibUtilities::TimeIntegrationSchemeOperators    m_ode;

    ALEDemo2(const LibUtilities::SessionReaderSharedPtr &pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph)
        : EquationSystem(pSession, pGraph),
          SolverUtils::ALEHelper()
    {
    }

    void v_InitObject()
    {
        EquationSystem::v_InitObject();
        ALEHelper::InitObject(m_graph, m_fields);

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
        UpdateGridVelocity(m_time);
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
        const int nvar = m_fields.size();
        const int nm = GetNcoeffs();

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

        m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);

        int step = 0;

        for (int i = 0; i < m_steps; ++i)
        {
            //std::cout << "SOLVING TIME = " << m_time << std::endl;
            fields = m_intScheme->TimeIntegrate(step, m_timestep, m_ode);
            //std::cout << "FINISHED SOLVING TIME" << std::endl;

            ++step;
            m_time += m_timestep;

            if (m_session->GetComm()->GetRank() == 0 &&
                !((step+1) % m_infosteps))
            {
                cout << "Steps: " << setw(8)  << left << step+1 << " "
                     << "Time: "  << setw(12) << left << m_time << endl;
            }

            m_fields[0]->Reset();
            m_fields[0]->GetTrace()->GetNormals(m_traceNormals);
            UpdateGridVelocity(m_time);
            SetBoundaryConditions(m_time);

            // Update m_fields with u^n by multiplying by inverse mass
            // matrix. That's then used in e.g. checkpoint output and L^2 error
            // calculation.
            //std::cout << "ELMT INV MASS TIME = " << m_time << std::endl;
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
        //std::cout << "EVAL RHS WITH TIME = " << time << std::endl;
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
        m_advObject->AdvectCoeffs(nVariables, m_fields, m_velocity, tmpin,
                                  outarray, m_time);

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
        //std::cout << "ODE PROJECT WITH TIME = " << time << std::endl;
        // Counter variable
        int i;

        // Number of fields (variables of the problem)
        int nVariables = inarray.size();

        // Get rid of multiplication by mass matrix
        Array<OneD, NekDouble> tmp(GetNcoeffs());

        // We need to update the geometry for the next stage, if necessary.
        if (time != m_prevStageTime)
        {
            //std::cout << "MOVING FIELDS TIME = " << time << std::endl;
            /*auto &ptsMap = m_graph->GetAllPointGeoms();
            for (auto &pt : ptsMap)
            {
                Array<OneD, NekDouble> coords(2, 0.0);
                pt.second->GetCoords(coords);

                Array<OneD, NekDouble> newLoc(3, 0.0);
                auto pnt = m_pts[pt.first];
                newLoc[0] = pnt(0) + 0.05 * sin(2*M_PI*time) * sin(2*M_PI*pnt(0)) * sin(2*M_PI*pnt(1));
                newLoc[1] = pnt(1) + 0.05 * sin(2*M_PI*time) * sin(2*M_PI*pnt(0)) * sin(2*M_PI*pnt(1));

                pt.second->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);

            }*/

            for (i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->GetMovement()->PerformMovement(time);
                m_fields[i]->Reset();
            }

            // Why does this break everything for rotating?!?
            m_fields[0]->SetUpPhysNormals();
            m_fields[0]->GetTrace()->GetNormals(m_traceNormals);

            // Recompute grid velocity.
            UpdateGridVelocity(m_time);

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

        for (i = 0; i < flux.size(); ++i)
        {
            for (j = 0; j < flux[0].size(); ++j)
            {
                // This is u * vel - u * gridvel
                //Vmath::Vvtvvtm(nq, physfield[i], 1, m_velocity[j], 1,
                //              physfield[i], 1, m_gridVelocity[j], 1,
                //               flux[i][j], 1);

                for (int k = 0; k < nq; ++k)
                {
                    flux[i][j][k] = physfield[i][k] * (
                        m_velocity[j][k] - m_gridVelocity[j][k]);
                }
            }
        }
    }

    virtual bool v_PostIntegrate(int step)
    {
        boost::ignore_unused(step);
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
        int nPts = m_velocity[0].size();

        // Auxiliary variable to compute the normal velocity
        Array<OneD, NekDouble> tmp(nPts), tmp2(nTracePts);

        // Reset the normal velocity
        Vmath::Zero(nTracePts, m_traceVn, 1);

        for (i = 0; i < m_velocity.size(); ++i)
        {
            // Subtract grid velocity here from velocity
            // velocity - grid velocity
            Vmath::Vsub(nPts, m_velocity[i], 1, m_gridVelocity[i], 1, tmp, 1);

            m_fields[0]->ExtractTracePhys(tmp, tmp2);

            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, tmp2, 1, m_traceVn,
                         1, m_traceVn, 1);
        }

        return m_traceVn;
    }

    void v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables)
    {
        int nCoeffs = m_fields[0]->GetNcoeffs();
        Array<OneD, NekDouble> gridVelFwdX(nCoeffs, 0.0);
        Array<OneD, NekDouble> gridVelFwdY(nCoeffs, 0.0), tmp;

        m_fields[0]->FwdTrans_IterPerExp(m_gridVelocity[0], gridVelFwdX);
        m_fields[0]->FwdTrans_IterPerExp(m_gridVelocity[1], gridVelFwdY);

        /*
        int nexp = m_fields[0]->GetExpSize();

        for (int i = 0; i < nexp; ++i)
        {
            auto exp = m_fields[0]->GetExp(i);
            auto offset_phys = m_fields[0]->GetPhys_Offset(i);
            auto offset_coeff = m_fields[0]->GetCoeff_Offset(i);

            auto stdexp = exp->GetStdExp();
            stdexp->FwdTrans(m_gridVelocity[0] + offset_phys,
                             tmp = gridVelFwdX + offset_coeff);
            stdexp->FwdTrans(m_gridVelocity[1] + offset_phys,
                             tmp = gridVelFwdY + offset_coeff);
        }
        */

        fieldcoeffs.push_back(gridVelFwdX);
        fieldcoeffs.push_back(gridVelFwdY);
        variables.push_back("gridVx");
        variables.push_back("gridVy");
    }
};

std::string ALEDemo2::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "ALEDemo2", ALEDemo2::create,
        "Arbitrary-Lagrangian-Eulerian advection system");
}