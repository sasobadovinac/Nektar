#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Advection/Advection.h>
#include <LocalRegions/MatrixKey.h>

namespace Nektar
{
class LaxFriedrichsSolver : public SolverUtils::RiemannSolver
{
public:
    static SolverUtils::RiemannSolverSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    {
        return SolverUtils::RiemannSolverSharedPtr(
            new LaxFriedrichsSolver(pSession));
    }

    static std::string solverName;

protected:
    LaxFriedrichsSolver(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : SolverUtils::RiemannSolver(pSession)
    {
    }

    void v_Solve(
        const int                                         nDim,
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
              Array<OneD,       Array<OneD, NekDouble> > &flux) final
    {
        boost::ignore_unused(nDim);
        const Array<OneD, NekDouble> &traceVel = m_scalars["Vn"]();

        for (int j = 0; j < traceVel.size(); ++j)
        {
            flux[0][j] = 0.5 * traceVel[j] * (Fwd[0][j] + Bwd[0][j]) +
                0.5 * abs(traceVel[j]) * (Fwd[0][j] - Bwd[0][j]);
        }
    }
};
std::string LaxFriedrichsSolver::solverName = SolverUtils::GetRiemannSolverFactory().
    RegisterCreatorFunction("LaxFriedrichs", LaxFriedrichsSolver::create,
                            "L-F solver");


    class ALEDemo : public SolverUtils::EquationSystem
    {
    public:
        /// Class may only be instantiated through the MemoryManager.
        friend class MemoryManager<ALEDemo>;

        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph)
        {
            SolverUtils::EquationSystemSharedPtr p =
                MemoryManager<ALEDemo>::AllocateSharedPtr(pSession, pGraph);
            p->InitObject();
            return p;
        }

        /// Name of class
        static std::string className;

        virtual ~ALEDemo() = default;

    protected:
        SolverUtils::RiemannSolverSharedPtr m_riemannSolver;
        Array<OneD, Array<OneD, NekDouble>> m_velocity;
        Array<OneD, Array<OneD, NekDouble>> m_gridVelocity;
        Array<OneD, NekDouble> m_traceVn;
        SolverUtils::AdvectionSharedPtr m_advObject;

        ALEDemo(const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &pGraph)
            : EquationSystem(pSession, pGraph)
        {
        }

        void v_InitObject()
        {
            EquationSystem::v_InitObject();

            std::vector<std::string> vel{"Vx", "Vy", "Vz"};

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

            m_advObject = SolverUtils::GetAdvectionFactory().CreateInstance(
                "WeakDG", "WeakDG");
            m_advObject->SetFluxVector(&ALEDemo::GetFluxVector, this);

            m_riemannSolver =
                SolverUtils::GetRiemannSolverFactory().CreateInstance(
                    "Upwind", m_session);
            m_riemannSolver->SetScalar("Vn", &ALEDemo::GetNormalVelocity, this);
            m_riemannSolver->SetVector("N", &ALEDemo::GetNormals, this);
            m_advObject->SetRiemannSolver(m_riemannSolver);
            m_advObject->InitObject(m_session, m_fields);
        }

        virtual void v_DoSolve()
        {
            NekDouble dt = m_session->GetParameter("TimeStep");
            int nSteps   = m_session->GetParameter("NumSteps");

            SetInitialConditions(0.0);
            SetBoundaryConditions(0.0);

            int nq      = m_fields[0]->GetNpoints();
            int nc      = m_fields[0]->GetNcoeffs();
            int nFields = m_fields.size();

            Array<OneD, Array<OneD, NekDouble>> un(nFields);
            Array<OneD, Array<OneD, NekDouble>> tmp(nFields);
            Array<OneD, Array<OneD, NekDouble>> tmpPhys(nFields);
            for (int i = 0; i < nFields; ++i)
            {
                un[i] = Array<OneD, NekDouble>(nc);
                Vmath::Vcopy(nc, m_fields[i]->GetCoeffs(), 1, un[i], 1);
                tmp[i]     = Array<OneD, NekDouble>(nc, 0.0);
                tmpPhys[i] = Array<OneD, NekDouble>(nq);
            }

            // Dump inital conditions
            Checkpoint_Output(m_nchk++);
            m_time = 0;

            for (int n = 0; n < nSteps; ++n)
            {
                SetBoundaryConditions(m_time);

                for (int i = 0; i < nFields; ++i)
                {
                    m_fields[i]->BwdTrans(un[i], tmpPhys[i]);
                }

                // Calculate M^n^{-1} f(u^n, t_n)
                m_advObject->AdvectCoeffs(nFields, m_fields, m_velocity,
                                          tmpPhys, tmp, m_time);

                // Compute FE step
                for (int i = 0; i < nFields; ++i)
                {
                    for (int j = 0; j < nc; ++j)
                    {
                        tmp[i][j] = un[i][j] - dt * tmp[i][j];
                    }

                    // Multiply by mass matrix M^n
                    MultiRegions::GlobalMatrixKey mkey(StdRegions::eMass);
                    m_fields[i]->GeneralMatrixOp_IterPerExp(mkey, tmp[i],
                                                            un[i]);
                }

                // Perform movement and reset the matrices
                for (auto &field : m_fields)
                {
                    field->GetMovement()->PerformMovement(m_time);
                    field->Reset();
                }

                // Why does this break everything for rotating?!?
                m_fields[0]->SetUpPhysNormals();
                m_fields[0]->GetTrace()->GetNormals(m_traceNormals);

                for (int i = 0; i < nFields; ++i)
                {
                    m_fields[i]->MultiplyByElmtInvMass(un[i], un[i]);
                    for (int j = 0; j < nc; ++j)
                    {
                        m_fields[i]->UpdateCoeffs()[j] = un[i][j];
                    }
                    m_fields[i]->BwdTrans(un[i], m_fields[i]->UpdatePhys());
                }

                // Timestep should be complete
                if (n % 10 == 0 && n != 0)
                {
                    Checkpoint_Output(m_nchk++);
                }

                m_time = (n + 1) * dt;
            }

            if (nSteps % 100 == 0)
            {
                Checkpoint_Output(m_nchk++);
            }
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
                    /*
                    for (int k = 0; k < nq; ++k)
                    {
                        flux[i][j][k] = physfield[i][k] * (
                            m_velocity[j][k] - m_gridVelocity[j][k]);
                    }
                    */
                }
            }
        }

        const Array<OneD, const Array<OneD, NekDouble>> &GetNormals()
        {
            return m_traceNormals;
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

            GetGridVelocity();

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

            auto zones = m_fields[0]->GetMovement()->GetZones();
            for (auto &zone : zones)
            {
                // If the zone domain is fixed then grid velocity is left at 0
                if (zone.second->GetMovementType() == SpatialDomains::MovementType::eFixed)
                {
                    continue;
                }
                else if (zone.second->GetMovementType() == SpatialDomains::MovementType::eRotate)
                {
                    SpatialDomains::ZoneRotateShPtr zoneRotate =
                        std::static_pointer_cast<
                            SpatialDomains::ZoneRotate>(
                            zone.second);
                    NekDouble angVel = zoneRotate->GetAngularVel();

                    auto ids = zone.second->GetElementIds();
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
                            m_gridVelocity[0][offset + i] =  angVel * -yc[i];
                            m_gridVelocity[1][offset + i] =  angVel * xc[i];
                        }
                    }
                }
                else if (zone.second->GetMovementType() == SpatialDomains::MovementType::eTranslate)
                {
                    std::cout << "translate get grid vel" << std::endl;
                    SpatialDomains::ZoneTranslateShPtr interfaceTranslate =
                        std::static_pointer_cast<
                            SpatialDomains::ZoneTranslate>(
                            zone.second);
                    std::vector<NekDouble> velocity =
                        interfaceTranslate->GetVel();

                    auto ids = zone.second->GetElementIds();
                    for (auto id : ids)
                    {
                        int indx       = elmtToExpId[id];
                        int offset     = m_fields[0]->GetPhys_Offset(indx);
                        auto expansion = (*exp)[indx];

                        int nq = expansion->GetTotPoints();
                        for (int i = 0; i < nq; ++i)
                        {
                            m_gridVelocity[0][offset + i] = 0;
                            m_gridVelocity[1][offset + i] = 0;
                        }
                    }
                }
                else if (zone.second->GetMovementType() == SpatialDomains::MovementType::ePrescribe)
                {
                    /*NekDouble Lx = 20, Ly = 20;         // Size of mesh
                    NekDouble nx = 1, ny = 1, nt = 1;   // Space and time period
                    NekDouble X0 = 0.5, Y0 = 0.5;       // Amplitude
                    NekDouble t0 = sqrt(5*5 + 5*5);     // Time domain*/

                    auto ids = zone.second->GetElementIds();
                    for (auto id : ids)
                    {
                        int indx       = elmtToExpId[id];
                        int offset     = m_fields[0]->GetPhys_Offset(indx);
                        auto expansion = (*exp)[indx];

                        int nq = expansion->GetTotPoints();
                        Array<OneD, NekDouble> xc(nq, 0.0), yc(nq, 0.0), zc(nq, 0.0);
                        expansion->GetCoords(xc, yc, zc);

                        for (int i = 0; i < nq; ++i)
                        {
                            //NekDouble distx, disty;
                            // @TODO: Put analytic solution in (partial derivative r/t time)
                            /*distx = X0 * sin((nt * 2 * M_PI * m_time) / t0)
                                                    * sin((nx * 2 * M_PI * xc[i]) / Lx)
                                                    * sin((ny * 2 * M_PI * yc[i]) / Ly);

                            disty = Y0 * sin((nt * 2 * M_PI * m_time) / t0)
                                                    * sin((nx * 2 * M_PI * xc[i]) / Lx)
                                                    * sin((ny * 2 * M_PI * yc[i]) / Ly);*/
                            /*
                            if (xc[i] < 1e-8 || fabs(xc[i] - 1) < 1e-8)
                            {
                                m_gridVelocity[0][offset + i] = 0; //2 * M_PI * cos(2 * M_PI * m_time) * xc[i] * (1 - xc[i]);
                                m_gridVelocity[1][offset + i] = 0;
                            }
                            else
                            {
                                m_gridVelocity[0][offset + i] = 0.5; //2 * M_PI * cos(2 * M_PI * m_time) * xc[i] * (1 - xc[i]);
                                m_gridVelocity[1][offset + i] = 0;
                            }
                            */

                            m_gridVelocity[0][offset + i] = 0.05 * 2 * M_PI * cos(2*M_PI*m_time) * sin(2*M_PI*xc[i]) * sin(2*M_PI*yc[i]);
                            m_gridVelocity[1][offset + i] = 0.05 * 2 * M_PI * cos(2*M_PI*m_time) * sin(2*M_PI*xc[i]) * sin(2*M_PI*yc[i]);
                            //m_gridVelocity[0][offset + i] = 0.0;
                            //m_gridVelocity[1][offset + i] = 0.0;
                        }
                    }
                }
            }

            // Print coords trace
            /*for (auto &traceExp : *m_fields[0]->GetTrace()->GetExp())
            {
                int nq = traceExp->GetTotPoints();
                Array<OneD, NekDouble> xc(nq, 0.0), yc(nq, 0.0), zc(nq, 0.0);
                traceExp->GetCoords(xc, yc, zc);
                if (traceExp->GetGeom()->GetGlobalID() == 115)
                {
                    std::cout << "Time = " << m_time << " | Element ID = " << traceExp->GetGeom()->GetGlobalID() << " | Coords = ";
                    for (int i = 0; i < nq; ++i)
                    {
                        std::cout << "(" << xc[i] << " " << yc[i] << ") ->";
                    }
                    std::cout << std::endl;
                }
            }*/

            return m_gridVelocity;
        }
    };

    std::string ALEDemo::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "ALEDemo", ALEDemo::create,
            "Arbitrary-Lagrangian-Eulerian advection system");
}
