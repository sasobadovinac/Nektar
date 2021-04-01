#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Advection/Advection.h>

namespace Nektar
{
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
                    field->GetInterfaces()->PerformMovement(m_timestep);
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
                // Subtract grid velocity here from trace velocity
                // tmp - grid velocity
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
                            m_gridVelocity[0][offset + i] = angVel * -yc[i];
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
    };

    std::string ALEDemo::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "ALEDemo", ALEDemo::create,
            "Arbitrary-Lagrangian-Eulerian advection system");
}