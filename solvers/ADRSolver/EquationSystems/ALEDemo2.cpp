#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Advection/Advection.h>
#include <SolverUtils/EquationSystem.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeOperators.h>
#include <iomanip>

#include <StdRegions/StdQuadExp.h>

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
    Array<OneD, NekDouble> m_traceVn, m_xc, m_yc;
    SolverUtils::AdvectionSharedPtr m_advObject;
    int m_infosteps;
    std::map<int, SpatialDomains::PointGeom> m_pts, m_ptsCurve;
    SpatialDomains::CurveMap m_curveEdge, m_curveFace;

    /// Wrapper to the time integration scheme
    LibUtilities::TimeIntegrationSchemeSharedPtr    m_intScheme;
    /// The time integration scheme operators to use.
    LibUtilities::TimeIntegrationSchemeOperators    m_ode;

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

        // Initialise grid velocity as 0s
        m_gridVelocity = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_gridVelocity[i] = Array<OneD, NekDouble>(m_fields[0]->GetTotPoints(), 0.0);
        }

        // Storage of points from original mesh.
        int nq = m_fields[0]->GetNpoints();
        m_xc = Array<OneD, NekDouble>(nq);
        m_yc = Array<OneD, NekDouble>(nq);
        m_fields[0]->GetCoords(m_xc, m_yc);
    }

    virtual void v_DoInitialise()
    {
        // Store all points
        auto &ptsMap = m_graph->GetAllPointGeoms();
        for (auto &pt : ptsMap)
        {
            m_pts[pt.first] = *(pt.second);
        }

        // Storage of curve points from original mesh
        m_curveEdge = m_graph->GetCurvedEdges();
        m_curveFace = m_graph->GetCurvedFaces();
        int cnt = 0;
        for (auto &edge : m_curveEdge)
        {
            for (auto &point : edge.second->m_points)
            {
                m_ptsCurve[cnt++] = *(point);
            }
        }

        SetBoundaryConditions(m_time);
        SetInitialConditions(m_time);
        GetGridVelocity(m_time);
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
            auto &ptsMap = m_graph->GetAllPointGeoms();
            for (auto &pt : ptsMap)
            {
                Array<OneD, NekDouble> newLoc(3, 0.0);
                auto pnt = m_pts[pt.first];

                newLoc[0] = pnt(0) + 0.05 * sin(2*M_PI*time) * sin(2*M_PI*pnt(0)) * sin(2*M_PI*pnt(1));
                newLoc[1] = pnt(1) + 0.05 * sin(2*M_PI*time) * sin(2*M_PI*pnt(0)) * sin(2*M_PI*pnt(1));
                pt.second->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
            }

            int cnt = 0;
            for (auto &edge : m_curveEdge)
            {
                for (auto &pt : edge.second->m_points)
                {
                    Array<OneD, NekDouble> newLoc(3, 0.0);
                    auto pnt = m_ptsCurve[cnt++];

                    newLoc[0] = pnt(0) + 0.05 * sin(2*M_PI*time) * sin(2*M_PI*pnt(0)) * sin(2*M_PI*pnt(1));
                    newLoc[1] = pnt(1) + 0.05 * sin(2*M_PI*time) * sin(2*M_PI*pnt(0)) * sin(2*M_PI*pnt(1));
                    pt->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
                }
            }

            // Attempt to reset element
            auto &quads = m_graph->GetAllQuadGeoms();
            for (auto &q : quads)
            {
                q.second->Reset(m_curveEdge, m_curveFace);
            }
            auto &tris = m_graph->GetAllTriGeoms();
            for (auto &t : tris)
            {
                t.second->Reset(m_curveEdge, m_curveFace);
            }

            for (int i = 0; i < m_fields.size(); ++i)
            {
                //m_fields[i]->GetInterfaces()->PerformMovement(time);
                m_fields[i]->Reset();
            }

            // Why does this break everything for rotating?!?
            m_fields[0]->SetUpPhysNormals();
            m_fields[0]->GetTrace()->GetNormals(m_traceNormals);

            // Recompute grid velocity.
            GetGridVelocity(time);

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

        int nq = physfield[0].size();

        for (int i = 0; i < flux.size(); ++i)
        {
            for (int j = 0; j < flux[0].size(); ++j)
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
        int nTracePts = GetTraceNpoints();
        int nPts = m_velocity[0].size();

        // Auxiliary variable to compute the normal velocity
        Array<OneD, NekDouble> tmp(nPts), tmp2(nTracePts);

        // Reset the normal velocity
        Vmath::Zero(nTracePts, m_traceVn, 1);

        for (int i = 0; i < m_velocity.size(); ++i)
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

    const Array<OneD, const Array<OneD, NekDouble>> &GetGridVelocity(NekDouble time)
    {
        //std::cout << "RECOMPUTE GRID VELOCITY TIME = " << time << std::endl;
        boost::ignore_unused(time);

        // First order approximation...
        //int nq = m_fields[0]->GetNpoints();
        //
        // Get updated coordinates.
        //Array<OneD, NekDouble> xc(nq), yc(nq);
        //m_fields[0]->GetCoords(xc, yc);
        //
        // for (int i = 0; i < nq; ++i)
        // {
        //     //m_gridVelocity[0][i] = (xc[i] - m_xc[i]) / (time - m_prevStageTime);
        //     //m_gridVelocity[1][i] = (yc[i] - m_yc[i]) / (time - m_prevStageTime);
        //     //m_xc[i] = xc[i];
        //     //m_yc[i] = yc[i];
        // }

        // My method for GLL_LAGRANGE_SEM same order as geometry order
        auto exp = m_fields[0]->GetExp();
        int nexp = m_fields[0]->GetExpSize();
        for (int i = 0; i < nexp; ++i)
        {
            auto elmt = m_fields[0]->GetExp(i);
            int offset = m_fields[0]->GetPhys_Offset(i);
            auto expansion = (*exp)[i];

            int nq = expansion->GetTotPoints();
            Array<OneD, NekDouble> xc(nq, 0.0), yc(nq, 0.0), zc(nq, 0.0);
            expansion->GetCoords(xc, yc, zc);

            for (int j = 0; j < nq; ++j)
            {
                m_gridVelocity[0][offset + j] = 0.05 * 2 * M_PI * cos(2*M_PI*time) * sin(2*M_PI*xc[j]) * sin(2*M_PI*yc[j]);
                m_gridVelocity[1][offset + j] = 0.05 * 2 * M_PI * cos(2*M_PI*time) * sin(2*M_PI*xc[j]) * sin(2*M_PI*yc[j]);
            }
        }

        // Daves method for interpolating the grid velocity
        //int nq = m_fields[0]->GetNpoints();
        //
        // Get updated coordinates.
        //Array<OneD, NekDouble> xc(nq), yc(nq);
        //m_fields[0]->GetCoords(xc, yc);
        /*int nexp = m_fields[0]->GetExpSize();
        for (int i = 0; i < nexp; ++i)
        {
            // Construct a standard expansion for each quadrilateral to evaluate
            // grid velocity inside deformed element.
            auto elmt = m_fields[0]->GetExp(i);
            auto pkey = elmt->GetBasis(0)->GetPointsKey();
            LibUtilities::BasisKey bkey(LibUtilities::eModified_A, 2, pkey);
            StdRegions::StdQuadExp stdexp(bkey, bkey);

            // Grid velocity of each quadrilateral vertex.
            Array<OneD, NekDouble> tmpx(4), tmpy(4);

            for (int j = 0; j < 4; ++j)
            {
                auto vert = elmt->GetGeom()->GetVertex(j);

                // original vertex from t=0 mesh
                auto pnt  = m_pts[vert->GetGlobalID()];

                // x/y velocity for each vertex
                tmpx[j] = 0.05 * 2 * M_PI * cos(2*M_PI*time) * sin(2*M_PI*pnt.x()) * sin(2*M_PI*pnt.y());
                tmpy[j] = 0.05 * 2 * M_PI * cos(2*M_PI*time) * sin(2*M_PI*pnt.x()) * sin(2*M_PI*pnt.y());
            }

            // swap to match tensor product order of coefficients
            // vs. anti-clockwise order of vertices.
            std::swap(tmpx[2], tmpx[3]);
            std::swap(tmpy[2], tmpy[3]);

            // Evaluate at all points in the deformed (time t) high-order
            // element.
            Array<OneD, NekDouble> tmp;
            stdexp.BwdTrans(
                tmpx, tmp = m_gridVelocity[0] + m_fields[0]->GetPhys_Offset(i));
            stdexp.BwdTrans(
                tmpy, tmp = m_gridVelocity[1] + m_fields[0]->GetPhys_Offset(i));
        }*/

        // Update
        /*
        // Do I need to do this for each field? I think it is consistent?
        auto exp = m_fields[0]->GetExp();

        // Create map from element ID to expansion ID
        std::map<int, int> elmtToExpId;
        for (int i = (*exp).size() - 1; i >= 0; --i)
        {
            elmtToExpId[(*exp)[i]->GetGeom()->GetGlobalID()] = i;
        }

        //std::cout << "COMPUTE GRID TIME = " << time << std::endl;

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
                        m_gridVelocity[0][offset + i] = -angVel*yc[i];
                        m_gridVelocity[1][offset + i] = angVel*xc[i];
                        //std::cout << "Coordinate: (" << xc[i] << ", " << yc[i] << ") has velocity = "
                        //          << m_gridVelocity[0][offset + i] << ", " << m_gridVelocity[1][offset + i] << std::endl;
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
            else if (interface->GetInterfaceType() == SpatialDomains::ePrescribed)
            {
                // This is hacky - as interface is set up for 2 sides usually, we only use the left side in this case
                if (interface->GetSide() == SpatialDomains::eLeft)
                {
                    auto ids = interface->GetElementIds();
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
                            m_gridVelocity[0][offset + i] = 0.05 * 2 * M_PI * cos(2*M_PI*time) * sin(2*M_PI*xc[i]) * sin(2*M_PI*yc[i]);
                            m_gridVelocity[1][offset + i] = 0.05 * 2 * M_PI * cos(2*M_PI*time) * sin(2*M_PI*xc[i]) * sin(2*M_PI*yc[i]);
                            //m_gridVelocity[0][offset + i] = -yc[i];
                            //m_gridVelocity[1][offset + i] = xc[i];
                        }
                    }
                }
            }
        }
*/

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

    void v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables)
    {
        int nCoeffs = m_fields[0]->GetNcoeffs();
        Array<OneD, NekDouble> gridVelFwdX(nCoeffs, 0.0);
        Array<OneD, NekDouble> gridVelFwdY(nCoeffs, 0.0), tmp;

        m_fields[0]->FwdTrans_IterPerExp(m_gridVelocity[0], gridVelFwdX);
        m_fields[0]->FwdTrans_IterPerExp(m_gridVelocity[1], gridVelFwdY);


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
