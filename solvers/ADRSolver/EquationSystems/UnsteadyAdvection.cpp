/////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyAdvection.cpp
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
// Description: Unsteady linear advection solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <ADRSolver/EquationSystems/UnsteadyAdvection.h>

using namespace std;

namespace Nektar
{
    string UnsteadyAdvection::className = SolverUtils::GetEquationSystemFactory().
        RegisterCreatorFunction("UnsteadyAdvection",
                                UnsteadyAdvection::create,
                                "Unsteady Advection equation.");

    UnsteadyAdvection::UnsteadyAdvection(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          AdvectionSystem(pSession, pGraph)
    {
        m_planeNumber = 0;
    }

    /**
     * @brief Initialisation object for the unsteady linear advection equation.
     */
    void UnsteadyAdvection::v_InitObject()
    {
        // Call to the initialisation object of UnsteadySystem
        AdvectionSystem::v_InitObject();

        m_session->LoadParameter("wavefreq",   m_waveFreq, 0.0);
        // Read the advection velocities from session file

        std::vector<std::string> vel;
        vel.push_back("Vx");
        vel.push_back("Vy");
        vel.push_back("Vz");

        // Resize the advection velocities vector to dimension of the problem
        vel.resize(m_spacedim);

        // Store in the global variable m_velocity the advection velocities
        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        GetFunction( "AdvectionVelocity")->Evaluate(vel,  m_velocity);

        // Type of advection class to be used
        switch(m_projectionType)
        {
            // Continuous field
            case MultiRegions::eGalerkin:
            {
                string advName;
                m_session->LoadSolverInfo(
                    "AdvectionType", advName, "NonConservative");
                m_advObject = SolverUtils::
                    GetAdvectionFactory().CreateInstance(advName, advName);
                if (m_specHP_dealiasing)
                {
                    m_advObject->SetFluxVector(
                        &UnsteadyAdvection::GetFluxVectorDeAlias, this);
                }
                else
                {
                    m_advObject->SetFluxVector(
                        &UnsteadyAdvection::GetFluxVector, this);
                }
                break;
            }
            // Discontinuous field
            case MultiRegions::eDiscontinuous:
            {
                // Do not forwards transform initial condition
                m_homoInitialFwd = false;
                
                // Define the normal velocity fields
                if (m_fields[0]->GetTrace())
                {
                    m_traceVn = Array<OneD, NekDouble>(GetTraceNpoints());
                }

                string advName;
                string riemName;
                m_session->LoadSolverInfo(
                    "AdvectionType", advName, "WeakDG");
                m_advObject = SolverUtils::
                    GetAdvectionFactory().CreateInstance(advName, advName);
                if (m_specHP_dealiasing)
                {
                    m_advObject->SetFluxVector(
                        &UnsteadyAdvection::GetFluxVectorDeAlias, this);
                }
                else
                {
                    m_advObject->SetFluxVector(
                        &UnsteadyAdvection::GetFluxVector, this);
                }
                m_session->LoadSolverInfo(
                    "UpwindType", riemName, "Upwind");
                m_riemannSolver = SolverUtils::
                    GetRiemannSolverFactory().CreateInstance(
                        riemName, m_session);
                m_riemannSolver->SetScalar(
                    "Vn", &UnsteadyAdvection::GetNormalVelocity, this);

                m_advObject->SetRiemannSolver(m_riemannSolver);
                m_advObject->InitObject(m_session, m_fields);
                break;
            }
            default:
            {
                ASSERTL0(false, "Unsupported projection type.");
                break;
            }
        }

        // If explicit it computes RHS and PROJECTION for the time integration
        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs     (&UnsteadyAdvection::DoOdeRhs,        this);
            m_ode.DefineProjection (&UnsteadyAdvection::DoOdeProjection, this);
        }
        // Otherwise it gives an error (no implicit integration)
        else
        {
            ASSERTL0(false, "Implicit unsteady Advection not set up.");
        }

        remove( "maxU.txt" );
        m_outfile.open("maxU.txt", std::ios_base::app);
    }

    /**
     * @brief Unsteady linear advection equation destructor.
     */
    UnsteadyAdvection::~UnsteadyAdvection()
    {
    }

    /**
     * @brief Get the normal velocity for the linear advection equation.
     */
    Array<OneD, NekDouble> &UnsteadyAdvection::GetNormalVelocity()
    {
        // Number of trace (interface) points
        int i;
        int nTracePts = GetTraceNpoints();

        // Auxiliary variable to compute the normal velocity
        Array<OneD, NekDouble> tmp(nTracePts);

        // Reset the normal velocity
        Vmath::Zero(nTracePts, m_traceVn, 1);

        MultiRegions::ExpListSharedPtr traceExp = m_fields[0]->GetTrace();
        Array<OneD, NekDouble> xc(nTracePts), yc(nTracePts);
        traceExp->GetCoords(xc, yc);

        for (i = 0; i < m_velocity.num_elements(); ++i)
        {
            //m_fields[0]->ExtractTracePhys(m_velocity[i], tmp);

            //Hack: force velocity in correct direction for mortars
            //for (auto &z : tmp)
            //{
             //z = tmp[0];
            //}

            if (i == 0)
            {
                for (int j = 0; j < nTracePts; ++j)
                {
                    tmp[j] = -yc[j];
                }
            }
            else
            {
                for (int j = 0; j < nTracePts; ++j)
                {
                    tmp[j] = xc[j];
                }
            }

            Vmath::Vvtvp(nTracePts,
                         m_traceNormals[i], 1,
                         tmp,               1,
                         m_traceVn,         1,
                         m_traceVn,         1);
        }

        return m_traceVn;
    }

    /**
     * @brief Compute the right-hand side for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyAdvection::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
              Array<OneD,        Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        // Counter variable
        int i;

        // Number of fields (variables of the problem)
        int nVariables = inarray.num_elements();

        // Number of solution points
        int nSolutionPts = GetNpoints();

        // RHS computation using the new advection base class
        m_advObject->Advect(nVariables, m_fields, m_velocity, inarray,
                            outarray, time);

        // Negate the RHS
        for (i = 0; i < nVariables; ++i)
        {
            Vmath::Neg(nSolutionPts, outarray[i], 1);
        }
    }

    /**
     * @brief Compute the projection for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void UnsteadyAdvection::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD,       Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        // Counter variable
        int i;

        // Number of fields (variables of the problem)
        int nVariables = inarray.num_elements();

        // Set the boundary conditions
        SetBoundaryConditions(time);

        // Switch on the projection type (Discontinuous or Continuous)
        switch(m_projectionType)
        {
            // Discontinuous projection
            case MultiRegions::eDiscontinuous:
            {
                // Number of quadrature points
                int nQuadraturePts = GetNpoints();

                // Just copy over array
                for(i = 0; i < nVariables; ++i)
                {
                    Vmath::Vcopy(nQuadraturePts, inarray[i], 1, outarray[i], 1);
                }
                break;
            }

            // Continuous projection
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs(),0.0);
                for(i = 0; i < nVariables; ++i)
                {
                    m_fields[i]->FwdTrans(inarray[i], coeffs);
                    m_fields[i]->BwdTrans_IterPerExp(coeffs, outarray[i]);
                }
                break;
            }

            default:
                ASSERTL0(false,"Unknown projection scheme");
                break;
        }
    }

    /**
     * @brief Return the flux vector for the linear advection equation.
     *
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void UnsteadyAdvection::GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        ASSERTL1(flux[0].num_elements() == m_velocity.num_elements(),
                 "Dimension of flux array and velocity array do not match");

        int i , j;
        int nq = physfield[0].num_elements();

        for (i = 0; i < flux.num_elements(); ++i)
        {
            for (j = 0; j < flux[0].num_elements(); ++j)
            {
                Vmath::Vmul(nq, physfield[i], 1, m_velocity[j], 1,
                            flux[i][j], 1);
            }
        }
    }

    /**
     * @brief Return the flux vector for the linear advection equation using
     * the dealiasing technique.
     *
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void UnsteadyAdvection::GetFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        ASSERTL1(flux[0].num_elements() == m_velocity.num_elements(),
                 "Dimension of flux array and velocity array do not match");

        int i, j;
        int nq = physfield[0].num_elements();
        int nVariables = physfield.num_elements();

        // Factor to rescale 1d points in dealiasing
        NekDouble OneDptscale = 2;

        Array<OneD, Array<OneD, NekDouble> >
            advVel_plane(m_velocity.num_elements());

        // Get number of points to dealias a cubic non-linearity
        nq = m_fields[0]->Get1DScaledTotPoints(OneDptscale);

        // Initialisation of higher-space variables
        Array<OneD, Array<OneD, NekDouble> >physfieldInterp(nVariables);
        Array<OneD, Array<OneD, NekDouble> >velocityInterp(m_expdim);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >fluxInterp(nVariables);

        // Interpolation to higher space of physfield
        for (i = 0; i < nVariables; ++i)
        {
            physfieldInterp[i] = Array<OneD, NekDouble>(nq);
            fluxInterp[i] = Array<OneD, Array<OneD, NekDouble> >(m_expdim);
            for (j = 0; j < m_expdim; ++j)
            {
                fluxInterp[i][j] = Array<OneD, NekDouble>(nq);
            }

            m_fields[0]->PhysInterp1DScaled(
                OneDptscale, physfield[i], physfieldInterp[i]);
        }

        // Interpolation to higher space of velocity
        for (j = 0; j < m_expdim; ++j)
        {
            velocityInterp[j] = Array<OneD, NekDouble>(nq);

            m_fields[0]->PhysInterp1DScaled(
                OneDptscale, m_velocity[j], velocityInterp[j]);
        }

        // Evaluation of flux vector in the higher space
        for (i = 0; i < flux.num_elements(); ++i)
        {
            for (j = 0; j < flux[0].num_elements(); ++j)
            {
                Vmath::Vmul(nq, physfieldInterp[i], 1, velocityInterp[j], 1,
                            fluxInterp[i][j], 1);
            }
        }

        // Galerkin project solution back to original space
        for (i = 0; i < nVariables; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                m_fields[0]->PhysGalerkinProjection1DScaled(
                    OneDptscale, fluxInterp[i][j], flux[i][j]);
            }
        }
    }

    void UnsteadyAdvection::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        AdvectionSystem::v_GenerateSummary(s);
    }

    bool UnsteadyAdvection::v_PostIntegrate(int step)
    {
        //Find maximum phys value for graphing of decay
        int n = 25;     //Number of points to check in each quad (n x n)
        int itAll = 50; //Maximum iterations in newton raphson loop

        if (step % m_checksteps == 0 || step == 0)
        {
            //Searches for element containing maximum value at quad point
            NekDouble maxQ = std::numeric_limits<NekDouble>::min();
            int maxEl = -1;
            for (int el = 0; el < GetExpSize(); ++el)
            {
                Array<OneD, NekDouble> physVals(m_fields[0]->GetExp(el)->GetTotPoints(),
                m_fields[0]->GetPhys() + m_fields[0]->GetPhys_Offset(el));

                for (NekDouble maxQPoint : physVals)
                {
                    if (maxQPoint > maxQ)
                    {
                        maxQ = maxQPoint;
                        maxEl = el;
                    }
                }
            }

            //Finds the stencil of elements surrounding that containing the maximum quad point
            //Just quads (could generalise for quads + tris)
            //Does not work across interface (need to add check for if segment is on interface,
            // then include elements across from it also, got mappings for this ...
            // m_leftEdgeToMortar then m_mortarToRightEdge or other way.)
            std::unordered_set<int> stencilElementIDs;
            SpatialDomains::QuadGeomSharedPtr quadMax =
                    std::dynamic_pointer_cast<SpatialDomains::QuadGeom>(
                            m_fields[0]->GetExp(maxEl)->GetGeom());

            for (int i = 0; i < 4; ++i)
            {
                SpatialDomains::PointGeomSharedPtr vertex = quadMax->GetVertex(i);
                SpatialDomains::GeometryLinkSharedPtr edgesLink =
                        m_graph->GetElementsFromVertex(vertex);

                for (std::pair<SpatialDomains::GeometrySharedPtr, int> edge
                        : *edgesLink)
                {
                    SpatialDomains::Geometry1DSharedPtr edge1DGeom =
                            std::dynamic_pointer_cast<SpatialDomains::Geometry1D>(
                                    edge.first);
                    SpatialDomains::GeometryLinkSharedPtr elementsLink =
                            m_graph->GetElementsFromEdge(edge1DGeom);

                    for (std::pair<SpatialDomains::GeometrySharedPtr, int> element
                            : *elementsLink)
                    {
                        stencilElementIDs.insert(element.first->GetGlobalID());
                    }
                }
            }

            //**********************************
            // Output debug for stencil
            std::cout << "Quad max U found: " <<  maxQ << " in element " << maxEl <<
                      " giving stencil for Newton-Raphson of " << stencilElementIDs.size() <<
                      " elements: ";
            for (auto i : stencilElementIDs)
            {
                std::cout << i << ", ";
            }
            std::cout << "\b\b." << std::endl;
            //**********************************

            int maxIt = -1;
            NekDouble max = std::numeric_limits<NekDouble>::min(); //reset max value to try and remove spikes
            for (int el : stencilElementIDs)
            {
                const int nq = m_fields[0]->GetExp(el)->GetTotPoints();
                Array<OneD, NekDouble> x(nq, 0.0), y(nq, 0.0), xder(nq, 0.0),
                        yder(nq, 0.0), xder2(nq, 0.0), yder2(nq, 0.0), xdery(nq, 0.0),
                        yderx(nq, 0.0);

                x = m_fields[0]->GetPhys() + GetPhys_Offset(el);
                m_fields[0]->GetExp(el)->StdPhysDeriv(x, xder, yder);
                m_fields[0]->GetExp(el)->StdPhysDeriv(xder, xder2, xdery);
                m_fields[0]->GetExp(el)->StdPhysDeriv(yder, yderx, yder2);

                for (int i = 0; i < n; ++i)
                {
                    for (int j = 0; j < n; ++j)
                    {
                        Array<OneD, NekDouble> xi(2, 0.0),
                                xi_prev(2, std::numeric_limits<NekDouble>::min());

                        xi[0] = (i * (2.0 / (n - 1)) - 1.0);
                        xi[1] = (j * (2.0 / (n - 1)) - 1.0);

                        for (int k = 0; k < itAll; ++k)
                        {
                            CopyArray(xi, xi_prev);

                            NekDouble xc =  m_fields[0]->GetExp(el)->StdPhysEvaluate(xi, x);
                            NekDouble xc_derx =  m_fields[0]->GetExp(el)->StdPhysEvaluate(xi, xder);
                            NekDouble xc_dery =  m_fields[0]->GetExp(el)->StdPhysEvaluate(xi, yder);
                            NekDouble xc_derxx =  m_fields[0]->GetExp(el)->StdPhysEvaluate(xi, xder2);
                            NekDouble xc_deryy =  m_fields[0]->GetExp(el)->StdPhysEvaluate(xi, yder2);
                            NekDouble xc_derxy =  m_fields[0]->GetExp(el)->StdPhysEvaluate(xi, xdery);
                            NekDouble xc_deryx =  m_fields[0]->GetExp(el)->StdPhysEvaluate(xi, yderx);

                            //Newton's method for 2 variables
                            // i.e. xi = xi_prev - J^-1 * [xc_derx, xc_dery]
                            NekDouble det =
                                    1 / (xc_derxx * xc_deryy - xc_derxy * xc_deryx);
                            xi[0] = xi_prev[0] - ((det * xc_deryy) * xc_derx +
                                                  (det * -xc_derxy) * xc_dery);
                            xi[1] = xi_prev[1] - ((det * -xc_deryx) * xc_derx +
                                                  (det * xc_derxx) * xc_dery);

                            xi[0] = (xi[0] < -1) ? -1 : xi[0];
                            xi[1] = (xi[1] < -1) ? -1 : xi[1];
                            xi[0] = (xi[0] > 1) ? 1 : xi[0];
                            xi[1] = (xi[1] > 1) ? 1 : xi[1];

                            if ((abs(xi[0] - xi_prev[0]) < 1e-10)
                                && (abs(xi[1] - xi_prev[1]) < 1e-10))
                            {
                                if (xc > max)
                                {
                                    max = xc;
                                    maxEl = el;
                                }
                                break;
                            }

                            WARNINGL0((k + 1) != itAll,
                                      "Maximum iterations of " + std::to_string(itAll) +
                                      " reached while finding Max U.");
                        }
                    }
                }
            }

            std::cout << "Exact max U found: " <<  max << " in element " << maxEl <<
                      " at iteration " << maxIt + 1  << "." << std::endl;

            m_outfile << max << ", " << maxIt << ", " << maxQ << std::endl;
        }

        return false;
    }
}
