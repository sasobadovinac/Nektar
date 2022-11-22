///////////////////////////////////////////////////////////////////////////////
//
// File: IncNavierStokes.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Incompressible Navier Stokes class definition built on
// ADRBase class
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/algorithm/string.hpp>
#include <iomanip>

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Communication/Comm.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <SolverUtils/Filters/Filter.h>

#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <SolverUtils/Forcing/ForcingMovingReferenceFrame.h>
#include <tinyxml.h>

using namespace std;

namespace Nektar
{

std::string IncNavierStokes::eqTypeLookupIds[6] = {
    LibUtilities::SessionReader::RegisterEnumValue("EqType", "SteadyStokes",
                                                   eSteadyStokes),
    LibUtilities::SessionReader::RegisterEnumValue("EqType", "SteadyOseen",
                                                   eSteadyOseen),
    LibUtilities::SessionReader::RegisterEnumValue(
        "EqType", "SteadyNavierStokes", eSteadyNavierStokes),
    LibUtilities::SessionReader::RegisterEnumValue(
        "EqType", "SteadyLinearisedNS", eSteadyLinearisedNS),
    LibUtilities::SessionReader::RegisterEnumValue(
        "EqType", "UnsteadyNavierStokes", eUnsteadyNavierStokes),
    LibUtilities::SessionReader::RegisterEnumValue("EqType", "UnsteadyStokes",
                                                   eUnsteadyStokes)};

/**
 * Constructor. Creates ...
 *
 * \param
 * \param
 */
IncNavierStokes::IncNavierStokes(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), AdvectionSystem(pSession, pGraph),
      m_SmoothAdvection(false)
{
}

void IncNavierStokes::v_InitObject(bool DeclareField)
{
    AdvectionSystem::v_InitObject(DeclareField);

    int i, j;
    int numfields        = m_fields.size();
    std::string velids[] = {"u", "v", "w"};

    // Set up Velocity field to point to the first m_expdim of m_fields;
    m_velocity = Array<OneD, int>(m_spacedim);

    for (i = 0; i < m_spacedim; ++i)
    {
        for (j = 0; j < numfields; ++j)
        {
            std::string var = m_boundaryConditions->GetVariable(j);
            if (boost::iequals(velids[i], var))
            {
                m_velocity[i] = j;
                break;
            }

            ASSERTL0(j != numfields, "Failed to find field: " + var);
        }
    }

    // Set up equation type enum using kEquationTypeStr
    for (i = 0; i < (int)eEquationTypeSize; ++i)
    {
        bool match;
        m_session->MatchSolverInfo("EQTYPE", kEquationTypeStr[i], match, false);
        if (match)
        {
            m_equationType = (EquationType)i;
            break;
        }
    }
    ASSERTL0(i != eEquationTypeSize, "EQTYPE not found in SOLVERINFO section");

    m_session->LoadParameter("Kinvis", m_kinvis);

    // Default advection type per solver
    std::string vConvectiveType;
    switch (m_equationType)
    {
        case eUnsteadyStokes:
        case eSteadyLinearisedNS:
            vConvectiveType = "NoAdvection";
            break;
        case eUnsteadyNavierStokes:
        case eSteadyNavierStokes:
            vConvectiveType = "Convective";
            break;
        case eUnsteadyLinearisedNS:
            vConvectiveType = "Linearised";
            break;
        default:
            break;
    }

    // Check if advection type overridden
    if (m_session->DefinesTag("AdvectiveType") &&
        m_equationType != eUnsteadyStokes &&
        m_equationType != eSteadyLinearisedNS)
    {
        vConvectiveType = m_session->GetTag("AdvectiveType");
    }

    // Initialise advection
    m_advObject = SolverUtils::GetAdvectionFactory().CreateInstance(
        vConvectiveType, vConvectiveType);
    m_advObject->InitObject(m_session, m_fields);

    // Set up arrays for moving reference frame
    // Note: this must be done before the forcing
    m_movingFrameProjMat = bnu::identity_matrix<NekDouble>(3, 3);
    if (DefinedForcing("MovingReferenceFrame"))
    {
        m_movingFrameVelsxyz = Array<OneD, NekDouble>(6, 0.0);
        m_movingFrameTheta   = Array<OneD, NekDouble>(3, 0.0);
        m_pivotPoint         = Array<OneD, NekDouble>(3, 0.0);
        GetPivotPoint(m_pivotPoint);
    }

    // Forcing terms
    m_forcing = SolverUtils::Forcing::Load(m_session, shared_from_this(),
                                           m_fields, v_GetForceDimension());

    // check to see if any Robin boundary conditions and if so set
    // up m_field to boundary condition maps;
    m_fieldsBCToElmtID      = Array<OneD, Array<OneD, int>>(numfields);
    m_fieldsBCToTraceID     = Array<OneD, Array<OneD, int>>(numfields);
    m_fieldsRadiationFactor = Array<OneD, Array<OneD, NekDouble>>(numfields);

    for (i = 0; i < m_fields.size(); ++i)
    {
        bool Set = false;

        Array<OneD, const SpatialDomains::BoundaryConditionShPtr> BndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr> BndExp;
        int radpts = 0;

        BndConds = m_fields[i]->GetBndConditions();
        BndExp   = m_fields[i]->GetBndCondExpansions();
        for (int n = 0; n < BndConds.size(); ++n)
        {
            if (boost::iequals(BndConds[n]->GetUserDefined(), "Radiation"))
            {
                ASSERTL0(
                    BndConds[n]->GetBoundaryConditionType() ==
                        SpatialDomains::eRobin,
                    "Radiation boundary condition must be of type Robin <R>");

                if (Set == false)
                {
                    m_fields[i]->GetBoundaryToElmtMap(m_fieldsBCToElmtID[i],
                                                      m_fieldsBCToTraceID[i]);
                    Set = true;
                }
                radpts += BndExp[n]->GetTotPoints();
            }
            if (boost::iequals(BndConds[n]->GetUserDefined(),
                               "ZeroNormalComponent"))
            {
                ASSERTL0(BndConds[n]->GetBoundaryConditionType() ==
                             SpatialDomains::eDirichlet,
                         "Zero Normal Component boundary condition option must "
                         "be of type Dirichlet <D>");

                if (Set == false)
                {
                    m_fields[i]->GetBoundaryToElmtMap(m_fieldsBCToElmtID[i],
                                                      m_fieldsBCToTraceID[i]);
                    Set = true;
                }
            }
        }

        m_fieldsRadiationFactor[i] = Array<OneD, NekDouble>(radpts);

        radpts = 0; // reset to use as a counter

        for (int n = 0; n < BndConds.size(); ++n)
        {
            if (boost::iequals(BndConds[n]->GetUserDefined(), "Radiation"))
            {

                int npoints = BndExp[n]->GetNpoints();
                Array<OneD, NekDouble> x0(npoints, 0.0);
                Array<OneD, NekDouble> x1(npoints, 0.0);
                Array<OneD, NekDouble> x2(npoints, 0.0);
                Array<OneD, NekDouble> tmpArray;

                BndExp[n]->GetCoords(x0, x1, x2);

                LibUtilities::Equation coeff =
                    std::static_pointer_cast<
                        SpatialDomains::RobinBoundaryCondition>(BndConds[n])
                        ->m_robinPrimitiveCoeff;

                coeff.Evaluate(x0, x1, x2, m_time,
                               tmpArray = m_fieldsRadiationFactor[i] + radpts);
                // Vmath::Neg(npoints,tmpArray = m_fieldsRadiationFactor[i]+
                // radpts,1);
                radpts += npoints;
            }
        }
    }

    // Set up maping for womersley BC - and load variables
    for (int i = 0; i < m_fields.size(); ++i)
    {
        for (int n = 0; n < m_fields[i]->GetBndConditions().size(); ++n)
        {
            if (boost::istarts_with(
                    m_fields[i]->GetBndConditions()[n]->GetUserDefined(),
                    "Womersley"))
            {
                // assumes that boundary condition is applied in normal
                // direction and is decomposed for each direction. There could
                // be a unique file for each direction
                m_womersleyParams[i][n] =
                    MemoryManager<WomersleyParams>::AllocateSharedPtr(
                        m_spacedim);
                // Read in fourier coeffs and precompute coefficients
                SetUpWomersley(
                    i, n, m_fields[i]->GetBndConditions()[n]->GetUserDefined());

                m_fields[i]->GetBoundaryToElmtMap(m_fieldsBCToElmtID[i],
                                                  m_fieldsBCToTraceID[i]);
            }
        }
    }

    // Check if we have moving reference frame boundary conditions,if so, all
    // velocity vectors at those boundaries should be Drichlet and have a same
    // tag loop over the boundaries
    for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
    {
        // loop over the velocities
        for (int i = 0; i < m_velocity.size(); ++i)
        {
            if (boost::iequals(
                    m_fields[i]->GetBndConditions()[n]->GetUserDefined(),
                    "MovingFrameWall"))
            {
                // all velocities must have the same userdefined tag and all
                // must be Dirichlet
                for (int j = 0; j < m_velocity.size(); ++j)
                {
                    ASSERTL0(m_fields[j]
                                     ->GetBndConditions()[n]
                                     ->GetBoundaryConditionType() ==
                                 SpatialDomains::eDirichlet,
                             "All velocities with userdefined tag: "
                             "\"MovingFrameWall\" must be Dirichlet boundary "
                             "condition");
                    ASSERTL0(boost::iequals(m_fields[j]
                                                ->GetBndConditions()[n]
                                                ->GetUserDefined(),
                                            "MovingFrameWall"),
                             "If any of the velocity components at a "
                             "boundary is defined as \"MovingFrameWall\", "
                             "the rest of velocity components are also "
                             "must be defined as \"MovingFrameWall\" ");
                }
            }

            if (boost::iequals(
                    m_fields[i]->GetBndConditions()[n]->GetUserDefined(),
                    "MovingFrameDomainVel"))
            {
                // all velocities must have the same userdefined tag and all
                // must be Dirichlet
                for (int j = 0; j < m_velocity.size(); ++j)
                {
                    ASSERTL0(
                        m_fields[j]
                                ->GetBndConditions()[n]
                                ->GetBoundaryConditionType() ==
                            SpatialDomains::eDirichlet,
                        "All velocities with userdefined tag: "
                        "\"MovingFrameDomainVel\" must be Dirichlet boundary "
                        "condition");
                    ASSERTL0(boost::iequals(m_fields[j]
                                                ->GetBndConditions()[n]
                                                ->GetUserDefined(),
                                            "MovingFrameDomainVel"),
                             "If any of the velocity components at a "
                             "boundary is defined as \"MovingFrameDomainVel\", "
                             "the rest of velocity components are also "
                             "must be defined as \"MovingFrameDomainVel\" ");
                }
            }
        }
    }

    // Set up Field Meta Data for output files
    m_fieldMetaDataMap["Kinvis"] = boost::lexical_cast<std::string>(m_kinvis);
    m_fieldMetaDataMap["TimeStep"] =
        boost::lexical_cast<std::string>(m_timestep);

    // Check to see if the metadata for moving reference frame is defined
    if (DefinedForcing("MovingReferenceFrame"))
    {
        if (m_fieldMetaDataMap != LibUtilities::NullFieldMetaDataMap)
        {
            std::vector<std::string> vSuffix = {"_x", "_y", "_z"};
            for (int i = 0; i < 3; ++i)
            {
                NekDouble dTheta{0};
                std::string sTheta = "Theta" + vSuffix[i];
                auto it            = m_fieldMetaDataMap.find(sTheta);
                if (it != m_fieldMetaDataMap.end())
                {
                    dTheta = boost::lexical_cast<NekDouble>(it->second);
                    m_movingFrameTheta[i] = dTheta;
                }
                else
                {
                    m_fieldMetaDataMap[sTheta] =
                        boost::lexical_cast<std::string>(m_movingFrameTheta[i]);
                }
            }
        }
        else
        {
            std::vector<std::string> vSuffix = {"_x", "_y", "_z"};
            for (int i = 0; i < 3; ++i)
            {
                std::string sTheta = "Theta" + vSuffix[i];
                m_fieldMetaDataMap[sTheta] =
                    boost::lexical_cast<std::string>(m_movingFrameTheta[i]);
            }
        }
    }
}

/**
 * Destructor
 */
IncNavierStokes::~IncNavierStokes(void)
{
}

/**
 * Evaluation -N(V) for all fields except pressure using m_velocity
 */
void IncNavierStokes::EvaluateAdvectionTerms(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int i;
    int VelDim = m_velocity.size();
    Array<OneD, Array<OneD, NekDouble>> velocity(VelDim);

    int npoints = m_fields[0]->GetNpoints();
    for (i = 0; i < VelDim; ++i)
    {
        velocity[i] = Array<OneD, NekDouble>(npoints);
        Vmath::Vcopy(npoints, inarray[m_velocity[i]], 1, velocity[i], 1);
    }
    for (auto &x : m_forcing)
    {
        x->PreApply(m_fields, velocity, velocity, time);
    }

    m_advObject->Advect(m_nConvectiveFields, m_fields, velocity, inarray,
                        outarray, time);
}

/**
 * Time dependent boundary conditions updating
 */
void IncNavierStokes::SetBoundaryConditions(NekDouble time)
{
    int i, n;
    std::string varName;
    int nvariables = m_fields.size();

    for (i = 0; i < nvariables; ++i)
    {
        for (n = 0; n < m_fields[i]->GetBndConditions().size(); ++n)
        {
            if (m_fields[i]->GetBndConditions()[n]->IsTimeDependent())
            {
                varName = m_session->GetVariable(i);
                m_fields[i]->EvaluateBoundaryConditions(time, varName);
            }
            else if (boost::istarts_with(
                         m_fields[i]->GetBndConditions()[n]->GetUserDefined(),
                         "Womersley"))
            {
                SetWomersleyBoundary(i, n);
            }
        }

        // Set Radiation conditions if required
        SetRadiationBoundaryForcing(i);
    }

    // Enforcing the boundary conditions (Inlet and wall) for the
    // Moving reference frame
    SetMovingReferenceFrameBCs(time);
    SetZeroNormalVelocity();
}

/**
 * Probably should be pushed back into ContField?
 */
void IncNavierStokes::SetRadiationBoundaryForcing(int fieldid)
{
    int i, n;

    Array<OneD, const SpatialDomains::BoundaryConditionShPtr> BndConds;
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp;

    BndConds = m_fields[fieldid]->GetBndConditions();
    BndExp   = m_fields[fieldid]->GetBndCondExpansions();

    LocalRegions::ExpansionSharedPtr elmt;
    StdRegions::StdExpansionSharedPtr Bc;

    int cnt;
    int elmtid, nq, offset, boundary;
    Array<OneD, NekDouble> Bvals, U;
    int cnt1 = 0;

    for (cnt = n = 0; n < BndConds.size(); ++n)
    {
        std::string type = BndConds[n]->GetUserDefined();

        if ((BndConds[n]->GetBoundaryConditionType() ==
             SpatialDomains::eRobin) &&
            (boost::iequals(type, "Radiation")))
        {
            for (i = 0; i < BndExp[n]->GetExpSize(); ++i, cnt++)
            {
                elmtid = m_fieldsBCToElmtID[m_velocity[fieldid]][cnt];
                elmt   = m_fields[fieldid]->GetExp(elmtid);
                offset = m_fields[fieldid]->GetPhys_Offset(elmtid);

                U  = m_fields[fieldid]->UpdatePhys() + offset;
                Bc = BndExp[n]->GetExp(i);

                boundary = m_fieldsBCToTraceID[fieldid][cnt];

                // Get edge values and put into ubc
                nq = Bc->GetTotPoints();
                Array<OneD, NekDouble> ubc(nq);
                elmt->GetTracePhysVals(boundary, Bc, U, ubc);

                Vmath::Vmul(nq,
                            &m_fieldsRadiationFactor
                                [fieldid][cnt1 + BndExp[n]->GetPhys_Offset(i)],
                            1, &ubc[0], 1, &ubc[0], 1);

                Bvals =
                    BndExp[n]->UpdateCoeffs() + BndExp[n]->GetCoeff_Offset(i);

                Bc->IProductWRTBase(ubc, Bvals);
            }
            cnt1 += BndExp[n]->GetTotPoints();
        }
        else
        {
            cnt += BndExp[n]->GetExpSize();
        }
    }
}

void IncNavierStokes::SetZeroNormalVelocity()
{
    // use static trip since cannot use UserDefinedTag for zero
    // velocity and have time dependent conditions
    static bool Setup = false;

    if (Setup == true)
    {
        return;
    }
    Setup = true;

    int i, n;

    Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr>>
        BndConds(m_spacedim);
    Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr>> BndExp(m_spacedim);

    for (i = 0; i < m_spacedim; ++i)
    {
        BndConds[i] = m_fields[m_velocity[i]]->GetBndConditions();
        BndExp[i]   = m_fields[m_velocity[i]]->GetBndCondExpansions();
    }

    LocalRegions::ExpansionSharedPtr elmt, Bc;

    int cnt;
    int elmtid, nq, boundary;

    Array<OneD, Array<OneD, NekDouble>> normals;
    Array<OneD, NekDouble> Bphys, Bcoeffs;

    int fldid = m_velocity[0];

    for (cnt = n = 0; n < BndConds[0].size(); ++n)
    {
        if ((BndConds[0][n]->GetBoundaryConditionType() ==
             SpatialDomains::eDirichlet) &&
            (boost::iequals(BndConds[0][n]->GetUserDefined(),
                            "ZeroNormalComponent")))
        {
            for (i = 0; i < BndExp[0][n]->GetExpSize(); ++i, cnt++)
            {
                elmtid   = m_fieldsBCToElmtID[fldid][cnt];
                elmt     = m_fields[0]->GetExp(elmtid);
                boundary = m_fieldsBCToTraceID[fldid][cnt];

                normals = elmt->GetTraceNormal(boundary);

                nq = BndExp[0][n]->GetExp(i)->GetTotPoints();
                Array<OneD, NekDouble> normvel(nq, 0.0);

                for (int k = 0; k < m_spacedim; ++k)
                {
                    Bphys = BndExp[k][n]->UpdatePhys() +
                            BndExp[k][n]->GetPhys_Offset(i);
                    Bc = BndExp[k][n]->GetExp(i);
                    Vmath::Vvtvp(nq, normals[k], 1, Bphys, 1, normvel, 1,
                                 normvel, 1);
                }

                // negate normvel for next step
                Vmath::Neg(nq, normvel, 1);

                for (int k = 0; k < m_spacedim; ++k)
                {
                    Bphys = BndExp[k][n]->UpdatePhys() +
                            BndExp[k][n]->GetPhys_Offset(i);
                    Bcoeffs = BndExp[k][n]->UpdateCoeffs() +
                              BndExp[k][n]->GetCoeff_Offset(i);
                    Bc = BndExp[k][n]->GetExp(i);
                    Vmath::Vvtvp(nq, normvel, 1, normals[k], 1, Bphys, 1, Bphys,
                                 1);
                    Bc->FwdTransBndConstrained(Bphys, Bcoeffs);
                }
            }
        }
        else
        {
            cnt += BndExp[0][n]->GetExpSize();
        }
    }
}

/**
 *  Womersley boundary condition defintion
 */
void IncNavierStokes::SetWomersleyBoundary(const int fldid, const int bndid)
{
    ASSERTL1(m_womersleyParams.count(bndid) == 1,
             "Womersley parameters for this boundary have not been set up");

    WomersleyParamsSharedPtr WomParam = m_womersleyParams[fldid][bndid];
    NekComplexDouble zvel;
    int i, j, k;

    int M_coeffs = WomParam->m_wom_vel.size();

    NekDouble T           = WomParam->m_period;
    NekDouble axis_normal = WomParam->m_axisnormal[fldid];

    // Womersley Number
    NekComplexDouble omega_c(2.0 * M_PI / T, 0.0);
    NekComplexDouble k_c(0.0, 0.0);
    NekComplexDouble m_time_c(m_time, 0.0);
    NekComplexDouble zi(0.0, 1.0);
    NekComplexDouble i_pow_3q2(-1.0 / sqrt(2.0), 1.0 / sqrt(2.0));

    MultiRegions::ExpListSharedPtr BndCondExp;
    BndCondExp = m_fields[fldid]->GetBndCondExpansions()[bndid];

    StdRegions::StdExpansionSharedPtr bc;
    int cnt = 0;
    int nfq;
    Array<OneD, NekDouble> Bvals;
    int exp_npts = BndCondExp->GetExpSize();
    Array<OneD, NekDouble> wbc(exp_npts, 0.0);

    Array<OneD, NekComplexDouble> zt(M_coeffs);

    // preallocate the exponent
    for (k = 1; k < M_coeffs; k++)
    {
        k_c   = NekComplexDouble((NekDouble)k, 0.0);
        zt[k] = std::exp(zi * omega_c * k_c * m_time_c);
    }

    // Loop over each element in an expansion
    for (i = 0; i < exp_npts; ++i, cnt++)
    {
        // Get Boundary and trace expansion
        bc  = BndCondExp->GetExp(i);
        nfq = bc->GetTotPoints();
        Array<OneD, NekDouble> wbc(nfq, 0.0);

        // Compute womersley solution
        for (j = 0; j < nfq; j++)
        {
            wbc[j] = WomParam->m_poiseuille[i][j];
            for (k = 1; k < M_coeffs; k++)
            {
                zvel   = WomParam->m_zvel[i][j][k] * zt[k];
                wbc[j] = wbc[j] + zvel.real();
            }
        }

        // Multiply w by normal to get u,v,w component of velocity
        Vmath::Smul(nfq, axis_normal, wbc, 1, wbc, 1);
        // get the offset
        Bvals = BndCondExp->UpdateCoeffs() + BndCondExp->GetCoeff_Offset(i);

        // Push back to Coeff space
        bc->FwdTrans(wbc, Bvals);
    }
}

void IncNavierStokes::SetUpWomersley(const int fldid, const int bndid,
                                     std::string womStr)
{
    std::string::size_type indxBeg = womStr.find_first_of(':') + 1;
    string filename                = womStr.substr(indxBeg, string::npos);

    TiXmlDocument doc(filename);

    bool loadOkay = doc.LoadFile();
    ASSERTL0(loadOkay,
             (std::string("Failed to load file: ") + filename).c_str());

    TiXmlHandle docHandle(&doc);

    int err; /// Error value returned by TinyXML.

    TiXmlElement *nektar = doc.FirstChildElement("NEKTAR");
    ASSERTL0(nektar, "Unable to find NEKTAR tag in file.");

    TiXmlElement *wombc = nektar->FirstChildElement("WOMERSLEYBC");
    ASSERTL0(wombc, "Unable to find WOMERSLEYBC tag in file.");

    // read womersley parameters
    TiXmlElement *womparam = wombc->FirstChildElement("WOMPARAMS");
    ASSERTL0(womparam, "Unable to find WOMPARAMS tag in file.");

    // Input coefficients
    TiXmlElement *params = womparam->FirstChildElement("W");
    map<std::string, std::string> Wparams;

    // read parameter list
    while (params)
    {

        std::string propstr;
        propstr = params->Attribute("PROPERTY");

        ASSERTL0(!propstr.empty(),
                 "Failed to read PROPERTY value Womersley BC Parameter");

        std::string valstr;
        valstr = params->Attribute("VALUE");

        ASSERTL0(!valstr.empty(),
                 "Failed to read VALUE value Womersley BC Parameter");

        std::transform(propstr.begin(), propstr.end(), propstr.begin(),
                       ::toupper);
        Wparams[propstr] = valstr;

        params = params->NextSiblingElement("W");
    }
    bool parseGood;

    // Read parameters

    ASSERTL0(
        Wparams.count("RADIUS") == 1,
        "Failed to find Radius parameter in Womersley boundary conditions");
    std::vector<NekDouble> rad;
    ParseUtils::GenerateVector(Wparams["RADIUS"], rad);
    m_womersleyParams[fldid][bndid]->m_radius = rad[0];

    ASSERTL0(
        Wparams.count("PERIOD") == 1,
        "Failed to find period parameter in Womersley boundary conditions");
    std::vector<NekDouble> period;
    parseGood = ParseUtils::GenerateVector(Wparams["PERIOD"], period);
    m_womersleyParams[fldid][bndid]->m_period = period[0];

    ASSERTL0(
        Wparams.count("AXISNORMAL") == 1,
        "Failed to find axisnormal parameter in Womersley boundary conditions");
    std::vector<NekDouble> anorm;
    parseGood = ParseUtils::GenerateVector(Wparams["AXISNORMAL"], anorm);
    m_womersleyParams[fldid][bndid]->m_axisnormal[0] = anorm[0];
    m_womersleyParams[fldid][bndid]->m_axisnormal[1] = anorm[1];
    m_womersleyParams[fldid][bndid]->m_axisnormal[2] = anorm[2];

    ASSERTL0(
        Wparams.count("AXISPOINT") == 1,
        "Failed to find axispoint parameter in Womersley boundary conditions");
    std::vector<NekDouble> apt;
    parseGood = ParseUtils::GenerateVector(Wparams["AXISPOINT"], apt);
    m_womersleyParams[fldid][bndid]->m_axispoint[0] = apt[0];
    m_womersleyParams[fldid][bndid]->m_axispoint[1] = apt[1];
    m_womersleyParams[fldid][bndid]->m_axispoint[2] = apt[2];

    // Read Temporal Fourier Coefficients.

    // Find the FourierCoeff tag
    TiXmlElement *coeff = wombc->FirstChildElement("FOURIERCOEFFS");

    // Input coefficients
    TiXmlElement *fval = coeff->FirstChildElement("F");

    int indx;
    int nextFourierCoeff = -1;

    while (fval)
    {
        nextFourierCoeff++;

        TiXmlAttribute *fvalAttr = fval->FirstAttribute();
        std::string attrName(fvalAttr->Name());

        ASSERTL0(attrName == "ID",
                 (std::string("Unknown attribute name: ") + attrName).c_str());

        err = fvalAttr->QueryIntValue(&indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");

        std::string coeffStr = fval->FirstChild()->ToText()->ValueStr();
        vector<NekDouble> coeffvals;

        parseGood = ParseUtils::GenerateVector(coeffStr, coeffvals);
        ASSERTL0(
            parseGood,
            (std::string("Problem reading value of fourier coefficient, ID=") +
             boost::lexical_cast<string>(indx))
                .c_str());
        ASSERTL1(
            coeffvals.size() == 2,
            (std::string(
                 "Have not read two entries of Fourier coefficicent from ID=" +
                 boost::lexical_cast<string>(indx))
                 .c_str()));

        m_womersleyParams[fldid][bndid]->m_wom_vel.push_back(
            NekComplexDouble(coeffvals[0], coeffvals[1]));

        fval = fval->NextSiblingElement("F");
    }

    // starting point of precalculation
    int i, j, k;
    // M fourier coefficients
    int M_coeffs = m_womersleyParams[fldid][bndid]->m_wom_vel.size();
    NekDouble R  = m_womersleyParams[fldid][bndid]->m_radius;
    NekDouble T  = m_womersleyParams[fldid][bndid]->m_period;
    Array<OneD, NekDouble> x0 = m_womersleyParams[fldid][bndid]->m_axispoint;

    NekComplexDouble rqR;
    // Womersley Number
    NekComplexDouble omega_c(2.0 * M_PI / T, 0.0);
    NekComplexDouble alpha_c(R * sqrt(omega_c.real() / m_kinvis), 0.0);
    NekComplexDouble z1(1.0, 0.0);
    NekComplexDouble i_pow_3q2(-1.0 / sqrt(2.0), 1.0 / sqrt(2.0));

    MultiRegions::ExpListSharedPtr BndCondExp;
    BndCondExp = m_fields[fldid]->GetBndCondExpansions()[bndid];

    StdRegions::StdExpansionSharedPtr bc;
    int cnt = 0;
    int nfq;
    Array<OneD, NekDouble> Bvals;

    int exp_npts = BndCondExp->GetExpSize();
    Array<OneD, NekDouble> wbc(exp_npts, 0.0);

    // allocate time indepedent variables
    m_womersleyParams[fldid][bndid]->m_poiseuille =
        Array<OneD, Array<OneD, NekDouble>>(exp_npts);
    m_womersleyParams[fldid][bndid]->m_zvel =
        Array<OneD, Array<OneD, Array<OneD, NekComplexDouble>>>(exp_npts);
    // could use M_coeffs - 1 but need to avoid complicating things
    Array<OneD, NekComplexDouble> zJ0(M_coeffs);
    Array<OneD, NekComplexDouble> lamda_n(M_coeffs);
    Array<OneD, NekComplexDouble> k_c(M_coeffs);
    NekComplexDouble zJ0r;

    for (k = 1; k < M_coeffs; k++)
    {
        k_c[k]     = NekComplexDouble((NekDouble)k, 0.0);
        lamda_n[k] = i_pow_3q2 * alpha_c * sqrt(k_c[k]);
        zJ0[k]     = Polylib::ImagBesselComp(0, lamda_n[k]);
    }

    // Loop over each element in an expansion
    for (i = 0; i < exp_npts; ++i, cnt++)
    {
        // Get Boundary and trace expansion
        bc  = BndCondExp->GetExp(i);
        nfq = bc->GetTotPoints();

        Array<OneD, NekDouble> x(nfq, 0.0);
        Array<OneD, NekDouble> y(nfq, 0.0);
        Array<OneD, NekDouble> z(nfq, 0.0);
        bc->GetCoords(x, y, z);

        m_womersleyParams[fldid][bndid]->m_poiseuille[i] =
            Array<OneD, NekDouble>(nfq);
        m_womersleyParams[fldid][bndid]->m_zvel[i] =
            Array<OneD, Array<OneD, NekComplexDouble>>(nfq);

        // Compute coefficients
        for (j = 0; j < nfq; j++)
        {
            rqR = NekComplexDouble(sqrt((x[j] - x0[0]) * (x[j] - x0[0]) +
                                        (y[j] - x0[1]) * (y[j] - x0[1]) +
                                        (z[j] - x0[2]) * (z[j] - x0[2])) /
                                       R,
                                   0.0);

            // Compute Poiseulle Flow
            m_womersleyParams[fldid][bndid]->m_poiseuille[i][j] =
                m_womersleyParams[fldid][bndid]->m_wom_vel[0].real() *
                (1. - rqR.real() * rqR.real());

            m_womersleyParams[fldid][bndid]->m_zvel[i][j] =
                Array<OneD, NekComplexDouble>(M_coeffs);

            // compute the velocity information
            for (k = 1; k < M_coeffs; k++)
            {
                zJ0r = Polylib::ImagBesselComp(0, rqR * lamda_n[k]);
                m_womersleyParams[fldid][bndid]->m_zvel[i][j][k] =
                    m_womersleyParams[fldid][bndid]->m_wom_vel[k] *
                    (z1 - (zJ0r / zJ0[k]));
            }
        }
    }
}

/**
 * Set boundary conditions for moving frame of reference
 */
void IncNavierStokes::SetMovingReferenceFrameBCs(const NekDouble &time)
{
    SetMRFWallBCs(time);
    SetMRFDomainVelBCs(time);
}

/**
 * Set Wall boundary conditions for moving frame of reference
 */
void IncNavierStokes::SetMRFWallBCs(const NekDouble &time)
{
    boost::ignore_unused(time);

    // for the wall we need to calculate:
    // [V_wall]_xyz = [V_frame]_xyz + [Omega X r]_xyz
    // Note all vectors must be in moving frame coordinates xyz
    // not in inertial frame XYZ

    Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr>>
        BndConds(m_spacedim);
    Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr>> BndExp(m_spacedim);

    for (int i = 0; i < m_spacedim; ++i)
    {
        BndConds[i] = m_fields[m_velocity[i]]->GetBndConditions();
        BndExp[i]   = m_fields[m_velocity[i]]->GetBndCondExpansions();
    }

    int npoints;
    Array<OneD, NekDouble> Bphys, Bcoeffs;

    // loop over the boundary regions
    for (int n = 0; n < BndExp[0].size(); ++n)
    {
        if (BndConds[0][n]->GetBoundaryConditionType() ==
                SpatialDomains::eDirichlet &&
            (boost::iequals(BndConds[0][n]->GetUserDefined(),
                            "MovingFrameWall")))
        {
            npoints = BndExp[0][n]->GetNpoints();
            //
            // define arrays to calculate the prescribed velocities and modifed
            // ones
            Array<OneD, Array<OneD, NekDouble>> velocities(m_velocity.size());
            Array<OneD, Array<OneD, NekDouble>> unitArray(m_velocity.size());
            Array<OneD, Array<OneD, NekDouble>> coords(3);
            for (int k = 0; k < m_velocity.size(); ++k)
            {
                velocities[k] = Array<OneD, NekDouble>(npoints, 0.0);
                unitArray[k]  = Array<OneD, NekDouble>(npoints, 1.0);
            }
            Array<OneD, NekDouble> tmp(npoints, 0.0);
            for (int k = 0; k < 3; ++k)
            {
                coords[k] = Array<OneD, NekDouble>(npoints, 0.0);
            }
            BndExp[0][n]->GetCoords(coords[0], coords[1], coords[2]);

            // move the centre to the location of pivot
            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Sadd(npoints, -m_pivotPoint[i], coords[i], 1, coords[i],
                            1);
            }
            /////////////////////////////////////
            // Note that both Omega and the absolute velocities have been
            // expressed in the moving reference frame unit vectors
            // therefore the result is in moving ref frame and no furhter
            // transformaton is required
            //
            // compute Omega X r = vx ex + vy ey + vz ez
            // Note OmegaX : movingFrameVelsxyz[3]
            // Note OmegaY : movingFrameVelsxyz[4]
            // Note OmegaZ : movingFrameVelsxyz[5]
            //
            // vx = OmegaY*z-OmegaZ*y
            Vmath::Smul(npoints, -1 * m_movingFrameVelsxyz[5], coords[1], 1,
                        velocities[0], 1);
            // vy = OmegaZ*x-OmegaX*z
            Vmath::Smul(npoints, m_movingFrameVelsxyz[5], coords[0], 1,
                        velocities[1], 1);
            if (m_spacedim == 3)
            {
                // add the OmegaY*z to vx
                Vmath::Svtvp(npoints, m_movingFrameVelsxyz[4], coords[2], 1,
                             velocities[0], 1, velocities[0], 1);
                // add the -OmegaX*z to vy
                Vmath::Svtvp(npoints, -1 * m_movingFrameVelsxyz[3], coords[2],
                             1, velocities[1], 1, velocities[1], 1);

                // vz = OmegaX*y-OmegaY*x
                Vmath::Svtsvtp(npoints, m_movingFrameVelsxyz[3], coords[1], 1,
                               -1.0 * m_movingFrameVelsxyz[4], coords[0], 1,
                               velocities[2], 1);
            }

            // add the translation velocity
            for (int k = 0; k < m_spacedim; ++k)
            {
                Vmath::Sadd(npoints, m_movingFrameVelsxyz[k], velocities[k], 1,
                            velocities[k], 1);
            }

            // update the boundary values
            for (int k = 0; k < m_spacedim; ++k)
            {
                Vmath::Vmul(npoints, unitArray[k], 1, velocities[k], 1,
                            BndExp[k][n]->UpdatePhys(), 1);
            }
            // update the coefficients
            for (int k = 0; k < m_spacedim; ++k)
            {
                if (m_fields[k]->GetExpType() == MultiRegions::e3DH1D)
                {
                    BndExp[k][n]->SetWaveSpace(false);
                }
                BndExp[k][n]->FwdTransBndConstrained(
                    BndExp[k][n]->GetPhys(), BndExp[k][n]->UpdateCoeffs());
            }
        }
    }
}

/**
 * Set inlet boundary conditions for moving frame of reference
 */
void IncNavierStokes::SetMRFDomainVelBCs(const NekDouble &time)
{
    // The inlet conditions for the velocities given in the session xml file
    // however, those are in the inertial reference coordinate XYZ and
    // needed to be converted to the moving reference coordinate xyz

    // define arrays to calculate the prescribed velocities and modifed ones
    Array<OneD, Array<OneD, NekDouble>> definedVels(m_velocity.size());
    Array<OneD, Array<OneD, NekDouble>> velocities(m_velocity.size());
    Array<OneD, Array<OneD, NekDouble>> unitArray(m_velocity.size());
    Array<OneD, Array<OneD, NekDouble>> coords(3);

    Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr>>
        BndConds(m_spacedim);
    Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr>> BndExp(m_spacedim);

    for (int i = 0; i < m_spacedim; ++i)
    {
        BndConds[i] = m_fields[m_velocity[i]]->GetBndConditions();
        BndExp[i]   = m_fields[m_velocity[i]]->GetBndCondExpansions();
    }

    int npoints;
    Array<OneD, NekDouble> Bphys, Bcoeffs;

    // loop over the boundary regions
    for (int n = 0; n < BndExp[0].size(); ++n)
    {
        if (BndConds[0][n]->GetBoundaryConditionType() ==
                SpatialDomains::eDirichlet &&
            (boost::iequals(BndConds[0][n]->GetUserDefined(),
                            "MovingFrameDomainVel")))
        {
            npoints = BndExp[0][n]->GetNpoints();
            for (int k = 0; k < m_velocity.size(); ++k)
            {
                definedVels[k] = Array<OneD, NekDouble>(npoints, 0.0);
                velocities[k]  = Array<OneD, NekDouble>(npoints, 0.0);
                unitArray[k]   = Array<OneD, NekDouble>(npoints, 1.0);
            }
            Array<OneD, NekDouble> tmp(npoints, 0.0);
            for (int k = 0; k < 3; ++k)
            {
                coords[k] = Array<OneD, NekDouble>(npoints, 0.0);
            }
            BndExp[0][n]->GetCoords(coords[0], coords[1], coords[2]);

            // loop over the velocity fields and compute the boundary
            // condition
            for (int k = 0; k < m_velocity.size(); ++k)
            {
                LibUtilities::Equation condition =
                    std::static_pointer_cast<
                        SpatialDomains::DirichletBoundaryCondition>(
                        BndConds[k][n])
                        ->m_dirichletCondition;
                // Evaluate
                condition.Evaluate(coords[0], coords[1], coords[2], time,
                                   definedVels[k]);
            }

            // We have all velocity components
            // transform them to the moving refernce frame
            for (int l = 0; l < m_spacedim; ++l)
            {
                Array<OneD, NekDouble> tmp0(npoints, 0.0);
                Array<OneD, NekDouble> tmp1(npoints, 0.0);
                Array<OneD, NekDouble> tmp2(npoints, 0.0);
                for (int m = 0; m < m_spacedim; ++m)
                {
                    Vmath::Svtvp(npoints, m_movingFrameProjMat(l, m),
                                 tmp0 = definedVels[m], 1, tmp1 = velocities[l],
                                 1, tmp2 = velocities[l], 1);
                }
            }

            // update the boundary values
            for (int k = 0; k < m_spacedim; ++k)
            {
                Vmath::Vmul(npoints, unitArray[k], 1, velocities[k], 1,
                            BndExp[k][n]->UpdatePhys(), 1);
            }
            // update the coefficients
            for (int k = 0; k < m_spacedim; ++k)
            {
                if (m_fields[k]->GetExpType() == MultiRegions::e3DH1D)
                {
                    BndExp[k][n]->SetWaveSpace(false);
                }
                BndExp[k][n]->FwdTransBndConstrained(
                    BndExp[k][n]->GetPhys(), BndExp[k][n]->UpdateCoeffs());
            }
        }
    }
}

/**
 * Add an additional forcing term programmatically.
 */
void IncNavierStokes::AddForcing(const SolverUtils::ForcingSharedPtr &pForce)
{
    m_forcing.push_back(pForce);
}

/**
 *
 */
Array<OneD, NekDouble> IncNavierStokes::v_GetMaxStdVelocity(
    const NekDouble SpeedSoundFactor)
{
    boost::ignore_unused(SpeedSoundFactor);
    int nvel  = m_velocity.size();
    int nelmt = m_fields[0]->GetExpSize();

    Array<OneD, NekDouble> stdVelocity(nelmt, 0.0);
    Array<OneD, Array<OneD, NekDouble>> velfields;

    if (m_HomogeneousType == eHomogeneous1D) // just do check on 2D info
    {
        velfields = Array<OneD, Array<OneD, NekDouble>>(2);

        for (int i = 0; i < 2; ++i)
        {
            velfields[i] = m_fields[m_velocity[i]]->UpdatePhys();
        }
    }
    else
    {
        velfields = Array<OneD, Array<OneD, NekDouble>>(nvel);

        for (int i = 0; i < nvel; ++i)
        {
            velfields[i] = m_fields[m_velocity[i]]->UpdatePhys();
        }
    }

    stdVelocity = m_extrapolation->GetMaxStdVelocity(velfields);

    return stdVelocity;
}

/**
 *
 */
void IncNavierStokes::v_GetPressure(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &pressure)
{
    pressure = physfield[m_nConvectiveFields];
}

/**
 *
 */
void IncNavierStokes::v_GetDensity(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &density)
{
    int nPts = physfield[0].size();
    Vmath::Fill(nPts, 1.0, density, 1);
}

/**
 *
 */
void IncNavierStokes::v_GetVelocity(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &velocity)
{
    for (int i = 0; i < m_spacedim; ++i)
    {
        velocity[i] = physfield[i];
    }
}

/**
 * Function to set the moving frame velocities calucated in the forcing
 * this gives access to the moving reference forcing to set the velocities
 * to be later used in enforcing the boundary condition in IncNavierStokes
 * class
 */
void IncNavierStokes::v_SetMovingFrameVelocities(
    const Array<OneD, NekDouble> &vFrameVels)
{
    ASSERTL0(vFrameVels.size() == m_movingFrameVelsxyz.size(),
             "Arrays have different dimensions, cannot set moving frame "
             "velocities");
    Vmath::Vcopy(vFrameVels.size(), vFrameVels, 1, m_movingFrameVelsxyz, 1);
}

void IncNavierStokes::v_GetMovingFrameVelocities(
    Array<OneD, NekDouble> &vFrameVels)
{
    ASSERTL0(vFrameVels.size() == m_movingFrameVelsxyz.size(),
             "Arrays have different dimensions, cannot get moving frame "
             "velocities");
    unsigned int size = m_movingFrameVelsxyz.size();
    Vmath::Vcopy(size, m_movingFrameVelsxyz, 1, vFrameVels, 1);
}

/**
 * Function to set the rotation matrix computed in the forcing
 * this gives access to the moving reference forcing to set the Projection
 * matrix to be used later in IncNavierStokes calss for enforcing the
 * boundary conditions
 */
void IncNavierStokes::v_SetMovingFrameProjectionMat(
    const bnu::matrix<NekDouble> &vProjMat)
{
    ASSERTL0(vProjMat.size1() == m_movingFrameProjMat.size1(),
             "Matrices have different numbers of rows, cannot Set the "
             "moving frame projection matrix");
    ASSERTL0(vProjMat.size2() == m_movingFrameProjMat.size2(),
             "Matrices have different numbers of columns, cannot Set the "
             "moving frame projection matrix");
    for (int i = 0; i < vProjMat.size1(); ++i)
    {
        for (int j = 0; j < vProjMat.size2(); ++j)
        {
            m_movingFrameProjMat(i, j) = vProjMat(i, j);
        }
    }
}

void IncNavierStokes::v_GetMovingFrameProjectionMat(
    bnu::matrix<NekDouble> &vProjMat)
{
    ASSERTL0(vProjMat.size1() == m_movingFrameProjMat.size1(),
             "Matrices have different numbers of rows, cannot Get the "
             "moving frame projection matrix");
    ASSERTL0(vProjMat.size2() == m_movingFrameProjMat.size2(),
             "Matrices have different numbers of columns, cannot Get the "
             "moving frame projection matrix");

    for (int i = 0; i < vProjMat.size1(); ++i)
    {
        for (int j = 0; j < vProjMat.size2(); ++j)
        {
            vProjMat(i, j) = m_movingFrameProjMat(i, j);
        }
    }
}

/**
 * Function to set the angles between the moving frame of reference and
 * stationary inertial reference frame
 **/
void IncNavierStokes::v_SetMovingFrameAngles(
    const Array<OneD, NekDouble> &vFrameTheta)
{
    ASSERTL0(vFrameTheta.size() == m_movingFrameTheta.size(),
             "Arrays have different size, cannot set moving frame angles");
    for (int i = 0; i < vFrameTheta.size(); ++i)
    {
        m_movingFrameTheta[i] = vFrameTheta[i];
    }
}

/**
 * Function to get the angles between the moving frame of reference and
 * stationary inertial reference frame
 **/
void IncNavierStokes::v_GetMovingFrameAngles(
    Array<OneD, NekDouble> &vFrameTheta)
{
    ASSERTL0(vFrameTheta.size() == m_movingFrameTheta.size(),
             "Arrays have different size, cannot get moving frame angles");
    for (int i = 0; i < m_movingFrameTheta.size(); ++i)
    {
        vFrameTheta[i] = m_movingFrameTheta[i];
    }
}

/**
 * Function to check the type of forcing
 **/
bool IncNavierStokes::DefinedForcing(const std::string &sForce)
{
    vector<std::string> vForceList;
    bool hasForce{false};

    if (!m_session->DefinesElement("Nektar/Forcing"))
    {
        return hasForce;
    }

    TiXmlElement *vForcing = m_session->GetElement("Nektar/Forcing");
    if (vForcing)
    {
        TiXmlElement *vForce = vForcing->FirstChildElement("FORCE");
        while (vForce)
        {
            string vType = vForce->Attribute("TYPE");

            vForceList.push_back(vType);
            vForce = vForce->NextSiblingElement("FORCE");
        }
    }

    for (auto &f : vForceList)
    {
        if (boost::iequals(f, sForce))
        {
            hasForce = true;
        }
    }

    return hasForce;
}

/**
 * Get the pivot point for moving reference
 **/
void IncNavierStokes::GetPivotPoint(Array<OneD, NekDouble> &vPivotPoint)
{
    std::string sMRFForcingType = "MovingReferenceFrame";

    // only if we use moving reference frame formulation
    if (!DefinedForcing(sMRFForcingType))
    {
        return;
    }

    TiXmlElement *vForcing = m_session->GetElement("Nektar/Forcing");
    if (vForcing)
    {
        TiXmlElement *vForce = vForcing->FirstChildElement("FORCE");
        while (vForce)
        {
            string vType = vForce->Attribute("TYPE");

            // if it is moving reference frame
            if (boost::iequals(vType, sMRFForcingType))
            {
                TiXmlElement *pivotElmt =
                    vForce->FirstChildElement("PivotPoint");

                // if not defined, zero would be the default
                if (pivotElmt)
                {
                    std::stringstream pivotPointStream;
                    std::string pivotPointStr = pivotElmt->GetText();
                    pivotPointStream.str(pivotPointStr);

                    for (int j = 0; j < m_spacedim; ++j)
                    {
                        pivotPointStream >> pivotPointStr;
                        if (!pivotPointStr.empty())
                        {
                            LibUtilities::Equation equ(
                                m_session->GetInterpreter(), pivotPointStr);
                            vPivotPoint[j] = equ.Evaluate();
                        }
                    }
                }

                break;
            }
        }
    }
}

/**
 * Perform the extrapolation.
 */
bool IncNavierStokes::v_PreIntegrate(int step)
{
    m_extrapolation->SubStepSaveFields(step);
    m_extrapolation->SubStepAdvance(step, m_time);
    SetBoundaryConditions(m_time + m_timestep);
    return false;
}

} // namespace Nektar
