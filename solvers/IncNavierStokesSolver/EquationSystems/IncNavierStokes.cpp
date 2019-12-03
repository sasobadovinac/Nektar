///////////////////////////////////////////////////////////////////////////////
//
// File IncNavierStokes.cpp
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
// Description: Incompressible Navier Stokes class definition built on
// ADRBase class
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <boost/algorithm/string.hpp>

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Communication/Comm.h>
#include <SolverUtils/Filters/Filter.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>

#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

#include <tinyxml.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>

using namespace std;

namespace Nektar
{

    /**
     * Constructor. Creates ...
     *
     * \param
     * \param
     */
    IncNavierStokes::IncNavierStokes(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph):
        UnsteadySystem(pSession, pGraph),
        AdvectionSystem(pSession, pGraph),
        m_SmoothAdvection(false)
    {
    }

    void IncNavierStokes::v_InitObject()
    {
        AdvectionSystem::v_InitObject();

        int i,j;
        int numfields = m_fields.num_elements();
        std::string velids[] = {"u","v","w"};

        // Set up Velocity field to point to the first m_expdim of m_fields;
        m_velocity = Array<OneD,int>(m_spacedim);

        for(i = 0; i < m_spacedim; ++i)
        {
            for(j = 0; j < numfields; ++j)
            {
                std::string var = m_boundaryConditions->GetVariable(j);
                if(boost::iequals(velids[i], var))
                {
                    m_velocity[i] = j;
                    break;
                }

                ASSERTL0(j != numfields, "Failed to find field: " + var);
            }
        }

        // Set up equation type enum using kEquationTypeStr
        for(i = 0; i < (int) eEquationTypeSize; ++i)
        {
            bool match;
            m_session->MatchSolverInfo("EQTYPE",
                            kEquationTypeStr[i],match,false);
            if(match)
            {
                m_equationType = (EquationType)i;
                break;
            }
        }
        ASSERTL0(i != eEquationTypeSize,
                "EQTYPE not found in SOLVERINFO section");

        // This probably should to into specific implementations
        // Equation specific Setups
        switch(m_equationType)
        {
        case eSteadyStokes:
        case eSteadyOseen:
        case eSteadyNavierStokes:
        case eSteadyLinearisedNS:
        case eUnsteadyNavierStokes:
        case eUnsteadyStokes:
            break;
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type");
        }

        m_session->LoadParameter("Kinvis", m_kinvis);
        m_session->LoadParameter("NumSteps", m_numSteps);

        // Default advection type per solver
        std::string vConvectiveType;
        switch(m_equationType)
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
        m_advObject->InitObject( m_session, m_fields);

        // Forcing terms
        m_forcing = SolverUtils::Forcing::Load(m_session, shared_from_this(),
                                            m_fields, v_GetForceDimension());

        // check to see if any Robin boundary conditions and if so set
        // up m_field to boundary condition maps;
        m_fieldsBCToElmtID  = Array<OneD, Array<OneD, int> >(numfields);
        m_fieldsBCToTraceID = Array<OneD, Array<OneD, int> >(numfields);
        m_fieldsRadiationFactor  = 
                Array<OneD, Array<OneD, NekDouble> > (numfields);

        for (i = 0; i < m_fields.num_elements(); ++i)
        {
            bool Set = false;

            Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
            Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
            int radpts = 0;

            BndConds = m_fields[i]->GetBndConditions();
            BndExp   = m_fields[i]->GetBndCondExpansions();
            for(int n = 0; n < BndConds.num_elements(); ++n)
            {
                if(boost::iequals(BndConds[n]->GetUserDefined(),"Radiation"))
                {
                    ASSERTL0(BndConds[n]->GetBoundaryConditionType() ==
                        SpatialDomains::eRobin,
                        "Radiation boundary condition must be of type Robin <R>");

                    if(Set == false)
                    {
                        m_fields[i]->GetBoundaryToElmtMap(
                            m_fieldsBCToElmtID[i],m_fieldsBCToTraceID[i]);
                        Set = true;
                    }
                    radpts += BndExp[n]->GetTotPoints();
                }
                if(boost::iequals(BndConds[n]->GetUserDefined(),
                        "ZeroNormalComponent") || boost::istarts_with
                        (m_fields[i]->GetBndConditions()[n]->GetUserDefined(),
                        "BlowingSuction"))
                {
                    ASSERTL0(BndConds[n]->GetBoundaryConditionType() ==
                            SpatialDomains::eDirichlet,
                            "Zero Normal Component boundary condition option must be of type Dirichlet <D>");

                    if(Set == false)
                    {
                        m_fields[i]->GetBoundaryToElmtMap(
                            m_fieldsBCToElmtID[i],m_fieldsBCToTraceID[i]);
                        Set = true;
                    }
                }
            }

            m_fieldsRadiationFactor[i] = Array<OneD, NekDouble>(radpts);

            radpts = 0; // reset to use as a counter

            for(int n = 0; n < BndConds.num_elements(); ++n)
            {
                if(boost::iequals(BndConds[n]->GetUserDefined(),"Radiation"))
                {
                    int npoints    = BndExp[n]->GetNpoints();
                    Array<OneD, NekDouble> x0(npoints,0.0);
                    Array<OneD, NekDouble> x1(npoints,0.0);
                    Array<OneD, NekDouble> x2(npoints,0.0);
                    Array<OneD, NekDouble> tmpArray;

                    BndExp[n]->GetCoords(x0,x1,x2);

                    LibUtilities::Equation coeff =
                        std::static_pointer_cast<
                    SpatialDomains::RobinBoundaryCondition
                        >(BndConds[n])->m_robinPrimitiveCoeff;

                    coeff.Evaluate(x0,x1,x2,m_time,
                                   tmpArray = m_fieldsRadiationFactor[i]+ radpts);
                    //Vmath::Neg(npoints,tmpArray = m_fieldsRadiationFactor[i]+ radpts,1);
                    radpts += npoints;
                }
            }
        }

        // Set up maping for womersley BC - and load variables
        int var;
        for (int i = 0; i < m_fields.num_elements(); ++i)
        {
            for(int n = 0; n < m_fields[i]->GetBndConditions().num_elements(); ++n)
            {
                if(boost::istarts_with(m_fields[i]->GetBndConditions()[n]->GetUserDefined(),"Womersley"))
                {
                    // assumes that boundary condition is applied in normal direction
                    // and is decomposed for each direction. There could be a
                    // unique file for each direction
                    m_womersleyParams[i][n] = MemoryManager<WomersleyParams>::AllocateSharedPtr(m_spacedim);
                    // Read in fourier coeffs and precompute coefficients
                    SetUpWomersley(i, n,
                                   m_fields[i]->GetBndConditions()[n]->GetUserDefined());

                    m_fields[i]->GetBoundaryToElmtMap(m_fieldsBCToElmtID[i],m_fieldsBCToTraceID[i]);

                }
                else if(boost::istarts_with(m_fields[i]->GetBndConditions()[n]->GetUserDefined(),"BlowingSuction"))
                {
                    m_BlowingSuction = true;
                    var = i;
                }
            }
        }

        if(m_BlowingSuction)
        {
            m_bsbcParams = MemoryManager<BlowingSuctionParams>::AllocateSharedPtr();
            SetUpBlowingSuction(m_fields[var]->GetBndConditions()[0]->GetUserDefined());
        }

        // Set up Field Meta Data for output files
        m_fieldMetaDataMap["Kinvis"]   = boost::lexical_cast<std::string>(m_kinvis);
        m_fieldMetaDataMap["TimeStep"] = boost::lexical_cast<std::string>(m_timestep);

        m_finTime = m_timestep * m_numSteps;
        m_iteration = 1;
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
                const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, NekDouble> > &outarray)
    {
        int i;
        int VelDim     = m_velocity.num_elements();
        Array<OneD, Array<OneD, NekDouble> > velocity(VelDim);

        for(i = 0; i < VelDim; ++i)
        {
            velocity[i] = inarray[m_velocity[i]];
        }

        m_advObject->Advect(m_nConvectiveFields, m_fields,
                            velocity, inarray, outarray, m_time);
    }

    /**
     * Time dependent boundary conditions updating
     */
    void IncNavierStokes::SetBoundaryConditions(NekDouble time)
    {
        int i, n;
        std::string varName;
        int nvariables = m_fields.num_elements();

        for (i = 0; i < nvariables; ++i)
        {
            for(n = 0; n < m_fields[i]->GetBndConditions().num_elements(); ++n)
            {
                if(m_fields[i]->GetBndConditions()[n]->IsTimeDependent())
                {
                    varName = m_session->GetVariable(i);
                    m_fields[i]->EvaluateBoundaryConditions(time, varName);
                }
                else if(boost::istarts_with(
                          m_fields[i]->GetBndConditions()[n]->GetUserDefined(),
                          "Womersley"))
                {
                    SetWomersleyBoundary(i,n);
                }
            }

            // Set Radiation conditions if required
            SetRadiationBoundaryForcing(i);
        }

        SetZeroNormalVelocity();
    }

    /**
     * Probably should be pushed back into ContField?
     */
    void IncNavierStokes::SetRadiationBoundaryForcing(int fieldid)
    {
        int  i,n;

        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>                BndExp;

        BndConds = m_fields[fieldid]->GetBndConditions();
        BndExp   = m_fields[fieldid]->GetBndCondExpansions();

        StdRegions::StdExpansionSharedPtr elmt;
        StdRegions::StdExpansionSharedPtr Bc;

        int cnt;
        int elmtid,nq,offset, boundary;
        Array<OneD, NekDouble> Bvals, U;
        int cnt1 = 0;

        for(cnt = n = 0; n < BndConds.num_elements(); ++n)
        {
            std::string type = BndConds[n]->GetUserDefined();

            if((BndConds[n]->GetBoundaryConditionType() ==
                    SpatialDomains::eRobin) &&
                (boost::iequals(type,"Radiation")))
            {
                for(i = 0; i < BndExp[n]->GetExpSize(); ++i,cnt++)
                {
                    elmtid = m_fieldsBCToElmtID[m_velocity[fieldid]][cnt];
                    elmt   = m_fields[fieldid]->GetExp(elmtid);
                    offset = m_fields[fieldid]->GetPhys_Offset(elmtid);

                    U = m_fields[fieldid]->UpdatePhys() + offset;
                    Bc = BndExp[n]->GetExp(i);

                    boundary = m_fieldsBCToTraceID[fieldid][cnt];

                    // Get edge values and put into ubc
                    nq = Bc->GetTotPoints();
                    Array<OneD, NekDouble> ubc(nq);
                    elmt->GetTracePhysVals(boundary,Bc,U,ubc);

                    Vmath::Vmul(nq,&m_fieldsRadiationFactor[fieldid][cnt1 +
                                BndExp[n]->GetPhys_Offset(i)],1,
                                &ubc[0],1,&ubc[0],1);

                    Bvals = BndExp[n]->UpdateCoeffs()+BndExp[n]->
                                GetCoeff_Offset(i);

                    Bc->IProductWRTBase(ubc,Bvals);
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
        static bool Setup  = false;

        if(Setup == true)
        {
            return;
        }
        Setup = true;

        int  i,n;

        Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr > >
            BndConds(m_spacedim);
        Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >
            BndExp(m_spacedim);


        for(i = 0; i < m_spacedim; ++i)
        {
            BndConds[i] = m_fields[m_velocity[i]]->GetBndConditions();
            BndExp[i]   = m_fields[m_velocity[i]]->GetBndCondExpansions();
        }

        StdRegions::StdExpansionSharedPtr elmt,Bc;

        int cnt;
        int elmtid,nq, boundary;

        Array<OneD, Array<OneD, NekDouble> > normals;
        Array<OneD, NekDouble> Bphys,Bcoeffs;

        int fldid = m_velocity[0];

        for(cnt = n = 0; n < BndConds[0].num_elements(); ++n)
        {
            if((BndConds[0][n]->GetBoundaryConditionType() ==
                    SpatialDomains::eDirichlet) &&
                (boost::iequals(BndConds[0][n]->GetUserDefined(),
                    "ZeroNormalComponent")))
            {
                for(i = 0; i < BndExp[0][n]->GetExpSize(); ++i,cnt++)
                {
                    elmtid   = m_fieldsBCToElmtID[fldid][cnt];
                    elmt     = m_fields[0]->GetExp(elmtid);
                    boundary = m_fieldsBCToTraceID[fldid][cnt];

                    normals = elmt->GetSurfaceNormal(boundary);

                    nq = BndExp[0][n]->GetExp(i)->GetTotPoints();
                    Array<OneD, NekDouble> normvel(nq,0.0);

                    for(int k = 0; k < m_spacedim; ++k)
                    {
                        Bphys  = BndExp[k][n]->UpdatePhys()+
                            BndExp[k][n]->GetPhys_Offset(i);
                        Bc  = BndExp[k][n]->GetExp(i);
                        Vmath::Vvtvp(nq,normals[k],1,Bphys,1,normvel,1,
                                     normvel,1);
                    }

                    // negate normvel for next step
                    Vmath::Neg(nq,normvel,1);

                    for(int k = 0; k < m_spacedim; ++k)
                    {
                        Bphys  = BndExp[k][n]->UpdatePhys()+
                            BndExp[k][n]->GetPhys_Offset(i);
                        Bcoeffs = BndExp[k][n]->UpdateCoeffs()+
                            BndExp[k][n]->GetCoeff_Offset(i);
                        Bc  = BndExp[k][n]->GetExp(i);
                        Vmath::Vvtvp(nq,normvel,1,normals[k],1,Bphys,1,
                                     Bphys,1);
                        Bc->FwdTrans_BndConstrained(Bphys,Bcoeffs);
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
        int  i,j,k;

        int M_coeffs = WomParam->m_wom_vel.size();

        NekDouble T = WomParam->m_period;
        NekDouble axis_normal = WomParam->m_axisnormal[fldid];

        // Womersley Number
        NekComplexDouble omega_c (2.0*M_PI/T, 0.0);
        NekComplexDouble k_c (0.0, 0.0);
        NekComplexDouble m_time_c (m_time, 0.0);
        NekComplexDouble zi (0.0,1.0);
        NekComplexDouble i_pow_3q2 (-1.0/sqrt(2.0),1.0/sqrt(2.0));

        MultiRegions::ExpListSharedPtr  BndCondExp;
        BndCondExp   = m_fields[fldid]->GetBndCondExpansions()[bndid];

        StdRegions::StdExpansionSharedPtr bc;
        int cnt=0;
        int nfq;
        Array<OneD, NekDouble> Bvals;
        int exp_npts = BndCondExp->GetExpSize();
        Array<OneD, NekDouble> wbc(exp_npts,0.0);

        Array<OneD, NekComplexDouble> zt(M_coeffs);

        // preallocate the exponent
        for (k=1; k < M_coeffs; k++)
        {
            k_c  = NekComplexDouble((NekDouble) k, 0.0);
            zt[k] = std::exp(zi * omega_c * k_c * m_time_c);
        }

        // Loop over each element in an expansion
        for(i = 0; i < exp_npts; ++i,cnt++)
        {
            // Get Boundary and trace expansion
            bc = BndCondExp->GetExp(i);
            nfq = bc->GetTotPoints();
            Array<OneD, NekDouble> wbc(nfq,0.0);

            // Compute womersley solution
            for (j=0; j < nfq; j++)
            {
                wbc[j] = WomParam->m_poiseuille[i][j];
                for (k=1; k < M_coeffs; k++)
                {
                    zvel =  WomParam->m_zvel[i][j][k] * zt[k];
                    wbc[j] = wbc[j] + zvel.real();
                }
            }

            // Multiply w by normal to get u,v,w component of velocity
            Vmath::Smul(nfq,axis_normal,wbc,1,wbc,1);
            // get the offset
            Bvals = BndCondExp->UpdateCoeffs()+
                    BndCondExp->GetCoeff_Offset(i);

            // Push back to Coeff space
            bc->FwdTrans(wbc,Bvals);
        }
    }

    /**
     *  Womersley boundary condition init
     */
    void IncNavierStokes::SetUpWomersley(const int fldid, const int bndid, std::string womStr)
    {
        std::string::size_type indxBeg = womStr.find_first_of(':') + 1;
        string filename = womStr.substr(indxBeg,string::npos);

        NekComplexDouble coef;

        TiXmlDocument doc(filename);

        bool loadOkay = doc.LoadFile();
        ASSERTL0(loadOkay,(std::string("Failed to load file: ") +
                           filename).c_str());

        TiXmlHandle docHandle(&doc);

        int err;    /// Error value returned by TinyXML.

        TiXmlElement *nektar = doc.FirstChildElement("NEKTAR");
        ASSERTL0(nektar, "Unable to find NEKTAR tag in file.");

        TiXmlElement *wombc = nektar->FirstChildElement("WOMERSLEYBC");
        ASSERTL0(wombc, "Unable to find WOMERSLEYBC tag in file.");

        // read womersley parameters
        TiXmlElement *womparam = wombc->FirstChildElement("WOMPARAMS");
        ASSERTL0(womparam, "Unable to find WOMPARAMS tag in file.");

        // Input coefficients
        TiXmlElement *params = womparam->FirstChildElement("W");
        map<std::string,std::string> Wparams;

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

            std::transform(propstr.begin(),propstr.end(),propstr.begin(),
                           ::toupper);
            Wparams[propstr] = valstr;

            params = params->NextSiblingElement("W");
        }
        bool parseGood;

        // Read parameters

        ASSERTL0(Wparams.count("RADIUS") == 1,
          "Failed to find Radius parameter in Womersley boundary conditions");
        std::vector<NekDouble> rad;
        ParseUtils::GenerateVector(Wparams["RADIUS"],rad);
        m_womersleyParams[fldid][bndid]->m_radius = rad[0];

        ASSERTL0(Wparams.count("PERIOD") == 1,
          "Failed to find period parameter in Womersley boundary conditions");
        std::vector<NekDouble> period;
        parseGood = ParseUtils::GenerateVector(Wparams["PERIOD"],period);
        m_womersleyParams[fldid][bndid]->m_period = period[0];


        ASSERTL0(Wparams.count("AXISNORMAL") == 1,
          "Failed to find axisnormal parameter in Womersley boundary conditions");
        std::vector<NekDouble> anorm;
        parseGood = ParseUtils::GenerateVector(Wparams["AXISNORMAL"],anorm);
        m_womersleyParams[fldid][bndid]->m_axisnormal[0] = anorm[0];
        m_womersleyParams[fldid][bndid]->m_axisnormal[1] = anorm[1];
        m_womersleyParams[fldid][bndid]->m_axisnormal[2] = anorm[2];


        ASSERTL0(Wparams.count("AXISPOINT") == 1,
          "Failed to find axispoint parameter in Womersley boundary conditions");
        std::vector<NekDouble> apt;
        parseGood = ParseUtils::GenerateVector(Wparams["AXISPOINT"],apt);
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
            ASSERTL0(parseGood,
                    (std::string("Problem reading value of fourier coefficient, ID=") +
                    boost::lexical_cast<string>(indx)).c_str());
            ASSERTL1(coeffvals.size() == 2,
                    (std::string("Have not read two entries of Fourier coefficicent from ID="+
                    boost::lexical_cast<string>(indx)).c_str()));

            m_womersleyParams[fldid][bndid]->m_wom_vel.push_back(NekComplexDouble (coeffvals[0], coeffvals[1]));

            fval = fval->NextSiblingElement("F");
        }

        // starting point of precalculation
        int  i,j,k;
        // M fourier coefficients
        int M_coeffs = m_womersleyParams[fldid][bndid]->m_wom_vel.size();
        NekDouble R = m_womersleyParams[fldid][bndid]->m_radius;
        NekDouble T = m_womersleyParams[fldid][bndid]->m_period;
        Array<OneD, NekDouble > x0 = m_womersleyParams[fldid][bndid]->m_axispoint;

        NekComplexDouble rqR;
        // Womersley Number
        NekComplexDouble omega_c (2.0*M_PI/T, 0.0);
        NekComplexDouble alpha_c (R*sqrt(omega_c.real()/m_kinvis), 0.0);
        NekComplexDouble z1 (1.0,0.0);
        NekComplexDouble i_pow_3q2 (-1.0/sqrt(2.0),1.0/sqrt(2.0));

        MultiRegions::ExpListSharedPtr  BndCondExp;
        BndCondExp   = m_fields[fldid]->GetBndCondExpansions()[bndid];

        StdRegions::StdExpansionSharedPtr bc;
        int cnt = 0;
        int nfq;
        Array<OneD, NekDouble> Bvals;

        int exp_npts = BndCondExp->GetExpSize();
        Array<OneD, NekDouble> wbc(exp_npts,0.0);

        // allocate time indepedent variables
        m_womersleyParams[fldid][bndid]->m_poiseuille = Array<OneD, Array<OneD, NekDouble> > (exp_npts);
        m_womersleyParams[fldid][bndid]->m_zvel = Array<OneD, Array<OneD, Array<OneD, NekComplexDouble> > > (exp_npts);
        // could use M_coeffs - 1 but need to avoid complicating things
        Array<OneD, NekComplexDouble> zJ0(M_coeffs);
        Array<OneD, NekComplexDouble> lamda_n(M_coeffs);
        Array<OneD, NekComplexDouble> k_c(M_coeffs);
        NekComplexDouble zJ0r;

        for (k=1; k < M_coeffs; k++)
        {
            k_c[k]  = NekComplexDouble((NekDouble) k, 0.0);
            lamda_n[k] = i_pow_3q2 * alpha_c * sqrt(k_c[k]);
            zJ0[k]  = Polylib::ImagBesselComp(0,lamda_n[k]);
        }

        // Loop over each element in an expansion
        for(i = 0; i < exp_npts; ++i,cnt++)
        {
            // Get Boundary and trace expansion
            bc = BndCondExp->GetExp(i);
            nfq = bc->GetTotPoints();

            Array<OneD, NekDouble> x(nfq,0.0);
            Array<OneD, NekDouble> y(nfq,0.0);
            Array<OneD, NekDouble> z(nfq,0.0);
            bc->GetCoords(x,y,z);

            m_womersleyParams[fldid][bndid]->m_poiseuille[i] =
                    Array<OneD, NekDouble> (nfq);
            m_womersleyParams[fldid][bndid]->m_zvel[i] =
                    Array<OneD, Array<OneD, NekComplexDouble> > (nfq);

            // Compute coefficients
            for (j=0; j < nfq; j++)
            {
                rqR = NekComplexDouble (sqrt((x[j]-x0[0])*(x[j]-x0[0]) +
                         (y[j]-x0[1])*(y[j]-x0[1]) +
                         (z[j]-x0[2])*(z[j]-x0[2]))/R, 0.0);

                // Compute Poiseulle Flow
                m_womersleyParams[fldid][bndid]->m_poiseuille[i][j] =
                        m_womersleyParams[fldid][bndid]->m_wom_vel[0].real() *
                        (1. - rqR.real()*rqR.real());


                m_womersleyParams[fldid][bndid]->m_zvel[i][j] =
                        Array<OneD, NekComplexDouble> (M_coeffs);

                // compute the velocity information
                for (k=1; k < M_coeffs; k++)
                {
                    zJ0r = Polylib::ImagBesselComp(0,rqR * lamda_n[k]);
                    m_womersleyParams[fldid][bndid]->m_zvel[i][j][k] =
                            m_womersleyParams[fldid][bndid]->m_wom_vel[k] *
                            (z1 - (zJ0r / zJ0[k]));
                }
            }
        }
    }

    /**
     *  Solve structural equation and get m_angle and m_angleVel
     */
    void IncNavierStokes::SolveStructuralDirect(NekDouble time)
    {
        static int totalCalls = -1;
        ++totalCalls;

        // Determine order for this time step
        static int intCalls = 0;
        ++intCalls;
        int order = min(intCalls, m_bsbcParams->m_intSteps);

        int expdim = m_fields[0]->GetGraph()->GetMeshDimension();

        // Get aerodynamic forces
        Array<OneD, NekDouble> moments(expdim, 0.0);
        m_bsbcParams->m_filterForces->GetTotalMoments(m_fields, moments, time);

        // Shift moment storage
        for(int n = m_bsbcParams->m_intSteps-1; n > 0; --n)
        {
            m_bsbcParams->m_moment[n] = m_bsbcParams->m_moment[n-1];
        }

        // Calculate total moment
        if( expdim == 2)
        {
            // Only have moment in z-direction
            m_bsbcParams->m_moment[0] = moments[0];
        }
        else
        {
            // Project to axis direction
            m_bsbcParams->m_moment[0] = 0;
            for(int j = 0; j < 3; ++j)
            {
                m_bsbcParams->m_moment[0] += moments[j] * m_bsbcParams->m_axis[j];
            }
        }

        // Account for torsional spring and damping contributions
        m_bsbcParams->m_moment[0] -= m_bsbcParams->m_K * m_bsbcParams->m_angle
         + m_bsbcParams->m_C * m_bsbcParams->m_angleVel[0];

        // Shift velocity storage
        for(int n = m_bsbcParams->m_intSteps-1; n > 0; --n)
        {
            m_bsbcParams->m_angleVel[n] = m_bsbcParams->m_angleVel[n-1];
        }

        // Update velocity
        for(int j = 0; j < order; ++j)
        {
            m_bsbcParams->m_angleVel[0] += m_timestep *
                m_bsbcParams->AdamsBashforth_coeffs[order-1][j] 
                * m_bsbcParams->m_moment[j] / m_bsbcParams->m_I;
        }

        // Update position
        m_bsbcParams->m_previousAngle = m_bsbcParams->m_angle;
        for(int j = 0; j < order; ++j)
        {
            m_bsbcParams->m_angle += m_timestep *
                m_bsbcParams->AdamsMoulton_coeffs[order-1][j] * m_bsbcParams->m_angleVel[j];
        }

        if( m_bsbcParams->m_doOutput )
        {
            m_bsbcParams->m_filterForces->Update(m_fields, time);
            ++m_bsbcParams->m_index;
            if ( !(m_bsbcParams->m_index % m_bsbcParams->m_outputFrequency) )
            {
                if ( m_fields[0]->GetComm()->GetRank() == 0)
                {
                    m_bsbcParams->m_outputStream.width(8);
                    m_bsbcParams->m_outputStream << setprecision(6) << time;
                    m_bsbcParams->m_outputStream.width(15);
                    m_bsbcParams->m_outputStream << setprecision(8) << m_bsbcParams->m_angleVel[0];
                    m_bsbcParams->m_outputStream.width(15);
                    m_bsbcParams->m_outputStream << setprecision(8) << m_bsbcParams->m_angle;
                    m_bsbcParams->m_outputStream << endl;
                }
            }
        }
    }

    /**
     *  Solve structural equation and get m_angle and m_angleVel
     */
    void IncNavierStokes::SolveStructuralAdjoint(NekDouble time)
    {
        static int totalCalls = -1;
        ++totalCalls;

        ////////// TEST angleVel = sigma * exp(sigma * t) ///////
        NekDouble sigma = 6.7256e-02;
        NekDouble timeCopy = time;
        m_bsbcParams->m_iter += 1;

        while (timeCopy > m_finTime)
        {
            timeCopy -= m_finTime;
        }
        if ((m_bsbcParams->m_iter % m_numSteps) == 0)
        {
            m_bsbcParams->m_previousAngle = m_bsbcParams->m_angle;
        }
        if ((m_bsbcParams->m_iter % m_numSteps) == 1)
        {
            m_bsbcParams->m_alpha = m_bsbcParams->m_previousAngle/m_bsbcParams->m_angle;
        }
        m_bsbcParams->m_angleVel[0] = (1/m_bsbcParams->m_alpha) * sigma * exp(sigma * time);
        m_bsbcParams->m_angle       = (1/m_bsbcParams->m_alpha) * exp(sigma * time);
        ////////////////////////////////////////////////////////

        if( m_bsbcParams->m_doOutput )
        {
            //m_bsbcParams->m_filterForces->Update(m_fields, time);
            ++m_bsbcParams->m_index;
            if ( !(m_bsbcParams->m_index % m_bsbcParams->m_outputFrequency) )
            {
                if ( m_fields[0]->GetComm()->GetRank() == 0)
                {
                    m_bsbcParams->m_outputStream.width(8);
                    m_bsbcParams->m_outputStream << setprecision(6) << time;
                    m_bsbcParams->m_outputStream.width(15);
                    m_bsbcParams->m_outputStream << setprecision(8) << m_bsbcParams->m_angleVel[0];
                    m_bsbcParams->m_outputStream.width(15);
                    m_bsbcParams->m_outputStream << setprecision(8) << m_bsbcParams->m_angle;
                    m_bsbcParams->m_outputStream << endl;
                }
            }
        }
    }

    /**
     *  Scale boundary conditions by m_angle and m_angleVel
     */
    void IncNavierStokes::ScaleBSBCDirect()
    {
        int nfields = m_fields.num_elements();
        int nvel    = nfields - 1;
        int nbnds   = m_fields[0]->GetBndConditions().num_elements();

        // Declare variables
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr> BndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;

        // Loop boundary conditions looking for BS-BC's
        for(int n = 0 ; n < nbnds ; ++n)
        {
            if (m_bsbcParams->m_isBlowingSuction[0])
            {
                // Now, project result to boundary
                for (int i = 0; i < nfields; ++i)
                {
                    BndConds   = m_fields[i]->GetBndConditions();
                    BndExp     = m_fields[i]->GetBndCondExpansions();
                    if(boost::istarts_with(BndConds[n]->GetUserDefined(),"BlowingSuction"))
                    {
                        // Apply BSBC correction
                        if (i<nvel)
                        {
                            // Scale file with new m_angle and m_angleVel
                            // Compute base flow gradient times angle:
                            Vmath::Smul(BndExp[n]->GetTotPoints(), m_bsbcParams->m_angle, 
                                m_bsbcParams->m_deltaGradBnd[i], 1, 
                                BndExp[n]->UpdatePhys(), 1);

                            // Add angular velocity:
                            Vmath::Svtvp(BndExp[n]->GetTotPoints(), m_bsbcParams->m_angleVel[0], 
                                m_bsbcParams->m_deltaGammaBnd[i], 1, 
                                BndExp[n]->GetPhys(), 1, 
                                BndExp[n]->UpdatePhys(), 1);
                        }
                    }
                }
            }
        }
        // Finally, perform FwdTrans in all fields
        for (int i = 0; i < m_fields.num_elements(); ++i)
        {
            // Get boundary condition information
            BndConds   = m_fields[i]->GetBndConditions();
            BndExp     = m_fields[i]->GetBndCondExpansions();
            for(int n = 0 ; n < BndConds.num_elements(); ++n)
            {   
                if ( BndConds[n]->GetBoundaryConditionType() == 
                        SpatialDomains::eDirichlet)
                {
                    if( boost::istarts_with(BndConds[n]->GetUserDefined(),"BlowingSuction"))
                    {
                        BndExp[n]->FwdTrans_BndConstrained(BndExp[n]->GetPhys(),
                                                    BndExp[n]->UpdateCoeffs());
                        if (m_fields[i]->GetExpType() == MultiRegions::e3DH1D)
                        {
                            BndExp[n]->HomogeneousFwdTrans(BndExp[n]->GetCoeffs(),
                                                    BndExp[n]->UpdateCoeffs());
                        }
                    }
                }
            }
        }
    }

    /**
     *  Scale boundary conditions for Adjoint evolution operator
     */
    void IncNavierStokes::ScaleBSBCAdjoint()
    {
        int nfields = m_fields.num_elements();
        int nvel    = nfields - 1;
        int nbnds   = m_fields[0]->GetBndConditions().num_elements();

        // Declare variables
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr> BndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;

        // Loop boundary conditions looking for BS-BC's
        for(int n = 0 ; n < nbnds ; ++n)
        {
            if (m_bsbcParams->m_isBlowingSuction[0])
            {
                // Now, project result to boundary
                for (int i = 0; i < nfields; ++i)
                {
                    BndConds   = m_fields[i]->GetBndConditions();
                    BndExp     = m_fields[i]->GetBndCondExpansions();
                    if(boost::istarts_with(BndConds[n]->GetUserDefined(),"BlowingSuction"))
                    {
                        // Apply BSBC correction
                        if (i<nvel)
                        {
                            // Scale BC with m_angleVel times deltaGamma:
                            Vmath::Smul(BndExp[n]->GetTotPoints(), m_bsbcParams->m_angleVel[0], 
                                m_bsbcParams->m_deltaGammaBnd[i], 1, 
                                BndExp[n]->UpdatePhys(), 1);
                        }
                    }
                }
            }
        }
        // Finally, perform FwdTrans in all fields
        for (int i = 0; i < m_fields.num_elements(); ++i)
        {
            // Get boundary condition information
            BndConds   = m_fields[i]->GetBndConditions();
            BndExp     = m_fields[i]->GetBndCondExpansions();
            for(int n = 0 ; n < BndConds.num_elements(); ++n)
            {   
                if ( BndConds[n]->GetBoundaryConditionType() == 
                        SpatialDomains::eDirichlet)
                {
                    if( boost::istarts_with(BndConds[n]->GetUserDefined(),"BlowingSuction"))
                    {
                        BndExp[n]->FwdTrans_BndConstrained(BndExp[n]->GetPhys(),
                                                    BndExp[n]->UpdateCoeffs());
                        if (m_fields[i]->GetExpType() == MultiRegions::e3DH1D)
                        {
                            BndExp[n]->HomogeneousFwdTrans(BndExp[n]->GetCoeffs(),
                                                    BndExp[n]->UpdateCoeffs());
                        }
                    }
                }
            }
        }
    }


    /**
     *  Blowing suction boundary condition initialisation
     */
    void IncNavierStokes::SetUpBlowingSuction(std::string BSBCStr)
    {
        // Parameter list for creating FilterAeroForces
        std::map<std::string, std::string> vParams;
        std::stringstream                  sStream;

        // Enter file containing BC parameters
        std::string::size_type indxBeg = BSBCStr.find_first_of(':') + 1;
        string filename = BSBCStr.substr(indxBeg,string::npos);

        TiXmlDocument doc(filename);

        bool loadOkay = doc.LoadFile();
        ASSERTL0(loadOkay,(std::string("Failed to load file: ") +
                           filename).c_str());

        TiXmlHandle docHandle(&doc);

        int err;    /// Error value returned by TinyXML.

        TiXmlElement *nektar = doc.FirstChildElement("NEKTAR");
        ASSERTL0(nektar, "Unable to find NEKTAR tag in file.");

        TiXmlElement *bsbc = nektar->FirstChildElement("BLOWINGSUCTIONBC");
        ASSERTL0(bsbc, "Unable to find BLOWINGSUCTIONBC tag in file.");

        // read blowing suction parameters
        TiXmlElement *bsparam = bsbc->FirstChildElement("BLOWINGSUCTIONPARAMS");
        ASSERTL0(bsparam, "Unable to find BLOWINGSUCTIONPARAMS tag in file.");

        // Input coefficients
        TiXmlElement *params = bsparam->FirstChildElement("P");
        map<std::string,std::string> BSparams;

        // read parameter list
        while (params)
        {
            std::string propstr;
            propstr = params->Attribute("PROPERTY");

            ASSERTL0(!propstr.empty(),
                    "Failed to read PROPERTY value Blowing Suction BC Parameter");

            std::string valstr;
            valstr = params->Attribute("VALUE");

            ASSERTL0(!valstr.empty(),
                    "Failed to read VALUE value Blowing Suction BC Parameter");

            std::transform(propstr.begin(),propstr.end(),propstr.begin(),
                           ::toupper);
            BSparams[propstr] = valstr;

            params = params->NextSiblingElement("P");
        }
        bool parseGood;

        // Read parameters
        // Inertia
        ASSERTL0(BSparams.count("I") == 1,
          "Failed to find Inertia parameter in Blowing Suction boundary conditions");
        std::vector<NekDouble> mass;
        ParseUtils::GenerateVector(BSparams["I"],mass);
        m_bsbcParams->m_I = mass[0];
        
        // Damping
        ASSERTL0(BSparams.count("C") == 1,
          "Failed to find Damping parameter in Blowing Suction boundary conditions");
        std::vector<NekDouble> damp;
        ParseUtils::GenerateVector(BSparams["C"],damp);
        m_bsbcParams->m_C = damp[0];

        // Stiffness
        ASSERTL0(BSparams.count("K") == 1,
          "Failed to find Stiffness parameter in Blowing Suction boundary conditions");
        std::vector<NekDouble> stiff;
        ParseUtils::GenerateVector(BSparams["K"],stiff);
        m_bsbcParams->m_K = stiff[0];


        // Rotation axis
        m_bsbcParams->m_axis = Array<OneD, NekDouble> (3, 0.0);
        ASSERTL0(BSparams.count("AXIS") == 1,
          "Failed to find axis parameter in Blowing Suction boundary conditions");
        std::vector<NekDouble> axis;
        parseGood = ParseUtils::GenerateVector(BSparams["AXIS"],axis);
        m_bsbcParams->m_axis[0] = axis[0];
        m_bsbcParams->m_axis[1] = axis[1];
        m_bsbcParams->m_axis[2] = axis[2];

        // Output frequency 1st vParams
        ASSERTL0(BSparams.count("OUTPUTFREQUENCY") == 1,
          "Failed to find output frequency parameter in Blowing Suction boundary conditions");
        std::vector<NekDouble> freq;
        ParseUtils::GenerateVector(BSparams["OUTPUTFREQUENCY"],freq);
        m_bsbcParams->m_outputFrequency = freq[0];
        vParams["OutputFrequency"] = m_bsbcParams->m_outputFrequency;

        // File name 2nd vParams
        ASSERTL0(BSparams.count("OUTPUTFILE") == 1,
          "Failed to find output file name in Blowing Suction boundary conditions");
        std::vector<std::string> outfile;
        ParseUtils::GenerateVector(BSparams["OUTPUTFILE"],outfile);
        m_bsbcParams->m_outputFile = outfile[0];
        vParams["OutputFile"] = m_bsbcParams->m_outputFile;
        m_bsbcParams->m_doOutput = true;

        if (!(m_bsbcParams->m_outputFile.length() >= 4
              && m_bsbcParams->m_outputFile.substr(m_bsbcParams->m_outputFile.length() - 4) == ".mot"))
        {
            m_bsbcParams->m_outputFile += ".mot";
        }

        // Hinge point 3rd vParams
        m_bsbcParams->m_hingePoint = Array<OneD, NekDouble> (3, 0.0);
        ASSERTL0(BSparams.count("HINGEPOINT") == 1,
          "Failed to find hinge point parameter in Blowing Suction boundary conditions");
        std::vector<NekDouble> hinge;
        parseGood = ParseUtils::GenerateVector(BSparams["HINGEPOINT"],hinge);
        m_bsbcParams->m_hingePoint[0] = hinge[0];
        m_bsbcParams->m_hingePoint[1] = hinge[1];
        m_bsbcParams->m_hingePoint[2] = hinge[2];

        // Convert all but the last element to avoid a trailing ","
        std::copy(hinge.begin(), hinge.end()-1,std::ostream_iterator<int>(sStream, " "));
        // Now add the last element with no delimiter
        sStream << hinge.back();
        vParams["MomentPoint"] = sStream.str();

        // Boundary ID (modify to be able to read numerous boundary regions)
        ASSERTL0(BSparams.count("BOUNDARYID") == 1,
          "Failed to find Boundary parameter in Blowing Suction boundary conditions");
        std::vector<std::string> bnd;
        ParseUtils::GenerateVector(BSparams["BOUNDARYID"],bnd);
        m_bsbcParams->m_boundaryList = bnd[0];
        vParams["Boundary"]  = m_bsbcParams->m_boundaryList;

        // Create FilterAeroForces object 
        m_bsbcParams->m_filterForces = MemoryManager<SolverUtils::FilterAeroForces>::
        AllocateSharedPtr(m_session, shared_from_this(), vParams);
        m_bsbcParams->m_filterForces->Initialise(m_fields, 0.0);

        // Determine time integration order
        std::string intMethod = m_session->GetSolverInfo("TIMEINTEGRATIONMETHOD");
        if(boost::iequals(intMethod, "IMEXOrder1"))
        {
            m_bsbcParams->m_intSteps = 1;
        }
        else if (boost::iequals(intMethod, "IMEXOrder2"))
        {
            m_bsbcParams->m_intSteps = 2;
        }
        else if (boost::iequals(intMethod, "IMEXOrder3"))
        {
            m_bsbcParams->m_intSteps = 3;
        }
        else
        {
            ASSERTL0(false, "Time integration method not supported.");
        }

        // Initialise variables
        m_bsbcParams->m_angle         = 0.001;
        m_bsbcParams->m_previousAngle = 0.001;
        m_bsbcParams->m_angleVel      = Array<OneD,NekDouble> (m_intSteps, 0.001);
        m_bsbcParams->m_moment        = Array<OneD,NekDouble> (m_intSteps, 0.0);
        m_bsbcParams->m_alpha         = 1;
        m_bsbcParams->m_iter          = 0;

        // Initialise constants for the BC sclaing
        SetStructConsts();

        // Open outputstream and write header
        if( m_bsbcParams->m_doOutput)
        {
            m_bsbcParams->m_index = 0;
            if ( m_fields[0]->GetComm()->GetRank() == 0)
            {
                // Open output stream
                bool adaptive;
                m_session->MatchSolverInfo("Driver", "Adaptive", adaptive, false);
                if (adaptive)
                {
                    m_bsbcParams->m_outputStream.open(m_bsbcParams->m_outputFile.c_str(), ofstream::app);
                }
                else
                {
                    m_bsbcParams->m_outputStream.open(m_bsbcParams->m_outputFile.c_str());
                }
                // Write header
                m_bsbcParams->m_outputStream << "# Angular velocity and angle of hinged bodies"
                               << endl;
                m_bsbcParams->m_outputStream << "#" << " Hinge Position = " << " (";
                m_bsbcParams->m_outputStream.width(8);
                m_bsbcParams->m_outputStream << setprecision(4) << m_bsbcParams->m_hingePoint[0];
                m_bsbcParams->m_outputStream.width(8);
                m_bsbcParams->m_outputStream << setprecision(4) << m_bsbcParams->m_hingePoint[1];
                m_bsbcParams->m_outputStream.width(8);
                m_bsbcParams->m_outputStream << setprecision(4) << m_bsbcParams->m_hingePoint[2];
                m_bsbcParams->m_outputStream << ")" << endl;
                m_bsbcParams->m_outputStream << "#" << " Rotation  Axis = " << " (";
                m_bsbcParams->m_outputStream.width(8);
                m_bsbcParams->m_outputStream << setprecision(4) << m_bsbcParams->m_axis[0];
                m_bsbcParams->m_outputStream.width(8);
                m_bsbcParams->m_outputStream << setprecision(4) << m_bsbcParams->m_axis[1];
                m_bsbcParams->m_outputStream.width(8);
                m_bsbcParams->m_outputStream << setprecision(4) << m_bsbcParams->m_axis[2];
                m_bsbcParams->m_outputStream << ")" << endl;
                m_bsbcParams->m_outputStream << "# Boundary regions: "
                               << m_bsbcParams->m_boundaryList.c_str() << endl;
                m_bsbcParams->m_outputStream << "#";
                m_bsbcParams->m_outputStream.width(7);
                m_bsbcParams->m_outputStream << "Time";
                m_bsbcParams->m_outputStream.width(15);
                m_bsbcParams->m_outputStream <<  "Angular_Vel";
                m_bsbcParams->m_outputStream.width(15);
                m_bsbcParams->m_outputStream <<  "Angle";
                m_bsbcParams->m_outputStream << endl;
            }
        }
    }

    /**
    * Define values on the boundary of delta Gamma and 
    * delta Gamma times base flow gradient
    * for boundary conditions scaling
    */
    void IncNavierStokes::SetStructConsts()
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nfields = m_fields.num_elements();
        int nvel = nfields - 1;
        int nbnds    = m_fields[0]->GetBndConditions().num_elements();

        // Declare variables
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr> BndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;

        Array<OneD, Array<OneD, NekDouble> > deltaGrad(nvel);
        Array<OneD, Array<OneD, NekDouble> > deltaGamma(nvel);

        m_bsbcParams->m_deltaGammaBnd = Array<OneD, Array<OneD, NekDouble> >(nvel);
        m_bsbcParams->m_deltaGradBnd = Array<OneD, Array<OneD, NekDouble> >(nvel);

        m_bsbcParams->m_isBlowingSuction = Array<OneD, bool> (nfields);

        for (int i = 0; i < nfields; ++i)
        {
            m_bsbcParams->m_isBlowingSuction[i] = false;
        }

        for (int i = 0; i< nvel; i++)
        {
            deltaGrad[i] = Array<OneD, NekDouble> (physTot, 0.0);
            deltaGamma[i] = Array<OneD, NekDouble> (physTot, 0.0);
        }

        // Get Cartesian coordinates for evaluating boundary conditions    
        Array<OneD, Array<OneD, NekDouble> > coords(3);
        for (int dir=0; dir < 3; dir++)
        {
            coords[dir] = Array<OneD, NekDouble> (physTot, 0.0);
        }
        m_fields[0]->GetCoords(coords[0],coords[1],coords[2]);

        for(int j = 0; j < 3; ++j)
            {
                // Subtract m_momPoint
                Vmath::Sadd (physTot, -1.0*m_bsbcParams->m_hingePoint[j],
                             coords[j], 1, coords[j], 1);
            }

        // Loop boundary conditions looking for BS-BC's
        for(int n = 0 ; n < nbnds ; ++n)
        {
            // Evaluate original Dirichlet boundary conditions in whole domain
            for (int i = 0; i < nfields; ++i)
            {
                BndConds   = m_fields[i]->GetBndConditions();
                BndExp     = m_fields[i]->GetBndCondExpansions();
                if ( BndConds[n]->GetBoundaryConditionType() == 
                                    SpatialDomains::eDirichlet && 
                                    boost::istarts_with(BndConds[n]->GetUserDefined(),"BlowingSuction"))
                {
                    m_bsbcParams->m_isBlowingSuction[i] = true;
                    // If we have the a velocity component
                    //      check if all vel bc's are also Dirichlet
                    if ( i<nvel )
                    {
                        for (int j = 0; j < nvel; ++j)
                        {
                            ASSERTL0(m_fields[j]->GetBndConditions()[n]->
                                                  GetBoundaryConditionType() == 
                                                    SpatialDomains::eDirichlet,
                                "BSBC only supported when all velocity components have the same type of boundary conditions");
                        }
                    }
                    // Check if bc is time-dependent
                    ASSERTL0( !BndConds[n]->IsTimeDependent(),
                        "Time-dependent Dirichlet boundary conditions not supported with BSBC yet.");
                }
            }

            if (m_bsbcParams->m_isBlowingSuction[0])
            {
                // Compute delta Gamma:
                Vmath::Vcopy(physTot, coords[0], 1, deltaGamma[1], 1);
                Vmath::Neg(physTot, coords[1], 1);
                Vmath::Vcopy(physTot, coords[1], 1, deltaGamma[0], 1);
                Vmath::Neg(physTot, coords[1], 1);

                // Compute delta Gamma times base-flow gradient:
                // u = y.du/dx - x.du/dy
                Vmath::Vmul(physTot, coords[1], 1, m_advObject->GetGradBase()[0],
                    1 , deltaGrad[0], 1);

                Vmath::Neg(physTot, coords[0], 1);

                Vmath::Vvtvp(physTot, coords[0], 1, m_advObject->GetGradBase()[1],
                    1, deltaGrad[0], 1, deltaGrad[0], 1);

                // v = y.dv/dx - x.dv/dy
                Vmath::Vmul(physTot, coords[1], 1, m_advObject->GetGradBase()[2],
                    1 , deltaGrad[1], 1);

                Vmath::Vvtvp(physTot, coords[0], 1, m_advObject->GetGradBase()[3],
                    1, deltaGrad[1], 1, deltaGrad[1], 1);

                // Now, project result to boundary
                for (int i = 0; i < nfields; ++i)
                {
                    BndConds   = m_fields[i]->GetBndConditions();
                    BndExp     = m_fields[i]->GetBndCondExpansions();

                    if(boost::istarts_with(BndConds[n]->GetUserDefined(),"BlowingSuction"))
                    {
                        // Apply BSBC correction
                        if (i<nvel)
                        {
                            m_bsbcParams->m_deltaGammaBnd[i] = Array<OneD, NekDouble> (BndExp[n]->GetTotPoints());
                            m_bsbcParams->m_deltaGradBnd[i] = Array<OneD, NekDouble> (BndExp[n]->GetTotPoints());

                            // Get coordinate and values on boundary
                            m_fields[i]->ExtractPhysToBnd(n, deltaGrad[i], m_bsbcParams->m_deltaGradBnd[i]);
                            m_fields[i]->ExtractPhysToBnd(n, deltaGamma[i], m_bsbcParams->m_deltaGammaBnd[i]);
                        }
                    }
                }
            }
        }
    }

    /**
     *  
     */
    void IncNavierStokes::v_GetStruct(NekDouble &angle, NekDouble &angleVel)
    {
        angle = m_bsbcParams->m_angle;
        angleVel = m_bsbcParams->m_angleVel[0];
    }

    /**
     *  
     */
    void IncNavierStokes::v_SetStruct(NekDouble &angle, NekDouble &angleVel)
    {
        m_bsbcParams->m_angle = angle;
        m_bsbcParams->m_angleVel[0] = angleVel;
    }

    /**
     *  
     */
    bool IncNavierStokes::v_CheckBSBC()
    {
        return m_BlowingSuction;
    }

    /**
    * Add an additional forcing term programmatically.
    */
    void IncNavierStokes::AddForcing(
        const SolverUtils::ForcingSharedPtr& pForce)
    {
        m_forcing.push_back(pForce);
    }

    /**
     *
     */
    Array<OneD, NekDouble> IncNavierStokes::v_GetMaxStdVelocity(void)
    {
        int nvel  = m_velocity.num_elements();
        int nelmt = m_fields[0]->GetExpSize();

        Array<OneD, NekDouble> stdVelocity(nelmt, 0.0);
        Array<OneD, Array<OneD, NekDouble> > velfields;

        if(m_HomogeneousType == eHomogeneous1D) // just do check on 2D info
        {
            velfields = Array<OneD, Array<OneD, NekDouble> >(2);

            for(int i = 0; i < 2; ++i)
            {
                velfields[i] = m_fields[m_velocity[i]]->UpdatePhys();
            }
        }
        else
        {
            velfields = Array<OneD, Array<OneD, NekDouble> >(nvel);

            for(int i = 0; i < nvel; ++i)
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
    void IncNavierStokes::GetPressure(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD, NekDouble>                     &pressure)
    {
        pressure = physfield[m_nConvectiveFields];
    }

    /**
     *
     */
    void IncNavierStokes::GetDensity(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD, NekDouble>                     &density)
    {
        int nPts  = physfield[0].num_elements();
        Vmath::Fill(nPts, 1.0, density, 1);
    }

    /**
     *
     */
    void IncNavierStokes::GetVelocity(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> >       &velocity)
    {
        for(int i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = physfield[i];
        }
    }

    /**
     * Perform the extrapolation.
     */
    bool IncNavierStokes::v_PreIntegrate(int step)
    {
        m_extrapolation->SubStepSaveFields(step);
        m_extrapolation->SubStepAdvance(m_intSoln,step,m_time);
        SetBoundaryConditions(m_time+m_timestep);
        return false;
    }

} //end of namespace

