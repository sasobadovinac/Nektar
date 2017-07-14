////////////////////////////////////////////////////////////////////////////////
//
// File: CouplingCwipi.cpp
//
// For more information, please see: http://www.nektar.info/
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: CWIPI Exchange class
//
////////////////////////////////////////////////////////////////////////////////

#include "CouplingCwipi.h"

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/CsvIO.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/Foundations/PhysGalerkinProject.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>

#include <boost/functional/hash.hpp>

#define OUTPUT_FREQ 0

namespace Nektar
{
namespace SolverUtils
{

using namespace std;

std::string CouplingCwipi::className =
    GetCouplingFactory().RegisterCreatorFunction(
        "Cwipi", CouplingCwipi::create, "Cwipi Coupling");

void CouplingCwipi::InterpCallback(
    const int entities_dim,
    const int n_local_vertex,
    const int n_local_element,
    const int n_local_polhyedra,
    const int n_distant_point,
    const double local_coordinates[],
    const int local_connectivity_index[],
    const int local_connectivity[],
    const int local_polyhedra_face_index[],
    const int local_polyhedra_cell_to_face_connectivity[],
    const int local_polyhedra_face_connectivity_index[],
    const int local_polyhedra_face_connectivity[],
    const double distant_points_coordinates[],
    const int distant_points_location[],
    const float distant_points_distance[],
    const int distant_points_barycentric_coordinates_index[],
    const double distant_points_barycentric_coordinates[],
    const int stride,
    const cwipi_solver_type_t solver_type,
    const void *local_field,
    void *distant_field)
{
    Array<OneD, Array<OneD, NekDouble> > interpField(stride);

    Array<OneD, Array<OneD, NekDouble> > distCoords(n_distant_point);
    for (int i = 0; i < n_distant_point; ++i)
    {
        distCoords[i] = Array<OneD, NekDouble>(3);
        for (int j = 0; j < 3; ++j)
        {
            distCoords[i][j] = distant_points_coordinates[3 * i + j];
        }
    }

    std::stringstream sst;
    sst << entities_dim << "," << n_local_vertex << "," << stride;
    SendCallbackMap[sst.str()](interpField, distCoords);

    ASSERTL0(interpField.num_elements() == stride, "size mismatch");
    ASSERTL0(interpField[0].num_elements() == n_distant_point, "size mismatch");

    for (int i = 0; i < n_distant_point; i++)
    {
        for (int j = 0; j < stride; ++j)
        {
            ((double *)distant_field)[i * stride + j] = interpField[j][i];
        }
    }
}

CouplingCwipi::CouplingCwipi(MultiRegions::ExpListSharedPtr field): Coupling(field), m_sendHandle(-1), m_recvHandle(-1), m_lastSend(-1E6),m_lastReceive(-1E6), m_points(NULL), m_coords(NULL), m_connecIdx(NULL),
      m_connec(NULL), m_rValsInterl(NULL), m_sValsInterl(NULL)
{
    // defaults
    m_config["GEOMTOL"]          = "0.1";
    m_config["LOCALNAME"]        = "nektar";
    m_config["REMOTENAME"]       = "precise";
    m_config["OVERSAMPLE"]       = "0";
    m_config["FILTERWIDTH"]      = "-1";
    m_config["DUMPRAW"]          = "0";
    m_config["SENDMETHOD"]       = "NEARESTNEIGHBOUR";
    m_config["NOTLOCMETHOD"]     = "NOTOUCH";
}

void CouplingCwipi::v_Init()
{
    Coupling::v_Init();

    ReadConfig(m_evalField->GetSession());

    cwipi_add_local_int_control_parameter("nSendVars", m_nSendVars);
    cwipi_add_local_int_control_parameter("nRecvVars", m_nRecvVars);
    cwipi_add_local_string_control_parameter(
        "recvFieldNames", m_config["RECEIVEVARIABLES"].c_str());
    cwipi_add_local_string_control_parameter("sendFieldNames",
                                             m_config["SENDVARIABLES"].c_str());
    m_recvTag = boost::hash<std::string>()(m_couplingName + m_config["REMOTENAME"] + m_config["LOCALNAME"]) % INT_MAX;
    cwipi_add_local_int_control_parameter("receiveTag", m_recvTag);

    m_spacedim = m_evalField->GetGraph()->GetSpaceDimension();

    if (m_filtWidth > 0)
    {
        m_filtWidth = 2 * M_PI / m_filtWidth;
        m_filtWidth = m_filtWidth * m_filtWidth;
    }

    //  Init Coupling
    cwipi_solver_type_t solver_type = CWIPI_SOLVER_CELL_VERTEX;
    NekDouble geom_tol = boost::lexical_cast<NekDouble>(m_config["GEOMTOL"]);
    cwipi_create_coupling(m_couplingName.c_str(),
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                          m_config["REMOTENAME"].c_str(),
                          m_spacedim,
                          geom_tol,
                          CWIPI_STATIC_MESH,
                          solver_type,
                          OUTPUT_FREQ,
                          "Ensight Gold",
                          "text");
    cwipi_synchronize_control_parameter(m_config["REMOTENAME"].c_str());

    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cwipi_dump_application_properties();
    }

    m_sendTag = cwipi_get_distant_int_control_parameter(m_config["REMOTENAME"].c_str(), "receiveTag");

    if (cwipi_has_int_parameter(m_config["REMOTENAME"].c_str(), "nRecvVars"))
    {
        int remoteNRecvVars = cwipi_get_distant_int_control_parameter(m_config["REMOTENAME"].c_str(), "nRecvVars");
        ASSERTL0(remoteNRecvVars == m_nSendVars, "Number of local send vars different to remote received vars");
    }

    if (cwipi_has_int_parameter(m_config["REMOTENAME"].c_str(), "nSendVars"))
    {
        int remoteNSendVars = cwipi_get_distant_int_control_parameter(m_config["REMOTENAME"].c_str(), "nSendVars");
        ASSERTL0(remoteNSendVars == m_nRecvVars, "Number of local receive vars different to remote sent vars");
    }

    AnnounceMesh();

    if (m_nRecvVars > 0 and m_recvSteps > 0)
    {
        SetupReceive();
    }

    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "locating..." << endl;
    }
    cwipi_locate(m_couplingName.c_str());

    if (m_nSendVars > 0 and m_sendSteps > 0)
    {
        SetupSend();
    }

    ReceiveStart();
}

CouplingCwipi::~CouplingCwipi()
{
    free(m_coords);
    free(m_points);
    free(m_connec);
    free(m_connecIdx);
    free(m_rValsInterl);
    free(m_sValsInterl);
}

void CouplingCwipi::ReadConfig(LibUtilities::SessionReaderSharedPtr session)
{
    ASSERTL0(session->DefinesElement("Nektar/Coupling"),
             "No Coupling config found");

    m_config["LOCALNAME"] = session->GetCmdLineArgument<std::string>("cwipi");

    TiXmlElement *vCoupling = session->GetElement("Nektar/Coupling");
    ASSERTL0(vCoupling, "Invalid Coupling config");

    m_filtWidth = boost::lexical_cast<NekDouble>(m_config["FILTERWIDTH"]);
}

void CouplingCwipi::SetupReceive()
{
    int oversamp = boost::lexical_cast<int>(m_config["OVERSAMPLE"]);

    SpatialDomains::MeshGraphSharedPtr recvGraph =
        SpatialDomains::MeshGraph::Read(m_evalField->GetSession());
    recvGraph->SetExpansionsToPointOrder(
        oversamp + m_evalField->GetExp(0)->GetNumPoints(0));

    // TODO: DeclareCoeffPhysArrays
    switch (m_spacedim)
    {
        case 1:
        {
            m_recvField =
                MemoryManager<MultiRegions::ContField1D>::AllocateSharedPtr(
                    m_evalField->GetSession(), recvGraph, "DefaultVar");
            break;
        }

        case 2:
        {

            m_recvField =
                MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(
                    m_evalField->GetSession(), recvGraph);
            break;
        }

        case 3:
        {

            m_recvField =
                MemoryManager<MultiRegions::ContField3D>::AllocateSharedPtr(
                    m_evalField->GetSession(), recvGraph);
            break;
        }

        default:
        {
            ASSERTL0(false, "Expansion dimension not recognised");
            break;
        }
    }

    m_oldFields = Array<OneD, Array<OneD, NekDouble> >(m_nRecvVars);
    m_newFields = Array<OneD, Array<OneD, NekDouble> >(m_nRecvVars);
    for (int i = 0; i < m_nRecvVars; ++i)
    {
        m_oldFields[i] = Array<OneD, NekDouble>(m_evalField->GetTotPoints());
        m_newFields[i] = Array<OneD, NekDouble>(m_evalField->GetTotPoints());
    }

    // define the quadrature points at which we want to receive data
    m_nPoints = m_recvField->GetTotPoints();
    Array<OneD, Array<OneD, NekDouble> > coords(3);
    coords[0] = Array<OneD, NekDouble>(m_nPoints);
    coords[1] = Array<OneD, NekDouble>(m_nPoints);
    coords[2] = Array<OneD, NekDouble>(m_nPoints);
    m_recvField->GetCoords(coords[0], coords[1], coords[2]);

    m_points = (double *)malloc(sizeof(double) * 3 * m_nPoints);
    ASSERTL1(m_points != NULL, "malloc failed for m_points");

    for (int i = 0; i < m_nPoints; ++i)
    {
        m_points[3 * i + 0] = double(coords[0][i]);

        if (m_spacedim > 1)
        {
            m_points[3 * i + 1] = double(coords[1][i]);
        }
        else
        {
            m_points[3 * i + 1] = 0.0;
        }

        if (m_spacedim > 2)
        {
            m_points[3 * i + 2] = double(coords[2][i]);
        }
        else
        {
            m_points[3 * i + 2] = 0.0;
        }
    }

    cwipi_set_points_to_locate(m_couplingName.c_str(), m_nPoints, m_points);

    m_rValsInterl = (double *)malloc(sizeof(double) * m_nPoints * m_nRecvVars);
    ASSERTL1(m_rValsInterl != NULL, "malloc failed for m_rValsInterl");
}

void CouplingCwipi::SetupSend()
{
    // this array is never used because of our send callback method
    m_sValsInterl = (double *)malloc(
        sizeof(double) * m_evalField->GetGraph()->GetNvertices() * m_nSendVars);
    ASSERTL1(m_sValsInterl != NULL, "malloc failed for m_sValsInterl");
    for (int i = 0; i < m_evalField->GetGraph()->GetNvertices() * m_nSendVars;
         ++i)
    {
        m_sValsInterl[i] = i;
    }

    // register this coupling as sender
    std::stringstream sst;
    sst << m_spacedim << "," << m_evalField->GetGraph()->GetNvertices() << ","
        << m_nSendVars;
    SendCallbackMap[sst.str()] = boost::bind(&CouplingCwipi::SendCallback, this, _1, _2);
    cwipi_set_interpolation_function(m_couplingName.c_str(), CouplingCwipi::InterpCallback);
}

void CouplingCwipi::EvaluateFields(
    Array<OneD, Array<OneD, NekDouble> > interpField,
    Array<OneD, Array<OneD, NekDouble> > distCoords)
{
    int nOutPts = distCoords.num_elements();

    Array<OneD, NekDouble> Lcoords(m_spacedim, 0.0);
    for (int i = 0; i < nOutPts; ++i)
    {
        // Obtain Element and LocalCoordinate to interpolate
        int elmtid    = -1;
        NekDouble tol = NekConstants::kNekZeroTol;
        while (elmtid < 0 and tol <= 1E3 * NekConstants::kNekZeroTol)
        {
            elmtid = m_evalField->GetExpIndex(distCoords[i], Lcoords, tol);
            tol *= 2;
        }
        if (tol > 2 * NekConstants::kNekZeroTol)
        {
            for (int j = 0; j < m_spacedim; ++j)
            {
                if (Lcoords[j] < -1 - 0.75 * NekConstants::kNekZeroTol)
                {
                    Lcoords[j] = -1;
                }
                if (Lcoords[j] > 1 + 0.75 * NekConstants::kNekZeroTol)
                {
                    Lcoords[j] = 1;
                }
            }
        }

        ASSERTL0(elmtid >= 0,
                 "no element found for (" +
                     boost::lexical_cast<string>(distCoords[i][0]) + ", " +
                     boost::lexical_cast<string>(distCoords[i][1]) + ", " +
                     boost::lexical_cast<string>(distCoords[i][2]) + ")");

        int offset =
            m_evalField->GetPhys_Offset(m_evalField->GetOffset_Elmt_Id(elmtid));

        for (int f = 0; f < m_nSendVars; ++f)
        {
            NekDouble value = m_evalField->GetExp(elmtid)->StdPhysEvaluate(
                Lcoords, m_sendField[f] + offset);

            ASSERTL0(!(boost::math::isnan)(value), "new value is not a number");
            interpField[f][i] = value;
        }
    }
}

void CouplingCwipi::SetupSendInterpolation()
{
    const double *distCoords =
        cwipi_get_distant_coordinates(m_couplingName.c_str());
    int nPts = cwipi_get_n_distant_points(m_couplingName.c_str());

    Array<OneD, Array<OneD, NekDouble> > local(3);
    for (int i = 0; i < 3; ++i)
    {
        local[i] = Array<OneD, NekDouble>(m_evalField->GetTotPoints(), 0.0);
    }
    m_evalField->GetCoords(local[0], local[1], local[2]);
    LibUtilities::PtsFieldSharedPtr locatPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(3, local);

    Array<OneD, Array<OneD, NekDouble> > dist(3);
    for (int i = 0; i < 3; ++i)
    {
        dist[i] = Array<OneD, NekDouble>(nPts);
        for (int j = 0; j < nPts; ++j)
        {
            dist[i][j] = distCoords[3 * j + i];
        }
    }
    LibUtilities::PtsFieldSharedPtr distPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(3, dist);

    FieldUtils::InterpMethod method = FieldUtils::eNearestNeighbour;
    if (boost::to_upper_copy(m_config["SENDMETHOD"]) == "SHEPARD")
    {
        method = FieldUtils::eShepard;
    }
    m_sendInterpolator =
        MemoryManager<FieldUtils::Interpolator>::AllocateSharedPtr(method);
    m_sendInterpolator->CalcWeights(locatPts, distPts);
    m_sendInterpolator->PrintStatistics();
}

void CouplingCwipi::AnnounceMesh()
{
    SpatialDomains::MeshGraphSharedPtr graph = m_evalField->GetGraph();

    // get Elements
    SpatialDomains::SegGeomMap seggeom;
    SpatialDomains::TriGeomMap trigeom;
    SpatialDomains::QuadGeomMap quadgeom;
    SpatialDomains::TetGeomMap tetgeom;
    SpatialDomains::PyrGeomMap pyrgeom;
    SpatialDomains::PrismGeomMap prismgeom;
    SpatialDomains::HexGeomMap hexgeom;
    if (m_spacedim == 1)
    {
        seggeom = graph->GetAllSegGeoms();
    }
    else if (m_spacedim == 2)
    {
        trigeom  = graph->GetAllTriGeoms();
        quadgeom = graph->GetAllQuadGeoms();
    }
    else if (m_spacedim == 3)
    {
        tetgeom   = graph->GetAllTetGeoms();
        pyrgeom   = graph->GetAllPyrGeoms();
        prismgeom = graph->GetAllPrismGeoms();
        hexgeom   = graph->GetAllHexGeoms();
    };

    int nVerts = graph->GetNvertices();
    int nElts = seggeom.size() + trigeom.size() + quadgeom.size() +
                tetgeom.size() + pyrgeom.size() + prismgeom.size() +
                hexgeom.size();

    // allocate CWIPI arrays
    m_coords = (double *)malloc(sizeof(double) * 3 * nVerts);
    ASSERTL1(m_coords != NULL, "malloc failed for m_coords");
    int tmp = 2 * seggeom.size() + 3 * trigeom.size() + 4 * quadgeom.size() +
              4 * tetgeom.size() + 5 * pyrgeom.size() + 6 * prismgeom.size() +
              8 * hexgeom.size();
    m_connec = (int *)malloc(sizeof(int) * tmp);
    ASSERTL1(m_connec != NULL, "malloc failed for m_connec");
    m_connecIdx = (int *)malloc(sizeof(int) * (nElts + 1));
    ASSERTL1(m_connecIdx != NULL, "malloc failed for m_connecIdx");

    m_connecIdx[0] = 0;
    int coordsPos  = 0;
    int connecPos  = 0;
    int conidxPos  = 0;

    AddElementsToMesh(seggeom, coordsPos, connecPos, conidxPos);
    AddElementsToMesh(trigeom, coordsPos, connecPos, conidxPos);
    AddElementsToMesh(quadgeom, coordsPos, connecPos, conidxPos);
    AddElementsToMesh(tetgeom, coordsPos, connecPos, conidxPos);
    AddElementsToMesh(pyrgeom, coordsPos, connecPos, conidxPos);
    AddElementsToMesh(prismgeom, coordsPos, connecPos, conidxPos);
    AddElementsToMesh(hexgeom, coordsPos, connecPos, conidxPos);

    // output the mesh in tecplot format. If this works, CWIPI will be able
    // to process it, too
    /*
    cout << "VARIABLES = \"X\", \"Y\", \"Z\", \"U\"" <<  endl;
    cout << "ZONE NODES=" <<  nVerts << ",ELEMENTS=" <<  nElts;
    cout << ",DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON" <<  endl;
    for (int i = 0; i < nVerts; ++i)
    {
        cout << m_coords[3*i + 0] << " " << m_coords[3*i + 1] << " ";
        cout << m_coords[3*i + 2] << " " <<  1.0 << endl;
    }
    for (int i = 0; i < nElts; ++i)
    {
        cout << m_connec[i*4 + 0] << " " << m_connec[i*4 + 1] << " ";
        cout << m_connec[i*4 + 2] <<  " " <<  m_connec[i*4 + 3] <<  endl;
    }
    */

    cwipi_define_mesh(
        m_couplingName.c_str(), nVerts, nElts, m_coords, m_connecIdx, m_connec);
}

void CouplingCwipi::v_Finalize(void)
{
    cwipi_delete_coupling(m_couplingName.c_str());
}

template <typename T>
void CouplingCwipi::AddElementsToMesh(T geom,
                                      int &coordsPos,
                                      int &connecPos,
                                      int &conidxPos)
{
    // helper variables
    Array<OneD, NekDouble> x(3);
    SpatialDomains::PointGeomSharedPtr vert;
    int vertID;

    int kNverts = T::mapped_type::element_type::kNverts;

    // iterate over all elements
    typename T::iterator it;

    for (it = geom.begin(); it != geom.end(); it++)
    {
        //  iterate over the elements vertices
        for (int j = 0; j < kNverts; ++j)
        {
            vert   = it->second->GetVertex(j);
            vertID = vert->GetVid();

            // check if we already stored the vertex
            if (m_vertMap.count(vertID) == 0)
            {
                //  store the vertex
                vert->GetCoords(x[0], x[1], x[2]);
                m_coords[3 * coordsPos + 0] = double(x[0]);
                m_coords[3 * coordsPos + 1] = double(x[1]);
                m_coords[3 * coordsPos + 2] = double(x[2]);

                // store the vertex position in the m_coords array
                m_vertMap[vertID] = coordsPos;
                coordsPos++;
            }

            m_connec[connecPos] = m_vertMap[vertID] + 1;
            connecPos++;
        }

        m_connecIdx[conidxPos + 1] = m_connecIdx[conidxPos] + kNverts;
        conidxPos++;
    }
}


void CouplingCwipi::SendCallback(
    Array<OneD, Array<OneD, NekDouble> > &interpField,
    Array<OneD, Array<OneD, NekDouble> > &distCoords)
{
    ASSERTL0(interpField.num_elements() == m_nSendVars, "size mismatch");

    for (int i = 0; i < m_nSendVars; ++i)
    {
        interpField[i] = Array<OneD, NekDouble>(distCoords.num_elements());
    }

    if (boost::to_upper_copy(m_config["SENDMETHOD"]) == "NEARESTNEIGHBOUR" ||
        boost::to_upper_copy(m_config["SENDMETHOD"]) == "SHEPARD"
    )
    {
        if (not m_sendInterpolator)
        {
            SetupSendInterpolation();
        }

        LibUtilities::PtsFieldSharedPtr ptsIn =
            MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
                0, m_sendField);
        LibUtilities::PtsFieldSharedPtr ptsOut =
            MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
                0, interpField);
        m_sendInterpolator->Interpolate(ptsIn, ptsOut);
    }
    else
    {
        EvaluateFields(interpField, distCoords);
    }
}

void CouplingCwipi::v_Send(
    const int step,
    const NekDouble time,
    const Array<OneD, const Array<OneD, NekDouble> > &field,
    LibUtilities::FieldMetaDataMap &fieldMetaDataMap)
{
    if (m_nSendVars < 1 or m_sendSteps < 1)
    {
        return;
    }

    if (step >= m_lastSend + m_sendSteps)
    {
        SendComplete();

        m_lastSend = step;

        string tmp = fieldMetaDataMap["Variables"] + fieldMetaDataMap["AuxVariables"];
        vector<string> vars;
        ParseUtils::GenerateOrderedStringVector(tmp.c_str(), vars);
        vector<int> sendVarsToVars = GenerateVariableMapping(vars, m_sendFieldNames);
        m_sendField = Array<OneD, Array<OneD, NekDouble> > (m_nSendVars);
        for (int i = 0; i < sendVarsToVars.size(); ++i)
        {
            m_sendField[i] = field[sendVarsToVars[i]];
        }

        if (m_evalField->GetComm()->GetRank() == 0 &&
            m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
        {
            cout << "sending fields at i = " << step << ", t = " << time
                 << endl;
        }

        Timer timer1;
        timer1.Start();

        char sendFN[10];
        strcpy(sendFN, "dummyName");

        cwipi_issend(m_couplingName.c_str(),
                       "ex1",
                       m_sendTag,
                       m_nSendVars,
                       step,
                       time,
                       sendFN,
                       m_sValsInterl,
                       &m_sendHandle);
        timer1.Stop();

        if (m_evalField->GetComm()->GetRank() == 0 &&
            m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
        {
            cout << "Send total time: " << timer1.TimePerTest(1) << endl;
        }
    }
}


void CouplingCwipi::SendComplete()
{
    if (m_sendHandle < 0)
    {
        return;
    }

    Timer timer1;
    timer1.Start();
    cwipi_wait_issend(m_couplingName.c_str(), m_sendHandle);
    timer1.Stop();

    // set to -1 so we dont try finishing a send before a new one was started
    m_sendHandle = -1;

    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "Send waiting time: " << timer1.TimePerTest(1) << endl;
    }
}


void CouplingCwipi::ReceiveStart()
{
    if (m_recvHandle >= 0)
    {
        return;
    }

    Timer timer1;
    timer1.Start();
    // workaround a bug in cwipi: receiving_field_name should be const char* but
    // is char*
    char recFN[10];
    strcpy(recFN, "dummyName");

    cwipi_irecv(m_couplingName.c_str(),
                "ex1",
                m_recvTag,
                m_nRecvVars,
                0,
                0.0,
                recFN,
                m_rValsInterl,
                &m_recvHandle);
    timer1.Stop();

    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "Receive start time: " << timer1.TimePerTest(1) << endl;
    }
}

void CouplingCwipi::v_Receive(const int step,
                                  const NekDouble time,
                                  Array<OneD, Array<OneD, NekDouble> > &field,
                                  LibUtilities::FieldMetaDataMap &fieldMetaDataMap)
{
    if (m_nRecvVars < 1 or m_recvSteps < 1)
    {
        return;
    }

    Array<OneD, Array<OneD, NekDouble> > recvFields(m_nRecvVars);
    string tmp = fieldMetaDataMap["Variables"] + fieldMetaDataMap["AuxVariables"];
    vector<string> vars;
    ParseUtils::GenerateOrderedStringVector(tmp.c_str(), vars);
    vector<int> recvVarsToVars = GenerateVariableMapping(vars, m_recvFieldNames);
    ASSERTL1(m_nRecvVars == recvVarsToVars.size(), "field size mismatch");
    for (int i = 0; i < recvVarsToVars.size(); ++i)
    {
        recvFields[i] = field[recvVarsToVars[i]];
    }

    int nq = m_evalField->GetTotPoints();

    // make sure we have sensible data in old/new recvFields the first time this
    // method is called
    if (m_lastReceive < 0)
    {
        for (int i = 0; i < m_nRecvVars; ++i)
        {
            Vmath::Vcopy(nq, recvFields[i], 1, m_oldFields[i], 1);
            Vmath::Vcopy(nq, recvFields[i], 1, m_newFields[i], 1);
        }
    }

    if (step >= m_lastReceive + m_recvSteps)
    {
        for (int i = 0; i < m_nRecvVars; ++i)
        {
            Vmath::Vcopy(nq, m_newFields[i], 1, m_oldFields[i], 1);
        }

        ReceiveCwipi(step, time, m_newFields);
    }

    NekDouble fact =
        NekDouble(step - m_lastReceive + 1) / NekDouble(m_recvSteps);
    for (int i = 0; i < m_nRecvVars; ++i)
    {
        Vmath::Svtsvtp(nq,
                       fact,
                       m_newFields[i],
                       1,
                       (1 - fact),
                       m_oldFields[i],
                       1,
                       recvFields[i],
                       1);
    }
}

void CouplingCwipi::ReceiveCwipi(const int step,
                            const NekDouble time,
                            Array<OneD, Array<OneD, NekDouble> > &field)
{
    ASSERTL1(m_nRecvVars == field.num_elements(), "field size mismatch");

    if (m_nRecvVars < 1 or m_recvSteps < 1)
    {
        return;
    }

    if (step >= m_lastReceive + m_recvSteps)
    {
        m_lastReceive = step;

        if (m_evalField->GetComm()->GetRank() == 0 &&
            m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
        {
            cout << "waiting for receive at i = " << step << ", t = " << time << endl;
        }

        Timer timer1, timer2, timer3;
        timer1.Start();

        Array<OneD, NekDouble> tmpC(m_recvField->GetNcoeffs());
        Array<OneD, Array<OneD, NekDouble> > rVals(m_nRecvVars);
        for (int i = 0; i < m_nRecvVars; ++i)
        {
            rVals[i] = Array<OneD, NekDouble>(m_recvField->GetTotPoints());
        }

        timer2.Start();
        cwipi_wait_irecv(m_couplingName.c_str(), m_recvHandle);
        timer2.Stop();

        // set to -1 so we know we can start receiving again
        m_recvHandle = -1;


        int nNotLoc = cwipi_get_n_not_located_points(m_couplingName.c_str());
        Array<OneD, int> notLoc;
        if (nNotLoc != 0)
        {
            cout << "WARNING: relocating " << nNotLoc << " of " << m_nPoints
                << " points" << endl;

            const int *tmp = cwipi_get_not_located_points(m_couplingName.c_str());
            notLoc = Array<OneD, int>(nNotLoc);
            for (int i = 0; i < nNotLoc; ++i)
            {
                notLoc[i] = tmp[i] - 1;
            }

            if (boost::to_upper_copy(m_config["NOTLOCMETHOD"]) == "NOTOUCH")
            {
                // interpolate from m_evalField to m_recvField
                for (int i = 0; i < m_nRecvVars; ++i)
                {
                    m_evalField->FwdTrans(field[i], tmpC);
                    m_recvField->BwdTrans(tmpC, rVals[i]);
                }
            }
        }


        for (int i = 0, locPos = 0, intPos = 0; i < m_nPoints; ++i)
        {
            if (locPos < nNotLoc && notLoc[locPos] == i)
            {
                // keep the original value of field[j][i]
                locPos++;
            }
            else
            {
                for (int j = 0; j < m_nRecvVars; ++j)
                {
                    rVals[j][i] = m_rValsInterl[intPos * m_nRecvVars + j];
                }
                intPos++;
            }
        }

        if (boost::to_upper_copy(m_config["NOTLOCMETHOD"]) == "EXTRAPOLATE")
        {
            int doExtrapolate = 0;
            if (nNotLoc != 0)
            {
                doExtrapolate = 1;
            }
            m_evalField->GetSession()->GetComm()->AllReduce(doExtrapolate, LibUtilities::ReduceMax);
            if (doExtrapolate > 0)
            {
                ExtrapolateFields(rVals, notLoc);
            }
        }

        OverrrideFields(rVals);

        if (m_config["DUMPRAW"] != "0")
        {
            DumpRawFields(time, rVals);
        }

        if (m_filtWidth > 0)
        {
            for (int i = 0; i < m_nRecvVars; ++i)
            {
                timer3.Start();

                Array<OneD, NekDouble> forcing(m_nPoints);

                Array<OneD, Array<OneD, NekDouble> > Velocity(m_spacedim);
                for (int j = 0; j < m_spacedim; ++j)
                {
                    Velocity[j] = Array<OneD, NekDouble>(m_nPoints, 0.0);
                }

                Vmath::Smul(m_nPoints, -m_filtWidth, rVals[i], 1, forcing, 1);

                // Note we are using the
                // LinearAdvectionDiffusionReaction solver here
                // instead of HelmSolve since m_filtWidth is negative and
                // so matrices are not positive definite. Ideally
                // should allow for negative m_filtWidth coefficient in
                // HelmSolve
                m_recvField->LinearAdvectionDiffusionReactionSolve(
                    Velocity, forcing, tmpC, -m_filtWidth);

                m_evalField->BwdTrans(tmpC, field[i]);
                timer3.Stop();

                if (m_evalField->GetComm()->GetRank() == 0 &&
                    m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
                {
                    cout << "Smoother time (" << m_recvFieldNames[i]
                        << "): " << timer3.TimePerTest(1) << endl;
                }
            }
        }
        else
        {
            for (int i = 0; i < m_nRecvVars; ++i)
            {
                m_recvField->FwdTrans(rVals[i], tmpC);
                m_evalField->BwdTrans(tmpC, field[i]);
            }
        }
        timer1.Stop();

        if (m_evalField->GetComm()->GetRank() == 0 &&
            m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
        {
            cout << "Receive total time: " << timer1.TimePerTest(1) << ", ";
            cout << "Receive waiting time: " << timer2.TimePerTest(1) << endl;
        }

        ReceiveStart();
    }
}

void CouplingCwipi::OverrrideFields(Array<OneD, Array<OneD, NekDouble> > &rVals)
{
    if (m_evalField->GetSession()->GetSessionName() == "CESAM-HP-single")
    {
        // HACK
        Array<OneD, Array<OneD, NekDouble> > x(3);
        x[0] = Array<OneD, NekDouble>(m_recvField->GetTotPoints(), 0.0);
        x[1] = Array<OneD, NekDouble>(m_recvField->GetTotPoints(), 0.0);
        x[2] = Array<OneD, NekDouble>(m_recvField->GetTotPoints(), 0.0);
        m_recvField->GetCoords(x[0], x[1], x[2]);
        NekDouble distsq = 1E23;
        NekDouble distsqNew = distsq;
        int refID = -1;
        for (int i = 0; i < m_recvField->GetTotPoints(); ++i)
        {
            distsqNew = pow(x[0][i] - 0.135, NekDouble(2)) + pow(x[1][i], NekDouble(2)) + pow(x[2][i], NekDouble(2));
            if (distsqNew < distsq)
            {
                distsq = distsqNew;
                refID = i;
            }
        }
        ASSERTL0(refID >= 0, "refID < 0");

        int ownerRank = -1;
        distsqNew = distsq;
        m_evalField->GetSession()->GetComm()->AllReduce(distsqNew, LibUtilities::ReduceMin);
        if(abs(distsqNew - distsq) <= NekConstants::kNekSqrtTol)
        {
            // this rank has the closest point
            ownerRank = m_evalField->GetSession()->GetComm()->GetRank();
        }
        m_evalField->GetSession()->GetComm()->AllReduce(ownerRank, LibUtilities::ReduceMax);
        ASSERTL0(ownerRank >= 0, "ownerRank < 0");


        Array<OneD, NekDouble> data(m_nRecvVars, 0.0);
        if (ownerRank == m_evalField->GetSession()->GetComm()->GetRank())
        {
            for (int j = 0; j < m_nRecvVars; ++j)
            {
                data[j] = rVals[j][refID];
            }
        }
        m_evalField->GetSession()->GetComm()->AllReduce(data, LibUtilities::ReduceSum);


        for (int i = 0; i < m_recvField->GetTotPoints(); ++i)
        {
            NekDouble tmp = pow(x[0][i] - 0.173, NekDouble(2)) + pow(x[1][i], NekDouble(2)) + pow(x[2][i], NekDouble(2));
            if (tmp <= pow(0.173-0.135, 2))
            {
                for (int j = 0; j < m_nRecvVars; ++j)
                {
                    rVals[j][i] = data[j];
                }
            }
        }
    }
    else if (m_evalField->GetSession()->GetSessionName() == "chakra")
    {
        // HACK
        Array<OneD, Array<OneD, NekDouble> > x(3);
        x[0] = Array<OneD, NekDouble>(m_recvField->GetTotPoints(), 0.0);
        x[1] = Array<OneD, NekDouble>(m_recvField->GetTotPoints(), 0.0);
        x[2] = Array<OneD, NekDouble>(m_recvField->GetTotPoints(), 0.0);
        m_recvField->GetCoords(x[0], x[1], x[2]);

        for (int i = 0; i < m_recvField->GetTotPoints(); ++i)
        {
            if (x[0][i] >= 0.67)
            {
                rVals[5][i] = 0.0;
            }
        }
    }
}

void CouplingCwipi::ExtrapolateFields(Array<OneD, Array<OneD, NekDouble> > &rVals, Array<OneD, int> &notLoc)
{
    Timer timer1, timer2;
    timer1.Start();

    int totvars = 3 + m_nRecvVars;
    int nNotLoc = notLoc.num_elements();

    Array<OneD, Array<OneD, NekDouble> > allVals(totvars);
    Array<OneD, Array<OneD, NekDouble> > scatterVals(totvars);
    Array<OneD, Array<OneD, NekDouble> > gatheredVals(totvars);
    for (int i = 0; i < 3; ++i)
    {
        allVals[i] = Array<OneD, NekDouble>(m_nPoints);
        scatterVals[i] = Array<OneD, NekDouble>(m_nPoints - nNotLoc);
    }
    m_recvField->GetCoords(allVals[0], allVals[1], allVals[2]);

    if (m_spacedim < 3)
    {
        Vmath::Zero(m_nPoints, allVals[2], 1);
    }
    if (m_spacedim < 2)
    {
        Vmath::Zero(m_nPoints, allVals[1], 1);
    }

    for (int i = 0; i < m_nRecvVars; ++i)
    {
        allVals[3 + i] = rVals[i];
        scatterVals[3 + i] = Array<OneD, NekDouble>(m_nPoints - nNotLoc);
    }

    // only copy points from allVals to scatterVals that were located
    for (int i = 0, intPos = 0, locPos = 0; i < m_nPoints; ++i)
    {
        if (locPos < nNotLoc && notLoc[locPos] == i)
        {
            // do nothing
            locPos++;
        }
        else
        {
            for (int k = 0; k < totvars; ++k)
            {
                scatterVals[k][intPos] = allVals[k][i];
            }
            intPos++;
        }
    }

    // send all located points to all ranks. This is probably horribly expensive
    int nranks = m_evalField->GetSession()->GetComm()->GetSize();
    Array<OneD, int> nThisNotLoc(1, m_nPoints - nNotLoc);
    Array<OneD, int> recvDataSizeMap(nranks);
    m_evalField->GetSession()->GetComm()->AllGather(nThisNotLoc, recvDataSizeMap);

    Array<OneD, int> recvDataOffsetMap(nranks);
    recvDataOffsetMap[0] = 0;
    for (int i = 1; i < nranks; ++i)
    {
        recvDataOffsetMap[i] = recvDataOffsetMap[i - 1] + recvDataSizeMap[i - 1];
    }
    int totRecvDataSize = recvDataOffsetMap[nranks - 1] + recvDataSizeMap[nranks - 1];

    timer2.Start();
    for (int i = 0; i < totvars; ++i)
    {
        gatheredVals[i] = Array<OneD, NekDouble>(totRecvDataSize);

        m_evalField->GetSession()->GetComm()->AllGatherv(
            scatterVals[i],
            gatheredVals[i],
            recvDataSizeMap,
            recvDataOffsetMap);
    }
    timer2.Stop();

    if (nNotLoc > 0)
    {
        LibUtilities::PtsFieldSharedPtr gatheredPts =
            MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(3, gatheredVals);

        Array<OneD, Array<OneD, NekDouble > > tmp(totvars);
        for (int j = 0;  j < totvars; ++j)
        {
            tmp[j] = Array<OneD, NekDouble>(nNotLoc);
            for (int i = 0; i < nNotLoc; ++i)
            {
                tmp[j][i] = allVals[j][notLoc[i]];
            }
        }
        LibUtilities::PtsFieldSharedPtr notlocPts =
            MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(3, tmp);

        // perform a nearest neighbour interpolation from gatheredVals to the not located rVals
        if (not m_extrapInterpolator)
        {
            m_extrapInterpolator =
                MemoryManager<FieldUtils::Interpolator>::AllocateSharedPtr(FieldUtils::eNearestNeighbour);
            m_extrapInterpolator->CalcWeights(gatheredPts, notlocPts);
            m_extrapInterpolator->PrintStatistics();
        }
        m_extrapInterpolator->Interpolate(gatheredPts, notlocPts);

        for (int j = 3;  j < totvars; ++j)
        {
            for (int i = 0; i < nNotLoc; ++i)
            {
                allVals[j][notLoc[i]] = tmp[j][i];
            }
        }
    }

    timer1.Stop();
    if (m_evalField->GetComm()->GetRank() == 0 &&
            m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "ExtrapolateFields total time: " << timer1.TimePerTest(1);
        cout << " (AllGatherv: " << timer2.TimePerTest(1) << ")" << endl;
    }
}

void CouplingCwipi::DumpRawFields(const NekDouble time,
                                  Array<OneD, Array<OneD, NekDouble> > rVals)
{
    Timer timer1;
    timer1.Start();

#ifdef _WIN32
    // We need this to make sure boost::format has always
    // two digits in the exponents of Scientific notation.
    unsigned int old_exponent_format;
    old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);
    filename = boost::str(boost::format(filename) % m_time);
    _set_output_format(old_exponent_format);
#else
    std::string filename =
        boost::str(boost::format(m_config["DUMPRAW"]) % time);
#endif

    Array<OneD, Array<OneD, NekDouble> > tmp(m_nRecvVars + m_spacedim);
    for (int i = 0; i < 3; ++i)
    {
        tmp[i] = Array<OneD, NekDouble>(m_recvField->GetTotPoints(), 0.0);
    }
    m_recvField->GetCoords(tmp[0], tmp[1], tmp[2]);

    for (int i = 0; i < m_nRecvVars; ++i)
    {
        tmp[m_spacedim + i] = rVals[i];
    }

    LibUtilities::CsvIO csvIO(m_evalField->GetSession()->GetComm());
    LibUtilities::PtsFieldSharedPtr rvPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
            m_spacedim, m_recvFieldNames, tmp);
    csvIO.Write(filename, rvPts);

    timer1.Stop();
    if (m_evalField->GetComm()->GetRank() == 0 &&
            m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "DumpRawFields total time: " << timer1.TimePerTest(1) << endl;
    }
}
}
}
