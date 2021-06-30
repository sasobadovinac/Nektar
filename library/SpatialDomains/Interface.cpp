////////////////////////////////////////////////////////////////////////////////
//
//  File: InterfaceConditions.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/Interface.h>
#include <tinyxml.h>

using namespace std;

namespace Nektar
{
namespace SpatialDomains
{
Movement::Movement(const LibUtilities::SessionReaderSharedPtr &pSession,
                       const MeshGraphSharedPtr &meshGraph)
    : m_meshGraph(meshGraph), m_session(pSession)
{
    TiXmlElement *xmlDoc = m_session->GetElement("NEKTAR");
    if(xmlDoc->FirstChild("CONDITIONS") != nullptr) // @TODO: Can I remove this line?
    {
        if(xmlDoc->FirstChild("CONDITIONS")->FirstChild("MOVEMENT")->FirstChild("ZONES") != nullptr)
        {
            ReadZones(m_session->GetElement("NEKTAR/CONDITIONS/MOVEMENT/ZONES"));
        }
        if(xmlDoc->FirstChild("CONDITIONS")->FirstChild("MOVEMENT")->FirstChild("INTERFACES") != nullptr)
        {
            ReadInterfaces(m_session->GetElement("NEKTAR/CONDITIONS/MOVEMENT/INTERFACES"));
        }
    }
}

ZoneBase::ZoneBase(MovementType type, int indx, CompositeMap domain)
: m_type(type), m_id(indx), m_domain(domain)
{
    // Fill element Ids
    for (auto &comp : domain)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            m_elementIds.emplace_back(geom->GetGlobalID());
            m_elements.emplace_back(geom);
        }
    }
}

ZoneRotate::ZoneRotate(int id, const CompositeMap &domain,
                                     const PointGeom &origin,
                                     const std::vector<NekDouble> &axis,
                                     const NekDouble angularVel)
    : ZoneBase(eRotating, id, domain), m_origin(origin), m_axis(axis),
      m_angularVel(angularVel)
{
    std::set<int> seenVerts, seenEdges;
    for (auto &comp : m_domain)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            for (int i = 0; i < geom->GetNumVerts(); ++i)
            {
                PointGeomSharedPtr vert = geom->GetVertex(i);

                if (seenVerts.find(vert->GetGlobalID()) != seenVerts.end())
                {
                    continue;
                }

                seenVerts.insert(vert->GetGlobalID());
                m_rotateVerts.emplace_back(vert);
            }

            for (int i = 0; i < geom->GetNumEdges(); ++i)
            {
                SegGeomSharedPtr edge =
                    std::static_pointer_cast<SegGeom>(geom->GetEdge(i));

                if (seenEdges.find(edge->GetGlobalID()) != seenEdges.end())
                {
                    continue;
                }

                seenEdges.insert(edge->GetGlobalID());

                CurveSharedPtr curve = edge->GetCurve();
                if (!curve)
                {
                    continue;
                }

                m_rotateCurves.emplace_back(curve);
            }
        }
    }

    // Copy points so we know original positions.
    for (auto &pt : m_rotateVerts)
    {
        m_origPosition.emplace_back(*pt);
    }

    for (auto &curve : m_rotateCurves)
    {
        for (auto &pt : curve->m_points)
        {
            m_origPosition.emplace_back(*pt);
        }
    }
}

ZoneTranslate::ZoneTranslate(int id, const CompositeMap &domain,
                                   const std::vector<NekDouble> &velocity)
    : ZoneBase(eSliding, id, domain), m_velocity(velocity)
{
    std::set<int> seenVerts, seenEdges;
    for (auto &comp : m_domain)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            for (int i = 0; i < geom->GetNumVerts(); ++i)
            {
                PointGeomSharedPtr vert = geom->GetVertex(i);

                if (seenVerts.find(vert->GetGlobalID()) != seenVerts.end())
                {
                    continue;
                }

                seenVerts.insert(vert->GetGlobalID());
                m_slideVerts.emplace_back(vert);
            }

            for (int i = 0; i < geom->GetNumEdges(); ++i)
            {
                SegGeomSharedPtr edge =
                    std::static_pointer_cast<SegGeom>(geom->GetEdge(i));

                if (seenEdges.find(edge->GetGlobalID()) != seenEdges.end())
                {
                    continue;
                }

                seenEdges.insert(edge->GetGlobalID());

                CurveSharedPtr curve = edge->GetCurve();
                if (!curve)
                {
                    continue;
                }

                m_slideCurves.emplace_back(curve);
            }
        }
    }

    // Copy points so we know original positions.
    for (auto &pt : m_slideVerts)
    {
        m_origPosition.emplace_back(*pt);
    }

    for (auto &curve : m_slideCurves)
    {
        for (auto &pt : curve->m_points)
        {
            m_origPosition.emplace_back(*pt);
        }
    }
}

ZonePrescribe::ZonePrescribe(int id,
                                         const CompositeMap &domain,
                                         LibUtilities::EquationSharedPtr xDeform,
                                         LibUtilities::EquationSharedPtr yDeform)
    : ZoneBase(ePrescribed, id, domain), m_xDeform(xDeform), m_yDeform(yDeform)
{
    std::set<int> seenVerts, seenEdges;
    for (auto &comp : m_domain)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            for (int i = 0; i < geom->GetNumVerts(); ++i)
            {
                PointGeomSharedPtr vert = geom->GetVertex(i);

                if (seenVerts.find(vert->GetGlobalID()) != seenVerts.end())
                {
                    continue;
                }

                seenVerts.insert(vert->GetGlobalID());
                m_interiorVerts.emplace_back(vert);
            }

            for (int i = 0; i < geom->GetNumEdges(); ++i)
            {
                SegGeomSharedPtr edge =
                    std::static_pointer_cast<SegGeom>(geom->GetEdge(i));

                if (seenEdges.find(edge->GetGlobalID()) != seenEdges.end())
                {
                    continue;
                }

                seenEdges.insert(edge->GetGlobalID());

                CurveSharedPtr curve = edge->GetCurve();
                if (!curve)
                {
                    continue;
                }

                for (auto &point : curve->m_points)
                {
                    if (seenVerts.find(point->GetGlobalID()) != seenVerts.end())
                    {
                        m_interiorVerts.emplace_back(point);
                    }
                }
            }
        }
    }

    // Copy points so we know original positions.
    for (auto &pt : m_interiorVerts)
    {
        m_origPosition.emplace_back(*pt);
    }
}



std::string ReadTag(std::string &tagStr)
{
    std::string::size_type indxBeg = tagStr.find_first_of('[') + 1;
    std::string::size_type indxEnd = tagStr.find_last_of(']') - 1;

    ASSERTL0(
        indxBeg <= indxEnd,
        (std::string("Error reading interface region definition:") + tagStr)
            .c_str());

    std::string indxStr = tagStr.substr(indxBeg, indxEnd - indxBeg + 1);

    return indxStr;
}

void Movement::ReadZones(TiXmlElement *zonesTag)
{
    ASSERTL0(zonesTag, "Unable to find ZONES tag in file.");

    TiXmlElement *zonesElement = zonesTag->FirstChildElement();

    while (zonesElement)
    {
        std::string zoneType = zonesElement->Value();

        int err;
        int indx;

        err = zonesElement->QueryIntAttribute("ID", &indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read zone ID.");

        std::string interfaceDomainStr;
        err = zonesElement->QueryStringAttribute("DOMAIN",
                                                     &interfaceDomainStr);
        ASSERTL0(err == TIXML_SUCCESS,
                 "Unable to read zone domain.");

        auto &domains = m_meshGraph->GetDomain();
        auto domFind = stoi(ReadTag(interfaceDomainStr));
        std::map<int, CompositeSharedPtr> domain;
        if (domains.find(domFind) != domains.end())
        {
            domain = domains.at(domFind);
        }

        /*map<int, CompositeSharedPtr> domainEdge;
        if (zoneType != "P")
        {
            std::string domainEdgeStr;
            int domainEdgeErr =
                zonesElement->QueryStringAttribute("EDGE", &domainEdgeStr);

            if (domainEdgeErr == TIXML_SUCCESS)
            {
                std::string indxStr = ReadTag(domainEdgeStr);
                m_meshGraph->GetCompositeList(indxStr, domainEdge);
            }
        }*/

        ZoneBaseShPtr zone;
        if (zoneType == "R" || zoneType == "ROTATE" || zoneType == "ROTATING")
        {
            std::string originStr;
            err = zonesElement->QueryStringAttribute("ORIGIN", &originStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read origin.");
            std::vector<NekDouble> originVec;
            ParseUtils::GenerateVector(originStr, originVec);
            auto origin =
                PointGeom(3, 0, originVec[0], originVec[1], originVec[2]);

            std::string axisStr;
            err = zonesElement->QueryStringAttribute("AXIS", &axisStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read axis.");
            std::vector<NekDouble> axis;
            ParseUtils::GenerateVector(axisStr, axis);

            std::string angularVelStr;
            err = zonesElement->QueryStringAttribute("ANGVEL", &angularVelStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read angular velocity.");

            LibUtilities::Equation angularVelEqn(
                m_session->GetInterpreter(), angularVelStr);
            NekDouble angularVel = angularVelEqn.Evaluate();

            zone = ZoneRotateShPtr(MemoryManager<ZoneRotate>::AllocateSharedPtr(indx, domain,  origin, axis, angularVel));
        }
        else if (zoneType == "T" || zoneType == "TRANSLATE" || zoneType == "TRANSLATING")
        {
            std::string velocityStr;
            err = zonesElement->QueryStringAttribute("VELOCITY", &velocityStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read direction.");
            std::vector<NekDouble> velocity;
            ParseUtils::GenerateVector(velocityStr, velocity);

            zone = ZoneTranslateShPtr(MemoryManager<ZoneTranslate>::AllocateSharedPtr(indx, domain, velocity));

        }
        else if (zoneType == "F" || zoneType == "FIXED")
        {
            zone = ZoneFixedShPtr(MemoryManager<ZoneFixed>::AllocateSharedPtr(indx, domain));
        }
        else if (zoneType == "P" || zoneType == "PRESCRIBED")
        {
            std::string xDeformStr;
            err = zonesElement->QueryStringAttribute("XDEFORM", &xDeformStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read x deform equation.");
            LibUtilities::EquationSharedPtr xDeformEqn =
                std::make_shared<LibUtilities::Equation>(m_session->GetInterpreter(), xDeformStr);

            std::string yDeformStr;
            err = zonesElement->QueryStringAttribute("YDEFORM", &yDeformStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read y deform equation.");
            LibUtilities::EquationSharedPtr yDeformEqn =
                std::make_shared<LibUtilities::Equation>(m_session->GetInterpreter(), yDeformStr);

            zone = ZonePrescribeShPtr(MemoryManager<ZonePrescribe>::AllocateSharedPtr(indx, domain, xDeformEqn, yDeformEqn));
        }

        m_zones[indx] = zone;

        zonesElement = zonesElement->NextSiblingElement();
    }
}

void Movement::ReadInterfaces(TiXmlElement *interfacesTag)
{
    ASSERTL0(interfacesTag, "Unable to find INTERFACES tag in file.");

    TiXmlElement *interfaceElement = interfacesTag->FirstChildElement();

    while (interfaceElement)
    {
        ASSERTL0("INTERFACE" == (std::string)interfaceElement->Value(),
                 "Only INTERFACE tags may be present inside the INTERFACES block.")

        int err;

        std::string name;
        err = interfaceElement->QueryStringAttribute("NAME", &name);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read interface name.");

        TiXmlElement *sideElement = interfaceElement->FirstChildElement();

        int cnt = 0;
        Array<OneD, InterfaceShPtr> interfaces(2);
        while(sideElement)
        {
            ASSERTL0(cnt < 2, "Only two sides may be present in each interface block.")
            std::string sideStr = interfaceElement->Value();

            InterfaceSide side;
            if(sideStr == "L" || sideStr == "LEFT")
            {
                side = eLeft;
            }
            else if (sideStr == "R" || sideStr == "RIGHT")
            {
                side = eRight;
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                         "Only LEFT or RIGHT tags may be present inside the INTERFACE block.")
            }

            int indx;
            err = sideElement->QueryIntAttribute("ID", &indx);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read interface ID.");

            std::string boundaryStr;
            int boundaryErr =
                sideElement->QueryStringAttribute("BOUNDARY", &boundaryStr);

            CompositeMap boundaryEdge;
            if (boundaryErr == TIXML_SUCCESS)
            {
                std::string indxStr = ReadTag(boundaryStr);
                m_meshGraph->GetCompositeList(indxStr, boundaryEdge);
            }

            interfaces[cnt] = InterfaceShPtr(MemoryManager<Interface>::AllocateSharedPtr(indx, boundaryEdge));

            sideElement = sideElement->NextSiblingElement();
        }

        m_interfaces[name] = InterfacePairShPtr(MemoryManager<InterfacePair>::AllocateSharedPtr(interfaces[0], interfaces[1]));
        interfaceElement = interfaceElement->NextSiblingElement();
    }
}

Interface::Interface(int indx, CompositeMap edge)
    : m_id(indx)
{
    // Fill element Ids
    for (auto &comp : edge)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            m_edgeIds.emplace_back(geom->GetGlobalID());
            m_edge[geom->GetGlobalID()] = geom;
        }
    }
}

void ZoneBase::SetEdge(const CompositeMap &edge)
{
    for (auto &compIt : edge)
    {
        for (auto &elmtIt : compIt.second->m_geomVec)
        {
            auto shapeType = elmtIt->GetShapeType();
            auto dim = elmtIt->GetCoordim();

            ASSERTL0((shapeType == LibUtilities::eSegment && dim == 2) ||
                     (shapeType == LibUtilities::eQuadrilateral && dim == 3) ||
                     (shapeType == LibUtilities::eTriangle && dim == 3),
                "Interface edge must be a segment for 2D or quad/tri for 3D.")

            size_t id = elmtIt->GetGlobalID();
            m_edge[id] = elmtIt;
            m_edgeIds.emplace_back(id);
        }
    }

    std::sort(m_edgeIds.begin(), m_edgeIds.end());
}

void Movement::PerformMovement(NekDouble timeStep)
{
    for (auto &interface : m_zones)
    {
        interface->Move(timeStep);
    }
}

void Movement::GenGeomFactors()
{
    for (auto &interface : m_zones)
    {
        auto edges = interface->GetEdge();
        for (auto &el : edges)
        {
            el.second->GenGeomFactors();
        }
    }

    // @TODO: Don't know if below is needed?
    /*for (auto &interface : m_zones)
    {
        auto elements = interface->GetElements();
        for (auto &el : elements)
        {
            int ne = el->GetNumEdges();
            for (int i = 0; i < ne; ++i)
            {
                el->GetEdge(i)->GenGeomFactors();
            }

            el->GenGeomFactors();
        }
    }*/
}

void ZoneRotate::v_Move(NekDouble time)
{
    NekDouble angle = m_angularVel * time;

    NekDouble ox, oy, oz;
    m_origin.GetCoords(ox, oy, oz);
    std::cout << "ox = " << ox << " oy = " << oy << " oz = " << oz << std::endl;
    std::cout << "angular vel = " << m_angularVel << " time = " << time << " angle = " << angle << std::endl;

    int cnt = 0;
    for (auto &vert : m_rotateVerts)
    {
        Array<OneD, NekDouble> newLoc(3, 0.0);
        auto pnt = m_origPosition[cnt];
        newLoc[0] = std::cos(angle) * pnt(0) - std::sin(angle) * pnt(1);
        newLoc[1] = std::sin(angle) * pnt(0) + std::cos(angle) * pnt(1);

        vert->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
        cnt++;
    }

    for (auto &curve : m_rotateCurves)
    {
        for (auto &vert : curve->m_points)
        {
            Array<OneD, NekDouble> newLoc(3, 0.0);
            auto pnt = m_origPosition[cnt];
            newLoc[0] = std::cos(angle) * pnt(0) - std::sin(angle) * pnt(1);
            newLoc[1] = std::sin(angle) * pnt(0) + std::cos(angle) * pnt(1);

            vert->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
            cnt++;
        }
    }

    m_moved = true;
    m_oppInterface->GetMoved() = true;
}

void ZoneTranslate::v_Move(NekDouble timeStep)
{
    int dim = 3; // @TODO: Think a way to get this for coorddim even if interface isn't present on rank.

    Array<OneD, NekDouble> dist(3, 0.0);
    for (int i = 0; i < dim; ++i)
    {
        dist[i] = m_velocity[i] * timeStep;
    }

    int cnt = 0;
    for (auto &vert : m_slideVerts)
    {
        Array<OneD, NekDouble> newLoc(3, 0.0);
        auto pnt = m_origPosition[cnt];

        for (int i = 0; i < dim; ++i)
        {
            newLoc[i] = pnt(i) + dist[i];
        }

        vert->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
        cnt++;
    }

    for (auto &curve : m_slideCurves)
    {
        for (auto &vert : curve->m_points)
        {
            Array<OneD, NekDouble> newLoc(3, 0.0);
            auto pnt = m_origPosition[cnt];

            for (int i = 0; i < dim; ++i)
            {
                newLoc[i] = pnt(i) + dist[i];
            }

            vert->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
            cnt++;
        }
    }

    m_moved = true;
    m_oppInterface->GetMoved() = true;
}

void ZoneFixed::v_Move(NekDouble time)
{
    boost::ignore_unused(time);
}

void ZonePrescribe::v_Move(NekDouble time)
{
    boost::ignore_unused(time);
    // This is hacky - as interface is set up for 2 sides usually, we only use the left side in this case
    if (m_side == eLeft)
    {
        int dim = 3;
        int cnt = 0;

        for (auto &vert : m_interiorVerts)
        {
            Array<OneD, NekDouble> coords(dim, 0.0);
            vert->GetCoords(coords);

            Array<OneD, NekDouble> newLoc(3, 0.0);

            // newLoc[0] = m_xDeform->Evaluate(coords[0], coords[1], coords[2], time);
            // newLoc[1] = m_yDeform->Evaluate(coords[0], coords[1], coords[2], time); newLoc[2] = coords[2];

            /*NekDouble Lx = 20, Ly = 20;         // Size of mesh
            NekDouble nx = 1, ny = 1, nt = 1;   // Space and time period
            NekDouble X0 = 0.5, Y0 = 0.5;       // Amplitude
            NekDouble t0 = sqrt(5*5 + 5*5);     // Time domain*/

            //newLoc[0] = coords[0] + X0 * sin((nt * 2 * M_PI * time) / t0)
                                       //* sin((nx * 2 * M_PI * coords[0]) / Lx)
                                       //* sin((ny * 2 * M_PI * coords[1]) / Ly);

            //newLoc[1] = coords[1] + Y0 * sin((nt * 2 * M_PI * time) / t0)
                                       //* sin((nx * 2 * M_PI * coords[0]) / Lx)
                                       //* sin((ny * 2 * M_PI * coords[1]) / Ly);
            /*
            if (coords[0] < 1e-8 || fabs(coords[0] - 1) < 1e-8)
            {
                newLoc[0] = coords[0]; // + 0.001 * sin(2 * M_PI * time) * coords[0] * (1 - coords[0]);
                newLoc[1] = coords[1];
                newLoc[2] = coords[2];
            }
            else
            {
                newLoc[0] = coords[0] + 0.001 * 0.5; // + 0.001 * sin(2 * M_PI * time) * coords[0] * (1 - coords[0]);
                newLoc[1] = coords[1];
                newLoc[2] = coords[2];
            }
            */

            auto pnt = m_origPosition[cnt];
            newLoc[0] = pnt(0) + 0.05 * sin(2*M_PI*time) * sin(2*M_PI*coords[0]) * sin(2*M_PI*coords[1]);
            newLoc[1] = pnt(1) + 0.05 * sin(2*M_PI*time) * sin(2*M_PI*coords[0]) * sin(2*M_PI*coords[1]);
            cnt++;

            vert->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
        }
    }

    for (auto &el : m_elements)
    {
        SpatialDomains::CurveMap edges;
        el->Reset(edges, edges);
        el->Setup();
        el->GenGeomFactors();
        el->FillGeom();
    }
}

}
}
