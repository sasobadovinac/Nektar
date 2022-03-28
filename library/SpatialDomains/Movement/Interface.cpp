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
#include <SpatialDomains/Movement/Interface.h>
#include <tinyxml.h>

namespace Nektar
{
namespace SpatialDomains
{
Movement::Movement(const LibUtilities::SessionReaderSharedPtr &pSession,
                   const MeshGraphSharedPtr &meshGraph)
    : m_meshGraph(meshGraph), m_session(pSession)
{
    TiXmlElement *xmlDoc     = m_session->GetElement("NEKTAR");
    TiXmlNode *conditionsXml = xmlDoc->FirstChild("CONDITIONS");
    if (conditionsXml != nullptr) // @TODO: Can I remove this line?
    {
        TiXmlNode *movement = conditionsXml->FirstChild("MOVEMENT");
        if (movement != nullptr)
        {
            if (movement->FirstChild("ZONES") != nullptr)
            {
                ReadZones(
                    m_session->GetElement("NEKTAR/CONDITIONS/MOVEMENT/ZONES"));
            }

            if (movement->FirstChild("INTERFACES") != nullptr)
            {
                ReadInterfaces(m_session->GetElement(
                    "NEKTAR/CONDITIONS/MOVEMENT/INTERFACES"));
            }

            // @TODO: Put in a check that both zones and interfaces are defined
        }
    }

    // DEBUG COMMENTS
    if (conditionsXml != nullptr) // Set if verbose/debug mode? to output rank interface information
    {
        if (m_session->GetComm()->GetRank() == 0 &&
                m_session->DefinesCmdLineArgument("verbose"))
        {
            std::cout << "Num zones: " << m_zones.size();
            std::cout << "\n-----------------------------\n";
            std::cout << "Zone ID "
                      << "Type \t"
                      << "# elmts \n";
            std::cout << "-----------------------------" << std::endl;

            std::array<std::string, 4> edgeName = {{"n", "l", "r", "b"}};

            for (auto &zone : m_zones)
            {
                std::cout << zone.first << "\t"
                          << MovementTypeStr[static_cast<int>(
                                 zone.second->GetMovementType())]
                          << "\t" << zone.second->GetElements().size()
                          << std::endl;
            }
            std::cout << "-----------------------------" << std::endl
                      << std::endl;

            std::cout << "Num interfaces: " << m_interfaces.size();
            std::cout << "\n-----------------------------\n";
            std::cout << "Name \t"
                      << "L # elmts \t"
                      << "R # elmts \n";
            std::cout << "-----------------------------" << std::endl;

            for (auto &interface : m_interfaces)
            {
                std::cout
                    << interface.first.second << "\t"
                    << interface.second->GetLeftInterface()->GetEdgeIds().size()
                    << "\t"
                    << interface.second->GetRightInterface()
                           ->GetEdgeIds()
                           .size()
                    << std::endl;
            }
            std::cout << "-----------------------------" << std::endl;
        }
    }
}

ZoneBase::ZoneBase(MovementType type, int indx, CompositeMap domain,
                   int coordDim)
    : m_type(type), m_id(indx), m_domain(domain), m_coordDim(coordDim)
{
    // Fill elements and element IDs from domain
    int shapeDim = coordDim; // @TODO: We are using coordim when this should be shapedim, find easy way in MPI to determine shapedim
    m_constituentElements = Array<OneD,std::set<GeometrySharedPtr>>(shapeDim);

    for (auto &comp : domain)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            m_elementIds.emplace_back(geom->GetGlobalID());
            m_elements.emplace_back(geom);

            // Fill consituent elements (i.e. faces and edges)
            switch (shapeDim)
            {
                case 3:
                    for(int i = 0; i < geom->GetNumFaces(); ++i)
                    {
                        m_constituentElements[2].insert(geom->GetFace(i));
                    }
                    /* fall through */
                case 2:
                    for(int i = 0; i < geom->GetNumEdges(); ++i)
                    {
                        m_constituentElements[1].insert(geom->GetEdge(i));
                    }
                    /* fall through */
                case 1:
                    m_constituentElements[0].insert(geom);
                    /* fall through */
                default:
                    break;
            }
        }
    }
}

ZoneRotate::ZoneRotate(int id, const CompositeMap &domain, const int coordDim,
                       const NekPoint<NekDouble> &origin, const DNekVec &axis,
                       const LibUtilities::EquationSharedPtr &angularVelEqn)
    : ZoneBase(MovementType::eRotate, id, domain, coordDim), m_origin(origin),
      m_axis(axis), m_angularVelEqn(angularVelEqn)
{
    std::set<int> seenVerts, seenEdges, seenFaces;
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

            for (int i = 0; i < geom->GetNumFaces(); ++i)
            {
                Geometry2DSharedPtr face =
                    std::static_pointer_cast<Geometry2D>(geom->GetFace(i));

                if (seenFaces.find(face->GetGlobalID()) != seenFaces.end())
                {
                    continue;
                }

                seenFaces.insert(face->GetGlobalID());

                CurveSharedPtr curve = face->GetCurve();
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

    // Construct rotation matrix
    m_W(0, 1) = -m_axis[2];
    m_W(0, 2) = m_axis[1];
    m_W(1, 0) = m_axis[2];
    m_W(1, 2) = -m_axis[0];
    m_W(2, 0) = -m_axis[1];
    m_W(2, 1) = m_axis[0];

    m_W2 = m_W * m_W;

    // std::cout << "W matrix: " << std::endl;
    // std::cout << m_W << std::endl << std::endl;
    // std::cout << "W^2 matrix: " << std::endl;
    // std::cout << m_W2 << std::endl << std::endl;
    // std::cout << "Origin: " << std::endl;
    // std::cout << m_origin << std::endl << std::endl;
}

ZoneTranslate::ZoneTranslate(int id, const CompositeMap &domain,
                             const int coordDim,
                             const std::vector<NekDouble> &velocity)
    : ZoneBase(MovementType::eTranslate, id, domain, coordDim),
      m_velocity(velocity)
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

ZonePrescribe::ZonePrescribe(int id, const CompositeMap &domain,
                             const int coordDim,
                             LibUtilities::EquationSharedPtr xDeform,
                             LibUtilities::EquationSharedPtr yDeform,
                             LibUtilities::EquationSharedPtr zDeform)
    : ZoneBase(MovementType::ePrescribe, id, domain, coordDim),
      m_xDeform(xDeform), m_yDeform(yDeform), m_zDeform(zDeform)
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
    int coordDim = m_meshGraph->GetSpaceDimension();

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
        err = zonesElement->QueryStringAttribute("DOMAIN", &interfaceDomainStr);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read zone domain.");

        auto &domains = m_meshGraph->GetDomain();
        auto domFind  = stoi(ReadTag(interfaceDomainStr));
        std::map<int, CompositeSharedPtr> domain;
        if (domains.find(domFind) != domains.end())
        {
            domain = domains.at(domFind);
        }

        ZoneBaseShPtr zone;
        if (zoneType == "R" || zoneType == "ROTATE" || zoneType == "ROTATING")
        {
            std::string originStr;
            err = zonesElement->QueryStringAttribute("ORIGIN", &originStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read origin.");
            std::vector<NekDouble> originVec;
            ParseUtils::GenerateVector(originStr, originVec);
            NekPoint<NekDouble> origin =
                NekPoint<NekDouble>(originVec[0], originVec[1], originVec[2]);

            std::string axisStr;
            err = zonesElement->QueryStringAttribute("AXIS", &axisStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read axis.");
            DNekVec axis(axisStr);
            axis.Normalize();

            std::string angularVelStr;
            err = zonesElement->QueryStringAttribute("ANGVEL", &angularVelStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read angular velocity.");

            m_session->SubstituteExpressions(angularVelStr);
            LibUtilities::EquationSharedPtr angularVelEqn =
                MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(
                        m_session->GetInterpreter(), angularVelStr);

            zone = ZoneRotateShPtr(MemoryManager<ZoneRotate>::AllocateSharedPtr(
                indx, domain, coordDim, origin, axis, angularVelEqn));

            m_moveFlag = true;
        }
        else if (zoneType == "T" || zoneType == "TRANSLATE" ||
                 zoneType == "TRANSLATING")
        {
            std::string velocityStr;
            err = zonesElement->QueryStringAttribute("VELOCITY", &velocityStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read direction.");
            std::vector<NekDouble> velocity;
            ParseUtils::GenerateVector(velocityStr, velocity);

            zone = ZoneTranslateShPtr(
                MemoryManager<ZoneTranslate>::AllocateSharedPtr(
                    indx, domain, coordDim, velocity));

            m_moveFlag = true;
        }
        else if (zoneType == "F" || zoneType == "FIXED")
        {
            zone = ZoneFixedShPtr(MemoryManager<ZoneFixed>::AllocateSharedPtr(
                indx, domain, coordDim));
        }
        else if (zoneType == "P" || zoneType == "PRESCRIBED")
        {
            std::string xDeformStr;
            err = zonesElement->QueryStringAttribute("XDEFORM", &xDeformStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read x deform equation.");
            LibUtilities::EquationSharedPtr xDeformEqn =
                std::make_shared<LibUtilities::Equation>(
                    m_session->GetInterpreter(), xDeformStr);

            std::string yDeformStr;
            err = zonesElement->QueryStringAttribute("YDEFORM", &yDeformStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read y deform equation.");
            LibUtilities::EquationSharedPtr yDeformEqn =
                std::make_shared<LibUtilities::Equation>(
                    m_session->GetInterpreter(), yDeformStr);

            std::string zDeformStr;
            err = zonesElement->QueryStringAttribute("ZDEFORM", &zDeformStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read z deform equation.");
            LibUtilities::EquationSharedPtr zDeformEqn =
                std::make_shared<LibUtilities::Equation>(
                    m_session->GetInterpreter(), zDeformStr);

            zone = ZonePrescribeShPtr(
                MemoryManager<ZonePrescribe>::AllocateSharedPtr(
                    indx, domain, coordDim, xDeformEqn, yDeformEqn, zDeformEqn));

            m_moveFlag = true;
        }
        else
        {
            WARNINGL0(false, "Zone type '" + zoneType +
                                 "' is unsupported. Valid types are: 'Fixed', 'Rotate', 'Translate', or 'Prescribe'.")
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
        ASSERTL0(
            "INTERFACE" == (std::string)interfaceElement->Value(),
            "Only INTERFACE tags may be present inside the INTERFACES block.")

        int err;

        std::string name;
        err = interfaceElement->QueryStringAttribute("NAME", &name);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read interface name.");
        TiXmlElement *sideElement = interfaceElement->FirstChildElement();

        int cnt = 0;
        Array<OneD, InterfaceShPtr> interfaces(2);
        while (sideElement)
        {
            ASSERTL0(cnt < 2,
                     "Only two sides may be present in each interface block.")
            std::string sideStr = sideElement->Value();

            InterfaceSide side = InterfaceSide::eNone; // @TODO: Currently don't use these sides. Change to make sure we define a left and right.
            if (sideStr == "L" || sideStr == "LEFT")
            {
                side = InterfaceSide::eLeft;
            }
            else if (sideStr == "R" || sideStr == "RIGHT")
            {
                side = InterfaceSide::eRight;
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                         "Only LEFT or RIGHT tags may be present inside the INTERFACE block.")
            }

            int indx; // @TODO: Do I need ID?
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

            interfaces[cnt++] =
                InterfaceShPtr(MemoryManager<Interface>::AllocateSharedPtr(
                    indx, side, boundaryEdge));

            sideElement = sideElement->NextSiblingElement();
        }

        m_interfaces[std::make_pair(m_interfaces.size(), name)] =
            InterfacePairShPtr(MemoryManager<InterfacePair>::AllocateSharedPtr(
                interfaces[0], interfaces[1]));
        interfaceElement = interfaceElement->NextSiblingElement();
    }
}

Interface::Interface(int indx, InterfaceSide side, CompositeMap edge)
    : m_id(indx), m_side(side)
{
    // Fill element Ids
    for (auto &comp : edge)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            m_edgeIds.emplace_back(geom->GetGlobalID());
            m_edge[geom->GetGlobalID()] = geom;
            m_edgeDeque.emplace_back(geom);
        }
    }
}

void Movement::PerformMovement(NekDouble timeStep)
{
    std::set<int> movedZoneIds;
    for (auto &zone : m_zones)
    {
        if (zone.second->Move(timeStep))
        {
            movedZoneIds.insert(zone.first);
        }
    }

    // If zone has moved, set all interfaces on that zone to moved.
    // @TODO: Probably better to save the moved flag on the interface pair obj?
    for (auto &interPair : m_interfaces)
    {
        int leftId  = interPair.second->GetLeftInterface()->GetId();
        int rightId = interPair.second->GetRightInterface()->GetId();

        if (movedZoneIds.find(leftId) != movedZoneIds.end() ||
            movedZoneIds.find(rightId) != movedZoneIds.end())
        {
            m_zones[leftId]->GetMoved()  = true;
            m_zones[rightId]->GetMoved() = true;
        }
    }
}

NekDouble ZoneRotate::GetAngularVel(NekDouble &time) const
{
    int spinUpTime = 10;
    if(time < spinUpTime)
    {
        return m_angularVelEqn->Evaluate(0,0,0,time) * spinUpTime / 10;
    }
    else
    {
        return m_angularVelEqn->Evaluate(0,0,0,time);
    }
}

// Calculate new location of points using Rodrigues formula
bool ZoneRotate::v_Move(NekDouble time)
{
    NekDouble angle = GetAngularVel(time) * time;
    // TODO: For none constant angular velocity this doesn't work ^^
    // @TODO: I need to take into account the total angle, summing timesteps works here
    // @TODO: but then it doesn't work for FieldConvert where only the checkpoint time is known
    // TODO: I want to integrate m_angularVelEqn up to current time
    //NekDouble period = 20.0;
    //NekDouble angle = ceil(1.0/period*time - 1.0) * time -  period / 2.0 * (ceil(1.0/period*time - 1.0) *  ceil(1.0/period*time));

    // Identity matrix
    DNekMat rot(3, 3, 0.0);
    rot(0, 0) = 1.0;
    rot(1, 1) = 1.0;
    rot(2, 2) = 1.0;

    rot = rot + sin(angle) * m_W + (1 - cos(angle)) * m_W2;

    //std::cout << "Rotation matrix: " << std::endl;
    //std::cout << rot << std::endl;
    int cnt = 0;
    for (auto &vert : m_rotateVerts)
    {
        NekPoint<NekDouble> pnt = m_origPosition[cnt] - m_origin;
        DNekVec pntVec          = {pnt[0], pnt[1], pnt[2]};

        DNekVec newLoc = rot * pntVec;

        vert->UpdatePosition(newLoc(0) + m_origin[0], newLoc(1) + m_origin[1],
                             newLoc(2) + m_origin[2]);
        cnt++;
    }

    for (auto &curve : m_rotateCurves)
    {
        for (auto &vert : curve->m_points)
        {
            NekPoint<NekDouble> pnt = m_origPosition[cnt] - m_origin;
            DNekVec pntVec          = {pnt[0], pnt[1], pnt[2]};

            DNekVec newLoc = rot * pntVec;

            vert->UpdatePosition(newLoc(0) + m_origin[0],
                                 newLoc(1) + m_origin[1],
                                 newLoc(2) + m_origin[2]);
            cnt++;
        }
    }

    // Clear bounding boxes (these will be regenerated next time GetBoundingBox is called)
    for (auto &el : m_elements)
    {
        el->ClearBoundingBox();

        int nfaces = el->GetNumFaces();
        for (int i = 0; i < nfaces; ++i)
        {
            el->GetFace(i)->ClearBoundingBox();
        }

        int nedges = el->GetNumEdges();
        for (int i = 0; i < nedges; ++i)
        {
            el->GetEdge(i)->ClearBoundingBox();
        }
    }

    return true;
}

bool ZoneTranslate::v_Move(NekDouble timeStep)
{
    Array<OneD, NekDouble> dist(3, 0.0);
    for (int i = 0; i < m_coordDim; ++i)
    {
        dist[i] = m_velocity[i] * timeStep;
    }

    int cnt = 0;
    for (auto &vert : m_slideVerts)
    {
        Array<OneD, NekDouble> newLoc(3, 0.0);
        auto pnt = m_origPosition[cnt];

        for (int i = 0; i < m_coordDim; ++i)
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

            for (int i = 0; i < m_coordDim; ++i)
            {
                newLoc[i] = pnt(i) + dist[i];
            }

            vert->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
            cnt++;
        }
    }

    // Clear bounding boxes (these will be regenerated next time GetBoundingBox is called)
    for (auto &el : m_elements)
    {
        el->ClearBoundingBox();

        int nfaces = el->GetNumFaces();
        for (int i = 0; i < nfaces; ++i)
        {
            el->GetFace(i)->ClearBoundingBox();
        }

        int nedges = el->GetNumEdges();
        for (int i = 0; i < nedges; ++i)
        {
            el->GetEdge(i)->ClearBoundingBox();
        }
    }

    return true;
}

bool ZoneFixed::v_Move(NekDouble time)
{
    boost::ignore_unused(time);
    return false;
}

bool ZonePrescribe::v_Move(NekDouble time)
{
    int cnt = 0;
    for (auto &vert : m_interiorVerts)
    {
        auto pnt = m_origPosition[cnt++];

        Array<OneD, NekDouble> coords(3, 0.0);
        vert->GetCoords(coords);

        Array<OneD, NekDouble> newLoc(3, 0.0);
        newLoc[0] = m_xDeform->Evaluate(coords[0], coords[1], coords[2], time) + pnt(0);
        newLoc[1] = m_yDeform->Evaluate(coords[0], coords[1], coords[2], time) + pnt(1);
        newLoc[2] = m_zDeform->Evaluate(coords[0], coords[1], coords[2], time) + pnt(2);

        vert->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
    }

    // Clear bounding boxes (these will be regenerated next time GetBoundingBox is called)
    for (auto &el : m_elements)
    {
        el->ClearBoundingBox();

        int nfaces = el->GetNumFaces();
        for (int i = 0; i < nfaces; ++i)
        {
            el->GetFace(i)->ClearBoundingBox();
        }

        int nedges = el->GetNumEdges();
        for (int i = 0; i < nedges; ++i)
        {
            el->GetEdge(i)->ClearBoundingBox();
        }
    }

    return true;
}

}
}