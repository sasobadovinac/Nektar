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
Interfaces::Interfaces(const LibUtilities::SessionReaderSharedPtr &pSession,
                       const MeshGraphSharedPtr &meshGraph)
    : m_meshGraph(meshGraph), m_session(pSession)
{
    TiXmlElement *xmlDoc = m_session->GetElement("NEKTAR");
    if(xmlDoc->FirstChild("CONDITIONS") != nullptr)
    {
        if(xmlDoc->FirstChild("CONDITIONS")->FirstChild("INTERFACES") != nullptr)
        {
            ReadInterfaces(
                m_session->GetElement("NEKTAR/CONDITIONS/INTERFACES"));
        }
    }
}

InterfaceBase::InterfaceBase(InterfaceType type, int indx, CompositeMap domain)
: m_type(type), m_id(indx), m_domain(domain)
{
    // Fill element Ids
    for (auto &comp : domain)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            m_elementIds.emplace_back(geom->GetGlobalID());
        }
    }
}

RotatingInterface::RotatingInterface(int id, const CompositeMap &domain,
                                     const PointGeom &origin,
                                     const std::vector<NekDouble> &axis,
                                     const NekDouble angularVel)
    : InterfaceBase(eRotating, id, domain), m_origin(origin), m_axis(axis),
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

void Interfaces::ReadInterfaces(TiXmlElement *interfaces)
{

    ASSERTL0(interfaces, "Unable to find INTERFACES tag in file.");

    TiXmlElement *interfaceElement = interfaces->FirstChildElement();

    std::map<int, std::vector<InterfaceBaseShPtr>> tmpInterfaceMap;

    while (interfaceElement)
    {
        std::string interfaceType = interfaceElement->Value();

        int err;
        int indx;

        err = interfaceElement->QueryIntAttribute("ID", &indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read interface ID.");

        std::string interfaceDomainStr;
        err = interfaceElement->QueryStringAttribute("DOMAIN",
                                                     &interfaceDomainStr);
        ASSERTL0(err == TIXML_SUCCESS,
                 "Unable to read interface domain.");
        auto domain =
            m_meshGraph->GetDomain(stoi(ReadTag(interfaceDomainStr)));

        std::string domainEdgeStr;
        int domainEdgeErr =
            interfaceElement->QueryStringAttribute("EDGE", &domainEdgeStr);
        map<int, CompositeSharedPtr> domainEdge;
        if (domainEdgeErr == TIXML_SUCCESS)
        {
            std::string indxStr = ReadTag(domainEdgeStr);
            m_meshGraph->GetCompositeList(indxStr, domainEdge);
        }

        InterfaceBaseShPtr interface;

        if (interfaceType == "R")
        {
            std::string originStr;
            err = interfaceElement->QueryStringAttribute("ORIGIN", &originStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read origin.");
            std::vector<NekDouble> originVec;
            ParseUtils::GenerateVector(originStr, originVec);
            auto origin =
                PointGeom(3, 0, originVec[0], originVec[1], originVec[2]);

            std::string axisStr;
            err = interfaceElement->QueryStringAttribute("AXIS", &axisStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read axis.");
            std::vector<NekDouble> axis;
            ParseUtils::GenerateVector(axisStr, axis);

            std::string angularVelStr;
            err = interfaceElement->QueryStringAttribute("ANGVEL", &angularVelStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read angular velocity.");
            LibUtilities::Equation angularVelEqn(
                m_session->GetInterpreter(), angularVelStr);
            NekDouble angularVel = angularVelEqn.Evaluate();

            interface = RotatingInterfaceShPtr(MemoryManager<RotatingInterface>::AllocateSharedPtr(indx, domain,  origin, axis, angularVel));
        }
        else if (interfaceType == "F")
        {
            interface = FixedInterfaceShPtr(MemoryManager<FixedInterface>::AllocateSharedPtr(indx, domain));
        }

        interface->SetEdge(domainEdge);

        tmpInterfaceMap[indx].emplace_back(interface);

        m_interfaceVector.emplace_back(interface);

        interfaceElement = interfaceElement->NextSiblingElement();
    }

    for (auto interfacePair : tmpInterfaceMap)
    {
        ASSERTL0(interfacePair.second.size() == 2,
                 "Every interface ID must have two domains associated with it")

        interfacePair.second[0]->SetOppInterface(interfacePair.second[1]);
        interfacePair.second[1]->SetOppInterface(interfacePair.second[0]);

        m_interfaces[interfacePair.first] =
            MemoryManager<SpatialDomains::InterfacePair>::AllocateSharedPtr(
                interfacePair.second[0], interfacePair.second[1]);
    }
}

void InterfaceBase::SetEdge(const CompositeMap &edge)
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

void Interfaces::PerformMovement(NekDouble timeStep)
{
    for (auto interface : m_interfaceVector)
    {
        interface->Move(timeStep);
    }
}

void RotatingInterface::v_Move(NekDouble timeStep)
{
    NekDouble angle = m_angularVel * timeStep;

    NekDouble ox, oy, oz;
    m_origin.GetCoords(ox, oy, oz);

    for (auto &vert : m_rotateVerts)
    {
        NekDouble x, y, z;
        vert->GetCoords(x, y, z);

        vert->UpdatePosition(
            std::cos(angle) * (x - ox) - std::sin(angle) * (y - oy) + ox,
            std::sin(angle) * (x - ox) + std::cos(angle) * (y - oy) + oy,
            0);
    }

    for (auto &curve : m_rotateCurves)
    {
        for (auto &vert : curve->m_points)
        {
            NekDouble x, y, z;
            vert->GetCoords(x, y, z);

            vert->UpdatePosition(
                std::cos(angle) * (x - ox) - std::sin(angle) * (y - oy) + ox,
                std::sin(angle) * (x - ox) + std::cos(angle) * (y - oy) + oy,
                0);
        }
    }
}

void FixedInterface::v_Move(NekDouble time)
{
    boost::ignore_unused(time);
}

}
}