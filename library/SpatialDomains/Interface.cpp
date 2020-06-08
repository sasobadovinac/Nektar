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

std::vector<GeometrySharedPtr> GetElementsFromVertex(CompositeMap &domain,
                                                     int vertId1, int vertId2)
{
    std::vector<GeometrySharedPtr> ret;
    for (auto &comp : domain)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            for (int i = 0; i < geom->GetNumVerts(); ++i)
            {
                if (geom->GetVid(i) == vertId1 || geom->GetVid(i) == vertId2)
                {
                    ret.push_back(geom);
                    break;
                }
            }
        }
    }
    return ret;
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


        std::string interfaceEdgeStr;
        int interfaceErr = interfaceElement->QueryStringAttribute(
            "INTERFACEEDGE", &interfaceEdgeStr);
        map<int, CompositeSharedPtr> interfaceEdge;
        if (interfaceErr == TIXML_SUCCESS)
        {
            std::string indxStr = ReadTag(interfaceEdgeStr);
            m_meshGraph->GetCompositeList(indxStr, interfaceEdge);
        }

        std::string domainEdgeStr;
        int domainEdgeErr =
            interfaceElement->QueryStringAttribute("EDGE", &domainEdgeStr);
        map<int, CompositeSharedPtr> domainEdge;
        if (domainEdgeErr == TIXML_SUCCESS)
        {
            std::string indxStr = ReadTag(domainEdgeStr);
            m_meshGraph->GetCompositeList(indxStr, domainEdge);
        }

        if (interfaceErr == TIXML_SUCCESS)
        {
            ASSERTL0(domainEdgeErr != TIXML_SUCCESS,
                     "Choose to define either INTERFACEEDGE or EDGE")
        }
        else if (domainEdgeErr == TIXML_SUCCESS)
        {
            ASSERTL0(interfaceErr != TIXML_SUCCESS,
                     "Choose to define either INTERFACEEDGE or EDGE")

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
            err = interfaceElement->QueryStringAttribute("ANGVEL",
                                                         &angularVelStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read angular velocity.");

            NekDouble angularVel = stod(angularVelStr);

            interface = RotatingInterfaceShPtr(MemoryManager<RotatingInterface>::AllocateSharedPtr(domain,  origin, axis, angularVel));
        }
        else if (interfaceType == "F")
        {
            interface = FixedInterfaceShPtr(MemoryManager<FixedInterface>::AllocateSharedPtr(domain));
        }

        if (interfaceErr == TIXML_SUCCESS)
        {
            interface->SetInterfaceEdge(interfaceEdge);
        }
        else
        {
            interface->SetEdge(domainEdge);
        }

        tmpInterfaceMap[indx].emplace_back(interface);

        interfaceElement = interfaceElement->NextSiblingElement();
    }

    for (auto interfacePair : tmpInterfaceMap)
    {
        ASSERTL0(interfacePair.second.size() == 2,
                 "Every interface ID must have two domains associated with it")

        m_interfaces[interfacePair.first] =
            MemoryManager<SpatialDomains::InterfacePair>::AllocateSharedPtr(
                interfacePair.second[0], interfacePair.second[1]);
    }
}

void Interfaces::CommSetup(LibUtilities::CommSharedPtr &comm)
{
    // myIndxLR contains the info about what interface edges are present on
    // current rank with each interface no,  i,  consisting of:
    // [i] = indx
    // [i + 1] = 0 (non), = 1 (left only), = 2 (right only), = 3 (both)

    Array<OneD, int> myIndxLR(m_interfaces.size() * 2, 0);
    std::map<int, int> myIndxLRMap;
    size_t cnt = 0;
    for (const auto &interface : m_interfaces)
    {
        myIndxLR[2 * cnt] = interface.first;

        if(!interface.second->GetLeftInterface()->IsEmpty())
        {
            myIndxLR[2 * cnt + 1] += 1;
        }
        if (!interface.second->GetRightInterface()->IsEmpty())
        {
            myIndxLR[2 * cnt + 1] += 2;
        }

        myIndxLRMap[interface.first] = myIndxLR[2 * cnt + 1];
        cnt++;
    }


    //Send num of interfaces size so all partitions can prepare buffers
    int nRanks = comm->GetSize();
    Array<OneD, int> rankNumInterfaces(nRanks);
    Array<OneD, int> localInterfaceSize(1, myIndxLR.size());
    comm->AllGather(localInterfaceSize, rankNumInterfaces);

    Array<OneD, int> rankLocalInterfaceDisp(nRanks, 0);
    for (size_t i = 1; i < nRanks; ++i)
    {
        rankLocalInterfaceDisp[i] = rankLocalInterfaceDisp[i - 1] + rankNumInterfaces[i - 1];
    }

    Array<OneD, int> rankLocalInterfaceIds(
        std::accumulate(rankNumInterfaces.begin(), rankNumInterfaces.end(), 0), 0);

    // Send all interface IDs to all partitions
    comm->AllGatherv(myIndxLR, rankLocalInterfaceIds, rankNumInterfaces,
                     rankLocalInterfaceDisp);

    // Find what interface Ids match with other ranks, then check if opposite edge
    size_t myRank = comm->GetRank();
    for (size_t i = 0; i < nRanks; ++i)
    {
        if (i == myRank)
        {
            continue;
        }

        for (size_t j = 0; j < rankNumInterfaces[i] / 2; ++j)
        {
            int otherId =
                rankLocalInterfaceIds[rankLocalInterfaceDisp[i] + 2 * j];
            int otherCode =
                rankLocalInterfaceIds[rankLocalInterfaceDisp[i] + 2 * j + 1];
            if (myIndxLRMap.find(otherId) != myIndxLRMap.end())
            {
                // Set interface opposite ranks (could probably simplify logic
                // here but this is easy to understand
                int myCode = myIndxLRMap[otherId];
                if ((myCode == 1 && otherCode == 2) ||
                    (myCode == 1 && otherCode == 3) ||
                    (myCode == 3 && otherCode == 2))
                {
                    m_interfaces[otherId]->GetLeftInterface()->AddOppRank(i);
                }
                else if ((myCode == 2 && otherCode == 1) ||
                         (myCode == 2 && otherCode == 3) ||
                         (myCode == 3 && otherCode == 1))
                {
                    m_interfaces[otherId]->GetRightInterface()->AddOppRank(i);
                }
                else if (myCode == 3 && otherCode == 3)
                {
                    m_interfaces[otherId]->GetLeftInterface()->AddOppRank(i);
                    m_interfaces[otherId]->GetRightInterface()->AddOppRank(i);
                }
            }
        }
    }
}

void InterfacePair::SeparateGraph(MeshGraphSharedPtr &graph) const
{
    auto rightDomain   = m_rightInterface->GetDomain();
    auto interfaceEdge = m_rightInterface->GetInterfaceEdge();

    int maxVertId = -1;
    for (auto &vert : graph->GetAllPointGeoms())
    {
        maxVertId = std::max(maxVertId, vert.first);
    }

    int maxEdgeId = -1;
    for (auto &edge : graph->GetAllSegGeoms())
    {
        maxEdgeId = std::max(maxEdgeId, edge.first);
    }

    ++maxVertId;
    ++maxEdgeId;

    // Map that stores existing renumbered geometry.
    std::map<int, int> vertDone;
    std::map<int, SegGeomSharedPtr> edgeDone;
    // Map that stores elements to process for renumbered edges and points
    std::map<int, GeometrySharedPtr> elementToDo;

    CurveMap &curvedEdges = graph->GetCurvedEdges();

    for (auto &comp : interfaceEdge)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            ASSERTL0(geom->GetShapeType() == LibUtilities::eSegment,
                     "Unexpected geometry type in composite");

            GeometryLinkSharedPtr elmtLink = graph->GetElementsFromEdge(
                std::static_pointer_cast<Geometry1D>(geom));

            size_t numElmts = elmtLink->size();
            if (numElmts == 1)
            {
                continue;
            }

            int vid[2] = {geom->GetVid(0), geom->GetVid(1)};
            PointGeomSharedPtr newVerts[2];

            for (int i = 0; i < 2; ++i)
            {
                auto it = vertDone.find(vid[i]);
                if (it == vertDone.end())
                {
                    // Create a new vertex
                    newVerts[i] = MemoryManager<PointGeom>::AllocateSharedPtr(
                        *geom->GetVertex(i));
                    newVerts[i]->SetGlobalID(maxVertId);
                    graph->GetAllPointGeoms()[maxVertId] = newVerts[i];
                    vertDone[vid[i]]                     = maxVertId++;
                }
                else
                {
                    newVerts[i] = graph->GetVertex(it->second);
                }
            }

            SegGeomSharedPtr oldEdge = std::static_pointer_cast<SegGeom>(geom);

            CurveSharedPtr newCurve;
            if (oldEdge->GetCurve())
            {
                newCurve = MemoryManager<Curve>::AllocateSharedPtr(
                    maxEdgeId, oldEdge->GetCurve()->m_ptype);

                for (auto &pt : oldEdge->GetCurve()->m_points)
                {
                    newCurve->m_points.push_back(
                        MemoryManager<PointGeom>::AllocateSharedPtr(*pt));
                }

                curvedEdges[maxEdgeId] = newCurve;
            }

            auto newEdge = MemoryManager<SegGeom>::AllocateSharedPtr(
                maxEdgeId, newVerts[0]->GetCoordim(), newVerts, newCurve);

            graph->GetAllSegGeoms()[maxEdgeId] = newEdge;
            edgeDone[geom->GetGlobalID()]      = newEdge;

            m_leftInterface->SetEdge(oldEdge);
            m_rightInterface->SetEdge(newEdge);
            maxEdgeId++;

            auto toProcess = GetElementsFromVertex(rightDomain, vid[0], vid[1]);
            for (auto &elementToProcess : toProcess)
            {
                elementToDo[elementToProcess->GetGlobalID()] = elementToProcess;
            }
        }
    }

    for (auto &elementMap : elementToDo)
    {
        auto rightGeom = elementMap.second;

        std::vector<SegGeomSharedPtr> newEdges(rightGeom->GetNumEdges());

        // Loop over edges
        for (int j = 0; j < newEdges.size(); ++j)
        {
            auto edge =
                std::static_pointer_cast<SegGeom>(rightGeom->GetEdge(j));
            auto edgeIt = edgeDone.find(edge->GetGlobalID());
            if (edgeIt != edgeDone.end())
            {
                newEdges[j] = edgeIt->second;
                continue;
            }

            int edgeVids[2] = {edge->GetVid(0), edge->GetVid(1)};

            PointGeomSharedPtr newEdgeVerts[2];
            bool create = false;

            for (int k = 0; k < 2; ++k)
            {
                auto vertIt = vertDone.find(edgeVids[k]);
                if (vertIt != vertDone.end())
                {
                    newEdgeVerts[k] = graph->GetVertex(vertIt->second);
                    create          = true;
                }
                else
                    newEdgeVerts[k] = graph->GetVertex(edgeVids[k]);
            }

            if (create)
            {
                auto newEdge = MemoryManager<SegGeom>::AllocateSharedPtr(
                    edge->GetGlobalID(), edge->GetVertex(0)->GetCoordim(),
                    newEdgeVerts, edge->GetCurve());
                graph->GetAllSegGeoms()[edge->GetGlobalID()] = newEdge;
                edgeDone[edge->GetGlobalID()]                = newEdge;
                newEdges[j]                                  = newEdge;
            }
            else
            {
                newEdges[j] = edge;
            }
        }

        GeometrySharedPtr newGeom;
        if (rightGeom->GetShapeType() == LibUtilities::eQuadrilateral)
        {
            // Create a new quad
            QuadGeomSharedPtr quad =
                std::static_pointer_cast<QuadGeom>(rightGeom);
            QuadGeomSharedPtr newQuad =
                MemoryManager<QuadGeom>::AllocateSharedPtr(
                    quad->GetGlobalID(), &newEdges[0], quad->GetCurve());
            graph->GetAllQuadGeoms()[quad->GetGlobalID()] = newQuad;
            newGeom                                       = newQuad;
        }
        else if (rightGeom->GetShapeType() == LibUtilities::eTriangle)
        {
            // Create a new tri
            TriGeomSharedPtr tri = std::static_pointer_cast<TriGeom>(rightGeom);
            TriGeomSharedPtr newTri = MemoryManager<TriGeom>::AllocateSharedPtr(
                tri->GetGlobalID(), &newEdges[0], tri->GetCurve());
            graph->GetAllTriGeoms()[tri->GetGlobalID()] = newTri;
            newGeom                                     = newTri;
        }

        // Replace this geometry in any composites.
        for (auto &comp : graph->GetComposites())
        {
            auto tmp = comp.second->m_geomVec.size();
            for (int n = 0; n < tmp; ++n)
            {
                if (comp.second->m_geomVec[n]->GetGlobalID() ==
                        newGeom->GetGlobalID() &&
                    comp.second->m_geomVec[n]->GetShapeType() ==
                        newGeom->GetShapeType())
                {
                    comp.second->m_geomVec[n] = newGeom;
                }
            }
        }

        //Offset right domain for visualisation
#if 0
        std::set<int> seenVerts, seenEdges;
        for (auto &comp : GetRightDomain())
        {
            for (auto &geom : comp.second->m_geomVec)
            {
                auto newGeom = graph->GetGeometry2D(geom->GetGlobalID());
                for (int i = 0; i < newGeom->GetNumVerts(); ++i)
                {
                    PointGeomSharedPtr vert = newGeom->GetVertex(i);

                    if (seenVerts.find(vert->GetGlobalID()) != seenVerts.end())
                    {
                        continue;
                    }

                    (*vert)(1) += 0.5;
                    seenVerts.insert(vert->GetGlobalID());

                }

                for (int i = 0; i < newGeom->GetNumEdges(); ++i)
                {
                    SegGeomSharedPtr edge =
                        std::static_pointer_cast<SegGeom>(newGeom->GetEdge(i));

                    // move curve points
                    if (seenEdges.find(edge->GetGlobalID()) != seenEdges.end())
                    {
                        continue;
                    }

                    CurveSharedPtr curve = edge->GetCurve();
                    if (!curve)
                    {
                        continue;
                    }

                    for (auto &pt : curve->m_points)
                    {
                        (*pt)(1) += 0.5;
                    }

                    seenEdges.insert(edge->GetGlobalID());
                }
            }
        }
#endif

    }
}

void InterfaceBase::SetEdge(const CompositeMap &edge)
{
    for (auto &compIt : edge)
    {
        for (auto &elmtIt : compIt.second->m_geomVec)
        {
            SegGeomSharedPtr elmt = std::dynamic_pointer_cast<SegGeom>(elmtIt);

            ASSERTL0(elmt,
                     "Composite for edge on the interface should only contain "
                     "segments");

            size_t id = elmt->GetGlobalID();
            m_edge[id] = elmt;
            m_edgeIds.emplace_back(id);
        }
    }

    std::sort(m_edgeIds.begin(), m_edgeIds.end());
}

}
}