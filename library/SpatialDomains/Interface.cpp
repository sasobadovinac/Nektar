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
    ReadInterfaces(m_session->GetElement("NEKTAR/CONDITIONS/INTERFACES"));
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

    while (interfaceElement)
    {
        std::string interfaceType = interfaceElement->Value();

        int err;
        int indx;

        err = interfaceElement->QueryIntAttribute("ID", &indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read interface ID.");

        std::string interfaceRightDomainStr;
        err = interfaceElement->QueryStringAttribute("RIGHTDOMAIN",
                                                     &interfaceRightDomainStr);
        ASSERTL0(err == TIXML_SUCCESS,
                 "Unable to read right interface domain.");
        auto rightDomain =
            m_meshGraph->GetDomain(stoi(ReadTag(interfaceRightDomainStr)));

        std::string interfaceLeftDomainStr;
        err = interfaceElement->QueryStringAttribute("LEFTDOMAIN",
                                                     &interfaceLeftDomainStr);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read left interface domain.");
        auto leftDomain =
            m_meshGraph->GetDomain(stoi(ReadTag(interfaceLeftDomainStr)));

        std::string interfaceEdgeStr;
        int interfaceErr = interfaceElement->QueryStringAttribute(
            "INTERFACE", &interfaceEdgeStr);
        map<int, CompositeSharedPtr> interfaceEdge;
        if (interfaceErr == TIXML_SUCCESS)
        {
            std::string indxStr = ReadTag(interfaceEdgeStr);
            m_meshGraph->GetCompositeList(indxStr, interfaceEdge);
        }

        std::string leftEdgeStr;
        int leftEdgeErr =
            interfaceElement->QueryStringAttribute("LEFTEDGE", &leftEdgeStr);
        map<int, CompositeSharedPtr> leftEdge;
        if (leftEdgeErr == TIXML_SUCCESS)
        {
            std::string indxStr = ReadTag(leftEdgeStr);
            m_meshGraph->GetCompositeList(indxStr, leftEdge);
        }

        std::string rightEdgeStr;
        int rightEdgeErr =
            interfaceElement->QueryStringAttribute("RIGHTEDGE", &rightEdgeStr);
        map<int, CompositeSharedPtr> rightEdge;
        if (rightEdgeErr == TIXML_SUCCESS)
        {
            std::string indxStr = ReadTag(rightEdgeStr);
            m_meshGraph->GetCompositeList(indxStr, rightEdge);
        }

        if (interfaceErr == TIXML_SUCCESS)
        {
            ASSERTL0(leftEdgeErr != TIXML_SUCCESS &&
                         rightEdgeErr != TIXML_SUCCESS,
                     "Choose to define either INTERFACE or both LEFTEDGE "
                     "and RIGHTEDGE.")
        }
        else if (leftEdgeErr == TIXML_SUCCESS && rightEdgeErr == TIXML_SUCCESS)
        {
            ASSERTL0(interfaceErr != TIXML_SUCCESS,
                     "Choose to define either INTERFACE or both LEFTEDGE "
                     "and RIGHTEDGE.")
        }
        else
        {
            ASSERTL0((interfaceErr + 1) * (leftEdgeErr + 1) *
                             (rightEdgeErr + 1) ==
                         1,
                     "Choose to define either INTERFACE or both LEFTEDGE "
                     "and RIGHTEDGE.")
        }

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

            InterfaceShPtr rotatingInterface(
                MemoryManager<RotatingInterface>::AllocateSharedPtr(
                    rightDomain, leftDomain, origin, axis, angularVel));

            if (interfaceErr == TIXML_SUCCESS)
            {
                rotatingInterface->SetInterfaceEdge(interfaceEdge);
            }
            else
            {
                rotatingInterface->SetEdgeRight(rightEdge);
                rotatingInterface->SetEdgeLeft(leftEdge);
            }

            m_interfaces[indx] = rotatingInterface;
        }

        interfaceElement = interfaceElement->NextSiblingElement();
    }
}

void InterfaceBase::SeparateGraph(MeshGraphSharedPtr &graph)
{

    auto rightDomain   = GetRightDomain();
    auto interfaceEdge = GetInterfaceEdge();

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

            m_leftEdge[oldEdge->GetGlobalID()] = oldEdge;
            m_rightEdge[maxEdgeId]             = newEdge;
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
            for (int n = 0; n < comp.second->m_geomVec.size(); ++n)
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
    }
};

void InterfaceBase::SetEdgeLeft(const CompositeMap &leftEdge)
{
    for (auto &compIt : leftEdge)
    {
        for (auto &elmtIt : compIt.second->m_geomVec)
        {
            SegGeomSharedPtr elmt = std::dynamic_pointer_cast<SegGeom>(elmtIt);
            ASSERTL0(elmt,
                     "Composite for left edge should only contain segments");
            m_leftEdge[elmt->GetGlobalID()] = elmt;
        }
    }
}

void InterfaceBase::SetEdgeRight(const CompositeMap &rightEdge)
{
    for (auto &compIt : rightEdge)
    {
        for (auto &elmtIt : compIt.second->m_geomVec)
        {
            SegGeomSharedPtr elmt = std::dynamic_pointer_cast<SegGeom>(elmtIt);
            ASSERTL0(elmt,
                     "Composite for right edge should only contain segments");
            m_rightEdge[elmt->GetGlobalID()] = elmt;
        }
    }
}

void InterfaceBase::FillInterfaceBoundingBoxTree()
{
    if (m_boundingInterfaceRight.empty())
    {
        for (auto &x : m_rightEdge)
        {
            BgBox b = x.second->GetBoundingBox();
            m_boundingInterfaceRight.insert(std::make_pair(b, x.first));
        }
    }
    if (m_boundingInterfaceLeft.empty())
    {
        for (auto &x : m_leftEdge)
        {
            BgBox b = x.second->GetBoundingBox();
            m_boundingInterfaceLeft.insert(std::make_pair(b, x.first));
        }
    }
}

std::vector<BgRtreeValue> InterfaceBase::GetLeftEdgesContainingPoint(
    NekDouble x, NekDouble y, NekDouble z)
{
    if (m_boundingInterfaceLeft.empty())
    {
        FillInterfaceBoundingBoxTree();
    }

    std::vector<BgRtreeValue> vals;

    BgBox b(BgPoint(x, y, z), BgPoint(x, y, z));

    m_boundingInterfaceLeft.query(bg::index::intersects(b),
                                  std::back_inserter(vals));

    return vals;
}

std::vector<BgRtreeValue> InterfaceBase::GetRightEdgesContainingPoint(
    NekDouble x, NekDouble y, NekDouble z)
{
    if (m_boundingInterfaceRight.empty())
    {
        FillInterfaceBoundingBoxTree();
    }

    std::vector<BgRtreeValue> vals;

    BgBox b(BgPoint(x, y, z), BgPoint(x, y, z));

    m_boundingInterfaceRight.query(bg::index::intersects(b),
                                   std::back_inserter(vals));

    return vals;
}
} // namespace SpatialDomains
} // namespace Nektar
