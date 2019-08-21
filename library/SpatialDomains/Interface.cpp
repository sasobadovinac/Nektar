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

    std::map<int, std::vector<InterfaceShPtr>> tmpInterfaceMap;

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
            "INTERFACE", &interfaceEdgeStr);
        map<int, CompositeSharedPtr> interfaceEdge;
        if (interfaceErr == TIXML_SUCCESS)
        {
            std::string indxStr = ReadTag(interfaceEdgeStr);
            m_meshGraph->GetCompositeList(indxStr, interfaceEdge);
        }

        std::string domainEdgeStr;
        int domainEdgeErr = interfaceElement->QueryStringAttribute(
                "EDGE", &domainEdgeStr);
        map<int, CompositeSharedPtr> domainEdge;
        if (domainEdgeErr == TIXML_SUCCESS)
        {
            std::string indxStr = ReadTag(domainEdgeStr);
            m_meshGraph->GetCompositeList(indxStr, domainEdge);
        }

        if (interfaceErr == TIXML_SUCCESS)
        {
            ASSERTL0(domainEdgeErr != TIXML_SUCCESS,
                     "Choose to define either INTERFACE or EDGE")
        }
        else if (domainEdgeErr == TIXML_SUCCESS)
        {
            ASSERTL0(interfaceErr != TIXML_SUCCESS,
                     "Choose to define either INTERFACE or EDGE")

        }
        else
        {
            ASSERTL0(false, "Choose to define either INTERFACE or EDGE")
        }

        InterfaceShPtr interface;

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

            interface = RotatingInterfaceShPtr(MemoryManager<RotatingInterface>::AllocateSharedPtr(domain, origin, axis, angularVel));
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
        ASSERTL0(interfacePair.second.size() == 2, "Every interface ID must have two domains associated with it")

           m_interfaces[interfacePair.first] = std::make_pair(
                   interfacePair.second[0], interfacePair.second[1]);
    }
}

void Interfaces::SeparateGraph(MeshGraphSharedPtr &graph, int indx)
{
    auto &keepInterface = m_interfaces[indx].first;
    auto &changeInterface = m_interfaces[indx].second;
    auto rightDomain   = changeInterface->GetDomain();
    auto interfaceEdge = changeInterface->GetInterfaceEdge();

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

            keepInterface->SetEdge(oldEdge);
            changeInterface->SetEdge(newEdge);
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
    }
};

void Interfaces::GenerateMortars(int indx)
{
    //set up 'right' edge points set
    auto &iterateEdge = m_interfaces[indx].second->GetEdge();
    std::map<int, PointGeomSharedPtr> points;

    for (auto &seg : iterateEdge)
    {
        auto vert0 = MemoryManager<PointGeom>::AllocateSharedPtr(
                *seg.second->GetVertex(0));
        auto vert1 = MemoryManager<PointGeom>::AllocateSharedPtr(
                *seg.second->GetVertex(1));

        points[vert0->GetVid()] = vert0;
        points[vert1->GetVid()] = vert1;
    }

    for (auto tmp : points)
    {
        cout << tmp.first << endl;
    }

    //mortars begin as copy of left side edges, change to vector and increment global ID
    auto &otherEdge = m_interfaces[indx].first->GetEdge();
    std::vector<SegGeomSharedPtr> mortars(otherEdge.size());
    int cnt = 0;
    for (auto iter : otherEdge)
    {
        auto edge = MemoryManager<SegGeom>::AllocateSharedPtr(*iter.second);
        edge->SetGlobalID(cnt);
        edge->GetGeomFactors()->GetGtype(); //populate Gtype for FindDistance
        mortars[cnt] = edge;
        cnt++;
    }

    //iterate over 'right' points and find matching seg in 'left' side
    Array<OneD, NekDouble> xs(3);
    for (auto &vert : points)
    {
        vert.second->GetCoords(xs);
        NekDouble foundPoint;
        for (int i = 0; i < mortars.size(); ++i)
        {
            auto geomSeg = mortars[i];
            NekDouble dist = geomSeg->FindDistance(xs, foundPoint);

            if (dist > 1e-8)
            {
                continue;
            }

            mortars.erase(mortars.begin() + i); //remove old split segment from mortar

            //create two new mortars splitting found seg around trial point
            bool indxFlag = false; //flag for if old mortar global ID has already been replaced
            NekDouble xOld = 0, yOld = 0, xNew = 0, yNew = 0, zOld = 0, zNew = 0;
            for (int k = 0; k < 2; ++k)
            {
                PointGeomSharedPtr newEdgeVerts[2];
                newEdgeVerts[0] = MemoryManager<PointGeom>::AllocateSharedPtr(
                        *vert.second);
                newEdgeVerts[1] = MemoryManager<PointGeom>::AllocateSharedPtr(
                        *geomSeg->GetVertex(k));

                newEdgeVerts[0]->GetCoords(xOld, yOld, zOld);
                newEdgeVerts[1]->GetCoords(xNew, yNew, zNew);

                if (yOld > yNew)
                {
                    std::swap(newEdgeVerts[0], newEdgeVerts[1]);
                }

                int id;
                if (!(xOld == xNew && yOld == yNew)) //check that verts aren't same
                {
                    if (!indxFlag)
                    {
                        id = geomSeg->GetGlobalID();
                        indxFlag = true;
                    }
                    else
                    {
                        id = cnt;
                        cnt++;
                    }

                    auto newMortar = MemoryManager<SegGeom>::AllocateSharedPtr(
                            id, 2, newEdgeVerts, geomSeg->GetCurve());
                    newMortar->GetGeomFactors()->GetGtype(); //populate new
                    mortars.emplace_back(newMortar);
                }
            }

            break;
        }
    }

    m_mortars.insert(m_mortars.end(), mortars.begin(), mortars.end());

    //Debugging output mortar coords
    for (auto it : mortars)
    {
        NekDouble x, y, z, a, b, c;
        it->GetVertex(0)->GetCoords(x, y, z);
        it->GetVertex(1)->GetCoords(a, b, c);

        cout << "Coords for mortar: " << it->GetGlobalID() <<  " | (" << x << ", " << y << ")" << " -> (" << a << ", " << b << ")" << endl;
    }

    //iterate over mortar and identify left and right corresponding edge
    //by evaluating center of standard element into Cartesian
    Array<OneD, NekDouble> xi(1, 0.0);
    cnt = 0;
    for (const auto it : mortars)
    {
        cout << endl << "Mortar: " << it->GetGlobalID() << endl;
        int nq = it->GetXmap()->GetTotPoints();
        Array<OneD, NekDouble> x(nq), y(nq);
        it->GetXmap()->BwdTrans(it->GetCoeffs(0),x);
        it->GetXmap()->BwdTrans(it->GetCoeffs(1),y);

        NekDouble xc = it->GetXmap()->PhysEvaluate(xi, x);
        NekDouble yc = it->GetXmap()->PhysEvaluate(xi, y);

        cout << "Centre: (" << xc << ", " << yc << ")" << endl;

        NekDouble foundPoint;
        Array<OneD, NekDouble> xs(3);
        xs[0] = xc;
        xs[1] = yc;
        xs[2] = 0.0;

        for (const auto edge: iterateEdge)
        {
            edge.second->GetGeomFactors()->GetGtype(); //populate Gtype for FindDistance
            NekDouble dist = edge.second->FindDistance(xs, foundPoint);
            if (dist > 1e-8)
            {
                continue;
            }

            m_mortarToRightEdgeMap.emplace_back(edge.second->GetGlobalID());
            m_rightEdgeToMortarMap[edge.second->GetGlobalID()].emplace_back(cnt);

            cout << "'Right' edge | Segment ID: " << edge.second->GetGlobalID() << endl;
        }

        for (const auto edge: otherEdge)
        {
            edge.second->GetGeomFactors()->GetGtype(); //populate Gtype for FindDistance
            NekDouble dist = edge.second->FindDistance(xs, foundPoint);
            if (dist > 1e-8)
            {
                continue;
            }

            m_mortarToLeftEdgeMap.emplace_back(edge.second->GetGlobalID());
            m_leftEdgeToMortarMap[edge.second->GetGlobalID()].emplace_back(cnt);

            cout << "'Left' edge | Segment ID: " << edge.second->GetGlobalID() << endl;
        }

        cnt++;
    }


    ASSERTL0(m_mortarToLeftEdgeMap.size() == m_mortars.size(),
            "Length of mortar to left edge map is not equal to number of mortars.");

    ASSERTL0(m_mortarToRightEdgeMap.size() == m_mortars.size(),
            "Length of mortar to right edge map is not equal to number of mortars.");

    cout << endl << "Number of mortars: " << m_mortars.size() << endl << endl;

    for (auto tmp : m_rightEdgeToMortarMap)
    {
        cout << "Right edge " << tmp.first << " -> ";
        for (auto tmp2: tmp.second)
        {
            cout << tmp2 << ' ';
        }
        cout << endl;
    }

    cout << endl;
    for (auto tmp : m_leftEdgeToMortarMap)
    {
        cout << "Left edge " << tmp.first << " -> ";
        for (auto tmp2: tmp.second)
        {
            cout << tmp2 << ' ';
        }
        cout << endl;
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
                     "Composite for left edge should only contain segments");
            m_edge[elmt->GetGlobalID()] = elmt;
        }
    }
}

void InterfaceBase::FillInterfaceBoundingBoxTree()
{
    if (m_boundingInterface.empty())
    {
        for (auto &x : m_edge)
        {
            BgBox b = x.second->GetBoundingBox();
            m_boundingInterface.insert(std::make_pair(b, x.first));
        }
    }
}

std::vector<BgRtreeValue> InterfaceBase::GetEdgesContainingPoint(NekDouble x,
        NekDouble y, NekDouble z)
{
    if (m_boundingInterface.empty())
    {
        FillInterfaceBoundingBoxTree();
    }

    std::vector<BgRtreeValue> vals;

    BgBox b(BgPoint(x, y, z), BgPoint(x, y, z));

    m_boundingInterface.query(bg::index::intersects(b),
                                  std::back_inserter(vals));

    return vals;
}
}
}
