#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/Interface.h>

using namespace Nektar;
using namespace Nektar::SpatialDomains;

bool ContainsElement(CompositeMap &domain, int edgeID)
{
    for (auto &comp : domain)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            if (edgeID == geom->GetGlobalID())
            {
                return true;
            }
        }
    }

    return false;
}

std::vector<GeometrySharedPtr>
VertexToElmt(CompositeMap &domain, int vertId1, int vertId2)
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

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session = LibUtilities::SessionReader::CreateInstance(
            argc, argv);

    MeshGraphSharedPtr graph = MeshGraph::Read(session);
    Interfaces interfaceCollection = Interfaces(session, graph);

    std::map<int, InterfaceShPtr> interfaces = interfaceCollection.GetInterfaces();

    for (auto &interface : interfaces)
    {
        int indx = interface.first;
        auto interfaceProperties = interface.second;

        auto movingDomain = interfaceProperties->GetMovingDomain();
        auto interfaceEdge = interfaceProperties->GetInterfaceEdge();

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
        //Map that stores elements to process for renumbered edges and points
        std::map<int, GeometrySharedPtr> elementToDo;

        for (auto &comp : interfaceEdge)
        {
            for (auto &geom : comp.second->m_geomVec)
            {
                std::cout << "Processing edge ID = " << geom->GetGlobalID()
                          << std::endl;
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
                        vertDone[vid[i]] = maxVertId++;
                        std::cout << "Replacing vert " << vid[i] << " with "
                                  << maxVertId - 1 << std::endl;
                    }
                    else
                    {
                        newVerts[i] = graph->GetVertex(it->second);
                    }
                }

                SegGeomSharedPtr oldEdge = std::static_pointer_cast<SegGeom>(
                        geom);

                CurveSharedPtr newCurve;
                if (oldEdge->GetCurve())
                {
                    newCurve = MemoryManager<Curve>::AllocateSharedPtr(
                            maxEdgeId,
                            oldEdge->GetCurve()->m_ptype);
                }

                auto newEdge = MemoryManager<SegGeom>::AllocateSharedPtr(
                        maxEdgeId, newVerts[0]->GetCoordim(), newVerts,
                        newCurve);
                std::cout << "Creating new edge " << maxEdgeId << ": "
                          << newVerts[0]->GetGlobalID() << " "
                          << newVerts[1]->GetGlobalID() << std::endl << endl;
                graph->GetAllSegGeoms()[maxEdgeId] = newEdge;
                edgeDone[geom->GetGlobalID()] = newEdge;
                maxEdgeId++;

                auto toProcess = VertexToElmt(movingDomain, vid[0], vid[1]);
                for (auto &elementToProcess : toProcess)
                {
                    elementToDo[elementToProcess->GetGlobalID()] = elementToProcess;
                }
            }
        }


        for (auto &elementMap : elementToDo)
        {
            auto movingGeom = elementMap.second;

            cout << endl << "Looking at element: " << movingGeom->GetGlobalID()
                 << endl;
            std::vector<SegGeomSharedPtr> newEdges(
                    movingGeom->GetNumEdges());

            // Loop over edges
            for (int j = 0; j < newEdges.size(); ++j)
            {
                auto edge = std::static_pointer_cast<SegGeom>(
                        movingGeom->GetEdge(j));
                cout << "Looking at edge: " << edge->GetGlobalID();
                auto edgeIt = edgeDone.find(edge->GetGlobalID());
                if (edgeIt != edgeDone.end())
                {
                    cout << " - already redefined to edge " << edgeIt->second->GetGlobalID() << endl;
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
                        create = true;
                        cout << " - redefine edge vertex " << vertIt->first << " to "<< vertIt->second << endl;
                    }
                    else newEdgeVerts[k] = graph->GetVertex(edgeVids[k]);
                }

                if (create)
                {
                    auto newEdge = MemoryManager<SegGeom>::AllocateSharedPtr(
                            edge->GetGlobalID(),
                            edge->GetVertex(0)->GetCoordim(),
                            newEdgeVerts, edge->GetCurve());
                    graph->GetAllSegGeoms()[edge->GetGlobalID()] = newEdge;
                    edgeDone[edge->GetGlobalID()] = newEdge;
                    newEdges[j] = newEdge;
                }
                else
                {
                    newEdges[j] = edge;
                    cout << " - keep old edge vertices" << endl;
                }
            }

            if (movingGeom->GetShapeType() == LibUtilities::eQuadrilateral)
            {
                // Create a new quad
                QuadGeomSharedPtr quad = std::static_pointer_cast<QuadGeom>(
                        movingGeom);
                QuadGeomSharedPtr newQuad = MemoryManager<QuadGeom>::AllocateSharedPtr(
                        quad->GetGlobalID(), &newEdges[0], quad->GetCurve());
                graph->GetAllQuadGeoms()[quad->GetGlobalID()] = newQuad;
                std::cout << "Redefining element: " << quad->GetGlobalID()
                          << std::endl;
            }
            else if (movingGeom->GetShapeType() == LibUtilities::eTriangle)
            {
                // Create a new tri
                TriGeomSharedPtr tri = std::static_pointer_cast<TriGeom>(
                        movingGeom);
                TriGeomSharedPtr newTri = MemoryManager<TriGeom>::AllocateSharedPtr(
                        tri->GetGlobalID(), &newEdges[0], tri->GetCurve());
                graph->GetAllTriGeoms()[tri->GetGlobalID()] = newTri;
            }
        }

        std::set<int> seenVerts, seenEdges;
        for (auto &comp : interfaceProperties->GetMovingDomain())
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

                    (*vert)(1) += 10.0;
                    seenVerts.insert(vert->GetGlobalID());
                }
            }
        }
    }

    std::string filename = "out.xml";
    graph->WriteGeometry(filename, true);
}
