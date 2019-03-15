#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/Interface.h>

using namespace Nektar;
using namespace Nektar::SpatialDomains;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session = LibUtilities::SessionReader::CreateInstance(argc, argv);

    MeshGraphSharedPtr graph = MeshGraph::Read(session);
    Interfaces interface = Interfaces(session, graph);
    std::map<int, InterfaceShPtr> interfaces = interface.GetInterfaces();

    for (auto &interface : interfaces)
    {
        int indx = interface.first;
        auto interfaceProperties = interface.second;

        InterfaceEdgeMap interfaceEdgeCompositeCollection = interfaceProperties->GetInterfaceEdge();
        std::vector<int> interfaceVertexIds;
        std::vector<int> interfaceEdgeIds;
        for (auto &interfaceEdgeCompositeMap : interfaceEdgeCompositeCollection)
        {
            CompositeSharedPtr interfaceEdgeComposite = interfaceEdgeCompositeMap.second;
            std::vector<GeometrySharedPtr> geomVec = interfaceEdgeComposite->m_geomVec;
            for (GeometrySharedPtr &geom : geomVec)
            {
                interfaceEdgeIds.emplace_back(geom->GetGlobalID());
                interfaceVertexIds.emplace_back (geom->GetVertex(0)->GetGlobalID());
                interfaceVertexIds.emplace_back (geom->GetVertex(1)->GetGlobalID());
            }
        }

        sort(interfaceVertexIds.begin(), interfaceVertexIds.end());
        interfaceVertexIds.erase(unique(interfaceVertexIds.begin(), interfaceVertexIds.end()), interfaceVertexIds.end());

        std::vector<std::map<int, CompositeSharedPtr>> domains = graph->GetDomain();
        std::map<int, CompositeSharedPtr> movingDomain;
        std::map<int, CompositeSharedPtr> fixedDomain;

        for (auto domain : domains)
        {
            if (domain == interfaceProperties->GetMovingDomain())
            {
                domain = movingDomain;
            }
            else if (domain == interfaceProperties->GetFixedDomain())
            {
                domain = fixedDomain;
            }
        }

        for (auto composite : movingDomain)
        {
            //
        }
    }
}
