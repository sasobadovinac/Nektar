//
// Created by Edward Laughton on 2019-03-12.
//

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/Interface.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session =
            LibUtilities::SessionReader::CreateInstance(argc, argv);

    SpatialDomains::MeshGraphSharedPtr graph =
            SpatialDomains::MeshGraph::Read(session);

    Interface conds = SpatialDomains::Interface(session, graph);
    auto regions = conds.GetInterfaceRegions();

    /*for (auto &region : regions)
    {
        std::cout << region->GetInterfaceConditionType() << std::endl;
        region->AdjustMeshGraph(graph);
    }*/


}