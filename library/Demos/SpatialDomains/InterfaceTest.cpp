//
// Created by Edward Laughton on 2019-03-12.
//

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/Interface.h>

using namespace Nektar;
using namespace SpatialDomains;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session = LibUtilities::SessionReader::CreateInstance(argc, argv);

    SpatialDomains::MeshGraphSharedPtr graph = SpatialDomains::MeshGraph::Read(session);
    Interfaces interfaceCollection = SpatialDomains::Interfaces(session, graph);

    auto regions = interfaceCollection.GetInterfaces();

    for (auto &region : regions)
    {
        cout << std::to_string(region.second->GetInterfaceType()) << endl; //<--- works
        // cout << std::to_string(region.second->GetAngularVel()) << endl; <--- doesn't work
    }
    cout << "All read" << endl;
}