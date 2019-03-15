//
// Created by Edward Laughton on 2019-03-12.
//

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
    Interfaces interfaceCollection = Interfaces(session, graph);

    auto interfaces = interfaceCollection.GetInterfaces();

    for (auto &interface : interfaces)
    {
        int indx = interface.first;
        cout << "INTERFACE INDEX: " + std::to_string(indx) << endl;

        auto interfaceProperties = interface.second;

        InterfaceType type = interfaceProperties->GetInterfaceType();
        cout << "Type: " << InterfaceTypeMap[type] << endl;

        CompositeMap movingDomainComposites = interfaceProperties->GetMovingDomain();
        for (auto &movingDomainComposite : movingDomainComposites)
        {
            cout << "Moving domain composite: " << movingDomainComposite.first << endl;
        }

        CompositeMap fixedDomainComposites = interfaceProperties->GetFixedDomain();
        for (auto &fixedDomainComposite : fixedDomainComposites)
        {
            cout << "Fixed domain composite: " << fixedDomainComposite.first << endl;
        }

        InterfaceEdgeMap interfaceEdgeComposites = interfaceProperties->GetInterfaceEdge();
        std::vector<unsigned> interfaceEdgeCompositeList;
        for (auto &interfaceEdgeComposite : interfaceEdgeComposites)
        {
            interfaceEdgeCompositeList.emplace_back(interfaceEdgeComposite.first);
        }
        std::string interfaceEdgeString = ParseUtils::GenerateSeqString(interfaceEdgeCompositeList);
        cout << "Interface edge composite: " << interfaceEdgeString << endl;

        if(type == eRotating)
        {
            cout << "Angular vel: " << std::to_string(interfaceProperties->GetAngularVel()) << endl;

            std::vector<NekDouble> origin(3,0);
            interfaceProperties->GetOrigin().GetCoords(origin[0], origin[1], origin[2]);
            cout << "Origin: (" << std::to_string(origin[0]) << ", " << std::to_string(origin[1]) << ", "<< std::to_string(origin[2]) << ")" << endl;

            std::vector<NekDouble> axis = interfaceProperties->GetAxis();
            cout << "Axis: (" << std::to_string(axis[0]) << ", " << std::to_string(axis[1]) << ", "<< std::to_string(axis[2]) << ")" << endl << endl;
        }
    }
}