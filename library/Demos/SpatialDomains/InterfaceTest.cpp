//
// Created by Edward Laughton on 2019-03-12.
//

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <SpatialDomains/MeshGraph.h>


using namespace Nektar;
using namespace Nektar::SpatialDomains;
using namespace std;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session = LibUtilities::SessionReader::CreateInstance(argc, argv);
    MeshGraphSharedPtr graph = MeshGraph::Read(session);
    for (int i =0; i < graph->GetNumElements(); ++i)
    {
        auto quad = graph->GetAllQuadGeoms()[i];
        auto bbox = quad->GetBoundingBox(true);

        cout << "Global ID: " << quad->GetGlobalID() << endl;
        cout << "Min X: " << bbox.at(0) << endl;
        cout << "Min Y: " << bbox.at(1) << endl;
        cout << "Min Z: " << bbox.at(2) << endl;
        cout << "Max X: " << bbox.at(3) << endl;
        cout << "Max Y: " << bbox.at(4) << endl;
        cout << "Max Z: " << bbox.at(5) << endl << endl;
    }

}