//
// Created by Edward Laughton on 2019-03-12.
//

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/Interface.h>
#include <SpatialDomains/SegGeom.h>


using namespace Nektar;
using namespace Nektar::SpatialDomains;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session = LibUtilities::SessionReader::CreateInstance(argc, argv);
    MeshGraphSharedPtr graph = MeshGraph::Read(session);

    auto bbox = graph->GetAllQuadGeoms()[0]->GetBoundingBox();

    Array<OneD, NekDouble> min(2), max(2);

    min[0] = boost::geometry::get<boost::geometry::min_corner, 0>(bbox);
    min[1] = boost::geometry::get<boost::geometry::min_corner, 1>(bbox);

    max[0] = boost::geometry::get<boost::geometry::max_corner, 0>(bbox);
    max[1] = boost::geometry::get<boost::geometry::max_corner, 1>(bbox);

    cout << "MINIMUM CORNER IS: (" << min[0] << ", " << min[1] << ")" << endl;
    cout << "MAXIMUM CORNER IS: (" << max[0] << ", " << max[1] << ")" << endl;

}