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
    PointGeomSharedPtr verts[2];
    verts[0]  = MemoryManager<PointGeom>::AllocateSharedPtr(2, 0, 2, 2, 0);
    verts[1]  = MemoryManager<PointGeom>::AllocateSharedPtr(2, 0, 3, 3, 0);
    SegGeomSharedPtr seg = MemoryManager<SegGeom>::AllocateSharedPtr(0, 2, verts);

    Array<OneD, NekDouble> xs(3);
    xs[0] = 2;
    xs[1] = 3;
    xs[2] = 0;

    NekDouble foundPoint;
    NekDouble dist = seg->FindDistance(xs, foundPoint);

    cout << "LOCAL POINT: "<< foundPoint << " | DISTANCE FROM SEG: " << dist << endl;

}