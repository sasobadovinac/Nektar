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
    boost::ignore_unused(argc, argv);


    PointGeomSharedPtr verts[2];
    verts[0]  = MemoryManager<PointGeom>::AllocateSharedPtr(3, 0, 0, 0, 0);
    verts[1]  = MemoryManager<PointGeom>::AllocateSharedPtr(3, 1, 1, 1, 1);

    // Create curve
    auto curve = MemoryManager<Curve>::AllocateSharedPtr(0, LibUtilities::ePolyEvenlySpaced);
    curve->m_points.emplace_back(verts[0]);
    auto curveVert  = MemoryManager<PointGeom>::AllocateSharedPtr(3, 2, 0.33, 0.66, 0.33);
    curve->m_points.emplace_back(curveVert);
    curve->m_points.emplace_back(verts[1]);

    SegGeomSharedPtr seg = MemoryManager<SegGeom>::AllocateSharedPtr(0, 3, verts, curve);
    seg->GenGeomFactors();

    Array<OneD, NekDouble> xs(3);
    xs[0] = strtod(argv[1], nullptr);
    xs[1] = strtod(argv[2], nullptr);
    xs[2] = strtod(argv[3], nullptr);

    std::cout << "LOOKING FOR POINT: " << xs[0] << " " << xs[1] << " " << xs[2] << std::endl;
    Array<OneD, NekDouble> foundPoint;
    NekDouble dist = seg->FindDistance(xs, foundPoint);

    std::cout << "LOCAL POINT ON SEG: "<< foundPoint[0] << " | DISTANCE FROM SEG: " << dist << std::endl;
    Array<OneD, NekDouble> globalCoord(3);
    globalCoord[0] = seg->GetCoord(0, foundPoint);
    globalCoord[1] = seg->GetCoord(1, foundPoint);
    globalCoord[2] = seg->GetCoord(2, foundPoint);
    std::cout << "GLOBAL POINT ON SEG: " << globalCoord[0] << " " << globalCoord[1] << " " << globalCoord[2] << std::endl;

}