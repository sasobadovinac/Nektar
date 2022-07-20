///////////////////////////////////////////////////////////////////////////////
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
/// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <Collections/Collection.h>
#include <Collections/CollectionOptimisation.h>
#include <LocalRegions/HexExp.h>
#include <SpatialDomains/MeshGraph.h>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

namespace Nektar
{
namespace FieldStorageTests
{

BOOST_AUTO_TEST_CASE(TestFieldStorageConstructor)
{
    // TODO :-)

    // SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u,
    // 0u, -1.0, -1.0, -1.0)); SpatialDomains::PointGeomSharedPtr v1(new
    // SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
    // SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u,
    // 2u, 1.0, 1.0, -1.0)); SpatialDomains::PointGeomSharedPtr v3(new
    // SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
    // SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u,
    // 4u, -1.0, -1.0, 1.0)); SpatialDomains::PointGeomSharedPtr v5(new
    // SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
    // SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u,
    // 6u, 1.0, 1.0, 1.0)); SpatialDomains::PointGeomSharedPtr v7(new
    // SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

    // SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4,
    // v5, v6, v7);

    // Nektar::LibUtilities::PointsType quadPointsTypeDir1 =
    // Nektar::LibUtilities::eGaussLobattoLegendre;
    // Nektar::LibUtilities::BasisType basisTypeDir1 =
    // Nektar::LibUtilities::eModified_A; unsigned int numQuadPoints = 6; const
    // Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints,
    // quadPointsTypeDir1); const Nektar::LibUtilities::BasisKey
    // basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

    // Nektar::LocalRegions::HexExpSharedPtr Exp =
    //     MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
    //     basisKeyDir1, basisKeyDir1, hexGeom);

    // std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
    // CollExp.push_back(Exp);

    // LibUtilities::SessionReaderSharedPtr dummySession;
    // Collections::CollectionOptimisation colOpt(dummySession, 3,
    // Collections::eStdMat); Collections::OperatorImpMap impTypes =
    // colOpt.GetOperatorImpMap(Exp); Collections::Collection     c(CollExp,
    // impTypes); c.Initialise(Collections::eBwdTrans);

    // Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
    // Array<OneD, NekDouble> phys1(Exp->GetTotPoints());
    // Array<OneD, NekDouble> phys2(Exp->GetTotPoints());

    // Exp->BwdTrans(coeffs, phys1);
    // c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

    // double epsilon = 1.0e-8;
    // for(int i = 0; i < phys1.size(); ++i)
    // {
    //     BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
    // }
}
} // namespace FieldStorageTests
} // namespace Nektar
