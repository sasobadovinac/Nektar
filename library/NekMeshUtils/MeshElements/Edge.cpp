////////////////////////////////////////////////////////////////////////////////
//
//  File: Edge.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Mesh manipulation objects.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/Edge.h>
#ifdef NEKTAR_USE_MESHGEN
#include <NekMeshUtils/CADSystem/CADCurve.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>
#endif

namespace Nektar
{
namespace NekMeshUtils
{

using namespace std;

string Edge::GetXmlCurveString()
{
    std::vector<NodeSharedPtr> nodeList;

    GetCurvedNodes(nodeList);

    std::stringstream s;
    std::string str;

    // put them into a string and return
    for (int k = 0; k < nodeList.size(); ++k)
    {
        s << std::scientific << std::setprecision(8) << "    "
          << nodeList[k]->m_x << "  " << nodeList[k]->m_y << "  "
          << nodeList[k]->m_z << "    ";
    }

    return s.str();
}

SpatialDomains::SegGeomSharedPtr Edge::GetGeom(int coordDim)
{
    static boost::mutex io_mutex;
    boost::mutex::scoped_lock lock(io_mutex);

    // Create edge vertices.
    SpatialDomains::PointGeomSharedPtr p[2];
    SpatialDomains::SegGeomSharedPtr ret;

    p[0] = m_n1->GetGeom(coordDim);
    p[1] = m_n2->GetGeom(coordDim);

    // Create a curve if high-order information exists.
    if (m_edgeNodes.size() > 0)
    {
        SpatialDomains::CurveSharedPtr c =
            MemoryManager<SpatialDomains::Curve>::AllocateSharedPtr(
                m_id, m_curveType);

        c->m_points.push_back(p[0]);
        for (int i = 0; i < m_edgeNodes.size(); ++i)
        {
            c->m_points.push_back(m_edgeNodes[i]->GetGeom(coordDim));
        }
        c->m_points.push_back(p[1]);

        ret = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(
            m_id, coordDim, p, c);
    }
    else
    {
        ret = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(
            m_id, coordDim, p);
    }

    return ret;
}


void Edge::MakeOrder(int                                order,
                     SpatialDomains::GeometrySharedPtr  geom,
                     LibUtilities::PointsType           edgeType,
                     int                                coordDim,
                     int                               &id)
{
    int nPoints = order + 1;
    StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

    Array<OneD, NekDouble> edgePoints;
    LibUtilities::PointsKey edgeKey(nPoints, edgeType);
    LibUtilities::PointsManager()[edgeKey]->GetPoints(edgePoints);

    Array<OneD, Array<OneD, NekDouble> > phys(coordDim);

    for (int i = 0; i < coordDim; ++i)
    {
        phys[i] = Array<OneD, NekDouble>(xmap->GetTotPoints());
        xmap->BwdTrans(geom->GetCoeffs(i), phys[i]);
    }

    m_edgeNodes.resize(nPoints - 2);

    for (int i = 1; i < nPoints - 1; ++i)
    {
        Array<OneD, NekDouble> x(3, 0.0);
        for (int j = 0; j < coordDim; ++j)
        {
            x[j] = xmap->PhysEvaluate(edgePoints + i, phys[j]);
        }

        m_edgeNodes[i-1] = boost::shared_ptr<Node>(
            new Node(id++, x[0], x[1], x[2]));
    }

    m_curveType = edgeType;

#ifdef NEKTAR_USE_MESHGEN
    if(onSurf || onCurve)
    {
        if(onCurve)
        {
            for(int i = 0; i < m_edgeNodes.size(); i++)
            {
                Array<OneD, NekDouble> loc(3);
                loc[0] = m_edgeNodes[i]->m_x;
                loc[1] = m_edgeNodes[i]->m_y;
                loc[2] = m_edgeNodes[i]->m_z;
                NekDouble t = CADCurve->loct(loc);
                m_edgeNodes[i]->SetCADCurve(CADCurveId,CADCurve,t);

                std::vector<CADSurfSharedPtr> s = CADCurve->GetAdjSurf();
                for(int j = 0; j < s.size(); j++)
                {
                    Array<OneD, NekDouble> uv = s[j]->locuv(loc);
                    m_edgeNodes[i]->SetCADSurf(s[j]->GetId(),s[j],uv);
                }
            }
        }
        else
        {
            for(int i = 0; i < m_edgeNodes.size(); i++)
            {
                Array<OneD, NekDouble> loc(3);
                loc[0] = m_edgeNodes[i]->m_x;
                loc[1] = m_edgeNodes[i]->m_y;
                loc[2] = m_edgeNodes[i]->m_z;
                Array<OneD, NekDouble> uv = CADSurf->locuv(loc);
                m_edgeNodes[i]->SetCADSurf(CADSurfId,CADSurf,uv);
            }
        }
    }
#endif
}

}
}
