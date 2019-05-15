////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry.cpp
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
//  Description:  This file contains the base class implementation for the
//                Geometry class.
//
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/Geometry1D.h>
#include <SpatialDomains/Geometry2D.h>

namespace Nektar
{
namespace SpatialDomains
{

// static class property
GeomFactorsVector Geometry::m_regGeomFactorsManager;

/**
 * @brief Default constructor.
 */
Geometry::Geometry()
    : m_coordim(0), m_geomFactorsState(eNotFilled), m_state(eNotFilled),
      m_setupState(false), m_shapeType(LibUtilities::eNoShapeType),
      m_globalID(-1)
{
}

/**
 * @brief Constructor when supplied a coordinate dimension.
 */
Geometry::Geometry(const int coordim)
    : m_coordim(coordim), m_geomFactorsState(eNotFilled), m_state(eNotFilled),
      m_setupState(false), m_shapeType(LibUtilities::eNoShapeType),
      m_globalID(-1)
{
}

/**
 * @brief Default destructor.
 */
Geometry::~Geometry()
{
}

/**
 * @brief Check to see if a geometric factor has already been created that
 * contains the same regular information.
 *
 * The principle behind this is that many regular (i.e. constant Jacobian)
 * elements have identicial geometric factors. Memory may therefore be reduced
 * by storing only the unique factors.
 *
 * @param geomFactor  The GeomFactor to check.
 *
 * @return Either the cached GeomFactor or @p geomFactor.
 *
 * @todo Currently this method is disabled since the lookup is very expensive.
 */
GeomFactorsSharedPtr Geometry::ValidateRegGeomFactor(
    GeomFactorsSharedPtr geomFactor)
{
    GeomFactorsSharedPtr returnval = geomFactor;

/// \todo should this '#if 0' statement be removed?
#if 0
    bool found = false;
    if (geomFactor->GetGtype() == eRegular)
    {
        for (GeomFactorsVectorIter iter = m_regGeomFactorsManager.begin();
             iter != m_regGeomFactorsManager.end();
             ++iter)
        {
            if (**iter == *geomFactor)
            {
                returnval = *iter;
                found = true;
                break;
            }
        }

        if (!found)
        {
            m_regGeomFactorsManager.push_back(geomFactor);
            returnval = geomFactor;
        }
    }
#endif
    return returnval;
}

bool SortByGlobalId(const boost::shared_ptr<Geometry> &lhs,
                    const boost::shared_ptr<Geometry> &rhs)
{
    return lhs->GetGlobalID() < rhs->GetGlobalID();
}

bool GlobalIdEquality(const boost::shared_ptr<Geometry> &lhs,
                      const boost::shared_ptr<Geometry> &rhs)
{
    return lhs->GetGlobalID() == rhs->GetGlobalID();
}

/**
 * @brief Get the ID of vertex @p i of this object.
 */
int Geometry::GetVid(int i) const
{
    return GetVertex(i)->GetGlobalID();
}

/**
 * @brief Get the ID of edge @p i of this object.
 */
int Geometry::GetEid(int i) const
{
    return GetEdge(i)->GetGlobalID();
}

/**
 * @brief Get the ID of face @p i of this object.
 */
int Geometry::GetFid(int i) const
{
    return GetFace(i)->GetGlobalID();
}

/**
 * @copydoc Geometry::GetEdge()
 */
Geometry1DSharedPtr Geometry::v_GetEdge(int i) const
{
    NEKERROR(ErrorUtil::efatal,
             "This function is only valid for shape type geometries");
    return Geometry1DSharedPtr();
}

/**
 * @copydoc Geometry::GetFace()
 */
Geometry2DSharedPtr Geometry::v_GetFace(int i) const
{
    NEKERROR(ErrorUtil::efatal,
             "This function is only valid for shape type geometries");
    return Geometry2DSharedPtr();
}

/**
 * @copydoc Geometry::GetNumVerts()
 */
int Geometry::v_GetNumVerts() const
{
    NEKERROR(ErrorUtil::efatal,
             "This function is only valid for shape type geometries");
    return 0;
}

/**
 * @copydoc Geometry::GetEorient()
 */
StdRegions::Orientation Geometry::v_GetEorient(const int i) const
{
    NEKERROR(ErrorUtil::efatal,
             "This function is not valid for this geometry.");
    return StdRegions::eForwards;
}

/**
 * @copydoc Geometry::GetForient()
 */
StdRegions::Orientation Geometry::v_GetForient(const int i) const
{
    NEKERROR(ErrorUtil::efatal,
             "This function is not valid for this geometry.");
    return StdRegions::eFwd;
}

/**
 * @copydoc Geometry::GetNumEdges()
 */
int Geometry::v_GetNumEdges() const
{
    NEKERROR(ErrorUtil::efatal,
             "This function is only valid for shape type geometries");
    return 0;
}

/**
 * @copydoc Geometry::GetNumFaces()
 */
int Geometry::v_GetNumFaces() const
{
    NEKERROR(ErrorUtil::efatal,
             "This function is only valid for shape type geometries");
    return 0;
}

/**
 * @copydoc Geometry::GetShapeDim()
 */
int Geometry::v_GetShapeDim() const
{
    NEKERROR(ErrorUtil::efatal,
             "This function is only valid for shape type geometries");
    return 0;
}

/**
 * @copydoc Geometry::GetXmap()
 */
StdRegions::StdExpansionSharedPtr Geometry::v_GetXmap() const
{
    return m_xmap;
}

/**
 * @copydoc Geometry::ContainsPoint(
 *     const Array<OneD, const NekDouble> &, Array<OneD, NekDouble> &,
 *     NekDouble, NekDouble&)
 */
bool Geometry::v_ContainsPoint(const Array<OneD, const NekDouble> &gloCoord,
                               Array<OneD, NekDouble> &locCoord, NekDouble tol,
                               NekDouble &resid)
{
    NEKERROR(ErrorUtil::efatal,
             "This function has not been defined for this geometry");
    return false;
}

NekDouble Geometry::v_FindDistance(const Array<OneD, const NekDouble> &xs,
                                   NekDouble &xi)
{
    NEKERROR(ErrorUtil::efatal,
             "This function has not been defined for this geometry");
    return false;
}

/**
 * @copydoc Geometry::GetVertexEdgeMap()
 */
int Geometry::v_GetVertexEdgeMap(const int i, const int j) const
{
    NEKERROR(ErrorUtil::efatal,
             "This function has not been defined for this geometry");
    return 0;
}

/**
 * @copydoc Geometry::GetVertexFaceMap()
 */
int Geometry::v_GetVertexFaceMap(const int i, const int j) const
{
    NEKERROR(ErrorUtil::efatal,
             "This function has not been defined for this geometry");
    return 0;
}

/**
 * @copydoc Geometry::GetEdgeFaceMap()
 */
int Geometry::v_GetEdgeFaceMap(const int i, const int j) const
{
    NEKERROR(ErrorUtil::efatal,
             "This function has not been defined for this geometry");
    return 0;
}

/**
 * @copydoc Geometry::GetCoord()
 */
NekDouble Geometry::v_GetCoord(const int i,
                               const Array<OneD, const NekDouble> &Lcoord)
{
    NEKERROR(ErrorUtil::efatal,
             "This function is only valid for expansion type geometries");
    return 0.0;
}

/**
 * @copydoc Geometry::GetLocCoords()
 */
NekDouble Geometry::v_GetLocCoords(const Array<OneD, const NekDouble> &coords,
                                   Array<OneD, NekDouble> &Lcoords)
{
    NEKERROR(ErrorUtil::efatal,
             "This function is only valid for expansion type geometries");
    return 0.0;
}

/**
 * @copydoc Geometry::FillGeom()
 */
void Geometry::v_FillGeom()
{
    NEKERROR(ErrorUtil::efatal,
             "This function is only valid for expansion type geometries");
}

/**
 * @copydoc Geometry::Reset()
 */
void Geometry::v_Reset(CurveMap &curvedEdges, CurveMap &curvedFaces)
{
    // Reset state
    m_state            = eNotFilled;
    m_geomFactorsState = eNotFilled;

    // Junk geometric factors
    m_geomFactors = GeomFactorsSharedPtr();
}

void Geometry::v_Setup()
{
    NEKERROR(ErrorUtil::efatal,
             "This function is only valid for expansion type geometries");
}

/**
 * @brief Generates the bounding box for the element.
 *
 * For regular elements, the vertices are sufficient to define the extent of
 * the bounding box. For non-regular elements, the extremes of the quadrature
 * point coordinates are used. A 10% margin is added around this computed
 * region to account for convex hull elements where the true extent of the
 * element may extend slightly beyond the quadrature points.
 */
void Geometry::GenBoundingBox()
{
    // NekDouble minx, miny, minz, maxx, maxy, maxz;
    Array<OneD, NekDouble> min(3, 0.0), max(3, 0.0);

    if (GetGeomFactors()->GetGtype() == eRegular)
    {
        PointGeomSharedPtr p = GetVertex(0);
        Array<OneD, NekDouble> x(3, 0.0);
        p->GetCoords(x[0], x[1], x[2]);
        for (int j = 0; j < 3; ++j)
        {
            min[j] = x[j] - NekConstants::kGeomFactorsTol;
            max[j] = x[j] + NekConstants::kGeomFactorsTol;
        }
        for (int i = 1; i < GetNumVerts(); ++i)
        {
            p = GetVertex(i);
            p->GetCoords(x[0], x[1], x[2]);
            for (int j = 0; j < 3; ++j)
            {
                min[j] = (x[j] < min[j] ? x[j] : min[j]);
                max[j] = (x[j] > max[j] ? x[j] : max[j]);
            }
        }
    }
    else
    {
        for (int j = 0; j < GetCoordim(); ++j)
        {
            auto minMax = FindMinMaxCoord(j);
            min[j]      = minMax.first - NekConstants::kGeomFactorsTol;
            max[j]      = minMax.second + NekConstants::kGeomFactorsTol;
        }
    }

    BgPoint pmin(min[0], min[1], min[2]);
    BgPoint pmax(max[0], max[1], max[2]);
    m_boundingBox = BgBox(pmin, pmax);
}

std::pair<NekDouble, NekDouble> Geometry::FindMinMaxCoord(int coordDir)
{
    switch (m_shapeType)
    {
        case LibUtilities::eSegment:
        {
            std::unordered_set<NekDouble> values;

            const int nq = m_xmap->GetTotPoints();
            Array<OneD, NekDouble> x(nq, 0.0), xder(nq, 0.0), xder2(nq, 0.0);

            const int n = 3;
            Array<OneD, NekDouble> points(n + 1, 0.0);

            m_xmap->BwdTrans(m_coeffs[coordDir], x);
            m_xmap->PhysDeriv(x, xder);
            m_xmap->PhysDeriv(xder, xder2);

            for (int i = 0; i < n + 1; ++i)
            {
                Array<OneD, NekDouble> xi(1, (i * (2.0 / n) - 1.0));
                NekDouble xi_prev = xi[0];
                for (int j = 0; j < 10; ++j)
                {
                    NekDouble xc       = m_xmap->PhysEvaluate(xi, x);
                    NekDouble xc_derx  = m_xmap->PhysEvaluate(xi, xder);
                    NekDouble xc_derxx = m_xmap->PhysEvaluate(xi, xder2);

                    xi[0] = xi_prev - xc_derx / xc_derxx;

                    if ((abs(xi[0] - xi_prev) < 1e-10) && (xi[0] >= -1) &&
                        (xi[0] <= 1))
                    {
                        values.insert(xc);
                        break;
                    }

                    xi_prev = xi[0];
                }
            }

            const Array<OneD, NekDouble> minusOne(1, -1);
            const Array<OneD, NekDouble> plusOne(1, 1);
            values.insert(m_xmap->PhysEvaluate(minusOne, x));
            values.insert(m_xmap->PhysEvaluate(plusOne, x));

            const auto res =
                std::minmax_element(std::begin(values), std::end(values));
            return std::make_pair(*res.first, *res.second);
        }

        case LibUtilities::eTriangle:
        case LibUtilities::eQuadrilateral:
        {
            if (m_coordim == 2)
            {
                std::unordered_set<NekDouble> values;

                for (int edgeID = 0; edgeID < GetNumEdges(); ++edgeID)
                {
                    auto edgeMinMax =
                        GetEdge(edgeID)->FindMinMaxCoord(coordDir);
                    values.insert(edgeMinMax.first);
                    values.insert(edgeMinMax.second);
                }

                const auto res =
                    std::minmax_element(std::begin(values), std::end(values));
                return std::make_pair(*res.first, *res.second);
            }
        }

        default: // @todo change default method so it solves numerically
        {
            std::unordered_set<NekDouble> values;

            const int nq = m_xmap->GetTotPoints();
            Array<OneD, NekDouble> x(nq, 0.0);

            m_xmap->BwdTrans(m_coeffs[coordDir], x);

            for (int i = 0; i < nq; ++i)
            {
                values.insert(x[i]);
            }

            const auto res =
                std::minmax_element(std::begin(values), std::end(values));
            auto minMax = std::make_pair(*res.first, *res.second);

            // Add 10% margin to bounding box in case elements have
            // convex boundaries.
            const int len = minMax.second - minMax.first;

            minMax.first -= 0.1 * len;
            minMax.second += 0.1 * len;

            return minMax;
        }
    }
}

} // namespace SpatialDomains
} // namespace Nektar
