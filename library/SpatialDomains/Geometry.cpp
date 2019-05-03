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
    Array<OneD, NekDouble> min(3), max(3);

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
        const int nq = GetXmap()->GetTotPoints();
        Array<OneD, Array<OneD, NekDouble>> x(3);
        for (int j = 0; j < 3; ++j)
        {
            x[j] = Array<OneD, NekDouble>(nq, 0.0);
        }
        for (int j = 0; j < GetCoordim(); ++j)
        {
            GetXmap()->BwdTrans(m_coeffs[j], x[j]);
        }
        for (int j = 0; j < 3; ++j)
        {
            min[j] = x[j][0] - NekConstants::kGeomFactorsTol;
            max[j] = x[j][0] + NekConstants::kGeomFactorsTol;

            for (int i = 1; i < nq; ++i)
            {
                min[j] = (x[j][i] < min[j] ? x[j][i] : min[j]);
                max[j] = (x[j][i] > max[j] ? x[j][i] : max[j]);
            }

            // Add 10% margin to bounding box in case elements have
            // convex boundaries.
            const int len = max[j] - min[j];
            max[j] += 0.1 * len;
            min[j] -= 0.1 * len;
        }

        for (int j = 0; j < GetCoordim(); ++j)
        {
            //cout << "Dimension: " << j <<  " | Estimate min coord:" << min[j] << endl;
            //cout << "Dimension: " << j <<  "| Estimate max coord: " << max[j] << endl;

            min[j] = FindMinCoord(j, min[j]);
            //max[j] = FindMinCoord(j, max[j]);

            //cout << "Dimension: " << j <<  " | Real min coord:" << min[j] << endl << endl;
            //cout << "Dimension: " << j <<  "| Real max coord: " << max[j] << endl;
        }
    }

    BgPoint pmin(min[0], min[1], min[2]);
    BgPoint pmax(max[0], max[1], max[2]);
    m_boundingBox = BgBox(pmin, pmax);
}

NekDouble Geometry::FindMinCoord(int coordDir, NekDouble xcEst)
{
    Array<OneD, NekDouble> xi(1, 0.5);
    const NekDouble c1 = 1e-4, c2 = 0.9;

    const int nq = GetXmap()->GetTotPoints();

    Array<OneD, NekDouble> x(nq);
    m_xmap->BwdTrans(m_coeffs[coordDir], x);

    Array<OneD, NekDouble> xder(nq, 0.0), yder(nq, 0.0);
    m_xmap->PhysDeriv(x, xder);//, yder);

    Array<OneD, NekDouble> xder2(nq, 0.0), xdery(nq, 0.0);
    Array<OneD, NekDouble> yder2(nq, 0.0), yderx(nq, 0.0);
    m_xmap->PhysDeriv(xder, xder2);//, xdery);
    //m_xmap->PhysDeriv(yder, yderx, yder2);

    // cout << "Shape type: " << LibUtilities::ShapeTypeMap[m_shapeType] << endl;
    // cout << "Shape dim: " << GetShapeDim() << endl;
    // cout << "xder: " << xder[0] << endl;
    // cout << "yder: " << yder[0] << endl;
    // cout << "xder2: " << xder2[0] << endl;
    // cout << "xdery: " << xdery[0] << endl;
    // cout << "yderx: " << yderx[0] << endl;
    // cout << "yder2: " << yder2[0] << endl;

    bool opt_succeed = false;

    NekDouble xc_prev = std::numeric_limits<NekDouble>::max();

    for (int i = 0; i < 100; ++i)
    {

        // Compute f(x_k) and its derivatives
        NekDouble xc       = m_xmap->PhysEvaluate(xi, x);
        NekDouble xc_derx  = m_xmap->PhysEvaluate(xi, xder);
        //NekDouble xc_dery  = m_xmap->PhysEvaluate(xi, yder);
        NekDouble xc_derxx = m_xmap->PhysEvaluate(xi, xder2);
        //NekDouble xc_derxy = m_xmap->PhysEvaluate(xi, xdery);
        //NekDouble xc_deryx = m_xmap->PhysEvaluate(xi, yderx);
        //NekDouble xc_deryy = m_xmap->PhysEvaluate(xi, yder2);

        cout << "Iteration: " << i << "| xc: " << xc << endl;
        // cout << "xc_derx: " << xc_derx << endl;
        // cout << "xc_dery: " << xc_dery << endl;
        // cout << "xc_derxx: " << xc_derxx << endl;
        // cout << "xc_deryy: " << xc_derxy << endl;
        // cout << "xc_deryx: " << xc_deryx << endl;
        // cout << "xc_deryy: " << xc_deryy << endl;

        // Check for convergence
        if (abs(xc - xc_prev) < 1e-16)
        {
            opt_succeed = true;
            xc_prev     = xc;
            break;
        }
        else
        {
            xc_prev = xc;
        }

        NekDouble gamma = 1.0;
        bool conv       = false;

        // Search direction: quasi-Newton, H^{-1}f . \nabla f
        NekDouble pk[1];
        NekDouble det;
        if (GetShapeDim() == 1)
        {
            pk[0] = -xc_derx / xc_derxx;
            //pk[1] = 0;
        }
        /*
        else if (GetShapeDim() == 2)
        {
            det   = xc_derxx * xc_deryy - xc_derxy * xc_deryx;
            pk[0] = (1 / det) * (xc_deryy * xc_derx - xc_derxy * xc_dery);
            pk[1] = (1 / det) * (xc_derxx * xc_dery - xc_deryx * xc_derx);
        }
        else
        {
            ASSERTL0(false, "Shape dimension is not supported")
        }
         */

        // Backtracking line search
        while (gamma > 1e-10)
        {
            Array<OneD, NekDouble> xi_pk(1);
            xi_pk[0] = xi[0] + pk[0] * gamma;
            //xi_pk[1] = xi[1] + pk[1] * gamma;

            if (xi_pk[0] < -1.0 || xi_pk[0] > 1.0)// || xi_pk[1] < -1.0 ||
                //xi_pk[1] > 1.0)
            {
                gamma /= 2.0;
                continue;
            }

            NekDouble xc_pk      = m_xmap->PhysEvaluate(xi_pk, x);
            NekDouble xc_pk_derx = m_xmap->PhysEvaluate(xi_pk, xder);
            //NekDouble xc_pk_dery = m_xmap->PhysEvaluate(xi_pk, yder);

            // Check Wolfe conditions
            //cout << "Condition 1: " << (xc_pk - (xc + c1 * gamma * pk[0] * xc_derx)) << endl;
            //cout << "Condition 2: " << (-pk[0] * xc_pk_derx + c2 * pk[0] * xc_derx) << endl;
            //if ((xc_pk - (xc + c1 * gamma * (pk[0] * xc_derx + pk[1] * xc_dery)) < std::numeric_limits<NekDouble>::epsilon()) &&
            //    (-(pk[0] * xc_pk_derx + pk[1] * xc_pk_dery) + c2 * (pk[0] * xc_derx + pk[1] * xc_dery) < std::numeric_limits<NekDouble>::epsilon()))
            if ((xc_pk - (xc + c1 * gamma * pk[0] * xc_derx) < std::numeric_limits<NekDouble>::epsilon()) &&
                (-pk[0] * xc_pk_derx + c2 * pk[0] * xc_derx) < std::numeric_limits<NekDouble>::epsilon())
            {
                conv = true;
                break;
            }

            gamma /= 2.0;
        }

        if (!conv)
        {
            opt_succeed = false;
            break;
        }

        xi[0] += gamma * pk[0];
        //xi[1] += gamma * pk[1];
    }

    if (opt_succeed)
    {
        std::cout << " success " << xc_prev << std::endl;
        return xc_prev;
    }
    else
    {
        std::cout << "failed " << std::endl;
        return std::numeric_limits<NekDouble>::max();
    }
}

} // namespace SpatialDomains
} // namespace Nektar
