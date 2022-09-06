////////////////////////////////////////////////////////////////////////////////
//
//  File:  SegGeom.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/SegGeom.h>

#include <LibUtilities/Foundations/ManagerAccess.h> // for PointsManager, etc
#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdSegExp.h>

namespace Nektar
{
namespace SpatialDomains
{
SegGeom::SegGeom()
{
    m_shapeType = LibUtilities::eSegment;
}

SegGeom::SegGeom(int id, const int coordim, const PointGeomSharedPtr vertex[],
                 const CurveSharedPtr curve)
    : Geometry1D(coordim)
{
    m_shapeType = LibUtilities::eSegment;
    m_globalID  = id;
    m_state     = eNotFilled;
    m_curve     = curve;

    m_verts[0] = vertex[0];
    m_verts[1] = vertex[1];
}

SegGeom::SegGeom(const SegGeom &in) : Geometry1D(in)
{
    // From Geometry class
    m_shapeType = in.m_shapeType;

    // info from EdgeComponent class
    m_globalID = in.m_globalID;
    m_xmap     = in.m_xmap;
    SetUpCoeffs(m_xmap->GetNcoeffs());

    // info from SegGeom class
    m_coordim  = in.m_coordim;
    m_verts[0] = in.m_verts[0];
    m_verts[1] = in.m_verts[1];

    m_state = in.m_state;
}

void SegGeom::SetUpXmap()
{
    if (m_curve)
    {
        int npts = m_curve->m_points.size();
        LibUtilities::PointsKey pkey(npts + 1,
                                     LibUtilities::eGaussLobattoLegendre);
        const LibUtilities::BasisKey B(LibUtilities::eModified_A, npts, pkey);
        m_xmap = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
    }
    else
    {
        const LibUtilities::BasisKey B(
            LibUtilities::eModified_A, 2,
            LibUtilities::PointsKey(2, LibUtilities::eGaussLobattoLegendre));
        m_xmap = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
    }
}

/**
 * \brief Generate a one dimensional space segment geometry where the vert[0]
 * has the same x value and vert[1] is set to vert[0] plus the length of the
 * original segment
 **/
SegGeomSharedPtr SegGeom::GenerateOneSpaceDimGeom(void)
{
    SegGeomSharedPtr returnval = MemoryManager<SegGeom>::AllocateSharedPtr();

    // info about numbering
    returnval->m_globalID = m_globalID;

    // geometric information.
    returnval->m_coordim     = 1;
    NekDouble x0             = (*m_verts[0])[0];
    PointGeomSharedPtr vert0 = MemoryManager<PointGeom>::AllocateSharedPtr(
        1, m_verts[0]->GetGlobalID(), x0, 0.0, 0.0);
    vert0->SetGlobalID(vert0->GetGlobalID());
    returnval->m_verts[0] = vert0;

    // Get information to calculate length.
    const Array<OneD, const LibUtilities::BasisSharedPtr> base =
        m_xmap->GetBase();
    LibUtilities::PointsKeyVector v;
    v.push_back(base[0]->GetPointsKey());
    v_GenGeomFactors();

    const Array<OneD, const NekDouble> jac = m_geomFactors->GetJac(v);

    NekDouble len = 0.0;
    if (jac.size() == 1)
    {
        len = jac[0] * 2.0;
    }
    else
    {
        Array<OneD, const NekDouble> w0 = base[0]->GetW();
        len                             = 0.0;

        for (int i = 0; i < jac.size(); ++i)
        {
            len += jac[i] * w0[i];
        }
    }
    // Set up second vertex.
    PointGeomSharedPtr vert1 = MemoryManager<PointGeom>::AllocateSharedPtr(
        1, m_verts[1]->GetGlobalID(), x0 + len, 0.0, 0.0);
    vert1->SetGlobalID(vert1->GetGlobalID());

    returnval->m_verts[1] = vert1;

    // at present just use previous m_xmap[0];
    returnval->m_xmap = m_xmap;
    returnval->SetUpCoeffs(m_xmap->GetNcoeffs());
    returnval->m_state = eNotFilled;

    return returnval;
}

SegGeom::~SegGeom()
{
}

LibUtilities::ShapeType SegGeom::v_GetShapeType() const
{
    return LibUtilities::eSegment;
}

NekDouble SegGeom::v_GetCoord(const int i,
                              const Array<OneD, const NekDouble> &Lcoord)
{
    ASSERTL1(m_state == ePtsFilled, "Geometry is not in physical space");

    Array<OneD, NekDouble> tmp(m_xmap->GetTotPoints());
    m_xmap->BwdTrans(m_coeffs[i], tmp);

    return m_xmap->PhysEvaluate(Lcoord, tmp);
}

/**
 * @brief Get the orientation of @p edge1.
 *
 * If @p edge1 is connected to @p edge2 in the same direction as the points
 * comprising @p edge1 then it is forward, otherwise it is backward.
 *
 * For example, assume @p edge1 is comprised of points 1 and 2, and @p edge2 is
 * comprised of points 2 and 3, then @p edge1 is forward.
 *
 * If @p edge1 is comprised of points 2 and 1 and @p edge2 is comprised of
 * points 3 and 2, then @p edge1 is backward.
 *
 * Since both edges are passed, it does not need any information from the
 * EdgeComponent instance.
 */
StdRegions::Orientation SegGeom::GetEdgeOrientation(const SegGeom &edge1,
                                                    const SegGeom &edge2)
{
    StdRegions::Orientation returnval = StdRegions::eForwards;

    if ((*edge1.GetVertex(0) == *edge2.GetVertex(0)) ||
        (*edge1.GetVertex(0) == *edge2.GetVertex(1)))
    {
        // Backward direction.  Vertex 0 is connected to edge 2.
        returnval = StdRegions::eBackwards;
    }
    else if ((*edge1.GetVertex(1) != *edge2.GetVertex(0)) &&
             (*edge1.GetVertex(1) != *edge2.GetVertex(1)))
    {
        // Not forward either, then we have a problem.
        std::ostringstream errstrm;
        errstrm << "Connected edges do not share a vertex. Edges ";
        errstrm << edge1.GetGlobalID() << ", " << edge2.GetGlobalID();
        ASSERTL0(false, errstrm.str());
    }

    return returnval;
}

void SegGeom::v_GenGeomFactors()
{
    if (!m_setupState)
    {
        SegGeom::v_Setup();
    }

    if (m_geomFactorsState != ePtsFilled)
    {
        SpatialDomains::GeomType gType = eRegular;
        SegGeom::v_FillGeom();

        if (m_xmap->GetBasisNumModes(0) != 2)
        {
            gType = eDeformed;
        }

        m_geomFactors = MemoryManager<GeomFactors>::AllocateSharedPtr(
            gType, m_coordim, m_xmap, m_coeffs);
        m_geomFactorsState = ePtsFilled;
    }
}

void SegGeom::v_FillGeom()
{
    if (m_state != ePtsFilled)
    {
        int i;

        if (m_coordim > 0 && m_curve)
        {
            int npts = m_curve->m_points.size();
            LibUtilities::PointsKey pkey(npts + 1,
                                         LibUtilities::eGaussLobattoLegendre);
            Array<OneD, NekDouble> tmp(npts);

            if (m_verts[0]->dist(*(m_curve->m_points[0])) >
                NekConstants::kVertexTheSameDouble)
            {
                std::string err =
                    "Vertex 0 is separated from first point by more than ";
                std::stringstream strstrm;
                strstrm << NekConstants::kVertexTheSameDouble << " in edge "
                        << m_globalID;
                err += strstrm.str();
                NEKERROR(ErrorUtil::ewarning, err.c_str());
            }

            if (m_verts[1]->dist(*(m_curve->m_points[npts - 1])) >
                NekConstants::kVertexTheSameDouble)
            {
                std::string err =
                    "Vertex 1 is separated from last point by more than ";
                std::stringstream strstrm;
                strstrm << NekConstants::kVertexTheSameDouble << " in edge "
                        << m_globalID;
                err += strstrm.str();
                NEKERROR(ErrorUtil::ewarning, err.c_str());
            }

            LibUtilities::PointsKey fkey(npts, m_curve->m_ptype);
            DNekMatSharedPtr I0 =
                LibUtilities::PointsManager()[fkey]->GetI(pkey);
            NekVector<NekDouble> out(npts + 1);

            for (int i = 0; i < m_coordim; ++i)
            {
                // Load up coordinate values into tmp
                for (int j = 0; j < npts; ++j)
                {
                    tmp[j] = (m_curve->m_points[j]->GetPtr())[i];
                }

                // Interpolate to GLL points
                NekVector<NekDouble> in(npts, tmp, eWrapper);
                out = (*I0) * in;

                m_xmap->FwdTrans(out.GetPtr(), m_coeffs[i]);
            }
        }

        for (i = 0; i < m_coordim; ++i)
        {
            m_coeffs[i][0] = (*m_verts[0])[i];
            m_coeffs[i][1] = (*m_verts[1])[i];
        }

        m_state = ePtsFilled;
    }
}

void SegGeom::v_Reset(CurveMap &curvedEdges, CurveMap &curvedFaces)
{
    Geometry::v_Reset(curvedEdges, curvedFaces);
    CurveMap::iterator it = curvedEdges.find(m_globalID);

    if (it != curvedEdges.end())
    {
        m_curve = it->second;
    }

    SetUpXmap();
    SetUpCoeffs(m_xmap->GetNcoeffs());
}

void SegGeom::v_Setup()
{
    if (!m_setupState)
    {
        SetUpXmap();
        SetUpCoeffs(m_xmap->GetNcoeffs());
        m_setupState = true;
    }
}

PointGeomSharedPtr SegGeom::v_GetVertex(const int i) const
{
    PointGeomSharedPtr returnval;

    if (i >= 0 && i < kNverts)
    {
        returnval = m_verts[i];
    }

    return returnval;
}

int SegGeom::v_GetNumVerts() const
{
    return kNverts;
}

NekDouble SegGeom::v_FindDistance(const Array<OneD, const NekDouble> &xs,
                                  Array<OneD, NekDouble> &xiOut)
{
    if (m_geomFactors->GetGtype() == eRegular)
    {
        xiOut = Array<OneD, NekDouble>(1, 0.0);

        GetLocCoords(xs, xiOut);
        ClampLocCoords(xiOut);

        Array<OneD, NekDouble> gloCoord(m_coordim);
        NekDouble tmp = 0;
        for (int i = 0; i < m_coordim; ++i)
        {
            gloCoord[i] = GetCoord(i, xiOut);
            tmp += (xs[i] - gloCoord[i]) * (xs[i] - gloCoord[i]);
        }

        return sqrt(tmp);
    }
    else if (m_geomFactors->GetGtype() == eDeformed)
    {
        Array<OneD, NekDouble> xi(1, 0.0);
        const NekDouble c1 = 1e-4, c2 = 0.9;

        int dim = GetCoordim();
        int nq  = m_xmap->GetTotPoints();

        Array<OneD, Array<OneD, NekDouble>> x(dim), xder(dim), xder2(dim);
        for (int i = 0; i < dim; ++i)
        {
            x[i]     = Array<OneD, NekDouble>(nq);
            xder[i]  = Array<OneD, NekDouble>(nq);
            xder2[i] = Array<OneD, NekDouble>(nq);

            m_xmap->BwdTrans(m_coeffs[i], x[i]);
        }

        bool opt_succeed = false;

        NekDouble fx_prev = std::numeric_limits<NekDouble>::max();

        for (int i = 0; i < 100; ++i)
        {
            // Compute f(x_k) and its derivatives
            Array<OneD, NekDouble> xc(dim), xc_der(dim), xc_der2(dim);
            NekDouble fx = 0, fxp = 0, fxp2 = 0, xcDiff = 0;
            for (int j = 0; j < dim; ++j)
            {
                xc[j] = m_xmap->PhysEvaluate2ndDeriv(xi, x[j], xc_der[j],
                                                     xc_der2[j]);

                xcDiff = xc[j] - xs[j];
                fx += xcDiff * xcDiff;
                fxp += xc_der[j] * xcDiff;
                fxp2 += xc_der2[j] * xcDiff + xc_der[j] * xc_der[j];
            }

            fxp *= 2;
            fxp2 *= 2;

            // Check for convergence
            if (std::abs(fx - fx_prev) < 1e-12)
            {
                opt_succeed = true;
                fx_prev     = fx;
                break;
            }
            else
            {
                fx_prev = fx;
            }

            NekDouble gamma = 1.0;
            bool conv       = false;

            // Search direction: Newton's method
            NekDouble pk = -fxp / fxp2;

            // Backtracking line search
            while (gamma > 1e-10)
            {
                Array<OneD, NekDouble> xi_pk(1);
                xi_pk[0] = xi[0] + pk * gamma;

                if (xi_pk[0] < -1.0 || xi_pk[0] > 1.0)
                {
                    gamma /= 2.0;
                    continue;
                }

                Array<OneD, NekDouble> xc_pk(dim), xc_der_pk(dim);
                NekDouble fx_pk = 0, fxp_pk = 0, xc_pkDiff = 0;
                for (int j = 0; j < dim; ++j)
                {
                    xc_pk[j] = m_xmap->PhysEvaluate(xi_pk, x[j], xc_der_pk[j]);

                    xc_pkDiff = xc_pk[j] - xs[j];
                    fx_pk += xc_pkDiff * xc_pkDiff;
                    fxp_pk += xc_der_pk[j] * xc_pkDiff;
                }

                fxp_pk *= 2;

                // Check Wolfe conditions
                if ((fx_pk - (fx + c1 * gamma * pk * fxp)) <
                        std::numeric_limits<NekDouble>::epsilon() &&
                    (-pk * fxp_pk + c2 * pk * fxp) <
                        std::numeric_limits<NekDouble>::epsilon())
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

            xi[0] += gamma * pk;
        }

        if (opt_succeed)
        {
            xiOut = xi;
            return sqrt(fx_prev);
        }
        else
        {
            xiOut = Array<OneD, NekDouble>(
                2, std::numeric_limits<NekDouble>::max());
            return std::numeric_limits<NekDouble>::max();
        }
    }
    else
    {
        ASSERTL0(false, "Geometry type unknown")
    }

    return -1.0;
}

} // namespace SpatialDomains
} // namespace Nektar
