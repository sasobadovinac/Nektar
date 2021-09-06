////////////////////////////////////////////////////////////////////////////////
//
//  File: TriGeom.cpp
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

#include <SpatialDomains/TriGeom.h>
#include <StdRegions/StdNodalTriExp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Interp.h>

#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/GeomFactors.h>

using namespace std;

namespace Nektar
{
namespace SpatialDomains
{

TriGeom::TriGeom()
{
    m_shapeType = LibUtilities::eTriangle;
}

TriGeom::TriGeom(const int id,
                 const SegGeomSharedPtr edges[],
                 const CurveSharedPtr curve)
    : Geometry2D(edges[0]->GetVertex(0)->GetCoordim(), curve)
{
    int j;

    m_shapeType = LibUtilities::eTriangle;
    m_globalID = id;

    // Copy the edge shared pointers.
    m_edges.insert(m_edges.begin(), edges, edges + TriGeom::kNedges);
    m_eorient.resize(kNedges);

    for (j = 0; j < kNedges; ++j)
    {
        m_eorient[j] =
            SegGeom::GetEdgeOrientation(*edges[j], *edges[(j + 1) % kNedges]);
        m_verts.push_back(
            edges[j]->GetVertex(m_eorient[j] == StdRegions::eForwards ? 0 : 1));
    }

    m_eorient[2] = m_eorient[2] == StdRegions::eBackwards ?
        StdRegions::eForwards : StdRegions::eBackwards;

    m_coordim = edges[0]->GetVertex(0)->GetCoordim();
    ASSERTL0(m_coordim > 1, "Cannot call function with dim == 1");
}

TriGeom::TriGeom(const TriGeom &in)
    : Geometry2D(in)
{
    // From Geometry
    m_shapeType = in.m_shapeType;

    // From TriFaceComponent
    m_globalID = in.m_globalID;

    // From TriGeom
    m_verts = in.m_verts;
    m_edges = in.m_edges;
    for (int i = 0; i < kNedges; i++)
    {
        m_eorient[i] = in.m_eorient[i];
    }
}

TriGeom::~TriGeom()
{
}

NekDouble TriGeom::v_GetCoord(const int i,
                              const Array<OneD, const NekDouble> &Lcoord)
{
    ASSERTL1(m_state == ePtsFilled, "Geometry is not in physical space");

    Array<OneD, NekDouble> tmp(m_xmap->GetTotPoints());
    m_xmap->BwdTrans(m_coeffs[i], tmp);

    return m_xmap->PhysEvaluate(Lcoord, tmp);
}

StdRegions::Orientation TriGeom::GetFaceOrientation(const TriGeom &face1,
                              const TriGeom &face2, bool doRot, int dir,
                              NekDouble angle, NekDouble tol)
{
    return GetFaceOrientation(face1.m_verts, face2.m_verts,
                              doRot, dir, angle, tol);
}

StdRegions::Orientation TriGeom::GetFaceOrientation(
              const PointGeomVector &face1, const PointGeomVector &face2,
              bool doRot, int dir, NekDouble angle, NekDouble tol)
{
    int i, j, vmap[3] = {-1, -1, -1};

    if(doRot)
    {
        PointGeom rotPt;

        for (i = 0; i < 3; ++i)
        {
            rotPt.Rotate((*face1[i]), dir, angle);
            for (j = 0; j < 3; ++j)
            {
                if (rotPt.dist(*face2[j]) < tol)
                {
                    vmap[j] = i;
                    break;
                }
            }
        }
    }
    else
    {

        NekDouble x, y, z, x1, y1, z1, cx = 0.0, cy = 0.0, cz = 0.0;

        // For periodic faces, we calculate the vector between the centre
        // points of the two faces. (For connected faces this will be
        // zero). We can then use this to determine alignment later in the
        // algorithm.
        for (i = 0; i < 3; ++i)
        {
            cx += (*face2[i])(0) - (*face1[i])(0);
            cy += (*face2[i])(1) - (*face1[i])(1);
            cz += (*face2[i])(2) - (*face1[i])(2);
        }
        cx /= 3;
        cy /= 3;
        cz /= 3;

        // Now construct a mapping which takes us from the vertices of one
        // face to the other. That is, vertex j of face2 corresponds to
        // vertex vmap[j] of face1.
        for (i = 0; i < 3; ++i)
        {
            x = (*face1[i])(0);
            y = (*face1[i])(1);
            z = (*face1[i])(2);
            for (j = 0; j < 3; ++j)
            {
                x1 = (*face2[j])(0) - cx;
                y1 = (*face2[j])(1) - cy;
                z1 = (*face2[j])(2) - cz;
                if (sqrt((x1 - x) * (x1 - x) + (y1 - y) * (y1 - y) +
                         (z1 - z) * (z1 - z)) < 1e-8)
                {
                    vmap[j] = i;
                break;
                }
            }
        }
    }

    if (vmap[1] == (vmap[0] + 1) % 3)
    {
        switch (vmap[0])
        {
        case 0:
            return StdRegions::eDir1FwdDir1_Dir2FwdDir2;
            break;
        case 1:
            return StdRegions::eDir1FwdDir2_Dir2BwdDir1;
            break;
        case 2:
            return StdRegions::eDir1BwdDir1_Dir2BwdDir2;
            break;
        }
    }
    else
    {
        switch (vmap[0])
        {
        case 0:
            return StdRegions::eDir1FwdDir2_Dir2FwdDir1;
            break;
        case 1:
            return StdRegions::eDir1BwdDir1_Dir2FwdDir2;
            break;
        case 2:
            return StdRegions::eDir1BwdDir2_Dir2BwdDir1;
            break;
        }
    }

    ASSERTL0(false, "Unable to determine triangle orientation");
    return StdRegions::eNoOrientation;
}

void TriGeom::v_GenGeomFactors()
{
    if(!m_setupState)
    {
        TriGeom::v_Setup();
    }

    if (m_geomFactorsState != ePtsFilled)
    {
        GeomType Gtype = eRegular;

        TriGeom::v_FillGeom();

        // check to see if expansions are linear
        for (int i = 0; i < m_coordim; ++i)
        {
            if (m_xmap->GetBasisNumModes(0) != 2 ||
                m_xmap->GetBasisNumModes(1) != 2)
            {
                Gtype = eDeformed;
            }
        }

        m_geomFactors = MemoryManager<GeomFactors>::AllocateSharedPtr(
            Gtype, m_coordim, m_xmap, m_coeffs);

        m_geomFactorsState = ePtsFilled;
    }
}

/**
 * Note verts and edges are listed according to anticlockwise
 * convention but points in _coeffs have to be in array format from
 * left to right.
 */
void TriGeom::v_FillGeom()
{
    // check to see if geometry structure is already filled
    if (m_state == ePtsFilled)
    {
        return;
    }

    int i, j, k;
    int nEdgeCoeffs = m_xmap->GetTraceNcoeffs(0);

    if (m_curve)
    {
        int pdim = LibUtilities::PointsManager()[LibUtilities::PointsKey(
                                                     2, m_curve->m_ptype)]
                       ->GetPointsDim();

        // Deal with 2D points type separately
        // (e.g. electrostatic or Fekete points) to 1D tensor
        // product.
        if (pdim == 2)
        {
            int N = m_curve->m_points.size();
            int nEdgePts =
                (-1 + (int)sqrt(static_cast<NekDouble>(8 * N + 1))) / 2;

            ASSERTL0(nEdgePts * (nEdgePts + 1) / 2 == N,
                     "NUMPOINTS should be a triangle number for"
                     " triangle curved face " +
                         boost::lexical_cast<string>(m_globalID));

            // Sanity check 1: are curved vertices consistent with
            // triangle vertices?
            for (i = 0; i < 3; ++i)
            {
                NekDouble dist = m_verts[i]->dist(*(m_curve->m_points[i]));
                if (dist > NekConstants::kVertexTheSameDouble)
                {
                    std::stringstream ss;
                    ss << "Curved vertex " << i << " of triangle " << m_globalID
                       << " is separated from expansion vertex by"
                       << " more than " << NekConstants::kVertexTheSameDouble
                       << " (dist = " << dist << ")";
                    NEKERROR(ErrorUtil::ewarning, ss.str().c_str());
                }
            }

            // Sanity check 2: are curved edges from the face curvature
            // consistent with curved edges?
            for (i = 0; i < kNedges; ++i)
            {
                CurveSharedPtr edgeCurve = m_edges[i]->GetCurve();

                ASSERTL0(edgeCurve->m_points.size() == nEdgePts,
                         "Number of edge points does not correspond "
                         "to number of face points in triangle " +
                             boost::lexical_cast<string>(m_globalID));

                const int offset = 3 + i * (nEdgePts - 2);
                NekDouble maxDist = 0.0;

                // Account for different ordering of nodal coordinates
                // vs. Cartesian ordering of element.
                StdRegions::Orientation orient = m_eorient[i];

                if (i == 2)
                {
                    orient = orient == StdRegions::eForwards ?
                        StdRegions::eBackwards : StdRegions::eForwards;
                }

                if (orient == StdRegions::eForwards)
                {
                    for (j = 0; j < nEdgePts - 2; ++j)
                    {
                        NekDouble dist = m_curve->m_points[offset + j]->dist(
                            *(edgeCurve->m_points[j + 1]));
                        maxDist = dist > maxDist ? dist : maxDist;
                    }
                }
                else
                {
                    for (j = 0; j < nEdgePts - 2; ++j)
                    {
                        NekDouble dist = m_curve->m_points[offset + j]->dist(
                            *(edgeCurve->m_points[nEdgePts - 2 - j]));
                        maxDist = dist > maxDist ? dist : maxDist;
                    }
                }

                if (maxDist > NekConstants::kVertexTheSameDouble)
                {
                    std::stringstream ss;
                    ss << "Curved edge " << i << " of triangle " << m_globalID
                       << " has a point separated from edge interior"
                       << " points by more than "
                       << NekConstants::kVertexTheSameDouble
                       << " (maxdist = " << maxDist << ")";
                    NEKERROR(ErrorUtil::ewarning, ss.str().c_str());
                }
            }

            const LibUtilities::PointsKey P0(
                nEdgePts, LibUtilities::eGaussLobattoLegendre);
            const LibUtilities::PointsKey P1(
                nEdgePts, LibUtilities::eGaussRadauMAlpha1Beta0);
            const LibUtilities::BasisKey T0(
                LibUtilities::eOrtho_A, nEdgePts, P0);
            const LibUtilities::BasisKey T1(
                LibUtilities::eOrtho_B, nEdgePts, P1);
            Array<OneD, NekDouble> phys(
                max(nEdgePts * nEdgePts, m_xmap->GetTotPoints()));
            Array<OneD, NekDouble> tmp(nEdgePts * nEdgePts);

            for (i = 0; i < m_coordim; ++i)
            {
                // Create a StdNodalTriExp.
                StdRegions::StdNodalTriExpSharedPtr t =
                    MemoryManager<StdRegions::StdNodalTriExp>::
                        AllocateSharedPtr(T0, T1, m_curve->m_ptype);

                for (j = 0; j < N; ++j)
                {
                    phys[j] = (m_curve->m_points[j]->GetPtr())[i];
                }

                t->BwdTrans(phys, tmp);

                // Interpolate points to standard region.
                LibUtilities::Interp2D(P0,
                                       P1,
                                       tmp,
                                       m_xmap->GetBasis(0)->GetPointsKey(),
                                       m_xmap->GetBasis(1)->GetPointsKey(),
                                       phys);

                // Forwards transform to get coefficient space.
                m_xmap->FwdTrans(phys, m_coeffs[i]);
            }
        }
        else if (pdim == 1)
        {
            int npts = m_curve->m_points.size();
            int nEdgePts = (int)sqrt(static_cast<NekDouble>(npts));
            Array<OneD, NekDouble> tmp(npts);
            Array<OneD, NekDouble> phys(m_xmap->GetTotPoints());
            LibUtilities::PointsKey curveKey(nEdgePts, m_curve->m_ptype);

            // Sanity checks:
            // - Curved faces should have square number of points;
            // - Each edge should have sqrt(npts) points.
            ASSERTL0(nEdgePts * nEdgePts == npts,
                     "NUMPOINTS should be a square number for"
                     " triangle " +
                         boost::lexical_cast<string>(m_globalID));

            for (i = 0; i < kNedges; ++i)
            {
                ASSERTL0(m_edges[i]->GetXmap()->GetNcoeffs() == nEdgePts,
                         "Number of edge points does not correspond to "
                         "number of face points in triangle " +
                             boost::lexical_cast<string>(m_globalID));
            }

            for (i = 0; i < m_coordim; ++i)
            {
                for (j = 0; j < npts; ++j)
                {
                    tmp[j] = (m_curve->m_points[j]->GetPtr())[i];
                }

                // Interpolate curve points to standard triangle
                // points.
                LibUtilities::Interp2D(curveKey,
                                       curveKey,
                                       tmp,
                                       m_xmap->GetBasis(0)->GetPointsKey(),
                                       m_xmap->GetBasis(1)->GetPointsKey(),
                                       phys);

                // Forwards transform to get coefficient space.
                m_xmap->FwdTrans(phys, m_coeffs[i]);
            }
        }
        else
        {
            ASSERTL0(false,
                     "Only 1D/2D points distributions "
                     "supported.");
        }
    }

    Array<OneD, unsigned int> mapArray(nEdgeCoeffs);
    Array<OneD, int> signArray(nEdgeCoeffs);

    for (i = 0; i < kNedges; i++)
    {
        m_edges[i]->FillGeom();
        m_xmap->GetTraceToElementMap(i,  mapArray, signArray, m_eorient[i]);

        nEdgeCoeffs = m_edges[i]->GetXmap()->GetNcoeffs();

        for (j = 0; j < m_coordim; j++)
        {
            for (k = 0; k < nEdgeCoeffs; k++)
            {
                m_coeffs[j][mapArray[k]] =
                    signArray[k] * m_edges[i]->GetCoeffs(j)[k];
            }
        }
    }

    m_state = ePtsFilled;
}

int TriGeom::v_GetDir(const int i, const int j) const
{
    boost::ignore_unused(j); // required in 3D shapes

    return i == 0 ? 0:1;
}

void TriGeom::v_Reset(CurveMap &curvedEdges, CurveMap &curvedFaces)
{
    Geometry::v_Reset(curvedEdges, curvedFaces);
    CurveMap::iterator it = curvedFaces.find(m_globalID);

    if (it != curvedFaces.end())
    {
        m_curve = it->second;
    }

    for (int i = 0; i < 3; ++i)
    {
        m_edges[i]->Reset(curvedEdges, curvedFaces);
    }

    SetUpXmap();
    SetUpCoeffs(m_xmap->GetNcoeffs());
}

void TriGeom::v_Setup()
{
    if(!m_setupState)
    {
        for (int i = 0; i < 3; ++i)
        {
            m_edges[i]->Setup();
        }
        SetUpXmap();
        SetUpCoeffs(m_xmap->GetNcoeffs());
        m_setupState = true;
    }
}

void TriGeom::SetUpXmap()
{
    int order0 = m_edges[0]->GetXmap()->GetBasis(0)->GetNumModes();
    int order1 = max(order0,
                     max(m_edges[1]->GetXmap()->GetBasis(0)->GetNumModes(),
                         m_edges[2]->GetXmap()->GetBasis(0)->GetNumModes()));

    const LibUtilities::BasisKey B0(
        LibUtilities::eModified_A,
        order0,
        LibUtilities::PointsKey(order0+1, LibUtilities::eGaussLobattoLegendre));
    const LibUtilities::BasisKey B1(
        LibUtilities::eModified_B,
        order1,
        LibUtilities::PointsKey(order1,
                                LibUtilities::eGaussRadauMAlpha1Beta0));

    m_xmap = MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(B0, B1);
}

NekDouble TriGeom::v_FindDistance(const Array<OneD, const NekDouble> &xs,
                                   Array<OneD, NekDouble> &xiOut)
{
    // Debug to print objective function values
    if (false)
    {
        std::cout << "Looking for point: " << xs[0] << ", " << xs[1] << std::endl;

        // triangle verts
        Array<OneD, NekDouble> xt(3), yt(3), zt(3);
        for (int i = 0; i < 3; ++i)
        {
            m_verts[i]->GetCoords(xt[i], yt[i], zt[i]);
        }

        std::cout << "In triangle ID " << GetGlobalID() << ": "  << xt[0] << ", " << yt[0] << ", " << zt[0] << " -> "
                  << xt[1] << ", " << yt[1] << ", " << zt[1] << " -> "
                  << xt[2] << ", " << yt[2] << ", " << zt[2] << std::endl;


        // Print curve point locations
        ofstream file_xcurve;
        ofstream file_ycurve;
        file_xcurve.open("xcurve.txt");
        file_ycurve.open("ycurve.txt");
        std::cout << "Curve points = " << std::endl;
        for (auto & m_point : m_curve->m_points)
        {
            file_xcurve << m_point->x() << std::endl;
            file_ycurve << m_point->y() << std::endl;

            Array<OneD, NekDouble> gloCoord(2), locCoord(2);
            gloCoord[0] = m_point->x();
            gloCoord[1] = m_point->y();
            GetLocCoords(gloCoord, locCoord);
            std::cout << locCoord[0] << " " << locCoord[1] << std::endl;
        }
        std::cout << std::endl;
        file_xcurve.close();
        file_ycurve.close();


        Array<OneD, NekDouble> xi(2, 0.0);
        int nq = m_xmap->GetTotPoints();

        Array<OneD, NekDouble> x(nq), y(nq);
        m_xmap->BwdTrans(m_coeffs[0], x);
        m_xmap->BwdTrans(m_coeffs[1], y);

        Array<OneD, NekDouble> xderxi1(nq, 0.0), yderxi1(nq, 0.0),
            xderxi2(nq, 0.0), yderxi2(nq, 0.0),
            xderxi1xi1(nq, 0.0), yderxi1xi1(nq, 0.0),
            xderxi1xi2(nq, 0.0), yderxi1xi2(nq, 0.0),
            xderxi2xi1(nq, 0.0), yderxi2xi1(nq, 0.0),
            xderxi2xi2(nq, 0.0), yderxi2xi2(nq, 0.0);

        m_xmap->PhysDeriv(x, xderxi1, xderxi2);
        m_xmap->PhysDeriv(y, yderxi1, yderxi2);

        m_xmap->PhysDeriv(xderxi1, xderxi1xi1, xderxi1xi2);
        m_xmap->PhysDeriv(yderxi1, yderxi1xi1, yderxi1xi2);

        m_xmap->PhysDeriv(yderxi2, yderxi2xi1, yderxi2xi2);
        m_xmap->PhysDeriv(xderxi2, xderxi2xi1, xderxi2xi2);

        // Open files
        ofstream file_xc;
        ofstream file_yc;
        ofstream file_xi_x;
        ofstream file_xi_y;
        ofstream file_fx;
        ofstream file_fx_derxi1;
        ofstream file_fx_derxi2;
        ofstream file_fx_derxi1xi1;
        ofstream file_fx_derxi1xi2;
        ofstream file_fx_derxi2xi2;

        file_xc.open("xc.txt");
        file_yc.open("yc.txt");
        file_xi_x.open("xi_x.txt");
        file_xi_y.open("xi_y.txt");
        file_fx.open("fx.txt");
        file_fx_derxi1.open("fx_derxi1.txt");
        file_fx_derxi2.open("fx_derxi2.txt");
        file_fx_derxi1xi1.open("fx_derxi1xi1.txt");
        file_fx_derxi1xi2.open("fx_derxi1xi2.txt");
        file_fx_derxi2xi2.open("fx_derxi2xi2.txt");

        int n = 20;
        for (int i = 0; i < n + 1; ++i)
        {
            xi[0] = -1.0 + 2.0 * i / n;
            for (int j = 0; j < n; ++j)
            {
                xi[1] = -1.0 + 2.0 * j / n;

                file_xi_x << xi[0] << std::endl;
                file_xi_y << xi[1] << std::endl;

                // Compute f(x_k) and its derivatives
                NekDouble xc = m_xmap->PhysEvaluate(xi, x);
                NekDouble yc = m_xmap->PhysEvaluate(xi, y);

                NekDouble xc_derxi1 = m_xmap->PhysEvaluate(xi, xderxi1);
                NekDouble yc_derxi1 = m_xmap->PhysEvaluate(xi, yderxi1);
                NekDouble xc_derxi2 = m_xmap->PhysEvaluate(xi, xderxi2);
                NekDouble yc_derxi2 = m_xmap->PhysEvaluate(xi, yderxi2);

                NekDouble xc_derxi1xi1 = m_xmap->PhysEvaluate(xi, xderxi1xi1);
                NekDouble yc_derxi1xi1 = m_xmap->PhysEvaluate(xi, yderxi1xi1);

                NekDouble xc_derxi1xi2 = m_xmap->PhysEvaluate(xi, xderxi1xi2);
                NekDouble yc_derxi1xi2 = m_xmap->PhysEvaluate(xi, yderxi1xi2);

                NekDouble xc_derxi2xi2 = m_xmap->PhysEvaluate(xi, xderxi2xi2);
                NekDouble yc_derxi2xi2 = m_xmap->PhysEvaluate(xi, yderxi2xi2);

                file_xc << xc << std::endl;
                file_yc << yc << std::endl;

                // Objective function
                NekDouble fx =
                    (xc - xs[0]) * (xc - xs[0]) + (yc - xs[1]) * (yc - xs[1]);

                NekDouble fx_derxi1 = 2.0 * (xc - xs[0]) * xc_derxi1 +
                                      2.0 * (yc - xs[1]) * yc_derxi1;

                NekDouble fx_derxi2 = 2.0 * (xc - xs[0]) * xc_derxi2 +
                                      2.0 * (yc - xs[1]) * yc_derxi2;

                NekDouble fx_derxi1xi1 = 2.0 * (xc - xs[0]) * xc_derxi1xi1 +
                                         2.0 * xc_derxi1 * xc_derxi1 +
                                         2.0 * (yc - xs[1]) * yc_derxi1xi1 +
                                         2.0 * yc_derxi1 * yc_derxi1;

                NekDouble fx_derxi1xi2 = 2.0 * (xc - xs[0]) * xc_derxi1xi2 +
                                         2.0 * xc_derxi2 * xc_derxi1 +
                                         2.0 * (yc - xs[1]) * yc_derxi1xi2 +
                                         2.0 * yc_derxi2 * yc_derxi1;

                NekDouble fx_derxi2xi2 = 2.0 * (xc - xs[0]) * xc_derxi2xi2 +
                                         2.0 * xc_derxi2 * xc_derxi2 +
                                         2.0 * (yc - xs[1]) * yc_derxi2xi2 +
                                         2.0 * yc_derxi2 * yc_derxi2;

                // Print to file
                file_fx << fx << std::endl;
                file_fx_derxi1 << fx_derxi1 << std::endl;
                file_fx_derxi2 << fx_derxi2 << std::endl;
                file_fx_derxi1xi1 << fx_derxi1xi1 << std::endl;
                file_fx_derxi1xi2 << fx_derxi1xi2 << std::endl;
                file_fx_derxi2xi2 << fx_derxi2xi2 << std::endl;
            }
        }

        file_xc.close();
        file_yc.close();
        file_xi_x.close();
        file_xi_y.close();
        file_fx.close();
        file_fx_derxi1.close();
        file_fx_derxi2.close();
        file_fx_derxi1xi1.close();
        file_fx_derxi1xi2.close();
        file_fx_derxi2xi2.close();
    }

    if (m_geomFactors->GetGtype() == eRegular && false)
    {
        xiOut = Array<OneD, NekDouble>(2, 0.0);

        GetLocCoords(xs, xiOut);
        ClampLocCoords(xiOut);

        Array<OneD, NekDouble> gloCoord(3);
        gloCoord[0] = GetCoord(0, xiOut);
        gloCoord[1] = GetCoord(1, xiOut);
        gloCoord[2] = GetCoord(2, xiOut);

        return sqrt((xs[0] - gloCoord[0]) * (xs[0] - gloCoord[0]) +
                    (xs[1] - gloCoord[1]) * (xs[1] - gloCoord[1]) +
                    (xs[2] - gloCoord[2]) * (xs[2] - gloCoord[2]));
    }
    else if (false)
    {
        int nq = m_xmap->GetTotPoints();

        Array<OneD, NekDouble> x(nq), y(nq);
        m_xmap->BwdTrans(m_coeffs[0], x);
        m_xmap->BwdTrans(m_coeffs[1], y);

        // Pick quad point closest to xs for starter xi
        Array<OneD, NekDouble> xi(2, 0.0);
        Array<OneD, NekDouble> mxc(nq), myc(nq);
        m_xmap->GetCoords(mxc, myc);
        NekDouble fx_test = std::numeric_limits<NekDouble>::max();
        for (int i = 0; i < nq; ++i)
        {
            Array<OneD, NekDouble> testxi(2);
            testxi[0] = mxc[i];
            testxi[1] = myc[i];

            NekDouble testxc = m_xmap->PhysEvaluate(testxi, x);
            NekDouble testyc = m_xmap->PhysEvaluate(testxi, y);

            NekDouble fx =
                (testxc - xs[0]) * (testxc - xs[0]) + (testyc - xs[1]) * (testyc - xs[1]);

            if (fx < fx_test)
            {
                fx_test = fx;
                xi = testxi;
            }
        }

        Array<OneD, NekDouble> xderxi1(nq, 0.0), yderxi1(nq, 0.0),
            xderxi2(nq, 0.0), yderxi2(nq, 0.0),
            xderxi1xi1(nq, 0.0), yderxi1xi1(nq, 0.0),
            xderxi1xi2(nq, 0.0), yderxi1xi2(nq, 0.0),
            xderxi2xi1(nq, 0.0), yderxi2xi1(nq, 0.0),
            xderxi2xi2(nq, 0.0), yderxi2xi2(nq, 0.0);

        m_xmap->PhysDeriv(x, xderxi1, xderxi2);
        m_xmap->PhysDeriv(y, yderxi1, yderxi2);

        m_xmap->PhysDeriv(xderxi1, xderxi1xi1, xderxi1xi2);
        m_xmap->PhysDeriv(yderxi1, yderxi1xi1, yderxi1xi2);

        m_xmap->PhysDeriv(yderxi2, yderxi2xi1, yderxi2xi2);
        m_xmap->PhysDeriv(xderxi2, xderxi2xi1, xderxi2xi2);

        bool opt_succeed = false;
        NekDouble fx_prev = std::numeric_limits<NekDouble>::max();

        ofstream file_xc_newton;
        ofstream file_yc_newton;
        file_xc_newton.open("xc_newton.txt");
        file_yc_newton.open("yc_newton.txt");
        for (int i = 0; i < 100; ++i)
        {
            // Compute f(x_k) and its derivatives
            NekDouble xc = m_xmap->PhysEvaluate(xi, x);
            NekDouble yc = m_xmap->PhysEvaluate(xi, y);
            file_xc_newton << xc << std::endl;
            file_yc_newton << yc << std::endl;

            NekDouble xc_derxi1 = m_xmap->PhysEvaluate(xi, xderxi1);
            NekDouble yc_derxi1 = m_xmap->PhysEvaluate(xi, yderxi1);
            NekDouble xc_derxi2 = m_xmap->PhysEvaluate(xi, xderxi2);
            NekDouble yc_derxi2 = m_xmap->PhysEvaluate(xi, yderxi2);

            NekDouble xc_derxi1xi1 = m_xmap->PhysEvaluate(xi, xderxi1xi1);
            NekDouble yc_derxi1xi1 = m_xmap->PhysEvaluate(xi, yderxi1xi1);

            NekDouble xc_derxi1xi2 = m_xmap->PhysEvaluate(xi, xderxi1xi2);
            NekDouble yc_derxi1xi2 = m_xmap->PhysEvaluate(xi, yderxi1xi2);

            NekDouble xc_derxi2xi2 = m_xmap->PhysEvaluate(xi, xderxi2xi2);
            NekDouble yc_derxi2xi2 = m_xmap->PhysEvaluate(xi, yderxi2xi2);

            // Objective function
            NekDouble fx =
                (xc - xs[0]) * (xc - xs[0]) + (yc - xs[1]) * (yc - xs[1]);

            //std::cout << fx << std::endl;

            NekDouble fx_derxi1 = 2.0 * (xc - xs[0]) * xc_derxi1 +
                                  2.0 * (yc - xs[1]) * yc_derxi1;

            NekDouble fx_derxi2 = 2.0 * (xc - xs[0]) * xc_derxi2 +
                                  2.0 * (yc - xs[1]) * yc_derxi2;

            NekDouble fx_derxi1xi1 = 2.0 * (xc - xs[0]) * xc_derxi1xi1 +
                                     2.0 * xc_derxi1 * xc_derxi1 +
                                     2.0 * (yc - xs[1]) * yc_derxi1xi1 +
                                     2.0 * yc_derxi1 * yc_derxi1;

            NekDouble fx_derxi1xi2 = 2.0 * (xc - xs[0]) * xc_derxi1xi2 +
                                     2.0 * xc_derxi2 * xc_derxi1 +
                                     2.0 * (yc - xs[1]) * yc_derxi1xi2 +
                                     2.0 * yc_derxi2 * yc_derxi1;

            NekDouble fx_derxi2xi2 = 2.0 * (xc - xs[0]) * xc_derxi2xi2 +
                                     2.0 * xc_derxi2 * xc_derxi2 +
                                     2.0 * (yc - xs[1]) * yc_derxi2xi2 +
                                     2.0 * yc_derxi2 * yc_derxi2;

            // Jacobian
            NekDouble jac[2];
            jac[0] = fx_derxi1;
            jac[1] = fx_derxi2;

            // Inverse of 2x2 hessian
            NekDouble hessInv[2][2];
            NekDouble det =
                1 / (fx_derxi1xi1 * fx_derxi2xi2 - fx_derxi1xi2 * fx_derxi1xi2);
            hessInv[0][0] = det * fx_derxi2xi2;
            hessInv[0][1] = det * -fx_derxi1xi2;
            hessInv[1][0] = det * -fx_derxi1xi2;
            hessInv[1][1] = det * fx_derxi1xi1;

            xi[0] = xi[0] - (hessInv[0][0] * jac[0] + hessInv[1][0] * jac[1]);
            xi[1] = xi[1] - (hessInv[0][1] * jac[1] + hessInv[1][1] * jac[1]);

            if (xi[0] < -1.1 || xi[0] > 1.1 ||
                xi[1] < -1.1 || xi[1] > 1.1)
            {
                fx_prev = fx;
                continue;
            }

            if (abs(fx - fx_prev) < 1e-16)
            {
                opt_succeed = true;
                fx_prev     = fx;
                break;
            }
            else
            {
                fx_prev = fx;
            }
        }

        file_xc_newton.close();
        file_yc_newton.close();

        if (opt_succeed)
        {
            xiOut = xi;
            return  sqrt(fx_prev);
        }
        else
        {
            xiOut = Array<OneD, NekDouble>(2, std::numeric_limits<NekDouble>::max());
            return std::numeric_limits<NekDouble>::max();
        }
    }
    else if (m_geomFactors->GetGtype() == eDeformed)
    {
        Array<OneD, NekDouble> xi(2, 0.0);
        const NekDouble c1 = 1e-4, c2 = 0.9;

        int nq = m_xmap->GetTotPoints();

        Array<OneD, NekDouble> x(nq), y(nq);
        m_xmap->BwdTrans(m_coeffs[0], x);
        m_xmap->BwdTrans(m_coeffs[1], y);

        Array<OneD, NekDouble> xderxi1(nq, 0.0), yderxi1(nq, 0.0),
                               xderxi2(nq, 0.0), yderxi2(nq, 0.0),
                               xderxi1xi1(nq, 0.0), yderxi1xi1(nq, 0.0),
                               xderxi1xi2(nq, 0.0), yderxi1xi2(nq, 0.0),
                               xderxi2xi1(nq, 0.0), yderxi2xi1(nq, 0.0),
                               xderxi2xi2(nq, 0.0), yderxi2xi2(nq, 0.0);

        m_xmap->PhysDeriv(x, xderxi1, xderxi2);
        m_xmap->PhysDeriv(y, yderxi1, yderxi2);

        m_xmap->PhysDeriv(xderxi1, xderxi1xi1, xderxi1xi2);
        m_xmap->PhysDeriv(yderxi1, yderxi1xi1, yderxi1xi2);

        m_xmap->PhysDeriv(yderxi2, yderxi2xi1, yderxi2xi2);
        m_xmap->PhysDeriv(xderxi2, xderxi2xi1, xderxi2xi2);

        bool opt_succeed = false;

        NekDouble fx_prev = std::numeric_limits<NekDouble>::max();
        for (int i = 0; i < 100; ++i)
        {
            // Compute f(x_k) and its derivatives
            NekDouble xc = m_xmap->PhysEvaluate(xi, x);
            NekDouble yc = m_xmap->PhysEvaluate(xi, y);

            NekDouble xc_derxi1 = m_xmap->PhysEvaluate(xi, xderxi1);
            NekDouble yc_derxi1 = m_xmap->PhysEvaluate(xi, yderxi1);
            NekDouble xc_derxi2 = m_xmap->PhysEvaluate(xi, xderxi2);
            NekDouble yc_derxi2 = m_xmap->PhysEvaluate(xi, yderxi2);

            NekDouble xc_derxi1xi1 = m_xmap->PhysEvaluate(xi, xderxi1xi1);
            NekDouble yc_derxi1xi1 = m_xmap->PhysEvaluate(xi, yderxi1xi1);

            NekDouble xc_derxi1xi2 = m_xmap->PhysEvaluate(xi, xderxi1xi2);
            NekDouble yc_derxi1xi2 = m_xmap->PhysEvaluate(xi, yderxi1xi2);

            NekDouble xc_derxi2xi2 = m_xmap->PhysEvaluate(xi, xderxi2xi2);
            NekDouble yc_derxi2xi2 = m_xmap->PhysEvaluate(xi, yderxi2xi2);

            // Objective function
            NekDouble fx = (xc - xs[0]) * (xc - xs[0]) +
                           (yc - xs[1]) * (yc - xs[1]);

            NekDouble fx_derxi1 = 2.0 * (xc - xs[0]) * xc_derxi1 +
                                  2.0 * (yc - xs[1]) * yc_derxi1;

            NekDouble fx_derxi2 = 2.0 * (xc - xs[0]) * xc_derxi2 +
                                  2.0 * (yc - xs[1]) * yc_derxi2;

            NekDouble fx_derxi1xi1 = 2.0 * (xc - xs[0]) * xc_derxi1xi1 +
                                     2.0 * xc_derxi1 * xc_derxi1 +
                                     2.0 * (yc - xs[1]) * yc_derxi1xi1 +
                                     2.0 * yc_derxi1 * yc_derxi1;

            NekDouble fx_derxi1xi2 =  2.0 * (xc - xs[0]) * xc_derxi1xi2 +
                                     2.0 * xc_derxi2 * xc_derxi1 +
                                      2.0 * (yc - xs[1]) * yc_derxi1xi2 +
                                     2.0 * yc_derxi2 * yc_derxi1;

            NekDouble fx_derxi2xi2 = 2.0 * (xc - xs[0]) * xc_derxi2xi2 +
                                     2.0 * xc_derxi2 * xc_derxi2 +
                                     2.0 * (yc - xs[1]) * yc_derxi2xi2 +
                                     2.0 * yc_derxi2 * yc_derxi2;

            //std::cout << "xi[0] = " << xi[0] << ", xi[1] = " << xi[1]
            //          << ", xc = " << xc << ", yc = " << yc
            //          << ", fx = " << fx << std::endl;

            // Jacobian
            NekDouble jac[2];
            jac[0] = fx_derxi1;
            jac[1] = fx_derxi2;

            //std::cout << "jac.." << std::endl;
            //std::cout << jac[0] << " " << jac[1] << std::endl;

            // Inverse of 2x2 hessian
            NekDouble hessInv[2][2];
            //std::cout << "fx_der..." << std::endl;
            //std::cout << fx_derxi1xi1 << " " << fx_derxi1xi2 << " " << fx_derxi2xi2 << std::endl;
            NekDouble det =
                1 / (fx_derxi1xi1 * fx_derxi2xi2 - fx_derxi1xi2 * fx_derxi1xi2);
            hessInv[0][0] = det * fx_derxi2xi2;
            hessInv[0][1] = det * -fx_derxi1xi2;
            hessInv[1][0] = det * -fx_derxi1xi2;
            hessInv[1][1] = det * fx_derxi1xi1;

            //std::cout << "hess inv" << std::endl;
            //std::cout << hessInv[0][0] << " " << hessInv[0][1] << " " << hessInv[1][1] << std::endl;

            // Check for convergence
            if (abs(fx - fx_prev) < 1e-16)
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
            NekDouble pk[2];
            pk[0] = -(hessInv[0][0] * jac[0] + hessInv[1][0] * jac[1]);
            pk[1] = -(hessInv[0][1] * jac[1] + hessInv[1][1] * jac[1]);

            //std::cout << "pk..." << std::endl;
            //std::cout << pk[0] << " " << pk[1] << std::endl;

            //std::cout << "pk[0] = " << pk[0] << ", pk[1] = " << pk[1] << std::endl;
            // Backtracking line search
            while (gamma > 1e-10)
            {
                Array<OneD, NekDouble> xi_pk(2);
                xi_pk[0] = xi[0] + pk[0] * gamma;
                xi_pk[1] = xi[1] + pk[1] * gamma;

                if (xi_pk[0] < -1.0 || xi_pk[0] > 1.0 ||
                    xi_pk[1] < -1.0 || xi_pk[1] > 1.0)
                {
                    gamma /= 2.0;
                    continue;
                }

                NekDouble xc_pk = m_xmap->PhysEvaluate(xi_pk, x);
                NekDouble yc_pk = m_xmap->PhysEvaluate(xi_pk, y);

                NekDouble xc_pk_derxi1 = m_xmap->PhysEvaluate(xi_pk, xderxi1);
                NekDouble yc_pk_derxi1 = m_xmap->PhysEvaluate(xi_pk, yderxi1);
                NekDouble xc_pk_derxi2 = m_xmap->PhysEvaluate(xi_pk, xderxi2);
                NekDouble yc_pk_derxi2 = m_xmap->PhysEvaluate(xi_pk, yderxi2);

                NekDouble fx_pk = (xc_pk - xs[0]) * (xc_pk - xs[0]) +
                                  (yc_pk - xs[1]) * (yc_pk - xs[1]);

                //std::cout << "xi_pk[0] = " << xi_pk[0]
                //        << ", xi_pk[1] = " << xi_pk[1]
                //        <<" xc_pk = " << xc_pk
                //        << ", yc_pk = " << yc_pk
                //        << ", fx_pk = " << fx_pk << std::endl;

                NekDouble fx_pk_derxi1 = 2.0 * (xc_pk - xs[0]) * xc_pk_derxi1 +
                                         2.0 * (yc_pk - xs[1]) * yc_pk_derxi1;

                NekDouble fx_pk_derxi2 = 2.0 * (xc_pk - xs[0]) * xc_pk_derxi2 +
                                         2.0 * (yc_pk - xs[1]) * yc_pk_derxi2;

                // Check Wolfe conditions
                // Armijo: fx_pk =< fx + c1 * gamma * pk * fx_der;
                // Curvature (weak): -pk * fx_pk_der =< -c2 * pk * fx_der;

                // pk^T * fx_der
                NekDouble tmp = pk[0] * fx_derxi1 + pk[1] * fx_derxi2;
                // pk^T * fx_pk_der;
                NekDouble tmp2 = pk[0] * fx_pk_derxi1 + pk[1] * fx_pk_derxi2;

                //std::cout << "Armijo condition: " << fx_pk << " < " << fx + c1 * gamma * tmp << std::endl;
                //std::cout << "Curvature condition: " << tmp2 << " < " <<  c2 * tmp << std::endl;
                // Armijo condition
                if ((fx_pk  - (fx + c1 * gamma * tmp))
                        < std::numeric_limits<NekDouble>::epsilon()
                    // Curvature condition (weak) @TODO: Should this be strong condition?
                    && (-tmp2 - (-c2 * tmp))
                           < std::numeric_limits<NekDouble>::epsilon())
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
            xi[1] += gamma * pk[1];
        }

        if (opt_succeed)
        {
            xiOut = xi;
            return  sqrt(fx_prev);
        }
        else
        {
            xiOut = Array<OneD, NekDouble>(2, std::numeric_limits<NekDouble>::max());
            return std::numeric_limits<NekDouble>::max();
        }
    }
    else
    {
        ASSERTL0(false, "Geometry type unknown")
    }

    return -1.0;
}

}
}
