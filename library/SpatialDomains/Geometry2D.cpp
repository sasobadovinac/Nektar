////////////////////////////////////////////////////////////////////////////////
//
//  File: Geometry2D.cpp
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
//  Description: 2D geometry information
//
//
////////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <SpatialDomains/Geometry2D.h>
#include <SpatialDomains/SegGeom.h>

#include <iomanip>

using namespace std;

namespace Nektar
{
namespace SpatialDomains
{

Geometry2D::Geometry2D()
{
}

Geometry2D::Geometry2D(const int coordim, CurveSharedPtr curve)
    : Geometry(coordim), m_curve(curve)
{
    ASSERTL0(m_coordim > 1,
             "Coordinate dimension should be at least 2 for a 2D geometry");
}

Geometry2D::~Geometry2D()
{
}

int Geometry2D::v_AllLeftCheck(const Array<OneD, const NekDouble> &gloCoord)
{
    int nc = 1, d0 = m_manifold[0], d1 = m_manifold[1];
    if (0 == m_edgeNormal.size())
    {
        m_edgeNormal = Array<OneD, Array<OneD, NekDouble>>(m_verts.size());
        Array<OneD, Array<OneD, NekDouble>> x(2);
        x[0] = Array<OneD, NekDouble>(3);
        x[1] = Array<OneD, NekDouble>(3);
        m_verts[m_verts.size() - 1]->GetCoords(x[0]);
        int i0 = 1, i1 = 0;
        for (size_t i = 0; i < m_verts.size(); ++i)
        {
            i0 ^= 1;
            i1 ^= 1;
            m_verts[i]->GetCoords(x[i1]);
            if (m_edges[i]->GetXmap()->GetBasis(0)->GetNumModes() > 2)
            {
                continue;
            }
            m_edgeNormal[i]    = Array<OneD, NekDouble>(2);
            m_edgeNormal[i][0] = x[i0][d1] - x[i1][d1];
            m_edgeNormal[i][1] = x[i1][d0] - x[i0][d0];
        }
    }

    Array<OneD, NekDouble> vertex(3);
    for (size_t i = 0; i < m_verts.size(); ++i)
    {
        m_verts[i]->GetCoords(vertex);
        if (m_edgeNormal[i].size() == 0)
        {
            nc = 0; // not sure
            continue;
        }
        NekDouble value = m_edgeNormal[i][0] * (gloCoord[d0] - vertex[d0]) +
                          m_edgeNormal[i][1] * (gloCoord[d1] - vertex[d1]);
        if (value < 0)
        {
            return -1; // outside
        }
    }
    // nc: 1 (side element), 0 (maybe inside), -1 (outside)
    return nc;
}

void Geometry2D::NewtonIterationForLocCoord(
    const Array<OneD, const NekDouble> &coords, Array<OneD, NekDouble> &Lcoords)
{
    // Maximum iterations for convergence
    const int MaxIterations = 51;
    // |x-xp|^2 < EPSILON  error    tolerance
    const NekDouble Tol = 1.e-8;
    // |r,s|    > LcoordDIV stop   the search
    const NekDouble LcoordDiv = 15.0;

    NekDouble ScaledTol = fabs(m_isoParameter[0][1] * m_isoParameter[1][2] -
                               m_isoParameter[1][1] * m_isoParameter[0][2]);
    ScaledTol *= Tol * Tol;

    NekDouble xmap, ymap, F1, F2;
    NekDouble derx_1, derx_2, dery_1, dery_2, jac;

    NekDouble res;
    int cnt = 0;
    while (cnt++ < MaxIterations)
    {
        NekDouble tmp = Lcoords[0] * Lcoords[1];
        // calculate the global point corresponding to Lcoords
        xmap = m_isoParameter[0][0] + m_isoParameter[0][1] * Lcoords[0] +
               m_isoParameter[0][2] * Lcoords[1] + m_isoParameter[0][3] * tmp;
        ymap = m_isoParameter[1][0] + m_isoParameter[1][1] * Lcoords[0] +
               m_isoParameter[1][2] * Lcoords[1] + m_isoParameter[1][3] * tmp;

        F1 = coords[0] - xmap;
        F2 = coords[1] - ymap;

        res = F1 * F1 + F2 * F2;
        if (res < ScaledTol)
        {
            break;
        }

        // Interpolate derivative metric at Lcoords
        derx_1 = m_isoParameter[0][1] + m_isoParameter[0][3] * Lcoords[1];
        derx_2 = m_isoParameter[0][2] + m_isoParameter[0][3] * Lcoords[0];
        dery_1 = m_isoParameter[1][1] + m_isoParameter[1][3] * Lcoords[1];
        dery_2 = m_isoParameter[1][2] + m_isoParameter[1][3] * Lcoords[0];

        jac = 1. / (dery_2 * derx_1 - dery_1 * derx_2);

        // use analytical inverse of derivitives which are
        // also similar to those of metric factors.
        Lcoords[0] += (dery_2 * F1 - derx_2 * F2) * jac;

        Lcoords[1] += (-dery_1 * F1 + derx_1 * F2) * jac;

        if (!(std::isfinite(Lcoords[0]) && std::isfinite(Lcoords[1])) ||
            fabs(Lcoords[0]) > LcoordDiv || fabs(Lcoords[1]) > LcoordDiv)
        {
            std::ostringstream ss;
            ss << "Iteration has diverged in NewtonIterationForLocCoord in "
                  "element "
               << GetGlobalID();
            WARNINGL1(false, ss.str());
            return;
        }
    }

    if (cnt >= MaxIterations)
    {
        std::ostringstream ss;

        ss << "Reached MaxIterations (" << MaxIterations
           << ") in Newton iteration ";

        WARNINGL1(cnt < MaxIterations, ss.str());
    }
}

void Geometry2D::NewtonIterationForLocCoord(
    const Array<OneD, const NekDouble> &coords,
    const Array<OneD, const NekDouble> &ptsx,
    const Array<OneD, const NekDouble> &ptsy, Array<OneD, NekDouble> &Lcoords,
    NekDouble &dist)
{
    // Maximum iterations for convergence
    const int MaxIterations = NekConstants::kNewtonIterations;
    // |x-xp|^2 < EPSILON  error    tolerance
    const NekDouble Tol = 1.e-8;
    // |r,s|    > LcoordDIV stop   the search
    const NekDouble LcoordDiv = 15.0;

    Array<OneD, const NekDouble> Jac =
        m_geomFactors->GetJac(m_xmap->GetPointsKeys());

    NekDouble ScaledTol =
        Vmath::Vsum(Jac.size(), Jac, 1) / ((NekDouble)Jac.size());
    ScaledTol *= Tol;

    NekDouble xmap, ymap, F1, F2;
    NekDouble derx_1, derx_2, dery_1, dery_2, jac;

    // save intiial guess for later reference if required.
    NekDouble init0 = Lcoords[0], init1 = Lcoords[1];

    Array<OneD, NekDouble> DxD1(ptsx.size());
    Array<OneD, NekDouble> DxD2(ptsx.size());
    Array<OneD, NekDouble> DyD1(ptsx.size());
    Array<OneD, NekDouble> DyD2(ptsx.size());

    // Ideally this will be stored in m_geomfactors
    m_xmap->PhysDeriv(ptsx, DxD1, DxD2);
    m_xmap->PhysDeriv(ptsy, DyD1, DyD2);

    int cnt = 0;
    Array<OneD, DNekMatSharedPtr> I(2);
    Array<OneD, NekDouble> eta(2);

    F1 = F2         = 2000; // Starting value of Function
    NekDouble resid = sqrt(F1 * F1 + F2 * F2);
    while (cnt++ < MaxIterations)
    {
        //  evaluate lagrange interpolant at Lcoords
        m_xmap->LocCoordToLocCollapsed(Lcoords, eta);
        I[0] = m_xmap->GetBasis(0)->GetI(eta);
        I[1] = m_xmap->GetBasis(1)->GetI(eta + 1);

        // calculate the global point `corresponding to Lcoords
        xmap = m_xmap->PhysEvaluate(I, ptsx);
        ymap = m_xmap->PhysEvaluate(I, ptsy);

        F1 = coords[0] - xmap;
        F2 = coords[1] - ymap;

        if (F1 * F1 + F2 * F2 < ScaledTol)
        {
            resid = sqrt(F1 * F1 + F2 * F2);
            break;
        }

        // Interpolate derivative metric at Lcoords
        derx_1 = m_xmap->PhysEvaluate(I, DxD1);
        derx_2 = m_xmap->PhysEvaluate(I, DxD2);
        dery_1 = m_xmap->PhysEvaluate(I, DyD1);
        dery_2 = m_xmap->PhysEvaluate(I, DyD2);

        jac = dery_2 * derx_1 - dery_1 * derx_2;

        // use analytical inverse of derivitives which are
        // also similar to those of metric factors.
        Lcoords[0] =
            Lcoords[0] +
            (dery_2 * (coords[0] - xmap) - derx_2 * (coords[1] - ymap)) / jac;

        Lcoords[1] =
            Lcoords[1] +
            (-dery_1 * (coords[0] - xmap) + derx_1 * (coords[1] - ymap)) / jac;

        if (!(std::isfinite(Lcoords[0]) && std::isfinite(Lcoords[1])))
        {
            dist = 1e16;
            std::ostringstream ss;
            ss << "nan or inf found in NewtonIterationForLocCoord in element "
               << GetGlobalID();
            WARNINGL1(false, ss.str());
            return;
        }
        if (fabs(Lcoords[0]) > LcoordDiv || fabs(Lcoords[1]) > LcoordDiv)
        {
            break; // lcoords have diverged so stop iteration
        }
    }

    m_xmap->LocCoordToLocCollapsed(Lcoords, eta);
    if (ClampLocCoords(eta, 0.))
    {
        I[0] = m_xmap->GetBasis(0)->GetI(eta);
        I[1] = m_xmap->GetBasis(1)->GetI(eta + 1);
        // calculate the global point corresponding to Lcoords
        xmap = m_xmap->PhysEvaluate(I, ptsx);
        ymap = m_xmap->PhysEvaluate(I, ptsy);
        F1   = coords[0] - xmap;
        F2   = coords[1] - ymap;
        dist = sqrt(F1 * F1 + F2 * F2);
    }
    else
    {
        dist = 0.;
    }

    if (cnt >= MaxIterations)
    {
        Array<OneD, NekDouble> collCoords(2);
        m_xmap->LocCoordToLocCollapsed(Lcoords, collCoords);

        // if coordinate is inside element dump error!
        if ((collCoords[0] >= -1.0 && collCoords[0] <= 1.0) &&
            (collCoords[1] >= -1.0 && collCoords[1] <= 1.0))
        {
            std::ostringstream ss;

            ss << "Reached MaxIterations (" << MaxIterations
               << ") in Newton iteration ";
            ss << "Init value (" << setprecision(4) << init0 << "," << init1
               << ","
               << ") ";
            ss << "Fin  value (" << Lcoords[0] << "," << Lcoords[1] << ","
               << ") ";
            ss << "Resid = " << resid << " Tolerance = " << sqrt(ScaledTol);

            WARNINGL1(cnt < MaxIterations, ss.str());
        }
    }
}

NekDouble Geometry2D::v_GetLocCoords(const Array<OneD, const NekDouble> &coords,
                                     Array<OneD, NekDouble> &Lcoords)
{
    NekDouble dist = std::numeric_limits<double>::max();
    Array<OneD, NekDouble> tmpcoords(2);
    tmpcoords[0] = coords[m_manifold[0]];
    tmpcoords[1] = coords[m_manifold[1]];
    if (GetMetricInfo()->GetGtype() == eRegular)
    {
        tmpcoords[0] -= m_isoParameter[0][0];
        tmpcoords[1] -= m_isoParameter[1][0];
        Lcoords[0] = m_invIsoParam[0][0] * tmpcoords[0] +
                     m_invIsoParam[0][1] * tmpcoords[1];
        Lcoords[1] = m_invIsoParam[1][0] * tmpcoords[0] +
                     m_invIsoParam[1][1] * tmpcoords[1];
    }
    else if (m_straightEdge)
    {
        ClampLocCoords(Lcoords, 0.);
        NewtonIterationForLocCoord(tmpcoords, Lcoords);
    }
    else if (GetMetricInfo()->GetGtype() == eDeformed)
    {
        v_FillGeom();
        // Determine nearest point of coords  to values in m_xmap
        int npts = m_xmap->GetTotPoints();
        Array<OneD, NekDouble> ptsx(npts), ptsy(npts);
        Array<OneD, NekDouble> tmpx(npts), tmpy(npts);

        // Determine 3D manifold orientation
        m_xmap->BwdTrans(m_coeffs[m_manifold[0]], ptsx);
        m_xmap->BwdTrans(m_coeffs[m_manifold[1]], ptsy);

        Array<OneD, NekDouble> eta(2, 0.);
        m_xmap->LocCoordToLocCollapsed(Lcoords, eta);
        ClampLocCoords(eta, 0.);

        m_xmap->LocCollapsedToLocCoord(eta, Lcoords);

        // Perform newton iteration to find local coordinates
        NewtonIterationForLocCoord(tmpcoords, ptsx, ptsy, Lcoords, dist);
    }
    return dist;
}

int Geometry2D::v_GetNumVerts() const
{
    return m_verts.size();
}

int Geometry2D::v_GetNumEdges() const
{
    return m_edges.size();
}

PointGeomSharedPtr Geometry2D::v_GetVertex(int i) const
{
    ASSERTL2(i >= 0 && i < m_verts.size(), "Index out of range");
    return m_verts[i];
}

Geometry1DSharedPtr Geometry2D::v_GetEdge(int i) const
{
    ASSERTL2(i >= 0 && i < m_edges.size(), "Index out of range");
    return m_edges[i];
}

StdRegions::Orientation Geometry2D::v_GetEorient(const int i) const
{
    ASSERTL2(i >= 0 && i < m_eorient.size(), "Index out of range");
    return m_eorient[i];
}

int Geometry2D::v_GetShapeDim() const
{
    return 2;
}

NekDouble Geometry2D::v_FindDistance(const Array<OneD, const NekDouble> &xs,
                                     Array<OneD, NekDouble> &xiOut)
{
    if (m_geomFactors->GetGtype() == eRegular)
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
    // If deformed edge then the inverse mapping is non-linear so need to
    // numerically solve for the local coordinate
    else if (m_geomFactors->GetGtype() == eDeformed)
    {
        // Choose starting based on closest quad
        Array<OneD, NekDouble> xi(2, 0.0), eta(2, 0.0);
        m_xmap->LocCollapsedToLocCoord(eta, xi);

        // Armijo constants:
        // https://en.wikipedia.org/wiki/Backtracking_line_search
        const NekDouble c1 = 1e-4, c2 = 0.9;

        int nq = m_xmap->GetTotPoints();

        Array<OneD, NekDouble> x(nq), y(nq), z(nq);
        m_xmap->BwdTrans(m_coeffs[0], x);
        m_xmap->BwdTrans(m_coeffs[1], y);
        m_xmap->BwdTrans(m_coeffs[2], z);

        Array<OneD, NekDouble> xderxi1(nq, 0.0), yderxi1(nq, 0.0),
            zderxi1(nq, 0.0), xderxi2(nq, 0.0), yderxi2(nq, 0.0),
            zderxi2(nq, 0.0), xderxi1xi1(nq, 0.0), yderxi1xi1(nq, 0.0),
            zderxi1xi1(nq, 0.0), xderxi1xi2(nq, 0.0), yderxi1xi2(nq, 0.0),
            zderxi1xi2(nq, 0.0), xderxi2xi1(nq, 0.0), yderxi2xi1(nq, 0.0),
            zderxi2xi1(nq, 0.0), xderxi2xi2(nq, 0.0), yderxi2xi2(nq, 0.0),
            zderxi2xi2(nq, 0.0);

        // Get first & second derivatives & partial derivatives of x,y,z values
        std::array<NekDouble, 3> xc_derxi, yc_derxi, zc_derxi;

        m_xmap->PhysDeriv(x, xderxi1, xderxi2);
        m_xmap->PhysDeriv(y, yderxi1, yderxi2);
        m_xmap->PhysDeriv(z, zderxi1, zderxi2);

        m_xmap->PhysDeriv(xderxi1, xderxi1xi1, xderxi1xi2);
        m_xmap->PhysDeriv(yderxi1, yderxi1xi1, yderxi1xi2);
        m_xmap->PhysDeriv(zderxi1, zderxi1xi1, zderxi1xi2);

        m_xmap->PhysDeriv(yderxi2, yderxi2xi1, yderxi2xi2);
        m_xmap->PhysDeriv(xderxi2, xderxi2xi1, xderxi2xi2);
        m_xmap->PhysDeriv(zderxi2, zderxi2xi1, zderxi2xi2);

        // Minimisation loop (Quasi-newton method)
        NekDouble fx_prev = std::numeric_limits<NekDouble>::max();
        for (int i = 0; i < NekConstants::kNewtonIterations; ++i)
        {
            // Compute the objective function, f(x_k) and its derivatives
            NekDouble xc = m_xmap->PhysEvaluate(xi, x, xc_derxi);
            NekDouble yc = m_xmap->PhysEvaluate(xi, y, yc_derxi);
            NekDouble zc = m_xmap->PhysEvaluate(xi, z, zc_derxi);

            NekDouble xc_derxi1xi1 = m_xmap->PhysEvaluate(xi, xderxi1xi1);
            NekDouble yc_derxi1xi1 = m_xmap->PhysEvaluate(xi, yderxi1xi1);
            NekDouble zc_derxi1xi1 = m_xmap->PhysEvaluate(xi, zderxi1xi1);

            NekDouble xc_derxi1xi2 = m_xmap->PhysEvaluate(xi, xderxi1xi2);
            NekDouble yc_derxi1xi2 = m_xmap->PhysEvaluate(xi, yderxi1xi2);
            NekDouble zc_derxi1xi2 = m_xmap->PhysEvaluate(xi, zderxi1xi2);

            NekDouble xc_derxi2xi2 = m_xmap->PhysEvaluate(xi, xderxi2xi2);
            NekDouble yc_derxi2xi2 = m_xmap->PhysEvaluate(xi, yderxi2xi2);
            NekDouble zc_derxi2xi2 = m_xmap->PhysEvaluate(xi, zderxi2xi2);

            // Objective function is the distance to the search point
            NekDouble xdiff = xc - xs[0];
            NekDouble ydiff = yc - xs[1];
            NekDouble zdiff = zc - xs[2];

            NekDouble fx = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;

            NekDouble fx_derxi1 = 2.0 * xdiff * xc_derxi[0] +
                                  2.0 * ydiff * yc_derxi[0] +
                                  2.0 * zdiff * zc_derxi[0];

            NekDouble fx_derxi2 = 2.0 * xdiff * xc_derxi[1] +
                                  2.0 * ydiff * yc_derxi[1] +
                                  2.0 * zdiff * zc_derxi[1];

            NekDouble fx_derxi1xi1 =
                2.0 * xdiff * xc_derxi1xi1 + 2.0 * xc_derxi[0] * xc_derxi[0] +
                2.0 * ydiff * yc_derxi1xi1 + 2.0 * yc_derxi[0] * yc_derxi[0] +
                2.0 * zdiff * zc_derxi1xi1 + 2.0 * zc_derxi[0] * zc_derxi[0];

            NekDouble fx_derxi1xi2 =
                2.0 * xdiff * xc_derxi1xi2 + 2.0 * xc_derxi[1] * xc_derxi[0] +
                2.0 * ydiff * yc_derxi1xi2 + 2.0 * yc_derxi[1] * yc_derxi[0] +
                2.0 * zdiff * zc_derxi1xi2 + 2.0 * zc_derxi[1] * zc_derxi[0];

            NekDouble fx_derxi2xi2 =
                2.0 * xdiff * xc_derxi2xi2 + 2.0 * xc_derxi[1] * xc_derxi[1] +
                2.0 * ydiff * yc_derxi2xi2 + 2.0 * yc_derxi[1] * yc_derxi[1] +
                2.0 * zdiff * zc_derxi2xi2 + 2.0 * zc_derxi[1] * zc_derxi[1];

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

            // Check for convergence
            if (abs(fx - fx_prev) < 1e-12)
            {
                fx_prev = fx;
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
            pk[1] = -(hessInv[0][1] * jac[0] + hessInv[1][1] * jac[1]);

            // Backtracking line search
            while (gamma > 1e-10)
            {
                Array<OneD, NekDouble> xi_pk(2);
                xi_pk[0] = xi[0] + pk[0] * gamma;
                xi_pk[1] = xi[1] + pk[1] * gamma;

                Array<OneD, NekDouble> eta_pk(2, 0.0);
                m_xmap->LocCoordToLocCollapsed(xi_pk, eta_pk);

                if (eta_pk[0] <
                        (-1 - std::numeric_limits<NekDouble>::epsilon()) ||
                    eta_pk[0] >
                        (1 + std::numeric_limits<NekDouble>::epsilon()) ||
                    eta_pk[1] <
                        (-1 - std::numeric_limits<NekDouble>::epsilon()) ||
                    eta_pk[1] > (1 + std::numeric_limits<NekDouble>::epsilon()))
                {
                    gamma /= 2.0;
                    continue;
                }

                std::array<NekDouble, 3> xc_pk_derxi, yc_pk_derxi, zc_pk_derxi;

                NekDouble xc_pk = m_xmap->PhysEvaluate(xi_pk, x, xc_pk_derxi);
                NekDouble yc_pk = m_xmap->PhysEvaluate(xi_pk, y, yc_pk_derxi);
                NekDouble zc_pk = m_xmap->PhysEvaluate(xi_pk, z, zc_pk_derxi);

                NekDouble xc_pk_diff = xc_pk - xs[0];
                NekDouble yc_pk_diff = yc_pk - xs[1];
                NekDouble zc_pk_diff = zc_pk - xs[2];

                NekDouble fx_pk = xc_pk_diff * xc_pk_diff +
                                  yc_pk_diff * yc_pk_diff +
                                  zc_pk_diff * zc_pk_diff;

                NekDouble fx_pk_derxi1 = 2.0 * xc_pk_diff * xc_pk_derxi[0] +
                                         2.0 * yc_pk_diff * yc_pk_derxi[0] +
                                         2.0 * zc_pk_diff * zc_pk_derxi[0];

                NekDouble fx_pk_derxi2 = 2.0 * xc_pk_diff * xc_pk_derxi[1] +
                                         2.0 * yc_pk_diff * yc_pk_derxi[1] +
                                         2.0 * zc_pk_diff * zc_pk_derxi[1];

                // Check Wolfe conditions using Armijo constants
                // https://en.wikipedia.org/wiki/Wolfe_conditions
                NekDouble tmp  = pk[0] * fx_derxi1 + pk[1] * fx_derxi2;
                NekDouble tmp2 = pk[0] * fx_pk_derxi1 + pk[1] * fx_pk_derxi2;
                if ((fx_pk - (fx + c1 * gamma * tmp)) <
                        std::numeric_limits<NekDouble>::epsilon() &&
                    (-tmp2 - (-c2 * tmp)) <
                        std::numeric_limits<NekDouble>::epsilon())
                {
                    conv = true;
                    break;
                }

                gamma /= 2.0;
            }

            if (!conv)
            {
                break;
            }

            xi[0] += gamma * pk[0];
            xi[1] += gamma * pk[1];
        }

        xiOut = xi;
        return sqrt(fx_prev);
    }
    else
    {
        ASSERTL0(false, "Geometry type unknown")
    }

    return -1.0;
}

void Geometry2D::v_CalculateInverseIsoParam()
{
    NekDouble Jac = m_isoParameter[0][1] * m_isoParameter[1][2] -
                    m_isoParameter[1][1] * m_isoParameter[0][2];
    Jac = 1. / Jac;
    // a12, -a02, -a11, a01
    m_invIsoParam       = Array<OneD, Array<OneD, NekDouble>>(2);
    m_invIsoParam[0]    = Array<OneD, NekDouble>(2);
    m_invIsoParam[1]    = Array<OneD, NekDouble>(2);
    m_invIsoParam[0][0] = m_isoParameter[1][2] * Jac;
    m_invIsoParam[0][1] = -m_isoParameter[0][2] * Jac;
    m_invIsoParam[1][0] = -m_isoParameter[1][1] * Jac;
    m_invIsoParam[1][1] = m_isoParameter[0][1] * Jac;
}

} // namespace SpatialDomains
} // namespace Nektar
