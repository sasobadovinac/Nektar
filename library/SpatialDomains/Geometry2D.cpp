////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry2D.cpp
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
//  Description:  2D geometry information
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

void Geometry2D::NewtonIterationForLocCoord(
    const Array<OneD, const NekDouble> &coords,
    const Array<OneD, const NekDouble> &ptsx,
    const Array<OneD, const NekDouble> &ptsy,
    Array<OneD, NekDouble> &Lcoords,
    NekDouble &dist)
{
    // Maximum iterations for convergence
    const int MaxIterations = 51;
    // |x-xp|^2 < EPSILON  error    tolerance
    const NekDouble Tol = 1.e-8;
    // |r,s|    > LcoordDIV stop   the search
    const NekDouble LcoordDiv = 15.0;

    Array<OneD, const NekDouble> Jac =
        m_geomFactors->GetJac(m_xmap->GetPointsKeys());

    NekDouble ScaledTol = Vmath::Vsum(Jac.size(), Jac, 1) /
                          ((NekDouble)Jac.size());
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

    F1 = F2 = 2000; // Starting value of Function
    NekDouble resid;
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

        if( !(std::isfinite(Lcoords[0]) && std::isfinite(Lcoords[1])) )
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
    if(ClampLocCoords(eta, 0.))
    {
        I[0] = m_xmap->GetBasis(0)->GetI(eta);
        I[1] = m_xmap->GetBasis(1)->GetI(eta + 1);
        // calculate the global point corresponding to Lcoords
        xmap = m_xmap->PhysEvaluate(I, ptsx);
        ymap = m_xmap->PhysEvaluate(I, ptsy);
        F1 = coords[0] - xmap;
        F2 = coords[1] - ymap;
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
    NekDouble dist = 0.;
    if (GetMetricInfo()->GetGtype() == eRegular)
    {
        int v2;
        if(m_shapeType == LibUtilities::eTriangle)
        {
            v2 = 2;
        }
        else if(m_shapeType == LibUtilities::eQuadrilateral)
        {
            v2 = 3;
        }
        else
        {
            v2 = 2;
            ASSERTL0(false, "unrecognized 2D element type");
        }

        NekDouble coords2 = (m_coordim == 3) ? coords[2] : 0.0;
        PointGeom r(m_coordim, 0, coords[0], coords[1], coords2);

        // Edges
        PointGeom er0, e10, e20;
        PointGeom norm, orth1, orth2;
        er0.Sub(r, *m_verts[0]);
        e10.Sub(*m_verts[1], *m_verts[0]);
        e20.Sub(*m_verts[v2], *m_verts[0]);

        // Obtain normal to element plane
        norm.Mult(e10, e20);
        // Obtain vector which are proportional to normal of e10 and e20.
        orth1.Mult(norm, e10);
        orth2.Mult(norm, e20);

        // Calculate length using L/|dv1| = (x-v0).n1/(dv1.n1) for coordiante 1
        // Then rescale to [-1,1].
        Lcoords[0] = er0.dot(orth2) / e10.dot(orth2);
        Lcoords[0] = 2. * Lcoords[0] - 1.;
        Lcoords[1] = er0.dot(orth1) / e20.dot(orth1);
        Lcoords[1] = 2. * Lcoords[1] - 1.;

        // Set distance
        Array<OneD, NekDouble> eta(2, 0.);
        m_xmap->LocCoordToLocCollapsed(Lcoords, eta);
        if(ClampLocCoords(eta, 0.))
        {
            Array<OneD, NekDouble> xi(2, 0.);
            m_xmap->LocCollapsedToLocCoord(eta, xi);
            xi[0] = (xi[0] + 1.) * 0.5; //re-scaled to ratio [0, 1]
            xi[1] = (xi[1] + 1.) * 0.5;
            for (int i = 0; i < m_coordim; ++i)
            {
                NekDouble tmp = xi[0]*e10[i] + xi[1]*e20[i] - er0[i];
                dist += tmp * tmp;
            }
            dist = sqrt(dist);
        }
    }
    else
    {
        v_FillGeom();
        // Determine nearest point of coords  to values in m_xmap
        int npts = m_xmap->GetTotPoints();
        Array<OneD, NekDouble> ptsx(npts), ptsy(npts);
        Array<OneD, NekDouble> tmpx(npts), tmpy(npts);

        // Determine 3D manifold orientation
        int idx = 0, idy = 1;
        if(m_coordim == 3)
        {
            PointGeom e01, e21, norm;
            e01.Sub(*m_verts[0], *m_verts[1]);
            e21.Sub(*m_verts[2], *m_verts[1]);
            norm.Mult(e01, e21);
            int tmpi = 0;
            double tmp = std::fabs(norm[0]);
            if(tmp < fabs(norm[1]))
            {
                tmp = fabs(norm[1]);
                tmpi = 1;
            }
            if(tmp < fabs(norm[2]))
            {
                tmpi = 2;
            }
            idx = (tmpi + 1) % 3;
            idy = (tmpi + 2) % 3;
        }
        Array<OneD, NekDouble> tmpcoords(2);
        tmpcoords[0] = coords[idx];
        tmpcoords[1] = coords[idy];

        m_xmap->BwdTrans(m_coeffs[idx], ptsx);
        m_xmap->BwdTrans(m_coeffs[idy], ptsy);

        const Array<OneD, const NekDouble> za = m_xmap->GetPoints(0);
        const Array<OneD, const NekDouble> zb = m_xmap->GetPoints(1);

        // guess the first local coords based on nearest point
        Vmath::Sadd(npts, -tmpcoords[0], ptsx, 1, tmpx, 1);
        Vmath::Sadd(npts, -tmpcoords[1], ptsy, 1, tmpy, 1);
        Vmath::Vmul(npts, tmpx, 1, tmpx, 1, tmpx, 1);
        Vmath::Vvtvp(npts, tmpy, 1, tmpy, 1, tmpx, 1, tmpx, 1);

        int min_i = Vmath::Imin(npts, tmpx, 1);

        Array<OneD, NekDouble> eta(2, 0.);
        eta[0] = za[min_i % za.size()];
        eta[1] = zb[min_i / za.size()];
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
    // Debug to print objective function values
    if (false)
    {
        std::cout << std::scientific << std::setprecision(10) << "Looking for point: " << xs[0] << ", " << xs[1] << std::endl;
        std::cout << "Triangle: " << m_globalID << " contains point: " << ContainsPoint(xs) << std::endl;
        // triangle verts
        Array<OneD, NekDouble> xt(3), yt(3), zt(3);
        for (int i = 0; i < 3; ++i)
        {
            m_verts[i]->GetCoords(xt[i], yt[i], zt[i]);
        }

        ofstream file_xtri;
        ofstream file_ytri;
        file_xtri.open("xtri.txt");
        file_ytri.open("ytri.txt");
        file_xtri << xt[0] << std::endl << xt[1] << std::endl << xt[2] << std::endl << xt[0] << std::endl;
        file_ytri << yt[0] << std::endl << yt[1] << std::endl << yt[2] << std::endl << yt[0] << std::endl;
        file_xtri.close();
        file_ytri.close();


        // Print curve point locations
        ofstream file_xcurve;
        ofstream file_ycurve;
        file_xcurve.open("xcurve.txt");
        file_ycurve.open("ycurve.txt");
        for (auto & m_point : m_curve->m_points)
        {
            file_xcurve << m_point->x() << std::endl;
            file_ycurve << m_point->y() << std::endl;
        }
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

        int n = 100;
        for (int i = 0; i < n + 1; ++i)
        {
            xi[0] = -1.0 + 2.0 * i / n;
            for (int j = 0; j < n + 1; ++j)
            {
                xi[1] = -1.0 + 2.0 * j / n;

                Array<OneD, NekDouble> eta(2);
                m_xmap->LocCollapsedToLocCoord(xi, eta);
                xi = eta;

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
    else if (false)
    {
        std::cout << "Triangle: " << m_globalID << " contains point: " << ContainsPoint(xs) << std::endl;

        int nq = m_xmap->GetTotPoints();

        Array<OneD, NekDouble> x(nq), y(nq);
        m_xmap->BwdTrans(m_coeffs[0], x);
        m_xmap->BwdTrans(m_coeffs[1], y);

        Array<OneD, NekDouble> xi(2, 0.0);
        // Pick quad point closest to xs for starter xi
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
        std::cout << "Closest quad point: " << xi[0] << " " << xi[1] << std::endl;

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

        //ofstream file_xc_newton;
        //ofstream file_yc_newton;
        //file_xc_newton.open("xc_newton.txt");
        //file_yc_newton.open("yc_newton.txt");
        for (int i = 0; i < 100; ++i)
        {
            // Compute f(x_k) and its derivatives
            NekDouble xc = m_xmap->PhysEvaluate(xi, x);
            NekDouble yc = m_xmap->PhysEvaluate(xi, y);
            //file_xc_newton << xc << std::endl;
            //file_yc_newton << yc << std::endl;

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

            xi[0] = xi[0] - (hessInv[0][0] * jac[0] + hessInv[0][1] * jac[1]);
            xi[1] = xi[1] - (hessInv[1][0] * jac[0] + hessInv[1][1] * jac[1]);

            // Gradient descent
            // xi[0] = xi[0] - jac[0];

            if (xi[0] < -1 || xi[0] > 1 ||
                xi[1] < -1 || xi[1] > 1)
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

        //file_xc_newton.close();
        //file_yc_newton.close();

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
    else if(m_geomFactors->GetGtype() == eDeformed && false)
    {
        xiOut = Array<OneD, NekDouble>(2,0.0);

        GetLocCoords(xs, xiOut);
        ClampLocCoords(xiOut);

        Array<OneD, NekDouble> gloCoord(3);
        gloCoord[0] = GetCoord(0, xiOut);
        gloCoord[1] = GetCoord(1, xiOut);
        gloCoord[2] = GetCoord(2, xiOut);

        return sqrt((xs[0] - gloCoord[0])*(xs[0] - gloCoord[0]) +
                    (xs[1] - gloCoord[1])*(xs[1] - gloCoord[1]) +
                    (xs[2] - gloCoord[2])*(xs[2] - gloCoord[2]));
    }
    else if (m_geomFactors->GetGtype() == eDeformed)
    {
        Array<OneD, NekDouble> xi(2, 0.0), eta(2, 0.0);
        m_xmap->LocCollapsedToLocCoord(eta, xi);

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

        //ofstream file_xc_newton;
        //ofstream file_yc_newton;
        //file_xc_newton.open("xc_newton.txt");
        //file_yc_newton.open("yc_newton.txt");
        NekDouble fx_prev = std::numeric_limits<NekDouble>::max();
        for (int i = 0; i < 100; ++i)
        {
            // Compute f(x_k) and its derivatives
            NekDouble xc = m_xmap->PhysEvaluate(xi, x);
            NekDouble yc = m_xmap->PhysEvaluate(xi, y);

            //file_xc_newton << xc << std::endl;
            //file_yc_newton << yc << std::endl;

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

            //std::cout << "xi[0] = " << xi[0] << ", xi[1] = " << xi[1] << ", xc = " << xc << ", yc = " << yc << ", fx = " << fx << ", fx_diff = " << abs(fx - fx_prev) << std::endl;

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
            pk[1] = -(hessInv[0][1] * jac[0] + hessInv[1][1] * jac[1]);

            //std::cout << "pk[0] = " << pk[0] << ", pk[1] = " << pk[1] << std::endl;

            // Backtracking line search
            while (gamma > 1e-16)
            {
                Array<OneD, NekDouble> xi_pk(2);
                xi_pk[0] = xi[0] + pk[0] * gamma;
                xi_pk[1] = xi[1] + pk[1] * gamma;

                Array<OneD, NekDouble> eta_pk(2, 0.0);
                m_xmap->LocCoordToLocCollapsed(xi_pk, eta_pk);

                if (eta_pk[0] < (-1 - std::numeric_limits<NekDouble>::epsilon()) ||
                    eta_pk[0] > ( 1 + std::numeric_limits<NekDouble>::epsilon()) ||
                    eta_pk[1] < (-1 - std::numeric_limits<NekDouble>::epsilon()) ||
                    eta_pk[1] > ( 1 + std::numeric_limits<NekDouble>::epsilon()))
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

                //std::cout << "xi_pk[0] = " << xi_pk[0] << ", xi_pk[1] = " << xi_pk[1] <<" xc_pk = " << xc_pk << ", yc_pk = " << yc_pk << ", fx_pk = " << fx_pk << std::endl;

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
                break;
            }

            xi[0] += gamma * pk[0];
            xi[1] += gamma * pk[1];

            Array<OneD, NekDouble> eta(2, 0.0);
            m_xmap->LocCoordToLocCollapsed(xi, eta);
            if (eta[0] < (-1 - std::numeric_limits<NekDouble>::epsilon()) ||
                eta[0] > ( 1 + std::numeric_limits<NekDouble>::epsilon()) ||
                eta[1] < (-1 - std::numeric_limits<NekDouble>::epsilon()) ||
                eta[1] > ( 1 + std::numeric_limits<NekDouble>::epsilon()))
            {
                std::cout << "BORKEN" << std::endl;
                exit(0);
            }
        }

        //file_xc_newton.close();
        //file_yc_newton.close();


        xiOut = xi;
        return  sqrt(fx_prev);

    }
    else
    {
        ASSERTL0(false, "Geometry type unknown")
    }

    return -1.0;
}

}
}
