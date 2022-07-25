///////////////////////////////////////////////////////////////////////////////
//
// File FilterStructurePres.cpp
//
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
// The above copyright notice and this permission notice shall be included
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
// Description: Output Structure preserved coeffs.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/Filters/FilterStructurePres.h>
#include <boost/core/ignore_unused.hpp>

namespace Nektar
{
namespace SolverUtils
{

std::string FilterStructurePres::className =
    GetFilterFactory().RegisterCreatorFunction("StructurePres",
                                               FilterStructurePres::create);

FilterStructurePres::FilterStructurePres(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem> &pEquation, const ParamMap &pParams)
    : Filter(pSession, pEquation)
{
    // setup filter
    auto it = pParams.find("c_gd");
    if (it != pParams.end())
    {
        demo.Setchold(std::stod(it->second));
    }

    it = pParams.find("gam_gd");
    if (it != pParams.end())
    {
        demo.Setgamhold(std::stod(it->second));
    }
    it = pParams.find("tol_gd");
    if (it != pParams.end())
    {
        demo.SettolGD(std::stod(it->second));
    }

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = 1;
    }
    roots1dtimehold = 0.0, roots2dtimehold = 0.0, roots3dtimehold = 0.0,
    totalroottime = 0.0;
}

FilterStructurePres::~FilterStructurePres()
{
    cout << "\n time spent inside str pres filter: \t" << timeStrPres;
    cout << "\n time spent checking for target elements: \t" << timeFindIds;
    cout << "\n time spent for orthonormalization and reverse : \t" << timeOrth
         << " ";
    cout << "\n number of iterations for GD (all) = " << retGDall;
    cout << "\n time spent inside basis eval = " << basiscalls;

    //      cout<<"\n time taken inside GD  (all) = "<<retGDtime;
    cout << "\n number of iterations for filter = " << tot_optiter << " ";
    cout << "\n number of flagged ele in simulation = " << total_flagged << " ";
    cout << "\n time taken by root finding (all) = " << totalroottime << " ";
    cout << "\n time taken by 1D rootfinding = " << roots1dtimehold << " ";
    cout << "\n time spent for basis eval = " << roots3dtimehold << " ";
    cout << "\n post-setup these many points are interpolated = "
         << tot_interp_points << " ";

    if (dimension > 2)
    {
        cout << "\n time spent inside 3d GD = " << roots3dtimehold << " ";
    }
    if (dimension > 1)
    {
        cout << "\n time spent inside 2d GD = " << roots2dtimehold << "\n";
    }
}

void FilterStructurePres::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(time);

    nfields = pFields.size(); // assumption: only 1 variable

    // loop this for every variable and initialize only the ptrs that are still
    // null

    nelmt = pFields[0]->GetExpSize();

    MultiRegions::ExpansionType exptype = pFields[0]->GetExpType();
    if (E3seg == nullptr)
    {
        Array<OneD, int> nmodesperel = pFields[0]->EvalBasisNumModesMaxPerExp();

        int orderseg3 = 3 * Vmath::Vmax(nmodesperel.size(), nmodesperel, 1);
        call_setup_seg(pFields[0]->GetExp(0), orderseg3);
    }

    // colleague matrix
    C = demo.formConf(E3seg->GetBasis(0)->GetNumModes());
    // uncomment this for using companion matrix approach
    // C = demo.formCompanion(E3seg->GetBasis(0)->GetNumModes());
    Timer tmr;
    tmr.Start();
    if (exptype > 1)
    {
        for (int k = 0; k < 1; k++)
        {
            listeletype = Array<OneD, StdExpansion *>(nelmt);
            // find out what is the max num of quad points in x-dir
            for (int i = 0; i < nelmt; i++)
            {

                LocalRegions::ExpansionSharedPtr exp = pFields[k]->GetExp(i);

                switch (exp->DetShapeType())
                {
                    case LibUtilities::eHexahedron:

                        if (Equad == nullptr)
                        {
                            call_setup_quad(exp);
                        }

                        if (Ehex == nullptr)
                        {
                            call_setup_hex(exp);
                        }
                        listeletype[i] = Ehex;
                        break;
                    case LibUtilities::eTetrahedron:

                        if (Etri == nullptr)
                        {
                            call_setup_tri(exp);
                        }
                        if (Etet == nullptr)
                        {
                            call_setup_tet(exp);
                        }

                        listeletype[i] = Etet;
                        break;
                    case LibUtilities::ePyramid:

                        if (Etri == nullptr)
                        {
                            call_setup_tri(exp);
                        }
                        if (Equad == nullptr)
                        {
                            call_setup_quad(exp);
                        }
                        if (Epyr == nullptr)
                        {
                            call_setup_pyr(exp);
                        }
                        listeletype[i] = Epyr;
                        break;

                    case LibUtilities::ePrism:
                        if (Etri == nullptr)
                        {
                            call_setup_tri(exp);
                        }
                        if (Equad == nullptr)
                        {
                            call_setup_quad(exp);
                        }
                        if (Epri == nullptr)
                        {
                            call_setup_pri(exp);
                        }

                        listeletype[i] = Epri;
                        break;
                    case LibUtilities::eQuadrilateral:
                        dimension = 2;
                        if (Equad == nullptr)
                        {
                            call_setup_quad(exp);
                        }
                        listeletype[i] = Equad;
                        break;

                    case LibUtilities::eTriangle:
                        if (Etri == nullptr)
                        {
                            call_setup_tri(exp);
                        }

                        listeletype[i] = Etri;
                        break;

                    default:
                        cout << "\n unsupported element for this filter!\n";
                        exit(0);
                        break;
                }
            }
        }
    }
    v_Update(pFields, time);
}

void FilterStructurePres::call_setup_tri(LocalRegions::ExpansionSharedPtr exp)
{
    int nmodes0 = exp->GetBasis(0)->GetNumModes();

    int nmodes1 = exp->GetBasis(1)->GetNumModes();
    int npts0   = exp->GetBasis(0)->GetNumPoints();
    int npts1   = exp->GetBasis(1)->GetNumPoints();
    PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
    PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
    BasisKey b0(LibUtilities::eOrtho_A, nmodes0, p0);
    BasisKey b1(LibUtilities::eOrtho_B, nmodes1, p1);
    Etri = new StdTriExp(b0, b1);

    demo.storage2dt  = Etri->GetPhysEvaluateStorage();
    demo.coordtri    = demo.GetCoords(Etri);
    demo.coordmidtri = demo.GetQuadratureMidCoords(demo.coordtri);
    demo.coordlatticetri =
        demo.GetLatticeCoords(demo.coordtri, demo.coordmidtri);

    demo.midptevaltri = Etri->PhysEvaluateBasis(
        demo.coordmidtri, demo.storage2dt, NullNekDouble1DArray,
        NullNekDouble1DArray, NullNekDouble1DArray);

    Array<OneD, Array<OneD, NekDouble>> edgeptsin(2), edgexy(2);
    for (int k = 0; k < 2; k++)
    {
        edgexy[k] = Array<OneD, NekDouble>(E3seg->GetBasis(0)->GetNumPoints());
        edgeptsin[k] = Array<OneD, NekDouble>(edgexy[k]);
    }
    Array<OneD, NekDouble> edgexytemp = E3seg->GetBasis(0)->GetZ();
    int totszedges = edgexytemp.size() * (Etri->GetNcoeffs());

    // left (x = -1)
    Vxm1t   = Array<OneD, NekDouble>(totszedges);
    Vdyxm1t = Array<OneD, NekDouble>(totszedges);
    Vdxxm1t = Array<OneD, NekDouble>(totszedges);

    // left x = -1
    edgexy[1] = edgexytemp;
    edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
    Vxm1t = Etri->PhysEvaluateBasis(edgexy, demo.storage2dt, Vdxxm1t, Vdyxm1t,
                                    NullNekDouble1DArray);

    // bot (y = -1)
    Vym1t   = Array<OneD, NekDouble>(totszedges);
    Vdxym1t = Array<OneD, NekDouble>(totszedges);
    Vdyym1t = Array<OneD, NekDouble>(totszedges);

    // bot y = -1
    edgexy[0] = edgexytemp;
    edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);

    Vym1t = Etri->PhysEvaluateBasis(edgexy, demo.storage2dt, Vdxym1t, Vdyym1t,
                                    NullNekDouble1DArray);

    // hypt tri (y = -x)
    Vxyhypt   = Array<OneD, NekDouble>(totszedges);
    Vdxxyhypt = Array<OneD, NekDouble>(totszedges);
    Vdyxyhypt = Array<OneD, NekDouble>(totszedges);

    edgexy[0] = edgexytemp;
    Vmath::Smul(edgexy[0].size(), -1.0, edgexy[0], 1, edgexy[1], 1);
    Vxyhypt = Etri->PhysEvaluateBasis(edgexy, demo.storage2dt, Vdxxyhypt,
                                      Vdyxyhypt, NullNekDouble1DArray);
}

void FilterStructurePres::call_setup_seg(LocalRegions::ExpansionSharedPtr exp,
                                         int orderseg3)
{

    LibUtilities::PointsKey pkeycheb(exp->GetBasis(0)->GetNumPoints(),
                                     LibUtilities::eGaussLobattoChebyshev);
    LibUtilities::BasisKey bkeyorth(LibUtilities::eOrtho_A, orderseg3,
                                    pkeycheb);
    E3seg = new StdSegExp(bkeyorth);
}

void FilterStructurePres::call_setup_quad(LocalRegions::ExpansionSharedPtr exp)
{
    // get # of quad pts in curr exp
    // get order of current exp
    // get ptypes and btypes

    int nmodes0 = exp->GetBasis(0)->GetNumModes();
    int nmodes1 = exp->GetBasis(1)->GetNumModes();
    int npts0   = exp->GetBasis(0)->GetNumPoints();
    int npts1   = exp->GetBasis(1)->GetNumPoints();

    PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
    PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
    BasisKey b0(LibUtilities::eOrtho_A, nmodes0, p0);
    BasisKey b1(LibUtilities::eOrtho_A, nmodes1, p1);
    Equad             = new StdQuadExp(b0, b1);
    demo.storage2dq   = Equad->GetPhysEvaluateStorage();
    demo.coordquad    = demo.GetCoords(Equad);
    demo.coordmidquad = demo.GetQuadratureMidCoords(demo.coordquad);
    demo.coordlatticequad =
        demo.GetLatticeCoords(demo.coordquad, demo.coordmidquad);

    demo.midptevalquad = Equad->PhysEvaluateBasis(
        demo.coordmidquad, demo.storage2dq, NullNekDouble1DArray,
        NullNekDouble1DArray, NullNekDouble1DArray);

    Array<OneD, Array<OneD, NekDouble>> edgeptsin(2), edgexy(2);
    for (int k = 0; k < 2; k++)
    {
        edgexy[k] = Array<OneD, NekDouble>(E3seg->GetBasis(0)->GetNumPoints());
        edgeptsin[k] = Array<OneD, NekDouble>(edgexy[k]);
    }
    Array<OneD, NekDouble> edgexytemp = E3seg->GetBasis(0)->GetZ();
    int totszedges = edgexytemp.size() * (Equad->GetNcoeffs());

    // left (x = -1)
    Vxm1q   = Array<OneD, NekDouble>(totszedges);
    Vdyxm1q = Array<OneD, NekDouble>(totszedges);
    Vdxxm1q = Array<OneD, NekDouble>(totszedges);

    // left x = -1
    edgexy[1] = edgexytemp;
    edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);
    Vxm1q = Equad->PhysEvaluateBasis(edgexy, demo.storage2dq, Vdxxm1q, Vdyxm1q,
                                     NullNekDouble1DArray);

    // bot (y = -1)
    Vym1q   = Array<OneD, NekDouble>(totszedges);
    Vdxym1q = Array<OneD, NekDouble>(totszedges);
    Vdyym1q = Array<OneD, NekDouble>(totszedges);

    // bot y = -1
    edgexy[0] = edgexytemp;
    edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), -1.0);

    Vym1q = Equad->PhysEvaluateBasis(edgexy, demo.storage2dq, Vdxym1q, Vdyym1q,
                                     NullNekDouble1DArray);

    // right quad (x = 1)
    Vx1q   = Array<OneD, NekDouble>(totszedges);
    Vdxx1q = Array<OneD, NekDouble>(totszedges);
    Vdyx1q = Array<OneD, NekDouble>(totszedges);

    // top quad (y = 1)
    Vdxy1q = Array<OneD, NekDouble>(totszedges);
    Vdyy1q = Array<OneD, NekDouble>(totszedges);
    Vy1q   = Array<OneD, NekDouble>(totszedges);

    // right x = 1
    edgexy[1] = edgexytemp;
    edgexy[0] = Array<OneD, NekDouble>(edgexytemp.size(), 1.0);
    Vx1q = Equad->PhysEvaluateBasis(edgexy, demo.storage2dq, Vdxx1q, Vdyx1q,
                                    NullNekDouble1DArray);

    // top y = 1
    edgexy[0] = edgexytemp;
    edgexy[1] = Array<OneD, NekDouble>(edgexytemp.size(), 1.0);
    Vy1q = Equad->PhysEvaluateBasis(edgexy, demo.storage2dq, Vdxy1q, Vdyy1q,
                                    NullNekDouble1DArray);
}

void FilterStructurePres::call_setup_tet(LocalRegions::ExpansionSharedPtr exp)
{
    dimension   = 3;
    int nmodes0 = exp->GetBasis(0)->GetNumModes();
    int nmodes1 = exp->GetBasis(1)->GetNumModes();
    int nmodes2 = exp->GetBasis(2)->GetNumModes();
    int npts0   = exp->GetBasis(0)->GetNumPoints();
    int npts1   = exp->GetBasis(1)->GetNumPoints();
    int npts2   = exp->GetBasis(2)->GetNumPoints();

    PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
    PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
    PointsKey p2(npts2, exp->GetBasis(2)->GetPointsType());
    BasisKey b0(LibUtilities::eOrtho_A, nmodes0, p0);
    BasisKey b1(LibUtilities::eOrtho_B, nmodes1, p1);
    BasisKey b2(LibUtilities::eOrtho_C, nmodes2, p2);

    Etet = new StdTetExp(b0, b1, b2);

    demo.storage3dtet = Etet->GetPhysEvaluateStorage();
    demo.coordtet     = demo.GetCoords(Etet);
    demo.coordmidtet  = demo.GetQuadratureMidCoords(demo.coordtet);
    demo.coordlatticetri =
        demo.GetLatticeCoords(demo.coordtet, demo.coordmidtet);

    demo.midptevaltet = Etet->PhysEvaluateBasis(
        demo.coordmidtet, demo.storage3dtet, NullNekDouble1DArray,
        NullNekDouble1DArray, NullNekDouble1DArray);

    Array<OneD, NekDouble> edgexyztemp = E3seg->GetBasis(0)->GetZ();
    int totszedges1d = edgexyztemp.size() * (exp->GetNcoeffs());
    Array<OneD, Array<OneD, NekDouble>> edgeptsin(dimension);
    for (int p = 0; p < dimension; p++)
    {
        edgeptsin[p] = Array<OneD, NekDouble>(edgexyztemp);
    }

    // edge front left (AD) (x = -1) (y = -1)
    Vxm1ym1ztet   = Array<OneD, NekDouble>(totszedges1d);
    Vdxxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
    Vdyxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
    Vdzxm1ym1ztet = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[0]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    edgeptsin[1]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

    Vxm1ym1ztet =
        Etet->PhysEvaluateBasis(edgeptsin, demo.storage3dtet, Vdxxm1ym1ztet,
                                Vdyxm1ym1ztet, Vdzxm1ym1ztet);

    // edge front hypt (DB) (y = -1) (z = -x)
    Vym1xmztet   = Array<OneD, NekDouble>(totszedges1d);
    Vdxym1xmztet = Array<OneD, NekDouble>(totszedges1d);
    Vdyym1xmztet = Array<OneD, NekDouble>(totszedges1d);
    Vdzym1xmztet = Array<OneD, NekDouble>(totszedges1d);

    Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[2][0], 1,
                &edgeptsin[0][0], 1);

    Vym1xmztet = Etet->PhysEvaluateBasis(
        edgeptsin, demo.storage3dtet, Vdxym1xmztet, Vdyym1xmztet, Vdzym1xmztet);

    // edge front bot (AB) (y = -1) (z = -1)
    Vym1xzm1tet   = Array<OneD, NekDouble>(totszedges1d);
    Vdxym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
    Vdyym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
    Vdzym1xzm1tet = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[0]  = edgexyztemp;
    edgeptsin[2]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

    Vym1xzm1tet =
        Etet->PhysEvaluateBasis(edgeptsin, demo.storage3dtet, Vdxym1xzm1tet,
                                Vdyym1xzm1tet, Vdzym1xzm1tet);

    // edge left hypt (DC) ( x = -1) (z = -y)
    Vxm1ymztet   = Array<OneD, NekDouble>(totszedges1d);
    Vdxxm1ymztet = Array<OneD, NekDouble>(totszedges1d);
    Vdyxm1ymztet = Array<OneD, NekDouble>(totszedges1d);
    Vdzxm1ymztet = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    edgeptsin[1] = edgexyztemp;
    Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[1][0], 1,
                &edgeptsin[2][0], 1);

    Vxm1ymztet = Etet->PhysEvaluateBasis(
        edgeptsin, demo.storage3dtet, Vdxxm1ymztet, Vdyxm1ymztet, Vdzxm1ymztet);

    // edge bot diag (BC) (z = -1) (y = -x)
    Vxmyzm1tet   = Array<OneD, NekDouble>(totszedges1d);
    Vdxxmyzm1tet = Array<OneD, NekDouble>(totszedges1d);
    Vdyxmyzm1tet = Array<OneD, NekDouble>(totszedges1d);
    Vdzxmyzm1tet = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[0] = edgexyztemp;
    Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[0][0], 1,
                &edgeptsin[1][0], 1);
    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[2].size(), -1.0);
    Vxmyzm1tet   = Etet->PhysEvaluateBasis(
          edgeptsin, demo.storage3dtet, Vdxxmyzm1tet, Vdyxmyzm1tet, Vdzxmyzm1tet);

    // edge CA bot left (x = -1) (z = -1)
    Vxm1yzm1tet   = Array<OneD, NekDouble>(totszedges1d);
    Vdxxm1yzm1tet = Array<OneD, NekDouble>(totszedges1d);
    Vdyxm1yzm1tet = Array<OneD, NekDouble>(totszedges1d);
    Vdzxm1yzm1tet = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[1]  = edgexyztemp;
    edgeptsin[2]  = Array<OneD, NekDouble>(edgeptsin[2].size(), -1.0);

    Vxm1yzm1tet =
        Etet->PhysEvaluateBasis(edgeptsin, demo.storage3dtet, Vdxxm1yzm1tet,
                                Vdyxm1yzm1tet, Vdzxm1yzm1tet);
    int totsurf2d = (demo.coordtri[0].size()) * Etet->GetNcoeffs();
    Array<OneD, Array<OneD, NekDouble>> surfptsin(dimension),
        surfptsintemp(dimension);

    for (int k = 0; k < dimension - 1; k++)
    {
        surfptsin[k]     = Array<OneD, NekDouble>(demo.coordtri[k]);
        surfptsintemp[k] = Array<OneD, NekDouble>(demo.coordtri[k]);
    }
    surfptsin[dimension - 1] = Array<OneD, NekDouble>(demo.coordtri[0].size());
    surfptsintemp[dimension - 1] =
        Array<OneD, NekDouble>(demo.coordtri[0].size());

    int totpt = surfptsin[0].size();

    // surface bot z = -1, (ABC)
    Vxyzm1tet        = Array<OneD, NekDouble>(totsurf2d);
    surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
    surfptsintemp[0] = surfptsin[0];
    surfptsintemp[1] = surfptsin[1];
    Vxyzm1tet        = Etet->PhysEvaluateBasis(
               surfptsintemp, demo.storage3dtet, NullNekDouble1DArray,
               NullNekDouble1DArray, NullNekDouble1DArray);

    // surface left x = -1  (DAC)
    Vxm1ypz0tet      = Array<OneD, NekDouble>(totsurf2d);
    surfptsintemp[0] = Array<OneD, NekDouble>(totpt, -1.0);
    surfptsintemp[1] = surfptsin[0];
    surfptsintemp[2] = surfptsin[1];

    Vxm1ypz0tet = Etet->PhysEvaluateBasis(
        surfptsintemp, demo.storage3dtet, NullNekDouble1DArray,
        NullNekDouble1DArray, NullNekDouble1DArray);

    // surf front y = -1,  (DAB)
    Vxpz0ym1tet      = Array<OneD, NekDouble>(totsurf2d);
    surfptsintemp[1] = Array<OneD, NekDouble>(totpt, -1.0);
    surfptsintemp[0] = surfptsin[0];
    surfptsintemp[2] = surfptsin[1];
    Vxpz0ym1tet      = Etet->PhysEvaluateBasis(
             surfptsintemp, demo.storage3dtet, NullNekDouble1DArray,
             NullNekDouble1DArray, NullNekDouble1DArray);

    // surf DCB (x + y + z = -1),
    Vxpypzm1tet      = Array<OneD, NekDouble>(totsurf2d);
    surfptsintemp[1] = surfptsin[1];
    ;

    Vmath::Vadd(totpt, &surfptsin[0][0], 1, &surfptsin[1][0], 1,
                &surfptsintemp[2][0], 1);
    Vmath::Smul(totpt, -1.0, &surfptsintemp[2][0], 1, &surfptsintemp[2][0], 1);
    Vmath::Sadd(totpt, -1.0, &surfptsintemp[2][0], 1, &surfptsintemp[2][0], 1);
    Vxpypzm1tet = Etet->PhysEvaluateBasis(
        surfptsintemp, demo.storage3dtet, NullNekDouble1DArray,
        NullNekDouble1DArray, NullNekDouble1DArray);
}

void FilterStructurePres::call_setup_hex(LocalRegions::ExpansionSharedPtr exp)
{
    dimension = 3;

    int nmodes0 = exp->GetBasis(0)->GetNumModes();
    int nmodes1 = exp->GetBasis(1)->GetNumModes();
    int nmodes2 = exp->GetBasis(2)->GetNumModes();
    int npts0   = exp->GetBasis(0)->GetNumPoints();
    int npts1   = exp->GetBasis(1)->GetNumPoints();
    int npts2   = exp->GetBasis(2)->GetNumPoints();
    PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
    PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
    PointsKey p2(npts2, exp->GetBasis(2)->GetPointsType());
    BasisKey b0(LibUtilities::eOrtho_A, nmodes0, p0);
    BasisKey b1(LibUtilities::eOrtho_A, nmodes1, p1);
    BasisKey b2(LibUtilities::eOrtho_A, nmodes2, p2);

    Ehex              = new StdHexExp(b0, b1, b2);
    demo.storage3dhex = Ehex->GetPhysEvaluateStorage();
    tot_interp_points += Ehex->GetTotPoints();
    demo.coordhex    = demo.GetCoords(Ehex);
    demo.coordmidhex = demo.GetQuadratureMidCoords(demo.coordhex);
    demo.coordlatticehex =
        demo.GetLatticeCoords(demo.coordhex, demo.coordmidhex);

    demo.midptevalhex = Ehex->PhysEvaluateBasis(
        demo.coordmidhex, demo.storage3dhex, NullNekDouble1DArray,
        NullNekDouble1DArray, NullNekDouble1DArray);

    Array<OneD, NekDouble> edgexyztemp = E3seg->GetBasis(0)->GetZ();

    int totszedges1d = edgexyztemp.size() * (Ehex->GetNcoeffs());
    Array<OneD, Array<OneD, NekDouble>> edgeptsin(dimension);
    for (int p = 0; p < dimension; p++)
    {
        edgeptsin[p] = Array<OneD, NekDouble>(edgexyztemp);
    }

    // edge front left (x = -1) (y = -1)
    Vxm1ym1z   = Array<OneD, NekDouble>(totszedges1d);
    Vdxxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
    Vdyxm1ym1z = Array<OneD, NekDouble>(totszedges1d);
    Vdzxm1ym1z = Array<OneD, NekDouble>(totszedges1d);

    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

    edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

    Vxm1ym1z = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxxm1ym1z,
                                       Vdyxm1ym1z, Vdzxm1ym1z);

    // edge front right (x = 1) (y = -1)
    Vx1ym1z      = Array<OneD, NekDouble>(totszedges1d);
    Vdxx1ym1z    = Array<OneD, NekDouble>(totszedges1d);
    Vdyx1ym1z    = Array<OneD, NekDouble>(totszedges1d);
    Vdzx1ym1z    = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
    Vx1ym1z = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxx1ym1z,
                                      Vdyx1ym1z, Vdzx1ym1z);

    // edge front top (y = -1) (z = 1)
    Vym1xz1      = Array<OneD, NekDouble>(totszedges1d);
    Vdxym1xz1    = Array<OneD, NekDouble>(totszedges1d);
    Vdyym1xz1    = Array<OneD, NekDouble>(totszedges1d);
    Vdzym1xz1    = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[0] = edgexyztemp;
    edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
    Vym1xz1 = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxym1xz1,
                                      Vdyym1xz1, Vdzym1xz1);

    // edge front bot (y = -1) (z = -1)
    Vym1xzm1     = Array<OneD, NekDouble>(totszedges1d);
    Vdxym1xzm1   = Array<OneD, NekDouble>(totszedges1d);
    Vdyym1xzm1   = Array<OneD, NekDouble>(totszedges1d);
    Vdzym1xzm1   = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    Vym1xzm1 = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxym1xzm1,
                                       Vdyym1xzm1, Vdzym1xzm1);

    // edge back left (y = 1), (x = -1)
    Vxm1y1z      = Array<OneD, NekDouble>(totszedges1d);
    Vdxxm1y1z    = Array<OneD, NekDouble>(totszedges1d);
    Vdyxm1y1z    = Array<OneD, NekDouble>(totszedges1d);
    Vdzxm1y1z    = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    edgeptsin[2] = edgexyztemp;

    Vxm1y1z = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxxm1y1z,
                                      Vdyxm1y1z, Vdzxm1y1z);

    // edge back right (x = 1), (y = 1))
    Vx1y1z       = Array<OneD, NekDouble>(totszedges1d);
    Vdxx1y1z     = Array<OneD, NekDouble>(totszedges1d);
    Vdyx1y1z     = Array<OneD, NekDouble>(totszedges1d);
    Vdzx1y1z     = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);

    Vx1y1z = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxx1y1z,
                                     Vdyx1y1z, Vdzx1y1z);

    // edge back top ( y = 1) (z = 1)
    Vy1xz1       = Array<OneD, NekDouble>(totszedges1d);
    Vdxy1xz1     = Array<OneD, NekDouble>(totszedges1d);
    Vdyy1xz1     = Array<OneD, NekDouble>(totszedges1d);
    Vdzy1xz1     = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[0] = edgexyztemp;
    edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
    Vy1xz1 = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxy1xz1,
                                     Vdyy1xz1, Vdzy1xz1);

    // edge back bot (y = 1) (z = -1))
    Vy1xzm1      = Array<OneD, NekDouble>(totszedges1d);
    Vdxy1xzm1    = Array<OneD, NekDouble>(totszedges1d);
    Vdyy1xzm1    = Array<OneD, NekDouble>(totszedges1d);
    Vdzy1xzm1    = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

    Vy1xzm1 = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxy1xzm1,
                                      Vdyy1xzm1, Vdzy1xzm1);

    // edge left bot (z = -1), (x = -1)
    Vxm1yzm1     = Array<OneD, NekDouble>(totszedges1d);
    Vdxxm1yzm1   = Array<OneD, NekDouble>(totszedges1d);
    Vdyxm1yzm1   = Array<OneD, NekDouble>(totszedges1d);
    Vdzxm1yzm1   = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[1] = edgexyztemp;
    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);

    Vxm1yzm1 = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxxm1yzm1,
                                       Vdyxm1yzm1, Vdzxm1yzm1);

    // edge left top (x = -1), (z = 1))
    Vxm1yz1   = Array<OneD, NekDouble>(totszedges1d);
    Vdxxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
    Vdyxm1yz1 = Array<OneD, NekDouble>(totszedges1d);
    Vdzxm1yz1 = Array<OneD, NekDouble>(totszedges1d);

    edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
    Vxm1yz1 = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxxm1yz1,
                                      Vdyxm1yz1, Vdzxm1yz1);

    // edge right bot ( z = -1) (x = 1)
    Vx1yzm1   = Array<OneD, NekDouble>(totszedges1d);
    Vdxx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
    Vdyx1yzm1 = Array<OneD, NekDouble>(totszedges1d);
    Vdzx1yzm1 = Array<OneD, NekDouble>(totszedges1d);

    edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
    Vx1yzm1 = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxx1yzm1,
                                      Vdyx1yzm1, Vdzx1yzm1);

    // edge right top (z  1) (x  1))
    Vx1yz1   = Array<OneD, NekDouble>(totszedges1d);
    Vdxx1yz1 = Array<OneD, NekDouble>(totszedges1d);
    Vdyx1yz1 = Array<OneD, NekDouble>(totszedges1d);
    Vdzx1yz1 = Array<OneD, NekDouble>(totszedges1d);

    edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);

    Vx1yz1 = Ehex->PhysEvaluateBasis(edgeptsin, demo.storage3dhex, Vdxx1yz1,
                                     Vdyx1yz1, Vdzx1yz1);
    int totszsurf2d = (demo.coordquad[0].size()) * Ehex->GetNcoeffs();

    Array<OneD, Array<OneD, NekDouble>> surfptsin(dimension),
        surfptsintemp(dimension);
    for (int k = 0; k < dimension - 1; k++)
    {
        surfptsin[k]     = Array<OneD, NekDouble>(demo.coordquad[k]);
        surfptsintemp[k] = Array<OneD, NekDouble>(demo.coordquad[k]);
    }

    surfptsin[dimension - 1] = Array<OneD, NekDouble>(demo.coordquad[0].size());
    surfptsintemp[dimension - 1] =
        Array<OneD, NekDouble>(demo.coordquad[0].size());

    // surface bot z = -1
    Vxyzm1           = Array<OneD, NekDouble>(totszsurf2d);
    surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
    surfptsintemp[0] = surfptsin[0];
    surfptsintemp[1] = surfptsin[1];
    Vxyzm1           = Ehex->PhysEvaluateBasis(surfptsintemp, demo.storage3dhex,
                                               NullNekDouble1DArray, NullNekDouble1DArray,
                                               NullNekDouble1DArray);

    // surface right x = 1
    Vx1yz            = Array<OneD, NekDouble>(totszsurf2d);
    surfptsintemp[0] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
    surfptsintemp[1] = surfptsin[0];
    surfptsintemp[2] = surfptsin[1];
    Vx1yz            = Ehex->PhysEvaluateBasis(surfptsintemp, demo.storage3dhex,
                                               NullNekDouble1DArray, NullNekDouble1DArray,
                                               NullNekDouble1DArray);

    // surface top z = 1
    Vxyz1            = Array<OneD, NekDouble>(totszsurf2d);
    surfptsintemp[0] = surfptsin[0];
    surfptsintemp[1] = surfptsin[1];
    surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
    Vxyz1            = Ehex->PhysEvaluateBasis(surfptsintemp, demo.storage3dhex,
                                               NullNekDouble1DArray, NullNekDouble1DArray,
                                               NullNekDouble1DArray);

    // surface left x = -1
    Vxm1yz           = Array<OneD, NekDouble>(totszsurf2d);
    surfptsintemp[0] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
    surfptsintemp[1] = surfptsin[0];
    surfptsintemp[2] = surfptsin[1];
    Vxm1yz           = Ehex->PhysEvaluateBasis(surfptsintemp, demo.storage3dhex,
                                               NullNekDouble1DArray, NullNekDouble1DArray,
                                               NullNekDouble1DArray);

    // surface front y = -1
    Vxym1z           = Array<OneD, NekDouble>(totszsurf2d);
    surfptsintemp[1] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
    surfptsintemp[0] = surfptsin[0];
    surfptsintemp[2] = surfptsin[1];
    Vxym1z           = Ehex->PhysEvaluateBasis(surfptsintemp, demo.storage3dhex,
                                               NullNekDouble1DArray, NullNekDouble1DArray,
                                               NullNekDouble1DArray);

    // surface back y = 1
    Vxy1z            = Array<OneD, NekDouble>(totszsurf2d);
    surfptsintemp[1] = Array<OneD, NekDouble>(surfptsin[0].size(), 1.0);
    surfptsintemp[0] = surfptsin[0];
    surfptsintemp[2] = surfptsin[1];
    Vxy1z            = Ehex->PhysEvaluateBasis(surfptsintemp, demo.storage3dhex,
                                               NullNekDouble1DArray, NullNekDouble1DArray,
                                               NullNekDouble1DArray);
}

void FilterStructurePres::call_setup_pri(LocalRegions::ExpansionSharedPtr exp)
{
    dimension   = 3;
    int nmodes0 = exp->GetBasis(0)->GetNumModes();
    int nmodes1 = exp->GetBasis(1)->GetNumModes();
    int nmodes2 = exp->GetBasis(2)->GetNumModes();
    int npts0   = exp->GetBasis(0)->GetNumPoints();
    int npts1   = exp->GetBasis(1)->GetNumPoints();
    int npts2   = exp->GetBasis(2)->GetNumPoints();
    PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
    PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
    PointsKey p2(npts2, exp->GetBasis(2)->GetPointsType());
    BasisKey b0(LibUtilities::eOrtho_A, nmodes0, p0);
    BasisKey b1(LibUtilities::eOrtho_A, nmodes1, p1);
    BasisKey b2(LibUtilities::eOrtho_B, nmodes2, p2);

    Epri              = new StdPrismExp(b0, b1, b2);
    demo.storage3dpri = Epri->GetPhysEvaluateStorage();
    demo.coordpri     = demo.GetCoords(Epri);
    demo.coordmidpri  = demo.GetQuadratureMidCoords(demo.coordpri);
    demo.coordlatticepri =
        demo.GetLatticeCoords(demo.coordpri, demo.coordmidpri);
    demo.midptevalpri = Epri->PhysEvaluateBasis(
        demo.coordmidpri, demo.storage3dpri, NullNekDouble1DArray,
        NullNekDouble1DArray, NullNekDouble1DArray);

    // we need qzin + qzinmid as input points
    Array<OneD, NekDouble> edgexyztemp = E3seg->GetBasis(0)->GetZ();
    int totszedges1d = edgexyztemp.size() * (Epri->GetNcoeffs());
    Array<OneD, Array<OneD, NekDouble>> edgeptsin(dimension);
    for (int p = 0; p < dimension; p++)
    {
        edgeptsin[p] = Array<OneD, NekDouble>(edgexyztemp);
    }
    // edge front left EA (x = -1) (y = -1)
    Vxm1ym1zpri   = Array<OneD, NekDouble>(totszedges1d);
    Vdxxm1ym1zpri = Array<OneD, NekDouble>(totszedges1d);
    Vdyxm1ym1zpri = Array<OneD, NekDouble>(totszedges1d);
    Vdzxm1ym1zpri = Array<OneD, NekDouble>(totszedges1d);

    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    Vxm1ym1zpri =
        Epri->PhysEvaluateBasis(edgeptsin, demo.storage3dpri, Vdxxm1ym1zpri,
                                Vdyxm1ym1zpri, Vdzxm1ym1zpri);

    // edge front bot (y = -1) (z = -1) AB
    //       Array<OneD, NekDouble> Vym1xzm1pri   ;
    Vym1xzm1pri   = Array<OneD, NekDouble>(totszedges1d);
    Vdxym1xzm1pri = Array<OneD, NekDouble>(totszedges1d);
    Vdyym1xzm1pri = Array<OneD, NekDouble>(totszedges1d);
    Vdzym1xzm1pri = Array<OneD, NekDouble>(totszedges1d);

    edgeptsin[0] = Array<OneD, NekDouble>(edgexyztemp);
    edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    Vym1xzm1pri =
        Epri->PhysEvaluateBasis(edgeptsin, demo.storage3dpri, Vdxym1xzm1pri,
                                Vdyym1xzm1pri, Vdzym1xzm1pri);

    // edge front hypt (y = -1) (z = -x) EB
    Vym1xmzpri   = Array<OneD, NekDouble>(totszedges1d);
    Vdxym1xmzpri = Array<OneD, NekDouble>(totszedges1d);
    Vdyym1xmzpri = Array<OneD, NekDouble>(totszedges1d);
    Vdzym1xmzpri = Array<OneD, NekDouble>(totszedges1d);
    Vmath::Smul(edgeptsin[0].size(), -1.0, edgeptsin[0], 1, edgeptsin[2], 1);
    Vym1xmzpri = Epri->PhysEvaluateBasis(
        edgeptsin, demo.storage3dpri, Vdxym1xmzpri, Vdyym1xmzpri, Vdzym1xmzpri);

    // edge back bot (y = 1) (z = -1)) DC
    Vy1xzm1pri   = Array<OneD, NekDouble>(totszedges1d);
    Vdxy1xzm1pri = Array<OneD, NekDouble>(totszedges1d);
    Vdyy1xzm1pri = Array<OneD, NekDouble>(totszedges1d);
    Vdzy1xzm1pri = Array<OneD, NekDouble>(totszedges1d);

    edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
    edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    Vy1xzm1pri   = Epri->PhysEvaluateBasis(
          edgeptsin, demo.storage3dpri, Vdxy1xzm1pri, Vdyy1xzm1pri, Vdzy1xzm1pri);

    // edge back left (y = 1) (x = -1))
    Vxm1y1zpri   = Array<OneD, NekDouble>(totszedges1d);
    Vdxxm1y1zpri = Array<OneD, NekDouble>(totszedges1d);
    Vdyxm1y1zpri = Array<OneD, NekDouble>(totszedges1d);
    Vdzxm1y1zpri = Array<OneD, NekDouble>(totszedges1d);

    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    edgeptsin[2] = Array<OneD, NekDouble>(edgexyztemp);
    Vxm1y1zpri   = Epri->PhysEvaluateBasis(
          edgeptsin, demo.storage3dpri, Vdxxm1y1zpri, Vdyxm1y1zpri, Vdzxm1y1zpri);

    // edge back hypt (y = 1) (z = -x)
    Vy1xmzpri    = Array<OneD, NekDouble>(totszedges1d);
    Vdxy1xmzpri  = Array<OneD, NekDouble>(totszedges1d);
    Vdyy1xmzpri  = Array<OneD, NekDouble>(totszedges1d);
    Vdzy1xmzpri  = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[0] = Array<OneD, NekDouble>(edgexyztemp);
    Vmath::Smul(edgeptsin[0].size(), -1.0, edgeptsin[0], 1, edgeptsin[2], 1);
    Vy1xmzpri = Epri->PhysEvaluateBasis(edgeptsin, demo.storage3dpri,
                                        Vdxy1xmzpri, Vdyy1xmzpri, Vdzy1xmzpri);
    // edge left top (x = -1) (z = 1)
    Vxm1yz1pri = Array<OneD, NekDouble>(totszedges1d);
    ;
    Vdxxm1yz1pri = Array<OneD, NekDouble>(totszedges1d);
    ;
    Vdyxm1yz1pri = Array<OneD, NekDouble>(totszedges1d);
    ;
    Vdzxm1yz1pri = Array<OneD, NekDouble>(totszedges1d);
    ;
    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
    edgeptsin[1] = Array<OneD, NekDouble>(edgexyztemp);
    Vxm1yz1pri   = Epri->PhysEvaluateBasis(
          edgeptsin, demo.storage3dpri, Vdxxm1yz1pri, Vdyxm1yz1pri, Vdzxm1yz1pri);

    // edge bot  right (x = 1) (z = -1)
    Vx1yzm1pri   = Array<OneD, NekDouble>(totszedges1d);
    Vdxx1yzm1pri = Array<OneD, NekDouble>(totszedges1d);
    Vdyx1yzm1pri = Array<OneD, NekDouble>(totszedges1d);
    Vdzx1yzm1pri = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
    edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    Vx1yzm1pri   = Epri->PhysEvaluateBasis(
          edgeptsin, demo.storage3dpri, Vdxx1yzm1pri, Vdyx1yzm1pri, Vdzx1yzm1pri);

    // edge bot left (x = -1) (z = -1)
    Vxm1yzm1pri   = Array<OneD, NekDouble>(totszedges1d);
    Vdxxm1yzm1pri = Array<OneD, NekDouble>(totszedges1d);
    Vdyxm1yzm1pri = Array<OneD, NekDouble>(totszedges1d);
    Vdzxm1yzm1pri = Array<OneD, NekDouble>(totszedges1d);
    edgeptsin[0]  = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    Vxm1yzm1pri =
        Epri->PhysEvaluateBasis(edgeptsin, demo.storage3dpri, Vdxxm1yzm1pri,
                                Vdyxm1yzm1pri, Vdzxm1yzm1pri);

    int totsurf2d = (demo.coordquad[0].size()) * Epri->GetNcoeffs();

    Array<OneD, Array<OneD, NekDouble>> surfptsin(dimension),
        surfptsintemp(dimension);
    for (int k = 0; k < dimension - 1; k++)
    {
        surfptsin[k]     = Array<OneD, NekDouble>(demo.coordquad[k]);
        surfptsintemp[k] = Array<OneD, NekDouble>(demo.coordquad[k]);
    }

    surfptsin[dimension - 1] = Array<OneD, NekDouble>(demo.coordquad[0].size());
    surfptsintemp[dimension - 1] =
        Array<OneD, NekDouble>(demo.coordquad[0].size());

    // surface bot z = -1, (ABC)
    Vxyzm1pri        = Array<OneD, NekDouble>(totsurf2d);
    surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
    surfptsintemp[0] = surfptsin[0];
    surfptsintemp[1] = surfptsin[1];
    Vxyzm1pri        = Epri->PhysEvaluateBasis(
               surfptsintemp, demo.storage3dpri, NullNekDouble1DArray,
               NullNekDouble1DArray, NullNekDouble1DArray);

    int totpt = surfptsin[0].size();

    // surf hypt x+z <=0
    Vxemzym1pri = Array<OneD, NekDouble>(totsurf2d);
    Vmath::Smul(totpt, -1.0, surfptsintemp[0], 1, surfptsintemp[2], 1);
    Vxemzym1pri = Epri->PhysEvaluateBasis(
        surfptsintemp, demo.storage3dpri, NullNekDouble1DArray,
        NullNekDouble1DArray, NullNekDouble1DArray);

    // left x = -1
    Vxm1yzpri        = Array<OneD, NekDouble>(totsurf2d);
    surfptsintemp[2] = surfptsin[2];
    surfptsintemp[0] = Array<OneD, NekDouble>(totpt, -1.0);
    Vxm1yzpri        = Epri->PhysEvaluateBasis(
               surfptsintemp, demo.storage3dpri, NullNekDouble1DArray,
               NullNekDouble1DArray, NullNekDouble1DArray);

    // rest surfaces are triangles:
    totsurf2d = (demo.coordtri[0].size()) * Epri->GetNcoeffs();

    for (int k = 0; k < dimension - 1; k++)
    {
        surfptsin[k]     = Array<OneD, NekDouble>(demo.coordtri[k]);
        surfptsintemp[k] = Array<OneD, NekDouble>(demo.coordtri[k]);
    }
    surfptsin[dimension - 1] = Array<OneD, NekDouble>(demo.coordtri[0].size());
    surfptsintemp[dimension - 1] =
        Array<OneD, NekDouble>(demo.coordtri[0].size());
    totpt = surfptsin[0].size();

    // surface front  x+z <= 0, y = -1

    Vxym1zpri = Array<OneD, NekDouble>(totsurf2d);
    Vmath::Smul(totpt, -1.0, &surfptsintemp[0][0], 1, &surfptsintemp[2][0], 1);
    surfptsintemp[1] = Array<OneD, NekDouble>(totpt, -1.0);
    Vxym1zpri        = Epri->PhysEvaluateBasis(
               surfptsintemp, demo.storage3dpri, NullNekDouble1DArray,
               NullNekDouble1DArray, NullNekDouble1DArray);

    // back y=1, x+z <= 0
    Vxy1zpri         = Array<OneD, NekDouble>(totsurf2d);
    surfptsintemp[1] = Array<OneD, NekDouble>(totpt, 1.0);
    Vxy1zpri         = Epri->PhysEvaluateBasis(
                surfptsintemp, demo.storage3dpri, NullNekDouble1DArray,
                NullNekDouble1DArray, NullNekDouble1DArray);
}

void FilterStructurePres::call_setup_pyr(LocalRegions::ExpansionSharedPtr exp)
{
    dimension   = 3;
    int nmodes0 = exp->GetBasis(0)->GetNumModes();
    int nmodes1 = exp->GetBasis(1)->GetNumModes();
    int nmodes2 = exp->GetBasis(2)->GetNumModes();
    int npts0   = exp->GetBasis(0)->GetNumPoints();
    int npts1   = exp->GetBasis(1)->GetNumPoints();
    int npts2   = exp->GetBasis(2)->GetNumPoints();
    PointsKey p0(npts0, exp->GetBasis(0)->GetPointsType());
    PointsKey p1(npts1, exp->GetBasis(1)->GetPointsType());
    PointsKey p2(npts2, exp->GetBasis(2)->GetPointsType());
    BasisKey b0(LibUtilities::eOrtho_A, nmodes0, p0);
    BasisKey b1(LibUtilities::eOrtho_A, nmodes1, p1);
    BasisKey b2(LibUtilities::eOrthoPyr_C, nmodes2, p2);
    Epyr              = new StdPyrExp(b0, b1, b2);
    demo.storage3dpyr = Epyr->GetPhysEvaluateStorage();
    demo.coordpyr     = demo.GetCoords(Epyr);
    demo.coordmidpyr  = demo.GetQuadratureMidCoords(demo.coordpyr);
    demo.coordlatticepyr =
        demo.GetLatticeCoords(demo.coordpyr, demo.coordmidpyr);

    demo.midptevalpyr = Epyr->PhysEvaluateBasis(
        demo.coordmidpyr, demo.storage3dpyr, NullNekDouble1DArray,
        NullNekDouble1DArray, NullNekDouble1DArray);

    // we need qzin + qzinmid as input points
    Array<OneD, NekDouble> edgexyztemp = E3seg->GetBasis(0)->GetZ();
    int totszedges1d = edgexyztemp.size() * (Epyr->GetNcoeffs());
    Array<OneD, Array<OneD, NekDouble>> edgeptsin(dimension);
    for (int p = 0; p < dimension; p++)
    {
        edgeptsin[p] = Array<OneD, NekDouble>(edgexyztemp);
    }

    // edge front left EA (x = -1) (y = -1)
    Vxm1ym1zpyr   = Array<OneD, NekDouble>(totszedges1d);
    Vdxxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);
    Vdyxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);
    Vdzxm1ym1zpyr = Array<OneD, NekDouble>(totszedges1d);

    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    Vxm1ym1zpyr =
        Epyr->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxxm1ym1zpyr,
                                Vdyxm1ym1zpyr, Vdzxm1ym1zpyr);

    // edge front hypt EB (y = -1) (z + x = 0)
    Vym1xmzpyr   = Array<OneD, NekDouble>(totszedges1d);
    Vdxym1xmzpyr = Array<OneD, NekDouble>(totszedges1d);
    Vdyym1xmzpyr = Array<OneD, NekDouble>(totszedges1d);
    Vdzym1xmzpyr = Array<OneD, NekDouble>(totszedges1d);
    Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[2][0], 1,
                &edgeptsin[0][0], 1);

    Vxm1ym1zpyr = Epyr->PhysEvaluateBasis(
        edgeptsin, demo.storage3dpyr, Vdxym1xmzpyr, Vdxym1xmzpyr, Vdzym1xmzpyr);
    // edge back hypt (z = -x  or x = -z) (x = y) EC
    Vxeyxmzpyr   = Array<OneD, NekDouble>(totszedges1d);
    Vdxxeyxmzpyr = Array<OneD, NekDouble>(totszedges1d);
    Vdyxeyxmzpyr = Array<OneD, NekDouble>(totszedges1d);
    Vdzxeyxmzpyr = Array<OneD, NekDouble>(totszedges1d);

    edgeptsin[1] = Array<OneD, NekDouble>(edgexyztemp);
    Vxeyxmzpyr   = Epyr->PhysEvaluateBasis(
          edgeptsin, demo.storage3dpyr, Vdxxeyxmzpyr, Vdyxeyxmzpyr, Vdzxeyxmzpyr);

    // edge back left (y = -z) (x = -1) ED
    Vx1ymzpyr   = Array<OneD, NekDouble>(totszedges1d);
    Vdxx1ymzpyr = Array<OneD, NekDouble>(totszedges1d);
    Vdyx1ymzpyr = Array<OneD, NekDouble>(totszedges1d);
    Vdzx1ymzpyr = Array<OneD, NekDouble>(totszedges1d);
    ;
    Vmath::Smul(edgeptsin[0].size(), -1.0, &edgeptsin[2][0], 1,
                &edgeptsin[1][0], 1);

    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    Vx1ymzpyr    = Epyr->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr,
                                           Vdxx1ymzpyr, Vdyx1ymzpyr, Vdzx1ymzpyr);

    // edge front bot (y = -1) (z = -1) AB
    Vym1xzm1pyr   = Array<OneD, NekDouble>(totszedges1d);
    Vdxym1xzm1pyr = Array<OneD, NekDouble>(totszedges1d);
    Vdyym1xzm1pyr = Array<OneD, NekDouble>(totszedges1d);
    Vdzym1xzm1pyr = Array<OneD, NekDouble>(totszedges1d);

    edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    edgeptsin[0] = Array<OneD, NekDouble>(edgexyztemp);
    Vym1xzm1pyr =
        Epyr->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxym1xzm1pyr,
                                Vdyym1xzm1pyr, Vdzym1xzm1pyr);

    // edge back bot (y = 1) (z = -1)) DC
    Vy1xzm1pyr   = Array<OneD, NekDouble>(totszedges1d);
    Vdyy1xzm1pyr = Array<OneD, NekDouble>(totszedges1d);
    Vdxy1xzm1pyr = Array<OneD, NekDouble>(totszedges1d);
    Vdzy1xzm1pyr = Array<OneD, NekDouble>(totszedges1d);

    edgeptsin[0] = edgexyztemp;
    edgeptsin[1] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
    Vy1xzm1pyr   = Epyr->PhysEvaluateBasis(
          edgeptsin, demo.storage3dpyr, Vdxy1xzm1pyr, Vdyy1xzm1pyr, Vdzy1xzm1pyr);

    // edge left bot (z = -1), (x = -1) AD
    Vxm1yzm1pyr   = Array<OneD, NekDouble>(totszedges1d);
    Vdyxm1yzm1pyr = Array<OneD, NekDouble>(totszedges1d);
    Vdxxm1yzm1pyr = Array<OneD, NekDouble>(totszedges1d);
    Vdzxm1yzm1pyr = Array<OneD, NekDouble>(totszedges1d);

    edgeptsin[1] = edgexyztemp;
    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    Vxm1yzm1pyr =
        Epyr->PhysEvaluateBasis(edgeptsin, demo.storage3dpyr, Vdxxm1yzm1pyr,
                                Vdyxm1yzm1pyr, Vdzxm1yzm1pyr);

    // edge right bot ( z = -1) (x = 1) BC
    Vx1yzm1pyr   = Array<OneD, NekDouble>(totszedges1d);
    Vdyx1yzm1pyr = Array<OneD, NekDouble>(totszedges1d);
    Vdxx1yzm1pyr = Array<OneD, NekDouble>(totszedges1d);
    Vdzx1yzm1pyr = Array<OneD, NekDouble>(totszedges1d);

    edgeptsin[1] = edgexyztemp;
    edgeptsin[0] = Array<OneD, NekDouble>(edgeptsin[0].size(), 1.0);
    edgeptsin[2] = Array<OneD, NekDouble>(edgeptsin[0].size(), -1.0);
    Vx1yzm1pyr   = Epyr->PhysEvaluateBasis(
          edgeptsin, demo.storage3dpyr, Vdxx1yzm1pyr, Vdyx1yzm1pyr, Vdzx1yzm1pyr);

    int totsurf2d = (demo.coordquad[0].size()) * Epyr->GetNcoeffs();

    Array<OneD, Array<OneD, NekDouble>> surfptsin(dimension),
        surfptsintemp(dimension);
    for (int k = 0; k < dimension - 1; k++)
    {
        surfptsin[k]     = Array<OneD, NekDouble>(demo.coordquad[k]);
        surfptsintemp[k] = Array<OneD, NekDouble>(demo.coordquad[k]);
    }

    surfptsin[dimension - 1] = Array<OneD, NekDouble>(demo.coordquad[0].size());
    surfptsintemp[dimension - 1] =
        Array<OneD, NekDouble>(demo.coordquad[0].size());

    // surface bot z = -1, (ABC)
    Vxyzm1pyr        = Array<OneD, NekDouble>(totsurf2d);
    surfptsintemp[2] = Array<OneD, NekDouble>(surfptsin[0].size(), -1.0);
    surfptsintemp[0] = surfptsin[0];
    surfptsintemp[1] = surfptsin[1];
    Vxyzm1pyr        = Epyr->PhysEvaluateBasis(
               surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray,
               NullNekDouble1DArray, NullNekDouble1DArray);

    totsurf2d = (demo.coordtri[0].size()) * Epyr->GetNcoeffs();

    for (int k = 0; k < dimension - 1; k++)
    {
        surfptsin[k]     = Array<OneD, NekDouble>(demo.coordtri[k]);
        surfptsintemp[k] = Array<OneD, NekDouble>(demo.coordtri[k]);
    }
    surfptsin[dimension - 1] = Array<OneD, NekDouble>(demo.coordtri[0].size());
    surfptsintemp[dimension - 1] =
        Array<OneD, NekDouble>(demo.coordtri[0].size());

    int totpt = surfptsin[0].size();

    // surface hypt (tri) x+z < 0,
    Vxmzypyr = Array<OneD, NekDouble>(totsurf2d);
    Vmath::Smul(totpt, -1.0, &surfptsintemp[0][0], 1, &surfptsintemp[2][0], 1);
    Vmath::Vcopy(totpt, &surfptsintemp[0][0], 1, &surfptsintemp[1][0], 1);
    Vxmzypyr = Epyr->PhysEvaluateBasis(
        surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray,
        NullNekDouble1DArray, NullNekDouble1DArray);

    // surface left (tri) x = -1, (y+z < 0)
    Vxm1ymzpyr       = Array<OneD, NekDouble>(totsurf2d);
    surfptsintemp[1] = surfptsin[1];
    surfptsintemp[2] = surfptsin[2];
    surfptsintemp[0] = Array<OneD, NekDouble>(totpt, -1.0);
    Vxm1ymzpyr       = Epyr->PhysEvaluateBasis(
              surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray,
              NullNekDouble1DArray, NullNekDouble1DArray);

    // surface front (tri) y = -1, x+z < 0
    Vxmzym1pyr       = Array<OneD, NekDouble>(totsurf2d);
    surfptsintemp[1] = Array<OneD, NekDouble>(totpt, -1.0);
    surfptsintemp[0] = surfptsin[0];
    Vxmzym1pyr       = Epyr->PhysEvaluateBasis(
              surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray,
              NullNekDouble1DArray, NullNekDouble1DArray);

    // surface back (tri) y = -z
    Vxymzpyr = Array<OneD, NekDouble>(totsurf2d);
    Vmath::Smul(totpt, -1.0, &surfptsintemp[2][0], 1, &surfptsintemp[1][0], 1);
    Vxymzpyr = Epyr->PhysEvaluateBasis(
        surfptsintemp, demo.storage3dpyr, NullNekDouble1DArray,
        NullNekDouble1DArray, NullNekDouble1DArray);
}

Array<OneD, Array<OneD, NekDouble>> FilterStructurePres::FindIds(
    Array<OneD, NekDouble> &uhatslocal, StdExpansion *Eorth, NekDouble tol)
{
    boost::ignore_unused(tol);
    Array<OneD, Array<OneD, NekDouble>> minvalandpthold;

    // // check on mid-quadrature pts of staggered lattice
    switch (Eorth->DetShapeType())
    {
        case LibUtilities::eQuadrilateral:
            minvalandpthold = demo.FindLatticeEval(
                uhatslocal, demo.storage2dq, demo.coordquad, demo.midptevalquad,
                demo.coordmidquad);
            break;
        case LibUtilities::eTriangle:
            minvalandpthold =
                demo.FindLatticeEval(uhatslocal, demo.storage2dt, demo.coordtri,
                                     demo.midptevaltri, demo.coordmidtri);
            break;
        case LibUtilities::eHexahedron:
            minvalandpthold = demo.FindLatticeEval(
                uhatslocal, demo.storage3dhex, demo.coordhex, demo.midptevalhex,
                demo.coordmidhex);
            break;
        case LibUtilities::eTetrahedron:
            minvalandpthold = demo.FindLatticeEval(
                uhatslocal, demo.storage3dtet, demo.coordtet, demo.midptevaltet,
                demo.coordmidtet);
            break;
        case LibUtilities::ePyramid:
            minvalandpthold = demo.FindLatticeEval(
                uhatslocal, demo.storage3dpyr, demo.coordpyr, demo.midptevalpyr,
                demo.coordmidpyr);
            break;
        case LibUtilities::ePrism:
            minvalandpthold = demo.FindLatticeEval(
                uhatslocal, demo.storage3dpri, demo.coordpri, demo.midptevalpri,
                demo.coordmidpri);
            break;
        default:
            cout << "\n not implemented for this shape type yet!";
            exit(0);
    }

    int sz = minvalandpthold[0].size();

    Array<OneD, Array<OneD, NekDouble>> retarr(2);
    //
    if (demo.details)
    {
        if (minvalandpthold[1][dimension] > 1 + 1e-9 && demo.ub == 1)
        {
            cout << " max = " << minvalandpthold[1][dimension] << " at ("
                 << minvalandpthold[1][0] << "," << minvalandpthold[1][1] << ")"
                 << " \n*****\n";
        }
    }

    if ((minvalandpthold[0][sz - 1] < 0 &&
         abs(minvalandpthold[0][sz - 1]) > 1e-6) ||
        (minvalandpthold[1][sz - 1] > 1 && demo.ub == 1))
    {
        if ((minvalandpthold[0][sz - 1] < 0 &&
             abs(minvalandpthold[0][sz - 1]) > 1e-6))
        {
            retarr[0] = minvalandpthold[0];
        }

        if (demo.ub == 1 && minvalandpthold[1][sz - 1] > 1 &&
            abs(minvalandpthold[1][sz - 1] - 1) > 1e-6)
        {
            retarr[1] = minvalandpthold[1];
        }

        return retarr;
    }
    return NullNekDoubleArrayOfArray;
}

void FilterStructurePres::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    Timer t;
    boost::ignore_unused(time);
    NekDouble avgiterall = 0;
    for (int i = 0; i < nfields; i++)

    {
        if (demo.details)
        {
            cout << "\t field = " << i << "\n";
        }
        int ctrcoeff = 0, ctrphys = 0;

        int nelmt = pFields[i]->GetExpSize();

        Array<OneD, NekDouble> uhatsallel = pFields[i]->GetCoeffs();
        Array<OneD, Array<OneD, NekDouble>> flaggedele(
            2); // both upper and lower bnds

        Array<OneD, NekDouble> physallel = pFields[i]->GetPhys();

        Array<OneD, NekDouble> saveuhatsallel(uhatsallel.size());
        Array<OneD, NekDouble> saveuhatsallel2(uhatsallel.size());
        Array<OneD, NekDouble> savephysallel(physallel.size());
        Array<OneD, NekDouble> savephysallel2(physallel.size());

        Vmath::Vcopy(uhatsallel.size(), uhatsallel, 1, saveuhatsallel, 1);
        Vmath::Vcopy(uhatsallel.size(), uhatsallel, 1, saveuhatsallel2, 1);
        Vmath::Vcopy(physallel.size(), physallel, 1, savephysallel, 1);
        Vmath::Vcopy(physallel.size(), physallel, 1, savephysallel2, 1);
        // int total_flagged = 0;
        int timestep_flagged = 0;
        for (int k = 0; k < nelmt; k++)
        {

            Array<OneD, NekDouble> uhatslocal(
                pFields[i]->GetExp(k)->GetNcoeffs());
            Array<OneD, NekDouble> physlocal(
                pFields[i]->GetExp(k)->GetTotPoints());
            Vmath::Vcopy(physlocal.size(), &physallel[ctrphys], 1,
                         &physlocal[0], 1);

            t.Start();

            listeletype[k]->FwdTrans(physlocal, uhatslocal);

            t.Stop();

            timeOrth += t.TimePerTest(1);
            t.Start();
            flaggedele = FindIds(uhatslocal, listeletype[k], demo.GettolGD());
            t.Stop();
            timeFindIds += t.TimePerTest(1);

            if (flaggedele != NullNekDoubleArrayOfArray)
            {
	        timestep_flagged++;
                total_flagged++;
                
                if (demo.details)
                {
                    if (flaggedele[1].size() > 0)
                    {
                        cout << "\n flaggedele max= " << k << " v(x,y)= ("
                             << flaggedele[1][0] << ", " << flaggedele[1][1]
                             << ", "
                             << ") val = " << flaggedele[1][2];
                    }
                }
                NekDouble retavgiter = 0, basiscallshold = 0;
                t.Start();
                Optimize(listeletype[k], uhatslocal, flaggedele, retavgiter,
                         basiscallshold);
                t.Stop();
                NekDouble tmptime = t.TimePerTest(1);
                
                timeStrPres += tmptime;
                retGDall += retavgiter;
                basiscalls += basiscallshold;

                avgiterall += retavgiter;

                flaggedele =
                    FindIds(uhatslocal, listeletype[k], demo.GettolGD());

                timeFindIds += t.TimePerTest(1);
                if (flaggedele.size() > 0 &&
                    (flaggedele[0].size() > 0 || flaggedele[1].size() > 0))
                {
                    cout << "\nStructure Preservation failed at : \n";
                    cout << "\n elid = " << k << " ";
                    for (int p = 0; p < flaggedele.size(); p++)
                    {
                        cout << "\n val" << p << "=";

                        for (int j = 0; j < flaggedele[p].size(); j++)
                        {
                            cout << flaggedele[p][j] << " ";
                        }
                    }

                    exit(0);
                }
                listeletype[k]->BwdTrans(uhatslocal, physlocal);
            }
            else
            {
            }
            Vmath::Vcopy(physlocal.size(), &physlocal[0], 1,
                         &savephysallel[ctrphys], 1);

            ctrcoeff += uhatslocal.size();

            ctrphys += physlocal.size();
        }
        pFields[i]->SetPhys(savephysallel);
    }
}

void FilterStructurePres::Optimize(
    StdExpansion *exp, Array<OneD, NekDouble> &coeffs,
    Array<OneD, Array<OneD, NekDouble>> flaggedelecoord,
    NekDouble &avgiterGDret, NekDouble &basiscallstotal)
{
    boost::ignore_unused(flaggedelecoord);

    Timer t1, t, totroottimer;
    ;
    t1.Start();

    StdExpansion *E = exp;
    int ns, dim; 
    int stype = exp->DetShapeType();
    Array<OneD, Array<OneD, NekDouble>> surfaceuhats, storage, Pf, tmpcoord;
    int N1 = coeffs.size();

    switch (stype)
    {
        case LibUtilities::eQuadrilateral:

            ns      = 0;
            dim     = 2;
            storage = demo.storage2dq;

            break;
        case LibUtilities::eTriangle:

            ns      = 0;
            dim     = 2;
            storage = demo.storage2dt;

            break;
        case LibUtilities::eTetrahedron:
            ns           = 4;
            dim          = 3;
            storage      = demo.storage3dtet;
            surfaceuhats = Array<OneD, Array<OneD, NekDouble>>(ns);
            for (int k = 0; k < ns; k++)
            {
                surfaceuhats[k] = Array<OneD, NekDouble>(Etri->GetNcoeffs());
            }

            break;
        case LibUtilities::eHexahedron:
            ns           = 6;
            dim          = 3;
            storage      = demo.storage3dhex;
            surfaceuhats = Array<OneD, Array<OneD, NekDouble>>(ns);
            break;
        case LibUtilities::ePyramid:
            ns              = 5;
            dim             = 3;
            storage         = demo.storage3dpyr;
            surfaceuhats    = Array<OneD, Array<OneD, NekDouble>>(ns);
            surfaceuhats[0] = Array<OneD, NekDouble>(Equad->GetNcoeffs());
            for (int k = 1; k < ns; k++)
            {
                surfaceuhats[k] = Array<OneD, NekDouble>(Etri->GetNcoeffs());
            }
            break;
        case LibUtilities::ePrism:
            ns           = 5;
            dim          = 3;
            storage      = demo.storage3dpri;
            surfaceuhats = Array<OneD, Array<OneD, NekDouble>>(ns);

            // surf seq = base, front, right, back, left, rop

            surfaceuhats[0] = Array<OneD, NekDouble>(Equad->GetNcoeffs());
            surfaceuhats[1] = Array<OneD, NekDouble>(Etri->GetNcoeffs());
            surfaceuhats[2] = Array<OneD, NekDouble>(Equad->GetNcoeffs());
            surfaceuhats[3] = Array<OneD, NekDouble>(Equad->GetNcoeffs());
            surfaceuhats[4] = Array<OneD, NekDouble>(Etri->GetNcoeffs());

            break;

        default:
            cout << "\n element type not supported yet";
            exit(0);
    }
    double inf = numeric_limits<double>::infinity();

    int niter = 1e3, counter = 0;
    NekDouble tol  = demo.GettolGD(),
              minv = inf; // constraint specific tolerances

    Array<OneD, NekDouble> optima(dim), pqvalcoords(dim + 1), xastarr(dim),
        utemp(N1), wsp1(1), Vtmp(N1);
    tmpcoord = Array<OneD, Array<OneD, NekDouble>>(dim);

    NekDouble pqval; // timeprojectedges = 0.0;//, timeprojectsurf = 0.0;

    NekDouble avgiterGD = 0.0, basiscalls1 = 0.0;

    NekDouble startcoordx, startcoordy, startcoordz;
    ;

    m_session->LoadParameter("UpperBound", demo.ub, inf);
    m_session->LoadParameter("details", demo.details, 0);

    // by default lower bound specified and always 0
    // Change this for more constraints
    //int numConstraints = 1;
    Array<OneD, NekDouble> fac(2);
    Array<OneD, NekDouble> cc1(N1, 0.0), hold(N1);
    vector<Array<OneD, NekDouble>> constrcoeffs;

    // change this if lower bound is not 0
    constrcoeffs.push_back(cc1);
    // nu for lower bound
    fac[0] = 1;
    if (demo.ub != inf)
    {
        // nu for upper bound
        fac[1]         = -1;
        // numConstraints = 2;
        // considering upperbound as 1: project fn of all 1s on E

        Array<OneD, NekDouble> fub(E->GetTotPoints(), 1.0);
        Array<OneD, NekDouble> cc2(N1);
        E->FwdTrans(fub, cc2);
        constrcoeffs.push_back(cc2);
    }
    int idxsave = 0, j = 0;
    while (counter <= niter)
    {
        utemp = coeffs;

        NekDouble roots1dtime = 0.0, roots2dtime = 0.0, roots3dtime = 0.0;
        pqval                 = inf;
        NekDouble avgiterhold = 0, callstobasishold = 0.0;
        ;
        idxsave = -1;
        //	  for( int j = 0; j < numConstraints; j++)
        //{
        Vmath::Vsub(N1, utemp, 1, cc1, 1, hold, 1);
        if (counter > 0)
        {

            avgiterhold = 0;
            totroottimer.Start();
            optima = call_find_roots(hold, avgiterhold, Pf, surfaceuhats, minv,
                                     E, 1e-7, roots1dtime, roots2dtime,
                                     roots3dtime, callstobasishold, fac[j]);
            totroottimer.Stop();
            totalroottime += totroottimer.TimePerTest(1);

            avgiterGDret += avgiterhold;
            basiscalls1 += callstobasishold;
            roots1dtimehold += roots1dtime;
            roots2dtimehold += roots2dtime;
            roots3dtimehold += roots3dtime;
            avgiterGD += avgiterhold;

            avgiterhold = 0;
            for (int k = 0; k < dim; k++)
            {
                tmpcoord[k] = Array<OneD, NekDouble>(1, optima[k]);
            }
            totroottimer.Start();
            Vtmp = E->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray,
                                        NullNekDouble1DArray,
                                        NullNekDouble1DArray);
            tot_interp_points += tmpcoord[0].size();

            tot_interp_points += tmpcoord[0].size();
            totroottimer.Stop();
            basiscalls1 += totroottimer.TimePerTest(1);
        }
        else
	{
            if (j == 0 || (j == 1 && flaggedelecoord.size() > 1))
            {
                startcoordx = flaggedelecoord[j][0];
                if (dim > 1)
                {
                    startcoordy = flaggedelecoord[j][1];
                }
                if (dim > 2)
                {
                    startcoordz = flaggedelecoord[j][2];
                }

                for (int k = 0; k < dim; k++)
                {
                    tmpcoord[k] = Array<OneD, NekDouble>(1);
                }

                tmpcoord[0][0] = startcoordx;

                if (dim > 1)
                    tmpcoord[1][0] = startcoordy;

                if (dim > 2)
                    tmpcoord[2][0] = startcoordz;

                optima[0] = startcoordx;
                if (dim > 1)
                    optima[1] = startcoordy;
                if (dim > 2)
                    optima[2] = startcoordz;
                totroottimer.Start();
                Vtmp = E->PhysEvaluateBasis(
                    tmpcoord, storage, NullNekDouble1DArray,
                    NullNekDouble1DArray, NullNekDouble1DArray);
                tot_interp_points += tmpcoord[0].size();

                tot_interp_points += tmpcoord[0].size();
                totroottimer.Stop();
                basiscalls1 += totroottimer.TimePerTest(1);
            }
        }

        demo.pq(hold, tmpcoord, Vtmp, pqvalcoords, wsp1, fac[j]);
        if (counter == 0 && abs(pqvalcoords[0]) < tol)
        {
            tol = abs(pqvalcoords[0]);
        }
        if (demo.details)
        {

            cout << "\n constr number = 0 pqvalcoords[0] = " << pqvalcoords[0]
                 << " at " << tmpcoord[0][0] << "," << tmpcoord[1][0] << ", "
                 << tmpcoord[2][0] << "  iter=" << counter << "\n";
        }

        if (pqvalcoords[0] < pqval)
        {
            xastarr[0] = tmpcoord[0][0];
            if (dim > 1)
                xastarr[1] = tmpcoord[1][0];

            if (dim > 2)
                xastarr[2] = tmpcoord[2][0];
            pqval   = pqvalcoords[0];
            idxsave = j;
        }
        // loop on number of constr ends

        // If minimum is non-negative, we're done

        if (pqval > -tol)
        {
            break;
        }
        counter = counter + 1;

        Array<OneD, NekDouble> Vastsq(N1);
        Array<OneD, NekDouble> Vast(N1);
        Array<OneD, Array<OneD, NekDouble>> xastarrofarr(dim);
        for (int k = 0; k < dim; k++)
        {
            xastarrofarr[k] = Array<OneD, NekDouble>(1, xastarr[k]);
        }
        Array<OneD, NekDouble> tmp;
        NekDouble vastsqsum;
        // basiscalls++;
        for (int i = 0; i < N1; i++)
        {
            totroottimer.Start();
            tmp = E->PhysEvaluateBasis(xastarrofarr, storage, i);
            tot_interp_points += 1;

            tot_interp_points += 1;
            totroottimer.Stop();
            basiscalls1 += totroottimer.TimePerTest(1);
            Vast[i]   = tmp[0];
            Vastsq[i] = (Vast[i] * Vast[i]);
        }

        vastsqsum = Vmath::Vsum(N1, &Vastsq[0], 1);
        Array<OneD, NekDouble> qast(N1);

        for (int i = 0; i < N1; i++)
        {
            qast[i] = ((1 / sqrt(vastsqsum)) * (Vast[i]));
        }

        Vmath::Smul(N1, fac[idxsave] * pqval, &qast[0], 1, &qast[0], 1);

        Vmath::Vsub(utemp.size(), utemp, 1, constrcoeffs[idxsave], 1, utemp, 1);
        Vmath::Vsub(utemp.size(), &utemp[0], 1, &qast[0], 1, &qast[0], 1);
        Vmath::Vadd(utemp.size(), qast, 1, constrcoeffs[idxsave], 1, utemp, 1);
    }
    coeffs = utemp;
    t1.Stop();
    roots1dtimehold = roots1dtimehold;
    roots2dtimehold = roots2dtimehold;
    roots3dtimehold = roots3dtimehold;

    tot_optiter += counter;
    basiscallstotal = basiscalls1;
    if (demo.details)
    {
        cout << " optimizer iterations = " << counter;
    }
}

void FilterStructurePres::project_edges(
    Array<OneD, NekDouble> uhats, Array<OneD, Array<OneD, NekDouble>> &ret,
    StdExpansion *E, NekDouble nu)
{
    if (E->DetShapeType() == LibUtilities::eQuadrilateral) // quad
    {
        // bot edge
        edgederpquhats(uhats, ret[0], E->GetNcoeffs(), Vym1q, Vdxym1q, Vdyym1q,
                       NullNekDouble1DArray, nu);
        // right edge
        edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vx1q, Vdxx1q, Vdyx1q,
                       NullNekDouble1DArray, nu);
        // top edge
        edgederpquhats(uhats, ret[2], E->GetNcoeffs(), Vy1q, Vdxy1q, Vdyy1q,
                       NullNekDouble1DArray, nu);

        // left edge
        edgederpquhats(uhats, ret[3], E->GetNcoeffs(), Vxm1q, Vdxxm1q, Vdyxm1q,
                       NullNekDouble1DArray, nu);

        return;
    }
    if (E->DetShapeType() == LibUtilities::eTriangle) // tri
    {

        // bot edge
        edgederpquhats(uhats, ret[0], E->GetNcoeffs(), Vym1t, Vdxym1t, Vdyym1t,
                       NullNekDouble1DArray, nu);
        // left edge
        edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vxm1t, Vdxxm1t, Vdyxm1t,
                       NullNekDouble1DArray, nu);

        // hypto
        edgederpquhats(uhats, ret[2], E->GetNcoeffs(), Vxyhypt, Vdxxyhypt,
                       Vdyxyhypt, NullNekDouble1DArray, nu);

        return;
    }

    if (E->DetShapeType() == LibUtilities::eHexahedron) // hex
    {

        // edge front left (x = -1) (y = -1)
        edgederpquhats(uhats, ret[0], E->GetNcoeffs(), Vxm1ym1z, Vdxxm1ym1z,
                       Vdyxm1ym1z, Vdzxm1ym1z, nu);
        // edge front right (x = 1) (y = -1)
        edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vx1ym1z, Vdxx1ym1z,
                       Vdyx1ym1z, Vdzx1ym1z, nu);

        // edge front top (y = -1) (z = 1)
        edgederpquhats(uhats, ret[2], E->GetNcoeffs(), Vym1xz1, Vdxym1xz1,
                       Vdyym1xz1, Vdzym1xz1, nu);

        // edge front bot (y = -1) (z = -1)
        edgederpquhats(uhats, ret[3], E->GetNcoeffs(), Vym1xzm1, Vdxym1xzm1,
                       Vdyym1xzm1, Vdzym1xzm1, nu);

        // edge back left (y = 1), (x = -1)
        edgederpquhats(uhats, ret[4], E->GetNcoeffs(), Vxm1y1z, Vdxxm1y1z,
                       Vdyxm1y1z, Vdzxm1y1z, nu);

        // edge back right (x = 1), (y = 1))
        edgederpquhats(uhats, ret[5], E->GetNcoeffs(), Vx1y1z, Vdxx1y1z,
                       Vdyx1y1z, Vdzx1y1z, nu);

        // edge back top ( y = 1) (z = 1)
        edgederpquhats(uhats, ret[6], E->GetNcoeffs(), Vy1xz1, Vdxy1xz1,
                       Vdyy1xz1, Vdzy1xz1, nu);

        // edge back bot (y = 1) (z = -1))
        edgederpquhats(uhats, ret[7], E->GetNcoeffs(), Vy1xzm1, Vdxy1xzm1,
                       Vdyy1xzm1, Vdzy1xzm1, nu);

        // edge left bot (z = -1), (x = -1)
        edgederpquhats(uhats, ret[8], E->GetNcoeffs(), Vxm1yzm1, Vdxxm1yzm1,
                       Vdyxm1yzm1, Vdzxm1yzm1, nu);

        // edge left top (x = -1), (z = 1))
        edgederpquhats(uhats, ret[9], E->GetNcoeffs(), Vxm1yz1, Vdxxm1yz1,
                       Vdyxm1yz1, Vdzxm1yz1, nu);

        // edge right bot ( z = -1) (x = 1)
        edgederpquhats(uhats, ret[10], E->GetNcoeffs(), Vx1yzm1, Vdxx1yzm1,
                       Vdyx1yzm1, Vdzx1yzm1, nu);

        // edge right top (z  1) (x  1))
        edgederpquhats(uhats, ret[11], E->GetNcoeffs(), Vx1yz1, Vdxx1yz1,
                       Vdyx1yz1, Vdzx1yz1, nu);
    }
    else if (E->DetShapeType() == LibUtilities::eTetrahedron) // tet
    {
        // edge front left (AD) (x = -1) (y = -1)
        edgederpquhats(uhats, ret[0], E->GetNcoeffs(), Vxm1ym1ztet,
                       Vdxxm1ym1ztet, Vdyxm1ym1ztet, Vdzxm1ym1ztet, nu);

        // edge front hypt (DB) (y = -1) (z = -x)
        edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vym1xmztet, Vdxym1xmztet,
                       Vdyym1xmztet, Vdzym1xmztet, nu);

        // edge front bot (AB) (y = -1) (z = -1)
        edgederpquhats(uhats, ret[2], E->GetNcoeffs(), Vym1xzm1tet,
                       Vdxym1xzm1tet, Vdyym1xzm1tet, Vdzym1xzm1tet, nu);

        // edge left hypt (DC) ( x = -1) (z = -y)
        edgederpquhats(uhats, ret[3], E->GetNcoeffs(), Vxm1ymztet, Vdxxm1ymztet,
                       Vdyxm1ymztet, Vdzxm1ymztet, nu);

        // edge bot diag (BC) (z = -1) (y = -x)
        edgederpquhats(uhats, ret[4], E->GetNcoeffs(), Vxmyzm1tet, Vdxxmyzm1tet,
                       Vdyxmyzm1tet, Vdzxmyzm1tet, nu);

        // edge CA bot left (x = -1) (z = -1)
        edgederpquhats(uhats, ret[5], E->GetNcoeffs(), Vxm1yzm1tet,
                       Vdxxm1yzm1tet, Vdyxm1yzm1tet, Vdzxm1yzm1tet, nu);
    }
    else if (E->DetShapeType() == LibUtilities::ePrism) // prism
    {
        // edge front left (x = -1) (y = -1) EA
        edgederpquhats(uhats, ret[0], E->GetNcoeffs(), Vxm1ym1zpri,
                       Vdxxm1ym1zpri, Vdyxm1ym1zpri, Vdzxm1ym1zpri, nu);

        // edge front hypt (y = -1) (z = -x) EB
        edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vym1xmzpri, Vdxym1xmzpri,
                       Vdyym1xmzpri, Vdzym1xmzpri, nu);

        // edge front bot (y = -1) (z = -1) AB
        edgederpquhats(uhats, ret[2], E->GetNcoeffs(), Vym1xzm1pri,
                       Vdxym1xzm1pri, Vdyym1xzm1pri, Vdzym1xzm1pri, nu);

        // edge back left (y = 1) (x = -1))
        edgederpquhats(uhats, ret[3], E->GetNcoeffs(), Vxm1y1zpri, Vdxxm1y1zpri,
                       Vdyxm1y1zpri, Vdzxm1y1zpri, nu);

        // edge back hypt (y = 1) (z = -x)
        edgederpquhats(uhats, ret[4], E->GetNcoeffs(), Vy1xmzpri, Vdxy1xmzpri,
                       Vdyy1xmzpri, Vdzy1xmzpri, nu);
        // edge back bot (y = 1) (z = -1)) DC
        edgederpquhats(uhats, ret[5], E->GetNcoeffs(), Vy1xzm1pri, Vdxy1xzm1pri,
                       Vdyy1xzm1pri, Vdzy1xzm1pri, nu);

        // edge left top (x = -1) (z = 1)
        edgederpquhats(uhats, ret[6], E->GetNcoeffs(), Vxm1yz1pri, Vdxxm1yz1pri,
                       Vdyxm1yz1pri, Vdzxm1yz1pri, nu);

        // edge bot  right (x = 1) (z = -1)
        edgederpquhats(uhats, ret[7], E->GetNcoeffs(), Vx1yzm1pri, Vdxx1yzm1pri,
                       Vdyx1yzm1pri, Vdzx1yzm1pri, nu);

        // edge bot left (x = -1) (z = -1)
        edgederpquhats(uhats, ret[8], E->GetNcoeffs(), Vxm1yzm1pri,
                       Vdxxm1yzm1pri, Vdyxm1yzm1pri, Vdzxm1yzm1pri, nu);
    }
    else if (E->DetShapeType() == LibUtilities::ePyramid) // pyr
    {

        // edge front left (x = -1) (y = -1) EA
        edgederpquhats(uhats, ret[0], E->GetNcoeffs(), Vxm1ym1zpyr,
                       Vdxxm1ym1zpyr, Vdyxm1ym1zpyr, Vdzxm1ym1zpyr, nu);

        // edge front hypt (y = -1) (z = -x) EB
        edgederpquhats(uhats, ret[1], E->GetNcoeffs(), Vym1xmzpyr, Vdxym1xmzpyr,
                       Vdyym1xmzpyr, Vdzym1xmzpyr, nu);

        // edge back hypt (z = -x  or x = -z) (x = y) EC
        edgederpquhats(uhats, ret[2], E->GetNcoeffs(), Vxeyxmzpyr, Vdxxeyxmzpyr,
                       Vdyxeyxmzpyr, Vdzxeyxmzpyr, nu);

        // edge back left (y = -z) (x = -1) ED
        edgederpquhats(uhats, ret[3], E->GetNcoeffs(), Vxeyxmzpyr, Vdxxeyxmzpyr,
                       Vdyxeyxmzpyr, Vdzxeyxmzpyr, nu);

        // edge front bot (y = -1) (z = -1) AB
        edgederpquhats(uhats, ret[4], E->GetNcoeffs(), Vym1xzm1pyr,
                       Vdxym1xzm1pyr, Vdyym1xzm1pyr, Vdzym1xzm1pyr, nu);

        // edge back bot (y = 1) (z = -1)) DC
        edgederpquhats(uhats, ret[5], E->GetNcoeffs(), Vy1xzm1pyr, Vdxy1xzm1pyr,
                       Vdyy1xzm1pyr, Vdzy1xzm1pyr, nu);

        // edge left bot (z = -1), (x = -1) AD
        edgederpquhats(uhats, ret[6], E->GetNcoeffs(), Vxm1yzm1pyr,
                       Vdxxm1yzm1pyr, Vdyxm1yzm1pyr, Vdzxm1yzm1pyr, nu);

        // edge right bot ( z = -1) (x = 1) BC
        edgederpquhats(uhats, ret[7], E->GetNcoeffs(), Vx1yzm1pyr, Vdxx1yzm1pyr,
                       Vdyx1yzm1pyr, Vdzx1yzm1pyr, nu);
    }
}

void FilterStructurePres::edgederpquhats(
    Array<OneD, NekDouble> &uhats, Array<OneD, NekDouble> &ret, int modes,
    Array<OneD, NekDouble> Vxy, Array<OneD, NekDouble> Vxyd0,
    Array<OneD, NekDouble> Vxyd1, Array<OneD, NekDouble> Vxyd2, NekDouble nu)
{
    boost::ignore_unused(nu);
    int uhatstot = uhats.size();
    int totpts   = E3seg->GetBasis(0)->GetZ().size();

    Array<OneD, NekDouble> temp(totpts), temp2(modes), temp3(uhatstot);
    Array<OneD, NekDouble> pqeval(totpts);
    NekDouble v1, v2;

    for (int i = 0; i < totpts; i++)
    {
        Vmath::Vmul(modes, &Vxy[0] + i, totpts, &Vxyd0[0] + i, totpts,
                    &temp2[0], 1);
        v1 = Vmath::Vsum(modes, temp2, 1);
        Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);
        v2 = Vmath::Vsum(uhatstot, temp3, 1);

        v1 = v2 * v1;

        // At this point,
        // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

        // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)

        Vmath::Vmul(uhatstot, &Vxy[0] + i, totpts, &Vxy[0] + i, totpts,
                    &temp2[0], 1);

        v2 = Vmath::Vsum(uhatstot, temp2, 1);

        Vmath::Vmul(uhats.size(), &Vxyd0[i], totpts, &uhats[0], 1, &temp2[0],
                    1);

        v1 = v2 * Vmath::Vsum(uhats.size(), temp2, 1) - v1;

        pqeval[i] = v1;

        if (dimension > 1)
        {
            Vmath::Vmul(modes, &Vxy[0] + i, totpts, &Vxyd1[0] + i, totpts,
                        &temp2[0], 1);
            v1 = Vmath::Vsum(modes, temp2, 1);
            Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);

            v2 = Vmath::Vsum(uhatstot, temp3, 1);
            v1 = v2 * v1;

            // At this point,
            // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

            // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
            Vmath::Vmul(uhatstot, &Vxy[0] + i, totpts, &Vxy[0] + i, totpts,
                        &temp2[0], 1);

            v2 = Vmath::Vsum(uhatstot, temp2, 1);
            Vmath::Vmul(uhats.size(), &Vxyd1[i], totpts, &uhats[0], 1,
                        &temp2[0], 1);

            v1 = v2 * Vmath::Vsum(uhats.size(), temp2, 1) - v1;
            pqeval[i] += v1;
        }
        if (dimension == 3)
        {
            Vmath::Vmul(modes, &Vxy[0] + i, totpts, &Vxyd2[0] + i, totpts,
                        &temp2[0], 1);
            v1 = Vmath::Vsum(modes, temp2, 1);
            Vmath::Vmul(uhatstot, &Vxy[i], totpts, &uhats[0], 1, &temp3[0], 1);

            v2 = Vmath::Vsum(uhatstot, temp3, 1);
            v1 = v2 * v1;

            // At this point,
            // v1 = (sum(Vxm1.*Vdxm1, 2)).*(Vxm1*c) for 1 point

            // calculate (sum(Vxm1.*Vxm1, 2)).*(Vd*c)
            Vmath::Vmul(uhatstot, &Vxy[0] + i, totpts, &Vxy[0] + i, totpts,
                        &temp2[0], 1);

            v2 = Vmath::Vsum(uhatstot, temp2, 1);
            Vmath::Vmul(uhats.size(), &Vxyd2[i], totpts, &uhats[0], 1,
                        &temp2[0], 1);

            v1 = v2 * Vmath::Vsum(uhats.size(), temp2, 1) - v1;
            pqeval[i] += v1;
        }
    }
    //      pqeval = nu*pqeval;
    // if nu = 1, lower bnd
    // if nu = -1, upper bnd
    Vmath::Smul(pqeval.size(), nu, pqeval, 1, pqeval, 1);
    E3seg->FwdTrans(pqeval, ret);
}

void FilterStructurePres::project_surfaces(
    Array<OneD, NekDouble> uhats, Array<OneD, Array<OneD, NekDouble>> &ret,
    StdExpansion *E)
{

    if (E->DetShapeType() == LibUtilities::eHexahedron) // hex
    {
        // bot surface  z = -1
        surfuhats(uhats, ret[0], Vxyzm1, Equad);

        // right surface x = 1
        surfuhats(uhats, ret[1], Vx1yz, Equad);

        // surface top z = 1
        surfuhats(uhats, ret[2], Vxyz1, Equad);

        // left surface x = -1
        surfuhats(uhats, ret[3], Vxm1yz, Equad);

        // surface front y = -1
        surfuhats(uhats, ret[4], Vxym1z, Equad);

        // surface back y = 1
        surfuhats(uhats, ret[5], Vxy1z, Equad);
    }
    else if (E->DetShapeType() == LibUtilities::eTetrahedron) // tets
    {
        // surface bot z = -1 restriction on GD: x+y = 0 (ABC)
        surfuhats(uhats, ret[0], Vxyzm1tet, Etri);

        // surface left x = -1 restriction on GD: y + z = 0 (DAC)
        surfuhats(uhats, ret[1], Vxm1ypz0tet, Etri);

        // surf front y = -1   restriction on GD: x + z = 0 (DAB)
        surfuhats(uhats, ret[2], Vxpz0ym1tet, Etri);

        // surf DCB restriction  on GD: (x + y + z = -1)
        surfuhats(uhats, ret[3], Vxpypzm1tet, Etri);
    }
    else if (E->DetShapeType() == LibUtilities::ePyramid)
    {
        // surface bot z = -1
        surfuhats(uhats, ret[0], Vxyzm1pyr, Equad);

        // surface hypt x+z = 0
        surfuhats(uhats, ret[1], Vxmzypyr, Etri);

        // surface left x = -1
        surfuhats(uhats, ret[2], Vxm1ymzpyr, Etri);

        // surface front y = -1
        surfuhats(uhats, ret[3], Vxmzym1pyr, Etri);

        // surface back y +z =  0
        surfuhats(uhats, ret[4], Vxymzpyr, Etri);
    }
    else if (E->DetShapeType() == LibUtilities::ePrism)
    {
        // surf bot z = -1
        surfuhats(uhats, ret[0], Vxyzm1pri, Equad);

        // surface front x+z <= 0, y = -1
        surfuhats(uhats, ret[1], Vxym1zpri, Etri);

        // surf right (hypt) x+z <=0

        surfuhats(uhats, ret[3], Vxemzym1pri, Equad);

        // surface back y +z <=  0, y = 1
        surfuhats(uhats, ret[4], Vxy1zpri, Etri);

        // surface left x = -1
        surfuhats(uhats, ret[2], Vxm1yzpri, Equad);
    }
}

void FilterStructurePres::surfuhats(Array<OneD, NekDouble> &uhats,
                                    Array<OneD, NekDouble> &ret,
                                    Array<OneD, NekDouble> Vxyz,
                                    StdExpansion *Etemp)
{
    int modes = uhats.size();

    Array<OneD, NekDouble> temp(uhats.size());
    int totpts = Etemp->GetTotPoints();
    Array<OneD, NekDouble> vals(totpts);
    // Vxyz*uhats -> project to -> E3tri or Equad
    for (int k = 0; k < totpts; k++)
    {

        Vmath::Vmul(uhats.size(), &Vxyz[k], totpts, &uhats[0], 1, &temp[0], 1);

        vals[k] = Vmath::Vsum(modes, temp, 1);
    }

    Etemp->FwdTrans(vals, ret);
}

Array<OneD, NekDouble> FilterStructurePres::call_find_roots(
    Array<OneD, NekDouble> &uhats, NekDouble &avgiterGD,
    Array<OneD, Array<OneD, NekDouble>> &uhatsedges,
    Array<OneD, Array<OneD, NekDouble>> &surfaceuhats, NekDouble &minv,
    StdExpansion *E, NekDouble tol, NekDouble &roots1dtime,
    NekDouble &roots2dtime, NekDouble &roots3dtime, NekDouble &callstobasis,
    NekDouble nu)
{
    boost::ignore_unused(nu);
    Timer t1;
    int tot_pts_interp = 0;

    int dimension = 0;
    Array<OneD, NekDouble> temp(uhats.size());
    roots1dtime = 0.0;
    roots2dtime = 0.0;
    roots3dtime = 0.0;

    NekDouble inf = numeric_limits<double>::infinity();
    minv          = inf;

    NekDouble avgiterGDhold = 0, callstobasishold = 0;

    Array<OneD, Array<OneD, NekDouble>> rethold, storage;
    vector<vector<NekDouble>> retall;
    // EDGES
    if (E->DetShapeType() == LibUtilities::eQuadrilateral)
    {
        dimension = 2;
        retall    = vector<vector<NekDouble>>(dimension);
        storage = demo.storage2dq;
        rethold = NullNekDoubleArrayOfArray;

        t1.Start();
        demo.steepestgradient_descent2D(uhats, Equad, rethold, avgiterGDhold,
                                        callstobasishold, 1e-7, nu);
        t1.Stop();
        callstobasis += callstobasishold;
        roots2dtime += t1.TimePerTest(1);
        avgiterGD += avgiterGDhold;

        if (rethold != NullNekDoubleArrayOfArray)
        {
            retall[0].push_back(rethold[0][0]);
            retall[1].push_back(rethold[1][0]);
        }
    }
    else if (E->DetShapeType() == LibUtilities::eTriangle)
    {
        dimension = 2;

        retall  = vector<vector<NekDouble>>(dimension);
        storage = demo.storage2dt;
        rethold = NullNekDoubleArrayOfArray;
        t1.Start();

        demo.steepestgradient_descent2D(uhats, Etri, rethold, avgiterGDhold,
                                        callstobasishold, 1e-7, nu);
        t1.Stop();
        roots2dtime += t1.TimePerTest(1);
        callstobasis += callstobasishold;
        avgiterGD += avgiterGDhold;

        if (rethold != NullNekDoubleArrayOfArray)
        {
            retall[0].push_back(rethold[0][0]);
            retall[1].push_back(rethold[1][0]);
        }
    }
    else if (E->DetShapeType() == LibUtilities::eTetrahedron) // tet
    {
        dimension = 3;

        retall = vector<vector<NekDouble>>(dimension);
        //	  int numedges = 6;
        storage = demo.storage3dtet;
        retall[0].push_back(-1.0);
        retall[1].push_back(-1.0);
        retall[2].push_back(-1.0);
    }
    else if (E->DetShapeType() == LibUtilities::eHexahedron) // hex
    {
        dimension = 3;
        retall    = vector<vector<NekDouble>>(dimension);
        storage   = demo.storage3dhex;
    }
    else if (E->DetShapeType() == LibUtilities::ePrism) // prism
    {
        dimension = 3;

        retall       = vector<vector<NekDouble>>(dimension);
        storage      = demo.storage3dpri;
        int numedges = 9;
        for (int i = 0; i < numedges; i++)
        {
            t1.Start();
            rethold = demo.call_companion_rf(uhatsedges[i], C);
            t1.Stop();
            roots1dtime += t1.TimePerTest(1);
            if (rethold != NullNekDoubleArrayOfArray)
            {
                for (int p = 0; p < rethold[0].size(); p++)
                {
                    switch (i)
                    {
                        case 0: // edge front left (x = -1) (y = -1) EA
                            retall[0].push_back(-1);
                            retall[1].push_back(-1);
                            retall[2].push_back(rethold[0][p]);
                            break;
                        case 1: // edge front hypt (y = -1) (z = -x) EB
                            retall[0].push_back(rethold[0][p]);
                            retall[1].push_back(-1);
                            retall[2].push_back(-rethold[0][p]);
                            break;
                        case 2: // edge front bot (y = -1) (z = -1) AB
                            retall[0].push_back(rethold[0][p]);
                            retall[1].push_back(-1);
                            retall[2].push_back(-1);
                            break;
                        case 3: // edge back left (y = 1) (x = -1))
                            retall[0].push_back(-1);
                            retall[1].push_back(1);
                            retall[2].push_back(rethold[0][p]);
                            break;
                        case 4: // edge back hypt (y = 1) (z = -x)
                            retall[0].push_back(rethold[0][p]);
                            retall[1].push_back(1);
                            retall[2].push_back(-rethold[0][p]);
                            retall[0].push_back(-rethold[0][p]);
                            retall[1].push_back(1);
                            retall[2].push_back(rethold[0][p]);
                            break;
                        case 5: // edge back bot (y = 1) (z = -1)) DC
                            retall[0].push_back(rethold[0][p]);
                            retall[1].push_back(1);
                            retall[2].push_back(-1);
                            break;
                        case 6: // edge left top (x = -1) (z = 1)
                            retall[0].push_back(-1);
                            retall[1].push_back(rethold[0][p]);
                            retall[2].push_back(1);
                            break;
                        case 7: // edge bot  right (x = 1) (z = -1)
                            retall[0].push_back(1);
                            retall[1].push_back(rethold[0][p]);
                            retall[2].push_back(-1);
                            break;
                        case 8: // edge bot left (x = -1) (z = -1)
                            retall[0].push_back(-1);
                            retall[2].push_back(-1);
                            retall[1].push_back(rethold[0][p]);
                            break;
                        default:
                            cout << "\n wrong edge number in prism!\n";
                            exit(0);
                    }
                }
            }
        }
    }
    else if (E->DetShapeType() == LibUtilities::ePyramid) // pyr
    {
        dimension = 3;

        retall       = vector<vector<NekDouble>>(dimension);
        storage      = demo.storage3dpyr;
        int numedges = 8;
        for (int i = 0; i < numedges; i++)
        {
            t1.Start();
            rethold = demo.call_companion_rf(uhatsedges[i], C);
            t1.Stop();

            roots1dtime += t1.TimePerTest(1);
            if (rethold != NullNekDoubleArrayOfArray)
            {
                for (int p = 0; p < rethold[0].size(); p++)
                {
                    if (i == 0)
                    {
                        retall[0].push_back(-1);
                        retall[1].push_back(-1);
                        retall[2].push_back(rethold[0][p]);
                    }
                    else if (i == 1) // edge front hypt (y = -1) (z = -x)
                    {
                        retall[0].push_back(rethold[0][p]);
                        retall[2].push_back(-rethold[0][p]);
                        retall[1].push_back(-1);

                        retall[0].push_back(-rethold[0][p]);
                        retall[2].push_back(rethold[0][p]);
                        retall[1].push_back(-1);
                    }
                    else if (i == 2) // edge back hypt ( y = -z) (z +x = 0)
                    {
                        retall[0].push_back(rethold[0][p]);
                        retall[1].push_back(rethold[0][p]);
                        retall[2].push_back(-rethold[0][p]);

                        retall[0].push_back(-rethold[0][p]);
                        retall[1].push_back(-rethold[0][p]);
                        retall[2].push_back(rethold[0][p]);
                    }
                    else if (i == 3) // edge back left (y = -z) (x = -1) ED
                    {
                        retall[0].push_back(-1);
                        retall[1].push_back(-rethold[0][p]);
                        retall[2].push_back(rethold[0][p]);

                        retall[0].push_back(-1);
                        retall[1].push_back(rethold[0][p]);
                        retall[2].push_back(-rethold[0][p]);
                    }
                    else if (i == 4) // edge front bot (y = -1) (z = -1)
                    {
                        retall[1].push_back(-1);
                        retall[0].push_back(rethold[0][p]);
                        retall[2].push_back(-1);
                    }
                    else if (i == 5) // edge back bot (y = 1) (z = -1)
                    {
                        retall[1].push_back(1);
                        retall[0].push_back(rethold[0][p]);
                        retall[2].push_back(-1);
                    }
                    else if (i == 6) // edge left bot (z = -1), (x = -1)
                    {
                        retall[0].push_back(-1);
                        retall[1].push_back(rethold[0][p]);
                        retall[2].push_back(-1);
                    }
                    else if (i == 7) // edge right bot ( z = -1) (x = 1)
                    {
                        retall[2].push_back(-1);
                        retall[1].push_back(rethold[0][p]);
                        retall[0].push_back(1);
                    }
                }
            }
        }
    }

    // SURFACES:
    if (dimension == 3)
    {
        if (E->DetShapeType() == LibUtilities::eTetrahedron) // tet
        {
            t1.Start();

            demo.steepestgradientdescent3D(
                uhats, Etet, demo.storage3dtet, demo.coordtet, demo.coordmidtet,
                demo.midptevaltet, rethold, avgiterGDhold, callstobasishold,
                1e-7, 1, tot_pts_interp);

            t1.Stop();
            roots3dtime += t1.TimePerTest(1);

            avgiterGD += avgiterGDhold;
            callstobasis += callstobasishold;
            if (rethold != NullNekDoubleArrayOfArray)
            {
                for (int k = 0; k < dimension; k++)
                {
                    if (rethold[k][0] < inf)
                    {
                        retall[k].push_back(rethold[k][0]);
                    }
                }
            }
        }
        else if (E->DetShapeType() == LibUtilities::eHexahedron) // hex
        {
            // add all vertices

            retall[0].push_back(1.0);
            retall[1].push_back(1.0);
            retall[2].push_back(1.0);

            retall[0].push_back(-1.0);
            retall[1].push_back(1.0);
            retall[2].push_back(1.0);

            retall[0].push_back(1.0);
            retall[1].push_back(-1.0);
            retall[2].push_back(1.0);

            retall[0].push_back(1.0);
            retall[1].push_back(1.0);
            retall[2].push_back(-1.0);

            retall[0].push_back(-1.0);
            retall[1].push_back(-1.0);
            retall[2].push_back(1.0);

            retall[0].push_back(1.0);
            retall[1].push_back(-1.0);
            retall[2].push_back(-1.0);

            retall[0].push_back(-1.0);
            retall[1].push_back(1.0);
            retall[2].push_back(-1.0);

            retall[0].push_back(-1.0);
            retall[1].push_back(-1.0);
            retall[2].push_back(-1.0);

            // 3d gd:

            t1.Start();
            demo.steepestgradientdescent3D(
                uhats, Ehex, demo.storage3dhex, demo.coordhex, demo.coordmidhex,
                demo.midptevalhex, rethold, avgiterGDhold, callstobasishold,
                1e-7, 1, tot_pts_interp);
            callstobasis += callstobasishold;
            tot_interp_points += tot_pts_interp;
            t1.Stop();
            roots3dtime += t1.TimePerTest(1);
            avgiterGD += avgiterGDhold;
            if (rethold != NullNekDoubleArrayOfArray)
            {
                for (int k = 0; k < dimension; k++)
                {
                    if (rethold[k][0] < inf)
                    {
                        retall[k].push_back(rethold[k][0]);
                    }
                }
            }
        }
        else if (E->DetShapeType() == LibUtilities::ePrism) // pri
        {
            NekDouble avgiterGDhold;

            // surface bot z = -1
            t1.Start();
            demo.steepestgradient_descent2D(surfaceuhats[0], Equad, rethold,
                                            avgiterGDhold, callstobasishold,
                                            tol, 1);
            t1.Stop();
            if (rethold != NullNekDoubleArrayOfArray)
            {
                roots2dtime += t1.TimePerTest(1);
                avgiterGD += avgiterGDhold;
                retall[2].push_back(-1);
                retall[1].push_back(rethold[1][0]);
                retall[0].push_back(rethold[0][0]);
            }

            // surface front x+z <= 0, y = -1
            t1.Start();
            demo.steepestgradient_descent2D(surfaceuhats[1], Etri, rethold,
                                            avgiterGDhold, callstobasishold,
                                            tol, 1);
            t1.Stop();
            if (rethold != NullNekDoubleArrayOfArray)
            {
                roots2dtime += t1.TimePerTest(1);
                avgiterGD += avgiterGDhold;
                retall[1].push_back(-1);
                retall[0].push_back(rethold[0][0]);
                retall[2].push_back(rethold[1][0]);
            }

            // surf right (hypt) x+z <=0
            t1.Start();
            demo.steepestgradient_descent2D(surfaceuhats[3], Equad, rethold,
                                            avgiterGDhold, callstobasishold,
                                            tol, 1);
            t1.Stop();
            if (rethold != NullNekDoubleArrayOfArray)
            {
                roots2dtime += t1.TimePerTest(1);
                avgiterGD += avgiterGDhold;
                retall[0].push_back(rethold[0][0]);
                retall[1].push_back(rethold[1][0]);
                retall[2].push_back(-rethold[0][0]);
            }

            // surface back y +z <=  0, y = 1
            t1.Start();
            demo.steepestgradient_descent2D(surfaceuhats[4], Etri, rethold,
                                            avgiterGDhold, callstobasishold,
                                            tol, 1);
            t1.Stop();
            if (rethold != NullNekDoubleArrayOfArray)
            {
                roots2dtime += t1.TimePerTest(1);
                avgiterGD += avgiterGDhold;
                retall[0].push_back(rethold[0][0]);
                retall[1].push_back(1);
                retall[2].push_back(rethold[1][0]);
            }

            // surface left x = -1
            t1.Start();
            demo.steepestgradient_descent2D(surfaceuhats[2], Equad, rethold,
                                            avgiterGDhold, callstobasishold,
                                            tol, 1);
            t1.Stop();
            if (rethold != NullNekDoubleArrayOfArray)
            {
                roots2dtime += t1.TimePerTest(1);
                avgiterGD += avgiterGDhold;
                retall[0].push_back(-1);
                retall[1].push_back(rethold[0][0]);
                retall[2].push_back(rethold[1][0]);
            }
            // 3d gd:
            t1.Start();
            demo.steepestgradientdescent3D(
                uhats, Epri, demo.storage3dpri, demo.coordpri, demo.coordmidpri,
                demo.midptevalpri, rethold, avgiterGDhold, callstobasis, 1e-7,
                1, tot_pts_interp);
            t1.Stop();
            roots3dtime += t1.TimePerTest(1);

            avgiterGD += avgiterGDhold;
            if (rethold != NullNekDoubleArrayOfArray)
            {
                for (int k = 0; k < dimension; k++)
                {
                    if (rethold[k][0] < inf)
                    {
                        retall[k].push_back(rethold[k][0]);
                    }
                }
            }
        }
        else if (E->DetShapeType() == LibUtilities::ePyramid) // pyr
        {
            int numsurfaces = 5;
            NekDouble avgiterGDhold;
            for (int i = 0; i < numsurfaces; i++)
            {
                if (i == 0) // bot surface: quad, z = -1
                {
                    t1.Start();
                    demo.steepestgradient_descent2D(surfaceuhats[i], Equad,
                                                    rethold, avgiterGDhold,
                                                    callstobasishold, tol, 1);
                    t1.Stop();
                    callstobasis += callstobasishold;
                    if (rethold != NullNekDoubleArrayOfArray)
                    {
                        roots2dtime += t1.TimePerTest(1);
                        avgiterGD += avgiterGDhold;

                        retall[2].push_back(-1);
                        retall[1].push_back(rethold[1][0]);
                        retall[0].push_back(rethold[0][0]);
                    }
                }
                else // rest surfaces: tri
                {
                    t1.Start();
                    demo.steepestgradient_descent2D(surfaceuhats[i], Etri,
                                                    rethold, avgiterGDhold,
                                                    callstobasishold, 1e-6, 1);

                    t1.Stop();
                    callstobasis += callstobasishold;
                    roots2dtime += t1.TimePerTest(1);
                    avgiterGD += avgiterGDhold;
                    if (rethold != NullNekDoubleArrayOfArray)
                    {
                        switch (i)
                        {
                            case 1: // x+z = 0
                                retall[0].push_back(rethold[0][0]);
                                retall[1].push_back(rethold[1][0]);
                                retall[2].push_back(-rethold[0][0]);
                                break;
                            case 2: // x = -1
                                retall[2].push_back(rethold[1][0]);
                                retall[1].push_back(rethold[0][0]);
                                retall[0].push_back(-1.0);
                                break;
                            case 3: // y = -1
                                retall[2].push_back(rethold[1][0]);
                                retall[1].push_back(-1.0);
                                retall[0].push_back(rethold[0][0]);
                                break;
                            case 4: // y+z = 0
                                retall[2].push_back(-rethold[1][0]);
                                retall[1].push_back(rethold[1][0]);
                                retall[0].push_back(rethold[0][0]);
                                break;
                        }
                    }
                }
            }

            // 3d gd:
            t1.Start();
            demo.steepestgradientdescent3D(
                uhats, Epyr, demo.storage3dpyr, demo.coordpyr, demo.coordmidpyr,
                demo.midptevalpyr, rethold, avgiterGDhold, callstobasis, 1e-7,
                1, tot_pts_interp);
            t1.Stop();
            roots3dtime += t1.TimePerTest(1);

            avgiterGD += avgiterGDhold;
            if (rethold.size() > 0)
            {
                for (int k = 0; k < dimension; k++)
                {
                    if (rethold[k][0] < inf)
                    {
                        retall[k].push_back(rethold[k][0]);
                    }
                }
            }
        }
    }
    Array<OneD, NekDouble> evalroots, ret(dimension);
    if (retall[0].size() > 0)
    {
        Array<OneD, Array<OneD, NekDouble>> tmpcoord(dimension);
        for (int p = 0; p < dimension; p++)
        {
            tmpcoord[p] = Array<OneD, NekDouble>(retall[0].size());
        }

        for (int p = 0; p < dimension; p++)
        {
            for (int q = 0; q < retall[0].size(); q++)
            {
                tmpcoord[p][q] = retall[p][q];
            }
        }

        NekDouble tempmin = inf;
        t1.Start();
        evalroots =
            E->PhysEvaluateBasis(tmpcoord, storage, NullNekDouble1DArray,
                                 NullNekDouble1DArray, NullNekDouble1DArray);

        tot_interp_points += tmpcoord[0].size();
        t1.Stop();
        roots3dtime += t1.TimePerTest(1);
        callstobasis += t1.TimePerTest(1);
        Array<OneD, NekDouble> minvandx(dimension + 1);
        Array<OneD, NekDouble> tmparr(retall[0].size());
        t1.Start();
        demo.pq(uhats, tmpcoord, evalroots, minvandx, tmparr, nu);
        t1.Stop();
        tempmin = minvandx[0];
        for (int k = 0; k < dimension; k++)
        {
            ret[k] = minvandx[k + 1];
        }
        minv = tempmin;
    }
    return ret;
}

void FilterStructurePres::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(pFields, time);
}

bool FilterStructurePres::v_IsTimeDependent()
{
    return false;
}
} // namespace SolverUtils
} // namespace Nektar
