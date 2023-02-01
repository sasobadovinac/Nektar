////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessBodyFittedVelocity.cpp
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
//  Description: Generate the body-fitted coordinate system and compute the
//               body-fitted velocity components.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include <FieldUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

#include "ProcessMapping.h"
#include <GlobalMapping/Mapping.h>

#include "ProcessBodyFittedVelocity.h"
#include "ProcessWallNormalData.h"

using namespace std;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessBodyFittedVelocity::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "bodyFittedVelocity"),
        ProcessBodyFittedVelocity::create,
        "Convert the velocity components from the Cartesian coordinate into "
        "the body-fitted coordinate.");

ProcessBodyFittedVelocity::ProcessBodyFittedVelocity(FieldSharedPtr f)
    : ProcessBoundaryExtract(f)
{
    f->m_writeBndFld =
        false; // turned on in the father class ProcessBoundaryExtract

    // bnd is read from father class ProcessBoundaryExtract
    m_config["assistDir"] =
        ConfigOption(false, "0.0,0.0,1.0",
                     "The normal direction of the assistant plane, where \
                             the first body-fitted coordinate direction is located. \
                             The assistance plane is defined in point normal form. \
                             Default = [0,0,1]");
    m_config["checkAngle"] =
        ConfigOption(true, "0", "The flag for checking the distance vector \
                             perpendicular to the nearest element or not. It \
                             does not need to be turned on if geometry is \
                             smooth and there is no special consideration of \
                             the body-fitted coordinate system. It would be \
                             helpful to set up a desired coordinate system \
                             based on the surface with with sharp corners, \
                             such as at the presence of a step or gap. \
                             Default = false.");
    m_config["distTol"] = ConfigOption(
        false, "1.0e-12", "The distance tolerence to properly find out the \
                             nearest boundary element. It works together with \
                             the CHECKANGLE option to set up the body-fitted \
                             coordinate system near the corner. \
                             Default = 1.0e-12.");
    m_config["bfcOut"] = ConfigOption(
        true, "1", "The flag control whether the body-fitted coordinate \
                             system will be output. Default = true.");
}

ProcessBodyFittedVelocity::~ProcessBodyFittedVelocity()
{
}

/*Note
 * This module converts the velocity components into the body-fitted coordinate.
 * It works for both incompressible and compressible cases.
 * The only requirement is the [u,v] fields the 2D case or [u,v,w] fields for
 * 2.5D/3D case exist in the input fld fild.
 * The input cases can be 2D, 2.5D and 3D.
 * The data will be exported with .fld extension.
 *
 * Mapping is not supported at the moment
 *
 * The user defined parameters are: bnd, fixedDir
 * - "bnd" for the wall reference id, which is defiend in the father class
 *   "ProcessBoundaryExtract".
 * - "fixedDir" is the dorection at which the velocity component is not rotated.
 */
void ProcessBodyFittedVelocity::Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::Process(vm);

    // Get dim to store data
    const int nFields    = m_f->m_variables.size();
    const int nCoordDim  = m_f->m_exp[0]->GetCoordim(0);
    const int nSpaceDim  = nCoordDim + m_f->m_numHomogeneousDir;
    const int nAddFields = nSpaceDim + 1; // distanceToWall, u_bfc, v_bfc, w_bfc

    // Tolerence setting
    // To include more elmts for further check
    const NekDouble distTol = m_config["distTol"].as<NekDouble>();
    const NekDouble geoTol  = 1.0e-12; // To check if two pnts are too close
    const NekDouble dirTol  = 1.0e-4;  // To check if two unit vecs are parallel
    const NekDouble iterTol = 1.0e-12; // To check if locCoord iter convergence

    const bool isCheckAngle = m_config["checkAngle"].as<bool>();
    const bool isOutputBfc  = m_config["bfcOut"].as<bool>();

    // Get the assist vector
    vector<NekDouble> assistDir;
    ASSERTL0(ParseUtils::GenerateVector(m_config["assistDir"].as<string>(),
                                        assistDir),
             "Failed to interpret assistance direction");

    Array<OneD, NekDouble> assistVec(3); // Save the assist vector
    for (int i = 0; i < 3; ++i)
    {
        assistVec[i] = (assistDir.size() > i) ? assistDir[i] : 0.0;
    }
    if (nCoordDim == 2)
    {
        assistVec[0] = 0.0;
        assistVec[1] = 0.0;
    }

    NekDouble norm = sqrt(Vmath::Dot(3, assistVec, 1, assistVec, 1));
    if (norm < geoTol)
    {
        ASSERTL0(false, "Error. The amplitude of assist vector is smaller than \
                         the geometry tolerence 1.0e-12.");
    }
    Vmath::Smul(3, 1.0 / norm, assistVec, 1, assistVec, 1);

    // Get the bnd id
    // We only use the first one [0], check ProcessBoundaryExtract for more info
    SpatialDomains::BoundaryConditions bcs(m_f->m_session,
                                           m_f->m_exp[0]->GetGraph());
    const SpatialDomains::BoundaryRegionCollection bregions =
        bcs.GetBoundaryRegions();
    map<int, int> BndRegionMap;
    int cnt = 0;
    for (auto &breg_it : bregions)
    {
        BndRegionMap[breg_it.first] = cnt++;
    }
    int bnd = BndRegionMap[m_f->m_bndRegionsToWrite[0]];

    // Generate the distance array and three arrays for directional vectors
    // bfcsDir[i][j][k]
    //   i - dir vec:   0-main tangential, 1-normal, (2-minor tangential)
    //   j - component: 0-x, 1-y, (2-z)
    //   k - point id
    int npoints = m_f->m_exp[0]->GetNpoints();
    Array<OneD, NekDouble> distance(npoints, 9999);
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> bfcsDir(nSpaceDim);
    for (int i = 0; i < nSpaceDim; ++i)
    {
        bfcsDir[i] = Array<OneD, Array<OneD, NekDouble>>(nSpaceDim);

        for (int j = 0; j < nSpaceDim; ++j)
        {
            bfcsDir[i][j] = Array<OneD, NekDouble>(npoints);
        }
    }

    // Key function. At each quadrature point inside the domian, compute the
    // body-fitted coordinate system with respect to the input bnd id
    GenPntwiseBodyFittedCoordSys(bnd, assistVec, distance, bfcsDir,
                                 isCheckAngle, distTol, iterTol, dirTol,
                                 geoTol);

    // Compute the velocity
    Array<OneD, Array<OneD, NekDouble>> vel_car(nSpaceDim); // Cartesian
    Array<OneD, Array<OneD, NekDouble>> vel_bfc(nSpaceDim); // body-fiitted
    Array<OneD, NekDouble> vel_tmp(npoints);

    for (int i = 0; i < nSpaceDim; ++i)
    {
        vel_bfc[i] = Array<OneD, NekDouble>(npoints, 0.0);
    }

    // Get the velocity in Cartesian coordinate system
    GetVelAndConvertToCartSys(vel_car);

    // Project the velocity into the body-fitted coordinate system
    for (int i = 0; i < nSpaceDim; ++i) // loop for bfc velocity
    {
        for (int j = 0; j < nSpaceDim; ++j)
        {
            Vmath::Vmul(npoints, vel_car[j], 1, bfcsDir[i][j], 1, vel_tmp, 1);
            Vmath::Vadd(npoints, vel_tmp, 1, vel_bfc[i], 1, vel_bfc[i], 1);
        }
    }

    // Add var names
    m_f->m_variables.push_back("distanceToWall");
    m_f->m_variables.push_back("u_bfc");
    m_f->m_variables.push_back("v_bfc");
    if (nSpaceDim == 3)
    {
        m_f->m_variables.push_back("w_bfc");
    }
    else if (nSpaceDim == 1)
    {
        ASSERTL0(false, "Velocity in 1D case is already in the body fitted \
                         coordinate. The dimension should be 2 or 3.");
    }

    // Add the new fields (distance + u,v,w_bfc)
    // copy distance
    m_f->m_exp.resize(nFields + nAddFields);

    m_f->m_exp[nFields] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
    Vmath::Vcopy(npoints, distance, 1, m_f->m_exp[nFields]->UpdatePhys(), 1);
    m_f->m_exp[nFields]->FwdTransLocalElmt(distance,
                                           m_f->m_exp[nFields]->UpdateCoeffs());

    // copy vel_bfc
    for (int i = 1; i < nAddFields; ++i)
    {
        m_f->m_exp[nFields + i] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
        Vmath::Vcopy(npoints, vel_bfc[i - 1], 1,
                     m_f->m_exp[nFields + i]->UpdatePhys(), 1);
        m_f->m_exp[nFields + i]->FwdTransLocalElmt(
            vel_bfc[i - 1], m_f->m_exp[nFields + i]->UpdateCoeffs());
    }

    // Output the the body-fitted coordinate into the file
    // The filter will read the coordiante through file
    // bfcsDir[i][j][k]
    //   i - dir vec:   0-main tangential, 1-normal, (2-minor tangential)
    //   j - component: 0-x, 1-y, (2-z)
    //   k - point id
    if (isOutputBfc)
    {
        m_f->m_variables.push_back("bfc_dir0_x");
        m_f->m_variables.push_back("bfc_dir0_y");
        if (nSpaceDim == 3)
        {
            m_f->m_variables.push_back("bfc_dir0_z");
        }
        m_f->m_variables.push_back("bfc_dir1_x");
        m_f->m_variables.push_back("bfc_dir1_y");
        if (nSpaceDim == 3)
        {
            m_f->m_variables.push_back("bfc_dir1_z");

            m_f->m_variables.push_back("bfc_dir2_x");
            m_f->m_variables.push_back("bfc_dir2_y");
            m_f->m_variables.push_back("bfc_dir2_z");
        }

        int nAddFields_coord = (nSpaceDim == 3) ? 9 : 4;
        m_f->m_exp.resize(nFields + nAddFields + nAddFields_coord);

        cnt = 0;
        for (int i = 0; i < nSpaceDim; ++i)
        {
            for (int j = 0; j < nSpaceDim; ++j)
            {

                m_f->m_exp[nFields + nAddFields + cnt] =
                    m_f->AppendExpList(m_f->m_numHomogeneousDir);

                Vmath::Vcopy(
                    npoints, bfcsDir[i][j], 1,
                    m_f->m_exp[nFields + nAddFields + cnt]->UpdatePhys(), 1);
                m_f->m_exp[nFields + nAddFields + cnt]->FwdTransLocalElmt(
                    bfcsDir[i][j],
                    m_f->m_exp[nFields + nAddFields + cnt]->UpdateCoeffs());

                ++cnt;
            }
        }
    }
}

/**
 * @brief Compute the point-to-point distance to find the estimated cloest
 * element.
 * @param pts    Global coordinate of the quadrature points in the inner
 * element.
 * @param pId    Id of the inner point of interest in the current loop.
 * @param bndPts Global coordinate of the quadrature points in the boundary
 * element.
 */
NekDouble ProcessBodyFittedVelocity::PntToBndElmtPntDistance(
    const Array<OneD, Array<OneD, NekDouble>> &pts, const int pId,
    const Array<OneD, Array<OneD, NekDouble>> &bndPts)
{
    const int nCoordDim = pts.size();

    NekDouble dist = 9999, dist_tmp;

    for (int i = 0; i < bndPts[0].size(); ++i)
    {
        // Compute distance
        dist_tmp = 0.0;
        for (int j = 0; j < nCoordDim; ++j)
        {
            dist_tmp +=
                (bndPts[j][i] - pts[j][pId]) * (bndPts[j][i] - pts[j][pId]);
        }
        dist_tmp = sqrt(dist_tmp);

        if (dist_tmp < dist)
        {
            dist = dist_tmp;
        }
    }

    return dist;
}

/**
 * @brief Compute the local coordinate for the nearest point on the given
 *        2D boundary element to the input point.
 * @param inGloCoord  Global coordinate for the input point.
 * @param bndGeom     Geometry of the boundary element to search from.
 * @param pts         Global coordinate of the quadrature points in the boundary
 * element.
 * @param locCoord    Local coordinate of the result.
 * @param gloCoord    Global coordinate of the result.
 * @param dist        Distance from the input point to the boundary element.
 * @param iterTol     Iteration tolerence.
 * @param iterMax     Max iteration steps.
 */
bool ProcessBodyFittedVelocity::LocCoordForNearestPntOnBndElmt_2D(
    const Array<OneD, const NekDouble> &inGloCoord,
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array<OneD, Array<OneD, NekDouble>> &pts,
    Array<OneD, NekDouble> &locCoord, Array<OneD, NekDouble> &gloCoord,
    NekDouble &dist, const NekDouble iterTol, const int iterMax)
{

    // Initial settings
    Array<OneD, NekDouble> etaLR(2);      // range [-1,1]
    etaLR[0] = -1.0;                      // left
    etaLR[1] = 1.0;                       // right
    NekDouble tmpLx, tmpLy, tmpRx, tmpRy; // tmp values for L/R
    NekDouble distL2, distR2;             // distance square

    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    locCoord[0]                               = -2.0;
    int cnt                                   = 0;
    bool isConverge                           = false;
    while (cnt < iterMax)
    {
        tmpLx = bndXmap->PhysEvaluate(etaLR, pts[0]);
        tmpLy = bndXmap->PhysEvaluate(etaLR, pts[1]);
        tmpRx = bndXmap->PhysEvaluate(etaLR + 1, pts[0]);
        tmpRy = bndXmap->PhysEvaluate(etaLR + 1, pts[1]);

        distL2 = (tmpLx - inGloCoord[0]) * (tmpLx - inGloCoord[0]) +
                 (tmpLy - inGloCoord[1]) * (tmpLy - inGloCoord[1]);
        distR2 = (tmpRx - inGloCoord[0]) * (tmpRx - inGloCoord[0]) +
                 (tmpRy - inGloCoord[1]) * (tmpRy - inGloCoord[1]);

        if (distL2 >= distR2)
        {
            etaLR[0] = 0.5 * (etaLR[0] + etaLR[1]);
        }
        else
        {
            etaLR[1] = 0.5 * (etaLR[0] + etaLR[1]);
        }

        if ((etaLR[1] - etaLR[0]) < iterTol)
        {
            locCoord[0] = 0.5 * (etaLR[0] + etaLR[1]);
            gloCoord[0] = 0.5 * (tmpLx + tmpRx);
            gloCoord[1] = 0.5 * (tmpLy + tmpRy);
            dist        = sqrt(0.5 * (distL2 + distR2));
            isConverge  = true;
            break;
        }

        ++cnt;
    }

    // Warning if failed
    if (cnt >= iterMax)
    {
        WARNINGL1(false, "Bisection iteration is not converged");
    }

    return isConverge;
}

/**
 * @brief Compute the local coordinate for the nearest point on the given
 *        boundary element to the input point. The locCoord is the position to
 *        set up the body-fitted coordinate. This function works as a driver.
 * @param inGloCoord  Global coordinate for the input point.
 * @param bndGeom     Geometry of the boundary element to search from.
 * @param locCoord    Local coordinate of the result.
 * @param gloCoord    Global coordinate of the result.
 * @param dist        Distance from the input point to the boundary element.
 * @param iterTol     Iteration tolerence.
 * @param iterMax     Max iteration steps.
 */
bool ProcessBodyFittedVelocity::LocCoordForNearestPntOnBndElmt(
    const Array<OneD, const NekDouble> &inGloCoord,
    SpatialDomains::GeometrySharedPtr bndGeom, Array<OneD, NekDouble> &locCoord,
    Array<OneD, NekDouble> &gloCoord, NekDouble &dist, const NekDouble iterTol,
    const int iterMax)
{
    bool isConverge     = false;
    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0); // =2 for 2.5D cases

    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    int nBndPts                               = bndXmap->GetTotPoints();

    Array<OneD, Array<OneD, const NekDouble>> bndCoeffs(nCoordDim);
    Array<OneD, Array<OneD, NekDouble>> bndPts(nCoordDim);
    for (int i = 0; i < nCoordDim; ++i)
    {
        bndPts[i]    = Array<OneD, NekDouble>(nBndPts);
        bndCoeffs[i] = bndGeom->GetCoeffs(i);
        bndXmap->BwdTrans(bndCoeffs[i], bndPts[i]);
    }

    // Compute the locCoord for the pnt on the bnd elmt
    if (nCoordDim == 2)
    {
        // Bisection search
        isConverge = LocCoordForNearestPntOnBndElmt_2D(
            inGloCoord, bndGeom, bndPts, locCoord, gloCoord, dist, iterTol,
            iterMax);
    }
    else
    {
        // Bi-bisection search
        ASSERTL0(false, "Not available at the moment.");
    }

    return isConverge;
}

/**
 * @brief Compute the normalized cross product for two 2D vectors.
 *        vec3 = vec1 x vec2
 */
void ProcessBodyFittedVelocity::ScaledCrosssProduct(
    const Array<OneD, NekDouble> &vec1, const Array<OneD, NekDouble> &vec2,
    Array<OneD, NekDouble> &vec3)
{
    vec3[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    vec3[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    vec3[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

    NekDouble coef;
    coef =
        1.0 / sqrt(vec3[0] * vec3[0] + vec3[1] * vec3[1] + vec3[2] * vec3[2]);

    vec3[0] = vec3[0] * coef;
    vec3[1] = vec3[1] * coef;
    vec3[2] = vec3[2] * coef;
}

/**
 * @brief At each quadrature point inside the domian, compute the body-fitted
 *        coordinate system with respect to the input boundary id.
 * @param targetBndId  Target boundary id.
 * @param assistVec    Unit assistant vector, the cross product of inward-
 *                     pointing wall normalId and wihch gives of the main
 *                     tangential direction of the body-fitted system.
 * @param bfcsDir      Pointwise body-fitted coordinate system.
 * @param isPerpendicularCondition Flag for using perpendicular check or not
 * @param distTol      Distance tolerence. Used to find the boundary elements
 *                     for local coordinate iteration.
 * @param iterTol      Iteration tolerence. Used to check iteration convergence.
 * @param dirTol       Direction tolerencce. Used to check if the inner product
 *                     of two unit vectors is cloes enough to 1.0.
 * @param geoTol       Geometry tolerence. Used as the relative tolerence for
 *                     local coord and distance absolute tolerence.
 */
void ProcessBodyFittedVelocity::GenPntwiseBodyFittedCoordSys(
    const int targetBndId, const Array<OneD, NekDouble> assistVec,
    Array<OneD, NekDouble> &distance,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &bfcsDir,
    const bool isCheckAngle, const NekDouble distTol, const NekDouble iterTol,
    const NekDouble dirTol, const NekDouble geoTol)
{

    const int nFields       = m_f->m_variables.size();
    const int nCoordDim     = m_f->m_exp[0]->GetCoordim(0);
    const int nSpaceDim     = nCoordDim + m_f->m_numHomogeneousDir;
    const int nBndLcoordDim = nCoordDim - 1;

    // Get expansion list for boundary
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp(nFields);
    for (int i = 0; i < nFields; ++i)
    {
        BndExp[i] = m_f->m_exp[i]->UpdateBndCondExpansion(targetBndId);
    }

    const int nElmts    = m_f->m_exp[0]->GetNumElmts(); // num of inner element
    const int nBndElmts = BndExp[0]->GetNumElmts();     // num of bnd element

    // Use the outward-pointing normal at quardrature pts on bnd elmt as dir ref
    Array<OneD, Array<OneD, NekDouble>> normalQ;
    m_f->m_exp[0]->GetBoundaryNormals(targetBndId, normalQ);

    // Vars in the loop
    SpatialDomains::GeometrySharedPtr bndGeom;
    StdRegions::StdExpansionSharedPtr bndXmap;
    Array<OneD, Array<OneD, NekDouble>> bndCoeffs(nCoordDim);
    Array<OneD, Array<OneD, NekDouble>> bndPts(nCoordDim);

    int nPts, nBndPts;
    NekDouble dist_tmp, dist_bak;
    bool isConverge;
    vector<int> bndElmtIds;
    int bndElmtId, bndElmtId_bak;
    int phys_offset, bnd_phys_offset;

    Array<OneD, NekDouble> inGloCoord(3, 0.0); // inner pnt
    Array<OneD, NekDouble> locCoord(nBndLcoordDim), locCoord_tmp(nBndLcoordDim);
    Array<OneD, NekDouble> gloCoord(3, 0.0), gloCoord_tmp(3, 0.0);
    Array<OneD, NekDouble> normal(3), normalRef(3);
    Array<OneD, NekDouble> tangential1(3),
        tangential2(3); // 1 - main, 2 - minor
    Array<OneD, NekDouble> normalChk(3);
    ProcessWallNormalData wnd(m_f); // wallNormalData object to use its routine

    // backup for angle check
    Array<OneD, NekDouble> locCoord_bak(nBndLcoordDim), gloCoord_bak(3, 0.0);

    //-------------------------------------------------------------------------
    // Loop all the quadrature points of the field, and find out the closest
    // point on the bnd for each of them. The bnd element contains this bnd
    // point is recorded.

    // Loop inner element
    for (int eId = 0; eId < nElmts; ++eId) // element id, nElmts
    {
        // Get inner points
        nPts = m_f->m_exp[0]->GetExp(eId)->GetTotPoints();
        Array<OneD, Array<OneD, NekDouble>> pts(nCoordDim);
        for (int i = 0; i < nCoordDim; ++i)
        {
            pts[i] = Array<OneD, NekDouble>(nPts, 0.0);
        }

        if (nCoordDim == 2)
        {
            m_f->m_exp[0]->GetExp(eId)->GetCoords(pts[0], pts[1]);
        }
        else if (nCoordDim == 3)
        {
            m_f->m_exp[0]->GetExp(eId)->GetCoords(pts[0], pts[1], pts[2]);
        }

        // Get the offset for the current elmt
        phys_offset = m_f->m_exp[0]->GetPhys_Offset(eId);

        // Loop points in the inner element
        for (int pId = 0; pId < nPts; ++pId)
        {

            // Compute the distance from the pnt (bnd elmt) to the inner pnt
            // Step 1 - Estimate search
            //  - Find the closest quadrature pnt among all the bnd elmt;
            //  - The bnd elmt which contains this pnt is saved for Step 2.
            //  * If the quadtature pnt is on the edge or corner of the bnd
            //    elmt, there will be more than one bnd elmt which contains
            //    this point.

            bndElmtIds.clear(); // clear the bnd elmt id vec for a new inner pnt

            // Loop bnd elmt to find the elmt to build body-fitted coord sys
            for (int beId = 0; beId < nBndElmts; ++beId)
            {
                // Get geomery and points
                bndGeom = BndExp[0]->GetExp(beId)->GetGeom();
                bndXmap = bndGeom->GetXmap();
                nBndPts = bndXmap->GetTotPoints();

                for (int i = 0; i < nCoordDim; ++i)
                {
                    bndPts[i]    = Array<OneD, NekDouble>(nBndPts);
                    bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
                    bndXmap->BwdTrans(bndCoeffs[i], bndPts[i]);
                }

                // Compute the distance (estimate)
                dist_tmp = PntToBndElmtPntDistance(pts, pId, bndPts);

                if ((dist_tmp < (distance[phys_offset + pId] - geoTol)) &&
                    (dist_tmp < (distance[phys_offset + pId] - distTol)))
                {
                    // Find a closer quadrature point from a different bnd elmt
                    // Update distance and bnd elmt id vector
                    distance[phys_offset + pId] = dist_tmp;
                    bndElmtIds.clear();
                    bndElmtIds.push_back(beId);
                }
                else if ((dist_tmp < (distance[phys_offset + pId] + geoTol)) ||
                         (dist_tmp < (distance[phys_offset + pId] + distTol)))
                {
                    // The same quadrature point from a different bnd elmt is
                    // found for the other time.
                    // Keep the distance and add the bnd elmt id to the vector
                    bndElmtIds.push_back(beId);
                }

            } // end of bnd elmt loop

            // Step 2 - Accurate search
            //  - Find the accurate distance and locCoord for the nearest pnt
            //    (to the inner point) on the bnd elmt
            //  - This process will be repeated on all the bnd elmts in the bnd
            //    elmt id vector.

            // Generate the coordinate for the inner point
            for (int i = 0; i < nCoordDim; ++i)
            {
                inGloCoord[i] = pts[i][pId];
            }

            // reset the dist to make sure the condition for if will be
            // satisfied
            distance[phys_offset + pId] = 9999;
            bndElmtId                   = -1;
            dist_bak                    = 9999;
            bndElmtId_bak               = -1;

            // Compute the accurate distance and locCoord
            for (int i = 0; i < bndElmtIds.size(); ++i)
            {
                // Iterate on each of the possible bnd elmt to get the locCoord
                // and the smallest dist
                bndGeom = BndExp[0]->GetExp(bndElmtIds[i])->GetGeom();

                isConverge = LocCoordForNearestPntOnBndElmt(
                    inGloCoord, bndGeom, locCoord_tmp, gloCoord_tmp, dist_tmp,
                    iterTol, 51);

                if (!isConverge)
                {
                    WARNINGL1(false, "Bisection iteration is not converged!!!");
                }

                // Perpendicular check (angle check)
                // When activated (distance > geoTol)
                // if the pnt-to-pnt vector is not parallel with the wall
                // normal, this bnd will be skipped. It is helpful to avoid
                // confusion when there are cornors
                if (isCheckAngle && (dist_tmp > geoTol)) // skip pts on the bnd
                {
                    // Generate the scaled vector
                    Vmath::Vcopy(nCoordDim, gloCoord_tmp, 1, normalChk, 1);
                    Vmath::Neg(nCoordDim, normalChk, 1);
                    Vmath::Vadd(nCoordDim, inGloCoord, 1, normalChk, 1,
                                normalChk, 1);
                    Vmath::Smul(nCoordDim,
                                1.0 / sqrt(Vmath::Dot(nCoordDim, normalChk, 1,
                                                      normalChk, 1)),
                                normalChk, 1, normalChk, 1);

                    // Get normal based on the locCoord_tmp
                    wnd.GetNormals(bndGeom, locCoord_tmp, normal);

                    // Reverse the direction to make it point into the field
                    bnd_phys_offset = BndExp[0]->GetPhys_Offset(bndElmtIds[i]);
                    for (int j = 0; j < nCoordDim; ++j)
                    {
                        normalRef[j] = normalQ[j][bnd_phys_offset];
                    }

                    if (Vmath::Dot(nCoordDim, normal, normalRef) > 0.0)
                    {
                        Vmath::Neg(nCoordDim, normal, 1);
                    }

                    // Check if the current bnd elmt convers this inner pnt
                    if (Vmath::Dot(nCoordDim, normalChk, 1, normal, 1) <
                        (1 - dirTol))
                    {
                        // If not covered, save the nearest result as backup
                        // Then skip
                        if (dist_tmp < dist_bak)
                        {
                            bndElmtId_bak = bndElmtIds[i];
                            dist_bak      = dist_tmp;

                            for (int j = 0; j < nBndLcoordDim; ++j)
                            {
                                locCoord_bak[j] = locCoord_tmp[j];
                            }

                            for (int j = 0; j < 3; ++j)
                            {
                                gloCoord_bak[j] = gloCoord_tmp[j];
                            }
                        }

                        continue; // Skip current bnd elmt
                    }
                }

                // No violation occurs, update the current result
                // This section must be executed at least once if there the
                // angle check is turned off. It means the following
                // (bndElmtId == -1) can only be caused by angle check.
                if (dist_tmp < distance[phys_offset + pId])
                {
                    bndElmtId                   = bndElmtIds[i];
                    distance[phys_offset + pId] = dist_tmp;

                    for (int j = 0; j < nBndLcoordDim; ++j)
                    {
                        locCoord[j] = locCoord_tmp[j];
                    }

                    for (int j = 0; j < 3; ++j)
                    {
                        gloCoord[j] = gloCoord_tmp[j];
                    }
                }
            }

            // Check if the bnd elmt is found. If not, use nearest one.
            if (bndElmtId == -1)
            {
                if (bndElmtIds.size() == 0)
                {
                    ASSERTL0(false, "No boundary element is found to be the \
                                     possible nearest element. Please check.");
                }

                distance[phys_offset + pId] = dist_bak;
                bndElmtId                   = bndElmtId_bak;

                for (int i = 0; i < nBndLcoordDim; ++i)
                {
                    locCoord[i] = locCoord_bak[i];
                }

                for (int i = 0; i < 3; ++i)
                {
                    gloCoord[i] = gloCoord_bak[i];
                }

                WARNINGL1(false,
                          "The boundary element is not found under given \
                                  tolerence. Please check and reset the \
                                  tolerence, or just turn off the angle check. \
                                  If you would like to continue, the nearest \
                                  boundary element is used.");
            }

            // Get the wall normal using the function in wallNormalData class
            bndGeom = BndExp[0]->GetExp(bndElmtId)->GetGeom();
            wnd.GetNormals(bndGeom, locCoord, normal);

            // Get ref normal for outward-point dir
            bnd_phys_offset = BndExp[0]->GetPhys_Offset(bndElmtId);
            for (int i = 0; i < nCoordDim; ++i)
            {
                normalRef[i] = normalQ[i][bnd_phys_offset];
            }

            // Correct the normal direction to make sure it points inward
            if (Vmath::Dot(nCoordDim, normal, normalRef) > 0.0)
            {
                Vmath::Neg(nCoordDim, normal, 1);
            }

            // The main tagential direction is defined to be overlapped with the
            // intersection line of the tangantial plane and the assistant
            // plane, whose normal is the input assistDir. Both plane is defined
            // using the point normal form, where the point is the nearest point
            // on the boundary (gloCoord), and the vectors are the normal
            // (pointing out) and the assistant direction (assistDir).
            // Therefore, the main tangantial direction can be obtained by the
            // cross product of the two vectors, since the two vectors has the
            // same starting point and the tangantial direction is perpendicular
            // to both of them.
            ScaledCrosssProduct(normal, assistVec, tangential1);

            // Save the direction vectors
            for (int i = 0; i < nSpaceDim; ++i)
            {
                bfcsDir[0][i][phys_offset + pId] = tangential1[i];
                bfcsDir[1][i][phys_offset + pId] = normal[i];
            }

            if (nSpaceDim == 3)
            {
                ScaledCrosssProduct(tangential1, normal, tangential2);
                for (int i = 0; i < 3; ++i)
                {
                    bfcsDir[2][i][phys_offset + pId] = tangential2[i];
                }
            }

        } // end of inner pts loop

    } // end of inner elmt loop
}

/**
 * @brief Get velocity and convert to Cartesian system, if it is still in
 *        transformed system. It is copied and modified from from
 *        ProcessGrad.cpp
 */
void ProcessBodyFittedVelocity::GetVelAndConvertToCartSys(
    Array<OneD, Array<OneD, NekDouble>> &vel)
{
    int spacedim   = m_f->m_exp[0]->GetCoordim(0) + m_f->m_numHomogeneousDir;
    int nfields    = m_f->m_variables.size();
    int npoints    = m_f->m_exp[0]->GetNpoints();
    int var_offset = 0;

    // Check type of the fields and set variable offset
    vector<string> vars = m_f->m_exp[0]->GetSession()->GetVariables();
    if (boost::iequals(vars[0], "rho") && boost::iequals(vars[1], "rhou"))
    {
        var_offset = spacedim + 2;
    }

    // Get mapping
    GlobalMapping::MappingSharedPtr mapping = ProcessMapping::GetMapping(m_f);

    // Get velocity and convert to Cartesian system,
    //      if it is still in transformed system
    if (m_f->m_fieldMetaDataMap.count("MappingCartesianVel"))
    {
        if (m_f->m_fieldMetaDataMap["MappingCartesianVel"] == "False")
        {
            // Initialize arrays and copy velocity
            for (int i = 0; i < spacedim; ++i)
            {
                vel[i] = Array<OneD, NekDouble>(npoints);
                if (m_f->m_exp[0]->GetWaveSpace())
                {
                    m_f->m_exp[0]->HomogeneousBwdTrans(
                        npoints, m_f->m_exp[var_offset + i]->GetPhys(), vel[i]);
                }
                else
                {
                    Vmath::Vcopy(npoints, m_f->m_exp[var_offset + i]->GetPhys(),
                                 1, vel[i], 1);
                }
            }
            // Convert velocity to cartesian system
            mapping->ContravarToCartesian(vel, vel);
            // Convert back to wavespace if necessary
            if (m_f->m_exp[0]->GetWaveSpace())
            {
                for (int i = 0; i < spacedim; ++i)
                {
                    m_f->m_exp[0]->HomogeneousFwdTrans(npoints, vel[i], vel[i]);
                }
            }
        }
        else
        {
            for (int i = 0; i < spacedim; ++i)
            {
                vel[i] = Array<OneD, NekDouble>(npoints);
                Vmath::Vcopy(npoints, m_f->m_exp[var_offset + i]->GetPhys(), 1,
                             vel[i], 1);
            }
        }
    }
    else
    {
        for (int i = 0; i < spacedim && i < nfields; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(npoints);
            Vmath::Vcopy(npoints, m_f->m_exp[var_offset + i]->GetPhys(), 1,
                         vel[i], 1);
        }
    }
}

} // namespace FieldUtils
} // namespace Nektar
