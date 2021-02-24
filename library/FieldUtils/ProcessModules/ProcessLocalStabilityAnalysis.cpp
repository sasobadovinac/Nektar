////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessLocalStabilityAnalysis.cpp
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
//  Description: Local linear stability analysis of compressible flow (for now).
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <FieldUtils/Interpolator.h>

#include "ProcessLocalStabilityAnalysis.h"

using namespace std;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessLocalStabilityAnalysis::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "localStabilityAnalysis"),
    ProcessLocalStabilityAnalysis::create,
    "Computes wall shear stress field.");

ProcessLocalStabilityAnalysis::ProcessLocalStabilityAnalysis(FieldSharedPtr f) : ProcessBoundaryExtract(f)
{
    f->m_writeBndFld = false; // turned on in upstream ProcessBoundaryExtract
    m_config["x"]    = ConfigOption(false, "0.1", 
                       "Sampling position given by x-coordinate.");
    m_config["y"]    = ConfigOption(false, "0.0",
                       "Sampling position given by y-coordinate.");
    m_config["z"]    = ConfigOption(false, "0.0",
                       "Sampling position given by z-coordinate.");
    m_config["tol"]  = ConfigOption(false, "0.1", 
                       "Relative tolerence to find the origin");
    m_config["useY"] = ConfigOption(true, "0", 
                       "Use y-coordinate to set up the position");
    m_config["h"]    = ConfigOption(false, "0.01", 
                       "Sampling distance along the wall normals.");
    m_config["nh"]   = ConfigOption(false, "11", 
                       "Number of sampling points along the wall normals.");
    m_config["d"]    = ConfigOption(false, "0.1", 
                       "Points distribution control in h direction, in (0,1)");
}

ProcessLocalStabilityAnalysis::~ProcessLocalStabilityAnalysis()
{
}

void ProcessLocalStabilityAnalysis::Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::Process(vm);

    /*Note
    * This module Should be able to process the LST for a given position
    * The x-direction is set to be the chordwise direction
    * The y-direction is set to be the vertical direction
    * The z-direction is set to be the spanwise direction
    *
    * For the cases to process the analysis on, they can be 2D, 2.5D and 3D. 
    * The 2D cases are normal. 
    * The 2.5D and 3D cases considered homogeneous in the z-direction.
    * However, we also allow the users to set origZ to be not zero
    */

    // Initialize sampling parameters
    Array<OneD, NekDouble> orig(3); // gloCoord of the origin
    orig[0]                   = m_config["x"].as<NekDouble>();
    orig[1]                   = m_config["y"].as<NekDouble>();
    orig[2]                   = m_config["z"].as<NekDouble>();
    const NekDouble relTol    = m_config["tol"].as<NekDouble>(); //
    const bool      isUseY    = m_config["useY"].as<bool>();     //
    const NekDouble distanceH = m_config["h"].as<NekDouble>();   //
    const int       nptsH     = m_config["nh"].as<int>();        //  
    const NekDouble delta     = m_config["d"].as<NekDouble>(); //0.1  
   
    // Get space dim
    const int nfields   = m_f->m_variables.size();              // number of fields
    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0);         // =2 for 2.5D cases
    m_spacedim          = nCoordDim + m_f->m_numHomogeneousDir; // overall dim of the input case
    const int nBndLcoordDim = 
        (m_spacedim==3 && m_f->m_numHomogeneousDir==0) ? 2 : 1; // local coordinate dim
    const int totVars   = m_spacedim + m_f->m_variables.size();

    // Get bnd info
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
    int bnd = BndRegionMap[m_f->m_bndRegionsToWrite[0]]; // Refer to wnd module


    // Get expansion list for boundary and the number of points and elements
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp(nfields); 
    for (int i = 0; i < nfields; ++i) {
        BndExp[i] = m_f->m_exp[i]->UpdateBndCondExpansion(bnd);
    }
    const int nqb    = BndExp[0]->GetTotPoints(); // points for all HomModesZ planes
    const int nElmts = BndExp[0]->GetNumElmts();  // number of elements on the boundary

    // remove this part later
    // Get inward-pointing wall-normal vectors for all quadrature points on bnd
    Array<OneD, Array<OneD, NekDouble> > normalsQ; 
    m_f->m_exp[0]->GetBoundaryNormals(bnd, normalsQ);
    for (int i = 0; i < m_spacedim; ++i) {
        Vmath::Neg(nqb, normalsQ[i], 1);
    }

    //-------------------------------------------------------------------------
    // Find the element that includes the origin
    SpatialDomains::GeometrySharedPtr bndGeom; 
    StdRegions::StdExpansionSharedPtr bndXmap;
    Array<OneD, const NekDouble> Jac;
    NekDouble scaledTol;
    int npts;
    Array<OneD, Array<OneD, NekDouble> >       pts(nCoordDim);
    Array<OneD, Array<OneD, const NekDouble> > bndCoeffs(nCoordDim);
    Array<OneD, NekDouble> locCoord(nCoordDim-1, -999.0);
    NekDouble resid;
    int elmtid;
    bool isInside = false;

    // Search and get precise Lcoord
    for (elmtid=0; elmtid<nElmts; ++elmtid) //nElmts
    {    
        bndGeom = BndExp[0]->GetExp(elmtid)->GetGeom(); 
        bndXmap = bndGeom->GetXmap();
        
        // Get the scaled tol by 2 times averaged Jacobian ([-1,1] -> [0,1]) 
        Jac = bndGeom->GetMetricInfo()->GetJac(bndXmap->GetPointsKeys());
        scaledTol = Vmath::Vsum(Jac.size(), Jac, 1)/((NekDouble)Jac.size());      //average Jac
        scaledTol = pow(Jac[0], 1.0/(static_cast<NekDouble>(nCoordDim)-1.0))*2.0; //length scale
        scaledTol *= relTol;

        // Get the coords for the element to check if origin is located inside
        /*Note
        * Currently (19 Feb 2021) 
        * SegGeom::v_GetLocCoords(...) can only get the Lcoord for 
        * straight-sided (eRegular) segment;
        * https://doc.nektar.info/doxygen/latest/_seg_geom_8cpp_source.html#l00372
        * QuadGeom::v_GetLocCoords(...) can only find Lcoord in the x-y plane
        * https://doc.nektar.info/doxygen/latest/_quad_geom_8cpp_source.html#l00481
        * TriGeom::v_GetLocCoords(...) can only find Lcoord in the x-y plane
        * https://doc.nektar.info/doxygen/latest/_tri_geom_8cpp_source.html#l00546
        * The above three possible types of boundary element do not include 
        * the needed routines to get Lcoord on curved boundary elements
        */

        npts = bndXmap->GetTotPoints();
        for (int i=0; i<nCoordDim; ++i) 
        {
            pts[i] = Array<OneD, NekDouble>(npts);
            bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
            bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
        }

        isInside = BndElmtContainsPoint(bndGeom, orig, locCoord, isUseY, scaledTol, resid);

        if (isInside) 
        {
            break;
        }

    }
    ASSERTL0(isInside, "Failed to find the sampling position."); 

    cout << "new gloCoord !!! " << endl;


    // Update the precise sampling position
    // x is precise, update y, or vice versa
    if(!isUseY)
    {
        orig[1] = bndXmap->PhysEvaluate(locCoord, pts[1]);                 
    }
    else
    {
        orig[0] = bndXmap->PhysEvaluate(locCoord, pts[0]);
    }

    // Update z according to the closest plane in 2.5D cases
    if(m_f->m_numHomogeneousDir==1)
    {
        int nPlanes    = m_f->m_exp[0]->GetHomogeneousBasis()->GetZ().size();
        NekDouble lHom = m_f->m_exp[0]->GetHomoLen();
        if (orig[2]<0.0 || orig[2]>lHom)
        {
            orig[2] = 0.0;
        }
        else
        {
            NekDouble dZ = lHom / nPlanes;
            NekDouble zTmp, zCur=0.0, distTmp, distCur = 999.0;
            for(int i=0; i<=nPlanes; ++i)
            {
                zTmp    = dZ * i;
                distTmp = abs(orig[2] - zTmp); 
                if(distTmp < distCur)
                {
                    distCur = distTmp;
                    zCur    = zTmp;
                }
            }
            orig[2] = zCur;
        }
    }
    
    cout << "Ox = " << orig[0] << ", Oy = " << orig[1] << ", Oz = " << orig[2] << endl;
    //-----------------------------------------------------------------------

    // Get the normals on the quadrature points in the element as reference
    // Set point key to get points
    Array<OneD, LibUtilities::PointsKey> from_key(nBndLcoordDim);
    int from_nPtsPerElmt = 1;
    for (int i=0; i<nBndLcoordDim; ++i)
    {
        from_key[i] = BndExp[0]->GetExp(elmtid)->GetBasis(i)->GetPointsKey();
        from_nPtsPerElmt *= from_key[i].GetNumPoints();
    }

    cout << "Ref normals:" << endl;
    for(int i=0;i<from_nPtsPerElmt;++i){
        cout << i << ", Q [nx,ny] = [" 
             << normalsQ[0][elmtid*from_nPtsPerElmt+i] << ", " 
             << normalsQ[1][elmtid*from_nPtsPerElmt+i] << "]" << endl;
    }

    // Correct the direction of the normals
    Array< OneD, NekDouble > normals(3, 0.0);
    GetNormals(bndGeom, locCoord, normals);
    Vmath::Neg(3, normals, 1);

    
    cout << "Final normals:" << endl;
    cout << "[nx, ny] = [" << normals[0] << ", " << normals[1] << "]" << endl;
    //cout << "Ox = " << orig[0] << ", Oy = " << orig[1] << ", Oz = " << orig[2] << endl;
    //---------------------------------------------------------------------
    cout << "========================================================" << endl;
    cout << "nqb = " << nqb << ", nElmts = " << nElmts << endl;
    cout << "Base points num on bnd geo= " << bndXmap->GetBasis(0)->GetNumPoints() << endl;
    cout << "Quadrature points num on bnd for computating = " << from_nPtsPerElmt << endl;
    cout << "scaledTol = " << scaledTol << endl;

    //-------------------------------------------------------------------------
    // Important block
    /*
    cout << "nqb = " << nqb << ", nElmts = " << nElmts << endl;
    cout << "bndExp[0] size = " << BndExp[0]->GetExpSize() << endl;    
    
    int elmtid = 0;
    SpatialDomains::GeometrySharedPtr bndGeom = BndExp[0]->GetExp(elmtid)->GetGeom(); 
    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    cout << "bnd numBase = " << bndXmap->GetNumBases() << endl;
  

    Array<OneD, NekDouble> ksi(2, 0.);
    Array<OneD, NekDouble> eta(2, 0.);
    ksi[0] = 1.0;
    ksi[1] = 2.0;
    bndXmap->LocCoordToLocCollapsed(ksi, eta);

    cout << "ksi: " << ksi[0] <<", " << ksi[1] << endl; 
    cout << "eta: " << eta[0] <<", " << eta[1] << endl; 


  
    int npts = bndXmap->GetTotPoints();
    cout << "npts = " << npts << endl;


    Array<OneD, const LibUtilities::BasisSharedPtr> base = bndXmap->GetBase();
    cout << "base type = " << base[0]->GetBasisType() << endl;
    cout << "base size = " << base.size() << endl;
    cout << "base npts = " << base[0]->GetNumPoints() << endl; 


    Array<OneD, NekDouble> ptsx(npts), ptsy(npts);

    Array<OneD, const NekDouble> bndCoeffsX = bndGeom->GetCoeffs(0);// 0 for x
    Array<OneD, const NekDouble> bndCoeffsY = bndGeom->GetCoeffs(1);// 1 for y 
   
    bndXmap->BwdTrans(bndCoeffsX, ptsx);
    bndXmap->BwdTrans(bndCoeffsY, ptsy);
 
    cout <<"xy = "<<ptsx[0]<< ", " <<ptsx[1]<< ", " <<ptsy[0]<< ", " <<ptsy[1]<< endl;
    */

    /*Note
    * According to the link
    * https://doc.nektar.info/doxygen/latest/class_nektar_1_1_std_regions_1_1_std_expansion1_d.html#ac73d10bd5d0fe32fefbb18bd20640a81
    * StdExpansion1D::v_PhysEvaluate uses Lcoord and calls Blas routine to compute. 
    * However, in the current master branch, StdExpansion1D::v_PhysEvaluate calls 
    * StdExpansion::BaryEvaluate<0>(Lcoord[0], &physvals[0]), which calls 
    * inline NekDouble BaryEvaluate(const NekDouble &coord, const NekDouble *physvals)
    * in the library/StdRegions/StdExpansion.h
    * The variables do not match (Lcoord vs. coord)
    * But the results seem to be correct? 
    */
    
    // Important
    /*
    Array<OneD, NekDouble> eta2(1);
    eta2[0] = 0; // range [-1,1]
    //Array<OneD, DNekMatSharedPtr> I(1); // 2
    //I[0] = bndXmap->GetBasis(0)->GetI(eta2); //GetI needs a pointer as input
    NekDouble x_tmp;
    x_tmp = bndXmap->PhysEvaluate(eta2, ptsx);
    //x_tmp = bndXmap->PhysEvaluate(I, ptsx);
    cout << x_tmp << endl;
    */


    // Set parametric distance of sampling points
    // Expression in Agrawal's paper:
    // h = 1- tanh((1-ksi)*atanh(sqrt(1-delta)))/sqrt(1-delta), ksi in [0,1]
    Array<OneD, NekDouble> h(nptsH);
    NekDouble tmp1;
    const NekDouble tmp2 = 1.0/(static_cast<NekDouble>(nptsH)-1.0);
    const NekDouble tmp3 = sqrt(1.0-delta);
    const NekDouble tmp4 = atanh(tmp3);
    const NekDouble tmp5 = 1.0/tmp3;
    for (int i=0; i<nptsH; ++i)
    {
        tmp1 = 1.0 - i * tmp2; // tmp1 = 1-ksi
        h[i] = 1 - tanh(tmp1*tmp4)*tmp5;
        h[i] *= distanceH; // physical distance in normal direction
    }


    // Set pts coordinates and interpoate the data
    Array<OneD, Array<OneD, NekDouble> > ptsH(totVars);
    for (int i=0; i<totVars; ++i)
    {
        ptsH[i] = Array<OneD, NekDouble>(nptsH, 0.0);
    }


    for(int i=0; i<m_spacedim; ++i)
    {
        for(int j=0; j<nptsH; ++j)
        {
            ptsH[i][j] = orig[i] + h[j] * normals[i]; // x0+dist*nx 
        }
    }

    m_f->m_fieldPts = MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
        m_spacedim, m_f->m_variables, ptsH, LibUtilities::NullPtsInfoMap);

    Interpolator interp;
    interp.Interpolate(m_f->m_exp, m_f->m_fieldPts, NekConstants::kNekUnsetDouble);

    // C2 Peojection? [in python]
    // Rotate the velocity in the tangential/normal direction [in python]
    // Re-scale the data as the CoPSE3D_LST needed [in python]
    //cout << ptsH[0][1] << ", " << ptsH[1][1] << ", " << ptsH[2][1] << ", " << ptsH[3][1] << ", "
    //     << ptsH[4][1] << ", " << ptsH[5][1] << ", " << ptsH[6][1] << ", " << ptsH[7][1] << ", "
    //     << ptsH[8][1] << ", " << ptsH[9][1] << ", " << ptsH[10][1]<< ", " << ptsH[11][1] << ", "
    //     << ptsH[12][1] << ", "<< ptsH[13][1]<< endl;
    //
    //0.0190167, 0.0208434, 0.755373, 0.771503, 0.392795, 2.14044, 1.02135, 
    //0.520001, 0.657731, 0.893228, -0.00169613, 1.1041, 1.03805, -3.46444




    //------------------------------------------------------------------------- 
    // call lst routine, do nothing for now
    call_hello();
   
    //call_lst();   

}

/*
* Check if the boundary element contain an input point
*/
bool ProcessLocalStabilityAnalysis::BndElmtContainsPoint(
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array< OneD, const NekDouble > & gloCoord,
        Array< OneD, NekDouble > & locCoord,
        const bool isUseY, 
        const NekDouble geomTol,
        NekDouble & resid)
{

    const int MaxIterations = 51;     //51
    const NekDouble iterTol  = 1.e-8;  // absolute tolerence of position

    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0);     // =2 for 2.5D cases
    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    Array<OneD, Array<OneD, NekDouble> >        pts(nCoordDim);
    Array<OneD, Array<OneD, const NekDouble> >  bndCoeffs(nCoordDim);
    int npts = bndXmap->GetTotPoints();

    for (int i=0; i<nCoordDim; ++i) 
    {
        pts[i] = Array<OneD, NekDouble>(npts);
        bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
        bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
    }


    const int dir1 = isUseY ? 1 : 0; // Use x or y to determine the location
    const int dir2 = 1 - dir1;       // The other Lcoord to be computed

    if(nCoordDim==2)
    {
        // Check if the point is in the range 
        if(pts[dir1][0] <= gloCoord[dir1] && pts[dir1][npts-1] >= gloCoord[dir1])
        {
            // In this range, start iteration of bisection
            Array<OneD, NekDouble> etaLR(2); // range [-1,1]
            etaLR[0] = -1.0; // left
            etaLR[1] =  1.0; // right
            NekDouble tmpL, tmpR;
            bool isConverge = false;
            bool isDesired  = false;
            int cnt = 0;

            while(cnt<MaxIterations)
            {
                tmpL = bndXmap->PhysEvaluate(etaLR  , pts[dir1]);
                tmpR = bndXmap->PhysEvaluate(etaLR+1, pts[dir1]);

                if (abs(gloCoord[dir1]-tmpL) >= abs(gloCoord[dir1]-tmpR))
                {
                    etaLR[0] = 0.5 * (etaLR[0]+etaLR[1]);
                }
                else
                {
                    etaLR[1] = 0.5 * (etaLR[0]+etaLR[1]);
                }

                if ( (etaLR[1]-etaLR[0]) < iterTol)
                {
                    locCoord[0] = 0.5 * (etaLR[0]+etaLR[1]);
                    resid = abs(0.5*(tmpL+tmpR)-gloCoord[dir1]);
                    isConverge = true;
                    break;
                }

                ++cnt;
            }

            tmpL = bndXmap->PhysEvaluate(locCoord, pts[dir2]);
            if (abs(gloCoord[dir2]-tmpL) < geomTol)
            {
                isDesired = true;
            }

            return isConverge && isDesired;

        }
        else
        {
            return false;
        }
        
    }
    else
    {
        NekDouble angleAbs, angleSign, angleSum = 0.0, angleCos;
        const NekDouble angleTol = 1e-6;
        Array<OneD, NekDouble> vec1(2, 0.0), vec2(2, 0.0);
        int id1, id2;

        // Generate a polygen, re-order the points on the edge
        int nptsEdge = sqrt(npts);  // num of points on edge for geo representation
        Array<OneD, Array<OneD, NekDouble> > ptsPolygon(3);
        for (int i=0; i<3; ++i) 
        {
            ptsPolygon[i] = Array<OneD, NekDouble>(4*nptsEdge-4);

            for(int j=0; j<nptsEdge; ++j)
            {
                ptsPolygon[i][j] = pts[i][j];
            }
            for(int j=0; j<nptsEdge-2; ++j)
            {
                ptsPolygon[i][nptsEdge+j] = pts[i][(j+2)*nptsEdge-1];
            }
            for(int j=0; j<nptsEdge; ++j)
            {
                ptsPolygon[i][2*nptsEdge-2+j] = pts[i][npts-1-j];
            }
            for(int j=0; j<nptsEdge-2; ++j)
            {
                ptsPolygon[i][3*nptsEdge-2+j] = pts[i][nptsEdge*(nptsEdge-j-2)];
            }
        }


        // Determine relation using z-x (isUseY==0) or y-z (isUseY==1) 
        for(int i=0; i<4*nptsEdge-4; ++i)
        {
            id1 = i;
            id2 = (id1==(4*nptsEdge-5)) ? 0 : id1+1;

            vec1[0] = ptsPolygon[dir1][id1] - gloCoord[dir1];
            vec1[1] = ptsPolygon[2][id1]    - gloCoord[2];
            vec2[0] = ptsPolygon[dir1][id2] - gloCoord[dir1];
            vec2[1] = ptsPolygon[2][id2]    - gloCoord[2];
        
            angleSign = ((vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0) ? 1.0 : -1.0;
            angleCos  = (vec1[0]*vec2[0]+vec1[1]*vec2[1])/
                         sqrt( (vec1[0]*vec1[0]+vec1[1]*vec1[1])
                               *(vec2[0]*vec2[0]+vec2[1]*vec2[1]));
            if(angleCos>1.0)
            {
                angleCos=1.0;
            }
            else if(angleCos<-1.0)
            {
                angleCos=-1.0;
            }

            angleAbs = acos(angleCos);
            angleSum += angleSign * angleAbs;
        }
        angleSum = abs(angleSum);        
        //cout <<"angle = " << angleSum << endl;


        if( abs(angleSum-2.0*M_PI) < angleTol )
        {
            /*
            cout <<"inside." << endl;
            cout << "npts = " << npts << endl;
            for (int i=0; i<npts; ++i)
            {
                cout << i <<", old [x,y,z]=[" << pts[0][i]<<", " << pts[1][i]<<", "<< pts[2][i]<<"]" << endl;
            }
            cout << "gloCoordXYZ = " << gloCoord[0] << ", " << gloCoord[1] << ", " << gloCoord[2] << endl;
            */
            
            // Inside the element projection, use Newton iteration to find the Lcoord
            bool isConverge = false;
            bool isDesired  = false;        
            NekDouble dist  = 999.0;
            NewtonIterationForLocCoordOnBndElmt(bndGeom, gloCoord, pts[0], pts[1], pts[2], dir2, locCoord, dist); //1

            if(dist < geomTol)
            {
                isConverge = true; // converged to get Lcoord 
            }

            NekDouble gloCoordTmp = bndXmap->PhysEvaluate(locCoord, pts[dir2]);
            if (abs(gloCoord[dir2]-gloCoordTmp) < geomTol)
            {
                isDesired = true;
            }

            return isConverge && isDesired;            
        }
        else
        {
            return false; // not in the projected area
        }
 
    }

    return false;
}


// 3D case -> bnd element projected to 2D plane
void ProcessLocalStabilityAnalysis::NewtonIterationForLocCoordOnBndElmt(
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array<OneD, const NekDouble> &coords,
    const Array<OneD, const NekDouble> &ptsx,
    const Array<OneD, const NekDouble> &ptsy,
    const Array<OneD, const NekDouble> &ptsz,
    const int projDir,
    Array<OneD, NekDouble> &Lcoords,
    NekDouble &dist)
{
    // Iteration settings 
    // Maximum iterations for convergence
    const int MaxIterations = 51;
    // |x-xp|^2 < EPSILON  error    tolerance
    const NekDouble Tol = 1.e-8;
    // |r,s|    > LcoordDIV stop   the search
    const NekDouble LcoordDiv = 15.0;

    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();

    Array<OneD, const NekDouble> Jac =
        bndGeom->GetMetricInfo()->GetJac(bndXmap->GetPointsKeys());
    NekDouble ScaledTol = Vmath::Vsum(Jac.size(), Jac, 1) /
                          ((NekDouble)Jac.size());
    ScaledTol *= Tol;

    
    // Initial the locCoord
    Lcoords[0] = 0.0;
    Lcoords[1] = 0.0;

    NekDouble x1map, x2map, F1, F2;
    NekDouble derx1_1, derx1_2, derx2_1, derx2_2, jac;

    Array<OneD, const NekDouble> ptsx1, ptsx2;
    Array<OneD, NekDouble> LcoordsTmp(2);
    int dir[2]={-1, -1};
    if(projDir==0)
    {
        dir[0] = 1;
        dir[1] = 2;
        ptsx1 = ptsy;
        ptsx2 = ptsz;
    }
    else if(projDir==1)
    {
        dir[0] = 2;
        dir[1] = 0;
        ptsx1 = ptsz;
        ptsx2 = ptsx;     
    }
    else if(projDir==2)
    {
        dir[0] = 0;
        dir[1] = 1;
        ptsx1 = ptsx;
        ptsx2 = ptsy;
    }
    else{
        ASSERTL0(false, "The projection direction needs to be 0 or 1 or 2.");  
    }
     
    Array<OneD, NekDouble> Dx1D1(ptsx1.size());
    Array<OneD, NekDouble> Dx1D2(ptsx1.size());
    Array<OneD, NekDouble> Dx2D1(ptsx1.size());
    Array<OneD, NekDouble> Dx2D2(ptsx1.size());

    // Ideally this will be stored in m_geomfactors
    bndXmap->PhysDeriv(ptsx1, Dx1D1, Dx1D2);
    bndXmap->PhysDeriv(ptsx2, Dx2D1, Dx2D2);
    
    int cnt = 0;

    F1 = F2 = 2000; // Starting value of Function
    NekDouble resid;
    while (cnt++ < MaxIterations)
    {
        x1map = bndXmap->PhysEvaluate(Lcoords, ptsx1);
        x2map = bndXmap->PhysEvaluate(Lcoords, ptsx2);

        F1 = coords[dir[0]] - x1map;
        F2 = coords[dir[1]] - x2map;

        if (F1 * F1 + F2 * F2 < ScaledTol)
        {
            resid = sqrt(F1 * F1 + F2 * F2);
            break;
        }

        // Interpolate derivative metric at Lcoords
        derx1_1 = bndXmap->PhysEvaluate(Lcoords, Dx1D1);
        derx1_2 = bndXmap->PhysEvaluate(Lcoords, Dx1D2);
        derx2_1 = bndXmap->PhysEvaluate(Lcoords, Dx2D1);
        derx2_2 = bndXmap->PhysEvaluate(Lcoords, Dx2D2);

        jac = derx2_2 * derx1_1 - derx2_1 * derx1_2;
        
        // use analytical inverse of derivitives which are
        // also similar to those of metric factors.
        Lcoords[0] =
            Lcoords[0] +
            ( derx2_2 * (coords[dir[0]] - x1map) - derx1_2 * (coords[dir[1]] - x2map)) / jac;

        Lcoords[1] =
            Lcoords[1] +
            (-derx2_1 * (coords[dir[0]] - x1map) + derx1_1 * (coords[dir[1]] - x2map)) / jac;

        
        if( !(std::isfinite(Lcoords[0]) && std::isfinite(Lcoords[1])) )
        {
            dist = 1e16;
            std::ostringstream ss;
            ss << "nan or inf found in NewtonIterationForLocCoord in element "
               << bndGeom->GetGlobalID();
            WARNINGL1(false, ss.str());
            return;
        }
        if (fabs(Lcoords[0]) > LcoordDiv || fabs(Lcoords[1]) > LcoordDiv)
        {
            break; // lcoords have diverged so stop iteration
        }
    }

    Array<OneD, DNekMatSharedPtr> I(2);
    Array<OneD, NekDouble> eta(2);

    bndXmap->LocCoordToLocCollapsed(Lcoords, eta);
    if(bndGeom->ClampLocCoords(eta, 0.0))
    {
        I[0] = bndXmap->GetBasis(0)->GetI(eta);
        I[1] = bndXmap->GetBasis(1)->GetI(eta + 1);
        // calculate the global point corresponding to Lcoords
        x1map = bndXmap->PhysEvaluate(I, ptsx1);
        x2map = bndXmap->PhysEvaluate(I, ptsx2);
        F1 = coords[dir[0]] - x1map;
        F2 = coords[dir[1]] - x2map;
        dist = sqrt(F1 * F1 + F2 * F2);
    }
    else
    {
        dist = 0.0;
    }

    if (cnt >= MaxIterations)
    {
        Array<OneD, NekDouble> collCoords(2);
        bndXmap->LocCoordToLocCollapsed(LcoordsTmp, collCoords);

        // if coordinate is inside element dump error!
        if ((collCoords[0] >= -1.0 && collCoords[0] <= 1.0) &&
            (collCoords[1] >= -1.0 && collCoords[1] <= 1.0))
        {
            std::ostringstream ss;

            ss << "Reached MaxIterations (" << MaxIterations
               << ") in Newton iteration ";
            ss << "Init value (" << setprecision(4) << 0 << "," << 0
               << ","
               << ") ";
            ss << "Fin  value (" << LcoordsTmp[0] << "," << LcoordsTmp[1] << ","
               << ") ";
            ss << "Resid = " << resid << " Tolerance = " << sqrt(ScaledTol);

            WARNINGL1(cnt < MaxIterations, ss.str());
        }
    }
    
    /*
    cout << "dir1, dir2 = [" << dir[0] << ", " << dir[1] << "]" << endl;
    cout << "Lcoords = [" << Lcoords[0] << ", " << Lcoords[1] << ", " << Lcoords[2] << "]" << endl;
    cout << "xyz = " << bndXmap->PhysEvaluate(Lcoords, ptsx1) << ", " 
                     << bndXmap->PhysEvaluate(Lcoords, ptsx2) << endl;
    */
}


/*
* Get the normals for a given locCoord
*/
void ProcessLocalStabilityAnalysis::GetNormals(
    SpatialDomains::GeometrySharedPtr bndGeom,
    Array< OneD, NekDouble > & locCoord, 
    Array< OneD, NekDouble > & normals)
{
    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0);     // =2 for 2.5D cases
    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    Array<OneD, Array<OneD, NekDouble> >        pts(m_spacedim);
    Array<OneD, Array<OneD, const NekDouble> >  bndCoeffs(m_spacedim);
    int npts = bndXmap->GetTotPoints();

    // Get pts in the element
    for (int i=0; i<nCoordDim; ++i) 
    {
        pts[i] = Array<OneD, NekDouble>(npts);
        bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
        bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
    }

    // Get the outward-pointing normals according to the given locCoord
    if(nCoordDim==2)
    {
        Array<OneD, NekDouble> DxD1(pts[0].size());
        Array<OneD, NekDouble> DyD1(pts[0].size());

        bndXmap->PhysDeriv(pts[0], DxD1);
        bndXmap->PhysDeriv(pts[1], DyD1);

        NekDouble dxd1, dyd1, tmp;
        dxd1 = bndXmap->PhysEvaluate(locCoord, DxD1);
        dyd1 = bndXmap->PhysEvaluate(locCoord, DyD1);
        tmp = sqrt(dxd1*dxd1 + dyd1*dyd1);

        normals[0] =  dyd1/tmp;
        normals[1] = -dxd1/tmp;
        normals[2] =  0.0;
    }
    else
    {
        // Will be developed soon
        ASSERTL0(false, "Analysis for full 3D cases is currently unavailable");
    }

}

/*
Case 1: 2D channel
/disk_two/Nek_Test/nektar++/build_f90/dist/bin/FieldConvert -f -m localStabilityAnalysis:bnd=2:x=0.10472:y=1 test_1.xml test_1.fld test_1.pts

Case 2: 2D airfoil
/disk_two/Nek_Test/nektar++/build_f90/dist/bin/FieldConvert -f -m localStabilityAnalysis:bnd=0:x=0.02:y=0.0195:useY=1:h=0.02:nh=5 test_2_mesh.xml test_2_session.xml test_2.fld test_2.pts

Case3: 2.5D flat plate
/disk_two/Nek_Test/nektar++/build_f90/dist/bin/FieldConvert -f -m localStabilityAnalysis:bnd=0:x=0.335:y=0.0001:z=0.007:h=0.002:nh=5 test_3_mesh.xml test_3_session.xml test_3.fld test_3.pts

Case4: 3D flat plate
/disk_two/Nek_Test/nektar++/build_f90/dist/bin/FieldConvert -f -m localStabilityAnalysis:bnd=0:x=0.601:y=0.001:h=0.02:nh=5 test_4_mesh.xml test_4_session.xml test_4.fld test_4.pts
*/



}
}
