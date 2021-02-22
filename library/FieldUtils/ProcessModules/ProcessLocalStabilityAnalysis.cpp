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

#include "ProcessLocalStabilityAnalysis.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>


#include <SpatialDomains/PointGeom.h>

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
                       "Sampling position given by x-coordinate.");
    m_config["tol"]  = ConfigOption(false, "0.1", 
                       "Relative tolerence to find the origin");

    m_config["useY"] = ConfigOption(true, "0", 
                       "Use y-coordinate to set up the position");
    m_config["h"]    = ConfigOption(false, "0.01", 
                       "Sampling distance along the wall normals.");
    m_config["nh"]   = ConfigOption(false, "11", 
                       "Number of sampling points along the wall normals.");

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
    //const NekDouble origX      = m_config["x"].as<NekDouble>();   //
    //const NekDouble origY      = m_config["y"].as<NekDouble>();   //
    //const NekDouble origZ      = m_config["z"].as<NekDouble>();   //
    Array<OneD, NekDouble> orig(3);
    orig[0] = m_config["x"].as<NekDouble>();
    orig[1] = m_config["y"].as<NekDouble>();
    orig[2] = m_config["z"].as<NekDouble>();
    const NekDouble relTol     = m_config["tol"].as<NekDouble>(); //
    const bool      isUseY     = m_config["useY"].as<bool>();     //
    const NekDouble distance_h = m_config["h"].as<NekDouble>();   //
    const int       npts_h     = m_config["nh"].as<int>();        //  
    
    cout << distance_h << ", " << npts_h << endl;

    
    // Get space dim
    const int nfields   = m_f->m_variables.size();              // number of fields
    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0);         // =2 for 2.5D cases
    m_spacedim          = nCoordDim + m_f->m_numHomogeneousDir; // overall dim of the input case
    const int nBndLcoordDim = 
        (m_spacedim==3 && m_f->m_numHomogeneousDir==0) ? 2 : 1; // local coordinate dim

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
    Array<OneD, Array<OneD, NekDouble> >       pts(m_spacedim);
    Array<OneD, Array<OneD, const NekDouble> > bndCoeffs(m_spacedim);
    Array<OneD, NekDouble> gloCoord(m_spacedim), locCoord(m_spacedim, -999.0);
    NekDouble resid;
    int elmtid;
    bool isInside = false;
    
    for (int i=0; i<m_spacedim; ++i)
    {
        gloCoord[i] = orig[i];
    }

    // Search and get precise Lcoord
    for (elmtid=0; elmtid<nElmts; ++elmtid)
    {    
        bndGeom = BndExp[0]->GetExp(elmtid)->GetGeom(); 
        bndXmap = bndGeom->GetXmap();
        
        // Get the scaled tol by 2 times averaged Jacobian ([-1,1] -> [0,1]) 
        Jac = bndGeom->GetMetricInfo()->GetJac(bndXmap->GetPointsKeys());
        scaledTol = Vmath::Vsum(Jac.size(), Jac, 1)*2./((NekDouble)Jac.size());
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
        for (int i=0; i<m_spacedim; ++i) 
        {
            pts[i] = Array<OneD, NekDouble>(npts);
            bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
            bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
        }
        //cout <<"xy = "<<pts[0][0]<< ", " <<pts[0][1]<< ", " <<pts[1][0]<< ", " <<pts[1][1]<< endl;
        //cout << "eid = " << elmtid << endl;
        isInside = BndElmtContainsPoint(bndGeom, gloCoord, locCoord, isUseY, scaledTol, resid);


        if (isInside) 
        {
            break;
        }

    }
    ASSERTL0(isInside, "Failed to find the sampling position."); 

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
    NekDouble normalDir = 0.0;
    for (int i=0; i<m_spacedim; ++i)
    {
        normalDir += normalsQ[i][elmtid*from_nPtsPerElmt]*normals[i];
    }

    if (normalDir<0)
    {
        Vmath::Neg(3, normals, 1);
    }
    
    cout << "Final normals:" << endl;
    cout << "[nx, ny] = [" << normals[0] << ", " << normals[1] << "]" << endl;
    cout << "Ox = " << orig[0] << ", Oy = " << orig[1] << ", Oz = " << orig[2] << endl;
    //---------------------------------------------------------------------
    cout << "========================================================" << endl;
    cout << "nqb = " << nqb << ", nElmts = " << nElmts << endl;
    cout << "Base points num on bnd geo= " << bndXmap->GetBasis(0)->GetNumPoints() << endl;
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
 
    cout << "exp size = " << m_f->m_exp[0]->GetExpSize() << endl;
    SpatialDomains::GeometrySharedPtr geom = m_f->m_exp[0]->GetExp(0)->GetGeom(); 
    StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();
    cout << "numBase = " << xmap->GetNumBases() << endl; 

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
        const NekDouble tol,
        NekDouble & resid)
{
    
    ASSERTL0(m_spacedim > 1,
             "Space dimension should be 2 or 3 for local stability analysis");

    const int MaxIterations = 51;     //51
    const NekDouble AbsTol  = 1.e-8;  // absolute tolerence of position

    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    Array<OneD, Array<OneD, NekDouble> >        pts(m_spacedim);
    Array<OneD, Array<OneD, const NekDouble> >  bndCoeffs(m_spacedim);
    int npts = bndXmap->GetTotPoints();

    for (int i=0; i<m_spacedim; ++i) 
    {
        pts[i] = Array<OneD, NekDouble>(npts);
        bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
        bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
    }


    const int dir1 = isUseY ? 1 : 0; // Use x or y to determine the location
    const int dir2 = 1 - dir1;       // The other Lcoord to be computed

    if(m_spacedim==2)
    {
        // Check if the point is in the range 
        if(pts[dir1][0] <= gloCoord[dir1] && pts[dir1][npts-1] >= gloCoord[dir1])
        {
            // In this range, start iteration
            Array<OneD, NekDouble> etaLR(2); // range [-1,1]
            etaLR[0] = -1; // left
            etaLR[1] =  1; // right
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

                if ( (etaLR[1]-etaLR[0]) < AbsTol)
                {
                    locCoord[0] = 0.5 * (etaLR[0]+etaLR[1]);
                    resid = abs(0.5*(tmpL+tmpR)-gloCoord[dir1]);
                    isConverge = true;
                    break;
                }

                ++cnt;
            }

            tmpL = bndXmap->PhysEvaluate(locCoord, pts[dir2]);
            if (abs(gloCoord[dir2]-tmpL) < tol)
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
    else if (m_spacedim==3 && m_f->m_numHomogeneousDir==1)
    {
        // Will be developed soon
        ASSERTL0(false, "Analysis for 2.5D cases is currently unavailable");  
    }
    else
    {
        // Will be developed soon
        ASSERTL0(false, "Analysis for full 3D cases is currently unavailable");  
    }



    return false;
}

/*
* Get the normals for a given locCoord
*/
void ProcessLocalStabilityAnalysis::GetNormals(
    SpatialDomains::GeometrySharedPtr bndGeom,
    Array< OneD, NekDouble > & locCoord, 
    Array< OneD, NekDouble > & normals)
{
    ASSERTL0(m_spacedim > 1,
             "The dimension for boundary element should be 2 or 3.");

    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    Array<OneD, Array<OneD, NekDouble> >        pts(m_spacedim);
    Array<OneD, Array<OneD, const NekDouble> >  bndCoeffs(m_spacedim);
    int npts = bndXmap->GetTotPoints();

    // Get pts in the element
    for (int i=0; i<m_spacedim; ++i) 
    {
        pts[i] = Array<OneD, NekDouble>(npts);
        bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
        bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
    }

    

    // Get the normals according to the given locCoord
    if(m_spacedim==2)
    {
        Array<OneD, NekDouble> DxD1(pts[0].size());
        Array<OneD, NekDouble> DyD1(pts[0].size());

        bndXmap->PhysDeriv(pts[0], DxD1);
        bndXmap->PhysDeriv(pts[1], DyD1);


        NekDouble dxd1, dyd1, tmp;
        dxd1 = bndXmap->PhysEvaluate(locCoord, DxD1);
        dyd1 = bndXmap->PhysEvaluate(locCoord, DyD1);
        tmp = sqrt(dxd1*dxd1 + dyd1*dyd1);
        
        normals[0] = -dyd1/tmp;
        normals[1] =  dxd1/tmp;
        //cout << "![nx, ny] = [" << normals[0] << ", " << normals[1] << "]" << endl;
    }
    else
    {
        // Will be developed soon
        ASSERTL0(false, "Analysis for full 3D cases is currently unavailable");
        /*
        Array<OneD, NekDouble> DxD1(ptsx.size());
        Array<OneD, NekDouble> DxD2(ptsx.size());
        Array<OneD, NekDouble> DyD1(ptsx.size());
        Array<OneD, NekDouble> DyD2(ptsx.size());

        bndXmap->PhysDeriv(ptsx, DxD1, DxD2);
        bndXmap->PhysDeriv(ptsy, DyD1, DyD2);

        m_xmap->LocCoordToLocCollapsed(Lcoords, eta);

        */
    }

}



}
}
