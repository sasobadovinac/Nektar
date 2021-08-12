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
//  Description: Get the wall-normal data at a given origin.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <FieldUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>

#include "ProcessWallNormalData.h"
#include "ProcessBodyFittedVelocity.h"


using namespace std;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessBodyFittedVelocity::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "bodyFittedVelocity"),
    ProcessBodyFittedVelocity::create,
    "Convert the velocity components from the Cartesian coordinate into the body-fitted coordinate.");

ProcessBodyFittedVelocity::ProcessBodyFittedVelocity(FieldSharedPtr f) : ProcessBoundaryExtract(f)
{
    f->m_writeBndFld = false; // turned on in the father class ProcessBoundaryExtract

    // bnd is read from father class ProcessBoundaryExtract
    m_config["fixedDir"]  = ConfigOption(false, "z",
                            "The fixed direction at which the velocity \
                            component does not rotated. default=z");
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
    
    string fixedDir  = m_config["fixedDir"].as<string>();

    if (!(boost::iequals(fixedDir, "x") || boost::iequals(fixedDir, "X") ||
          boost::iequals(fixedDir, "y") || boost::iequals(fixedDir, "Y") ||
          boost::iequals(fixedDir, "z") || boost::iequals(fixedDir, "Z") ) )
    {
        ASSERTL0(false, "fixedDir must be set as x, y or z, \
                         or the uppercase.");
    }

    // Get dim to store data
    const int nFields   = m_f->m_variables.size();
    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0);
    m_spacedim          = nCoordDim + m_f->m_numHomogeneousDir;
    const int nBndLcoordDim = nCoordDim - 1;
    const int nTotVars   = m_spacedim + nFields;

    cout <<nFields<<nCoordDim<<m_spacedim<<nBndLcoordDim<<nTotVars << endl;


    // Add var names
    m_f->m_variables.push_back("distance2Wall");
    /*
    if (m_spacedim == 2)
    {
        
        //m_f->m_variables.push_back("u_bfc");
        //m_f->m_variables.push_back("v_bfc");
    }
    else if (m_spacedim == 3)
    {
        //m_f->m_variables.push_back("u_bfc");
        //m_f->m_variables.push_back("v_bfc");
        //m_f->m_variables.push_back("w_bfc");
    }
    else
    {
        ASSERTL0(false, "Velocity in 1D case is already in the body fitted \
                         coordinate. The dimension should be 2 or 3.");
    }
    */



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

    // Get expansion list for boundary
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp(nFields); 
    for (int i = 0; i < nFields; ++i) {
        BndExp[i] = m_f->m_exp[i]->UpdateBndCondExpansion(bnd);
    }
    

    //-------------------------------------------------------------------------
    // Loop all the Q points of the field, and find out the closest point on
    // the bnd, and the bnd element contains this bnd point

    // Add the new fields (distance + u,v,w_bfc)
    m_f->m_exp.resize(nFields + 1);
    int npoints = m_f->m_exp[0]->GetNpoints();
    Array<OneD, NekDouble> distance(npoints, 9999);
    /*
    Array<OneD, Array<OneD, NekDouble> > vel_bfc(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        vel_bfc[i] = Array<OneD, NekDouble>(npoints, 0.0);
    }
    */



    const int nElmts    = m_f->m_exp[0]->GetNumElmts(); // num for domain element
    const int nBndElmts = BndExp[0]->GetNumElmts();     // num for bnd element

    cout << "nElmts = " << nElmts <<", nBndElmts = " << nBndElmts << endl;

    SpatialDomains::GeometrySharedPtr bndGeom; //geom
    StdRegions::StdExpansionSharedPtr bndXmap; //xmap
    int nPts, nBndPts;
    NekDouble dist_tmp;
    bool isConverge;
    int bndElmtId=-1;
    int phys_offset; // coef_offset = m_field->GetCoeff_Offset(i);

    // Loop inner element
    for (int eId=0; eId<nElmts; ++eId) //element id, nElmts
    {
        // Get inner points
        nPts = m_f->m_exp[0]->GetExp(eId)->GetTotPoints();
        Array<OneD, Array<OneD, NekDouble> > pts(nCoordDim);
        for (int i=0; i<nCoordDim; ++i) 
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

        // -------Debug------
        if (eId==0){
            for (int i=0; i<nPts; ++i){
                cout << "pId = " << i <<", [x,y]=["<<pts[0][i]<<", "<<pts[1][i]<<"]"<<endl;
            }
        }

        
        phys_offset = m_f->m_exp[0]->GetPhys_Offset(eId);
        bndElmtId   = -1;
        // Loop points in the inner element
        for (int pId=0; pId<nPts; ++pId)
        {

            // Loop bnd element to find the element to build body-fitted coordinate 
            for (int beId=0; beId<nBndElmts; ++beId)
            {
                // Get geomery and points
                bndGeom = BndExp[0]->GetExp(beId)->GetGeom();   // get geom of a bnd element
                bndXmap = bndGeom->GetXmap();
                nBndPts = bndXmap->GetTotPoints();

                Array<OneD, Array<OneD, const NekDouble> > bndCoeffs(nCoordDim);
                Array<OneD, Array<OneD, NekDouble> >       bndPts(nCoordDim);
                for (int i=0; i<nCoordDim; ++i) 
                {
                    bndPts[i]    = Array<OneD, NekDouble>(nBndPts);
                    bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
                    bndXmap->BwdTrans(bndCoeffs[i], bndPts[i]);
                }

                // Compute the distance from the pnt (bnd elmt) to the inner pnt
                // The elemt with the smallest distance is the one to place the 
                // body-fitted coordinate.  
                dist_tmp = PntToBndElmtPntDistance(pts, pId, bndPts);

                if (dist_tmp < distance[phys_offset+pId])
                {
                    distance[phys_offset+pId] = dist_tmp;
                    bndElmtId = beId;
                }

            } // end of bnd elmt loop  

            //cout << "pnt " << pId <<", dist = " << dist <<", beId = " << bndElmtId << endl;
            // Use iteration to find out the local coordinate and wall normal
            bndGeom = BndExp[0]->GetExp(bndElmtId)->GetGeom(); // Get the geom for target bnd elmt
            Array<OneD, NekDouble > normals(3, 0.0);
            Array<OneD, NekDouble > locCoord(nBndLcoordDim, -999.0);
            Array<OneD, NekDouble > gloCoord(3, 0.0);
            for (int i=0; i<nCoordDim; ++i) 
            {
                gloCoord[i] = pts[i][pId];
            }

            isConverge = LocCoordForNearestPntOnBndElmt( gloCoord, bndGeom,
                                    locCoord, normals, dist_tmp, 1.0e-8, 51);
            
            if (!isConverge)
            {
                WARNINGL1(false, "Bisection iteration is not converged!!!");
            }

            distance[phys_offset+pId] = dist_tmp; // Update a more accurate distance field

            // -------Debug------
            if (eId==0 && pId==12){
                cout << "  - pId = " << pId <<", [x,y]=["<<pts[0][pId]<<", "<<pts[1][pId]<<"]"<<endl;
                cout << "    - dist     = " << dist_tmp << ", locCoord = "<<locCoord[0]<<", bndElmtId = "<< bndElmtId <<endl;
                cout << "    - dist_fld = " << distance[phys_offset+pId] << endl;
            }

        } // end of inner pts loop
             
    } // end of inner elmt loop

    // ==========================================================================================
    

    
    
    int number = 0;
    for (int i=0;i<m_f->m_exp[0]->GetNumElmts(); ++i)
    {
        //number += m_f->m_exp[0]->GetExp(i)->GetGeom()->GetXmap()->GetTotPoints();
        number += m_f->m_exp[0]->GetExp(i)->GetTotPoints();
    }
    cout << "totNumPts_1 = " << number<<endl;
    cout << "totNumPts_2 = " << m_f->m_exp[0]->GetNpoints()<<endl;
    cout << " " << bndElmtId << endl;
    
    cout << (m_f->m_exp[0]->GetExp(1)->GetBase())[0]->GetNumPoints() << endl;        // 4
    cout << (m_f->m_exp[0]->GetExp(1)->GetBase())[1]->GetNumPoints() << endl;        // 4
    cout <<  m_f->m_exp[0]->GetExp(1)->GetTotPoints() << endl;                       // 16
    cout <<  m_f->m_exp[0]->GetExp(1)->GetGeom()->GetXmap()->GetTotPoints() << endl; // 15


    // copy distance
    m_f->m_exp[nFields] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
    Vmath::Vcopy(npoints, distance, 1, 
                     m_f->m_exp[nFields]->UpdatePhys(), 1);
    m_f->m_exp[nFields]->FwdTrans_IterPerExp(
            distance, m_f->m_exp[nFields]->UpdateCoeffs());
    /*
    // copy vel_bfc
    for (int i = 1; i < (m_spacedim+1); ++i)
    {
        m_f->m_exp[nFields + i] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
        Vmath::Vcopy(npoints, vel_bfc[i], 1, 
                     m_f->m_exp[nFields + i]->UpdatePhys(), 1);
        m_f->m_exp[nFields + i]->FwdTrans_IterPerExp(
            vel_bfc[i], m_f->m_exp[nFields + i]->UpdateCoeffs());
    }
    */



    /*
    SpatialDomains::GeometrySharedPtr geom, bndGeom;
    geom     = m_f->m_exp[0]->GetExp(140)->GetGeom(); // id=0
    bndGeom  = BndExp[0]->GetExp(0)->GetGeom();     //  bndElmtId=0

    StdRegions::StdExpansionSharedPtr xmap    = geom->GetXmap();
    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();

    const int nPts     = xmap->GetTotPoints();
    const int nBndPts  = bndXmap->GetTotPoints();

    cout << "nPts0 = " << nPts <<", nBndPts0 = " << nBndPts << endl;


    Array<OneD, Array<OneD, const NekDouble> > coeffs(nCoordDim);
    Array<OneD, Array<OneD, NekDouble> >       pts(nCoordDim);
    for (int i=0; i<nCoordDim; ++i) 
    {
        pts[i]     = Array<OneD, NekDouble>(nPts);
        coeffs[i]  = geom->GetCoeffs(i); // 0/1/2 for x/y/z
        xmap->BwdTrans(coeffs[i], pts[i]);
    }

    Array<OneD, Array<OneD, const NekDouble> > bndCoeffs(nCoordDim);
    Array<OneD, Array<OneD, NekDouble> >       bndPts(nCoordDim);
    for (int i=0; i<nCoordDim; ++i) 
    {
        bndPts[i]     = Array<OneD, NekDouble>(nBndPts);
        bndCoeffs[i]  = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
        bndXmap->BwdTrans(bndCoeffs[i], bndPts[i]);
    }



    for (int i=0;i<nPts;++i){
        cout << pts[0][i] << ", " << pts[1][i] << endl;
    }
    cout << "---"<<endl;
    for (int i=0;i<nBndPts;++i){
        cout << bndPts[0][i] << ", " << bndPts[1][i] << endl;
    }
    */



    /*
    //-------------------------------------------------------------------------
    // Find the element that contains the origin, and give the precise origin
    SpatialDomains::GeometrySharedPtr bndGeom;
    Array<OneD, NekDouble> locCoord(nBndLcoordDim, -999.0);
    NekDouble projDist;
    int elmtid;
    bool isInside = false;

    // Search and get precise locCoord
    for (elmtid=0; elmtid<nBndElmts; ++elmtid) //nBndElmts
    {    
        bndGeom  = BndExp[0]->GetExp(elmtid)->GetGeom();         
        isInside = BndElmtContainsPoint(bndGeom, orig, projDir, locCoord,
                                        projDist, maxDist);
        if (isInside) 
        {
            break;
        }

    }
    ASSERTL0(isInside, "Failed to find the sampling origin on the boundary.");
    */
 
    /*
    // Then Update the precise sampling position
    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    const int npts = bndXmap->GetTotPoints();
    Array<OneD, Array<OneD, NekDouble> >       pts(nCoordDim);
    Array<OneD, Array<OneD, const NekDouble> > bndCoeffs(nCoordDim);

    for (int i=0; i<nCoordDim; ++i) 
    {
        pts[i] = Array<OneD, NekDouble>(npts);
        bndCoeffs[i] = bndGeom->GetCoeffs(i);              // 0/1/2 for x/y/z
        bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
        orig[i] = bndXmap->PhysEvaluate(locCoord, pts[i]); // Update wall point
    }
 
    
    // Get outward-pointing normal vectors for all quadrature points on bnd
    // Use these normals as diretion references
    Array<OneD, Array<OneD, NekDouble> > normalsQ; 
    m_f->m_exp[0]->GetBoundaryNormals(bnd, normalsQ);

    // Get the normals in the element
    // Set point key to get point id
    Array<OneD, LibUtilities::PointsKey> from_key(nBndLcoordDim);
    int from_nPtsPerElmt = 1;
    for (int i=0; i<nBndLcoordDim; ++i)
    {
        from_key[i] = BndExp[0]->GetExp(elmtid)->GetBasis(i)->GetPointsKey();
        from_nPtsPerElmt *= from_key[i].GetNumPoints();
    }

    //Get ref normal
    Array< OneD, NekDouble > normalsRef(3, 0.0);
    int refId = elmtid*from_nPtsPerElmt + 0; // the 1st normal in the bnd elmt
    for(int i=0;i<m_spacedim; ++i)
    {
        normalsRef[i] = normalsQ[i][refId];
    }
    
    // Get the precise normals and correct the direction to be inward
    Array< OneD, NekDouble > normals(3, 0.0);
    GetNormals(bndGeom, locCoord, normals);

    if(Vmath::Dot(3, normals, normalsRef) > 0.0)
    {
        Vmath::Neg(3, normals, 1);
    }
    */
    
    

}


// Compute point-to-point distance to find the cloest element
NekDouble ProcessBodyFittedVelocity::PntToBndElmtPntDistance(
    const Array<OneD, Array<OneD, NekDouble> > & pts,
    const int pId,
    const Array<OneD, Array<OneD, NekDouble> > & bndPts)
{
    const int nCoordDim = pts.size();

    NekDouble dist = 9999, dist_tmp;

    for (int i=0; i<bndPts.size(); ++i)
    {
        // Compute distance
        dist_tmp = 0.0;
        for (int j=0; j<nCoordDim; ++j)
        {
            dist_tmp += (bndPts[j][i]-pts[j][pId])*(bndPts[j][i]-pts[j][pId]);
        }
        dist_tmp = sqrt(dist_tmp);

        if (dist_tmp < dist)
        {
            dist = dist_tmp;
        }
    }

    return dist;
}


bool ProcessBodyFittedVelocity::LocCoordForNearestPntOnBndElmt_2D(
    const Array<OneD, const NekDouble > & gloCoord,
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array<OneD, Array<OneD, NekDouble> > & pts,
    Array<OneD, NekDouble > & locCoord,
    NekDouble & dist,
    const NekDouble iterTol,
    const int iterMax)
{

    // Initial settings
    Array<OneD, NekDouble> etaLR(2);        // range [-1,1]
    etaLR[0] = -1.0;                        // left
    etaLR[1] =  1.0;                        // right
    NekDouble tmpLx, tmpLy, tmpRx, tmpRy;   // tmp values for L/R
    NekDouble distL2, distR2;               // distance square

    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    locCoord[0] = -2.0;
    int cnt = 0;
    bool isConverge = false;
    while(cnt<iterMax)
    {
        tmpLx = bndXmap->PhysEvaluate(etaLR  , pts[0]);
        tmpLy = bndXmap->PhysEvaluate(etaLR  , pts[1]);
        tmpRx = bndXmap->PhysEvaluate(etaLR+1, pts[0]);
        tmpRy = bndXmap->PhysEvaluate(etaLR+1, pts[1]);
        
        distL2 = (tmpLx-gloCoord[0])*(tmpLx-gloCoord[0])
               + (tmpLy-gloCoord[1])*(tmpLy-gloCoord[1]);
        distR2 = (tmpRx-gloCoord[0])*(tmpRx-gloCoord[0])
               + (tmpRy-gloCoord[1])*(tmpRy-gloCoord[1]);
        
        if (distL2 >= distR2)
        {
            etaLR[0]  = 0.5 * (etaLR[0]+etaLR[1]);
        }
        else
        {
            etaLR[1]  = 0.5 * (etaLR[0]+etaLR[1]);
        }


        if ((etaLR[1]-etaLR[0]) < iterTol)
        {
            locCoord[0] = 0.5 * (etaLR[0]+etaLR[1]);
            dist        = sqrt(0.5*(distL2+distR2));
            isConverge  = true;
            break;
        }

        ++cnt;
    }

    // Warning if failed
    if(cnt >= iterMax)
    {
        WARNINGL1(false, "Bisection iteration is not converged");
    }

    return isConverge;
}



// Compute the locCoord for the nearest point on the bnd elmt
// The locCoord is the position to set up the body-fitted coordinate
bool ProcessBodyFittedVelocity::LocCoordForNearestPntOnBndElmt(
    const Array<OneD, const NekDouble > & gloCoord,
    SpatialDomains::GeometrySharedPtr bndGeom,
    Array<OneD, NekDouble > & locCoord,
    Array<OneD, NekDouble > & normals,
    NekDouble & dist,
    const NekDouble iterTol,
    const int iterMax)
{
    bool isConverge = false;
    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0);     // =2 for 2.5D cases

    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    int nBndPts = bndXmap->GetTotPoints();

    Array<OneD, Array<OneD, const NekDouble> > bndCoeffs(nCoordDim);
    Array<OneD, Array<OneD, NekDouble> >       bndPts(nCoordDim);
    for (int i=0; i<nCoordDim; ++i) 
    {
        bndPts[i]    = Array<OneD, NekDouble>(nBndPts);
        bndCoeffs[i] = bndGeom->GetCoeffs(i);
        bndXmap->BwdTrans(bndCoeffs[i], bndPts[i]);
    }


    // Compute the locCoord for the pnt on the bnd elmt
    if (nCoordDim==2)
    {
        // Bisection search
        isConverge = LocCoordForNearestPntOnBndElmt_2D(
                        gloCoord, bndGeom, bndPts, locCoord, dist, iterTol, iterMax);
    }
    else
    {
        // Bi-bisection search
        ASSERTL0(false, "Not available at the moment.");
    }

    // Compute the wall normal
    // Using the function in wallNormalData class
    ProcessWallNormalData wnd(m_f);
    wnd.GetNormals(bndGeom, locCoord, normals);

    return isConverge;
}




//========================================================================


/**
 * @brief Use iteration to get the locCoord. This routine should be used after
 *        we have checked the projected point is inside the projected element.
 * @param bndGeom      Geometry to get the xmap.
 * @param gloCoord     Global coordinate of the point. size=3.
 * @param pts          Global coordinate of the vertices of the elmt. size=2/3.
 * @param dieUse       The main direction(s) used to compute local coordinate
 * @param locCoord     Iteration results for local coordinate(s)
 * @param dist         Returned distance in physical space if the collapsed 
 *                     locCoord is out of range [-1,1].
 * @param iterTol      Tolerence for iteration.
 * @param iterMax      Maximum iteration steps
 * @return             Converged (true) or not (false)
 */
bool ProcessBodyFittedVelocity::BisectionForLocCoordOnBndElmt(
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array<OneD, const NekDouble > & gloCoord,
    const Array<OneD, const Array<OneD, NekDouble> > & pts,
    const Array<OneD, const int > & dirUse,
    Array<OneD, NekDouble > & locCoord,
    const NekDouble iterTol,
    const int iterMax)
{
    // Initial settings
    Array<OneD, NekDouble> etaLR(2); // range [-1,1]
    etaLR[0] = -1.0;                 // left
    etaLR[1] =  1.0;                 // right
    NekDouble tmpL, tmpR;            // tmp values for L/R

    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    locCoord[0] = -2.0;
    int cnt = 0;
    bool isConverge = false;
    while(cnt<iterMax)
    {
        tmpL = bndXmap->PhysEvaluate(etaLR  , pts[dirUse[0]]);
        tmpR = bndXmap->PhysEvaluate(etaLR+1, pts[dirUse[0]]);

        if (fabs(gloCoord[dirUse[0]]-tmpL) >= fabs(gloCoord[dirUse[0]]-tmpR))
        {
            etaLR[0] = 0.5 * (etaLR[0]+etaLR[1]);
        }
        else
        {
            etaLR[1] = 0.5 * (etaLR[0]+etaLR[1]);
        }

        if ((etaLR[1]-etaLR[0]) < iterTol)
        {
            locCoord[0] = 0.5 * (etaLR[0]+etaLR[1]);
            isConverge = true;
            break;
        }

        ++cnt;
    }

    // Warning if failed
    if(cnt >= iterMax)
    {
        WARNINGL1(false, "Bisection iteration is not converged");
    }

    return isConverge;
}


bool ProcessBodyFittedVelocity::NewtonIterForLocCoordOnBndElmt(
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array<OneD, const NekDouble> & gloCoord,
    const Array<OneD, const Array<OneD, NekDouble> > & pts,
    const Array<OneD, const int > & dirUse,
    Array<OneD, NekDouble> & locCoord,
    NekDouble & dist,
    const NekDouble iterTol,
    const int iterMax)
{

    const NekDouble LcoordDiv = 15.0;

    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();

    Array<OneD, const NekDouble> Jac =
        bndGeom->GetMetricInfo()->GetJac(bndXmap->GetPointsKeys());
    NekDouble scaledTol = Vmath::Vsum(Jac.size(), Jac, 1) /
                          ((NekDouble)Jac.size());
    scaledTol *= iterTol;


    // Set the gloCoord used to compute locCoord
    const int dir1 = dirUse[0]; 
    const int dir2 = dirUse[1];

    Array<OneD, NekDouble> Dx1D1(pts[dir1].size());
    Array<OneD, NekDouble> Dx1D2(pts[dir1].size());
    Array<OneD, NekDouble> Dx2D1(pts[dir1].size());
    Array<OneD, NekDouble> Dx2D2(pts[dir1].size());

    // Ideally this will be stored in m_geomfactors
    bndXmap->PhysDeriv(pts[dir1], Dx1D1, Dx1D2);
    bndXmap->PhysDeriv(pts[dir2], Dx2D1, Dx2D2);
       
    // Initial the locCoord, in [-1,1]
    locCoord[0] = 0.0;
    locCoord[1] = 0.0;

    NekDouble x1map, x2map, F1, F2;
    NekDouble derx1_1, derx1_2, derx2_1, derx2_2, jac;
    NekDouble resid;
    int cnt = 0;
    bool isConverge = false;
    

    F1 = F2 = 2000; // Starting value of Function
    while (cnt++ < iterMax)
    {
        x1map = bndXmap->PhysEvaluate(locCoord, pts[dir1]);
        x2map = bndXmap->PhysEvaluate(locCoord, pts[dir2]);

        F1 = gloCoord[dir1] - x1map;
        F2 = gloCoord[dir2] - x2map;

        if (F1 * F1 + F2 * F2 < scaledTol)
        {
            resid = sqrt(F1 * F1 + F2 * F2);
            isConverge = true;
            break;
        }

        // Interpolate derivative metric at locCoord
        derx1_1 = bndXmap->PhysEvaluate(locCoord, Dx1D1);
        derx1_2 = bndXmap->PhysEvaluate(locCoord, Dx1D2);
        derx2_1 = bndXmap->PhysEvaluate(locCoord, Dx2D1);
        derx2_2 = bndXmap->PhysEvaluate(locCoord, Dx2D2);

        jac = derx2_2 * derx1_1 - derx2_1 * derx1_2;
        
        // Use analytical inverse of derivitives which are
        // also similar to those of metric factors.
        locCoord[0] = locCoord[0] + ( derx2_2 * (gloCoord[dir1] - x1map) -
                                      derx1_2 * (gloCoord[dir2] - x2map)) / jac;

        locCoord[1] = locCoord[1] + (-derx2_1 * (gloCoord[dir1] - x1map) +
                                      derx1_1 * (gloCoord[dir2] - x2map)) / jac;

        
        // locCoord have diverged so stop iteration
        if( !(std::isfinite(locCoord[0]) && std::isfinite(locCoord[1])) )
        {
            dist = 1e16;
            std::ostringstream ss;
            ss << "nan or inf found in NewtonIterForLocCoordOnProjBndElmt in element "
               << bndGeom->GetGlobalID();
            WARNINGL1(false, ss.str());
            return false;
        }
        if (fabs(locCoord[0]) > LcoordDiv || fabs(locCoord[1]) > LcoordDiv)
        {
            break; 
        }
    }

    // Check distance for collapsed coordinate 
    Array<OneD, NekDouble> eta(2);
    bndXmap->LocCoordToLocCollapsed(locCoord, eta);

    if(bndGeom->ClampLocCoords(eta, 0.0))
    {
        // calculate the global point corresponding to locCoord
        x1map = bndXmap->PhysEvaluate(eta, pts[dir1]);
        x2map = bndXmap->PhysEvaluate(eta, pts[dir2]);

        F1 = gloCoord[dir1] - x1map;
        F2 = gloCoord[dir2] - x2map;

        dist = sqrt(F1 * F1 + F2 * F2);
    }
    else
    {
        dist = 0.0;
    }

    // Warning if failed
    if (cnt >= iterMax)
    {
        Array<OneD, NekDouble> collCoords(2);
        bndXmap->LocCoordToLocCollapsed(locCoord, collCoords);

        // if coordinate is inside element dump error!
        if ((collCoords[0] >= -1.0 && collCoords[0] <= 1.0) &&
            (collCoords[1] >= -1.0 && collCoords[1] <= 1.0))
        {
            std::ostringstream ss;

            ss << "Reached MaxIterations (" << iterMax
               << ") in Newton iteration ";
            ss << "Init value (" << setprecision(4) << 0 << "," << 0
               << ","
               << ") ";
            ss << "Fin  value (" << locCoord[0] << "," << locCoord[1] << ","
               << ") ";
            ss << "Resid = " << resid << " Tolerance = " << sqrt(scaledTol);

            WARNINGL1(cnt < iterMax, ss.str());
        }
    }

    return isConverge;

}


/**
 * @brief Check if a point can be projected onto an oundary element in a given
 *        direction. If yes, give the local coordinates of the projected point.
 *        we have checked the projected point is inside the projected element.
 * @param bndGeom      Pointer to the geometry of the boundary element.
 * @param gloCoord     Global coordinate of the point. size=3.
 * @param projDir      Projection direction, which is used as the reference
 *                     direction in the 3D routine. size=3, norm=1. 
 * @param locCoord     Iteration results for local coordinates (if inside).
 * @param projDist     Projection distance betweem the point to the wall point.
 * @param maxDist      Disntance to check if the wall point is desired.
 * @param iterTol      Tolerence for iteration.
 * @return             Inside (true) or not (false)
 */
/*
bool ProcessBodyFittedVelocity::BndElmtContainsPoint(
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array<OneD, const NekDouble > & gloCoord,
    const Array<OneD, const NekDouble > & projDir,
    Array< OneD, NekDouble > & locCoord,
    NekDouble & projDist,
    const NekDouble maxDist,
    const NekDouble iterTol)
{
    // Get variables
    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    const int npts      = bndXmap->GetTotPoints();
    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0);     // =2 for 2.5D cases

    Array<OneD, Array<OneD, const NekDouble> > bndCoeffs(nCoordDim);
    Array<OneD, Array<OneD, NekDouble> >       pts(nCoordDim);
    Array<OneD, Array<OneD, NekDouble> >       projPts(nCoordDim);
    Array<OneD, NekDouble >                    projGloCoord(3, 0.0);
    
    for (int i=0; i<nCoordDim; ++i) 
    {
        pts[i]     = Array<OneD, NekDouble>(npts);
        projPts[i] = Array<OneD, NekDouble>(npts);
        bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
        bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
    }

    // Project the point and vertices of the element in the input direction
    ProjectPoint(gloCoord, projDir, 0.0, projGloCoord);
    ProjectVertices(pts, projDir, 0.0, projPts);


    // Set the main direction(s) and the minor direction
    // The gloCoord for minor direction will not be used for locCoord iteration
    // dirUse[0] is the main dir for 2D/2.5D cases, dirUse[1]/[2] is the minor one
    // dirUse[0] and dirUse[1] are the main dir for 3D cases, dirUse[2] is hte minor one
    int dirMaxId = 0; // id to get the dir with largest projDir component
    for (int i=1; i<nCoordDim; ++i)
    {
        if (fabs(projDir[i])>fabs(projDir[dirMaxId]))
        {
            dirMaxId = i;
        }
    }

    Array<OneD, int > dirUse(3, 0); 
    if (nCoordDim==2)
    {
        // 2D or 2.5D cases
        if (dirMaxId==0)
        {
            dirUse[0] = 1;
            dirUse[1] = 0;
            dirUse[2] = 2;
        }
        else 
        {
            dirUse[0] = 0;
            dirUse[1] = 1;
            dirUse[2] = 2;
        }
    }
    else
    {
        // 3D cases
        if (dirMaxId==0)
        {
            dirUse[0] = 1;
            dirUse[1] = 2;
            dirUse[2] = 0;
        }
        else if (dirMaxId==1)
        {
            dirUse[0] = 2;
            dirUse[1] = 0;
            dirUse[2] = 1;
        }
        else
        {
            dirUse[0] = 0;
            dirUse[1] = 1;
            dirUse[2] = 2;
        }

    }


    // Check if the projected point is in the projected elmt
    // If yes, then compute the locCoord and check if desired point is found
    if(nCoordDim==2)
    {
        if (isInProjectedArea2D(projGloCoord, projPts, 1.0e-12))
        {
            bool isConverge, isDesired;
            
            isConverge = BisectionForLocCoordOnBndElmt(bndGeom, projGloCoord,
                             projPts, dirUse, locCoord, iterTol);

            Array<OneD, NekDouble > tmp(2, 0.0);
            tmp[0] = bndXmap->PhysEvaluate(locCoord, pts[0]) - gloCoord[0];
            tmp[1] = bndXmap->PhysEvaluate(locCoord, pts[1]) - gloCoord[1];
            projDist = Vmath::Dot(2, tmp, 1, projDir, 1);  // can be negative
            
            isDesired = (projDist > 0.0) && (projDist < maxDist);

            return isConverge && isDesired;
        }
        else
        {
            return false;
        }
        
    }
    else
    {
        if (isInProjectedArea3D(projGloCoord, projPts, projDir, 1.0e-12, 1.0e-6))
        {
            NekDouble dist; 
            bool isConverge, isDesired;

            isConverge = NewtonIterForLocCoordOnBndElmt(bndGeom, projGloCoord,
                             projPts, dirUse, locCoord, dist, iterTol);
            
            if (dist>iterTol)
            {
                std::ostringstream ss;
                ss << "Collapsed locCoord out of range.\n"
                   << "Newton iteration gives the distance: " << dist;
                WARNINGL1(false, ss.str());
            }
            
            Array<OneD, NekDouble > tmp(3, 0.0);
            tmp[0] = bndXmap->PhysEvaluate(locCoord, pts[0]) - gloCoord[0];
            tmp[1] = bndXmap->PhysEvaluate(locCoord, pts[1]) - gloCoord[1];
            tmp[2] = bndXmap->PhysEvaluate(locCoord, pts[2]) - gloCoord[2];
            projDist = Vmath::Dot(3, tmp, 1, projDir, 1);  // can be negative

            isDesired = (projDist > 0.0) && (projDist < maxDist);

            return isConverge && isDesired;
        }
        else
        {
            return false;
        }
        
    }
    
}
*/




}
}
