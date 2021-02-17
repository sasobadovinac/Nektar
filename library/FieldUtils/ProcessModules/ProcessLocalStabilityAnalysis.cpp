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
    m_config["x"]  = ConfigOption(false, "0.1", 
                     "Sampling position given by x-coordinate.");
    m_config["y"]  = ConfigOption(false, "0.0",
                     "Sampling position given by y-coordinate.");
    m_config["z"]  = ConfigOption(false, "0.0",
                     "Sampling position given by x-coordinate.");

    m_config["h"]  = ConfigOption(false, "0.01", 
                     "Sampling distance along the wall normals.");
    m_config["nh"]  = ConfigOption(false, "11", 
                     "Number of sampling points along the wall normals.");

}

ProcessLocalStabilityAnalysis::~ProcessLocalStabilityAnalysis()
{
}

void ProcessLocalStabilityAnalysis::Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::Process(vm);


    // Initialize sampling parameters
    const NekDouble orig_x     = m_config["x"].as<NekDouble>(); //
    const NekDouble orig_y     = m_config["y"].as<NekDouble>(); //
    const NekDouble orig_z     = m_config["z"].as<NekDouble>(); //
    const NekDouble distance_h = m_config["h"].as<NekDouble>(); //
    const int       npts_h     = m_config["nh"].as<int>();      //  
    
    cout << orig_x << ", " << orig_y << ", " << orig_z << ", " << distance_h << ", " << npts_h << endl;

    
    // 
    int nfields   = m_f->m_variables.size();
    int nCoordDim = m_f->m_exp[0]->GetCoordim(0);
    m_spacedim    = nCoordDim + m_f->m_numHomogeneousDir;


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
    int bnd = BndRegionMap[m_f->m_bndRegionsToWrite[0]]; // Ref wnd module


    // Get expansion list for boundary and the number of points
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp(nfields); 
    for (int i = 0; i < nfields; ++i) {
        BndExp[i] = m_f->m_exp[i]->UpdateBndCondExpansion(bnd);
    }
    int nqb = BndExp[0]->GetTotPoints(); // points for all HomModesZ planes



    //-------------------------------------------------------------------------
    cout << "nqb = " << nqb << endl;
    cout << "bndExp[0] size = " << BndExp[0]->GetExpSize() << endl;    
    
    SpatialDomains::GeometrySharedPtr bndGeom = BndExp[0]->GetExp(0)->GetGeom(); 
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

    Array<OneD, const NekDouble> bndCoeffsX = bndGeom->GetCoeffs(0);// 0 for x and 1 for y
    Array<OneD, const NekDouble> bndCoeffsY = bndGeom->GetCoeffs(1); 
   
    bndXmap->BwdTrans(bndCoeffsX, ptsx);
    bndXmap->BwdTrans(bndCoeffsY, ptsy);
 
    cout <<"xy = "<<ptsx[0]<< ", " <<ptsx[1]<< ", " <<ptsy[0]<< ", " <<ptsy[1]<< endl;
    
    Array<OneD, DNekMatSharedPtr> I(1); // 2

    Array<OneD, NekDouble> eta2(1);
    eta2[0] = 0.5;
    I[0] = bndXmap->GetBasis(0)->GetI(eta2); //GetI needs a pointer as input
    
    NekDouble x_tmp;//, y_tmp;
    x_tmp = bndXmap->PhysEvaluate(I, ptsx);
    cout << x_tmp << endl;
 
    /*
    // Get inward-pointing wall-normal vectors for all points on bnd
    Array<OneD, Array<OneD, NekDouble> > normals; 
    m_f->m_exp[0]->GetBoundaryNormals(bnd, normals);
    for (int i = 0; i < m_spacedim; ++i) {
        Vmath::Neg(nqb, normals[i], 1);
    }


    // loop the element on the bnd
    Array<OneD, Array<OneD, NekDouble> > from_ptsInElmt(3);
    for (int i=0; i<3; ++i)
    {
        from_ptsInElmt[i]   = Array<OneD, NekDouble>(from_nPtsPerElmt, 0.0);
    }
    */




    //exp = m_f->m_exp[0];
    //cout << "exp size = " << exp->GetExpSize() << endl;
    //LocalRegions::ExpansionSharedPtr Elmt = exp->GetExp(i);
    //SpatialDomains::GeometrySharedPtr   geom = m_geom;
 
    cout << "exp size = " << m_f->m_exp[0]->GetExpSize() << endl;
    SpatialDomains::GeometrySharedPtr geom = m_f->m_exp[0]->GetExp(0)->GetGeom(); 
    StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

    cout << "numBase = " << xmap->GetNumBases() << endl; 

    //------------------------------------------------------------------------- 

    // CurveSharedPtr
    //vector<PointGeomSharedPtr> pts = m_f->m_graph->GetSegGeom(1)->GetCurve()->m_points;
    //SpatialDomains::MeshGraphSharedPtr graph = m_f->m_graph;
    //SpatialDomains::SegGeomSharedPtr segGeom = graph->GetSegGeom(1);
    //SpatialDomains::CurveSharedPtr curve = segGeom->GetCurve();
    //std::vector<SpatialDomains::PointGeomSharedPtr> pts = curve->m_points; 




    // call lst routine, do nothing for now
    call_hello();
   
    call_lst();   

    

}



}
}
