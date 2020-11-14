////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessWear.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
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
//  Description: Computes Erosion Wear field.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include <boost/geometry.hpp>
#include "ProcessWear.h"

#include <FieldUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/BasicUtils/CsvIO.h>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessWear::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "wear"),
        ProcessWear::create,
       "Computes Erosion Wear field.");

ProcessWear::ProcessWear(FieldSharedPtr f)
    : ProcessBoundaryExtract(f)
{
        m_config["frompts"] = ConfigOption(
        false, "NotSet", "Pts file from which to interpolate field");

    m_config["interpcoord"] =
        ConfigOption(false, "-1", "coordinate id to use for interpolation");
}

ProcessWear::~ProcessWear()
{
}

void ProcessWear::Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::Process(vm);

    LibUtilities::PtsFieldSharedPtr fieldPts;
    // Load pts file
    ASSERTL0( m_config["frompts"].as<string>().compare("NotSet") != 0,
            "ProcessWear requires frompts parameter");
    string inFile = m_config["frompts"].as<string>().c_str();

    if (boost::filesystem::path(inFile).extension() == ".pts")
    {
        LibUtilities::PtsIOSharedPtr ptsIO =
            MemoryManager<LibUtilities::PtsIO>::AllocateSharedPtr(m_f->m_comm);

        ptsIO->Import(inFile, fieldPts);
    }
    else if (boost::filesystem::path(inFile).extension() == ".csv")
    {
        LibUtilities::CsvIOSharedPtr csvIO =
            MemoryManager<LibUtilities::CsvIO>::AllocateSharedPtr(m_f->m_comm);

        csvIO->Import(inFile, fieldPts);
    }
    else
    {
        ASSERTL0(false, "unknown frompts file type");
    }

    int nFields = fieldPts->GetNFields();
    ASSERTL0( nFields == 3, 
        "Velocity, collision angle and diameter are required  input fields");
       
    // Define new expansions.
    ASSERTL0(m_f->m_numHomogeneousDir == 0,
        "ProcessWear does not support homogeneous expansion");

    m_f->m_exp.resize(1);
    m_f->m_variables.push_back("Wear");
    //m_f->m_exp[1] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
 
    // Create map of boundary ids for partitioned domains
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

    int dim = fieldPts->GetDim();
    Array<OneD, Array<OneD, NekDouble> > pts; fieldPts->GetPts(pts);
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp(1);
    NekDouble ONE_OVER_SQRT_2PI = 0.39894228040143267793994605993438;
    NekDouble dist2 = 0.0, Sigma = 0.39894228040143267793994605993438, Wear = 0.0;

#if 1  // Gauss counting
    for (int b = 0; b < m_f->m_bndRegionsToWrite.size(); ++b)
    {
        if (BndRegionMap.count(m_f->m_bndRegionsToWrite[b]) == 1)
        {
           int bnd = BndRegionMap[m_f->m_bndRegionsToWrite[b]];
           // Get expansion list for boundary and for elements containing this bnd
           BndExp[0] = m_f->m_exp[0]->UpdateBndCondExpansion(bnd);

           //Create variables for the data generation
           int totpoints = BndExp[0]->GetTotPoints();
           Array<OneD, Array<OneD, NekDouble> > intFields(4);
           for (int i = 0; i < 4 ; ++i)
           {
              intFields[i] = Array<OneD, NekDouble>(totpoints,0.0);
           }
           //Load the quadrature point information
           BndExp[0]->GetCoords(intFields[0], intFields[1], intFields[2]);

           for (int  k = 0; k < pts[0].size() ; ++k)
           {
              /* Wear  = Model1wear(pts[dim][k],fabs(pts[dim+1][k])); */
              /* Wear  = ECRCwear(pts[dim][k],(pts[dim+1][k])); */
              Wear  = TulsaAnsys(pts[dim][k],(pts[dim+1][k]));

              // Sum Gauss function
              for ( int i = 0; i < totpoints; ++i)
              {
                 dist2 = 0.0;
                 for (int j = 0; j < dim; ++j)
                 {
                    dist2 += pow(intFields[j][i] - pts[j][k],2);
                 }

                 intFields[3][i] += Wear * ONE_OVER_SQRT_2PI / Sigma
                    * exp(- 0.5 * dist2 / pow(Sigma,2));
              
              }
           }    

           //Interpolate data
           for (int i = 0; i < totpoints; ++i)
           {
           			BndExp[0]->SetPhys(i, intFields[3][i]);
           }

           // forward transform fields
           BndExp[0]->FwdTrans_IterPerExp(BndExp[0]->GetPhys(),
           										   BndExp[0]->UpdateCoeffs());

        }
    }
#endif

#if 0 // elemental counting
    for (int b = 0; b < m_f->m_bndRegionsToWrite.size(); ++b)
    {
        if (BndRegionMap.count(m_f->m_bndRegionsToWrite[b]) == 1)
        {
           int bnd = BndRegionMap[m_f->m_bndRegionsToWrite[b]];
           // Get expansion list for boundary and for elements containing this bnd
           BndExp[0] = m_f->m_exp[0]->UpdateBndCondExpansion(bnd);
           
           // Count the collisions in each elements of each boundary
           Array<OneD, int> elmtWear(BndExp[0]->GetExpSize(),0.0);
           for (int k = 0; k < pts[0].num_elements(); ++k)
           {
              Array<OneD, NekDouble> newCoord(3,0.0), locCoord(3,0.0);
              for (int i = 0; i < dim; ++i)
              {
                 newCoord[i] = pts[i][k];

              }
              int eId = BndExp[0]->GetExpIndex(newCoord,locCoord,
                                       NekConstants::kNekZeroTol);
              if (eId > -1)
              {
                 elmtWear[eId] +=  ECRCwear(pts[dim][k],(pts[dim+1][k]));
                 /* elmtWear[eId] += Model1wear(pts[dim][k],fabs(pts[dim+1][k])); */
              }
           }

           const int nBndPts = BndExp[0]->GetNpoints();
           Array<OneD, NekDouble> wearBnd(nBndPts, 0.0);

           //Set erosion to each element 
           for (int i = 0; i < BndExp[0]->GetExpSize(); ++i)
           {
              LocalRegions::ExpansionSharedPtr elmt = BndExp[0]->GetExp(i);

              int nElmtPts = elmt->GetTotPoints();
              Array<OneD, NekDouble> wear(nElmtPts, 0.0);

              // Calculate elemental area or length 
              Array<OneD, NekDouble> ones(nElmtPts, 1.0);
              NekDouble scale = elmt->Integral(ones);

              for (int j = 0; j < nElmtPts; ++j)
              {
                 wear[j] = elmtWear[i]/scale;
              }

              Vmath::Vcopy(nElmtPts, &wear[0], 1, &wearBnd[BndExp[0]->GetPhys_Offset(i)], 1);
           }

           BndExp[0]->FwdTrans_IterPerExp(wearBnd, BndExp[0]->UpdateCoeffs());
        }
    }
#endif
    
}

NekDouble ProcessWear::ECRCwear(NekDouble Vel, NekDouble angle)
{
    // Angle function variables
    NekDouble n1 = 0.15;
    NekDouble n2 = 0.85;
    NekDouble n3 = 0.65;
    NekDouble Hv = 1.83;
    NekDouble f  = 1.53;
    NekDouble K  = 2.16E-8;
//    NekDouble C  = 4.62E-7;
    //NekDouble BH = 178.9;
    NekDouble fs = 1.0;
    
    return    K * fs * pow(Vel,2.41)
            * (1/f * pow(sin(angle),n1)) 
            * pow(1 + pow(Hv,n3) * (1-sin(angle)),n2); 
}

NekDouble ProcessWear::Model1wear(NekDouble Vel, NekDouble angle)
{
    //Angle function variables
    NekDouble A = -0.396;
    NekDouble B = 8.380;
    NekDouble C = -16.92;
    NekDouble D = 10.747;
    NekDouble E = -1.765;
    NekDouble F = 0.434; 
    
    return fabs(pow(Vel,2)*(A*pow(sin(angle),4)+B*pow(sin(angle),3)+
                    C*pow(sin(angle),2)+D*sin(angle)+E)*F);
    
}

NekDouble ProcessWear::TulsaAnsys(NekDouble Vel, NekDouble angle)
{
    // Angle function variables
   // NekDouble mp = 1; // mass flow
    NekDouble Fs = 1; 
    NekDouble B = 200; // Brinel hardness 
    NekDouble Fa = 0.0; 
    if (angle<=0.267)
    {
      Fa = 22.7*angle-38.4*pow(angle,2); 
    }
    else
    {
      Fa = 2+6.8*angle-7.5*pow(angle,2)+2.25*pow(angle,3); 
    }
    return 1.559E-6*pow(B,-0.59)*Fs*pow(Vel,1.73)*Fa; 
}

}
}
