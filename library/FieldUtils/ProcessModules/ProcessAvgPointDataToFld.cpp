////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessAvgPointDataToFld.cpp
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
//  Description: Average of point values given data to a fld file
//
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>
#include <boost/geometry.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <FieldUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/BasicUtils/CsvIO.h>

#include "ProcessAvgPointDataToFld.h"

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessAvgPointDataToFld::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "avgpointdatatofld"),
        ProcessAvgPointDataToFld::create,
        "Evaluate discrete data using average to a fld file given a xml file");

ProcessAvgPointDataToFld::ProcessAvgPointDataToFld(FieldSharedPtr f)
    : ProcessModule(f)
{
    m_config["frompts"] = ConfigOption(
        false, "NotSet", "Pts file from which to interpolate field");

    m_config["interpcoord"] =
        ConfigOption(false, "-1", "coordinate id to use for interpolation");
}

ProcessAvgPointDataToFld::~ProcessAvgPointDataToFld()
{
}

void ProcessAvgPointDataToFld::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

//    int i, j;
    LibUtilities::PtsFieldSharedPtr fieldPts;
    // Load pts file
    ASSERTL0( m_config["frompts"].as<string>().compare("NotSet") != 0,
            "ProcessAvgPointDataToFld requires frompts parameter");
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
    ASSERTL0(nFields > 0, "No field values provided in input");

    // Define new expansions.
    ASSERTL0(m_f->m_numHomogeneousDir == 0,
        "ProcessAvgPointDataToFld does not support homogeneous expansion");

    m_f->m_exp.resize(nFields);
    for (int i = 1; i < nFields; ++i)
    {
        m_f->m_exp[i] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
    }

    int totpoints = m_f->m_exp[0]->GetTotPoints();
    Array<OneD, Array<OneD, NekDouble> > intFields(3 + nFields);
    for (int i = 0; i < 3 + nFields; ++i)
    {
        intFields[i] = Array<OneD, NekDouble>(totpoints);
    }
    m_f->m_exp[0]->GetCoords(intFields[0], intFields[1], intFields[2]);
    LibUtilities::PtsFieldSharedPtr outPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(3, intFields);

    int coord_id = m_config["interpcoord"].as<int>();
    ASSERTL0(coord_id <= static_cast<int>(outPts->GetDim()) - 1,
             "Avgcoord is bigger than the Pts files dimension");

    int dim = fieldPts->GetDim();
    int numelemts = m_f->m_exp[0]->GetExpSize();  
    Array<OneD, Array<OneD, NekDouble> > pts; fieldPts->GetPts(pts);
    
   for (int j = 0; j < nFields; ++j)
    {
      // Count the collisions in each elements of each boundary
      Array<OneD, NekDouble> elmtAvg(numelemts,0.0);
      Array<OneD, NekDouble> PtsAvg(totpoints,0.0);
      Array<OneD, int> elmtCount(numelemts,0);
      
     for (int k = 0; k < pts[0].size(); ++k)
      {
         Array<OneD, NekDouble> newCoord(3,0.0), locCoord(3,0.0);
         for (int i = 0; i < dim; ++i)
         {
            newCoord[i] = pts[i][k];
         }
         int eId = m_f->m_exp[0]->GetExpIndex(newCoord,locCoord,
                                          NekConstants::kNekZeroTol);
         if (eId > -1)
         {
            elmtAvg[eId] +=  pts[dim+j][k];
            elmtCount[eId]++;
         }
       }
              
      for (int i = 0; i < numelemts; ++i)
      {
        if (elmtCount[i] != 0)
        {
           LocalRegions::ExpansionSharedPtr elmt = m_f->m_exp[0]->GetExp(i);
           int nElmtPts = elmt->GetTotPoints();
           Array<OneD, NekDouble> avg(nElmtPts, 0.0);

           for (int k = 0; k < nElmtPts; ++k)
           {
              avg[k] = elmtAvg[i] / elmtCount[i];
           }
           Vmath::Vcopy(nElmtPts, &avg[0], 1, &PtsAvg[m_f->m_exp[0]->GetPhys_Offset(i)], 1);
        }
      }
      m_f->m_exp[j]->FwdTrans_IterPerExp(PtsAvg,m_f->m_exp[j]->UpdateCoeffs());
     }

    // save field names
    for (int j = 0; j < nFields; ++j)
    {
        m_f->m_variables.push_back(fieldPts->GetFieldName(j));
    }
}
}
}
