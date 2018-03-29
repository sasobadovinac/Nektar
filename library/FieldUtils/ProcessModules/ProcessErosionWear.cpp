////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessErosionWear.cpp
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
//  Description: Computes erosion wear field.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

using namespace std;

#include <boost/geometry.hpp>
#include "ProcessInterpPointDataToFld.h"

#include <FieldUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/BasicUtils/CsvIO.h>
#include <boost/math/special_functions/fpclassify.hpp>

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessInterpPointDataToFld::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "interppointdatatofld"),
        ProcessInterpPointDataToFld::create,
        "Interpolates given discrete data using a finite difference "
        "approximation to a fld file given a xml file");

ProcessInterpPointDataToFld::ProcessInterpPointDataToFld(FieldSharedPtr f)
    : ProcessModule(f)
{
    m_config["frompts"] = ConfigOption(
        false, "NotSet", "Pts file from which to interpolate field");

    m_config["interpcoord"] =
        ConfigOption(false, "-1", "coordinate id to use for interpolation");
}

ProcessInterpPointDataToFld::~ProcessInterpPointDataToFld()
{
}

void ProcessInterpPointDataToFld::Process(po::variables_map &vm)
{
    int i, j;
    LibUtilities::PtsFieldSharedPtr fieldPts;
    // Load pts file
    ASSERTL0( m_config["frompts"].as<string>().compare("NotSet") != 0,
            "ProcessInterpPointDataToFld requires frompts parameter");
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
        "ProcessInterpPointDataToFld does not support homogeneous expansion");

    m_f->m_exp.resize(nFields);
    
