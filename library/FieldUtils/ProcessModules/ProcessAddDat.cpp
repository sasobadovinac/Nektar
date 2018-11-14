////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessAddDat.cpp
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
//  Description: Add a field to the intput field
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include "ProcessAddDat.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessAddDat::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "adddat"),
    ProcessAddDat::create,
    "add two dat fields together with optional scaling. Must specify fromfld and "
    "scaling is optionally specified with input option scale.");

ProcessAddDat::ProcessAddDat(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["scale"] = ConfigOption(false, "1.0", "scale factor");

    m_config["fromdat"] =
        ConfigOption(false, "NotSet", "Dat file form which to add dat");

    m_priority = eModifyPts;
}

ProcessAddDat::~ProcessAddDat()
{
}

void ProcessAddDat::Process(po::variables_map &vm)
{
    string scalestr = m_config["scale"].as<string>();
    NekDouble scale = boost::lexical_cast<NekDouble>(scalestr);

    ASSERTL0(m_config["fromdat"].as<string>().compare("NotSet") != 0,
             "Need to specify fromdat=file.dat ");
    string fromdat = m_config["fromdat"].as<string>();
    
    LibUtilities::PtsFieldSharedPtr fromFieldPts;
    
    InputDatFile(fromdat, fromFieldPts);

    ASSERTL0(fromFieldPts->GetNpoints() == m_f->m_fieldPts->GetNpoints(),
             "dat files are of different length");

    int dim = m_f->m_fieldPts->GetDim();
    for( int n = 0; n < m_f->m_fieldPts->GetNFields(); ++n)
    {
        for (int i = 0; i < m_f->m_fieldPts->GetNpoints(); ++i)
        {
            m_f->m_fieldPts->SetPointVal(dim+n,i,
                                         m_f->m_fieldPts->GetPointVal(dim+n,i)
                                         + scale*fromFieldPts->GetPointVal(dim+n,i));
        }
    }
}
}
}
