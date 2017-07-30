////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessAddConstToFld.cpp
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

#include "ProcessAddConstToFld.h"

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessAddConstToFld::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "addconsttofld"),
    ProcessAddConstToFld::create,
    "add  a constant value to a field file.  Must specify constant value and "
    "optionally specify which field id to modify. ");

ProcessAddConstToFld::ProcessAddConstToFld(FieldSharedPtr f) : ProcessModule(f)
{
    ASSERTL0(f->m_inputfiles.count("xml"),"Must provide an xml field to use this module");

    m_priority = eModifyExp;

    m_config["value"]   = ConfigOption(false, "NotSet",
                                         "constant to add to fld ");
    m_config["fieldid"]   = ConfigOption(false, "NotSet",
                                         "field id to modify");

}

ProcessAddConstToFld::~ProcessAddConstToFld()
{
}

void ProcessAddConstToFld::Process(po::variables_map &vm)
{
    NekDouble value = 0.0;
    
    if(m_config["value"].m_beenSet)
    {
        value = m_config["value"].as<NekDouble>();
    }
    else
    {
        ASSERTL0(false,"Must specify the constant \"value\" to use the AddConstToFld module");
    }

    int fieldid_start = 0;
    int fieldid_stop = m_f->m_variables.size();
    
    if(m_config["fieldid"].m_beenSet)
    {
        fieldid_start = m_config["fieldid"].as<int>();
        
        ASSERTL1(fieldid_start > 0, "Specified fieldid must be positivie");
        
        ASSERTL1(fieldid_start < fieldid_stop,"Specified fieldid is larger than the number of field files in input fld file");
        
        fieldid_stop = fieldid_start+1;
    }

    
    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }
    
    int nphys = m_f->m_exp[0]->GetTotPoints();

    for (int j = fieldid_start; j < fieldid_stop; ++j)
    {
        m_f->m_exp[j]->BwdTrans(m_f->m_exp[j]->GetCoeffs(),
                                m_f->m_exp[j]->UpdatePhys());

        // add constant
        Vmath::Sadd(nphys, value, m_f->m_exp[j]->GetPhys(),1,
                    m_f->m_exp[j]->UpdatePhys(),1);

        m_f->m_exp[j]->FwdTrans_IterPerExp(m_f->m_exp[j]->GetPhys(),
                                           m_f->m_exp[j]->UpdateCoeffs());
    }
}
}
}
