////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessVariableProduct.cpp
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
//  Description: Computes product of all variables with others: uu, uv, uw, ...
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <GlobalMapping/Mapping.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessVariableProduct.h"
#include "ProcessMapping.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessVariableProduct::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "varprod"),
    ProcessVariableProduct::create,
    "Computes product of all fields with all other: uu, uv, uw, vv, vw, ww ... ");

ProcessVariableProduct::ProcessVariableProduct(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessVariableProduct::~ProcessVariableProduct()
{
}

void ProcessVariableProduct::Process(po::variables_map &vm)
{
    m_f->SetUpExp(vm);

    int i, j;
    int expdim    = m_f->m_graph->GetMeshDimension();
    int spacedim  = m_f->m_numHomogeneousDir + expdim;
    int nfields   = m_f->m_variables.size();
    int addfields = 0;

    for (i = 0; i < nfields; ++i)
        for (j = i; j < nfields; ++j)
        {
            m_f->m_variables.push_back(m_f->m_variables[i]+m_f->m_variables[j]);
            addfields++;
        }

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    int npoints = m_f->m_exp[0]->GetNpoints();
    m_f->m_exp.resize(nfields + addfields);

    addfields = 0;
    for (i = 0; i < nfields; ++i)
        for (j = i; j < nfields; ++j)
        {
            m_f->m_exp[nfields + addfields] =
                m_f->AppendExpList(m_f->m_numHomogeneousDir);

            Vmath::Vmul(npoints, m_f->m_exp[i]->GetPhys(), 1, 
                            m_f->m_exp[j]->GetPhys(), 1,
                            m_f->m_exp[nfields + addfields]->UpdatePhys(), 1);

            m_f->m_exp[nfields + addfields]->FwdTrans_IterPerExp(
                m_f->m_exp[nfields + addfields]->GetPhys(),
                m_f->m_exp[nfields + addfields]->UpdateCoeffs());
            addfields++;
        }
}
}
}
