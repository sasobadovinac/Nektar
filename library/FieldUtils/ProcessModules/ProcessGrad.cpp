////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessGrad.cpp
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
//  Description: Computes gradient of fields.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <GlobalMapping/Mapping.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessGrad.h"
#include "ProcessMapping.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessGrad::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "gradient"), ProcessGrad::create,
    "Computes gradient of fields.");

ProcessGrad::ProcessGrad(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["vars"] = ConfigOption(false, "NotSet", "Select variables");
    m_config["dirs"] = ConfigOption(false, "NotSet", "Select directions");
}

ProcessGrad::~ProcessGrad()
{
}

void ProcessGrad::ParserOptions(std::set<int> &variables,
                                std::set<int> &directions)
{
    int expdim   = m_f->m_graph->GetMeshDimension();
    int spacedim = m_f->m_numHomogeneousDir + expdim;
    variables.clear();
    if (m_config["vars"].as<string>().compare("NotSet"))
    {
        vector<string> vars;
        ASSERTL0(
            ParseUtils::GenerateVector(m_config["vars"].as<string>(), vars),
            "Failed to interpret variable numbers or names");
        for (string s : vars)
        {
            int v = -1;
            try
            {
                v = boost::lexical_cast<int>(s);
            }
            catch (const boost::bad_lexical_cast &)
            {
                auto index =
                    find(m_f->m_variables.begin(), m_f->m_variables.end(), s);
                v = index - m_f->m_variables.begin();
            }
            if (v < 0 || v >= m_f->m_variables.size())
            {
                WARNINGL0(false, "Warning: variable " + s + " not found");
            }
            else
            {
                variables.insert(v);
            }
        }
    }
    else
    {
        for (int v = 0; v < m_f->m_variables.size(); ++v)
        {
            variables.insert(v);
        }
    }
    directions.clear();
    vector<string> coords = {"x", "y", "z"};
    if (m_config["dirs"].as<string>().compare("NotSet"))
    {
        vector<string> vars;
        ASSERTL0(
            ParseUtils::GenerateVector(m_config["dirs"].as<string>(), vars),
            "Failed to interpret coordinate numbers or names");
        for (string s : vars)
        {
            int d = -1;
            try
            {
                d = boost::lexical_cast<int>(s);
            }
            catch (const boost::bad_lexical_cast &)
            {
                auto index =
                    find(coords.begin(), coords.end(), boost::to_lower_copy(s));
                d = index - coords.begin();
            }
            if (d < 0 || d >= spacedim)
            {
                WARNINGL0(false,
                          "Warning: space direction " + s + " not found");
            }
            else
            {
                directions.insert(d);
            }
        }
    }
    else
    {
        for (int d = 0; d < spacedim; ++d)
        {
            directions.insert(d);
        }
    }
}

void ProcessGrad::ProcessCartesianFld(Array<OneD, Array<OneD, NekDouble>> &grad)
{
    // Calculate Gradient
    int n = 0;
    for (int i : m_selectedVars)
    {
        for (int j : m_directions)
        {
            m_f->m_exp[i]->PhysDeriv(MultiRegions::DirCartesianMap[j],
                                     m_f->m_exp[i]->GetPhys(), grad[n]);
            ++n;
        }
    }
}

void ProcessGrad::ProcessMappingFld(Array<OneD, Array<OneD, NekDouble>> &grad)
{
    int expdim   = m_f->m_graph->GetMeshDimension();
    int spacedim = m_f->m_numHomogeneousDir + expdim;
    int nfields  = m_f->m_variables.size();
    bool hasvel = m_selectedVars.size() && (*m_selectedVars.begin() < spacedim);

    int npoints = m_f->m_exp[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> tmp(spacedim);
    for (int i = 0; i < spacedim; i++)
    {
        tmp[i] = Array<OneD, NekDouble>(npoints);
    }

    // Get mapping
    GlobalMapping::MappingSharedPtr mapping = ProcessMapping::GetMapping(m_f);

    // Get velocity and convert to Cartesian system,
    //      if it is still in transformed system
    Array<OneD, Array<OneD, NekDouble>> vel(spacedim);
    if (hasvel && m_f->m_fieldMetaDataMap.count("MappingCartesianVel"))
    {
        if (m_f->m_fieldMetaDataMap["MappingCartesianVel"] == "False")
        {
            // Initialize arrays and copy velocity
            for (int i = 0; i < spacedim; ++i)
            {
                vel[i] = Array<OneD, NekDouble>(npoints);
                if (m_f->m_exp[0]->GetWaveSpace())
                {
                    m_f->m_exp[0]->HomogeneousBwdTrans(m_f->m_exp[i]->GetPhys(),
                                                       vel[i]);
                }
                else
                {
                    Vmath::Vcopy(npoints, m_f->m_exp[i]->GetPhys(), 1, vel[i],
                                 1);
                }
            }
            // Convert velocity to cartesian system
            mapping->ContravarToCartesian(vel, vel);
            // Convert back to wavespace if necessary
            if (m_f->m_exp[0]->GetWaveSpace())
            {
                for (int i = 0; i < spacedim; ++i)
                {
                    m_f->m_exp[0]->HomogeneousFwdTrans(vel[i], vel[i]);
                }
            }
        }
        else
        {
            for (int i = 0; i < spacedim; ++i)
            {
                vel[i] = Array<OneD, NekDouble>(npoints);
                Vmath::Vcopy(npoints, m_f->m_exp[i]->GetPhys(), 1, vel[i], 1);
            }
        }
    }
    else if (hasvel)
    {
        for (int i = 0; i < spacedim && i < nfields; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(npoints);
            Vmath::Vcopy(npoints, m_f->m_exp[i]->GetPhys(), 1, vel[i], 1);
        }
    }

    // Calculate Gradient
    int n = 0;
    for (int i : m_selectedVars)
    {
        for (int j = 0; j < spacedim; ++j)
        {
            if (i < spacedim)
            {
                m_f->m_exp[i]->PhysDeriv(MultiRegions::DirCartesianMap[j],
                                         vel[i], tmp[j]);
            }
            else
            {
                m_f->m_exp[i]->PhysDeriv(MultiRegions::DirCartesianMap[j],
                                         m_f->m_exp[i]->GetPhys(), tmp[j]);
            }
        }
        mapping->CovarToCartesian(tmp, tmp);
        for (int j : m_directions)
        {
            Vmath::Vcopy(npoints, tmp[j], 1, grad[n], 1);
            ++n;
        }
    }
}

void ProcessGrad::Process(po::variables_map &vm)
{
    m_f->SetUpExp(vm);
    ParserOptions(m_selectedVars, m_directions);

    int nfields   = m_f->m_variables.size();
    int addfields = m_selectedVars.size() * m_directions.size();
    m_f->m_exp.resize(nfields + addfields);

    vector<string> coords = {"x", "y", "z"};
    for (int i : m_selectedVars)
    {
        for (int j : m_directions)
        {
            m_f->m_variables.push_back(m_f->m_variables[i] + "_" + coords[j]);
        }
    }

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    int npoints = m_f->m_exp[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble>> grad(addfields);
    for (int i = 0; i < addfields; ++i)
    {
        grad[i] = Array<OneD, NekDouble>(npoints);
    }

    if (m_f->m_fieldMetaDataMap.count("MappingType"))
    {
        ProcessMappingFld(grad);
    }
    else
    {
        ProcessCartesianFld(grad);
    }

    for (int i = 0; i < addfields; ++i)
    {
        m_f->m_exp[nfields + i] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
        Vmath::Vcopy(npoints, grad[i], 1, m_f->m_exp[nfields + i]->UpdatePhys(),
                     1);
        m_f->m_exp[nfields + i]->FwdTransLocalElmt(
            grad[i], m_f->m_exp[nfields + i]->UpdateCoeffs());
    }
}
} // namespace FieldUtils
} // namespace Nektar
