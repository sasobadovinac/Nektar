////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessVelocityDivergence.cpp
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
//  Description: Computes velocity divergence field.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <GlobalMapping/Mapping.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessMapping.h"
#include "ProcessVelocityDivergence.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessVelocityDivergence::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "divergence"),
        ProcessVelocityDivergence::create,
        "Computes divergence of the velocity field.");

ProcessVelocityDivergence::ProcessVelocityDivergence(FieldSharedPtr f)
    : ProcessModule(f)
{
}

ProcessVelocityDivergence::~ProcessVelocityDivergence()
{
}

void ProcessVelocityDivergence::v_Process(po::variables_map &vm)
{
    m_f->SetUpExp(vm);

    int i, s;
    int expdim = m_f->m_graph->GetMeshDimension();
    m_spacedim = expdim;
    if ((m_f->m_numHomogeneousDir) == 1 || (m_f->m_numHomogeneousDir) == 2)
    {
        m_spacedim = 3;
    }
    int nfields = m_f->m_variables.size();
    ASSERTL0(m_spacedim != 1,
             "Error: Divergence for a 1D problem cannot be computed");

    // Append field names
    m_f->m_variables.push_back("divV");

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }
    int npoints = m_f->m_exp[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble>> grad(m_spacedim * m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> outfield(1);

    int nstrips;

    m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

    for (i = 0; i < m_spacedim * m_spacedim; ++i)
    {
        grad[i] = Array<OneD, NekDouble>(npoints);
    }

    outfield[0] = Array<OneD, NekDouble>(npoints);

    Array<OneD, Array<OneD, NekDouble>> tmp(m_spacedim);
    for (int i = 0; i < m_spacedim; i++)
    {
        tmp[i] = Array<OneD, NekDouble>(npoints);
    }

    vector<MultiRegions::ExpListSharedPtr> Exp(nstrips);

    // Get mapping
    GlobalMapping::MappingSharedPtr mapping = ProcessMapping::GetMapping(m_f);

    for (s = 0; s < nstrips; ++s) // homogeneous strip varient
    {
        // Get velocity and convert to Cartesian system,
        //      if it is still in transformed system
        Array<OneD, Array<OneD, NekDouble>> vel(m_spacedim);
        GetVelocity(vel, s);
        if (m_f->m_fieldMetaDataMap.count("MappingCartesianVel"))
        {
            if (m_f->m_fieldMetaDataMap["MappingCartesianVel"] == "False")
            {
                // Initialize arrays and copy velocity
                if (m_f->m_exp[0]->GetWaveSpace())
                {
                    for (int i = 0; i < m_spacedim; ++i)
                    {
                        m_f->m_exp[0]->HomogeneousBwdTrans(npoints, vel[i],
                                                           vel[i]);
                    }
                }
                // Convert velocity to cartesian system
                mapping->ContravarToCartesian(vel, vel);
                // Convert back to wavespace if necessary
                if (m_f->m_exp[0]->GetWaveSpace())
                {
                    for (int i = 0; i < m_spacedim; ++i)
                    {
                        m_f->m_exp[0]->HomogeneousFwdTrans(npoints, vel[i],
                                                           vel[i]);
                    }
                }
            }
        }

        // Calculate Gradient
        if (m_spacedim == 2)
        {
            for (i = 0; i < m_spacedim; ++i)
            {
                m_f->m_exp[s * nfields + i]->PhysDeriv(vel[i], tmp[0], tmp[1]);
                mapping->CovarToCartesian(tmp, tmp);
                for (int j = 0; j < m_spacedim; j++)
                {
                    Vmath::Vcopy(npoints, tmp[j], 1, grad[i * m_spacedim + j],
                                 1);
                }
            }
            // diV = Ux + Vy
            Vmath::Vadd(npoints, grad[0 * m_spacedim + 0], 1,
                        grad[1 * m_spacedim + 1], 1, outfield[0], 1);
        }
        else
        {
            for (i = 0; i < m_spacedim; ++i)
            {
                m_f->m_exp[s * nfields + i]->PhysDeriv(vel[i], tmp[0], tmp[1],
                                                       tmp[2]);
                mapping->CovarToCartesian(tmp, tmp);
                for (int j = 0; j < m_spacedim; j++)
                {
                    Vmath::Vcopy(npoints, tmp[j], 1, grad[i * m_spacedim + j],
                                 1);
                }
            }

            // diV = Ux + Vy + Wz
            Vmath::Vadd(npoints, grad[0 * m_spacedim + 0], 1,
                        grad[1 * m_spacedim + 1], 1, outfield[0], 1);
            Vmath::Vadd(npoints, outfield[0], 1, grad[2 * m_spacedim + 2], 1,
                        outfield[0], 1);
        }

        Exp[s] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
        Vmath::Vcopy(npoints, outfield[0], 1, Exp[s]->UpdatePhys(), 1);
        Exp[s]->FwdTransLocalElmt(outfield[0], Exp[s]->UpdateCoeffs());
    }

    for (s = 0; s < nstrips; ++s)
    {
        m_f->m_exp.insert(m_f->m_exp.begin() + s * (nfields + 1) + nfields,
                          Exp[s]);
    }
}

void ProcessVelocityDivergence::GetVelocity(
    Array<OneD, Array<OneD, NekDouble>> &vel, int strip)
{
    int nfields = m_f->m_variables.size();
    int npoints = m_f->m_exp[0]->GetNpoints();
    if (boost::iequals(m_f->m_variables[0], "u"))
    {
        // IncNavierStokesSolver
        for (int i = 0; i < m_spacedim; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(npoints);
            Vmath::Vcopy(npoints, m_f->m_exp[strip * nfields + i]->GetPhys(), 1,
                         vel[i], 1);
        }
    }
    else if (boost::iequals(m_f->m_variables[0], "rho") &&
             boost::iequals(m_f->m_variables[1], "rhou"))
    {
        // CompressibleFlowSolver
        for (int i = 0; i < m_spacedim; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(npoints);
            Vmath::Vdiv(npoints, m_f->m_exp[strip * nfields + i + 1]->GetPhys(),
                        1, m_f->m_exp[strip * nfields + 0]->GetPhys(), 1,
                        vel[i], 1);
        }
    }
    else
    {
        // Unknown
        ASSERTL0(false,
                 "Could not identify velocity for ProcessVelocityDivergence");
    }
}

} // namespace FieldUtils
} // namespace Nektar
