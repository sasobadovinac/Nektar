////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessWSS.cpp
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
//  Description: Computes wss field.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include "ProcessWSS.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

using namespace std;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessWSS::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "wss"),
    ProcessWSS::create,
    "Computes wall shear stress field.");

ProcessWSS::ProcessWSS(FieldSharedPtr f) : ProcessBoundaryExtract(f)
{
}

ProcessWSS::~ProcessWSS()
{
}

void ProcessWSS::Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::Process(vm);

    int i, j;
    int nfields = m_f->m_variables.size();
    int expdim  = m_f->m_graph->GetSpaceDimension();
    m_spacedim  = expdim + m_f->m_numHomogeneousDir;

    m_f->m_variables.push_back("Wear");
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    // Resize m_exp
    m_f->m_exp.resize(nfields + 1);
    m_f->m_exp[nfields] = m_f->AppendExpList(m_f->m_numHomogeneousDir);


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

    // Loop over boundaries to Write
    for (int b = 0; b < m_f->m_bndRegionsToWrite.size(); ++b)
    {
        if (BndRegionMap.count(m_f->m_bndRegionsToWrite[b]) == 1)
        {
            int bnd = BndRegionMap[m_f->m_bndRegionsToWrite[b]];
            // Get expansion list for boundary
            MultiRegions::ExpListSharedPtr BndExp =
                m_f->m_exp[nfields + i]->UpdateBndCondExpansion(bnd);

            // Get number of points in expansions
            int nqb = BndExp->GetTotPoints();
            Array<OneD, NekDouble> wear(nqb, 0.0);

            // Calculate wear


            BndExp->FwdTrans_IterPerExp(wear, BndExp->UpdateCoeffs());
        }
    }
}

void ProcessWSS::GetViscosity(NekDouble &kinvis, NekDouble &lambda)
{
    if(boost::iequals(m_f->m_variables[0], "u"))
    {
        // IncNavierStokesSolver
        kinvis = m_f->m_session->GetParameter("Kinvis");
        lambda = 0;
    }
    else if(boost::iequals(m_f->m_variables[0], "rho") &&
            boost::iequals(m_f->m_variables[1], "rhou"))
    {
        // CompressibleFlowSolver
        m_f->m_session->LoadParameter ("mu",     kinvis, 1.78e-05);
        m_f->m_session->LoadParameter ("lambda", lambda, -2.0/3.0);
    }
    else
    {
        // Unknown
        ASSERTL0(false, "Invalid variables for WSS");
    }
}

void ProcessWSS::GetVelocity(
        const Array<OneD, MultiRegions::ExpListSharedPtr> exp,
              Array<OneD, Array<OneD, NekDouble> > &vel)
{
    int npoints = exp[0]->GetNpoints();
    if(boost::iequals(m_f->m_variables[0], "u"))
    {
        // IncNavierStokesSolver
        for (int i = 0; i < m_spacedim; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(npoints);
            Vmath::Vcopy(npoints,
                         exp[i]->GetPhys(), 1,
                         vel[i], 1);
        }
    }
    else if(boost::iequals(m_f->m_variables[0], "rho") &&
            boost::iequals(m_f->m_variables[1], "rhou"))
    {
        // CompressibleFlowSolver
        for (int i = 0; i < m_spacedim; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(npoints);
            Vmath::Vdiv(npoints,
                         exp[i + 1]->GetPhys(), 1,
                         exp[0]->GetPhys(), 1,
                         vel[i], 1);
        }
    }
    else
    {
        // Unknown
        ASSERTL0(false, "Could not identify velocity for WSS");
    }
}

}
}
