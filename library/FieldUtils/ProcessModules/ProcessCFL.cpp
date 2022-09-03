////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessCFL.cpp
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
//  Description: Computes CFL number over the entire domain for the
//  incompressible flow simulaiton. This is helpful in terms of debugging
//  and tracing the evolution of CFL in time over the domain.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <GlobalMapping/Mapping.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessCFL.h"
#include "ProcessMapping.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessCFL::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "CFL"), ProcessCFL::create,
    "Computes CFL number for the entire domain for Incompressible flow.");

ProcessCFL::ProcessCFL(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessCFL::~ProcessCFL()
{
}

void ProcessCFL::Process(po::variables_map &vm)
{
    m_f->SetUpExp(vm);

    int expdim    = m_f->m_graph->GetMeshDimension();
    int nelmt     = m_f->m_exp[0]->GetExpSize();
    int nfields   = m_f->m_variables.size();
    m_spacedim    = expdim;

    NekDouble timeStep = m_f->m_session->GetParameter("TimeStep");
    NekDouble cLambda  = 0.2; // Spencer's book

    if (m_f->m_numHomogeneousDir == 1)
    {
        m_spacedim = 3;
    }
    ASSERTL0(m_f->m_numHomogeneousDir != 2,
             "CFL for 3DH2D simulations is not supported");
    ASSERTL0(m_spacedim != 1, "Error: CFL for a 1D problem is not supported");

    // Append field names
    m_f->m_variables.push_back("CFL");

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }
    int npoints = m_f->m_exp[0]->GetNpoints();
    Array<OneD, NekDouble> outfield(npoints);

    int nstrips;
    m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

    MultiRegions::ExpListSharedPtr Exp;
#if EXPLISTDATA
    //add in new fields 
    for (int s = 0; s < nstrips; ++s)
    {
        Exp = m_f->AppendExpList(m_f->m_numHomogeneousDir);
        m_f->m_exp.insert(m_f->m_exp.begin() + s*(nfields + 1) + nfields, Exp);
    }
#else
    std::vector<MultiRegions::ExpListSharedPtr> varExp;
    // add in new fields 
    for (int s = 0; s < nstrips; ++s)
    {
        Exp = m_f->AppendExpList(m_f->m_numHomogeneousDir);
            
        m_f->m_exp.insert(m_f->m_exp.begin() + s*(nfields + 1)+ nfields, Exp);
            varExp.push_back(Exp);
    }

    m_f->m_fieldPhys  ->AddVariable(varExp);
    m_f->m_fieldCoeffs->AddVariable(varExp);

    // Reshuffle data for strip case before filling new variables.
    for (int s = nstrips-1; s > 0; --s) // homogeneous strip varient
    {    
        int ncoeffs = m_f->m_exp[0]->GetNcoeffs();
        for(int n = nfields; n > 0; --n)
        {
            int fid  = s*(nfields)+n;
            int fid1 = s*(nfields+1)+n;

            Vmath::Vcopy(ncoeffs,m_f->m_fieldCoeffs->GetArray1D(fid-1),1,
                         m_f->m_fieldCoeffs->UpdateArray1D(fid1-1),1);

            Vmath::Vcopy(npoints,m_f->m_fieldPhys->GetArray1D(fid-1),1,
                         m_f->m_fieldPhys->UpdateArray1D(fid1-1),1);
        }
    }
#endif
    
    
    for (int s = 0; s < nstrips; ++s) // homogeneous strip varient
    {
        Array<OneD, Array<OneD, NekDouble>> velocityField(expdim);

        // Get the velocity field
        GetVelocity(velocityField, s);

        // compute the max velocity in the std regions
        Array<OneD, NekDouble> stdVel = GetMaxStdVelocity(velocityField);

        // get the maximum expansion order in each element
        Array<OneD, int> expOrder =
            m_f->m_exp[s * nfields + 0]->EvalBasisNumModesMaxPerExp();

        // compute the CFL number
        Array<OneD, NekDouble> cfl(nelmt);
        for (int el = 0; el < nelmt; ++el)
        {
            int order = std::max(expOrder[el] - 1, 1);
            cfl[el]   = timeStep * stdVel[el] * cLambda * order * order;
        }

        int cnt = 0;
        for (int el = 0; el < nelmt; ++el)
        {
            // using the field[0]==m_exp[s*nfields + 0]
            int nquad = m_f->m_exp[s * nfields + 0]->GetExp(el)->GetTotPoints();
            Vmath::Fill(nquad, cfl[el], &outfield[cnt], 1);
            cnt += nquad;
        }

        // temporary store the CFL number field for each strip
#if EXPLISTDATA
        Vmath::Vcopy(npoints, outfield, 1, m_f->m_exp[s*(nfields+1) + nfields]->
                     UpdatePhys(), 1);
        m_f->m_exp[0]->FwdTransLocalElmt(outfield,
                       m_f->m_exp[s*(nfields+1) + nfields]->UpdateCoeffs());
#else
        Vmath::Vcopy(npoints, outfield, 1, m_f->m_fieldPhys->UpdateArray1D
                     (s*(nfields+1) + nfields), 1); 
        
        m_f->m_exp[0]->FwdTransLocalElmt(outfield,
                          m_f->m_fieldCoeffs->UpdateArray1D(s*(nfields+1)
                                                            + nfields));
#endif
    }
}

void ProcessCFL::GetVelocity(Array<OneD, Array<OneD, NekDouble>> &vel,
                             int strip)
{
    int expdim  = m_f->m_graph->GetMeshDimension();
    int nfields = m_f->m_variables.size();
    int npoints = m_f->m_exp[0]->GetNpoints();
    if (boost::iequals(m_f->m_variables[0], "u"))
    {
        // IncNavierStokesSolver
        // Using expdim instead of spacedim
        // This is because for 3DH1D, only a 2D plane will be considered
        for (int i = 0; i < expdim; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(npoints);
#if EXPLISTDATA
            Vmath::Vcopy(npoints, m_f->m_exp[strip * nfields + i]->GetPhys(), 1,
                         vel[i], 1);
#else
            Vmath::Vcopy(npoints, m_f->m_fieldPhys->GetArray1D
                         (strip * nfields + i), 1, vel[i], 1);
#endif
        }
    }
    else if (boost::iequals(m_f->m_variables[0], "rho") &&
             boost::iequals(m_f->m_variables[1], "rhou"))
    {
        // CompressibleFlowSolver
        ASSERTL0(false, "CFL calculation is not supported for the compressible "
                        "flow simulations at the moment");
    }
    else
    {
        // Unknown
        ASSERTL0(false, "Could not identify velocity for ProcessCFL");
    }
}

/**
 *
 */
Array<OneD, NekDouble> ProcessCFL::GetMaxStdVelocity(
    const Array<OneD, Array<OneD, NekDouble>> &vel, int strip)
{
    int nfields    = m_f->m_variables.size();
    int n_points_0 = m_f->m_exp[0]->GetExp(0)->GetTotPoints();
    int n_element  = m_f->m_exp[0]->GetExpSize();
    int nvel       = vel.size();
    int cnt;

    NekDouble pntVelocity;

    // Getting the standard velocity vector
    Array<OneD, Array<OneD, NekDouble>> stdVelocity(nvel);
    Array<OneD, NekDouble> tmp;
    Array<OneD, NekDouble> maxV(n_element, 0.0);
    LibUtilities::PointsKeyVector ptsKeys;

    for (int i = 0; i < nvel; ++i)
    {
        stdVelocity[i] = Array<OneD, NekDouble>(n_points_0);
    }

    cnt = 0.0;
    for (int el = 0; el < n_element; ++el)
    {
        int n_points = m_f->m_exp[0]->GetExp(el)->GetTotPoints();
        ptsKeys      = m_f->m_exp[0]->GetExp(el)->GetPointsKeys();

        // reset local space
        if (n_points != n_points_0)
        {
            for (int j = 0; j < nvel; ++j)
            {
                stdVelocity[j] = Array<OneD, NekDouble>(n_points, 0.0);
            }
            n_points_0 = n_points;
        }
        else
        {
            for (int j = 0; j < nvel; ++j)
            {
                Vmath::Zero(n_points, stdVelocity[j], 1);
            }
        }

        Array<TwoD, const NekDouble> gmat = m_f->m_exp[strip * nfields + 0]
                                                ->GetExp(el)
                                                ->GetGeom()
                                                ->GetMetricInfo()
                                                ->GetDerivFactors(ptsKeys);

        if (m_f->m_exp[strip * nfields + 0]
                ->GetExp(el)
                ->GetGeom()
                ->GetMetricInfo()
                ->GetGtype() == SpatialDomains::eDeformed)
        {
            for (int j = 0; j < nvel; ++j)
            {
                for (int k = 0; k < nvel; ++k)
                {
                    Vmath::Vvtvp(n_points, gmat[k * nvel + j], 1,
                                 tmp = vel[k] + cnt, 1, stdVelocity[j], 1,
                                 stdVelocity[j], 1);
                }
            }
        }
        else
        {
            for (int j = 0; j < nvel; ++j)
            {
                for (int k = 0; k < nvel; ++k)
                {
                    Vmath::Svtvp(n_points, gmat[k * nvel + j][0],
                                 tmp = vel[k] + cnt, 1, stdVelocity[j], 1,
                                 stdVelocity[j], 1);
                }
            }
        }
        cnt += n_points;

        // Calculate total velocity in stdVelocity[0]
        Vmath::Vmul(n_points, stdVelocity[0], 1, stdVelocity[0], 1,
                    stdVelocity[0], 1);
        for (int k = 1; k < nvel; ++k)
        {
            Vmath::Vvtvp(n_points, stdVelocity[k], 1, stdVelocity[k], 1,
                         stdVelocity[0], 1, stdVelocity[0], 1);
        }
        pntVelocity = Vmath::Vmax(n_points, stdVelocity[0], 1);
        maxV[el]    = sqrt(pntVelocity);
    }

    return maxV;
}
} // namespace FieldUtils
} // namespace Nektar
