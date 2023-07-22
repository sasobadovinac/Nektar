///////////////////////////////////////////////////////////////////////////////
//
// File: AdvectionWeakDG.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Weak DG advection class.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Advection/AdvectionWeakDG.h>
#include <iomanip>
#include <iostream>

#include <LibUtilities/BasicUtils/Timer.h>

namespace Nektar
{
namespace SolverUtils
{
std::string AdvectionWeakDG::type =
    GetAdvectionFactory().RegisterCreatorFunction(
        "WeakDG", AdvectionWeakDG::create, "Weak DG");

AdvectionWeakDG::AdvectionWeakDG()
{
}

/**
 * @brief Initialise AdvectionWeakDG objects and store them before
 * starting the time-stepping.
 *
 * @param pSession  Pointer to session reader.
 * @param pFields   Pointer to fields.
 */
void AdvectionWeakDG::v_InitObject(
    LibUtilities::SessionReaderSharedPtr pSession,
    Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
{
    Advection::v_InitObject(pSession, pFields);
}

/**
 * @brief Compute the advection term at each time-step using the
 * Discontinuous Galerkin approach (DG).
 *
 * @param nConvectiveFields   Number of fields.
 * @param fields              Pointer to fields.
 * @param advVel              Advection velocities.
 * @param inarray             Solution at the previous time-step.
 * @param outarray            Advection term to be passed at the
 *                            time integration class.
 */
void AdvectionWeakDG::v_Advect(
    const int nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &advVel,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble &time,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    LibUtilities::Timer timer1;
    timer1.Start();

    boost::ignore_unused(advVel, time);
    size_t nCoeffs = fields[0]->GetNcoeffs();

    Array<OneD, Array<OneD, NekDouble>> tmp{size_t(nConvectiveFields)};
    for (int i = 0; i < nConvectiveFields; ++i)
    {
        tmp[i] = Array<OneD, NekDouble>{nCoeffs, 0.0};
    }
    LibUtilities::Timer timer;
    timer.Start();
    AdvectCoeffs(nConvectiveFields, fields, advVel, inarray, tmp, time, pFwd,
                 pBwd);
    timer.Stop();
    timer.AccumulateRegion("AdvectCoeffs", 2);

    // Multiply by inverse mass matrix
    LibUtilities::Timer timer2;
    for (int i = 0; i < nConvectiveFields; ++i)
    {
        timer2.Start();
        fields[i]->MultiplyByElmtInvMass(tmp[i], tmp[i]);
        timer2.Stop();
        timer2.AccumulateRegion("AdvWeakDG:_MultiplyByElmtInvMass", 1);
    }

    LibUtilities::Timer timer3;
    timer3.Start();
    for (int i = 0; i < nConvectiveFields; ++i)
    {
        fields[i]->BwdTrans(tmp[i], outarray[i]);
    }
    timer3.Stop();
    timer3.AccumulateRegion("AdvWeakDG:_BwdTrans", 1);
    timer1.Stop();
    timer1.AccumulateRegion("AdvWeakDG:All", 10);
}

/**
 *
 */
void AdvectionWeakDG::v_AdvectCoeffs(
    const int nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &advVel,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble &time,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    LibUtilities::Timer timer1;
    timer1.Start();

    size_t nPointsTot      = fields[0]->GetTotPoints();
    size_t nCoeffs         = fields[0]->GetNcoeffs();
    size_t nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> fluxvector(
        nConvectiveFields);
    // Allocate storage for flux vector F(u).
    for (int i = 0; i < nConvectiveFields; ++i)
    {
        fluxvector[i] = Array<OneD, Array<OneD, NekDouble>>{size_t(m_spaceDim)};
        for (int j = 0; j < m_spaceDim; ++j)
        {
            fluxvector[i][j] = Array<OneD, NekDouble>(nPointsTot, 0.0);
        }
    }

    LibUtilities::Timer timer;
    timer.Start();
    AdvectVolumeFlux(inarray, fluxvector);
    timer.Stop();
    timer.AccumulateRegion("AdvWeakDG:_fluxVector", 3);

    timer.Start();
    // Get the advection part (without numerical flux)
    for (int i = 0; i < nConvectiveFields; ++i)
    {
        Vmath::Fill(outarray[i].size(), 0.0, outarray[i], 1);
        fields[i]->IProductWRTDerivBase(fluxvector[i], outarray[i]);
    }
    timer.Stop();
    timer.AccumulateRegion("AdvWeakDG:_IProductWRTDerivBase", 10);

    Array<OneD, Array<OneD, NekDouble>> numflux{size_t(nConvectiveFields)};

    for (int i = 0; i < nConvectiveFields; ++i)
    {
        numflux[i] = Array<OneD, NekDouble>{nTracePointsTot, 0.0};
    }

    AdvectTraceFlux(nConvectiveFields, fields, advVel, inarray, numflux, time,
                    pFwd, pBwd);

    // Evaulate <\phi, \hat{F}\cdot n> - OutField[i]
    for (int i = 0; i < nConvectiveFields; ++i)
    {
        Vmath::Neg(nCoeffs, outarray[i], 1);

        timer.Start();
        fields[i]->AddTraceIntegral(numflux[i], outarray[i]);
        timer.Stop();
        timer.AccumulateRegion("AdvWeakDG:_AddTraceIntegral", 1);
    }

    timer1.Stop();
    timer1.AccumulateRegion("AdvWeakDG: Coeff All");
}

/**
 *
 */
void AdvectionWeakDG::AdvectTraceFlux(
    const int nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &advVel,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &TraceFlux, const NekDouble &time,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    boost::ignore_unused(advVel, time);
    int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();

    ASSERTL1(m_riemann, "Riemann solver must be provided for AdvectionWeakDG.");

    // Store forwards/backwards space along trace space
    Array<OneD, Array<OneD, NekDouble>> Fwd(nConvectiveFields);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nConvectiveFields);

    if (pFwd == NullNekDoubleArrayOfArray || pBwd == NullNekDoubleArrayOfArray)
    {
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
            Bwd[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
            fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
        }
    }
    else
    {
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            Fwd[i] = pFwd[i];
            Bwd[i] = pBwd[i];
        }
    }

    LibUtilities::Timer timer;
    timer.Start();
    m_riemann->Solve(m_spaceDim, Fwd, Bwd, TraceFlux);
    timer.Stop();
    timer.AccumulateRegion("AdvWeakDG:_Riemann", 10);
}

} // end of namespace SolverUtils
} // end of namespace Nektar
