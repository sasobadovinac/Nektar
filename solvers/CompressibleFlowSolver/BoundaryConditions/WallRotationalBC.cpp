///////////////////////////////////////////////////////////////////////////////
//
// File: WallRotationalBC.cpp
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
// Description: No-slip wall boundary condition
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include "WallRotationalBC.h"

using namespace std;

namespace Nektar
{

std::string WallRotationalBC::classNameRotational = GetCFSBndCondFactory().
    RegisterCreatorFunction("WallRotational",
                            WallRotationalBC::create,
                            "Adiabatic rotational wall boundary condition.");

WallRotationalBC::WallRotationalBC(const LibUtilities::SessionReaderSharedPtr& pSession,
           const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
           const Array<OneD, Array<OneD, NekDouble> >& pTraceNormals,
           const Array<OneD, Array<OneD, NekDouble> >& pGridVelocity,
           const int pSpaceDim,
           const int bcRegion,
           const int cnt)
    : CFSBndCond(pSession, pFields, pTraceNormals, pGridVelocity, pSpaceDim, bcRegion, cnt)
{
    m_diffusionAveWeight = 0.5;

    // Set up rotational boundary edge velocities
    const Array<OneD, const int> &traceBndMap
        = m_fields[0]->GetTraceBndMap();

    int eMax = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetExpSize();
    for (int e = 0; e < eMax; ++e)
    {
        int nBCEdgePts = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
                     GetExp(e)->GetTotPoints();
        int id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[m_offset+e]);

        Array<OneD, NekDouble> x(nBCEdgePts, 0.0);
        Array<OneD, NekDouble> y(nBCEdgePts, 0.0);
        m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetExp(e)->GetCoords(x, y);
        for (int i = 0; i < m_spacedim; i++)
        {
            for (int j = 0; j < nBCEdgePts; ++j)
            {
                m_gridVelocity[i][id2 + j] = (i == 0) ? m_angVel * y[j] : m_angVel * x[j];
            }
        }
    }
}

void WallRotationalBC::v_Apply(
        Array<OneD, Array<OneD, NekDouble> >               &Fwd,
        Array<OneD, Array<OneD, NekDouble> >               &physarray,
        const NekDouble                                    &time)
{
    boost::ignore_unused(time);
    int nVariables = physarray.size();

    // @TODO: ALE we subtract the grid velocity ? "Set u = to ug for this one" - Dave

    int i;
    const Array<OneD, const int> &traceBndMap
        = m_fields[0]->GetTraceBndMap();

    // Take into account that for PDE based shock capturing, eps = 0 at the
    // wall. Adjust the physical values of the trace to take user defined
    // boundaries into account
    int e, id1, id2, nBCEdgePts, eMax;

    eMax = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetExpSize();

    for (e = 0; e < eMax; ++e)
    {
        nBCEdgePts = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
            GetExp(e)->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
            GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[m_offset+e]);

        // Boundary condition for epsilon term. @TODO: Is this correct, or should I do E = p/(gamma -1) + 1/2*rho(u^2 +v^2 + w^2)... or
        if (nVariables == m_spacedim+3)
        {
            Vmath::Zero(nBCEdgePts, &Fwd[nVariables-1][id2], 1);
        }

        for (i = 0; i < m_spacedim; i++)
        {
            // V = -Vin
            //Vmath::Neg(nBCEdgePts, &Fwd[i+1][id2], 1);

            // This now does Vg * rho + Vin
            //Vmath::Vvtvp(nBCEdgePts, &m_gridVelocity[i][id2], 1, &Fwd[0][id2], 1, &Fwd[i+1][id2], 1, &Fwd[i+1][id2], 1);

            for (int j = 0; j < nBCEdgePts; ++j)
            {
                Fwd[i+1][id2 + j] = 2 * m_gridVelocity[i][id2 + j] * Fwd[0][id2 + j] - Fwd[i+1][id2 + j];
            }
        }

        // Copy boundary adjusted values into the boundary expansion
        for (i = 0; i < nVariables; ++i)
        {
            Vmath::Vcopy(nBCEdgePts, &Fwd[i][id2], 1,
                         &(m_fields[i]->GetBndCondExpansions()[m_bcRegion]->
                           UpdatePhys())[id1], 1);
        }
    }
}

}
