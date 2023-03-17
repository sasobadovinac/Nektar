///////////////////////////////////////////////////////////////////////////////
//
// File: StagnationInflowBC.cpp
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
// Description: Stagnation conditions inflow boundary condition
//
///////////////////////////////////////////////////////////////////////////////

#include "StagnationInflowBC.h"
#include <boost/core/ignore_unused.hpp>

using namespace std;

namespace Nektar
{

std::string StagnationInflowBC::className =
    GetCFSBndCondFactory().RegisterCreatorFunction(
        "StagnationInflow", StagnationInflowBC::create,
        "Stagnation conditions inflow boundary condition.");

StagnationInflowBC::StagnationInflowBC(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
    const int pSpaceDim, const int bcRegion, const int cnt)
    : CFSBndCond(pSession, pFields, pTraceNormals, pSpaceDim, bcRegion, cnt)
{
    const size_t nvariables = m_fields.size();
    const int expdim        = m_fields[0]->GetGraph()->GetMeshDimension();
    const int spacedim      = m_fields[0]->GetGraph()->GetSpaceDimension();
    const int numBCPts =
        m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetNpoints();

    // Allocate internal storage for boundary condition values
    m_fieldStorage = Array<OneD, Array<OneD, NekDouble>>(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        m_fieldStorage[i] = Array<OneD, NekDouble>(numBCPts, 0.0);
        Vmath::Vcopy(numBCPts,
                     m_fields[i]->GetBndCondExpansions()[m_bcRegion]->GetPhys(),
                     1, m_fieldStorage[i], 1);
    }

    // Compute the flow direction that will be imposed
    for (int j = 0; j < numBCPts; ++j)
    {
        // Compute norm of user-specified flow direction
        NekDouble dirNorm = 0.;
        for (int i = 0; i < spacedim; ++i)
        {
            dirNorm += std::pow(m_fieldStorage[i + 1][j], 2);
        }
        dirNorm = std::sqrt(dirNorm);

        // Use the user-specified flow direction if it's nonzero
        if (dirNorm > 1.E-8)
        {
            for (int i = 0; i < spacedim; ++i)
            {
                m_fieldStorage[i + 1][j] /= dirNorm;
            }
        }
        else
        {
            for (int i = 0; i < expdim; ++i)
            {
                m_fieldStorage[i + 1][j] = -m_traceNormals[i][j];
            }
            for (int i = expdim; i < spacedim; ++i)
            {
                m_fieldStorage[i + 1][j] = 0.;
            }
        }
    }
}

void StagnationInflowBC::v_Apply(Array<OneD, Array<OneD, NekDouble>> &Fwd,
                                 Array<OneD, Array<OneD, NekDouble>> &physarray,
                                 const NekDouble &time)
{
    boost::ignore_unused(time);

    const size_t nTracePts  = Fwd[0].size();
    const size_t nVariables = physarray.size();

    ASSERTL0(nTracePts == m_fields[0]->GetTrace()->GetNpoints(),
             "Number of trace points does not match in "
             "StagnationInflowBC::v_Apply()");

    const Array<OneD, const int> &traceBndMap = m_fields[0]->GetTraceBndMap();

    // Get velocity magnitude from Fwd
    Array<OneD, NekDouble> absVel(nTracePts, 0.0);
    m_varConv->GetAbsoluteVelocity(Fwd, absVel);

    // Loop over all elements on the boundary
    for (int e = 0;
         e < m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetExpSize(); ++e)
    {
        // Number of quadrature points in this element
        const int npts = m_fields[0]
                             ->GetBndCondExpansions()[m_bcRegion]
                             ->GetExp(e)
                             ->GetTotPoints();
        const int id1 =
            m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetPhys_Offset(e);
        const int id2 =
            m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[m_offset + e]);

        // Loop over all quadrature points on this element
        for (int i = 0; i < npts; ++i)
        {
            // Copy conserved variables at stagnation state from BC
            NekDouble rhoStag = m_fieldStorage[0][id1 + i];
            NekDouble EStag   = m_fieldStorage[nVariables - 1][id1 + i];

            // Compute stagnation enthalpy
            NekDouble hStag = m_gamma * EStag / rhoStag;

            // Compute static enthalpy by solving 1D energy balance
            NekDouble hStat = hStag - 0.5 * pow(absVel[id2 + i], 2);

            // Compute density from isentropic flow relation
            NekDouble rho = rhoStag * pow(rhoStag * hStat / (EStag * m_gamma),
                                          1. / (m_gamma - 1.));

            // Update density for BC
            (m_fields[0]
                 ->GetBndCondExpansions()[m_bcRegion]
                 ->UpdatePhys())[id1 + i] = rho;

            // Update momentum for BC
            for (int n = 0; n < m_spacedim; ++n)
            {
                (m_fields[1 + n]
                     ->GetBndCondExpansions()[m_bcRegion]
                     ->UpdatePhys())[id1 + i] =
                    rho * absVel[id2 + i] * m_fieldStorage[1 + n][id1 + i];
            }

            // Update total energy for BC
            (m_fields[nVariables - 1]
                 ->GetBndCondExpansions()[m_bcRegion]
                 ->UpdatePhys())[id1 + i] =
                rho * (hStat / m_gamma + 0.5 * pow(absVel[id2 + i], 2));
        }
    }
}

} // namespace Nektar
