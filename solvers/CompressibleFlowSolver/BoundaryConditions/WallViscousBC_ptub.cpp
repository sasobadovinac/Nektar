///////////////////////////////////////////////////////////////////////////////
//
// File: WallViscousBC_ptub.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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

#include "WallViscousBC_ptub.h"

using namespace std;

namespace Nektar
{

std::string WallViscousBC_ptub::classNameViscous = GetCFSBndCondFactory().
    RegisterCreatorFunction("WallViscous_ptub",
                            WallViscousBC_ptub::create,
                            "No-slip (viscous) wall boundary condition.");

std::string WallViscousBC_ptub::classNameAdiabatic = GetCFSBndCondFactory().
    RegisterCreatorFunction("WallAdiabatic_ptub",
                            WallViscousBC_ptub::create,
                            "Adiabatic wall boundary condition.");

WallViscousBC_ptub::WallViscousBC_ptub(const LibUtilities::SessionReaderSharedPtr& pSession,
           const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
           const Array<OneD, Array<OneD, NekDouble> >& pTraceNormals,
           const int pSpaceDim,
           const int bcRegion,
           const int cnt)
    : CFSBndCond(pSession, pFields, pTraceNormals, pSpaceDim, bcRegion, cnt)
{
    m_diffusionAveWeight = 0.5;

    const MultiRegions::ExpListSharedPtr bndexp =
        m_fields[0]->GetBndCondExpansions()[m_bcRegion];

    m_npts    = bndexp->GetTotPoints();
    m_bndPhys = Array<OneD, Array<OneD, NekDouble> > (m_fields.num_elements());
}

void WallViscousBC_ptub::v_Apply(
        Array<OneD, Array<OneD, NekDouble> >               &Fwd,
        Array<OneD, Array<OneD, NekDouble> >               &physarray,
        const NekDouble                                    &time)
{
    int i;
    int nVariables = physarray.num_elements();

    //----------------------------------------------------
    // Update variables on the BC for its time dependence
    std::string varName;
    for (i = 0; i < nVariables; ++i)
    {
        varName = m_session->GetVariable(i);
        m_fields[i]->EvaluateBoundaryConditions(time, varName);
    }

    // Get the variables on the boundary 
    // Merge the two for-loop
    for(i = 0; i < nVariables; ++i) // number of fields
    {
        m_bndPhys[i] = m_fields[i]->GetBndCondExpansions()[m_bcRegion]
            ->UpdatePhys();
    }
    //----------------------------------------------------


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
        id1  = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
            GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[m_offset+e]);

        for (i = 0; i < m_spacedim; i++)
        {
            Vmath::Neg(nBCEdgePts, &Fwd[i+1][id2], 1);
        }

        //------------------------------
        // Super-impose the perturbation 
        // Fwd is created in CompressibleFlowSystem.cpp -> CompressibleFlowSystem::SetBoundaryConditions(...)
        // So it can be modified as we need
        // Fwd is the the array for pyhysical variables, so we can directly manipulate it.    
        for (int i = 0; i < m_spacedim; i++)
        {
            Vmath::Vadd(nBCEdgePts, &m_bndPhys[i+1][id2],1,&Fwd[i+1][id2],1,&Fwd[i+1][id2],1); //id1 or id2 ?
        }
        //------------------------------

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
