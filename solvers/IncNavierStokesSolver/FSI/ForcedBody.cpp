///////////////////////////////////////////////////////////////////////////////
//
// File: ForcedBody.cpp
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
// Description: Forced body motion for FSI problems
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/FSI/ForcedBody.h>

namespace Nektar
{

std::string ForcedBody::className =
    GetFSIBodyFactory().RegisterCreatorFunction("Forced",
    ForcedBody::create, "Forced body motion");

/**
 * @class ForcedBody
 */
ForcedBody::ForcedBody(
        const LibUtilities::SessionReaderSharedPtr          &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields)
    : FSIBody(pSession, pFields)
{
}

/**
 *
 */
void ForcedBody::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields,
        const std::map<std::string, std::string>            &pParams)
{
    FSIBody::v_InitObject(pFields, pParams);

    // Read displacement parameter
    auto it = pParams.find("DisplacementFCN");
    ASSERTL0(it != pParams.end(),     "Missing parameter 'DisplacementFCN'.");
    ASSERTL0(it->second.length() > 0, "Empty parameter 'DisplacementFCN'.");
    ASSERTL0(m_session->DefinesFunction(it->second),
            "Function '" + it->second + "' not defined.");
    m_funcName = it->second;
}

void ForcedBody::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pDisplFields,
        const NekDouble                                      &time)
{
    int    dim           = pDisplFields.num_elements();
    string fieldNames[3] = {"x", "y", "z"};

    // Loop coordinates
    for( int i = 0; i < dim; ++i)
    {
        // Skip this direction if function is not defined for it
        if( !m_session->DefinesFunction(m_funcName, fieldNames[i]))
        {
            continue;
        }

        // Get boundary expansions
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr> bndConds =
                pDisplFields[i]->GetBndConditions();
        Array<OneD, MultiRegions::ExpListSharedPtr> bndExp =
                pDisplFields[i]->GetBndCondExpansions();

        // Get function for displacement in this direction
        LibUtilities::EquationSharedPtr ffunc =
            m_session->GetFunction(m_funcName, fieldNames[i]);

        // Loop on boundary regions
        for( int n = 0; n < bndConds.num_elements(); ++n)
        {
            // Only modify boundary regions corresponding to this body
            if(m_boundaryRegionIsInList[n] == 1)
            {
                // Number of points on this boundary
                int nPts = bndExp[n]->GetTotPoints();

                // Get coordinates
                Array<OneD, Array<OneD, NekDouble>> coords(3);
                for (int j = 0; j < 3; ++j)
                {
                    coords[j] = Array<OneD, NekDouble>(nPts, 0.0);
                }
                bndExp[n]->GetCoords(coords[0], coords[1], coords[2]);

                // Evaluate function
                ffunc->Evaluate(coords[0], coords[1], coords[2], time,
                        bndExp[n]->UpdatePhys());

                // Update coefficients
                bndExp[n]->FwdTrans_BndConstrained(bndExp[n]->GetPhys(),
                                            bndExp[n]->UpdateCoeffs());
                if (pDisplFields[i]->GetExpType() == MultiRegions::e3DH1D)
                {
                    bndExp[n]->HomogeneousFwdTrans(bndExp[n]->GetCoeffs(),
                                            bndExp[n]->UpdateCoeffs());
                }
            }
        }
    }
}

}
