///////////////////////////////////////////////////////////////////////////////
//
// File: GeneralCoupler.cpp
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
// Description: FSICoupler with general domain deformation
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/FSI/GeneralCoupler.h>

namespace Nektar
{

std::string GeneralCoupler::className =
    GetFSICouplerFactory().RegisterCreatorFunction("General",
    GeneralCoupler::create, "General domain deformation");

/**
 * @class GeneralCoupler
 */
GeneralCoupler::GeneralCoupler(
        const LibUtilities::SessionReaderSharedPtr          &pSession,
        const std::weak_ptr<SolverUtils::EquationSystem>   &pEquation)
    : FSICoupler(pSession, pEquation)
{
}

void GeneralCoupler::v_CalculateDisplacement()
{
    int nPts = m_coords[0].num_elements();
    // Solve Laplace equation to obtain displacements
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 0.0;
    Array<OneD, NekDouble> forcing(nPts, 0.0);
    for(int i = 0; i < m_expDim; ++i)
    {
        m_displFields[i]->HelmSolve(forcing, m_displFields[i]->UpdateCoeffs(),
                NullFlagList, factors);
        // Bwd transform the result
        m_displFields[i]->BwdTrans(m_displFields[i]->GetCoeffs(),
                                   m_displFields[i]->UpdatePhys());
    }
}

}
