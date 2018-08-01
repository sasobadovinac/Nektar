///////////////////////////////////////////////////////////////////////////////
//
// File CrossField.cpp
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
// Description: CrossField solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <ADRSolver/EquationSystems/CrossField.h>

using namespace std;

namespace Nektar
{
string CrossField::className =
    GetEquationSystemFactory().RegisterCreatorFunction("CrossField",
                                                       CrossField::create);

CrossField::CrossField(const LibUtilities::SessionReaderSharedPtr &pSession,
                       const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : Laplace(pSession, pGraph)
{
}

CrossField::~CrossField()
{
}

void CrossField::v_InitObject()
{
    Laplace::v_InitObject();

    ASSERTL0(m_expdim == 2, "Only 2D supported with cross fields.")

    for (int bcRegion = 0;
         bcRegion < m_fields[0]->GetBndConditions().num_elements(); ++bcRegion)
    {
        Array<OneD, Array<OneD, NekDouble>> normals;
        m_fields[0]->GetBoundaryNormals(bcRegion, normals);
        int npts = normals[0].num_elements();

        Array<OneD, Array<OneD, NekDouble>> u(2);
        u[0] = Array<OneD, NekDouble>(npts);
        u[1] = Array<OneD, NekDouble>(npts);

        for (int i = 0; i < normals[0].num_elements(); ++i)
        {
            NekDouble theta = atan2(normals[1][i], normals[0][i]);
            NekDouble t     = fmod(4 * theta, 2 * M_PI);

            u[0][i] = cos(t);
            u[1][i] = sin(t);
        }

        m_fields[0]->GetBndCondExpansions()[bcRegion]->SetPhys(u[0]);
        m_fields[1]->GetBndCondExpansions()[bcRegion]->SetPhys(u[1]);

        m_fields[0]->GetBndCondExpansions()[bcRegion]->FwdTrans_BndConstrained(
            m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys(),
            m_fields[0]->GetBndCondExpansions()[bcRegion]->UpdateCoeffs());
        m_fields[1]->GetBndCondExpansions()[bcRegion]->FwdTrans_BndConstrained(
            m_fields[1]->GetBndCondExpansions()[bcRegion]->GetPhys(),
            m_fields[1]->GetBndCondExpansions()[bcRegion]->UpdateCoeffs());
    }
}

void CrossField::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    const int nPhys   = m_fields[0]->GetNpoints();
    const int nCoeffs = m_fields[0]->GetNcoeffs();

    Array<OneD, Array<OneD, NekDouble>> tmp(m_fields.num_elements());
    for (int i = 0; i < m_fields.num_elements(); ++i)
    {
        tmp[i] = m_fields[i]->GetPhys();
    }

    Array<OneD, NekDouble> psi(nPhys);
    for (int i = 0; i < psi.num_elements(); ++i)
    {
        psi[i] = fmod(atan2(tmp[1][i], tmp[0][i]) / 4.0, 2.0 * M_PI);
    }

    Array<OneD, NekDouble> psiFwd(nCoeffs);
    m_fields[0]->FwdTrans_IterPerExp(psi, psiFwd);

    variables.push_back("psi");
    fieldcoeffs.push_back(psiFwd);
}
}
