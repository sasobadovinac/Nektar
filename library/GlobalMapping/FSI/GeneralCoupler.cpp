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

#include <GlobalMapping/FSI/GeneralCoupler.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>

namespace Nektar
{
namespace GlobalMapping
{

std::string GeneralCoupler::className =
    GetFSICouplerFactory().RegisterCreatorFunction("General",
    GeneralCoupler::create, "General domain deformation");

/**
 * @class GeneralCoupler
 */
GeneralCoupler::GeneralCoupler(
        const LibUtilities::SessionReaderSharedPtr          &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields)
    : FSICoupler(pSession, pFields)
{
}

/**
 *
 */
void GeneralCoupler::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields,
              TiXmlElement* pFSI)
{
    FSICoupler::v_InitObject(pFields, pFSI);

    // Create m_displFields
    m_displFields = Array<OneD, MultiRegions::ExpListSharedPtr> (m_expDim);
    const SpatialDomains::MeshGraphSharedPtr graph = pFields[0]->GetGraph();
    string fieldNames[3] = {"x", "y", "z"};
    switch (pFields[0]->GetExpType())
    {
        case MultiRegions::e2D:
        {
            MultiRegions::ContField2DSharedPtr tmp =
                std::dynamic_pointer_cast<
                    MultiRegions::ContField2D>(pFields[0]);

            for(int i = 0; i < m_expDim; ++i)
            {
                m_displFields[i] =
                    MemoryManager<MultiRegions::ContField2D>::
                        AllocateSharedPtr(*tmp, graph, fieldNames[i]);
            }
        }
        break;
        case MultiRegions::e3D:
        {
            MultiRegions::ContField3DSharedPtr tmp =
                std::dynamic_pointer_cast<
                    MultiRegions::ContField3D>(pFields[0]);

            for(int i = 0; i < m_expDim; ++i)
            {
                m_displFields[i] =
                    MemoryManager<MultiRegions::ContField3D>::
                        AllocateSharedPtr(*tmp, graph, fieldNames[i]);
            }
        }
        break;
        case MultiRegions::e3DH1D:
        {
            MultiRegions::ContField3DHomogeneous1DSharedPtr tmp =
                std::dynamic_pointer_cast<
                    MultiRegions::ContField3DHomogeneous1D>(pFields[0]);

            for(int i = 0; i < m_expDim; ++i)
            {
                m_displFields[i] =
                    MemoryManager<MultiRegions::ContField3DHomogeneous1D>::
                        AllocateSharedPtr(*tmp, graph, fieldNames[i]);
            }
        }
        break;
        default:
            ASSERTL0(0,"Dimension not supported");
        break;
    }
}

void GeneralCoupler::v_CalculateDisplacement()
{
    // TO DO: solve Laplace equation to obtain displacements
}

}
}
