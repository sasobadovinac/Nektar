///////////////////////////////////////////////////////////////////////////////
//
// File: RigidPlaneCoupler.cpp
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
// Description: FSICoupler with rigid 3DH1D planes deformation
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/FSI/RigidPlaneCoupler.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>

namespace Nektar
{

std::string RigidPlaneCoupler::className =
    GetFSICouplerFactory().RegisterCreatorFunction("RigidPlane",
    RigidPlaneCoupler::create, "Rigid 3DH1D planes deformation");

/**
 * @class RigidPlaneCoupler
 */
RigidPlaneCoupler::RigidPlaneCoupler(
        const LibUtilities::SessionReaderSharedPtr          &pSession,
         const std::weak_ptr<SolverUtils::EquationSystem>   &pEquation)
    : FSICoupler(pSession, pEquation)
{
}

/**
 * 
 */
void RigidPlaneCoupler::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
              TiXmlElement                                *pFSI)
{
    FSICoupler::v_InitObject(pFields, pFSI);

    // Check if using homogeneous expansion
    ASSERTL0( pFields[0]->GetExpType() == MultiRegions::e3DH1D,
            "RigidPlane require 3DH1D expansion");

    // Check if we have only one body
    ASSERTL0( m_bodies.size() == 1, "Need exactly one body for RigidPlane.");

    // Get boundary parameter for the body
    TiXmlElement *body = pFSI->FirstChildElement("BODY");
    std::map<std::string, std::string> vParams;
    TiXmlElement *param = body->FirstChildElement("PARAM");
    while (param)
    {
        ASSERTL0(param->Attribute("NAME"),
                "Missing attribute 'NAME' for parameter in body ");
        std::string nameStr = param->Attribute("NAME");

        ASSERTL0(param->GetText(), "Empty value string for param.");
        std::string valueStr = param->GetText();

        vParams[nameStr] = valueStr;

        param = param->NextSiblingElement("PARAM");
    }
    auto it = vParams.find("Boundary");
    ASSERTL0(it != vParams.end(), "Missing parameter 'Boundary'.");
    ASSERTL0(it->second.length() > 0, "Empty parameter 'Boundary'.");
    std::string boundaryString = it->second;

    // Parse the boundary regions into a list.
    std::string::size_type firstInd =
                            boundaryString.find_first_of('[') + 1;
    std::string::size_type lastInd =
                            boundaryString.find_last_of(']') - 1;

    ASSERTL0(firstInd <= lastInd,
            (std::string("Error reading boundary region definition:") +
             boundaryString).c_str());

    std::string indString =
            boundaryString.substr(firstInd, lastInd - firstInd + 1);
    std::vector<unsigned int> boundaryRegionsIdList;
    ParseUtils::GenerateSeqVector(indString, boundaryRegionsIdList);

    // Determine boundary condition corresponding to first bnd region in list
    SpatialDomains::BoundaryConditions bcs(m_session, pFields[0]->GetGraph());
    const SpatialDomains::BoundaryRegionCollection &bregions =
                                            bcs.GetBoundaryRegions();
    int cnt = 0;
    m_bndId = -1;
    for (auto &it : bregions)
    {
        if ( it.first ==  boundaryRegionsIdList[0])
        {
            m_bndId = cnt;
        }
        cnt++;
    }

    // In parallel, not all processes will have this boundary
    //    Therefore, pick one to broadcast the result
    LibUtilities::CommSharedPtr comm =
            m_displFields[0]->GetComm()->GetRowComm();
    m_bcastRank = -1;
    if (m_bndId != -1)
    {
        m_bcastRank = comm->GetRank();
    }
    comm->AllReduce(m_bcastRank, LibUtilities::ReduceMax);
    ASSERTL0(m_bcastRank >= 0, "Boundary not found.");
}

void RigidPlaneCoupler::v_CalculateDisplacement()
{
    LibUtilities::CommSharedPtr comm =
            m_displFields[0]->GetComm()->GetRowComm();

    // Process m_bcastRank gets the displacements and then broadcasts them
    int nPlanes = m_displFields[0]->GetZIDs().num_elements();
    Array<OneD, NekDouble> displ(m_expDim * nPlanes, 0.0);
    if( comm->GetRank() == m_bcastRank)
    {
        for(int i = 0; i < m_expDim; ++i)
        {
            for(int p = 0; p < nPlanes; ++p)
            {
                displ[i + p * m_expDim] = m_displFields[i]->
                    GetBndCondExpansions()[m_bndId]->GetPlane(p)->GetPhys()[0];
            }
        }
    }
    comm->Bcast(displ, m_bcastRank);

    // Fill m_displFields with constant displacement per plane
    for(int i = 0; i < m_expDim; ++i)
    {
        for(int p = 0; p < nPlanes; ++p)
        {
            int nPts = m_displFields[i]->GetPlane(p)->GetTotPoints();
            Vmath::Fill(nPts, displ[i + p * m_expDim],
                    m_displFields[i]->GetPlane(p)->UpdatePhys(), 1);
        }
    }
}

}
