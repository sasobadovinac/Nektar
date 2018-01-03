///////////////////////////////////////////////////////////////////////////////
//
// File: FSICoupler.cpp
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
// Description: Abstract base class for fluid-structure interaction coupler.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/FSI/FSICoupler.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>

using namespace std;

namespace Nektar
{

NekDouble FSICoupler::BDF_Alpha_Coeffs[3][3] = {
    { 1.0,  0.0, 0.0},{ 2.0, -0.5, 0.0},{ 3.0, -1.5, 1.0/3.0}};
NekDouble FSICoupler::BDF_Gamma0_Coeffs[3] = {
      1.0,  1.5, 11.0/6.0};

FSICouplerFactory& GetFSICouplerFactory()
{
    static FSICouplerFactory instance;
    return instance;
}

FSICoupler::FSICoupler(
        const LibUtilities::SessionReaderSharedPtr         &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr>  &pFields)
    : m_session(pSession)
{
    // Get timestep
    m_session->LoadParameter("TimeStep", m_timestep);

    switch (pFields[0]->GetExpType())
    {
        case MultiRegions::e2D:
        {
            m_expDim   = 2;
            m_spaceDim = 2;
        }
        break;
        case MultiRegions::e3D:
        {
            m_expDim   = 3;
            m_spaceDim = 3;
        }
        break;
        case MultiRegions::e3DH1D:
        {
            m_expDim   = 2;
            m_spaceDim = 3;
        }
        break;
        default:
            ASSERTL0(0,"Dimension not supported");
        break;
    }

    std::string intMethod = m_session->GetSolverInfo("TIMEINTEGRATIONMETHOD");
    if(boost::iequals(intMethod, "IMEXOrder1"))
    {
        m_intSteps = 1;
    }
    else if (boost::iequals(intMethod, "IMEXOrder2"))
    {
        m_intSteps = 2;
    }
    else if (boost::iequals(intMethod, "IMEXOrder3"))
    {
        m_intSteps = 3;
    }
    else
    {
        ASSERTL0(false, "Time integration method not supported.");
    }
}

/**
 * 
 */
void FSICoupler::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
              TiXmlElement                                *pFSI)
{
    // Adjust mapping
    GlobalMapping::MappingSharedPtr mapping =
            GlobalMapping::Mapping::Load(m_session, pFields);
    mapping->SetTimeDependent(true);
    mapping->SetFromFunction (false);

    // Create m_displFields
    CreateDisplacementFields(pFields);

    // Create entries for m_bodies
    ReadBodies(pFields, pFSI);

    // Allocate storage
    int nPts     = pFields[0]->GetTotPoints();
    m_coords     = Array<OneD, Array<OneD, NekDouble>> (3);
    m_meshCoords = Array<OneD, Array<OneD, NekDouble>> (3);
    m_coordsVel  = Array<OneD, Array<OneD, NekDouble>> (3);
    for( int i = 0; i < 3; ++i)
    {
        m_coords[i]     = Array<OneD, NekDouble> (nPts, 0.0);
        m_meshCoords[i] = Array<OneD, NekDouble> (nPts, 0.0);
        m_coordsVel[i]  = Array<OneD, NekDouble> (nPts, 0.0);
    }

    m_oldCoords = Array<OneD, Array<OneD, Array<OneD, NekDouble>>> (m_intSteps);
    for( int j = 0; j < m_intSteps; ++j)
    {
        m_oldCoords[j] = Array<OneD, Array<OneD, NekDouble>> (m_expDim);
        for( int i = 0; i < m_expDim; ++i)
        {
            m_oldCoords[j][i] = Array<OneD, NekDouble> (nPts, 0.0);
        }
    }

    // Get mesh coordinates
    m_meshCoords = Array<OneD, Array<OneD, NekDouble>> (3);
    for( int i = 0; i < 3; ++i)
    {
        m_meshCoords[i] = Array<OneD, NekDouble> (nPts, 0.0);
    }
    pFields[0]->GetCoords(m_meshCoords[0],m_meshCoords[1],m_meshCoords[2]);

    // Initialise m_coords to m_meshCoords (zero displacement in IC)
    for( int i = 0; i < 3; ++i)
    {
        Vmath::Vcopy(nPts, m_meshCoords[i], 1, m_coords[i], 1);
    }

    // Set m_oldCoords[0] to m_coords for calculating coordinates velocity
    for( int i = 0; i < m_expDim; ++i)
    {
        Vmath::Vcopy(nPts, m_coords[i], 1, m_oldCoords[0][i], 1);
    }
}

void FSICoupler::CreateDisplacementFields(
    const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields)
{
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

    // Initialise to zero
    for(int i = 0; i < m_expDim; ++i)
    {
        Vmath::Zero ( m_displFields[i]->GetTotPoints(),
                        m_displFields[i]->UpdatePhys(), 1);
        Vmath::Zero ( m_displFields[i]->GetNcoeffs(),
                        m_displFields[i]->UpdateCoeffs(), 1);
    }
}

/**
 * 
 */
void FSICoupler::v_UpdateMappingCoordsVel(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        GlobalMapping::MappingSharedPtr                   &mapping,
        const NekDouble                                   &time)
{
    m_time = time;
    // Call m_bodies to update the boundary conditions of m_displFields
    for (auto &x : m_bodies)
    {
        x->Apply(pFields, m_displFields, time);
    }

    CalculateDisplacement();

    UpdateCoordinates();

    CalculateCoordVel();

    mapping->UpdateMappingCoordsVel(time, m_coordsVel);
}

/**
 * 
 */
void FSICoupler::v_UpdateMappingCoords(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        GlobalMapping::MappingSharedPtr                   &mapping,
        const NekDouble                                   &time)
{
    ASSERTL1(time == m_time, "UpdateMappingCoordsVel should be called"
            " before UpdateMappingCoords");
    mapping->UpdateMappingCoords(time, m_coords);
}

void FSICoupler::UpdateCoordinates()
{
    for( int i = 0; i < m_expDim; ++i)
    {
        Vmath::Vadd(m_coords[i].num_elements(),
                m_displFields[i]->GetPhys(), 1,
                m_meshCoords[i], 1,
                m_coords[i], 1);
    }
}

void FSICoupler::CalculateCoordVel()
{
    int nPts = m_coords[0].num_elements();

    // Determine correct order (for initial time-steps)
    static int nCalls = 0;
    ++nCalls;
    int order = min(nCalls,m_intSteps);

    for( int i = 0; i < m_expDim; ++i)
    {
        // Calculate m_coordsVel
        Vmath::Smul(nPts,
                    BDF_Gamma0_Coeffs[order-1],
                    m_coords[i], 1,
                    m_coordsVel[i],  1);

        for(int j = 0; j < order; j++)
        {
            Vmath::Svtvp(nPts,
                        -1*BDF_Alpha_Coeffs[order-1][j],
                         m_oldCoords[j][i], 1,
                         m_coordsVel[i],    1,
                         m_coordsVel[i],    1);
        }
        // Divide by time-step
        Vmath::Smul(nPts, 1.0/m_timestep, m_coordsVel[i], 1,
                                          m_coordsVel[i], 1);

        // Update m_oldCoords
        Array<OneD, NekDouble> tmp = m_oldCoords[m_intSteps-1][i];
        for(int n = m_intSteps-1; n > 0; --n)
        {
            m_oldCoords[n][i] = m_oldCoords[n-1][i];
        }
        m_oldCoords[0][i] = tmp;
        Vmath::Vcopy(nPts, m_coords[i], 1, m_oldCoords[0][i], 1);
    }
}

void FSICoupler::ReadBodies(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
          TiXmlElement* pFSI)
{
    m_bodies.clear();

    TiXmlElement *body = pFSI->FirstChildElement("BODY");
    while (body)
    {
        ASSERTL0(body->Attribute("TYPE"),
                "Missing attribute 'TYPE' for body.");
        std::string typeStr = body->Attribute("TYPE");

        std::map<std::string, std::string> vParams;

        TiXmlElement *param = body->FirstChildElement("PARAM");
        while (param)
        {
            ASSERTL0(param->Attribute("NAME"),
                    "Missing attribute 'NAME' for parameter in body "
                    + typeStr + "'.");
            std::string nameStr = param->Attribute("NAME");

            ASSERTL0(param->GetText(), "Empty value string for param.");
            std::string valueStr = param->GetText();

            vParams[nameStr] = valueStr;

            param = param->NextSiblingElement("PARAM");
        }

        m_bodies.push_back(
            GetFSIBodyFactory().CreateInstance(
                        typeStr, m_session, pFields, vParams));

        body = body->NextSiblingElement("BODY");
    }
}

}
