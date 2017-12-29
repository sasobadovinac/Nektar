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

#include <GlobalMapping/FSI/FSICoupler.h>

using namespace std;

namespace Nektar
{
namespace GlobalMapping
{

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
    if(intMethod == "IMEXOrder1")
    {
        m_intSteps = 1;
    }
    else if (intMethod == "IMEXOrder2")
    {
        m_intSteps = 2;
    }
    else if (intMethod == "IMEXOrder3")
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

    // Create entries for m_bodies
    ReadBodies(pFSI);

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
        m_oldCoords[j] = Array<OneD, Array<OneD, NekDouble>> (3);
        for( int i = 0; i < 3; ++i)
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
}

/**
 * 
 */
void FSICoupler::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        GlobalMapping::MappingSharedPtr                   &mapping,
        const NekDouble                                   &time)
{
    // Call m_bodies to update the boundary conditions of m_displFields
    for (auto &x : m_bodies)
    {
        x->Apply(pFields, m_displFields, time);
    }

    CalculateDisplacement();

    UpdateCoordinates();

    CalculateCoordVel();

    mapping->UpdateMapping(time, m_coords, m_coordsVel);
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
    for( int i = m_expDim; i < m_spaceDim; ++i)
    {
        m_coords[i] = m_meshCoords[i];
    }
}

void FSICoupler::CalculateCoordVel()
{
    // TO DO: Calculate m_coordsVel and update m_oldCoords
}

void FSICoupler::ReadBodies(TiXmlElement* pFSI)
{
    m_bodies.clear();

    TiXmlElement *body = pFSI->FirstChildElement("BODY");
    while (body)
    {
        ASSERTL0(body->Attribute("TYPE"),
                "Missing attribute 'TYPE' for filter.");
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
                        typeStr, m_session, m_displFields, vParams));

        body = body->NextSiblingElement("BODY");
    }
}

}
}
