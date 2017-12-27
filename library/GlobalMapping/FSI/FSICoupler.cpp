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
}

/**
 * 
 */
void FSICoupler::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
              TiXmlElement                                *pFSI)
{
    ReadBodies(pFSI);
}

/**
 * 
 */
void FSICoupler::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        GlobalMapping::MappingSharedPtr                   &mapping,
        const NekDouble                                   &time)
{
    for (auto &x : m_bodies)
    {
        x->Apply(pFields, m_displFields, time);
    }
}

void FSICoupler::CalculateCoordVel()
{

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
