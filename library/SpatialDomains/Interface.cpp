////////////////////////////////////////////////////////////////////////////////
//
//  File: InterfaceConditions.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <SpatialDomains/Interface.h>
#include <tinyxml.h>

using namespace std;

namespace Nektar
{
    namespace SpatialDomains
    {
        Interfaces::Interfaces(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const MeshGraphSharedPtr &meshGraph) :
                m_meshGraph(meshGraph), m_session(pSession)
        {
            Read(m_session->GetElement("NEKTAR/CONDITIONS/INTERFACES"));
        }

        Interfaces::Interfaces(void)
        {
        }

        Interfaces::~Interfaces(void)
        {
        }

        void Interfaces::Read(TiXmlElement *interfacesTag)
        {
            ASSERTL0(interfacesTag, "Unable to find INTERFACES tag in file.");

            if (interfacesTag)
            {
                ReadInterface(interfacesTag);
            }
        }

        std::string ReadTag(std::string domainStr)
        {
            std::string::size_type indxBeg = domainStr.find_first_of('[') + 1;
            std::string::size_type indxEnd = domainStr.find_last_of(']') - 1;

            ASSERTL0(indxBeg <= indxEnd,
                     (std::string("Error reading boundary region definition:") + domainStr).c_str());

            std::string indxStr = domainStr.substr(indxBeg, indxEnd - indxBeg + 1);

            return indxStr;
        }

        void Interfaces::ReadInterfaces(TiXmlElement *interfacesTag)
        {
            TiXmlElement *interfaceElementTag = interfacesTag->FirstChildElement();
            std::string interfaceType = interfaceElementTag->Value();

            int err; //variable to check attributes are read correctly
            int indx; //value that holds interface ID
            err = interfaceElementTag->QueryIntAttribute("ID", &indx);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read interface ID.");

            std::string interfaceMovingDomainStr;
            err = interfaceElementTag->QueryStringAttribute("DOMAIN", &interfaceMovingDomainStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read moving interface domain.");
            CompositeMap movingDomain = m_meshGraph->GetDomain(stoi(ReadTag(interfaceMovingDomainStr)));

            std::string interfaceFixedDomainStr;
            err = interfaceElementTag->QueryStringAttribute("FIXED", &interfaceFixedDomainStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read fixed interface domain.");
            CompositeMap fixedDomain = m_meshGraph->GetDomain(stoi(ReadTag(interfaceFixedDomainStr)));

            std::string interfaceEdgeStr;
            err = interfaceElementTag->QueryStringAttribute("INTERFACE", &interfaceEdgeStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read interface composite.");

            InterfaceEdgeShPtr interfaceEdge;
            m_meshGraph->GetCompositeList(ReadTag(interfaceEdgeStr), *interfaceEdge);

            TiXmlAttribute *attr = interfaceElementTag->FirstAttribute();

            if (interfaceType == "R")
            {

            }
        }
    }
}
