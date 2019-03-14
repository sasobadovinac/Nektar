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
                ReadInterfaces(interfacesTag);
            }
        }

        std::string ReadTag(std::string tagStr)
        {
            std::string::size_type indxBeg = tagStr.find_first_of('[') + 1;
            std::string::size_type indxEnd = tagStr.find_last_of(']') - 1;

            ASSERTL0(indxBeg <= indxEnd,
                     (std::string("Error reading boundary region definition:")
                     + tagStr).c_str());

            std::string indxStr = tagStr.substr(indxBeg, indxEnd - indxBeg + 1);

            return indxStr;
        }

        void Interfaces::ReadInterfaces(TiXmlElement *interfacesTag)
        {
            TiXmlElement *interfaceElementTag = interfacesTag->FirstChildElement();
            while(interfaceElementTag)
            {
                std::vector<std::string>::iterator iter;
                std::string interfaceType = interfaceElementTag->Value();

                int err; //variable to check attributes are read correctly
                int indx; //value that holds interface ID
                err = interfaceElementTag->QueryIntAttribute("ID", &indx);
                ASSERTL0(err == TIXML_SUCCESS, "Unable to read interface ID.");

                std::string interfaceMovingDomainStr;
                err = interfaceElementTag->QueryStringAttribute("DOMAIN",
                                                                &interfaceMovingDomainStr);
                ASSERTL0(err == TIXML_SUCCESS,
                         "Unable to read moving interface domain.");
                CompositeMap movingDomain = m_meshGraph->GetDomain(
                        stoi(ReadTag(interfaceMovingDomainStr)));

                std::string interfaceFixedDomainStr;
                err = interfaceElementTag->QueryStringAttribute("FIXED",
                                                                &interfaceFixedDomainStr);
                ASSERTL0(err == TIXML_SUCCESS,
                         "Unable to read fixed interface domain.");
                CompositeMap fixedDomain = m_meshGraph->GetDomain(
                        stoi(ReadTag(interfaceFixedDomainStr)));

                std::string interfaceEdgeStr;
                err = interfaceElementTag->QueryStringAttribute("INTERFACE",
                                                                &interfaceEdgeStr);
                ASSERTL0(err == TIXML_SUCCESS,
                         "Unable to read interface composite.");

                InterfaceEdgeShPtr interfaceEdge(
                        MemoryManager<InterfaceEdge>::AllocateSharedPtr());

                interfaceEdgeStr = ReadTag(interfaceEdgeStr);
                m_meshGraph->GetCompositeList(interfaceEdgeStr, *interfaceEdge);

                if (interfaceType == "R")
                {

                    std::string originStr;
                    err = interfaceElementTag->QueryStringAttribute("ORIGIN",
                                                                    &originStr);

                    ASSERTL0(err == TIXML_SUCCESS,
                             "Unable to read origin.");

                    std::vector<NekDouble> originVec;
                    ParseUtils::GenerateVector(originStr, originVec);
                    PointGeomSharedPtr origin = MemoryManager<PointGeom>::
                            AllocateSharedPtr(PointGeom(3, 0, originVec[0], originVec[1], originVec[2]));

                    std::string axisStr;
                    err = interfaceElementTag->QueryStringAttribute("AXIS",
                                                                    &axisStr);

                    ASSERTL0(err == TIXML_SUCCESS,
                             "Unable to read axis.");

                    std::vector<NekDouble> axis;
                    ParseUtils::GenerateVector(originStr, axis);

                    std::string angularVelStr;
                    err = interfaceElementTag->QueryStringAttribute("ANGVEL",
                                                                    &angularVelStr);

                    ASSERTL0(err == TIXML_SUCCESS,
                             "Unable to read angular velocity.");

                    NekDouble angularVel = stod(angularVelStr);

                    InterfaceShPtr rotatingInterface(
                            MemoryManager<RotatingInterface>::AllocateSharedPtr(
                                    m_session, movingDomain, fixedDomain,
                                    interfaceEdge, origin, axis, angularVel));

                    m_interfaces[indx] = rotatingInterface;
                }

                interfaceElementTag = interfaceElementTag->NextSiblingElement();
            }
            //
        }
    }
}
