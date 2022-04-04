////////////////////////////////////////////////////////////////////////////////
//
//  File: Movement.cpp
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
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/Movement/Movement.h>
#include <tinyxml.h>

std::string ReadTag(std::string &tagStr)
{
    std::string::size_type indxBeg = tagStr.find_first_of('[') + 1;
    std::string::size_type indxEnd = tagStr.find_last_of(']') - 1;

    ASSERTL0(
        indxBeg <= indxEnd,
        (std::string("Error reading interface region definition:") + tagStr)
            .c_str());

    std::string indxStr = tagStr.substr(indxBeg, indxEnd - indxBeg + 1);

    return indxStr;
}

namespace Nektar
{
namespace SpatialDomains
{

Movement::Movement(const LibUtilities::SessionReaderSharedPtr &pSession,
                   const MeshGraphSharedPtr &meshGraph)
    : m_meshGraph(meshGraph), m_session(pSession)
{
    TiXmlElement *xmlDoc     = m_session->GetElement("NEKTAR");
    TiXmlNode *conditionsXml = xmlDoc->FirstChild("CONDITIONS");
    if (conditionsXml != nullptr) // @TODO: Can I remove this line?
    {
        TiXmlNode *movement = conditionsXml->FirstChild("MOVEMENT");
        if (movement != nullptr)
        {
            if (movement->FirstChild("ZONES") != nullptr)
            {
                ReadZones(
                    m_session->GetElement("NEKTAR/CONDITIONS/MOVEMENT/ZONES"));
            }

            if (movement->FirstChild("INTERFACES") != nullptr)
            {
                ReadInterfaces(m_session->GetElement(
                    "NEKTAR/CONDITIONS/MOVEMENT/INTERFACES"));
            }

            // @TODO: Put in a check that both zones and interfaces are defined
        }
    }

    // DEBUG COMMENTS
    if (conditionsXml != nullptr) // Set if verbose/debug mode? to output rank interface information
    {
        if (m_session->GetComm()->GetRank() == 0 &&
            m_session->DefinesCmdLineArgument("verbose"))
        {
            std::cout << "Num zones: " << m_zones.size();
            std::cout << "\n-----------------------------\n";
            std::cout << "Zone ID "
                      << "Type \t"
                      << "# elmts \n";
            std::cout << "-----------------------------" << std::endl;

            std::array<std::string, 4> edgeName = {{"n", "l", "r", "b"}};

            for (auto &zone : m_zones)
            {
                std::cout << zone.first << "\t"
                          << MovementTypeStr[static_cast<int>(
                                 zone.second->GetMovementType())]
                          << "\t" << zone.second->GetElements().size()
                          << std::endl;
            }
            std::cout << "-----------------------------" << std::endl
                      << std::endl;

            std::cout << "Num interfaces: " << m_interfaces.size();
            std::cout << "\n-----------------------------\n";
            std::cout << "Name \t"
                      << "L # elmts \t"
                      << "R # elmts \n";
            std::cout << "-----------------------------" << std::endl;

            for (auto &interface : m_interfaces)
            {
                std::cout
                    << interface.first.second << "\t"
                    << interface.second->GetLeftInterface()->GetEdgeIds().size()
                    << "\t"
                    << interface.second->GetRightInterface()
                           ->GetEdgeIds()
                           .size()
                    << std::endl;
            }
            std::cout << "-----------------------------" << std::endl;
        }
    }
}

void Movement::ReadZones(TiXmlElement *zonesTag)
{
    int coordDim = m_meshGraph->GetSpaceDimension();

    ASSERTL0(zonesTag, "Unable to find ZONES tag in file.");
    TiXmlElement *zonesElement = zonesTag->FirstChildElement();
    while (zonesElement)
    {
        std::string zoneType = zonesElement->Value();

        int err;
        int indx;

        err = zonesElement->QueryIntAttribute("ID", &indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read zone ID.");

        std::string interfaceDomainStr;
        err = zonesElement->QueryStringAttribute("DOMAIN", &interfaceDomainStr);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read zone domain.");

        auto &domains = m_meshGraph->GetDomain();
        auto domFind  = stoi(ReadTag(interfaceDomainStr));
        std::map<int, CompositeSharedPtr> domain;
        if (domains.find(domFind) != domains.end())
        {
            domain = domains.at(domFind);
        }

        ZoneBaseShPtr zone;
        if (zoneType == "R" || zoneType == "ROTATE" || zoneType == "ROTATING")
        {
            std::string originStr;
            err = zonesElement->QueryStringAttribute("ORIGIN", &originStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read origin.");
            std::vector<NekDouble> originVec;
            ParseUtils::GenerateVector(originStr, originVec);
            NekPoint<NekDouble> origin =
                NekPoint<NekDouble>(originVec[0], originVec[1], originVec[2]);

            std::string axisStr;
            err = zonesElement->QueryStringAttribute("AXIS", &axisStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read axis.");
            DNekVec axis(axisStr);
            axis.Normalize();

            std::string angularVelStr;
            err = zonesElement->QueryStringAttribute("ANGVEL", &angularVelStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read angular velocity.");

            m_session->SubstituteExpressions(angularVelStr);
            LibUtilities::EquationSharedPtr angularVelEqn =
                MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(
                    m_session->GetInterpreter(), angularVelStr);

            zone = ZoneRotateShPtr(MemoryManager<ZoneRotate>::AllocateSharedPtr(
                indx, domain, coordDim, origin, axis, angularVelEqn));

            m_moveFlag = true;
        }
        else if (zoneType == "T" || zoneType == "TRANSLATE" ||
                 zoneType == "TRANSLATING")
        {
            std::string velocityStr;
            err = zonesElement->QueryStringAttribute("VELOCITY", &velocityStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read direction.");
            std::vector<NekDouble> velocity;
            ParseUtils::GenerateVector(velocityStr, velocity);

            zone = ZoneTranslateShPtr(
                MemoryManager<ZoneTranslate>::AllocateSharedPtr(
                    indx, domain, coordDim, velocity));

            m_moveFlag = true;
        }
        else if (zoneType == "F" || zoneType == "FIXED")
        {
            zone = ZoneFixedShPtr(MemoryManager<ZoneFixed>::AllocateSharedPtr(
                indx, domain, coordDim));
        }
        else if (zoneType == "P" || zoneType == "PRESCRIBED")
        {
            std::string xDeformStr;
            err = zonesElement->QueryStringAttribute("XDEFORM", &xDeformStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read x deform equation.");
            LibUtilities::EquationSharedPtr xDeformEqn =
                std::make_shared<LibUtilities::Equation>(
                    m_session->GetInterpreter(), xDeformStr);

            std::string yDeformStr;
            err = zonesElement->QueryStringAttribute("YDEFORM", &yDeformStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read y deform equation.");
            LibUtilities::EquationSharedPtr yDeformEqn =
                std::make_shared<LibUtilities::Equation>(
                    m_session->GetInterpreter(), yDeformStr);

            std::string zDeformStr;
            err = zonesElement->QueryStringAttribute("ZDEFORM", &zDeformStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read z deform equation.");
            LibUtilities::EquationSharedPtr zDeformEqn =
                std::make_shared<LibUtilities::Equation>(
                    m_session->GetInterpreter(), zDeformStr);

            zone = ZonePrescribeShPtr(
                MemoryManager<ZonePrescribe>::AllocateSharedPtr(
                    indx, domain, coordDim, xDeformEqn, yDeformEqn, zDeformEqn));

            m_moveFlag = true;
        }
        else
        {
            WARNINGL0(false, "Zone type '" + zoneType +
                                 "' is unsupported. Valid types are: 'Fixed', 'Rotate', 'Translate', or 'Prescribe'.")
        }

        m_zones[indx] = zone;

        zonesElement = zonesElement->NextSiblingElement();
    }
}

void Movement::ReadInterfaces(TiXmlElement *interfacesTag)
{
    ASSERTL0(interfacesTag, "Unable to find INTERFACES tag in file.");

    TiXmlElement *interfaceElement = interfacesTag->FirstChildElement();

    while (interfaceElement)
    {
        ASSERTL0(
            "INTERFACE" == (std::string)interfaceElement->Value(),
            "Only INTERFACE tags may be present inside the INTERFACES block.")

        int err;

        std::string name;
        err = interfaceElement->QueryStringAttribute("NAME", &name);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read interface name.");
        TiXmlElement *sideElement = interfaceElement->FirstChildElement();

        int cnt = 0;
        Array<OneD, InterfaceShPtr> interfaces(2);
        while (sideElement)
        {
            ASSERTL0(cnt < 2,
                     "Only two sides may be present in each interface block.")
            std::string sideStr = sideElement->Value();

            int indx; // @TODO: Do I need ID?
            err = sideElement->QueryIntAttribute("ID", &indx);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read interface ID.");

            std::string boundaryStr;
            int boundaryErr =
                sideElement->QueryStringAttribute("BOUNDARY", &boundaryStr);

            CompositeMap boundaryEdge;
            if (boundaryErr == TIXML_SUCCESS)
            {
                std::string indxStr = ReadTag(boundaryStr);
                m_meshGraph->GetCompositeList(indxStr, boundaryEdge);
            }

            interfaces[cnt++] =
                InterfaceShPtr(MemoryManager<Interface>::AllocateSharedPtr(
                    indx, boundaryEdge));

            sideElement = sideElement->NextSiblingElement();
        }

        m_interfaces[std::make_pair(m_interfaces.size(), name)] =
            InterfacePairShPtr(MemoryManager<InterfacePair>::AllocateSharedPtr(
                interfaces[0], interfaces[1]));
        interfaceElement = interfaceElement->NextSiblingElement();
    }
}

void Movement::PerformMovement(NekDouble timeStep)
{
    std::set<int> movedZoneIds;
    for (auto &zone : m_zones)
    {
        if (zone.second->Move(timeStep))
        {
            movedZoneIds.insert(zone.first);
        }
    }

    // If zone has moved, set all interfaces on that zone to moved.
    // @TODO: Probably better to save the moved flag on the interface pair obj?
    for (auto &interPair : m_interfaces)
    {
        int leftId  = interPair.second->GetLeftInterface()->GetId();
        int rightId = interPair.second->GetRightInterface()->GetId();

        if (movedZoneIds.find(leftId) != movedZoneIds.end() ||
            movedZoneIds.find(rightId) != movedZoneIds.end())
        {
            m_zones[leftId]->GetMoved()  = true;
            m_zones[rightId]->GetMoved() = true;
        }
    }
}

}
}