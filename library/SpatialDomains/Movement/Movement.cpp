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
//  Description: This file contains the base class for implementing
//               non-conformal geometry using the Movement object
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/Movement/Movement.h>
#include <tinyxml.h>

namespace Nektar
{
namespace SpatialDomains
{

std::string static inline ReadTag(std::string &tagStr)
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

Movement::Movement(const LibUtilities::SessionReaderSharedPtr &pSession,
                   MeshGraph* meshGraph)
{
    TiXmlNode *conditions = pSession->GetElement("NEKTAR")->FirstChild("CONDITIONS");
    if(conditions == nullptr)
    {
        return;
    }

    TiXmlNode *movement = conditions->FirstChild("MOVEMENT");
    if (movement != nullptr)
    {
        bool zones = movement->FirstChild("ZONES") != nullptr;
        if (zones)
        {
            ReadZones(pSession->GetElement(
                "NEKTAR/CONDITIONS/MOVEMENT/ZONES"), meshGraph, pSession);
        }

        bool interfaces = movement->FirstChild("INTERFACES") != nullptr;
        if (interfaces)
        {
            ReadInterfaces(pSession->GetElement(
                "NEKTAR/CONDITIONS/MOVEMENT/INTERFACES"), meshGraph);
        }

        ASSERTL0(zones == interfaces,
                 "Only one of ZONES or INTERFACES present in the MOVEMENT "
                 "block.")
    }

    // DEBUG COMMENTS
    if (movement != nullptr && pSession->DefinesCmdLineArgument("verbose"))
    {
        auto comm = pSession->GetComm();
        if (comm->TreatAsRankZero())
        {
            std::cout << "Movement Info:\n";
            std::cout << "\tNum zones: " << m_zones.size() << "\n";
        }

        for (auto &zone : m_zones)
        {
            auto numEl = zone.second->GetElements().size();
            comm->AllReduce(numEl, LibUtilities::ReduceSum);

            // Find shape type if not on this proc
            int shapeType = !zone.second->GetElements().empty()
                    ? zone.second->GetElements().front()->GetShapeType() : -1;
            comm->AllReduce(shapeType, LibUtilities::ReduceMax);

            if (comm->TreatAsRankZero())
            {
                std::cout
                    << "\t- " << zone.first << " "
                    << MovementTypeStr[static_cast<int>(
                           zone.second->GetMovementType())]
                    << ": " << numEl << " "
                    << LibUtilities::ShapeTypeMap
                           [shapeType]
                    << "s\n";
            }
        }

        if (comm->TreatAsRankZero())
        {
            std::cout << "\tNum interfaces: " << m_interfaces.size() << "\n";
        }

        for (auto &interface : m_interfaces)
        {
            auto numLeft =
                interface.second->GetLeftInterface()->GetEdge().size();
            auto numRight =
                interface.second->GetRightInterface()->GetEdge().size();
            comm->AllReduce(numLeft, LibUtilities::ReduceSum);
            comm->AllReduce(numRight, LibUtilities::ReduceSum);

            // Find shape type if not on this proc
            int shapeTypeLeft = ! interface.second->GetLeftInterface()
                                        ->GetEdge().empty()
                                ? interface.second->GetLeftInterface()
                                          ->GetEdge().begin()
                                          ->second->GetShapeType() : -1;
            comm->AllReduce(shapeTypeLeft, LibUtilities::ReduceMax);
            int shapeTypeRight = ! interface.second->GetRightInterface()
                                        ->GetEdge().empty()
                                    ? interface.second->GetRightInterface()
                                          ->GetEdge().begin()
                                          ->second->GetShapeType() : -1;
            comm->AllReduce(shapeTypeRight, LibUtilities::ReduceMax);

            if (comm->TreatAsRankZero())
            {
                std::cout
                    << "\t- \"" << interface.first.second << "\": "
                    << interface.second->GetLeftInterface()->GetId() << " ("
                    << numLeft << " "
                    << LibUtilities::ShapeTypeMap[shapeTypeLeft]
                    << "s) <-> "
                    << interface.second->GetRightInterface()->GetId() << " ("
                    << numRight << " "
                    << LibUtilities::ShapeTypeMap[shapeTypeRight]
                    << "s)\n";
            }
        }

        comm->Block();
        if (comm->TreatAsRankZero())
        {
            std::cout << std::endl;
        }
    }
}

void Movement::ReadZones(TiXmlElement *zonesTag,
                         MeshGraph* meshGraph,
                         const LibUtilities::SessionReaderSharedPtr &pSession)
{
    int coordDim = meshGraph->GetSpaceDimension();

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

        auto &domains = meshGraph->GetDomain();
        auto domFind  = stoi(ReadTag(interfaceDomainStr));
        std::map<int, CompositeSharedPtr> domain;
        if (domains.find(domFind) != domains.end())
        {
            domain = domains.at(domFind);
        }

        ZoneBaseShPtr zone;
        if (zoneType == "F" ||
            zoneType == "FIXED")
        {
            zone = ZoneFixedShPtr(MemoryManager<ZoneFixed>::AllocateSharedPtr(
                indx, domain, coordDim));
        }
        else if (zoneType == "R" ||
                 zoneType == "ROTATE" ||
                 zoneType == "ROTATING")
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

            pSession->SubstituteExpressions(angularVelStr);
            LibUtilities::EquationSharedPtr angularVelEqn =
                MemoryManager<LibUtilities::Equation>::AllocateSharedPtr(
                    pSession->GetInterpreter(), angularVelStr);

            zone = ZoneRotateShPtr(MemoryManager<ZoneRotate>::AllocateSharedPtr(
                indx, domain, coordDim, origin, axis, angularVelEqn));

            m_moveFlag = true;
        }
        else if (zoneType == "T" ||
                 zoneType == "TRANSLATE" ||
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
        else if (zoneType == "P" ||
                 zoneType == "PRESCRIBED")
        {
            std::string xDeformStr;
            err = zonesElement->QueryStringAttribute("XDEFORM", &xDeformStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read x deform equation.");
            LibUtilities::EquationSharedPtr xDeformEqn =
                std::make_shared<LibUtilities::Equation>(
                    pSession->GetInterpreter(), xDeformStr);

            std::string yDeformStr;
            err = zonesElement->QueryStringAttribute("YDEFORM", &yDeformStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read y deform equation.");
            LibUtilities::EquationSharedPtr yDeformEqn =
                std::make_shared<LibUtilities::Equation>(
                    pSession->GetInterpreter(), yDeformStr);

            std::string zDeformStr;
            err = zonesElement->QueryStringAttribute("ZDEFORM", &zDeformStr);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read z deform equation.");
            LibUtilities::EquationSharedPtr zDeformEqn =
                std::make_shared<LibUtilities::Equation>(
                    pSession->GetInterpreter(), zDeformStr);

            zone = ZonePrescribeShPtr(
                MemoryManager<ZonePrescribe>::AllocateSharedPtr(
                    indx, domain, coordDim, xDeformEqn, yDeformEqn, zDeformEqn));

            m_moveFlag = true;
        }
        else
        {
            WARNINGL0(false, "Zone type '" + zoneType +
                                 "' is unsupported. Valid types are: 'Fixed', "
                                 "'Rotate', 'Translate', or 'Prescribe'.")
        }

        m_zones[indx] = zone;
        zonesElement = zonesElement->NextSiblingElement();
    }
}

void Movement::ReadInterfaces(TiXmlElement *interfacesTag,
                              MeshGraph* meshGraph)
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

        // @TODO: For different interface types have a string attribute type in
        // @TODO: the INTERFACE element like for NAME above

        Array<OneD, InterfaceShPtr> interfaces(2);
        std::vector<int> cnt;
        while (sideElement)
        {
            ASSERTL0(cnt.size() < 2,
                "In INTERFACE NAME " + name +
                    ", only two sides may be present in each INTERFACE block.")

            int indx;
            err = sideElement->QueryIntAttribute("ID", &indx);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read interface ID.");

            std::string boundaryStr;
            int boundaryErr =
                sideElement->QueryStringAttribute("BOUNDARY", &boundaryStr);

            CompositeMap boundaryEdge;
            if (boundaryErr == TIXML_SUCCESS)
            {
                std::string indxStr = ReadTag(boundaryStr);
                meshGraph->GetCompositeList(indxStr, boundaryEdge);
            }

            // Sets location in interface pair to 0 for left, and 1 for right
            auto sideElVal = sideElement->ValueStr();
            if (sideElVal == "LEFT" || sideElVal == "L")
            {
                cnt.emplace_back(0);
            }
            else if (sideElVal == "RIGHT" || sideElVal == "R")
            {
                cnt.emplace_back(1);
            }
            else
            {
                NEKERROR(ErrorUtil::efatal, sideElement->ValueStr() +
                             " is not a valid interface side for interface "
                             "NAME " + name + ". Please only use LEFT or RIGHT.")
            }

            interfaces[cnt[cnt.size()-1]] =
                InterfaceShPtr(MemoryManager<Interface>::AllocateSharedPtr(
                    indx, boundaryEdge));

            sideElement = sideElement->NextSiblingElement();
        }

        ASSERTL0(std::accumulate(cnt.begin(), cnt.end(), 0) == 1,
                 "You must have only one LEFT and one RIGHT side"
                 " present in interface NAME " + name)

        m_interfaces[std::make_pair(m_interfaces.size(), name)] =
            InterfacePairShPtr(MemoryManager<InterfacePair>::AllocateSharedPtr(
                interfaces[0], interfaces[1]));
        interfaceElement = interfaceElement->NextSiblingElement();
    }
}

// Acts as a placeholder for when ALE function and moving geometry capability
// is added. Currently unused.
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