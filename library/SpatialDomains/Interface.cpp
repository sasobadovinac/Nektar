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
        /**
         * Constructor - collective on the session's communicator.
         */
        InterfaceConditions::InterfaceConditions(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const MeshGraphSharedPtr &meshGraph) :
                m_meshGraph(meshGraph), m_session(pSession)
        {
            Read(m_session->GetElement("Nektar/Conditions"));
        }

        InterfaceConditions::InterfaceConditions(void)
        {
        }

        InterfaceConditions::~InterfaceConditions(void)
        {
        }

        /**
         * Create a new communicator for each boundary region.
         * Collective on the session's communicator.
         */
        //void InterfaceConditions::CreateInterfaceComms()
        /*{
            LibUtilities::CommSharedPtr comm = m_session->GetComm();

            std::set<int> allids = ShareAllBoundaryIDs(m_boundaryRegions, comm);

            for (auto &it : allids)
            {
                auto reg_it = m_boundaryRegions.find(it);
                int this_rank_participates = (reg_it != m_boundaryRegions.end());
                LibUtilities::CommSharedPtr comm_region = comm->CommCreateIf(
                        this_rank_participates);

                ASSERTL0(bool(comm_region) == bool(this_rank_participates),
                         "Rank should be in communicator but wasn't or is in "
                         "communicator but shouldn't be.");

                if (this_rank_participates)
                {
                    m_boundaryCommunicators[reg_it->first] = comm_region;
                }
            }
        }*/

        /**
         * Collective on the session's communicator.
         */
        void InterfaceConditions::Read(TiXmlElement *conditions)
        {
            ASSERTL0(conditions, "Unable to find CONDITIONS tag in file.");

            TiXmlElement *interfaceRegions = conditions->FirstChildElement(
                    "INTERFACEREGIONS");

            if (interfaceRegions)
            {
                ReadInterfaceRegions(conditions);
                CreateInterfaceComms();
                ReadInterfaceConditions(conditions);
            }
        }

        /**
         *
         */
        void InterfaceConditions::ReadInterfaceRegions(TiXmlElement *conditions)
        {
            // ensure interface regions only read once per class definition
            if (m_interfaceRegions.size() != 0)
            {
                return;
            }

            TiXmlElement *interfaceRegions = conditions->FirstChildElement(
                    "INTERFACEREGIONS");
            ASSERTL0(interfaceRegions,
                     "Unable to find INTERFACEREGIONS block.");

            // See if we have boundary regions defined.
            TiXmlElement *interfaceRegionsElement =
                    interfaceRegions->FirstChildElement("I");

            while (interfaceRegionsElement)
            {
                /// All elements are of the form: "<I ID="#"> ... </B>", with
                /// ? being the element type.
                int indx;
                int err = interfaceRegionsElement->QueryIntAttribute("ID",
                                                                     &indx);
                ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");

                TiXmlNode *interfaceRegionChild =
                        interfaceRegionsElement->FirstChild();
                // This is primarily to skip comments that may be present.
                // Comments appear as nodes just like elements.
                // We are specifically looking for text in the body
                // of the definition.
                while (interfaceRegionChild
                       && interfaceRegionChild->Type()
                          != TiXmlNode::TINYXML_TEXT)
                {
                    interfaceRegionChild = interfaceRegionChild->NextSibling();
                }

                ASSERTL0(interfaceRegionChild,
                         "Unable to read variable definition body.");
                std::string interfaceRegionStr =
                        interfaceRegionChild->ToText()->ValueStr();

                std::string::size_type indxBeg =
                        interfaceRegionStr.find_first_of('[') + 1;
                std::string::size_type indxEnd =
                        interfaceRegionStr.find_last_of(
                                ']') - 1;

                ASSERTL0(indxBeg <= indxEnd,
                         (std::string(
                                 "Error reading interface region definition:")
                          + interfaceRegionStr).c_str());

                std::string indxStr = interfaceRegionStr.substr(indxBeg,
                                                                indxEnd -
                                                                indxBeg + 1);

                if (!indxStr.empty())
                {
                    // Extract the composites from the string and return them in a list.
                    InterfaceRegionShPtr interfaceRegion(
                            MemoryManager<InterfaceRegion>::AllocateSharedPtr());

                    ASSERTL0(m_interfaceRegions.count(indx) == 0,
                             "Interface region " + indxStr +
                             " defined more than "
                             "once!");

                    m_meshGraph->GetCompositeList(indxStr, *interfaceRegion);
                    m_interfaceRegions[indx] = interfaceRegion;
                }

                interfaceRegionsElement =
                        interfaceRegionsElement->NextSiblingElement("I");
            }
        }

        /**
        *
        */
        void InterfaceConditions::ReadInterfaceConditions(
                TiXmlElement *conditions)
        {
            // Protect against multiple reads.
            if (!m_interfaceConditions.empty()) //before was .size()!=0;
            {
                return;
            }

            // Read REGION tags
            TiXmlElement *interfaceConditionsElement =
                    conditions->FirstChildElement("INTERFACECONDITIONS");
            ASSERTL0(interfaceConditionsElement,
                     "Boundary conditions must be specified.");

            TiXmlElement *regionElement =
                    interfaceConditionsElement->FirstChildElement("REGION");

            // Read F (Fixed), R (Rotating), S (Sliding) tags
            while (regionElement)
            {
                InterfaceConditionMapShPtr interfaceConditions = MemoryManager<
                        InterfaceConditionMap>::AllocateSharedPtr();

                int interfaceRegionID;
                int err = regionElement->QueryIntAttribute("REF",
                                                           &interfaceRegionID);
                ASSERTL0(err == TIXML_SUCCESS,
                         "Error reading interface region reference.");

                ASSERTL0(m_interfaceConditions.count(interfaceRegionID) == 0,
                         "Interface region '" + boost::lexical_cast<std::string
                                                                   >(
                                 interfaceRegionID)
                         + "' appears multiple times.");

                // Find the interface region corresponding to this ID.
                std::string interfaceRegionIDStr;
                std::ostringstream interfaceRegionIDStrm(interfaceRegionIDStr);
                interfaceRegionIDStrm << interfaceRegionID;

                ASSERTL0(m_interfaceRegions.count(interfaceRegionID) == 1,
                         "Interface region " + boost::lexical_cast<string
                                                                  >(
                                 interfaceRegionID) + " not found");

                // Find the communicator that belongs to this ID
                LibUtilities::CommSharedPtr interfaceRegionComm =
                        m_interfaceCommunicators[interfaceRegionID];

                TiXmlElement *conditionElement =
                        regionElement->FirstChildElement();
                //Below should get list of all domains defined in geometry
                //std::vector <std::map<int, CompositeSharedPtr>> domain = m_meshGraph->GetDomain();


                while (conditionElement)
                {
                    // Check type.
                    std::string conditionType = conditionElement->Value();
                    std::string attrData;
                    bool isTimeDependent = false;

                    // All have var specified, or else all variables are zero.
                    TiXmlAttribute *attr = conditionElement->FirstAttribute();

                    std::vector<std::string>::iterator iter;
                    std::string attrName;

                    attrData = conditionElement->Attribute("DOMAIN");

                    //CHANGE THIS: Currently checks all domains are in every interface region
                    //Will be wrong with multiple interfaces in future
                    /*
                    if (!attrData.empty())
                    {
                        iter = std::find(domain.begin(), domain.end(), attrData);
                        ASSERTL0(iter != domain.end(),
                                 (std::string("Cannot find domain: ")
                                  + attrData).c_str());
                    }
                    */

                    if (conditionType == "R") //Rotating domain
                    {
                        if (attrData.empty())
                        {
                            ASSERTL0(!attrData.empty(),
                                     "Rotating region must have associated attributes.");
                        }
                        else
                        {
                            // Use the iterator from above, which must point to the variable.
                            attr = attr->Next();

                            if (attr)
                            {
                                vector<unsigned int> rotatingIntRegionIndex;
                                std::vector<NekDouble> origin, axis;
                                NekDouble angularVel;
                                while (attr)
                                {

                                    attrName = attr->Name();

                                    if (attrName == "ORIGIN")
                                    {
                                        attrData = attr->Value();
                                        ASSERTL0(!attrData.empty(),
                                                 "ORIGIN attribute must have associated value.");

                                        // Suppose to go here?
                                        m_session->SubstituteExpressions(
                                                attrData);

                                        ParseUtils::GenerateVector(attrData,
                                                                   origin);

                                        ASSERTL0(origin.size() == 3,
                                                 "ORIGIN attribute must have length of 3.");

                                    }
                                    else if (attrName == "AXIS")
                                    {
                                        attrData = attr->Value();
                                        ASSERTL0(!attrData.empty(),
                                                 "AXIS attribute must have associated value.");

                                        m_session->SubstituteExpressions(
                                                attrData);

                                        ParseUtils::GenerateVector(attrData,
                                                                   axis);

                                        ASSERTL0(axis.size() == 3,
                                                 "AXIS attribute must have length of 3.");

                                    }
                                    else if (attrName == "ANGVEL")
                                    {
                                        attrData = attr->Value();
                                        ASSERTL0(!attrData.empty(),
                                                 "ANGVEL attribute must be specified.");

                                        m_session->SubstituteExpressions(
                                                attrData);

                                        //Convert from string to NekDouble
                                        //@todo add PI parsing so expressions such as 2*PI can be understood
                                        angularVel = static_cast<NekDouble>(stod(
                                                attrData));
                                    }
                                    else
                                    {
                                        ASSERTL0(false,
                                                 (std::string(
                                                         "Unknown rotating interface condition attribute: ") +
                                                  attrName).c_str());
                                    }
                                    attr = attr->Next();
                                }
                                InterfaceConditionShPtr rotatingCondition(
                                        MemoryManager<RotatingInterfaceCondition>::AllocateSharedPtr(
                                                rotatingIntRegionIndex[0],
                                                m_session, origin,
                                                axis, angularVel,
                                                interfaceRegionComm));
                                (*interfaceConditions)[*iter] = rotatingCondition;
                            }
                        }
                    }
                    else if (conditionType == "F")
                    {

                    }
                    else if (conditionType == "S") // Sliding
                    {

                    }

                    conditionElement = conditionElement->NextSiblingElement();

                }

                m_interfaceConditions[interfaceRegionID] = interfaceConditions;
                regionElement = regionElement->NextSiblingElement("REGION");
            }
        }
    }
}
