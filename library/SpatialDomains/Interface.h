////////////////////////////////////////////////////////////////////////////////
//
//  File:  Conditions.h
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
#ifndef NEKTAR_SPATIALDOMAINS_INTERFACECONDITIONS_H
#define NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H

#include <string>
#include <map>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/Equation.h>
#include <SpatialDomains/MeshGraph.h>

namespace Nektar
{
    struct OneD;

    namespace SpatialDomains
    {
        enum InterfaceConditionType
        {
            eFixed,
            eRotating,
            eSliding,
            eNotDefinedInterface
        };

        struct InterfaceConditionBase
        {
            InterfaceConditionBase(
                    InterfaceConditionType type,
                    LibUtilities::CommSharedPtr comm = LibUtilities::CommSharedPtr()):
                    m_interfaceConditionType(type),
                    m_comm(comm)
            {
            }

            virtual ~InterfaceConditionBase()
            {};

            InterfaceConditionType GetInterfaceConditionType() const
            {
                return m_interfaceConditionType;
            }

            void SetInterfaceConditionType(InterfaceConditionType interfaceType)
            {
                m_interfaceConditionType = interfaceType;
            }

            LibUtilities::CommSharedPtr GetComm()
            {
                return m_comm;
            }

        protected:
            InterfaceConditionType      m_interfaceConditionType;
            std::vector<NekDouble>      m_origin;
            std::vector<NekDouble>      m_axis;
            NekDouble                   m_angularVel;
            LibUtilities::CommSharedPtr m_comm;
        };

        struct RotatingInterfaceCondition : public InterfaceConditionBase
        {
            RotatingInterfaceCondition(
                    const unsigned int n,
                    const LibUtilities::SessionReaderSharedPtr &pSession,
                    const std::vector<NekDouble> origin,
                    const std::vector<NekDouble> axis,
                    const NekDouble angularVel,
                    const LibUtilities::CommSharedPtr comm=LibUtilities::CommSharedPtr()):
                    InterfaceConditionBase(eRotating, comm),
                    m_connectedInterfaceRegion(n)
            {
            }

            unsigned int m_connectedInterfaceRegion;

        };

        typedef std::map<int, CompositeSharedPtr> InterfaceRegion;
        typedef std::shared_ptr<InterfaceRegion> InterfaceRegionShPtr;
        typedef std::shared_ptr<const InterfaceRegion> ConstInterfaceRegionShPtr;
        typedef std::map<int, InterfaceRegionShPtr> InterfaceRegionCollection;

        typedef std::shared_ptr<InterfaceConditionBase> InterfaceConditionShPtr;
        typedef std::shared_ptr<RotatingInterfaceCondition> RotatingInterfaceShPtr;

        typedef std::map<std::string,InterfaceConditionShPtr>  InterfaceConditionMap;
        typedef std::shared_ptr<InterfaceConditionMap>  InterfaceConditionMapShPtr;
        typedef std::map<int, InterfaceConditionMapShPtr> InterfaceConditionCollection;

        const static Array<OneD, InterfaceConditionShPtr> NullInterfaceConditionShPtrArray;

        class InterfaceConditions
        {
        public:
            SPATIAL_DOMAINS_EXPORT InterfaceConditions(const LibUtilities::SessionReaderSharedPtr &pSession, const MeshGraphSharedPtr &meshGraph);

            SPATIAL_DOMAINS_EXPORT InterfaceConditions(void);
            SPATIAL_DOMAINS_EXPORT ~InterfaceConditions(void);

            const InterfaceRegionCollection &GetInterfaceRegions(void) const
            {
                return m_interfaceRegions;
            }

            void AddInterfaceRegions(const int regionID, InterfaceRegionShPtr &iRegion)
            {
                m_interfaceRegions[regionID] = iRegion;
            }

            const InterfaceConditionCollection &GetInterfaceConditions(void) const
            {
                return m_interfaceConditions;
            }


            void AddInterfaceConditions(const int regionID, InterfaceConditionMapShPtr &iCond)
            {
                m_interfaceConditions[regionID] = iCond;
            }

            const std::string GetVariable(unsigned int indx)
            {
                return m_session->GetVariable(indx);
            }

        protected:
            /// The mesh graph to use for referencing geometry info.
            MeshGraphSharedPtr                      m_meshGraph;
            LibUtilities::SessionReaderSharedPtr    m_session;

            InterfaceRegionCollection                m_interfaceRegions;
            InterfaceConditionCollection             m_interfaceConditions;
            std::map<int, LibUtilities::CommSharedPtr> m_interfaceCommunicators;

        private:

            /// Read segments (and general MeshGraph) given TiXmlDocument.
            void Read(TiXmlElement *conditions);
            void ReadInterfaceRegions(TiXmlElement *regions);
            void ReadInterfaceConditions(TiXmlElement *conditions);
            void CreateInterfaceComms();
        };
    }
}

#endif //NEKTAR_SPATIALDOMAINS_INTERFACECONDITIONS_H
