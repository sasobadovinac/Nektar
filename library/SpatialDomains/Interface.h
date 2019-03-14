////////////////////////////////////////////////////////////////////////////////
//
//  File:  s.h
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
//  Software is furnished to do so, subject to the following s:
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
#ifndef NEKTAR_SPATIALDOMAINS_INTERFACES_H
#define NEKTAR_SPATIALDOMAINS_INTERFACES_H

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
        enum InterfaceType
        {
            eFixed,
            eRotating,
            eSliding,
            eNotDefined
        };

        typedef std::map<int, CompositeSharedPtr> InterfaceEdge;
        typedef std::shared_ptr<InterfaceEdge> InterfaceEdgeShPtr;

        struct InterfaceBase
        {
            InterfaceBase(
                    InterfaceType type,
                    CompositeMap movingDomain,
                    CompositeMap fixedDomain,
                    InterfaceEdgeShPtr interfaceEdge):
                    m_interfaceType(type),
                    m_movingDomain(movingDomain),
                    m_fixedDomain(fixedDomain),
                    m_interfaceEdge(interfaceEdge)

            {
            }

            virtual ~InterfaceBase()
            {};

            InterfaceType GetInterfaceType() const
            {
                return m_interfaceType;
            }

            CompositeMap GetMovingDomain() const
            {
                return m_movingDomain;
            }

            CompositeMap GetFixedDomain() const
            {
                return m_fixedDomain;
            }

            InterfaceEdgeShPtr GetInterfaceEdge() const
            {
                return m_interfaceEdge;
            }

        protected:
            InterfaceType      m_interfaceType;
            CompositeMap       m_movingDomain;
            CompositeMap       m_fixedDomain;
            InterfaceEdgeShPtr m_interfaceEdge;
        };

        struct RotatingInterface : public InterfaceBase
        {
            RotatingInterface(
                    const LibUtilities::SessionReaderSharedPtr &pSession,
                    const CompositeMap movingDomain,
                    const CompositeMap fixedDomain,
                    const InterfaceEdgeShPtr interfaceEdge,
                    const PointGeomSharedPtr origin,
                    const std::vector<NekDouble> axis,
                    const NekDouble angularVel)
                    : InterfaceBase(eRotating, movingDomain, fixedDomain, interfaceEdge),
                      m_origin(origin),
                      m_axis(axis),
                      m_angularVel(angularVel)
            {
            }

            PointGeomSharedPtr GetOrigin() const
            {
                return m_origin;
            }

            std::vector<NekDouble> GetAxis() const
            {
                return m_axis;
            }

            NekDouble GetAngularVel() const
            {
                return m_angularVel;
            }

            PointGeomSharedPtr      m_origin;
            std::vector<NekDouble>  m_axis;
            NekDouble               m_angularVel;
        };

        typedef std::shared_ptr<InterfaceBase> InterfaceShPtr;
        typedef std::shared_ptr<RotatingInterface> RotatingInterfaceShPtr;

        typedef std::map<int, InterfaceShPtr> InterfaceCollection;



        class Interfaces
        {
        public:
            SPATIAL_DOMAINS_EXPORT Interfaces(
                    const LibUtilities::SessionReaderSharedPtr &pSession,
                    const MeshGraphSharedPtr &meshGraph);

            SPATIAL_DOMAINS_EXPORT Interfaces(void);
            SPATIAL_DOMAINS_EXPORT ~Interfaces(void);

            const InterfaceCollection &GetInterfaces(void) const
            {
                return m_interfaces;
            }


        protected:
            /// The mesh graph to use for referencing geometry info.
            MeshGraphSharedPtr                      m_meshGraph;
            LibUtilities::SessionReaderSharedPtr    m_session;

            InterfaceCollection                     m_interfaces;
        private:

            /// Read segments (and general MeshGraph) given TiXmlDocument.
            void Read(TiXmlElement *interfaceTag);
            void ReadInterfaces(TiXmlElement *interfaceTag);
        };
    }
}

#endif //NEKTAR_SPATIALDOMAINS_INTERFACES_H
