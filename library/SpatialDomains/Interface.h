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
            eSliding
        };

        const char* const InterfaceTypeMap[] =
        {
            "Fixed",
            "Rotating",
            "Sliding",
            "NotDefined"
        };

        struct InterfaceBase
        {
            InterfaceBase(
                    InterfaceType type,
                    CompositeMap leftDomain,
                    CompositeMap rightDomain,
                    std::map<int, CompositeSharedPtr> interfaceEdge):
                    m_interfaceType(type),
                    m_leftDomain(leftDomain),
                    m_rightDomain(rightDomain),
                    m_interfaceEdge(interfaceEdge)

            {
            }

            InterfaceBase(
                InterfaceType type,
                CompositeMap leftDomain,
                CompositeMap rightDomain,
                CompositeMap leftEdge,
                CompositeMap rightEdge):
                m_interfaceType(type),
                m_leftDomain(leftDomain),
                m_rightDomain(rightDomain),
                m_leftEdge(leftEdge),
                m_rightEdge(rightEdge)

            {
                for (auto &it : m_leftEdge)
                {
                    m_leftEdgeVector.push_back(it.first);
                }

                for (auto &it : m_rightEdge)
                {
                    m_rightEdgeVector.push_back(it.first);
                }
            }

            InterfaceType GetInterfaceType() const
            {
                return m_interfaceType;
            }

            CompositeMap GetRightDomain() const
            {
                return m_rightDomain;
            }

            CompositeMap GetLeftDomain() const
            {
                return m_leftDomain;
            }

            CompositeMap GetInterfaceEdge() const
            {
                return m_interfaceEdge;
            }

            std::vector<int> const &GetEdgeLeftVector() const
            {
                return m_leftEdgeVector;
            }

            std::vector<int> const &GetEdgeRightVector() const
            {
                return m_rightEdgeVector;
            }

            void SeparateGraph(MeshGraphSharedPtr &graph);

        protected:
            InterfaceType                       m_interfaceType;
            CompositeMap                        m_leftDomain;
            CompositeMap                        m_rightDomain;
            CompositeMap                        m_leftEdge;
            CompositeMap                        m_rightEdge;
            std::vector<int>                    m_leftEdgeVector;
            std::vector<int>                    m_rightEdgeVector;
            CompositeMap                        m_interfaceEdge;

        };

        struct RotatingInterface : public InterfaceBase
        {
            RotatingInterface(
                    const CompositeMap leftDomain,
                    const CompositeMap rightDomain,
                    const std::map<int, CompositeSharedPtr> interfaceEdge,
                    const PointGeom origin,
                    const std::vector<NekDouble> axis,
                    const NekDouble angularVel)
                    : InterfaceBase(eRotating, leftDomain, rightDomain, interfaceEdge),
                      m_origin(origin),
                      m_axis(axis),
                      m_angularVel(angularVel)


            {
            }

            RotatingInterface(
                const CompositeMap leftDomain,
                const CompositeMap rightDomain,
                const std::map<int, CompositeSharedPtr> leftEdge,
                const std::map<int, CompositeSharedPtr> rightEdge,
                const PointGeom origin,
                const std::vector<NekDouble> axis,
                const NekDouble angularVel)
                : InterfaceBase(eRotating, leftDomain, rightDomain, leftEdge, rightEdge),
                  m_origin(origin),
                  m_axis(axis),
                  m_angularVel(angularVel)


            {
            }

            PointGeom GetOrigin() const
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

        protected:
            PointGeom               m_origin;
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

            SPATIAL_DOMAINS_EXPORT Interfaces() = default;

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
