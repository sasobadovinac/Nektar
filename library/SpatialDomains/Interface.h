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
                    CompositeMap rightDomain,
                    CompositeMap leftDomain):
                    m_interfaceType(type),
                    m_rightDomain(rightDomain),
                    m_leftDomain(leftDomain)
            {
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

            void SetInterfaceEdge(const CompositeMap &interfaceEdge)
            {
                m_interfaceEdge = interfaceEdge;
            }

            std::map<int, SegGeomSharedPtr> const &GetEdgeLeft() const
            {
                return m_leftEdge;
            }

            void SetEdgeLeft(const std::map<int, SegGeomSharedPtr> &leftEdge)
            {
                m_leftEdge = leftEdge;
            }

            void SetEdgeLeft(const CompositeMap &leftEdge)
            {
                for (auto &compIt : leftEdge)
                {
                    for (auto &elmtIt : compIt.second->m_geomVec)
                    {
                        SegGeomSharedPtr elmt = std::dynamic_pointer_cast<SegGeom>(elmtIt);
                        ASSERTL0(elmt, "Composite for left edge should only contain segments");
                        m_leftEdge[elmt->GetGlobalID()] = elmt;
                    }
                }
            }

            std::map<int, SegGeomSharedPtr> const &GetEdgeRight() const
            {
                return m_rightEdge;
            }

            void SetEdgeRight(const std::map<int, SegGeomSharedPtr> &rightEdge)
            {
                m_rightEdge = rightEdge;
            }

            void SetEdgeRight(const CompositeMap &rightEdge)
            {
                for (auto &compIt : rightEdge)
                {
                    for (auto &elmtIt : compIt.second->m_geomVec)
                    {
                        SegGeomSharedPtr elmt = std::dynamic_pointer_cast<SegGeom>(elmtIt);
                        ASSERTL0(elmt, "Composite for right edge should only contain segments");
                        m_rightEdge[elmt->GetGlobalID()] = elmt;
                    }
                }
            }

            void SeparateGraph(MeshGraphSharedPtr &graph);

        protected:
            InterfaceType      m_interfaceType;
            CompositeMap       m_rightDomain;
            CompositeMap       m_leftDomain;
            std::map<int, SegGeomSharedPtr>   m_leftEdge;
            std::map<int, SegGeomSharedPtr>   m_rightEdge;
            CompositeMap       m_interfaceEdge;

        };

        struct RotatingInterface : public InterfaceBase
        {
            RotatingInterface(
                    const CompositeMap rightDomain,
                    const CompositeMap leftDomain,
                    const PointGeom origin,
                    const std::vector<NekDouble> axis,
                    const NekDouble angularVel)
                    : InterfaceBase(eRotating, rightDomain, leftDomain),
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
