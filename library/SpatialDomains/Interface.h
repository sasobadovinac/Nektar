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

#include <map>
#include <string>

#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

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

struct Composite;
typedef std::map<int, std::shared_ptr<Composite>> CompositeMap;

struct InterfaceBase
{
    InterfaceBase(InterfaceType type, CompositeMap domain)
        : m_type(type), m_domain(domain)
    {
    }

    InterfaceType GetInterfaceType() const
    {
        return m_type;
    }

    CompositeMap GetDomain() const
    {
        return m_domain;
    }

    CompositeMap GetInterfaceEdge() const
    {
        return m_interfaceEdge;
    }

    void SetInterfaceEdge(const CompositeMap &interfaceEdge)
    {
        m_interfaceEdge = interfaceEdge;
    }

    std::map<int, SegGeomSharedPtr> const &GetEdge() const
    {
        return m_edge;
    }

    std::vector<int> const &GetEdgeIds() const
    {
        return m_edgeIds;
    }

    void SetEdge(const SegGeomSharedPtr &edge)
    {
        m_edge[edge->GetGlobalID()] = edge;
    }

    void SetEdge(const CompositeMap &edge);

protected:
    InterfaceType m_type;
    CompositeMap m_domain;
    std::map<int, SegGeomSharedPtr> m_edge;
    std::vector<int> m_edgeIds;
    CompositeMap m_interfaceEdge;
};

struct RotatingInterface : public InterfaceBase
{
    RotatingInterface(const CompositeMap &domain, const PointGeom &origin,
                      const std::vector<NekDouble> &axis,
                      const NekDouble angularVel)
        : InterfaceBase(eRotating, domain), m_origin(origin),
          m_axis(axis), m_angularVel(angularVel)
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
    PointGeom m_origin;
    std::vector<NekDouble> m_axis;
    NekDouble m_angularVel;
};

struct FixedInterface : public InterfaceBase
{
    FixedInterface(const CompositeMap &domain)
            : InterfaceBase(eFixed, domain)
    {
    }
};

typedef std::shared_ptr<InterfaceBase> InterfaceBaseShPtr;
typedef std::shared_ptr<RotatingInterface> RotatingInterfaceShPtr;
typedef std::shared_ptr<FixedInterface> FixedInterfaceShPtr;

struct InterfacePair
{
    InterfacePair(InterfaceBaseShPtr leftInterface,
                  InterfaceBaseShPtr rightInterface)
                 : m_leftInterface(leftInterface),
                   m_rightInterface(rightInterface)
    {
    }

    InterfaceBaseShPtr m_leftInterface;
    InterfaceBaseShPtr m_rightInterface;
    bool m_calcFlag = true;

public:
    const InterfaceBaseShPtr &GetLeftInterface() const
    {
        return m_leftInterface;
    }

    const InterfaceBaseShPtr &GetRightInterface() const
    {
        return m_rightInterface;
    }

    bool GetCalcFlag()
    {
        return m_calcFlag;
    }

    void SetCalcFlag(bool flag)
    {
        m_calcFlag = flag;
    }

    void SeparateGraph(MeshGraphSharedPtr &graph);
};


typedef std::shared_ptr<InterfacePair> InterfacePairShPtr;
typedef std::map<int, InterfacePairShPtr> InterfaceCollection;
typedef std::shared_ptr<InterfaceCollection> InterfaceCollectionShPtr;

class Interfaces
{
public:
    SPATIAL_DOMAINS_EXPORT Interfaces(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const MeshGraphSharedPtr &meshGraph);

    SPATIAL_DOMAINS_EXPORT Interfaces() = default;

    const InterfaceCollection &GetInterfaces() const
    {
        return m_interfaces;
    }

protected:
    /// The mesh graph to use for referencing geometry info.
    MeshGraphSharedPtr m_meshGraph;
    LibUtilities::SessionReaderSharedPtr m_session;
    InterfaceCollection m_interfaces;

private:
    /// Read segments (and general MeshGraph) given TiXmlDocument.
    void Read(TiXmlElement *interfaceTag);
    void ReadInterfaces(TiXmlElement *interfaceTag);
};

typedef std::shared_ptr<Interfaces> InterfacesSharedPtr;

} // namespace SpatialDomains
} // namespace Nektar

#endif // NEKTAR_SPATIALDOMAINS_INTERFACES_H
