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

enum InterfaceSide
{
    eNone,
    eLeft,
    eRight
};

const std::string InterfaceSideStr[] =
{
    "NONE",
    "LEFT",
    "RIGHT"
};

struct Composite;
typedef std::map<int, std::shared_ptr<Composite>> CompositeMap;

struct InterfaceBase;
typedef std::shared_ptr<InterfaceBase> InterfaceBaseShPtr;

struct InterfaceBase
{
    InterfaceBase(InterfaceType type, int indx, CompositeMap domain)
        : m_type(type), m_id(indx), m_domain(domain)
    {
    }

    inline InterfaceType GetInterfaceType() const
    {
        return m_type;
    }

    inline CompositeMap GetDomain() const
    {
        return m_domain;
    }

    inline CompositeMap GetInterfaceEdge() const
    {
        return m_interfaceEdge;
    }

    inline void SetInterfaceEdge(const CompositeMap &interfaceEdge)
    {
        m_interfaceEdge = interfaceEdge;
    }

    inline std::map<int, SegGeomSharedPtr> const &GetEdge() const
    {
        return m_edge;
    }

    inline std::vector<int> const &GetEdgeIds() const
    {
        return m_edgeIds;
    }

    inline bool IsEmpty() const
    {
        return m_edge.empty();
    }

    inline void SetEdge(const SegGeomSharedPtr &edge)
    {
        m_edge[edge->GetGlobalID()] = edge;
    }

    void SetEdge(const CompositeMap &edge);

    inline InterfaceSide GetSide() const
    {
        return m_side;
    }

    inline void SetSide(const InterfaceSide &side)
    {
        m_side = side;
    }

    inline void SetTotPoints(const int &points)
    {
        m_totPoints = points;
    }

    inline int &GetTotPoints()
    {
        return m_totPoints;
    }

    inline void AddOppRank(const int &rank)
    {
        m_oppRanks.emplace_back(rank);
    }

    inline std::vector<int> &GetOppRanks()
    {
        return m_oppRanks;
    }

    inline void SetOppInterface(InterfaceBaseShPtr oppInterface)
    {
        m_oppInterface = oppInterface;
    }

    inline InterfaceBaseShPtr GetOppInterface()
    {
        return m_oppInterface;
    }

    inline int &GetId()
    {
        return m_id;
    }

    inline std::map<int, std::tuple<NekDouble, NekDouble, NekDouble>> GetMissing()
    {
        return m_missingCoords;
    }

    inline void SetMissing(std::map<int, std::tuple<NekDouble, NekDouble, NekDouble>> missingCoords)
    {
        m_missingCoords = missingCoords;
    }

protected:
    InterfaceBaseShPtr m_oppInterface;
    InterfaceType m_type;
    int m_id;
    InterfaceSide m_side = eNone;
    CompositeMap m_domain;
    std::map<int, SegGeomSharedPtr> m_edge;
    std::vector<int> m_edgeIds;
    CompositeMap m_interfaceEdge;
    int m_totPoints;
    std::vector<int> m_oppRanks;
    bool m_checkLocal = false;
    std::map<int, std::tuple<NekDouble, NekDouble, NekDouble>> m_missingCoords;
};

struct RotatingInterface : public InterfaceBase
{
    RotatingInterface(int id, const CompositeMap &domain, const PointGeom &origin,
                      const std::vector<NekDouble> &axis,
                      const NekDouble angularVel)
        : InterfaceBase(eRotating, id, domain), m_origin(origin),
          m_axis(axis), m_angularVel(angularVel)
    {
    }

    inline PointGeom GetOrigin() const
    {
        return m_origin;
    }

    inline std::vector<NekDouble> GetAxis() const
    {
        return m_axis;
    }

    inline NekDouble GetAngularVel() const
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
    FixedInterface(int id, const CompositeMap &domain)
            : InterfaceBase(eFixed, id, domain)
    {
    }
};

typedef std::shared_ptr<RotatingInterface> RotatingInterfaceShPtr;
typedef std::shared_ptr<FixedInterface> FixedInterfaceShPtr;

struct InterfacePair
{
    InterfacePair(InterfaceBaseShPtr leftInterface,
                  InterfaceBaseShPtr rightInterface)
                 : m_leftInterface(leftInterface),
                   m_rightInterface(rightInterface)
    {
        leftInterface->SetSide(eLeft);
        rightInterface->SetSide(eRight);
    }

    InterfaceBaseShPtr m_leftInterface;
    InterfaceBaseShPtr m_rightInterface;
    bool m_calcFlag = true;
    bool m_checkLocal = false;

public:
    inline const InterfaceBaseShPtr &GetLeftInterface() const
    {
        return m_leftInterface;
    }

    inline const InterfaceBaseShPtr &GetRightInterface() const
    {
        return m_rightInterface;
    }

    inline bool GetCalcFlag()
    {
        return m_calcFlag;
    }

    inline void SetCalcFlag(bool flag)
    {
        m_calcFlag = flag;
    }
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

    inline const InterfaceCollection &GetInterfaces() const
    {
        return m_interfaces;
    }

    void CommSetup(LibUtilities::CommSharedPtr &comm);

protected:
    /// The mesh graph to use for referencing geometry info.
    MeshGraphSharedPtr m_meshGraph;
    LibUtilities::SessionReaderSharedPtr m_session;
    InterfaceCollection m_interfaces;
    /// Map of interface ID to ranks sharing that interface (on the other side)
    std::map<int, std::vector<int>> m_interfacesSharedRank;


private:
    /// Read segments (and general MeshGraph) given TiXmlDocument.
    void Read(TiXmlElement *interfaceTag);
    void ReadInterfaces(TiXmlElement *interfaceTag);
};

typedef std::shared_ptr<Interfaces> InterfacesSharedPtr;

} // namespace SpatialDomains
} // namespace Nektar

#endif // NEKTAR_SPATIALDOMAINS_INTERFACES_H
