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
#include <deque>

#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

namespace Nektar
{

namespace SpatialDomains
{

enum class MovementType
{
    eNone,
    eFixed,
    eRotate,
    eTranslate,
    ePrescribe
};

const std::string MovementTypeStr[] =
{
    "None",
    "Fixed",
    "Rotating",
    "Sliding",
    "Prescribed"
};

enum class InterfaceSide
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

struct ZoneBase;
typedef std::shared_ptr<ZoneBase> ZoneBaseShPtr;

struct ZoneBase
{
    ZoneBase(MovementType type, int indx, CompositeMap domain, int coordDim);

    virtual ~ZoneBase() = default;

    inline MovementType GetMovementType() const
    {
        return m_type;
    }

    inline CompositeMap GetDomain() const
    {
        return m_domain;
    }

    inline int &GetId()
    {
        return m_id;
    }

    inline virtual bool v_Move(NekDouble time)
    {
        boost::ignore_unused(time);
        return false;
    }

    inline bool Move(NekDouble time)
    {
        return v_Move(time);
    }

    inline std::vector<int> const &GetElementIds() const
    {
        return m_elementIds;
    }

    inline std::vector<GeometrySharedPtr> const &GetElements() const
    {
        return m_elements;
    }

    inline bool &GetMoved()
    {
        return m_moved;
    }

    inline Array<OneD,std::set<GeometrySharedPtr>> &GetConstituentElements()
    {
        return m_constituentElements;
    }

protected:
    //ZoneBaseShPtr m_oppInterface;
    MovementType m_type = MovementType::eNone;
    int m_id;
    //InterfaceSide m_side = eNone;
    CompositeMap m_domain;
    //std::map<int, GeometrySharedPtr> m_edge;
    //std::vector<int> m_edgeIds;
    std::vector<int> m_elementIds;
    std::vector<GeometrySharedPtr> m_elements;
    Array<OneD,std::set<GeometrySharedPtr>> m_constituentElements;
    bool m_moved = true;
    int m_coordDim;
};

struct ZoneRotate final: public ZoneBase
{
    ZoneRotate(int id,
               const CompositeMap &domain,
               const int coordDim,
               const NekPoint<NekDouble> &origin,
               const DNekVec &axis,
               const LibUtilities::EquationSharedPtr &angularVelEqn);

    virtual ~ZoneRotate() = default;

    inline NekPoint<NekDouble> GetOrigin() const
    {
        return m_origin;
    }

    inline NekVector<NekDouble> GetAxis() const
    {
        return m_axis;
    }

    NekDouble GetAngularVel(NekDouble &time) const;

    virtual bool v_Move(NekDouble timeStep) final;

protected:
    NekPoint<NekDouble> m_origin;
    DNekVec m_axis;
    LibUtilities::EquationSharedPtr m_angularVelEqn;
    std::vector<PointGeomSharedPtr> m_rotateVerts;
    std::vector<CurveSharedPtr> m_rotateCurves;
    std::vector<PointGeom> m_origPosition;
    DNekMat m_W  = DNekMat(3, 3, 0.0);
    DNekMat m_W2 = DNekMat(3, 3, 0.0);

};

struct ZoneTranslate final: public ZoneBase
{
    ZoneTranslate(int id,
                  const CompositeMap &domain,
                  const int coordDim,
                  const std::vector<NekDouble> &velocity);

    virtual ~ZoneTranslate() = default;

    inline std::vector<NekDouble> GetVel() const
    {
        return m_velocity;
    }

    virtual bool v_Move(NekDouble timeStep) final;

protected:
    std::vector<NekDouble> m_velocity;
    std::vector<PointGeomSharedPtr> m_slideVerts;
    std::vector<CurveSharedPtr> m_slideCurves;
    std::vector<PointGeom> m_origPosition;
};

struct ZonePrescribe final: public ZoneBase
{
    ZonePrescribe(int id,
                  const CompositeMap &domain,
                  const int coordDim,
                  LibUtilities::EquationSharedPtr xDeform,
                  LibUtilities::EquationSharedPtr yDeform);

    virtual ~ZonePrescribe() = default;

    inline NekDouble GetXDeform(NekDouble x, NekDouble y,  NekDouble z, NekDouble t) const
    {
        return m_xDeform->Evaluate(x, y, z, t);
    }

    inline NekDouble GetYDeform(NekDouble x, NekDouble y,  NekDouble z, NekDouble t) const
    {
        return m_yDeform->Evaluate(x, y, z, t);
    }

    virtual bool v_Move(NekDouble timeStep) final;

protected:
    LibUtilities::EquationSharedPtr m_xDeform;
    LibUtilities::EquationSharedPtr m_yDeform;
    std::vector<PointGeomSharedPtr> m_interiorVerts;
    std::vector<PointGeom> m_origPosition;
};

struct ZoneFixed final: public ZoneBase
{
    ZoneFixed(int id,
              const CompositeMap &domain,
              const int coordDim)
            : ZoneBase(MovementType::eFixed, id, domain, coordDim)
    {
    }

    virtual ~ZoneFixed() = default;

    virtual bool v_Move(NekDouble timeStep) final;

};

typedef std::shared_ptr<ZoneRotate> ZoneRotateShPtr;
typedef std::shared_ptr<ZoneTranslate> ZoneTranslateShPtr;
typedef std::shared_ptr<ZonePrescribe> ZonePrescribeShPtr;
typedef std::shared_ptr<ZoneFixed> ZoneFixedShPtr;

struct Interface;
typedef std::shared_ptr<Interface> InterfaceShPtr;
struct Interface
{
    Interface(int indx, InterfaceSide side, CompositeMap edge);

    virtual ~Interface() = default;

    inline std::map<int, GeometrySharedPtr> const &GetEdge() const
    {
        return m_edge;
    }

    inline GeometrySharedPtr const &GetEdge(int id)
    {
        return m_edge[id];
    }

    inline std::deque<GeometrySharedPtr> const &GetEdgeDeque() const
    {
        return m_edgeDeque;
    }

    inline std::vector<int> const &GetEdgeIds() const
    {
        return m_edgeIds;
    }

    inline bool IsEmpty() const
    {
        return m_edge.empty();
    }

    inline void SetEdge(const GeometrySharedPtr &edge)
    {
        m_edge[edge->GetGlobalID()] = edge;
    }

    void SetEdge(const CompositeMap &edge);

    inline void SetOppInterface(const InterfaceShPtr &oppInterface)
    {
        m_oppInterface = oppInterface;
    }

    inline InterfaceShPtr GetOppInterface()
    {
        return m_oppInterface;
    }

    inline int &GetId()
    {
        return m_id;
    }

    inline InterfaceSide GetSide() const
    {
        return m_side;
    }

    inline void SetSide(const InterfaceSide &side)
    {
        m_side = side;
    }

protected:
    InterfaceShPtr m_oppInterface;
    int m_id;
    InterfaceSide m_side = InterfaceSide::eNone;
    std::map<int, GeometrySharedPtr> m_edge;
    std::deque<GeometrySharedPtr> m_edgeDeque;
    std::vector<int> m_edgeIds;
};

struct InterfacePair
{
    InterfacePair(const InterfaceShPtr &leftInterface,
                  const InterfaceShPtr &rightInterface)
                 : m_leftInterface(leftInterface),
                   m_rightInterface(rightInterface)
    {
        leftInterface->SetOppInterface(rightInterface);
        rightInterface->SetOppInterface(leftInterface);
    }

    InterfaceShPtr m_leftInterface;
    InterfaceShPtr m_rightInterface;

public:
    inline const InterfaceShPtr &GetLeftInterface() const
    {
        return m_leftInterface;
    }

    inline const InterfaceShPtr &GetRightInterface() const
    {
        return m_rightInterface;
    }
};

typedef std::shared_ptr<InterfacePair> InterfacePairShPtr;
typedef std::map<std::pair<int, std::string>, InterfacePairShPtr> InterfaceCollection;
typedef std::shared_ptr<InterfaceCollection> InterfaceCollectionShPtr;

class Movement
{
public:
    SPATIAL_DOMAINS_EXPORT Movement(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const MeshGraphSharedPtr &meshGraph);

    SPATIAL_DOMAINS_EXPORT Movement() = default;

    inline const InterfaceCollection &GetInterfaces() const
    {
        return m_interfaces;
    }

    inline const std::map<int, ZoneBaseShPtr> &GetZones() const
    {
        return m_zones;
    }

    void PerformMovement(NekDouble timeStep);

    void GenGeomFactors();

    inline const bool &GetMoveFlag() const
    {
        return m_moveFlag;
    }

protected:
    /// The mesh graph to use for referencing geometry info.
    MeshGraphSharedPtr m_meshGraph;
    LibUtilities::SessionReaderSharedPtr m_session;
    InterfaceCollection m_interfaces;
    std::map<int, ZoneBaseShPtr> m_zones;
    bool m_moveFlag = false;


private:
    /// Read interfaces (and general MeshGraph) given TiXmlDocument.
    void Read(TiXmlElement *movementTag);
    void ReadZones(TiXmlElement *zonesTag);
    void ReadInterfaces(TiXmlElement *interfacesTag);
};

typedef std::shared_ptr<Movement> MovementSharedPtr;

} // namespace SpatialDomains
} // namespace Nektar

#endif // NEKTAR_SPATIALDOMAINS_INTERFACES_H
