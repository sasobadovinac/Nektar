////////////////////////////////////////////////////////////////////////////////
//
//  File: Zones.h
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
#ifndef NEKTAR_SPATIALDOMAINS_ZONES_H
#define NEKTAR_SPATIALDOMAINS_ZONES_H

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

const std::string MovementTypeStr[] = {"None",
                                       "Fixed",
                                       "Rotated",
                                       "Translated",
                                       "Prescribed"};

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

typedef std::shared_ptr<ZoneBase> ZoneBaseShPtr;

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
                  LibUtilities::EquationSharedPtr yDeform,
                  LibUtilities::EquationSharedPtr zDeform);

    virtual ~ZonePrescribe() = default;

    inline NekDouble GetXDeform(NekDouble x, NekDouble y,  NekDouble z, NekDouble t) const
    {
        return m_xDeform->Evaluate(x, y, z, t);
    }

    inline NekDouble GetYDeform(NekDouble x, NekDouble y,  NekDouble z, NekDouble t) const
    {
        return m_yDeform->Evaluate(x, y, z, t);
    }

    inline NekDouble GetZDeform(NekDouble x, NekDouble y,  NekDouble z, NekDouble t) const
    {
        return m_zDeform->Evaluate(x, y, z, t);
    }

    virtual bool v_Move(NekDouble timeStep) final;

protected:
    LibUtilities::EquationSharedPtr m_xDeform;
    LibUtilities::EquationSharedPtr m_yDeform;
    LibUtilities::EquationSharedPtr m_zDeform;
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

} // namespace SpatialDomains
} // namespace Nektar

#endif // NEKTAR_SPATIALDOMAINS_ZONES_H