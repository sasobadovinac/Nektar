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
//  Description: Zones used in the non-conformal interfaces
//               and ALE implementations
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_ZONES_H
#define NEKTAR_SPATIALDOMAINS_ZONES_H

#include <LibUtilities/BasicUtils/Equation.h>

namespace Nektar
{

namespace SpatialDomains
{

/// Enum of zone movement type
enum class MovementType
{
    eNone,
    eFixed,
    eRotate,
    eTranslate,
    ePrescribe
};

/// Map of zone movement type to movement type string
const std::string MovementTypeStr[] = {"None", "Fixed", "Rotated", "Translated",
                                       "Prescribed"};

/// Zone base: Contains the shared functions and variables
struct ZoneBase
{
    /// Constructor
    ZoneBase(MovementType type, int indx, CompositeMap domain, int coordDim);

    /// Default destructor
    virtual ~ZoneBase() = default;

    /// Returns the type of movement
    inline MovementType GetMovementType() const
    {
        return m_type;
    }

    /// Returns the domain the zone is on
    inline CompositeMap GetDomain() const
    {
        return m_domain;
    }

    /// Returns the zone ID
    inline int &GetId()
    {
        return m_id;
    }

    /// Virtual function for movement of the zone at @param time
    inline virtual bool v_Move(NekDouble time)
    {
        boost::ignore_unused(time);
        return false;
    }

    /// Performs the movement of the zone at @param time
    inline bool Move(NekDouble time)
    {
        return v_Move(time);
    }

    /// Returns all highest dimension elements in the zone
    inline std::vector<GeometrySharedPtr> const &GetElements() const
    {
        return m_elements;
    }

    /// Returns the flag which states if the zone has moved in this timestep
    inline bool &GetMoved()
    {
        return m_moved;
    }

    /// Clears all bounding boxes associated with the zones elements
    void ClearBoundingBoxes();

    /// Returns constituent elements, i.e. faces + edges
    inline std::array<std::set<GeometrySharedPtr>, 3> &GetConstituentElements()
    {
        return m_constituentElements;
    }

protected:
    /// Type of zone movement
    MovementType m_type = MovementType::eNone;
    /// Zone ID
    int m_id;
    /// Zone domain
    CompositeMap m_domain;
    /// Vector of highest dimension zone elements
    std::vector<GeometrySharedPtr> m_elements;
    /// Array of all dimension elements i.e. faces = [2], edges = [1], geom =
    /// [0]
    std::array<std::set<GeometrySharedPtr>, 3> m_constituentElements;
    /// Moved flag
    bool m_moved = true;
    /// Coordinate dimension
    int m_coordDim;
    /// Vector of all points in the zone
    std::vector<PointGeomSharedPtr> m_verts;
    /// Vector of all curves in the zone
    std::vector<CurveSharedPtr> m_curves;
    /// Vector of all points in the zone at initialisation
    std::vector<PointGeom> m_origVerts;
};

typedef std::shared_ptr<ZoneBase> ZoneBaseShPtr;

/// Rotating zone: Motion of every point around a given axis on an origin
struct ZoneRotate final : public ZoneBase
{
    /**
     * Constructor for rotating zones
     *
     * @param id Zone ID
     * @param domain Domain that the zone consists of
     * @param coordDim Coordinate dimension
     * @param origin Origin that the zone rotates about
     * @param axis Axis that the zone rotates about
     * @param angularVelEqn Equation for the angular velocity of rotation
     */
    ZoneRotate(int id, const CompositeMap &domain, const int coordDim,
               const NekPoint<NekDouble> &origin, const DNekVec &axis,
               const LibUtilities::EquationSharedPtr &angularVelEqn);

    /// Default destructor
    virtual ~ZoneRotate() = default;

    /// Return the angular velocity of the zone at @param time
    NekDouble GetAngularVel(NekDouble &time) const;

    /// Virtual function for movement of the zone at @param time
    virtual bool v_Move(NekDouble time) final;

    inline NekPoint<NekDouble> GetOrigin() const
    {
        return m_origin;
    }

    inline DNekVec GetAxis() const
    {
        return m_axis;
    }

protected:
    ///  Origin point rotation is performed around
    NekPoint<NekDouble> m_origin;
    /// Axis rotation is performed around
    DNekVec m_axis;
    /// Equation defining angular velocity as a function of time
    LibUtilities::EquationSharedPtr m_angularVelEqn;
    /// W matrix Rodrigues' rotation formula, cross product of axis
    DNekMat m_W = DNekMat(3, 3, 0.0);
    /// W^2 matrix Rodrigues' rotation formula, cross product of axis squared
    DNekMat m_W2 = DNekMat(3, 3, 0.0);
};

/// Translating zone: addition of a constant vector to every point
struct ZoneTranslate final : public ZoneBase
{
    /**
     * Constructor for translating zone
     *
     * @param id Zone ID
     * @param domain Domain that the zone consists of
     * @param coordDim Coordinate dimension
     * @param velocity Vector of translation velocity in x,y,z direction
     */
    ZoneTranslate(int id, const CompositeMap &domain, const int coordDim,
        const Array<OneD, LibUtilities::EquationSharedPtr> &velocityEqns,
        const Array<OneD, LibUtilities::EquationSharedPtr> &displacementEqns)
        : ZoneBase(MovementType::eTranslate, id, domain, coordDim),
          m_velocityEqns(velocityEqns), m_displacementEqns(displacementEqns)
    {
    }

    /// Default destructor
    virtual ~ZoneTranslate() = default;

    /// Returns the velocity of the zone
    std::vector<NekDouble> GetVel(NekDouble &time) const;

    /// Returns the displacement of the zone
    std::vector<NekDouble> GetDisp(NekDouble &time) const;

    /// Virtual function for movement of the zone at @param time
    virtual bool v_Move(NekDouble time) final;

protected:
    Array<OneD, LibUtilities::EquationSharedPtr> m_velocityEqns;
    Array<OneD, LibUtilities::EquationSharedPtr> m_displacementEqns;
};

/// Prescribed zone: applies equation to every point
struct ZonePrescribe final : public ZoneBase
{
    /**
     * Constructor for prescribed zone
     *
     * @param id Zone ID
     * @param domain Domain that the zone consists of
     * @param coordDim Coordinate dimension
     * @param xDeform Equation for prescribed motion of x-coordinate
     * @param yDeform Equation for prescribed motion of y-coordinate
     * @param zDeform Equation for prescribed motion of z-coordinate
     */
    ZonePrescribe(int id, const CompositeMap &domain, const int coordDim,
                  LibUtilities::EquationSharedPtr xDeform,
                  LibUtilities::EquationSharedPtr yDeform,
                  LibUtilities::EquationSharedPtr zDeform)
        : ZoneBase(MovementType::ePrescribe, id, domain, coordDim),
          m_xDeform(xDeform), m_yDeform(yDeform), m_zDeform(zDeform)
    {
    }

    /// Default destructor
    virtual ~ZonePrescribe() = default;

    /**
     * Returns point @param x @param y @param z deformation in the x direction
     * at time @param t
     * @param x x-coordinate
     * @param y y-coordinate
     * @param z z-coordinate
     * @param t time
     * @return deformation in x direction
     */
    inline NekDouble GetXDeform(NekDouble x, NekDouble y, NekDouble z,
                                NekDouble t) const
    {
        return m_xDeform->Evaluate(x, y, z, t);
    }

    /**
     * Returns point @param x @param y @param z deformation in the y direction
     * at time @param t
     * @param x x-coordinate
     * @param y y-coordinate
     * @param z z-coordinate
     * @param t time
     * @return deformation in y direction
     */
    inline NekDouble GetYDeform(NekDouble x, NekDouble y, NekDouble z,
                                NekDouble t) const
    {
        return m_yDeform->Evaluate(x, y, z, t);
    }

    /**
     * Returns point @param x @param y @param z deformation in the z direction
     * at time @param t
     * @param x x-coordinate
     * @param y y-coordinate
     * @param z z-coordinate
     * @param t time
     * @return deformation in z direction
     */
    inline NekDouble GetZDeform(NekDouble x, NekDouble y, NekDouble z,
                                NekDouble t) const
    {
        return m_zDeform->Evaluate(x, y, z, t);
    }

    /// Virtual function for movement of the zone at @param time
    virtual bool v_Move(NekDouble time) final;

protected:
    /// Equation specifying prescribed motion in x-direction
    LibUtilities::EquationSharedPtr m_xDeform;
    /// Equation specifying prescribed motion in y-direction
    LibUtilities::EquationSharedPtr m_yDeform;
    /// Equation specifying prescribed motion in z-direction
    LibUtilities::EquationSharedPtr m_zDeform;
};

/// Fixed zone: does not move
struct ZoneFixed final : public ZoneBase
{
    /// Constructor
    ZoneFixed(int id, const CompositeMap &domain, const int coordDim)
        : ZoneBase(MovementType::eFixed, id, domain, coordDim)
    {
    }

    /// Default destructor
    virtual ~ZoneFixed() = default;

    /// Virtual function for movement of the zone at @param time
    virtual bool v_Move(NekDouble time) final;
};

typedef std::shared_ptr<ZoneRotate> ZoneRotateShPtr;
typedef std::shared_ptr<ZoneTranslate> ZoneTranslateShPtr;
typedef std::shared_ptr<ZonePrescribe> ZonePrescribeShPtr;
typedef std::shared_ptr<ZoneFixed> ZoneFixedShPtr;

} // namespace SpatialDomains
} // namespace Nektar

#endif // NEKTAR_SPATIALDOMAINS_ZONES_H
