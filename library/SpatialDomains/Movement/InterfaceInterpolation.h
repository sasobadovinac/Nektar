////////////////////////////////////////////////////////////////////////////////
//
//  File: InterfaceInterpolation.h
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
#ifndef NEKTAR_SPATIALDOMAINS_INTERFACEINTERPOLATION_H
#define NEKTAR_SPATIALDOMAINS_INTERFACEINTERPOLATION_H

#include <SpatialDomains/GeomFactors.h>

namespace Nektar
{

namespace SpatialDomains
{

// Fwd def to allow for inclusion in meshgraph
struct Composite;
typedef std::shared_ptr<Composite> CompositeSharedPtr;
typedef std::map<int, CompositeSharedPtr> CompositeMap;

/// A interface which is a single edge on a zone for handling non-conformality
struct Interface
{
    /// Constructor
    Interface(int indx, const CompositeMap &edge);

    /// Default destructor
    virtual ~Interface() = default;

    /// Returns map of global ID to geometry of the interface edge
    inline std::map<int, GeometrySharedPtr> const &GetEdge() const
    {
        return m_edge;
    }

    /// Returns geometry of the interface edge with global ID @param id
    inline GeometrySharedPtr const &GetEdge(int id)
    {
        return m_edge[id];
    }


    /// Checks if the interface edge is empty (used for parallelisation)
    inline bool IsEmpty() const
    {
        return m_edge.empty();
    }

    /// Returns the matching opposite interface from the interface pair
    inline std::shared_ptr<Interface> &GetOppInterface()
    {
        return m_oppInterface;
    }

    /// Returns the interface ID
    inline int &GetId()
    {
        return m_id;
    }

protected:
    /// Matching opposite interface of the interface pair
    std::shared_ptr<Interface> m_oppInterface;
    /// Interface ID
    int m_id;
    /// Map of global ID to geometry of the interface edge
    std::map<int, GeometrySharedPtr> m_edge;
};

typedef std::shared_ptr<Interface> InterfaceShPtr;

/// Interface pair consisting of a 'left' and 'right' interface
struct InterfacePair
{
    /// Constructor
    InterfacePair(const InterfaceShPtr &leftInterface,
                  const InterfaceShPtr &rightInterface)
        : m_leftInterface(leftInterface),
          m_rightInterface(rightInterface)
    {
        // Sets the opposite interfaces
        leftInterface->GetOppInterface() = rightInterface;
        rightInterface->GetOppInterface() = leftInterface;
    }

    /// 'Left' interface of the interface pair
    InterfaceShPtr m_leftInterface;
    /// 'Right' interface of the interface pair
    InterfaceShPtr m_rightInterface;

public:
    /// Return the 'left' interface from the interface pair
    inline const InterfaceShPtr &GetLeftInterface() const
    {
        return m_leftInterface;
    }

    /// Return the 'right' interface from the interface pair
    inline const InterfaceShPtr &GetRightInterface() const
    {
        return m_rightInterface;
    }
};

typedef std::shared_ptr<InterfacePair> InterfacePairShPtr;

} // namespace SpatialDomains
} // namespace Nektar

#endif // NEKTAR_INTERFACEINTERPOLATION_MOVEMENT_H