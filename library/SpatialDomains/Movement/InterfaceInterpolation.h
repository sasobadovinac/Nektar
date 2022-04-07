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

#include <SpatialDomains/MeshGraph.h>

namespace Nektar
{

namespace SpatialDomains
{


struct Interface;
typedef std::shared_ptr<Interface> InterfaceShPtr;

struct Interface
{
    Interface(int indx, const CompositeMap &edge);

    virtual ~Interface() = default;

    inline std::map<int, GeometrySharedPtr> const &GetEdge() const
    {
        return m_edge;
    }

    inline GeometrySharedPtr const &GetEdge(int id)
    {
        return m_edge[id];
    }

    inline bool IsEmpty() const
    {
        return m_edge.empty();
    }

    inline InterfaceShPtr &GetOppInterface()
    {
        return m_oppInterface;
    }

    inline int &GetId()
    {
        return m_id;
    }

protected:
    InterfaceShPtr m_oppInterface;
    int m_id;
    std::map<int, GeometrySharedPtr> m_edge;
};

struct InterfacePair
{
    InterfacePair(const InterfaceShPtr &leftInterface,
                  const InterfaceShPtr &rightInterface)
        : m_leftInterface(leftInterface),
          m_rightInterface(rightInterface)
    {
        // Sets the opposite interfaces
        leftInterface->GetOppInterface() = rightInterface;
        rightInterface->GetOppInterface() = leftInterface;
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

} // namespace SpatialDomains
} // namespace Nektar

#endif // NEKTAR_INTERFACEINTERPOLATION_MOVEMENT_H