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

#include <map>
#include <string>
#include <deque>

#include <SpatialDomains/MeshGraph.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

namespace Nektar
{

namespace SpatialDomains
{

enum class InterfaceSide
{
    eNone,
    eLeft,
    eRight
};

const std::string InterfaceSideStr[] = {"NONE", "LEFT", "RIGHT"};

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

    // Using these as various caches for the InterfaceMapDG class
    // This allows caching of calculated values across different fields when searching
    // for missing coordinates locally on own rank
    std::vector<Array<OneD, NekDouble>> m_missingCoords;
    std::map<int, std::pair<int, Array<OneD, NekDouble>>> m_foundLocalCoords;
    std::vector<int> m_mapMissingCoordToTrace;
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

} // namespace SpatialDomains
} // namespace Nektar

#endif // NEKTAR_INTERFACEINTERPOLATION_MOVEMENT_H