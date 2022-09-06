////////////////////////////////////////////////////////////////////////////////
//
//  File: Movement.h
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
#ifndef NEKTAR_SPATIALDOMAINS_MOVEMENT_H
#define NEKTAR_SPATIALDOMAINS_MOVEMENT_H

#include <SpatialDomains/Movement/InterfaceInterpolation.h>
#include <SpatialDomains/Movement/Zones.h>

namespace Nektar
{

namespace SpatialDomains
{

typedef std::map<std::pair<int, std::string>, InterfacePairShPtr>
    InterfaceCollection;

class Movement
{
public:
    /// Constructor
    SPATIAL_DOMAINS_EXPORT Movement(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        MeshGraph *meshGraph);

    /// Default destructor
    SPATIAL_DOMAINS_EXPORT ~Movement() = default;

    inline const InterfaceCollection &GetInterfaces() const
    {
        return m_interfaces;
    }

    inline const std::map<int, ZoneBaseShPtr> &GetZones() const
    {
        return m_zones;
    }

    void PerformMovement(NekDouble timeStep);

    inline const bool &GetMoveFlag() const
    {
        return m_moveFlag;
    }

    inline bool &GetCoordExchangeFlag()
    {
        return m_coordExchangeFlag;
    }

protected:
    InterfaceCollection m_interfaces;
    std::map<int, ZoneBaseShPtr> m_zones;
    bool m_moveFlag = false; // Flags presence of moving zones
    bool m_coordExchangeFlag =
        true; // Flags if missing coordinates need to be calculated

private:
    /// Read zones given TiXmlDocument
    void ReadZones(TiXmlElement *zonesTag, MeshGraph *meshGraph,
                   const LibUtilities::SessionReaderSharedPtr &pSession);
    /// Read interfaces given TiXmlDocument
    void ReadInterfaces(TiXmlElement *interfacesTag, MeshGraph *meshGraph);
};

typedef std::shared_ptr<Movement> MovementSharedPtr;

} // namespace SpatialDomains
} // namespace Nektar

#endif // NEKTAR_SPATIALDOMAINS_MOVEMENT_H