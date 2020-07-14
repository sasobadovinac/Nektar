///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyCommDG.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Communication for interfaces
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_INTERFACEMAPDG_H
#define NEKTAR_INTERFACEMAPDG_H

#include <MultiRegions/ExpList.h>
#include <SpatialDomains/Interface.h>

namespace Nektar
{
namespace MultiRegions
{

class InterfaceExchange
{

public:
    /// Default destructor
    MULTI_REGIONS_EXPORT virtual ~InterfaceExchange() = default;

    /// Default constructor
    MULTI_REGIONS_EXPORT InterfaceExchange(
        const LibUtilities::CommSharedPtr &comm,
        std::pair<int, std::vector<SpatialDomains::InterfaceBaseShPtr>> rankPair)
        : m_comm(comm),
          m_rank(rankPair.first),
          m_interfaces(rankPair.second)
    {
    }

    MULTI_REGIONS_EXPORT void RankFillSizes(LibUtilities::CommRequestSharedPtr request, int requestNum);
    MULTI_REGIONS_EXPORT void RankCoordCalc(LibUtilities::CommRequestSharedPtr request, int requestNum);
private:
    const LibUtilities::CommSharedPtr m_comm;
    int m_rank;
    const std::vector<SpatialDomains::InterfaceBaseShPtr> m_interfaces;
    Array<OneD, int> m_sendSize;
    Array<OneD, int> m_recvSize;
    int m_totSendSize = 0;
    int m_totRecvSize = 0;
    Array<OneD, NekDouble> m_send;
    Array<OneD, NekDouble> m_recv;
    std::map<int, bool> m_checkLocal;
};

typedef std::shared_ptr<InterfaceExchange>  InterfaceExchangeSharedPtr;

class InterfaceMapDG
{
public:
    /// Default destructor
    MULTI_REGIONS_EXPORT ~InterfaceMapDG() = default;

    // Constructor for interface communication
    MULTI_REGIONS_EXPORT InterfaceMapDG(
        const SpatialDomains::InterfacesSharedPtr &interfaces,
        const ExpListSharedPtr &trace,
        const std::map<int, int> geomIdToTraceId);

    void CalcLocalCoords(const ExpListSharedPtr &trace, std::map<int, int> geomIdToTraceId);

private:
    SpatialDomains::InterfacesSharedPtr m_interfaces;
    std::vector<SpatialDomains::InterfaceBaseShPtr> m_localInterfaces;
    const ExpListSharedPtr m_trace;
    std::vector<InterfaceExchangeSharedPtr> m_exchange;
    std::map<int, int> m_geomIdToTraceId;
};

typedef std::shared_ptr<InterfaceMapDG>  InterfaceMapDGSharedPtr;


} // namespace MultiRegions
} // namespace Nektar

#endif