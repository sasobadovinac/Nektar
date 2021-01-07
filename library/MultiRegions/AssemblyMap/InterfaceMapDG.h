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

class InterfaceTrace
{
public:
    /// Default constructor
    MULTI_REGIONS_EXPORT InterfaceTrace(
        const ExpListSharedPtr &trace,
        const SpatialDomains::InterfaceBaseShPtr &interfaceBaseShPtr,
        const std::map<int, int> &geomIdToTraceId);

    /// Default destructor
    MULTI_REGIONS_EXPORT virtual ~InterfaceTrace() = default;

    inline void SetCheckLocal(bool flag)
    {
        m_checkLocal = flag;
    }

    inline bool &RecalcCoords()
    {
        return m_moved;
    }

    inline std::vector<Array<OneD, NekDouble>> GetMissingCoords()
    {
        return m_missingCoords;
    }

    SpatialDomains::InterfaceBaseShPtr GetInterface()
    {
        return m_interfaceBase;
    }

    void CalcLocalMissing();

    void FillLocalBwdTrace(Array<OneD, NekDouble> &Fwd,
                           Array<OneD, NekDouble> &Bwd);
    void FillRankBwdTrace(Array<OneD, NekDouble> &trace,
                          Array<OneD, NekDouble> &Bwd);

private:
    ExpListSharedPtr m_trace;
    SpatialDomains::InterfaceBaseShPtr m_interfaceBase;
    std::map<int, int> m_geomIdToTraceId;
    int m_totQuadPnts = 0;
    bool m_checkLocal = false;
    bool m_moved = true;
    std::vector<Array<OneD, NekDouble>> m_missingCoords;
    std::vector<std::pair<int, Array<OneD, NekDouble>>> m_foundLocalCoords;
    std::vector<int> m_mapMissingCoordToTrace;
    std::vector<int> m_mapFoundCoordToTrace;
};

typedef std::shared_ptr<InterfaceTrace> InterfaceTraceSharedPtr;

class InterfaceExchange
{
public:
    /// Default destructor
    MULTI_REGIONS_EXPORT virtual ~InterfaceExchange() = default;

    /// Default constructor
    MULTI_REGIONS_EXPORT InterfaceExchange(
        const ExpListSharedPtr &trace, const LibUtilities::CommSharedPtr &comm,
        std::pair<int, std::vector<InterfaceTraceSharedPtr>> rankPair,
        const std::map<int, int> &geomIdToTraceId)
        : m_trace(trace), m_comm(comm), m_rank(rankPair.first),
          m_interfaces(rankPair.second), m_geomIdToTraceId(geomIdToTraceId)
    {
    }

    MULTI_REGIONS_EXPORT void RankFillSizes(
        LibUtilities::CommRequestSharedPtr request, int requestNum);
    MULTI_REGIONS_EXPORT void SendMissing(
        LibUtilities::CommRequestSharedPtr request, int requestNum);
    MULTI_REGIONS_EXPORT void CalcRankDistances();
    MULTI_REGIONS_EXPORT void SendFwdTrace(
        LibUtilities::CommRequestSharedPtr request, int requestNum,
        Array<OneD, NekDouble> &Fwd);
    MULTI_REGIONS_EXPORT void FillRankBwdTraceExchange(
        Array<OneD, NekDouble> &Bwd);

private:
    const ExpListSharedPtr m_trace;
    const LibUtilities::CommSharedPtr m_comm;
    int m_rank;
    const std::vector<InterfaceTraceSharedPtr> m_interfaces;
    Array<OneD, int> m_sendSize;
    Array<OneD, int> m_recvSize;
    int m_totSendSize = 0;
    int m_totRecvSize = 0;
    Array<OneD, NekDouble> m_send;
    Array<OneD, NekDouble> m_recv;
    std::map<int, std::pair<int, Array<OneD, NekDouble>>> m_foundRankCoords;
    Array<OneD, NekDouble> m_recvTrace;
    std::map<int, int> m_geomIdToTraceId;
};

typedef std::shared_ptr<InterfaceExchange> InterfaceExchangeSharedPtr;

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

    MULTI_REGIONS_EXPORT void ExchangeTrace(Array<OneD, NekDouble> &Fwd,
                                            Array<OneD, NekDouble> &Bwd);
    MULTI_REGIONS_EXPORT void ExchangeCoords();

private:
    SpatialDomains::InterfacesSharedPtr m_interfaces;
    std::vector<InterfaceTraceSharedPtr> m_localInterfaces;
    const ExpListSharedPtr m_trace;
    std::vector<InterfaceExchangeSharedPtr> m_exchange;
    std::map<int, int> m_geomIdToTraceId;
};

typedef std::shared_ptr<InterfaceMapDG> InterfaceMapDGSharedPtr;

} // namespace MultiRegions
} // namespace Nektar

#endif