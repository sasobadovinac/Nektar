///////////////////////////////////////////////////////////////////////////////
//
// File InterfaceMapDG.h
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
// Description: MPI communication for interfaces, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_INTERFACEMAPDG_H
#define NEKTAR_INTERFACEMAPDG_H

#include <MultiRegions/ExpList.h>
#include <SpatialDomains/Movement/Movement.h>

namespace Nektar
{
namespace MultiRegions
{

/**
 * Object for each interface present between two ranks, which are held in the
 * InterfaceExchange object.
 */

class InterfaceTrace
{
public:
    /// Default constructor
    MULTI_REGIONS_EXPORT InterfaceTrace(
        const ExpListSharedPtr &trace,
        const SpatialDomains::InterfaceShPtr &interfaceShPtr);

    /// Default destructor
    MULTI_REGIONS_EXPORT virtual ~InterfaceTrace() = default;

    inline void SetCheckLocal(bool flag)
    {
        m_checkLocal = flag;
    }

    /// Returns the missing coordinates vector
    inline std::vector<Array<OneD, NekDouble>> GetMissingCoords()
    {
        return m_missingCoords;
    }

    /// Returns the interface object
    SpatialDomains::InterfaceShPtr GetInterface()
    {
        return m_interface;
    }

    /// Calculates what coordinates on the interface are missing locally
    void CalcLocalMissing();
    /// Fills the Bwd trace by interpolating from the Fwd for local interfaces
    void FillLocalBwdTrace(Array<OneD, NekDouble> &Fwd,
                           Array<OneD, NekDouble> &Bwd);
    /// Fills the Bwd trace from partitioned trace
    void FillRankBwdTrace(Array<OneD, NekDouble> &trace,
                          Array<OneD, NekDouble> &Bwd);

private:
    /// Trace expansion list
    ExpListSharedPtr m_trace;
    /// Local interface object
    SpatialDomains::InterfaceShPtr m_interface;
    /// Flag whether the opposite side of the interface is present locally
    bool m_checkLocal = false;
    /// Vector of coordinates on interface not present on opposite side local
    std::vector<Array<OneD, NekDouble>> m_missingCoords;
    /// Map of found coordinates present locally
    std::map<int, std::pair<int, Array<OneD, NekDouble>>> m_foundLocalCoords;
    /// Vector of indices corresponding to m_missingCoord locations in trace
    std::vector<int> m_mapMissingCoordToTrace;
};

typedef std::shared_ptr<InterfaceTrace> InterfaceTraceSharedPtr;

/**
 * Object for one rank-to-rank communication for all interfaces shared
 * between those ranks. e.g. if on rank 1 and there are interfaces shared with
 * rank 2 and 3, there will be two InterfaceExchange objects in m_exchange in
 * InterfaceMapDG. This holds the InterfaceTrace objects in m_interfaceTraces.
 */
class InterfaceExchange
{
public:
    /// Default destructor
    MULTI_REGIONS_EXPORT virtual ~InterfaceExchange() = default;

    /// Constructor
    MULTI_REGIONS_EXPORT InterfaceExchange(
        SpatialDomains::MovementSharedPtr movement,
        const ExpListSharedPtr &trace, const LibUtilities::CommSharedPtr &comm,
        std::pair<int, std::vector<InterfaceTraceSharedPtr>> rankPair)
        : m_movement(movement), m_zones(movement->GetZones()), m_trace(trace),
          m_comm(comm), m_rank(rankPair.first),
          m_interfaceTraces(rankPair.second)
    {
    }

    /**
     * Communicates with other ranks how many missing coordinates for each
     * interface to expect
     *
     * @param requestSend List of send requests
     * @param requestRecv List of receive requests
     * @param requestNum Index of request in list to use
     */
    MULTI_REGIONS_EXPORT void RankFillSizes(
        LibUtilities::CommRequestSharedPtr &requestSend,
        LibUtilities::CommRequestSharedPtr &requestRecv, int requestNum);

    /**
     * Sends/receives the missing coordinates to/from other ranks
     *
     * @param requestSend List of send requests
     * @param requestRecv List of receive requests
     * @param requestNum Index of request in list to use
     */
    MULTI_REGIONS_EXPORT void SendMissing(
        LibUtilities::CommRequestSharedPtr &requestSend,
        LibUtilities::CommRequestSharedPtr &requestRecv, int requestNum);

    /// Populates m_foundRankCoords using the FindDistance function
    MULTI_REGIONS_EXPORT void CalcRankDistances();

    /**
     * Calculates and sends the trace to other rank from the m_foundRankCoords
     * structure using non-blocking pairwise communication i.e. Isend & Irecv
     *
     * @param requestSend List of send requests
     * @param requestRecv List of receive requests
     * @param requestNum Index of request in list to use
     * @param Fwd The values to send across the interface
     */
    MULTI_REGIONS_EXPORT void SendFwdTrace(
        LibUtilities::CommRequestSharedPtr &requestSend,
        LibUtilities::CommRequestSharedPtr &requestRecv, int requestNum,
        Array<OneD, NekDouble> &Fwd);

    /**
     * Loops over interfaces and partitions out the received trace from the
     * other ranks for insertion into Bwd using FillRankBwdTrace
     *
     * @param Bwd The Bwd trace to be filled from across the interface
     */
    MULTI_REGIONS_EXPORT void FillRankBwdTraceExchange(
        Array<OneD, NekDouble> &Bwd);

private:
    /// Movement object associated with the non-conformal interfaces
    SpatialDomains::MovementSharedPtr m_movement;
    /// Map of zone IDs to zone bases
    std::map<int, SpatialDomains::ZoneBaseShPtr> m_zones;
    /// Trace expansion list
    const ExpListSharedPtr m_trace;
    /// Communicator
    const LibUtilities::CommSharedPtr m_comm;
    /// Process rank
    int m_rank;
    /// Vector of interface traces i.e. every interface side present on m_rank
    const std::vector<InterfaceTraceSharedPtr> m_interfaceTraces;
    /// Send buffer for coord exchange
    Array<OneD, NekDouble> m_send;
    /// Receive buffer for coord exchange
    Array<OneD, NekDouble> m_recv;
    /// Receive buffer for trace exchange
    Array<OneD, NekDouble> m_recvTrace;
    /// Send buffer for trace exchange
    Array<OneD, NekDouble> m_sendTrace;
    /// Map of rank to total size of send buffer for all interfaces
    std::map<int, int> m_totSendSize;
    /// Map of rank to total size of receive buffer for all interfaces
    std::map<int, int> m_totRecvSize;
    /// Map of rank to array of size of send buffer for each interface
    std::map<int, Array<OneD, int>> m_sendSize;
    /// Map of rank to array of size of receive buffer for each interface
    std::map<int, Array<OneD, int>> m_recvSize;

    /**
     *  Caches the found coordinates to reuse when exchanging the trace in a
     *  map of integer rank to a map of integer missing coordinate location to a
     *  pair of local edge ID and found local coordinate
     */
    std::map<int, std::map<int, std::pair<int, Array<OneD, NekDouble>>>>
        m_foundRankCoords;
};

typedef std::shared_ptr<InterfaceExchange> InterfaceExchangeSharedPtr;

/**
 * Implements the communication patterns to allow for exchange of information
 * across non-conformal interfaces and across different partitions. Holds all
 * the InterfaceExchange objects in m_exchange.
 */
class InterfaceMapDG
{
public:
    /// Default destructor
    MULTI_REGIONS_EXPORT ~InterfaceMapDG() = default;

    /**
     * Sets up the InterfaceExchange objects stored in m_exchange, each object
     * is rank -> rank and contains a vector of InterfaceTrace objects
     * corresponding to shared interfaces between those ranks.
     */
    MULTI_REGIONS_EXPORT InterfaceMapDG(
        const SpatialDomains::MeshGraphSharedPtr &graph,
        const ExpListSharedPtr &trace);

    /**
     * @brief Perform the trace exchange between processors, given the forwards
     * and backwards spaces.
     *
     * @param Fwd  Local forwards space of the trace (which will be sent)
     * @param Bwd  Local backwards space of the trace (which will receive
     *                 contributions)
     */
    MULTI_REGIONS_EXPORT void ExchangeTrace(Array<OneD, NekDouble> &Fwd,
                                            Array<OneD, NekDouble> &Bwd);
    /**
     * @brief Perform the coordinate exchange between processors. This is where
     * the missing coordinates on the interface are found and sent to all other
     * processors on the other side of that interface so the matching ranks can
     * be found.
     */
    MULTI_REGIONS_EXPORT void ExchangeCoords();

private:
    /// Mesh associated with this expansion list.
    SpatialDomains::MeshGraphSharedPtr m_graph;
    /// Movement object associated with the non-conformal interfaces
    SpatialDomains::MovementSharedPtr m_movement;
    /// Interface sides present on current process
    std::vector<InterfaceTraceSharedPtr> m_localInterfaces;
    /// Trace expansion list
    const ExpListSharedPtr m_trace;
    /// Vector of interface exchanges, i.e. every rank-to-rank comm needed
    std::vector<InterfaceExchangeSharedPtr> m_exchange;
};

typedef std::shared_ptr<InterfaceMapDG> InterfaceMapDGSharedPtr;

} // namespace MultiRegions
} // namespace Nektar

#endif