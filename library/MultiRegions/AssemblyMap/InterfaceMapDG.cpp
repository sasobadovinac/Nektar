///////////////////////////////////////////////////////////////////////////////
//
// File InterfaceMapDG.cpp
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
// Description: MPI communication for interfaces
//
///////////////////////////////////////////////////////////////////////////////

#include "InterfaceMapDG.h"

namespace Nektar
{
namespace MultiRegions
{

InterfaceTrace::InterfaceTrace(
    const ExpListSharedPtr &trace,
    const SpatialDomains::InterfaceShPtr &interfaceShPtr)
    : m_trace(trace), m_interface(interfaceShPtr)
{
}

// Creates the objects that communicates trace values across the interface in
// serial and parallel. Does this by checking all interfaces to see if they
// are contained locally (on the current rank), and checks all other ranks for
// the same interface ID and if contains the opposite interface side i.e.
// 'left' or 'right'
InterfaceMapDG::InterfaceMapDG(
    const SpatialDomains::MeshGraphSharedPtr &meshGraph,
    const ExpListSharedPtr &trace)
    : m_graph(meshGraph), m_movement(meshGraph->GetMovement()), m_trace(trace)
{
    auto comm                = m_trace->GetComm();
    auto interfaceCollection = m_movement->GetInterfaces();

    // myIndxLR contains the info about what interface edges are present on
    // current rank with each interface no, i,  consisting of:
    // [i] = indx
    // [i + 1] = 0 (non), = 1 (left only), = 2 (right only), = 3 (both)
    std::map<int, int> myIndxLRMap;

    // localInterfaces contains a map of interface ID to a pair of interface
    // traces, this pair is 'left' and 'right' interface and is used to
    // construct and store the traces for interfaces present on the current rank
    std::map<int, std::pair<InterfaceTraceSharedPtr, InterfaceTraceSharedPtr>>
        localInterfaces;

    // Map of m_localInterfaces vector to interface ID
    Array<OneD, int> indxToInterfaceID(interfaceCollection.size());

    // Loops over all interfaces and check if either side, 'left' or 'right'
    // are empty on the current rank, if not empty then populate the local
    // interface data structures and keep track using myIndxLRMap
    size_t cnt = 0;
    for (const auto &interface : interfaceCollection)
    {
        indxToInterfaceID[cnt]             = interface.first.first;
        myIndxLRMap[interface.first.first] = 0;

        if (!interface.second->GetLeftInterface()->IsEmpty())
        {
            myIndxLRMap[interface.first.first] += 1;

            localInterfaces[interface.first.first].first =
                MemoryManager<InterfaceTrace>::AllocateSharedPtr(
                    trace, interface.second->GetLeftInterface());
            m_localInterfaces.emplace_back(
                localInterfaces[interface.first.first].first);
        }

        if (!interface.second->GetRightInterface()->IsEmpty())
        {
            myIndxLRMap[interface.first.first] += 2;

            localInterfaces[interface.first.first].second =
                MemoryManager<InterfaceTrace>::AllocateSharedPtr(
                    trace, interface.second->GetRightInterface());
            m_localInterfaces.emplace_back(
                localInterfaces[interface.first.first].second);
        }

        cnt++;
    }

    // Send num of interfaces size so all partitions can prepare buffers
    int nRanks = comm->GetSize();

    // Send all interface edges present to all partitions
    Array<OneD, int> interfaceEdges(myIndxLRMap.size());
    cnt = 0;
    for (auto pres : myIndxLRMap)
    {
        interfaceEdges[cnt++] = pres.second;
    }
    Array<OneD, int> rankLocalInterfaceIds(myIndxLRMap.size() * nRanks, 0);
    comm->AllGather(interfaceEdges, rankLocalInterfaceIds);

    // Map of rank to vector of interface traces
    std::map<int, std::vector<InterfaceTraceSharedPtr>> oppRankSharedInterface;

    // Find what interface Ids match with other ranks, then check if opposite
    // edge
    size_t myRank     = comm->GetRank();
    size_t numInterfaces = interfaceCollection.size();
    for (int i = 0; i < nRanks; ++i)
    {
        for (size_t j = 0; j < numInterfaces; ++j)
        {
            int otherId   = indxToInterfaceID[j];

            // otherCode represents for a specific interface ID what sides are
            // present on rank i, 0=non, 1=left, 2=right, 3=both
            int otherCode = rankLocalInterfaceIds[i * numInterfaces + j];

            // myCode represents for a specific interface ID what sides are
            // present on this rank/process 0=non, 1=left, 2=right, 3=both
            int myCode    = myIndxLRMap[otherId];

            // Special case if checking current rank (this process)
            if (i == myRank)
            {
                // If contains both edges locally then set check local to true
                if (myCode == 3)
                {
                    localInterfaces[otherId].first->SetCheckLocal(true);
                    localInterfaces[otherId].second->SetCheckLocal(true);
                }

                continue;
            }

            // Checks if this ranks 'left' matches any 'right' on other rank
            if ((myCode == 1 && otherCode == 2) ||
                (myCode == 1 && otherCode == 3) ||
                (myCode == 3 && otherCode == 2))
            {
                oppRankSharedInterface[i].emplace_back(
                    localInterfaces[otherId].first);
            }
            // Checks if this ranks 'right' matches any 'left' on other rank
            else if ((myCode == 2 && otherCode == 1) ||
                     (myCode == 2 && otherCode == 3) ||
                     (myCode == 3 && otherCode == 1))
            {
                oppRankSharedInterface[i].emplace_back(
                    localInterfaces[otherId].second);
            }
            // Checks if this ranks 'both' matches any 'both' on other rank
            else if (myCode == 3 && otherCode == 3)
            {
                oppRankSharedInterface[i].emplace_back(
                    localInterfaces[otherId].first);
                oppRankSharedInterface[i].emplace_back(
                    localInterfaces[otherId].second);
            }
        }
    }

    // Create individual interface exchange objects (each object is rank ->
    // rank) and contains a vector of interfaceTrace objects
    for (auto &rank : oppRankSharedInterface)
    {
        m_exchange.emplace_back(
            MemoryManager<InterfaceExchange>::AllocateSharedPtr(
                m_movement, m_trace, comm, rank));
    }

    // Find missing coordinates on interface from other side
    ExchangeCoords();
}

void InterfaceMapDG::ExchangeCoords()
{
    auto comm  = m_trace->GetComm();
    auto zones = m_movement->GetZones();

    for (auto &interfaceTrace : m_localInterfaces)
    {
        interfaceTrace->CalcLocalMissing();
    }

    // If parallel communication is needed
    if (!m_exchange.empty())
    {
        auto requestSend = comm->CreateRequest(m_exchange.size());
        auto requestRecv = comm->CreateRequest(m_exchange.size());

        for (int i = 0; i < m_exchange.size(); ++i)
        {
            m_exchange[i]->RankFillSizes(requestSend, requestRecv, i);
        }
        comm->WaitAll(requestSend);
        comm->WaitAll(requestRecv);

        for (int i = 0; i < m_exchange.size(); ++i)
        {
            m_exchange[i]->SendMissing(requestSend, requestRecv, i);
        }
        comm->WaitAll(requestSend);
        comm->WaitAll(requestRecv);

        for (auto &i : m_exchange)
        {
            i->CalcRankDistances();
        }
    }
}

/**
 * Calculates what coordinates on the interface are missing; i.e. aren't located
 * in an edge from the other side of the interface on this partition. These are
 * stored in the vector #m_missingCoords, and another vector of the same index
 * layout contains the location of that coordinate in the trace i.e. the
 * 'offset + i' value. It checks first if the interface is flagged local, which
 * denotes whether the opposite side is also present on this rank, if not then
 * all coordinates are missing. If local then the same structure as
 * CalcRankDistances is used to minimise computational cost.
 */
void InterfaceTrace::CalcLocalMissing()
{
    auto childEdge = m_interface->GetEdge();
    // If not flagged 'check local' then all points are missing
    if (!m_checkLocal)
    {
        int cnt = 0;
        for (auto &childId : childEdge)
        {
            auto childElmt =
                m_trace->GetExpFromGeomId(childId.first);
            size_t nq = childElmt->GetTotPoints();
            Array<OneD, NekDouble> xc(nq, 0.0), yc(nq, 0.0), zc(nq, 0.0);
            childElmt->GetCoords(xc, yc, zc);
            int offset =
                m_trace->GetPhys_Offset(m_trace->GetElmtToExpId(childId.first));

            for (int i = 0; i < nq; ++i, ++cnt)
            {
                Array<OneD, NekDouble> xs(3);
                xs[0] = xc[i];
                xs[1] = yc[i];
                xs[2] = zc[i];

                m_missingCoords.emplace_back(xs);
                m_mapMissingCoordToTrace.emplace_back(offset + i);
            }
        }
    }
    // Otherwise we need to check to see what points the other side of the
    // interface on this rank contains. This is done by looping over each
    // quadrature point in all the edges on this side of the interface, then
    // check if this quadrature point is also in any edges on the other side
    // of the interface. First use MinMaxCheck as a cheap method to eliminate
    // far away edges. If it passes MinMaxCheck perform the FindDistance search
    // which numerically searches for the closest local coordinate on the other
    // interface side edge and returns the cartesian distance from the search
    // coordinate to the found local coordinate. If determined to be close
    // enough then add to found local coordinates.
    else
    {
        for (auto &childId : childEdge)
        {
            auto childElmt =
                m_trace->GetExpFromGeomId(childId.first);
            size_t nq = childElmt->GetTotPoints();
            Array<OneD, NekDouble> xc(nq, 0.0), yc(nq, 0.0), zc(nq, 0.0);
            childElmt->GetCoords(xc, yc, zc);
            int offset =
                m_trace->GetPhys_Offset(m_trace->GetElmtToExpId(childId.first));

            for (int i = 0; i < nq; ++i)
            {
                bool found = false;
                Array<OneD, NekDouble> foundLocCoord;
                Array<OneD, NekDouble> xs(3);
                xs[0] = xc[i];
                xs[1] = yc[i];
                xs[2] = zc[i];

                auto parentEdge = m_interface->GetOppInterface()->GetEdge();
                for (auto &edge : parentEdge)
                {
                    // First check if inside the edge bounding box
                    if (!edge.second->MinMaxCheck(xs))
                    {
                        continue;
                    }

                    NekDouble dist =
                        edge.second->FindDistance(xs, foundLocCoord);
                    if (dist < NekConstants::kFindDistanceMin)
                    {
                        found                          = true;
                        m_foundLocalCoords[offset + i] = std::make_pair(
                            edge.second->GetGlobalID(), foundLocCoord);
                        break;
                    }
                }

                if (!found)
                {
                    m_missingCoords.emplace_back(xs);
                    m_mapMissingCoordToTrace.emplace_back(offset + i);
                }
            }
        }
    }

    // If running in serial there shouldn't be any missing coordinates.
    if (m_trace->GetComm()->IsSerial())
    {
        ASSERTL0(m_missingCoords.empty(),
                 "Missing " + std::to_string(m_missingCoords.size()) +
                     " coordinates on interface ID " +
                     std::to_string(m_interface->GetId()) +
                     " linked to interface ID " +
                     std::to_string(m_interface->GetOppInterface()->GetId()) +
                     ". Check both sides of the interface line up.");
    }
}

// This fills m_sendSize and m_recvSize with the correct sizes for the
// communication lengths for each rank when communicating the traces. I.e.
// however many quadrature points on the opposite side of the interface to
// send, and however many on this side to receive.
void InterfaceExchange::RankFillSizes(
    LibUtilities::CommRequestSharedPtr &requestSend,
    LibUtilities::CommRequestSharedPtr &requestRecv, int requestNum)
{

    // Recalc size is the number of interfaces missing points that need to be
    // communicated.
    int recalcSize = 0;
    for (auto &localInterface : m_interfaceTraces)
    {
        recalcSize +=
            m_zones[localInterface->GetInterface()->GetId()]->GetMoved();
    }

    m_sendSize[m_rank] = Array<OneD, int>(recalcSize);
    m_recvSize[m_rank] = Array<OneD, int>(recalcSize);

    int cnt = 0;
    for (auto &interface : m_interfaceTraces)
    {
        if (m_zones[interface->GetInterface()->GetId()]->GetMoved())
        {
            // Calculates the number of missing coordinates to be communicated
            // This size is 3 times the number of missing coordinates to allow
            // for 3D points
            m_sendSize[m_rank][cnt++] =
                interface->GetMissingCoords().size() * 3;
        }
    }

    // Uses non-blocking send/recvs to send communication sizes to other ranks
    m_comm->Isend(m_rank, m_sendSize[m_rank], m_interfaceTraces.size(),
                  requestSend, requestNum);
    m_comm->Irecv(m_rank, m_recvSize[m_rank], m_interfaceTraces.size(),
                  requestRecv, requestNum);
}

// Uses the calculated send/recv sizes to send all missing coordinates to
// other ranks.
void InterfaceExchange::SendMissing(
    LibUtilities::CommRequestSharedPtr &requestSend,
    LibUtilities::CommRequestSharedPtr &requestRecv, int requestNum)
{
    m_totSendSize[m_rank] = std::accumulate(m_sendSize[m_rank].begin(),
                                            m_sendSize[m_rank].end(), 0);
    m_totRecvSize[m_rank] = std::accumulate(m_recvSize[m_rank].begin(),
                                            m_recvSize[m_rank].end(), 0);

    m_send = Array<OneD, NekDouble>(m_totSendSize[m_rank]);
    m_recv = Array<OneD, NekDouble>(m_totRecvSize[m_rank]);

    int cnt = 0;
    for (auto &interface : m_interfaceTraces)
    {
        if (m_zones[interface->GetInterface()->GetId()]->GetMoved())
        {
            auto missing = interface->GetMissingCoords();
            for (auto coord : missing)
            {
                // Points are organised in blocks of three
                // x coord, then y coord, then z coord.
                for (int k = 0; k < 3; ++k, ++cnt)
                {
                    m_send[cnt] = coord[k];
                }
            }
        }
    }

    // Uses non-blocking send/recvs to send missing points to other ranks
    m_comm->Isend(m_rank, m_send, m_totSendSize[m_rank], requestSend,
                  requestNum);
    m_comm->Irecv(m_rank, m_recv, m_totRecvSize[m_rank], requestRecv,
                  requestNum);
}

/**
 * Performs the trace exchange across interfaces
 * 1) Calculate and send the Fwd trace to other ranks with pairwise comms
 * 2) While waiting for the send/recv to complete we fill the local interfaces
 *    Bwd trace
 * 3) Fill the remaining Bwd trace with the received trace data from other ranks
 */
void InterfaceMapDG::ExchangeTrace(Array<OneD, NekDouble> &Fwd,
                                   Array<OneD, NekDouble> &Bwd)
{
    auto comm = m_trace->GetComm();

    // If no parallel exchange needed we only fill the local traces
    if (m_exchange.empty())
    {
        // Fill local interface traces
        for (auto &m_localInterface : m_localInterfaces)
        {
            m_localInterface->FillLocalBwdTrace(Fwd, Bwd);
        }
    }
    else
    {
        auto requestSend = comm->CreateRequest(m_exchange.size());
        auto requestRecv = comm->CreateRequest(m_exchange.size());
        for (int i = 0; i < m_exchange.size(); ++i)
        {
            m_exchange[i]->SendFwdTrace(requestSend, requestRecv, i, Fwd);
        }

        // Fill local interface traces
        for (auto &m_localInterface : m_localInterfaces)
        {
            m_localInterface->FillLocalBwdTrace(Fwd, Bwd);
        }

        comm->WaitAll(requestSend);
        comm->WaitAll(requestRecv);

        // Fill communicated interface traces
        for (auto &i : m_exchange)
        {
            i->FillRankBwdTraceExchange(Bwd);
        }
    }
}

void InterfaceTrace::FillLocalBwdTrace(Array<OneD, NekDouble> &Fwd,
                                       Array<OneD, NekDouble> &Bwd)
{
    // If flagged then fill trace from local coords
    if (m_checkLocal)
    {
        for (auto &foundLocCoord : m_foundLocalCoords)
        {
            int traceId = m_trace->GetElmtToExpId(foundLocCoord.second.first);
            Array<OneD, NekDouble> locCoord = foundLocCoord.second.second;

            Array<OneD, NekDouble> edgePhys =
                Fwd + m_trace->GetPhys_Offset(traceId);

            Bwd[foundLocCoord.first] =
                m_trace->GetExp(traceId)->StdPhysEvaluate(locCoord, edgePhys);
        }
    }
}

void InterfaceExchange::FillRankBwdTraceExchange(Array<OneD, NekDouble> &Bwd)
{
    int cnt = 0;
    for (int i = 0; i < m_interfaceTraces.size(); ++i)
    {
        Array<OneD, NekDouble> traceTmp(m_sendSize[m_rank][i] / 3, 0.0);
        for (int j = 0; j < m_sendSize[m_rank][i] / 3; ++j, ++cnt)
        {
            traceTmp[j] = m_recvTrace[cnt];
        }

        m_interfaceTraces[i]->FillRankBwdTrace(traceTmp, Bwd);
    }
}

void InterfaceTrace::FillRankBwdTrace(Array<OneD, NekDouble> &trace,
                                      Array<OneD, NekDouble> &Bwd)
{
    for (int i = 0; i < m_mapMissingCoordToTrace.size(); ++i)
    {
        if (!std::isnan(trace[i]))
        {
            Bwd[m_mapMissingCoordToTrace[i]] = trace[i];
        }
    }
}

void InterfaceExchange::SendFwdTrace(
    LibUtilities::CommRequestSharedPtr &requestSend,
    LibUtilities::CommRequestSharedPtr &requestRecv, int requestNum,
    Array<OneD, NekDouble> &Fwd)
{
    m_recvTrace =
        Array<OneD, NekDouble>(m_totSendSize[m_rank] / 3, std::nan(""));
    m_sendTrace =
        Array<OneD, NekDouble>(m_totRecvSize[m_rank] / 3, std::nan(""));

    for (auto &i : m_foundRankCoords[m_rank])
    {
        int traceId                     = m_trace->GetElmtToExpId(i.second.first);
        Array<OneD, NekDouble> locCoord = i.second.second;

        Array<OneD, NekDouble> edgePhys =
            Fwd + m_trace->GetPhys_Offset(traceId);

        m_sendTrace[i.first] =
            m_trace->GetExp(traceId)->StdPhysEvaluate(locCoord, edgePhys);
    }

    m_comm->Isend(m_rank, m_sendTrace, m_sendTrace.size(), requestSend,
                  requestNum);
    m_comm->Irecv(m_rank, m_recvTrace, m_recvTrace.size(), requestRecv,
                  requestNum);
}

/**
 * Check coords in m_recv from other rank to see if present on this rank, and
 * then populates m_foundRankCoords if found.
 * - Loops over all interface traces present in the exchange, i.e. every
 *   interface needed on this rank-to-rank communication object
 * - Then loops over all the quadrature points missing from the other side in
 *   that interface using global coordinates
 * - Loop over each local edge and use MinMaxCheck on the point as a fast
 *   comparison on whether the point lies on the edge
 * - If true then use the more expensive FindDistance function to find the
 *   closest local coordinate of the global 'missing' point on the local edge
 *   and the distance from the edge that corresponds to
 * - If distance is less than 5e-5 save found coordinate in m_foundRankCoords
 */

void InterfaceExchange::CalcRankDistances()
{
    Array<OneD, int> disp(m_recvSize[m_rank].size() + 1, 0.0);
    std::partial_sum(m_recvSize[m_rank].begin(), m_recvSize[m_rank].end(),
                     &disp[1]);

    for (int i = 0; i < m_interfaceTraces.size(); ++i)
    {
        auto localEdge = m_interfaceTraces[i]->GetInterface()->GetEdge();

        for (int j = disp[i]; j < disp[i + 1]; j += 3)
        {
            Array<OneD, NekDouble> foundLocCoord;
            Array<OneD, NekDouble> xs(3);
            xs[0] = m_recv[j];
            xs[1] = m_recv[j + 1];
            xs[2] = m_recv[j + 2];

            for (auto &edge : localEdge)
            {
                // First check if inside the edge bounding box
                if (!edge.second->MinMaxCheck(xs))
                {
                    continue;
                }

                NekDouble dist = edge.second->FindDistance(xs, foundLocCoord);

                if (dist < NekConstants::kFindDistanceMin)
                {
                    m_foundRankCoords[m_rank][j / 3] = std::make_pair(
                        edge.second->GetGlobalID(), foundLocCoord);
                    break;
                }
            }
        }
    }
}

} // namespace MultiRegions
} // namespace Nektar