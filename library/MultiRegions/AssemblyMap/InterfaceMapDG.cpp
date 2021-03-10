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

#include "InterfaceMapDG.h"

namespace Nektar
{
namespace MultiRegions
{

InterfaceTrace::InterfaceTrace(
    const ExpListSharedPtr &trace,
    const SpatialDomains::InterfaceBaseShPtr &interfaceBaseShPtr,
    const std::map<int, int> &geomIdToTraceId)
    : m_trace(trace), m_interfaceBase(interfaceBaseShPtr),
      m_geomIdToTraceId(geomIdToTraceId)
{
    // Calc total quad points
    for (auto id : m_interfaceBase->GetEdgeIds())
    {
        m_totQuadPnts +=
            m_trace->GetExp(m_geomIdToTraceId.at(id))->GetTotPoints();
    }
}

InterfaceMapDG::InterfaceMapDG(
    const SpatialDomains::InterfacesSharedPtr &interfaces,
    const ExpListSharedPtr &trace, const std::map<int, int> geomIdToTraceId)
    : m_interfaces(interfaces), m_trace(trace),
      m_geomIdToTraceId(geomIdToTraceId)
{
    auto comm                = m_trace->GetComm();
    auto interfaceCollection = m_interfaces->GetInterfaces();

    // myIndxLR contains the info about what interface edges are present on
    // current rank with each interface no,  i,  consisting of:
    // [i] = indx
    // [i + 1] = 0 (non), = 1 (left only), = 2 (right only), = 3 (both)
    std::map<int, int> myIndxLRMap;
    std::map<int,
             std::map<SpatialDomains::InterfaceSide, InterfaceTraceSharedPtr>>
        localInterfaces;
    Array<OneD, int> indxToInterfaceID(interfaceCollection.size());
    size_t cnt = 0;
    for (const auto &interface : interfaceCollection)
    {
        indxToInterfaceID[cnt]       = interface.first;
        myIndxLRMap[interface.first] = 0;

        if (!interface.second->GetLeftInterface()->IsEmpty())
        {
            myIndxLRMap[interface.first] += 1;
            localInterfaces[interface.first][SpatialDomains::eLeft] =
                MemoryManager<InterfaceTrace>::AllocateSharedPtr(
                    trace, interface.second->GetLeftInterface(),
                    geomIdToTraceId);
            m_localInterfaces.emplace_back(
                localInterfaces[interface.first][SpatialDomains::eLeft]);
        }
        if (!interface.second->GetRightInterface()->IsEmpty())
        {
            myIndxLRMap[interface.first] += 2;
            localInterfaces[interface.first][SpatialDomains::eRight] =
                MemoryManager<InterfaceTrace>::AllocateSharedPtr(
                    trace, interface.second->GetRightInterface(),
                    geomIdToTraceId);
            m_localInterfaces.emplace_back(
                localInterfaces[interface.first][SpatialDomains::eRight]);
        }

        cnt++;
    }

    // DEBUG COMMENTS
    if (false) // Set if verbose/debug mode? to output rank interface information
    {
        if (comm->GetRank() == 0)
        {
            std::cout << "\n-----------------------------\n";
            std::cout << "Rank \t" << "Interfaces\n";
            std::cout << "-----------------------------" << std::endl;
        }

        comm->Block();

        std::array<std::string, 4> edgeName = {"n", "l", "r", "b"};
        for (int rank = 0; rank < comm->GetSize(); ++rank)
        {
            if (comm->GetRank() == rank)
            {
                std::cout << rank << "\t";
                for (auto set : myIndxLRMap)
                {
                    std::cout << set.first << "[" << edgeName[set.second]
                              << "],\t";
                }
                std::cout << "\b\b  " << std::endl;
            }

            comm->Block();
        }

        if (comm->GetRank() == 0)
        {
            std::cout << "-----------------------------" << std::endl;
        }
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

    // Find what interface Ids match with other ranks, then check if opposite
    // edge
    std::map<int, std::vector<InterfaceTraceSharedPtr>>
        oppRankSharedInterface; // Map of rank to vector of interface traces

    size_t myRank     = comm->GetRank();
    int numInterfaces = interfaceCollection.size();
    for (size_t i = 0; i < nRanks; ++i)
    {
        for (size_t j = 0; j < numInterfaces; ++j)
        {
            int otherId   = indxToInterfaceID[j];
            int otherCode = rankLocalInterfaceIds[i * numInterfaces + j];

            InterfaceExchangeSharedPtr exchange;
            if (i == myRank)
            {
                // If contains opposite edge locally then set true
                if (otherCode == 3)
                {
                    localInterfaces[otherId][SpatialDomains::eLeft]
                        ->SetCheckLocal(true);
                    localInterfaces[otherId][SpatialDomains::eRight]
                        ->SetCheckLocal(true);
                }

                continue;
            }

            // Interface opposite ranks
            int myCode = myIndxLRMap[otherId];
            if ((myCode == 1 && otherCode == 2) ||
                (myCode == 1 && otherCode == 3) ||
                (myCode == 3 && otherCode == 2))
            {
                oppRankSharedInterface[i].emplace_back(
                    localInterfaces[otherId][SpatialDomains::eLeft]);
            }
            else if ((myCode == 2 && otherCode == 1) ||
                     (myCode == 2 && otherCode == 3) ||
                     (myCode == 3 && otherCode == 1))
            {
                oppRankSharedInterface[i].emplace_back(
                    localInterfaces[otherId][SpatialDomains::eRight]);
            }
            else if (myCode == 3 && otherCode == 3)
            {
                oppRankSharedInterface[i].emplace_back(
                    localInterfaces[otherId][SpatialDomains::eLeft]);
                oppRankSharedInterface[i].emplace_back(
                    localInterfaces[otherId][SpatialDomains::eRight]);
            }
        }
    }

    // Create individual interface exchange objects (each object is rank ->
    // rank) and contains a vector of interfaceTrace objects
    for (auto &rank : oppRankSharedInterface)
    {
        m_exchange.emplace_back(
            MemoryManager<InterfaceExchange>::AllocateSharedPtr(
                m_trace, comm, rank, geomIdToTraceId));
    }
}

void InterfaceMapDG::ExchangeCoords()
{
    auto comm = m_trace->GetComm();

    int skipExchange = 0;
    for (auto &interfaceTrace : m_localInterfaces)
    {
        if (interfaceTrace->GetInterface()->GetMoved())
        {
            interfaceTrace->CalcLocalMissing();
            skipExchange += 1;
        }
    }

    // If parallel communication is needed
    if (!m_exchange.empty() && skipExchange > 0)
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

    for (auto &interfaceTrace : m_localInterfaces)
    {
        interfaceTrace->GetInterface()->GetMoved() = false;
    }
}

/**
 * Calculates what coordinates on the interface are missing; i.e. aren't located
 * in an edge from the other side of the interface on this partition. These are
 * stored in the vector m_missingCoords, and another vector of the same index
 * layout contains the location of that coordinate in the trace i.e. the
 * 'offset + i' value.
 */
void InterfaceTrace::CalcLocalMissing()
{
    // Nuke old missing/found
    m_missingCoords.clear();
    m_mapMissingCoordToTrace.clear();
    m_foundLocalCoords.clear();
    m_mapFoundCoordToTrace.clear();

    auto childEdgeIds = m_interfaceBase->GetEdgeIds();

    // If not flagged 'check local' then all points are missing
    if (!m_checkLocal)
    {
        int cnt = 0;
        for (auto childId : childEdgeIds)
        {
            auto childElmt = m_trace->GetExp(m_geomIdToTraceId.at(childId));
            size_t nq      = childElmt->GetTotPoints();
            Array<OneD, NekDouble> xc(nq, 0.0), yc(nq, 0.0), zc(nq, 0.0);
            childElmt->GetCoords(xc, yc, zc);
            int offset = m_trace->GetPhys_Offset(m_geomIdToTraceId.at(childId));

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
    // interface on this rank contains
    else
    {
        auto parentEdge = m_interfaceBase->GetOppInterface()->GetEdge();

        int cnt = 0;
        for (auto childId : childEdgeIds)
        {
            auto childElmt = m_trace->GetExp(m_geomIdToTraceId.at(childId));
            size_t nq      = childElmt->GetTotPoints();
            Array<OneD, NekDouble> xc(nq, 0.0), yc(nq, 0.0), zc(nq, 0.0);
            childElmt->GetCoords(xc, yc, zc);
            int offset = m_trace->GetPhys_Offset(m_geomIdToTraceId.at(childId));

            for (int i = 0; i < nq; ++i, ++cnt)
            {
                bool found = false;
                Array<OneD, NekDouble> foundLocCoord;
                Array<OneD, NekDouble> xs(3);
                xs[0] = xc[i];
                xs[1] = yc[i];
                xs[2] = zc[i];

                for (auto &edge : parentEdge)
                {
                    NekDouble dist =
                        edge.second->FindDistance(xs, foundLocCoord);

                    if (dist < 1e-8)
                    {
                        found = true;
                        m_foundLocalCoords.emplace_back(
                            edge.second->GetGlobalID(), foundLocCoord);
                        m_mapFoundCoordToTrace.emplace_back(offset + i);
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
}

/**
 * Communicates with other ranks how many missing coordinates for each interface
 * to expect
 */
void InterfaceExchange::RankFillSizes(
    LibUtilities::CommRequestSharedPtr &requestSend,
    LibUtilities::CommRequestSharedPtr &requestRecv, int requestNum)
{
    // Get size of all interfaces missing and moved to communicate
    int recalcSize = 0;
    for (auto &localInterface : m_interfaces)
    {
        recalcSize += localInterface->GetInterface()->GetMoved();
    }

    m_sendSize = Array<OneD, int>(recalcSize);
    m_recvSize = Array<OneD, int>(recalcSize);

    int cnt = 0;
    for (auto &m_interface : m_interfaces)
    {
        if (m_interface->GetInterface()->GetMoved())
        {
            m_sendSize[cnt++] = m_interface->GetMissingCoords().size() * 3;
        }
    }

    m_comm->Isend(m_rank, m_sendSize, m_interfaces.size(), requestSend,
                  requestNum);
    m_comm->Irecv(m_rank, m_recvSize, m_interfaces.size(), requestRecv,
                  requestNum);
}

/**
 * Sends/receives the missing coordinates to/from other ranks
 */
void InterfaceExchange::SendMissing(
    LibUtilities::CommRequestSharedPtr &requestSend,
    LibUtilities::CommRequestSharedPtr &requestRecv, int requestNum)
{
    m_totSendSize = std::accumulate(m_sendSize.begin(), m_sendSize.end(), 0);
    m_totRecvSize = std::accumulate(m_recvSize.begin(), m_recvSize.end(), 0);
    m_send        = Array<OneD, NekDouble>(m_totSendSize);
    m_recv        = Array<OneD, NekDouble>(m_totRecvSize);

    int cnt = 0;
    for (auto &m_interface : m_interfaces)
    {
        if (m_interface->GetInterface()->GetMoved())
        {
            auto missing = m_interface->GetMissingCoords();
            for (auto coord : missing)
            {
                for (int k = 0; k < 3; ++k, ++cnt)
                {
                    m_send[cnt] = coord[k];
                }
            }
        }
    }

    m_comm->Isend(m_rank, m_send, m_totSendSize, requestSend, requestNum);
    m_comm->Irecv(m_rank, m_recv, m_totRecvSize, requestRecv, requestNum);
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
    ExchangeCoords();
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

    // LDG needs a consistent flux definition so swap Fwd/Bwd for right-hand side
    // of the interface
    for (auto &localInterface : m_localInterfaces)
    {
        if (localInterface->GetInterface()->GetSide() == SpatialDomains::eRight && false)
        {
            localInterface->SwapFwdBwdTrace(Fwd, Bwd);
        }
    }
}

/**
 * Swap the Fwd trace with the Bwd trace
 */
void InterfaceTrace::SwapFwdBwdTrace(Array<OneD, NekDouble> &Fwd,
                                     Array<OneD, NekDouble> &Bwd)
{
    // Flips interface edges
    for (auto &id : m_interfaceBase->GetEdgeIds())
    {
        int traceId = m_geomIdToTraceId.at(id);
        int offset = m_trace->GetPhys_Offset(traceId);

        int nq = m_trace->GetExp(traceId)->GetTotPoints();
        Array<OneD, NekDouble> tmp(nq);
        for (int i = 0; i < nq; ++i)
        {
            std::swap(Fwd[offset + i], Bwd[offset + i]);
        }
    }
}

/**
 * Fills the Bwd trace by interpolating from the Fwd for local interfaces
 */
void InterfaceTrace::FillLocalBwdTrace(Array<OneD, NekDouble> &Fwd,
                                       Array<OneD, NekDouble> &Bwd)
{
    // If flagged then fill trace from local coords
    if (m_checkLocal)
    {
        for (int i = 0; i < m_foundLocalCoords.size(); ++i)
        {
            int traceId = m_geomIdToTraceId[m_foundLocalCoords[i].first];
            Array<OneD, NekDouble> locCoord = m_foundLocalCoords[i].second;

            Array<OneD, NekDouble> edgePhys =
                Fwd + m_trace->GetPhys_Offset(traceId);
            Bwd[m_mapFoundCoordToTrace[i]] =
                m_trace->GetExp(traceId)->StdPhysEvaluate(locCoord, edgePhys);
        }
    }
}

/**
 * Loops over interfaces and partitions out the received trace from the other
 * ranks for insertion into Bwd
 */
void InterfaceExchange::FillRankBwdTraceExchange(Array<OneD, NekDouble> &Bwd)
{
    int cnt = 0;
    for (int i = 0; i < m_interfaces.size(); ++i)
    {
        Array<OneD, NekDouble> traceTmp(m_sendSize[i] / 3, 0.0);
        for (int j = 0; j < m_sendSize[i] / 3; ++j, ++cnt)
        {
            traceTmp[j] = m_recvTrace[cnt];
        }

        m_interfaces[i]->FillRankBwdTrace(traceTmp, Bwd);
    }
}

/**
 * Fills the Bwd trace from partitioned trace
 */
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

/**
 * Calculates and sends the trace to other rank from m_foundRankCoords
 */
void InterfaceExchange::SendFwdTrace(
    LibUtilities::CommRequestSharedPtr &requestSend,
    LibUtilities::CommRequestSharedPtr &requestRecv, int requestNum,
    Array<OneD, NekDouble> &Fwd)
{
    m_recvTrace = Array<OneD, NekDouble>(m_totSendSize / 3, std::nan(""));
    m_sendTrace = Array<OneD, NekDouble>(m_totRecvSize / 3, std::nan(""));

    for (auto &i : m_foundRankCoords)
    {
        Array<OneD, NekDouble> locCoord = i.second.second;
        int edgeId                      = i.second.first;

        Array<OneD, NekDouble> edgePhys =
            Fwd + m_trace->GetPhys_Offset(m_geomIdToTraceId.at(edgeId));
        m_sendTrace[i.first] = m_trace->GetExp(m_geomIdToTraceId.at(edgeId))
                                   ->StdPhysEvaluate(locCoord, edgePhys);
    }

    m_comm->Isend(m_rank, m_sendTrace, m_sendTrace.size(), requestSend,
                  requestNum);
    m_comm->Irecv(m_rank, m_recvTrace, m_recvTrace.size(), requestRecv,
                  requestNum);
}

/**
 * Check coords in m_recv from other rank to see if present on this rank, and
 * populates m_foundRankCoords
 */
// @TODO: Probably want to take this down into InterfaceTrace and use a
// GetFound() to communicate ?
void InterfaceExchange::CalcRankDistances()
{
    Array<OneD, NekDouble> disp(m_recvSize.size() + 1, 0.0);
    std::partial_sum(m_recvSize.begin(), m_recvSize.end(), &disp[1]);

    Array<OneD, int> foundNum(m_interfaces.size(), 0);
    for (int i = 0; i < m_interfaces.size(); ++i)
    {
        if (m_interfaces[i]->GetInterface()->GetMoved())
        {
            auto parentEdge = m_interfaces[i]->GetInterface()->GetEdge();

            for (int j = disp[i]; j < disp[i + 1]; j += 3)
            {
                Array<OneD, NekDouble> foundLocCoord;
                Array<OneD, NekDouble> xs(3);
                xs[0] = m_recv[j];
                xs[1] = m_recv[j + 1];
                xs[2] = m_recv[j + 2];

                for (auto &edge : parentEdge)
                {
                    NekDouble dist =
                        edge.second->FindDistance(xs, foundLocCoord);

                    if (dist < 1e-8)
                    {
                        m_foundRankCoords[j / 3] = std::make_pair(
                            edge.second->GetGlobalID(), foundLocCoord);
                        foundNum[i] += 1;
                        break;
                    }
                }
            }
        }
    }
}

} // namespace MultiRegions
} // namespace Nektar