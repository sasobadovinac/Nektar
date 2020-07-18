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

InterfaceMapDG::InterfaceMapDG(
    const SpatialDomains::InterfacesSharedPtr &interfaces,
    const ExpListSharedPtr &trace,
    const std::map<int, int> geomIdToTraceId)
    : m_interfaces(interfaces),
      m_trace(trace),
      m_geomIdToTraceId(geomIdToTraceId)
{
    auto comm = m_trace->GetComm();
    auto interfaceCollection = m_interfaces->GetInterfaces();

    // myIndxLR contains the info about what interface edges are present on
    // current rank with each interface no,  i,  consisting of:
    // [i] = indx
    // [i + 1] = 0 (non), = 1 (left only), = 2 (right only), = 3 (both)
    Array<OneD, int> myIndxLR(interfaceCollection.size() * 2, 0);
    std::map<int, int> myIndxLRMap;
    std::map<int, std::map<SpatialDomains::InterfaceSide, InterfaceTraceSharedPtr>> localInterfaces;
    size_t cnt = 0;
    for (const auto &interface : interfaceCollection)
    {
        myIndxLR[2 * cnt] = interface.first;

        if(!interface.second->GetLeftInterface()->IsEmpty())
        {
            myIndxLR[2 * cnt + 1] += 1;
            localInterfaces[interface.first][SpatialDomains::eLeft] = MemoryManager<InterfaceTrace>::AllocateSharedPtr(trace, interface.second->GetLeftInterface(), geomIdToTraceId);
            m_localInterfaces.emplace_back(localInterfaces[interface.first][SpatialDomains::eLeft]);
        }
        if (!interface.second->GetRightInterface()->IsEmpty())
        {
            myIndxLR[2 * cnt + 1] += 2;
            localInterfaces[interface.first][SpatialDomains::eRight] = MemoryManager<InterfaceTrace>::AllocateSharedPtr(trace, interface.second->GetRightInterface(), geomIdToTraceId);
            m_localInterfaces.emplace_back(localInterfaces[interface.first][SpatialDomains::eRight]);
        }

        myIndxLRMap[interface.first] = myIndxLR[2 * cnt + 1];
        cnt++;
    }


    //Send num of interfaces size so all partitions can prepare buffers
    int nRanks = comm->GetSize();
    Array<OneD, int> rankNumInterfaces(nRanks);
    Array<OneD, int> localInterfaceSize(1, myIndxLR.size());
    comm->AllGather(localInterfaceSize, rankNumInterfaces);

    Array<OneD, int> rankLocalInterfaceDisp(nRanks, 0);
    for (size_t i = 1; i < nRanks; ++i)
    {
        rankLocalInterfaceDisp[i] = rankLocalInterfaceDisp[i - 1] + rankNumInterfaces[i - 1];
    }

    Array<OneD, int> rankLocalInterfaceIds(
        std::accumulate(rankNumInterfaces.begin(), rankNumInterfaces.end(), 0), 0);

    // Send all interface IDs to all partitions
    comm->AllGatherv(myIndxLR, rankLocalInterfaceIds, rankNumInterfaces,rankLocalInterfaceDisp);

    // Find what interface Ids match with other ranks, then check if opposite edge
    std::map<int, std::vector<InterfaceTraceSharedPtr>> oppRankSharedInterface; // Map of rank to vector of interface traces

    size_t myRank = comm->GetRank();
    for (size_t i = 0; i < nRanks; ++i)
    {
        for (size_t j = 0; j < rankNumInterfaces[ i] / 2; ++j)
        {
            int otherId =
                rankLocalInterfaceIds[rankLocalInterfaceDisp[i] + 2 * j];
            int otherCode =
                rankLocalInterfaceIds[rankLocalInterfaceDisp[i] + 2 * j + 1];

            InterfaceExchangeSharedPtr exchange;
            if (myIndxLRMap.find(otherId) != myIndxLRMap.end())
            {
                if (i == myRank)
                {
                    //If contains opposite edge locally then set true
                    if (otherCode == 3)
                    {
                        localInterfaces[otherId][SpatialDomains::eLeft]->SetCheckLocal(true);
                        localInterfaces[otherId][SpatialDomains::eRight]->SetCheckLocal(true);
                    }

                    continue;
                }

                // Interface opposite ranks
                int myCode = myIndxLRMap[otherId];
                if ((myCode == 1 && otherCode == 2) ||
                    (myCode == 1 && otherCode == 3) ||
                    (myCode == 3 && otherCode == 2))
                {
                    oppRankSharedInterface[i].emplace_back(localInterfaces[otherId][SpatialDomains::eLeft]);
                }
                else if ((myCode == 2 && otherCode == 1) ||
                         (myCode == 2 && otherCode == 3) ||
                         (myCode == 3 && otherCode == 1))
                {
                    oppRankSharedInterface[i].emplace_back(localInterfaces[otherId][SpatialDomains::eRight]);
                }
                else if (myCode == 3 && otherCode == 3)
                {
                    oppRankSharedInterface[i].emplace_back(localInterfaces[otherId][SpatialDomains::eLeft]);
                    oppRankSharedInterface[i].emplace_back(localInterfaces[otherId][SpatialDomains::eRight]);
                }
            }
        }
    }

    // Create individual interface exchange objects (each object is rank -> rank)
    // and contains a vector of interfaceTrace objects
    for (auto &rank : oppRankSharedInterface)
    {
        m_exchange.emplace_back(MemoryManager<InterfaceExchange>::AllocateSharedPtr(m_trace, comm, rank));
    }

    for (auto &interfaceTrace : m_localInterfaces)
    {
        interfaceTrace->CalcLocalMissing();
    }

    auto request = comm->CreateRequest(2 * m_exchange.size());
    for (int i = 0; i < m_exchange.size(); ++i)
    {
        m_exchange[i]->RankFillSizes(request, i);
    }
    comm->WaitAll(request);

    for (int i = 0; i < m_exchange.size(); ++i)
    {
        m_exchange[i]->SendMissing(request, i);
    }
    comm->WaitAll(request);

    for (int i = 0; i < m_exchange.size(); ++i)
    {
        m_exchange[i]->CalcRankDistances();
    }
}

void InterfaceMapDG::ExchangeTrace(Array<OneD, NekDouble> &Fwd, Array<OneD, NekDouble> &Bwd)
{
    auto comm = m_trace->GetComm();
    auto request = comm->CreateRequest(2 * m_exchange.size());
    for (int i = 0; i < m_exchange.size(); ++i)
    {
        m_exchange[i]->SendFwdTrace(request, i, Fwd);
    }

    // Fill local interface traces
    for (int i = 0; i < m_localInterfaces.size(); ++i)
    {
        m_localInterfaces[i]->FillLocalBwdTrace(Fwd, Bwd);
    }

    comm->WaitAll(request);

    for (int i = 0; i < m_exchange.size(); ++i)
    {
        m_exchange[i]->FillRankBwdTrace(request, i, Bwd);
    }
    comm->WaitAll(request);
}

void InterfaceTrace::FillLocalBwdTrace(Array<OneD, NekDouble> &Fwd, Array<OneD, NekDouble> &Bwd)
{
    auto childEdgeIds = m_interfaceBase->GetEdgeIds();

    // If flagged then fill trace from local coords
    if(m_checkLocal)
    {
        int cnt = 0;
        for (auto childId : childEdgeIds)
        {
            auto childElmt = m_trace->GetExp(m_geomIdToTraceId.at(childId));
            int nq = childElmt->GetTotPoints();
            int offset = m_trace->GetPhys_Offset(m_geomIdToTraceId.at(childId));

            for (int i = 0; i < nq; ++i, ++cnt)
            {
                if(m_foundLocalCoords.find(cnt) != m_foundLocalCoords.end())
                {
                    Array<OneD, NekDouble> edgePhys = Fwd + m_trace->GetPhys_Offset(m_foundLocalCoords[cnt].first);
                    Array<OneD, NekDouble> foundPointArray(1, m_foundLocalCoords[cnt].second);
                    Bwd[offset + i] = m_trace->GetExp(m_foundLocalCoords[cnt].first)->StdPhysEvaluate(foundPointArray, edgePhys);
                }
            }
        }
    }
}

void InterfaceExchange::FillRankBwdTrace(LibUtilities::CommRequestSharedPtr request, int requestNum, Array<OneD, NekDouble> &Bwd)
{
    boost::ignore_unused(requestNum, request, Bwd);
    /*int cnt = 0;
    for (int i = 0; i < m_interfaces.size(); ++i)
    {
        Array<OneD, NekDouble> traceTmp(m_sendSize[i]/3);
        for (int j = 0; j < m_sendSize[i]/3; ++j)
        {
            traceTmp[j] = m_recvTrace[cnt];
        }

        auto childIds = m_interfaces[i]->GetInterface()->GetEdgeIds();
        for (auto id : childIds)
        {

        }

        auto missing = m_interfaces[i]->GetMissingCoords();
        for (int j = 0; j < missing.size(); ++j)
        {
            auto coord = missing[j];
            for (int k = 0; k < 3; ++k, ++cnt)
            {
                m_send[cnt] = coord[k];
            }
        }
    }*/
}

void InterfaceExchange::SendFwdTrace(LibUtilities::CommRequestSharedPtr request, int requestNum, Array<OneD, NekDouble> &Fwd)
{
    m_recvTrace = Array<OneD, NekDouble>(m_totSendSize/3);
    Array<OneD, NekDouble> sendTrace(m_totSendSize/3, std::numeric_limits<NekDouble>::epsilon());
    for (auto i : m_foundRankCoords)
    {
        Array<OneD, NekDouble> locCoord(1, i.second.second);
        int edgeId = i.second.first;

        Array<OneD, NekDouble> edgePhys =  Fwd + m_trace->GetPhys_Offset(edgeId);

        std::cout << m_comm->GetRank() << " AFTER GETPHYSOFFSET " << std::endl;

        sendTrace[i.first] = m_trace->GetExp(edgeId)->StdPhysEvaluate(locCoord, edgePhys);
    }

    m_comm->Isend(m_rank, sendTrace, m_totRecvSize, request, 2 * requestNum);
    m_comm->Irecv(m_rank, m_recvTrace, m_totSendSize, request, 2 * requestNum + 1);
}

void InterfaceExchange::CalcRankDistances()
{
    int cnt = 0;
    for (int i = 0; i < m_interfaces.size(); ++i)
    {
        auto parentEdgeIds = m_interfaces[i]->GetInterface()->GetEdgeIds();

        Array<OneD, NekDouble> recvTmp(m_recvSize[i]);
        for (int j = 0; j < m_recvSize[i]; j+=3, ++cnt)
        {
            NekDouble foundLocCoord;
            Array<OneD, NekDouble> xs(3);
            xs[0] = m_recv[j];
            xs[1] = m_recv[j + 1];
            xs[2] = m_recv[j + 2];

            for (auto id : parentEdgeIds)
            {
                SpatialDomains::SegGeomSharedPtr searchSeg = std::static_pointer_cast<SpatialDomains::SegGeom>(m_trace->GetGraph()->GetSegGeom(id));
                NekDouble dist = searchSeg->FindDistance(xs, foundLocCoord);

                if (dist < 1e-8)
                {
                    std::cout << "RANK: " << m_trace->GetComm()->GetRank() << " INTERFACE " << m_interfaces[i]->GetInterface()->GetId() << " I FOUND A RANK COORDINATE! " << xs[0] << " " << xs[1] << " " << xs[2] << std::endl;
                    m_foundRankCoords[cnt] = std::make_pair(searchSeg->GetGlobalID(), foundLocCoord);
                    break;
                }
            }
        }
    }
}

InterfaceTrace::InterfaceTrace(
    const ExpListSharedPtr &trace,
    const SpatialDomains::InterfaceBaseShPtr &interfaceBaseShPtr,
    const std::map<int, int> &geomIdToTraceId)
    : m_trace(trace),
      m_interfaceBase(interfaceBaseShPtr),
      m_geomIdToTraceId(geomIdToTraceId)
{
    // Calc total quad points
    for (auto id : m_interfaceBase->GetEdgeIds())
    {
        m_totQuadPnts += m_trace->GetExp(m_geomIdToTraceId.at(id))->GetTotPoints();
    }
}

void InterfaceTrace::CalcLocalMissing()
{
    auto childEdgeIds = m_interfaceBase->GetEdgeIds();

    // If not flagged check local then all points are missing
    if(!m_checkLocal)
    {
        int cnt = 0;
        for (auto childId : childEdgeIds)
        {
            auto childElmt = m_trace->GetExp(m_geomIdToTraceId.at(childId));
            size_t nq = childElmt->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);
            childElmt->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i, ++cnt)
            {
                Array<OneD, NekDouble> xs(3);
                xs[0] = xc[i];
                xs[1] = yc[i];
                xs[2] = zc[i];

                m_missingCoords[cnt] = xs;
            }
        }
    }
    // Otherwise we need to check to see what points the other side of the
    // interface on this rank contains
    else
    {
        auto parentEdgeIds = m_interfaceBase->GetOppInterface()->GetEdgeIds();

        int cnt = 0;
        for (auto childId : childEdgeIds)
        {
            auto childElmt = m_trace->GetExp(m_geomIdToTraceId.at(childId));
            size_t nq = childElmt->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);
            childElmt->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i, ++cnt)
            {
                bool found = false;
                NekDouble foundLocCoord;
                Array<OneD, NekDouble> xs(3);
                xs[0] = xc[i];
                xs[1] = yc[i];
                xs[2] = zc[i];

                for (auto id : parentEdgeIds)
                {
                    SpatialDomains::SegGeomSharedPtr searchSeg = std::static_pointer_cast<SpatialDomains::SegGeom>(m_trace->GetGraph()->GetSegGeom(id));
                    NekDouble dist = searchSeg->FindDistance(xs, foundLocCoord);

                    if (dist < 1e-8)
                    {
                        std::cout << "RANK: " << m_trace->GetComm()->GetRank() << " INTERFACE " << m_interfaceBase->GetId() << " I FOUND A LOCAL COORDINATE! " << xs[0] << " " << xs[1] << " " << xs[2] << std::endl;
                        found = true;
                        m_foundLocalCoords[cnt] = std::make_pair(m_geomIdToTraceId[searchSeg->GetGlobalID()], foundLocCoord);
                        break;
                    }
                }

                if (!found)
                {
                    ASSERTL0(m_trace->GetComm()->IsSerial(), "Serial interface and missing a coordinate.")
                    m_missingCoords[cnt] = xs;
                }
            }
        }
    }
}

void InterfaceExchange::RankFillSizes(LibUtilities::CommRequestSharedPtr request, int requestNum)
{
    // Get size of all interfaces missing to communicate
    m_sendSize = Array<OneD, int>(m_interfaces.size());
    m_recvSize = Array<OneD, int>(m_interfaces.size());
    for (int i = 0; i < m_interfaces.size(); ++i)
    {
        m_sendSize[i] = m_interfaces[i]->GetMissingCoords().size()*3;
    }

    m_comm->Isend(m_rank, m_sendSize, m_interfaces.size(), request, 2 * requestNum);
    m_comm->Irecv(m_rank, m_recvSize, m_interfaces.size(), request, 2 * requestNum + 1);
}

void InterfaceExchange::SendMissing(LibUtilities::CommRequestSharedPtr request, int requestNum)
{
    m_totSendSize = std::accumulate(m_sendSize.begin(), m_sendSize.end(), 0);
    m_totRecvSize = std::accumulate(m_recvSize.begin(), m_recvSize.end(), 0);
    m_send = Array<OneD, NekDouble>(m_totSendSize);
    m_recv = Array<OneD, NekDouble>(m_totRecvSize);

    int cnt = 0;
    for (int i = 0; i < m_interfaces.size(); ++i)
    {
        auto missing = m_interfaces[i]->GetMissingCoords();
        for (int j = 0; j < missing.size(); ++j)
        {
            auto coord = missing[j];
            for (int k = 0; k < 3; ++k, ++cnt)
            {
                m_send[cnt] = coord[k];
            }
        }
    }

    // Debug output
    std::ostringstream output;
    output << "MYRANK: " << m_comm->GetRank() << " TO RANK: " << m_rank << "\n";
    for (int i = 0; i < m_interfaces.size(); ++i)
    {
        output << "\tINTERFACE: " << m_interfaces[i]->GetInterface()->GetId() << " SEND SIZE: " << m_sendSize[i] << " RECV SIZE: " << m_recvSize[i] << "\n";
    }
    output << "\tTOTAL: " << " SEND SIZE: " << m_totSendSize << " RECV SIZE: " << m_totRecvSize << "\n";
    std::cout << output.str() << std::endl;

    m_comm->Isend(m_rank, m_send, m_totSendSize, request, 2 * requestNum);
    m_comm->Irecv(m_rank, m_recv, m_totRecvSize, request, 2 * requestNum + 1);
}

/*std::vector<std::tuple<NekDouble, int, int>> CalcCoordsOneWay(
    const SpatialDomains::InterfaceBaseShPtr &child,
    const SpatialDomains::InterfaceBaseShPtr &parent,
    const ExpListSharedPtr &trace,
    const SpatialDomains::MeshGraphSharedPtr &graph,
    const std::map<int, int> &geomToTrace,
    LibUtilities::CommSharedPtr &comm)
{

    int cnt = 0;
    //Vector of quad point matching coords {local coord, edge ID, rank}
    std::vector<std::tuple<NekDouble, int, int>> tmp(child->GetTotPoints());
    std::map<int, Array<OneD, NekDouble>> missingCoords;

    auto childEdge = child->GetEdgeIds();
    auto parentEdge = parent->GetEdgeIds();

    for (auto childId : childEdge)
    {
        int myRank = comm->GetRank();
        std::cout << "My rank: " << myRank << " child ID" << childId << std::endl;
        LocalRegions::ExpansionSharedPtr childElmt = trace->GetExp(geomToTrace.at(childId));
        size_t nq = childElmt->GetTotPoints();
        Array<OneD, NekDouble> xc(nq), yc(nq);
        childElmt->GetCoords(xc, yc);

        // Check local interface
        for (int i = 0; i < nq; ++i)
        {
            Array<OneD, NekDouble> xs(3);
            xs[0] = xc[i];
            xs[1] = yc[i];
            xs[2] = 0;

            bool found = false;
            int foundTraceId = -1;
            NekDouble foundLocCoord = -1;

            for (auto id : parentEdge)
            {
                SpatialDomains::SegGeomSharedPtr searchSeg =
                    std::static_pointer_cast<SpatialDomains::SegGeom>(graph->GetSegGeom(id));
                NekDouble dist = searchSeg->FindDistance(xs, foundLocCoord);

                if (dist < 1e-8)
                {
                    foundTraceId = geomToTrace.at(id);
                    //std::cout << "Found : "<< xc << " " << yc << " in trace " << foundTraceId << " loc coord " << foundLocCoord << " @ dist: " << dist << std::endl;
                    found = true;
                    break;
                }
            }

            if(!found)
            {
                missingCoords[cnt] = xs;
            }

            tmp[cnt] = std::make_tuple(foundLocCoord, foundTraceId, myRank);

            cnt++;
        }
    }

    std::cout << "Couldn't find locally: ";
    for (auto i : missingCoords)
    {
        std::cout << i.first << " @ " << i.second[0] << ", " << i.second[1] << " | ";
    }
    std::cout << std::endl;

    // Share how many unknown points to expect on each rank
    std::vector<int> oppRanks(child->GetOppRank().begin(), child->GetOppRank().end());

    Array<OneD, NekDouble> sendBuff(1, missingCoords.size());
    Array<OneD, NekDouble> recvBuff(oppRanks.size(), -1);

    LibUtilities::CommRequestSharedPtr recvRequest;
    LibUtilities::CommRequestSharedPtr sendRequest;
    for (int i = 0; i < oppRanks.size(); ++i)
    {
        comm->Irecv(oppRanks[i], recvBuff[i], 1, recvRequest, i);
    }

    for (int i = 0; i < oppRanks.size(); ++i)
    {
        comm->Isend(oppRanks[i], sendBuff, 1, sendRequest, i);
    }

    comm->WaitAll(sendRequest);
    comm->WaitAll(recvRequest);

    return tmp;
}

void InterfaceExchange::CalcLocalInterfaceCoords()
{
    for (auto &interface : m_interfaces->GetInterfaces())
    {
        size_t indx = interface.first;
        auto pair = interface.second;
        if (pair->GetCalcFlag())
        {
            auto left = pair->GetLeftInterface();
            auto right = pair->GetRightInterface();

            auto tmp = CalcCoordsOneWay(left, right, m_trace, m_graph, m_geomIdToTraceId, m_comm);
            auto tmp2 = CalcCoordsOneWay(right, left, m_trace, m_graph, m_geomIdToTraceId, m_comm);
            tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());

            //m_locCoordSegIdPair[indx] = tmp;
            pair->SetCalcFlag(false);
        }
    }
}*/

} // namespace MultiRegions
} // namespace Nektar