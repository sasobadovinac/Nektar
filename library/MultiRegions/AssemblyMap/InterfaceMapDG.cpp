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

    std::cout << "RANK: " << comm->GetRank() << " BEFORE InterfaceMapDG" << std::endl;

    // myIndxLR contains the info about what interface edges are present on
    // current rank with each interface no,  i,  consisting of:
    // [i] = indx
    // [i + 1] = 0 (non), = 1 (left only), = 2 (right only), = 3 (both)
    Array<OneD, int> myIndxLR(interfaceCollection.size() * 2, 0);
    std::map<int, int> myIndxLRMap;
    size_t cnt = 0;
    for (const auto &interface : interfaceCollection)
    {
        myIndxLR[2 * cnt] = interface.first;

        if(!interface.second->GetLeftInterface()->IsEmpty())
        {
            myIndxLR[2 * cnt + 1] += 1;
            m_localInterfaces.emplace_back(interface.second->GetLeftInterface());
        }
        if (!interface.second->GetRightInterface()->IsEmpty())
        {
            myIndxLR[2 * cnt + 1] += 2;
            m_localInterfaces.emplace_back(interface.second->GetRightInterface());
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
    std::cout << "RANK: " << comm->GetRank() << " BEFORE Allgatherv for interface IDs" << std::endl;
    comm->AllGatherv(myIndxLR, rankLocalInterfaceIds, rankNumInterfaces,
                     rankLocalInterfaceDisp);
    std::cout << "RANK: " << comm->GetRank() << " AFTER Allgatherv for interface IDs" << std::endl;

    // Find what interface Ids match with other ranks, then check if opposite edge
    std::map<int, std::vector<SpatialDomains::InterfaceBaseShPtr>> oppRankSharedInterface; // Map of rank to map of interfaces to opposite edge

    size_t myRank = comm->GetRank();
    for (size_t i = 0; i < nRanks; ++i)
    {
        for (size_t j = 0; j < rankNumInterfaces[i] / 2; ++j)
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
                        interfaceCollection[otherId]->GetLeftInterface()->SetCheckLocal(true);
                        interfaceCollection[otherId]->GetRightInterface()->SetCheckLocal(true);
                    }

                    continue;
                }

                // Interface opposite ranks (could probably simplify logic
                // here but this is easy to understand
                int myCode = myIndxLRMap[otherId];
                if ((myCode == 1 && otherCode == 2) ||
                    (myCode == 1 && otherCode == 3) ||
                    (myCode == 3 && otherCode == 2))
                {
                    oppRankSharedInterface[i].emplace_back(interfaceCollection[otherId]->GetLeftInterface());
                }
                else if ((myCode == 2 && otherCode == 1) ||
                         (myCode == 2 && otherCode == 3) ||
                         (myCode == 3 && otherCode == 1))
                {
                    oppRankSharedInterface[i].emplace_back(interfaceCollection[otherId]->GetRightInterface());
                }
                else if (myCode == 3 && otherCode == 3)
                {
                    oppRankSharedInterface[i].emplace_back(interfaceCollection[otherId]->GetLeftInterface());
                    oppRankSharedInterface[i].emplace_back(interfaceCollection[otherId]->GetRightInterface());
                }
            }
        }
    }

    /*std::cout << "MY RANK: " << myRank;
    for (auto i : leftEdgeOppRanks)
    {
        std::cout << "\n\t shares interface (left) " << i.first << " with rank: ";
        for (auto j : i.second)
        {
            std::cout << j << " ";
        }
    }
    for (auto i : rightEdgeOppRanks)
    {
        std::cout << "\n\t shares interface (right) " << i.first << " with rank: ";
        for (auto j : i.second)
        {
            std::cout << j << " ";
        }
    }
    std::cout << std::endl;*/

    // Create individual interface exchange objects
    for (auto &rank : oppRankSharedInterface)
    {
        m_exchange.emplace_back(MemoryManager<InterfaceExchange>::AllocateSharedPtr(comm, rank));
    }

    // Calculate total quadrature points on each interface edge
    for (const auto &interface : m_interfaces->GetInterfaces())
    {
        int tmp = 0;
        auto leftInterface = interface.second->GetLeftInterface();
        for (auto id : leftInterface->GetEdgeIds())
        {
            tmp += m_trace->GetExp(geomIdToTraceId.at(id))->GetTotPoints();
        }
        leftInterface->SetTotPoints(tmp);

        tmp = 0;
        auto rightInterface = interface.second->GetRightInterface();
        for (auto id : leftInterface->GetEdgeIds())
        {
            tmp += m_trace->GetExp(geomIdToTraceId.at(id))->GetTotPoints();
        }
        rightInterface->SetTotPoints(tmp);
    }

    // Find local seg-coord pairs and missing coords to communicate to other ranks
    std::cout << "RANK: " << comm->GetRank() << " BEFORE CalcLocalCoords" << std::endl;
    CalcLocalCoords(trace, geomIdToTraceId);
    std::cout << "RANK: " << comm->GetRank() << " AFTER CalcLocalCoords" << std::endl;

    // Combine multiple interfaces on each rank for communication
    auto request = comm->CreateRequest(2 * m_exchange.size());
    std::cout << "RANK: " << comm->GetRank() << " BEFORE RankFillSizes()" << std::endl;
    for (int i = 0; i < m_exchange.size(); ++i)
    {
        m_exchange[i]->RankFillSizes(request, i);
    }
    std::cout << "RANK: " << comm->GetRank() << " AT WAITALL AFTER RankFillSizes()" << std::endl;

    comm->WaitAll(request);

    std::cout << "RANK: " << comm->GetRank() << " BEFORE RankCoordCalc()" << std::endl;
    for (int i = 0; i < m_exchange.size(); ++i)
    {
        m_exchange[i]->RankCoordCalc(request, i);
    }
    std::cout << "RANK: " << comm->GetRank() << " AT WAITALL AFTER RankCoordCalc()" << std::endl;
    comm->WaitAll(request);

    comm->Block();
    exit(0);
}

void InterfaceExchange::RankFillSizes(LibUtilities::CommRequestSharedPtr request, int requestNum)
{
    // Get size of all interfaces missing to communicate
    m_sendSize = Array<OneD, int>(m_interfaces.size());
    m_recvSize = Array<OneD, int>(m_interfaces.size());
    for (int i = 0; i < m_interfaces.size(); ++i)
    {
        m_sendSize[i] = m_interfaces[i]->GetMissing().size();
    }

    m_comm->Isend(m_rank, m_sendSize, m_interfaces.size(), request, 2 * requestNum);
    m_comm->Irecv(m_rank, m_recvSize, m_interfaces.size(), request, 2 * requestNum + 1);
}

void InterfaceExchange::RankCoordCalc(LibUtilities::CommRequestSharedPtr request, int requestNum)
{
    m_totSendSize = std::accumulate(m_sendSize.begin(), m_sendSize.end(), 0);
    m_totRecvSize = std::accumulate(m_recvSize.begin(), m_recvSize.end(), 0);
    m_send = Array<OneD, NekDouble>(m_totSendSize);
    m_recv = Array<OneD, NekDouble>(m_totRecvSize);

    int cnt = 0;
    for (int i = 0; i < m_interfaces.size(); ++i)
    {
        auto missing = m_interfaces[i]->GetMissing();
        for (auto coord : missing)
        {
            m_send[3 * cnt] = std::get<0>(coord.second);
            m_send[3 * cnt + 1] = std::get<1>(coord.second);
            m_send[3 * cnt + 2] = std::get<2>(coord.second);
        }
    }

    // Debug output
    boost::ignore_unused(request, requestNum);
    std::ostringstream output;
    output << "MYRANK: " << m_comm->GetRank() << " TO RANK: " << m_rank << "\n";
    for (int i = 0; i < m_interfaces.size(); ++i)
    {
        output << "\tINTERFACE: " << m_interfaces[i]->GetId() << " SEND SIZE: " << m_sendSize[i] << " RECV SIZE: " << m_recvSize[i] << "\n";
    }
    output << "\tTOTAL: " << " SEND SIZE: " << m_totSendSize << " RECV SIZE: " << m_totRecvSize << "\n";
    std::cout << output.str() << std::endl;

    m_comm->Isend(m_rank, m_send, m_totSendSize, request, 2 * requestNum);
    m_comm->Irecv(m_rank, m_recv, m_totRecvSize, request, 2 * requestNum + 1);
}

void InterfaceMapDG::CalcLocalCoords(const ExpListSharedPtr &trace, std::map<int, int> geomIdToTraceId)
{
    for (auto &interface : m_localInterfaces)
    {
        std::map<int, std::tuple<NekDouble, NekDouble, NekDouble>> missingCoords;
        std::map<int, std::pair<int, NekDouble>> foundEdgeLocalCoordPair;

        if (interface->GetCheckLocal())
        {
            auto graph = trace->GetGraph();
            Array<OneD, std::pair<int, NekDouble>> geomEdgeIdLocalCoordPair(interface->GetTotPoints());

            int cnt            = 0;
            auto childEdgeIds  = interface->GetEdgeIds();
            auto parentEdgeIds = interface->GetOppInterface()->GetEdgeIds();
            for (auto childId : childEdgeIds)
            {
                auto childElmt = trace->GetExp(geomIdToTraceId.at(childId));
                size_t nq      = childElmt->GetTotPoints();
                Array<OneD, NekDouble> xc(nq), yc(nq);
                childElmt->GetCoords(xc, yc);

                // Check local interface
                for (int i = 0; i < nq; ++i)
                {
                    Array<OneD, NekDouble> xs(3);
                    xs[0] = xc[i];
                    xs[1] = yc[i];
                    xs[2] = 0;

                    bool found              = false;
                    int foundEdgeId         = -1;
                    NekDouble foundLocCoord = -1;

                    for (auto parentId : parentEdgeIds)
                    {
                        SpatialDomains::SegGeomSharedPtr searchSeg =
                            std::static_pointer_cast<SpatialDomains::SegGeom>(
                                graph->GetSegGeom(parentId));
                        NekDouble dist =
                            searchSeg->FindDistance(xs, foundLocCoord);

                        if (dist < 1e-8)
                        {
                            foundEdgeId = parentId;
                            // std::cout << "Found : "<< xc[i] << " " << yc[i] << " in trace " << foundEdgeId << " loc coord " << foundLocCoord << " @ dist: " << dist << std::endl;
                            found = true;
                            break;
                        }
                    }

                    if (!found)
                    {
                        missingCoords[cnt] = std::make_tuple(xs[0], xs[1], xs[2]);
                    }

                    foundEdgeLocalCoordPair[cnt] = std::make_pair(foundEdgeId, foundLocCoord);
                    cnt++;
                }
            }
        }
        else
        {
            int cnt           = 0;
            auto childEdgeIds = interface->GetEdgeIds();
            for (auto childId : childEdgeIds)
            {
                auto childElmt = trace->GetExp(geomIdToTraceId.at(childId));
                size_t nq      = childElmt->GetTotPoints();
                Array<OneD, NekDouble> xc(nq), yc(nq);
                childElmt->GetCoords(xc, yc);

                for (int i = 0; i < nq; ++i)
                {
                    missingCoords[cnt] = std::make_tuple(xc[i], yc[i], 0);
                    cnt++;
                }
            }
        }

        interface->SetMissing(missingCoords);

        // Cout missing coords (debug)
        /*std::ostringstream output;
        output << "MY RANK: " << trace->GetComm()->GetRank() << " INTERFACE: " << interface->GetId() << " " << SpatialDomains::InterfaceSideStr[interface->GetSide()] << " MISSING: \n";
        for (auto i : missingCoords)
        {
            output << "\t" << i.first << "\t" << std::get<0>(i.second) << " "
                   << std::get<1>(i.second) << " " << std::get<2>(i.second)
                   << "\n";
        }
        output << std::endl;
        std::cout << output.str();*/
    }
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