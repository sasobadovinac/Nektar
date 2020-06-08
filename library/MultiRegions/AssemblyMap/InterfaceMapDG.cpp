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
    const ExpList &locExp, const ExpListSharedPtr &trace,
    const std::map<int, int> geomIdToTraceId)
{
    auto comm = trace->GetComm();
    auto interfaceCollection = interfaces->GetInterfaces();

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
        }
        if (!interface.second->GetRightInterface()->IsEmpty())
        {
            myIndxLR[2 * cnt + 1] += 2;
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
    comm->AllGatherv(myIndxLR, rankLocalInterfaceIds, rankNumInterfaces,
                     rankLocalInterfaceDisp);

    // Find what interface Ids match with other ranks, then check if opposite edge
    size_t myRank = comm->GetRank();
    for (size_t i = 0; i < nRanks; ++i)
    {
        if (i == myRank)
        {
            continue;
        }

        for (size_t j = 0; j < rankNumInterfaces[i] / 2; ++j)
        {
            int otherId =
                rankLocalInterfaceIds[rankLocalInterfaceDisp[i] + 2 * j];
            int otherCode =
                rankLocalInterfaceIds[rankLocalInterfaceDisp[i] + 2 * j + 1];

            InterfaceExchangeSharedPtr exchange;
            if (myIndxLRMap.find(otherId) != myIndxLRMap.end())
            {
                // `Set` interface opposite ranks (could probably simplify logic
                // here but this is easy to understand
                int myCode = myIndxLRMap[otherId];
                if ((myCode == 1 && otherCode == 2) ||
                    (myCode == 1 && otherCode == 3) ||
                    (myCode == 3 && otherCode == 2))
                {
                    m_exchange.emplace_back(
                        MemoryManager<InterfaceExchange>::AllocateSharedPtr(
                            interfaceCollection[otherId]->GetLeftInterface(), i));
                }
                else if ((myCode == 2 && otherCode == 1) ||
                         (myCode == 2 && otherCode == 3) ||
                         (myCode == 3 && otherCode == 1))
                {
                    m_exchange.emplace_back(
                        MemoryManager<InterfaceExchange>::AllocateSharedPtr(
                            interfaceCollection[otherId]->GetRightInterface(),
                            i));
                }
                else if (myCode == 3 && otherCode == 3)
                {
                    m_exchange.emplace_back(
                        MemoryManager<InterfaceExchange>::AllocateSharedPtr(
                            interfaceCollection[otherId]->GetLeftInterface(), i));
                    m_exchange.emplace_back(
                        MemoryManager<InterfaceExchange>::AllocateSharedPtr(
                            interfaceCollection[otherId]->GetRightInterface(),
                            i));
                }
            }
        }
    }
}

std::vector<std::tuple<NekDouble, int, int>> CalcCoordsOneWay(
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
}

} // namespace MultiRegions
} // namespace Nektar