///////////////////////////////////////////////////////////////////////////////
//
// File: CommSerial.h
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
// Description: CommSerial header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_UTILITIES_COMMSERIAL_H
#define NEKTAR_LIB_UTILITIES_COMMSERIAL_H

#include <string>

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar
{
namespace LibUtilities
{
// Forward declarations
class CommSerial;

/// Pointer to a Communicator object.
typedef std::shared_ptr<CommSerial> CommSerialSharedPtr;

/// A global linear system.
class CommSerial : public Comm
{
public:
    /// Creates an instance of this class
    LIB_UTILITIES_EXPORT static CommSharedPtr create(int narg, char *arg[])
    {
        return MemoryManager<CommSerial>::AllocateSharedPtr(narg, arg);
    }

    /// Name of class
    LIB_UTILITIES_EXPORT static std::string className;

    LIB_UTILITIES_EXPORT CommSerial(int argc, char *argv[]);
    LIB_UTILITIES_EXPORT virtual ~CommSerial() override;

protected:
    LIB_UTILITIES_EXPORT virtual void v_Finalise() override final;
    LIB_UTILITIES_EXPORT virtual int v_GetRank() override;
    LIB_UTILITIES_EXPORT virtual bool v_TreatAsRankZero() override;
    LIB_UTILITIES_EXPORT virtual bool v_IsSerial() override;
    LIB_UTILITIES_EXPORT virtual std::tuple<int, int, int> v_GetVersion()
        override final;

    LIB_UTILITIES_EXPORT virtual void v_Block() override final;
    LIB_UTILITIES_EXPORT virtual NekDouble v_Wtime() override final;
    LIB_UTILITIES_EXPORT virtual void v_Send(void *buf, int count,
                                             CommDataType dt,
                                             int dest) override final;
    LIB_UTILITIES_EXPORT virtual void v_Recv(void *buf, int count,
                                             CommDataType dt,
                                             int source) override final;
    LIB_UTILITIES_EXPORT virtual void v_SendRecv(void *sendbuf, int sendcount,
                                                 CommDataType sendtype,
                                                 int dest, void *recvbuf,
                                                 int recvcount,
                                                 CommDataType recvtype,
                                                 int source) override final;
    LIB_UTILITIES_EXPORT virtual void v_AllReduce(
        void *buf, int count, CommDataType dt,
        enum ReduceOperator pOp) override final;
    LIB_UTILITIES_EXPORT virtual void v_AlltoAll(
        void *sendbuf, int sendcount, CommDataType sendtype, void *recvbuf,
        int recvcount, CommDataType recvtype) override final;
    LIB_UTILITIES_EXPORT virtual void v_AlltoAllv(
        void *sendbuf, int sendcounts[], int sensdispls[],
        CommDataType sendtype, void *recvbuf, int recvcounts[], int rdispls[],
        CommDataType recvtype) override final;
    LIB_UTILITIES_EXPORT virtual void v_AllGather(
        void *sendbuf, int sendcount, CommDataType sendtype, void *recvbuf,
        int recvcount, CommDataType recvtype) override final;
    LIB_UTILITIES_EXPORT virtual void v_AllGatherv(
        void *sendbuf, int sendcount, CommDataType sendtype, void *recvbuf,
        int recvcounts[], int rdispls[], CommDataType recvtype) override final;
    LIB_UTILITIES_EXPORT virtual void v_AllGatherv(
        void *recvbuf, int recvcounts[], int rdispls[],
        CommDataType recvtype) override final;
    LIB_UTILITIES_EXPORT virtual void v_Bcast(void *buffer, int count,
                                              CommDataType dt,
                                              int root) override final;
    LIB_UTILITIES_EXPORT virtual void v_Gather(void *sendbuf, int sendcount,
                                               CommDataType sendtype,
                                               void *recvbuf, int recvcount,
                                               CommDataType recvtype,
                                               int root) override final;
    LIB_UTILITIES_EXPORT virtual void v_Scatter(void *sendbuf, int sendcount,
                                                CommDataType sendtype,
                                                void *recvbuf, int recvcount,
                                                CommDataType recvtype,
                                                int root) override final;

    LIB_UTILITIES_EXPORT virtual void v_DistGraphCreateAdjacent(
        int indegree, const int sources[], const int sourceweights[],
        int reorder) override final;

    LIB_UTILITIES_EXPORT virtual void v_NeighborAlltoAllv(
        void *sendbuf, int sendcounts[], int sdispls[], CommDataType sendtype,
        void *recvbuf, int recvcounts[], int rdispls[],
        CommDataType recvtype) override final;

    LIB_UTILITIES_EXPORT virtual void v_Irsend(void *buf, int count,
                                               CommDataType dt, int dest,
                                               CommRequestSharedPtr request,
                                               int loc) override final;

    LIB_UTILITIES_EXPORT virtual void v_Isend(void *buf, int count,
                                              CommDataType dt, int dest,
                                              CommRequestSharedPtr request,
                                              int loc) final;

    LIB_UTILITIES_EXPORT virtual void v_SendInit(void *buf, int count,
                                                 CommDataType dt, int dest,
                                                 CommRequestSharedPtr request,
                                                 int loc) override final;

    LIB_UTILITIES_EXPORT virtual void v_Irecv(void *buf, int count,
                                              CommDataType dt, int source,
                                              CommRequestSharedPtr request,
                                              int loc) override final;

    LIB_UTILITIES_EXPORT virtual void v_RecvInit(void *buf, int count,
                                                 CommDataType dt, int source,
                                                 CommRequestSharedPtr request,
                                                 int loc) override final;

    LIB_UTILITIES_EXPORT virtual void v_StartAll(
        CommRequestSharedPtr request) override final;
    LIB_UTILITIES_EXPORT virtual void v_WaitAll(
        CommRequestSharedPtr request) override final;
    LIB_UTILITIES_EXPORT virtual CommRequestSharedPtr v_CreateRequest(
        int num) override final;

    LIB_UTILITIES_EXPORT virtual void v_SplitComm(int pRows, int pColumns,
                                                  int pTime) override;
    LIB_UTILITIES_EXPORT virtual CommSharedPtr v_CommCreateIf(
        int flag) override final;
};
} // namespace LibUtilities
} // namespace Nektar

#endif
