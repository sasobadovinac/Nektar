///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyCommDG.h
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
// Description: Local to Global DG mapping routines, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_ASSEMBLY_COMM_DG_H
#define MULTIREGIONS_ASSEMBLY_COMM_DG_H

#include <MultiRegions/ExpList.h>
#include <LibUtilities/BasicUtils/Timer.h>

namespace Nektar
{
namespace MultiRegions
{

class ExchangeMethod
{
public:
    /// Default constructor
    ExchangeMethod() = default;
    /// Default deconstructor
    virtual  ~ExchangeMethod() = default;
    /// Perform MPI comm exchange
    virtual void PerformExchange(
        const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd) = 0;
};

typedef std::shared_ptr<ExchangeMethod>  ExchangeMethodSharedPtr;

class AllToAll: public ExchangeMethod
{
public:
    /// Default constructor.
    AllToAll(
            const std::map<int, std::vector<int>> &rankSharedEdges,
            const std::map<int, std::vector<int>> &edgeToTrace,
            const int &nRanks,
            const int &maxQuad,
            const LibUtilities::CommSharedPtr &comm);

    virtual void PerformExchange(
        const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd) override;

private:
    /// List of trace map indices of the quad points to exchange
    std::vector<int> m_allEdgeIndex;
    /// Largest shared partition edge
    int m_maxCount = 0;
    /// Max number of quadrature points in an element
    int m_maxQuad = 0;
    /// Number of ranks/processes/partitions
    int m_nRanks = 0;
    LibUtilities::CommSharedPtr m_comm;
};

class AllToAllV: public ExchangeMethod
{
public:
    /// Default constructor.
    AllToAllV(
            const std::map<int, std::vector<int>> &rankSharedEdges,
            const std::map<int, std::vector<int>> &edgeToTrace,
            const int &nRanks,
            const LibUtilities::CommSharedPtr &comm);

    virtual void PerformExchange(
        const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd) override;

private:
    /// List of trace map indices of the quad points to exchange
    std::vector<int> m_allVEdgeIndex;
    /// List of counts for MPI_alltoallv
    Array<OneD, int> m_allVSendCount;
    /// List of displacements for MPI_alltoallv
    Array<OneD, int> m_allVSendDisp;
    LibUtilities::CommSharedPtr m_comm;
};

class NeighborAllToAllV: public ExchangeMethod
{
public:
    /// Default constructor.
    NeighborAllToAllV(
            const std::map<int, std::vector<int>> &rankSharedEdges,
            const std::map<int, std::vector<int>> &edgeToTrace);

    virtual void PerformExchange(
        const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd) override;

private:
    /// List of displacements
    Array<OneD, int> m_sendDisp;
    ///Distributed graph communicator
    MPI_Comm m_commGraph;
    /// List of trace map indices of the quad points to exchange
    std::vector<int> m_edgeTraceIndex;
    /// List of counts
    Array<OneD,int> m_sendCount;
};

class Pairwise: public ExchangeMethod
{
public:
    Pairwise(
            const std::map<int, std::vector<int>> &rankSharedEdges,
            const std::map<int, std::vector<int>> &edgeToTrace);

    virtual void PerformExchange(
        const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd) override;

private:
    /// List of partition to trace map indices of the quad points to exchange
    std::vector<std::pair<int, std::vector<int>>> m_vecPairPartitionTrace;
    /// Total quadrature points to send/recv
    int m_totSends = 0;
    /// List of displacements
    Array<OneD,int> m_sendDisp;
};

class AssemblyCommDG
{
    public:
        /// Default deconstructor
        MULTI_REGIONS_EXPORT ~AssemblyCommDG() = default;
        // Constructor for MPI communication methods

        MULTI_REGIONS_EXPORT AssemblyCommDG(
            const ExpList &locExp, const ExpListSharedPtr &trace,
            const Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> > elmtToTrace,
            const Array<OneD, const ExpListSharedPtr> &bndCondExp,
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndCond,
            const PeriodicMap &perMap);

        /// Main function that performs MPI exchange using timed fastest method
        MULTI_REGIONS_EXPORT inline void PerformExchange(
            const Array<OneD, double> &testFwd, Array<OneD, double> &testBwd)
        {
            m_exchange->PerformExchange(testFwd, testBwd);
        }

    private:
        ExchangeMethodSharedPtr m_exchange;
        /// Max number of quadrature points in an element
        int m_maxQuad = 0;
        /// Number of ranks/processes/partitions
        int m_nRanks = 0;
        /// Map of process to shared edge IDs
        std::map<int, std::vector<int>> m_rankSharedEdges;
        /// Map of edge ID to quad point trace indices
        std::map<int, std::vector<int>> m_edgeToTrace;

        /// Initalises the structure for the MPI communication
        void InitialiseStructure(
            const ExpList &locExp,
            const ExpListSharedPtr &trace,
            const Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> > elmtToTrace,
            const Array<OneD, const ExpListSharedPtr> &bndCondExp,
            const Array< OneD, const SpatialDomains::BoundaryConditionShPtr> &bndCond,
            const PeriodicMap &perMap,
            const LibUtilities::CommSharedPtr &comm);

        /// Timing of the MPI exchange method f
        static std::tuple<NekDouble, NekDouble, NekDouble> Timing(
            const LibUtilities::CommSharedPtr &comm,
            const int &count,
            const int &num,
            ExchangeMethodSharedPtr f);
};

typedef std::shared_ptr<AssemblyCommDG>  AssemblyCommDGSharedPtr;

} // end of namespace
} // end of namespace

#endif