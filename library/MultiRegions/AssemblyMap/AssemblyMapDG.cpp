///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMapDG.cpp
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
// Description: Local to Global Base Class mapping routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/HashUtils.hpp>
#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/ExpList.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/HexExp.h>
#include <LocalRegions/TetExp.h>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/PyrExp.h>
#include <LocalRegions/PointExp.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>

#include <boost/core/ignore_unused.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        AssemblyMapDG::AssemblyMapDG():
            m_numDirichletBndPhys(0)
        {
        }

        AssemblyMapDG::~AssemblyMapDG()
        {
        }

        AssemblyMapDG::AssemblyMapDG(
            const LibUtilities::SessionReaderSharedPtr                &pSession,
            const SpatialDomains::MeshGraphSharedPtr                  &graph,
            const ExpListSharedPtr                                    &trace,
            const ExpList                                             &locExp,
            const Array<OneD, const MultiRegions::ExpListSharedPtr>         &bndCondExp,
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndCond,
            const PeriodicMap                                         &periodicTrace,
            const std::string variable):
            AssemblyMap(pSession,variable)
        {
            boost::ignore_unused(graph);

            int i, j, k, cnt, eid, id, id1, gid;
            int order_e   = 0;
            int nTraceExp = trace->GetExpSize();
            int nbnd      = bndCondExp.num_elements();

            LocalRegions::ExpansionSharedPtr  exp;
            LocalRegions::ExpansionSharedPtr  bndExp;
            SpatialDomains::GeometrySharedPtr traceGeom;

            const LocalRegions::ExpansionVector expList = *(locExp.GetExp());
            int nel = expList.size();

            map<int, int> meshTraceId;

            m_signChange = true;

            // determine mapping from geometry edges to trace
            for(i = 0; i < nTraceExp; ++i)
            {
                meshTraceId[trace->GetExp(i)->GetGeom()->GetGlobalID()] = i;
            }

            // Count total number of trace elements
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                cnt += expList[i]->GetNtrace();
            }

            Array<OneD, LocalRegions::ExpansionSharedPtr> tracemap(cnt);
            m_elmtToTrace = Array<
                OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >(nel);

            // set up trace expansions links;
            for(cnt = i = 0; i < nel; ++i)
            {
                m_elmtToTrace[i] = tracemap + cnt;

                for(j = 0; j < expList[i]->GetNtrace(); ++j)
                {
                    id = expList[i]->GetGeom()->GetTid(j);

                    if(meshTraceId.count(id) > 0)
                    {
                        m_elmtToTrace[i][j] =
                            trace->GetExp(meshTraceId.find(id)->second);
                    }
                    else
                    {
                        ASSERTL0(false, "Failed to find trace map");
                    }
                }

                cnt += expList[i]->GetNtrace();
            }

            // Set up boundary mapping
            cnt = 0;
            for(i = 0; i < nbnd; ++i)
            {
                if (bndCond[i]->GetBoundaryConditionType() ==
                        SpatialDomains::ePeriodic)
                {
                    continue;
                }
                cnt += bndCondExp[i]->GetExpSize();
            }

            set<int> dirTrace;

            m_numLocalDirBndCoeffs = 0;
            m_numDirichletBndPhys  = 0;

            for(i = 0; i < bndCondExp.num_elements(); ++i)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    bndExp    = bndCondExp[i]->GetExp(j);
                    traceGeom = bndExp->GetGeom();
                    id        = traceGeom->GetGlobalID();

                    if(bndCond[i]->GetBoundaryConditionType() ==
                           SpatialDomains::eDirichlet)
                    {
                        m_numLocalDirBndCoeffs += bndExp->GetNcoeffs();
                        m_numDirichletBndPhys  += bndExp->GetTotPoints();
                        dirTrace.insert(id);
                    }
                }
            }

            // Set up integer mapping array and sign change for each degree of
            // freedom + initialise some more data members.
            m_staticCondLevel           = 0;
            m_lowestStaticCondLevel     = 0;
            m_numPatches                = nel;
            m_numLocalBndCoeffsPerPatch = Array<OneD, unsigned int>(nel);
            m_numLocalIntCoeffsPerPatch = Array<OneD, unsigned int>(nel);

            int nbndry = 0;
            for(i = 0; i < nel; ++i) // count number of elements in array
            {
                eid     = i;
                nbndry += expList[eid]->NumDGBndryCoeffs();
                m_numLocalIntCoeffsPerPatch[i] = 0;
                m_numLocalBndCoeffsPerPatch[i] =
                    (unsigned int) expList[eid]->NumDGBndryCoeffs();
            }

            m_numGlobalDirBndCoeffs = m_numLocalDirBndCoeffs;
            m_numLocalBndCoeffs     = nbndry;
            m_numLocalCoeffs        = nbndry;
            m_localToGlobalBndMap   = Array<OneD, int>       (nbndry);
            m_localToGlobalBndSign  = Array<OneD, NekDouble> (nbndry,1);

            // Set up array for potential mesh optimsation
            Array<OneD,int> traceElmtGid(nTraceExp, -1);
            int nDir = 0;
            cnt = 0;

            // We are now going to construct a graph of the mesh which can be
            // reordered depending on the type of solver we would like to use.
            typedef boost::adjacency_list<
                boost::setS, boost::vecS, boost::undirectedS> BoostGraph;

            BoostGraph boostGraphObj;
            int trace_id, trace_id1;
            int dirOffset = 0;

            // make trace trace renumbering map where first solved trace starts
            // at 0 so we can set up graph.
            for(i = 0; i < nTraceExp; ++i)
            {
                id = trace->GetExp(i)->GetGeom()->GetGlobalID();

                if (dirTrace.count(id) == 0)
                {
                    // Initial put in element ordering (starting from zero) into
                    // traceElmtGid
                    boost::add_vertex(boostGraphObj);
                    traceElmtGid[i] = cnt++;
                }
                else
                {
                    // Use existing offset for Dirichlet edges
                    traceElmtGid[i] = dirOffset;
                    dirOffset      += trace->GetExp(i)->GetNcoeffs();
                    nDir++;
                }
            }

            // Set up boost Graph
            for(i = 0; i < nel; ++i)
            {
                eid = i;

                for(j = 0; j < expList[eid]->GetNtrace(); ++j)
                {
                    // Add trace to boost graph for non-Dirichlet Boundary
                    traceGeom = m_elmtToTrace[eid][j]->GetGeom();
                    id        = traceGeom->GetGlobalID();
                    trace_id  = meshTraceId.find(id)->second;

                    if(dirTrace.count(id) == 0)
                    {
                        for(k = j+1; k < expList[eid]->GetNtrace(); ++k)
                        {
                            traceGeom = m_elmtToTrace[eid][k]->GetGeom();
                            id1       = traceGeom->GetGlobalID();
                            trace_id1 = meshTraceId.find(id1)->second;

                            if(dirTrace.count(id1) == 0)
                            {
                                boost::add_edge((size_t)traceElmtGid[trace_id],
                                                (size_t)traceElmtGid[trace_id1],
                                                boostGraphObj);
                            }
                        }
                    }
                }
            }

            int                                 nGraphVerts = nTraceExp - nDir;
            Array<OneD, int>                    perm (nGraphVerts);
            Array<OneD, int>                    iperm(nGraphVerts);
            Array<OneD, int>                    vwgts(nGraphVerts);
            BottomUpSubStructuredGraphSharedPtr bottomUpGraph;

            for(i = 0; i < nGraphVerts; ++i)
            {
                vwgts[i] = trace->GetExp(i+nDir)->GetNcoeffs();
            }

            if(nGraphVerts)
            {
                switch(m_solnType)
                {
                    case eDirectFullMatrix:
                    case eIterativeFull:
                    case eIterativeStaticCond:
                    case eXxtFullMatrix:
                    case eXxtStaticCond:
                    case ePETScFullMatrix:
                    case ePETScStaticCond:
                    {
                        NoReordering(boostGraphObj,perm,iperm);
                        break;
                    }
                    case eDirectStaticCond:
                    {
                        CuthillMckeeReordering(boostGraphObj,perm,iperm);
                        break;
                    }
                    case eDirectMultiLevelStaticCond:
                    case eIterativeMultiLevelStaticCond:
                    case eXxtMultiLevelStaticCond:
                    case ePETScMultiLevelStaticCond:
                    {
                        MultiLevelBisectionReordering(boostGraphObj,perm,iperm,
                                                      bottomUpGraph);
                        break;
                    }
                    default:
                    {
                        ASSERTL0(false,"Unrecognised solution type");
                    }
                }
            }

            // Recast the permutation so that it can be used as a map from old
            // trace ID to new trace ID
            cnt = m_numLocalDirBndCoeffs;
            for(i = 0; i < nTraceExp - nDir; ++i)
            {
                traceElmtGid[perm[i]+nDir] = cnt;
                cnt += trace->GetExp(perm[i]+nDir)->GetNcoeffs();
            }

            // Now have trace edges Gid position

            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                eid = i;
                exp = expList[eid];

                for(j = 0; j < exp->GetNtrace(); ++j)
                {
                    traceGeom = m_elmtToTrace[eid][j]->GetGeom();
                    id        = traceGeom->GetGlobalID();
                    gid       = traceElmtGid[meshTraceId.find(id)->second];

                    const int nDim = expList[eid]->GetNumBases();

                    if (nDim == 1)
                    {
                        order_e = 1;
                        m_localToGlobalBndMap[cnt] = gid;
                    }
                    else if (nDim == 2)
                    {
                        order_e = expList[eid]->GetEdgeNcoeffs(j);
                    
                        if(expList[eid]->GetEorient(j) == StdRegions::eForwards)
                        {
                            for(k = 0; k < order_e; ++k)
                            {
                                m_localToGlobalBndMap[k+cnt] = gid + k;
                            }
                        }
                        else
                        {
                            switch(m_elmtToTrace[eid][j]->GetBasisType(0))
                            {
                                case LibUtilities::eModified_A:
                                {
                                    // reverse vertex order
                                    m_localToGlobalBndMap[cnt]   = gid + 1;
                                    m_localToGlobalBndMap[cnt+1] = gid;
                                    for (k = 2; k < order_e; ++k)
                                    {
                                        m_localToGlobalBndMap[k+cnt] = gid + k;
                                    }

                                    // negate odd modes
                                    for(k = 3; k < order_e; k+=2)
                                    {
                                        m_localToGlobalBndSign[cnt+k] = -1.0;
                                    }
                                    break;
                                }
                                case LibUtilities::eGLL_Lagrange:
                                {
                                    // reverse  order
                                    for(k = 0; k < order_e; ++k)
                                    {
                                        m_localToGlobalBndMap[cnt+order_e-k-1] = gid + k;
                                    }
                                    break;
                                }
                                case LibUtilities::eGauss_Lagrange:
                                {
                                    // reverse  order
                                    for(k = 0; k < order_e; ++k)
                                    {
                                        m_localToGlobalBndMap[cnt+order_e-k-1] = gid + k;
                                    }
                                    break;
                                }
                                default:
                                {
                                    ASSERTL0(false,"Boundary type not permitted");
                                }
                            }
                        }
                    }
                    else if (nDim == 3)
                    {
                        order_e = expList[eid]->GetFaceNcoeffs(j);

                        std::map<int, int> orientMap;

                        Array<OneD, unsigned int> elmMap1 (order_e);
                        Array<OneD,          int> elmSign1(order_e);
                        Array<OneD, unsigned int> elmMap2 (order_e);
                        Array<OneD,          int> elmSign2(order_e);

                        StdRegions::Orientation fo = expList[eid]->GetForient(j);

                        // Construct mapping which will permute global IDs
                        // according to face orientations.
                        expList[eid]->GetFaceToElementMap(j,fo,elmMap1,elmSign1);
                        expList[eid]->GetFaceToElementMap(
                            j,StdRegions::eDir1FwdDir1_Dir2FwdDir2,elmMap2,elmSign2);

                        for (k = 0; k < elmMap1.num_elements(); ++k)
                        {
                            // Find the elemental co-efficient in the original
                            // mapping.
                            int idx = -1;
                            for (int l = 0; l < elmMap2.num_elements(); ++l)
                            {
                                if (elmMap1[k] == elmMap2[l])
                                {
                                    idx = l;
                                    break;
                                }
                            }

                            ASSERTL2(idx != -1, "Problem with face to element map!");
                            orientMap[k] = idx;
                        }

                        for(k = 0; k < order_e; ++k)
                        {
                            m_localToGlobalBndMap [k+cnt] = gid + orientMap[k];
                            m_localToGlobalBndSign[k+cnt] = elmSign2[orientMap[k]];
                        }
                    }

                    cnt += order_e;
                }
            }

            // set up m_bndCondCoeffsToGlobalCoeffsMap to align with map
            cnt = 0;
            for(i = 0; i < nbnd; ++i)
            {
                if (bndCond[i]->GetBoundaryConditionType() ==
                    SpatialDomains::ePeriodic)
                {
                    continue;
                }
                cnt += bndCondExp[i]->GetNcoeffs();
            }

            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD, int>(cnt);

            // Number of boundary expansions
            int nbndexp = 0, bndOffset, bndTotal = 0;
            for(cnt = i = 0; i < nbnd; ++i)
            {
                if (bndCond[i]->GetBoundaryConditionType() ==
                    SpatialDomains::ePeriodic)
                {
                    continue;
                }

                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    bndExp    = bndCondExp[i]->GetExp(j);
                    id        = bndExp->GetGeom()->GetGlobalID();
                    gid       = traceElmtGid[meshTraceId.find(id)->second];
                    bndOffset = bndCondExp[i]->GetCoeff_Offset(j) + bndTotal;

                    // Since boundary information is defined to be aligned with
                    // the geometry just use forward/forward (both coordinate
                    // directions) defintiion for gid's
                    for(k = 0; k < bndExp->GetNcoeffs(); ++k)
                    {
                        m_bndCondCoeffsToGlobalCoeffsMap[bndOffset+k] = gid + k;
                    }
                }

                nbndexp  += bndCondExp[i]->GetExpSize();
                bndTotal += bndCondExp[i]->GetNcoeffs();
            }

            m_numGlobalBndCoeffs = trace->GetNcoeffs();
            m_numGlobalCoeffs    = m_numGlobalBndCoeffs;

            CalculateBndSystemBandWidth();

            cnt = 0;
            m_bndCondTraceToGlobalTraceMap = Array<OneD, int>(nbndexp);
            for(i = 0; i < bndCondExp.num_elements(); ++i)
            {
                if (bndCond[i]->GetBoundaryConditionType() ==
                    SpatialDomains::ePeriodic)
                {
                    continue;
                }

                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    bndExp    = bndCondExp[i]->GetExp(j);
                    traceGeom = bndExp->GetGeom();
                    id        = traceGeom->GetGlobalID();
                    m_bndCondTraceToGlobalTraceMap[cnt++] =
                        meshTraceId.find(id)->second;
                }
            }

            // Now set up mapping from global coefficients to universal.
            ExpListSharedPtr tr = std::dynamic_pointer_cast<ExpList>(trace);
            SetUpUniversalDGMap   (locExp);
            SetUpUniversalTraceMap(locExp, tr, bndCondExp, bndCond, periodicTrace);

            if ((m_solnType == eDirectMultiLevelStaticCond ||
                 m_solnType == eIterativeMultiLevelStaticCond ||
                 m_solnType == eXxtMultiLevelStaticCond ||
                 m_solnType == ePETScMultiLevelStaticCond)
                && nGraphVerts)
            {
                if (m_staticCondLevel < (bottomUpGraph->GetNlevels() - 1))
                {
                    Array<OneD, int> vwgts_perm(nGraphVerts);

                    for (int i = 0; i < nGraphVerts; i++)
                    {
                        vwgts_perm[i] = vwgts[perm[i]];
                    }

                    bottomUpGraph->ExpandGraphWithVertexWeights(vwgts_perm);
                    m_nextLevelLocalToGlobalMap = MemoryManager<AssemblyMap>::
                        AllocateSharedPtr(this, bottomUpGraph);
                }
            }

            m_hash = hash_range(m_localToGlobalBndMap.begin(),
                                m_localToGlobalBndMap.end());
        }

        /**
         * Constructs a mapping between the process-local global numbering and
         * a universal numbering of the trace space expansion. The universal
         * numbering is defined by the mesh edge IDs to enforce consistency
         * across processes.
         *
         * @param       locExp  List of local elemental expansions.
         */
        void AssemblyMapDG::SetUpUniversalDGMap(const ExpList &locExp)
        {
            LocalRegions::ExpansionSharedPtr locExpansion;
            int eid       = 0;
            int cnt       = 0;
            int id        = 0;
            int order_e   = 0;
            int vGlobalId = 0;
            int maxDof    = 0;
            int dof       = 0;
            int nDim      = 0;
            int i,j,k;

            const LocalRegions::ExpansionVector &locExpVector = *(locExp.GetExp());

            // Initialise the global to universal maps.
            m_globalToUniversalBndMap = Nektar::Array<OneD, int>(m_numGlobalBndCoeffs, -1);
            m_globalToUniversalBndMapUnique = Nektar::Array<OneD, int>(m_numGlobalBndCoeffs, -1);

            // Loop over all the elements in the domain and compute max
            // DOF. Reduce across all processes to get universal maximum.
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = locExpVector[i];
                nDim = locExpansion->GetShapeDimension();

                // Loop over all edges of element i
                if (nDim == 1)
                {
                    maxDof = (1 > maxDof ? 1 : maxDof);
                }
                else if (nDim == 2)
                {
                    for (j = 0; j < locExpansion->GetNedges(); ++j)
                    {
                        dof    = locExpansion->GetEdgeNcoeffs(j);
                        maxDof = (dof > maxDof ? dof : maxDof);
                    }
                }
                else if (nDim == 3)
                {
                    for (j = 0; j < locExpansion->GetNfaces(); ++j)
                    {
                        dof    = locExpansion->GetFaceNcoeffs(j);
                        maxDof = (dof > maxDof ? dof : maxDof);
                    }
                }
            }
            m_comm->AllReduce(maxDof, LibUtilities::ReduceMax);

            // Now have trace edges Gid position
            cnt = 0;
            for(i = 0; i < locExpVector.size(); ++i)
            {
                eid = i;
                locExpansion = locExpVector[eid];
                nDim = locExpansion->GetShapeDimension();

                // Populate mapping for each edge of the element.
                if (nDim == 1)
                {
                    int nverts = locExpansion->GetNverts();
                    for(j = 0; j < nverts; ++j)
                    {
                        LocalRegions::PointExpSharedPtr locPointExp =
                            m_elmtToTrace[eid][j]->as<LocalRegions::PointExp>();
                        id = locPointExp->GetGeom()->GetGlobalID();
                        vGlobalId = m_localToGlobalBndMap[cnt+j];
                        m_globalToUniversalBndMap[vGlobalId]
                            = id * maxDof + j + 1;
                    }
                    cnt += nverts;
                }
                else if (nDim == 2)
                {
                    for(j = 0; j < locExpansion->GetNedges(); ++j)
                    {
                        LocalRegions::SegExpSharedPtr locSegExp =
                            m_elmtToTrace[eid][j]->as<LocalRegions::SegExp>();

                        id  = locSegExp->GetGeom()->GetGlobalID();
                        order_e = locExpansion->GetEdgeNcoeffs(j);

                        map<int,int> orientMap;
                        Array<OneD, unsigned int> map1(order_e), map2(order_e);
                        Array<OneD, int> sign1(order_e), sign2(order_e);

                        locExpansion->GetEdgeToElementMap(j, StdRegions::eForwards, map1, sign1);
                        locExpansion->GetEdgeToElementMap(j, locExpansion->GetEorient(j), map2, sign2);

                        for (k = 0; k < map1.num_elements(); ++k)
                        {
                            // Find the elemental co-efficient in the original
                            // mapping.
                            int idx = -1;
                            for (int l = 0; l < map2.num_elements(); ++l)
                            {
                                if (map1[k] == map2[l])
                                {
                                    idx = l;
                                    break;
                                }
                            }

                            ASSERTL2(idx != -1, "Problem with face to element map!");
                            orientMap[k] = idx;
                        }

                        for(k = 0; k < order_e; ++k)
                        {
                            vGlobalId = m_localToGlobalBndMap[k+cnt];
                            m_globalToUniversalBndMap[vGlobalId]
                                = id * maxDof + orientMap[k] + 1;
                        }
                        cnt += order_e;
                    }
                }
                else if (nDim == 3)
                {
                    for(j = 0; j < locExpansion->GetNfaces(); ++j)
                    {
                        LocalRegions::Expansion2DSharedPtr locFaceExp =
                                m_elmtToTrace[eid][j]
                                           ->as<LocalRegions::Expansion2D>();

                        id  = locFaceExp->GetGeom()->GetGlobalID();
                        order_e = locExpansion->GetFaceNcoeffs(j);

                        map<int,int> orientMap;
                        Array<OneD, unsigned int> map1(order_e), map2(order_e);
                        Array<OneD, int> sign1(order_e), sign2(order_e);

                        locExpansion->GetFaceToElementMap(j, StdRegions::eDir1FwdDir1_Dir2FwdDir2, map1, sign1);
                        locExpansion->GetFaceToElementMap(j, locExpansion->GetForient(j), map2, sign2);

                        for (k = 0; k < map1.num_elements(); ++k)
                        {
                            // Find the elemental co-efficient in the original
                            // mapping.
                            int idx = -1;
                            for (int l = 0; l < map2.num_elements(); ++l)
                            {
                                if (map1[k] == map2[l])
                                {
                                    idx = l;
                                    break;
                                }
                            }

                            ASSERTL2(idx != -1, "Problem with face to element map!");
                            orientMap[k] = idx;
                        }

                        for(k = 0; k < order_e; ++k)
                        {
                            vGlobalId = m_localToGlobalBndMap[k+cnt];
                            m_globalToUniversalBndMap[vGlobalId]
                                = id * maxDof + orientMap[k] + 1;
                        }
                        cnt += order_e;
                    }
                }
            }

            // Initialise GSlib and populate the unique map.
            Array<OneD, long> tmp(m_globalToUniversalBndMap.num_elements());
            for (i = 0; i < m_globalToUniversalBndMap.num_elements(); ++i)
            {
                tmp[i] = m_globalToUniversalBndMap[i];
            }
            m_bndGsh = m_gsh = Gs::Init(tmp, m_comm);
            Gs::Unique(tmp, m_comm);
            for (i = 0; i < m_globalToUniversalBndMap.num_elements(); ++i)
            {
                m_globalToUniversalBndMapUnique[i] = (tmp[i] >= 0 ? 1 : 0);
            }
        }

        void AssemblyMapDG::SetUpUniversalTraceMap(
            const ExpList         &locExp,
            const ExpListSharedPtr trace,
            const Array<OneD, const MultiRegions::ExpListSharedPtr>         &bndCondExp,
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>  &bndCond,
            const PeriodicMap     &perMap)
        {
            Array<OneD, int> tmp;
            LocalRegions::ExpansionSharedPtr locExpansion;
            int i;
            int maxQuad = 0, quad = 0, nDim = 0, eid = 0, offset = 0;

            const LocalRegions::ExpansionVector &locExpVector = *(locExp.GetExp());

            int nTracePhys = trace->GetTotPoints();

            // Initialise the trace to universal maps.
            m_traceToUniversalMap       =
                Nektar::Array<OneD, int>(nTracePhys, -1);
            m_traceToUniversalMapUnique =
                Nektar::Array<OneD, int>(nTracePhys, -1);

            // Assume that each element of the expansion is of the same
            // dimension.
            nDim = locExpVector[0]->GetShapeDimension();

            if (nDim == 1)
            {
                maxQuad = (1 > maxQuad ? 1 : maxQuad);
            }
            else
            {
                for (i = 0; i < trace->GetExpSize(); ++i)
                {
                    quad = trace->GetExp(i)->GetTotPoints();
                    if (quad > maxQuad)
                    {
                        maxQuad = quad;
                    }
                }
            }
            m_comm->AllReduce(maxQuad, LibUtilities::ReduceMax);

            if (nDim == 1)
            {
                for (int i = 0; i < trace->GetExpSize(); ++i)
                {
                    eid = trace->GetExp(i)->GetGeom()->GetGlobalID();
                    offset = trace->GetPhys_Offset(i);

                    // Check to see if this vert is periodic. If it is, then we
                    // need use the unique eid of the two points
                    auto it = perMap.find(eid);
                    if (perMap.count(eid) > 0)
                    {
                        PeriodicEntity ent = it->second[0];
                        if (ent.isLocal == false) // Not sure if true in 1D
                        {
                            eid = min(eid, ent.id);
                        }
                    }

                    m_traceToUniversalMap[offset] = eid*maxQuad+1;
                    m_edgeToTrace[eid].emplace_back(offset);
                }
            }
            else
            {
                for (int i = 0; i < trace->GetExpSize(); ++i)
                {
                    eid    = trace->GetExp(i)->GetGeom()->GetGlobalID();
                    offset = trace->GetPhys_Offset(i);
                    quad   = trace->GetExp(i)->GetTotPoints();

                    // Check to see if this edge is periodic. If it is, then we
                    // need to reverse the trace order of one edge only in the
                    // universal map so that the data are reversed w.r.t each
                    // other. We do this by using the minimum of the two IDs.
                    auto it = perMap.find(eid);
                    bool realign = false;
                    if (perMap.count(eid) > 0)
                    {
                        PeriodicEntity ent = it->second[0];
                        if (ent.isLocal == false)
                        {
                            realign = eid == min(eid, ent.id);
                            eid = min(eid, ent.id);
                        }
                    }

                    for (int j = 0; j < quad; ++j)
                    {
                        m_traceToUniversalMap[j+offset] = eid*maxQuad+j+1;
                        m_edgeToTrace[eid].emplace_back(offset + j);
                    }

                    if (realign)
                    {
                        if (nDim == 2)
                        {
                            RealignTraceElement(
                                tmp = m_traceToUniversalMap+offset,
                                it->second[0].orient, quad);
                        }
                        else
                        {
                            RealignTraceElement(
                                tmp = m_traceToUniversalMap+offset,
                                it->second[0].orient,
                                trace->GetExp(i)->GetNumPoints(0),
                                trace->GetExp(i)->GetNumPoints(1));
                        }
                    }
                }
            }

            //Set up graph topology for MPI
            //This creates a list of all geometry of problem dimension - 1
            std::vector<int> localEdgeIds; //
            for (int i = 0; i < locExpVector.size(); ++i)
            {
                eid = i;
                locExpansion = locExpVector[eid];
                nDim = locExpansion->GetShapeDimension();
                if (nDim == 1)
                {
                    int nVerts = locExpansion->GetNverts();
                    for (int j = 0; j < nVerts; ++j)
                    {
                        LocalRegions::PointExpSharedPtr locPointExp =
                                m_elmtToTrace[eid][j]->as<LocalRegions::PointExp>();
                        int id = locPointExp->GetGeom()->GetGlobalID();
                        localEdgeIds.emplace_back(id);
                    }
                }
                else if (nDim == 2)
                {
                    for (int j = 0; j < locExpansion->GetNedges(); ++j)
                    {
                        LocalRegions::SegExpSharedPtr locSegExp =
                                m_elmtToTrace[eid][j]->as<LocalRegions::SegExp>();
                        int id  = locSegExp->GetGeom()->GetGlobalID();
                        localEdgeIds.emplace_back(id);
                    }
                }
                else if (nDim == 3)
                {
                    for(int j = 0; j < locExpansion->GetNfaces(); ++j)
                    {
                        LocalRegions::Expansion2DSharedPtr locFaceExp =
                                m_elmtToTrace[eid][j]
                                        ->as<LocalRegions::Expansion2D>();
                        int id  = locFaceExp->GetGeom()->GetGlobalID();
                        localEdgeIds.emplace_back(id);
                    }
                }
            }

            //Create list of boundary edge IDs
            std::set<int> bndIdList;
            for (int i = 0; i < bndCond.num_elements(); ++i)
            {
                for (int j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    bndIdList.insert(bndCondExp[i]->GetExp(j)->GetGeom()->GetGlobalID());
                }
            }

            //Get unique edges to send (could be a better way of doing this? or is it even worth doing?)
            std::vector<int> uniqueEdgeIds;
            std::vector<bool> duplicated(localEdgeIds.size(), false);
            for (int i = 0; i < localEdgeIds.size(); ++i)
            {
                for (int j = i + 1; j < localEdgeIds.size(); ++j)
                {
                    if (localEdgeIds[i] == localEdgeIds[j])
                    {
                        duplicated[i] = duplicated[j] = true;
                    }
                }

                if (!duplicated[i])
                {
                    if (bndIdList.find(localEdgeIds[i]) == bndIdList.end())
                    {
                        uniqueEdgeIds.emplace_back(localEdgeIds[i]);
                    }
                }
            }

            // Send uniqueEdgeIds size so all partitions can prepare buffers
            int nRanks =  m_comm->GetSize();
            Array<OneD, int> rankNumEdges(nRanks);
            Array<OneD, int> localEdgeSize(1, uniqueEdgeIds.size());
            m_comm->AllGather(localEdgeSize, rankNumEdges);

            Array<OneD, int> rankLocalEdgeDisp(nRanks, 0);
            for (i = 1; i < nRanks; ++i)
            {
                rankLocalEdgeDisp[i] = rankLocalEdgeDisp[i-1] + rankNumEdges[i-1];
            }

            Array<OneD, int> localEdgeIdsArray(uniqueEdgeIds.size());
            for (i = 0; i < uniqueEdgeIds.size(); ++i)
            {
                localEdgeIdsArray[i] = uniqueEdgeIds[i];
            }

            //Sort localEdgeIdsArray before sending (this is important!)
            std::sort(localEdgeIdsArray.begin(), localEdgeIdsArray.end());

            Array<OneD, int> rankLocalEdgeIds(std::accumulate(
                    rankNumEdges.begin(), rankNumEdges.end(), 0), 0);

            //Send all unique edge IDs to all partitions
            m_comm->AllGatherv(localEdgeIdsArray, rankLocalEdgeIds, rankNumEdges, rankLocalEdgeDisp);

            //Create an array of other rank IDs (all except mine)
            int myRank = m_comm->GetRank();
            int cnt = 0;
            Array<OneD, int> otherRanks(nRanks - 1);
            for(i = 0; i < nRanks; ++i)
            {
                if (i != myRank)
                {
                    otherRanks[cnt] = i;
                    ++cnt;
                }
            }

            //Find what edge Ids match with other ranks
            std::map<int, std::vector<int>> rankSharedEdges;
            for (auto &rank : otherRanks)
            {
                for (int j = 0; j < rankNumEdges[rank]; ++j)
                {
                    int edgeId = rankLocalEdgeIds[rankLocalEdgeDisp[rank] + j];
                    int *found = std::find(localEdgeIdsArray.begin(), localEdgeIdsArray.end(), edgeId);
                    if (found != localEdgeIdsArray.end())
                    {
                        rankSharedEdges[rank].emplace_back(edgeId);
                    }
                }
            }
            // Setting up MPI neighbor_alltoallv
            //Prepare destinations & weights for creation of new distributed
            // MPI communicator for the neighbour communication
            int nNeighbours = rankSharedEdges.size();
            Array<OneD,int> destinations(nNeighbours, 0);
            Array<OneD,int> weights(nNeighbours, 0);
            cnt = 0;
            for (auto &rankEdgeVec : rankSharedEdges)
            {
                destinations[cnt] = rankEdgeVec.first;
                weights[cnt] = rankEdgeVec.second.size();
                ++cnt;
            }

            int retval = MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD,
                                                        nNeighbours, destinations.get(), weights.get(),  // Sources
                                                        nNeighbours, destinations.get(), weights.get(),  // Destinations
                                                        MPI_INFO_NULL, 1, &m_commGraph);

            ASSERTL0(retval == MPI_SUCCESS, "MPI error creating the distributed graph.");

            //Setting up indices
            m_sendCount = Array<OneD, int>(nNeighbours, 0);
            cnt = 0;
            for (auto &rankEdgeSet : rankSharedEdges)
            {
                for (int i : rankEdgeSet.second)
                {
                    std::vector<int> edgeIndex = m_edgeToTrace[i];
                    m_edgeTraceIndex.insert(m_edgeTraceIndex.end(), edgeIndex.begin(), edgeIndex.end());
                    m_sendCount[cnt] += edgeIndex.size();
                }

                ++cnt;
            }

            m_sendDisp = Nektar::Array<OneD, int>(nNeighbours, 0);
            for (i = 1; i < nNeighbours; ++i)
            {
                m_sendDisp[i] = m_sendDisp[i-1] + m_sendCount[i-1];
            }

            // Setting up standard MPI alltoall
            // Get maxCount which is the largest shared partition edge
            int maxCount = 0;
            for (int i = 0; i < nRanks; ++i)
            {
                if (rankSharedEdges.find(i) != rankSharedEdges.end())
                {
                    maxCount = (rankSharedEdges[i].size() > maxCount) ? rankSharedEdges[i].size() : maxCount;
                }
            }

            m_comm->AllReduce(maxCount, LibUtilities::ReduceMax);

            // Creates the edge index vector where value -1 indicates
            // padding of value 0 to be inserted instead of value from Fwd
            std::vector<int> allEdgeIndex;
            for (int i = 0; i < nRanks; ++i)
            {
                if (rankSharedEdges.find(i) != rankSharedEdges.end())
                {
                    for (int j = 0; j < rankSharedEdges[i].size(); ++j)
                    {
                        std::vector<int> edgeIndex = m_edgeToTrace[rankSharedEdges[i][j]];
                        if (edgeIndex.size() < maxQuad)
                        {
                            std::vector<int> diff(maxQuad - edgeIndex.size(), -1);
                            edgeIndex.insert(edgeIndex.end(), diff.begin(), diff.end());
                        }

                        allEdgeIndex.insert(allEdgeIndex.end(), edgeIndex.begin(), edgeIndex.end());
                    }

                    if (rankSharedEdges[i].size() < maxCount)
                    {
                        std::vector<int> edgeIndex(maxQuad * (maxCount - rankSharedEdges[i].size()), -1);
                        allEdgeIndex.insert(allEdgeIndex.end(), edgeIndex.begin(), edgeIndex.end());
                    }
                }
                else
                {
                    std::vector<int> edgeIndex(maxQuad * maxCount, -1);
                    allEdgeIndex.insert(allEdgeIndex.end(), edgeIndex.begin(), edgeIndex.end());
                }
            }

            // Setting up standard MPI alltoallv
            std::vector<int> allVEdgeIndex;
            Array<OneD, int> allVSendCount(nRanks, 0);
            for (int i = 0; i < nRanks; ++i)
            {
                if (rankSharedEdges.find(i) != rankSharedEdges.end())
                {
                    for (int j = 0; j < rankSharedEdges[i].size(); ++j)
                    {
                        std::vector<int> edgeIndex = m_edgeToTrace[rankSharedEdges[i][j]];
                        allVEdgeIndex.insert(allVEdgeIndex.end(), edgeIndex.begin(), edgeIndex.end());
                        allVSendCount[i] += edgeIndex.size();
                    }
                }
                else
                {
                    allVSendCount[i] = 0;
                }
            }

            Array<OneD, int> allVSendDisp(nRanks, 0);
            for (i = 1; i < nRanks; ++i)
            {
                allVSendDisp[i] = allVSendDisp[i-1] + allVSendCount[i-1];
            }

            //Setting up pairwise send/recv between individual partitions
            std::vector<std::pair<int, std::vector<int>>> vecPairPartitionTrace;
            int totSends = 0;
            for (const auto& rankEdgeSet : rankSharedEdges)
            {
                std::vector<int> edgeTraceIndex;
                for (int i : rankEdgeSet.second)
                {
                    std::vector<int> edgeIndex = m_edgeToTrace[i];
                    edgeTraceIndex.insert(edgeTraceIndex.end(), edgeIndex.begin(), edgeIndex.end());

                    totSends += m_edgeToTrace[i].size();
                }

                vecPairPartitionTrace.emplace_back(std::make_pair(rankEdgeSet.first, edgeTraceIndex));
            }

            //Timing MPI comm methods
            int testLoopCount = 100;
            Array<OneD, double> testFwd(trace->GetNpoints(), 1);
            Array<OneD, double> testBwd(trace->GetNpoints(), -2);

            NekDouble timeAvgAll = -1;
            NekDouble timeAvgAllv = -1;
            NekDouble timeAvgNeighbourV = -1;
            NekDouble timeAvgPairwise = -1;

            //Alltoall
            LibUtilities::Timer tAll;
            tAll.Start();

            for (int i = 0; i < testLoopCount; ++i)
            {
                int size = maxQuad * maxCount * nRanks;
                Array<OneD, double> sendBuff(size, -1);
                Array<OneD, double> recvBuff(size, -1);

                for (int j = 0; j < size; ++j)
                {
                    if (allEdgeIndex[j] == -1)
                    {
                        sendBuff[j] = 0;
                    }
                    else
                    {
                        sendBuff[j] = testFwd[allEdgeIndex[j]];
                    }
                }

                m_comm->AlltoAll(sendBuff, recvBuff);

                for (int j = 0; j < size; ++j)
                {
                    if (allEdgeIndex[j] != -1)
                    {
                        testBwd[allEdgeIndex[j]] = recvBuff[j];
                    }
                }
            }

            tAll.Stop();
            Array<OneD, NekDouble> minTimeAll(1, tAll.TimePerTest(testLoopCount));
            m_comm->AllReduce(minTimeAll, LibUtilities::ReduceMin);

            Array<OneD, NekDouble> maxTimeAll(1, tAll.TimePerTest(testLoopCount));
            m_comm->AllReduce(maxTimeAll, LibUtilities::ReduceMax);

            Array<OneD, NekDouble> sumTimeAll(1, tAll.TimePerTest(testLoopCount));
            m_comm->AllReduce(sumTimeAll, LibUtilities::ReduceSum);

            if (myRank == 0)
            {
                timeAvgAll = sumTimeAll[0]/nRanks;
                std::cout << "alltoall times (avg, min, max):   " << timeAvgAll << " " << minTimeAll[0] << " " << maxTimeAll[0] << std::endl;
            }

            //Alltoallv
            LibUtilities::Timer tAllv;
            tAllv.Start();

            for (int i = 0; i < testLoopCount; ++i)
            {
                Array<OneD, double> sendBuff(allVEdgeIndex.size(), -1);
                Array<OneD, double> recvBuff(allVEdgeIndex.size(), -1);

                for (int i = 0; i < allVEdgeIndex.size(); ++i)
                {
                    sendBuff[i] = testFwd[allVEdgeIndex[i]];
                }

                m_comm->AlltoAllv(sendBuff, allVSendCount, allVSendDisp,
                                  recvBuff, allVSendCount, allVSendDisp);

                for (int i = 0; i < allVEdgeIndex.size(); ++i)
                {
                    testBwd[allVEdgeIndex[i]] = recvBuff[i];
                }
            }

            tAllv.Stop();

            Array<OneD, NekDouble> minTimeAllv(1, tAllv.TimePerTest(testLoopCount));
            m_comm->AllReduce(minTimeAllv, LibUtilities::ReduceMin);

            Array<OneD, NekDouble> maxTimeAllv(1, tAllv.TimePerTest(testLoopCount));
            m_comm->AllReduce(maxTimeAllv, LibUtilities::ReduceMax);

            Array<OneD, NekDouble> sumTimeAllv(1, tAllv.TimePerTest(testLoopCount));
            m_comm->AllReduce(sumTimeAllv, LibUtilities::ReduceSum);

            if (myRank == 0)
            {
                timeAvgAllv = sumTimeAllv[0]/nRanks;
                std::cout << "alltoallv times (avg, min, max):  " << timeAvgAllv << " " << minTimeAllv[0] << " " << maxTimeAllv[0] << std::endl;
            }

            //neighbor_alltoallv
            LibUtilities::Timer tNeighbourAllv;
            tNeighbourAllv.Start();

            for (int i = 0; i < testLoopCount; ++i)
            {
                Array<OneD, double> sendBuff(m_edgeTraceIndex.size(), -1);
                Array<OneD, double> recvBuff(m_edgeTraceIndex.size(), -1);
                for (int i = 0; i < m_edgeTraceIndex.size(); ++i)
                {
                    sendBuff[i] = testFwd[m_edgeTraceIndex[i]];
                }


                MPI_Neighbor_alltoallv(sendBuff.get(), m_sendCount.get(), m_sendDisp.get(), MPI_DOUBLE,
                                       recvBuff.get(), m_sendCount.get(), m_sendDisp.get(), MPI_DOUBLE,
                                       m_commGraph);

                for (int i = 0; i < m_edgeTraceIndex.size(); ++i)
                {
                    testBwd[m_edgeTraceIndex[i]] = recvBuff[i];
                }
            }

            tNeighbourAllv.Stop();
            NekDouble timeNeighbourAllvLocal = tNeighbourAllv.TimePerTest(testLoopCount);

            Array<OneD, NekDouble> timeNeighbourAllv(nRanks, 0.0);
            MPI_Gather(&timeNeighbourAllvLocal, 1, MPI_DOUBLE, timeNeighbourAllv.get(), 1, MPI_DOUBLE, 0, m_commGraph);

            if (myRank == 0)
            {
                NekDouble totalTime = timeNeighbourAllv[0];
                NekDouble minTime = timeNeighbourAllv[0];
                NekDouble maxTime = timeNeighbourAllv[0];

                for (int i = 1; i < nRanks; ++i)
                {
                    totalTime += timeNeighbourAllv[i];
                    minTime = (timeNeighbourAllv[i] < minTime) ? timeNeighbourAllv[i] : minTime;
                    maxTime = (timeNeighbourAllv[i] > maxTime) ? timeNeighbourAllv[i] : maxTime;
                }

                timeAvgNeighbourV = totalTime/nRanks;

                std::cout << "n_alltoall times (avg, min, max): " << timeAvgNeighbourV << " " << minTime << " " << maxTime << std::endl;
            }

            //pairwise send/recv
            LibUtilities::Timer tPairwise;
            tPairwise.Start();

            for (int i = 0; i < testLoopCount; ++i)
            {
                Array<OneD, MPI_Request> request(vecPairPartitionTrace.size() * 2);
                Array<OneD, MPI_Status> status(vecPairPartitionTrace.size() * 2);
                int count = 0, count2 = 0;

                Array<OneD, NekDouble> recvBuff(totSends, -1);
                for (auto pairPartitionTrace : vecPairPartitionTrace)
                {
                    Array<OneD, NekDouble> sendBuff(pairPartitionTrace.second.size(), -1);
                    for (int i = 0; i < pairPartitionTrace.second.size(); ++i)
                    {
                        sendBuff[i] = testFwd[pairPartitionTrace.second[i]];
                    }

                    MPI_Isend(sendBuff.get(),
                            sendBuff.num_elements(),
                            MPI_DOUBLE,
                            pairPartitionTrace.first,  // rank of destination
                            0,                         // message tag
                            m_commGraph,
                            &request[count2++]);

                    MPI_Irecv(static_cast<void *>(recvBuff.get() + m_sendDisp[count++]),
                            sendBuff.num_elements(),
                            MPI_DOUBLE,
                            pairPartitionTrace.first,  // rank of source
                            0,                         // message tag
                            m_commGraph,
                            &request[count2++]);

                }

                MPI_Waitall(vecPairPartitionTrace.size() * 2, request.get(), status.get());

                count = 0;
                for (auto pairPartitionTrace : vecPairPartitionTrace)
                {
                    for (int i = 0; i < pairPartitionTrace.second.size(); ++i)
                    {
                        testBwd[pairPartitionTrace.second[i]] = recvBuff[m_sendDisp[count] + i];
                    }
                    count++;
                }
            }

            tPairwise.Stop();
            Array<OneD, NekDouble> minTimePairwise(1, tPairwise.TimePerTest(testLoopCount));
            m_comm->AllReduce(minTimePairwise, LibUtilities::ReduceMin);

            Array<OneD, NekDouble> maxTimePairwise(1, tPairwise.TimePerTest(testLoopCount));
            m_comm->AllReduce(maxTimePairwise, LibUtilities::ReduceMax);

            Array<OneD, NekDouble> sumTimePairwise(1, tPairwise.TimePerTest(testLoopCount));
            m_comm->AllReduce(sumTimePairwise, LibUtilities::ReduceSum);

            if (myRank == 0)
            {
                timeAvgPairwise = sumTimePairwise[0]/nRanks;
                std::cout << "pairwise times (avg, min, max):   " << timeAvgPairwise << " " << minTimePairwise[0] << " " << maxTimePairwise[0] << std::endl;
            }

            Array<OneD, long> tmp2(nTracePhys);
            for (int i = 0; i < nTracePhys; ++i)
            {
                tmp2[i] = m_traceToUniversalMap[i];
            }
            m_traceGsh = Gs::Init(tmp2, m_comm);
            Gs::Unique(tmp2, m_comm);
            for (int i = 0; i < nTracePhys; ++i)
            {
                m_traceToUniversalMapUnique[i] = tmp2[i];
            }
        }

        void AssemblyMapDG::RealignTraceElement(
            Array<OneD, int>        &toAlign,
            StdRegions::Orientation  orient,
            int                      nquad1,
            int                      nquad2)
        {
            if (orient == StdRegions::eBackwards)
            {
                ASSERTL1(nquad2 == 0, "nquad2 != 0 for reorienation");
                for (int i = 0; i < nquad1/2; ++i)
                {
                    swap(toAlign[i], toAlign[nquad1-1-i]);
                }
            }
            else if (orient != StdRegions::eForwards)
            {
                ASSERTL1(nquad2 != 0, "nquad2 == 0 for reorienation");

                Array<OneD, int> tmp(nquad1*nquad2);

                // Copy transpose.
                if (orient == StdRegions::eDir1FwdDir2_Dir2FwdDir1 ||
                    orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
                    orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1 ||
                    orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
                {
                    for (int i = 0; i < nquad2; ++i)
                    {
                        for (int j = 0; j < nquad1; ++j)
                        {
                            tmp[i*nquad1 + j] = toAlign[j*nquad2 + i];
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < nquad2; ++i)
                    {
                        for (int j = 0; j < nquad1; ++j)
                        {
                            tmp[i*nquad1 + j] = toAlign[i*nquad1 + j];
                        }
                    }
                }

                if (orient == StdRegions::eDir1BwdDir1_Dir2FwdDir2 ||
                    orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
                    orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
                    orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
                {
                    // Reverse x direction
                    for (int i = 0; i < nquad2; ++i)
                    {
                        for (int j = 0; j < nquad1/2; ++j)
                        {
                            swap(tmp[i*nquad1 + j],
                                 tmp[i*nquad1 + nquad1-j-1]);
                        }
                    }
                }

                if (orient == StdRegions::eDir1FwdDir1_Dir2BwdDir2 ||
                    orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
                    orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1 ||
                    orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
                {
                    // Reverse y direction
                    for (int j = 0; j < nquad1; ++j)
                    {
                        for (int i = 0; i < nquad2/2; ++i)
                        {
                            swap(tmp[i*nquad1 + j],
                                 tmp[(nquad2-i-1)*nquad1 + j]);
                        }
                    }
                }
                Vmath::Vcopy(nquad1*nquad2, tmp, 1, toAlign, 1);
            }
        }

        void AssemblyMapDG::UniversalTraceAssemble(Array<OneD, NekDouble> &Fwd, Array<OneD, NekDouble> &Bwd)
        {
            Array<OneD, double> sendBuff(m_edgeTraceIndex.size(), -1);
            Array<OneD, double> recvBuff(m_edgeTraceIndex.size(), -1);
            for (int i = 0; i < m_edgeTraceIndex.size(); ++i)
            {
                sendBuff[i] = Fwd[m_edgeTraceIndex[i]];
            }

            MPI_Neighbor_alltoallv(sendBuff.get(), m_sendCount.get(), m_sendDisp.get(), MPI_DOUBLE,
                                   recvBuff.get(), m_sendCount.get(), m_sendDisp.get(), MPI_DOUBLE,
                                   m_commGraph);

            for (int i = 0; i < m_edgeTraceIndex.size(); ++i)
            {
                Bwd[m_edgeTraceIndex[i]] = recvBuff[i];
            }
        }

        void AssemblyMapDG::UniversalTraceAssembleGS(Array<OneD, NekDouble> &pGlobal) const
        {
            Gs::Gather(pGlobal, Gs::gs_add, m_traceGsh);

        }
        int AssemblyMapDG::v_GetLocalToGlobalMap(const int i) const
        {
            return m_localToGlobalBndMap[i];
        }

        int AssemblyMapDG::v_GetGlobalToUniversalMap(const int i) const
        {
            return m_globalToUniversalBndMap[i];
        }

        int AssemblyMapDG::v_GetGlobalToUniversalMapUnique(const int i) const
        {
            return m_globalToUniversalBndMapUnique[i];
        }

        const Array<OneD,const int>& AssemblyMapDG::v_GetLocalToGlobalMap()
        {
            return m_localToGlobalBndMap;
        }

        const Array<OneD,const int>& AssemblyMapDG::v_GetGlobalToUniversalMap()
        {
            return m_globalToUniversalBndMap;
        }

        const Array<OneD,const int>& AssemblyMapDG::v_GetGlobalToUniversalMapUnique()
        {
            return m_globalToUniversalBndMapUnique;
        }

        NekDouble AssemblyMapDG::v_GetLocalToGlobalSign(
                    const int i) const
        {
            return GetLocalToGlobalBndSign(i);
        }

        void AssemblyMapDG::v_LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                    Array<OneD,       NekDouble>& global,
                    bool useComm ) const
        {
            boost::ignore_unused(useComm);
            AssembleBnd(loc,global);
        }

        void AssemblyMapDG::v_LocalToGlobal(
                    const NekVector<NekDouble>& loc,
                    NekVector<      NekDouble>& global,
                    bool useComm) const
        {
            boost::ignore_unused(useComm);
            AssembleBnd(loc,global);
        }

        void AssemblyMapDG::v_GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const
        {
            GlobalToLocalBnd(global,loc);
        }

        void AssemblyMapDG::v_GlobalToLocal(
                    const NekVector<NekDouble>& global,
                          NekVector<      NekDouble>& loc) const
        {
            GlobalToLocalBnd(global,loc);
        }

        void AssemblyMapDG::v_Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const
        {
            AssembleBnd(loc,global);
        }

        void AssemblyMapDG::v_Assemble(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const
        {
            AssembleBnd(loc,global);
        }

        void AssemblyMapDG::v_UniversalAssemble(
                      Array<OneD,     NekDouble>& pGlobal) const
        {
            Gs::Gather(pGlobal, Gs::gs_add, m_gsh);
        }

        void AssemblyMapDG::v_UniversalAssemble(
                      NekVector<      NekDouble>& pGlobal) const
        {
            UniversalAssemble(pGlobal.GetPtr());
        }

        int AssemblyMapDG::v_GetFullSystemBandWidth() const
        {
            return GetBndSystemBandWidth();
        }

        int AssemblyMapDG::GetTraceToUniversalMap(int i)
        {
            return m_traceToUniversalMap[i];
        }

        int AssemblyMapDG::GetTraceToUniversalMapUnique(int i)
        {
            return m_traceToUniversalMapUnique[i];
        }

        int AssemblyMapDG::GetNumDirichletBndPhys()
        {
            return m_numDirichletBndPhys;
        }

        Array<OneD, LocalRegions::ExpansionSharedPtr>&
                    AssemblyMapDG::GetElmtToTrace(const int i)
        {
            ASSERTL1(i >= 0 && i < m_elmtToTrace.num_elements(),
                     "i is out of range");
            return m_elmtToTrace[i];
        }

        Array<OneD, Array< OneD, LocalRegions::ExpansionSharedPtr> >&
                    AssemblyMapDG::GetElmtToTrace()
        {
            return m_elmtToTrace;
        }
    } //namespace
} // namespace
