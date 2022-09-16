///////////////////////////////////////////////////////////////////////////////
//
// File PulseWaveSystem.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Generic timestepping for Pulse Wave Solver
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Communication/CommSerial.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <MultiRegions/ContField.h>
#include <PulseWaveSolver/EquationSystems/PulseWaveSystem.h>
using namespace std;

namespace Nektar
{

/**
 *  @class PulseWaveSystem
 *
 *  Initialises the arterial subdomains in m_vessels and sets up
 *  all domain-linking conditions (bifurcations, intefaces/junctions,
 *  merging flows). Detects the network structure and assigns
 *  boundary conditons. Also provides the underlying timestepping
 *  framework for pulse wave solvers including the general
 *  timestepping routines.
 */

/**
 *  Processes SolverInfo parameters from the session file and sets
 *  up timestepping-specific code.
 *
 *  @param   m_Session        Session object to read parameters from.
 */
PulseWaveSystem::PulseWaveSystem(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph)
{
}

/**
 *  Destructor
 */
PulseWaveSystem::~PulseWaveSystem()
{
}

/**
 *  Initialisation routine for multidomain solver. Sets up the
 *  expansions for every arterial segment (m_vessels) and for one
 *  complete field m_outfield which is needed to write the
 *  postprocessing output. Also determines which upwind strategy
 *  is used (currently only upwinding scheme available) and reads
 *  blodd flow specific parameters from the inputfile
 *
 */
void PulseWaveSystem::v_InitObject(bool DeclareField)
{
    // Initialise base class
    UnsteadySystem::v_InitObject(DeclareField);

    // Read the geometry and the expansion information
    m_domain   = m_graph->GetDomain();
    m_nDomains = m_domain.size();

    // Determine projectiontype
    ASSERTL0(m_session->MatchSolverInfo("Projection", "DisContinuous"),
             "Pulse solver only set up for Discontinuous projections");
    m_projectionType = MultiRegions::eDiscontinuous;
    ASSERTL0(
        m_graph->GetMeshDimension() == 1,
        "Pulse wave solver only set up for expansion dimension equal to 1");

    int i;
    m_nVariables = m_session->GetVariables().size();

    m_fields = Array<OneD, MultiRegions::ExpListSharedPtr>(m_nVariables);
    m_vessels =
        Array<OneD, MultiRegions::ExpListSharedPtr>(m_nVariables * m_nDomains);

    /* If the PressureArea property is not specified, default to the Beta
     * law; it's the most well known and this way previous examples that did
     * not specify the tube law will still work.
     */
    if (m_session->DefinesSolverInfo("PressureArea"))
    {
        m_pressureArea = GetPressureAreaFactory().CreateInstance(
            m_session->GetSolverInfo("PressureArea"), m_vessels, m_session);
    }
    else
    {
        m_pressureArea = GetPressureAreaFactory().CreateInstance(
            "Beta", m_vessels, m_session);
    }

    m_pressure = Array<OneD, Array<OneD, NekDouble>>(m_nDomains);

    m_PWV = Array<OneD, Array<OneD, NekDouble>>(m_nDomains);
    m_W1  = Array<OneD, Array<OneD, NekDouble>>(m_nDomains);
    m_W2  = Array<OneD, Array<OneD, NekDouble>>(m_nDomains);

    SpatialDomains::BoundaryConditions Allbcs(m_session, m_graph);

    // Set up domains and put geometry to be only one space dimension.
    int cnt                     = 0;
    bool SetToOneSpaceDimension = true;

    if (m_session->DefinesCmdLineArgument("SetToOneSpaceDimension"))
    {
        std::string cmdline = m_session->GetCmdLineArgument<std::string>(
            "SetToOneSpaceDimension");
        if (boost::to_upper_copy(cmdline) == "FALSE")
        {
            SetToOneSpaceDimension = false;
        }
    }

    // get parallel communicators as required
    map<int, LibUtilities::CommSharedPtr> domComm;
    GetCommArray(domComm);

    SetUpDomainInterfaceBCs(Allbcs, domComm);

    for (auto &d : m_domOrder)
    {
        for (int j = 0; j < m_nVariables; ++j)
        {
            m_vessels[cnt++] =
                MemoryManager<MultiRegions::DisContField>::AllocateSharedPtr(
                    m_session, m_graph, m_domain[d], Allbcs,
                    m_session->GetVariable(j), domComm[d],
                    SetToOneSpaceDimension);
        }
    }

    // Reset coeff and phys space to be continuous over all domains
    int totcoeffs     = 0;
    int totphys       = 0;
    m_fieldPhysOffset = Array<OneD, int>(m_nDomains + 1, 0);
    for (i = 0; i < m_nDomains; ++i)
    {
        totcoeffs += m_vessels[i * m_nVariables]->GetNcoeffs();

        m_fieldPhysOffset[i] = totphys;
        totphys += m_vessels[i * m_nVariables]->GetTotPoints();
    }
    m_fieldPhysOffset[m_nDomains] = totphys;

    for (int n = 0; n < m_nVariables; ++n)
    {
        Array<OneD, NekDouble> coeffs(totcoeffs, 0.0);
        Array<OneD, NekDouble> phys(totphys, 0.0);
        Array<OneD, NekDouble> tmpcoeffs, tmpphys;

        m_vessels[n]->SetCoeffsArray(coeffs);
        m_vessels[n]->SetPhysArray(phys);

        int cnt  = m_vessels[n]->GetNcoeffs();
        int cnt1 = m_vessels[n]->GetTotPoints();

        for (i = 1; i < m_nDomains; ++i)
        {
            m_vessels[i * m_nVariables + n]->SetCoeffsArray(tmpcoeffs =
                                                                coeffs + cnt);
            m_vessels[i * m_nVariables + n]->SetPhysArray(tmpphys =
                                                              phys + cnt1);
            cnt += m_vessels[i * m_nVariables + n]->GetNcoeffs();
            cnt1 += m_vessels[i * m_nVariables + n]->GetTotPoints();
        }
    }

    // If Discontinuous Galerkin determine upwinding method to use
    for (int i = 0; i < (int)SIZE_UpwindTypePulse; ++i)
    {
        bool match;
        m_session->MatchSolverInfo("UPWINDTYPEPULSE", UpwindTypeMapPulse[i],
                                   match, false);
        if (match)
        {
            m_upwindTypePulse = (UpwindTypePulse)i;
            break;
        }
    }

    // Load blood density and external pressure
    m_session->LoadParameter("rho", m_rho, 0.5);
    m_session->LoadParameter("pext", m_pext, 0.0);

    int nq = 0;
    /*
     *  Gets the Material Properties of each arterial segment
     *  specified in the inputfile from section MaterialProperties
     *  Also gets the Area at static equilibrium A_0 specified in the
     *  inputfile.
     *
     * Having found these points also extract the values at the
     * trace points and the normal direction consistent with the
     * left adjacent definition of Fwd and Bwd
     */
    m_beta             = Array<OneD, Array<OneD, NekDouble>>(m_nDomains);
    m_beta_trace       = Array<OneD, Array<OneD, NekDouble>>(m_nDomains);
    m_gamma            = Array<OneD, Array<OneD, NekDouble>>(m_nDomains);
    m_alpha            = Array<OneD, Array<OneD, NekDouble>>(m_nDomains);
    m_alpha_trace      = Array<OneD, Array<OneD, NekDouble>>(m_nDomains);
    m_A_0              = Array<OneD, Array<OneD, NekDouble>>(m_nDomains);
    m_A_0_trace        = Array<OneD, Array<OneD, NekDouble>>(m_nDomains);
    m_trace_fwd_normal = Array<OneD, Array<OneD, NekDouble>>(m_nDomains);

    int omega = 0;
    for (auto &d : m_domOrder)
    {
        nq          = m_vessels[2 * omega]->GetNpoints();
        m_fields[0] = m_vessels[2 * omega];

        m_beta[omega] = Array<OneD, NekDouble>(nq);
        GetFunction("MaterialProperties")
            ->Evaluate("beta", m_beta[omega], m_time, d);

        // If the input file doesn't specify viscoelasticity, set it to zero.
        m_gamma[omega] = Array<OneD, NekDouble>(nq);

        if (m_session->DefinesFunction("Viscoelasticity"))
        {
            GetFunction("Viscoelasticity")
                ->Evaluate("gamma", m_gamma[omega], m_time, d);
        }
        else
        {
            for (int j = 0; j < nq; ++j)
            {
                m_gamma[omega][j] = 0;
            }
        }

        /* If the input file doesn't specify strain-stiffening, set it to
         * 0.5 for elastic behaviour.
         */
        m_alpha[omega] = Array<OneD, NekDouble>(nq);

        if (m_session->DefinesFunction("StrainStiffening"))
        {
            GetFunction("StrainStiffening")
                ->Evaluate("alpha", m_alpha[omega], m_time, d);
        }
        else
        {
            for (int j = 0; j < nq; ++j)
            {
                m_alpha[omega][j] = 0.5;
            }
        }

        m_A_0[omega] = Array<OneD, NekDouble>(nq);
        GetFunction("A_0")->Evaluate("A_0", m_A_0[omega], m_time, d);

        int nqTrace = GetTraceTotPoints();

        m_beta_trace[omega] = Array<OneD, NekDouble>(nqTrace);
        m_fields[0]->ExtractTracePhys(m_beta[omega], m_beta_trace[omega]);

        m_A_0_trace[omega] = Array<OneD, NekDouble>(nqTrace);
        m_fields[0]->ExtractTracePhys(m_A_0[omega], m_A_0_trace[omega]);

        m_alpha_trace[omega] = Array<OneD, NekDouble>(nqTrace);
        m_fields[0]->ExtractTracePhys(m_alpha[omega], m_alpha_trace[omega]);

        if (SetToOneSpaceDimension)
        {
            m_trace_fwd_normal[omega] = Array<OneD, NekDouble>(nqTrace, 0.0);

            MultiRegions::ExpListSharedPtr trace = m_fields[0]->GetTrace();
            int nelmt_trace                      = trace->GetExpSize();

            Array<OneD, Array<OneD, NekDouble>> normals(nelmt_trace);

            for (int i = 0; i < nelmt_trace; ++i)
            {
                normals[i] = m_trace_fwd_normal[omega] + i;
            }

            // need to set to 1 for consistency since boundary
            // conditions may not have coordim = 1
            trace->GetExp(0)->GetGeom()->SetCoordim(1);

            trace->GetNormals(normals);
        }
        omega++;
    }

    SetUpDomainInterfaces();
}

void PulseWaveSystem::GetCommArray(
    map<int, LibUtilities::CommSharedPtr> &retval)
{

    // set up defualt serial communicator
    std::string def("default");
    char *argv = new char[def.length() + 1];
    std::strcpy(argv, def.c_str());
    LibUtilities::CommSharedPtr serialComm =
        MemoryManager<LibUtilities::CommSerial>::AllocateSharedPtr(1, &argv);

    int nprocs = m_comm->GetSize();

    if (nprocs == 1) // serial case
    {
        for (auto &d : m_domain)
        {
            // serial communicator
            retval[d.first] = serialComm;
            // fill out m_domOrder set for serial case
            m_domOrder.push_back(d.first);
        }
    }
    else // parallel case
    {
        int rank = m_comm->GetRank();

        // Fill array with domain details
        int dmax = 0;
        for (auto &d : m_domain)
        {
            dmax = (d.first > dmax) ? d.first : dmax;
        }
        dmax += 1;

        // get domain maximum id + 1;
        m_comm->AllReduce(dmax, LibUtilities::ReduceMax);

        Array<OneD, int> commtmp(dmax, 0);
        for (auto &d : m_domain)
        {
            commtmp[d.first] = 1;
        }

        // identfy which domains are split between two procs
        m_comm->AllReduce(commtmp, LibUtilities::ReduceSum);

        Array<OneD, bool> DoneCreate(nprocs, false);

        // setup communicators for domain partitions greater than 2
        for (int d = 0; d < dmax; ++d)
        {
            // set up communicator for cases where more than 2 procs
            // share domain
            if (commtmp[d] > 2)
            {
                LibUtilities::CommSharedPtr newcomm;
                int flag = 0;
                // set flag to 10 if d is on this processors domain
                for (auto &d1 : m_domain)
                {
                    if (d1.first == d)
                    {
                        flag = 10;
                    }
                }

                newcomm = m_comm->CommCreateIf(flag);

                // set up communicator if domain on this processor
                for (auto &d1 : m_domain)
                {
                    if (d1.first == d)
                    {
                        retval[d] = newcomm;
                    }
                }

                // turn off commtmp flag since should now be setup
                commtmp[d] = -1;
            }
        }

        // we now should have a list of the number of each
        // domains/vessel that is split over two processors

        // set up serial communicators for domains on only one
        // processor and set up map of domains that require parallel
        // communicator and number of processors in parallel

        // communicator
        map<int, int> SharedProc;
        Array<OneD, int> commtmp1(dmax, 0);

        for (auto &d : m_domain)
        {
            if (commtmp[d.first] == 1) // set to serial communicator
            {
                retval[d.first] = serialComm;
            }
            else if (commtmp[d.first] == 2) // shared by 2 procs
            {
                SharedProc[d.first] = 1;
                commtmp1[d.first]   = rank;
            }
            else // do nothing
            {
            }
        }

        // find the maximum number of communicator on all processors
        int nShrProc = SharedProc.size();
        int commMax  = nShrProc;
        m_comm->AllReduce(commMax, LibUtilities::ReduceMax);

        // find out processor id of each communicator by adding rank
        // and then subtracting local rank

        m_comm->AllReduce(commtmp1, LibUtilities::ReduceSum);

        // reset SharedProc to hold shared proc id.
        for (auto &d : SharedProc)
        {
            SharedProc[d.first] = commtmp1[d.first] - rank;
        }

        // now work out a comm ordering;
        Array<OneD, int> procs(2);
        int commcall = 0;
        map<int, pair<int, int>> CreateComm;

        bool search = true;
        Array<OneD, int> setdone(3);

        while (search)
        {
            int dorank = 0; // search index over processors
            set<int>
                proclist;   // has proc been identified for a comm at this level
            int flag = 100; // colour for communicators in this search level

            while (dorank < nprocs)
            {
                int sharedproc = -1;
                int ncomms     = 0;

                if (DoneCreate[dorank] == false)
                {
                    // check to see if proc even has a shared domain
                    ncomms = (rank == dorank) ? nShrProc : 0;
                    m_comm->AllReduce(ncomms, LibUtilities::ReduceMax);

                    // set DoneCreate to true if n comms required on
                    // proc
                    DoneCreate[dorank] = (ncomms == 0) ? true : false;
                }

                // requires communicator and needs to create one.
                if (ncomms)
                {
                    bool doneCreateOnDoRank     = false;
                    bool doneCreateOnSharedProc = false;
                    bool createdComm            = false;

                    // check not already callng communicator on
                    // this processor for this search sweep
                    if ((proclist.count(dorank) == 0))
                    {

                        // give all processors rank id of shared proc
                        if (rank == dorank)
                        {
                            for (auto &d : SharedProc)
                            {
                                if (CreateComm.count(d.first) == 0)
                                {
                                    sharedproc = SharedProc[d.first];
                                    break;
                                }
                            }

                            ASSERTL1(sharedproc != -1,
                                     "Failed to fine a new "
                                     "processor to set up another "
                                     "communicator");

                            if (proclist.count(sharedproc) == 0)
                            {
                                // save all communicators on this for
                                // split comm setup
                                for (auto &d : SharedProc)
                                {
                                    if (d.second == sharedproc)
                                    {
                                        CreateComm[d.first] = pair<int, int>(
                                            commcall, flag + d.first);
                                        createdComm = true;
                                    }
                                }
                            }

                            // set bool to mark as done on dorank
                            if (CreateComm.size() == nShrProc)
                            {
                                doneCreateOnDoRank = true;
                            }
                        }
                        else
                        {
                            sharedproc = -1;
                        }

                        m_comm->AllReduce(sharedproc, LibUtilities::ReduceMax);
                        ASSERTL1(
                            sharedproc < nprocs,
                            "Error in obtaining rank "
                            "shared by domain. Sharedproc value is greater "
                            "than number of procs");

                        if (proclist.count(sharedproc) == 0)
                        {
                            if (rank == sharedproc)
                            {
                                // save all communicators on this for
                                // split comm setup
                                for (auto &d : SharedProc)
                                {
                                    if (d.second == dorank)
                                    {
                                        CreateComm[d.first] = pair<int, int>(
                                            commcall, flag + d.first);
                                        createdComm = true;
                                    }
                                }

                                if (CreateComm.size() == nShrProc)
                                {
                                    doneCreateOnSharedProc = true;
                                }
                            }
                        }
                    }

                    // mark DoneCreate as if complete
                    setdone[0] = (doneCreateOnDoRank) ? 1 : 0;
                    setdone[1] = (doneCreateOnSharedProc) ? 1 : 0;
                    setdone[2] = (createdComm) ? 1 : 0;
                    m_comm->AllReduce(setdone, LibUtilities::ReduceMax);
                    DoneCreate[dorank] = (bool)setdone[0];
                    if (sharedproc != -1)
                    {
                        DoneCreate[sharedproc] = (bool)setdone[1];
                    }
                    if (setdone[2]) // created new communicator
                    {
                        proclist.insert(dorank);
                        proclist.insert(sharedproc);
                    }
                }

                dorank++;
            }

            // have we found all comms on all processors.
            int foundallcomms = 0;
            int i;
            for (i = 0; i < nprocs; ++i)
            {
                if (DoneCreate[i] == false)
                {
                    break;
                }
            }

            if (i == nprocs)
            {
                search = false;
            }
            else
            {
                commcall++;
                ASSERTL0(commcall < 4 * commMax,
                         "Failed to find sub communicators "
                         "pattern in 4*commMax searches");
            }
        }

        ASSERTL1(CreateComm.size() == nShrProc,
                 "Have not created communicators for all shared procs");

        // determine maxmimum number of CreateComm size
        int maxCreateComm = CreateComm.size();
        m_comm->AllReduce(maxCreateComm, LibUtilities::ReduceMax);

        // loop over CreateComm list
        for (int i = 0; i < maxCreateComm; ++i)
        {
            LibUtilities::CommSharedPtr newcomm;
            int flag = 0;

            for (auto &d : CreateComm)
            {
                if (d.second.first == i)
                {
                    flag = d.second.second;
                }
            }

            newcomm = m_comm->CommCreateIf(flag);

            for (auto &d : CreateComm)
            {
                if (d.second.first == i)
                {
                    retval[d.first] = newcomm;
                }
            }
        }

        // Finally need to re-order domain list so that we call the
        // communicators in a non-locking manner. This is done by
        // uniquely numbering each shared proc/comm interface and
        // ensuring the domains are ordered so that we call this
        // numbering in an increasing order

        Array<OneD, NekDouble> numShared(nprocs, 0.0);
        // determine the number of communicators on this rank where
        // the share proc has a higher rank id
        numShared[rank] = nShrProc;
        m_comm->AllReduce(numShared, LibUtilities::ReduceSum);
        int totShared = (int)Vmath::Vsum(nprocs, numShared, 1);

        if (totShared)
        {
            int cnt = 0;
            Array<OneD, NekDouble> numOffset(nprocs, 0.0);
            for (auto &s : SharedProc)
            {
                if (s.second > rank)
                {
                    cnt++;
                }
            }

            numOffset[rank] = cnt;
            m_comm->AllReduce(numOffset, LibUtilities::ReduceMax);

            // make numShared into a cumulative list
            for (int i = 1; i < nprocs; ++i)
            {
                numShared[i] += numShared[i - 1];
                numOffset[i] += numOffset[i - 1];
            }
            for (int i = nprocs - 1; i > 0; --i)
            {
                numShared[i] = numShared[i - 1];
                numOffset[i] = numOffset[i - 1];
            }
            numShared[0] = 0;
            numOffset[0] = 0;

            Array<OneD, NekDouble> shareddom(totShared, -1.0);
            cnt = 0; // collect a list of domains that each shared commm is
                     // attached to
            for (auto &s : SharedProc)
            {
                shareddom[numShared[rank] + cnt] = s.first;
                ++cnt;
            }
            m_comm->AllReduce(shareddom, LibUtilities::ReduceMax);

            // define a numbering scheme
            Array<OneD, NekDouble> sharedid(totShared, -1.0);
            cnt      = 0;
            int cnt1 = 0;
            for (auto &s : SharedProc)
            {
                if (s.second > rank)
                {
                    sharedid[numShared[rank] + cnt] = numOffset[rank] + cnt1;

                    // find shared proc offset by matching the domain ids.
                    int j;
                    for (j = 0; j < maxCreateComm; ++j)
                    {
                        if ((numShared[s.second] + j < totShared) &&
                            (shareddom[numShared[s.second] + j] == s.first))
                        {
                            break;
                        }
                    }

                    sharedid[numShared[s.second] + j] = numOffset[rank] + cnt1;
                    cnt1++;
                }
                cnt++;
            }

            // now communicate ids to all procesors.
            m_comm->AllReduce(sharedid, LibUtilities::ReduceMax);

            if (rank == 0)
            {
                for (int i = 0; i < totShared; ++i)
                {
                    ASSERTL1(sharedid[i] != -1.0,
                             "Failed to number shared proc uniquely");
                }
            }

            cnt = 0;
            int maxoffset =
                (int)Vmath::Vmax(nShrProc, &sharedid[numShared[rank]], 1);
            m_comm->AllReduce(maxoffset, LibUtilities::ReduceMax);
            maxoffset++;

            // make a map relating the order of the SharedProc to the domain id
            map<int, int> ShrToDom;
            cnt = 0;
            for (auto &s : SharedProc)
            {
                ShrToDom[cnt] = s.first;
                ++cnt;
            }

            // Set up a set the domain ids listed in the ordering of
            // the SharedProc unique numbering
            set<int> doneDom;
            NekDouble maxdom = Vmath::Vmax(totShared, shareddom, 1);
            maxdom++;
            for (int i = 0; i < nShrProc; ++i)
            {
                int minId =
                    Vmath::Imin(nShrProc, &sharedid[numShared[rank]], 1);
                sharedid[numShared[rank] + minId] += maxoffset;

                m_domOrder.push_back(ShrToDom[minId]);
                doneDom.insert(ShrToDom[minId]);
            }

            // add all the other
            for (auto &d : m_domain)
            {
                if (doneDom.count(d.first) == 0)
                {
                    m_domOrder.push_back(d.first);
                }
            }
        }
    }
}

void PulseWaveSystem::SetUpDomainInterfaces(void)
{
    map<int, std::vector<InterfacePointShPtr>> VidToDomain;

    // First set up domain to determine how many vertices meet at a
    // vertex.  This allows us to

    /* Loop over domain and find out if we have any undefined
     * boundary conditions representing interfaces. If so make a
     * map based around vid and storing the domains that are
     * part of interfaces.
     */
    for (int omega = 0; omega < m_nDomains; ++omega)
    {
        int vesselID = omega * m_nVariables;

        for (int i = 0; i < (m_vessels[vesselID]->GetBndConditions()).size();
             ++i)
        {
            if (m_vessels[vesselID]->GetBndConditions()[i]->GetUserDefined() ==
                "Interface")
            {
                // Get Vid of interface
                int vid = m_vessels[vesselID]
                              ->UpdateBndCondExpansion(i)
                              ->GetExp(0)
                              ->GetGeom()
                              ->GetVid(0);

                MultiRegions::ExpListSharedPtr trace =
                    m_vessels[vesselID]->GetTrace();
                InterfacePointShPtr Ipt;

                bool finish = false;
                // find which elmt, the local vertex and the data
                // offset of point
                for (int n = 0; n < m_vessels[vesselID]->GetExpSize(); ++n)
                {
                    for (int p = 0; p < 2; ++p)
                    {

                        if (m_vessels[vesselID]
                                ->GetTraceMap()
                                ->GetElmtToTrace()[n][p]
                                ->as<LocalRegions::Expansion>()
                                ->GetGeom()
                                ->GetVid(0) == vid)
                        {
                            int eid = m_vessels[vesselID]
                                          ->GetTraceMap()
                                          ->GetElmtToTrace()[n][p]
                                          ->GetElmtId();

                            int tid = m_vessels[vesselID]
                                          ->GetTrace()
                                          ->GetCoeff_Offset(eid);

                            Ipt = MemoryManager<
                                InterfacePoint>::AllocateSharedPtr(vid, omega,
                                                                   n, p, tid,
                                                                   i);

                            finish = true;
                            break;
                        }
                    }
                    if (finish == true)
                    {
                        break;
                    }
                }

                VidToDomain[vid].push_back(Ipt);
            }
        }
    }

    // Set up comms for parallel cases
    map<int, int> domId;

    int cnt = 0;
    for (auto &d : m_domain)
    {
        domId[cnt] = d.first;
        ++cnt;
    }
    ASSERTL1(domId.size() == m_nDomains, "Number of domains do not match");

    Gs::gs_data *gscomm;
    int nvid2dom = VidToDomain.size();
    Array<OneD, long> tmp(nvid2dom);
    Array<OneD, NekDouble> nvid(nvid2dom, 0.0);

    cnt = 0;
    for (auto &v : VidToDomain)
    {
        tmp[cnt] = v.first + 1;
        for (int i = 0; i < v.second.size(); ++i)
        {
            nvid[cnt] += 1.0;
        }
        cnt++;
    }

    const bool verbose = m_session->DefinesCmdLineArgument("verbose");

    gscomm = Gs::Init(tmp, m_comm->GetRowComm(), verbose);
    Gs::Gather(nvid, Gs::gs_add, gscomm);

    // Loop over domains and identify how many vessels at a
    // bifurcation use the beginning of the element at the
    // junction. This distinguishes a bifurcation and merging junction
    // since we implicitly assume elements are ordered left ot right
    // and so at a bifurcation we have two elements meeting this point
    // at the first vertex
    Array<OneD, NekDouble> nbeg(nvid2dom, 0.0);
    cnt = 0;
    for (auto &v : VidToDomain)
    {
        if (nvid[cnt] == 3.0)
        {
            for (int i = 0; i < v.second.size(); ++i)
            {
                // for bifurcations and merging junctions store how
                // many points at interface are at the beginning or
                // end
                if (v.second[i]->m_elmtVert == 0)
                {
                    nbeg[cnt] += 1.0;
                }
            }
        }
        ++cnt;
    }
    Gs::Gather(nbeg, Gs::gs_add, gscomm);

    // Finally determine the maximum domain id between the two vessels
    // at an interfaces or the maximum domain id of the daughter
    // vessels for a bifurcation or the two parent id's for a merging
    // junction
    Array<OneD, NekDouble> dom(VidToDomain.size(), 0.0);
    cnt = 0;
    for (auto &v : VidToDomain)
    {
        if (nvid[cnt] == 2.0)
        {
            // store the domain id of the two domains
            for (int i = 0; i < v.second.size(); ++i)
            {
                dom[cnt] =
                    max(dom[cnt], (NekDouble)domId[v.second[i]->m_domain]);
            }
        }
        else if (nvid[cnt] == 3.0)
        {
            // store the domain id of the two daughter vessels if a
            // bifurcation or the two parent vessels if a merging
            // junction

            int val = (nbeg[cnt] == 2.0) ? 0 : 1;

            for (int i = 0; i < v.second.size(); ++i)
            {
                if (v.second[i]->m_elmtVert == val)
                {
                    dom[cnt] =
                        max(dom[cnt], (NekDouble)domId[v.second[i]->m_domain]);
                }
            }
        }
        ++cnt;
    }
    Gs::Gather(dom, Gs::gs_max, gscomm);

    // loop over map and set up Interface information;
    int v = 0;
    for (auto &iter : VidToDomain)
    {
        ASSERTL1(nvid[v] != 1.0, "Found an interface wth only "
                                 "one domain which should not happen");

        if (nvid[v] == 2.0) // Vessel jump interface
        {
            for (int i = 0; i < iter.second.size(); ++i)
            {
                if (domId[iter.second[i]->m_domain] == dom[v])
                {
                    iter.second[i]->m_riemannOrd = 1;
                }
                else
                {
                    iter.second[i]->m_riemannOrd = 0;
                }
            }
            m_vesselIntfcs.push_back(iter.second);
        }
        else if (nvid[v] == 3.0) // Bifurcation or Merging junction.
        {
            // Set up Bifurcation information
            int val = (nbeg[v] == 2.0) ? 1 : 0;

            for (int i = 0; i < iter.second.size(); ++i)
            {
                if (iter.second[i]->m_elmtVert == val)
                {
                    iter.second[i]->m_riemannOrd = 0;
                }
                else
                {
                    if (domId[iter.second[i]->m_domain] == dom[v])
                    {
                        iter.second[i]->m_riemannOrd = 2;
                    }
                    else
                    {
                        iter.second[i]->m_riemannOrd = 1;
                    }
                }
            }

            if (nbeg[v] == 2.0)
            {
                m_bifurcations.push_back(iter.second);
            }
            else
            {
                m_mergingJcts.push_back(iter.second);
            }
        }
        else
        {
            ASSERTL0(false, "Unknown junction type");
        }
        ++v; // incrementing loop over vertices
    }

    // finally set up a communicator for interface data collection
    Array<OneD, long> tmp1(3 * nvid2dom);
    if (nvid2dom)
    {
        Vmath::Zero(3 * nvid2dom, tmp1, 1);
    }

    cnt = 0;
    for (int n = 0; n < m_bifurcations.size(); ++n, ++cnt)
    {
        int vid = m_bifurcations[n][0]->m_vid;
        for (int i = 0; i < 3; ++i)
        {
            tmp1[3 * cnt + i] = (NekDouble)(3 * vid + i + 1);
        }
    }

    for (int n = 0; n < m_mergingJcts.size(); ++n, ++cnt)
    {
        int vid = m_mergingJcts[n][0]->m_vid;
        for (int i = 0; i < 3; ++i)
        {
            int offset        = m_mergingJcts[n][i]->m_riemannOrd;
            tmp1[3 * cnt + i] = (NekDouble)(3 * vid + i + 1);
        }
    }

    for (int n = 0; n < m_vesselIntfcs.size(); ++n, ++cnt)
    {
        int vid = m_vesselIntfcs[n][0]->m_vid;
        for (int i = 0; i < 3; ++i)
        {
            tmp1[3 * cnt + i] = (NekDouble)(3 * vid + i + 1);
        }
    }

    m_intComm = Gs::Init(tmp1, m_comm->GetRowComm(), verbose);
}

void PulseWaveSystem::SetUpDomainInterfaceBCs(
    SpatialDomains::BoundaryConditions &Allbcs,
    map<int, LibUtilities::CommSharedPtr> &domComm)
{

    const SpatialDomains::BoundaryRegionCollection &bregions =
        Allbcs.GetBoundaryRegions();

    /* Loop over domain and find out if we have any undefined
     * boundary conditions representing interfaces. If so make a
     * map based around vid and storing the domains that are
     * part of interfaces.
     */
    int numNewBc = 1;
    for (auto &d : m_domOrder)
    {

        map<int, SpatialDomains::GeometrySharedPtr> domvids;

        // Loop over each domain and find which vids are only touched
        // by one element
        for (auto &compIt : m_domain[d])
        {
            for (int j = 0; j < compIt.second->m_geomVec.size(); ++j)
            {

                // get hold of vids of each segment
                for (int p = 0; p < 2; ++p)
                {
                    SpatialDomains::GeometrySharedPtr vert =
                        compIt.second->m_geomVec[j]->GetVertex(p);
                    int vid = vert->GetGlobalID();

                    // if vid has already been added remove it otherwise add
                    // into set
                    if (domvids.count(vid))
                    {
                        domvids.erase(vid);
                    }
                    else
                    {
                        domvids[vid] = vert;
                    }
                }
            }
        }

        ASSERTL1(domvids.size() == 2,
                 "Failed to find two end points "
                 "of a domain (domvids = " +
                     boost::lexical_cast<std::string>(domvids.size()) + ")");

        int nprocs = domComm[d]->GetSize();

        if (nprocs > 1) // Remove parallel interfaces
        {
            int rank = domComm[d]->GetRank();
            Array<OneD, int> nvids(nprocs, 0);
            nvids[rank] = domvids.size();
            domComm[d]->AllReduce(nvids, LibUtilities::ReduceSum);

            int totvids = Vmath::Vsum(nprocs, nvids, 1);

            Array<OneD, int> locids(totvids, -1);

            int cnt = 0;
            for (int i = 0; i < rank; ++i)
            {
                cnt += nvids[i];
            }

            for (auto &i : domvids)
            {
                locids[cnt++] = i.first;
            }

            // collect lcoal vids on this domain on all processor
            domComm[d]->AllReduce(locids, LibUtilities::ReduceMax);

            set<int> chkvids;
            for (int i = 0; i < totvids; ++i)
            {
                if (chkvids.count(locids[i]))
                {
                    // if this id is on local domain then remove
                    if (domvids.count(locids[i]))
                    {
                        domvids.erase(locids[i]);
                    }
                }
                else
                {
                    chkvids.insert(locids[i]);
                }
            }
        }

        // Finally check if remaining domvids are already declared, if
        // not add them

        // Erase from domvids any existing BCs
        for (auto &it : bregions)
        {
            for (auto &bregionIt : *it.second)
            {
                // can assume that all regions only contain one point in 1D
                // Really do not need loop above
                int id = bregionIt.second->m_geomVec[0]->GetGlobalID();

                if (domvids.count(id))
                {
                    domvids.erase(id);
                }
            }
        }

        // Finally add interface conditions to Allbcs
        std::vector<std::string> variables = m_session->GetVariables();

        // determine the maxmimum bregion index
        int maxbregind = 0;
        for (auto &b : bregions)
        {
            maxbregind = (maxbregind > b.first) ? maxbregind : b.first;
        }

        for (auto &b : domvids)
        {
            SpatialDomains::BoundaryRegionShPtr breg(
                MemoryManager<
                    SpatialDomains::BoundaryRegion>::AllocateSharedPtr());

            // Set up Composite (GemetryVector) to contain vertex and put into
            // bRegion
            SpatialDomains::CompositeSharedPtr gvec =
                MemoryManager<SpatialDomains::Composite>::AllocateSharedPtr();
            gvec->m_geomVec.push_back(b.second);
            (*breg)[b.first] = gvec;

            Allbcs.AddBoundaryRegions(maxbregind + numNewBc, breg);

            SpatialDomains::BoundaryConditionMapShPtr bCondition =
                MemoryManager<
                    SpatialDomains::BoundaryConditionMap>::AllocateSharedPtr();

            std::string userDefined = "Interface";
            // Set up just boundary condition for this variable.
            SpatialDomains::BoundaryConditionShPtr DirichletInterface(
                MemoryManager<SpatialDomains::DirichletBoundaryCondition>::
                    AllocateSharedPtr(m_session, "0", userDefined));

            for (int i = 0; i < variables.size(); ++i)
            {
                (*bCondition)[variables[i]] = DirichletInterface;
            }

            Allbcs.AddBoundaryConditions(maxbregind + numNewBc, bCondition);
            ++numNewBc;
        }
    }
}

/**
 *  Initialisation routine for multiple subdomain case. Sets the
 *  initial conditions for all arterial subdomains read from the
 *  inputfile. Sets the material properties and the A_0 area for
 *  all subdomains and fills the domain-linking boundary
 *  conditions with the initial values of their domain.
 */
void PulseWaveSystem::v_DoInitialise()
{

    if (m_session->GetComm()->GetRank() == 0)
    {
        cout << "Initial Conditions: " << endl;
    }

    /* Loop over all subdomains to initialize all with the Initial
     * Conditions read from the inputfile*/
    int omega = 0;
    for (auto &d : m_domOrder)
    {
        for (int i = 0; i < 2; ++i)
        {
            m_fields[i] = m_vessels[m_nVariables * omega + i];
        }

        if (m_session->GetComm()->GetRank() == 0)
        {
            cout << "Subdomain: " << omega << endl;
        }

        SetInitialConditions(0.0, 0, d);
        omega++;
    }

    // Reset 2 variables to first vessels
    for (int i = 0; i < 2; ++i)
    {
        m_fields[i] = m_vessels[i];
    }
}

/**
 * NEEDS Updating:
 *
 *  DoSolve routine for PulseWavePropagation with multiple
 *  subdomains taken from UnsteadySystem and modified for
 *  multidomain case. Initialises the time integration scheme (as
 *  specified in the session file), and perform the time
 *  integration. Within the timestepping loop the following is
 *  done: 1. Link all arterial segments according to the network
 *  structure, solve the Riemann problem between different
 *  arterial segments and assign the values to the boundary
 *  conditions (LinkSubdomains) 2. Every arterial segment is
 *  solved independently for this timestep. This is done by handing
 *  the solution vector \f$ \mathbf{u} \f$ and the right hand side
 *  m_ode, which is the PulseWavePropagation class in this example
 *  over to the time integration scheme
 */
void PulseWaveSystem::v_DoSolve()
{
    NekDouble IntegrationTime = 0.0;
    int i;
    int n;
    int nchk = 1;

    Array<OneD, Array<OneD, NekDouble>> fields(m_nVariables);

    for (int i = 0; i < m_nVariables; ++i)
    {
        fields[i] = m_vessels[i]->UpdatePhys();
        m_fields[i]->SetPhysState(false);
    }

    m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);

    // Time loop
    for (n = 0; n < m_steps; ++n)
    {
        LibUtilities::Timer time_v_IntStep;
        time_v_IntStep.Start();
        fields = m_intScheme->TimeIntegrate(n, m_timestep, m_ode);
        m_time += m_timestep;
        time_v_IntStep.Stop();
        time_v_IntStep.AccumulateRegion("PulseWaveSystem::TimeIntegrationStep",
                                        0);
        // IntegrationTime += timer.TimePerTest(1);

        // Write out status information.
        if (m_session->GetComm()->GetRank() == 0 && !((n + 1) % m_infosteps))
        {
            cout << "Steps: " << n + 1 << "\t Time: " << m_time
                 << "\t Time-step: " << m_timestep << "\t" << endl;
        }

        // Transform data if needed
        if (!((n + 1) % m_checksteps))
        {
            for (i = 0; i < m_nVariables; ++i)
            {
                int cnt = 0;
                for (int omega = 0; omega < m_nDomains; omega++)
                {
                    m_vessels[omega * m_nVariables + i]->FwdTrans(
                        fields[i] + cnt,
                        m_vessels[omega * m_nVariables + i]->UpdateCoeffs());
                    cnt += m_vessels[omega * m_nVariables + i]->GetTotPoints();
                }
            }
            CheckPoint_Output(nchk++);
        }

    } // end of timeintegration

    // Copy Array To Vessel Phys Fields
    for (int i = 0; i < m_nVariables; ++i)
    {
        Vmath::Vcopy(fields[i].size(), fields[i], 1, m_vessels[i]->UpdatePhys(),
                     1);
    }

    /**  if(m_session->GetComm()->GetRank() == 0)
       {
           cout << "Time-integration timing: " << IntegrationTime << " s" <<
       endl
                << endl;
       } **/
}

void PulseWaveSystem::FillDataFromInterfacePoint(
    InterfacePointShPtr &I,
    const Array<OneD, const Array<OneD, NekDouble>> &fields, NekDouble &A,
    NekDouble &u, NekDouble &beta, NekDouble &A_0, NekDouble &alpha)
{
    LibUtilities::Timer time_FillDataFromInterfacePoint;
    time_FillDataFromInterfacePoint.Start();

    int omega       = I->m_domain;
    int traceId     = I->m_traceId;
    int eid         = I->m_elmt;
    int vert        = I->m_elmtVert;
    int vesselID    = omega * m_nVariables;
    int phys_offset = m_vessels[vesselID]->GetPhys_Offset(eid);
    LocalRegions::ExpansionSharedPtr dumExp;
    Array<OneD, NekDouble> Tmp(1);

    m_vessels[vesselID]->GetExp(eid)->GetTracePhysVals(
        vert, dumExp, fields[0] + m_fieldPhysOffset[omega] + phys_offset, Tmp);
    A = Tmp[0];
    m_vessels[vesselID]->GetExp(eid)->GetTracePhysVals(
        vert, dumExp, fields[1] + m_fieldPhysOffset[omega] + phys_offset, Tmp);
    u = Tmp[0];

    beta  = m_beta_trace[omega][traceId];
    A_0   = m_A_0_trace[omega][traceId];
    alpha = m_alpha_trace[omega][traceId];
    time_FillDataFromInterfacePoint.Stop();
    time_FillDataFromInterfacePoint.AccumulateRegion(
        "PulseWaveSystem::FillDataFromInterfacePoint", 2);
}

void PulseWaveSystem::EnforceInterfaceConditions(
    const Array<OneD, const Array<OneD, NekDouble>> &fields)
{
    LibUtilities::Timer time_EnforceInterfaceConditions;
    time_EnforceInterfaceConditions.Start();
    int dom, bcpos;

    int totif =
        m_bifurcations.size() + m_mergingJcts.size() + m_vesselIntfcs.size();

    Array<OneD, NekDouble> Aut, Au(3 * totif, 0.0);
    Array<OneD, NekDouble> uut, uu(3 * totif, 0.0);
    Array<OneD, NekDouble> betat, beta(3 * totif, 0.0);
    Array<OneD, NekDouble> A_0t, A_0(3 * totif, 0.0);
    Array<OneD, NekDouble> alphat, alpha(3 * totif, 0.0);

    // Bifurcations Data:
    int cnt = 0;
    for (int n = 0; n < m_bifurcations.size(); ++n, ++cnt)
    {
        for (int i = 0; i < m_bifurcations[n].size(); ++i)
        {
            int l = m_bifurcations[n][i]->m_riemannOrd;
            FillDataFromInterfacePoint(
                m_bifurcations[n][i], fields, Au[l + 3 * cnt], uu[l + 3 * cnt],
                beta[l + 3 * cnt], A_0[l + 3 * cnt], alpha[l + 3 * cnt]);
        }
    }

    // Enforce Merging vessles Data:
    for (int n = 0; n < m_mergingJcts.size(); ++n, ++cnt)
    {
        // Merged vessel
        for (int i = 0; i < m_mergingJcts.size(); ++i)
        {
            int l = m_mergingJcts[n][i]->m_riemannOrd;
            FillDataFromInterfacePoint(
                m_mergingJcts[n][i], fields, Au[l + 3 * cnt], uu[l + 3 * cnt],
                beta[l + 3 * cnt], A_0[l + 3 * cnt], alpha[l + 3 * cnt]);
        }
    }

    // Enforce interface junction between two vessesls
    for (int n = 0; n < m_vesselIntfcs.size(); ++n, ++cnt)
    {
        for (int i = 0; i < m_vesselIntfcs[n].size(); ++i)
        {
            int l = m_vesselIntfcs[n][i]->m_riemannOrd;
            FillDataFromInterfacePoint(
                m_vesselIntfcs[n][i], fields, Au[l + 3 * cnt], uu[l + 3 * cnt],
                beta[l + 3 * cnt], A_0[l + 3 * cnt], alpha[l + 3 * cnt]);
        }
    }

    // Gather data if running in parallel
    Gs::Gather(Au, Gs::gs_add, m_intComm);
    Gs::Gather(uu, Gs::gs_add, m_intComm);
    Gs::Gather(beta, Gs::gs_add, m_intComm);
    Gs::Gather(A_0, Gs::gs_add, m_intComm);
    Gs::Gather(alpha, Gs::gs_add, m_intComm);

    // Enforce Bifurcations:
    cnt = 0;
    for (int n = 0; n < m_bifurcations.size(); ++n, ++cnt)
    {
        // Solve the Riemann problem for a bifurcation
        BifurcationRiemann(Aut = Au + 3 * cnt, uut = uu + 3 * cnt,
                           betat = beta + 3 * cnt, A_0t = A_0 + 3 * cnt,
                           alphat = alpha + 3 * cnt);

        // Store the values into the right positions:
        for (int i = 0; i < m_bifurcations[n].size(); ++i)
        {
            dom   = m_bifurcations[n][i]->m_domain;
            bcpos = m_bifurcations[n][i]->m_bcPosition;
            int l = m_bifurcations[n][i]->m_riemannOrd;
            m_vessels[dom * m_nVariables]
                ->UpdateBndCondExpansion(bcpos)
                ->UpdatePhys()[0] = Au[l + 3 * cnt];
            m_vessels[dom * m_nVariables + 1]
                ->UpdateBndCondExpansion(bcpos)
                ->UpdatePhys()[0] = uu[l + 3 * cnt];
        }
    }

    // Enforce Merging vessles;
    for (int n = 0; n < m_mergingJcts.size(); ++n, ++cnt)
    {
        // Solve the Riemann problem for a merging vessel
        MergingRiemann(Aut = Au + 3 * cnt, uut = uu + 3 * cnt,
                       betat = beta + 3 * cnt, A_0t = A_0 + 3 * cnt,
                       alphat = alpha + 3 * cnt);

        // Store the values into the right positions:
        for (int i = 0; i < m_mergingJcts[n].size(); ++i)
        {
            int dom   = m_mergingJcts[n][i]->m_domain;
            int bcpos = m_mergingJcts[n][i]->m_bcPosition;
            int l     = m_mergingJcts[n][i]->m_riemannOrd;
            m_vessels[dom * m_nVariables]
                ->UpdateBndCondExpansion(bcpos)
                ->UpdatePhys()[0] = Au[l + 3 * cnt];
            m_vessels[dom * m_nVariables + 1]
                ->UpdateBndCondExpansion(bcpos)
                ->UpdatePhys()[0] = uu[l + 3 * cnt];
        }
    }

    // Enforce interface junction between two vessesls
    for (int n = 0; n < m_vesselIntfcs.size(); ++n, ++cnt)
    {
        InterfaceRiemann(Aut = Au + 3 * cnt, uut = uu + 3 * cnt,
                         betat = beta + 3 * cnt, A_0t = A_0 + 3 * cnt,
                         alphat = alpha + 3 * cnt);

        // Store the values into the right positions:
        for (int i = 0; i < m_vesselIntfcs[n].size(); ++i)
        {
            int dom   = m_vesselIntfcs[n][i]->m_domain;
            int bcpos = m_vesselIntfcs[n][i]->m_bcPosition;
            int l     = m_vesselIntfcs[n][i]->m_riemannOrd;
            m_vessels[dom * m_nVariables]
                ->UpdateBndCondExpansion(bcpos)
                ->UpdatePhys()[0] = Au[l + 3 * cnt];
            m_vessels[dom * m_nVariables + 1]
                ->UpdateBndCondExpansion(bcpos)
                ->UpdatePhys()[0] = uu[l + 3 * cnt];
        }
    }
    time_EnforceInterfaceConditions.Stop();
    time_EnforceInterfaceConditions.AccumulateRegion(
        "PulseWaveSystem::EnforceInterfaceConditions", 1);
}

/**
 *  Solves the Riemann problem at a bifurcation by assuming
 *  subsonic flow at both sides of the boundary and by
 *  applying conservation of mass and continuity of the total
 *  pressure \f$ \frac{p}{rho} + \frac{u^{2}}{2}. \f$ The
 *  other 3 missing equations come from the characteristic
 *  variables. For further information see "Pulse
 *  WavePropagation in the human vascular system" Section
 *  3.4.4
 */
void PulseWaveSystem::BifurcationRiemann(Array<OneD, NekDouble> &Au,
                                         Array<OneD, NekDouble> &uu,
                                         Array<OneD, NekDouble> &beta,
                                         Array<OneD, NekDouble> &A_0,
                                         Array<OneD, NekDouble> &alpha)
{

    NekDouble rho = m_rho;
    Array<OneD, NekDouble> W(3);
    Array<OneD, NekDouble> P_Au(3);
    Array<OneD, NekDouble> W_Au(3);
    NekMatrix<NekDouble> invJ(6, 6);
    NekVector<NekDouble> f(6);
    NekVector<NekDouble> dx(6);

    int proceed  = 1;
    int iter     = 0;
    int MAX_ITER = 100;

    // Forward characteristic
    m_pressureArea->GetW1(W[0], uu[0], beta[0], Au[0], A_0[0], alpha[0]);

    // Backward characteristics
    for (int i = 1; i < 3; ++i)
    {
        m_pressureArea->GetW2(W[i], uu[i], beta[i], Au[i], A_0[i], alpha[i]);
    }

    // Tolerances for the algorithm
    NekDouble Tol = 1.0E-10;

    // Newton Iteration
    while ((proceed) && (iter < MAX_ITER))
    {
        LibUtilities::Timer time_BifurcationRiemann;
        time_BifurcationRiemann.Start();
        iter += 1;

        /*
         * We solve the six constraint equations via a multivariate Newton
         * iteration. Equations are:
         * 1. Forward characteristic:          W1(A_L,  U_L)  = W1(Au_L,  Uu_L)
         * 2. Backward characteristic 1:       W2(A_R1, U_R1) = W2(Au_R1, Uu_R1)
         * 3. Backward characteristic 2:       W2(A_R2, U_R2) = W2(Au_R2, Uu_R2)
         * 4. Conservation of mass:            Au_L * Uu_L    = Au_R1 * Uu_R1 +
         *                                                      Au_R2 * Uu_R2
         * 5. Continuity of total pressure 1:  rho * Uu_L  * Uu_L  / 2 + p(Au_L)
         * = rho * Uu_R1 * Uu_R1 / 2 + p(Au_R1)
         * 6. Continuity of total pressure 2:  rho * Uu_L  * Uu_L  / 2 + p(Au_L)
         * = rho * Uu_R2 * Uu_R2 / 2 + p(Au_R2)
         */
        for (int i = 0; i < 3; ++i)
        {
            m_pressureArea->GetPressure(P_Au[i], beta[i], Au[i], A_0[i], 0, 0,
                                        alpha[i]);
        }

        m_pressureArea->GetW1(W_Au[0], uu[0], beta[0], Au[0], A_0[0], alpha[0]);
        m_pressureArea->GetW2(W_Au[1], uu[1], beta[1], Au[1], A_0[1], alpha[1]);
        m_pressureArea->GetW2(W_Au[2], uu[2], beta[2], Au[2], A_0[2], alpha[2]);

        // Constraint equations set to zero
        f[0] = W_Au[0] - W[0];
        f[1] = W_Au[1] - W[1];
        f[2] = W_Au[2] - W[2];
        f[3] = Au[0] * uu[0] - Au[1] * uu[1] - Au[2] * uu[2];
        f[4] = uu[0] * uu[0] + 2 * P_Au[0] / rho - uu[1] * uu[1] -
               2 * P_Au[1] / rho;
        f[5] = uu[0] * uu[0] + 2 * P_Au[0] / rho - uu[2] * uu[2] -
               2 * P_Au[2] / rho;

        // Inverse Jacobian for x + dx = x - J^(-1) * f(x)
        m_pressureArea->GetJacobianInverse(invJ, Au, uu, beta, A_0, alpha,
                                           "Bifurcation");

        Multiply(dx, invJ, f);

        // Update the solution: x_new = x_old - dx
        for (int i = 0; i < 3; ++i)
        {
            uu[i] -= dx[i];
            Au[i] -= dx[i + 3];
        }

        // Check if the error of the solution is smaller than Tol
        if (Dot(dx, dx) < Tol)
        {
            proceed = 0;
        }

        // Check if solver converges
        if (iter >= MAX_ITER)
        {
            ASSERTL0(false, "Riemann solver for Bifurcation did not converge.");
        }
        time_BifurcationRiemann.Stop();
        time_BifurcationRiemann.AccumulateRegion(
            "PulseWaveSystem::Bifurcation Riemann", 2);
    }
}

/**
 *  Solves the Riemann problem at an merging flow condition by
 *  assuming subsonic flow at both sides of the boundary and by
 *  applying conservation of mass and continuity of the total
 *  pressure \f$ \frac{p}{rho} + \frac{u^{2}}{2}. \f$ The other 3
 *  missing equations come from the characteristic variables. For
 *  further information see "Pulse WavePropagation in the human
 *  vascular system" Section 3.4.4
 */
void PulseWaveSystem::MergingRiemann(Array<OneD, NekDouble> &Au,
                                     Array<OneD, NekDouble> &uu,
                                     Array<OneD, NekDouble> &beta,
                                     Array<OneD, NekDouble> &A_0,
                                     Array<OneD, NekDouble> &alpha)
{

    NekDouble rho = m_rho;
    Array<OneD, NekDouble> W(3);
    Array<OneD, NekDouble> W_Au(3);
    Array<OneD, NekDouble> P_Au(3);
    NekMatrix<NekDouble> invJ(6, 6);
    NekVector<NekDouble> f(6);
    NekVector<NekDouble> dx(6);

    int proceed  = 1;
    int iter     = 0;
    int MAX_ITER = 15;

    // Backward characteristic
    m_pressureArea->GetW2(W[0], uu[0], beta[0], Au[0], A_0[0], alpha[0]);

    // Forward characteristics
    for (int i = 1; i < 3; ++i)
    {
        m_pressureArea->GetW1(W[i], uu[i], beta[i], Au[i], A_0[i], alpha[i]);
    }

    // Tolerances for the algorithm
    NekDouble Tol = 1.0E-10;

    // Newton Iteration
    while ((proceed) && (iter < MAX_ITER))
    {
        LibUtilities::Timer time_MergingRiemann;
        time_MergingRiemann.Start();
        iter += 1;

        /*
         * We solve the six constraint equations via a multivariate Newton
         * iteration. Equations are:
         * 1. Backward characteristic:          W2(A_R,  U_R)  = W1(Au_R, Uu_R)
         * 2. Forward characteristic 1:         W1(A_L1, U_L1) = W1(Au_L1,
         * Uu_R1)
         * 3. Forward characteristic 2:         W1(A_L2, U_L2) = W1(Au_L2,
         * Uu_L2)
         * 4. Conservation of mass:             Au_R * Uu_R    = Au_L1 * Uu_L1 +
         *                                                       Au_L2 * Uu_L2
         * 5. Continuity of total pressure 1:  rho * Uu_R  * Uu_R  / 2 + p(Au_R)
         * = rho * Uu_L1 * Uu_L1 / 2 + p(Au_L1)
         * 6. Continuity of total pressure 2:  rho * Uu_R  * Uu_R  / 2 + p(Au_R)
         * = rho * Uu_L2 * Uu_L2 / 2 + p(Au_L2)
         */
        for (int i = 0; i < 3; ++i)
        {
            m_pressureArea->GetPressure(P_Au[i], beta[i], Au[i], A_0[i], 0, 0,
                                        alpha[i]);
        }

        m_pressureArea->GetW2(W_Au[0], uu[0], beta[0], Au[0], A_0[0], alpha[0]);
        m_pressureArea->GetW1(W_Au[1], uu[1], beta[1], Au[1], A_0[1], alpha[1]);
        m_pressureArea->GetW1(W_Au[2], uu[2], beta[2], Au[2], A_0[2], alpha[2]);

        // Constraint equations set to zero
        f[0] = W_Au[0] - W[0];
        f[1] = W_Au[1] - W[1];
        f[2] = W_Au[2] - W[2];
        f[3] = Au[0] * uu[0] - Au[1] * uu[1] - Au[2] * uu[2];
        f[4] = uu[0] * uu[0] + 2 * P_Au[0] / rho - uu[1] * uu[1] -
               2 * P_Au[1] / rho;
        f[5] = uu[0] * uu[0] + 2 * P_Au[0] / rho - uu[2] * uu[2] -
               2 * P_Au[2] / rho;

        // Inverse Jacobian for x + dx = x - J^(-1) * f(x)
        m_pressureArea->GetJacobianInverse(invJ, Au, uu, beta, A_0, alpha,
                                           "Merge");

        Multiply(dx, invJ, f);

        // Update the solution: x_new = x_old - dx
        for (int i = 0; i < 3; ++i)
        {
            uu[i] -= dx[i];
            Au[i] -= dx[i + 3];
        }

        // Check if the error of the solution is smaller than Tol
        if (Dot(dx, dx) < Tol)
        {
            proceed = 0;
        }

        // Check if solver converges
        if (iter >= MAX_ITER)
        {
            ASSERTL0(false, "Riemann solver for Merging Flow did not converge");
        }
        time_MergingRiemann.Stop();
        time_MergingRiemann.AccumulateRegion("PulseWaveSystem::MergingRiemann",
                                             2);
    }
}

/**
 *  Solves the Riemann problem at an interdomain
 *  junction/Interface by assuming subsonic flow at both sides of
 *  the boundary and by applying conservation of mass and
 *  continuity of the total pressure \f$ \frac{p}{rho} +
 *  \frac{u^{2}}{2}. \f$ The other 2 missing equations come from
 *  the characteristic variables. For further information see
 *  "Pulse WavePropagation in the human vascular system" Section
 *  3.4.
 */
void PulseWaveSystem::InterfaceRiemann(Array<OneD, NekDouble> &Au,
                                       Array<OneD, NekDouble> &uu,
                                       Array<OneD, NekDouble> &beta,
                                       Array<OneD, NekDouble> &A_0,
                                       Array<OneD, NekDouble> &alpha)
{

    NekDouble rho = m_rho;
    Array<OneD, NekDouble> W(2);
    Array<OneD, NekDouble> W_Au(2);
    Array<OneD, NekDouble> P_Au(2);
    NekMatrix<NekDouble> invJ(4, 4);
    NekVector<NekDouble> f(4);
    NekVector<NekDouble> dx(4);

    int proceed   = 1;
    int iter      = 0;
    int MAX_ITER  = 15;
    NekDouble Tol = 1.0E-10;

    // Forward and backward characteristics
    m_pressureArea->GetW1(W[0], uu[0], beta[0], Au[0], A_0[0], alpha[0]);
    m_pressureArea->GetW2(W[1], uu[1], beta[1], Au[1], A_0[1], alpha[1]);

    while ((proceed) && (iter < MAX_ITER))
    {
        LibUtilities::Timer time_InterfaceRiemann;
        time_InterfaceRiemann.Start();
        iter += 1;

        /*
         * We solve the four constraint equations via a multivariate Newton
         * iteration. Equations are:
         * 1. Forward characteristic:        W1(A_L, U_L) = W1(Au_L, Uu_L)
         * 2. Backward characteristic:       W2(A_R, U_R) = W2(Au_R, Uu_R)
         * 3. Conservation of mass:          Au_L * Uu_L = Au_R * Uu_R
         * 4. Continuity of total pressure:  rho * Uu_L * Uu_L / 2 + p(Au_L) =
         *                                   rho * Uu_R * Uu_R / 2 + p(Au_R)
         */
        for (int i = 0; i < 2; ++i)
        {
            m_pressureArea->GetPressure(P_Au[i], beta[i], Au[i], A_0[i], 0, 0,
                                        alpha[i]);
        }

        m_pressureArea->GetW1(W_Au[0], uu[0], beta[0], Au[0], A_0[0], alpha[0]);
        m_pressureArea->GetW2(W_Au[1], uu[1], beta[1], Au[1], A_0[1], alpha[1]);

        // Constraint equations set to zero
        f[0] = W_Au[0] - W[0];
        f[1] = W_Au[1] - W[1];
        f[2] = Au[0] * uu[0] - Au[1] * uu[1];
        f[3] = uu[0] * uu[0] + 2 * P_Au[0] / rho - uu[1] * uu[1] -
               2 * P_Au[1] / rho;

        // Inverse Jacobian for x + dx = x - J^(-1) * f(x)
        m_pressureArea->GetJacobianInverse(invJ, Au, uu, beta, A_0, alpha,
                                           "Interface");

        Multiply(dx, invJ, f);

        // Update solution: x_new = x_old - dx
        for (int i = 0; i < 2; ++i)
        {
            uu[i] -= dx[i];
            Au[i] -= dx[i + 2];
        }

        // Check if the error of the solution is smaller than Tol.
        if (Dot(dx, dx) < Tol)
        {
            proceed = 0;
        }
        time_InterfaceRiemann.Stop();
        time_InterfaceRiemann.AccumulateRegion(
            "PulseWaveSystem::InterfaceRiemann", 2);
    }

    if (iter >= MAX_ITER)
    {
        ASSERTL0(false, "Riemann solver for Interface did not converge");
    }
}

/**
 *  Writes the .fld file at the end of the simulation. Similar to the normal
 *  v_Output however the Multidomain output has to be prepared.
 */
void PulseWaveSystem::v_Output(void)
{
    /**
     * Write the field data to file. The file is named according to the session
     * name with the extension .fld appended.
     */
    std::string outname = m_sessionName + ".fld";

    WriteVessels(outname);
}

/**
 *  Writes the .fld file at the end of the simulation. Similar to the normal
 *  v_Output however the Multidomain output has to be prepared.
 */
void PulseWaveSystem::CheckPoint_Output(const int n)
{
    std::stringstream outname;
    outname << m_sessionName << "_" << n << ".chk";

    WriteVessels(outname.str());
}

/**
 * Writes the field data to a file with the given filename.
 * @param   outname         Filename to write to.
 */
void PulseWaveSystem::WriteVessels(const std::string &outname)
{
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
    std::vector<std::string> variables = m_session->GetVariables();

    for (int n = 0; n < m_nDomains; ++n)
    {
        m_vessels[n * m_nVariables]->GetFieldDefinitions(FieldDef);
    }

    std::vector<std::vector<NekDouble>> FieldData(FieldDef.size());

    int nFieldDefPerDomain = FieldDef.size() / m_nDomains;
    int cnt;

    // Copy Data into FieldData and set variable
    for (int n = 0; n < m_nDomains; ++n)
    {
        // Outputs area and velocity
        for (int j = 0; j < m_nVariables; ++j)
        {
            for (int i = 0; i < nFieldDefPerDomain; ++i)
            {
                cnt = n * nFieldDefPerDomain + i;

                FieldDef[cnt]->m_fields.push_back(variables[j]);

                m_vessels[n * m_nVariables]->AppendFieldData(
                    FieldDef[cnt], FieldData[cnt],
                    m_vessels[n * m_nVariables + j]->UpdateCoeffs());
            }
        }

        // Outputs pressure
        Array<OneD, NekDouble> PFwd(m_vessels[n * m_nVariables]->GetNcoeffs());

        m_vessels[n * m_nVariables]->FwdTransLocalElmt(m_pressure[n], PFwd);

        for (int i = 0; i < nFieldDefPerDomain; ++i)
        {
            cnt = n * nFieldDefPerDomain + i;

            FieldDef[cnt]->m_fields.push_back("P");

            m_vessels[n * m_nVariables]->AppendFieldData(FieldDef[cnt],
                                                         FieldData[cnt], PFwd);
        }

        if (extraFields)
        {
            Array<OneD, NekDouble> PWVFwd(
                m_vessels[n * m_nVariables]->GetNcoeffs());
            Array<OneD, NekDouble> W1Fwd(
                m_vessels[n * m_nVariables]->GetNcoeffs());
            Array<OneD, NekDouble> W2Fwd(
                m_vessels[n * m_nVariables]->GetNcoeffs());

            m_vessels[n * m_nVariables]->FwdTransLocalElmt(m_PWV[n], PWVFwd);
            m_vessels[n * m_nVariables]->FwdTransLocalElmt(m_W1[n], W1Fwd);
            m_vessels[n * m_nVariables]->FwdTransLocalElmt(m_W2[n], W2Fwd);

            for (int i = 0; i < nFieldDefPerDomain; ++i)
            {
                cnt = n * nFieldDefPerDomain + i;

                FieldDef[cnt]->m_fields.push_back("c");
                FieldDef[cnt]->m_fields.push_back("W1");
                FieldDef[cnt]->m_fields.push_back("W2");

                m_vessels[n * m_nVariables]->AppendFieldData(
                    FieldDef[cnt], FieldData[cnt], PWVFwd);
                m_vessels[n * m_nVariables]->AppendFieldData(
                    FieldDef[cnt], FieldData[cnt], W1Fwd);
                m_vessels[n * m_nVariables]->AppendFieldData(
                    FieldDef[cnt], FieldData[cnt], W2Fwd);
            }
        }
    }

    // Update time in field info if required
    if (m_fieldMetaDataMap.find("Time") != m_fieldMetaDataMap.end())
    {
        m_fieldMetaDataMap["Time"] = boost::lexical_cast<std::string>(m_time);
    }

    m_fld->Write(outname, FieldDef, FieldData, m_fieldMetaDataMap);
}

/* Compute the error in the L2-norm
 * @param   field           The field to compare.
 * @param   exactsoln       The exact solution to compare with.
 * @param   Normalised      Normalise L2-error.
 * @returns                 Error in the L2-norm.
 */
NekDouble PulseWaveSystem::v_L2Error(unsigned int field,
                                     const Array<OneD, NekDouble> &exact,
                                     bool Normalised)
{
    NekDouble L2error = 0.0;
    NekDouble L2error_dom;
    NekDouble Vol = 0.0;

    if (m_NumQuadPointsError == 0)
    {
        for (int omega = 0; omega < m_nDomains; omega++)
        {
            int vesselid = field + omega * m_nVariables;

            // since exactsoln is passed for just the first field size we need
            // to reset it each domain
            Array<OneD, NekDouble> exactsoln(m_vessels[vesselid]->GetNpoints());
            m_fields[field] = m_vessels[omega * m_nVariables];
            EvaluateExactSolution(field, exactsoln, m_time);

            if (m_vessels[vesselid]->GetPhysState() == false)
            {
                m_vessels[vesselid]->BwdTrans(
                    m_vessels[vesselid]->GetCoeffs(),
                    m_vessels[vesselid]->UpdatePhys());
            }

            if (exactsoln.size())
            {
                L2error_dom = m_vessels[vesselid]->L2(
                    m_vessels[vesselid]->GetPhys(), exactsoln);
            }
            else if (m_session->DefinesFunction("ExactSolution"))
            {
                Array<OneD, NekDouble> exactsoln(
                    m_vessels[vesselid]->GetNpoints());

                LibUtilities::EquationSharedPtr vEqu =
                    m_session->GetFunction("ExactSolution", field, omega);
                GetFunction("ExactSolution")
                    ->Evaluate(m_session->GetVariable(field), exactsoln,
                               m_time);

                L2error_dom = m_vessels[vesselid]->L2(
                    m_vessels[vesselid]->GetPhys(), exactsoln);
            }
            else
            {
                L2error_dom =
                    m_vessels[vesselid]->L2(m_vessels[vesselid]->GetPhys());
            }

            if (m_vessels[vesselid]->GetComm()->GetRank())
            {
                // ensure domains are only summed on local root of
                // domain
                L2error_dom = 0.0;
            }

            L2error += L2error_dom * L2error_dom;

            if (Normalised == true)
            {
                Array<OneD, NekDouble> one(m_vessels[vesselid]->GetNpoints(),
                                           1.0);

                Vol += m_vessels[vesselid]->Integral(one);
            }
        }
    }
    else
    {
        ASSERTL0(false, "Not set up");
    }

    m_comm->AllReduce(L2error, LibUtilities::ReduceSum);

    if (Normalised == true)
    {
        m_comm->AllReduce(Vol, LibUtilities::ReduceSum);

        L2error = sqrt(L2error / Vol);
    }
    else
    {
        L2error = sqrt(L2error);
    }

    return L2error;
}

/**
 * Compute the error in the L_inf-norm
 * @param   field           The field to compare.
 * @param   exactsoln       The exact solution to compare with.
 * @returns                 Error in the L_inft-norm.
 */
NekDouble PulseWaveSystem::v_LinfError(unsigned int field,
                                       const Array<OneD, NekDouble> &exact)
{
    NekDouble LinferrorDom, Linferror = -1.0;

    for (int omega = 0; omega < m_nDomains; ++omega)
    {
        int vesselid = field + omega * m_nVariables;

        // since exactsoln is passed for just the first field size we need
        // to reset it each domain
        Array<OneD, NekDouble> exactsoln(m_vessels[vesselid]->GetNpoints());
        m_fields[field] = m_vessels[omega * m_nVariables];
        EvaluateExactSolution(field, exactsoln, m_time);

        if (m_NumQuadPointsError == 0)
        {
            if (m_vessels[vesselid]->GetPhysState() == false)
            {
                m_vessels[vesselid]->BwdTrans(
                    m_vessels[vesselid]->GetCoeffs(),
                    m_vessels[vesselid]->UpdatePhys());
            }

            if (exactsoln.size())
            {
                LinferrorDom = m_vessels[vesselid]->Linf(
                    m_vessels[vesselid]->GetPhys(), exactsoln);
            }
            else if (m_session->DefinesFunction("ExactSolution"))
            {
                Array<OneD, NekDouble> exactsoln(
                    m_vessels[vesselid]->GetNpoints());

                GetFunction("ExactSolution")
                    ->Evaluate(m_session->GetVariable(field), exactsoln,
                               m_time);

                LinferrorDom = m_vessels[vesselid]->Linf(
                    m_vessels[vesselid]->GetPhys(), exactsoln);
            }
            else
            {
                LinferrorDom = 0.0;
            }

            Linferror = (Linferror > LinferrorDom) ? Linferror : LinferrorDom;
        }
        else
        {
            ASSERTL0(false, "ErrorExtraPoints not allowed for this solver");
        }
    }
    m_comm->AllReduce(Linferror, LibUtilities::ReduceMax);
    return Linferror;
}

} // namespace Nektar
