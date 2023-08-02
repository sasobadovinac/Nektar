///////////////////////////////////////////////////////////////////////////////
//
// File DriverParallelInTime.cpp
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
// Description: Driver class for the parallel-in-time solver
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/DriverParallelInTime.h>
#include <boost/format.hpp>

namespace Nektar
{
namespace SolverUtils
{

/**
 *
 */
DriverParallelInTime::DriverParallelInTime(
    const LibUtilities::SessionReaderSharedPtr pSession,
    const SpatialDomains::MeshGraphSharedPtr pGraph)
    : Driver(pSession, pGraph)
{
}

/**
 *
 */
DriverParallelInTime::~DriverParallelInTime()
{
}

/**
 *
 */
void DriverParallelInTime::v_InitObject(std::ostream &out)
{
    try
    {
        // Retrieve the type of evolution operator to use
        m_EvolutionOperator =
            m_session->GetSolverInfoAsEnum<EvolutionOperatorType>(
                "EvolutionOperator");

        m_nequ = 2; // Only two time levels currently implemented.

        m_equ = Array<OneD, EquationSystemSharedPtr>(m_nequ);

        // Set the AdvectiveType tag and create EquationSystem objects.
        switch (m_EvolutionOperator)
        {
            case eNonlinear:
                SetParallelInTimeEquationSystem("Convective");
                break;
            case eDirect:
                SetParallelInTimeEquationSystem("Linearised");
                break;
            case eAdjoint:
                SetParallelInTimeEquationSystem("Adjoint");
                break;
            case eSkewSymmetric:
                SetParallelInTimeEquationSystem("SkewSymmetric");
                break;
            default:
                ASSERTL0(false, "Unrecognised evolution operator.");
        }

        // Set pointers.
        m_fineEqSys =
            std::dynamic_pointer_cast<SolverUtils::UnsteadySystem>(m_equ[0]);
        m_coarseEqSys =
            std::dynamic_pointer_cast<SolverUtils::UnsteadySystem>(m_equ[1]);
    }
    catch (int e)
    {
        ASSERTL0(e == -1, "No such class class defined.");
        out << "An error occurred during driver initialisation." << std::endl;
    }
}

/**
 *
 */
void DriverParallelInTime::v_Execute(std::ostream &out)
{
    boost::ignore_unused(out);

    ASSERTL0(false, "Specific version of Parallel-in-Time not implemented");
}

/**
 *
 */
NekDouble DriverParallelInTime::v_EstimateCommunicationTime(void)
{
    ASSERTL0(false, "Specific version of Parallel-in-Time not implemented");
    return 0.0;
}

/**
 *
 */
NekDouble DriverParallelInTime::v_EstimateRestrictionTime(void)
{
    ASSERTL0(false, "Specific version of Parallel-in-Time not implemented");
    return 0.0;
}

/**
 *
 */
NekDouble DriverParallelInTime::v_EstimateInterpolationTime(void)
{
    ASSERTL0(false, "Specific version of Parallel-in-Time not implemented");
    return 0.0;
}

/**
 *
 */
NekDouble DriverParallelInTime::v_EstimateCoarseSolverTime(void)
{
    ASSERTL0(false, "Specific version of Parallel-in-Time not implemented");
    return 0.0;
}

/**
 *
 */
NekDouble DriverParallelInTime::v_EstimateFineSolverTime(void)
{
    ASSERTL0(false, "Specific version of Parallel-in-Time not implemented");
    return 0.0;
}

/**
 *
 */
NekDouble DriverParallelInTime::v_EstimatePredictorTime(void)
{
    ASSERTL0(false, "Specific version of Parallel-in-Time not implemented");
    return 0.0;
}

/**
 *
 */
NekDouble DriverParallelInTime::v_EstimateOverheadTime(void)
{
    ASSERTL0(false, "Specific version of Parallel-in-Time not implemented");
    return 0.0;
}

/**
 *
 */
NekDouble DriverParallelInTime::v_ComputeSpeedUp(
    const size_t iter, NekDouble fineSolveTime, NekDouble coarseSolveTime,
    NekDouble restTime, NekDouble interTime, NekDouble commTime,
    NekDouble predictorOverheadTime, NekDouble overheadTime)
{
    boost::ignore_unused(iter, fineSolveTime, coarseSolveTime, restTime,
                         interTime, commTime, predictorOverheadTime,
                         overheadTime);
    ASSERTL0(false, "Specific version of Parallel-in-Time not implemented");
    return 0.0;
}

/**
 * Set the ParallelInTime (coarse solver) session file
 */
void DriverParallelInTime::SetParallelInTimeEquationSystem(
    std::string AdvectiveType)
{
    // Retrieve the equation system to solve.
    ASSERTL0(m_session->DefinesSolverInfo("EqType"),
             "EqType SolverInfo tag must be defined.");
    std::string vEquation = m_session->DefinesSolverInfo("SolverType")
                                ? m_session->GetSolverInfo("SolverType")
                                : m_session->GetSolverInfo("EqType");

    // Check such a module exists for this equation.
    ASSERTL0(GetEquationSystemFactory().ModuleExists(vEquation),
             "EquationSystem '" + vEquation +
                 "' is not defined.\n"
                 "Ensure equation name is correct and module is compiled.\n");

    // Set fine parallel-in-time solver.
    m_session->SetTag("AdvectiveType", AdvectiveType);
    m_session->SetTag("ParallelInTimeSolver", "TimeLevel0");
    m_equ[0] = GetEquationSystemFactory().CreateInstance(vEquation, m_session,
                                                         m_graph);

    // Define argument for the coarse parallel-in-time solver.
    int npx = m_session->DefinesCmdLineArgument("npx")
                  ? m_session->GetCmdLineArgument<int>("npx")
                  : 1;
    int npy = m_session->DefinesCmdLineArgument("npy")
                  ? m_session->GetCmdLineArgument<int>("npy")
                  : 1;
    int npz = m_session->DefinesCmdLineArgument("npz")
                  ? m_session->GetCmdLineArgument<int>("npz")
                  : 1;
    int nsz = m_session->DefinesCmdLineArgument("nsz")
                  ? m_session->GetCmdLineArgument<int>("nsz")
                  : 1;
    int npt = m_session->DefinesCmdLineArgument("npt")
                  ? m_session->GetCmdLineArgument<int>("npt")
                  : 1;

    // Convert into string.
    std::string npx_string = std::to_string(npx);
    std::string npy_string = std::to_string(npy);
    std::string npz_string = std::to_string(npz);
    std::string nsz_string = std::to_string(nsz);
    std::string npt_string = std::to_string(npt);
    std::string optfile    = m_session->GetSessionName() + ".opt";

    char *argv[] = {const_cast<char *>("Solver"), // this is just a place holder
                    const_cast<char *>("--npx"),
                    const_cast<char *>(npx_string.c_str()),
                    const_cast<char *>("--npy"),
                    const_cast<char *>(npy_string.c_str()),
                    const_cast<char *>("--npz"),
                    const_cast<char *>(npz_string.c_str()),
                    const_cast<char *>("--nsz"),
                    const_cast<char *>(nsz_string.c_str()),
                    const_cast<char *>("--npt"),
                    const_cast<char *>(npt_string.c_str()),
                    const_cast<char *>("--useoptfile"),
                    const_cast<char *>(optfile.c_str()),
                    nullptr};

    int argc = m_session->DefinesCmdLineArgument("useoptfile") ? 11 : 13;

    // Set session for coarse solver.
    std::vector<std::string> sessionFileNames(m_session->GetFilenames());
    for (size_t timeLevel = 1; timeLevel < m_nequ; timeLevel++)
    {
        auto session = LibUtilities::SessionReader::CreateInstance(
            argc, argv, sessionFileNames, m_session->GetComm(), timeLevel);

        // Set graph for coarse solver.
        auto graph = SpatialDomains::MeshGraph::Read(session);

        // Set BndRegionOrdering (necessary for DG with periodic BC) FIXME
        graph->SetBndRegionOrdering(m_graph->GetBndRegionOrdering());

        // Set CompositeOrdering (necessary for DG with periodic BC) FIXME
        graph->SetCompositeOrdering(m_graph->GetCompositeOrdering());

        // Retrieve the equation system to solve.
        ASSERTL0(session->DefinesSolverInfo("EqType"),
                 "EqType SolverInfo tag must be defined.");
        auto vEquation = session->DefinesSolverInfo("SolverType")
                             ? session->GetSolverInfo("SolverType")
                             : session->GetSolverInfo("EqType");

        // Check such a module exists for this equation.
        ASSERTL0(
            GetEquationSystemFactory().ModuleExists(vEquation),
            "EquationSystem '" + vEquation +
                "' is not defined.\n"
                "Ensure equation name is correct and module is compiled.\n");

        // Set coarse parallel-in-time solver.
        session->SetTag("AdvectiveType", AdvectiveType);
        session->SetTag("ParallelInTimeSolver",
                        "TimeLevel" + std::to_string(timeLevel));
        m_equ[timeLevel] = GetEquationSystemFactory().CreateInstance(
            vEquation, session, graph);
    }
}

/**
 *
 */
void DriverParallelInTime::GetParametersFromSession(void)
{
    // Parallel-in-Time iteration parameters.
    m_tolerPIT      = m_session->GetParameter("PITToler");
    m_iterMaxPIT    = m_session->GetParameter("PITIterMax");
    m_numWindowsPIT = m_session->DefinesParameter("NumWindows")
                          ? m_session->GetParameter("NumWindows")
                          : m_numWindowsPIT;

    // Time stepping parameters.
    m_fineTimeStep   = m_equ[0]->GetTimeStep();
    m_coarseTimeStep = m_equ[1]->GetTimeStep();
    m_fineSteps      = m_equ[0]->GetSteps();
    m_coarseSteps    = m_equ[1]->GetSteps();

    // I/O parameters.
    m_infoSteps  = m_fineEqSys->GetInfoSteps();
    m_checkSteps = m_fineEqSys->GetCheckpointSteps();

    // Other parameters.
    m_exactSolution = m_session->DefinesParameter("ExactSolution")
                          ? m_session->GetParameter("ExactSolution")
                          : m_exactSolution;
}

/**
 *
 */
void DriverParallelInTime::InitialiseEqSystem(bool turnoff_output)
{
    // Initialize fine solver.
    if (turnoff_output)
    {
        m_fineEqSys->SetInfoSteps(0);
        m_fineEqSys->SetCheckpointSteps(0);
    }
    m_fineEqSys->DoInitialise(true);

    // Initialize coarse solver.
    m_coarseEqSys->SetInfoSteps(0);
    m_coarseEqSys->SetCheckpointSteps(0);
    m_coarseEqSys->DoInitialise(false);
}

/**
 *
 */
void DriverParallelInTime::AllocateMemory(void)
{
    // Set some member variables.
    m_nVar       = m_fineEqSys->GetNvariables();
    m_fineNpts   = m_fineEqSys->GetNpoints();
    m_coarseNpts = m_coarseEqSys->GetNpoints();

    // Allocate memory.
    m_exactsoln = Array<OneD, Array<OneD, NekDouble>>(m_nVar);
    m_tmpfine   = Array<OneD, Array<OneD, NekDouble>>(m_nVar);
    m_tmpcoarse = Array<OneD, Array<OneD, NekDouble>>(m_nVar);
    for (size_t i = 0; i < m_nVar; ++i)
    {
        m_exactsoln[i] = Array<OneD, NekDouble>(m_fineNpts, 0.0);
        m_tmpfine[i]   = Array<OneD, NekDouble>(m_fineNpts, 0.0);
        m_tmpcoarse[i] = Array<OneD, NekDouble>(m_coarseNpts, 0.0);
    }
    m_vL2Errors   = Array<OneD, NekDouble>(m_nVar, 0.0);
    m_vLinfErrors = Array<OneD, NekDouble>(m_nVar, 0.0);
}

/**
 *
 */
void DriverParallelInTime::InitialiseInterpolationField(void)
{
    m_fineFields   = Array<OneD, MultiRegions::ExpListSharedPtr>(m_nVar);
    m_coarseFields = Array<OneD, MultiRegions::ExpListSharedPtr>(m_nVar);
    for (size_t i = 0; i < m_nVar; ++i)
    {
        m_fineFields[i] =
            MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(
                m_fineEqSys->UpdateFields()[i], true, false);
        m_coarseFields[i] =
            MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(
                m_coarseEqSys->UpdateFields()[i], true, false);
    }
}

/**
 *
 */
void DriverParallelInTime::PrintCoarseSolverInfo(std::ostream &out)
{
    if (m_chunkRank == 0 && m_comm->GetSpaceComm()->GetRank() == 0)
    {
        std::cout << "========================================================="
                     "=============="
                  << std::endl
                  << std::flush;
        std::cout << "======================== COARSE PROPAGATOR INFO "
                     "======================="
                  << std::endl
                  << std::flush;

        m_coarseEqSys->PrintSummary(out);

        std::cout << std::endl << std::flush;
    }
}

/**
 *
 */
void DriverParallelInTime::PrintFineSolverInfo(std::ostream &out)
{
    if (m_chunkRank == 0 && m_comm->GetSpaceComm()->GetRank() == 0)
    {
        std::cout << "========================================================="
                     "=============="
                  << std::endl
                  << std::flush;
        std::cout << "========================= FINE PROPAGATOR INFO "
                     "========================"
                  << std::endl
                  << std::flush;

        m_fineEqSys->PrintSummary(out);

        std::cout << std::endl << std::flush;
    }
}

/**
 *
 */
void DriverParallelInTime::PrintHeaderTitle1(const std::string &title)
{
    if (m_chunkRank == m_numChunks - 1 &&
        m_comm->GetSpaceComm()->GetRank() == 0)
    {
        std::cout << std::endl;
        std::cout << "*******************************************" << std::endl
                  << std::flush;
        std::cout << title << std::endl << std::flush;
        std::cout << "*******************************************" << std::endl
                  << std::flush;
    }
}

/**
 *
 */
void DriverParallelInTime::PrintHeaderTitle2(const std::string &title)
{
    if (m_chunkRank == m_numChunks - 1 &&
        m_comm->GetSpaceComm()->GetRank() == 0)
    {
        std::cout << "-------------------------------------------" << std::endl
                  << std::flush;
        std::cout << title << std::endl << std::flush;
        std::cout << "-------------------------------------------" << std::endl
                  << std::flush;
    }
}

/**
 *
 */
void DriverParallelInTime::PrintComputationalTime(const NekDouble time)
{
    if (m_chunkRank == m_numChunks - 1 &&
        m_comm->GetSpaceComm()->GetRank() == 0)
    {
        std::cout << "Total Computation Time : " << time << "s" << std::endl
                  << std::flush;
    }
}

/**
 *
 */
void DriverParallelInTime::RecvInitialConditionFromPreviousProc(
    Array<OneD, Array<OneD, NekDouble>> &array, int &convergence)
{
    LibUtilities::CommSharedPtr tComm = m_session->GetComm()->GetTimeComm();

    if (!convergence)
    {
        if (m_chunkRank > 0)
        {
            tComm->Recv(m_chunkRank - 1, convergence);
            RecvInitialConditionFromPreviousProc(array);
        }
    }
}

/**
 *
 */
void DriverParallelInTime::RecvInitialConditionFromPreviousProc(
    Array<OneD, Array<OneD, NekDouble>> &array)
{
    LibUtilities::CommSharedPtr tComm = m_session->GetComm()->GetTimeComm();

    if (m_chunkRank > 0)
    {
        for (size_t i = 0; i < array.size(); ++i)
        {
            tComm->Recv(m_chunkRank - 1, array[i]);
        }
    }
}

/**
 *
 */
void DriverParallelInTime::SendSolutionToNextProc(
    Array<OneD, Array<OneD, NekDouble>> &array, int &convergence)
{
    LibUtilities::CommSharedPtr tComm = m_session->GetComm()->GetTimeComm();

    if (m_chunkRank < m_numChunks - 1)
    {
        tComm->Send(m_chunkRank + 1, convergence);
    }
    SendSolutionToNextProc(array);
}

/**
 *
 */
void DriverParallelInTime::SendSolutionToNextProc(
    Array<OneD, Array<OneD, NekDouble>> &array)
{
    LibUtilities::CommSharedPtr tComm = m_session->GetComm()->GetTimeComm();

    if (m_chunkRank < m_numChunks - 1)
    {
        for (size_t i = 0; i < array.size(); ++i)
        {
            tComm->Send(m_chunkRank + 1, array[i]);
        }
    }
}

/**
 *
 */
void DriverParallelInTime::CopySolutionVector(
    const Array<OneD, const Array<OneD, NekDouble>> &in,
    Array<OneD, Array<OneD, NekDouble>> &out)
{
    for (size_t i = 0; i < in.size(); ++i)
    {
        Vmath::Vcopy(in[i].size(), in[i], 1, out[i], 1);
    }
}

/**
 *
 */
void DriverParallelInTime::CopyFromFinePhysField(
    Array<OneD, Array<OneD, NekDouble>> &out)
{
    for (size_t i = 0; i < out.size(); ++i)
    {
        m_fineEqSys->CopyFromPhysField(i, out[i]);
    }
}

/**
 *
 */
void DriverParallelInTime::CopyFromCoarsePhysField(
    Array<OneD, Array<OneD, NekDouble>> &out)
{
    for (size_t i = 0; i < out.size(); ++i)
    {
        m_coarseEqSys->CopyFromPhysField(i, out[i]);
    }
}

/**
 *
 */
void DriverParallelInTime::CopyToFinePhysField(
    const Array<OneD, const Array<OneD, NekDouble>> &in)
{
    for (size_t i = 0; i < in.size(); ++i)
    {
        m_fineEqSys->CopyToPhysField(i, in[i]);
    }
}

/**
 *
 */
void DriverParallelInTime::CopyToCoarsePhysField(
    const Array<OneD, const Array<OneD, NekDouble>> &in)
{
    for (size_t i = 0; i < in.size(); ++i)
    {
        m_coarseEqSys->CopyToPhysField(i, in[i]);
    }
}

/**
 *
 */
void DriverParallelInTime::UpdateSolution(
    const Array<OneD, const Array<OneD, NekDouble>> &in)
{
    for (size_t i = 0; i < in.size(); ++i)
    {
        m_fineEqSys->CopyToPhysField(i, in[i]);
        m_fineEqSys->UpdateFields()[i]->FwdTrans(
            m_fineEqSys->UpdateFields()[i]->GetPhys(),
            m_fineEqSys->UpdateFields()[i]->UpdateCoeffs());
    }
}

/**
 *
 */
void DriverParallelInTime::EvaluateExactSolution(const NekDouble &time)
{
    for (size_t i = 0; i < m_exactsoln.size(); ++i)
    {
        m_fineEqSys->EvaluateExactSolution(i, m_exactsoln[i], time);
    }
}

/**
 *
 */
void DriverParallelInTime::SolutionConvergenceMonitoring(
    const NekDouble &CPUtime)
{
    UpdateErrorNorm(true);
    PrintErrorNorm(true);
    PrintComputationalTime(CPUtime);
}

/**
 *
 */
void DriverParallelInTime::SolutionConvergenceSummary(const NekDouble &CPUtime)
{
    UpdateErrorNorm(false);
    PrintErrorNorm(false);
    PrintComputationalTime(CPUtime);
}

/**
 *
 */
void DriverParallelInTime::UpdateErrorNorm(const bool normalized)
{
    for (size_t i = 0; i < m_vL2Errors.size(); ++i)
    {
        m_vL2Errors[i]   = m_fineEqSys->L2Error(i, m_exactsoln[i], normalized);
        m_vLinfErrors[i] = m_fineEqSys->LinfError(i, m_exactsoln[i]);
    }
}

/**
 *
 */
void DriverParallelInTime::PrintErrorNorm(const bool normalized)
{
    if (m_chunkRank == m_numChunks - 1 &&
        m_comm->GetSpaceComm()->GetRank() == 0)
    {
        for (size_t i = 0; i < m_vL2Errors.size(); ++i)
        {
            if (normalized)
            {
                std::cout << "L2 error (variable "
                          << m_fineEqSys->GetVariable(i)
                          << ") : " << m_vL2Errors[i] << std::endl
                          << std::flush;
                std::cout << "Linf error (variable "
                          << m_fineEqSys->GetVariable(i)
                          << ") : " << m_vLinfErrors[i] << std::endl
                          << std::flush;
            }
            else
            {
                std::cout << "L 2 error (variable "
                          << m_fineEqSys->GetVariable(i)
                          << ") : " << m_vL2Errors[i] << std::endl
                          << std::flush;
                std::cout << "L inf error (variable "
                          << m_fineEqSys->GetVariable(i)
                          << ") : " << m_vLinfErrors[i] << std::endl
                          << std::flush;
            }
        }
    }
}

/**
 *
 */
NekDouble DriverParallelInTime::vL2ErrorMax(void)
{
    NekDouble L2Error = 0.0;
    for (size_t i = 0; i < m_vL2Errors.size(); ++i)
    {
        L2Error = std::max(L2Error, m_vL2Errors[i]);
    }
    return L2Error;
}

/**
 *
 */
void DriverParallelInTime::Interpolator(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    size_t nVar = inarray.size();

    // Interpolate from fine to coarse
    if (inarray[0].size() > outarray[0].size())
    {
        for (size_t i = 0; i < nVar; ++i)
        {
            m_fineFields[i]->UpdatePhys()   = inarray[i];
            m_coarseFields[i]->UpdatePhys() = outarray[i];
        }
        InterpExp1ToExp2(m_fineFields, m_coarseFields);
    }
    // Interpolate from coarse to fine
    else
    {
        for (size_t i = 0; i < nVar; ++i)
        {
            m_coarseFields[i]->UpdatePhys() = inarray[i];
            m_fineFields[i]->UpdatePhys()   = outarray[i];
        }
        InterpExp1ToExp2(m_coarseFields, m_fineFields);
    }
}

/**
 *
 */
void DriverParallelInTime::SpeedUpAnalysis()
{
    // Print header.
    PrintHeaderTitle1("PARAREAL SPEED-UP ANALYSIS");

    // Mean communication time.
    NekDouble commTime = v_EstimateCommunicationTime();
    PrintHeaderTitle2("Mean Communication Time = " +
                      (boost::format("%1$.6e") % commTime).str() + "s");

    // Mean restriction time.
    NekDouble restTime = v_EstimateRestrictionTime();
    PrintHeaderTitle2("Mean Restriction Time = " +
                      (boost::format("%1$.6f") % restTime).str() + "s");

    // Mean interpolation time.
    NekDouble interTime = v_EstimateInterpolationTime();
    PrintHeaderTitle2("Mean Interpolation Time = " +
                      (boost::format("%1$.6f") % interTime).str() + "s");

    // Mean coarse solver time.
    NekDouble coarseSolveTime = v_EstimateCoarseSolverTime();
    PrintHeaderTitle2("Mean Coarse Solve Time = " +
                      (boost::format("%1$.6f") % coarseSolveTime).str() + "s");

    // Mean fine solver time.
    NekDouble fineSolveTime = v_EstimateFineSolverTime();
    PrintHeaderTitle2("Mean Fine Solve Time = " +
                      (boost::format("%1$.6f") % fineSolveTime).str() + "s");

    // Mean predictor time.
    NekDouble predictorTime = v_EstimatePredictorTime();
    PrintHeaderTitle2("Mean Predictor Time = " +
                      (boost::format("%1$.6f") % predictorTime).str() + "s");

    // Mean overhead time.
    NekDouble OverheadTime = v_EstimateOverheadTime();
    PrintHeaderTitle2("Mean Overhead Time = " +
                      (boost::format("%1$.6f") % OverheadTime).str() + "s");

    // Print speedup time.
    PrintSpeedUp(fineSolveTime, coarseSolveTime, restTime, interTime, commTime,
                 predictorTime, OverheadTime);
}

/**
 *
 */
void DriverParallelInTime::PrintSpeedUp(NekDouble fineSolveTime,
                                        NekDouble coarseSolveTime,
                                        NekDouble restTime, NekDouble interTime,
                                        NekDouble commTime,
                                        NekDouble predictOverheadTime,
                                        NekDouble OverheadTime)
{
    if (m_chunkRank == m_numChunks - 1 &&
        m_comm->GetSpaceComm()->GetRank() == 0)
    {
        // Print maximum theoretical speed-up
        PrintHeaderTitle2("Maximum Speed-up");
        for (size_t k = 1; k <= m_numChunks; k++)
        {
            NekDouble speedup =
                v_ComputeSpeedUp(k, fineSolveTime, coarseSolveTime, 0.0, 0.0,
                                 0.0, predictOverheadTime, OverheadTime);
            std::cout << "Speed-up (" << k << ") = " << speedup << std::endl
                      << std::flush;
        }

        // Print speed-up with interpolation
        PrintHeaderTitle2("Speed-up with interp.");
        for (size_t k = 1; k <= m_numChunks; k++)
        {
            NekDouble speedup = v_ComputeSpeedUp(
                k, fineSolveTime, coarseSolveTime, restTime, interTime, 0.0,
                predictOverheadTime, OverheadTime);
            std::cout << "Speed-up (" << k << ") = " << speedup << std::endl
                      << std::flush;
        }

        // Print speed-up with interpolation and communication
        PrintHeaderTitle2("Speed-up with comm. and interp.");
        for (size_t k = 1; k <= m_numChunks; k++)
        {
            NekDouble speedup = v_ComputeSpeedUp(
                k, fineSolveTime, coarseSolveTime, restTime, interTime,
                commTime, predictOverheadTime, OverheadTime);
            std::cout << "Speed-up (" << k << ") = " << speedup << std::endl
                      << std::flush;
        }
        std::cout << "-------------------------------------------" << std::endl
                  << std::flush;
    }
}

/**
 *
 */
NekDouble DriverParallelInTime::EstimateCommunicationTime(
    Array<OneD, Array<OneD, NekDouble>> &buffer1,
    Array<OneD, Array<OneD, NekDouble>> &buffer2)
{
    if (m_numChunks == 1)
    {
        return 0.0;
    }
    else
    {
        // Average communication time over niter iteration.
        size_t niter = 20;
        Nektar::LibUtilities::Timer timer;
        LibUtilities::CommSharedPtr tComm = m_session->GetComm()->GetTimeComm();
        for (size_t n = 0; n <= niter; n++)
        {
            if (n == 1)
            {
                timer.Start(); // Ignore the first iteration
            }

            if (m_chunkRank == 0)
            {
                for (size_t i = 0; i < buffer1.size(); ++i)
                {
                    tComm->Send(m_numChunks - 1, buffer1[i]);
                }
            }

            if (m_chunkRank == m_numChunks - 1)
            {
                for (size_t i = 0; i < buffer2.size(); ++i)
                {
                    tComm->Recv(0, buffer2[i]);
                }
            }
        }
        timer.Stop();
        return timer.Elapsed().count() / niter;
    }
}

/**
 *
 */
void InterpExp1ToExp2(const Array<OneD, MultiRegions::ExpListSharedPtr> exp1,
                      Array<OneD, MultiRegions::ExpListSharedPtr> &exp2)
{

    // Interpolation from exp1 -> exp2 assuming that exp1 and exp2 are the
    // same explists, but at potentially different polynomial orders.
    if (exp1.size() != exp2.size())
    {
        NEKERROR(ErrorUtil::efatal, "not the same mesh")
    }

    for (int n = 0; n < exp1.size(); ++n)
    {
        // Interpolation from exp1 -> exp2 assuming that exp1 and exp2 are the
        // same explists, but at potentially different polynomial orders.
        if (exp1[n]->GetExpSize() != exp2[n]->GetExpSize())
        {
            NEKERROR(ErrorUtil::efatal, "not the same mesh")
        }

        // If same polynomial orders, simply copy solution
        if (exp1[n]->GetTotPoints() == exp2[n]->GetTotPoints())
        {
            Vmath::Vcopy(exp1[n]->GetTotPoints(), exp1[n]->GetPhys(), 1,
                         exp2[n]->UpdatePhys(), 1);
        }
        // If different polynomial orders, interpolate solution
        else
        {
            // Transform solution from physical to coefficient space
            exp1[n]->FwdTrans(exp1[n]->GetPhys(), exp1[n]->UpdateCoeffs());

            for (int i = 0; i < exp1[n]->GetExpSize(); ++i)
            {
                // Get the elements
                LocalRegions::ExpansionSharedPtr elmt1 = exp1[n]->GetExp(i),
                                                 elmt2 = exp2[n]->GetExp(i);

                // Get the offset of elements in the storage arrays.
                int offset1 = exp1[n]->GetCoeff_Offset(i);
                int offset2 = exp2[n]->GetCoeff_Offset(i);

                // Get number of modes
                Array<OneD, LibUtilities::BasisSharedPtr> base1 =
                    elmt1->GetBase();
                std::vector<unsigned int> nummodes1(base1.size());
                std::vector<LibUtilities::BasisType> btype1(base1.size());
                for (int j = 0; j < nummodes1.size(); ++j)
                {
                    nummodes1[j] = base1[j]->GetNumModes();
                    btype1[j]    = base1[j]->GetBasisType();
                }

                // Extract data from exp1 -> exp2.
                elmt2->ExtractDataToCoeffs(
                    &exp1[n]->GetCoeffs()[offset1], nummodes1, 0,
                    &exp2[n]->UpdateCoeffs()[offset2], btype1);
            }

            // Transform solution back to physical space
            exp2[n]->BwdTrans(exp2[n]->GetCoeffs(), exp2[n]->UpdatePhys());
        }
    }
}

} // namespace SolverUtils
} // namespace Nektar
