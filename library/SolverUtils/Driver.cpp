///////////////////////////////////////////////////////////////////////////////
//
// File Driver.cpp
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
// Description: Base class for Drivers
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Driver.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

std::string Driver::evolutionOperatorLookupIds[6] = {
    LibUtilities::SessionReader::RegisterEnumValue(
            "EvolutionOperator","Nonlinear"      ,eNonlinear),
    LibUtilities::SessionReader::RegisterEnumValue(
            "EvolutionOperator","Direct"         ,eDirect),
    LibUtilities::SessionReader::RegisterEnumValue(
            "EvolutionOperator","Adjoint"        ,eAdjoint),
    LibUtilities::SessionReader::RegisterEnumValue(
            "EvolutionOperator","TransientGrowth",eTransientGrowth),
    LibUtilities::SessionReader::RegisterEnumValue(
            "EvolutionOperator","SkewSymmetric"  ,eSkewSymmetric),
    LibUtilities::SessionReader::RegisterEnumValue(
            "EvolutionOperator","AdaptiveSFD"    ,eAdaptiveSFD)
};
std::string Driver::evolutionOperatorDef =
    LibUtilities::SessionReader::RegisterDefaultSolverInfo(
            "EvolutionOperator","Nonlinear");
std::string Driver::driverDefault =
    LibUtilities::SessionReader::RegisterDefaultSolverInfo(
            "Driver",           "Standard");

DriverFactory& GetDriverFactory()
{
    static DriverFactory instance;
    return instance;
}

/**
 *
 */
Driver::Driver(const LibUtilities::SessionReaderSharedPtr pSession,
               const SpatialDomains::MeshGraphSharedPtr pGraph)
    : m_comm(pSession->GetComm()),
      m_session(pSession),
      m_graph(pGraph)
{
}

Driver::~Driver()

{
}

/**
 *
 */
void Driver::v_InitObject(ostream &out)
{
    try
    {
        // Retrieve the equation system to solve.
        ASSERTL0(m_session->DefinesSolverInfo("EqType"),
                 "EqType SolverInfo tag must be defined.");
        std::string vEquation = m_session->GetSolverInfo("EqType");
        if (m_session->DefinesSolverInfo("SolverType"))
        {
            vEquation = m_session->GetSolverInfo("SolverType");
        }

        // Check such a module exists for this equation.
        ASSERTL0(GetEquationSystemFactory().ModuleExists(vEquation),
                 "EquationSystem '" + vEquation + "' is not defined.\n"
                 "Ensure equation name is correct and module is compiled.\n");

        // Retrieve the type of evolution operator to use
        /// @todo At the moment this is Navier-Stokes specific - generalise?
        m_EvolutionOperator =
            m_session->GetSolverInfoAsEnum<EvolutionOperatorType>(
                    "EvolutionOperator");

        m_nequ = ((m_EvolutionOperator == eTransientGrowth ||
                   m_EvolutionOperator == eAdaptiveSFD) ? 2 : 1);

        // Determine if parareal solver is used
        std::string DriverType = m_session->GetSolverInfo("Driver");
        if (DriverType == "Parareal")
        {
            m_nequ += 1; // Add coarse parareal solver
        }

        m_equ = Array<OneD, EquationSystemSharedPtr>(m_nequ);

        // Set the AdvectiveType tag and create EquationSystem objects.
        switch (m_EvolutionOperator)
        {
            case eNonlinear:
                m_session->SetTag("AdvectiveType","Convective");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);
                if (DriverType == "Parareal")
                {
                    // Set coarse parareal session file
                    setPararealSessionFile();

                    // Set coarse parareal solver
                    m_equ[1] = GetEquationSystemFactory().CreateInstance(
                        vEquation, m_session_coarse, m_graph_coarse);
                }
                break;
            case eDirect:
                m_session->SetTag("AdvectiveType","Linearised");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);
                if (DriverType == "Parareal")
                {
                    // Set coarse parareal session file
                    setPararealSessionFile();

                    // Set coarse parareal solver
                    m_equ[1] = GetEquationSystemFactory().CreateInstance(
                        vEquation, m_session_coarse, m_graph_coarse);
                }
                break;
            case eAdjoint:
                m_session->SetTag("AdvectiveType","Adjoint");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);
                if (DriverType == "Parareal")
                {
                    // Set coarse parareal session file
                    setPararealSessionFile();

                    // Set coarse parareal solver
                    m_equ[1] = GetEquationSystemFactory().CreateInstance(
                        vEquation, m_session_coarse, m_graph_coarse);
                }
                break;
            case eTransientGrowth:
                // forward timestepping
                m_session->SetTag("AdvectiveType","Linearised");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                                    vEquation, m_session, m_graph);

                // backward timestepping
                m_session->SetTag("AdvectiveType","Adjoint");
                m_equ[1] = GetEquationSystemFactory().CreateInstance(
                                    vEquation, m_session, m_graph);
                break;
            case eSkewSymmetric:
                m_session->SetTag("AdvectiveType","SkewSymmetric");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);
                if (DriverType == "Parareal")
                {
                    // Set coarse parareal session file
                    setPararealSessionFile();

                    // Set coarse parareal solver
                    m_equ[1] = GetEquationSystemFactory().CreateInstance(
                        vEquation, m_session_coarse, m_graph_coarse);
                }
                break;
            case eAdaptiveSFD:
            {
                // Coupling SFD method and Arnoldi algorithm
                // For having 2 equation systems defined into 2 different
                // session files (with the mesh into a file named 'session'.gz)
                string          meshfile;
                string          LinNSCondFile;
                vector<string>  LinNSFilename;
                meshfile = m_session->GetFilenames()[0];
                LinNSCondFile = m_session->GetSessionName();
                LinNSCondFile += "_LinNS.xml";
                LinNSFilename.push_back(meshfile);
                LinNSFilename.push_back(LinNSCondFile);

                char *argv[] = {
                    const_cast<char*>("IncNavierStokesSolver"), nullptr };
                session_LinNS = LibUtilities::SessionReader::CreateInstance(
                    1, argv, LinNSFilename, m_session->GetComm());

                SpatialDomains::MeshGraphSharedPtr graph_linns =
                    SpatialDomains::MeshGraph::Read(session_LinNS);

                //For running stability analysis
                session_LinNS->SetTag("AdvectiveType","Linearised");
                m_equ[0] = GetEquationSystemFactory().CreateInstance(
                    vEquation, session_LinNS, graph_linns);

                //For running the SFD method on the nonlinear problem
                m_session->SetTag("AdvectiveType","Convective");
                m_equ[1] = GetEquationSystemFactory().CreateInstance(
                    vEquation, m_session, m_graph);
            }
                break;
            default:
                ASSERTL0(false, "Unrecognised evolution operator.");

        }
    }
    catch (int e)
    {
        ASSERTL0(e == -1, "No such class class defined.");
        out << "An error occurred during driver initialisation." << endl;
    }
}

Array<OneD, NekDouble> Driver::v_GetRealEvl(void)
{
    ASSERTL0(false,"This routine is not valid in this class");
    return NullNekDouble1DArray;
}

Array<OneD, NekDouble> Driver::v_GetImagEvl(void)
{
    ASSERTL0(false,"This routine is not valid in this class");
    return NullNekDouble1DArray;
}

/// Set the Parareal (coarse solver) session file
void Driver::setPararealSessionFile(void)
{
    // Get the coarse solver session file
    string meshFile;
    string coarseSolverFile;
    vector<string> coarseSolverFilenames;
    meshFile         = m_session->GetFilenames()[0];
    coarseSolverFile = m_session->GetFilenames()[1];
    coarseSolverFile = coarseSolverFile.substr(0, coarseSolverFile.size() - 4);
    coarseSolverFile += "_coarseSolver.xml";
    std::ifstream f(coarseSolverFile);
    if (f.good())
    {
        coarseSolverFilenames.push_back(meshFile);
        coarseSolverFilenames.push_back(coarseSolverFile);
    }
    else
    {
        // if _coarseSolver.xml does not exit, use original session file
        coarseSolverFilenames.push_back(m_session->GetFilenames()[0]);
        coarseSolverFilenames.push_back(m_session->GetFilenames()[1]);
    }

    // Define argument for the coarse parareal solver
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

    char *argv[] = {const_cast<char *>("Solver"), // this is just a place holder
                    const_cast<char *>("--npx"),
                    const_cast<char *>(std::to_string(npx).c_str()),
                    const_cast<char *>("--npy"),
                    const_cast<char *>(std::to_string(npy).c_str()),
                    const_cast<char *>("--npz"),
                    const_cast<char *>(std::to_string(npz).c_str()),
                    const_cast<char *>("--nsz"),
                    const_cast<char *>(std::to_string(nsz).c_str()),
                    const_cast<char *>("--npt"),
                    const_cast<char *>(std::to_string(npt).c_str()),
                    nullptr};

    // Set session for coarse solver
    m_session_coarse = LibUtilities::SessionReader::CreateInstance(
        11, argv, coarseSolverFilenames, m_session->GetComm());

    // Set graph for coarse solver
    m_graph_coarse = SpatialDomains::MeshGraph::Read(m_session_coarse);

    // If a coarse solver session file is not specified, use
    // m_coarseSolveFactor to determine the timestep of the coarse solver
    if (!f.good())
    {
        double TimeStep =
            m_session->GetParameter("TimeStep") * m_coarseSolveFactor;
        int NumSteps =
            m_session->GetParameter("NumSteps") / m_coarseSolveFactor;
        m_session_coarse->SetParameter("TimeStep", TimeStep);
        m_session_coarse->SetParameter("NumSteps", NumSteps);
    }
}

}
}
