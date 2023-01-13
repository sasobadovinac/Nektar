///////////////////////////////////////////////////////////////////////////////
//
// File: TestVarcoeffHashing.cpp
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
// Description: Unit test for varcoeff hashing
//
///////////////////////////////////////////////////////////////////////////////
#include <MultiRegions/ContField.h>
#include <MultiRegions/MultiRegions.hpp>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/test/unit_test.hpp>

namespace Nektar
{
namespace VarcoeffHashingTest
{
using namespace MultiRegions;

// Forward declare global variables for all tests
boost::filesystem::path pathGlobal;
LibUtilities::SessionReaderSharedPtr lsession;
MultiRegions::ContFieldSharedPtr contfield, forcefield;
Array<OneD, Array<OneD, NekDouble>> velocity1, velocity2;

StdRegions::ConstFactorMap factors;
StdRegions::VarCoeffMap varcoeffs1, varcoeffs2;

// Define session file for the tests
// Ideally use a session-free approach in the future
static const std::string sessionfile =
    R"END(<?xml version="1.0" encoding="UTF-8"?>
<NEKTAR>
    <GEOMETRY DIM="2" SPACE="2">

        <VERTEX>
            <V ID="0"> 0.0 0.0 0.0 </V>
            <V ID="1"> 0.5 0.0 0.0 </V>
            <V ID="2"> 1.0 0.0 0.0 </V>
            <V ID="3"> 0.0 0.5 0.0 </V>
            <V ID="4"> 0.5 0.5 0.0 </V>
            <V ID="5"> 1.0 0.5 0.0 </V>
            <V ID="6"> 0.0 1.0 0.0 </V>
            <V ID="7"> 0.5 1.0 0.0 </V>
            <V ID="8"> 1.0 1.0 0.0 </V>
        </VERTEX>

        <EDGE>
            <E ID="0"> 0 1 </E>
            <E ID="1"> 1 2 </E>
            <E ID="2"> 0 3 </E>
            <E ID="3"> 1 4 </E>
            <E ID="4"> 2 5 </E>
            <E ID="5"> 3 4 </E>
            <E ID="6"> 4 5 </E>
            <E ID="7"> 3 6 </E>
            <E ID="8"> 4 7 </E>
            <E ID="9"> 5 8 </E>
            <E ID="10"> 6 7 </E>
            <E ID="11"> 7 8 </E>
        </EDGE>

        <ELEMENT>
            <Q ID="0"> 0 3 5 2 </Q>
            <Q ID="1"> 1 4 6 3 </Q>
            <Q ID="2"> 5 8 10 7 </Q>
            <Q ID="3"> 6 9 11 8 </Q>
        </ELEMENT>

        <COMPOSITE>
            <C ID="0"> Q[0-3] </C>
            <C ID="1"> E[0,1] </C>          <!-- Lower wall -->
            <C ID="2"> E[2,7] </C>          <!-- Inflow -->
            <C ID="3"> E[4,9] </C>          <!-- Outflow -->
            <C ID="4"> E[10,11] </C>        <!-- Upper wall -->
        </COMPOSITE>

        <DOMAIN> C[0] </DOMAIN>

    </GEOMETRY>
    
    <EXPANSIONS>
        <E COMPOSITE="C[0]" FIELDS="u" TYPE="MODIFIED" NUMMODES="3"/>
    </EXPANSIONS>
    
    <CONDITIONS>

        <PARAMETERS>
            <P> lambda = 1.0  </P>
        </PARAMETERS>
        <VARIABLES>
            <V ID="0"> u </V>
        </VARIABLES>

        <BOUNDARYREGIONS>
            <B ID="0"> C[1] </B>
            <B ID="1"> C[2] </B>
            <B ID="2"> C[3] </B>
            <B ID="3"> C[4] </B>
        </BOUNDARYREGIONS>

        <BOUNDARYCONDITIONS>
            <REGION REF="0"> <!-- Stationary -->
                <D VAR="u" VALUE="0" />
            </REGION>
            <REGION REF="1"> <!-- Inflow -->
                <D VAR="u" VALUE="y" />
            </REGION>
            <REGION REF="2"> <!-- Outflow -->
                <N VAR="u" VALUE="0" />
            </REGION>
            <REGION REF="3"> <!-- Moving -->
                <D VAR="u" VALUE="1" />
            </REGION>
        </BOUNDARYCONDITIONS>

        <FUNCTION NAME="d00">
            <E VAR="u" VALUE="0.1*y" />
        </FUNCTION>

        <FUNCTION NAME="d00B">
            <E VAR="u" VALUE="0.5*y" />
        </FUNCTION>

        <FUNCTION NAME="Forcing">
            <E VAR="u" VALUE="sin(y)" />
        </FUNCTION>

        <FUNCTION NAME="ExactSolution">
            <E VAR="u" VALUE="y" />
        </FUNCTION>

    </CONDITIONS>
    
</NEKTAR>
)END";

/*
 * Create the session file at working directory
 * TODO Use session-free approach to save effort on this test
 */
void createSessionFile(boost::filesystem::path &ph)
{
    // Create temporary directory
    using namespace boost::filesystem;
    ph = temp_directory_path() / unique_path();
    create_directories(ph);

    // Create file in working directory
    ph /= "TestVarcoeffHashing.xml"; // append filename
    boost::filesystem::ofstream sfile(ph);
    sfile << sessionfile;
    sfile.close();
}

/*
 * Setup for generic ContField solve
 * Read session file, setup contfield and
 * auxilliary functions for forcing and varcoeffs
 */
void setupContFieldSolve(boost::filesystem::path &ph,
                         LibUtilities::SessionReaderSharedPtr &Session,
                         MultiRegions::ContFieldSharedPtr &Exp,
                         MultiRegions::ContFieldSharedPtr &Fce,
                         Array<OneD, Array<OneD, NekDouble>> &vel1,
                         Array<OneD, Array<OneD, NekDouble>> &vel2)
{
#ifdef _WIN32
    // Prepare input file name - convert wstring to string to support Windows
    const std::wstring input_file_wstr(ph.c_str());
    const std::string input_file_str(input_file_wstr.begin(),
                                     input_file_wstr.end());
    char *input_file = (char *)input_file_str.c_str();
#else
    char *input_file = (char *)ph.c_str();
#endif

    // Setup parameters for call to IncNavierStokesSolver
    int argc     = 2;
    char *argv[] = {(char *)("IncNavierStokesSolver"), input_file, nullptr};

    // Read Session, MeshGraph and create ContField
    Session = LibUtilities::SessionReader::CreateInstance(argc, argv);
    SpatialDomains::MeshGraphSharedPtr Graph =
        SpatialDomains::MeshGraph::Read(Session);
    Exp = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
        Session, Graph, Session->GetVariable(0));
    Fce = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(*Exp);

    // Read functons for varcoeffs and Forcing
    LibUtilities::EquationSharedPtr d00func  = Session->GetFunction("d00", 0);
    LibUtilities::EquationSharedPtr d00Bfunc = Session->GetFunction("d00B", 0);
    LibUtilities::EquationSharedPtr ffunc = Session->GetFunction("Forcing", 0);

    // Get coordinates to evaluate functions
    unsigned int npoints = Exp->GetNpoints();
    unsigned int coordim = Exp->GetCoordim(0);
    Array<OneD, NekDouble> x0(npoints), x1(npoints), x2(npoints);
    Array<OneD, NekDouble> d00(npoints, 0.0), d00B(npoints, 0.0),
        fcePhys(npoints);
    vel1 = Array<OneD, Array<OneD, NekDouble>>(coordim);
    vel2 = Array<OneD, Array<OneD, NekDouble>>(coordim);

    // Define lambda
    factors[StdRegions::eFactorLambda] = Session->GetParameter("lambda");
    Exp->GetCoords(x0, x1, x2);

    // Define 1st varcoeff
    d00func->Evaluate(x0, x1, x2, d00);
    varcoeffs1[StdRegions::eVarCoeffD00] = d00;

    // Define 2nd varcoeff
    d00Bfunc->Evaluate(x0, x1, x2, d00B);
    varcoeffs2[StdRegions::eVarCoeffD00] = d00B;

    // Define Forcing
    ffunc->Evaluate(x0, x1, x2, fcePhys);
    Fce->SetPhys(fcePhys);

    // Define Advection velocities
    for (int i = 0; i < coordim; i++)
    {
        vel1[i] = Array<OneD, NekDouble>(npoints, 0.1);
        vel2[i] = Array<OneD, NekDouble>(npoints, 0.2);
    }

    // Clear solution
    Vmath::Zero(Exp->GetNcoeffs(), Exp->UpdateCoeffs(), 1);
}

BOOST_AUTO_TEST_CASE(TestVarcoeffHashing)
{
    // Create session file and setup Contfield
    createSessionFile(pathGlobal);
    setupContFieldSolve(pathGlobal, lsession, contfield, forcefield, velocity1,
                        velocity2);

    // Check no initial GlobalLinSys
    BOOST_CHECK_EQUAL(contfield->GetPoolCount("GlobalLinSys"),
                      0); // Check count == 0

    // Create GlobalLinSys(varcoeff1)
    contfield->HelmSolve(forcefield->GetPhys(), contfield->UpdateCoeffs(),
                         factors, varcoeffs1);
    BOOST_CHECK_EQUAL(contfield->GetPoolCount("GlobalLinSys"),
                      1); // Check count == 1

    // Create new GlobalLinSys(varcoeff2)
    contfield->HelmSolve(forcefield->GetPhys(), contfield->UpdateCoeffs(),
                         factors, varcoeffs2);
    BOOST_CHECK_EQUAL(contfield->GetPoolCount("GlobalLinSys"),
                      2); // Check count == 2

    // Again use 2nd GlobalLinSys(varcoeff2)
    contfield->HelmSolve(forcefield->GetPhys(), contfield->UpdateCoeffs(),
                         factors, varcoeffs2);
    BOOST_CHECK_EQUAL(contfield->GetPoolCount("GlobalLinSys"),
                      2); // Check count still == 2

    // Again use 1st GlobalLinSys(varcoeff1)
    contfield->HelmSolve(forcefield->GetPhys(), contfield->UpdateCoeffs(),
                         factors, varcoeffs1);
    BOOST_CHECK_EQUAL(contfield->GetPoolCount("GlobalLinSys"),
                      2); // Check count still == 2

    // Clear for next tests
    contfield->ClearGlobalLinSysManager();
}

BOOST_AUTO_TEST_CASE(TestUnsetGlobalLinSys)
{
    // Check that the test setup is clear
    BOOST_CHECK_EQUAL(contfield->GetPoolCount("GlobalLinSys"),
                      0); // Check count == 0

    // Create 1st GlobalLinSys
    auto gkey1 = contfield->LinearAdvectionDiffusionReactionSolve(
        velocity1, forcefield->GetPhys(), contfield->UpdateCoeffs(),
        factors[StdRegions::eFactorLambda]);
    BOOST_CHECK_EQUAL(contfield->GetPoolCount("GlobalLinSys"),
                      1); // Check count == 1

    // Create 2nd GlobalLinSys
    auto gkey2 = contfield->LinearAdvectionDiffusionReactionSolve(
        velocity2, forcefield->GetPhys(), contfield->UpdateCoeffs(),
        factors[StdRegions::eFactorLambda]);
    BOOST_CHECK_EQUAL(contfield->GetPoolCount("GlobalLinSys"),
                      2); // Check count == 2

    // Unset GlobalLinSys1
    contfield->UnsetGlobalLinSys(gkey1, true);
    BOOST_CHECK_EQUAL(contfield->GetPoolCount("GlobalLinSys"),
                      1); // Check count == 1

    // Unset GlobalLinSys2
    contfield->UnsetGlobalLinSys(gkey2, true);
    BOOST_CHECK_EQUAL(contfield->GetPoolCount("GlobalLinSys"),
                      0); // Check count == 0

    // Clear for next tests
    contfield->ClearGlobalLinSysManager();
}
} // namespace VarcoeffHashingTest
} // namespace Nektar
