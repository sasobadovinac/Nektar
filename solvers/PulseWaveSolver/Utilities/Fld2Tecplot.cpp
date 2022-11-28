///////////////////////////////////////////////////////////////////////////////
//
// File: Fld2Tecplot.cpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SolverUtils/Driver.h>

#include <PulseWaveSolver/EquationSystems/PulseWaveSystem.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::SolverUtils;

static std::string SetToOneD =
    LibUtilities::SessionReader::RegisterCmdLineArgument(
        "SetToOneSpaceDimension", "1", "Redefine mesh to be aligned to x-axis");

int main(int argc, char *argv[])
{
    if ((argc < 3) || (argc > 4))
    {
        fprintf(stderr, "Usage: ./Fld2Tecplot [-c] file.xml file.fld\n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr session;
    SpatialDomains::MeshGraphSharedPtr graph;
    string vDriverModule;
    DriverSharedPtr drv;

    try
    {

        // Define new input with extra argument to intialisae -OneD=false
        int newargc    = argc + 1;
        char **newargv = new char *[newargc];

        newargv[0] = argv[0];
        newargv[1] = new char[31];
        strcpy(newargv[1], "--SetToOneSpaceDimension=false");

        for (int i = 1; i < argc; ++i)
        {
            newargv[i + 1] = argv[i];
        }

        // Create session reader and MeshGraph.
        session = LibUtilities::SessionReader::CreateInstance(newargc, newargv);
        graph   = SpatialDomains::MeshGraph::Read(session);
        delete[] newargv;

        // Create driver
        session->LoadSolverInfo("Driver", vDriverModule, "Standard");
        drv = GetDriverFactory().CreateInstance(vDriverModule, session, graph);

        EquationSystemSharedPtr EqSys = drv->GetEqu()[0];

        PulseWaveSystemSharedPtr PulseWave;
        if (!(PulseWave = std::dynamic_pointer_cast<PulseWaveSystem>(EqSys)))
        {
            ASSERTL0(false,
                     "Failed to dynamically cast to PulseWaveSystemOutput");
        }

        std::string fname(argv[argc - 1]);
        Array<OneD, MultiRegions::ExpListSharedPtr> Vessels;

        int ndomains = PulseWave->GetNdomains();

        PulseWave->ImportFldToMultiDomains(
            fname, Vessels = PulseWave->UpdateVessels(), ndomains);
        int fdot = fname.find_last_of('.');

        if (fdot != std::string::npos)
        {
            string ending = fname.substr(fdot);

            // If .chk or .fld we exchange the extension in the output file.
            // For all other files (e.g. .bse) we append the extension to avoid
            // conflicts.
            if (ending == ".chk" || ending == ".fld")
            {
                fname = fname.substr(0, fdot);
            }
        }

        fname += ".dat";

        ofstream outfile(fname.c_str());
        int nvariables  = session->GetVariables().size();
        std::string var = "";
        int j;
        for (j = 0; j < nvariables - 1; ++j)
        {
            var += session->GetVariable(j) + ", ";
        }
        var += session->GetVariable(j);

        Vessels[0]->WriteTecplotHeader(outfile, var);

        for (int n = 0; n < ndomains; ++n)
        {
            Vessels[n * nvariables]->WriteTecplotZone(outfile);
            for (int j = 0; j < nvariables; ++j)
            {
                Vessels[n * nvariables + j]->WriteTecplotField(outfile);
            }

            Vessels[n * nvariables]->WriteTecplotConnectivity(outfile);
        }
    }

    catch (const std::runtime_error &)
    {
        return 1;
    }
    catch (const std::string &eStr)
    {
        cout << "Error: " << eStr << endl;
    }

    return 0;
}
