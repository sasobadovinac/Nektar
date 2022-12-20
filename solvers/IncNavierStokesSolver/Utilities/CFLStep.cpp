///////////////////////////////////////////////////////////////////////////////
//
// File: CFLStep.cpp
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

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::SolverUtils;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        fprintf(stderr, "Usage: ./CflStep file.xml \n");
        fprintf(stderr, "\t Method will read intiial conditions section of "
                        ".xml file for input \n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr session;
    SpatialDomains::MeshGraphSharedPtr graph;

    string vDriverModule;
    DriverSharedPtr drv;
    try
    {
        // Create session reader.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv);

        // Create MeshGraph.
        graph = SpatialDomains::MeshGraph::Read(session);

        // Create driver
        session->LoadSolverInfo("Driver", vDriverModule, "Standard");
        drv = GetDriverFactory().CreateInstance(vDriverModule, session, graph);

        EquationSystemSharedPtr EqSys   = drv->GetEqu()[0];
        IncNavierStokesSharedPtr IncNav = EqSys->as<IncNavierStokes>();

        IncNav->SetInitialConditions(0.0, false);
        Array<OneD, NekDouble> cfl = IncNav->GetElmtCFLVals();

        // Reset Pressure field with CFL values
        Array<OneD, MultiRegions::ExpListSharedPtr> fields =
            IncNav->UpdateFields();
        int i, n, nquad, cnt;
        int nfields = fields.size();
        int nexp    = fields[0]->GetExpSize();

        int elmtid = Vmath::Imax(nexp, cfl, 1);

        cout << "Max CFL: " << cfl[elmtid] << " In element " << elmtid << endl;

        for (n = 0; n < nfields; ++n)
        {
            if (session->GetVariable(n) == "p")
            {
                break;
            }
        }

        ASSERTL0(n != nfields, "Could not find field named p in m_fields");

        Array<OneD, NekDouble> phys = fields[n]->UpdatePhys();

        cnt = 0;
        for (i = 0; i < fields[n]->GetExpSize(); ++i)
        {
            nquad = fields[n]->GetExp(i)->GetTotPoints();
            Vmath::Fill(nquad, cfl[i], &phys[cnt], 1);
            cnt += nquad;
        }

        fields[n]->FwdTransLocalElmt(fields[n]->GetPhys(),
                                     fields[n]->UpdateCoeffs());

        // Need to reset varibale name for output
        session->SetVariable(n, "CFL");

        // Reset session name for output file
        std::string outname = IncNav->GetSessionName();

        outname += "_CFLStep";
        IncNav->ResetSessionName(outname);
        IncNav->Output();
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
