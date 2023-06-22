///////////////////////////////////////////////////////////////////////////////
//
// File: Prepacing.cpp
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

#include <CardiacEPSolver/CellModels/CellModel.h>
#include <CardiacEPSolver/Stimuli/Stimulus.h>
#include <MultiRegions/ExpList.h>
#include <SpatialDomains/MeshComponents.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    SpatialDomains::PointGeomSharedPtr vPoint;
    MultiRegions::ExpListSharedPtr vExp;
    LibUtilities::SessionReaderSharedPtr vSession;
    std::string vCellModel;
    CellModelSharedPtr vCell;
    std::vector<StimulusSharedPtr> vStimulus;
    Array<OneD, Array<OneD, NekDouble>> vWsp(1);
    Array<OneD, Array<OneD, NekDouble>> vSol(1);
    NekDouble vDeltaT;
    NekDouble vTime;
    unsigned int nSteps;

    // Create a session reader to read pacing parameters
    vSession = LibUtilities::SessionReader::CreateInstance(argc, argv);
    vSession->InitSession();

    try
    {
        // Construct a field consisting of a single vertex
        vPoint = MemoryManager<SpatialDomains::PointGeom>::AllocateSharedPtr(
            3, 0, 0.0, 0.0, 0.0);
        vExp = MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(vPoint);

        // Get cell model name and create it
        vSession->LoadSolverInfo("CELLMODEL", vCellModel, "");
        ASSERTL0(vCellModel != "", "Cell Model not specified.");

        vCell =
            GetCellModelFactory().CreateInstance(vCellModel, vSession, vExp);
        vCell->Initialise();

        // Load the stimuli
        vStimulus = Stimulus::LoadStimuli(vSession, vExp);

        // Set up solution arrays, workspace and read in parameters
        vSol[0] = Array<OneD, NekDouble>(1, 0.0);
        vWsp[0] = Array<OneD, NekDouble>(1, 0.0);
        vDeltaT = vSession->GetParameter("TimeStep");
        vTime   = 0.0;
        nSteps  = vSession->GetParameter("NumSteps");

        LibUtilities::EquationSharedPtr e =
            vSession->GetFunction("InitialConditions", "u");
        vSol[0][0] = e->Evaluate(0.0, 0.0, 0.0, 0.0);

        cout << "#";
        for (unsigned int i = 0; i < vCell->GetNumCellVariables(); ++i)
        {
            cout << "   " << vCell->GetCellVarName(i);
        }
        cout << endl;

        // Time integrate cell model
        for (unsigned int i = 0; i < nSteps; ++i)
        {
            // Compute J_ion
            vCell->TimeIntegrate(vSol, vWsp, vTime);

            // Add stimuli J_stim
            for (unsigned int i = 0; i < vStimulus.size(); ++i)
            {
                vStimulus[i]->Update(vWsp, vTime);
            }

            // Time-step with forward Euler
            Vmath::Svtvp(1, vDeltaT, vWsp[0], 1, vSol[0], 1, vSol[0], 1);

            // Increment time
            vTime += vDeltaT;

            // Output current solution to stdout
            cout << vTime << "   " << vSol[0][0];
            for (unsigned int j = 0; j < vCell->GetNumCellVariables(); ++j)
            {
                cout << "   " << vCell->GetCellSolution(j)[0];
            }
            cout << endl;
        }

        for (unsigned int i = 0; i < vCell->GetNumCellVariables(); ++i)
        {
            cout << "# " << vCell->GetCellVarName(i) << "  "
                 << vCell->GetCellSolution(i)[0] << endl;
        }
    }
    catch (...)
    {
        cerr << "An error occured" << endl;
    }

    return 0;
}
