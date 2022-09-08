#include <cstdio>
#include <cstdlib>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <MultiRegions/NekField/NekField.hpp>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

int NoCaseStringCompare(const string &s1, const string &s2);

int main(int argc, char *argv[])
{

    if (argc < 2)
    {
        fprintf(stderr, "Usage: Helmholtz3DHomo2D meshfile [SysSolnType]   \n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();
    string meshfile(argv[1]);

    try
    {
        //----------------------------------------------
        // Read in mesh from input file
        SpatialDomains::MeshGraphSharedPtr graph =
            SpatialDomains::MeshGraph::Read(vSession);
        //----------------------------------------------

        //----------------------------------------------
        // Define Expansion
        int nypoints;
        int nzpoints;
        NekDouble ly;
        NekDouble lz;

        vSession->LoadParameter("HomModesY", nypoints);
        vSession->LoadParameter("HomModesZ", nzpoints);
        vSession->LoadParameter("LY", ly);
        vSession->LoadParameter("LZ", lz);
    
        bool useFFT  = false;
        bool deal    = false;

        const LibUtilities::PointsKey PkeyY(nypoints,
                                        LibUtilities::eFourierEvenlySpaced);
        const LibUtilities::BasisKey BkeyY(LibUtilities::eFourier, nypoints, PkeyY);

        const LibUtilities::PointsKey PkeyZ(nzpoints,
                                        LibUtilities::eFourierEvenlySpaced);
        const LibUtilities::BasisKey BkeyZ(LibUtilities::eFourier, nzpoints, PkeyZ);
 
        int nvars = vSession->GetVariables().size();
        Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nvars);
        for(int v = 0; v < nvars; ++v)
        {
            Exp[v] = MemoryManager<MultiRegions::ContField3DHomogeneous2D>::
                AllocateSharedPtr(vSession, BkeyY, BkeyZ, lz, lz, useFFT, deal,
                                  graph,  vSession->GetVariable(v));
        }

        //----------------------------------------------
        
        //----------------------------------------------
        // Print summary of solution details
        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
        const SpatialDomains::ExpansionInfoMap &expansions =
            graph->GetExpansionInfo();
        LibUtilities::BasisKey bkey0 =
            expansions.begin()->second->m_basisKeyVector[0];
        
        cout << "Solving 3D Helmholtz (Homogeneous in yz-plane):" << endl;
        cout << "         Lambda          : " << factors[StdRegions::eFactorLambda]
             << endl;
        cout << "         Ly              : " << ly << endl;
        cout << "         Lz              : " << lz << endl;
        cout << "         N.modes         : " << bkey0.GetNumModes() << endl;
        cout << "         N.Y homo modes  : " << BkeyY.GetNumModes() << endl;
        cout << "         N.Z homo modes  : " << BkeyZ.GetNumModes() << endl;
        cout << endl;
        //----------------------------------------------

        //----------------------------------------------
        // Set up coordinates of mesh for Forcing function evaluation
        NekField<NekDouble, ePhys> xc(Exp[0], 0.0, 3);
        Exp[0]->GetCoords(xc);
        //----------------------------------------------

        //----------------------------------------------
        // Define forcing function for first variable defined in file
        NekField<NekDouble, ePhys> Fce(Exp, 0.0);
        LibUtilities::EquationSharedPtr ffunc;
        Array<OneD, NekDouble> tmp;
        for(int v = 0; v < nvars; ++v)
        {
            ffunc = vSession->GetFunction("Forcing",  vSession->GetVariable(v));
            ffunc->Evaluate(xc, tmp = Fce.UpdateArray1D(v));
        }

        //----------------------------------------------
        // Helmholtz solution taking physical forcing after setting
        // initial condition to zero
        NekField<NekDouble, eCoeff> Coeffs(Exp, 0.0);
        for(int v = 0; v < nvars; ++v)
        {
            Exp[v]->HelmSolve(v, Fce, Coeffs, factors);
        }
        //----------------------------------------------

        //----------------------------------------------
        // Backward Transform Solution to get solved values at
        NekField<NekDouble, ePhys> Phys(Exp, 0.0);
        Exp[0]->BwdTrans(Coeffs, Phys);
        //----------------------------------------------

        //-----------------------------------------------
        // Write solution to file
        LibUtilities::FieldIOSharedPtr fld =
            LibUtilities::FieldIO::CreateDefault(vSession);
        string out = vSession->GetSessionName() + ".fld";
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
            Exp[0]->GetFieldDefinitions();
        std::vector<std::vector<NekDouble>> FieldData(FieldDef.size());
        for (int v = 0; v < nvars; ++v)
        {
            for (int i = 0; i < FieldDef.size(); ++i)
            {
                FieldDef[i]->m_fields.push_back(vSession->GetVariable(v));
                Exp[v]->AppendFieldData(FieldDef[i], FieldData[i],
                                        tmp = Coeffs.UpdateArray1D(v));
            }
        }
        fld->Write(out, FieldDef, FieldData);
        //-----------------------------------------------

        //----------------------------------------------
        // See if there is an exact solution, if so
        // evaluate and plot errors
        for(int v = 0; v < nvars; ++v)
        {
            LibUtilities::EquationSharedPtr ex_sol =
                vSession->GetFunction("ExactSolution", v);

            if (ex_sol)
            {
                //----------------------------------------------
                // evaluate exact solution
                ex_sol->Evaluate(xc, tmp = Fce.UpdateArray1D(v));
                //----------------------------------------------

                //--------------------------------------------
                // Calculate errors
                NekDouble vLinfError = Exp[0]->Linf(Phys, Fce, v);
                NekDouble vL2Error   = Exp[0]->L2  (Phys, Fce, v);
                if (vComm->GetRank() == 0)
                {
                    cout << "L inf error (variable "<< vSession->GetVariable(v) <<
                        ") : " << vLinfError << endl;
                    cout << "L 2 error   (variable "<< vSession->GetVariable(v) <<
                        ") : " << vL2Error << endl;
                }
                //--------------------------------------------
            }
        }    
    }
    catch (const std::runtime_error &)
    {
        cerr << "Caught exception." << endl;
        return 1;
    }

    vSession->Finalise();

    return 0;
}

