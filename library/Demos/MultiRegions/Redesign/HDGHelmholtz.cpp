#include <cstdio>
#include <cstdlib>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/DisContField.h>
#include <MultiRegions/NekField/NekField.hpp>
#include <SpatialDomains/MeshGraph.h>

//#define TIMING
#ifdef TIMING
#include <time.h>
#define Timing(s)                                                              \
    fprintf(stdout, "%s Took %g seconds\n", s,                                 \
            (clock() - st) / (double)CLOCKS_PER_SEC);                          \
    st = clock();
#else
#define Timing(s) /* Nothing */
#endif

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: HDGHelmholtz  meshfile \n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();


    try
    {
        //----------------------------------------------
        // Read in mesh from input file
        SpatialDomains::MeshGraphSharedPtr graph =
            SpatialDomains::MeshGraph::Read(vSession);
        //----------------------------------------------


        //----------------------------------------------
        // Print summary of solution details
        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
        factors[StdRegions::eFactorTau]    = 1.0;
        const SpatialDomains::ExpansionInfoMap &expansions =
            graph->GetExpansionInfo();
        LibUtilities::BasisKey bkey0 =
            expansions.begin()->second->m_basisKeyVector[0];
        
        if (vComm->GetRank() == 0)
        {
            cout << "Solving Helmholtz (HDG): " << endl;
            cout << "         Communication: " << vSession->GetComm()->GetType()
                 << endl;
            cout << "         Solver type  : "
                 << vSession->GetSolverInfo("GlobalSysSoln") << endl;
            cout << "         Lambda       : " << factors[StdRegions::eFactorLambda]
                 << endl;
            cout << "         No. modes    : " << bkey0.GetNumModes() << endl;
            cout << endl;
        }
        
        //----------------------------------------------
        

        //----------------------------------------------
        // Define Expansion
        int nvars = vSession->GetVariables().size();
        Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nvars);
        for(int v = 0; v < nvars; ++v)
        {
            Exp[v] = MemoryManager<MultiRegions::DisContField>::AllocateSharedPtr(
                                    vSession, graph, vSession->GetVariable(v));
        }
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

        //-----------------------------------------------
        // See if there is an exact solution, if so
        // evaluate and plot errors
        for(int v = 0; v < nvars; ++v)
        {
            LibUtilities::EquationSharedPtr ex_sol =
                vSession->GetFunction("ExactSolution", v);

            if (ex_sol)
            {
                Array<OneD, NekDouble> tmp =
                    Array<OneD, NekDouble>(Exp[0]->GetTotPoints());
                
                //----------------------------------------------
                // evaluate exact solution
                ex_sol->Evaluate(xc, tmp = Fce.UpdateArray1D(v));
                //----------------------------------------------
                
                //--------------------------------------------
                // Calculate error
                NekDouble vLinfError = Exp[0]->Linf(Phys, Fce, v);
                NekDouble vL2Error   = Exp[0]->L2  (Phys, Fce, v);
                NekDouble vH1Error   = Exp[0]->H1  (Phys, Fce, v);

                if (vSession->GetComm()->GetRank() == 0)
                {
                    std::string var = vSession->GetVariable(v); 
                    cout << "L inf error (variable "<< var << ") : " <<
                        vLinfError << endl;
                    cout << "L 2 error   (variable "<< var << ") : " <<
                        vL2Error << endl;
                    cout << "H 1 error   (variable "<< var << ") : " <<
                        vH1Error << endl;
                }
            }
        }
        //--------------------------------------------
    }
    catch (const std::runtime_error &)
    {
        cerr << "Caught exception." << endl;
        return 1;
    }
    
    vSession->Finalise();

    return 0;
}
