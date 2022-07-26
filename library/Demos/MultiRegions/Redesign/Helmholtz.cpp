#include <cstdio>
#include <cstdlib>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ContField.h>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/Field/Field.hpp>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

void SetVariableCoeffs(LibUtilities::SessionReaderSharedPtr &vSession,
                       Field<NekDouble, ePhys> &xc,
                       StdRegions::VarCoeffMap &varcoeffs);

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: Helmholtz  meshfile \n");
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
        const SpatialDomains::ExpansionInfoMap &expansions =
            graph->GetExpansionInfo();
        LibUtilities::BasisKey bkey =
            expansions.begin()->second->m_basisKeyVector[0];

        if (vComm->GetRank() == 0)
        {
            cout << "Solving    Helmholtz: " << endl;
            cout << "       Communication: " << vComm->GetType() << endl;
            cout << "       Solver type  : "
                 << vSession->GetSolverInfo("GlobalSysSoln") << endl;
            cout << "       Lambda       : "
                 << factors[StdRegions::eFactorLambda] << endl;
            cout << "       No. modes    : " << bkey.GetNumModes() << endl;
        }
        //----------------------------------------------

        //----------------------------------------------
        // Define Expansion
        MultiRegions::ContFieldSharedPtr Exp =
            MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
                vSession, graph, vSession->GetVariable(0));
        //----------------------------------------------

        //----------------------------------------------
        // Set up coordinates of mesh for Forcing function evaluation
        Field<NekDouble, ePhys> xc(Exp, 0.0, 3);
        Exp->GetCoords(xc);

        //----------------------------------------------

        //----------------------------------------------
        // Set up variable coefficients if defined
        StdRegions::VarCoeffMap varcoeffs;
        SetVariableCoeffs(vSession, xc, varcoeffs);
        //----------------------------------------------

        //----------------------------------------------
        // Define forcing function for first variable defined in file
        LibUtilities::EquationSharedPtr ffunc =
            vSession->GetFunction("Forcing", 0);
        Field<NekDouble, ePhys> Fce(Exp);
        ffunc->Evaluate(xc, Fce);
        //----------------------------------------------

        //----------------------------------------------
        // Helmholtz solution taking physical forcing after setting
        // initial condition to zero
        Field<NekDouble, eCoeff> Coeffs(Exp);
        Vmath::Zero(Exp->GetNcoeffs(), Coeffs.UpdateData(), 1);
        Exp->HelmSolve(Fce, Coeffs, factors, varcoeffs);
        //----------------------------------------------

        //----------------------------------------------
        // Backward Transform Solution to get solved values at
        Field<NekDouble, ePhys> Phys(Exp);
        Exp->BwdTrans(Coeffs, Phys);
        //----------------------------------------------

        //-----------------------------------------------
        // Write solution to file
        LibUtilities::FieldIOSharedPtr fld =
            LibUtilities::FieldIO::CreateDefault(vSession);
        string out = vSession->GetSessionName() + ".fld";
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
            Exp->GetFieldDefinitions();
        std::vector<std::vector<NekDouble>> FieldData(FieldDef.size());
        for (int i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back("u");
            Exp->AppendFieldData(FieldDef[i], FieldData[i]);
        }
        fld->Write(out, FieldDef, FieldData);
        //-----------------------------------------------

        //----------------------------------------------
        // See if there is an exact solution, if so
        // evaluate and plot errors
        LibUtilities::EquationSharedPtr ex_sol =
            vSession->GetFunction("ExactSolution", 0);

        if (ex_sol)
        {
            //----------------------------------------------
            // evaluate exact solution
            ex_sol->Evaluate(xc, Fce);
            //----------------------------------------------

            //--------------------------------------------
            // Calculate errors
            NekDouble vLinfError = Exp->Linf(Phys, Fce);
            NekDouble vL2Error   = Exp->L2(Phys, Fce);
            NekDouble vH1Error   = Exp->H1(Phys, Fce);
            if (vComm->GetRank() == 0)
            {
                cout << "L infinity error: " << vLinfError << endl;
                cout << "L 2 error:        " << vL2Error << endl;
                cout << "H 1 error:        " << vH1Error << endl;
            }
            //--------------------------------------------
        }
        //----------------------------------------------
    }
    catch (const std::runtime_error &)
    {
        cerr << "Caught exception." << endl;
        return 1;
    }

    vComm->Finalise();

    return 0;
}

void SetVariableCoeffs(LibUtilities::SessionReaderSharedPtr &vSession,
                       Field<NekDouble, ePhys> &xc,
                       StdRegions::VarCoeffMap &varcoeffs)
{
    int varsize = xc.GetArray1D(xc.GetNumVariables()-1).size();
        
    if (vSession->DefinesFunction("d00"))
    {
        Array<OneD, NekDouble> d00(varsize);
        LibUtilities::EquationSharedPtr d00func =
            vSession->GetFunction("d00", 0);
        d00func->Evaluate(xc, d00);
        varcoeffs[StdRegions::eVarCoeffD00] = d00;
    }

    if (vSession->DefinesFunction("d01"))
    {
        Array<OneD, NekDouble> d01(varsize);
        LibUtilities::EquationSharedPtr d01func =
            vSession->GetFunction("d01", 0);
        d01func->Evaluate(xc, d01);
        varcoeffs[StdRegions::eVarCoeffD01] = d01;
    }

    if (vSession->DefinesFunction("d02"))
    {
        Array<OneD, NekDouble> d02(varsize);
        LibUtilities::EquationSharedPtr d02func =
            vSession->GetFunction("d02", 0);
        d02func->Evaluate(xc, d02);
        varcoeffs[StdRegions::eVarCoeffD02] = d02;
    }

    if (vSession->DefinesFunction("d11"))
    {
        Array<OneD, NekDouble> d11(varsize);
        LibUtilities::EquationSharedPtr d11func =
            vSession->GetFunction("d11", 0);
        d11func->Evaluate(xc, d11);
        varcoeffs[StdRegions::eVarCoeffD11] = d11;
    }

    if (vSession->DefinesFunction("d12"))
    {
        Array<OneD, NekDouble> d12(varsize);
        LibUtilities::EquationSharedPtr d12func =
            vSession->GetFunction("d12", 0);
        d12func->Evaluate(xc, d12);
        varcoeffs[StdRegions::eVarCoeffD12] = d12;
    }

    if (vSession->DefinesFunction("d22"))
    {
        Array<OneD, NekDouble> d22(varsize);
        LibUtilities::EquationSharedPtr d22func =
            vSession->GetFunction("d22", 0);
        d22func->Evaluate(xc, d22);
        varcoeffs[StdRegions::eVarCoeffD22] = d22;
    }
}
