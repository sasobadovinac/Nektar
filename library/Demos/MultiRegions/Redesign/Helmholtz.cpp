#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ContField.h>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/FieldStorage/FieldStorage.hpp>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

void SetVariableCoeffs(LibUtilities::SessionReaderSharedPtr &vSession,
		       Array<OneD, NekDouble> &xc0,
		       Array<OneD, NekDouble> &xc1,
		       Array<OneD, NekDouble> &xc2,
		       StdRegions::VarCoeffMap &varcoeffs);

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();
    MultiRegions::ContFieldSharedPtr Exp;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    StdRegions::ConstFactorMap factors;
    StdRegions::VarCoeffMap varcoeffs;

    if( argc < 2 )
    {
        fprintf(stderr,"Usage: Helmholtz  meshfile \n");
        exit(1);
    }

    try
    {
        LibUtilities::FieldIOSharedPtr fld =
            LibUtilities::FieldIO::CreateDefault(vSession);

        //----------------------------------------------
        // Read in mesh from input file
        SpatialDomains::MeshGraphSharedPtr graph =
            SpatialDomains::MeshGraph::Read(vSession);
        //----------------------------------------------

        //----------------------------------------------
        // Print summary of solution details
        factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
        const SpatialDomains::ExpansionInfoMap &expansions = graph->GetExpansionInfo();
        LibUtilities::BasisKey bkey = expansions.begin()->second->m_basisKeyVector[0];

        if (vComm->GetRank() ==0)
        {
            cout << "Solving    Helmholtz: "  << endl;
            cout << "       Communication: " << vComm->GetType() << endl;
            cout << "       Solver type  : " << vSession->GetSolverInfo("GlobalSysSoln") << endl;
            cout << "       Lambda       : " << factors[StdRegions::eFactorLambda] << endl;
            cout << "       No. modes    : " << bkey.GetNumModes() << endl;
        }
        //----------------------------------------------

        //----------------------------------------------
        // Define Expansion
        Exp = MemoryManager<MultiRegions::ContField>::
            AllocateSharedPtr(vSession,graph,vSession->GetVariable(0));
        //----------------------------------------------

        //----------------------------------------------
        // Set up coordinates of mesh for Forcing function evaluation
        coordim = Exp->GetCoordim(0);
        nq      = Exp->GetTotPoints();

        xc0 = Array<OneD,NekDouble>(nq);
        xc1 = Array<OneD,NekDouble>(nq);
        xc2 = Array<OneD,NekDouble>(nq);

        switch(coordim)
        {
        case 1:
            Exp->GetCoords(xc0);
            Vmath::Zero(nq,&xc1[0],1);
            Vmath::Zero(nq,&xc2[0],1);
            break;
        case 2:
            Exp->GetCoords(xc0,xc1);
            Vmath::Zero(nq,&xc2[0],1);
            break;
        case 3:
            Exp->GetCoords(xc0,xc1,xc2);
            break;
        }
        //----------------------------------------------

        //----------------------------------------------
        // Set up variable coefficients if defined

	SetVariableCoeffs(vSession, xc0, xc1, xc2, varcoeffs);
        //----------------------------------------------
        
        //----------------------------------------------
        // Define forcing function for first variable defined in file
        LibUtilities::EquationSharedPtr ffunc
                                        = vSession->GetFunction("Forcing", 0);
	MultiRegions::FieldStorage<NekDouble,MultiRegions::ePhys> Fce(Exp); 
        ffunc->Evaluate(xc0,xc1,xc2, Fce.UpdateData());
        //----------------------------------------------

        //----------------------------------------------
        //Helmholtz solution taking physical forcing after setting
        //initial condition to zero
	MultiRegions::FieldStorage<NekDouble,MultiRegions::eCoeff> Coeffs(Exp); 
        Vmath::Zero(Exp->GetNcoeffs(),Coeffs.UpdateData(),1);
        Exp->HelmSolve(Fce, Coeffs, factors, varcoeffs);
        //----------------------------------------------

        //----------------------------------------------
        // Backward Transform Solution to get solved values at
	MultiRegions::FieldStorage<NekDouble,MultiRegions::ePhys> Phys(Exp); 
        Exp->BwdTrans(Coeffs, Phys);
        //Exp->FwdTrans(Phys,Coeffs);
        //----------------------------------------------

        //-----------------------------------------------
        // Write solution to file
        string out = vSession->GetSessionName() + ".fld";
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                                                    = Exp->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        for(i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back("u");
            Exp->AppendFieldData(FieldDef[i], FieldData[i]);
        }
        fld->Write(out, FieldDef, FieldData);
        //-----------------------------------------------

        //----------------------------------------------
        // See if there is an exact solution, if so
        // evaluate and plot errors
        LibUtilities::EquationSharedPtr ex_sol
                                = vSession->GetFunction("ExactSolution", 0);

        if(ex_sol)
        {
            //----------------------------------------------
            // evaluate exact solution

            ex_sol->Evaluate(xc0,xc1,xc2, Fce.UpdateData());
            //----------------------------------------------

            //--------------------------------------------
            // Calculate errors
            NekDouble vLinfError = Exp->Linf(Phys.GetData(), Fce.GetData());
            NekDouble vL2Error   = Exp->L2(Phys.GetData(), Fce.GetData());
            NekDouble vH1Error   = Exp->H1(Phys.GetData(), Fce.GetData());
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
    catch (const std::runtime_error&)
    {
        cerr << "Caught exception." << endl;
        return 1;
    }

    vComm->Finalise();

    return 0;
}

void SetVariableCoeffs(LibUtilities::SessionReaderSharedPtr &vSession,
		       Array<OneD, NekDouble> &xc0,
		       Array<OneD, NekDouble> &xc1,
		       Array<OneD, NekDouble> &xc2,
		       StdRegions::VarCoeffMap &varcoeffs)
{
    int nq = xc0.size();
    if (vSession->DefinesFunction("d00"))
    {
         Array<OneD, NekDouble> d00(nq,0.0);
	 LibUtilities::EquationSharedPtr d00func =
	   vSession->GetFunction("d00",0);
	 d00func->Evaluate(xc0, xc1, xc2, d00);
	 varcoeffs[StdRegions::eVarCoeffD00] = d00;
    }
        
    if (vSession->DefinesFunction("d01"))
    {
        Array<OneD, NekDouble> d01(nq,0.0);
	LibUtilities::EquationSharedPtr d01func =
	  vSession->GetFunction("d01",0);
	d01func->Evaluate(xc0, xc1, xc2, d01);
	varcoeffs[StdRegions::eVarCoeffD01] = d01;
    }
    
    if (vSession->DefinesFunction("d02"))
    {
        Array<OneD, NekDouble> d02(nq,0.0);
	LibUtilities::EquationSharedPtr d02func =
	  vSession->GetFunction("d02",0);
	d02func->Evaluate(xc0, xc1, xc2, d02);
	varcoeffs[StdRegions::eVarCoeffD02] = d02;
    }
	
    if (vSession->DefinesFunction("d11"))
    {
        Array<OneD, NekDouble> d11(nq,0.0);
	LibUtilities::EquationSharedPtr d11func =
	  vSession->GetFunction("d11",0);
	d11func->Evaluate(xc0, xc1, xc2, d11);
	varcoeffs[StdRegions::eVarCoeffD11] = d11;
    }

    if (vSession->DefinesFunction("d12"))
    {
        Array<OneD, NekDouble> d12(nq,0.0);
	LibUtilities::EquationSharedPtr d12func =
	  vSession->GetFunction("d12",0);
	d12func->Evaluate(xc0, xc1, xc2, d12);
	varcoeffs[StdRegions::eVarCoeffD12] = d12;
    }
	
    if (vSession->DefinesFunction("d22"))
    {
        Array<OneD, NekDouble> d22(nq,0.0);
	LibUtilities::EquationSharedPtr d22func =
	  vSession->GetFunction("d22",0);
	d22func->Evaluate(xc0, xc1, xc2, d22);
	varcoeffs[StdRegions::eVarCoeffD22] = d22;
    }
}
