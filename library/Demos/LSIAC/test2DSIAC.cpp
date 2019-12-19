#include <cstdio>
#include <cstdlib>
#include <iostream>

#include <LSIAC/CentralBSplines.h>
#include <LSIAC/GeneralBSplines.h>
#include <LSIAC/HandleNekMesh.h>
#include <LSIAC/HandleNekMesh2D.h>
#include <LSIAC/SmoothieSIAC2D.h>
#include <LSIAC/SymmetricSIAC.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <cmath>

using namespace std;
using namespace Nektar;
using namespace LSIAC;
using namespace SIACUtilities;

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        cout << "please enter the following information" << endl;
        cout << "1st arg xml file." << endl;
        cout << "2nd arg polynomial degree filter you want to apply" << endl;
        return 0;
    }

    argc = 2;
    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(argc, argv);
    HandleNekMesh2D *HNM2D = new HandleNekMesh2D(vSession);
    vector<string> var     = vSession->GetVariables();
    HNM2D->LoadMesh(var[0]);
    string fname = vSession->GetSessionName();

    SpatialDomains::MeshGraphSharedPtr graph2D =
        SpatialDomains::MeshGraph::Read(vSession);
    //	graph2D->SetExpansionsToPolyOrder( 2*(atoi(argv[2])-1)+2);
    MultiRegions::ExpListSharedPtr Exp_u =
        MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(
            vSession, graph2D, vSession->GetVariable(0));

    HNM2D->CalculateDynamicScaling();

    int tNquadPts = HNM2D->m_expansions[0]->GetTotPoints();
    Array<OneD, NekDouble> xc0(tNquadPts);
    Array<OneD, NekDouble> xc1(tNquadPts);
    Array<OneD, NekDouble> xc2(tNquadPts);

    switch (HNM2D->m_expansions[0]->GetCoordim(0))
    {
        case 1:
            HNM2D->m_expansions[0]->GetCoords(xc0);
            Vmath::Zero(tNquadPts, &xc1[0], 1);
            Vmath::Zero(tNquadPts, &xc2[0], 1);
            break;
        case 2:
            HNM2D->m_expansions[0]->GetCoords(xc0, xc1);
            Vmath::Zero(tNquadPts, &xc2[0], 1);
            break;
        case 3:
            HNM2D->m_expansions[0]->GetCoords(xc0, xc1, xc2);
            break;
        default:
            assert(false && "looks dim not taken into account");
            cout << "opps did not plan for this" << endl;
    }

    // Define forcing function for first variable defined in file
    Array<OneD, NekDouble> fce, sce, ece, temp0, temp1;
    fce = Array<OneD, NekDouble>(tNquadPts);
    //	LibUtilities::EquationSharedPtr ffunc
    //					= vSession->GetFunction("ExactSolution",
    // 0); 	ffunc->Evaluate(xc0,xc1,xc2, fce);

    for (int i = 0; i < tNquadPts; i++)
    {
        fce[i] = std::cos(2.0 * M_PI * (xc0[i] + xc1[i]));
    }

    HNM2D->m_expansions[0]->FwdTrans(fce,
                                     HNM2D->m_expansions[0]->UpdateCoeffs());
    HNM2D->m_expansions[0]->BwdTrans(HNM2D->m_expansions[0]->GetCoeffs(),
                                     HNM2D->m_expansions[0]->UpdatePhys());
    ece                               = HNM2D->m_expansions[0]->GetPhys();
    Array<OneD, NekDouble> ece_Coeffs = HNM2D->m_expansions[0]->GetCoeffs();

    HNM2D->m_Arrays.push_back(ece);
	// Not loading Rtree since, it is already included in code.
 	//HNM2D->LoadExpListIntoRTree();  

    // Second Target.
    int new_tNquadPts = Exp_u->GetTotPoints();
    Array<OneD, NekDouble> direction(3, 0.0), glCoord(3, 0.0);
    direction[0] = 1.0;
    Array<OneD, NekDouble> sceN(new_tNquadPts);
    Array<OneD, NekDouble> xcN0(new_tNquadPts), xcN1(new_tNquadPts),
        xcN2(new_tNquadPts);
    Array<OneD, NekDouble> fceN(new_tNquadPts), eceN(new_tNquadPts),
        eIDN(new_tNquadPts);
    SmoothieSIAC2D sm(SIACUtilities::eSYM_2kp1_1SIDED_2kp2, HNM2D,
                      atoi(argv[2]), 0.1);
    NekDouble valY, valZ;
    Exp_u->GetCoords(xcN0, xcN1, xcN2);

    for (int elmID = 0; elmID < Exp_u->GetExpSize(); elmID++)
    {
        int physOffsetN = Exp_u->GetPhys_Offset(elmID);
        LocalRegions::ExpansionSharedPtr lexp =
            HNM2D->m_expansions[0]->GetExp(elmID);
        Array<OneD, NekDouble> u_phys_old = ece.CreateWithOffset(
            ece, HNM2D->m_expansions[0]->GetPhys_Offset(elmID));
        for (int j = 0; j < Exp_u->GetExp(elmID)->GetTotPoints(); j++)
        {
            int index  = j + physOffsetN;
            glCoord[0] = xcN0[index];
			NekDouble dynScaling = HNM2D->GetDynamicScaling(glCoord,elmID);
            sm.EvaluateAt(xcN0[index], xcN1[index], 0.0, sceN[index], valY,
                          valZ, direction, dynScaling, 0);
            fceN[index] = std::cos(2.0 * M_PI * (xcN0[index] + xcN1[index]));
            eceN[index] = lexp->PhysEvaluate(glCoord, u_phys_old);
            eIDN[index] = elmID;
        }
    }

    // Calculate errors
    NekDouble vLinfError = Exp_u->Linf(sceN, fceN);
    NekDouble vL2Error   = Exp_u->L2(sceN, fceN);
    cout << "L infinity error: " << vLinfError << endl;
    cout << "L 2 error:        " << vL2Error << endl;

    vSession->Finalise();

    return 0;
}
