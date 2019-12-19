#include <cstdio>
#include <cstdlib>
#include <iostream>

#include <LSIAC/CentralBSplines.h>
#include <LSIAC/GeneralBSplines.h>
#include <LSIAC/HandleNekMesh.h>
#include <LSIAC/HandleNekMesh1D.h>
#include <LSIAC/SmoothieSIAC1D.h>
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
    HandleNekMesh1D *HNM1D = new HandleNekMesh1D(vSession);
    vector<string> var     = vSession->GetVariables();
    HNM1D->LoadMesh(var[0]);
    string fname = vSession->GetSessionName();

    SpatialDomains::MeshGraphSharedPtr graph1D =
        SpatialDomains::MeshGraph::Read(vSession);
    graph1D->SetExpansionsToPolyOrder(2 * (atoi(argv[2]) - 1) + 2);
    MultiRegions::ExpListSharedPtr Exp_u =
        MemoryManager<MultiRegions::DisContField1D>::AllocateSharedPtr(
            vSession, graph1D, vSession->GetVariable(0));

    HNM1D->CalculateDynamicScaling();

    int tNquadPts = HNM1D->m_expansions[0]->GetTotPoints();
    Array<OneD, NekDouble> xc0(tNquadPts);
    Array<OneD, NekDouble> xc1(tNquadPts);
    Array<OneD, NekDouble> xc2(tNquadPts);

    switch (HNM1D->m_expansions[0]->GetCoordim(0))
    {
        case 1:
            HNM1D->m_expansions[0]->GetCoords(xc0);
            Vmath::Zero(tNquadPts, &xc1[0], 1);
            Vmath::Zero(tNquadPts, &xc2[0], 1);
            break;
        case 2:
            HNM1D->m_expansions[0]->GetCoords(xc0, xc1);
            Vmath::Zero(tNquadPts, &xc2[0], 1);
            break;
        case 3:
            HNM1D->m_expansions[0]->GetCoords(xc0, xc1, xc2);
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
        fce[i] = std::cos(2.0 * M_PI * (xc0[i]));
    }

    HNM1D->m_expansions[0]->FwdTrans(fce,
                                     HNM1D->m_expansions[0]->UpdateCoeffs());
    HNM1D->m_expansions[0]->BwdTrans(HNM1D->m_expansions[0]->GetCoeffs(),
                                     HNM1D->m_expansions[0]->UpdatePhys());
    ece                               = HNM1D->m_expansions[0]->GetPhys();
    Array<OneD, NekDouble> ece_Coeffs = HNM1D->m_expansions[0]->GetCoeffs();

    HNM1D->m_Arrays.push_back(ece);
    // HNM1D->LoadExpListIntoRTree();

    // Second Target.
    int new_tNquadPts = Exp_u->GetTotPoints();
    Array<OneD, NekDouble> direction(3, 0.0), glCoord(3, 0.0);
    direction[0] = 1.0;
    Array<OneD, NekDouble> sceN(new_tNquadPts);
    Array<OneD, NekDouble> xcN0(new_tNquadPts), xcN1(new_tNquadPts),
        xcN2(new_tNquadPts);
    Array<OneD, NekDouble> fceN(new_tNquadPts), eceN(new_tNquadPts),
        eIDN(new_tNquadPts);
    SmoothieSIAC1D sm(SIACUtilities::eSYM_2kp1_1SIDED_2kp2, HNM1D,
                      atoi(argv[2]), 0.1);
    NekDouble valY, valZ;
    Exp_u->GetCoords(xcN0, xcN1, xcN2);

    for (int elmID = 0; elmID < Exp_u->GetExpSize(); elmID++)
    {
        int physOffsetN = Exp_u->GetPhys_Offset(elmID);
        LocalRegions::ExpansionSharedPtr lexp =
            HNM1D->m_expansions[0]->GetExp(elmID);
        Array<OneD, NekDouble> u_phys_old = ece.CreateWithOffset(
            ece, HNM1D->m_expansions[0]->GetPhys_Offset(elmID));
        for (int j = 0; j < Exp_u->GetExp(elmID)->GetTotPoints(); j++)
        {
            int index  = j + physOffsetN;
            glCoord[0] = xcN0[index];
            sm.EvaluateAt(xcN0[index], 0.0, 0.0, sceN[index], valY, valZ,
                          direction, HNM1D->GetDynamicScaling(glCoord, elmID),
                          0);
            fceN[index] = std::cos(2.0 * M_PI * (xcN0[index]));
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
