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

void printNekArray(Array<OneD, NekDouble> &ar);
void printNekArray(Array<OneD, NekDouble> &ar, int del);

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        cout << "please enter the following information" << endl;
        cout << "1st arg xml file." << endl;
        cout << "2nd arg polynomial degree filter you want to apply" << endl;
        cout << "3rd arg meshscaling you want to use." << endl;
        cout << "4th arg Resolution of output." << endl;
        return 0;
    }

    argc = 2;
    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(argc, argv);
    HandleNekMesh1D *HNM1D = new HandleNekMesh1D(vSession);
    vector<string> var     = vSession->GetVariables();
    HNM1D->LoadMesh(var[0]);
    string fname = vSession->GetSessionName();

    HNM1D->CalculateDynamicScaling();

    int tNquadPts = HNM1D->m_expansions[0]->GetTotPoints();
    cout << "fc:" << tNquadPts << endl;
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
    //0); 	ffunc->Evaluate(xc0,xc1,xc2, fce);

    for (int i = 0; i < tNquadPts; i++)
    {
        // fce[i] = std::cos(2.0*M_PI*(xc0[i]));
        fce[i] = exp(-0.5 * (xc0[i] - 0.5) * (xc0[i] - 0.5) / 0.025 / 0.025);
        // fce[i] = exp(-0.5* (xc0[i]-0.5)*(xc0[i]-0.5)/0.05/0.05);
    }

    HNM1D->m_expansions[0]->FwdTrans(fce,
                                     HNM1D->m_expansions[0]->UpdateCoeffs());
    HNM1D->m_expansions[0]->BwdTrans(HNM1D->m_expansions[0]->GetCoeffs(),
                                     HNM1D->m_expansions[0]->UpdatePhys());
    ece                               = HNM1D->m_expansions[0]->GetPhys();
    Array<OneD, NekDouble> ece_Coeffs = HNM1D->m_expansions[0]->GetCoeffs();

    HNM1D->m_Arrays.push_back(ece);
    //	HNM1D->LoadExpListIntoRTree(); // CHECK THIS
    //	HNM1D->CalculateDynamicScaling(); // CHECK THIS

    // Second Target.

    // Evaluate on a new equal space grid mesh.
    int gPts     = atoi(argv[4]);
    int totPts   = gPts;
    int Nx       = gPts;
    NekDouble sx = 1.0 / (Nx - 1.0);

    SmoothieSIAC1D sm(SIACUtilities::eSYM_2kp1_1SIDED_2kp2, HNM1D,
                      atoi(argv[2]), atof(argv[3]));
    SmoothieSIAC1D sm_NONSYM(SIACUtilities::eSYM_UNEVEN_2kp1, HNM1D,
                             atoi(argv[2]), atof(argv[3]));
    NekDouble valY, valZ;
    vector<NekDouble> pos_x;
    pos_x.clear();
    vector<NekDouble> values_x;
    values_x.clear();
    Array<OneD, NekDouble> direction(3, 0.0), coord(3, 0.0);
    direction[0] = 1.0;
    Array<OneD, NekDouble> pX(totPts), pY(totPts), pV(totPts), pP(totPts),
        pS(totPts), pDyn(totPts), pSDyn(totPts), pE(totPts);
    Array<OneD, NekDouble> glCoords(3, 0.0), pS_NSYM(totPts);
    for (int i = 0; i < Nx; i++)
    {
        int index = i;
        pX[index] = i * sx;
        // pV[index] = std::cos(2.0*M_PI*(pX[index]) );
        pV[index] =
            exp(-0.5 * (pX[index] - 0.5) * (pX[index] - 0.5) / 0.025 / 0.025);
        // pV[index] = exp(-0.5* (pX[index]-0.5)*(pX[index]-0.5)/0.05/0.05);
        glCoords[0] = pX[index];
        //			NekDouble dynScaling =
        //HNM1D->GetDynamicScaling(glCoords); // CHECK THIS 			pDyn[index] =
        //dynScaling; // CHECK THIS
        sm.EvaluateAt(pX[index], valY, valZ, pS[index], valY, valZ);
        sm_NONSYM.EvaluateNonSymAt(pX[index], 0.0, 0.0, pS_NSYM[index], valY,
                                   valZ);
        // sm.EvaluateAt(pX[index],pY[index],0.0,pSDyn[index],valY,valZ,direction,
        // dynScaling ,0);  // CHECK THIS
        coord[0]           = pX[index];
        coord[1]           = 0.0;
        coord[2]           = 0.0;
        int elid           = HNM1D->m_expansions[0]->GetExpIndex(coord);
        NekDouble dScaling = HNM1D->GetDynamicScaling(coord, elid, 1.0);
        pDyn[index]        = dScaling;
        sm.EvaluateAt(pX[index], 0.0, 0.0, pSDyn[index], valY, valZ, direction,
                      dScaling, 0);
        LocalRegions::ExpansionSharedPtr lexp =
            HNM1D->m_expansions[0]->GetExp(elid);
        // int elId = HNM1D->m_expansions[0]->GetExpIndex(coord);
        int physOffset = HNM1D->m_expansions[0]->GetPhys_Offset(elid);
        const Array<OneD, NekDouble> el_phys =
            ece.CreateWithOffset(ece, physOffset);
        pP[index] = lexp->PhysEvaluate(coord, el_phys);
        pE[index] = elid;
        //			cout << j << endl;

        // cout << coord[0] << "\t"<<coord[1] << "\t"<<coord[2] <<endl;
        //			cout << pV[index] << "\t"<<pP[index] <<
        //"\t"<<pS[index] << "\t"<<pS_NSYM[index] <<endl;
    }

    LSIACPostProcessor k;
    k.writeNekArray(pX, fname + "_" + argv[2] + "_SExp_" + argv[3] + "_R_" +
                            argv[4] + "_pX_2DDyn.txt");
    k.writeNekArray(pE, fname + "_" + argv[2] + "_SExp_" + argv[3] + "_R_" +
                            argv[4] + "_pE_2DDyn.txt");
    k.writeNekArray(pV, fname + "_" + argv[2] + "_SExp_" + argv[3] + "_R_" +
                            argv[4] + "_pV_2DDyn.txt");
    k.writeNekArray(pP, fname + "_" + argv[2] + "_SExp_" + argv[3] + "_R_" +
                            argv[4] + "_pP_2DDyn.txt");
    k.writeNekArray(pS, fname + "_" + argv[2] + "_SExp_" + argv[3] + "_R_" +
                            argv[4] + "_pS_2DDyn.txt");
    k.writeNekArray(pS_NSYM, fname + "_" + argv[2] + "_SExp_" + argv[3] +
                                 "_R_" + argv[4] + "_pS_NSYM_2DDyn.txt");
    k.writeNekArray(pSDyn, fname + "_" + argv[2] + "_SExp_" + argv[3] + "_R_" +
                               argv[4] + "_pSDyn_2DDyn.txt");
    k.writeNekArray(pDyn, fname + "_" + argv[2] + "_SExp_" + argv[3] + "_R_" +
                              argv[4] + "_pDyn_2DDyn.txt");
    return 0;
}
