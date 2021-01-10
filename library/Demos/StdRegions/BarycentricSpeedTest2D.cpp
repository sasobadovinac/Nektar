//
// Created by Edward Laughton on 30/09/2020.
//

#include <iostream>
#include <LibUtilities/BasicUtils/Timer.h>

#include "StdDemoSupport.hpp"

// polynomial = x^2 + y^2 - z^2
Array<OneD, NekDouble> EvalPoly(const Array<OneD, const Array<OneD, NekDouble>> &pts)
{
    Array<OneD, NekDouble> ret(pts.size());
    for (int i = 0; i < pts.size(); i++)
    {
        ret[i] = pow(pts[i][0],2) + pow(pts[i][1], 2) - pow(pts[i][2], 2);
    }

    return ret;
}

// derivative in x = 2x
Array<OneD, NekDouble> EvalPolyDerivx(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts.size());
    for (int i = 0; i < pts.size(); i++)
    {
        ret[i] =  2*pts[i][0];
    }

    return ret;
}

// derivative in y = 2y
Array<OneD, NekDouble> EvalPolyDerivy(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts.size());
    for (int i = 0; i < pts.size(); i++)
    {
        ret[i] =  2*pts[i][1];
    }

    return ret;
}

// derivative in z = -2z
Array<OneD, NekDouble> EvalPolyDerivz(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts.size());
    for (int i = 0; i < pts.size(); i++)
    {
        ret[i] =  -2*pts[i][2];
    }

    return ret;
}

int main(int argc, char *argv[])
{
    DemoSupport demo;
    demo.ParseArguments(argc, argv);
    StdExpansion *E = demo.CreateStdExpansion();

    Array<OneD, Array<OneD, NekDouble> > coordsE = demo.GetCoords(E);

    // Get physical values from coordsE
    Array<OneD, Array<OneD, NekDouble>> tmpCoords(coordsE[0].size());
    for (int i = 0; i < coordsE[0].size(); ++i)
    {
        tmpCoords[i] = Array<OneD, NekDouble>(3);
        tmpCoords[i][0] = coordsE[0][i];
        tmpCoords[i][1] = coordsE[1][i];
        tmpCoords[i][2] = 0;
    }
    Array<OneD, NekDouble> physIn = EvalPoly(tmpCoords);

    // Fill coords array for interpolation with nInt x nInt evenly spaced values
    int nInt = 10;
    Array<OneD, Array<OneD, NekDouble>> coordIn(nInt * nInt);
    for (NekDouble i = 0, cnt = 0; i < nInt; ++i)
    {
        for (NekDouble j = 0; j < nInt; ++j, ++cnt)
        {
            coordIn[cnt] = Array<OneD,NekDouble>(3);
            coordIn[cnt][0] = i * 2 / nInt - 1;
            coordIn[cnt][1] = j * 2 / nInt - 1;
            coordIn[cnt][2] = 0;
        }
    }

    const auto totPoints = coordIn.size();

    // Do tests of cycles depending on order
    int nModes = E->GetBasisNumModes(0);
    int totCyc = 1000;
    std::cout << "Num of modes is " << nModes << " therefore running for " << totCyc << " cycles." << std::endl;

    // Calc interpolation matrix every call
    LibUtilities::Timer t;

    Array<OneD, NekDouble> Ephys(totPoints), Ederiv0(totPoints), Ederiv1(totPoints);
    Array<OneD, NekDouble> Sphys(totPoints), Sderiv0(totPoints), Sderiv1(totPoints);
    Array<OneD, NekDouble> Bphys(totPoints), Bderiv0(totPoints), Bderiv1(totPoints);


    Array<OneD, NekDouble> EphysDeriv0(coordsE[0].size()), EphysDeriv1(coordsE[0].size());
    std::cout << "Testing old method, calc interp matrix every cycle" << std::endl;
    t.Start();
    for (int cyc = 0; cyc < totCyc; ++cyc)
    {
        Array<OneD, Array<OneD, DNekMatSharedPtr>>  I(totPoints);
        for (int i = 0; i < totPoints; ++i)
        {
            // Obtain local collapsed coordinate from
            // cartesian coordinate.
            Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(2);
            E->LocCoordToLocCollapsed(coordIn[i],eta);

            // Get Lagrange interpolants.
            I[i] = Array<OneD, DNekMatSharedPtr>(2);
            I[i][0] = E->GetBase()[0]->GetI(eta);
            I[i][1] = E->GetBase()[1]->GetI(eta+1);
        }

        E->PhysDeriv(physIn, EphysDeriv0, EphysDeriv1);
        for (int i = 0; i < totPoints; ++i)
        {
            Ephys[i]   = E->PhysEvaluate(I[i], physIn);
            Ederiv0[i] = E->PhysEvaluate(I[i], EphysDeriv0);
            Ederiv1[i] = E->PhysEvaluate(I[i], EphysDeriv1);
        }
    }
    t.Stop();
    NekDouble timeOld = t.TimePerTest(totCyc);
    std::cout << "Old method: " << t.Elapsed().count() << "s - > " << timeOld << " per cycle (" << totCyc << " cycles)." << std::endl;

    Array<OneD, NekDouble> SphysDeriv0(coordsE[0].size()), SphysDeriv1(coordsE[0].size());
    Array<OneD, Array<OneD, DNekMatSharedPtr>>  I(totPoints);
    for (int i = 0; i < totPoints; ++i)
    {
        // Obtain local collapsed coordinate from
        // cartesian coordinate.
        Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(2);
        E->LocCoordToLocCollapsed(coordIn[i],eta);

        // Get Lagrange interpolants.
        I[i] = Array<OneD, DNekMatSharedPtr>(2);
        I[i][0] = E->GetBase()[0]->GetI(eta);
        I[i][1] = E->GetBase()[1]->GetI(eta+1);
    }

    t.Start();
    for (int cyc = 0; cyc < totCyc; ++cyc)
    {
        E->PhysDeriv(physIn, SphysDeriv0, SphysDeriv1);
        for (int i = 0; i < totPoints; ++i)
        {
            Sphys[i]   = E->PhysEvaluate(I[i], physIn);
            Sderiv0[i] = E->PhysEvaluate(I[i], SphysDeriv0);
            Sderiv1[i] = E->PhysEvaluate(I[i], SphysDeriv1);
        }
    }
    t.Stop();
    NekDouble timePrecalc = t.TimePerTest(totCyc);
    std::cout << "Precalc method: " << t.Elapsed().count() << "s - > " << timePrecalc << " per cycle (" << totCyc << " cycles)." << std::endl;

    Array<OneD, NekDouble> BphysDeriv0(totPoints), BphysDeriv1(totPoints);
    t.Start();
    for (int cyc = 0; cyc < totCyc; ++cyc)
    {
        for (int i = 0; i < totPoints; ++i)
        {
            Bphys[i] = E->PhysEvaluate(coordIn[i], physIn, Bderiv0[i], Bderiv1[i]);
        }
    }
    t.Stop();
    NekDouble timeBary = t.TimePerTest(totCyc);
    std::cout << "New method: " << t.Elapsed().count() << "s - > " << timeBary << " per cycle (" << totCyc << " cycles)." << std::endl;

    Array<OneD, NekDouble> sol(totPoints), sol0(totPoints), sol1(totPoints), sol2(totPoints);
    sol1 = EvalPolyDerivy(coordIn);
    sol0 = EvalPolyDerivx(coordIn);
    sol = EvalPoly(coordIn);

    std::ofstream outfile;
    std::string fileName = LibUtilities::ShapeTypeMap[E->DetShapeType()];
    outfile.open(fileName + "AllDeriv.txt", std::ios_base::app); // append instead of overwrite
    outfile << nModes << " " << timeOld << " " << timePrecalc << " " << timeBary << std::endl;
    outfile.close();
    std::cout << "Saved to file: " << fileName + "AllDeriv.txt" << std::endl;

    //Check error
    std::cout << "\nBarycentric \t\t L2 Error \t Linf Error" << std::endl;
    std::cout << "\tPhys: "   << "\t\t" << E->L2(Bphys, sol) << "\t" << E->Linf(Bphys, sol)   << std::endl;
    std::cout << "\tDeriv0: " << "\t" << E->L2(Bderiv0, sol0) << "\t" << E->Linf(Bderiv0, sol0)   << std::endl;
    std::cout << "\tDeriv1: " << "\t" << E->L2(Bderiv1, sol1) << "\t" << E->Linf(Bderiv1, sol1)   << std::endl;

    std::cout << "\nStore matrix \t\t L2 Error \t Linf Error" << std::endl;
    std::cout << "\tPhys: "   << "\t\t" <<E->L2(Sphys, sol) << "\t" << E->Linf(Sphys, sol)   << std::endl;
    std::cout << "\tDeriv0: " << "\t" << E->L2(Sderiv0, sol0) << "\t" << E->Linf(Sderiv0, sol0)   << std::endl;
    std::cout << "\tDeriv1: " << "\t" << E->L2(Sderiv1, sol1) << "\t" << E->Linf(Sderiv1, sol1)   << std::endl;

    std::cout << "\nCalc every cycle \t L2 Error \t Linf Error" << std::endl;
    std::cout << "\tPhys: "   << "\t\t" << E->L2(Ephys, sol) << "\t" << E->Linf(Ephys, sol)   << std::endl;
    std::cout << "\tDeriv0: " << "\t" << E->L2(Ederiv0, sol0) << "\t" << E->Linf(Ederiv0, sol0)   << std::endl;
    std::cout << "\tDeriv1: " << "\t" << E->L2(Ederiv1, sol1) << "\t" << E->Linf(Ederiv1, sol1)   << std::endl;
}
