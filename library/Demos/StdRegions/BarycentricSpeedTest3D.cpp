//
// Created by Edward Laughton on 30/09/2020.
//

#include <iostream>
#include <LibUtilities/BasicUtils/Timer.h>

#include "StdDemoSupport.hpp"

// polynomial = x^2 + y^2 - z^2
Array<OneD, NekDouble> EvalPoly(const Array<OneD, const Array<OneD, NekDouble>> &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    unsigned dim = pts.size();
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] = pow(pts[0][i],2)
                 + (dim >= 2 ? pow(pts[1][i], 2) : 0.0)
                 - (dim >= 3 ? pow(pts[2][i], 2) : 0.0);
    }

    return ret;
}

// derivative in x = 2x
Array<OneD, NekDouble> EvalPolyDerivx(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  2*pts[0][i];
    }

    return ret;
}

// derivative in y = 2y
Array<OneD, NekDouble> EvalPolyDerivy(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  2*pts[1][i];
    }

    return ret;
}

// derivative in z = -2z
Array<OneD, NekDouble> EvalPolyDerivz(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    Array<OneD, NekDouble> ret(pts[0].size());
    for (int i = 0; i < pts[0].size(); i++)
    {
        ret[i] =  -2*pts[2][i] ;
    }

    return ret;
}

int main(int argc, char *argv[])
{
    DemoSupport demo;
    demo.ParseArguments(argc, argv);
    StdExpansion *E = demo.CreateStdExpansion();

    const auto dimension = (unsigned) E->GetShapeDimension();
    const auto totPoints = (unsigned) E->GetTotPoints();

    // Create a new element but with the evenly-spaced points type, so that we
    // perform a PhysEvaluateDeriv at a different set of nodal points
    // (i.e. non-collocated interpolation).
    vector<string> &ptypes = demo.GetPointsType();
    ptypes[0] = "GaussRadauPLegendre";
    ptypes[1] = "GaussRadauMLegendre";
    ptypes[2] = "GaussRadauMLegendre";
    StdExpansion *F = demo.CreateStdExpansion();

    Array<OneD, Array<OneD, NekDouble> > coordsF = demo.GetCoords(F);
    Array<OneD, Array<OneD, NekDouble> > coordsE = demo.GetCoords(E);

    // Get physical values
    Array<OneD, NekDouble> physIn = EvalPoly(coordsE);

    // Fill coords array
    Array<OneD, Array<OneD, NekDouble>> coordIn(totPoints);
    for (int i = 0; i < totPoints; ++i)
    {
        coordIn[i] = Array<OneD, NekDouble>(dimension, 0.0);
        for (int j =0; j < dimension; ++j)
        {
            coordIn[i][j] = coordsF[j][i];
        }
    }

    // Do tests of cycles depending on order
    int nModes = E->GetBasisNumModes(0);
    int totCyc;
    if (nModes < 6)
    {
        totCyc = 4333504.20338048 / pow (nModes, 2.98768871615519);
    }
    else
    {
        totCyc = 123492872.476171 / pow(nModes, 4.95727898146688);
    }

    if (totCyc < 200) {totCyc = 200;} // Set minimum cycles to 200
    std::cout << "Num of modes is " << nModes << " therefore running for " << totCyc << " cycles." << std::endl;

    // Calc interpolation matrix every call
    LibUtilities::Timer t;

    Array<OneD, NekDouble> Ephys(totPoints), Ederiv0(totPoints), Ederiv1(totPoints), Ederiv2(totPoints);
    Array<OneD, NekDouble> Sphys(totPoints), Sderiv0(totPoints), Sderiv1(totPoints), Sderiv2(totPoints);
    Array<OneD, NekDouble> Bphys(totPoints), Bderiv0(totPoints), Bderiv1(totPoints), Bderiv2(totPoints);


    Array<OneD, NekDouble> EphysDeriv0(totPoints), EphysDeriv1(totPoints), EphysDeriv2(totPoints);
    t.Start();
    for (int cyc = 0; cyc < totCyc; ++cyc)
    {
        Array<OneD, Array<OneD, DNekMatSharedPtr>>  I(totPoints);
        for (int i = 0; i < totPoints; ++i)
        {
            // Obtain local collapsed coordinate from
            // cartesian coordinate.
            Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(3);
            E->LocCoordToLocCollapsed(coordIn[i],eta);

            // Get Lagrange interpolants.
            I[i] = Array<OneD, DNekMatSharedPtr>(3);
            I[i][0] = E->GetBase()[0]->GetI(eta);
            I[i][1] = E->GetBase()[1]->GetI(eta+1);
            I[i][2] = E->GetBase()[2]->GetI(eta+2);
        }

        E->PhysDeriv(physIn, EphysDeriv0, EphysDeriv1, EphysDeriv2);
        for (int i = 0; i < totPoints; ++i)
        {
            Ephys[i]   = E->PhysEvaluate(I[i], physIn);
            Ederiv0[i] = E->PhysEvaluate(I[i], EphysDeriv0);
            Ederiv1[i] = E->PhysEvaluate(I[i], EphysDeriv1);
            Ederiv2[i] = E->PhysEvaluate(I[i], EphysDeriv2);
        }
    }
    t.Stop();
    NekDouble timeOld = t.TimePerTest(totCyc);
    std::cout << "Old method: " << t.Elapsed().count() << "s - > " << timeOld << " per cycle (" << totCyc << " cycles)." << std::endl;

    Array<OneD, NekDouble> SphysDeriv0(totPoints), SphysDeriv1(totPoints), SphysDeriv2(totPoints);
    Array<OneD, Array<OneD, DNekMatSharedPtr>>  I(totPoints);
    for (int i = 0; i < totPoints; ++i)
    {
        // Obtain local collapsed coordinate from
        // cartesian coordinate.
        Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(3);
        E->LocCoordToLocCollapsed(coordIn[i],eta);

        // Get Lagrange interpolants.
        I[i] = Array<OneD, DNekMatSharedPtr>(3);
        I[i][0] = E->GetBase()[0]->GetI(eta);
        I[i][1] = E->GetBase()[1]->GetI(eta+1);
        I[i][2] = E->GetBase()[2]->GetI(eta+2);
    }

    t.Start();
    for (int cyc = 0; cyc < totCyc; ++cyc)
    {
        E->PhysDeriv(physIn, SphysDeriv0, SphysDeriv1, SphysDeriv2);
        for (int i = 0; i < totPoints; ++i)
        {
            Sphys[i]   = E->PhysEvaluate(I[i], physIn);
            Sderiv0[i] = E->PhysEvaluate(I[i], SphysDeriv0);
            Sderiv1[i] = E->PhysEvaluate(I[i], SphysDeriv1);
            Sderiv2[i] = E->PhysEvaluate(I[i], SphysDeriv2);
        }
    }
    t.Stop();

    NekDouble timePrecalc = t.TimePerTest(totCyc);
    std::cout << "Precalc method: " << t.Elapsed().count() << "s - > " << timePrecalc << " per cycle (" << totCyc << " cycles)." << std::endl;

    Array<OneD, NekDouble> BphysDeriv0(totPoints), BphysDeriv1(totPoints), BphysDeriv2(totPoints);
    t.Start();
    for (int cyc = 0; cyc < totCyc; ++cyc)
    {
        for (int i = 0; i < totPoints; ++i)
        {
            Bphys[i] = E->PhysEvaluate(coordIn[i], physIn, Bderiv0[i], Bderiv1[i], Bderiv2[i]);
        }
    }
    t.Stop();
    NekDouble timeBary = t.TimePerTest(totCyc);
    std::cout << "New method: " << t.Elapsed().count() << "s - > " << timeBary << " per cycle (" << totCyc << " cycles)." << std::endl;

    Array<OneD, NekDouble> sol(totPoints), sol0(totPoints), sol1(totPoints), sol2(totPoints);
    sol2 = EvalPolyDerivz(coordsF);
    sol1 = EvalPolyDerivy(coordsF);
    sol0 = EvalPolyDerivx(coordsF);
    sol = EvalPoly(coordsF);

    //Check error
    std::cout << "\nBarycentric \t\t L2 Error \t Linf Error" << std::endl;
    std::cout << "\tPhys: "   << "\t\t" << F->L2(Bphys, sol) << "\t" << F->Linf(Bphys, sol)   << std::endl;
    std::cout << "\tDeriv0: " << "\t" << F->L2(Bderiv0, sol0) << "\t" << F->Linf(Bderiv0, sol0)   << std::endl;
    std::cout << "\tDeriv1: " << "\t" << F->L2(Bderiv1, sol1) << "\t" << F->Linf(Bderiv1, sol1)   << std::endl;
    std::cout << "\tDeriv2: " << "\t" << F->L2(Bderiv2, sol2) << "\t" << F->Linf(Bderiv2, sol2)   << std::endl;

    std::cout << "\nStore matrix \t\t L2 Error \t Linf Error" << std::endl;
    std::cout << "\tPhys: "   << "\t\t" << F->L2(Sphys, sol) << "\t" << F->Linf(Sphys, sol)   << std::endl;
    std::cout << "\tDeriv0: " << "\t" << F->L2(Sderiv0, sol0) << "\t" << F->Linf(Sderiv0, sol0)   << std::endl;
    std::cout << "\tDeriv1: " << "\t" << F->L2(Sderiv1, sol1) << "\t" << F->Linf(Sderiv1, sol1)   << std::endl;
    std::cout << "\tDeriv2: " << "\t" << F->L2(Sderiv2, sol2) << "\t" << F->Linf(Sderiv2, sol2)   << std::endl;

    std::cout << "\nCalc every cycle \t L2 Error \t Linf Error" << std::endl;
    std::cout << "\tPhys: "   << "\t\t" << F->L2(Ephys, sol) << "\t" << F->Linf(Ephys, sol)   << std::endl;
    std::cout << "\tDeriv0: " << "\t" << F->L2(Ederiv0, sol0) << "\t" << F->Linf(Ederiv0, sol0)   << std::endl;
    std::cout << "\tDeriv1: " << "\t" << F->L2(Ederiv1, sol1) << "\t" << F->Linf(Ederiv1, sol1)   << std::endl;
    std::cout << "\tDeriv2: " << "\t" << F->L2(Ederiv2, sol2) << "\t" << F->Linf(Ederiv2, sol2)   << std::endl;

    std::ofstream outfile;
    std::string fileName = LibUtilities::ShapeTypeMap[E->DetShapeType()];
    outfile.open(fileName + "All.txt", std::ios_base::app); // append instead of overwrite
    outfile << nModes << " " << timeOld << " " << timePrecalc << " " << timeBary << std::endl;
    outfile.close();
    std::cout << "Saved to file: " << fileName + "All.txt" << std::endl;
}
