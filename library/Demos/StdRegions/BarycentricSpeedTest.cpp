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
    int dim = coordsE.size();

    std::cout << "ptypes[0] = " << LibUtilities::kPointsTypeStr[E->GetPointsType(0)] << std::endl;
    if (dim > 1)
    {
        std::cout << "ptypes[1] = " << LibUtilities::kPointsTypeStr[E->GetPointsType(1)] << std::endl;
    }
    if (dim > 2)
    {
        std::cout << "ptypes[2] = " << LibUtilities::kPointsTypeStr[E->GetPointsType(2)] << std::endl;
    }

    // Get physical values from coordsE
    Array<OneD, Array<OneD, NekDouble>> tmpCoords(coordsE[0].size());
    for (int i = 0; i < coordsE[0].size(); ++i)
    {
        tmpCoords[i] = Array<OneD, NekDouble>(3);
        tmpCoords[i][0] = coordsE[0][i];
        tmpCoords[i][1] = (dim > 1) ? coordsE[1][i] : 0;
        tmpCoords[i][2] = (dim > 2) ? coordsE[2][i] : 0;
    }
    Array<OneD, NekDouble> physIn = EvalPoly(tmpCoords);

    // Create a new element but with the evenly-spaced points type, so that we
    // perform a PhysEvaluateDeriv at a different set of nodal points
    // (i.e. non-collocated interpolation).
    vector<string> &ptypes = demo.GetPointsType();
    ptypes[0] = "PolyEvenlySpaced";
    ptypes[1] = "PolyEvenlySpaced";
    ptypes[2] = "PolyEvenlySpaced";

    vector<int> &points = demo.GetPoints();
    if (dim == 1)
    {
        points[0] = 64;
        std::cout << "num points to sample = 64^1" << std::endl;

    }
    else if (dim == 2)
    {
        points[0] = 8;
        points[1] = 8;
        std::cout << "num points to sample = 8^2" << std::endl;
    }
    else if (dim == 3)
    {
        points[0] = 4;
        points[1] = 4;
        points[2] = 4;
        std::cout << "num points to sample = 4^3" << std::endl;
    }

    StdExpansion *F = demo.CreateStdExpansion();
    auto tmpCoords2 = demo.GetCoords(F);

    Array<OneD, Array<OneD, NekDouble>> coordIn(tmpCoords2[0].size());
    for (int i = 0; i < tmpCoords2[0].size(); ++i)
    {
        coordIn[i] = Array<OneD, NekDouble>(3);
        coordIn[i][0] = tmpCoords2[0][i];
        coordIn[i][1] = (dim > 1) ? tmpCoords2[1][i] : 0;
        coordIn[i][2] = (dim > 2) ? tmpCoords2[2][i] : 0;
    }

    int totPoints = coordIn.size();

    int nModes = E->GetBasisNumModes(0);
    std::cout << "Num of modes is: " << nModes << std::endl;
    int totCyc = 100000;
    if (dim == 1)
    {
        totCyc = 1000000;
    }

    std::cout << "Running timings for " << totCyc << " cycles." << std::endl;

    // Calc interpolation matrix every call
    LibUtilities::Timer t;

    Array<OneD, NekDouble> Ephys(totPoints), Ederiv0(totPoints), E2deriv0(totPoints), Ederiv1(totPoints), Ederiv2(totPoints);
    Array<OneD, NekDouble> Sphys(totPoints), Sderiv0(totPoints), S2deriv0(totPoints), Sderiv1(totPoints), Sderiv2(totPoints);
    Array<OneD, NekDouble> Bphys(totPoints), Bderiv0(totPoints), B2deriv0(totPoints), Bderiv1(totPoints), Bderiv2(totPoints);

    Array<OneD, NekDouble> EphysDeriv0(coordsE[0].size()), Ephys2Deriv0(coordsE[0].size()), EphysDeriv1(coordsE[0].size()), EphysDeriv2(coordsE[0].size());
    std::cout << "Testing old method, calc interp matrix every cycle" << std::endl;

    if (dim == 1)
    {
        t.Start();
        for (int cyc = 0; cyc < totCyc; ++cyc)
        {
            Array<OneD, Array<OneD, DNekMatSharedPtr>> I(totPoints);
            for (int i = 0; i < totPoints; ++i)
            {
                // Obtain local collapsed coordinate from
                // cartesian coordinate.
                Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(3);
                eta[0] = coordIn[i][0];
                eta[1] = coordIn[i][1];
                eta[2] = coordIn[i][2];

                // Get Lagrange interpolants.
                I[i]    = Array<OneD, DNekMatSharedPtr>(1);
                I[i][0] = E->GetBase()[0]->GetI(eta);
            }

            E->PhysDeriv(physIn, EphysDeriv0);
            E->PhysDeriv(EphysDeriv0, Ephys2Deriv0);

            for (int i = 0; i < totPoints; ++i)
            {
                Ephys[i]   = E->PhysEvaluate(I[i], physIn);
                Ederiv0[i] = E->PhysEvaluate(I[i], EphysDeriv0);
                E2deriv0[i] = E->PhysEvaluate(I[i], Ephys2Deriv0);

            }
        }
        t.Stop();
    }
    else if (dim == 2)
    {
        t.Start();
        for (int cyc = 0; cyc < totCyc; ++cyc)
        {
            Array<OneD, Array<OneD, DNekMatSharedPtr>> I(totPoints);
            for (int i = 0; i < totPoints; ++i)
            {
                // Obtain local collapsed coordinate from
                // cartesian coordinate.
                Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(3);
                E->LocCoordToLocCollapsed(coordIn[i], eta);

                // Get Lagrange interpolants.
                I[i]    = Array<OneD, DNekMatSharedPtr>(2);
                I[i][0] = E->GetBase()[0]->GetI(eta);
                I[i][1] = E->GetBase()[1]->GetI(eta + 1);
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
    }
    else if (dim == 3)
    {
        t.Start();
        for (int cyc = 0; cyc < totCyc; ++cyc)
        {
            Array<OneD, Array<OneD, DNekMatSharedPtr>> I(totPoints);
            for (int i = 0; i < totPoints; ++i)
            {
                // Obtain local collapsed coordinate from
                // cartesian coordinate.
                Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(3);
                E->LocCoordToLocCollapsed(coordIn[i], eta);

                // Get Lagrange interpolants.
                I[i]    = Array<OneD, DNekMatSharedPtr>(3);
                I[i][0] = E->GetBase()[0]->GetI(eta);
                I[i][1] = E->GetBase()[1]->GetI(eta + 1);
                I[i][2] = E->GetBase()[2]->GetI(eta + 2);
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
    }

    NekDouble timeOld = t.TimePerTest(totCyc);
    std::cout << "Old method: " << t.Elapsed().count() << "s - > " << timeOld << " per cycle (" << totCyc << " cycles)." << std::endl;

    Array<OneD, NekDouble> SphysDeriv0(coordsE[0].size()), Sphys2Deriv0(coordsE[0].size()), SphysDeriv1(coordsE[0].size()), SphysDeriv2(coordsE[0].size());
    Array<OneD, Array<OneD, DNekMatSharedPtr>>  I(totPoints);

    if (dim == 1)
    {
        for (int i = 0; i < totPoints; ++i)
        {
            // Obtain local collapsed coordinate from
            // cartesian coordinate.
            Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(3);
            eta[0] = coordIn[i][0];
            eta[1] = coordIn[i][1];
            eta[2] = coordIn[i][2];

            // Get Lagrange interpolants.
            I[i]    = Array<OneD, DNekMatSharedPtr>(1);
            I[i][0] = E->GetBase()[0]->GetI(eta);
        }
    }
    else if (dim == 2)
    {
        for (int i = 0; i < totPoints; ++i)
        {
            // Obtain local collapsed coordinate from
            // cartesian coordinate.
            Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(3);
            E->LocCoordToLocCollapsed(coordIn[i], eta);

            // Get Lagrange interpolants.
            I[i]    = Array<OneD, DNekMatSharedPtr>(2);
            I[i][0] = E->GetBase()[0]->GetI(eta);
            I[i][1] = E->GetBase()[1]->GetI(eta + 1);
        }
    }
    else if (dim == 3)
    {
        for (int i = 0; i < totPoints; ++i)
        {
            // Obtain local collapsed coordinate from
            // cartesian coordinate.
            Array<OneD, NekDouble> eta = Array<OneD, NekDouble>(3);
            E->LocCoordToLocCollapsed(coordIn[i], eta);

            // Get Lagrange interpolants.
            I[i]    = Array<OneD, DNekMatSharedPtr>(3);
            I[i][0] = E->GetBase()[0]->GetI(eta);
            I[i][1] = E->GetBase()[1]->GetI(eta + 1);
            I[i][2] = E->GetBase()[2]->GetI(eta + 2);
        }
    }

    if(dim == 1)
    {
        t.Start();
        for (int cyc = 0; cyc < totCyc; ++cyc)
        {
            E->PhysDeriv(physIn, SphysDeriv0);
            E->PhysDeriv(SphysDeriv0, Sphys2Deriv0);
            for (int i = 0; i < totPoints; ++i)
            {
                Sphys[i]   = E->PhysEvaluate(I[i], physIn);
                Sderiv0[i] = E->PhysEvaluate(I[i], SphysDeriv0);
                S2deriv0[i] = E->PhysEvaluate(I[i], Sphys2Deriv0);
            }
        }
        t.Stop();
    }
    else if(dim == 2)
    {
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
    }
    else if (dim == 3)
    {
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
    }

    NekDouble timePrecalc = t.TimePerTest(totCyc);
    std::cout << "Precalc method: " << t.Elapsed().count() << "s - > " << timePrecalc << " per cycle (" << totCyc << " cycles)." << std::endl;

    if (dim == 1)
    {
        t.Start();
        for (int cyc = 0; cyc < totCyc; ++cyc)
        {
            for (int i = 0; i < totPoints; ++i)
            {
                Bphys[i] = E->PhysEvaluate2ndDeriv(coordIn[i], physIn, Bderiv0[i], B2deriv0[i]);
            }
        }
        t.Stop();
    }
    else if (dim == 2)
    {
        t.Start();
        for (int cyc = 0; cyc < totCyc; ++cyc)
        {
            for (int i = 0; i < totPoints; ++i)
            {
                Bphys[i] = E->PhysEvaluate(coordIn[i], physIn, Bderiv0[i], Bderiv1[i]);
            }
        }
        t.Stop();
    }
    else if (dim == 3)
    {
        t.Start();
        for (int cyc = 0; cyc < totCyc; ++cyc)
        {
            for (int i = 0; i < totPoints; ++i)
            {
                Bphys[i] = E->PhysEvaluate(coordIn[i], physIn, Bderiv0[i], Bderiv1[i], Bderiv2[i]);
            }
        }
        t.Stop();
    }

    NekDouble timeBary = t.TimePerTest(totCyc);
    std::cout << "New method: " << t.Elapsed().count() << "s - > " << timeBary << " per cycle (" << totCyc << " cycles)." << std::endl;

    Array<OneD, NekDouble> sol(totPoints), sol0(totPoints), sol1(totPoints), sol2(totPoints), sol2_0(totPoints, 2);
    sol2 = EvalPolyDerivz(coordIn);
    sol1 = EvalPolyDerivy(coordIn);
    sol0 = EvalPolyDerivx(coordIn);
    sol = EvalPoly(coordIn);

    std::ofstream outfile;
    std::string fileName = LibUtilities::ShapeTypeMap[E->DetShapeType()];
    outfile.open(fileName + "AllDeriv.txt", std::ios_base::app); // append instead of overwrite
    outfile << nModes << " " << timeOld << " " << timePrecalc << " " << timeBary << std::endl;
    outfile.close();
    std::cout << "Saved to file: " << fileName + "AllDeriv.txt" << std::endl;

    // Debug output
    for (int i = 0; i < totPoints; ++i)
    {
        //std::cout << sol[i] << " " << Bphys[i] << " " << Sphys[i] << " " << Ephys[i] << " -> Coord = (" << coordIn[i][0] << ", " << coordIn[i][1] << ", " << coordIn[i][2] << ")" << std::endl;
        //std::cout << sol0[i] << " " << Bderiv0[i] << " " << Sderiv0[i] << " " << Ederiv0[i] << " -> Coord = (" << coordIn[i][0] << ", " << coordIn[i][1] << ", " << coordIn[i][2] << ")" << std::endl;
        //std::cout << sol1[i] << " " << Bderiv1[i] << " " << Sderiv1[i] << " " << Ederiv1[i] << " -> Coord = (" << coordIn[i][0] << ", " << coordIn[i][1] << ", " << coordIn[i][2] << ")" << std::endl;
        //std::cout << sol2[i] << " " << Bderiv2[i] << " " << Sderiv2[i] << " " << Ederiv2[i] << " -> Coord = (" << coordIn[i][0] << ", " << coordIn[i][1] << ", " << coordIn[i][2] << ")" << std::endl;
        //std::cout << sol2_0[i] << " " << B2deriv0[i] << " " << S2deriv0[i] << " " << E2deriv0[i] << " -> Coord = (" << coordIn[i][0] << ", " << coordIn[i][1] << ", " << coordIn[i][2] << ")" << std::endl;
        //std::cout << sol0[i] << " " << Bderiv0[i] << "\t -> Coord = (" << coordIn[i][0] << ", " << coordIn[i][1] << ", " << coordIn[i][2] << ")" << std::endl;
    }

    //Check error
    std::cout << "\nBarycentric \t\t L2 Error \t Linf Error" << std::endl;
    std::cout << "\tPhys: "   << "\t\t" << F->L2(Bphys, sol) << "\t" << F->Linf(Bphys, sol)   << std::endl;
    if (dim > 0)
    {
        std::cout << "\tDeriv0: " << "\t" << F->L2(Bderiv0, sol0) << "\t" << F->Linf(Bderiv0, sol0)   << std::endl;
    }
    if (dim == 1)
    {
        std::cout << "\t2Deriv0: " << "\t" << F->L2(B2deriv0, sol2_0) << "\t" << F->Linf(B2deriv0, sol2_0)   << std::endl;
    }
    if (dim > 1)
    {
        std::cout << "\tDeriv1: " << "\t" << F->L2(Bderiv1, sol1) << "\t" << F->Linf(Bderiv1, sol1)   << std::endl;
    }
    if (dim > 2)
    {
        std::cout << "\tDeriv2: " << "\t" << F->L2(Bderiv2, sol2) << "\t" << F->Linf(Bderiv2, sol2)   << std::endl;
    }

    std::cout << "\nStore matrix \t\t L2 Error \t Linf Error" << std::endl;
    std::cout << "\tPhys: "   << "\t\t" <<F->L2(Sphys, sol) << "\t" << F->Linf(Sphys, sol)   << std::endl;
    if (dim > 0)
    {
        std::cout << "\tDeriv0: " << "\t" << F->L2(Sderiv0, sol0) << "\t" << F->Linf(Sderiv0, sol0) << std::endl;
    }
    if (dim == 1)
    {
        std::cout << "\t2Deriv0: " << "\t" << F->L2(S2deriv0, sol2_0) << "\t" << F->Linf(S2deriv0, sol2_0)   << std::endl;
    }
    if (dim > 1)
    {
        std::cout << "\tDeriv1: " << "\t" << F->L2(Sderiv1, sol1) << "\t" << F->Linf(Sderiv1, sol1) << std::endl;
    }
    if (dim > 2)
    {
        std::cout << "\tDeriv2: " << "\t" << F->L2(Sderiv2, sol2) << "\t" << F->Linf(Sderiv2, sol2) << std::endl;
    }

    std::cout << "\nCalc every cycle \t L2 Error \t Linf Error" << std::endl;
    std::cout << "\tPhys: "   << "\t\t" << F->L2(Ephys, sol) << "\t" << F->Linf(Ephys, sol)   << std::endl;
    if (dim > 0)
    {
        std::cout << "\tDeriv0: " << "\t" << F->L2(Ederiv0, sol0) << "\t" << F->Linf(Ederiv0, sol0) << std::endl;
    }
    if (dim == 1)
    {
        std::cout << "\t2Deriv0: " << "\t" << F->L2(E2deriv0, sol2_0) << "\t" << F->Linf(E2deriv0, sol2_0)   << std::endl;
    }
    if (dim > 1)
    {
        std::cout << "\tDeriv1: " << "\t" << F->L2(Ederiv1, sol1) << "\t" << F->Linf(Ederiv1, sol1) << std::endl;
    }
    if (dim > 2)
    {
        std::cout << "\tDeriv2: " << "\t" << F->L2(Ederiv2, sol2) << "\t" << F->Linf(Ederiv2, sol2) << std::endl;
    }
}
