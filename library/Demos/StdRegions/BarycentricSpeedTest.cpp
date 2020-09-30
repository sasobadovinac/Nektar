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
            coordIn[i][j] = coordsE[j][i];
        }
    }


    // Do tests for 1000 cycles
    int totCyc = 1000;
    // Calc interpolation matrix every call
    LibUtilities::Timer t;

    NekDouble phys, deriv0, deriv1, deriv2;

    Array<OneD, NekDouble> physDeriv0(totPoints), physDeriv1(totPoints), physDeriv2(totPoints);
    std::cout << "Testing old method, calc deriv matrix every cycle" << std::endl;
    t.Start();
    for (int cyc = 0; cyc < totCyc; ++cyc)
    {
        E->PhysDeriv(physIn, physDeriv0, physDeriv1, physDeriv2);
        for (int i = 0; i < totPoints; ++i)
        {
            phys   = E->PhysEvaluateOld(coordIn[i], physIn);
            deriv0 = E->PhysEvaluateOld(coordIn[i], physDeriv0);
            deriv1 = E->PhysEvaluateOld(coordIn[i], physDeriv1);
            deriv2 = E->PhysEvaluateOld(coordIn[i], physDeriv1);
        }
    }
    t.Stop();
    std::cout << "Old method: " << t.TimePerTest(totCyc) << " per cycle (" << totCyc << " cycles)." << std::endl;

    std::cout << "Testing old method, store deriv matrix" << std::endl;
    t.Start();
    E->PhysDeriv(physIn, physDeriv0, physDeriv1, physDeriv2);
    for (int cyc = 0; cyc < totCyc; ++cyc)
    {
        for (int i = 0; i < totPoints; ++i)
        {
            phys   = E->PhysEvaluateOld(coordIn[i], physIn);
            deriv0 = E->PhysEvaluateOld(coordIn[i], physDeriv0);
            deriv1 = E->PhysEvaluateOld(coordIn[i], physDeriv1);
            deriv2 = E->PhysEvaluateOld(coordIn[i], physDeriv1);
        }
    }
    t.Stop();
    std::cout << "Old method: " << t.TimePerTest(totCyc) << " per cycle (" << totCyc << " cycles)." << std::endl;


    std::cout << "Testing barycentric method" << std::endl;
    t.Start();
    for (int cyc = 0; cyc < totCyc; ++cyc)
    {
        for (int i = 0; i < totPoints; ++i)
        {
            phys = E->PhysEvaluate(coordIn[i], physIn, deriv0, deriv1, deriv2);
        }
    }
    t.Stop();
    std::cout << "New method: " << t.TimePerTest(totCyc) << " per cycle (" << totCyc << " cycles)." << std::endl;


    // Iterative calc interpolation matrix


    // Barycentric interpolation

}