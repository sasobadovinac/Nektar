///////////////////////////////////////////////////////////////////////////////
//
// File: StdInterpBasis.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Demo for testing functionality of PhysEvaluateBasis
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Timer.h>
#include <iostream>

#include "StdDemoSupport.hpp"

int main(int argc, char *argv[])
{
    DemoSupport demo;
    demo.ParseArguments(argc, argv);
    StdExpansion *E = demo.CreateStdExpansion();

    if (E == nullptr)
    {
        return 1;
    }

    int nCoeffs = E->GetNcoeffs(), totPoints = E->GetTotPoints();
    int nTot = nCoeffs * totPoints, dimension = E->GetShapeDimension();

    Array<OneD, Array<OneD, NekDouble>> coordsE = demo.GetCoords(E);
    Array<OneD, Array<OneD, NekDouble>> coordsIn(totPoints);

    for (int i = 0; i < totPoints; ++i)
    {
        coordsIn[i] = Array<OneD, NekDouble>(dimension, 0.0);
        for (int j = 0; j < dimension; ++j)
        {
            coordsIn[i][j] = coordsE[j][i];
        }
    }

    Array<OneD, Array<OneD, NekDouble>> coordsInmid(3);

    int totCyc = 5;
    std::cout << "Testing new method, bary eval  every cycle" << std::endl;
    LibUtilities::Timer t;
    Array<OneD, Array<OneD, NekDouble>> tempIn(3);
    Array<OneD, Array<OneD, NekDouble>> tempmid(3);
    int ctr                   = 0;
    Array<OneD, NekDouble> qZ = E->GetBasis(0)->GetZ();

    int coordsInmidsz = pow(qZ.size() - 1, 3);
    Array<OneD, Array<OneD, NekDouble>> m_physevalall =
        E->GetPhysEvaluateStorage();

    for (int k = 0; k < 3; k++)
    {
        tempmid[k] = Array<OneD, NekDouble>(coordsInmidsz);
        tempIn[k]  = Array<OneD, NekDouble>(coordsIn.size());
        for (int p = 0; p < coordsIn.size(); p++)
        {
            tempIn[k][p] = coordsIn[p][k];
        }
        ctr++;
    }
    int ctmid = 0;
    for (int k = 0; k < qZ.size() - 1; k++)
    {
        for (int p = 0; p < qZ.size() - 1; p++)
        {
            for (int y = 0; y < qZ.size() - 1; y++)
            {

                tempmid[0][ctmid] = (qZ[k] + qZ[k + 1]) / 2.0;
                tempmid[1][ctmid] = (qZ[p] + qZ[p + 1]) / 2.0;
                tempmid[2][ctmid] = (qZ[y] + qZ[y + 1]) / 2.0;
                ctr++;
                ctmid++;
            }
        }
    }
    // for basis evaluation on quadrature points using Lagrange
    Array<OneD, NekDouble> sol1(nTot);
    // for basis evaluation on midpoint grid using Lagrange
    Array<OneD, NekDouble> sol2(tempmid[0].size() * nCoeffs);
    // for basis evaluation on quadrature points using Bary
    Array<OneD, NekDouble> phys1(nTot);
    // for basis evaluation on midpt grid using Bary
    Array<OneD, NekDouble> phys2(tempmid[0].size() * nCoeffs);
    Array<OneD, NekDouble> ar1(totPoints);

    // evaluation on quadrature points (coulkd just copy from storage)
    for (int i = 0; i < nCoeffs; ++i)
    {
        Vmath::Vcopy(totPoints, &m_physevalall[0][i * totPoints], 1, &ar1[0],
                     1);

        for (int k = 0; k < tempIn[0].size(); ++k)
        {
            Array<OneD, NekDouble> tmpcoord(3);
            tmpcoord[0]                     = tempIn[0][k];
            tmpcoord[1]                     = tempIn[1][k];
            tmpcoord[2]                     = tempIn[2][k];
            phys1[tempIn[0].size() * i + k] = E->PhysEvaluate(tmpcoord, ar1);
        }
    }
    t.Start();
    for (int cyc = 0; cyc < totCyc; ++cyc)

    {
        for (int i = 0; i < nCoeffs; ++i)
        {
            Vmath::Vcopy(totPoints, &m_physevalall[0][i * totPoints], 1,
                         &ar1[0], 1);
            for (int k = 0; k < tempmid[0].size(); ++k)
            {

                Array<OneD, NekDouble> tmpcoord(3);
                tmpcoord[0] = tempmid[0][k];
                tmpcoord[1] = tempmid[1][k];
                tmpcoord[2] = tempmid[2][k];
                // evaluation on midpt grid using barycentric
                phys2[tempmid[0].size() * i + k] =
                    E->PhysEvaluate(tmpcoord, ar1);
            }
        }
    }

    t.Stop();

    std::cout << "Bary method: " << t.TimePerTest(1) / 5.0 << " per cycle ("
              << totCyc << " cycles)." << std::endl;
    // Another approach: Use Nektar++'s approach using interpolation matrices
    for (int cyc = 0; cyc < totCyc; ++cyc)
    {

        for (int i = 0; i < nCoeffs; ++i)
        {
            Vmath::Vcopy(totPoints, &m_physevalall[0][i * totPoints], 1,
                         &ar1[0], 1);

            for (int k = 0; k < tempIn[0].size(); ++k)
            {
                Array<OneD, NekDouble> tmpcoord(3);
                tmpcoord[0] = tempIn[0][k];
                tmpcoord[1] = tempIn[1][k];
                tmpcoord[2] = tempIn[2][k];
                sol1[i * (tempIn[0].size()) + k] =
                    E->PhysEvaluateOld(tmpcoord, ar1);
            }
        }
    }

    t.Start();
    for (int cyc = 0; cyc < totCyc; ++cyc)
    {

        for (int i = 0; i < nCoeffs; ++i)
        {
            Vmath::Vcopy(totPoints, &m_physevalall[0][i * totPoints], 1,
                         &ar1[0], 1);

            for (int k = 0; k < tempmid[0].size(); ++k)
            {
                Array<OneD, NekDouble> tmpcoord(3);
                tmpcoord[0] = tempmid[0][k];
                tmpcoord[1] = tempmid[1][k];
                tmpcoord[2] = tempmid[2][k];
                sol2[i * (tempmid[0].size()) + k] =
                    E->PhysEvaluateOld(tmpcoord, ar1);
            }
        }
    }

    t.Stop();
    std::cout << "Lagrange method: " << t.TimePerTest(totCyc) << " per cycle ("
              << totCyc << " cycles)." << std::endl;

    NekDouble errL2 = 0, errLinf = 0;
    for (int k = 0; k < nCoeffs * (totPoints); k++)
    {
        //      if(sol[k] - phys[k] > 1e-10)
        //{
        //      cout<<" sol :"<<sol1[k]<<" phys:"<<phys1[k]<<"\n\n";
        // exit(0);
        //	}
    }
    cout << "\n midpt interp vals:\n";
    for (int k = 0; k < nCoeffs * (tempmid[0].size()); k++)
    {
        //      if(sol[k] - phys[k] > 1e-10)
        //{
        //    cout<<" sol :"<<sol2[k]<<" phys:"<<phys2[k]<<"\n\n";
        // exit(0);
        //	}
    }

    cout << "L infinity error : " << scientific << errLinf << endl;
    cout << "L 2 error        : " << scientific << errL2 << endl;

    return 0;
}
