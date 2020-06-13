////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessSensor.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Computes gradient of fields.
//
////////////////////////////////////////////////////////////////////////////////

// #include <iostream>
// #include <string>
using namespace std;

// #include <boost/core/ignore_unused.hpp>

// #include <GlobalMapping/Mapping.h>
// #include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessSensor.h"
// #include "ProcessMapping.h"

#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdHexExp.h>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessSensor::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "sensor"),
    ProcessSensor::create,
    "Computes the shock sensor.");

ProcessSensor::ProcessSensor(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["var"] = ConfigOption(false, "0", "Variable to process");
}

ProcessSensor::~ProcessSensor()
{
}

void ProcessSensor::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    int expdim    = m_f->m_graph->GetMeshDimension();
    // int spacedim  = m_f->m_numHomogeneousDir + expdim;
    int nfields   = m_f->m_variables.size();

    int var = m_config["var"].as<int>();

    int nExp = m_f->m_exp[0]->GetExpSize();

    vector<vector<NekDouble>> sensors;
    int maxSize = 0;

    for (int i = 0; i < nExp; ++i)
    {
        int offset = m_f->m_exp[var]->GetPhys_Offset(i);
        LocalRegions::ExpansionSharedPtr Exp = m_f->m_exp[var]->GetExp(i);

        Array<OneD, int> P(expdim);
        Array<OneD, int> numPoints(expdim);
        Array<OneD, LibUtilities::PointsKey> ptsKey(expdim);

        int minP = numeric_limits<int>::max();

        for (int j = 0; j < expdim; ++j)
        {
            P[j]         = Exp->GetBasis(j)->GetNumModes();
            numPoints[j] = Exp->GetBasis(j)->GetNumPoints();
            ptsKey[j]    = LibUtilities::PointsKey (numPoints[j],
                            Exp->GetBasis(j)->GetPointsType());
            
            minP = min(minP, P[j]);
        }

        vector<NekDouble> sensor;

        NekDouble den = 0;
        Array<OneD, NekDouble> phys = m_f->m_exp[var]->GetPhys() + offset;

        for (int k = 1; k < minP; ++k)
        {
            // Declare orthogonal basis.
            StdRegions::StdExpansionSharedPtr OrthoExp;
            switch (Exp->GetGeom()->GetShapeType())
            {
                case LibUtilities::eQuadrilateral:
                {
                    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, P[0] - k,
                                                ptsKey[0]);
                    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, P[1] - k,
                                                ptsKey[1]);
                    OrthoExp = MemoryManager<
                        StdRegions::StdQuadExp>::AllocateSharedPtr(Ba, Bb);
                    break;
                }
                case LibUtilities::eTriangle:
                {
                    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, P[0] - k,
                                                ptsKey[0]);
                    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_B, P[1] - k,
                                                ptsKey[1]);
                    OrthoExp =
                        MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(
                            Ba, Bb);
                    break;
                }
                case LibUtilities::eTetrahedron:
                {
                    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, P[0] - k,
                                                ptsKey[0]);
                    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_B, P[1] - k,
                                                ptsKey[1]);
                    LibUtilities::BasisKey Bc(LibUtilities::eOrtho_C, P[2] - k,
                                                ptsKey[2]);
                    OrthoExp =
                        MemoryManager<StdRegions::StdTetExp>::AllocateSharedPtr(
                            Ba, Bb, Bc);
                    break;
                }
                case LibUtilities::ePyramid:
                {
                    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, P[0] - k,
                                                ptsKey[0]);
                    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, P[1] - k,
                                                ptsKey[1]);
                    LibUtilities::BasisKey Bc(LibUtilities::eOrtho_C, P[2] - k,
                                                ptsKey[2]);
                    OrthoExp =
                        MemoryManager<StdRegions::StdPyrExp>::AllocateSharedPtr(
                            Ba, Bb, Bc);
                    break;
                }
                case LibUtilities::ePrism:
                {
                    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, P[0] - k,
                                                ptsKey[0]);
                    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, P[1] - k,
                                                ptsKey[1]);
                    LibUtilities::BasisKey Bc(LibUtilities::eOrtho_B, P[2] - k,
                                                ptsKey[2]);
                    OrthoExp =
                        MemoryManager<StdRegions::StdPrismExp>::AllocateSharedPtr(
                            Ba, Bb, Bc);
                    break;
                }
                case LibUtilities::eHexahedron:
                {
                    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, P[0] - k,
                                                ptsKey[0]);
                    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, P[1] - k,
                                                ptsKey[1]);
                    LibUtilities::BasisKey Bc(LibUtilities::eOrtho_A, P[2] - k,
                                                ptsKey[2]);
                    OrthoExp =
                        MemoryManager<StdRegions::StdHexExp>::AllocateSharedPtr(
                            Ba, Bb, Bc);
                    break;
                }
                default:
                    ASSERTL0(false, "Shape not supported.");
                    break;
            }

            int nq = OrthoExp->GetTotPoints();

            Array<OneD, NekDouble> coeffs      = Array<OneD, NekDouble>(OrthoExp->GetNcoeffs());
            Array<OneD, NekDouble> physReduced = Array<OneD, NekDouble>(OrthoExp->GetTotPoints());
            Array<OneD, NekDouble> tmpArray    = Array<OneD, NekDouble>(OrthoExp->GetTotPoints(), 0.0);

            // Project solution to lower order
            OrthoExp->FwdTrans(phys, coeffs);
            OrthoExp->BwdTrans(coeffs, physReduced);

            // ||phys-physReduced||^2
            Vmath::Vsub(nq, phys,     1, physReduced, 1, tmpArray, 1);
            Vmath::Vmul(nq, tmpArray, 1, tmpArray,    1, tmpArray, 1);
            NekDouble num = Exp->Integral(tmpArray);

            if (k == 1)
            {
                // ||phys||^2
                Vmath::Vmul(nq, phys, 1, phys, 1, tmpArray, 1);
                den = Exp->Integral(tmpArray);
            }

            // Calculate error =||phys-physReduced||^2 / ||phys||^2
            sensor.push_back(fabs(num / den));

            if (k + 1 == minP)
            {
                // ||physReduced||^2
                Vmath::Vmul(nq, physReduced, 1, physReduced, 1, tmpArray, 1);
                num = Exp->Integral(tmpArray);

                // ||physReduced||^2 / ||phys||^2
                sensor.push_back(fabs(num / den));
            }
            else
            {
                phys = physReduced;
            }
        }

        maxSize = max(maxSize, int(sensor.size()));
        sensors.push_back(sensor);
    }

    m_f->m_exp.resize(nfields + maxSize);

    for (int i = 0; i < maxSize; ++i)
    {
        m_f->m_variables.push_back("Sensor" + to_string(i));
        m_f->m_exp[nfields + i] = m_f->AppendExpList(0);
    }

    for (int i = 0; i < sensors.size(); ++i)
    {
        int npts = m_f->m_exp[var]->GetExp(i)->GetTotPoints();
        int offset = m_f->m_exp[var]->GetPhys_Offset(i);

        for (int j = 0; j < sensors[i].size(); j ++)
        {
            // Add sensor to all points in the element
            Array<OneD, NekDouble> tmp = m_f->m_exp[nfields + j]->UpdatePhys() + offset;
            Vmath::Fill(npts, sensors[i][j], tmp, 1);
        }
    }
}
}
}
