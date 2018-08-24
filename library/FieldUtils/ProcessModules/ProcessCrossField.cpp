////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessCrossField.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
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
//  Description: Post-processes cross field simulation results.
//
////////////////////////////////////////////////////////////////////////////////

using namespace std;

#include "ProcessCrossField.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessCrossField::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "crossfield"), ProcessCrossField::create,
        "Post-processes cross field simulation results.");

ProcessCrossField::ProcessCrossField(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessCrossField::~ProcessCrossField()
{
}

void ProcessCrossField::Process(po::variables_map &vm)
{
    ASSERTL0(m_f->m_graph->GetSpaceDimension() == 2,
             "Cross field post-processing only possible in 2D.")

    // Find all elements containing u = 0 and v = 0
    vector<set<int>> isoElmts(2);

    for (int i = 0; i < 2; ++i)
    {
        Array<OneD, NekDouble> phys = m_f->m_exp[i]->GetPhys();

        for (int elmt = 0; elmt < m_f->m_exp[i]->GetNumElmts(); ++elmt)
        {
            int offset  = m_f->m_exp[i]->GetPhys_Offset(elmt);
            int npoints = m_f->m_exp[i]->GetTotPoints(elmt);

            auto mm = minmax_element(phys.begin() + offset,
                                     phys.begin() + offset + npoints);

            LocalRegions::ExpansionSharedPtr expansion =
                m_f->m_exp[i]->GetExp(elmt);

            NekDouble minV = *(mm.first);
            NekDouble maxV = *(mm.second);

            // Check value at collapsed vertex of triangles
            if (expansion->DetShapeType() == LibUtilities::eTriangle)
            {
                if (minV > 0.0 || maxV < 0.0)
                {
                    // Is it possible to do aggregate initialisation?
                    Array<OneD, NekDouble> lcoords(3);
                    lcoords[0] = 0.0;
                    lcoords[1] = 1.0;
                    lcoords[2] = 0.0;

                    NekDouble colVal = expansion->StdPhysEvaluate(
                        lcoords, m_f->m_exp[i]->GetPhys() + offset);

                    minV = min(minV, colVal);
                    maxV = max(maxV, colVal);
                }
            }

            if (minV <= 0.0 && maxV >= 0.0)
            {
                isoElmts[i].insert(elmt);
            }
        }
    }

    // Find intersection of u = 0 and v = 0 element sets
    vector<int> intElmts;
    set_intersection(isoElmts[0].begin(), isoElmts[0].end(),
                     isoElmts[1].begin(), isoElmts[1].end(),
                     back_inserter(intElmts));

    // Count the number of singularities
    int cnt = 0;

    // Find exact location of each singularity
    // We assume there is at most 1 singularity per element
    for (auto &elmt : intElmts)
    {
        LocalRegions::ExpansionSharedPtr expansion =
            m_f->m_exp[0]->GetExp(elmt);

        int offset  = m_f->m_exp[0]->GetPhys_Offset(elmt);
        int npoints = m_f->m_exp[0]->GetTotPoints(elmt);

        /*
        // This commented part is not guaranteed to work if starting point is on
        // the boundary
        // I might need to re-use it if I have to find multiple signularities
        // inside a single element

        // Find nearest quadrature point to u = v = 0
        Array<OneD, NekDouble> costFun(npoints, 0.0);
        Array<OneD, NekDouble> tmp(npoints);
        for (int i = 0; i < 2; ++i)
        {
            Vmath::Vabs(npoints, m_f->m_exp[i]->GetPhys() + offset, 1, tmp, 1);
            Vmath::Vadd(npoints, costFun, 1, tmp, 1, costFun, 1);
        }

        auto point = min_element(costFun.begin(), costFun.end());

        Array<OneD, Array<OneD, NekDouble>> quadEta(2);
        for (int i = 0; i < 2; ++i)
        {
            quadEta[i] = Array<OneD, NekDouble>(npoints);
        }
        expansion->StdExpansion::GetCoords(quadEta[0], quadEta[1]);

        // Starting point for Newton's method
        Array<OneD, NekDouble> eta(2, 0.0);
        eta[0] = quadEta[0][point - costFun.begin()];
        eta[1] = quadEta[1][point - costFun.begin()];
        */

        // Starting point for Newton's method
        Array<OneD, NekDouble> eta(2, expansion->GetGeom()->GetShapeType() ==
                                              LibUtilities::eTriangle
                                          ? -1.0 / 3.0
                                          : 0.0);

        // Variables for iterations
        Array<OneD, NekDouble> u(2);
        Array<OneD, NekDouble> u_eta(4);
        Array<OneD, NekDouble> delta(2, numeric_limits<NekDouble>::max());

        // Hold u_eta at quadrature points
        Array<OneD, Array<OneD, NekDouble>> deriv(4);
        for (int i = 0; i < 4; ++i)
        {
            deriv[i] = Array<OneD, NekDouble>(npoints);
        }
        for (int i = 0; i < 2; ++i)
        {
            expansion->StdPhysDeriv(m_f->m_exp[i]->GetPhys() + offset,
                                    deriv[2 * i], deriv[2 * i + 1]);
        }

        // Newton's method to find point where u = v = 0
        while (sqrt(pow(delta[0], 2) + pow(delta[1], 2)) > 1.0e-9)
        {
            // Evaluate u and u_eta at current point
            for (int i = 0; i < 2; ++i)
            {
                u[i] = expansion->StdPhysEvaluate(
                    eta, m_f->m_exp[i]->GetPhys() + offset);
            }
            for (int i = 0; i < 4; ++i)
            {
                u_eta[i] = expansion->StdPhysEvaluate(eta, deriv[i]);
            }

            // u
            // 0 -> u
            // 1 -> v

            // u_eta
            // 0 -> u_eta0
            // 1 -> u_eta1
            // 2 -> v_eta0
            // 3 -> v_eta1

            // Invert the jacobian
            iter_swap(u_eta.begin() + 0, u_eta.begin() + 3);
            Vmath::Neg(2, &u_eta[1], 1);
            Vmath::Smul(4, 1.0 / (u_eta[0] * u_eta[3] - u_eta[1] * u_eta[2]),
                        &u_eta[0], 1, &u_eta[0], 1);

            // u_eta
            // 2 -> eta0_u
            // 3 -> eta0_v
            // 4 -> eta0_u
            // 5 -> eta0_v

            // delta = J^(-1) * u(x_n)
            delta[0] = u_eta[0] * u[0] + u_eta[1] * u[1];
            delta[1] = u_eta[2] * u[0] + u_eta[3] * u[1];

            // x_(n+1) = x_n - delta
            eta[0] -= delta[0];
            eta[1] -= delta[1];
        }

        NekDouble tol = 0.0;

        // Check if point is still inside the element
        // If not, we dismiss it;
        // it's most likely been detected by adjacent element
        switch (expansion->GetGeom()->GetShapeType())
        {
            case LibUtilities::eTriangle:
                if (!(eta[0] >= -(1.0 + tol) && eta[1] >= -(1.0 + tol) &&
                      eta[0] + eta[1] <= tol))
                {
                    // cout << "Outside!\t\t" << eta[0] << "\t\t" << eta[1] <<
                    // endl;
                    // cout << "Element #" << elmt << endl;
                    continue;
                }
                break;

            case LibUtilities::eQuadrilateral:
                if (!(eta[0] >= -(1.0 + tol) && eta[1] >= -(1.0 + tol) &&
                      eta[0] <= (1.0 + tol) && eta[1] <= (1.0 + tol)))
                {
                    continue;
                }
                break;

            default:
                ASSERTL0(false, "Unexpected shape type");
        }

        Array<OneD, NekDouble> x(2);
        expansion->GetCoord(eta, x);

        cout << ++cnt << "\t\t" << x[0] << "\t\t" << x[1] << endl;
    }
}
}
}
