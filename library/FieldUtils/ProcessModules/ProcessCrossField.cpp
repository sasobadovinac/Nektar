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

ProcessCrossField::ProcessCrossField(FieldSharedPtr f) : ProcessGrad(f)
{
}

ProcessCrossField::~ProcessCrossField()
{
}

void ProcessCrossField::Process(po::variables_map &vm)
{
    ASSERTL0(m_f->m_graph->GetSpaceDimension() == 2,
             "Cross field post-processing only possible in 2D.")

    // Calculate the Jacobian of (u, v)
    ProcessGrad::Process(vm);

    // Find all elements containing u = 0 and v = 0
    vector<set<LocalRegions::ExpansionSharedPtr>> isoElmts(2);

    for (int i = 0; i < 2; ++i)
    {
        Array<OneD, NekDouble> phys = m_f->m_exp[i]->GetPhys();

        for (int elmt = 0; elmt < m_f->m_exp[i]->GetNumElmts(); ++elmt)
        {
            int offset  = m_f->m_exp[i]->GetPhys_Offset(elmt);
            int npoints = m_f->m_exp[i]->GetTotPoints(elmt);

            auto mm = minmax_element(phys.begin() + offset,
                                     phys.begin() + offset + npoints);

            if (*(mm.first) <= 0.0 && *(mm.second) >= 0.0)
            {
                isoElmts[i].insert(m_f->m_exp[i]->GetExp(elmt));
            }
            else
            {
                // TODO: check collapsed vertex
            }
        }
    }

    // Find intersection of u = 0 and v = 0 element sets
    vector<LocalRegions::ExpansionSharedPtr> intElmts;
    set_intersection(isoElmts[0].begin(), isoElmts[0].end(),
                     isoElmts[1].begin(), isoElmts[1].end(),
                     back_inserter(intElmts));

    // Find exact location of each singularity
    // We assume there is at most 1 singularity per element
    for (auto &elmt : intElmts)
    {
        int id = elmt->GetElmtId();

        int offset  = m_f->m_exp[0]->GetPhys_Offset(id);
        int npoints = m_f->m_exp[0]->GetTotPoints(id);

        // Find nearest quadrature point to u = v = 0
        Array<OneD, NekDouble> costFun(npoints, 0.0);
        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < npoints; ++j)
            {
                costFun[j] += fabs(m_f->m_exp[i]->GetPhys()[offset + j]);
            }
        }

        auto point = min_element(costFun.begin(), costFun.end());

        Array<OneD, NekDouble> x(npoints);
        Array<OneD, NekDouble> y(npoints);
        elmt->GetCoords(x, y);

        // Starting point for Newton's method
        Array<OneD, NekDouble> coords(3, 0.0);
        coords[0] = x[point - costFun.begin()];
        coords[1] = y[point - costFun.begin()];

        Array<OneD, NekDouble> vals(m_f->m_exp.size());
        Array<OneD, NekDouble> delta(2, numeric_limits<NekDouble>::max());

        // Copy of local phys for evaluation at point
        Array<OneD, Array<OneD, NekDouble>> locPhys(m_f->m_exp.size());
        for (int i = 0; i < m_f->m_exp.size(); ++i)
        {
            locPhys[i] = Array<OneD, NekDouble>(npoints);
            Vmath::Vcopy(npoints, &m_f->m_exp[i]->GetPhys()[offset], 1,
                         &locPhys[i][0], 1);
        }

        // Newton's method to find point where u = v = 0
        while (sqrt(delta[0] * delta[0] + delta[1] * delta[1]) > 1.0e-6)
        {
            // Evaluate all fields at current point
            for (int i = 0; i < m_f->m_exp.size(); ++i)
            {
                vals[i] = elmt->PhysEvaluate(coords, locPhys[i]);
            }

            // vals
            // 0 -> u
            // 1 -> v
            // 2 -> u_x
            // 3 -> u_y
            // 4 -> v_x
            // 5 -> v_y

            // Invert the jacobian
            iter_swap(vals.begin() + 2, vals.begin() + 5);
            Vmath::Neg(2, &vals[3], 1);
            Vmath::Smul(4, 1.0 / (vals[2] * vals[5] - vals[3] * vals[4]),
                        &vals[2], 1, &vals[2], 1);

            // vals
            // 2 -> x_u
            // 3 -> x_v
            // 4 -> y_u
            // 5 -> y_v

            // J^(-1) * u(x_n)
            delta[0] = vals[2] * vals[0] + vals[3] * vals[1];
            delta[1] = vals[4] * vals[0] + vals[5] * vals[1];

            // x_(n+1) = x_n - J^(-1) * u(x_n)
            coords[0] -= delta[0];
            coords[1] -= delta[1];
        }

        // Check if point is inside element
        // If not, we dismiss it
        // It's most likely been detected by adjacent element
        if (m_f->m_exp[0]->GetExpIndex(coords) != id)
        {
            continue;
        }

        cout << coords[0] << "\t" << coords[1] << endl;
    }
}
}
}
