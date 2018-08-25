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

NekDouble ProcessCrossField::AdamsBashforth_coeffs[4][4] = {
    {1.0, 0.0, 0.0, 0.0},
    {3.0 / 2.0, -1.0 / 2.0, 0.0, 0.0},
    {23.0 / 12.0, -4.0 / 3.0, 5.0 / 12.0, 0.0},
    {55.0 / 24.0, -59.0 / 24.0, 37.0 / 24.0, -3.0 / 8.0}};

ProcessCrossField::ProcessCrossField(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["outfile"] = ConfigOption(false, "", "Output filename.");
}

ProcessCrossField::~ProcessCrossField()
{
}

void ProcessCrossField::Process(po::variables_map &vm)
{
    ofstream file;
    file.open(m_config["outfile"].as<string>());

    int dim = m_f->m_graph->GetSpaceDimension();

    ASSERTL0(dim == 2, "Cross field post-processing only possible in 2D.")

    // Find all elements containing u = 0 and v = 0
    vector<set<int>> isoElmts(dim);

    for (int i = 0; i < dim; ++i)
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
                    Array<OneD, NekDouble> lcoords(dim);
                    lcoords[0] = 0.0;
                    lcoords[1] = 1.0;

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

    // Store singularities and their starting element
    map<int, pair<int, Array<OneD, NekDouble>>> singularities;

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
        for (int i = 0; i < dim; ++i)
        {
            Vmath::Vabs(npoints, m_f->m_exp[i]->GetPhys() + offset, 1, tmp, 1);
            Vmath::Vadd(npoints, costFun, 1, tmp, 1, costFun, 1);
        }

        auto point = min_element(costFun.begin(), costFun.end());

        Array<OneD, Array<OneD, NekDouble>> quadEta(2);
        for (int i = 0; i < dim; ++i)
        {
            quadEta[i] = Array<OneD, NekDouble>(npoints);
        }
        expansion->StdExpansion::GetCoords(quadEta[0], quadEta[1]);

        // Starting point for Newton's method
        Array<OneD, NekDouble> eta(dim, 0.0);
        for (int i = 0; i < dim; ++i)
        {
            eta[i] = quadEta[i][point - costFun.begin()];
        }
        */

        // Starting point for Newton's method
        Array<OneD, NekDouble> eta(dim, expansion->GetGeom()->GetShapeType() ==
                                                LibUtilities::eTriangle
                                            ? -1.0 / 3.0
                                            : 0.0);

        // Variables for iterations
        Array<OneD, NekDouble> u(dim);
        Array<OneD, NekDouble> u_eta(pow(dim, 2));
        Array<OneD, NekDouble> delta(dim);

        // Hold u_eta at quadrature points
        Array<OneD, Array<OneD, NekDouble>> deriv(pow(dim, 2));
        for (int i = 0; i < pow(dim, 2); ++i)
        {
            deriv[i] = Array<OneD, NekDouble>(npoints);
        }
        for (int i = 0; i < dim; ++i)
        {
            expansion->StdPhysDeriv(m_f->m_exp[i]->GetPhys() + offset,
                                    deriv[dim * i], deriv[dim * i + 1]);
        }

        // Newton's method to find point where u = v = 0
        do
        {
            // Evaluate u and u_eta at current point
            for (int i = 0; i < dim; ++i)
            {
                u[i] = expansion->StdPhysEvaluate(
                    eta, m_f->m_exp[i]->GetPhys() + offset);
            }
            for (int i = 0; i < pow(dim, 2); ++i)
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

            // Invert the 2x2 jacobian
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

        } while (sqrt(pow(delta[0], 2) + pow(delta[1], 2)) > 1.0e-9);

        // Check if point is still inside the element
        // If not, we dismiss it;
        // it's most likely been detected by adjacent element
        switch (expansion->GetGeom()->GetShapeType())
        {
            case LibUtilities::eTriangle:
                if (!(eta[0] >= -1.0 && eta[1] >= -1.0 &&
                      eta[0] + eta[1] <= 0.0))
                {
                    continue;
                }
                break;

            case LibUtilities::eQuadrilateral:
                if (!(eta[0] >= -1.0 && eta[1] >= -1.0 && eta[0] <= 1.0 &&
                      eta[1] <= 1.0))
                {
                    continue;
                }
                break;

            default:
                ASSERTL0(false, "Unexpected shape type");
        }

        // Determine the number of branches
        // In theory, we should be doing:
        // \( nbranches = 4 - \int_\pi {\psi'(\theta)} \)
        // In practice, we check for jumps of atan2(v, u) at +-\pi
        int ncquads    = 100;
        NekDouble dist = 1.0e-6;
        Array<OneD, NekDouble> cquad(dim);
        int nbranches  = 4;
        NekDouble oldV = 0.0;

        for (int i = 0; i <= ncquads; ++i)
        {
            cquad[0] = eta[0] + dist * cos(2.0 * M_PI * i / ncquads);
            cquad[1] = eta[1] + dist * sin(2.0 * M_PI * i / ncquads);

            NekDouble u = expansion->StdPhysEvaluate(
                cquad, m_f->m_exp[0]->GetPhys() + offset);
            NekDouble v = expansion->StdPhysEvaluate(
                cquad, m_f->m_exp[1]->GetPhys() + offset);

            // Condition: u < 0 and v changes sign, i.e. +-\pi
            if (u < 0.0)
            {
                if (oldV * v < 0.0)
                {
                    if (v < 0.0)
                    {
                        --nbranches;
                    }
                    else
                    {
                        ++nbranches;
                    }
                }
            }

            oldV = v;
        }

        singularities[elmt] = make_pair(nbranches, eta);

        Array<OneD, NekDouble> x(dim);
        expansion->GetCoord(eta, x);

        cout << "Singularity #" << singularities.size() << " with " << nbranches
             << " branches at (" << x[0] << ", " << x[1] << ")" << endl;

        file << x[0] << "," << x[1] << endl;
    }

    for (auto &s : singularities)
    {
        int dim = m_f->m_graph->GetSpaceDimension();

        // I can make dist the first stage of RK based on |U|=1
        NekDouble dist = 1.0e-2;

        int elmt      = s.first;
        int nbranches = s.second.first;
        Array<OneD, NekDouble> eta = s.second.second;

        Array<OneD, NekDouble> x(dim);
        m_f->m_exp[0]->GetExp(elmt)->GetCoord(eta, x);

        vector<vector<pair<Array<OneD, NekDouble>, NekDouble>>> fullHist;

        for (int i = 0; i < nbranches; ++i)
        {
            vector<pair<Array<OneD, NekDouble>, NekDouble>> history;

            NekDouble dir;
            if (i == 0)
            {
                dir = 0.0;
            }
            else
            {
                dir = fullHist.back().front().second + 2.0 * M_PI / nbranches;
            }

            NekDouble delta;
            Array<OneD, NekDouble> point(dim);
            Array<OneD, NekDouble> lpoint(dim);
            int rot;
            bool uneg, vneg;

            do
            {
                point[0] = x[0] + dist * cos(dir);
                point[1] = x[1] + dist * sin(dir);

                int id     = m_f->m_exp[0]->GetExpIndex(point, lpoint);
                int offset = m_f->m_exp[0]->GetPhys_Offset(id);

                NekDouble u = m_f->m_exp[0]->GetExp(id)->StdPhysEvaluate(
                    lpoint, m_f->m_exp[0]->GetPhys() + offset);
                NekDouble v = m_f->m_exp[0]->GetExp(id)->StdPhysEvaluate(
                    lpoint, m_f->m_exp[1]->GetPhys() + offset);

                NekDouble psi = atan2(v, u) / 4.0;

                delta =
                    dir -
                    (psi + round((dir - psi) / (M_PI / 2.0)) * (M_PI / 2.0));
                dir -= delta;

                rot  = round((dir - psi) / (M_PI / 2.0));
                uneg = u < 0.0;
                vneg = v < 0.0;

            } while (fabs(delta) < 1.0e-9);

            file << point[0] << "," << point[1] << endl;

            history.push_back(make_pair(x, dir));     // singularity itself
            history.push_back(make_pair(point, dir)); // first point

            while (true)
            {
                int numIte = history.size();
                int order  = min(numIte, 4);

                Array<OneD, NekDouble> newPoint(2);
                newPoint[0] = history.back().first[0];
                newPoint[1] = history.back().first[1];

                for (int j = 0; j < order; ++j)
                {
                    newPoint[0] += dist * AdamsBashforth_coeffs[order - 1][j] *
                                   cos(history[history.size() - 1 - j].second);
                    newPoint[1] += dist * AdamsBashforth_coeffs[order - 1][j] *
                                   sin(history[history.size() - 1 - j].second);
                }

                file << newPoint[0] << "," << newPoint[1] << endl;

                int id = m_f->m_exp[0]->GetExpIndex(newPoint, lpoint);

                if (id == -1)
                {
                    break;
                }

                int offset = m_f->m_exp[0]->GetPhys_Offset(id);

                NekDouble u = m_f->m_exp[0]->GetExp(id)->StdPhysEvaluate(
                    lpoint, m_f->m_exp[0]->GetPhys() + offset);
                NekDouble v = m_f->m_exp[0]->GetExp(id)->StdPhysEvaluate(
                    lpoint, m_f->m_exp[1]->GetPhys() + offset);

                if (uneg)
                {
                    if (vneg && v > 0.0)
                    {
                        --rot;
                    }
                    else if (!vneg && v < 0.0)
                    {
                        ++rot;
                    }
                }

                uneg = u < 0.0;
                vneg = v < 0.0;

                history.push_back(
                    make_pair(newPoint, atan2(v, u) / 4.0 + rot * M_PI / 2.0));
            }

            fullHist.push_back(history);
        }
    }

    file.close();
}
}
}
