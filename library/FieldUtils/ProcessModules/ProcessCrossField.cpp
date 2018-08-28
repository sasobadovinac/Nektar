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
    m_config["outcsv"] = ConfigOption(false, "", "Output filename.");
    m_config["step"]   = ConfigOption(false, "", "Streamline step size.");
}

ProcessCrossField::~ProcessCrossField()
{
}

void ProcessCrossField::Process(po::variables_map &vm)
{
    m_dim = m_f->m_graph->GetSpaceDimension();
    ASSERTL0(m_dim == 2, "Cross field post-processing only possible in 2D.")

    ASSERTL0(m_config["outcsv"].m_beenSet, "An output CSV file must be set.")
    m_csvfile.open(m_config["outcsv"].as<string>());

    vector<pair<Array<OneD, NekDouble>, int>> singularities =
        FindAllSingularities();

    m_step =
        m_config["step"].m_beenSet ? m_config["step"].as<NekDouble>() : 1.0e-2;

    for (auto &s : singularities)
    {
        Array<OneD, NekDouble> x = s.first;
        int nbranches = s.second;

        NekDouble dir = 0.0;

        for (int i = 0; i < nbranches; ++i)
        {
            Streamline sl = Streamline(m_f, x, m_step, m_csvfile);
            sl.Initialise(dir);

            bool cont = true;

            while (cont)
            {
                cont = cont && sl.Advance();
            }

            dir += 2.0 * M_PI / nbranches;
        }
    }

    m_csvfile.close();
}

vector<pair<Array<OneD, NekDouble>, int>> ProcessCrossField::
    FindAllSingularities()
{
    vector<set<int>> isoElmts = FindIsocontourElements();

    // Find intersection of isocontour elements
    vector<int> intElmts;
    set_intersection(isoElmts[0].begin(), isoElmts[0].end(),
                     isoElmts[1].begin(), isoElmts[1].end(),
                     back_inserter(intElmts));

    // Store singularities and their number of branches
    vector<pair<Array<OneD, NekDouble>, int>> ret;

    // Find exact location of each singularity
    // We assume there is at most 1 singularity per element
    for (auto &elmt : intElmts)
    {
        Array<OneD, NekDouble> eta = FindSingularityInElmt(elmt);
        if (!eta.num_elements())
        {
            // No singularity found
            continue;
        }

        int nbranches = CalculateNumberOfBranches(elmt, eta);

        Array<OneD, NekDouble> x(m_dim);
        m_f->m_exp[0]->GetExp(elmt)->GetCoord(eta, x);

        ret.push_back(make_pair(x, nbranches));

        cout << "Singularity #" << ret.size() << " with " << nbranches
             << " branches at (" << x[0] << ", " << x[1] << ")" << endl;

        m_csvfile << x[0] << "," << x[1] << endl;
    }

    return ret;
}

vector<set<int>> ProcessCrossField::FindIsocontourElements()
{
    vector<set<int>> ret(m_dim);

    for (int i = 0; i < m_dim; ++i)
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
                    Array<OneD, NekDouble> lcoords(m_dim);
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
                ret[i].insert(elmt);
            }
        }
    }

    return ret;
}

Array<OneD, NekDouble> ProcessCrossField::FindSingularityInElmt(int id)
{
    LocalRegions::ExpansionSharedPtr expansion = m_f->m_exp[0]->GetExp(id);

    int offset  = m_f->m_exp[0]->GetPhys_Offset(id);
    int npoints = m_f->m_exp[0]->GetTotPoints(id);

    /*
    // This commented part is not guaranteed to work if starting point is on
    // the boundary
    // I might need to re-use it if I have to find multiple signularities
    // inside a single element

    // Find nearest quadrature point to u = v = 0
    Array<OneD, NekDouble> costFun(npoints, 0.0);
    Array<OneD, NekDouble> tmp(npoints);
    for (int i = 0; i < m_dim; ++i)
    {
        Vmath::Vabs(npoints, m_f->m_exp[i]->GetPhys() + offset, 1, tmp, 1);
        Vmath::Vadd(npoints, costFun, 1, tmp, 1, costFun, 1);
    }

    auto point = min_element(costFun.begin(), costFun.end());

    Array<OneD, Array<OneD, NekDouble>> quadEta(2);
    for (int i = 0; i < m_dim; ++i)
    {
        quadEta[i] = Array<OneD, NekDouble>(npoints);
    }
    expansion->StdExpansion::GetCoords(quadEta[0], quadEta[1]);

    // Starting point for Newton's method
    Array<OneD, NekDouble> eta(m_dim, 0.0);
    for (int i = 0; i < m_dim; ++i)
    {
        eta[i] = quadEta[i][point - costFun.begin()];
    }
    */

    // Starting point for Newton's method
    Array<OneD, NekDouble> eta(m_dim, expansion->GetGeom()->GetShapeType() ==
                                              LibUtilities::eTriangle
                                          ? -1.0 / 3.0
                                          : 0.0);

    // Variables for iterations
    Array<OneD, NekDouble> u(m_dim);
    Array<OneD, NekDouble> u_eta(pow(m_dim, 2));
    Array<OneD, NekDouble> delta(m_dim);

    // Hold u_eta at quadrature points
    Array<OneD, Array<OneD, NekDouble>> deriv(pow(m_dim, 2));
    for (int i = 0; i < pow(m_dim, 2); ++i)
    {
        deriv[i] = Array<OneD, NekDouble>(npoints);
    }
    for (int i = 0; i < m_dim; ++i)
    {
        expansion->StdPhysDeriv(m_f->m_exp[i]->GetPhys() + offset,
                                deriv[m_dim * i], deriv[m_dim * i + 1]);
    }

    // Newton's method to find point where u = v = 0
    do
    {
        // Evaluate u and u_eta at current point
        for (int i = 0; i < m_dim; ++i)
        {
            u[i] = expansion->StdPhysEvaluate(eta, m_f->m_exp[i]->GetPhys() +
                                                       offset);
        }
        for (int i = 0; i < pow(m_dim, 2); ++i)
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
            if (!(eta[0] >= -1.0 && eta[1] >= -1.0 && eta[0] + eta[1] <= 0.0))
            {
                eta = Array<OneD, NekDouble>(0);
                return eta;
            }
            break;

        case LibUtilities::eQuadrilateral:
            if (!(eta[0] >= -1.0 && eta[1] >= -1.0 && eta[0] <= 1.0 &&
                  eta[1] <= 1.0))
            {
                eta = Array<OneD, NekDouble>(0);
                return eta;
            }
            break;

        default:
            ASSERTL0(false, "Unexpected shape type");
    }

    return eta;
}

int ProcessCrossField::CalculateNumberOfBranches(int id,
                                                 Array<OneD, NekDouble> eta)
{
    // Determine the number of branches
    // In theory, we should be doing:
    // \( nbranches = 4 - \int_\pi {\psi'(\theta)} \)
    // In practice, we check for jumps of atan2(v, u) at +-\pi
    int ncquads    = 100;
    NekDouble dist = 1.0e-6;
    Array<OneD, NekDouble> cquad(m_dim);
    int nbranches  = 4;
    NekDouble oldV = 0.0;

    LocalRegions::ExpansionSharedPtr expansion = m_f->m_exp[0]->GetExp(id);
    int offset = m_f->m_exp[0]->GetPhys_Offset(id);

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

    return nbranches;
}

NekDouble Streamline::AdamsBashforth_coeffs[4][4] = {
    {1.0, 0.0, 0.0, 0.0},
    {3.0 / 2.0, -1.0 / 2.0, 0.0, 0.0},
    {23.0 / 12.0, -4.0 / 3.0, 5.0 / 12.0, 0.0},
    {55.0 / 24.0, -59.0 / 24.0, 37.0 / 24.0, -3.0 / 8.0}};

void Streamline::Initialise(NekDouble &angle)
{
    ASSERTL0(m_points.size() == 1, "Streamline cannot be intialised twice.")

    NekDouble delta;
    Array<OneD, NekDouble> point(m_dim);
    Array<OneD, NekDouble> lpoint(m_dim);
    NekDouble u, v;

    do
    {
        // Maybe use a smaller value than m_step?
        point[0] = m_points[0][0] + m_step * cos(angle);
        point[1] = m_points[0][1] + m_step * sin(angle);

        int id     = m_f->m_exp[0]->GetExpIndex(point, lpoint);
        int offset = m_f->m_exp[0]->GetPhys_Offset(id);

        u = m_f->m_exp[0]->GetExp(id)->StdPhysEvaluate(
            lpoint, m_f->m_exp[0]->GetPhys() + offset);
        v = m_f->m_exp[0]->GetExp(id)->StdPhysEvaluate(
            lpoint, m_f->m_exp[1]->GetPhys() + offset);

        NekDouble psi = atan2(v, u) / 4.0;

        m_rot = round((angle - psi) / (M_PI / 2.0));

        delta = angle - (psi + m_rot * (M_PI / 2.0));
        angle -= delta;

    } while (fabs(delta) < 1.0e-9);

    // Over-write angle of singularity
    m_angles.back() = angle;

    // First point
    AddPoint(point, angle);

    m_neg.first  = u < 0.0;
    m_neg.second = v < 0.0;
}

bool Streamline::Advance()
{
    int numPoints = m_points.size();
    int order     = min(numPoints, 4);

    Array<OneD, NekDouble> point(2);
    point[0] = m_points.back()[0];
    point[1] = m_points.back()[1];

    for (int j = 0; j < order; ++j)
    {
        point[0] += m_step * AdamsBashforth_coeffs[order - 1][j] *
                    cos(m_angles[numPoints - 1 - j]);
        point[1] += m_step * AdamsBashforth_coeffs[order - 1][j] *
                    sin(m_angles[numPoints - 1 - j]);
    }

    Array<OneD, NekDouble> lpoint(2);
    int id = m_f->m_exp[0]->GetExpIndex(point, lpoint);

    // We got out of the domain
    if (id == -1)
    {
        AddPoint(point);
        return false;
    }

    // We are still inside the domain and can continue
    int offset = m_f->m_exp[0]->GetPhys_Offset(id);

    NekDouble u = m_f->m_exp[0]->GetExp(id)->StdPhysEvaluate(
        lpoint, m_f->m_exp[0]->GetPhys() + offset);
    NekDouble v = m_f->m_exp[0]->GetExp(id)->StdPhysEvaluate(
        lpoint, m_f->m_exp[1]->GetPhys() + offset);

    if (m_neg.first)
    {
        if (m_neg.second && v > 0.0)
        {
            --m_rot;
        }
        else if (!m_neg.second && v < 0.0)
        {
            ++m_rot;
        }
    }

    m_neg.first  = u < 0.0;
    m_neg.second = v < 0.0;

    AddPoint(point, atan2(v, u) / 4.0 + m_rot * M_PI / 2.0);

    return true;
}
}
}
