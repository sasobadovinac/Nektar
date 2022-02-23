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

#include <iostream>
#include <string>
using namespace std;

#include "ProcessCrossField.h"

#include <LibUtilities/BasicUtils/ParseUtils.h>

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
    m_config["vertices"] = ConfigOption(false, "", "Vertices filename.");
    m_config["outcsv"]   = ConfigOption(false, "", "Output filename.");
    m_config["step"]     = ConfigOption(false, "", "Streamline step size.");
    m_config["mergetol"] = ConfigOption(false, "", "Tolerance for merging.");
}

ProcessCrossField::~ProcessCrossField()
{
}

void ProcessCrossField::Process(po::variables_map &vm)
{
    m_f->SetUpExp(vm);

    Initialise();

    vector<pair<Array<OneD, NekDouble>, int>> vertices = AnalyseVertices();

    vector<pair<Array<OneD, NekDouble>, int>> singularities =
        FindAllSingularities();

    vector<Streamline> allSls = CreateStreamlines(vertices, singularities);

    AdvanceStreamlines(allSls);

    Finalise(allSls);
}

void ProcessCrossField::Initialise()
{
    m_dim = m_f->m_graph->GetSpaceDimension();
    ASSERTL0(m_dim == 2, "Cross field post-processing only possible in 2D.")

    ASSERTL0(m_config["vertices"].m_beenSet, "A vertices file must be set.")
    ASSERTL0(m_config["outcsv"].m_beenSet, "An output CSV file must be set.")

    m_step =
        m_config["step"].m_beenSet ? m_config["step"].as<NekDouble>() : 1.0e-2;

    m_mergeTol = m_config["mergetol"].m_beenSet
                     ? m_config["mergetol"].as<NekDouble>()
                     : 1.0;
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

    // Store singularities and their number of quadrants
    vector<pair<Array<OneD, NekDouble>, int>> ret;

    // Find exact location of each singularity
    // We assume there is at most 1 singularity per element
    for (auto &elmt : intElmts)
    {
        Array<OneD, NekDouble> eta = FindSingularityInElmt(elmt);
        if (!eta.size())
        {
            // No singularity found
            continue;
        }

        int nQuadrants = CalculateNumberOfQuadrants(elmt, eta);

        Array<OneD, NekDouble> x(m_dim);
        m_f->m_exp[0]->GetExp(elmt)->GetCoord(eta, x);

        ret.push_back(make_pair(x, nQuadrants));

        cout << "Singularity #" << ret.size() << " has " << nQuadrants
             << " quadrants at (" << x[0] << ", " << x[1] << ")" << endl;
    }

    return ret;
}

vector<pair<Array<OneD, NekDouble>, int>> ProcessCrossField::AnalyseVertices()
{
    ifstream file;
    file.open(m_config["vertices"].as<string>());
    ASSERTL0(file.is_open(), "Could not open the vertices file.")

    vector<pair<Array<OneD, NekDouble>, int>> ret;

    string line;
    while (getline(file, line))
    {
        vector<NekDouble> data;
        ParseUtils::GenerateVector(line, data);

        NekDouble delta = numeric_limits<NekDouble>::max();
        Array<OneD, NekDouble> point(m_dim + 2);
        SpatialDomains::PointGeomSharedPtr geom;

        // Find closest point to the coordinates
        for (int i = 0; i < m_f->m_graph->GetNvertices(); ++i)
        {
            m_f->m_graph->GetVertex(i)->GetCoords(point);
            NekDouble newDelta =
                sqrt(pow(data[0] - point[0], 2) + pow(data[1] - point[1], 2));

            if (newDelta < delta)
            {
                delta = newDelta;
                geom  = m_f->m_graph->GetVertex(i);
            }
        }

        geom->GetCoords(point);

        // Find starting i for sweeping angle inside the domain ccw
        int nquads     = 100;
        NekDouble dist = m_step / 100;
        Array<OneD, NekDouble> quad(m_dim);
        Array<OneD, NekDouble> lcoords(m_dim);
        int i = 0;

        bool wasInside = false;
        bool first     = true;

        while (true)
        {
            quad[0] = point[0] + dist * cos(2.0 * M_PI * i / nquads);
            quad[1] = point[1] + dist * sin(2.0 * M_PI * i / nquads);

            int id = m_f->m_exp[0]->GetExpIndex(quad, lcoords);

            // outside
            if (id == -1)
            {
                ++i;

                if (!first && wasInside)
                {
                    break;
                }

                wasInside = false;
            }
            // inside
            else
            {
                if (!first && !wasInside)
                {
                    break;
                }

                --i;
                wasInside = true;
            }

            first = false;
        }

        // Determine the number of quadrants
        // We adjust the total angle of the corner by the difference in psi
        // We take into account jumps of atan2(v, u) at +-\pi too
        NekDouble oldV   = 0.0;
        NekDouble anglei = 0.0;
        NekDouble anglef = 0.0;
        NekDouble psii   = 0.0;
        NekDouble psif   = 0.0;
        int adj          = 0;

        first = true;

        while (true)
        {
            quad[0] = point[0] + dist * cos(2.0 * M_PI * i / nquads);
            quad[1] = point[1] + dist * sin(2.0 * M_PI * i / nquads);

            ++i;

            int id = m_f->m_exp[0]->GetExpIndex(quad, lcoords);

            if (id == -1)
            {
                break;
            }

            int offset = m_f->m_exp[0]->GetPhys_Offset(id);

            NekDouble u = m_f->m_exp[0]->GetExp(id)->StdPhysEvaluate(
                lcoords, m_f->m_exp[0]->GetPhys() + offset);
            NekDouble v = m_f->m_exp[0]->GetExp(id)->StdPhysEvaluate(
                lcoords, m_f->m_exp[1]->GetPhys() + offset);

            psif   = atan2(v, u) / 4.0;
            anglef = 2.0 * M_PI * (i - 1) / nquads;

            if (u < 0.0)
            {
                if (oldV * v < 0.0)
                {
                    if (v < 0.0)
                    {
                        --adj;
                    }
                    else
                    {
                        ++adj;
                    }
                }
            }

            oldV = v;

            if (first)
            {
                psii   = psif;
                anglei = anglef;
                first  = false;
            }
        }

        point[2] = anglei;
        point[3] = anglef;

        NekDouble fullAngle = point[3] - point[2];
        NekDouble dpsi      = psif - psii - adj * M_PI / 2.0;

        // Number of pi/2 quadrants
        int nQuadrants = int(round((fullAngle - dpsi) / (M_PI / 2.0)));

        ret.push_back(make_pair(point, nQuadrants));

        cout << "Vertex #" << ret.size() << " has " << nQuadrants
             << " quadrants at (" << point[0] << ", " << point[1] << ")"
             << endl;
    }

    return ret;
}

vector<Streamline> ProcessCrossField::CreateStreamlines(
    vector<pair<Array<OneD, NekDouble>, int>> &vertices,
    vector<pair<Array<OneD, NekDouble>, int>> &singularities)
{
    vector<Streamline> ret;

    for (auto &v : vertices)
    {
        Array<OneD, NekDouble> x(m_dim);
        x[0] = v.first[0];
        x[1] = v.first[1];

        NekDouble anglei = v.first[2];
        NekDouble anglef = v.first[3];

        // Number of streamlines to trace
        // No need to trace exterior streamlines: they correspond to CAD curves
        int nbranches = v.second - 1;

        NekDouble dir = anglei;

        for (int i = 0; i < nbranches; ++i)
        {
            dir += (anglef - anglei) / (nbranches + 1);

            Streamline sl = Streamline(m_f, x, m_step);
            sl.Initialise(dir);
            ret.push_back(sl);
        }

        // Collapsed cross-field at this vertex: need a splitting streamline
        if (nbranches == -1)
        {
            dir += anglef;
            dir *= 0.5;

            Streamline sl = Streamline(m_f, x, m_step, true);
            sl.Initialise(dir);
            ret.push_back(sl);
        }
    }

    for (auto &s : singularities)
    {
        Array<OneD, NekDouble> x = s.first;

        // Number of streamlines to trace
        int nbranches = s.second;

        NekDouble dir = 0.0;

        for (int i = 0; i < nbranches; ++i)
        {
            Streamline sl = Streamline(m_f, x, m_step);
            sl.Initialise(dir);
            ret.push_back(sl);

            dir += 2.0 * M_PI / nbranches;
        }
    }

    return ret;
}

void ProcessCrossField::AdvanceStreamlines(vector<Streamline> &sls)
{
    bool allLocked = false;

    while (!allLocked)
    {
        allLocked = true;

        for (auto &sl : sls)
        {
            allLocked = !sl.Advance() && allLocked;
        }

        for (int i = 0; i < sls.size(); ++i)
        {
            for (int j = 0; j < sls.size();)
            {
                if (i == j || sls[i].GetFirstPoint() == sls[j].GetFirstPoint())
                {
                    ++j;
                    continue;
                }

                int n = sls[i].CheckMerge(
                    sls[j],
                    m_step * ((sls[i].IsSplitter() || sls[j].IsSplitter())
                                  ? 1.0
                                  : m_mergeTol));

                if (n == -1)
                {
                    ++j;
                    continue;
                }

                if (!sls[i].IsSplitter() && !sls[j].IsSplitter())
                {
                    int n2 = sls[j].CheckMerge(sls[i], m_step * m_mergeTol);
                    if (n2 != -1)
                    {
                        n = min(n, n2);
                    }

                    for (int k = 0; k < n; ++k)
                    {
                        // We assume we won't get out of the domain
                        sls[i].Advance();
                        sls[j].Advance();
                    }

                    sls.push_back(sls[i].MergeWith(sls[j]));

                    if (i > j)
                    {
                        sls.erase(sls.begin() + i--);
                        sls.erase(sls.begin() + j);
                    }
                    else
                    {
                        sls.erase(sls.begin() + j);
                        sls.erase(sls.begin() + i);
                    }

                    j = 0;
                }
                else
                {
                    vector<Streamline> tmp;

                    if (sls[j].IsSplitter())
                    {
                        tmp = sls[j].ConvertSplitter(
                            sls[j].GetAllPoints().size());
                    }
                    else
                    {
                        tmp = sls[i].ConvertSplitter(n);
                    }

                    sls.insert(sls.end(), tmp.begin(), tmp.end());
                }
            }
        }
    }

    cout << "Final number of streamlines (after merging): " << sls.size()
         << endl;
}

void ProcessCrossField::Finalise(vector<Streamline> &sls)
{
    m_csvfile.open(m_config["outcsv"].as<string>());

    for (auto &sl : sls)
    {
        sl.WritePoints(m_csvfile);
    }

    m_csvfile.close();
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
                    lcoords[0] = -1.0;
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

        // Check if point is still inside the element
        // If not, we dismiss it;
        // it's most likely been detected by adjacent element
        // This probably should be outside of the loop because there is no
        // guarantee that the search path will not cross the boundary even
        // though there is indeed a singularity inside the element However,
        // StdExpansion2D has an ASSERTL2 on coordinates in v_PhysEvaluate that
        // blocks GitLab tests
        switch (expansion->GetGeom()->GetShapeType())
        {
            case LibUtilities::eTriangle:
                if (!(eta[0] >= -1.0 && eta[1] >= -1.0 &&
                      eta[0] + eta[1] <= 0.0))
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

    } while (sqrt(pow(delta[0], 2) + pow(delta[1], 2)) > 1.0e-9);

    return eta;
}

int ProcessCrossField::CalculateNumberOfQuadrants(int id,
                                                  Array<OneD, NekDouble> eta)
{
    // Determine the number of quadrants
    // In theory, we should be doing:
    // \( nQuadrants = 4 - \int_\pi {\psi'(\theta)} \)
    // In practice, we check for jumps of atan2(v, u) at +-\pi
    int ncquads    = 100;
    NekDouble dist = 1.0e-6;
    Array<OneD, NekDouble> cquad(m_dim);
    int nQuadrants = 4;
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
                    --nQuadrants;
                }
                else
                {
                    ++nQuadrants;
                }
            }
        }

        oldV = v;
    }

    return nQuadrants;
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
        // 1000 arbitrarily chosen; must just be a small step
        point[0] = m_points[0][0] + m_step * cos(angle) / 1000;
        point[1] = m_points[0][1] + m_step * sin(angle) / 1000;

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

    } while (fabs(delta) > 1.0e-9);

    // Over-write angle of singularity
    m_angles.back() = angle;

    m_neg.first  = u < 0.0;
    m_neg.second = v < 0.0;
}

bool Streamline::Advance()
{
    if (m_locked)
    {
        return false;
    }

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
        m_locked = true;
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

int Streamline::CheckMerge(Streamline sl, NekDouble tol)
{
    if ((m_locked || sl.IsLocked()) && !(m_splitter || sl.IsSplitter()))
    {
        return -1;
    }

    const Array<OneD, NekDouble> point = sl.GetLastPoint();
    const NekDouble angle = sl.GetLastAngle();

    int i;

    for (i = 0; i < m_points.size(); ++i)
    {
        // Check distance
        if (sqrt(pow(m_points[i][0] - point[0], 2) +
                 pow(m_points[i][1] - point[1], 2)) < tol)
        {
            if (m_splitter || sl.IsSplitter())
            {
                break;
            }

            // Check opposite directions
            if (abs(int(round((m_angles[i] - angle) / (M_PI / 2.0)))) % 4 == 2)
            {
                break;
            }
        }
    }

    return i == m_points.size() ? -1 : i;
}

Streamline Streamline::MergeWith(Streamline sl)
{
    vector<Array<OneD, NekDouble>> points = sl.GetAllPoints();

    ASSERTL0(m_points.size() == points.size(),
             "Problem with the number of points to merge.");

    int npoints = m_points.size();
    vector<Array<OneD, NekDouble>> mergedPoints(npoints);

    int nweights = npoints - 1;

    for (int i = 0; i < npoints; ++i)
    {
        mergedPoints[i] = Array<OneD, NekDouble>(m_dim, 0.0);

        /*
        // Linear interpolation
        mergedPoints[i][0] =
            (m_points[i][0] * (nweights - i) + points[nweights - i][0] * i) /
            nweights;
        mergedPoints[i][1] =
            (m_points[i][1] * (nweights - i) + points[nweights - i][1] * i) /
            nweights;
        */
       
        // Trigonometric interpolation
        mergedPoints[i][0] =
            m_points[i][0] * pow(cos(M_PI / 2.0 * i / nweights), 2) +
            points[nweights - i][0] * pow(sin(M_PI / 2.0 * i / nweights), 2);
        mergedPoints[i][1] =
            m_points[i][1] * pow(cos(M_PI / 2.0 * i / nweights), 2) +
            points[nweights - i][1] * pow(sin(M_PI / 2.0 * i / nweights), 2);
    }

    return Streamline(mergedPoints);
}

vector<Streamline> Streamline::ConvertSplitter(int n)
{
    m_splitter = false;

    vector<Streamline> ret;

    m_points.erase(m_points.begin(), m_points.begin() + 2 * n / 3);
    m_angles.erase(m_angles.begin(), m_angles.begin() + 2 * n / 3);

    Array<OneD, NekDouble> p0 = m_points[0];
    Array<OneD, NekDouble> p1 = m_points[1];

    NekDouble angle = atan2(p1[1] - p0[1], p1[0] - p0[0]);

    for (int i = 0; i < 2; ++i)
    {
        vector<Array<OneD, NekDouble>> points(2);

        points[0]    = Array<OneD, NekDouble>(m_dim);
        points[0][0] = p0[0];
        points[0][1] = p0[1];

        points[1]    = Array<OneD, NekDouble>(m_dim);
        points[1][0] = p0[0];
        points[1][1] = p0[1];

        Array<OneD, NekDouble> vec(m_dim);
        vec[0] = m_step * cos(angle + (i + 1) / 3.0 * 2.0 * M_PI);
        vec[1] = m_step * sin(angle + (i + 1) / 3.0 * 2.0 * M_PI);

        do
        {
            Vmath::Vadd(m_dim, points[1], 1, vec, 1, points[1], 1);
        } while (m_f->m_exp[0]->GetExpIndex(points[1]) != -1);

        ret.emplace_back(points);
    }

    return ret;
}

void Streamline::WritePoints(ofstream &csvfile)
{
    for (auto &p : m_points)
    {
        csvfile << p[0] << "," << p[1] << endl;
    }

    csvfile << endl;
}
}
}
