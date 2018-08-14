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

#include <SpatialDomains/Geometry.h>

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

    vector<set<shared_ptr<SpatialDomains::Geometry>>> isoElmts(2);

    for (int i = 0; i < 2; ++i)
    {
        Array<OneD, NekDouble> phys = m_f->m_exp[i]->GetPhys();

        for (int elmt = 0; elmt < m_f->m_exp[i]->GetNumElmts(); ++elmt)
        {
            int offset  = m_f->m_exp[i]->GetPhys_Offset(elmt);
            int npoints = m_f->m_exp[i]->GetTotPoints(elmt);

            auto mm = minmax_element(phys.begin() + offset,
                                     phys.begin() + offset + npoints);

            if (*(mm.first) * *(mm.second) < 0.0)
            {
                isoElmts[i].insert(m_f->m_exp[i]->GetExp(elmt)->GetGeom());
            }
        }
    }

    vector<shared_ptr<SpatialDomains::Geometry>> intElmts(
        min(isoElmts[0].size(), isoElmts[1].size()));
    auto it = set_intersection(isoElmts[0].begin(), isoElmts[0].end(),
                               isoElmts[1].begin(), isoElmts[1].end(),
                               intElmts.begin());
    intElmts.resize(it - intElmts.begin());
}
}
}
