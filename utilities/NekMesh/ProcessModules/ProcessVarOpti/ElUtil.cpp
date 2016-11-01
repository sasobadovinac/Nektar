////////////////////////////////////////////////////////////////////////////////
//
//  File: ElUtil.cpp
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
//  Description: Calculate jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include "ElUtil.h"
#include "ProcessVarOpti.h"

#include <LibUtilities/Foundations/ManagerAccess.h>

#include <Kokkos_Core.hpp>

using namespace std;

namespace Nektar
{
namespace Utilities
{


ElUtil::ElUtil(ElementSharedPtr e, DerivUtilSharedPtr d, int n)
{
    m_el = e;
    derivUtil = d;
    m_mode = n;
    m_dim = m_el->GetDim();
    vector<NodeSharedPtr> ns;
    m_el->GetCurvedNodes(ns);
    nodes.resize(ns.size());
    nodeIds.resize(ns.size());
    for (int i = 0; i < ns.size(); ++i)
    {
        nodes[i].resize(m_dim+1);
        nodes[i][0] = &ns[i]->m_x;

        if (m_dim >= 2)
        {
            nodes[i][1] = &ns[i]->m_y;
        }

        if (m_dim >= 3)
        {
            nodes[i][2] = &ns[i]->m_z;
        }
        nodeIds[i] = &ns[i]->m_id;
    }
    maps = MappingIdealToRef();
}

vector<Array<OneD, NekDouble> > ElUtil::MappingIdealToRef()
{
    vector<Array<OneD, NekDouble> > ret;

    if(m_el->GetConf().m_e == LibUtilities::eQuadrilateral)
    {
        ASSERTL0(false,"Not coded");
        
    }
    else if(m_el->GetConf().m_e == LibUtilities::eTriangle)
    {
        DNekMat J(2,2,0.0);
        J(0,0) = (*nodes[1][0] - *nodes[0][0]);
        J(1,0) = (*nodes[1][1] - *nodes[0][1]);
        J(0,1) = (*nodes[2][0] - *nodes[0][0]);
        J(1,1) = (*nodes[2][1] - *nodes[0][1]);

        J.Invert();

        DNekMat R(2,2,0.0);
        R(0,0) = 2.0;
        R(1,1) = 2.0;

        J = J * R;

        for(int i = 0 ; i < derivUtil->ptsHigh; i++)
        {
            Array<OneD, NekDouble> r(10,0.0); //store det in 10th entry

            r[9] = 1.0 / (J(0,0) * J(1,1) - J(0,1) * J(1,0));
            r[0] = J(0,0);
            r[1] = J(1,0);
            r[2] = 0.0;
            r[3] = J(0,1);
            r[4] = J(1,1);
            r[5] = 0.0;
            r[6] = 0.0;
            r[7] = 0.0;
            r[8] = 0.0;
            ret.push_back(r);
        }
    }
    else if(m_el->GetConf().m_e == LibUtilities::eTetrahedron)
    {
        DNekMat J(3,3,0.0);
        J(0,0) = (*nodes[1][0] - *nodes[0][0]);
        J(1,0) = (*nodes[1][1] - *nodes[0][1]);
        J(2,0) = (*nodes[1][2] - *nodes[0][2]);
        J(0,1) = (*nodes[2][0] - *nodes[0][0]);
        J(1,1) = (*nodes[2][1] - *nodes[0][1]);
        J(2,1) = (*nodes[2][2] - *nodes[0][2]);
        J(0,2) = (*nodes[3][0] - *nodes[0][0]);
        J(1,2) = (*nodes[3][1] - *nodes[0][1]);
        J(2,2) = (*nodes[3][2] - *nodes[0][2]);

        J.Invert();

        DNekMat R(3,3,0.0);
        R(0,0) = 2.0;
        R(1,1) = 2.0;
        R(2,2) = 2.0;

        J = J * R;

        for(int i = 0 ; i < derivUtil->ptsHigh; i++)
        {
            Array<OneD, NekDouble> r(10,0.0); //store det in 10th entry

            r[9] = 1.0/(J(0,0)*(J(1,1)*J(2,2)-J(2,1)*J(1,2))
                       -J(0,1)*(J(1,0)*J(2,2)-J(2,0)*J(1,2))
                       +J(0,2)*(J(1,0)*J(2,1)-J(2,0)*J(1,1)));

            r[0] = J(0,0);
            r[1] = J(1,0);
            r[2] = J(2,0);
            r[3] = J(0,1);
            r[4] = J(1,1);
            r[5] = J(2,1);
            r[6] = J(0,2);
            r[7] = J(1,2);
            r[8] = J(2,2);
            ret.push_back(r);
        }
    }
    else if(m_el->GetConf().m_e == LibUtilities::ePrism)
    {
        ASSERTL0(false, "not coded");
        
    }
    else
    {
        ASSERTL0(false,"not coded");
    }

    return ret;
}



}
}
