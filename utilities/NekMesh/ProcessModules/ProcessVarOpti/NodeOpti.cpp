////////////////////////////////////////////////////////////////////////////////
//
//  File: NodeOpti.cpp
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

#include <boost/multi_array.hpp>

#include <limits>

#include "NodeOpti.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{


NodeOptiFactory &GetNodeOptiFactory()
{
    
    static NodeOptiFactory asd;
    return asd;
}

int NodeOpti3D3D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    33, NodeOpti3D3D::create, "3D3D");


void NodeOpti::Node_indexing(NodesGPU &nodes, 
		std::vector<std::vector<boost::shared_ptr<NodeOpti>> > &optiNodes)
{
	// construct vector that contains the number of nodes in each colourset
    int nColoursets = optiNodes.size();
    nodes.coloursetSize = Kokkos::View<int*> ("coloursetSize",nColoursets);
    nodes.h_coloursetSize = Kokkos::create_mirror_view(nodes.coloursetSize);
    int maxColoursetSize = 0;
    for(int cs = 0; cs < nColoursets; cs++)
    {            
        nodes.h_coloursetSize(cs) = optiNodes[cs].size();
        maxColoursetSize = (maxColoursetSize < nodes.h_coloursetSize(cs)
                ? nodes.h_coloursetSize(cs) : maxColoursetSize);
    }
    //construct Vector that contains the number of elements that each node in the colorset is associated with
    nodes.nElmtArray = Kokkos::View<int**> ("nElmtArray",nColoursets,maxColoursetSize);
    nodes.h_nElmtArray = Kokkos::create_mirror_view(nodes.nElmtArray);

    int maxnElmt = 0;
    for(int cs = 0; cs < nColoursets; cs++)
    {
        for(int node = 0; node < nodes.h_coloursetSize(cs); node++)
        {
            nodes.h_nElmtArray(cs,node) = optiNodes[cs][node]->data.size();
            maxnElmt = (maxnElmt < nodes.h_nElmtArray(cs,node)
                    ? nodes.h_nElmtArray(cs,node) : maxnElmt);
        }
    }
    // construct Array that contains the element Ids of all elements that each node in the colorset is associated with
    nodes.elIdArray = Kokkos::View<int***> ("elIdArray",nColoursets,maxColoursetSize, maxnElmt);
    nodes.h_elIdArray = Kokkos::create_mirror_view(nodes.elIdArray);

    // construct Array that contains the local node Ids in all elements that each node in the colorset is associated with
    nodes.localNodeIdArray = Kokkos::View<int***> ("localNodeIdArray",nColoursets,maxColoursetSize, maxnElmt);
    nodes.h_localNodeIdArray = Kokkos::create_mirror_view(nodes.localNodeIdArray);

    for(int cs = 0; cs < nColoursets; cs++)
    {
        for(int node = 0; node < nodes.h_coloursetSize(cs); node++)            
        {
            const int globalNodeId = optiNodes[cs][node]->node->m_id;

            for (int el = 0; el < nodes.h_nElmtArray(cs,node); ++el)
            {                
                int globalEl = optiNodes[cs][node]->data[el]->GetId();
                nodes.h_elIdArray(cs,node,el) = globalEl;

                for(int ln = 0; ln < nodes.nodes_size; ln++)
                {
                    if(nodes.h_Id(globalEl,ln) == globalNodeId)
                    {
                        nodes.h_localNodeIdArray(cs,node,el) = ln;
                    }
                }
            }
        }
    }

    Kokkos::deep_copy(nodes.localNodeIdArray,nodes.h_localNodeIdArray);
    Kokkos::deep_copy(nodes.elIdArray,nodes.h_elIdArray);
    Kokkos::deep_copy(nodes.nElmtArray,nodes.h_nElmtArray);
    Kokkos::deep_copy(nodes.coloursetSize,nodes.h_coloursetSize);
}



}
}
