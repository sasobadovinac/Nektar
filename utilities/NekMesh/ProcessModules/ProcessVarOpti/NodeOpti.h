////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.h
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

#ifndef UTILITIES_NEKMESH_NODEOPTI
#define UTILITIES_NEKMESH_NODEOPTI

#include <iostream>

#include "../../Module.h"
#include "ProcessVarOpti.h"

#include <LibUtilities/BasicUtils/Thread.h>

namespace Nektar
{
namespace Utilities
{


class NodeOptiJob;

class NodeOpti
{
public:
    NodeOpti(NodeSharedPtr n,
             std::pair<std::vector<int>, std::vector<ElUtilSharedPtr> > e,
             std::map<LibUtilities::ShapeType,DerivUtilSharedPtr> d,
             optimiser o)
        : node(n), nodeIds(e.first), data(e.second), derivUtil(d), opti(o)
    {
    }

    virtual ~NodeOpti(){};
        
    std::vector<ElUtilSharedPtr> data;
    NodeSharedPtr node;

    static void Node_indexing(NodesGPU &nodes, 
            std::vector<std::vector<boost::shared_ptr<NodeOpti>> > &optiNodes);
    

protected:

    std::vector<int> nodeIds;
    
    bool Linear();
     

    NekDouble dx;
    NekDouble minJac;
    Residual res;
    std::map<LibUtilities::ShapeType,DerivUtilSharedPtr> derivUtil;
    optimiser opti;

    static const NekDouble gam;

};

typedef boost::shared_ptr<NodeOpti> NodeOptiSharedPtr;
typedef LibUtilities::NekFactory<int,
                                 NodeOpti,
                                 NodeSharedPtr,
                                 std::pair<std::vector<int>,
                                           std::vector<ElUtilSharedPtr> >,
                                 std::map<LibUtilities::ShapeType,DerivUtilSharedPtr>,
                                 optimiser> NodeOptiFactory;

NodeOptiFactory &GetNodeOptiFactory();


class NodeOpti2D2D : public NodeOpti // 1D optimsation in 3D space
{
public:
    NodeOpti2D2D(NodeSharedPtr n,
                 std::pair<std::vector<int>, std::vector<ElUtilSharedPtr> > e,
                 std::map<LibUtilities::ShapeType,DerivUtilSharedPtr> d,
                 optimiser o)
                 : NodeOpti(n,e,d,o)
    {
    }

    ~NodeOpti2D2D(){};
    
    static int m_type;
    static NodeOptiSharedPtr create(
        NodeSharedPtr n,
        std::pair<std::vector<int>, std::vector<ElUtilSharedPtr> > e,
        std::map<LibUtilities::ShapeType,DerivUtilSharedPtr> d,
        optimiser o)
    {
        return NodeOptiSharedPtr(new NodeOpti2D2D(n, e, d, o));
    }

private:
};

class NodeOpti3D3D : public NodeOpti //1D optimsation in 3D space
{
public:
    NodeOpti3D3D(NodeSharedPtr n,
                 std::pair<std::vector<int>, std::vector<ElUtilSharedPtr> > e,
                 std::map<LibUtilities::ShapeType,DerivUtilSharedPtr> d,
                 optimiser o)
                 : NodeOpti(n,e,d,o)
    {
    }

    ~NodeOpti3D3D(){};    

    static int m_type;
    static NodeOptiSharedPtr create(
        NodeSharedPtr n,
        std::pair<std::vector<int>, std::vector<ElUtilSharedPtr> > e,
        std::map<LibUtilities::ShapeType,DerivUtilSharedPtr> d,
        optimiser o)
    {
        return NodeOptiSharedPtr(new NodeOpti3D3D(n, e, d, o));
    }

private:

};

}
}

#endif
