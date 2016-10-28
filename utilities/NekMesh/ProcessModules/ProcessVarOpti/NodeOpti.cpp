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
#include "Evaluator.hxx"
#include "Hessian.hxx"

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

const NekDouble NodeOpti::gam = numeric_limits<float>::epsilon();

KOKKOS_INLINE_FUNCTION
double NodeOpti::CalcMinJac( ElUtilGPU &elUtil, int nElmt, typename Kokkos::View<int*>::HostMirror elIdArray)
{
    double minJac = DBL_MAX;
    int el;

    for(int i = 0; i < nElmt; i++)
    {
        el = elIdArray[i];
        minJac = (minJac > elUtil.h_minJac(el) ? elUtil.h_minJac(el) : minJac);
            
    }
    return minJac;
}



KOKKOS_INLINE_FUNCTION
void NodeOpti::GetNodeCoord(double (&X)[3], int id, NodesGPU &nodes,
            typename Kokkos::View<int*>::HostMirror elIdArray, typename Kokkos::View<int*>::HostMirror localNodeIdArray)
{   
    int elmt = elIdArray[0];
    int node = localNodeIdArray[0];
    Kokkos::View<double[3]> Xt("Xt");
    typename Kokkos::View< double[3]>::HostMirror h_Xt = Kokkos::create_mirror_view(Xt);
    Kokkos::parallel_for(range_policy(0,1), KOKKOS_LAMBDA (const int k)
    {
        Xt(0) = nodes.X(elmt,node);
        Xt(1) = nodes.Y(elmt,node);
        Xt(2) = nodes.Z(elmt,node);
    });    
    Kokkos::deep_copy(h_Xt,Xt);
    X[0] = h_Xt(0);
    X[1] = h_Xt(1);
    X[2] = h_Xt(2); 
}

KOKKOS_INLINE_FUNCTION
void NodeOpti::SetNodeCoord(double (&X)[3], int id, NodesGPU &nodes,
             typename Kokkos::View<int*>::HostMirror elIdArray, typename Kokkos::View<int*>::HostMirror localNodeIdArray, int nElmt)
{
    Kokkos::View<double[3]> Xt("Xt");
    typename Kokkos::View< double[3]>::HostMirror h_Xt = Kokkos::create_mirror_view(Xt);
    h_Xt(0) = X[0];
    h_Xt(1) = X[1];
    h_Xt(2) = X[2];
    Kokkos::deep_copy(Xt,h_Xt);
    for(int el = 0; el < nElmt; el++)
    {
        int elmt = elIdArray[el];
        int node = localNodeIdArray[el];
        Kokkos::parallel_for(range_policy(0,1), KOKKOS_LAMBDA (const int k)
        {
            nodes.X(elmt,node) = Xt(0);
            nodes.Y(elmt,node) = Xt(1);
            nodes.Z(elmt,node) = Xt(2);
        });
    }   
}




/*int NodeOpti2D2D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    22, NodeOpti2D2D::create, "2D2D");

void NodeOpti2D2D::Optimise(DerivUtilGPU &derivUtil,NodesGPU &nodes, 
        NodeMap &nodeMap, ElUtilGPU &elUtil, Residual &res, 
        int nElmt, int globalNodeId)
{

}*/




int NodeOpti3D3D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    33, NodeOpti3D3D::create, "3D3D");

KOKKOS_INLINE_FUNCTION
void NodeOpti3D3D::Optimise(DerivUtilGPU &derivUtil,NodesGPU &nodes, 
        NodeMap &nodeMap, ElUtilGPU &elUtil, Residual &res, 
        int nElmt)
{      
    /*Grad grad;
    grad.G = Kokkos::View<double[9]> ("G");
    grad.h_G = Kokkos::create_mirror_view(grad.G);
    grad.integral = Kokkos::View<double[1]> ("integral");
    grad.h_integral = Kokkos::create_mirror_view(grad.integral);    
    
    Kokkos::parallel_for( team_policy( 1, Kokkos::AUTO ), KOKKOS_LAMBDA ( const member_type& teamMember)
    {
        //const int el = teamMember.league_rank();
        
        double currentW, newVal;

        double minJac = CalcMinJacGPU(elUtil, nElmt, nodes.elIdArray);
        double ep = minJac < 0.0 ? sqrt(1e-9 + 0.04*minJac*minJac) : sqrt(1e-9); 
        
        

        currentW = NodeOpti::GetFunctional<3>(derivUtil, nodes, elUtil, grad, nElmt, ep, teamMember);
    
    
        if(grad.G[0]*grad.G[0] + grad.G[1]*grad.G[1] + grad.G[2]*grad.G[2] > gradTol())
        {        
            double h_Xc[3];        
            double h_Xn[3];        
            GetNodeCoordGPU(h_Xc, nodes, nodes.elIdArray, nodes.localNodeIdArray); 

            double sk[3];
            CalcSK(grad.G, sk);
            double pg = (grad.G[0]*sk[0]+grad.G[1]*sk[1]+grad.G[2]*sk[2]);
            double hes    = sk[0] * (sk[0]*grad.G[3] + sk[1]*grad.G[4] + sk[2]*grad.G[5]) +
                            sk[1] * (sk[0]*grad.G[4] + sk[1]*grad.G[6] + sk[2]*grad.G[7]) +
                            sk[2] * (sk[0]*grad.G[5] + sk[1]*grad.G[7] + sk[2]*grad.G[8]);
            hes = (hes > 0.0 ? 0.0 : hes);

            double alpha  = 1.0;
            bool found  = false;
            //while (alpha > alphaTol())
            for(int c  = 0; c < 34; c++)
            {
            // Update node                
                h_Xn[0] = h_Xc[0] + alpha * sk[0];
                h_Xn[1] = h_Xc[1] + alpha * sk[1];
                h_Xn[2] = h_Xc[2] + alpha * sk[2];
                SetNodeCoordGPU(h_Xn, nodes, nodes.elIdArray, nodes.localNodeIdArray, nElmt);
                
                
                newVal = NodeOpti::GetFunctional<3>(derivUtil, nodes, elUtil, grad, nElmt, ep, teamMember);
           

                if (newVal <= currentW + c1() * (alpha*pg+ 0.5*alpha*alpha*hes))
                {
                    found = true;
                    break;
                }
                alpha /= 2.0;        
            }        

            if(!found)
            {
                h_Xn[0] = h_Xc[0];
                h_Xn[1] = h_Xc[1];
                h_Xn[2] = h_Xc[2];
                SetNodeCoordGPU(h_Xn, nodes, nodes.elIdArray, nodes.localNodeIdArray, nElmt);

                //mtx.lock();
                Kokkos::single(Kokkos::PerTeam(teamMember),[&] ()
                {
                    res.nReset[0]++;
                    printf("%s\n", "3D reset");
                });            
                //mtx.unlock();
            }

            // update residual value
            double thisval = sqrt( (h_Xn[0]-h_Xc[0])*(h_Xn[0]-h_Xc[0])
                                  +(h_Xn[1]-h_Xc[1])*(h_Xn[1]-h_Xc[1])
                                  +(h_Xn[2]-h_Xc[2])*(h_Xn[2]-h_Xc[2]) );

            //mtx.lock();
            Kokkos::single(Kokkos::PerTeam(teamMember),[&] ()
            {
                //res.val[0] = max(thisval, res.val[0] );
                res.val[0] = (res.val[0] < thisval ? thisval : res.val[0]);
                res.func[0] +=newVal;
            });
            //mtx.unlock();

        }

    });*/

}



}
}
