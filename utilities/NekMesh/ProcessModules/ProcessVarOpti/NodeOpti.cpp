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
    /*
    typedef Loki::SingletonHolder<NodeOptiFactory, Loki::CreateUsingNew,
                                  Loki::NoDestroy, Loki::ClassLevelLockable> Type;
    return Type::Instance();
    */
    //typedef Loki::SingletonHolder<NodeOptiFactory, Loki::CreateUsingNew,
    //                              Loki::NoDestroy, Loki::ClassLevelLockable> Type;
    //return Type::Instance();

    static NodeOptiFactory asd;
    return asd;
}

const NekDouble NodeOpti::gam = numeric_limits<float>::epsilon();

void NodeOpti::CalcMinJac()
{
    minJac = numeric_limits<double>::max();

    for(int i = 0; i < data.size(); i++)
    {
        minJac = min(minJac, data[i]->minJac);
    }
}
double NodeOpti::CalcMinJac(ElUtilGPU &elUtil, int nElmt, int * elIdArray)
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

void NodeOpti::GetNodeCoord(double (&X)[3], int id, NodesGPU &nodes,
            int * elIdArray, int * localNodeIdArray)
{
    //x = node->m_x;
    //y = node->m_y;
    //z = node->m_z;

    // it is sufficient to pull the coordinates from the first instance of the node
    ////NodeMap::const_iterator coeffs;
    ////coeffs = nodeMap.find(id);    
    ////int elmt = std::get<0>(coeffs->second);
    ////int node = std::get<1>(coeffs->second);

    //x = nodes.h_X(elmt,node);
    //y = nodes.h_Y(elmt,node);
    //z = nodes.h_Z(elmt,node);

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

void NodeOpti::SetNodeCoord(double (&X)[3], int id, NodesGPU &nodes,
             int * elIdArray, int * localNodeIdArray, int nElmt)
{
    //node->m_x = x;
    //node->m_y = y;
    //node->m_z = z;
    
    Kokkos::View<double[3]> Xt("Xt");
    typename Kokkos::View< double[3]>::HostMirror h_Xt = Kokkos::create_mirror_view(Xt);
    h_Xt(0) = X[0];
    h_Xt(1) = X[1];
    h_Xt(2) = X[2];
    Kokkos::deep_copy(Xt,h_Xt);

    ////NodeMap::const_iterator coeffs;
    ////coeffs = nodeMap.find(id);    
    for(int el = 0; el < nElmt; el++)
    {
        ////int elmt = std::get<0>(coeffs->second);
        ////int node = std::get<1>(coeffs->second);
        int elmt = elIdArray[el];
        int node = localNodeIdArray[el];
        //nodes.h_X(elmt,node) = x;
        //nodes.h_Y(elmt,node) = y;
        //nodes.h_Z(elmt,node) = z;
        Kokkos::parallel_for(range_policy(0,1), KOKKOS_LAMBDA (const int k)
        {
            nodes.X(elmt,node) = Xt(0);
            nodes.Y(elmt,node) = Xt(1);
            nodes.Z(elmt,node) = Xt(2);
        });
        ////coeffs++;
    }    

}


int NodeOpti2D2D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    22, NodeOpti2D2D::create, "2D2D");

void NodeOpti2D2D::Optimise(DerivUtilGPU &derivUtil,NodesGPU &nodes, 
        NodeMap &nodeMap, ElUtilGPU &elUtil, Residual &res, 
        int nElmt, int globalNodeId, int * elIdArray, int * localNodeIdArray)
{

}

int NodeOpti3D3D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    33, NodeOpti3D3D::create, "3D3D");




void NodeOpti3D3D::Optimise(DerivUtilGPU &derivUtil,NodesGPU &nodes, 
        NodeMap &nodeMap, ElUtilGPU &elUtil, Residual &res, 
        int nElmt, int globalNodeId, int * elIdArray, int * localNodeIdArray)
{
       
    // using the node ID and the node map find the local coordinates of the node,
                // depending on the considered element
    /*const int nElmt = data.size();
    const int globalNodeId = node->m_id;
    
    int elIdArray[nElmt];    
    int localNodeIdArray[nElmt];

    NodeMap::const_iterator coeffs;
    coeffs = nodeMap.find(globalNodeId); 
    for (int el = 0; el < nElmt; ++el)
    {
        elIdArray[el] = data[el]->GetId();

        int elmt = std::get<0>(coeffs->second);
        if (elmt == elIdArray[el])
        {
            localNodeIdArray[el] = std::get<1>(coeffs->second);
        }            
        coeffs++;
    }*/


    double minJac = CalcMinJac(elUtil, nElmt, elIdArray);
    double ep = minJac < 0.0 ? sqrt(1e-9 + 0.04*minJac*minJac) : sqrt(1e-9);    

    Grad grad;
    grad.G = Kokkos::View<double[9]> ("G");
    grad.h_G = Kokkos::create_mirror_view(grad.G);
    grad.integral = Kokkos::View<double[1]> ("integral");
    grad.h_integral = Kokkos::create_mirror_view(grad.integral);

    //Kokkos::parallel_for(range_policy_host(0,1), KOKKOS_LAMBDA (const int dum)
    //{

    double currentW, newVal, dbVal;
    int elId, localNodeId;    
    grad.h_integral[0] = 0.0;
    Kokkos::deep_copy(grad.integral,grad.h_integral);
    for (int el = 0; el < nElmt; ++el)
    { 
        elId = elIdArray[el];
        localNodeId = localNodeIdArray[el];
        GetFunctional<3>(derivUtil, nodes, elUtil, nodeMap, grad,
            elId, localNodeId, ep);
    }
    Kokkos::deep_copy(grad.h_integral, grad.integral);
    currentW = grad.h_integral[0];
    Kokkos::deep_copy(grad.h_G, grad.G);

    
    
    
    if(grad.h_G[0]*grad.h_G[0] + grad.h_G[1]*grad.h_G[1] + grad.h_G[2]*grad.h_G[2] > gradTol())
    {
        
        double h_Xc[3];        
        double h_Xn[3];

        //GetNodeCoord(h_Xc, globalNodeId, nodes, nodeMap); 

        //calculate sk
        double sk[3];
        double det = grad.h_G[3]*(grad.h_G[6]*grad.h_G[8]-grad.h_G[7]*grad.h_G[7])
                       -grad.h_G[4]*(grad.h_G[4]*grad.h_G[8]-grad.h_G[5]*grad.h_G[7])
                       +grad.h_G[5]*(grad.h_G[4]*grad.h_G[7]-grad.h_G[5]*grad.h_G[6]);

        sk[0] = grad.h_G[0]*(grad.h_G[6]*grad.h_G[8]-grad.h_G[7]*grad.h_G[7]) +
                grad.h_G[1]*(grad.h_G[5]*grad.h_G[7]-grad.h_G[4]*grad.h_G[8]) +
                grad.h_G[2]*(grad.h_G[4]*grad.h_G[7]-grad.h_G[3]*grad.h_G[7]);
        sk[1] = grad.h_G[0]*(grad.h_G[7]*grad.h_G[5]-grad.h_G[4]*grad.h_G[5]) +
                grad.h_G[1]*(grad.h_G[3]*grad.h_G[8]-grad.h_G[5]*grad.h_G[5]) +
                grad.h_G[2]*(grad.h_G[4]*grad.h_G[5]-grad.h_G[3]*grad.h_G[7]);
        sk[2] = grad.h_G[0]*(grad.h_G[4]*grad.h_G[7]-grad.h_G[6]*grad.h_G[5]) +
                grad.h_G[1]*(grad.h_G[4]*grad.h_G[5]-grad.h_G[3]*grad.h_G[7]) +
                grad.h_G[2]*(grad.h_G[3]*grad.h_G[6]-grad.h_G[4]*grad.h_G[4]);

        sk[0] /= det * -1.0;
        sk[1] /= det * -1.0;
        sk[2] /= det * -1.0;

        double pg = (grad.h_G[0]*sk[0]+grad.h_G[1]*sk[1]+grad.h_G[2]*sk[2]);

        bool runDNC = false; //so want to make this varible runDMC
        bool found  = false;        
        
        //checks for DNC
        /*double dk[3];
        double lhs;
        int def = IsIndefinite<3>(grad);
        if(def) // def != 0
        {
            //the dk vector needs calculating
            double val;
            MinEigen<3>(val,dk, grad);

            if(dk[0]*grad.h_G[0] + dk[1]*grad.h_G[1] + dk[2]*grad.h_G[2] > 0.0)
            {
                for(int i = 0; i < 3; i++)
                {
                    dk[i] *= -1.0;
                }
            }

            lhs = dk[0] * (dk[0]*grad.h_G[3] + dk[1]*grad.h_G[4] + dk[2]*grad.h_G[5]) +
                  dk[1] * (dk[0]*grad.h_G[4] + dk[1]*grad.h_G[6] + dk[2]*grad.h_G[7]) +
                  dk[2] * (dk[0]*grad.h_G[5] + dk[1]*grad.h_G[7] + dk[2]*grad.h_G[8]);

            ASSERTL0(lhs < 0.0 , "weirdness");
        
            NekDouble skmag = sqrt(sk[0]*sk[0] + sk[1]*sk[1] + sk[2]*sk[2]);
            runDNC = !((grad.h_G[0]*sk[0]+grad.h_G[1]*sk[1]+grad.h_G[2]*sk[2])/skmag <=
                        2.0*(0.5*lhs + grad.h_G[0]*dk[0]+grad.h_G[1]*dk[1]+grad.h_G[2]*dk[2]));
        }*/

        

        if(!runDNC)
        {
            //normal gradient line Search
            GetNodeCoord(h_Xc, globalNodeId, nodes, elIdArray, localNodeIdArray); 

            double alpha  = 1.0;
            double hes    = sk[0] * (sk[0]*grad.h_G[3] + sk[1]*grad.h_G[4] + sk[2]*grad.h_G[5]) +
                            sk[1] * (sk[0]*grad.h_G[4] + sk[1]*grad.h_G[6] + sk[2]*grad.h_G[7]) +
                            sk[2] * (sk[0]*grad.h_G[5] + sk[1]*grad.h_G[7] + sk[2]*grad.h_G[8]);
            //hes = min(hes,0.0);
            hes = (hes > 0.0 ? 0.0 : hes);
            
            while (alpha > alphaTol())
            {
            // Update node                
                h_Xn[0] = h_Xc[0] + alpha * sk[0];
                h_Xn[1] = h_Xc[1] + alpha * sk[1];
                h_Xn[2] = h_Xc[2] + alpha * sk[2];
                SetNodeCoord(h_Xn, globalNodeId, nodes, elIdArray, localNodeIdArray, nElmt);
                
                grad.h_integral[0] = 0.0;
                Kokkos::deep_copy(grad.integral,grad.h_integral);
                for (int el = 0; el < nElmt; ++el)
                {
                    elId = elIdArray[el];
                    localNodeId = localNodeIdArray[el];
                    GetFunctional<3>(derivUtil, nodes, elUtil, nodeMap, grad,
                            elId, localNodeId, ep, false,false);
                }
                Kokkos::deep_copy(grad.h_integral, grad.integral);
                newVal = grad.h_integral[0];
            // end evaluate node

                //dont need the hessian again this function updates G to be the new
                //location
                //
                // Wolfe conditions
                if (newVal <= currentW + c1() * (alpha*pg+ 0.5*alpha*alpha*hes))
                {
                    found = true;
                    break;
                }

                alpha /= 2.0;
            }
        }
        /*else
        {
            GetNodeCoord(h_Xc, globalNodeId, nodes, elIdArray, localNodeIdArray); 

            NekDouble beta = 0.5;
            int l = 0;
            NekDouble alpha = pow(beta,l);

            NekDouble hes = lhs;

            pg = (grad.h_G[0]*dk[0]+grad.h_G[1]*dk[1]+grad.h_G[2]*dk[2]);

        //choose whether to do forward or reverse line search
            h_Xn[0] = h_Xc[0] + dk[0];
            h_Xn[1] = h_Xc[1] + dk[1];
            h_Xn[2] = h_Xc[2] + dk[2];
            SetNodeCoord(h_Xn, globalNodeId, nodes, elIdArray, localNodeIdArray, nElmt);

            grad.h_integral[0] = 0.0;
            Kokkos::deep_copy(grad.integral,grad.h_integral);
            for (int el = 0; el < nElmt; ++el)
            {
                elId = elIdArray[el];
                localNodeId = localNodeIdArray[el];
                GetFunctional<3>(derivUtil, nodes, elUtil, nodeMap, grad,
                            elId, localNodeId, ep, false,false);
            }
            Kokkos::deep_copy(grad.h_integral, grad.integral);
            newVal = grad.h_integral[0];
        // end evaluate node

            if(newVal <= currentW + c1() * (
                pg + 0.5*hes))
            {
                //this is a minimser so see if we can extend further
                while (l > -10)
                {
                // Update node
                    h_Xn[0] = h_Xc[0] + alpha * dk[0];
                    h_Xn[1] = h_Xc[1] + alpha * dk[1];
                    h_Xn[2] = h_Xc[2] + alpha * dk[2];
                    SetNodeCoord(h_Xn, globalNodeId, nodes, elIdArray, localNodeIdArray, nElmt);

                    grad.h_integral[0] = 0.0;
                    Kokkos::deep_copy(grad.integral,grad.h_integral);
                    for (int el = 0; el < nElmt; ++el)
                    {
                        GetFunctional<3>(derivUtil, nodes, elUtil, nodeMap, grad,
                            elId, localNodeId, ep, false,false);
                    }
                    Kokkos::deep_copy(grad.h_integral, grad.integral);
                    newVal = grad.h_integral[0];
                // end evaluate node

                // update node, but scaled
                    h_Xn[0] = h_Xc[0] + alpha/beta * dk[0];
                    h_Xn[1] = h_Xc[1] + alpha/beta * dk[1];
                    h_Xn[2] = h_Xc[2] + alpha/beta * dk[2];
                    SetNodeCoord(h_Xn, globalNodeId, nodes, elIdArray, localNodeIdArray, nElmt);

                    grad.h_integral[0] = 0.0;
                    Kokkos::deep_copy(grad.integral,grad.h_integral);
                    for (int el = 0; el < nElmt; ++el)
                    {
                        elId = elIdArray[el];
                        localNodeId = localNodeIdArray[el];
                        GetFunctional<3>(derivUtil, nodes, elUtil, nodeMap, grad,
                            elId, localNodeId, ep, false,false);
                    }
                    Kokkos::deep_copy(grad.h_integral, grad.integral);
                    dbVal = grad.h_integral[0];
                // end evaluate node

                    if (newVal <= currentW + c1() * (
                        alpha*pg + 0.5*alpha*alpha*hes) &&
                        dbVal > currentW + c1() *(
                        alpha/beta*pg + 0.5*alpha*alpha*hes/beta/beta))
                    {
                        found = true;
                        break;
                    }

                    l--;
                    alpha = pow(beta,l);
                }
            }
            else
            {
                //this is not a minimser so reverse line search
                while (alpha > alphaTol())
                {
                // Update node
                    h_Xn[0] = h_Xc[0] + alpha * dk[0];
                    h_Xn[1] = h_Xc[1] + alpha * dk[1];
                    h_Xn[2] = h_Xc[2] + alpha * dk[2];
                    SetNodeCoord(h_Xn, globalNodeId, nodes, elIdArray, localNodeIdArray, nElmt);

                    grad.h_integral[0] = 0.0;
                    Kokkos::deep_copy(grad.integral,grad.h_integral);
                    for (int el = 0; el < nElmt; ++el)
                    {
                        elId = elIdArray[el];
                        localNodeId = localNodeIdArray[el];
                        GetFunctional<3>(derivUtil, nodes, elUtil, nodeMap, grad,
                            elId, localNodeId, ep, false,false);
                    }
                    Kokkos::deep_copy(grad.h_integral, grad.integral);
                    newVal = grad.h_integral[0];
                // end evaluate node


                    if (newVal <= currentW + c1() * (
                        alpha*pg + 0.5*alpha*alpha*hes))
                    {
                        found = true;
                        break;
                    }

                    l++;
                    alpha = pow(beta,l);
                }
            }
        }*/

        if(!found)
        {
            h_Xn[0] = h_Xc[0];
            h_Xn[1] = h_Xc[1];
            h_Xn[2] = h_Xc[2];
            SetNodeCoord(h_Xn, globalNodeId, nodes, elIdArray, localNodeIdArray, nElmt);

            mtx.lock();
            res.nReset++;
            cout << "3d reset " << runDNC << endl;
            mtx.unlock();
        }

        // update residual value
        mtx.lock();        
        //GetNodeCoord(h_Xn, id, nodes, nodeMap);
        res.val = max(sqrt(  (h_Xn[0]-h_Xc[0])*(h_Xn[0]-h_Xc[0])
                            +(h_Xn[1]-h_Xc[1])*(h_Xn[1]-h_Xc[1])
                            +(h_Xn[2]-h_Xc[2])*(h_Xn[2]-h_Xc[2]) ),res.val );
        res.func +=newVal;
        mtx.unlock();
    }
    //});
}


/*NodeOptiJob* NodeOpti::GetJob()
{
    return new NodeOptiJob(this);
}*/


}
}
