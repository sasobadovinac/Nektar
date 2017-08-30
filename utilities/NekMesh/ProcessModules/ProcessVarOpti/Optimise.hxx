////////////////////////////////////////////////////////////////////////////////
//
//  File: OptimiseGPU.hxx
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

#ifndef UTILITIES_NEKMESH_PROCESSVAROPTI_OPTMISE
#define UTILITIES_NEKMESH_PROCESSVAROPTI_OPTMISE

#include "Hessian.hxx"
#include "GetFunctional_eHypEl.hxx"
#include "GetFunctional_eLinEl.hxx"
#include "GetFunctional_eRoca.hxx"
#include "GetFunctional_eWins.hxx"

namespace Nektar
{
namespace Utilities
{


KOKKOS_INLINE_FUNCTION
void ProcessVarOpti::CalcMinJacGPU (const NodesGPU &nodes, int nElmt, int node, int cs)
{
    double minJac = DBL_MAX;
    for(int el = 0; el < nElmt; el++)
    {
        int elmt = nodes.elIdArray(cs, node, el);
        int nodeId = nodes.localNodeIdArray(cs, node, el);

        minJac = (minJac > nodes.minJac(elmt,nodeId) ? nodes.minJac(elmt,nodeId) : minJac);        
    }
    for(int el = 0; el < nElmt; el++)
    {
        int elmt = nodes.elIdArray(cs, node, el);
        int nodeId = nodes.localNodeIdArray(cs, node, el);

        nodes.minJac(elmt,nodeId) = minJac;        
    }  
}

KOKKOS_INLINE_FUNCTION
double ProcessVarOpti::GetMinJacGPU (const NodesGPU &nodes, int nElmt, int node, int cs)
{
    
    int elmt = nodes.elIdArray(cs,node,0);
    int nodeId = nodes.localNodeIdArray(cs,node, 0);

    return nodes.minJac(elmt,nodeId);     
       
}

KOKKOS_INLINE_FUNCTION
void ProcessVarOpti::SetMinJacGPU (const Grad &grad, const NodesGPU &nodes, int nElmt, int node, int cs)
{
    for(int el = 0; el < nElmt; el++)
    {
        int elmt = nodes.elIdArray(cs,node,el);
        int nodeId = nodes.localNodeIdArray(cs,node, el);

        nodes.minJac(elmt,nodeId) = grad.minJacNew(node);     
    }   
}

template<int DIM>
KOKKOS_INLINE_FUNCTION
void ProcessVarOpti::GetNodeCoordGPU ( double (&X)[DIM], const NodesGPU &nodes, int node, int cs)
{
}
template<>
KOKKOS_INLINE_FUNCTION
void ProcessVarOpti::GetNodeCoordGPU<2> ( double (&X)[2], const NodesGPU &nodes, int node, int cs)
{   
    int elmt = nodes.elIdArray(cs,node,0);
    int nodeId = nodes.localNodeIdArray(cs,node,0);
    
    X[0] = nodes.X(elmt,nodeId);
    X[1] = nodes.Y(elmt,nodeId); 
}
template<>
KOKKOS_INLINE_FUNCTION
void ProcessVarOpti::GetNodeCoordGPU<3> ( double (&X)[3], const NodesGPU &nodes, int node, int cs)
{   
    int elmt = nodes.elIdArray(cs,node,0);
    int nodeId = nodes.localNodeIdArray(cs,node,0);
    
    X[0] = nodes.X(elmt,nodeId);
    X[1] = nodes.Y(elmt,nodeId);
    X[2] = nodes.Z(elmt,nodeId); 
}

template<int DIM>
KOKKOS_INLINE_FUNCTION
void ProcessVarOpti::SetNodeCoordGPU (const double (&X)[DIM], const NodesGPU &nodes, int nElmt, int node, int cs)
{
}
template<>
KOKKOS_INLINE_FUNCTION
void ProcessVarOpti::SetNodeCoordGPU<2> (const double (&X)[2], const NodesGPU &nodes, int nElmt, int node, int cs)
{
    for(int el = 0; el < nElmt; el++)
    {
        int elmt = nodes.elIdArray(cs,node,el);
        int nodeId = nodes.localNodeIdArray(cs,node, el);

        nodes.X(elmt,nodeId) = X[0];
        nodes.Y(elmt,nodeId) = X[1];       
    }   
}
template<>
KOKKOS_INLINE_FUNCTION
void ProcessVarOpti::SetNodeCoordGPU<3> (const double (&X)[3], const NodesGPU &nodes, int nElmt, int node, int cs)
{
    for(int el = 0; el < nElmt; el++)
    {
        int elmt = nodes.elIdArray(cs,node,el);
        int nodeId = nodes.localNodeIdArray(cs,node, el);

        nodes.X(elmt,nodeId) = X[0];
        nodes.Y(elmt,nodeId) = X[1];
        nodes.Z(elmt,nodeId) = X[2];        
    }   
}


void ProcessVarOpti::Optimise3D3D(DerivUtilGPU &derivUtil, NodesGPU &nodes, 
        ElUtilGPU &elUtil, Residual &res, int cs, optimiser opti)
{

    //printf("%s\n", "in OptimiseGPU");
    const int coloursetSize = nodes.h_coloursetSize(cs);
    Grad grad;
    grad.G = Kokkos::View<double*[9]> ("G",coloursetSize);
    grad.integral = Kokkos::View<double*> ("integral", coloursetSize);
    grad.minJacNew = Kokkos::View<double*> ("minJacNew", coloursetSize);

    int scratch_size = ScratchViewType::shmem_size(3*derivUtil.ptsLow); // for all three dimensions

    Kokkos::parallel_for( team_policy( coloursetSize, Kokkos::AUTO).set_scratch_size(0,Kokkos::PerTeam(scratch_size))
        , KOKKOS_LAMBDA ( const member_type& teamMember)
    {
        const int node = teamMember.league_rank();
        const int nElmt = nodes.nElmtArray(cs,node);        
        double currentW, newVal;

        double minJac = GetMinJacGPU(nodes, nElmt, node, cs);
        double ep = minJac < 0.0 ? sqrt(1e-8 + 0.04*minJac*minJac) : 1e-4; 
        
        // Timings
        /*#ifdef __CUDA_ARCH__
        long long int start = clock64();
        #endif*/
        
        if (opti == eHypEl)
        {   
            currentW = GetFunctional<3,true,eHypEl>()(derivUtil, nodes, elUtil,
                    grad, nElmt, node, cs, ep, teamMember);
        }
        else if (opti == eLinEl)
        {   
            currentW = GetFunctional<3,true,eLinEl>()(derivUtil, nodes, elUtil,
                    grad, nElmt, node, cs, ep, teamMember);
        }
        else if (opti == eRoca)
        {   
            currentW = GetFunctional<3,true,eRoca>()(derivUtil, nodes, elUtil,
                    grad, nElmt, node, cs, ep, teamMember);
        }
        else if (opti == eWins)
        {   
            currentW = GetFunctional<3,true,eWins>()(derivUtil, nodes, elUtil,
                    grad, nElmt, node, cs, ep, teamMember);
        }
        
        newVal = currentW;

        /*#ifdef __CUDA_ARCH__
        long long int stop = clock64();
        long long int cycles = stop - start;
        int thread = threadIdx.y;
        //if (thread ==0){ // its the same for the whole block!
        printf("GPUTimer Node %i: %lli\n", node, cycles); //}
        #endif*/


        double G[9];
        for (int i = 0; i < 9; ++i)
        {
            G[i] = grad.G(node,i);
        }                

        if(G[0]*G[0] + G[1]*G[1] + G[2]*G[2] > gradTol())
        {        
            double h_Xc[3];        // current
            double h_Xn[3];        // new       
            GetNodeCoordGPU<3>(h_Xc, nodes, node, cs);
           
            double sk[3];
            double eval[3]; // the eigenvalues
            CalcEValues<3>(G, eval); //eval[2] is the minimum EigenValue
            
            if(eval[2] < 1e-6)
            {
                G[3] += 1e-6 - eval[2];
                G[6] += 1e-6 - eval[2];
                G[8] += 1e-6 - eval[2];
            }
            
            //calculate sk
            NekDouble det = G[3]*(G[6]*G[8]-G[7]*G[7])
                   -G[4]*(G[4]*G[8]-G[5]*G[7])
                   +G[5]*(G[4]*G[7]-G[5]*G[6]);

            sk[0] = G[0]*(G[6]*G[8]-G[7]*G[7]) +
                    G[1]*(G[5]*G[7]-G[4]*G[8]) +
                    G[2]*(G[4]*G[7]-G[3]*G[7]);
            sk[1] = G[0]*(G[7]*G[5]-G[4]*G[5]) +
                    G[1]*(G[3]*G[8]-G[5]*G[5]) +
                    G[2]*(G[4]*G[5]-G[3]*G[7]);
            sk[2] = G[0]*(G[4]*G[7]-G[6]*G[5]) +
                    G[1]*(G[4]*G[5]-G[3]*G[7]) +
                    G[2]*(G[3]*G[6]-G[4]*G[4]);

            sk[0] /= det * -1.0;
            sk[1] /= det * -1.0;
            sk[2] /= det * -1.0;

            bool found  = false;

            double pg  = (G[0]*sk[0] + G[1]*sk[1] + G[2]*sk[2]);
            
            double alpha  = 1.0;
            
            while (alpha > alphaTol())
            {
                // Update node                
                h_Xn[0] = h_Xc[0] + alpha * sk[0];
                h_Xn[1] = h_Xc[1] + alpha * sk[1];
                h_Xn[2] = h_Xc[2] + alpha * sk[2];
                SetNodeCoordGPU<3>(h_Xn, nodes, nElmt, node, cs);
                
                //newVal = GetFunctional<3,false,eHypEl>()(derivUtil, nodes, elUtil, grad, 
                //        nElmt, node, cs, ep, teamMember);
                if (opti == eHypEl)
                {   
                    newVal = GetFunctional<3,false,eHypEl>()(derivUtil, nodes, elUtil,
                            grad, nElmt, node, cs, ep, teamMember);
                }
                else if (opti == eLinEl)
                {   
                    newVal = GetFunctional<3,false,eLinEl>()(derivUtil, nodes, elUtil,
                            grad, nElmt, node, cs, ep, teamMember);
                }
                else if (opti == eRoca)
                {   
                    newVal = GetFunctional<3,false,eRoca>()(derivUtil, nodes, elUtil,
                            grad, nElmt, node, cs, ep, teamMember);
                }
                else if (opti == eWins)
                {   
                    newVal = GetFunctional<3,false,eWins>()(derivUtil, nodes, elUtil,
                            grad, nElmt, node, cs, ep, teamMember);
                }

           
                if (newVal <= currentW + c1() * alpha * pg)
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
                SetNodeCoordGPU<3>(h_Xn, nodes, nElmt, node, cs);

                Kokkos::single(Kokkos::PerTeam(teamMember),[&] ()
                {
                    Kokkos::atomic_add(&res.nReset[0], 1); 
                    //printf("%s\n", "3d reset");                   
                });   
            }
            else
            {
                SetMinJacGPU(grad, nodes, nElmt, node, cs);
            }

            // store residual values of each node
            double resid = sqrt(   (h_Xn[0]-h_Xc[0])*(h_Xn[0]-h_Xc[0])
                                  +(h_Xn[1]-h_Xc[1])*(h_Xn[1]-h_Xc[1])
                                  +(h_Xn[2]-h_Xc[2])*(h_Xn[2]-h_Xc[2]) );            

            Kokkos::single(Kokkos::PerTeam(teamMember),[&] ()
            {
                Kokkos::atomic_add(&res.func[0], newVal);
                Kokkos::atomic_fetch_max(&res.val[0], resid);
            });            

        }        
        //printf("node %i finished\n", node);
    
    });
    //printf("colorset %i finished\n", cs);

}



void ProcessVarOpti::Optimise2D2D(DerivUtilGPU &derivUtil, NodesGPU &nodes, 
        ElUtilGPU &elUtil, Residual &res, int cs, optimiser opti)
{
    const int coloursetSize = nodes.h_coloursetSize(cs);
    Grad grad;
    grad.G = Kokkos::View<double*[9]> ("G",coloursetSize);
    grad.integral = Kokkos::View<double*> ("integral", coloursetSize);
    grad.minJacNew = Kokkos::View<double*> ("minJacNew", coloursetSize);

    int scratch_size = ScratchViewType::shmem_size(2*derivUtil.ptsLow); // for all two dimensions

    Kokkos::parallel_for( team_policy( coloursetSize, Kokkos::AUTO ).set_scratch_size(0,Kokkos::PerTeam(scratch_size))
        , KOKKOS_LAMBDA ( const member_type& teamMember)
    {
        const int node = teamMember.league_rank();
        const int nElmt = nodes.nElmtArray(cs,node);        
        double currentW, newVal;

        double minJac = GetMinJacGPU(nodes, nElmt, node, cs);
        double ep = minJac < 0.0 ? sqrt(1e-8 + 0.04*minJac*minJac) : 1e-4; 
        
        //currentW = GetFunctional<2,true,eHypEl>()(derivUtil, nodes, elUtil,
        //        grad, nElmt, node, cs, ep, teamMember);
        if (opti == eHypEl)
        {   
            currentW = GetFunctional<2,true,eHypEl>()(derivUtil, nodes, elUtil,
                    grad, nElmt, node, cs, ep, teamMember);
        }
        else if (opti == eLinEl)
        {   
            currentW = GetFunctional<2,true,eLinEl>()(derivUtil, nodes, elUtil,
                    grad, nElmt, node, cs, ep, teamMember);
        }
        else if (opti == eRoca)
        {   
            currentW = GetFunctional<2,true,eRoca>()(derivUtil, nodes, elUtil,
                    grad, nElmt, node, cs, ep, teamMember);
        }
        else if (opti == eWins)
        {   
            currentW = GetFunctional<2,true,eWins>()(derivUtil, nodes, elUtil,
                    grad, nElmt, node, cs, ep, teamMember);
        }

        newVal = currentW;

        double G[4];
        for (int i = 0; i < 4; ++i)
        {
            G[i] = grad.G(node,i);
        }    

        if(G[0]*G[0] + G[1]*G[1] > gradTol())
        {        
            double h_Xc[2];        // current
            double h_Xn[2];        // new       
            GetNodeCoordGPU<2>(h_Xc, nodes, node, cs);
           
            double sk[2];
            double eval[2]; // the eigenvalues
            CalcEValues<2>(G, eval); //eval[1] is the minimum EigenValue
            
            if(eval[1] < 1e-6)
            {
                G[2] += 1e-6 - eval[1];
                G[4] += 1e-6 - eval[1];
            }
            
            // calculate sk
            sk[0] = -1.0 / (G[2] * G[4] - G[3] * G[3]) *
                (G[4] * G[0] - G[3] * G[1]);
            sk[1] = -1.0 / (G[2] * G[4] - G[3] * G[3]) *
                (G[2] * G[1] - G[3] * G[0]);

            bool found  = false;

            double pg  = (G[0]*sk[0] + G[1]*sk[1]);
            
            double alpha  = 1.0;
            
            while (alpha > alphaTol())
            {
                // Update node                
                h_Xn[0] = h_Xc[0] + alpha * sk[0];
                h_Xn[1] = h_Xc[1] + alpha * sk[1];
                SetNodeCoordGPU<2>(h_Xn, nodes, nElmt, node, cs);
                
                //newVal = GetFunctional<2,false,eHypEl>()(derivUtil, nodes, elUtil, grad, 
                //        nElmt, node, cs, ep, teamMember);
                if (opti == eHypEl)
                {   
                    newVal = GetFunctional<2,false,eHypEl>()(derivUtil, nodes, elUtil,
                            grad, nElmt, node, cs, ep, teamMember);
                }
                else if (opti == eLinEl)
                {   
                    newVal = GetFunctional<2,false,eLinEl>()(derivUtil, nodes, elUtil,
                            grad, nElmt, node, cs, ep, teamMember);
                }
                else if (opti == eRoca)
                {   
                    newVal = GetFunctional<2,false,eRoca>()(derivUtil, nodes, elUtil,
                            grad, nElmt, node, cs, ep, teamMember);
                }
                else if (opti == eWins)
                {   
                    newVal = GetFunctional<2,false,eWins>()(derivUtil, nodes, elUtil,
                            grad, nElmt, node, cs, ep, teamMember);
                }
           
                if (newVal <= currentW + c1() * alpha * pg)
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
                SetNodeCoordGPU<2>(h_Xn, nodes, nElmt, node, cs);

                Kokkos::single(Kokkos::PerTeam(teamMember),[&] ()
                {
                    Kokkos::atomic_add(&res.nReset[0], 1); 
                    //printf("%s\n", "3d reset");                   
                });   
            }
            else
            {
                SetMinJacGPU(grad, nodes, nElmt, node, cs);
            }

            // store residual values of each node
            double resid = sqrt(   (h_Xn[0]-h_Xc[0])*(h_Xn[0]-h_Xc[0])
                                  +(h_Xn[1]-h_Xc[1])*(h_Xn[1]-h_Xc[1]) );            

            Kokkos::single(Kokkos::PerTeam(teamMember),[&] ()
            {
                Kokkos::atomic_add(&res.func[0], newVal);
                Kokkos::atomic_fetch_max(&res.val[0], resid);
            });            

        }        
        //printf("node %i finished\n", node);
    
    });
    //printf("colorset %i finished\n", cs);

}





} // Utilities

} // Nektar

# endif