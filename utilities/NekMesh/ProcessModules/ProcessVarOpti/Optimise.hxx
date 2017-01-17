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
#include "GetFunctional.hxx"

namespace Nektar
{
namespace Utilities
{


KOKKOS_INLINE_FUNCTION
double ProcessVarOpti::CalcMinJacGPU(const ElUtilGPU &elUtil, int nElmt, int node, int cs, Kokkos::View<int***> elIdArray)
{
    double minJac = DBL_MAX;
    int el;

    for(int i = 0; i < nElmt; i++)
    {
        el = elIdArray(cs,node,i);
        minJac = (minJac > elUtil.minJac(el) ? elUtil.minJac(el) : minJac);
            
    }
    return minJac;
}


KOKKOS_INLINE_FUNCTION
void ProcessVarOpti::GetNodeCoordGPU ( double (&X)[3], const NodesGPU &nodes,
            Kokkos::View<int***> elIdArray, Kokkos::View<int***> localNodeIdArray, int node, int cs)
{   
    int elmt = elIdArray(cs,node,0);
    int nodeId = localNodeIdArray(cs,node,0);
    
    X[0] = nodes.X(elmt,nodeId);
    X[1] = nodes.Y(elmt,nodeId);
    X[2] = nodes.Z(elmt,nodeId); 
}

KOKKOS_INLINE_FUNCTION
void ProcessVarOpti::SetNodeCoordGPU (const double (&X)[3], const NodesGPU &nodes,
             Kokkos::View<int***> elIdArray, Kokkos::View<int***> localNodeIdArray, int nElmt, int node, int cs)
{
    for(int el = 0; el < nElmt; el++)
    {
        int elmt = elIdArray(cs,node,el);
        int nodeId = localNodeIdArray(cs,node, el);

        nodes.X(elmt,nodeId) = X[0];
        nodes.Y(elmt,nodeId) = X[1];
        nodes.Z(elmt,nodeId) = X[2];        
    }   
}

KOKKOS_INLINE_FUNCTION
void ProcessVarOpti::GetNodeCoord(double (&X)[3], int id, NodesGPU &nodes,
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
void ProcessVarOpti::SetNodeCoord(double (&X)[3], int id, NodesGPU &nodes,
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


inline void CalcSK(double G[9], double sk[3])
{

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
}



void ProcessVarOpti::Optimise(DerivUtilGPU &derivUtil, NodesGPU &nodes, 
        ElUtilGPU &elUtil, Residual &res, int cs, optimiser opti)
{

    //printf("%s\n", "in OptimiseGPU");
    const int coloursetSize = nodes.h_coloursetSize(cs);
    Grad grad;
    grad.G = Kokkos::View<double*[9]> ("G",coloursetSize);
    grad.integral = Kokkos::View<double*> ("integral", coloursetSize);

    //int optiint = opti;
    //printf("optiint %i\n", optiint);
    //KOKKOS_INLINE_FUNCTION
    //__launch_bounds__ (32,24)
    Kokkos::parallel_for( team_policy( coloursetSize, 32 ), KOKKOS_LAMBDA ( const member_type& teamMember)
    {
        const int node = teamMember.league_rank();
        const int nElmt = nodes.nElmtArray(cs,node);        
        double currentW, newVal;

        double minJac = CalcMinJacGPU(elUtil, nElmt, node, cs, nodes.elIdArray);
        double ep = minJac < 0.0 ? sqrt(1e-9 + 0.04*minJac*minJac) : sqrt(1e-9); 
        
        currentW = GetFunctional<3>(derivUtil, nodes, elUtil,
                grad, nElmt, node, cs, ep, teamMember, opti, true);       

        double G[9];
        for (int i = 0; i < 9; ++i)
        {
            G[i] = grad.G(node,i);
        }
        
        /*if (node == 1)
        {
            double eval3[3];
            double eval2[2];
            CalcEValues<3>(G, eval3);
            double G2[4];
            for (int i = 0; i < 4; ++i)
            {
                G2[i] = G[i];
            }
            CalcEValues<2>(G2, eval2);
            IsIndefinite<3>(eval3);
            IsIndefinite<2>(eval2);
            double evec3[3];
            double evec2[2];
            CalcEVector<3>(G, eval3[2], evec3);
            CalcEVector<2>(G2, eval2[1], evec2);
        }*/        

        if(G[0]*G[0] + G[1]*G[1] + G[2]*G[2] > gradTol())
        {        
            double h_Xc[3];        
            double h_Xn[3];        
            GetNodeCoordGPU(h_Xc, nodes, nodes.elIdArray, nodes.localNodeIdArray, node, cs); 
        
            double sk[3];
            CalcSK(G, sk);
            double pg  = (G[0]*sk[0] + G[1]*sk[1] + G[2]*sk[2]);
            double lhs;
            double dk[3];
            
            bool runDNC = false; //so want to make this varible runDMC
            bool found  = false;

            double eval[3]; // the eigenvalues
            CalcEValues<3>(G, eval);
            
            int def = IsIndefinite<3>(eval);
            if(def)
            {               
                Kokkos::single(Kokkos::PerTeam(teamMember),[&] ()
                {   
                    printf("%s", "matrix is indefinite");
                    printf("lowest eigenvalue %16.16e\n", eval[2]);
                });
                //the dk vector needs calculating                
                CalcEVector<3>(G, eval[2], dk); //eval[2] is the minimum EigenValue

                if(dk[0]*G[0] + dk[1]*G[1] + dk[2]*G[2] > 0.0)
                {
                    for(int i = 0; i < 3; i++)
                    {
                        dk[i] *= -1.0;
                    }
                }

                lhs = dk[0] * (dk[0]*G[3] + dk[1]*G[4] + dk[2]*G[5]) +
                      dk[1] * (dk[0]*G[4] + dk[1]*G[6] + dk[2]*G[7]) +
                      dk[2] * (dk[0]*G[5] + dk[1]*G[7] + dk[2]*G[8]);

                ASSERTL0(lhs < 0.0 , "weirdness");

                double skmag = sqrt(sk[0]*sk[0] + sk[1]*sk[1] + sk[2]*sk[2]);
                runDNC =  !((G[0]*sk[0]+G[1]*sk[1]+G[2]*sk[2])/skmag <=
                            2.0*(0.5*lhs + G[0]*dk[0]+G[1]*dk[1]+G[2]*dk[2]));                
            }

            if(!runDNC)
            {
                //printf("%s\n", "run normal opti");
                double hes =sk[0] * (sk[0]*G[3] + sk[1]*G[4] + sk[2]*G[5]) +
                            sk[1] * (sk[0]*G[4] + sk[1]*G[6] + sk[2]*G[7]) +
                            sk[2] * (sk[0]*G[5] + sk[1]*G[7] + sk[2]*G[8]);
                hes = (hes > 0.0 ? 0.0 : hes);
                double alpha  = 1.0;
                
                //for(int c  = 0; c < 34; c++)
                while (alpha > alphaTol())
                {
                // Update node                
                    h_Xn[0] = h_Xc[0] + alpha * sk[0];
                    h_Xn[1] = h_Xc[1] + alpha * sk[1];
                    h_Xn[2] = h_Xc[2] + alpha * sk[2];
                    SetNodeCoordGPU(h_Xn, nodes, nodes.elIdArray, 
                            nodes.localNodeIdArray, nElmt, node, cs);
                    
                    newVal = GetFunctional<3>(derivUtil, nodes, elUtil, grad, 
                            nElmt, node, cs, ep, teamMember, opti, false);
               
                    if (newVal <= currentW + c1() * (alpha*pg+ 0.5*alpha*alpha*hes))
                    {
                        found = true;
                        break;
                    }
                    alpha /= 2.0;        
                }
            }
            else
            {
                //printf("%s\n", "run DNC opti");
                double beta = 0.5;
                int l = 0;
                double alpha = pow(beta,l);

                double hes = lhs;

                pg = (G[0]*dk[0]+G[1]*dk[1]+G[2]*dk[2]);

                //choose whether to do forward or reverse line search
                h_Xn[0] = h_Xc[0] + dk[0];
                h_Xn[1] = h_Xc[1] + dk[1];
                h_Xn[2] = h_Xc[2] + dk[2];
                SetNodeCoordGPU(h_Xn, nodes, nodes.elIdArray, 
                        nodes.localNodeIdArray, nElmt, node, cs);

                newVal = GetFunctional<3>(derivUtil, nodes, elUtil, grad, 
                        nElmt, node, cs, ep, teamMember, opti, false);

                if(newVal <= currentW + c1() * (pg + 0.5*hes))
                {
                    //this is a minimser so see if we can extend further
                    while (l > -10)
                    {
                        // Update node
                        h_Xn[0] = h_Xc[0] + alpha * dk[0];
                        h_Xn[1] = h_Xc[1] + alpha * dk[1];
                        h_Xn[2] = h_Xc[2] + alpha * dk[2];
                        SetNodeCoordGPU(h_Xn, nodes, nodes.elIdArray, 
                                nodes.localNodeIdArray, nElmt, node, cs);

                        newVal = GetFunctional<3>(derivUtil, nodes, elUtil, grad, 
                                nElmt, node, cs, ep, teamMember, opti, false);

                        h_Xn[0] = h_Xc[0] + alpha/beta * dk[0];
                        h_Xn[1] = h_Xc[1] + alpha/beta * dk[1];
                        h_Xn[2] = h_Xc[2] + alpha/beta * dk[2];
                        SetNodeCoordGPU(h_Xn, nodes, nodes.elIdArray, 
                                nodes.localNodeIdArray, nElmt, node, cs);

                        double dbVal = GetFunctional<3>(derivUtil, nodes, elUtil, grad, 
                                nElmt, node, cs, ep, teamMember, opti, false);

                        if (newVal <= currentW + c1() * (
                            alpha*pg + 0.5*alpha*alpha*hes) &&
                            dbVal > currentW + c1() * (
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
                        SetNodeCoordGPU(h_Xn, nodes, nodes.elIdArray, 
                                nodes.localNodeIdArray, nElmt, node, cs);

                        newVal = GetFunctional<3>(derivUtil, nodes, elUtil, grad, 
                                nElmt, node, cs, ep, teamMember, opti, false);

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
            }       

            if(!found)
            {
                h_Xn[0] = h_Xc[0];
                h_Xn[1] = h_Xc[1];
                h_Xn[2] = h_Xc[2];
                SetNodeCoordGPU(h_Xn, nodes, nodes.elIdArray, 
                        nodes.localNodeIdArray, nElmt, node, cs);

                Kokkos::single(Kokkos::PerTeam(teamMember),[&] ()
                {
                    Kokkos::atomic_add(&res.nReset[0], 1); 
                    printf("%s\n", "3d reset");                   
                });   
            }

            // store residual values of each node
            double resid = sqrt( (h_Xn[0]-h_Xc[0])*(h_Xn[0]-h_Xc[0])
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



} // Utilities

} // Nektar

# endif