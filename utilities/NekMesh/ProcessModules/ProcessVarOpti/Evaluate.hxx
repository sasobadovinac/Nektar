////////////////////////////////////////////////////////////////////////////////
//
//  File: Evaluate.hxx
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

#ifndef UTILITIES_NEKMESH_PROCESSVAROPTI_EVALUATE
#define UTILITIES_NEKMESH_PROCESSVAROPTI_EVALUATE


namespace Nektar
{
namespace Utilities
{

void ProcessVarOpti::Evaluate(DerivUtilGPU &derivUtil,NodesGPU &nodes, ElUtilGPU &elUtil, Residual &res)
{
    int m_dim = dataSet[0]->m_dim;
    if(m_dim == 2)
    {
        ASSERTL0(false,"not coded");
    }
    else if (m_dim == 3)
    {
        int nodes_size = nodes.nodes_size;
        int nElmt = nodes.globalnElmt;   

        // declare and initialise min and max Jacobian
        Kokkos::View<double*> mx("mx",nElmt);
        Kokkos::View<double*> mn("mn",nElmt);
        Kokkos::parallel_for(range_policy(0,nElmt), KOKKOS_LAMBDA (const int el)
        {
            mx(el) = -1.0 * DBL_MAX;
            mn(el) = DBL_MAX;
        });
        
        // do the matrix vector multiplication on the GPU
        Kokkos::View<double[3][3]> dxdz("dxdz");
        Kokkos::View<double**> jacDet("jacDet", nElmt, nodes_size);
        
        Kokkos::parallel_for( team_policy( nElmt , Kokkos::AUTO ), KOKKOS_LAMBDA ( const member_type& teamMember)
        {
            const int el = teamMember.league_rank();
        
            Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , nodes_size), [&] ( const int i)
            {
                double x1i = 0.0;
                double y1i = 0.0;
                double z1i = 0.0;
                double x2i = 0.0;
                double y2i = 0.0;
                double z2i = 0.0;
                double x3i = 0.0;
                double y3i = 0.0;
                double z3i = 0.0;
                            
                for (int j = 0; j < nodes_size; ++j)
                {
                    //x1i += derivUtil.VdmDL_0( i , j ) * nodes.X( el , j );
                    y1i += derivUtil.VdmDL_0( i , j ) * nodes.Y( el , j );
                    z1i += derivUtil.VdmDL_0( i , j ) * nodes.Z( el , j );
                    x2i += derivUtil.VdmDL_1( i , j ) * nodes.X( el , j );
                    y2i += derivUtil.VdmDL_1( i , j ) * nodes.Y( el , j );
                    z2i += derivUtil.VdmDL_1( i , j ) * nodes.Z( el , j );
                    x3i += derivUtil.VdmDL_2( i , j ) * nodes.X( el , j );
                    y3i += derivUtil.VdmDL_2( i , j ) * nodes.Y( el , j );
                    z3i += derivUtil.VdmDL_2( i , j ) * nodes.Z( el , j );
                }
                  
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, nodes_size), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_0( i , j ) * nodes.X( el , j );
                }, x1i);
                /*Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, nodes_size), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_0( i , j ) * Y( el , j );
                }, y1i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, nodes_size), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_0( i , j ) * Z( el , j );
                }, z1i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, nodes_size), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_1( i , j ) * X( el , j );
                }, x2i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, nodes_size), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_1( i , j ) * Y( el , j );
                }, y2i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, nodes_size), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_1( i , j ) * Z( el , j );
                }, z2i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, nodes_size), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_2( i , j ) * X( el , j );
                }, x3i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, nodes_size), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_2( i , j ) * Y( el , j );
                }, y3i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, nodes_size), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_2( i , j ) * Z( el , j );
                }, z3i);*/
                  
                // calculate the Jacobian
                dxdz(0,0) = x1i;
                dxdz(0,1) = y1i;
                dxdz(0,2) = z1i;
                dxdz(1,0) = x2i;
                dxdz(1,1) = y2i;
                dxdz(1,2) = z2i;
                dxdz(2,0) = x3i;
                dxdz(2,1) = y3i;
                dxdz(2,2) = z3i;

                jacDet(el,i) = dxdz(0,0)*(dxdz(1,1)*dxdz(2,2)-dxdz(2,1)*dxdz(1,2))
                                    -dxdz(0,1)*(dxdz(1,0)*dxdz(2,2)-dxdz(2,0)*dxdz(1,2))
                                    +dxdz(0,2)*(dxdz(1,0)*dxdz(2,1)-dxdz(2,0)*dxdz(1,1));           
            });

        });
        
        // initialisation
        Kokkos::parallel_for("output", range_policy(0,1), KOKKOS_LAMBDA (const int& el)
        {
            res.worstJac[0] = DBL_MAX;
            res.startInv[0] = 0;
        });
        
        // compute minimum and scaled Jacobian of each element 
        Kokkos::parallel_for(range_policy(0,nElmt), KOKKOS_LAMBDA (const int el)
        {
            for(int j = 0; j < nodes_size; j++)
            {
                //  mx = max(mx,jacDet);
                mx(el) = (mx(el) < jacDet(el,j) ? jacDet(el,j) : mx(el));
                //  mn = min(mn,jacDet);
                mn(el) = (mn(el) > jacDet(el,j) ? jacDet(el,j) : mn(el));
            }
            //dataSet[el]->minJac = h_mn(el);
            //elUtil.minJac(el) = mn(el);            
            //dataSet[el]->scaledJac = h_mn(el)/h_mx(el);
            elUtil.scaledJac(el) = mn(el)/mx(el);

            // compute the smallest (worst) Jacobian of all elements
            Kokkos::atomic_fetch_min(&res.worstJac[0], elUtil.scaledJac(el));

            // compute number of invalied elements (Jacobian < 0)            
            if(mn(el) < 0)
            {
                Kokkos::atomic_add(&res.startInv[0],1);
            }
        });

        // compute the smallest (worst) Jacobian of all elements
        // doesnt work yet due to Kokkos bug !!!!

        /*Kokkos::parallel_for( team_policy( 1 , Kokkos::AUTO ), KOKKOS_LAMBDA ( const member_type& teamMember)
        {  
            double JacJac;
            MinFunctor <double> mnfunctor(elUtil.scaledJac);
            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(teamMember, nElmt) , mnfunctor, JacJac);//res.worstJac[0]);
            Kokkos::single( Kokkos::PerTeam(teamMember), [&] ()
            {
                printf("Worst JacJac: %e  ", JacJac);
            });           
        });
        MinFunctor <double> mnfunctor(elUtil.scaledJac);
        Kokkos::parallel_reduce(range_policy(0, nElmt) , mnfunctor, res.h_worstJac[0]);
        printf("Worst Jacobian: %e  ", res.h_worstJac[0]);*/
		
        Kokkos::parallel_for("output", range_policy(0,1), KOKKOS_LAMBDA (const int& el)
        {
            // print the residual as calculated in the OptimiseGPU function
            printf("Residual: %e  ", res.val[0]);

            // print worst Jacobian and number of invalid elements
            printf("Worst Jacobian: %e  ", res.worstJac[0]);            
            printf("Invalid Elements: %i  ", res.startInv[0]);
            
            // print and reset functional and reset-nodes for next run of OptimiseGPU
            printf("Reset Nodes: %i  ", res.nReset[0]);
            printf("Functional: %e\n", res.func[0]);
            res.func[0] = 0.0;
            res.nReset[0] = 0;
        });
        
    }
}



void ProcessVarOpti::InitialMinJac(DerivUtilGPU &derivUtil,NodesGPU &nodes)
{
    int m_dim = dataSet[0]->m_dim;
    if(m_dim == 2)
    {
        ASSERTL0(false,"not coded");
    }
    else if (m_dim == 3)
    {
        int nodes_size = nodes.nodes_size;
        int nElmt = nodes.globalnElmt;
        int ptsHigh = derivUtil.ptsHigh;   

        // declare and initialise min and max Jacobian
        Kokkos::View<double*> mn("mn",nElmt);
        Kokkos::parallel_for(range_policy(0,nElmt), KOKKOS_LAMBDA (const int el)
        {
            mn(el) = DBL_MAX;
        });
        
        // do the matrix vector multiplication on the GPU
        Kokkos::View<double**> jacDet("jacDet", nElmt, ptsHigh);
        
        Kokkos::parallel_for( team_policy( nElmt , Kokkos::AUTO ), KOKKOS_LAMBDA ( const member_type& teamMember)
        {
            const int el = teamMember.league_rank();
        
            Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , ptsHigh), [&] ( const int i)
            {
                double x1i = 0.0;
                double y1i = 0.0;
                double z1i = 0.0;
                double x2i = 0.0;
                double y2i = 0.0;
                double z2i = 0.0;
                double x3i = 0.0;
                double y3i = 0.0;
                double z3i = 0.0;
                            
                for (int j = 0; j < nodes_size; ++j)
                {
                    x1i += derivUtil.VdmD_0( i , j ) * nodes.X( el , j );
                    y1i += derivUtil.VdmD_0( i , j ) * nodes.Y( el , j );
                    z1i += derivUtil.VdmD_0( i , j ) * nodes.Z( el , j );
                    x2i += derivUtil.VdmD_1( i , j ) * nodes.X( el , j );
                    y2i += derivUtil.VdmD_1( i , j ) * nodes.Y( el , j );
                    z2i += derivUtil.VdmD_1( i , j ) * nodes.Z( el , j );
                    x3i += derivUtil.VdmD_2( i , j ) * nodes.X( el , j );
                    y3i += derivUtil.VdmD_2( i , j ) * nodes.Y( el , j );
                    z3i += derivUtil.VdmD_2( i , j ) * nodes.Z( el , j );
                }
                // calculate the Jacobian  
                double dxdz[3][3];
                dxdz[0][0] = x1i;
                dxdz[0][1] = y1i;
                dxdz[0][2] = z1i;
                dxdz[1][0] = x2i;
                dxdz[1][1] = y2i;
                dxdz[1][2] = z2i;
                dxdz[2][0] = x3i;
                dxdz[2][1] = y3i;
                dxdz[2][2] = z3i;
                jacDet(el,i) = dxdz[0][0]*(dxdz[1][1]*dxdz[2][2]-dxdz[2][1]*dxdz[1][2])
                                    -dxdz[0][1]*(dxdz[1][0]*dxdz[2][2]-dxdz[2][0]*dxdz[1][2])
                                    +dxdz[0][2]*(dxdz[1][0]*dxdz[2][1]-dxdz[2][0]*dxdz[1][1]);
                           
            });

            for(int i = 0; i < ptsHigh; i++)
            {
                //  mn = min(mn,jacDet);
                mn(el) = (mn(el) > jacDet(el,i) ? jacDet(el,i) : mn(el));
            }

            // copy the minimum Jacobian of each element to each node of that element
            for (int j = 0; j < nodes_size; j++)
            {
                nodes.minJac(el,j) = mn(el);
            }            


        });

    }
         
}



} // Utilities

} // Nektar

# endif