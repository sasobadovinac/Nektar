////////////////////////////////////////////////////////////////////////////////
//
//  File: EvaluateGPU.hxx
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

#ifndef UTILITIES_NEKMESH_PROCESSVAROPTI_EVALUATEGPU
#define UTILITIES_NEKMESH_PROCESSVAROPTI_EVALUATEGPU

#include <map>

namespace Nektar
{
namespace Utilities
{

boost::mutex mtx3;

void ProcessVarOpti::Load_derivUtil(DerivUtilGPU &derivUtil)
{
    //typedef Kokkos::MemoryTraits<Kokkos::RandomAccess> random_memory; //?
    int nodes_size = derivUtil.nodes_size;
    int ptsHigh = derivUtil.ptsHigh;
    derivUtil.VdmDL_0 = Kokkos::View<double**>("VdmDL_0",nodes_size,nodes_size);
    derivUtil.h_VdmDL_0 = Kokkos::create_mirror_view(derivUtil.VdmDL_0);
    derivUtil.VdmDL_1 = Kokkos::View<double**>("VdmDL_1",nodes_size,nodes_size);
    derivUtil.h_VdmDL_1 = Kokkos::create_mirror_view(derivUtil.VdmDL_1);
    int N1 = nodes_size;
    int M1 = nodes_size;
    Kokkos::parallel_for(team_policy_host(N1,M1), KOKKOS_LAMBDA (const member_type_host& teamMember)
    {         
        const int i = teamMember.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, M1 ), [&] (const int j)
        {                
            derivUtil.h_VdmDL_0(i,j) = dataSet[0]->derivUtil->VdmDL[0](i,j);
            derivUtil.h_VdmDL_1(i,j) = dataSet[0]->derivUtil->VdmDL[1](i,j);
        }); 
    });
    Kokkos::deep_copy(derivUtil.VdmDL_0,derivUtil.h_VdmDL_0);
    Kokkos::deep_copy(derivUtil.VdmDL_1,derivUtil.h_VdmDL_1);
    

    derivUtil.VdmD_0 = Kokkos::View<double**>("VdmD_0",ptsHigh,nodes_size);
    derivUtil.h_VdmD_0 = Kokkos::create_mirror_view(derivUtil.VdmD_0);
    derivUtil.VdmD_1 = Kokkos::View<double**>("VdmD_1",ptsHigh,nodes_size);
    derivUtil.h_VdmD_1 = Kokkos::create_mirror_view(derivUtil.VdmD_1);   
    int N2 = ptsHigh;
    int M2 = nodes_size;
    Kokkos::parallel_for(team_policy_host(N2,M2), KOKKOS_LAMBDA (const member_type_host& teamMember)
    {         
        const int i = teamMember.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, M2 ), [&] (const int j)
        {                
            derivUtil.h_VdmD_0(i,j) = dataSet[0]->derivUtil->VdmD[0](i,j);
            derivUtil.h_VdmD_1(i,j) = dataSet[0]->derivUtil->VdmD[1](i,j);                      
        }); 
    });   
    Kokkos::deep_copy(derivUtil.VdmD_0,derivUtil.h_VdmD_0);
    Kokkos::deep_copy(derivUtil.VdmD_1,derivUtil.h_VdmD_1);

    
    int m_dim = dataSet[0]->m_dim;
    if(m_dim == 3)
    {
        derivUtil.VdmDL_2 = Kokkos::View<double**>("VdmDL_2",nodes_size,nodes_size);
        derivUtil.h_VdmDL_2 = Kokkos::create_mirror_view(derivUtil.VdmDL_2);
        int N1 = nodes_size;
	    int M1 = nodes_size;
	    Kokkos::parallel_for(team_policy_host(N1,M1), KOKKOS_LAMBDA (const member_type_host& teamMember)
	    {         
	        const int i = teamMember.league_rank();
	        Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, M1 ), [&] (const int j)
	        {                
	            derivUtil.h_VdmDL_2(i,j) = dataSet[0]->derivUtil->VdmDL[2](i,j);
	        }); 
	    });
	    Kokkos::deep_copy(derivUtil.VdmDL_2,derivUtil.h_VdmDL_2);


        derivUtil.VdmD_2 = Kokkos::View<double**>("VdmD_2",ptsHigh,nodes_size);
        derivUtil.h_VdmD_2 = Kokkos::create_mirror_view(derivUtil.VdmD_2);
        int N2 = ptsHigh;
	    int M2 = nodes_size;
	    Kokkos::parallel_for(team_policy_host(N2,M2), KOKKOS_LAMBDA (const member_type_host& teamMember)
	    {         
	        const int i = teamMember.league_rank();
	        Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, M2 ), [&] (const int j)
	        {                
	            derivUtil.h_VdmD_2(i,j) = dataSet[0]->derivUtil->VdmD[2](i,j);                     
	        }); 
	    });        
        Kokkos::deep_copy(derivUtil.VdmD_2,derivUtil.h_VdmD_2);
    }
}

template<int DIM>
void ProcessVarOpti::Load_elUtils(ElUtilGPU &elUtil)
{
	elUtil.idealMap = Kokkos::View<double**[10]> ("idealMap", elUtil.nElmt, elUtil.ptsHigh);
	elUtil.h_idealMap = Kokkos::create_mirror_view(elUtil.idealMap);

	elUtil.minJac = Kokkos::View<double*> ("minJac", elUtil.nElmt);
	elUtil.h_minJac = Kokkos::create_mirror_view(elUtil.minJac);
	elUtil.scaledJac = Kokkos::View<double*> ("scaledJac", elUtil.nElmt);
	elUtil.h_scaledJac = Kokkos::create_mirror_view(elUtil.scaledJac);

	int N1 = elUtil.nElmt;
    int M1 = elUtil.ptsHigh;
    Kokkos::parallel_for(team_policy_host(N1,M1), KOKKOS_LAMBDA (const member_type_host& teamMember)
    {         
        const int el = teamMember.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, M1 ), [&] (const int node)
        {                
            for (int i = 0; i < 10; i++)
        	{
        		elUtil.h_idealMap(el,node,i) = dataSet[el]->maps[node][i];
        	}
        	 
        });
        elUtil.h_minJac(el) = DBL_MAX;
        elUtil.h_scaledJac(el) = DBL_MAX;   

        //const int ElmtId = dataSet[el]->GetId();
        //elUtil.h_ElmtOffset(ElmtId) = el;
    });
    
    Kokkos::deep_copy(elUtil.idealMap,elUtil.h_idealMap);

}


void ProcessVarOpti::Create_nodes_view(NodesGPU &nodes)
{
    int nodes_size = nodes.nodes_size;
    int nElmt = nodes.nElmt;    
    nodes.X = Kokkos::View<double**> ("X",nElmt, nodes_size);
    nodes.h_X = Kokkos::create_mirror_view(nodes.X);
    nodes.Y = Kokkos::View<double**> ("Y",nElmt, nodes_size);
    nodes.h_Y = Kokkos::create_mirror_view(nodes.Y);
    nodes.Z = Kokkos::View<double**> ("Z",nElmt, nodes_size);
    nodes.h_Z = Kokkos::create_mirror_view(nodes.Z);
    nodes.Id = Kokkos::View<int**> ("Id",nElmt, nodes_size);
    nodes.h_Id = Kokkos::create_mirror_view(nodes.Id);

    nodes.ElmtOffset = Kokkos::View<int*> ("ElmtOffset",nElmt);
    nodes.h_ElmtOffset = Kokkos::create_mirror_view(nodes.ElmtOffset);
}


void ProcessVarOpti::Load_nodes(NodesGPU &nodes)
{
    int N1 = nodes.nElmt;
    int M1 = nodes.nodes_size;
    Kokkos::parallel_for(team_policy_host(N1,M1), KOKKOS_LAMBDA (const member_type_host& teamMember)
    {         
        const int el = teamMember.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, M1 ), [&] (const int node)
        {                
            nodes.h_X(el,node) = *dataSet[el]->nodes[node][0];
            nodes.h_Y(el,node) = *dataSet[el]->nodes[node][1];
            nodes.h_Z(el,node) = *dataSet[el]->nodes[node][2];
            nodes.h_Id(el,node) = *dataSet[el]->nodeIds[node];
        });
        const int ElmtId = dataSet[el]->GetId();
        nodes.h_ElmtOffset(ElmtId) = el;
    });
    Kokkos::deep_copy(nodes.X,nodes.h_X);
    Kokkos::deep_copy(nodes.Y,nodes.h_Y);
    Kokkos::deep_copy(nodes.Z,nodes.h_Z);
    Kokkos::deep_copy(nodes.Id,nodes.h_Id);

    Kokkos::deep_copy(nodes.ElmtOffset,nodes.h_ElmtOffset);
}

void ProcessVarOpti::Create_NodeMap(NodesGPU &nodes, 
			std::vector<std::vector<NodeSharedPtr> > &freenodes, NodeMap &nodeMap)
{
	int N1 = nodes.nElmt;
    int M1 = nodes.nodes_size;
    int id;
    //std::multimap<int,std::pair<int,int>> nodeMap;
    for(int ii = 0; ii < freenodes.size(); ii++)
    {
        for(int jj = 0; jj < freenodes[ii].size(); jj++)
        {
        	id = freenodes[ii][jj]->m_id;
		    for(int k = 0; k < N1; k++)
		    {         
		        for(int j = 0; j < M1; j++)
		        {                
		            if(nodes.h_Id(k,j) == id)
		            {
		                nodeMap.insert(std::pair<int,std::pair<int,int>> (id, {k,j} ));
		            }
		        }
		    }
		}
	}
}

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
        int nElmt = nodes.nElmt;   

        // declare min and max Jacobian
        Kokkos::View<double*,Kokkos::DefaultHostExecutionSpace> h_mx("mx",nElmt);
        Kokkos::View<double*,Kokkos::DefaultHostExecutionSpace> h_mn("mn",nElmt);   
        Kokkos::parallel_for(range_policy_host(0,nElmt), KOKKOS_LAMBDA (const int k)
        {
            h_mx(k) = -1.0 * DBL_MAX;
            h_mn(k) = DBL_MAX;
        });    

        // do the matrix vector multiplication on the GPU
        Kokkos::View<double[3][3]> dxdz("dxdz");
        Kokkos::View<double**> jacDet("jacDet", nElmt, nodes_size);
        typename Kokkos::View< double**>::HostMirror h_jacDet = Kokkos::create_mirror_view(jacDet);

        Kokkos::parallel_for( team_policy( nElmt , Kokkos::AUTO ), KOKKOS_LAMBDA ( const member_type& teamMember)
        {
            const int k = teamMember.league_rank();
        
            int N = nodes_size;
            int M = nodes_size;     
            /*Kokkos::parallel_for( team_policy( N , Kokkos::AUTO ), KOKKOS_LAMBDA ( const member_type& teamMember)
            {
                const int i = teamMember.league_rank();
                double x1i_i = 0.0;    
                Kokkos::parallel_reduce( Kokkos::TeamThreadRange( teamMember, M ), [&] (const int& j, double &update )
                {
                    update += derivUtil.VdmDL_0( i , j ) * X( k , j );              
                } ,x1i_i);
                Kokkos::single( Kokkos::PerTeam(teamMember), [&] ()
                {
                  x1i = x1i_i;
                });
            });*/
            Kokkos::parallel_for( Kokkos::TeamThreadRange( teamMember , N ), [&] ( const int i)
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
                            
                for (int j = 0; j < M; ++j)
                {
                    x1i += derivUtil.VdmDL_0( i , j ) * nodes.X( k , j );
                    y1i += derivUtil.VdmDL_0( i , j ) * nodes.Y( k , j );
                    z1i += derivUtil.VdmDL_0( i , j ) * nodes.Z( k , j );
                    x2i += derivUtil.VdmDL_1( i , j ) * nodes.X( k , j );
                    y2i += derivUtil.VdmDL_1( i , j ) * nodes.Y( k , j );
                    z2i += derivUtil.VdmDL_1( i , j ) * nodes.Z( k , j );
                    x3i += derivUtil.VdmDL_2( i , j ) * nodes.X( k , j );
                    y3i += derivUtil.VdmDL_2( i , j ) * nodes.Y( k , j );
                    z3i += derivUtil.VdmDL_2( i , j ) * nodes.Z( k , j );
                }
                  
                /*Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, M), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_0( i , j ) * X( k , j );
                }, x1i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, M), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_0( i , j ) * Y( k , j );
                }, y1i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, M), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_0( i , j ) * Z( k , j );
                }, z1i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, M), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_1( i , j ) * X( k , j );
                }, x2i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, M), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_1( i , j ) * Y( k , j );
                }, y2i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, M), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_1( i , j ) * Z( k , j );
                }, z2i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, M), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_2( i , j ) * X( k , j );
                }, x3i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, M), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_2( i , j ) * Y( k , j );
                }, y3i);
                Kokkos::parallel_reduce( Kokkos::ThreadVectorRange(teamMember, M), [&] (const int j, double &update )
                {
                    update += derivUtil.VdmDL_2( i , j ) * Z( k , j );
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

                jacDet(k,i) = dxdz(0,0)*(dxdz(1,1)*dxdz(2,2)-dxdz(2,1)*dxdz(1,2))
                                    -dxdz(0,1)*(dxdz(1,0)*dxdz(2,2)-dxdz(2,0)*dxdz(1,2))
                                    +dxdz(0,2)*(dxdz(1,0)*dxdz(2,1)-dxdz(2,0)*dxdz(1,1));             
            });

        });

        Kokkos::deep_copy(h_jacDet,jacDet);

        
        Kokkos::parallel_for(range_policy_host(0,nElmt), KOKKOS_LAMBDA (const int k)
        {
            for(int j = 0; j < nodes_size; j++)
            {
                //  mx = max(mx,jacDet);
                h_mx(k) = (h_mx(k) < h_jacDet(k,j) ? h_jacDet(k,j) : h_mx(k));
                //  mn = min(mn,jacDet);
                h_mn(k) = (h_mn(k) > h_jacDet(k,j) ? h_jacDet(k,j) : h_mn(k));
            }
            /*
            auto jacDet_sub = Kokkos::subview(jacDet, k, Kokkos::ALL);
            MaxFunctor <double> mxfunctor(jacDet_sub);
            Kokkos::parallel_reduce(range_policy(0, nodes_size) , mxfunctor, h_mx(k));

            MinFunctor <double> mnfunctor(jacDet_sub);
            Kokkos::parallel_reduce(range_policy(0, nodes_size) , mnfunctor, h_mn(k));      
            */   

            /*mtx3.lock();
            if(h_mn(k) < 0)
            {
                res.startInv++;
            }
            //  res.worstJac = min(res.worstJac,mn/mx);
            res.worstJac = (res.worstJac > h_mn(k)/h_mx(k) ? h_mn(k)/h_mx(k) : res.worstJac);
            mtx3.unlock();*/

            
            //dataSet[k]->minJac = h_mn(k);
            elUtil.h_minJac(k) = h_mn(k);     

            //dataSet[k]->scaledJac = h_mn(k)/h_mx(k);
            elUtil.h_scaledJac(k) = h_mn(k)/h_mx(k); 
        });

        Kokkos::deep_copy(elUtil.minJac,elUtil.h_minJac);
        Kokkos::deep_copy(elUtil.scaledJac,elUtil.h_scaledJac);
    }
}


void ProcessVarOpti::Optimise()
{

}


} // Utilities

} // Nektar

# endif