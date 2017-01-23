////////////////////////////////////////////////////////////////////////////////
//
//  File: LoadData.cpp
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

#ifndef UTILITIES_NEKMESH_PROCESSVAROPTI_LOADDATA
#define UTILITIES_NEKMESH_PROCESSVAROPTI_LOADDATA


#include "ProcessVarOpti.h"


namespace Nektar
{
namespace Utilities
{

void ProcessVarOpti::Load_derivUtil(DerivUtilGPU &derivUtil)
{
    int m_dim = dataSet[0]->m_dim;
    derivUtil.nodes_size = dataSet[0]->nodes.size();
    derivUtil.ptsLow = dataSet[0]->derivUtil->ptsLow;
    derivUtil.ptsHigh = dataSet[0]->derivUtil->ptsHigh;
    int nodes_size = derivUtil.nodes_size;
    int ptsHigh = derivUtil.ptsHigh;
    int ptsLow = derivUtil.ptsLow;

    derivUtil.VdmDL_0 = Kokkos::View<double**,random_memory>("VdmDL_0",ptsLow,nodes_size);
    derivUtil.h_VdmDL_0 = Kokkos::create_mirror_view(derivUtil.VdmDL_0);
    derivUtil.VdmDL_1 = Kokkos::View<double**,random_memory>("VdmDL_1",ptsLow,nodes_size);
    derivUtil.h_VdmDL_1 = Kokkos::create_mirror_view(derivUtil.VdmDL_1);
    Kokkos::parallel_for(team_policy_host(ptsLow,nodes_size), KOKKOS_LAMBDA (const member_type_host& teamMember)
    {         
        const int i = teamMember.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, nodes_size ), [&] (const int j)
        {                
            derivUtil.h_VdmDL_0(i,j) = dataSet[0]->derivUtil->VdmDL[0](i,j);
            derivUtil.h_VdmDL_1(i,j) = dataSet[0]->derivUtil->VdmDL[1](i,j);
            //printf("VdmDL_0(%i,%i) = %e\n", i,j, derivUtil.h_VdmDL_0(i,j));
        }); 
    });
    Kokkos::deep_copy(derivUtil.VdmDL_0,derivUtil.h_VdmDL_0);
    Kokkos::deep_copy(derivUtil.VdmDL_1,derivUtil.h_VdmDL_1);
    

    derivUtil.VdmD_0 = Kokkos::View<double**,random_memory>("VdmD_0",ptsHigh,nodes_size);
    derivUtil.h_VdmD_0 = Kokkos::create_mirror_view(derivUtil.VdmD_0);
    derivUtil.VdmD_1 = Kokkos::View<double**,random_memory>("VdmD_1",ptsHigh,nodes_size);
    derivUtil.h_VdmD_1 = Kokkos::create_mirror_view(derivUtil.VdmD_1);   
    Kokkos::parallel_for(team_policy_host(ptsHigh,nodes_size), KOKKOS_LAMBDA (const member_type_host& teamMember)
    {         
        const int i = teamMember.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, nodes_size ), [&] (const int j)
        {                
            derivUtil.h_VdmD_0(i,j) = dataSet[0]->derivUtil->VdmD[0](i,j);
            derivUtil.h_VdmD_1(i,j) = dataSet[0]->derivUtil->VdmD[1](i,j); 
            //printf("VdmD_0(%i,%i) = %e\n", i,j, derivUtil.h_VdmD_0(i,j));                     
        }); 
    });   
    Kokkos::deep_copy(derivUtil.VdmD_0,derivUtil.h_VdmD_0);
    Kokkos::deep_copy(derivUtil.VdmD_1,derivUtil.h_VdmD_1);


    // Quadrature Weights
    derivUtil.quadW = Kokkos::View<double*,random_memory>("quadW",ptsHigh);
    derivUtil.h_quadW = Kokkos::create_mirror_view(derivUtil.quadW);
    Kokkos::parallel_for(range_policy_host(0,ptsHigh), KOKKOS_LAMBDA (const int pts)
    {
    	derivUtil.h_quadW(pts) = dataSet[0]->derivUtil->quadW[pts];
    });
    Kokkos::deep_copy(derivUtil.quadW,derivUtil.h_quadW);

    
    if(m_dim == 3)
    {
        derivUtil.VdmDL_2 = Kokkos::View<double**,random_memory>("VdmDL_2",ptsLow,nodes_size);
        derivUtil.h_VdmDL_2 = Kokkos::create_mirror_view(derivUtil.VdmDL_2);
        Kokkos::parallel_for(team_policy_host(ptsLow,nodes_size), KOKKOS_LAMBDA (const member_type_host& teamMember)
	    {         
	        const int i = teamMember.league_rank();
	        Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, nodes_size ), [&] (const int j)
	        {                
	            derivUtil.h_VdmDL_2(i,j) = dataSet[0]->derivUtil->VdmDL[2](i,j);
	        }); 
	    });
	    Kokkos::deep_copy(derivUtil.VdmDL_2,derivUtil.h_VdmDL_2);


        derivUtil.VdmD_2 = Kokkos::View<double**,random_memory>("VdmD_2",ptsHigh,nodes_size);
        derivUtil.h_VdmD_2 = Kokkos::create_mirror_view(derivUtil.VdmD_2);
        Kokkos::parallel_for(team_policy_host(ptsHigh,nodes_size), KOKKOS_LAMBDA (const member_type_host& teamMember)
	    {         
	        const int i = teamMember.league_rank();
	        Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, nodes_size ), [&] (const int j)
	        {                
	            derivUtil.h_VdmD_2(i,j) = dataSet[0]->derivUtil->VdmD[2](i,j);                     
	        }); 
	    });        
        Kokkos::deep_copy(derivUtil.VdmD_2,derivUtil.h_VdmD_2);
    }
}


void ProcessVarOpti::Load_elUtils(ElUtilGPU &elUtil)
{
	elUtil.ptsHigh = dataSet[0]->derivUtil->ptsHigh;
    elUtil.globalnElmt = dataSet.size();
    int ptsHigh = elUtil.ptsHigh;
    int nElmt = elUtil.globalnElmt;

    elUtil.idealMap = Kokkos::View<double**[10],random_memory> ("idealMap", nElmt, ptsHigh);
	elUtil.h_idealMap = Kokkos::create_mirror_view(elUtil.idealMap);

	elUtil.minJac = Kokkos::View<double*> ("minJac", nElmt);
	elUtil.h_minJac = Kokkos::create_mirror_view(elUtil.minJac);
	elUtil.scaledJac = Kokkos::View<double*> ("scaledJac", nElmt);
	elUtil.h_scaledJac = Kokkos::create_mirror_view(elUtil.scaledJac);

	Kokkos::parallel_for(team_policy_host(nElmt,ptsHigh), KOKKOS_LAMBDA (const member_type_host& teamMember)
    {         
        const int el = teamMember.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, ptsHigh ), [&] (const int node)
        {                
            for (int i = 0; i < 10; i++)
        	{
        		elUtil.h_idealMap(el,node,i) = dataSet[el]->maps[node][i];
        	}
        	 
        });        
    });
    
    Kokkos::deep_copy(elUtil.idealMap,elUtil.h_idealMap);
}




void ProcessVarOpti::Load_nodes(NodesGPU &nodes)
{
    nodes.nodes_size = dataSet[0]->nodes.size();
    nodes.globalnElmt = dataSet.size();
    int nodes_size = nodes.nodes_size;
    int nElmt = nodes.globalnElmt;

    nodes.X = Kokkos::View<double**> ("X",nElmt, nodes_size);
    nodes.h_X = Kokkos::create_mirror_view(nodes.X);
    nodes.Y = Kokkos::View<double**> ("Y",nElmt, nodes_size);
    nodes.h_Y = Kokkos::create_mirror_view(nodes.Y);
    nodes.Z = Kokkos::View<double**> ("Z",nElmt, nodes_size);
    nodes.h_Z = Kokkos::create_mirror_view(nodes.Z);
    nodes.Id = Kokkos::View<int**> ("Id",nElmt, nodes_size);
    nodes.h_Id = Kokkos::create_mirror_view(nodes.Id);

    //nodes.ElmtOffset = Kokkos::View<int*> ("ElmtOffset",nElmt);
    //nodes.h_ElmtOffset = Kokkos::create_mirror_view(nodes.ElmtOffset);

    Kokkos::parallel_for(team_policy_host(nElmt,nodes_size), KOKKOS_LAMBDA (const member_type_host& teamMember)
    {         
        const int el = teamMember.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, nodes_size ), [&] (const int node)
        {                
            nodes.h_X(el,node) = *dataSet[el]->nodes[node][0];
            nodes.h_Y(el,node) = *dataSet[el]->nodes[node][1];
            nodes.h_Z(el,node) = *dataSet[el]->nodes[node][2];
            nodes.h_Id(el,node) = *dataSet[el]->nodeIds[node];
        });
        //const int ElmtId = dataSet[el]->GetId();
        //nodes.h_ElmtOffset(ElmtId) = el;
    });
    Kokkos::deep_copy(nodes.X,nodes.h_X);
    Kokkos::deep_copy(nodes.Y,nodes.h_Y);
    Kokkos::deep_copy(nodes.Z,nodes.h_Z);
    Kokkos::deep_copy(nodes.Id,nodes.h_Id);

    //Kokkos::deep_copy(nodes.ElmtOffset,nodes.h_ElmtOffset);
}

void ProcessVarOpti::Load_residual(Residual &res)
{
    res.startInv = Kokkos::View<int[1]>("startInv");
    res.h_startInv = Kokkos::create_mirror_view(res.startInv);
    
    res.worstJac = Kokkos::View<double[1]>("worstJac");
    res.h_worstJac = Kokkos::create_mirror_view(res.worstJac);    
    
    res.func = Kokkos::View<double[1]>("func");
    res.h_func = Kokkos::create_mirror_view(res.func);

    res.nReset = Kokkos::View<int[1]>("nReset");

    res.val = Kokkos::View<double[1]>("val");
    res.h_val = Kokkos::create_mirror_view(res.val); 
}





} // Utilities

} // Nektar

# endif