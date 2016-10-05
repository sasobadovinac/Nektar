////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.cpp
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
//  Description: Calculate Jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/ManagerAccess.h>

#include <NekMeshUtils/MeshElements/Element.h>
#include "ProcessVarOpti.h"
#include "NodeOpti.h"
#include "ElUtil.h"

#include <boost/thread/mutex.hpp>

#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdPrismExp.h>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Foundations/NodalUtil.h>

#include <Kokkos_Core.hpp>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>


using namespace std;
using namespace Nektar::NekMeshUtils;


void p2rocess2()
{

  int num_vectors = 1000; // number of vectors
  int len  = 4000;       // length of vectors 
  int nrepeat = 10;       // number of repeats of the test

    // allocate space for vectors to do num_vectors dot products of length len
  Kokkos::View<double**,Kokkos::LayoutRight> a("A",num_vectors,len);
  Kokkos::View<double**,Kokkos::LayoutRight> b("B",num_vectors,len);
  Kokkos::View<double*>  c("C",num_vectors);
  auto h_c = Kokkos::create_mirror_view(c);
  
  //int team_size = std::is_same<Kokkos::DefaultExecutionSpace,Kokkos::DefaultHostExecutionSpace>::value?1:256;

  typedef Kokkos::TeamPolicy<>::member_type team_member;

  // Initialize vectors
  Kokkos::parallel_for( Kokkos::TeamPolicy<>(num_vectors,Kokkos::AUTO), 
                        KOKKOS_LAMBDA (const team_member& thread) {
    const int i = thread.league_rank();
    Kokkos::parallel_for( Kokkos::TeamThreadRange(thread,len), [&] (const int& j) {
      a(i,j) = i+1;
      b(i,j) = j+1;
    });
    //Kokkos::single( Kokkos::PerTeam(thread), [&] () {
    //  c(i) = 0.0;
    //});
  });


  // Time dot products
  struct timeval begin,end;

  gettimeofday(&begin,NULL);

  for(int repeat = 0; repeat < nrepeat; repeat++)
  {
    Kokkos::parallel_for( Kokkos::TeamPolicy<>(num_vectors,Kokkos::AUTO), 
                          KOKKOS_LAMBDA (const team_member& thread) {
      double c_i = 0.0;
      const int i = thread.league_rank();
      Kokkos::parallel_reduce( Kokkos::TeamThreadRange(thread,len), [&] (const int& j, double& ctmp) {
        ctmp += a(i,j) * b(i,j);
      },c_i);
      Kokkos::single( Kokkos::PerTeam(thread), [&] () {
        c(i) = c_i;
      });
    });
    /*Kokkos::parallel_for( num_vectors, KOKKOS_LAMBDA (const int i)
    {
      double c_i = 0.0;
      for( int j = 0; j < len; j++)
      {
        c_i += a(i,j) * b(i,j);
      }
      c(i) = c_i;      
    });*/
  }

  gettimeofday(&end,NULL);

  // Calculate time
  double time = 1.0*(end.tv_sec-begin.tv_sec) + 1.0e-6*(end.tv_usec-begin.tv_usec);

  // Error check
  Kokkos::deep_copy(h_c,c);
  int error = 0;
  for(int i = 0; i < num_vectors; i++)
  {
    double diff = ((h_c(i) - 1.0*(i+1)*len*(len+1)/2))/((i+1)*len*(len+1)/2);
    if ( diff*diff>1e-20 )
    { 
      error = 1;
      printf("Error: %i %i %i %lf %lf %e %lf\n",i,num_vectors,len,h_c(i),1.0*(i+1)*len*(len+1)/2,h_c(i) - 1.0*(i+1)*len*(len+1)/2,diff);
    }
  }

  // Print results (problem size, time and bandwidth in GB/s)

  if(error==0) { 
    printf("#NumVector Length Time(s) ProblemSize(MB) Bandwidth(GB/s)\n");
    printf("%i %i %e %lf %lf\n",num_vectors,len,time,1.0e-6*num_vectors*len*2*8,1.0e-9*num_vectors*len*2*8*nrepeat/time);
  }
  else printf("Error\n");
  
}




namespace Nektar
{
namespace Utilities
{

boost::mutex mtx3;

ModuleKey ProcessVarOpti::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "varopti"),
    ProcessVarOpti::create,
    "Optimise mesh locations.");

ProcessVarOpti::ProcessVarOpti(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["linearelastic"] =
        ConfigOption(true, "", "Optimise for linear elasticity");
    m_config["winslow"] =
        ConfigOption(true, "", "Optimise for winslow");
    m_config["roca"] =
        ConfigOption(true, "", "Optimise for roca method");
    m_config["hyperelastic"] =
        ConfigOption(true, "", "Optimise for hyper elasticity");
    m_config["numthreads"] =
        ConfigOption(false, "1", "Number of threads");
    m_config["restol"] =
        ConfigOption(false, "1e-6", "Tolerance criterion");
    m_config["maxiter"] =
        ConfigOption(false, "500", "Maximum number of iterations");
    m_config["nq"] =
        ConfigOption(false, "-1", "Order of mesh");
    m_config["Boost"] =
        ConfigOption(true, "", "Parallelise with Boost");
    m_config["Kokkos"] =
        ConfigOption(true, "", "Parallelise with Kokkos");
    m_config["region"] =
        ConfigOption(false, "0.0", "create regions based on target");
    m_config["resfile"] =
        ConfigOption(false, "", "writes residual values to file");
    m_config["histfile"] =
        ConfigOption(false, "", "histogram of scaled jac");
}

ProcessVarOpti::~ProcessVarOpti()
{
}

void ProcessVarOpti::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "ProcessVarOpti: Optimising... " << endl;
    }

    if(m_config["linearelastic"].beenSet)
    {
        opti = eLinEl;
    }
    else if(m_config["winslow"].beenSet)
    {
        opti = eWins;
    }
    else if(m_config["roca"].beenSet)
    {
        opti = eRoca;
    }
    else if(m_config["hyperelastic"].beenSet)
    {
        opti = eHypEl;
    }
    else
    {
        ASSERTL0(false,"no opti type set");
    }

    if(m_config["Kokkos"].beenSet)
    {
        ThreadManagerType = "ThreadManagerKokkos";        
    }
    else if(m_config["Boost"].beenSet)
    {
        ThreadManagerType = "ThreadManagerBoost";
    }
    else
    {
        ThreadManagerType = "ThreadManagerBoost";
        cout << endl << "Default parallelisation using Boost" << endl;
    }

    const int maxIter = m_config["maxiter"].as<int>();
    const NekDouble restol = m_config["restol"].as<NekDouble>();

    
    EdgeSet::iterator eit;
    m_mesh->m_nummode = -1;
    if (m_mesh->m_nummode == -1)
    {
        bool fd = false;
        for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
        {
            if((*eit)->m_edgeNodes.size() > 0)
            {
                m_mesh->m_nummode = (*eit)->m_edgeNodes.size() + 2;
                fd = true;
                break;
            }
        }

        if(!fd && m_config["nq"].beenSet)
        {
            m_mesh->m_nummode = m_config["nq"].as<int>();
            fd = true;
        }

        ASSERTL0(fd,"failed to find order of mesh");
    }

    if(m_mesh->m_verbose)
    {
        cout << "Identified mesh order as: " << m_mesh->m_nummode - 1 << endl;
    }

    if (m_mesh->m_expDim == 2 && m_mesh->m_spaceDim == 3)
    {
        ASSERTL0(false,"cannot deal with manifolds");
    }

    res = boost::shared_ptr<Residual>(new Residual);
    res->val = 1.0;

    m_mesh->MakeOrder(m_mesh->m_nummode-1,LibUtilities::eGaussLobattoLegendre);
    BuildDerivUtil();
    GetElementMap();

    vector<ElementSharedPtr> elLock;

    if(m_config["region"].beenSet)
    {
        elLock = GetLockedElements(m_config["region"].as<NekDouble>());
    }

    vector<vector<NodeSharedPtr> > freenodes = GetColouredNodes(elLock);
    vector<vector<NodeOptiSharedPtr> > optiNodes;

    //turn the free nodes into optimisable objects with all required data
    for(int i = 0; i < freenodes.size(); i++)
    {
        vector<NodeOptiSharedPtr> ns;
        for(int j = 0; j < freenodes[i].size(); j++)
        {
            NodeElMap::iterator it = nodeElMap.find(freenodes[i][j]->m_id);
            ASSERTL0(it != nodeElMap.end(), "could not find");

            int optiType = m_mesh->m_spaceDim;

            if (freenodes[i][j]->GetNumCadCurve() == 1)
            {
                optiType += 10;
            }
            else if (freenodes[i][j]->GetNumCADSurf() == 1)
            {
                optiType += 20;
            }
            else
            {
                optiType += 10*m_mesh->m_expDim;
            }

            ns.push_back(
                GetNodeOptiFactory().CreateInstance(
                    optiType, freenodes[i][j], it->second, res, derivUtil, opti));
        }
        optiNodes.push_back(ns);
    }

    int nset = optiNodes.size();
    int p = 0;
    int mn = numeric_limits<int>::max();
    int mx = 0;
    for(int i = 0; i < nset; i++)
    {
        p += optiNodes[i].size();
        mn = min(mn, int(optiNodes[i].size()));
        mx = max(mx, int(optiNodes[i].size()));
    }

    res->startInv =0;
    res->worstJac = numeric_limits<double>::max();

    if(m_config["Kokkos"].beenSet)
    {
        int nThreads = m_config["numthreads"].as<int>();
        Kokkos::InitArguments args;
        args.num_threads = nThreads;
        Kokkos::DefaultHostExecutionSpace::initialize(args.num_threads);
        Kokkos::DefaultExecutionSpace::initialize();

        std::cout << " Template using int" << std::endl;  
        process<int>();
        std::cout << " Template using double" << std::endl; 
        process<double>();

        //p2rocess2(); 
    }
    
    // initialise and load derivUtils onto GPU
    DerivUtilGPU derivUtil;
    derivUtil.nodes_size = dataSet[0]->nodes.size();
    Load_derivUtil(derivUtil);

    // initialise views for node coordinates on GPU
    NodesGPU nodes;
    nodes.nodes_size = dataSet[0]->nodes.size();
    nodes.dataSet_size = dataSet.size();
    Create_nodes_view(nodes);    

    // copy nodes onto GPU and evaluate the elements on GPU
    Load_nodes(nodes);
    Evaluate(derivUtil, nodes);

    // Evaluate elements on CPU
    /*Kokkos::parallel_for(range_policy_host(0,dataSet.size()), KOKKOS_LAMBDA (const int i)
    {
        dataSet[i]->Evaluate();
    });*/ 


    if(m_config["histfile"].beenSet)
    {
        ofstream histFile;
        string name = m_config["histfile"].as<string>() + "_start.txt";
        histFile.open(name.c_str());

        for(int i = 0; i < dataSet.size(); i++)
        {
            histFile << dataSet[i]->scaledJac << endl;
        }
        histFile.close();
    }

    cout << scientific << endl;
    cout << "N elements:\t\t" << m_mesh->m_element[m_mesh->m_expDim].size() - elLock.size() << endl
         << "N elements invalid:\t" << res->startInv << endl
         << "Worst jacobian:\t\t" << res->worstJac << endl
         << "N free nodes:\t\t" << res->n << endl
         << "N Dof:\t\t\t" << res->nDoF << endl
         << "N color sets:\t\t" << nset << endl
         << "Avg set colors:\t\t" << p/nset << endl
         << "Min set:\t\t" << mn << endl
         << "Max set:\t\t" << mx << endl
         << "Residual tolerance:\t" << restol << endl;   

    Timer t;
    t.Start();

    int ctr = 0;
    ofstream resFile;
    if(m_config["resfile"].beenSet)
    {
        resFile.open(m_config["resfile"].as<string>().c_str());
    }

    while (ctr < maxIter && res->val > restol)
    {
        ctr++;
        res->val = 0.0;
        res->func = 0.0;
        res->nReset = 0;
        for(int i = 0; i < optiNodes.size(); i++)
        {
            Kokkos::parallel_for(range_policy_host(0,optiNodes[i].size()),KOKKOS_LAMBDA (const int j)
            {
                optiNodes[i][j]->Optimise();
            });
        }

        res->startInv = 0;
        res->worstJac = numeric_limits<double>::max();

        // copy nodes onto GPU and evaluate the elements on GPU
        Load_nodes(nodes);
        Evaluate(derivUtil, nodes);

        // Evaluate elements on CPU
        /*Kokkos::parallel_for(range_policy_host(0,dataSet.size()), KOKKOS_LAMBDA (const int i)
        {
            dataSet[i]->Evaluate();
        });*/

        if(m_config["resfile"].beenSet)
        {
            resFile << res->val << " " << res->worstJac << " " << res->func << endl;
        }

        cout << ctr << "\tResidual: " << res->val
                    << "\tMin Jac: " << res->worstJac
                    << "\tInvalid: " << res->startInv
                    << "\tReset nodes: " << res->nReset
                    << "\tFunctional: " << res->func
                    << endl;
    }

    if(m_config["histfile"].beenSet)
    {
        ofstream histFile;
        string name = m_config["histfile"].as<string>() + "_end.txt";
        histFile.open(name.c_str());

        for(int i = 0; i < dataSet.size(); i++)
        {
            histFile << dataSet[i]->scaledJac << endl;
        }
        histFile.close();
    }
    if(m_config["resfile"].beenSet)
    {
        resFile.close();
    }

    t.Stop();
    cout << "Time to compute: " << t.TimePerTest(1) << endl;

    cout << "Invalid at end:\t\t" << res->startInv << endl;
    cout << "Worst at end:\t\t" << res->worstJac << endl;

    if(m_config["Kokkos"].beenSet)
    {
        Kokkos::DefaultExecutionSpace::finalize();
        Kokkos::DefaultHostExecutionSpace::finalize();        
    }
}


void ProcessVarOpti::Load_derivUtil(DerivUtilGPU &derivUtil)
{
    //typedef Kokkos::MemoryTraits<Kokkos::RandomAccess> random_memory; //?
    int nodes_size = derivUtil.nodes_size;
    derivUtil.VdmDL_0 = Kokkos::View<double**>("VdmDL_0",nodes_size,nodes_size);
    derivUtil.h_VdmDL_0 = Kokkos::create_mirror_view(derivUtil.VdmDL_0);
    derivUtil.VdmDL_1 = Kokkos::View<double**>("VdmDL_1",nodes_size,nodes_size);
    derivUtil.h_VdmDL_1 = Kokkos::create_mirror_view(derivUtil.VdmDL_1);
    
    NekMatrix<NekDouble> derivUtil_VdmDL_0 = dataSet[0]->derivUtil->VdmDL[0];
    NekMatrix<NekDouble> derivUtil_VdmDL_1 = dataSet[0]->derivUtil->VdmDL[1];
    int N2 = nodes_size;
    int M2 = nodes_size;
    Kokkos::parallel_for(team_policy_host(N2,M2), KOKKOS_LAMBDA (const member_type_host& teamMember)
    {         
        const int i = teamMember.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, M2 ), [&] (const int j)
        {                
            derivUtil.h_VdmDL_0(i,j) = derivUtil_VdmDL_0(i,j);
            derivUtil.h_VdmDL_1(i,j) = derivUtil_VdmDL_1(i,j);                      
        }); 
    });        
   
    Kokkos::deep_copy(derivUtil.VdmDL_0,derivUtil.h_VdmDL_0);
    Kokkos::deep_copy(derivUtil.VdmDL_1,derivUtil.h_VdmDL_1);


    int m_dim = dataSet[0]->m_dim;
    if(m_dim == 3)
    {
        derivUtil.VdmDL_2 = Kokkos::View<double**>("VdmDL_2",nodes_size,nodes_size);
        derivUtil.h_VdmDL_2 = Kokkos::create_mirror_view(derivUtil.VdmDL_2);
        
        NekMatrix<NekDouble> derivUtil_VdmDL_2 = dataSet[0]->derivUtil->VdmDL[2];
        int N2 = nodes_size;
        int M2 = nodes_size;
        Kokkos::parallel_for(team_policy_host(N2,M2), KOKKOS_LAMBDA (const member_type_host& teamMember)
        {         
            const int i = teamMember.league_rank();
            Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, M2 ), [&] (const int j)
            {                
                derivUtil.h_VdmDL_2(i,j) = derivUtil_VdmDL_2(i,j);                
            }); 
        }); 
        
        Kokkos::deep_copy(derivUtil.VdmDL_2,derivUtil.h_VdmDL_2);
    }
}


void ProcessVarOpti::Create_nodes_view(NodesGPU &nodes)
{
    int nodes_size = nodes.nodes_size;
    int dataSet_size = nodes.dataSet_size;    
    nodes.X = Kokkos::View<double**> ("X",dataSet_size, nodes_size);
    nodes.h_X = Kokkos::create_mirror_view(nodes.X);
    nodes.Y = Kokkos::View<double**> ("Y",dataSet_size, nodes_size);
    nodes.h_Y = Kokkos::create_mirror_view(nodes.Y);
    nodes.Z = Kokkos::View<double**> ("Z",dataSet_size, nodes_size);
    nodes.h_Z = Kokkos::create_mirror_view(nodes.Z);
}


void ProcessVarOpti::Load_nodes(NodesGPU &nodes)
{
    int N1 = nodes.dataSet_size;
    int M1 = nodes.nodes_size;
    Kokkos::parallel_for(team_policy_host(N1,M1), KOKKOS_LAMBDA (const member_type_host& teamMember)
    {         
        const int k = teamMember.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange( teamMember, M1 ), [&] (const int j)
        {                
            nodes.h_X(k,j) = *dataSet[k]->nodes[j][0];
            nodes.h_Y(k,j) = *dataSet[k]->nodes[j][1];
            nodes.h_Z(k,j) = *dataSet[k]->nodes[j][2];              
        }); 
    });
    Kokkos::deep_copy(nodes.X,nodes.h_X);
    Kokkos::deep_copy(nodes.Y,nodes.h_Y);
    Kokkos::deep_copy(nodes.Z,nodes.h_Z);
}

void ProcessVarOpti::Evaluate(DerivUtilGPU &derivUtil,NodesGPU &nodes)
{
    int m_dim = dataSet[0]->m_dim;
    if(m_dim == 2)
    {
        ASSERTL0(false,"not coded");
    }
    else if (m_dim == 3)
    {
        int nodes_size = nodes.nodes_size;
        int dataSet_size = nodes.dataSet_size;   

        // declare min and max Jacobian
        Kokkos::View<double*> mx("mx",dataSet_size);
        typename Kokkos::View< double*>::HostMirror h_mx = Kokkos::create_mirror_view(mx);
        Kokkos::View<double*> mn("mn",dataSet_size);
        typename Kokkos::View< double*>::HostMirror h_mn = Kokkos::create_mirror_view(mn);    
        Kokkos::parallel_for(range_policy_host(0,dataSet_size), KOKKOS_LAMBDA (const int k)
        {
          h_mx(k) = -1.0 * DBL_MAX;
          h_mn(k) = DBL_MAX;
        });    

        // do the matrix vector multiplication on the GPU
        Kokkos::View<double[3][3]> dxdz("dxdz");
        Kokkos::View<double**> jacDet("jacDet", dataSet_size, nodes_size);
        typename Kokkos::View< double**>::HostMirror h_jacDet = Kokkos::create_mirror_view(jacDet);

        Kokkos::parallel_for( team_policy( dataSet_size , Kokkos::AUTO ), KOKKOS_LAMBDA ( const member_type& teamMember)
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

        
        Kokkos::parallel_for(range_policy_host(0,dataSet_size), KOKKOS_LAMBDA (const int k)
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

            mtx3.lock();
            if(h_mn(k) < 0)
            {
                dataSet[k]->res->startInv++;
            }
            //  res->worstJac = min(res->worstJac,mn/mx);
            dataSet[k]->res->worstJac = (dataSet[k]->res->worstJac > h_mn(k)/h_mx(k) ? 
                                          h_mn(k)/h_mx(k) : dataSet[k]->res->worstJac);
            mtx3.unlock();

            dataSet[k]->maps = dataSet[k]->MappingIdealToRef();
            
            dataSet[k]->minJac = h_mn(k);
            dataSet[k]->scaledJac = h_mn(k)/h_mx(k); 
        });
    }

}

} // Utilities
} // Nektar

