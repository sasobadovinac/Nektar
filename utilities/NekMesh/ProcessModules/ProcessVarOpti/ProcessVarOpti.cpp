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
#include "EvaluateGPU.hxx"

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


    if(m_config["Kokkos"].beenSet)
    {
        int nThreads = m_config["numthreads"].as<int>();
        Kokkos::InitArguments args;
        args.num_threads = nThreads;
        Kokkos::DefaultHostExecutionSpace::initialize(args.num_threads);
        Kokkos::DefaultExecutionSpace::initialize();

        // Test functions
        std::cout << " Template using int" << std::endl;  
        process<int>();
        std::cout << " Template using double" << std::endl; 
        process<double>();

        //p2rocess2(); 
    }

    Residual res;
    res.val = Kokkos::View<double[1]>("val");
    res.h_val = Kokkos::create_mirror_view(res.val);
    res.h_val[0] = 1.0;
    
    m_mesh->MakeOrder(m_mesh->m_nummode-1,LibUtilities::eGaussLobattoLegendre);
    BuildDerivUtil();
    GetElementMap();

    vector<ElementSharedPtr> elLock;

    if(m_config["region"].beenSet)
    {
        elLock = GetLockedElements(m_config["region"].as<NekDouble>());
    }

    vector<vector<NodeSharedPtr> > freenodes = GetColouredNodes(elLock, res);
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
                    optiType, freenodes[i][j], it->second, derivUtil, opti));
        }
        optiNodes.push_back(ns);
    }

    // determine statistics of element colouring
    int nset = optiNodes.size();
    int p = 0;
    int mn = INT_MAX;
    int mx = 0;
    for(int i = 0; i < nset; i++)
    {
        p += optiNodes[i].size();
        mn = min(mn, int(optiNodes[i].size()));
        mx = max(mx, int(optiNodes[i].size()));
    }

    // initialise variable for number of invalid elements
    res.startInv = Kokkos::View<int[1]>("startInv");
    res.h_startInv = Kokkos::create_mirror_view(res.startInv);
    res.h_startInv[0] = 0;
    Kokkos::deep_copy(res.startInv,res.h_startInv);

    // initialise variable for smallest Jacobian of all elements
    res.worstJac = Kokkos::View<double[1]>("worstJac");
    res.h_worstJac = Kokkos::create_mirror_view(res.worstJac);    
    res.h_worstJac[0] = DBL_MAX;

    res.func = Kokkos::View<double[1]>("func");
    res.h_func = Kokkos::create_mirror_view(res.func);
    res.nReset = Kokkos::View<int[1]>("nReset");
    res.h_nReset = Kokkos::create_mirror_view(res.nReset);
    
    
    // initialise and load derivUtils onto GPU
    DerivUtilGPU derivUtil;
    derivUtil.nodes_size = dataSet[0]->nodes.size();
    derivUtil.ptsLow = dataSet[0]->derivUtil->ptsLow;
    derivUtil.ptsHigh = dataSet[0]->derivUtil->ptsHigh;
    Load_derivUtil(derivUtil);//
    

    // load mapping ideal to ref element
    ElUtilGPU elUtil;
    elUtil.ptsHigh = dataSet[0]->derivUtil->ptsHigh;
    elUtil.globalnElmt = dataSet.size();
    Load_elUtils<3>(elUtil);
    

    // initialise views for node coordinates on GPU
    NodesGPU nodes;
    nodes.nodes_size = dataSet[0]->nodes.size();
    nodes.globalnElmt = dataSet.size();
    Create_nodes_view(nodes); 

    // copy nodes onto GPU and evaluate the elements on GPU
    Load_nodes(nodes);

    // create mapping between node ids and their position on node array
    NodeMap nodeMap;
    Create_NodeMap(nodes, freenodes, nodeMap);

    // Evaluate the Jacobian of all elements on GPU
    Evaluate(derivUtil, nodes, elUtil, res);

    
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
         << "N elements invalid:\t" << res.h_startInv[0] << endl
         << "Worst jacobian:\t\t" << res.h_worstJac[0] << endl
         << "N free nodes:\t\t" << res.n << endl
         << "N Dof:\t\t\t" << res.nDoF << endl
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

    while (ctr < maxIter && res.h_val[0] > restol)
    {
        ctr++;
        res.h_val[0] = 0.0;
        res.h_func[0] = 0.0;
        res.h_nReset[0] = 0;
        Kokkos::deep_copy(res.val,res.h_val);
        Kokkos::deep_copy(res.func,res.h_func);
        Kokkos::deep_copy(res.nReset,res.h_nReset);
        
        for(int cs = 0; cs < optiNodes.size(); cs++)
        {
            nodes.coloursetSize = optiNodes[cs].size();

            //construct Vector that contains the number of elements that each node in the colorset is associated with
            nodes.nElmtArray = Kokkos::View<int*> ("nElmtArray",nodes.coloursetSize);
            nodes.h_nElmtArray = Kokkos::create_mirror_view(nodes.nElmtArray);
            int maxnElmt = 0;
            for(int node = 0; node < nodes.coloursetSize; node++)
            {
                //const int nElmt = optiNodes[cs][node]->data.size();
                nodes.h_nElmtArray(node) = optiNodes[cs][node]->data.size();
                maxnElmt = (maxnElmt < nodes.h_nElmtArray(node) ? nodes.h_nElmtArray(node) : maxnElmt);
            }

            // construct Array that contains the element Ids of all elements that each node in the colorset is associated with
            nodes.elIdArray = Kokkos::View<int**> ("elIdArray",nodes.coloursetSize, maxnElmt);
            nodes.h_elIdArray = Kokkos::create_mirror_view(nodes.elIdArray);

            // construct Array that contains the local node Ids in all elements that each node in the colorset is associated with
            nodes.localNodeIdArray = Kokkos::View<int**> ("localNodeIdArray",nodes.coloursetSize, maxnElmt);
            nodes.h_localNodeIdArray = Kokkos::create_mirror_view(nodes.localNodeIdArray);

            for(int node = 0; node < nodes.coloursetSize; node++)            
            {
                const int globalNodeId = optiNodes[cs][node]->node->m_id;
                NodeMap::const_iterator coeffs;
                coeffs = nodeMap.find(globalNodeId); 
                for (int el = 0; el < nodes.h_nElmtArray(node); ++el)
                {
                    nodes.h_elIdArray(node,el) = optiNodes[cs][node]->data[el]->GetId();

                    int elmt = std::get<0>(coeffs->second);
                    if (elmt == nodes.h_elIdArray(node,el))
                    {
                        nodes.h_localNodeIdArray(node,el) = std::get<1>(coeffs->second);
                    }            
                    coeffs++;
                }
            }
            Kokkos::deep_copy(nodes.localNodeIdArray,nodes.h_localNodeIdArray);
            Kokkos::deep_copy(nodes.elIdArray,nodes.h_elIdArray);
            Kokkos::deep_copy(nodes.nElmtArray,nodes.h_nElmtArray);


            
            for(int node = 0; node < nodes.coloursetSize; node++)
            {   
                //optiNodes[cs][node]->Optimise(derivUtil, nodes, nodeMap, elUtil, res, nElmt);
                //printf("node %cs finished\n", node);
            }
            //printf("colorset %cs finished\n", cs);
            
            OptimiseGPU(derivUtil, nodes, nodeMap, elUtil, res);
            
        }

        res.h_startInv[0] = 0;
        Kokkos::deep_copy(res.startInv,res.h_startInv);
        res.h_worstJac[0] = DBL_MAX;

        // Evaluate the elements on GPU
        Evaluate(derivUtil, nodes, elUtil, res);
        
        Kokkos::deep_copy(res.h_val,res.val);
        Kokkos::deep_copy(res.h_func,res.func);
        Kokkos::deep_copy(res.h_nReset,res.nReset);

        if(m_config["resfile"].beenSet)
        {
            resFile << res.h_val[0] << " " << res.h_worstJac[0] << " " << res.h_func[0] << endl;
        }

        cout << ctr << "\tResidual: " << res.h_val[0]
                    << "\tMin Jac: " << res.h_worstJac[0]
                    << "\tInvalid: " << res.h_startInv[0]
                    << "\tReset nodes: " << res.h_nReset[0]
                    << "\tFunctional: " << res.h_func[0]
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

    cout << "Invalid at end:\t\t" << res.h_startInv[0] << endl;
    cout << "Worst at end:\t\t" << res.h_worstJac[0] << endl;

    if(m_config["Kokkos"].beenSet)
    {
        Kokkos::DefaultExecutionSpace::finalize();
        Kokkos::DefaultHostExecutionSpace::finalize();        
    }
}

} // Utilities
} // Nektar

