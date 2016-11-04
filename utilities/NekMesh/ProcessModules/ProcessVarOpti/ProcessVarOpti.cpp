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


//#include <boost/thread/mutex.hpp>

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
        //Kokkos::DefaultHostExecutionSpace::initialize(args.num_threads);
        //Kokkos::DefaultExecutionSpace::initialize();
        Kokkos::initialize(args);

    }

    
     
    m_mesh->MakeOrder(m_mesh->m_nummode-1,LibUtilities::eGaussLobattoLegendre);
    BuildDerivUtil();
    GetElementMap();

    vector<ElementSharedPtr> elLock;

    if(m_config["region"].beenSet)
    {
        elLock = GetLockedElements(m_config["region"].as<NekDouble>());
    }

    Residual res;
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

    cout << scientific << endl;
    cout << "N elements:\t\t" << m_mesh->m_element[m_mesh->m_expDim].size() - elLock.size() << endl
         << "N free nodes:\t\t" << res.n << endl
         << "N Dof:\t\t\t" << res.nDoF << endl
         << "N color sets:\t\t" << nset << endl
         << "Avg set colors:\t\t" << p/nset << endl
         << "Min set:\t\t" << mn << endl
         << "Max set:\t\t" << mx << endl
         << "Residual tolerance:\t" << restol << endl;


    // initialise residual variables
    res.startInv = Kokkos::View<int[1]>("startInv");
    res.h_startInv = Kokkos::create_mirror_view(res.startInv);
    
    res.worstJac = Kokkos::View<double[1]>("worstJac");
    res.h_worstJac = Kokkos::create_mirror_view(res.worstJac);    
    
    res.func = Kokkos::View<double[1]>("func");
    res.h_func = Kokkos::create_mirror_view(res.func);

    res.nReset = Kokkos::View<int[1]>("nReset");

    res.val = Kokkos::View<double[1]>("val");
    res.h_val = Kokkos::create_mirror_view(res.val);   
    
    
    // initialise and load derivUtils onto GPU
    DerivUtilGPU derivUtil;
    derivUtil.nodes_size = dataSet[0]->nodes.size();
    derivUtil.ptsLow = dataSet[0]->derivUtil->ptsLow;
    derivUtil.ptsHigh = dataSet[0]->derivUtil->ptsHigh;
    Load_derivUtil(derivUtil);
    

    // load mapping ideal to ref element
    ElUtilGPU elUtil;
    elUtil.ptsHigh = dataSet[0]->derivUtil->ptsHigh;
    elUtil.globalnElmt = dataSet.size();
    Load_elUtils(elUtil);
    

    // initialise views for node coordinates on GPU and copy nodes onto GPU
    NodesGPU nodes;
    nodes.nodes_size = dataSet[0]->nodes.size();
    nodes.globalnElmt = dataSet.size();
    Load_nodes(nodes);
    

    // create mapping between node ids and their position on node array
    
    // construct vector that contains the number of nodes in each colourset
    int nColoursets = optiNodes.size();
    nodes.coloursetSize = Kokkos::View<int*> ("coloursetSize",nColoursets);
    nodes.h_coloursetSize = Kokkos::create_mirror_view(nodes.coloursetSize);
    int maxColoursetSize = 0;
    for(int cs = 0; cs < nColoursets; cs++)
    {            
        nodes.h_coloursetSize(cs) = optiNodes[cs].size();
        maxColoursetSize = (maxColoursetSize < nodes.h_coloursetSize(cs)
                ? nodes.h_coloursetSize(cs) : maxColoursetSize);
    }
    //construct Vector that contains the number of elements that each node in the colorset is associated with
    nodes.nElmtArray = Kokkos::View<int**> ("nElmtArray",nColoursets,maxColoursetSize);
    nodes.h_nElmtArray = Kokkos::create_mirror_view(nodes.nElmtArray);

    int maxnElmt = 0;
    for(int cs = 0; cs < nColoursets; cs++)
    {
        for(int node = 0; node < nodes.h_coloursetSize(cs); node++)
        {
            nodes.h_nElmtArray(cs,node) = optiNodes[cs][node]->data.size();
            maxnElmt = (maxnElmt < nodes.h_nElmtArray(cs,node)
                    ? nodes.h_nElmtArray(cs,node) : maxnElmt);
        }
    }
    // construct Array that contains the element Ids of all elements that each node in the colorset is associated with
    nodes.elIdArray = Kokkos::View<int***> ("elIdArray",nColoursets,maxColoursetSize, maxnElmt);
    nodes.h_elIdArray = Kokkos::create_mirror_view(nodes.elIdArray);

    // construct Array that contains the local node Ids in all elements that each node in the colorset is associated with
    nodes.localNodeIdArray = Kokkos::View<int***> ("localNodeIdArray",nColoursets,maxColoursetSize, maxnElmt);
    nodes.h_localNodeIdArray = Kokkos::create_mirror_view(nodes.localNodeIdArray);

    for(int cs = 0; cs < nColoursets; cs++)
    {
        for(int node = 0; node < nodes.h_coloursetSize(cs); node++)            
        {
            const int globalNodeId = optiNodes[cs][node]->node->m_id;

            for (int el = 0; el < nodes.h_nElmtArray(cs,node); ++el)
            {                
                int globalEl = optiNodes[cs][node]->data[el]->GetId();
                nodes.h_elIdArray(cs,node,el) = globalEl;

                for(int ln = 0; ln < nodes.nodes_size; ln++)
                {
                    if(nodes.h_Id(globalEl,ln) == globalNodeId)
                    {
                        nodes.h_localNodeIdArray(cs,node,el) = ln;
                    }
                }
            }
        }
    }

    Kokkos::deep_copy(nodes.localNodeIdArray,nodes.h_localNodeIdArray);
    Kokkos::deep_copy(nodes.elIdArray,nodes.h_elIdArray);
    Kokkos::deep_copy(nodes.nElmtArray,nodes.h_nElmtArray);
    Kokkos::deep_copy(nodes.coloursetSize,nodes.h_coloursetSize);
    // END INDEXING ----------------------------------------------


    // Evaluate the Jacobian of all elements on GPU
    Evaluate(derivUtil, nodes, elUtil, res);

    
    if(m_config["histfile"].beenSet)
    {
        ofstream histFile;
        string name = m_config["histfile"].as<string>() + "_start.txt";
        histFile.open(name.c_str());
        Kokkos::deep_copy(elUtil.h_scaledJac,elUtil.scaledJac);

        for(int i = 0; i < dataSet.size(); i++)
        {
            //histFile << dataSet[i]->scaledJac << endl;
            histFile << elUtil.h_scaledJac(i) << endl;
        }
        histFile.close();
    }  


    Timer t;
    t.Start();

    
    ofstream resFile;
    if(m_config["resfile"].beenSet)
    {
        resFile.open(m_config["resfile"].as<string>().c_str());
    }

    int ctr = 0;
    res.h_val[0] = 1.0;
    while (ctr < maxIter && res.h_val[0] > restol)
    {
        ctr++;
        res.h_val[0] = 0.0;
        Kokkos::deep_copy(res.val,res.h_val);
        
        for(int cs = 0; cs < optiNodes.size(); cs++)
        {              
            OptimiseGPU(derivUtil, nodes, elUtil, res, cs);
            //printf("colorset %cs finished\n", cs);            
        }
                
        Kokkos::deep_copy(res.h_val,res.val);

        Evaluate(derivUtil, nodes, elUtil, res);        
        
        if(m_config["resfile"].beenSet)
        {
            Kokkos::deep_copy(res.h_func,res.func);
            resFile << res.h_val[0] << " " << res.h_worstJac[0] << " " << res.h_func[0] << endl;
        }
        
    }

    if(m_config["histfile"].beenSet)
    {
        ofstream histFile;
        string name = m_config["histfile"].as<string>() + "_end.txt";
        histFile.open(name.c_str());

        Kokkos::deep_copy(elUtil.h_scaledJac,elUtil.scaledJac);

        for(int i = 0; i < dataSet.size(); i++)
        {
            //histFile << dataSet[i]->scaledJac << endl;
            histFile << elUtil.h_scaledJac(i) << endl;
        }
        histFile.close();
    }
    if(m_config["resfile"].beenSet)
    {
        resFile.close();
    }

    t.Stop();
    cout << "Time to compute: " << t.TimePerTest(1) << endl;
    Kokkos::deep_copy(res.h_startInv,res.startInv);
    Kokkos::deep_copy(res.h_worstJac,res.worstJac);
    cout << "Invalid at end:\t\t" << res.h_startInv[0] << endl;
    cout << "Worst at end:\t\t" << res.h_worstJac[0] << endl;

    if(m_config["Kokkos"].beenSet)
    {
        //Kokkos::DefaultExecutionSpace::finalize();
        //Kokkos::DefaultHostExecutionSpace::finalize(); 
        Kokkos::finalize();       
    }
}

} // Utilities
} // Nektar

