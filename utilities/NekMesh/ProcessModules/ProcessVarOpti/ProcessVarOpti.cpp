

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
        
        
    }
    typedef Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace> range_policy_host;

    Kokkos::parallel_for(range_policy_host(0,dataSet.size()), KOKKOS_LAMBDA (const int i)
    {
        dataSet[i]->Evaluate();
    });

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


    /*
    if(m_config["Kokkos"].beenSet)
    {
        // Initialise Kokkos
        Kokkos::InitArguments args;
        args.num_threads = m_config["numthreads"].as<int>();
        Kokkos::DefaultHostExecutionSpace::initialize(args.num_threads);       
        Kokkos::DefaultExecutionSpace::initialize();

        int N = pow(2,12) ;       // number of rows 2^12
        int M = pow(2,10);       // number of columns 2^10
        int S = pow(2,22) ;      // total size 2^22
        int nrepeat = 100 ;    // number of repeats of the test
        
        // typedef Kokkos::Serial   ExecSpace ;
        // typedef Kokkos::Threads  ExecSpace ;
        // typedef Kokkos::OpenMP   ExecSpace ;
        typedef Kokkos::Cuda        ExecSpace ;
        
        // typedef Kokkos::HostSpace    MemSpace; 
        // typedef Kokkos::OpenMP       MemSpace; 
        typedef Kokkos::CudaSpace       MemSpace;
        // typedef Kokkos::CudaUVMSpace MemSpace;

        typedef Kokkos::LayoutLeft   Layout ;
        // typedef Kokkos::LayoutRight  Layout ;

        typedef Kokkos::RangePolicy<ExecSpace> range_policy ;

        // Allocate y, x vectors and Matrix A:
        // Device
        typedef Kokkos::View<double*, Layout, MemSpace>   ViewVectorType;
        typedef Kokkos::View<double**, Layout, MemSpace>   ViewMatrixType;
        ViewVectorType devy("devy", N);
        ViewVectorType devx("devx", M);
        ViewMatrixType devA("devA", N, M);

        //Host mirror
        ViewVectorType::HostMirror y =  Kokkos::create_mirror_view(devy);
        ViewVectorType::HostMirror x =  Kokkos::create_mirror_view(devx);
        ViewMatrixType::HostMirror A =  Kokkos::create_mirror_view(devA);


        // Initialize y vector on host
        for (int i = 0; i < N; ++i) {
          y( i ) = 1;
        }

        // Initialize x vector on host
        for (int i = 0; i < M; ++i) {
          x( i ) = 1;
        }

        // Initialize A matrix, note 2D indexing computation on host
        for (int j = 0; j < N; ++j) {
          for ( int i = 0 ; i < M ; ++i ) {
            A( j , i ) = 1;
          }
        }

        //Deep copy host view to device views
        Kokkos::deep_copy(devy, y);
        Kokkos::deep_copy(devx, x);
        Kokkos::deep_copy(devA, A);

        // Timer products
        struct timeval begin,end;
        gettimeofday(&begin,NULL);

        for ( int repeat = 0; repeat < nrepeat; repeat++)
        {
            //Application: <y,Ax> = y^T*A*x
            double result = 0;
            Kokkos::parallel_reduce( range_policy( 0, N ), KOKKOS_LAMBDA ( int j, double &update ) {
              double temp2 = 0;
              for ( int i = 0 ; i < M ; ++i ) {
                temp2 += devA( j , i ) * devx( i );
              }
              update += devy( j ) * temp2;
            }, result );

            //Output result
            if ( repeat == (nrepeat - 1) )
              printf("  Computed result for %d x %d is %lf\n", N, M, result);
            const double solution = (double)N *(double)M;

            if ( result != solution ) {
              printf("  Error: result( %lf ) != solution( %lf )\n",result,solution);
            }
        }

        gettimeofday(&end,NULL);

        // Calculate time
        double time = 1.0*(end.tv_sec-begin.tv_sec) +
                        1.0e-6*(end.tv_usec-begin.tv_usec);

        // Calculate bandwidth.
        double Gbytes = 1.0e-9 * double(sizeof(double) * ( M + M * N + N )) ;

        // Print results (problem size, time and bandwidth in GB/s)
        printf("  M( %d ) N( %d ) nrepeat ( %d ) problem( %g MB ) time( %g s ) bandwidth( %g GB/s )\n",
                 M , N, nrepeat, Gbytes * 1000, time, Gbytes * nrepeat / time );

        
    }
    */

    Timer t;
    t.Start();

    int ctr = 0;
    ofstream resFile;
    if(m_config["resfile"].beenSet)
    {
        resFile.open(m_config["resfile"].as<string>().c_str());
    }

    while (res->val > restol)
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

        
        Kokkos::parallel_for(range_policy_host(0,dataSet.size()), KOKKOS_LAMBDA (const int i)
        {
            dataSet[i]->Evaluate();
        });


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
        if(ctr >= maxIter)
            break;
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

}
}

