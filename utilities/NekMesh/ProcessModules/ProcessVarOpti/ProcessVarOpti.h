

////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.h
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

#ifndef UTILITIES_NEKMESH_PROCESSVAROPTI
#define UTILITIES_NEKMESH_PROCESSVAROPTI

#include "../../Module.h"

#include "ElUtil.h"
//#include "NodeOpti.h"

namespace Nektar
{
namespace Utilities
{

typedef Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace> range_policy;
typedef Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace> range_policy_host;
typedef Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> team_policy;
typedef Kokkos::TeamPolicy<Kokkos::DefaultHostExecutionSpace> team_policy_host;
typedef Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type  member_type;
typedef Kokkos::TeamPolicy<Kokkos::DefaultHostExecutionSpace>::member_type  member_type_host;

typedef Kokkos::DefaultHostExecutionSpace host_space;
typedef Kokkos::DefaultExecutionSpace exe_space;

typedef Kokkos::MemoryTraits<Kokkos::RandomAccess> random_memory;



struct DerivUtil
{
    NekMatrix<NekDouble> VdmD[3];
    NekMatrix<NekDouble> VdmDL[3]; //deriv matrix without interp
    NekVector<NekDouble> quadW;

    Array<OneD, Array<OneD, NekDouble> > basisDeriv;

    int ptsHigh;
    int ptsLow;
};

struct DerivUtilGPU
{
    Kokkos::View<double**,random_memory> VdmDL_0;
    typename Kokkos::View< double**>::HostMirror h_VdmDL_0;
    Kokkos::View<double**,random_memory> VdmDL_1;
    typename Kokkos::View< double**>::HostMirror h_VdmDL_1;        
    Kokkos::View<double**,random_memory> VdmDL_2;
    typename Kokkos::View< double**>::HostMirror h_VdmDL_2;

    Kokkos::View<double**,random_memory> VdmD_0;
    typename Kokkos::View< double**>::HostMirror h_VdmD_0;
    Kokkos::View<double**,random_memory> VdmD_1;
    typename Kokkos::View< double**>::HostMirror h_VdmD_1;        
    Kokkos::View<double**,random_memory> VdmD_2;
    typename Kokkos::View< double**>::HostMirror h_VdmD_2;

    Kokkos::View<double*,random_memory> quadW;
    typename Kokkos::View< double*>::HostMirror h_quadW;

    int nodes_size;

    int ptsHigh;
    int ptsLow;
};


struct ElUtilGPU
{
    Kokkos::View<double**[10],random_memory> idealMap;
    typename Kokkos::View< double**[10]>::HostMirror h_idealMap;
    int ptsHigh;
    int globalnElmt;

    Kokkos::View<double*> minJac;
    typename Kokkos::View< double*>::HostMirror h_minJac;
    Kokkos::View<double*> scaledJac;
    typename Kokkos::View< double*>::HostMirror h_scaledJac;
};

struct NodesGPU
{
    Kokkos::View<double**> X;
    typename Kokkos::View< double**>::HostMirror h_X;
    Kokkos::View<double**> Y;
    typename Kokkos::View< double**>::HostMirror h_Y;        
    Kokkos::View<double**> Z;
    typename Kokkos::View< double**>::HostMirror h_Z;
    Kokkos::View<int**> Id;
    typename Kokkos::View< int**>::HostMirror h_Id;
    int nodes_size;
    int globalnElmt;
    Kokkos::View<int*> ElmtOffset;
    typename Kokkos::View< int*>::HostMirror h_ElmtOffset;

    Kokkos::View<int*> elIdArray;
    typename Kokkos::View<int*>::HostMirror h_elIdArray;
    Kokkos::View<int*> localNodeIdArray;
    typename Kokkos::View<int*>::HostMirror h_localNodeIdArray;
};

typedef boost::shared_ptr<DerivUtil> DerivUtilSharedPtr;

enum optimiser
{
    eLinEl,
    eWins,
    eRoca,
    eHypEl
};

struct Residual
{
    //NekDouble val;
    Kokkos::View<double[1]> val;
    typename Kokkos::View< double[1]>::HostMirror h_val;
    //NekDouble func;
    Kokkos::View<double[1]> func;
    typename Kokkos::View< double[1]>::HostMirror h_func;
    //NekDouble worstJac;
    Kokkos::View<double[1]> worstJac;
    typename Kokkos::View< double[1]>::HostMirror h_worstJac;
    
    int n;
    int nDoF;

    //int startInv;
    Kokkos::View<int[1]> startInv;
    typename Kokkos::View< int[1]>::HostMirror h_startInv;
    //int nReset;
    Kokkos::View<int[1]> nReset;
    typename Kokkos::View< int[1]>::HostMirror h_nReset;     
};

struct Grad
{
    Kokkos::View<double[9]> G;
    typename Kokkos::View< double[9]>::HostMirror h_G;
    Kokkos::View<double[1]> integral;
    typename Kokkos::View< double[1]>::HostMirror h_integral;
};

//typedef boost::shared_ptr<Residual> ResidualSharedPtr;
typedef std::multimap<int,std::pair<int,int>> NodeMap;



class ProcessVarOpti : public ProcessModule
{
public:
    /// Creates an instance of this class
    static boost::shared_ptr<Module> create(MeshSharedPtr m)
    {
        return MemoryManager<ProcessVarOpti>::AllocateSharedPtr(m);
    }
    static ModuleKey className;

    ProcessVarOpti(MeshSharedPtr m);
    virtual ~ProcessVarOpti();

    virtual void Process();

    void Load_derivUtil(DerivUtilGPU &derivUtil);
    template<int DIM> void Load_elUtils(ElUtilGPU &elUtil);
    void Create_nodes_view(NodesGPU &nodes);
    void Load_nodes(NodesGPU &nodes);
    void Evaluate(DerivUtilGPU &derivUtil,NodesGPU &nodes, ElUtilGPU &elUtil, Residual &res);
    void Create_NodeMap(NodesGPU &nodes, std::vector<std::vector<NodeSharedPtr> > &freenodes, NodeMap &nodeMap);

    

    
    
private:
    typedef std::map<int, std::pair<std::vector<int>,
                                    std::vector<ElUtilSharedPtr> > > NodeElMap;
    
    void BuildDerivUtil();
    void GetElementMap();


    std::vector<ElementSharedPtr> GetLockedElements(NekDouble thres);
    std::vector<Array<OneD, NekDouble> > MappingIdealToRef(ElementSharedPtr el);
    std::vector<std::vector<NodeSharedPtr> > GetColouredNodes(std::vector<ElementSharedPtr> elLock, Residual &res);

    NodeElMap nodeElMap;
    std::vector<ElUtilSharedPtr> dataSet;

    NodeMap nodeMap; //maps node_ids to GPU node array

    ResidualSharedPtr res;
    std::map<LibUtilities::ShapeType,DerivUtilSharedPtr> derivUtil;
    optimiser opti;

    std::string ThreadManagerType;
};




}
}


KOKKOS_INLINE_FUNCTION
int num_min(int & mn){return INT_MIN;}
KOKKOS_INLINE_FUNCTION
int num_max(int & mx){return INT_MAX;}
KOKKOS_INLINE_FUNCTION
double num_min(double & mn){return DBL_MIN;}
KOKKOS_INLINE_FUNCTION
double num_max(double & mx){return DBL_MAX;}
KOKKOS_INLINE_FUNCTION
float num_min(float & mn){return FLT_MIN;}
KOKKOS_INLINE_FUNCTION
float num_max(float & mx){return FLT_MAX;}


template <typename T>
struct MinFunctor {
  Kokkos::View<T*> vect;  
  MinFunctor(const Kokkos::View<T*> vect_):
    vect(vect_) {}
  KOKKOS_INLINE_FUNCTION
  void init(T& mn )const {
    mn = num_max(mn);
  }
  KOKKOS_INLINE_FUNCTION
  void join(volatile T& mn, const volatile T& update) const {
    if(update < mn) {
      mn = update;
    }
  }
  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i, T& mn) const {
    const T value = vect(i);
    if(value < mn) {
       mn = value;
    }
  }
};

template <typename T>
struct MaxFunctor {
  Kokkos::View<T*> vect;
  MaxFunctor(const Kokkos::View<T*> vect_):
    vect(vect_) {}
  KOKKOS_INLINE_FUNCTION
  void init(T& mx) const {
    mx = num_min(mx);
  }
  KOKKOS_INLINE_FUNCTION
  void join(volatile T& mx, const volatile T& update) const {
    if(update > mx) {
      mx = update;
    }
  }
  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i, T& mx) const {
    const T value = vect(i);
    if(value > mx) {
       mx = value;
    }
  }
};

template <typename T>
void process()
{
  // initialise some vector
  T nset = 1024;
  Kokkos::View< T*> points("Points",nset);
  typename Kokkos::View<T*>::HostMirror h_points = Kokkos::create_mirror_view(points);
  srand(56779);
  for(int i = 0; i < nset; i++)
  {
      //h_points(i) = 1.0*i;
      h_points(i) = rand()%(1024);
  }
  // compute sum min and max of vector on CPU for comparison
  T h_sm = 0;
  T h_mn = num_max(h_mn);
  T h_mx = num_min(h_mx);
  for(int i = 0; i < nset; i++)
  {
      h_sm += h_points[i];
      h_mn = (h_mn > h_points[i] ? h_points[i] : h_mn);
      h_mx = (h_mx < h_points[i] ? h_points[i] : h_mx);
  }
  // copy vector to GPU
  Kokkos::deep_copy(points,h_points);
  typedef Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace> range_policy;

  // calculate sum of vector on GPU
  T sm = 0.0;
  Kokkos::parallel_reduce(range_policy(0, nset), KOKKOS_LAMBDA (const int& i, T& sum)
  {
      sum += points(i);
  }, sm);  

  // atomic add of vector on GPU
  Kokkos::View< T[1]> sm_atomic("atomic");
  typename Kokkos::View< T*>::HostMirror h_sm_atomic = Kokkos::create_mirror_view(sm_atomic);  
  h_sm_atomic[0] = 0.0;
  Kokkos::deep_copy(sm_atomic,h_sm_atomic);
  
  Kokkos::parallel_for("summation", range_policy(0,nset), KOKKOS_LAMBDA (const int& i)  {
      Kokkos::atomic_add(&sm_atomic[0], points(i));
  });
  Kokkos::deep_copy(h_sm_atomic,sm_atomic);
  std::cout << "h_sm = " << h_sm << ", sm = " << sm << ", h_sm_atomic = " << h_sm_atomic[0] <<std::endl;

  //calculate max of vector on GPU
  T mx;    
  MaxFunctor <T> mxfunctor(points);
  Kokkos::parallel_reduce(range_policy(0, nset) , mxfunctor, mx);
  std::cout << "h_mx = " << h_mx << ", mx = " << mx <<std::endl;

  //calculate min of vector on GPU
  T mn;
  MinFunctor <T> mnfunctor(points);
  Kokkos::parallel_reduce( range_policy(0, nset) , mnfunctor, mn);
  std::cout << "h_mn = " << h_mn << ", mn = " << mn <<std::endl;  
}


#endif

