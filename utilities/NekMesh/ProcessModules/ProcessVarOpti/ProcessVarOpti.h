

////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessVarOpti.h
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

typedef Kokkos::View<double*,
    Kokkos::DefaultExecutionSpace::scratch_memory_space ,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > ScratchViewType;

enum optimiser
{
    eLinEl,
    eWins,
    eRoca,
    eHypEl
};

struct DerivUtil
{
    NekMatrix<NekDouble> VdmD[3];
    NekMatrix<NekDouble> VdmDL[3]; //deriv matrix without interp
    NekVector<NekDouble> quadW;

    Array<OneD, Array<OneD, NekDouble> > basisDeriv;

    int ptsHigh;
    int ptsLow;
};

typedef boost::shared_ptr<DerivUtil> DerivUtilSharedPtr;

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
    //Kokkos::View<int*> ElmtOffset;
    //typename Kokkos::View< int*>::HostMirror h_ElmtOffset;

    Kokkos::View<int*> coloursetSize;
    typename Kokkos::View<int*>::HostMirror h_coloursetSize;
    Kokkos::View<int**> nElmtArray;
    typename Kokkos::View<int**>::HostMirror h_nElmtArray;
    Kokkos::View<int***> elIdArray;
    typename Kokkos::View<int***>::HostMirror h_elIdArray;
    Kokkos::View<int***> localNodeIdArray;
    typename Kokkos::View<int***>::HostMirror h_localNodeIdArray;
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
    //int startInv;
    Kokkos::View<int[1]> startInv;
    typename Kokkos::View< int[1]>::HostMirror h_startInv;
    //int nReset;
    Kokkos::View<int[1]> nReset;
    typename Kokkos::View< int[1]>::HostMirror h_nReset;

    int n;
    int nDoF;    

    //Kokkos::View<double*> resid;  // for residual of every node in a colourset
};

struct Grad
{
    Kokkos::View<double*[9]> G;
    typename Kokkos::View< double*[9]>::HostMirror h_G;
    Kokkos::View<double*> integral;
    typename Kokkos::View< double*>::HostMirror h_integral;
};


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
    
    // in LoadData.cpp
    void Load_derivUtil(DerivUtilGPU &derivUtil);
    void Load_elUtils(ElUtilGPU &elUtil);
    void Load_nodes(NodesGPU &nodes);
    void Load_residual(Residual &res);

    // in Evaluate.hxx
    void Evaluate(DerivUtilGPU &derivUtil,NodesGPU &nodes, ElUtilGPU &elUtil, Residual &res);
   
    // in Optimise.hxx
    void GetNodeCoordGPU(double (&X)[3], const NodesGPU &nodes,
            Kokkos::View<int***> elIdArray, Kokkos::View<int***> localNodeIdArray, int node, int cs);
    void SetNodeCoordGPU(const double (&X)[3], const NodesGPU &nodes,
            Kokkos::View<int***> elIdArray, Kokkos::View<int***> localNodeIdArray, int nElmt, int node, int cs);
    void GetNodeCoord(double (&X)[3], int id,NodesGPU &nodes,
            typename Kokkos::View<int*>::HostMirror elIdArray, typename Kokkos::View<int*>::HostMirror localNodeIdArray);
    void SetNodeCoord(double (&X)[3], int id,NodesGPU &nodes,
            typename Kokkos::View<int*>::HostMirror elIdArray, typename Kokkos::View<int*>::HostMirror localNodeIdArray, int nElmt);
    double CalcMinJacGPU(const ElUtilGPU &elUtil, int nElmt, int node, int cs, Kokkos::View<int***> elIdArray);
    void Optimise(DerivUtilGPU &derivUtil,NodesGPU &nodes, 
        ElUtilGPU &elUtil, Residual &res, int cs, optimiser opti);

    // in GetFunctional.hxx
    /*template<const int DIM, const bool gradient, const optimiser opti> 
    struct  GetFunctional
    {
         NekDouble operator() (const DerivUtilGPU &derivUtilGPU,
         const NodesGPU &nodes, const ElUtilGPU &elUtil, 
         const Grad &grad, int nElmt, int node, int cs,//const int elId, const int localNodeId,
         const double ep, const member_type &teamMember);
    };*/
    template<const int DIM, const bool gradient> 
    NekDouble GetFunctional (const DerivUtilGPU &derivUtilGPU,
         const NodesGPU &nodes, const ElUtilGPU &elUtil, 
         const Grad &grad, int nElmt, int node, int cs,//const int elId, const int localNodeId,
         const double ep, const member_type &teamMember);
    
    
    // in Hessian.hxx
    template<int DIM> void CalcEValues(const double (&G)[DIM*DIM], double (&eval)[DIM]);
    template<int DIM> int IsIndefinite(const double (&eval)[DIM]);
    template<int DIM> void CalcEVector(const double (&G)[DIM*DIM], const double &eval, double (&evec)[DIM]);

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

    std::map<LibUtilities::ShapeType,DerivUtilSharedPtr> derivUtil;
    optimiser opti;

    static NekDouble c1() {return 1e-3;}
    static NekDouble c2() {return 0.9;}
    static NekDouble gradTol() {return 1e-8;}
    static NekDouble alphaTol() {return 1e-8;}

    
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
  KOKKOS_INLINE_FUNCTION 
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
    if(vect(i) < mn) {
       mn = vect(i);
    }
  }
};

template <typename T>
struct MaxFunctor {
  Kokkos::View<T*> vect;
  KOKKOS_INLINE_FUNCTION 
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




#endif

