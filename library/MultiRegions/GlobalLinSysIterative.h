///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysIterative.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: GlobalLinSysIterative header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSITERATIVE_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSITERATIVE_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/Preconditioner.h>

#include <boost/circular_buffer.hpp>

#include <Kokkos_Core.hpp>


namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;

        /// A global linear system.
        class GlobalLinSysIterative : virtual public GlobalLinSys
        {
        public:
            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysIterative(
                    const GlobalLinSysKey                &pKey,
                    const boost::weak_ptr<ExpList>       &pExpList,
                    const boost::shared_ptr<AssemblyMap> &pLocToGloMap);

            MULTI_REGIONS_EXPORT virtual ~GlobalLinSysIterative();


                        //Kokkos
            typedef Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace> range_policy_host;
            typedef Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace> range_policy;
            typedef Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> team_policy;
            typedef Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type  member_type;
            typedef Kokkos::MemoryTraits<Kokkos::RandomAccess> random_memory;
            typedef Kokkos::View<double*,
                Kokkos::DefaultExecutionSpace::scratch_memory_space ,
                Kokkos::MemoryTraits<Kokkos::Unmanaged> > ScratchViewType;


            std::vector<std::vector<int> > CreateColourSets(
                Array<OneD, const int> &localToGlobalMap,
                int ncoeffs, int elmts);


            // functions for plain parallel Conjugate Gradient
            void GeneralMatrixOp_plain(
                    const Array<OneD,const NekDouble> &inarray,
                    Array<OneD,      NekDouble> &outarray,
                    const NekDouble &lambda,
                    const Array<OneD, const NekDouble> &quadMetricGlo,                
                    const Array<OneD, const NekDouble> &laplacian00Glo,
                    const Array<OneD, const NekDouble> &laplacian01Glo,
                    const Array<OneD, const NekDouble> &laplacian11Glo,
                    const int &nquad0, const int &nquad1, 
                    const int &nmodes0, const int &nmodes1, const int &ncoeffs, 
                    const Array<OneD, const int>  &coeff_offset, const int &elmts,
                    const Array<OneD, const NekDouble> &base0,
                    const Array<OneD, const NekDouble> &base1,
                    const Array<OneD, const NekDouble> &dbase0,
                    const Array<OneD, const NekDouble> &dbase1,
                    const DNekMatSharedPtr &D0, const DNekMatSharedPtr &D1,
                    const int &numLocalCoeffs, const int &numGlobalCoeffs,
                    const Array<OneD, const int> &localToGlobalMap,
                    const Array<OneD, const NekDouble> &localToGlobalSign,
                    const int iteration);


            void GeneralMatrixOp_IterPerExp_plain(
                    const Array<OneD,const NekDouble> &inarray,
                    //Array<OneD,      NekDouble> &outarray,
                    Kokkos::View<double*, Kokkos::HostSpace> transfer_out,                     
                    const NekDouble &lambda,
                    const Array<OneD, const NekDouble> &quadMetricGlo,                
                    const Array<OneD, const NekDouble> &laplacian00Glo,
                    const Array<OneD, const NekDouble> &laplacian01Glo,
                    const Array<OneD, const NekDouble> &laplacian11Glo,
                    const int &nquad0, const int &nquad1, 
                    const int &nmodes0, const int &nmodes1, const int &ncoeffs, 
                    const Array<OneD, const int>  &coeff_offset, const int &elmts,
                    const Array<OneD, const NekDouble> &base0,
                    const Array<OneD, const NekDouble> &base1,
                    const Array<OneD, const NekDouble> &dbase0,
                    const Array<OneD, const NekDouble> &dbase1,
                    const DNekMatSharedPtr &D0, const DNekMatSharedPtr &D1);

            void HelmholtzMatrixOp_MatFree_plain(
                    const Array<OneD, const NekDouble> &inarray,
                    //      Array<OneD,       NekDouble> &outarray,
                    Kokkos::View<double*, Kokkos::HostSpace> transfer_out, 
                    const int &el, const Array<OneD, const int>  &coeff_offset,
                    const NekDouble &lambda,
                    const Array<OneD, const NekDouble> &quadMetric,
                    const Array<OneD, const NekDouble> &laplacian00,
                    const Array<OneD, const NekDouble> &laplacian01,
                    const Array<OneD, const NekDouble> &laplacian11,
                    const int &nquad0, const int &nquad1,
                    const int &nmodes0, const int &nmodes1, const int &ncoeffs,
                    const Array<OneD, const NekDouble> &base0,
                    const Array<OneD, const NekDouble> &base1,
                    const Array<OneD, const NekDouble> &dbase0,
                    const Array<OneD, const NekDouble> &dbase1,
                    const DNekMatSharedPtr &D0, const DNekMatSharedPtr &D1);

            void IProductWRTBase_SumFacKernel_plain(
                    const Array<OneD, const NekDouble>& base0,
                    const Array<OneD, const NekDouble>& base1,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD,       NekDouble> &outarray,
                    Array<OneD, NekDouble> &wsp,
                    const int &nmodes0, const int &nmodes1,
                    const int &nquad0, const int &nquad1);


            void BwdTrans_SumFacKernel_plain(
                    const Array<OneD, const NekDouble>& base0,
                    const Array<OneD, const NekDouble>& base1,
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &outarray,
                    Array<OneD, NekDouble> &wsp,
                    const int &nmodes0, const int &nmodes1,
                    const int &nquad0, const int &nquad1);

            void PhysTensorDeriv_plain(
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &outarray_d0,
                    Array<OneD, NekDouble> &outarray_d1,
                    const int &nquad0, const int &nquad1,
                    const DNekMatSharedPtr &D0, const DNekMatSharedPtr &D1);




            // functions for full Kokkos Conjugate Gradient
            void DoConjugateGradient_Kokkos(
                    const int                          nGlobal,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr        &plocToGloMap,
                    const int                          nDir);


            void GeneralMatrixOp_Kokkos(
                const Kokkos::View<double*> inarray,
                Kokkos::View<double*> outarray,
                const Kokkos::View<double[1]> lambda,
                const Kokkos::View<double*> quadMetricGlo,
                const Kokkos::View<double*> laplacian00Glo,
                const Kokkos::View<double*> laplacian01Glo,
                const Kokkos::View<double*> laplacian11Glo,
                const int &nquad0, const int &nquad1,
                const int &nmodes0, const int &nmodes1, const int &ncoeffs, 
                const Kokkos::View<int*> coeff_offset,
                const int &elmts,
                const Kokkos::View<double*> base0,
                const Kokkos::View<double*> base1,
                const Kokkos::View<double*> dbase0,
                const Kokkos::View<double*> dbase1,
                const Kokkos::View<double*> D0,
                const Kokkos::View<double*> D1,
                const int &numLocalCoeffs, const int &numGlobalCoeffs,
                const Kokkos::View<int*> localToGlobalMap,
                const Kokkos::View<double*> localToGlobalSign,
                const int iteration,
                Kokkos::View<double*> transfer_in,
                Kokkos::View<double*> transfer_out);


            void GeneralMatrixOp_IterPerExp_Kokkos(
                const Kokkos::View<double*> transfer_in,
                Kokkos::View<double*> transfer_out,
                const Kokkos::View<double[1]> lambda,
                const Kokkos::View<double*> quadMetricGlo,
                const Kokkos::View<double*> laplacian00Glo,
                const Kokkos::View<double*> laplacian01Glo,
                const Kokkos::View<double*> laplacian11Glo,
                const int &nquad0, const int &nquad1,
                const int &nmodes0, const int &nmodes1, const int &ncoeffs, 
                const Kokkos::View<int*> coeff_offset,
                const int &elmts,
                const Kokkos::View<double*> base0,
                const Kokkos::View<double*> base1,
                const Kokkos::View<double*> dbase0,
                const Kokkos::View<double*> dbase1,
                const Kokkos::View<double*> D0,
                const Kokkos::View<double*> D1);

            KOKKOS_INLINE_FUNCTION
            void HelmholtzMatrixOp_MatFree_Kokkos(
                const ScratchViewType s_tmp_inarray,
                Kokkos::View<double*> outarray,
                const int &el,
                const Kokkos::View<int*>  coeff_offset,
                const Kokkos::View<double[1]> lambda,
                const Kokkos::View<double*>  quadMetricGlo,
                const Kokkos::View<double*> laplacian00Glo,
                const Kokkos::View<double*> laplacian01Glo,
                const Kokkos::View<double*> laplacian11Glo,
                const int &nquad0, const int &nquad1,
                const int &nmodes0, const int &nmodes1, const int &ncoeffs,
                const Kokkos::View<double*> base0,
                const Kokkos::View<double*> base1,
                const Kokkos::View<double*> dbase0,
                const Kokkos::View<double*> dbase1,
                const Kokkos::View<double*> D0,
                const Kokkos::View<double*> D1,
                const member_type &teamMember, const int &wspsize);

            KOKKOS_INLINE_FUNCTION
            void IProductWRTBase_SumFacKernel_Kokkos(
                const Kokkos::View<double*> base0,
                const Kokkos::View<double*> base1,
                const ScratchViewType inarray,
                ScratchViewType outarray,
                ScratchViewType wsp,
                const int &nmodes0, const int &nmodes1,
                const int &nquad0, const int &nquad1);

            KOKKOS_INLINE_FUNCTION
            void BwdTrans_SumFacKernel_Kokkos(
                const Kokkos::View<double*> base0,
                const Kokkos::View<double*> base1,
                const ScratchViewType inarray,
                ScratchViewType outarray,
                ScratchViewType wsp,
                const int &nmodes0, const int &nmodes1,
                const int &nquad0, const int &nquad1);

            KOKKOS_INLINE_FUNCTION
            void PhysTensorDeriv_Kokkos(
                const ScratchViewType inarray,
                ScratchViewType outarray_d0,
                ScratchViewType outarray_d1,
                const int &nquad0, const int &nquad1,
                const Kokkos::View<double*> D0,
                const Kokkos::View<double*> D1);


            /*KOKKOS_INLINE_FUNCTION
            double plainDdot(int n, const double *dx, int incx,
                    const double *dy, int incy);

            KOKKOS_INLINE_FUNCTION
            int plainDaxpy(int n, const double da, const double *dx,
                    int incx, double *dy, int incy);

            KOKKOS_INLINE_FUNCTION
            int plainDgemm(char transa, char transb, int m, int n, int k,
                    const double alpha, const double *a, int lda, const double *b,
                    int ldb, const double beta, double *c, int ldc);

            KOKKOS_INLINE_FUNCTION
            int plainDgemv(char trans, int m, int n, const double alpha,
                    const double *a, int lda, const double *x, int incx,
                    const double beta, double *y, int incy);*/



        protected:
            /// Global to universal unique map
            Array<OneD, int>                            m_map;

            /// maximum iterations
            int                                         m_maxiter;

            /// Tolerance of iterative solver.
            NekDouble                                   m_tolerance;

            /// dot product of rhs to normalise stopping criterion
            NekDouble                                   m_rhs_magnitude;

            /// cnt to how many times rhs_magnitude is called 
            NekDouble                                   m_rhs_mag_sm; 
            
            PreconditionerSharedPtr                     m_precon;

            MultiRegions::PreconditionerType            m_precontype;
            
            int                                         m_totalIterations;

            /// Whether to apply projection technique
            bool                                        m_useProjection;

            /// Root if parallel
            bool                                        m_root;

            /// Storage for solutions to previous linear problems
            boost::circular_buffer<Array<OneD, NekDouble> > m_prevLinSol;

            /// Total counter of previous solutions
            int m_numPrevSols;

            /// A-conjugate projection technique
            void DoAconjugateProjection(
                    const int pNumRows,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &locToGloMap,
                    const int pNumDir);

            /// Actual iterative solve
            void DoConjugateGradient(
                    const int pNumRows,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &locToGloMap,
                    const int pNumDir);

            void DoConjugateGradient_plain(
                    const int                          nGlobal,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr        &plocToGloMap,
                    const int                          nDir);

            void DoConjugateGradient_OpenMP(
                    const int                          nGlobal,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr        &plocToGloMap,
                    const int                          nDir);


            void Set_Rhs_Magnitude(const NekVector<NekDouble> &pIn);

            virtual void v_UniqueMap() = 0;
            
        private:
            void UpdateKnownSolutions(
                    const int pGlobalBndDofs,
                    const Array<OneD,const NekDouble> &pSolution,
                    const int pNumDirBndDofs);

            NekDouble CalculateAnorm(
                    const int nGlobal,
                    const Array<OneD,const NekDouble> &in,
                    const int nDir);


            /// Solve the matrix system
            virtual void v_SolveLinearSystem(
                    const int pNumRows,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &locToGloMap,
                    const int pNumDir);

            virtual void v_DoMatrixMultiply(
                    const Array<OneD, NekDouble>& pInput,
                          Array<OneD, NekDouble>& pOutput) = 0;

            // functions for parallel Conjugate Gradient
            void GetMatrixMultiplyMetrics(
                Array<OneD, NekDouble> &quadMetricGlo,                
                Array<OneD, NekDouble> &laplacian00Glo,
                Array<OneD, NekDouble> &laplacian01Glo,
                Array<OneD, NekDouble> &laplacian11Glo,
                int &nquad0, int &nquad1, int &elmts,
                int &numLocalCoeffs, int &numGlobalCoeffs,
                Array<OneD, const int> &localToGlobalMap,
                Array<OneD, const NekDouble> &localToGlobalSign);






            // functions for OpenMP Conjugate Gradient
            void GeneralMatrixOp_OpenMP(
                    const Array<OneD,const NekDouble> &inarray,
                    Array<OneD,      NekDouble> &outarray,
                    const NekDouble &lambda,
                    const Array<OneD, const NekDouble> &quadMetricGlo,                
                    const Array<OneD, const NekDouble> &laplacian00Glo,
                    const Array<OneD, const NekDouble> &laplacian01Glo,
                    const Array<OneD, const NekDouble> &laplacian11Glo,
                    const int &nquad0, const int &nquad1, 
                    const int &nmodes0, const int &nmodes1, const int &ncoeffs, 
                    const Array<OneD, const int>  &coeff_offset, const int &elmts,
                    const Array<OneD, const NekDouble> &base0,
                    const Array<OneD, const NekDouble> &base1,
                    const Array<OneD, const NekDouble> &dbase0,
                    const Array<OneD, const NekDouble> &dbase1,
                    const DNekMatSharedPtr &D0, const DNekMatSharedPtr &D1,
                    const int &numLocalCoeffs, const int &numGlobalCoeffs,
                    const Array<OneD, const int> &localToGlobalMap,
                    const Array<OneD, const NekDouble> &localToGlobalSign);


            void GeneralMatrixOp_IterPerExp_OpenMP(
                    const Array<OneD,const NekDouble> &inarray,
                    Array<OneD,      NekDouble> &outarray,
                    const NekDouble &lambda,
                    const Array<OneD, const NekDouble> &quadMetricGlo,                
                    const Array<OneD, const NekDouble> &laplacian00Glo,
                    const Array<OneD, const NekDouble> &laplacian01Glo,
                    const Array<OneD, const NekDouble> &laplacian11Glo,
                    const int &nquad0, const int &nquad1, 
                    const int &nmodes0, const int &nmodes1, const int &ncoeffs, 
                    const Array<OneD, const int>  &coeff_offset, const int &elmts,
                    const Array<OneD, const NekDouble> &base0,
                    const Array<OneD, const NekDouble> &base1,
                    const Array<OneD, const NekDouble> &dbase0,
                    const Array<OneD, const NekDouble> &dbase1,
                    const DNekMatSharedPtr &D0, const DNekMatSharedPtr &D1);

            void HelmholtzMatrixOp_MatFree_OpenMP(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const NekDouble &lambda,
                    const Array<OneD, const NekDouble> &quadMetric,
                    const Array<OneD, const NekDouble> &laplacian00,
                    const Array<OneD, const NekDouble> &laplacian01,
                    const Array<OneD, const NekDouble> &laplacian11,
                    const int &nquad0, const int &nquad1,
                    const int &nmodes0, const int &nmodes1, const int &ncoeffs,
                    const Array<OneD, const NekDouble> &base0,
                    const Array<OneD, const NekDouble> &base1,
                    const Array<OneD, const NekDouble> &dbase0,
                    const Array<OneD, const NekDouble> &dbase1,
                    const DNekMatSharedPtr &D0, const DNekMatSharedPtr &D1);

            void IProductWRTBase_SumFacKernel_OpenMP(
                    const Array<OneD, const NekDouble>& base0,
                    const Array<OneD, const NekDouble>& base1,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD,       NekDouble> &outarray,
                    Array<OneD, NekDouble> &wsp,
                    const int &nmodes0, const int &nmodes1,
                    const int &nquad0, const int &nquad1);


            void BwdTrans_SumFacKernel_OpenMP(
                    const Array<OneD, const NekDouble>& base0,
                    const Array<OneD, const NekDouble>& base1,
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &outarray,
                    Array<OneD, NekDouble> &wsp,
                    const int &nmodes0, const int &nmodes1,
                    const int &nquad0, const int &nquad1);

            void PhysTensorDeriv_OpenMP(
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &outarray_d0,
                    Array<OneD, NekDouble> &outarray_d1,
                    const int &nquad0, const int &nquad1,
                    const DNekMatSharedPtr &D0, const DNekMatSharedPtr &D1);
            

        };
    }
}

#endif
