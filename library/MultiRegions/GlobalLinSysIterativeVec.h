///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysIterativeVec.h
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
// Description: GlobalLinSysIterativeVec header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSITERATIVEVEC_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSITERATIVEVEC_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/Preconditioner.h>

#include <boost/circular_buffer.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;

        /// A global linear system.
        class GlobalLinSysIterativeVec : virtual public GlobalLinSys
        {
        public:
            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysIterativeVec(
                    const GlobalLinSysKey                &pKey,
                    const Array<OneD, std::weak_ptr<ExpList> >&pExpList,
                    const Array<OneD, std::shared_ptr<AssemblyMap> >
                    &pLocToGloMap);

            MULTI_REGIONS_EXPORT virtual ~GlobalLinSysIterativeVec();

        protected:
            /// Local Matrix System/Expandiosn
            const Array<OneD, std::weak_ptr<ExpList>>   m_expListVec;

            /// Global to universal unique map
            Array<OneD, Array<OneD, int> >              m_mapVec;

            /// maximum iterations
            int                                         m_maxiter;

            /// Tolerance of iterative solver.
            NekDouble                                   m_tolerance;

            /// dot product of rhs to normalise stopping criterion
            NekDouble                                   m_rhs_magnitude;

            /// cnt to how many times rhs_magnitude is called 
            NekDouble                                   m_rhs_mag_sm; 
            
            Array<OneD, PreconditionerSharedPtr>        m_preconVec;

            MultiRegions::PreconditionerType            m_precontype;
            
            int                                         m_totalIterations;

            /// Whether to apply projection technique
            bool                                        m_useProjection;

            /// Root if parallel
            bool                                        m_root;

            /// Actual iterative solve
            void DoConjugateGradient(
                const Array<OneD, int >                    &nGlobal,
                const Array<OneD, Array<OneD, NekDouble> > &pvecInput,
                      Array<OneD, Array<OneD, NekDouble> > &pvecOutput,
                const Array<OneD, AssemblyMapSharedPtr >   &pvecLocToGloMap,
                const Array<OneD, int >                    &nDir);
            
            void Set_Rhs_Magnitude( const int ntest,
                                    const Array<OneD, Array<OneD, NekDouble> >
                                    &pIn);

            virtual void v_UniqueMap() = 0;
            
        private:
            /// Solve the matrix system
            virtual void v_SolveVecLinearSystem(
                    const Array<OneD, int >& nGlobal,
                    const Array<OneD, const Array<OneD, NekDouble> > &pvecInput,
                    Array <OneD, Array<OneD, NekDouble> > &pvecOutput,
                    const Array<OneD, AssemblyMapSharedPtr > &pvecLocToGloMap,
                    const Array<OneD, int> & nDir);

            virtual void v_DoMatrixMultiply(
                            const Array<OneD, int> nglobal,
                            const Array<OneD, NekDouble>& pInput,
                            Array<OneD, NekDouble>& pOutput) = 0;
        };
    }
}

#endif
