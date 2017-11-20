///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysIterativeStaticCondVec.h
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
// Description: GlobalLinSysIterativeStaticCondVec header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSITERATIVESTATICCONDVEC_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSITERATIVESTATICCONDVEC_H

#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/GlobalLinSysIterative.h>
#include <MultiRegions/GlobalLinSysStaticCondVec.h>
#include <LibUtilities/LinearAlgebra/SparseMatrixFwd.hpp>


namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;
        class GlobalLinSysIterativeStaticCondVec;

        typedef std::shared_ptr<GlobalLinSysIterativeStaticCondVec>
            GlobalLinSysIterativeStaticCondVecSharedPtr;

        enum LocalMatrixStorageStrategy
        {
            eNoStrategy,
            eContiguous,
            eNonContiguous,
            eSparse
        };

        const char* const LocalMatrixStorageStrategyMap[] =
        {
            "Contiguous",
            "Non-contiguous",
            "Sparse"
        };


        /// A global linear system.
        class GlobalLinSysIterativeStaticCondVec : virtual public GlobalLinSysIterativeVec,
            virtual public GlobalLinSysStaticCondVec
        {
        public:
            typedef NekSparseDiagBlkMatrix<StorageSmvBsr<NekDouble> >
                                            DNekSmvBsrDiagBlkMat;
            typedef std::shared_ptr<DNekSmvBsrDiagBlkMat>
                                            DNekSmvBsrDiagBlkMatSharedPtr;

            /// Name of class
            static std::string className;
            static std::string className2;

            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysIterativeStaticCondVec(
                const GlobalLinSysKey                &mkey,
                const Array<OneD, std::weak_ptr<ExpList> > &pVecExpList,
                const Array<OneD, std::shared_ptr<AssemblyMap> > &pVecLocToGloMap);

            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysIterativeStaticCondVec(
                const GlobalLinSysKey                &mkey,
                const Array<OneD, std::weak_ptr<ExpList> > &pVecExpList,
                const DNekScalBlkMatSharedPtr         pSchurCompl,
                const DNekScalBlkMatSharedPtr         pBinvD,
                const DNekScalBlkMatSharedPtr         pC,
                const DNekScalBlkMatSharedPtr         pInvD,
                const Array<OneD, std::shared_ptr<AssemblyMap> > &pVecLocToGloMap,
                const Array<OneD, PreconditionerSharedPtr>       &pPreconVec);

            virtual ~GlobalLinSysIterativeStaticCondVec();

        protected:
            virtual DNekScalBlkMatSharedPtr v_GetStaticCondBlock(unsigned int n);

            virtual DNekScalBlkMatSharedPtr v_PreSolve(
                int                     scLevel,
                Array<OneD, Array<OneD, NekDouble> >  &F);
            virtual void v_BasisTransformLoc(Array<OneD, NekDouble>& pInOut);
            virtual void v_BasisInvTransformLoc(Array<OneD, NekDouble>& pInOut);

        private:
            DNekScalBlkMatSharedPtr                  m_S1Blk;
            /// Dense storage for block Schur complement matrix
            std::vector<double>                      m_storage;
            /// Vector of pointers to local matrix data
            std::vector<const double*>               m_denseBlocks;
            /// Ranks of local matrices
            Array<OneD, unsigned int>                m_rows;
            /// Scaling factors for local matrices
            Array<OneD, NekDouble>                   m_scale;
            /// Sparse representation of Schur complement matrix at this level
            DNekSmvBsrDiagBlkMatSharedPtr            m_sparseSchurCompl;
            /// Utility strings
            static std::string                       storagedef;
            static std::string                       storagelookupIds[];

            virtual void v_InitObject();

            /// Assemble the Schur complement matrix.
            void v_AssembleSchurComplement(
                const std::shared_ptr<AssemblyMap> locToGloMap);

            /// Prepares local representation of Schur complement
            /// stored as a sparse block-diagonal matrix.
            void PrepareLocalSchurComplement();

            /// Perform a Shur-complement matrix multiply operation.
            virtual void v_DoMatrixMultiply(
                            const Array<OneD, int> nglobal,
                            const Array<OneD, NekDouble>& pInput,
                            Array<OneD, NekDouble>& pOutput);

            virtual void v_UniqueMap();
        };
    }
}

#endif
