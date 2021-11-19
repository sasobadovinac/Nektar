///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysIterativeStaticCondVec.cpp
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
// Description: Implementation to linear solver using single-
//              or multi-level static condensation
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysIterativeStaticCondVec.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/StorageSmvBsr.hpp>
#include <LibUtilities/LinearAlgebra/SparseDiagBlkMatrix.hpp>
#include <LibUtilities/LinearAlgebra/SparseUtils.hpp>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysIterativeStaticCondVec
         *
         * Solves a linear system iteratively using single- or multi-level
         * static condensation.
         */


        std::string GlobalLinSysIterativeStaticCondVec::storagedef = 
            LibUtilities::SessionReader::RegisterDefaultSolverInfo(
                "LocalMatrixStorageStrategy",
                "Sparse");
        std::string GlobalLinSysIterativeStaticCondVec::storagelookupIds[3] = {
            LibUtilities::SessionReader::RegisterEnumValue(
                "LocalMatrixStorageStrategy",
                "Contiguous",
                MultiRegions::eContiguous),
            LibUtilities::SessionReader::RegisterEnumValue(
                "LocalMatrixStorageStrategy",
                "Non-contiguous",
                MultiRegions::eNonContiguous),
            LibUtilities::SessionReader::RegisterEnumValue(
                "LocalMatrixStorageStrategy",
                "Sparse",
                MultiRegions::eSparse),
        };

        /**
         * For a matrix system of the form @f[
         * \left[ \begin{array}{cc}
         * \boldsymbol{A} & \boldsymbol{B}\\
         * \boldsymbol{C} & \boldsymbol{D}
         * \end{array} \right]
         * \left[ \begin{array}{c} \boldsymbol{x_1}\\ \boldsymbol{x_2}
         * \end{array}\right]
         * = \left[ \begin{array}{c} \boldsymbol{y_1}\\ \boldsymbol{y_2}
         * \end{array}\right],
         * @f]
         * where @f$\boldsymbol{D}@f$ and
         * @f$(\boldsymbol{A-BD^{-1}C})@f$ are invertible, store and assemble
         * a static condensation system, according to a given local to global
         * mapping. #m_linSys is constructed by AssembleSchurComplement().
         * @param   mKey        Associated matrix key.
         * @param   pLocMatSys  LocalMatrixSystem
         * @param   locToGloMap Local to global mapping.
         */
        GlobalLinSysIterativeStaticCondVec::GlobalLinSysIterativeStaticCondVec(
            const GlobalLinSysKey                &pKey,
            const Array<OneD, std::weak_ptr<ExpList> > &pVecExpList,
            const Array<OneD, std::shared_ptr<AssemblyMap> > &pVecLocToGloMap)
                : GlobalLinSys(pKey, pVecExpList, pVecLocToGloMap[0]),
                  GlobalLinSysIterativeVec(pKey,pVecExpList,pVecLocToGloMap),
                  GlobalLinSysStaticCondVec(pKey,pVecExpList,pVecLocToGloMap)
        {
            ASSERTL1((pKey.GetGlobalSysSolnType()==eIterativeStaticCond),
                     "This constructor is only valid when using static "
                     "condensation");
            ASSERTL1(pKey.GetGlobalSysSolnType()
                        == pVecLocToGloMap[0]->GetGlobalSysSolnType(),
                     "The local to global map is not set up for the requested "
                     "solution type");
        }

        
        /**
         *
         */
        GlobalLinSysIterativeStaticCondVec::GlobalLinSysIterativeStaticCondVec(
            const GlobalLinSysKey                &pKey,
            const Array<OneD, std::weak_ptr<ExpList> > &pVecExpList,
            const DNekScalBlkMatSharedPtr         pSchurCompl,
            const DNekScalBlkMatSharedPtr         pBinvD,
            const DNekScalBlkMatSharedPtr         pC,
            const DNekScalBlkMatSharedPtr         pInvD,
            const Array<OneD, std::shared_ptr<AssemblyMap> > &pVecLocToGloMap,
            const Array<OneD, PreconditionerSharedPtr>       &pPreconVec)
            : GlobalLinSys          (pKey, pVecExpList, pVecLocToGloMap[0]),
              GlobalLinSysIterativeVec (pKey, pVecExpList, pVecLocToGloMap),
              GlobalLinSysStaticCondVec(pKey, pVecExpList, pVecLocToGloMap)
        {
            m_schurCompl  = pSchurCompl;
            m_S1Blk       = pSchurCompl;
            m_BinvD       = pBinvD;
            m_C           = pC;
            m_invD        = pInvD;
            m_preconVec   = pPreconVec;
        }


        void GlobalLinSysIterativeStaticCondVec::v_InitObject()
        {
            int nvec = m_expListVec.num_elements();
            
            m_preconVec = Array<OneD, PreconditionerSharedPtr>(nvec);
            for(int i = 0; i < nvec; ++i)
            {
                m_preconVec[i] = CreatePrecon(m_locToGloMapVec[i]);
            }

            // Allocate memory for top-level structure
            SetupTopLevel(m_locToGloMapVec[0]);

            // Setup Block Matrix systems
            int n, n_exp = m_expListVec[0].lock()->GetNumElmts();

            MatrixStorage blkmatStorage = eDIAGONAL;
            const Array<OneD,const unsigned int>& nbdry_size
                    = m_locToGloMapVec[0]->GetNumLocalBndCoeffsPerPatch();

            m_S1Blk      = MemoryManager<DNekScalBlkMat>
                ::AllocateSharedPtr(nbdry_size, nbdry_size, blkmatStorage);

            // Preserve original matrix in m_S1Blk
            for (n = 0; n < n_exp; ++n)
            {
                DNekScalMatSharedPtr mat = m_schurCompl->GetBlock(n, n);
                m_S1Blk->SetBlock(n, n, mat);
            }

            // Build preconditioner
            for(int i =0; i < nvec; ++i)
            {
                m_preconVec[i]->BuildPreconditioner();
            }

            // Do transform of Schur complement matrix
            for (n = 0; n < n_exp; ++n)
            {
                if (m_linSysKey.GetMatrixType() !=
                        StdRegions::eHybridDGHelmBndLam)
                {
                    DNekScalMatSharedPtr mat = m_S1Blk->GetBlock(n, n);
                    DNekScalMatSharedPtr t = m_preconVec[0]->TransformedSchurCompl(n, mat);
                    m_schurCompl->SetBlock(n, n, t);
                }
            }

            // Construct this level
            Initialise(m_locToGloMapVec[0]);
        }
        
        /**
         *
         */
        GlobalLinSysIterativeStaticCondVec::~GlobalLinSysIterativeStaticCondVec()
        {
            
        }

        DNekScalBlkMatSharedPtr GlobalLinSysIterativeStaticCondVec::
            v_GetStaticCondBlock(unsigned int n)
        {
            DNekScalBlkMatSharedPtr schurComplBlock;
            int  scLevel           = m_locToGloMapVec[0]->GetStaticCondLevel();
            DNekScalBlkMatSharedPtr sc = scLevel == 0 ? m_S1Blk : m_schurCompl;
            DNekScalMatSharedPtr    localMat = sc->GetBlock(n,n);
            unsigned int nbdry    = localMat->GetRows();
            unsigned int nblks    = 1;
            unsigned int esize[1] = {nbdry};

            schurComplBlock = MemoryManager<DNekScalBlkMat>
                ::AllocateSharedPtr(nblks, nblks, esize, esize);
            schurComplBlock->SetBlock(0, 0, localMat);

            return schurComplBlock;
        }

        /**
         * Assemble the schur complement matrix from the block matrices stored
         * in #m_blkMatrices and the given local to global mapping information.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysIterativeStaticCondVec::v_AssembleSchurComplement(
            const AssemblyMapSharedPtr pLocToGloMap)
        {
            // Set up unique map
            v_UniqueMap();

            PrepareLocalSchurComplement();
            return;
        }


        /**
         * Populates sparse block-diagonal schur complement matrix from
         * the block matrices stored in #m_blkMatrices.
         */
        void GlobalLinSysIterativeStaticCondVec::PrepareLocalSchurComplement()
        {
            LocalMatrixStorageStrategy storageStrategy =
                m_expListVec[0].lock()->GetSession()->
                    GetSolverInfoAsEnum<LocalMatrixStorageStrategy>(
                                       "LocalMatrixStorageStrategy");

            switch(storageStrategy)
            {
                case MultiRegions::eContiguous:
                case MultiRegions::eNonContiguous:
                {
                    size_t storageSize = 0;
                    int nBlk           = m_schurCompl->GetNumberOfBlockRows();

                    m_scale = Array<OneD, NekDouble> (nBlk, 1.0);
                    m_rows  = Array<OneD, unsigned int> (nBlk, 0U);

                    // Determine storage requirements for dense blocks.
                    for (int i = 0; i < nBlk; ++i)
                    {
                        m_rows[i]    = m_schurCompl->GetBlock(i,i)->GetRows();
                        m_scale[i]   = m_schurCompl->GetBlock(i,i)->Scale();
                        storageSize += m_rows[i] * m_rows[i];
                    }

                    // Assemble dense storage blocks.
                    DNekScalMatSharedPtr loc_mat;
                    m_denseBlocks.resize(nBlk);
                    double *ptr = 0;

                    if (MultiRegions::eContiguous == storageStrategy)
                    {
                        m_storage.resize    (storageSize);
                        ptr = &m_storage[0];
                    }

                    for (unsigned int n = 0; n < nBlk; ++n)
                    {
                        loc_mat = m_schurCompl->GetBlock(n,n);

                        if (MultiRegions::eContiguous == storageStrategy)
                        {
                            int loc_lda      = loc_mat->GetRows();
                            int blockSize    = loc_lda * loc_lda;
                            m_denseBlocks[n] = ptr;
                            for(int i = 0; i < loc_lda; ++i)
                            {
                                for(int j = 0; j < loc_lda; ++j)
                                {
                                    ptr[j*loc_lda+i] = (*loc_mat)(i,j);
                                }
                            }
                            ptr += blockSize;
                            GlobalLinSys::v_DropStaticCondBlock(n);
                        }
                        else
                        {
                            m_denseBlocks[n] = loc_mat->GetRawPtr();
                        }
                    }
                    break;
                }
                case MultiRegions::eSparse:
                {
                    ASSERTL0(false,"Not set up for sparse storage");
                    break;
                }
                default:
                    ErrorUtil::NekError("Solver info property \
                        LocalMatrixStorageStrategy takes values \
                        Contiguous, Non-contiguous and Sparse");
            }
        }

        /**
         *
         */
        void GlobalLinSysIterativeStaticCondVec::v_DoMatrixMultiply(
                            const Array<OneD, int> nGlobal,
                            const Array<OneD, NekDouble>& pInput,
                            Array<OneD, NekDouble>& pOutput)
        {
            int nLocal = m_locToGloMapVec[0]->GetNumLocalBndCoeffs();
            bool doGlobalOp = m_expList.lock()->GetGlobalOptParam()->
                    DoGlobalMatOp(m_linSysKey.GetMatrixType());

            if(doGlobalOp)
            {
                ASSERTL0(false,"Option doGlobalOp not set up");
            }
            else if (m_sparseSchurCompl)
            {
                ASSERTL0(false,"Option spaceSchurCompl not set up");
            }
            else
            {
                int i,n, cnt;
                Array<OneD, NekDouble> tmp,tmp1,tmp2;

                int nvec = nGlobal.num_elements();
                
                RotPeriodicInfoSharedPtr perRotInfo = m_locToGloMapVec[0]->GetPerRotInfo();
                Array<OneD, int> periodicRotBndMap = m_locToGloMapVec[0]->GetPeriodicRotBndMap();

                // Do matrix multiply locally, using direct BLAS calls
                for(n = cnt = 0; n < nvec; cnt += nGlobal[n], ++n)
                {                    
                    m_locToGloMapVec[n]->GlobalToLocalBnd(pInput+cnt,
                                                          tmp = m_wsp+n*nLocal);
                }

                // put in bwd rotation term here.
                if(perRotInfo.get())
                {
                    if(nvec == 1)
                        perRotInfo->RotateBwd(periodicRotBndMap,m_wsp,
                                          tmp, tmp1);
                    else if(nvec == 2)
                        perRotInfo->RotateBwd(periodicRotBndMap,m_wsp,
                                          tmp  = m_wsp + nLocal,
                                          tmp1);
                    else
                        perRotInfo->RotateBwd(periodicRotBndMap,m_wsp,
                                          tmp  = m_wsp + nLocal,
                                          tmp1 = m_wsp+ 2*nLocal);
                }
            
                Array<OneD, NekDouble> tmpout = m_wsp + nLocal*nvec;
                
                cnt = 0; 
                for(n = 0; n < nvec; ++n)
                {
                    for (i = 0; i < m_denseBlocks.size(); cnt += m_rows[i], ++i)
                    {
                        const int rows = m_rows[i];
                        Blas::Dgemv('N', rows, rows,
                                    m_scale[i], m_denseBlocks[i], rows,
                                    m_wsp.get()+cnt, 1,
                                    0.0, tmpout.get()+cnt, 1);
                    }
                }
                
                // put in fwd rotation term here.
                if(perRotInfo.get())
                {
                    if(nvec == 1)
                        perRotInfo->RotateFwd(periodicRotBndMap,tmp = m_wsp + nLocal*nvec,
                                          tmp1, tmp2);
                    else if(nvec == 2)
                        perRotInfo->RotateFwd(periodicRotBndMap,tmp = m_wsp + nLocal*nvec,
                                          tmp1 = m_wsp + nLocal*(nvec+1), tmp2);
                    else
                        perRotInfo->RotateFwd(periodicRotBndMap,tmp = m_wsp + nLocal*nvec,
                                          tmp1 = m_wsp + nLocal*(nvec+1),
                                          tmp2 = m_wsp + nLocal*(nvec+2));

                }

                for(n = cnt = 0; n < nvec; cnt += nGlobal[n], ++n)
                {                    
                    m_locToGloMapVec[n]->AssembleBnd(tmpout+n*nLocal, tmp = pOutput+cnt);
                }
            }
        }

        void GlobalLinSysIterativeStaticCondVec::v_UniqueMap()
        {

            m_mapVec = Array<OneD, Array<OneD, int> > (m_locToGloMapVec.num_elements());
            
            for(int n = 0; n < m_locToGloMapVec.num_elements(); ++n)
            {
                m_mapVec[n] = m_locToGloMapVec[n]->GetGlobalToUniversalBndMapUnique();
            }
        }

        DNekScalBlkMatSharedPtr GlobalLinSysIterativeStaticCondVec::v_PreSolve(
            int                     scLevel,
            Array<OneD, Array<OneD, NekDouble> >  &F)
        {

            if (scLevel == 0)
            {
                // When matrices are supplied to the constructor at the top
                // level, the preconditioner is never set up.

                if (m_preconVec.num_elements() == 0)
                {
                    
                    int nvec = m_locToGloMapVec.num_elements();
                    m_preconVec = Array<OneD, PreconditionerSharedPtr>(nvec);
                    
                    for(int n = 0; n < nvec; ++n)
                    {
                        m_preconVec[n] = CreatePrecon(m_locToGloMapVec[n]);
                        m_preconVec[n]->BuildPreconditioner();
                    }
                }
                
                int nLocBndDofs   = m_locToGloMapVec[0]->GetNumLocalBndCoeffs();

                ASSERTL1(F[0].num_elements() >= nLocBndDofs,"Wrong size array");
                
                //Set_Rhs_Magnitude - version using local array

                Array<OneD, NekDouble> vExchange(1, 0.0);

                for(int i = 0;  i < F.num_elements(); ++i)
                {
                    vExchange[0] += Blas::Ddot(nLocBndDofs, F[i],1,F[i],1);
                }

                m_expListVec[0].lock()->GetComm()->GetRowComm()->AllReduce(
                vExchange, Nektar::LibUtilities::ReduceSum);

                // To ensure that very different rhs values are not being
                // used in subsequent solvers such as the velocit solve in
                // INC NS. If this works we then need to work out a better
                // way to control this.
                NekDouble new_rhs_mag = (vExchange[0] > 1e-6)? vExchange[0] : 1.0;
                
                if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
                {
                    m_rhs_magnitude = new_rhs_mag;
                }
                else
                {
                    m_rhs_magnitude = (m_rhs_mag_sm*(m_rhs_magnitude) + 
                                       (1.0-m_rhs_mag_sm)*new_rhs_mag); 
                }
 
                return m_S1Blk;
            }
            else
            {
                // for multilevel iterative solver always use rhs
                // vector value with no weighting
                m_rhs_magnitude = NekConstants::kNekUnsetDouble;

                return m_schurCompl;
            }
        }

        void GlobalLinSysIterativeStaticCondVec::v_BasisTransformLoc(
                       Array<OneD, NekDouble>& pInOut)
        {
            m_preconVec[0]->DoTransformToLowEnergyLoc(pInOut);
        }

        void GlobalLinSysIterativeStaticCondVec::v_BasisInvTransformLoc(
            Array<OneD, NekDouble>& pInOut)
        {
            m_preconVec[0]->DoTransformFromLowEnergyLoc(pInOut);
        }
    }
}
