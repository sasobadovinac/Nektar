///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysSaenaStaticCond.cpp
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
// Description: GlobalLinSysSaenaStaticCond definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysSaenaStaticCond.h>

//#include <petscsys.h>
//#include <petscksp.h>
//#include <petscmat.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysSaena
         *
         * Solves a linear system using single- or multi-level static
         * condensation.
         */

        /**
         * Registers the class with the Factory.
         */
        string GlobalLinSysSaenaStaticCond::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "SaenaStaticCond",
                    GlobalLinSysSaenaStaticCond::create,
                    "Saena static condensation.");

        string GlobalLinSysSaenaStaticCond::className2
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "SaenaMultiLevelStaticCond",
                    GlobalLinSysSaenaStaticCond::create,
                    "Saena multi-level static condensation.");

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
        GlobalLinSysSaenaStaticCond::GlobalLinSysSaenaStaticCond(
                     const GlobalLinSysKey                &pKey,
                     const std::weak_ptr<ExpList>         &pExpList,
                     const std::shared_ptr<AssemblyMap>   &pLocToGloMap)
            : GlobalLinSys          (pKey, pExpList, pLocToGloMap),
              GlobalLinSysSaena     (pKey, pExpList, pLocToGloMap),
              GlobalLinSysStaticCond(pKey, pExpList, pLocToGloMap)
        {
            std::cout << __func__ << std::endl;

            ASSERTL1((pKey.GetGlobalSysSolnType()==eSaenaStaticCond)||
                     (pKey.GetGlobalSysSolnType()==eSaenaMultiLevelStaticCond),
                     "This constructor is only valid when using static "
                     "condensation");
            ASSERTL1(pKey.GetGlobalSysSolnType()
                        == pLocToGloMap->GetGlobalSysSolnType(),
                     "The local to global map is not set up for the requested "
                     "solution type");
        }

        /**
         *
         */
        GlobalLinSysSaenaStaticCond::GlobalLinSysSaenaStaticCond(
                     const GlobalLinSysKey                &pKey,
                     const std::weak_ptr<ExpList>         &pExpList,
                     const DNekScalBlkMatSharedPtr         pSchurCompl,
                     const DNekScalBlkMatSharedPtr         pBinvD,
                     const DNekScalBlkMatSharedPtr         pC,
                     const DNekScalBlkMatSharedPtr         pInvD,
                     const std::shared_ptr<AssemblyMap>   &pLocToGloMap,
                     const PreconditionerSharedPtr         pPrecon)
            : GlobalLinSys          (pKey, pExpList, pLocToGloMap),
              GlobalLinSysSaena     (pKey, pExpList, pLocToGloMap),
              GlobalLinSysStaticCond(pKey, pExpList, pLocToGloMap)
        {
            std::cout << __func__ << std::endl;

            m_schurCompl = pSchurCompl;
            m_BinvD      = pBinvD;
            m_C          = pC;
            m_invD       = pInvD;
            m_precon     = pPrecon;
        }

        /**
         *
         */
        GlobalLinSysSaenaStaticCond::~GlobalLinSysSaenaStaticCond()
        {

        }

        /**
         * Assemble the schur complement matrix from the block matrices stored
         * in #m_blkMatrices and the given local to global mapping information.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysSaenaStaticCond::v_AssembleSchurComplement(
            AssemblyMapSharedPtr pLocToGloMap)
        {
            std::cout << __func__ << std::endl;

            int i, j, n, cnt, gid1, gid2, loc_lda;
            NekDouble sign1, sign2, value;

            const int nDirDofs = pLocToGloMap->GetNumGlobalDirBndCoeffs();

            DNekScalBlkMatSharedPtr SchurCompl = m_schurCompl;
            DNekScalBlkMatSharedPtr BinvD      = m_BinvD;
            DNekScalBlkMatSharedPtr C          = m_C;
            DNekScalBlkMatSharedPtr invD       = m_invD;
            DNekScalMatSharedPtr    loc_mat;

            // CALCULATE REORDERING MAPPING
            CalculateReordering(pLocToGloMap->GetGlobalToUniversalBndMap(),
                                pLocToGloMap->GetGlobalToUniversalBndMapUnique(),
                                pLocToGloMap);

            // SET UP VECTORS AND MATRIX
//            SetUpMatVec(pLocToGloMap->GetNumGlobalBndCoeffs(), nDirDofs);
            SetUpMatVec();

            // CONSTRUCT KSP OBJECT
            SetUpSolver(pLocToGloMap->GetIterativeTolerance());

            // POPULATE MATRIX
            for(n = cnt = 0; n < m_schurCompl->GetNumberOfBlockRows(); ++n)
            {
                loc_mat = m_schurCompl->GetBlock(n,n);
                loc_lda = loc_mat->GetRows();

                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalBndMap(cnt + i)-nDirDofs;
                    sign1 = pLocToGloMap->GetLocalToGlobalBndSign(cnt + i);
                    if(gid1 >= 0)
                    {
                        int gid1ro = m_reorderedMap[gid1];
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = pLocToGloMap->GetLocalToGlobalBndMap(cnt + j)
                                                                    - nDirDofs;
                            sign2 = pLocToGloMap->GetLocalToGlobalBndSign(cnt + j);
                            if(gid2 >= 0)
                            {
                                int gid2ro = m_reorderedMap[gid2];
                                value = sign1*sign2*(*loc_mat)(i,j);
                                m_matrix.set(gid1ro, gid2ro, value);
                            }
                        }
                    }
                }
                cnt += loc_lda;
            }

            m_matrix.assemble();
        }

        GlobalLinSysStaticCondSharedPtr GlobalLinSysSaenaStaticCond::v_Recurse(
            const GlobalLinSysKey                &mkey,
            const std::weak_ptr<ExpList>         &pExpList,
            const DNekScalBlkMatSharedPtr         pSchurCompl,
            const DNekScalBlkMatSharedPtr         pBinvD,
            const DNekScalBlkMatSharedPtr         pC,
            const DNekScalBlkMatSharedPtr         pInvD,
            const std::shared_ptr<AssemblyMap>   &l2gMap)
        {
//            GlobalLinSysSaenaStaticCondSharedPtr sys = MemoryManager<
//                GlobalLinSysSaenaStaticCond>::AllocateSharedPtr(
//                    mkey, pExpList, pSchurCompl, pBinvD, pC, pInvD, l2gMap,
//                    m_precon);
            GlobalLinSysSaenaStaticCondSharedPtr sys = MemoryManager<
                    GlobalLinSysSaenaStaticCond>::AllocateSharedPtr(
                    mkey, pExpList, pSchurCompl, pBinvD, pC, pInvD, l2gMap);

            std::cout << __func__ << std::endl;

            sys->Initialise(l2gMap);
            return sys;
        }
    }
}
