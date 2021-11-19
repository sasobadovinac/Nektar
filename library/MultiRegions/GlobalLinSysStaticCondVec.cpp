///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysStaticCondVec.cpp
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

#include <MultiRegions/GlobalLinSysStaticCondVec.h>
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
         * @class GlobalLinSysStaticCondVec
         *
         * Solves a linear system using single- or multi-level static
         * condensation.
         */

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
        GlobalLinSysStaticCondVec::GlobalLinSysStaticCondVec(
            const GlobalLinSysKey                &pKey,
            const Array<OneD, std::weak_ptr<ExpList> > &pVecExpList,
            const Array<OneD, std::shared_ptr<AssemblyMap> > &pVecLocToGloMap)
                : GlobalLinSys(pKey, pVecExpList, pVecLocToGloMap[0]),
                  m_locToGloMapVec (pVecLocToGloMap)
        {
        }

        void GlobalLinSysStaticCondVec::v_InitObject()
        {
            // Allocate memory for top-level structures of Static Condensed systems
            SetupTopLevel(m_locToGloMapVec[0]);

            // Construct this level
            Initialise(m_locToGloMapVec[0]);
        }
        
        /**
         *
         */
        GlobalLinSysStaticCondVec::~GlobalLinSysStaticCondVec()
        {

        }
        
        
        /**
         *
         */
        void GlobalLinSysStaticCondVec::v_SolveVec(
                      const Array<OneD, Array<OneD,  NekDouble>  >&in,
                      Array<OneD,  Array<OneD,       NekDouble> > &out)
        {
            int nvec = m_locToGloMapVec.num_elements();
            
            Array<OneD,int> nGlobBndDofs(nvec);
            Array<OneD,int> nDirBndDofs(nvec);

            // Can evaluate these sizes here since are assuming they
            // are constant over each vector components.
            int nGlobDofs   = m_locToGloMapVec[0]->GetNumGlobalCoeffs();
            int nLocBndDofs = m_locToGloMapVec[0]->GetNumLocalBndCoeffs();
            int nIntDofs    = nGlobDofs - m_locToGloMapVec[0]->GetNumGlobalBndCoeffs();
            
            Array<OneD, Array<OneD, NekDouble > > F(nvec), V_locbnd(nvec), F_bnd(nvec);
            Array<OneD, NekDouble > tmp;
                        
            int n;
            
            // Set up definitions and get local bnd vectoår 
            for(n = 0; n < nvec; ++n)
            {
                // This value is actually the same but needed for LinearSolve call 
                nGlobBndDofs[n]  = m_locToGloMapVec[n]->GetNumGlobalBndCoeffs();
                // This value can change if number of Dirichlet BCs
                // different in each component
                nDirBndDofs[n]   = m_locToGloMapVec[n]->GetNumGlobalDirBndCoeffs();

                // These first two arrays need to be untouched until end
                // First space is used in Multiply_a operation 
                V_locbnd[n] = m_wsp + 2*nvec*nLocBndDofs + n*nLocBndDofs;
                F[n]        = m_wsp + 3*nvec*nLocBndDofs + n*nGlobDofs;
                F_bnd[n]    = m_wsp + n*nLocBndDofs;

                m_locToGloMapVec[n]->LocalToLocalBnd(in[n], F_bnd[n]);
            }
            
            DNekScalBlkMatSharedPtr sc = v_PreSolve(0, F_bnd);
            
            for(n = 0; n < nvec; ++n)
            {
                // Gather boundary expansison into locbnd 
                m_locToGloMapVec[n]->LocalToLocalBnd(out[n],V_locbnd[n]);
                NekVector<NekDouble> V_LocBnd(nLocBndDofs,V_locbnd[n],eWrapper);
                NekVector<NekDouble> F_Bnd(nLocBndDofs,m_wsp +nvec*nLocBndDofs,eWrapper);

                // construct boundary forcing
                if(nIntDofs)
                {
                    m_locToGloMapVec[n]->LocalToLocalInt(in[n],
                                          tmp = F[n]+nGlobBndDofs[n]);
                    NekVector<NekDouble> F_Int(nIntDofs,
                                          tmp = F[n]+nGlobBndDofs[n],eWrapper);

                    DNekScalBlkMat &BinvD      = *m_BinvD;
                    DNekScalBlkMat &SchurCompl = *sc;

                    // include dirichlet boundary forcing
                    F_Bnd = BinvD*F_Int + SchurCompl*V_LocBnd;
                }
                else
                {
                    // include dirichlet boundary forcing
                    DNekScalBlkMat &SchurCompl = *sc;
                    F_Bnd = SchurCompl*V_LocBnd;
                }

                Vmath::Vsub(nLocBndDofs, &F_bnd[n][0],1, &F_Bnd[0], 1,
                            &F_bnd[n][0],1);

                v_BasisTransformLoc(F_bnd[n]);
            }


            RotPeriodicInfoSharedPtr perRotInfo = m_locToGloMapVec[0]->GetPerRotInfo();
            Array<OneD, int> periodicRotBndMap = m_locToGloMapVec[0]->GetPeriodicRotBndMap();
            // for(auto &it: periodicRotBndMap)
            // {
            //     cout << "\t\t\t id=" << it << endl;
            //     cout << m_locToGloMapVec[0]->GetLocalToGlobalMap(it) << endl;
            //     // cout << m_locToGloMapVec[0]->GetGlobalToUniversalMap(it) << endl;
            //     // cout << m_locToGloMapVec[0]->GetGlobalToUniversalMapUnique(it) << endl;
            //     // cout << m_locToGloMapVec[0]->GetLocalToGlobalBndMap(it) << endl;
            //     // cout << m_locToGloMapVec[0]->GetBndCondCoeffsToGlobalCoeffsMap(it) << endl;
            // }
            
            // put in fwd rotation term here.
            if(perRotInfo.get())
            {
                // cout << "---------------------------------------VCS.cpp RotateFwd "  << F_bnd[0].num_elements() <<
                // " " << periodicRotBndMap.num_elements() << " -------------------------------" << endl;
                // for(auto &it : periodicRotBndMap) cout << it << endl;
                
                // for(auto &it : F_bnd[0]) cout << it << ", ";
                // cout << endl;
                // for(auto &it : F_bnd[1]) cout << it << ", ";
                // cout << endl;

                Array<OneD, NekDouble> tmp;
                if(nvec == 1)
                    perRotInfo->RotateFwd(periodicRotBndMap,F_bnd[0],tmp, tmp);
                else if(nvec == 2)
                    perRotInfo->RotateFwd(periodicRotBndMap,F_bnd[0],F_bnd[1],
                                        tmp);
                else
                    perRotInfo->RotateFwd(periodicRotBndMap,F_bnd[0],F_bnd[1],
                                        F_bnd[2]);

                // for(auto &it : F_bnd[0]) cout << it << ", ";
                // cout << endl;
                // for(auto &it : F_bnd[1]) cout << it << ", ";
                // cout << endl;
            }

            // calculate globally  condensed forcing
            for(n = 0; n < nvec; ++n)
            {
                m_locToGloMapVec[n]->AssembleBnd(F_bnd[n], F[n]);
            }

            // solve boundary system
            // Solve for difference from initial solution given inout;
            SolveVecLinearSystem(nGlobBndDofs, F, out, m_locToGloMapVec,
                                 nDirBndDofs);
                
            
            Array<OneD, Array<OneD, NekDouble> > outloc(3); 

            for(n = 0; n < nvec; ++n)
            {
                outloc[n] = m_wsp + n*nLocBndDofs;
                
                // put solution into local format
                m_locToGloMapVec[n]->GlobalToLocalBnd(out[n],outloc[n]);
            }

            // periodic rotate bwd
            if(perRotInfo.get())
            {
                // cout << "---------------------------------------VCS.cpp RotateBwd "  << outloc[0].num_elements() <<
                // " " << periodicRotBndMap.num_elements() << " -------------------------------" << endl;
                // for(auto &it : periodicRotBndMap) cout << it << endl;
                
                // for(auto &it : outloc[0]) cout << it << ", ";
                // cout << endl;
                // for(auto &it : outloc[1]) cout << it << ", ";
                // cout << endl;

                Array<OneD, NekDouble> tmp;
                if(nvec == 1)
                    perRotInfo->RotateBwd(periodicRotBndMap,outloc[0],tmp, tmp);
                else if(nvec == 2)
                    perRotInfo->RotateBwd(periodicRotBndMap,outloc[0],outloc[1],
                                      tmp);
                else
                    perRotInfo->RotateBwd(periodicRotBndMap,outloc[0],outloc[1],
                                      outloc[2]);

                // for(auto &it : outloc[0]) cout << it << ", ";
                // cout << endl;
                // for(auto &it : outloc[1]) cout << it << ", ";
                // cout << endl;
            }
            
            Array<OneD, NekDouble> V_int = m_wsp + 3*nvec*nLocBndDofs + nvec*nGlobDofs;
            for(n = 0; n < nvec; ++n)
            {
                // Transform back to original basis 
                v_BasisInvTransformLoc(outloc[n]);
                
                Vmath::Vadd(nLocBndDofs, V_locbnd[n], 1, outloc[n], 1, V_locbnd[n],1);

                // put final solution back in out array
                m_locToGloMapVec[n]->LocalBndToLocal(V_locbnd[n],out[n]);

                // solve interior system
                if(nIntDofs)
                {
                    // get array of local solutions
                    DNekScalBlkMat &invD  = *m_invD;
                    
                    NekVector<NekDouble> V_Int(nIntDofs, V_int ,eWrapper);
                    NekVector<NekDouble> F_Int(nIntDofs, tmp = F[n]+nGlobBndDofs[n],
                                               eWrapper);

                    if(nGlobBndDofs[n] - nDirBndDofs[n] || nDirBndDofs[n])
                    {
                        NekVector<NekDouble> V_LocBnd(nLocBndDofs, V_locbnd[n],eWrapper);
                        
                        DNekScalBlkMat &C     = *m_C;
                        
                        F_Int = F_Int - C*V_LocBnd;
                    }
                    
                    Multiply(V_Int, invD, F_Int);

                    m_locToGloMapVec[n]->LocalIntToLocal(V_int,out[n]);
                }

                
            }
        }


        /**
         * If at the last level of recursion (or the only level in the case of
         * single-level static condensation), assemble the Schur complement.
         * For other levels, in the case of multi-level static condensation,
         * the next level of the condensed system is computed.
         * @param   pLocToGloMap    Local to global mapping.
         */
        void GlobalLinSysStaticCondVec::v_Initialise(
                                     const std::shared_ptr<AssemblyMap>& pLocToGloMap)
        {
            int nvec = m_locToGloMapVec.num_elements();

            int nLocalBnd  = m_locToGloMapVec[0]->GetNumLocalBndCoeffs();
            int nGlobal = m_locToGloMapVec[0]->GetNumGlobalCoeffs();
            
            m_wsp = Array<OneD, NekDouble>(nvec*3*nLocalBnd + (nvec+1)*nGlobal, 0.0);

            v_AssembleSchurComplement(m_locToGloMapVec[0]);
        }

        int GlobalLinSysStaticCondVec::v_GetNumBlocks()
        {
            return m_schurCompl->GetNumberOfBlockRows();
        }

        /**
         * For the first level in multi-level static condensation, or the only
         * level in the case of single-level static condensation, allocate the
         * condensed matrices and populate them with the local matrices
         * retrieved from the expansion list.
         * @param
         */
        void GlobalLinSysStaticCondVec::SetupTopLevel(
                const std::shared_ptr<AssemblyMap>& pLocToGloMap)
        {
            int n;
            int n_exp = m_expList.lock()->GetNumElmts();

            const Array<OneD,const unsigned int>& nbdry_size
                    = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
            const Array<OneD,const unsigned int>& nint_size
                    = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();

            // Setup Block Matrix systems
            MatrixStorage blkmatStorage = eDIAGONAL;
            m_schurCompl = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nbdry_size, blkmatStorage);
            m_BinvD      = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nint_size , blkmatStorage);
            m_C          = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nbdry_size, blkmatStorage);
            m_invD       = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nint_size , blkmatStorage);

            for(n = 0; n < n_exp; ++n)
            {
                if (m_linSysKey.GetMatrixType() ==
                        StdRegions::eHybridDGHelmBndLam)
                {
                    DNekScalMatSharedPtr loc_mat
                        = GlobalLinSys::v_GetBlock(n);
                    m_schurCompl->SetBlock(n,n,loc_mat);
                }
                else
                {
                    DNekScalBlkMatSharedPtr loc_schur
                        = GlobalLinSys::v_GetStaticCondBlock(n);
                    DNekScalMatSharedPtr t;
                    m_schurCompl->SetBlock(n, n, t = loc_schur->GetBlock(0,0));
                    m_BinvD     ->SetBlock(n, n, t = loc_schur->GetBlock(0,1));
                    m_C         ->SetBlock(n, n, t = loc_schur->GetBlock(1,0));
                    m_invD      ->SetBlock(n, n, t = loc_schur->GetBlock(1,1));
                }
            }
        }

    }
}
