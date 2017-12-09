///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysStaticCondVec.h
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
// Description: A collection of routines common to statically condensed systems.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSSTATICCONDVEC_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSSTATICCONDVEC_H

#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/GlobalLinSysIterativeVec.h>
#include <LibUtilities/LinearAlgebra/SparseMatrixFwd.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;
        class GlobalLinSysStaticCondVec;

        typedef std::shared_ptr<GlobalLinSysStaticCondVec>
            GlobalLinSysStaticCondVecSharedPtr;

        /// A global linear system.
        class GlobalLinSysStaticCondVec : virtual public GlobalLinSys
        {
        public:
            /// Constructor for full direct matrix solve.
            GlobalLinSysStaticCondVec(
                const GlobalLinSysKey                &mkey,
                const Array<OneD, std::weak_ptr<ExpList> > &pVecExpList,
                const Array<OneD, std::shared_ptr<AssemblyMap> > &pVecLocToGloMap);

            virtual ~GlobalLinSysStaticCondVec();

        protected:

            virtual DNekScalBlkMatSharedPtr v_PreSolve(
                int                   scLevel,
                Array<OneD, Array<OneD, NekDouble> >  &F)
            {
                return m_schurCompl;
            }

            virtual void v_BasisTransformLoc(
                Array<OneD, NekDouble>& pInOut)
            {

            }

            virtual void v_BasisInvTransformLoc(
                Array<OneD, NekDouble>& pInOut)
            {
                
            }

            virtual void v_AssembleSchurComplement(
                std::shared_ptr<AssemblyMap> pLoctoGloMap)
            {
                
            }

            virtual int v_GetNumBlocks();

            /// Block Schur complement matrix.
            DNekScalBlkMatSharedPtr                  m_schurCompl;
            /// Block \f$ BD^{-1} \f$ matrix.
            DNekScalBlkMatSharedPtr                  m_BinvD;
            /// Block \f$ C \f$ matrix.
            DNekScalBlkMatSharedPtr                  m_C;
            /// Block \f$ D^{-1} \f$ matrix.
            DNekScalBlkMatSharedPtr                  m_invD;
            /// Local to global map.
            Array<OneD, std::shared_ptr<AssemblyMap> > m_locToGloMapVec;
            /// Workspace array for matrix multiplication
            Array<OneD, NekDouble>                   m_wsp;

            /// Rotational information for periodic non-planar cases
            //shared_ptr<RotPeriodicInfo>  m_perRotInfo;
            /// Rotational local coefficients for periodic non-planar cases
            //Array<OneD, int>           m_periodicRotMap;
            
            /// using a specified local to global map.
            virtual void v_SolveVec(
                      const Array<OneD, Array<OneD,  NekDouble> >&in,
                      Array<OneD,  Array<OneD,       NekDouble> > &out);

            virtual void v_InitObject();

            /// Initialise this object
            virtual void v_Initialise(
                    const std::shared_ptr<AssemblyMap>& locToGloMap);

            /// Set up the storage for the Schur complement or the top level
            /// of the multi-level Schur complement.
            void SetupTopLevel(
                    const std::shared_ptr<AssemblyMap>& locToGloMap);

            /// Set up rotational information for local coefficients
            /// if these are defined in mesh
            void SetupPeriodicRotation(void);
        };
    }
}

#endif
