///////////////////////////////////////////////////////////////////////////////
//
// File StdNodalTriExp.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Header field for Nodal triangle routines built upon
// StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDNODALTRIEXP_H
#define STDNODALTRIEXP_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar {
    namespace StdRegions {
        class StdNodalTriExp : virtual public StdTriExp {
        public:
            STD_REGIONS_EXPORT StdNodalTriExp() = default;

            STD_REGIONS_EXPORT StdNodalTriExp(
                    const LibUtilities::BasisKey &Ba,
                    const LibUtilities::BasisKey &Bb,
                    LibUtilities::PointsType Ntype);

            STD_REGIONS_EXPORT StdNodalTriExp(const StdNodalTriExp &T);

            STD_REGIONS_EXPORT ~StdNodalTriExp() override = default;

            //-------------------------------
            // Nodal basis specific routines
            //-------------------------------
            STD_REGIONS_EXPORT void NodalToModal(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT void NodalToModalTranspose(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT void ModalToNodal(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT void GetNodalPoints(
                    Array<OneD, const NekDouble> &x,
                    Array<OneD, const NekDouble> &y);

            STD_REGIONS_EXPORT DNekMatSharedPtr GenNBasisTransMatrix();


        protected:
            LibUtilities::PointsKey m_nodalPointsKey;

            STD_REGIONS_EXPORT const LibUtilities::PointsKey
            v_GetNodalPointsKey() const override
            {
                return m_nodalPointsKey;
            };

            STD_REGIONS_EXPORT bool v_IsNodalNonTensorialExp() override;

            //---------------------------------------
            // Transforms
            //---------------------------------------
            STD_REGIONS_EXPORT void v_BwdTrans(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray) override;

            STD_REGIONS_EXPORT void v_BwdTrans_SumFac(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray) override;

            STD_REGIONS_EXPORT void v_FwdTrans(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray) override;


            //---------------------------------------
            // Inner product functions
            //---------------------------------------
            STD_REGIONS_EXPORT void v_IProductWRTBase(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray) override;

            STD_REGIONS_EXPORT void v_IProductWRTBase_SumFac(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray,
                    bool multiplybyweights = true) override;

            STD_REGIONS_EXPORT void v_IProductWRTDerivBase(
                    int dir,
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray) override;

            STD_REGIONS_EXPORT void v_IProductWRTDerivBase_SumFac(
                    int dir,
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray) override;


            //---------------------------------------
            // Evaluation functions
            //---------------------------------------
            STD_REGIONS_EXPORT void v_FillMode(
                    int mode,
                    Array<OneD, NekDouble> &outarray) override;


            //---------------------------
            // Helper functions
            //---------------------------
            STD_REGIONS_EXPORT int v_NumBndryCoeffs() const override;

            //--------------------------
            // Mappings
            //--------------------------
            STD_REGIONS_EXPORT int v_GetVertexMap(int localVertexId,
                                                  bool useCoeffPacking = false) override;

            STD_REGIONS_EXPORT void v_GetTraceToElementMap(
                    int eid,
                    Array<OneD, unsigned int> &maparray,
                    Array<OneD, int> &signarray,
                    Orientation edgeOrient = eForwards,
                    int P = -1,
                    int Q = -1) override;

            STD_REGIONS_EXPORT void v_GetTraceInteriorToElementMap(
                    int eid,
                    Array<OneD, unsigned int> &maparray,
                    Array<OneD, int> &signarray,
                    Orientation edgeOrient = eForwards) override;

            STD_REGIONS_EXPORT void v_GetInteriorMap(
                    Array<OneD, unsigned int> &outarray) override;

            STD_REGIONS_EXPORT void v_GetBoundaryMap(
                    Array<OneD, unsigned int> &outarray) override;


            //---------------------------------------
            // Operator evaluation functions
            //---------------------------------------
            STD_REGIONS_EXPORT void v_MassMatrixOp(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray,
                    const StdMatrixKey &mkey) override;

            STD_REGIONS_EXPORT void v_LaplacianMatrixOp(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray,
                    const StdMatrixKey &mkey) override;

            STD_REGIONS_EXPORT void v_LaplacianMatrixOp(
                    int k1,
                    int k2,
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray,
                    const StdMatrixKey &mkey) override;

            STD_REGIONS_EXPORT void v_WeakDerivMatrixOp(
                    int i,
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray,
                    const StdMatrixKey &mkey) override;

            STD_REGIONS_EXPORT void v_HelmholtzMatrixOp(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray,
                    const StdMatrixKey &mkey) override;


            //---------------------------------------
            // Wrapper functions
            //---------------------------------------
            STD_REGIONS_EXPORT DNekMatSharedPtr v_GenMatrix(
                    const StdMatrixKey &mkey) override;

            STD_REGIONS_EXPORT DNekMatSharedPtr v_CreateStdMatrix(
                    const StdMatrixKey &mkey) override;
        };

        typedef std::shared_ptr<StdNodalTriExp> StdNodalTriExpSharedPtr;
    } //end of namespace
} //end of namespace

#endif //STDNODALTRIEXP_H
