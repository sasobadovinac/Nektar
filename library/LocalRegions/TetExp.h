///////////////////////////////////////////////////////////////////////////////
//
// File TetExp.h
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef TETEXP_H
#define TETEXP_H

#include <StdRegions/StdTetExp.h>
#include <SpatialDomains/TetGeom.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/LocalRegionsDeclspec.h>

namespace Nektar {
    namespace LocalRegions {

        class TetExp
                : virtual public StdRegions::StdTetExp,
                  virtual public Expansion3D {

        public:
            LOCAL_REGIONS_EXPORT TetExp(const LibUtilities::BasisKey &Ba,
                                        const LibUtilities::BasisKey &Bb,
                                        const LibUtilities::BasisKey &Bc,
                                        const SpatialDomains::TetGeomSharedPtr &geom);

            LOCAL_REGIONS_EXPORT TetExp(const TetExp &T);

            LOCAL_REGIONS_EXPORT ~TetExp() override = default;

        protected:
            //-----------------------------
            // Integration Methods
            //-----------------------------
            LOCAL_REGIONS_EXPORT NekDouble
            v_Integral(const Array<OneD, const NekDouble> &inarray) override;


            //-----------------------------
            // Differentiation Methods
            //-----------------------------
            LOCAL_REGIONS_EXPORT void v_PhysDeriv(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &out_d0,
                    Array<OneD, NekDouble> &out_d1,
                    Array<OneD, NekDouble> &out_d2) override;

            //-----------------------------
            // Transforms
            //-----------------------------
            LOCAL_REGIONS_EXPORT void v_FwdTrans(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray) override;

            //-----------------------------
            // Inner product functions
            //-----------------------------
            LOCAL_REGIONS_EXPORT void v_IProductWRTBase(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray) override;

            LOCAL_REGIONS_EXPORT void v_IProductWRTBase_SumFac(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray,
                    bool multiplybyweights = true) override;

            LOCAL_REGIONS_EXPORT void v_IProductWRTDerivBase(
                    int dir,
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray) override;

            LOCAL_REGIONS_EXPORT void v_AlignVectorToCollapsedDir(
                    int dir,
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray) override;


            //-----------------------------
            // Evaluation functions
            //-----------------------------
            LOCAL_REGIONS_EXPORT NekDouble v_StdPhysEvaluate(
                    const Array<OneD, const NekDouble> &Lcoord,
                    const Array<OneD, const NekDouble> &physvals) override;

            LOCAL_REGIONS_EXPORT NekDouble v_PhysEvaluate(
                    const Array<OneD, const NekDouble> &coords,
                    const Array<OneD, const NekDouble> &physvals) override;

            LOCAL_REGIONS_EXPORT NekDouble v_PhysEvaluate(
                    const Array<OneD, NekDouble> &coord,
                    const Array<OneD, const NekDouble> &inarray,
                    std::array<NekDouble, 3> &firstOrderDerivs) final;

            LOCAL_REGIONS_EXPORT void v_GetCoord(
                    const Array<OneD, const NekDouble> &Lcoords,
                    Array<OneD, NekDouble> &coords) override;

            LOCAL_REGIONS_EXPORT void v_GetCoords(
                    Array<OneD, NekDouble> &coords_1,
                    Array<OneD, NekDouble> &coords_2,
                    Array<OneD, NekDouble> &coords_3) override;

            //-----------------------------
            // Helper functions
            //-----------------------------
            LOCAL_REGIONS_EXPORT
            LibUtilities::ShapeType v_DetShapeType() const override;

            LOCAL_REGIONS_EXPORT
            StdRegions::StdExpansionSharedPtr v_GetStdExp() const override;

            LOCAL_REGIONS_EXPORT
            StdRegions::StdExpansionSharedPtr
            v_GetLinStdExp() const override;

            LOCAL_REGIONS_EXPORT int v_GetCoordim() override;

            LOCAL_REGIONS_EXPORT void v_ExtractDataToCoeffs(
                    const NekDouble *data,
                    const std::vector<unsigned int> &nummodes,
                    int mode_offset,
                    NekDouble *coeffs,
                    std::vector<LibUtilities::BasisType> &fromType) override;


            LOCAL_REGIONS_EXPORT void v_GetTracePhysMap(
                    int face,
                    Array<OneD, int> &outarray) override;

            LOCAL_REGIONS_EXPORT void
            v_ComputeTraceNormal(int face) override;
            //-----------------------------
            // Operator creation functions
            //-----------------------------
            LOCAL_REGIONS_EXPORT void v_HelmholtzMatrixOp(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray,
                    const StdRegions::StdMatrixKey &mkey) override;

            LOCAL_REGIONS_EXPORT void v_LaplacianMatrixOp(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray,
                    const StdRegions::StdMatrixKey &mkey) override;

            LOCAL_REGIONS_EXPORT void v_LaplacianMatrixOp(
                    int k1,
                    int k2,
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray,
                    const StdRegions::StdMatrixKey &mkey) override;

            LOCAL_REGIONS_EXPORT void v_SVVLaplacianFilter(
                    Array<OneD, NekDouble> &array,
                    const StdRegions::StdMatrixKey &mkey) override;

            //-----------------------------
            // Matrix creation functions
            //-----------------------------
            LOCAL_REGIONS_EXPORT DNekMatSharedPtr v_GenMatrix(
                    const StdRegions::StdMatrixKey &mkey) override;

            LOCAL_REGIONS_EXPORT DNekMatSharedPtr v_CreateStdMatrix(
                    const StdRegions::StdMatrixKey &mkey) override;

            LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr v_GetLocMatrix(
                    const MatrixKey &mkey) override;

            LOCAL_REGIONS_EXPORT
            DNekScalBlkMatSharedPtr v_GetLocStaticCondMatrix(
                    const MatrixKey &mkey) override;

            LOCAL_REGIONS_EXPORT void v_DropLocStaticCondMatrix(
                    const MatrixKey &mkey) override;

            LOCAL_REGIONS_EXPORT void SetUpInverseTransformationMatrix(
                    const DNekMatSharedPtr &m_transformationmatrix,
                    DNekMatSharedPtr m_inversetransformationmatrix,
                    DNekMatSharedPtr m_inversetransposedtransformationmatrix);

            LOCAL_REGIONS_EXPORT void v_ComputeConditionNumberOfMatrix(
                    const DNekScalMatSharedPtr &mat);

            LOCAL_REGIONS_EXPORT void v_ComputeLaplacianMetric() override;

            LOCAL_REGIONS_EXPORT void v_NormalTraceDerivFactors
                    (Array<OneD, Array<OneD, NekDouble> > &d0factors,
                     Array<OneD, Array<OneD, NekDouble> > &d1factors,
                     Array<OneD, Array<OneD, NekDouble> > &d2factors) override;

        private:
            LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
            LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;

            LOCAL_REGIONS_EXPORT void GeneralMatrixOp_MatOp(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray,
                    const StdRegions::StdMatrixKey &mkey);

            void v_LaplacianMatrixOp_MatFree_Kernel(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray,
                    Array<OneD, NekDouble> &wsp) override;
        };

        typedef std::shared_ptr<TetExp> TetExpSharedPtr;
        typedef std::vector<TetExpSharedPtr> TetExpVector;
    } //end of namespace
} //end of namespace

#endif // TETEXP_H
