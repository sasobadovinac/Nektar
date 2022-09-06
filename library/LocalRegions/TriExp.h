///////////////////////////////////////////////////////////////////////////////
//
// File TriExp.h
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
// Description: Expansion for triangular elements.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef TRIEXP_H
#define TRIEXP_H

#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
#include <LocalRegions/MatrixKey.h>
#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/TriGeom.h>
#include <StdRegions/StdTriExp.h>

namespace Nektar
{
namespace LocalRegions
{

class TriExp : virtual public StdRegions::StdTriExp, virtual public Expansion2D
{

public:
    /**
     * @brief Constructor using BasisKey class for quadrature
     * points and order definition
     */
    LOCAL_REGIONS_EXPORT TriExp(
        const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb,
        const SpatialDomains::Geometry2DSharedPtr &geom);

    LOCAL_REGIONS_EXPORT TriExp(const TriExp &T);

    LOCAL_REGIONS_EXPORT ~TriExp() override = default;

protected:
    //-------------------------------
    // Integration Methods
    //-------------------------------
    LOCAL_REGIONS_EXPORT NekDouble
    v_Integral(const Array<OneD, const NekDouble> &inarray) override;

    //----------------------------
    // Differentiation Methods
    //----------------------------
    LOCAL_REGIONS_EXPORT void v_PhysDeriv(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &out_d0, Array<OneD, NekDouble> &out_d1,
        Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray) override;
    LOCAL_REGIONS_EXPORT void v_PhysDeriv(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;
    LOCAL_REGIONS_EXPORT void v_PhysDirectionalDeriv(
        const Array<OneD, const NekDouble> &inarray,
        const Array<OneD, const NekDouble> &direction,
        Array<OneD, NekDouble> &out) override;

    //---------------------------------------
    // Transforms
    //---------------------------------------
    LOCAL_REGIONS_EXPORT void v_FwdTrans(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;
    LOCAL_REGIONS_EXPORT void v_FwdTransBndConstrained(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    //---------------------------------------
    // Inner product functions
    //---------------------------------------
    LOCAL_REGIONS_EXPORT void v_IProductWRTBase(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;
    LOCAL_REGIONS_EXPORT void v_IProductWRTDerivBase(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;
    LOCAL_REGIONS_EXPORT void v_IProductWRTBase_SumFac(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        bool multiplybyweights = true) override;
    LOCAL_REGIONS_EXPORT void v_IProductWRTBase_MatOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;
    LOCAL_REGIONS_EXPORT void v_IProductWRTDerivBase_SumFac(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;
    LOCAL_REGIONS_EXPORT void v_AlignVectorToCollapsedDir(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray) override;

    LOCAL_REGIONS_EXPORT void v_IProductWRTDerivBase_MatOp(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    LOCAL_REGIONS_EXPORT void v_IProductWRTDirectionalDerivBase(
        const Array<OneD, const NekDouble> &direction,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    LOCAL_REGIONS_EXPORT void v_IProductWRTDirectionalDerivBase_SumFac(
        const Array<OneD, const NekDouble> &direction,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    LOCAL_REGIONS_EXPORT void v_NormVectorIProductWRTBase(
        const Array<OneD, const NekDouble> &Fx,
        const Array<OneD, const NekDouble> &Fy,
        const Array<OneD, const NekDouble> &Fz,
        Array<OneD, NekDouble> &outarray) override;

    LOCAL_REGIONS_EXPORT void v_NormVectorIProductWRTBase(
        const Array<OneD, const Array<OneD, NekDouble>> &Fvec,
        Array<OneD, NekDouble> &outarray) override;

    //---------------------------------------
    // Evaluation functions
    //---------------------------------------
    LOCAL_REGIONS_EXPORT
    StdRegions::StdExpansionSharedPtr v_GetStdExp(void) const override;

    LOCAL_REGIONS_EXPORT
    StdRegions::StdExpansionSharedPtr v_GetLinStdExp(void) const override;

    LOCAL_REGIONS_EXPORT void v_GetCoord(
        const Array<OneD, const NekDouble> &Lcoords,
        Array<OneD, NekDouble> &coords) override;
    LOCAL_REGIONS_EXPORT void v_GetCoords(
        Array<OneD, NekDouble> &coords_1, Array<OneD, NekDouble> &coords_2,
        Array<OneD, NekDouble> &coords_3) override;
    LOCAL_REGIONS_EXPORT NekDouble
    v_StdPhysEvaluate(const Array<OneD, const NekDouble> &Lcoord,
                      const Array<OneD, const NekDouble> &physvals) override;
    LOCAL_REGIONS_EXPORT NekDouble
    v_PhysEvaluate(const Array<OneD, const NekDouble> &coord,
                   const Array<OneD, const NekDouble> &physvals) override;
    LOCAL_REGIONS_EXPORT NekDouble
    v_PhysEvaluate(const Array<OneD, NekDouble> &coord,
                   const Array<OneD, const NekDouble> &inarray,
                   std::array<NekDouble, 3> &firstOrderDerivs) final;
    LOCAL_REGIONS_EXPORT void v_GetTracePhysVals(
        const int edge, const StdRegions::StdExpansionSharedPtr &EdgeExp,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        StdRegions::Orientation orient) override;
    LOCAL_REGIONS_EXPORT void v_GetEdgeInterpVals(
        const int edge, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray);
    LOCAL_REGIONS_EXPORT void v_GetTraceQFactors(
        const int edge, Array<OneD, NekDouble> &outarray) override;
    LOCAL_REGIONS_EXPORT void v_ComputeTraceNormal(const int edge) override;

    //---------------------------------------
    // Helper functions
    //---------------------------------------
    LOCAL_REGIONS_EXPORT int v_GetCoordim() override;
    LOCAL_REGIONS_EXPORT void v_ExtractDataToCoeffs(
        const NekDouble *data, const std::vector<unsigned int> &nummodes,
        const int mode_offset, NekDouble *coeffs,
        std::vector<LibUtilities::BasisType> &fromType) override;
    LOCAL_REGIONS_EXPORT
    StdRegions::Orientation v_GetTraceOrient(int edge) override;
    LOCAL_REGIONS_EXPORT virtual const LibUtilities::BasisSharedPtr &v_GetBasis(
        int dir) const;
    LOCAL_REGIONS_EXPORT virtual int v_GetNumPoints(const int dir) const;
    LOCAL_REGIONS_EXPORT void v_GetTracePhysMap(
        const int edge, Array<OneD, int> &outarray) override;

    //---------------------------------------
    // Matrix creation functions
    //---------------------------------------
    LOCAL_REGIONS_EXPORT
    DNekMatSharedPtr v_GenMatrix(const StdRegions::StdMatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT
    DNekMatSharedPtr v_CreateStdMatrix(
        const StdRegions::StdMatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT
    DNekScalMatSharedPtr v_GetLocMatrix(const MatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT
    DNekScalBlkMatSharedPtr v_GetLocStaticCondMatrix(
        const MatrixKey &mkey) override;

    LOCAL_REGIONS_EXPORT void v_DropLocStaticCondMatrix(
        const MatrixKey &mkey) override;

    LOCAL_REGIONS_EXPORT void v_MassMatrixOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::StdMatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT void v_LaplacianMatrixOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::StdMatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT void v_LaplacianMatrixOp(
        const int k1, const int k2, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::StdMatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT void v_WeakDerivMatrixOp(
        const int i, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::StdMatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT void v_WeakDirectionalDerivMatrixOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::StdMatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT void v_MassLevelCurvatureMatrixOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::StdMatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT void v_HelmholtzMatrixOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::StdMatrixKey &mkey) override;
    LOCAL_REGIONS_EXPORT void v_GeneralMatrixOp_MatOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::StdMatrixKey &mkey) final;
    LOCAL_REGIONS_EXPORT void v_LaplacianMatrixOp_MatFree_Kernel(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, Array<OneD, NekDouble> &wsp) override;
    LOCAL_REGIONS_EXPORT void v_ReduceOrderCoeffs(
        int numMin, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    LOCAL_REGIONS_EXPORT void v_ComputeLaplacianMetric() override;

    LOCAL_REGIONS_EXPORT void v_SVVLaplacianFilter(
        Array<OneD, NekDouble> &array,
        const StdRegions::StdMatrixKey &mkey) override;

    LOCAL_REGIONS_EXPORT void v_NormalTraceDerivFactors(
        Array<OneD, Array<OneD, NekDouble>> &factors,
        Array<OneD, Array<OneD, NekDouble>> &d0factors,
        Array<OneD, Array<OneD, NekDouble>> &d1factors) override;

private:
    LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess>
        m_matrixManager;
    LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess>
        m_staticCondMatrixManager;
};

typedef std::shared_ptr<TriExp> TriExpSharedPtr;
typedef std::vector<TriExpSharedPtr> TriExpVector;
} // namespace LocalRegions
} // namespace Nektar

#endif
