///////////////////////////////////////////////////////////////////////////////
//
// File PyrExp.h
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
// Description: Header file for PyrExp routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef PYREXP_H
#define PYREXP_H

#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
#include <LocalRegions/MatrixKey.h>
#include <SpatialDomains/PyrGeom.h>
#include <StdRegions/StdPyrExp.h>

namespace Nektar
{
namespace LocalRegions
{

class PyrExp final : virtual public StdRegions::StdPyrExp,
                     virtual public Expansion3D
{
public:
    /** \brief Constructor using BasisKey class for quadrature points and order
        definition */
    LOCAL_REGIONS_EXPORT PyrExp(const LibUtilities::BasisKey &Ba,
                                const LibUtilities::BasisKey &Bb,
                                const LibUtilities::BasisKey &Bc,
                                const SpatialDomains::PyrGeomSharedPtr &geom);

    LOCAL_REGIONS_EXPORT PyrExp(const PyrExp &T);

    LOCAL_REGIONS_EXPORT ~PyrExp() final = default;

protected:
    //-------------------------------
    // Integration Methods
    //-------------------------------
    LOCAL_REGIONS_EXPORT NekDouble
    v_Integral(const Array<OneD, const NekDouble> &inarray) final;

    //----------------------------
    // Differentiation Methods
    //----------------------------
    LOCAL_REGIONS_EXPORT void v_PhysDeriv(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &out_d0, Array<OneD, NekDouble> &out_d1,
        Array<OneD, NekDouble> &out_d2) final;

    //---------------------------------------
    // Transforms
    //---------------------------------------
    LOCAL_REGIONS_EXPORT void v_FwdTrans(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) final;

    //---------------------------------------
    // Inner product functions
    //---------------------------------------
    LOCAL_REGIONS_EXPORT void v_IProductWRTBase(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) final;
    LOCAL_REGIONS_EXPORT void v_IProductWRTBase_SumFac(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, bool multiplybyweights = true) final;
    LOCAL_REGIONS_EXPORT void v_IProductWRTDerivBase(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) final;
    LOCAL_REGIONS_EXPORT void v_IProductWRTDerivBase_SumFac(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) final;

    LOCAL_REGIONS_EXPORT void v_AlignVectorToCollapsedDir(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray) final;

    //---------------------------------------
    // Evaluation functions
    //---------------------------------------
    LOCAL_REGIONS_EXPORT
    StdRegions::StdExpansionSharedPtr v_GetStdExp(void) const final;

    LOCAL_REGIONS_EXPORT
    StdRegions::StdExpansionSharedPtr v_GetLinStdExp(void) const final;

    LOCAL_REGIONS_EXPORT void v_GetCoord(
        const Array<OneD, const NekDouble> &Lcoords,
        Array<OneD, NekDouble> &coords) final;

    LOCAL_REGIONS_EXPORT void v_GetCoords(
        Array<OneD, NekDouble> &coords_1, Array<OneD, NekDouble> &coords_2,
        Array<OneD, NekDouble> &coords_3) final;

    LOCAL_REGIONS_EXPORT void v_ExtractDataToCoeffs(
        const NekDouble *data, const std::vector<unsigned int> &nummodes,
        const int mode_offset, NekDouble *coeffs,
        std::vector<LibUtilities::BasisType> &fromType) final;

    LOCAL_REGIONS_EXPORT NekDouble
    v_StdPhysEvaluate(const Array<OneD, const NekDouble> &Lcoord,
                      const Array<OneD, const NekDouble> &physvals);

    LOCAL_REGIONS_EXPORT NekDouble
    v_PhysEvaluate(const Array<OneD, const NekDouble> &coord,
                   const Array<OneD, const NekDouble> &physvals) final;

    LOCAL_REGIONS_EXPORT NekDouble
    v_PhysEvaluate(const Array<OneD, NekDouble> &coord,
                   const Array<OneD, const NekDouble> &inarray,
                   std::array<NekDouble, 3> &firstOrderDerivs) final;

    //---------------------------------------
    // Helper functions
    //---------------------------------------
    LOCAL_REGIONS_EXPORT int v_GetCoordim() final;

    LOCAL_REGIONS_EXPORT void v_GetTracePhysMap(
        const int face, Array<OneD, int> &outarray) final;

    LOCAL_REGIONS_EXPORT void v_ComputeTraceNormal(const int face) final;

    LOCAL_REGIONS_EXPORT void v_SVVLaplacianFilter(
        Array<OneD, NekDouble> &array,
        const StdRegions::StdMatrixKey &mkey) final;

    //---------------------------------------
    // Matrix creation functions
    //---------------------------------------
    LOCAL_REGIONS_EXPORT DNekMatSharedPtr
    v_GenMatrix(const StdRegions::StdMatrixKey &mkey) final;
    LOCAL_REGIONS_EXPORT DNekMatSharedPtr
    v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey) final;
    LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr
    v_GetLocMatrix(const MatrixKey &mkey) final;
    LOCAL_REGIONS_EXPORT DNekScalBlkMatSharedPtr
    v_GetLocStaticCondMatrix(const MatrixKey &mkey) final;
    LOCAL_REGIONS_EXPORT void v_DropLocStaticCondMatrix(
        const MatrixKey &mkey) final;
    LOCAL_REGIONS_EXPORT void v_ComputeLaplacianMetric() final;
    LOCAL_REGIONS_EXPORT void v_NormalTraceDerivFactors(
        Array<OneD, Array<OneD, NekDouble>> &d0factors,
        Array<OneD, Array<OneD, NekDouble>> &d1factors,
        Array<OneD, Array<OneD, NekDouble>> &d2factors) final;

private:
    LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess>
        m_matrixManager;
    LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess>
        m_staticCondMatrixManager;

    void v_LaplacianMatrixOp_MatFree_Kernel(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, Array<OneD, NekDouble> &wsp) final;
};

typedef std::shared_ptr<PyrExp> PyrExpSharedPtr;
typedef std::vector<PyrExpSharedPtr> PyrExpVector;
} // namespace LocalRegions
} // namespace Nektar

#endif // PYREXP_H
