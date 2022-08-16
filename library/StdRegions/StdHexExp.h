///////////////////////////////////////////////////////////////////////////////
//
// File StdHexExp.h
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
// Description: Header file for hexahedron routines.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBRARY_STDREGIONS_STDHEXEXP_H
#define NEKTAR_LIBRARY_STDREGIONS_STDHEXEXP_H

//#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion3D.h>
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar
{
namespace StdRegions
{
/// Class representing a hexehedral element in reference space.
class StdHexExp : virtual public StdExpansion3D
{

public:
    STD_REGIONS_EXPORT StdHexExp() = default;
    STD_REGIONS_EXPORT StdHexExp(const LibUtilities::BasisKey &Ba,
                                 const LibUtilities::BasisKey &Bb,
                                 const LibUtilities::BasisKey &Bc);
    STD_REGIONS_EXPORT StdHexExp(const StdHexExp &T);
    STD_REGIONS_EXPORT ~StdHexExp() override = default;

protected:
    //----------------------------
    // Differentiation Methods
    //----------------------------
    STD_REGIONS_EXPORT void v_PhysDeriv(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &out_d0, Array<OneD, NekDouble> &out_d1,
        Array<OneD, NekDouble> &out_d2) override;
    STD_REGIONS_EXPORT void v_PhysDeriv(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;
    STD_REGIONS_EXPORT void v_StdPhysDeriv(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &out_d0, Array<OneD, NekDouble> &out_d1,
        Array<OneD, NekDouble> &out_d2) override;
    STD_REGIONS_EXPORT void v_StdPhysDeriv(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    //---------------------------------------
    // Transforms
    //---------------------------------------
    STD_REGIONS_EXPORT void v_BwdTrans(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    STD_REGIONS_EXPORT void v_BwdTrans_SumFac(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    STD_REGIONS_EXPORT void v_BwdTrans_SumFacKernel(
        const Array<OneD, const NekDouble> &base0,
        const Array<OneD, const NekDouble> &base1,
        const Array<OneD, const NekDouble> &base2,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, Array<OneD, NekDouble> &wsp,
        bool doCheckCollDir0, bool doCheckCollDir1,
        bool doCheckCollDir2) override;

    STD_REGIONS_EXPORT void v_FwdTrans(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    //---------------------------------------
    // Inner product functions
    //---------------------------------------
    STD_REGIONS_EXPORT void v_IProductWRTBase(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    STD_REGIONS_EXPORT virtual void v_IProductWRTBase_MatOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray);

    STD_REGIONS_EXPORT void v_IProductWRTBase_SumFac(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, bool multbyweights = true) override;

    STD_REGIONS_EXPORT void v_IProductWRTBase_SumFacKernel(
        const Array<OneD, const NekDouble> &base0,
        const Array<OneD, const NekDouble> &base1,
        const Array<OneD, const NekDouble> &base2,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, Array<OneD, NekDouble> &wsp,
        bool doCheckCollDir0, bool doCheckCollDir1,
        bool doCheckCollDir2) override;

    STD_REGIONS_EXPORT void v_IProductWRTDerivBase(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    STD_REGIONS_EXPORT virtual void v_IProductWRTDerivBase_MatOp(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray);

    STD_REGIONS_EXPORT void v_IProductWRTDerivBase_SumFac(
        const int dir, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    //---------------------------------------
    // Evaluation functions
    //---------------------------------------
    STD_REGIONS_EXPORT void v_LocCoordToLocCollapsed(
        const Array<OneD, const NekDouble> &xi,
        Array<OneD, NekDouble> &eta) override;
    STD_REGIONS_EXPORT void v_LocCollapsedToLocCoord(
        const Array<OneD, const NekDouble> &eta,
        Array<OneD, NekDouble> &xi) override;
    STD_REGIONS_EXPORT void v_FillMode(
        const int mode, Array<OneD, NekDouble> &outarray) override;
    STD_REGIONS_EXPORT NekDouble v_PhysEvaluateBasis(
        const Array<OneD, const NekDouble> &coords, int mode) final;
    STD_REGIONS_EXPORT NekDouble
    v_PhysEvaluate(const Array<OneD, NekDouble> &coord,
                   const Array<OneD, const NekDouble> &inarray,
                   std::array<NekDouble, 3> &firstOrderDerivs) override;

    //---------------------------
    // Helper functions
    //---------------------------
    STD_REGIONS_EXPORT int v_GetNverts() const override;
    STD_REGIONS_EXPORT int v_GetNedges() const override;
    STD_REGIONS_EXPORT int v_GetNtraces() const override;
    STD_REGIONS_EXPORT LibUtilities::ShapeType v_DetShapeType() const override;
    STD_REGIONS_EXPORT int v_NumBndryCoeffs() const override;
    STD_REGIONS_EXPORT int v_NumDGBndryCoeffs() const override;
    STD_REGIONS_EXPORT int v_GetTraceNcoeffs(const int i) const override;
    STD_REGIONS_EXPORT int v_GetTotalTraceIntNcoeffs() const override;
    STD_REGIONS_EXPORT int v_GetTraceIntNcoeffs(const int i) const override;
    STD_REGIONS_EXPORT int v_GetTraceNumPoints(const int i) const override;

    STD_REGIONS_EXPORT LibUtilities::PointsKey v_GetTracePointsKey(
        const int i, const int j) const override;
    STD_REGIONS_EXPORT int v_CalcNumberOfCoefficients(
        const std::vector<unsigned int> &nummodes, int &modes_offset) override;
    STD_REGIONS_EXPORT const LibUtilities::BasisKey v_GetTraceBasisKey(
        const int i, const int k) const override;
    STD_REGIONS_EXPORT bool v_IsBoundaryInteriorExpansion() override;
    STD_REGIONS_EXPORT void v_GetCoords(
        Array<OneD, NekDouble> &coords_x, Array<OneD, NekDouble> &coords_y,
        Array<OneD, NekDouble> &coords_z) override;
    STD_REGIONS_EXPORT void v_GetTraceNumModes(
        const int fid, int &numModes0, int &numModes1,
        Orientation faceOrient = eDir1FwdDir1_Dir2FwdDir2) override;

    STD_REGIONS_EXPORT int v_GetEdgeNcoeffs(const int i) const override;

    //--------------------------
    // Mappings
    //--------------------------
    STD_REGIONS_EXPORT int v_GetVertexMap(
        int localVertexId, bool useCoeffPacking = false) override;
    STD_REGIONS_EXPORT void v_GetInteriorMap(
        Array<OneD, unsigned int> &outarray) override;
    STD_REGIONS_EXPORT void v_GetBoundaryMap(
        Array<OneD, unsigned int> &outarray) override;
    STD_REGIONS_EXPORT void v_GetTraceCoeffMap(
        const unsigned int fid, Array<OneD, unsigned int> &maparray) override;
    STD_REGIONS_EXPORT void v_GetElmtTraceToTraceMap(
        const unsigned int fid, Array<OneD, unsigned int> &maparray,
        Array<OneD, int> &signarray, Orientation faceOrient, int P,
        int Q) override;
    STD_REGIONS_EXPORT void v_GetEdgeInteriorToElementMap(
        const int tid, Array<OneD, unsigned int> &maparray,
        Array<OneD, int> &signarray,
        const Orientation traceOrient = eDir1FwdDir1_Dir2FwdDir2) override;

    STD_REGIONS_EXPORT void v_GetTraceInteriorToElementMap(
        const int tid, Array<OneD, unsigned int> &maparray,
        Array<OneD, int> &signarray,
        const Orientation traceOrient = eDir1FwdDir1_Dir2FwdDir2) override;

    //---------------------------------------
    // Wrapper functions
    //---------------------------------------
    STD_REGIONS_EXPORT DNekMatSharedPtr
    v_GenMatrix(const StdMatrixKey &mkey) override;
    STD_REGIONS_EXPORT DNekMatSharedPtr
    v_CreateStdMatrix(const StdMatrixKey &mkey) override;

    //---------------------------------------
    // Operator evaluation functions
    //---------------------------------------
    STD_REGIONS_EXPORT void v_MassMatrixOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey) override;

    STD_REGIONS_EXPORT void v_LaplacianMatrixOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey) override;

    STD_REGIONS_EXPORT void v_LaplacianMatrixOp(
        const int k1, const int k2, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey) override;

    STD_REGIONS_EXPORT void v_WeakDerivMatrixOp(
        const int i, const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey) override;

    STD_REGIONS_EXPORT void v_HelmholtzMatrixOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey) override;

    STD_REGIONS_EXPORT virtual void v_GeneralMatrixOp_MatOp(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray, const StdMatrixKey &mkey);

    STD_REGIONS_EXPORT void v_MultiplyByStdQuadratureMetric(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    STD_REGIONS_EXPORT void v_SVVLaplacianFilter(
        Array<OneD, NekDouble> &array, const StdMatrixKey &mkey) override;

    STD_REGIONS_EXPORT void v_ExponentialFilter(
        Array<OneD, NekDouble> &array, const NekDouble alpha,
        const NekDouble exponent, const NekDouble cutoff) override;
};

typedef std::shared_ptr<StdHexExp> StdHexExpSharedPtr;

} // namespace StdRegions
} // namespace Nektar

#endif
