///////////////////////////////////////////////////////////////////////////////
//
// File SmoothieSIAC.h
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
// Description: SmoothieSIAC definition
//
///////////////////////////////////////////////////////////////////////////////
#pragma once

#include "HandleNekMesh.h"
#include "NonSymmetricSIAC.h"
#include "SIACFilter.h"
#include "Smoothie.h"
#include "SymmetricSIAC.h"
#include <LibUtilities/Foundations/ManagerAccess.h> // for Points Manager, etc
#include <StdRegions/StdSegExp.h>
#include <memory>

#include <algorithm>
#include <boost/timer.hpp>

using namespace Nektar::LSIAC::SIACUtilities;

namespace Nektar
{
namespace LSIAC
{
/** @brief High level class that can apply pre-determined filters.
 */
class SmoothieSIAC : public Smoothie
{
public:
    enum filter_flags
    {
        DERIVATIVE_AFTER_SIAC =
            0x001, /**< Flag for Derivative to be applied after SIAC filter is
                      applied at each point */
        DERIVATIVE_USING_CONV_SIAC =
            0x002, /**< Flag to apply derivative on filter itself. */
        DERIVATIVE_RESERVED_FUTURE1 =
            0x004, /**< Flag reserved for future use. */
        DERIVATIVE_RESERVED_FUTURE2 =
            0x008, /**< Flag reserved for future use. */
        SYMMETRIC_2kp1 =
            0x010, /**< Flag to apply Symmetric filter with 2kp1 order. */
        SYMMETRIC_4kp1 =
            0x020, /**< Flag to apply Symmetric filter with 4kp1 order. */
        SYMMETRIC_FUTURE_1 = 0x040, /**< Flag reserved for future use. */
        SYMMETRIC_FUTURE_2 = 0x080, /**< Flag reserved for future use. */
        ONESIDED_2kPlus1 =
            0x100, /**< Flag to apply Symmetric filter with 2kp1 order. */
        ONESIDED_2kPlus2 = 0x200, /**< Flag to apply Symmetric filter with 2kp2
                                     order with one last general B-Spline. */
        ONESIDED_4kPlus1 =
            0x400, /**< Flag to apply Symmetric filter with 4kp1 order. */
        ONESIDED_FUTURE_1 = 0x800 /**< Flag reserved for future use. */
    };

private:
protected:
public:
    FilterType m_filterType;
    int m_order;
    unique_ptr<SymmetricSIAC> m_symSIACptr;
    vector<std::unique_ptr<SIACFilter>> m_siacFilterPtrs;
    HandleNekMesh *m_meshHandlePtr;
    int m_SymID = 0;
    int m_OneID = 1; // These should not be set until intitalization. Ideal for
                     // extra check.
public:
    SmoothieSIAC(FilterType filter, int order)
        : m_filterType(filter), m_order(order),
          m_symSIACptr(new SymmetricSIAC(order)){};

    SmoothieSIAC(const FilterType filter, const HandleMesh &meshHandle,
                 const int Order, Array<OneD, NekDouble> &direction);

    /**
     * @brief Apply L-SIAC at a given point in the mesh.
     *
     * This function takes the input location using ptsX, ptsY, ptsZ and output
     * variables in valX, valY, valZ.
     *
     * Note:
     * 1) Currently, the L-SIAC filter can only be applied to scalar fields and
     * produce scalar outputs; hence only valX, is being used, and valY and valZ
     * always result in 0.0.
     */
    bool EvaluateAt(const NekDouble ptsX, const NekDouble ptsY,
                    const NekDouble ptsZ, NekDouble &valX, NekDouble &valY,
                    NekDouble &valZ)
    {
        return v_EvaluateAt(ptsX, ptsY, ptsZ, valX, valY, valZ);
    }

    /**
     * @brief Research Phase.
     *
     * In the research phase, the behavior may change in time.
     *
     */
    bool EvaluateNonSymAt(const NekDouble ptsX, const NekDouble ptsY,
                          const NekDouble ptsZ, NekDouble &valX,
                          NekDouble &valY, NekDouble &valZ)
    {
        return v_EvaluateNonSymAt(ptsX, ptsY, ptsZ, valX, valY, valZ);
    }

    /**
     * @brief Research Phase.
     *
     * In the research phase, the behavior may change in time.
     *
     */
    bool EvaluateRecursiveAt(
        const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        const vector<std::shared_ptr<SmoothieSIAC>> &Sms,
        vector<Array<OneD, NekDouble>> directions,
        const vector<NekDouble> &meshSpacings = vector<NekDouble>(),
        const vector<int> &varNums = vector<int>(), const int curLevel = 0)
    {
        return v_EvaluateRecursiveAt(ptsX, ptsY, ptsZ, valX, valY, valZ, Sms,
                                     directions, meshSpacings, varNums,
                                     curLevel);
    }

    /**
     * @brief Research Phase.
     *
     * In the research phase, the behavior may change in time.
     *
     */
    bool EvaluateAt_NSK_GivenFilterLength_v1(
        const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        Array<OneD, NekDouble> &direction, NekDouble meshSpacing = -1.0,
        int varNum = 0)
    {
        return v_EvaluateAt_NSK_GivenFilterLength_v1(
            ptsX, ptsY, ptsZ, valX, valY, valZ, direction, meshSpacing, varNum);
    }

    /**
     * @brief Apply L-SIAC using non uniform knot sequence.
     *
     * In the research phase, the behavior may change in time.
     */
    bool EvaluateAt_NSK_GivenFilterLength_v2(
        const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        Array<OneD, NekDouble> &direction, NekDouble meshSpacing = -1.0,
        int varNum = 0)
    {
        return v_EvaluateAt_NSK_GivenFilterLength_v2(
            ptsX, ptsY, ptsZ, valX, valY, valZ, direction, meshSpacing, varNum);
    }

    /**
     * @brief Apply L-SIAC at a given point the mesh.
     *
     * This function takes the input location using ptsX, ptsY, ptsZ and output
     * variables in valX, valY, valZ.
     * The L-SIAC parameters direction and characteristic lengths can be set
     * using direction and meshSpacing.
     *
     * Note:
     * 1) Currently L-SIAC filter can only be applied to scalar fields and
     * produce scalar outputs; hence only valX is being used and valY and valZ
     * always result 0.0.
     * 2) If meshSpacing is negative, the meshSpacing specified along with the
     * constructor is used.
     */
    bool EvaluateAt(const NekDouble ptsX, const NekDouble ptsY,
                    const NekDouble ptsZ, NekDouble &valX, NekDouble &valY,
                    NekDouble &valZ, Array<OneD, NekDouble> &direction,
                    NekDouble meshSpacing = -1.0, int varNum = 0)
    {
        return v_EvaluateAt(ptsX, ptsY, ptsZ, valX, valY, valZ, direction,
                            meshSpacing, varNum);
    }

    /**
     * @brief Apply L-SIAC at a given point in the mesh.
     *
     * This function takes the input location using ptsX, ptsY, ptsZ and output
     * variables in valX, valY, valZ.
     * The L-SIAC parameters direction and characteristic lengths can be set
     * using directions and meshSpacing. This function
     * loops through an array of directions to find a direction in which the
     * symmetric L-SIAC filter can be applied. If the symmetric filter cannot be
     * applied, then the last direction is used to apply a one-sided filter.
     *
     * Note:
     * 1) Currently, the L-SIAC filter can only be applied to scalar fields and
     * produce scalar outputs; hence only valX is being used and valY and valZ
     * always result in 0.0.
     * 2) If meshSpacing is negative, the meshSpacing specified along with the
     * constructor is used.
     */
    bool EvaluateAt(const NekDouble ptsX, const NekDouble ptsY,
                    const NekDouble ptsZ, NekDouble &valX, NekDouble &valY,
                    NekDouble &valZ, vector<Array<OneD, NekDouble>> &directions,
                    NekDouble meshSpacing = -1.0, int varNum = 0)
    {
        return v_EvaluateAt(ptsX, ptsY, ptsZ, valX, valY, valZ, directions,
                            meshSpacing, varNum);
    }

    /**
     * @brief Apply L-SIAC using non uniform knot sequence.
     *
     * In the research phase, the behavior may change in time.
     */
    bool EvaluateAt_NSK_FixedNumBSpl(const NekDouble ptsX, const NekDouble ptsY,
                                     const NekDouble ptsZ, NekDouble &valX,
                                     NekDouble &valY, NekDouble &valZ,
                                     Array<OneD, NekDouble> &direction,
                                     NekDouble meshSpacing = -1.0,
                                     int varNum            = 0)
    {
        return v_EvaluateAt_NSK_FixedNumBSpl(ptsX, ptsY, ptsZ, valX, valY, valZ,
                                             direction, meshSpacing, varNum);
    }

    /**
     * @brief Apply L-SIAC using non uniform knot sequence.
     *
     * In the research phase, the behavior may change in time.
     */
    bool EvaluateAt_NUK_MetricTensor(const NekDouble ptsX, const NekDouble ptsY,
                                     const NekDouble ptsZ, NekDouble &valX,
                                     NekDouble &valY, NekDouble &valZ,
                                     Array<OneD, NekDouble> &direction,
                                     NekDouble meshSpacing = -1.0,
                                     int varNum            = 0)
    {
        return v_EvaluateAt_NUK_MetricTensor(ptsX, ptsY, ptsZ, valX, valY, valZ,
                                             direction, meshSpacing, varNum);
    }

    /**
     * @brief Apply L-SIAC using non uniform knot sequence.
     *
     * In the research phase, the behavior may change in time.
     */
    bool EvaluateAt_SymY(const NekDouble ptsX, const NekDouble ptsY,
                         const NekDouble ptsZ, NekDouble &valX, NekDouble &valY,
                         NekDouble &valZ, Array<OneD, NekDouble> &direction,
                         NekDouble meshSpacing = -1.0, int varNum = 0)
    {
        return v_EvaluateAt_SymY(ptsX, ptsY, ptsZ, valX, valY, valZ, direction,
                                 meshSpacing, varNum);
    }

    /**
     * @brief Research Phase. Post process at list of specified locations.
     *
     * In the research phase, the behavior may change in time.
     * This function is useful to design a better algorithm at post
     * processing all locations at once, instead of calculating one point at a
     * time.
     */
    bool EvaluateAt(const Array<OneD, NekDouble> &listPtsX,
                    const Array<OneD, NekDouble> &listPtsY,
                    const Array<OneD, NekDouble> &listPtsZ,
                    Array<OneD, NekDouble> &valX, Array<OneD, NekDouble> &valY,
                    Array<OneD, NekDouble> &valZ);

    /**
     * @brief Construct a 1D mesh using a line and the given mesh.
     *
     * Applicable only for 2D and 3D meshes. This function constructs a 1D mesh
     * across the line specified by the parameters (stPoint, direction) over
     * the higher dimensional mesh. The vertices of the 1D mesh are the
     * intersections of the line with the higher dimensional mesh.
     *
     * Output:
     * The HvalT has the locations of the vertices for the 1D mesh.
     * The element ids of the higher dimensional mesh are stored in t_EIDs.
     * The global element ids of the higher dimensional mesh are stored in
     * t_GIDs (currently t_EIDs are t_GIDs are the same).
     * The location of the quadrature points in each 1D element is stored in
     * t_LineElm and each element as t_quadPts.
     */
    bool SetupLineForLSIAC(const Array<OneD, NekDouble> &direction,
                           const vector<NekDouble> &stPoint,
                           const NekDouble tmin, const NekDouble tmax,
                           const int n_quadPts, vector<NekDouble> &HvalT,
                           vector<int> &t_GIDs, vector<int> &t_EIDs,
                           Array<OneD, NekDouble> &t_LineElm)
    {
        return v_SetupLineForLSIAC(direction, stPoint, tmin, tmax, n_quadPts,
                                   HvalT, t_GIDs, t_EIDs, t_LineElm);
    }

    /**
     * @brief Construct a 1D mesh using a line and the given mesh.
     *
     * Applicable only for 2D and 3D meshes. This function constructs a 1D mesh
     * across the line specified by the parameters (stPoint, direction) over
     * the higher dimensional mesh. The vertices of the 1D mesh are the
     * intersections of the line with the higher dimensional mesh.
     *
     * Output:
     * The HvalT has the locations of the vertices for the 1D mesh.
     * The element ids of the higher dimensional mesh are stored in t_EIDs.
     * The global element ids of the higher dimensional mesh are stored in
     * t_GIDs (currently t_EIDs are t_GIDs are the same).
     * The locations of quadrature points with t_quadPts and t_quadPts_Resample
     * points for each element are stored in t_LineElm, and t_LineElm_Resample,
     * respectively.
     */
    bool SetupLineForLSIAC_ReSamp(const Array<OneD, NekDouble> &direction,
                                  const vector<NekDouble> &stPoint,
                                  const NekDouble tmin, const NekDouble tmax,
                                  const int n_quadPts,
                                  const int n_quadPts_Resample,
                                  vector<NekDouble> &HvalT, vector<int> &t_GIDs,
                                  vector<int> &t_EIDs,
                                  Array<OneD, NekDouble> &t_LineElm,
                                  Array<OneD, NekDouble> &t_LineElm_Resample)
    {
        return v_SetupLineForLSIAC_ReSamp(
            direction, stPoint, tmin, tmax, n_quadPts, n_quadPts_Resample,
            HvalT, t_GIDs, t_EIDs, t_LineElm, t_LineElm_Resample);
    }

    /**
     * @brief Project the higher dimension field on the 1D mesh.
     *
     * Applicable only for 2D and 3D meshes. This function projects the fields
     * (varNUM) across the line specified by the parameters (stPoint,
     * direction).
     *
     * Output:
     * The value of the field at the locations specified by t_LineElm are
     * stored in tv_LineELm.
     */
    bool GetVLineForLSIAC(const int n_quadPts, const vector<NekDouble> &stPoint,
                          const Array<OneD, NekDouble> &direction,
                          const vector<NekDouble> &HvalT,
                          const vector<int> &t_GIDs, const vector<int> &t_EIDs,
                          const Array<OneD, NekDouble> &t_LineElm,
                          Array<OneD, NekDouble> tv_LineElm, int varNum)
    {
        return v_GetVLineForLSIAC(n_quadPts, stPoint, direction, HvalT, t_GIDs,
                                  t_EIDs, t_LineElm, tv_LineElm, varNum);
    }

    /**
     * @brief Project the higher dimension field on the 1D mesh.
     *
     * Applicable only for 2D and 3D meshes. This function projects the fields
     * (varNUM) across the line specified by the parameters (stPoint,
     * direction).
     *
     * Output:
     * The value of the field at the locations specified by t_LineElm_resample
     * are stored in tv_LineELm_resample_resample.
     */
    bool GetVLineForLSIAC_resample(
        const int n_quadPts, const int n_quadPts_resample,
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const vector<NekDouble> &HvalT,
        const vector<int> &t_GIDs, const vector<int> &t_EIDs,
        const Array<OneD, NekDouble> &t_LineElm,
        const Array<OneD, NekDouble> tv_LineElm,
        const Array<OneD, NekDouble> &t_LineElm_resample,
        Array<OneD, NekDouble> tv_LineElm_resample, int varNum)
    {
        return v_GetVLineForLSIAC_resample(
            n_quadPts, n_quadPts_resample, stPoint, direction, HvalT, t_GIDs,
            t_EIDs, t_LineElm, tv_LineElm, t_LineElm_resample,
            tv_LineElm_resample, varNum);
    }

    /**
     *  @brief Return the adaptive characteristic lengths for a set of ponits.
     */
    bool GetDynScalingForLSIAC(const vector<NekDouble> &stPoint,
                               const Array<OneD, NekDouble> &direction,
                               const vector<NekDouble> &tparams,
                               vector<NekDouble> &t_dynScaling,
                               const NekDouble meshSpacing, const NekDouble mu)
    {
        return v_GetDynScalingForLSIAC(stPoint, direction, tparams,
                                       t_dynScaling, meshSpacing, mu);
    }

    /**
     *  @brief Evaluate LSIAC filter at a set of points across a line.
     *
     *  \param n_quadPts Number of quadrature points used for integration.
     *  \param tparams Locations of the points for evaluation.
     *  \param t_dynScaling The adaptive characteristic lenghts at the points.
     *  \param t_mesh_min The minimum value of line parameter within the mesh
     *  boundary.
     *  \param t_mesh_max The maximum value of line parameter within the mesh
     *  boundary.
     *  \param t_LineElm The locations of quadrature points across the line for
     *  each element.
     *  \param t_vLineElm The value of the field at the quadrature points
     *  specified by t_LineElm.
     *  \param HvalT The locations of vertices across the line.
     *  \param tvals The postprocessed values at locations tparams are written
     *  to this variable.
     */
    bool EvaluateLineForLSIAC_v1_dynScaling(
        const int n_quadPts, const vector<NekDouble> &tparams,
        const vector<NekDouble> &t_dynScaling, const NekDouble t_mesh_min,
        const NekDouble t_mesh_max, const Array<OneD, NekDouble> &t_LineElm,
        const Array<OneD, NekDouble> &tv_LineElm,
        const vector<NekDouble> &HvalT, vector<NekDouble> &tvals)
    {
        return v_EvaluateLineForLSIAC_v1_dynScaling(
            n_quadPts, tparams, t_dynScaling, t_mesh_min, t_mesh_max, t_LineElm,
            tv_LineElm, HvalT, tvals);
    }

    /**
     *  @brief Evaluate LSIAC filter at a set of points across a line.
     *
     *  \param n_quadPts Number of quadrature points used for integration.
     *  \param tparams Locations of the points for evaluation.
     *  \param meshSpacing This characteristic length is used for all the
     *  points .
     *  \param t_mesh_min The minimum value of line parameter within the mesh
     *  boundary.
     *  \param t_mesh_max The maximum value of line parameter within the mesh
     *  boundary.
     *  \param t_LineElm The locations of quadrature points across the line for
     *  each element.
     *  \param t_vLineElm The value of the field at the quadrature points
     *  specified by t_LineElm.
     *  \param HvalT The locations of vertices across the line.
     *  \param tvals The postprocessed values at locations tparams are written
     *  to this variable.
     */
    bool EvaluateLineForLSIAC_v1(
        const int n_quadPts, const vector<NekDouble> &tparams,
        const NekDouble meshSpacing, const NekDouble t_mesh_min,
        const NekDouble t_mesh_max, const Array<OneD, NekDouble> &t_LineElm,
        const Array<OneD, NekDouble> &tv_LineElm,
        const vector<NekDouble> &HvalT, vector<NekDouble> &tvals)
    {
        return v_EvaluateLineForLSIAC_v1(n_quadPts, tparams, meshSpacing,
                                         t_mesh_min, t_mesh_max, t_LineElm,
                                         tv_LineElm, HvalT, tvals);
    }

    /**
     *  @brief Evaluate LSIAC filter at a set of points across a line.
     *
     *  \param n_quadPts Number of quadrature points at which the field is
     *  defined for element.
     *  \param n_quadPts_resample Number of quadrature points used for
     *  integration.
     *  \param tparams Locations of the points for evaluation.
     *  \param meshSpacing This characteristic length is used for all the
     *  points .
     *  \param t_mesh_min The minimum value of line parameter within the mesh
     *  boundary.
     *  \param t_mesh_max The maximum value of line parameter within the mesh
     *  boundary.
     *  \param t_LineElm The locations of quadrature points across the line for
     *  each element.
     *  \param t_vLineElm The value of the field at the quadrature points
     *  specified by t_LineElm.
     *  \param HvalT The locations of vertices across the line.
     *  \param tvals The postprocessed values at locations tparams are written
     *  to this variable.
     */
    bool EvaluateLineForLSIAC_v3(
        const int n_quadPts, const int n_quadPts_resample,
        const vector<NekDouble> &tparams, const NekDouble meshSpacing,
        const NekDouble t_mesh_min, const NekDouble t_mesh_max,
        const Array<OneD, NekDouble> &t_LineElm,
        const Array<OneD, NekDouble> &tv_LineElm,
        const vector<NekDouble> &HvalT, vector<NekDouble> &tvals)
    {
        return v_EvaluateLineForLSIAC_v3(n_quadPts, n_quadPts_resample, tparams,
                                         meshSpacing, t_mesh_min, t_mesh_max,
                                         t_LineElm, tv_LineElm, HvalT, tvals);
    }

    /**
     *  @brief Evaluate LSIAC filter at a set of points across a line.
     *
     *  \param n_quadPts Number of quadrature points at which the field is
     *  defined for element.
     *  \param n_quadPts_resample Number of quadrature points used for
     *  integration.
     *  \param tparams Locations of the points for evaluation.
     *  \param t_dynScaling The adaptive characteristic lenghts at the points.
     *  \param t_mesh_min The minimum value of line parameter within the mesh
     *  boundary.
     *  \param t_mesh_max The maximum value of line parameter within the mesh
     *  boundary.
     *  \param t_LineElm The locations of quadrature points across the line for
     *  each element.
     *  \param t_vLineElm The value of the field at the quadrature points
     *  specified by t_LineElm.
     *  \param HvalT The locations of vertices across the line.
     *  \param tvals The postprocessed values at locations tparams are written
     *  to this variable.
     */
    bool EvaluateLineForLSIAC_v3_dynScaling(
        const int n_quadPts, const int n_quadPts_resample,
        const vector<NekDouble> &tparams, const vector<NekDouble> &t_dynScaling,
        const NekDouble t_mesh_min, const NekDouble t_mesh_max,
        const Array<OneD, NekDouble> &t_LineElm,
        const Array<OneD, NekDouble> &tv_LineElm,
        const vector<NekDouble> &HvalT, vector<NekDouble> &tvals)
    {
        return v_EvaluateLineForLSIAC_v3_dynScaling(
            n_quadPts, n_quadPts_resample, tparams, t_dynScaling, t_mesh_min,
            t_mesh_max, t_LineElm, tv_LineElm, HvalT, tvals);
    }

    /**
     * @brief Research Phase.
     *
     * In the research phase, the behavior may change in time.
     */
    bool EvaluateUsingLineAt(const vector<NekDouble> &stPoint,
                             const Array<OneD, NekDouble> &direction,
                             const int n_quadPts, const NekDouble meshScaling,
                             const vector<NekDouble> &tparams,
                             vector<NekDouble> &tvals, int varNum)
    {
        return v_EvaluateUsingLineAt(stPoint, direction, n_quadPts, meshScaling,
                                     tparams, tvals, varNum);
    }

    /**
     * @brief Research Phase.
     *
     * In the research phase, the behavior may change in time.
     */
    bool EvaluateUsingLineAt_v1DynScaling(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const NekDouble meshScaling, const vector<NekDouble> &tparams,
        vector<NekDouble> &tvals, int varNum)
    {
        return v_EvaluateUsingLineAt_v1DynScaling(
            stPoint, direction, n_quadPts, meshScaling, tparams, tvals, varNum);
    }

    /**
     * @brief Research Phase.
     *
     * In the research phase, the behavior may change in time.
     */
    bool EvaluateL2UsingLineAt(const vector<NekDouble> &stPoint,
                               const Array<OneD, NekDouble> &direction,
                               const int n_quadPts, const NekDouble meshScaling,
                               const vector<NekDouble> &tparams,
                               vector<NekDouble> &tvals, int varNum)
    {
        return v_EvaluateL2UsingLineAt(stPoint, direction, n_quadPts,
                                       meshScaling, tparams, tvals, varNum);
    }

    /**
     * @brief Research Phase.
     *
     * In the research phase, the behavior may change in time.
     */
    bool EvaluateUsingLineAt_v2DynScaling(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const int n_quadPts_resample, const NekDouble meshScaling,
        const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int varNum)
    {
        return v_EvaluateUsingLineAt_v2DynScaling(
            stPoint, direction, n_quadPts, n_quadPts_resample, meshScaling,
            tparams, tvals, varNum);
    }

    /**
     * @brief Research Phase.
     *
     * In the research phase, the behavior may change in time.
     */
    bool EvaluateUsingLineAt_v2(const vector<NekDouble> &stPoint,
                                const Array<OneD, NekDouble> &direction,
                                const int n_quadPts,
                                const int n_quadPts_resample,
                                const NekDouble meshScaling,
                                const vector<NekDouble> &tparams,
                                vector<NekDouble> &tvals, int varNum)
    {
        return v_EvaluateUsingLineAt_v2(stPoint, direction, n_quadPts,
                                        n_quadPts_resample, meshScaling,
                                        tparams, tvals, varNum);
    }

    /**
     * @brief Research Phase.
     *
     * In the research phase, the behavior may change in time.
     */
    bool EvaluateUsingLineAt_v3DynScaling(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const int n_quadPts_resample, const NekDouble meshScaling,
        const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int varNum)
    {
        return v_EvaluateUsingLineAt_v3DynScaling(
            stPoint, direction, n_quadPts, n_quadPts_resample, meshScaling,
            tparams, tvals, varNum);
    }

    /**
     * @brief Research Phase.
     *
     * In the research phase, the behavior may change in time.
     */
    bool EvaluateUsingLineAt_v3(const vector<NekDouble> &stPoint,
                                const Array<OneD, NekDouble> &direction,
                                const int n_quadPts,
                                const int n_quadPts_resample,
                                const NekDouble meshScaling,
                                const vector<NekDouble> &tparams,
                                vector<NekDouble> &tvals, int varNum)
    {
        return v_EvaluateUsingLineAt_v3(stPoint, direction, n_quadPts,
                                        n_quadPts_resample, meshScaling,
                                        tparams, tvals, varNum);
    }

    /**
     * @brief Research Phase.
     *
     * In the research phase, the behavior may change in time.
     */
    bool EvaluateL2UsingLineAt_v3(const vector<NekDouble> &stPoint,
                                  const Array<OneD, NekDouble> &direction,
                                  const int n_quadPts,
                                  const int n_quadPts_resample,
                                  const NekDouble meshScaling,
                                  const vector<NekDouble> &tparams,
                                  vector<NekDouble> &tvals, int varNum)
    {
        return v_EvaluateL2UsingLineAt_v3(stPoint, direction, n_quadPts,
                                          n_quadPts_resample, meshScaling,
                                          tparams, tvals, varNum);
    }

    /**
     * @brief Research Phase.
     *
     * In the research phase, the behavior may change in time.
     */
    bool EvaluateUsingLineAt_vNonSymKnots(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const NekDouble meshSpacing, const vector<NekDouble> &tparams,
        vector<NekDouble> &tvals, int varNum)
    {
        return v_EvaluateUsingLineAt_vNonSymKnots(
            stPoint, direction, n_quadPts, meshSpacing, tparams, tvals, varNum);
    }

    //		bool b_symMesh = m_meshHandlePtr->GetKnotVec(m_order-1,
    // t,HvaT,knotVec, nonSymShift); TBD
    bool CalculateKnotVec(NekDouble t, vector<NekDouble> &HvalT,
                          Array<OneD, NekDouble> &knotVec,
                          NekDouble &nonSymShift);

    /**
     * @brief Research Phase.
     *
     * In the research phase, the behavior may change in time.
     */
    bool EvaluatePt_vNonSymKnots(const vector<NekDouble> &stPoint,
                                 const Array<OneD, NekDouble> &direction,
                                 const int n_quadPts,
                                 const NekDouble meshSpacing,
                                 vector<NekDouble> &tvals, int varNum)
    {
        return v_EvaluatePt_vNonSymKnots(stPoint, direction, n_quadPts,
                                         meshSpacing, tvals, varNum);
    }

    /**
     * @brief Research Phase.
     *
     * In the research phase, the behavior may change in time.
     */
    bool Cal_NUK_ConstMetricTensor(const NekDouble ptsX, const NekDouble ptsY,
                                   const NekDouble ptsZ,
                                   const NekDouble meshSpacing,
                                   Array<OneD, NekDouble> &direction,
                                   Array<OneD, NekDouble> &knotVec)
    {
        return v_Cal_NUK_ConstMetricTensor(ptsX, ptsY, ptsZ, meshSpacing,
                                           direction, knotVec);
    }

protected:
    //! Post process at a specified location.
    /*!
            \param ptsX x-coordinate.
            \param ptsY y-coordinate.
            \param ptsZ z-coordinate.
            \param valX scalar/vector1st(X) result.
            \param valY vector2nd(Y)  result.
            \param valZ vector3rd(Z)  result.
     */
    virtual bool v_EvaluateAt(const NekDouble ptsX, const NekDouble ptsY,
                              const NekDouble ptsZ, NekDouble &valX,
                              NekDouble &valY, NekDouble &valZ) = 0;

    virtual bool v_EvaluateNonSymAt(const NekDouble ptsX, const NekDouble ptsY,
                                    const NekDouble ptsZ, NekDouble &valX,
                                    NekDouble &valY, NekDouble &valZ)
    {
        boost::ignore_unused(ptsX, ptsY, ptsZ, valX, valY, valZ);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }

    virtual bool v_EvaluateAt(const NekDouble ptsX, const NekDouble ptsY,
                              const NekDouble ptsZ, NekDouble &valX,
                              NekDouble &valY, NekDouble &valZ,
                              Array<OneD, NekDouble> &direction,
                              NekDouble meshSpacing = -1.0, int varNum = 0)
    {
        boost::ignore_unused(ptsX, ptsY, ptsZ, valX, valY, valZ, direction,
                             meshSpacing, varNum);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }

    virtual bool v_EvaluateAt(const NekDouble ptsX, const NekDouble ptsY,
                              const NekDouble ptsZ, NekDouble &valX,
                              NekDouble &valY, NekDouble &valZ,
                              vector<Array<OneD, NekDouble>> &directions,
                              NekDouble meshSpacing = -1.0, int varNum = 0)
    {
        boost::ignore_unused(ptsX, ptsY, ptsZ, valX, valY, valZ, directions,
                             meshSpacing, varNum);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }
    virtual bool v_EvaluateAt_NSK_FixedNumBSpl(
        const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        Array<OneD, NekDouble> &direction, NekDouble meshSpacing = -1.0,
        int varNum = 0)
    {
        boost::ignore_unused(ptsX, ptsY, ptsZ, valX, valY, valZ, direction,
                             meshSpacing, varNum);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }

    virtual bool v_EvaluateAt_SymY(const NekDouble ptsX, const NekDouble ptsY,
                                   const NekDouble ptsZ, NekDouble &valX,
                                   NekDouble &valY, NekDouble &valZ,
                                   Array<OneD, NekDouble> &direction,
                                   NekDouble meshSpacing = -1.0, int varNum = 0)
    {
        boost::ignore_unused(ptsX, ptsY, ptsZ, valX, valY, valZ, direction,
                             meshSpacing, varNum);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }

    virtual bool v_EvaluateAt_NUK_MetricTensor(
        const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        Array<OneD, NekDouble> &direction, NekDouble meshSpacing = -1.0,
        int varNum = 0)
    {
        boost::ignore_unused(ptsX, ptsY, ptsZ, valX, valY, valZ, direction,
                             meshSpacing, varNum);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }

    virtual bool v_EvaluateRecursiveAt(
        const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        const vector<std::shared_ptr<SmoothieSIAC>> &Sms,
        vector<Array<OneD, NekDouble>> directions,
        const vector<NekDouble> &meshSpacing = vector<NekDouble>(),
        const vector<int> &varNum = vector<int>(), const int curLevel = 0)
    {
        boost::ignore_unused(ptsX, ptsY, ptsZ, valX, valY, valZ, Sms,
                             directions, meshSpacing, varNum, curLevel);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }

    virtual bool v_SetupLineForLSIAC(const Array<OneD, NekDouble> &direction,
                                     const vector<NekDouble> &stPoint,
                                     const NekDouble tmin, const NekDouble tmax,
                                     const int n_quadPts,
                                     vector<NekDouble> &HvalT,
                                     vector<int> &t_GIDs, vector<int> &t_EIDs,
                                     Array<OneD, NekDouble> &t_LineElm);
    virtual bool v_SetupLineForLSIAC_ReSamp(
        const Array<OneD, NekDouble> &direction,
        const vector<NekDouble> &stPoint, const NekDouble tmin,
        const NekDouble tmax, const int n_quadPts, const int n_quadPts_Resample,
        vector<NekDouble> &HvalT, vector<int> &t_GIDs, vector<int> &t_EIDs,
        Array<OneD, NekDouble> &t_LineElm,
        Array<OneD, NekDouble> &t_LineElm_Resample);
    virtual bool v_GetVLineForLSIAC(
        const int n_quadPts, const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const vector<NekDouble> &HvalT,
        const vector<int> &t_GIDs, const vector<int> &t_EIDs,
        const Array<OneD, NekDouble> &t_LineElm,
        Array<OneD, NekDouble> tv_LineElm, int varNum);
    virtual bool v_GetVLineForLSIAC_resample(
        const int n_quadPts, const int n_quadPts_resample,
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const vector<NekDouble> &HvalT,
        const vector<int> &t_GIDs, const vector<int> &t_EIDs,
        const Array<OneD, NekDouble> &t_LineElm,
        const Array<OneD, NekDouble> tv_LineElm,
        const Array<OneD, NekDouble> &t_LineElm_resample,
        Array<OneD, NekDouble> tv_LineElm_resample, int varNum);
    virtual bool v_GetDynScalingForLSIAC(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction,
        const vector<NekDouble> &tparams, vector<NekDouble> &t_dynScaling,
        const NekDouble meshSpacing, const NekDouble mu);
    bool mergeBreakPts(const vector<NekDouble> &HvalT,
                       const vector<NekDouble> &SvalT,
                       vector<NekDouble> &TvalT);
    bool mergeBreakPtsAdv(const vector<NekDouble> &HvalT,
                          const vector<NekDouble> &SvalT,
                          vector<NekDouble> &TvalT, const vector<int> &t_h_GIDs,
                          const vector<int> &t_h_EIDs, vector<int> &t_GIDs,
                          vector<int> &t_EIDs);
    virtual bool v_EvaluateLineForLSIAC_v1(
        const int n_quadPts, const vector<NekDouble> &tparams,
        const NekDouble meshSpacing, const NekDouble t_mesh_min,
        const NekDouble t_mesh_max, const Array<OneD, NekDouble> &t_LineElm,
        const Array<OneD, NekDouble> &tv_LineElm,
        const vector<NekDouble> &HvalT, vector<NekDouble> &tvals);
    virtual bool v_EvaluateLineForLSIAC_v1_dynScaling(
        const int n_quadPts, const vector<NekDouble> &tparams,
        const vector<NekDouble> &t_dynScaling, const NekDouble t_mesh_min,
        const NekDouble t_mesh_max, const Array<OneD, NekDouble> &t_LineElm,
        const Array<OneD, NekDouble> &tv_LineElm,
        const vector<NekDouble> &HvalT, vector<NekDouble> &tvals);

    virtual bool v_EvaluateLineForLSIAC_v3(
        const int n_quadPts, const int n_quadPts_resample,
        const vector<NekDouble> &tparams, const NekDouble meshSpacing,
        const NekDouble t_mesh_min, const NekDouble t_mesh_max,
        const Array<OneD, NekDouble> &t_LineElm,
        const Array<OneD, NekDouble> &tv_LineElm,
        const vector<NekDouble> &HvalT, vector<NekDouble> &tvals);

    virtual bool v_EvaluateLineForLSIAC_v3_dynScaling(
        const int n_quadPts, const int n_quadPts_resample,
        const vector<NekDouble> &tparams, const vector<NekDouble> &t_dynScaling,
        const NekDouble t_mesh_min, const NekDouble t_mesh_max,
        const Array<OneD, NekDouble> &t_LineElm,
        const Array<OneD, NekDouble> &tv_LineElm,
        const vector<NekDouble> &HvalT, vector<NekDouble> &tvals);

    virtual bool v_EvaluateUsingLineAt(const vector<NekDouble> &stPoint,
                                       const Array<OneD, NekDouble> &direction,
                                       const int n_quadPts,
                                       const NekDouble meshScaling,
                                       const vector<NekDouble> &tparams,
                                       vector<NekDouble> &tvals, int varNum);

    virtual bool v_EvaluateUsingLineAt_v1DynScaling(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const NekDouble meshScaling, const vector<NekDouble> &tparams,
        vector<NekDouble> &tvals, int varNum);

    virtual bool v_EvaluateL2UsingLineAt(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const NekDouble meshScaling, const vector<NekDouble> &tparams,
        vector<NekDouble> &tvals, int varNum)
    {
        boost::ignore_unused(stPoint, direction, n_quadPts, meshScaling,
                             tparams, tvals, varNum);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }

    virtual bool v_EvaluateUsingLineAt_v2DynScaling(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const int n_quadPts_resample, const NekDouble meshScaling,
        const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int varNum);

    virtual bool v_EvaluateUsingLineAt_v2(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const int n_quadPts_resample, const NekDouble meshScaling,
        const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int varNum);

    virtual bool v_EvaluateUsingLineAt_v3(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const int n_quadPts_resample, const NekDouble meshScaling,
        const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int varNum);

    virtual bool v_EvaluateUsingLineAt_v3DynScaling(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const int n_quadPts_resample, const NekDouble meshScaling,
        const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int varNum);

    virtual bool v_EvaluateL2UsingLineAt_v3(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const int n_quadPts_resample, const NekDouble meshScaling,
        const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int varNum)
    {
        boost::ignore_unused(stPoint, direction, n_quadPts, n_quadPts_resample,
                             meshScaling, tparams, tvals, varNum);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }

    virtual bool v_EvaluateUsingLineAt_vNonSymKnots(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const NekDouble meshSpacing, const vector<NekDouble> &tparams,
        vector<NekDouble> &tvals, int varNum);

    virtual bool v_EvaluatePt_vNonSymKnots(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const NekDouble meshSpacing, vector<NekDouble> &tvals, int varNum);

    virtual bool v_EvaluateAt_NSK_GivenFilterLength_v1(
        const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        Array<OneD, NekDouble> &direction, NekDouble meshSpacing = -1.0,
        int varNum = 0)
    {
        boost::ignore_unused(ptsX, ptsY, ptsZ, valX, valY, valZ, direction,
                             meshSpacing, varNum);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }

    virtual bool v_EvaluateAt_NSK_GivenFilterLength_v2(
        const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        Array<OneD, NekDouble> &direction, NekDouble meshSpacing, int varNum)
    {
        boost::ignore_unused(ptsX, ptsY, ptsZ, valX, valY, valZ, direction,
                             meshSpacing, varNum);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }

    virtual bool v_Cal_NUK_ConstMetricTensor(const NekDouble ptsX,
                                             const NekDouble ptsY,
                                             const NekDouble ptsZ,
                                             const NekDouble meshSpacing,
                                             Array<OneD, NekDouble> &direction,
                                             Array<OneD, NekDouble> &knotVec)
    {
        boost::ignore_unused(ptsX, ptsY, ptsZ, meshSpacing, direction, knotVec);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }
};
} // namespace LSIAC
} // namespace Nektar
