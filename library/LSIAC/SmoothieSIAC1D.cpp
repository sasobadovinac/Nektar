///////////////////////////////////////////////////////////////////////////////
//
// File SmoothieSIAC1D.cpp
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
// Description: SmoothieSIAC1D definition
//
///////////////////////////////////////////////////////////////////////////////
#include "SmoothieSIAC1D.h"
#include "OneSidedSIAC.h"
#include <LibUtilities/Foundations/ManagerAccess.h> // for Points Manager, etc

namespace Nektar
{
namespace LSIAC
{
SmoothieSIAC1D::SmoothieSIAC1D(const FilterType filter,
                               HandleNekMesh *meshHandle, const int order,
                               NekDouble meshSpacing, const int derivative)
    : SmoothieSIAC(filter, order), m_meshSpacing(meshSpacing)
{
    // Need to Create SIAC filter depending on filterType.
    // Set the filterType flags to be used in the code.
    m_meshHandlePtr = meshHandle;
    switch (filter)
    {
        case eNONE:
            m_OneID = -1;
            m_SymID = -1;
            break;
        case eSYM_2kp1:
            // Set up symmetric siac filter and any other parameters needed for
            // case scenario.
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(order));
            m_OneID = -1;
            break;
        case eSYM_2kp1_1SIDED_2kp1:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(order));
            m_siacFilterPtrs.emplace_back(new OneSidedSIAC(
                order, OneSidedSIAC::OneSidedFilterType::BASIC_SIAC_2kp1));
            break;
        case eSYM_2kp1_1SIDED_4kp1:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(order));
            m_siacFilterPtrs.emplace_back(new OneSidedSIAC(
                order, OneSidedSIAC::OneSidedFilterType::VAN_SIAC_4kp1));
            break;
        case eSYM_4kp1:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(order, 4 * order));
            m_OneID = -1;
            ASSERTL1(false, "symmetric 4k+1 filter is some how screwed up.");
            break;
        case eSYM_DER_2kp1_1SIDED_2kp1:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(
                order,
                SymmetricSIAC::SymFilterType::CUSTOM_SMOOTH_Derivative_SIAC,
                derivative));
            m_siacFilterPtrs.emplace_back(new OneSidedSIAC(
                order,
                OneSidedSIAC::OneSidedFilterType::Der_SMOOTH_BASIC_SIAC_2kp1,
                derivative));
            break;
        case eSYM_DER_2kp1_1SIDED_4kp1:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(
                order,
                SymmetricSIAC::SymFilterType::CUSTOM_SMOOTH_Derivative_SIAC,
                derivative));
            m_siacFilterPtrs.emplace_back(new OneSidedSIAC(
                order, OneSidedSIAC::OneSidedFilterType::Der_BASIC_SIAC_4kp1,
                derivative));
            break;
        case eSYM_2kp1_1SIDED_2kp2:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(order));
            m_siacFilterPtrs.emplace_back(new OneSidedSIAC(
                order, OneSidedSIAC::OneSidedFilterType::XLi_SIAC_2kp2));
            break;
        case eSYM_DER_2kp1_1SIDED_2kp2:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(
                order,
                SymmetricSIAC::SymFilterType::CUSTOM_SMOOTH_Derivative_SIAC,
                derivative));
            m_siacFilterPtrs.emplace_back(new OneSidedSIAC(
                order, OneSidedSIAC::OneSidedFilterType::Der_XLi_SIAC_2kp2,
                derivative));
            break;
        case eSYM_UNEVEN_2kp1:
            m_siacFilterPtrs.emplace_back(new NonSymmetricSIAC(order));
            m_OneID = -1;
            break;
        default:
            NEKERROR(ErrorUtil::efatal, "Filter not defined");
    }
}

bool SmoothieSIAC1D::v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                                  const NekDouble PtsZ, NekDouble &valX,
                                  NekDouble &valY, NekDouble &valZ,
                                  Array<OneD, NekDouble> &direction,
                                  NekDouble meshSpacing, int varNum)
{
    boost::ignore_unused(valY, valZ); // reserved for derivative filter.
    // Direction by default is set in x-direction.
    direction[0]         = 1.0;
    NekDouble meshTShift = 0.0;
    NekDouble tmin, tmax;
    vector<NekDouble> HvalX, HvalY, HvalZ, HvalT, SvalT, TvalT;
    // HNM to get list of break points.
    // check if a symmetric fitler can be applied.
    // in position, symmetric range

    // Get filter range given scaling.
    m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, tmin, tmax);

    bool b_symMesh = m_meshHandlePtr->CanTRangebeApplied(
        PtsX, PtsY, PtsZ, direction, tmin, tmax, meshTShift);
    if (b_symMesh)
    {
        m_siacFilterPtrs[m_SymID]->GetBreakPts(
            meshSpacing, SvalT); //"??? specify parameters"
    }
    else
    {
        SvalT.clear();
        if (m_OneID >= 0)
        { // OneSided Filter is defined and given by .
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, tmin, tmax);
            m_meshHandlePtr->CanTRangebeApplied(PtsX, PtsY, PtsZ, direction,
                                                tmin, tmax, meshTShift);
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, tmin, tmax,
                                                      meshTShift);
            m_siacFilterPtrs[m_OneID]->GetBreakPts(meshSpacing, SvalT,
                                                   meshTShift);
        }
        else
        { // OneSided Filter is not defined.
            valX = -1;
            return false;
        }
    }

    m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ, direction, tmin, tmax, HvalX,
                                 HvalY, HvalZ, HvalT);
    mergeBreakPts(HvalT, SvalT, TvalT);

    // Loop for sumation.
    NekDouble sum = 0.0;
    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.
    int quad_npoints =
        ceil((m_order + 3) / 2) +
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetTotPoints() + 2;
    LibUtilities::PointsKey quadPointsKey(
        quad_npoints, Nektar::LibUtilities::eGaussGaussLegendre);
    Array<OneD, NekDouble> quad_points =
        LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    Array<OneD, NekDouble> quad_weights =
        LibUtilities::PointsManager()[quadPointsKey]->GetW();

    Array<OneD, NekDouble> t_quad(quad_npoints), t_x(quad_npoints),
        t_y(quad_npoints), t_z(quad_npoints), t_quad_vals(quad_npoints),
        t_xyz_vals(quad_npoints);
    vector<int> t_GIDs, t_EIDs;
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, TvalT, t_GIDs,
                                   t_EIDs);

    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        // Create a seg element
        // Actual quadrature points needed are for order = polyorder+Bspline
        // order.
        NekDouble a = TvalT[i];
        NekDouble b = TvalT[i + 1];
        int gID     = t_GIDs[i];
        int eID     = t_EIDs[i];
        for (int j = 0; j < quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (quad_points[j] + 1.0) / 2.0 + a;
            t_x[j]    = PtsX + t_quad[j] * direction[0];
            t_y[j]    = PtsY + t_quad[j] * direction[1];
            t_z[j]    = PtsZ + t_quad[j] * direction[2];
        }
        // Evaluate filter  at these points.
        if (b_symMesh)
        {
            m_siacFilterPtrs[m_SymID]->EvaluateFilter(t_quad, t_quad_vals,
                                                      meshSpacing);
        }
        else
        {
            m_siacFilterPtrs[m_OneID]->EvaluateFilter(
                t_quad, t_quad_vals, meshSpacing, meshTShift, true);
        }
        m_meshHandlePtr->EvaluateAt(t_x, t_y, t_z, gID, eID, t_xyz_vals,
                                    varNum);
        NekDouble integral = 0;
        // integral
        for (int k = 0; k < quad_npoints; k++)
        {
            integral += 0.5 * quad_weights[k] * t_quad_vals[k] *
                        std::abs((b - a) * direction[0]) * t_xyz_vals[k];
        }
        sum += integral;
    }
    valX = sum;
    return true;
}

bool SmoothieSIAC1D::v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                                  const NekDouble PtsZ, NekDouble &valX,
                                  NekDouble &valY, NekDouble &valZ)
{
    boost::ignore_unused(valY, valZ); // reserved for derivative filter.
    Array<OneD, NekDouble> direction(3, 0.0);
    direction[0]         = 1.0;
    NekDouble meshTShift = 0.0;
    NekDouble tmin, tmax;
    vector<NekDouble> HvalX, HvalY, HvalZ, HvalT, SvalT, TvalT;
    // HNM to get list of break points.
    // check if a symmetric fitler can be applied.
    // in position, symmetric range

    // Get filter range given scaling.
    m_siacFilterPtrs[m_SymID]->GetFilterRange(m_meshSpacing, tmin, tmax);

    bool b_symMesh = m_meshHandlePtr->CanTRangebeApplied(
        PtsX, PtsY, PtsZ, direction, tmin, tmax, meshTShift);
    // Made an assumtion that a symmetric filter can be applied ???
    // return List of mesh breakpoints in that range.
    // call SSS/SSO
    if (b_symMesh)
    {
        m_siacFilterPtrs[m_SymID]->GetBreakPts(
            m_meshSpacing, SvalT); //"??? specify parameters"
    }
    else
    {
        SvalT.clear();
        if (m_OneID >= 0)
        { // OneSided Filter is defined and given by .
            m_siacFilterPtrs[m_OneID]->GetFilterRange(m_meshSpacing, tmin,
                                                      tmax);
            m_meshHandlePtr->CanTRangebeApplied(PtsX, PtsY, PtsZ, direction,
                                                tmin, tmax, meshTShift);
            m_siacFilterPtrs[m_OneID]->GetFilterRange(m_meshSpacing, tmin, tmax,
                                                      meshTShift);
            m_siacFilterPtrs[m_OneID]->GetBreakPts(
                m_meshSpacing, SvalT, meshTShift); //"??? specify parameters"
        }
        else
        { // OneSided Filter is not defined.
            valX = -1;
            return false;
        }
    }

    m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ, direction, tmin, tmax, HvalX,
                                 HvalY, HvalZ, HvalT);
    mergeBreakPts(HvalT, SvalT, TvalT);

    // Loop for sumation.
    NekDouble sum = 0.0;
    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.

    int quad_npoints =
        ceil((m_order + 3) / 2) +
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetTotPoints() + 2;
    LibUtilities::PointsKey quadPointsKey(
        quad_npoints, Nektar::LibUtilities::eGaussGaussLegendre);
    Array<OneD, NekDouble> quad_points =
        LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    Array<OneD, NekDouble> quad_weights =
        LibUtilities::PointsManager()[quadPointsKey]->GetW();

    Array<OneD, NekDouble> t_quad(quad_npoints), t_x(quad_npoints),
        t_y(quad_npoints), t_z(quad_npoints), t_quad_vals(quad_npoints),
        t_xyz_vals(quad_npoints);
    vector<int> t_GIDs, t_EIDs;
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, TvalT, t_GIDs,
                                   t_EIDs);

    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        // Create a seg element
        // Actual quadrature points needed are for order = polyorder+Bspline
        // order.
        NekDouble a = TvalT[i];
        NekDouble b = TvalT[i + 1];
        int gID     = t_GIDs[i];
        int eID     = t_EIDs[i];
        for (int j = 0; j < quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (quad_points[j] + 1.0) / 2.0 + a;
            t_x[j]    = PtsX + t_quad[j] * direction[0];
            t_y[j]    = PtsY + t_quad[j] * direction[1];
            t_z[j]    = PtsZ + t_quad[j] * direction[2];
        }
        // Evaluate filter  at these points.
        if (b_symMesh)
        {
            m_siacFilterPtrs[m_SymID]->EvaluateFilter(t_quad, t_quad_vals,
                                                      m_meshSpacing);
        }
        else
        {
            m_siacFilterPtrs[m_OneID]->EvaluateFilter(
                t_quad, t_quad_vals, m_meshSpacing, meshTShift, true);
        }
        m_meshHandlePtr->EvaluateAt(t_x, t_y, t_z, gID, eID, t_xyz_vals, 0);
        NekDouble integral = 0;
        // integral
        for (int k = 0; k < quad_npoints; k++)
        {
            integral += 0.5 * quad_weights[k] * t_quad_vals[k] *
                        std::abs((b - a) * direction[0]) * t_xyz_vals[k];
        }
        sum += integral;
    }
    valX = sum;
    return true;
}

bool SmoothieSIAC1D::v_EvaluateNonSymAt(const NekDouble PtsX,
                                        const NekDouble PtsY,
                                        const NekDouble PtsZ, NekDouble &valX,
                                        NekDouble &valY, NekDouble &valZ)

{
    boost::ignore_unused(valY, valZ); // reserved for derivative filter.
    // call HNM
    Array<OneD, NekDouble> direction(3, 0.0), coord(3, 0.0);
    coord[0] = PtsX;
    coord[1] = PtsY;
    coord[2] = PtsZ;
    // TO do . Direction needs to picked up from meshptr. ???
    direction[0] = 1.0;
    NekDouble tmin, tmax;
    vector<NekDouble> HvalX, HvalY, HvalZ, HvalT, SvalT, TvalT;

    NekDouble nonSymShift = 0.0;
    Array<OneD, NekDouble> knotVec(3 * (m_order - 1) + 2, 0.0);
    vector<NekDouble> coords1D;
    m_meshHandlePtr->Get1DVec(coords1D);
    bool b_symMesh = CalculateKnotVec(coord[0], coords1D, knotVec, nonSymShift);

    if (b_symMesh)
    {
        SvalT.resize(knotVec.num_elements());
        for (int i = 0; i < knotVec.num_elements(); i++)
        {
            SvalT[i] = knotVec[i];
        }
        tmin = knotVec[0];
        tmax = knotVec[knotVec.num_elements() - 1];
    }
    else
    {
        valX = -1;
        return false;
    }

    m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ, direction, tmin, tmax, HvalX,
                                 HvalY, HvalZ, HvalT);
    mergeBreakPts(HvalT, SvalT, TvalT);

    // Loop forsumation.
    NekDouble sum = 0.0;
    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.

    int quad_npoints =
        ceil((m_order + 3) / 2) +
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetTotPoints() + 2;
    LibUtilities::PointsKey quadPointsKey(
        quad_npoints, Nektar::LibUtilities::eGaussGaussLegendre);
    Array<OneD, NekDouble> quad_points =
        LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    Array<OneD, NekDouble> quad_weights =
        LibUtilities::PointsManager()[quadPointsKey]->GetW();

    Array<OneD, NekDouble> t_quad(quad_npoints), t_x(quad_npoints),
        t_y(quad_npoints), t_z(quad_npoints), t_quad_vals(quad_npoints),
        t_xyz_vals(quad_npoints);
    vector<int> t_GIDs, t_EIDs;
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, TvalT, t_GIDs,
                                   t_EIDs);

    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        // Create a seg element
        // Actual quadrature points needed are for order = polyorder+Bspline
        // order.
        NekDouble a = TvalT[i];
        NekDouble b = TvalT[i + 1];
        int gID     = t_GIDs[i];
        int eID     = t_EIDs[i];
        for (int j = 0; j < quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (quad_points[j] + 1.0) / 2.0 + a;
            t_x[j]    = PtsX + t_quad[j] * direction[0];
            t_y[j]    = PtsY + t_quad[j] * direction[1];
            t_z[j]    = PtsZ + t_quad[j] * direction[2];
        }
        // Evaluate filter  at these points.
        if (b_symMesh)
        {
            m_siacFilterPtrs[m_SymID]->EvaluateFilterWknots(t_quad, t_quad_vals,
                                                            knotVec, 1.0);
        }
        else
        {
            m_siacFilterPtrs[m_SymID]->EvaluateFilterWknots(t_quad, t_quad_vals,
                                                            knotVec, 1.0);
        }
        m_meshHandlePtr->EvaluateAt(t_x, t_y, t_z, gID, eID, t_xyz_vals, 0);
        NekDouble integral = 0;
        // integral
        for (int k = 0; k < quad_npoints; k++)
        {
            integral += 0.5 * quad_weights[k] * t_quad_vals[k] *
                        std::abs((b - a) * direction[0]) * t_xyz_vals[k];
        }
        sum += integral;
    }
    valX = sum;
    return true;
}

bool SmoothieSIAC1D::v_EvaluateAt_NUK_MetricTensor(
    const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
    NekDouble &valX, NekDouble &valY, NekDouble &valZ,
    Array<OneD, NekDouble> &direction, NekDouble meshSpacing, int varNum)
{
    boost::ignore_unused(valY, valZ); // reserved for derivative filter.
    // call HNM
    Array<OneD, NekDouble> coord(3, 0.0);
    coord[0] = PtsX;
    coord[1] = PtsY;
    coord[2] = PtsZ;
    direction[0] = 1.0;
    NekDouble tmin, tmax;
    vector<NekDouble> HvalX, HvalY, HvalZ, HvalT, SvalT, TvalT;

    Array<OneD, NekDouble> knotVec_BeforeReverse(3 * (m_order - 1) + 2, 0.0);
    Array<OneD, NekDouble> knotVec(3 * (m_order - 1) + 2, 0.0);

    bool b_symMesh = Cal_NUK_ConstMetricTensor(PtsX, PtsY, PtsZ, meshSpacing,
                                               direction, knotVec);
    if (b_symMesh)
    {
        SvalT.resize(knotVec.num_elements());
        for (int i = 0; i < knotVec.num_elements(); i++)
        {
            SvalT[i] = knotVec[i];
        }
        tmin = knotVec[0];
        tmax = knotVec[knotVec.num_elements() - 1];
    }
    else
    {
        valX = -1;
        return false;
    }

    m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ, direction, tmin, tmax, HvalX,
                                 HvalY, HvalZ, HvalT);
    mergeBreakPts(HvalT, SvalT, TvalT);

    // Loop forsumation.
    NekDouble sum = 0.0;
    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.
    int quad_npoints =
        ceil((m_order + 3) / 2) +
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetTotPoints() + 2;
    LibUtilities::PointsKey quadPointsKey(
        quad_npoints, Nektar::LibUtilities::eGaussGaussLegendre);
    Array<OneD, NekDouble> quad_points =
        LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    Array<OneD, NekDouble> quad_weights =
        LibUtilities::PointsManager()[quadPointsKey]->GetW();

    Array<OneD, NekDouble> t_quad(quad_npoints), t_x(quad_npoints),
        t_y(quad_npoints), t_z(quad_npoints), t_quad_vals(quad_npoints),
        t_xyz_vals(quad_npoints);
    vector<int> t_GIDs, t_EIDs;
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, TvalT, t_GIDs,
                                   t_EIDs);

    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        // Create a seg element
        // Actual quadrature points needed are for order = polyorder+Bspline
        // order.
        NekDouble a = TvalT[i];
        NekDouble b = TvalT[i + 1];
        int gID     = t_GIDs[i];
        int eID     = t_EIDs[i];
        for (int j = 0; j < quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (quad_points[j] + 1.0) / 2.0 + a;
            t_x[j]    = PtsX + t_quad[j] * direction[0];
            t_y[j]    = PtsY + t_quad[j] * direction[1];
            t_z[j]    = PtsZ + t_quad[j] * direction[2];
        }
        // Evaluate filter  at these points.
        if (b_symMesh)
        {
            m_siacFilterPtrs[m_SymID]->EvaluateFilterWknots(t_quad, t_quad_vals,
                                                            knotVec, 1.0);
        }
        else
        {
            m_siacFilterPtrs[m_SymID]->EvaluateFilterWknots(t_quad, t_quad_vals,
                                                            knotVec, 1.0);
        }
        m_meshHandlePtr->EvaluateAt(t_x, t_y, t_z, gID, eID, t_xyz_vals,
                                    varNum);
        NekDouble integral = 0;
        // integral
        for (int k = 0; k < quad_npoints; k++)
        {
            integral += 0.5 * quad_weights[k] * t_quad_vals[k] *
                        std::abs((b - a) * direction[0]) * t_xyz_vals[k];
        }
        sum += integral;
    }
    valX = sum;
    return true;
}

bool SmoothieSIAC1D::v_Cal_NUK_ConstMetricTensor(
    const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
    const NekDouble meshSpacing, Array<OneD, NekDouble> &direction,
    Array<OneD, NekDouble> &knotVec)
{
    vector<NekDouble> HvalX, HvalY, HvalZ, HvalT, SvalT;
    NekDouble tmin, tmax;
    m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, tmin, tmax);
    m_siacFilterPtrs[m_SymID]->GetBreakPts(1.0, SvalT);

    // UK_pos
    vector<NekDouble>::iterator UK_pos_iter;
    for (std::vector<NekDouble>::iterator it = SvalT.begin(); it != SvalT.end();
         ++it)
    {
        if (*it > 0.0 + 1e-9)
        {
            UK_pos_iter = it;
            break;
        }
    }

    // Find NUK on the positive side
    m_meshHandlePtr->GetBreakPts_Without_Tmin_Tmax(
        PtsX, PtsY, PtsZ, direction, 0.0, tmax, HvalX, HvalY, HvalZ, HvalT);
    if (HvalT.size() == 0)
    {
        return false;
    }
    ASSERTL0(HvalT.size() > 0, "Make sure more neightbourhood is searched");
    vector<int> t_GIDs, t_EIDs;

    HvalT.insert(HvalT.begin(),
                 0.0); // Adding zero to make sure GetListOfGIDs work.
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, HvalT, t_GIDs,
                                   t_EIDs);
    t_GIDs.pop_back();
    t_EIDs.pop_back();
    HvalT.erase(
        HvalT.begin()); // Removing the first element to maintain consistency.

    // Design a method which return direction for a list of element ratios for
    // the line segments. Use Metric Tensor here.
    vector<NekDouble> scale; // initializing them to 1.0 for test case. (Can do
                             // mesh scaling also.)
    m_meshHandlePtr->GetMTScalingOfGIDs(t_EIDs, direction, scale);

    // Need ratios also.
    vector<NekDouble> tvalues = HvalT;
    vector<NekDouble> NUKpos, NUKneg;

    // Determine zero knot
    if ((m_order % 2) == 0)
    {
        NUKpos.push_back(0.0);
    }

    // Initialize loop
    NekDouble rem_KL, rem_EL, accum_KL, step;
    int elmIndex = 0; // location of the element needs to calculated.
    rem_KL       = *UK_pos_iter;
    accum_KL     = 0.0;
    rem_EL       = tvalues[elmIndex];
    bool b_inMesh;
    int multi = 0;

    while (true)
    {
        step = rem_EL; // Here not req to mul with scale[elmIndex];
        if (step > rem_KL * scale[elmIndex])
        {
            accum_KL += rem_KL * scale[elmIndex];
            rem_EL -= rem_KL * scale[elmIndex];
            NUKpos.push_back(accum_KL);
            UK_pos_iter++;
            if (UK_pos_iter != SvalT.end())
            {
                accum_KL = 0.0;
                rem_KL   = *UK_pos_iter - *(UK_pos_iter - 1);
            }
            else
            {
                break;
            }
        }
        else
        {
            accum_KL += rem_EL;
            rem_KL = rem_KL - (rem_EL) / scale[elmIndex];
            // Check if the next element there or not.
            if (elmIndex + 1 >= tvalues.size())
            {
                b_inMesh = true;
                HvalT.clear();
                while (HvalT.size() == 0 && b_inMesh == true)
                {
                    multi = multi + 1;
                    m_meshHandlePtr->GetBreakPts_Without_Tmin_Tmax(
                        PtsX, PtsY, PtsZ, direction, tmax * multi,
                        tmax * (multi + 1), HvalX, HvalY, HvalZ, HvalT);
                    b_inMesh = m_meshHandlePtr->CanTRangebeAppliedWOMeshShift(
                        PtsX, PtsY, PtsZ, direction, tmax * multi,
                        tmax * (multi + 1));
                }
                if (HvalT.size() == 0 && b_inMesh == false)
                {
                    return false;
                }
                else
                {
                    HvalT.insert(HvalT.begin(),
                                 tmax * multi); // Adding zero to make sure
                                                // GetListOfGIDs work.
                    t_GIDs.clear();
                    t_EIDs.clear();
                    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction,
                                                   HvalT, t_GIDs, t_EIDs);
                    t_GIDs.pop_back();
                    t_EIDs.pop_back();
                    HvalT.erase(HvalT.begin()); // Removing the first element to
                                                // maintain consistency.
                    vector<NekDouble>
                        scale_append; // initializing them to 1.0 for test case.
                                      // (Can do mesh scaling also.)
                    m_meshHandlePtr->GetMTScalingOfGIDs(t_EIDs, direction,
                                                        scale_append);

                    tvalues.insert(tvalues.end(), HvalT.begin(), HvalT.end());
                    scale.insert(scale.end(), scale_append.begin(),
                                 scale_append.end());
                }
            }
            elmIndex++;
            rem_EL = tvalues[elmIndex] - tvalues[elmIndex - 1];
        }
    }

    // Find NUK on the negative side
    m_meshHandlePtr->GetBreakPts_Without_Tmin_Tmax(
        PtsX, PtsY, PtsZ, direction, tmin, 0.0, HvalX, HvalY, HvalZ, HvalT);
    if (HvalT.size() == 0)
    {
        return false;
    }
    ASSERTL0(HvalT.size() > 0, "Make sure more neightbourhood is searched");
    t_GIDs.clear();
    t_EIDs.clear();

    HvalT.push_back(0.0); // Adding zero to make sure GetListOfGIDs work.
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, HvalT, t_GIDs,
                                   t_EIDs);
    t_GIDs.pop_back();
    t_EIDs.pop_back();
    HvalT.pop_back(); // Removing the last element to maintain consistency.

    // Design a method which return direction for a list of element ratios for
    // the line segments. Use Metric Tensor here.
    scale.clear(); // initializing them to 1.0 for test case. (Can do mesh
                   // scaling also.)
    m_meshHandlePtr->GetMTScalingOfGIDs(t_EIDs, direction, scale);

    // Need ratios also.
    tvalues = HvalT;
    std::reverse(std::begin(tvalues), std::end(tvalues));
    std::reverse(std::begin(scale), std::end(scale));

    // Initialize loop
    rem_KL   = 0;
    rem_EL   = 0;
    accum_KL = 0;
    step     = 0;
    elmIndex = 0; // location of the element needs to calculated.
                  // UK_neg_iter
    vector<NekDouble>::reverse_iterator UK_neg_iter;
    for (std::vector<NekDouble>::reverse_iterator it = SvalT.rbegin();
         it != SvalT.rend(); ++it)
    {
        if (*it < 0.0 - 1e-9)
        {
            UK_neg_iter = it;
            break;
        }
    }
    rem_KL   = -1.0 * *UK_neg_iter;
    accum_KL = 0.0;
    rem_EL   = -1 * tvalues[elmIndex]; // tvalues are negative.

    multi = 0;
    while (true)
    {
        step = rem_EL; // Here not req to mul with scale[elmIndex];
        if (step > rem_KL * scale[elmIndex])
        {
            accum_KL += rem_KL * scale[elmIndex];
            rem_EL -= rem_KL * scale[elmIndex];
            NUKneg.push_back(-1.0 * accum_KL);
            UK_neg_iter++;
            if (UK_neg_iter != SvalT.rend())
            {
                accum_KL = 0.0;
                rem_KL   = *(UK_neg_iter - 1) - *(UK_neg_iter);
            }
            else
            {
                break;
            }
        }
        else
        {
            accum_KL += rem_EL;
            rem_KL = rem_KL - (rem_EL) / scale[elmIndex];
            if (elmIndex + 1 >= tvalues.size())
            {
                b_inMesh = true;
                HvalT.clear();
                while (HvalT.size() == 0 && b_inMesh == true)
                {
                    multi = multi + 1;
                    m_meshHandlePtr->GetBreakPts_Without_Tmin_Tmax(
                        PtsX, PtsY, PtsZ, direction, tmin * (multi + 1),
                        tmin * (multi), HvalX, HvalY, HvalZ, HvalT);
                    b_inMesh = m_meshHandlePtr->CanTRangebeAppliedWOMeshShift(
                        PtsX, PtsY, PtsZ, direction, tmin * (multi + 1),
                        tmin * (multi));
                }
                if (HvalT.size() == 0 && b_inMesh == false)
                {
                    // Debug info
                    for (auto it = NUKneg.begin(); it != NUKneg.end(); ++it)
                    {
                    }
                    return false;
                }
                else
                {
                    // HvalT.insert(HvalT.rbegin(),tmin*multi); // Adding zero
                    // to make sure GetListOfGIDs work.
                    HvalT.push_back(
                        tmin *
                        multi); // Adding zero to make sure GetListOfGIDs work.
                    t_GIDs.clear();
                    t_EIDs.clear();
                    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction,
                                                   HvalT, t_GIDs, t_EIDs);
                    t_GIDs.pop_back();
                    t_EIDs.pop_back();
                    HvalT.pop_back(); // erase(HvalT.rbegin()); //Removing the
                                      // first element to maintain consistency.
                    vector<NekDouble>
                        scale_append; // initializing them to 1.0 for test case.
                                      // (Can do mesh scaling also.)
                    m_meshHandlePtr->GetMTScalingOfGIDs(t_EIDs, direction,
                                                        scale_append);

                    tvalues.insert(tvalues.end(), HvalT.rbegin(), HvalT.rend());
                    scale.insert(scale.end(), scale_append.rbegin(),
                                 scale_append.rend());
                }
            }
            elmIndex++;
            rem_EL = tvalues[elmIndex - 1] - tvalues[elmIndex];
        }
    }

    vector<NekDouble> NUK(NUKpos.size());
    vector<NekDouble> NUKNEG(NUKneg.size());
    std::partial_sum(NUKpos.begin(), NUKpos.end(), NUK.begin());
    std::partial_sum(NUKneg.begin(), NUKneg.end(), NUKNEG.begin());
    for (auto it = NUKNEG.begin(); it != NUKNEG.end(); ++it)
    {
        NUK.insert(NUK.begin(), *it);
    }
    for (int i = 0; i < NUK.size(); i++)
    {
        knotVec[i] = NUK[i];
    }
    return true;
}
} // namespace LSIAC
} // namespace Nektar
