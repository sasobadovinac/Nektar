///////////////////////////////////////////////////////////////////////////////
//
// File SmoothieSIAC3D.cpp
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
// Description: SmoothieSIAC3D definition
//
///////////////////////////////////////////////////////////////////////////////
#include "SmoothieSIAC3D.h"
#include "OneSidedSIAC.h"
#include <LibUtilities/Foundations/ManagerAccess.h> // for Points Manager, etc
#include <boost/timer.hpp>

namespace Nektar
{
namespace LSIAC
{
SmoothieSIAC3D::SmoothieSIAC3D(const FilterType filter,
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
            NEKERROR(ErrorUtil
                     : efatal, "symmetric 4k+1 filter is some how screwed up.");
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
    }
}

bool SmoothieSIAC3D::v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                                  const NekDouble PtsZ, NekDouble &valX,
                                  NekDouble &valY, NekDouble &valZ)
{
    Array<OneD, NekDouble> direction(3, 0.0);
    direction[0] = 1.0;
    return EvaluateAt(PtsX, PtsY, PtsZ, valX, valY, valZ, direction,
                      m_meshSpacing);
}

bool SmoothieSIAC3D::v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                                  const NekDouble PtsZ, NekDouble &valX,
                                  NekDouble &valY, NekDouble &valZ,
                                  Array<OneD, NekDouble> &direction,
                                  NekDouble meshSpacing, int varNum)
{
    boost::ignore_unused(valY, valZ);
    if (meshSpacing < 0)
    { // No parameter was specified.
        meshSpacing = m_meshSpacing;
    }
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
    // Made an assumtion that a symmetric filter can be applied ???
    // return List of mesh breakpoints in that range.
    if (b_symMesh)
    {
        m_siacFilterPtrs[m_SymID]->GetBreakPts(meshSpacing, SvalT);
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

    // Loop forsumation.
    NekDouble sum = 0.0;
    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.
    if (m_calculateQuadrature)
    {
        int quadratureOfMesh =
            m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(0) +
            m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(1) +
            m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(2);

        m_quad_npoints = ceil((m_order + quadratureOfMesh + 1) / 2);

        LibUtilities::PointsKey quadPointsKey(
            m_quad_npoints, LibUtilities::PointsType::eGaussGaussLegendre);

        m_quad_points  = LibUtilities::PointsManager()[quadPointsKey]->GetZ();
        m_quad_weights = LibUtilities::PointsManager()[quadPointsKey]->GetW();
        m_calculateQuadrature = false;
    }

    Array<OneD, NekDouble> t_quad(m_quad_npoints), t_x(m_quad_npoints),
        t_y(m_quad_npoints), t_z(m_quad_npoints), t_quad_vals(m_quad_npoints),
        t_xyz_vals(m_quad_npoints);
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
        for (int j = 0; j < m_quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (m_quad_points[j] + 1.0) / 2.0 + a;
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
        for (int k = 0; k < m_quad_npoints; k++)
        {
            integral += 0.5 * m_quad_weights[k] * t_quad_vals[k] *
                        std::abs((b - a)) * t_xyz_vals[k];
        }
        sum += integral;
    }
    valX = sum;
    return true;
}

bool SmoothieSIAC3D::v_EvaluateRecursiveAt(
    const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
    NekDouble &valX, NekDouble &valY, NekDouble &valZ,
    const vector<std::shared_ptr<SmoothieSIAC>> &Sms,
    vector<Array<OneD, NekDouble>> directions,
    const vector<NekDouble> &meshSpacings, const vector<int> &varNums,
    const int curLevel)
{
    int totLevels = directions.size();
    ASSERTL0((totLevels > curLevel) && (curLevel >= 0) &&
             "Some parameters are not right.");
    NekDouble meshSpacing            = meshSpacings[curLevel];
    Array<OneD, NekDouble> direction = directions[curLevel];
    if (totLevels - 1 == curLevel)
    {
        return Sms[curLevel]->EvaluateAt(
            PtsX, PtsY, PtsZ, valX, valY, valZ, directions[curLevel],
            meshSpacings[curLevel], varNums[curLevel]);
    }

    if (meshSpacing < 0)
    { // No parameter was specified.
        meshSpacing = m_meshSpacing;
    }
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
    // Made an assumtion that a symmetric filter can be applied ???
    // return List of mesh breakpoints in that range.
    if (b_symMesh)
    {
        m_siacFilterPtrs[m_SymID]->GetBreakPts(meshSpacing, SvalT);
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

    // Loop forsumation.
    NekDouble sum = 0.0;
    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.
    if (m_calculateQuadrature)
    {
        int quadratureOfMesh =
            m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(0) +
            m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(1) +
            m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(2);

        m_quad_npoints = ceil((m_order + quadratureOfMesh + 1.0) / 2.0);

        LibUtilities::PointsKey quadPointsKey(
            m_quad_npoints, LibUtilities::PointsType::eGaussGaussLegendre);

        m_quad_points  = LibUtilities::PointsManager()[quadPointsKey]->GetZ();
        m_quad_weights = LibUtilities::PointsManager()[quadPointsKey]->GetW();
        m_calculateQuadrature = false;
    }

    Array<OneD, NekDouble> t_quad(m_quad_npoints), t_x(m_quad_npoints),
        t_y(m_quad_npoints), t_z(m_quad_npoints), t_quad_vals(m_quad_npoints),
        t_xyz_vals(m_quad_npoints);
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
        for (int j = 0; j < m_quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (m_quad_points[j] + 1.0) / 2.0 + a;
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

        for (int j = 0; j < m_quad_npoints; j++)
        {
            Sms[curLevel]->EvaluateRecursiveAt(
                t_x[j], t_y[j], t_z[j], t_xyz_vals[j], valY, valZ, Sms,
                directions, meshSpacings, varNums, curLevel + 1);
        }

        NekDouble integral = 0;
        // integral
        for (int k = 0; k < m_quad_npoints; k++)
        {
            integral += 0.5 * m_quad_weights[k] * t_quad_vals[k] *
                        std::abs((b - a)) * t_xyz_vals[k];
        }
        sum += integral;
    }
    valX = sum;
    return true;
}

} // namespace LSIAC
} // namespace Nektar
