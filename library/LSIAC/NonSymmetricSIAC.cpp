///////////////////////////////////////////////////////////////////////////////
//
// File NonSymmetricSIAC.cpp
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
// Description: NonSymmetricSIAC definition
//
///////////////////////////////////////////////////////////////////////////////
#include "NonSymmetricSIAC.h"
#include <Eigen/Dense>

namespace Nektar
{
namespace LSIAC
{
NonSymmetricSIAC::NonSymmetricSIAC(int order)
    : m_coeffCalculated(false), m_cenBSpline(order), m_splines(order - 1)
{
    m_nthDer = 0;
    // initialize number of coeffecients array using order.
    // If only order is specified. It is by default BASIC_SIAC_4kp1.
    m_filterType = SymFilterType::BASIC_SIAC_2kp1;

    switch (m_filterType)
    {
        case (SymFilterType::BASIC_SIAC_2kp1):
            m_order  = order;
            m_nBSpl  = 2 * (m_order - 1) + 1;
            m_coeffs = Array<OneD, NekDouble>((2 * (order - 1) + 1), 0.0);
            m_genBSplinePtr = std::make_shared<GeneralBSplines>(order);
            break;
        case (SymFilterType::SIAC_GIVEN_R_SPLINES):
            m_order  = order;
            m_nBSpl  = 2 * (m_order - 1) + 1;
            m_coeffs = Array<OneD, NekDouble>((2 * (order - 1) + 1), 0.0);
            m_genBSplinePtr = std::make_shared<GeneralBSplines>(order);
            break;
        default:
            NEKERROR(ErrorUtil : efatal, "Filter not defined");
            break;
    }
}

bool NonSymmetricSIAC::v_EvaluateCoefficients_GivenNumSplines(
    const Array<OneD, NekDouble> &knotVec, const NekDouble kernelShift)
{
    boost::ignore_unused(kernelShift);
    // Changing knot vector to knot Matrix.
    int t_nBSpl = knotVec.num_elements() - m_order;
    vector<vector<NekDouble>> kv_mat;
    for (int nb = 0; nb < t_nBSpl; nb++)
    {
        vector<NekDouble> kv(m_order + 1);
        for (int i = 0; i < m_order + 1; i++)
        {
            kv[i] = knotVec[nb + i];
        }
        kv_mat.push_back(kv);
    }

    switch (m_filterType)
    {
        case (SymFilterType::BASIC_SIAC_2kp1):
            CalCoeffForWideSymKernel(m_order - 1, knotVec, m_coeffs);
            m_splines.Initialize(m_order - 1, t_nBSpl, m_coeffs);
            CalCoeffForKnotMatrixVec_Hanieh(m_order - 1, kv_mat, m_coeffs);
            break;
        default:
            NEKERROR(ErrorUtil : efatal, " Filter type not defined");
    }
    return true;
}

bool NonSymmetricSIAC::v_EvaluateCoefficients(
    const Array<OneD, NekDouble> &knotVec, const NekDouble kernelShift)
{
    boost::ignore_unused(kernelShift);
    // Changing knot vector to knot Matrix.
    vector<vector<NekDouble>> kv_mat;
    for (int nb = 0; nb < m_nBSpl; nb++)
    {
        vector<NekDouble> kv(m_order + 1);
        for (int i = 0; i < m_order + 1; i++)
        {
            kv[i] = knotVec[nb + i];
        }
        kv_mat.push_back(kv);
    }

    switch (m_filterType)
    {
        case (SymFilterType::BASIC_SIAC_2kp1):
            CalCoeffForWideSymKernel(m_order - 1, knotVec, m_coeffs);
            m_splines.Initialize(m_order - 1, m_nBSpl, m_coeffs);
            CalCoeffForKnotMatrixVec_Hanieh(m_order - 1, kv_mat, m_coeffs);
            break;
        default:
            NEKERROR(ErrorUtil
                     : efatal, "Assert or add code for all the other cases");
    }
    return true;
}

bool NonSymmetricSIAC::v_EvaluateCoefficients(const NekDouble kernelShift)
{
    boost::ignore_unused(kernelShift);
    switch (m_filterType)
    {
        case (SymFilterType::BASIC_SIAC_2kp1):
            CalCoeffForWideSymKernel(m_order - 1, m_nBSpl, m_coeffs);
            m_splines.Initialize(m_order - 1, m_nBSpl, m_coeffs);
            break;
        default:
            NEKERROR(ErrorUtil
                     : efatal, "Assert or add code for all the other cases");
    }
    return true;
}

bool NonSymmetricSIAC::EvaluateFilterUsingSplines(
    const Array<OneD, NekDouble> &x_pos, Array<OneD, NekDouble> &t_vals,
    const NekDouble meshScaling, const NekDouble meshShift,
    const bool evalCoeff)
{
    boost::ignore_unused(meshShift, evalCoeff);
    // Always check if coeffecients have already been calculated.
    // For symmetric case. They are always calculated during initilization.
    int nq;
    nq = x_pos.num_elements();
    Vmath::Fill(nq, 0.0, t_vals, 1);
    Array<OneD, NekDouble> t_valTemp(nq);
    for (int i = 0; i < m_nBSpl; i++)
    {
        //	NekDouble it = -m_R*0.5 + i;
        NekDouble it = -(m_nBSpl - 1.0) / 2.0 + i;
        m_cenBSpline.EvaluateBSplines(
            x_pos, t_valTemp, ((NekDouble)it) * meshScaling, meshScaling);
        Vmath::Smul(
            nq, m_coeffs[m_nBSpl - 1 - i] / std::pow(meshScaling, 1 + m_nthDer),
            t_valTemp, 1, t_valTemp, 1);
        Vmath::Vadd(nq, t_valTemp, 1, t_vals, 1, t_vals, 1);
    }
    return true;
}

bool NonSymmetricSIAC::v_EvaluateFilter(const Array<OneD, NekDouble> &x_pos,
                                        Array<OneD, NekDouble> &t_vals,
                                        const NekDouble meshScaling,
                                        const NekDouble meshShift,
                                        const bool evalCoeff)
{
    boost::ignore_unused(meshShift, evalCoeff);
    vector<NekDouble> m_tempVector(10); // Need to fix this.
    vector<NekDouble> v_x_pos(x_pos.num_elements());
    memcpy(&v_x_pos[0], &x_pos[0], x_pos.num_elements() * sizeof(NekDouble));
    vector<NekDouble> v_t_vals(t_vals.num_elements());
    // call the funtion here.
    m_splines.EvaluateUArr(v_x_pos, v_t_vals, m_tempVector, meshScaling, 0.0,
                           m_nthDer);
    memcpy(&t_vals[0], &v_t_vals[0], t_vals.num_elements() * sizeof(NekDouble));
    return true;
}

bool NonSymmetricSIAC::v_EvaluateFilterWknots(
    const Array<OneD, NekDouble> &x_pos, Array<OneD, NekDouble> &t_vals,
    const Array<OneD, NekDouble> &knotVec, const NekDouble meshScaling,
    const NekDouble meshShift, const bool evalCoeff)
{
    boost::ignore_unused(evalCoeff);
    int nq = x_pos.num_elements();
    Vmath::Fill(nq, 0.0, t_vals, 1);
    Array<OneD, NekDouble> t_valTemp(nq);

    // need to create knot vector or matrix
    // Changing knot vector to knot Matrix.
    vector<vector<NekDouble>> kv_mat;
    for (int nb = 0; nb < m_nBSpl; nb++)
    {
        vector<NekDouble> kv(m_order + 1);
        for (int i = 0; i < m_order + 1; i++)
        {
            kv[i] = knotVec[nb + i];
        }
        kv_mat.push_back(kv);
    }
    // Evaluate coefficients
    EvaluateCoefficients(knotVec, meshShift);

    // Evaluate filter
    for (int i = 0; i < m_nBSpl; i++)
    {
        m_genBSplinePtr->EvaluateBSplines(x_pos, kv_mat[i], 0, t_valTemp, 0.0,
                                          meshScaling);
        Vmath::Smul(nq, m_coeffs[i] / meshScaling, t_valTemp, 1, t_valTemp, 1);
        Vmath::Vadd(nq, t_valTemp, 1, t_vals, 1, t_vals, 1);
    }

    return true;
}

bool NonSymmetricSIAC::v_EvaluateFilterWknots_GivenNumSplines(
    const Array<OneD, NekDouble> &x_pos, Array<OneD, NekDouble> &t_vals,
    const Array<OneD, NekDouble> &knotVec, const NekDouble meshScaling,
    const NekDouble meshShift, const bool evalCoeff)
{
    boost::ignore_unused(evalCoeff);
    int nq = x_pos.num_elements();
    Vmath::Fill(nq, 0.0, t_vals, 1);
    Array<OneD, NekDouble> t_valTemp(nq);

    // need to create knot vector or matrix
    // Changing knot vector to knot Matrix.
    int t_nBSpl = knotVec.num_elements() - m_order;
    m_coeffs    = Array<OneD, NekDouble>(t_nBSpl);
    vector<vector<NekDouble>> kv_mat;
    for (int nb = 0; nb < t_nBSpl; nb++)
    {
        vector<NekDouble> kv(m_order + 1);
        for (int i = 0; i < m_order + 1; i++)
        {
            kv[i] = knotVec[nb + i];
        }
        kv_mat.push_back(kv);
    }
    // Evaluate coefficients
    EvaluateCoefficients_GivenNumSplines(knotVec, meshShift);

    // Evaluate filter
    for (int i = 0; i < t_nBSpl; i++)
    {
        m_genBSplinePtr->EvaluateBSplines(x_pos, kv_mat[i], 0, t_valTemp, 0.0,
                                          meshScaling);
        Vmath::Smul(nq, m_coeffs[i] / meshScaling, t_valTemp, 1, t_valTemp, 1);
        Vmath::Vadd(nq, t_valTemp, 1, t_vals, 1, t_vals, 1);
    }

    return true;
}

/// This function return all break point of filter scaled appropriately.
// This function should not be calculating break points everytime.
// They should have been saved at initilization (That poses other problems :( )
bool NonSymmetricSIAC::v_GetBreakPts(const NekDouble scaling,
                                     vector<NekDouble> &valT,
                                     const NekDouble shift)
{
    boost::ignore_unused(shift);
    NekDouble tmin = -1.0, tmax = -1.0;
    int t_nBSpl = 0;
    switch (m_filterType)
    {
        case (SymFilterType::BASIC_SIAC_2kp1):
        case (SymFilterType::CUSTOM_SIAC):
        case (SymFilterType::CUSTOM_SMOOTH_Derivative_SIAC):
        case (SymFilterType::CUSTOM_Derivative_SIAC):
        case (SymFilterType::CUSTOM_SMOOTH_Derivative_SIAC_WOUT_DIVDIFF):
            tmin = -((m_order) / 2.0 + (m_nBSpl - 1.0) / 2.0) * scaling;
            tmax = ((m_order) / 2.0 + (m_nBSpl - 1.0) / 2.0) * scaling;
            break;
        case (SymFilterType::SIAC_GIVEN_R_SPLINES):
            t_nBSpl = 2 * (m_order - 1) + 1;
            tmin    = -((m_order) / 2.0 + (t_nBSpl - 1.0) / 2.0) * scaling;
            tmax    = ((m_order) / 2.0 + (t_nBSpl - 1.0) / 2.0) * scaling;
            break;
        default:
            NEKERROR(ErrorUtil : efatal, "Filter not defined yet");
    }

    valT.clear();
    for (NekDouble t = tmin; t <= tmax; t += scaling)
    {
        valT.push_back(t);
    }
    return true;
}

bool NonSymmetricSIAC::v_GetFilterRange(NekDouble scaling, NekDouble &tmin,
                                        NekDouble &tmax, const NekDouble shift)
{
    boost::ignore_unused(shift);
    int t_nBSpl = 0;
    switch (m_filterType)
    {
        case (SymFilterType::BASIC_SIAC_2kp1):
        case (SymFilterType::CUSTOM_SIAC):
        case (SymFilterType::CUSTOM_SMOOTH_Derivative_SIAC):
        case (SymFilterType::CUSTOM_Derivative_SIAC):
        case (SymFilterType::CUSTOM_SMOOTH_Derivative_SIAC_WOUT_DIVDIFF):
            tmin = -((m_order) / 2.0 + (m_nBSpl - 1.0) / 2.0) * scaling;
            tmax = ((m_order) / 2.0 + (m_nBSpl - 1.0) / 2.0) * scaling;
            break;
        case (SymFilterType::SIAC_GIVEN_R_SPLINES):
            t_nBSpl = 2 * (m_order - 1) + 1;
            tmin    = -((m_order) / 2.0 + (t_nBSpl - 1.0) / 2.0) * scaling;
            tmax    = ((m_order) / 2.0 + (t_nBSpl - 1.0) / 2.0) * scaling;
            break;
        default:
            NEKERROR(ErrorUtil : efatal, "Filter not defined yet");
    }
    return true;
}
} // namespace LSIAC
} // namespace Nektar
