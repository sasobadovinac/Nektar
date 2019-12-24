///////////////////////////////////////////////////////////////////////////////
//
// File SymmetricSIAC.cpp
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
// Description: SymmetricSIAC definition
//
///////////////////////////////////////////////////////////////////////////////
#include "SymmetricSIAC.h"
#include <Eigen/Dense>

namespace Nektar
{
namespace LSIAC
{

SymmetricSIAC::SymmetricSIAC(int order)
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
            break;
        default:
            ASSERTL0(false, "something is wrong");
            break;
    }
    EvaluateCoefficients();
}

SymmetricSIAC::SymmetricSIAC(int Order, int nBSpl, int nthDerivative)
    : m_cenBSpline(Order), m_splines(Order - 1)
{
    m_nthDer = nthDerivative;
    if (m_nthDer == 0)
    {
        m_filterType = SymFilterType::CUSTOM_SIAC;
    }
    switch (m_filterType)
    {
        case (SymFilterType::CUSTOM_SIAC):
            m_order  = Order;
            m_nBSpl  = nBSpl;
            m_coeffs = Array<OneD, NekDouble>(nBSpl, 0.0);
            break;
        default:
            NEKERROR(ErrorUtil::efatal, "Symmetric SIAC Constructor");
    }
    EvaluateCoefficients();
}

SymmetricSIAC::SymmetricSIAC(int Order, SymFilterType filterType, int nthDer)
    : m_filterType(filterType), m_cenBSpline(Order), m_splines(Order - 1)
{
    m_nthDer = nthDer;
    switch (m_filterType)
    {
        case (SymFilterType::CUSTOM_SMOOTH_Derivative_SIAC):
            m_order  = Order;
            m_nBSpl  = 2 * (m_order - 1) + 1 + nthDer;
            m_coeffs = Array<OneD, NekDouble>(m_nBSpl, 0.0);
            break;
        case (SymFilterType::CUSTOM_Derivative_SIAC):
            m_order  = Order;
            m_nBSpl  = 2 * (m_order - 1) + 1;
            m_coeffs = Array<OneD, NekDouble>(m_nBSpl, 0.0);
            break;
        case (SymFilterType::CUSTOM_SMOOTH_Derivative_SIAC_WOUT_DIVDIFF):
            m_order      = Order + nthDer;
            m_nthDer     = 0;
            m_nBSpl      = 2 * (Order - 1) + 1;
            m_coeffs     = Array<OneD, NekDouble>(m_nBSpl, 0.0);
            m_cenBSpline = CentralBSplines(m_order);
            break;
        default:
            NEKERROR(ErrorUtil::efatal, "Filter type not defined");
    }
    EvaluateCoefficients();
}

bool SymmetricSIAC::v_EvaluateCoefficients(const NekDouble kernelShift)
{
    boost::ignore_unused(kernelShift); // no effect in case of symmetric kernel;
    switch (m_filterType)
    {
        case (SymFilterType::BASIC_SIAC_2kp1):
            CalCoeffForWideSymKernel(m_order - 1, m_nBSpl, m_coeffs);
            m_splines.Initialize(m_order - 1, m_nBSpl, m_coeffs);
            break;
        case (SymFilterType::CUSTOM_SIAC):
            CalCoeffForWideSymKernel(m_order - 1, m_nBSpl, m_coeffs);
            m_splines.Initialize(m_order - 1, m_nBSpl, m_coeffs);
            break;
        case (SymFilterType::CUSTOM_SMOOTH_Derivative_SIAC):
            CalCoeffForWideSymKernel(m_order - 1 + m_nthDer, m_nthDer, m_nBSpl,
                                     m_coeffs);
            CalCoeffForCenBSplDerivatives(m_order - 1 + m_nthDer, m_nthDer,
                                          m_nBSpl, m_coeffs);
            break;
        case (SymFilterType::CUSTOM_Derivative_SIAC):
            CalCoeffForWideSymKernel(m_order - 1 + m_nthDer, m_nthDer, m_nBSpl,
                                     m_coeffs);
            CalCoeffForCenBSplDerivatives(m_order - 1 + m_nthDer, m_nthDer,
                                          m_nBSpl, m_coeffs);
            m_splines.Initialize(m_order - 1, m_nBSpl, m_coeffs);
            break;
        case (SymFilterType::CUSTOM_SMOOTH_Derivative_SIAC_WOUT_DIVDIFF):
            CalCoeffForWideSymKernel(m_order - 1, m_nBSpl, m_coeffs);
            m_splines.Initialize(m_order - 1, m_nBSpl, m_coeffs);
            break;
        default:
            NEKERROR(ErrorUtil::efatal, "New filter type. Not accounted for");
    }

    return true;
}

bool SymmetricSIAC::EvaluateFilterUsingSplines(
    const Array<OneD, NekDouble> &x_pos, Array<OneD, NekDouble> &t_vals,
    const NekDouble meshScaling, const NekDouble meshShift,
    const bool evalCoeff)
{
    boost::ignore_unused(
        meshShift,
        evalCoeff); // this parameter has no meaning in symmetric kernel
    // Always check if coeffecients have already been calculated.
    // For symmetric case. They are always calculated during initilization.
    int nq;
    nq = x_pos.num_elements();
    Vmath::Fill(nq, 0.0, t_vals, 1);
    Array<OneD, NekDouble> t_valTemp(nq);
    for (int i = 0; i < m_nBSpl; i++)
    {
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

bool SymmetricSIAC::v_EvaluateFilter(const Array<OneD, NekDouble> &x_pos,
                                     Array<OneD, NekDouble> &t_vals,
                                     const NekDouble meshScaling,
                                     const NekDouble meshShift,
                                     const bool evalCoeff)
{
    boost::ignore_unused(
        meshShift,
        evalCoeff); // this parameter has no meaning in symmetric kernel
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

/// This function return all break point of filter scaled appropriately.
// This function should not be calculating break points everytime.
// They should have been saved at initilization (That poses other problems :( )
bool SymmetricSIAC::v_GetBreakPts(const NekDouble scaling,
                                  vector<NekDouble> &valT,
                                  const NekDouble shift)
{
    boost::ignore_unused(
        shift); // this parameter has no meaning in symmetric kernel
    NekDouble tmin = -1.0, tmax = -1.0;
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
        default:
            NEKERROR(ErrorUtil::efatal, "Filter out of scope");
    }

    valT.clear();
    for (NekDouble t = tmin; t <= tmax; t += scaling)
    {
        valT.push_back(t);
    }
    return true;
}

bool SymmetricSIAC::v_GetFilterRange(NekDouble scaling, NekDouble &tmin,
                                     NekDouble &tmax, const NekDouble shift)
{
    boost::ignore_unused(
        shift); // this parameter has no meaning in symmetric kernel
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
        default:
            NEKERROR(ErrorUtil::efatal, "Filter out of scope");
    }
    return true;
}
} // namespace LSIAC
} // namespace Nektar
