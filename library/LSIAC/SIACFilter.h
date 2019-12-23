///////////////////////////////////////////////////////////////////////////////
//
// File SIACFilter.h
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
// Description: SIACFilter definition
//
///////////////////////////////////////////////////////////////////////////////
#pragma once

#include "LSIACPostProcessor.h"
#include <iostream>
#include <vector>
using namespace std;

namespace Nektar
{
namespace LSIAC
{
/**
 * @brief This class has different SIAC filters as its subclasses.
 * Currently all the filters have been classified into
 *       Symmetric and OneSided.
 */
class SIACFilter : public LSIACPostProcessor
{
public:
    int m_order;
    int m_nBSpl;
    Array<OneD, NekDouble> m_coeffs;
    int m_nthDer;

protected:
    SIACFilter(int order, int nBSpl, int nthDer)
        : m_order(order), m_nBSpl(nBSpl), m_nthDer(nthDer){};
    SIACFilter(){};

public:
    /**
     * @brief Returns the span of the filter after scaling it.
     */
    bool GetFilterRange(NekDouble scaling, NekDouble &tmin, NekDouble &tmax,
                        const NekDouble shift = 0.0)
    {
        return v_GetFilterRange(scaling, tmin, tmax, shift);
    }

    /**
     * @brief Returns the knot positions across the SIAC kernel.
     *
     * The knot positions are scaled and shifted based on the inputs.
     */
    bool GetBreakPts(
        const NekDouble scaling, vector<NekDouble> &valT,
        const NekDouble shift = 0.0) // This works for symmetric filter
    {
        return v_GetBreakPts(scaling, valT,
                             shift); // This works for symmetric filter.
    }

    /**
     * @brief Evaluates the value of the L-SIAC filter at given locations.
     *
     * The input should be in kernel space coordinates.
     * If evalCoeff is true, the coefficients need to be recalculated before
     * evaluating the kernel.
     */
    bool EvaluateFilter(const Array<OneD, NekDouble> &x_pos,
                        Array<OneD, NekDouble> &vals,
                        const NekDouble meshScaling = 1.0,
                        const NekDouble meshShift   = 0.0,
                        const bool evalCoeff        = false)
    {
        return v_EvaluateFilter(x_pos, vals, meshScaling, meshShift, evalCoeff);
    }
    /**
     * @brief This function evaluates the SIAC kernel coefficients.
     */
    bool EvaluateCoefficients(const NekDouble kernelShift = 0.0)
    {
        return v_EvaluateCoefficients(kernelShift);
    }

    /**
     * @brief This function evaluates the SIAC kernel coefficients.
     *
     * The knotVec specifies the location of the knots w.r.t point of
     * evaluation.
     */
    bool EvaluateCoefficients(const Array<OneD, NekDouble> &knotVec,
                              const NekDouble kernelShift = 0.0)
    {
        return v_EvaluateCoefficients(knotVec, kernelShift);
    }

    /**
     * @brief This function evaluates the SIAC kernel coefficients, using the
     * given knots in knotVec.
     */
    bool EvaluateCoefficients_GivenNumSplines(
        const Array<OneD, NekDouble> &knotVec,
        const NekDouble kernelShift = 0.0)
    {
        return v_EvaluateCoefficients_GivenNumSplines(knotVec, kernelShift);
    }

    /**
     * @brief This function evaluates the SIAC kernel coefficients.
     *
     * The coefficients are calculated using the knotvec. The evaluated kernel
     * is evaluated at x_pos and written to vals.
     */
    bool EvaluateFilterWknots(const Array<OneD, NekDouble> &x_pos,
                              Array<OneD, NekDouble> &vals,
                              const Array<OneD, NekDouble> &knotVec,
                              const NekDouble meshScaling = 1.0,
                              const NekDouble meshShift   = 0.0,
                              const bool evalCoeff        = false)
    {
        return v_EvaluateFilterWknots(x_pos, vals, knotVec, meshScaling,
                                      meshShift, evalCoeff);
    }

    /**
     * @brief This function evaluates the SIAC kernel coefficients.
     *
     * The coefficients are calculated using the knotvec. The evaluated kernel
     * is evaluated at x_pos and written to vals.
     */
    bool EvaluateFilterWknots_GivenNumSplines(
        const Array<OneD, NekDouble> &x_pos, Array<OneD, NekDouble> &vals,
        const Array<OneD, NekDouble> &knotVec,
        const NekDouble meshScaling = 1.0, const NekDouble meshShift = 0.0,
        const bool evalCoeff = false)
    {
        return v_EvaluateFilterWknots_GivenNumSplines(
            x_pos, vals, knotVec, meshScaling, meshShift, evalCoeff);
    }

protected:
    virtual bool v_EvaluateFilterWknots(const Array<OneD, NekDouble> &x_pos,
                                        Array<OneD, NekDouble> &vals,
                                        const Array<OneD, NekDouble> &knotVec,
                                        const NekDouble meshScaling,
                                        const NekDouble meshShift,
                                        const bool evalCoeff)
    {
        boost::ignore_unused(x_pos, vals, knotVec, meshScaling, meshShift,
                             evalCoeff);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }
    virtual bool v_EvaluateFilterWknots_GivenNumSplines(
        const Array<OneD, NekDouble> &x_pos, Array<OneD, NekDouble> &vals,
        const Array<OneD, NekDouble> &knotVec,
        const NekDouble meshScaling = 1.0, const NekDouble meshShift = 0.0,
        const bool evalCoeff = false)
    {
        boost::ignore_unused(x_pos, vals, knotVec, meshScaling, meshShift,
                             evalCoeff);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }

    virtual bool v_EvaluateCoefficients(const NekDouble kernelShift = 0.0)
    {
        boost::ignore_unused(kernelShift);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }

    virtual bool v_EvaluateCoefficients(const Array<OneD, NekDouble> &knotVec,
                                        const NekDouble kernelShift = 0.0)
    {
        boost::ignore_unused(knotVec, kernelShift);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }

    virtual bool v_EvaluateCoefficients_GivenNumSplines(
        const Array<OneD, NekDouble> &knotVec,
        const NekDouble kernelShift = 0.0)
    {
        boost::ignore_unused(knotVec, kernelShift);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }

    virtual bool v_GetBreakPts(const NekDouble scaling, vector<NekDouble> &valT,
                               const NekDouble shift = 0.0)
    {
        boost::ignore_unused(scaling, valT, shift);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    }; // This works for symmetric filter.
    virtual bool v_GetFilterRange(NekDouble scaling, NekDouble &tmin,
                                  NekDouble &tmax, const NekDouble shift = 0.0)
    {
        boost::ignore_unused(scaling, tmin, tmax, shift);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    };
    virtual bool v_EvaluateFilter(const Array<OneD, NekDouble> &x_pos,
                                  Array<OneD, NekDouble> &vals,
                                  const NekDouble meshScaling = 1.0,
                                  const NekDouble meshShift   = 0.0,
                                  const bool evalCoeff        = false)
    {
        boost::ignore_unused(x_pos, vals, meshScaling, meshShift, evalCoeff);
        NEKERROR(ErrorUtil::efatal, "Not yet coded");
        return false;
    };

public:
    void CalCoeffForKnotMatrixVec_Hanieh(
        int deg, const std::vector<std::vector<NekDouble>> &kMatrix,
        Nektar::Array<OneD, NekDouble> &coeffs);

protected:
    int factorial(int jj);
    void CalCoeffForStandardKernel(int deg,
                                   Nektar::Array<OneD, NekDouble> &coeffs,
                                   const NekDouble shift = 0.0);
    NekDouble dividedDiff(int ii, int jj,
                          const Nektar::Array<OneD, NekDouble> &tknots,
                          NekDouble x, NekDouble D);
    NekDouble dividedDiff(int ii, int jj, const std::vector<NekDouble> &tknots,
                          NekDouble x, NekDouble D);
    NekDouble Monomial(NekDouble t, NekDouble x, int der, int degree);
    void CalCoeffForWideSymKernel(int deg,
                                  const Nektar::Array<OneD, NekDouble> &tknots,
                                  Nektar::Array<OneD, NekDouble> &coeffs,
                                  const NekDouble shift = 0.0);
    void CalCoeffForWideSymKernel(const int deg, const int nSpl,
                                  Nektar::Array<OneD, NekDouble> &coeffs,
                                  const NekDouble shift = 0.0);
    void CalCoeffForCenBSplDerivatives(const int degree, const int alpha,
                                       int &nBSpl,
                                       Array<OneD, NekDouble> &coeffs);
    void CalCoeffForWideSymKernel(const int deg, const int der, const int nSpl,
                                  Nektar::Array<OneD, NekDouble> &coeffs,
                                  const NekDouble shift = 0.0);
    void CalCoeffForKnotMatrixVec(
        int deg, const std::vector<std::vector<NekDouble>> &kMatrix,
        Nektar::Array<OneD, NekDouble> &coeffs);
    void CalDerivativesForKnotMatrixVec_Hanieh(
        const int deg, const int derivative,
        std::vector<std::vector<NekDouble>> &kMatrix,
        Nektar::Array<OneD, NekDouble> &coeffs);
};
} // namespace LSIAC
} // namespace Nektar
