#pragma once

#include "LSIACPostProcessor.h"
#include <iostream>
#include <vector>
using namespace std;

/// This class has different SIAC filters as its subclasses.
/** Currently all the filters have been classified into
        Symmetric and OneSided.
*/
namespace Nektar
{
namespace LSIAC
{
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
    bool GetFilterRange(NekDouble scaling, NekDouble &tmin, NekDouble &tmax,
                        const NekDouble shift = 0.0)
    {
        return v_GetFilterRange(scaling, tmin, tmax, shift);
    }
    bool GetBreakPts(
        const NekDouble scaling, vector<NekDouble> &valT,
        const NekDouble shift = 0.0) // This works for symmetric filter
    {
        return v_GetBreakPts(scaling, valT,
                             shift); // This works for symmetric filter.
    }

    bool EvaluateFilter(const Array<OneD, NekDouble> &x_pos,
                        Array<OneD, NekDouble> &vals,
                        const NekDouble meshScaling = 1.0,
                        const NekDouble meshShift   = 0.0,
                        const bool evalCoeff        = false)
    {
        return v_EvaluateFilter(x_pos, vals, meshScaling, meshShift, evalCoeff);
    }
    bool EvaluateCoefficients(const NekDouble kernelShift = 0.0)
    {
        return v_EvaluateCoefficients(kernelShift);
    }

    bool EvaluateCoefficients(const Array<OneD, NekDouble> &knotVec,
                              const NekDouble kernelShift = 0.0)
    {
        return v_EvaluateCoefficients(knotVec, kernelShift);
    }

    bool EvaluateCoefficients_GivenNumSplines(
        const Array<OneD, NekDouble> &knotVec,
        const NekDouble kernelShift = 0.0)
    {
        return v_EvaluateCoefficients_GivenNumSplines(knotVec, kernelShift);
    }

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
