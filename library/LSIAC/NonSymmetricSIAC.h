#pragma once
#include "CentralBSplines.h"
#include "SIACFilter.h"
#include "Splines.h"
/// This class evaluates symmetric siac filter and saves coeffecients for future
/// use.
/**
        - This class is useful when post processing parts of the mesh where
   Symmetric coeffecients are need multiple times.
        - Not calculating coeffecients and filter knots everytime would save lot
   of computation time.
        - Usage: The user would initiate a class with the order,numberof Bspline
   basis he needs.
              - Since this is a symmetric filter knots can be calculated
   directly from number of Bsplines.
                  - The coefficient calculation for normal filters is different
   from derivative filter.
                  - There are different way to calcualte coeffecients. (Matrix
   slove and direct.). We need to direct method in final version.
        - This class assumes the use of only centralBSplines as basis.
*/

namespace Nektar
{
namespace LSIAC
{
class NonSymmetricSIAC : public SIACFilter
{

public:
    enum SymFilterType
    {
        BASIC_SIAC_2kp1,
        SIAC_GIVEN_R_SPLINES,
        EXTENDED_SIAC_4kp1,
        CUSTOM_SIAC,
        CUSTOM_SMOOTH_Derivative_SIAC,
        CUSTOM_Derivative_SIAC,
        CUSTOM_SMOOTH_Derivative_SIAC_WOUT_DIVDIFF
    };
    bool m_coeffCalculated;
    Array<OneD, NekDouble> m_coeffs;
    SymFilterType m_filterType;

protected:
private:
    std::shared_ptr<GeneralBSplines> m_genBSplinePtr;

public:
    CentralBSplines m_cenBSpline;
    NekDouble m_scaling, m_shift;
    Splines m_splines;

    NonSymmetricSIAC(int Order);

    bool SetCoefficients(const Array<OneD, NekDouble> &coeff);
    bool GetCoefficients(const Array<OneD, NekDouble> &coeff) const;
    bool GetFilterBreakPoints(Array<OneD, NekDouble> &bPts) const;
    //		bool GetFilterRange(const NekDouble scaling, NekDouble &tmin,
    // NekDouble &tmax) const;
    bool EvaluateFilterUsingSplines(const Array<OneD, NekDouble> &x_pos,
                                    Array<OneD, NekDouble> &vals,
                                    const NekDouble meshScaling = 1.0,
                                    const NekDouble meshShift   = 0.0,
                                    const bool evalCoeff        = false);

    // This functions takes the standard break points multiply by this scaling
    // factor and return them.
protected:
    virtual bool v_EvaluateFilterWknots(const Array<OneD, NekDouble> &x_pos,
                                        Array<OneD, NekDouble> &vals,
                                        const Array<OneD, NekDouble> &knotVec,
                                        const NekDouble meshScaling = 1.0,
                                        const NekDouble meshShift   = 0.0,
                                        const bool evalCoeff        = false);

    virtual bool v_EvaluateFilterWknots_GivenNumSplines(
        const Array<OneD, NekDouble> &x_pos, Array<OneD, NekDouble> &vals,
        const Array<OneD, NekDouble> &knotVec,
        const NekDouble meshScaling = 1.0, const NekDouble meshShift = 0.0,
        const bool evalCoeff = false);

    virtual bool v_GetBreakPts(const NekDouble scaling, vector<NekDouble> &valT,
                               const NekDouble shift = 0.0);
    virtual bool v_GetFilterRange(NekDouble scaling, NekDouble &tmin,
                                  NekDouble &tmax, const NekDouble shift = 0.0);

    virtual bool v_EvaluateFilter(const Array<OneD, NekDouble> &x_pos,
                                  Array<OneD, NekDouble> &vals,
                                  const NekDouble meshScaling = 1.0,
                                  const NekDouble meshShift   = 0.0,
                                  const bool evalCoeff        = false);

    //        virtual bool v_EvaluateFilter(const Array<OneD,NekDouble>
    //        &x_pos,Array<OneD,NekDouble> &vals,
    //                   const NekDouble meshScaling=1.0 ) const;
    virtual bool v_EvaluateCoefficients(const NekDouble kernelShift = 0.0);
    virtual bool v_EvaluateCoefficients(const Array<OneD, NekDouble> &knotVec,
                                        const NekDouble kernelShift);
    virtual bool v_EvaluateCoefficients_GivenNumSplines(
        const Array<OneD, NekDouble> &knotVec, const NekDouble kernelShift);
};
} // namespace LSIAC
} // namespace Nektar
