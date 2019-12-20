///////////////////////////////////////////////////////////////////////////////
//
// File SymmetricSIAC.h
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
class SymmetricSIAC : public SIACFilter
{

public:
    enum SymFilterType
    {
        BASIC_SIAC_2kp1,
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
public:
    CentralBSplines m_cenBSpline;
    NekDouble m_scaling, m_shift;
    Splines m_splines;

    SymmetricSIAC(int Order);
    SymmetricSIAC(int Order, int nBspl, int nthDerivative = 0);
    SymmetricSIAC(int Order, SymFilterType filterType, int aplha = 0);
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
};
} // namespace LSIAC
} // namespace Nektar
