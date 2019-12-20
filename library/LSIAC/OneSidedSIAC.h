///////////////////////////////////////////////////////////////////////////////
//
// File OneSidedSIAC.h
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
// Description: OneSidedSIAC definition
//
///////////////////////////////////////////////////////////////////////////////
#pragma once

#include "CentralBSplines.h"
#include "SIACFilter.h"
#include <vector>

/// This class implements OneSided SIAC at different locations.
/** OneSided SIAC is useful when post processing near the boundaries of the
   mesh. The coefficients need to be re-evaluated everytime a new point is been
   requested. Re-usability of coefficients is possible in 2D and 3D meshes but
   not been implented in this version yet.
*/
namespace Nektar
{
namespace LSIAC
{
class OneSidedSIAC : public SIACFilter
{
public:
    enum OneSidedFilterType
    {
        BASIC_SIAC_2kp1,
        XLi_SIAC_2kp2,
        VAN_SIAC_4kp1,
        Der_SMOOTH_BASIC_SIAC_2kp1,
        Der_BASIC_SIAC_2kp1,
        Der_SMOOTH_BASIC_SIAC_4kp1,
        Der_BASIC_SIAC_4kp1,
        Der_XLi_SIAC_2kp2,
        N_Der_SMOOTH_BASIC_SIAC_2kp1
    };

    std::vector<std::vector<NekDouble>> m_knotMatrix;

private:
    OneSidedFilterType m_filter;
    std::shared_ptr<CentralBSplines> m_cenBSplinePtr;
    // std::vector<<std::vector<NekDouble>> KnotMatrix;
    // std::vector<NekDouble> GenBSPKnotMatrix;
    std::shared_ptr<GeneralBSplines> m_genBSplinePtr;

public:
    OneSidedSIAC(const int Order,
                 OneSidedSIAC::OneSidedFilterType filter =
                     OneSidedFilterType::BASIC_SIAC_2kp1,
                 const int Derivative = 0);
    OneSidedSIAC(int Order, int R, OneSidedSIAC::OneSidedFilterType filter);

protected:
    virtual bool v_GetBreakPts(const NekDouble scaling, vector<NekDouble> &valT,
                               const NekDouble shift = 0.0);
    virtual bool v_GetFilterRange(NekDouble scaling, NekDouble &tmin,
                                  NekDouble &tmax, const NekDouble shift);

    virtual bool v_EvaluateFilter(const Array<OneD, NekDouble> &x_pos,
                                  Array<OneD, NekDouble> &vals,
                                  const NekDouble meshScaling = 1.0,
                                  const NekDouble meshShift   = 0.0,
                                  const bool evalCoeff        = true);

    virtual bool v_EvaluateCoefficients(const NekDouble kernelShift = 0.0);
};
} // namespace LSIAC
} // namespace Nektar
