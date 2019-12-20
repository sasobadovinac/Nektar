///////////////////////////////////////////////////////////////////////////////
//
// File CentralBSplines.cpp
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
// Description: CentralBSplines definition
//
///////////////////////////////////////////////////////////////////////////////
#include "CentralBSplines.h"

namespace Nektar
{
namespace LSIAC
{

CentralBSplines::CentralBSplines(int Order, NekDouble Shift, NekDouble scale)
    : GeneralBSplines(Order), m_shift(Shift), m_scale(scale)
{
    Array<OneD, NekDouble> knots(Order + 1, 0.0);
    for (int i = 0; i < Order + 1; i++)
    {
        knots[i] = -Order / 2.0 + i;
    }
    this->SetKnotVector(knots);
    GeneralBSplines::SetKnotVector(knots);
}

CentralBSplines::CentralBSplines(int Order)
    : GeneralBSplines(Order), m_shift(0.0), m_scale(1.0)
{
    Array<OneD, NekDouble> knots(Order + 1, 0.0);
    for (int i = 0; i < Order + 1; i++)
    {
        knots[i] = -Order / 2.0 + i;
    }
    this->SetKnotVector(knots);
    GeneralBSplines::SetKnotVector(knots);
}

CentralBSplines::CentralBSplines(int Order, NekDouble shift)
    : GeneralBSplines(Order), m_shift(shift), m_scale(1.0)
{
    Array<OneD, NekDouble> knots(Order + 1, 0.0);
    for (int i = 0; i < Order + 1; i++)
    {
        knots[i] = -Order / 2.0 + i;
    }
    this->SetKnotVector(knots);
    GeneralBSplines::SetKnotVector(knots);
}

bool CentralBSplines::EvaluateBSplines(const Array<OneD, NekDouble> &t_pos,
                                       Array<OneD, NekDouble> &t_val,
                                       const NekDouble shift,
                                       const NekDouble meshScaling) const
{
    // For central BSplines. The j is always zero.
    // The spline needs to be shifted and evaluated at x_pos locations.
    GeneralBSplines::EvaluateBSplines(t_pos, 0, t_val, shift, meshScaling);
    return true;
}

} // namespace LSIAC
} // namespace Nektar
