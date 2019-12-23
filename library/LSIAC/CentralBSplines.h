///////////////////////////////////////////////////////////////////////////////
//
// File CentralBSplines.h
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
#pragma once
#include "GeneralBSplines.h"

/// The class evaluates CentralB-Splines at given locations.
/**	This calss automatically calculates knot positions, when using the order
   specified. All the knot positions calcualted will range -(k+1)/2 to k+1)/2
*/
namespace Nektar
{
namespace LSIAC
{
class CentralBSplines : public GeneralBSplines
{

private:
protected:
public:
    NekDouble m_shift;
    NekDouble m_scale;

    CentralBSplines(int Order);
    CentralBSplines(int Order, NekDouble shift);
    CentralBSplines(int Order, NekDouble shift, NekDouble scale);

    bool GetShift(NekDouble &shift) const;
    bool SetShift(NekDouble shift);

    bool GetScale(NekDouble &scale) const;
    bool SetScale(NekDouble scale);

    bool SetOrder(int Order);
    int GetOrder() const;

    bool EvaluateBSplines(const Array<OneD, NekDouble> &t_pos,
                          Array<OneD, NekDouble> &t_vals,
                          const NekDouble shift       = 0.0,
                          const NekDouble meshScaling = 1.0) const;
};
} // namespace LSIAC
} // namespace Nektar
