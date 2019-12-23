///////////////////////////////////////////////////////////////////////////////
//
// File GeneralBSplines.h
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
// Description: GeneralBSplines definition
//
///////////////////////////////////////////////////////////////////////////////
#pragma once
#include "BSplines.h"
#include <iostream>
#include <vector>



namespace Nektar
{
namespace LSIAC
{

/**
 * @brief This class evaluates constructs and evaluates the general BSplines
 */
class GeneralBSplines : public BSplines
{
    // data
public:
    Array<OneD, NekDouble> m_knotVector;
    int m_order;
    // functions
private:
    bool BSplinesBasis(const Array<OneD, NekDouble> &t_pos, const int k,
                       const int j, Array<OneD, NekDouble> &t_val,
                       const NekDouble scaling     = 0.0,
                       const NekDouble meshScaling = 1.0) const;

    bool BSplinesBasis(const Array<OneD, NekDouble> &t_pos,
                       const vector<NekDouble> &kvector, const int k,
                       const int j, Array<OneD, NekDouble> &t_val,
                       const NekDouble scaling     = 0.0,
                       const NekDouble meshScaling = 1.0) const;

protected:
public:
    /**
     *
     */
    GeneralBSplines(const int order);
    GeneralBSplines(const Array<OneD, NekDouble> &knots, const int order);

    bool SetKnotVector(const Array<OneD, NekDouble> &knots);

    bool GetKnotVector(Array<OneD, NekDouble> &knots);

    bool SetOrder(const int order);

    int GetOrder() const;

    bool EvaluateBSplines(const Array<OneD, NekDouble> &t_pos, const int j_th,
                          Array<OneD, NekDouble> &t_values,
                          const NekDouble shift       = 0.0,
                          const NekDouble meshScaling = 1.0) const;

    bool EvaluateBSplines(const Array<OneD, NekDouble> &t_pos,
                          const std::vector<NekDouble> &kvec, const int j_th,
                          Array<OneD, NekDouble> &t_values,
                          const NekDouble shift       = 0.0,
                          const NekDouble meshScaling = 1.0) const;
};
} // namespace LSIAC
} // namespace Nektar
