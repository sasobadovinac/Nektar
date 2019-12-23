///////////////////////////////////////////////////////////////////////////////
//
// File Splines.h
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
// Description: Splines definition
//
///////////////////////////////////////////////////////////////////////////////
#pragma once
#include "LSIACPostProcessor.h"


namespace Nektar
{
namespace LSIAC
{

/** @This class evaluates Splines.
 *
 *This class implements De Boor's algorithm to calculate splines faster.
 */
class Splines : LSIACPostProcessor
{
public:
    int m_deg;
    vector<NekDouble> m_cpts;
    vector<NekDouble> m_knots;

    Splines();
    Splines(int deg);
    void findks(const NekDouble u, int &k, int &s) const;
    void expandSupport();
    bool validk(int k, int s);

    void Initialize(int deg, int n_Bspl, const Array<OneD, NekDouble> &coeffs);
    void EvaluateUArr(const vector<NekDouble> &uAr, vector<NekDouble> &solAr,
                      vector<NekDouble> &pts_eval, NekDouble meshScale = 1.0,
                      NekDouble meshShift = 0.0, int m_nthDer = 0);
    void EvaluateUA(const NekDouble u, NekDouble &sol,
                    vector<NekDouble> &pts_eval);
    void EvaluateU(const NekDouble u, NekDouble &sol);
};

} // namespace LSIAC
} // namespace Nektar
