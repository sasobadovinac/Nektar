///////////////////////////////////////////////////////////////////////////////
//
// File SmoothieSIAC1D.h
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
// Description: SmoothieSIAC1D definition
//
///////////////////////////////////////////////////////////////////////////////
#pragma once
#include "HandleNekMesh.h"
#include "HandleNekMesh1D.h"
#include "SmoothieSIAC.h"
#include "SymmetricSIAC.h"


namespace Nektar
{
namespace LSIAC
{
class SmoothieSIAC1D : public SmoothieSIAC
{
private:
protected:
public:
    SmoothieSIAC1D(const FilterType filter, HandleNekMesh *meshHandle,
                   const int Order, NekDouble meshSpacing = 1.0,
                   const int derivative = 0);

private:
    NekDouble m_meshSpacing;

protected:
    virtual bool v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                              const NekDouble PtsZ, NekDouble &valX,
                              NekDouble &valY, NekDouble &valZ);

    virtual bool v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                              const NekDouble PtsZ, NekDouble &valX,
                              NekDouble &valY, NekDouble &valZ,
                              Array<OneD, NekDouble> &direction,
                              NekDouble meshSpacing = -1.0, int varNum = 0);

    virtual bool v_EvaluateNonSymAt(const NekDouble PtsX, const NekDouble PtsY,
                                    const NekDouble PtsZ, NekDouble &valX,
                                    NekDouble &valY, NekDouble &valZ);

    virtual bool v_EvaluateAt_NUK_MetricTensor(
        const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        Array<OneD, NekDouble> &direction, NekDouble meshSpacing = -1.0,
        int varNum = 0);

    virtual bool v_Cal_NUK_ConstMetricTensor(const NekDouble PtsX,
                                             const NekDouble PtsY,
                                             const NekDouble PtsZ,
                                             const NekDouble meshSpacing,
                                             Array<OneD, NekDouble> &direction,
                                             Array<OneD, NekDouble> &knotVec);
};
} // namespace LSIAC
} // namespace Nektar
