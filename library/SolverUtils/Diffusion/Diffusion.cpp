///////////////////////////////////////////////////////////////////////////////
//
// File: Diffusion.cpp
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
// Description: Abstract base class for diffusion.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar
{
namespace SolverUtils
{
DiffusionFactory &GetDiffusionFactory()
{
    static DiffusionFactory instance;
    return instance;
}

void Diffusion::InitObject(const LibUtilities::SessionReaderSharedPtr pSession,
                           Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
{
    v_InitObject(pSession, pFields);

    // Div curl storage
    int nPts        = pFields[0]->GetTotPoints();
    m_divVel        = Array<OneD, NekDouble>(nPts, 0.0);
    m_divVelSquare  = Array<OneD, NekDouble>(nPts, 0.0);
    m_curlVelSquare = Array<OneD, NekDouble>(nPts, 0.0);
}

void Diffusion::v_DiffuseCoeffs(
    const std::size_t nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    boost::ignore_unused(nConvectiveFields, fields, inarray, outarray, pFwd,
                         pBwd);
    NEKERROR(ErrorUtil::efatal, "v_DiffuseCoeffs not defined");
}

const Array<OneD, const Array<OneD, NekDouble>> &Diffusion::v_GetTraceNormal()
{
    NEKERROR(ErrorUtil::efatal, "v_GetTraceNormal not defined");
    return NullNekDoubleArrayOfArray;
}

void Diffusion::v_DiffuseCalcDerivative(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    TensorOfArray3D<NekDouble> &qfields,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    boost::ignore_unused(fields, inarray, qfields, pFwd, pBwd);
    NEKERROR(ErrorUtil::efatal, "Not defined for function DiffuseVolumeFLux.");
}

void Diffusion::v_DiffuseVolumeFlux(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    TensorOfArray3D<NekDouble> &qfields, TensorOfArray3D<NekDouble> &VolumeFlux,
    Array<OneD, int> &nonZeroIndex)
{
    boost::ignore_unused(fields, inarray, qfields, VolumeFlux, nonZeroIndex);
    NEKERROR(ErrorUtil::efatal, "Not defined for function DiffuseVolumeFLux.");
}

void Diffusion::v_DiffuseTraceFlux(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    TensorOfArray3D<NekDouble> &qfields, TensorOfArray3D<NekDouble> &VolumeFlux,
    Array<OneD, Array<OneD, NekDouble>> &TraceFlux,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd,
    Array<OneD, int> &nonZeroIndex)
{
    boost::ignore_unused(fields, inarray, qfields, VolumeFlux, TraceFlux, pFwd,
                         pBwd, nonZeroIndex);
    NEKERROR(ErrorUtil::efatal, "Not defined function DiffuseTraceFLux.");
}

} // namespace SolverUtils
} // namespace Nektar
