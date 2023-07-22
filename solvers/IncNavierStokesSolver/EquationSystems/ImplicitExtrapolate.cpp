///////////////////////////////////////////////////////////////////////////////
//
// File: ImplicitExtrapolate.cpp
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
// Description: Abstract base class for ImplicitExtrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/ImplicitExtrapolate.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{
/**
 * Registers the class with the Factory.
 */
std::string ImplicitExtrapolate::className =
    GetExtrapolateFactory().RegisterCreatorFunction(
        "Implicit", ImplicitExtrapolate::create, "Implicit");

std::string ImplicitExtrapolate::solverTypeLookupId =
    LibUtilities::SessionReader::RegisterEnumValue("SolverType", "Implicit",
                                                   eImplicit);
ImplicitExtrapolate::ImplicitExtrapolate(
    const LibUtilities::SessionReaderSharedPtr pSession,
    Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
    MultiRegions::ExpListSharedPtr pPressure, const Array<OneD, int> pVel,
    const SolverUtils::AdvectionSharedPtr advObject)
    : WeakPressureExtrapolate(pSession, pFields, pPressure, pVel, advObject)
{
}

ImplicitExtrapolate::~ImplicitExtrapolate()
{
}

/**
 * Analogous to WeakPressureExtrapolate Implementation
 * However, does not use Extrapolation of CurlCurl-term
 * on boundary
 */
void ImplicitExtrapolate::v_EvaluatePressureBCs(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    const Array<OneD, const Array<OneD, NekDouble>> &N, NekDouble kinvis)
{
    m_pressureCalls++;
    if (m_HBCnumber > 0)
    {
        // Calculate just viscous BCs at current level and put in
        // m_pressureHBCs[nlevels-1]
        CalcNeumannPressureBCs(fields, N, kinvis);

        // \int_bnd q  n.u^{n} ds update current normal of field
        // add m_pressureHBCs to gamma_0/Dt * m_acceleration[0]
        AddVelBC();

        // Copy m_pressureHBCs to m_PbndExp
        CopyPressureHBCsToPbndExp();
    }

    // Evaluate High order outflow conditions
    CalcOutflowBCs(fields, kinvis);
}
} // namespace Nektar
