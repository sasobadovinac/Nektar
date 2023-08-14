///////////////////////////////////////////////////////////////////////////////
//
// File: SchemeInitializer.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: This file isused to add each of the Time Integration Schemes
//              to the NekFactory.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/TimeIntegration/AdamsBashforthTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/AdamsMoultonTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/BDFImplicitTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/CNABTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/DIRKTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/EulerExponentialTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/EulerTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/IMEXGearTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/IMEXdirkTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/MCNABTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/RungeKuttaTimeIntegrationSchemes.h>

#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeFIT.h>

#include <LibUtilities/TimeIntegration/ExplicitTimeIntegrationSchemeSDC.h>
#include <LibUtilities/TimeIntegration/IMEXTimeIntegrationSchemeSDC.h>
#include <LibUtilities/TimeIntegration/ImplicitTimeIntegrationSchemeSDC.h>

#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeGEM.h>

#include <LibUtilities/BasicUtils/SessionReader.h>

namespace Nektar
{
namespace LibUtilities
{

// Register all the schemes with the Time Integration Scheme Factory...
//
#define FACTORYREGISTER(scheme)                                                \
    std::string scheme##TimeIntegrationScheme::className =                     \
        GetTimeIntegrationSchemeFactory().RegisterCreatorFunction(             \
            #scheme, scheme##TimeIntegrationScheme::create)
#define SESSIONREGISTER(scheme)                                                \
    std::string scheme##TimeIntegrationScheme::TimeIntegrationMethodLookupId = \
        SessionReader::RegisterEnumValue("TimeIntegrationMethod", #scheme, 0)

// AdamsBashforthTimeIntegrationSchemes.h
FACTORYREGISTER(AdamsBashforth);
FACTORYREGISTER(AdamsBashforthOrder1);
SESSIONREGISTER(AdamsBashforthOrder1);
FACTORYREGISTER(AdamsBashforthOrder2);
SESSIONREGISTER(AdamsBashforthOrder2);
FACTORYREGISTER(AdamsBashforthOrder3);
SESSIONREGISTER(AdamsBashforthOrder3);
FACTORYREGISTER(AdamsBashforthOrder4);
SESSIONREGISTER(AdamsBashforthOrder4);

// AdamsMoultonTimeIntegrationSchemes.h
FACTORYREGISTER(AdamsMoulton);
FACTORYREGISTER(AdamsMoultonOrder1);
SESSIONREGISTER(AdamsMoultonOrder1);
FACTORYREGISTER(AdamsMoultonOrder2);
SESSIONREGISTER(AdamsMoultonOrder2);
FACTORYREGISTER(AdamsMoultonOrder3);
SESSIONREGISTER(AdamsMoultonOrder3);
FACTORYREGISTER(AdamsMoultonOrder4);
SESSIONREGISTER(AdamsMoultonOrder4);

// BDFImplicitTimeIntegrationSchemes.h
FACTORYREGISTER(BDFImplicit);
FACTORYREGISTER(BDFImplicitOrder1);
SESSIONREGISTER(BDFImplicitOrder1);
FACTORYREGISTER(BDFImplicitOrder2);
SESSIONREGISTER(BDFImplicitOrder2);
FACTORYREGISTER(BDFImplicitOrder3);
SESSIONREGISTER(BDFImplicitOrder3);
FACTORYREGISTER(BDFImplicitOrder4);
SESSIONREGISTER(BDFImplicitOrder4);

// EulerTimeIntegrationSchemes.h
FACTORYREGISTER(Euler);
FACTORYREGISTER(BackwardEuler);
SESSIONREGISTER(BackwardEuler);
FACTORYREGISTER(ForwardEuler);
SESSIONREGISTER(ForwardEuler);

// EulerExponentialTimeIntegrationSchemes.h
FACTORYREGISTER(EulerExponential);

// TimeIntegrationSchemesFIT.h
std::string FractionalInTimeIntegrationScheme::className =
    GetTimeIntegrationSchemeFactory().RegisterCreatorFunction(
        "FractionalInTime", FractionalInTimeIntegrationScheme::create);

// CNABTimeIntegrationScheme.h
FACTORYREGISTER(CNAB);
SESSIONREGISTER(CNAB);

// DIRKTimeIntegrationSchemes.h
FACTORYREGISTER(DIRK);
FACTORYREGISTER(DIRKOrder1);
SESSIONREGISTER(DIRKOrder1);
FACTORYREGISTER(DIRKOrder2);
SESSIONREGISTER(DIRKOrder2);
FACTORYREGISTER(DIRKOrder3);
SESSIONREGISTER(DIRKOrder3);
FACTORYREGISTER(DIRKOrder3_ES5);
SESSIONREGISTER(DIRKOrder3_ES5);
FACTORYREGISTER(DIRKOrder4_ES6);
SESSIONREGISTER(DIRKOrder4_ES6);

// IMEXdirkTimeIntegrationSchemes.h
FACTORYREGISTER(IMEXdirk);
FACTORYREGISTER(IMEXdirk_1_1_1);
SESSIONREGISTER(IMEXdirk_1_1_1);
FACTORYREGISTER(IMEXdirk_1_2_1);
SESSIONREGISTER(IMEXdirk_1_2_1);
FACTORYREGISTER(IMEXdirk_1_2_2);
SESSIONREGISTER(IMEXdirk_1_2_2);
FACTORYREGISTER(IMEXdirk_2_2_2);
SESSIONREGISTER(IMEXdirk_2_2_2);
FACTORYREGISTER(IMEXdirk_2_3_2);
SESSIONREGISTER(IMEXdirk_2_3_2);
FACTORYREGISTER(IMEXdirk_2_3_3);
SESSIONREGISTER(IMEXdirk_2_3_3);
FACTORYREGISTER(IMEXdirk_3_4_3);
SESSIONREGISTER(IMEXdirk_3_4_3);
FACTORYREGISTER(IMEXdirk_4_4_3);
SESSIONREGISTER(IMEXdirk_4_4_3);

// IMEXGearTimeIntegrationScheme.h
FACTORYREGISTER(IMEXGear);
SESSIONREGISTER(IMEXGear);

// IMEXTimeIntegrationSchemes.h
FACTORYREGISTER(IMEX);
FACTORYREGISTER(IMEXOrder1);
SESSIONREGISTER(IMEXOrder1);
FACTORYREGISTER(IMEXOrder2);
SESSIONREGISTER(IMEXOrder2);
FACTORYREGISTER(IMEXOrder3);
SESSIONREGISTER(IMEXOrder3);
FACTORYREGISTER(IMEXOrder4);
SESSIONREGISTER(IMEXOrder4);

// MCNABTimeIntegrationScheme.h
FACTORYREGISTER(MCNAB);
SESSIONREGISTER(MCNAB);

// RungeKuttaTimeIntegrationSchemes.h
FACTORYREGISTER(RungeKutta);
FACTORYREGISTER(RungeKutta1);
SESSIONREGISTER(RungeKutta1);
FACTORYREGISTER(RungeKutta2);
SESSIONREGISTER(RungeKutta2);
FACTORYREGISTER(RungeKutta2_ImprovedEuler);
SESSIONREGISTER(RungeKutta2_ImprovedEuler);
FACTORYREGISTER(RungeKutta2_SSP);
SESSIONREGISTER(RungeKutta2_SSP);
FACTORYREGISTER(RungeKutta3);
SESSIONREGISTER(RungeKutta3);
FACTORYREGISTER(RungeKutta3_SSP);
SESSIONREGISTER(RungeKutta3_SSP);
FACTORYREGISTER(ClassicalRungeKutta4);
SESSIONREGISTER(ClassicalRungeKutta4);
FACTORYREGISTER(RungeKutta4);
SESSIONREGISTER(RungeKutta4);
FACTORYREGISTER(RungeKutta5);
SESSIONREGISTER(RungeKutta5);

// TimeIntegrationSchemesSDC.h
std::string ExplicitTimeIntegrationSchemeSDC::className =
    GetTimeIntegrationSchemeFactory().RegisterCreatorFunction(
        "ExplicitSDC", ExplicitTimeIntegrationSchemeSDC::create);
std::string ImplicitTimeIntegrationSchemeSDC::className =
    GetTimeIntegrationSchemeFactory().RegisterCreatorFunction(
        "ImplicitSDC", ImplicitTimeIntegrationSchemeSDC::create);
std::string IMEXTimeIntegrationSchemeSDC::className =
    GetTimeIntegrationSchemeFactory().RegisterCreatorFunction(
        "IMEXSDC", IMEXTimeIntegrationSchemeSDC::create);

// TimeIntegrationSchemesGEM.h
std::string TimeIntegrationSchemeGEM::className =
    GetTimeIntegrationSchemeFactory().RegisterCreatorFunction(
        "ExtrapolationMethod", TimeIntegrationSchemeGEM::create);

} // end namespace LibUtilities
} // namespace Nektar
