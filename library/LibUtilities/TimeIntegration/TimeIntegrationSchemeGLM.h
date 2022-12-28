///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationSchemeGLM.h
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
// Description: Header file of time integration scheme GLM base class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_SCHEME_GLM
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_SCHEME_GLM

#define LUE LIB_UTILITIES_EXPORT

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSolutionGLM.h>

namespace Nektar
{
namespace LibUtilities
{

/**
 * @brief Base class for GLM time integration schemes.
 */
class TimeIntegrationSchemeGLM : public TimeIntegrationScheme
{
public:
    LUE void InitializeSecondaryData(TimeIntegrationAlgorithmGLM *phase,
                                     NekDouble deltaT) const
    {
        v_InitializeSecondaryData(phase, deltaT);
    }

    // Friend classes
    LUE friend std::ostream &operator<<(std::ostream &os,
                                        const TimeIntegrationSchemeGLM &rhs);
    LUE friend std::ostream &operator<<(
        std::ostream &os, const TimeIntegrationSchemeGLMSharedPtr &rhs);

protected:
    // Values stored by each integration phase.
    LUE virtual std::string v_GetVariant() const override;
    LUE virtual unsigned int v_GetOrder() const override;
    LUE virtual std::vector<NekDouble> v_GetFreeParams() const override;
    LUE virtual TimeIntegrationSchemeType v_GetIntegrationSchemeType()
        const override;
    LUE virtual unsigned int v_GetNumIntegrationPhases() const override;

    // Gets the solution Vector
    virtual const TripleArray &v_GetSolutionVector() const override
    {
        return m_solVector->GetSolutionVector();
    }

    // Sets the solution Vector
    virtual void v_SetSolutionVector(const int Offset,
                                     const DoubleArray &y) override
    {
        m_solVector->SetSolutionVector(Offset, y);
    }

    // The worker methods
    LUE virtual void v_InitializeScheme(
        const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
        const TimeIntegrationSchemeOperators &op) override;

    LUE virtual ConstDoubleArray &v_TimeIntegrate(
        const int timestep, const NekDouble delta_t,
        const TimeIntegrationSchemeOperators &op) override;

    LUE virtual void v_InitializeSecondaryData(
        TimeIntegrationAlgorithmGLM *phase, NekDouble deltaT) const;

    LUE virtual void v_print(std::ostream &os) const override;
    LUE virtual void v_printFull(std::ostream &os) const override;

    // These methods should never be used directly, only used by child classes.
    LUE TimeIntegrationSchemeGLM(std::string variant, unsigned int order,
                                 std::vector<NekDouble> freeParams)
        : TimeIntegrationScheme(variant, order, freeParams)
    {
        boost::ignore_unused(variant, order, freeParams);
    }

    virtual ~TimeIntegrationSchemeGLM()
    {
    }

    TimeIntegrationAlgorithmGLMVector m_integration_phases;

    TimeIntegrationSolutionGLMSharedPtr m_solVector;

}; // end class TimeIntegrationScheme

LUE std::ostream &operator<<(std::ostream &os,
                             const TimeIntegrationSchemeGLM &rhs);
LUE std::ostream &operator<<(std::ostream &os,
                             const TimeIntegrationSchemeGLMSharedPtr &rhs);

} // end of namespace LibUtilities
} // end of namespace Nektar

#endif
