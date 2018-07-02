///////////////////////////////////////////////////////////////////////////////
//
// File: RigidBody.h
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
// Description: Rigid body motion for FSI problems
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_FSI_RIGIDBODY
#define NEKTAR_SOLVERS_FSI_RIGIDBODY

#include <string>
#include <IncNavierStokesSolver/FSI/Body/FSIBody.h>
#include <SolverUtils/Filters/FilterAeroForces.h>

namespace Nektar
{

class RigidBody: public FSIBody
{
public:

    friend class MemoryManager<RigidBody>;

    /// Creates an instance of this class
    static FSIBodySharedPtr create(
        const LibUtilities::SessionReaderSharedPtr        &pSession,
        const std::weak_ptr<SolverUtils::EquationSystem>   &pEquation,
        const std::map<std::string, std::string>          &pParams)
    {
        FSIBodySharedPtr p =
                MemoryManager<RigidBody>::AllocateSharedPtr(pSession,
                                                            pEquation);
        p->InitObject(pEquation.lock()->UpdateFields(), pParams);
        return p;
    }

    ///Name of the class
    static std::string className;

protected:
    //Mass
    NekDouble                               m_M;
    // Spring coefficients
    Array<OneD, NekDouble>                  m_K;
    // Damping coefficients
    Array<OneD, NekDouble>                  m_C;
    // Number of degrees of freedom
    int                                     m_nDof;
    // Movement directions 
    Array<OneD, Array<OneD, NekDouble> >    m_directions;
    // AeroForces filter
    SolverUtils::FilterAeroForcesSharedPtr  m_filterForces;

    /// Time step
    NekDouble                               m_timestep;
    /// Number of sub steps (i.e. ratio of structure timestep to fluid timestep)
    int                                     m_subSteps;
    /// Time integration order
    int                                     m_intSteps;
    /// Initial time when the body is fixed, to prevent instability in startup
    NekDouble                               m_startTime;

    // Output information
    unsigned int                            m_outputFrequency;
    std::string                             m_outputFile;
    std::ofstream                           m_outputStream;
    bool                                    m_doOutput;
    unsigned int                            m_index;


    // Variables for time integration
    Array<OneD, NekDouble>                  m_displacement;
    Array<OneD, NekDouble>                  m_previousDisp;
    Array<OneD, Array<OneD, NekDouble>>     m_velocity;
    Array<OneD, Array<OneD, NekDouble>>     m_force;


    // Coefficients for Adams time-integration
    static NekDouble AdamsBashforth_coeffs[3][3];
    static NekDouble AdamsMoulton_coeffs[3][3];


    // Constructor
    RigidBody(const LibUtilities::SessionReaderSharedPtr        &pSession,
              const std::weak_ptr<SolverUtils::EquationSystem>  &pEquation);

    // Virtual functions
    virtual void v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pFields,
        const std::map<std::string, std::string>             &pParams);

    virtual void v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pDisplFields,
        const NekDouble &time);

    void GetInitialCondition(
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pFields);

};

}

#endif
