///////////////////////////////////////////////////////////////////////////////
//
// File: FlexibleCylinderBody.h
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

#ifndef NEKTAR_SOLVERS_FSI_FLEXIBLECYLINDERBODY
#define NEKTAR_SOLVERS_FSI_FLEXIBLECYLINDERBODY

#include <string>
#include <IncNavierStokesSolver/FSI/Body/FSIBody.h>
#include <LibUtilities/FFT/NektarFFT.h>
#include <SolverUtils/Filters/FilterAeroForces.h>

namespace Nektar
{

class FlexibleCylinderBody: public FSIBody
{
public:

    friend class MemoryManager<FlexibleCylinderBody>;

    /// Creates an instance of this class
    static FSIBodySharedPtr create(
        const LibUtilities::SessionReaderSharedPtr        &pSession,
        const std::weak_ptr<SolverUtils::EquationSystem>   &pEquation,
        const std::map<std::string, std::string>          &pParams)
    {
        FSIBodySharedPtr p =
                MemoryManager<FlexibleCylinderBody>::AllocateSharedPtr(pSession,
                                                            pEquation);

        p->InitObject(pEquation.lock()->UpdateFields(), pParams);
        return p;
    }

    ///Name of the class
    static std::string className;

protected:
    /// Time step
    NekDouble                               m_timestep;

    // Output information
    unsigned int                            m_outputFrequency;
    std::string                             m_outputFile;
    std::ofstream                           m_outputStream;
    bool                                    m_doOutput;
    unsigned int                            m_index;

    // AeroForces filter
    SolverUtils::FilterAeroForcesSharedPtr  m_filterForces;

    /// Number of planes per processors
    int m_np;
    /// Vibration dimension
    int m_vdim;

    /// Structural dynamics parameters
    NekDouble m_structstiff;
    NekDouble m_cabletension;
    NekDouble m_bendingstiff;
    /// Mass of the cable per unit length
    NekDouble m_structrho;
    /// Damping ratio of the cable
    NekDouble m_structdamp;
    /// Flag marking if using fictitious mass
    bool      m_fictmass;
    /// Fictitious mass
    NekDouble m_fictrho;
    /// Fictitious damping
    NekDouble m_fictdamp;
    /// Length ratio of the cable
    NekDouble m_lhom;
    ///
    LibUtilities::NektarFFTSharedPtr m_FFT;
    /// Support type
    std::string m_supptype;

    /// Storage for the cable's motion(x,y) variables
    Array<OneD, Array<OneD, NekDouble> >m_MotionVars;
    /// Fictitious velocity storage
    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_fV;
    /// Fictitious acceleration storage
    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_fA;
    /// Matrices in Newmart-beta method
    Array<OneD, DNekMatSharedPtr> m_CoeffMat_A;
    /// Matrices in Newmart-beta method
    Array<OneD, DNekMatSharedPtr> m_CoeffMat_B;
    /// [0] is displacements, [1] is velocities, [2] is accelerations
    Array<OneD, std::string> m_funcName;
    /// Motion direction: [0] is 'x' and [1] is 'y'
    Array<OneD, std::string> m_motion;
    /// Determine if the the body motion come from an external file
    Array<OneD, bool>        m_IsFromFile;
    /// Store the derivatives of motion variables in x-direction
    Array<OneD, Array< OneD, NekDouble> > m_zta;
    /// Store the derivatives of motion variables in y-direction
    Array<OneD, Array< OneD, NekDouble> > m_eta;

    // Constructor
    FlexibleCylinderBody(
              const LibUtilities::SessionReaderSharedPtr &pSession,
              const std::weak_ptr<SolverUtils::EquationSystem>   &pEquation);


    // Virtual functions
    virtual void v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pFields,
        const std::map<std::string, std::string>             &pParams);

    virtual void v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pDisplFields,
        const NekDouble &time);

private:
    void CheckIsFromFile(const std::map<std::string, std::string> &pParams);

    void InitialiseCableModel(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);

    void Newmark_betaSolver(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
              Array<OneD, NekDouble> &FcePhysinArray,
              Array<OneD, NekDouble> &MotPhysinArray);

    void EvaluateStructDynModel(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
              Array<OneD, NekDouble> &Hydroforces,
              NekDouble time );

    void SetDynEqCoeffMatrix(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);

    void RollOver(Array<OneD, Array<OneD, NekDouble> > &input);

};

}

#endif
