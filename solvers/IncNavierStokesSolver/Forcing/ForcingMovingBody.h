///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingBody.h
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
// Description: Moving Body (Wavyness and acceleration)
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCINGMOVINGBODY
#define NEKTAR_SOLVERUTILS_FORCINGMOVINGBODY

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/FFT/NektarFFT.h>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <IncNavierStokesSolver/Filters/FilterMovingBody.h>
#include <GlobalMapping/Mapping.h>
#include <IncNavierStokesSolver/Forcing/ExecuteSharpy.h>

namespace Nektar
{

class ForcingMovingBody : public SolverUtils::Forcing
{
    public:

        friend class MemoryManager<ForcingMovingBody>;

        /// Creates an instance of this class
        static SolverUtils::ForcingSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const unsigned int& pNumForcingFields,
                const TiXmlElement* pForce)
        {
            SolverUtils::ForcingSharedPtr p =
                                    MemoryManager<ForcingMovingBody>::
                                            AllocateSharedPtr(pSession);
            p->InitObject(pFields, pNumForcingFields, pForce);
            return p;
        }

        ///Name of the class
        static std::string className;

    protected:
        // Mapping object
        GlobalMapping::MappingSharedPtr               m_mapping;

        virtual void v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
            const unsigned int&                         pNumForcingFields,
            const TiXmlElement*                         pForce);

        virtual void v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr>& fields,
            const Array<OneD, Array<OneD, NekDouble> >& inarray,
                  Array<OneD, Array<OneD, NekDouble> >& outarray,
            const NekDouble&                            time);

    private:

        ForcingMovingBody(
            const LibUtilities::SessionReaderSharedPtr& pSession);

        void CheckIsFromFile(const TiXmlElement* pForce);

        void InitialiseVibrationModel(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);

        void InitialiseFilter(
			const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const TiXmlElement* pForce);

        void DFT(
            const Array<OneD, NekDouble> &Hydroforces,
                  Array<OneD, Array<OneD, NekDouble> > &motions);

        void EvaluateVibrationModel(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble                                   &time);

        void SetModeMatrix();

        void RollOver(Array<OneD, Array<OneD, NekDouble> > &input);

        void AverageForceCoefficients(const Array<OneD, Array <OneD, NekDouble> >& fces,
            Array<OneD, Array<OneD, NekDouble> >& AverageHydFCoeffs);
        
        int m_movingBodyCalls;     ///< number of times the movbody have been called
        int m_np;                  ///< number of planes per processors
        int m_vdim;                ///< vibration dimension
        int m_nstrips;             ///< number of strips
        int m_nv;                  ///< number of vibration modes
        int m_nz;                  ///< number of homo modes

        NekDouble m_structrho;     ///< mass of the cable per unit length
        NekDouble m_structdamp;    ///< damping ratio of the cable
        NekDouble m_wave_number;   ///< damping ratio of the cable
        NekDouble m_omega;		   ///< damping ratio of the cable
        NekDouble m_A;    		   ///< damping ratio of the cable
        NekDouble m_length;        ///< length ratio of the cable
        NekDouble m_timestep;      ///< time step
        ///
        LibUtilities::NektarFFTSharedPtr m_FFT;
        ///
        Nektar::FilterMovingBodySharedPtr m_MovBodyfilter;
        /// storage for the cable's motion(x,y) variables
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_motions;
        /// matrices in mode decomposition method
        Array<OneD, DNekMatSharedPtr> m_CoeffMat_A;
        /// matrices in mode decomposition method
        Array<OneD, DNekMatSharedPtr> m_CoeffMat_B;
        /// [0] is displacements, [1] is velocities, [2] is accelerations
        Array<OneD, std::string> m_funcName;
        Array<OneD, std::string> m_mapfuncName;
        /// motion direction: [0] is 'x' and [1] is 'y'
        Array<OneD, std::string> m_motion;
        /// do determine if the the body motion come from an extern file
        Array<OneD, bool>        m_IsFromFile;
        Array<OneD, bool>        m_IsMapFromFile;
        /// Store the derivatives of motion variables in x-direction
        Array<OneD, Array< OneD, NekDouble> > m_zta;
        /// Store the derivatives of motion variables in y-direction
        Array<OneD, Array< OneD, NekDouble> > m_eta;

        unsigned int                    m_outputFrequency;
        Array<OneD, std::ofstream>      m_outputStream;
        std::string                     m_outputFile_fce;
        std::string                     m_outputFile_mot;
        bool                            m_IsHomostrip;
        bool                            m_IsFictmass;


        Sharpy::Beam m_beam;
};

}

#endif
