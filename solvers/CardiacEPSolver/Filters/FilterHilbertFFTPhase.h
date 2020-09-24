///////////////////////////////////////////////////////////////////////////////
//
// File FilterHilbertFFTPhase.h
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
// Description: Outputs phase fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FilterHilbertFFTPhase_H
#define NEKTAR_SOLVERUTILS_FILTERS_FilterHilbertFFTPhase_H

#include <CardiacEPSolver/CellModels/CellModel.h>
#include <SolverUtils/Filters/Filter.h>
#include <fftw3.h>
//using namespace SolverUtils;

namespace Nektar
{
    //namespace SolverUtils
    //{
        class FilterHilbertFFTPhase : public SolverUtils::Filter
        {
        public:
            friend class MemoryManager<FilterHilbertFFTPhase>;

            /// Creates an instance of this class
            static SolverUtils::FilterSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr &pSession,
				const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
                const ParamMap &pParams) {
                SolverUtils::FilterSharedPtr p = MemoryManager<FilterHilbertFFTPhase>::AllocateSharedPtr(pSession, pEquation, pParams);
                //p->InitObject();
                return p;
            }

            ///Name of the class
            static std::string className;

            FilterHilbertFFTPhase(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
				const ParamMap &pParams);
            ~FilterHilbertFFTPhase();


        protected:
            virtual void v_Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
            virtual void v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
            virtual void v_Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
            virtual bool v_IsTimeDependent();

        private:
			//Counters
            unsigned int m_index;
            unsigned int m_outputIndex;		
			unsigned int vCounter;
			
			//Inputs
            unsigned int m_outputFrequency;
			int m_window;
			unsigned int m_overlap;
			bool linear_trans = true;
			bool provided_mean = false;			
            std::string m_outputFile;
            LibUtilities::FieldIOSharedPtr m_fld;
			
			//Section details
			unsigned int npoints;
			NekDouble m_TimeStep;			
			
			//Variables
			double m_mean;
            Array<OneD, NekDouble> oldV; //first dimension is space, second dimension is time
			Array<OneD, NekDouble> oldV_zero_mean;
			Array<OneD, NekDouble> fcoef;
			Array<OneD, NekDouble> h;
			Array<OneD, NekDouble> overlap_phase;
			
			fftw_plan plan_forward;
			fftw_plan plan_backward;
        };
    //}
}

#endif