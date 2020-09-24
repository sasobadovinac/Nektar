///////////////////////////////////////////////////////////////////////////////
//
// File FilterHilbertFFTPhase.cpp
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

#include <CardiacEPSolver/Filters/FilterHilbertFFTPhase.h>

namespace Nektar
{
    std::string FilterHilbertFFTPhase::className = SolverUtils::GetFilterFactory().RegisterCreatorFunction("HilbertPhase", FilterHilbertFFTPhase::create);

    FilterHilbertFFTPhase::FilterHilbertFFTPhase(
        const LibUtilities::SessionReaderSharedPtr &pSession,
		const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
        const ParamMap &pParams) :
        Filter(pSession, pEquation)
    {
        if (pParams.find("OutputFile") == pParams.end())
        {
            m_outputFile = m_session->GetSessionName();
        }
        else
        {
            ASSERTL0(!(pParams.find("OutputFile")->second.empty()),
                     "Missing parameter 'OutputFile'.");
            m_outputFile = pParams.find("OutputFile")->second;
        }
		
        ASSERTL0(pParams.find("OutputFrequency") != pParams.end(), "Missing parameter 'OutputFrequency'."); //Sampling frequency = per number of timesteps
        LibUtilities::Equation equ(m_session->GetInterpreter(), pParams.find("OutputFrequency")->second);
		m_outputFrequency = floor(equ.Evaluate());
		
		ASSERTL0(pParams.find("WindowSize") != pParams.end(), "Missing parameter 'WindowSize'.");
        LibUtilities::Equation equ1(m_session->GetInterpreter(), pParams.find("WindowSize")->second);
		m_window = floor(equ1.Evaluate()); //No. of output steps taken in each FFT
		ASSERTL0(m_window%2 == 0, "'WindowSize' must be even.");
		
		ASSERTL0(pParams.find("OverlapSize") != pParams.end(), "Missing parameter 'OverlapSize'.");
        LibUtilities::Equation equ2(m_session->GetInterpreter(), pParams.find("OverlapSize")->second);
		m_overlap = floor(equ2.Evaluate());
		ASSERTL0(m_window>m_overlap, "'OverlapSize' must be smaller than 'WindowSize'.");
		
		if (pParams.find("MeanV") != pParams.end())
        {
            provided_mean = true;
			LibUtilities::Equation equ3(m_session->GetInterpreter(), pParams.find("MeanV")->second);
			m_mean = equ3.Evaluate();
        }
		
		if (pParams.find("LinearTransitionOverlap") != pParams.end())
        {
            linear_trans = ( pParams.find("LinearTransitionOverlap")->second == "true" );
        }
				
		m_outputIndex = 0;
        m_index = 0;
		
		m_fld = LibUtilities::FieldIO::CreateDefault(pSession);
    }

    FilterHilbertFFTPhase::~FilterHilbertFFTPhase()
    {

    }

    void FilterHilbertFFTPhase::v_Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
    {
        m_index = 0;
        m_outputIndex = 0;
		vCounter = 0;
		
		m_TimeStep = pFields[0]->GetSession()->GetParameter("TimeStep");
		
		npoints = pFields[0]->GetPhys().size();
		oldV           = Array<OneD, NekDouble> (npoints*m_window);
		oldV_zero_mean = Array<OneD, NekDouble> (npoints*m_window);
		fcoef          = Array<OneD, NekDouble> (npoints*m_window);
		h              = Array<OneD, NekDouble> (npoints*m_window);
		overlap_phase  = Array<OneD, NekDouble> (npoints*m_overlap);		
		std::fill_n(&overlap_phase[0], npoints*m_overlap, 0.0);		
		
		// Initialise FFTW
		fftw_r2r_kind fwd_kind[npoints];
		std::fill_n(fwd_kind, npoints, FFTW_R2HC);		
		plan_forward = fftw_plan_many_r2r(1, 			&m_window,		npoints,
										  &oldV_zero_mean[0], 	&m_window,
										  npoints, 		1,
										  &fcoef[0], 	&m_window,
										  npoints, 		1,
										  &fwd_kind[0], FFTW_ESTIMATE);	
		
		fftw_r2r_kind bwd_kind[npoints];
		std::fill_n(bwd_kind, npoints, FFTW_HC2R);
		plan_backward = fftw_plan_many_r2r(1, 			&m_window,		npoints,
										  &fcoef[0], 	&m_window,
										  npoints, 		1,
										  &h[0], 		&m_window,
										  npoints, 		1,
										  &bwd_kind[0], FFTW_ESTIMATE);	
									
		
		v_Update(pFields, 0.0);
    }

    void FilterHilbertFFTPhase::v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
    {						
		if (m_index++ % m_outputFrequency > 0)
        {
            return;
        }
		
		
		// Store old data
		Array<OneD, NekDouble> currentV = pFields[0]->GetPhys();
		memcpy(&oldV[vCounter*npoints], &currentV[0], npoints*sizeof(NekDouble));
		vCounter++;
		
		// Do Hilbert Transform and output phase
		if (vCounter == m_window) 
		{
			// Zero mean
			if (!provided_mean)
			{
				for (int j = 0; j < npoints; ++j)
				{
					m_mean = 0.0;
					for (int i = 0; i < m_window; ++i)
					{
						m_mean += oldV[i*npoints + j];
					}
					m_mean = m_mean/(double)m_window; //time-average
				
					for (int i = 0; i < m_window; ++i)
					{
						oldV_zero_mean[i*npoints + j] = oldV[i*npoints + j] - m_mean;
					}
				}
			}
			else //minus mean directly
			{
				for (int i = 0; i < npoints*m_window; ++i)
				{
					oldV_zero_mean[i] = oldV[i] - m_mean;
				}
			}
			
			
			// Do Hilbert Transform:
			fftw_execute(plan_forward);
			// Rearrange out: imag->real, -real->imag
			double tmp_fcoef;
			for (int i = 1; i < m_window/2; ++i){
				for (int j = 0; j < npoints; ++j) {					
					tmp_fcoef = -fcoef[i*npoints + j]/m_window;
					fcoef[i*npoints + j] = fcoef[(m_window-i)*npoints + j]/m_window;
					fcoef[(m_window-i)*npoints + j] = tmp_fcoef;
				}
			}
			std::fill_n(&fcoef[0], npoints, 0.0);
			std::fill_n(&fcoef[(m_window/2)*npoints], npoints, 0.0);
			fftw_execute(plan_backward);			
			
			// Output batch		
			Array<OneD, NekDouble> phase_coeff(pFields[0]->GetNcoeffs());
			Array<OneD, NekDouble> phase(npoints);
			
			for (int t = 0; t < m_window-m_overlap; ++t)
			{
				std::stringstream vOutputFilename;
				vOutputFilename << m_outputFile << "_" << m_outputIndex << ".chk";
				
				/* 
				Get FieldDef
				Allocate FieldData
				*/ 
				std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef = pFields[0]->GetFieldDefinitions();		
				std::vector<std::vector<NekDouble>> FieldData(FieldDef.size());
				
				// Calculate phase
				if (t < m_overlap) 
				{
					if (linear_trans)
					{
						for (int i = 0; i < npoints; ++i)
						{
							phase[i] = ((double)(m_overlap-t)/m_overlap)*overlap_phase[t*npoints + i]
							           + ((double)t/m_overlap)*atan2(h[t*npoints + i], oldV_zero_mean[t*npoints + i]);
						}
					}
					else //simple average
					{
						for (int i = 0; i < npoints; ++i)
						{
							phase[i] = 0.5*overlap_phase[t*npoints + i] + 0.5*atan2(h[t*npoints + i], oldV_zero_mean[t*npoints + i]);
						}
					}
				}
				else
				{
					for (int i = 0; i < npoints; ++i)
					{
						phase[i] = atan2(h[t*npoints + i], oldV_zero_mean[t*npoints + i]);
					}
				}
				
				// Convert to coefficient space
				pFields[0]->FwdTrans_IterPerExp(phase, phase_coeff);
				
				/*
				for each item in FieldDef
				push "phase" to FieldDef[i]->m_fields
				call AppendFieldData(FieldDef[i], FieldData[i], data)
				*/
				for(int i = 0; i < FieldDef.size(); ++i)
				{
					FieldDef[i]->m_fields.push_back("phase");
					pFields[0]->AppendFieldData(FieldDef[i], FieldData[i], phase_coeff);
				}		

				// Update time in field info if required
				LibUtilities::FieldMetaDataMap fieldMetaDataMap;
				fieldMetaDataMap["Time"] =  boost::lexical_cast<std::string>(m_outputIndex*m_outputFrequency*m_TimeStep);

				m_fld->Write(vOutputFilename.str(),FieldDef,FieldData,fieldMetaDataMap);
				m_outputIndex++;
			}
			
			// Save phase for next iter
			for (int t = m_window-m_overlap; t < m_window; ++t)
			{
				for (int i = 0; i < npoints; ++i)
				{
					overlap_phase[(t-m_window+m_overlap)*npoints + i] = atan2(h[t*npoints + i], oldV_zero_mean[t*npoints + i]);
				}
			}
			
			// Move overlap part to front
			memcpy(&oldV[0], &oldV[(m_window-m_overlap)*npoints], npoints*m_overlap*sizeof(NekDouble));	
			vCounter = m_overlap;
		}
		
	}


    void FilterHilbertFFTPhase::v_Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
    {

    }

    bool FilterHilbertFFTPhase::v_IsTimeDependent()
    {
        return true;
    }
}
