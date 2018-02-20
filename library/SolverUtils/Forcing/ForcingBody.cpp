///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingBody.cpp
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
// Description: Body forcing
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Forcing/ForcingBody.h>
#include <MultiRegions/ExpList.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

    std::string ForcingBody::classNameBody = GetForcingFactory().
        RegisterCreatorFunction("Body",
                                ForcingBody::create,
                                "Body Forcing");
    std::string ForcingBody::classNameField = GetForcingFactory().
        RegisterCreatorFunction("Field",
                                ForcingBody::create,
                                "Field Forcing");

    ForcingBody::ForcingBody(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : Forcing(pSession),
          m_hasTimeFcnScaling(false)
    {
    }

    void ForcingBody::v_InitObject(
                                   const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                                   const unsigned int& pNumForcingFields,
                                   const TiXmlElement* pForce)
    {
        m_NumVariable   = pNumForcingFields;
	m_bodyforcing   = false;
	m_fieldforcing  = false;
	m_oubodyforcing = false;
        const TiXmlElement* funcNameElmt = pForce->FirstChildElement("BODYFORCE");

	const TiXmlElement* ouparameters = NULL; 

	m_iter = 0;
	// get time step size
	m_session->LoadParameter("TimeStep", m_timestep,   0.01);
	m_dist = boost::random::normal_distribution<double>(0.0,1.0);


	if(funcNameElmt) m_bodyforcing = true;

        if(!m_bodyforcing && !m_oubodyforcing)
	{
		funcNameElmt = pForce->FirstChildElement("FIELDFORCE");
		if(funcNameElmt) m_fieldforcing = true;
	}
	if(!m_bodyforcing && !m_fieldforcing)
	{
		funcNameElmt = pForce->FirstChildElement("OUBODYFORCE");
		if(funcNameElmt) m_oubodyforcing = true;
	}
        if(!m_bodyforcing && !m_fieldforcing && !m_oubodyforcing)
        {
            ASSERTL0(funcNameElmt, "Requires BODYFORCE or FIELDFORCE or OUBODYFORCE tag "
                     "specifying function name which prescribes body force.");
        }

	// Get parameters for OU process
	if (m_oubodyforcing)
	{
		ouparameters = pForce->FirstChildElement("OUTAU");
		if (ouparameters)
		{	
			m_tau = boost::lexical_cast<NekDouble>(ouparameters->GetText());
		}
		else
		{
			ASSERTL0(ouparameters, "Requires OUTAU tag, specifying the relaxation parameter.");
		}
		ouparameters = pForce->FirstChildElement("OUDIFF");
		if (ouparameters)
		{
			m_diff = boost::lexical_cast<NekDouble>(ouparameters->GetText());
		}
		else
		{
			ASSERTL0(ouparameters, "Requires OUDIFF tag, specifying the diffusion parameter.");
		}
		ouparameters = pForce->FirstChildElement("OUINITX");
		if (ouparameters)
		{
			m_initx = boost::lexical_cast<NekDouble>(ouparameters->GetText());   // default value of the initial X
		}
		else
		{
			m_initx = 0.0;
		}
		ouparameters = pForce->FirstChildElement("CHECKOU");
		if (ouparameters)
		{
			m_checkou  = boost::lexical_cast<int>(ouparameters->GetText());    // output Xi for each m_checkou steps for restart
		}
		else
		{
			m_checkou  = 100;
		}
		m_xi          = m_initx;
		m_standardvar = sqrt(m_diff*m_tau/2.0); 
	}


        m_funcName = funcNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(m_funcName),
                 "Function '" + m_funcName + "' not defined.");

        bool singleMode, halfMode;
        m_session->MatchSolverInfo("ModeType","SingleMode",singleMode,false);
        m_session->MatchSolverInfo("ModeType","HalfMode",  halfMode,  false);
        bool homogeneous = pFields[0]->GetExpType() == MultiRegions::e3DH1D ||
                           pFields[0]->GetExpType() == MultiRegions::e3DH2D;
        m_transform = (singleMode || halfMode || homogeneous);
        // Time function is optional
        if (m_bodyforcing)   funcNameElmt = pForce->FirstChildElement("BODYFORCETIMEFCN");
        if (m_fieldforcing)  funcNameElmt = pForce->FirstChildElement("FIELDFORCETIMEFCN");
        if (m_oubodyforcing) funcNameElmt = pForce->FirstChildElement("OUBODYFORCETIMEFCN");
	

        // Load time function if specified
        if(funcNameElmt)
        {
            std::string funcNameTime = funcNameElmt->GetText();

            ASSERTL0(!funcNameTime.empty(),
                     "Expression must be given in BODYFORCETIMEFCN or "
                     "FIELDFORCETIMEFCN or OUFORCETIMEFCN");

            m_session->SubstituteExpressions(funcNameTime);
            m_timeFcnEqn = MemoryManager<LibUtilities::Equation>
                ::AllocateSharedPtr(m_session->GetExpressionEvaluator(),funcNameTime);

            m_hasTimeFcnScaling = true;
        }

        m_Forcing = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
        for (int i = 0; i < m_NumVariable; ++i)
        {
            m_Forcing[i] = Array<OneD, NekDouble> (pFields[0]->GetTotPoints(), 0.0);
        }


        Update(pFields, 0.0);
    }


    void ForcingBody::Update(
            const Array< OneD, MultiRegions::ExpListSharedPtr > &pFields,
            const NekDouble &time)
    {
        for (int i = 0; i < m_NumVariable; ++i)
        {
            std::string  s_FieldStr   = m_session->GetVariable(i);
            ASSERTL0(m_session->DefinesFunction(m_funcName, s_FieldStr),
                     "Variable '" + s_FieldStr + "' not defined.");
            GetFunction(pFields, m_session, m_funcName, true)->Evaluate(s_FieldStr, m_Forcing[i], time);
        }

        // If singleMode or halfMode, transform the forcing term to be in
        // physical space in the plane, but Fourier space in the homogeneous
        // direction
        if (m_transform)
        {
            for (int i = 0; i < m_NumVariable; ++i)
            {
                pFields[0]->HomogeneousFwdTrans(m_Forcing[i], m_Forcing[i]);
            }
        }
    }

    void ForcingBody::OUProcess(NekDouble& xi)
    {	
	xi      = m_initx*exp(-m_timestep/m_tau)+sqrt((m_diff*m_tau*0.5)*(1.0-exp(-2.0*m_timestep/m_tau)))*m_dist(m_rng);
	m_initx = xi;
    }
    void ForcingBody::v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble &time)
    {
        if(m_hasTimeFcnScaling)
        {
            Array<OneD, NekDouble>  TimeFcn(1);
            for (int i = 0; i < m_NumVariable; i++)
            {
                EvaluateTimeFunction(time, m_timeFcnEqn, TimeFcn);
		if (m_oubodyforcing)
		{
			// modulate TimeFcn
			TimeFcn[0] =  TimeFcn[0]*m_xi/m_standardvar;
		}	
                Vmath::Svtvp(outarray[i].num_elements(), TimeFcn[0],
                             m_Forcing[i], 1,
                             outarray[i],  1,
                             outarray[i],  1);
            }
	    if (m_oubodyforcing)
	    {
		OUProcess(m_xi);
		if(m_session->GetComm()->GetRank()==0)
		{
			if ((m_iter+1)%m_checkou==0)
			{
				cout <<"OU Process normalised X at the step "<< m_iter+1 << "\t : " << m_xi << endl;  
			}
		}
	    }   
        }
        else
        {
            Update(fields, time);

            for (int i = 0; i < m_NumVariable; i++)
            {
                Vmath::Vadd(outarray[i].num_elements(), outarray[i], 1,
                            m_Forcing[i], 1, outarray[i], 1);
            }
        }
	++m_iter;
    }

}
}

// Add the OU process
