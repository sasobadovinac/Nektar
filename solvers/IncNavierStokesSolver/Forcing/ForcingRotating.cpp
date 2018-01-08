///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingRotating.cpp
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
// Description: Rotating forcing
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Forcing/ForcingRotating.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

    std::string ForcingRotating::className = GetForcingFactory().
                                RegisterCreatorFunction("Rotating",
                                                        ForcingRotating::create,
                                                        "Forcing Rotating");

    ForcingRotating::ForcingRotating(const LibUtilities::SessionReaderSharedPtr& pSession)
            : Forcing(pSession),
              m_hasAngularVelocity(false)	
    {
    }

    void ForcingRotating::v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
            const unsigned int& pNumForcingFields,
            const TiXmlElement* pForce)
    {
        m_NumVariable = pNumForcingFields;
        int npts       = pFields[0]->GetTotPoints();
        const TiXmlElement* funcNameElmtx;
        const TiXmlElement* funcNameElmty;
        const TiXmlElement* funcNameElmtz;
        funcNameElmtx = pForce->FirstChildElement("OMEGAX");
        funcNameElmty = pForce->FirstChildElement("OMEGAY");
        funcNameElmtz = pForce->FirstChildElement("OMEGAZ");
        ASSERTL0(funcNameElmtx, "Requires OMEGAX tag, specifying the x-rotating angular veloicty.");
        ASSERTL0(funcNameElmty, "Requires OMEGAY tag, specifying the y-rotating angular veloicty.");
        ASSERTL0(funcNameElmtz, "Requires OMEGAZ tag, specifying the z-rotating angular veloicty.");
        if (funcNameElmtx && funcNameElmty && funcNameElmtz)
	{
		m_omegax = 2.0*boost::lexical_cast<NekDouble>(funcNameElmtx->GetText());
		m_omegay = 2.0*boost::lexical_cast<NekDouble>(funcNameElmty->GetText());
		m_omegaz = 2.0*boost::lexical_cast<NekDouble>(funcNameElmtz->GetText());
		if (3==m_NumVariable)
		{
			m_Forcing  = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
			for (int i = 0; i < m_NumVariable; ++i)
                        { 
                        	m_Forcing[i]  = Array<OneD, NekDouble> (npts, 0.0);
			}
		}
		m_hasAngularVelocity =true;
	}
    }

    void ForcingRotating::v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble &time)
    {
        int nq = m_Forcing[0].num_elements();
        if (m_hasAngularVelocity)
	{
		if (2==m_NumVariable) //2D
		{
			// 2*oz*v+out
			Vmath::Svtvp(nq,      m_omegaz, inarray[1], 1, outarray[0],  1, outarray[0], 1);
			//-2*oz*u+out	
			Vmath::Svtvp(nq, -1.0*m_omegaz, inarray[0], 1, outarray[1],  1, outarray[1], 1);	
		}
		else //3D
		{
			// 2*oy*w - 2*oz*v
			Vmath::Svtsvtp(nq, m_omegay,  inarray[2], 1,   m_omegaz, inarray[1], 1, m_Forcing[0], 1);
			Vmath::Svtvp(nq, -1.0, m_Forcing[0], 1, outarray[0], 1, outarray[0],  1);	
			// 2*oz*u - 2*ox*w
			Vmath::Svtsvtp(nq, m_omegaz,  inarray[0], 1,   m_omegax, inarray[2], 1, m_Forcing[1], 1);
			Vmath::Svtvp(nq, -1.0, m_Forcing[1], 1, outarray[1], 1, outarray[1],  1);	
			// 2*ox*v - 2*oy*u
			Vmath::Svtsvtp(nq, m_omegax,  inarray[1], 1,   m_omegay, inarray[0], 1, m_Forcing[2], 1);
			Vmath::Svtvp(nq, -1.0, m_Forcing[2], 1, outarray[2], 1, outarray[2],  1);	
		}
	}
    }
        
}
}
