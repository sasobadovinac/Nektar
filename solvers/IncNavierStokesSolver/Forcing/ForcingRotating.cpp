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
            m_omegax = boost::lexical_cast<NekDouble>(funcNameElmtx->GetText());
            m_omegay = boost::lexical_cast<NekDouble>(funcNameElmty->GetText());
            m_omegaz = boost::lexical_cast<NekDouble>(funcNameElmtz->GetText());
            if (m_NumVariable == 3)
            {
                m_Forcing  = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
                for (int i = 0; i < m_NumVariable; ++i)
                { 
                    	m_Forcing[i]  = Array<OneD, NekDouble> (npts, 0.0);
                }
		m_tmp = Array<OneD, NekDouble> (npts, 0.0);
            }
	    if ((m_NumVariable == 2) || (m_NumVariable == 3))
	    {
                m_xyz      = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
                for (int i = 0; i < m_NumVariable; ++i)
                { 
                    	m_xyz[i]      = Array<OneD, NekDouble> (npts, 0.0);
                }
		// initialise the coordinates
		if (m_NumVariable == 2)
		{
			pFields[0]->GetCoords( m_xyz[0], m_xyz[1]);
		}
		if (m_NumVariable == 3)
		{	
			pFields[0]->GetCoords( m_xyz[0], m_xyz[1], m_xyz[2]);
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
		/// Coriolis forcing
                //-2*oz*v+out
                Vmath::Svtvp(nq, -2.0*m_omegaz, inarray[1], 1,
                             outarray[0],  1, outarray[0], 1);
                //+2*oz*u+out	
                Vmath::Svtvp(nq,  2.0*m_omegaz, inarray[0], 1,
                             outarray[1],  1, outarray[1], 1);	
		/// Centrifugal forcing
            }
            else //3D - 2 \Omega x u 
            {
                /// Coriolis forcing 
                // 2*oy*w - 2*oz*v
                Vmath::Svtsvtp(nq, 2.0*m_omegay,  inarray[2], 1,
                               -2.0*m_omegaz, inarray[1], 1, m_Forcing[0], 1);
                Vmath::Vsub(nq,outarray[0],1, m_Forcing[0], 1, outarray[0], 1);
                // 2*oz*u - 2*ox*w
                Vmath::Svtsvtp(nq, 2.0*m_omegaz,  inarray[0], 1,
                               -2.0*m_omegax, inarray[2], 1, m_Forcing[1], 1);
                Vmath::Vsub(nq,outarray[1],1, m_Forcing[1], 1, outarray[1], 1);	
                // 2*ox*v - 2*oy*u
                Vmath::Svtsvtp(nq, 2.0*m_omegax,  inarray[1], 1,
                               -2.0*m_omegay, inarray[0], 1, m_Forcing[2], 1);
                Vmath::Vsub(nq,outarray[2],1, m_Forcing[2], 1, outarray[2], 1);	


		/// Centrifugal forcing
		// - x*oy^2 + ox*y*oy - x*oz^2 + ox*z*oz
		Vmath::Svtsvtp(nq, -m_omegay*m_omegay,  m_xyz[0], 1,
                                    m_omegax*m_omegay,  m_xyz[1], 1, m_tmp, 1);
		Vmath::Svtsvtp(nq, -m_omegaz*m_omegaz,  m_xyz[0], 1,
                                    m_omegax*m_omegaz,  m_xyz[2], 1, m_Forcing[0], 1);
		Vmath::Vadd(nq,m_tmp,1, m_Forcing[0], 1, m_Forcing[0], 1);
		Vmath::Vsub(nq,outarray[0],1, m_Forcing[0], 1, outarray[0], 1);
		// - y*ox^2 + oy*x*ox - y*oz^2 + oy*z*oz
		Vmath::Svtsvtp(nq, -m_omegax*m_omegax,  m_xyz[1], 1,
                                    m_omegay*m_omegax,  m_xyz[0], 1, m_tmp, 1);
                Vmath::Svtsvtp(nq, -m_omegaz*m_omegaz,  m_xyz[1], 1,
                                    m_omegay*m_omegaz,  m_xyz[2], 1, m_Forcing[1], 1);
                Vmath::Vadd(nq,m_tmp,1, m_Forcing[1], 1, m_Forcing[1], 1);
                Vmath::Vsub(nq,outarray[1],1, m_Forcing[1], 1, outarray[1], 1);
		// - z*ox^2 + oz*x*ox - z*oy^2 + oz*y*oy
		Vmath::Svtsvtp(nq, -m_omegax*m_omegax,  m_xyz[2], 1,
                                    m_omegaz*m_omegax,  m_xyz[0], 1, m_tmp, 1);
                Vmath::Svtsvtp(nq, -m_omegay*m_omegay,  m_xyz[2], 1,
                                    m_omegaz*m_omegay,  m_xyz[1], 1, m_Forcing[2], 1);
                Vmath::Vadd(nq,m_tmp,1, m_Forcing[2], 1, m_Forcing[2], 1);
                Vmath::Vsub(nq,outarray[2],1, m_Forcing[2], 1, outarray[2], 1);
            }
        }
    }   
}
}
