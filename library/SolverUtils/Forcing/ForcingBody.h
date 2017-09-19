///////////////////////////////////////////////////////////////////////////////
//ost::variate_generator<boost::mt19937&,boost::normal_distribution<> >  m_normal(m_rng, boost::normal_distribution<>(0, 1.0) )
//
// File: ForcingBody.h
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
// Description: Body or FIeld forcing
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCINGBODY
#define NEKTAR_SOLVERUTILS_FORCINGBODY

#include <string>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <SolverUtils/Forcing/Forcing.h>

#include <boost/random/mersenne_twister.hpp>  // for mt19937
#include <boost/random/variate_generator.hpp>  // for variate_generator
#include <boost/random/normal_distribution.hpp>

namespace Nektar
{
namespace SolverUtils
{
    class ForcingBody : public Forcing
    {
        public:

            friend class MemoryManager<ForcingBody>;

            /// Creates an instance of this class
            SOLVER_UTILS_EXPORT static ForcingSharedPtr create(
                    const LibUtilities::SessionReaderSharedPtr& pSession,
                    const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                    const unsigned int& pNumForcingFields,
                    const TiXmlElement* pForce)
            {
                ForcingSharedPtr p = MemoryManager<ForcingBody>::
                                                AllocateSharedPtr(pSession);
                p->InitObject(pFields, pNumForcingFields, pForce);
                return p;
            }

            ///Name of the class
            static std::string classNameBody;
            static std::string classNameField;

        protected:
            SOLVER_UTILS_EXPORT virtual void v_InitObject(
                    const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                    const unsigned int& pNumForcingFields,
                    const TiXmlElement* pForce);

            SOLVER_UTILS_EXPORT virtual void v_Apply(
                    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                    const Array<OneD, Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble &time);

        private:
	    int                             m_iter;
	    NekDouble                       m_timestep;
            std::string                     m_funcName;
            bool                            m_hasTimeFcnScaling;
            LibUtilities::EquationSharedPtr m_timeFcnEqn;
            bool                            m_transform;
            bool                            m_bodyforcing;
            bool                            m_fieldforcing;
	    boost::mt19937  		    m_rng;
	    boost::random::normal_distribution<double>  m_dist;

            bool                            m_oubodyforcing;    // mark OU process
	    NekDouble                       m_tau;              // relaxation time
	    NekDouble                       m_diff;             // diffusion constant
	    NekDouble                       m_initx;            // initial position 
	    NekDouble                       m_xi;               // random variable
	    int                             m_checkou;          // random variable
	    NekDouble                       m_standardvar;      // standard variation
            ForcingBody(const LibUtilities::SessionReaderSharedPtr& pSession);

            virtual ~ForcingBody(void){};

            void Update(
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const NekDouble &time);
	    void OUProcess(NekDouble&  xi);
    };

}
}

#endif
