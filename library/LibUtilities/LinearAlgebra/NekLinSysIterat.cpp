///////////////////////////////////////////////////////////////////////////////
//
// File:  NekLinSysIterat.cpp
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
// Description:  NekLinSysIterat definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekLinSysIterat.h>
#include <LibUtilities/BasicUtils/Timer.h>

using namespace std;

namespace Nektar
{   
    namespace LibUtilities
    {    
        /**
         * @class  NekLinSysIterat
         *
         * Solves a linear system using iterative methods.
         */
        NekLinSysIteratFactory& GetNekLinSysIteratFactory()
        {
            static NekLinSysIteratFactory instance;
            return instance;
        }

        /// Constructor for full direct matrix solve.
        NekLinSysIterat::NekLinSysIterat(
            const LibUtilities::SessionReaderSharedPtr  &pSession,
            const LibUtilities::CommSharedPtr           &vComm,
            const int                                   nDimen)
            : NekNonlinLinSys(pSession,vComm, nDimen)
        {
            std::vector<std::string>  variables(1);
            variables[0] =  pSession->GetVariable(0);
            string variable = variables[0];

            if(pSession->DefinesGlobalSysSolnInfo(variable,
                                              "LinIteratSolverTolerance"))
            {
                m_tolerance = boost::lexical_cast<NekDouble>(
                        pSession->GetGlobalSysSolnInfo(variable,
                                "LinIteratSolverTolerance").c_str());
            }
            else
            {
                pSession->LoadParameter("LinIteratSolverTolerance",
                                        m_tolerance,
                                        NekConstants::kNekIterativeTol);
            }

            if(pSession->DefinesGlobalSysSolnInfo(variable,
                                                  "MaxIterations"))
            {
                m_maxiter = boost::lexical_cast<int>(
                        pSession->GetGlobalSysSolnInfo(variable,
                                "MaxIterations").c_str());
            }
            else
            {
                pSession->LoadParameter("MaxIterations",
                                        m_maxiter,
                                        5000);
            }
        }

        void NekLinSysIterat::v_InitObject()
        {
            NekNonlinLinSys::v_InitObject();
            setUniversalUniqueMap();
        }

        NekLinSysIterat::~NekLinSysIterat()
        {
        }

        void NekLinSysIterat::setUniversalUniqueMap(Array<OneD, int> &map)
        {
            int nmap = map.num_elements();
            if(m_map.num_elements()!=nmap)
            {
                m_map   =   Array<OneD, int>(nmap,0);
            }
            Vmath::Vcopy(nmap,map,1,m_map,1);
        }

        void NekLinSysIterat::setUniversalUniqueMap()
        {
            m_map   =   Array<OneD, int>(m_SysDimen,1);
        }

        void  NekLinSysIterat::Set_Rhs_Magnitude(
            const NekVector<NekDouble> &pIn)
        {
            Array<OneD, NekDouble> vExchange(1, 0.0);
            if (m_map.num_elements() > 0)
            {
                vExchange[0] = Vmath::Dot2(pIn.GetDimension(),
                                        &pIn[0], &pIn[0], &m_map[0]);
            }
            m_Comm->AllReduce(vExchange, LibUtilities::ReduceSum);

            // To ensure that very different rhs values are not being
            // used in subsequent solvers such as the velocit solve in
            // INC NS. If this works we then need to work out a better
            // way to control this.
            NekDouble new_rhs_mag = (vExchange[0] > 1e-6) ? vExchange[0] : 1.0;

            if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
            {
                m_rhs_magnitude = new_rhs_mag;
            }
            else
            {
                m_rhs_magnitude = (m_rhs_mag_sm * (m_rhs_magnitude) +
                                (1.0 - m_rhs_mag_sm) * new_rhs_mag);
            }
        }
    }
}

