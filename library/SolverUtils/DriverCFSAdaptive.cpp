///////////////////////////////////////////////////////////////////////////////
//
// File DriverCFSAdaptive.cpp
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
// Description: Incompressible Navier Stokes solver
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <SolverUtils/DriverCFSAdaptive.h>

using namespace std;

namespace Nektar
{
    namespace SolverUtils
    {
        string DriverCFSAdaptive::className = GetDriverCFSFactory().RegisterCreatorFunction("CFSAdaptive", DriverCFSAdaptive::create);
        string DriverCFSAdaptive::driverLookupId = LibUtilities::SessionReader::RegisterEnumValue("Driver","CFSAdaptive",0);

        /**
	 *
         */
        DriverCFSAdaptive::DriverCFSAdaptive(
            const LibUtilities::SessionReaderSharedPtr pSession,
            const SpatialDomains::MeshGraphSharedPtr pGraph,
            const SpatialDomains::MeshGraphSharedPtr pHigherOrderGraph)
            : DriverCFS(pSession, pGraph,pHigherOrderGraph)
        {
        }
    
    
        /**
         *
         */
        DriverCFSAdaptive:: ~DriverCFSAdaptive()
        {
        }
    
    
        /**
         *
         */
        void DriverCFSAdaptive::v_InitObject(ostream &out)
        {
            DriverCFS::v_InitObject(out);
        }
    
    
        void DriverCFSAdaptive::v_Execute(ostream &out)
        
        {
            time_t starttime, endtime;
            NekDouble CPUtime;

            m_equ[0]->PrintSummary(out);

            time(&starttime);
            ASSERTL0(m_equ[0]->GetNpoints()==m_equ[1]->GetNpoints(),"Quad Points Not the Same");
            m_equ[0]->DoInitialise();   
            m_equ[0]->DoSolve();
            //Not sure if need initialize m_equ[1] becasue boundary conditions is also included in this function
            //m_equ[1]->DoInitialise();

            time(&endtime);

            m_equ[0]->Output();
        
            if (m_comm->GetRank() == 0)
            {
                CPUtime = difftime(endtime, starttime);
                cout << "-------------------------------------------" << endl;
                cout << "Total Computation Time = " << CPUtime << "s" << endl;
                cout << "-------------------------------------------" << endl;
            }

            int nwidthcolm = 7+6; // the second value determines the number of sigificant digits

            // Evaluate and output computation time and solution accuracy.
            // The specific format of the error output is essential for the
            // regression tests to work.
            // Evaluate L2 Error
            for(int i = 0; i < m_equ[0]->GetNvariables(); ++i)
            {
                Array<OneD, NekDouble> exactsoln(m_equ[0]->GetTotPoints(), 0.0);

                // Evaluate "ExactSolution" function, or zero array
                m_equ[0]->EvaluateExactSolution(i, exactsoln, 
                                                    m_equ[0]->GetFinalTime());

                NekDouble vL2Error   = m_equ[0]->L2Error  (i, exactsoln);
                NekDouble vLinfError = m_equ[0]->LinfError(i, exactsoln);

                if (m_comm->GetRank() == 0)
                {
                    out << "L 2 error (variable " << m_equ[0]->GetVariable(i) 
                        << ") : " ;
                    out <<std::scientific<<std::setw(nwidthcolm)<<std::setprecision(nwidthcolm-8) 
                        << vL2Error << endl;
                    out << "L inf error (variable " << m_equ[0]->GetVariable(i) 
                        << ") : " ;
                    out <<std::scientific<<std::setw(nwidthcolm)<<std::setprecision(nwidthcolm-8) 
                        << vLinfError << endl;
                }
            }
        }
    }
}

