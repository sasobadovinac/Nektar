///////////////////////////////////////////////////////////////////////////////
//
// File DriverMultiLevel.cpp
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

#include <SolverUtils/DriverMultiLevel.h>

using namespace std;

namespace Nektar
{
    namespace SolverUtils
    {
        string DriverMultiLevel::className = GetDriverFactory().RegisterCreatorFunction("MultiLevel", DriverMultiLevel::create);
        string DriverMultiLevel::driverLookupId = LibUtilities::SessionReader::RegisterEnumValue("Driver","MultiLevel",0);

        /**
	 *
         */
        DriverMultiLevel::DriverMultiLevel(
            const LibUtilities::SessionReaderSharedPtr pSession,
            const SpatialDomains::MeshGraphSharedPtr pGraph)
            : Driver(pSession, pGraph)
        {
        }
    
    
        /**
         *
         */
        DriverMultiLevel:: ~DriverMultiLevel()
        {
        }

        void DriverMultiLevel::v_DoMultiOrderProjection(const Array<OneD,const Array<OneD, NekDouble>> &inarray,
                                                           Array<OneD, Array<OneD,NekDouble>> &outarray,
                                                                        const NekDouble time)
        {
            ASSERTL1(m_equ[1],"Need to define EquationSystem[1]");
            m_equ[1]->DoOdeProjection1(inarray,outarray,time);
        }

        void DriverMultiLevel::v_DoMultiOrderOdeRhs(const Array<OneD,const Array<OneD, NekDouble>> &inarray,
                                                           Array<OneD, Array<OneD,NekDouble>> &outarray,
                                                                        const NekDouble time)
        {
            ASSERTL1(m_equ[1],"Need to define EquationSystem[1]");
            m_equ[1]->DoOdeRhs1(inarray,outarray,time);
        }
    
    
        /**
         *
         */
        void DriverMultiLevel::v_InitObject(ostream &out)
        {
            Driver::v_InitObject(out);
            int nCycles=m_nLevels-1;

            Array<OneD,Array<OneD,DNekMatSharedPtr>> m_RestrictionMatrix(nCycles);
            Array<OneD,Array<OneD,DNekMatSharedPtr>>m_ProlongationMatrix(nCycles);

            for(int k=0;k<nCycles;k++)
            {

                //To Do: operate in StdElement, to save memory, need to only store nType Element
                int nElmts=m_equ[k]->GetExpSize();
                m_RestrictionMatrix[k]=Array<OneD,DNekMatSharedPtr>(nElmts);
                m_ProlongationMatrix[k]=Array<OneD,DNekMatSharedPtr>(nElmts);
                for (int i=0;i<nElmts;i++)
                {
                    LocalRegions::ExpansionSharedPtr LowerOrderExpansion=m_equ[k+1]->GetExp(i);
                    LocalRegions::ExpansionSharedPtr HigherOrderExpansion=m_equ[k]->GetExp(i);
                    int nHighOrderCoeffs=m_equ[k]->GetNcoeffs(i);
                    int nLowOrderCoeffs=m_equ[k+1]->GetNcoeffs(i);
                    m_RestrictionMatrix[k][i]=MemoryManager<DNekMat>::AllocateSharedPtr(nLowOrderCoeffs,nHighOrderCoeffs, eFULL, 0.0);
                    m_ProlongationMatrix[k][i]=MemoryManager<DNekMat>::AllocateSharedPtr(nHighOrderCoeffs,nLowOrderCoeffs, eFULL, 0.0);
                    LibUtilities::PointsKeyVector HigherOrderExpansionKeys,LowerOrderExpansionKeys;
                    LowerOrderExpansionKeys=LowerOrderExpansion->GetPointsKeys();
                    HigherOrderExpansionKeys=HigherOrderExpansion->GetPointsKeys();
                    StdRegions::StdMatrixKey LowerOrderMatKey(StdRegions::eBwdTrans,
                                            LowerOrderExpansion->DetShapeType(),
                                            *(LowerOrderExpansion));
                    DNekMatSharedPtr LowerOrderBwdMat = LowerOrderExpansion->GetStdMatrix(LowerOrderMatKey);
                    StdRegions::StdMatrixKey HigherOrderMatKey(StdRegions::eBwdTrans,
                                            HigherOrderExpansion->DetShapeType(),
                                            *(HigherOrderExpansion));
                    DNekMatSharedPtr HigherOrderBwdMat = HigherOrderExpansion->GetStdMatrix(HigherOrderMatKey);
                    LowerOrderExpansion->CreateInterpolationMatrix(HigherOrderExpansionKeys,HigherOrderBwdMat,m_RestrictionMatrix[k][i]);
                    HigherOrderExpansion->CreateInterpolationMatrix(LowerOrderExpansionKeys,LowerOrderBwdMat,m_ProlongationMatrix[k][i]);
                    // cout<<"RestrictionMatrix["<<k<<"]["<<i<<"]"<<endl;
                    // OutputMatrix(m_RestrictionMatrix[k][i]);
                    // cout<<"ProlongationMatrix["<<k<<"]["<<i<<"]"<<endl;
                    // OutputMatrix(m_ProlongationMatrix[k][i]);
                
                }
            }


            m_driverOperator.DefineMultiOrderProjection(&Driver::DoMultiOrderProjection, this);
            m_driverOperator.DefineMultiOrderOdeRhs(&Driver::DoMultiOrderOdeRhs, this);
            m_equ[0]->SetdriverOperator(m_driverOperator);
        }
    
    
        void DriverMultiLevel::v_Execute(ostream &out)
        
        {
            time_t starttime, endtime;
            NekDouble CPUtime;

            m_equ[0]->PrintSummary(out);

            time(&starttime);
            m_equ[0]->DoInitialise();   
            // m_driver.DefineMultiOrderOdeRhs(&m_equ[1]->GetDoOdeObj().DoOdeRhs,this);
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
        
        void DriverMultiLevel::PrintMatrix(DNekMatSharedPtr &Matrix)
        {
            int nrows                = Matrix->GetRows();
            int ncols                = Matrix->GetColumns();
            MatrixStorage matStorage = Matrix->GetStorageType();
            for (int i = 0; i < nrows; i++)
            {

                for (int j = 0; j < ncols; j++)
                {
                    cout << setprecision(1)  << i  << "    " << j 
                            << "     " << setprecision(16) << (*Matrix)(i, j) << endl;
                }
            
            }
        }

        void DriverMultiLevel::OutputMatrix(DNekMatSharedPtr &Matrix)
        {
            int rows = Matrix->GetRows();
            int cols = Matrix->GetColumns();
            ofstream outfile1;
            outfile1.open("./Matrix.txt");
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    outfile1 << i+1<< "     " << j+1  << "    "
                            << std::setprecision(16) << (*Matrix)(i, j) << endl;
                }
            }
            outfile1.close();
        }
    }
}

