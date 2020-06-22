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

        void DriverMultiLevel::v_MultiLevel(const Array<OneD, NekDouble> &inarray,
                                                Array<OneD,NekDouble> &outarray,
                                                const bool UpDateOperatorflag,
                                                const int Level)
        {
            ASSERTL1(m_equ[Level],"Need to define EquationSystem[level]");
            
            m_equ[Level]->MultiLevel(inarray, outarray, Level, m_MultiLevelCoeffs[Level],
                       m_MultiLevelCoeffs[Level+1], UpDateOperatorflag);
        } 

        void DriverMultiLevel::v_MultiLvlJacMultiplyMatFree(
            const int                         Level,
            const TensorOfArray1D<NekDouble>  &inarray, 
            TensorOfArray1D<NekDouble>        &out, 
            const NekDouble                   time, 
            const NekDouble                   dtlamda, 
            const TensorOfArray2D<NekDouble>  &refFields, 
            const bool                        flagUpdateJac)
        {
            ASSERTL1(m_equ[Level],"Need to define EquationSystem[level]");
            m_equ[Level]->MatrixMultiply_MatrixFree_coeff(
                inarray,
                out,
                time,
                dtlamda,
                refFields,
                flagUpdateJac);
        }    

        void DriverMultiLevel::v_CalculateNextLevelPreconditioner(
                    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                    const NekDouble                                 time,
                    const NekDouble                                 lambda,
                    const int                                       Level)
        {
            ASSERTL1(m_equ[Level],"Need to define EquationSystem[level]");
            
            m_equ[Level]->CalculateNextLevelPreconditioner(inarray, time, lambda);
        }    
    
        /**
         *
         */
        void DriverMultiLevel::v_InitObject(ostream &out)
        {
            Driver::v_InitObject(out);

            ASSERTL0(m_nLevels > 1,"No need use MultiLevel");

            CalcMultiLevelMatrix();

            std::string vEquation = m_session->GetSolverInfo("EqType");
            if (m_session->DefinesSolverInfo("SolverType"))
            {
                vEquation = m_session->GetSolverInfo("SolverType");
            }

            int nCycles = m_nLevels - 1;
            
            for (int k = 0; k < nCycles; ++k)
            {
                vector<string>  MultiLevelFilename;
                for(int i=0;i<m_session->GetFilenames().size(); ++i)
                {
                    std::string TmpInputFile=m_session->GetFilenames()[i];
                    MultiLevelFilename.push_back(TmpInputFile);
                }

                LibUtilities::SessionReaderSharedPtr MultiLevelSession = 
                    LibUtilities::SessionReader::CreateInstance(
                                    0, NULL, MultiLevelFilename, 
                                    m_session->GetComm());
                SpatialDomains::MeshGraphSharedPtr MultiLevelGraph =
                    SpatialDomains::MeshGraph::Read(MultiLevelSession,
                            SpatialDomains::NullDomainRangeShPtr,
                            true,
                            m_graph->GetCompositeOrdering(),
                            m_graph->GetBndRegionOrdering(), 
                            m_multilvlncoeffOffset[k], 
                            m_multilvlnphyscOffset[k]);

                MultiLevelSession->SetTag("AdvectiveType","Convective");
                m_equ[k+1] = GetEquationSystemFactory().CreateInstance(
                    vEquation, MultiLevelSession, MultiLevelGraph);
            }  

            m_MultiLevelCoeffs = Array<OneD, int>(m_nLevels+1, 0);
            for(int k = 0; k < m_nLevels; ++k)
            {
                m_MultiLevelCoeffs[k] = m_equ[k]->GetNcoeffs();
            }
            //m_MultiLevelCoeffs[m_nLevels]=-1 used in MultiLevel's LowLevelCoeffs 
            //so that avoid repeated setting flag of lowest level
            m_MultiLevelCoeffs[m_nLevels] = -1;
            
            m_driverOperator.DefineCalculateNextLevelPreconditioner
                (&Driver::CalculateNextLevelPreconditioner, this);
            m_driverOperator.DefineMultiLevel(&Driver::MultiLevel, this);
            for (int k = 0;k < nCycles; ++k)
            {
                m_equ[k]->SetdriverOperator(m_driverOperator);
            }
        }  
        void DriverMultiLevel::CalcMultiLevelMatrix()
        {
            std::string vEquation = m_session->GetSolverInfo("EqType");
            if (m_session->DefinesSolverInfo("SolverType"))
            {
                vEquation = m_session->GetSolverInfo("SolverType");
            }

            int nCycles = m_nLevels - 1;
            m_RestrictionResidualMatrix  = 
                Array<OneD,Array<OneD,DNekMatSharedPtr>> (nCycles);
            m_RestrictionMatrix = 
                Array<OneD,Array<OneD,DNekMatSharedPtr>> (nCycles);
            m_ProlongationMatrix = 
                Array<OneD,Array<OneD,DNekMatSharedPtr>>(nCycles);
            for (int k = 0; k < nCycles; ++k)
            {
                vector<string>  MultiLevelFilename;
                //Notice, after partition, FileName_xml/P0000000.xml, 
                for(int i = 0; i < m_session->GetFilenames().size(); ++i)
                {
                    string TmpInputFile=m_session->GetFilenames()[i];
                    MultiLevelFilename.push_back(TmpInputFile);
                }

                LibUtilities::SessionReaderSharedPtr MultiLevelSession= 
                LibUtilities::SessionReader::CreateInstance(
                                0, NULL, MultiLevelFilename, 
                                m_session->GetComm());
                
                //Temporary Graph uses highest Quad number to Calculate Restriction and prolongation
                //This way can help jump levels
                SpatialDomains::MeshGraphSharedPtr MultiLevelGraph=
                SpatialDomains::MeshGraph::Read(MultiLevelSession,
                        SpatialDomains::NullDomainRangeShPtr,
                        true,
                        m_graph->GetCompositeOrdering(),
                        m_graph->GetBndRegionOrdering(), 
                        m_multilvlncoeffOffset[k], 0);

                MultiLevelSession->SetTag("AdvectiveType","Convective");
                m_equ[k+1] = GetEquationSystemFactory().CreateInstance(
                    vEquation, MultiLevelSession, MultiLevelGraph);
                
                int nElmts = m_equ[k]->GetExpSize();
                m_RestrictionResidualMatrix[k] = 
                    Array<OneD,DNekMatSharedPtr>(nElmts);
                m_RestrictionMatrix[k] = 
                    Array<OneD,DNekMatSharedPtr>(nElmts);
                m_ProlongationMatrix[k] = 
                    Array<OneD,DNekMatSharedPtr>(nElmts);
                for (int i = 0; i < nElmts; ++i)
                {
                    LocalRegions::ExpansionSharedPtr LowerLevelExpansion = 
                        m_equ[k+1]->GetExp(i);
                    LocalRegions::ExpansionSharedPtr HigherLevelExpansion = 
                        m_equ[k]->GetExp(i);
                    int nHighLevelCoeffs = m_equ[k]->GetNcoeffs(i);
                    int nLowLevelCoeffs = m_equ[k+1]->GetNcoeffs(i);
                    m_RestrictionResidualMatrix[k][i] = 
                        MemoryManager<DNekMat>::AllocateSharedPtr(
                            nLowLevelCoeffs,nHighLevelCoeffs, eFULL, 0.0);
                    m_RestrictionMatrix[k][i] = 
                        MemoryManager<DNekMat>::AllocateSharedPtr(
                            nLowLevelCoeffs,nHighLevelCoeffs, eFULL, 0.0);
                    m_ProlongationMatrix[k][i] = 
                        MemoryManager<DNekMat>::AllocateSharedPtr(
                            nHighLevelCoeffs,nLowLevelCoeffs, eFULL, 0.0);
                    LibUtilities::PointsKeyVector LowerLevelExpansionKeys = 
                        LowerLevelExpansion->GetPointsKeys();
                    LibUtilities::PointsKeyVector HigherLevelExpansionKeys =
                        HigherLevelExpansion->GetPointsKeys();
                    StdRegions::StdMatrixKey LowerLevelMatKey(
                        StdRegions::eBwdTrans,
                        LowerLevelExpansion->DetShapeType(),
                        *(LowerLevelExpansion));
                    DNekMatSharedPtr LowerLevelBwdMat = 
                        LowerLevelExpansion->GetStdMatrix(LowerLevelMatKey);
                    StdRegions::StdMatrixKey HigherLevelMatKey(
                        StdRegions::eBwdTrans,
                        HigherLevelExpansion->DetShapeType(),
                        *(HigherLevelExpansion));
                    DNekMatSharedPtr HigherLevelBwdMat = 
                        HigherLevelExpansion->GetStdMatrix(HigherLevelMatKey);
                    LowerLevelExpansion->CreateInterpolationMatrix(
                        HigherLevelExpansionKeys,HigherLevelBwdMat,
                        m_RestrictionMatrix[k][i]);
                    HigherLevelExpansion->CreateInterpolationMatrix(
                        LowerLevelExpansionKeys,LowerLevelBwdMat,
                        m_ProlongationMatrix[k][i]);
                    //Since have not found related DNekMatSharedPtr Transpose
                    for(int row=0;row<nLowLevelCoeffs;row++)
                    {
                        for(int col=0;col<nHighLevelCoeffs;col++)
                        {
                            (*m_RestrictionResidualMatrix[k][i])(row,col) = 
                                (*m_ProlongationMatrix[k][i])(col,row);
                        }
                    }
                }
                m_equ[k]->SetRestrictionResidualMatrix(
                    m_RestrictionResidualMatrix[k]);
                m_equ[k]->SetRestrictionMatrix(m_RestrictionMatrix[k]);
                m_equ[k]->SetProlongationMatrix(m_ProlongationMatrix[k]);
            }
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

