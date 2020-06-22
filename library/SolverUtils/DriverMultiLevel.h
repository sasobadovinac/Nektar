///////////////////////////////////////////////////////////////////////////////
//
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
// Description: Driver class for the stability solver
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DRIVERMULTILEVEL_H
#define NEKTAR_SOLVERUTILS_DRIVERMULTILEVEL_H

#include <SolverUtils/Driver.h>

namespace Nektar
{
    namespace SolverUtils
    {
        /// Base class for the development of solvers.
        class DriverMultiLevel: public Driver
        {
        public:
            friend class MemoryManager<DriverMultiLevel>;

            /// Creates an instance of this class
            static DriverSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const SpatialDomains::MeshGraphSharedPtr& pGraph)
            {
                DriverSharedPtr p = MemoryManager<DriverMultiLevel>
                    ::AllocateSharedPtr(pSession, pGraph);
                p->InitObject();
                return p;
            }

            void CalcMultiLevelMatrix();

            void PrintMatrix(DNekMatSharedPtr &Matrix);
    
            void OutputMatrix(DNekMatSharedPtr &Matrix);
	
            ///Name of the class
            static std::string className;

        protected:
            /// Constructor
            SOLVER_UTILS_EXPORT DriverMultiLevel(
                const LibUtilities::SessionReaderSharedPtr pSession,
                const SpatialDomains::MeshGraphSharedPtr pGraph);

            /// Destructor
            SOLVER_UTILS_EXPORT virtual ~DriverMultiLevel();

            SOLVER_UTILS_EXPORT virtual void v_MultiLevel(
                    const Array<OneD, NekDouble> &inarray,
                    Array<OneD,NekDouble>  &outarray,
                    const bool UpDateOperatorflag,
                    const int Level);

            SOLVER_UTILS_EXPORT virtual void  v_MultiLvlJacMultiplyMatFree(
                const int                         Level,
                const TensorOfArray1D<NekDouble>  &inarray, 
                TensorOfArray1D<NekDouble>        &out, 
                const NekDouble                   time, 
                const NekDouble                   dtlamda, 
                const TensorOfArray2D<NekDouble>  &refFields, 
                const bool                        flagUpdateJac);

            SOLVER_UTILS_EXPORT virtual void  v_CalculateNextLevelPreconditioner(
                    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                    const NekDouble                                 time,
                    const NekDouble                                 lambda,
                    const int                                       Level);
        
            /// Second-stage initialisation
            SOLVER_UTILS_EXPORT virtual void v_InitObject(std::ostream &out = std::cout);

            /// Virtual function for solve implementation.
            SOLVER_UTILS_EXPORT virtual void v_Execute(std::ostream &out = std::cout);
		
            static std::string driverLookupId;
	};
    }	
} //end of namespace

#endif //NEKTAR_SOLVERUTILS_DriverMultiLevel_H

