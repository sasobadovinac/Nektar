///////////////////////////////////////////////////////////////////////////////
//
// File: RoeSolver.h
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
// Description: Roe Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_RIEMANNSOLVER_ROESOLVER
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_RIEMANNSOLVER_ROESOLVER

#include <CompressibleFlowSolver/RiemannSolvers/CompressibleSolver.h>

namespace Nektar
{
    class RoeSolver : public CompressibleSolver
    {
    public:
        static RiemannSolverSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        {
            return RiemannSolverSharedPtr(
                new RoeSolver(pSession));
        }
        
        static std::string solverName;
        
    protected:
        RoeSolver(const LibUtilities::SessionReaderSharedPtr& pSession);
        
        virtual void v_PointSolve(
            double  rhoL, double  rhouL, double  rhovL, double  rhowL, double  EL,
            double  rhoR, double  rhouR, double  rhovR, double  rhowR, double  ER,
            double &rhof, double &rhouf, double &rhovf, double &rhowf, double &Ef);
        virtual void v_PointSolve(
            const Array<OneD, NekDouble> &Fwd, 
            const Array<OneD, NekDouble> &Bwd, 
            Array<OneD, NekDouble>       &flux, 
            const Array<OneD, NekDouble> &FwdJ, 
            const Array<OneD, NekDouble> &BwdJ);
        void RoeAverage(
            double rhoL, double rhouL, double rhovL, double rhowL, double EL,
            double rhoR, double rhouR, double rhovR, double rhowR, double ER,
            double &uRoe, double &vRoe, double &wRoe, double &hRoe, double &URoe,
            double &cRoe, double &ocRoe);
    };
}

#endif
