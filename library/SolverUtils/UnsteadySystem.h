///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadySystem.h
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
// Description: Generic timestepping for Unsteady solvers
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_UNSTEADYSYSTEM_H
#define NEKTAR_SOLVERUTILS_UNSTEADYSYSTEM_H

#include <LibUtilities/TimeIntegration/TimeIntegrationWrapper.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Filters/Filter.h> 

namespace Nektar
{
    namespace SolverUtils
    {
        /// Base class for unsteady solvers.
        class UnsteadySystem : public EquationSystem
        {
        public:
            /// Destructor
            SOLVER_UTILS_EXPORT virtual ~UnsteadySystem();

            /// Calculate the larger time-step mantaining the problem stable.
            SOLVER_UTILS_EXPORT NekDouble GetTimeStep(
                const Array<OneD, const Array<OneD, NekDouble> > &inarray);

            SOLVER_UTILS_EXPORT void SteadyStateResidual(
                int                         step, 
                Array<OneD, NekDouble>      &L2)
            {
                v_SteadyStateResidual(step,L2);
            }

            /// CFL safety factor (comprise between 0 to 1)(may be larger than 1 for implicit solvers).
            NekDouble m_cflSafetyFactor;
            NekDouble m_cflNonAcoustic;
            /// CFL growth rate
            NekDouble m_CFLGrowth;
            /// maximun cfl in cfl growth
            NekDouble m_CFLEnd;


        protected:
            int  m_ExtractRhsCalculator=0;
            
            /// Number of time steps between outputting status information.
            int                                             m_infosteps;

            int                                             m_nanSteps;
            /// Wrapper to the time integration scheme
            LibUtilities::TimeIntegrationWrapperSharedPtr   m_intScheme;
            /// The time integration scheme operators to use.
            LibUtilities::TimeIntegrationSchemeOperators    m_ode;

            // SolverUtils::DriverOperators                    m_driver;
            ///
            LibUtilities::TimeIntegrationSolutionSharedPtr  m_intSoln;
            ///
            NekDouble                                       m_epsilon;
            /// Indicates if explicit or implicit treatment of diffusion is used.
            bool                                            m_explicitDiffusion;
            /// Indicates if explicit or implicit treatment of advection is used.
            bool                                            m_explicitAdvection;
            /// Indicates if explicit or implicit treatment of reaction is used.
            bool                                            m_explicitReaction;
            /// Flag to determine if simulation should start in homogeneous
            /// forward transformed state.
            bool                                            m_homoInitialFwd;

            /// Tolerance to which steady state should be evaluated at
            NekDouble                                       m_steadyStateTol;
            /// Check for steady state at step interval
            int                                             m_steadyStateSteps;

            NekDouble                                       m_steadyStateRes    =1.0;
            NekDouble                                       m_steadyStateRes0   =1.0;

            /// Storage for previous solution for steady-state check
            Array<OneD, Array<OneD, NekDouble> >            m_previousSolution;
            // Steady-state residual file
            std::ofstream                                   m_errFile;

            std::vector<int>                                m_intVariables;

            std::vector<FilterSharedPtr>                    m_filters;

            /// at which time to evaluate the boundary conditions(used in unsteady time integrations)
            NekDouble                                       m_BndEvaluateTime;

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
            /// coefff of spacial derivatives(rhs or m_F in GLM) in calculating the residual of the whole equation(used in unsteady time integrations)
            NekDouble                                       m_TimeIntegLambda=0.0;

            NekDouble                                       m_TimeIntegLambdaPrcMat=0.0;

            NekDouble                                       m_Res0PreviousStep=-1.0;

            ///Solution of The kth iteration in the Newton method(Nonlinear iteration)
            Array<OneD,       Array<OneD, NekDouble> >      m_TimeIntegtSol_k;

            /// Solution at time step n(input valure from timeintegration)
            Array<OneD,       Array<OneD, NekDouble> >      m_TimeIntegtSol_n;
            /// Residual of the nonlinear system at the kth iteration in the Newton method(Nonlinear iteration)
            /// also the b of linearsys(Ax=b) stored to compute Jacobian_
            Array<OneD,       Array<OneD, NekDouble> >      m_SysEquatResid_k;

            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >  m_PrecMatVars;

            Array<OneD, Array<OneD, NekDouble> >            m_PrecMatVarsOffDiag;

            DNekBlkMatSharedPtr                             m_PrecMat;
            Array<OneD, Array<OneD, SNekBlkMatSharedPtr> >  m_PrecMatVarsSingle;
            SNekBlkMatSharedPtr                             m_PrecMatSingle;

            bool                                            m_flagPrecMatVarsSingle;
            bool                                            m_flagPrecondCacheOptmis;
            bool                                            m_flagImplItsStatistcs;

            Array<OneD, DNekBlkMatSharedPtr >               m_TraceJac;
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,NekDouble >>>>  m_TraceJacArray;

            Array<OneD, DNekBlkMatSharedPtr >               m_TraceJacDeriv;
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,NekDouble >>>>  m_TraceJacDerivArray;

            Array<OneD,       Array<OneD, NekDouble> >      m_TraceJacDerivSign;

            Array<OneD, SNekBlkMatSharedPtr >               m_TraceJacSingle;
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,NekSingle >>>>  m_TraceJacArraySingle;

            Array<OneD, SNekBlkMatSharedPtr >               m_TraceJacDerivSingle;
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,NekSingle >>>>  m_TraceJacDerivArraySingle;

            Array<OneD,       Array<OneD, NekSingle> >      m_TraceJacDerivSignSingle;
            
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,Array<OneD,NekSingle >>>>>  m_TraceIPSymJacArraySingle;
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,Array<OneD,NekDouble >>>>>  m_TraceIPSymJacArray;

            /// estimate the magnitude of each conserved varibles
            Array<OneD, NekDouble>                          m_magnitdEstimat;

            /// local time step(notice only for jfnk other see m_cflSafetyFactor)
            Array<OneD, NekDouble>                          m_locTimeStep;

            NekDouble   m_inArrayNorm=-1.0;

            bool m_CalcuPrecMatFlag     = true;
            ////////////////////////////////////////////////////////////////////////////
            //Paremeters to control the Direct Error
            int m_CalculateTemporalErrorCounter=0;
            NekDouble m_TemporalErrorNorm;
            Array<OneD,NekDouble> m_TemporalErrorNormArray;
            Array<OneD,Array<OneD,NekDouble>> m_TemporalError;
            int m_CalculateSpatialErrorCounter=0;
            bool m_CalculateSpatialErrorFlag=false;
            NekDouble m_SpatialErrorNorm;
            Array<OneD,NekDouble> m_SpatialErrorNormArray;
            Array<OneD,Array<OneD,NekDouble>> m_SpatialError;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            int m_CalcuPrecMatCounter  = std::numeric_limits<int>::max();

            int m_TotLinItePerStep=0;
            int m_StagesPerStep=1;

            int m_maxLinItePerNewton;

            int m_TotNewtonIts  =0;
            int m_TotGMRESIts   =0;
            int m_TotOdeRHS     =0;
            int m_TotImpStages  =0;

            /// flag to update artificial viscosity
            bool m_calcuPhysicalAV = true;

#endif
            /// Initialises UnsteadySystem class members.
            SOLVER_UTILS_EXPORT UnsteadySystem(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const SpatialDomains::MeshGraphSharedPtr& pGraph);

            /// Init object for UnsteadySystem class.
            SOLVER_UTILS_EXPORT virtual void v_InitObject();

            /// Get the maximum timestep estimator for cfl control.
            SOLVER_UTILS_EXPORT NekDouble MaxTimeStepEstimator();

            /// Solves an unsteady problem.
            SOLVER_UTILS_EXPORT virtual void v_DoSolve();

            /// Sets up initial conditions.
            SOLVER_UTILS_EXPORT virtual void v_DoInitialise();

            /// Print a summary of time stepping parameters.
            SOLVER_UTILS_EXPORT virtual void v_GenerateSummary(SummaryList& s);

            /// Print the solution at each solution point in a txt file
            SOLVER_UTILS_EXPORT virtual void v_AppendOutput1D(
                Array<OneD, Array<OneD, NekDouble> > &solution1D);

            SOLVER_UTILS_EXPORT virtual NekDouble v_GetTimeStep(
                const Array<OneD, const Array<OneD, NekDouble> > &inarray);

            SOLVER_UTILS_EXPORT virtual bool v_PreIntegrate(int step);
            SOLVER_UTILS_EXPORT virtual bool v_PostIntegrate(int step);

            SOLVER_UTILS_EXPORT virtual bool v_RequireFwdTrans()
            {
                return true;
            }

            SOLVER_UTILS_EXPORT virtual void v_SteadyStateResidual(
                int                         step, 
                Array<OneD, NekDouble>      &L2);

            SOLVER_UTILS_EXPORT void CheckForRestartTime(NekDouble &time, int &nchk);

            /// \brief Evaluate the SVV diffusion coefficient
            /// according to Moura's paper where it should
            /// proportional to h time velocity
            SOLVER_UTILS_EXPORT void SVVVarDiffCoeff(const Array<OneD, Array<OneD, NekDouble> >
                                                     vel,
                                                     StdRegions::VarCoeffMap &varCoeffMap);


        private:

            void InitializeSteadyState();

            bool CheckSteadyState(int step);
            bool CheckSteadyState(int step, NekDouble totCPUTime);
            NekDouble CalculateToleranceSafetyFactor();

            NekDouble CalculateToleranceErrorNorm();
        };

    }
}

#endif
