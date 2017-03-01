///////////////////////////////////////////////////////////////////////////////
//
// File BS
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
// Description: B-S solve routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_OPTIONPRICINGSOLVER_EQUATIONSYSTEMS_BS_H
#define NEKTAR_SOLVERS_OPTIONPRICINGSOLVER_EQUATIONSYSTEMS_BS_H

#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Diffusion/Diffusion.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>


namespace Nektar
{
    class BS : public SolverUtils::AdvectionSystem
    {
    public:
        friend class MemoryManager<BS>;

        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession){
            SolverUtils::EquationSystemSharedPtr p = MemoryManager<BS>::
                AllocateSharedPtr(pSession);
            
            p->InitObject();
            return p;
        }
        
        /// Name of class
        static std::string className;

        /// Destructor
        virtual ~BS();

    protected:
        int       m_intSteps;
        NekDouble m_cflSafetyFactor;
        int       m_infosteps;
        int       m_minsubsteps;
        bool      m_subSteppingScheme;
        bool      m_useSpecVanVisc;
        
        // Cut-off ratio from which to start decayhing modes
        NekDouble m_sVVCutoffRatio;
        
        // Diffusion coefficient of SVV modes
        NekDouble m_sVVDiffCoeff;
        
        SolverUtils::RiemannSolverSharedPtr  m_riemannSolver;
        SolverUtils::DiffusionSharedPtr      m_diffusion;
        Array<OneD, Array<OneD, NekDouble> > m_stockPrice;
        Array<OneD, NekDouble>               m_traceSn;

        // Plane (used only for Discontinous projection with 3DHomogenoeus1D)
        int m_planeNumber;

        /// Session reader
        BS(const LibUtilities::SessionReaderSharedPtr& pSession);

        /// Evaluate the flux at each solution point for the advection part
        void GetFluxVectorAdv(
            const Array<OneD, Array<OneD, NekDouble> >               &physfield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);

        /// Evaluate the flux at each solution point for the diffusion part
        void GetFluxVectorDiff(
            const int i, 
            const int j,
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
                  Array<OneD, Array<OneD, NekDouble> > &derivatives,
                  Array<OneD, Array<OneD, NekDouble> > &flux);

        /// Compute the RHS
        virtual void DoOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time);

        /// Perform the projection
        void DoOdeProjection(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const NekDouble time);

        /// Solve implicitly the diffusion term
        virtual void DoImplicitSolve(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                  Array<OneD,       Array<OneD, NekDouble> >&outarray,
            NekDouble time,
            NekDouble lambda);
        
        /// Get the normal stock price based on m_stockPrice
        Array<OneD, NekDouble> &GetNormalStockPrice();

        /// Get the normal stock based on input stock price
        Array<OneD, NekDouble> &GetNormalStock(
            const Array<OneD, const Array<OneD, NekDouble> > &stockField);
        
        /// Initialise the object
        virtual void v_InitObject();

        /// Print Summary
        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);

        /// PreIntegration step for substepping. 
        virtual bool v_PreIntegrate(int step);

        // SubsStepping methods -> Probably could be set up in separate class
        void SubStepAdvance(
            const LibUtilities::TimeIntegrationSolutionSharedPtr
                      &integrationSoln,
            int       nstep,
            NekDouble time);
        
        NekDouble GetSubstepTimeStep();
        
        void SetUpSubSteppingTimeIntegration(
            int intMethod,
            const LibUtilities::TimeIntegrationWrapperSharedPtr
                      &IntegrationScheme);

        void SubStepAdvection(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const NekDouble                                  time);

        void SubStepProjection(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const NekDouble                                  time);

        void AddAdvectionPenaltyFlux(
            const Array<OneD, const Array<OneD, NekDouble> > &velfield,
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD,       Array<OneD, NekDouble> > &Outarray);
        
        Array<OneD, NekDouble> GetMaxStdVelocity(
            const Array<OneD, Array<OneD,NekDouble> > inarray);

        LibUtilities::TimeIntegrationWrapperSharedPtr
            m_subStepIntegrationScheme;
        
        LibUtilities::TimeIntegrationSchemeOperators
            m_subStepIntegrationOps;


        
    private:
        
        // Risk-free interest rate and volatility
        NekDouble m_interest;
        NekDouble m_volatility;

        /// Variable diffusivity
        StdRegions::VarCoeffMap m_vardiff;

    };
}

#endif
