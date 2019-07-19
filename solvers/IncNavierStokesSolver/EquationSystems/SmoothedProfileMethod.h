///////////////////////////////////////////////////////////////////////////////
//
// File SmoothedProfileMethod.h
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
// Description: Smoothed Profile Method header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SMOOTHEDPROFILEMETHOD_H
#define NEKTAR_SOLVERS_SMOOTHEDPROFILEMETHOD_H

#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.h>

namespace Nektar
{
    class SmoothedProfileMethod: public VelocityCorrectionScheme
    {
    public:
        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph)
        {
            SolverUtils::EquationSystemSharedPtr p =
                MemoryManager<SmoothedProfileMethod>::AllocateSharedPtr(
                    pSession, pGraph);
            p->InitObject();
            return p;
        }

        /// Name of class
        static std::string className;

        //Constructor
        SmoothedProfileMethod(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph);

        // Destructor
        virtual ~SmoothedProfileMethod();

        virtual void v_InitObject();

        // Solves the linear part of the velocity correction scheme incluiding
        // the SPM method calculation for 'fs'
        void SolveUnsteadyStokesSystem(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble time,
                    const NekDouble a_iixDt)
        {
            v_SolveUnsteadyStokesSystem(inarray, outarray, time, a_iixDt);
        }

    protected:
        /// Correction pressure field for SPM
        MultiRegions::ExpListSharedPtr m_pressureP;
        /// DEBUG: u_p constant and equal to 0
        Array<OneD, NekDouble> m_up;
        /// Stiffly-stable scheme \gamma_0 coefficient
        NekDouble m_gamma0;

        // Calculates the shape function values
        // (only for non-moving boundaries)
        void CalcPhi(const MultiRegions::ExpListSharedPtr &expansion,
                    Array<OneD, NekDouble> &phi);
        // Calculates the virtual force 'fs'
        void IBForcing(const Array<OneD, const Array<OneD, NekDouble> > &fields,
                    NekDouble dt,
                    Array<OneD, Array<OneD, NekDouble> > &f_s);
        // Calculates the virtual force 'fs' in the boundary 'BndExp'
        void IBForcingBC(int bndInd,
                         const MultiRegions::ExpListSharedPtr &BndExp,
                         NekDouble dt,
                         Array<OneD, Array<OneD, NekDouble> > &f_s);
        // Interface for 'v_SolveUnsteadyStokesSystem'
        virtual void v_SolveUnsteadyStokesSystem(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble time,
                    const NekDouble a_iixDt);
        // Sets the parameters and BCs for the Poisson equation
        virtual void v_SetUpCorrectionPressure(
                    const Array<OneD, const Array<OneD, NekDouble> > &fields,
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    const NekDouble aii_Dt);
        // Solves the Poisson equation for the correction pressure
        virtual void v_SolveCorrectionPressure(
                    const Array<OneD, NekDouble> &Forcing);
        // Explicitly corrects the velocity by using the force 'fs'
        virtual void v_SolveCorrectedVelocity(
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    Array<OneD, Array<OneD, NekDouble> > &fields,
                    const NekDouble dt);
        // Set proper BCs for the corrected pressure 'p_p'
        virtual void v_SetCorrectionPressureBCs(NekDouble dt);

    private:
    };

    typedef std::shared_ptr<SmoothedProfileMethod>
            SmoothedProfileMethodSharedPtr;

} // end of namespace

#endif // SMOOTHEDPROFILEMETHOD_H
