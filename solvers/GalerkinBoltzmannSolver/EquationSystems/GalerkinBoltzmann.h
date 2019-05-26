///////////////////////////////////////////////////////////////////////////////
//
// File GalkerinBoltzmann.h
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
// Description: Unsteady advection solve routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_GALERKINBOLTZMANN_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_GALERKINBOLTZMANN_H

#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/AdvectionSystem.h>

namespace Nektar
{
    class GalerkinBoltzmann : public SolverUtils::AdvectionSystem
    {
    public:
        friend class MemoryManager<GalerkinBoltzmann>;

        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr& pGraph)
        {
            SolverUtils::EquationSystemSharedPtr p = MemoryManager<
                GalerkinBoltzmann>::AllocateSharedPtr(pSession, pGraph);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;

        /// Destructor
        virtual ~GalerkinBoltzmann();

    protected:

        // Plane (used only for Discontinous projection
        //        with 3DHomogenoeus1D expansion)
        int                                     m_planeNumber;

        /// Session reader
        GalerkinBoltzmann(const LibUtilities::SessionReaderSharedPtr& pSession,
                          const SpatialDomains::MeshGraphSharedPtr& pGraph);

        /// Compute the RHS
        void DoOdeRhs(
            const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                  Array<OneD,         Array<OneD, NekDouble> > &outarray,
            const NekDouble time);

        /// Compute the projection
        void DoOdeProjection(
            const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                  Array<OneD,         Array<OneD, NekDouble> > &outarray,
            const NekDouble time);

        void EvaluateIProductWRTDerivBaseVolFlux
              ( const Array<OneD, const  Array<OneD, NekDouble> >&physfield,
                Array<OneD,        Array<OneD, NekDouble>        >&volflux);

        void EvaluateNormalNumFlux(const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
                                   const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                                   Array<OneD, Array<OneD, NekDouble> > &numflux);

        void AddSourceTerms( const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                             Array<OneD,        Array<OneD, NekDouble>        >&outarray);

        /// Initialise the object
        virtual void v_InitObject();

        /// Print Summary
        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);

    private:
        /// Square Rood of RT
        NekDouble m_sqrtRT;
        NekDouble m_tau;
    };
}

#endif
