///////////////////////////////////////////////////////////////////////////////
//
// File: AdvectionWeakDG.h
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
// Description: Weak DG advection class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_ADVECTIONWEAKDG
#define NEKTAR_SOLVERUTILS_ADVECTIONWEAKDG

#include <SolverUtils/Advection/Advection.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <StdRegions/StdSegExp.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/DisContField1D.h>
#include <boost/math/special_functions/gamma.hpp>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>

#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        class AdvectionWeakDG : public Advection
        {
        public:
            static AdvectionSharedPtr create(std::string advType)
            {
                return AdvectionSharedPtr(new AdvectionWeakDG());
            }

            static std::string type;
            
            Array<OneD, NekDouble>  m_dx;
            Array<OneD, NekDouble>  m_dxFwd;
            Array<OneD, NekDouble>  m_dxBwd;

            Array<OneD, Array<OneD, NekDouble> > m_gmat;
            DNekMatSharedPtr                     m_Ixm;
            DNekMatSharedPtr                     m_Ixp;
            
        protected:
            AdvectionWeakDG();

            Array<OneD, Array<OneD, NekDouble> >               m_traceNormals;

            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);

            virtual void v_Advect(
                const int                                          nConvective,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &advVel,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const NekDouble                                   &time,
                const Array<OneD, Array<OneD, NekDouble> >
                                  &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> >
                                  &pBwd = NullNekDoubleArrayofArray);
            
            virtual void v_SetupFluxLength(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields);
        };
    }
}

#endif
