///////////////////////////////////////////////////////////////////////////////
//
// File: EnforceSU.h
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
// Description: Modified Riemann invariant boundary condition.
//              Enforcing the entropy and normal velocity at the inflow boundary;
//              Enforcing the Riemann invariant at the outflow boundary.
//              The input can be either VALUE or FILE. 
// 
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_BNDCOND_ENFORCESU
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_BNDCOND_ENFORCESU

#include "CFSBndCond.h"


namespace Nektar
{

/**
* @brief Outflow characteristic boundary conditions for compressible
* flow problems.
*/
class EnforceSU : public CFSBndCond
{
    public:

        friend class MemoryManager<EnforceSU>;

        /// Creates an instance of this class
        static CFSBndCondSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const Array<OneD, Array<OneD, NekDouble> >& pTraceNormals,
                const int pSpaceDim, const int bcRegion, const int cnt)
        {
            CFSBndCondSharedPtr p = MemoryManager<EnforceSU>::
                                    AllocateSharedPtr(pSession, pFields,
                                    pTraceNormals, pSpaceDim, bcRegion, cnt);
            return p;
        }

        ///Name of the class
        static std::string className;

    protected:
    int m_npts; 
    // Reference rho on boundary 
    Array<OneD, NekDouble>               m_rhoBC;
    // Reference Velocity on BC
    Array<OneD, Array<OneD, NekDouble> > m_velBC;
    // Reference pressure on boundary 
    Array<OneD, NekDouble>               m_pBC;
    /// Reference normal velocity
    Array<OneD, NekDouble>               m_VnInf;
    // Mapping from boundary to Trace values
    Array<OneD, int> m_bndToTraceMap; 
    
    // Arrays of arrays pointing to the boundary condition physical
    // space for the specified region.
    Array<OneD, Array<OneD, NekDouble> > m_bndPhys;

    virtual void v_Apply(
            Array<OneD, Array<OneD, NekDouble> >               &Fwd,
            Array<OneD, Array<OneD, NekDouble> >               &physarray,
            const NekDouble                                    &time);

    void GenerateRotationMatrices(
            const Array<OneD, const Array<OneD, NekDouble> > &normals);

    void FromToRotation(
                Array<OneD, const NekDouble> &from,
                Array<OneD, const NekDouble> &to,
                NekDouble                    *mat);

    SOLVER_UTILS_EXPORT void rotateToNormal  (
                const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                const Array<OneD, const Array<OneD, NekDouble> > &normals,
                const Array<OneD, const Array<OneD, NekDouble> > &vecLocs,
                      Array<OneD,       Array<OneD, NekDouble> > &outarray);

    SOLVER_UTILS_EXPORT void rotateFromNormal(
                const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                const Array<OneD, const Array<OneD, NekDouble> > &normals,
                const Array<OneD, const Array<OneD, NekDouble> > &vecLocs,
                      Array<OneD,       Array<OneD, NekDouble> > &outarray);


    private:
        EnforceSU(const LibUtilities::SessionReaderSharedPtr& pSession,
               const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
               const Array<OneD, Array<OneD, NekDouble> >& pTraceNormals,
               const int pSpaceDim,
               const int bcRegion,
               const int cnt);
        
        virtual ~EnforceSU(void){};
};

}

#endif
