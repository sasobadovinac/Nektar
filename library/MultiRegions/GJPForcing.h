///////////////////////////////////////////////////////////////////////////////
//
// File GJPForcing.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: GJP data 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_GJPFORCING_H
#define NEKTAR_LIBS_MULTIREGIONS_GJPFORCING_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/DisContField.h>
#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>

namespace Nektar
{
    namespace MultiRegions
    {

        class GJPForcing
        {
        public:
            MULTI_REGIONS_EXPORT GJPForcing(ExpListSharedPtr field);

            MULTI_REGIONS_EXPORT ~GJPForcing() {};

            void Apply(const Array<OneD, NekDouble>  &inarray,
                       Array<OneD, NekDouble> &outarray,
                       bool OutArrayInCoeffSpace,
                       const Array<OneD, NekDouble> &pUnorm = NullNekDouble1DArray,
                       const NekDouble scale = 1.0) const;
            

            Array<OneD, Array<OneD, NekDouble>> &GetTraceNormals(void)
            {
                return m_traceNormals;
            }

            int GetNumTracePts(void) const
            {
                return m_dgfield->GetTrace()->GetTotPoints();
            }

            bool IsSemiImplicit() const 
            {
                return m_useGJPSemiImplicit;
            }
        private:
            int m_coordDim; 
            int m_traceDim; 

            /// Number of planes in expansion to be stabilised for
            /// Homgoeneous expansion
            int m_nplanes;
            
            bool m_useGJPSemiImplicit; 

            // Trace normals 
            Array<OneD, Array<OneD, NekDouble> > m_traceNormals; 
            
            /// DG expansion for projection evalaution along trace
            MultiRegions::ExpListSharedPtr m_dgfield;    
            /// LocaTraceToTraceMap 
            MultiRegions::LocTraceToTraceMapSharedPtr  m_locTraceToTraceMap;
            /// Local Elemental trace expansions
            MultiRegions::ExpListSharedPtr m_locElmtTrace;

            /// Scale factor for phys values along trace involving the lcoal
            /// normals and tangent geometric factors n
            Array<OneD, Array<OneD, NekDouble>> m_scalTrace;
            
            std::vector<std::pair<int,Array<OneD, DNekMatSharedPtr>>>
                m_StdDBaseOnTraceMat; 

            void SetUpExpansionInfoMapForGJP
                             (SpatialDomains::MeshGraphSharedPtr graph,
                              std::string variable);

            void MultiplyByStdDerivBaseOnTraceMat(int i,
                                           Array<OneD, NekDouble> &in,
                                           Array<OneD, NekDouble> &out) const;
        };

        typedef std::shared_ptr<GJPForcing>  GJPForcingSharedPtr;
    }
}
#endif
