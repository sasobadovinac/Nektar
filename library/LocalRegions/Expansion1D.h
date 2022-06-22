///////////////////////////////////////////////////////////////////////////////
//
// File Expansion1D.h
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
// Description: Header file for Expansion1D routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPANSION1D_H
#define EXPANSION1D_H

#include <LocalRegions/Expansion.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
#include <StdRegions/StdExpansion1D.h>
#include <SpatialDomains/Geometry1D.h>

namespace Nektar {
    namespace LocalRegions {
        class Expansion2D;

        typedef std::shared_ptr<Expansion2D> Expansion2DSharedPtr;
        typedef std::weak_ptr<Expansion2D> Expansion2DWeakPtr;

        class Expansion1D;

        typedef std::shared_ptr<Expansion1D> Expansion1DSharedPtr;
        typedef std::weak_ptr<Expansion1D> Expansion1DWeakPtr;
        typedef std::vector<Expansion1DSharedPtr> Expansion1DVector;

        class Expansion1D : virtual public Expansion,
                            virtual public StdRegions::StdExpansion1D {
        public:
            LOCAL_REGIONS_EXPORT Expansion1D(SpatialDomains::
                                             Geometry1DSharedPtr pGeom)
                    : Expansion(pGeom),
                      StdExpansion1D()
            {
            }

            LOCAL_REGIONS_EXPORT ~Expansion1D() override = default;

            LOCAL_REGIONS_EXPORT void AddNormTraceInt(
                    int dir,
                    Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray);


            void AddHDGHelmholtzTraceTerms(
                    NekDouble tau,
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray);

            inline SpatialDomains::Geometry1DSharedPtr GetGeom1D() const;

        protected:
            DNekMatSharedPtr v_GenMatrix(
                    const StdRegions::StdMatrixKey &mkey) override;

            void v_AddRobinMassMatrix(
                    int vert,
                    const Array<OneD, const NekDouble> &primCoeffs,
                    DNekMatSharedPtr &inoutmat) override;

            void v_AddRobinEdgeContribution(
                    int vert,
                    const Array<OneD, const NekDouble> &primCoeffs,
                    const Array<OneD, NekDouble> &incoeffs,
                    Array<OneD, NekDouble> &coeffs);

            NekDouble v_VectorFlux(
                    const Array<OneD, Array<OneD, NekDouble> > &vec) override;

            void v_NormalTraceDerivFactors(
                    Array<OneD, Array<OneD, NekDouble> > &factors,
                    Array<OneD, Array<OneD, NekDouble> > &d0factors,
                    Array<OneD, Array<OneD, NekDouble> > &d1factors) override;

            virtual const NormalVector &
            v_GetTraceNormal(int edge) const final;

            void v_ReOrientTracePhysMap(
                    StdRegions::Orientation orient,
                    Array<OneD, int> &idmap,
                    int nq0, int nq1) override;

            void
            v_TraceNormLen(int traceid, NekDouble &h, NekDouble &p) override;

        private:
        };

        inline SpatialDomains::Geometry1DSharedPtr Expansion1D
        ::GetGeom1D() const
        {
            return std::dynamic_pointer_cast<SpatialDomains
            ::Geometry1D>(m_geom);
        }
    } //end of namespace
} //end of namespace

#endif //EXPANSION1D_H

