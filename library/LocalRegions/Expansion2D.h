///////////////////////////////////////////////////////////////////////////////
//
// File Expansion2D.h
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
// Description: Header file for Expansion2D routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPANSION2D_H
#define EXPANSION2D_H

#include <boost/core/ignore_unused.hpp>

#include <LocalRegions/Expansion1D.h>
#include <StdRegions/StdExpansion2D.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
#include <SpatialDomains/Geometry2D.h>
#include <LocalRegions/MatrixKey.h>

namespace Nektar {
    namespace LocalRegions {
        class Expansion3D;

        typedef std::shared_ptr<Expansion3D> Expansion3DSharedPtr;
        typedef std::weak_ptr<Expansion3D> Expansion3DWeakPtr;

        class Expansion2D;

        typedef std::shared_ptr<Expansion2D> Expansion2DSharedPtr;
        typedef std::weak_ptr<Expansion2D> Expansion2DWeakPtr;
        typedef std::vector<Expansion2DSharedPtr> Expansion2DVector;

        class Expansion2D : virtual public Expansion,
                            virtual public StdRegions::StdExpansion2D {
        public:
            LOCAL_REGIONS_EXPORT Expansion2D(SpatialDomains::
                                             Geometry2DSharedPtr pGeom);

            LOCAL_REGIONS_EXPORT ~Expansion2D() override = default;

            LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr
            CreateMatrix(const MatrixKey &mkey);

            LOCAL_REGIONS_EXPORT void SetTraceToGeomOrientation(
                    Array<OneD, ExpansionSharedPtr> &EdgeExp,
                    Array<OneD, NekDouble> &inout);

            LOCAL_REGIONS_EXPORT Array<OneD, unsigned int>

            GetTraceInverseBoundaryMap(int eid);

            inline void AddNormTraceInt(
                    int dir,
                    Array<OneD, ExpansionSharedPtr> &EdgeExp,
                    Array<OneD, Array<OneD, NekDouble> > &edgeCoeffs,
                    Array<OneD, NekDouble> &outarray);

            inline void AddNormTraceInt(
                    int dir,
                    Array<OneD, const NekDouble> &inarray,
                    Array<OneD, ExpansionSharedPtr> &EdgeExp,
                    Array<OneD, NekDouble> &outarray,
                    const StdRegions::VarCoeffMap &varcoeffs);

            inline void AddEdgeBoundaryInt(
                    int edge,
                    ExpansionSharedPtr &EdgeExp,
                    Array<OneD, NekDouble> &edgePhys,
                    Array<OneD, NekDouble> &outarray,
                    const StdRegions::VarCoeffMap &varcoeffs = StdRegions::
                    NullVarCoeffMap);

            inline void AddHDGHelmholtzEdgeTerms(
                    NekDouble tau,
                    int edge,
                    Array<OneD, ExpansionSharedPtr> &EdgeExp,
                    Array<OneD, NekDouble> &edgePhys,
                    const StdRegions::VarCoeffMap &dirForcing,
                    Array<OneD, NekDouble> &outarray);

            inline void AddHDGHelmholtzTraceTerms(
                    NekDouble tau,
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, ExpansionSharedPtr> &EdgeExp,
                    const StdRegions::VarCoeffMap &dirForcing,
                    Array<OneD, NekDouble> &outarray);

            inline SpatialDomains::Geometry2DSharedPtr GetGeom2D() const;

            DNekMatSharedPtr
            v_GenMatrix(const StdRegions::StdMatrixKey &mkey) override;

            void v_GenTraceExp(int traceid,
                               ExpansionSharedPtr &exp) override;

        protected:
            std::vector<bool> m_requireNeg;

            LOCAL_REGIONS_EXPORT Array<OneD, NekDouble> v_GetMF(
                    int dir,
                    int shapedim,
                    const StdRegions::VarCoeffMap &varcoeffs) override;

            LOCAL_REGIONS_EXPORT Array<OneD, NekDouble> v_GetMFDiv(
                    int dir,
                    const StdRegions::VarCoeffMap &varcoeffs) override;

            LOCAL_REGIONS_EXPORT Array<OneD, NekDouble>
            v_GetMFMag(int dir,
                       const StdRegions::VarCoeffMap &varcoeffs) override;


            // Hybridized DG routines
            void v_DGDeriv(
                    int dir,
                    const Array<OneD, const NekDouble> &incoeffs,
                    Array<OneD, ExpansionSharedPtr> &EdgeExp,
                    Array<OneD, Array<OneD, NekDouble> > &edgeCoeffs,
                    Array<OneD, NekDouble> &out_d) override;

            void v_AddEdgeNormBoundaryInt(
                    int edge,
                    const ExpansionSharedPtr &EdgeExp,
                    const Array<OneD, const NekDouble> &Fx,
                    const Array<OneD, const NekDouble> &Fy,
                    Array<OneD, NekDouble> &outarray) override;

            void v_AddEdgeNormBoundaryInt(
                    int edge,
                    const ExpansionSharedPtr &EdgeExp,
                    const Array<OneD, const NekDouble> &Fn,
                    Array<OneD, NekDouble> &outarray) override;

            void v_AddRobinMassMatrix(
                    int edgeid,
                    const Array<OneD, const NekDouble> &primCoeffs,
                    DNekMatSharedPtr &inoutmat) override;

            void v_AddRobinTraceContribution(
                    int traceid,
                    const Array<OneD, const NekDouble> &primCoeffs,
                    const Array<OneD, NekDouble> &incoeffs,
                    Array<OneD, NekDouble> &coeffs) override;

            DNekMatSharedPtr v_BuildVertexMatrix(
                    const DNekScalMatSharedPtr &r_bnd) override;

            void GetPhysEdgeVarCoeffsFromElement(
                    int edge,
                    ExpansionSharedPtr &EdgeExp,
                    const Array<OneD, const NekDouble> &varcoeff,
                    Array<OneD, NekDouble> &outarray);

            Array<OneD, NekDouble> v_GetnEdgecdotMF(
                    int dir,
                    int edge,
                    ExpansionSharedPtr &EdgeExp_e,
                    const Array<OneD, const Array<OneD, NekDouble> > &normals,
                    const StdRegions::VarCoeffMap &varcoeffs);

            void v_ReOrientTracePhysMap
                    (StdRegions::Orientation orient,
                     Array<OneD, int> &idmap,
                     int nq0, int nq1) override;

            void v_SetUpPhysNormals(int edge) override;

            NekDouble v_VectorFlux(
                    const Array<OneD, Array<OneD, NekDouble> > &vec) override;

            void
            v_TraceNormLen(int traceid, NekDouble &h, NekDouble &p) override;

        };


        inline SpatialDomains::Geometry2DSharedPtr
        Expansion2D::GetGeom2D() const
        {
            return std::dynamic_pointer_cast<SpatialDomains::
            Geometry2D>(m_geom);
        }
    } //end of namespace
} //end of namespace

#endif //EXPANSION2D_H
