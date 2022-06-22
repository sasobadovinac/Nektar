///////////////////////////////////////////////////////////////////////////////
//
// File: Expansion3D.h
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
// Description: Header file for Expansion3D routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPANSION3D_H
#define EXPANSION3D_H

#include <LocalRegions/Expansion.h>
#include <StdRegions/StdExpansion3D.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
#include <SpatialDomains/Geometry3D.h>

namespace Nektar {
    namespace LocalRegions {

        class Expansion2D;

        typedef std::shared_ptr<Expansion2D> Expansion2DSharedPtr;
        typedef std::weak_ptr<Expansion2D> Expansion2DWeakPtr;

        class Expansion3D;

        typedef std::shared_ptr<Expansion3D> Expansion3DSharedPtr;
        typedef std::weak_ptr<Expansion3D> Expansion3DWeakPtr;
        typedef std::vector<Expansion3DSharedPtr> Expansion3DVector;

        class Expansion3D : virtual public Expansion,
                            virtual public StdRegions::StdExpansion3D {
        public:
            LOCAL_REGIONS_EXPORT Expansion3D(
                    SpatialDomains::Geometry3DSharedPtr pGeom) : Expansion(
                    pGeom), StdExpansion3D(), m_requireNeg() {}

            LOCAL_REGIONS_EXPORT ~Expansion3D() override = default;

            LOCAL_REGIONS_EXPORT void SetTraceToGeomOrientation
                    (Array<OneD, NekDouble> &inout);

            LOCAL_REGIONS_EXPORT void SetFaceToGeomOrientation
                    (int face, Array<OneD, NekDouble> &inout);

            inline void AddHDGHelmholtzFaceTerms(
                    NekDouble tau,
                    int edge,
                    Array<OneD, NekDouble> &facePhys,
                    const StdRegions::VarCoeffMap &dirForcing,
                    Array<OneD, NekDouble> &outarray);

            inline void AddNormTraceInt(
                    int dir,
                    Array<OneD, ExpansionSharedPtr> &FaceExp,
                    Array<OneD, Array<OneD, NekDouble> > &faceCoeffs,
                    Array<OneD, NekDouble> &outarray);

            inline void AddNormTraceInt(
                    int dir,
                    Array<OneD, const NekDouble> &inarray,
                    Array<OneD, ExpansionSharedPtr> &FaceExp,
                    Array<OneD, NekDouble> &outarray,
                    const StdRegions::VarCoeffMap &varcoeffs);

            inline void AddFaceBoundaryInt(
                    int face,
                    ExpansionSharedPtr &FaceExp,
                    Array<OneD, NekDouble> &facePhys,
                    Array<OneD, NekDouble> &outarray,
                    const StdRegions::VarCoeffMap &varcoeffs
                    = StdRegions::NullVarCoeffMap);

            inline SpatialDomains::Geometry3DSharedPtr GetGeom3D() const;

            LOCAL_REGIONS_EXPORT void v_ReOrientTracePhysMap(
                    StdRegions::Orientation orient,
                    Array<OneD, int> &idmap,
                    int nq0,
                    int nq1) override;

            LOCAL_REGIONS_EXPORT void v_NormVectorIProductWRTBase(
                    const Array<OneD, const Array<OneD, NekDouble> > &Fvec,
                    Array<OneD, NekDouble> &outarray) override;

            LOCAL_REGIONS_EXPORT Array<OneD, unsigned int>
            GetEdgeInverseBoundaryMap(int eid);

            LOCAL_REGIONS_EXPORT Array<OneD, unsigned int>
            GetTraceInverseBoundaryMap(int fid,
                                       StdRegions::Orientation
                                       faceOrient = StdRegions::eNoOrientation,
                                       int P1 = -1, int P2 = -1);

            LOCAL_REGIONS_EXPORT void GetInverseBoundaryMaps(
                    Array<OneD, unsigned int> &vmap,
                    Array<OneD, Array<OneD, unsigned int> > &emap,
                    Array<OneD, Array<OneD, unsigned int> > &fmap);

            LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr CreateMatrix(
                    const MatrixKey &mkey);

        protected:
            std::map<int, NormalVector> m_faceNormals;

            void v_DGDeriv(
                    int dir,
                    const Array<OneD, const NekDouble> &incoeffs,
                    Array<OneD, ExpansionSharedPtr> &FaceExp,
                    Array<OneD, Array<OneD, NekDouble> > &faceCoeffs,
                    Array<OneD, NekDouble> &out_d) override;

            DNekMatSharedPtr v_GenMatrix(
                    const StdRegions::StdMatrixKey &mkey) override;

            void v_AddFaceNormBoundaryInt(
                    int face,
                    const ExpansionSharedPtr &FaceExp,
                    const Array<OneD, const NekDouble> &Fn,
                    Array<OneD, NekDouble> &outarray) override;

            void v_AddRobinMassMatrix(
                    int face,
                    const Array<OneD, const NekDouble> &primCoeffs,
                    DNekMatSharedPtr &inoutmat) override;

            StdRegions::Orientation v_GetTraceOrient(int face) override;

            void v_GetTracePhysVals(
                    int face,
                    const StdRegions::StdExpansionSharedPtr &FaceExp,
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray,
                    StdRegions::Orientation orient) override;


            void v_GenTraceExp(int traceid,
                               ExpansionSharedPtr &exp) override;

            void GetPhysFaceVarCoeffsFromElement(
                    int face,
                    ExpansionSharedPtr &FaceExp,
                    const Array<OneD, const NekDouble> &varcoeff,
                    Array<OneD, NekDouble> &outarray);

            virtual Array<OneD, NekDouble> v_GetnFacecdotMF(
                    int dir,
                    int face,
                    ExpansionSharedPtr &FaceExp_f,
                    const Array<OneD, const Array<OneD, NekDouble> > &normals,
                    const StdRegions::VarCoeffMap &varcoeffs);

            //-----------------------------
            // Low Energy Basis functions
            //-----------------------------

            LOCAL_REGIONS_EXPORT DNekMatSharedPtr
            v_BuildTransformationMatrix(
                    const DNekScalMatSharedPtr &r_bnd,
                    StdRegions::MatrixType matrixType) override;

            LOCAL_REGIONS_EXPORT DNekMatSharedPtr
            v_BuildInverseTransformationMatrix(
                    const DNekScalMatSharedPtr &transformationmatrix) override;

            LOCAL_REGIONS_EXPORT DNekMatSharedPtr v_BuildVertexMatrix(
                    const DNekScalMatSharedPtr &r_bnd) override;

            LOCAL_REGIONS_EXPORT void v_TraceNormLen(int traceid,
                                                     NekDouble &h,
                                                     NekDouble &p) override;

        private:
            // Do not add members here since it may lead to conflicts.
            // Only use this class for member functions

            std::vector<bool> m_requireNeg;
        };

        inline SpatialDomains::Geometry3DSharedPtr
        Expansion3D::GetGeom3D() const
        {
            return std::dynamic_pointer_cast<SpatialDomains::Geometry3D>(
                    m_geom);
        }
    } //end of namespace
} //end of namespace

#endif //EXPANSION3D_H
