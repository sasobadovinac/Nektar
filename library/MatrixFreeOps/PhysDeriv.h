///////////////////////////////////////////////////////////////////////////////
//
// File: PhysDeriv.h
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MF_PHYSDERIV_H
#define MF_PHYSDERIV_H

#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "Operator.hpp"
#include "PhysDerivKernels.hpp"

namespace Nektar
{
namespace MatrixFree
{

// As each opertor has seven shapes over three dimension to get to the
// "work" each operator uses a series of preprocessor directives based
// on the dimension and shape type so to limit the code
// inclusion. This keeps the library size as minimal as possible while
// removing runtime conditionals.

// The preprocessor directives, SHAPE_DIMENSION_?D and SHAPE_TYPE_*
// are constructed by CMake in the CMakeLists.txt file which uses an
// implementation file.  See the CMakeLists.txt files for more
// details.
template <LibUtilities::ShapeType SHAPE_TYPE, bool DEFORMED = false>
struct PhysDerivTemplate
    : public PhysDeriv,
      public Helper<LibUtilities::ShapeTypeDimMap[SHAPE_TYPE], DEFORMED>
{
    PhysDerivTemplate(std::vector<LibUtilities::BasisSharedPtr> basis,
                      int nElmt)
        : PhysDeriv(basis, nElmt),
          Helper<LibUtilities::ShapeTypeDimMap[SHAPE_TYPE], DEFORMED>(basis,
                                                                      nElmt)
    {
        constexpr auto DIM = LibUtilities::ShapeTypeDimMap[SHAPE_TYPE];

        if (DIM == 1)
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients(SHAPE_TYPE,
                                                            this->m_nm[0]);
        }
        else if (DIM == 2)
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients(
                SHAPE_TYPE, this->m_nm[0], this->m_nm[1]);
        }
        else if (DIM == 3)
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients(
                SHAPE_TYPE, this->m_nm[0], this->m_nm[1], this->m_nm[2]);
        }
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
    {
        return std::make_shared<PhysDerivTemplate<SHAPE_TYPE, DEFORMED>>(basis,
                                                                         nElmt);
    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, Array<OneD, NekDouble>> &output) final
    {
#include "SwitchPoints.h"
    }

    // There must be separate 1D, 2D, and 3D operators because the
    // helper base class is based on the shape type dim. Thus indices
    // for the basis must be respected.

    // Further there are duplicate 1D, 2D, and 3D operators so to
    // allow for compile time array sizes to be part of the template
    // and thus gain loop unrolling. The only difference is the
    // location of the size value, in the template or the function:
    // foo<int bar>() {} vs foo() { int bar = ... }

#if defined(SHAPE_DIMENSION_1D)

    // Non-size based operator.
    void operator1D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, Array<OneD, NekDouble>> &output)
    {
        ASSERTL1(output.size() > 0, "PhysDerivTemplate::Operator1D: Cannot "
                                    "call 1D routine with no output.");
        ASSERTL1(output.size() <= 3, "PhysDerivTemplate::Operator1D: Operator "
                                     "not set up for other dimensions.")

        const auto nq0 = m_basis[0]->GetNumPoints();

        const auto nqTot    = nq0;
        const auto nqBlocks = nqTot * vec_t::width;

        const auto ndf        = output.size();
        constexpr int max_ndf = 3;

        auto *inptr = &input[0];

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        // Workspace for kernels - also checks preconditions
        // PhysDeriv1DWorkspace<SHAPE_TYPE>(nq0);

        // Call 1D kernel
        const vec_t *df_ptr;
        vec_t df_tmp[max_ndf]; // max_ndf is a constexpr

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut[max_ndf];
        NekDouble *out_ptr[max_ndf]; // max_ndf is a constexpr

        for (int d = 0; d < ndf; ++d)
        {
            tmpOut[d].resize(nqTot);
            out_ptr[d] = &output[d][0];
        }

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            df_ptr = &((*this->m_df)[dfSize * e]);

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            PhysDeriv1DKernel(nq0, ndf, tmpIn, this->m_D[0], df_ptr, df_tmp,
                              tmpOut);

            // De-interleave and store data
            for (int d = 0; d < ndf; ++d)
            {
                deinterleave_store(tmpOut[d], nqTot, out_ptr[d]);
                out_ptr[d] += nqBlocks;
            }

            inptr += nqBlocks;
        }
    }

    // Size based template version.
    template <int nq0>
    void operator1D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, Array<OneD, NekDouble>> &output)
    {
        ASSERTL1(output.size() > 0, "PhysDerivTemplate::Operator1D: Cannot "
                                    "call 1D routine with no output.");
        ASSERTL1(output.size() <= 3, "PhysDerivTemplate::Operator1D: Operator "
                                     "not set up for other dimensions.")

        constexpr auto nqTot    = nq0;
        constexpr auto nqBlocks = nqTot * vec_t::width;

        const auto ndf        = output.size();
        constexpr int max_ndf = 3;

        auto *inptr = &input[0];

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        // Workspace for kernels - also checks preconditions
        PhysDeriv1DWorkspace<SHAPE_TYPE>(nq0);

        // Call 1D kernel
        const vec_t *df_ptr;
        vec_t df_tmp[max_ndf]; // max_ndf is a constexpr

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut[max_ndf];
        NekDouble *out_ptr[max_ndf]; // max_ndf is a constexpr

        for (int d = 0; d < ndf; ++d)
        {
            tmpOut[d].resize(nqTot);
            out_ptr[d] = &output[d][0];
        }

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            df_ptr = &((*this->m_df)[dfSize * e]);

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            PhysDeriv1DKernel(nq0, ndf, tmpIn, this->m_D[0], df_ptr, df_tmp,
                              tmpOut);

            // De-interleave and store data
            for (int d = 0; d < ndf; ++d)
            {
                deinterleave_store(tmpOut[d], nqTot, out_ptr[d]);
                out_ptr[d] += nqBlocks;
            }

            inptr += nqBlocks;
        }
    }

#elif defined(SHAPE_DIMENSION_2D)

    // Non-size based operator.
    void operator2D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, Array<OneD, NekDouble>> &output)
    {
        ASSERTL1(output.size() > 1, "PhysDerivTemplate::Operator2D: Cannot "
                                    "call 2D routine with one output.");
        ASSERTL1(output.size() <= 3, "PhysDerivTemplate::Operator2D: Operator "
                                     "not set up for 3D coordinates.");

        const auto nq0 = m_basis[0]->GetNumPoints();
        const auto nq1 = m_basis[1]->GetNumPoints();

        const auto nqTot    = nq0 * nq1;
        const auto nqBlocks = nqTot * vec_t::width;

        const auto outdim        = output.size();
        constexpr int max_outdim = 3;

        const auto ndf        = 2 * outdim;
        constexpr int max_ndf = 2 * max_outdim;

        auto *inptr = &input[0];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut[max_ndf];
        NekDouble *out_ptr[max_outdim]; // max_outdim is a constexpr

        for (int d = 0; d < outdim; ++d)
        {
            tmpOut[d].resize(nqTot);
            out_ptr[d] = &output[d][0];
        }

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        // Workspace for kernels - also checks preconditions
        // PhysDeriv2DWorkspace<SHAPE_TYPE>(nq0, nq1);

        const vec_t *df_ptr;
        vec_t df_tmp[max_ndf]; // max_ndf is a constexpr

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            df_ptr = &((*this->m_df)[dfSize * e]);

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            PhysDeriv2DKernel(nq0, nq1, outdim, tmpIn, this->m_Z[0],
                              this->m_Z[1], this->m_D[0], this->m_D[1], df_ptr,
                              df_tmp, tmpOut);

            inptr += nqBlocks;

            // de-interleave and store data
            for (int d = 0; d < outdim; ++d)
            {
                deinterleave_store(tmpOut[d], nqTot, out_ptr[d]);
                out_ptr[d] += nqBlocks;
            }
        }
    }

    // Size based template version.
    template <int nq0, int nq1>
    void operator2D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, Array<OneD, NekDouble>> &output)
    {
        ASSERTL1(output.size() > 1, "PhysDerivTemplate::Operator2D: Cannot "
                                    "call 2D routine with one output.");
        ASSERTL1(output.size() <= 3, "PhysDerivTemplate::Operator2D: Operator "
                                     "not set up for 3D coordinates.");

        constexpr auto nqTot    = nq0 * nq1;
        constexpr auto nqBlocks = nqTot * vec_t::width;

        const auto outdim        = output.size();
        constexpr int max_outdim = 3;

        const auto ndf        = 2 * outdim;
        constexpr int max_ndf = 2 * max_outdim;

        auto *inptr = &input[0];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut[max_ndf];
        NekDouble *out_ptr[max_outdim]; // max_outdim is a constexpr

        for (int d = 0; d < outdim; ++d)
        {
            tmpOut[d].resize(nqTot);
            out_ptr[d] = &output[d][0];
        }

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        // Workspace for kernels - also checks preconditions
        // PhysDeriv2DWorkspace<SHAPE_TYPE>(nq0, nq1);

        const vec_t *df_ptr;
        vec_t df_tmp[max_ndf]; // max_ndf is a constexpr

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            df_ptr = &((*this->m_df)[dfSize * e]);

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            PhysDeriv2DKernel(nq0, nq1, outdim, tmpIn, this->m_Z[0],
                              this->m_Z[1], this->m_D[0], this->m_D[1], df_ptr,
                              df_tmp, tmpOut);

            inptr += nqBlocks;

            // de-interleave and store data
            for (int d = 0; d < outdim; ++d)
            {
                deinterleave_store(tmpOut[d], nqTot, out_ptr[d]);
                out_ptr[d] += nqBlocks;
            }
        }
    }

#elif defined(SHAPE_DIMENSION_3D)

    // Non-size based operator.
    void operator3D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, Array<OneD, NekDouble>> &output)
    {
        ASSERTL1(output.size() == 3, "PhysDerivTemplate::Operator3D: Cannot "
                                     "call 3D routine with 1 or 2 outputs.");

        const auto nq0 = m_basis[0]->GetNumPoints();
        const auto nq1 = m_basis[1]->GetNumPoints();
        const auto nq2 = m_basis[2]->GetNumPoints();

        const auto nqTot    = nq0 * nq1 * nq2;
        const auto nqBlocks = nqTot * vec_t::width;

        constexpr auto ndf = 9;

        const auto *inptr = input.data();
        auto *outptr_d0   = &output[0][0];
        auto *outptr_d1   = &output[1][0];
        auto *outptr_d2   = &output[2][0];

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0, wsp1Size = 0;
        PhysDeriv3DWorkspace<SHAPE_TYPE>(nq0, nq1, nq2, wsp0Size, wsp1Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), wsp1(wsp1Size),
            tmpIn(nqTot), tmpd0(nqTot), tmpd1(nqTot), tmpd2(nqTot);

        const vec_t *df_ptr;
        vec_t df_tmp[ndf]; // ndf is a constexpr

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            df_ptr = &((*this->m_df)[dfSize * e]);

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            PhysDeriv3DKernel(nq0, nq1, nq2, tmpIn, this->m_Z[0], this->m_Z[1],
                              this->m_Z[2], this->m_D[0], this->m_D[1],
                              this->m_D[2], df_ptr, df_tmp, wsp0, wsp1, tmpd0,
                              tmpd1, tmpd2);

            // de-interleave and store data
            deinterleave_store(tmpd0, nqTot, outptr_d0);
            deinterleave_store(tmpd1, nqTot, outptr_d1);
            deinterleave_store(tmpd2, nqTot, outptr_d2);

            inptr += nqBlocks;
            outptr_d0 += nqBlocks;
            outptr_d1 += nqBlocks;
            outptr_d2 += nqBlocks;
        }
    }

    // Size based template version.
    template <int nq0, int nq1, int nq2>
    void operator3D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, Array<OneD, NekDouble>> &output)
    {
        ASSERTL1(output.size() == 3, "PhysDerivTemplate::Operator3D: Cannot "
                                     "call 3D routine with 1 or 2 outputs.");

        constexpr auto nqTot    = nq0 * nq1 * nq2;
        constexpr auto nqBlocks = nqTot * vec_t::width;

        constexpr auto ndf = 9;

        const auto *inptr = input.data();
        auto *outptr_d0   = &output[0][0];
        auto *outptr_d1   = &output[1][0];
        auto *outptr_d2   = &output[2][0];

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0, wsp1Size = 0;
        PhysDeriv3DWorkspace<SHAPE_TYPE>(nq0, nq1, nq2, wsp0Size, wsp1Size);

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpd0(nqTot),
            tmpd1(nqTot), tmpd2(nqTot), wsp0(wsp0Size), wsp1(wsp1Size);

        const vec_t *df_ptr;
        vec_t df_tmp[ndf]; // ndf is a constexpr

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            df_ptr = &((*this->m_df)[dfSize * e]);

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            PhysDeriv3DKernel(nq0, nq1, nq2, tmpIn, this->m_Z[0], this->m_Z[1],
                              this->m_Z[2], this->m_D[0], this->m_D[1],
                              this->m_D[2], df_ptr, df_tmp, wsp0, wsp1, tmpd0,
                              tmpd1, tmpd2);

            // de-interleave and store data
            deinterleave_store(tmpd0, nqTot, outptr_d0);
            deinterleave_store(tmpd1, nqTot, outptr_d1);
            deinterleave_store(tmpd2, nqTot, outptr_d2);

            inptr += nqBlocks;
            outptr_d0 += nqBlocks;
            outptr_d1 += nqBlocks;
            outptr_d2 += nqBlocks;
        }
    }

#endif

#if defined(SHAPE_DIMENSION_1D)

    // This code is the kernel and has shaped based conditionals. For
    // consistancy it should be in the PhysDerivKernel.hpp
    // code. However, currently it is only used in this class. So for
    // simplicity it is part of the class.
    NEK_FORCE_INLINE void PhysDeriv1DKernel(
        const int nq0, const size_t ndf,
        const std::vector<vec_t, allocator<vec_t>> &in,
        const std::vector<vec_t, allocator<vec_t>> &D0, const vec_t *df_ptr,
        vec_t *df_tmp, std::vector<vec_t, allocator<vec_t>> *out)
    {
        // Get the basic derivative
        PhysDerivTensor1DKernel(nq0, in, D0, out[0]);

        if (!DEFORMED)
        {
            // if( ndf >= 1 )
            df_tmp[0] = df_ptr[0];
            if (ndf >= 2)
                df_tmp[1] = df_ptr[1];
            if (ndf == 3)
                df_tmp[2] = df_ptr[2];
        }

        for (int j = 0; j < nq0; ++j)
        {
            if (DEFORMED)
            {
                // if( ndf >= 1 )
                df_tmp[0] = df_ptr[j * ndf]; // load 1x
                if (ndf >= 2)
                    df_tmp[1] = df_ptr[j * ndf + 1]; // load 1x
                if (ndf == 3)
                    df_tmp[2] = df_ptr[j * ndf + 2]; // load 1x
            }

            // Multiply by derivative factors
            if (ndf == 3)
                out[2][j] = out[0][j] * df_tmp[2]; // Store 1x
            if (ndf >= 2)
                out[1][j] = out[0][j] * df_tmp[1]; // Store 1x
            out[0][j] *= df_tmp[0];                // Store 1x
        }
    }

#elif defined(SHAPE_DIMENSION_2D)

    // This code is the kernel and has shaped based conditionals. For
    // consistancy it should be in the PhysDerivKernel.hpp
    // code. However, currently it is only used in this class. So for
    // simplicity it is part of the class.
    NEK_FORCE_INLINE void PhysDeriv2DKernel(
        const int nq0, const int nq1, const size_t outdim,
        const std::vector<vec_t, allocator<vec_t>> &in,
        const std::vector<vec_t, allocator<vec_t>> &Z0,
        const std::vector<vec_t, allocator<vec_t>> &Z1,
        const std::vector<vec_t, allocator<vec_t>> &D0,
        const std::vector<vec_t, allocator<vec_t>> &D1, const vec_t *df_ptr,
        vec_t *df_tmp, std::vector<vec_t, allocator<vec_t>> *out)
    {
        boost::ignore_unused(Z0, Z1);

        // Results written to out_d0, out_d1
        PhysDerivTensor2DKernel(nq0, nq1, in, D0, D1, out[0], out[1]);

        auto ndf = 2 * outdim;

        if (!DEFORMED)
        {
            // if( outdim >= 2 )
            {
                df_tmp[0] = df_ptr[0];
                df_tmp[1] = df_ptr[1];
                df_tmp[2] = df_ptr[2];
                df_tmp[3] = df_ptr[3];
            }

            if (outdim == 3)
            {
                df_tmp[4] = df_ptr[4];
                df_tmp[5] = df_ptr[5];
            }
        }

        for (int j = 0, cnt_ji = 0; j < nq1; ++j)
        {
#if defined(SHAPE_TYPE_TRI)
            vec_t xfrm0 = 2.0 / (1.0 - Z1[j]); // Load 1x
#endif

            for (int i = 0; i < nq0; ++i, ++cnt_ji)
            {
                vec_t d0 = out[0][cnt_ji]; // Load 1x
                vec_t d1 = out[1][cnt_ji]; // Load 1x

#if defined(SHAPE_TYPE_TRI)
                {
                    // Moving from standard to collapsed coordinates
                    vec_t xfrm1 = 0.5 * (1.0 + Z0[i]); // Load 1x

                    d0 *= xfrm0;
                    d1.fma(d0, xfrm1);
                }
#elif defined(SHAPE_TYPE_QUAD)
                // Nothing to do.
#endif

                if (DEFORMED)
                {
                    // if( outdim >= 2 )
                    {
                        df_tmp[0] = df_ptr[cnt_ji * ndf];
                        df_tmp[1] = df_ptr[cnt_ji * ndf + 1];
                        df_tmp[2] = df_ptr[cnt_ji * ndf + 2];
                        df_tmp[3] = df_ptr[cnt_ji * ndf + 3];
                    }
                    if (outdim == 3)
                    {
                        df_tmp[4] = df_ptr[cnt_ji * ndf + 4];
                        df_tmp[5] = df_ptr[cnt_ji * ndf + 5];
                    }
                }

                // Multiply by derivative factors
                vec_t out0 = d0 * df_tmp[0]; // d0 * df0 + d1 * df1
                out0.fma(d1, df_tmp[1]);
                out[0][cnt_ji] = out0; // Store 1x

                vec_t out1 = d0 * df_tmp[2]; // d0 * df2 + d1 * df3
                out1.fma(d1, df_tmp[3]);
                out[1][cnt_ji] = out1; // Store 1x

                if (outdim == 3)
                {
                    vec_t out2 = d0 * df_tmp[4]; // d0 * df4 + d1 * df5
                    out2.fma(d1, df_tmp[5]);
                    out[2][cnt_ji] = out2; // Store 1x
                }
            }
        }
    }

#elif defined(SHAPE_DIMENSION_3D)

    // This code is the kernel and has shaped based conditionals. For
    // consistancy it should be in the PhysDerivKernel.hpp
    // code. However, currently it is only used herein this class. So
    // for simplicity it is part of the class.
    NEK_FORCE_INLINE void PhysDeriv3DKernel(
        const int nq0, const int nq1, const int nq2,
        const std::vector<vec_t, allocator<vec_t>> &in,
        const std::vector<vec_t, allocator<vec_t>> &Z0,
        const std::vector<vec_t, allocator<vec_t>> &Z1,
        const std::vector<vec_t, allocator<vec_t>> &Z2,
        const std::vector<vec_t, allocator<vec_t>> &D0,
        const std::vector<vec_t, allocator<vec_t>> &D1,
        const std::vector<vec_t, allocator<vec_t>> &D2, const vec_t *df_ptr,
        vec_t *df_tmp, std::vector<vec_t, allocator<vec_t>> &wsp0, // Tets only
        std::vector<vec_t, allocator<vec_t>> &wsp1,                // Tets only
        std::vector<vec_t, allocator<vec_t>> &out_d0,
        std::vector<vec_t, allocator<vec_t>> &out_d1,
        std::vector<vec_t, allocator<vec_t>> &out_d2)
    {
#ifndef SHAPE_TYPE_TET
        boost::ignore_unused(Z0, Z1, Z2, wsp0, wsp1);
#endif
        // Results written to out_d0, out_d1, out_d2
        PhysDerivTensor3DKernel(nq0, nq1, nq2, in, D0, D1, D2, out_d0, out_d1,
                                out_d2);

        constexpr auto ndf = 9;

#if defined(SHAPE_TYPE_TET)
        // Tets get special handling.
        {
            for (int k = 0, eta0 = 0; k < nq2; ++k)
            {
                vec_t xfrm_eta2 = 2.0 / (1.0 - Z2[k]); // Load 1x

                for (int j = 0; j < nq1; ++j)
                {
                    vec_t xfrm_eta1 = 2.0 / (1.0 - Z1[j]); // Load 1x
                    vec_t xfrm      = xfrm_eta1 * xfrm_eta2;

                    for (int i = 0; i < nq0; ++i, ++eta0)
                    {
                        vec_t d0 = xfrm * out_d0[eta0]; // Load 1x

                        out_d0[eta0] = d0; // Store 1x
                        wsp0[eta0]   = d0; // Store 1x partial form for reuse
                    }
                }
            }

            for (int k = 0, eta0 = 0; k < nq2; ++k)
            {
                vec_t xfrm_eta2 = 2.0 / (1.0 - Z2[k]); // Load 1x

                for (int j = 0; j < nq1; ++j)
                {
                    for (int i = 0; i < nq0; ++i, ++eta0)
                    {
                        vec_t xfrm_eta0 = 0.5 * (1.0 + Z0[i]); // Load 1x

                        vec_t out0 = xfrm_eta0 * wsp0[eta0]; // Load 1x
                        wsp0[eta0] = out0; // 2 * (1 + eta_0) / (1 -
                                           // eta_1)(1-eta2) | store 1x

                        vec_t d1 = out_d1[eta0]; // Load 1x
                        d1       = xfrm_eta2 * d1;

                        out_d1[eta0] = out0 + d1; // Store 1x
                        wsp1[eta0]   = d1;        // Store 1x
                    }
                }
            }

            for (int k = 0, eta0 = 0; k < nq2; ++k)
            {
                for (int j = 0; j < nq1; ++j)
                {
                    vec_t xfrm_eta1 = 0.5 * (1.0 + Z1[j]); // Load 1x

                    for (int i = 0; i < nq0; ++i, ++eta0)
                    {
                        vec_t out = wsp0[eta0]; // Load 1x
                        vec_t d1  = wsp1[eta0]; // Load 1x
                        out.fma(d1, xfrm_eta1);
                        out          = out + out_d2[eta0]; // Load 1x
                        out_d2[eta0] = out;                // Store 1x
                    }
                }
            }
        }
#endif
        if (!DEFORMED)
        {
            df_tmp[0] = df_ptr[0];
            df_tmp[1] = df_ptr[1];
            df_tmp[2] = df_ptr[2];
            df_tmp[3] = df_ptr[3];
            df_tmp[4] = df_ptr[4];
            df_tmp[5] = df_ptr[5];
            df_tmp[6] = df_ptr[6];
            df_tmp[7] = df_ptr[7];
            df_tmp[8] = df_ptr[8];
        }

        for (int k = 0, cnt_ijk = 0; k < nq2; ++k)
        {
#if defined(SHAPE_TYPE_PRISM) || defined(SHAPE_TYPE_PYR)
            vec_t xfrm_eta2 = 2.0 / (1.0 - Z2[k]); // Load 1x
#endif

            for (int j = 0; j < nq1; ++j)
            {
#if defined(SHAPE_TYPE_PYR)
                vec_t xfrm_eta1 = 0.5 * (1.0 + Z1[j]); // Load 1x
#endif
                for (int i = 0; i < nq0; ++i, ++cnt_ijk)
                {
                    vec_t d0 = out_d0[cnt_ijk]; // Load 1x
                    vec_t d1 = out_d1[cnt_ijk]; // Load 1x
                    vec_t d2 = out_d2[cnt_ijk]; // Load 1x

#if defined(SHAPE_TYPE_HEX) || defined(SHAPE_TYPE_TET)
                    //     Nothing to do.
#elif defined(SHAPE_TYPE_PRISM)
                    {
                        // Chain-rule for eta_0 and eta_2
                        d0 *= xfrm_eta2; // Load 1x

                        vec_t xfrm_eta0 = 0.5 * (1.0 + Z0[i]); // Load 1x
                        d2.fma(xfrm_eta0, d0);
                    }
#elif defined(SHAPE_TYPE_PYR)
                    {
                        // Chain-rule for eta_0 and eta_2
                        d0 *= xfrm_eta2; // Load 1x
                        d1 *= xfrm_eta2; // Load 1x

                        vec_t xfrm_eta0 = 0.5 * (1.0 + Z0[i]); // Load 1x
                        d2.fma(xfrm_eta0, d0);
                        d2.fma(xfrm_eta1, d1);
                    }
#endif

                    if (DEFORMED)
                    {
                        df_tmp[0] = df_ptr[cnt_ijk * ndf];
                        df_tmp[1] = df_ptr[cnt_ijk * ndf + 1];
                        df_tmp[2] = df_ptr[cnt_ijk * ndf + 2];
                        df_tmp[3] = df_ptr[cnt_ijk * ndf + 3];
                        df_tmp[4] = df_ptr[cnt_ijk * ndf + 4];
                        df_tmp[5] = df_ptr[cnt_ijk * ndf + 5];
                        df_tmp[6] = df_ptr[cnt_ijk * ndf + 6];
                        df_tmp[7] = df_ptr[cnt_ijk * ndf + 7];
                        df_tmp[8] = df_ptr[cnt_ijk * ndf + 8];
                    }

                    // Metric for eta_0, xi_1, eta_2
                    vec_t out0 = d0 * df_tmp[0];
                    out0.fma(d1, df_tmp[1]);
                    out0.fma(d2, df_tmp[2]);
                    out_d0[cnt_ijk] = out0; // Store 1x

                    vec_t out1 = d0 * df_tmp[3];
                    out1.fma(d1, df_tmp[4]);
                    out1.fma(d2, df_tmp[5]);
                    out_d1[cnt_ijk] = out1; // Store 1x

                    vec_t out2 = d0 * df_tmp[6];
                    out2.fma(d1, df_tmp[7]);
                    out2.fma(d2, df_tmp[8]);
                    out_d2[cnt_ijk] = out2; // Store 1x
                }
            }
        }
    }

#endif // SHAPE_DIMENSION

private:
    int m_nmTot;
};

} // namespace MatrixFree
} // namespace Nektar

#endif
