///////////////////////////////////////////////////////////////////////////////
//
// File: PhysDerivKernels.hpp
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

#ifndef NEKTAR_LIBRARY_MF_PHYS_DERIV_KERNELS_HPP
#define NEKTAR_LIBRARY_MF_PHYS_DERIV_KERNELS_HPP

namespace Nektar
{
namespace MatrixFree
{

using namespace tinysimd;
using vec_t = simd<NekDouble>;

// The dimension kernels. NOTE: They are NOT duplicate
// templated version based on the array size like the
// operators. HOWEVER, they are forced to be INLINED. The inlining is
// critical so that when used in the templated version of the operator
// that loop unrolling occurs.

// The three shape kernels where the work gets done.
#if defined(SHAPE_DIMENSION_1D)

NEK_FORCE_INLINE static void PhysDerivTensor1DKernel(
    const size_t nq0, const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    std::vector<vec_t, allocator<vec_t>> &out_d0)
{
    // All matricies are column major ordered since operators used to
    // be computed via BLAS.

    // D0 * in
    for (int i = 0; i < nq0; ++i)
    { // Row index of D0 matrix

        vec_t prod_sum = 0.0;
        for (int k = 0; k < nq0; ++k)
        {                               // Col index of D0, row index of IN
            vec_t v1 = D0[k * nq0 + i]; // Load 1x
            vec_t v2 = in[k];           // Load 1x

            prod_sum.fma(v1, v2);
        }

        out_d0[i] = prod_sum; // Store 1x
    }
}

#elif defined(SHAPE_DIMENSION_2D)

NEK_FORCE_INLINE static void PhysDerivTensor2DKernel(
    const size_t nq0, const size_t nq1,
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const std::vector<vec_t, allocator<vec_t>> &D1,
    std::vector<vec_t, allocator<vec_t>> &out_d0,
    std::vector<vec_t, allocator<vec_t>> &out_d1)
{
    // All matricies are column major ordered since operators used to
    // be computed via BLAS.

    // D0 * in
    for (int i = 0; i < nq0; ++i)
    { // Row index of D0 matrix
        for (int j = 0; j < nq1; ++j)
        { // Col index of IN matrix

            vec_t prod_sum = 0.0;
            for (int k = 0; k < nq0; ++k)
            {                               // Col index of D0, row index of IN
                vec_t v1 = D0[k * nq0 + i]; // Load 1x
                vec_t v2 = in[j * nq0 + k]; // Load 1x

                prod_sum.fma(v1, v2);
            }

            out_d0[j * nq0 + i] = prod_sum; // Store 1x
        }
    }

    // in * D1^T
    for (int i = 0; i < nq0; ++i)
    { // row index for grid
        for (int j = 0; j < nq1; ++j)
        { // Column index for D1^T (row idx for D1)

            vec_t prod_sum = 0.0;
            for (int k = 0; k < nq1; ++k)
            {
                vec_t v1 = in[k * nq0 + i]; // Load 1x
                vec_t v2 = D1[k * nq1 + j]; // Load 1x

                prod_sum.fma(v1, v2);
            }

            out_d1[j * nq0 + i] = prod_sum; // Store 1x
        }
    }
}

#elif defined(SHAPE_DIMENSION_3D)

NEK_FORCE_INLINE static void PhysDerivTensor3DKernel(
    const size_t nq0, const size_t nq1, const size_t nq2,
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const std::vector<vec_t, allocator<vec_t>> &D1,
    const std::vector<vec_t, allocator<vec_t>> &D2,
    std::vector<vec_t, allocator<vec_t>> &out_d0,
    std::vector<vec_t, allocator<vec_t>> &out_d1,
    std::vector<vec_t, allocator<vec_t>> &out_d2)
{
    // All matricies are column major ordered since operators used to
    // be computed via BLAS.

    // Direction 1
    for (int i = 0; i < nq0; ++i)
    {
        for (int j = 0; j < nq1 * nq2; ++j)
        {
            vec_t prod_sum = 0.0;
            for (int k = 0; k < nq0; ++k)
            {
                vec_t v1 = D0[k * nq0 + i]; // Load 1x
                vec_t v2 = in[j * nq0 + k]; // Load 1x

                prod_sum.fma(v1, v2);
            }

            out_d0[j * nq0 + i] = prod_sum; // Store 1x
        }
    }

    // Direction 2
    for (int block = 0; block < nq2; ++block)
    {
        int start = block * nq0 * nq1;

        for (int i = 0; i < nq0; ++i)
        {
            for (int j = 0; j < nq1; ++j)
            {
                vec_t prod_sum = 0.0;
                for (int k = 0; k < nq1; ++k)
                {
                    vec_t v1 = in[start + k * nq0 + i]; // Load 1x
                    vec_t v2 = D1[k * nq1 + j];         // Load 1x

                    prod_sum.fma(v1, v2);
                }

                out_d1[start + j * nq0 + i] = prod_sum; // Store 1x
            }
        }
    }

    // Direction 3
    for (int i = 0; i < nq0 * nq1; ++i)
    {
        for (int j = 0; j < nq2; ++j)
        {
            vec_t prod_sum = 0.0;
            for (int k = 0; k < nq2; ++k)
            {
                vec_t v1 = in[k * nq0 * nq1 + i]; // Load 1x
                vec_t v2 = D2[k * nq2 + j];       // Load 1x

                prod_sum.fma(v1, v2);
            }

            out_d2[j * nq0 * nq1 + i] = prod_sum; // Store 1x
        }
    }
}

#endif

// Workspace - used to dynamically get the workspace size needed for
// temporary memory.
#if defined(SHAPE_DIMENSION_1D)

template <LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void PhysDeriv1DWorkspace(const size_t nq0)
{
    boost::ignore_unused(SHAPE_TYPE, nq0);

    // Check preconditions
    // None
}

#elif defined(SHAPE_DIMENSION_2D)

template <LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void PhysDeriv2DWorkspace(const size_t nq0,
                                                  const size_t nq1)
{
    boost::ignore_unused(SHAPE_TYPE, nq0, nq1);

    // Check preconditions
    ASSERTL1((SHAPE_TYPE == LibUtilities::ShapeType::Tri && nq0 == nq1 + 1) ||
                 (SHAPE_TYPE == LibUtilities::ShapeType::Quad && nq0 == nq1),
             "PhysDeriv2DWorkspace: Requires homogenous points.");
}

#elif defined(SHAPE_DIMENSION_3D)

template <LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void PhysDeriv3DWorkspace(const size_t nq0,
                                                  const size_t nq1,
                                                  const size_t nq2,
                                                  size_t &wsp1Size,
                                                  size_t &wsp2Size)
{
    boost::ignore_unused(SHAPE_TYPE, nq0, nq1, nq2, wsp1Size, wsp2Size);

    // Check preconditions
    ASSERTL1((SHAPE_TYPE == LibUtilities::ShapeType::Hex && nq0 == nq1 &&
              nq0 == nq2) ||
                 (SHAPE_TYPE == LibUtilities::ShapeType::Tet &&
                  nq0 == nq1 + 1 && nq0 == nq2 + 1) ||
                 (SHAPE_TYPE == LibUtilities::ShapeType::Prism && nq0 == nq1 &&
                  nq0 == nq2 + 1) ||
                 (SHAPE_TYPE == LibUtilities::ShapeType::Pyr && nq0 == nq1 &&
                  nq0 == nq2 + 1),
             "PhysDeriv3DWorkspace: Requires homogenous points.");

#if defined(SHAPE_TYPE_TET)
    wsp1Size = std::max(wsp1Size, nq0 * nq1 * nq2);
    wsp2Size = std::max(wsp2Size, nq0 * nq1 * nq2);
#endif
}

#endif // SHAPE_DIMENSION

} // namespace MatrixFree
} // namespace Nektar

#endif
