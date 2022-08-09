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

// Overloaded with diffusion
NEK_FORCE_INLINE static void PhysDerivTensor2DKernel(
    const size_t nq0, const size_t nq1,
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const std::vector<vec_t, allocator<vec_t>> &D1,
    std::vector<vec_t, allocator<vec_t>> &out_d0,
    std::vector<vec_t, allocator<vec_t>> &out_d1,
    const Array<OneD, Array<OneD, NekDouble>> &diff)
{
    const auto nqTot = nq0 * nq1;

    vec_t d00 = diff[0][0], d01 = diff[0][1], d11 = diff[1][1];

    std::vector<vec_t, allocator<vec_t>> temp_d0(nqTot), temp_d1(nqTot);

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

            temp_d0[j * nq0 + i] = prod_sum; // Store 1x
        }
    }

    // in * D1^T
    for (int i = 0; i < nq0; ++i)
    { // Row index for grid
        for (int j = 0; j < nq1; ++j)
        { // Column index for D1^T (row idx for D1)

            vec_t prod_sum = 0.0;
            for (int k = 0; k < nq1; ++k)
            {
                vec_t v1 = in[k * nq0 + i]; // Load 1x
                vec_t v2 = D1[k * nq1 + j]; // Load 1x

                prod_sum.fma(v1, v2);
            }

            temp_d1[j * nq0 + i] = prod_sum; // Store 1x
        }
    }

    for (int ij = 0; ij < nqTot; ++ij)
    {
        vec_t prod_sum = 0.0;
        vec_t temp0    = temp_d0[ij];
        prod_sum.fma(d00, temp0);
        vec_t temp1 = temp_d1[ij];
        prod_sum.fma(d01, temp1);
        out_d0[ij] = prod_sum;
    }

    for (int ij = 0; ij < nqTot; ++ij)
    {
        vec_t prod_sum = 0.0;
        vec_t temp0    = temp_d0[ij];
        prod_sum.fma(d01, temp0);
        vec_t temp1 = temp_d1[ij];
        prod_sum.fma(d11, temp1);
        out_d1[ij] = prod_sum;
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
