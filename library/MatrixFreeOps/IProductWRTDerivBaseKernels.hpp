#ifndef NEKTAR_LIBRARY_MF_IPRODUCTWRTDERIVBASE_KERNELS_H
#define NEKTAR_LIBRARY_MF_IPRODUCTWRTDERIVBASE_KERNELS_H

namespace Nektar
{
namespace MatrixFree
{

// The dimension and shape kernels. NOTE: They are NOT duplicate
// templated version based on the array size like the
// operators. HOWEVER, they are forced to be INLINED. The inlining is
// critical so that when used in the templated version of the operator
// that loop unrolling occurs.

// The seven shape kernels where the work gets done.
#if defined(SHAPE_TYPE_SEG)

template <bool DEFORMED>
NEK_FORCE_INLINE static void IProductWRTDerivBaseSegKernel(
    const size_t nq0, const size_t ndf, const vec_t *df_ptr, vec_t *df_tmp,
    const std::vector<vec_t, allocator<vec_t>> *tmpIn,
    std::vector<vec_t, allocator<vec_t>> &tmp0)
{
    // Calculate dxi/dx in[0] + dxi/dy in[1] + dxi/dz in[2]
    if (!DEFORMED)
    {
        for (int d = 0; d < ndf; ++d)
            df_tmp[d] = df_ptr[d];
    }

    for (int i = 0; i < nq0; ++i)
    {
        if (DEFORMED)
        {
            for (int d = 0; d < ndf; ++d)
                df_tmp[d] = df_ptr[i * ndf + d];
        }

        tmp0[i] = df_tmp[0] * tmpIn[0][i];
        for (int d = 1; d < ndf; ++d)
            tmp0[i] += (df_tmp[d] * tmpIn[d][i]);
    }
}

#elif defined(SHAPE_TYPE_TRI)

template <bool DEFORMED>
NEK_FORCE_INLINE static void IProductWRTDerivBaseTriKernel(
    const size_t nq0, const size_t nq1, const size_t indim, const vec_t *df_ptr,
    vec_t *df_tmp, const std::vector<vec_t, allocator<vec_t>> &Z0,
    const std::vector<vec_t, allocator<vec_t>> &Z1,
    const std::vector<vec_t, allocator<vec_t>> *tmpIn,
    std::vector<vec_t, allocator<vec_t>> &tmp0,
    std::vector<vec_t, allocator<vec_t>> &tmp1)
{
    const auto ndf = 2 * indim;

    // Calculate dxi/dx in[0] + dxi/dy in[1] + dxi/dz in[2]
    if (!DEFORMED)
    {
        for (int d = 0; d < ndf; ++d)
        {
            df_tmp[d] = df_ptr[d];
        }
    }

    size_t cnt_ji = 0;
    for (size_t j = 0; j < nq1; ++j)
    {
        vec_t f0 = 2.0 / (1.0 - Z1[j]); // Load 1x

        for (size_t i = 0; i < nq0; ++i, ++cnt_ji)
        {
            if (DEFORMED)
            {
                for (int d = 0; d < ndf; ++d)
                {
                    df_tmp[d] = df_ptr[cnt_ji * ndf + d];
                }
            }

            // Calculate dx/dxi in[0] + dy/dxi in[1]
            vec_t tI0 = tmpIn[0][cnt_ji]; // load 1x
            vec_t tI1 = tmpIn[1][cnt_ji]; // load 1x

            vec_t t0 = df_tmp[0] * tI0 + df_tmp[2] * tI1;
            vec_t t1 = df_tmp[1] * tI0 + df_tmp[3] * tI1;

            if (indim == 3)
            {
                vec_t tI2 = tmpIn[2][cnt_ji]; // load 1x

                t0 += df_tmp[4] * tI2;
                t1 += df_tmp[5] * tI2;
            }

            // Multiply by geometric factors
            vec_t hf1 = 0.5 * (1.0 + Z0[i]); // Load 1x
            // Scale by geometric factor 2/(1-z1)
            t0 *= f0;
            // Scale by geometric factor (1+z0)/(1-z1)
            vec_t c1 = hf1 * t1;
            t0.fma(c1, f0);

            // Store
            tmp0[cnt_ji] = t0; // store 1x
            tmp1[cnt_ji] = t1; // store 1x
        }
    }
}

#elif defined(SHAPE_TYPE_QUAD)

template <bool DEFORMED>
NEK_FORCE_INLINE static void IProductWRTDerivBaseQuadKernel(
    const size_t nq0, const size_t nq1, const size_t indim, const vec_t *df_ptr,
    vec_t *df_tmp, const std::vector<vec_t, allocator<vec_t>> *tmpIn,
    std::vector<vec_t, allocator<vec_t>> &tmp0,
    std::vector<vec_t, allocator<vec_t>> &tmp1)
{
    const auto ndf   = 2 * indim;
    const auto nqTot = nq0 * nq1;

    // Calculate dxi/dx in[0] + dxi/dy in[1] + dxi/dz in[2]
    if (!DEFORMED)
    {
        for (int d = 0; d < ndf; ++d)
        {
            df_tmp[d] = df_ptr[d];
        }
    }

    for (int i = 0; i < nqTot; ++i)
    {
        if (DEFORMED)
        {
            for (int d = 0; d < ndf; ++d)
            {
                df_tmp[d] = df_ptr[i * ndf + d];
            }
        }

        tmp0[i] = df_tmp[0] * tmpIn[0][i] + df_tmp[2] * tmpIn[1][i];
        tmp1[i] = df_tmp[1] * tmpIn[0][i] + df_tmp[3] * tmpIn[1][i];

        if (indim == 3)
        {
            tmp0[i] += df_tmp[4] * tmpIn[2][i];
            tmp1[i] += df_tmp[5] * tmpIn[2][i];
        }
    }
}

#elif defined(SHAPE_TYPE_HEX)

template <bool DEFORMED>
NEK_FORCE_INLINE static void IProductWRTDerivBaseHexKernel(
    const size_t nq0, const size_t nq1, const size_t nq2, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn0,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn1,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn2,
    std::vector<vec_t, allocator<vec_t>> &tmp0,
    std::vector<vec_t, allocator<vec_t>> &tmp1,
    std::vector<vec_t, allocator<vec_t>> &tmp2)
{
    constexpr auto ndf = 9u;
    const auto nqTot   = nq0 * nq1 * nq2;

    // Calculate dx/dxi in[0] + dy/dxi in[1] + dz/dxi in[2]
    vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;

    if (!DEFORMED)
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
        df4 = df_ptr[4];
        df5 = df_ptr[5];
        df6 = df_ptr[6];
        df7 = df_ptr[7];
        df8 = df_ptr[8];
    }

    for (int i = 0; i < nqTot; ++i)
    {
        if (DEFORMED)
        {
            df0 = df_ptr[i * ndf];
            df1 = df_ptr[i * ndf + 1];
            df2 = df_ptr[i * ndf + 2];
            df3 = df_ptr[i * ndf + 3];
            df4 = df_ptr[i * ndf + 4];
            df5 = df_ptr[i * ndf + 5];
            df6 = df_ptr[i * ndf + 6];
            df7 = df_ptr[i * ndf + 7];
            df8 = df_ptr[i * ndf + 8];
        }

        vec_t tI0 = tmpIn0[i]; // load 1x
        vec_t tI1 = tmpIn1[i]; // load 1x
        vec_t tI2 = tmpIn2[i]; // load 1x
        tmp0[i]   = df0 * tI0 + df3 * tI1 + df6 * tI2;
        tmp1[i]   = df1 * tI0 + df4 * tI1 + df7 * tI2;
        tmp2[i]   = df2 * tI0 + df5 * tI1 + df8 * tI2;
    }
}

#elif defined(SHAPE_TYPE_TET)

template <bool DEFORMED>
NEK_FORCE_INLINE static void IProductWRTDerivBaseTetKernel(
    const size_t nq0, const size_t nq1, const size_t nq2, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &Z0,
    const std::vector<vec_t, allocator<vec_t>> &Z1,
    const std::vector<vec_t, allocator<vec_t>> &Z2,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn0,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn1,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn2,
    std::vector<vec_t, allocator<vec_t>> &tmp0,
    std::vector<vec_t, allocator<vec_t>> &tmp1,
    std::vector<vec_t, allocator<vec_t>> &tmp2)
{
    constexpr auto ndf = 9u;

    // Calculate dx/dxi in[0] + dy/dxi in[1] + dz/dxi in[2]
    vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;

    if (!DEFORMED)
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
        df4 = df_ptr[4];
        df5 = df_ptr[5];
        df6 = df_ptr[6];
        df7 = df_ptr[7];
        df8 = df_ptr[8];
    }

    for (size_t k = 0, cnt_kji = 0; k < nq2; ++k)
    {
        vec_t f2 = 2.0 / (1.0 - Z2[k]); // Load 1x

        for (size_t j = 0; j < nq1; ++j)
        {
            vec_t Z1Load = Z1[j]; // Load1x
            vec_t f3     = 0.5 * (1.0 + Z1Load);
            vec_t f0     = 2.0 * f2 / (1.0 - Z1Load);

            for (size_t i = 0; i < nq0; ++i, ++cnt_kji)
            {
                if (DEFORMED)
                {
                    df0 = df_ptr[cnt_kji * ndf];
                    df1 = df_ptr[cnt_kji * ndf + 1];
                    df2 = df_ptr[cnt_kji * ndf + 2];
                    df3 = df_ptr[cnt_kji * ndf + 3];
                    df4 = df_ptr[cnt_kji * ndf + 4];
                    df5 = df_ptr[cnt_kji * ndf + 5];
                    df6 = df_ptr[cnt_kji * ndf + 6];
                    df7 = df_ptr[cnt_kji * ndf + 7];
                    df8 = df_ptr[cnt_kji * ndf + 8];
                }

                // Calculate dx/dxi in[0] + dy/dxi in[1] + dz/dxi in[2]
                vec_t tI0 = tmpIn0[cnt_kji]; // load 1x
                vec_t tI1 = tmpIn1[cnt_kji]; // load 1x
                vec_t tI2 = tmpIn2[cnt_kji]; // load 1x
                vec_t t0  = df0 * tI0 + df3 * tI1 + df6 * tI2;
                vec_t t1  = df1 * tI0 + df4 * tI1 + df7 * tI2;
                vec_t t2  = df2 * tI0 + df5 * tI1 + df8 * tI2;

                // Multiply by geometric factors
                vec_t f1 = 0.5 * (1.0 + Z0[i]); // Load 1x

                // Scale by geometric factor 1 and add to t0
                t0.fma(t1 + t2, f1);
                // Scale by geometric factor 0
                t0 *= f0;

                // Scale by geometric factor 3 and add to t1
                t1.fma(t2, f3);
                // Scale by geometric factor 2
                t1 *= f2;

                // Store
                tmp0[cnt_kji] = t0; // store 1x
                tmp1[cnt_kji] = t1; // store 1x
                tmp2[cnt_kji] = t2; // store 1x
            }
        }
    }
}

#elif defined(SHAPE_TYPE_PRISM)

template <bool DEFORMED>
NEK_FORCE_INLINE static void IProductWRTDerivBasePrismKernel(
    const size_t nq0, const size_t nq1, const size_t nq2, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &Z0,
    const std::vector<vec_t, allocator<vec_t>> &Z2,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn0,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn1,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn2,
    std::vector<vec_t, allocator<vec_t>> &tmp0,
    std::vector<vec_t, allocator<vec_t>> &tmp1,
    std::vector<vec_t, allocator<vec_t>> &tmp2)
{
    constexpr auto ndf = 9u;

    // Calculate dx/dxi in[0] + dy/dxi in[1] + dz/dxi in[2]
    vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;

    if (!DEFORMED)
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
        df4 = df_ptr[4];
        df5 = df_ptr[5];
        df6 = df_ptr[6];
        df7 = df_ptr[7];
        df8 = df_ptr[8];
    }

    for (size_t k = 0, cnt_kji = 0; k < nq2; ++k)
    {
        // div in most external loop
        vec_t f0 = 2.0 / (1.0 - Z2[k]); // Load 1x

        for (size_t j = 0; j < nq1; ++j)
        {
            for (size_t i = 0; i < nq0; ++i, ++cnt_kji)
            {
                if (DEFORMED)
                {
                    df0 = df_ptr[cnt_kji * ndf];
                    df1 = df_ptr[cnt_kji * ndf + 1];
                    df2 = df_ptr[cnt_kji * ndf + 2];
                    df3 = df_ptr[cnt_kji * ndf + 3];
                    df4 = df_ptr[cnt_kji * ndf + 4];
                    df5 = df_ptr[cnt_kji * ndf + 5];
                    df6 = df_ptr[cnt_kji * ndf + 6];
                    df7 = df_ptr[cnt_kji * ndf + 7];
                    df8 = df_ptr[cnt_kji * ndf + 8];
                }

                // Calculate dx/dxi in[0] + dy/dxi in[1] + dz/dxi in[2]
                vec_t tI0 = tmpIn0[cnt_kji]; // load 1x
                vec_t tI1 = tmpIn1[cnt_kji]; // load 1x
                vec_t tI2 = tmpIn2[cnt_kji]; // load 1x
                vec_t t0  = df0 * tI0 + df3 * tI1 + df6 * tI2;
                vec_t t1  = df1 * tI0 + df4 * tI1 + df7 * tI2;
                vec_t t2  = df2 * tI0 + df5 * tI1 + df8 * tI2;

                // Multiply by geometric factors
                vec_t hf1 = 0.5 * (1.0 + Z0[i]); // Load 1x
                // Scale by geometric factor 2/(1-z2)
                t0 *= f0;
                // Scale by geometric factor (1+z0)/(1-z2)
                vec_t f1t2 = hf1 * t2;
                t0.fma(f1t2, f0);

                // Store
                tmp0[cnt_kji] = t0; // store 1x
                tmp1[cnt_kji] = t1; // store 1x
                tmp2[cnt_kji] = t2; // store 1x
            }
        }
    }
}

#elif defined(SHAPE_TYPE_PYR)

template <bool DEFORMED>
NEK_FORCE_INLINE static void IProductWRTDerivBasePyrKernel(
    const size_t nq0, const size_t nq1, const size_t nq2, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &Z0,
    const std::vector<vec_t, allocator<vec_t>> &Z1,
    const std::vector<vec_t, allocator<vec_t>> &Z2,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn0,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn1,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn2,
    std::vector<vec_t, allocator<vec_t>> &tmp0,
    std::vector<vec_t, allocator<vec_t>> &tmp1,
    std::vector<vec_t, allocator<vec_t>> &tmp2)
{
    constexpr auto ndf = 9u;

    // Calculate dx/dxi in[0] + dy/dxi in[1] + dz/dxi in[2]
    vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;

    if (!DEFORMED)
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
        df4 = df_ptr[4];
        df5 = df_ptr[5];
        df6 = df_ptr[6];
        df7 = df_ptr[7];
        df8 = df_ptr[8];
    }

    for (size_t k = 0, cnt_kji = 0; k < nq2; ++k)
    {
        // div in most external loop
        vec_t f0 = 2.0 / (1.0 - Z2[k]); // Load 1x

        for (size_t j = 0; j < nq1; ++j)
        {
            vec_t hf2 = 0.5 * (1.0 + Z1[j]); // Load 1x

            for (size_t i = 0; i < nq0; ++i, ++cnt_kji)
            {
                if (DEFORMED)
                {
                    df0 = df_ptr[cnt_kji * ndf];
                    df1 = df_ptr[cnt_kji * ndf + 1];
                    df2 = df_ptr[cnt_kji * ndf + 2];
                    df3 = df_ptr[cnt_kji * ndf + 3];
                    df4 = df_ptr[cnt_kji * ndf + 4];
                    df5 = df_ptr[cnt_kji * ndf + 5];
                    df6 = df_ptr[cnt_kji * ndf + 6];
                    df7 = df_ptr[cnt_kji * ndf + 7];
                    df8 = df_ptr[cnt_kji * ndf + 8];
                }

                // Calculate dx/dxi in[0] + dy/dxi in[1] + dz/dxi in[2]
                vec_t tI0 = tmpIn0[cnt_kji]; // load 1x
                vec_t tI1 = tmpIn1[cnt_kji]; // load 1x
                vec_t tI2 = tmpIn2[cnt_kji]; // load 1x
                vec_t t0  = df0 * tI0 + df3 * tI1 + df6 * tI2;
                vec_t t1  = df1 * tI0 + df4 * tI1 + df7 * tI2;
                vec_t t2  = df2 * tI0 + df5 * tI1 + df8 * tI2;

                // Scale by geometric factor 2/(1-z2)
                t0 *= f0;
                vec_t hf1 = 0.5 * (1.0 + Z0[i]); // Load 1x
                // Scale by geometric factor (1+z0)/(1-z2)
                vec_t f1t2 = hf1 * t2;
                t0.fma(f1t2, f0);

                // Scale by geometric factor 2/(1-z2)
                t1 *= f0;
                // Scale by geometric factor (1+z1)/(1-z2)
                f1t2 = hf2 * t2;
                t1.fma(f1t2, f0);

                // Store
                tmp0[cnt_kji] = t0; // store 1x
                tmp1[cnt_kji] = t1; // store 1x
                tmp2[cnt_kji] = t2; // store 1x
            }
        }
    }
}

#endif // SHAPE_TYPE

// The dimension kernels which select the shape kernel.
#if defined(SHAPE_DIMENSION_1D)

template <LibUtilities::ShapeType SHAPE_TYPE, bool DEFORMED>
NEK_FORCE_INLINE static void IProductWRTDerivBase1DKernel(
    const size_t nq0, const size_t ndf, const vec_t *df_ptr, vec_t *df_tmp,
    const std::vector<vec_t, allocator<vec_t>> *tmpIn,
    std::vector<vec_t, allocator<vec_t>> &tmp0)
{
#if defined(SHAPE_TYPE_SEG)
    IProductWRTDerivBaseSegKernel<DEFORMED>(nq0, ndf, df_ptr, df_tmp, tmpIn,
                                            tmp0);
#endif
}

#elif defined(SHAPE_DIMENSION_2D)

template <LibUtilities::ShapeType SHAPE_TYPE, bool DEFORMED>
NEK_FORCE_INLINE static void IProductWRTDerivBase2DKernel(
    const size_t nq0, const size_t nq1, const size_t indim, const vec_t *df_ptr,
    vec_t *df_tmp, const std::vector<vec_t, allocator<vec_t>> &Z0,
    const std::vector<vec_t, allocator<vec_t>> &Z1,
    const std::vector<vec_t, allocator<vec_t>> *tmpIn,
    std::vector<vec_t, allocator<vec_t>> &tmp0,
    std::vector<vec_t, allocator<vec_t>> &tmp1)
{
#if defined(SHAPE_TYPE_TRI)
    IProductWRTDerivBaseTriKernel<DEFORMED>(nq0, nq1, indim, df_ptr, df_tmp, Z0,
                                            Z1, tmpIn, tmp0, tmp1);
#elif defined(SHAPE_TYPE_QUAD)
    boost::ignore_unused(Z0, Z1);
    IProductWRTDerivBaseQuadKernel<DEFORMED>(nq0, nq1, indim, df_ptr, df_tmp,
                                             tmpIn, tmp0, tmp1);
#endif
}

#elif defined(SHAPE_DIMENSION_3D)

template <LibUtilities::ShapeType SHAPE_TYPE, bool DEFORMED>
NEK_FORCE_INLINE static void IProductWRTDerivBase3DKernel(
    const size_t nq0, const size_t nq1, const size_t nq2, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &Z0,
    const std::vector<vec_t, allocator<vec_t>> &Z1,
    const std::vector<vec_t, allocator<vec_t>> &Z2,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn0,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn1,
    const std::vector<vec_t, allocator<vec_t>> &tmpIn2,
    std::vector<vec_t, allocator<vec_t>> &tmp0,
    std::vector<vec_t, allocator<vec_t>> &tmp1,
    std::vector<vec_t, allocator<vec_t>> &tmp2)
{
#if defined(SHAPE_TYPE_HEX)
    boost::ignore_unused(Z0, Z1, Z2);
    IProductWRTDerivBaseHexKernel<DEFORMED>(nq0, nq1, nq2, df_ptr, tmpIn0,
                                            tmpIn1, tmpIn2, tmp0, tmp1, tmp2);
#elif defined(SHAPE_TYPE_TET)
    IProductWRTDerivBaseTetKernel<DEFORMED>(nq0, nq1, nq2, df_ptr, Z0, Z1, Z2,
                                            tmpIn0, tmpIn1, tmpIn2, tmp0, tmp1,
                                            tmp2);
#elif defined(SHAPE_TYPE_PRISM)
    boost::ignore_unused(Z1);
    IProductWRTDerivBasePrismKernel<DEFORMED>(nq0, nq1, nq2, df_ptr, Z0, Z2,
                                              tmpIn0, tmpIn1, tmpIn2, tmp0,
                                              tmp1, tmp2);
#elif defined(SHAPE_TYPE_PYR)
    IProductWRTDerivBasePyrKernel<DEFORMED>(nq0, nq1, nq2, df_ptr, Z0, Z1, Z2,
                                            tmpIn0, tmpIn1, tmpIn2, tmp0, tmp1,
                                            tmp2);
#endif
}

#endif // SHAPE_DIMENSION

} // namespace MatrixFree
} // namespace Nektar

#endif
