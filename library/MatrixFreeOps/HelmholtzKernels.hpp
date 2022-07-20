#ifndef NEKTAR_LIBRARY_MF_HELMHOLTZ_KERNELS_H
#define NEKTAR_LIBRARY_MF_HELMHOLTZ_KERNELS_H

namespace Nektar
{
namespace MatrixFree
{

using namespace tinysimd;
using vec_t = simd<NekDouble>;

// The dimension and shape kernels. NOTE: They are NOT duplicate
// templated version based on the array size like the
// operators. HOWEVER, they are forced to be INLINED. The inlining is
// critical so that when used in the templated version of the operator
// that loop unrolling occurs.

// The seven shape kernels where the work gets done.
#if defined(SHAPE_TYPE_TRI)

template <bool DEFORMED>
NEK_FORCE_INLINE static void HelmholtzTriKernel(
    const size_t nq0, const size_t nq1, const bool isConstVarDiff,
    const Array<OneD, NekDouble> &constVarDiff, const bool isVarDiff,
    const Array<OneD, NekDouble> &varD00, const Array<OneD, NekDouble> &varD01,
    const Array<OneD, NekDouble> &varD11, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &m_h0,
    const std::vector<vec_t, allocator<vec_t>> &m_h1,
    std::vector<vec_t, allocator<vec_t>> &bwd,
    std::vector<vec_t, allocator<vec_t>> &deriv0,
    std::vector<vec_t, allocator<vec_t>> &deriv1)
{
    constexpr auto ndf = 4;

    vec_t df0, df1, df2, df3;
    vec_t metric00, metric01, metric11;

    vec_t d00 = {1.0};
    vec_t d01 = {0.0};
    vec_t d11 = {1.0};                // var diffusion terms
    vec_t dtmp0, dtmp1, dtmp2, dtmp3; // metric products

    if (isConstVarDiff)
    {
        d00 = constVarDiff[0];
        d01 = constVarDiff[1];
        d11 = constVarDiff[2];
    }

    // Precompute Laplacian metricsp
    if (!DEFORMED)
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
    }

    // Step 4a: Construct Laplacian metrics
    for (size_t j = 0, cnt = 0; j < nq1; ++j)
    {
        vec_t h1j = m_h1[j];
        for (size_t i = 0; i < nq0; ++i, ++cnt)
        {
            if (DEFORMED)
            {
                df0 = df_ptr[cnt * ndf];
                df1 = df_ptr[cnt * ndf + 1];
                df2 = df_ptr[cnt * ndf + 2];
                df3 = df_ptr[cnt * ndf + 3];
            }

            vec_t h0i = m_h0[i];

            // M = [M_00, df1; M_10; df3]
            metric00  = h1j * (df0 + h0i * df1); // M_00
            vec_t tmp = h1j * (df2 + h0i * df3); // M_10

            if (!isConstVarDiff && !isVarDiff)
            {
                metric01 = metric00 * df1;
                metric00 = metric00 * metric00;

                metric01.fma(tmp, df3);
                metric00.fma(tmp, tmp);

                metric11 = df1 * df1;
                metric11.fma(df3, df3);
            }
            else
            {
                if (isVarDiff)
                {
                    d00 = varD00[cnt];
                    d01 = varD01[cnt];
                    d11 = varD11[cnt];
                }
                // M = [M_00, df1; M_10; df3]
                dtmp0 = metric00 * d00;
                dtmp0.fma(tmp, d01);
                dtmp1 = metric00 * d01;
                dtmp1.fma(tmp, d11);
                dtmp2 = df1 * d00;
                dtmp2.fma(df3, d01);
                dtmp3 = df1 * d01;
                dtmp3.fma(df3, d11);

                metric00 = metric00 * dtmp0;
                metric00.fma(tmp, dtmp1);
                metric01 = df1 * dtmp0;
                metric01.fma(df3, dtmp1);
                metric11 = df1 * dtmp2;
                metric11.fma(df3, dtmp3);
            }

            vec_t d0 = deriv0[cnt];
            vec_t d1 = deriv1[cnt];

            tmp = metric00 * d0;
            tmp.fma(metric01, d1);
            bwd[cnt] = tmp;

            tmp = metric01 * d0;
            tmp.fma(metric11, d1);
            deriv0[cnt] = tmp;
        }
    }
}

#elif defined(SHAPE_TYPE_QUAD)

template <bool DEFORMED>
NEK_FORCE_INLINE static void HelmholtzQuadKernel(
    const size_t nq0, const size_t nq1, const bool isConstVarDiff,
    const Array<OneD, NekDouble> &constVarDiff, const bool isVarDiff,
    const Array<OneD, NekDouble> &varD00, const Array<OneD, NekDouble> &varD01,
    const Array<OneD, NekDouble> &varD11, const vec_t *df_ptr,
    std::vector<vec_t, allocator<vec_t>> &bwd,
    std::vector<vec_t, allocator<vec_t>> &deriv0,
    std::vector<vec_t, allocator<vec_t>> &deriv1)
{
    constexpr auto ndf = 4;

    const auto nqTot = nq0 * nq1;

    vec_t d00 = {1.0};
    vec_t d01 = {0.0};
    vec_t d11 = {1.0};                // var diffusion terms
    vec_t dtmp0, dtmp1, dtmp2, dtmp3; // temp for vardiff
    vec_t df0, df1, df2, df3;
    vec_t metric00, metric01, metric11;

    if (isConstVarDiff)
    {
        d00 = constVarDiff[0];
        d01 = constVarDiff[1];
        d11 = constVarDiff[2];
    }

    // Precompute Laplacian metricsp
    if (!DEFORMED)
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];

        if (!isConstVarDiff && !isVarDiff)
        {
            metric00 = df0 * df0;
            metric00.fma(df2, df2);
            metric01 = df0 * df1;
            metric01.fma(df2, df3);
            metric11 = df1 * df1;
            metric11.fma(df3, df3);
        }
        else if (isConstVarDiff)
        {
            dtmp0 = df0 * d00;
            dtmp0.fma(df2, d01);
            dtmp1 = df0 * d01;
            dtmp1.fma(df2, d11);
            dtmp2 = df1 * d00;
            dtmp2.fma(df3, d01);
            dtmp3 = df1 * d01;
            dtmp3.fma(df3, d11);

            metric00 = df0 * dtmp0;
            metric00.fma(df2, dtmp1);
            metric01 = df1 * dtmp0;
            metric01.fma(df3, dtmp1);
            metric11 = df1 * dtmp2;
            metric11.fma(df3, dtmp3);
        }
    }

    // Step 4: Apply Laplacian metrics & inner product
    if (!isVarDiff)
    {
        if (DEFORMED)
        {
            for (size_t j = 0, cnt = 0; j < nq1; ++j)
            {
                for (size_t i = 0; i < nq0; ++i, ++cnt)
                {
                    df0 = df_ptr[cnt * ndf];
                    df1 = df_ptr[cnt * ndf + 1];
                    df2 = df_ptr[cnt * ndf + 2];
                    df3 = df_ptr[cnt * ndf + 3];

                    if (!isConstVarDiff)
                    {
                        metric00 = df0 * df0;
                        metric00.fma(df2, df2);
                        metric01 = df0 * df1;
                        metric01.fma(df2, df3);
                        metric11 = df1 * df1;
                        metric11.fma(df3, df3);
                    }
                    else
                    {
                        dtmp0 = df0 * d00;
                        dtmp0.fma(df2, d01);
                        dtmp1 = df0 * d01;
                        dtmp1.fma(df2, d11);
                        dtmp2 = df1 * d00;
                        dtmp2.fma(df3, d01);
                        dtmp3 = df1 * d01;
                        dtmp3.fma(df3, d11);

                        metric00 = df0 * dtmp0;
                        metric00.fma(df2, dtmp1);
                        metric01 = df1 * dtmp0;
                        metric01.fma(df3, dtmp1);
                        metric11 = df1 * dtmp2;
                        metric11.fma(df3, dtmp3);
                    }

                    vec_t d0 = deriv0[cnt];
                    vec_t d1 = deriv1[cnt];

                    vec_t tmp = metric00 * d0;
                    tmp.fma(metric01, d1);
                    bwd[cnt] = tmp;

                    tmp = metric01 * d0;
                    tmp.fma(metric11, d1);
                    deriv0[cnt] = tmp;
                }
            }
        }
        else
        {
            for (int i = 0; i < nqTot; ++i)
            {
                vec_t d0 = deriv0[i];
                vec_t d1 = deriv1[i];

                vec_t tmp = metric00 * d0;
                tmp.fma(metric01, d1);
                bwd[i] = tmp;

                tmp = metric01 * d0;
                tmp.fma(metric11, d1);
                deriv0[i] = tmp;
            }
        }
    }
    else
    {
        if (DEFORMED)
        {
            for (size_t j = 0, cnt = 0; j < nq1; ++j)
            {
                for (size_t i = 0; i < nq0; ++i, ++cnt)
                {
                    df0 = df_ptr[cnt * ndf];
                    df1 = df_ptr[cnt * ndf + 1];
                    df2 = df_ptr[cnt * ndf + 2];
                    df3 = df_ptr[cnt * ndf + 3];

                    d00 = varD00[cnt];
                    d01 = varD01[cnt];
                    d11 = varD11[cnt];

                    dtmp0 = df0 * d00;
                    dtmp0.fma(df2, d01);
                    dtmp1 = df0 * d01;
                    dtmp1.fma(df2, d11);
                    dtmp2 = df1 * d00;
                    dtmp2.fma(df3, d01);
                    dtmp3 = df1 * d01;
                    dtmp3.fma(df3, d11);

                    metric00 = df0 * dtmp0;
                    metric00.fma(df2, dtmp1);
                    metric01 = df1 * dtmp0;
                    metric01.fma(df3, dtmp1);
                    metric11 = df1 * dtmp2;
                    metric11.fma(df3, dtmp3);

                    vec_t d0 = deriv0[cnt];
                    vec_t d1 = deriv1[cnt];

                    vec_t tmp = metric00 * d0;
                    tmp.fma(metric01, d1);
                    bwd[cnt] = tmp;

                    tmp = metric01 * d0;
                    tmp.fma(metric11, d1);
                    deriv0[cnt] = tmp;
                }
            }
        }
        else
        {
            for (size_t j = 0, cnt = 0; j < nq1; ++j)
            {
                for (size_t i = 0; i < nq0; ++i, ++cnt)
                {
                    d00 = varD00[cnt];
                    d01 = varD01[cnt];
                    d11 = varD11[cnt];

                    dtmp0 = df0 * d00;
                    dtmp0.fma(df2, d01);
                    dtmp1 = df0 * d01;
                    dtmp1.fma(df2, d11);
                    dtmp2 = df1 * d00;
                    dtmp2.fma(df3, d01);
                    dtmp3 = df1 * d01;
                    dtmp3.fma(df3, d11);

                    metric00 = df0 * dtmp0;
                    metric00.fma(df2, dtmp1);
                    metric01 = df1 * dtmp0;
                    metric01.fma(df3, dtmp1);
                    metric11 = df1 * dtmp2;
                    metric11.fma(df3, dtmp3);

                    vec_t d0 = deriv0[cnt];
                    vec_t d1 = deriv1[cnt];

                    vec_t tmp = metric00 * d0;
                    tmp.fma(metric01, d1);
                    bwd[cnt] = tmp;

                    tmp = metric01 * d0;
                    tmp.fma(metric11, d1);
                    deriv0[cnt] = tmp;
                }
            }
        }
    }
}

#elif defined(SHAPE_TYPE_HEX)

template <bool DEFORMED>
NEK_FORCE_INLINE static void HelmholtzHexKernel(
    const size_t nq0, const size_t nq1, const size_t nq2,
    const bool isConstVarDiff, const Array<OneD, NekDouble> &constVarDiff,
    const bool isVarDiff, const Array<OneD, NekDouble> &varD00,
    const Array<OneD, NekDouble> &varD01, const Array<OneD, NekDouble> &varD11,
    const Array<OneD, NekDouble> &varD02, const Array<OneD, NekDouble> &varD12,
    const Array<OneD, NekDouble> &varD22, const vec_t *df_ptr,
    std::vector<vec_t, allocator<vec_t>> &deriv0,
    std::vector<vec_t, allocator<vec_t>> &deriv1,
    std::vector<vec_t, allocator<vec_t>> &deriv2)
{
    constexpr auto ndf = 9;
    const auto nqTot   = nq0 * nq1 * nq2;

    vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;
    vec_t metric00, metric01, metric02, metric11, metric12, metric22;

    vec_t d00 = {1.0};
    vec_t d01 = {0.0};
    vec_t d11 = {1.0};
    vec_t d02 = {0.0};
    vec_t d12 = {0.0};
    vec_t d22 = {1.0};                                 // var diffusion terms
    vec_t td0, td1, td2, td3, td4, td5, td6, td7, td8; // temp terms for vardiff

    if (isConstVarDiff)
    {
        d00 = constVarDiff[0];
        d01 = constVarDiff[1];
        d11 = constVarDiff[2];
        d02 = constVarDiff[3];
        d12 = constVarDiff[4];
        d22 = constVarDiff[5];
    }

    // Precompute Laplacian metricsp
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

        if (!isConstVarDiff && !isConstVarDiff)
        {
            metric00 = df0 * df0;
            metric00.fma(df3, df3);
            metric00.fma(df6, df6);

            metric01 = df0 * df1;
            metric01.fma(df3, df4);
            metric01.fma(df6, df7);

            metric02 = df0 * df2;
            metric02.fma(df3, df5);
            metric02.fma(df6, df8);

            metric11 = df1 * df1;
            metric11.fma(df4, df4);
            metric11.fma(df7, df7);

            metric12 = df1 * df2;
            metric12.fma(df4, df5);
            metric12.fma(df7, df8);

            metric22 = df2 * df2;
            metric22.fma(df5, df5);
            metric22.fma(df8, df8);
        }
        else if (isConstVarDiff)
        { // with vardiff
            td0 = df0 * d00;
            td0.fma(df3, d01);
            td0.fma(df6, d02);

            td1 = df0 * d01;
            td1.fma(df3, d11);
            td1.fma(df6, d12);

            td2 = df0 * d02;
            td2.fma(df3, d12);
            td2.fma(df6, d22);

            td3 = df1 * d00;
            td3.fma(df4, d01);
            td3.fma(df7, d02);

            td4 = df1 * d01;
            td4.fma(df4, d11);
            td4.fma(df7, d12);

            td5 = df1 * d02;
            td5.fma(df4, d12);
            td5.fma(df7, d22);

            td6 = df2 * d00;
            td6.fma(df5, d01);
            td6.fma(df8, d02);

            td7 = df2 * d01;
            td7.fma(df5, d11);
            td7.fma(df8, d12);

            td8 = df2 * d02;
            td8.fma(df5, d12);
            td8.fma(df8, d22);

            metric00 = td0 * df0;
            metric00.fma(td1, df3);
            metric00.fma(td2, df6);

            metric01 = td0 * df1;
            metric01.fma(td1, df4);
            metric01.fma(td2, df7);

            metric02 = td0 * df2;
            metric02.fma(td1, df5);
            metric02.fma(td2, df8);

            metric11 = td3 * df1;
            metric11.fma(td4, df4);
            metric11.fma(td5, df7);

            metric12 = td3 * df2;
            metric12.fma(td4, df5);
            metric12.fma(td5, df8);

            metric22 = td6 * df2;
            metric22.fma(td7, df5);
            metric22.fma(td8, df8);
        }
    }

    // Step 4: Apply Laplacian metrics & inner product
    if (DEFORMED || isVarDiff)
    {
        for (size_t k = 0, cnt = 0; k < nq2; k++)
        {
            for (size_t j = 0; j < nq1; j++)
            {
                for (size_t i = 0; i < nq0; i++, ++cnt)
                {
                    if (DEFORMED)
                    {
                        df0 = df_ptr[cnt * ndf];
                        df1 = df_ptr[cnt * ndf + 1];
                        df2 = df_ptr[cnt * ndf + 2];
                        df3 = df_ptr[cnt * ndf + 3];
                        df4 = df_ptr[cnt * ndf + 4];
                        df5 = df_ptr[cnt * ndf + 5];
                        df6 = df_ptr[cnt * ndf + 6];
                        df7 = df_ptr[cnt * ndf + 7];
                        df8 = df_ptr[cnt * ndf + 8];
                    }

                    if (!isConstVarDiff && !isVarDiff)
                    {
                        metric00 = df0 * df0;
                        metric00.fma(df3, df3);
                        metric00.fma(df6, df6);

                        metric01 = df0 * df1;
                        metric01.fma(df3, df4);
                        metric01.fma(df6, df7);

                        metric02 = df0 * df2;
                        metric02.fma(df3, df5);
                        metric02.fma(df6, df8);

                        metric11 = df1 * df1;
                        metric11.fma(df4, df4);
                        metric11.fma(df7, df7);

                        metric12 = df1 * df2;
                        metric12.fma(df4, df5);
                        metric12.fma(df7, df8);

                        metric22 = df2 * df2;
                        metric22.fma(df5, df5);
                        metric22.fma(df8, df8);
                    }
                    else
                    { // with vardiff

                        if (isVarDiff)
                        {
                            d00 = varD00[cnt];
                            d01 = varD01[cnt];
                            d11 = varD11[cnt];
                            d02 = varD02[cnt];
                            d12 = varD12[cnt];
                            d22 = varD22[cnt];
                        }

                        td0 = df0 * d00;
                        td0.fma(df3, d01);
                        td0.fma(df6, d02);

                        td1 = df0 * d01;
                        td1.fma(df3, d11);
                        td1.fma(df6, d12);

                        td2 = df0 * d02;
                        td2.fma(df3, d12);
                        td2.fma(df6, d22);

                        td3 = df1 * d00;
                        td3.fma(df4, d01);
                        td3.fma(df7, d02);

                        td4 = df1 * d01;
                        td4.fma(df4, d11);
                        td4.fma(df7, d12);

                        td5 = df1 * d02;
                        td5.fma(df4, d12);
                        td5.fma(df7, d22);

                        td6 = df2 * d00;
                        td6.fma(df5, d01);
                        td6.fma(df8, d02);

                        td7 = df2 * d01;
                        td7.fma(df5, d11);
                        td7.fma(df8, d12);

                        td8 = df2 * d02;
                        td8.fma(df5, d12);
                        td8.fma(df8, d22);

                        metric00 = td0 * df0;
                        metric00.fma(td1, df3);
                        metric00.fma(td2, df6);

                        metric01 = td0 * df1;
                        metric01.fma(td1, df4);
                        metric01.fma(td2, df7);

                        metric02 = td0 * df2;
                        metric02.fma(td1, df5);
                        metric02.fma(td2, df8);

                        metric11 = td3 * df1;
                        metric11.fma(td4, df4);
                        metric11.fma(td5, df7);

                        metric12 = td3 * df2;
                        metric12.fma(td4, df5);
                        metric12.fma(td5, df8);

                        metric22 = td6 * df2;
                        metric22.fma(td7, df5);
                        metric22.fma(td8, df8);
                    }

                    vec_t d0 = deriv0[cnt];
                    vec_t d1 = deriv1[cnt];
                    vec_t d2 = deriv2[cnt];

                    vec_t tmp = metric00 * d0;
                    tmp.fma(metric01, d1);
                    tmp.fma(metric02, d2);
                    deriv0[cnt] = tmp;

                    tmp = metric01 * d0;
                    tmp.fma(metric11, d1);
                    tmp.fma(metric12, d2);
                    deriv1[cnt] = tmp;

                    tmp = metric02 * d0;
                    tmp.fma(metric12, d1);
                    tmp.fma(metric22, d2);
                    deriv2[cnt] = tmp;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < nqTot; ++i)
        {
            vec_t d0 = deriv0[i];
            vec_t d1 = deriv1[i];
            vec_t d2 = deriv2[i];

            vec_t tmp = metric00 * d0;
            tmp.fma(metric01, d1);
            tmp.fma(metric02, d2);
            deriv0[i] = tmp;

            tmp = metric01 * d0;
            tmp.fma(metric11, d1);
            tmp.fma(metric12, d2);
            deriv1[i] = tmp;

            tmp = metric02 * d0;
            tmp.fma(metric12, d1);
            tmp.fma(metric22, d2);
            deriv2[i] = tmp;
        }
    }
}

#elif defined(SHAPE_TYPE_TET)

template <bool DEFORMED>
NEK_FORCE_INLINE static void HelmholtzTetKernel(
    const size_t nq0, const size_t nq1, const size_t nq2,
    const bool isConstVarDiff, const Array<OneD, NekDouble> &constVarDiff,
    const bool isVarDiff, const Array<OneD, NekDouble> &varD00,
    const Array<OneD, NekDouble> &varD01, const Array<OneD, NekDouble> &varD11,
    const Array<OneD, NekDouble> &varD02, const Array<OneD, NekDouble> &varD12,
    const Array<OneD, NekDouble> &varD22, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &m_h0,
    const std::vector<vec_t, allocator<vec_t>> &m_h1,
    const std::vector<vec_t, allocator<vec_t>> &m_h2,
    const std::vector<vec_t, allocator<vec_t>> &m_h3,
    std::vector<vec_t, allocator<vec_t>> &deriv0,
    std::vector<vec_t, allocator<vec_t>> &deriv1,
    std::vector<vec_t, allocator<vec_t>> &deriv2)
{
    boost::ignore_unused(varD11);

    constexpr auto ndf = 9;

    vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;
    vec_t g0, g1, g2, g3, g4, g5;
    vec_t d00 = {1.0};
    vec_t d01 = {0.0};
    vec_t d11 = {1.0};
    vec_t d02 = {0.0};
    vec_t d12 = {0.0};
    vec_t d22 = {1.0};                                 // var diffusion terms
    vec_t td0, td1, td2, td3, td4, td5, td6, td7, td8; // temp terms for vardiff

    if (isConstVarDiff)
    {
        d00 = constVarDiff[0];
        d01 = constVarDiff[1];
        d11 = constVarDiff[2];
        d02 = constVarDiff[3];
        d12 = constVarDiff[4];
        d22 = constVarDiff[5];
    }

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

    // Step 4a: Construct Laplacian metrics
    for (size_t k = 0, cnt = 0; k < nq2; ++k)
    {
        vec_t h3 = m_h3[k];
        for (size_t j = 0; j < nq1; ++j)
        {
            vec_t h1   = m_h1[j];
            vec_t h2   = m_h2[j];
            vec_t h2h3 = h2 * h3;
            vec_t h1h3 = h1 * h3;

            for (int i = 0; i < nq0; ++i, ++cnt)
            {
                if (DEFORMED)
                {
                    df0 = df_ptr[cnt * ndf];
                    df1 = df_ptr[cnt * ndf + 1];
                    df2 = df_ptr[cnt * ndf + 2];
                    df3 = df_ptr[cnt * ndf + 3];
                    df4 = df_ptr[cnt * ndf + 4];
                    df5 = df_ptr[cnt * ndf + 5];
                    df6 = df_ptr[cnt * ndf + 6];
                    df7 = df_ptr[cnt * ndf + 7];
                    df8 = df_ptr[cnt * ndf + 8];
                }

                vec_t h0h2h3 = m_h0[i] * h2h3;

                vec_t tmp1 = h0h2h3 * (df1 + df2);
                tmp1.fma(df0, h2h3);
                vec_t tmp2 = h0h2h3 * (df4 + df5);
                tmp2.fma(df3, h2h3);
                vec_t tmp3 = h0h2h3 * (df7 + df8);
                tmp3.fma(df6, h2h3);

                vec_t tmp4 = df1 * h3;
                tmp4.fma(df2, h1h3);
                vec_t tmp5 = df4 * h3;
                tmp5.fma(df5, h1h3);
                vec_t tmp6 = df7 * h3;
                tmp6.fma(df8, h1h3);

                if (!isConstVarDiff && !isVarDiff)
                {
                    g0 = tmp1 * tmp1;
                    g0.fma(tmp2, tmp2);
                    g0.fma(tmp3, tmp3);

                    g4 = df2 * tmp1;
                    g4.fma(df5, tmp2);
                    g4.fma(df8, tmp3);

                    g3 = tmp1 * tmp4;
                    g3.fma(tmp2, tmp5);
                    g3.fma(tmp3, tmp6);

                    g1 = tmp4 * tmp4;
                    g1.fma(tmp5, tmp5);
                    g1.fma(tmp6, tmp6);

                    g5 = df2 * tmp4;
                    g5.fma(df5, tmp5);
                    g5.fma(df8, tmp6);

                    g2 = df2 * df2;
                    g2.fma(df5, df5);
                    g2.fma(df8, df8);
                }
                else
                {
                    if (isVarDiff)
                    {
                        d00 = varD00[cnt];
                        d01 = varD01[cnt];
                        d02 = varD02[cnt];
                        d12 = varD12[cnt];
                        d22 = varD22[cnt];
                    }

                    td0 = tmp1 * d00;
                    td0.fma(tmp2, d01);
                    td0.fma(tmp3, d02);

                    td1 = tmp1 * d01;
                    td1.fma(tmp2, d11);
                    td1.fma(tmp3, d12);

                    td2 = tmp1 * d02;
                    td2.fma(tmp2, d12);
                    td2.fma(tmp3, d22);

                    td3 = tmp4 * d00;
                    td3.fma(tmp5, d01);
                    td3.fma(tmp6, d02);

                    td4 = tmp4 * d01;
                    td4.fma(tmp5, d11);
                    td4.fma(tmp6, d12);

                    td5 = tmp4 * d02;
                    td5.fma(tmp5, d12);
                    td5.fma(tmp6, d22);

                    td6 = df2 * d00;
                    td6.fma(df5, d01);
                    td6.fma(df8, d02);

                    td7 = df2 * d01;
                    td7.fma(df5, d11);
                    td7.fma(df8, d12);

                    td8 = df2 * d02;
                    td8.fma(df5, d12);
                    td8.fma(df8, d22);

                    g0 = td0 * tmp1;
                    g0.fma(td1, tmp2);
                    g0.fma(td2, tmp3);

                    g3 = td0 * tmp4;
                    g3.fma(td1, tmp5);
                    g3.fma(td2, tmp6);

                    g4 = td0 * df2;
                    g4.fma(td1, df5);
                    g4.fma(td2, df8);

                    g1 = td3 * tmp4;
                    g1.fma(td4, tmp5);
                    g1.fma(td5, tmp6);

                    g5 = td3 * df2;
                    g5.fma(td4, df5);
                    g5.fma(td5, df8);

                    g2 = td6 * df2;
                    g2.fma(td7, df5);
                    g2.fma(td8, df8);
                }

                vec_t d0 = deriv0[cnt];
                vec_t d1 = deriv1[cnt];
                vec_t d2 = deriv2[cnt];

                tmp1 = g0 * d0;
                tmp1.fma(g3, d1);
                tmp1.fma(g4, d2);
                deriv0[cnt] = tmp1;

                tmp2 = g3 * d0;
                tmp2.fma(g1, d1);
                tmp2.fma(g5, d2);
                deriv1[cnt] = tmp2;

                tmp3 = g4 * d0;
                tmp3.fma(g5, d1);
                tmp3.fma(g2, d2);
                deriv2[cnt] = tmp3;
            }
        }
    }
}

#elif defined(SHAPE_TYPE_PRISM)

template <bool DEFORMED>
NEK_FORCE_INLINE static void HelmholtzPrismKernel(
    const size_t nq0, const size_t nq1, const size_t nq2,
    const bool isConstVarDiff, const Array<OneD, NekDouble> &constVarDiff,
    const bool isVarDiff, const Array<OneD, NekDouble> &varD00,
    const Array<OneD, NekDouble> &varD01, const Array<OneD, NekDouble> &varD11,
    const Array<OneD, NekDouble> &varD02, const Array<OneD, NekDouble> &varD12,
    const Array<OneD, NekDouble> &varD22, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &m_h0,
    const std::vector<vec_t, allocator<vec_t>> &m_h1,
    std::vector<vec_t, allocator<vec_t>> &deriv0,
    std::vector<vec_t, allocator<vec_t>> &deriv1,
    std::vector<vec_t, allocator<vec_t>> &deriv2)
{
    constexpr auto ndf = 9;

    vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;
    vec_t g0, g1, g2, g3, g4, g5; // metrics
    vec_t d00 = {1.0};
    vec_t d01 = {0.0};
    vec_t d11 = {1.0};
    vec_t d02 = {0.0};
    vec_t d12 = {0.0};
    vec_t d22 = {1.0};                                 // var diffusion terms
    vec_t td0, td1, td2, td3, td4, td5, td6, td7, td8; // temp terms for vardiff
    boost::ignore_unused(d00, d01, d11, d02, d12, d22, td0, td1, td2, td3, td4,
                         td5, td6, td7, td8);

    if (isConstVarDiff)
    {
        d00 = constVarDiff[0];
        d01 = constVarDiff[1];
        d11 = constVarDiff[2];
        d02 = constVarDiff[3];
        d12 = constVarDiff[4];
        d22 = constVarDiff[5];
    }

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

    // Step 4a: Construct Laplacian metrics
    for (size_t k = 0, cnt = 0; k < nq2; ++k)
    {
        vec_t h1 = m_h1[k];
        for (size_t j = 0; j < nq1; ++j)
        {
            for (size_t i = 0; i < nq0; ++i, cnt++)
            {
                vec_t h0 = m_h0[i];

                if (DEFORMED)
                {
                    df0 = df_ptr[cnt * ndf];
                    df1 = df_ptr[cnt * ndf + 1];
                    df2 = df_ptr[cnt * ndf + 2];
                    df3 = df_ptr[cnt * ndf + 3];
                    df4 = df_ptr[cnt * ndf + 4];
                    df5 = df_ptr[cnt * ndf + 5];
                    df6 = df_ptr[cnt * ndf + 6];
                    df7 = df_ptr[cnt * ndf + 7];
                    df8 = df_ptr[cnt * ndf + 8];
                }

                vec_t tmp1 = h1 * (h0 * df2 + df0);
                vec_t tmp2 = h1 * (h0 * df5 + df3);
                vec_t tmp3 = h1 * (h0 * df8 + df6);

                if (!isConstVarDiff && !isVarDiff)
                {
                    g0 = tmp1 * tmp1;
                    g0.fma(tmp2, tmp2);
                    g0.fma(tmp3, tmp3);

                    g3 = df1 * tmp1;
                    g3.fma(df4, tmp2);
                    g3.fma(df7, tmp3);

                    g4 = df2 * tmp1;
                    g4.fma(df5, tmp2);
                    g4.fma(df8, tmp3);

                    g1 = df1 * df1;
                    g1.fma(df4, df4);
                    g1.fma(df7, df7);

                    g2 = df2 * df2;
                    g2.fma(df5, df5);
                    g2.fma(df8, df8);

                    g5 = df1 * df2;
                    g5.fma(df4, df5);
                    g5.fma(df7, df8);
                }
                else
                {
                    // vardiff
                    if (isVarDiff)
                    {
                        d00 = varD00[cnt];
                        d01 = varD01[cnt];
                        d11 = varD11[cnt];
                        d02 = varD02[cnt];
                        d12 = varD12[cnt];
                        d22 = varD22[cnt];
                    }

                    td0 = tmp1 * d00;
                    td0.fma(tmp2, d01);
                    td0.fma(tmp3, d02);

                    td1 = tmp1 * d01;
                    td1.fma(tmp2, d11);
                    td1.fma(tmp3, d12);

                    td2 = tmp1 * d02;
                    td2.fma(tmp2, d12);
                    td2.fma(tmp3, d22);

                    td3 = df1 * d00;
                    td3.fma(df4, d01);
                    td3.fma(df7, d02);

                    td4 = df1 * d01;
                    td4.fma(df4, d11);
                    td4.fma(df7, d12);

                    td5 = df1 * d02;
                    td5.fma(df4, d12);
                    td5.fma(df7, d22);

                    td6 = df2 * d00;
                    td6.fma(df5, d01);
                    td6.fma(df8, d02);

                    td7 = df2 * d01;
                    td7.fma(df5, d11);
                    td7.fma(df8, d12);

                    td8 = df2 * d02;
                    td8.fma(df5, d12);
                    td8.fma(df8, d22);

                    g0 = td0 * tmp1;
                    g0.fma(td1, tmp2);
                    g0.fma(td2, tmp3);

                    g3 = td0 * df1;
                    g3.fma(td1, df4);
                    g3.fma(td2, df7);

                    g4 = td0 * df2;
                    g4.fma(td1, df5);
                    g4.fma(td2, df8);

                    g1 = td3 * df1;
                    g1.fma(td4, df4);
                    g1.fma(td5, df7);

                    g5 = td3 * df2;
                    g5.fma(td4, df5);
                    g5.fma(td5, df8);

                    g2 = td6 * df2;
                    g2.fma(td7, df5);
                    g2.fma(td8, df8);
                }

                vec_t d0 = deriv0[cnt];
                vec_t d1 = deriv1[cnt];
                vec_t d2 = deriv2[cnt];

                tmp1 = g0 * d0;
                tmp1.fma(g3, d1);
                tmp1.fma(g4, d2);
                deriv0[cnt] = tmp1;

                tmp2 = g3 * d0;
                tmp2.fma(g1, d1);
                tmp2.fma(g5, d2);
                deriv1[cnt] = tmp2;

                tmp3 = g4 * d0;
                tmp3.fma(g5, d1);
                tmp3.fma(g2, d2);
                deriv2[cnt] = tmp3;
            }
        }
    }
}

#elif defined(SHAPE_TYPE_PYR)

template <bool DEFORMED>
NEK_FORCE_INLINE static void HelmholtzPyrKernel(
    const size_t nq0, const size_t nq1, const size_t nq2,
    const bool isConstVarDiff, const Array<OneD, NekDouble> &constVarDiff,
    const bool isVarDiff, const Array<OneD, NekDouble> &varD00,
    const Array<OneD, NekDouble> &varD01, const Array<OneD, NekDouble> &varD11,
    const Array<OneD, NekDouble> &varD02, const Array<OneD, NekDouble> &varD12,
    const Array<OneD, NekDouble> &varD22, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &m_h0,
    const std::vector<vec_t, allocator<vec_t>> &m_h1,
    const std::vector<vec_t, allocator<vec_t>> &m_h2,
    std::vector<vec_t, allocator<vec_t>> &deriv0,
    std::vector<vec_t, allocator<vec_t>> &deriv1,
    std::vector<vec_t, allocator<vec_t>> &deriv2)
{
    constexpr auto ndf = 9;

    vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;
    vec_t g0, g1, g2, g3, g4, g5; // metrics
    vec_t d00 = {1.0};
    vec_t d01 = {0.0};
    vec_t d11 = {1.0};
    vec_t d02 = {0.0};
    vec_t d12 = {0.0};
    vec_t d22 = {1.0};                                 // var diffusion terms
    vec_t td0, td1, td2, td3, td4, td5, td6, td7, td8; // temp terms for vardiff

    if (isConstVarDiff)
    {
        d00 = constVarDiff[0];
        d01 = constVarDiff[1];
        d11 = constVarDiff[2];
        d02 = constVarDiff[3];
        d12 = constVarDiff[4];
        d22 = constVarDiff[5];
    }

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

    // Step 4a: Construct Laplacian metrics
    for (size_t k = 0, cnt = 0; k < nq2; ++k)
    {
        vec_t h2 = m_h2[k];
        for (size_t j = 0; j < nq1; ++j)
        {
            vec_t h1   = m_h1[j];
            vec_t h1h2 = h1 * h2;
            for (size_t i = 0; i < nq0; ++i, cnt++)
            {
                vec_t h0   = m_h0[i];
                vec_t h0h2 = h0 * h2;

                if (DEFORMED)
                {
                    df0 = df_ptr[cnt * ndf];
                    df1 = df_ptr[cnt * ndf + 1];
                    df2 = df_ptr[cnt * ndf + 2];
                    df3 = df_ptr[cnt * ndf + 3];
                    df4 = df_ptr[cnt * ndf + 4];
                    df5 = df_ptr[cnt * ndf + 5];
                    df6 = df_ptr[cnt * ndf + 6];
                    df7 = df_ptr[cnt * ndf + 7];
                    df8 = df_ptr[cnt * ndf + 8];
                }

                vec_t tmp0 = h2 * df0;
                tmp0.fma(h0h2, df2);
                vec_t tmp1 = h2 * df3;
                tmp1.fma(h0h2, df5);
                vec_t tmp2 = h2 * df6;
                tmp2.fma(h0h2, df8);

                vec_t tmp3 = h2 * df1;
                tmp3.fma(h1h2, df2);
                vec_t tmp4 = h2 * df4;
                tmp4.fma(h1h2, df5);
                vec_t tmp5 = h2 * df7;
                tmp5.fma(h1h2, df8);

                if (!isConstVarDiff && !isVarDiff)
                {
                    g0 = tmp0 * tmp0;
                    g0.fma(tmp1, tmp1);
                    g0.fma(tmp2, tmp2);

                    g1 = tmp3 * tmp3;
                    g1.fma(tmp4, tmp4);
                    g1.fma(tmp5, tmp5);

                    g2 = df2 * df2;
                    g2.fma(df5, df5);
                    g2.fma(df8, df8);

                    g3 = tmp0 * tmp3;
                    g3.fma(tmp1, tmp4);
                    g3.fma(tmp2, tmp5);

                    g4 = df2 * tmp0;
                    g4.fma(df5, tmp1);
                    g4.fma(df8, tmp2);

                    g5 = df2 * tmp3;
                    g5.fma(df5, tmp4);
                    g5.fma(df8, tmp5);
                }
                else
                {
                    if (isVarDiff)
                    {
                        d00 = varD00[cnt];
                        d01 = varD01[cnt];
                        d11 = varD11[cnt];
                        d02 = varD02[cnt];
                        d12 = varD12[cnt];
                        d22 = varD22[cnt];
                    }

                    td0 = tmp0 * d00;
                    td0.fma(tmp1, d01);
                    td0.fma(tmp2, d02);

                    td1 = tmp0 * d01;
                    td1.fma(tmp1, d11);
                    td1.fma(tmp2, d12);

                    td2 = tmp0 * d02;
                    td2.fma(tmp1, d12);
                    td2.fma(tmp2, d22);

                    td3 = tmp3 * d00;
                    td3.fma(tmp4, d01);
                    td3.fma(tmp5, d02);

                    td4 = tmp3 * d01;
                    td4.fma(tmp4, d11);
                    td4.fma(tmp5, d12);

                    td5 = tmp3 * d02;
                    td5.fma(tmp4, d12);
                    td5.fma(tmp5, d22);

                    td6 = df2 * d00;
                    td6.fma(df5, d01);
                    td6.fma(df8, d02);

                    td7 = df2 * d01;
                    td7.fma(df5, d11);
                    td7.fma(df8, d12);

                    td8 = df2 * d02;
                    td8.fma(df5, d12);
                    td8.fma(df8, d22);

                    g0 = td0 * tmp0;
                    g0.fma(td1, tmp1);
                    g0.fma(td2, tmp2);

                    g3 = td0 * tmp3;
                    g3.fma(td1, tmp4);
                    g3.fma(td2, tmp5);

                    g4 = td0 * df2;
                    g4.fma(td1, df5);
                    g4.fma(td2, df8);

                    g1 = td3 * tmp3;
                    g1.fma(td4, tmp4);
                    g1.fma(td5, tmp5);

                    g5 = td3 * df2;
                    g5.fma(td4, df5);
                    g5.fma(td5, df8);

                    g2 = td6 * df2;
                    g2.fma(td7, df5);
                    g2.fma(td8, df8);
                }

                vec_t d0 = deriv0[cnt];
                vec_t d1 = deriv1[cnt];
                vec_t d2 = deriv2[cnt];

                tmp1 = g0 * d0;
                tmp1.fma(g3, d1);
                tmp1.fma(g4, d2);
                deriv0[cnt] = tmp1;

                tmp2 = g3 * d0;
                tmp2.fma(g1, d1);
                tmp2.fma(g5, d2);
                deriv1[cnt] = tmp2;

                tmp3 = g4 * d0;
                tmp3.fma(g5, d1);
                tmp3.fma(g2, d2);
                deriv2[cnt] = tmp3;
            }
        }
    }
}

#endif // SHAPE_TYPE

#if defined(SHAPE_TYPE_TRI)
template <LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void GetHelmholtz2DHalfSpace(
    const size_t nq0, const size_t nq1, const Array<OneD, const NekDouble> &z0,
    const Array<OneD, const NekDouble> &z1,
    std::vector<vec_t, allocator<vec_t>> &h0,
    std::vector<vec_t, allocator<vec_t>> &h1)
{
    boost::ignore_unused(SHAPE_TYPE);

    h0.resize(nq0);
    h1.resize(nq1);

    for (int i = 0; i < nq0; ++i)
    {
        h0[i] = 0.5 * (1 + z0[i]);
    }

    for (int j = 0; j < nq1; ++j)
    {
        h1[j] = 2.0 / (1 - z1[j]);
    }
}

#elif defined(SHAPE_TYPE_TET) || defined(SHAPE_TYPE_PRISM) ||                  \
    defined(SHAPE_TYPE_PYR)
template <LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void GetHelmholtz3DHalfSpace(
    const size_t nq0, const size_t nq1, const size_t nq2,
    const Array<OneD, const NekDouble> &z0,
    const Array<OneD, const NekDouble> &z1,
    const Array<OneD, const NekDouble> &z2,
    std::vector<vec_t, allocator<vec_t>> &h0,
    std::vector<vec_t, allocator<vec_t>> &h1,
    std::vector<vec_t, allocator<vec_t>> &h2,
    std::vector<vec_t, allocator<vec_t>> &h3)
{
    boost::ignore_unused(SHAPE_TYPE, nq1, z1, h2, h3);

#if defined(SHAPE_TYPE_TET)
    h0.resize(nq0);
    h1.resize(nq1);
    h2.resize(nq1);
    h3.resize(nq2);

    for (int i = 0; i < nq0; ++i)
    {
        h0[i] = 0.5 * (1 + z0[i]);
    }

    for (int j = 0; j < nq1; ++j)
    {
        h1[j] = 0.5 * (1 + z1[j]);
        h2[j] = 2.0 / (1 - z1[j]);
    }

    for (int k = 0; k < nq2; ++k)
    {
        h3[k] = 2.0 / (1 - z2[k]);
    }

#elif defined(SHAPE_TYPE_PRISM)

    h0.resize(nq0);
    h1.resize(nq2);

    for (int i = 0; i < nq0; ++i)
    {
        h0[i] = 0.5 * (1 + z0[i]);
    }

    for (int k = 0; k < nq2; ++k)
    {
        h1[k] = 2.0 / (1 - z2[k]);
    }

#elif defined(SHAPE_TYPE_PYR)

    h0.resize(nq0);
    h1.resize(nq1);
    h2.resize(nq2);

    for (int i = 0; i < nq0; ++i)
    {
        h0[i] = 0.5 * (1 + z0[i]);
    }

    for (int j = 0; j < nq1; ++j)
    {
        h1[j] = 0.5 * (1 + z1[j]);
    }

    for (int k = 0; k < nq2; ++k)
    {
        h2[k] = 2.0 / (1 - z2[k]);
    }
#endif
}
#endif

// The dimension kernels which select the shape kernel.
#if defined(SHAPE_DIMENSION_2D)

template <LibUtilities::ShapeType SHAPE_TYPE, bool DEFORMED>
NEK_FORCE_INLINE static void Helmholtz2DKernel(
    const size_t nq0, const size_t nq1, const bool isConstVarDiff,
    const Array<OneD, NekDouble> &constVarDiff, const bool isVarDiff,
    const Array<OneD, NekDouble> &varD00, const Array<OneD, NekDouble> &varD01,
    const Array<OneD, NekDouble> &varD11, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &h0,
    const std::vector<vec_t, allocator<vec_t>> &h1,
    std::vector<vec_t, allocator<vec_t>> &bwd,
    std::vector<vec_t, allocator<vec_t>> &deriv0,
    std::vector<vec_t, allocator<vec_t>> &deriv1)
{
#if defined(SHAPE_TYPE_TRI)
    HelmholtzTriKernel<DEFORMED>(nq0, nq1, isConstVarDiff, constVarDiff,
                                 isVarDiff, varD00, varD01, varD11, df_ptr, h0,
                                 h1, bwd, deriv0, deriv1);
#elif defined(SHAPE_TYPE_QUAD)
    boost::ignore_unused(h0, h1);
    HelmholtzQuadKernel<DEFORMED>(nq0, nq1, isConstVarDiff, constVarDiff,
                                  isVarDiff, varD00, varD01, varD11, df_ptr,
                                  bwd, deriv0, deriv1);
#endif
}

#elif defined(SHAPE_DIMENSION_3D)

template <LibUtilities::ShapeType SHAPE_TYPE, bool DEFORMED>
NEK_FORCE_INLINE static void Helmholtz3DKernel(
    const size_t nq0, const size_t nq1, const size_t nq2,
    const bool isConstVarDiff, const Array<OneD, NekDouble> &constVarDiff,
    const bool isVarDiff, const Array<OneD, NekDouble> &varD00,
    const Array<OneD, NekDouble> &varD01, const Array<OneD, NekDouble> &varD11,
    const Array<OneD, NekDouble> &varD02, const Array<OneD, NekDouble> &varD12,
    const Array<OneD, NekDouble> &varD22, const vec_t *df_ptr,
    const std::vector<vec_t, allocator<vec_t>> &h0,
    const std::vector<vec_t, allocator<vec_t>> &h1,
    const std::vector<vec_t, allocator<vec_t>> &h2,
    const std::vector<vec_t, allocator<vec_t>> &h3,
    std::vector<vec_t, allocator<vec_t>> &deriv0,
    std::vector<vec_t, allocator<vec_t>> &deriv1,
    std::vector<vec_t, allocator<vec_t>> &deriv2)
{
#if defined(SHAPE_TYPE_HEX)
    boost::ignore_unused(h0, h1, h2, h3);
    HelmholtzHexKernel<DEFORMED>(
        nq0, nq1, nq2, isConstVarDiff, constVarDiff, isVarDiff, varD00, varD01,
        varD11, varD02, varD12, varD22, df_ptr, deriv0, deriv1, deriv2);
#elif defined(SHAPE_TYPE_TET)
    HelmholtzTetKernel<DEFORMED>(nq0, nq1, nq2, isConstVarDiff, constVarDiff,
                                 isVarDiff, varD00, varD01, varD11, varD02,
                                 varD12, varD22, df_ptr, h0, h1, h2, h3, deriv0,
                                 deriv1, deriv2);
#elif defined(SHAPE_TYPE_PRISM)
    boost::ignore_unused(h2, h3);
    HelmholtzPrismKernel<DEFORMED>(
        nq0, nq1, nq2, isConstVarDiff, constVarDiff, isVarDiff, varD00, varD01,
        varD11, varD02, varD12, varD22, df_ptr, h0, h1, deriv0, deriv1, deriv2);
#elif defined(SHAPE_TYPE_PYR)
    boost::ignore_unused(h3);
    HelmholtzPyrKernel<DEFORMED>(nq0, nq1, nq2, isConstVarDiff, constVarDiff,
                                 isVarDiff, varD00, varD01, varD11, varD02,
                                 varD12, varD22, df_ptr, h0, h1, h2, deriv0,
                                 deriv1, deriv2);
#endif
}

#endif // SHAPE_DIMENSION

} // namespace MatrixFree
} // namespace Nektar

#endif
