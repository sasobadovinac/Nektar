#ifndef NEKTAR_LIBRARY_MF_IPRODUCT_KERNELS_HPP
#define NEKTAR_LIBRARY_MF_IPRODUCT_KERNELS_HPP

namespace Nektar
{
namespace MatrixFree
{

using namespace tinysimd;
using vec_t = simd<NekDouble>;

template<bool SCALE, bool APPEND>
NEK_FORCE_INLINE static void ScaleAppend(vec_t &store,
                                         vec_t &pos,
                                         NekDouble scale)
{
    if (SCALE && APPEND)
    {
        store.fma(pos, vec_t(scale));
    }
    else if (APPEND)
    {
        store = store + pos;
    }
    else if (SCALE)
    {
        store = pos * scale;
    }
    else
    {
        store = pos;
    }
}

// The dimension and shape kernels. NOTE: They are NOT duplicate
// templated version based on the array size like the
// operators. HOWEVER, they are forced to be INLINED. The inlining is
// critical so that when used in the templated version of the operator
// that loop unrolling occurs.

// The seven shape kernels where the work gets done.
#if defined(SHAPE_TYPE_SEG)

template<bool SCALE, bool APPEND, bool DEFORMED>
NEK_FORCE_INLINE static void IProductSegKernel(
    const size_t nm0,
    const size_t nq0,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& w0,
    const vec_t* jac,
          std::vector<vec_t, allocator<vec_t>> &out,
          NekDouble scale = 1.0)
{
    for (int p = 0; p < nm0; ++p)
    {
        vec_t sum = 0.0;

        for (int i = 0; i < nq0; ++i)
        {
            vec_t jac_val;

            if(DEFORMED)
            {
                jac_val = jac[i]; //J for each quadrature point.
            }
            else
            {
                jac_val = jac[0];
            }

            vec_t prod = in[i] * basis0[p*nq0 + i] * jac_val; //Load 2x
            sum.fma(prod, w0[i]); //Load 1x
        }

        // Modes are reversed from what they normally are for tris, tets etc.
        ScaleAppend<SCALE, APPEND>(out[p], sum, scale); // Store x1
    }
}

#elif defined(SHAPE_TYPE_TRI)

template<bool SCALE, bool APPEND, bool DEFORMED>
NEK_FORCE_INLINE static void IProductTriKernel(
    const size_t nm0, const size_t nm1,
    const size_t nq0, const size_t nq1,
    const bool correct,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
    const std::vector<vec_t, allocator<vec_t>>& w0,
    const std::vector<vec_t, allocator<vec_t>>& w1,
    const vec_t* jac,
          std::vector<vec_t, allocator<vec_t>>& eta0_sums, // nq1
          std::vector<vec_t, allocator<vec_t>>& out,
          NekDouble scale = 1.0)
{
    int mode = 0;

    for (int p = 0; p < nm0; ++p)
    {
        int eta_idx = 0;

        // Our inner loop is phi_p not phi_pq since we want to put as
        // much work as we can in the p-only loop instead of the full
        // pq loop.
        for (int eta1 = 0; eta1 < nq1; ++eta1)
        {
            vec_t eta0_sum = 0.0;

            for (int eta0 = 0; eta0 < nq0; ++eta0, ++eta_idx)
            {
                //eta0_sum += phi_p(eta0) * fn(eta0, eta1) * J * w0(eta0)
                vec_t jac_val;
                if (DEFORMED)
                {
                    jac_val = jac[eta1*nq0 + eta0];
                }
                else
                {
                    jac_val = jac[0];
                }

                vec_t prod = in[eta_idx] * basis0[p*nq0 + eta0] * jac_val; //Load 2x
                eta0_sum.fma(prod, w0[eta0]); //Load 1x
            }

            eta0_sums[eta1] = eta0_sum;
        }

        for (int q = 0; q < nm1 - p; ++q, ++mode)
        {
            vec_t sum_eta1 = 0.0;
            for (int eta1 = 0; eta1 < nq1; ++eta1)
            {
                vec_t prod = eta0_sums[eta1] * basis1[mode*nq1 + eta1]; //Load 2x
                sum_eta1.fma(prod, w1[eta1]); //Load 1x
            }

            ScaleAppend<SCALE, APPEND>(out[mode], sum_eta1, scale); //Store x1
        }
    }

    // Correction for singular vertex in collpased coordinates.
    // Basically we add phi_1 * phi_01 * (weighting, etc) to mode 00
    // With contributions from every quadrature point
    if (correct)
    {
        int eta_idx = 0;
        vec_t iprod_01 = 0.0; //T(outptr + VW); //Load 1x

        for (int eta1 = 0; eta1 < nq1; ++eta1)
        {
            vec_t preweight_eta1;

            if (DEFORMED)
            {
                preweight_eta1 = w1[eta1] * basis1[nq1 + eta1];
            }
            else
            {
                preweight_eta1 = w1[eta1] * jac[0] * basis1[nq1 + eta1];
            }

            for (int eta0 = 0; eta0 < nq0; ++eta0, ++eta_idx)
            {
                vec_t prod = in[eta_idx] * preweight_eta1 * w0[eta0];

                if (DEFORMED)
                {
                    prod = prod * jac[eta1*nq0 + eta0];
                }

                vec_t basis_val1 = basis0[nq0 + eta0];
                iprod_01.fma(prod, basis_val1);
            }
        }

        ScaleAppend<SCALE, true>(out[1], iprod_01, scale);
    }
}

#elif defined(SHAPE_TYPE_QUAD)

template<bool SCALE, bool APPEND, bool DEFORMED>
NEK_FORCE_INLINE static void IProductQuadKernel(
    const size_t nm0, const size_t nm1,
    const size_t nq0, const size_t nq1,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
    const std::vector<vec_t, allocator<vec_t>>& w0,
    const std::vector<vec_t, allocator<vec_t>>& w1,
    const vec_t* jac,
          std::vector<vec_t, allocator<vec_t>>& sums_j, // nq1
          std::vector<vec_t, allocator<vec_t>>& out,
          NekDouble scale = 1.0)
{
    for (int p = 0; p < nm0; ++p)
    {
        int cnt_ji = 0;
        for (int j = 0; j < nq1; ++j)
        {
            vec_t sum_j = 0.0;

            for (int i = 0; i < nq0; ++i, ++cnt_ji)
            {
                vec_t jac_val;

                if(DEFORMED)
                {
                    jac_val = jac[j*nq0 + i]; //J for each quadrature point.
                }
                else
                {
                    jac_val = jac[0];
                }

                vec_t prod = in[cnt_ji] * basis0[p*nq0 + i] * jac_val; //Load 2x
                sum_j.fma(prod, w0[i]); //Load 1x
            }

            sums_j[j] = sum_j; //Store 1
        }

        for (int q = 0; q < nm1; ++q)
        {
            vec_t sum = 0.0;

            for (int j = 0; j < nq1; ++j)
            {
                vec_t prod = sums_j[j] * basis1[q*nq1 + j]; //Load 2x
                sum.fma(prod, w1[j]); // Load 1x
            }

            //Modes are reversed from what they normally are for tris, tets etc.
            ScaleAppend<SCALE, APPEND>(out[q*nm1 + p], sum, scale); // Store x1
        }
    }
}

#elif defined(SHAPE_TYPE_HEX)

template<bool SCALE, bool APPEND, bool DEFORMED>
NEK_FORCE_INLINE static void IProductHexKernel(
    const size_t nm0, const size_t nm1, const size_t nm2,
    const size_t nq0, const size_t nq1, const size_t nq2,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
    const std::vector<vec_t, allocator<vec_t>>& basis2,
    const std::vector<vec_t, allocator<vec_t>>& w0,
    const std::vector<vec_t, allocator<vec_t>>& w1,
    const std::vector<vec_t, allocator<vec_t>>& w2,
    const vec_t* jac,
          std::vector<vec_t, allocator<vec_t>>& sums_kj, // nq2 * nq1
          std::vector<vec_t, allocator<vec_t>>& sums_k,  // nq2
          std::vector<vec_t, allocator<vec_t>>& out,
          NekDouble scale = 1.0)
{
    for (int p = 0; p < nm0; ++p)
    {
        int cnt_kji = 0, cnt_kj = 0;

        for (int k = 0; k < nq2; ++k)
        {
            for (int j = 0; j < nq1; ++j, ++cnt_kj)
            {
                vec_t sum_kj = 0.0;

                for (int i = 0; i < nq0; ++i, ++cnt_kji)
                {

                    vec_t jac_val;

                    if(DEFORMED)
                    {
                        jac_val = jac[nq0*nq1*k + nq0*j + i];
                    }
                    else
                    {
                        jac_val = jac[0];
                    }

                    vec_t prod = in[cnt_kji] * basis0[i + nq0*p] * jac_val; // load 2x
                    sum_kj.fma(prod, w0[i]); //Load 1x
                }

                sums_kj[cnt_kj] = sum_kj;
            }
        }

        for (int q = 0; q < nm1; ++q)
        {
            cnt_kj = 0;
            for (int k = 0; k < nq2; ++k)
            {
                vec_t sum_k = 0.0;

                for (int j = 0; j < nq1; ++j, ++cnt_kj)
                {
                    vec_t prod = sums_kj[cnt_kj] * basis1[q*nq1 + j]; //Load 2x
                    sum_k.fma(prod, w1[j]); //Load 1x
                }

                sums_k[k] = sum_k;
            }

            for (int r = 0; r < nm2; ++r)
            {
                vec_t sum = 0.0;

                for(int k = 0; k < nq2; ++k)
                {
                    vec_t prod = sums_k[k] * basis2[r*nq2 + k]; //Load 2x
                    sum.fma(prod,w2[k]); //Load 1x
                }

                ScaleAppend<SCALE, APPEND>(out[r*nm0*nm1 + q*nm0 + p], sum, scale); // Store x1
            }
        }
    }
}

#elif defined(SHAPE_TYPE_TET)

template<bool SCALE, bool APPEND, bool DEFORMED>
NEK_FORCE_INLINE static void IProductTetKernel(
    const size_t nm0, const size_t nm1, const size_t nm2,
    const size_t nq0, const size_t nq1, const size_t nq2,
    const bool correct,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
    const std::vector<vec_t, allocator<vec_t>>& basis2,
    const std::vector<vec_t, allocator<vec_t>>& w0,
    const std::vector<vec_t, allocator<vec_t>>& w1,
    const std::vector<vec_t, allocator<vec_t>>& w2,
    const vec_t* jac,
          std::vector<vec_t, allocator<vec_t>>& sums_kj, // nq2 * nq1
          std::vector<vec_t, allocator<vec_t>>& sums_k,  // nq2
          std::vector<vec_t, allocator<vec_t>>& out,
          NekDouble scale = 1.0)
{
    for (int p = 0, mode = 0, mode2 = 0, cnt_pqr = 0; p < nm0; ++p)
    {
        int cnt_kji = 0;

        for (int k = 0, cnt_kj = 0; k < nq2; ++k)
        {
            for (int j = 0; j < nq1; ++j, ++cnt_kj)
            {
                // Unroll first entry of for loop below and multiply by
                // quadrature weights in dir0 & jacobian.
                vec_t jac_val;

                if (DEFORMED)
                {
                    jac_val = jac[nq0*nq1*k + nq0*j];
                }
                else
                {
                    jac_val = jac[0];
                }

                vec_t sum_kj = in[cnt_kji] * basis0[nq0*p] * jac_val * w0[0]; //Load 3x
                ++cnt_kji;

                for (int i = 1; i < nq0; ++i, ++cnt_kji)
                {
                    if (DEFORMED)
                    {
                        jac_val = jac[nq0*nq1*k + nq0*j + i];
                    }
                    else
                    {
                        jac_val = jac[0];
                    }

                    vec_t inxmm = in[cnt_kji] * basis0[i + nq0*p] * jac_val; //Load 2x
                    sum_kj.fma(inxmm, w0[i]); //Load 1x
                }

                sums_kj[cnt_kj] =  sum_kj; //Store 1x
            }
        }

        for (int q = 0; q < nm1-p; ++q, ++mode)
        {
            int cnt_kj = 0;

            for (int k = 0; k < nq2; ++k)
            {
                vec_t sum_k = basis1[mode*nq1] * sums_kj[cnt_kj] * w1[0]; //Load 3x
                ++cnt_kj;

                for (int j = 1; j < nq1; ++j, ++cnt_kj)
                {
                    vec_t tmp2 = basis1[mode*nq1 + j] * sums_kj[cnt_kj]; //Load 2x
                    sum_k.fma(tmp2, w1[j]); //Load 1x
                }

                sums_k[k] = sum_k; //Store 1x
            }

            for (int r = 0; r < nm2-p-q; ++r, ++mode2, ++cnt_pqr)
            {
                vec_t tmp = sums_k[0] * basis2[mode2*nq2] * w2[0]; //Load 3x

                for (int k = 1; k < nq2; ++k)
                {
                    vec_t tmp2 = sums_k[k] * basis2[mode2*nq2 + k]; //Load 2x
                    tmp.fma(tmp2, w2[k]); //Load 1x
                }

                ScaleAppend<SCALE, APPEND>(out[cnt_pqr], tmp, scale); //Store 1x
            }
        }
    }

    if (correct)
    {
        for (int k = 0, cnt = 0; k < nq2; ++k)
        {
            vec_t tmpQ2 = w2[k]; //Load 1x
            if (!DEFORMED)
            {
                tmpQ2 = tmpQ2 * jac[0];
            }

            for (int j = 0; j < nq1; ++j)
            {
                vec_t tmpQ1 = tmpQ2 * w1[j]; //Load 1x

                for (int i = 0; i < nq0; ++i, ++cnt)
                {
                    // Store jac * quadrature weight
                    vec_t tmpQ = tmpQ1 * w0[i]; //Load 1x
                    vec_t tmpIn = in[cnt]; //Load 1x

                    if (DEFORMED)
                    {
                        tmpQ = tmpQ * jac[k*nq0*nq1 + j*nq0 + i];
                    }

                    // top vertex
                    //
                    // outarray[1] += inarray[cnt] * basis2[nq2 + k] * (
                    //     basis0[i]*basis1[nq1+j] + basis0[nq0+i]*basis1[j] +
                    //     basis0[nq0+i]*basis1[nq1+j]);

                    vec_t tmp = basis0[i] * basis1[nq1+j]; //Load 2x
                    tmp.fma(basis0[nq0+i], basis1[j]); //Load 2x
                    tmp.fma(basis0[nq0+i], basis1[nq1+j]); //Load 2x
                    tmp = tmp * basis2[nq2+k]; //Load 1x
                    tmp = tmp * tmpIn;

                    // add to existing entry
                    vec_t tmpOut = tmp * tmpQ;
                    ScaleAppend<SCALE, true>(out[1], tmpOut, scale); //Store 1x

                    // bottom vertex
                    //
                    // outarray[nm2] += inarray[cnt] * basis2[k] * (
                    //    basis0[nq0+i] * basis1[nq1+j]);

                    tmp = basis0[nq0+i] * basis1[nq1+j] * basis2[k] * tmpIn; //Load 3x
                    tmpOut = tmp * tmpQ;
                    ScaleAppend<SCALE, true>(out[nm2], tmpOut, scale); //Store 1x

                    // singular edge
                    for (int r = 1; r < nm2-1; ++r)
                    {
                        // outarray[nm2+r] += inarray[cnt] *
                        //     basis2[(r+1)*nq2+k] * basis1[nq1+j] * basis0[nq0+i];
                        tmp = basis2[(r+1)*nq2+k] * basis1[nq1+j] * basis0[nq0+i] * tmpIn; //Load 3x
                        tmpOut = tmp * tmpQ;
                        ScaleAppend<SCALE, true>(out[nm2+r], tmpOut, scale); //Store 1x
                    }
                }
            }
        }
    }
}

#elif defined(SHAPE_TYPE_PRISM)

template<bool SCALE, bool APPEND, bool DEFORMED>
NEK_FORCE_INLINE static void IProductPrismKernel(
    const size_t nm0, const size_t nm1, const size_t nm2,
    const size_t nq0, const size_t nq1, const size_t nq2,
    const bool correct,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
    const std::vector<vec_t, allocator<vec_t>>& basis2,
    const std::vector<vec_t, allocator<vec_t>>& w0,
    const std::vector<vec_t, allocator<vec_t>>& w1,
    const std::vector<vec_t, allocator<vec_t>>& w2,
    const vec_t* jac,
          std::vector<vec_t, allocator<vec_t>>& sums_kj, // nq2 * nq1
          std::vector<vec_t, allocator<vec_t>>& sums_k,  // nq2
          std::vector<vec_t, allocator<vec_t>>& corr_q,  // nm1
          std::vector<vec_t, allocator<vec_t>>& out,
          NekDouble scale = 1.0)
{
    int mode_pr = 0, mode_pqr = 0;

    for (int p = 0; p < nm0; ++p)
    {
        int cnt_kji = 0, cnt_kj = 0;

        for (int k = 0; k < nq2; ++k)
        {
            for (int j = 0; j < nq1; ++j, ++cnt_kj)
            {
                vec_t sum_kj = 0.0;

                for (int i = 0; i < nq0; ++i, ++cnt_kji)
                {
                    vec_t jac_val;

                    if (DEFORMED)
                    {
                        jac_val = jac[nq0*nq1*k + nq0*j + i];
                    }
                    else{
                        jac_val = jac[0];
                    }

                    vec_t prod = basis0[nq0*p + i] * jac_val * w0[i]; // load 2x
                    vec_t fn = in[cnt_kji]; //load 1x
                    sum_kj.fma(prod, fn);
                }

                sums_kj[cnt_kj] = sum_kj; //store 1x
            }
        }

        for (int q = 0; q < nm1; ++q)
        {
            cnt_kj = 0;

            for (int k = 0; k < nq2; ++k)
            {
                vec_t sum_k = 0.0;

                for (int j = 0; j < nq1; ++j, ++cnt_kj)
                {
                    sum_k.fma(basis1[q*nq1 + j] * w1[j], sums_kj[cnt_kj]); //Load 3x
                }

                sums_k[k] = sum_k; //Store 1x
            }

            // Start with nesting. Should be able to move out of q
            // loop and sotre identical copies...
            for (int r = 0; r < nm2-p; ++r, ++mode_pqr)
            {
                vec_t sum_k = 0.0;

                for (int k = 0; k < nq2; ++k)
                {
                    sum_k.fma(basis2[(mode_pr + r)*nq2 + k] * w2[k], sums_k[k]); //Load 3x
                }

                ScaleAppend<SCALE, APPEND>(out[mode_pqr], sum_k, scale); //Store 1x
            }
        }

        mode_pr += nm2 - p;
    }

    if (correct)
    {
        // Corrections for singular edge
        for (int q = 0; q < nm1; ++q)
        {
            corr_q[q] = 0.0; //T(outptr + (nm2*q + 1)*VW);
        }

        int cnt_kji = 0;
        for (int k = 0; k < nq2; ++k)
        {
            vec_t k_weight = w2[k];
            if (!DEFORMED)
            {
                k_weight = k_weight * jac[0];
            }

            for (int j = 0; j < nq1; ++j)
            {
                vec_t kj_weight = k_weight * w1[j];
                for (int i = 0; i < nq0; ++i, ++cnt_kji)
                {

                    vec_t kji_weight = kj_weight * w0[i];
                    vec_t prod = kji_weight * in[cnt_kji];

                    if (DEFORMED)
                    {
                        prod *= jac[k*nq1*nq0 + j*nq0 + i];
                    }

                    vec_t basis_2 = basis2[nq2 + k];
                    vec_t basis_0 = basis0[nq0 + i];
                    //Add phi_1q1 to phi_0q1
                    for (int q = 0; q < nm1; ++q)
                    {
                        vec_t basis_1 = basis1[q*nq1 + j];

                        corr_q[q].fma(basis_2*basis_1, basis_0*prod);
                    }
                }
            }
        }

        for (int q = 0; q < nm1; ++q)
        {
            ScaleAppend<SCALE, true>(out[nm2*q + 1], corr_q[q], scale);
        }
    }
}

#elif defined(SHAPE_TYPE_PYR)

template<bool SCALE, bool APPEND, bool DEFORMED>
NEK_FORCE_INLINE static void IProductPyrKernel(
    const size_t nm0, const size_t nm1, const size_t nm2,
    const size_t nq0, const size_t nq1, const size_t nq2,
    const bool correct,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
    const std::vector<vec_t, allocator<vec_t>>& basis2,
    const std::vector<vec_t, allocator<vec_t>>& w0,
    const std::vector<vec_t, allocator<vec_t>>& w1,
    const std::vector<vec_t, allocator<vec_t>>& w2,
    const vec_t* jac,
          std::vector<vec_t, allocator<vec_t>>& sums_kj,
          std::vector<vec_t, allocator<vec_t>>& sums_k,
          std::vector<vec_t, allocator<vec_t>>& out,
          NekDouble scale = 1.0)
{
    int mode_pqr = 0;

    for (int p = 0; p < nm0; ++p)
    {
        int cnt_kji = 0, cnt_kj = 0;

        for (int k = 0; k < nq2; ++k)
        {
            for (int j = 0; j < nq1; ++j, ++cnt_kj)
            {
                vec_t sum_kj = 0.0;

                for (int i = 0; i < nq0; ++i, ++cnt_kji)
                {
                    vec_t jac_val;

                    if (DEFORMED)
                    {
                        jac_val = jac[nq0*nq1*k + nq0*j + i];
                    }
                    else
                    {
                        jac_val = jac[0];
                    }

                    vec_t prod = basis0[nq0*p + i] * jac_val * w0[i]; // load 2x
                    vec_t fn = in[cnt_kji]; //load 1x
                    sum_kj.fma(prod, fn);
                }

                sums_kj[cnt_kj] = sum_kj; //store 1x
            }
        }

        for (int q = 0; q < p; ++q)
        {
            cnt_kj = 0;

            for (int k = 0; k < nq2; ++k)
            {
                vec_t sum_k = 0.0;

                for (int j = 0; j < nq1; ++j, ++cnt_kj)
                {
                    sum_k.fma(basis1[q*nq1+j]*w1[j], sums_kj[cnt_kj]); //Load 3x
                }

                sums_k[k] = sum_k; //Store 1x
            }

            for (int r = 0; r < nm2-p; ++r, ++mode_pqr)
            {
                vec_t sum_k = 0.0;

                for (int k = 0; k < nq2; ++k)
                {
                    sum_k.fma(basis2[mode_pqr*nq2+k]*w2[k], sums_k[k]); //Load 3x
                }

                ScaleAppend<SCALE, APPEND>(out[mode_pqr], sum_k, scale); //Store 1x
            }
        }

        for (int q = p; q < nm1; ++q)
        {
            cnt_kj = 0;

            for (int k = 0; k < nq2; ++k)
            {
                vec_t sum_k = 0.0;

                for (int j = 0; j < nq1; ++j, ++cnt_kj)
                {
                    sum_k.fma(basis1[q*nq1+j]*w1[j], sums_kj[cnt_kj]); //Load 3x
                }

                sums_k[k] = sum_k; //Store 1x
            }

            for (int r = 0; r < nm2-q; ++r, ++mode_pqr)
            {
                vec_t sum_k = 0.0;

                for (int k = 0; k < nq2; ++k)
                {
                    sum_k.fma(basis2[mode_pqr*nq2+k]*w2[k], sums_k[k]); //Load 3x
                }

                ScaleAppend<SCALE, APPEND>(out[mode_pqr], sum_k, scale); //Store 1x
            }
        }
    }

    if (correct)
    {
        for (int k = 0, cnt = 0; k < nq2; ++k)
        {
            vec_t tmpQ2 = w2[k]; //Load 1x

            if (!DEFORMED)
            {
                tmpQ2 = tmpQ2 * jac[0];
            }

            for (int j = 0; j < nq1; ++j)
            {
                vec_t tmpQ1 = tmpQ2 * w1[j]; //Load 1x

                for (int i = 0; i < nq0; ++i, ++cnt)
                {
                    // Store jac * quadrature weight
                    vec_t tmpQ = tmpQ1 * w0[i]; //Load 1x
                    vec_t tmpIn = in[cnt]; //Load 1x

                    if (DEFORMED)
                    {
                        tmpQ = tmpQ * jac[k*nq0*nq1 + j*nq0 + i];
                    }

                    // top vertex
                    //
                    // outarray[1] += inarray[cnt] * basis2[nq2 + k] * (
                    //     basis0[i]*basis1[nq1+j] + basis0[nq0+i]*basis1[j] +
                    //     basis0[nq0+i]*basis1[nq1+j]);

                    vec_t tmp = basis0[i] * basis1[nq1+j]; //Load 2x
                    tmp.fma(basis0[nq0+i], basis1[j]); //Load 2x
                    tmp.fma(basis0[nq0+i], basis1[nq1+j]); //Load 2x
                    tmp = tmp * basis2[nq2+k]; //Load 1x
                    tmp = tmp * tmpIn;

                    // add to existing entry
                    vec_t tmpOut = tmp * tmpQ;
                    ScaleAppend<SCALE, true>(out[1], tmpOut, scale); //Store 1x
                }
            }
        }
    }
}

#endif

// Workspace - used to dynamically get the workspace size needed for
// temporary memory.
#if defined(SHAPE_DIMENSION_1D)

template<LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void IProduct1DWorkspace(
    const size_t nm0, const size_t nq0)

{
    boost::ignore_unused(SHAPE_TYPE, nm0, nq0);

    // Check preconditions
    // None
}

#elif defined(SHAPE_DIMENSION_2D)

template<LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void IProduct2DWorkspace(
    const size_t nm0, const size_t nm1,
    const size_t nq0, const size_t nq1,
          size_t& wsp0Size)
{
    boost::ignore_unused(SHAPE_TYPE, nm0, nm1, nq1);

    // Check preconditions
    ASSERTL1((SHAPE_TYPE == LibUtilities::ShapeType::Tri  && nm0 == nm1 && nq0 == nq1+1) ||
             (SHAPE_TYPE == LibUtilities::ShapeType::Quad && nm0 == nm1 && nq0 == nq1),
             "IProduct2DWorkspace: Requires homogenous points.");

    wsp0Size = std::max(wsp0Size, nq0);
}

#elif defined(SHAPE_DIMENSION_3D)

template<LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void IProduct3DWorkspace(
    const size_t nm0, const size_t nm1, const size_t nm2,
    const size_t nq0, const size_t nq1, const size_t nq2,
          size_t& wsp0Size, size_t& wsp1Size, size_t& wsp2Size)
{
    boost::ignore_unused(SHAPE_TYPE, nm0, nm1, nm2, nq0, wsp2Size);

    // Check preconditions
    ASSERTL1((SHAPE_TYPE == LibUtilities::ShapeType::Hex   &&
              nm0 == nm1 && nm0 == nm2 && nq0 == nq1   && nq0 == nq2  ) ||
             (SHAPE_TYPE == LibUtilities::ShapeType::Tet   &&
              nm0 == nm1 && nm0 == nm2 && nq0 == nq1+1 && nq0 == nq2+1) ||
             (SHAPE_TYPE == LibUtilities::ShapeType::Prism &&
              nm0 == nm1 && nm0 == nm2 && nq0 == nq1   && nq0 == nq2+1) ||
             (SHAPE_TYPE == LibUtilities::ShapeType::Pyr   &&
              nm0 == nm1 && nm0 == nm2 && nq0 == nq1   && nq0 == nq2+1),
             "IProduct3DWorkspace: Requires homogenous points.");

    wsp0Size = std::max(wsp0Size, nq1 * nq2);
    wsp1Size = std::max(wsp1Size, nq2);

#if defined(SHAPE_TYPE_PRISM)
    wsp2Size = std::max(wsp0Size, nm1);
#endif
}
#endif // SHAPE_DIMENSION

// The dimension kernels which select the shape kernel.
#if defined(SHAPE_DIMENSION_1D)

template<LibUtilities::ShapeType SHAPE_TYPE,
         bool SCALE, bool APPEND, bool DEFORMED>
NEK_FORCE_INLINE static void IProduct1DKernel(
    const size_t nm0,
    const size_t nq0,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& w0,
    const vec_t *jac,
          std::vector<vec_t, allocator<vec_t>>& out,
          NekDouble scale = 1.0)
{
#if defined(SHAPE_TYPE_SEG)
  IProductSegKernel<SCALE, APPEND, DEFORMED>(nm0, nq0,
                    in, basis0, w0, jac, out, scale);
#endif
}

#elif defined(SHAPE_DIMENSION_2D)

template<LibUtilities::ShapeType SHAPE_TYPE,
         bool SCALE, bool APPEND, bool DEFORMED>
NEK_FORCE_INLINE static void IProduct2DKernel(
    const size_t nm0, const size_t nm1,
    const size_t nq0, const size_t nq1,
    const bool correct,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
    const std::vector<vec_t, allocator<vec_t>>& w0,
    const std::vector<vec_t, allocator<vec_t>>& w1,
    const vec_t* jac,
          std::vector<vec_t, allocator<vec_t>>& wsp0,
          std::vector<vec_t, allocator<vec_t>>& out,
          NekDouble scale = 1.0)
{
#if defined(SHAPE_TYPE_TRI)
    IProductTriKernel<SCALE, APPEND, DEFORMED>(
                     nm0, nm1, nq0, nq1, correct,
                      in, basis0, basis1, w0, w1, jac, wsp0, out, scale);
#elif defined(SHAPE_TYPE_QUAD)
    boost::ignore_unused(correct);
    IProductQuadKernel<SCALE, APPEND, DEFORMED>(
                       nm0, nm1, nq0, nq1,
                       in, basis0, basis1, w0, w1, jac, wsp0, out, scale);
#endif
}

#elif defined(SHAPE_DIMENSION_3D)

template<LibUtilities::ShapeType SHAPE_TYPE,
         bool SCALE, bool APPEND, bool DEFORMED>
NEK_FORCE_INLINE static void IProduct3DKernel(
    const size_t nm0, const size_t nm1, const size_t nm2,
    const size_t nq0, const size_t nq1, const size_t nq2,
    const bool correct,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
    const std::vector<vec_t, allocator<vec_t>>& basis2,
    const std::vector<vec_t, allocator<vec_t>>& w0,
    const std::vector<vec_t, allocator<vec_t>>& w1,
    const std::vector<vec_t, allocator<vec_t>>& w2,
    const vec_t* jac,
          std::vector<vec_t, allocator<vec_t>>& wsp0,
          std::vector<vec_t, allocator<vec_t>>& wsp1,
          std::vector<vec_t, allocator<vec_t>>& wsp2,
          std::vector<vec_t, allocator<vec_t>>& out,
          NekDouble scale = 1.0)
{
#if defined(SHAPE_TYPE_HEX)
    boost::ignore_unused(correct, wsp2);
    IProductHexKernel<SCALE, APPEND, DEFORMED>(
                      nm0, nm1, nm2, nq0, nq1, nq2,
                      in, basis0, basis1, basis2, w0, w1, w2, jac,
                      wsp0, wsp1, out, scale);
#elif defined(SHAPE_TYPE_TET)
    boost::ignore_unused(wsp2);
    IProductTetKernel<SCALE, APPEND, DEFORMED>(
                     nm0, nm1, nm2, nq0, nq1, nq2,
                     correct,
                     in, basis0, basis1, basis2, w0, w1, w2, jac,
                     wsp0, wsp1, out, scale);
#elif defined(SHAPE_TYPE_PRISM)
    IProductPrismKernel<SCALE, APPEND, DEFORMED>(
                        nm0, nm1, nm2, nq0, nq1, nq2,
                        correct,
                        in, basis0, basis1, basis2, w0, w1, w2, jac,
                        wsp0, wsp1, wsp2, out, scale);
#elif defined(SHAPE_TYPE_PYR)
    boost::ignore_unused(wsp2);
    IProductPyrKernel<SCALE, APPEND, DEFORMED>(
                      nm0, nm1, nm2, nq0, nq1, nq2,
                      correct,
                      in, basis0, basis1, basis2, w0, w1, w2, jac,
                      wsp0, wsp1, out, scale);
#endif
}

#endif // SHAPE_DIMENSION

} // namespace MatrixFree
} // namespace Nektar

#endif
