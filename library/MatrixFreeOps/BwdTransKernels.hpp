#ifndef NEKTAR_LIBRARY_MF_BWD_TRANS_KERNELS_HPP
#define NEKTAR_LIBRARY_MF_BWD_TRANS_KERNELS_HPP

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
#if defined(SHAPE_TYPE_SEG)

NEK_FORCE_INLINE static void BwdTransSegKernel(
    const size_t nm0, const size_t nq0,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
          std::vector<vec_t, allocator<vec_t>>& out)
{
    for (int i = 0; i < nq0; ++i)
    {
        vec_t tmp =  in[0] * basis0[i]; //Load 2x

        for (int p = 1; p < nm0; ++p)
        {
            tmp.fma(in[p], basis0[p * nq0 + i]); //Load 2x
        }

        out[i] = tmp; //Store 1x
    }
}

#elif defined(SHAPE_TYPE_TRI)

NEK_FORCE_INLINE static void BwdTransTriKernel(
    const size_t nm0, const size_t nm1,
    const size_t nq0, const size_t nq1,
    const bool correct,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
          std::vector<vec_t, allocator<vec_t>>& p_sums, // nm0
          std::vector<vec_t, allocator<vec_t>>& out)
{
    for (int eta1 = 0, eta_idx = 0; eta1 < nq1; ++eta1)
    {
        for (int p = 0, mode = 0; p < nm0; ++p)
        {
            vec_t p_sum = 0.0;

            for (int q = 0; q < (nm1-p); ++q, ++mode)
            {
               p_sum.fma(basis1[mode * nq1 + eta1], in[mode]);
            }

            p_sums[p] = p_sum; //Store 1x
        }

        // Already have q_sum at each quadrature point in eta 1 for
        // each mode, p.  From this assemble the tensor produce of
        // each quadrature point, eta1
        for (int eta0 = 0; eta0 < nq0; ++eta0, ++eta_idx)
        {
            vec_t p_sum = 0.0;
            for (int p = 0; p < nm0; ++p)
            {
                p_sum.fma(p_sums[p], basis0[p*nq0 + eta0]); //Load 2x
            }

            if (correct)
            {
                //p_sum += coef * basis0 * basis1
                p_sum.fma(in[1] * basis0[nq0 + eta0], basis1[nq1 + eta1]);
            }

            out[eta_idx] = p_sum;
        }
    }
}

#elif defined(SHAPE_TYPE_QUAD)

NEK_FORCE_INLINE static void BwdTransQuadKernel(
    const size_t nm0, const size_t nm1,
    const size_t nq0, const size_t nq1,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
          std::vector<vec_t, allocator<vec_t>>& wsp, // nq0 * nm1
          std::vector<vec_t, allocator<vec_t>>& out)
{
    for (int i = 0, cnt_iq = 0; i < nq0; ++i)
    {
        for (int q = 0, cnt_pq = 0; q < nm1; ++q, ++cnt_iq)
        {
            vec_t tmp = in[cnt_pq] * basis0[i]; //Load 2x
            ++cnt_pq;
            for (int p = 1; p < nm0; ++p, ++cnt_pq)
            {
                tmp.fma(in[cnt_pq], basis0[p * nq0 + i]); //Load 2x
            }
            wsp[cnt_iq] = tmp; //Store 1x
        }
    }

    for (int j = 0, cnt_ij = 0; j < nq1; ++j)
    {
        for (int i = 0, cnt_iq = 0; i < nq0; ++i, ++cnt_ij)
        {
            vec_t tmp = wsp[cnt_iq] * basis1[j]; //Load 2x
            ++cnt_iq;
            for (int q = 1; q < nm1; ++q, ++cnt_iq)
            {
                tmp.fma(wsp[cnt_iq], basis1[q * nq1 + j]); //Load 2x
            }
            out[cnt_ij] = tmp; //Store 1x
        }
    }
}

#elif defined(SHAPE_TYPE_HEX)

NEK_FORCE_INLINE static void BwdTransHexKernel(
    const size_t nm0, const size_t nm1, const size_t nm2,
    const size_t nq0, const size_t nq1, const size_t nq2,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
    const std::vector<vec_t, allocator<vec_t>>& basis2,
          std::vector<vec_t, allocator<vec_t>>& sum_irq, // nq0 * nm2 * nm1
          std::vector<vec_t, allocator<vec_t>>& sum_jir, // nq1 * nq0 * nm2
          std::vector<vec_t, allocator<vec_t>>& out)
{
    for (int i = 0, cnt_irq = 0; i < nq0; ++i)
    {
        for (int r = 0, cnt_rqp = 0; r < nm2; ++r)
        {
            for (int q = 0; q < nm1; ++q, ++cnt_irq)
            {
                vec_t tmp = in[cnt_rqp] * basis0[i];
                ++cnt_rqp;

                for (int p = 1; p < nm0; ++p, ++cnt_rqp)
                {
                    tmp.fma(in[cnt_rqp], basis0[p*nq0+i]);
                }

                sum_irq[cnt_irq] = tmp;
            }
        }
    }

    for (int j = 0, cnt_jir = 0; j < nq1; ++j)
    {
        for (int i = 0, cnt_irq = 0; i < nq0; ++i)
        {
            for (int r = 0; r < nm2; ++r, ++cnt_jir)
            {
                vec_t tmp = sum_irq[cnt_irq]* basis1[j];
                ++cnt_irq;

                for (int q = 1; q < nm1; ++q)
                {
                    tmp.fma(sum_irq[cnt_irq++], basis1[q*nq1+j]);
                }

                sum_jir[cnt_jir] = tmp;
            }
        }
    }

    for (int k = 0, cnt_kji = 0; k < nq2; ++k)
    {
        for (int j = 0, cnt_jir = 0; j < nq1; ++j)
        {
            for (int i = 0; i < nq0; ++i, ++cnt_kji)
            {
                vec_t tmp = sum_jir[cnt_jir] * basis2[k];
                ++cnt_jir;

                for (int r = 1; r < nm2; ++r)
                {
                    tmp.fma(sum_jir[cnt_jir++], basis2[r*nq2+k]);
                }

                out[cnt_kji] = tmp;
            }
        }
    }
}


#elif defined(SHAPE_TYPE_TET)

NEK_FORCE_INLINE static void BwdTransTetKernel(
    const size_t nm0, const size_t nm1, const size_t nm2,
    const size_t nq0, const size_t nq1, const size_t nq2,
    const bool correct,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
    const std::vector<vec_t, allocator<vec_t>>& basis2,
          std::vector<vec_t, allocator<vec_t>>& fpq, // nm0 * nm1
          std::vector<vec_t, allocator<vec_t>>& fp,  // nm0
          std::vector<vec_t, allocator<vec_t>>& out)
{
    for (int k = 0, cnt_kji = 0; k < nq2; ++k)
    {
        int cnt_pq = 0, mode = 0;
        for (int p = 0; p < nm0; ++p)
        {
            for (int q = 0; q < nm1-p; ++q, ++cnt_pq)
            {
                vec_t prod = in[mode] * basis2[k + nq2*mode]; //Load 2x
                ++mode;

                for (int r = 1; r < nm2-p-q; ++r, ++mode)
                {
                    vec_t inxmm = in[mode]; //Load 1x
                    prod.fma(inxmm, basis2[k + nq2*mode]); //Load 1x
                }

                fpq[cnt_pq] = prod; //Store 1x
            }

            //increment mode in case order1!=order2
            for(int q = nm1-p; q < nm2-p; ++q)
            {
                mode += nm2-p-q;
            }
        }

        for (int j = 0; j < nq1; ++j)
        {
            mode = cnt_pq = 0;
            for (int p = 0; p < nm0; ++p)
            {
                vec_t prod = fpq[cnt_pq] * basis1[mode*nq1+j]; //Load 2x
                ++cnt_pq;

                for (int q = 1; q < nm1 - p; ++q, ++cnt_pq)
                {
                    prod.fma(fpq[cnt_pq], basis1[(mode+q)*nq1+j]); //Load 2x
                }

                fp[p] = prod; //Store 1x
                mode += nm1 - p;
            }

            for (int i = 0; i < nq0; ++i, ++cnt_kji)
            {
                vec_t tmp = basis0[i] * fp[0]; //Load 2x

                for (int p = 1; p < nm0; ++p)
                {
                    tmp.fma(basis0[p*nq0+i], fp[p]); //Load 2x
                }

                if (correct)
                {
                    // top vertex
                    //
                    // sum += inarray[1] * base2[nquad2 + k] * (
                    //     base0[i] * base1[nquad1+j] +
                    //     base0[nquad0+i] * base1[j] +
                    //     base0[nquad0+i] * base1[nquad1+j]);

                    vec_t tmp1 = basis0[i] * basis1[nq1+j]; //Load 2x
                    tmp1.fma(basis0[nq0+i], basis1[j]); //Load 2x
                    tmp1.fma(basis0[nq0+i], basis1[nq1+j]); //Load 2x
                    tmp1 = tmp1 * basis2[nq2+k]; //Load 1x

                    vec_t inarray1 = in[1]; //Load 1x
                    tmp.fma(tmp1, inarray1);

                    // bottom vertex
                    //
                    // sum += inarray[order2] * base2[k] * (
                    //     base0[nquad0+i] * base1[nquad1+j]);
                    tmp1 = basis0[nq0+i] * basis1[nq1+j]; //Load 2x
                    tmp1 = tmp1 * basis2[k]; //Load 1x
                    inarray1 = in[nm2]; //Load 1x
                    tmp.fma(inarray1, tmp1);

                    // singular edge
                    for (int r = 1; r < nm2-1; ++r)
                    {
                        // sum += inarray[order2+r] * base2[(r+1)*nquad2+k] *
                        //     base1[nquad1+j] * base0[nquad0+i];
                        tmp1 = basis1[nq1+j] * basis0[nq0+i]; //Load 2x
                        tmp1 = tmp1 * basis2[(r+1)*nq2+k]; //Load 1x
                        inarray1 = in[nm2+r]; //Load 1x
                        tmp.fma(inarray1, tmp1);
                        // multiply by (1-a)/2
                    }
                }

                out[cnt_kji] = tmp; //Store 1x
            }
        }
    }
}

#elif defined(SHAPE_TYPE_PRISM)

NEK_FORCE_INLINE static void BwdTransPrismKernel(
    const size_t nm0, const size_t nm1, const size_t nm2,
    const size_t nq0, const size_t nq1, const size_t nq2,
    const bool correct,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
    const std::vector<vec_t, allocator<vec_t>>& basis2,
          std::vector<vec_t, allocator<vec_t>>& fpq, // nm0 * nm1
          std::vector<vec_t, allocator<vec_t>>& fp,  // nm0
          std::vector<vec_t, allocator<vec_t>>& out)
{
    for (int k = 0, cnt_kji = 0; k < nq2; ++k)
    {
        int mode_pqr = 0, mode_pq = 0, mode_pr = 0;
        for (int p = 0; p < nm0; ++p)
        {
            for (int q = 0; q < nm1; ++q, ++mode_pq)
            {
                vec_t prod = 0.0;
                for (int r = 0; r < nm2-p; ++r, ++mode_pqr)
                {
                    vec_t coef = in[mode_pqr]; //Load 1x
                    prod.fma(coef, basis2[(mode_pr + r)*nq2 + k]); //Load 1x
                }

                fpq[mode_pq] = prod; //Store 1x
            }

            mode_pr += nm2 - p;
        }

        for (int j = 0; j < nq1; ++j)
        {
            mode_pq = 0;
            for (int p = 0; p < nm0; ++p)
            {
                vec_t prod = 0.0;
                for (int q = 0; q < nm1; ++q, ++mode_pq)
                {
                    prod.fma(fpq[mode_pq], basis1[q*nq1 + j]); //Load 2x
                }
                fp[p] = prod; //Store 1x
            }

            for (int i = 0; i < nq0; ++i, ++cnt_kji)
            {

                vec_t val_kji = 0.0;
                for (int p = 0; p < nm0; ++p)
                {
                    val_kji.fma(fp[p], basis0[p*nq0 + i]); //Load 2x
                }

                if (correct)
                {
                    vec_t basis_2 = basis2[nq2 + k]; //Load 1x
                    vec_t basis_0 = basis0[nq0 + i]; //Load 1x

                    for (int q = 0; q < nm1; ++q)
                    {
                        vec_t coef_0q1 = in[q*nm2 + 1]; //Load 1x
                        vec_t basis_1 = basis1[q*nq1 + j]; //Load 1x
                        val_kji.fma(basis_2*basis_1, basis_0*coef_0q1);
                    }
                }
                out[cnt_kji] = val_kji; //store 1x
            }
        }
    }
}

#elif defined(SHAPE_TYPE_PYR)

NEK_FORCE_INLINE static void BwdTransPyrKernel(
    const size_t nm0, const size_t nm1, const size_t nm2,
    const size_t nq0, const size_t nq1, const size_t nq2,
    const bool correct,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
    const std::vector<vec_t, allocator<vec_t>>& basis2,
          std::vector<vec_t, allocator<vec_t>>& fpq, // nm0 * nm1
          std::vector<vec_t, allocator<vec_t>>& fp,  // nm0
          std::vector<vec_t, allocator<vec_t>>& out)
{
    for (int k = 0, cnt_kji = 0; k < nq2; ++k)
    {
        int mode_pqr = 0, mode_pq = 0;
        for (int p = 0; p < nm0; ++p)
        {
            for (int q = 0; q < p; ++q, ++mode_pq)
            {
                vec_t prod = 0.0;
                for (int r = 0; r < nm2-p; ++r, ++mode_pqr)
                {
                    vec_t coef = in[mode_pqr]; //Load 1x
                    prod.fma(coef, basis2[mode_pqr*nq2 + k]); //Load 1x
                }
                fpq[mode_pq] = prod; //Store 1x
            }

            for (int q = p; q < nm1; ++q, ++mode_pq)
            {
                vec_t prod = 0.0;
                for (int r = 0; r < nm2-q; ++r, ++mode_pqr)
                {
                    vec_t coef = in[mode_pqr]; //Load 1x
                    prod.fma(coef, basis2[mode_pqr*nq2 + k]); //Load 1x
                }

                fpq[mode_pq] = prod; //Store 1x
            }

            //increment mode in case nm2>nm1
            for(int q = nm1; q < nm2-p; ++q)
            {
                mode_pqr += nm2-q;
            }
        }

        for (int j = 0; j < nq1; ++j)
        {
            mode_pq = 0;
            for (int p = 0; p < nm0; ++p)
            {
                vec_t prod = 0.0;
                for (int q = 0; q < nm1; ++q, ++mode_pq)
                {
                    prod.fma(fpq[mode_pq], basis1[q*nq1 + j]); //Load 2x
                }
                fp[p] = prod; //Store 1x
            }

            for (int i = 0; i < nq0; ++i, ++cnt_kji)
            {
                vec_t val_kji = 0.0;
                for (int p = 0; p < nm0; ++p)
                {
                    val_kji.fma(fp[p], basis0[p*nq0 + i]); //Load 2x
                }

                if (correct)
                {
                    // top vertex
                    //
                    // sum += inarray[1] * base2[nquad2 + k] * (
                    //     base0[i] * base1[nquad1+j] +
                    //     base0[nquad0+i] * base1[j] +
                    //     base0[nquad0+i] * base1[nquad1+j]);
                    vec_t tmp1 = basis0[i] * basis1[nq1+j]; //Load 2x
                    tmp1.fma(basis0[nq0+i], basis1[j]); //Load 2x
                    tmp1.fma(basis0[nq0+i], basis1[nq1+j]); //Load 2x
                    tmp1 = tmp1 * basis2[nq2+k]; //Load 1x

                    vec_t inarray1 = in[1]; //Load 1x
                    val_kji.fma(tmp1, inarray1);
                }
                out[cnt_kji] = val_kji; //store 1x
            }
        }
    }
}

#endif // SHAPE_TYPE

// Workspace - used to dynamically get the workspace size needed for
// temporary memory.
#if defined(SHAPE_DIMENSION_1D)

template<LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void BwdTrans1DWorkspace(
    const size_t nm0, const size_t nq0)

{
    boost::ignore_unused(SHAPE_TYPE, nm0, nq0);
}

#elif defined(SHAPE_DIMENSION_2D)

template<LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void BwdTrans2DWorkspace(
    const size_t nm0, const size_t nm1,
    const size_t nq0, const size_t nq1,
          size_t& wsp0Size)
{
    boost::ignore_unused(SHAPE_TYPE, nm0, nm1, nq0, nq1);

    // Check preconditions
    ASSERTL1((SHAPE_TYPE == LibUtilities::ShapeType::Tri  && nm0 == nm1 && nq0 == nq1+1) ||
             (SHAPE_TYPE == LibUtilities::ShapeType::Quad && nm0 == nm1 && nq0 == nq1),
             "BwdTrans2DWorkspace: Requires homogenous points.");

#if defined(SHAPE_TYPE_TRI)
    wsp0Size = std::max(wsp0Size, nm0);
#elif defined(SHAPE_TYPE_QUAD)
    wsp0Size = std::max(wsp0Size, nm0 * nq0);
#endif
}

#elif defined(SHAPE_DIMENSION_3D)

template<LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void BwdTrans3DWorkspace(
    const size_t nm0, const size_t nm1, const size_t nm2,
    const size_t nq0, const size_t nq1, const size_t nq2,
          size_t& wsp0Size, size_t& wsp1Size)
{
    boost::ignore_unused(SHAPE_TYPE, nm0, nm1, nm2, nq0, nq1, nq2);

    // Check preconditions
    ASSERTL1((SHAPE_TYPE == LibUtilities::ShapeType::Hex   &&
              nm0 == nm1 && nm0 == nm2 && nq0 == nq1   && nq0 == nq2  ) ||
             (SHAPE_TYPE == LibUtilities::ShapeType::Tet   &&
              nm0 == nm1 && nm0 == nm2 && nq0 == nq1+1 && nq0 == nq2+1) ||
             (SHAPE_TYPE == LibUtilities::ShapeType::Pyr   &&
              nm0 == nm1 && nm0 == nm2 && nq0 == nq1   && nq0 == nq2+1) ||
             (SHAPE_TYPE == LibUtilities::ShapeType::Prism &&
              nm0 == nm1 && nm0 == nm2 && nq0 == nq1   && nq0 == nq2+1),
             "BwdTrans3DWorkspace: Requires homogenous points.");

#if defined(SHAPE_TYPE_HEX)
    wsp0Size = std::max(wsp0Size, nq0 * nm1 * nm2); // nm1 == nm2
    wsp1Size = std::max(wsp1Size, nq0 * nq1 * nm2); // nq0 == nq1
#elif defined(SHAPE_TYPE_TET) || defined(SHAPE_TYPE_PRISM) || defined(SHAPE_TYPE_PYR)
    wsp0Size = std::max(wsp0Size, nm0 * nm1); // nm0 == nm1 == nm2
    wsp1Size = std::max(wsp1Size, nm0);
#endif
}

#endif // SHAPE_DIMENSION

// The dimension kernels which select the shape kernel.
#if defined(SHAPE_DIMENSION_1D)

template<LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void BwdTrans1DKernel(
    const size_t nm0, const size_t nq0,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
          std::vector<vec_t, allocator<vec_t>>& out)
{
#if defined(SHAPE_TYPE_SEG)
    BwdTransSegKernel(nm0, nq0, in, basis0, out);
#endif
}

#elif defined(SHAPE_DIMENSION_2D)

template<LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void BwdTrans2DKernel(
    const size_t nm0, const size_t nm1,
    const size_t nq0, const size_t nq1,
    const bool correct,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
          std::vector<vec_t, allocator<vec_t>>& wsp0,
          std::vector<vec_t, allocator<vec_t>>& out)
{
#if defined(SHAPE_TYPE_TRI)
    BwdTransTriKernel(nm0, nm1, nq0, nq1, correct,
                      in, basis0, basis1,
                      wsp0, out);

#elif defined(SHAPE_TYPE_QUAD)
    boost::ignore_unused(correct);
    BwdTransQuadKernel(nm0, nm1, nq0, nq1,
                       in, basis0, basis1,
                       wsp0, out);
#endif
}

#elif defined(SHAPE_DIMENSION_3D)

template<LibUtilities::ShapeType SHAPE_TYPE>
NEK_FORCE_INLINE static void BwdTrans3DKernel(
    const size_t nm0, const size_t nm1, const size_t nm2,
    const size_t nq0, const size_t nq1, const size_t nq2,
    const bool correct,
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& basis0,
    const std::vector<vec_t, allocator<vec_t>>& basis1,
    const std::vector<vec_t, allocator<vec_t>>& basis2,
          std::vector<vec_t, allocator<vec_t>>& wsp0,
          std::vector<vec_t, allocator<vec_t>>& wsp1,
          std::vector<vec_t, allocator<vec_t>>& out)
{
#if defined(SHAPE_TYPE_HEX)
    boost::ignore_unused(correct);
    BwdTransHexKernel(nm0, nm1, nm2, nq0, nq1, nq2,
                      in, basis0, basis1, basis2,
                      wsp0, wsp1, out);
#elif defined(SHAPE_TYPE_TET)
    BwdTransTetKernel(nm0, nm1, nm2, nq0, nq1, nq2, correct,
                      in, basis0, basis1, basis2,
                      wsp0, wsp1, out);
#elif defined(SHAPE_TYPE_PRISM)
    BwdTransPrismKernel(nm0, nm1, nm2, nq0, nq1, nq2, correct,
                        in, basis0, basis1, basis2,
                        wsp0, wsp1, out);
#elif defined(SHAPE_TYPE_PYR)
    BwdTransPyrKernel(nm0, nm1, nm2, nq0, nq1, nq2, correct,
                      in, basis0, basis1, basis2,
                      wsp0, wsp1, out);
#endif
}

#endif // SHAPE_DIMENSION

} // namespace MatrixFree
} // namespace Nektar

#endif
