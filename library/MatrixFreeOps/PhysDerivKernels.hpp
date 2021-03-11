#ifndef MF_PHYSDERIVKERNELS_HPP
#define MF_PHYSDERIVKERNELS_HPP

#include <LibUtilities/SimdLib/tinysimd.hpp>

namespace Nektar
{
namespace MatrixFree
{
using namespace tinysimd;
using vec_t = simd<NekDouble>;

template<int NUMQUAD0>
inline static void PhysDerivTensor1DKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    std::vector<vec_t, allocator<vec_t>> &outptr_d0)
{
    constexpr auto nq0 = NUMQUAD0;

    //All matricies are column major ordered since operators used to
    //be computed via BLAS.

    //D0 * in
    for (int i = 0; i < nq0; ++i)
    { //row index of D0 matrix
        vec_t prod_sum = 0.0;
        for (int k = 0; k < nq0; ++k)
        { //Col index of D0, row index of IN
            vec_t v1 = D0[k * nq0 + i]; //Load 1x
            vec_t v2 = in[k]; //Load 1x
            
            prod_sum.fma(v1, v2);
        }
        outptr_d0[i] = prod_sum; //Store 1x
    }
}

template<int NUMQUAD0, int NUMQUAD1>
inline static void PhysDerivTensor2DKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const std::vector<vec_t, allocator<vec_t>> &D1,
          std::vector<vec_t, allocator<vec_t>> &outptr_d0,
          std::vector<vec_t, allocator<vec_t>> &outptr_d1)
{
    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1;

    //All matricies are column major ordered since operators used to be computed via BLAS.

    //D0 * in
    for (int i = 0; i < nq0; ++i)
    { //row index of D0 matrix
        for (int j = 0; j < nq1; ++j)
        { //col index of IN matrix

            vec_t prod_sum = 0.0;
            for (int k = 0; k < nq0; ++k)
            { //Col index of D0, row index of IN
                vec_t v1 = D0[k * nq0 + i]; //Load 1x
                vec_t v2 = in[j * nq0 + k]; //Load 1x

                prod_sum.fma(v1, v2);
            }

            outptr_d0[j * nq0 + i] = prod_sum; //Store 1x
        }
    }

    //in * D1^T
    for (int i = 0; i < nq0; ++i)
    { //row index for grid
        for (int j = 0; j < nq1; ++j)
        { //Column index for D1^T (row idx for D1)

            vec_t prod_sum = 0.0;
            for(int k = 0; k < nq1; ++k)
            {
                vec_t v1 = in[k * nq0 + i]; //Load 1x
                vec_t v2 = D1[k * nq1 + j]; //Load 1x

                prod_sum.fma(v1, v2);
            }

            outptr_d1[j * nq0 + i] = prod_sum; //Store 1x

        }
    }

}

// Overloaded with diffusion
template<int NUMQUAD0, int NUMQUAD1>
inline static void PhysDerivTensor2DKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const std::vector<vec_t, allocator<vec_t>> &D1,
          std::vector<vec_t, allocator<vec_t>> &outptr_d0,
          std::vector<vec_t, allocator<vec_t>> &outptr_d1,
    const Array<OneD, Array<OneD, NekDouble> > &diff)
{
    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1, nqTot = NUMQUAD0 * NUMQUAD1;
    vec_t d00 = diff[0][0], d01 = diff[0][1], d11 = diff[1][1];
    
    std::vector<vec_t, allocator<vec_t>> temp_d0 (nqTot), temp_d1 (nqTot);
    

    //All matricies are column major ordered since operators used to be computed via BLAS.

    //D0 * in
    for (int i = 0; i < nq0; ++i)
    { //row index of D0 matrix
        for (int j = 0; j < nq1; ++j)
        { //col index of IN matrix

            vec_t prod_sum = 0.0;
            for (int k = 0; k < nq0; ++k)
            { //Col index of D0, row index of IN
                vec_t v1 = D0[k * nq0 + i]; //Load 1x
                vec_t v2 = in[j * nq0 + k]; //Load 1x

                prod_sum.fma(v1, v2);
            }

            temp_d0[j * nq0 + i] = prod_sum; //Store 1x
        }
    }

    //in * D1^T
    for (int i = 0; i < nq0; ++i)
    { //row index for grid
        for (int j = 0; j < nq1; ++j)
        { //Column index for D1^T (row idx for D1)

            vec_t prod_sum = 0.0;
            for(int k = 0; k < nq1; ++k)
            {
                vec_t v1 = in[k * nq0 + i]; //Load 1x
                vec_t v2 = D1[k * nq1 + j]; //Load 1x

                prod_sum.fma(v1, v2);
            }

            temp_d1[j * nq0 + i] = prod_sum; //Store 1x

        }
    }
    
    for (int ij = 0; ij < nqTot; ++ij) {
        vec_t prod_sum = 0.0;
        vec_t temp0 = temp_d0[ij];
        prod_sum.fma(d00,temp0);
        vec_t temp1 = temp_d1[ij];
        prod_sum.fma(d01,temp1);
        outptr_d0[ij] = prod_sum;
    }

    for (int ij = 0; ij < nqTot; ++ij) {
        vec_t prod_sum = 0.0;
        vec_t temp0 = temp_d0[ij];
        prod_sum.fma(d01,temp0);
        vec_t temp1 = temp_d1[ij];
        prod_sum.fma(d11,temp1);
        outptr_d1[ij] = prod_sum;
    }

}

// Diffusion kernel 2D
# if 0
template<int NUMQUAD0, int NUMQUAD1>
inline static void PhysDiffusion2DKernel(
    const std::vector<vec_t, allocator<vec_t>> &inptr_d0,
    const std::vector<vec_t, allocator<vec_t>> &inptr_d1,
          std::vector<vec_t, allocator<vec_t>> &outptr_d0,
          std::vector<vec_t, allocator<vec_t>> &outptr_d1,
          Array<OneD, Array<OneD, NekDouble> > diff)
{
    constexpr auto nqTot = NUMQUAD0 * NUMQUAD1;
    NekDouble d00 = diff[0][0], d01 = diff[0][1], d11 = diff[1][1];
    for (int i; i < nqTot; i++) {
        vec_t D0 = inptr_d0[i];
        vec_t D1 = inptr_d1[i];
        vec_t prod_sum0, prod_sum1;
        prod_sum0.fma(d00,D0);
        prod_sum0.fma(d01,D1);
        prod_sum1.fma(d01,D0);
        prod_sum1.fma(d11,D1);
        outptr_d0[i] = prod_sum0;
        outptr_d1[i] = prod_sum1;
    }

}
#endif


#if 0 
template<int NUMQUAD0, bool DEFORMED>
static void PhysDerivSegKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const vec_t* df_ptr,
    std::vector<vec_t, allocator<vec_t>> &out_d0)
{

    constexpr auto nq0 = NUMQUAD0;

    PhysDerivTensor1DKernel<NUMQUAD0>(in, D0, out_d0);
    //Results written to out_d0

    vec_t df0;
    if (!DEFORMED)
    {
        df0 = df_ptr[0];
    }

    for (int j = 0; j < nq0; ++j)
    {
        if (DEFORMED)
        {
            df0 = df_ptr[j];
        }
        //Multiply by derivative factors
        out_d0[j] *= df0; //Store 1x
    }
}
#endif
    
template<int NUMQUAD0, int NUMQUAD1, bool DEFORMED>
static void PhysDerivQuadKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &Z0,
    const std::vector<vec_t, allocator<vec_t>> &Z1,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const std::vector<vec_t, allocator<vec_t>> &D1,
    const vec_t* df_ptr,
          std::vector<vec_t, allocator<vec_t>> &out_d0,
          std::vector<vec_t, allocator<vec_t>> &out_d1)
{
    boost::ignore_unused(Z0, Z1);

    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1;
    constexpr auto ndf = 4;

    PhysDerivTensor2DKernel<NUMQUAD0, NUMQUAD1>(
        in,
        D0, D1,
        out_d0, out_d1); //Results written to out_d0, out_d1

    vec_t df0, df1, df2, df3;
    if (!DEFORMED)
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
    }

    int cnt_ji = 0;
    for (int j = 0; j < nq0; ++j)
    {
        for (int i = 0; i < nq1; ++i, ++cnt_ji)
        {

            vec_t d0 = out_d0[cnt_ji];// Load 1x
            vec_t d1 = out_d1[cnt_ji]; //Load 1x

            if (DEFORMED)
            {
                df0 = df_ptr[cnt_ji * ndf];
                df1 = df_ptr[cnt_ji * ndf + 1];
                df2 = df_ptr[cnt_ji * ndf + 2];
                df3 = df_ptr[cnt_ji * ndf + 3];
            }

            //Multiply by derivative factors
            vec_t out0 = d0 * df0; //d0 * df0 + d1 * df1
            out0.fma(d1, df1);
            out_d0[cnt_ji] = out0; //Store 1x

            vec_t out1 = d0 * df2; //d0 * df2 + d1 * df3
            out1.fma(d1, df3);
            out_d1[cnt_ji] = out1; //Store 1x
        }
    }

}

template<int NUMQUAD0, int NUMQUAD1, bool DEFORMED>
static void PhysDerivTriKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &Z0,
    const std::vector<vec_t, allocator<vec_t>> &Z1,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const std::vector<vec_t, allocator<vec_t>> &D1,
    const vec_t* df_ptr,
    std::vector<vec_t, allocator<vec_t>> &out_d0,
    std::vector<vec_t, allocator<vec_t>> &out_d1)
{

    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1;
    constexpr auto ndf = 4;

    PhysDerivTensor2DKernel<NUMQUAD0, NUMQUAD1>(
        in,
        D0, D1,
        out_d0, out_d1); //Results written to out_d0, out_d1

    vec_t df0, df1, df2, df3;
    if (!DEFORMED) //Optimized out by compiler
    {
        df0 = df_ptr[0];
        df1 = df_ptr[1];
        df2 = df_ptr[2];
        df3 = df_ptr[3];
    }

    int cnt_ji = 0;
    for (int j = 0; j < nq1; ++j)
    {

        vec_t xfrm0 = 2.0 / (1.0 - Z1[j]); //Load 1x

        for (int i = 0; i < nq0; ++i, ++cnt_ji)
        {

            if (DEFORMED)
            {
                df0 = df_ptr[cnt_ji * ndf];
                df1 = df_ptr[cnt_ji * ndf + 1];
                df2 = df_ptr[cnt_ji * ndf + 2];
                df3 = df_ptr[cnt_ji * ndf + 3];
            }

            //moving from standard to collapsed coordinates
            vec_t d0 = xfrm0 * out_d0[cnt_ji]; //Load 1x
            vec_t d1 = out_d1[cnt_ji]; //Load 1x
            vec_t xfrm1 = 0.5 * (1.0 + Z0[i]); //Load 1x
            d1.fma(d0, xfrm1);

            //Multiply by derivative factors
            vec_t out0 = d0 * df0; //d0 * df0 + d1 * df10
            out0.fma(d1, df1);
            out_d0[cnt_ji] = out0; // store 1x

            vec_t out1 = d0 * df2; //d0 * df2 + d1 * df3
            out1.fma(d1, df3);
            out_d1[cnt_ji] = out1; //Store 1x

        }
    }

}

template<int NUMQUAD0, int NUMQUAD1, int NUMQUAD2>
inline static void PhysDerivTensor3DKernel(
    const std::vector<vec_t, allocator<vec_t>> &in,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const std::vector<vec_t, allocator<vec_t>> &D1,
    const std::vector<vec_t, allocator<vec_t>> &D2,
          std::vector<vec_t, allocator<vec_t>> &out_d0,
          std::vector<vec_t, allocator<vec_t>> &out_d1,
          std::vector<vec_t, allocator<vec_t>> &out_d2)
{
    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;

    //Direction 1
    for (int i = 0; i < nq0; ++i)
    {
        for (int j = 0; j < nq1*nq2; ++j)
        {
            vec_t prod_sum = 0.0;
            for (int k = 0; k < nq0; ++k)
            {
                vec_t v1 = D0[k*nq0+i];
                vec_t v2 = in[j*nq0+k];

                prod_sum.fma(v1, v2);
            }

            out_d0[j*nq0+i] = prod_sum;

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
                    vec_t v1 = in[start + k*nq0+i];
                    vec_t v2 = D1[k*nq1 + j];

                    prod_sum.fma(v1,v2);
                }

                out_d1[ start + j*nq0+i] = prod_sum;
            }
        }

    }

    //Direction 3
    for (int i = 0; i < nq0 * nq1; ++i)
    {
        for (int j = 0; j < nq2; ++j)
        {

            vec_t prod_sum = 0.0;
            for (int k = 0; k < nq2; ++k)
            {
                vec_t v1 = in[k*nq0*nq1 + i];
                vec_t v2 = D2[k*nq2 + j];

                prod_sum.fma(v1,v2);
            }

            out_d2[j*nq0*nq1 + i] = prod_sum;
        }
    }

}

template<int NUMQUAD0, int NUMQUAD1, int NUMQUAD2, bool DEFORMED>
static void PhysDerivHexKernel(
    const std::vector<vec_t, allocator<vec_t>> &inptr,
    const std::vector<vec_t, allocator<vec_t>> &Z0,
    const std::vector<vec_t, allocator<vec_t>> &Z1,
    const std::vector<vec_t, allocator<vec_t>> &Z2,
    const std::vector<vec_t, allocator<vec_t>> &D0,
    const std::vector<vec_t, allocator<vec_t>> &D1,
    const std::vector<vec_t, allocator<vec_t>> &D2,
    const vec_t* df_ptr,
    std::vector<vec_t, allocator<vec_t>> &outptr_d0,
    std::vector<vec_t, allocator<vec_t>> &outptr_d1,
    std::vector<vec_t, allocator<vec_t>> &outptr_d2)
{
    boost::ignore_unused(Z0, Z1, Z2);
    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;
    constexpr auto ndf = 9;

    PhysDerivTensor3DKernel<NUMQUAD0, NUMQUAD1, NUMQUAD2>(
        inptr,
        D0, D1, D2,
        outptr_d0, outptr_d1, outptr_d2);

    vec_t df0, df1, df2, df3, df4, df5, df6, df7, df8;
    if(!DEFORMED){
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

    int cnt_ijk = 0;
    for (int k = 0; k < nq2; ++k)
    {
        for (int j = 0; j < nq1; ++j)
        {
            for (int i = 0; i < nq0; ++i, ++cnt_ijk)
            {
                vec_t d0 = outptr_d0[cnt_ijk]; //Load 1x
                vec_t d1 = outptr_d1[cnt_ijk]; //Load 1x
                vec_t d2 = outptr_d2[cnt_ijk]; //Load 1x

                if(DEFORMED){
                    df0 = df_ptr[cnt_ijk*ndf + 0];
                    df1 = df_ptr[cnt_ijk*ndf + 1];
                    df2 = df_ptr[cnt_ijk*ndf + 2];
                    df3 = df_ptr[cnt_ijk*ndf + 3];
                    df4 = df_ptr[cnt_ijk*ndf + 4];
                    df5 = df_ptr[cnt_ijk*ndf + 5];
                    df6 = df_ptr[cnt_ijk*ndf + 6];
                    df7 = df_ptr[cnt_ijk*ndf + 7];
                    df8 = df_ptr[cnt_ijk*ndf + 8];
                }

                vec_t out0 = d0 * df0;
                out0.fma(d1, df1);
                out0.fma(d2, df2);
                outptr_d0[cnt_ijk] = out0; //Store 1x

                vec_t out1 = d0 * df3;
                out1.fma(d1, df4);
                out1.fma(d2, df5);
                outptr_d1[cnt_ijk] = out1; //Store 1x

                vec_t out2 = d0 * df6;
                out2.fma(d1, df7);
                out2.fma(d2, df8);
                outptr_d2[cnt_ijk] = out2; //Store 1x

            }
        }
    }

}

template<int NUMQUAD0, int NUMQUAD1, int NUMQUAD2, bool DEFORMED>
static void PhysDerivPrismKernel(
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& Z0,
    const std::vector<vec_t, allocator<vec_t>>& Z1,
    const std::vector<vec_t, allocator<vec_t>>& Z2,
    const std::vector<vec_t, allocator<vec_t>>& D0,
    const std::vector<vec_t, allocator<vec_t>>& D1,
    const std::vector<vec_t, allocator<vec_t>>& D2,
    const vec_t* df_ptr,
    std::vector<vec_t, allocator<vec_t>>& out_d0,
    std::vector<vec_t, allocator<vec_t>>& out_d1,
    std::vector<vec_t, allocator<vec_t>>& out_d2)
{
    boost::ignore_unused(Z1);

    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;
    constexpr auto ndf = 9;

    PhysDerivTensor3DKernel<NUMQUAD0, NUMQUAD1, NUMQUAD2>(
        in,
        D0, D1, D2,
        out_d0, out_d1, out_d2);

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

    int cnt_ijk = 0;
    for (int k = 0; k < nq2; ++k)
    {
        vec_t xfrm_eta2 = 2.0 / (1.0 - Z2[k]); //Load 1x
        for (int j = 0; j < nq1; ++j)
        {
            for (int i = 0; i < nq0; ++i, ++cnt_ijk)
            {
                //Chain-rule for eta_0 and eta_2
                vec_t d0 = out_d0[cnt_ijk] * xfrm_eta2; //Load 1x
                vec_t d2 = out_d2[cnt_ijk];//Load 1x
                vec_t xfrm_eta0 = 0.5 * (1.0 + Z0[i]); //Load 1x
                d2.fma(xfrm_eta0, d0);

                vec_t d1 = out_d1[cnt_ijk]; //Load 1x

                if (DEFORMED)
                {
                    df0 = df_ptr[cnt_ijk*ndf];
                    df1 = df_ptr[cnt_ijk*ndf + 1];
                    df2 = df_ptr[cnt_ijk*ndf + 2];
                    df3 = df_ptr[cnt_ijk*ndf + 3];
                    df4 = df_ptr[cnt_ijk*ndf + 4];
                    df5 = df_ptr[cnt_ijk*ndf + 5];
                    df6 = df_ptr[cnt_ijk*ndf + 6];
                    df7 = df_ptr[cnt_ijk*ndf + 7];
                    df8 = df_ptr[cnt_ijk*ndf + 8];
                }

                //Metric for eta_0, xi_1, eta_2
                vec_t out0 = d0 * df0;
                out0.fma(d1, df1);
                out0.fma(d2, df2);
                out_d0[cnt_ijk] = out0; //Store 1x

                vec_t out1 = d0 * df3;
                out1.fma(d1, df4);
                out1.fma(d2, df5);
                out_d1[cnt_ijk] = out1; //Store 1x

                vec_t out2 = d0 * df6;
                out2.fma(d1, df7);
                out2.fma(d2, df8);
                out_d2[cnt_ijk] = out2; //Store 1x
            }
        }
    }
}

template<int NUMQUAD0, int NUMQUAD1, int NUMQUAD2, bool DEFORMED>
static void PhysDerivPyrKernel(
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& Z0,
    const std::vector<vec_t, allocator<vec_t>>& Z1,
    const std::vector<vec_t, allocator<vec_t>>& Z2,
    const std::vector<vec_t, allocator<vec_t>>& D0,
    const std::vector<vec_t, allocator<vec_t>>& D1,
    const std::vector<vec_t, allocator<vec_t>>& D2,
    const vec_t* df_ptr,
    std::vector<vec_t, allocator<vec_t>>& out_d0,
    std::vector<vec_t, allocator<vec_t>>& out_d1,
    std::vector<vec_t, allocator<vec_t>>& out_d2)
{
    boost::ignore_unused(Z1);

    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;
    constexpr auto ndf = 9;

    PhysDerivTensor3DKernel<NUMQUAD0, NUMQUAD1, NUMQUAD2>(
        in,
        D0, D1, D2,
        out_d0, out_d1, out_d2);

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

    int cnt_ijk = 0;
    for (int k = 0; k < nq2; ++k)
    {
        vec_t xfrm_eta2 = 2.0 / (1.0 - Z2[k]); //Load 1x

        for (int j = 0; j < nq1; ++j)
        {
            vec_t xfrm_eta1 = 0.5*(1+Z1[j]); //Load 1x

            for (int i = 0; i < nq0; ++i, ++cnt_ijk)
            {
                
                //Chain-rule for eta_0 and eta_2
                vec_t d0 = out_d0[cnt_ijk] * xfrm_eta2; //Load 1x
                vec_t d1 = out_d1[cnt_ijk] * xfrm_eta2; //Load 1x
                
                vec_t d2 = out_d2[cnt_ijk];//Load 1x
                vec_t xfrm_eta0 = 0.5*(1+Z0[i]); //Load 1x
                d2.fma(xfrm_eta0, d0);
                d2.fma(xfrm_eta1, d1); 
                
                if (DEFORMED)
                {
                    df0 = df_ptr[cnt_ijk*ndf];
                    df1 = df_ptr[cnt_ijk*ndf + 1];
                    df2 = df_ptr[cnt_ijk*ndf + 2];
                    df3 = df_ptr[cnt_ijk*ndf + 3];
                    df4 = df_ptr[cnt_ijk*ndf + 4];
                    df5 = df_ptr[cnt_ijk*ndf + 5];
                    df6 = df_ptr[cnt_ijk*ndf + 6];
                    df7 = df_ptr[cnt_ijk*ndf + 7];
                    df8 = df_ptr[cnt_ijk*ndf + 8];
                }
                
                //Metric for eta_0, xi_1, eta_2
                vec_t out0 = d0 * df0;
                out0.fma(d1, df1);
                out0.fma(d2, df2);
                out_d0[cnt_ijk] = out0; //Store 1x
                
                vec_t out1 = d0 * df3;
                out1.fma(d1, df4);
                out1.fma(d2, df5);
                out_d1[cnt_ijk] = out1; //Store 1x
                
                vec_t out2 = d0 * df6;
                out2.fma(d1, df7);
                out2.fma(d2, df8);
                out_d2[cnt_ijk] = out2; //Store 1x
            }
        }
    }
}

template<int NUMQUAD0, int NUMQUAD1, int NUMQUAD2, bool DEFORMED>
static void PhysDerivTetKernel(
    const std::vector<vec_t, allocator<vec_t>>& in,
    const std::vector<vec_t, allocator<vec_t>>& Z0,
    const std::vector<vec_t, allocator<vec_t>>& Z1,
    const std::vector<vec_t, allocator<vec_t>>& Z2,
    const std::vector<vec_t, allocator<vec_t>>& D0,
    const std::vector<vec_t, allocator<vec_t>>& D1,
    const std::vector<vec_t, allocator<vec_t>>& D2,
    const vec_t* df_ptr,
    std::vector<vec_t, allocator<vec_t>>& diff0, // tmp store
    std::vector<vec_t, allocator<vec_t>>& diff1, // tmp store
    std::vector<vec_t, allocator<vec_t>>& diff2, // tmp store
    std::vector<vec_t, allocator<vec_t>>& out_d0,
    std::vector<vec_t, allocator<vec_t>>& out_d1,
    std::vector<vec_t, allocator<vec_t>>& out_d2)
{

    constexpr auto nq0 = NUMQUAD0, nq1 = NUMQUAD1, nq2 = NUMQUAD2;
    constexpr auto ndf = 9;

    PhysDerivTensor3DKernel<NUMQUAD0, NUMQUAD1, NUMQUAD2>(
        in,
        D0, D1, D2,
        diff0, diff1, diff2);

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

    for (int k = 0, eta0 = 0; k < nq2; ++k)
    {
        vec_t xfrm_eta2 = 2.0 / (1.0 - Z2[k]); //Load 1x

        for (int j = 0; j < nq1; ++j)
        {

            vec_t xfrm_eta1 = 2.0 / (1.0 -  Z1[j]); //Load 1x
            vec_t xfrm = xfrm_eta1 * xfrm_eta2;

            for (int i = 0; i < nq0; ++i, ++eta0)
            {
                vec_t d0 = xfrm * diff0[eta0]; //Load 1x

                out_d0[eta0] = d0; //Store 1x
                diff0[eta0] = d0; //Store 1x partial form for reuse

            }
        }
    }

    for (int k = 0, eta0 = 0; k < nq2; ++k)
    {
        vec_t xfrm_eta2 = 2.0 / (1.0 - Z2[k]); //Load 1x

        for (int j = 0; j < nq1; ++j)
        {
            for (int i = 0; i < nq0; ++i, ++eta0)
            {
                vec_t xfrm_eta0 = 0.5 * (1.0 + Z0[i]); //Load 1x

                vec_t out0 = xfrm_eta0 * diff0[eta0]; //Load 1x
                diff0[eta0] = out0; //2 * (1 + eta_0) / (1 - eta_1)(1-eta2) | store 1x

                vec_t d1 = diff1[eta0]; //Load 1x
                d1 = xfrm_eta2 * d1;

                out_d1[eta0] = out0 + d1; //store 1x
                diff1[eta0] = d1; //store 1x
            }
        }
    }

    for (int k = 0, eta0 = 0; k < nq2; ++k)
    {
        for (int j = 0; j < nq1; ++j)
        {
            vec_t xfrm_eta1 = 0.5 * (1.0 + Z1[j]); //Load 1x

            for (int i = 0; i < nq0; ++i, ++eta0)
            {
                vec_t out = diff0[eta0]; //Load 1x
                vec_t d1 = diff1[eta0];  //Load 1x
                out.fma(d1, xfrm_eta1);
                out = out + diff2[eta0]; //Load 1x
                out_d2[eta0] = out;      //Store 1x
            }
        }
    }

    for (int k = 0, cnt_kji = 0; k < nq2; ++k)
    {
        for (int j = 0; j < nq1; ++j)
        {
            for (int i = 0; i < nq0; ++i, ++cnt_kji)
            {
                vec_t d0 = out_d0[cnt_kji]; //Load 1x
                vec_t d1 = out_d1[cnt_kji]; //Load 1x
                vec_t d2 = out_d2[cnt_kji]; //Load 1x

                if (DEFORMED)
                {
                    df0 = df_ptr[cnt_kji*ndf + 0];
                    df1 = df_ptr[cnt_kji*ndf + 1];
                    df2 = df_ptr[cnt_kji*ndf + 2];
                    df3 = df_ptr[cnt_kji*ndf + 3];
                    df4 = df_ptr[cnt_kji*ndf + 4];
                    df5 = df_ptr[cnt_kji*ndf + 5];
                    df6 = df_ptr[cnt_kji*ndf + 6];
                    df7 = df_ptr[cnt_kji*ndf + 7];
                    df8 = df_ptr[cnt_kji*ndf + 8];
                }

                vec_t out0 = d0 * df0;
                out0.fma(d1, df1);
                out0.fma(d2, df2);
                out_d0[cnt_kji] = out0; //Store 1x

                vec_t out1 = d0 * df3;
                out1.fma(d1, df4);
                out1.fma(d2, df5);
                out_d1[cnt_kji] = out1; //Store 1x

                vec_t out2 = d0 * df6;
                out2.fma(d1, df7);
                out2.fma(d2, df8);
                out_d2[cnt_kji] = out2; //Store 1x
            }
        }
    }

}

} // namespace MatrixFree
} // namespace Nektar

#endif
