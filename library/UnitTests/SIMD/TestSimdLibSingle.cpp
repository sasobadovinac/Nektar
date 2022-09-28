///////////////////////////////////////////////////////////////////////////////
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
/// The above copyright notice and this permission notice shall be included
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

#include <LibUtilities/SimdLib/io.hpp>
#include <LibUtilities/SimdLib/tinysimd.hpp>

#include <boost/core/ignore_unused.hpp>
#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <iostream>

// define types in use and number of lanes
#define NUM_LANES_32BITS 1
#define ALIGNMENT 4
#define USING_SCALAR
#if defined(__x86_64__)
#if defined(__AVX512F__) && defined(NEKTAR_ENABLE_SIMD_AVX512)
#define USING_AVX512
#undef NUM_LANES_32BITS
#define NUM_LANES_32BITS 16
#undef ALIGNMENT
#define ALIGNMENT 64
#undef USING_SCALAR
#elif defined(__AVX2__) && defined(NEKTAR_ENABLE_SIMD_AVX2)
#define USING_AVX2
#undef NUM_LANES_32BITS
#define NUM_LANES_32BITS 8
#undef ALIGNMENT
#define ALIGNMENT 32
#undef USING_SCALAR
#elif defined(__SSE2__) && defined(NEKTAR_ENABLE_SIMD_SSE2)
#define USING_SSE2
#undef NUM_LANES_32BITS
#define NUM_LANES_32BITS 4
#undef ALIGNMENT
#define ALIGNMENT 64
#undef USING_SCALAR
#endif
#endif
#if defined(__ARM_FEATURE_SVE) && defined(NEKTAR_ENABLE_SIMD_SVE)
#define USING_SVE
#undef NUM_LANES_32BITS
#define NUM_LANES_32BITS __ARM_FEATURE_SVE_BITS / 32
#undef USING_SCALAR
#endif

namespace Nektar
{
namespace SimdLibTests
{
using namespace tinysimd;
using vec_t = simd<float>;

BOOST_AUTO_TEST_CASE(SimdLibSingle_width_alignment)
{
    std::size_t width, alignment;

#if defined(USING_SCALAR)
    std::cout << "scalar float" << std::endl;
#endif
#if defined(USING_AVX2)
    std::cout << "avx2 float" << std::endl;
#endif
#if defined(USING_AVX512)
    std::cout << "avx512 float" << std::endl;
#endif

    // float
    width     = simd<float>::width;
    alignment = simd<float>::alignment;
    BOOST_CHECK_EQUAL(width, NUM_LANES_32BITS);
    BOOST_CHECK_EQUAL(alignment, ALIGNMENT);
    // std::int32_t index forcing # of lanes
    width     = simd<std::int32_t, vec_t::width>::width;
    alignment = simd<std::int32_t, vec_t::width>::alignment;
    BOOST_CHECK_EQUAL(width, NUM_LANES_32BITS);
    BOOST_CHECK_EQUAL(alignment, ALIGNMENT);
    // std::int32_t default index
    width     = simd<std::int32_t>::width;
    alignment = simd<std::int32_t>::alignment;
    BOOST_CHECK_EQUAL(width, NUM_LANES_32BITS);
    BOOST_CHECK_EQUAL(alignment, ALIGNMENT);

#if defined(USING_SVE)
    std::cout << "sve float" << std::endl;
    //
    // these are going to be machine/compilation dependent
    // we are forcing VLA -> VLST
    // so we can still check consistency with __ARM_FEATURE_SVE_BITS
    //
    // // int 32 bit on number of lanes for int 64 bit
    // width = simd<int>::width;
    // alignment = simd<int>::alignment;
    // BOOST_CHECK_EQUAL(width, __ARM_FEATURE_SVE_BITS/sizeof(std::int32_t));
    // BOOST_CHECK_EQUAL(alignment, 32);
    // int
    width     = simd<int>::width;
    alignment = simd<int>::alignment;
    BOOST_CHECK_EQUAL(width, NUM_LANES_32BITS);
    BOOST_CHECK_EQUAL(alignment, __ARM_FEATURE_SVE_BITS / sizeof(std::int64_t));
    // float
    width     = simd<float>::width;
    alignment = simd<float>::alignment;
    BOOST_CHECK_EQUAL(width, NUM_LANES_32BITS);
    BOOST_CHECK_EQUAL(alignment, __ARM_FEATURE_SVE_BITS / sizeof(float));
#endif
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_type_traits)
{
    using namespace details;

    BOOST_CHECK_EQUAL(has_width<float>::value, false);
    BOOST_CHECK_EQUAL(has_width<vec_t>::value, true);

    BOOST_CHECK_EQUAL(has_alignment<float>::value, false);
    BOOST_CHECK_EQUAL(has_alignment<vec_t>::value, true);

    BOOST_CHECK_EQUAL(has_scalarType<float>::value, false);
    BOOST_CHECK_EQUAL(has_scalarType<vec_t>::value, true);

    BOOST_CHECK_EQUAL(is_vector<float>::value, false);
    BOOST_CHECK_EQUAL(is_vector<vec_t>::value, true);

    BOOST_CHECK_EQUAL(is_vector_floating_point<float>::value, false);
    BOOST_CHECK_EQUAL(is_vector_floating_point<simd<int>>::value, false);
    BOOST_CHECK_EQUAL(is_vector_floating_point<vec_t>::value, true);
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_mem_size)
{
    BOOST_CHECK_EQUAL(sizeof(vec_t), sizeof(float) * vec_t::width);
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_ctors)
{
    vec_t avec1;

    vec_t::scalarType ascalar = 0;
    vec_t avec2(ascalar);
    vec_t avec3{ascalar};
    vec_t avec4 = ascalar;

    vec_t avec5(avec2);
    vec_t avec6{avec4};

    vec_t avec7(avec2._data);
    vec_t avec8{avec2._data};

    vec_t::vectorType anative;
    vec_t avec9(anative);
    vec_t avec10{anative};

    boost::ignore_unused(avec1, avec3, avec5, avec6, avec7, avec8, avec9,
                         avec10);
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_load)
{
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // double brace to deal with gcc 4.8.5 ...
    vec_t avec;
    avec.load(ascalararr.data());
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_load_implicit)
{
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // double brace to deal with gcc 4.8.5 ...
    vec_t avec;
    avec = *(reinterpret_cast<vec_t *>(ascalararr.data()));
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_load_aligned)
{
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // double brace to deal with gcc 4.8.5 ...
    vec_t avec;
    avec.load(ascalararr.data(), is_aligned);
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_load_unaligned)
{
    std::array<float, vec_t::width> ascalararr{
        {}}; // double brace to deal with gcc 4.8.5 ...
    vec_t avec;
    avec.load(ascalararr.data(), is_not_aligned);
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_store)
{
    float val = 4.0;
    vec_t avec(val);
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // float brace to deal with gcc 4.8.5 ...
    avec.store(ascalararr.data());

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], val);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_store_aligned)
{
    float val = 4.0;
    vec_t avec(val);
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // float brace to deal with gcc 4.8.5 ...
    avec.store(ascalararr.data(), is_aligned);

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], val);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_store_unaligned)
{
    float val = 4.0;
    vec_t avec(val);
    std::array<float, vec_t::width> ascalararr{
        {}}; // float brace to deal with gcc 4.8.5 ...
    avec.store(ascalararr.data(), is_not_aligned);

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], val);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_store_non_temporal)
{
    float val = 4.0;
    vec_t avec(val);
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // float brace to deal with gcc 4.8.5 ...
    avec.store(ascalararr.data(), is_not_reused);

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], val);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibSingle_broadcast)
{
    // explictly float literal for
    // Visual Studio 14 2015 Win64 compiler bug
    vec_t::scalarType ascalar{3.333f};
    vec_t avec;
    avec.broadcast(ascalar);
}

BOOST_AUTO_TEST_CASE(SimdLibSingle_subscript_assign_read)
{
    vec_t avec;
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // double brace to deal with gcc 4.8.5 ...

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        ascalararr[i] = i;
    }

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        avec[i] = ascalararr[i];
    }

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], avec[i]);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_gather32)
{
    vec_t avec;
    using index_t = simd<std::uint32_t>;
    index_t aindexvec;

    // create and fill index
    std::array<std::uint32_t, vec_t::width> aindex;
    aindex[0] = 0;
    if (vec_t::width > 2)
    {
        aindex[1] = 3;
        aindex[2] = 5;
        aindex[3] = 6;
    }
    if (vec_t::width > 4)
    {
        aindex[4] = 8;
        aindex[5] = 15;
        aindex[6] = 20;
        aindex[7] = 23;
    }
    if (vec_t::width > 8)
    {
        aindex[8]  = 24;
        aindex[9]  = 28;
        aindex[10] = 33;
        aindex[11] = 40;
        aindex[12] = 41;
        aindex[13] = 45;
        aindex[14] = 60;
        aindex[15] = 61;
    }

    // load index
    aindexvec.load(aindex.data(), is_not_aligned);

    // create and fill scalar array
    constexpr size_t scalarArraySize = 64;
    std::array<float, scalarArraySize> ascalararr;
    for (size_t i = 0; i < scalarArraySize; ++i)
    {
        ascalararr[i] = i;
    }

    avec.gather(ascalararr.data(), aindexvec);

    // check
    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[aindex[i]], avec[i]);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_scatter32)
{
    vec_t avec;
    using index_t = simd<std::uint32_t>;
    index_t aindexvec;

    // create and fill index
    std::array<std::uint32_t, vec_t::width> aindex;
    aindex[0] = 1;
    if (vec_t::width > 1)
    {
        aindex[1] = 3;
    }
    if (vec_t::width > 2)
    {
        aindex[2] = 5;
        aindex[3] = 6;
    }
    if (vec_t::width > 4)
    {
        aindex[4] = 8;
        aindex[5] = 15;
        aindex[6] = 20;
        aindex[7] = 30;
    }
    if (vec_t::width > 8)
    {
        aindex[8]  = 31;
        aindex[9]  = 32;
        aindex[10] = 35;
        aindex[11] = 40;
        aindex[12] = 41;
        aindex[13] = 45;
        aindex[14] = 60;
        aindex[15] = 61;
    }

    // load index
    aindexvec.load(aindex.data(), is_not_aligned);

    // create scalar array
    constexpr size_t scalarArraySize = 64;
    std::array<float, scalarArraySize> ascalararr;

    // fill vector
    alignas(vec_t::alignment) std::array<float, vec_t::width> avecarr{{}};
    avecarr[0] = 10;
    if (vec_t::width > 1)
    {
        avecarr[1] = 9;
    }

    if (vec_t::width > 2)
    {
        avec[2] = 8;
        avec[3] = 7;
    }
    if (vec_t::width > 4)
    {
        avec[4] = 4;
        avec[5] = 3;
        avec[6] = 2;
        avec[7] = 1;
    }
    avec.load(avecarr.data());

    avec.scatter(ascalararr.data(), aindexvec);

    // check
    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(avec[i], ascalararr[aindex[i]]);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_add_unary)
{
    float val1 = -4.0;
    float val2 = 2.0;
    vec_t res(val1);
    vec_t avec(val2);
    res += avec;
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // float brace to deal with gcc 4.8.5 ...
    res.store(ascalararr.data());

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], val1 + val2);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_sub_unary)
{
    float val1 = -4.0;
    float val2 = 2.0;
    vec_t res(val1);
    vec_t avec(val2);
    res -= avec;
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // float brace to deal with gcc 4.8.5 ...
    res.store(ascalararr.data());

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], val1 - val2);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_mul_unary)
{
    float val1 = -4.0;
    float val2 = 2.0;
    vec_t res(val1);
    vec_t avec(val2);
    res *= avec;
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // float brace to deal with gcc 4.8.5 ...
    res.store(ascalararr.data());

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], val1 * val2);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_div_unary)
{
    float val1 = -4.0;
    float val2 = 2.0;
    vec_t res(val1);
    vec_t avec(val2);
    res /= avec;
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // float brace to deal with gcc 4.8.5 ...
    res.store(ascalararr.data());

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], val1 / val2);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_add_binary)
{
    float val1 = -4.0;
    float val2 = 2.5;
    vec_t avec1(val1);
    vec_t avec2(val2);
    vec_t res = avec1 + avec2;
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // double brace to deal with gcc 4.8.5 ...
    res.store(ascalararr.data());

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], val1 + val2);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_sub_binary)
{
    float val1 = -4.0;
    float val2 = 2.5;
    vec_t avec1(val1);
    vec_t avec2(val2);
    vec_t res = avec1 - avec2;
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // double brace to deal with gcc 4.8.5 ...
    res.store(ascalararr.data());

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], val1 - val2);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_mul_binary)
{
    float val1 = -4.0;
    float val2 = 2.5;
    vec_t avec1(val1);
    vec_t avec2(val2);
    vec_t res = avec1 * avec2;
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // double brace to deal with gcc 4.8.5 ...
    res.store(ascalararr.data());

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], val1 * val2);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_div_binary)
{
    float val1 = -4.0;
    float val2 = 2.5;
    vec_t avec1(val1);
    vec_t avec2(val2);
    vec_t res = avec1 / avec2;
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // double brace to deal with gcc 4.8.5 ...
    res.store(ascalararr.data());

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], val1 / val2);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_add_mul)
{
    float val1 = -4.0;
    float val2 = 2.0;
    float val3 = 2.0;
    vec_t avec1(val1);
    vec_t avec2(val2);
    vec_t avec3(val3);
    vec_t res = avec1 + avec2 * avec3;
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // float brace to deal with gcc 4.8.5 ...
    res.store(ascalararr.data());

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], val1 + val2 + val3);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_fused_add_mul)
{
    float val1 = -4.0;
    float val2 = 1.5;
    float val3 = 5.0;
    vec_t avec1(val1);
    vec_t avec2(val2);
    vec_t avec3(val3);
    avec1.fma(avec2, avec3);
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // double brace to deal with gcc 4.8.5 ...
    avec1.store(ascalararr.data());

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], val1 + val2 * val3);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_sqrt)
{
    float val = 4.0;
    vec_t avec(val);
    vec_t asqrt = sqrt(avec);
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // float brace to deal with gcc 4.8.5 ...
    asqrt.store(ascalararr.data());

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], std::sqrt(val));
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_abs)
{
    float val = -4.0;
    vec_t avec(val);
    vec_t aabs = abs(avec);
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // float brace to deal with gcc 4.8.5 ...
    aabs.store(ascalararr.data());

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], std::abs(val));
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_log)
{
    float val = 4.0;
    vec_t avec(val);
    vec_t alog = log(avec);
    alignas(vec_t::alignment) std::array<float, vec_t::width> ascalararr{
        {}}; // float brace to deal with gcc 4.8.5 ...
    alog.store(ascalararr.data());

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        BOOST_CHECK_EQUAL(ascalararr[i], std::log(val));
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_greater)
{
    float aval = 4.0;
    vec_t avec(aval);
    using mask_t = simd<bool, vec_t::width>;
    mask_t amask;

    amask = avec > avec;
    // check
    alignas(vec_t::alignment) std::array<std::uint32_t, vec_t::width>
        ascalararr{{}}; // float brace to deal with gcc 4.8.5 ...
    amask.store(ascalararr.data());
    for (size_t i = 0; i < vec_t::width; ++i)
    {
        // type conversion make lvalue in rvalue, needed pre-c++17
        BOOST_CHECK_EQUAL(ascalararr[i], (std::uint32_t)mask_t::false_v);
    }

    float bval = 3.0;
    vec_t bvec(bval);

    amask = avec > bvec;
    // check
    amask.store(ascalararr.data());
    for (size_t i = 0; i < vec_t::width; ++i)
    {
        // type conversion make lvalue in rvalue, needed pre-c++17
        BOOST_CHECK_EQUAL(ascalararr[i], (std::uint32_t)mask_t::true_v);
    }

    float cval = 5.0;
    vec_t cvec(cval);

    amask = avec > cvec;
    // check
    amask.store(ascalararr.data());
    for (size_t i = 0; i < vec_t::width; ++i)
    {
        // type conversion make lvalue in rvalue, needed pre-c++17
        BOOST_CHECK_EQUAL(ascalararr[i], (std::uint32_t)mask_t::false_v);
    }

    if (vec_t::width == 4)
    {
        alignas(vec_t::alignment) std::array<float, 4> ascalararr2{
            {1.0, 2.0, 3.0, 4.0}}; // float brace to deal with gcc 4.8.5 ...
        float dval = 2.0;
        vec_t dvec(dval);
        vec_t evec;
        evec.load(ascalararr2.data());

        amask = dvec > evec;
        // check
        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(static_cast<bool>(amask[i]), dvec[i] > evec[i]);
        }
    }

    if (vec_t::width == 8)
    {
        alignas(vec_t::alignment) std::array<float, 8> ascalararr2{
            {1.0, 2.0, 3.0, 4.0, 3.0, 2.0,
             1.0}}; // double brace to deal with gcc 4.8.5 ...
        float dval = 2.0;
        vec_t dvec(dval);
        vec_t evec;
        evec.load(ascalararr2.data());

        amask = dvec > evec;
        // check
        for (size_t i = 0; i < vec_t::width; ++i)
        {
            BOOST_CHECK_EQUAL(static_cast<bool>(amask[i]),
                              dval > ascalararr2[i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_logic_and)
{
    float aval = 4.0;
    vec_t avec(aval);
    using mask_t = simd<bool, vec_t::width>;
    mask_t amask;

    alignas(vec_t::alignment) std::array<std::uint32_t, vec_t::width>
        ascalararr{{}}; // float brace to deal with gcc 4.8.5 ...
    for (size_t i = 0; i < vec_t::width; ++i)
    {
        ascalararr[i] = mask_t::true_v;
    }
    amask.load(ascalararr.data());

    // check
    BOOST_CHECK_EQUAL(amask && false, false);
    BOOST_CHECK_EQUAL(amask && true, true);

    for (size_t i = 0; i < vec_t::width; ++i)
    {
        ascalararr[i] = mask_t::false_v;
    }
    amask.load(ascalararr.data());

    // check
    BOOST_CHECK_EQUAL(amask && false, false);
    BOOST_CHECK_EQUAL(amask && true, false);

    if (vec_t::width > 1)
    {
        ascalararr[0] = mask_t::true_v;
        // check
        BOOST_CHECK_EQUAL(amask && false, false);
        BOOST_CHECK_EQUAL(amask && true, false);
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_load_interleave_unload)
{
    constexpr size_t nDof{5};
    // no padding in load_interleave deinterleave_store
    constexpr size_t nEle{vec_t::width * 5};
    constexpr size_t nDofBlock = nDof * vec_t::width;

    constexpr size_t size{nDof * nEle};
    std::array<float, size> dofScalarArr{{}};
    for (size_t i = 0; i < size; ++i)
    {
        dofScalarArr[i] = i;
    }

    // number of blocks
    size_t nBlock = nEle / vec_t::width;

    // aligned vector
    std::vector<vec_t, allocator<vec_t>> dofVectorArr(nDof);

    float *dataPtr = dofScalarArr.data();
    // loop over blocks vec_t::width elements at the time
    for (size_t b = 0; b < nBlock; ++b)
    {
        // load
        load_interleave(dataPtr, nDof, dofVectorArr);

        // manipulate each block
        for (size_t j = 0; j < nDof; ++j)
        {
            dofVectorArr[j] = dofVectorArr[j] + j;
        }

        // store
        deinterleave_store(dofVectorArr, nDof, dataPtr);
        dataPtr += nDofBlock;
    }

    // check
    for (size_t b = 0, i = 0; b < nBlock; ++b)
    {
        for (size_t j = 0; j < nDof; ++j, ++i)
        {
            BOOST_CHECK_EQUAL(dofScalarArr[i], i + j);
        }
    }
}

BOOST_AUTO_TEST_CASE(SimdLibFloat_io)
{
    vec_t avec(3.14);
    std::cout << avec << std::endl;
}

} // namespace SimdLibTests
} // namespace Nektar
