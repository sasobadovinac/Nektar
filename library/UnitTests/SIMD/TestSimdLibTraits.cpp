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

#include <LibUtilities/SimdLib/tinysimd.hpp>
#include <LibUtilities/SimdLib/io.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/core/ignore_unused.hpp>


#include <array>
#include <cmath>
#include <iostream>

// type in use
#if defined(__SSE2__) && defined(NEKTAR_ENABLE_SIMD_SSE2)
    #define USING_SSE2
#endif
#if defined(__AVX2__) && defined(NEKTAR_ENABLE_SIMD_AVX2)
    #define USING_AVX2
#endif
#if defined(__AVX512__) && defined(NEKTAR_ENABLE_SIMD_AVX512)
    #define USING_AVX512
#endif


namespace Nektar
{
namespace SimdLibTests
{
    using namespace tinysimd;

    BOOST_AUTO_TEST_CASE(SimdLib_width_alignment)
    {
        std::size_t width, alignment;

        #if !defined(USING_SSE2) && !defined(USING_AVX2) && !defined(USING_AVX512)
        std::cout << "scalar" << std::endl;
        // std::int32_t aka (usually) int
        width = simd<std::int32_t>::width;
        alignment = simd<std::int32_t>::alignment;
        BOOST_CHECK_EQUAL(width, 1);
        BOOST_CHECK_EQUAL(alignment, 4);
        // std::int64_t aka (usually) long
        width = simd<std::int64_t>::width;
        alignment = simd<std::int64_t>::alignment;
        BOOST_CHECK_EQUAL(width, 1);
        BOOST_CHECK_EQUAL(alignment, 8);
        // double
        width = simd<double>::width;
        alignment = simd<double>::alignment;
        BOOST_CHECK_EQUAL(width, 1);
        BOOST_CHECK_EQUAL(alignment, 8);
        #endif

        #if defined(USING_SSE2) && !defined(USING_AVX2) && !defined(USING_AVX512)
        std::cout << "sse2" << std::endl;
        // std::int32_t
        width = simd<std::int32_t>::width;
        alignment = simd<std::int32_t>::alignment;
        BOOST_CHECK_EQUAL(width, 4);
        BOOST_CHECK_EQUAL(alignment, 16);
        #endif

        #if defined(USING_AVX2) && !defined(USING_AVX512)
        std::cout << "avx2" << std::endl;
        // int (sse2int4)
        width = simd<int>::width;
        alignment = simd<int>::alignment;
        BOOST_CHECK_EQUAL(width, 4);
        BOOST_CHECK_EQUAL(alignment, 16);
        // long
        width = simd<long>::width;
        alignment = simd<long>::alignment;
        BOOST_CHECK_EQUAL(width, 4);
        BOOST_CHECK_EQUAL(alignment, 32);
        // double
        width = simd<double>::width;
        alignment = simd<double>::alignment;
        BOOST_CHECK_EQUAL(width, 4);
        BOOST_CHECK_EQUAL(alignment, 32);
        // float
        width = simd<float>::width;
        alignment = simd<float>::alignment;
        BOOST_CHECK_EQUAL(width, 8);
        BOOST_CHECK_EQUAL(alignment, 32);
        #endif

        #if defined(USING_AVX512)
        std::cout << "avx512" << std::endl;
        // int (avx2int8)
        width = simd<int>::width;
        alignment = simd<int>::alignment;
        BOOST_CHECK_EQUAL(width, 8);
        BOOST_CHECK_EQUAL(alignment, 32);
        // long
        width = simd<long>::width;
        alignment = simd<long>::alignment;
        BOOST_CHECK_EQUAL(width, 8);
        BOOST_CHECK_EQUAL(alignment, 64);
        // double
        width = simd<double>::width;
        alignment = simd<double>::alignment;
        BOOST_CHECK_EQUAL(width, 8);
        BOOST_CHECK_EQUAL(alignment, 64);
        #endif

    }

}
}
