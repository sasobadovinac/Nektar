///////////////////////////////////////////////////////////////////////////////
//
// File: avx2.hpp
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
// Description: Vector type using avx2 extension.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITES_SIMDLIB_AVX2_H
#define NEKTAR_LIB_LIBUTILITES_SIMDLIB_AVX2_H

#if defined(__x86_64__)
#include <immintrin.h>
#if defined(__INTEL_COMPILER) && !defined(TINYSIMD_HAS_SVML)
#define TINYSIMD_HAS_SVML
#endif
#endif
#include "allocator.hpp"
#include "sse2.hpp"
#include "traits.hpp"
#include <cmath>
#include <vector>

namespace tinysimd
{

namespace abi
{

template <typename scalarType, int width = 0> struct avx2
{
    using type = void;
};

} // namespace abi

#if defined(__AVX2__) && defined(NEKTAR_ENABLE_SIMD_AVX2)

// forward declaration of concrete types
template <typename T> struct avx2Int8;
template <typename T> struct avx2Long4;
struct avx2Double4;
struct avx2Float8;
struct avx2Mask4;
struct avx2Mask8;

namespace abi
{

// mapping between abstract types and concrete types
template <> struct avx2<double>
{
    using type = avx2Double4;
};
template <> struct avx2<float>
{
    using type = avx2Float8;
};
template <> struct avx2<std::int64_t>
{
    using type = avx2Long4<std::int64_t>;
};
template <> struct avx2<std::uint64_t>
{
    using type = avx2Long4<std::uint64_t>;
};
template <> struct avx2<std::int32_t>
{
    using type = avx2Int8<std::int32_t>;
};
template <> struct avx2<std::uint32_t>
{
    using type = avx2Int8<std::uint32_t>;
};
template <> struct avx2<bool, 4>
{
    using type = avx2Mask4;
};
template <> struct avx2<bool, 8>
{
    using type = avx2Mask8;
};

} // namespace abi

// concrete types
template <typename T> struct avx2Int8
{
    static_assert(std::is_integral<T>::value && sizeof(T) == 4,
                  "4 bytes Integral required.");

    static constexpr unsigned int width     = 8;
    static constexpr unsigned int alignment = 32;

    using scalarType  = T;
    using vectorType  = __m256i;
    using scalarArray = scalarType[width];

    // storage
    vectorType _data;

    // ctors
    inline avx2Int8()                    = default;
    inline avx2Int8(const avx2Int8 &rhs) = default;
    inline avx2Int8(const vectorType &rhs) : _data(rhs)
    {
    }
    inline avx2Int8(const scalarType rhs)
    {
        _data = _mm256_set1_epi32(rhs);
    }
    explicit inline avx2Int8(scalarArray &rhs)
    {
        _data = _mm256_load_si256(reinterpret_cast<vectorType *>(rhs));
    }

    // store
    inline void store(scalarType *p) const
    {
        _mm256_store_si256(reinterpret_cast<vectorType *>(p), _data);
    }

    template <class flag,
              typename std::enable_if<is_requiring_alignment<flag>::value &&
                                          !is_streaming<flag>::value,
                                      bool>::type = 0>
    inline void store(scalarType *p, flag) const
    {
        _mm256_store_si256(reinterpret_cast<vectorType *>(p), _data);
    }

    template <class flag,
              typename std::enable_if<!is_requiring_alignment<flag>::value,
                                      bool>::type = 0>
    inline void store(scalarType *p, flag) const
    {
        _mm256_storeu_si256(reinterpret_cast<vectorType *>(p), _data);
    }

    inline void load(const scalarType *p)
    {
        _data = _mm256_load_si256(reinterpret_cast<const vectorType *>(p));
    }

    template <class flag,
              typename std::enable_if<is_requiring_alignment<flag>::value &&
                                          !is_streaming<flag>::value,
                                      bool>::type = 0>
    inline void load(const scalarType *p, flag)
    {
        _data = _mm256_load_si256(reinterpret_cast<const vectorType *>(p));
    }

    template <class flag,
              typename std::enable_if<!is_requiring_alignment<flag>::value,
                                      bool>::type = 0>
    inline void load(const scalarType *p, flag)
    {
        _data = _mm256_loadu_si256(reinterpret_cast<const vectorType *>(p));
    }

    inline void broadcast(const scalarType rhs)
    {
        _data = _mm256_set1_epi32(rhs);
    }

    // subscript
    // subscriptsoperators are convienient but expensive
    // should not be used in optimized kernels
    inline scalarType operator[](size_t i) const
    {
        alignas(alignment) scalarArray tmp;
        store(tmp, is_aligned);
        return tmp[i];
    }
};

template <typename T>
inline avx2Int8<T> operator+(avx2Int8<T> lhs, avx2Int8<T> rhs)
{
    return _mm256_add_epi32(lhs._data, rhs._data);
}

template <
    typename T, typename U,
    typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
inline avx2Int8<T> operator+(avx2Int8<T> lhs, U rhs)
{
    return _mm256_add_epi32(lhs._data, _mm256_set1_epi32(rhs));
}

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct avx2Long4
{
    static_assert(std::is_integral<T>::value && sizeof(T) == 8,
                  "8 bytes Integral required.");

    static constexpr unsigned int width     = 4;
    static constexpr unsigned int alignment = 32;

    using scalarType  = T;
    using vectorType  = __m256i;
    using scalarArray = scalarType[width];

    // storage
    vectorType _data;

    // ctors
    inline avx2Long4()                     = default;
    inline avx2Long4(const avx2Long4 &rhs) = default;
    inline avx2Long4(const vectorType &rhs) : _data(rhs)
    {
    }
    inline avx2Long4(const scalarType rhs)
    {
        _data = _mm256_set1_epi64x(rhs);
    }
    explicit inline avx2Long4(scalarArray &rhs)
    {
        _data = _mm256_load_si256(reinterpret_cast<vectorType *>(rhs));
    }

    // store
    inline void store(scalarType *p) const
    {
        _mm256_store_si256(reinterpret_cast<vectorType *>(p), _data);
    }

    template <class flag,
              typename std::enable_if<is_requiring_alignment<flag>::value &&
                                          !is_streaming<flag>::value,
                                      bool>::type = 0>
    inline void store(scalarType *p, flag) const
    {
        _mm256_store_si256(reinterpret_cast<vectorType *>(p), _data);
    }

    template <class flag,
              typename std::enable_if<!is_requiring_alignment<flag>::value,
                                      bool>::type = 0>
    inline void store(scalarType *p, flag) const
    {
        _mm256_storeu_si256(reinterpret_cast<vectorType *>(p), _data);
    }

    inline void load(const scalarType *p)
    {
        _data = _mm256_load_si256(reinterpret_cast<const vectorType *>(p));
    }

    template <class flag,
              typename std::enable_if<is_requiring_alignment<flag>::value &&
                                          !is_streaming<flag>::value,
                                      bool>::type = 0>
    inline void load(const scalarType *p, flag)
    {
        _data = _mm256_load_si256(reinterpret_cast<const vectorType *>(p));
    }

    template <class flag,
              typename std::enable_if<!is_requiring_alignment<flag>::value,
                                      bool>::type = 0>
    inline void load(const scalarType *p, flag)
    {
        _data = _mm256_loadu_si256(reinterpret_cast<const vectorType *>(p));
    }

    inline void broadcast(const scalarType rhs)
    {
        _data = _mm256_set1_epi64x(rhs);
    }

    // subscript
    // subscript operators are convienient but expensive
    // should not be used in optimized kernels
    inline scalarType operator[](size_t i) const
    {
        alignas(alignment) scalarArray tmp;
        store(tmp, is_aligned);
        return tmp[i];
    }
};

template <typename T>
inline avx2Long4<T> operator+(avx2Long4<T> lhs, avx2Long4<T> rhs)
{
    return _mm256_add_epi64(lhs._data, rhs._data);
}

template <
    typename T, typename U,
    typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
inline avx2Long4<T> operator+(avx2Long4<T> lhs, U rhs)
{
    return _mm256_add_epi64(lhs._data, _mm256_set1_epi64x(rhs));
}

////////////////////////////////////////////////////////////////////////////////

struct avx2Double4
{
    static constexpr unsigned width     = 4;
    static constexpr unsigned alignment = 32;

    using scalarType      = double;
    using scalarIndexType = std::uint64_t;
    using vectorType      = __m256d;
    using scalarArray     = scalarType[width];

    // storage
    vectorType _data;

    // ctors
    inline avx2Double4()                       = default;
    inline avx2Double4(const avx2Double4 &rhs) = default;
    inline avx2Double4(const vectorType &rhs) : _data(rhs)
    {
    }
    inline avx2Double4(const scalarType rhs)
    {
        _data = _mm256_set1_pd(rhs);
    }

    // store
    inline void store(scalarType *p) const
    {
        _mm256_store_pd(p, _data);
    }

    template <class flag,
              typename std::enable_if<is_requiring_alignment<flag>::value &&
                                          !is_streaming<flag>::value,
                                      bool>::type = 0>
    inline void store(scalarType *p, flag) const
    {
        _mm256_store_pd(p, _data);
    }

    template <class flag,
              typename std::enable_if<!is_requiring_alignment<flag>::value,
                                      bool>::type = 0>
    inline void store(scalarType *p, flag) const
    {
        _mm256_storeu_pd(p, _data);
    }

    template <class flag, typename std::enable_if<is_streaming<flag>::value,
                                                  bool>::type = 0>
    inline void store(scalarType *p, flag) const
    {
        _mm256_stream_pd(p, _data);
    }

    // load packed
    inline void load(const scalarType *p)
    {
        _data = _mm256_load_pd(p);
    }

    template <class flag,
              typename std::enable_if<is_requiring_alignment<flag>::value,
                                      bool>::type = 0>
    inline void load(const scalarType *p, flag)
    {
        _data = _mm256_load_pd(p);
    }

    template <class flag,
              typename std::enable_if<!is_requiring_alignment<flag>::value,
                                      bool>::type = 0>
    inline void load(const scalarType *p, flag)
    {
        _data = _mm256_loadu_pd(p);
    }

    // broadcast
    inline void broadcast(const scalarType rhs)
    {
        _data = _mm256_set1_pd(rhs);
    }

#if defined(__SSE2__) && defined(NEKTAR_ENABLE_SIMD_SSE2)
    // gather/scatter with sse2
    template <typename T>
    inline void gather(scalarType const *p, const sse2Int4<T> &indices)
    {
        _data = _mm256_i32gather_pd(p, indices._data, 8);
    }

    template <typename T>
    inline void scatter(scalarType *out, const sse2Int4<T> &indices) const
    {
        // no scatter intrinsics for AVX2
        alignas(alignment) scalarArray tmp;
        _mm256_store_pd(tmp, _data);

        out[_mm_extract_epi32(indices._data, 0)] = tmp[0]; // SSE4.1
        out[_mm_extract_epi32(indices._data, 1)] = tmp[1];
        out[_mm_extract_epi32(indices._data, 2)] = tmp[2];
        out[_mm_extract_epi32(indices._data, 3)] = tmp[3];
    }
#endif

    // gather scatter with avx2
    template <typename T>
    inline void gather(scalarType const *p, const avx2Long4<T> &indices)
    {
        _data = _mm256_i64gather_pd(p, indices._data, 8);
    }

    template <typename T>
    inline void scatter(scalarType *out, const avx2Long4<T> &indices) const
    {
        // no scatter intrinsics for AVX2
        alignas(alignment) scalarArray tmp;
        _mm256_store_pd(tmp, _data);

        out[_mm256_extract_epi64(indices._data, 0)] = tmp[0];
        out[_mm256_extract_epi64(indices._data, 1)] = tmp[1];
        out[_mm256_extract_epi64(indices._data, 2)] = tmp[2];
        out[_mm256_extract_epi64(indices._data, 3)] = tmp[3];
    }

    // fma
    // this = this + a * b
    inline void fma(const avx2Double4 &a, const avx2Double4 &b)
    {
        _data = _mm256_fmadd_pd(a._data, b._data, _data);
    }

    // subscript
    // subscript operators are convienient but expensive
    // should not be used in optimized kernels
    inline scalarType operator[](size_t i) const
    {
        alignas(alignment) scalarArray tmp;
        store(tmp, is_aligned);
        return tmp[i];
    }

    // unary ops
    inline void operator+=(avx2Double4 rhs)
    {
        _data = _mm256_add_pd(_data, rhs._data);
    }

    inline void operator-=(avx2Double4 rhs)
    {
        _data = _mm256_sub_pd(_data, rhs._data);
    }

    inline void operator*=(avx2Double4 rhs)
    {
        _data = _mm256_mul_pd(_data, rhs._data);
    }

    inline void operator/=(avx2Double4 rhs)
    {
        _data = _mm256_div_pd(_data, rhs._data);
    }
};

inline avx2Double4 operator+(avx2Double4 lhs, avx2Double4 rhs)
{
    return _mm256_add_pd(lhs._data, rhs._data);
}

inline avx2Double4 operator-(avx2Double4 lhs, avx2Double4 rhs)
{
    return _mm256_sub_pd(lhs._data, rhs._data);
}

inline avx2Double4 operator*(avx2Double4 lhs, avx2Double4 rhs)
{
    return _mm256_mul_pd(lhs._data, rhs._data);
}

inline avx2Double4 operator/(avx2Double4 lhs, avx2Double4 rhs)
{
    return _mm256_div_pd(lhs._data, rhs._data);
}

inline avx2Double4 sqrt(avx2Double4 in)
{
    return _mm256_sqrt_pd(in._data);
}

inline avx2Double4 abs(avx2Double4 in)
{
    // there is no avx2 _mm256_abs_pd intrinsic
    static const __m256d sign_mask = _mm256_set1_pd(-0.); // -0. = 1 << 63
    return _mm256_andnot_pd(sign_mask, in._data);         // !sign_mask & x
}

inline avx2Double4 log(avx2Double4 in)
{
#if defined(TINYSIMD_HAS_SVML)
    return _mm256_log_pd(in._data);
#else
    // there is no avx2 log intrinsic
    // this is a dreadful implementation and is simply a stop gap measure
    alignas(avx2Double4::alignment) avx2Double4::scalarArray tmp;
    in.store(tmp);
    tmp[0] = std::log(tmp[0]);
    tmp[1] = std::log(tmp[1]);
    tmp[2] = std::log(tmp[2]);
    tmp[3] = std::log(tmp[3]);
    avx2Double4 ret;
    ret.load(tmp);
    return ret;
#endif
}

inline void load_interleave(
    const double *in, size_t dataLen,
    std::vector<avx2Double4, allocator<avx2Double4>> &out)
{
    alignas(avx2Double4::alignment)
        size_t tmp[avx2Double4::width] = {0, dataLen, 2 * dataLen, 3 * dataLen};
    using index_t                      = avx2Long4<size_t>;
    index_t index0(tmp);
    index_t index1 = index0 + 1;
    index_t index2 = index0 + 2;
    index_t index3 = index0 + 3;

    // 4x unrolled loop
    constexpr uint16_t unrl = 4;
    size_t nBlocks          = dataLen / unrl;
    for (size_t i = 0; i < nBlocks; ++i)
    {
        out[unrl * i + 0].gather(in, index0);
        out[unrl * i + 1].gather(in, index1);
        out[unrl * i + 2].gather(in, index2);
        out[unrl * i + 3].gather(in, index3);
        index0 = index0 + unrl;
        index1 = index1 + unrl;
        index2 = index2 + unrl;
        index3 = index3 + unrl;
    }

    // spillover loop
    for (size_t i = unrl * nBlocks; i < dataLen; ++i)
    {
        out[i].gather(in, index0);
        index0 = index0 + 1;
    }
}

inline void deinterleave_store(
    const std::vector<avx2Double4, allocator<avx2Double4>> &in, size_t dataLen,
    double *out)
{
    alignas(avx2Double4::alignment)
        size_t tmp[avx2Double4::width] = {0, dataLen, 2 * dataLen, 3 * dataLen};
    using index_t                      = avx2Long4<size_t>;
    index_t index0(tmp);

    for (size_t i = 0; i < dataLen; ++i)
    {
        in[i].scatter(out, index0);
        index0 = index0 + 1;
    }
}

//////////////////////////////////////////////////////////////////////////////

struct avx2Float8
{
    static constexpr unsigned width     = 8;
    static constexpr unsigned alignment = 32;

    using scalarType      = float;
    using scalarIndexType = std::uint32_t;
    using vectorType      = __m256;
    using scalarArray     = scalarType[width];

    // storage
    vectorType _data;

    // ctors
    inline avx2Float8()                      = default;
    inline avx2Float8(const avx2Float8 &rhs) = default;
    inline avx2Float8(const vectorType &rhs) : _data(rhs)
    {
    }
    inline avx2Float8(const scalarType rhs)
    {
        _data = _mm256_set1_ps(rhs);
    }

    // store
    inline void store(scalarType *p) const
    {
        _mm256_store_ps(p, _data);
    }

    template <class flag,
              typename std::enable_if<is_requiring_alignment<flag>::value &&
                                          !is_streaming<flag>::value,
                                      bool>::type = 0>
    inline void store(scalarType *p, flag) const
    {
        _mm256_store_ps(p, _data);
    }

    template <class flag,
              typename std::enable_if<!is_requiring_alignment<flag>::value,
                                      bool>::type = 0>
    inline void store(scalarType *p, flag) const
    {
        _mm256_storeu_ps(p, _data);
    }

    template <class flag, typename std::enable_if<is_streaming<flag>::value,
                                                  bool>::type = 0>
    inline void store(scalarType *p, flag) const
    {
        _mm256_stream_ps(p, _data);
    }

    // load packed
    inline void load(const scalarType *p)
    {
        _data = _mm256_load_ps(p);
    }

    template <class flag,
              typename std::enable_if<is_requiring_alignment<flag>::value,
                                      bool>::type = 0>
    inline void load(const scalarType *p, flag)
    {
        _data = _mm256_load_ps(p);
    }

    template <class flag,
              typename std::enable_if<!is_requiring_alignment<flag>::value,
                                      bool>::type = 0>
    inline void load(const scalarType *p, flag)
    {
        _data = _mm256_loadu_ps(p);
    }

    // broadcast
    inline void broadcast(const scalarType rhs)
    {
        _data = _mm256_set1_ps(rhs);
    }

    // gather scatter with avx2
    template <typename T>
    inline void gather(scalarType const *p, const avx2Int8<T> &indices)
    {
        _data = _mm256_i32gather_ps(p, indices._data, 4);
    }

    template <typename T>
    inline void scatter(scalarType *out, const avx2Int8<T> &indices) const
    {
        // no scatter intrinsics for AVX2
        alignas(alignment) scalarArray tmp;
        _mm256_store_ps(tmp, _data);

        out[_mm256_extract_epi32(indices._data, 0)] = tmp[0];
        out[_mm256_extract_epi32(indices._data, 1)] = tmp[1];
        out[_mm256_extract_epi32(indices._data, 2)] = tmp[2];
        out[_mm256_extract_epi32(indices._data, 3)] = tmp[3];
        out[_mm256_extract_epi32(indices._data, 4)] = tmp[4];
        out[_mm256_extract_epi32(indices._data, 5)] = tmp[5];
        out[_mm256_extract_epi32(indices._data, 6)] = tmp[6];
        out[_mm256_extract_epi32(indices._data, 7)] = tmp[7];
    }

    // fma
    // this = this + a * b
    inline void fma(const avx2Float8 &a, const avx2Float8 &b)
    {
        _data = _mm256_fmadd_ps(a._data, b._data, _data);
    }

    // subscript
    // subscript operators are convienient but expensive
    // should not be used in optimized kernels
    inline scalarType operator[](size_t i) const
    {
        alignas(alignment) scalarArray tmp;
        store(tmp, is_aligned);
        return tmp[i];
    }

    inline scalarType &operator[](size_t i)
    {
        scalarType *tmp = reinterpret_cast<scalarType *>(&_data);
        return tmp[i];
    }

    // unary ops
    inline void operator+=(avx2Float8 rhs)
    {
        _data = _mm256_add_ps(_data, rhs._data);
    }

    inline void operator-=(avx2Float8 rhs)
    {
        _data = _mm256_sub_ps(_data, rhs._data);
    }

    inline void operator*=(avx2Float8 rhs)
    {
        _data = _mm256_mul_ps(_data, rhs._data);
    }

    inline void operator/=(avx2Float8 rhs)
    {
        _data = _mm256_div_ps(_data, rhs._data);
    }
};

inline avx2Float8 operator+(avx2Float8 lhs, avx2Float8 rhs)
{
    return _mm256_add_ps(lhs._data, rhs._data);
}

inline avx2Float8 operator-(avx2Float8 lhs, avx2Float8 rhs)
{
    return _mm256_sub_ps(lhs._data, rhs._data);
}

inline avx2Float8 operator*(avx2Float8 lhs, avx2Float8 rhs)
{
    return _mm256_mul_ps(lhs._data, rhs._data);
}

inline avx2Float8 operator/(avx2Float8 lhs, avx2Float8 rhs)
{
    return _mm256_div_ps(lhs._data, rhs._data);
}

inline avx2Float8 sqrt(avx2Float8 in)
{
    return _mm256_sqrt_ps(in._data);
}

inline avx2Float8 abs(avx2Float8 in)
{
    // there is no avx2 _mm256_abs_ps intrinsic
    static const __m256 sign_mask = _mm256_set1_ps(-0.); // -0. = 1 << 63
    return _mm256_andnot_ps(sign_mask, in._data);        // !sign_mask & x
}

inline avx2Float8 log(avx2Float8 in)
{
    // there is no avx2 log intrinsic
    // this is a dreadful implementation and is simply a stop gap measure
    alignas(avx2Float8::alignment) avx2Float8::scalarArray tmp;
    in.store(tmp);
    tmp[0] = std::log(tmp[0]);
    tmp[1] = std::log(tmp[1]);
    tmp[2] = std::log(tmp[2]);
    tmp[3] = std::log(tmp[3]);
    tmp[4] = std::log(tmp[4]);
    tmp[5] = std::log(tmp[5]);
    tmp[6] = std::log(tmp[6]);
    tmp[7] = std::log(tmp[7]);
    avx2Float8 ret;
    ret.load(tmp);
    return ret;
}

inline void load_interleave(const float *in, std::uint32_t dataLen,
                            std::vector<avx2Float8, allocator<avx2Float8>> &out)
{

    alignas(avx2Float8::alignment) avx2Float8::scalarIndexType tmp[8] = {
        0,           dataLen,     2 * dataLen, 3 * dataLen,
        4 * dataLen, 5 * dataLen, 6 * dataLen, 7 * dataLen};

    using index_t = avx2Int8<avx2Float8::scalarIndexType>;
    index_t index0(tmp);
    index_t index1 = index0 + 1;
    index_t index2 = index0 + 2;
    index_t index3 = index0 + 3;

    // 4x unrolled loop
    size_t nBlocks = dataLen / 4;
    for (size_t i = 0; i < nBlocks; ++i)
    {
        out[4 * i + 0].gather(in, index0);
        out[4 * i + 1].gather(in, index1);
        out[4 * i + 2].gather(in, index2);
        out[4 * i + 3].gather(in, index3);
        index0 = index0 + 4;
        index1 = index1 + 4;
        index2 = index2 + 4;
        index3 = index3 + 4;
    }

    // spillover loop
    for (size_t i = 4 * nBlocks; i < dataLen; ++i)
    {
        out[i].gather(in, index0);
        index0 = index0 + 1;
    }
}

inline void deinterleave_store(
    const std::vector<avx2Float8, allocator<avx2Float8>> &in,
    std::uint32_t dataLen, float *out)
{
    alignas(avx2Float8::alignment) avx2Float8::scalarIndexType tmp[8] = {
        0,           dataLen,     2 * dataLen, 3 * dataLen,
        4 * dataLen, 5 * dataLen, 6 * dataLen, 7 * dataLen};
    using index_t = avx2Int8<avx2Float8::scalarIndexType>;
    index_t index0(tmp);

    for (size_t i = 0; i < dataLen; ++i)
    {
        in[i].scatter(out, index0);
        index0 = index0 + 1;
    }
}

////////////////////////////////////////////////////////////////////////////////

// mask type
// mask is a int type with special properties (broad boolean vector)
// broad boolean vectors defined and allowed values are:
// false=0x0 and true=0xFFFFFFFF
//
// VERY LIMITED SUPPORT...just enough to make cubic eos work...
//
struct avx2Mask4 : avx2Long4<std::uint64_t>
{
    // bring in ctors
    using avx2Long4::avx2Long4;

    static constexpr scalarType true_v  = -1;
    static constexpr scalarType false_v = 0;
};

inline avx2Mask4 operator>(avx2Double4 lhs, avx2Double4 rhs)
{
    return reinterpret_cast<__m256i>(
        _mm256_cmp_pd(lhs._data, rhs._data, _CMP_GT_OQ));
}

inline bool operator&&(avx2Mask4 lhs, bool rhs)
{
    bool tmp =
        _mm256_testc_si256(lhs._data, _mm256_set1_epi64x(avx2Mask4::true_v));

    return tmp && rhs;
}

struct avx2Mask8 : avx2Int8<std::uint32_t>
{
    // bring in ctors
    using avx2Int8::avx2Int8;

    static constexpr scalarType true_v  = -1;
    static constexpr scalarType false_v = 0;
};

inline avx2Mask8 operator>(avx2Float8 lhs, avx2Float8 rhs)
{
    return reinterpret_cast<__m256i>(_mm256_cmp_ps(rhs._data, lhs._data, 1));
}

inline bool operator&&(avx2Mask8 lhs, bool rhs)
{
    bool tmp =
        _mm256_testc_si256(lhs._data, _mm256_set1_epi64x(avx2Mask8::true_v));

    return tmp && rhs;
}

#endif // defined(__AVX2__)

} // namespace tinysimd
#endif
