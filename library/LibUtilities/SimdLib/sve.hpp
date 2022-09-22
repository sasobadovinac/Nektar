///////////////////////////////////////////////////////////////////////////////
//
// File: sve.hpp
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
// Description: Vector type using Armv8 Scalable Vector Extension (SVE).
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITES_SIMDLIB_SVE_H
#define NEKTAR_LIB_LIBUTILITES_SIMDLIB_SVE_H

#if defined(__ARM_FEATURE_SVE)
#include <arm_acle.h>
#include <arm_sve.h>
#endif

#include "allocator.hpp"
#include "traits.hpp"
#include <vector>

namespace tinysimd
{

namespace abi
{

template <typename scalarType> struct sve
{
    using type = void;
};

} // namespace abi

// requires clang >= 12.0.0 or gcc >= 10
// requires -msve-vector-bits=<length>
#if __ARM_FEATURE_SVE_BITS > 0 && defined(NEKTAR_ENABLE_SIMD_SVE)
// from VLA to VLST
// C++ does not allow for incomplete class member types
// to get around that we force a known size at compile time
typedef svfloat64_t svfloat64_vlst_t
    __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));
typedef svint64_t svint64_vlst_t
    __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));
typedef svuint64_t svuint64_vlst_t
    __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));
typedef svbool_t svbool_vlst_t
    __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS)));

// forward declaration of concrete types
template <typename T> struct sveLong;
struct sveDouble;
struct sveMask;

namespace abi
{

// mapping between abstract types and concrete types
template <> struct sve<double>
{
    using type = sveDouble;
};
template <> struct sve<std::int64_t>
{
    using type = sveLong<std::int64_t>;
};
template <> struct sve<std::uint64_t>
{
    using type = sveLong<std::uint64_t>;
};
template <> struct sve<bool>
{
    using type = sveMask;
};

} // namespace abi

// concrete types, could add enable if to allow only unsigned long and long...
template <typename T> struct sveLong
{
    static_assert(std::is_integral<T>::value && sizeof(T) == 8,
                  "8 bytes Integral required.");

    static constexpr unsigned int alignment =
        __ARM_FEATURE_SVE_BITS / sizeof(T);
    static constexpr unsigned int width = alignment / 8;

    using scalarType = T;
    using vectorType =
        typename std::conditional<std::is_signed<T>::value, svint64_vlst_t,
                                  svuint64_vlst_t>::type;
    using scalarArray = scalarType[width];

    // storage
    vectorType _data;

    // ctors
    inline sveLong()                   = default;
    inline sveLong(const sveLong &rhs) = default;
    inline sveLong(const vectorType &rhs) : _data(rhs)
    {
    }
    inline sveLong(const scalarType rhs)
    {
        _data = svdup_s64(rhs);
    }
    explicit inline sveLong(scalarArray &rhs)
    {
        _data = svld1(svptrue_b64(), rhs);
    }

    // store packed
    inline void store(scalarType *p) const
    {
        svst1(svptrue_b64(), p, _data);
    }
    // refer to x86_64 implementations
    // sve has no requirements on alignment
    // nevertheless we should accept valid tags for compatibility
    template <typename TAG,
              typename std::enable_if<is_load_tag<TAG>::value, bool>::type = 0>
    inline void store(scalarType *p, TAG) const
    {
        svst1(svptrue_b64(), p, _data);
    }

    // load packed
    inline void load(const scalarType *p)
    {
        _data = svld1(svptrue_b64(), p);
    }
    // refer to x86_64 implementations
    // sve has no requirements on alignment
    // nevertheless we should accept valid tags for compatibility
    template <typename TAG,
              typename std::enable_if<is_load_tag<TAG>::value, bool>::type = 0>
    inline void load(const scalarType *p, TAG)
    {
        _data = svld1(svptrue_b64(), p);
    }

    // broadcast
    inline void broadcast(const scalarType rhs)
    {
        _data = svdup(rhs);
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
    inline void operator+=(sveLong rhs)
    {
        _data = svadd_x(svptrue_b64(), _data, rhs._data);
    }

    inline void operator-=(sveLong rhs)
    {
        _data = svsub_x(svptrue_b64(), _data, rhs._data);
    }

    inline void operator*=(sveLong rhs)
    {
        _data = svmul_x(svptrue_b64(), _data, rhs._data);
    }

    inline void operator/=(sveLong rhs)
    {
        _data = svdiv_x(svptrue_b64(), _data, rhs._data);
    }
};

template <typename T>
inline sveLong<T> operator+(sveLong<T> lhs, sveLong<T> rhs)
{
    return svadd_x(svptrue_b64(), lhs._data, rhs._data);
}

template <typename T> inline sveLong<T> operator+(sveLong<T> lhs, T rhs)
{
    return svadd_x(svptrue_b64(), lhs._data, sveLong<T>(rhs)._data);
}

template <typename T>
inline sveLong<T> operator-(sveLong<T> lhs, sveLong<T> rhs)
{
    return svsub_x(svptrue_b64(), lhs._data, rhs._data);
}

template <typename T>
inline sveLong<T> operator*(sveLong<T> lhs, sveLong<T> rhs)
{
    return svmul_x(svptrue_b64(), lhs._data, rhs._data);
}

template <typename T>
inline sveLong<T> operator/(sveLong<T> lhs, sveLong<T> rhs)
{
    return svdiv_x(svptrue_b64(), lhs._data, rhs._data);
}

template <typename T> inline sveLong<T> abs(sveLong<T> in)
{
    return svabs_x(svptrue_b64(), in._data);
}

////////////////////////////////////////////////////////////////////////////////

struct sveDouble
{
    static constexpr unsigned int alignment =
        __ARM_FEATURE_SVE_BITS / sizeof(double);
    static constexpr unsigned int width = alignment / 8;

    using scalarType  = double;
    using vectorType  = svfloat64_vlst_t;
    using scalarArray = scalarType[width];

    // storage
    vectorType _data;

    // ctors
    inline sveDouble()                     = default;
    inline sveDouble(const sveDouble &rhs) = default;
    inline sveDouble(const vectorType &rhs) : _data(rhs)
    {
    }
    inline sveDouble(const scalarType rhs)
    {
        _data = svdup_f64(rhs);
    }

    // store packed
    inline void store(scalarType *p) const
    {
        svst1_f64(svptrue_b64(), p, _data);
    }
    // refer to x86_64 implementations
    // sve has no requirements on alignment
    // nevertheless we should accept valid tags for compatibility
    template <typename T,
              typename std::enable_if<is_load_tag<T>::value, bool>::type = 0>
    inline void store(scalarType *p, T) const
    {
        svst1_f64(svptrue_b64(), p, _data);
    }

    // load packed
    inline void load(const scalarType *p)
    {
        _data = svld1_f64(svptrue_b64(), p);
    }
    // refer to x86_64 implementations
    // sve has no requirements on alignment
    // nevertheless we should accept valid tags for compatibility
    template <typename T,
              typename std::enable_if<is_load_tag<T>::value, bool>::type = 0>
    inline void load(const scalarType *p, T)
    {
        _data = svld1_f64(svptrue_b64(), p);
    }

    // broadcast
    inline void broadcast(const scalarType rhs)
    {
        _data = svdup_f64(rhs);
    }

    // gather/scatter
    template <typename T>
    inline void gather(scalarType const *p, const sveLong<T> &indices)
    {
        _data = svld1_gather_index(svptrue_b64(), p, indices._data);
    }

    template <typename T>
    inline void scatter(scalarType *out, const sveLong<T> &indices) const
    {
        svst1_scatter_index(svptrue_b64(), out, indices._data, _data);
    }

    // fma
    // this = this + a * b
    inline void fma(const sveDouble &a, const sveDouble &b)
    {
        _data = svmad_x(svptrue_b64(), a._data, b._data, _data);
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
    inline void operator+=(sveDouble rhs)
    {
        _data = svadd_x(svptrue_b64(), _data, rhs._data);
    }

    inline void operator-=(sveDouble rhs)
    {
        _data = svsub_x(svptrue_b64(), _data, rhs._data);
    }

    inline void operator*=(sveDouble rhs)
    {
        _data = svmul_x(svptrue_b64(), _data, rhs._data);
    }

    inline void operator/=(sveDouble rhs)
    {
        _data = svdiv_x(svptrue_b64(), _data, rhs._data);
    }
};

inline sveDouble operator+(sveDouble lhs, sveDouble rhs)
{
    return svadd_x(svptrue_b64(), lhs._data, rhs._data);
}

inline sveDouble operator-(sveDouble lhs, sveDouble rhs)
{
    return svsub_x(svptrue_b64(), lhs._data, rhs._data);
}

inline sveDouble operator*(sveDouble lhs, sveDouble rhs)
{
    return svmul_x(svptrue_b64(), lhs._data, rhs._data);
}

inline sveDouble operator/(sveDouble lhs, sveDouble rhs)
{
    return svdiv_x(svptrue_b64(), lhs._data, rhs._data);
}

inline sveDouble sqrt(sveDouble in)
{
    return svsqrt_x(svptrue_b64(), in._data);
}

inline sveDouble abs(sveDouble in)
{
    return svabs_x(svptrue_b64(), in._data);
}

inline sveDouble log(sveDouble in)
{
    // there is no sve log intrinsic
    // this is a dreadful implementation and is simply a stop gap measure
    alignas(sveDouble::alignment) sveDouble::scalarArray tmp;
    in.store(tmp);
    for (size_t i = 0; i < sveDouble::width; ++i)
    {
        tmp[i] = std::log(tmp[i]);
    }
    sveDouble ret;
    ret.load(tmp);
    return ret;
}

inline void load_interleave(const double *in, size_t dataLen,
                            std::vector<sveDouble, allocator<sveDouble>> &out)
{

    alignas(sveDouble::alignment) size_t tmp[sveDouble::width] = {};

    // populate scalar index of unknown size
    // (known at compile time)
    for (size_t i = 0; i < sveDouble::width; ++i)
    {
        tmp[i] = i * dataLen;
    }

    using index_t = sveLong<size_t>;
    index_t index0(tmp);
    index_t index1 = index0 + 1ul;

    // 2x unrolled loop -- minimun width is 2
    size_t nBlocks = dataLen / 2;
    for (size_t i = 0; i < nBlocks; ++i)
    {
        out[2 * i + 0].gather(in, index0);
        out[2 * i + 1].gather(in, index1);
        index0 = index0 + 2ul;
        index1 = index1 + 2ul;
    }

    // spillover loop
    for (size_t i = 2 * nBlocks; i < dataLen; ++i)
    {
        out[i].gather(in, index0);
        index0 = index0 + 1ul;
    }
}

inline void deinterleave_store(
    const std::vector<sveDouble, allocator<sveDouble>> &in, size_t dataLen,
    double *out)
{
    alignas(sveDouble::alignment) size_t tmp[sveDouble::width] = {};

    // populate scalar index of unknown size
    // (known at compile time)
    for (size_t i = 0; i < sveDouble::width; ++i)
    {
        tmp[i] = i * dataLen;
    }

    using index_t = sveLong<size_t>;
    index_t index0(tmp);

    for (size_t i = 0; i < dataLen; ++i)
    {
        in[i].scatter(out, index0);
        index0 = index0 + 1ul;
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
struct sveMask : sveLong<std::uint64_t>
{
    // bring in ctors
    using sveLong::sveLong;

    static constexpr scalarType true_v  = -1;
    static constexpr scalarType false_v = 0;
};

inline sveMask operator>(sveDouble lhs, sveDouble rhs)
{
    // set mask
    svbool_vlst_t mask = svcmpgt(svptrue_b64(), lhs._data, rhs._data);
    // abuse set inactive to zero to convert
    sveMask::vectorType sveTrue_v = svdup_u64(sveMask::true_v);
    return svand_z(mask, sveTrue_v, sveTrue_v);
}

// logical and
inline bool operator&&(sveMask lhs, bool rhs)
{
    // set mask
    sveMask::vectorType sveFalse_v = svdup_u64(sveMask::false_v);
    svbool_vlst_t mask = svcmpne(svptrue_b64(), lhs._data, sveFalse_v);
    // is any equal to false (zero)?
    bool tmp = svptest_any(svptrue_b64(), mask);
    return tmp && rhs;
}

#endif // defined(__ARM_FEATURE_SVE_BITS)

} // namespace tinysimd
#endif
