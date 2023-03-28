///////////////////////////////////////////////////////////////////////////////
//
// File: TestPolylib.cpp
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
// Description: Unit tests for the Polylib library.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Polylib/Polylib.h>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <tuple>
#include <vector>

namespace Nektar
{
namespace UnitTests
{

/// Lower range of number of quadrature points to use for unit tests.
#define NPLOWER 5
/// Upper range of number of quadrature points to use for unit tests.
#define NPUPPER 15
/// Tolerance to be used for floating point comparisons.
#define EPS 1e-11

/// Dot product two vectors.
double ddot(int n, double *x, int incx, double *y, int incy)
{
    double sum = 0.0;

    while (n--)
    {
        sum += (*x) * (*y);
        x += incx;
        y += incy;
    }

    return sum;
}

// Storage
std::vector<double> z(NPUPPER), w(NPUPPER), p(NPUPPER);
std::vector<double> d(NPUPPER *NPUPPER), q(NPUPPER *NPUPPER);

/**
 * @brief Test integrals of Jacobi polynomials.
 *
 * This routine evaluates the integral \f[ \int_{-1}^1 (1-x)^\alpha (1+x)^\beta
 * P_n^{\alpha,\beta} dx = 0 \f] using \f$ -0.5 \leq \alpha,\beta \leq 5 \f$ and
 * using \f$ N \f$ points where \f$ N \f$ lies between #NPLOWER and #NPUPPER.
 * Tolerance is taken using the #EPS value.
 */
inline void TestIntegral(
    std::function<void(double *, double *, int, double, double)> func, int o)
{
    double alpha = -0.5;
    while (alpha <= 5.0)
    {
        double beta = -0.5;
        while (beta <= 5.0)
        {
            for (int np = NPLOWER; np <= NPUPPER; ++np)
            {
                func(&z[0], &w[0], np, alpha, beta);
                for (int n = 2; n < 2 * np - 1 - o; ++n)
                {
                    Polylib::jacobfd(np, &z[0], &p[0], NULL, n, alpha, beta);
                    double sum = ddot(np, &w[0], 1, &p[0], 1);
                    BOOST_CHECK_SMALL(sum, EPS);
                }
            }

            beta += 0.5;
        }

        alpha += 0.5;
    }
}

/**
 * @brief Test derivatives using Gaussian quadrature.
 *
 * This routine evaluates the deriatives \f[ \frac{d}{dx}(x^n) = nx^{n-1} \f]
 * using \f$ -0.5 \leq \alpha,\beta \leq 5 \f$ and using \f$ N \f$ points where
 * \f$ N \f$ lies between #NPLOWER and #NPUPPER.  Tolerance is taken using the
 * #EPS value.
 */
inline void TestDifferentiation(
    std::function<void(double *, double *, int, double, double)> func,
    std::function<void(double *, double *, int, double, double)> funcD)
{
    double alpha = -0.5;
    while (alpha <= 5.0)
    {
        double beta = -0.5;
        while (beta <= 5.0)
        {
            for (int np = NPLOWER; np <= NPUPPER; ++np)
            {
                func(&z[0], &w[0], np, alpha, beta);
                funcD(&d[0], &z[0], np, alpha, beta);
                for (int n = 2; n < np - 1; ++n)
                {
                    for (int i = 0; i < np; ++i)
                    {
                        p[i] = pow(z[i], n);
                    }

                    double sum = 0.0;

                    for (int i = 0; i < np; ++i)
                    {
                        sum += fabs(ddot(np, &d[0] + i, np, &p[0], 1) -
                                    n * pow(z[i], n - 1));
                    }
                    sum /= np;

                    BOOST_CHECK_SMALL(sum, EPS);
                }
            }
            beta += 0.5;
        }
        alpha += 0.5;
    }
}

/**
 * @brief Evaluate interpolation using Gaussian quadrature.
 *
 * This routine evaluates the interpolation of \f$ z^n \f$ to \f$ x^n \f$, where
 * \f$ z \f$ are the quadrature zeros and \f$ x \f$ are the equispaced points
 * \f[ x_n = \frac{2n}{N-1} - 1, \quad 0\leq n\leq N, \f] using \f$ -0.5 \leq
 * \alpha,\beta \leq 5 \f$ and using \f$ N \f$ points where \f$ N \f$ lies
 * between #NPLOWER and #NPUPPER.  Tolerance is taken using the
 * #EPS value.
 */
inline void TestInterpolation(
    std::function<void(double *, double *, int, double, double)> func,
    std::function<void(double *, double *, double *, int, int, double, double)>
        funcI)
{
    double alpha = -0.5;
    while (alpha <= 5.0)
    {
        double beta = -0.5;
        while (beta <= 5.0)
        {
            for (int np = NPLOWER; np <= NPUPPER; ++np)
            {
                func(&z[0], &w[0], np, alpha, beta);
                for (int n = 2; n < np - 1; ++n)
                {
                    for (int i = 0; i < np; ++i)
                    {
                        w[i] = 2.0 * i / (double)(np - 1) - 1.0;
                        p[i] = pow(z[i], n);
                    }
                    funcI(&d[0], &z[0], &w[0], np, np, alpha, beta);

                    double sum = 0.0;
                    for (int i = 0; i < np; ++i)
                    {
                        sum += fabs(ddot(np, &d[0] + i, np, &p[0], 1) -
                                    pow(w[i], n));
                    }
                    sum /= np;

                    BOOST_CHECK_SMALL(sum, EPS);
                }
            }
            beta += 0.5;
        }
        alpha += 0.5;
    }
}

inline void TestIntegralMatrix(
    std::function<void(double *, double *, int, double, double)> func,
    std::function<void(double *, double *, int, int)> funcQ)
{
    double alpha = -0.5;
    while (alpha <= 5.0)
    {
        double beta = -0.5;
        while (beta <= 5.0)
        {
            for (int np = NPLOWER; np <= NPUPPER; ++np)
            {
                func(&z[0], &w[0], np, alpha, beta);
                funcQ(&q[0], &z[0], np, 0);
                for (int n = 2; n < np - 1; ++n)
                {
                    for (int i = 0; i < np; ++i)
                    {
                        p[i] = pow(z[i], n);
                    }

                    double sum = 0.0;

                    for (int i = 0; i < np; ++i)
                    {
                        sum += fabs(ddot(np, &q[0] + i * np, 1, &p[0], 1) -
                                    pow(z[i], n + 1) / (n + 1));
                    }
                    sum /= np;

                    BOOST_CHECK_SMALL(sum, EPS);
                }
            }
            beta += 0.5;
        }
        alpha += 0.5;
    }
}

BOOST_AUTO_TEST_CASE(TestGaussInt)
{
    TestIntegral(Polylib::zwgj, 0);
}

BOOST_AUTO_TEST_CASE(TestGaussRadauM)
{
    TestIntegral(Polylib::zwgrjm, 1);
}

BOOST_AUTO_TEST_CASE(TestGaussRadauP)
{
    TestIntegral(Polylib::zwgrjp, 1);
}

BOOST_AUTO_TEST_CASE(TestGaussLobatto)
{
    TestIntegral(Polylib::zwglj, 2);
}

BOOST_AUTO_TEST_CASE(TestGaussDiff)
{
    TestDifferentiation(Polylib::zwgj, Polylib::Dgj);
}

BOOST_AUTO_TEST_CASE(TestGaussRadauMDiff)
{
    TestDifferentiation(Polylib::zwgrjm, Polylib::Dgrjm);
}

BOOST_AUTO_TEST_CASE(TestGaussRadauPDiff)
{
    TestDifferentiation(Polylib::zwgrjp, Polylib::Dgrjp);
}

BOOST_AUTO_TEST_CASE(TestGaussLobattoDiff)
{
    TestDifferentiation(Polylib::zwglj, Polylib::Dglj);
}

BOOST_AUTO_TEST_CASE(TestGaussInterp)
{
    TestInterpolation(Polylib::zwgj, Polylib::Imgj);
}

BOOST_AUTO_TEST_CASE(TestGaussRadauMInterp)
{
    TestInterpolation(Polylib::zwgrjm, Polylib::Imgrjm);
}

BOOST_AUTO_TEST_CASE(TestGaussRadauPInterp)
{
    TestInterpolation(Polylib::zwgrjp, Polylib::Imgrjp);
}

BOOST_AUTO_TEST_CASE(TestGaussLobattoInterp)
{
    TestInterpolation(Polylib::zwglj, Polylib::Imglj);
}

BOOST_AUTO_TEST_CASE(TestGaussIntMatrix)
{
    TestIntegralMatrix(Polylib::zwgj, Polylib::Qg);
}

BOOST_AUTO_TEST_CASE(TestGaussRadauMIntMatrix)
{
    TestIntegralMatrix(Polylib::zwgrjm, Polylib::Qg);
}

BOOST_AUTO_TEST_CASE(TestGaussRadauPIntMatrix)
{
    TestIntegralMatrix(Polylib::zwgrjp, Polylib::Qg);
}

BOOST_AUTO_TEST_CASE(TestGaussLobattoIntMatrix)
{
    TestIntegralMatrix(Polylib::zwglj, Polylib::Qg);
}

BOOST_AUTO_TEST_CASE(TestGammaFraction)
{
    double a, c;

    using arrT = std::tuple<int, double, int, double>;
    std::vector<std::tuple<double, arrT, arrT>> ac = {
        {362880.0, std::make_tuple(10, 0.0, 0, 1.0),
         std::make_tuple(1, 0.0, 12, -2.0)},
        {30.0, std::make_tuple(4, 3., 5, 0.), std::make_tuple(5, 0., 7, 0.)},
        {21651.09375, std::make_tuple(7, 3.5, 5, 0.5),
         std::make_tuple(5, 0.5, 10, 0.5)},
        {97429.921875, std::make_tuple(10, 0.5, 5, -0.5),
         std::make_tuple(5, -0.5, 11, -0.5)},
        {2279.0625, std::make_tuple(10, -0.5, 5, 0.5),
         std::make_tuple(5, 0.5, 10, -0.5)},
        {639383.8623046875, std::make_tuple(10, 0.5, 0, 0.5),
         std::make_tuple(0, 0.5, 10, 0.5)},
        {5967498288235982848., std::make_tuple(100, 0., 90, 0.5),
         std::make_tuple(90, 0.5, 100, 0.)},
        {5618641603777298169856., std::make_tuple(200, 0., 191, -0.5),
         std::make_tuple(190, 0.5, 200, 0.)},
        {77396694214720029196288., std::make_tuple(200, 0., 190, 0.),
         std::make_tuple(190, 0., 200, 0.)}};

    for (auto &test : ac)
    {
        double a   = std::get<0>(test);
        arrT c1    = std::get<1>(test);
        arrT c2    = std::get<2>(test);
        double c1e = Polylib::gammaFracGammaF(std::get<0>(c1), std::get<1>(c1),
                                              std::get<2>(c1), std::get<3>(c1));
        double c2e = Polylib::gammaFracGammaF(std::get<0>(c2), std::get<1>(c2),
                                              std::get<2>(c2), std::get<3>(c2));

        BOOST_CHECK_SMALL(c1e / a - 1.0, EPS);
        BOOST_CHECK_SMALL(c2e * a - 1.0, EPS);
    }

    int Q = 2;
    for (double alpha = -0.5; alpha <= 150.; alpha += 0.5)
    {
        for (double beta = -0.5; beta <= 150.; beta += 0.5)
        {
            a = Polylib::gammaF(Q + alpha) / Polylib::gammaF(Q + beta);
            c = Polylib::gammaFracGammaF(Q, alpha, Q, beta) / a - 1.;
            BOOST_CHECK_SMALL(c, EPS);
        }
    }
}

} // namespace UnitTests
} // namespace Nektar
