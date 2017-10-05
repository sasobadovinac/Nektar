///////////////////////////////////////////////////////////////////////////////
//
// File: ConjugateGradient_BLAS.hxx
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_CONJUGATEGRADIENT_BLAS_HXX
#define NEKTAR_LIB_MULTIREGIONS_CONJUGATEGRADIENT_BLAS_HXX

#include <MultiRegions/GlobalLinSysIterative.h>

namespace Nektar
{
namespace MultiRegions
{

    KOKKOS_INLINE_FUNCTION
    double plainDdot(int n, const double *dx, int incx,
            const double *dy, int incy)
    {
        /* Local variables */
        int i, m, ix, iy, mp1;
        double dtemp;

        /* Parameter adjustments */
        --dy;
        --dx;

        /* Function Body */
        dtemp = 0.0;
        if (n <= 0)
        {
            return dtemp;
        }
        if (incx == 1 && incy == 1)
        {
            m = n % 4;
            if (m != 0)
            {
                for (i = 1; i <= m; ++i)
                {
                    dtemp += dx[i] * dy[i];
                }
                if (n < 4)
                {
                    return dtemp;
                }
            }
            mp1 = m + 1;
            for (i = mp1; i <= n; i += 4)
            {
                dtemp = dtemp + dx[i] * dy[i]
                              + dx[i + 1] * dy[i + 1]
                              + dx[i + 2] * dy[i + 2]
                              + dx[i + 3] * dy[i + 3];
            }
            return dtemp;
        }
        else
        {                
            ix = 1;
            iy = 1;
            if (incx < 0)
            {
                ix = (-n + 1) * incx + 1;
            }
            if (incy < 0)
            {
                iy = (-n + 1) * incy + 1;
            }
            for (i = 1; i <= n; ++i)
            {
                dtemp += dx[ix] * dy[iy];
                ix += incx;
                iy += incy;
            }
            return dtemp;
        }

    } /* ddot_ */


    KOKKOS_INLINE_FUNCTION
    int plainDaxpy(int n, const double da, const double *dx,
            int incx, double *dy, int incy)
    {
        /* Local variables */
        int i, m, ix, iy, mp1;

        /* Parameter adjustments */
        --dy;
        --dx;

        /* Function Body */
        if (n <= 0)
        {
            return 0;
        }
        if (da == 0.0)
        {
            return 0;
        }
        if (incx == 1 && incy == 1)
        {
            m = n % 4;
            if (m != 0)
            {
                for (i = 1; i <= m; ++i)
                {
                    dy[i] += da * dx[i];
                }
                if (n < 4)
                {
                    return 0;
                }
            }
            mp1 = m + 1;
            for (i = mp1; i <= n; i += 4)
            {
                dy[i] += da * dx[i];
                dy[i + 1] += da * dx[i + 1];
                dy[i + 2] += da * dx[i + 2];
                dy[i + 3] += da * dx[i + 3];
            }
            return 0;
        }
        else
        {
            ix = 1;
            iy = 1;
            if (incx < 0)
            {
                ix = (-n + 1) * incx + 1;
            }
            if (incy < 0) 
            {
                iy = (-n + 1) * incy + 1;
            }
            for (i = 1; i <= n; ++i) 
            {
                dy[iy] += da * dx[ix];
                ix += incx;
                iy += incy;
            }
            return 0;
        }
    } /* daxpy_ */


    KOKKOS_INLINE_FUNCTION
    int plainDgemv(char trans, int m, int n, const double alpha,
            const double *a, int lda, const double *x, int incx,
            const double beta, double *y, int incy)
    {
        /* Local variables */
        int i, j, ix, iy, jx, jy, kx, ky;
        double temp;
        int lenx, leny;

        /* Parameter adjustments */
        a -= 1 + lda;
        --x;
        --y;

        /* Function Body */
    /*     Quick return if possible. */

        if (m == 0 || n == 0 || alpha == 0. && beta == 1.) {
        return 0;
        }

    /*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
    /*     up the start points in  X  and  Y. */

        if (trans == 'N'){
        lenx = n;
        leny = m;
        } else {
        lenx = m;
        leny = n;
        }
        if (incx > 0) {
        kx = 1;
        } else {
        kx = 1 - (lenx - 1) * incx;
        }
        if (incy > 0) {
        ky = 1;
        } else {
        ky = 1 - (leny - 1) * incy;
        }

    /*     Start the operations. In this version the elements of A are */
    /*     accessed sequentially with one pass through A. */

    /*     First form  y := beta*y. */

        if (beta != 1.) {
        if (incy == 1) {
            if (beta == 0.) {
            for (i = 1; i <= leny; ++i) {
                y[i] = 0.;
    /* L10: */
            }
            } else {
            for (i = 1; i <= leny; ++i) {
                y[i] = beta * y[i];
    /* L20: */
            }
            }
        } else {
            iy = ky;
            if (beta == 0.) {
            for (i = 1; i <= leny; ++i) {
                y[iy] = 0.;
                iy += incy;
    /* L30: */
            }
            } else {
            for (i = 1; i <= leny; ++i) {
                y[iy] = beta * y[iy];
                iy += incy;
    /* L40: */
            }
            }
        }
        }
        if (alpha == 0.) {
        return 0;
        }
        if ((trans == 'N')) {

    /*        Form  y := alpha*A*x + y. */

        jx = kx;
        if (incy == 1) {
            for (j = 1; j <= n; ++j) {
            if (x[jx] != 0.) {
                temp = alpha * x[jx];
                for (i = 1; i <= m; ++i) {
                y[i] += temp * a[i + j * lda];
    /* L50: */
                }
            }
            jx += incx;
    /* L60: */
            }
        } else {
            for (j = 1; j <= n; ++j) {
            if (x[jx] != 0.) {
                temp = alpha * x[jx];
                iy = ky;
                for (i = 1; i <= m; ++i) {
                y[iy] += temp * a[i + j * lda];
                iy += incy;
    /* L70: */
                }
            }
            jx += incx;
    /* L80: */
            }
        }
        } else {

    /*        Form  y := alpha*A'*x + y. */

        jy = ky;
        if (incx == 1) {
            for (j = 1; j <= n; ++j) {
            temp = 0.;
            for (i = 1; i <= m; ++i) {
                temp += a[i + j * lda] * x[i];
    /* L90: */
            }
            y[jy] += alpha * temp;
            jy += incy;
    /* L100: */
            }
        } else {
            for (j = 1; j <= n; ++j) {
            temp = 0.;
            ix = kx;
            for (i = 1; i <= m; ++i) {
                temp += a[i + j * lda] * x[ix];
                ix += incx;
    /* L110: */
            }
            y[jy] += alpha * temp;
            jy += incy;
    /* L120: */
            }
        }
        }

        return 0;
    } /* dgemv_ */


    KOKKOS_INLINE_FUNCTION
    int plainDgemm(char transa, char transb, int m, int n, int k,
        const double alpha, const double *a, int lda, const double *b,
        int ldb, const double beta, double *c, int ldc)
    {
        /* Local variables */
        int i, j, l;
        bool nota, notb;
        double temp;

        /* Parameter adjustments */
        a -= 1 + lda;
        b -= 1 + ldb;
        c -= 1 + ldc;

        /* Function Body */
        nota = (transa == 'N');
        notb = (transb == 'N');       
    

    /*     Quick return if possible. */

        if (m == 0 || n == 0 || (alpha == 0. || k == 0) && beta == 1.) {
        return 0;
        }

    /*     And if  alpha.eq.zero. */

        if (alpha == 0.) {
        if (beta == 0.) {
            for (j = 1; j <= n; ++j) {            
            for (i = 1; i <= m; ++i) {
                c[i + j * ldc] = 0.;
    /* L10: */
            }
    /* L20: */
            }
        } else {
            for (j = 1; j <= n; ++j) {            
            for (i = 1; i <= m; ++i) {
                c[i + j * ldc] = beta * c[i + j * ldc];
    /* L30: */
            }
    /* L40: */
            }
        }
        return 0;
        }

    /*     Start the operations. */

        if (notb) {
        if (nota) {

    /*           Form  C := alpha*A*B + beta*C. */

            for (j = 1; j <= n; ++j) {
            if (beta == 0.) {                
                for (i = 1; i <= m; ++i) {
                c[i + j * ldc] = 0.;
    /* L50: */
                }
            } else if (beta != 1.) {                
                for (i = 1; i <= m; ++i) {
                c[i + j * ldc] = beta * c[i + j * ldc];
    /* L60: */
                }
            }
            for (l = 1; l <= k; ++l) {
                if (b[l + j * ldb] != 0.) {
                temp = alpha * b[l + j * ldb];
                for (i = 1; i <= m; ++i) {
                    c[i + j * ldc] += temp * a[i + l * 
                        lda];
    /* L70: */
                }
                }
    /* L80: */
            }
    /* L90: */
            }
        } else {

    /*           Form  C := alpha*A'*B + beta*C */

            for (j = 1; j <= n; ++j) {            
            for (i = 1; i <= m; ++i) {
                temp = 0.;
                for (l = 1; l <= k; ++l) {
                temp += a[l + i * lda] * b[l + j * ldb];
    /* L100: */
                }
                if (beta == 0.) {
                c[i + j * ldc] = alpha * temp;
                } else {
                c[i + j * ldc] = alpha * temp + beta * c[
                    i + j * ldc];
                }
    /* L110: */
            }
    /* L120: */
            }
        }
        } else {
        if (nota) {

    /*           Form  C := alpha*A*B' + beta*C */

            for (j = 1; j <= n; ++j) {
            if (beta == 0.) {                
                for (i = 1; i <= m; ++i) {
                c[i + j * ldc] = 0.;
    /* L130: */
                }
            } else if (beta != 1.) {                
                for (i = 1; i <= m; ++i) {
                c[i + j * ldc] = beta * c[i + j * ldc];
    /* L140: */
                }
            }
            for (l = 1; l <= k; ++l) {
                if (b[j + l * ldb] != 0.) {
                temp = alpha * b[j + l * ldb];
                for (i = 1; i <= m; ++i) {
                    c[i + j * ldc] += temp * a[i + l * 
                        lda];
    /* L150: */
                }
                }
    /* L160: */
            }
    /* L170: */
            }
        } else {

    /*           Form  C := alpha*A'*B' + beta*C */

            for (j = 1; j <= n; ++j) {            
            for (i = 1; i <= m; ++i) {
                temp = 0.;
                for (l = 1; l <= k; ++l) {
                temp += a[l + i * lda] * b[j + l * ldb];
    /* L180: */
                }
                if (beta == 0.) {
                c[i + j * ldc] = alpha * temp;
                } else {
                c[i + j * ldc] = alpha * temp + beta * c[
                    i + j * ldc];
                }
    /* L190: */
            }
    /* L200: */
            }
        }
        }

        return 0;

    /*     End of DGEMM . */

    } /* dgemm_ */


}
}

#endif