///////////////////////////////////////////////////////////////////////////////
//
// File: Polylib.h
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef H_PLYLIB

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <complex>

/*
 *  LIBRARY ROUTINES FOR POLYNOMIAL CALCULUS AND INTERPOLATION
 */

/** \brief The namespace associated with the the Polylib library
 * (\ref pagePolylib "Polylib introduction")
 */
namespace Polylib
{

/**
   \page pagePolylib The Polylib library
   \section sectionPolyLib Routines For Orthogonal Polynomial Calculus and
   Interpolation

   Spencer Sherwin,
   Aeronautics, Imperial College London

   Based on codes by Einar Ronquist and Ron Henderson

   Abbreviations
   - z    -   Set of collocation/quadrature points
   - w    -   Set of quadrature weights
   - D    -   Derivative matrix
   - h    -   Lagrange Interpolant
   - I    -   Interpolation matrix
   - g    -   Gauss
   - k    -   Kronrod
   - gr   -   Gauss-Radau
   - gl   -   Gauss-Lobatto
   - j    -   Jacobi
   - m    -   point at minus 1 in Radau rules
   - p    -   point at plus  1 in Radau rules

   -----------------------------------------------------------------------\n
   MAIN     ROUTINES\n
   -----------------------------------------------------------------------\n

   Points and Weights:

   - zwgj        Compute Gauss-Jacobi         points and weights
   - zwgrjm      Compute Gauss-Radau-Jacobi   points and weights (z=-1)
   - zwgrjp      Compute Gauss-Radau-Jacobi   points and weights (z= 1)
   - zwglj       Compute Gauss-Lobatto-Jacobi points and weights
       - zwgk        Compute Gauss-Kronrod-Jacobi points and weights
       - zwrk        Compute Radau-Kronrod        points and weights
       - zwlk        Compute Lobatto-Kronrod      points and weights
       - JacZeros    Compute Gauss-Jacobi         points and weights

   Integration Matrices:

   - Qg          Compute Gauss               integration matrix

   Derivative Matrices:

   - Dgj         Compute Gauss-Jacobi         derivative matrix
   - Dgrjm       Compute Gauss-Radau-Jacobi   derivative matrix (z=-1)
   - Dgrjp       Compute Gauss-Radau-Jacobi   derivative matrix (z= 1)
   - Dglj        Compute Gauss-Lobatto-Jacobi derivative matrix

   Lagrange Interpolants:

   - hgj         Compute Gauss-Jacobi         Lagrange interpolants
   - hgrjm       Compute Gauss-Radau-Jacobi   Lagrange interpolants (z=-1)
   - hgrjp       Compute Gauss-Radau-Jacobi   Lagrange interpolants (z= 1)
   - hglj        Compute Gauss-Lobatto-Jacobi Lagrange interpolants

   Interpolation Operators:

   - Imgj        Compute interpolation operator gj->m
   - Imgrjm      Compute interpolation operator grj->m (z=-1)
   - Imgrjp      Compute interpolation operator grj->m (z= 1)
   - Imglj       Compute interpolation operator glj->m

   Polynomial Evaluation:

   - polycoeff   Returns value and derivative of Jacobi poly.
   - jacobfd     Returns value and derivative of Jacobi poly. at point z
   - jacobd      Returns derivative of Jacobi poly. at point z (valid at z=-1,1)

   -----------------------------------------------------------------------\n
   LOCAL      ROUTINES\n
   -----------------------------------------------------------------------\n

   - jacobz      Returns Jacobi polynomial zeros
   - gammaf      Gamma function for integer values and halves
       - RecCoeff    Calculates the recurrence coefficients for orthogonal poly
       - TriQL		 QL algorithm for symmetrix tridiagonal matrix
       - JKMatrix	 Generates the Jacobi-Kronrod matrix

   ------------------------------------------------------------------------\n

   Useful references:

   - [1] Gabor Szego: Orthogonal Polynomials, American Mathematical Society,
   Providence, Rhode Island, 1939.
   - [2] Abramowitz \& Stegun: Handbook of Mathematical Functions,
   Dover, New York, 1972.
   - [3] Canuto, Hussaini, Quarteroni \& Zang: Spectral Methods in Fluid
   Dynamics, Springer-Verlag, 1988.
   - [4] Ghizzetti \& Ossicini: Quadrature Formulae, Academic Press, 1970.
   - [5] Karniadakis \& Sherwin: Spectral/hp element methods for CFD, 1999


   NOTES

   -# Legendre  polynomial \f$ \alpha = \beta = 0 \f$
   -# Chebychev polynomial \f$ \alpha = \beta = -0.5 \f$
   -# All routines are double precision.
   -# All array subscripts start from zero, i.e. vector[0..N-1]
*/

/*-----------------------------------------------------------------------
  M A I N     R O U T I N E S
  -----------------------------------------------------------------------*/

/* Points and weights */
LIB_UTILITIES_EXPORT void zwgj(double *, double *, const int, const double,
                               const double);
LIB_UTILITIES_EXPORT void zwgrjm(double *, double *, const int, const double,
                                 const double);
LIB_UTILITIES_EXPORT void zwgrjp(double *, double *, const int, const double,
                                 const double);
LIB_UTILITIES_EXPORT void zwglj(double *, double *, const int, const double,
                                const double);
LIB_UTILITIES_EXPORT void zwgk(double *, double *, const int, const double,
                               const double);
LIB_UTILITIES_EXPORT void zwrk(double *, double *, const int, const double,
                               const double);
LIB_UTILITIES_EXPORT void zwlk(double *, double *, const int, const double,
                               const double);
LIB_UTILITIES_EXPORT void JacZeros(const int, double *, double *, const double,
                                   const double);

/* Integration operators */
LIB_UTILITIES_EXPORT void Qg(double *, const double *, const int);

/* Derivative operators */
LIB_UTILITIES_EXPORT void Dgj(double *, const double *, const int, const double,
                              const double);
LIB_UTILITIES_EXPORT void Dgrjm(double *, const double *, const int,
                                const double, const double);
LIB_UTILITIES_EXPORT void Dgrjp(double *, const double *, const int,
                                const double, const double);
LIB_UTILITIES_EXPORT void Dglj(double *, const double *, const int,
                               const double, const double);

/* Lagrangian interpolants */
LIB_UTILITIES_EXPORT double laginterp(double, int, const double *, int);
LIB_UTILITIES_EXPORT double laginterpderiv(double, int, const double *, int);
LIB_UTILITIES_EXPORT double hgj(const int, const double, const double *,
                                const int, const double, const double);
LIB_UTILITIES_EXPORT double hgrjm(const int, const double, const double *,
                                  const int, const double, const double);
LIB_UTILITIES_EXPORT double hgrjp(const int, const double, const double *,
                                  const int, const double, const double);
LIB_UTILITIES_EXPORT double hglj(const int, const double, const double *,
                                 const int, const double, const double);

/* Interpolation operators */
LIB_UTILITIES_EXPORT void Imgj(double *, const double *, const double *,
                               const int, const int, const double,
                               const double);
LIB_UTILITIES_EXPORT void Imgrjm(double *, const double *, const double *,
                                 const int, const int, const double,
                                 const double);
LIB_UTILITIES_EXPORT void Imgrjp(double *, const double *, const double *,
                                 const int, const int, const double,
                                 const double);
LIB_UTILITIES_EXPORT void Imglj(double *, const double *, const double *,
                                const int, const int, const double,
                                const double);

/* Polynomial functions */
LIB_UTILITIES_EXPORT void polycoeffs(double *, const double *, const int,
                                     const int);
LIB_UTILITIES_EXPORT void jacobfd(const int, const double *, double *, double *,
                                  const int, const double, const double);
LIB_UTILITIES_EXPORT void jacobd(const int, const double *, double *, const int,
                                 const double, const double);

/* Gamma function routines */
LIB_UTILITIES_EXPORT double gammaF(const double);
LIB_UTILITIES_EXPORT double gammaFracGammaF(const int, const double, const int,
                                            const double);

/* Bessel function routines */
LIB_UTILITIES_EXPORT std::complex<Nektar::NekDouble> ImagBesselComp(
    const int, std::complex<Nektar::NekDouble>);

} // namespace Polylib

#define H_PLYLIB
#endif /* END OF POLYLIB.H DECLARATIONS */
