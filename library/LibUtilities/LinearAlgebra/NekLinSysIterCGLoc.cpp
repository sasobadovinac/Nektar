///////////////////////////////////////////////////////////////////////////////
//
// File: NekLinSysIterCGLoc.cpp
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
// Description: NekLinSysIterCGLoc definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/LinearAlgebra/NekLinSysIterCGLoc.h>

using namespace std;

namespace Nektar
{
namespace LibUtilities
{
/**
 * @class  NekLinSysIterCGLoc
 *
 * Solves a linear system using iterative methods.
 */
string NekLinSysIterCGLoc::className =
    LibUtilities::GetNekLinSysIterFactory().RegisterCreatorFunction(
        "ConjugateGradientLoc", NekLinSysIterCGLoc::create,
        "NekLinSysIterCG solver in Local Space.");

NekLinSysIterCGLoc::NekLinSysIterCGLoc(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::CommSharedPtr &vComm, const int nDimen,
    const NekSysKey &pKey)
    : NekLinSysIter(pSession, vComm, nDimen, pKey)
{
    m_isLocal = true;
}

void NekLinSysIterCGLoc::v_InitObject()
{
    NekLinSysIter::v_InitObject();
}

NekLinSysIterCGLoc::~NekLinSysIterCGLoc()
{
}

/**
 *
 */
int NekLinSysIterCGLoc::v_SolveSystem(
    const int nLocal, const Array<OneD, const NekDouble> &pInput,
    Array<OneD, NekDouble> &pOutput, const int nDir, const NekDouble tol,
    const NekDouble factor)
{
    boost::ignore_unused(tol, nDir);

    m_tolerance   = max(tol, 1.0E-16);
    m_prec_factor = factor;

    DoConjugateGradient(nLocal, pInput, pOutput);

    return m_totalIterations;
}

/**  
 * Solve a global linear system using the conjugate gradient method.  
 * We solve only for the non-Dirichlet modes. The operator is evaluated  
 * using an auxiliary function m_operator.DoNekSysLhsEval defined by the  
 * specific solver. Distributed math routines are used to support  
 * parallel execution of the solver.  
 *  
 * The implemented algorithm uses a reduced-communication reordering of  
 * the standard PCG method (Demmel, Heath and Vorst, 1993)  
 *  
 * @param       pInput      Input residual  of all DOFs.  
 * @param       pOutput     Solution vector of all DOFs.  
 */
void NekLinSysIterCGLoc::DoConjugateGradient(
    const int nLocal, const Array<OneD, const NekDouble> &pInput,
    Array<OneD, NekDouble> &pOutput)
{
    // Allocate array storage
    Array<OneD, NekDouble> w_A(nLocal, 0.0);
    Array<OneD, NekDouble> s_A(nLocal, 0.0);
    Array<OneD, NekDouble> p_A(nLocal, 0.0);
    Array<OneD, NekDouble> r_A(nLocal, 0.0);
    Array<OneD, NekDouble> q_A(nLocal, 0.0);
    Array<OneD, NekDouble> wk(nLocal, 0.0);

    int k;
    NekDouble alpha;
    NekDouble beta;
    NekDouble rho;
    NekDouble rho_new;
    NekDouble mu;
    NekDouble eps;
    NekDouble min_resid;
    Array<OneD, NekDouble> vExchange(3, 0.0);

    // Copy initial residual from input
    Vmath::Vcopy(nLocal, pInput, 1, r_A, 1);

    // zero homogeneous out array ready for solution updates
    // Should not be earlier in case input vector is same as
    // output and above copy has been peformed
    Vmath::Zero(nLocal, pOutput, 1);

    // evaluate initial residual error for exit check
    m_operator.DoAssembleLoc(r_A, wk, true);
    vExchange[2] = Vmath::Dot(nLocal, wk, r_A);
    m_Comm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

    eps = vExchange[2];

    if (m_rhs_magnitude == NekConstants::kNekUnsetDouble)
    {
        Set_Rhs_Magnitude(pInput);
    }

    m_totalIterations = 0;

    // If input residual is less than tolerance skip solve.
    if (eps < m_tolerance * m_tolerance * m_rhs_magnitude)
    {
        if (m_verbose && m_root)
        {
            cout << "CG iterations made = " << m_totalIterations
                 << " using tolerance of " << m_tolerance
                 << " (error = " << sqrt(eps / m_rhs_magnitude)
                 << ", rhs_mag = " << sqrt(m_rhs_magnitude) << ")" << endl;
        }
        return;
    }

    m_operator.DoNekSysPrecon(r_A, w_A, true);
    m_operator.DoNekSysLhsEval(w_A, s_A);

    k = 0;

    vExchange[0] = Vmath::Dot(nLocal, r_A, w_A);
    vExchange[1] = Vmath::Dot(nLocal, s_A, w_A);
    m_Comm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

    rho               = vExchange[0];
    mu                = vExchange[1];
    min_resid         = m_rhs_magnitude;
    beta              = 0.0;
    alpha             = rho / mu;
    m_totalIterations = 1;

    // Continue until convergence
    while (true)
    {
        if (k >= m_maxiter)
        {
            if (m_root)
            {
                cout << "CG iterations made = " << m_totalIterations
                     << " using tolerance of " << m_tolerance
                     << " (error = " << sqrt(eps / m_rhs_magnitude)
                     << ", rhs_mag = " << sqrt(m_rhs_magnitude) << ")" << endl;
            }
            ROOTONLY_NEKERROR(ErrorUtil::efatal,
                              "Exceeded maximum number of iterations");
        }

        // Compute new search direction p_k, q_k
        Vmath::Svtvp(nLocal, beta, p_A, 1, w_A, 1, p_A, 1);
        Vmath::Svtvp(nLocal, beta, q_A, 1, s_A, 1, q_A, 1);

        // Update solution x_{k+1}
        Vmath::Svtvp(nLocal, alpha, p_A, 1, pOutput, 1, pOutput, 1);

        // Update residual vector r_{k+1}
        Vmath::Svtvp(nLocal, -alpha, q_A, 1, r_A, 1, r_A, 1);

        // Apply preconditioner
        m_operator.DoNekSysPrecon(r_A, w_A, true);

        // Perform the method-specific matrix-vector multiply operation.
        m_operator.DoNekSysLhsEval(w_A, s_A);

        // <r_{k+1}, w_{k+1}>
        vExchange[0] = Vmath::Dot(nLocal, r_A, w_A);
        // <s_{k+1}, w_{k+1}>
        vExchange[1] = Vmath::Dot(nLocal, s_A, w_A);
        // <r_{k+1}, r_{k+1}>
        m_operator.DoAssembleLoc(r_A, wk, true);
        vExchange[2] = Vmath::Dot(nLocal, wk, r_A);

        // Perform inner-product exchanges
        m_Comm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);

        rho_new = vExchange[0];
        mu      = vExchange[1];
        eps     = vExchange[2];

        m_totalIterations++;

        // Test if norm is within tolerance
        if (eps < m_tolerance * m_tolerance * m_rhs_magnitude)
        {
            if (m_verbose && m_root)
            {
                cout << "CG iterations made = " << m_totalIterations
                     << " using tolerance of " << m_tolerance
                     << " (error = " << sqrt(eps / m_rhs_magnitude)
                     << ", rhs_mag = " << sqrt(m_rhs_magnitude) << ")" << endl;
            }
            break;
        }
        min_resid = min(min_resid, eps);

        // Compute search direction and solution coefficients
        beta  = rho_new / rho;
        alpha = rho_new / (mu - rho_new * beta / alpha);
        rho   = rho_new;
        k++;
    }
}
} // namespace LibUtilities
} // namespace Nektar
