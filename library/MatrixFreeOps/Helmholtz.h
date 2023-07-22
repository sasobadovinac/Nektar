///////////////////////////////////////////////////////////////////////////////
//
// File: Helmholtz.h
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

#ifndef NEKTAR_LIBRARY_MF_HELMHOLTZ_H
#define NEKTAR_LIBRARY_MF_HELMHOLTZ_H

#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "Operator.hpp"

#include "BwdTransKernels.hpp"
#include "HelmholtzKernels.hpp"
#include "IProductKernels.hpp"
#include "PhysDerivKernels.hpp"

namespace Nektar
{
namespace MatrixFree
{

// As each opertor has seven shapes over three dimension to get to the
// "work" each operator uses a series of preprocessor directives based
// on the dimension and shape type so to limit the code
// inclusion. This keeps the library size as minimal as possible while
// removing runtime conditionals.

// The preprocessor directives, SHAPE_DIMENSION_?D and SHAPE_TYPE_*
// are constructed by CMake in the CMakeLists.txt file which uses an
// implementation file.  See the CMakeLists.txt files for more
// details.
template <LibUtilities::ShapeType SHAPE_TYPE, bool DEFORMED = false>
struct HelmholtzTemplate
    : public Helmholtz,
      public Helper<LibUtilities::ShapeTypeDimMap[SHAPE_TYPE], DEFORMED>
{
    HelmholtzTemplate(std::vector<LibUtilities::BasisSharedPtr> basis,
                      int nElmt)
        : Helmholtz(basis, nElmt),
          Helper<LibUtilities::ShapeTypeDimMap[SHAPE_TYPE], DEFORMED>(basis,
                                                                      nElmt)
    {
        constexpr auto DIM = LibUtilities::ShapeTypeDimMap[SHAPE_TYPE];

        if (DIM == 1)
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients(SHAPE_TYPE,
                                                            this->m_nm[0]);
        }
        else if (DIM == 2)
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients(
                SHAPE_TYPE, this->m_nm[0], this->m_nm[1]);

#if defined(SHAPE_TYPE_TRI)
            const auto nq0 = m_basis[0]->GetNumPoints();
            const auto nq1 = m_basis[1]->GetNumPoints();

            const Array<OneD, const NekDouble> &z0 = m_basis[0]->GetZ();
            const Array<OneD, const NekDouble> &z1 = m_basis[1]->GetZ();

            GetHelmholtz2DHalfSpace<SHAPE_TYPE>(nq0, nq1, z0, z1, m_h0, m_h1);
#endif
        }
        else if (DIM == 3)
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients(
                SHAPE_TYPE, this->m_nm[0], this->m_nm[1], this->m_nm[2]);

#if defined(SHAPE_TYPE_TET) || defined(SHAPE_TYPE_PRISM) ||                    \
    defined(SHAPE_TYPE_PYR)
            const auto nq0 = m_basis[0]->GetNumPoints();
            const auto nq1 = m_basis[1]->GetNumPoints();
            const auto nq2 = m_basis[2]->GetNumPoints();

            const Array<OneD, const NekDouble> &z0 = m_basis[0]->GetZ();
            const Array<OneD, const NekDouble> &z1 = m_basis[1]->GetZ();
            const Array<OneD, const NekDouble> &z2 = m_basis[2]->GetZ();

            GetHelmholtz3DHalfSpace<SHAPE_TYPE>(nq0, nq1, nq2, z0, z1, z2, m_h0,
                                                m_h1, m_h2, m_h3);
#endif
        }
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
    {
        return std::make_shared<HelmholtzTemplate<SHAPE_TYPE, DEFORMED>>(basis,
                                                                         nElmt);
    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output) final
    {
#include "SwitchNodesPoints.h"
    }

    // There must be separate 1D, 2D, and 3D operators because the
    // helper base class is based on the shape type dim. Thus indices
    // for the basis must be respected.

    // Further there are duplicate 1D, 2D, and 3D operators so to
    // allow for compile time array sizes to be part of the template
    // and thus gain loop unrolling. The only difference is the
    // location of the size value, in the template or the function:
    // foo<int bar>() {} vs foo() { int bar = ... }

#if defined(SHAPE_DIMENSION_1D)

    // Non-size based operator.
    void operator1D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        boost::ignore_unused(input);
        boost::ignore_unused(output);

        ASSERTL0(false, "HelmholtzTemplate::operator1D: Not Impelented.");
    }

    // Size based template version.
    template <int nm0, int nq0>
    void operator1D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        boost::ignore_unused(input);
        boost::ignore_unused(output);

        ASSERTL0(false, "HelmholtzTemplate::operator1D: Not Impelented.");
    }

#elif defined(SHAPE_DIMENSION_2D)

    // Non-size based operator.
    void operator2D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        const auto nm0 = m_basis[0]->GetNumModes();
        const auto nm1 = m_basis[1]->GetNumModes();

        const auto nq0 = m_basis[0]->GetNumPoints();
        const auto nq1 = m_basis[1]->GetNumPoints();

        constexpr auto ndf = 4;

        const auto nqTot    = nq0 * nq1;
        const auto nmBlocks = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Allocate sufficient workspace for backwards transform and inner
        // product kernels.
        size_t wsp0Size = 0;
        BwdTrans2DWorkspace<SHAPE_TYPE>(nm0, nm1, nq0, nq1, wsp0Size);
        IProduct2DWorkspace<SHAPE_TYPE>(nm0, nm1, nq0, nq1, wsp0Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), tmpIn(m_nmTot),
            tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot), deriv0(nqTot),
            deriv1(nqTot);

        const vec_t *jac_ptr;
        const vec_t *df_ptr;

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            df_ptr = &((*this->m_df)[e * dfSize]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            if (!DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e * nqTot]);
            }

            // Step 1: BwdTrans
            BwdTrans2DKernel<SHAPE_TYPE>(nm0, nm1, nq0, nq1, correct, tmpIn,
                                         this->m_bdata[0], this->m_bdata[1],
                                         wsp0, bwd);

            // Step 2: inner product for mass matrix operation
            IProduct2DKernel<SHAPE_TYPE, true, false, DEFORMED>(
                nm0, nm1, nq0, nq1, correct, bwd, this->m_bdata[0],
                this->m_bdata[1], this->m_w[0], this->m_w[1], jac_ptr, wsp0,
                tmpOut, m_lambda);

            // Step 3: take derivatives in collapsed coordinate space
            PhysDerivTensor2DKernel(nq0, nq1, bwd, this->m_D[0], this->m_D[1],
                                    deriv0, deriv1);

            Helmholtz2DKernel<SHAPE_TYPE, DEFORMED>(
                nq0, nq1, this->m_isConstVarDiff, this->m_constVarDiff,
                this->m_isVarDiff, this->m_varD00, this->m_varD01,
                this->m_varD11, df_ptr, this->m_h0, this->m_h1, bwd, deriv0,
                deriv1);

            // Step 4: Apply Laplacian metrics & inner product
            IProduct2DKernel<SHAPE_TYPE, false, true, DEFORMED>(
                nm0, nm1, nq0, nq1, correct, bwd, this->m_dbdata[0],
                this->m_bdata[1], this->m_w[0], this->m_w[1], jac_ptr, wsp0,
                tmpOut);

            IProduct2DKernel<SHAPE_TYPE, false, true, DEFORMED>(
                nm0, nm1, nq0, nq1, correct, deriv0, this->m_bdata[0],
                this->m_dbdata[1], this->m_w[0], this->m_w[1], jac_ptr, wsp0,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nmBlocks;
            outptr += nmBlocks;
        }
    }

    // Size based template version.
    template <int nm0, int nm1, int nq0, int nq1>
    void operator2D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        constexpr auto ndf = 4;

        constexpr auto nqTot = nq0 * nq1;
        const auto nmBlocks  = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Allocate sufficient workspace for backwards transform and inner
        // product kernels.
        size_t wsp0Size = 0;
        BwdTrans2DWorkspace<SHAPE_TYPE>(nm0, nm1, nq0, nq1, wsp0Size);
        IProduct2DWorkspace<SHAPE_TYPE>(nm0, nm1, nq0, nq1, wsp0Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), tmpIn(m_nmTot),
            tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot), deriv0(nqTot),
            deriv1(nqTot);

        const vec_t *jac_ptr;
        const vec_t *df_ptr;

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            df_ptr = &((*this->m_df)[e * dfSize]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            if (!DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e * nqTot]);
            }

            // Step 1: BwdTrans
            BwdTrans2DKernel<SHAPE_TYPE>(nm0, nm1, nq0, nq1, correct, tmpIn,
                                         this->m_bdata[0], this->m_bdata[1],
                                         wsp0, bwd);

            // Step 2: inner product for mass matrix operation
            IProduct2DKernel<SHAPE_TYPE, true, false, DEFORMED>(
                nm0, nm1, nq0, nq1, correct, bwd, this->m_bdata[0],
                this->m_bdata[1], this->m_w[0], this->m_w[1], jac_ptr, wsp0,
                tmpOut, m_lambda);

            // Step 3: take derivatives in collapsed coordinate space
            PhysDerivTensor2DKernel(nq0, nq1, bwd, this->m_D[0], this->m_D[1],
                                    deriv0, deriv1);

            Helmholtz2DKernel<SHAPE_TYPE, DEFORMED>(
                nq0, nq1, this->m_isConstVarDiff, this->m_constVarDiff,
                this->m_isVarDiff, this->m_varD00, this->m_varD01,
                this->m_varD11, df_ptr, this->m_h0, this->m_h1, bwd, deriv0,
                deriv1);

            // Step 4: Apply Laplacian metrics & inner product
            IProduct2DKernel<SHAPE_TYPE, false, true, DEFORMED>(
                nm0, nm1, nq0, nq1, correct, bwd, this->m_dbdata[0],
                this->m_bdata[1], this->m_w[0], this->m_w[1], jac_ptr, wsp0,
                tmpOut);

            IProduct2DKernel<SHAPE_TYPE, false, true, DEFORMED>(
                nm0, nm1, nq0, nq1, correct, deriv0, this->m_bdata[0],
                this->m_dbdata[1], this->m_w[0], this->m_w[1], jac_ptr, wsp0,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nmBlocks;
            outptr += nmBlocks;
        }
    }

#elif defined(SHAPE_DIMENSION_3D)

    // Non-size based operator.
    void operator3D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        const auto nm0 = m_basis[0]->GetNumModes();
        const auto nm1 = m_basis[1]->GetNumModes();
        const auto nm2 = m_basis[2]->GetNumModes();

        const auto nq0 = m_basis[0]->GetNumPoints();
        const auto nq1 = m_basis[1]->GetNumPoints();
        const auto nq2 = m_basis[2]->GetNumPoints();

        constexpr auto ndf  = 9;
        const auto nqTot    = nq0 * nq1 * nq2;
        const auto nmBlocks = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Allocate sufficient workspace for backwards transform and inner
        // product kernels.
        size_t wsp0Size = 0, wsp1Size = 0, wsp2Size = 0;
        BwdTrans3DWorkspace<SHAPE_TYPE>(nm0, nm1, nm2, nq0, nq1, nq2, wsp0Size,
                                        wsp1Size);
        IProduct3DWorkspace<SHAPE_TYPE>(nm0, nm1, nm2, nq0, nq1, nq2, wsp0Size,
                                        wsp1Size, wsp2Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), wsp1(wsp1Size),
            wsp2(wsp2Size);
        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot);
        std::vector<vec_t, allocator<vec_t>> deriv0(nqTot), deriv1(nqTot),
            deriv2(nqTot);

        const vec_t *jac_ptr;
        const vec_t *df_ptr;

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            df_ptr = &((*this->m_df)[e * dfSize]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            if (!DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e * nqTot]);
            }

            // Step 1: BwdTrans
            BwdTrans3DKernel<SHAPE_TYPE>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, tmpIn, this->m_bdata[0],
                this->m_bdata[1], this->m_bdata[2], wsp0, wsp1, bwd);

            // // Step 2: inner product for mass matrix operation
            IProduct3DKernel<SHAPE_TYPE, true, false, DEFORMED>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, bwd, this->m_bdata[0],
                this->m_bdata[1], this->m_bdata[2], this->m_w[0], this->m_w[1],
                this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut, m_lambda);

            // Step 3: take derivatives in standard space
            PhysDerivTensor3DKernel(nq0, nq1, nq2, bwd, this->m_D[0],
                                    this->m_D[1], this->m_D[2], deriv0, deriv1,
                                    deriv2);

            // Step 4: Apply Laplacian metrics & inner product
            Helmholtz3DKernel<SHAPE_TYPE, DEFORMED>(
                nq0, nq1, nq2, this->m_isConstVarDiff, this->m_constVarDiff,
                this->m_isVarDiff, this->m_varD00, this->m_varD01,
                this->m_varD11, this->m_varD02, this->m_varD12, this->m_varD22,
                df_ptr, this->m_h0, this->m_h1, this->m_h2, this->m_h3, deriv0,
                deriv1, deriv2);

            IProduct3DKernel<SHAPE_TYPE, false, true, DEFORMED>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, deriv0,
                this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr, wsp0, wsp1,
                wsp2, tmpOut);

            IProduct3DKernel<SHAPE_TYPE, false, true, DEFORMED>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, deriv1, this->m_bdata[0],
                this->m_dbdata[1], this->m_bdata[2], this->m_w[0], this->m_w[1],
                this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut);

            IProduct3DKernel<SHAPE_TYPE, false, true, DEFORMED>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, deriv2, this->m_bdata[0],
                this->m_bdata[1], this->m_dbdata[2], this->m_w[0], this->m_w[1],
                this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nmBlocks;
            outptr += nmBlocks;
        }
    }

    // Size based template version.
    template <int nm0, int nm1, int nm2, int nq0, int nq1, int nq2>
    void operator3D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        constexpr auto ndf   = 9;
        constexpr auto nqTot = nq0 * nq1 * nq2;
        const auto nmBlocks  = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Allocate sufficient workspace for backwards transform and inner
        // product kernels.
        size_t wsp0Size = 0, wsp1Size = 0, wsp2Size = 0;
        BwdTrans3DWorkspace<SHAPE_TYPE>(nm0, nm1, nm2, nq0, nq1, nq2, wsp0Size,
                                        wsp1Size);
        IProduct3DWorkspace<SHAPE_TYPE>(nm0, nm1, nm2, nq0, nq1, nq2, wsp0Size,
                                        wsp1Size, wsp2Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), wsp1(wsp1Size),
            wsp2(wsp2Size);
        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot);
        std::vector<vec_t, allocator<vec_t>> deriv0(nqTot), deriv1(nqTot),
            deriv2(nqTot);

        const vec_t *jac_ptr;
        const vec_t *df_ptr;

        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            df_ptr = &((*this->m_df)[e * dfSize]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            if (!DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e * nqTot]);
            }

            // Step 1: BwdTrans
            BwdTrans3DKernel<SHAPE_TYPE>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, tmpIn, this->m_bdata[0],
                this->m_bdata[1], this->m_bdata[2], wsp0, wsp1, bwd);

            // // Step 2: inner product for mass matrix operation
            IProduct3DKernel<SHAPE_TYPE, true, false, DEFORMED>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, bwd, this->m_bdata[0],
                this->m_bdata[1], this->m_bdata[2], this->m_w[0], this->m_w[1],
                this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut, m_lambda);

            // Step 3: take derivatives in standard space
            PhysDerivTensor3DKernel(nq0, nq1, nq2, bwd, this->m_D[0],
                                    this->m_D[1], this->m_D[2], deriv0, deriv1,
                                    deriv2);

            // Step 4: Apply Laplacian metrics & inner product
            Helmholtz3DKernel<SHAPE_TYPE, DEFORMED>(
                nq0, nq1, nq2, this->m_isConstVarDiff, this->m_constVarDiff,
                this->m_isVarDiff, this->m_varD00, this->m_varD01,
                this->m_varD11, this->m_varD02, this->m_varD12, this->m_varD22,
                df_ptr, this->m_h0, this->m_h1, this->m_h2, this->m_h3, deriv0,
                deriv1, deriv2);

            IProduct3DKernel<SHAPE_TYPE, false, true, DEFORMED>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, deriv0,
                this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr, wsp0, wsp1,
                wsp2, tmpOut);

            IProduct3DKernel<SHAPE_TYPE, false, true, DEFORMED>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, deriv1, this->m_bdata[0],
                this->m_dbdata[1], this->m_bdata[2], this->m_w[0], this->m_w[1],
                this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut);

            IProduct3DKernel<SHAPE_TYPE, false, true, DEFORMED>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, deriv2, this->m_bdata[0],
                this->m_bdata[1], this->m_dbdata[2], this->m_w[0], this->m_w[1],
                this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nmBlocks;
            outptr += nmBlocks;
        }
    }

#endif // SHAPE_DIMENSION

private:
    int m_nmTot;

    std::vector<vec_t, allocator<vec_t>> m_h0, m_h1, m_h2, m_h3;
};

} // namespace MatrixFree
} // namespace Nektar

#endif
