#ifndef NEKTAR_LIBRARY_MF_IPRODUCT_H
#define NEKTAR_LIBRARY_MF_IPRODUCT_H

#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "Operator.hpp"

#include "IProductKernels.hpp"

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
struct IProductTemplate
    : public IProduct,
      public Helper<LibUtilities::ShapeTypeDimMap[SHAPE_TYPE], DEFORMED>
{
    IProductTemplate(std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
        : IProduct(basis, nElmt),
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
        }
        else if (DIM == 3)
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients(
                SHAPE_TYPE, this->m_nm[0], this->m_nm[1], this->m_nm[2]);
        }
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
    {
        return std::make_shared<IProductTemplate<SHAPE_TYPE, DEFORMED>>(basis,
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
        const auto nm0 = m_basis[0]->GetNumModes();
        const auto nq0 = m_basis[0]->GetNumPoints();

        const auto nqTot    = nq0;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        // Workspace for kernels - also checks preconditions
        // IProduct1DWorkspace<SHAPE_TYPE>(nm0, nq0);

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);

        vec_t *jac_ptr;

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[e * nqTot]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProduct1DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nq0, tmpIn, this->m_bdata[0], this->m_w[0], jac_ptr,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }

    // Size based template version.
    template <int nm0, int nq0>
    void operator1D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        constexpr auto nqTot    = nq0;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks     = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        // Workspace for kernels - also checks preconditions
        // IProduct1DWorkspace<SHAPE_TYPE>(nm0, nq0);

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);

        vec_t *jac_ptr;

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[e * nqTot]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProduct1DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nq0, tmpIn, this->m_bdata[0], this->m_w[0], jac_ptr,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
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

        const auto nqTot    = nq0 * nq1;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0;
        IProduct2DWorkspace<SHAPE_TYPE>(nm0, nm1, nq0, nq1, wsp0Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), tmpIn(nqTot),
            tmpOut(m_nmTot);

        vec_t *jac_ptr;

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot * e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProduct2DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nm1, nq0, nq1, correct, tmpIn, this->m_bdata[0],
                this->m_bdata[1], this->m_w[0], this->m_w[1], jac_ptr, wsp0,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }

    // Size based template version.
    template <int nm0, int nm1, int nq0, int nq1>
    void operator2D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        constexpr auto nqTot    = nq0 * nq1;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks     = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0;
        IProduct2DWorkspace<SHAPE_TYPE>(nm0, nm1, nq0, nq1, wsp0Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), tmpIn(nqTot),
            tmpOut(m_nmTot);

        vec_t *jac_ptr;

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot * e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProduct2DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nm1, nq0, nq1, correct, tmpIn, this->m_bdata[0],
                this->m_bdata[1], this->m_w[0], this->m_w[1], jac_ptr, wsp0,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
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

        const auto nqTot    = nq0 * nq1 * nq2;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0, wsp1Size = 0, wsp2Size = 0;
        IProduct3DWorkspace<SHAPE_TYPE>(nm0, nm1, nm2, nq0, nq1, nq2, wsp0Size,
                                        wsp1Size, wsp2Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), wsp1(wsp1Size),
            wsp2(wsp2Size), tmpIn(nqTot), tmpOut(m_nmTot);

        vec_t *jac_ptr;

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot * e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProduct3DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, tmpIn, this->m_bdata[0],
                this->m_bdata[1], this->m_bdata[2], this->m_w[0], this->m_w[1],
                this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }

    // Size based template version.
    template <int nm0, int nm1, int nm2, int nq0, int nq1, int nq2>
    void operator3D(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output)
    {
        constexpr auto nqTot    = nq0 * nq1 * nq2;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks     = m_nmTot * vec_t::width;

        auto *inptr  = &input[0];
        auto *outptr = &output[0];

        const bool correct =
            (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0, wsp1Size = 0, wsp2Size = 0;
        IProduct3DWorkspace<SHAPE_TYPE>(nm0, nm1, nm2, nq0, nq1, nq2, wsp0Size,
                                        wsp1Size, wsp2Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), wsp1(wsp1Size),
            wsp2(wsp2Size), tmpIn(nqTot), tmpOut(m_nmTot);

        vec_t *jac_ptr;

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot * e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProduct3DKernel<SHAPE_TYPE, false, false, DEFORMED>(
                nm0, nm1, nm2, nq0, nq1, nq2, correct, tmpIn, this->m_bdata[0],
                this->m_bdata[1], this->m_bdata[2], this->m_w[0], this->m_w[1],
                this->m_w[2], jac_ptr, wsp0, wsp1, wsp2, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }

#endif // SHAPE_DIMENSION

private:
    int m_nmTot;
};

} // namespace MatrixFree
} // namespace Nektar

#endif
