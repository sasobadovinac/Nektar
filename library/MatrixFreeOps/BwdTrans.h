#ifndef NEKTAR_LIBRARY_MF_BWDTRANS_H
#define NEKTAR_LIBRARY_MF_BWDTRANS_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "Operator.hpp"
#include "BwdTransKernels.hpp"


// As each opertor has seven shapes over three dimension to get to the
// "work" each operator uses a series of preprocessor directives based
// on the dimension and shape type so to limit the code
// inclusion. This keeps the library size as minimal as possible while
// removing runtime conditionals.

// The preprocessor directives, SHAPE_DIMENSION_?D and SHAPE_TYPE_*
// are constructed by CMake in the CMakeLists.txt file which uses an
// implementation file.  See the CMakeLists.txt files for more
// details.
namespace Nektar
{
namespace MatrixFree
{

template<LibUtilities::ShapeType SHAPE_TYPE, bool DEFORMED = false>
struct BwdTransTemplate : public BwdTrans,
                          public Helper<LibUtilities::ShapeTypeDimMap[SHAPE_TYPE]>
{
    BwdTransTemplate(std::vector<LibUtilities::BasisSharedPtr> basis,
                  int nElmt)
    : BwdTrans(basis, nElmt),
      Helper<LibUtilities::ShapeTypeDimMap[SHAPE_TYPE]>(basis, nElmt)
    {
        constexpr auto DIM = LibUtilities::ShapeTypeDimMap[SHAPE_TYPE];

        if( DIM == 1 )
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients( SHAPE_TYPE, this->m_nm[0] );
        }
        else if( DIM == 2 )
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients( SHAPE_TYPE, this->m_nm[0], this->m_nm[1] );
        }
        else if( DIM == 3 )
        {
            m_nmTot = LibUtilities::GetNumberOfCoefficients( SHAPE_TYPE, this->m_nm[0], this->m_nm[1], this->m_nm[2] );
        }
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
      return std::make_shared<BwdTransTemplate<SHAPE_TYPE>>(basis, nElmt);
    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

    void operator()(const Array<OneD, const NekDouble>& input,
                          Array<OneD,       NekDouble>& output) final
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
    void operator1D(const Array<OneD, const NekDouble>& input,
                          Array<OneD,       NekDouble>& output)
    {
        const auto nm0 = m_basis[0]->GetNumModes();
        const auto nq0 = m_basis[0]->GetNumPoints();

        const auto nqTot = nq0;
        const auto nqBlocks =   nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        auto* inptr  = &input[0];
        auto* outptr = &output[0];

        // Workspace for kernels - also checks preconditions
        BwdTrans1DWorkspace<SHAPE_TYPE>(nm0, nq0);

        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(nqTot);

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            BwdTrans1DKernel<SHAPE_TYPE>(nm0, nq0,
                                         tmpIn, this->m_bdata[0], tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nqTot, outptr);

            inptr  += nmBlocks;
            outptr += nqBlocks;
        }
    }

    // Size based template version.
    template<int nm0, int nq0>
    void operator1D(const Array<OneD, const NekDouble>& input,
                          Array<OneD,       NekDouble>& output)
    {
        constexpr auto nqTot = nq0;
        constexpr auto nqBlocks =   nqTot * vec_t::width;
        const     auto nmBlocks = m_nmTot * vec_t::width;

        auto* inptr  = &input[0];
        auto* outptr = &output[0];

        // Workspace for kernels - also checks preconditions
        BwdTrans1DWorkspace<SHAPE_TYPE>(nm0, nq0);

        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(nqTot);

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            BwdTrans1DKernel<SHAPE_TYPE>(nm0, nq0,
                                         tmpIn, this->m_bdata[0], tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nqTot, outptr);

            inptr  += nmBlocks;
            outptr += nqBlocks;
        }
    }

#elif defined(SHAPE_DIMENSION_2D)

    // Non-size based operator.
    void operator2D(const Array<OneD, const NekDouble>& input,
                          Array<OneD,       NekDouble>& output)
    {
        const auto nm0 = m_basis[0]->GetNumModes();
        const auto nm1 = m_basis[1]->GetNumModes();

        const auto nq0 = m_basis[0]->GetNumPoints();
        const auto nq1 = m_basis[1]->GetNumPoints();

        const auto nqTot = nq0 * nq1;
        const auto nqBlocks =   nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;


        auto* inptr  = &input[0];
        auto* outptr = &output[0];
        const bool correct = (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0;
        BwdTrans2DWorkspace<SHAPE_TYPE>(nm0, nm1, nq0, nq1, wsp0Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), tmpIn(m_nmTot), tmpOut(nqTot);

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            BwdTrans2DKernel<SHAPE_TYPE>(nm0, nm1, nq0, nq1, correct,
                                         tmpIn, this->m_bdata[0], this->m_bdata[1],
                                         wsp0, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nqTot, outptr);

            inptr  += nmBlocks;
            outptr += nqBlocks;
        }
    }

    // Size based template version.
    template<int nm0, int nm1, int nq0, int nq1>
    void operator2D(const Array<OneD, const NekDouble>& input,
                          Array<OneD,       NekDouble>& output)
    {
        constexpr auto nqTot = nq0 * nq1;
        constexpr auto nqBlocks =   nqTot * vec_t::width;
        const     auto nmBlocks = m_nmTot * vec_t::width;

        auto* inptr  = &input[0];
        auto* outptr = &output[0];
        const bool correct = (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0;
        BwdTrans2DWorkspace<SHAPE_TYPE>(nm0, nm1, nq0, nq1, wsp0Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), tmpIn(m_nmTot), tmpOut(nqTot);

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            BwdTrans2DKernel<SHAPE_TYPE>(nm0, nm1, nq0, nq1, correct,
                                         tmpIn, this->m_bdata[0], this->m_bdata[1],
                                         wsp0, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nqTot, outptr);

            inptr  += nmBlocks;
            outptr += nqBlocks;
        }
    }

#elif defined(SHAPE_DIMENSION_3D)

    // Non-size based operator.
    void operator3D(const Array<OneD, const NekDouble> &input,
                          Array<OneD,       NekDouble> &output)
    {
        const auto nm0 = m_basis[0]->GetNumModes();
        const auto nm1 = m_basis[1]->GetNumModes();
        const auto nm2 = m_basis[2]->GetNumModes();

        const auto nq0 = m_basis[0]->GetNumPoints();
        const auto nq1 = m_basis[1]->GetNumPoints();
        const auto nq2 = m_basis[2]->GetNumPoints();

        const auto nqTot = nq0 * nq1 * nq2;
        const auto nqBlocks =   nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        auto* inptr  = &input[0];
        auto* outptr = &output[0];

        const bool correct = (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0, wsp1Size = 0;
        BwdTrans3DWorkspace<SHAPE_TYPE>(nm0, nm1, nm2, nq0, nq1, nq2,
                            wsp0Size, wsp1Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), wsp1(wsp1Size),
            tmpIn(m_nmTot), tmpOut(nqTot);

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            BwdTrans3DKernel<SHAPE_TYPE>(nm0, nm1, nm2, nq0, nq1, nq2, correct,
                             tmpIn,
                             this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                             wsp0, wsp1, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nqTot, outptr);

            inptr  += nmBlocks;
            outptr += nqBlocks;
        }
    }

    // Size based template version.
    template<int nm0, int nm1, int nm2, int nq0, int nq1, int nq2>
    void operator3D(const Array<OneD, const NekDouble> &input,
                          Array<OneD,       NekDouble> &output)
    {
        constexpr auto nqTot = nq0 * nq1 * nq2;
        constexpr auto nqBlocks =   nqTot * vec_t::width;
        const     auto nmBlocks = m_nmTot * vec_t::width;

        auto* inptr  = &input[0];
        auto* outptr = &output[0];
        const bool correct = (m_basis[0]->GetBasisType() == LibUtilities::eModified_A);

        // Workspace for kernels - also checks preconditions
        size_t wsp0Size = 0, wsp1Size = 0;
        BwdTrans3DWorkspace<SHAPE_TYPE>(nm0, nm1, nm2, nq0, nq1, nq2,
                            wsp0Size, wsp1Size);

        std::vector<vec_t, allocator<vec_t>> wsp0(wsp0Size), wsp1(wsp1Size),
            tmpIn(m_nmTot), tmpOut(nqTot);

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            BwdTrans3DKernel<SHAPE_TYPE>(nm0, nm1, nm2, nq0, nq1, nq2, correct,
                             tmpIn,
                             this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                             wsp0, wsp1, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, nqTot, outptr);

            inptr  += nmBlocks;
            outptr += nqBlocks;
        }
    }

#endif // SHAPE_DIMENSION

private:
    int m_nmTot;
};

} // namespace MatrixFree
} // namespace Nektar

#endif
