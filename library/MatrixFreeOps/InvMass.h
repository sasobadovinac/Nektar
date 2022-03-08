#ifndef NEKTAR_LIBRARY_MF_INVMASS_H
#define NEKTAR_LIBRARY_MF_INVMASS_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

#include "Operator.hpp"
#include "BwdTransKernels.hpp"

namespace Nektar
{
namespace MatrixFree
{

template<bool DEFORMED = false>
struct InvMassQuad : public InvMass, public Helper<2>
{
private:
    std::vector<vec_t, tinysimd::allocator<vec_t>> m_bdata_inv;
    std::vector<vec_t, tinysimd::allocator<vec_t>> m_bdata_inv_trans;
    std::vector<vec_t, tinysimd::allocator<vec_t>> m_w_inv;

public:
    InvMassQuad(std::vector<LibUtilities::BasisSharedPtr> basis,
                    int nElmt)
        : InvMass(basis, nElmt),
          Helper<2>(basis, nElmt),
          m_nmTot(LibUtilities::StdQuadData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1]))
    {
        ASSERTL0(this->m_nm[0] == this->m_nm[1], "must be equal");

        const int nm = this->m_nm[0];

        // Create a matrix consisting of _inverse_ basis data from our existing
        // basis. For this we need to also adjust number of quadrature points.
        LibUtilities::PointsKey newpkey(nm, basis[0]->GetPointsType());
        LibUtilities::BasisKey  newbkey(basis[0]->GetBasisType(),
                                        nm, newpkey);

        auto newbasis = LibUtilities::BasisManager()[newbkey];

        NekMatrix<NekDouble> bmat(nm, nm);
        for (int i = 0; i < nm; ++i)
        {
            for (int j = 0; j < nm; ++j)
            {
                bmat(i,j) = newbasis->GetBdata()[i + j*nm];
            }
        }

        for (int i = 0; i < nm; ++i)
        {
            for (int j = 0; j < nm; ++j)
            {
                std::cout << bmat(i, j) << " ";
            }
            std::cout << std::endl;
        }

        bmat.Invert();

        for (int i = 0; i < nm; ++i)
        {
            for (int j = 0; j < nm; ++j)
            {
                std::cout << bmat(i, j) << " ";
            }
            std::cout << std::endl;
        }

        m_w_inv.resize(nm);
        m_bdata_inv.resize(nm * nm);
        m_bdata_inv_trans.resize(nm * nm);

        for (int i = 0; i < nm * nm; ++i)
        {
            m_bdata_inv[i] = bmat.GetPtr()[i];
        }

        for (int i = 0; i < nm; ++i)
        {
            m_w_inv[i] = 1.0 / newbasis->GetW()[i];

            for (int j = 0; j < nm; ++j)
            {
                NekDouble hmm = bmat(i,j);
                m_bdata_inv_trans[i*nm+j] = hmm;
            }
        }
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<InvMassQuad>(basis, nElmt);
    }

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

    void operator()(const Array<OneD, const NekDouble> &in,
        Array<OneD,       NekDouble> &out) final
    {
        // Check preconditions
        ASSERTL0(m_basis[0]->GetNumModes() == m_basis[1]->GetNumModes() &&
            m_basis[0]->GetNumPoints() == m_basis[1]->GetNumPoints(),
            "MatrixFree requires homogenous modes/points");

        switch(m_basis[0]->GetNumModes())
        {
            case 2: InvMassQuadImpl<2 ,2 >(in, out); break;
            case 3: InvMassQuadImpl<3 ,3 >(in, out); break;
            case 4: InvMassQuadImpl<4 ,4 >(in, out); break;
            case 5: InvMassQuadImpl<5 ,5 >(in, out); break;
            case 6: InvMassQuadImpl<6 ,6 >(in, out); break;
            case 7: InvMassQuadImpl<7 ,7 >(in, out); break;
            case 8: InvMassQuadImpl<8 ,8 >(in, out); break;
            case 9: InvMassQuadImpl<9 ,9 >(in, out); break;
            default: NEKERROR(ErrorUtil::efatal,
                "InvMassQuad: # of modes / points combo not yet initialised.");
        }
    }

    template<int NM0, int NM1>
    void InvMassQuadImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using namespace tinysimd;
        using vec_t = simd<NekDouble>;

        auto *inptr = &input[0];
        auto *outptr = &output[0];

        const auto nmBlocks = m_nmTot * vec_t::width;
        constexpr auto nqTot = NM0 * NM1;

        vec_t p_sums[NM0 * NM1]; //Sums over q for each quadpt p_i
        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(m_nmTot);

        const vec_t* jac_ptr;

        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            // TODO: DEFORMED would be broken here since m_jac is stored at quad
            // pts of original element. Could be fixed via an interpolation
            jac_ptr = DEFORMED ?
                &((*this->m_jac)[e * nqTot]) : &((*this->m_jac)[e]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            // Backwards transform with inv_bdata_trans
            BwdTransQuadKernel<NM0, NM1, NM0, NM1>(
                tmpIn, this->m_bdata_inv_trans, this->m_bdata_inv_trans,
                p_sums, tmpOut);

            // Divide by jac & quadrature weight term
            for (int i = 0; i < NM0; ++i)
            {
                for (int j = 0; j < NM0; ++j)
                {
                    tmpIn[i * NM0 + j] = m_w_inv[j] * tmpOut[i * NM0 + j];
                }
            }

            for (int i = 0; i < NM0; ++i)
            {
                for (int j = 0; j < NM0; ++j)
                {
                    tmpIn[i + j * NM0] = m_w_inv[j] * tmpIn[i + j * NM0];
                }
            }

            for (int i = 0; i < NM0 * NM0; ++i)
            {
                tmpIn[i] = tmpIn[i] / jac_ptr[0];
            }

            // Backwards transform with inv_bdata
            BwdTransQuadKernel<NM0, NM1, NM0, NM1>(
                tmpIn, this->m_bdata_inv, this->m_bdata_inv,
                p_sums, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, NM0*NM1, outptr);

            inptr += nmBlocks;
            outptr += nmBlocks;
        }
    }

private:
    int m_nmTot;
};

}
}

#endif
