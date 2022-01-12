#ifndef NEKTAR_LIBRARY_MF_IPRODUCT_H
#define NEKTAR_LIBRARY_MF_IPRODUCT_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "Operator.hpp"
#include "IProductKernels.hpp"

namespace Nektar
{
namespace MatrixFree
{

template<bool DEFORMED = false>
struct IProductSeg : public IProduct, public Helper<1, DEFORMED>
{
    IProductSeg(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProduct(basis, nElmt),
          Helper<1, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdSegData::getNumberOfCoefficients(
                                                 this->m_nm[0]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductSeg<DEFORMED>>(basis, nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out) final
    {
        const int nm0 = m_basis[0]->GetNumModes(); 
        const int nq0 = m_basis[0]->GetNumPoints(); 

        switch(nm0)
        {
        case 2:
            switch(nq0)
            {
            case 2: IProductSegImpl<2 ,2 >(in, out); break;
            case 3: IProductSegImpl<2 ,3 >(in, out); break;
            case 4: IProductSegImpl<2 ,4 >(in, out); break;
            default: IProductSegImpl(nm0, nq0, in, out); break;
            } break;
        case 3:
            switch(nq0)
            {
            case 3: IProductSegImpl<3 ,3 >(in, out); break;
            case 4: IProductSegImpl<3 ,4 >(in, out); break;
            case 5: IProductSegImpl<3 ,5 >(in, out); break;
            case 6: IProductSegImpl<3 ,6 >(in, out); break;
            default: IProductSegImpl(nm0, nq0, in, out); break;
            } break;
        case 4:
            switch(nq0)
            {
            case 4: IProductSegImpl<4 ,4 >(in, out); break;
            case 5: IProductSegImpl<4 ,5 >(in, out); break;
            case 6: IProductSegImpl<4 ,6 >(in, out); break;
            case 7: IProductSegImpl<4 ,7 >(in, out); break;
            case 8: IProductSegImpl<4 ,8 >(in, out); break;
            default: IProductSegImpl(nm0, nq0, in, out); break;
            } break;
        case 5:
            switch(nq0)
            {
            case 5: IProductSegImpl<5 ,5 >(in, out); break;
            case 6: IProductSegImpl<5 ,6 >(in, out); break;
            case 7: IProductSegImpl<5 ,7 >(in, out); break;
            case 8: IProductSegImpl<5 ,8 >(in, out); break;
            case 9: IProductSegImpl<5 ,9 >(in, out); break;
            case 10: IProductSegImpl<5 ,10 >(in, out); break;
            default: IProductSegImpl(nm0, nq0, in, out); break;           
            } break;
        case 6:
            switch(nq0)
            {
            case 6: IProductSegImpl<6 ,6 >(in, out); break;
            case 7: IProductSegImpl<6 ,7 >(in, out); break;
            case 8: IProductSegImpl<6 ,8 >(in, out); break;
            case 9: IProductSegImpl<6 ,9 >(in, out); break;
            case 10: IProductSegImpl<6 ,10 >(in, out); break;
            case 11: IProductSegImpl<6 ,11 >(in, out); break;
            case 12: IProductSegImpl<6 ,12 >(in, out); break;
            default: IProductSegImpl(nm0, nq0, in, out); break;
            } break;
        case 7:
            switch(nq0)
            {
            case 7: IProductSegImpl<7 ,7 >(in, out); break;
            case 8: IProductSegImpl<7 ,8 >(in, out); break;
            case 9: IProductSegImpl<7 ,9 >(in, out); break;
            case 10: IProductSegImpl<7 ,10 >(in, out); break;
            case 11: IProductSegImpl<7 ,11 >(in, out); break;
            case 12: IProductSegImpl<7 ,12 >(in, out); break;
            case 13: IProductSegImpl<7 ,13 >(in, out); break;
            case 14: IProductSegImpl<7 ,14 >(in, out); break;
            default: IProductSegImpl(nm0, nq0, in, out); break;
            } break;
        case 8:
            switch(nq0)
            {
            case 8: IProductSegImpl<8 ,8 >(in, out); break;
            case 9: IProductSegImpl<8 ,9 >(in, out); break;
            case 10: IProductSegImpl<8 ,10 >(in, out); break;
            case 11: IProductSegImpl<8 ,11 >(in, out); break;
            case 12: IProductSegImpl<8 ,12 >(in, out); break;
            case 13: IProductSegImpl<8 ,13 >(in, out); break;
            case 14: IProductSegImpl<8 ,14 >(in, out); break;
            case 15: IProductSegImpl<8 ,15 >(in, out); break;
            case 16: IProductSegImpl<8 ,16 >(in, out); break;
            default: IProductSegImpl(nm0, nq0, in, out); break;
            } break;;
        default: IProductSegImpl(nm0, nq0, in, out); break;
        }
    }

 template<int NM0, int NQ0>
    void IProductSegImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto *inptr = &input[0];
        auto *outptr = &output[0];

        constexpr auto nqTot = NQ0;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[e*nqTot]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductSegKernel(NM0, NQ0, false, false, DEFORMED,
                tmpIn, this->m_bdata[0], this->m_w[0], jac_ptr, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }
    
    void IProductSegImpl(
        const int nm0, const int nq0, 
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto *inptr = &input[0];
        auto *outptr = &output[0];

        const auto nqTot = nq0;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[e*nqTot]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductSegKernel(nm0, nq0, false, false, DEFORMED, 
                tmpIn, this->m_bdata[0], this->m_w[0], jac_ptr, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }
public:

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }
                  
private:
    int m_nmTot;
};

template<bool DEFORMED = false>
struct IProductQuad : public IProduct, public Helper<2, DEFORMED>
{
    IProductQuad(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProduct(basis, nElmt),
          Helper<2, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdQuadData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductQuad<DEFORMED>>(basis, nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out) final
    {
        const int nm0 = m_basis[0]->GetNumModes(); 
        const int nm1 = m_basis[1]->GetNumModes(); 
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints(); 

        if((nm0 == nm1)&&(nq0 == nq1))
        {
            switch(nm0)
            {
            case 2:
                switch(nq0)
                {
                case 2: IProductQuadImpl<2 ,2 ,2 ,2 >(in, out); break;
                case 3: IProductQuadImpl<2 ,2 ,3 ,3 >(in, out); break;
                case 4: IProductQuadImpl<2 ,2 ,4 ,4 >(in, out); break;
                default: IProductQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
                } break;
            case 3:
                switch(nq0)
                {
                case 3: IProductQuadImpl<3 ,3 ,3 ,3 >(in, out); break;
                case 4: IProductQuadImpl<3 ,3 ,4 ,4 >(in, out); break;
                case 5: IProductQuadImpl<3 ,3 ,5 ,5 >(in, out); break;
                case 6: IProductQuadImpl<3 ,3 ,6 ,6 >(in, out); break;
                default: IProductQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
                } break;
            case 4:
                switch(nq0)
                {
                case 4: IProductQuadImpl<4 ,4 ,4 ,4 >(in, out); break;
                case 5: IProductQuadImpl<4 ,4 ,5 ,5 >(in, out); break;
                case 6: IProductQuadImpl<4 ,4 ,6 ,6 >(in, out); break;
                case 7: IProductQuadImpl<4 ,4 ,7 ,7 >(in, out); break;
                case 8: IProductQuadImpl<4 ,4 ,8 ,8 >(in, out); break;
                default: IProductQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
                } break;
            case 5:
                switch(nq0)
                {
                case 5: IProductQuadImpl<5 ,5 ,5 ,5 >(in, out); break;
                case 6: IProductQuadImpl<5 ,5 ,6 ,6 >(in, out); break;
                case 7: IProductQuadImpl<5 ,5 ,7 ,7 >(in, out); break;
                case 8: IProductQuadImpl<5 ,5 ,8 ,8 >(in, out); break;
                case 9: IProductQuadImpl<5 ,5 ,9 ,9 >(in, out); break;
                case 10: IProductQuadImpl<5 ,5 ,10 ,10 >(in, out); break;
                default: IProductQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
                } break;
            case 6:
                switch(nq0)
                {
                case 6: IProductQuadImpl<6 ,6 ,6 ,6 >(in, out); break;
                case 7: IProductQuadImpl<6 ,6 ,7 ,7 >(in, out); break;
                case 8: IProductQuadImpl<6 ,6 ,8 ,8 >(in, out); break;
                case 9: IProductQuadImpl<6 ,6 ,9 ,9 >(in, out); break;
                case 10: IProductQuadImpl<6 ,6 ,10 ,10 >(in, out); break;
                case 11: IProductQuadImpl<6 ,6 ,11 ,11 >(in, out); break;
                case 12: IProductQuadImpl<6 ,6 ,12 ,12 >(in, out); break;
                } break;
            case 7:
                switch(nq0)
                {
                case 7: IProductQuadImpl<7 ,7 ,7 ,7 >(in, out); break;
                case 8: IProductQuadImpl<7 ,7 ,8 ,8 >(in, out); break;
                case 9: IProductQuadImpl<7 ,7 ,9 ,9 >(in, out); break;
                case 10: IProductQuadImpl<7 ,7 ,10 ,10 >(in, out); break;
                case 11: IProductQuadImpl<7 ,7 ,11 ,11 >(in, out); break;
                case 12: IProductQuadImpl<7 ,7 ,12 ,12 >(in, out); break;
                case 13: IProductQuadImpl<7 ,7 ,13 ,13 >(in, out); break;
                case 14: IProductQuadImpl<7 ,7 ,14 ,14 >(in, out); break;
                default: IProductQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
                } break;
            case 8:
                switch(nq0)
                {
                case 8: IProductQuadImpl<8 ,8 ,8 ,8 >(in, out); break;
                case 9: IProductQuadImpl<8 ,8 ,9 ,9 >(in, out); break;
                case 10: IProductQuadImpl<8 ,8 ,10 ,10 >(in, out); break;
                case 11: IProductQuadImpl<8 ,8 ,11 ,11 >(in, out); break;
                case 12: IProductQuadImpl<8 ,8 ,12 ,12 >(in, out); break;
                case 13: IProductQuadImpl<8 ,8 ,13 ,13 >(in, out); break;
                case 14: IProductQuadImpl<8 ,8 ,14 ,14 >(in, out); break;
                case 15: IProductQuadImpl<8 ,8 ,15 ,15 >(in, out); break;
                case 16: IProductQuadImpl<8 ,8 ,16 ,16 >(in, out); break;
                default: IProductQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
                } break;
            default: IProductQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
            }
        }
    }

    template<int NM0, int NM1, int NQ0, int NQ1>
    void IProductQuadImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr = input.data();
        auto* outptr = output.data();

        constexpr auto nqTot = NQ0 * NQ1;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t sums_j[NQ1]; //Sums over eta0 for each value of eta1;

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductQuadKernel(NM0, NM1, NQ0, NQ1, false, false, DEFORMED,
                tmpIn, this->m_bdata[0], this->m_bdata[1],
                this->m_w[0], this->m_w[1], jac_ptr,
                sums_j, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }

    void IProductQuadImpl(
        const int nm0, const int nm1,
        const int nq0, const int nq1, 
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr = input.data();
        auto* outptr = output.data();

        const auto nqTot = nq0 * nq1;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t sums_j[nq1]; //Sums over eta0 for each value of eta1;

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductQuadKernel(nm0, nm1, nq0, nq1, false, false, DEFORMED,
                tmpIn, this->m_bdata[0], this->m_bdata[1],
                this->m_w[0], this->m_w[1], jac_ptr,
                sums_j, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }
public:
    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

private:
    int m_nmTot;
};

template<bool DEFORMED = false>
struct IProductTri : public IProduct, public Helper<2, DEFORMED>
{
    IProductTri(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProduct(basis, nElmt),
          Helper<2, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdTriData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductTri<DEFORMED>>(basis, nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
        Array<OneD, NekDouble> &out) final
    {
        const int nm0 = m_basis[0]->GetNumModes(); 
        const int nm1 = m_basis[1]->GetNumModes(); 
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints(); 

        if((nm0 == nm1)&&(nq0 == nq1+1))
        {
            if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                switch(nm0)
                {
                case 2:
                    switch(nq0)
                    {
                    case 3: IProductTriImpl<2 ,2 ,3 ,2 ,true>(in, out); break;
                    case 4: IProductTriImpl<2 ,2 ,4 ,3 ,true>(in, out); break;
                    default: IProductTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                    } break;
                case 3:
                    switch(nq0)
                    {
                    case 4: IProductTriImpl<3 ,3 ,4 ,3 ,true>(in, out); break;
                    case 5: IProductTriImpl<3 ,3 ,5 ,4 ,true>(in, out); break;
                    case 6: IProductTriImpl<3 ,3 ,6 ,5 ,true>(in, out); break;
                    case 7: IProductTriImpl<3 ,3 ,7 ,6 ,true>(in, out); break;
                    default: IProductTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                    } break;
                case 4:
                    switch(nq0)
                    {
                    case 5: IProductTriImpl<4 ,4 ,5 ,4 ,true>(in, out); break;
                    case 6: IProductTriImpl<4 ,4 ,6 ,5 ,true>(in, out); break;
                    case 7: IProductTriImpl<4 ,4 ,7 ,6 ,true>(in, out); break;
                    case 8: IProductTriImpl<4 ,4 ,8 ,7 ,true>(in, out); break;
                    default: IProductTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                    } break;
                case 5:
                    switch(nq0)
                    {
                    case 6: IProductTriImpl<5 ,5 ,6 ,5 ,true>(in, out); break;
                    case 7: IProductTriImpl<5 ,5 ,7 ,6 ,true>(in, out); break;
                    case 8: IProductTriImpl<5 ,5 ,8 ,7 ,true>(in, out); break;
                    case 9: IProductTriImpl<5 ,5 ,9 ,8 ,true>(in, out); break;
                    case 10: IProductTriImpl<5 ,5 ,10 ,9 ,true>(in, out); break;
                    default: IProductTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                    } break;
                case 6:
                    switch(nq0)
                    {
                    case 7: IProductTriImpl<6 ,6 ,7 ,6 ,true>(in, out); break;
                    case 8: IProductTriImpl<6 ,6 ,8 ,7 ,true>(in, out); break;
                    case 9: IProductTriImpl<6 ,6 ,9 ,8 ,true>(in, out); break;
                    case 10: IProductTriImpl<6 ,6 ,10 ,9 ,true>(in, out); break;
                    case 11: IProductTriImpl<6 ,6 ,11 ,10 ,true>(in, out); break;
                    case 12: IProductTriImpl<6 ,6 ,12 ,11 ,true>(in, out); break;
                    default: IProductTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                    } break;
                case 7:
                    switch(nq0)
                    {
                    case 8: IProductTriImpl<7 ,7 ,8 ,7 ,true>(in, out); break;
                    case 9: IProductTriImpl<7 ,7 ,9 ,8 ,true>(in, out); break;
                    case 10: IProductTriImpl<7 ,7 ,10 ,9 ,true>(in, out); break;
                    case 11: IProductTriImpl<7 ,7 ,11 ,10 ,true>(in, out); break;
                    case 12: IProductTriImpl<7 ,7 ,12 ,11 ,true>(in, out); break;
                    case 13: IProductTriImpl<7 ,7 ,13 ,12 ,true>(in, out); break;
                    case 14: IProductTriImpl<7 ,7 ,14 ,13 ,true>(in, out); break;
                    default: IProductTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                    } break;
                case 8:
                    switch(nq0)
                    {
                    case 9: IProductTriImpl<8 ,8 ,9 ,8 ,true>(in, out); break;
                    case 10: IProductTriImpl<8 ,8 ,10 ,9 ,true>(in, out); break;
                    case 11: IProductTriImpl<8 ,8 ,11 ,10 ,true>(in, out); break;
                    case 12: IProductTriImpl<8 ,8 ,12 ,11 ,true>(in, out); break;
                    case 13: IProductTriImpl<8 ,8 ,13 ,12 ,true>(in, out); break;
                    case 14: IProductTriImpl<8 ,8 ,14 ,13 ,true>(in, out); break;
                    case 15: IProductTriImpl<8 ,8 ,15 ,14 ,true>(in, out); break;
                    case 16: IProductTriImpl<8 ,8 ,16 ,15 ,true>(in, out); break;
                    default: IProductTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                    } break;
                default: IProductTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                }
            }
            else
            {
                switch(nm0)
                {
                case 2:
                    switch(nq0)
                    {
                    case 3: IProductTriImpl<2 ,2 ,3 ,2 ,false>(in, out); break;
                    case 4: IProductTriImpl<2 ,2 ,4 ,3 ,false>(in, out); break;
                    default: IProductTriImpl(nm0,nm1,nq0,nq1,false,in,out); break;
                    } break;
                case 3:
                    switch(nq0)
                    {
                    case 4: IProductTriImpl<3 ,3 ,4 ,3 ,false>(in, out); break;
                    case 5: IProductTriImpl<3 ,3 ,5 ,4 ,false>(in, out); break;
                    case 6: IProductTriImpl<3 ,3 ,6 ,5 ,false>(in, out); break;
                    case 7: IProductTriImpl<3 ,3 ,7 ,6 ,false>(in, out); break;
                    default: IProductTriImpl(nm0,nm1,nq0,nq1,false,in,out); break;
                    } break;
                case 4:
                    switch(nq0)
                    {
                    case 5: IProductTriImpl<4 ,4 ,5 ,4 ,false>(in, out); break;
                    case 6: IProductTriImpl<4 ,4 ,6 ,5 ,false>(in, out); break;
                    case 7: IProductTriImpl<4 ,4 ,7 ,6 ,false>(in, out); break;
                    case 8: IProductTriImpl<4 ,4 ,8 ,7 ,false>(in, out); break;
                    default: IProductTriImpl(nm0,nm1,nq0,nq1,false,in,out); break;
                    } break;
                case 5:
                    switch(nq0)
                    {
                    case 6: IProductTriImpl<5 ,5 ,6 ,5 ,false>(in, out); break;
                    case 7: IProductTriImpl<5 ,5 ,7 ,6 ,false>(in, out); break;
                    case 8: IProductTriImpl<5 ,5 ,8 ,7 ,false>(in, out); break;
                    case 9: IProductTriImpl<5 ,5 ,9 ,8 ,false>(in, out); break;
                    case 10: IProductTriImpl<5 ,5 ,10 ,9 ,false>(in, out); break;
                    default: IProductTriImpl(nm0,nm1,nq0,nq1,false,in,out); break;
                    } break;
                case 6:
                    switch(nq0)
                    {
                    case 7: IProductTriImpl<6 ,6 ,7 ,6 ,false>(in, out); break;
                    case 8: IProductTriImpl<6 ,6 ,8 ,7 ,false>(in, out); break;
                    case 9: IProductTriImpl<6 ,6 ,9 ,8 ,false>(in, out); break;
                    case 10: IProductTriImpl<6 ,6 ,10 ,9 ,false>(in, out); break;
                    case 11: IProductTriImpl<6 ,6 ,11 ,10 ,false>(in, out); break;
                    case 12: IProductTriImpl<6 ,6 ,12 ,11 ,false>(in, out); break;
                    default: IProductTriImpl(nm0,nm1,nq0,nq1,false,in,out); break;
                    } break;
                case 7:
                    switch(nq0)
                    {
                    case 8: IProductTriImpl<7 ,7 ,8 ,7 ,false>(in, out); break;
                    case 9: IProductTriImpl<7 ,7 ,9 ,8 ,false>(in, out); break;
                    case 10: IProductTriImpl<7 ,7 ,10 ,9 ,false>(in, out); break;
                    case 11: IProductTriImpl<7 ,7 ,11 ,10 ,false>(in, out); break;
                    case 12: IProductTriImpl<7 ,7 ,12 ,11 ,false>(in, out); break;
                    case 13: IProductTriImpl<7 ,7 ,13 ,12 ,false>(in, out); break;
                    case 14: IProductTriImpl<7 ,7 ,14 ,13 ,false>(in, out); break;
                    default: IProductTriImpl(nm0,nm1,nq0,nq1,false,in,out); break;
                    } break;
                case 8:
                    switch(nq0)
                    {
                    case 9: IProductTriImpl<8 ,8 ,9 ,8 ,false>(in, out); break;
                    case 10: IProductTriImpl<8 ,8 ,10 ,9 ,false>(in, out); break;
                    case 11: IProductTriImpl<8 ,8 ,11 ,10 ,false>(in, out); break;
                    case 12: IProductTriImpl<8 ,8 ,12 ,11 ,false>(in, out); break;
                    case 13: IProductTriImpl<8 ,8 ,13 ,12 ,false>(in, out); break;
                    case 14: IProductTriImpl<8 ,8 ,14 ,13 ,false>(in, out); break;
                    case 15: IProductTriImpl<8 ,8 ,15 ,14 ,false>(in, out); break;
                    case 16: IProductTriImpl<8 ,8 ,16 ,15 ,false>(in, out); break;
                    default: IProductTriImpl(nm0,nm1,nq0,nq1,false,in,out); break;
                    } break;
                default: IProductTriImpl(nm0,nm1,nq0,nq1,false,in,out); break;
                }
            }
        }
        else
        {
            if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                IProductTriImpl(nm0,nm1,nq0,nq1,true,in,out);
            }
            else
            {
                IProductTriImpl(nm0,nm1,nq0,nq1,false,in,out);
            }
        }
    }

    template<int NM0, int NM1, int NQ0, int NQ1, bool CORRECT>
    void IProductTriImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto *inptr = input.data();
        auto *outptr = output.data();

        constexpr auto nqTot = NQ0 * NQ1;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t eta0_sums[NQ1]; //Sums over eta0 for each value of eta1;

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e =0; e < this->m_nBlocks; ++e)
        {

            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductTriKernel(NM0, NM1, NQ0, NQ1, CORRECT, false,
                              false, DEFORMED,
                tmpIn, this->m_bdata[0], this->m_bdata[1],
                this->m_w[0], this->m_w[1], jac_ptr,
                eta0_sums, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }

    void IProductTriImpl(
        const int nm0, const int nm1,
        const int nq0, const int nq1,
        const bool CORRECT,
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto *inptr = input.data();
        auto *outptr = output.data();

        const auto nqTot = nq0 * nq1;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t eta0_sums[nq1]; //Sums over eta0 for each value of eta1;

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e =0; e < this->m_nBlocks; ++e)
        {

            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductTriKernel(nm0,nm1,nq0,nq1, CORRECT, false,
                              false, DEFORMED,
                tmpIn, this->m_bdata[0], this->m_bdata[1],
                this->m_w[0], this->m_w[1], jac_ptr,
                eta0_sums, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }

public:

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

private:
    int m_nmTot;
};

template<bool DEFORMED = false>
struct IProductHex : public IProduct, public Helper<3, DEFORMED>
{
    IProductHex(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProduct(basis, nElmt),
          Helper<3, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdHexData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductHex<DEFORMED>>(basis, nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out)
    {
        auto nm0 = m_basis[0]->GetNumModes();
        auto nm1 = m_basis[1]->GetNumModes();
        auto nm2 = m_basis[2]->GetNumModes();

        auto nq0 = m_basis[0]->GetNumPoints();
        auto nq1 = m_basis[1]->GetNumPoints();
        auto nq2 = m_basis[2]->GetNumPoints();

        if((nm0 == nm1) && (nm0 == nm2)&&(nq0 == nq1) &&(nq0 == nq2))
        {
            switch(nm0)
            {
            case 2:
                switch(nq0)
                {
                case 2: IProductHexImpl<2, 2, 2, 2, 2, 2>(in, out); break;
                case 3: IProductHexImpl<2, 2, 2, 3, 3, 3>(in, out); break;
                case 4: IProductHexImpl<2, 2, 2, 4, 4, 4>(in, out); break;
                default: IProductHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
                } break;
            case 3:
                switch(nq0)
                {
                case 3: IProductHexImpl<3, 3, 3, 3, 3, 3>(in, out); break;
                case 4: IProductHexImpl<3, 3, 3, 4, 4, 4>(in, out); break;
                case 5: IProductHexImpl<3, 3, 3, 5, 5, 5>(in, out); break;
                case 6: IProductHexImpl<3, 3, 3, 6, 6, 6>(in, out); break;
                default: IProductHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
                } break;
            case 4:
                switch(nq0)
                {
                case 4: IProductHexImpl<4, 4, 4, 4, 4 ,4>(in, out); break;
                case 5: IProductHexImpl<4, 4, 4, 5, 5 ,5>(in, out); break;
                case 6: IProductHexImpl<4, 4, 4, 6, 6 ,6>(in, out); break;
                case 7: IProductHexImpl<4, 4, 4, 7, 7 ,7>(in, out); break;
                case 8: IProductHexImpl<4, 4, 4, 8, 8 ,8>(in, out); break;
                default: IProductHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
                } break;
            case 5:
                switch(nq0)
                {
                case 5: IProductHexImpl<5, 5, 5, 5, 5, 5>(in, out); break;
                case 6: IProductHexImpl<5, 5, 5, 6, 6, 6>(in, out); break;
                case 7: IProductHexImpl<5, 5, 5, 7, 7, 7>(in, out); break;
                case 8: IProductHexImpl<5, 5, 5, 8, 8, 8>(in, out); break;
                case 9: IProductHexImpl<5, 5, 5, 9, 9, 9>(in, out); break;
                case 10: IProductHexImpl<5, 5, 5, 10, 10, 10>(in, out); break;
                default: IProductHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
                } break;
            case 6:
                switch(nq0)
                {
                case 6: IProductHexImpl<6, 6, 6, 6, 6, 6>(in, out); break;
                case 7: IProductHexImpl<6, 6, 6, 7, 7, 7>(in, out); break;
                case 8: IProductHexImpl<6, 6, 6, 8, 8, 8>(in, out); break;
                case 9: IProductHexImpl<6, 6, 6, 9, 9, 9>(in, out); break;
                case 10: IProductHexImpl<6, 6, 6, 10, 10, 10>(in, out); break;
                case 11: IProductHexImpl<6, 6, 6, 11, 11, 11>(in, out); break;
                case 12: IProductHexImpl<6, 6, 6, 12, 12, 12>(in, out); break;
                default: IProductHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
                } break;
            case 7:
                switch(nq0)
                {
                case 7: IProductHexImpl<7, 7, 7, 7, 7, 7>(in, out); break;
                case 8: IProductHexImpl<7, 7, 7, 8, 8, 8>(in, out); break;
                case 9: IProductHexImpl<7, 7, 7, 9, 9, 9>(in, out); break;
                case 10: IProductHexImpl<7, 7, 7, 10, 10, 10>(in, out); break;
                case 11: IProductHexImpl<7, 7, 7, 11, 11, 11>(in, out); break;
                case 12: IProductHexImpl<7, 7, 7, 12, 12, 12>(in, out); break;
                case 13: IProductHexImpl<7, 7, 7, 13, 13, 13>(in, out); break;
                case 14: IProductHexImpl<7, 7, 7, 14, 14, 14>(in, out); break;
                default: IProductHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
                } break;
            case 8:
                switch(nq0)
                {
                case 8: IProductHexImpl<8, 8, 8, 8, 8, 8>(in, out); break;
                case 9: IProductHexImpl<8, 8, 8, 9, 9, 9>(in, out); break;
                case 10: IProductHexImpl<8, 8, 8, 10, 10, 10>(in, out); break;
                case 11: IProductHexImpl<8, 8, 8, 11, 11, 11>(in, out); break;
                case 12: IProductHexImpl<8, 8, 8, 12, 12, 12>(in, out); break;
                case 13: IProductHexImpl<8, 8, 8, 13, 13, 13>(in, out); break;
                case 14: IProductHexImpl<8, 8, 8, 14, 14, 14>(in, out); break;
                case 15: IProductHexImpl<8, 8, 8, 15, 15, 15>(in, out); break;
                case 16: IProductHexImpl<8, 8, 8, 16, 16, 16>(in, out); break;
                default: IProductHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
                } break;
            default: IProductHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
            }
        }
        else
        {
            IProductHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out);
        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2>
    void IProductHexImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using namespace tinysimd;
        using vec_t = simd<NekDouble>;

        auto* inptr = input.data();
        auto* outptr = output.data();

        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t sums_kj[NQ1 * NQ2];
        vec_t sums_k[NQ2];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductHexKernel(NM0, NM1, NM2, NQ0, NQ1, NQ2, false, false,
                              DEFORMED,
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                sums_kj, sums_k,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }

    void IProductHexImpl(
        const int nm0, const int nm1, const int nm2,
        const int nq0, const int nq1, const int nq2,
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        using namespace tinysimd;
        using vec_t = simd<NekDouble>;

        auto* inptr = input.data();
        auto* outptr = output.data();

        const auto nqTot = nq0 * nq1 * nq2;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t sums_kj[nq1 * nq2];
        vec_t sums_k[nq2];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductHexKernel(nm0,nm1,nm2,nq0,nq1,nq2, false, false, DEFORMED,
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                sums_kj, sums_k,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }

public:

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }


private:
    /// Padded basis
    int m_nmTot;
};

template<bool DEFORMED = false>
struct IProductPrism : public IProduct, public Helper<3, DEFORMED>
{
    IProductPrism(std::vector<LibUtilities::BasisSharedPtr> basis,
                     int nElmt)
        : IProduct(basis, nElmt),
          Helper<3, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdPrismData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductPrism<DEFORMED>>(basis, nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
        Array<OneD,       NekDouble> &out) final
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes(); 
        const int nm2 = m_basis[2]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        if((nm0 == nm1)&&(nm0 == nm2)&&(nq0 == nq1)&&(nq0 == nq2 + 1))
        {
            if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                switch(nm0)
                {
                case 2: switch(nq0)
                    {
                    case 3: IProductPrismImpl<2, 2, 2, 3, 3, 2, true>
                            (in, out); break;
                    case 4: IProductPrismImpl<2, 2, 2, 4, 4, 3, true>
                            (in, out); break;
                    default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                               in, out); break;
                    } break;
                case 3:
                    switch(nq0)
                    {
                    case 4: IProductPrismImpl<3, 3, 3, 4, 4, 3, true>
                            (in, out); break;
                    case 5: IProductPrismImpl<3, 3, 3, 5, 5, 4, true>
                            (in, out); break;
                    case 6: IProductPrismImpl<3, 3, 3, 6, 6, 5, true>
                            (in, out); break;
                    default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                               in, out); break;
                    } break;
                case 4:
                    switch(nq0)
                    {
                    case 5: IProductPrismImpl<4, 4, 4, 5, 5, 4, true>
                            (in, out); break;
                    case 6: IProductPrismImpl<4, 4, 4, 6, 6, 5, true>
                            (in, out); break;
                    case 7: IProductPrismImpl<4, 4, 4, 7, 7, 6, true>
                            (in, out); break;
                    case 8: IProductPrismImpl<4, 4, 4, 8, 8, 7, true>
                            (in, out); break;
                    default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                               in, out); break;
                    } break;
                case 5:
                    switch(nq0)
                    {
                    case 6: IProductPrismImpl<5, 5, 5, 6, 6, 5, true>
                            (in, out); break;
                    case 7: IProductPrismImpl<5, 5, 5, 7, 7, 6, true>
                            (in, out); break;
                    case 8: IProductPrismImpl<5, 5, 5, 8, 8, 7, true>
                            (in, out); break;
                    case 9: IProductPrismImpl<5, 5, 5, 9, 9, 8, true>
                            (in, out); break;
                    case 10: IProductPrismImpl<5, 5, 5, 10, 10, 9, true>
                            (in, out); break;
                    default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                               in, out); break;
                    } break;
                case 6:
                    switch(nq0)
                    {
                    case 7: IProductPrismImpl<6, 6, 6, 7, 7, 6, true>
                            (in, out); break;
                    case 8: IProductPrismImpl<6, 6, 6, 8, 8, 7, true>
                            (in, out); break;
                    case 9: IProductPrismImpl<6, 6, 6, 9, 9, 8, true>
                            (in, out); break;
                    case 10: IProductPrismImpl<6, 6, 6, 10, 10, 9, true>
                            (in, out); break;
                    case 11: IProductPrismImpl<6, 6, 6, 11, 11, 10, true>
                            (in, out); break;
                    case 12: IProductPrismImpl<6, 6, 6, 12, 12, 11, true>
                            (in, out); break;
                    default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                               in, out); break;
                    } break;
                case 7:
                    switch(nq0)
                    {
                    case 8: IProductPrismImpl<7, 7, 7, 8, 8, 7, true>
                            (in, out); break;
                    case 9: IProductPrismImpl<7, 7, 7, 9, 9, 8, true>
                            (in, out); break;
                    case 10: IProductPrismImpl<7, 7, 7, 10, 10, 9, true>
                            (in, out); break;
                    case 11: IProductPrismImpl<7, 7, 7, 11, 11, 10, true>
                            (in, out); break;
                    case 12: IProductPrismImpl<7, 7, 7, 12, 12, 11, true>
                            (in, out); break;
                    case 13: IProductPrismImpl<7, 7, 7, 13, 13, 12, true>
                            (in, out); break;
                    case 14: IProductPrismImpl<7, 7, 7, 14, 14, 13, true>
                            (in, out); break;
                    default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                               in, out); break;
                    } break;
                case 8:
                    switch(nq0)
                    {
                    case 9: IProductPrismImpl<8, 8, 8, 9, 9, 8, true>
                            (in, out); break;
                    case 10: IProductPrismImpl<8, 8, 8, 10, 10, 9, true>
                            (in, out); break;
                    case 11: IProductPrismImpl<8, 8, 8, 11, 11, 10, true>
                            (in, out); break;
                    case 12: IProductPrismImpl<8, 8, 8, 12, 12, 11, true>
                            (in, out); break;
                    case 13: IProductPrismImpl<8, 8, 8, 13, 13, 12, true>
                            (in, out); break;
                    case 14: IProductPrismImpl<8, 8, 8, 14, 14, 13, true>
                            (in, out); break;
                    case 15: IProductPrismImpl<8, 8, 8, 15, 15, 14, true>
                            (in, out); break;
                    case 16: IProductPrismImpl<8, 8, 8, 16, 16, 15, true>
                            (in, out); break;
                    default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                               in, out); break;
                    } break;
                default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                           in, out); break;
                }
            }
            else
            {
                switch(nm0)
                {
                case 2:
                    switch(nq0)
                    {
                    case 3: IProductPrismImpl<2, 2, 2, 3, 3, 2, false>
                            (in, out); break;
                    case 4: IProductPrismImpl<2, 2, 2, 4, 4, 3, false>
                            (in, out); break;
                    default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, false,
                                               in, out); break;
                    } break;
                case 3:
                    switch(m_basis[0]->GetNumPoints())
                    {
                    case 4: IProductPrismImpl<3, 3, 3, 4, 4, 3, false>
                            (in, out); break;
                    case 5: IProductPrismImpl<3, 3, 3, 5, 5, 4, false>
                            (in, out); break;
                    case 6: IProductPrismImpl<3, 3, 3, 6, 6, 5, false>
                            (in, out); break;
                    default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, false,
                                               in, out); break;
                    } break;
                case 4:
                    switch(m_basis[0]->GetNumPoints())
                    {
                    case 5: IProductPrismImpl<4, 4, 4, 5, 5, 4, false>
                            (in, out); break;
                    case 6: IProductPrismImpl<4, 4, 4, 6, 6, 5, false>
                            (in, out); break;
                    case 7: IProductPrismImpl<4, 4, 4, 7, 7, 6, false>
                            (in, out); break;
                    case 8: IProductPrismImpl<4, 4, 4, 8, 8, 7, false>
                            (in, out); break;
                    default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, false,
                                               in, out); break;
                    } break;
                case 5:
                    switch(m_basis[0]->GetNumPoints())
                    {
                    case 6: IProductPrismImpl<5, 5, 5, 6, 6, 5, false>
                            (in, out); break;
                    case 7: IProductPrismImpl<5, 5, 5, 7, 7, 6, false>
                            (in, out); break;
                    case 8: IProductPrismImpl<5, 5, 5, 8, 8, 7, false>
                            (in, out); break;
                    case 9: IProductPrismImpl<5, 5, 5, 9, 9, 8, false>
                            (in, out); break;
                    case 10: IProductPrismImpl<5, 5, 5, 10, 10, 9, false>
                            (in, out); break;
                    default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, false,
                                               in, out); break;
                    } break;
                case 6:
                    switch(m_basis[0]->GetNumPoints())
                    {
                    case 7: IProductPrismImpl<6, 6, 6, 7, 7, 6, false>
                            (in, out); break;
                    case 8: IProductPrismImpl<6, 6, 6, 8, 8, 7, false>
                            (in, out); break;
                    case 9: IProductPrismImpl<6, 6, 6, 9, 9, 8, false>
                            (in, out); break;
                    case 10: IProductPrismImpl<6, 6, 6, 10, 10, 9, false>
                            (in, out); break;
                    case 11: IProductPrismImpl<6, 6, 6, 11, 11, 10, false>
                            (in, out); break;
                    case 12: IProductPrismImpl<6, 6, 6, 12, 12, 11, false>
                            (in, out); break;
                    default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, false,
                                               in, out); break;
                    } break;
                case 7:
                    switch(m_basis[0]->GetNumPoints())
                    {
                    case 8: IProductPrismImpl<7, 7, 7, 8, 8, 7, false>
                            (in, out); break;
                    case 9: IProductPrismImpl<7, 7, 7, 9, 9, 8, false>
                            (in, out); break;
                    case 10: IProductPrismImpl<7, 7, 7, 10, 10, 9, false>
                            (in, out); break;
                    case 11: IProductPrismImpl<7, 7, 7, 11, 11, 10, false>
                            (in, out); break;
                    case 12: IProductPrismImpl<7, 7, 7, 12, 12, 11, false>
                            (in, out); break;
                    case 13: IProductPrismImpl<7, 7, 7, 13, 13, 12, false>
                            (in, out); break;
                    case 14: IProductPrismImpl<7, 7, 7, 14, 14, 13, false>
                            (in, out); break;
                    default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, false,
                                               in, out); break;
                    } break;
                case 8:
                    switch(m_basis[0]->GetNumPoints())
                    {
                    case 9: IProductPrismImpl<8, 8, 8, 9, 9, 8, false>
                            (in, out); break;
                    case 10: IProductPrismImpl<8, 8, 8, 10, 10, 9, false>
                            (in, out); break;
                    case 11: IProductPrismImpl<8, 8, 8, 11, 11, 10, false>
                            (in, out); break;
                    case 12: IProductPrismImpl<8, 8, 8, 12, 12, 11, false>
                            (in, out); break;
                    case 13: IProductPrismImpl<8, 8, 8, 13, 13, 12, false>
                            (in, out); break;
                    case 14: IProductPrismImpl<8, 8, 8, 14, 14, 13, false>
                            (in, out); break;
                    case 15: IProductPrismImpl<8, 8, 8, 15, 15, 14, false>
                            (in, out); break;
                    case 16: IProductPrismImpl<8, 8, 8, 16, 16, 15, false>
                            (in, out); break;
                    default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, false,
                                               in, out); break;
                    } break;
                default: IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, false,
                                           in, out); break;
                }
            }
        }
        else
        {
            if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                  in, out);
            }
            else
            {
                IProductPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, false,
                                  in, out);
            }
        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void IProductPrismImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr = input.data();
        auto* outptr = output.data();

        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = NQ0 * NQ1 * NQ2 * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t sums_kj[NQ1 * NQ2];
        vec_t sums_k[NQ2];
        vec_t corr_q[NM1];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductPrismKernel(NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT, false,
                                false, DEFORMED,
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                sums_kj, sums_k,
                corr_q,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }

    }

    void IProductPrismImpl(
        const int nm0, const int nm1, const int nm2,
        const int nq0, const int nq1, const int nq2,
        const bool CORRECT,
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr = input.data();
        auto* outptr = output.data();

        const auto nqTot = nq0 * nq1 * nq2; 
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t sums_kj[nq1 * nq2];
        vec_t sums_k[nq2];
        vec_t corr_q[nm1];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductPrismKernel(nm0,nm1,nm2,nq0,nq1,nq2, CORRECT, false,
                                false, DEFORMED,
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                sums_kj, sums_k,
                corr_q,
                tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }

    }
public:

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

private:
    /// Padded basis
    int m_nmTot;
};

template<bool DEFORMED = false>
struct IProductPyr : public IProduct, public Helper<3, DEFORMED>
{
    IProductPyr(std::vector<LibUtilities::BasisSharedPtr> basis,
                     int nElmt)
        : IProduct(basis, nElmt),
          Helper<3, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdPyrData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductPyr<DEFORMED>>(basis, nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
        Array<OneD,       NekDouble> &out) final
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes(); 
        const int nm2 = m_basis[2]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        if((nm0 == nm1)&&(nm0 == nm2)&&(nq0 == nq1)&&(nq0 == nq2 + 1))
        {
            if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                switch(nm0)
                {
                case 2:
                    switch(nq0)
                    {
                    case 3: IProductPyrImpl<2, 2, 2, 3, 3, 2, true>
                            (in, out); break;
                    case 4: IProductPyrImpl<2, 2, 2, 4, 4, 3, true>
                            (in, out); break;
                    default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                             in, out); break;
                    } break;
                case 3:
                    switch(nq0)
                    {
                    case 4: IProductPyrImpl<3, 3, 3, 4, 4, 3, true>
                            (in, out); break;
                    case 5: IProductPyrImpl<3, 3, 3, 5, 5, 4, true>
                            (in, out); break;
                    case 6: IProductPyrImpl<3, 3, 3, 6, 6, 5, true>
                            (in, out); break;
                    default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                             in, out); break;
                    } break;
                case 4:
                    switch(nq0)
                    {
                    case 5: IProductPyrImpl<4, 4, 4, 5, 5, 4, true>
                            (in, out); break;
                    case 6: IProductPyrImpl<4, 4, 4, 6, 6, 5, true>
                            (in, out); break;
                    case 7: IProductPyrImpl<4, 4, 4, 7, 7, 6, true>
                            (in, out); break;
                    case 8: IProductPyrImpl<4, 4, 4, 8, 8, 7, true>
                            (in, out); break;
                    default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                             in, out); break;
                    } break;
                case 5:
                    switch(nq0)
                    {
                    case 6: IProductPyrImpl<5, 5, 5, 6, 6, 5, true>
                            (in, out); break;
                    case 7: IProductPyrImpl<5, 5, 5, 7, 7, 6, true>
                            (in, out); break;
                    case 8: IProductPyrImpl<5, 5, 5, 8, 8, 7, true>
                            (in, out); break;
                    case 9: IProductPyrImpl<5, 5, 5, 9, 9, 8, true>
                            (in, out); break;
                    case 10: IProductPyrImpl<5, 5, 5, 10, 10, 9, true>
                            (in, out); break;
                    default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                             in, out); break;
                    } break;
                case 6:
                    switch(nq0)
                    {
                    case 7: IProductPyrImpl<6, 6, 6, 7, 7, 6, true>
                            (in, out); break;
                    case 8: IProductPyrImpl<6, 6, 6, 8, 8, 7, true>
                            (in, out); break;
                    case 9: IProductPyrImpl<6, 6, 6, 9, 9, 8, true>
                            (in, out); break;
                    case 10: IProductPyrImpl<6, 6, 6, 10, 10, 9, true>
                            (in, out); break;
                    case 11: IProductPyrImpl<6, 6, 6, 11, 11, 10, true>
                            (in, out); break;
                    case 12: IProductPyrImpl<6, 6, 6, 12, 12, 11, true>
                            (in, out); break;
                    default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                             in, out); break;
                    } break;
                case 7:
                    switch(nq0)
                    {
                    case 8: IProductPyrImpl<7, 7, 7, 8, 8, 7, true>
                            (in, out); break;
                    case 9: IProductPyrImpl<7, 7, 7, 9, 9, 8, true>
                            (in, out); break;
                    case 10: IProductPyrImpl<7, 7, 7, 10, 10, 9, true>
                            (in, out); break;
                    case 11: IProductPyrImpl<7, 7, 7, 11, 11, 10, true>
                            (in, out); break;
                    case 12: IProductPyrImpl<7, 7, 7, 12, 12, 11, true>
                            (in, out); break;
                    case 13: IProductPyrImpl<7, 7, 7, 13, 13, 12, true>
                            (in, out); break;
                    case 14: IProductPyrImpl<7, 7, 7, 14, 14, 13, true>
                            (in, out); break;
                    default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                             in, out); break;
                    } break;
                case 8:
                    switch(nq0)
                    {
                    case 9: IProductPyrImpl<8, 8, 8, 9, 9, 8, true>
                            (in, out); break;
                    case 10: IProductPyrImpl<8, 8, 8, 10, 10, 9, true>
                            (in, out); break;
                    case 11: IProductPyrImpl<8, 8, 8, 11, 11, 10, true>
                            (in, out); break;
                    case 12: IProductPyrImpl<8, 8, 8, 12, 12, 11, true>
                            (in, out); break;
                    case 13: IProductPyrImpl<8, 8, 8, 13, 13, 12, true>
                            (in, out); break;
                    case 14: IProductPyrImpl<8, 8, 8, 14, 14, 13, true>
                            (in, out); break;
                    case 15: IProductPyrImpl<8, 8, 8, 15, 15, 14, true>
                            (in, out); break;
                    case 16: IProductPyrImpl<8, 8, 8, 16, 16, 15, true>
                            (in, out); break;
                    default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                             in, out); break;
                    } break;
                default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                         in, out); break;
                }
            }
            else
            {
                switch(nm0)
                {
                case 2:
                    switch(nq0)
                    {
                    case 3: IProductPyrImpl<2, 2, 2, 3, 3, 2, false>
                            (in, out); break;
                    case 4: IProductPyrImpl<2, 2, 2, 4, 4, 3, false>
                            (in, out); break;
                    default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                             in, out); break;
                    } break;
                case 3:
                    switch(nq0)
                    {
                    case 4: IProductPyrImpl<3, 3, 3, 4, 4, 3, false>
                            (in, out); break;
                    case 5: IProductPyrImpl<3, 3, 3, 5, 5, 4, false>
                            (in, out); break;
                    case 6: IProductPyrImpl<3, 3, 3, 6, 6, 5, false>
                            (in, out); break;
                    default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                             in, out); break;
                    } break;
                case 4:
                    switch(nq0)
                    {
                    case 5: IProductPyrImpl<4, 4, 4, 5, 5, 4, false>
                            (in, out); break;
                    case 6: IProductPyrImpl<4, 4, 4, 6, 6, 5, false>
                            (in, out); break;
                    case 7: IProductPyrImpl<4, 4, 4, 7, 7, 6, false>
                            (in, out); break;
                    case 8: IProductPyrImpl<4, 4, 4, 8, 8, 7, false>
                            (in, out); break;
                    default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                             in, out); break;
                    } break;
                case 5:
                    switch(nq0)
                    {
                    case 6: IProductPyrImpl<5, 5, 5, 6, 6, 5, false>
                            (in, out); break;
                    case 7: IProductPyrImpl<5, 5, 5, 7, 7, 6, false>
                            (in, out); break;
                    case 8: IProductPyrImpl<5, 5, 5, 8, 8, 7, false>
                            (in, out); break;
                    case 9: IProductPyrImpl<5, 5, 5, 9, 9, 8, false>
                            (in, out); break;
                    case 10: IProductPyrImpl<5, 5, 5, 10, 10, 9, false>
                            (in, out); break;
                    default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                             in, out); break;
                    } break;
                case 6:
                    switch(nq0)
                    {
                    case 7: IProductPyrImpl<6, 6, 6, 7, 7, 6, false>
                            (in, out); break;
                    case 8: IProductPyrImpl<6, 6, 6, 8, 8, 7, false>
                            (in, out); break;
                    case 9: IProductPyrImpl<6, 6, 6, 9, 9, 8, false>
                            (in, out); break;
                    case 10: IProductPyrImpl<6, 6, 6, 10, 10, 9, false>
                            (in, out); break;
                    case 11: IProductPyrImpl<6, 6, 6, 11, 11, 10, false>
                            (in, out); break;
                    case 12: IProductPyrImpl<6, 6, 6, 12, 12, 11, false>
                            (in, out); break;
                    default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                             in, out); break;
                    } break;
                case 7:
                    switch(nq0)
                    {
                    case 8: IProductPyrImpl<7, 7, 7, 8, 8, 7, false>
                            (in, out); break;
                    case 9: IProductPyrImpl<7, 7, 7, 9, 9, 8, false>
                            (in, out); break;
                    case 10: IProductPyrImpl<7, 7, 7, 10, 10, 9, false>
                            (in, out); break;
                    case 11: IProductPyrImpl<7, 7, 7, 11, 11, 10, false>
                            (in, out); break;
                    case 12: IProductPyrImpl<7, 7, 7, 12, 12, 11, false>
                            (in, out); break;
                    case 13: IProductPyrImpl<7, 7, 7, 13, 13, 12, false>
                            (in, out); break;
                    case 14: IProductPyrImpl<7, 7, 7, 14, 14, 13, false>
                            (in, out); break;
                    default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                             in, out); break;
                    } break;
                case 8:
                    switch(nq0)
                    {
                    case 9: IProductPyrImpl<8, 8, 8, 9, 9, 8, false>
                            (in, out); break;
                    case 10: IProductPyrImpl<8, 8, 8, 10, 10, 9, false>
                            (in, out); break;
                    case 11: IProductPyrImpl<8, 8, 8, 11, 11, 10, false>
                            (in, out); break;
                    case 12: IProductPyrImpl<8, 8, 8, 12, 12, 11, false>
                            (in, out); break;
                    case 13: IProductPyrImpl<8, 8, 8, 13, 13, 12, false>
                            (in, out); break;
                    case 14: IProductPyrImpl<8, 8, 8, 14, 14, 13, false>
                            (in, out); break;
                    case 15: IProductPyrImpl<8, 8, 8, 15, 15, 14, false>
                            (in, out); break;
                    case 16: IProductPyrImpl<8, 8, 8, 16, 16, 15, false>
                            (in, out); break;
                    default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                             in, out); break;
                    } break;
                default: IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                         in, out); break;
                }
            }
        }
        else
        {
            if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,true, in, out);
            }
            else
            {
                IProductPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,in, out);
            }
        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void IProductPyrImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr  = &input[0];
        auto* outptr = &output[0];

        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = NQ0 * NQ1 * NQ2 * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t sums_kj[NQ1 * NQ2];
        vec_t sums_k[NQ2];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductPyrKernel(NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT, false,
                              false, DEFORMED,
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                sums_kj, sums_k, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }

    }

    void IProductPyrImpl(
        const int nm0, const int nm1, const int nm2,
        const int nq0, const int nq1, const int nq2,
        const bool CORRECT,
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr  = &input[0];
        auto* outptr = &output[0];

        const auto nqTot = nq0 * nq1 * nq2;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t sums_kj[nq1 * nq2];
        vec_t sums_k[nq2];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductPyrKernel(nm0, nm1, nm2, nq0, nq1, nq2, CORRECT, false,
                              false, DEFORMED,
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                sums_kj, sums_k, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }

    }
public:

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

private:
    /// Padded basis
    int m_nmTot;
};



template<bool DEFORMED = false>
struct IProductTet : public IProduct, public Helper<3, DEFORMED>
{
    IProductTet(std::vector<LibUtilities::BasisSharedPtr> basis,
                   int nElmt)
        : IProduct(basis, nElmt),
          Helper<3, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdTetData::getNumberOfCoefficients(
                      this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<IProductTet<DEFORMED>>(basis, nElmt);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
        Array<OneD, NekDouble> &out) final
    {
        const int nm0 = m_basis[0]->GetNumModes();
        const int nm1 = m_basis[1]->GetNumModes(); 
        const int nm2 = m_basis[2]->GetNumModes();
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

        if((nm0 == nm1)&&(nm0 == nm2)&&(nq0 == nq1+1)&&(nq0 == nq2 + 1))
        {
            if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                switch(nm0)
                {
                case 2: switch(nq0)
                    {
                    case 3: IProductTetImpl<2, 2, 2, 3, 2, 2, true>
                            (in, out); break;
                    case 4: IProductTetImpl<2, 2, 2, 4, 3, 3, true>
                            (in, out); break;
                    default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2, true,
                                             in, out); break;
                    } break;
                case 3:
                    switch(nq0)
                    {
                    case 4: IProductTetImpl<3, 3, 3, 4, 3, 3, true>
                            (in, out); break;
                    case 5: IProductTetImpl<3, 3, 3, 5, 4, 4, true>
                            (in, out); break;
                    case 6: IProductTetImpl<3, 3, 3, 6, 5, 5, true>
                            (in, out); break;
                    default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2, true,
                                             in, out); break;
                    } break;
                case 4:
                    switch(nq0)
                    {
                    case 5: IProductTetImpl<4, 4, 4, 5, 4, 4, true>
                            (in, out); break;
                    case 6: IProductTetImpl<4, 4, 4, 6, 5, 5, true>
                            (in, out); break;
                    case 7: IProductTetImpl<4, 4, 4, 7, 6, 6, true>
                            (in, out); break;
                    case 8: IProductTetImpl<4, 4, 4, 8, 7, 7, true>
                            (in, out); break;
                    default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2, true,
                                             in, out); break;
                    } break;
                case 5:
                    switch(nq0)
                    {
                    case 6: IProductTetImpl<5, 5, 5, 6, 5, 5, true>
                            (in, out); break;
                    case 7: IProductTetImpl<5, 5, 5, 7, 6, 6, true>
                            (in, out); break;
                    case 8: IProductTetImpl<5, 5, 5, 8, 7, 7, true>
                            (in, out); break;
                    case 9: IProductTetImpl<5, 5, 5, 9, 8, 8, true>
                            (in, out); break;
                    case 10: IProductTetImpl<5, 5, 5, 10, 9, 9, true>
                            (in, out); break;
                    default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2, true,
                                             in, out); break;
                    } break;
                case 6:
                    switch(nq0)
                    {
                    case 7: IProductTetImpl<6, 6, 6, 7, 6, 6, true>
                        (in, out); break;
                    case 8: IProductTetImpl<6, 6, 6, 8, 7, 7, true>
                            (in, out); break;
                    case 9: IProductTetImpl<6, 6, 6, 9, 8, 8, true>
                            (in, out); break;
                    case 10: IProductTetImpl<6, 6, 6, 10, 9, 9, true>
                            (in, out); break;
                    case 11: IProductTetImpl<6, 6, 6, 11, 10, 10, true>
                            (in, out); break;
                    case 12: IProductTetImpl<6, 6, 6, 12, 11, 11, true>
                            (in, out); break;
                    default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2, true,
                                             in, out); break;
                    } break;
                case 7:
                    switch(nq0)
                    {
                    case 8: IProductTetImpl<7, 7, 7, 8, 7, 7, true>
                            (in, out); break;
                    case 9: IProductTetImpl<7, 7, 7, 9, 8, 8, true>
                            (in, out); break;
                    case 10: IProductTetImpl<7, 7, 7, 10, 9, 9, true>
                            (in, out); break;
                    case 11: IProductTetImpl<7, 7, 7, 11, 10, 10, true>
                            (in, out); break;
                    case 12: IProductTetImpl<7, 7, 7, 12, 11, 11, true>
                        (in, out); break;
                    case 13: IProductTetImpl<7, 7, 7, 13, 12, 12, true>
                            (in, out); break;
                    case 14: IProductTetImpl<7, 7, 7, 14, 13, 13, true>
                            (in, out); break;
                    default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2, true,
                                             in, out); break;
                    } break;
                case 8:
                    switch(nq0)
                    {
                    case 9: IProductTetImpl<8, 8, 8, 9, 8, 8, true>
                            (in, out); break;
                    case 10: IProductTetImpl<8, 8, 8, 10, 9, 9, true>
                            (in, out); break;
                    case 11: IProductTetImpl<8, 8, 8, 11, 10, 10, true>
                            (in, out); break;
                    case 12: IProductTetImpl<8, 8, 8, 12, 11, 11, true>
                            (in, out); break;
                    case 13: IProductTetImpl<8, 8, 8, 13, 12, 12, true>
                            (in, out); break;
                    case 14: IProductTetImpl<8, 8, 8, 14, 13, 13, true>
                            (in, out); break;
                    case 15: IProductTetImpl<8, 8, 8, 15, 14, 14, true>
                            (in, out); break;
                    case 16: IProductTetImpl<8, 8, 8, 16, 15, 15, true>
                            (in, out); break;
                    default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2, true,
                                             in, out); break;
                    } break;
                default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2, true,
                                         in, out); break;
                    
                }
            }
            else
            {
                switch(nm0)
                {
                case 2: switch(nq0)
                    {
                    case 3: IProductTetImpl<2, 2, 2, 3, 2, 2, false>
                            (in, out); break;
                    case 4: IProductTetImpl<2, 2, 2, 4, 3, 3, false>
                            (in, out); break;
                    default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,false,
                                             in, out); break;
                    } break;
                case 3:
                    switch(nq0)
                    {
                    case 4: IProductTetImpl<3, 3, 3, 4, 3, 3, false>
                            (in, out); break;
                    case 5: IProductTetImpl<3, 3, 3, 5, 4, 4, false>
                            (in, out); break;
                    case 6: IProductTetImpl<3, 3, 3, 6, 5, 5, false>
                            (in, out); break;
                    default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,false,
                                             in, out); break;
                } break;
                case 4:
                    switch(nq0)
                    {
                    case 5: IProductTetImpl<4, 4, 4, 5, 4, 4, false>
                            (in, out); break;
                    case 6: IProductTetImpl<4, 4, 4, 6, 5, 5, false>
                            (in, out); break;
                    case 7: IProductTetImpl<4, 4, 4, 7, 6, 6, false>
                            (in, out); break;
                    case 8: IProductTetImpl<4, 4, 4, 8, 7, 7, false>
                            (in, out); break;
                    default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,false,
                                             in, out); break;
                    } break;
                case 5:
                    switch(nq0)
                    {
                    case 6: IProductTetImpl<5, 5, 5, 6, 5, 5, false>
                            (in, out); break;
                    case 7: IProductTetImpl<5, 5, 5, 7, 6, 6, false>
                            (in, out); break;
                    case 8: IProductTetImpl<5, 5, 5, 8, 7, 7, false>
                            (in, out); break;
                    case 9: IProductTetImpl<5, 5, 5, 9, 8, 8, false>
                            (in, out); break;
                    case 10: IProductTetImpl<5, 5, 5, 10, 9, 9, false>
                            (in, out); break;
                    default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,false,
                                             in, out); break;
                    } break;
                case 6:
                    switch(nq0)
                    {
                    case 7: IProductTetImpl<6, 6, 6, 7, 6, 6, false>
                            (in, out); break;
                    case 8: IProductTetImpl<6, 6, 6, 8, 7, 7, false>
                            (in, out); break;
                    case 9: IProductTetImpl<6, 6, 6, 9, 8, 8, false>
                            (in, out); break;
                    case 10: IProductTetImpl<6, 6, 6, 10, 9, 9, false>
                            (in, out); break;
                    case 11: IProductTetImpl<6, 6, 6, 11, 10, 10, false>
                            (in, out); break;
                    case 12: IProductTetImpl<6, 6, 6, 12, 11, 11, false>
                            (in, out); break;
                    default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,false,
                                             in, out); break;
                    } break;
                case 7:
                    switch(nq0)
                    {
                    case 8: IProductTetImpl<7, 7, 7, 8, 7, 7, false>
                            (in, out); break;
                    case 9: IProductTetImpl<7, 7, 7, 9, 8, 8, false>
                            (in, out); break;
                    case 10: IProductTetImpl<7, 7, 7, 10, 9, 9, false>
                            (in, out); break;
                    case 11: IProductTetImpl<7, 7, 7, 11, 10, 10, false>
                            (in, out); break;
                    case 12: IProductTetImpl<7, 7, 7, 12, 11, 11, false>
                            (in, out); break;
                    case 13: IProductTetImpl<7, 7, 7, 13, 12, 12, false>
                            (in, out); break;
                    case 14: IProductTetImpl<7, 7, 7, 14, 13, 13, false>
                            (in, out); break;
                    default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,false,
                                             in, out); break;
                    } break;
                case 8:
                    switch(nq0)
                    {
                    case 9: IProductTetImpl<8, 8, 8, 9, 8, 8, false>
                            (in, out); break;
                    case 10: IProductTetImpl<8, 8, 8, 10, 9, 9, false>
                            (in, out); break;
                    case 11: IProductTetImpl<8, 8, 8, 11, 10, 10, false>
                            (in, out); break;
                    case 12: IProductTetImpl<8, 8, 8, 12, 11, 11, false>
                            (in, out); break;
                    case 13: IProductTetImpl<8, 8, 8, 13, 12, 12, false>
                            (in, out); break;
                    case 14: IProductTetImpl<8, 8, 8, 14, 13, 13, false>
                            (in, out); break;
                    case 15: IProductTetImpl<8, 8, 8, 15, 14, 14, false>
                            (in, out); break;
                    case 16: IProductTetImpl<8, 8, 8, 16, 15, 15, false>
                            (in, out); break;
                    default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,false,
                                             in, out); break;
                    } break;
                default: IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,false,
                                         in, out); break;
                }
            }
        }
        else
        {
            if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,true,in, out);
            }
            else
            {
                IProductTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,false,in, out);
            }
        }
    }
    
    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void IProductTetImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr = input.data();
        auto* outptr = output.data();

        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        constexpr auto nqBlocks = NQ0 * NQ1 * NQ2 * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t wsp[NQ1 * NQ2 + NQ2];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductTetKernel(NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT,
                              false, false, DEFORMED,
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                wsp, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }
        
    void IProductTetImpl(
        const int nm0, const int nm1, const int nm2,
        const int nq0, const int nq1, const int nq2,
        const bool CORRECT, 
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        auto* inptr = input.data();
        auto* outptr = output.data();

        const auto nqTot = nq0 * nq1 * nq2;
        const auto nqBlocks = nqTot * vec_t::width;
        const auto nmBlocks = m_nmTot * vec_t::width;

        vec_t wsp[nq1 * nq2 + nq2];

        std::vector<vec_t, allocator<vec_t>> tmpIn(nqTot), tmpOut(m_nmTot);
        vec_t* jac_ptr;
        for (int e = 0; e < this->m_nBlocks; ++e)
        {
            if (DEFORMED)
            {
                jac_ptr = &((*this->m_jac)[nqTot*e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e]);
            }

            // Load and transpose data
            load_interleave(inptr, nqTot, tmpIn);

            IProductTetKernel(nm0, nm1, nm2, nq0, nq1, nq2, CORRECT,
                              false, false, DEFORMED,
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                wsp, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr += nqBlocks;
            outptr += nmBlocks;
        }
    }
public:

    NekDouble Ndof() final
    {
        return m_nmTot * this->m_nElmt;
    }

private:
    /// Padded basis
    int m_nmTot;
};

} // namespace MatrixFree
} // namespace Nektar

#endif
