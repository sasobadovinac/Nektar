#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include "MatrixFreeDeclspec.h"

#include <iostream>

#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <Collections/Operator.h>
// #include "AVXAssembly.h"

#include <LibUtilities/SimdLib/tinysimd.hpp>

namespace Nektar
{
namespace MatrixFree
{

using vec_t = tinysimd::simd<NekDouble>;

/// Operator base class
class Operator
{
public:
    virtual ~Operator()
    {
    }

    /// Number of gigflops required to compute the operator.
    virtual NekDouble GFlops()
    {
        return 0.0;
    }

    /// Number of degrees of freedom that this operator will process.
    virtual NekDouble Ndof()
    {
        return 0.0;
    }

    /// Number of stores to memory
    virtual NekDouble NStores()
    {
        return 0.0;
    }

    /// Number of loads from memory
    virtual NekDouble NLoads()
    {
        return 0.0;
    }

    /// This operator requires derivative factors.
    MATRIXFREE_EXPORT virtual bool NeedsDF()
    {
        return false;
    }

    /// This operator requires Jacobian.
    MATRIXFREE_EXPORT virtual bool NeedsJac()
    {
        return false;
    }

    MATRIXFREE_EXPORT virtual void SetJac(const Array<OneD, const NekDouble> &jac) = 0;

    MATRIXFREE_EXPORT virtual void SetDF(const Array<TwoD, const NekDouble> &df) = 0;

    // virtual void set_asmMap(MultiRegions::AssemblyMapSharedPtr asmMap)
    // {
    // }

    // /// Provides a reference function for operators which can be used to
    // /// validate correctness compared to Nektar++.
    // ///
    // /// In 2D this is taken as \f$ \sin x \cos y \f$, and in 3D \f$ \sin x \cos
    // /// y \sin z \f$.
    // static void RefFn(MultiRegions::ExpListSharedPtr expList,
    //                   Array<OneD, NekDouble> &in)
    // {
    //     const int nq = expList->GetNpoints();
    //     const int dim = expList->GetExp(0)->GetShapeDimension();

    //     if(dim == 2){
    //         Array<OneD, NekDouble> xc(nq), yc(nq);

    //         expList->GetCoords(xc,yc);
    //         for(int i = 0; i < nq; i++){
    //             in[i] = sin(xc[i]) * cos(yc[i]);
    //         }
    //     }
    //     else if(dim == 3){
    //         Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

    //         expList->GetCoords(xc,yc,zc);
    //         for(int i = 0; i < nq; i++){
    //             in[i] = sin(xc[i]) * cos(yc[i]) * sin(zc[i]);
    //         }
    //     }
    // }

};

typedef std::shared_ptr<Operator> OperatorSharedPtr;

/// Base class for backwards transform operator.
class BwdTrans : virtual public Operator
{
public:
    BwdTrans(std::vector<LibUtilities::BasisSharedPtr> basis,
             int nElmt) :
        m_basis(basis), m_nElmt(nElmt)
    {
    }

    virtual ~BwdTrans()
    {
    }

    MATRIXFREE_EXPORT virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output) = 0; //Abstract Method

    // void Ref(MultiRegions::ExpListSharedPtr expList,
    //          Array<OneD, NekDouble> &ref_exp,
    //          Array<OneD, NekDouble> &ref_bwd)
    // {
    //     Array<OneD, NekDouble> ref_fn(expList->GetNpoints());
    //     this->RefFn(expList, ref_fn);
    //     expList->FwdTrans_IterPerExp(ref_fn, ref_exp);
    //     expList->BwdTrans(ref_exp, ref_bwd);
    // }

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
};

class IProduct : virtual public Operator
{
public:
    IProduct(std::vector<LibUtilities::BasisSharedPtr> basis,
             int nElmt) :
        m_basis(basis), m_nElmt(nElmt)
    {
    }

    virtual ~IProduct()
    {
    }

    bool NeedsJac() final
    {
        return true;
    }

    MATRIXFREE_EXPORT virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output) = 0;

protected:
    /// Vector of tensor product basis directions
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
};

class PhysDeriv : virtual public Operator
{
public:
    PhysDeriv(std::vector<LibUtilities::BasisSharedPtr> basis,
             int nElmt) :
        m_basis(basis), m_nElmt(nElmt)
    {
    }

    virtual ~PhysDeriv()
    {
    }

    bool NeedsDF() final
    {
        return true;
    }

    MATRIXFREE_EXPORT virtual void operator()(const Array<OneD, const NekDouble> &in,
                                Array<OneD,       NekDouble> &out_d0,
                                Array<OneD,       NekDouble> &out_d1) = 0;

    MATRIXFREE_EXPORT virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output0,
        Array<OneD, NekDouble> &output1,
        Array<OneD, NekDouble> &output2) = 0;

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
};

class IProductWRTDerivBase : virtual public Operator
{
public:
    IProductWRTDerivBase(std::vector<LibUtilities::BasisSharedPtr> basis,
              int nElmt) :
        m_basis(basis), m_nElmt(nElmt)
    {
    }

    virtual ~IProductWRTDerivBase()
    {
    }

    MATRIXFREE_EXPORT bool NeedsJac() final
    {
        return true;
    }

    MATRIXFREE_EXPORT bool NeedsDF() final
    {
        return true;
    }

    MATRIXFREE_EXPORT virtual void operator()(
        const Array<OneD, Array<OneD, NekDouble>> &in,
        Array<OneD, NekDouble> &out) = 0;

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
};

class Helmholtz : virtual public Operator
{
public:
    Helmholtz(std::vector<LibUtilities::BasisSharedPtr> basis,
              int nElmt) :
        m_basis(basis), m_nElmt(nElmt)
    {
    }

    virtual ~Helmholtz()
    {
    }

    bool NeedsJac() final
    {
        return true;
    }

    bool NeedsDF() final
    {
        return true;
    }

    MATRIXFREE_EXPORT virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output) = 0;

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
};

// template<int VW>
// class HelmholtzGlobal : virtual public Operator
// {
//     public:
//     HelmholtzGlobal(std::vector<LibUtilities::BasisSharedPtr> basis,
//             int nElmt) :
//             m_basis(basis), m_nElmt(nElmt)
//     {
//     }

//     virtual ~HelmholtzGlobal()
//     {
//     }

//     virtual bool NeedsJac() override
//     {
//         return true;
//     }

//     virtual bool NeedsDF() override
//     {
//         return true;
//     }

//     virtual void operator()(
//         const Array<OneD, const NekDouble> &globalIn,
//         Array<OneD,       NekDouble> &globalOut,
//         AVXAssembly<VW>              &l2g) = 0;


//     // void Ref(MultiRegions::ExpListSharedPtr expList,
//     //          Array<OneD, NekDouble> &ref_exp,
//     //          Array<OneD, NekDouble> &ref_helm)
//     // {
//     //     Array<OneD, NekDouble> ref_fn(expList->GetNpoints());
//     //     this->RefFn(expList, ref_fn);
//     //     expList->FwdTrans_IterPerExp(ref_fn, ref_exp);

//     //     StdRegions::ConstFactorMap factors;
//     //     factors[StdRegions::eFactorLambda] = 1.0;

//     //     MultiRegions::GlobalMatrixKey mkey(
//     //         StdRegions::eHelmholtz,
//     //         MultiRegions::NullAssemblyMapSharedPtr,
//     //         factors);

//     //     expList->GeneralMatrixOp_IterPerExp(mkey, ref_exp, ref_helm);
//     // }

//     // void Ref2(MultiRegions::ExpListSharedPtr expList,
//     //         Array<OneD, NekDouble> &ref_exp,
//     //         Array<OneD, NekDouble> &ref_helm,
//     //         const MultiRegions::AssemblyMapSharedPtr asmMap)
//     // {
//     //     Array<OneD, NekDouble> ref_fn(expList->GetNpoints());
//     //     Array<OneD, NekDouble> ref_exp_local (expList->GetNcoeffs());
//     //     this->RefFn(expList, ref_fn);
//     //     expList->FwdTrans_IterPerExp(ref_fn, ref_exp_local);

//     //     StdRegions::ConstFactorMap factors;
//     //     factors[StdRegions::eFactorLambda] = 1.0;

//     //     MultiRegions::GlobalMatrixKey mkey(
//     //         StdRegions::eHelmholtz,
//     //         MultiRegions::NullAssemblyMapSharedPtr,
//     //         factors);

//     //     Array<OneD, NekDouble> ref_helm_local(expList->GetNcoeffs());
//     //     expList->GeneralMatrixOp_IterPerExp(mkey, ref_exp_local, ref_helm_local);

//     //     asmMap->Assemble(ref_helm_local, ref_helm);
//     // }

// protected:
//     std::vector<LibUtilities::BasisSharedPtr> m_basis;
//     int m_nElmt;

// };

template <int DIM, bool DEFORMED = false>
class Helper : virtual public Operator
{
protected:
    Helper(std::vector<LibUtilities::BasisSharedPtr> basis,
              int nElmt)
        : Operator()
    {
        // Sanity check: no padding yet!
        ASSERTL0(nElmt % vec_t::width == 0, "Number of elements not divisible by vector "
                 "width, padding not yet implemented.");

        // Calculate number of 'blocks', i.e. meta-elements
        m_nBlocks = nElmt / vec_t::width;

        // Depending on element dimension, set up basis information, quadrature,
        // etc, inside vectorised environment.
        for (int i = 0; i < DIM; ++i)
        {
            const Array<OneD, const NekDouble> bdata = basis[i]->GetBdata();
            const Array<OneD, const NekDouble> dbdata = basis[i]->GetDbdata();
            const Array<OneD, const NekDouble> w = basis[i]->GetW();

            m_nm[i] = basis[i]->GetNumModes();
            m_nq[i] = basis[i]->GetNumPoints();

            m_bdata[i].resize(bdata.size());
            for (auto j = 0; j < bdata.size(); ++j)
            {
                m_bdata[i][j] = bdata[j];
            }

            m_dbdata[i].resize(dbdata.size());
            for (auto j = 0; j < dbdata.size(); ++j)
            {
                m_dbdata[i][j] = dbdata[j];
            }

            NekDouble fac = 1.0;
            if (basis[i]->GetPointsType() == LibUtilities::eGaussRadauMAlpha1Beta0)
            {
                fac = 0.5;
            }
            else if (basis[i]->GetPointsType() == LibUtilities::eGaussRadauMAlpha2Beta0)
            {
                fac = 0.25;
            }

            m_w[i].resize(w.size());
            for (auto j = 0; j < w.size(); ++j)
            {
                m_w[i][j] = fac * w[j];
            }

            auto D = basis[i]->GetD()->GetPtr();
            m_D[i].resize(D.size());
            for (int j = 0; j < D.size(); ++j)
            {
                m_D[i][j] = D[j];
            }

            auto Z = basis[i]->GetZ();
            m_Z[i].resize(Z.size());
            for (int j = 0; j < Z.size(); ++j)
            {
                m_Z[i][j] = Z[j];
            }

        }
    }

    /// Set up Jacobian array for those operators that require geometric
    /// information.
    void SetJac(const Array<OneD, const NekDouble> &jac) final
    {
        if (DEFORMED)
        {
            int nq = m_nq[0];
            for (int i = 1; i < DIM; i++)
            {
                nq *= m_nq[i];
            }

            m_jac.resize(m_nBlocks*nq);

            alignas(vec_t::alignment) NekDouble tmp[vec_t::width];
            for (size_t block = 0; block < m_nBlocks; ++block)
            {
                for(size_t q = 0; q < nq; q++)
                {
                    for (int j = 0; j < vec_t::width; ++j)
                    {
                        tmp[j] = jac[block*nq*vec_t::width + nq*j + q]; //Unvalidated until I can get an actual deformed mesh.
                    }

                    //Order is [block][quadpt]
                    m_jac[block*nq + q].load(&tmp[0]);
                }
            }
        }
        else{
            m_jac.resize(m_nBlocks);

            alignas(vec_t::alignment) NekDouble tmp[vec_t::width];
            for (size_t i = 0; i < m_nBlocks; ++i)
            {
                for (int j = 0; j < vec_t::width; ++j)
                {
                    tmp[j] = jac[vec_t::width*i+j];
                }
                m_jac[i].load(&tmp[0]);
            }
        }
    }

    void SetDF(const Array<TwoD, const NekDouble> &df) final
    {
        constexpr unsigned int n_df = DIM * DIM;
        alignas(vec_t::alignment) NekDouble vec[vec_t::width];

        if (DEFORMED)
        {
            int nq = m_nq[0];
            for (int i = 1; i < DIM; ++i)
            {
                nq *= m_nq[i];
            }

            m_df.resize(m_nBlocks * n_df*nq);
            auto *df_ptr = &m_df[0];
            for (int e = 0; e < m_nBlocks; ++e)
            {
                for (int q = 0; q < nq; q++)
                {
                    for (int dir = 0; dir < n_df; ++dir, ++df_ptr)
                    {
                        for (int j = 0; j < vec_t::width; ++j)
                        {
                            vec[j] = df[dir][(vec_t::width*e + j)*nq + q];
                        }
                        (*df_ptr).load(&vec[0]);
                    }
                }
            }
        }
        else
        {
            m_df.resize(m_nBlocks * n_df);
            for (int e = 0; e < m_nBlocks; ++e)
            {
                for (int dir = 0; dir < n_df; ++dir)
                {
                    for (int j = 0; j < vec_t::width; ++j)
                    {
                        vec[j] = df[dir][vec_t::width*e + j];
                    }
                    // Must have all vec_t::width elemnts aligned to do a load.
                    m_df[e*n_df + dir].load(&vec[0]);
                }
            }
        }

    }

    int m_nBlocks;
    std::array<int, DIM> m_nm, m_nq;
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM> m_bdata;
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM> m_dbdata;
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM> m_D; //Derivatives
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM> m_Z; //Zeroes
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM> m_w; //Weights
    std::vector<vec_t, tinysimd::allocator<vec_t>> m_df; //Chain rule function deriviatives for each element (00, 10, 20, 30...)
    std::vector<vec_t, tinysimd::allocator<vec_t>> m_jac;

};

using OperatorFactory = LibUtilities::NekFactory<
    std::string,
    Operator,
    std::vector<LibUtilities::BasisSharedPtr>,
    int
    >;

OperatorFactory &GetOperatorFactory();

/// Helper function, get operator string
std::string GetOpstring(LibUtilities::ShapeType shape, bool deformed=false);

} // namespace MatrixFree
} // namespace Nektar

#endif