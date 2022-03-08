#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include "MatrixFreeDeclspec.h"

#include <iostream>

#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <Collections/Operator.h>
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

    /// Number of degrees of freedom that this operator will process.
    MATRIXFREE_EXPORT virtual NekDouble Ndof()
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

    MATRIXFREE_EXPORT virtual void SetJac(
     const std::shared_ptr<std::vector<vec_t, tinysimd::allocator<vec_t>>>  &jac) = 0;
    
    MATRIXFREE_EXPORT virtual void SetDF(
     const std::shared_ptr<std::vector<vec_t, tinysimd::allocator<vec_t>>>  &df) = 0;
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
                            Array<OneD, Array<OneD,   NekDouble> > &out) = 0;

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

    bool NeedsJac() final
    {
        return true;
    }

    bool NeedsDF() final
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
        m_basis(basis), m_nElmt(nElmt),
        m_lambda(1.0), m_isConstVarDiff(false),
        m_isVarDiff(false)
    {
        int n = m_basis.size();
        m_constVarDiff = Array<OneD, NekDouble>(n*(n+1)/2);
        int tp = 1;
        for (int bn = 0; bn < n; ++bn)
        {
            tp *= m_basis[bn]->GetNumPoints();
        }
        switch (n)
        {
            case 2:
                m_varD00 = Array<OneD, NekDouble>(tp);
                m_varD01 = Array<OneD, NekDouble>(tp);
                m_varD11 = Array<OneD, NekDouble>(tp);
                break;
            case 3:
                m_varD00 = Array<OneD, NekDouble>(tp);
                m_varD01 = Array<OneD, NekDouble>(tp);
                m_varD11 = Array<OneD, NekDouble>(tp);
                m_varD02 = Array<OneD, NekDouble>(tp);
                m_varD12 = Array<OneD, NekDouble>(tp);
                m_varD22 = Array<OneD, NekDouble>(tp);
                break;
            default:
                break;
        }
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

    inline void SetLambda(NekDouble lambda)
    {
        m_lambda = lambda;
    }
    

    inline void SetConstVarDiffusion(Array<OneD, NekDouble> diff)
    {
        m_isConstVarDiff = true;
        
        int n = m_basis.size();
        
        for(int i = 0; i < n*(n+1)/2; ++i)
        {
            m_constVarDiff[i] = diff[i];
        }
    }
    
    inline void SetVarDiffusion(Array<OneD, NekDouble> diff)
    {
        boost::ignore_unused(diff);
        m_isVarDiff = true;
        m_isConstVarDiff = false;
        
        int n = m_basis.size();
        int tp = 1;
        
        for (int bn = 0; bn < n; ++bn)
        {
            tp *= m_basis[bn]->GetNumPoints();
        }
        
        // fixed values for testing!
        for(int i = 0; i < tp; ++i)
        {
            switch (n)
            {
                case 2:
                    m_varD00[i] = diff[0];
                    m_varD01[i] = diff[1];
                    m_varD11[i] = diff[2];
                    break;
                case 3:
                    m_varD00[i] = diff[0];
                    m_varD01[i] = diff[1];
                    m_varD11[i] = diff[2];
                    m_varD02[i] = diff[3];
                    m_varD12[i] = diff[4];
                    m_varD22[i] = diff[5];
                    break;
                default:
                    break;       
            }
        }
    }
    
    

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
    NekDouble m_lambda;
    bool m_isConstVarDiff;
    Array<OneD, NekDouble> m_constVarDiff;
    bool m_isVarDiff;
    Array<OneD, NekDouble> m_varD00;
    Array<OneD, NekDouble> m_varD01;
    Array<OneD, NekDouble> m_varD11;
    Array<OneD, NekDouble> m_varD02;
    Array<OneD, NekDouble> m_varD12;
    Array<OneD, NekDouble> m_varD22;
};


template <int DIM, bool DEFORMED = false>
class Helper : virtual public Operator
{
protected:
    Helper(std::vector<LibUtilities::BasisSharedPtr> basis,
              int nElmt)
        : Operator()
    {
        // Sanity check: no padding yet!
        ASSERTL1(nElmt % vec_t::width == 0,
                 "Number of elements not divisible by vector "
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

    void SetJac(
     const std::shared_ptr<std::vector<vec_t, tinysimd::allocator<vec_t>>>  &jac)
        final
    {
        m_jac = jac; 
    }
    
    void SetDF(
     const std::shared_ptr<std::vector<vec_t, tinysimd::allocator<vec_t>>>  &df)
        final
    {
        m_df = df; 
    }
    
    int m_nBlocks;
    std::array<int, DIM> m_nm, m_nq;
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM> m_bdata;
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM> m_dbdata;
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM> m_D; //Derivatives
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM> m_Z; //Zeroes
    std::array<std::vector<vec_t, tinysimd::allocator<vec_t>>, DIM> m_w; //Weights
    std::shared_ptr<std::vector<vec_t, tinysimd::allocator<vec_t>>>  m_df;//Chain rule function deriviatives for each element (00, 10, 20, 30...)
    std::shared_ptr<std::vector<vec_t, tinysimd::allocator<vec_t>>>  m_jac;
    
};

using OperatorFactory = LibUtilities::NekFactory<
    std::string,
    Operator,
    std::vector<LibUtilities::BasisSharedPtr>,
    int
    >;

MATRIXFREE_EXPORT OperatorFactory &GetOperatorFactory();

/// Helper function, get operator string
MATRIXFREE_EXPORT std::string GetOpstring(LibUtilities::ShapeType shape, bool deformed=false);

} // namespace MatrixFree
} // namespace Nektar

#endif
