///////////////////////////////////////////////////////////////////////////////
//
// File: Helmholtz.cpp
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
// Description: Helmholtz operator implementations
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <MatrixFreeOps/Operator.hpp>

#include <Collections/Operator.h>
#include <Collections/MatrixFreeBase.h>
#include <Collections/Collection.h>
#include <Collections/IProduct.h>

using namespace std;

namespace Nektar {
namespace Collections {

using LibUtilities::eSegment;
using LibUtilities::eQuadrilateral;
using LibUtilities::eTriangle;
using LibUtilities::eHexahedron;
using LibUtilities::eTetrahedron;
using LibUtilities::ePrism;
using LibUtilities::ePyramid;

/**
 * @brief Helmholtz operator using LocalRegions implementation.
 */
class Helmholtz_NoCollection : public Operator
{
    public:
        OPERATOR_CREATE(Helmholtz_NoCollection)

        ~Helmholtz_NoCollection() final
        {
        }

        void operator()(
                const Array<OneD, const NekDouble> &entry0,
                      Array<OneD, NekDouble> &entry1,
                      Array<OneD, NekDouble> &entry2,
                      Array<OneD, NekDouble> &entry3,
                      Array<OneD, NekDouble> &wsp) final
        {
            boost::ignore_unused(entry2,entry3,wsp);

            unsigned int nmodes = m_expList[0]->GetNcoeffs();
            Array<OneD, NekDouble> tmp;

            for(int n = 0; n < m_numElmt; ++n)
            {
                StdRegions::StdMatrixKey mkey(StdRegions::eHelmholtz,
                                              (m_expList)[n]->DetShapeType(),
                                              *(m_expList)[n], m_factors);
                m_expList[n]->GeneralMatrixOp(entry0 + n *nmodes,
                                              tmp = entry1 + n * nmodes,
                                              mkey);
            }
        }

         void operator()(int dir,
                         const Array<OneD, const NekDouble> &input,
                         Array<OneD, NekDouble> &output,
                         Array<OneD, NekDouble> &wsp) final
        {
            boost::ignore_unused(dir, input, output, wsp);
            NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
        }

        virtual void CheckFactors(StdRegions::FactorMap factors,
                                  int coll_phys_offset)
        {
            boost::ignore_unused(coll_phys_offset);
            m_factors = factors;
        }


    protected:
        int                                         m_dim;
        int                                         m_coordim;
        vector<StdRegions::StdExpansionSharedPtr>   m_expList;
        StdRegions::FactorMap                       m_factors;

    private:
        Helmholtz_NoCollection(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData,
                StdRegions::FactorMap                     factors)
            : Operator(pCollExp, pGeomData, factors)
        {
            m_expList = pCollExp;
            m_dim     = pCollExp[0]->GetNumBases();
            m_coordim = pCollExp[0]->GetCoordim();

            m_factors = factors;
        }
};

/// Factory initialisation for the Helmholtz_NoCollection operators
OperatorKey Helmholtz_NoCollection::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment,       eHelmholtz, eNoCollection,false),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      eHelmholtz, eNoCollection,false),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      eHelmholtz, eNoCollection,true),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eHelmholtz, eNoCollection,false),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   eHelmholtz, eNoCollection,false),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   eHelmholtz, eNoCollection,true),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid,       eHelmholtz, eNoCollection,false),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         eHelmholtz, eNoCollection,false),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         eHelmholtz, eNoCollection,true),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron,    eHelmholtz, eNoCollection,false),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_Hex")
};


/**
 * @brief Helmholtz operator using LocalRegions implementation.
 */
class Helmholtz_IterPerExp : public Operator
{
    public:
        OPERATOR_CREATE(Helmholtz_IterPerExp)

        ~Helmholtz_IterPerExp() final
        {
        }

        void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD, NekDouble> &output,
                      Array<OneD, NekDouble> &output1,
                      Array<OneD, NekDouble> &output2,
                      Array<OneD, NekDouble> &wsp) final
        {
            boost::ignore_unused(output1,output2);

            const int nCoeffs = m_stdExp->GetNcoeffs();
            const int nPhys   = m_stdExp->GetTotPoints();

            ASSERTL1(input.size() >= m_numElmt*nCoeffs,
                     "input array size is insufficient");
            ASSERTL1(output.size() >= m_numElmt*nCoeffs,
                     "output array size is insufficient");

            Array<OneD, NekDouble> tmpphys, t1;
            Array<OneD, Array<OneD, NekDouble> > dtmp(3);
            Array<OneD, Array<OneD, NekDouble> > tmp(3);

            tmpphys = wsp;
            for(int i = 1; i < m_coordim+1; ++i)
            {
                dtmp[i-1] = wsp + i*nPhys;
                tmp [i-1] = wsp + (i+m_coordim)*nPhys;
            }

            for (int i = 0; i < m_numElmt; ++i)
            {
                m_stdExp->BwdTrans(input + i*nCoeffs, tmpphys);

                // local derivative
                m_stdExp->PhysDeriv(tmpphys, dtmp[0], dtmp[1], dtmp[2]);

                // determine mass matrix term
                if(m_isDeformed)
                {
                    Vmath::Vmul(nPhys,m_jac+i*nPhys,1,tmpphys,1,tmpphys,1);
                }
                else
                {
                    Vmath::Smul(nPhys,m_jac[i],tmpphys,1,tmpphys,1);
                }

                m_stdExp->IProductWRTBase(tmpphys,t1 = output + i*nCoeffs);
                Vmath::Smul(nCoeffs,m_lambda,output + i*nCoeffs,1,
                            t1 = output+i*nCoeffs,1);

                if(m_isDeformed)
                {
                    // calculate full derivative
                    for(int j = 0; j < m_coordim; ++j)
                    {
                        Vmath::Vmul(nPhys,
                                    m_derivFac[j*m_dim].origin() + i*nPhys, 1,
                                    &dtmp[0][0], 1, &tmp[j][0], 1);

                        for(int k = 1; k < m_dim; ++k)
                        {
                            Vmath::Vvtvp (nPhys, m_derivFac[j*m_dim+k].origin()
                                          + i*nPhys, 1, &dtmp[k][0], 1,
                                          &tmp[j][0],   1,  &tmp[j][0],   1);
                        }
                    }

                    if(m_HasVarCoeffDiff)
                    {
                        // calculate dtmp[i] = dx/dxi sum_j diff[0][j] tmp[j]
                        //                   + dy/dxi sum_j diff[1][j] tmp[j]
                        //                   + dz/dxi sum_j diff[2][j] tmp[j]

                        // First term
                        Vmath::Smul(nPhys, m_diff[0][0], &tmp[0][0],  1,
                                                         &tmpphys[0], 1);
                        for(int l = 1; l < m_coordim; ++l)
                        {
                            Vmath::Svtvp(nPhys, m_diff[0][l], &tmp[l][0], 1,
                                                &tmpphys[0], 1, &tmpphys[0], 1);
                        }

                        for(int j = 0; j < m_dim; ++j)
                        {
                            Vmath::Vmul(nPhys,
                                        m_derivFac[j].origin() + i*nPhys, 1,
                                        &tmpphys[0], 1, &dtmp[j][0], 1);
                        }

                        // Second and third terms
                        for(int k = 1; k < m_coordim; ++k)
                        {
                            Vmath::Smul(nPhys,m_diff[k][0], &tmp[0][0],  1,
                                                            &tmpphys[0], 1);
                            for(int l = 1; l < m_coordim; ++l)
                            {
                                Vmath::Svtvp(nPhys, m_diff[k][l], &tmp[l][0], 1,
                                             &tmpphys[0], 1, &tmpphys[0], 1);
                            }

                            for(int j = 0; j < m_dim; ++j)
                            {
                                Vmath::Vvtvp (nPhys,
                                              m_derivFac[j +k*m_dim].origin()
                                                + i*nPhys, 1,
                                              &tmpphys[0], 1,
                                              &dtmp[j][0], 1, &dtmp[j][0], 1);
                            }
                        }
                    }
                    else
                    {
                        // calculate dx/dxi tmp[0] + dy/dxi tmp[1]
                        //                         + dz/dxi tmp[2]
                        for(int j = 0; j < m_dim; ++j)
                        {
                            Vmath::Vmul (nPhys,
                                         m_derivFac[j].origin() + i*nPhys, 1,
                                         &tmp[0][0], 1, &dtmp[j][0], 1);

                            for(int k = 1; k < m_coordim; ++k)
                            {
                                Vmath::Vvtvp (nPhys,
                                              m_derivFac[j +k*m_dim].origin()
                                                + i*nPhys, 1, &tmp[k][0], 1,
                                              &dtmp[j][0], 1, &dtmp[j][0], 1);
                            }
                        }
                    }

                    // calculate Iproduct WRT Std Deriv
                    for(int j = 0; j < m_dim; ++j)
                    {

                        // multiply by Jacobian
                        Vmath::Vmul(nPhys,m_jac+i*nPhys,1,dtmp[j],1,dtmp[j],1);

                        m_stdExp->IProductWRTDerivBase(j,dtmp[j],tmp[0]);
                        Vmath::Vadd(nCoeffs,tmp[0],1,output+i*nCoeffs,1,
                                    t1 = output+i*nCoeffs,1);
                    }
                }
                else
                {
                    // calculate full derivative
                    for(int j = 0; j < m_coordim; ++j)
                    {
                        Vmath::Smul(nPhys, m_derivFac[j*m_dim][i],
                                    &dtmp[0][0], 1, &tmp[j][0], 1);

                        for(int k = 1; k < m_dim; ++k)
                        {
                            Vmath::Svtvp (nPhys,  m_derivFac[j*m_dim+k][i],
                                                  &dtmp[k][0], 1,
                                                  &tmp[j][0],  1,
                                                  &tmp[j][0],  1);
                        }

                    }

                    if(m_HasVarCoeffDiff)
                    {
                        // calculate dtmp[i] = dx/dxi sum_j diff[0][j] tmp[j]
                        //                   + dy/dxi sum_j diff[1][j] tmp[j]
                        //                   + dz/dxi sum_j diff[2][j] tmp[j]

                        // First term
                        Vmath::Smul(nPhys, m_diff[0][0], &tmp[0][0],  1,
                                                         &tmpphys[0], 1);
                        for(int l = 1; l < m_coordim; ++l)
                        {
                            Vmath::Svtvp(nPhys, m_diff[0][l], &tmp[l][0], 1,
                                         &tmpphys[0], 1, &tmpphys[0], 1);
                        }

                        for(int j = 0; j < m_dim; ++j)
                        {
                            Vmath::Smul (nPhys, m_derivFac[j][i],
                                                &tmpphys[0], 1, &dtmp[j][0], 1);
                        }

                        // Second and third terms
                        for(int k = 1; k < m_coordim; ++k)
                        {
                            Vmath::Smul(nPhys, m_diff[k][0], &tmp[0][0],  1,
                                                             &tmpphys[0], 1);
                            for(int l = 1; l < m_coordim; ++l)
                            {
                                Vmath::Svtvp(nPhys, m_diff[k][l], &tmp[l][0], 1,
                                                &tmpphys[0], 1, &tmpphys[0], 1);
                            }

                            for(int j = 0; j < m_dim; ++j)
                            {
                                Vmath::Svtvp (nPhys, m_derivFac[j +k*m_dim][i],
                                              &tmpphys[0], 1, &dtmp[j][0], 1,
                                              &dtmp[j][0], 1);
                            }
                        }
                    }
                    else
                    {
                        // calculate dx/dxi tmp[0] + dy/dxi tmp[2]
                        //                         + dz/dxi tmp[3]
                        for(int j = 0; j < m_dim; ++j)
                        {
                            Vmath::Smul (nPhys,m_derivFac[j][i],
                                         &tmp[0][0], 1, &dtmp[j][0],1);

                            for(int k = 1; k < m_coordim; ++k)
                            {
                                Vmath::Svtvp (nPhys, m_derivFac[j +k*m_dim][i],
                                              &tmp[k][0], 1, &dtmp[j][0], 1,
                                              &dtmp[j][0], 1);
                            }
                        }
                    }

                    // calculate Iproduct WRT Std Deriv
                    for(int j = 0; j < m_dim; ++j)
                    {
                        // multiply by Jacobian
                        Vmath::Smul(nPhys,m_jac[i],dtmp[j],1,dtmp[j],1);

                        m_stdExp->IProductWRTDerivBase(j,dtmp[j],tmp[0]);
                        Vmath::Vadd(nCoeffs,tmp[0],1,output+i*nCoeffs,1,
                                    t1 = output+i*nCoeffs,1);
                    }
                }
            }
        }

        void operator()(int dir,
                        const Array<OneD, const NekDouble> &input,
                        Array<OneD, NekDouble> &output,
                        Array<OneD, NekDouble> &wsp) final
        {
            boost::ignore_unused(dir, input, output, wsp);
            NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
        }


        /**
         * @brief Check the validity of supplied constant factors.
         *
         * @param factors Map of factors
         * @param coll_phys_offset Unused
         */
        virtual void CheckFactors(StdRegions::FactorMap factors,
                                  int coll_phys_offset)
        {
            boost::ignore_unused(coll_phys_offset);

            // If match previous factors, nothing to do.
            if (m_factors == factors)
            {
                return;
            }

            m_factors = factors;

            // Check Lambda constant of Helmholtz operator
            auto x = factors.find(StdRegions::eFactorLambda);
            ASSERTL1(x != factors.end(),
                     "Constant factor not defined: "
                     + std::string(StdRegions::ConstFactorTypeMap
                                   [StdRegions::eFactorLambda]));
            m_lambda = x->second;

            // If varcoeffs not supplied, nothing else to do.
            m_HasVarCoeffDiff = false;
            auto d = factors.find(StdRegions::eFactorCoeffD00);
            if (d == factors.end())
            {
                return;
            }

            m_diff = Array<OneD, Array<OneD, NekDouble> >(m_coordim);
            for(int i = 0; i < m_coordim; ++i)
            {
                m_diff[i] = Array<OneD, NekDouble>(m_coordim, 0.0);
            }

            for(int i = 0; i < m_coordim; ++i)
            {
                d = factors.find(m_factorCoeffDef[i][i]);
                if (d != factors.end())
                {
                    m_diff[i][i] = d->second;
                }
                else
                {
                    m_diff[i][i] = 1.0;
                }

                for(int j = i+1; j < m_coordim; ++j)
                {
                    d = factors.find(m_factorCoeffDef[i][j]);
                    if (d != factors.end())
                    {
                        m_diff[i][j] = m_diff[j][i] = d->second;
                    }
                }
            }
            m_HasVarCoeffDiff = true;
        }


    protected:
        Array<TwoD, const NekDouble>    m_derivFac;
        Array<OneD, const NekDouble>    m_jac;
        int                             m_dim;
        int                             m_coordim;
        StdRegions::FactorMap           m_factors;
        NekDouble                       m_lambda;
        bool                            m_HasVarCoeffDiff;
        Array<OneD, Array<OneD, NekDouble>>   m_diff;
        const StdRegions::ConstFactorType     m_factorCoeffDef[3][3] =
            {{StdRegions::eFactorCoeffD00,StdRegions::eFactorCoeffD01,
                    StdRegions::eFactorCoeffD02},
             {StdRegions::eFactorCoeffD01,StdRegions::eFactorCoeffD11,
                    StdRegions::eFactorCoeffD12},
             {StdRegions::eFactorCoeffD02,StdRegions::eFactorCoeffD12,
                    StdRegions::eFactorCoeffD22}};

    private:
        Helmholtz_IterPerExp(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData,
                StdRegions::FactorMap                     factors)
            : Operator(pCollExp, pGeomData, factors)
        {
            LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();
            m_dim      = PtsKey.size();
            m_coordim  = pCollExp[0]->GetCoordim();
            int nqtot  = m_stdExp->GetTotPoints();

            m_derivFac = pGeomData->GetDerivFactors(pCollExp);
            m_jac      = pGeomData->GetJac(pCollExp);
            m_wspSize = (2*m_coordim+1)*nqtot;

            m_lambda = 1.0;
            m_HasVarCoeffDiff = false;
            m_factors = StdRegions::NullFactorMap;
            this->CheckFactors(factors, 0);
        }
};

/// Factory initialisation for the Helmholtz_IterPerExp operators
OperatorKey Helmholtz_IterPerExp::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment,       eHelmholtz, eIterPerExp,false),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      eHelmholtz, eIterPerExp,false),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      eHelmholtz, eIterPerExp,true),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eHelmholtz, eIterPerExp,false),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   eHelmholtz, eIterPerExp,false),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   eHelmholtz, eIterPerExp,true),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid,       eHelmholtz, eIterPerExp,false),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         eHelmholtz, eIterPerExp,false),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         eHelmholtz, eIterPerExp,true),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron,    eHelmholtz, eIterPerExp,false),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_Hex")
};


/**
 * @brief Helmholtz operator using matrix free operators.
 */
class Helmholtz_MatrixFree : public Operator, MatrixFreeOneInOneOut
{
public:
    OPERATOR_CREATE(Helmholtz_MatrixFree)

    ~Helmholtz_MatrixFree() final
    {
    }

    void operator()(
            const Array<OneD, const NekDouble> &input,
                  Array<OneD,       NekDouble> &output0,
                  Array<OneD,       NekDouble> &output1,
                  Array<OneD,       NekDouble> &output2,
                  Array<OneD,       NekDouble> &wsp) final
    {
        boost::ignore_unused(output1,output2,wsp);

        if (m_isPadded)
        {
            // copy into padded vector
            Vmath::Vcopy(m_nmtot, input, 1, m_input, 1);
            (*m_oper)(m_input, m_output);
            Vmath::Vcopy(m_nmtot, m_output, 1, output0, 1);
        }
        else
        {
            (*m_oper)(input, output0);
        }
    }

    void operator()(int dir,
                    const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output,
                    Array<OneD, NekDouble> &wsp) final
    {
        boost::ignore_unused(dir,input,output,wsp);
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }


    /**
     *
     */
    virtual void CheckFactors(StdRegions::FactorMap factors,
                              int coll_phys_offset)
    {
        boost::ignore_unused(coll_phys_offset);

        if (factors == m_factors)
        {
            return;
        }

        m_factors = factors;

        // Set lambda for this call
        auto x = factors.find(StdRegions::eFactorLambda);
        ASSERTL1(x != factors.end(),
                     "Constant factor not defined: "
                     + std::string(StdRegions::ConstFactorTypeMap[
                                            StdRegions::eFactorLambda]));
        m_oper->SetLambda(x->second);

        // set constant diffusion coefficients
        bool isConstVarDiff = false;
        Array<OneD, NekDouble> diff = Array<OneD, NekDouble> (6, 0.0);
        diff[0] = diff[2] = diff[5] = 1.0;

        auto xd00 = factors.find(StdRegions::eFactorCoeffD00);
        if (xd00 != factors.end() && xd00->second != 1.0) {
            isConstVarDiff = true;
            diff[0] = xd00->second;
        }

        auto xd01 = factors.find(StdRegions::eFactorCoeffD01);
        if (xd01 != factors.end() && xd01->second != 0.0) {
            isConstVarDiff = true;
            diff[1] = xd01->second;
        }

        auto xd11 = factors.find(StdRegions::eFactorCoeffD11);
        if (xd11 != factors.end() && xd11->second != 1.0) {
            isConstVarDiff = true;
            diff[2] = xd11->second;
        }

        auto xd02 = factors.find(StdRegions::eFactorCoeffD02);
        if (xd02 != factors.end() && xd02->second != 0.0) {
            isConstVarDiff = true;
            diff[3] = xd02->second;
        }

        auto xd12 = factors.find(StdRegions::eFactorCoeffD12);
        if (xd12 != factors.end() && xd12->second != 0.0) {
            isConstVarDiff = true;
            diff[4] = xd12->second;
        }

        auto xd22 = factors.find(StdRegions::eFactorCoeffD22);
        if (xd22 != factors.end() && xd22->second != 1.0) {
            isConstVarDiff = true;
            diff[5] = xd22->second;
        }

        if (isConstVarDiff) {
            m_oper->SetConstVarDiffusion(diff);
        }

        // set random here for fn variable diffusion
        auto k = factors.find(StdRegions::eFactorTau);
        if (k != factors.end() && k->second != 0.0) {
            m_oper->SetVarDiffusion(diff);
        }
    }

private:
    std::shared_ptr<MatrixFree::Helmholtz> m_oper;
    unsigned int m_nmtot;
    StdRegions::FactorMap m_factors;

    Helmholtz_MatrixFree(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                         CoalescedGeomDataSharedPtr                pGeomData,
                         StdRegions::FactorMap                     factors)
        : Operator(pCollExp, pGeomData, factors),
          MatrixFreeOneInOneOut(pCollExp[0]->GetStdExp()->GetNcoeffs(),
                                pCollExp[0]->GetStdExp()->GetNcoeffs(),
                                pCollExp.size())
    {

        m_nmtot = m_numElmt*pCollExp[0]->GetStdExp()->GetNcoeffs();

        const auto dim = pCollExp[0]->GetStdExp()->GetShapeDimension();

        // Basis vector.
        std::vector<LibUtilities::BasisSharedPtr> basis(dim);
        for (auto i = 0; i < dim; ++i)
        {
            basis[i] = pCollExp[0]->GetBasis(i);
        }

        // Get shape type
        auto shapeType = pCollExp[0]->GetStdExp()->DetShapeType();

        // Generate operator string and create operator.
        std::string op_string = "Helmholtz";
        op_string += MatrixFree::GetOpstring(shapeType, m_isDeformed);
        auto oper = MatrixFree::GetOperatorFactory().
            CreateInstance(op_string, basis, m_nElmtPad);

        // Set Jacobian
        oper->SetJac(pGeomData->GetJacInterLeave(pCollExp,m_nElmtPad));

        // Store derivative factor
        oper->SetDF(pGeomData->GetDerivFactorsInterLeave
                    (pCollExp,m_nElmtPad));

        m_oper = std::dynamic_pointer_cast<MatrixFree::Helmholtz>(oper);
        ASSERTL0(m_oper, "Failed to cast pointer.");

        // Set factors
        m_factors = StdRegions::NullFactorMap;
        this->CheckFactors(factors, 0);
    }
};

/// Factory initialisation for the Helmholtz_MatrixFree operators
OperatorKey Helmholtz_MatrixFree::m_typeArr[] =
{
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Hex"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Tet"),
};

}
}

