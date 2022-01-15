#ifndef NEKTAR_LIBRARY_MF_HELMHOLTZ_H
#define NEKTAR_LIBRARY_MF_HELMHOLTZ_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include "IProduct.h"
#include "IProductKernels.hpp"
#include "BwdTransKernels.hpp"
#include "BwdTrans.h"
#include "PhysDeriv.h"

namespace Nektar
{
namespace MatrixFree
{

template<bool DEFORMED = false>
struct HelmholtzQuad : public Helmholtz, public Helper<2, DEFORMED>
{
    HelmholtzQuad(std::vector<LibUtilities::BasisSharedPtr> basis,
                     int nElmt)
    : Helmholtz(basis, nElmt),
      Helper<2, DEFORMED>(basis, nElmt),
      m_nmTot(LibUtilities::StdQuadData::getNumberOfCoefficients(
                                 this->m_nm[0], this->m_nm[1]))
    {
    }

    static std::shared_ptr<Operator> Create(
                     std::vector<LibUtilities::BasisSharedPtr> basis,
                     int nElmt)
    {
        return std::make_shared<HelmholtzQuad<DEFORMED>>(basis, nElmt);
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
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
                case 2: HelmholtzQuadImpl<2 ,2 ,2 ,2 >(in, out); break;
                case 3: HelmholtzQuadImpl<2 ,2 ,3 ,3 >(in, out); break;
                case 4: HelmholtzQuadImpl<2 ,2 ,4 ,4 >(in, out); break;
                default: HelmholtzQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
                } break;
            case 3:
                switch(nq0)
                {
                case 3: HelmholtzQuadImpl<3 ,3 ,3 ,3 >(in, out); break;
                case 4: HelmholtzQuadImpl<3 ,3 ,4 ,4 >(in, out); break;
                case 5: HelmholtzQuadImpl<3 ,3 ,5 ,5 >(in, out); break;
                case 6: HelmholtzQuadImpl<3 ,3 ,6 ,6 >(in, out); break;
                default: HelmholtzQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
                } break;
            case 4:
                switch(nq0)
                {
                case 4: HelmholtzQuadImpl<4 ,4 ,4 ,4 >(in, out); break;
                case 5: HelmholtzQuadImpl<4 ,4 ,5 ,5 >(in, out); break;
                case 6: HelmholtzQuadImpl<4 ,4 ,6 ,6 >(in, out); break;
                case 7: HelmholtzQuadImpl<4 ,4 ,7 ,7 >(in, out); break;
                case 8: HelmholtzQuadImpl<4 ,4 ,8 ,8 >(in, out); break;
                default: HelmholtzQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
                } break;
            case 5:
                switch(nq0)
                {
                case 5: HelmholtzQuadImpl<5 ,5 ,5 ,5 >(in, out); break;
                case 6: HelmholtzQuadImpl<5 ,5 ,6 ,6 >(in, out); break;
                case 7: HelmholtzQuadImpl<5 ,5 ,7 ,7 >(in, out); break;
                case 8: HelmholtzQuadImpl<5 ,5 ,8 ,8 >(in, out); break;
                case 9: HelmholtzQuadImpl<5 ,5 ,9 ,9 >(in, out); break;
                case 10: HelmholtzQuadImpl<5 ,5 ,10 ,10 >(in, out); break;
                default: HelmholtzQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
                } break;
            case 6:
                switch(nq0)
                {
                case 6: HelmholtzQuadImpl<6 ,6 ,6 ,6 >(in, out); break;
                case 7: HelmholtzQuadImpl<6 ,6 ,7 ,7 >(in, out); break;
                case 8: HelmholtzQuadImpl<6 ,6 ,8 ,8 >(in, out); break;
                case 9: HelmholtzQuadImpl<6 ,6 ,9 ,9 >(in, out); break;
                case 10: HelmholtzQuadImpl<6 ,6 ,10 ,10 >(in, out); break;
                case 11: HelmholtzQuadImpl<6 ,6 ,11 ,11 >(in, out); break;
                case 12: HelmholtzQuadImpl<6 ,6 ,12 ,12 >(in, out); break;
                default: HelmholtzQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
                } break;
            case 7:
                switch(nq0)
                {
                case 7: HelmholtzQuadImpl<7 ,7 ,7 ,7 >(in, out); break;
                case 8: HelmholtzQuadImpl<7 ,7 ,8 ,8 >(in, out); break;
                case 9: HelmholtzQuadImpl<7 ,7 ,9 ,9 >(in, out); break;
                case 10: HelmholtzQuadImpl<7 ,7 ,10 ,10 >(in, out); break;
                case 11: HelmholtzQuadImpl<7 ,7 ,11 ,11 >(in, out); break;
                case 12: HelmholtzQuadImpl<7 ,7 ,12 ,12 >(in, out); break;
                case 13: HelmholtzQuadImpl<7 ,7 ,13 ,13 >(in, out); break;
                case 14: HelmholtzQuadImpl<7 ,7 ,14 ,14 >(in, out); break;
                default: HelmholtzQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
                } break;
            case 8:
                switch(nq0)
                {
                case 8: HelmholtzQuadImpl<8 ,8 ,8 ,8 >(in, out); break;
                case 9: HelmholtzQuadImpl<8 ,8 ,9 ,9 >(in, out); break;
                case 10: HelmholtzQuadImpl<8 ,8 ,10 ,10 >(in, out); break;
                case 11: HelmholtzQuadImpl<8 ,8 ,11 ,11 >(in, out); break;
                case 12: HelmholtzQuadImpl<8 ,8 ,12 ,12 >(in, out); break;
                case 13: HelmholtzQuadImpl<8 ,8 ,13 ,13 >(in, out); break;
                case 14: HelmholtzQuadImpl<8 ,8 ,14 ,14 >(in, out); break;
                case 15: HelmholtzQuadImpl<8 ,8 ,15 ,15 >(in, out); break;
                case 16: HelmholtzQuadImpl<8 ,8 ,16 ,16 >(in, out); break;
                default: HelmholtzQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
                } break;
            default: HelmholtzQuadImpl(nm0,nm1,nq0,nq1,in,out); break;
            }
        }
        else
        {
            HelmholtzQuadImpl(nm0,nm1,nq0,nq1,in,out);
        }
    }    

    template<int NM0, int NM1, int NQ0, int NQ1>
    void HelmholtzQuadImpl(
          const Array<OneD, const NekDouble> &input, Array<OneD, NekDouble> &out)
    {
        auto *inptr  = &input[0];
        auto *outptr = &out[0];

        constexpr auto ndf = 4;
        constexpr auto nqTot = NQ0 * NQ1;
        const auto nmBlocks = m_nmTot * vec_t::width;

        // Allocate sufficient workspace for backwards transform and inner
        // product kernels.
        constexpr auto wspInnerProd = NQ1;
        constexpr auto wspBwdTrans = NQ0 * NM0;
        constexpr auto wspSize = wspInnerProd > wspBwdTrans ?
            wspInnerProd : wspBwdTrans;

        vec_t wsp[wspSize]; // workspace for kernels

        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot),  deriv0(nqTot),
            deriv1(nqTot);

        const vec_t* jac_ptr;
        const vec_t* df_ptr;

        vec_t df0,df1,df2,df3;
        vec_t metric00,metric01,metric11;
        
        vec_t d00 = {1.0};
        vec_t d01 = {0.0};
        vec_t d11 = {1.0}; // var diffusion terms
        vec_t dtmp0,dtmp1,dtmp2,dtmp3; // temp for vardiff
        boost::ignore_unused(d00,d01,d11,dtmp0,dtmp1,dtmp2,dtmp3);
        
        if (this->m_isConstVarDiff)
        {
            d00 = this->m_constVarDiff[0];
            d01 = this->m_constVarDiff[1];
            d11 = this->m_constVarDiff[2];
        }
                
        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            df_ptr = &((*this->m_df)[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            // Precompute Laplacian metricsp
            if(!DEFORMED)
            {
                df0 = df_ptr[0];  df1 = df_ptr[1];
                df2 = df_ptr[2];  df3 = df_ptr[3];

                if (!this->m_isConstVarDiff && !this->m_isVarDiff)
                {
                    metric00 = df0*df0;
                    metric00.fma(df2,df2);
                    metric01 = df0*df1; 
                    metric01.fma(df2,df3);
                    metric11 = df1*df1;
                    metric11.fma(df3,df3);
                }
                else if (this->m_isConstVarDiff)
                {
                    dtmp0 = df0 * d00;
                    dtmp0.fma(df2,d01);
                    dtmp1 = df0 * d01;
                    dtmp1.fma(df2,d11);
                    dtmp2 = df1 * d00;
                    dtmp2.fma(df3,d01);
                    dtmp3 = df1 * d01;
                    dtmp3.fma(df3,d11);
                    
                    metric00 = df0 * dtmp0;
                    metric00.fma(df2,dtmp1);
                    metric01 = df1 * dtmp0;
                    metric01.fma(df3,dtmp1);
                    metric11 = df1 * dtmp2;
                    metric11.fma(df3,dtmp3);
                }

                jac_ptr = &((*this->m_jac)[e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e*nqTot]);
            }

            // Step 1: BwdTrans
            BwdTransQuadKernel(NM0, NM1, NQ0, NQ1,
                    tmpIn, this->m_bdata[0], this->m_bdata[1], wsp, bwd);

            // Step 2: inner product for mass matrix operation
            IProductQuadKernel<NM0, NM1, NQ0, NQ1, true, false, DEFORMED>
                (bwd, this->m_bdata[0], this->m_bdata[1], this->m_w[0],
                 this->m_w[1], jac_ptr, wsp, tmpOut, m_lambda);

            // Step 3: take derivatives in quadrature space
            PhysDerivTensor2DKernel(NQ0, NQ1, 
                                    bwd, this->m_D[0], this->m_D[1], deriv0, deriv1);
            
            // Step 4: Apply Laplacian metrics & inner product
            if (!this->m_isVarDiff) 
            {
                if(DEFORMED)
                {
                    for(size_t j = 0, cnt = 0; j < NQ1; ++j)
                    {
                        for (size_t i = 0; i < NQ0; ++i, ++cnt)
                        {
                            df0 = df_ptr[cnt * ndf];
                            df1 = df_ptr[cnt * ndf + 1];
                            df2 = df_ptr[cnt * ndf + 2];
                            df3 = df_ptr[cnt * ndf + 3];
                            
                            if (!this->m_isConstVarDiff)
                            {
                                metric00 = df0*df0;
                                metric00.fma(df2,df2);
                                metric01 = df0*df1; 
                                metric01.fma(df2,df3);
                                metric11 = df1*df1;
                                metric11.fma(df3,df3);
                            }
                            else
                            {
                                dtmp0 = df0 * d00;
                                dtmp0.fma(df2,d01);
                                dtmp1 = df0 * d01;
                                dtmp1.fma(df2,d11);
                                dtmp2 = df1 * d00;
                                dtmp2.fma(df3,d01);
                                dtmp3 = df1 * d01;
                                dtmp3.fma(df3,d11);
                                
                                metric00 = df0 * dtmp0;
                                metric00.fma(df2,dtmp1);
                                metric01 = df1 * dtmp0;
                                metric01.fma(df3,dtmp1);
                                metric11 = df1 * dtmp2;
                                metric11.fma(df3,dtmp3);
                            }

                            vec_t d0 = deriv0[cnt];
                            vec_t d1 = deriv1[cnt];

                            vec_t tmp = metric00 * d0;
                            tmp.fma(metric01, d1);
                            bwd[cnt]  = tmp;
                            
                            tmp = metric01 * d0;
                            tmp.fma(metric11, d1);
                            deriv0[cnt] = tmp;
                        }
                    }
                }
                else
                {   
                    for (int i = 0; i < nqTot; ++i)
                    {
                        vec_t d0 = deriv0[i];
                        vec_t d1 = deriv1[i];

                        vec_t tmp = metric00 * d0;
                        tmp.fma(metric01, d1);
                        bwd[i] = tmp;
                        
                        tmp = metric01 * d0;
                        tmp.fma(metric11, d1);
                        deriv0[i] = tmp; 
                    }
                }
            }
            else
            {
                if(DEFORMED)
                {
                    for(size_t j = 0, cnt = 0; j < NQ1; ++j)
                    {
                        for (size_t i = 0; i < NQ0; ++i, ++cnt)
                        {
                            df0 = df_ptr[cnt * ndf];
                            df1 = df_ptr[cnt * ndf + 1];
                            df2 = df_ptr[cnt * ndf + 2];
                            df3 = df_ptr[cnt * ndf + 3];
                            
                            d00 = m_varD00[cnt];
                            d01 = m_varD01[cnt];
                            d11 = m_varD11[cnt];
                            
                            dtmp0 = df0 * d00;
                            dtmp0.fma(df2,d01);
                            dtmp1 = df0 * d01;
                            dtmp1.fma(df2,d11);
                            dtmp2 = df1 * d00;
                            dtmp2.fma(df3,d01);
                            dtmp3 = df1 * d01;
                            dtmp3.fma(df3,d11);
                            
                            metric00 = df0 * dtmp0;
                            metric00.fma(df2,dtmp1);
                            metric01 = df1 * dtmp0;
                            metric01.fma(df3,dtmp1);
                            metric11 = df1 * dtmp2;
                            metric11.fma(df3,dtmp3);

                            vec_t d0 = deriv0[cnt];
                            vec_t d1 = deriv1[cnt];

                            vec_t tmp = metric00 * d0;
                            tmp.fma(metric01, d1);
                            bwd[cnt]  = tmp;
                            
                            tmp = metric01 * d0;
                            tmp.fma(metric11, d1);
                            deriv0[cnt] = tmp;
                        }
                    }
                }
                else
                {   
                    for(size_t j = 0, cnt = 0; j < NQ1; ++j)
                    {
                        for (size_t i = 0; i < NQ0; ++i, ++cnt)
                        {
                            d00 = m_varD00[cnt];
                            d01 = m_varD01[cnt];
                            d11 = m_varD11[cnt];
                            
                            dtmp0 = df0 * d00;
                            dtmp0.fma(df2,d01);
                            dtmp1 = df0 * d01;
                            dtmp1.fma(df2,d11);
                            dtmp2 = df1 * d00;
                            dtmp2.fma(df3,d01);
                            dtmp3 = df1 * d01;
                            dtmp3.fma(df3,d11);
                            
                            metric00 = df0 * dtmp0;
                            metric00.fma(df2,dtmp1);
                            metric01 = df1 * dtmp0;
                            metric01.fma(df3,dtmp1);
                            metric11 = df1 * dtmp2;
                            metric11.fma(df3,dtmp3);

                            vec_t d0 = deriv0[cnt];
                            vec_t d1 = deriv1[cnt];

                            vec_t tmp = metric00 * d0;
                            tmp.fma(metric01, d1);
                            bwd[cnt]  = tmp;
                            
                            tmp = metric01 * d0;
                            tmp.fma(metric11, d1);
                            deriv0[cnt] = tmp;
                        }
                    }
                }
            }
            
            IProductQuadKernel<NM0, NM1, NQ0, NQ1, false, true, DEFORMED>
                (bwd, this->m_dbdata[0], this->m_bdata[1], this->m_w[0],
                 this->m_w[1], jac_ptr, wsp, tmpOut);
            
            IProductQuadKernel<NM0, NM1, NQ0, NQ1, false, true, DEFORMED>
                (deriv0, this->m_bdata[0], this->m_dbdata[1], this->m_w[0],
                 this->m_w[1], jac_ptr, wsp, tmpOut);
           
            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);
            
            inptr  += nmBlocks;
            outptr += nmBlocks;
        }
    }

    void HelmholtzQuadImpl(
          const int nm0, const int nm1,
          const int nq0, const int nq1, 
          const Array<OneD, const NekDouble> &input, Array<OneD, NekDouble> &out)
    {
        auto *inptr  = &input[0];
        auto *outptr = &out[0];

        const auto ndf = 4;
        const auto nqTot = nq0 * nq1;
        const auto nmBlocks = m_nmTot * vec_t::width;

        // Allocate sufficient workspace for backwards transform and inner
        // product kernels.
        const auto wspInnerProd = nq1;
        const auto wspBwdTrans = nq0 * nm0;
        const auto wspSize = wspInnerProd > wspBwdTrans ?
            wspInnerProd : wspBwdTrans;

        std::vector<vec_t, allocator<vec_t>> wsp(wspSize), tmpIn(m_nmTot),
            tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot),  deriv0(nqTot),
            deriv1(nqTot);

        const vec_t* jac_ptr;
        const vec_t* df_ptr;

        vec_t df0,df1,df2,df3;
        vec_t metric00,metric01,metric11;
        
        vec_t d00 = {1.0};
        vec_t d01 = {0.0};
        vec_t d11 = {1.0}; // var diffusion terms
        vec_t dtmp0,dtmp1,dtmp2,dtmp3; // temp for vardiff
        boost::ignore_unused(d00,d01,d11,dtmp0,dtmp1,dtmp2,dtmp3);
        
        if (this->m_isConstVarDiff)
        {
            d00 = this->m_constVarDiff[0];
            d01 = this->m_constVarDiff[1];
            d11 = this->m_constVarDiff[2];
        }
                
        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            df_ptr = &((*this->m_df)[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            // Precompute Laplacian metricsp
            if(!DEFORMED)
            {
                df0 = df_ptr[0];  df1 = df_ptr[1];
                df2 = df_ptr[2];  df3 = df_ptr[3];

                if (!this->m_isConstVarDiff && !this->m_isVarDiff)
                {
                    metric00 = df0*df0;
                    metric00.fma(df2,df2);
                    metric01 = df0*df1; 
                    metric01.fma(df2,df3);
                    metric11 = df1*df1;
                    metric11.fma(df3,df3);
                }
                else if (this->m_isConstVarDiff)
                {
                    dtmp0 = df0 * d00;
                    dtmp0.fma(df2,d01);
                    dtmp1 = df0 * d01;
                    dtmp1.fma(df2,d11);
                    dtmp2 = df1 * d00;
                    dtmp2.fma(df3,d01);
                    dtmp3 = df1 * d01;
                    dtmp3.fma(df3,d11);
                    
                    metric00 = df0 * dtmp0;
                    metric00.fma(df2,dtmp1);
                    metric01 = df1 * dtmp0;
                    metric01.fma(df3,dtmp1);
                    metric11 = df1 * dtmp2;
                    metric11.fma(df3,dtmp3);
                }

                jac_ptr = &((*this->m_jac)[e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e*nqTot]);
            }

            // Step 1: BwdTrans
            BwdTransQuadKernel(nm0, nm1, nq0, nq1,
                    tmpIn, this->m_bdata[0], this->m_bdata[1], &wsp[0], bwd);

            // Step 2: inner product for mass matrix operation
            IProductQuadKernel(nm0, nm1, nq0, nq1, true, false, DEFORMED,
                 bwd, this->m_bdata[0], this->m_bdata[1], this->m_w[0],
                 this->m_w[1], jac_ptr, &wsp[0], tmpOut, m_lambda);

            // Step 3: take derivatives in quadrature space
            PhysDerivTensor2DKernel(nq0, nq1, bwd, this->m_D[0],
                                    this->m_D[1], deriv0, deriv1);
            
            // Step 4: Apply Laplacian metrics & inner product
            if (!this->m_isVarDiff) 
            {
                if(DEFORMED)
                {
                    for(size_t j = 0, cnt = 0; j < nq1; ++j)
                    {
                        for (size_t i = 0; i < nq0; ++i, ++cnt)
                        {
                            df0 = df_ptr[cnt * ndf];
                            df1 = df_ptr[cnt * ndf + 1];
                            df2 = df_ptr[cnt * ndf + 2];
                            df3 = df_ptr[cnt * ndf + 3];
                            
                            if (!this->m_isConstVarDiff)
                            {
                                metric00 = df0*df0;
                                metric00.fma(df2,df2);
                                metric01 = df0*df1; 
                                metric01.fma(df2,df3);
                                metric11 = df1*df1;
                                metric11.fma(df3,df3);
                            }
                            else
                            {
                                dtmp0 = df0 * d00;
                                dtmp0.fma(df2,d01);
                                dtmp1 = df0 * d01;
                                dtmp1.fma(df2,d11);
                                dtmp2 = df1 * d00;
                                dtmp2.fma(df3,d01);
                                dtmp3 = df1 * d01;
                                dtmp3.fma(df3,d11);
                                
                                metric00 = df0 * dtmp0;
                                metric00.fma(df2,dtmp1);
                                metric01 = df1 * dtmp0;
                                metric01.fma(df3,dtmp1);
                                metric11 = df1 * dtmp2;
                                metric11.fma(df3,dtmp3);
                            }

                            vec_t d0 = deriv0[cnt];
                            vec_t d1 = deriv1[cnt];

                            vec_t tmp = metric00 * d0;
                            tmp.fma(metric01, d1);
                            bwd[cnt]  = tmp;
                            
                            tmp = metric01 * d0;
                            tmp.fma(metric11, d1);
                            deriv0[cnt] = tmp;
                        }
                    }
                }
                else
                {   
                    for (int i = 0; i < nqTot; ++i)
                    {
                        vec_t d0 = deriv0[i];
                        vec_t d1 = deriv1[i];

                        vec_t tmp = metric00 * d0;
                        tmp.fma(metric01, d1);
                        bwd[i] = tmp;
                        
                        tmp = metric01 * d0;
                        tmp.fma(metric11, d1);
                        deriv0[i] = tmp; 
                    }
                }
            }
            else
            {
                if(DEFORMED)
                {
                    for(size_t j = 0, cnt = 0; j < nq1; ++j)
                    {
                        for (size_t i = 0; i < nq0; ++i, ++cnt)
                        {
                            df0 = df_ptr[cnt * ndf];
                            df1 = df_ptr[cnt * ndf + 1];
                            df2 = df_ptr[cnt * ndf + 2];
                            df3 = df_ptr[cnt * ndf + 3];
                            
                            d00 = m_varD00[cnt];
                            d01 = m_varD01[cnt];
                            d11 = m_varD11[cnt];
                            
                            dtmp0 = df0 * d00;
                            dtmp0.fma(df2,d01);
                            dtmp1 = df0 * d01;
                            dtmp1.fma(df2,d11);
                            dtmp2 = df1 * d00;
                            dtmp2.fma(df3,d01);
                            dtmp3 = df1 * d01;
                            dtmp3.fma(df3,d11);
                            
                            metric00 = df0 * dtmp0;
                            metric00.fma(df2,dtmp1);
                            metric01 = df1 * dtmp0;
                            metric01.fma(df3,dtmp1);
                            metric11 = df1 * dtmp2;
                            metric11.fma(df3,dtmp3);

                            vec_t d0 = deriv0[cnt];
                            vec_t d1 = deriv1[cnt];

                            vec_t tmp = metric00 * d0;
                            tmp.fma(metric01, d1);
                            bwd[cnt]  = tmp;
                            
                            tmp = metric01 * d0;
                            tmp.fma(metric11, d1);
                            deriv0[cnt] = tmp;
                        }
                    }
                }
                else
                {   
                    for(size_t j = 0, cnt = 0; j < nq1; ++j)
                    {
                        for (size_t i = 0; i < nq0; ++i, ++cnt)
                        {
                            d00 = m_varD00[cnt];
                            d01 = m_varD01[cnt];
                            d11 = m_varD11[cnt];
                            
                            dtmp0 = df0 * d00;
                            dtmp0.fma(df2,d01);
                            dtmp1 = df0 * d01;
                            dtmp1.fma(df2,d11);
                            dtmp2 = df1 * d00;
                            dtmp2.fma(df3,d01);
                            dtmp3 = df1 * d01;
                            dtmp3.fma(df3,d11);
                            
                            metric00 = df0 * dtmp0;
                            metric00.fma(df2,dtmp1);
                            metric01 = df1 * dtmp0;
                            metric01.fma(df3,dtmp1);
                            metric11 = df1 * dtmp2;
                            metric11.fma(df3,dtmp3);

                            vec_t d0 = deriv0[cnt];
                            vec_t d1 = deriv1[cnt];

                            vec_t tmp = metric00 * d0;
                            tmp.fma(metric01, d1);
                            bwd[cnt]  = tmp;
                            
                            tmp = metric01 * d0;
                            tmp.fma(metric11, d1);
                            deriv0[cnt] = tmp;
                        }
                    }
                }
            }
            
            IProductQuadKernel(nm0, nm1, nq0, nq1, false, true, DEFORMED,
                bwd, this->m_dbdata[0], this->m_bdata[1], this->m_w[0],
                 this->m_w[1], jac_ptr, &wsp[0], tmpOut);
            
            IProductQuadKernel(nm0, nm1, nq0, nq1, false, true, DEFORMED,
                 deriv0, this->m_bdata[0], this->m_dbdata[1], this->m_w[0],
                 this->m_w[1], jac_ptr, &wsp[0], tmpOut);
           
            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);
            
            inptr  += nmBlocks;
            outptr += nmBlocks;
        }
    }

public:

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }

private:
    const int m_nmTot;
};


    template<bool DEFORMED = false>
    struct HelmholtzTri : public Helmholtz, public Helper<2, DEFORMED>
    {
        HelmholtzTri(std::vector<LibUtilities::BasisSharedPtr> basis,
                     int nElmt)
            : Helmholtz(basis, nElmt),
              Helper<2, DEFORMED>(basis, nElmt),
              m_nmTot(LibUtilities::StdTriData::getNumberOfCoefficients(
                                                                        this->m_nm[0], this->m_nm[1])),
              m_h0(basis[0]->GetNumPoints()),
              m_h1(basis[1]->GetNumPoints())
        {
            const int nq0 = basis[0]->GetNumPoints();
            const int nq1 = basis[1]->GetNumPoints();

            const Array<OneD, const NekDouble> &z0 = basis[0]->GetZ();
            const Array<OneD, const NekDouble> &z1 = basis[1]->GetZ();

            for (int i = 0; i < nq0; ++i)
            {
                m_h0[i] = 0.5 * (1 + z0[i]);
            }

            for (int j = 0; j < nq1; ++j)
            {
                m_h1[j] = 2.0 / (1 - z1[j]);
            }
        }

            static std::shared_ptr<Operator> Create(
                                                    std::vector<LibUtilities::BasisSharedPtr> basis,
                                                    int nElmt)
        {
            return std::make_shared<HelmholtzTri<DEFORMED>>(basis, nElmt);
        }

        virtual void operator()(const Array<OneD, const NekDouble> &in,
                                Array<OneD,       NekDouble> &out) final
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
                    case 2: switch(nq0)
                        {
                        case 3: HelmholtzTriImpl<2, 2, 3, 2, true>(in, out); break;
                        case 4: HelmholtzTriImpl<2, 2, 4, 3, true>(in, out); break;
                        default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                        } break;
                    case 3:
                        switch(nq0)
                        {
                        case 4: HelmholtzTriImpl<3, 3, 4, 3, true>(in, out); break;
                        case 5: HelmholtzTriImpl<3, 3, 5, 4, true>(in, out); break;
                        case 6: HelmholtzTriImpl<3, 3, 6, 5, true>(in, out); break;
                        default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                        } break;
                    case 4:
                        switch(nq0)
                        {
                        case 5: HelmholtzTriImpl<4, 4, 5, 4, true>(in, out); break;
                        case 6: HelmholtzTriImpl<4, 4, 6, 5, true>(in, out); break;
                        case 7: HelmholtzTriImpl<4, 4, 7, 6, true>(in, out); break;
                        case 8: HelmholtzTriImpl<4, 4, 8, 7, true>(in, out); break;
                        default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                        } break;
                    case 5:
                        switch(nq0)
                        {
                        case 6: HelmholtzTriImpl<5, 5, 6, 5, true>(in, out); break;
                        case 7: HelmholtzTriImpl<5, 5, 7, 6, true>(in, out); break;
                        case 8: HelmholtzTriImpl<5, 5, 8, 7, true>(in, out); break;
                        case 9: HelmholtzTriImpl<5, 5, 9, 8, true>(in, out); break;
                        case 10: HelmholtzTriImpl<5, 5, 10, 9, true>(in, out); break;
                        default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                        } break;
                    case 6:
                        switch(nq0)
                        {
                        case 7: HelmholtzTriImpl<6, 6, 7, 6, true>(in, out); break;
                        case 8: HelmholtzTriImpl<6, 6, 8, 7, true>(in, out); break;
                        case 9: HelmholtzTriImpl<6, 6, 9, 8, true>(in, out); break;
                        case 10: HelmholtzTriImpl<6, 6, 10, 9, true>(in, out); break;
                        case 11: HelmholtzTriImpl<6, 6, 11, 10, true>(in, out); break;
                        case 12: HelmholtzTriImpl<6, 6, 12, 11, true>(in, out); break;
                        default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                        } break;
                    case 7:
                        switch(nq0)
                        {
                        case 8: HelmholtzTriImpl<7, 7, 8, 7, true>(in, out); break;
                        case 9: HelmholtzTriImpl<7, 7, 9, 8, true>(in, out); break;
                        case 10: HelmholtzTriImpl<7, 7, 10, 9, true>(in, out); break;
                        case 11: HelmholtzTriImpl<7, 7, 11, 10, true>(in, out); break;
                        case 12: HelmholtzTriImpl<7, 7, 12, 11, true>(in, out); break;
                        case 13: HelmholtzTriImpl<7, 7, 13, 12, true>(in, out); break;
                        case 14: HelmholtzTriImpl<7, 7, 14, 13, true>(in, out); break;
                        default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                        } break;
                    case 8:
                        switch(nq0)
                        {
                        case 9: HelmholtzTriImpl<8, 8, 9, 8, true>(in, out); break;
                        case 10: HelmholtzTriImpl<8, 8, 10, 9, true>(in, out); break;
                        case 11: HelmholtzTriImpl<8, 8, 11, 10, true>(in, out); break;
                        case 12: HelmholtzTriImpl<8, 8, 12, 11, true>(in, out); break;
                        case 13: HelmholtzTriImpl<8, 8, 13, 12, true>(in, out); break;
                        case 14: HelmholtzTriImpl<8, 8, 14, 13, true>(in, out); break;
                        case 15: HelmholtzTriImpl<8, 8, 15, 14, true>(in, out); break;
                        case 16: HelmholtzTriImpl<8, 8, 16, 15, true>(in, out); break;
                        default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                        } break;
                    default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,true,in, out); break;
                    }
                }
                else
                {
                    switch(nm0)
                    {
                    case 2:
                        switch(nq0)
                        {
                        case 3: HelmholtzTriImpl<2 ,2 ,3 ,2, false >(in, out); break;
                        case 4: HelmholtzTriImpl<2 ,2 ,4 ,3, false >(in, out); break;
                        default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,false,in, out); break;
                        } break;
                    case 3:
                        switch(nq0)
                        {
                        case 4: HelmholtzTriImpl<3 ,3 ,4 ,3, false >(in, out); break;
                        case 5: HelmholtzTriImpl<3 ,3 ,5 ,4, false >(in, out); break;
                        case 6: HelmholtzTriImpl<3 ,3 ,6 ,5, false >(in, out); break;
                        default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,false,in, out); break;
                    } break;
                    case 4:
                        switch(nq0)
                        {
                        case 5: HelmholtzTriImpl<4, 4, 5, 4, false>(in, out); break;
                        case 6: HelmholtzTriImpl<4, 4, 6, 5, false>(in, out); break;
                        case 7: HelmholtzTriImpl<4, 4, 7, 6, false>(in, out); break;
                        case 8: HelmholtzTriImpl<4, 4, 8, 7, false>(in, out); break;
                        default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,false,in, out); break;
                        } break;
                    case 5:
                        switch(nq0)
                        {
                        case 6: HelmholtzTriImpl<5, 5, 6, 5, false>(in, out); break;
                        case 7: HelmholtzTriImpl<5, 5, 7, 6, false>(in, out); break;
                        case 8: HelmholtzTriImpl<5, 5, 8, 7, false>(in, out); break;
                        case 9: HelmholtzTriImpl<5, 5, 9, 8, false>(in, out); break;
                        case 10: HelmholtzTriImpl<5, 5, 10, 9, false>(in, out); break;
                        default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,false,in, out); break;
                        } break;
                    case 6:
                        switch(nq0)
                        {
                        case 7: HelmholtzTriImpl<6, 6, 7, 6, false>(in, out); break;
                        case 8: HelmholtzTriImpl<6, 6, 8, 7, false>(in, out); break;
                        case 9: HelmholtzTriImpl<6, 6, 9, 8, false>(in, out); break;
                        case 10: HelmholtzTriImpl<6, 6, 10, 9, false>(in, out); break;
                        case 11: HelmholtzTriImpl<6, 6, 11, 10, false>(in, out); break;
                        case 12: HelmholtzTriImpl<6, 6, 12, 11, false>(in, out); break;
                        default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,false,in, out); break;
                        } break;
                    case 7:
                        switch(nq0)
                        {
                        case 8: HelmholtzTriImpl<7, 7, 8, 7, false>(in, out); break;
                        case 9: HelmholtzTriImpl<7, 7, 9, 8, false>(in, out); break;
                        case 10: HelmholtzTriImpl<7, 7, 10, 9, false>(in, out); break;
                        case 11: HelmholtzTriImpl<7, 7, 11, 10, false>(in, out); break;
                        case 12: HelmholtzTriImpl<7, 7, 12, 11, false>(in, out); break;
                        case 13: HelmholtzTriImpl<7, 7, 13, 12, false>(in, out); break;
                        case 14: HelmholtzTriImpl<7, 7, 14, 13, false>(in, out); break;
                        default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,false,in, out); break;
                        } break;
                    case 8:
                        switch(nq0)
                        {
                        case 9: HelmholtzTriImpl<8, 8, 9, 8, false>(in, out); break;
                        case 10: HelmholtzTriImpl<8, 8, 10, 9, false>(in, out); break;
                        case 11: HelmholtzTriImpl<8, 8, 11, 10, false>(in, out); break;
                        case 12: HelmholtzTriImpl<8, 8, 12, 11, false>(in, out); break;
                        case 13: HelmholtzTriImpl<8, 8, 13, 12, false>(in, out); break;
                        case 14: HelmholtzTriImpl<8, 8, 14, 13, false>(in, out); break;
                        case 15: HelmholtzTriImpl<8, 8, 15, 14, false>(in, out); break;
                        case 16: HelmholtzTriImpl<8, 8, 16, 15, false>(in, out); break;
                        default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,false,in, out); break;
                        } break;
                    default: HelmholtzTriImpl(nm0,nm1,nq0,nq1,false,in, out); break;
                    }
                }
            }
            else
            {
                if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
                {
                    HelmholtzTriImpl(nm0,nm1,nq0,nq1,true,in, out);
                }
                else
                {
                    HelmholtzTriImpl(nm0,nm1,nq0,nq1,false,in, out);
                }
            }
        }
        
        template<int NM0, int NM1, int NQ0, int NQ1, bool CORRECT>
        void HelmholtzTriImpl(const Array<OneD, const NekDouble> &input,
                              Array<OneD,NekDouble> &out)
        {
            auto *inptr  = &input[0];
            auto *outptr = &out[0];

            constexpr auto ndf = 4;
            constexpr auto nqTot = NQ0 * NQ1;
            const auto nmBlocks = m_nmTot * vec_t::width;


            // Allocate sufficient workspace for backwards transform and inner
            // product kernels.
            constexpr auto wspInnerProd = NQ1;
            constexpr auto wspBwdTrans = NM0;
            constexpr auto wspSize = wspInnerProd > wspBwdTrans ? wspInnerProd : wspBwdTrans;

            vec_t wsp[wspSize]; // workspace for kernels

            std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(m_nmTot);
            std::vector<vec_t, allocator<vec_t>> bwd(nqTot),  deriv0(nqTot), deriv1(nqTot);

            const vec_t* jac_ptr;
            const vec_t* df_ptr;

            vec_t df0,df1,df2,df3;
            vec_t metric00,metric01,metric11; 
        
            vec_t d00 = {1.0};
            vec_t d01 = {0.0};
            vec_t d11 = {1.0}; // var diffusion terms
            vec_t dtmp0,dtmp1,dtmp2,dtmp3; // metric products
            boost::ignore_unused(d00,d01,d11,dtmp0,dtmp1,dtmp2,dtmp3);
        
            if (this->m_isConstVarDiff)
            {
                d00 = this->m_constVarDiff[0];
                d01 = this->m_constVarDiff[1];
                d11 = this->m_constVarDiff[2];
            }
        
            // Get size of derivative factor block
            auto dfSize = ndf;
            if (DEFORMED)
            {
                dfSize *= nqTot;
            }

            for (int e = 0; e < this->m_nBlocks; e++)
            {
                df_ptr = &((*this->m_df)[e*dfSize]);

                // Load and transpose data
                load_interleave(inptr, m_nmTot, tmpIn);

                // Precompute Laplacian metricsp
                if(!DEFORMED)
                {
                    df0 = df_ptr[0];  df1 = df_ptr[1];
                    df2 = df_ptr[2];  df3 = df_ptr[3];

                    jac_ptr = &((*this->m_jac)[e]);
                }
                else
                {
                    jac_ptr = &((*this->m_jac)[e*nqTot]);
                }

                // Step 1: BwdTrans
                BwdTransTriKernel(NM0, NM1, NQ0, NQ1, CORRECT,
                  tmpIn, this->m_bdata[0], this->m_bdata[1], wsp, bwd);

                // Step 2: inner product for mass matrix operation
                IProductTriKernel<NM0, NM1, NQ0, NQ1, CORRECT, true, false,
                                  DEFORMED>
                    (bwd, this->m_bdata[0], this->m_bdata[1], this->m_w[0],
                     this->m_w[1], jac_ptr, wsp, tmpOut, m_lambda);
                
                // Step 3: take derivatives in collapsed coordinate space
                PhysDerivTensor2DKernel(NQ0, NQ1, 
                                        bwd, this->m_D[0], this->m_D[1], deriv0, deriv1);
                
                // Step 4a: Construct Laplacian metrics
                for (size_t j = 0, cnt = 0; j < NQ1; ++j)
                {
                    vec_t h1j = m_h1[j];
                    for (size_t i = 0; i < NQ0; ++i, ++cnt)
                    {
                        if(DEFORMED)
                        {
                            df0 = df_ptr[cnt * ndf];
                            df1 = df_ptr[cnt * ndf + 1];
                            df2 = df_ptr[cnt * ndf + 2];
                            df3 = df_ptr[cnt * ndf + 3];
                        }
                        
                        vec_t h0i = m_h0[i];
                        
                        // M = [M_00, df1; M_10; df3]
                        metric00 = h1j * (df0 + h0i * df1);  // M_00
                        vec_t tmp = h1j * (df2 + h0i * df3); // M_10
                        
                        if (!this->m_isConstVarDiff && !this->m_isVarDiff)
                        {
                            metric01 = metric00 * df1;
                            metric00 = metric00 * metric00;
                            
                            metric01.fma(tmp, df3);
                            metric00.fma(tmp, tmp);
                            
                            metric11 = df1 * df1;
                            metric11.fma(df3, df3);
                        }
                        else
                        {
                            if (this->m_isVarDiff)
                            {
                                d00 = m_varD00[cnt];
                                d01 = m_varD01[cnt];
                                d11 = m_varD11[cnt];
                            }
                            // M = [M_00, df1; M_10; df3]
                            dtmp0 = metric00 * d00;
                            dtmp0.fma(tmp,d01);
                            dtmp1 = metric00 * d01;
                            dtmp1.fma(tmp,d11);
                            dtmp2 = df1 * d00;
                            dtmp2.fma(df3,d01);
                            dtmp3 = df1 * d01;
                            dtmp3.fma(df3,d11);
                            
                            metric00 = metric00 * dtmp0;
                            metric00.fma(tmp,dtmp1);
                            metric01 = df1 * dtmp0;
                            metric01.fma(df3,dtmp1);
                            metric11 = df1 * dtmp2;
                            metric11.fma(df3,dtmp3);
                        }
                        
                        vec_t d0 = deriv0[cnt];
                        vec_t d1 = deriv1[cnt];
                        
                        tmp = metric00 * d0;
                        tmp.fma(metric01, d1);
                        bwd[cnt]  = tmp;
                        
                        tmp = metric01 * d0;
                        tmp.fma(metric11, d1);
                        deriv0[cnt] = tmp;
                    }
                }
                                  
                // Step 4b: Take inner products
                IProductTriKernel<NM0, NM1, NQ0, NQ1, CORRECT,
                                  false, true, DEFORMED>
                    (bwd, this->m_dbdata[0], this->m_bdata[1], this->m_w[0],
                     this->m_w[1], jac_ptr, wsp, tmpOut);
                
                IProductTriKernel<NM0, NM1, NQ0, NQ1, CORRECT,
                                  false, true, DEFORMED>
                    (deriv0, this->m_bdata[0], this->m_dbdata[1], this->m_w[0],
                     this->m_w[1], jac_ptr, wsp, tmpOut);

                // de-interleave and store data
                deinterleave_store(tmpOut, m_nmTot, outptr);
                
                inptr  += nmBlocks;
                outptr += nmBlocks;
            }
        }
    
        void HelmholtzTriImpl(
          const int nm0, const int nm1,
          const int nq0, const int nq1,
          const bool CORRECT,
          const Array<OneD, const NekDouble> &input, Array<OneD,NekDouble> &out)
        {
            auto *inptr  = &input[0];
            auto *outptr = &out[0];

            constexpr auto ndf = 4;
            const auto nqTot = nq0 * nq1;
            const auto nmBlocks = m_nmTot * vec_t::width;


            // Allocate sufficient workspace for backwards transform and inner
            // product kernels.
            const auto wspInnerProd = nq1;
            const auto wspBwdTrans = nm0;
            const auto wspSize = wspInnerProd > wspBwdTrans ? wspInnerProd : wspBwdTrans;

            std::vector<vec_t, allocator<vec_t>> wsp(wspSize), tmpIn(m_nmTot),
                tmpOut(m_nmTot);
            std::vector<vec_t, allocator<vec_t>> bwd(nqTot),  deriv0(nqTot),
                deriv1(nqTot);

            const vec_t* jac_ptr;
            const vec_t* df_ptr;

            vec_t df0,df1,df2,df3;
            vec_t metric00,metric01,metric11; 
        
            vec_t d00 = {1.0};
            vec_t d01 = {0.0};
            vec_t d11 = {1.0}; // var diffusion terms
            vec_t dtmp0,dtmp1,dtmp2,dtmp3; // metric products
            boost::ignore_unused(d00,d01,d11,dtmp0,dtmp1,dtmp2,dtmp3);
        
            if (this->m_isConstVarDiff)
            {
                d00 = this->m_constVarDiff[0];
                d01 = this->m_constVarDiff[1];
                d11 = this->m_constVarDiff[2];
            }
        
            // Get size of derivative factor block
            auto dfSize = ndf;
            if (DEFORMED)
            {
                dfSize *= nqTot;
            }

            for (int e = 0; e < this->m_nBlocks; e++)
            {
                df_ptr = &((*this->m_df)[e*dfSize]);

                // Load and transpose data
                load_interleave(inptr, m_nmTot, tmpIn);

                // Precompute Laplacian metricsp
                if(!DEFORMED)
                {
                    df0 = df_ptr[0];  df1 = df_ptr[1];
                    df2 = df_ptr[2];  df3 = df_ptr[3];

                    jac_ptr = &((*this->m_jac)[e]);
                }
                else
                {
                    jac_ptr = &((*this->m_jac)[e*nqTot]);
                }

                // Step 1: BwdTrans
                BwdTransTriKernel(nm0,nm1,nq0,nq1, CORRECT,
                  tmpIn, this->m_bdata[0], this->m_bdata[1], &wsp[0], bwd);

                // Step 2: inner product for mass matrix operation
                IProductTriKernel(nm0,nm1,nq0,nq1, CORRECT, true, false,
                                  DEFORMED,
                     bwd, this->m_bdata[0], this->m_bdata[1], this->m_w[0],
                     this->m_w[1], jac_ptr, &wsp[0], tmpOut, m_lambda);
                
                // Step 3: take derivatives in collapsed coordinate space
                PhysDerivTensor2DKernel(nq0,nq1, bwd, this->m_D[0],
                                        this->m_D[1], deriv0, deriv1);
                
                // Step 4a: Construct Laplacian metrics
                for (size_t j = 0, cnt = 0; j < nq1; ++j)
                {
                    vec_t h1j = m_h1[j];
                    for (size_t i = 0; i < nq0; ++i, ++cnt)
                    {
                        if(DEFORMED)
                        {
                            df0 = df_ptr[cnt * ndf];
                            df1 = df_ptr[cnt * ndf + 1];
                            df2 = df_ptr[cnt * ndf + 2];
                            df3 = df_ptr[cnt * ndf + 3];
                        }
                        
                        vec_t h0i = m_h0[i];
                        
                        // M = [M_00, df1; M_10; df3]
                        metric00 = h1j * (df0 + h0i * df1);  // M_00
                        vec_t tmp = h1j * (df2 + h0i * df3); // M_10
                        
                        if (!this->m_isConstVarDiff && !this->m_isVarDiff)
                        {
                            metric01 = metric00 * df1;
                            metric00 = metric00 * metric00;
                            
                            metric01.fma(tmp, df3);
                            metric00.fma(tmp, tmp);
                            
                            metric11 = df1 * df1;
                            metric11.fma(df3, df3);
                        }
                        else
                        {
                            if (this->m_isVarDiff)
                            {
                                d00 = m_varD00[cnt];
                                d01 = m_varD01[cnt];
                                d11 = m_varD11[cnt];
                            }
                            // M = [M_00, df1; M_10; df3]
                            dtmp0 = metric00 * d00;
                            dtmp0.fma(tmp,d01);
                            dtmp1 = metric00 * d01;
                            dtmp1.fma(tmp,d11);
                            dtmp2 = df1 * d00;
                            dtmp2.fma(df3,d01);
                            dtmp3 = df1 * d01;
                            dtmp3.fma(df3,d11);
                            
                            metric00 = metric00 * dtmp0;
                            metric00.fma(tmp,dtmp1);
                            metric01 = df1 * dtmp0;
                            metric01.fma(df3,dtmp1);
                            metric11 = df1 * dtmp2;
                            metric11.fma(df3,dtmp3);
                        }
                        
                        vec_t d0 = deriv0[cnt];
                        vec_t d1 = deriv1[cnt];
                        
                        tmp = metric00 * d0;
                        tmp.fma(metric01, d1);
                        bwd[cnt]  = tmp;
                        
                        tmp = metric01 * d0;
                        tmp.fma(metric11, d1);
                        deriv0[cnt] = tmp;
                    }
                }
                                  
                // Step 4b: Take inner products
                IProductTriKernel(nm0,nm1,nq0,nq1, CORRECT,
                                  false, true, DEFORMED,
                bwd, this->m_dbdata[0], this->m_bdata[1], this->m_w[0],
                this->m_w[1], jac_ptr, &wsp[0], tmpOut);
                
                IProductTriKernel(nm0,nm1,nq0,nq1, CORRECT,
                                  false, true, DEFORMED,
                deriv0, this->m_bdata[0], this->m_dbdata[1], this->m_w[0],
                this->m_w[1], jac_ptr, &wsp[0], tmpOut);

                // de-interleave and store data
                deinterleave_store(tmpOut, m_nmTot, outptr);
                
                inptr  += nmBlocks;
                outptr += nmBlocks;
            }
        }
        
    
    public:
        

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }
    
private:
    const int m_nmTot;
    std::vector<vec_t, allocator<vec_t>> m_h0;
    std::vector<vec_t, allocator<vec_t>> m_h1;
};

template<bool DEFORMED = false>
struct HelmholtzHex : public Helmholtz, public Helper<3, DEFORMED>
{
    HelmholtzHex(std::vector<LibUtilities::BasisSharedPtr> basis,
                    int nElmt)
    : Helmholtz(basis, nElmt),
      Helper<3, DEFORMED>(basis, nElmt),
      m_nmTot(LibUtilities::StdHexData::getNumberOfCoefficients(
                  this->m_nm[0], this->m_nm[1], this->m_nm[2]))
    {
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<HelmholtzHex<DEFORMED>>(basis, nElmt);
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out) final
    {

        const int nm0 = m_basis[0]->GetNumModes();  
        const int nm1 = m_basis[1]->GetNumModes();  
        const int nm2 = m_basis[2]->GetNumModes();  
        const int nq0 = m_basis[0]->GetNumPoints();  
        const int nq1 = m_basis[1]->GetNumPoints();  
        const int nq2 = m_basis[2]->GetNumPoints();  

        if((nm0 == nm1)&&(nm0 == nm2)&&(nq0 == nq1)&&(nq0 == nq2))
        {
            switch(nm0)
            {
            case 2:
                switch(nq0)
                {
                case 2: HelmholtzHexImpl<2, 2, 2, 2, 2, 2>(in, out); break;
                case 3: HelmholtzHexImpl<2, 2, 2, 3, 3, 3>(in, out); break;
                case 4: HelmholtzHexImpl<2, 2, 2, 4, 4, 4>(in, out); break;
                default: HelmholtzHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
                } break;
            case 3:
                switch(nq0)
                {
                case 3: HelmholtzHexImpl<3, 3, 3, 3, 3, 3>(in, out); break;
                case 4: HelmholtzHexImpl<3, 3, 3, 4, 4, 4>(in, out); break;
                case 5: HelmholtzHexImpl<3, 3, 3, 5, 5, 5>(in, out); break;
                case 6: HelmholtzHexImpl<3, 3, 3, 6, 6, 6>(in, out); break;
                default: HelmholtzHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
                } break;
            case 4:
                switch(nq0)
                {
                case 4: HelmholtzHexImpl<4, 4, 4, 4, 4, 4>(in, out); break;
                case 5: HelmholtzHexImpl<4, 4, 4, 5, 5, 5>(in, out); break;
                case 6: HelmholtzHexImpl<4, 4, 4, 6, 6, 6>(in, out); break;
                case 7: HelmholtzHexImpl<4, 4, 4, 7, 7, 7>(in, out); break;
                case 8: HelmholtzHexImpl<4, 4, 4, 8, 8, 8>(in, out); break;
                default: HelmholtzHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
                } break;
            case 5:
                switch(nq0)
                {
                case 5: HelmholtzHexImpl<5, 5, 5, 5, 5, 5>(in, out); break;
                case 6: HelmholtzHexImpl<5, 5, 5, 6, 6, 6>(in, out); break;
                case 7: HelmholtzHexImpl<5, 5, 5, 7, 7, 7>(in, out); break;
                case 8: HelmholtzHexImpl<5, 5, 5, 8, 8, 8>(in, out); break;
                case 9: HelmholtzHexImpl<5, 5, 5, 9, 9, 9>(in, out); break;
                case 10: HelmholtzHexImpl<5, 5, 5, 10, 10, 10>(in, out); break;
                default: HelmholtzHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
                } break;
            case 6:
                switch(nq0)
                {
                case 6: HelmholtzHexImpl<6, 6, 6, 6, 6, 6>(in, out); break;
                case 7: HelmholtzHexImpl<6, 6, 6, 7, 7, 7>(in, out); break;
                case 8: HelmholtzHexImpl<6, 6, 6, 8, 8, 8>(in, out); break;
                case 9: HelmholtzHexImpl<6, 6, 6, 9, 9, 9>(in, out); break;
                case 10: HelmholtzHexImpl<6, 6, 6, 10, 10, 10>(in, out); break;
                case 11: HelmholtzHexImpl<6, 6, 6, 11, 11, 11>(in, out); break;
                case 12: HelmholtzHexImpl<6, 6, 6, 12, 12, 12>(in, out); break;
                default: HelmholtzHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
                } break;
            case 7:
                switch(nq0)
                {
                case 7: HelmholtzHexImpl<7, 7, 7, 7, 7, 7>(in, out); break;
                case 8: HelmholtzHexImpl<7, 7, 7, 8, 8, 8>(in, out); break;
                case 9: HelmholtzHexImpl<7, 7, 7, 9, 9, 9>(in, out); break;
                case 10: HelmholtzHexImpl<7, 7, 7, 10, 10, 10>(in, out); break;
                case 11: HelmholtzHexImpl<7, 7, 7, 11, 11, 11>(in, out); break;
                case 12: HelmholtzHexImpl<7, 7, 7, 12, 12, 12>(in, out); break;
                case 13: HelmholtzHexImpl<7, 7, 7, 13, 13, 13>(in, out); break;
                case 14: HelmholtzHexImpl<7, 7, 7, 14, 14, 14>(in, out); break;
                default: HelmholtzHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
                } break;
            case 8:
                switch(nq0)
                {
                case 8: HelmholtzHexImpl<8, 8, 8, 8, 8, 8>(in, out); break;
                case 9: HelmholtzHexImpl<8, 8, 8, 9, 9, 9>(in, out); break;
                case 10: HelmholtzHexImpl<8, 8, 8, 10, 10, 10>(in, out); break;
                case 11: HelmholtzHexImpl<8, 8, 8, 11, 11, 11>(in, out); break;
                case 12: HelmholtzHexImpl<8, 8, 8, 12, 12, 12>(in, out); break;
                case 13: HelmholtzHexImpl<8, 8, 8, 13, 13, 13>(in, out); break;
                case 14: HelmholtzHexImpl<8, 8, 8, 14, 14, 14>(in, out); break;
                case 15: HelmholtzHexImpl<8, 8, 8, 15, 15, 15>(in, out); break;
                case 16: HelmholtzHexImpl<8, 8, 8, 16, 16, 16>(in, out); break;
                default: HelmholtzHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
                } break;
            default: HelmholtzHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out); break;
            }
        }
        else
        {
            HelmholtzHexImpl(nm0,nm1,nm2,nq0,nq1,nq2,in,out);
        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2>
    void HelmholtzHexImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out)
    {
        auto *inptr  = &input[0];
        auto *outptr = &out[0];

        constexpr auto ndf = 9;
        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        const auto nmBlocks = m_nmTot * vec_t::width;

        // Workspace for kernels
        vec_t wsp1[NQ0 * NQ1 * NQ2], wsp2[NQ0 * NQ1 * NQ2];
        

        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot),  deriv0(nqTot),
            deriv1(nqTot), deriv2(nqTot);

        const vec_t* jac_ptr;
        const vec_t* df_ptr;

        vec_t df0,df1,df2,df3,df4,df5,df6,df7,df8;
        vec_t metric00,metric01,metric02,metric11,metric12, metric22;
        vec_t d00 = {1.0};
        vec_t d01 = {0.0};
        vec_t d11 = {1.0};
        vec_t d02 = {0.0};
        vec_t d12 = {0.0};
        vec_t d22 = {1.0}; // var diffusion terms
        vec_t td0,td1,td2,td3,td4,td5,td6,td7,td8; // temp terms for vardiff
        boost::ignore_unused(d00,d01,d11,d02,d12,d22,
                             td0,td1,td2,td3,td4,td5,td6,td7,td8);
        
        if (this->m_isConstVarDiff)
        {
            d00 = this->m_constVarDiff[0];
            d01 = this->m_constVarDiff[1];
            d11 = this->m_constVarDiff[2];
            d02 = this->m_constVarDiff[3];
            d12 = this->m_constVarDiff[4];
            d22 = this->m_constVarDiff[5];
        } 
        
        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            // Precompute Laplacian metrics
            df_ptr = &((*this->m_df)[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            // Precompute Laplacian metricsp
            if(!DEFORMED)
            {
                df0 = df_ptr[0];  df1 = df_ptr[1]; df2 = df_ptr[2];
                df3 = df_ptr[3];  df4 = df_ptr[4]; df5 = df_ptr[5];
                df6 = df_ptr[6];  df7 = df_ptr[7]; df8 = df_ptr[8];
                
                if (!this->m_isConstVarDiff && !this->m_isConstVarDiff)
                { 
                    metric00 = df0 * df0;
                    metric00.fma(df3, df3);
                    metric00.fma(df6, df6);
                    
                    metric01 = df0 * df1;
                    metric01.fma(df3, df4);
                    metric01.fma(df6, df7);
                    
                    metric02 = df0 * df2;
                    metric02.fma(df3, df5);
                    metric02.fma(df6, df8);
                    
                    metric11 = df1 * df1;
                    metric11.fma(df4, df4);
                    metric11.fma(df7, df7);
                    
                    metric12 = df1 * df2;
                    metric12.fma(df4, df5);
                    metric12.fma(df7, df8);
                    
                    metric22 = df2 * df2;
                    metric22.fma(df5, df5);
                    metric22.fma(df8, df8);
                }
                else if (this->m_isConstVarDiff)
                {   // with vardiff
                    td0 = df0 * d00;
                    td0.fma(df3,d01);
                    td0.fma(df6,d02);
                    
                    td1 = df0 * d01;
                    td1.fma(df3,d11);
                    td1.fma(df6,d12);
                    
                    td2 = df0 * d02;
                    td2.fma(df3,d12);
                    td2.fma(df6,d22);
                
                    td3 = df1 * d00;
                    td3.fma(df4,d01);
                    td3.fma(df7,d02);
                    
                    td4 = df1 * d01;
                    td4.fma(df4,d11);
                    td4.fma(df7,d12);
                    
                    td5 = df1 * d02;
                    td5.fma(df4,d12);
                    td5.fma(df7,d22);
                                        
                    td6 = df2 * d00;
                    td6.fma(df5,d01);
                    td6.fma(df8,d02);
                    
                    td7 = df2 * d01;
                    td7.fma(df5,d11);
                    td7.fma(df8,d12);
                    
                    td8 = df2 * d02;
                    td8.fma(df5,d12);
                    td8.fma(df8,d22);
                    
                    metric00 = td0 * df0;
                    metric00.fma(td1,df3);
                    metric00.fma(td2,df6);
                    
                    metric01 = td0 * df1;
                    metric01.fma(td1,df4);
                    metric01.fma(td2,df7);
                    
                    metric02 = td0 * df2;
                    metric02.fma(td1,df5);
                    metric02.fma(td2,df8);
                    
                    metric11 = td3 * df1;
                    metric11.fma(td4,df4);
                    metric11.fma(td5,df7);
                    
                    metric12 = td3 * df2;
                    metric12.fma(td4,df5);
                    metric12.fma(td5,df8);
                    
                    metric22 = td6 * df2;
                    metric22.fma(td7,df5);
                    metric22.fma(td8,df8);
                    
                }
                
                jac_ptr = &((*this->m_jac)[e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e*nqTot]);
            }

            // Step 1: BwdTrans
            BwdTransHexKernel(NM0, NM1, NM2, NQ0, NQ1, NQ2,
                              tmpIn, this->m_bdata[0],
                              this->m_bdata[0], this->m_bdata[0],
                              wsp1, wsp2, bwd);

            // Step 2: inner product for mass matrix operation
            IProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, true, false, DEFORMED>
                (bwd, this->m_bdata[0], this->m_bdata[0], this->m_bdata[0],
                 this->m_w[0], this->m_w[0], this->m_w[0], jac_ptr,
                 wsp1, wsp2, tmpOut, m_lambda);

            // Step 3: take derivatives in standard space
            PhysDerivTensor3DKernel(NQ0, NQ1, NQ2,
                 bwd, this->m_D[0], this->m_D[0], this->m_D[0],
                 deriv0, deriv1, deriv2);

            // Step 4: Apply Laplacian metrics & inner product
            if (DEFORMED || this->m_isVarDiff)
            {
                for (size_t k = 0,cnt=0; k < NQ2; k++)
                {
                    for (size_t j = 0; j < NQ1; j++)
                    {
                        for (size_t i = 0; i < NQ0; i++, ++cnt)
                        {
                            if (DEFORMED)
                            {
                                df0 = df_ptr[cnt * ndf];
                                df1 = df_ptr[cnt * ndf + 1];
                                df2 = df_ptr[cnt * ndf + 2];
                                df3 = df_ptr[cnt * ndf + 3];
                                df4 = df_ptr[cnt * ndf + 4];
                                df5 = df_ptr[cnt * ndf + 5];
                                df6 = df_ptr[cnt * ndf + 6];
                                df7 = df_ptr[cnt * ndf + 7];
                                df8 = df_ptr[cnt * ndf + 8];
                            }
                            
                            if (!this->m_isConstVarDiff && !this->m_isVarDiff)
                            {
                                metric00 = df0 * df0;
                                metric00.fma(df3, df3);
                                metric00.fma(df6, df6);

                                metric01 = df0 * df1;
                                metric01.fma(df3, df4);
                                metric01.fma(df6, df7);

                                metric02 = df0 * df2;
                                metric02.fma(df3, df5);
                                metric02.fma(df6, df8);

                                metric11 = df1 * df1;
                                metric11.fma(df4, df4);
                                metric11.fma(df7, df7);

                                metric12 = df1 * df2;
                                metric12.fma(df4, df5);
                                metric12.fma(df7, df8);

                                metric22 = df2 * df2;
                                metric22.fma(df5, df5);
                                metric22.fma(df8, df8);
                            }
                            else
                            {   // with vardiff
                                
                                if (this->m_isVarDiff)
                                {
                                    d00 = m_varD00[cnt];
                                    d01 = m_varD01[cnt];
                                    d11 = m_varD11[cnt];
                                    d02 = m_varD02[cnt];
                                    d12 = m_varD12[cnt];
                                    d22 = m_varD22[cnt];
                                }
                                
                                td0 = df0 * d00;
                                td0.fma(df3,d01);
                                td0.fma(df6,d02);
                                
                                td1 = df0 * d01;
                                td1.fma(df3,d11);
                                td1.fma(df6,d12);
                                
                                td2 = df0 * d02;
                                td2.fma(df3,d12);
                                td2.fma(df6,d22);
                            
                                td3 = df1 * d00;
                                td3.fma(df4,d01);
                                td3.fma(df7,d02);
                                
                                td4 = df1 * d01;
                                td4.fma(df4,d11);
                                td4.fma(df7,d12);
                                
                                td5 = df1 * d02;
                                td5.fma(df4,d12);
                                td5.fma(df7,d22);
                                                    
                                td6 = df2 * d00;
                                td6.fma(df5,d01);
                                td6.fma(df8,d02);
                                
                                td7 = df2 * d01;
                                td7.fma(df5,d11);
                                td7.fma(df8,d12);
                                
                                td8 = df2 * d02;
                                td8.fma(df5,d12);
                                td8.fma(df8,d22);
                                
                                metric00 = td0 * df0;
                                metric00.fma(td1,df3);
                                metric00.fma(td2,df6);
                                
                                metric01 = td0 * df1;
                                metric01.fma(td1,df4);
                                metric01.fma(td2,df7);
                                
                                metric02 = td0 * df2;
                                metric02.fma(td1,df5);
                                metric02.fma(td2,df8);
                                
                                metric11 = td3 * df1;
                                metric11.fma(td4,df4);
                                metric11.fma(td5,df7);
                                
                                metric12 = td3 * df2;
                                metric12.fma(td4,df5);
                                metric12.fma(td5,df8);
                                
                                metric22 = td6 * df2;
                                metric22.fma(td7,df5);
                                metric22.fma(td8,df8);
                                
                            }

                            vec_t d0 = deriv0[cnt];
                            vec_t d1 = deriv1[cnt];
                            vec_t d2 = deriv2[cnt];

                            vec_t tmp = metric00 * d0;
                            tmp.fma(metric01, d1);
                            tmp.fma(metric02, d2);
                            deriv0[cnt] = tmp; 

                            tmp = metric01 * d0;
                            tmp.fma(metric11, d1);
                            tmp.fma(metric12, d2);
                            deriv1[cnt] = tmp;

                            tmp = metric02 * d0;
                            tmp.fma(metric12, d1);
                            tmp.fma(metric22, d2);
                            deriv2[cnt] = tmp; 
                        }
                    }
                }
            }
            else
            {
                /*
                metric00 = df0 * df0;
                metric00.fma(df3, df3);
                metric00.fma(df6, df6);
                
                metric01 = df0 * df1;
                metric01.fma(df3, df4);
                metric01.fma(df6, df7);

                metric02 = df0 * df2;
                metric02.fma(df3, df5);
                metric02.fma(df6, df8);
                            
                metric11 = df1 * df1;
                metric11.fma(df4, df4);
                metric11.fma(df7, df7);
                
                metric12 = df1 * df2;
                metric12.fma(df4, df5);
                metric12.fma(df7, df8);
                
                metric22 = df2 * df2;
                metric22.fma(df5, df5);
                metric22.fma(df8, df8);
                */
                
                for (int i = 0; i < nqTot; ++i)
                {
                    vec_t d0 = deriv0[i];
                    vec_t d1 = deriv1[i];
                    vec_t d2 = deriv2[i];
                    
                    vec_t tmp = metric00 * d0;
                    tmp.fma(metric01, d1);
                    tmp.fma(metric02, d2);
                    deriv0[i] = tmp; 
                    
                    tmp = metric01 * d0;
                    tmp.fma(metric11, d1);
                    tmp.fma(metric12, d2);
                    deriv1[i] = tmp;
                    
                    tmp = metric02 * d0;
                    tmp.fma(metric12, d1);
                    tmp.fma(metric22, d2);
                    deriv2[i] = tmp; 
                }
            }

            IProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, false, true, DEFORMED>
                (deriv0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 wsp1, wsp2, tmpOut);
            
            IProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, false, true, DEFORMED>
                (deriv1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 wsp1, wsp2, tmpOut);

            IProductHexKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, false, true, DEFORMED>
                (deriv2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 wsp1, wsp2, tmpOut);
            
            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr  += nmBlocks;
            outptr += nmBlocks;
        }
    }

    void HelmholtzHexImpl(
        const int nm0, const int nm1, const int nm2,
        const int nq0, const int nq1, const int nq2,
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out)
    {
        auto *inptr  = &input[0];
        auto *outptr = &out[0];

        constexpr auto ndf = 9;
        const auto nqTot = nq0 * nq1 * nq2;
        const auto nmBlocks = m_nmTot * vec_t::width;

        // Workspace for kernels
        std::vector<vec_t, allocator<vec_t>> wsp1(nqTot), wsp2(nqTot),
            tmpIn(m_nmTot), tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot),  deriv0(nqTot),
            deriv1(nqTot), deriv2(nqTot);

        const vec_t* jac_ptr;
        const vec_t* df_ptr;

        vec_t df0,df1,df2,df3,df4,df5,df6,df7,df8;
        vec_t metric00,metric01,metric02,metric11,metric12, metric22;
        vec_t d00 = {1.0};
        vec_t d01 = {0.0};
        vec_t d11 = {1.0};
        vec_t d02 = {0.0};
        vec_t d12 = {0.0};
        vec_t d22 = {1.0}; // var diffusion terms
        vec_t td0,td1,td2,td3,td4,td5,td6,td7,td8; // temp terms for vardiff
        boost::ignore_unused(d00,d01,d11,d02,d12,d22,
                             td0,td1,td2,td3,td4,td5,td6,td7,td8);
        
        if (this->m_isConstVarDiff)
        {
            d00 = this->m_constVarDiff[0];
            d01 = this->m_constVarDiff[1];
            d11 = this->m_constVarDiff[2];
            d02 = this->m_constVarDiff[3];
            d12 = this->m_constVarDiff[4];
            d22 = this->m_constVarDiff[5];
        } 
        
        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            // Precompute Laplacian metrics
            df_ptr = &((*this->m_df)[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            // Precompute Laplacian metricsp
            if(!DEFORMED)
            {
                df0 = df_ptr[0];  df1 = df_ptr[1]; df2 = df_ptr[2];
                df3 = df_ptr[3];  df4 = df_ptr[4]; df5 = df_ptr[5];
                df6 = df_ptr[6];  df7 = df_ptr[7]; df8 = df_ptr[8];
                
                if (!this->m_isConstVarDiff && !this->m_isConstVarDiff)
                { 
                    metric00 = df0 * df0;
                    metric00.fma(df3, df3);
                    metric00.fma(df6, df6);
                    
                    metric01 = df0 * df1;
                    metric01.fma(df3, df4);
                    metric01.fma(df6, df7);
                    
                    metric02 = df0 * df2;
                    metric02.fma(df3, df5);
                    metric02.fma(df6, df8);
                    
                    metric11 = df1 * df1;
                    metric11.fma(df4, df4);
                    metric11.fma(df7, df7);
                    
                    metric12 = df1 * df2;
                    metric12.fma(df4, df5);
                    metric12.fma(df7, df8);
                    
                    metric22 = df2 * df2;
                    metric22.fma(df5, df5);
                    metric22.fma(df8, df8);
                }
                else if (this->m_isConstVarDiff)
                {   // with vardiff
                    td0 = df0 * d00;
                    td0.fma(df3,d01);
                    td0.fma(df6,d02);
                    
                    td1 = df0 * d01;
                    td1.fma(df3,d11);
                    td1.fma(df6,d12);
                    
                    td2 = df0 * d02;
                    td2.fma(df3,d12);
                    td2.fma(df6,d22);
                
                    td3 = df1 * d00;
                    td3.fma(df4,d01);
                    td3.fma(df7,d02);
                    
                    td4 = df1 * d01;
                    td4.fma(df4,d11);
                    td4.fma(df7,d12);
                    
                    td5 = df1 * d02;
                    td5.fma(df4,d12);
                    td5.fma(df7,d22);
                                        
                    td6 = df2 * d00;
                    td6.fma(df5,d01);
                    td6.fma(df8,d02);
                    
                    td7 = df2 * d01;
                    td7.fma(df5,d11);
                    td7.fma(df8,d12);
                    
                    td8 = df2 * d02;
                    td8.fma(df5,d12);
                    td8.fma(df8,d22);
                    
                    metric00 = td0 * df0;
                    metric00.fma(td1,df3);
                    metric00.fma(td2,df6);
                    
                    metric01 = td0 * df1;
                    metric01.fma(td1,df4);
                    metric01.fma(td2,df7);
                    
                    metric02 = td0 * df2;
                    metric02.fma(td1,df5);
                    metric02.fma(td2,df8);
                    
                    metric11 = td3 * df1;
                    metric11.fma(td4,df4);
                    metric11.fma(td5,df7);
                    
                    metric12 = td3 * df2;
                    metric12.fma(td4,df5);
                    metric12.fma(td5,df8);
                    
                    metric22 = td6 * df2;
                    metric22.fma(td7,df5);
                    metric22.fma(td8,df8);
                    
                }
                
                jac_ptr = &((*this->m_jac)[e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e*nqTot]);
            }

            // Step 1: BwdTrans
            BwdTransHexKernel(nm0,nm1,nm2,nq0,nq1,nq2,
                              tmpIn, this->m_bdata[0],
                              this->m_bdata[0], this->m_bdata[0],
                              &wsp1[0], &wsp2[0], bwd);

            // Step 2: inner product for mass matrix operation
            IProductHexKernel(nm0,nm1,nm2,nq0,nq1,nq2, true, false, DEFORMED,
                 bwd, this->m_bdata[0], this->m_bdata[0], this->m_bdata[0],
                 this->m_w[0], this->m_w[0], this->m_w[0], jac_ptr,
                              &wsp1[0], &wsp2[0], tmpOut, m_lambda);

            // Step 3: take derivatives in standard space
            PhysDerivTensor3DKernel(nq0,nq1,nq2,
                 bwd, this->m_D[0], this->m_D[0], this->m_D[0],
                 deriv0, deriv1, deriv2);

            // Step 4: Apply Laplacian metrics & inner product
            if (DEFORMED || this->m_isVarDiff)
            {
                for (size_t k = 0,cnt=0; k < nq2; k++)
                {
                    for (size_t j = 0; j < nq1; j++)
                    {
                        for (size_t i = 0; i < nq0; i++, ++cnt)
                        {
                            if (DEFORMED)
                            {
                                df0 = df_ptr[cnt * ndf];
                                df1 = df_ptr[cnt * ndf + 1];
                                df2 = df_ptr[cnt * ndf + 2];
                                df3 = df_ptr[cnt * ndf + 3];
                                df4 = df_ptr[cnt * ndf + 4];
                                df5 = df_ptr[cnt * ndf + 5];
                                df6 = df_ptr[cnt * ndf + 6];
                                df7 = df_ptr[cnt * ndf + 7];
                                df8 = df_ptr[cnt * ndf + 8];
                            }
                            
                            if (!this->m_isConstVarDiff && !this->m_isVarDiff)
                            {
                                metric00 = df0 * df0;
                                metric00.fma(df3, df3);
                                metric00.fma(df6, df6);

                                metric01 = df0 * df1;
                                metric01.fma(df3, df4);
                                metric01.fma(df6, df7);

                                metric02 = df0 * df2;
                                metric02.fma(df3, df5);
                                metric02.fma(df6, df8);

                                metric11 = df1 * df1;
                                metric11.fma(df4, df4);
                                metric11.fma(df7, df7);

                                metric12 = df1 * df2;
                                metric12.fma(df4, df5);
                                metric12.fma(df7, df8);

                                metric22 = df2 * df2;
                                metric22.fma(df5, df5);
                                metric22.fma(df8, df8);
                            }
                            else
                            {   // with vardiff
                                
                                if (this->m_isVarDiff)
                                {
                                    d00 = m_varD00[cnt];
                                    d01 = m_varD01[cnt];
                                    d11 = m_varD11[cnt];
                                    d02 = m_varD02[cnt];
                                    d12 = m_varD12[cnt];
                                    d22 = m_varD22[cnt];
                                }
                                
                                td0 = df0 * d00;
                                td0.fma(df3,d01);
                                td0.fma(df6,d02);
                                
                                td1 = df0 * d01;
                                td1.fma(df3,d11);
                                td1.fma(df6,d12);
                                
                                td2 = df0 * d02;
                                td2.fma(df3,d12);
                                td2.fma(df6,d22);
                            
                                td3 = df1 * d00;
                                td3.fma(df4,d01);
                                td3.fma(df7,d02);
                                
                                td4 = df1 * d01;
                                td4.fma(df4,d11);
                                td4.fma(df7,d12);
                                
                                td5 = df1 * d02;
                                td5.fma(df4,d12);
                                td5.fma(df7,d22);
                                                    
                                td6 = df2 * d00;
                                td6.fma(df5,d01);
                                td6.fma(df8,d02);
                                
                                td7 = df2 * d01;
                                td7.fma(df5,d11);
                                td7.fma(df8,d12);
                                
                                td8 = df2 * d02;
                                td8.fma(df5,d12);
                                td8.fma(df8,d22);
                                
                                metric00 = td0 * df0;
                                metric00.fma(td1,df3);
                                metric00.fma(td2,df6);
                                
                                metric01 = td0 * df1;
                                metric01.fma(td1,df4);
                                metric01.fma(td2,df7);
                                
                                metric02 = td0 * df2;
                                metric02.fma(td1,df5);
                                metric02.fma(td2,df8);
                                
                                metric11 = td3 * df1;
                                metric11.fma(td4,df4);
                                metric11.fma(td5,df7);
                                
                                metric12 = td3 * df2;
                                metric12.fma(td4,df5);
                                metric12.fma(td5,df8);
                                
                                metric22 = td6 * df2;
                                metric22.fma(td7,df5);
                                metric22.fma(td8,df8);
                                
                            }

                            vec_t d0 = deriv0[cnt];
                            vec_t d1 = deriv1[cnt];
                            vec_t d2 = deriv2[cnt];

                            vec_t tmp = metric00 * d0;
                            tmp.fma(metric01, d1);
                            tmp.fma(metric02, d2);
                            deriv0[cnt] = tmp; 

                            tmp = metric01 * d0;
                            tmp.fma(metric11, d1);
                            tmp.fma(metric12, d2);
                            deriv1[cnt] = tmp;

                            tmp = metric02 * d0;
                            tmp.fma(metric12, d1);
                            tmp.fma(metric22, d2);
                            deriv2[cnt] = tmp; 
                        }
                    }
                }
            }
            else
            {
                /*
                metric00 = df0 * df0;
                metric00.fma(df3, df3);
                metric00.fma(df6, df6);
                
                metric01 = df0 * df1;
                metric01.fma(df3, df4);
                metric01.fma(df6, df7);

                metric02 = df0 * df2;
                metric02.fma(df3, df5);
                metric02.fma(df6, df8);
                            
                metric11 = df1 * df1;
                metric11.fma(df4, df4);
                metric11.fma(df7, df7);
                
                metric12 = df1 * df2;
                metric12.fma(df4, df5);
                metric12.fma(df7, df8);
                
                metric22 = df2 * df2;
                metric22.fma(df5, df5);
                metric22.fma(df8, df8);
                */
                
                for (int i = 0; i < nqTot; ++i)
                {
                    vec_t d0 = deriv0[i];
                    vec_t d1 = deriv1[i];
                    vec_t d2 = deriv2[i];
                    
                    vec_t tmp = metric00 * d0;
                    tmp.fma(metric01, d1);
                    tmp.fma(metric02, d2);
                    deriv0[i] = tmp; 
                    
                    tmp = metric01 * d0;
                    tmp.fma(metric11, d1);
                    tmp.fma(metric12, d2);
                    deriv1[i] = tmp;
                    
                    tmp = metric02 * d0;
                    tmp.fma(metric12, d1);
                    tmp.fma(metric22, d2);
                    deriv2[i] = tmp; 
                }
            }

            IProductHexKernel(nm0,nm1,nm2,nq0,nq1,nq2, false, true, DEFORMED,
                 deriv0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 &wsp1[0], &wsp2[0], tmpOut);
            
            IProductHexKernel(nm0,nm1,nm2,nq0,nq1,nq2, false, true, DEFORMED,
                 deriv1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                              &wsp1[0], &wsp2[0], tmpOut);

            IProductHexKernel(nm0,nm1,nm2,nq0,nq1,nq2, false, true, DEFORMED,
                 deriv2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                              &wsp1[0], &wsp2[0], tmpOut);
            
            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr  += nmBlocks;
            outptr += nmBlocks;
        }
    }
    

public:

private:
    const int m_nmTot;
};

template<bool DEFORMED = false>
struct HelmholtzPrism : public Helmholtz, public Helper<3, DEFORMED>
{
    HelmholtzPrism(std::vector<LibUtilities::BasisSharedPtr> basis,
                      int nElmt)
    : Helmholtz(basis, nElmt),
      Helper<3, DEFORMED>(basis, nElmt),
      m_nmTot(LibUtilities::StdPrismData::getNumberOfCoefficients(
                  this->m_nm[0], this->m_nm[1], this->m_nm[2])),
      m_h0(basis[0]->GetNumPoints()),
      m_h1(basis[2]->GetNumPoints())
    {
        const int nq0 = basis[0]->GetNumPoints();
        const int nq2 = basis[2]->GetNumPoints();

        const Array<OneD, const NekDouble> &z0 = basis[0]->GetZ();
        const Array<OneD, const NekDouble> &z2 = basis[2]->GetZ();

        for (int i = 0; i < nq0; ++i)
        {
            m_h0[i] = 0.5 * (1 + z0[i]);
        }

        for (int k = 0; k < nq2; ++k)
        {
            m_h1[k] = 2.0 / (1 - z2[k]);
        }
    }

    static std::shared_ptr<Operator> Create(
                                            std::vector<LibUtilities::BasisSharedPtr> basis,
                                            int nElmt)
    {
        return std::make_shared<HelmholtzPrism<DEFORMED>>(basis, nElmt);
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                            Array<OneD,       NekDouble> &out) final
    {
        const int nm0 = m_basis[0]->GetNumModes();  
        const int nm1 = m_basis[1]->GetNumModes();  
        const int nm2 = m_basis[2]->GetNumModes();  
        const int nq0 = m_basis[0]->GetNumPoints();  
        const int nq1 = m_basis[1]->GetNumPoints();  
        const int nq2 = m_basis[2]->GetNumPoints();  

        if((nm0 == nm1)&&(nm0 == nm2)&&(nq0 == nq1)&&(nq0 == nq2+1))
        {

            if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                switch(nm0)
                {
                case 2:
                    switch(nq0)
                    {
                    case 3: HelmholtzPrismImpl<2, 2, 2, 3, 3, 2, true>
                            (in, out); break;
                    case 4: HelmholtzPrismImpl<2, 2, 2, 4, 4, 3, true>
                            (in, out); break;
                    default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                                in, out); break;
                    } break;
                case 3:
                    switch(nq0)
                    {
                    case 4: HelmholtzPrismImpl<3, 3, 3, 4, 4, 3, true>
                            (in, out); break;
                    case 5: HelmholtzPrismImpl<3, 3, 3, 5, 5, 4, true>
                            (in, out); break;
                    case 6: HelmholtzPrismImpl<3, 3, 3, 6, 6, 5, true>
                            (in, out); break;
                    default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                                in, out); break;
                    } break;
                case 4:
                    switch(nq0)
                    {
                    case 5: HelmholtzPrismImpl<4, 4, 4, 5, 5, 4, true>
                            (in, out); break;
                    case 6: HelmholtzPrismImpl<4, 4, 4, 6, 6, 5, true>
                            (in, out); break;
                    case 7: HelmholtzPrismImpl<4, 4, 4, 7, 7, 6, true>
                            (in, out); break;
                    case 8: HelmholtzPrismImpl<4, 4, 4, 8, 8, 7, true>
                            (in, out); break;
                    default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                                in, out); break;
                    } break;
                case 5:
                    switch(nq0)
                    {
                    case 6: HelmholtzPrismImpl<5, 5, 5, 6, 6, 5, true>
                            (in, out); break;
                    case 7: HelmholtzPrismImpl<5, 5, 5, 7, 7, 6, true>
                            (in, out); break;
                    case 8: HelmholtzPrismImpl<5, 5, 5, 8, 8, 7, true>
                            (in, out); break;
                    case 9: HelmholtzPrismImpl<5, 5, 5, 9, 9, 8, true>
                            (in, out); break;
                    case 10: HelmholtzPrismImpl<5, 5, 5, 10, 10, 9, true>
                            (in, out); break;
                    default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                                in, out); break;
                    } break;
                case 6:
                    switch(nq0)
                    {
                    case 7: HelmholtzPrismImpl<6, 6, 6, 7, 7, 6, true>
                            (in, out); break;
                    case 8: HelmholtzPrismImpl<6, 6, 6, 8, 8, 7, true>
                            (in, out); break;
                    case 9: HelmholtzPrismImpl<6, 6, 6, 9, 9, 8, true>
                            (in, out); break;
                    case 10: HelmholtzPrismImpl<6, 6, 6, 10, 10, 9, true>
                            (in, out); break;
                    case 11: HelmholtzPrismImpl<6, 6, 6, 11, 11, 10, true>
                            (in, out); break;
                    case 12: HelmholtzPrismImpl<6, 6, 6, 12, 12, 11, true>
                            (in, out); break;
                    default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                                in, out); break;
                    } break;
                case 7:
                    switch(nq0)
                    {
                    case 8: HelmholtzPrismImpl<7, 7, 7, 8, 8, 7, true>
                            (in, out); break;
                    case 9: HelmholtzPrismImpl<7, 7, 7, 9, 9, 8, true>
                            (in, out); break;
                    case 10: HelmholtzPrismImpl<7, 7, 7, 10, 10, 9, true>
                            (in, out); break;
                    case 11: HelmholtzPrismImpl<7, 7, 7, 11, 11, 10, true>
                            (in, out); break;
                    case 12: HelmholtzPrismImpl<7, 7, 7, 12, 12, 11, true>
                            (in, out); break;
                    case 13: HelmholtzPrismImpl<7, 7, 7, 13, 13, 12, true>
                            (in, out); break;
                    case 14: HelmholtzPrismImpl<7, 7, 7, 14, 14, 13, true>
                            (in, out); break;
                    default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                                in, out); break;
                    } break;
                case 8:
                    switch(nq0)
                    {
                    case 9: HelmholtzPrismImpl<8, 8, 8, 9, 9, 8, true>
                            (in, out); break;
                    case 10: HelmholtzPrismImpl<8, 8, 8, 10, 10, 9, true>
                            (in, out); break;
                    case 11: HelmholtzPrismImpl<8, 8, 8, 11, 11, 10, true>
                            (in, out); break;
                    case 12: HelmholtzPrismImpl<8, 8, 8, 12, 12, 11, true>
                            (in, out); break;
                    case 13: HelmholtzPrismImpl<8, 8, 8, 13, 13, 12, true>
                            (in, out); break;
                    case 14: HelmholtzPrismImpl<8, 8, 8, 14, 14, 13, true>
                            (in, out); break;
                    case 15: HelmholtzPrismImpl<8, 8, 8, 15, 15, 14, true>
                            (in, out); break;
                    case 16: HelmholtzPrismImpl<8, 8, 8, 16, 16, 15, true>
                            (in, out); break;
                    default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
                                                in, out); break;
                    } break;
                default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2, true,
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
                    case 3: HelmholtzPrismImpl<2, 2, 2, 3, 3, 2, false>
                            (in, out); break;
                    case 4: HelmholtzPrismImpl<2, 2, 2, 4, 4, 3, false>
                            (in, out); break;
                    default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                                in, out); break;
                    } break;
                case 3:
                    switch(nq0)
                    {
                    case 4: HelmholtzPrismImpl<3, 3, 3, 4, 4, 3, false>
                            (in, out); break;
                    case 5: HelmholtzPrismImpl<3, 3, 3, 5, 5, 4, false>
                            (in, out); break;
                    case 6: HelmholtzPrismImpl<3, 3, 3, 6, 6, 5, false>
                            (in, out); break;
                    default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                                in, out); break;
                    } break;
                case 4:
                    switch(nq0)
                    {
                    case 5: HelmholtzPrismImpl<4, 4, 4, 5, 5, 4, false>
                            (in, out); break;
                    case 6: HelmholtzPrismImpl<4, 4, 4, 6, 6, 5, false>
                            (in, out); break;
                    case 7: HelmholtzPrismImpl<4, 4, 4, 7, 7, 6, false>
                            (in, out); break;
                    case 8: HelmholtzPrismImpl<4, 4, 4, 8, 8, 7, false>
                            (in, out); break;
                    default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                                in, out); break;
                    } break;
                case 5:
                    switch(nq0)
                    {
                    case 6: HelmholtzPrismImpl<5, 5, 5, 6, 6, 5, false>
                            (in, out); break;
                    case 7: HelmholtzPrismImpl<5, 5, 5, 7, 7, 6, false>
                            (in, out); break;
                    case 8: HelmholtzPrismImpl<5, 5, 5, 8, 8, 7, false>
                            (in, out); break;
                    case 9: HelmholtzPrismImpl<5, 5, 5, 9, 9, 8, false>
                            (in, out); break;
                    case 10: HelmholtzPrismImpl<5, 5, 5, 10, 10, 9, false>
                            (in, out); break;
                    default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                                in, out); break;
                    } break;
                case 6:
                    switch(nq0)
                    {
                    case 7: HelmholtzPrismImpl<6, 6, 6, 7, 7, 6, false>
                            (in, out); break;
                    case 8: HelmholtzPrismImpl<6, 6, 6, 8, 8, 7, false>
                            (in, out); break;
                    case 9: HelmholtzPrismImpl<6, 6, 6, 9, 9, 8, false>
                            (in, out); break;
                    case 10: HelmholtzPrismImpl<6, 6, 6, 10, 10, 9, false>
                            (in, out); break;
                    case 11: HelmholtzPrismImpl<6, 6, 6, 11, 11, 10, false>
                            (in, out); break;
                    case 12: HelmholtzPrismImpl<6, 6, 6, 12, 12, 11, false>
                            (in, out); break;
                    default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                                in, out); break;
                    } break;
                case 7:
                    switch(nq0)
                    {
                    case 8: HelmholtzPrismImpl<7, 7, 7, 8, 8, 7, false>
                            (in, out); break;
                    case 9: HelmholtzPrismImpl<7, 7, 7, 9, 9, 8, false>
                            (in, out); break;
                    case 10: HelmholtzPrismImpl<7, 7, 7, 10, 10, 9, false>
                            (in, out); break;
                    case 11: HelmholtzPrismImpl<7, 7, 7, 11, 11, 10, false>
                            (in, out); break;
                    case 12: HelmholtzPrismImpl<7, 7, 7, 12, 12, 11, false>
                            (in, out); break;
                    case 13: HelmholtzPrismImpl<7, 7, 7, 13, 13, 12, false>
                            (in, out); break;
                    case 14: HelmholtzPrismImpl<7, 7, 7, 14, 14, 13, false>
                            (in, out); break;
                    default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                                in, out); break;
                    } break;
                case 8:
                    switch(nq0)
                    {
                    case 9: HelmholtzPrismImpl<8, 8, 8, 9, 9, 8, false>
                            (in, out); break;
                    case 10: HelmholtzPrismImpl<8, 8, 8, 10, 10, 9, false>
                            (in, out); break;
                    case 11: HelmholtzPrismImpl<8, 8, 8, 11, 11, 10, false>
                            (in, out); break;
                    case 12: HelmholtzPrismImpl<8, 8, 8, 12, 12, 11, false>
                            (in, out); break;
                    case 13: HelmholtzPrismImpl<8, 8, 8, 13, 13, 12, false>
                            (in, out); break;
                    case 14: HelmholtzPrismImpl<8, 8, 8, 14, 14, 13, false>
                            (in, out); break;
                    case 15: HelmholtzPrismImpl<8, 8, 8, 15, 15, 14, false>
                            (in, out); break;
                    case 16: HelmholtzPrismImpl<8, 8, 8, 16, 16, 15, false>
                            (in, out); break;
                    default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                                in, out); break;
                    } break;
                default: HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                            in, out); break;
                }
            }
        }
        else
        {
            if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2,true,in,out);
            }
            else
            {
                HelmholtzPrismImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,in,out);
            }
        }
    }
    
    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void HelmholtzPrismImpl(const Array<OneD, const NekDouble> &input,
                            Array<OneD,       NekDouble> &out)
    {
        auto *inptr  = &input[0];
        auto *outptr = &out[0];

        constexpr auto ndf = 9;
        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        const auto nmBlocks = m_nmTot * vec_t::width;

        // Workspace for kernels
        vec_t wsp1[NQ1 * NQ2], wsp2[NQ2], wsp3[NM1];

        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot),  deriv0(nqTot),
            deriv1(nqTot), deriv2(nqTot);

        const vec_t* jac_ptr;
        const vec_t* df_ptr;

        vec_t df0,df1,df2,df3,df4,df5,df6,df7,df8;
        vec_t g0,g1,g2,g3,g4,g5; // metrics
        vec_t d00 = {1.0};
        vec_t d01 = {0.0};
        vec_t d11 = {1.0};
        vec_t d02 = {0.0};
        vec_t d12 = {0.0};
        vec_t d22 = {1.0}; // var diffusion terms
        vec_t td0,td1,td2,td3,td4,td5,td6,td7,td8; // temp terms for vardiff
        boost::ignore_unused(d00,d01,d11,d02,d12,d22,
                             td0,td1,td2,td3,td4,td5,td6,td7,td8);
                             
        if (this->m_isConstVarDiff)
        {
            d00 = this->m_constVarDiff[0];
            d01 = this->m_constVarDiff[1];
            d11 = this->m_constVarDiff[2];
            d02 = this->m_constVarDiff[3];
            d12 = this->m_constVarDiff[4];
            d22 = this->m_constVarDiff[5];
        }
        
        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            // Precompute Laplacian metrics
            df_ptr = &((*this->m_df)[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            if(!DEFORMED)
            {
                df0 = df_ptr[0];  df1 = df_ptr[1]; df2 = df_ptr[2];
                df3 = df_ptr[3];  df4 = df_ptr[4]; df5 = df_ptr[5];
                df6 = df_ptr[6];  df7 = df_ptr[7]; df8 = df_ptr[8];

                jac_ptr = &((*this->m_jac)[e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e*nqTot]);
            }

            // Step 1: BwdTrans
            BwdTransPrismKernel(NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT,
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                wsp1, wsp2, bwd);

            // Step 2: inner product for mass matrix operation
            IProductPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT,
                                true, false, DEFORMED>
                (bwd, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 wsp1, wsp2, wsp3, tmpOut, m_lambda);

            // Step 3: take derivatives in standard space
            PhysDerivTensor3DKernel(NQ0, NQ1, NQ2,
                 bwd, this->m_D[0], this->m_D[1], this->m_D[2], deriv0, deriv1,
                 deriv2);

            // Step 4: Apply Laplacian metrics & inner product

            // Step 4a: Construct Laplacian metrics
            for (size_t k = 0, cnt = 0; k < NQ2; ++k)
            {
                vec_t h1 = m_h1[k];
                for (size_t j = 0; j < NQ1; ++j)
                {
                    for (size_t i = 0; i < NQ0; ++i, cnt++)
                    {
                       vec_t  h0 = m_h0[i];

                        if(DEFORMED)
                        {
                            df0 = df_ptr[cnt * ndf];
                            df1 = df_ptr[cnt * ndf + 1];
                            df2 = df_ptr[cnt * ndf + 2];
                            df3 = df_ptr[cnt * ndf + 3];
                            df4 = df_ptr[cnt * ndf + 4];
                            df5 = df_ptr[cnt * ndf + 5];
                            df6 = df_ptr[cnt * ndf + 6];
                            df7 = df_ptr[cnt * ndf + 7];
                            df8 = df_ptr[cnt * ndf + 8];
                        }

                        vec_t tmp1 = h1 * (h0 * df2 + df0);
                        vec_t tmp2 = h1 * (h0 * df5 + df3);
                        vec_t tmp3 = h1 * (h0 * df8 + df6);
                        
                        if (!this->m_isConstVarDiff && !this->m_isVarDiff)
                        {
                            g0 = tmp1 * tmp1;
                            g0.fma(tmp2, tmp2);
                            g0.fma(tmp3, tmp3);

                            g3 = df1 * tmp1;
                            g3.fma(df4, tmp2);
                            g3.fma(df7, tmp3);
                            
                            g4 = df2 * tmp1;
                            g4.fma(df5, tmp2);
                            g4.fma(df8, tmp3);

                            g1 = df1 * df1;
                            g1.fma(df4, df4);
                            g1.fma(df7, df7);

                            g2 = df2 * df2;
                            g2.fma(df5, df5);
                            g2.fma(df8, df8);

                            g5 = df1 * df2;
                            g5.fma(df4, df5);
                            g5.fma(df7, df8);
                        }
                        else 
                        {
                            // vardiff
                            if (this->m_isVarDiff)
                            {
                                d00 = m_varD00[cnt];
                                d01 = m_varD01[cnt];
                                d11 = m_varD11[cnt];
                                d02 = m_varD02[cnt];
                                d12 = m_varD12[cnt];
                                d22 = m_varD22[cnt];
                            }
                            
                            td0 = tmp1 * d00;
                            td0.fma(tmp2,d01);
                            td0.fma(tmp3,d02);
                            
                            td1 = tmp1 * d01;
                            td1.fma(tmp2,d11);
                            td1.fma(tmp3,d12);
                            
                            td2 = tmp1 * d02;
                            td2.fma(tmp2,d12);
                            td2.fma(tmp3,d22);
                        
                            td3 = df1 * d00;
                            td3.fma(df4,d01);
                            td3.fma(df7,d02);
                            
                            td4 = df1 * d01;
                            td4.fma(df4,d11);
                            td4.fma(df7,d12);
                            
                            td5 = df1 * d02;
                            td5.fma(df4,d12);
                            td5.fma(df7,d22);
                                                
                            td6 = df2 * d00;
                            td6.fma(df5,d01);
                            td6.fma(df8,d02);
                            
                            td7 = df2 * d01;
                            td7.fma(df5,d11);
                            td7.fma(df8,d12);
                            
                            td8 = df2 * d02;
                            td8.fma(df5,d12);
                            td8.fma(df8,d22);
                            
                            g0 = td0 * tmp1;
                            g0.fma(td1,tmp2);
                            g0.fma(td2,tmp3);
                            
                            g3 = td0 * df1;
                            g3.fma(td1,df4);
                            g3.fma(td2,df7);
                            
                            g4 = td0 * df2;
                            g4.fma(td1,df5);
                            g4.fma(td2,df8);
                            
                            g1 = td3 * df1;
                            g1.fma(td4,df4);
                            g1.fma(td5,df7);
                            
                            g5 = td3 * df2;
                            g5.fma(td4,df5);
                            g5.fma(td5,df8);
                            
                            g2 = td6 * df2;
                            g2.fma(td7,df5);
                            g2.fma(td8,df8);
                        
                        }

                        vec_t d0 = deriv0[cnt];
                        vec_t d1 = deriv1[cnt];
                        vec_t d2 = deriv2[cnt];

                        tmp1 = g0 * d0;
                        tmp1.fma(g3, d1);
                        tmp1.fma(g4, d2);
                        deriv0[cnt] = tmp1;

                        tmp2 = g3 * d0;
                        tmp2.fma(g1, d1);
                        tmp2.fma(g5, d2);
                        deriv1[cnt] = tmp2;

                        tmp3 = g4 * d0;
                        tmp3.fma(g5, d1);
                        tmp3.fma(g2, d2);
                        deriv2[cnt] = tmp3;
                    }
                }
            }

            IProductPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2,
                                CORRECT, false, true, DEFORMED>
                (deriv0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                  this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                  wsp1, wsp2, wsp3, tmpOut);

            IProductPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2,
                                CORRECT, false, true, DEFORMED>
                (deriv1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 wsp1, wsp2, wsp3, tmpOut);

            IProductPrismKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2,
                                CORRECT, false, true, DEFORMED>
                (deriv2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 wsp1, wsp2, wsp3, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr  += nmBlocks;
            outptr += nmBlocks;
        }
    }

    void HelmholtzPrismImpl(
        const int nm0, const int nm1, const int nm2,
        const int nq0, const int nq1, const int nq2,
        const bool CORRECT,
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out)
    {
        auto *inptr  = &input[0];
        auto *outptr = &out[0];

        constexpr auto ndf = 9;
        const auto nqTot = nq0 * nq1 * nq2; 
        const auto nmBlocks = m_nmTot * vec_t::width;

        // Workspace for kernels
        std::vector<vec_t, allocator<vec_t>> wsp1(nq1 * nq2), wsp2(nq2),
            wsp3(nm1), tmpIn(m_nmTot), tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot),  deriv0(nqTot),
            deriv1(nqTot), deriv2(nqTot);

        const vec_t* jac_ptr;
        const vec_t* df_ptr;

        vec_t df0,df1,df2,df3,df4,df5,df6,df7,df8;
        vec_t g0,g1,g2,g3,g4,g5; // metrics
        vec_t d00 = {1.0};
        vec_t d01 = {0.0};
        vec_t d11 = {1.0};
        vec_t d02 = {0.0};
        vec_t d12 = {0.0};
        vec_t d22 = {1.0}; // var diffusion terms
        vec_t td0,td1,td2,td3,td4,td5,td6,td7,td8; // temp terms for vardiff
        boost::ignore_unused(d00,d01,d11,d02,d12,d22,
                             td0,td1,td2,td3,td4,td5,td6,td7,td8);
                             
        if (this->m_isConstVarDiff)
        {
            d00 = this->m_constVarDiff[0];
            d01 = this->m_constVarDiff[1];
            d11 = this->m_constVarDiff[2];
            d02 = this->m_constVarDiff[3];
            d12 = this->m_constVarDiff[4];
            d22 = this->m_constVarDiff[5];
        }
        
        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            // Precompute Laplacian metrics
            df_ptr = &((*this->m_df)[e*dfSize]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            if(!DEFORMED)
            {
                df0 = df_ptr[0];  df1 = df_ptr[1]; df2 = df_ptr[2];
                df3 = df_ptr[3];  df4 = df_ptr[4]; df5 = df_ptr[5];
                df6 = df_ptr[6];  df7 = df_ptr[7]; df8 = df_ptr[8];

                jac_ptr = &((*this->m_jac)[e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e*nqTot]);
            }

            // Step 1: BwdTrans
            BwdTransPrismKernel(nm0,nm1,nm2,nq0,nq1,nq2, CORRECT,
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                &wsp1[0], &wsp2[0], bwd);

            // Step 2: inner product for mass matrix operation
            IProductPrismKernel(nm0,nm1,nm2,nq0,nq1,nq2, CORRECT,
                true, false, DEFORMED,
                bwd, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                &wsp1[0], &wsp2[0], &wsp3[0], tmpOut, m_lambda);
        
            // Step 3: take derivatives in standard space
            PhysDerivTensor3DKernel(nq0,nq1,nq2,
                bwd, this->m_D[0], this->m_D[1], this->m_D[2], deriv0, deriv1,
                                    deriv2);

            // Step 4: Apply Laplacian metrics & inner product
            
            // Step 4a: Construct Laplacian metrics
            for (size_t k = 0, cnt = 0; k < nq2; ++k)
            {
                vec_t h1 = m_h1[k];
                for (size_t j = 0; j < nq1; ++j)
                {
                    for (size_t i = 0; i < nq0; ++i, cnt++)
                    {
                       vec_t  h0 = m_h0[i];

                        if(DEFORMED)
                        {
                            df0 = df_ptr[cnt * ndf];
                            df1 = df_ptr[cnt * ndf + 1];
                            df2 = df_ptr[cnt * ndf + 2];
                            df3 = df_ptr[cnt * ndf + 3];
                            df4 = df_ptr[cnt * ndf + 4];
                            df5 = df_ptr[cnt * ndf + 5];
                            df6 = df_ptr[cnt * ndf + 6];
                            df7 = df_ptr[cnt * ndf + 7];
                            df8 = df_ptr[cnt * ndf + 8];
                        }

                        vec_t tmp1 = h1 * (h0 * df2 + df0);
                        vec_t tmp2 = h1 * (h0 * df5 + df3);
                        vec_t tmp3 = h1 * (h0 * df8 + df6);
                        
                        if (!this->m_isConstVarDiff && !this->m_isVarDiff)
                        {
                            g0 = tmp1 * tmp1;
                            g0.fma(tmp2, tmp2);
                            g0.fma(tmp3, tmp3);

                            g3 = df1 * tmp1;
                            g3.fma(df4, tmp2);
                            g3.fma(df7, tmp3);
                            
                            g4 = df2 * tmp1;
                            g4.fma(df5, tmp2);
                            g4.fma(df8, tmp3);

                            g1 = df1 * df1;
                            g1.fma(df4, df4);
                            g1.fma(df7, df7);

                            g2 = df2 * df2;
                            g2.fma(df5, df5);
                            g2.fma(df8, df8);

                            g5 = df1 * df2;
                            g5.fma(df4, df5);
                            g5.fma(df7, df8);
                        }
                        else 
                        {
                            // vardiff
                            if (this->m_isVarDiff)
                            {
                                d00 = m_varD00[cnt];
                                d01 = m_varD01[cnt];
                                d11 = m_varD11[cnt];
                                d02 = m_varD02[cnt];
                                d12 = m_varD12[cnt];
                                d22 = m_varD22[cnt];
                            }
                            
                            td0 = tmp1 * d00;
                            td0.fma(tmp2,d01);
                            td0.fma(tmp3,d02);
                            
                            td1 = tmp1 * d01;
                            td1.fma(tmp2,d11);
                            td1.fma(tmp3,d12);
                            
                            td2 = tmp1 * d02;
                            td2.fma(tmp2,d12);
                            td2.fma(tmp3,d22);
                        
                            td3 = df1 * d00;
                            td3.fma(df4,d01);
                            td3.fma(df7,d02);
                            
                            td4 = df1 * d01;
                            td4.fma(df4,d11);
                            td4.fma(df7,d12);
                            
                            td5 = df1 * d02;
                            td5.fma(df4,d12);
                            td5.fma(df7,d22);
                                                
                            td6 = df2 * d00;
                            td6.fma(df5,d01);
                            td6.fma(df8,d02);
                            
                            td7 = df2 * d01;
                            td7.fma(df5,d11);
                            td7.fma(df8,d12);
                            
                            td8 = df2 * d02;
                            td8.fma(df5,d12);
                            td8.fma(df8,d22);
                            
                            g0 = td0 * tmp1;
                            g0.fma(td1,tmp2);
                            g0.fma(td2,tmp3);
                            
                            g3 = td0 * df1;
                            g3.fma(td1,df4);
                            g3.fma(td2,df7);
                            
                            g4 = td0 * df2;
                            g4.fma(td1,df5);
                            g4.fma(td2,df8);
                            
                            g1 = td3 * df1;
                            g1.fma(td4,df4);
                            g1.fma(td5,df7);
                            
                            g5 = td3 * df2;
                            g5.fma(td4,df5);
                            g5.fma(td5,df8);
                            
                            g2 = td6 * df2;
                            g2.fma(td7,df5);
                            g2.fma(td8,df8);
                        
                        }

                        vec_t d0 = deriv0[cnt];
                        vec_t d1 = deriv1[cnt];
                        vec_t d2 = deriv2[cnt];

                        tmp1 = g0 * d0;
                        tmp1.fma(g3, d1);
                        tmp1.fma(g4, d2);
                        deriv0[cnt] = tmp1;

                        tmp2 = g3 * d0;
                        tmp2.fma(g1, d1);
                        tmp2.fma(g5, d2);
                        deriv1[cnt] = tmp2;

                        tmp3 = g4 * d0;
                        tmp3.fma(g5, d1);
                        tmp3.fma(g2, d2);
                        deriv2[cnt] = tmp3;
                    }
                }
            }

            IProductPrismKernel(nm0,nm1,nm2,nq0,nq1,nq2,
                                CORRECT, false, true, DEFORMED,
                  deriv0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                  this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                  &wsp1[0], &wsp2[0], &wsp3[0], tmpOut);

            IProductPrismKernel(nm0,nm1,nm2,nq0,nq1,nq2,
                                CORRECT, false, true, DEFORMED,
                 deriv1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 &wsp1[0], &wsp2[0], &wsp3[0], tmpOut);

            IProductPrismKernel(nm0,nm1,nm2,nq0,nq1,nq2,
                                CORRECT, false, true, DEFORMED,
                 deriv2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 &wsp1[0], &wsp2[0], &wsp3[0], tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr  += nmBlocks;
            outptr += nmBlocks;
        }
    }


public:

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }

private:
    int m_nmTot;
    std::vector<vec_t, allocator<vec_t>> m_h0;
    std::vector<vec_t, allocator<vec_t>> m_h1;
    std::vector<vec_t, allocator<vec_t>> m_h2;
};

    
template<bool DEFORMED = false>
struct HelmholtzPyr : public Helmholtz, public Helper<3, DEFORMED>
{
    HelmholtzPyr(std::vector<LibUtilities::BasisSharedPtr> basis,
                      int nElmt)
    : Helmholtz(basis, nElmt),
      Helper<3, DEFORMED>(basis, nElmt),
      m_nmTot(LibUtilities::StdPyrData::getNumberOfCoefficients(
                  this->m_nm[0], this->m_nm[1], this->m_nm[2])),
      m_h0(basis[0]->GetNumPoints()),
      m_h1(basis[1]->GetNumPoints()),
      m_h2(basis[2]->GetNumPoints())
    {
        const int nq0 = basis[0]->GetNumPoints();
        const int nq1 = basis[1]->GetNumPoints();
        const int nq2 = basis[2]->GetNumPoints();

        const Array<OneD, const NekDouble> &z0 = basis[0]->GetZ();
        const Array<OneD, const NekDouble> &z1 = basis[1]->GetZ();
        const Array<OneD, const NekDouble> &z2 = basis[2]->GetZ();

        for (int i = 0; i < nq0; ++i)
        {
            m_h0[i] = 0.5 * (1 + z0[i]);
        }

        for (int j = 0; j< nq1; ++j)
        {
            m_h1[j] = 0.5 * (1 + z1[j]);
        }

        for (int k = 0; k< nq2; ++k)
        {
            m_h2[k] = 2.0 / (1 - z2[k]);
        }
    }

    static std::shared_ptr<Operator> Create(
        std::vector<LibUtilities::BasisSharedPtr> basis,
        int nElmt)
    {
        return std::make_shared<HelmholtzPyr<DEFORMED>>(basis, nElmt);
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out) final
    {
        const int nm0 = m_basis[0]->GetNumModes();  
        const int nm1 = m_basis[1]->GetNumModes();  
        const int nm2 = m_basis[2]->GetNumModes();  
        const int nq0 = m_basis[0]->GetNumPoints();  
        const int nq1 = m_basis[1]->GetNumPoints();  
        const int nq2 = m_basis[2]->GetNumPoints();  

        if((nm0 == nm1)&&(nm0 == nm2)&&(nq0 == nq1)&&(nq0 == nq2+1))
        {

            if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                switch(nm0)
                {
                case 2:
                    switch(nq0)
                    {
                    case 3: HelmholtzPyrImpl<2, 2, 2, 3, 3, 2, true>
                            (in, out); break;
                    case 4: HelmholtzPyrImpl<2, 2, 2, 4, 4, 3, true>
                            (in, out); break;
                    default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,true,
                                             in, out); break;
                    } break;
                case 3:
                    switch(nq0)
                    {
                    case 4: HelmholtzPyrImpl<3, 3, 3, 4, 4, 3, true>
                            (in, out); break;
                    case 5: HelmholtzPyrImpl<3, 3, 3, 5, 5, 4, true>
                            (in, out); break;
                    case 6: HelmholtzPyrImpl<3, 3, 3, 6, 6, 5, true>
                            (in, out); break;
                    default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,true,
                                              in, out); break;
                    } break;
                case 4:
                    switch(nq0)
                    {
                    case 5: HelmholtzPyrImpl<4, 4, 4, 5, 5, 4, true>
                            (in, out); break;
                    case 6: HelmholtzPyrImpl<4, 4, 4, 6, 6, 5, true>
                            (in, out); break;
                    case 7: HelmholtzPyrImpl<4, 4, 4, 7, 7, 6, true>
                            (in, out); break;
                    case 8: HelmholtzPyrImpl<4, 4, 4, 8, 8, 7, true>
                            (in, out); break;
                    default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,true,
                                              in, out); break;
                    } break;
                case 5:
                    switch(nq0)
                    {
                    case 6: HelmholtzPyrImpl<5, 5, 5, 6, 6, 5, true>
                            (in, out); break;
                    case 7: HelmholtzPyrImpl<5, 5, 5, 7, 7, 6, true>
                            (in, out); break;
                    case 8: HelmholtzPyrImpl<5, 5, 5, 8, 8, 7, true>
                            (in, out); break;
                    case 9: HelmholtzPyrImpl<5, 5, 5, 9, 9, 8, true>
                            (in, out); break;
                    case 10: HelmholtzPyrImpl<5, 5, 5, 10, 10, 9, true>
                            (in, out); break;
                    default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,true,
                                              in, out); break;
                    } break;
                case 6:
                    switch(nq0)
                    {
                    case 7: HelmholtzPyrImpl<6, 6, 6, 7, 7, 6, true>
                            (in, out); break;
                    case 8: HelmholtzPyrImpl<6, 6, 6, 8, 8, 7, true>
                            (in, out); break;
                    case 9: HelmholtzPyrImpl<6, 6, 6, 9, 9, 8, true>
                            (in, out); break;
                    case 10: HelmholtzPyrImpl<6, 6, 6, 10, 10, 9, true>
                            (in, out); break;
                    case 11: HelmholtzPyrImpl<6, 6, 6, 11, 11, 10, true>
                            (in, out); break;
                    case 12: HelmholtzPyrImpl<6, 6, 6, 12, 12, 11, true>
                            (in, out); break;
                    default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,true,
                                              in, out); break;
                    } break;
                case 7:
                    switch(nq0)
                    {
                    case 8: HelmholtzPyrImpl<7, 7, 7, 8, 8, 7, true>
                            (in, out); break;
                    case 9: HelmholtzPyrImpl<7, 7, 7, 9, 9, 8, true>
                            (in, out); break;
                    case 10: HelmholtzPyrImpl<7, 7, 7, 10, 10, 9, true>
                            (in, out); break;
                    case 11: HelmholtzPyrImpl<7, 7, 7, 11, 11, 10, true>
                            (in, out); break;
                    case 12: HelmholtzPyrImpl<7, 7, 7, 12, 12, 11, true>
                            (in, out); break;
                    case 13: HelmholtzPyrImpl<7, 7, 7, 13, 13, 12, true>
                            (in, out); break;
                    case 14: HelmholtzPyrImpl<7, 7, 7, 14, 14, 13, true>
                            (in, out); break;
                    default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,true,
                                              in, out); break;
                    } break;
                case 8:
                    switch(nq0)
                    {
                    case 9: HelmholtzPyrImpl<8, 8, 8, 9, 9, 8, true>
                            (in, out); break;
                    case 10: HelmholtzPyrImpl<8, 8, 8, 10, 10, 9, true>
                            (in, out); break;
                    case 11: HelmholtzPyrImpl<8, 8, 8, 11, 11, 10, true>
                            (in, out); break;
                    case 12: HelmholtzPyrImpl<8, 8, 8, 12, 12, 11, true>
                            (in, out); break;
                    case 13: HelmholtzPyrImpl<8, 8, 8, 13, 13, 12, true>
                            (in, out); break;
                    case 14: HelmholtzPyrImpl<8, 8, 8, 14, 14, 13, true>
                            (in, out); break;
                    case 15: HelmholtzPyrImpl<8, 8, 8, 15, 15, 14, true>
                            (in, out); break;
                    case 16: HelmholtzPyrImpl<8, 8, 8, 16, 16, 15, true>
                            (in, out); break;
                    default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,true,
                                              in, out); break;
                    } break;
                default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,true,
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
                    case 3: HelmholtzPyrImpl<2, 2, 2, 3, 3, 2, false>
                            (in, out); break;
                    case 4: HelmholtzPyrImpl<2, 2, 2, 4, 4, 3, false>
                            (in, out); break;
                    default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                              in, out); break;
                    }
                    break;
                case 3:
                    switch(nq0)
                    {
                    case 4: HelmholtzPyrImpl<3, 3, 3, 4, 4, 3, false>
                            (in, out); break;
                    case 5: HelmholtzPyrImpl<3, 3, 3, 5, 5, 4, false>
                            (in, out); break;
                    case 6: HelmholtzPyrImpl<3, 3, 3, 6, 6, 5, false>
                            (in, out); break;
                    default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                              in, out); break;
                    } break;
                case 4:
                    switch(nq0)
                    {
                    case 5: HelmholtzPyrImpl<4, 4, 4, 5, 5, 4, false>
                            (in, out); break;
                    case 6: HelmholtzPyrImpl<4, 4, 4, 6, 6, 5, false>
                            (in, out); break;
                    case 7: HelmholtzPyrImpl<4, 4, 4, 7, 7, 6, false>
                            (in, out); break;
                    case 8: HelmholtzPyrImpl<4, 4, 4, 8, 8, 7, false>
                            (in, out); break;
                    default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                              in, out); break;
                    } break;
                case 5:
                    switch(nq0)
                    {
                    case 6: HelmholtzPyrImpl<5, 5, 5, 6, 6, 5, false>
                            (in, out); break;
                    case 7: HelmholtzPyrImpl<5, 5, 5, 7, 7, 6, false>
                            (in, out); break;
                    case 8: HelmholtzPyrImpl<5, 5, 5, 8, 8, 7, false>
                            (in, out); break;
                    case 9: HelmholtzPyrImpl<5, 5, 5, 9, 9, 8, false>
                            (in, out); break;
                    case 10: HelmholtzPyrImpl<5, 5, 5, 10, 10, 9, false>
                            (in, out); break;
                    default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                              in, out); break;
                    } break;
                case 6:
                    switch(nq0)
                    {
                    case 7: HelmholtzPyrImpl<6, 6, 6, 7, 7, 6, false>
                            (in, out); break;
                    case 8: HelmholtzPyrImpl<6, 6, 6, 8, 8, 7, false>
                            (in, out); break;
                    case 9: HelmholtzPyrImpl<6, 6, 6, 9, 9, 8, false>
                            (in, out); break;
                    case 10: HelmholtzPyrImpl<6, 6, 6, 10, 10, 9, false>
                            (in, out); break;
                    case 11: HelmholtzPyrImpl<6, 6, 6, 11, 11, 10, false>
                            (in, out); break;
                    case 12: HelmholtzPyrImpl<6, 6, 6, 12, 12, 11, false>
                            (in, out); break;
                    default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                              in, out); break;
                    } break;
                case 7:
                    switch(nq0)
                    {
                    case 8: HelmholtzPyrImpl<7, 7, 7, 8, 8, 7, false>
                            (in, out); break;
                    case 9: HelmholtzPyrImpl<7, 7, 7, 9, 9, 8, false>
                            (in, out); break;
                    case 10: HelmholtzPyrImpl<7, 7, 7, 10, 10, 9, false>
                            (in, out); break;
                    case 11: HelmholtzPyrImpl<7, 7, 7, 11, 11, 10, false>
                            (in, out); break;
                    case 12: HelmholtzPyrImpl<7, 7, 7, 12, 12, 11, false>
                            (in, out); break;
                    case 13: HelmholtzPyrImpl<7, 7, 7, 13, 13, 12, false>
                            (in, out); break;
                    case 14: HelmholtzPyrImpl<7, 7, 7, 14, 14, 13, false>
                            (in, out); break;
                    default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                              in, out); break;
                    } break;
                case 8:
                    switch(nq0)
                    {
                    case 9: HelmholtzPyrImpl<8, 8, 8, 9, 9, 8, false>
                            (in, out); break;
                    case 10: HelmholtzPyrImpl<8, 8, 8, 10, 10, 9, false>
                            (in, out); break;
                    case 11: HelmholtzPyrImpl<8, 8, 8, 11, 11, 10, false>
                            (in, out); break;
                    case 12: HelmholtzPyrImpl<8, 8, 8, 12, 12, 11, false>
                            (in, out); break;
                    case 13: HelmholtzPyrImpl<8, 8, 8, 13, 13, 12, false>
                            (in, out); break;
                    case 14: HelmholtzPyrImpl<8, 8, 8, 14, 14, 13, false>
                            (in, out); break;
                    case 15: HelmholtzPyrImpl<8, 8, 8, 15, 15, 14, false>
                            (in, out); break;
                    case 16: HelmholtzPyrImpl<8, 8, 8, 16, 16, 15, false>
                            (in, out); break;
                    default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                              in, out); break;
                    } break;
                default: HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false,
                                          in, out); break;
                    
                }
            }
        }
        else
        {
            if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,true,in, out);
            }
            else
            {
                HelmholtzPyrImpl(nm0,nm1,nm2,nq0,nq1,nq2,false, in, out);
            }
        }
    }
                
    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void HelmholtzPyrImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out)
    {
        auto *inptr  = &input[0];
        auto *outptr = &out[0];

        constexpr auto ndf = 9;
        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        const auto nmBlocks = m_nmTot * vec_t::width;

        // Workspace for kernels
        vec_t wsp1[NQ1 * NQ2], wsp2[NQ2]; 

        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot),  deriv0(nqTot),
            deriv1(nqTot), deriv2(nqTot);

        const vec_t* jac_ptr;
        const vec_t* df_ptr;

        vec_t df0,df1,df2,df3,df4,df5,df6,df7,df8;
        vec_t g0,g1,g2,g3,g4,g5; // metrics
        vec_t d00 = {1.0};
        vec_t d01 = {0.0};
        vec_t d11 = {1.0};
        vec_t d02 = {0.0};
        vec_t d12 = {0.0};
        vec_t d22 = {1.0}; // var diffusion terms
        vec_t td0,td1,td2,td3,td4,td5,td6,td7,td8; // temp terms for vardiff
        boost::ignore_unused(d00,d01,d11,d02,d12,d22,
                             td0,td1,td2,td3,td4,td5,td6,td7,td8);
                             
        if (this->m_isConstVarDiff)
        {
            d00 = this->m_constVarDiff[0];
            d01 = this->m_constVarDiff[1];
            d11 = this->m_constVarDiff[2];
            d02 = this->m_constVarDiff[3];
            d12 = this->m_constVarDiff[4];
            d22 = this->m_constVarDiff[5];
        }
        
        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            // Precompute Laplacian metrics
            df_ptr = &((*this->m_df)[dfSize*e]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            if(!DEFORMED)
            {
                df0 = df_ptr[0];  df1 = df_ptr[1]; df2 = df_ptr[2];
                df3 = df_ptr[3];  df4 = df_ptr[4]; df5 = df_ptr[5];
                df6 = df_ptr[6];  df7 = df_ptr[7]; df8 = df_ptr[8];

                jac_ptr = &((*this->m_jac)[e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e*nqTot]);
            }

            // Step 1: BwdTrans
            BwdTransPyrKernel(NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT,
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                wsp1, wsp2, bwd);

            // Step 2: inner product for mass matrix operation
            IProductPyrKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT,
                              true, false, DEFORMED>
                (bwd, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 wsp1, wsp2, tmpOut, m_lambda);

            // Step 3: take derivatives in standard space
            PhysDerivTensor3DKernel(NQ0, NQ1, NQ2,
                 bwd, this->m_D[0], this->m_D[1], this->m_D[2], deriv0, deriv1,
                 deriv2);

            // Step 4: Apply Laplacian metrics & inner product

            // Step 4a: Construct Laplacian metrics
            for (size_t k = 0, cnt = 0; k < NQ2; ++k)
            {
                vec_t h2 = m_h2[k];
                for (size_t j = 0; j < NQ1; ++j)
                {
                    vec_t h1 = m_h1[j];
                    vec_t h1h2 = h1 * h2;
                    for (size_t i = 0; i < NQ0; ++i, cnt++)
                    {
                        vec_t  h0 = m_h0[i];
                        vec_t  h0h2 = h0 * h2;

                        if(DEFORMED)
                        {
                            df0 = df_ptr[cnt * ndf];
                            df1 = df_ptr[cnt * ndf + 1];
                            df2 = df_ptr[cnt * ndf + 2];
                            df3 = df_ptr[cnt * ndf + 3];
                            df4 = df_ptr[cnt * ndf + 4];
                            df5 = df_ptr[cnt * ndf + 5];
                            df6 = df_ptr[cnt * ndf + 6];
                            df7 = df_ptr[cnt * ndf + 7];
                            df8 = df_ptr[cnt * ndf + 8];
                        }

                        vec_t tmp0 = h2*df0;
                        tmp0.fma(h0h2,df2);
                        vec_t tmp1 = h2*df3;
                        tmp1.fma(h0h2,df5);
                        vec_t tmp2 = h2*df6;
                        tmp2.fma(h0h2,df8);

                        vec_t tmp3 = h2*df1;
                        tmp3.fma(h1h2,df2);
                        vec_t tmp4 = h2*df4;
                        tmp4.fma(h1h2,df5);
                        vec_t tmp5 = h2*df7;
                        tmp5.fma(h1h2,df8);
                        
                        if (!this->m_isConstVarDiff && !this->m_isVarDiff)
                        {
                            g0 = tmp0 * tmp0;
                            g0.fma(tmp1, tmp1);
                            g0.fma(tmp2, tmp2);

                            g1 = tmp3 * tmp3;
                            g1.fma(tmp4, tmp4);
                            g1.fma(tmp5, tmp5);
                            
                            g2 = df2 * df2;
                            g2.fma(df5, df5);
                            g2.fma(df8, df8);
                            
                            g3 = tmp0 * tmp3;
                            g3.fma(tmp1, tmp4);
                            g3.fma(tmp2, tmp5);
                            
                            g4 = df2 * tmp0;
                            g4.fma(df5, tmp1);
                            g4.fma(df8, tmp2);

                            g5 = df2 * tmp3;
                            g5.fma(df5, tmp4);
                            g5.fma(df8, tmp5);
                        }
                        else
                        {
                            if (this->m_isVarDiff)
                            {
                                d00 = m_varD00[cnt];
                                d01 = m_varD01[cnt];
                                d11 = m_varD11[cnt];
                                d02 = m_varD02[cnt];
                                d12 = m_varD12[cnt];
                                d22 = m_varD22[cnt];
                            }
                            
                            td0 = tmp0 * d00;
                            td0.fma(tmp1,d01);
                            td0.fma(tmp2,d02);
                            
                            td1 = tmp0 * d01;
                            td1.fma(tmp1,d11);
                            td1.fma(tmp2,d12);
                            
                            td2 = tmp0 * d02;
                            td2.fma(tmp1,d12);
                            td2.fma(tmp2,d22);
                        
                            td3 = tmp3 * d00;
                            td3.fma(tmp4,d01);
                            td3.fma(tmp5,d02);
                            
                            td4 = tmp3 * d01;
                            td4.fma(tmp4,d11);
                            td4.fma(tmp5,d12);
                            
                            td5 = tmp3 * d02;
                            td5.fma(tmp4,d12);
                            td5.fma(tmp5,d22);
                                                
                            td6 = df2 * d00;
                            td6.fma(df5,d01);
                            td6.fma(df8,d02);
                            
                            td7 = df2 * d01;
                            td7.fma(df5,d11);
                            td7.fma(df8,d12);
                            
                            td8 = df2 * d02;
                            td8.fma(df5,d12);
                            td8.fma(df8,d22);
                            
                            g0 = td0 * tmp0;
                            g0.fma(td1,tmp1);
                            g0.fma(td2,tmp2);
                            
                            g3 = td0 * tmp3;
                            g3.fma(td1,tmp4);
                            g3.fma(td2,tmp5);
                            
                            g4 = td0 * df2;
                            g4.fma(td1,df5);
                            g4.fma(td2,df8);
                            
                            g1 = td3 * tmp3;
                            g1.fma(td4,tmp4);
                            g1.fma(td5,tmp5);
                            
                            g5 = td3 * df2;
                            g5.fma(td4,df5);
                            g5.fma(td5,df8);
                            
                            g2 = td6 * df2;
                            g2.fma(td7,df5);
                            g2.fma(td8,df8);
                        
                        }

                        vec_t d0 = deriv0[cnt];
                        vec_t d1 = deriv1[cnt];
                        vec_t d2 = deriv2[cnt];

                        tmp1 = g0 * d0;
                        tmp1.fma(g3, d1);
                        tmp1.fma(g4, d2);
                        deriv0[cnt] = tmp1;

                        tmp2 = g3 * d0;
                        tmp2.fma(g1, d1);
                        tmp2.fma(g5, d2);
                        deriv1[cnt] = tmp2;

                        tmp3 = g4 * d0;
                        tmp3.fma(g5, d1);
                        tmp3.fma(g2, d2);
                        deriv2[cnt] = tmp3;
                    }
                }
            }

            IProductPyrKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2,
                              CORRECT, false, true, DEFORMED>
                (deriv0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                  this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                  wsp1, wsp2, tmpOut);

            IProductPyrKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2,
                              CORRECT, false, true, DEFORMED>
                (deriv1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 wsp1, wsp2, tmpOut);

            IProductPyrKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2,
                              CORRECT, false, true, DEFORMED>
                (deriv2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 wsp1, wsp2, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr  += nmBlocks;
            outptr += nmBlocks;
        }
    }

    void HelmholtzPyrImpl(
        const int nm0, const int nm1, const int nm2,
        const int nq0, const int nq1, const int nq2,
        const bool CORRECT,
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out)
    {
        auto *inptr  = &input[0];
        auto *outptr = &out[0];

        constexpr auto ndf = 9;
        const auto nqTot = nq0 * nq1 * nq2;
        const auto nmBlocks = m_nmTot * vec_t::width;

        std::vector<vec_t, allocator<vec_t>>  wsp1(nq1 * nq2), wsp2(nq2),
            tmpIn(m_nmTot), tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot),  deriv0(nqTot),
            deriv1(nqTot), deriv2(nqTot);

        const vec_t* jac_ptr;
        const vec_t* df_ptr;

        vec_t df0,df1,df2,df3,df4,df5,df6,df7,df8;
        vec_t g0,g1,g2,g3,g4,g5; // metrics
        vec_t d00 = {1.0};
        vec_t d01 = {0.0};
        vec_t d11 = {1.0};
        vec_t d02 = {0.0};
        vec_t d12 = {0.0};
        vec_t d22 = {1.0}; // var diffusion terms
        vec_t td0,td1,td2,td3,td4,td5,td6,td7,td8; // temp terms for vardiff
        boost::ignore_unused(d00,d01,d11,d02,d12,d22,
                             td0,td1,td2,td3,td4,td5,td6,td7,td8);
                             
        if (this->m_isConstVarDiff)
        {
            d00 = this->m_constVarDiff[0];
            d01 = this->m_constVarDiff[1];
            d11 = this->m_constVarDiff[2];
            d02 = this->m_constVarDiff[3];
            d12 = this->m_constVarDiff[4];
            d22 = this->m_constVarDiff[5];
        }
        
        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            // Precompute Laplacian metrics
            df_ptr = &((*this->m_df)[dfSize*e]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            if(!DEFORMED)
            {
                df0 = df_ptr[0];  df1 = df_ptr[1]; df2 = df_ptr[2];
                df3 = df_ptr[3];  df4 = df_ptr[4]; df5 = df_ptr[5];
                df6 = df_ptr[6];  df7 = df_ptr[7]; df8 = df_ptr[8];

                jac_ptr = &((*this->m_jac)[e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e*nqTot]);
            }

            // Step 1: BwdTrans
            BwdTransPyrKernel(nm0,nm1,nm2,nq0,nq1,nq2, CORRECT,
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                              &wsp1[0], &wsp2[0], bwd);

            // Step 2: inner product for mass matrix operation
            IProductPyrKernel(nm0,nm1,nm2,nq0,nq1,nq2, CORRECT,
                              true, false, DEFORMED,
                 bwd, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 &wsp1[0], &wsp2[0], tmpOut, m_lambda);

            // Step 3: take derivatives in standard space
            PhysDerivTensor3DKernel(nq0, nq1, nq2,
                 bwd, this->m_D[0], this->m_D[1], this->m_D[2], deriv0, deriv1,
                 deriv2);

            // Step 4: Apply Laplacian metrics & inner product

            // Step 4a: Construct Laplacian metrics
            for (size_t k = 0, cnt = 0; k < nq2; ++k)
            {
                vec_t h2 = m_h2[k];
                for (size_t j = 0; j < nq1; ++j)
                {
                    vec_t h1 = m_h1[j];
                    vec_t h1h2 = h1 * h2;
                    for (size_t i = 0; i < nq0; ++i, cnt++)
                    {
                        vec_t  h0 = m_h0[i];
                        vec_t  h0h2 = h0 * h2;

                        if(DEFORMED)
                        {
                            df0 = df_ptr[cnt * ndf];
                            df1 = df_ptr[cnt * ndf + 1];
                            df2 = df_ptr[cnt * ndf + 2];
                            df3 = df_ptr[cnt * ndf + 3];
                            df4 = df_ptr[cnt * ndf + 4];
                            df5 = df_ptr[cnt * ndf + 5];
                            df6 = df_ptr[cnt * ndf + 6];
                            df7 = df_ptr[cnt * ndf + 7];
                            df8 = df_ptr[cnt * ndf + 8];
                        }

                        vec_t tmp0 = h2*df0;
                        tmp0.fma(h0h2,df2);
                        vec_t tmp1 = h2*df3;
                        tmp1.fma(h0h2,df5);
                        vec_t tmp2 = h2*df6;
                        tmp2.fma(h0h2,df8);

                        vec_t tmp3 = h2*df1;
                        tmp3.fma(h1h2,df2);
                        vec_t tmp4 = h2*df4;
                        tmp4.fma(h1h2,df5);
                        vec_t tmp5 = h2*df7;
                        tmp5.fma(h1h2,df8);
                        
                        if (!this->m_isConstVarDiff && !this->m_isVarDiff)
                        {
                            g0 = tmp0 * tmp0;
                            g0.fma(tmp1, tmp1);
                            g0.fma(tmp2, tmp2);

                            g1 = tmp3 * tmp3;
                            g1.fma(tmp4, tmp4);
                            g1.fma(tmp5, tmp5);
                            
                            g2 = df2 * df2;
                            g2.fma(df5, df5);
                            g2.fma(df8, df8);
                            
                            g3 = tmp0 * tmp3;
                            g3.fma(tmp1, tmp4);
                            g3.fma(tmp2, tmp5);
                            
                            g4 = df2 * tmp0;
                            g4.fma(df5, tmp1);
                            g4.fma(df8, tmp2);

                            g5 = df2 * tmp3;
                            g5.fma(df5, tmp4);
                            g5.fma(df8, tmp5);
                        }
                        else
                        {
                            if (this->m_isVarDiff)
                            {
                                d00 = m_varD00[cnt];
                                d01 = m_varD01[cnt];
                                d11 = m_varD11[cnt];
                                d02 = m_varD02[cnt];
                                d12 = m_varD12[cnt];
                                d22 = m_varD22[cnt];
                            }
                            
                            td0 = tmp0 * d00;
                            td0.fma(tmp1,d01);
                            td0.fma(tmp2,d02);
                            
                            td1 = tmp0 * d01;
                            td1.fma(tmp1,d11);
                            td1.fma(tmp2,d12);
                            
                            td2 = tmp0 * d02;
                            td2.fma(tmp1,d12);
                            td2.fma(tmp2,d22);
                        
                            td3 = tmp3 * d00;
                            td3.fma(tmp4,d01);
                            td3.fma(tmp5,d02);
                            
                            td4 = tmp3 * d01;
                            td4.fma(tmp4,d11);
                            td4.fma(tmp5,d12);
                            
                            td5 = tmp3 * d02;
                            td5.fma(tmp4,d12);
                            td5.fma(tmp5,d22);
                                                
                            td6 = df2 * d00;
                            td6.fma(df5,d01);
                            td6.fma(df8,d02);
                            
                            td7 = df2 * d01;
                            td7.fma(df5,d11);
                            td7.fma(df8,d12);
                            
                            td8 = df2 * d02;
                            td8.fma(df5,d12);
                            td8.fma(df8,d22);
                            
                            g0 = td0 * tmp0;
                            g0.fma(td1,tmp1);
                            g0.fma(td2,tmp2);
                            
                            g3 = td0 * tmp3;
                            g3.fma(td1,tmp4);
                            g3.fma(td2,tmp5);
                            
                            g4 = td0 * df2;
                            g4.fma(td1,df5);
                            g4.fma(td2,df8);
                            
                            g1 = td3 * tmp3;
                            g1.fma(td4,tmp4);
                            g1.fma(td5,tmp5);
                            
                            g5 = td3 * df2;
                            g5.fma(td4,df5);
                            g5.fma(td5,df8);
                            
                            g2 = td6 * df2;
                            g2.fma(td7,df5);
                            g2.fma(td8,df8);
                        
                        }

                        vec_t d0 = deriv0[cnt];
                        vec_t d1 = deriv1[cnt];
                        vec_t d2 = deriv2[cnt];

                        tmp1 = g0 * d0;
                        tmp1.fma(g3, d1);
                        tmp1.fma(g4, d2);
                        deriv0[cnt] = tmp1;

                        tmp2 = g3 * d0;
                        tmp2.fma(g1, d1);
                        tmp2.fma(g5, d2);
                        deriv1[cnt] = tmp2;

                        tmp3 = g4 * d0;
                        tmp3.fma(g5, d1);
                        tmp3.fma(g2, d2);
                        deriv2[cnt] = tmp3;
                    }
                }
            }

            IProductPyrKernel(nm0,nm1,nm2,nq0,nq1,nq2,
                              CORRECT, false, true, DEFORMED,
                 deriv0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                  this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                  &wsp1[0], &wsp2[0], tmpOut);

            IProductPyrKernel(nm0,nm1,nm2,nq0,nq1,nq2,
                              CORRECT, false, true, DEFORMED,
                 deriv1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 &wsp1[0], &wsp2[0], tmpOut);

            IProductPyrKernel(nm0,nm1,nm2,nq0,nq1,nq2,
                              CORRECT, false, true, DEFORMED,
                 deriv2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 &wsp1[0], &wsp2[0], tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr  += nmBlocks;
            outptr += nmBlocks;
        }
    }
    
public:

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }

private:
    int m_nmTot;
    std::vector<vec_t, allocator<vec_t>> m_h0;
    std::vector<vec_t, allocator<vec_t>> m_h1;
    std::vector<vec_t, allocator<vec_t>> m_h2;
};


template<bool DEFORMED = false>
struct HelmholtzTet : public Helmholtz, public Helper<3, DEFORMED>
{
    HelmholtzTet(std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
        : Helmholtz(basis, nElmt),
          Helper<3, DEFORMED>(basis, nElmt),
          m_nmTot(LibUtilities::StdTetData::getNumberOfCoefficients
                  (this->m_nm[0], this->m_nm[1], this->m_nm[2])),
          m_h0(basis[0]->GetNumPoints()),
          m_h1(basis[1]->GetNumPoints()),
          m_h2(basis[1]->GetNumPoints()),
          m_h3(basis[2]->GetNumPoints())
    {
        const int nq0 = basis[0]->GetNumPoints();
        const int nq1 = basis[1]->GetNumPoints();
        const int nq2 = basis[2]->GetNumPoints();

        const Array<OneD, const NekDouble> &z0 = basis[0]->GetZ();
        const Array<OneD, const NekDouble> &z1 = basis[1]->GetZ();
        const Array<OneD, const NekDouble> &z2 = basis[2]->GetZ();

        for (int i = 0; i < nq0; ++i)
        {
            m_h0[i] = 0.5 * (1 + z0[i]);
        }

        for (int j = 0; j < nq1; ++j)
        {
            m_h1[j] = 0.5 * (1 + z1[j]);
            m_h2[j] = 2.0 / (1 - z1[j]);
        }

        for (int k = 0; k < nq2; ++k)
        {
            m_h3[k] = 2.0 / (1 - z2[k]);
        }
    }
    
    static std::shared_ptr<Operator> Create
    (std::vector<LibUtilities::BasisSharedPtr> basis, int nElmt)
    {
        return std::make_shared<HelmholtzTet<DEFORMED>>(basis, nElmt);
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                  Array<OneD,       NekDouble> &out) final
    {
        const int nm0 = m_basis[0]->GetNumModes();  
        const int nm1 = m_basis[1]->GetNumModes();  
        const int nm2 = m_basis[2]->GetNumModes();  
        const int nq0 = m_basis[0]->GetNumPoints();  
        const int nq1 = m_basis[1]->GetNumPoints();  
        const int nq2 = m_basis[2]->GetNumPoints();  

        if((nm0 == nm1)&&(nm0 == nm2)&&(nq0 == nq1+1)&&(nq0 == nq2+1))
        {
            if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                switch(nm0)
                {
                case 2:
                    switch(nq0)
                    {
                    case 3: HelmholtzTetImpl<2, 2, 2, 3, 2, 2, true>
                            (in, out); break;
                    case 4: HelmholtzTetImpl<2, 2, 2, 4, 3, 3, true>
                            (in, out); break;
                    default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                              true, in, out); break;
                    } break;
                case 3:
                    switch(nq0)
                    {
                    case 4: HelmholtzTetImpl<3, 3, 3, 4, 3, 3, true>
                            (in, out); break;
                    case 5: HelmholtzTetImpl<3, 3, 3, 5, 4, 4, true>
                            (in, out); break;
                    case 6: HelmholtzTetImpl<3, 3, 3, 6, 5, 5, true>
                            (in, out); break;
                    default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                              true, in, out); break;
                    } break;
                case 4:
                    switch(nq0)
                    {
                    case 5: HelmholtzTetImpl<4, 4, 4, 5, 4, 4, true>
                            (in, out); break;
                    case 6: HelmholtzTetImpl<4, 4, 4, 6, 5, 5, true>
                            (in, out); break;
                    case 7: HelmholtzTetImpl<4, 4, 4, 7, 6, 6, true>
                            (in, out); break;
                    case 8: HelmholtzTetImpl<4, 4, 4, 8, 7, 7, true>
                            (in, out); break;
                    default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                              true, in, out); break;
                    } break;
                case 5:
                    switch(nq0)
                    {
                    case 6: HelmholtzTetImpl<5, 5, 5, 6, 5, 5, true>
                            (in, out); break;
                    case 7: HelmholtzTetImpl<5, 5, 5, 7, 6, 6, true>
                            (in, out); break;
                    case 8: HelmholtzTetImpl<5, 5, 5, 8, 7, 7, true>
                            (in, out); break;
                    case 9: HelmholtzTetImpl<5, 5, 5, 9, 8, 8, true>
                            (in, out); break;
                    case 10: HelmholtzTetImpl<5, 5, 5, 10, 9, 9, true>
                            (in, out); break;
                    default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                              true, in, out); break;
                    } break;
                case 6:
                    switch(nq0)
                    {
                    case 7: HelmholtzTetImpl<6, 6, 6, 7, 6, 6, true>
                            (in, out); break;
                    case 8: HelmholtzTetImpl<6, 6, 6, 8, 7, 7, true>
                            (in, out); break;
                    case 9: HelmholtzTetImpl<6, 6, 6, 9, 8, 8, true>
                            (in, out); break;
                    case 10: HelmholtzTetImpl<6, 6, 6, 10, 9, 9, true>
                            (in, out); break;
                    case 11: HelmholtzTetImpl<6, 6, 6, 11, 10, 10, true>
                            (in, out); break;
                    case 12: HelmholtzTetImpl<6, 6, 6, 12, 11, 11, true>
                            (in, out); break;
                    default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                              true, in, out); break;
                    } break;
                case 7:
                    switch(nq0)
                    {
                    case 8: HelmholtzTetImpl<7, 7, 7, 8, 7, 7, true>
                            (in, out); break;
                    case 9: HelmholtzTetImpl<7, 7, 7, 9, 8, 8, true>
                            (in, out); break;
                    case 10: HelmholtzTetImpl<7, 7, 7, 10, 9, 9, true>
                            (in, out); break;
                    case 11: HelmholtzTetImpl<7, 7, 7, 11, 10, 10, true>
                            (in, out); break;
                    case 12: HelmholtzTetImpl<7, 7, 7, 12, 11, 11, true>
                            (in, out); break;
                    case 13: HelmholtzTetImpl<7, 7, 7, 13, 12, 12, true>
                            (in, out); break;
                    case 14: HelmholtzTetImpl<7, 7, 7, 14, 13, 13, true>
                            (in, out); break;
                    default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                              true, in, out); break;
                    } break;
                case 8:
                    switch(nq0)
                    {
                    case 9: HelmholtzTetImpl<8, 8, 8, 9, 8, 8, true>
                            (in, out); break;
                    case 10: HelmholtzTetImpl<8, 8, 8, 10, 9, 9, true>
                            (in, out); break;
                    case 11: HelmholtzTetImpl<8, 8, 8, 11, 10, 10, true>
                            (in, out); break;
                    case 12: HelmholtzTetImpl<8, 8, 8, 12, 11, 11, true>
                            (in, out); break;
                    case 13: HelmholtzTetImpl<8, 8, 8, 13, 12, 12, true>
                            (in, out); break;
                    case 14: HelmholtzTetImpl<8, 8, 8, 14, 13, 13, true>
                            (in, out); break;
                    case 15: HelmholtzTetImpl<8, 8, 8, 15, 14, 14, true>
                            (in, out); break;
                    case 16: HelmholtzTetImpl<8, 8, 8, 16, 15, 15, true>
                            (in, out); break;
                    default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                              true, in, out); break;
                    } break;
                default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                          true, in, out); break;
                }
            }
            else
            {
                switch(nm0)
                {
                case 2:
                    switch(nq0)
                    {
                    case 3: HelmholtzTetImpl<2, 2, 2, 3, 2, 2, false>
                            (in, out); break;
                    case 4: HelmholtzTetImpl<2, 2, 2, 4, 3, 3, false>
                            (in, out); break;
                    default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                              false, in, out); break;
                    } break;
                case 3:
                    switch(nq0)
                    {
                    case 4: HelmholtzTetImpl<3, 3, 3, 4, 3, 3, false>
                            (in, out); break;
                    case 5: HelmholtzTetImpl<3, 3, 3, 5, 4, 4, false>
                            (in, out); break;
                    case 6: HelmholtzTetImpl<3, 3, 3, 6, 5, 5, false>
                            (in, out); break;
                    default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                              false, in, out); break;
                    } break;
                case 4:
                    switch(nq0)
                    {
                    case 5: HelmholtzTetImpl<4, 4, 4, 5, 4, 4, false>
                            (in, out); break;
                    case 6: HelmholtzTetImpl<4, 4, 4, 6, 5, 5, false>
                            (in, out); break;
                    case 7: HelmholtzTetImpl<4, 4, 4, 7, 6, 6, false>
                            (in, out); break;
                    case 8: HelmholtzTetImpl<4, 4, 4, 8, 7, 7, false>
                            (in, out); break;
                    default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                              false, in, out); break;
                    } break;
                case 5:
                    switch(nq0)
                    {
                    case 6: HelmholtzTetImpl<5, 5, 5, 6, 5, 5, false>
                            (in, out); break;
                    case 7: HelmholtzTetImpl<5, 5, 5, 7, 6, 6, false>
                            (in, out); break;
                    case 8: HelmholtzTetImpl<5, 5, 5, 8, 7, 7, false>
                            (in, out); break;
                    case 9: HelmholtzTetImpl<5, 5, 5, 9, 8, 8, false>
                            (in, out); break;
                    case 10: HelmholtzTetImpl<5, 5, 5, 10, 9, 9, false>
                            (in, out); break;
                    default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                              false, in, out); break;
                    } break;
                case 6:
                    switch(nq0)
                    {
                    case 7: HelmholtzTetImpl<6, 6, 6, 7, 6, 6, false>
                            (in, out); break;
                    case 8: HelmholtzTetImpl<6, 6, 6, 8, 7, 7, false>
                            (in, out); break;
                    case 9: HelmholtzTetImpl<6, 6, 6, 9, 8, 8, false>
                            (in, out); break;
                    case 10: HelmholtzTetImpl<6, 6, 6, 10, 9, 9, false>
                            (in, out); break;
                    case 11: HelmholtzTetImpl<6, 6, 6, 11, 10, 10, false>
                            (in, out); break;
                    case 12: HelmholtzTetImpl<6, 6, 6, 12, 11, 11, false>
                            (in, out); break;
                    default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                              false, in, out); break;
                    } break;
                case 7:
                    switch(nq0)
                    {
                    case 8: HelmholtzTetImpl<7, 7, 7, 8, 7, 7, false>
                            (in, out); break;
                    case 9: HelmholtzTetImpl<7, 7, 7, 9, 8, 8, false>
                            (in, out); break;
                    case 10: HelmholtzTetImpl<7, 7, 7, 10, 9, 9, false>
                            (in, out); break;
                    case 11: HelmholtzTetImpl<7, 7, 7, 11, 10, 10, false>
                            (in, out); break;
                    case 12: HelmholtzTetImpl<7, 7, 7, 12, 11, 11, false>
                            (in, out); break;
                    case 13: HelmholtzTetImpl<7, 7, 7, 13, 12, 12, false>
                            (in, out); break;
                    case 14: HelmholtzTetImpl<7, 7, 7, 14, 13, 13, false>
                            (in, out); break;
                    default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                              false, in, out); break;
                    } break;
                case 8:
                    switch(nq0)
                    {
                    case 9: HelmholtzTetImpl<8, 8, 8, 9, 8, 8, false>
                            (in, out); break;
                    case 10: HelmholtzTetImpl<8, 8, 8, 10, 9, 9, false>
                            (in, out); break;
                    case 11: HelmholtzTetImpl<8, 8, 8, 11, 10, 10, false>
                            (in, out); break;
                    case 12: HelmholtzTetImpl<8, 8, 8, 12, 11, 11, false>
                            (in, out); break;
                    case 13: HelmholtzTetImpl<8, 8, 8, 13, 12, 12, false>
                            (in, out); break;
                    case 14: HelmholtzTetImpl<8, 8, 8, 14, 13, 13, false>
                            (in, out); break;
                    case 15: HelmholtzTetImpl<8, 8, 8, 15, 14, 14, false>
                            (in, out); break;
                    case 16: HelmholtzTetImpl<8, 8, 8, 16, 15, 15, false>
                            (in, out); break;
                    default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                              false, in, out); break;
                    } break;
                default: HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2,
                                          false, in, out); break;
                }
            }
        }
        else
        {
            if (m_basis[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2, true, in, out);
            }
            else
            {
                HelmholtzTetImpl(nm0, nm1, nm2, nq0, nq1, nq2, false, in, out);
            }
        }
    }

    template<int NM0, int NM1, int NM2, int NQ0, int NQ1, int NQ2, bool CORRECT>
    void HelmholtzTetImpl(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out)
    {
        auto *inptr  = &input[0];
        auto *outptr = &out[0];

        constexpr auto ndf = 9;
        constexpr auto nqTot = NQ0 * NQ1 * NQ2;
        const auto nmBlocks = m_nmTot * vec_t::width;

        // Workspace for kernels
        vec_t wsp1[NQ1 * NQ2 + NQ2], wsp2[NM0];

        std::vector<vec_t, allocator<vec_t>> tmpIn(m_nmTot), tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot),  deriv0(nqTot),
            deriv1(nqTot), deriv2(nqTot);

        const vec_t* jac_ptr;
        const vec_t* df_ptr;
        
        vec_t df0,df1,df2,df3,df4,df5,df6,df7,df8;
        vec_t g0,g1,g2,g3,g4,g5;
        vec_t d00 = {1.0};
        vec_t d01 = {0.0};
        vec_t d11 = {1.0};
        vec_t d02 = {0.0};
        vec_t d12 = {0.0};
        vec_t d22 = {1.0}; // var diffusion terms
        vec_t td0,td1,td2,td3,td4,td5,td6,td7,td8; // temp terms for vardiff
        boost::ignore_unused(d00,d01,d11,d02,d12,d22,
                             td0,td1,td2,td3,td4,td5,td6,td7,td8);
                             
        if (this->m_isConstVarDiff)
        {
            d00 = this->m_constVarDiff[0];
            d01 = this->m_constVarDiff[1];
            d11 = this->m_constVarDiff[2];
            d02 = this->m_constVarDiff[3];
            d12 = this->m_constVarDiff[4];
            d22 = this->m_constVarDiff[5];
        }
        
        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            // Precompute Laplacian metrics
            df_ptr = &((*this->m_df)[dfSize*e]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            if(!DEFORMED)
            {
                df0 = df_ptr[0];  df1 = df_ptr[1]; df2 = df_ptr[2];
                df3 = df_ptr[3];  df4 = df_ptr[4]; df5 = df_ptr[5];
                df6 = df_ptr[6];  df7 = df_ptr[7]; df8 = df_ptr[8];

                jac_ptr = &((*this->m_jac)[e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e*nqTot]);
            }

            // Step 1: BwdTrans
            BwdTransTetKernel(NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT,
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                wsp1, wsp2, bwd);

            // Step 2: inner product for mass matrix operation
            IProductTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT,
                              true, false, DEFORMED>
                (bwd, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 wsp1, tmpOut, m_lambda);

            // Step 3: take derivatives in standard space
            PhysDerivTensor3DKernel(NQ0, NQ1, NQ2,
                bwd, this->m_D[0], this->m_D[1], this->m_D[2],
                deriv0, deriv1, deriv2);

            // Step 4: Apply Laplacian metrics & inner product

            // Step 4a: Construct Laplacian metrics
            for (size_t k = 0, cnt = 0; k < NQ2; ++k)
            {
                vec_t h3 = m_h3[k];
                for (size_t j = 0; j < NQ1; ++j)
                {
                    vec_t h1 = m_h1[j];
                    vec_t h2 = m_h2[j];
                    vec_t h2h3 = h2 * h3;
                    vec_t h1h3 = h1 * h3;

                    for (int i = 0; i < NQ0; ++i, ++cnt)
                    {
                        if(DEFORMED)
                        {
                            df0 = df_ptr[cnt * ndf];
                            df1 = df_ptr[cnt * ndf + 1];
                            df2 = df_ptr[cnt * ndf + 2];
                            df3 = df_ptr[cnt * ndf + 3];
                            df4 = df_ptr[cnt * ndf + 4];
                            df5 = df_ptr[cnt * ndf + 5];
                            df6 = df_ptr[cnt * ndf + 6];
                            df7 = df_ptr[cnt * ndf + 7];
                            df8 = df_ptr[cnt * ndf + 8];
                        }

                        vec_t h0h2h3 = m_h0[i] * h2h3;

                        vec_t tmp1 = h0h2h3 * (df1 + df2);
                        tmp1.fma(df0, h2h3);
                        vec_t tmp2 = h0h2h3 * (df4 + df5);
                        tmp2.fma(df3, h2h3);
                        vec_t tmp3 = h0h2h3 * (df7 + df8);
                        tmp3.fma(df6, h2h3);

                        vec_t tmp4 = df1 * h3;
                        tmp4.fma(df2, h1h3);
                        vec_t tmp5 = df4 * h3;
                        tmp5.fma(df5, h1h3);
                        vec_t tmp6 = df7 * h3;
                        tmp6.fma(df8, h1h3);
                        
                        if (!this->m_isConstVarDiff && !this->m_isVarDiff)
                        {
                            g0 = tmp1 * tmp1;
                            g0.fma(tmp2, tmp2);
                            g0.fma(tmp3, tmp3);

                            g4 = df2 * tmp1;
                            g4.fma(df5, tmp2);
                            g4.fma(df8, tmp3);

                            g3 = tmp1 * tmp4;
                            g3.fma(tmp2, tmp5);
                            g3.fma(tmp3, tmp6);

                            g1 = tmp4 * tmp4;
                            g1.fma(tmp5, tmp5);
                            g1.fma(tmp6, tmp6);

                            g5 = df2 * tmp4;
                            g5.fma(df5, tmp5);
                            g5.fma(df8, tmp6);

                            g2 = df2 * df2;
                            g2.fma(df5, df5);
                            g2.fma(df8, df8);
                        }
                        else
                        {
                            if (this->m_isVarDiff)
                            {
                                d00 = m_varD00[cnt];
                                d01 = m_varD01[cnt];
                                d02 = m_varD02[cnt];
                                d12 = m_varD12[cnt];
                                d22 = m_varD22[cnt];
                            }
                                
                            td0 = tmp1 * d00;
                            td0.fma(tmp2,d01);
                            td0.fma(tmp3,d02);
                            
                            td1 = tmp1 * d01;
                            td1.fma(tmp2,d11);
                            td1.fma(tmp3,d12);
                            
                            td2 = tmp1 * d02;
                            td2.fma(tmp2,d12);
                            td2.fma(tmp3,d22);
                        
                            td3 = tmp4 * d00;
                            td3.fma(tmp5,d01);
                            td3.fma(tmp6,d02);
                            
                            td4 = tmp4 * d01;
                            td4.fma(tmp5,d11);
                            td4.fma(tmp6,d12);
                            
                            td5 = tmp4 * d02;
                            td5.fma(tmp5,d12);
                            td5.fma(tmp6,d22);
                                                
                            td6 = df2 * d00;
                            td6.fma(df5,d01);
                            td6.fma(df8,d02);
                            
                            td7 = df2 * d01;
                            td7.fma(df5,d11);
                            td7.fma(df8,d12);
                            
                            td8 = df2 * d02;
                            td8.fma(df5,d12);
                            td8.fma(df8,d22);
                            
                            g0 = td0 * tmp1;
                            g0.fma(td1,tmp2);
                            g0.fma(td2,tmp3);
                            
                            g3 = td0 * tmp4;
                            g3.fma(td1,tmp5);
                            g3.fma(td2,tmp6);
                            
                            g4 = td0 * df2;
                            g4.fma(td1,df5);
                            g4.fma(td2,df8);
                            
                            g1 = td3 * tmp4;
                            g1.fma(td4,tmp5);
                            g1.fma(td5,tmp6);
                            
                            g5 = td3 * df2;
                            g5.fma(td4,df5);
                            g5.fma(td5,df8);
                            
                            g2 = td6 * df2;
                            g2.fma(td7,df5);
                            g2.fma(td8,df8);
                        
                        }

                        vec_t d0 = deriv0[cnt];
                        vec_t d1 = deriv1[cnt];
                        vec_t d2 = deriv2[cnt];

                        tmp1 = g0 * d0;
                        tmp1.fma(g3, d1);
                        tmp1.fma(g4, d2);
                        deriv0[cnt] = tmp1; 

                        tmp2 = g3 * d0;
                        tmp2.fma(g1, d1);
                        tmp2.fma(g5, d2);
                        deriv1[cnt] = tmp2;
                        
                        tmp3 = g4 * d0;
                        tmp3.fma(g5, d1);
                        tmp3.fma(g2, d2);
                        deriv2[cnt] = tmp3;
                    }
                }
            }

            IProductTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT,
                              false, true, DEFORMED>
                (deriv0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 wsp1, tmpOut);

            IProductTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT,
                              false, true, DEFORMED>
                (deriv1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 wsp1, tmpOut);

            IProductTetKernel<NM0, NM1, NM2, NQ0, NQ1, NQ2, CORRECT,
                              false, true, DEFORMED>
                (deriv2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 wsp1, tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr  += nmBlocks;
            outptr += nmBlocks;
        }
    }

    void HelmholtzTetImpl(
        const int nm0, const int nm1, const int nm2,
        const int nq0, const int nq1, const int nq2,
        const bool CORRECT,
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &out)
    {
        auto *inptr  = &input[0];
        auto *outptr = &out[0];

        constexpr auto ndf = 9;
        const auto nqTot = nq0 * nq1 * nq2;
        const auto nmBlocks = m_nmTot * vec_t::width;

        std::vector<vec_t, allocator<vec_t>> wsp1(nqTot), wsp2(nm0),
            tmpIn(m_nmTot), tmpOut(m_nmTot);
        std::vector<vec_t, allocator<vec_t>> bwd(nqTot),  deriv0(nqTot),
            deriv1(nqTot), deriv2(nqTot);

        const vec_t* jac_ptr;
        const vec_t* df_ptr;
        
        vec_t df0,df1,df2,df3,df4,df5,df6,df7,df8;
        vec_t g0,g1,g2,g3,g4,g5;
        vec_t d00 = {1.0};
        vec_t d01 = {0.0};
        vec_t d11 = {1.0};
        vec_t d02 = {0.0};
        vec_t d12 = {0.0};
        vec_t d22 = {1.0}; // var diffusion terms
        vec_t td0,td1,td2,td3,td4,td5,td6,td7,td8; // temp terms for vardiff
        boost::ignore_unused(d00,d01,d11,d02,d12,d22,
                             td0,td1,td2,td3,td4,td5,td6,td7,td8);
                             
        if (this->m_isConstVarDiff)
        {
            d00 = this->m_constVarDiff[0];
            d01 = this->m_constVarDiff[1];
            d11 = this->m_constVarDiff[2];
            d02 = this->m_constVarDiff[3];
            d12 = this->m_constVarDiff[4];
            d22 = this->m_constVarDiff[5];
        }
        
        // Get size of derivative factor block
        auto dfSize = ndf;
        if (DEFORMED)
        {
            dfSize *= nqTot;
        }

        for (int e = 0; e < this->m_nBlocks; e++)
        {
            // Precompute Laplacian metrics
            df_ptr = &((*this->m_df)[dfSize*e]);

            // Load and transpose data
            load_interleave(inptr, m_nmTot, tmpIn);

            if(!DEFORMED)
            {
                df0 = df_ptr[0];  df1 = df_ptr[1]; df2 = df_ptr[2];
                df3 = df_ptr[3];  df4 = df_ptr[4]; df5 = df_ptr[5];
                df6 = df_ptr[6];  df7 = df_ptr[7]; df8 = df_ptr[8];

                jac_ptr = &((*this->m_jac)[e]);
            }
            else
            {
                jac_ptr = &((*this->m_jac)[e*nqTot]);
            }

            // Step 1: BwdTrans
            BwdTransTetKernel(nm0, nm1, nm2, nq0, nq1, nq2, CORRECT,
                tmpIn, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                &wsp1[0], &wsp2[0], bwd);

            // Step 2: inner product for mass matrix operation
            IProductTetKernel(nm0, nm1, nm2, nq0, nq1, nq2, CORRECT,
                              true, false, DEFORMED,
                 bwd, this->m_bdata[0], this->m_bdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 &wsp1[0], tmpOut, m_lambda);

            // Step 3: take derivatives in standard space
            PhysDerivTensor3DKernel(nq0, nq1, nq2,
                bwd, this->m_D[0], this->m_D[1], this->m_D[2],
                deriv0, deriv1, deriv2);

            // Step 4: Apply Laplacian metrics & inner product

            // Step 4a: Construct Laplacian metrics
            for (size_t k = 0, cnt = 0; k < nq2; ++k)
            {
                vec_t h3 = m_h3[k];
                for (size_t j = 0; j < nq1; ++j)
                {
                    vec_t h1 = m_h1[j];
                    vec_t h2 = m_h2[j];
                    vec_t h2h3 = h2 * h3;
                    vec_t h1h3 = h1 * h3;

                    for (int i = 0; i < nq0; ++i, ++cnt)
                    {
                        if(DEFORMED)
                        {
                            df0 = df_ptr[cnt * ndf];
                            df1 = df_ptr[cnt * ndf + 1];
                            df2 = df_ptr[cnt * ndf + 2];
                            df3 = df_ptr[cnt * ndf + 3];
                            df4 = df_ptr[cnt * ndf + 4];
                            df5 = df_ptr[cnt * ndf + 5];
                            df6 = df_ptr[cnt * ndf + 6];
                            df7 = df_ptr[cnt * ndf + 7];
                            df8 = df_ptr[cnt * ndf + 8];
                        }

                        vec_t h0h2h3 = m_h0[i] * h2h3;

                        vec_t tmp1 = h0h2h3 * (df1 + df2);
                        tmp1.fma(df0, h2h3);
                        vec_t tmp2 = h0h2h3 * (df4 + df5);
                        tmp2.fma(df3, h2h3);
                        vec_t tmp3 = h0h2h3 * (df7 + df8);
                        tmp3.fma(df6, h2h3);

                        vec_t tmp4 = df1 * h3;
                        tmp4.fma(df2, h1h3);
                        vec_t tmp5 = df4 * h3;
                        tmp5.fma(df5, h1h3);
                        vec_t tmp6 = df7 * h3;
                        tmp6.fma(df8, h1h3);
                        
                        if (!this->m_isConstVarDiff && !this->m_isVarDiff)
                        {
                            g0 = tmp1 * tmp1;
                            g0.fma(tmp2, tmp2);
                            g0.fma(tmp3, tmp3);

                            g4 = df2 * tmp1;
                            g4.fma(df5, tmp2);
                            g4.fma(df8, tmp3);

                            g3 = tmp1 * tmp4;
                            g3.fma(tmp2, tmp5);
                            g3.fma(tmp3, tmp6);

                            g1 = tmp4 * tmp4;
                            g1.fma(tmp5, tmp5);
                            g1.fma(tmp6, tmp6);

                            g5 = df2 * tmp4;
                            g5.fma(df5, tmp5);
                            g5.fma(df8, tmp6);

                            g2 = df2 * df2;
                            g2.fma(df5, df5);
                            g2.fma(df8, df8);
                        }
                        else
                        {
                            if (this->m_isVarDiff)
                            {
                                d00 = m_varD00[cnt];
                                d01 = m_varD01[cnt];
                                d02 = m_varD02[cnt];
                                d12 = m_varD12[cnt];
                                d22 = m_varD22[cnt];
                            }
                                
                            td0 = tmp1 * d00;
                            td0.fma(tmp2,d01);
                            td0.fma(tmp3,d02);
                            
                            td1 = tmp1 * d01;
                            td1.fma(tmp2,d11);
                            td1.fma(tmp3,d12);
                            
                            td2 = tmp1 * d02;
                            td2.fma(tmp2,d12);
                            td2.fma(tmp3,d22);
                        
                            td3 = tmp4 * d00;
                            td3.fma(tmp5,d01);
                            td3.fma(tmp6,d02);
                            
                            td4 = tmp4 * d01;
                            td4.fma(tmp5,d11);
                            td4.fma(tmp6,d12);
                            
                            td5 = tmp4 * d02;
                            td5.fma(tmp5,d12);
                            td5.fma(tmp6,d22);
                                                
                            td6 = df2 * d00;
                            td6.fma(df5,d01);
                            td6.fma(df8,d02);
                            
                            td7 = df2 * d01;
                            td7.fma(df5,d11);
                            td7.fma(df8,d12);
                            
                            td8 = df2 * d02;
                            td8.fma(df5,d12);
                            td8.fma(df8,d22);
                            
                            g0 = td0 * tmp1;
                            g0.fma(td1,tmp2);
                            g0.fma(td2,tmp3);
                            
                            g3 = td0 * tmp4;
                            g3.fma(td1,tmp5);
                            g3.fma(td2,tmp6);
                            
                            g4 = td0 * df2;
                            g4.fma(td1,df5);
                            g4.fma(td2,df8);
                            
                            g1 = td3 * tmp4;
                            g1.fma(td4,tmp5);
                            g1.fma(td5,tmp6);
                            
                            g5 = td3 * df2;
                            g5.fma(td4,df5);
                            g5.fma(td5,df8);
                            
                            g2 = td6 * df2;
                            g2.fma(td7,df5);
                            g2.fma(td8,df8);
                        
                        }

                        vec_t d0 = deriv0[cnt];
                        vec_t d1 = deriv1[cnt];
                        vec_t d2 = deriv2[cnt];

                        tmp1 = g0 * d0;
                        tmp1.fma(g3, d1);
                        tmp1.fma(g4, d2);
                        deriv0[cnt] = tmp1; 

                        tmp2 = g3 * d0;
                        tmp2.fma(g1, d1);
                        tmp2.fma(g5, d2);
                        deriv1[cnt] = tmp2;
                        
                        tmp3 = g4 * d0;
                        tmp3.fma(g5, d1);
                        tmp3.fma(g2, d2);
                        deriv2[cnt] = tmp3;
                    }
                }
            }

            IProductTetKernel(nm0, nm1, nm2, nq0, nq1, nq2, CORRECT,
                              false, true, DEFORMED,
                deriv0, this->m_dbdata[0], this->m_bdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 &wsp1[0], tmpOut);

            IProductTetKernel(nm0, nm1, nm2, nq0, nq1, nq2, CORRECT,
                              false, true, DEFORMED,
                 deriv1, this->m_bdata[0], this->m_dbdata[1], this->m_bdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 &wsp1[0], tmpOut);

            IProductTetKernel(nm0, nm1, nm2, nq0, nq1, nq2, CORRECT,
                              false, true, DEFORMED,
                 deriv2, this->m_bdata[0], this->m_bdata[1], this->m_dbdata[2],
                 this->m_w[0], this->m_w[1], this->m_w[2], jac_ptr,
                 &wsp1[0], tmpOut);

            // de-interleave and store data
            deinterleave_store(tmpOut, m_nmTot, outptr);

            inptr  += nmBlocks;
            outptr += nmBlocks;
        }
    }
    
public:

    virtual NekDouble Ndof() override
    {
        return m_nmTot * this->m_nElmt;
    }
      
private:
    int m_nmTot;
    std::vector<vec_t, allocator<vec_t>> m_h0, m_h1, m_h2, m_h3;
};

} // namespace MatrixFree
} // namespace Nektar

#endif
