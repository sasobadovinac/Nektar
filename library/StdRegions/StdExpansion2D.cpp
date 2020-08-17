///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion2D.cpp
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
// Description: Daughter of StdExpansion. This class contains routine
// which are common to 2D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion2D.h>

#ifdef max
#undef max
#endif

namespace Nektar
{
    namespace StdRegions
    {

        StdExpansion2D::StdExpansion2D()
        {
        }

        StdExpansion2D::StdExpansion2D(int numcoeffs,
                                       const LibUtilities::BasisKey &Ba,
                                       const LibUtilities::BasisKey &Bb):
            StdExpansion(numcoeffs,2, Ba, Bb)
        {
        }

        StdExpansion2D::StdExpansion2D(const StdExpansion2D &T):
            StdExpansion(T)
        {
        }

        StdExpansion2D::~StdExpansion2D()
        {
        }

        //----------------------------
        // Differentiation Methods
        //----------------------------
        
        // find derivative of u (inarray) at all quad points
        void StdExpansion2D::PhysTensorDeriv(const Array<OneD, const NekDouble>& inarray,
                                             Array<OneD, NekDouble> &outarray_d0,
                                             Array<OneD, NekDouble> &outarray_d1)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();

            if (outarray_d0.size() > 0) // calculate du/dx_0
            {
                DNekMatSharedPtr D0 = m_base[0]->GetD();
                if(inarray.data() == outarray_d0.data())
                {
                    Array<OneD, NekDouble> wsp(nquad0 * nquad1);
                    Vmath::Vcopy(nquad0 * nquad1,inarray.get(),1,wsp.get(),1);
                    Blas::Dgemm('N', 'N', nquad0, nquad1, nquad0, 1.0,
                                &(D0->GetPtr())[0], nquad0, &wsp[0], nquad0, 0.0,
                                &outarray_d0[0], nquad0);
                }
                else
                {
                    Blas::Dgemm('N', 'N', nquad0, nquad1, nquad0, 1.0,
                                &(D0->GetPtr())[0], nquad0, &inarray[0], nquad0, 0.0,
                                &outarray_d0[0], nquad0);
                }
            }

            if (outarray_d1.size() > 0) // calculate du/dx_1
            {
                DNekMatSharedPtr D1 = m_base[1]->GetD();
                if(inarray.data() == outarray_d1.data())
                {
                    Array<OneD, NekDouble> wsp(nquad0 * nquad1);
                    Vmath::Vcopy(nquad0 * nquad1,inarray.get(),1,wsp.get(),1);
                    Blas:: Dgemm('N', 'T', nquad0, nquad1, nquad1, 1.0, &wsp[0], nquad0,
                                 &(D1->GetPtr())[0], nquad1, 0.0, &outarray_d1[0], nquad0);
                }
                else
                {
                    Blas:: Dgemm('N', 'T', nquad0, nquad1, nquad1, 1.0, &inarray[0], nquad0,
                                 &(D1->GetPtr())[0], nquad1, 0.0, &outarray_d1[0], nquad0);
                }
            }
        }
       
        
        
        // find derivative of u (inarray) at all coords points
        void StdExpansion2D::PhysTensorDerivFast(
                                                 const Array<OneD, const Array<OneD, NekDouble> >& coords,
                                                 const Array<OneD, const NekDouble>& inarray,
                                                 Array<OneD, NekDouble> &out_d0,
                                                 Array<OneD, NekDouble> &out_d1)

        {        
            //int sz = GetTotPoints();
            const int nq0 = m_base[0]->GetNumPoints();
            const int nq1 = m_base[1]->GetNumPoints();
  
            // collapse coords;
            Array<OneD, NekDouble>  collcoords(2);

            if(out_d0.size() > 0)
            {    
                for(int i = 0; i < coords[0].size(); i++)
                {
                    Array<OneD, NekDouble> tmp(2);
                    tmp[0] = coords[0][i];
                    tmp[1] = coords[1][i];
                    LocCoordToLocCollapsed(tmp, collcoords);                    
                    Array<OneD, NekDouble> wsp(nq1);
                    for(int k = 0; k < nq0; k++)
                    {
                        for (int j = 0; j < nq1; ++j)
                        {
                            wsp[j] = StdExpansion::BaryEvaluateDeriv<0>(
                                                                        collcoords[0], &inarray[0] + j * nq0);
                            
                        }
                        
                        out_d0[i] =  StdExpansion::BaryEvaluate<1>( collcoords[1], &wsp[0]);
                    } 
                }
            }
            if(out_d1.size()>0)
            {
                
                for(int i = 0; i < coords[0].size(); i++)
                {

                    Array<OneD, NekDouble> tmp(2);
                    tmp[0] = coords[0][i];
                    tmp[1] = coords[1][i];
                    LocCoordToLocCollapsed(tmp, collcoords);           
                
                    Array<OneD, NekDouble> wsp(nq1);
                    for(int k = 0; k < nq0; k++)
                    {
                        for (int j = 0; j < nq1; ++j)
                        {
                            wsp[j] = StdExpansion::BaryEvaluate<0>(
                                                                   collcoords[0], &inarray[0] + j * nq0);
                            
                        }
                        
                        out_d1[i] =  StdExpansion::BaryEvaluateDeriv<1>(collcoords[1], &wsp[0]);
                    } 
                }            
                
            }
        }
              
        //slow version
        // fast one impl as v_PhysEvalBasisGradFast() -> does not use storage space
        void StdExpansion2D::v_PhysEvalBasisGrad(
                                                 const Array<OneD, const Array<OneD, NekDouble> >coords,
                                                 Array<OneD, NekDouble> &out_eval,                    
                                                 Array<OneD, NekDouble> &out_d0,
                                                 Array<OneD, NekDouble> &out_d1,
                                                 Array<OneD, NekDouble> &out_d2)
        {
            boost::ignore_unused(out_d2);

            int tot = GetTotPoints();
                
            Array<OneD, NekDouble> physvals(tot);
            Array<OneD, NekDouble> physvalsder(tot);

            const int nq0 = m_base[0]->GetNumPoints();
            const int nq1 = m_base[1]->GetNumPoints();

            int neq = m_ncoeffs;
             

            if(out_eval.size() > 0)
            {    
                
                Array<OneD, NekDouble> wsp(nq1);
                for(int k = 0; k < neq; k++)
                {
                    Vmath::Vcopy(tot, &m_physevalall[0][k*tot], 1, &physvals[0], 1);       
                    for(int i = 0; i < tot; i++)
                    {
                        Array<OneD, NekDouble> coll1(2);
                        Array<OneD, NekDouble> coll2(2);
                        coll1[0] = coords[0][i];
                        coll1[1] = coords[1][i];
                        LocCoordToLocCollapsed(coll1, coll2);
                     
                        for (int j = 0; j < nq1; ++j)
                        {
                            wsp[j] = StdExpansion::BaryEvaluate<0>(
                                                                   coll2[0], &physvals[0] + j * nq0);
                        }
                        out_eval[i+k*tot] =  StdExpansion::BaryEvaluate<1>(coll2[1], &wsp[0]); 
                    }
                }
                
            } 

            if(out_d0.size() > 0)
            {    
                
                Array<OneD, NekDouble> wsp(nq1);
                for(int k = 0; k < neq; k++)
                {
                    Vmath::Vcopy(tot, &m_physevalall[1][k*tot], 1, &physvals[0], 1);       
                    for(int i = 0; i < tot; i++)
                    {
                        Array<OneD, NekDouble> coll1(2);
                        Array<OneD, NekDouble> coll2(2);
                        coll1[0] = coords[0][i];
                        coll1[1] = coords[1][i];
                        LocCoordToLocCollapsed(coll1, coll2);
                     
                        for (int j = 0; j < nq1; ++j)
                        {
                            wsp[j] = StdExpansion::BaryEvaluate<0>(
                                                                   coll2[0], &physvals[0] + j * nq0);
                        }
                        out_d0[i+k*tot] =  StdExpansion::BaryEvaluate<1>(coll2[1], &wsp[0]); 
                    }
                }
                
            }
            if(out_d1.size() > 0)
            {    
                
                Array<OneD, NekDouble> wsp(nq1);
                for(int k = 0; k < neq; k++)
                {
                    Vmath::Vcopy(tot, &m_physevalall[2][k*tot], 1, &physvals[0], 1);       
                    for(int i = 0; i < tot; i++)
                    {
                        Array<OneD, NekDouble> coll1(2);
                        Array<OneD, NekDouble> coll2(2);
                        coll1[0] = coords[0][i];
                        coll1[1] = coords[1][i];
                        LocCoordToLocCollapsed(coll1, coll2);
                        
                        for (int j = 0; j < nq1; ++j)
                        {
                            wsp[j] = StdExpansion::BaryEvaluate<0>(
                                                                   coll2[0], &physvals[0] + j * nq0);
                        }
                        out_d1[i+k*tot] =  StdExpansion::BaryEvaluate<1>(coll2[1], &wsp[0]); 
                    }
                }
                
            }
        }

        // create and populate storage for slow versions of physderiv
        // and physbasisderiv
        Array<OneD, Array<OneD, NekDouble> >StdExpansion2D::v_GetPhysEvalALL()
        {
            
            Array<OneD, Array<OneD, NekDouble> > ret(3);
            NekDouble nq = GetTotPoints();
           
            ret[0] = Array<OneD, NekDouble>(m_ncoeffs*nq);
            ret[1] = Array<OneD, NekDouble>(m_ncoeffs*nq);
            ret[2] = Array<OneD, NekDouble>(m_ncoeffs*nq);
            for(int i = 0; i < m_ncoeffs; i++)
            {
                Array<OneD, NekDouble> tmp(nq);
                           
                Array<OneD, NekDouble> tmp2(nq);

                Array<OneD, NekDouble> tmp3(nq);
                             
                FillMode(i, tmp);
                Vmath::Vcopy(nq, &tmp[0], 1, &ret[0][i*nq], 1);  

                PhysDeriv(0, tmp, tmp2);
                Vmath::Vcopy(nq, &tmp2[0], 1, &ret[1][i*nq], 1);  

                PhysDeriv(1, tmp, tmp3);
                Vmath::Vcopy(nq, &tmp3[0], 1, &ret[2][i*nq], 1);  

            }
            return ret;
        
        }        

        //evaluates der of multiple points given in coords(2D array) with x-coords in coords[0] and ycoords in coords[1]
        
        NekDouble StdExpansion2D::v_PhysEvaluate(
                                                 const Array<OneD, const NekDouble> &coords,
                                                 const Array<OneD, const NekDouble> &physvals)
        {
            ASSERTL2(coords[0] > -1 - NekConstants::kNekZeroTol,
                     "coord[0] < -1");
            ASSERTL2(coords[0] <  1 + NekConstants::kNekZeroTol,
                     "coord[0] >  1");
            ASSERTL2(coords[1] > -1 - NekConstants::kNekZeroTol,
                     "coord[1] < -1");
            ASSERTL2(coords[1] <  1 + NekConstants::kNekZeroTol,
                     "coord[1] >  1");

            Array<OneD, NekDouble> coll(2);
            LocCoordToLocCollapsed(coords,coll);
            //coll=coords;
            const int nq0 = m_base[0]->GetNumPoints();
            const int nq1 = m_base[1]->GetNumPoints();
            
            Array<OneD, NekDouble> wsp(nq1);
            for (int i = 0; i < nq1; ++i)
            {
                wsp[i] = StdExpansion::BaryEvaluate<0>(
                                                       coll[0], &physvals[0] + i * nq0);
            }

            return StdExpansion::BaryEvaluate<1>(coll[1], &wsp[0]);
        }

    
        NekDouble StdExpansion2D::v_PhysEvaluate(
                                                 const Array<OneD, DNekMatSharedPtr > &I,
                                                 const Array<OneD, const NekDouble> &physvals)
        {
            NekDouble val;
            int i;
            int nq0 = m_base[0]->GetNumPoints();
            int nq1 = m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> wsp1(nq1);

            // interpolate first coordinate direction
            for (i = 0; i < nq1;++i)
            {
                wsp1[i] = Blas::Ddot(nq0, &(I[0]->GetPtr())[0], 1,
                                     &physvals[i * nq0], 1);
            }

            // interpolate in second coordinate direction
            val = Blas::Ddot(nq1, I[1]->GetPtr(), 1, wsp1, 1);

            return val;
        }

        //////////////////////////////
        // Integration Methods
        //////////////////////////////

        NekDouble StdExpansion2D::Integral(const Array<OneD, const NekDouble>& inarray,
                                           const Array<OneD, const NekDouble>& w0,
                                           const Array<OneD, const NekDouble>& w1)
        {
            int i;
            NekDouble Int = 0.0;
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> tmp(nquad0 * nquad1);

            // multiply by integration constants
            for (i = 0; i < nquad1; ++i)
            {
                Vmath::Vmul(nquad0, &inarray[0] + i*nquad0, 1, w0.get(),
                            1, &tmp[0] + i*nquad0, 1);
            }

            for (i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1, &tmp[0]+ i, nquad0, w1.get(), 1,
                            &tmp[0] + i, nquad0);
            }
            Int = Vmath::Vsum(nquad0 * nquad1, tmp, 1);

            return Int;
        }

        void StdExpansion2D::BwdTrans_SumFacKernel(
                                                   const Array<OneD, const NekDouble>& base0,
                                                   const Array<OneD, const NekDouble>& base1,
                                                   const Array<OneD, const NekDouble>& inarray,
                                                   Array<OneD, NekDouble> &outarray,
                                                   Array<OneD, NekDouble> &wsp,
                                                   bool doCheckCollDir0,
                                                   bool doCheckCollDir1)
        {
            v_BwdTrans_SumFacKernel(base0, base1, inarray, outarray, wsp, doCheckCollDir0, doCheckCollDir1);
        }

        void StdExpansion2D::IProductWRTBase_SumFacKernel(
                                                          const Array<OneD, const NekDouble>& base0,
                                                          const Array<OneD, const NekDouble>& base1,
                                                          const Array<OneD, const NekDouble>& inarray,
                                                          Array<OneD, NekDouble> &outarray,
                                                          Array<OneD, NekDouble> &wsp,
                                                          bool doCheckCollDir0,
                                                          bool doCheckCollDir1)
        {
            v_IProductWRTBase_SumFacKernel(base0, base1, inarray, outarray, wsp, doCheckCollDir0, doCheckCollDir1);
        }

        void StdExpansion2D::v_LaplacianMatrixOp_MatFree(
                                                         const Array<OneD, const NekDouble> &inarray,
                                                         Array<OneD,NekDouble> &outarray,
                                                         const StdRegions::StdMatrixKey &mkey)
        {
            if (mkey.GetNVarCoeff() == 0
                &&!mkey.ConstFactorExists(StdRegions::eFactorSVVCutoffRatio))
            {
                using std::max;

                // This implementation is only valid when there are no
                // coefficients associated to the Laplacian operator
                int       nquad0  = m_base[0]->GetNumPoints();
                int       nquad1  = m_base[1]->GetNumPoints();
                int       nqtot   = nquad0*nquad1;
                int       nmodes0 = m_base[0]->GetNumModes();
                int       nmodes1 = m_base[1]->GetNumModes();
                int       wspsize = max(max(max(nqtot,m_ncoeffs),nquad1*nmodes0),nquad0*nmodes1);

                const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
                const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();

                // Allocate temporary storage
                Array<OneD,NekDouble> wsp0(4*wspsize);      // size wspsize
                Array<OneD,NekDouble> wsp1(wsp0+wspsize);   // size 3*wspsize

                if(!(m_base[0]->Collocation() && m_base[1]->Collocation()))
                {
                    // LAPLACIAN MATRIX OPERATION
                    // wsp0 = u       = B   * u_hat
                    // wsp1 = du_dxi1 = D_xi1 * wsp0 = D_xi1 * u
                    // wsp2 = du_dxi2 = D_xi2 * wsp0 = D_xi2 * u
                    BwdTrans_SumFacKernel(base0,base1,inarray,wsp0,wsp1,true,true);
                    LaplacianMatrixOp_MatFree_Kernel(wsp0, outarray, wsp1);
                }
                else
                {
                    LaplacianMatrixOp_MatFree_Kernel(inarray, outarray, wsp1);
                }
            }
            else
            {
                StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(
                                                                    inarray,outarray,mkey);
            }
        }



        void StdExpansion2D::v_HelmholtzMatrixOp_MatFree(
                                                         const Array<OneD, const NekDouble> &inarray,
                                                         Array<OneD,NekDouble> &outarray,
                                                         const StdRegions::StdMatrixKey &mkey)
        {
            if (mkey.GetNVarCoeff() == 0
                &&!mkey.ConstFactorExists(StdRegions::eFactorSVVCutoffRatio))
            {
                using std::max;

                int       nquad0  = m_base[0]->GetNumPoints();
                int       nquad1  = m_base[1]->GetNumPoints();
                int       nqtot   = nquad0*nquad1;
                int       nmodes0 = m_base[0]->GetNumModes();
                int       nmodes1 = m_base[1]->GetNumModes();
                int       wspsize = max(max(max(nqtot,m_ncoeffs),nquad1*nmodes0), nquad0*nmodes1);
                NekDouble lambda  =
                    mkey.GetConstFactor(StdRegions::eFactorLambda);

                const Array<OneD, const NekDouble>& base0 = m_base[0]->GetBdata();
                const Array<OneD, const NekDouble>& base1 = m_base[1]->GetBdata();

                // Allocate temporary storage
                Array<OneD,NekDouble> wsp0(5*wspsize);      // size wspsize
                Array<OneD,NekDouble> wsp1(wsp0 + wspsize);  // size wspsize
                Array<OneD,NekDouble> wsp2(wsp0 + 2*wspsize);// size 3*wspsize

                if (!(m_base[0]->Collocation() && m_base[1]->Collocation()))
                {
                    // MASS MATRIX OPERATION
                    // The following is being calculated:
                    // wsp0     = B   * u_hat = u
                    // wsp1     = W   * wsp0
                    // outarray = B^T * wsp1  = B^T * W * B * u_hat = M * u_hat
                    BwdTrans_SumFacKernel       (base0, base1, inarray,
                                                 wsp0, wsp2,true,true);
                    MultiplyByQuadratureMetric  (wsp0, wsp1);
                    IProductWRTBase_SumFacKernel(base0, base1, wsp1, outarray,
                                                 wsp2, true, true);

                    LaplacianMatrixOp_MatFree_Kernel(wsp0, wsp1, wsp2);
                }
                else
                {
                    MultiplyByQuadratureMetric(inarray,outarray);
                    LaplacianMatrixOp_MatFree_Kernel(inarray, wsp1, wsp2);
                }

                // outarray = lambda * outarray + wsp1
                //          = (lambda * M + L ) * u_hat
                Vmath::Svtvp(m_ncoeffs, lambda, &outarray[0], 1,
                             &wsp1[0], 1, &outarray[0], 1);
            }
            else
            {
                StdExpansion::HelmholtzMatrixOp_MatFree_GenericImpl(
                                                                    inarray,outarray,mkey);
            }
        }

    } //end namespace
} //end namespace
