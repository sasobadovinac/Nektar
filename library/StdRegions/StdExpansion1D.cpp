///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion1D.cpp
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
// which are common to 1d expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion1D.h>

namespace Nektar
{
    namespace StdRegions
    {

        StdExpansion1D::StdExpansion1D()
        {
        }

        StdExpansion1D::StdExpansion1D(int numcoeffs, const LibUtilities::BasisKey &Ba):
            StdExpansion(numcoeffs,1,Ba)
        {
        }

        StdExpansion1D::StdExpansion1D(const StdExpansion1D &T):StdExpansion(T)
        {
        }

        StdExpansion1D::~StdExpansion1D()
        {
        }


        //----------------------------
        // Differentiation Methods
        //-----------------------------
        // find derivative of u (inarray) at all coords points
        NekDouble StdExpansion1D::PhysTensorDerivFast(
            const Array<OneD, NekDouble> &coord,
            const Array<OneD, const NekDouble> &inarray,
            Array<OneD, NekDouble> &out_d0)
        {
            return StdExpansion::BaryEvaluate<0, true>(coord[0], &inarray[0], out_d0[0]);
        }

        // find derivative of u (inarray) at all quad points
        void StdExpansion1D::PhysTensorDeriv(const Array<OneD, const NekDouble>& inarray,
                                             Array<OneD, NekDouble>& outarray)
        {
            int nquad = GetTotPoints();
            DNekMatSharedPtr D = m_base[0]->GetD();

            if( inarray.data() == outarray.data())
            {
                Array<OneD, NekDouble> wsp(nquad);
                CopyArray(inarray, wsp);
                Blas::Dgemv('N',nquad,nquad,1.0,&(D->GetPtr())[0],nquad,
                            &wsp[0],1,0.0,&outarray[0],1);
            }
            else
            {
                Blas::Dgemv('N',nquad,nquad,1.0,&(D->GetPtr())[0],nquad,
                            &inarray[0],1,0.0,&outarray[0],1);
            }
        }


        Array<OneD, Array<OneD, NekDouble> >StdExpansion1D::v_GetPhysEvalALL()
        {
            Array<OneD, Array<OneD, NekDouble> > ret(2);
            NekDouble nq = GetTotPoints();
            ret[0] = Array<OneD, NekDouble>(m_ncoeffs*nq);
            ret[1] = Array<OneD, NekDouble>(m_ncoeffs*nq);
            for(int i = 0; i < m_ncoeffs; i++)
            {
                Array<OneD, NekDouble> tmp(nq);
                           
                Array<OneD, NekDouble> tmp2(nq);

                FillMode(i, tmp);
                Vmath::Vcopy(nq, &tmp[0], 1, &ret[0][i*nq], 1);  

                PhysDeriv(0, tmp, tmp2);
                Vmath::Vcopy(nq, &tmp2[0], 1, &ret[1][i*nq], 1);  

            }
            return ret;

        }       



        NekDouble StdExpansion1D::v_PhysEvaluate(
                                                 const Array<OneD, const NekDouble>& Lcoord,
                                                 const Array<OneD, const NekDouble>& physvals)
        {
            ASSERTL2(Lcoord[0] >= -1 - NekConstants::kNekZeroTol,"Lcoord[0] < -1");
            ASSERTL2(Lcoord[0] <=  1 + NekConstants::kNekZeroTol,"Lcoord[0] >  1");

            return StdExpansion::BaryEvaluate<0>(Lcoord[0], &physvals[0]);
        }

    NekDouble StdExpansion1D::v_PhysEvaluate(
        const Array<OneD, NekDouble> coord,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &out_d0,
        Array<OneD, NekDouble> &out_d1,
        Array<OneD, NekDouble> &out_d2)
    {
        boost::ignore_unused(coord, inarray, out_d0, out_d1, out_d2);

        return 0;
    }

        Array< OneD, NekDouble> StdExpansion1D::v_PhysEvaluateBasis(
                                                 const Array<OneD, const Array<OneD, NekDouble> >coords,
                                                 const Array<OneD, Array<OneD, NekDouble> > storage,
                                                 int mode)
        {
            int tot = GetTotPoints();
   
            Array<OneD, NekDouble> physall(tot), ret(coords[0].size());
            
            Vmath::Vcopy(tot, &storage[0][mode*tot], 1, &physall[0], 1);
            for(int i = 0; i < coords[0].size(); i++)
            {
                Array<OneD, NekDouble> ctemp(1);
                ctemp[0] = coords[0][i];
                
                ret[i] = v_PhysEvaluate(ctemp, physall);
            }
            return ret;
        }
  
        void StdExpansion1D::v_PhysEvalBasisGrad(
                                                 const Array<OneD, const Array<OneD, NekDouble> >coords,
                                                 Array<OneD, Array<OneD, NekDouble> > storage,
                                                 Array<OneD, NekDouble> &out_eval,                    
                                                 Array<OneD, NekDouble> &out_d0,
                                                 Array<OneD, NekDouble> &out_d1,
                                                 Array<OneD, NekDouble> &out_d2)
        {
            boost::ignore_unused(out_d1, out_d2);
            int tot = GetTotPoints();
   
            Array<OneD, NekDouble> physall(tot);
            int neq = m_ncoeffs;


            if(out_eval.size() > 0)
            {    
                for(int k = 0; k < neq; k++)
                {
                    Vmath::Vcopy(tot, &storage[0][k*tot], 1, &physall[0], 1);
                    for(int i = 0; i < coords[0].size(); i++)
                    {
                        Array<OneD, NekDouble> ctemp(1);
                        ctemp[0] = coords[0][i];
                          
                        out_eval[i+k*(coords[0].size())] = v_PhysEvaluate(ctemp, physall);              }
                }
                
            } 


            if(out_d0.size() > 0)
            {    
                
                std::cout<<"\n here\n\n";
                for(int k = 0; k < neq; k++)
                {
                    Vmath::Vcopy(tot, &storage[1][k*tot], 1, &physall[0], 1);
                    for(int i = 0; i < tot; i++)
                    {
                        Array<OneD, NekDouble> ctemp(1);
                        ctemp[0] = coords[0][i];
                        
                        out_d0[i+k*tot] = v_PhysEvaluate(ctemp, physall);              
                    }
                }
              
                
            }            

        }

    }//end namespace
}//end namespace
