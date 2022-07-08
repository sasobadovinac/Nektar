///////////////////////////////////////////////////////////////////////////////
//
// File:  PreconCfsBRJ.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description:  PreconCfsBRJ definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Timer.h>
#include <CompressibleFlowSolver/Preconditioner/PreconCfsBRJ.h>
#include <LibUtilities/BasicUtils/Timer.h>

using namespace std;

namespace Nektar
{
/**
 * @class  PreconCfsBRJ
 *
 * Solves a linear system using iterative methods.
 */
std::string PreconCfsBRJ::className =
    GetPreconCfsOpFactory().RegisterCreatorFunction(
        "PreconCfsBRJ", PreconCfsBRJ::create,
        "Block Relaxed Jacobi Preconditioner for CFS.");

PreconCfsBRJ::PreconCfsBRJ(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::CommSharedPtr &vComm)
    : PreconCfsOp(pFields, pSession, vComm)
{
    pSession->LoadParameter("PreconItsStep", m_PreconItsStep, 7);
    pSession->LoadParameter("BRJRelaxParam", m_BRJRelaxParam, 1.0);

    int nvariables     = pFields.size();
    m_PreconMatStorage = eDiagonal;

    m_PreconMatVarsSingle = TensorOfArray2D<SNekBlkMatSharedPtr>(nvariables);
    for (int i = 0; i < nvariables; i++)
    {
        m_PreconMatVarsSingle[i] = Array<OneD, SNekBlkMatSharedPtr>(nvariables);
    }
    AllocatePreconBlkDiagCoeff(pFields, m_PreconMatVarsSingle);

#ifdef SIMD
    AllocateSIMDPreconBlkMatDiag(pFields);
#else
    int nelmts = pFields[0]->GetNumElmts();
    int nelmtcoef;
    Array<OneD, unsigned int> nelmtmatdim(nelmts);
    for (int i = 0; i < nelmts; i++)
    {
        nelmtcoef      = pFields[0]->GetExp(i)->GetNcoeffs();
        nelmtmatdim[i] = nelmtcoef * nvariables;
    }
    AllocateNekBlkMatDig(m_PreconMatSingle, nelmtmatdim, nelmtmatdim);
#endif
}

void PreconCfsBRJ::v_InitObject()
{
    PreconCfsOp::v_InitObject();
}

void PreconCfsBRJ::v_DoPreconCfs(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray,
    const bool &flag)
{
    boost::ignore_unused(flag);

    int nBRJIterTot = m_PreconItsStep;
    if (0 == nBRJIterTot)
    {
        DoNullPrecon(inarray, outarray, flag);
    }
    else
    {
        const NekDouble BRJParam   = m_BRJRelaxParam;
        const NekDouble OmBRJParam = 1.0 - BRJParam;

        unsigned int nvariables = pFields.size();
        unsigned int npoints    = pFields[0]->GetNcoeffs();
        unsigned int ntotpnt    = inarray.size();

        ASSERTL0(nvariables * npoints == ntotpnt,
                 "nvariables*npoints!=ntotpnt in PreconCoeff");

        Array<OneD, NekDouble> rhs(ntotpnt);

        Array<OneD, NekDouble> outN(ntotpnt);
        Array<OneD, NekDouble> outTmp(ntotpnt);
        Array<OneD, Array<OneD, NekDouble>> rhs2d(nvariables);
        Array<OneD, Array<OneD, NekDouble>> out_2d(nvariables);
        Array<OneD, Array<OneD, NekDouble>> outTmp_2d(nvariables);
        for (int m = 0; m < nvariables; m++)
        {
            int moffset  = m * npoints;
            rhs2d[m]     = rhs + moffset;
            out_2d[m]    = outarray + moffset;
            outTmp_2d[m] = outTmp + moffset;
            pFields[m]->MultiplyByMassMatrix(inarray + moffset, rhs2d[m]);
        }

        int nphysic   = pFields[0]->GetNpoints();
        int nTracePts = pFields[0]->GetTrace()->GetNpoints();
        TensorOfArray3D<NekDouble> qfield(m_spacedim);
        for (int i = 0; i < m_spacedim; i++)
        {
            qfield[i] = Array<OneD, Array<OneD, NekDouble>>(nvariables);
            for (int j = 0; j < nvariables; j++)
            {
                qfield[i][j] = Array<OneD, NekDouble>(nphysic);
            }
        }
        int ntmpTrace = 4 + 2 * m_spacedim;
        TensorOfArray3D<NekDouble> tmpTrace(ntmpTrace);
        for (int i = 0; i < ntmpTrace; i++)
        {
            tmpTrace[i] = Array<OneD, Array<OneD, NekDouble>>(nvariables);
            for (int j = 0; j < nvariables; j++)
            {
                tmpTrace[i][j] = Array<OneD, NekDouble>(nTracePts);
            }
        }
        Array<OneD, Array<OneD, NekDouble>> FwdFluxDeriv(nvariables);
        Array<OneD, Array<OneD, NekDouble>> BwdFluxDeriv(nvariables);
        for (int j = 0; j < nvariables; j++)
        {
            FwdFluxDeriv[j] = Array<OneD, NekDouble>(nTracePts);
            BwdFluxDeriv[j] = Array<OneD, NekDouble>(nTracePts);
        }

        bool flagUpdateDervFlux = false;

        const int nwspTraceDataType = nvariables + 1;
        Array<OneD, Array<OneD, NekSingle>> wspTraceDataType(nwspTraceDataType);
        for (int m = 0; m < nwspTraceDataType; m++)
        {
            wspTraceDataType[m] = Array<OneD, NekSingle>(nTracePts);
        }

        LibUtilities::Timer timer;
        timer.Start();
        PreconBlkDiag(pFields, rhs, outarray);
        timer.Stop();
        timer.AccumulateRegion("PreconCfsBRJ::PreconBlkDiag", 2);

        for (int nrelax = 0; nrelax < nBRJIterTot - 1; nrelax++)
        {
            Vmath::Smul(ntotpnt, OmBRJParam, outarray, 1, outN, 1);

            timer.Start();
            MinusOffDiag2Rhs(
                pFields, nvariables, npoints, rhs2d, out_2d, flagUpdateDervFlux,
                FwdFluxDeriv, BwdFluxDeriv, qfield, tmpTrace, wspTraceDataType,
                m_TraceJacArraySingle, m_TraceJacDerivArraySingle,
                m_TraceJacDerivSignSingle, m_TraceIPSymJacArraySingle);
            timer.Stop();
            timer.AccumulateRegion("PreconCfsBRJ::MinusOffDiag2Rhs", 2);
            
            timer.Start();
            PreconBlkDiag(pFields, outarray, outTmp);
            timer.Stop();
            timer.AccumulateRegion("PreconCfsBRJ::PreconBlkDiag", 2);

            Vmath::Svtvp(ntotpnt, BRJParam, outTmp, 1, outN, 1, outarray, 1);
        }
    }
}

void PreconCfsBRJ::v_BuildPreconCfs(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, const Array<OneD, NekDouble>> &intmp,
    const NekDouble time, const NekDouble lambda)
{
    if (0 < m_PreconItsStep)
    {
        SNekBlkMatSharedPtr PreconMatSingle;
#ifdef SIMD
        using vec_t = simd<NekSingle>;
        int nvariables = pFields.size();
        int nelmts = pFields[0]->GetNumElmts();
        Array<OneD, unsigned int> matdim(nelmts);
        for (int i = 0; i < nelmts; i++)
        {
            matdim[i] = pFields[0]->GetExp(i)->GetNcoeffs() * nvariables;
        }
#ifdef OLDINIT
        AllocateNekBlkMatDig(m_PreconMatSingle, matdim, matdim);
        PreconMatSingle = m_PreconMatSingle;
#else
        AllocateNekBlkMatDig(PreconMatSingle, matdim, matdim);
#endif
#else
        PreconMatSingle = m_PreconMatSingle;
#endif
        m_operator.DoCalcPreconMatBRJCoeff(
            intmp, m_PreconMatVarsSingle, PreconMatSingle, m_TraceJacSingle,
            m_TraceJacDerivSingle, m_TraceJacDerivSignSingle,
            m_TraceJacArraySingle, m_TraceJacDerivArraySingle,
            m_TraceIPSymJacArraySingle);

        if (m_verbose && m_root)
        {
            cout << "     ## CalcuPreconMat " << endl;
        }

#ifdef SIMD // copy matrix to simd layout
#ifndef OLDINIT
        // load matrix 
        int cnt  = 0;
        int cnt1 = 0;
        const auto vecwidth = vec_t::width;

        alignas(vec_t::alignment) std::array<NekSingle, vec_t::width> tmp;

        for (int ne = 0; ne < nelmts; ne++)
        {
            const auto nElmtDof    = matdim[ne]; 
            const auto nblocks     = nElmtDof/vecwidth;

            const NekSingle *mmat = PreconMatSingle->
                GetBlockPtr(ne,ne)->GetRawPtr();
            /// Copy array into column major blocks of vector width
            for(int i1 = 0; i1 < nblocks; ++i1)
            {
                for(int j = 0; j < nElmtDof; ++j)
                {
                    for(int i = 0; i < vecwidth; ++i)
                    {
                        tmp[i] = mmat[j + (i1*vecwidth + i)*nElmtDof];
                    }
                    // store contiguous vec_t array. 
                    m_sBlkDiagMat[cnt++].load(tmp.data()); 
                }
            }

            const auto endwidth = nElmtDof - nblocks*vecwidth; 

            // end rows that do not fit into vector widths
            if(endwidth)
            {
                for(int j = 0; j < nElmtDof; ++j)
                {
                    for(int i = 0; i < endwidth; ++i)
                    {
                        tmp[i] = mmat[j + (nblocks*vecwidth + i)*nElmtDof];
                    }

                    for(int i = endwidth; i < vecwidth; ++i)
                    {
                        tmp[i] = 0.0;
                    }
                    m_sBlkDiagMat[cnt++].load(tmp.data()); 
                }
            }
        }
#endif
#endif
    }

    m_BndEvaluateTime   = time;
    m_DtLambdaPreconMat = lambda;

    m_CalcPreconMatFlag  = false;
    m_PreconTimesCounter = 1;

}

bool PreconCfsBRJ::UpdatePreconMatCheck(const Array<OneD, const NekDouble> &res,
                                        const NekDouble dtLambda)
{
    boost::ignore_unused(res);

    bool flag = false;

    if (m_CalcPreconMatFlag || (m_DtLambdaPreconMat != dtLambda))
    {
        flag = true;
    }

    if (m_PreconMatFreezNumb < m_PreconTimesCounter)
    {
        flag = true;
    }

    m_CalcPreconMatFlag = flag;
    return flag;
}

void PreconCfsBRJ::PreconBlkDiag(
                                 const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                                 const Array<OneD, NekDouble> &inarray,
                                 Array<OneD, NekDouble> &outarray)
{
    unsigned int nvariables = pFields.size();

    int nTotElmt = pFields[0]->GetNumElmts();
                                                          
    unsigned int ncoeffs    = pFields[0]->GetNcoeffs();

#define NTIME  1
    //#define NTIME1 1

    LibUtilities::Timer  timer1;

#ifdef SIMD
    using vec_t = simd<NekSingle>;
    const auto vecwidth = vec_t::width;

#ifdef OLDINIT // older intialisation 
    static unsigned int m_max_nblocks = 0;
    static unsigned int m_max_nElmtDof = 0;
    
    // remapped Precon Matrix
    //static std::vector<vec_t, tinysimd::allocator<vec_t>> m_sBlkDiagMat;
    //static std::vector<int> m_inputIdx; 


    if(m_max_nblocks == 0)
    {
        // will need to be a float in the end. 
        alignas(vec_t::alignment) std::array<NekSingle, vec_t::width> tmp;
        int TotMatLen = 0;
        int TotLen = 0;
        
        for (int ne = 0; ne < nTotElmt; ne++)
        {
            const auto nElmtDof    = pFields[0]->GetNcoeffs(ne)*nvariables;
            const auto nblocks     = nElmtDof/vecwidth;
            const auto totblocks   = (nElmtDof%vecwidth)? nblocks+1: nblocks; 

            m_max_nblocks  = max(m_max_nblocks,totblocks);
            m_max_nElmtDof = max(m_max_nElmtDof,nElmtDof);

            TotLen    += totblocks*vecwidth; 
            TotMatLen += nElmtDof*totblocks;
        }

        m_sBlkDiagMat.resize(TotMatLen);
        m_inputIdx.resize(TotLen);
       
        m_sBlkDiagMat[0].load(tmp.data()); 	
        // load matrix 
        int cnt  = 0;
        int cnt1 = 0;
        for (int ne = 0; ne < nTotElmt; ne++)
        {
            const auto nElmtCoeff  = pFields[0]->GetNcoeffs(ne); 
            const auto nElmtDof    = nElmtCoeff*nvariables;
            const auto nblocks     = nElmtDof/vecwidth;

            const NekSingle *mmat = m_PreconMatSingle->
                                GetBlockPtr(ne,ne)->GetRawPtr();

            /// Copy array into column major blocks of vector width
            for(int i1 = 0; i1 < nblocks; ++i1)
            {
                for(int j = 0; j < nElmtDof; ++j)
                {
                    for(int i = 0; i < vecwidth; ++i)
                    {
                        tmp[i] = mmat[j + (i1*vecwidth + i)*nElmtDof];
                    }
                    // store contiguous vec_t array. 
                    m_sBlkDiagMat[cnt++].load(tmp.data()); 
                }
            }

            unsigned int endwidth = nElmtDof - nblocks*vecwidth; 

            // end rows that do not fit into vector widths
            if(endwidth)
            {
                for(int j = 0; j < nElmtDof; ++j)
                {
                    for(int i = 0; i < endwidth; ++i)
                    {
                        tmp[i] = mmat[j + (nblocks*vecwidth + i)*nElmtDof];
                    }

                    for(int i = endwidth; i < vecwidth; ++i)
                    {
                        tmp[i] = 0.0;
                    }
                    m_sBlkDiagMat[cnt++].load(tmp.data()); 
                }
            }

            const auto nCoefOffset = pFields[0]->GetCoeff_Offset(ne);
            int i  = 0;
            int i0 = 0;
            int inOffset,j; 

            for (int m = 0; m < nvariables; m++)
            {
                inOffset  = m*ncoeffs + nCoefOffset;

                if(m)
                {
                    // May need to add entries from later variables to
                    // remainder of last variable if the vector width
                    // was not exact multiple of number of elemental
                    // coeffs
                    for(i = 0; i0 < vecwidth; ++i, ++i0)
                    {
                        m_inputIdx[cnt1++]= inOffset + i;
                    }
                }
                    
                // load up other vectors in varaible that fit into vector
                // width
                for (j = 0; j < (nElmtCoeff-i)/vecwidth; j += vecwidth)
                {
                    for(i0 = 0; i0 < vecwidth; ++i0)
                    {
                        m_inputIdx[cnt1++] = inOffset + i + j + i0;
                    }
                }
                
                // load up any residaul data for this varaible
                for(i0 = 0 ; j < nElmtCoeff-i; ++j, ++i0)
                {
                    m_inputIdx[cnt1++] = inOffset + i + j; 
                }
            }
            
            endwidth = nElmtDof - nblocks*vecwidth; 

            // fill out rest of index to match vector width with last entry
            if(endwidth)
            {
                for( ; i0 < vecwidth; ++i0)
                {
                    m_inputIdx[cnt1++] = inOffset + i + j - 1;
                }
            }
        }
    }
#endif

    // vectorized matrix multiply 
    for(int t = 0; t < NTIME; ++t) // timing loop
    {
    std::vector<vec_t, tinysimd::allocator<vec_t>> Sinarray (m_max_nblocks); 
    std::vector<vec_t, tinysimd::allocator<vec_t>> Soutarray(m_max_nElmtDof);
    //std::vector<vec_t, tinysimd::allocator<vec_t>> tmp;

    alignas(vec_t::alignment)  std::array<NekSingle, vec_t::width> tmp;

    for (int ne = 0, cnt = 0, icnt = 0, icnt1 = 0; ne < nTotElmt; ne++)
    {
        const auto nElmtCoef   = pFields[0]->GetNcoeffs(ne);
        const auto nElmtDof    = nElmtCoef*nvariables;
        const auto nblocks     = (nElmtDof%vecwidth)?
                nElmtDof/vecwidth + 1:  nElmtDof/vecwidth;

#ifdef NTIME1 // inner timing test
        timer1.Start(); 
        int icnt_sav = icnt; 
        for(int t = 0; t < NTIME1; ++t) // timing loop
        {
        icnt  = icnt_sav; 
#endif
        // gather data into blocks - could probably be done with a
        // gather call? can be replaced with a gather op when working
        for (int j = 0; j < nblocks; ++j, icnt += vecwidth)
        {
            for(int i = 0; i < vecwidth; ++i)
            {
                tmp[i] = inarray[m_inputIdx[icnt + i]]; 
            }
            
            Sinarray[j].load(tmp.data());
        }
#ifdef NTIME1 // inner timing test
        }
        timer1.Stop();
        timer1.AccumulateRegion("PreconCfsBRJ: Load");
        

        timer1.Start(); 
        int cnt_sav = cnt; 
        for(int t = 0; t < NTIME1; ++t) // timing loop
        {
            cnt  = cnt_sav; 
#endif
        // Do matrix multiply
        // first block just needs multiplying
        vec_t in = Sinarray[0];
        for (int i = 0; i < nElmtDof; ++i)
        {
            Soutarray[i] = m_sBlkDiagMat[cnt++] * in;
        }

        // rest of blocks are multiply add operations;
        for(int n = 1; n < nblocks; ++n)
        {
            in = Sinarray[n];
            for (int i = 0; i < nElmtDof; ++i)
            {
                Soutarray[i].fma(m_sBlkDiagMat[cnt++],in);
            }
        }
#ifdef NTIME1 // inner timing test
        }
        timer1.Stop();
        timer1.AccumulateRegion("PreconCfsBRJ: Mult");       
#endif
        
#if 1
#ifdef NTIME1 // inner timing test
        timer1.Start(); 
         // get block aligned index for this expansion

        int icnt1_sav = icnt1; 
        for(int t = 0; t < NTIME1; ++t) // timing loop
        {
        icnt1 = icnt1_sav; 
#endif
        NekSingle val; 
        for (int i = 0; i < nElmtDof; ++i)
        {
             // Get hold of datak
            Soutarray[i].store(tmp.data());
            
            // Sum vector width
            val = tmp[0];
            for(int j = 1; j < vecwidth; ++j)
            {
                val += tmp[j]; 
            }
            // put data into outarray 
            outarray[m_inputIdx[icnt1+i]] = NekDouble(val); 
        }
#ifdef NTIME1 // inner timing test
        }
        timer1.Stop();
        timer1.AccumulateRegion("PreconCfsBRJ: Unpack");
#endif
        icnt1 += nblocks*vecwidth;

#else

#ifdef NTIME1 // inner timing test
        timer1.Start(); 
        for(int t = 0; t < NTIME1; ++t) // timing loop
        {
#endif
        // sum vector and unpack data
        NekSingle val; 
        const auto nCoefOffset = pFields[0]->GetCoeff_Offset(ne);
        for (int m = 0, cnt1 = 0; m < nvariables; m++)
        {
            int inOffset = m*ncoeffs + nCoefOffset;

            for (int i = 0; i < nElmtCoef; ++i)
            {
                Soutarray[cnt1++].store(tmp.data());
                
                val = tmp[0];
                for(int j = 1; j < vecwidth; ++j)
                {
                    val += tmp[j]; 
                }
                outarray[inOffset + i] = NekDouble(val); 
            }
        }
#ifdef NTIME1 // inner timing test
        }
        timer1.Stop();
        timer1.AccumulateRegion("PreconCfsBRJ: Unpack");
#endif
#endif
    }
    }
#else // master implementation 

    for(int t = 0; t < NTIME; ++t)
    {
    unsigned int ncoeffsVar = nvariables * ncoeffs;
    Array<OneD, NekSingle> Sinarray(ncoeffsVar);
    Array<OneD, NekSingle> Soutarray(ncoeffsVar);


    NekVector<NekSingle> tmpVect(ncoeffsVar, Sinarray, eWrapper);
    NekVector<NekSingle> outVect(ncoeffsVar, Soutarray, eWrapper);
    for (int m = 0; m < nvariables; m++)
    {
        int nVarOffset = m * ncoeffs;
        for (int ne = 0; ne < nTotElmt; ne++)
        {
            int nCoefOffset = pFields[0]->GetCoeff_Offset(ne);
            int nElmtCoef   = pFields[0]->GetNcoeffs(ne);
            int inOffset    = nVarOffset + nCoefOffset;
            int outOffset   = nCoefOffset * nvariables + m * nElmtCoef;
            for (int i = 0; i < nElmtCoef; i++)
            {
                Sinarray[outOffset + i] = NekSingle(inarray[inOffset + i]);
            }
        }
    }

    outVect = (*m_PreconMatSingle) * tmpVect;

    for (int m = 0; m < nvariables; m++)
    {
        int nVarOffset = m * ncoeffs;
        for (int ne = 0; ne < nTotElmt; ne++)
        {
            int nCoefOffset = pFields[0]->GetCoeff_Offset(ne);
            int nElmtCoef   = pFields[0]->GetNcoeffs(ne);
            int inOffset    = nVarOffset + nCoefOffset;
            int outOffset   = nCoefOffset * nvariables + m * nElmtCoef;
            for (int i = 0; i < nElmtCoef; i++)
            {
                outarray[inOffset + i] = NekDouble(Soutarray[outOffset + i]);
            }
        }
    }
    }
#endif
}

template <typename DataType>
void PreconCfsBRJ::MinusOffDiag2Rhs(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const int nvariables, const int nCoeffs,
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, bool flagUpdateDervFlux,
    Array<OneD, Array<OneD, NekDouble>> &FwdFluxDeriv,
    Array<OneD, Array<OneD, NekDouble>> &BwdFluxDeriv,
    TensorOfArray3D<NekDouble> &qfield, TensorOfArray3D<NekDouble> &wspTrace,
    Array<OneD, Array<OneD, DataType>> &wspTraceDataType,
    const TensorOfArray4D<DataType> &TraceJacArray,
    const TensorOfArray4D<DataType> &TraceJacDerivArray,
    const Array<OneD, const Array<OneD, DataType>> &TraceJacDerivSign,
    const TensorOfArray5D<DataType> &TraceIPSymJacArray)
{
    boost::ignore_unused(flagUpdateDervFlux, qfield, TraceJacDerivArray,
                         TraceJacDerivSign, FwdFluxDeriv, BwdFluxDeriv,
                         TraceIPSymJacArray);

    int nTracePts = pFields[0]->GetTrace()->GetNpoints();
    int npoints   = pFields[0]->GetNpoints();
    int nDim      = m_spacedim;

    Array<OneD, Array<OneD, NekDouble>> outpnts(nvariables);
    for (int i = 0; i < nvariables; i++)
    {
        outpnts[i] = Array<OneD, NekDouble>(npoints, 0.0);
        pFields[i]->BwdTrans(outarray[i], outpnts[i]);
    }

    // Store forwards/backwards space along trace space
    Array<OneD, Array<OneD, NekDouble>> Fwd;
    Array<OneD, Array<OneD, NekDouble>> Bwd;
    Array<OneD, Array<OneD, NekDouble>> FwdFlux;
    Array<OneD, Array<OneD, NekDouble>> BwdFlux;
    TensorOfArray3D<NekDouble> numDerivBwd(nDim);
    TensorOfArray3D<NekDouble> numDerivFwd(nDim);
    int indexwspTrace = 0;
    Fwd               = wspTrace[indexwspTrace], indexwspTrace++;
    Bwd               = wspTrace[indexwspTrace], indexwspTrace++;
    FwdFlux           = wspTrace[indexwspTrace], indexwspTrace++;
    BwdFlux           = wspTrace[indexwspTrace], indexwspTrace++;

    LibUtilities::Timer timer;
    for (int i = 0; i < nvariables; ++i)
    {
        timer.Start(); 
        pFields[i]->GetFwdBwdTracePhys(outpnts[i], Fwd[i], Bwd[i]);
        timer.Stop();
        timer.AccumulateRegion("ExpList::GetFwdBwdTracePhys", 10);
    }

    int indexwspTraceDataType = 0;
    Array<OneD, Array<OneD, DataType>> Fwdarray(nvariables);
    for (int m = 0; m < nvariables; ++m)
    {
        Fwdarray[m] = wspTraceDataType[indexwspTraceDataType],
        indexwspTraceDataType++;
    }
    Array<OneD, DataType> Fwdreslt;
    Fwdreslt = wspTraceDataType[indexwspTraceDataType], indexwspTraceDataType++;

    for (int m = 0; m < nvariables; ++m)
    {
        for (int i = 0; i < nTracePts; ++i)
        {
            Fwdarray[m][i] = DataType(Fwd[m][i]);
        }
    }
    for (int m = 0; m < nvariables; ++m)
    {
        Vmath::Zero(nTracePts, Fwdreslt, 1);
        for (int n = 0; n < nvariables; ++n)
        {
            Vmath::Vvtvp(nTracePts, TraceJacArray[0][m][n], 1, Fwdarray[n], 1,
                         Fwdreslt, 1, Fwdreslt, 1);
        }

        for (int i = 0; i < nTracePts; ++i)
        {
            FwdFlux[m][i] = NekDouble(Fwdreslt[i]);
        }
    }

    for (int m = 0; m < nvariables; ++m)
    {
        for (int i = 0; i < nTracePts; ++i)
        {
            Fwdarray[m][i] = DataType(Bwd[m][i]);
        }
    }
    for (int m = 0; m < nvariables; ++m)
    {
        Vmath::Zero(nTracePts, Fwdreslt, 1);
        for (int n = 0; n < nvariables; ++n)
        {
            Vmath::Vvtvp(nTracePts, TraceJacArray[1][m][n], 1, Fwdarray[n], 1,
                         Fwdreslt, 1, Fwdreslt, 1);
        }
        for (int i = 0; i < nTracePts; ++i)
        {
            BwdFlux[m][i] = NekDouble(Fwdreslt[i]);
        }
    }

     
    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Fill(nCoeffs, 0.0, outarray[i], 1);
        timer.Start();
        pFields[i]->AddTraceIntegralToOffDiag(FwdFlux[i], BwdFlux[i],
                                              outarray[i]);
        timer.Stop();
        timer.AccumulateRegion("ExpList::AddTraceIntegralToOffDiag", 10);
    }
    

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Svtvp(nCoeffs, -m_DtLambdaPreconMat, outarray[i], 1, inarray[i],
                     1, outarray[i], 1);
    }
}


template <typename TypeNekBlkMatSharedPtr>
void PreconCfsBRJ::AllocatePreconBlkDiagCoeff(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr>> &gmtxarray,
    const int &nscale)
{

    int nvars  = pFields.size();
    int nelmts = pFields[0]->GetNumElmts();
    int nelmtcoef;
    Array<OneD, unsigned int> nelmtmatdim(nelmts);
    for (int i = 0; i < nelmts; i++)
    {
        nelmtcoef      = pFields[0]->GetExp(i)->GetNcoeffs();
        nelmtmatdim[i] = nelmtcoef * nscale;
    }

    for (int i = 0; i < nvars; i++)
    {
        for (int j = 0; j < nvars; j++)
        {
            AllocateNekBlkMatDig(gmtxarray[i][j], nelmtmatdim, nelmtmatdim);
        }
    }
}

} // namespace Nektar
