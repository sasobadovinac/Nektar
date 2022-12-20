///////////////////////////////////////////////////////////////////////////////
//
// File: PreconCfsBRJ.h
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
// Description: PreconCfsBRJ header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_PRECONCFSBRJ
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_PRECONCFSBRJ

#include <CompressibleFlowSolver/Preconditioner/PreconCfsOp.h>

namespace Nektar
{

using namespace tinysimd;

/**
 * Block Relaxed(weighted) Jacobi iterative (BRJ) Preconditioner for CFS
 *
 */
class PreconCfsBRJ : public PreconCfsOp
{
public:
    friend class MemoryManager<PreconCfsBRJ>;

    /// Creates an instance of this class
    static PreconCfsOpSharedPtr create(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm)
    {
        PreconCfsOpSharedPtr p = MemoryManager<PreconCfsBRJ>::AllocateSharedPtr(
            pFields, pSession, vComm);
        return p;
    }

    /// Name of the class
    static std::string className;

    PreconCfsBRJ(const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                 const LibUtilities::SessionReaderSharedPtr &pSession,
                 const LibUtilities::CommSharedPtr &vComm);
    ~PreconCfsBRJ(){};

    virtual bool v_UpdatePreconMatCheck(const Array<OneD, const NekDouble> &res,
                                        const NekDouble dtLambda) override;

protected:
    int m_PreconItsStep;
    int m_BRJRelaxParam;

    Array<OneD, Array<OneD, SNekBlkMatSharedPtr>> m_PreconMatVarsSingle;
    
    unsigned int m_max_nblocks;
    unsigned int m_max_nElmtDof;
    std::vector<simd<NekSingle>, tinysimd::allocator<simd<NekSingle>>>
        m_sBlkDiagMat;
    std::vector<int> m_inputIdx;
    
    Array<OneD, SNekBlkMatSharedPtr> m_TraceJacSingle;
    TensorOfArray4D<NekSingle> m_TraceJacArraySingle;
    Array<OneD, SNekBlkMatSharedPtr> m_TraceJacDerivSingle;
    TensorOfArray4D<NekSingle> m_TraceJacDerivArraySingle;
    Array<OneD, Array<OneD, NekSingle>> m_TraceJacDerivSignSingle;
    TensorOfArray5D<NekSingle> m_TraceIPSymJacArraySingle;

    PrecType m_PreconMatStorage;

    virtual void v_InitObject() override;

private:
    virtual void v_DoPreconCfs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, NekDouble> &pInput, Array<OneD, NekDouble> &pOutput,
        const bool &flag) override;

    virtual void v_BuildPreconCfs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, const Array<OneD, NekDouble>> &intmp,
        const NekDouble time, const NekDouble lambda) override;

    void PreconBlkDiag(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray);

    template <typename DataType>
    void MinusOffDiag2Rhs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const size_t nvariables, const size_t nCoeffs,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, bool flagUpdateDervFlux,
        Array<OneD, Array<OneD, NekDouble>> &FwdFluxDeriv,
        Array<OneD, Array<OneD, NekDouble>> &BwdFluxDeriv,
        TensorOfArray3D<NekDouble> &qfield,
        TensorOfArray3D<NekDouble> &wspTrace,
        Array<OneD, Array<OneD, DataType>> &wspTraceDataType,
        const TensorOfArray4D<DataType> &TraceJacArray,
        const TensorOfArray4D<DataType> &TraceJacDerivArray,
        const Array<OneD, const Array<OneD, DataType>> &TraceJacDerivSign,
        const TensorOfArray5D<DataType> &TraceIPSymJacArray);

    template <typename TypeNekBlkMatSharedPtr>
    void AllocatePreconBlkDiagCoeff(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr>> &gmtxarray,
        const int &nscale = 1);

    /// This function creates the matrix structure for the block
    /// diagonal operator. It organizes the way that the matrix
    // is loaded, multiplied and unpacked in order to be operated
    /// by SIMD instructions. The degrees of freedom of each element
    /// are reorganized, so they are placed in the correct location
    /// to perfom SIMD instructions.
    inline void AllocateSIMDPreconBlkMatDiag(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
    {
        using vec_t = simd<NekSingle>;

        int TotMatLen         = 0;
        int TotLen            = 0;
        const auto nTotElmt   = pFields[0]->GetNumElmts();
        const auto nvariables = pFields.size();
        const auto vecwidth   = vec_t::width;

        m_max_nblocks  = 0;
        m_max_nElmtDof = 0;

        for (int ne = 0; ne < nTotElmt; ne++)
        {
            const auto nElmtDof = pFields[0]->GetNcoeffs(ne) * nvariables;
            const auto nblocks  = nElmtDof / vecwidth;
            unsigned int totblocks =
                (nElmtDof % vecwidth) ? nblocks + 1 : nblocks;

            m_max_nblocks =
                (m_max_nblocks > totblocks) ? m_max_nblocks : totblocks;
            m_max_nElmtDof =
                (m_max_nElmtDof > nElmtDof) ? m_max_nElmtDof : nElmtDof;

            TotLen += totblocks * vecwidth;
            TotMatLen += nElmtDof * totblocks;
        }

        m_sBlkDiagMat.resize(TotMatLen);
        m_inputIdx.resize(TotLen);

        // generate a index list of vector width aware mapping from
        // local coeff storage over all variables to elemental storage
        // over variables
        unsigned int ncoeffs = pFields[0]->GetNcoeffs();
        for (int ne = 0, cnt1 = 0; ne < nTotElmt; ne++)
        {
            const auto nElmtCoeff  = pFields[0]->GetNcoeffs(ne);
            const auto nElmtDof    = nElmtCoeff * nvariables;
            const auto nblocks     = nElmtDof / vecwidth;
            const auto nCoefOffset = pFields[0]->GetCoeff_Offset(ne);
            int i                  = 0;
            int i0                 = 0;
            int inOffset, j;

            for (int m = 0; m < nvariables; m++)
            {
                inOffset = m * ncoeffs + nCoefOffset;

                if (m && (vecwidth - i0 < nElmtCoeff))
                {
                    // May need to add entries from later variables to
                    // remainder of last variable if the vector width
                    // was not exact multiple of number of elemental
                    // coeffs
                    for (i = 0; i0 < vecwidth; ++i, ++i0)
                    {
                        m_inputIdx[cnt1++] = inOffset + i;
                    }
                }
		else
		{
		     i = 0; 
		}

                // load up other vectors in variable that fit into vector
                // width
                for (j = 0; (j + 1) * vecwidth < nElmtCoeff - i; ++j)
                {
                    for (i0 = 0; i0 < vecwidth; ++i0)
                    {
                        m_inputIdx[cnt1++] = inOffset + i + j * vecwidth + i0;
                    }
                }

                // load up any residual data for this variable
                for (i0 = 0, j = j * vecwidth; j < nElmtCoeff - i; ++j, ++i0)
                {
                    m_inputIdx[cnt1++] = inOffset + i + j;
                }
            }

            const auto endwidth = nElmtDof - nblocks * vecwidth;

            // fill out rest of index to match vector width with last entry
            if (endwidth)
            {
                for (i0 = endwidth; i0 < vecwidth; ++i0)
                {
                    m_inputIdx[cnt1++] = inOffset + i + j - 1;
                }
            }
            ASSERTL1(cnt1 <= TotLen, "m_inputIdx over extended");
        }
    }

    inline void AllocateNekBlkMatDig(SNekBlkMatSharedPtr &mat,
                                     const Array<OneD, unsigned int> nrow,
                                     const Array<OneD, unsigned int> ncol)
    {
        mat =
            MemoryManager<SNekBlkMat>::AllocateSharedPtr(nrow, ncol, eDIAGONAL);
        SNekMatSharedPtr loc_matNvar;
        for (size_t nelm = 0; nelm < nrow.size(); ++nelm)
        {
            int nrowsVars = nrow[nelm];
            int ncolsVars = ncol[nelm];

            loc_matNvar = MemoryManager<SNekMat>::AllocateSharedPtr(
                nrowsVars, ncolsVars, 0.0);
            mat->SetBlock(nelm, nelm, loc_matNvar);
        }
    }
};
} // namespace Nektar

#endif
