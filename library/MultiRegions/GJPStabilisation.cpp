///////////////////////////////////////////////////////////////////////////////
//
// File GJPStabilisation.cpp
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
// Description: GJP data
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GJPStabilisation.h>

namespace Nektar
{
namespace MultiRegions
{
GJPStabilisation::GJPStabilisation(ExpListSharedPtr pField)
{
    LibUtilities::SessionReaderSharedPtr session = pField->GetSession();

    session->MatchSolverInfo("GJPStabilisation", "SemiImplicit",
                             m_useGJPSemiImplicit, false);

    // Call GetTrace on the initialising field will set up
    // DG. Store a copoy so that if we make a soft copy of
    // this class we can re-used this field for operators.
    pField->GetTrace();
    m_dgfield = pField;

    m_coordDim = m_dgfield->GetCoordim(0);
    m_traceDim = m_dgfield->GetShapeDimension() - 1;

    // set up trace normals but would be better if could use
    // equation system definition
    m_traceNormals = Array<OneD, Array<OneD, NekDouble>>(m_coordDim);
    for (int i = 0; i < m_coordDim; ++i)
    {
        m_traceNormals[i] =
            Array<OneD, NekDouble>(m_dgfield->GetTrace()->GetNpoints());
    }
    m_dgfield->GetTrace()->GetNormals(m_traceNormals);

    SetUpExpansionInfoMapForGJP(pField->GetGraph(), session->GetVariable(0));

    MultiRegions::DisContFieldSharedPtr dgfield;

    dgfield = MemoryManager<MultiRegions::DisContField>::AllocateSharedPtr(
        session, pField->GetGraph(), "GJP", true, false,
        Collections::eNoImpType, session->GetVariable(0));
    dgfield->GetLocTraceToTraceMap(m_locTraceToTraceMap);

    m_locElmtTrace = MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(
        session, *(dgfield->GetExp()), dgfield->GetGraph(), true, "GJP");

    m_scalTrace = Array<OneD, Array<OneD, NekDouble>>(m_traceDim + 1);

    const std::shared_ptr<LocalRegions::ExpansionVector> exp =
        dgfield->GetExp();

    Array<OneD, Array<OneD, NekDouble>> dfactors[3];
    Array<OneD, Array<OneD, NekDouble>> LocTrace(m_traceDim + 1);
    Array<OneD, NekDouble> e_tmp;

    for (int i = 0; i < m_traceDim + 1; ++i)
    {
        LocTrace[i] =
            Array<OneD, NekDouble>(m_locTraceToTraceMap->GetNLocTracePts());
    }

    int cnt         = 0;
    int offset_phys = 0;
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> dbasis;
    Array<OneD, Array<OneD, Array<OneD, unsigned int>>> traceToCoeffMap;

    Array<OneD, unsigned int> map, map1;
    Array<OneD, int> sign, sign1;
    NekDouble h, p;

    for (int e = 0; e < m_dgfield->GetExpSize(); ++e)
    {
        LocalRegions::ExpansionSharedPtr elmt = (*exp)[e];

        elmt->NormalTraceDerivFactors(dfactors[0], dfactors[1], dfactors[2]);

        for (int n = 0; n < elmt->GetNtraces(); ++n, ++cnt)
        {
            NekDouble jumpScal;
            elmt->TraceNormLen(n, h, p);
            ASSERTL0(boost::math::isnan(h) == false,
                     "h has a nan value when e = " +
                         boost::lexical_cast<std::string>(e) +
                         " n =" + boost::lexical_cast<std::string>(n));

            if (p == 1)
            {
                jumpScal = 0.02 * h * h;
            }
            else
            {
                jumpScal = 0.8 * pow(p + 1, -4.0) * h * h;
            }

            int nptrace = elmt->GetTraceNumPoints(n);
            elmt->GetTraceCoeffMap(n, map);

            for (int i = 0; i < m_traceDim + 1; ++i)
            {
                Vmath::Smul(nptrace, jumpScal, dfactors[i][n], 1,
                            e_tmp = LocTrace[i] + offset_phys, 1);
            }

            offset_phys += nptrace;
        }
    }

    for (int i = 0; i < m_traceDim + 1; ++i)
    {
        m_scalTrace[i] = LocTrace[i];

        if (m_traceDim > 0)
        {
            // multiply by Jacobian and quadrature points.
            m_locElmtTrace->MultiplyByQuadratureMetric(m_scalTrace[i],
                                                       m_scalTrace[i]);
        }
    }

    // Assemble list of Matrix Product
    Array<OneD, DNekMatSharedPtr> TraceMat;

    int nelmt = 1;
    Array<OneD, const LibUtilities::BasisSharedPtr> base_sav =
        dgfield->GetExp(0)->GetBase();

    dgfield->GetExp(0)->StdDerivBaseOnTraceMat(TraceMat);

    for (int n = 1; n < dgfield->GetExpSize(); ++n)
    {
        const Array<OneD, const LibUtilities::BasisSharedPtr> &base =
            dgfield->GetExp(n)->GetBase();

        int i;
        for (i = 0; i < base.size(); ++i)
        {
            if (base[i] != base_sav[i])
            {
                break;
            }
        }

        if (i == base.size())
        {
            nelmt++;
        }
        else
        {
            // save previous block of data.
            m_StdDBaseOnTraceMat.push_back(
                std::pair<int, Array<OneD, DNekMatSharedPtr>>(nelmt, TraceMat));

            // start new block
            dgfield->GetExp(n)->StdDerivBaseOnTraceMat(TraceMat);
            nelmt    = 1;
            base_sav = dgfield->GetExp(n)->GetBase();
        }
    }
    // save latest block of data.
    m_StdDBaseOnTraceMat.push_back(
        std::pair<int, Array<OneD, DNekMatSharedPtr>>(nelmt, TraceMat));
}

void GJPStabilisation::Apply(const Array<OneD, NekDouble> &inarray,
                             Array<OneD, NekDouble> &outarray,
                             const Array<OneD, NekDouble> &pUnorm,
                             NekDouble scale) const
{
    int ncoeffs    = m_dgfield->GetNcoeffs();
    int nphys      = m_dgfield->GetNpoints();
    int nTracePts  = m_dgfield->GetTrace()->GetTotPoints();
    int nLocETrace = m_locElmtTrace->GetTotPoints();

    ASSERTL1(nLocETrace == m_scalTrace[0].size(), "expect these to be similar");
    ASSERTL1(m_locElmtTrace->GetNcoeffs() <= nphys,
             "storage assumptions presume "
             "that nLocETraceCoeffs < nphys");

    Array<OneD, Array<OneD, NekDouble>> deriv(3, NullNekDouble1DArray);
    for (int i = 0; i < m_coordDim; ++i)
    {
        deriv[i] = Array<OneD, NekDouble>(nphys);
    }

    int nmax = std::max(ncoeffs, nphys);
    Array<OneD, NekDouble> FilterCoeffs(nmax, 0.0);
    Array<OneD, NekDouble> GradJumpOnTrace(nTracePts, 0.0);
    Array<OneD, NekDouble> Fwd(nTracePts), Bwd(nTracePts);

    Array<OneD, NekDouble> wsp(nLocETrace), tmp;
#if EXPLISTDATA
    Array<OneD, NekDouble> LocElmtTracePhys   = m_locElmtTrace->UpdatePhys();
    Array<OneD, NekDouble> LocElmtTraceCoeffs = m_locElmtTrace->UpdateCoeffs();
#else
    Array<OneD, NekDouble> LocElmtTracePhys(m_locElmtTrace->GetNpoints());
    Array<OneD, NekDouble> LocElmtTraceCoefs(m_locElmtTrace->GetNcoeffs()); 
#endif

    ASSERTL1(LocElmtTracePhys.size() <= nLocETrace,
             "expect this vector to be at least of size nLocETrace");

    Array<OneD, NekDouble> unorm;
    if (pUnorm == NullNekDouble1DArray)
    {
        unorm = Array<OneD, NekDouble>(nTracePts, 1.0);
    }
    else
    {
        unorm = pUnorm;
    }

    Array<OneD, NekDouble> GradJumpOnTraceBwd;
    if (m_useGJPSemiImplicit)
    {
        GradJumpOnTraceBwd = Array<OneD, NekDouble>(nTracePts);
    }

    if (m_useGJPSemiImplicit)
    {
        Vmath::Zero(nTracePts, GradJumpOnTraceBwd, 1);
    }

    // calculate derivative
    m_dgfield->PhysDeriv(inarray, deriv[0], deriv[1], deriv[2]);

    // Evaluate the  normal derivative jump on the trace
    for (int n = 0; n < m_coordDim; ++n)
    {
        m_dgfield->GetFwdBwdTracePhys(deriv[n], Fwd, Bwd, true, true);

        if (m_useGJPSemiImplicit)
        {
            // want to put Fwd vals on bwd trace and vice versa
            Vmath::Vvtvp(nTracePts, Bwd, 1, m_traceNormals[n], 1,
                         GradJumpOnTrace, 1, GradJumpOnTrace, 1);
            Vmath::Vvtvp(nTracePts, Fwd, 1, m_traceNormals[n], 1,
                         GradJumpOnTraceBwd, 1, GradJumpOnTraceBwd, 1);
        }
        else
        {
            // Multiply by normal and add to trace evaluation
            Vmath::Vsub(nTracePts, Fwd, 1, Bwd, 1, Fwd, 1);
            Vmath::Vvtvp(nTracePts, Fwd, 1, m_traceNormals[n], 1,
                         GradJumpOnTrace, 1, GradJumpOnTrace, 1);
        }
    }

    if (m_useGJPSemiImplicit)
    {
        // Need to negate Bwd case when  using Fwd normal
        Vmath::Neg(nTracePts, GradJumpOnTrace, 1);
    }

    Vmath::Vmul(nTracePts, unorm, 1, GradJumpOnTrace, 1, GradJumpOnTrace, 1);

    // Interpolate GradJumpOnTrace to Local elemental traces.
    m_locTraceToTraceMap->InterpTraceToLocTrace(0, GradJumpOnTrace, wsp);
    m_locTraceToTraceMap->UnshuffleLocTraces(0, wsp, LocElmtTracePhys);

    if (m_useGJPSemiImplicit)
    {
        // Vmath::Neg(nTracePts,GradJumpOnTraceBwd,1);
        Vmath::Vmul(nTracePts, unorm, 1, GradJumpOnTraceBwd, 1,
                    GradJumpOnTraceBwd, 1);
        m_locTraceToTraceMap->InterpTraceToLocTrace(1, GradJumpOnTraceBwd, wsp);
        m_locTraceToTraceMap->UnshuffleLocTraces(1, wsp, LocElmtTracePhys);
    }
    else
    {
        m_locTraceToTraceMap->InterpTraceToLocTrace(1, GradJumpOnTrace, wsp);
        m_locTraceToTraceMap->UnshuffleLocTraces(1, wsp, LocElmtTracePhys);
    }

    // Scale jump on trace
    Vmath::Vmul(nLocETrace, m_scalTrace[0], 1, LocElmtTracePhys, 1, wsp, 1);
    MultiplyByStdDerivBaseOnTraceMat(0, wsp, FilterCoeffs);

    for (int i = 0; i < m_traceDim; ++i)
    {
        // Scale jump on trace
        Vmath::Vmul(nLocETrace, m_scalTrace[i + 1], 1, LocElmtTracePhys, 1, wsp,
                    1);
        MultiplyByStdDerivBaseOnTraceMat(i + 1, wsp, deriv[0]);
        Vmath::Vadd(ncoeffs, deriv[0], 1, FilterCoeffs, 1, FilterCoeffs, 1);
    }

    Vmath::Svtvp(ncoeffs, scale, FilterCoeffs, 1, outarray, 1, outarray, 1);
}

void GJPStabilisation::SetUpExpansionInfoMapForGJP(
    SpatialDomains::MeshGraphSharedPtr graph, std::string variable)
{

    // check to see if already deifned and if so return
    if (graph->ExpansionInfoDefined("GJP"))
    {
        return;
    }

    const SpatialDomains::ExpansionInfoMap expInfo =
        graph->GetExpansionInfo(variable);

    SpatialDomains::ExpansionInfoMapShPtr newInfo =
        MemoryManager<SpatialDomains::ExpansionInfoMap>::AllocateSharedPtr();

    // loop over epxansion info
    for (auto expIt = expInfo.begin(); expIt != expInfo.end(); ++expIt)
    {
        std::vector<LibUtilities::BasisKey> BKeyVector;

        for (int i = 0; i < expIt->second->m_basisKeyVector.size(); ++i)
        {
            LibUtilities::BasisKey bkeyold = expIt->second->m_basisKeyVector[i];

            // Reset radauM alpha non-zero cases to radauM
            // Legendre at one order higher

            switch (bkeyold.GetPointsType())
            {
                case LibUtilities::eGaussRadauMAlpha1Beta0:
                case LibUtilities::eGaussRadauMAlpha2Beta0:
                {
                    int npts = bkeyold.GetNumPoints();

                    // const LibUtilities::PointsKey pkey(npts+1,
                    // LibUtilities::eGaussRadauMLegendre); trying
                    // npts to be consistent for tri faces
                    const LibUtilities::PointsKey pkey(
                        npts, LibUtilities::eGaussRadauMLegendre);
                    LibUtilities::BasisKey bkeynew(bkeyold.GetBasisType(),
                                                   bkeyold.GetNumModes(), pkey);
                    BKeyVector.push_back(bkeynew);
                }
                break;
                default:
                    BKeyVector.push_back(bkeyold);
                    break;
            }
        }

        (*newInfo)[expIt->first] =
            MemoryManager<SpatialDomains::ExpansionInfo>::AllocateSharedPtr(
                expIt->second->m_geomShPtr, BKeyVector);
    }

    graph->SetExpansionInfo("GJP", newInfo);
}

void GJPStabilisation::MultiplyByStdDerivBaseOnTraceMat(
    int i, Array<OneD, NekDouble> &in, Array<OneD, NekDouble> &out) const
{
    // Should probably be vectorised

    int cnt  = 0;
    int cnt1 = 0;
    for (auto &it : m_StdDBaseOnTraceMat)
    {
        int rows = it.second[i]->GetRows();
        int cols = it.second[i]->GetColumns();

        Blas::Dgemm('N', 'N', rows, it.first, cols, 1.0,
                    &(it.second[i]->GetPtr())[0], rows, &in[0] + cnt, cols, 0.0,
                    &out[0] + cnt1, rows);

        cnt += cols * it.first;
        cnt1 += rows * it.first;
    }
}
} // namespace MultiRegions
} // namespace Nektar
