///////////////////////////////////////////////////////////////////////////////
//
// File FilterGJPSmoothing.cpp
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
// Description: Gradient Jump Penalty Smoothing. 
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Filters/FilterGJPSmoothing.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/AssemblyMap/LocTraceToTraceMap.h>
                                           
namespace Nektar
{
namespace SolverUtils
{
std::string FilterGJPSmoothing::className =
        GetFilterFactory().RegisterCreatorFunction(
                "GJPSmoothing", FilterGJPSmoothing::create);

FilterGJPSmoothing::FilterGJPSmoothing(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem>        &pEquation,
    const ParamMap &pParams)
    : Filter(pSession, pEquation),
      m_index(0)
{
    // Load sampling frequency
    auto it = pParams.find("SmoothingFrequency");
    if (it == pParams.end())
    {
        m_smoothingFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(
            m_session->GetInterpreter(), it->second);
        m_smoothingFrequency = round(equ.Evaluate());
    }

    m_traceNormals = pEquation.lock()->GetTraceNormals();
}
    

FilterGJPSmoothing::~FilterGJPSmoothing()
{
}

void FilterGJPSmoothing::v_Initialise
(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
 const NekDouble &time)
{
    boost::ignore_unused(pFields,time);

    // Check this is a Mixed CG-DG expansion
    ASSERTL0(boost::iequals(m_session->GetSolverInfo("Projection"),
                            "Mixed_CG_Discontinuous"),
             "Currently the GJP is only set up for mixed CG/DG expansions");

    MultiRegions::ExpListSharedPtr      trace = pFields[0]->GetTrace(); 
    MultiRegions::AssemblyMapSharedPtr traceMap = pFields[0]->GetTraceMap(); 

    m_shapeDim = pFields[0]->GetShapeDimension();
    m_scalFwd = Array<OneD, NekDouble> (trace->GetNpoints());
    m_scalBwd = Array<OneD, NekDouble> (trace->GetNpoints());

    Array<OneD, Array<OneD, NekDouble> > factors;
    Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
        &elmtToTrace = std::dynamic_pointer_cast<MultiRegions::AssemblyMapDG>
                (traceMap)->GetElmtToTrace();
    
    MultiRegions::LocTraceToTraceMapSharedPtr locTraceToTraceMap
        = pFields[0]->GetLocTraceToTraceMap();

    std::vector<bool> leftAdjacentTrace = pFields[0]->GetLeftAdjacentTraces();

    Array<OneD, NekDouble> FwdLocTrace(locTraceToTraceMap->GetNLocTracePts());
    Array<OneD, NekDouble> BwdLocTrace(locTraceToTraceMap->GetNLocTracePts());
    Array<OneD, NekDouble> e_tmp;

    const std::shared_ptr<LocalRegions::ExpansionVector> exp = pFields[0]->GetExp();
    
    int cnt    = 0;
    int toffset_coeff = 0;
    int toffset_phys  = 0;
    int coeff_offset = 0; 
    Array<OneD, Array<OneD, NekDouble> > dbasis;
    Array<OneD, StdRegions::TraceCoeffSet >traceToCoeffMap;

    for(int e = 0; e < pFields[0]->GetExpSize(); ++e)
    {
        (*exp)[e]->NormalTraceDerivFactors(factors);

        (*exp)[e]->DerivNormalBasisOnTrace(dbasis,
                                           traceToCoeffMap);
        
        for(int n = 0; n < (*exp)[e]->GetNtraces(); ++n, ++cnt)
        {
            int nptrace = (*exp)[e]->GetTraceExp(n)->GetTotPoints();

            int eid = elmtToTrace[e][n]->GetElmtId();
            toffset_coeff = trace->GetCoeff_Offset(eid);

            e_tmp = (leftAdjacentTrace[cnt])? FwdLocTrace + toffset_phys:
                BwdLocTrace + toffset_phys; 
            
            Vmath::Vcopy(nptrace,factors[n],1,e_tmp,1);

            int nctrace = elmtToTrace[e][n]->GetNcoeffs();
            if(leftAdjacentTrace[cnt])
            {
                for(int i = 0; i < nctrace; ++i)
                {
                    std::set<int> dbset; 
                    for(auto &it : traceToCoeffMap[n][i])
                    {
                        dbset.insert(coeff_offset+ it);
                    }
                    // make a pair containing the dervative of the
                    // basis for a trace mode and teh set of
                    // coefficients it should add
                    std::pair<NekDouble, std::set<int> >  dbaseinfo(dbasis[n][i],dbset); 
                    m_fwdTraceToCoeffMap[i+toffset_coeff] = dbaseinfo;
                }
            }
            else
            {
                for(int i = 0; i < nctrace; ++i)
                {
                    std::set<int> dbset; 
                    for(auto &it : traceToCoeffMap[n][i])
                    {
                        dbset.insert(coeff_offset+it);
                    }
                    // make a pair containing the dervative of the
                    // basis for a trace mode and teh set of
                    // coefficients it should add
                    //
                    // Minus sign is to to normal need to be negated on bwd trace
                    std::pair<NekDouble, std::set<int> >  dbaseinfo(-1.0*dbasis[n][i],dbset); 
                    m_bwdTraceToCoeffMap[i+toffset_coeff] = dbaseinfo;
                }

            }
            toffset_phys  += nptrace; 
        }

        coeff_offset += (*exp)[e]->GetNcoeffs();
    }

    locTraceToTraceMap->InterpLocTracesToTrace(0,FwdLocTrace, m_scalFwd);
    locTraceToTraceMap->InterpLocTracesToTrace(0,BwdLocTrace, m_scalBwd);
}


void FilterGJPSmoothing::v_Update
(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
 const NekDouble &time)
{
    boost::ignore_unused(pFields,time);
    
    m_index++;
    if (m_index % m_smoothingFrequency > 0)
    {
        return;
    }

    int ncoeffs   = pFields[0]->GetNcoeffs();
    int nphys     = pFields[0]->GetNpoints();
    int nTracePts = pFields[0]->GetTrace()->GetTotPoints();
    int nTraceCoeffs = pFields[0]->GetTrace()->GetNcoeffs();

    Array<OneD, NekDouble> FilterCoeffs(ncoeffs);
    Array<OneD, Array<OneD, NekDouble> > deriv(3,NullNekDouble1DArray);
    for(int i = 0; i < m_shapeDim; ++i)
    {
        deriv[i] = Array<OneD, NekDouble> (nphys); 
    }

    Array<OneD, NekDouble> Fwd(nTracePts), Bwd(nTracePts); 
    Array<OneD, NekDouble> GradJumpOnTrace(nTracePts,0.0); 
    
    for(int f = 0; f < pFields.num_elements(); ++f)
    {
        // calculate derivative 
        pFields[f]->PhysDeriv(pFields[f]->GetPhys(),
                              deriv[0],deriv[1],deriv[2]);

        // Evaluate the  normal derivative jump on the trace
        for(int n = 0; n < m_shapeDim; ++n)
        {
            pFields[f]->GetFwdBwdTracePhys(deriv[n],Fwd,Bwd,true);

            // Multiply by normal and add to trace evaluation
            Vmath::Vsub(nTracePts,Fwd,1,Bwd,1,Fwd,1);
            Vmath::Vvtvp(nTracePts,Fwd,1,m_traceNormals[n],1,
                         GradJumpOnTrace,1,GradJumpOnTrace,1);
        }

        // Add filter to coeffs to inner product of phys values and
        // then bwdtrans
        pFields[f]->IProductWRTBase(pFields[f]->GetPhys(),FilterCoeffs);

        // Scale jump on fwd trace
        Vmath::Vmul(nTracePts,m_scalFwd,1,GradJumpOnTrace,1,Fwd,1);

        // Take inner product and put result into Fwd array
        pFields[f]->GetTrace()->IProductWRTBase(Fwd,Bwd); 

        // Add trace values to coeffs scaled by factor and sign. 
        for(int i = 0; i < nTraceCoeffs; ++i)
        {
            NekDouble scal = m_fwdTraceToCoeffMap[i].first; 
            for(auto &it:  m_fwdTraceToCoeffMap[i].second)
            {
                FilterCoeffs[it] += Bwd[i]*scal;
            }
        }

        // Scale jump on fwd trace
        Vmath::Vmul(nTracePts,m_scalBwd,1,GradJumpOnTrace,1,Fwd,1);

        // Take inner product and put result into Fwd array
        pFields[f]->GetTrace()->IProductWRTBase(Fwd,Bwd); 

        // Add trace values to coeffs scaled by factor and sign. 
        for(int i = 0; i < nTraceCoeffs; ++i)
        {
            NekDouble scal = m_bwdTraceToCoeffMap[i].first; 
            for(auto &it:  m_bwdTraceToCoeffMap[i].second)
            {
                FilterCoeffs[it] += Bwd[i]*scal;
            }
        }

        pFields[f]->MultiplyByInvMassMatrix(FilterCoeffs,deriv[0]);
        pFields[f]->BwdTrans(deriv[0],pFields[f]->UpdatePhys());
    }
}

void FilterGJPSmoothing::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(pFields,time);
}

bool FilterGJPSmoothing::v_IsTimeDependent()
{
    return true;
}
}
}
