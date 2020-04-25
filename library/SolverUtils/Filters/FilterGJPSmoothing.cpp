//////////////////////////////////////////////////////////////////////////////
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

    // Not sure what variable we shoudl set this up based on? 
    m_dgfield = MemoryManager<MultiRegions::DisContField>::AllocateSharedPtr
        (m_session,pFields[0]->GetGraph(), m_session->GetVariable(0),
         true,false);
    
    MultiRegions::ExpListSharedPtr      trace   = m_dgfield->GetTrace(); 
    MultiRegions::AssemblyMapSharedPtr traceMap = m_dgfield->GetTraceMap(); 

    m_expType  = m_dgfield->GetExpType();
    m_coordDim = m_dgfield->GetCoordim(0);

    // check to see normals are declared. For CG this does not happen. 
    if(m_traceNormals.size() == 0)
    {
        m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(m_coordDim);
        for(int i=0; i < m_coordDim; ++i)
        {
            m_traceNormals[i] = Array<OneD, NekDouble>(trace->GetNpoints()); 
        }
        
        m_dgfield->GetTrace()->GetNormals(m_traceNormals);
    }

    m_scalFwd = Array<OneD, NekDouble> (trace->GetNpoints());
    m_scalBwd = Array<OneD, NekDouble> (trace->GetNpoints());

    m_fwdTraceToCoeffMap = Array<OneD,
         std::set< std::pair<unsigned int, NekDouble> > >(trace->GetTotPoints());
    m_bwdTraceToCoeffMap = Array<OneD,
         std::set< std::pair<unsigned int, NekDouble> > >(trace->GetTotPoints());
    
    Array<OneD, Array<OneD, NekDouble> > factors;
    Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
        &elmtToTrace = std::dynamic_pointer_cast<MultiRegions::AssemblyMapDG>
                (traceMap)->GetElmtToTrace();
    
    MultiRegions::LocTraceToTraceMapSharedPtr locTraceToTraceMap
        = m_dgfield->GetLocTraceToTraceMap();

    std::vector<bool> leftAdjacentTrace = m_dgfield->GetLeftAdjacentTraces();

    Array<OneD, NekDouble> FwdLocTrace(locTraceToTraceMap->GetNLocTracePts());
    Array<OneD, NekDouble> BwdLocTrace(locTraceToTraceMap->GetNLocTracePts());
    Array<OneD, NekDouble> e_tmp;

    const std::shared_ptr<LocalRegions::ExpansionVector>
        exp = m_dgfield->GetExp();
    
    int cnt           = 0;
    int toffset_coeff = 0;
    int toffset_phys  = 0;
    int coeff_offset  = 0; 
    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >dbasis;
    Array<OneD, Array<OneD, Array<OneD, unsigned int> > >traceToCoeffMap;

    int  expdim = m_dgfield->GetShapeDimension();
    
    for(int e = 0; e < m_dgfield->GetExpSize(); ++e)
    {
        (*exp)[e]->NormalTraceDerivFactors(factors);

        (*exp)[e]->DerivNormalBasisOnTrace(dbasis,
                                           traceToCoeffMap);
        
        for(int n = 0; n < (*exp)[e]->GetNtraces(); ++n, ++cnt)
        {
            int nptrace = (*exp)[e]->GetTraceExp(n)->GetTotPoints();

            // Estimate h as average of adjacent elmeents
            ASSERTL1((*exp)[e]->GetTraceExp(n)->GetLeftAdjacentElementTrace()
                     != -1,
                     "Left adjacent expansion not setup");
            LocalRegions::ExpansionSharedPtr texp  = (*exp)[e]->GetTraceExp(n)->
                GetLeftAdjacentElementExp();
            int np = texp->GetTotPoints();
            NekDouble p  = texp->GetBasisNumModes(0);
            Array<OneD, NekDouble> one(np,1.0);
            NekDouble h = texp->Integral(one);
            h = pow(h,expdim);
            
            if((*exp)[e]->GetTraceExp(n)->GetRightAdjacentElementTrace() != -1)
            {
                texp  = (*exp)[e]->GetTraceExp(n)->GetRightAdjacentElementExp();
                np = texp->GetTotPoints();
                one = Array<OneD, NekDouble>(np,1.0);
                NekDouble  h2 = texp->Integral(one);
                h2 = pow(h,expdim);
                h = (h + h2)*0.5;
                NekDouble p2  = texp->GetBasisNumModes(0);
                p = (p + p2)*0.5;
            }            
            NekDouble jumpScal = 0.7*pow(p,-4.0)*h*h; 
            
            int eid = elmtToTrace[e][n]->GetElmtId();
            toffset_phys  =  trace->GetPhys_Offset(eid);
            toffset_coeff = trace->GetCoeff_Offset(eid);

            int nctrace = elmtToTrace[e][n]->GetNcoeffs();
            // Note curretly doing the same thing so can liklely remove this if
            if(leftAdjacentTrace[cnt])
            {
                Vmath::Smul(nptrace,jumpScal, factors[n], 1,
                            e_tmp = FwdLocTrace + toffset_phys,1);
                
                ASSERTL1(dbasis[n].size() == nctrace,"The dimension of the "
                    "normalderivtrace method do not match"
                         " with the size of the trace"); 
                for(int i = 0; i < nctrace; ++i)
                {
                    for(int j = 0; j < dbasis[n][i].size(); ++j)
                    {
                        std::pair<int, NekDouble>
                            dbaseinfo(traceToCoeffMap[n][i][j] + coeff_offset,
                                      dbasis[n][i][j]);
                            
                        m_fwdTraceToCoeffMap[i+toffset_coeff].insert(dbaseinfo);
                    }
                }
            }
            else
            {
                Vmath::Smul(nptrace,jumpScal, factors[n], 1,
                            e_tmp = BwdLocTrace + toffset_phys,1);

                for(int i = 0; i < nctrace; ++i)
                {
                    for(int j = 0; j < dbasis[n][i].size(); ++j)
                    {
                        std::pair<int, NekDouble> dbaseinfo
                            (traceToCoeffMap[n][i][j] + coeff_offset,
                             dbasis[n][i][j]);

                        m_bwdTraceToCoeffMap[i+toffset_coeff].insert(dbaseinfo);
                    }
                }

            }
            toffset_phys  += nptrace; 
        }

        coeff_offset += (*exp)[e]->GetNcoeffs();
    }

    locTraceToTraceMap->InterpLocTracesToTrace(0,FwdLocTrace, m_scalFwd);
    locTraceToTraceMap->InterpLocTracesToTrace(1,BwdLocTrace, m_scalBwd);
}


void FilterGJPSmoothing::v_Update
(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
 const NekDouble &time)
{
    boost::ignore_unused(time);
    
    m_index++;
    if (m_index % m_smoothingFrequency > 0)
    {
        return;
    }

    int ncoeffs   = m_dgfield->GetNcoeffs();
    int nphys     = m_dgfield->GetNpoints();
    int nTracePts = m_dgfield->GetTrace()->GetTotPoints();
    int nTraceCoeffs = m_dgfield->GetTrace()->GetNcoeffs();
    
    Array<OneD, NekDouble> FilterCoeffs(ncoeffs);
    Array<OneD, Array<OneD, NekDouble> > deriv(3,NullNekDouble1DArray);
    for(int i = 0; i < m_coordDim; ++i)
    {
        deriv[i] = Array<OneD, NekDouble> (nphys); 
    }
    
    Array<OneD, NekDouble> Fwd(nTracePts), Bwd(nTracePts); 
    Array<OneD, NekDouble> GradJumpOnTrace(nTracePts,0.0); 
    
    for(int f = 0; f < pFields.size(); ++f)
    {
        // calculate derivative 
        pFields[f]->PhysDeriv(pFields[f]->GetPhys(),
                              deriv[0],deriv[1],deriv[2]);
        
        // Evaluate the  normal derivative jump on the trace
        for(int n = 0; n < m_coordDim; ++n)
        {
            m_dgfield->GetFwdBwdTracePhys(deriv[n],Fwd,Bwd,true);
            
            // Multiply by normal and add to trace evaluation
            Vmath::Vsub(nTracePts,Fwd,1,Bwd,1,Fwd,1);
            Vmath::Vvtvp(nTracePts,Fwd,1,m_traceNormals[n],1,
                         GradJumpOnTrace,1,GradJumpOnTrace,1);
        }
        
        // Add filter to coeffs to inner product of phys values and
        // then bwdtrans
        pFields[f]->IProductWRTBase(pFields[f]->GetPhys(),FilterCoeffs);
        
        if(m_expType == MultiRegions::e1D)
        {
            // Scale jump on fwd trace
            Vmath::Vmul(nTracePts,m_scalFwd,1,GradJumpOnTrace,1,Bwd,1);
        }
        else
        {
            // Scale jump on fwd trace
            Vmath::Vmul(nTracePts,m_scalFwd,1,GradJumpOnTrace,1,Fwd,1);
            // Take inner product and put result into Fwd array
            m_dgfield->GetTrace()->IProductWRTBase(Fwd,Bwd);
        }
        
        // Add trace values to coeffs scaled by factor and sign. 
        for(int i = 0; i < nTraceCoeffs; ++i)
        {
            for(auto &it:  m_fwdTraceToCoeffMap[i])
            {
                FilterCoeffs[it.first] -= Bwd[i]*(it.second);
            }
        }
        
        if(m_expType == MultiRegions::e1D)
        {
            // Scale jump on fwd trace
            Vmath::Vmul(nTracePts,m_scalBwd,1,GradJumpOnTrace,1,Bwd,1);
        }
        else
        {
            // Scale jump on fwd trace
            Vmath::Vmul(nTracePts,m_scalBwd,1,GradJumpOnTrace,1,Fwd,1);
            
            // Take inner product and put result into Fwd array
            m_dgfield->GetTrace()->IProductWRTBase(Fwd,Bwd);
        }
        
        // Add trace values to coeffs scaled by factor and sign. 
        for(int i = 0; i < nTraceCoeffs; ++i)
        {
            for(auto &it:  m_bwdTraceToCoeffMap[i])
            {
                FilterCoeffs[it.first] -= Bwd[i]*(it.second);
            }
        }
        
        pFields[f]->MultiplyByElmtInvMass(FilterCoeffs,deriv[0]);
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
