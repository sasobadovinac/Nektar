///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingGJPStabilisaton.cpp
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
// Description: Body forcing
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Forcing/ForcingGJPStabilisation.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/AssemblyMap/LocTraceToTraceMap.h>
#include <SolverUtils/EquationSystem.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
    std::string ForcingGJPStabilisaton::classNameBody = GetForcingFactory().
        RegisterCreatorFunction("GJPStabilisation",
                                ForcingGJPStabilisaton::create,
                                "Graient Jump Penalty Stablisation");
    std::string ForcingGJPStabilisaton::classNameField = GetForcingFactory().
        RegisterCreatorFunction("GJPStabilization",
                                ForcingGJPStabilisaton::create,
                                "Graient Jump Penalty Stablization");

    ForcingGJPStabilisaton::ForcingGJPStabilisaton(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::weak_ptr<EquationSystem>      &pEquation)
        : Forcing(pSession, pEquation)
    {
    }

    void ForcingGJPStabilisaton::v_InitObject
         (const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
          const unsigned int& pNumForcingFields,
          const TiXmlElement* pForce)
    {
        boost::ignore_unused(pForce);
        
        m_numForcingFields = pNumForcingFields;

        bool isHomogeneous1D;
        m_session->MatchSolverInfo("Homogeneous", "1D", isHomogeneous1D, false);
        if(isHomogeneous1D)
        {
            m_nplanes = pFields[0]->GetZIDs().size();
        }
        else
        {
            m_nplanes = 1;
        }
        
        m_traceNormals = m_equ.lock()->GetTraceNormals();
            
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
        
        m_scalFwd = Array<OneD, NekDouble> (trace->GetNpoints(),0.0);
        m_scalBwd = Array<OneD, NekDouble> (trace->GetNpoints(),0.0);
        
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
        
        Array<OneD, NekDouble> FwdLocTrace(locTraceToTraceMap->
                                           GetNLocTracePts());
        Array<OneD, NekDouble> BwdLocTrace = FwdLocTrace +
            locTraceToTraceMap->GetNFwdLocTracePts();
        Array<OneD, NekDouble> e_tmp;
        
        const std::shared_ptr<LocalRegions::ExpansionVector>
            exp = m_dgfield->GetExp();
        
        int cnt              = 0;
        int toffset_coeff    = 0;
        int fwd_offset_phys  = 0;
        int bwd_offset_phys  = 0;
        int coeff_offset     = 0; 
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >dbasis;
        Array<OneD, Array<OneD, Array<OneD, unsigned int> > >traceToCoeffMap;
        
        int expdim = m_dgfield->GetShapeDimension();
        int ncoeffs = m_dgfield->GetNcoeffs();
        
        for(int e = 0; e < m_dgfield->GetExpSize(); ++e)
        {
            (*exp)[e]->NormalTraceDerivFactors(factors);
            
            (*exp)[e]->DerivNormalBasisOnTrace(dbasis,
                                               traceToCoeffMap);
            
            for(int n = 0; n < (*exp)[e]->GetNtraces(); ++n, ++cnt)
            {
                
                // Estimate h as average of adjacent elmeents
                ASSERTL1((*exp)[e]->GetTraceExp(n)->GetLeftAdjacentElementTrace()
                         != -1,
                         "Left adjacent expansion not setup");
                LocalRegions::ExpansionSharedPtr texp  =
                    (*exp)[e]->GetTraceExp(n)->
                    GetLeftAdjacentElementExp();
                int np = texp->GetTotPoints();
                NekDouble p  = texp->GetBasisNumModes(0);
                Array<OneD, NekDouble> one(np,1.0);
                NekDouble h = texp->Integral(one);
                h = pow(h,1.0/(NekDouble)expdim);
                
                if((*exp)[e]->GetTraceExp(n)->GetRightAdjacentElementTrace()!=-1)
                {
                    texp  = (*exp)[e]->GetTraceExp(n)->
                        GetRightAdjacentElementExp();
                    np = texp->GetTotPoints();
                    one = Array<OneD, NekDouble>(np,1.0);
                    NekDouble  h2 = texp->Integral(one);
                    h2 = pow(h2,1.0/(NekDouble)expdim);
                    h = (h + h2)*0.5;
                    NekDouble p2  = texp->GetBasisNumModes(0);
                    p = (p + p2)*0.5;
                }            
                NekDouble jumpScal = 0.7*pow(p,-4.0)*h*h; 
                
                unsigned long  nctrace = elmtToTrace[e][n]->GetNcoeffs();
                int eid = elmtToTrace[e][n]->GetElmtId();
                toffset_coeff = trace->GetCoeff_Offset(eid);

                int nptrace = (*exp)[e]->GetTraceNumPoints(n);
                
                // Note curretly doing the same thing so can liklely
                // remove this if
                if(leftAdjacentTrace[cnt])
                {
                    Vmath::Smul(nptrace,jumpScal, factors[n], 1,
                                e_tmp = FwdLocTrace + fwd_offset_phys,1);
                    fwd_offset_phys += nptrace;
                    
                    // note the min is for variable p expansions
                    for(int i = 0; i < min(nctrace,dbasis[n].size()) ; ++i)
                    {
                        for(int j = 0; j < dbasis[n][i].size(); ++j)
                        {
                            int ncoeffid =traceToCoeffMap[n][i][j] +
                                coeff_offset;
                            if(ncoeffid >= ncoeffs)
                            {
                                cout << "here" << endl;
                            }
                            ASSERTL1(ncoeffid < ncoeffs, "Error in evaluating "
                                     "which coefficient is being updated");
                            std::pair<int, NekDouble> dbaseinfo(ncoeffid,
                                                               dbasis[n][i][j]);
                            m_fwdTraceToCoeffMap[i+toffset_coeff].
                                insert(dbaseinfo);
                        }
                    }
                }
                else
                {
                    Vmath::Smul(nptrace,jumpScal, factors[n], 1,
                                e_tmp = BwdLocTrace + bwd_offset_phys,1);
                    bwd_offset_phys += nptrace;
                    
                    for(int i = 0; i < min(nctrace,dbasis[n].size()) ; ++i)
                    {
                        for(int j = 0; j < dbasis[n][i].size(); ++j)
                        {
                            int ncoeffid =traceToCoeffMap[n][i][j] + coeff_offset;
                            ASSERTL1(ncoeffid < ncoeffs, "Error in evaluating "
                                     "which coefficient is being updated");
                            std::pair<int, NekDouble> dbaseinfo(ncoeffid, dbasis[n][i][j]);
                            m_bwdTraceToCoeffMap[i+toffset_coeff].
                                insert(dbaseinfo);
                        }
                    }
                    
                }
            }
            
            coeff_offset += (*exp)[e]->GetNcoeffs();
        }
        
        locTraceToTraceMap->InterpLocTracesToTrace(0,FwdLocTrace, m_scalFwd);
        locTraceToTraceMap->InterpLocTracesToTrace(1,BwdLocTrace, m_scalBwd);
    }
    
    void ForcingGJPStabilisaton::v_Apply
         (const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
          const Array<OneD, Array<OneD, NekDouble> > &inarray,
          Array<OneD, Array<OneD, NekDouble> > &outarray,
          const NekDouble &time)
    {
        boost::ignore_unused(time,pFields);
        
        int ncoeffs   = m_dgfield->GetNcoeffs();
        int nphys     = m_dgfield->GetNpoints();
        int nTracePts = m_dgfield->GetTrace()->GetTotPoints();
        int nTraceCoeffs = m_dgfield->GetTrace()->GetNcoeffs();
        
        Array<OneD, Array<OneD, NekDouble> > deriv(3,NullNekDouble1DArray);
        for(int i = 0; i < m_coordDim; ++i)
        {
            deriv[i] = Array<OneD, NekDouble> (nphys); 
        }
        
        Array<OneD, NekDouble> Fwd(nTracePts), Bwd(nTracePts), tmp; 
        Array<OneD, NekDouble> GradJumpOnTrace(nTracePts,0.0); 

        int nmax = max(ncoeffs,nphys);
        
        for(int f = 0; f < m_numForcingFields; ++f)
        {
            for(int p = 0; p < m_nplanes; ++p)
            {
                Array<OneD, NekDouble> FilterCoeffs(nmax,0.0);
                
                // calculate derivative 
                m_dgfield->PhysDeriv(inarray[f] + p*nphys,deriv[0],deriv[1],deriv[2]);
        
                // Evaluate the  normal derivative jump on the trace
                for(int n = 0; n < m_coordDim; ++n)
                {
                    m_dgfield->GetFwdBwdTracePhys(deriv[n],Fwd,Bwd,true);
                    
                    // Multiply by normal and add to trace evaluation
                    Vmath::Vsub(nTracePts,Fwd,1,Bwd,1,Fwd,1);
                    Vmath::Vvtvp(nTracePts,Fwd,1,m_traceNormals[n],1,
                                 GradJumpOnTrace,1,GradJumpOnTrace,1);
                }
                
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
                
                m_dgfield->MultiplyByElmtInvMass(FilterCoeffs,deriv[0]);
                m_dgfield->BwdTrans(deriv[0],FilterCoeffs);
                
                Vmath::Vadd(nphys,outarray[f]+p*nphys,1,FilterCoeffs,1,
                            tmp = outarray[f]+p*nphys,1);
            }
        }
    }
}
}
