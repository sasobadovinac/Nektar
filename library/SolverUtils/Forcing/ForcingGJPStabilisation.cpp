///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingGJPStabilisation.cpp
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
    std::string ForcingGJPStabilisation::classNameBody = GetForcingFactory().
        RegisterCreatorFunction("GJPStabilisation",
                                ForcingGJPStabilisation::create,
                                "Graient Jump Penalty Stablisation");
    std::string ForcingGJPStabilisation::classNameField = GetForcingFactory().
        RegisterCreatorFunction("GJPStabilization",
                                ForcingGJPStabilisation::create,
                                "Graient Jump Penalty Stablization");

    ForcingGJPStabilisation::ForcingGJPStabilisation(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::weak_ptr<EquationSystem>      &pEquation)
        : Forcing(pSession, pEquation)
    {
    }

    void ForcingGJPStabilisation::v_InitObject
         (const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
          const unsigned int& pNumForcingFields,
          const TiXmlElement* pForce)
    {
        const TiXmlElement* funcNameElmt = pForce->FirstChildElement("HSCALING");

        if (funcNameElmt)
        {
            m_hScalingStr = funcNameElmt->GetText();
        }
        else
        {
            m_hScalingStr = std::string("DEFAULT");
        }

        funcNameElmt = pForce->FirstChildElement("VELSCALING");

        if (funcNameElmt)
        {
            m_velScalingStr = funcNameElmt->GetText();
        }
        else
        {
            m_velScalingStr = std::string("DEFAULT");
        }

        if(m_session->GetComm()->GetRank() == 0)
        {
            cout << "GJP Stabilisation:" << endl;
            cout << "\t H-Scaling:        " << m_hScalingStr << endl;
            cout << "\t Velocity-Scaling: " << m_velScalingStr << endl;
        }
        
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
            
        // Not sure what variable we should set this up based on? 
        m_dgfield = MemoryManager<MultiRegions::DisContField>::AllocateSharedPtr
            (m_session,pFields[0]->GetGraph(), m_session->GetVariable(0),
             true,false);
        MultiRegions::ExpListSharedPtr      trace   = m_dgfield->GetTrace(); 
        MultiRegions::AssemblyMapSharedPtr traceMap = m_dgfield->GetTraceMap(); 
        
        m_expType  = m_dgfield->GetExpType();
        m_coordDim = m_dgfield->GetCoordim(0);
        m_traceDim = m_dgfield->GetShapeDimension() - 1;
        
        // check to see normals are declared. For CG this does not happen. 
        if(m_traceNormals.size() == 0)
        {
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(m_coordDim);
            for(int i=0; i < m_coordDim; ++i)
            {
                m_traceNormals[i] = Array<OneD,NekDouble>(trace->GetNpoints()); 
            }
            
            m_dgfield->GetTrace()->GetNormals(m_traceNormals);
        }
        
        m_scalFwd = Array<OneD, Array<OneD, NekDouble>>(m_traceDim+1);
        m_scalBwd = Array<OneD, Array<OneD, NekDouble>>(m_traceDim+1);
        for(int i = 0; i < m_traceDim +1; ++i)
        {
            m_scalFwd[i] = Array<OneD, NekDouble> (trace->GetNpoints(),0.0);
            m_scalBwd[i] = Array<OneD, NekDouble> (trace->GetNpoints(),0.0);
        }

        m_fwdTraceToCoeffMap = Array<OneD,
            std::set< std::pair<unsigned int, NekDouble> > >(trace->GetNcoeffs());
        m_bwdTraceToCoeffMap = Array<OneD,
            std::set< std::pair<unsigned int, NekDouble> > >(trace->GetNcoeffs());
        m_fwdTraceToETraceCoeffMap = Array<OneD,
            std::set< std::pair<unsigned int, int> > >(trace->GetNcoeffs());
        m_bwdTraceToETraceCoeffMap = Array<OneD,
            std::set< std::pair<unsigned int, int> > >(trace->GetNcoeffs());

        Array<OneD, Array<OneD, NekDouble> > factors, dfactors[2];
        Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
            &elmtToTrace = std::dynamic_pointer_cast<MultiRegions::AssemblyMapDG>
            (traceMap)->GetElmtToTrace();
    
        MultiRegions::LocTraceToTraceMapSharedPtr locTraceToTraceMap
            = m_dgfield->GetLocTraceToTraceMap();
        
        std::vector<bool> leftAdjacentTrace = m_dgfield->GetLeftAdjacentTraces();

        Array<OneD, Array<OneD, NekDouble>> FwdLocTrace(m_traceDim+1); 
        Array<OneD, Array<OneD, NekDouble>> BwdLocTrace(m_traceDim+1); 
        Array<OneD, NekDouble> e_tmp;

        for(int i = 0; i < m_traceDim +1; ++i)
        {
            FwdLocTrace[i] = Array<OneD, NekDouble> (locTraceToTraceMap->
                                                    GetNLocTracePts());
            BwdLocTrace[i] = FwdLocTrace[i] +
                locTraceToTraceMap->GetNFwdLocTracePts();
        }
        
        const std::shared_ptr<LocalRegions::ExpansionVector>
            exp = m_dgfield->GetExp();
        
        int cnt              = 0;
        int toffset_coeff    = 0;
        int fwd_offset_phys  = 0;
        int bwd_offset_phys  = 0;
        int coeff_offset     = 0; 
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >dbasis;
        Array<OneD, Array<OneD, Array<OneD, unsigned int> > >traceToCoeffMap;
        
        Array<OneD, unsigned int> map, map1;
        Array<OneD, int> sign,sign1;
        NekDouble h,p;

        for(int e = 0; e < m_dgfield->GetExpSize(); ++e)
        {
            LocalRegions::ExpansionSharedPtr elmt = (*exp)[e];

            elmt->NormalTraceDerivFactors(factors,dfactors[0],dfactors[1]);
            elmt->DerivNormalBasisOnTrace(dbasis, traceToCoeffMap);
            
            for(int n = 0; n < elmt->GetNtraces(); ++n, ++cnt)
            {
                NekDouble jumpScal; 
                eval_h(elmt,n,h,p);
                ASSERTL0(boost::math::isnan(h) == false,
                         "h has a nan value when e = " + 
                         boost::lexical_cast<std::string>(e) + " n =" +
                         boost::lexical_cast<std::string>(n));

                if(boost::iequals(m_hScalingStr,"H-cubed"))
                {
                    if(p==1)
                    {
                        jumpScal = 0.02*h*h*h;
                    }
                    else
                    {
                        jumpScal = 0.8*pow(p+1,-4.0)*h*h*h;
                    }
                }
                else
                {
                    if(p==1)
                    {
                        jumpScal = 0.02*h*h;
                    }
                    else
                    {
                        jumpScal = 0.8*pow(p+1,-4.0)*h*h;
                    }
                }
                
                LocalRegions::ExpansionSharedPtr telmt = elmtToTrace[e][n];
                unsigned long  nctrace = telmt->GetNcoeffs(); 
                int eid = elmtToTrace[e][n]->GetElmtId();
                toffset_coeff = trace->GetCoeff_Offset(eid);

                int nptrace = elmt->GetTraceNumPoints(n);

                int P = (m_expType == MultiRegions::e1D)? -1:
                    telmt->GetBasisNumModes(0);
                int Q = (m_expType == MultiRegions::e3D)?
                    telmt->GetBasisNumModes(1): -1;
                
                elmt->GetElmtTraceToTraceMap(n,map,sign,
                                             elmt->GetTraceOrient(n),
                                             P,Q);                 

                elmt->GetTraceToElementMap(n,map1,sign1,elmt->GetTraceOrient(n),
                                           P,Q);                 

                if(leftAdjacentTrace[cnt])
                {
                    Vmath::Smul(nptrace, jumpScal, factors[n], 1,
                    e_tmp = FwdLocTrace[0] + fwd_offset_phys,1);

                    for(int i = 0; i < m_traceDim; ++i)
                    {
                        Vmath::Smul(nptrace, jumpScal, dfactors[i][n], 1,
                                 e_tmp = FwdLocTrace[i+1] + fwd_offset_phys,1);
                    }
                    
                    fwd_offset_phys += nptrace;
                    
                    // note the min is for variable p expansions
                    for(int i = 0; i < min(nctrace,
                                           dbasis[n].size()) ; ++i)
                    {
                        NekDouble Sign = sign[i];
                        int loc = map[i]; 
                        for(int j = 0; j < dbasis[n][i].size(); ++j)
                        {
                            int ncoeffid =traceToCoeffMap[n][i][j] +
                                coeff_offset;
                            ASSERTL1(ncoeffid < m_dgfield->GetNcoeffs(),
                                     "Error in evaluating "
                                     "which coefficient is being updated");
                            std::pair<int, NekDouble>
                                dbaseinfo(ncoeffid,Sign*dbasis[n][i][j]);
                            m_fwdTraceToCoeffMap[loc+toffset_coeff].
                                insert(dbaseinfo);
                        }
                    }

                    for(int i = 0; i<min((int)nctrace,elmt->GetTraceNcoeffs(n));
                        ++i)
                    {
                        std::pair<int, NekDouble> info(map1[i] + coeff_offset,
                                                       sign1[i]);
                        m_fwdTraceToETraceCoeffMap[i+toffset_coeff].insert(info);
                    }
                    
                }
                else
                {
                    Vmath::Smul(nptrace, jumpScal, factors[n], 1,
                                e_tmp = BwdLocTrace[0] + bwd_offset_phys,1);

                    for(int i = 0; i < m_traceDim; ++i)
                    {
                        Vmath::Smul(nptrace, jumpScal, dfactors[i][n], 1,
                                 e_tmp = BwdLocTrace[i+1] + bwd_offset_phys,1);
                    }

                    bwd_offset_phys += nptrace;
                    
                    for(int i = 0; i < min(nctrace,dbasis[n].size()) ; ++i)
                    {
                        NekDouble Sign = sign[i];
                        int loc = map[i]; 
                        for(int j = 0; j < dbasis[n][i].size(); ++j)
                        {
                            int ncoeffid =traceToCoeffMap[n][i][j] + coeff_offset;
                            ASSERTL1(ncoeffid < m_dgfield->GetNcoeffs(),
                                     "Error in evaluating "
                                     "which coefficient is being updated");
                            std::pair<int, NekDouble> dbaseinfo(ncoeffid,
                                                        Sign*dbasis[n][i][j]);
                            m_bwdTraceToCoeffMap[loc+toffset_coeff].
                                insert(dbaseinfo);
                        }
                    }

                    for(int i = 0; i<min((int)nctrace,elmt->GetTraceNcoeffs(n));
                        ++i)
                    {
                        std::pair<int, NekDouble> info(map1[i] + coeff_offset,
                                                       sign1[i]);
                        m_bwdTraceToETraceCoeffMap[i+toffset_coeff].insert(info);
                    }
                }                
            }
            
            coeff_offset += elmt->GetNcoeffs();
        }

        for(int i = 0; i < m_traceDim+1; ++i)
        {
            locTraceToTraceMap->InterpLocTracesToTrace(0,FwdLocTrace[i],
                                                       m_scalFwd[i]);
            locTraceToTraceMap->InterpLocTracesToTrace(1,BwdLocTrace[i],
                                                       m_scalBwd[i]);
        }
    }
    
    void ForcingGJPStabilisation::v_Apply
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
        
        int nmax = max(ncoeffs,nphys);
        Array<OneD, NekDouble> Fwd(nTracePts), Bwd(nTracePts), tmp; 
        Array<OneD, NekDouble> FilterCoeffs(nmax);
        Array<OneD, NekDouble> GradJumpOnTrace(nTracePts); 
        Array<OneD, NekDouble> unorm(nTracePts,1.0);

        if(boost::iequals(m_velScalingStr,"NormalVelocity"))
        {
            m_dgfield->GetFwdBwdTracePhys(inarray[0],Fwd,Bwd,true,true);
            Vmath::Vmul(nTracePts,Fwd,1,m_traceNormals[0],1,unorm,1);

            // Evaluate u.n on trace
            for(int f = 1; f < m_coordDim; ++f)
            {
                m_dgfield->GetFwdBwdTracePhys(inarray[f],Fwd,Bwd,true,true);
                Vmath::Vvtvp(nTracePts,Fwd,1,m_traceNormals[f],1,
                             unorm,1,unorm,1);
            }
            
            Vmath::Vabs(nTracePts,unorm,1,unorm,1);
        }
            
        for(int f = 0; f < m_numForcingFields; ++f)
        {
            for(int p = 0; p < m_nplanes; ++p)
            {
                Vmath::Zero(nmax,FilterCoeffs,1);
                Vmath::Zero(nTracePts,GradJumpOnTrace,1);
                
                // calculate derivative 
                m_dgfield->PhysDeriv(inarray[f] + p*nphys,deriv[0],
                                     deriv[1], deriv[2]);
        
                // Evaluate the  normal derivative jump on the trace
                for(int n = 0; n < m_coordDim; ++n)
                {
                    m_dgfield->GetFwdBwdTracePhys(deriv[n],Fwd,Bwd,true,true);
                    
                    // Multiply by normal and add to trace evaluation
                    Vmath::Vsub(nTracePts,Fwd,1,Bwd,1,Fwd,1);
                    Vmath::Vvtvp(nTracePts,Fwd,1,m_traceNormals[n],1,
                                 GradJumpOnTrace,1,GradJumpOnTrace,1);
                }
                
                // Terms related to derivative normal to fwd trace
                if(m_expType == MultiRegions::e1D)
                {
                    // Scale jump on fwd trace
                    Vmath::Vmul(nTracePts,m_scalFwd[0],1,GradJumpOnTrace,1,Bwd,1);
                }
                else
                {
                    // Scale jump on fwd trace
                    Vmath::Vmul(nTracePts,m_scalFwd[0],1,GradJumpOnTrace,1,Fwd,1);
                    Vmath::Vmul(nTracePts,unorm,1,Fwd,1,Fwd,1);

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
                
                // Terms related to derivative within fwd trace
                for(int i = 0; i < m_traceDim; ++i)
                {
                    // New approach
                    // Scale jump on bwd trace
                    Vmath::Vmul(nTracePts,m_scalFwd[i+1],1,
                                GradJumpOnTrace,1,Fwd,1);
                    Vmath::Vmul(nTracePts,unorm,1,Fwd,1,deriv[0],1);

                    m_dgfield->GetTrace()->
                        IProductWRTStdDerivBase(i,deriv[0],Fwd);

                    // map to expansion
                    for(int i = 0; i < nTraceCoeffs; ++i)
                    {
                        for(auto &it:  m_bwdTraceToETraceCoeffMap[i])
                        {
                            FilterCoeffs[it.first] -= Fwd[i]*(it.second);
                        }
                    }
                }
                
                // Terms related to derivative normal to bwd trace
                if(m_expType == MultiRegions::e1D)
                {
                    // Scale jump on bwd trace
                    Vmath::Vmul(nTracePts,m_scalBwd[0],1,GradJumpOnTrace,1,Bwd,1);
                }
                else
                {
                    // Scale jump on bwd trace
                    Vmath::Vmul(nTracePts,m_scalBwd[0],1,GradJumpOnTrace,1,Fwd,1);
                    Vmath::Vmul(nTracePts,unorm,1,Fwd,1,Fwd,1);

                    // Take inner product and put result into Fwd array
                    m_dgfield->GetTrace()->IProductWRTBase(Fwd,Bwd);
                }
                
                // Add trace values to coeffs scaled by factor and
                // sign.  Note although the nornal is negated in this case
	 	// the bwd-fwd is also reversed/negated so same sign used
                for(int i = 0; i < nTraceCoeffs; ++i)
                {
                    for(auto &it:  m_bwdTraceToCoeffMap[i])
                    {
                        FilterCoeffs[it.first] -= Bwd[i]*(it.second);
                    }
                }

                // Terms related to derivative within bwd trace
                for(int i = 0; i < m_traceDim; ++i)
                {
                    // New approach
                    // Scale jump on bwd trace
                    Vmath::Vmul(nTracePts,m_scalBwd[i+1],1,
                                GradJumpOnTrace,1,Fwd,1);
                    Vmath::Vmul(nTracePts,unorm,1,Fwd,1,deriv[0],1);

                    m_dgfield->GetTrace()->
                        IProductWRTStdDerivBase(i,deriv[0],Fwd);

                    // map to expansion
                    for(int i = 0; i < nTraceCoeffs; ++i)
                    {
                        for(auto &it:  m_bwdTraceToETraceCoeffMap[i])
                        {
                            FilterCoeffs[it.first] -= Fwd[i]*(it.second);
                        }
                    }
                }
                
                m_dgfield->MultiplyByElmtInvMass(FilterCoeffs,deriv[0]);
                
                m_dgfield->BwdTrans(deriv[0],FilterCoeffs);
                Vmath::Vadd(nphys,outarray[f]+p*nphys,1,FilterCoeffs,1,
                            tmp = outarray[f]+p*nphys,1);
            }
        }
    }

    void ForcingGJPStabilisation::eval_h(LocalRegions::ExpansionSharedPtr elmt,
                                      int traceid, NekDouble &h, NekDouble &p)
    {
        SpatialDomains::GeometrySharedPtr geom = elmt->GetGeom();

        h = 0.0;
        
        switch(geom->GetCoordim())
        {
        case 1:
            {
                h = geom->GetVertex(1)->dist(*geom->GetVertex(0));
                p = elmt->GetNcoeffs(); 
            }
            break;
        case 2:
            {
                int nverts = geom->GetNumVerts();
                int pe;
                //vertices on edges
                SpatialDomains::PointGeom ev0 = *geom->GetVertex(traceid);
                SpatialDomains::PointGeom ev1 = *geom->GetVertex((traceid+1)%
                                                                 nverts);

                //vertex on adjacent edge to ev0 
                SpatialDomains::PointGeom vadj = *geom->GetVertex
                    ((traceid+(nverts-1))%nverts);
                
                // calculate perpendicular distance of normal length
                // from first vertex
                NekDouble h1 = ev0.dist(vadj);
                SpatialDomains::PointGeom Dx, Dx1; 

                Dx.Sub(ev1,ev0);
                Dx1.Sub(vadj,ev0);
                
                NekDouble d1  = Dx.dot(Dx1); 
                NekDouble lenDx = Dx.dot(Dx);
                h = sqrt(h1*h1-d1*d1/lenDx);
                pe = elmt->GetTraceNcoeffs((traceid+nverts-1)%nverts)-1;
                p = pe; 
                    
                // perpendicular distanace from second vertex 
                vadj = *geom->GetVertex((traceid+2)%nverts);

                h1  = ev1.dist(vadj);
                Dx1.Sub(vadj,ev1);
                d1 = Dx.dot(Dx1); 

                h = (h+sqrt(h1*h1-d1*d1/lenDx))*0.5;
                pe = elmt->GetTraceNcoeffs((traceid+1)%nverts)-1;
                p = (p+pe)*0.5;
            }
        break;
        case 3:
        {
            int nverts = geom->GetFace(traceid)->GetNumVerts();

            SpatialDomains::PointGeom tn1,tn2, normal;
            tn1.Sub(*(geom->GetFace(traceid)->GetVertex(1)),
                    *(geom->GetFace(traceid)->GetVertex(0)));
            tn2.Sub(*(geom->GetFace(traceid)->GetVertex(nverts-1)),
                    *(geom->GetFace(traceid)->GetVertex(0)));

            normal.Mult(tn1,tn2);

            //normalise normal
            NekDouble mag = normal.dot(normal);
            mag = 1.0/sqrt(mag); 
            normal.UpdatePosition(normal.x()*mag,
                                  normal.y()*mag,
                                  normal.z()*mag);

            SpatialDomains::PointGeom Dx;
            for(int i = 0; i < nverts; ++i)
            {
                //vertices on edges
                int edgid = geom->GetEdgeNormalToFaceVert(traceid,i); 
                    
                //vector along noramal edge to each vertex 
                Dx.Sub(*(geom->GetEdge(edgid)->GetVertex(0)),
                       *(geom->GetEdge(edgid)->GetVertex(1)));

                // calculate perpendicular distance of normal length
                // from first vertex
                h  += fabs(normal.dot(Dx));
            }
            
            h /= (NekDouble)(nverts);

            // find normal basis direction
            int dir0 = geom->GetDir(traceid,0);
            int dir1 = geom->GetDir(traceid,1);
            int dirn;
            for(dirn = 0; dirn < 3; ++dirn)
            {
                if((dirn != dir0)&&(dirn != dir1))
                {
                    break;
                }
            }
            p = (NekDouble) (elmt->GetBasisNumModes(dirn)-1);
        }
        break;
        default:
        break;
        }
    }
}
}
