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


        funcNameElmt = pForce->FirstChildElement("JUMPSCALING");

        if (funcNameElmt)
        {
            std::string jumpscal = funcNameElmt->GetText();
            m_jumpScal = boost::lexical_cast<NekDouble>(jumpscal);
        }
        else
        {
            m_jumpScal = 1.0;
        }
        
        if(m_session->GetComm()->GetRank() == 0)
        {
            cout << "GJP Stabilisation:" << endl;
            cout << "\t H-Scaling:        " << m_hScalingStr << endl;
            cout << "\t Velocity-Scaling: " << m_velScalingStr << endl;
            cout << "\t jump-Scaling:    " << m_jumpScal << endl;
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
            
        m_dgfield = MemoryManager<MultiRegions::DisContField>::AllocateSharedPtr
            (m_session,pFields[0]->GetGraph(), m_session->GetVariable(0), true, false);

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
        
        SetUpExpansionInfoMapForGJP(pFields[0]->GetGraph());

        MultiRegions::DisContFieldSharedPtr dgfield; 

        dgfield = MemoryManager<MultiRegions::DisContField>::AllocateSharedPtr
            (m_session,pFields[0]->GetGraph(), "GJP", true,
             false, Collections::eNoImpType, m_session->GetVariable(0));
        dgfield->GetLocTraceToTraceMap(m_locTraceToTraceMap); 

        m_locElmtTrace = MemoryManager<MultiRegions::ExpList>::
            AllocateSharedPtr(m_session, *(dgfield->GetExp()),
                              dgfield->GetGraph(),true,"GJP"); 

        //Run a check to see that the trace of the dgfield is the same
        //as the trace of the m_dgfield
        bool sametrace = (m_dgfield->GetTrace()->GetExpSize() == dgfield->GetTrace()->GetExpSize())?
            true: false;

        for(int i = 0; i < m_dgfield->GetTrace()->GetExpSize(); ++i)
        {
            for(int j = 0; j < m_traceDim; ++j)
            {
                if(trace->GetExp(i)->GetBasis(j) !=
                   dgfield->GetTrace()->GetExp(i)->GetBasis(j))
                {
                    cout << "trace base in elmt " << i << " direction " << j <<
                        " has type " <<trace->GetExp(i)->GetBasis(j)->
                        GetPointsType() << " versus " << 
                        dgfield->GetTrace()->GetExp(i)->GetBasis(j)->
                        GetPointsType() << std::endl;
                    cout << "num points is " << trace->GetExp(i)->GetBasis(j)->
                        GetNumPoints() << " versus " << 
                        dgfield->GetTrace()->GetExp(i)->GetBasis(j)->
                        GetNumPoints() << std::endl;
                    sametrace = false;
                }
            }
            if(!sametrace)
            {
                break;
            }
        }
        ASSERTL1(sametrace,"GJP Stabilisation: The trace of the dg field is "
                 "not the same as the trace for the LocTraceToTraceMap. Has "
                 "something changed in the intiialisaiton of the expasions?");

        m_scalTrace = Array<OneD, Array<OneD, NekDouble>>(m_traceDim+1);

        const std::shared_ptr<LocalRegions::ExpansionVector>
            exp = dgfield->GetExp();
        
        Array<OneD, Array<OneD, NekDouble> > dfactors[3];
        Array<OneD, Array<OneD, NekDouble>> LocTrace(m_traceDim+1); 
        Array<OneD, NekDouble> e_tmp;

        for(int i = 0; i < m_traceDim +1; ++i)
        {
            LocTrace[i] = Array<OneD, NekDouble> (m_locTraceToTraceMap->
                                                  GetNLocTracePts());
        }
                
        int cnt              = 0;
        int offset           = 0;
        int offset_phys      = 0;
        int coeff_offset     = 0; 
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >dbasis;
        Array<OneD, Array<OneD, Array<OneD, unsigned int> > >traceToCoeffMap;
        
        Array<OneD, unsigned int> map, map1;
        Array<OneD, int> sign,sign1;
        NekDouble h,p;

        for(int e = 0; e < m_dgfield->GetExpSize(); ++e)
        {
            LocalRegions::ExpansionSharedPtr elmt = (*exp)[e];

            elmt->NormalTraceDerivFactors(dfactors[0],dfactors[1],dfactors[2]);
            
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

                jumpScal *= m_jumpScal;

                //#define GJPDEBUG 1
#if GJPDEBUG
                jumpScal = 1.0;
#endif

                int nptrace = elmt->GetTraceNumPoints(n);
                elmt->GetTraceCoeffMap(n,map);               
                int traceNcoeffs = elmt->GetTraceNcoeffs(n);

                for(int i = 0; i < m_traceDim+1; ++i)
                {
                    Vmath::Smul(nptrace, jumpScal, dfactors[i][n], 1,
                                e_tmp = LocTrace[i] + offset_phys,1);
                }
                
                offset_phys += nptrace;
                offset      += traceNcoeffs; 
            }
            coeff_offset += elmt->GetNcoeffs();
        }


        for(int i = 0; i < m_traceDim+1; ++i)
        {
            m_scalTrace[i] = LocTrace[i];

            if(m_traceDim > 0)
            {
                //multiply by Jacobian and quadrature points. 
                m_locElmtTrace->MultiplyByQuadratureMetric(m_scalTrace[i],m_scalTrace[i]);
            }
        }


        // Assemble list of Matrix Product
        Array<OneD, DNekMatSharedPtr> TraceMat;
        
        int nelmt = 1;
        Array<OneD, const LibUtilities::BasisSharedPtr>
            base_sav = dgfield->GetExp(0)->GetBase();

        dgfield->GetExp(0)->IProductWRTTensorDerivBaseOnTraceMat(TraceMat);

        for(int n = 1; n < dgfield->GetExpSize(); ++n)
        {
            const Array<OneD, const LibUtilities::BasisSharedPtr>&
                base = dgfield->GetExp(n)->GetBase();

            int i; 
            for(i = 0; i < base.size(); ++i)
            {
                if(base[i] != base_sav[i])
                {
                    break;
                }
            }

            if(i == base.size())
            {
                nelmt++;
            }
            else
            {
                // save previous block of data. 
                m_IPWRTDBOnTraceMat.push_back(std::pair<int,Array<OneD, DNekMatSharedPtr>>(nelmt,TraceMat));

                // start new block 
                dgfield->GetExp(n)->IProductWRTTensorDerivBaseOnTraceMat(TraceMat);
                nelmt = 1;
                base_sav = dgfield->GetExp(n)->GetBase();
            }
        }
        // save latest block of data. 
        m_IPWRTDBOnTraceMat.push_back(std::pair<int,Array<OneD, DNekMatSharedPtr>>(nelmt,TraceMat));
    }
    
    void ForcingGJPStabilisation::v_Apply(const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                                          const Array<OneD, Array<OneD, NekDouble> > &inarray,
                                          Array<OneD, Array<OneD, NekDouble> > &outarray,
                                          const NekDouble &time)
    {
        boost::ignore_unused(time,pFields);
        
        int ncoeffs   = m_dgfield->GetNcoeffs();
        int nphys     = m_dgfield->GetNpoints();
        int nTracePts = m_dgfield->GetTrace()->GetTotPoints();
        int nLocETrace = m_locElmtTrace->GetTotPoints();
        int nLocETraceCoeffs = m_locElmtTrace->GetNcoeffs();

        ASSERTL1(nLocETrace == m_scalTrace[0].size(),"expect these to be similar");
        ASSERTL1(nLocETraceCoeffs <= nphys,"storage assumptions presume that nLocETraceCoeffs < nphys");
        
        Array<OneD, Array<OneD, NekDouble> > deriv(3,NullNekDouble1DArray);
        for(int i = 0; i < m_coordDim; ++i)
        {
            deriv[i] = Array<OneD, NekDouble> (nphys); 
        }
        
        int nmax = max(ncoeffs,nphys);
        Array<OneD, NekDouble> FilterCoeffs(nmax);
        Array<OneD, NekDouble> GradJumpOnTrace(nTracePts); 
        Array<OneD, NekDouble> Fwd(nTracePts), Bwd(nTracePts); 
        Array<OneD, NekDouble> unorm(nTracePts,1.0);

        Array<OneD, NekDouble> tmp(nLocETrace);
        Array<OneD, NekDouble> LocElmtTracePhys   = m_locElmtTrace->UpdatePhys();
        Array<OneD, NekDouble> LocElmtTraceCoeffs = m_locElmtTrace->UpdateCoeffs();
        ASSERTL1(LocElmtTracePhys.size() <= nLocETrace,"expect this vector to be at least of size nLocETrace");

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
                Vmath::Vmul(nTracePts,unorm,1,GradJumpOnTrace,1,GradJumpOnTrace,1);

#if GJPDEBUG // debugging
                Vmath::Fill(nTracePts,1.0,GradJumpOnTrace,1); 
#endif
                
                // Interpolate GradJumpOnTrace to Local elemental traces.
                m_locTraceToTraceMap->InterpTraceToLocTrace(0,GradJumpOnTrace, tmp);
                m_locTraceToTraceMap->UnshuffleLocTraces(0,tmp,LocElmtTracePhys);
                m_locTraceToTraceMap->InterpTraceToLocTrace(1,GradJumpOnTrace, tmp);
                m_locTraceToTraceMap->UnshuffleLocTraces(1,tmp,LocElmtTracePhys);

                // Scale jump on trace
                Vmath::Vmul(nLocETrace,m_scalTrace[0],1,LocElmtTracePhys,1,tmp,1);
                MultiplyByIProductWRTDerivOnTraceMat(0,tmp,FilterCoeffs);
                
                for(int i = 0; i < m_traceDim; ++i)
                {
                    // Scale jump on trace
                    Vmath::Vmul(nLocETrace,m_scalTrace[i+1],1,LocElmtTracePhys,1,tmp,1);
                    MultiplyByIProductWRTDerivOnTraceMat(i+1,tmp,deriv[0]);
                    Vmath::Vadd(ncoeffs,deriv[0],1,FilterCoeffs,1,FilterCoeffs,1);
                }
                Vmath::Neg(ncoeffs,FilterCoeffs,1);
                m_dgfield->MultiplyByElmtInvMass(FilterCoeffs,deriv[0]);

#if GJPDEBUG
                {
                    Array<OneD, NekDouble> FilterCoeffs1(ncoeffs,0.0);
                    Array<OneD, NekDouble> GradFwdOnTrace(nTracePts);
                    Array<OneD, NekDouble> GradBwdOnTrace(nTracePts);
                    Array<OneD, NekDouble> tmp1(nphys);
                    
                    for(int i = 0; i < ncoeffs; ++i)
                    {
                        Vmath::Zero(ncoeffs,deriv[0],1);
                        deriv[0][i] = 1.0;
                        m_dgfield->BwdTrans(deriv[0],tmp1);

                        // calculate derivative 
                        m_dgfield->PhysDeriv(tmp1,deriv[0],
                                             deriv[1], deriv[2]);

                        Vmath::Zero(nTracePts,GradFwdOnTrace,1);
                        Vmath::Zero(nTracePts,GradBwdOnTrace,1);
                        for(int n = 0; n < m_coordDim; ++n)
                        {
                            // can't use with Periodic BCs
                            m_dgfield->GetFwdBwdTracePhys(deriv[n],Fwd,Bwd,
                                                          false,false,false);
                 
                            //Vmath::Vsub(nTracePts,Fwd,1,Bwd,1,Fwd,1);
                            Vmath::Vvtvp(nTracePts,Fwd,1,m_traceNormals[n],1,
                                         GradFwdOnTrace,1,GradFwdOnTrace,1);
                            Vmath::Vvtvp(nTracePts,Bwd,1,m_traceNormals[n],1,
                                         GradBwdOnTrace,1,GradBwdOnTrace,1);
                        }
                        Vmath::Neg(nTracePts,GradBwdOnTrace,1);
                        
                        Vmath::Zero(nLocETrace,LocElmtTracePhys,1);

                        // Interpolate GradJumpOnTrace to Local elemental traces.
                        Vmath::Zero(nLocETrace,tmp,1);
                        m_locTraceToTraceMap->InterpTraceToLocTrace
                                                   (0,GradFwdOnTrace, tmp);
                        m_locTraceToTraceMap->UnshuffleLocTraces
                                                   (0,tmp,LocElmtTracePhys);
                        Vmath::Zero(nLocETrace,tmp,1);
                        m_locTraceToTraceMap->InterpTraceToLocTrace
                                                   (1,GradBwdOnTrace, tmp);
                        m_locTraceToTraceMap->UnshuffleLocTraces
                                                   (1,tmp,LocElmtTracePhys);
                         
                        // Integrate 
                        FilterCoeffs1[i] = m_locElmtTrace->PhysIntegral(LocElmtTracePhys); 
                    }

                    for(int i = 0; i < ncoeffs; ++i)
                    {
                        if(fabs(FilterCoeffs1[i]+FilterCoeffs[i]) > 1e-10)
                        {
                            cout << "i= " << i << " diff is " <<
                                fabs(FilterCoeffs1[i]+FilterCoeffs[i])
                                 << " Direct Coeff: " << FilterCoeffs1[i]
                                 << " Indirect Coeff: " << FilterCoeffs[i]
                                 << std::endl;
                            
                        }
                    }
                    cout << "Done testing " << ncoeffs << " coefficients" <<
                        std::endl;
                    exit(1);
                }
#endif
                
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


    void ForcingGJPStabilisation::SetUpExpansionInfoMapForGJP(SpatialDomains::MeshGraphSharedPtr graph)
    {
        const SpatialDomains::ExpansionInfoMap  expInfo = graph->GetExpansionInfo(m_session->GetVariable(0));

        SpatialDomains::ExpansionInfoMapShPtr newInfo = MemoryManager<SpatialDomains::ExpansionInfoMap>::AllocateSharedPtr();

        // loop over epxansion info
        for (auto expIt = expInfo.begin(); expIt != expInfo.end(); ++expIt)
        {
            std::vector<LibUtilities::BasisKey> BKeyVector;
            

            for (int i = 0; i < expIt->second->m_basisKeyVector.size(); ++i)
            {
                LibUtilities::BasisKey bkeyold =
                    expIt->second->m_basisKeyVector[i];

                // reset radauM alpha non-zero cases to radauM Lenendre at one order higher
                switch(bkeyold.GetPointsType())
                {
                case LibUtilities::eGaussRadauMAlpha1Beta0:
                case LibUtilities::eGaussRadauMAlpha2Beta0:
                    {
                        int npts = bkeyold.GetNumPoints();
                        
                        //const LibUtilities::PointsKey pkey(npts+1, LibUtilities::eGaussRadauMLegendre);
                        // trying npts to be consistent for tri faces 
                        const LibUtilities::PointsKey pkey(npts, LibUtilities::eGaussRadauMLegendre);
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
            
            (*newInfo)[expIt->first] = MemoryManager<SpatialDomains::ExpansionInfo>::
                AllocateSharedPtr(expIt->second->m_geomShPtr, BKeyVector);
        }
                

        graph->SetExpansionInfo("GJP",newInfo);
    }

    void ForcingGJPStabilisation::MultiplyByIProductWRTDerivOnTraceMat(int i, Array<OneD, NekDouble> &in,
                                                                       Array<OneD, NekDouble> &out)
    {
        // Should be vectorised

        int cnt = 0; 
        int cnt1 = 0; 
        for(auto &it: m_IPWRTDBOnTraceMat)
        {
            int rows = it.second[i]->GetRows();
            int cols = it.second[i]->GetColumns();
            
            Blas::Dgemm('N','N', rows, it.first, cols,  1.0, &(it.second[i]->GetPtr())[0],
                        rows, &in[0]+cnt,  cols, 
                        0.0, &out[0]+cnt1, rows);

            cnt  += cols*it.first;
            cnt1 += rows*it.first; 
        }
        
    }
}
}
