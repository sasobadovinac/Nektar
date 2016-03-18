///////////////////////////////////////////////////////////////////////////////
//
// File: Extrapolate.cpp
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
// Description: Abstract base class for Extrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/Extrapolate.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{
    NekDouble Extrapolate::StifflyStable_Betaq_Coeffs[3][3] = {
        { 1.0,  0.0, 0.0},{ 2.0, -1.0, 0.0},{ 3.0, -3.0, 1.0}};
    NekDouble Extrapolate::StifflyStable_Alpha_Coeffs[3][3] = {
        { 1.0,  0.0, 0.0},{ 2.0, -0.5, 0.0},{ 3.0, -1.5, 1.0/3.0}};
    NekDouble Extrapolate::StifflyStable_Gamma0_Coeffs[3] = {
          1.0,  1.5, 11.0/6.0};

    ExtrapolateFactory& GetExtrapolateFactory()
    {
        typedef Loki::SingletonHolder<ExtrapolateFactory,
                                      Loki::CreateUsingNew,
                                      Loki::NoDestroy,
                                      Loki::SingleThreaded > Type;
        return Type::Instance();
    }

    Extrapolate::Extrapolate(
        const LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
        MultiRegions::ExpListSharedPtr pPressure,
        const Array<OneD, int> pVel,
        const SolverUtils::AdvectionSharedPtr advObject)
        : m_session(pSession),
          m_fields(pFields),
          m_pressure(pPressure),
          m_velocity(pVel),
          m_advObject(advObject)
    {      
        m_session->LoadParameter("TimeStep", m_timestep,   0.01);
        m_comm = m_session->GetComm();
    }
    
    Extrapolate::~Extrapolate()
    {
    }

        
    std::string Extrapolate::def =
        LibUtilities::SessionReader::RegisterDefaultSolverInfo(
            "StandardExtrapolate", "StandardExtrapolate");
    
    /** 
     * Function to extrapolate the new pressure boundary condition.
     * Based on the velocity field and on the advection term.
     * Acceleration term is also computed.
     * This routine is a general one for 2d and 3D application and it can be called
     * directly from velocity correction scheme. Specialisation on dimensionality is
     * redirected to the CalcPressureBCs method.
     */
    void Extrapolate::EvaluatePressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        const Array<OneD, const Array<OneD, NekDouble> >  &N,
        NekDouble kinvis)
    {
        m_pressureCalls++;
        if(m_HBCnumber > 0)
        {
            int  n,cnt;

            // Calculate Neumann BCs at current level
            CalcNeumannPressureBCs(fields, N, kinvis);

            //Calculate acceleration term at level n based on previous steps
            AccelerationBDF(m_acceleration);

            // Adding acceleration term to HOPBCs
            Vmath::Svtvp(m_numHBCDof, -1.0/m_timestep,
                         m_acceleration[m_intSteps],  1,
                         m_pressureHBCs[m_intSteps-1], 1,
                         m_pressureHBCs[m_intSteps-1], 1);

            // Extrapolate to n+1
            ExtrapolateArray(m_pressureHBCs);

            // Copy values of [dP/dn]^{n+1} in the pressure bcs storage.
            // m_pressureHBCS[nlevels-1] will be cancelled at next time step
            for(cnt = n = 0; n < m_PBndConds.num_elements(); ++n)
            {
                if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"H"))
                {
                    int nq = m_PBndExp[n]->GetNcoeffs();
                    Vmath::Vcopy(nq, &(m_pressureHBCs[m_intSteps-1])[cnt],  1,
                                     &(m_PBndExp[n]->UpdateCoeffs()[0]), 1);
                    cnt += nq;
                }
            }

        }

        CalcOutflowBCs(fields, N, kinvis);
    }


    /**
     * Unified routine for calculation high-oder terms
     */
    void Extrapolate::v_CalcNeumannPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        const Array<OneD, const Array<OneD, NekDouble> >  &N,
        NekDouble kinvis)
    {
        int n, cnt;

        Array<OneD, NekDouble> Pvals;
        Array<OneD, NekDouble> Uvals;

        Array<OneD, Array<OneD, NekDouble> > Velocity(m_curl_dim);
        Array<OneD, Array<OneD, NekDouble> > Advection(m_bnd_dim);

        Array<OneD, Array<OneD, NekDouble> > BndValues(m_bnd_dim);
        Array<OneD, Array<OneD, NekDouble> > Q(m_curl_dim);

        MultiRegions::ExpListSharedPtr BndElmtExp;
        for(n = cnt = 0; n < m_PBndConds.num_elements(); ++n)
        {
            // High order boundary condition;
            if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"H"))
            {
                m_fields[0]->GetBndElmtExpansion(n, BndElmtExp);
                int nqb = m_PBndExp[n]->GetTotPoints();
                int nq  = BndElmtExp->GetTotPoints();

                for(int i = 0; i < m_bnd_dim; i++)
                {
                    BndValues[i] = Array<OneD, NekDouble> (nqb,0.0);
                }

                for(int i = 0; i < m_curl_dim; i++)
                {
                    Q[i]         = Array<OneD, NekDouble> (nq,0.0);
                }

                // Obtaining fields on BndElmtExp
                for(int i = 0; i < m_curl_dim; i++)
                {
                    m_fields[0]->ExtractPhysToBndElmt(n, fields[i],Velocity[i]);
                }
                for(int i = 0; i < m_bnd_dim; i++)
                {
                    m_fields[0]->ExtractPhysToBndElmt(n, N[i],Advection[i]);
                }

                // CurlCurl
                BndElmtExp->CurlCurl(Velocity, Q);

                // Mounting advection component into the high-order condition
                for(int i = 0; i < m_bnd_dim; i++)
                {
                    MountHOPBCs(nq, kinvis,Q[i],Advection[i]);
                }

                Pvals = (m_pressureHBCs[m_intSteps-1]) + cnt;
                Uvals = (m_acceleration[m_intSteps]) + cnt;

                // Getting values on the edge and filling the pressure boundary
                // expansion and the acceleration term. Multiplication by the
                // normal is required
                for(int i = 0; i < m_bnd_dim; i++)
                {
                    m_fields[0]->ExtractElmtToBndPhys(n, Q[i],BndValues[i]);
                }
                m_PBndExp[n]->NormVectorIProductWRTBase(BndValues, Pvals);

                for(int i = 0; i < m_bnd_dim; i++)
                {
                    m_fields[0]->ExtractElmtToBndPhys(n, Velocity[i],BndValues[i]);
                }
                m_PBndExp[n]->NormVectorIProductWRTBase(BndValues, Uvals);

                // Get offset for next terms
                cnt += m_PBndExp[n]->GetNcoeffs();
            }
        }
    }



    void Extrapolate::CalcOutflowBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        const Array<OneD, const Array<OneD, NekDouble> >  &N,
        NekDouble kinvis)
    {
        static bool init = true;
        static bool noHOBC = false;

        if(noHOBC == true)
        {
           return;
        }
        
        if(init) // set up storage for boundary velocity at outflow
        {
            init = false;
            int totbndpts = 0;
            for(int n = 0; n < m_PBndConds.num_elements(); ++n)
            {
                if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"HOutflow"))
                {
                    totbndpts += m_PBndExp[n]->GetTotPoints();
                }
            }
          
            if(totbndpts == 0)
            {
                noHOBC = true;
                return;
            }
  
            m_outflowVel = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (m_curl_dim);
            for(int i = 0; i < m_curl_dim; ++i)
            {
                m_outflowVel[i] = Array<OneD, Array<OneD, NekDouble> >(m_curl_dim);
                for(int j = 0; j < m_curl_dim; ++j)
                {             
                    // currently just set up for 2nd order extrapolation 
                    m_outflowVel[i][j] = Array<OneD, NekDouble>(totbndpts,0.0);
                }
            }

            if (m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
            {
                m_PhyoutfVel = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (m_curl_dim);

                for(int i = 0; i < m_curl_dim; ++i)
                {
                    m_PhyoutfVel[i] = Array<OneD, Array<OneD, NekDouble> > (m_curl_dim);
                    for(int j = 0; j < m_curl_dim; ++j)
                    {
                        // currently just set up for 2nd order extrapolation
                        m_PhyoutfVel[i][j] = Array<OneD, NekDouble> (totbndpts,0.0);
                    }
                }

                m_nonlinearterm_phys   = Array<OneD, NekDouble> (totbndpts,0.0);
                m_nonlinearterm_coeffs = Array<OneD, NekDouble> (totbndpts,0.0);

                m_PBndCoeffs = Array<OneD, NekDouble> (totbndpts,0.0);
                m_UBndCoeffs = Array<OneD, Array<OneD, NekDouble> > (m_curl_dim);
                for(int i = 0; i < m_curl_dim; ++i)
                {
                    m_UBndCoeffs[i] = Array<OneD, NekDouble> (totbndpts);   
                }
                Array<OneD, unsigned int> planes;
                planes = m_pressure->GetZIDs();
                int num_planes = planes.num_elements();
                m_expsize_per_plane = Array<OneD, unsigned int> (m_PBndConds.num_elements());
                for(int n = 0; n < m_PBndConds.num_elements(); ++n)
                {
                    int exp_size = m_PBndExp[n]->GetExpSize();
                    m_expsize_per_plane[n] = exp_size/num_planes;
                }
                m_totexps_per_plane = 0;
                for(int n = 0; n < m_PBndConds.num_elements(); ++n)
                {
                    m_totexps_per_plane += m_PBndExp[n]->GetExpSize()/num_planes;
                }
            }
        }
        
        StdRegions::StdExpansionSharedPtr Bc,Pbc; 
        Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr> >
                                                        UBndConds(m_curl_dim);
        Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >
                                                        UBndExp(m_curl_dim);

        for (int i = 0; i < m_curl_dim; ++i)
        {
            UBndConds[i] = m_fields[m_velocity[i]]->GetBndConditions();
            UBndExp[i]   = m_fields[m_velocity[i]]->GetBndCondExpansions();
        }

        Array<OneD, Array<OneD, NekDouble> > BndValues(m_curl_dim);
        Array<OneD, Array<OneD, NekDouble> > BndElmt  (m_curl_dim);
        Array<OneD, Array<OneD, NekDouble> > nGradu   (m_curl_dim);
        Array<OneD, NekDouble > gradtmp (m_pressureBCsElmtMaxPts),
                                fgradtmp(m_pressureBCsElmtMaxPts);
        
        nGradu[0] = Array<OneD, NekDouble>(m_curl_dim*m_pressureBCsMaxPts);
        for(int i = 0; i < m_curl_dim; ++i)
        {
            BndElmt[i]   = Array<OneD, NekDouble> (m_pressureBCsElmtMaxPts,
                                                   0.0);
            BndValues[i] = Array<OneD, NekDouble> (m_pressureBCsMaxPts,0.0);
            nGradu[i]  = nGradu[0] + i*m_pressureBCsMaxPts;
            RollOver(m_outflowVel[i]);

            if (m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
            {
                RollOver(m_PhyoutfVel[i]);
            }
        }
            
        int nbc,cnt,cnt_start;
        int veloffset = 0;
        int  nint    = min(m_pressureCalls,m_intSteps);

        StdRegions::StdExpansionSharedPtr elmt;
        Array<OneD, NekDouble> PBCvals, UBCvals;
        Array<OneD, Array<OneD, NekDouble> > ubc(m_curl_dim);
        Array<OneD, Array<OneD, NekDouble> > normals;

        cnt = 0;
        for(int n = 0; n < m_PBndConds.num_elements(); ++n)
        {
            // Do outflow boundary conditions if they exist
            if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"HOutflow"))
            {
                for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i,cnt++)
                {
                    cnt = max(cnt,m_PBndExp[n]->GetTotPoints());
                }
            }
        }

        for(int i =0; i < m_curl_dim; ++i)
        {
            ubc[i] = Array<OneD, NekDouble>(cnt);
        }

        NekDouble U0,delta;
        m_session->LoadParameter("U0_HighOrderBC",U0,1.0);
        m_session->LoadParameter("Delta_HighOrderBC",delta,1/20.0);

        cnt = 0;
        for(int n = 0; n < m_PBndConds.num_elements(); ++n)
        {
            cnt_start = cnt;

            // Do outflow boundary conditions if they exist
            if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"HOutflow"))
            {

                if (m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
                {
                    int veloffset = 0;
                    for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i, cnt++)
                    {
                        // find element and edge of this expansion.
                        Bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion>
                            (m_PBndExp[n]->GetExp(i));

                        int elmtid = m_pressureBCtoElmtID[cnt];
                        elmt       = m_fields[0]->GetExp(elmtid);
                        int offset = m_fields[0]->GetPhys_Offset(elmtid);

                        int boundary = m_pressureBCtoTraceID[cnt];

                        // Determine extrapolated U,V values
                        int nq = elmt->GetTotPoints();
                        int nbc = m_PBndExp[n]->GetExp(i)->GetTotPoints();
                        // currently just using first order approximation here.
                        // previously have obtained value from m_integrationSoln
                        Array<OneD, NekDouble> veltmp;

                        for(int j = 0; j < m_curl_dim; ++j)
                        {
                            Vmath::Vcopy(nq, &fields[m_velocity[j]][offset], 1,
                                         &BndElmt[j][0],                 1);
                            elmt->GetTracePhysVals(boundary,Bc,BndElmt[j],
                                         veltmp = m_outflowVel[j][0] + veloffset);
                        }
                        veloffset += nbc;
                    }

                    // for velocity on the outflow boundary in e3DH1D,
                    // we need to make a backward fourier transformation
                    // to get the physical coeffs at the outflow BCs.
                    for(int j = 0; j < m_curl_dim; ++j)
                    {
                        m_PBndExp[n]->HomogeneousBwdTrans(
                                m_outflowVel[j][0],
                                m_PhyoutfVel[j][0]);
                    }

                    cnt = cnt_start;
                    veloffset = 0;
                    for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i, cnt++)
                    {

                        int elmtid = m_pressureBCtoElmtID[cnt];
                        elmt       = m_fields[0]->GetExp(elmtid);
                        int nbc = m_PBndExp[n]->GetExp(i)->GetTotPoints();

                        Array<OneD, NekDouble>  veltmp(nbc,0.0),
                                                normDotu(nbc,0.0), utot(nbc,0.0);
                        int boundary = m_pressureBCtoTraceID[cnt];
                        normals=elmt->GetSurfaceNormal(boundary);

                        // extrapolate velocity
                        if(nint <= 1)
                        {
                            for(int j = 0; j < m_curl_dim; ++j)
                            {
                                Vmath::Vcopy(nbc,
                                        veltmp = m_PhyoutfVel[j][0] +veloffset, 1,
                                        BndValues[j],                           1);
                            }
                        }
                        else // only set up for 2nd order extrapolation
                        {
                            for(int j = 0; j < m_curl_dim; ++j)
                            {
                                Vmath::Smul(nbc, 2.0,
                                        veltmp = m_PhyoutfVel[j][0] + veloffset, 1,
                                        BndValues[j],                            1);
                                Vmath::Svtvp(nbc, -1.0,
                                        veltmp = m_PhyoutfVel[j][1] + veloffset, 1,
                                        BndValues[j],                            1,
                                        BndValues[j],                            1);
                            }    
                        }

                        // Set up |u|^2, n.u in physical space
                        for(int j = 0; j < m_curl_dim; ++j)
                        {
                            Vmath::Vvtvp(nbc, BndValues[j], 1, BndValues[j], 1,
                                              utot,         1, utot,         1);
                        }
                        for(int j = 0; j < m_bnd_dim; ++j)
                        {
                            Vmath::Vvtvp(nbc, normals[j],   1, BndValues[j], 1,
                                              normDotu,     1, normDotu,     1);
                        }
                        
                        int Offset = m_PBndExp[n]->GetPhys_Offset(i);

                        for(int k = 0; k < nbc; ++k)
                        {
                            // calculate the nonlinear term (kinetic energy
                            //  multiplies step function) in physical space
                            NekDouble fac = 0.5*(1.0-tanh(normDotu[k]/(U0*delta)));
                            m_nonlinearterm_phys[k + Offset] =  0.5 * utot[k] * fac;
                        }

                        veloffset += nbc;
                    }

                    // for e3DH1D, we need to make a forward fourier transformation
                    // for the kinetic energy term (nonlinear)
                    UBndExp[0][n]->HomogeneousFwdTrans(
                            m_nonlinearterm_phys,
                            m_nonlinearterm_coeffs);

                    // for e3DH1D, we need to make a forward fourier transformation
                    // for Dirichlet pressure boundary condition that is from input file
                    m_PBndExp[n]->HomogeneousFwdTrans(
                            m_PBndExp[n]->UpdatePhys(),
                            m_PBndCoeffs);
                    // for e3DH1D, we need to make a forward fourier transformation
                    // for Neumann velocity boundary condition that is from input file
                    for (int j = 0; j < m_curl_dim; ++j)
                    {
                        UBndExp[j][n]->HomogeneousFwdTrans(
                            UBndExp[j][n]->UpdatePhys(),
                            m_UBndCoeffs[j]);
                    }
                }

                cnt = cnt_start;
                veloffset = 0;
                for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i,cnt++)
                {
                    // find element and edge of this expansion. 
                    Bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion> 
                        (m_PBndExp[n]->GetExp(i));

                    int elmtid = m_pressureBCtoElmtID[cnt];
                    elmt       = m_fields[0]->GetExp(elmtid);
                    int offset = m_fields[0]->GetPhys_Offset(elmtid);

                    // Determine extrapolated U,V values
                    int nq = elmt->GetTotPoints();

                    // currently just using first order approximation here. 
                    // previously have obtained value from m_integrationSoln
                    for(int j = 0; j < m_bnd_dim; ++j)
                    {
                        Vmath::Vcopy(nq, &fields[m_velocity[j]][offset], 1,
                                         &BndElmt[j][0],                 1);
                    }

                    int nbc      = m_PBndExp[n]->GetExp(i)->GetTotPoints();
                    int boundary = m_pressureBCtoTraceID[cnt];

                    Array<OneD, NekDouble>  ptmp(nbc,0.0),
                                            divU(nbc,0.0);

                    normals=elmt->GetSurfaceNormal(boundary);
                    Vmath::Zero(m_bnd_dim*m_pressureBCsMaxPts,nGradu[0],1);

                    for (int j = 0; j < m_bnd_dim; j++)
                    {
                        // Calculate Grad u =  du/dx, du/dy, du/dz, etc. 
                        for (int k = 0; k< m_bnd_dim; k++)
                        {
                            elmt->PhysDeriv(MultiRegions::DirCartesianMap[k],
                                            BndElmt[j], gradtmp);
                            elmt->GetTracePhysVals(boundary, Bc, gradtmp,
                                                   fgradtmp);
                            Vmath::Vvtvp(nbc,normals[k], 1, fgradtmp, 1,
                                             nGradu[j],  1, nGradu[j],1);
                            if(j == k)
                            {
                                Vmath::Vadd(nbc,fgradtmp, 1, divU, 1, divU, 1);
                            }
                        }
                    }

                    if (m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
                    {
                        // Set up |u|^2, n.u, div(u), and (n.grad(u) . n) for
                        // pressure condition
                        for(int j = 0; j < m_bnd_dim; ++j)
                        {
                            Vmath::Vvtvp(nbc, normals[j],   1, nGradu[j],    1,
                                              ptmp,         1, ptmp,         1);
                        }
                        int p_offset = m_PBndExp[n]->GetPhys_Offset(i);

                        for(int k = 0; k < nbc; ++k)
                        {
                            // Set up Dirichlet pressure condition and
                            // store in ptmp (m_UBndCoeffs contains Fourier Coeffs of the
                            // function from the input file )

                            ptmp[k] =  kinvis * ptmp[k] 
                                        - m_nonlinearterm_coeffs[k + p_offset]
                                                  - m_PBndCoeffs[k + p_offset];
                        }

                        int u_offset = UBndExp[0][n]->GetPhys_Offset(i);

                        for(int j = 0; j < m_bnd_dim; ++j)
                        {
                            for(int k = 0; k < nbc; ++k)
                            {
                                ubc[j][k + u_offset] = (1.0 / kinvis)
                                                * (m_UBndCoeffs[j][k + u_offset]
                                          + m_nonlinearterm_coeffs[k + u_offset]
                                                                * normals[j][k]);
                            }
                        }

                        // boundary condition for velocity in homogenous direction
                        for(int k = 0; k < nbc; ++k)
                        {
                            ubc[m_bnd_dim][k + u_offset] = (1.0 / kinvis)
                                                * m_UBndCoeffs[m_bnd_dim][k + u_offset];
                        }

                        u_offset = UBndExp[m_bnd_dim][n]->GetPhys_Offset(i);
                        UBCvals  = UBndExp[m_bnd_dim][n]->UpdateCoeffs()
                                    + UBndExp[m_bnd_dim][n]->GetCoeff_Offset(i);
                        Bc->IProductWRTBase(ubc[m_bnd_dim] + u_offset, UBCvals);
                    }
                    else
                    {

                        Array<OneD, NekDouble>  veltmp, utot(nbc,0.0),
                                                normDotu(nbc,0.0);
                        // extract velocity and store
                        for(int j = 0; j < m_bnd_dim; ++j)
                        {
                            elmt->GetTracePhysVals(boundary,Bc,BndElmt[j],
                                           veltmp = m_outflowVel[j][0] + veloffset);
                        }

                        // extrapolate velocity
                        if(nint <= 1)
                        {
                            for(int j = 0; j < m_bnd_dim; ++j)
                            {
                                Vmath::Vcopy(nbc,
                                        veltmp = m_outflowVel[j][0] +veloffset, 1,
                                        BndValues[j],                           1);
                            }
                        }
                        else // only set up for 2nd order extrapolation
                        {
                            for(int j = 0; j < m_bnd_dim; ++j)
                            {
                                Vmath::Smul(nbc, 2.0,
                                        veltmp = m_outflowVel[j][0] + veloffset, 1,
                                        BndValues[j],                            1);
                                Vmath::Svtvp(nbc, -1.0,
                                        veltmp = m_outflowVel[j][1] + veloffset, 1,
                                        BndValues[j],                            1,
                                        BndValues[j],                            1);
                            }
                        }

                        // Set up |u|^2, n.u, div(u), and (n.grad(u) . n) for
                        // pressure condition
                        for(int j = 0; j < m_bnd_dim; ++j)
                        {
                            Vmath::Vvtvp(nbc, BndValues[j], 1, BndValues[j], 1,
                                              utot,         1, utot,         1);
                            Vmath::Vvtvp(nbc, normals[j],   1, BndValues[j], 1,
                                              normDotu,     1, normDotu,     1);
                            Vmath::Vvtvp(nbc, normals[j],   1, nGradu[j],    1,
                                              ptmp,         1, ptmp,         1);
                        }

                        PBCvals = m_PBndExp[n]->GetPhys() +
                                        m_PBndExp[n]->GetPhys_Offset(i);

                        for(int k = 0; k < nbc; ++k)
                        {
                            NekDouble fac = 0.5*(1.0-tanh(normDotu[k]/(U0*delta)));

                            // Set up Dirichlet pressure condition and
                            // store in ptmp (PBCvals contains a
                            // function from the input file )
                            ptmp[k] =  kinvis * ptmp[k] - 0.5 * utot[k] * fac
                                                        - PBCvals[k];
                        }

                        int u_offset = UBndExp[0][n]->GetPhys_Offset(i);

                        for(int j = 0; j < m_bnd_dim; ++j)
                        {
                            UBCvals = UBndExp[j][n]->GetPhys()
                                        + UBndExp[j][n]->GetPhys_Offset(i);

                            for(int k = 0; k < nbc; ++k)
                            {
                                NekDouble fac        = 0.5 * (1.0 - tanh(normDotu[k]
                                                                / (U0 * delta)));
                                ubc[j][k + u_offset] = (1.0 / kinvis)
                                                * (UBCvals[k] + 0.5 * utot[k] * fac
                                                                * normals[j][k]);
                            }
                        }
                    }

                    // set up pressure boundary condition
                    PBCvals = m_PBndExp[n]->UpdateCoeffs()
                            + m_PBndExp[n]->GetCoeff_Offset(i);
                    m_PBndExp[n]->GetExp(i)->FwdTrans(ptmp,PBCvals);

                    veloffset += nbc;
                }

                // Now set up Velocity conditions.
                for(int j = 0; j < m_bnd_dim; j++)
                {
                    if(boost::iequals(UBndConds[j][n]->GetUserDefined(),"HOutflow"))
                    {
                        cnt = cnt_start;
                        for(int i = 0; i < UBndExp[0][n]->GetExpSize();
                                       ++i, cnt++)
                        {
                            Pbc =  StdRegions::StdExpansionSharedPtr
                                            (m_PBndExp[n]->GetExp(i));
                            Bc  =  StdRegions::StdExpansionSharedPtr
                                            (UBndExp[0][n]->GetExp(i));

                            nbc = UBndExp[0][n]->GetExp(i)->GetTotPoints();
                            int boundary = m_pressureBCtoTraceID[cnt];

                            Array<OneD, NekDouble> pb(nbc), ub(nbc);
                            int elmtid = m_pressureBCtoElmtID[cnt];

                            elmt   = m_fields[0]->GetExp(elmtid);

                            normals = elmt->GetSurfaceNormal(boundary);

                            // Get p from projected boundary condition
                            PBCvals = m_PBndExp[n]->UpdateCoeffs()
                                    + m_PBndExp[n]->GetCoeff_Offset(i);
                            Pbc->BwdTrans(PBCvals,pb);

                            int u_offset = UBndExp[j][n]->GetPhys_Offset(i);

                            for(int k = 0; k < nbc; ++k)
                            {
                                ub[k] = ubc[j][k + u_offset]
                                      + pb[k] * normals[j][k] / kinvis;
                            }

                            UBCvals = UBndExp[j][n]->UpdateCoeffs()
                                    + UBndExp[j][n]->GetCoeff_Offset(i);
                            Bc->IProductWRTBase(ub,UBCvals);
                        }
                    }
                }
            }
            else
            {
                cnt += m_PBndExp[n]->GetExpSize();
            }
        }
    }

    /** 
     * Function to roll time-level storages to the next step layout.
     * The stored data associated with the oldest time-level 
     * (not required anymore) are moved to the top, where they will
     * be overwritten as the solution process progresses.
     */
    void Extrapolate::RollOver(Array<OneD, Array<OneD, NekDouble> > &input)
    {
        int  nlevels = input.num_elements();
        
        Array<OneD, NekDouble> tmp;
    
        tmp = input[nlevels-1];
    
        for(int n = nlevels-1; n > 0; --n)
        {
            input[n] = input[n-1];
        }
    
        input[0] = tmp;
    }
    
    
    /**
     * Map to directly locate HOPBCs position and offsets in all scenarios
     */
    void Extrapolate::GenerateHOPBCMap()
    {
        m_PBndConds   = m_pressure->GetBndConditions();
        m_PBndExp     = m_pressure->GetBndCondExpansions();
    
        // Set up mapping from pressure boundary condition to pressure element
        // details.
        m_pressure->GetBoundaryToElmtMap(m_pressureBCtoElmtID,
                                         m_pressureBCtoTraceID);

        // find the maximum values of points  for pressure BC evaluation
        m_pressureBCsMaxPts = 0; 
        m_pressureBCsElmtMaxPts = 0; 
        int cnt, n;
        for(cnt = n = 0; n < m_PBndConds.num_elements(); ++n)
        {
            for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i)
            {
                m_pressureBCsMaxPts = max(m_pressureBCsMaxPts,
                                m_PBndExp[n]->GetExp(i)->GetTotPoints());
                m_pressureBCsElmtMaxPts = max(m_pressureBCsElmtMaxPts,
                                m_pressure->GetExp(m_pressureBCtoElmtID[cnt++])
                                                            ->GetTotPoints());
            }
        }
    
        // Storage array for high order pressure BCs
        m_pressureHBCs = Array<OneD, Array<OneD, NekDouble> > (m_intSteps);
        m_acceleration = Array<OneD, Array<OneD, NekDouble> > (m_intSteps + 1);
    
        m_HBCnumber = 0;
        m_numHBCDof = 0;
        for( n = 0; n < m_PBndConds.num_elements(); ++n)
        {
            // High order boundary condition;
            if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"H"))
            {
                m_numHBCDof += m_PBndExp[n]->GetNcoeffs();
                m_HBCnumber += m_PBndExp[n]->GetExpSize();
            }
        }

        //int checkHBC = m_HBCnumber;
        //m_comm->AllReduce(checkHBC,LibUtilities::ReduceSum);
        //ASSERTL0(checkHBC > 0 ,"At least one high-order pressure boundary "
        //                       "condition is required for scheme "
        //                       "consistency");

        m_acceleration[0] = Array<OneD, NekDouble>(m_numHBCDof, 0.0);
        for(n = 0; n < m_intSteps; ++n)
        {
            m_pressureHBCs[n]   = Array<OneD, NekDouble>(m_numHBCDof, 0.0);
            m_acceleration[n+1] = Array<OneD, NekDouble>(m_numHBCDof, 0.0);
        }

        m_pressureCalls = 0;
        
        switch(m_pressure->GetExpType())
        {
            case MultiRegions::e2D:
            {
                m_curl_dim = 2;
                m_bnd_dim  = 2;
            }
            break;
            case MultiRegions::e3DH1D:
            {
                m_curl_dim = 3;
                m_bnd_dim  = 2;
            }
            break;
            case MultiRegions::e3DH2D:
            {
                m_curl_dim = 3;
                m_bnd_dim  = 1;
            }
            break;
            case MultiRegions::e3D:
            {
                m_curl_dim = 3;
                m_bnd_dim  = 3;
            }
            break;
            default:
                ASSERTL0(0,"Dimension not supported");
                break;
        }
    }

    /**
     *
     */
    Array<OneD, NekDouble> Extrapolate::GetMaxStdVelocity(
        const Array<OneD, Array<OneD,NekDouble> > inarray)
    {
        // Checking if the problem is 2D
        ASSERTL0(m_curl_dim >= 2, "Method not implemented for 1D");
        
        int n_points_0      = m_fields[0]->GetExp(0)->GetTotPoints();
        int n_element       = m_fields[0]->GetExpSize();       
        int nvel            = inarray.num_elements();
        int cnt; 
        
        NekDouble pntVelocity;
        
        // Getting the standard velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(nvel);
        Array<OneD, NekDouble> maxV(n_element, 0.0);
        LibUtilities::PointsKeyVector ptsKeys;
        
        for (int i = 0; i < nvel; ++i)
        {
            stdVelocity[i] = Array<OneD, NekDouble>(n_points_0);
        }
        
        if (nvel == 2)
        {
            cnt = 0.0;
            for (int el = 0; el < n_element; ++el)
            { 
                int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();
                
                // reset local space if necessary
                if(n_points != n_points_0)
                {
                    for (int j = 0; j < nvel; ++j)
                    {
                        stdVelocity[j] = Array<OneD, NekDouble>(n_points);
                    }
                    n_points_0 = n_points;
                }        
                
                Array<TwoD, const NekDouble> gmat = 
                    m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetDerivFactors(ptsKeys);
                
                if (m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetGtype()
                    == SpatialDomains::eDeformed)
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][i]*inarray[0][i+cnt] 
                            + gmat[2][i]*inarray[1][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][i]*inarray[0][i+cnt] 
                            + gmat[3][i]*inarray[1][i+cnt];
                    }
                }
                else
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][0]*inarray[0][i+cnt] 
                            + gmat[2][0]*inarray[1][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][0]*inarray[0][i+cnt] 
                            + gmat[3][0]*inarray[1][i+cnt];
                    }
                }
                
                cnt += n_points;
                
                
                for (int i = 0; i < n_points; i++)
                {
                    pntVelocity = stdVelocity[0][i]*stdVelocity[0][i] 
                        + stdVelocity[1][i]*stdVelocity[1][i];
                    
                    if (pntVelocity>maxV[el])
                    {
                        maxV[el] = pntVelocity;
                    }
                }
                maxV[el] = sqrt(maxV[el]);
            }
        }
        else
        {
            cnt = 0;
            for (int el = 0; el < n_element; ++el)
            { 
                
                int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();
                
                // reset local space if necessary
                if(n_points != n_points_0)
                {
                    for (int j = 0; j < nvel; ++j)
                    {
                        stdVelocity[j] = Array<OneD, NekDouble>(n_points);
                    }
                    n_points_0 = n_points;
                }        
                
                Array<TwoD, const NekDouble> gmat =
                    m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetDerivFactors(ptsKeys);
                
                if (m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetGtype()
                    == SpatialDomains::eDeformed)
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][i]*inarray[0][i+cnt] 
                            + gmat[3][i]*inarray[1][i+cnt] 
                            + gmat[6][i]*inarray[2][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][i]*inarray[0][i+cnt] 
                            + gmat[4][i]*inarray[1][i+cnt] 
                            + gmat[7][i]*inarray[2][i+cnt];
                        
                        stdVelocity[2][i] = gmat[2][i]*inarray[0][i+cnt] 
                            + gmat[5][i]*inarray[1][i+cnt] 
                            + gmat[8][i]*inarray[2][i+cnt];
                    }
                }
                else
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][0]*inarray[0][i+cnt] 
                            + gmat[3][0]*inarray[1][i+cnt] 
                            + gmat[6][0]*inarray[2][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][0]*inarray[0][i+cnt] 
                            + gmat[4][0]*inarray[1][i+cnt] 
                            + gmat[7][0]*inarray[2][i+cnt];
                        
                        stdVelocity[2][i] = gmat[2][0]*inarray[0][i+cnt] 
                            + gmat[5][0]*inarray[1][i+cnt] 
                            + gmat[8][0]*inarray[2][i+cnt];
                    }
                }
                
                cnt += n_points;
                
                for (int i = 0; i < n_points; i++)
                {
                    pntVelocity = stdVelocity[0][i]*stdVelocity[0][i] 
                        + stdVelocity[1][i]*stdVelocity[1][i] 
                        + stdVelocity[2][i]*stdVelocity[2][i];
                    
                    if (pntVelocity > maxV[el])
                    {
                        maxV[el] = pntVelocity;
                    }
                }

                maxV[el] = sqrt(maxV[el]);
                //cout << maxV[el]*maxV[el] << endl;
            }
        }
        
        return maxV;
    }

    /**
     *    At the start, the newest value is stored in array[nlevels-1]
     *        and the previous values in the first positions
     *    At the end, the extrapolated value is stored in array[nlevels-1]
     *        and the storage has been updated to included the new value
     */    
    void Extrapolate::ExtrapolateArray(
            Array<OneD, Array<OneD, NekDouble> > &array)
    {
        int nint     = min(m_pressureCalls,m_intSteps);
        int nlevels  = array.num_elements();
        int nPts     = array[0].num_elements();

        // Update array
        RollOver(array);

        // Extrapolate to outarray
        Vmath::Smul(nPts, StifflyStable_Betaq_Coeffs[nint-1][nint-1],
                         array[nint-1],    1,
                         array[nlevels-1], 1);

        for(int n = 0; n < nint-1; ++n)
        {
            Vmath::Svtvp(nPts, StifflyStable_Betaq_Coeffs[nint-1][n],
                         array[n],1, array[nlevels-1],1,
                         array[nlevels-1],1);
        }
    }

    /**
     *    At the start, the newest value is stored in array[nlevels-1]
     *        and the previous values in the first positions
     *    At the end, the acceleration from BDF is stored in array[nlevels-1]
     *        and the storage has been updated to included the new value
     */
    void Extrapolate::AccelerationBDF(
            Array<OneD, Array<OneD, NekDouble> > &array)
    {
        int nlevels  = array.num_elements();
        int nPts     = array[0].num_elements();

        // Update array
        RollOver(array);

        // Calculate acceleration using Backward Differentiation Formula
        Array<OneD, NekDouble> accelerationTerm (nPts, 0.0);
        if (m_pressureCalls > 2)
        {
            int acc_order = min(m_pressureCalls-2,m_intSteps);
            Vmath::Smul(nPts,
                             StifflyStable_Gamma0_Coeffs[acc_order-1],
                             m_acceleration[0], 1,
                             accelerationTerm,  1);

            for(int i = 0; i < acc_order; i++)
            {
                Vmath::Svtvp(nPts,
                            -1*StifflyStable_Alpha_Coeffs[acc_order-1][i],
                             m_acceleration[i+1], 1,
                             accelerationTerm,    1,
                             accelerationTerm,    1);
            }
        }
        m_acceleration[nlevels-1] = accelerationTerm;
    }

}

