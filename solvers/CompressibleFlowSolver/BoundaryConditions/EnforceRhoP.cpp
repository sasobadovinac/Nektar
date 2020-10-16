///////////////////////////////////////////////////////////////////////////////
//
// File: EnforceRhoP.cpp
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
// Description: Modified Riemann invariant boundary condition.
//              Enforcing the density and pressure at the inflow boundary;
//              Enforcing the pressure at the outflow boundary. 
//              The input can be either VALUE or FILE.  
// 
///////////////////////////////////////////////////////////////////////////////


#include "EnforceRhoP.h"
using namespace std;

namespace Nektar
{

std::string EnforceRhoP::className = GetCFSBndCondFactory().
RegisterCreatorFunction("EnforceRhoP",
                        EnforceRhoP::create,
                        "Riemann invariant boundary condition, \
                        fixing Rho and P.");


EnforceRhoP::EnforceRhoP(
           const LibUtilities::SessionReaderSharedPtr& pSession,
           const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
           const Array<OneD, Array<OneD, NekDouble> >& pTraceNormals,
           const int pSpaceDim,
           const int bcRegion,
           const int cnt)
    : CFSBndCond(pSession, pFields, pTraceNormals, pSpaceDim, bcRegion, cnt)
{

    const MultiRegions::ExpListSharedPtr bndexp =
        m_fields[0]->GetBndCondExpansions()[m_bcRegion];

    //-> Gather a list of index from trace to this boundary
    m_npts = bndexp->GetTotPoints();

    m_bndToTraceMap = Array<OneD, int> (m_npts,-1);

    const Array<OneD, const int> &traceBndMap = m_fields[0]->GetTraceBndMap();

    // Construct a map for the boundary to trace map for easy acess to
    // phys space points
    int cnt1 = 0; 
    for (int e = 0; e < bndexp->GetNumElmts(); ++e)
    {
        int nTracePts = bndexp->GetExp(e)->GetTotPoints();
        
        int id = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[m_offset+e]);
        
        // Loop on the points of the m_bcRegion
        for (int i = 0; i < nTracePts; i++)
        {
            // the ith point in region e
            m_bndToTraceMap[cnt1++] = id+i;
        }
    }

    Array<OneD, Array<OneD, NekDouble> > BCvals(m_fields.num_elements());
    m_bndPhys = Array<OneD, Array<OneD, NekDouble> > (m_fields.num_elements());
    
    for(int i = 0; i < m_fields.num_elements(); ++i)
    {
        m_bndPhys[i] = m_fields[i]->GetBndCondExpansions()[m_bcRegion]
            ->UpdatePhys();
        
        BCvals[i] = Array<OneD, NekDouble>(m_npts);
        Vmath::Vcopy(m_npts, m_bndPhys[i], 1, BCvals[i], 1);
    }
    
    // Set up boudnary required BCs
    m_rhoBC = Array<OneD, NekDouble>(m_npts);
    Vmath::Vcopy(m_npts,BCvals[0],1,m_rhoBC,1);

    m_velBC = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
    // Evaluate velocity on boundary
    for(int i = 0; i < m_spacedim; ++i)
    {
        m_velBC[i] = Array<OneD, NekDouble>(m_npts);
        Vmath::Vcopy(m_npts,BCvals[i+1],1,m_velBC[i],1);
        Vmath::Vdiv(m_npts,m_velBC[i],1,m_rhoBC,1,m_velBC[i],1);
    }
    m_pBC = Array<OneD, NekDouble> (m_npts);
    m_varConv->GetPressure(BCvals, m_pBC);

    // Computing the normal velocity for characteristics coming
    // from outside the computational domain
    m_VnInf = Array<OneD, NekDouble> (m_npts, 0.0);
    for(int i = 0; i < m_spacedim; i++)
    {
        for(int j = 0; j < m_npts; ++j)
        {
            m_VnInf[j] += m_traceNormals[i][m_bndToTraceMap[j]]*m_velBC[i][j];
        }
    }

}


void EnforceRhoP::v_Apply
(
        Array<OneD, Array<OneD, NekDouble> >               &Fwd,
        Array<OneD, Array<OneD, NekDouble> >               &physarray,
        const NekDouble                                    &time)
{
    boost::ignore_unused(physarray,time);
    
    int i, j;
    int nDimensions = m_spacedim;

    Array<OneD, Array<OneD, NekDouble> > FwdBnd(Fwd.num_elements());
    Array<OneD, Array<OneD, NekDouble> > bndPhys(Fwd.num_elements());

    // make a local copy of Fwd along boundary of interest
    for(i = 0; i < Fwd.num_elements(); ++i)
    {
        FwdBnd[i] = Array<OneD, NekDouble> (m_npts);
        for(j =0; j < m_npts; ++j)
        {
            FwdBnd[i][j] = Fwd[i][m_bndToTraceMap[j]];
        }
    }

    // Computing the normal velocity for characteristics coming
    // from inside the computational domain
    Array<OneD, NekDouble > Vn (m_npts, 0.0);
    
    for (i = 0; i < nDimensions; ++i)
    {
        for(j = 0; j < m_npts; ++j)
        {
            Vn[j] += m_traceNormals[i][m_bndToTraceMap[j]]*FwdBnd[i+1][j];
        }
    }
    // divide by density. 
    Vmath::Vdiv(m_npts,Vn,1,FwdBnd[0],1,Vn,1);
    
    // Get speed of sound
    Array<OneD, NekDouble > pressure  (m_npts);
    Array<OneD, NekDouble > soundSpeed(m_npts);
    
    m_varConv->GetPressure(FwdBnd, pressure);
    m_varConv->GetSoundSpeed(FwdBnd, soundSpeed);
    
    // Get Mach. Note: it is computed by Vn/c
    Array<OneD, NekDouble > Mach(m_npts, 0.0);
    Vmath::Vdiv(m_npts, Vn, 1, soundSpeed, 1, Mach, 1);
    Vmath::Vabs(m_npts, Mach, 1, Mach, 1);
    
    // Auxiliary variables
    Array<OneD, NekDouble> velBC(nDimensions, 0.0);

    // L represents properties outside boundary
    // R represents properties inside boundary (numerical state)
    NekDouble rhoL, uL, pL;
    NekDouble EBC, vn, rR, cstar, pstar, rhostar, ustar; 

    NekDouble gamMinOne         = m_gamma - 1.0;
    NekDouble twoOverGamMinOne  = 2.0 / gamMinOne;
    NekDouble gamInv            = 1.0/m_gamma;
    
    // Loop on m_bcRegions
    for (int pnt = 0; pnt < m_npts; ++pnt )
    {
        // Impose inflow Riemann invariant
        if (Vn[pnt] <= 0.0)
        {
            // Subsonic flows
            if (Mach[pnt] < 1.00)
            {
                // right characteristic
                rR = -Vn[pnt] - sqrt(m_gamma * pressure[pnt] /
                                     FwdBnd[0][pnt])*twoOverGamMinOne; 
                vn = -m_VnInf[pnt];

                // fix rhostar and pstar to be the input values
		// compute ustar using left-pointing characteristic line IR^-
		pstar = m_pBC[pnt];
		rhostar = m_rhoBC[pnt];
		cstar = sqrt(m_gamma*pstar/rhostar);
		ustar = rR + cstar * twoOverGamMinOne;
               
                // add supplement equation that rhoL=rhostar
                // then pL=pstar, according to IL^0
                // and  uL=ustar, according to IL^+
		rhoL = rhostar;
		pL = pstar;
		uL = ustar;
            }
            else  // Supersonic inflow
            {
                // all characteristics are from left so just impose
                // star state to left values
                // Note: m_vnInf is the negative of the normal velocity
                // across boundary
                rhoL = m_rhoBC[pnt];
                uL =  -m_VnInf[pnt];
                pL =   m_pBC[pnt];
            }
            
            // Boundary energy
            EBC = pL *twoOverGamMinOne *0.5;

            // evaluate the different between the left state normal
            // velocity and that from the desired condition (note
            // m_VnInf is using an outwards normal definition.
            NekDouble VnDiff = uL + m_VnInf[pnt];

            // Boundary velocities & Kinite energy 
            // Note: normals are negated since they point outwards in
            // the domain

            // Note: Can just use the BC values directly!!
            for ( j = 0; j < nDimensions; ++j)
            {
                // Set velocity to the desired conditions modified to
                // take account of the normal state for Riemann
                // problem. (Negative accounts for outwards normal definition)
                velBC[j] = m_velBC[j][pnt] - VnDiff
                    *m_traceNormals[j][m_bndToTraceMap[pnt]];
                
                EBC += 0.5 * rhoL * velBC[j]*velBC[j];
            }
                
            // Impose Left hand Riemann Invariant boundary conditions
            m_bndPhys[0][pnt] = rhoL;
            for (j = 0; j < nDimensions; ++j)
            {
                m_bndPhys[j+1][pnt] = rhoL * velBC[j];
            }
            m_bndPhys[nDimensions+1][pnt] = EBC;
        }
        else // Outflow
        {

            // Note: Allowing the switch can cause worse convergence in this type BC. 
            //       So improve it later.
            if (Mach[pnt] < 1.00)
            {             
                // subsonic outflow: fix pstar
                rR = -Vn[pnt] - sqrt(m_gamma * pressure[pnt] /
                                     FwdBnd[0][pnt])*twoOverGamMinOne;
                vn = -m_VnInf[pnt];

                pstar   = m_pBC[pnt];
                rhostar = FwdBnd[0][pnt] * pow( (pstar/pressure[pnt]) , gamInv ); 
                cstar   = sqrt( m_gamma*pstar/rhostar );
                ustar   = rR + cstar * twoOverGamMinOne;            

                rhoL = rhostar;
                uL   = ustar;
                pL   = pstar;
            }
            else
            {
                // supersonic outflow 
                // Just set to imposed state and let Riemann BC dictate values
                rhoL = m_rhoBC[pnt];
                uL =  -m_VnInf[pnt];
                pL =   m_pBC[pnt]; 
            }

            // Boundary energy
            EBC = pL *twoOverGamMinOne *0.5;

            // Boundary velocities & Kinite energy 
            // Note: normals are negated since they point outwards in
            // the domain
            for ( j = 0; j < nDimensions; ++j)
            {
                velBC[j] = -1.0 * uL * m_traceNormals[j][m_bndToTraceMap[pnt]];
                EBC += 0.5 * rhoL * velBC[j]*velBC[j];
            }
                
            // Impose Left hand Riemann Invariant boundary conditions
            m_bndPhys[0][pnt] = rhoL;
            for (j = 0; j < nDimensions; ++j)
            {
                m_bndPhys[j+1][pnt] = rhoL * velBC[j];
            }
            m_bndPhys[nDimensions+1][pnt] = EBC;

        }

    }

}

}
