///////////////////////////////////////////////////////////////////////////////
//
// File ExecuateSharpy.cpp
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
// Description: *****.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Forcing/ExecuteSharpy.h>
#include <iostream> 

  
using namespace std;

namespace Nektar
{
namespace Sharpy
{
void Beam::InitialiseStatic(const LibUtilities::SessionReaderSharedPtr& m_session)
{
    //----------------------------Initial Displacement: ImpStart vs Static Solution
    // Initialise static beam data.
    // 1. Read input data.
    m_OutFile = Array<OneD, char>(25,"BeamSol.txt");                  // The file to store the solution of Beam Dynamics
    m_session->LoadParameter("BeamLength", m_BeamLength);			  // Length of the beam
    m_session->LoadParameter("BeamNumNodesElem", m_NumNodesElem,2);
    m_session->LoadParameter("BeamMaxElNod", m_MaxElNod,3);           // Max number of nodes per element.
    m_session->LoadParameter("BeamNumElems", m_NumElems,20);          // Elements discreising the beam
    m_session->LoadParameter("BeamNumForcedDisp", m_NForcedDisp,0);   //
    m_session->LoadParameter("BeamElemProj", m_ElemProj,0);           // Element info computed in (1) global frame (2) fixed element frame (3) moving element frame.
    m_session->LoadParameter("BeamMaxIterations", m_MaxIterations,20);// Newton-Raphson iterations for nonlinear solution
    m_session->LoadParameter("BeamNumLoadSteps", m_NumLoadSteps,60);  // Number of load increments
    m_session->LoadParameter("BeamNumGauss", m_NumGauss,3);           // Gauss points in element integration
    m_session->LoadParameter("BeamSolution", m_Solution,912);         // 902/912: cbeam3 linear/nonlinear flexible-body dynamic 
    m_session->LoadParameter("BeamDimMat", m_DimMat,24);
    m_session->LoadParameter("BeamDeltaCurved", m_DeltaCurved,1.e-4); // Min. angle for two vectors to be parallel
    m_session->LoadParameter("BeamMinDelta", m_MinDelta,1.e-6);       // Convergence param for Newton-Rhaphson iterations
    m_session->LoadParameter("BeamDiameter", m_D,1.0);				  // Diameter of the cylinder
    m_session->LoadParameter("Velocity_Infinity", m_U_inf,1.0);		  // Velocity of free stream
    m_session->LoadParameter("Density_Fluid", m_rho_f,1.0);			  // Density of the fluid

    m_session->MatchSolverInfo("BeamFollowerForce","True",m_FollowerForce,false);
    m_session->MatchSolverInfo("BeamFollowerForceRig","True",m_FollowerForceRig,false);
    m_session->MatchSolverInfo("BeamPrintInfo","True",m_PrintInfo,false);
    m_session->MatchSolverInfo("BeamOutInBframe","True",m_OutInBframe,false);
    m_session->MatchSolverInfo("BeamOutInaframe","True",m_OutInaframe,false);

    m_iStep             = 0;
    m_NumNodes          = Array<OneD, int> (m_NumElems);
    m_MemNo             = Array<OneD, int> (m_NumElems);
    m_Conn              = Array<OneD, int> (m_MaxElNod*m_NumElems);
    m_Master_Array      = Array<OneD, int> (m_MaxElNod*m_NumElems*2);
    m_Length            = Array<OneD, NekDouble> (m_NumElems,0.0);
    m_PreCurv           = Array<OneD, NekDouble> (3*m_NumElems,0.0);
    m_Psi               = Array<OneD, NekDouble> (3*m_NumElems,0.0);
    m_Vector            = Array<OneD, NekDouble> (3*m_NumElems,0.0);
    m_Mass_Array        = Array<OneD, NekDouble> (6*m_NumElems*6,0.0);
    m_Stiff_Array       = Array<OneD, NekDouble> (6*m_NumElems*6,0.0);
    m_InvStiff_Array    = Array<OneD, NekDouble> (6*m_NumElems*6,0.0);
    m_RBMass_Array      = Array<OneD, NekDouble> (m_MaxElNod*m_NumElems*6*6,0.0);

    // Read element information of the beam
    Input_elem(m_session);

    m_BoundConds        = Array<OneD, int> (m_NumNodes_tot);
    m_Nod_Sflag         = Array<OneD, int> (m_NumNodes_tot);
    m_Master            = Array<OneD, int> (2*m_NumNodes_tot);
    m_Vdof              = Array<OneD, int> (m_NumNodes_tot);
    m_Fdof              = Array<OneD, int> (m_NumNodes_tot);
    m_ListIN            = Array<OneD, int> (m_NumNodes_tot);
    m_NodeForcedDisp_Array  = Array<OneD, int> (m_NForcedDisp);
    m_PosForcedDisp_Array   = Array<OneD, NekDouble> (m_NForcedDisp*3,0.0);
    m_PosIni_Array          = Array<OneD, NekDouble> (m_NumNodes_tot*3,0.0);
    m_ForceStatic_Array     = Array<OneD, NekDouble> (m_NumNodes_tot*6,0.0);
    m_PhiNodes              = Array<OneD, NekDouble> (m_NumNodes_tot,0.0);
    m_PsiIni_Array          = Array<OneD, NekDouble> (m_NumElems*m_MaxElNod*3,0.0);
    m_PosDefor_Array        = Array<OneD, NekDouble> (m_NumNodes_tot*3,0.0);
    m_PsiDefor_Array        = Array<OneD, NekDouble> (m_NumElems*m_MaxElNod*3,0.0);
    m_PosDeforDot_Array     = Array<OneD, NekDouble> (m_NumNodes_tot*3,0.0);
    m_PsiDeforDot_Array     = Array<OneD, NekDouble> (m_NumElems*m_MaxElNod*3,0.0);
    m_PosDeforDDot_Array    = Array<OneD, NekDouble> (m_NumNodes_tot*3,0.0);
    m_PsiDeforDDot_Array    = Array<OneD, NekDouble> (m_NumElems*m_MaxElNod*3,0.0);

    m_Quat_Array            = Array<OneD, NekDouble> (4,0.0);
    m_PsiA_G_Array          = Array<OneD, NekDouble> (3,0.0);
    m_Cao_Array             = Array<OneD, NekDouble> (3*3,0.0);

    m_Quat_Array[0]         = 1.0;
    m_Cao_Array[0]          = 1.0;
    m_Cao_Array[4]          = 1.0;
    m_Cao_Array[8]          = 1.0;

    
    // Read node informatioin of the beam
    Input_node(m_session);

    // 2. Compute initial (undeformed) geometry.
    SHARPy::Wrap_xbeam_undef_geom(
        m_NumElems,m_NumNodes_tot,
        &m_NumNodes[0],&m_MemNo[0],
        &m_Conn[0],&m_Master_Array[0],
        &m_Length[0],&m_PreCurv[0],
        &m_Psi[0],&m_Vector[0],&m_Mass_Array[0],
        &m_Stiff_Array[0],&m_InvStiff_Array[0],
        &m_RBMass_Array[0],&m_PosIni_Array[0],
        &m_PhiNodes[0],&m_PsiIni_Array[0],
        m_FollowerForce, m_FollowerForceRig,
        m_PrintInfo, m_OutInBframe, m_OutInaframe,
        m_ElemProj, m_MaxIterations, m_NumLoadSteps,
        m_NumGauss, m_Solution, m_DeltaCurved,
        m_MinDelta, m_NewmarkDamp);

    // 3. Identify nodal degrees of freedom.
    SHARPy::Wrap_xbeam_undef_dofs(
        m_NumElems,m_NumNodes_tot,
        &m_NumNodes[0],&m_MemNo[0],
        &m_Conn[0],&m_Master_Array[0],
        &m_Length[0],&m_PreCurv[0],
        &m_Psi[0],&m_Vector[0],&m_Mass_Array[0],
        &m_Stiff_Array[0],&m_InvStiff_Array[0],
        &m_RBMass_Array[0],&m_BoundConds[0],
        &m_Master[0],&m_Vdof[0],&m_Fdof[0],m_NumDof,
        &m_Nod_Sflag[0]);

    // Add RBMasses (Rigid-Body masses)
    ////

    //Calculate initial displacements.
    Vmath::Vcopy(m_NumNodes_tot*3, m_PosIni_Array, 1, m_PosDefor_Array, 1);
    Vmath::Vcopy(m_NumElems*m_MaxElNod*3, m_PsiIni_Array, 1, m_PsiDefor_Array, 1);
}

/**
 *
 **/
void Beam::InitialiseDynamic(const LibUtilities::SessionReaderSharedPtr& m_session,
        const Array<OneD, Array<OneD, NekDouble> >& HydCoeffs)
{
     
    // 1. Input data for transient dynamic solution.    
    m_session->LoadParameter("BeamNewmarkDamp", m_NewmarkDamp,0.01);
    m_session->LoadParameter("BeamTime0", m_t0,0.0);                

    m_session->LoadParameter("BeamTimefinal", m_tf,99999);          
    m_session->LoadParameter("BeamDT", m_dt);                     

    m_NumDof2       = m_NumDof*m_NumDof;
    m_X             = Array<OneD, NekDouble> (m_NumDof,0.0);
    m_DX            = Array<OneD, NekDouble> (m_NumDof,0.0);
    m_DXdt          = Array<OneD, NekDouble> (m_NumDof,0.0);
    m_DXddt         = Array<OneD, NekDouble> (m_NumDof,0.0);
    m_Force_Array   = Array<OneD, NekDouble> (m_NumNodes_tot*6,0.0);
    m_Vrel          = Array<OneD, NekDouble> (6,0.0);
    m_VrelDot       = Array<OneD, NekDouble> (6,0.0);
    m_Mglobal_Array = Array<OneD, NekDouble> (m_NumDof2,0.0);
    m_Mvel_Array    = Array<OneD, NekDouble> (m_NumDof*6,0.0);
    m_Cglobal_Array = Array<OneD, NekDouble> (m_NumDof2,0.0);
    m_Cvel_Array    = Array<OneD, NekDouble> (m_NumDof*6,0.0);
    m_Kglobal_Array = Array<OneD, NekDouble> (m_NumDof2,0.0);
    m_Fglobal_Array = Array<OneD, NekDouble> (m_NumDof2,0.0);
    m_Qglobal       = Array<OneD, NekDouble> (m_NumDof,0.0);
    m_Asys_Array    = Array<OneD, NekDouble>(m_NumDof2,0.0);

    m_gamma         = 1.0/2.0+m_NewmarkDamp;
    m_beta          = 1.0/4.0*(m_gamma+0.5)*(m_gamma+0.5);
    m_NumSteps      = (m_tf - m_t0)/m_dt;

    m_Time          = Array<OneD, NekDouble> (m_NumSteps);
    for(int i = 0; i < m_NumSteps; i++)
    {
        m_Time[i] = m_t0 + m_dt*NekDouble(i);
    }

    //Initialize Local variables.
    for(int k = 0; k < m_NumNodes_tot; k++)
    {
        m_ListIN[k] = m_Vdof[k];
    }
   
    Array<OneD, NekDouble> Temp(m_NumNodes_tot,0.0);

    Vmath::Vcopy(m_NumNodes_tot,HydCoeffs[0],1,Temp = m_Force_Array +   m_NumNodes_tot,1);
    Vmath::Vcopy(m_NumNodes_tot,HydCoeffs[1],1,Temp = m_Force_Array + 2*m_NumNodes_tot,1);

    // dimensionalize the hydroforces applying on the cylinder.
    NekDouble C = m_rho_f*m_U_inf*m_U_inf*m_D;
    C = m_BeamLength/(m_NumNodes_tot-1)*C;
    Vmath::Smul(6*m_NumNodes_tot,C,m_Force_Array,1,m_Force_Array,1);
    for(int i = 0; i < 3; i++)
    {
        m_Force_Array[i*m_NumNodes_tot]         = 0.5*m_Force_Array[i*m_NumNodes_tot];
        m_Force_Array[(i+1)*m_NumNodes_tot-1]   = 0.5*m_Force_Array[(i+1)*m_NumNodes_tot-1];
    }


    Vmath::Zero(m_NumDof2,m_Asys_Array,1);
    Vmath::Zero(m_NumDof2,m_Mglobal_Array,1);
    Vmath::Zero(m_NumDof2,m_Cglobal_Array,1);
    Vmath::Zero(m_NumDof2,m_Kglobal_Array,1);
    Vmath::Zero(m_NumDof2,m_Fglobal_Array,1);
    Vmath::Zero(m_NumDof,m_Qglobal,1);

    //Extract initial displacements and velocities.
    SHARPy::Wrap_cbeam3_solv_disp2state(
        m_NumNodes_tot,m_NumDof,m_NumElems,
        &m_Master[0],&m_Vdof[0],&m_Fdof[0],
        &m_PosDefor_Array[0],&m_PsiDefor_Array[0],
        &m_PosDeforDot_Array[0],&m_PsiDeforDot_Array[0],
        &m_X[0],&m_DXdt[0]);

    //Compute initial acceleration (we are neglecting qdotdot in Kmass).
    Vmath::Smul(m_NumNodes_tot*3,0.0,m_PosDeforDot_Array,1,m_PosDeforDDot_Array,1);
    Vmath::Smul(m_NumElems*m_MaxElNod*3,0.0,m_PsiDeforDot_Array,1,m_PsiDeforDDot_Array,1);
    Vmath::Vadd(m_NumNodes_tot*6,m_Force_Array,1,m_ForceStatic_Array,1,m_Force_Array,1);

    SHARPy::Wrap_cbeam3_asbly_dynamic(
        m_NumElems,m_NumNodes_tot,&m_NumNodes[0],
        &m_MemNo[0],&m_Conn[0],&m_Master_Array[0],&m_Length[0],
        &m_PreCurv[0],&m_Psi[0],&m_Vector[0],&m_Mass_Array[0],
        &m_Stiff_Array[0],&m_InvStiff_Array[0],&m_RBMass_Array[0],
        &m_Master[0],&m_Vdof[0],&m_Fdof[0],
        &m_PosIni_Array[0],&m_PsiIni_Array[0],
        &m_PosDefor_Array[0],&m_PsiDefor_Array[0],
        &m_PosDeforDot_Array[0],&m_PsiDeforDot_Array[0],
        &m_PosDeforDDot_Array[0],&m_PsiDeforDDot_Array[0],
        &m_Force_Array[0],&m_Vrel[0],&m_VrelDot[0],
        m_NumDof,m_DimMat,
        m_ms,&m_Mglobal_Array[0],&m_Mvel_Array[0],
        m_cs,&m_Cglobal_Array[0],&m_Cvel_Array[0],
        m_ks,&m_Kglobal_Array[0],
        m_fs,&m_Fglobal_Array[0],&m_Qglobal[0],
        m_FollowerForce, m_FollowerForceRig,
        m_PrintInfo, m_OutInBframe, m_OutInaframe,
        m_ElemProj, m_MaxIterations, m_NumLoadSteps,
        m_NumGauss, m_Solution, m_DeltaCurved,
        m_MinDelta, m_NewmarkDamp,
        &m_Cao_Array[0]);

    Array<OneD, NekDouble> Temp0_Array(m_NumDof,0.0);
    Array<OneD, NekDouble> Temp1_Array(m_NumDof,0.0);

    SHARPy::Wrap_fem_m2v(
        m_NumNodes_tot,6,
        &m_Force_Array[0],m_NumDof,
        &Temp0_Array[0],&m_ListIN[0]);

    DNekMatSharedPtr fglobal
        = MemoryManager<DNekMat>::AllocateSharedPtr(
            m_NumDof,m_NumDof);

    for(int i = 0; i < m_NumDof; i++)
    {
        for(int j = 0; j < m_NumDof; j++)
        {
            (*fglobal)(i,j)= m_Fglobal_Array[i*m_NumDof+j];
        }
    }

    Blas::Dgemv('N', m_NumDof, m_NumDof,
                1.0, &(fglobal->GetPtr())[0],
                m_NumDof, &Temp0_Array[0], 1,
                0.0,  &Temp1_Array[0],1);

    Vmath::Vsub(m_NumDof,m_Qglobal,1,Temp1_Array,1,m_Qglobal,1);
    Vmath::Smul(m_NumDof,-1.0,m_Qglobal,1,m_Qglobal,1);
    Vmath::Vadd(m_NumDof2,m_Mglobal_Array,1,m_Asys_Array,1,m_Asys_Array,1);

    int N = m_NumDof;
    Array<OneD, int> ipivot (N);
    int info =0;
    Lapack::Dgetrf(N, N, m_Asys_Array.get(), N, ipivot.get(), info);
    if( info < 0 )
    {
        std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) +
                                "th parameter had an illegal parameter for dgetrf";
    }
    else if( info > 0 )
    {
        std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
        boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
        ASSERTL0(false, message.c_str());
    }

    // N means no transponse (direct matrix)
    int ncolumns_b =1;
    Vmath::Vcopy(m_NumDof, m_Qglobal, 1, m_DXddt, 1);
    Lapack::Dgetrs('N', N, ncolumns_b, m_Asys_Array.get(), N, ipivot.get(), m_DXddt.get(), N, info);
    if( info < 0 )
    {
        std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) +
        "th parameter had an illegal parameter for dgetrf";
        ASSERTL0(false, message.c_str());
    }
    else if( info > 0 )
    {
        std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
        boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
        ASSERTL0(false, message.c_str());
    }
}

/**
 *
 **/
void Beam::SolvenlnStatic(const LibUtilities::SessionReaderSharedPtr& m_session,
        const Array<OneD, Array<OneD, NekDouble> >& HydCoeffs)
{
    Array<OneD, NekDouble> Temp(m_NumNodes_tot,0.0);
	
	// Store the forces in x and y directions in the force array.
    Vmath::Vcopy(m_NumNodes_tot,HydCoeffs[0],1,Temp = m_ForceStatic_Array +   m_NumNodes_tot,1);
    Vmath::Vcopy(m_NumNodes_tot,HydCoeffs[1],1,Temp = m_ForceStatic_Array + 2*m_NumNodes_tot,1);

    // Dimensionalize the hydroforces.
    NekDouble C = m_rho_f*m_U_inf*m_U_inf*m_D*m_BeamLength/(m_NumNodes_tot-1);
    Vmath::Smul(6*m_NumNodes_tot,C,m_ForceStatic_Array,1,m_ForceStatic_Array,1);

    for(int i = 0; i < 3; i++)
    {
        m_ForceStatic_Array[i*m_NumNodes_tot]         = 0.5*m_ForceStatic_Array[i*m_NumNodes_tot];
        m_ForceStatic_Array[(i+1)*m_NumNodes_tot-1]   = 0.5*m_ForceStatic_Array[(i+1)*m_NumNodes_tot-1];
    }

    ///*chaplin_test_#3
    /*
    // Cable length [meters]
    // coordinates of support points [meters]
    NekDouble posB[3] = {-1.8,0.0,0.028}; // support point bottom
    NekDouble posT[3] = {1.8,0.0,2.2}   ; // support point top

    // translate reference system origin to B
    for(int i = 0; i < 3; i++)
    {
        posT[i]=posT[i]-posB[i];
        posB[i]=0.0*posB[i];
    }

    //# distance between support points
    NekDouble Dist_T_B;
    Dist_T_B = sqrt((posT[0]-posB[0])*(posT[0]-posB[0])+
                    (posT[1]-posB[1])*(posT[1]-posB[1])+
                    (posT[2]-posB[2])*(posT[2]-posB[2]));

    //# direction B to T
    NekDouble fi      = -atan((posT[2]-posB[2])/(posT[0]-posB[0]));
    m_PsiA_G_Array[1] = 1;
    Vmath::Smul(3, fi, m_PsiA_G_Array, 1, m_PsiA_G_Array, 1);

    Array<OneD, int>  ForcedDisp_node(m_NForcedDisp);
    Array<OneD, char> ForcedDisp_FoR (m_NForcedDisp);

    Array<OneD, Array<OneD, NekDouble> > ForcedDisp_pos(m_NForcedDisp);
    for(int i = 0; i < m_NForcedDisp; i++)
    {
        ForcedDisp_pos [i] = Array<OneD, NekDouble>(3);
    }

    for(int i = 0; i < m_NForcedDisp; i++)
    {
        ForcedDisp_node[i]    = m_NumNodes_tot - 1;
        ForcedDisp_pos [i][0] = Dist_T_B;
        ForcedDisp_pos [i][1] = 0.0;
        ForcedDisp_pos [i][2] = 0.0;
        ForcedDisp_FoR [i]    = 'A';
        if (ForcedDisp_FoR[i] == 'G')
        {
            ASSERTL0(false, "Only the FoR='A' option is "
                "available to specify forced displacements")
        }
        m_NodeForcedDisp_Array[i]   = ForcedDisp_node[i]  ;
        m_PosForcedDisp_Array[3*i]  = ForcedDisp_pos[i][0];
        m_PosForcedDisp_Array[3*i+1]= ForcedDisp_pos[i][1];
        m_PosForcedDisp_Array[3*i+2]= ForcedDisp_pos[i][2];
    }
    */
    //Add Gravity Force
    bool IsGravity;
    m_session->MatchSolverInfo("IsBeamGravity","True",IsGravity,false);
    if(IsGravity)
    {
        AddGravityLoad(m_session);
    }

    SHARPy::Wrap_cbeam3_solv_nlnstatic (
        m_NumDof,m_NumElems,
        &m_NumNodes[0], &m_MemNo[0], &m_Conn[0],
        &m_Master_Array[0],
        &m_Length[0], &m_PreCurv[0],
        &m_Psi[0], &m_Vector[0], &m_Mass_Array[0],
        &m_Stiff_Array[0],
        &m_InvStiff_Array[0], &m_RBMass_Array[0],
        m_NumNodes_tot,
        &m_Master[0], &m_Vdof[0], &m_Fdof[0],
        &m_ForceStatic_Array[0],
        &m_PosIni_Array[0],&m_PsiIni_Array[0],
        &m_PosDefor_Array[0],&m_PsiDefor_Array[0],
        m_FollowerForce, m_FollowerForceRig,
        m_PrintInfo, m_OutInBframe, m_OutInaframe,
        m_ElemProj, m_MaxIterations, m_NumLoadSteps,
        m_NumGauss, m_Solution, m_DeltaCurved,
        m_MinDelta, m_NewmarkDamp,
        m_NForcedDisp,&m_NodeForcedDisp_Array[0],
        &m_PosForcedDisp_Array[0],&m_PsiA_G_Array[0]);

    /*
    Array<OneD, Array<OneD, NekDouble> > Cao(3);
    Cao[0]= Array<OneD, NekDouble>(3,0.0);
    Cao[1]= Array<OneD, NekDouble>(3,0.0);
    Cao[2]= Array<OneD, NekDouble>(3,0.0);

    RotCRV(m_PsiA_G_Array,Cao);
    
    Array<OneD, Array<OneD, NekDouble> > PosDefG(m_NumNodes_tot);
    for(int i = 0; i < m_NumNodes_tot; i++)
    {
        PosDefG[i] = Array<OneD, NekDouble>(3,0.0);
    }

    for(int i = 0; i < m_NumNodes_tot; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                PosDefG[i][j] += Cao[j][k]*m_PosDefor_Array[3*i+k];
            }
        }
    }
    */
}


/**
 *
 **/
void Beam::SolvenlnDynamic(const LibUtilities::SessionReaderSharedPtr& m_session,
                const Array<OneD, Array<OneD, NekDouble> > &HydroForces,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &motions)
{
    //
    int nvibmodes;
    m_session->LoadParameter("VibModesZ", nvibmodes);
	Array<OneD, NekDouble> tp0(nvibmodes+1,0.0);
	Array<OneD, NekDouble> tp1(nvibmodes+1,0.0);
	Vmath::Vcopy(nvibmodes+1,motions[0][1],1,tp0,1);
	Vmath::Vcopy(nvibmodes+1,motions[1][1],1,tp1,1);

    Array<OneD, NekDouble> Temp(m_NumNodes_tot,0.0);
    Vmath::Vcopy(m_NumNodes_tot,HydroForces[0],1,Temp = m_Force_Array +   m_NumNodes_tot,1);
    Vmath::Vcopy(m_NumNodes_tot,HydroForces[1],1,Temp = m_Force_Array + 2*m_NumNodes_tot,1);
	
	
    // dimensionalize the hydroforces.
    NekDouble C = m_rho_f*m_U_inf*m_U_inf*m_D*m_BeamLength/(m_NumNodes_tot-1);
    Vmath::Smul(6*m_NumNodes_tot,C,m_Force_Array,1,m_Force_Array,1);

    for(int i = 0; i < 3; i++)
    {
        m_Force_Array[i*m_NumNodes_tot]         = 0.5*m_Force_Array[i*m_NumNodes_tot];
        m_Force_Array[(i+1)*m_NumNodes_tot-1]   = 0.5*m_Force_Array[(i+1)*m_NumNodes_tot-1];
    }
	

    m_dt = m_Time[m_iStep+1]-m_Time[m_iStep];

    int size = 4;
    Array<OneD, NekDouble> Temp0_Array(size*size,0.0);
    Array<OneD, NekDouble> Temp1_Array(size*size,0.0);
    Array<OneD, NekDouble> Temp2_Array(size*size,0.0);
    Array<OneD, NekDouble> Temp3_Array(size*size,0.0);
    Array<OneD, NekDouble> Unit4_Array(size*size,0.0);

    for(int i = 0; i < 4; i++)
    {
        Unit4_Array[i*4+i] = 1.0;
    }

    NekDouble tmp = 0.25*(m_Time[m_iStep+1]-m_Time[m_iStep]);

    // Update transformation matrix for given angular velocity
    SHARPy::Wrap_xbeam_QuadSkew(4, &m_Vrel[3], &Temp0_Array[0]);
    Vmath::Svtvp(size*size, tmp, Temp0_Array, 1, Unit4_Array, 1, Temp1_Array, 1);
    SHARPy::Wrap_lu_invers(size, &Temp1_Array[0], &Temp2_Array[0]);
    Vmath::Svtvp(size, -1.0*tmp, Temp0_Array, 1, Unit4_Array, 1,  Temp1_Array, 1);
    SHARPy::Wrap_matmul(size, size, size, &Temp1_Array[0], &m_Quat_Array[0], &Temp3_Array[0]);
    SHARPy::Wrap_matmul(size, size, size, &Temp2_Array[0], &Temp3_Array[0], &m_Quat_Array[0]);
    SHARPy::Wrap_xbeam_Rot(&m_Quat_Array[0], &m_Cao_Array[0]);

    // Predictor step.
    Array<OneD,NekDouble> Temp_Array(m_NumDof, 0.0);
    Vmath::Smul(m_NumDof, (0.5-m_beta)*m_dt*m_dt, m_DXddt, 1, Temp_Array, 1);
    Vmath::Svtvp(m_NumDof, m_dt, m_DXdt, 1, Temp_Array, 1, Temp_Array, 1);
    Vmath::Vadd(m_NumDof, m_X, 1, Temp_Array, 1, m_X, 1);
    Vmath::Svtvp(m_NumDof, (1.0-m_gamma)*m_dt, m_DXddt, 1, m_DXdt, 1, m_DXdt, 1);
    Vmath::Zero(m_NumDof, m_DXddt, 1);

    int Iter;
    // Iteration until convergence.
    for(Iter = 0; Iter < m_MaxIterations; Iter++)
    {
        if (Iter == m_MaxIterations)
        {
            ASSERTL0(false,"Solution did not converge in SHARPy");
        }

        // Update nodal positions and velocities .
        SHARPy::Wrap_cbeam3_solv_state2disp(
            m_NumElems,m_NumNodes_tot,&m_NumNodes[0],
            &m_MemNo[0],&m_Conn[0],&m_Master_Array[0],
            &m_Length[0],&m_PreCurv[0],&m_Psi[0],
            &m_Vector[0],&m_Mass_Array[0],&m_Stiff_Array[0],
            &m_InvStiff_Array[0],&m_RBMass_Array[0],
            &m_Master[0],&m_Vdof[0],&m_Fdof[0],
            &m_PosIni_Array[0],&m_PsiIni_Array[0],
            m_NumDof,&m_X[0],&m_DXdt[0],&m_PosDefor_Array[0],
            &m_PsiDefor_Array[0],&m_PosDeforDot_Array[0],
            &m_PsiDeforDot_Array[0]);

        // Compute system functionals and matrices. (Use initial accelerations for Kgyr).
        Vmath::Zero(m_NumDof, m_Qglobal, 1);
        Vmath::Zero(m_NumDof*6, m_Mvel_Array, 1);
        Vmath::Zero(m_NumDof*6, m_Cvel_Array, 1);
        Vmath::Zero(m_NumDof2,m_Mglobal_Array,1);
        Vmath::Zero(m_NumDof2,m_Cglobal_Array,1);
        Vmath::Zero(m_NumDof2,m_Kglobal_Array,1);
        Vmath::Zero(m_NumDof2,m_Fglobal_Array,1);
        Vmath::Smul(m_NumNodes_tot*3,0.0,m_PosDefor_Array,1,m_PosDeforDDot_Array,1);
        Vmath::Smul(m_NumElems*m_MaxElNod*3,0.0,m_PsiDefor_Array,1,m_PsiDeforDDot_Array,1);
        Vmath::Vadd(m_NumNodes_tot*6,m_Force_Array,1,m_ForceStatic_Array,1,m_Force_Array,1);

        SHARPy::Wrap_cbeam3_asbly_dynamic(
            m_NumElems,m_NumNodes_tot,&m_NumNodes[0],
            &m_MemNo[0],&m_Conn[0],&m_Master_Array[0],&m_Length[0],
            &m_PreCurv[0],&m_Psi[0],&m_Vector[0],&m_Mass_Array[0],
            &m_Stiff_Array[0],&m_InvStiff_Array[0],&m_RBMass_Array[0],
            &m_Master[0],&m_Vdof[0],&m_Fdof[0],
            &m_PosIni_Array[0],&m_PsiIni_Array[0],
            &m_PosDefor_Array[0],&m_PsiDefor_Array[0],
            &m_PosDeforDot_Array[0],&m_PsiDeforDot_Array[0],
            &m_PosDeforDDot_Array[0],&m_PsiDeforDDot_Array[0],
            &m_Force_Array[0],&m_Vrel[0],&m_VrelDot[0],m_NumDof,m_DimMat,
            m_ms,&m_Mglobal_Array[0],&m_Mvel_Array[0],
            m_cs,&m_Cglobal_Array[0],&m_Cvel_Array[0],
            m_ks,&m_Kglobal_Array[0],
            m_fs,&m_Fglobal_Array[0],&m_Qglobal[0],
            m_FollowerForce, m_FollowerForceRig,
            m_PrintInfo, m_OutInBframe, m_OutInaframe,
            m_ElemProj, m_MaxIterations, m_NumLoadSteps,
            m_NumGauss, m_Solution, m_DeltaCurved,
            m_MinDelta, m_NewmarkDamp,
            &m_Cao_Array[0]);

        // Compute admissible error.
        Array<OneD, NekDouble>tmp0(m_NumDof,0.0);
        Array<OneD, NekDouble>tmp1(m_NumDof,0.0);
        Vmath::Vabs(m_NumDof, m_Qglobal, 1, tmp0, 1);
        NekDouble MinDelta;
        MinDelta = m_MinDelta*std::max(1.0, Vmath::Vmax(m_NumDof,tmp0,1));

        // Compute the residual.
        Array<OneD,NekDouble> Temp4_Array(m_NumDof,0.0);
        Array<OneD,NekDouble> Temp5_Array(m_NumDof,0.0);
        Array<OneD,NekDouble> Temp6_Array(m_NumDof,0.0);
        Array<OneD,NekDouble> Temp7_Array(m_NumDof,0.0);

        DNekMatSharedPtr mglobal =
            MemoryManager<DNekMat>::AllocateSharedPtr(
                m_NumDof,m_NumDof);

        for(int i = 0; i < m_NumDof; i++)
        {
            for(int j = 0; j < m_NumDof; j++)
            {
                (*mglobal)(i,j) = m_Mglobal_Array[i*m_NumDof+j];
            }
        }

        Blas::Dgemv('N', m_NumDof, m_NumDof,
                    1.0, &(mglobal->GetPtr())[0],
                    m_NumDof, &m_DXddt[0], 1,
                    0.0, &Temp4_Array[0], 1);

        SHARPy::Wrap_matmul(
            m_NumDof,6,1,
            &m_Mvel_Array[0], &m_VrelDot[0],
            &Temp5_Array[0]);

        Vmath::Vadd(m_NumDof, Temp4_Array, 1, Temp5_Array, 1, Temp5_Array, 1);
        Vmath::Vadd(m_NumDof, Temp5_Array, 1, m_Qglobal,   1, m_Qglobal,   1);

        SHARPy::Wrap_fem_m2v(
            m_NumNodes_tot,6,
            &m_Force_Array[0],
            m_NumDof, &Temp6_Array[0],
            &m_ListIN[0]);

        DNekMatSharedPtr fglobal =
            MemoryManager<DNekMat>::AllocateSharedPtr(
                m_NumDof,m_NumDof);

        for(int i = 0; i < m_NumDof; i++)
        {
            for(int j = 0; j < m_NumDof; j++)
            {
                (*fglobal)(i,j)= m_Fglobal_Array[i*m_NumDof+j];
            }
        }

        Blas::Dgemv('N', m_NumDof, m_NumDof,
                    1.0, &(fglobal->GetPtr())[0],
                    m_NumDof, &Temp6_Array[0], 1,
                    0.0, &Temp7_Array[0], 1);

        Vmath::Vsub(m_NumDof, m_Qglobal, 1, Temp7_Array, 1, m_Qglobal, 1);

        // Check convergence.
        Vmath::Vabs(m_NumDof, m_Qglobal, 1, tmp0, 1);
        Vmath::Vabs(m_NumDof, m_DX,      1, tmp1, 1);
        Vmath::Vadd(m_NumDof, tmp0,      1, tmp1, 1, tmp1, 1);

        if (Vmath::Vmax(m_NumDof,tmp1,1) < MinDelta)
        {
            if (m_PrintInfo)
            {
                cout<<"Subiteration = "<<Iter<<"; Delta = "
                    <<Vmath::Vmax(m_NumDof,tmp0,1)<<endl;
            }
            break;
        }

        // Calculate Jacobian
        NekDouble factor1, factor2;
        factor1 = m_gamma/(m_beta*m_dt);
        factor2 = 1.0/(m_beta*m_dt*m_dt);

        Vmath::Zero(m_NumDof2,m_Asys_Array,1);
        Vmath::Vadd(m_NumDof2,m_Kglobal_Array,1,m_Asys_Array,1,m_Asys_Array,1);
        Vmath::Svtvp(m_NumDof2,factor1,m_Cglobal_Array,1,m_Asys_Array,1,m_Asys_Array,1);
        Vmath::Svtvp(m_NumDof2,factor2,m_Mglobal_Array,1,m_Asys_Array,1,m_Asys_Array,1);

        // Calculation of the correction.
        Vmath::Smul(m_NumDof, -1.0, m_Qglobal, 1, m_Qglobal, 1);

        Array<OneD, int> ipivot (m_NumDof);
        int info = 0;
        Lapack::Dgetrf(m_NumDof, m_NumDof, m_Asys_Array.get(), m_NumDof, ipivot.get(), info);
        if( info < 0 )
        {
            std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) +
                        "th parameter had an illegal parameter for dgetrf";
        }
        else if( info > 0 )
        {
            std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
                boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
                ASSERTL0(false, message.c_str());
        }

        // N means no transponse (direct matrix)
        int ncolumns_b = 1;
        Vmath::Vcopy(m_NumDof, m_Qglobal,1,m_DX,1);
        Lapack::Dgetrs( 'N', m_NumDof, ncolumns_b, m_Asys_Array.get(),
                m_NumDof, ipivot.get(), m_DX.get(), m_NumDof, info);
        if( info < 0 )
        {
            std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) +
                "th parameter had an illegal parameter for dgetrf";
                ASSERTL0(false, message.c_str());
        }
        else if( info > 0 )
        {
            std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
                boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
                ASSERTL0(false, message.c_str());
        }

        //Corrector step
        Vmath::Vadd(m_NumDof,  m_X, 1,  m_DX, 1, m_X, 1);
        Vmath::Svtvp(m_NumDof, factor1, m_DX, 1, m_DXdt,  1, m_DXdt,  1);
        Vmath::Svtvp(m_NumDof, factor2, m_DX, 1, m_DXddt, 1, m_DXddt, 1);
    }

    // Update nodal positions and velocities on the current converged time step.
    SHARPy::Wrap_cbeam3_solv_state2disp(
        m_NumElems,m_NumNodes_tot,&m_NumNodes[0],
        &m_MemNo[0],&m_Conn[0],&m_Master_Array[0],
        &m_Length[0],&m_PreCurv[0],&m_Psi[0],
        &m_Vector[0],&m_Mass_Array[0],
        &m_Stiff_Array[0],&m_InvStiff_Array[0],
        &m_RBMass_Array[0],&m_Master[0],&m_Vdof[0],&m_Fdof[0],
        &m_PosIni_Array[0],&m_PsiIni_Array[0],
        m_NumDof,&m_X[0],&m_DXdt[0],
        &m_PosDefor_Array[0],&m_PsiDefor_Array[0],
        &m_PosDeforDot_Array[0],&m_PsiDeforDot_Array[0]);
    
    //copy the displacement and velocity for current step
    Array<OneD, NekDouble> temp(3*m_NumNodes_tot,0.0);

    Vmath::Vsub(3*m_NumNodes_tot,  m_PosDefor_Array, 1,m_PosIni_Array,1,temp,1);
     
    for(int i = 0; i < nvibmodes+1; i++)
    {
        int inode;
        if(m_NumNodesElem == 2)
        {
            inode = i;
        }
        else //m_NumNodesElem == 3 
        {
            inode = 2*i;
        }
        motions[0][0][i] =  temp[  m_NumNodes_tot + inode];
        motions[1][0][i] =  temp[2*m_NumNodes_tot + inode];
        motions[0][1][i] =  m_PosDeforDot_Array[  m_NumNodes_tot + inode];
        motions[1][1][i] =  m_PosDeforDot_Array[2*m_NumNodes_tot + inode];
    }

    //scale the displacement and velocity
    Vmath::Smul(nvibmodes+1,1.0/m_D,motions[0][0],1,motions[0][0],1);
    Vmath::Smul(nvibmodes+1,1.0/m_D,motions[1][0],1,motions[1][0],1);
    Vmath::Smul(nvibmodes+1,1.0/m_U_inf,motions[0][1],1,motions[0][1],1);
    Vmath::Smul(nvibmodes+1,1.0/m_U_inf,motions[1][1],1,motions[1][1],1);
	
	//calculate accel.
	Vmath::Vsub(nvibmodes+1,motions[0][1],1,tp0,1,tp0,1);
	Vmath::Vsub(nvibmodes+1,motions[1][1],1,tp1,1,tp1,1);
	Vmath::Smul(nvibmodes+1,1.0/m_dt,tp0,1,motions[0][2],1);
	Vmath::Smul(nvibmodes+1,1.0/m_dt,tp1,1,motions[1][2],1);

    m_iStep++;
}

/**
 *
 **/
void Beam::Input_elem(const LibUtilities::SessionReaderSharedPtr& m_session)
{
    // Connectivies.
    Vmath::Zero(m_MaxElNod*m_NumElems,m_Conn,1);
    switch(m_NumNodesElem)
    {
        case 2:
            {
                for(int i = 0; i < m_NumElems; i++)
                {
                    m_Conn[m_MaxElNod*i]   = i+1;
                    m_Conn[m_MaxElNod*i+1] = i+2;
                    m_NumNodes[i] = 2;
                }
                m_NumNodes_tot    = m_NumElems+1;
            }
            break;
        case 3:
            {
                for(int i = 0; i < m_NumElems; i++)
                {
                    m_Conn[m_MaxElNod*i]   = 2*i+1;
                    m_Conn[m_MaxElNod*i+1] = 2*i+3;
                    m_Conn[m_MaxElNod*i+2] = 2*i+2;
                    m_NumNodes[i] = 3;
                }
                m_NumNodes_tot    = 2*m_NumElems+1;
            }
            break;
        default:
            ASSERTL0(false,"Number of Nodes in each element is not correct in SHARPy");
    }

    // Store element stiffness/mass (constant)
    NekDouble I,A,mu,E,G,rho;

    m_session->LoadParameter("BeamI", I);
    m_session->LoadParameter("BeamA", A);
    m_session->LoadParameter("BeamE", E);
    m_session->LoadParameter("BeamRho", rho);
    m_session->LoadParameter("BeamMu", mu);

    G = E/(2.0*(1.0+mu));

    NekDouble BeamMass[6][6];
    NekDouble BeamMassMatrix[6*m_NumElems][6];
    NekDouble BeamStiffness[6][6];
    NekDouble BeamStiffnessMatrix[6*m_NumElems][6];
    NekDouble BeamInvStiffnessMatrix[6*m_NumElems][6];

    for(int i = 0; i < 6; i++)
    {
        for(int j = 0; j < 6; j++)
        {
            BeamStiffness[i][j] = 0.0;
            BeamMass[i][j]      = 0.0;
        }
    }

    for(int i = 0; i < 6; i++)
    {
        for(int j = 0; j < 6; j++)
        {
            BeamMass[i][j]      = 0.0;
            BeamStiffness[i][j] = 0.0;
        }
    }

    BeamMass[0][0]      = rho*A;
    BeamMass[1][1]      = BeamMass[0][0];
    BeamMass[2][2]      = BeamMass[0][0];
    BeamMass[3][3]      = G*I*2.0;
    BeamMass[4][4]      = G*I;
    BeamMass[5][5]      = G*I;
    BeamStiffness[0][0] = E*A;
    BeamStiffness[1][1] = G*A;
    BeamStiffness[2][2] = G*A;
    BeamStiffness[3][3] = 2.0*G*I;
    BeamStiffness[4][4] = E*I;
    BeamStiffness[5][5] = E*I;
    
    ///*chaplin_test#3
    /*
    NekDouble d         = 56.e-3;
    A       = 0.25*M_PI*d*d;
    NekDouble m         = 4.0;
    NekDouble rho_water = 1e3;
    NekDouble m_water   = A*rho_water;
    NekDouble rho_cable = m/A;
    I       = rho_cable*0.25*M_PI*(0.5*d)*(0.5*d)*(0.5*d)*(0.5*d);
    NekDouble EI        = 1.3;
    NekDouble EA        = 790.0e3;
    NekDouble shear_factor = 0.4;
    BeamMass[0][0]      = m;
    BeamMass[1][1]      = m;
    BeamMass[2][2]      = m;
    BeamMass[3][3]      = 2.0*I;
    BeamMass[4][4]      = I;
    BeamMass[5][5]      = I;
    BeamStiffness[0][0] = EA;
    BeamStiffness[1][1] = shear_factor*EA;
    BeamStiffness[2][2] = shear_factor*EA;
    BeamStiffness[3][3] = 2.0*EI;
    BeamStiffness[4][4] = EI;
    BeamStiffness[5][5] = EI;
    */

    //*/ 
    int size = 6;
    Array<OneD, NekDouble> temp1(size*size,0.0);
    Array<OneD, NekDouble> temp2(size*size,0.0);

    for(int i = 0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
        {
            temp1[size*i+j] = BeamStiffness[i][j];
        }
    }

    SHARPy::Wrap_lu_invers(size, &temp1[0], &temp2[0]);

    for (int i = 0; i < m_NumElems; i++)
    {
        for(int j = 0; j < size; j++)
        {
            for(int k = 0; k < size; k++)
            {
                BeamStiffnessMatrix[i*size+j][k]     = BeamStiffness[j][k];
                BeamMassMatrix[i*size+j][k]          = BeamMass[j][k];
                BeamInvStiffnessMatrix[i*size+j][k]  = temp2[size*j+k];
            }
        }
    }

    for(int i = 0; i < m_NumElems*size; i++)
    {
        for(int j = 0; j < size; j++)
        {
            int cnt = i+m_NumElems*size*j;
            m_Stiff_Array[cnt]    = BeamStiffnessMatrix[i][j];
            m_Mass_Array[cnt]     = BeamMassMatrix[i][j];
            m_InvStiff_Array[cnt] = BeamInvStiffnessMatrix[i][j];
        }
    }

    // Define lumped masses at element nodes.
    for (int i = 0; i < m_NumElems; i++)
    {
         m_RBMass_Array[i] = 0.0;
    }
    // Element orientation.
    for (int i = 0; i < m_NumElems; i++)
    {
        m_Vector[3*i+1] = 1.0;
    }
    // Define element types.
    for (int i = 0; i < m_NumElems; i++)
    {
        m_MemNo[i] = 0;
    }
}

/**
 *
 **/
void Beam::Input_node(const LibUtilities::SessionReaderSharedPtr& m_session)
{
    //Initial position vector of grid points.
    for(int i = 0; i < m_NumNodes_tot; i++)
    {
        m_PosIni_Array[i] = m_BeamLength*i/(m_NumNodes_tot-1);
    }
    //Initial pretwist angle.
    Vmath::Zero(m_NumNodes_tot,m_PhiNodes,1);

    //Static point forces.
    Vmath::Zero(6*m_NumNodes_tot,m_ForceStatic_Array,1);

    //Boundary conditions.
    Vmath::Zero(m_NumNodes_tot,m_BoundConds,1);
    std::string BConds =
            m_session->GetSolverInfo("BeamBConds");

    if(BConds[0] == 'M')
    {
            m_BoundConds[0] = -1;
            m_BoundConds[m_NumNodes_tot-1] = -1;
            int MidNode = (m_NumNodes_tot-1)/2;
            if (BConds[1] == 'C') {m_BoundConds[MidNode] = 1;}
            if (BConds[1] == 'S') {m_BoundConds[MidNode] = 2;}
    }
    else
    {
        if(BConds[0] == 'C'){m_BoundConds[0] = 1;}
        if(BConds[0] == 'F'){m_BoundConds[0] =-1;}
        if(BConds[0] == 'S'){m_BoundConds[0] = 2;}
        if(BConds[0] == 'T'){m_BoundConds[0] = 3;}
        if(BConds[1] == 'C'){m_BoundConds[m_NumNodes_tot-1] = 1;}
        if(BConds[1] == 'F'){m_BoundConds[m_NumNodes_tot-1] =-1;}
        if(BConds[1] == 'S'){m_BoundConds[m_NumNodes_tot-1] = 2;}
        if(BConds[1] == 'T'){m_BoundConds[m_NumNodes_tot-1] = 3;}
    }
}

/**
 *
 **/
void Beam::AddGravityLoad(const LibUtilities::SessionReaderSharedPtr& m_session)
{
    NekDouble A, rho;

    m_session->LoadParameter("BeamA", A);
    m_session->LoadParameter("Beamrho", rho);

    NekDouble gravity       = 9.81;
    NekDouble MassPerLength = 4;
    NekDouble NodeSpacing   = m_BeamLength/(m_NumNodes_tot-1);  
    NekDouble ForcePerNode  = -1.0 * NodeSpacing * MassPerLength * gravity;

    for(int i = 1; i < m_NumNodes_tot-1; i++)
    {
        m_ForceStatic_Array[2*m_NumNodes_tot+i] += ForcePerNode;
    }

    m_ForceStatic_Array[2*m_NumNodes_tot]                   += 0.5 * ForcePerNode;
    m_ForceStatic_Array[2*m_NumNodes_tot+m_NumNodes_tot-1]  += 0.5 * ForcePerNode;
}

/**
 *  Calculates the rotation matrix from CRV.
 *  - Psi is the CRV that defines the rotation required to move a frame A over a
 *    frame B.
 *  - C will be the rotation matrix from A to B or, alternatively, the projection
 *    matrix from B to A
 *    @warning: for the same rotation, this function returns the transpose of
 *    Rot. i.e.
 *    RotCRV(Psi).T = Rot(psi2quat(Psi))
 *    @warning: equivalent to Psi2TransMat !!
 **/
void Beam::RotCRV(const Array<OneD, NekDouble> &Psi,
                                     Array<OneD, Array<OneD, NekDouble> > &C)
{
    double I[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
           C[i][j] = 0.0;
        }
    }

    Array<OneD, Array<OneD, NekDouble> > PsiSkew(3);
    PsiSkew[0] = Array<OneD, NekDouble>(3,0.0);
    PsiSkew[1] = Array<OneD, NekDouble>(3,0.0);
    PsiSkew[2] = Array<OneD, NekDouble>(3,0.0);
    Skew(Psi,PsiSkew);

    //norm
    NekDouble Psival = 0.0;
    Psival = sqrt(Psi[0]*Psi[0] + Psi[1]*Psi[1] +Psi[2]*Psi[2]);
    double dot[3][3];
    if(Psival > 1.0e-6)
    {
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                dot[i][j] = 0.0;
                for(int k = 0; k < 3; k++)
                {
                    dot[i][j] += PsiSkew[i][k]*PsiSkew[k][j];
                }
            }
        }
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                C[i][j] = I[i][j] + sin(Psival)/Psival*PsiSkew[i][j]
                    +(1.0-cos(Psival))/(Psival*Psival) * dot[i][j];
            }
        }
    }
    else
    {
        double Cadd[3][3];
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0 ; j < 3; j++)
            {
                Cadd[i][j]=I[i][j];
            }

        }

        double tol = 1e-8;
        double PsiPower[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

        int kk=0;
        NekDouble factorial = 1.0;
        NekDouble max_abs_Cadd = 0.0;
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                if(max_abs_Cadd < abs(Cadd[i][j]))
                {
                    max_abs_Cadd = abs(Cadd[i][j]);
                }
            }
        }
        while (max_abs_Cadd > tol)
        {
            //# add current
            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < 3; j++)
                {
                    C[i][j] += Cadd[i][j];
                }
            }

            //# new term
            kk += 1;
            double dot2[3][3];
            max_abs_Cadd = 0.0;
            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < 3; j++)
                {
                    dot2[i][j] = 0.0;
                    for(int k = 0; k < 3; k++)
                    {
                        dot2[i][j] += PsiPower[i][k]*PsiSkew[k][j];
                    }
                    factorial = factorial*kk;
                    Cadd[i][j] = 1.0/factorial*dot2[i][j];
                    if(max_abs_Cadd < abs(Cadd[i][j]))
                    {
                        max_abs_Cadd = abs(Cadd[i][j]);
                    }
                }
            }
        }
    }
}

/**
 *@brief Returns the skew-symmetric Matrix associated to Vector
 **/
void Beam::Skew(const Array<OneD, NekDouble> &Vector,
                                   Array<OneD, Array<OneD, NekDouble> > &SkewMat)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            SkewMat[i][j] = 0.0;
        }
    }
    SkewMat[0][1]   =-Vector[2];
    SkewMat[0][2]   = Vector[1];
    SkewMat[1][0]   = Vector[2];
    SkewMat[1][2]   =-Vector[0];
    SkewMat[2][0]   =-Vector[1];
    SkewMat[2][1]   = Vector[0];
}

}
}
