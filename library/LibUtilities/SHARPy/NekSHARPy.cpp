///////////////////////////////////////////////////////////////////////////////
//
// File NekSHARPy.cpp
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
// Description: SHARPy class in Nektar++
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/SHARPy/NekSHARPy.h>
#include <LibUtilities/SHARPy/SHARPy.hpp>
#include <LocalRegions/MatrixKey.h>

namespace Nektar
{
	namespace LibUtilities
	{
		/**
		 * @class NektarSHARPy
         * The NektarSHARPy class manages the use of the third-party code SHARPy to solve nonlinear structural dynamics.
         * The function here defined will link to a proper implementation of the SHARPy routine.
         * Depending on the user definition the functions can link to a class which is a wrapper around the SHARPy
         * library or to a specific SHARPy implementation.
         */
		
		/**
		 * This constructor is protected as the objects of this class are never
		 * instantiated directly.
		 */
        string NekSHARPy::className
            = GetNektarSHARPyFactory().RegisterCreatorFunction("NekSHARPy",
                                                            NekSHARPy::create);

		NekSHARPy::NekSHARPy()
				:NektarSHARPy()
		{
			
		}
		
		NekSHARPy::~NekSHARPy()
		{
			
		}
		
		/**
		 * This allows initialisation of the class which cannot be completed
		 * during object construction (such as setting of initial conditions).
		 *
		 * Public interface routine to virtual function implementation.
		 */
	
		void NekSHARPy::v_InitialiseStatic(const LibUtilities::SessionReaderSharedPtr &pSession)
		{
		    //----------------------------Initial Displacement: ImpStart vs Static Solution
    		// Initialise static beam data.
    		// 1. Read input data.
    		m_OutFile = Array<OneD, char>(25,"BeamSol.txt");                  // The file to store the solution of Beam Dynamics
   			pSession->LoadParameter("BeamLength", m_BeamLength);             // Length of the beam
    		pSession->LoadParameter("BeamNumNodesElem", m_NumNodesElem,2);
    		pSession->LoadParameter("BeamMaxElNod", m_MaxElNod,3);           // Max number of nodes per element.
    		pSession->LoadParameter("BeamNumElems", m_NElems);          // Elements discreising the beam
    		pSession->LoadParameter("BeamNumForcedDisp", m_NForcedDisp,0);   //
    		pSession->LoadParameter("BeamElemProj", m_ElemProj,0);           // Element info computed in (1) global frame (2) fixed element frame (3) moving element frame.
    		pSession->LoadParameter("BeamMaxIterations", m_MaxIterations,20);// Newton-Raphson iterations for nonlinear solution
    		pSession->LoadParameter("BeamNumLoadSteps", m_NumLoadSteps,60);  // Number of load increments
    		pSession->LoadParameter("BeamNumGauss", m_NumGauss,3);           // Gauss points in element integration
    		pSession->LoadParameter("BeamSolution", m_Solution,112);         // 902/912: cbeam3 linear/nonlinear flexible-body dynamic 
    		pSession->LoadParameter("BeamDimMat", m_DimMat,24);
    		pSession->LoadParameter("BeamDeltaCurved", m_DeltaCurved,1.e-4); // Min. angle for two vectors to be parallel
    		pSession->LoadParameter("BeamMinDelta", m_MinDelta,1.e-6);       // Convergence param for Newton-Rhaphson iterations
    		pSession->LoadParameter("BeamDiameter", m_D,1.0);                // Diameter of the cylinder
    		pSession->LoadParameter("Velocity_Infinity", m_U_inf,1.0);       // Velocity of free stream
    		pSession->LoadParameter("Density_Fluid", m_rho_f,1.0);           // Density of the fluid

    		pSession->MatchSolverInfo("BeamFollowerForce","True",m_FollowerForce,false);
    		pSession->MatchSolverInfo("BeamFollowerForceRig","True",m_FollowerForceRig,false);
    		pSession->MatchSolverInfo("BeamPrintInfo","True",m_PrintInfo,false);
    		pSession->MatchSolverInfo("BeamOutInBframe","True",m_OutInBframe,false);
    		pSession->MatchSolverInfo("BeamOutInaframe","True",m_OutInaframe,false);

		    m_iStep             = 0;
    		m_NumNodes          = Array<OneD, int> (m_NElems);
    		m_MemNo             = Array<OneD, int> (m_NElems);
    		m_Conn              = Array<OneD, int> (m_MaxElNod*m_NElems);
    		m_Master_Elem       = Array<OneD, int> (m_MaxElNod*m_NElems*2);
    		m_Length            = Array<OneD, NekDouble> (m_NElems,0.0);
    		m_PreCurv           = Array<OneD, NekDouble> (3*m_NElems,0.0);
    		m_Psi               = Array<OneD, NekDouble> (3*m_NElems,0.0);
    		m_Vector            = Array<OneD, NekDouble> (3*m_NElems,0.0);
    		m_Mass              = Array<OneD, NekDouble> (6*m_NElems*6,0.0);
    		m_Stiff             = Array<OneD, NekDouble> (6*m_NElems*6,0.0);
    		m_InvStiff          = Array<OneD, NekDouble> (6*m_NElems*6,0.0);
    		m_RBMass            = Array<OneD, NekDouble> (6*m_MaxElNod*m_NElems*6,0.0);

			// Read element information of the beam
    		Input_elem(pSession);

    		m_BoundConds        = Array<OneD, int> (m_tNNodes);
    		m_Master_Node       = Array<OneD, int> (2*m_tNNodes,0);
    		m_Vdof              = Array<OneD, int> (m_tNNodes,0);
    		m_Fdof              = Array<OneD, int> (m_tNNodes,0);
    		m_Sdof              = Array<OneD, int> (m_tNNodes);
    		m_ListIN            = Array<OneD, int> (m_tNNodes);
    		m_NodeForcedDisp  = Array<OneD, int> (m_NForcedDisp);
    		m_PosForcedDisp   = Array<OneD, NekDouble> (3*m_NForcedDisp,0.0);
    		m_PosIni          = Array<OneD, NekDouble> (3*m_tNNodes,0.0);
    		m_PhiNodes        = Array<OneD, NekDouble> (m_tNNodes,0.0);
    		m_PsiIni          = Array<OneD, NekDouble> (3*m_NElems*m_MaxElNod,0.0);
    		m_PosDefor        = Array<OneD, NekDouble> (3*m_tNNodes,0.0);
    		m_PsiDefor        = Array<OneD, NekDouble> (3*m_NElems*m_MaxElNod,0.0);
    		m_PosDeforDot     = Array<OneD, NekDouble> (3*m_tNNodes,0.0);
    		m_PsiDeforDot     = Array<OneD, NekDouble> (3*m_NElems*m_MaxElNod,0.0);
    		m_PosDeforDDot    = Array<OneD, NekDouble> (3*m_tNNodes,0.0);
    		m_PsiDeforDDot    = Array<OneD, NekDouble> (3*m_NElems*m_MaxElNod,0.0);

    		m_Quat            = Array<OneD, NekDouble> (4,0.0);
    		m_PsiA_G          = Array<OneD, NekDouble> (3,0.0);
    		m_Cao             = Array<OneD, NekDouble> (3*3,0.0);

    		m_Quat[0]         = 1.0;
    		m_Cao[0]          = 1.0;
    		m_Cao[4]          = 1.0;
    		m_Cao[8]          = 1.0;

    		// Read node informatioin of the beam
    		Input_node(pSession);

    		// 2. Compute initial (undeformed) geometry.
    		SHARPy::Wrap_xbeam_undef_geom(
        		m_NElems,m_tNNodes,
        		&m_NumNodes[0],&m_MemNo[0],
        		&m_Conn[0],&m_Master_Elem[0],
        		&m_Length[0],&m_PreCurv[0],
        		&m_Psi[0],&m_Vector[0],&m_Mass[0],
        		&m_Stiff[0],&m_InvStiff[0],
        		&m_RBMass[0],&m_PosIni[0],
        		&m_PhiNodes[0],&m_PsiIni[0],
        		m_FollowerForce, m_FollowerForceRig,
        		m_PrintInfo, m_OutInBframe, m_OutInaframe,
        		m_ElemProj, m_MaxIterations, m_NumLoadSteps,
        		m_NumGauss, m_Solution, m_DeltaCurved,
        		m_MinDelta, m_NewmarkDamp);

    		// 3. Identify nodal degrees of freedom.
    		SHARPy::Wrap_xbeam_undef_dofs(
        		m_NElems,m_tNNodes,
        		&m_NumNodes[0],&m_MemNo[0],
        		&m_Conn[0],&m_Master_Elem[0],
        		&m_Length[0],&m_PreCurv[0],
        		&m_Psi[0],&m_Vector[0],&m_Mass[0],
        		&m_Stiff[0],&m_InvStiff[0],
        		&m_RBMass[0],&m_BoundConds[0],
        		&m_Master_Node[0],&m_Vdof[0],&m_Fdof[0],m_NumDof,
        		&m_Sdof[0]);
    		// Add RBMasses (Rigid-Body masses)
    		////

    		//Copy initial displacements.
    		Vmath::Vcopy(m_tNNodes*3, m_PosIni, 1, m_PosDefor, 1);
    		Vmath::Vcopy(m_NElems*m_MaxElNod*3, m_PsiIni, 1, m_PsiDefor, 1);


		}

		/**
 		 *
 		**/
		void NekSHARPy::v_InitialiseDynamic(const LibUtilities::SessionReaderSharedPtr& pSession,
        		const Array<OneD, Array<OneD, NekDouble> >& CdCl)
		{
    		// 1. Input data for transient dynamic solution.    
    		pSession->LoadParameter("BeamTime0", m_t0,0.0);
    		pSession->LoadParameter("BeamNewmarkDamp", m_NewmarkDamp,0.01);
    		pSession->LoadParameter("BeamTimefinal", m_tf,9999999);
    		pSession->LoadParameter("BeamDT", m_dt);

    		m_NumDof2       = m_NumDof*m_NumDof;
    		m_X             = Array<OneD, NekDouble> (m_NumDof,0.0);
    		m_DX            = Array<OneD, NekDouble> (m_NumDof,0.0);
    		m_DXdt          = Array<OneD, NekDouble> (m_NumDof,0.0);
    		m_DXddt         = Array<OneD, NekDouble> (m_NumDof,0.0);
    		m_Vrel          = Array<OneD, NekDouble> (6,0.0);
    		m_VrelDot       = Array<OneD, NekDouble> (6,0.0);
    		m_Mglobal = Array<OneD, NekDouble> (m_NumDof2,0.0);
    		m_Mvel    = Array<OneD, NekDouble> (m_NumDof*6,0.0);
    		m_Cglobal = Array<OneD, NekDouble> (m_NumDof2,0.0);
    		m_Cvel    = Array<OneD, NekDouble> (m_NumDof*6,0.0);
    		m_Kglobal = Array<OneD, NekDouble> (m_NumDof2,0.0);
    		m_Fglobal = Array<OneD, NekDouble> (m_NumDof2,0.0);
    		m_Qglobal = Array<OneD, NekDouble> (m_NumDof,0.0);
    		m_Asys    = Array<OneD, NekDouble>(m_NumDof2,0.0);

    		m_gamma         = 1.0/2.0+m_NewmarkDamp;
    		m_beta          = 1.0/4.0*(m_gamma+0.5)*(m_gamma+0.5);
    		m_NumSteps      = (m_tf - m_t0)/m_dt;

    		m_Time          = Array<OneD, NekDouble> (m_NumSteps);
    		for(int i = 0; i < m_NumSteps; i++)
    		{
        		m_Time[i] = m_t0 + m_dt*NekDouble(i);
    		}

    		//Initialize Local variables.
    		for(int k = 0; k < m_tNNodes; k++)
    		{
        		m_ListIN[k] = m_Vdof[k];
    		}

    		Array<OneD, NekDouble> Temp(m_tNNodes,0.0);

    		// store in an 1D array
    		Array<OneD, NekDouble> hforces(6*m_tNNodes);
    		Vmath::Vcopy(m_tNNodes,CdCl[0],1,Temp = hforces +   m_tNNodes,1);
    		Vmath::Vcopy(m_tNNodes,CdCl[1],1,Temp = hforces + 2*m_tNNodes,1);

    		// calculate the hydroforces applying on the body.
    		NekDouble amp = m_rho_f*m_U_inf*m_U_inf*m_D*m_BeamLength/(m_tNNodes-1);

    		Vmath::Smul(6*m_tNNodes,amp,hforces,1,hforces,1);

    		// The hydroforces applying at the two ends
    		for(int i = 0; i < 3; i++)
    		{
        		hforces[i*m_tNNodes]         = 0.5*hforces[i*m_tNNodes];
        		hforces[(i+1)*m_tNNodes-1]   = 0.5*hforces[(i+1)*m_tNNodes-1];
    		}

    		// Add Gravity Forces
    		bool IsGravity;
    		pSession->MatchSolverInfo("IsBeamGravity","True",IsGravity,false);
    		Array<OneD, NekDouble>gforces(6*m_tNNodes,0.0);
    		Array<OneD, NekDouble>dforces(6*m_tNNodes,0.0);

    		if(IsGravity)
    		{
        		AddGravityLoad(pSession,gforces);
        		Vmath::Vadd(6*m_tNNodes,gforces,1,hforces,1,dforces,1);
    		}

    		Vmath::Zero(m_NumDof2,m_Asys,1);
    		Vmath::Zero(m_NumDof2,m_Mglobal,1);
    		Vmath::Zero(m_NumDof2,m_Cglobal,1);
    		Vmath::Zero(m_NumDof2,m_Kglobal,1);
    		Vmath::Zero(m_NumDof2,m_Fglobal,1);
    		Vmath::Zero(m_NumDof, m_Qglobal,1);

    		//Extract initial displacements and velocities.
    		SHARPy::Wrap_cbeam3_solv_disp2state(
        		m_tNNodes,m_NumDof,m_NElems,
        		&m_Master_Node[0],&m_Vdof[0],&m_Fdof[0],
        		&m_PosDefor[0],&m_PsiDefor[0],
        		&m_PosDeforDot[0],&m_PsiDeforDot[0],
        		&m_X[0],&m_DXdt[0]);

    		//Compute initial acceleration (we are neglecting qdotdot in Kmass).
    		Vmath::Zero(m_tNNodes*3,m_PosDeforDDot,1);
    		Vmath::Zero(m_NElems*m_MaxElNod*3,m_PsiDeforDDot,1);
    		SHARPy::Wrap_cbeam3_asbly_dynamic(
        		m_NElems,m_tNNodes,&m_NumNodes[0],
        		&m_MemNo[0],&m_Conn[0],&m_Master_Elem[0],&m_Length[0],
        		&m_PreCurv[0],&m_Psi[0],&m_Vector[0],&m_Mass[0],
        		&m_Stiff[0],&m_InvStiff[0],&m_RBMass[0],
        		&m_Master_Node[0],&m_Vdof[0],&m_Fdof[0],&m_Sdof[0],
        		&m_PosIni[0],&m_PsiIni[0],
        		&m_PosDefor[0],&m_PsiDefor[0],
        		&m_PosDeforDot[0],&m_PsiDeforDot[0],
        		&m_PosDeforDDot[0],&m_PsiDeforDDot[0],
        		&dforces[0],&m_Vrel[0],&m_VrelDot[0],
        		m_NumDof,m_DimMat,
        		m_ms,&m_Mglobal[0],&m_Mvel[0],
        		m_cs,&m_Cglobal[0],&m_Cvel[0],
        		m_ks,&m_Kglobal[0],
        		m_fs,&m_Fglobal[0],&m_Qglobal[0],
        		m_FollowerForce, m_FollowerForceRig,
        		m_PrintInfo, m_OutInBframe, m_OutInaframe,
        		m_ElemProj, m_MaxIterations, m_NumLoadSteps,
        		m_NumGauss, m_Solution, m_DeltaCurved,
        		m_MinDelta, m_NewmarkDamp,
        		&m_Cao[0]);

    		Array<OneD, NekDouble> Temp0(m_NumDof,0.0);
    		Array<OneD, NekDouble> Temp1(m_NumDof,0.0);

    		SHARPy::Wrap_fem_m2v(
        		m_tNNodes,6,
        		&dforces[0],m_NumDof,
        		&Temp0[0],&m_ListIN[0]);

    		DNekMatSharedPtr fglobal
        		= MemoryManager<DNekMat>::AllocateSharedPtr(
            	m_NumDof,m_NumDof);

    		for(int i = 0; i < m_NumDof; i++)
    		{
        		for(int j = 0; j < m_NumDof; j++)
        		{
            		(*fglobal)(i,j)= m_Fglobal[i*m_NumDof+j];
        		}
    		}

    		Blas::Dgemv('N', m_NumDof, m_NumDof,
            	    1.0, &(fglobal->GetPtr())[0],
                	m_NumDof, &Temp0[0], 1,
                	0.0,  &Temp1[0],1);

    		Vmath::Vsub(m_NumDof,m_Qglobal,1,Temp1,1,m_Qglobal,1);
    		Vmath::Smul(m_NumDof,-1.0,m_Qglobal,1,m_Qglobal,1);
    		Vmath::Vadd(m_NumDof2,m_Mglobal,1,m_Asys,1,m_Asys,1);

    		int N = m_NumDof;
    		Array<OneD, int> ipivot (N);
    		int info = 0;
    		Lapack::Dgetrf(N, N, m_Asys.get(), N, ipivot.get(), info);
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
    		Vmath::Vcopy(m_NumDof, m_Qglobal, 1, m_DXddt, 1);
    		Lapack::Dgetrs('N', N, ncolumns_b, m_Asys.get(), N, ipivot.get(), m_DXddt.get(), N, info);
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
		void NekSHARPy::Input_elem(const LibUtilities::SessionReaderSharedPtr& pSession)
		{
    		// Connectivies.
    		Vmath::Zero(m_MaxElNod*m_NElems,m_Conn,1);
    		switch(m_NumNodesElem)
    		{
        		case 2:
            	{
                	for(int i = 0; i < m_NElems; i++)
                	{
                    	m_Conn[m_MaxElNod*i]   = i+1;
                    	m_Conn[m_MaxElNod*i+1] = i+2;
                    	m_NumNodes[i] = 2;
                	}
                	m_tNNodes    = m_NElems+1;
            	}
            	break;
        		case 3:
            	{
                	for(int i = 0; i < m_NElems; i++)
                	{
                    	m_Conn[m_MaxElNod*i]   = 2*i+1;
                    	m_Conn[m_MaxElNod*i+1] = 2*i+3;
                    	m_Conn[m_MaxElNod*i+2] = 2*i+2;
                    	m_NumNodes[i] = 3;
                	}
                	m_tNNodes    = 2*m_NElems+1;
            	}
            	break;
        		default:
            		ASSERTL0(false,"Number of Nodes in each element is not correct in SHARPy");
    		}

    		// Store element stiffness/mass (constant)
    		NekDouble I,A,mu,E,G,rho;

    		pSession->LoadParameter("BeamI", I);
    		pSession->LoadParameter("BeamA", A);
    		pSession->LoadParameter("BeamE", E);
    		pSession->LoadParameter("BeamRho", rho);
    		pSession->LoadParameter("BeamMu", mu);

    		G = E/(2.0*(1.0+mu));

    		NekDouble BeamMass[6][6];
    		NekDouble BeamMassMatrix[6*m_NElems][6];
    		NekDouble BeamStiffness[6][6];
    		NekDouble BeamStiffnessMatrix[6*m_NElems][6];
    		NekDouble BeamInvStiffnessMatrix[6*m_NElems][6];

    		for(int i = 0; i < 6; i++)
    		{
        		for(int j = 0; j < 6; j++)
        		{
            		BeamStiffness[i][j] = 0.0;
            		BeamMass[i][j]      = 0.0;
        		}
    		}

    		for(int i = 0; i < 6*m_NElems; i++)
    		{
        		for(int j = 0; j < 6; j++)
        		{
            		BeamMassMatrix[i][j]        = 0.0;
            		BeamStiffnessMatrix[i][j]   = 0.0;
            		BeamInvStiffnessMatrix[i][j]= 0.0;
        		}
    		}
	
    		BeamMass[0][0]      = rho*A;
    		BeamMass[1][1]      = BeamMass[0][0];
    		BeamMass[2][2]      = BeamMass[0][0];
    		BeamMass[3][3]      = rho*I;
    		BeamMass[4][4]      = rho*I;
    		BeamMass[5][5]      = 2.0*rho*I; //???
    		BeamStiffness[0][0] = E*A;
    		BeamStiffness[1][1] = G*A;
    		BeamStiffness[2][2] = G*A;
    		BeamStiffness[3][3] = 2.0*G*I;
    		BeamStiffness[4][4] = E*I;
    		BeamStiffness[5][5] = E*I;

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

    		for (int i = 0; i < m_NElems; i++)
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

    		for(int i = 0; i < m_NElems*size; i++)
    		{
        		for(int j = 0; j < size; j++)
        		{
            		int cnt = i+m_NElems*size*j;
            		m_Stiff[cnt]    = BeamStiffnessMatrix[i][j];
            		m_Mass[cnt]     = BeamMassMatrix[i][j];
            		m_InvStiff[cnt] = BeamInvStiffnessMatrix[i][j];
       			}
    		}

    		// Define lumped masses at element nodes.
    		for (int i = 0; i < m_NElems; i++)
    		{
         		m_RBMass[i] = 0.0;
  			}
   		 	// Element orientation.
    		for (int i = 0; i < m_NElems; i++)
   		 	{
        		m_Vector[3*i+1] = 1.0;
   		 	}
    		// Define element types.
    		for (int i = 0; i < m_NElems; i++)
    		{
        		m_MemNo[i] = 0;
    		}
		}


		/**
		 *
 		**/
		void NekSHARPy::v_SolvenlnStatic(const LibUtilities::SessionReaderSharedPtr& pSession,
        	const Array<OneD, Array<OneD, NekDouble> >& CdCl)
		{
    		// Store the forces in x and y directions in an 1D array.
    		Array<OneD, NekDouble> hforces(6*m_tNNodes,0.0);
    		Array<OneD, NekDouble> Temp(m_tNNodes,0.0);
    		Vmath::Vcopy(m_tNNodes,CdCl[0],1,Temp = hforces +   m_tNNodes,1);
    		Vmath::Vcopy(m_tNNodes,CdCl[1],1,Temp = hforces + 2*m_tNNodes,1);

    		// Dimensionalize the hydroforces.
    		NekDouble amp = m_rho_f*m_U_inf*m_U_inf*m_D*m_BeamLength/(m_tNNodes-1);
    		Vmath::Smul(6*m_tNNodes,amp,hforces,1,hforces,1);

    		for(int i = 0; i < 3; i++)
    		{
        		hforces[i*m_tNNodes]         = 0.5*hforces[i*m_tNNodes];
        		hforces[(i+1)*m_tNNodes-1]   = 0.5*hforces[(i+1)*m_tNNodes-1];
    		}

    		//Add Gravity Forces
    		bool IsGravity;
    		pSession->MatchSolverInfo("IsBeamGravity","True",IsGravity,false);
    		Array<OneD, NekDouble>  gforces(6*m_tNNodes,0.0);
    		Array<OneD, NekDouble>  sforces(6*m_tNNodes,0.0);
    		if(IsGravity)
    		{
        		AddGravityLoad(pSession,gforces);
        		Vmath::Vadd(6*m_tNNodes,gforces,1,hforces,1,sforces,1);
    		}


    		SHARPy::Wrap_cbeam3_solv_nlnstatic (
        		m_NumDof,m_NElems,
        		&m_NumNodes[0], &m_MemNo[0], &m_Conn[0],
        		&m_Master_Elem[0],
        		&m_Length[0], &m_PreCurv[0],
        		&m_Psi[0], &m_Vector[0], &m_Mass[0],
        		&m_Stiff[0],
        		&m_InvStiff[0], &m_RBMass[0],
        		m_tNNodes,
        		&m_Master_Node[0], &m_Vdof[0], &m_Fdof[0],&m_Sdof[0],
        		&sforces[0],
        		&m_PosIni[0],&m_PsiIni[0],
        		&m_PosDefor[0],&m_PsiDefor[0],
        		m_FollowerForce, m_FollowerForceRig,
        		m_PrintInfo, m_OutInBframe, m_OutInaframe,
        		m_ElemProj, m_MaxIterations, m_NumLoadSteps,
        		m_NumGauss, m_Solution, m_DeltaCurved,
        		m_MinDelta, m_NewmarkDamp,
        		m_NForcedDisp,&m_NodeForcedDisp[0],
        		&m_PosForcedDisp[0],&m_PsiA_G[0]);
		}

		/**
 	     *
 		**/
		void NekSHARPy::v_SolvenlnDynamic(const LibUtilities::SessionReaderSharedPtr& pSession,
        	        const Array<OneD, Array<OneD, NekDouble> > &CdCl,
            	    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &Vmotions)
		{
    		//
    		Array<OneD, NekDouble> vel0(m_tNNodes,0.0);
    		Array<OneD, NekDouble> vel1(m_tNNodes,0.0);
    		Vmath::Vcopy(m_tNNodes,Vmotions[0][1],1,vel0,1);
    		Vmath::Vcopy(m_tNNodes,Vmotions[1][1],1,vel1,1);

    		Array<OneD, NekDouble> Temp(m_tNNodes,0.0);
    		Array<OneD, NekDouble> hforces(6*m_tNNodes,0.0);
    		Vmath::Vcopy(m_tNNodes,CdCl[0],1,Temp = hforces +   m_tNNodes,1);
    		Vmath::Vcopy(m_tNNodes,CdCl[1],1,Temp = hforces + 2*m_tNNodes,1);

    		// dimensionalize the hydroforces.
    		NekDouble amp = m_rho_f*m_U_inf*m_U_inf*m_D*m_BeamLength/(m_tNNodes-1);
    		Vmath::Smul(6*m_tNNodes,amp,hforces,1,hforces,1);
    		for(int i = 0; i < 3; i++)
    		{
        		hforces[i*m_tNNodes]         = 0.5*hforces[i*m_tNNodes];
        		hforces[(i+1)*m_tNNodes-1]   = 0.5*hforces[(i+1)*m_tNNodes-1];
    		}

    		//Add Gravity Forces
    		bool IsGravity;
    		pSession->MatchSolverInfo("IsBeamGravity","True",IsGravity,false);
    		Array<OneD, NekDouble>gforces(6*m_tNNodes,0.0);
    		Array<OneD, NekDouble>dforces(6*m_tNNodes,0.0);
    		if(IsGravity)
    		{
        		AddGravityLoad(pSession,gforces);
        		Vmath::Vadd(6*m_tNNodes,gforces,1,hforces,1,dforces,1);
    		}

    		int size = 4;
    		Array<OneD, NekDouble> Temp0(size*size,0.0);
    		Array<OneD, NekDouble> Temp1(size*size,0.0);
    		Array<OneD, NekDouble> Temp2(size*size,0.0);
    		Array<OneD, NekDouble> Temp3(size*size,0.0);
    		Array<OneD, NekDouble> Unit4(size*size,0.0);

    		for(int i = 0; i < 4; i++)
    		{
        		Unit4[i*4+i] = 1.0;
    		}

    		NekDouble tmp = 0.25*m_dt;

    		// Update transformation matrix for given angular velocity
    		SHARPy::Wrap_xbeam_QuadSkew(4, &m_Vrel[3], &Temp0[0]);
    		Vmath::Svtvp(size*size, tmp, Temp0, 1, Unit4, 1, Temp1, 1);
    		SHARPy::Wrap_lu_invers(size, &Temp1[0], &Temp2[0]);
    		Vmath::Svtvp(size, -1.0*tmp, Temp0, 1, Unit4, 1,  Temp1, 1);
    		SHARPy::Wrap_matmul(size, size, size, &Temp1[0], &m_Quat[0], &Temp3[0]);
    		SHARPy::Wrap_matmul(size, size, size, &Temp2[0], &Temp3[0], &m_Quat[0]);
    		SHARPy::Wrap_xbeam_Rot(&m_Quat[0], &m_Cao[0]);

    		// Predictor step.
    		Array<OneD,NekDouble> Temp_(m_NumDof, 0.0);
    		Vmath::Smul(m_NumDof, (0.5-m_beta)*m_dt*m_dt, m_DXddt, 1, Temp_, 1);
    		Vmath::Svtvp(m_NumDof, m_dt, m_DXdt, 1, Temp_, 1, Temp_, 1);
    		Vmath::Vadd(m_NumDof, m_X, 1, Temp_, 1, m_X, 1);
    		Vmath::Svtvp(m_NumDof, (1.0-m_gamma)*m_dt, m_DXddt, 1, m_DXdt, 1, m_DXdt, 1);
    		Vmath::Zero(m_NumDof, m_DXddt, 1);

    		int Iter;
    		// Iteration until convergence.
    		for(Iter = 0; Iter < m_MaxIterations+1; Iter++)
    		{
        		if (Iter == m_MaxIterations)
        		{
            		ASSERTL0(false,"Solution did not converge in SHARPy");
        		}

        		// Update nodal positions and velocities .
        		SHARPy::Wrap_cbeam3_solv_state2disp(
            		m_NElems,m_tNNodes,&m_NumNodes[0],
            		&m_MemNo[0],&m_Conn[0],&m_Master_Elem[0],
            		&m_Length[0],&m_PreCurv[0],&m_Psi[0],
            		&m_Vector[0],&m_Mass[0],&m_Stiff[0],
            		&m_InvStiff[0],&m_RBMass[0],
            		&m_Master_Node[0],&m_Vdof[0],&m_Fdof[0],
            		&m_PosIni[0],&m_PsiIni[0],
            		m_NumDof,&m_X[0],&m_DXdt[0],&m_PosDefor[0],
            		&m_PsiDefor[0],&m_PosDeforDot[0],
            		&m_PsiDeforDot[0]);

        		// Compute system functionals and matrices. (Use initial accelerations for Kgyr).
        		Vmath::Zero(m_NumDof, m_Qglobal, 1);
        		Vmath::Zero(m_NumDof*6, m_Mvel, 1);
        		Vmath::Zero(m_NumDof*6, m_Cvel, 1);
        		Vmath::Zero(m_NumDof2,m_Mglobal,1);
        		Vmath::Zero(m_NumDof2,m_Cglobal,1);
        		Vmath::Zero(m_NumDof2,m_Kglobal,1);
        		Vmath::Zero(m_NumDof2,m_Fglobal,1);
        		Vmath::Smul(m_tNNodes*3,0.0,m_PosDefor,1,m_PosDeforDDot,1);
        		Vmath::Smul(m_NElems*m_MaxElNod*3,0.0,m_PsiDefor,1,m_PsiDeforDDot,1);

        		SHARPy::Wrap_cbeam3_asbly_dynamic(
            		m_NElems,m_tNNodes,&m_NumNodes[0],
            		&m_MemNo[0],&m_Conn[0],&m_Master_Elem[0],&m_Length[0],
            		&m_PreCurv[0],&m_Psi[0],&m_Vector[0],&m_Mass[0],
            		&m_Stiff[0],&m_InvStiff[0],&m_RBMass[0],
            		&m_Master_Node[0],&m_Vdof[0],&m_Fdof[0],&m_Sdof[0],
            		&m_PosIni[0],&m_PsiIni[0],
            		&m_PosDefor[0],&m_PsiDefor[0],
            		&m_PosDeforDot[0],&m_PsiDeforDot[0],
            		&m_PosDeforDDot[0],&m_PsiDeforDDot[0],
            		&dforces[0],&m_Vrel[0],&m_VrelDot[0],m_NumDof,m_DimMat,
            		m_ms,&m_Mglobal[0],&m_Mvel[0],
            		m_cs,&m_Cglobal[0],&m_Cvel[0],
            		m_ks,&m_Kglobal[0],
            		m_fs,&m_Fglobal[0],&m_Qglobal[0],
            		m_FollowerForce, m_FollowerForceRig,
            		m_PrintInfo, m_OutInBframe, m_OutInaframe,
            		m_ElemProj, m_MaxIterations, m_NumLoadSteps,
            		m_NumGauss, m_Solution, m_DeltaCurved,
            		m_MinDelta, m_NewmarkDamp,
            		&m_Cao[0]);

        		// Compute admissible error.
        		Array<OneD, NekDouble>tmp0(m_NumDof,0.0);
        		Array<OneD, NekDouble>tmp1(m_NumDof,0.0);
        		Vmath::Vabs(m_NumDof, m_Qglobal, 1, tmp0, 1);
        		NekDouble MinDelta;
        		MinDelta = m_MinDelta*std::max(1.0, Vmath::Vmax(m_NumDof,tmp0,1));

        		// Compute the residual.
        		Array<OneD,NekDouble> Temp4(m_NumDof,0.0);
        		Array<OneD,NekDouble> Temp5(m_NumDof,0.0);
        		Array<OneD,NekDouble> Temp6(m_NumDof,0.0);
        		Array<OneD,NekDouble> Temp7(m_NumDof,0.0);

        		DNekMatSharedPtr mglobal =
            		MemoryManager<DNekMat>::AllocateSharedPtr(
                	m_NumDof,m_NumDof);

        		for(int i = 0; i < m_NumDof; i++)
        		{
            		for(int j = 0; j < m_NumDof; j++)
            		{
                		(*mglobal)(i,j) = m_Mglobal[i*m_NumDof+j];
            		}
        		}
        		Blas::Dgemv('N', m_NumDof, m_NumDof,
                	    1.0, &(mglobal->GetPtr())[0],
                    	m_NumDof, &m_DXddt[0], 1,
                    	0.0, &Temp4[0], 1);

        		SHARPy::Wrap_matmul(
            		m_NumDof,6,1,
            		&m_Mvel[0], &m_VrelDot[0],
            		&Temp5[0]);

        		Vmath::Vadd(m_NumDof, Temp4, 1, Temp5, 1, Temp5, 1);
        		Vmath::Vadd(m_NumDof, Temp5, 1, m_Qglobal, 1, m_Qglobal, 1);

        		SHARPy::Wrap_fem_m2v(
            		m_tNNodes,6,
            		&dforces[0],
            		m_NumDof, &Temp6[0],
            		&m_ListIN[0]);

        		DNekMatSharedPtr fglobal =
            		MemoryManager<DNekMat>::AllocateSharedPtr(
                		m_NumDof,m_NumDof);

        		for(int i = 0; i < m_NumDof; i++)
        		{
            		for(int j = 0; j < m_NumDof; j++)
            		{
                		(*fglobal)(i,j)= m_Fglobal[i*m_NumDof+j];
            		}
        		}

        		Blas::Dgemv('N', m_NumDof, m_NumDof,
            		        1.0, &(fglobal->GetPtr())[0],
                		    m_NumDof, &Temp6[0], 1,
                    		0.0, &Temp7[0], 1);

        		Vmath::Vsub(m_NumDof, m_Qglobal, 1, Temp7, 1, m_Qglobal, 1);

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

        		Vmath::Zero(m_NumDof2,m_Asys,1);
        		Vmath::Vadd(m_NumDof2,m_Kglobal,1,m_Asys,1,m_Asys,1);
        		Vmath::Svtvp(m_NumDof2,factor1,m_Cglobal,1,m_Asys,1,m_Asys,1);
        		Vmath::Svtvp(m_NumDof2,factor2,m_Mglobal,1,m_Asys,1,m_Asys,1);

        		// Calculation of the correction.
        		Vmath::Smul(m_NumDof, -1.0, m_Qglobal, 1, m_Qglobal, 1);

        		Array<OneD, int> ipivot (m_NumDof);
        		int info = 0;
        		Lapack::Dgetrf(m_NumDof, m_NumDof, m_Asys.get(), m_NumDof, ipivot.get(), info);
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
        			Lapack::Dgetrs( 'N', m_NumDof, ncolumns_b, m_Asys.get(),
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
        		m_NElems,m_tNNodes,&m_NumNodes[0],
        		&m_MemNo[0],&m_Conn[0],&m_Master_Elem[0],
        		&m_Length[0],&m_PreCurv[0],&m_Psi[0],
        		&m_Vector[0],&m_Mass[0],
        		&m_Stiff[0],&m_InvStiff[0],
        		&m_RBMass[0],&m_Master_Node[0],&m_Vdof[0],&m_Fdof[0],
        		&m_PosIni[0],&m_PsiIni[0],
        		m_NumDof,&m_X[0],&m_DXdt[0],
        		&m_PosDefor[0],&m_PsiDefor[0],
        		&m_PosDeforDot[0],&m_PsiDeforDot[0]);

    		//copy the displacement and velocity for current step
    		Array<OneD, NekDouble> temp(3*m_tNNodes,0.0);

    		Vmath::Vsub(3*m_tNNodes,  m_PosDefor, 1,m_PosIni,1,temp,1);

    		for(int i = 0; i < m_tNNodes; i++)
    		{
        		Vmotions[0][0][i] =  temp[  m_tNNodes + i];
        		Vmotions[1][0][i] =  temp[2*m_tNNodes + i];
        		Vmotions[0][1][i] =  m_PosDeforDot[  m_tNNodes + i];
        		Vmotions[1][1][i] =  m_PosDeforDot[2*m_tNNodes + i];
    		}

    		//scale the displacement and velocity for output file
    		Vmath::Smul(m_tNNodes,1.0/m_D,Vmotions[0][0],1,Vmotions[0][0],1);
    		Vmath::Smul(m_tNNodes,1.0/m_D,Vmotions[1][0],1,Vmotions[1][0],1);
    		Vmath::Smul(m_tNNodes,1.0/m_U_inf,Vmotions[0][1],1,Vmotions[0][1],1);
    		Vmath::Smul(m_tNNodes,1.0/m_U_inf,Vmotions[1][1],1,Vmotions[1][1],1);

    		//calculate accel for output file.
    		Vmath::Vsub(m_tNNodes,Vmotions[0][1],1,vel0,1,vel0,1);
    		Vmath::Vsub(m_tNNodes,Vmotions[1][1],1,vel1,1,vel1,1);
    		Vmath::Smul(m_tNNodes,1.0/m_dt,vel0,1,Vmotions[0][2],1);
    		Vmath::Smul(m_tNNodes,1.0/m_dt,vel1,1,Vmotions[1][2],1);

    		m_iStep++;
		}


		/**
		 *
		**/
        void NekSHARPy::Input_node(const LibUtilities::SessionReaderSharedPtr& pSession)
		{
			//Initial position vector of grid points.
   			for(int i = 0; i < m_tNNodes; i++)
    		{
        		m_PosIni[i] = m_BeamLength*i/(m_tNNodes-1);
 		    }
    		//Initial pretwist angle.
    		Vmath::Zero(m_tNNodes,m_PhiNodes,1);

    		//Boundary conditions.
    		Vmath::Zero(m_tNNodes,m_BoundConds,1);
    		std::string BConds =
            	pSession->GetSolverInfo("BeamBConds");

    		if(BConds[0] == 'M')
    		{
            	m_BoundConds[0] = -1;
            	m_BoundConds[m_tNNodes-1] = -1;
            	int MidNode = (m_tNNodes-1)/2;
            	if (BConds[1] == 'C') {m_BoundConds[MidNode] = 1;}
            	if (BConds[1] == 'S') {m_BoundConds[MidNode] = 2;}
    		}
    		else
    		{
        		if(BConds[0] == 'C'){m_BoundConds[0] = 1;}
        		if(BConds[0] == 'F'){m_BoundConds[0] =-1;}
        		if(BConds[0] == 'S'){m_BoundConds[0] = 2;}
        		if(BConds[0] == 'T'){m_BoundConds[0] = 3;}
        		if(BConds[1] == 'C'){m_BoundConds[m_tNNodes-1] = 1;}
        		if(BConds[1] == 'F'){m_BoundConds[m_tNNodes-1] =-1;}
        		if(BConds[1] == 'S'){m_BoundConds[m_tNNodes-1] = 2;}
        		if(BConds[1] == 'T'){m_BoundConds[m_tNNodes-1] = 3;}
    		}

		}

		/**
		 *
 		**/
		void NekSHARPy::AddGravityLoad(const LibUtilities::SessionReaderSharedPtr& pSession, Array<OneD,NekDouble> &gForces)
		{
    		NekDouble MPerLength, g;
    		pSession->LoadParameter("BeamRho",MPerLength);
    		pSession->LoadParameter("Gravity", g, 9.81);
    		pSession->LoadParameter("Gravity", g, 9.81);

    		NekDouble NodeSpacing   = m_BeamLength/(m_tNNodes-1);
    		NekDouble a3  = -1.0 * NodeSpacing * MPerLength * g;
    		NekDouble ForcePerNode[3] = {0.0,0.0,a3};

	
    		// Obtain transformation from Earth to a-fram.
    		// if PsiA_G == None: CGa = Psi2TransMat(PsiA_G)
    		// else CGa = Psi2TransMat(PsiA_G)
    		NekDouble CaG[3][3];
    		Array<OneD, NekDouble> PsiA_G(3,0.0);
    		Array<OneD, Array<OneD, NekDouble> > CGa(3);
    		CGa[0]= Array<OneD, NekDouble>(3,0.0);
    		CGa[1]= Array<OneD, NekDouble>(3,0.0);
    		CGa[2]= Array<OneD, NekDouble>(3,0.0);

    		if(m_FollowerForceRig)
   		 	{
        		// rotate gravity loads
        		Vmath::Vcopy(3, m_PsiA_G, 1, PsiA_G, 1);
    		}
   		 	else
    		{
        		Vmath::Zero(3,PsiA_G,1);
    		}

    		RotCRV(PsiA_G,CGa);

    		//CaG = CGa.T
    		for(int i = 0; i < 3; i++)
    		{
        		for(int j = 0; j < 3; j++)
        		{
            		CaG[i][j] = CGa[j][i];
        		}
    		}	
    		//  Blas::Dgemm('N', 'N', m_derivMat[i]->GetRows(), m_numElmt,
    		//  m_derivMat[i]->GetColumns(), 1.0,
    		//  m_derivMat[i]->GetRawPtr(),
    		//  m_derivMat[i]->GetRows(), input.get(), nPhys,
    		//  0.0, &Diff[i][0],nPhys);

    		// Force in a-frame
    		NekDouble Force_a[3];
    		for(int i = 0; i < 3; i++)
    		{
        		Force_a[i] = 0.0;
        		for(int k = 0; k < 3; k++)
        		{
            		Force_a[i] += CaG[i][k]*ForcePerNode[k];
        		}
    		}

    		// Indices for boundary nodes.
    		int Root = 0;
    		int Tip = m_tNNodes -1;

   			// Apply forces. 
    		Array<OneD, Array<OneD, NekDouble>> BeamForces(3);
    		BeamForces[0] =  Array<OneD, NekDouble>(m_tNNodes,0.0);
    		BeamForces[1] =  Array<OneD, NekDouble>(m_tNNodes,0.0);
    		BeamForces[2] =  Array<OneD, NekDouble>(m_tNNodes,0.0);
    		for(int i = 0; i < 3; i++)
    		{
        		for(int j = Root+1; j < m_tNNodes-1; j++)
        		{
            		BeamForces[i][j] = Force_a[i];
        		}
        		BeamForces[i][Root]  = 0.5*Force_a[i];
        		BeamForces[i][Tip]   = 0.5*Force_a[i];
    		}
    		// Add RB masses ***
    		// Loop through nodes to get moment arm at each.***
    		for(int i = 0; i < m_tNNodes; i++)
    		{
        		gForces[0*m_tNNodes+i] = BeamForces[0][i];
        		gForces[1*m_tNNodes+i] = BeamForces[1][i];
       			gForces[2*m_tNNodes+i] = BeamForces[2][i];
    		}
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
		void NekSHARPy::RotCRV(const Array<OneD, NekDouble> &Psi,
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
                    		PsiPower[i][j] = dot2[i][j];
                    		factorial = factorial*kk;
                    		Cadd[i][j] = 1.0/factorial*PsiPower[i][j];
                		}
            		}
            		max_abs_Cadd = 0.0;
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
        		}
    		}
		}


		/**
 		*@brief Returns the skew-symmetric Matrix associated to Vector
 		**/
		void NekSHARPy::Skew(const Array<OneD, NekDouble> &Vector,
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
		
	}//end namespace LibUtilities
}//end of namespace Nektar
