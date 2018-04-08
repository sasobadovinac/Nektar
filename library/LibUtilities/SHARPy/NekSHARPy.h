///////////////////////////////////////////////////////////////////////////////
//
// File NekSHARPy.h
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
// Description: Header file for the SHARPy class in Nektar++
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILIITIES_SHARPY_NEKSHARPY_H
#define NEKTAR_LIB_UTILIITIES_SHARPY_NEKSHARPY_H

#include <LibUtilities/SHARPy/NektarSHARPy.h>

namespace Nektar
{
    template <typename Dim, typename DataType>
    class Array;
    
	namespace LibUtilities
	{
		/**
		 * The NektarSHARPy class manages the use of the third-party code SHARPy to solve nonlinear structural dynamics.
		 * The function here defined will link to a proper implementation of the SHARPy routine.
		 * Depending on the user definition the functions can link to a class which is a wrapper around the SHARPy
		 * library or to a specific SHARPy implementation.
		 */
		class NekSHARPy;
		
		// A shared pointer to the NektarSHARPy object
		typedef std::shared_ptr<NekSHARPy>  NekSHARPySharedPtr;
		
		class NekSHARPy: public NektarSHARPy
		{
		public:
			/// Creates an instance of this class
            static NektarSHARPySharedPtr create()
            {
                return MemoryManager<NekSHARPy>::AllocateSharedPtr();
            }

            /// Name of class
            static std::string className;

			/// Initialises NektarFFT class members.
			NekSHARPy();
			
			// Distructor
			virtual  ~NekSHARPy();
			
		protected:
			
			virtual void v_InitialiseStatic(const LibUtilities::SessionReaderSharedPtr &pSession);

			virtual void v_InitialiseDynamic(const LibUtilities::SessionReaderSharedPtr& pSession, const Array<OneD, Array<OneD, NekDouble> >& CdCl);

			virtual void v_SolvenlnStatic(const LibUtilities::SessionReaderSharedPtr& pSession, const Array<OneD, Array<OneD, NekDouble> >& CdCl);

	        virtual void v_SolvenlnDynamic(const LibUtilities::SessionReaderSharedPtr& pSession,const Array<OneD, Array<OneD, NekDouble> > &CdCl, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &Vmotions);
			
		private:
		
            void Input_elem(const LibUtilities::SessionReaderSharedPtr& pSession);

            void Input_node(const LibUtilities::SessionReaderSharedPtr& pSession);

			void AddGravityLoad(const LibUtilities::SessionReaderSharedPtr& pSession, Array<OneD,NekDouble> &gForces);

			void RotCRV(const Array<OneD, NekDouble> &Psi, Array<OneD, Array<OneD, NekDouble> > &C);

			void Skew(const Array<OneD, NekDouble> &Vector, Array<OneD, Array<OneD, NekDouble> > &SkewMat);


			bool m_FollowerForce;        // Forces follow local deflections
            bool m_FollowerForceRig;     // Forces follow the body-fixed frame
            bool m_PrintInfo;            // Print informations about the solving process
            bool m_OutInBframe;          // Print velocities in local-frame (if not, use body-attached frame)
            bool m_OutInaframe;          // Print velocities in body-fixed-frame (if not, use inertial frame)

            int m_MaxElNod;              // Max number of nodes per element.
            int m_NElems;              	 // Elements discreising the beam
            int m_ElemProj;              // Element info computed in (1) global frame (2) fixed element frame (3) moving element frame.
            int m_MaxIterations;         // Newton-Raphson iterations for nonlinear solution
            int m_NumLoadSteps;          // Number of load increments
            int m_NumGauss;              // Gauss points in element integration
            int m_Solution;              // 902/912: cbeam3 linear/nonlinear flexible-body dynamic
            int m_NumNodesElem;
            int m_tNNodes;
            int m_NumDof;
            int m_DimMat;
            int m_ms;
            int m_cs;
            int m_ks;
            int m_fs;
            int m_as;
            int m_NumDof2;
            int m_NumSteps;
            int m_iStep;
            int m_NForcedDisp;

            NekDouble m_t0;
            NekDouble m_tf;
            NekDouble m_dt;
            NekDouble m_gamma;
            NekDouble m_beta;
            NekDouble m_DeltaCurved;
            NekDouble m_MinDelta;
            NekDouble m_NewmarkDamp;
            NekDouble m_BeamLength;
            NekDouble m_D;

            NekDouble m_U_inf;
            NekDouble m_rho_f;

            Array<OneD, int> m_NumNodes;
            Array<OneD, int> m_MemNo;
            Array<OneD, int> m_Conn;
            Array<OneD, int> m_Master_Elem;
            Array<OneD, int> m_BoundConds;
            Array<OneD, int> m_Master_Node;
            Array<OneD, int> m_Vdof;
            Array<OneD, int> m_Fdof;
            Array<OneD, int> m_Sdof;
            Array<OneD, int> m_ListIN;
            Array<OneD, int> m_NodeForcedDisp;

            Array<OneD, char>m_OutFile;

            Array<OneD, NekDouble> m_Length;
            Array<OneD, NekDouble> m_PreCurv;
            Array<OneD, NekDouble> m_Psi;
            Array<OneD, NekDouble> m_Vector;
            Array<OneD, NekDouble> m_Mass;
            Array<OneD, NekDouble> m_Stiff;
            Array<OneD, NekDouble> m_InvStiff;
            Array<OneD, NekDouble> m_RBMass;
            Array<OneD, NekDouble> m_PosIni;
            Array<OneD, NekDouble> m_ForceStatic;
            Array<OneD, NekDouble> m_PhiNodes;
            Array<OneD, NekDouble> m_PsiIni;
            Array<OneD, NekDouble> m_AppForces;
            Array<OneD, NekDouble> m_PosDefor;
            Array<OneD, NekDouble> m_PsiDefor;
            Array<OneD, NekDouble> m_PosDeforDot;
            Array<OneD, NekDouble> m_PsiDeforDot;
            Array<OneD, NekDouble> m_X;
            Array<OneD, NekDouble> m_DX;
            Array<OneD, NekDouble> m_DXdt;
            Array<OneD, NekDouble> m_DXddt;
            Array<OneD, NekDouble> m_Force;
            Array<OneD, NekDouble> m_Vrel;
            Array<OneD, NekDouble> m_VrelDot;
            Array<OneD, NekDouble> m_Mglobal;
            Array<OneD, NekDouble> m_Mvel;
            Array<OneD, NekDouble> m_Cglobal;
            Array<OneD, NekDouble> m_Cvel;
            Array<OneD, NekDouble> m_Kglobal;
            Array<OneD, NekDouble> m_Fglobal;
            Array<OneD, NekDouble> m_Qglobal;
            Array<OneD, NekDouble> m_PosDeforDDot;
            Array<OneD, NekDouble> m_PsiDeforDDot;
            Array<OneD, NekDouble> m_Quat;
            Array<OneD, NekDouble> m_Cao;
            Array<OneD, NekDouble> m_Asys;
            Array<OneD, NekDouble> m_F0;
            Array<OneD, NekDouble> m_Fa;
            Array<OneD, NekDouble> m_FTime;
            Array<OneD, NekDouble> m_Time;
            Array<OneD, NekDouble> m_ForcedVel;
            Array<OneD, NekDouble> m_ForcedVelDot;
            Array<OneD, NekDouble> m_PosForcedDisp;
            Array<OneD, NekDouble> m_PsiA_G;
		};
	}//end namespace LibUtilities
}//end of namespace Nektar
#endif //NEKTAR_LIB_UTILIITIES_SHARPY_NEKSHARPY_H
