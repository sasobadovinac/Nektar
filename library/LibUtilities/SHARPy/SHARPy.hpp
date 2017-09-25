///////////////////////////////////////////////////////////////////////////////
//
// File SHARPy.hpp
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
// Description: wrapper of functions around SHARPy routines.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_SHARPY_HPP
#define NEKTAR_LIB_UTILITIES_SHARPY_HPP


// Translations for using Fortran version of SHARPy
namespace SHARPy
{
    extern "C"
    {
		//SHARPY ROUTINES
        void __test_MOD_wrap_input_setup (const int& NumElems,       char*  OutFile,
                                          const bool& FollowerForce, const bool& FollowerForceRig,
                                          const bool& PrintInfo,     const bool& OutInBframe,
                                          const bool& OutInaframe,   const int& ElemProj,
                                          const int& MaxIterations,  const int& NumLoadSteps,
                                          const int& NumGauss,       const int& Solution,
                                          const double& DeltaCurved, const double &MinDelta,
                                          const double& NewmarkDamp);

		//
		void __test_MOD_wrap_input_elem  (const int& NumElems,		int& NumNodes_tot,
										  int* NumNodes,			int* MemNo,
										  int* Conn,				int* Master_Array,
										  double* Length,			double* PreCurv,
										  double* Psi,				double* Vector,
										  double* Mass_Array,		double* Stiff_Array,
										  double* InvStiff_Array, 	double* RBMass_Array);

		//
		void __test_MOD_wrap_input_node  (const int& NumElems,            const int& NumNodes_tot,
                                          const int* NumNodes,            const int* MemNo,
                                          const int* Conn,                const int* Master_Array,
                                          const double* Length,           const double* PreCurv,
                                          const double* Psi,              const double* Vector,
                                          const double* Mass_Array,       const double* Stiff_Array,
                                          const double* InvStiff_Array,   const double* RBMass_Array,
                                          int* BoundConds,	              double* PosIni_Array,
                                          double* ForceStatic_Array,      double*PhiNodes);

		//
		void __test_MOD_wrap_xbeam_undef_geom (const int& NumElems,          int& NumNodes_tot,
                                               int* NumNodes,                int* MemNo,
                                               int* Conn,                    int* Master_Array,
                                               double* Length,               double* PreCurv,
                                               double* Psi,                  double* Vector,
                                               double* Mass_Array,           double* Stiff_Array,
                                               double* InvStiff_Array,       double* RBMass_Array,
                                               const double* PosIni_Array,   const double* PhiNodes,
                                               double* PsiIni_Array,         const bool& FollowerForce,
                                               const bool& FollowerForceRig, const bool& PrintInfo,
                                               const bool& OutInBframe,      const bool& OutInaframe,
                                               const int& ElemProj,          const int& MaxIterations,
                                               const int& NumLoadSteps,      const int& NumGauss,
                                               const int& Solution,          const double& DeltaCurved,
                                               const double& MinDelta,       const double& NewmarkDamp);
 
		//
		void __test_MOD_wrap_xbeam_undef_dofs (const int& NumElems,              const int& NumNodes_tot,
                                               const int* NumNodes,              const int* MemNo,
                                               const int* Conn,                  const int* Master_Array,
                                               const double* Length,             const double* PreCurv,
                                               const double* Psi,                const double* Vector,
                                               const double* Mass_Array,         const double* Stiff_Array,
                                               const double* InvStiff_Array,     const double* RBMass_Array,
                                               const int* BoundConds,            int* Nod_Master,
                                               int* Nod_Vdof,                    int* Nod_Fdof,
                                               int& NumDof,                      int* Nod_Sflag);

		//
		void __test_MOD_wrap_cbeam3_solv_nlnstatic (const int& NumDof, 				const int& NumElems, 
												    const int* NumNodes, 			const int* MemNo, 
												    const int* Conn,        		const int* Master_Array,    
                       							    const double* Length, 			const double* PreCurv,    
                    							    const double* Psi, 				const double* Vector, 
												    const double* Mass_Array,   	const double* Stiff_Array,      
                      							    const double* InvStiff_Array, 	const double* RBMass_Array,           
												    const int& NumNodes_tot, 		const int* Master, 
												    const int* Vdof, 				const int* Fdof,       //for pack_xbnode
                   								    const double* AppForces_Array,	const double*  Coords_Array,
												    const double* Psi0_Array,		double* PosDefor_Array, 
												    double* PsiDefor_Array,		    const bool& FollowerForce, 
												    const bool& FollowerForceRig,   const bool& PrintInfo, 
												    const bool& OutInBframe, 	    const bool& OutInaframe,    
												    const int& ElemProj, 			const int& MaxIterations, 
												    const int& NumLoadSteps,  	    const int& NumGauss, 
												    const int& Solution, 			const double& DeltaCurved,	
												    const double& MinDelta, 		const double& NewmarkDamp,
													const int& NForcedDisp,			int* NodeForcedDisp_Array,
													double* PosForcedDisp_Array,	double* PsiA_G_Array);

		//
		void __test_MOD_wrap_input_dynsetup (int& NumSteps,                  double& t0,
                                             double& dt,                     const bool& FollowerForce,
                                             const bool& FollowerForceRig,   const bool& PrintInfo,
                                             const bool& OutInBframe,        const bool& OutInaframe,
                                             const int& ElemProj,            const int& MaxIterations,
                                             const int& NumLoadSteps,        const int& NumGauss,
                                             const int& Solution,            const double& DeltaCurved,
                                             const double& MinDelta,         const double& NewmarkDamp);

		//
		void __test_MOD_wrap_cbeam3_solv_disp2state (const int& NumNodes_tot,        const int& NumDof,
                                                     const int& NumElems,            const int* Nod_Master,
                                                     const int* Nod_Vdof,            const int* Nod_Fdof,
                                                     const double* PosDefor_Array,   const double* PsiDefor_Array,
                                                     const double* PosDeforDot_Array,const double* PsiDeforDot_Array,
                                                     double* X,                      double*DXdt);

		//
		void __test_MOD_wrap_cbeam3_asbly_dynamic (const int& NumElems,              const int& NumNodes_tot,
                                                   const int* NumNodes,              const int* MemNo,
                                                   const int* Conn,                  const int* Master_Array,
                                                   const double* Length,             const double* PreCurv,
                                                   const double* Psi,                const double*Vector,
                                                   const double* Mass_Array,         const double* Stiff_Array,
                                                   const double* InvStiff_Array,     const double* RBMass_Array,
                                                   const int* Nod_Master,            const int* Nod_Vdof,
                                                   const int* Nod_Fdof,              const double* Coords_Array,
                                                   const double* Psi0_Array,         const double* PosDefor_Array,
                                                   const double* PsiDefor_Array,     const double* PosDeforDot_Array,
                                                   const double* PsiDeforDot_Array,  const double* PosDeforDDot_Array,
                                                   const double* PsiDeforDDot_Array, const double* Force_Array,
                                                   const double* Vrel,               const double* VrelDot,
                                                   const int& NumDof,                const int& DimMat,
                                                   int& ms,                          double* Mglobal_Array,
                                                   double* Mvel_Array,               int& cs,
                                                   double* Cglobal_Array,            double* Cvel_Array,
                                                   int& ks,                          double* Kglobal_Array,
                                                   int& fs,                          double* Fglobal_Array,
                                                   double* Qglobal,                  const bool& FollowerForce,
                                                   const bool& FollowerForceRig,     const bool& PrintInfo,
                                                   const bool& OutInBframe,          const bool& OutInaframe,
                                                   const int& ElemProj,              const int& MaxIterations,
                                                   const int& NumLoadSteps,          const int& NumGauss,
                                                   const int& Solution,              const double& DeltaCurved,
                                                   const double& MinDelta,           const double& NewmarkDamp,
                                                   const double* Cao_Array);
 
		//
		void __test_MOD_wrap_cbeam3_solv_state2disp  (const int& NumElems,            const int& NumNodes_tot,
                                                      const int* NumNodes,            const int* MemNo,
                                                      const int* Conn,                const int* Master_Array,
                                                      const double* Length,           const double* PreCurv,
                                                      const double* Psi,              const double* Vector,
                                                      const double* Mass_Array,       const double* Stiff_Array,
                                                      const double* InvStiff_Array,   const double* RBMass_Array,
                                                      const int* Nod_Master,          const int* Nod_Vdof,
                                                      const int* Nod_Fdof,            const double* Coords_Array,
                                                      const double* Psi0_Array,       const int& NumDof,
                                                      const double* X,                const double* DXdt,
                                                      double* PosDefor_Array,         double* PsiDefor_Array,
                                                      double* PosDeforDot_Array,      double* PsiDeforDot_Array);
	
		//
		void __test_MOD_wrap_fem_m2v (const int& N1, const int& N2, const double* Mat_Array, const int& N3, double* Array_out, const int* FilterIN);

		//
		void __test_MOD_wrap_sparse_matvmul (const int& ms, const int& DimMat, const int& NumDof, const double* Mat_Array, const double* Vector, double* Array_out);

		//
		void __test_MOD_wrap_matmul(const int& n0, const int& n1, const int& n2, const double*Array0, const double*Array1, const double*Array2);

		//
		void __test_MOD_wrap_lu_sparse(const int& ms, const int& DimMat, const int& NumDof, const double*Mat_Array, const double* Vec, double*X);

		//
		void __test_MOD_wrap_lu_invers (const int& n, const double* Matrix_vec, double* InvMatrix_vec);

		//
		void __test_MOD_wrap_full_rank2sparse(const int& DimSprMat, const int& n1, const int& n2, const double* FulMat, double* SprMat);

		//
	 	void __test_MOD_wrap_xbeam_rot(const double* q, double* R_Array);

		//
	 	void __test_MOD_wrap_xbeam_quadskew(const int& n, const double* Array0, double* Array1);

		//
		void  __test_MOD_wrap_sparse_zero(int& ms, const int& DimSprMat, const int& NumDof, double* Mat_Array);

		//
		void __test_MOD_wrap_sparse_addsparse(const int& i1,             const int& j1, 
                                              const int& NumDof,         const int& DimMat, 
											  const int& DSubMat, 		 const double* Submat_Array,
											  const int& DMat,			 double* Mat_Array,         
											  double& Factor); 

		//
    	void __test_MOD_wrap_input_dynforce(const int& NumNodes,           const int& nrt,
                                            const double* Time, 	       const double* ForceStatic,
                                            double* ForceDynAmp,    	   double* ForceTime);
    	//
    	void __test_MOD_wrap_input_forcedvel(const int& NumNodes,        const int& nrt,
                                             const double* Time,         double* ForcedVel,
                                             double* ForcedVelDot);
 	}  

//#ifdef NEKTAR_USING_SHARPY

	///
	static inline void Wrap_input_setup (const int& NumElems,		char*  OutFile,
               							 const bool& FollowerForce, const bool& FollowerForceRig,
               							 const bool& PrintInfo, 	const bool& OutInBframe, 
										 const bool& OutInaframe,	const int& ElemProj, 
										 const int& MaxIterations, 	const int& NumLoadSteps,
               							 const int& NumGauss, 		const int& Solution, 
										 const double& DeltaCurved,	const double &MinDelta, 
										 const double& NewmarkDamp) 
	{
		__test_MOD_wrap_input_setup (NumElems,OutFile,
               FollowerForce, FollowerForceRig,        
               PrintInfo, OutInBframe, OutInaframe,    
               ElemProj, MaxIterations, NumLoadSteps,  
               NumGauss, Solution, DeltaCurved,        
               MinDelta, NewmarkDamp);
	}
	
	///
    static inline void Wrap_input_elem (const int& NumElems,        int& NumNodes_tot,
                                        int* NumNodes,              int* MemNo,
                                        int* Conn,                  int* Master_Array,
                                        double* Length,             double* PreCurv,
                                        double* Psi,                double* Vector,
                                        double* Mass_Array,         double* Stiff_Array,
                                        double* InvStiff_Array,     double* RBMass_Array)
	{
		__test_MOD_wrap_input_elem (NumElems,NumNodes_tot,NumNodes,MemNo,Conn,
               Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,
               Stiff_Array,InvStiff_Array,RBMass_Array);
	}

	///
	static inline void Wrap_input_node (const int& NumElems, 			const int& NumNodes_tot,
										const int* NumNodes,			const int* MemNo,
										const int* Conn,				const int* Master_Array,
										const double* Length,			const double* PreCurv,
										const double* Psi,				const double* Vector,
										const double* Mass_Array,		const double* Stiff_Array,
										const double* InvStiff_Array,	const double* RBMass_Array,
										int* BoundConds,				double* PosIni_Array,
										double* ForceStatic_Array,		double*PhiNodes)
	{
		__test_MOD_wrap_input_node (NumElems,NumNodes_tot,NumNodes,MemNo,Conn,
               Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,
               Stiff_Array,InvStiff_Array,RBMass_Array,BoundConds,
               PosIni_Array,ForceStatic_Array,PhiNodes);
	}

	///
	static inline void Wrap_xbeam_undef_geom (const int&NumElems,			int& NumNodes_tot,
											  int* NumNodes,				int* MemNo,
               								  int* Conn,					int* Master_Array,
											  double* Length,				double* PreCurv,
											  double* Psi,					double* Vector,
											  double* Mass_Array,			double* Stiff_Array,
											  double* InvStiff_Array,		double* RBMass_Array,
											  const double* PosIni_Array, 	const double* PhiNodes,
											  double* PsiIni_Array,			const bool& FollowerForce, 
											  const bool& FollowerForceRig,	const bool& PrintInfo, 
											  const bool& OutInBframe, 		const bool& OutInaframe,
               								  const int& ElemProj, 			const int& MaxIterations, 
											  const int& NumLoadSteps,		const int& NumGauss, 
											  const int& Solution, 			const double& DeltaCurved,
               								  const double& MinDelta, 		const double& NewmarkDamp)
	{
		__test_MOD_wrap_xbeam_undef_geom (NumElems,NumNodes_tot,NumNodes,MemNo,
               Conn,Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,
               Stiff_Array,InvStiff_Array,RBMass_Array,PosIni_Array,PhiNodes,PsiIni_Array,
               FollowerForce, FollowerForceRig,PrintInfo, OutInBframe, OutInaframe,    
               ElemProj, MaxIterations, NumLoadSteps,NumGauss, Solution, DeltaCurved,       
               MinDelta, NewmarkDamp); 
	}

	///
	static inline void Wrap_xbeam_undef_dofs (const int& NumElems, 				const int& NumNodes_tot,
											  const int* NumNodes, 				const int* MemNo,
               								  const int* Conn,					const int* Master_Array,
											  const double* Length,				const double* PreCurv,
											  const double* Psi,				const double* Vector,
											  const double* Mass_Array, 		const double* Stiff_Array,
											  const double* InvStiff_Array,		const double* RBMass_Array,
											  const int* BoundConds,			int* Nod_Master,
											  int* Nod_Vdof,					int* Nod_Fdof,
											  int& NumDof, 						int* Nod_Sflag)
	{
		__test_MOD_wrap_xbeam_undef_dofs (NumElems,NumNodes_tot,NumNodes,MemNo,
               Conn,Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,
               Stiff_Array,InvStiff_Array,RBMass_Array,BoundConds,
               Nod_Master,Nod_Vdof,Nod_Fdof,NumDof,Nod_Sflag);   
	}   

	///
	static inline void Wrap_cbeam3_solv_nlnstatic(const int& NumDof,            const int& NumElems,
                                                  const int* NumNodes,          const int* MemNo,
                                                  const int* Conn,              const int* Master_Array,   
                                                  const double* Length,         const double* PreCurv,    
                                                  const double* Psi,            const double* Vector,
                                                  const double* Mass_Array,     const double* Stiff_Array,   
                                                  const double* InvStiff_Array, const double* RBMass_Array,
                                                  const int& NumNodes_tot,      const int* Master,
                                                  const int* Vdof,              const int* Fdof,       //for pack_xbnode
                                                  const double* AppForces_Array,const double*  Coords_Array,
                                                  const double* Psi0_Array,     double* PosDefor_Array,
                                                  double* PsiDefor_Array,       const bool& FollowerForce,
                                                  const bool& FollowerForceRig, const bool& PrintInfo,
                                                  const bool& OutInBframe,      const bool& OutInaframe,
                                                  const int& ElemProj,          const int& MaxIterations,
                                                  const int& NumLoadSteps,      const int& NumGauss,
                                                  const int& Solution,          const double& DeltaCurved,        
                                                  const double& MinDelta,       const double& NewmarkDamp,
                                                  const int& NForcedDisp,       int* NodeForcedDisp_Array,
                                                  double* PosForcedDisp_Array,  double* PsiA_G_Array)	
	{
		__test_MOD_wrap_cbeam3_solv_nlnstatic (NumDof,NumElems, NumNodes, MemNo, Conn,       
                   Master_Array,Length, PreCurv,Psi, Vector, Mass_Array,               
                   Stiff_Array,InvStiff_Array, RBMass_Array,NumNodes_tot, Master, Vdof, Fdof,  
                   AppForces_Array,Coords_Array,Psi0_Array,PosDefor_Array,PsiDefor_Array,
                   FollowerForce, FollowerForceRig,PrintInfo, OutInBframe, OutInaframe,  
                   ElemProj, MaxIterations, NumLoadSteps,NumGauss, Solution, DeltaCurved, 
                   MinDelta, NewmarkDamp, NForcedDisp, NodeForcedDisp_Array, 
				   PosForcedDisp_Array, PsiA_G_Array);
	}
   
	///
	static inline void Wrap_input_dynsetup (int& NumSteps,					double& t0,
											double& dt,						const bool& FollowerForce, 
											const bool& FollowerForceRig,	const bool& PrintInfo, 
											const bool& OutInBframe, 		const bool& OutInaframe,
											const int& ElemProj, 			const int& MaxIterations, 
											const int& NumLoadSteps,		const int& NumGauss, 
											const int& Solution, 			const double& DeltaCurved, 
											const double& MinDelta, 		const double& NewmarkDamp)
	{
		__test_MOD_wrap_input_dynsetup (NumSteps,t0,dt,FollowerForce, FollowerForceRig,       
               		PrintInfo, OutInBframe, OutInaframe,ElemProj, MaxIterations, NumLoadSteps, 
               		NumGauss, Solution, DeltaCurved, MinDelta, NewmarkDamp);       
	}

	///
	static inline void Wrap_cbeam3_solv_disp2state (const int& NumNodes_tot,		const int& NumDof,
													const int& NumElems,			const int* Nod_Master,
													const int* Nod_Vdof,			const int* Nod_Fdof,
               										const double* PosDefor_Array,	const double* PsiDefor_Array,
													const double* PosDeforDot_Array,const double* PsiDeforDot_Array,
               										double* X, 						double*DXdt)
	{
		__test_MOD_wrap_cbeam3_solv_disp2state (
               NumNodes_tot,NumDof,NumElems,Nod_Master,Nod_Vdof,Nod_Fdof,
               PosDefor_Array,PsiDefor_Array,PosDeforDot_Array,PsiDeforDot_Array,
               X,DXdt);
	}

	///
	static inline void Wrap_cbeam3_asbly_dynamic (const int& NumElems,				const int& NumNodes_tot,
												  const int* NumNodes,				const int* MemNo,
												  const int* Conn,					const int* Master_Array,
												  const double* Length,				const double* PreCurv,
												  const double* Psi,				const double*Vector,
												  const double* Mass_Array,			const double* Stiff_Array,
												  const double* InvStiff_Array,		const double* RBMass_Array,
               									  const int* Nod_Master,			const int* Nod_Vdof,
												  const int* Nod_Fdof,				const double* Coords_Array,
												  const double* Psi0_Array,			const double* PosDefor_Array,
												  const double* PsiDefor_Array,		const double* PosDeforDot_Array,
												  const double* PsiDeforDot_Array, 	const double* PosDeforDDot_Array,
												  const double* PsiDeforDDot_Array, const double* Force_Array,
												  const double* Vrel,				const double* VrelDot,
												  const int& NumDof,				const int& DimMat,
               									  int& ms,							double* Mglobal_Array,
												  double* Mvel_Array,				int& cs,
												  double* Cglobal_Array,			double* Cvel_Array,
               									  int& ks,							double* Kglobal_Array,
               									  int& fs,							double* Fglobal_Array,
												  double* Qglobal,					const bool& FollowerForce,
										    	  const bool& FollowerForceRig,		const bool& PrintInfo, 
												  const bool& OutInBframe, 			const bool& OutInaframe,
               									  const int& ElemProj, 				const int& MaxIterations, 
												  const int& NumLoadSteps,			const int& NumGauss, 
												  const int& Solution, 				const double& DeltaCurved,
               									  const double& MinDelta, 			const double& NewmarkDamp,
               									  const double* Cao_Array)
	{
		__test_MOD_wrap_cbeam3_asbly_dynamic (NumElems,NumNodes_tot,NumNodes,
               MemNo,Conn,Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,
               Stiff_Array,InvStiff_Array,RBMass_Array,Nod_Master,Nod_Vdof,Nod_Fdof,
               Coords_Array,Psi0_Array,PosDefor_Array,PsiDefor_Array,
               PosDeforDot_Array,PsiDeforDot_Array,PosDeforDDot_Array,PsiDeforDDot_Array,
               Force_Array,Vrel,VrelDot,NumDof,DimMat,ms,Mglobal_Array,Mvel_Array,
               cs,Cglobal_Array,Cvel_Array,ks,Kglobal_Array,
               fs,Fglobal_Array,Qglobal,FollowerForce, FollowerForceRig,        
               PrintInfo, OutInBframe, OutInaframe,ElemProj, MaxIterations, NumLoadSteps,  
               NumGauss, Solution, DeltaCurved,MinDelta, NewmarkDamp,Cao_Array);
	}

	///
	static inline void Wrap_cbeam3_solv_state2disp (const int& NumElems, 			const int& NumNodes_tot,
													const int* NumNodes, 			const int* MemNo,
													const int* Conn,				const int* Master_Array,
													const double* Length,			const double* PreCurv,
													const double* Psi,				const double* Vector,
													const double* Mass_Array, 		const double* Stiff_Array, 
													const double* InvStiff_Array,	const double* RBMass_Array,
               										const int* Nod_Master,			const int* Nod_Vdof,
													const int* Nod_Fdof,			const double* Coords_Array,
													const double* Psi0_Array,		const int& NumDof,
													const double* X,				const double* DXdt,
               										double* PosDefor_Array,			double* PsiDefor_Array,
													double* PosDeforDot_Array,		double* PsiDeforDot_Array)
	{
		__test_MOD_wrap_cbeam3_solv_state2disp (
               NumElems,NumNodes_tot,NumNodes,
               MemNo,Conn,Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,
               Stiff_Array,InvStiff_Array,RBMass_Array,Nod_Master,Nod_Vdof,Nod_Fdof,
               Coords_Array,Psi0_Array,NumDof,X,DXdt,
               PosDefor_Array,PsiDefor_Array,PosDeforDot_Array,PsiDeforDot_Array);
	}

	///
	static inline void Wrap_fem_m2v(const int& N1, const int& N2, const double* Mat_Array, const int& N3, double* Array_out, const int* FilterIN)
	{
		__test_MOD_wrap_fem_m2v(N1, N2, Mat_Array, N3, Array_out, FilterIN);
	}

	///
	static inline void Wrap_full_rank2sparse(const int& DimSprMat, const int& n1, const int& n2, const double* FulMat, double* SprMat)
	{
		__test_MOD_wrap_full_rank2sparse(DimSprMat,n1,n2,FulMat,SprMat);
	}

	///
	static inline void Wrap_lu_invers (const int& n, const double* Matrix_vec, double* InvMatrix_vec)
	{
		__test_MOD_wrap_lu_invers (n, Matrix_vec, InvMatrix_vec);
	}

	///
	static inline void Wrap_lu_sparse (	const int& ms,  	const int& DimMat, 
										const int& NumDof, 	const double*Mat_Array, 
										const double* Vec, 	double*X)
	{
		__test_MOD_wrap_lu_sparse (ms, DimMat, NumDof, Mat_Array, Vec, X);
	}

	///
	static inline void Wrap_xbeam_QuadSkew(const int& n, const double* Array0, double* Array1)
	{
		__test_MOD_wrap_xbeam_quadskew(n, Array0, Array1);
	}

	///
	static inline void Wrap_xbeam_Rot(const double* q, double* R_Array)
	{
		__test_MOD_wrap_xbeam_rot(q, R_Array);
	}

	///
	static inline void Wrap_matmul (const int& n0, 			const int& n1, 
									const int& n2, 			const double*Array0, 
									const double*Array1, 	const double*Array2)
	{
		__test_MOD_wrap_matmul (n0, n1, n2, Array0, Array1, Array2);
	}

	///
	static inline void Wrap_sparse_matvmul (const int& ms, 				const int&DimMat,          
											const int& NumDof,			const double* Mat_Array,   	
											const double* Vector, 		double* Array_out)
	{
		__test_MOD_wrap_sparse_matvmul (ms, DimMat,	NumDof, Mat_Array, Vector, Array_out);
	}

	///
	static inline void Wrap_sparse_zero (int& ms, const int& DimSprMat, const int& NumDof, double* Mat_Array)
	{
		__test_MOD_wrap_sparse_zero(ms, DimSprMat, NumDof, Mat_Array);
	}

	///
	static inline void Wrap_sparse_addsparse(const int& i1, 			const int& j1, 
											 const int& NumDof,			const int& DimMat,
											 const int& DSubMat, 	    const double* Submat_Array,
											 const int& DMat, 		    double* Mat_Array,			
										     double& Factor)
	{
		__test_MOD_wrap_sparse_addsparse(i1,j1,NumDof,DimMat, DSubMat,Submat_Array,DMat,Mat_Array,Factor);
	}

	///
	static inline void Wrap_input_dynforce(const int& NumNodes, 		const int& nrt, 
											 const double* Time,		const double* ForceStatic,
											 double* ForceDynAmp,		double* ForceTime)
	{
		 __test_MOD_wrap_input_dynforce(NumNodes,nrt,Time,ForceStatic,ForceDynAmp,ForceTime);
	}

	///
	static inline void Wrap_input_forcedvel(const int& NumNodes,		const int& nrt,
											const double* Time,			double* ForcedVel,
											double* ForcedVelDot)
	{
		__test_MOD_wrap_input_forcedvel(NumNodes,nrt,Time,ForcedVel,ForcedVelDot);
	}
//#endif //NEKTAR_USING_SHARPY
}
#endif //NEKTAR_LIB_UTILITIES_SHARPY_HPP

