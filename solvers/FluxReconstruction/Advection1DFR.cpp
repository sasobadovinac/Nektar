///////////////////////////////////////////////////////////////////////////////
//
// File Advection1DFR.cpp
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
// Description: Demo to test the Flux Reconstruction 1D
// Unsteady advection (linear or spatially non linear advection) 
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <FluxReconstruction/Advection1DFR.h>

namespace Nektar
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// CONSTRUCTOR
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Advection1DFR::Advection1DFR(LibUtilities::SessionReaderSharedPtr &vSession)
	{
		// Read in mesh from input file and create an object of class MeshGraph1D
		// to encaplusate the mesh
		graph1D = MemoryManager<SpatialDomains::MeshGraph1D>::AllocateSharedPtr(vSession);
		
		// Feed our spatial discretisation object with the information coming from the session file
		// i.e. initialise all the memory
		Domain = MemoryManager<MultiRegions::DisContField1D>::AllocateSharedPtr(vSession,graph1D,vSession->GetVariable(0));
		
		// Get the total number of physical points (quadrature points)
		// For nodal expansion the number of physical points is equal to the number of coefficients 
		nq = Domain->GetTotPoints();
		
		// Create an array "x" to hold the coordinates value (x_i) and fill it up
		// with the the function GetCoords (y and z are required by default but they are zero if your mesh is
		// a line in a 1D space, if it is a line in a 3D space you may need y and z too)
		x = Array<OneD,double>(nq);
		y = Array<OneD,double>(nq);
		z = Array<OneD,double>(nq);
		
		Domain->GetCoords(x,y,z);
		
		// interface points
		ni = graph1D->GetNvertices();
		ne = Domain->GetExpSize();
		
		K = Array<OneD,double>(ni);
		
		xi = Array<OneD,double>(ni);
		yi = Array<OneD,double>(ni);
		zi = Array<OneD,double>(ni);
		
		for (int i = 0; i < ni; i++)
		{
			graph1D->GetVertex(i)->GetCoords(xi[i],yi[i],zi[i]);
		}
		
		//filling the correction function array
		dgL = Array<OneD,double>(nq/ne+2);
		dgR = Array<OneD,double>(nq/ne+2);
		GFunctionsGrad(dgL,dgR);
		
		// Loading from the session file the functions describing the advection term,
		// the initial condition and the exact solution
		AdveFunc = vSession->GetFunction("Advection",0);
		InitCond = vSession->GetFunction("InitialConditions",0);
		ExSol = vSession->GetFunction("ExactSolution",0);
		
		// Loading time-stepping parameters
		vSession->LoadParameter("InitialTime",Time,0.0);
		vSession->LoadParameter("TimeStep",TimeStep,0.01);
		vSession->LoadParameter("NumSteps",NumSteps,0);
		
		RiemSol = vSession->GetSolverInfo("ReimannSolver");
		
		// Filling two vectors with the value of the prescibed functions at the quadrature points
		Adv  = Array<OneD, Array<OneD, double> >(1);
		Ini  = Array<OneD, Array<OneD, double> >(1);
		Exac = Array<OneD, Array<OneD, double> >(1);
		
		Adv[0]  = Array<OneD,double>(nq);
		Ini[0]  = Array<OneD,double>(nq);
		Exac[0] = Array<OneD,double>(nq);
		
		for(int i = 0; i < nq; ++i)
		{
			Adv[0][i] = AdveFunc->Evaluate(x[i],y[i],z[i],Time);
			Ini[0][i] = InitCond->Evaluate(x[i],y[i],z[i],Time);
		}
		
		// Setting the physical value of the initial condition inside the physical space
		Domain->UpdatePhys() = Ini[0];
		// Setting also the coefficients consequently
		Domain->FwdTrans(Ini[0],Domain->UpdateCoeffs());
	}
	
	//disctructor
	Advection1DFR::~Advection1DFR()
	{
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// FUNCTIONs IMPLEMENTATION
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	void Advection1DFR::EvaluateAdvectionTerm(const Array<OneD, Array<OneD, double> > & inarray,
											  Array<OneD, Array<OneD, double> > & outarray,
											  const double time)
	{
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// inarray is intended to be the solution U at time n, so this RHS calculation we will be used to work out U at time n+1
		// which will then used to feed this fuction at the next step
		// outarray is our advection term (A*dU/dx) corrected with the flux recostruction technique
		// Note: The value of the Advection coefficients are stored in the vector Adv[i].
		// This is because we can have non-constant advection (some spatial non lineraity).
		// The advection coefficients could also be function of time, then the evaluation of Adv[i]
		// must be redone here as we did in the constructor with the function "Evaluate".
		// For a full non-linear advection where A=U, we need just to substitute in the Vmath function "Adv" with "inarray"
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// switching to Peter nomenclature where inarray = Ud
		
		Array<OneD,double> gradUd(nq,0.0);
		Array<OneD,double> fd(nq,0.0);
		Array<OneD,double> udi(2*ni-2,0.0);
		Array<OneD,double> fdi(2*ni-2,0.0);
		Array<OneD,double> fi(ni,0.0);
		
		Array<OneD,double> tmp1,tmp2;
		
		// Taking the gradient of gradUd = dUd/dx elemetally
		//for(int i=0; i<ne; i++)
		//{
		//	Domain->GetExp(i)->PhysDeriv(tmp1 = inarray[0]+i*nq/ne, tmp2 = gradUd+i*nq/ne);
		//}
		
		// Taking the gradient of gradUd = dUd/dx
		Domain->PhysDeriv(inarray[0],gradUd);
		
		for(int i = 0; i < nq; ++i)
		{
			Adv[0][i] = AdveFunc->Evaluate(x[i],y[i],z[i],time);
		}
		
		//calculating the discontinous flux fd = Adv*gradUd
		Vmath::Vmul(nq,Adv[0],1,gradUd,1,fd,1);
		
		//calculate the value of the fd at the interface points
		InterpToInterface(fd,fdi);
		
		//calculate the value of the ud at the interface points
		InterpToInterface(inarray[0],udi);
		
		//calculate the interface fluxes fi
		ReimannSolver(udi,fi);
		
		//calculate correction fluxes
		FluxesReconstruction(fd,fdi,fi,outarray[0]);
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::Projection(const Array<OneD, Array<OneD, double> > & inarray,
								   Array<OneD, Array<OneD, double> > & outarray,
								   const double time) const
	{
		// For DG it is just a copy
		Vmath::Vcopy(nq,inarray[0],1,outarray[0],1);
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::ReimannSolver(const Array<OneD, double> & inarray, 
									  Array<OneD, double> & F)
	{
		Array<OneD,double> a,sumInt,difInt;
		a = Array<OneD,double>(ni);
		sumInt = Array<OneD,double>(ni);
		difInt = Array<OneD,double>(ni);
		
		sumInt[0] = inarray[0]+inarray[2*ni-3];
		sumInt[ni-1] = inarray[0]+inarray[2*ni-3];
		difInt[0] = inarray[0]-inarray[2*ni-3];
		difInt[ni-1] = inarray[0]-inarray[2*ni-3];
		
		for(int i = 1; i < ni-1; i++)
		{
			sumInt[i] = inarray[2*i] + inarray[2*i-1];
			difInt[i] = inarray[2*i] - inarray[2*i-1];
		}
	
		for(int i = 0; i < ni; ++i)
		{
			a[i] = AdveFunc->Evaluate(xi[i],yi[i],zi[i],Time);
			
			if(RiemSol == "Up-Wind"){K[i] = 0.0;}
			else if(RiemSol == "Euler-Centered"){K[i] = 1.0;}
			else if(RiemSol == "Lax-Friedrichs"){K[i] = 1.0 - 1.0/(TimeStep*abs(a[i]));}
			else if(RiemSol == "Lax-Wendroff"){K[i] = 1.0 - (TimeStep*abs(a[i]));}
			else{ASSERTL0(false,"Reimann solver not implemented");}
			
			F[i] = 0.5*a[i]*(sumInt[i]) - 0.5*abs(a[i])*(1-K[i])*(difInt[i]);
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::FluxesReconstruction(const Array<OneD, double> & fd,
											 const Array<OneD, double> & fdi,
											 const Array<OneD, double> & fi,
											 Array<OneD, double> & outarray)
	{
		Array<OneD,double> jL(ni-1,0.0);  // jumps on the left interfaces
		Array<OneD,double> jR(ni-1,0.0);  // jumps on the righ interfaces
		
		Array<OneD,double> tmpDGL(nq/ne,0.0); // temp arrays containing j*dg/dx
		Array<OneD,double> tmpDGR(nq/ne,0.0);
		Array<OneD,double> tmp,tmparray;
		
		// calculating the jumps on the left and right side
		for(int i = 0; i< ni-1; i++)
		{
			jL[i] = fi[i] - fdi[2*i];
			jR[i] = fi[i+1] - fdi[2*i+1];
		}

		for(int i = 0; i < ne; i++)
		{
			Vmath::Smul(nq/ne,jL[i],tmp = dgL + 1,1,tmpDGL,1);
			Vmath::Smul(nq/ne,jR[i],tmp = dgR+1,1,tmpDGR,1);
			Vmath::Vadd(nq/ne,tmp = fd+i*nq/ne,1,tmpDGL,1,tmparray = outarray +i*nq/ne,1);
			Vmath::Vadd(nq/ne,tmparray = outarray +i*nq/ne,1,tmpDGR,1,tmparray = outarray +i*nq/ne,1);
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::GFunctionsGrad(Array<OneD, double> & dGL,
									   Array<OneD, double> & dGR)
	{	
		int np = nq/ne; // number of points per element
		
		int k = np-1;  // expansion order
		
		double sign = pow(-1.0,double(k));
		
		Array<OneD,double> zeros(np,0.0);
		Array<OneD,double> dump(np,0.0);
		//Gauss points
		Polylib::zwgj(&zeros[0],&dump[0],np,0.0,0.0);
		
		// adding the ends point z = -1 & +1
		Array<OneD,double> points(np+2,0.0);
		
		points[0] = -1.0;
		points[np+1] = 1.0;
		
		for(int i=1; i<np+1; i++)
		{
			points[i] = zeros[i-1]; 
		}
		
		Array<OneD,double> Lk(np+2,0.0);
		Array<OneD,double> dLk(np+2,0.0);
		Array<OneD,double> Lk1(np+2,0.0);
		Array<OneD,double> dLk1(np+2,0.0);

		Array<OneD,double> GL(np+2,0.0);
		Array<OneD,double> GR(np+2,0.0);
		
		Polylib::jacobfd(np+2,&(points[0]),&(Lk[0]),&(dLk[0]),k,0.0,0.0);
		Polylib::jacobfd(np+2,&(points[0]),&(Lk1[0]),&(dLk1[0]),k+1,0.0,0.0);
		
		Vmath::Vsub(np+2,Lk,1,Lk1,1,GL,1);
		Vmath::Smul(np+2,0.5*sign,GL,1,GL,1);
		
		Vmath::Vadd(np+2,Lk,1,Lk1,1,GR,1);
		Vmath::Smul(np+2,0.5,GR,1,GR,1);
		
		Vmath::Vsub(np+2,dLk,1,dLk1,1,dGL,1);
		Vmath::Smul(np+2,0.5*sign,dGL,1,dGL,1);
		
		Vmath::Vadd(np+2,dLk,1,dLk1,1,dGR,1);
		Vmath::Smul(np+2,0.5,dGR,1,dGR,1);
		//-------------------------------------------------------------------------------
		// Plot gL and gR
		
		// Plotting G functions
		int numpoints = 1000;
		
		Array<OneD,double> Lkplot(numpoints,0.0);
		Array<OneD,double> dLkplot(numpoints,0.0);
		
		Array<OneD,double> Lk1plot(numpoints,0.0);
		Array<OneD,double> dLk1plot(numpoints,0.0);
		
		Array<OneD,double> GLplot(numpoints,0.0);
		Array<OneD,double> GRplot(numpoints,0.0);
		
		double dx = 2.0/(numpoints-1);
		
		Array<OneD,double> plot(numpoints,0.0);
		plot[0] = -1.0;
		for(int i=1; i< numpoints; i++)
		{
			plot[i] = plot[i-1] + dx; 
		}
		
		
		Polylib::jacobfd(numpoints,&(plot[0]),&(Lkplot[0]),&(dLkplot[0]),k,0.0,0.0);
		Polylib::jacobfd(numpoints,&(plot[0]),&(Lk1plot[0]),&(dLk1plot[0]),k+1,0.0,0.0);
		
		Vmath::Vsub(numpoints,Lkplot,1,Lk1plot,1,GLplot,1);
		Vmath::Smul(numpoints,0.5*sign,GLplot,1,GLplot,1);
		
		Vmath::Vadd(numpoints,Lkplot,1,Lk1plot,1,GRplot,1);
		Vmath::Smul(numpoints,0.5,GRplot,1,GRplot,1);
			
		ofstream outfile1;
		outfile1.open("GLGR.dat");
		for(int i=0; i < numpoints; i++)
		{
			outfile1 << scientific
			<< setw (17) 
			<< setprecision(10) 
			<< plot[i]
			<< "  " 
			<< GLplot[i] 
			<< "  " 
			<< GRplot[i]
			<< endl;
		}
		outfile1 << endl << endl;
		outfile1.close();
		ofstream outfile2;
		outfile2.open("PlotGfunc.p");
		outfile2 << "# Gnuplot script file" << endl;
		outfile2 << "set   autoscale" << endl;                       
		outfile2 << "unset log" << endl;                           
		outfile2 << "unset label" << endl;                          
		outfile2 << "set xtic auto" << endl;                    
		outfile2 << "set ytic auto" << endl;                        
		outfile2 << "set xlabel \"x\"" << endl;
		outfile2 << "set ylabel \"g\"" << endl;
		outfile2 << "plot    \"GLGR.dat\" ";
		outfile2 << "using 1:2 index 0 ";
		outfile2 << " title 'gL' with linespoints lt 3, ";
		outfile2 << "\"GLGR.dat\" ";
		outfile2 << "using 1:3 index 0 ";
		outfile2 << " title 'gR' with linespoints lt 1" << endl;
		outfile2.close();
		
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::InterpToInterface(const Array<OneD, double> & total,
										  Array<OneD, double> & interface)
	{
		LibUtilities::PointsKey InteriorPoints(nq/ne,LibUtilities::eGaussGaussLegendre);
		LibUtilities::BasisKey  ProblemBase(LibUtilities::eGLL_Lagrange,nq/ne-1,InteriorPoints);
		
		LibUtilities::PointsKey InterfacePoints(2,LibUtilities::eGaussLobattoLegendre);
		LibUtilities::BasisKey  LocalBase(LibUtilities::eModified_A,2,InterfacePoints);
		
		for(int i=0; i < ne; i++)
		{
			LibUtilities::Interp1D(ProblemBase,&(total[i*nq/ne]),LocalBase,&(interface[2*i]));
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::AppendOutput(const Array<OneD, double> & approx,
									 const Array<OneD, double> & exact) const
	{
		ofstream outfile;
		outfile.open("Solution.dat");
		for(int i = 0; i < nq; i++)
		{
			outfile << scientific 
			<< setw (17) 
			<< setprecision(10) 
			<< x[i]
			<< "  " 
			<< approx[i] 
			<< "  " 
			<< exact[i] 
			<< endl;
		}
		outfile << endl << endl;
		outfile.close();
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::GenerateGnuplotScript() const
	{
		ofstream outfile;
		outfile.open("FRSolution.p");
			
		outfile << "# Gnuplot script file" << endl;
		outfile << "set   autoscale" << endl;                       
		outfile << "unset log" << endl;                           
		outfile << "unset label" << endl;                          
		outfile << "set xtic auto" << endl;                    
		outfile << "set ytic auto" << endl;                        
		outfile << "set xlabel \"x\"" << endl;
		outfile << "set ylabel \"u\"" << endl;

		outfile << "plot    \"Solution.dat\" ";
		outfile << "using 1:2 index 0 ";
		outfile << " title 'FR solution (t = " << Time << ")' with linespoints lt 3, ";
		outfile << "\"Solution.dat\" ";
		outfile << "using 1:3 index 0 ";
		outfile << " title 'Exact Solution (t = " << Time << ")' with linespoints lt 1" << endl;

		outfile.close();
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::SolutionPrint()
	{
		/////////////////////////////////////////////////////////////
		// Evaluating the error
		double L2, Linf;
		for(int i = 0; i < nq; ++i)
		{
			Exac[0][i] = ExSol->Evaluate(x[i],y[i],z[i],Time);
		}
		
		L2   = Domain->L2(Exac[0]);
		Linf = Domain->Linf(Exac[0]);
		
		/////////////////////////////////////////////////////////////
		// Write solution to file
		AppendOutput(Domain->GetPhys(),Exac[0]);
		GenerateGnuplotScript();
		/////////////////////////////////////////////////////////////
		// Print summary of solution details
		const SpatialDomains::ExpansionMap &expansions = graph1D->GetExpansions(); 
		LibUtilities::BasisKey bkey0 = expansions.begin()->second->m_basisKeyVector[0];
		
		cout << endl;
		cout << "Solving Unsteady 1D Advection with Flux Recostruction"  << endl;
		
		cout << endl;
		cout << "Expansion : " << LibUtilities::BasisTypeMap[bkey0.GetBasisType()] << endl;
		cout << "No. modes : " << bkey0.GetNumModes() << endl;
		cout << "Advection : " << AdveFunc->GetEquation() << endl;
		
		cout << endl;
		
		cout << "              Time Step : " << TimeStep << endl;
		cout << "        Number of Steps : " << NumSteps << endl;
		cout << "      Initial Condition : " << InitCond->GetEquation() << endl;
		
		cout << endl;
		cout << "Exact Solution : " << ExSol->GetEquation() << endl;
		cout << "    Linf error : " << L2 << endl;
		cout << "    L2   error : " << Linf << endl;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	MultiRegions::DisContField1DSharedPtr Advection1DFR::GetDomain() const
	{
		return Domain;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double Advection1DFR::GetTime() const
	{
		return Time;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::UpdateTime()
	{
		Time = Time + TimeStep;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double Advection1DFR::GetTimeStep() const
	{
		return TimeStep;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int Advection1DFR::GetNumSteps() const
	{
		return NumSteps;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int Advection1DFR::GetNumPoints() const
	{
		return nq;
	}
	//////////////////////////////////////////////////////////////////////////////	
}