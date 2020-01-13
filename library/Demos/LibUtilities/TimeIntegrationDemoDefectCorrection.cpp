///////////////////////////////////////////////////////////////////////////////
//
// File TimeIntegrationDemoDefectCorrection.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
// Description: Example of using time-integration schemes.
//
///////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------
// This file illustrates the use of the time-stepping framework.
//
// For more information, please consult the following reference:
//
// Vos, P.E.J., Eskilsson, C., Bolis, A., Chun, S., Kirby, R.M. and Sherwin, S.J.
// "A Generic Framework for Time-Stepping PDEs: general linear methods,
//  object-oriented implementation and application to fluid problems"
// International Journal of Computational Fluid Dynamics, to appear

// It solves the one-dimensional advection-diffusion problem, defined as
//
//  |    du     du     d^2 u
//  |    -- + U -- = D -----,
//  |    dt     dx     d x^2
//  |
//  |    subject to:
//  |    - periodic boundary conditions
//  |    - the initial condition
//  |        u(x,0) = sin(2*pi*k*x)
//  |
//  |    and with  x = [0,1]
//  |              t = [0,1]
//  |              U = 1
//  |              D = 0.05
//  |              k = 1
//  |
//
// using the finite difference method.
// The exact solution of the equation above is given by
//
//  u(x,t) = exp(-D * (2*pi*k)^2 * t) * sin(2*pi*k * (x - U*t))
//
// The output is written out to the files
//
//   - OneDfinDiffAdvDiffSolverOutput.dat (containing the data)
//   - OneDfinDiffAdvDiffSolverOutput.p   (containing a gnuplot script)
//
// and can be visualised by gnuplot using the command
//
//    gnuplot OneDfinDiffAdvDiffSolverOutput.p
//
// To run this code in perfect setting type, e.g., in the terminal:
// > TimeIntegrationOrderDemoDefectCorrection -p 100 -t 100  -m 2
//-------------------------------------------------------------------------------------------------
#include <fstream>
#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/program_options.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeOperators.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSolution.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;

namespace po = boost::program_options;

//Global array containing the defect-correction values to be recognised everywhere. Merely a hack!!
Array<OneD, double> null(1);
Array<OneD, double> dc2(1); 
Array<OneD, double> dc3(1); 
Array<OneD, double> anstpred(1);
Array<OneD, double> anstcorr(1);

// We first implement a class that represents
// the 1D finite difference solver
class OneDfinDiffAdvDiffSolver
{
public:
    // constructor based upon the discretisation details
    OneDfinDiffAdvDiffSolver(int nPoints, int nTimeSteps):
        m_x0(0.0),
        m_xend(1.0),
        m_nPoints(nPoints),
        m_dx((m_xend-m_x0)/((double) m_nPoints-1.0)),
        m_t0(0.0),
        m_tend(1.0),
        m_nTimeSteps(nTimeSteps),
        m_dt((m_tend-m_t0)/(double) m_nTimeSteps),
        m_wavenumber(1.0),
        m_U(1.0),
        m_D(0.05)
    {
    }

    // ---------------------------------------------------------------------------------------------
    // ---- These functions/methods below are the routines which will be used by
    // ---- the TimeIntegration framework (and are required for using it)...
    // ---- The implementation of these functions can be found at the end of the file
    
    void HelmSolveDc(const Array<OneD, const Array<OneD, double> >& inarray,
                           Array<OneD,       Array<OneD, double> >& outarray,
                     const NekDouble time,
                     const NekDouble lambda) const;

    void HelmSolveImex(const Array<OneD, const Array<OneD, double> >& inarray,
                             Array<OneD,       Array<OneD, double> >& outarray,
                       const NekDouble time,
                       const NekDouble lambda) const;

    void EvaluateAdvectionDc(const Array<OneD, const Array<OneD, double> >& inarray,
                                   Array<OneD,       Array<OneD, double> >& outarray,
                             const NekDouble time) const;   

    void EvaluateAdvectionImex(const Array<OneD, const  Array<OneD, double> >& inarray,                                                     
                                     Array<OneD,        Array<OneD, double> >& outarray,
                               const NekDouble time) const; 

    void EvaluateAdvectionSp(const Array<OneD, const double>& predictor,
                             const Array<OneD, const double>& corrector,
                                   Array<OneD, double>& outarray) const;
    
    int GetNpoints() const;        

    void EvaluateExactSolution(Array<OneD, double>& outarray,
                               const NekDouble time) const;    

    double EvaluateL2Error(const Array<OneD, const double>& approx,
                           const Array<OneD, const double>& exact) const;    
 
    void AppendOutput(ofstream& outfile,
                      const Array<OneD, const double>& approx,
                      const Array<OneD, const double>& exact) const;    

    void AsympExpansion(const Array<OneD, const double>& inarray0,
                        const Array<OneD, const double>& inarray1,
                        const Array<OneD, const double>& inarray2,
                              Array<OneD, double>& outarray,
                        const NekDouble dt) const;   
   
    void GenerateGnuplotScript(NekDouble levelref) const;    

    double GetInitialTime() const;

    double GetFinalTime() const;
    
    double GetTimeStep() const;
    
    // --------------------------------------------------------------------------------------------

private:
    // spatial discretisation
    double m_x0;         // the left boundary of the domain
    double m_xend;       // the right boundary of the domain
    int    m_nPoints;    // the number of grid-points used in the finite difference method
    double m_dx;         // the distance between 2 grid points

    // temporal discretisation
    double m_t0;         // the initial time
    double m_tend;       // the end time
    int    m_nTimeSteps; // the number of time-steps
    double m_dt;         // the size of a time-step

    // value of the coefficients
    double m_wavenumber; // wave number
    double m_U;          // advection speed
    double m_D;          // diffusion coefficient

    void solveTriDiagMatrix (int n, double alpha, double beta,
                             const Array<OneD, const double>& inarray,
                             Array<OneD,       double>& outarray) const;

    int factorial(int n) const;
        
};

int main(int argc, char *argv[])
{
    po::options_description desc("Usage:");
    desc.add_options()("help,h", "Produce this help message.")(
        "points,p", po::value<int>(), "Number of grid points to be used.")(
        "timesteps,t", po::value<int>(), "Number of timesteps to be used.")(
        "method,m", po::value<int>(),
        "TimeIntegrationMethod is a number in the range [1,9].\n"
        "It defines the time-integration method to be used:\n"
        "- 1: 1st order multi-step IMEX scheme\n"
        "     (Euler Backwards/Euler Forwards)\n"
        "- 2: 2nd order multi-step IMEX scheme\n"
        "- 3: 3rd order multi-step IMEX scheme\n"
        "- 4: 4th order multi-step IMEX scheme\n"
        "- 5: 2nd order multi-stage DIRK IMEX scheme\n"
        "- 6: 3nd order multi-stage DIRK IMEX scheme\n"
        "- 7: 2nd order IMEX Gear (Extrapolated Gear/SBDF-2)\n"
        "- 8: 2nd order Crank-Nicolson/Adams-Bashforth (CNAB)\n"
        "- 9: 2nd order Modified Crank-Nicolson/Adams-Bashforth\n"
        "     (MCNAB)"
        "-10: 2nd order Defect Correction IMEX scheme\n"
        "-11: 3rd order Defect Correction IMEX scheme\n");

    po::variables_map vm;
    try
    {
        po::store(po::command_line_parser(argc, argv).
                  options(desc).run(), vm);
        po::notify(vm);
    }
    catch (const std::exception& e)
    {
        cerr << e.what() << endl;
        cerr << desc;
        return 1;
    }

    if (!vm.count("points") || !vm.count("timesteps") || !vm.count("method")
            || vm.count("help"))
    {
        cout << "Please specify points, timesteps and method." << endl << endl;
        cout << desc;
        return 1;
    }

    int nPoints = vm["points"].as<int>();
    int nTimesteps = vm["timesteps"].as<int>();
    int nMethod = vm["method"].as<int>();

    // Open a file for writing the solution
    ofstream outfile;
    outfile.open("OneDfinDiffAdvDiffSolverOutput.dat");

    // --------------------------------------------------------------------------------------------
    // THE IMPLEMENTATION BELOW SHOWS HOW THE TIME-STEPPING FRAMEWORK CAN BE
    // USED FOR TIME-INTEGRATION PDEs

    // 1. THE SPATIAL DISCRETISATION
    //    Create an object of the OneDfinDiffAdvDiffSolver class.
    //    This class can be thought of as representing the
    //    spatial (finite difference) discretisation.
    OneDfinDiffAdvDiffSolver* solver = new OneDfinDiffAdvDiffSolver(nPoints,nTimesteps);
    //    After this spatial discretisation, the PDE has actually been
    //    reduced (through the method-of-lines) to an ODE.
    //    In order to use the time-stepping framework, we need to give it the necessary
    //    information about this ODE.
    //    Therefore, we create an object of the class TimeIntegrationSchemeOperators that
    //    contains a 'function pointer' (in fact a 'functor') to the
    //    - explicit term of the ODE (i.e. the advection term)
    //    - implicit solve routine (i.e. the Helmholtz solver)
    //    - projection operator (i.e. the identity operator in this case)
    LibUtilities::TimeIntegrationSchemeOperators ode;    
    ode.DefineOdeRhs        (&OneDfinDiffAdvDiffSolver::EvaluateAdvectionImex,solver);    
    ode.DefineImplicitSolve (&OneDfinDiffAdvDiffSolver::HelmSolveImex,        solver);

    LibUtilities::TimeIntegrationSchemeOperators odedc1;    
    odedc1.DefineOdeRhs        (&OneDfinDiffAdvDiffSolver::EvaluateAdvectionDc,solver);    
    odedc1.DefineImplicitSolve (&OneDfinDiffAdvDiffSolver::HelmSolveDc,        solver);
    
    LibUtilities::TimeIntegrationSchemeOperators odedc2;
    odedc2.DefineOdeRhs        (&OneDfinDiffAdvDiffSolver::EvaluateAdvectionDc,solver);    
    odedc2.DefineImplicitSolve (&OneDfinDiffAdvDiffSolver::HelmSolveDc,        solver);

    LibUtilities::TimeIntegrationSchemeOperators odedc3;
    odedc3.DefineOdeRhs        (&OneDfinDiffAdvDiffSolver::EvaluateAdvectionDc,solver);    
    odedc3.DefineImplicitSolve (&OneDfinDiffAdvDiffSolver::HelmSolveDc,        solver);
    
    // 2. THE TEMPORAL DISCRETISATION
    // 2.1 Read in which method should be used.
    //     For a multi-step scheme, also set up
    //     which method should be used for appropriately
    //     starting up the system

    LibUtilities::TimeIntegrationSchemeSharedPtr tiScheme;

    TimeIntegrationSchemeFactory &factory = LibUtilities::GetTimeIntegrationSchemeFactory();
    
    switch (nMethod)
    {
        case 1:
            tiScheme = factory.CreateInstance("IMEXOrder1");
            break;
        case 2:
            tiScheme = factory.CreateInstance("IMEXOrder2");
            break;
        case 3:
            tiScheme = factory.CreateInstance("IMEXOrder3");
            break;
        case 4:
            tiScheme = factory.CreateInstance("IMEXOrder4");
            break;
        case 5:
            tiScheme = factory.CreateInstance("IMEXdirk_2_3_2");
            break;
        case 6:
            tiScheme = factory.CreateInstance("IMEXdirk_3_4_3");
            break;
        case 7:
            tiScheme = factory.CreateInstance("IMEXGear");
            break;
        case 8:
            tiScheme = factory.CreateInstance("CNAB");
            break;
        case 9:
            tiScheme = factory.CreateInstance("MCNAB");
            break;
        case 10:
        case 11:
            tiScheme = factory.CreateInstance("IMEXOrder1");
            break;    
        default :
        {
            cout << "Invalid method." << endl << endl;
            cout << desc;
            exit(1);
        }
    }

    // 2.2 Create objects of the time-integration framework.
    //     These can later be used for actually doing the time-integration
    
    double l2error_coarser = 1.0;
    double l2error_finer = 1.0;
    double l2orderconv;

    int level = 5; //Set the level of refinement for the reference solution
        
    // 2.31 Initialise some arrays that contain the numerical and
    //      analytical solutions

        Array<OneD, Array<OneD, double> > fidifsol(1);   // Array containing the numerical solution 
        Array<OneD, Array<OneD, double> > fidifsolsp1(1);// Array containing the numerical solution (subproblem 1)
        Array<OneD, Array<OneD, double> > fidifsolsp2(1);// Array containing the numerical solution (subproblem 2)
        Array<OneD, Array<OneD, double> > fidifsolsp3(1);// Array containing the numerical solution (subproblem 3)
        fidifsol[0] = Array<OneD, double>(nPoints);        
        fidifsolsp1[0] = Array<OneD, double>(nPoints);   
        fidifsolsp2[0] = Array<OneD, double>(nPoints);
        fidifsolsp3[0] = Array<OneD, double>(nPoints); 

        Array<OneD, double> fidifrefsol(1);// Array containing the reference numerical solution         
        Array<OneD, double> exactsol(1);   // Array containing the exact solution  
        fidifrefsol = Array<OneD, double>(nPoints);   
        exactsol = Array<OneD, double>(nPoints);
          
        solver->AppendOutput(outfile,fidifsol[0],exactsol); // Write the initial condition to a file        
        
        // Auxiliary variables for DC-3 method (zeroth-order)
        Array<OneD, Array<OneD, double> > fidifsol00(1); // Array containing the numerical solution u00
        Array<OneD, Array<OneD, double> > fidifsol01(1); // Array containing the numerical solution u01
        Array<OneD, Array<OneD, double> > fidifsol02(1); // Array containing the numerical solution u02
        Array<OneD, Array<OneD, double> > fidifsol03(1); // Array containing the numerical solution u03
        Array<OneD, Array<OneD, double> > fidifsol10(1); // Array containing the numerical solution u10
        Array<OneD, Array<OneD, double> > fidifsol11(1); // Array containing the numerical solution u11
        Array<OneD, Array<OneD, double> > fidifsol12(1); // Array containing the numerical solution u12
        Array<OneD, Array<OneD, double> > fidifsol20(1); // Array containing the numerical solution u20
        Array<OneD, Array<OneD, double> > fidifsol21(1); // Array containing the numerical solution u21
        fidifsol00[0] = Array<OneD, double>(nPoints);  
        fidifsol01[0] = Array<OneD, double>(nPoints);  
        fidifsol02[0] = Array<OneD, double>(nPoints);  
        fidifsol03[0] = Array<OneD, double>(nPoints); 
        fidifsol10[0] = Array<OneD, double>(nPoints);  
        fidifsol11[0] = Array<OneD, double>(nPoints);  
        fidifsol12[0] = Array<OneD, double>(nPoints); 
        fidifsol20[0] = Array<OneD, double>(nPoints);
        fidifsol21[0] = Array<OneD, double>(nPoints);

        // Auxiliary variables for DC-3 method (first-order)
        Array<OneD, double> fidifsold01(nPoints); // Array containing the numerical solution du01
        Array<OneD, double> fidifsold02(nPoints); // Array containing the numerical solution du02 
        Array<OneD, double> fidifsold03(nPoints); // Array containing the numerical solution du03 
        Array<OneD, double> fidifsold11(nPoints); // Array containing the numerical solution du11
        Array<OneD, double> fidifsold12(nPoints); // Array containing the numerical solution du12        

        // Auxiliary variables for DC-3 method (second-order)        
        Array<OneD, double> fidifsold202(nPoints);// Array containing the numerical solution d2u02
        Array<OneD, double> fidifsold203(nPoints);// Array containing the numerical solution d2u03
        Array<OneD, double> fidifsold212(nPoints);// Array containing the numerical solution d2u12             

        // Auxiliary variables for DC-3 method (third-order)
        Array<OneD, double> fidifsold303(nPoints);// Array containing the numerical solution d3u03
              
        // Auxiliary variables for nonlinear ansatz
        Array<OneD, double> fidifsolnl01(nPoints);// Array containing the nonlinear term nlu01
        Array<OneD, double> fidifsolnl02(nPoints);// Array containing the nonlinear term nlu02
        Array<OneD, double> fidifsolnl03(nPoints);// Array containing the nonlinear term nlu03
        Array<OneD, double> fidifsolnl11(nPoints);// Array containing the nonlinear term nlu11
        Array<OneD, double> fidifsolnl12(nPoints);// Array containing the nonlinear term nlu12
        Array<OneD, double> fidifsolnl21(nPoints);// Array containing the nonlinear term nlu21
        null = Array<OneD, double>(nPoints,0.0);  // Array containing the null vector                 
        dc2 = Array<OneD, double>(nPoints);
        dc3 = Array<OneD, double>(nPoints); 
        anstpred = Array<OneD, double>(nPoints);  
        anstcorr = Array<OneD, double>(nPoints);        
    
    //  2.3 Looping over the time-integration routines for the convergence analysis (by refining time-steps)
    for(int j = level ; j >= 0; j--)
    {  
        double dt = solver->GetTimeStep();  
        double t0 = solver->GetInitialTime();
        double tend = solver->GetFinalTime(); 
        dt = dt / pow(2,j);          
               
        solver->EvaluateExactSolution(exactsol,t0); 
        Vmath::Vcopy(nPoints,exactsol,1,fidifsol[0],1);  
        Vmath::Vcopy(nPoints,exactsol,1,fidifsolsp1[0],1);        
           
        // 2.32 Initialize the time-integration scheme 
        LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr sol;      
        LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr soldc1;
        LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr soldc2;  
        LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr soldc3;  
        sol= tiScheme->InitializeScheme(dt,fidifsol,t0,ode);             
        soldc1 = tiScheme->InitializeScheme(dt,fidifsolsp1,t0,odedc1);
        soldc2 = tiScheme->InitializeScheme(dt,fidifsolsp2,t0,odedc2); 
        soldc3 = tiScheme->InitializeScheme(dt,fidifsolsp3,t0,odedc3); 
        // 2.33 Do the time-integration        
    
        int dcstart;
        nMethod == 11 ? dcstart = -1 : dcstart = 0;
        
        Vmath::Vcopy(nPoints,fidifsolsp1[0],1,fidifsol00[0],1);// Assign the initial condition to u00 
        Vmath::Vcopy(nPoints,fidifsolsp2[0],1,fidifsol10[0],1);// Assign the initial condition to u10  
        Vmath::Vcopy(nPoints,fidifsolsp3[0],1,fidifsol20[0],1);// Assign the initial condition to u20  

        nTimesteps = (tend - t0)/dt; // Calculating the number of time-discretization for each refining dt        
        double t = t0;

        for(int i = dcstart; i < nTimesteps + 1; i++)
        {

            int h = ((i - dcstart) < 1) ? (i - dcstart) : 0; // For multi-step schemes and DC schemes                
                                                           
            //==================================== Main body of DC-2 methods ====================================
            //Main body of DC-3 methods
            if( ( i == 0 && nMethod == 10) || (i == -1 &&  nMethod == 11) )
            {
                // Subproblem 1                    
                t = t + dt;    
                Vmath::Vcopy(nPoints,fidifsol00[0],1,fidifsolnl01,1); 
                Vmath::Vcopy(nPoints,fidifsolnl01,1,anstpred,1);                           
                fidifsol01 = tiScheme->TimeIntegrate(h, dt, soldc1, odedc1);  
                Vmath::Svtsvtp(nPoints,1/dt,fidifsol01[0],1,-1/dt,fidifsol00[0],1,fidifsold01,1);// du01 = (u01 - u00)/dt                                  
            }

            else if( (i > 0 && nMethod == 10) || (i == 0 && nMethod == 11) )
            {
                // Subproblem 1                
                t = t + dt;                
                solver->AsympExpansion(fidifsol01[0], fidifsol10[0], null, fidifsolnl02, dt);// nl03 = u02 + dt*u11  
                Vmath::Vcopy(nPoints,fidifsolnl02,1,anstpred,1);                                                   
                fidifsol02 = tiScheme->TimeIntegrate(h, dt, soldc1, odedc1);                   
                Vmath::Svtsvtp(nPoints,1/dt,fidifsol02[0],1,-1/dt,fidifsol01[0],1,fidifsold02,1);// du02 = (u02 - u01)/dt

                // Subproblem 2                
                t = t - dt;                
                Vmath::Svtsvtp(nPoints,1/dt,fidifsold02,1,-1/dt,fidifsold01,1,fidifsold202,1);// du202 = (du02 - du01)/dt              
                Vmath::Smul(nPoints,-dt/2,fidifsold202,1,dc2,1);                 
                solver->AsympExpansion(fidifsol01[0], fidifsol10[0], null, fidifsolnl11, dt);// nl03 = u02 + dt*u11 
                Vmath::Vcopy(nPoints,fidifsolnl11,1,anstpred,1); 
                Vmath::Vcopy(nPoints,fidifsolnl01,1,anstcorr,1);
                fidifsol11 = tiScheme->TimeIntegrate(h, dt,soldc2,odedc2);             
                Vmath::Svtsvtp(nPoints,1/dt,fidifsol11[0],1,-1/dt,fidifsol10[0],1,fidifsold11,1);// du11 = (u11 - u10)/dt

                if(nMethod == 10)
                {                     
                    solver->AsympExpansion(fidifsol01[0], fidifsol11[0], null, fidifsol[0], dt);// nl03 = u02 + dt*u11            

                    // Updating auxiliary variables
                    Vmath::Vcopy(nPoints,fidifsol02[0],1,fidifsol01[0],1);
                    Vmath::Vcopy(nPoints,fidifsol11[0],1,fidifsol10[0],1);                                
                    Vmath::Vcopy(nPoints,fidifsold02,1,fidifsold01,1);                    
                    Vmath::Vcopy(nPoints,fidifsolnl02,1,fidifsolnl01,1);  

                    solver->EvaluateExactSolution(exactsol,t);      // Calculate the exact solution                                
                    solver->AppendOutput(outfile,fidifsol[0],exactsol);// Dump the output to a file           
                    //Storing the reference solution at smallest dt/j^level                  
                    j == level ? Vmath::Vcopy(nPoints,fidifsol[0],1,fidifrefsol,1):Vmath::Vcopy(nPoints,fidifrefsol,1,fidifrefsol,1); 
                }
            }

            else if( (i > 0) && (nMethod == 11) )
            {
                // Subproblem 1                 
                t = t + 2.0 * dt;
                solver->AsympExpansion(fidifsol02[0], fidifsol11[0], null, fidifsolnl03, dt);// nl03 = u02 + dt*u11 
                Vmath::Vcopy(nPoints,fidifsolnl03,1,anstpred,1);                                          
                fidifsol03 = tiScheme->TimeIntegrate(h,dt,soldc1,odedc1);                    
                Vmath::Svtsvtp(nPoints,1/dt,fidifsol03[0],1,-1/dt,fidifsol02[0],1,fidifsold03,1);// du03 = (u03 - u02)/dt      

                // Subproblem 2                
                t = t - dt;                  
                Vmath::Svtsvtp(nPoints,1/dt,fidifsold03,1,-1/dt,fidifsold02,1,fidifsold203,1);// du203 = (du03 - du02)/dt
                Vmath::Smul(nPoints,-dt/2,fidifsold203,1,dc2,1);             
                solver->AsympExpansion(fidifsol02[0], fidifsol11[0],fidifsol20[0],fidifsolnl12,dt); // nl11 = u02 + dt*u11 + (dt^2)*u20
                Vmath::Vcopy(nPoints,fidifsolnl12,1,anstpred,1); 
                Vmath::Vcopy(nPoints,fidifsolnl02,1,anstcorr,1);                
                fidifsol12 = tiScheme->TimeIntegrate(h,dt,soldc2,odedc2);                    
                Vmath::Svtsvtp(nPoints,1/dt,fidifsol12[0],1,-1/dt,fidifsol11[0],1,fidifsold12,1);// du12 = (u12 - u11)/dt

                // Subproblem 3                
                t = t - dt;                   
                Vmath::Svtsvtp(nPoints,1/dt,fidifsold12,1,-1/dt,fidifsold11,1,fidifsold212,1);// du212 = (du12 - du11)/dt                 
                Vmath::Svtsvtp(nPoints,1/dt,fidifsold203,1,-1/dt,fidifsold202,1,fidifsold303,1);// du303 = (d203 - d202)/dt 
                Vmath::Smul(nPoints,-dt/2,fidifsold212,1,dc2,1); 
                Vmath::Smul(nPoints, dt/6,fidifsold303,1,dc3,1);
                
                solver->AsympExpansion(fidifsol01[0],fidifsol11[0],fidifsol20[0],fidifsolnl21,dt);// nl11 = u01 + dt*u11 + (dt^2)*u20
                Vmath::Vcopy(nPoints,fidifsolnl21,1,anstpred,1); 
                Vmath::Vcopy(nPoints,fidifsolnl11,1,anstcorr,1);                 
                fidifsol21 = tiScheme->TimeIntegrate(h, dt, soldc3, odedc3);               
                solver->AsympExpansion(fidifsol01[0],fidifsol11[0],fidifsol21[0],fidifsol[0],dt);// u = u01 + dt*u11 + (dt^2)*u21            

                // Updating auxiliary variables            
                Vmath::Vcopy(nPoints,fidifsol02[0],1,fidifsol01[0],1); 
                Vmath::Vcopy(nPoints,fidifsol03[0],1,fidifsol02[0],1);            
                Vmath::Vcopy(nPoints,fidifsol12[0],1,fidifsol11[0],1); 
                Vmath::Vcopy(nPoints,fidifsol21[0],1,fidifsol20[0],1);             
                Vmath::Vcopy(nPoints,fidifsold03,1,fidifsold02,1);
                Vmath::Vcopy(nPoints,fidifsold12,1,fidifsold11,1);   
                Vmath::Vcopy(nPoints,fidifsold203,1,fidifsold202,1);  
                Vmath::Vcopy(nPoints,fidifsolnl03,1,fidifsolnl02,1);  
                Vmath::Vcopy(nPoints,fidifsolnl12,1,fidifsolnl11,1);                
                solver->EvaluateExactSolution(exactsol,t);      // Calculate the exact solution                                
                solver->AppendOutput(outfile,fidifsol[0],exactsol);// Dump the output to a file           
                //Storing the reference solution at smallest dt/j^level                  
                j == level ? Vmath::Vcopy(nPoints,fidifsol[0],1,fidifrefsol,1):Vmath::Vcopy(nPoints,fidifrefsol,1,fidifrefsol,1);            
            }  

            else if( (i >= 0) && (i < nTimesteps) && (nMethod != 10) && (nMethod != 11))
            {
                t = t0 + (i+1) * dt;            
                fidifsol = tiScheme->TimeIntegrate(i, dt,sol,ode); // Time-integration for 1 time-step
                solver->EvaluateExactSolution(exactsol,t);       // Calculate the exact solution
                solver->AppendOutput(outfile,fidifsol[0],exactsol); // Dump the output to a file
                j == level ? Vmath::Vcopy(nPoints,fidifsol[0],1,fidifrefsol,1):Vmath::Vcopy(nPoints,fidifrefsol,1,fidifrefsol,1);                  
            }
            
            // Updating the local time step
            t += (i > 0) ? dt : 0;
                
        }                
            l2error_coarser = l2error_finer;
            l2error_finer = solver->EvaluateL2Error(fidifrefsol,fidifsol[0]);        
            l2orderconv = log(l2error_finer/l2error_coarser)/log(2.0);         
        
            // Screen dump **************************************************************************************            
            cout << "j :"    << setw(5) << j 
                 << "| dt :" << setw(10) << setprecision(8) << dt
                 << "| dx :"  << setw(10) << 1.0/nPoints 
                 << "| L 2 error (exact) :" <<  setw(10) 
                 << solver->EvaluateL2Error(exactsol,fidifsol[0])
                 << "| L 2 error (ref) :" <<  setw(10)
                 << solver->EvaluateL2Error(fidifrefsol,fidifsol[0]);
            cout << "| order (ref) : ";
            //Skipped the first two level of convergence analysis 
            (j < level-1)? cout << setw(5) << l2orderconv << endl: cout << " -- " << endl;               
            //***************************************************************************************************       
            solver->GenerateGnuplotScript(level);
            outfile.close();
            delete solver;               
    }   
    return 0;
}


void OneDfinDiffAdvDiffSolver::HelmSolveDc(const Array<OneD, const Array<OneD, double> >& inarray,
                                                 Array<OneD,       Array<OneD, double> >& outarray,
                                           const NekDouble time,
                                           const NekDouble lambda) const
{
    boost::ignore_unused(time);
    // This function implements a 1D finite difference helmholtz solver.
    // The 1D Helmholtz equation leads to a cyclic triadiagonal matrix to be
    // solved.
    // In this function, this is solved as follows:
    // - Based upon Ahlberg-Nilson-Walsh algorithm, we first solve
    //   for the periodic grid-points (i.e. the boundary nodes)
    // - Next, we solve for the interior grid-points. The associated tridiagonal
    //   system is solved based upon the Thomas algorithm
    double a = - m_D * lambda / (m_dx * m_dx); //off diagonal term
    double b = 1.0 + 2.0 * lambda * m_D / (m_dx * m_dx); // diagonal term of triadiagonal matrix
    double lambp = 1.0; //to control the power of lambda
    int nIntPoints = m_nPoints-2; //number of internal points minus the boundary points
    
    Array<OneD, double> invD_f(nIntPoints);    
    Array<OneD, double> rhs(m_nPoints);    

    //Evaluating the correction ansatz for the nonlinear term    

    if(dc2[1] == 0.0 && dc3[1] == 0.0)
    {
        lambp = 1.0; 
        EvaluateAdvectionSp(anstpred, null, rhs);
    }
    else if(dc2[1] != 0.0 && dc3[1] == 0.0)
    {
        lambp = 0.0;
        EvaluateAdvectionSp(anstpred, anstcorr, rhs);
    }
    else 
    {
        lambp = -1.0;
        EvaluateAdvectionSp(anstpred, anstcorr, rhs);
    }    

    for (int i = 0; i < m_nPoints; i++)
    {
        rhs[i] = rhs[i] * pow(lambda, lambp) + inarray[0][i] + dc2[i] + dc3[i];
        //Resetting the defect-correction value to zero
        dc2[i] = 0.0;
        dc3[i] = 0.0;                
    }    
    
    solveTriDiagMatrix(nIntPoints,a,b,rhs+1,invD_f);
    Array<OneD, double> C(nIntPoints,0.0);
    Array<OneD, double> invD_C(nIntPoints,0.0);    
    C[0]     = a;
    C[nIntPoints-1] = a;
    solveTriDiagMatrix(nIntPoints,a,b,C,invD_C);
    
    outarray[0][0] = ( rhs[0] - a * (invD_f[0] + invD_f[nIntPoints-1]) ) / (b - a*(invD_C[0] + invD_C[nIntPoints-1]) );
    outarray[0][m_nPoints-1] = outarray[0][0];

    //Solve for the internal points
    for(int i = 1; i < m_nPoints-1; i++)
    {
        outarray[0][i] = invD_f[i-1] - invD_C[i-1] * outarray[0][0];
    }       
    
}

void OneDfinDiffAdvDiffSolver::HelmSolveImex(const Array<OneD, const Array<OneD, double> >& inarray,
                                                   Array<OneD,       Array<OneD, double> >& outarray,
                                             const NekDouble time,
                                             const NekDouble lambda) const
{
    boost::ignore_unused(time);
    // This function implements a 1D finite difference helmholtz solver.
    // The 1D Helmholtz equation leads to a cyclic triadiagonal matrix to be
    // solved.
    // In this function, this is solved as follows:
    // - Based upon Ahlberg-Nilson-Walsh algorithm, we first solve
    //   for the periodic grid-points (i.e. the boundary nodes)
    // - Next, we solve for the interior grid-points. The associated tridiagonal
    //   system is solved based upon the Thomas algorithm
    double a = - m_D * lambda / (m_dx * m_dx); //off diagonal term
    double b = 1.0 + 2.0 * lambda * m_D / (m_dx * m_dx); // diagonal term of triadiagonal matrix
    
    int nIntPoints = m_nPoints-2; //number of internal points minus the boundary points
    
    Array<OneD, double> invD_f(nIntPoints);    
    solveTriDiagMatrix(nIntPoints,a,b,inarray[0]+1,invD_f);
    Array<OneD, double> C(nIntPoints,0.0);
    Array<OneD, double> invD_C(nIntPoints,0.0);    
    C[0]     = a;
    C[nIntPoints-1] = a;
    solveTriDiagMatrix(nIntPoints,a,b,C,invD_C);
    
    outarray[0][0] = ( inarray[0][0] - a * (invD_f[0] + invD_f[nIntPoints-1]) ) / (b - a*(invD_C[0] + invD_C[nIntPoints-1]) );
    outarray[0][m_nPoints-1] = outarray[0][0];

    //Solve for the internal points
    for(int i = 1; i < m_nPoints-1; i++)
    {
        outarray[0][i] = invD_f[i-1] - invD_C[i-1] * outarray[0][0];
    }       
    
}

void OneDfinDiffAdvDiffSolver::EvaluateAdvectionDc(const Array<OneD, const  Array<OneD, double> >& inarray,                                                     
                                                         Array<OneD,        Array<OneD, double> >& outarray,
                                                   const NekDouble time) const
{   
    boost::ignore_unused(time);
    boost::ignore_unused(inarray);     
    for(int i = 0; i < m_nPoints; i++)
    {
        outarray[0][i] = 0.0;
    }        
    
}

void OneDfinDiffAdvDiffSolver::EvaluateAdvectionImex(const Array<OneD, const  Array<OneD, double> >& inarray,                                                     
                                                           Array<OneD,        Array<OneD, double> >& outarray,
                                                     const NekDouble time) const
{   
    boost::ignore_unused(time); 
    EvaluateAdvectionSp(inarray[0], null, outarray[0]);
}

void OneDfinDiffAdvDiffSolver::EvaluateAdvectionSp(const Array<OneD, const double>& predictor,
                                                   const Array<OneD, const double>& corrector,                                                                                                      
                                                         Array<OneD, double>& outarray) const
{   
    // second-order central differences        
    outarray[0]           = - m_U * ( (predictor[1] - predictor[m_nPoints-2]) - (corrector[1] - corrector[m_nPoints-2]) ) / ( 2.0 * m_dx );
    outarray[m_nPoints-1] = outarray[0];

    for(int i = 1; i < m_nPoints-1; i++)
    {
        outarray[i] = - m_U * ( (predictor[i+1] - predictor[i-1]) - (corrector[i+1] - corrector[i-1]) ) / ( 2.0 * m_dx );
    }    
        
}


void OneDfinDiffAdvDiffSolver::solveTriDiagMatrix (int n, double a, double b,
                                                   const Array<OneD, const double>& inarray,
                                                         Array<OneD,       double>& outarray) const
{
    // Implementation of the Thomas algorithm for Tridiagonal systems
    Array<OneD, double> cprime(n);
    Array<OneD, double> dprime(n);

    cprime[0] = a / b;
    dprime[0] = inarray[0] / b;

    for(int i = 1; i < n; i++)
    {
        double id = 1.0 / (b - cprime[i-1] * a);
        if(i <= n - 2) 
        {
            cprime[i] = a * id;
        }
        dprime[i] = (inarray[i] - dprime[i-1] * a) * id;
    }

    outarray[n-1] = dprime[n-1];
    for (int i = n - 2; i >= 0; i--)
    {
        outarray[i] = dprime[i] - cprime[i] * outarray[i+1];
    }

}


int OneDfinDiffAdvDiffSolver::factorial(int n) const
{
    return (n == 0 || n == 1) ? 1 : n * factorial(n-1);        
}


int OneDfinDiffAdvDiffSolver::GetNpoints() const
{
    return m_nPoints;
}


void OneDfinDiffAdvDiffSolver::EvaluateExactSolution(Array<OneD, double>& outarray,
                                                     const NekDouble time) const
{
    boost::ignore_unused(time);
    double x;
    for(int i = 0; i < m_nPoints; i++)
    {
        x = m_x0 + i*m_dx;
        outarray[i] = exp(-m_D * 2.0 * 2.0 * M_PI * M_PI * m_wavenumber*m_wavenumber*time) *
                      sin( 2.0 * m_wavenumber * M_PI * (x - m_U*time) );
    }
}


double OneDfinDiffAdvDiffSolver::EvaluateL2Error(const Array<OneD, const double>& approx,
                                                 const Array<OneD, const double>& exact) const
{
    //Computing the Relative L2-Error = sum_{i=0}^{m_nPoints} |approx(i) - exact(i)|^2/|exact(i)|^2
    double a = 0.0;
    double b = 0.0;

    for(int i = 0; i < m_nPoints; i++)
    {
        a += ( approx[i] - exact[i] ) * ( approx[i] - exact[i] );
        b += exact[i] * exact[i];
    }
    //Because only one time frame (at t=1) is used to compute L2(L2)-norm
    return sqrt(( m_tend - m_t0 ) / 1 * a / b);
    //return sqrt(( m_tend - m_t0 ) / 1 * a);
}


void OneDfinDiffAdvDiffSolver::AppendOutput(ofstream& outfile,
                                            const Array<OneD, const double>& approx,
                                            const Array<OneD, const double>& exact) const
{
    //Writing the output "x(j)"  "approx(j)"  "exact(j)" for each TimeStep(i)
    for(int i = 0; i < m_nPoints; i++)
    {
        outfile << scientific
                << setw (17)
                << setprecision(10)
                << m_x0 + i*m_dx
                << "  "
                << approx[i]
                << "  "
                << exact[i]                
                << endl;
    }
    outfile << endl << endl;
}


void OneDfinDiffAdvDiffSolver::GenerateGnuplotScript(NekDouble levelref) const
{
    //Writing out gnuplot formatting file (*.p) to plot the generated output (*.dat) 
    ofstream outfile;
    outfile.open("OneDfinDiffAdvDiffSolverOutput.p");

    outfile << "# Gnuplot script file" << endl;
    outfile << "set   autoscale" << endl;
    outfile << "unset log" << endl;
    outfile << "unset label" << endl;
    outfile << "set xtic auto" << endl;
    outfile << "set ytic auto" << endl;
    outfile << "set title \"Finite Difference Solution to the 1D advection-diffusion equation\"" << endl;
    outfile << "set xlabel \"x\"" << endl;
    outfile << "set ylabel \"u\"" << endl;
    outfile << "set xr [" << m_x0 << ":" << m_xend << "]" << endl;
    outfile << "set yr [-1.0:1.0]" << endl;

    double t;
    for(int i=0; i <= m_nTimeSteps * pow(2,levelref); i++)
    {
        t = m_t0+i*m_dt/pow(2,levelref);
        outfile << "plot    \"OneDfinDiffAdvDiffSolverOutput.dat\" ";
        outfile << "using 1:2 index ";
        outfile << i << " title 'Finite Difference Ref Solution (t=" << t << ")' with linespoints , ";
        outfile << "\"OneDfinDiffAdvDiffSolverOutput.dat\" ";
        outfile << "using 1:3 index ";
        outfile << i << " title 'Exact Solution (t=" << t << ")' with linespoints" << endl;
        outfile << "pause " << 4.0/m_nTimeSteps << endl;
    }

    outfile.close();
}


double OneDfinDiffAdvDiffSolver::GetInitialTime() const
{
    return m_t0;
}


double OneDfinDiffAdvDiffSolver::GetFinalTime() const
{
    return m_tend;
}


double OneDfinDiffAdvDiffSolver::GetTimeStep() const
{
    return m_dt;
}


void OneDfinDiffAdvDiffSolver::AsympExpansion(const Array<OneD, const double>& inarray0,
                                              const Array<OneD, const double>& inarray1,
                                              const Array<OneD, const double>& inarray2,      
                                                    Array<OneD, double>& outarray,
                                              const NekDouble dt) const 
{  
    for(int i = 0; i < m_nPoints; i++)
    {
        outarray[i] = inarray0[i] + inarray1[i] * dt + inarray2[i] * dt * dt;        
    }
}
