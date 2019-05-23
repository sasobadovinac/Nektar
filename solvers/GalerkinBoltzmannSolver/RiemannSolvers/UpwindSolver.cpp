///////////////////////////////////////////////////////////////////////////////
//
// File: UpwindSolver.cpp
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
// Description: HLL Riemann solver for the Linear Shallow Water Equations.
//              Only valid for constant depth
//
///////////////////////////////////////////////////////////////////////////////

#include <GalerkinBoltzmannSolver/RiemannSolvers/UpwindSolver.h>

namespace Nektar
{
    std::string UpwindSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "BoltzmannUpwind",
            UpwindSolver::create,
            "Upwind Riemann solver");

    UpwindSolver::UpwindSolver(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : BoltzmannSolver(pSession)
    {

    }

    /**
     * @brief Upwind Riemann solver for the  Boltzmann Equations
     * @param nx     normal direction x-component
     * @param ny     normal direction y-component
     * @param nz     normal direction z-component
     * @param etaL   Free surface elevation left state.
     * @param etaR   Free surface elevation right state.  
     * @param uL     x-velocity  left state.  
     * @param uR     x-velocity  right state.  
     * @param vL     y-velocity  left state.  
     * @param vR     y-velocity  right state. 
     * @param dL     still water depth component left state.  
     * @param dR     still water depth component right state. 
     * @param uf     Computed Riemann flux for x-momentum component 
     * @param vf     Computed Riemann flux for y-momentum component 
     */
    void UpwindSolver::v_PointSolve(
         const int nDim, const NekDouble sqrt_rt,
         NekDouble   nx, NekDouble   ny, NekDouble nz,
         NekDouble  a0L, NekDouble  a1L, NekDouble a2L,
         NekDouble  a3L, NekDouble  a4L, NekDouble a5L,
         NekDouble  a0R, NekDouble  a1R, NekDouble a2R,
         NekDouble  a3R, NekDouble  a4R, NekDouble a5R,
         NekDouble  &f0, NekDouble  &f1, NekDouble &f2,
         NekDouble  &f3, NekDouble  &f4, NekDouble &f5)
    { 
	int i; 
	ASSERTL0(m_spacedim == 2, "Space Dimension should be 2."); // Currently support   
        NekDouble Interior[6][6]; // matrix V*Lambda^+*V^(-1)  
        NekDouble Exterior[6][6]; // matrix V*Lambda^-*V^(-1) 
 
        // Interior matrix
        //[ -3^(1/2)/6, nx/2, ny/2, -(3^(1/2)*nx*ny)/3, -(6^(1/2)*nx^2)/6, -(6^(1/2)*ny^2)/6]
        Interior[0][0] = -sqrt(3.0)/6.;          
        Interior[0][1] =  nx/2.;              
        Interior[0][2] =  ny/2.; 
        Interior[0][3] = -sqrt(3.0)*nx*ny/3.0;   
        Interior[0][4] = -sqrt(6.)*nx*nx/6.0; 
        Interior[0][5] = -sqrt(6.)*ny*ny/6.0;
        
        //[ nx/2, -(3^(1/2)*nx^2)/2 - ny^2/2, -(nx*ny*(3^(1/2) - 1))/2, ny/6, (2^(1/2)*nx)/6, 0]
        Interior[1][0] =  nx/2.;  
        Interior[1][1] = -sqrt(3.)*nx*nx/2. - ny*ny/2;  
        Interior[1][2] = -(nx*ny*(sqrt(3.) - 1))/2; 
        Interior[1][3] =  ny/6.; 
        Interior[1][4] =  sqrt(2.)*nx/6.0;               
        Interior[1][5] =  0;

        //[ ny/2, -(nx*ny*(3^(1/2) - 1))/2, -nx^2/2 - (3^(1/2)*ny^2)/2, nx/6, 0, (2^(1/2)*ny)/6]
        Interior[2][0] = ny/2.;  
        Interior[2][1] = -(nx*ny*(sqrt(3.) - 1))/2.;  
        Interior[2][2] = -nx*nx/2.-(sqrt(3.)*ny*ny)/2.;    
        Interior[2][3] = nx/6.0; 
        Interior[2][4] = 0;                           
        Interior[2][5] = sqrt(2.)*ny/6.0;
        // [ -(3^(1/2)*nx*ny)/3, ny/6, nx/6, - (nx^2 - ny^2)^2/2 - (2*3^(1/2)*nx^2*ny^2)/3, -(2^(1/2)*nx*ny*(2*3^(1/2)*nx^2 - 3*nx^2 + 3*ny^2))/6, -(2^(1/2)*nx*ny*(2*3^(1/2)*ny^2 + 3*nx^2 - 3*ny^2))/6]
        Interior[3][0] = -(sqrt(3.)*nx*ny)/3.;  
        Interior[3][1] = ny/6.;  
        Interior[3][2] = nx/6.;   
        Interior[3][3] = -(nx*nx-ny*ny)*(nx*nx-ny*ny)/2.-(2*sqrt(3.)*nx*nx*ny*ny)/3.; 
        Interior[3][4] = -(sqrt(2.)*nx*ny*(2*sqrt(3.)*nx*nx-3*nx*nx+3*ny*ny))/6.;
        Interior[3][5] = -(sqrt(2.)*nx*ny*(2*sqrt(3.)*ny*ny+3*nx*nx-3*ny*ny))/6.;
        // [ -(6^(1/2)*nx^2)/6, (2^(1/2)*nx)/6, 0, -(2^(1/2)*nx*ny*(2*3^(1/2)*nx^2 - 3*nx^2 + 3*ny^2))/6, - (3^(1/2)*nx^4)/3 - nx^2*ny^2, -(nx^2*ny^2*(3^(1/2) - 3))/3]
        Interior[4][0] = -(sqrt(6.)*nx*nx)/6.;
        Interior[4][1] =  (sqrt(2.)*nx)/6.;
        Interior[4][2] =  0;
        Interior[4][3] = -(sqrt(2.)*nx*ny*(2*sqrt(3.)*nx*nx-3*nx*nx+3*ny*ny))/6.;
        Interior[4][4] = -(sqrt(3.)*nx*nx*nx*nx)/3.-nx*nx*ny*ny;
        Interior[4][5] = -(nx*nx*ny*ny*(sqrt(3.)-3))/3.;
        // [ -(6^(1/2)*ny^2)/6, 0, (2^(1/2)*ny)/6, -(2^(1/2)*nx*ny*(2*3^(1/2)*ny^2 + 3*nx^2 - 3*ny^2))/6, -(nx^2*ny^2*(3^(1/2) - 3))/3, - nx^2*ny^2 - (3^(1/2)*ny^4)/3]
        Interior[5][0] = -(sqrt(6.)*ny*ny)/6.;
        Interior[5][1] =  0;
        Interior[5][2] =  (sqrt(2.)*ny)/6.;
        Interior[5][3] = -(sqrt(2.)*nx*ny*(2*sqrt(3.)*ny*ny+3*nx*nx-3*ny*ny))/6.;
        Interior[5][4] = -(nx*nx*ny*ny*(sqrt(3.)-3.))/3.;
        Interior[5][5] = - nx*nx*ny*ny-(sqrt(3)*ny*ny*ny*ny)/3.;

        // Exterior matrix
        // [ 3^(1/2)/6, nx/2, ny/2, (3^(1/2)*nx*ny)/3, (6^(1/2)*nx^2)/6, (6^(1/2)*ny^2)/6]
        Exterior[0][0] = sqrt(3.)/6.;   
        Exterior[0][1] = nx/2.; 
	Exterior[0][2] = ny/2.;
	Exterior[0][3] = sqrt(3.)*nx*ny/3.0;
	Exterior[0][4] = sqrt(6.)*nx*nx/6.0;
	Exterior[0][5] = sqrt(6.)*ny*ny/6.0;
        // [ nx/2, (3^(1/2)*nx^2)/2 + ny^2/2, (nx*ny*(3^(1/2) - 1))/2, ny/6, (2^(1/2)*nx)/6, 0]
        Exterior[1][0] = nx/2.;
	Exterior[1][1] = sqrt(3.)*nx*nx/2.0+ny*ny/2;
	Exterior[1][2] = (nx*ny*(sqrt(3.) - 1))/2;
	Exterior[1][3] = ny/6.0;
	Exterior[1][4] = sqrt(2.)*nx/6.0;
	Exterior[1][5] = 0;
        // [ ny/2, (nx*ny*(3^(1/2) - 1))/2, nx^2/2 + (3^(1/2)*ny^2)/2, nx/6, 0, (2^(1/2)*ny)/6]
        Exterior[2][0] = ny/2.;
	Exterior[2][1] = (nx*ny*(sqrt(3.)-1))/2.;
	Exterior[2][2] = nx*nx/2.-(sqrt(3.)*ny*ny)/2.;
	Exterior[2][3] = nx/6.;
	Exterior[2][4] = 0;
	Exterior[2][5] = sqrt(2.)*ny/6.;
        // [ (3^(1/2)*nx*ny)/3, ny/6, nx/6, (nx^2 - ny^2)^2/2 + (2*3^(1/2)*nx^2*ny^2)/3, (2^(1/2)*nx*ny*(2*3^(1/2)*nx^2 - 3*nx^2 + 3*ny^2))/6, (2^(1/2)*nx*ny*(2*3^(1/2)*ny^2 + 3*nx^2 - 3*ny^2))/6]
 	Exterior[3][0] = (sqrt(3.)*nx*ny)/3.; 
	Exterior[3][1] = ny/6.;
	Exterior[3][2] = nx/6.;
	Exterior[3][3] = (nx*nx-ny*ny)*(nx*nx-ny*ny)/2.-(2*sqrt(3.)*nx*nx*ny*ny)/3.;
	Exterior[3][4] = (sqrt(2.)*nx*ny*(2*sqrt(3.)*nx*nx-3*nx*nx+3*ny*ny))/6.;
	Exterior[3][5] = (sqrt(2.)*nx*ny*(2*sqrt(3.)*ny*ny+3*nx*nx-3*ny*ny))/6.;
        // [ (6^(1/2)*nx^2)/6, (2^(1/2)*nx)/6, 0, (2^(1/2)*nx*ny*(2*3^(1/2)*nx^2 - 3*nx^2 + 3*ny^2))/6, (3^(1/2)*nx^4)/3 + nx^2*ny^2, (nx^2*ny^2*(3^(1/2) - 3))/3]
        Exterior[4][0] = (sqrt(6.)*nx*nx)/6.;
	Exterior[4][1] = (sqrt(2.)*nx)/6.;
	Exterior[4][2] = 0; 
	Exterior[4][3] = (sqrt(2.)*nx*ny*(2*sqrt(3.)*nx*nx-3*nx*nx+3*ny*ny))/6.;
	Exterior[4][4] = (sqrt(3.)*nx*nx*nx*nx)/3.+nx*nx*ny*ny;
	Exterior[4][5] = (nx*nx*ny*ny*(sqrt(3.)-3))/3.;
        // [ (6^(1/2)*ny^2)/6, 0, (2^(1/2)*ny)/6, (2^(1/2)*nx*ny*(2*3^(1/2)*ny^2 + 3*nx^2 - 3*ny^2))/6, (nx^2*ny^2*(3^(1/2) - 3))/3, nx^2*ny^2 + (3^(1/2)*ny^4)/3]
    Exterior[5][0] = (sqrt(6.)*ny*ny)/6.;
	Exterior[5][1] = 0;
	Exterior[5][2] = (sqrt(2.)*ny)/6.;
	Exterior[5][3] = (sqrt(2.)*nx*ny*(2*sqrt(3.)*ny*ny+3*nx*nx-3*ny*ny))/6.;
	Exterior[5][4] = (nx*nx*ny*ny*(sqrt(3.)-3.))/3.;
	Exterior[5][5] =  nx*nx*ny*ny+(sqrt(3)*ny*ny*ny*ny)/3.;
        // Splitting scheme
	i=0;
	f0=Interior[i][0]*a0L+Interior[i][1]*a1L+Interior[i][2]*a2L+Interior[i][3]*a3L+Interior[i][4]*a4L+Interior[i][5]*a5L+
           Exterior[i][0]*a0R+Exterior[i][1]*a1R+Exterior[i][2]*a2R+Exterior[i][3]*a3R+Exterior[i][4]*a4R+Exterior[i][5]*a5R;
	f0=-f0*sqrt_rt;
        i=1;
        f1=Interior[i][0]*a0L+Interior[i][1]*a1L+Interior[i][2]*a2L+Interior[i][3]*a3L+Interior[i][4]*a4L+Interior[i][5]*a5L+
           Exterior[i][0]*a0R+Exterior[i][1]*a1R+Exterior[i][2]*a2R+Exterior[i][3]*a3R+Exterior[i][4]*a4R+Exterior[i][5]*a5R;
	f1=-f1*sqrt_rt;
        i=2;
        f2=Interior[i][0]*a0L+Interior[i][1]*a1L+Interior[i][2]*a2L+Interior[i][3]*a3L+Interior[i][4]*a4L+Interior[i][5]*a5L+
           Exterior[i][0]*a0R+Exterior[i][1]*a1R+Exterior[i][2]*a2R+Exterior[i][3]*a3R+Exterior[i][4]*a4R+Exterior[i][5]*a5R;
	f2=-f2*sqrt_rt;
        i=3;
        f3=Interior[i][0]*a0L+Interior[i][1]*a1L+Interior[i][2]*a2L+Interior[i][3]*a3L+Interior[i][4]*a4L+Interior[i][5]*a5L+
           Exterior[i][0]*a0R+Exterior[i][1]*a1R+Exterior[i][2]*a2R+Exterior[i][3]*a3R+Exterior[i][4]*a4R+Exterior[i][5]*a5R;
	f3=-f3*sqrt_rt;
        i=4;
        f4=Interior[i][0]*a0L+Interior[i][1]*a1L+Interior[i][2]*a2L+Interior[i][3]*a3L+Interior[i][4]*a4L+Interior[i][5]*a5L+
           Exterior[i][0]*a0R+Exterior[i][1]*a1R+Exterior[i][2]*a2R+Exterior[i][3]*a3R+Exterior[i][4]*a4R+Exterior[i][5]*a5R;
	f4=-f4*sqrt_rt;
        i=5;
        f5=Interior[i][0]*a0L+Interior[i][1]*a1L+Interior[i][2]*a2L+Interior[i][3]*a3L+Interior[i][4]*a4L+Interior[i][5]*a5L+
           Exterior[i][0]*a0R+Exterior[i][1]*a1R+Exterior[i][2]*a2R+Exterior[i][3]*a3R+Exterior[i][4]*a4R+Exterior[i][5]*a5R;
	f5=-f5*sqrt_rt;
    }

}
