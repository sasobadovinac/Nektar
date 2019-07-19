/////////////////////////////////////////////////////////////////////////////
//
// File GalerkinBoltzmann.cpp
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
// Description: Unsteady linear advection solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <GalerkinBoltzmannSolver/EquationSystems/GalerkinBoltzmann.h>

using namespace std;

namespace Nektar
{
    string GalerkinBoltzmann::className =
        SolverUtils::GetEquationSystemFactory().
        RegisterCreatorFunction("GalerkinBoltzmann",
                                GalerkinBoltzmann::create,
                                "Galerkin Boltzmann equation.");

    GalerkinBoltzmann::GalerkinBoltzmann(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          AdvectionSystem(pSession, pGraph)
    {
        m_planeNumber = 0;
    }

    /**
     * @brief Initialisation object for the Galerkin Boltzmann equation.
     */
    void GalerkinBoltzmann::v_InitObject()
    {
        // Call to the initialisation object of UnsteadySystem
        AdvectionSystem::v_InitObject();

        ASSERTL0(m_fields.num_elements() >= 5,
                 "Expected there to be five variables in Galerkin "
                 "Boltzmann solver (in 2-dimensions)");
        
        ASSERTL0(m_session->DefinesParameter("RT"),"Must define parameter RT");
        m_sqrtRT = sqrt(m_session->GetParameter("RT"));

        ASSERTL0(m_session->DefinesParameter("tau"),"Must define parameter tau");
        m_tau    = m_session->GetParameter("tau");
        m_traceNormals =  Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        int nTracePointsTot = m_fields[0]->GetTrace()->GetTotPoints();
        //int nTracePointsTot = m_fields[0]->GetTraceNpoints();
        for(int i = 0; i < m_spacedim; ++i)
        {
            m_traceNormals[i] = Array<OneD, NekDouble> (nTracePointsTot);
        }
	m_fields[0]->GetTrace()->GetNormals(m_traceNormals);
        // Type of advection class to be used
        switch(m_projectionType)
        {
            // Discontinuous field
            case MultiRegions::eDiscontinuous:
            {
                // Do not forwards transform initial condition
                m_homoInitialFwd = false;
                break;
            }
            default:
            {
                ASSERTL0(false, "Unsupported projection type.");
                break;
            }
        }

        // If explicit it computes RHS and PROJECTION for the time integration
        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs     (&GalerkinBoltzmann::DoOdeRhs,        this);
            m_ode.DefineProjection (&GalerkinBoltzmann::DoOdeProjection, this);
        }
        // Otherwise it gives an error (no implicit integration)
        else
        {
            ASSERTL0(false, "Implicit unsteady Advection not set up.");
        }
    }

    /**
     * @brief Galerkin Boltzmann destructor.
     */
    GalerkinBoltzmann::~GalerkinBoltzmann()
    {
    }

    /**
     * @brief Compute the right-hand side for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void GalerkinBoltzmann::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
              Array<OneD,        Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        int ndim       = m_spacedim;

        ASSERTL0(ndim == 2,"Method is not set up for 3D (or 1D) ");

        int nFields = inarray.num_elements();
        int ncoeffs = m_fields[0]->GetNcoeffs();

        // Calculate inner product with respect to derivative of basis of volume flux
        Array<OneD, Array<OneD, NekDouble> > volflux(nFields);
        for(int i = 0; i < nFields; ++i)
        {
            volflux[i] = Array<OneD, NekDouble> (ncoeffs);
        }
        
        EvaluateIProductWRTDerivBaseVolFlux(inarray,volflux);
        // Calculate the normal numerical flux on trace boudaries
        // Add upwind normal flux component
        int nTracePointsTot = m_fields[0]->GetTrace()->GetTotPoints();

        // Store forwards/backwards space along trace space
        Array<OneD, Array<OneD, NekDouble> > Fwd    (nFields);
        Array<OneD, Array<OneD, NekDouble> > Bwd    (nFields);
        Array<OneD, Array<OneD, NekDouble> > numflux(nFields);
        for(int i = 0; i < nFields; ++i)
        {
            Fwd[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
            Bwd[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
            numflux[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
            m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
        }

        EvaluateNormalNumFlux(Fwd,Bwd,numflux);
        // Project back to physical space
        for(int i = 0; i < nFields; ++i)
        {
            Vmath::Neg                        (ncoeffs,    volflux[i], 1);
            m_fields[i]->AddTraceIntegral     (numflux[i], volflux[i]);
            m_fields[i]->MultiplyByElmtInvMass(volflux[i], volflux[i]);
            m_fields[i]->BwdTrans             (volflux[i], outarray[i]);
        }
        //-------------------------------------------------------
        // Negate the RHS
        int nq = m_fields[0]->GetTotPoints();
        for (int i = 0; i < nFields; ++i)
        {
            Vmath::Neg(nq, outarray[i], 1);
            //Vmath::Zero(nq, outarray[i], 1);
        }
        AddSourceTerms(inarray,outarray);
       
    }

    /**
     * @brief Compute the projection for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void GalerkinBoltzmann::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD,       Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        // Counter variable
        int i;

        // Number of fields (variables of the problem)
        int nVariables = inarray.num_elements();


        // Number of quadrature points
        int nQuadraturePts = GetNpoints();
        
        // Just copy over array
        for(i = 0; i < nVariables; ++i)
        {
            Vmath::Vcopy(nQuadraturePts, inarray[i], 1, outarray[i], 1);
        } 

        // Set the boundary conditions
        SetBoundaryConditions(time);
   }

    /**
     * @brief Return the flux vector for the linear advection equation.
     *
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void GalerkinBoltzmann::EvaluateIProductWRTDerivBaseVolFlux(
        const Array<OneD, const  Array<OneD, NekDouble> >&physfield,
        Array<OneD,        Array<OneD, NekDouble> >&volflux)

    {

        NekDouble sqrt_2RT = sqrt(2)*m_sqrtRT;
        int nq = physfield[0].num_elements();

        Array<OneD, Array<OneD, NekDouble> > flux(2);
        flux[0] = Array<OneD, NekDouble>(nq);
        flux[1] = Array<OneD, NekDouble>(nq);

        //  Volume Flux for first component = sqrt(RT) [a1, a2]
        Vmath::Smul(nq, m_sqrtRT, physfield[1], 1, flux[0], 1);
        Vmath::Smul(nq, m_sqrtRT, physfield[2], 1, flux[1], 1);
        m_fields[0]->IProductWRTDerivBase(flux,volflux[0]);
            
        //  Volume Flux for second component =  sqrt(RT)[a0 + sqrt(2) a4, a3]
        Vmath::Smul(nq, m_sqrtRT, physfield[0], 1, flux[0], 1);
        Vmath::Svtvp(nq,sqrt_2RT,physfield[4],1,flux[0],1,flux[0],1);
        Vmath::Smul(nq, m_sqrtRT, physfield[3], 1, flux[1], 1);
        m_fields[1]->IProductWRTDerivBase(flux,volflux[1]);

        //  Volume Flux for third component = sqrt(RT)[a3, a0 + sqrt(2) a5]
        Vmath::Smul(nq, m_sqrtRT, physfield[3], 1, flux[0], 1);
        Vmath::Smul(nq, m_sqrtRT, physfield[0], 1, flux[1], 1);
        Vmath::Svtvp(nq,sqrt_2RT,physfield[5],1,flux[1],1,flux[1],1);
        m_fields[2]->IProductWRTDerivBase(flux,volflux[2]);

        //  Volume Flux for fourth component = sqrt(RT)[a2, a1]
        Vmath::Smul(nq, m_sqrtRT, physfield[2], 1, flux[0], 1);
        Vmath::Smul(nq, m_sqrtRT, physfield[1], 1, flux[1], 1);
        m_fields[3]->IProductWRTDerivBase(flux,volflux[3]);

        //  Volume Flux for fifth component = sqrt(RT)[sqrt(2) a1, 0]
        Vmath::Smul(nq, sqrt_2RT, physfield[1], 1, flux[0], 1);
        Vmath::Zero(nq,  flux[1], 1);
        m_fields[4]->IProductWRTDerivBase(flux,volflux[4]);

        //  Volume Flux for sixth component = sqrt(RT)[0, sqrt(2) a2]
        Vmath::Zero(nq,  flux[0], 1);
        Vmath::Smul(nq, sqrt_2RT, physfield[2], 1, flux[1], 1);
        m_fields[5]->IProductWRTDerivBase(flux,volflux[5]);

    }


    void GalerkinBoltzmann::EvaluateNormalNumFlux(const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
                                                  const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                                                  Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int nTracePointsTot = m_fields[0]->GetTrace()->GetTotPoints();
	NekDouble Interior[6][6]; // matrix V*Lambda^+*V^(-1)  along the normal
        NekDouble Exterior[6][6]; // matrix V*Lambda^-*V^(-1)  opposite to the normal
        // get hold of normals

        NekDouble sqrt2 = sqrt(2);
        NekDouble sqrt3 = sqrt(3);
        NekDouble sqrt6 = sqrt(6);

        for(int i = 0 ; i < nTracePointsTot; ++i)
        {
            NekDouble nx = m_traceNormals[0][i];
            NekDouble ny = m_traceNormals[1][i];
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

       
            // Fwd exter
            // Bwd inter

            for(int j = 0; j < 6; ++j)
            {
                 
		numflux[j][i] = Interior[j][0]*Bwd[0][i]+Interior[j][1]*Bwd[1][i]+Interior[j][2]*Bwd[2][i]+Interior[j][3]*Bwd[3][i]+Interior[j][4]*Bwd[4][i]+Interior[j][5]*Bwd[5][i]+
           		        Exterior[j][0]*Fwd[0][i]+Exterior[j][1]*Fwd[1][i]+Exterior[j][2]*Fwd[2][i]+Exterior[j][3]*Fwd[3][i]+Exterior[j][4]*Fwd[4][i]+Exterior[j][5]*Fwd[5][i];
                numflux[j][i] *= (m_sqrtRT);
            }
            //cout << nx << "\t" << ny <<"\t" <<endl; 
        }
    }

    void GalerkinBoltzmann::AddSourceTerms( const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                                                  Array<OneD,        Array<OneD, NekDouble> >&outarray)
    {

        int nq = m_fields[0]->GetTotPoints();
        
        //-------------------------------------------------
        // Add "source terms"
        // Input and output in physical space so using a collocation projection 
        //-------------------------------------------------
        Array<OneD, NekDouble> force(nq);
        
        NekDouble invsqrt_2 =  1.0/sqrt(2);
        NekDouble invtau    = -1.0/m_tau;
	/*
	Vmath::Zero(nq, outarray[0], 1);
	Vmath::Zero(nq, outarray[1], 1);
	Vmath::Zero(nq, outarray[2], 1);
	Vmath::Zero(nq, outarray[3], 1);
	Vmath::Zero(nq, outarray[4], 1);
	Vmath::Zero(nq, outarray[5], 1);
	*/
        // outarray[3] = -1/tau (a[3] - a[1] a[2]/a[0])
        Vmath::Vmul(nq,inarray[1],1,inarray[2],1,force,1);
        Vmath::Vdiv(nq,force,1,inarray[0],1,force,1);
        Vmath::Vsub(nq,inarray[3],1,force,1,force,1);
        Vmath::Smul(nq,invtau,force,1,force,1);
        Vmath::Vadd(nq,outarray[3],1,force,1,outarray[3],1);
        // outarray[4] = -1/tau (a[4] - 1/sqrt(2) x a[1] a[1]/a[0])
        Vmath::Vmul(nq,inarray[1],1,inarray[1],1,force,1);
        Vmath::Vdiv(nq,force,1,inarray[0],1,force,1);
        Vmath::Smul(nq,invsqrt_2,force,1,force,1);
        Vmath::Vsub(nq,inarray[4],1,force,1,force,1);
        Vmath::Smul(nq,invtau,force,1,force,1);
        Vmath::Vadd(nq,outarray[4],1,force,1,outarray[4],1);
        // outarray[5] = -1/tau (a[5] - 1/sqrt(2) x a[2] a[2]/a[0])
        Vmath::Vmul(nq,inarray[2],1,inarray[2],1,force,1);
        Vmath::Vdiv(nq,force,1,inarray[0],1,force,1);
        Vmath::Smul(nq,invsqrt_2,force,1,force,1);
        Vmath::Vsub(nq,inarray[5],1,force,1,force,1);
        Vmath::Smul(nq,invtau,force,1,force,1);
        Vmath::Vadd(nq,outarray[5],1,force,1,outarray[5],1);
    }
    
    void GalerkinBoltzmann::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        AdvectionSystem::v_GenerateSummary(s);
    }
}
