///////////////////////////////////////////////////////////////////////////////
//
// File: Extrapolate.h
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

#ifndef NEKTAR_SOLVERS_EXTRAPOLATE_H
#define NEKTAR_SOLVERS_EXTRAPOLATE_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationWrapper.h>
#include <IncNavierStokesSolver/AdvectionTerms/AdvectionTerm.h>

namespace Nektar
{
	
	// Velocity correction scheme coefficient required for extrapolation.
	
    static NekDouble StifflyStable_Betaq_Coeffs[][3] = {{ 1.0,  0.0, 0.0},
		{ 2.0, -1.0, 0.0},
		{ 3.0, -3.0, 1.0}};
    
    static NekDouble StifflyStable_Alpha_Coeffs[][3] = {{ 1.0,  0.0, 0.0},
		{ 2.0, -0.5, 0.0},
		{ 3.0, -1.5, 1.0/3.0}};
    
    static NekDouble StifflyStable_Gamma0_Coeffs[3]  = {1.0,  1.5, 11.0/6.0};
	
	struct HBCInfo
    {
        int m_globalElmtID;  // elements ID in the global ordering
        
        int m_ptsInElmt;     // number of physical points of the element
        
        int m_physOffset;    // elmenent physical offset in the global list
        
        int m_bndElmtOffset; // element offset in the boundary expansion
        
        int m_elmtTraceID;   // trace ID on the element
        
        int m_bndryElmtID;   // pressure boundary condition ID
        
        int m_assPhysOffset; // associated elments physical offset (k and k_c
		// are the real and the complex plane)
        
        int m_coeffOffset;   // coefficients offset used to locate the
		// acceleration term in the general m_pressureHBC
    };
	
	
    // Forward declaration
    class Extrapolate;

    typedef boost::shared_ptr<Extrapolate> ExtrapolateSharedPtr;
    
    typedef LibUtilities::NekFactory< std::string, Extrapolate,
        const LibUtilities::SessionReaderSharedPtr& ,
        Array<OneD, MultiRegions::ExpListSharedPtr>& ,
        Array<OneD, int> ,
        AdvectionTermSharedPtr > ExtrapolateFactory; 

    ExtrapolateFactory& GetExtrapolateFactory();


    class Extrapolate
    {
    public:
        Extrapolate(        
            const LibUtilities::SessionReaderSharedPtr pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
            Array<OneD, int> pVel,
            AdvectionTermSharedPtr advObject);
        
        virtual ~Extrapolate();
		
		void GenerateHOPBCMap();

        inline void SubSteppingTimeIntegration(
            int intMethod);

        inline void SubStepSaveFields(
            const int nstep);
        
        inline void SubStepSetPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            const NekDouble Aii_DT,
            NekDouble kinvis);

        inline void SubStepAdvance(
            LibUtilities::TimeIntegrationSolutionSharedPtr integrationSoln, 
            const int nstep, 
            NekDouble time);
		
        inline void MountHOPBCs(
            int HBCdata, 
            NekDouble kinvis, 
            Array<OneD, NekDouble> &Q, 
            Array<OneD, const NekDouble> &Advection);

        void EvaluatePressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &fields,
            const Array<OneD, const Array<OneD, NekDouble> >  &N,
            const int kinvis);
        
    protected:
        virtual void v_SubSteppingTimeIntegration(
            int intMethod)=0;

        virtual void v_SubStepSaveFields(
            const int nstep)=0;

        virtual void v_SubStepSetPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            const NekDouble Aii_DT,
            NekDouble kinvis)=0;

        virtual void v_SubStepAdvance(
            LibUtilities::TimeIntegrationSolutionSharedPtr integrationSoln, 
            const int nstep, 
            NekDouble time)=0;

        virtual void v_MountHOPBCs(
            int HBCdata, 
            NekDouble kinvis, 
            Array<OneD, NekDouble> &Q, 
            Array<OneD, const NekDouble> &Advection)=0;
        
        void CalcPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &fields,
            const Array<OneD, const Array<OneD, NekDouble> >  &N,
            const int kinvis);
        
        void RollOver(
            Array<OneD, Array<OneD, NekDouble> > &input);
		
        void CurlCurl(
            Array<OneD, Array<OneD, const NekDouble> > &Vel,
            Array<OneD, Array<OneD, NekDouble> > &Q,
            const int j);
        
        LibUtilities::SessionReaderSharedPtr m_session;

        LibUtilities::CommSharedPtr m_comm;

        Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;

        Array<OneD, int> m_velocity;

        AdvectionTermSharedPtr m_advObject;

        Array<OneD, Array<OneD, NekDouble> > m_previousVelFields;
		
        /// Curl-curl dimensionality
        int m_curl_dim;
        
        /// bounday dimensionality
        int m_bnd_dim;
		
        /// pressure boundary conditions container
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr> m_PBndConds;
		
        /// pressure boundary conditions expansion container
        Array<OneD, MultiRegions::ExpListSharedPtr>  m_PBndExp;
        
        /// number of times the high-order pressure BCs have been called
        int m_pressureCalls;
        
        /// Maximum points used in pressure BC evaluation
        int m_pressureBCsMaxPts;
		
        /// Maximum points used in pressure BC evaluation
        int m_intSteps;

        NekDouble m_timestep;

        /// Flag to determine if single homogeneous mode is used.
        bool m_SingleMode;
        /// Flag to determine if half homogeneous mode is used.
        bool m_HalfMode;
        /// Flag to determine if use multiple homogenenous modes are used.
        bool m_MultipleModes;

        NekDouble m_LhomZ;  ///< physical length in Z direction (if homogeneous)
        
        int m_npointsX;     ///< number of points in X direction (if homogeneous)
        int m_npointsY;     ///< number of points in Y direction (if homogeneous)
        int m_npointsZ;     ///< number of points in Z direction (if homogeneous)


		
        /// Id of element to which pressure  boundary condition belongs
        Array<OneD, int> m_pressureBCtoElmtID;
        
        /// Id of edge (2D) or face (3D) to which pressure boundary condition belongs
        Array<OneD, int> m_pressureBCtoTraceID;
        
        /// Storage for current and previous levels of high order pressure boundary conditions.
        Array<OneD, Array<OneD, NekDouble> >  m_pressureHBCs;
		
        /// Storage for current and previous levels of the acceleration term.
        Array<OneD, Array<OneD, NekDouble> >  m_acceleration;
        
        /// data structure to old all the information regarding High order pressure BCs
        Array<OneD, HBCInfo > m_HBCdata;
		
        /// general standard element used to deal with HOPBC calculations
        StdRegions::StdExpansionSharedPtr m_elmt;
        
        /// wave number 2 pi k /Lz
        Array<OneD, NekDouble>  m_wavenumber;
        
        /// minus Square of wavenumber
        Array<OneD, NekDouble>  m_negWavenumberSq;

    private:
        static std::string def;
    };

    /**
     *
     */
    inline void Extrapolate::SubSteppingTimeIntegration(
        int intMethod)
    {
        v_SubSteppingTimeIntegration(intMethod);
    }
    
    /**
     *
     */
    inline void Extrapolate::SubStepSaveFields(
        const int nstep)
    {
        v_SubStepSaveFields(nstep);
    }

    /**
     *
     */
    inline void Extrapolate::SubStepSetPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            const NekDouble Aii_DT,
            NekDouble kinvis)
    {
        v_SubStepSetPressureBCs(inarray,Aii_DT,kinvis);
    }
    
    /**
     *
     */
    inline void Extrapolate::SubStepAdvance(
        LibUtilities::TimeIntegrationSolutionSharedPtr integrationSoln, 
        const int nstep, 
        NekDouble time)
    {
        v_SubStepAdvance(integrationSoln,nstep, time);
    }
	
    /**
     *
     */
    inline void Extrapolate::MountHOPBCs(
        int HBCdata, 
        NekDouble kinvis, 
        Array<OneD, NekDouble> &Q, 
        Array<OneD, const NekDouble> &Advection)
    {
        v_MountHOPBCs(HBCdata,kinvis,Q,Advection);
    }
}

#endif

