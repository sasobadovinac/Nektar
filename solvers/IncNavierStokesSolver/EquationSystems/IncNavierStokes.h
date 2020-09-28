///////////////////////////////////////////////////////////////////////////////
//
// File IncNavierStokes.h
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
// Description: Basic Advection Diffusion Reaction Field definition in two-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_INCNAVIERSTOKES_H
#define NEKTAR_SOLVERS_INCNAVIERSTOKES_H

#include <LibUtilities/TimeIntegration/TimeIntegrationWrapper.h>
#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/AdvectionSystem.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <IncNavierStokesSolver/EquationSystems/Extrapolate.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <complex>

#include <SolverUtils/Filters/FilterAeroForces.h>
#include <SolverUtils/EquationSystem.h>

namespace Nektar
{
    enum EquationType
    {
        eNoEquationType,
        eSteadyStokes,
        eSteadyOseen,
        eSteadyLinearisedNS,
        eUnsteadyStokes,
        eUnsteadyLinearisedNS,
        eUnsteadyNavierStokes,
        eSteadyNavierStokes,
        eEquationTypeSize
    };

    // Keep this consistent with the enums in EquationType.
    const std::string kEquationTypeStr[] =
    {
        "NoType",
        "SteadyStokes",
        "SteadyOseen",
        "SteadyLinearisedNS",
        "UnsteadyStokes",
        "UnsteadyLinearisedNS",
        "UnsteadyNavierStokes",
        "SteadyNavierStokes",
    };


    enum AdvectionForm
    {
        eNoAdvectionForm,
        eConvective,
        eNonConservative,
        eLinearised,
        eAdjoint,
        eSkewSymmetric,
        eNoAdvection,
        eAdvectionFormSize
    };

    // Keep this consistent with the enums in EquationType.
    const std::string kAdvectionFormStr[] =
    {
        "NoType",
        "Convective",
        "NonConservative",
        "Linearised",
        "Adjoint",
        "SkewSymmetric"
        "NoAdvection"
    };


    typedef std::complex<double> NekComplexDouble;

    struct WomersleyParams
    {
        WomersleyParams(int dim)
        {
            m_axisnormal = Array<OneD, NekDouble>(dim,0.0);
            m_axispoint  = Array<OneD, NekDouble>(dim,0.0);
        };

        virtual ~WomersleyParams()
        {};

        // Real and imaginary velocity comp. of wom
        std::vector<NekComplexDouble> m_wom_vel;

        // Womersley  BC constants
        NekDouble m_radius;
        NekDouble m_period;
        Array<OneD, NekDouble> m_axisnormal;
        // currently this needs to be the point in the middle of the
        // axis but should be generalised to be any point on the axis
        Array<OneD, NekDouble> m_axispoint;

        // poiseuille flow and fourier coefficients
        Array<OneD, Array<OneD, NekDouble> > m_poiseuille;
        Array<OneD, Array<OneD, Array<OneD, NekComplexDouble> > > m_zvel;

    };
    typedef std::shared_ptr<WomersleyParams> WomersleyParamsSharedPtr;

    struct BlowingSuctionParams
    {
        virtual ~BlowingSuctionParams()
        {};

        //DOF [x-trans ; y-trans ; z-rot]
        Array<OneD, int>                                 m_DOF;
        //Mass
        NekDouble                                        m_M;
        //Inertia
        NekDouble                                        m_I;
        // Spring coefficient
        Array<OneD, NekDouble>                           m_K;
        // Added Spring coefficient
        NekDouble                                        m_Kadd;
        // Damping coefficient
        Array<OneD, NekDouble>                           m_C;
        // Hinge point 
        Array<OneD, NekDouble>                           m_hingePoint;
        // Rotation axis 
        Array<OneD, NekDouble>                           m_axis;
        // AeroForces filter
        SolverUtils::FilterAeroForcesSharedPtr           m_filterForces;
        // AeroForces filter
        SolverUtils::FilterAeroForcesSharedPtr           m_filterForcesA;
        // AeroForces filter
        SolverUtils::FilterAeroForcesSharedPtr           m_filterForcesB;
        // AeroForces filter
        SolverUtils::FilterAeroForcesSharedPtr           m_filterForcesAS;
        /// Soft pointer to original equation system 
        const std::weak_ptr<SolverUtils::EquationSystem> m_equ;

        // Boundary regions
        std::string                             m_boundaryList;

        /// Time integration order
        int                                     m_intSteps;
        /// Initial time when the body is fixed, to prevent instability in startup
        NekDouble                               m_startTime;

        // Output information
        unsigned int                            m_outputFrequency;
        std::string                             m_outputFile;
        std::ofstream                           m_outputStream;
        bool                                    m_doOutput;
        unsigned int                            m_index;


        // Variables for time integration
        NekDouble                               m_angle;
        NekDouble                               m_previousAngle;
        Array<OneD, NekDouble>                  m_angleVel;
        Array<OneD, NekDouble>                  m_angleAdj;
        Array<OneD, NekDouble>                  m_angleAdjPrev;

        Array<OneD, NekDouble>                  m_dof;
        Array<OneD, Array<OneD, NekDouble> >    m_dofVel;
        Array<OneD, Array<OneD, NekDouble> >    m_dofAdj;
        Array<OneD, NekDouble>                  m_moment;
        Array<OneD, Array<OneD, NekDouble> >    m_force;
        Array<OneD, Array<OneD, NekDouble> >    m_forceA;
        Array<OneD, Array<OneD, NekDouble> >    m_forceB;
        Array<OneD, NekDouble>                  m_momentA;
        Array<OneD, NekDouble>                  m_momentB;

        /// Storage for constants of scaling
        Array<OneD, Array<OneD, NekDouble> >    m_deltaGradBnd;
        Array<OneD, Array<OneD, NekDouble> >    m_deltaGrad;
        Array<OneD, Array<OneD, NekDouble> >    m_deltaGamma;
        Array<OneD, Array<OneD, NekDouble> >    m_deltaGammaBnd;
        Array<OneD, Array<OneD, NekDouble> >    m_GradBaseBnd;


        Array<OneD, bool>                       m_isBlowingSuction;
        bool                                    m_isMomentA;
        bool                                    m_isMomentB;

        bool                                    m_isModified;
        bool                                    m_isPitch;
        bool                                    m_isSway;

        NekDouble                               m_alpha;

        // Coefficients for Adams time-integration
        // NekDouble AdamsBashforth_coeffs[3][3];
        // NekDouble AdamsMoulton_coeffs[3][3];

        NekDouble AdamsBashforth_coeffs[3][3] = {
        { 1.0       , 0.0       , 0.0     },
        { 3.0/2.0   ,-1.0/2.0   , 0.0     },
        { 23.0/12.0 ,-4.0/3.0   , 5.0/12.0}};
        NekDouble AdamsMoulton_coeffs[3][3] = {
        { 1.0       ,  0.0      , 0.0     },
        { 1.0/2.0   ,  1.0/2.0  , 0.0     },
        { 5.0/12.0  ,  2.0/3.0  ,-1.0/12.0}};
    };
    typedef std::shared_ptr<BlowingSuctionParams> BlowingSuctionParamsSharedPtr;

    /**
     * \brief This class is the base class for Navier Stokes problems
     *
     */
    class IncNavierStokes: public SolverUtils::AdvectionSystem,
                           public SolverUtils::FluidInterface
    {
    public:
        // Destructor
        virtual ~IncNavierStokes();

        virtual void v_InitObject();

        int GetNConvectiveFields(void)
        {
            return m_nConvectiveFields;
        }

        Array<OneD, int> &GetVelocity(void)
        {
            return  m_velocity;
        }

        void AddForcing(const SolverUtils::ForcingSharedPtr& pForce);

        virtual void GetPressure(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble>                     &pressure);

        virtual void GetDensity(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble>                     &density);

        virtual bool HasConstantDensity()
        {
            return true;
        }

        virtual void GetVelocity(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD, Array<OneD, NekDouble> >       &velocity);

        /// Return theta and theta dot
        virtual void v_GetStruct(NekDouble &angle, NekDouble &angleVel);

        /// Return value of deltagrad on boundary
        virtual void v_GetDeltaGrad(Array<OneD, Array<OneD, NekDouble> > &deltaGrad);

        /// Return value of base flow gradient
        virtual void v_GetGradBase(Array<OneD, Array<OneD, NekDouble> > &gradBase);

        /// Return value of deltagamma on boundary
        virtual void v_GetDeltaGamma(Array<OneD, Array<OneD, NekDouble> > &deltaGamma);

        /// Return theta and theta dot
        virtual void v_SetStruct(NekDouble &angle, NekDouble &angleVel);

        /// Return scaling factor
        virtual void v_SetScalingFactor(NekDouble &alpha);

        /// Set BSBC flag
        virtual bool v_CheckBSBC();

        /// Get number of structural DOF for BSBC
        virtual void v_nDOF(int &nDOF); 

        /// Set moments flags
        virtual bool v_SetMoment(bool &isMomentA, bool &isMomentB);

        /// Set added stiffness flags
        virtual bool v_SetAddedStiff(bool &isModified, bool &isPitch, bool &isSway);

        // bool to check if BSBC needed
        bool m_BlowingSuction = false;

        /// Blowing suction parameters if required
        BlowingSuctionParamsSharedPtr                          m_bsbcParams;

        /// Number of steps for Arnoldi solver
        int                                                    m_numSteps;

        /// Iteration for Arnoldi solver
        int                                                    m_iteration;

        /// Period for Arnoldi solver
        NekDouble                                              m_finTime;

    protected:

        // pointer to the extrapolation class for sub-stepping and HOPBS
        ExtrapolateSharedPtr m_extrapolation;

        /// modal energy file
        std::ofstream m_mdlFile;

        /// bool to identify if advection term smoothing is requested
        bool m_SmoothAdvection;

        /// Forcing terms
        std::vector<SolverUtils::ForcingSharedPtr>               m_forcing;

        /// Number of fields to be convected;
        int   m_nConvectiveFields;

        /// int which identifies which components of m_fields contains the
        /// velocity (u,v,w);
        Array<OneD, int> m_velocity;

        /// Pointer to field holding pressure field
        MultiRegions::ExpListSharedPtr m_pressure;
        /// Kinematic viscosity
        NekDouble   m_kinvis;
        /// dump energy to file at steps time
        int         m_energysteps;

        /// equation type;
        EquationType  m_equationType;

        /// Mapping from BCs to Elmt IDs
        Array<OneD, Array<OneD, int> > m_fieldsBCToElmtID;
        /// Mapping from BCs to Elmt Edge IDs
        Array<OneD, Array<OneD, int> > m_fieldsBCToTraceID;
        /// RHS Factor for Radiation Condition
        Array<OneD, Array<OneD, NekDouble> > m_fieldsRadiationFactor;

        /// Number of time integration steps AND Order of extrapolation for
        /// pressure boundary conditions.
        int m_intSteps;

        /// Modified fields to use in FSI stability
        Array<OneD, MultiRegions::ExpListSharedPtr> m_modifiedBase;

        /// Constructor.
        IncNavierStokes(const LibUtilities::SessionReaderSharedPtr& pSession,
                        const SpatialDomains::MeshGraphSharedPtr &pGraph);

        EquationType GetEquationType(void)
        {
            return m_equationType;
        }

        void EvaluateAdvectionTerms(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
									Array<OneD, Array<OneD, NekDouble> > &outarray);

        void WriteModalEnergy(void);

        /// time dependent boundary conditions updating
        void SetBoundaryConditions(NekDouble time);

        /// Set Radiation forcing term
        void SetRadiationBoundaryForcing(int fieldid);

        /// Set Normal Velocity Component to Zero
        void SetZeroNormalVelocity();

        /// Set Womersley Profile if specified
        void SetWomersleyBoundary(const int fldid, const int bndid);

        /// Solve structural equation and obtain theta and theta dot
        void SolveStructuralDirect(NekDouble time);

        /// Solve structural equation and obtain theta and theta dot
        void SolveStructuralAdjoint(NekDouble time);

        /// Force structural equation and obtain theta and theta dot
        void ForceStructuralAdjoint(NekDouble time);

        /// Scale BCs by theta and theta dot
        void ScaleBSBCDirect();

        /// Scale BCs fo adjoint
        void ScaleBSBCAdjoint();

        /// Set Up Womersley details
        void SetUpWomersley(const int fldid, const int bndid, std::string womstr);

        /// Set up Blowing suction boundary conditions
        void SetUpBlowingSuction(std::string BSBCStr);

        /// Initialise structural scaling constants
        void SetStructConsts();

        /// Womersley parameters if required
        std::map<int, std::map<int,WomersleyParamsSharedPtr> > m_womersleyParams;

        virtual MultiRegions::ExpListSharedPtr v_GetPressure()
        {
            return m_pressure;
        }

        virtual void v_TransCoeffToPhys(void)
        {
            ASSERTL0(false,"This method is not defined in this class");
        }

        virtual void v_TransPhysToCoeff(void)
        {
            ASSERTL0(false,"This method is not defined in this class");
        }

        virtual int v_GetForceDimension()=0;

        virtual Array<OneD, NekDouble> v_GetMaxStdVelocity();

        virtual bool v_PreIntegrate(int step);

    private:

    };

    typedef std::shared_ptr<IncNavierStokes> IncNavierStokesSharedPtr;

} //end of namespace

#endif //NEKTAR_SOLVERS_INCNAVIERSTOKES_H
