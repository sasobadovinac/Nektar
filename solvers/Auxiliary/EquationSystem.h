///////////////////////////////////////////////////////////////////////////////
//
// File EquationSystem.h
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
// Description: Base class for individual solvers.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEM_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEM_H

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <SpatialDomains/HistoryPoints.h>
#include <SpatialDomains/SpatialData.h>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>

namespace Nektar
{
    class EquationSystem;

    /// A shared pointer to an EquationSystem object
    typedef boost::shared_ptr<EquationSystem> EquationSystemSharedPtr;
    /// Datatype of the NekFactory used to instantiate classes derived from
    /// the EquationSystem class.
    typedef LibUtilities::NekFactory<
                std::string, EquationSystem,
                LibUtilities::CommSharedPtr&,
                LibUtilities::SessionReaderSharedPtr&
            > EquationSystemFactory;
    EquationSystemFactory& GetEquationSystemFactory();

    /// A base class for describing how to solve specific equations.
    class EquationSystem
    {
    public:
        EquationSystem(void);

        /// Destructor
        virtual ~EquationSystem();

        void InitObject();

        virtual void v_InitObject();

        /// Perform any initialisation necessary before solving the problem.
        void DoInitialise(void);

        /// Solve the problem.
        void DoSolve(void);
		
        /// Print a summary of parameters and solver characteristics.
        void PrintSummary(std::ostream &out);

        /// Perform initialisation of the base flow.
        void InitialiseBaseFlow(Array<OneD, Array<OneD, NekDouble> > &base);

        /// Initialise the data in the dependent fields.
        inline void SetInitialConditions(NekDouble initialtime = 0.0,
                                  bool dumpInitialConditions = true);

        /// Evaluates an exact solution
        inline void EvaluateExactSolution(int field,
                Array<OneD, NekDouble> &outfield,
                const NekDouble time);

        /// Decide the kind of forcing functions
        void SetInitialForce(NekDouble initialtime=0.0) ;

        ///Initialise the dimendion of forcing functions
        void InitialiseForcingFunctions(
                         Array<OneD, MultiRegions::ExpListSharedPtr> &force);

        /// Perform output operations after solve.
        void Output();

        /// Compute the L2 error between fields and a given exact solution.
        NekDouble L2Error(int field,
                          const Array<OneD,NekDouble> &exactsoln,
                          bool Normalised = false);

        /// Compute the L2 error of the fields
        NekDouble L2Error(int field, bool Normalised = false)
        {
            return L2Error(field,NullNekDouble1DArray,Normalised);
        }

        /// Compute the L_inf error between fields and a given exact solution.
        NekDouble LinfError(int field,
                const Array<OneD,NekDouble> &exactsoln = NullNekDouble1DArray);

        ///Compute error (L2 and L_inf) over an larger set of quadrature points return [L2 Linf]
        Array<OneD,NekDouble> ErrorExtraPoints(int field);

         /// Compute the inner product \f$ (\nabla \phi \cdot F) \f$.
         void WeakAdvectionGreensDivergenceForm(
                 const Array<OneD, Array<OneD, NekDouble> > &F,
                 Array<OneD, NekDouble> &outarray);

         /// Compute the inner product \f$ (\phi, \nabla \cdot F) \f$.
         void WeakAdvectionDivergenceForm(
                 const Array<OneD, Array<OneD, NekDouble> > &F,
                 Array<OneD, NekDouble> &outarray);

         /// Compute the inner product \f$ (\phi, V\cdot \nabla u) \f$.
         void WeakAdvectionNonConservativeForm(
                 const Array<OneD, Array<OneD, NekDouble> > &V,
                 const Array<OneD, const NekDouble> &u,
                 Array<OneD, NekDouble> &outarray);

         /// Compute the non-conservative advection \f$ (V \cdot \nabla u) \f$.
         void AdvectionNonConservativeForm(
                 const Array<OneD, Array<OneD, NekDouble> > &V,
                 const Array<OneD, const NekDouble> &u,
                 Array<OneD, NekDouble> &outarray,
                 Array<OneD, NekDouble> &wk = NullNekDouble1DArray);

         /// Calculate the weak discontinuous Galerkin advection.
         void WeakDGAdvection(
                 const Array<OneD, Array<OneD, NekDouble> >& InField,
                 Array<OneD, Array<OneD, NekDouble> >& OutField,
                 bool NumericalFluxIncludesNormal = true,
                 bool InFieldIsInPhysSpace = false,
                 int nvariables = 0);

         /// Calculate weak DG Diffusion in the LDG form.
         void WeakDGDiffusion(
                 const Array<OneD, Array<OneD, NekDouble> >& InField,
                 Array<OneD, Array<OneD, NekDouble> >& OutField,
                 bool NumericalFluxIncludesNormal = true,
                 bool InFieldIsInPhysSpace = false);


         /// Write checkpoint file.
         void Checkpoint_Output(const int n);
         /// Write checkpoint file.
         void Checkpoint_Output(const int n, MultiRegions::ExpListSharedPtr &field, Array< OneD, Array<OneD, NekDouble> > &fieldcoeffs, Array<OneD, std::string> &variables);

         /// Write field data to the given filename.
         void WriteFld(std::string &outname);

         /// Write input fields to the given filename.
         void WriteFld(std::string &outname, MultiRegions::ExpListSharedPtr &field, Array<OneD, Array<OneD, NekDouble> > &fieldcoeffs, Array<OneD, std::string> &variables);

         /// Input field data from the given file.
         void ImportFld(std::string &infile, Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);

         /// Output a field.
         void Array_Output(const int n, std::string name,
                           const Array<OneD, const NekDouble>&inarray,
                           bool IsInPhysicalSpace);

         void WriteTecplotFile(const int n, std::string name, bool IsInPhysicalSpace);

         /// Builds map of which element holds each history point.
         void ScanForHistoryPoints();

         /// Probe each history point and write to file.
         void WriteHistoryData (std::ostream &out);

         /// Write out a full summary.
         void Summary          (std::ostream &out);

         /// Write out a session summary.
         void SessionSummary   (std::ostream &out);

         /// Write out a summary of the time parameters.
         void TimeParamSummary (std::ostream &out);

         inline Array<OneD, MultiRegions::ExpListSharedPtr> &UpdateFields(void)
         {
             return m_fields;
         }

         inline Array<OneD, MultiRegions::ExpListSharedPtr> &UpdateForces(void)
         {
             return m_forces;
         }

         /// Return final time
         inline NekDouble GetFinalTime()
         {
             return m_time;
         }

         inline int GetNcoeffs(void)
         {
             return m_fields[0]->GetNcoeffs();
         }

         inline int GetNcoeffs(const int eid)
         {
             return m_fields[0]->GetNcoeffs(eid);
         }

         inline int GetNumExpModes(void)
         {
             return m_graph->GetExpansions().begin()->second->m_basisKeyVector[0]
                                                         .GetNumModes();
         }

         inline const Array<OneD,int> GetNumExpModesPerExp(void)
         {
             return m_fields[0]->EvalBasisNumModesMaxPerExp();
         }

         inline int GetNvariables(void)
         {
             return m_fields.num_elements();
         }

         inline const std::string &GetVariable(unsigned int i)
         {
             return m_boundaryConditions->GetVariable(i);
         }

         inline int GetTraceTotPoints(void)
         {
             return GetTraceNpoints();
         }

         inline int GetTraceNpoints(void)
         {
             switch(m_expdim)
             {
                 case 1:
                     // can't have two &GetTrace in ExpList.h hmm...
                     //return m_fields[0]->GetTrace().num_elements();
                     break;
                 case 2:
                 case 3:
                     return m_fields[0]->GetTrace()->GetNpoints();
                     break;
                 default:
                     ASSERTL0(false,"illegal expansion dimension");
             }
         }

         inline int GetExpSize(void)
         {
           return m_fields[0]->GetExpSize();
         }

         inline int GetPhys_Offset(int n)
         {
           return m_fields[0]->GetPhys_Offset(n);
         }

         inline int GetCoeff_Offset(int n)
         {
             return m_fields[0]->GetCoeff_Offset(n);
         }

         inline int GetTotPoints(void)
         {
             return m_fields[0]->GetNpoints();
         }

         inline int GetTotPoints(int n)
         {
             return m_fields[0]->GetTotPoints(n);
         }

         inline int GetNpoints(void)
         {
             return m_fields[0]->GetNpoints();
         }

         inline int GetSteps(void)
         {
             return m_steps;
         }

         void ZeroPhysFields(void);

         void FwdTransFields(void);

         /// Type of Galerkin projection.
         enum ProjectionType
         {
             eGalerkin,
             eDiscontinuousGalerkin
         };

         //-----------------------------------------------------------
         // virtual functions wrappers

         void GetFluxVector(const int i,
                             Array<OneD, Array<OneD, NekDouble> >&physfield,
                             Array<OneD, Array<OneD, NekDouble> >&flux)
         {
             v_GetFluxVector(i,physfield, flux);
         }

         void GetFluxVector(const int i,
                             Array<OneD, Array<OneD, NekDouble> >&physfield,
                             Array<OneD, Array<OneD, NekDouble> >&fluxX,
                             Array<OneD, Array<OneD, NekDouble> > &fluxY)
         {
             v_GetFluxVector(i,physfield, fluxX, fluxY);
         }

         virtual void GetFluxVector(const int i, const int j,
                             Array<OneD, Array<OneD, NekDouble> > &physfield,
                             Array<OneD, Array<OneD, NekDouble> > &flux)
         {
             v_GetFluxVector(i,j,physfield,flux);
         }

         void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                            Array<OneD, Array<OneD, NekDouble> > &numflux)
         {
             v_NumericalFlux(physfield, numflux);
         }

         void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                            Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                            Array<OneD, Array<OneD, NekDouble> > &numfluxY)
         {
             v_NumericalFlux(physfield, numfluxX, numfluxY);
         }

         void NumFluxforScalar(Array<OneD, Array<OneD, NekDouble> > &ufield,
                     Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
         {
             v_NumFluxforScalar(ufield, uflux);
         }

         void NumFluxforVector(Array<OneD, Array<OneD, NekDouble> > &ufield,
                   Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                   Array<OneD, Array<OneD, NekDouble> >  &qflux)
         {
             v_NumFluxforVector(ufield,qfield, qflux);
         }

         NekDouble AdvectionSphere(const NekDouble x0j, const NekDouble x1j,
                                   const NekDouble x2j, const NekDouble time)
         {
             return v_AdvectionSphere(x0j, x1j, x2j, time);
         }

         NekDouble Morphogenesis(const int field, const NekDouble x0j,
                                 const NekDouble x1j, const NekDouble x2j,
                                 const NekDouble time)
         {
             return v_Morphogenesis(field, x0j, x1j, x2j, time);
         }

         /// Number of Quadrature points used to work out the error
         int  m_NumQuadPointsError;
         bool m_UseContCoeff;

         ///Parameter for homogeneous expansions
         enum HomogeneousType
         {
             eHomogeneous1D,
             eHomogeneous2D,
             eHomogeneous3D,
             eNotHomogeneous
         };

         bool m_useFFT;               ///< flag to determine if use or not the FFT for transformations

         enum HomogeneousType m_HomogeneousType;

         NekDouble m_LhomX; ///< physical length in X direction (if homogeneous)
         NekDouble m_LhomY; ///< physical length in Y direction (if homogeneous)
         NekDouble m_LhomZ; ///< physical length in Z direction (if homogeneous)

         int m_npointsX;    ///< number of points in X direction (if homogeneous)
         int m_npointsY;    ///< number of points in Y direction (if homogeneous)
         int m_npointsZ;    ///< number of points in Z direction (if homogeneous)

         int m_HomoDirec;   ///< number of homogenous directions

    protected:
         /// Communicator
         LibUtilities::CommSharedPtr                 m_comm;
         /// The session reader
         LibUtilities::SessionReaderSharedPtr        m_session;
         /// Array holding all dependent variables.
         Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
         /// Base fields.
         Array<OneD, MultiRegions::ExpListSharedPtr> m_base;
         /// Array holding force values.
         Array<OneD, MultiRegions::ExpListSharedPtr> m_forces;
         /// variable that determine if the force is necessary or not.
         bool                                        m_bforce;
         /// Array holding all dependent variables.
         Array<OneD, MultiRegions::ExpListSharedPtr> m_derivedfields;
         /// dimension force array
         int                                         m_FDim;
         /// Pointer to boundary conditions object.
         SpatialDomains::BoundaryConditionsSharedPtr m_boundaryConditions;
         /// Pointer to history data object.
         SpatialDomains::HistorySharedPtr            m_historyPoints;
         /// Pointer to graph defining mesh.
         SpatialDomains::MeshGraphSharedPtr          m_graph;

         std::list<std::pair<SpatialDomains::VertexComponentSharedPtr, int> >
                                                     m_historyList;

         SpatialDomains::SpatialParametersSharedPtr  m_spatialParameters;

         std::string m_filename;      ///< Filename
         std::string m_sessionName;   ///< Name of the sessions
         NekDouble m_time;            ///< Continous time
         NekDouble m_fintime;         ///< time to be taken during the simulation
         NekDouble m_timestep;        ///< Time step size
         int m_steps;                 ///< Number of steps to take
         int m_checksteps;            ///< Number of steps between checkpoints
         int m_spacedim;              ///< Spatial dimension (> expansion dim)
         int m_expdim;                ///< Dimension of the expansion

         MultiRegions::GlobalSysSolnType m_solnType;

         /// Type of projection, i.e. Galerkin or DG.
         enum ProjectionType m_projectionType;

         /// Array holding the forward normals.
         Array<OneD, Array<OneD, NekDouble> > m_traceNormals;
         /// 1 x nvariable x nq
         Array<OneD, Array<OneD, Array<OneD,NekDouble> > > m_gradtan;
         /// 2 x m_spacedim x nq
         Array<OneD, Array<OneD, Array<OneD,NekDouble> > > m_tanbasis;

         /// Flag to indicate if the fields should be checked for singularity.
         Array<OneD, bool>                               m_checkIfSystemSingular;


         /// Initialises EquationSystem class members.
         EquationSystem( LibUtilities::CommSharedPtr& pComm,
                         LibUtilities::SessionReaderSharedPtr& pSession);

         /// Perform a case-insensitive string comparison.
         int NoCaseStringCompare(const string & s1, const string& s2) ;

         // Here for consistency purposes with old version
         int nocase_cmp(const string & s1, const string& s2)
         {
             return NoCaseStringCompare(s1,s2);
         }

         /// Check for and load an integer parameter
         void LoadParameter(std::string name, int &var, int def = 0);

         /// Check for and load a double precision parameter
         void LoadParameter(std::string name, NekDouble &var, NekDouble def= 0.0);

        /// Evaluates a function as specified in the session file.
        void EvaluateFunction(Array<OneD, NekDouble>& pArray,
                SpatialDomains::ConstUserDefinedEqnShPtr pEqn,
                const NekDouble pTime = 0.0);

        /// Evaluates a function as specified in the session file.
        void EvaluateFunction(
                Array<OneD, Array<OneD, NekDouble> >& pArray,
                std::string pFunctionName,
                const NekDouble pTime = 0.0);

        /// Populate given fields with the forcing function from session.
        void EvaluateFunction(
                Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                const std::string& pName);

        /// Evaluates the boundary conditions at the given time.
        void SetBoundaryConditions(NekDouble time);

        /// Virtual function for initialisation implementation.
        virtual void v_DoInitialise();

        /// Virtual function for solve implementation.
        virtual void v_DoSolve();

        /// Virtual function for printing summary information.
        virtual void v_PrintSummary(std::ostream &out);

        virtual void v_SetInitialConditions(NekDouble initialtime = 0.0,
                            bool dumpInitialConditions = true);

        virtual void v_EvaluateExactSolution(int field,
                        Array<OneD, NekDouble> &outfield,
                        const NekDouble time);

        //Initialise m_base in order to store the base flow from a file 
        void SetUpBaseFields( SpatialDomains::MeshGraphSharedPtr &mesh);
        
        // Fill m_base with the values stored in a fld file
        void ImportFldBase(std::string pInfile, SpatialDomains::MeshGraphSharedPtr 
    	    pGraph);

        // Ouptut field information
        virtual void v_Output(void);

    private:

        virtual Array<OneD, bool> v_GetSystemSingularChecks()
        {
            return Array<OneD, bool>(m_boundaryConditions->GetNumVariables(), false);
        }

        virtual void v_GetFluxVector(const int i, Array<OneD,
                            Array<OneD, NekDouble> >&physfield,
                            Array<OneD, Array<OneD, NekDouble> >&flux)
        {
            ASSERTL0(false, "v_GetFluxVector: This function is not valid "
                            "for the Base class");
        }

        virtual void v_GetFluxVector(const int i, const int j,
                            Array<OneD, Array<OneD, NekDouble> >&physfield,
                            Array<OneD, Array<OneD, NekDouble> >&flux)
        {
            ASSERTL0(false, "v_GetqFluxVector: This function is not valid "
                            "for the Base class");
        }

        virtual void v_GetFluxVector(const int i, Array<OneD,
                            Array<OneD, NekDouble> >&physfield,
                            Array<OneD, Array<OneD, NekDouble> >&fluxX,
                            Array<OneD, Array<OneD, NekDouble> > &fluxY)
        {
            ASSERTL0(false, "v_GetFluxVector: This function is not valid "
                            "for the Base class");
        }

        virtual void v_NumericalFlux(
                            Array<OneD, Array<OneD, NekDouble> > &physfield,
                            Array<OneD, Array<OneD, NekDouble> > &numflux)
        {
            ASSERTL0(false, "v_NumericalFlux: This function is not valid "
                            "for the Base class");
        }

        virtual void v_NumericalFlux(
                            Array<OneD, Array<OneD, NekDouble> > &physfield,
                            Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                            Array<OneD, Array<OneD, NekDouble> > &numfluxY )
        {
            ASSERTL0(false, "v_NumericalFlux: This function is not valid "
                            "for the Base class");
        }

        virtual void v_NumFluxforScalar(
                    Array<OneD, Array<OneD, NekDouble> > &ufield,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
        {
            ASSERTL0(false, "v_NumFluxforScalar: This function is not valid "
                            "for the Base class");
        }

        virtual void v_NumFluxforVector(
                    Array<OneD, Array<OneD, NekDouble> > &ufield,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                    Array<OneD, Array<OneD, NekDouble > >  &qflux)
        {
            ASSERTL0(false, "v_NumFluxforVector: This function is not valid "
                            "for the Base class");
        }

        virtual NekDouble v_AdvectionSphere(const NekDouble x0j,
                            const NekDouble x1j, const NekDouble x2j,
                            const NekDouble time)
        {
            ASSERTL0(false, "v_AdvectionSphere: This function is not valid "
                            "for the Base class");
            return 0.0;
        }

        virtual NekDouble v_Morphogenesis(const int field, const NekDouble x0j,
                            const NekDouble x1j, const NekDouble x2j,
                            const NekDouble time)
        {
            ASSERTL0(false, "v_Morphogenesis: This function is not valid "
                            "for the Base class");
            return 0.0;
        }

    };

    inline void EquationSystem::SetInitialConditions(NekDouble initialtime,
                              bool dumpInitialConditions)
    {
        v_SetInitialConditions(initialtime,dumpInitialConditions);
    }

    /// Evaluates an exact solution
    inline void EquationSystem::EvaluateExactSolution(int field,
            Array<OneD, NekDouble> &outfield,
            const NekDouble time)
    {
        v_EvaluateExactSolution(field, outfield, time);
    }


}

#endif
