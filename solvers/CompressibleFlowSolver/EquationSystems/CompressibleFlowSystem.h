///////////////////////////////////////////////////////////////////////////////
//
// File CompressibleFlowSystem.h
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
// Description: Auxiliary functions for the compressible flow system
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_COMPRESSIBLEFLOWSYSTEM_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_COMPRESSIBLEFLOWSYSTEM_H

#include <CompressibleFlowSolver/ArtificialDiffusion/ArtificialDiffusion.h>
#include <CompressibleFlowSolver/Misc/VariableConverter.h>
#include <CompressibleFlowSolver/BoundaryConditions/CFSBndCond.h>
#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Diffusion/Diffusion.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/Expansion2D.h>

#define DEMO_IMPLICITSOLVER_JFNK_COEFF
#define CFS_DEBUGMODE
namespace Nektar
{
    /**
     *
     */
    class CompressibleFlowSystem: public SolverUtils::AdvectionSystem
    {
    public:

        friend class MemoryManager<CompressibleFlowSystem>;

        virtual ~CompressibleFlowSystem();
 
        /// Function to calculate the stability limit for DG/CG.
        NekDouble GetStabilityLimit(int n);

        /// Function to calculate the stability limit for DG/CG
        /// (a vector of them).
        Array<OneD, NekDouble> GetStabilityLimitVector(
            const Array<OneD,int> &ExpOrder);

        /// Function to get estimate of min h/p factor per element
        Array<OneD, NekDouble>  GetElmtMinHP(void);
    protected:
        SolverUtils::DiffusionSharedPtr     m_diffusion;
        ArtificialDiffusionSharedPtr        m_artificialDiffusion;
        Array<OneD, Array<OneD, NekDouble> >m_vecLocs;
        NekDouble                           m_gamma;
        std::string                         m_shockCaptureType;

        // Parameters for exponential filtering
        NekDouble                           m_filterAlpha;
        NekDouble                           m_filterExponent;
        NekDouble                           m_filterCutoff;

        NekDouble                           m_JFEps;

        NekDouble                        m_ForcingGama  = 1.0;
        NekDouble                        m_ForcingAlpha = 0.5 * (1.0 + sqrt(5));

        TensorOfArray3D<NekDouble>                     m_qfield;
        Array<OneD, Array<OneD, NekDouble> >           m_MatrixFreeRefFields;
        Array<OneD, Array<OneD, NekDouble> >           m_MatrixFreeRefFwd;
        Array<OneD, Array<OneD, NekDouble> >           m_MatrixFreeRefBwd;
        TensorOfArray2D<DNekBlkMatSharedPtr>           m_ElmtFluxJacArray;
        
        bool                                m_useFiltering;

        /// Store physical artificial viscosity
        Array<OneD, NekDouble>              m_muav;

        /// Store physical artificial viscosity
        Array<OneD, NekDouble>              m_muavTrace;

        // Parameters for local time-stepping
        bool                                m_useLocalTimeStep;

        bool                                m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT;
        bool                                m_DEBUG_VISCOUS_JAC_MAT;
        bool                                m_DEBUG_ADVECTION_JAC_MAT;
        bool                                m_centralDiffTracJac;

        bool                                m_CalcTracePartFlag;
        bool                                m_CalcVolumPartFlag;

        bool                                m_flagPrecMatFree = false;

#ifdef CFS_DEBUGMODE
       // 1: Adv; 2: Dif; Default: all
        int                                 m_DebugAdvDiffSwitch; 
       // 1: Vol; 2: Trace; Default: all
        int                                 m_DebugVolTraceSwitch; 
       // 1: Con; 2: Deriv; Default: all
        int                                 m_DebugConsDerivSwitch; 


        int                                 m_DebugNumJacMatSwitch;
        int                                 m_DebugOutputJacMatSwitch;

        int                                 m_DebugInvMassSwitch   ; 
        int                                 m_DebugPlusSourceSwitch; 

        int                                 m_DebugIPSymmFluxJacSwitch; 
        int                                 m_DebugNumJacBSOR;
#endif

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
        TensorOfArray4D<NekDouble>          m_StdDMatDataDBB;
        TensorOfArray5D<NekDouble>          m_StdDMatDataDBDB;

        TensorOfArray4D<NekSingle>          m_StdSMatDataDBB;
        TensorOfArray5D<NekSingle>          m_StdSMatDataDBDB;
        int                                 m_nPadding = 1;
#endif
        int                                 m_LiniearizationMethod;

        // Auxiliary object to convert variables
        VariableConverterSharedPtr          m_varConv;

        // User defined boundary conditions
        std::vector<CFSBndCondSharedPtr>    m_bndConds;

        // Forcing term
        std::vector<SolverUtils::ForcingSharedPtr> m_forcing;

        NekLinSysIterativeSharedPtr                 m_linsolBRJ;

        enum PreconditionerType
        {
            eNull,    ///< No Solution type specified
            eDiagonal,
            eSparse,
        };
        PreconditionerType                  m_PrecMatStorage;

        enum PrecDiagOffDiagTag
        {
            ePrecDiagOffDiagTagAll,    
            ePrecDiagOffDiagTagDiag,
            ePrecDiagOffDiagTagOffdiag
        };

        CompressibleFlowSystem(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr& pGraph);
        
        void v_DoOdeProjection1(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
              const NekDouble                              time);
        
        void v_DoOdeRhs1(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
              const NekDouble                              time);

        virtual void v_InitObject();

        void InitialiseParameters();

        void InitAdvection();

        void DoOdeRhs(
            const TensorOfArray2D<NekDouble>        &inarray,
            Array<OneD, Array<OneD, NekDouble> >    &outarray,
            const NekDouble                         time);
        void DoOdeProjection(
            const TensorOfArray2D<NekDouble>        &inarray,
            Array<OneD, Array<OneD, NekDouble> >    &outarray,
            const NekDouble                         time);

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
        void preconditioner(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&out);
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void preconditioner_BlkDiag(
            const Array<OneD, NekDouble>                  &inarray,
            Array<OneD, NekDouble >                       &outarray,
            const TensorOfArray2D<TypeNekBlkMatSharedPtr> &PrecMatVars,
            const DataType                                &tmpDataType);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void preconditioner_BlkDiag(
            const Array<OneD, NekDouble>  &inarray,
            Array<OneD, NekDouble >       &outarray,
            const TypeNekBlkMatSharedPtr  &PrecMatVars,
            const DataType                &tmpDataType);

        template<typename DataType>
        void preconditioner_BlkDiagMatFree(
            const TensorOfArray1D<NekDouble>     &inarray,
            TensorOfArray1D<NekDouble>           &outarray,
            const DataType                       &tmpDataType);

        void PrecBJRDiagMult(
            const  Array<OneD, NekDouble>   &inarray,
            Array<OneD, NekDouble >         &out,
            const  bool                     &controlFlag);

        void PrecBJRDiagPrec(
            const Array<OneD, NekDouble>    &inarray,
            Array<OneD, NekDouble >         &outarray,
            const bool                      &flag);

        void preconditioner_NumJac(
            const Array<OneD, NekDouble>                  &inarray,
            Array<OneD, NekDouble >                       &outarray,
            const TensorOfArray2D<DNekBlkMatSharedPtr>    &PrecMatVars,
            const Array<OneD, Array<OneD, NekDouble > >   &PrecMatVarsOffDiag);
        void MinusOffDiag2RhsNumJac(
            const int                                     nvariables,
            const int                                     nCoeffs,
            const Array<OneD, NekDouble>                  &inarray,
            Array<OneD, NekDouble>                        &outarray,
            const Array<OneD, Array<OneD, NekDouble > >   &PrecMatVarsOffDiag);
            
        void preconditioner_BlkSOR_coeff(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&outarray,
            const bool                   &Continueflag);

        void preconditioner_MultiLevel_coeff(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&outarray,
                  const bool             &UnusedFlag);

        // void MinusOffDiag2Rhs(
        //     const int nvariables,
        //     const int nCoeffs,
        //     const TensorOfArray2D<NekDouble>    &inarray,
        //           Array<OneD,       Array<OneD, NekDouble> >    &outarray,
        //     bool                                                flagUpdateDervFlux,
        //           Array<OneD,       Array<OneD, NekDouble> >    &FwdFluxDeriv,
        //           Array<OneD,       Array<OneD, NekDouble> >    &BwdFluxDeriv,
        //     TensorOfArray3D<NekDouble>  &qfield,
        //     TensorOfArray3D<NekDouble>  &tmpTrace);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void MinusOffDiag2Rhs(
            const int                                   nvariables,
            const int                                   nCoeffs,
            const TensorOfArray2D<NekDouble>            &inarray,
            Array<OneD, Array<OneD, NekDouble> >        &outarray,
            bool                                        flagUpdateDervFlux,
            Array<OneD, Array<OneD, NekDouble> >        &FwdFluxDeriv,
            Array<OneD, Array<OneD, NekDouble> >        &BwdFluxDeriv,
            TensorOfArray3D<NekDouble>                  &qfield,
            TensorOfArray3D<NekDouble>                  &tmpTrace,
            const Array<OneD, TypeNekBlkMatSharedPtr >  &TraceJac,
            const Array<OneD, TypeNekBlkMatSharedPtr >  &TraceJacDeriv,
            const Array<OneD, Array<OneD, DataType> >   &TraceJacDerivSign);

        template<typename DataType>
        void AddOrMinusPrecTraceFlux(
            const int                            &nswitchDiagOffdiag,
            const NekDouble                      AddOrMinusSign,
            const TensorOfArray2D<NekDouble>     &inarray,
            TensorOfArray2D<NekDouble>           &outarray,
            bool                                 flagUpdateDervFlux,
            TensorOfArray2D<NekDouble>           &FwdFluxDeriv,
            TensorOfArray2D<NekDouble>           &BwdFluxDeriv,
            TensorOfArray3D<NekDouble>           &qfield,
            TensorOfArray3D<NekDouble>           &wspTrace,
            TensorOfArray2D<DataType>            &wspTraceDataType,
            const TensorOfArray4D<DataType>      &TraceJacArray,
            const TensorOfArray4D<DataType>      &TraceJacDerivArray,
            const TensorOfArray2D<DataType>      &TraceJacDerivSign,
            const TensorOfArray5D<DataType>      &TraceIPSymJacArray);
        
        template<typename DataType>
        void CalTraceFwdBwdFlux(
            const TensorOfArray2D<NekDouble>         &Fwd,
            const TensorOfArray2D<NekDouble>         &Bwd,
            TensorOfArray2D<NekDouble>               &FwdFlux,
            TensorOfArray2D<NekDouble>               &BwdFlux,
            TensorOfArray2D<DataType>                &wspTraceDataType,
            const TensorOfArray4D<DataType>          &TraceJacArray);

        template<typename DataType>
        void AddTraceFwdBwdDerivFlux(
            TensorOfArray2D<NekDouble>               &FwdFlux,
            TensorOfArray2D<NekDouble>               &BwdFlux,
            TensorOfArray3D<NekDouble>               &numDerivFwd,
            TensorOfArray3D<NekDouble>               &numDerivBwd,
            TensorOfArray3D<NekDouble>               &qfield,
            TensorOfArray2D<DataType>                &wspTraceDataType,
            const TensorOfArray4D<DataType>          &TraceJacDerivArray,
            const TensorOfArray2D<DataType>          &TraceJacDerivSign);
        
        template<typename DataType>
        void CalTraceFwdBwdSymFlux(
            const TensorOfArray2D<NekDouble>         &Fwd,
            const TensorOfArray2D<NekDouble>         &Bwd,
            TensorOfArray3D<NekDouble>               &numDerivFwd,
            TensorOfArray3D<NekDouble>               &numDerivBwd,
            TensorOfArray2D<DataType>                &wspTraceDataType,
            const TensorOfArray2D<DataType>          &TraceJacDerivSign,
            const TensorOfArray5D<DataType>          &TraceIPSymJacArray);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void AddMatNSBlkDiag_volume(
            const TensorOfArray2D<NekDouble>         &inarray,
            const TensorOfArray3D<NekDouble>         &qfield,
            TensorOfArray2D<TypeNekBlkMatSharedPtr>  &gmtxarray,
            TensorOfArray4D<DataType>                &StdMatDataDBB,
            TensorOfArray5D<DataType>                &StdMatDataDBDB);

        template<typename DataType>
        void CalcVolJacStdMat(
            TensorOfArray4D<DataType>     &StdMatDataDBB,
            TensorOfArray5D<DataType>     &StdMatDataDBDB);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void AddMatNSBlkDiag_boundary(
            const TensorOfArray2D<NekDouble>         &inarray,
            TensorOfArray3D<NekDouble>               &qfield,
            TensorOfArray2D<TypeNekBlkMatSharedPtr>  &gmtxarray,
            Array<OneD, TypeNekBlkMatSharedPtr >     &TraceJac,
            Array<OneD, TypeNekBlkMatSharedPtr >     &TraceJacDeriv,
            Array<OneD, Array<OneD, DataType> >      &TraceJacDerivSign,
            TensorOfArray5D<DataType>                &TraceIPSymJacArray);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void ElmtVarInvMtrx(
            TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
            TypeNekBlkMatSharedPtr                  &gmtVar,
            const DataType                          &tmpDatatype);
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void ElmtVarInvMtrx(
            TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
            const DataType                          &tmpDatatype);
        
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void GetTraceJac(
            const TensorOfArray2D<NekDouble>      &inarray,
            const TensorOfArray3D<NekDouble>      &qfield,
            Array<OneD, TypeNekBlkMatSharedPtr >  &TraceJac,
            Array<OneD, TypeNekBlkMatSharedPtr >  &TraceJacDeriv,
            Array<OneD, Array<OneD, DataType> >   &TraceJacDerivSign,
            TensorOfArray5D<DataType>             &TraceIPSymJacArray);
       
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void NumCalRiemFluxJac(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &AdvVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            const TensorOfArray3D<NekDouble>                  &qfield,
            const NekDouble                                   &time,
            const Array<OneD, Array<OneD, NekDouble> >        &Fwd,
            const Array<OneD, Array<OneD, NekDouble> >        &Bwd,
            TypeNekBlkMatSharedPtr                            &FJac,
            TypeNekBlkMatSharedPtr                            &BJac,
            TensorOfArray5D<DataType>                      &TraceIPSymJacArray);

        void PointFluxJacobian_pn(
            const Array<OneD, NekDouble> &Fwd,
            const Array<OneD, NekDouble> &normals,
                  DNekMatSharedPtr       &FJac,
            const NekDouble efix,   const NekDouble fsw);

        void PointFluxJacobian_pn2D(
            const Array<OneD, NekDouble> &Fwd,
            const Array<OneD, NekDouble> &normals,
                  DNekMatSharedPtr       &FJac,
            const NekDouble efix,   const NekDouble fsw);

#ifdef CFS_DEBUGMODE

        void NumJacElemental(
            DNekMatSharedPtr    &NumericalJacobianMatrix,
            const int                 RowElementID,
            const int                 ColElementID);
        
        void CalOffDiagJacByMinusOffDiagElemental(
            DNekMatSharedPtr    &MinusoffJacobianMatrix,
            const int                 RowElementID,
            const int                 ColElementID);
        void DebugCheckJac(
            const int                 RowElementID,
            const int                 ColElementID);
        
        void DebugNumCalJac_coeff(
            TensorOfArray2D<DNekBlkMatSharedPtr>   &gmtxarray,
            Array<OneD, Array<OneD, NekDouble > >  &JacOffDiagArray =
                NullNekDoubleArrayofArray);
            
        void DebugNumCalElmtJac_coeff(
            Array<OneD, Array<OneD, DNekMatSharedPtr> > &ElmtPrecMatVars ,
            const int                                   nelmt,
            Array<OneD, Array<OneD, NekDouble > >       &JacOffDiagArray);
        void DebugNumCalElmtJac_coeff(
            Array<OneD, Array<OneD, DNekMatSharedPtr> > &ElmtPrecMatVars ,
            const int                                   nelmt);

        void NonlinSysEvaluator_coeff_out(
                Array<OneD, Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, NekDouble> > &out);

        template<typename TypeNekBlkMatSharedPtr>
        void CoutBlkMat(
            TypeNekBlkMatSharedPtr &gmtx, 
            const unsigned int nwidthcolm=12);

        template<typename TypeNekBlkMatSharedPtr>
        void CoutStandardMat(
            TypeNekBlkMatSharedPtr &loc_matNvar,
            const unsigned int nwidthcolm=12);

        template<typename TypeNekBlkMatSharedPtr>
        void Cout1DArrayBlkMat(
            Array<OneD, TypeNekBlkMatSharedPtr> &gmtxarray,
            const unsigned int nwidthcolm=12);
        
        template<typename TypeNekBlkMatSharedPtr>
        void Cout2DArrayBlkMat(
            TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
            const unsigned int                      nwidthcolm = 12);

        template<typename TypeNekBlkMatSharedPtr>
        void Cout2DArrayStdMat(
            TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
            const unsigned int                      nwidthcolm=12);
#endif
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void TranSamesizeBlkDiagMatIntoArray(
            const TypeNekBlkMatSharedPtr    &BlkMat,
            TensorOfArray3D<DataType>       &MatArray);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void TransTraceJacMatToArray(
            const Array<OneD, TypeNekBlkMatSharedPtr > &TraceJac,
            const Array<OneD, TypeNekBlkMatSharedPtr > &TraceJacDeriv,
            TensorOfArray4D<DataType>                  &TraceJacArray,
            TensorOfArray4D<DataType>                  &TraceJacDerivArray);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void Fill2DArrayOfBlkDiagonalMat(
            TensorOfArray2D<TypeNekBlkMatSharedPtr>   &gmtxarray,
            const DataType                            valu);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void Fill1DArrayOfBlkDiagonalMat( 
            Array<OneD, TypeNekBlkMatSharedPtr >    &gmtxarray,
            const DataType                          valu);

        void DoImplicitSolve_phy2coeff(
            const TensorOfArray2D<NekDouble>        &inarray,
            Array<OneD, Array<OneD, NekDouble> >    &out,
            const NekDouble                         time,
            const NekDouble                         lambda);

        void DoImplicitSolve_coeff(
            const TensorOfArray2D<NekDouble>        &inpnts,
            const TensorOfArray2D<NekDouble>        &inarray,
            Array<OneD, Array<OneD, NekDouble> >    &out,
            const NekDouble                         time,
            const NekDouble                         lambda);
        void UpdateSoltnRefNorms(
            const Array<OneD, Array<OneD, NekDouble>>   &inarray,
            const NekDouble                             &ototalDOF);

        bool NewtonStopCriteria(
            const Array<OneD, NekDouble>    &NonlinSysRes_1D,
            const int                       &k,
            const int                       &totalDOF,
            const NekDouble                 &resnorm,
            NekDouble                       &resnorm0,
            NekDouble                       &resmaxm,
            NekDouble                       &resratio);

        void CalPrecMat(
            const TensorOfArray2D<NekDouble> &inpnts,
            const NekDouble                  time,
            const NekDouble                  lambda);

        template<typename TypeNekBlkMatSharedPtr>
        void AllocatePrecondBlkDiag_coeff(
            TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
            const int                               &nscale=1 );

        inline void AllocateNekBlkMatDig(
            SNekBlkMatSharedPtr               &mat,
            const Array<OneD, unsigned int >  nrow,
            const Array<OneD, unsigned int >  ncol)
        {
            mat = MemoryManager<SNekBlkMat>
                ::AllocateSharedPtr(nrow, ncol, eDIAGONAL);
            SNekMatSharedPtr loc_matNvar;
            for(int nelm = 0; nelm < nrow.num_elements(); ++nelm)
            {
                int nrowsVars = nrow[nelm];
                int ncolsVars = ncol[nelm];
                
                loc_matNvar = MemoryManager<SNekMat>::
                        AllocateSharedPtr(nrowsVars,ncolsVars,0.0);
                mat->SetBlock(nelm,nelm,loc_matNvar);
            }
        }

        inline void AllocateNekBlkMatDig(
            DNekBlkMatSharedPtr              &mat,
            const Array<OneD, unsigned int > nrow,
            const Array<OneD, unsigned int > ncol)
        {
            mat = MemoryManager<DNekBlkMat>
                ::AllocateSharedPtr(nrow, ncol, eDIAGONAL);
            DNekMatSharedPtr loc_matNvar;
            for(int nelm = 0; nelm < nrow.num_elements(); ++nelm)
            {
                int nrowsVars = nrow[nelm];
                int ncolsVars = ncol[nelm];
                
                loc_matNvar = MemoryManager<DNekMat>::
                        AllocateSharedPtr(nrowsVars,ncolsVars,0.0);
                mat->SetBlock(nelm,nelm,loc_matNvar);
            }
        }
        void CalcInexactNewtonForcing(
            const int       &k,
            const NekDouble &resnormOld,
            const NekDouble &resnorm,
            NekDouble       &forcing);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void GetpreconditionerNSBlkDiag_coeff(
            const TensorOfArray2D<NekDouble>        &inarray,
            TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
            TypeNekBlkMatSharedPtr                  &gmtVar,
            Array<OneD, TypeNekBlkMatSharedPtr >    &TraceJac,
            Array<OneD, TypeNekBlkMatSharedPtr >    &TraceJacDeriv,
            Array<OneD, Array<OneD, DataType> >     &TraceJacDerivSign,
            TensorOfArray4D<DataType>               &TraceJacArray,
            TensorOfArray4D<DataType>               &TraceJacDerivArray,
            TensorOfArray5D<DataType>               &TraceIPSymJacArray,
            TensorOfArray4D<DataType>               &StdMatDataDBB,
            TensorOfArray5D<DataType>               &StdMatDataDBDB);
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void GetpreconditionerNSBlkDiag_coeff(
            const TensorOfArray2D<NekDouble>          &inarray,
            TensorOfArray2D<TypeNekBlkMatSharedPtr>   &gmtxarray,
            Array<OneD, TypeNekBlkMatSharedPtr >      &TraceJac,
            Array<OneD, TypeNekBlkMatSharedPtr >      &TraceJacDeriv,
            Array<OneD, Array<OneD, DataType> >       &TraceJacDerivSign,
            TensorOfArray5D<DataType>                 &TraceIPSymJacArray,
            TensorOfArray4D<DataType>                 &StdMatDataDBB,
            TensorOfArray5D<DataType>                 &StdMatDataDBDB);

        void MatrixMultiply_JacobianFree_coeff(
            const  Array<OneD, NekDouble> &inarray,
                   Array<OneD, NekDouble >&out);
        void MatrixMultiply_MatrixFree_coeff_noSource(
            const  Array<OneD, NekDouble> &inarray,
                   Array<OneD, NekDouble >&out);
        void MatrixMultiply_JacobianFree_coeff_central(
            const  Array<OneD, NekDouble> &inarray,
                Array<OneD, NekDouble >&out);
        void MatrixMultiply_MatrixFree_coeff(
            const  Array<OneD, NekDouble> &inarray,
                Array<OneD, NekDouble >&out);
        void MatrixMultiply_MatrixFree_coeff_FourthCentral(
            const  Array<OneD, NekDouble> &inarray,
                Array<OneD, NekDouble >&out);

        void MatrixMultiply_JacobianFree_coeff_dualtimestep(
            const  Array<OneD, NekDouble> &inarray,
                Array<OneD, NekDouble >&out,
            const  bool                   &controlFlag);

        void NonlinSysEvaluator_coeff(
                Array<OneD, Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, NekDouble> > &out);
        
        void NonlinSysEvaluator_coeff_noSource(
                Array<OneD, Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, NekDouble> > &out);

        void DoOdeRhs_coeff(
            const TensorOfArray2D<NekDouble>     &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble                      time,
            const bool                           flagFreezeJac = false);
        
        void DoOdeRhs_coeffVol(
            const TensorOfArray2D<NekDouble>        &inarray,
            TensorOfArray2D<NekDouble>              &outarray,
            const NekDouble                         time,
            const bool                              flagFreezeJac=false);

        void DoAdvection_coeff(
            const TensorOfArray2D<NekDouble>            &inarray,
            Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const NekDouble                             time,
            const Array<OneD, Array<OneD, NekDouble> >  &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >  &pBwd,
            const bool                                  flagFreezeJac=false);
        void DoAdvection_coeffVol(
            const TensorOfArray2D<NekDouble>      &inarray,
            Array<OneD, Array<OneD, NekDouble> >  &outarray,
            const NekDouble                       time,
            const bool                            flagFreezeJac);
        
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void MultiplyElmtInvMass_PlusSource(
            TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
            const NekDouble                         dtlamda,
            const DataType                          tmpDataType);

        void CalcFluxJacVolBnd(
            const TensorOfArray2D<NekDouble>  &inarray,
            const TensorOfArray3D<NekDouble>  &qfield);

        void GetFluxVectorMF(
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
            TensorOfArray3D<NekDouble>                 &flux);
        
        void GetFluxVectorTraceMF(
            const Array<OneD, Array<OneD, NekDouble> >              &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >              &pBwd,
            Array<OneD, Array<OneD, NekDouble> >                    &flux);
        void GetFluxVectorTraceMFNum(
            const Array<OneD, Array<OneD, NekDouble> >              &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >              &pBwd,
            Array<OneD, Array<OneD, NekDouble> >                    &flux);
        void GetFluxVectorTraceMFAna(
            const Array<OneD, Array<OneD, NekDouble> >              &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >              &pBwd,
            Array<OneD, Array<OneD, NekDouble> >                    &flux);

        void GetFluxVectorJacDirctn(
            const int                           nDirctn,
            const TensorOfArray2D<NekDouble>    &inarray,
            TensorOfArray5D<NekDouble>          &ElmtJacArray);
        
        void GetFluxVectorJacDirctnMat(
            const int                            nDirctn,
            const TensorOfArray2D<NekDouble>     &inarray,
            TensorOfArray2D<DNekBlkMatSharedPtr> &ElmtFluxJacArray);

        void GetFluxVectorJacDirctnElmt(
            const int                                    nConvectiveFields,
            const int                                    nElmtPnt,
            const Array<OneD, Array<OneD, NekDouble> >   &locVars,
            const Array<OneD, NekDouble>                 &normals,
            DNekMatSharedPtr                             &wspMat,
            Array<OneD, Array<OneD, NekDouble> >         &PntJacArray);

        void GetFluxVectorJacPoint(
            const int                                   nConvectiveFields,
            const Array<OneD, NekDouble>                &conservVar, 
            const Array<OneD, NekDouble>                &normals, 
            DNekMatSharedPtr                            &fluxJac);
        
        void CalTraceNumericalFlux(
            const int                                         nConvectiveFields,
            const int                                         nDim,
            const int                                         nPts,
            const int                                         nTracePts,
            const NekDouble                                   PenaltyFactor2,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &AdvVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            const NekDouble                                   time,
            const TensorOfArray3D<NekDouble>                  &qfield,
            const Array<OneD, Array<OneD, NekDouble> >        &vFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &vBwd,
            const TensorOfArray3D<NekDouble>                  &qFwd,
            const TensorOfArray3D<NekDouble>                  &qBwd,
            const Array<OneD, NekDouble >                     &MuVarTrace,
            Array<OneD, int >                                 &nonZeroIndex,
            Array<OneD, Array<OneD, NekDouble> >              &traceflux);
        
        void CalTraceIPSymmFlux(
            const int                                         nConvectiveFields,
            const int                                         nTracePts,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            const NekDouble                                   time,
            const TensorOfArray3D<NekDouble>                  &qfield,
            const Array<OneD, Array<OneD, NekDouble> >        &vFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &vBwd,
            const Array<OneD, NekDouble >                     &MuVarTrace,
            Array<OneD, int >                                 &nonZeroIndex,
            TensorOfArray3D<NekDouble>                        &traceflux);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void CalVisFluxDerivJac(
            const int                                         nConvectiveFields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            const Array<OneD, Array<OneD, NekDouble> >        &Fwd,
            const Array<OneD, Array<OneD, NekDouble> >        &Bwd,
            TypeNekBlkMatSharedPtr                            &BJac,
            DataType                                          &tmpDataType);

        void MinusDiffusionFluxJacDirctn(
            const int                          nDirctn,
            const TensorOfArray2D<NekDouble>   &inarray,
            const TensorOfArray3D<NekDouble>   &qfields,
            TensorOfArray5D<NekDouble>         &ElmtJacArray)
        {
            v_MinusDiffusionFluxJacDirctn(nDirctn,inarray, qfields,
                ElmtJacArray);
        }

        void MinusDiffusionFluxJacDirctnElmt(
            const int                                   nConvectiveFields,
            const int                                   nElmtPnt,
            const Array<OneD, Array<OneD, NekDouble> >  &locVars,
            const TensorOfArray3D<NekDouble>            &locDerv,
            const Array<OneD, NekDouble>                &locmu,
            const Array<OneD, NekDouble>                &locDmuDT,
            const Array<OneD, NekDouble>                &normals,
            DNekMatSharedPtr                            &wspMat,
            Array<OneD, Array<OneD, NekDouble> >        &PntJacArray)
        {
            v_MinusDiffusionFluxJacDirctnElmt(nConvectiveFields,nElmtPnt,
                locVars,locDerv,locmu,locDmuDT,normals,wspMat,PntJacArray);
        }
        void MinusDiffusionFluxJacDirctnMat(
            const int                            nDirctn,
            const TensorOfArray2D<NekDouble>     &inarray,
            const TensorOfArray3D<NekDouble>     &qfields,
            TensorOfArray2D<DNekBlkMatSharedPtr> &ElmtFluxJacArray)
        {
            v_MinusDiffusionFluxJacDirctnMat(nDirctn,inarray, qfields,
                ElmtFluxJacArray);
        }

        void GetFluxDerivJacDirctn(
            const MultiRegions::ExpListSharedPtr  &explist,
            const TensorOfArray2D<NekDouble>      &normals,
            const int                             nDervDir,
            const TensorOfArray2D<NekDouble>      &inarray,
            TensorOfArray5D<NekDouble>            &ElmtJacArray,
            const int                             nfluxDir)
        {
            v_GetFluxDerivJacDirctn(explist,normals,nDervDir,inarray,
                ElmtJacArray,nfluxDir);
        }

        void GetFluxDerivJacDirctnElmt(
            const int                                   nConvectiveFields,
            const int                                   nElmtPnt,
            const int                                   nDervDir,
            const Array<OneD, Array<OneD, NekDouble> >  &locVars,
            const Array<OneD, NekDouble>                &locmu,
            const Array<OneD, Array<OneD, NekDouble> >  &locnormal,
            DNekMatSharedPtr                            &wspMat,
            Array<OneD, Array<OneD, NekDouble> >        &PntJacArray)
        {
            v_GetFluxDerivJacDirctnElmt(nConvectiveFields,nElmtPnt,nDervDir,
                    locVars,locmu,locnormal,wspMat,PntJacArray);
        }
            
        void GetFluxDerivJacDirctn(
            const MultiRegions::ExpListSharedPtr        &explist,
            const TensorOfArray2D<NekDouble>            &normals,
            const int                                   nDervDir,
            const TensorOfArray2D<NekDouble>            &inarray,
            Array<OneD, Array<OneD, DNekMatSharedPtr> > &ElmtJac)
        {
            v_GetFluxDerivJacDirctn(explist,normals,nDervDir,inarray,ElmtJac);
        }
        // void GetFluxDerivJacDirctn(
        //     const MultiRegions::ExpListSharedPtr                            &explist,
        //     const int                                                       nFluxDir,
        //     const int                                                       nDervDir,
        //     const TensorOfArray2D<NekDouble>                &inarray,
        //           Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac)
        // {
        //     v_GetFluxDerivJacDirctn(explist,nFluxDir,nDervDir,inarray,ElmtJac);
        // }
        void CalphysDeriv(
            const TensorOfArray2D<NekDouble> &inarray,
            TensorOfArray3D<NekDouble>       &qfield)
        {
            v_CalphysDeriv(inarray, qfield);
        }
        void GetDiffusionFluxJacDirctn(
            const int                                   nDirctn,
            const TensorOfArray2D<NekDouble>            &inarray,
            const TensorOfArray3D<NekDouble>            &qfields,
            Array<OneD, Array<OneD, DNekMatSharedPtr> > &ElmtJac);
        void GetDiffusionFluxJacPoint(
            const Array<OneD, NekDouble>       &conservVar, 
            const TensorOfArray2D<NekDouble>   &conseDeriv, 
            const NekDouble                    mu,
            const NekDouble                    DmuDT,
            const Array<OneD, NekDouble>       &normals,
            DNekMatSharedPtr                   &fluxJac)
        {
            v_GetDiffusionFluxJacPoint(conservVar,conseDeriv,mu,DmuDT,normals,
                fluxJac);
        }
#endif

        void DoAdvection(
            const TensorOfArray2D<NekDouble>            &inarray,
            Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const NekDouble                             time,
            const Array<OneD, Array<OneD, NekDouble> >  &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >  &pBwd);

        void DoDiffusion(
            const TensorOfArray2D<NekDouble>            &inarray,
            Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >  &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >  &pBwd);
        void DoDiffusion_coeff(
            const TensorOfArray2D<NekDouble>            &inarray,
            Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >  &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >  &pBwd,
            const bool                                  flagFreezeJac = false);
        void DoDiffusion_coeffVol(
            const TensorOfArray2D<NekDouble>        &inarray,
            TensorOfArray2D<NekDouble>              &outarray,
            const bool                              flagFreezeJac = false)
        {
            v_DoDiffusion_coeffVol(inarray, outarray, flagFreezeJac);
        }

        void GetFluxVector(
            const Array<OneD, Array<OneD, NekDouble> >  &physfield,
            TensorOfArray3D<NekDouble>                  &flux);
        void GetFluxVectorDeAlias(
            const Array<OneD, Array<OneD, NekDouble> >  &physfield,
            TensorOfArray3D<NekDouble>                  &flux);

        void SetBoundaryConditions(
            Array<OneD, Array<OneD, NekDouble> >        &physarray,
            NekDouble                                    time);

        void SetBoundaryConditionsBwdWeight();

        void SetBoundaryConditionsDeriv(
            const TensorOfArray2D<NekDouble>  &physarray,
            const TensorOfArray3D<NekDouble>  &dervarray,
            NekDouble                         time,
            const TensorOfArray2D<NekDouble>  &pFwd     = 
                NullNekDoubleArrayofArray,
            const TensorOfArray3D<NekDouble>  &pDervFwd = 
                NullNekDoubleArrayofArrayofArray);

        void GetElmtTimeStep(
            const TensorOfArray2D<NekDouble> &inarray,
            Array<OneD, NekDouble>           &tstep);

        void GetViscousSymmtrFluxConservVar(
            const int                                  nConvectiveFields,
            const int                                  nSpaceDim,
            const Array<OneD, Array<OneD, NekDouble> > &inaverg,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            TensorOfArray3D<NekDouble>                 &outarray,
            Array< OneD, int >                         &nonZeroIndex,
            const Array<OneD, Array<OneD, NekDouble> > &normals)
        {
            v_GetViscousSymmtrFluxConservVar(nConvectiveFields,nSpaceDim,
                    inaverg,inarray,outarray,nonZeroIndex,normals);
        }

        void CalcMuDmuDT(
            const TensorOfArray2D<NekDouble> &inarray,
            Array<OneD, NekDouble>           &mu,
            Array<OneD, NekDouble>           &DmuDT)
        {
            v_CalcMuDmuDT(inarray,mu,DmuDT);
        }

        virtual NekDouble v_GetTimeStep(
            const TensorOfArray2D<NekDouble> &inarray);
        virtual void v_SetInitialConditions(
            NekDouble initialtime           = 0.0,
            bool      dumpInitialConditions = true,
            const int domain                = 0);

        NekDouble GetGamma()
        {
            return m_gamma;
        }

        const TensorOfArray2D<NekDouble> &GetVecLocs()
        {
            return m_vecLocs;
        }

        const TensorOfArray2D<NekDouble> &GetNormals()
        {
            return m_traceNormals;
        }

        virtual void v_ExtraFldOutput(
            std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
            std::vector<std::string>             &variables);

        virtual void v_DoDiffusion(
            const TensorOfArray2D<NekDouble>           &inarray,
            Array<OneD, Array<OneD, NekDouble> >       &outarray,
            const Array<OneD, Array<OneD, NekDouble> > &pFwd,
            const Array<OneD, Array<OneD, NekDouble> > &pBwd)
        {
            // Do nothing by default
        }

        virtual void v_DoDiffusion_coeff(
            const TensorOfArray2D<NekDouble>            &inarray,
            Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >  &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >  &pBwd,
            const bool                                  flagFreezeJac)
        {
            // Do nothing by default
        }

        virtual void v_DoDiffusion_coeffVol(
            const TensorOfArray2D<NekDouble>        &inarray,
            Array<OneD, Array<OneD, NekDouble> >    &outarray,
            const bool                              flagFreezeJac)
        {
            // Do nothing by default
        }

        virtual void v_DoDiffusionFlux(
            const TensorOfArray2D<NekDouble>            &inarray,
            TensorOfArray3D<NekDouble>                  &VolumeFlux,
            Array<OneD, Array<OneD, NekDouble>>         &TraceFlux,
            const Array<OneD, Array<OneD, NekDouble> >  &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >  &pBwd)
        {
            //Artificial Diffusion need to implement
            if (m_shockCaptureType != "Off")
            {
                m_artificialDiffusion->DoArtificialDiffusionFlux(inarray, 
                        VolumeFlux,TraceFlux);
            }
        }

        virtual Array<OneD, NekDouble> v_GetMaxStdVelocity(
            const NekDouble SpeedSoundFactor=1.0);

        virtual void v_GetViscousSymmtrFluxConservVar(
            const int                                   nConvectiveFields,
            const int                                   nSpaceDim,
            const Array<OneD, Array<OneD, NekDouble> >  &inaverg,
            const Array<OneD, Array<OneD, NekDouble> >  &inarray,
            TensorOfArray3D<NekDouble>                  &outarray,
            Array< OneD, int >                          &nonZeroIndex,
            const Array<OneD, Array<OneD, NekDouble> >  &normals);
        
        virtual void v_SteadyStateResidual(
                int                         step, 
                Array<OneD, NekDouble>      &L2);
        virtual void v_CalcMuDmuDT(
            const TensorOfArray2D<NekDouble> &inarray,
            Array<OneD, NekDouble>           &mu,
            Array<OneD, NekDouble>           &DmuDT)
        {
        }
                
#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
        virtual void v_GetDiffusionFluxJacPoint(
            const Array<OneD, NekDouble>        &conservVar, 
            const TensorOfArray2D<NekDouble>    &conseDeriv, 
            const NekDouble                     mu,
            const NekDouble                     DmuDT,
            const Array<OneD, NekDouble>        &normals,
            DNekMatSharedPtr                    &fluxJac);
        virtual void v_CalphysDeriv(
            const TensorOfArray2D<NekDouble>    &inarray,
            TensorOfArray3D<NekDouble>          &qfield)
        {}

        virtual void v_MinusDiffusionFluxJacDirctn(
            const int                          nDirctn,
            const TensorOfArray2D<NekDouble>   &inarray,
            const TensorOfArray3D<NekDouble>   &qfields,
            TensorOfArray5D<NekDouble>         &ElmtJacArray);
        virtual void v_MinusDiffusionFluxJacDirctnElmt(
            const int                                  nConvectiveFields,
            const int                                  nElmtPnt,
            const Array<OneD, Array<OneD, NekDouble> > &locVars,
            const TensorOfArray3D<NekDouble>           &locDerv,
            const Array<OneD, NekDouble>               &locmu,
            const Array<OneD, NekDouble>               &locDmuDT,
            const Array<OneD, NekDouble>               &normals,
            DNekMatSharedPtr                           &wspMat,
            Array<OneD, Array<OneD, NekDouble> >       &PntJacArray);

        virtual void v_MinusDiffusionFluxJacDirctnMat(
            const int                              nDirctn,
            const TensorOfArray2D<NekDouble>       &inarray,
            const TensorOfArray3D<NekDouble>       &qfields,
            TensorOfArray2D<DNekBlkMatSharedPtr>   &ElmtFluxJacArray);

        virtual void v_GetFluxDerivJacDirctn(
            const MultiRegions::ExpListSharedPtr &explist,
            const TensorOfArray2D<NekDouble>     &normals,
            const int                            nDervDir,
            const TensorOfArray2D<NekDouble>     &inarray,
            TensorOfArray5D<NekDouble>           &ElmtJacArray,
            const int                            nfluxDir);

        virtual void v_GetFluxDerivJacDirctnElmt(
            const int                                   nConvectiveFields,
            const int                                   nElmtPnt,
            const int                                   nDervDir,
            const Array<OneD, Array<OneD, NekDouble> >  &locVars,
            const Array<OneD, NekDouble>                &locmu,
            const Array<OneD, Array<OneD, NekDouble> >  &locnormal,
            DNekMatSharedPtr                            &wspMat,
            Array<OneD, Array<OneD, NekDouble> >        &PntJacArray);

        virtual void v_GetFluxDerivJacDirctn(
            const MultiRegions::ExpListSharedPtr            &explist,
            const TensorOfArray2D<NekDouble>                &normals,
            const int                                       nDervDir,
            const TensorOfArray2D<NekDouble>                &inarray,
            Array<OneD, Array<OneD, DNekMatSharedPtr> >     &ElmtJac);
        

        void v_MultiLevel(           
            const Array<OneD, NekDouble>                                    &inarray,
            Array<OneD, NekDouble>                                          &outarray, 
            const int                                                       Level,
            const int                                                       CurrentLevelCoeff,    
            const int                                                       LowLevelCoeff,     
            const bool                                                      UpDateOperatorflag);
        
        void v_CalculateNextLevelPreconditioner(
            const Array<OneD, const Array<OneD, NekDouble>>                 &inarrayCoeff,
            const NekDouble                                                 time,
            const NekDouble                                                 lambda);

        void RestrictResidual(
        const Array<OneD,DNekMatSharedPtr>                                  &RestrictionMatrix,
        const Array<OneD, NekDouble>                                        &inarray,
              Array<OneD, NekDouble>                                        &outarray);
        
        void RestrictSolution(
        const Array<OneD,DNekMatSharedPtr>                                  &RestrictionMatrix,
        const Array<OneD, const Array<OneD, NekDouble>>                     &inarray,
              Array<OneD, Array<OneD, NekDouble>>                           &outarray);
        
        void ProlongateSolution(
        const Array<OneD, DNekMatSharedPtr>                                 &ProlongationMatrix,
        const Array<OneD, NekDouble>                                        &inarray,
              Array<OneD, NekDouble>                                        &outarray);

        void PrintArray(Array<OneD, NekDouble> &Array);

        void OutputArray(Array<OneD, NekDouble> &Array);

        void OutputConstArray(const Array<OneD, NekDouble> &Array);

        void PrintMatrix(DNekMatSharedPtr &Matrix);

        void OutputMatrix(DNekMatSharedPtr &Matrix);

        void OutputConstMatrix(const DNekMatSharedPtr &Matrix);

        // virtual void v_GetFluxDerivJacDirctn(
        //     const MultiRegions::ExpListSharedPtr                            &explist,
        //     const int                                                       nFluxDir,
        //     const int                                                       nDervDir,
        //     const TensorOfArray2D<NekDouble>                &inarray,
        //           Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac);
#endif
    };
}
#endif
