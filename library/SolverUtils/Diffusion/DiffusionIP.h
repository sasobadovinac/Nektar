///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionIP.h
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
// Description: IP diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DIFFUSIONWEAKDG
#define NEKTAR_SOLVERUTILS_DIFFUSIONWEAKDG

#include <SolverUtils/Diffusion/Diffusion.h>
#define CFS_DEBUGMODE

namespace Nektar
{
    namespace SolverUtils
    {
        class DiffusionIP : public Diffusion
        {
        public:
            static DiffusionSharedPtr create(std::string diffType)
            {
                return DiffusionSharedPtr(new DiffusionIP());
            }
            
            static std::string type;

            void ConsVarAve(
                const int                                           nConvectiveFields,
                const int                                           npnts,
                const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
                const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                      Array<OneD,       Array<OneD, NekDouble> >    &aver);
            
        protected:
            DiffusionIP();
   		
            std::string                             m_shockCaptureType;
            NekDouble                               m_IPSymmFtluxCoeff;
            NekDouble                               m_IP2ndDervCoeff;
            NekDouble                               m_IPPenaltyCoeff;
	         
            Array<OneD, NekDouble>                            m_MuVarTrace;
            Array<OneD, Array<OneD, NekDouble> >              m_traceNormals;
            Array<OneD, Array<OneD, NekDouble> >              m_traceAver;
            Array<OneD, Array<OneD, NekDouble> >              m_traceJump;
            Array<OneD, NekDouble>                            m_tracBwdWeightAver;
            Array<OneD, NekDouble>                            m_tracBwdWeightJump;
            Array<OneD, NekDouble>                            m_traceNormDirctnElmtLength;
            Array<OneD, NekDouble>                            m_oIPPenaltyLength;
            LibUtilities::SessionReaderSharedPtr              m_session;

#ifdef CFS_DEBUGMODE
            int                                 m_DebugVolTraceSwitch; 
            int                                 m_DebugIP_DDGSwitch; 
#endif

            void GetPenaltyFactor(
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                Array<OneD, NekDouble >                             &factor); 
            void GetPenaltyFactor_const(
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                Array<OneD, NekDouble >                             &factor); 

            virtual void v_AddSymmFluxIntegralToCoeff(
                const int                                         nvariables,
                const int                                         nDim,
                const int                                         nPts,
                const int                                         nTracePts,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, const int >                     &nonZeroIndex,
                TensorOfArray3D<NekDouble>                        &tracflux,
                TensorOfArray2D<NekDouble>                        &outarray);

            void AddSymmFluxIntegralToPhys(
                const int                                                           nConvectiveFields,
                const int                                                           nDim,
                const int                                                           nPts,
                const int                                                           nTracePts,
                const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
                const Array<OneD, const int >                                       &nonZeroIndex,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &tracflux,
                      Array<OneD, Array<OneD, NekDouble> >                          &outarray);
            
            void CalTraceSymFlux(
                const int                                                           nConvectiveFields,
                const int                                                           nDim,
                const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
                const Array<OneD, Array<OneD, NekDouble> >                          &solution_Aver,
                const Array<OneD, Array<OneD, NekDouble> >                          &solution_jump,
                    Array<OneD, int >                                             &nonZeroIndexsymm,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &traceSymflux);
            
            virtual void v_DiffuseTraceSymmFlux(
                const int                                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
                const Array<OneD, Array<OneD, NekDouble>>                   &inarray,
                const Array<OneD,Array<OneD, Array<OneD, NekDouble> > >     &qfield,
                const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    &VolumeFlux,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >          &SymmFlux,
                const Array<OneD, Array<OneD, NekDouble>>                   &pFwd,
                const Array<OneD, Array<OneD, NekDouble>>                   &pBwd,
                Array< OneD, int >                                          &nonZeroIndex,
                Array<OneD, Array<OneD, NekDouble> >                        &solution_Aver,
                Array<OneD, Array<OneD, NekDouble> >                        &solution_jump);
            
            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);
            
            virtual void v_Diffuse(
                const int                                          nConvective,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);
            
            virtual void v_DiffuseCoeff(
                const int                                          nConvective,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);

            virtual void v_DiffuseCoeff(
                const int                                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
                const Array<OneD, Array<OneD, NekDouble> >                  &inarray,
                Array<OneD, Array<OneD, NekDouble> >                        &outarray,
                const Array<OneD, Array<OneD, NekDouble> >                  &vFwd,
                const Array<OneD, Array<OneD, NekDouble> >                  &vBwd,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >          &qfield,
                Array< OneD, int >                                          &nonZeroIndex);

            virtual void v_DiffuseCoeffVol(
                const int                                          nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>  &fields,
                const Array<OneD, Array<OneD, NekDouble> >         &inarray,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                Array<OneD, Array<OneD, NekDouble> >               &outarray,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &elmtFlux,
                Array< OneD, int >                                 &nonZeroIndex);

            virtual void v_DiffuseCoeffTrac(
                const int                                          nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>  &fields,
                const Array<OneD, Array<OneD, NekDouble> >         &inarray,
                Array<OneD, Array<OneD, NekDouble> >               &outarray,
                const Array<OneD, Array<OneD, NekDouble> >         &vFwd,
                const Array<OneD, Array<OneD, NekDouble> >         &vBwd,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &elmtFlux,
                Array< OneD, int >                                 &nonZeroIndex);

            virtual void v_DiffuseVolumeFlux(
                const int                                           nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array< OneD, int >                                  &nonZeroIndex) ;
            
            virtual void v_DiffuseTraceFlux(
                const int                                           nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &TraceFlux,
                const Array<OneD, Array<OneD, NekDouble>>           &pFwd,
                const Array<OneD, Array<OneD, NekDouble>>           &pBwd,
                Array< OneD, int >                                  &nonZeroIndex);
            void v_DiffuseTraceFlux(
                const int                                                       nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>               &fields,
                const Array<OneD, Array<OneD, NekDouble>>                       &inarray,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                Array<OneD, Array<OneD, NekDouble> >                            &TraceFlux,
                const Array<OneD, Array<OneD, NekDouble>>                       &pFwd,
                const Array<OneD, Array<OneD, NekDouble>>                       &pBwd,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qFwd,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qBwd,
                const Array<OneD, NekDouble>                                    &MuAVTrace,
                Array< OneD, int >                                              &nonZeroIndex  ,
                const Array<OneD, Array<OneD, NekDouble>>                       &Aver          ,
                const Array<OneD, Array<OneD, NekDouble>>                       &Jump          );
            virtual void v_AddDiffusionSymmFluxToCoeff(
                const int                                           nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble> >          &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &outarray,
                const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >          &pBwd);

            // virtual void v_AddDiffusionSymmFluxToPhys(
            //     const int                                           nConvectiveFields,
            //     const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
            //     const Array<OneD, Array<OneD, NekDouble> >          &inarray,
            //     Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
            //     Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
            //     Array<OneD, Array<OneD, NekDouble> >                &outarray,
            //     const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
            //     const Array<OneD, Array<OneD, NekDouble> >          &pBwd);

            virtual void v_DiffuseCalculateDerivative(
                const int                                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
                const Array<OneD, Array<OneD, NekDouble> >                  &inarray,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    &qfield,
                const Array<OneD, Array<OneD, NekDouble> >                  &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >                  &pBwd);
           
            virtual void v_ConsVarAveJump(
                const int                                           nConvectiveFields,
                const int                                           npnts,
                const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
                const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                      Array<OneD,       Array<OneD, NekDouble> >    &aver,
                      Array<OneD,       Array<OneD, NekDouble> >    &jump);
            
            virtual const Array<OneD, const Array<OneD, NekDouble> > &v_GetTraceNormal()
            {
                return m_traceNormals;
            }

            void CalTraceNumFlux_ReduceComm(
                const int                                                           nConvectiveFields,
                const int                                                           nDim,
                const int                                                           nPts,
                const int                                                           nTracePts,
                const NekDouble                                                     PenaltyFactor2,
                const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
                const Array<OneD, Array<OneD, NekDouble> >                          &inarray,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qfield,
                const Array<OneD, Array<OneD, NekDouble> >                          &vFwd,
                const Array<OneD, Array<OneD, NekDouble> >                          &vBwd,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qFwd,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qBwd,
                const Array<OneD, NekDouble >                                       &MuVarTrace,
                      Array<OneD, int >                                             &nonZeroIndexflux,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &traceflux,
                      Array<OneD, Array<OneD, NekDouble> >                          &solution_Aver,
                      Array<OneD, Array<OneD, NekDouble> >                          &solution_jump);

            void ApplyFluxBndConds(
                const int                                               nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>       &fields,
                Array<OneD,       Array<OneD, NekDouble> >              &flux);
            
           void AddSecondDerivTOTrace_ReduceComm(
                const int                                                           nConvectiveFields,
                const int                                                           nDim,
                const int                                                           nPts,
                const int                                                           nTracePts,
                const NekDouble                                                     PenaltyFactor2,
                const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qfield,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >              &numDerivFwd,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >              &numDerivBwd);

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
            virtual void v_MinusVolumDerivJacToMat( 
                const int                                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &pFields,
                const Array<OneD, const Array<OneD,  Array<OneD, 
                    Array<OneD,  Array<OneD,  NekDouble> > > > >            &ElmtJacArray,
                const int                                                   nDervDir, 
                Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >              &gmtxarray);
            virtual void v_MinusVolumDerivJacToMat( 
                const int                                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &pFields,
                const Array<OneD, const Array<OneD,  Array<OneD, 
                    Array<OneD,  Array<OneD,  NekDouble> > > > >            &ElmtJacArray,
                const int                                                   nDervDir, 
                Array<OneD, Array<OneD, SNekBlkMatSharedPtr> >              &gmtxarray);
            
#endif
        }; 
    }
}
    
#endif
