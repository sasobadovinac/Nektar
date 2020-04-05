///////////////////////////////////////////////////////////////////////////////
//
// File: Diffusion.cpp
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
// Description: Abstract base class for diffusion.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar
{
    namespace SolverUtils
    {
        DiffusionFactory& GetDiffusionFactory()
        {
            static DiffusionFactory instance;
            return instance;
        }
        
        void Diffusion::InitObject(
            const LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
        {
            v_InitObject(pSession, pFields);
        }
        
        void Diffusion::Diffuse(
            const int nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            v_Diffuse(nConvectiveFields, fields, inarray, outarray, pFwd, pBwd);
        }

        void Diffusion::Diffuse_coeff(
            const int nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            v_Diffuse_coeff(nConvectiveFields, fields, inarray, outarray, pFwd, pBwd);
        }

        void Diffusion::Diffuse(
            const int nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray,
            NekDouble                                         time,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            m_time  =    time;
            v_Diffuse(nConvectiveFields, fields, inarray, outarray, pFwd, pBwd);
        }

        void Diffusion::Diffuse_coeff(
            const int nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray,
            NekDouble                                         time,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd,
            const bool                                        flagFreezeJac)
        {
            m_time  =    time;
            m_flagFreezeJac = flagFreezeJac;
            v_Diffuse_coeff(nConvectiveFields, fields, inarray, outarray, pFwd, pBwd);
            m_flagFreezeJac = false;
        }
        void Diffusion::GetAVmu(
            const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
            const Array<OneD, Array<OneD, NekDouble> >                  &inarray,
                  Array<OneD, NekDouble >                               &muvar,
                  Array<OneD, NekDouble >                               &MuVarTrace)
        {
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();

            Array<OneD, NekDouble> Fwd(nTracePts,0.0);
            Array<OneD, NekDouble> Bwd(nTracePts,0.0);
            
            m_ArtificialDiffusionVector(inarray, muvar);

            // BwdMuvar is left to be 0.0 according to DiffusionLDG.cpp
            fields[0]->GetFwdBwdTracePhysNoBndFill(muvar,Fwd,Bwd);

            for(int k = 0; k < nTracePts; ++k)
            {
                MuVarTrace[k] = 0.5 * (Fwd[k] + Bwd[k]) ;
            }
        }
        
        void Diffusion::v_Diffuse_coeff(
            const int nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> > &pFwd,
            const Array<OneD, Array<OneD, NekDouble> > &pBwd)
        {
            ASSERTL0(false,"v_Diffuse_coeff not defined");
        }

        void Diffusion::v_Diffuse_coeff(
            const int                                                   nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
            const Array<OneD, Array<OneD, NekDouble> >                  &inarray,
            Array<OneD, Array<OneD, NekDouble> >                        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >                  &vFwd,
            const Array<OneD, Array<OneD, NekDouble> >                  &vBwd,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >          &qfield,
            Array< OneD, int >                                          &nonZeroIndex)
        {
            ASSERTL0(false,"v_Diffuse_coeff not defined");
        }

        const Array<OneD, const Array<OneD, NekDouble> > &Diffusion::v_GetTraceNormal()
        {
            ASSERTL0(false," not defined");
        }

        void Diffusion::v_ConsVarAveJump(
                const int                                           nConvectiveFields,
                const int                                           npnts,
                const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
                const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                      Array<OneD,       Array<OneD, NekDouble> >    &aver,
                      Array<OneD,       Array<OneD, NekDouble> >    &jump)
        {
            ASSERTL0(false," not defined");
        }
        
        // No multiply(check if diffsionVolume difined)
        void Diffusion::v_DiffuseCalculateDerivative(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble>>         &inarray,
            Array<OneD,Array<OneD, Array<OneD, NekDouble> > > &inarrayderivative,
            const Array<OneD, Array<OneD, NekDouble>>         &pFwd,
            const Array<OneD, Array<OneD, NekDouble>>         &pBwd)
        {
            ASSERTL0(false, "Not defined for function DiffuseVolumeFLux.");
        }

        // No multiply(check if diffsionVolume difined)
        void Diffusion::v_DiffuseVolumeFlux(
            const int                                           nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
            const Array<OneD, Array<OneD, NekDouble>>           &inarray,
            Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
            Array< OneD, int >                                  &nonZeroIndex)       
        {
            ASSERTL0(false, "Not defined for function DiffuseVolumeFLux.");
        }

        // No multiply(check if diffsionTraceFlux difined)
        void Diffusion::v_DiffuseTraceFlux(
            const int                                           nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
            const Array<OneD, Array<OneD, NekDouble>>           &inarray,
            Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
            Array<OneD, Array<OneD, NekDouble> >                &TraceFlux,
            const Array<OneD, Array<OneD, NekDouble>>           &pFwd,
            const Array<OneD, Array<OneD, NekDouble>>           &pBwd,
            Array< OneD, int >                                  &nonZeroIndex)     
        {
            ASSERTL0(false, "Not defined function DiffuseTraceFLux.");
        }

        void Diffusion::v_DiffuseTraceFlux(
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
                const Array<OneD, Array<OneD, NekDouble>>                       &Jump          )
            {
                ASSERTL0(false, "Not defined function DiffuseTraceFLux.");
            }
        
#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
        void Diffusion::v_MinusVolumDerivJacToMat( 
            const int                                                   nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>           &pFields,
            const Array<OneD, const Array<OneD,  Array<OneD, 
                Array<OneD,  Array<OneD,  NekDouble> > > > >            &ElmtJacArray,
            const int                                                   nDervDir, 
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >              &gmtxarray)
        {
            ASSERTL0(false," not defined");
        }
#endif
    }
}
