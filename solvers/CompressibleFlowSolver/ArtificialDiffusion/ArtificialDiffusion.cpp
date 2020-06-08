///////////////////////////////////////////////////////////////////////////////
//
// File: ArtificialDiffusion.cpp
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
// Description: Abstract base class for compressible solver artificial diffusion
//              used for shock capturing artificial diffusion.
//
///////////////////////////////////////////////////////////////////////////////

#include "ArtificialDiffusion.h"

using namespace std;

namespace Nektar
{
ArtificialDiffusionFactory& GetArtificialDiffusionFactory()
{
    static ArtificialDiffusionFactory instance;
    return instance;
}

ArtificialDiffusion::ArtificialDiffusion(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const int spacedim)
        : m_session(pSession),
          m_fields(pFields)
{
    // Create auxiliary object to convert variables
    m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(
                m_session, spacedim);

    m_diffusion = SolverUtils::GetDiffusionFactory()
                                    .CreateInstance("LDG", "LDG");
    m_diffusion->SetArtificialDiffusionVector(
                        &ArtificialDiffusion::GetArtificialViscosity, this);
    m_diffusion->InitObject (m_session, m_fields);
}

/**
 *
 */
void ArtificialDiffusion::DoArtificialDiffusion(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD,       Array<OneD, NekDouble> > &outarray)
{
    v_DoArtificialDiffusion(inarray, outarray);
}

/**
 *
 */
void ArtificialDiffusion::v_DoArtificialDiffusion(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD,       Array<OneD, NekDouble> > &outarray)
{
    int i;
    int nvariables = inarray.num_elements();
    int npoints    = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);

    for (i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
    }

    m_diffusion->Diffuse(nvariables, m_fields, inarray, outarrayDiff);

    for (i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(npoints,
                    outarray[i], 1,
                    outarrayDiff[i], 1,
                    outarray[i], 1);
    }
}



/**
 *
 */
void ArtificialDiffusion::DoArtificialDiffusion_coeff(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD,       Array<OneD, NekDouble> > &outarray)
{
    v_DoArtificialDiffusion_coeff(inarray, outarray);
}


/**
 *
 */
void ArtificialDiffusion::v_DoArtificialDiffusion_coeff(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD,       Array<OneD, NekDouble> > &outarray)
{
    int i;
    int nvariables = inarray.num_elements();
    int npoints    = m_fields[0]->GetNpoints();
    int ncoeffs    = m_fields[0]->GetNcoeffs();

    Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);

    for (i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>(ncoeffs, 0.0);
    }

    m_diffusion->DiffuseCoeff(nvariables, m_fields, inarray, outarrayDiff);

    for (i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(ncoeffs,
                    outarray[i], 1,
                    outarrayDiff[i], 1,
                    outarray[i], 1);
    }
}
/**
 *
 */
//To DO, need to judge whether conservative/primal derivatives!
//To Do, SmoothCapture has not been modified because it includes  force.
void ArtificialDiffusion::v_DoArtificialDiffusionFlux(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>>&VolumeFlux,
    Array<OneD, Array<OneD, NekDouble>>             &TraceFlux)
{
    int nvariables = inarray.num_elements();
    int npoints    = m_fields[0]->GetNpoints();
    int nTracePts  = m_fields[0]->GetTrace()->GetTotPoints();
    int nDim       = m_fields[0]->GetCoordim(0);

    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> VolumeDiff(nDim);
    Array<OneD, Array<OneD, NekDouble>> TraceDiff(nvariables);
    Array<OneD,Array<OneD, Array<OneD, NekDouble>>> inarrayDiffderivative(nDim);
    
    for (int j = 0; j < nDim; ++j)
    {
        VolumeDiff[j] = Array<OneD, Array<OneD, NekDouble>>(nvariables);
        inarrayDiffderivative[j]=Array<OneD, Array<OneD, NekDouble>> (nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            VolumeDiff[j][i] = Array<OneD, NekDouble>(npoints, 0.0);
            inarrayDiffderivative[j][i]=Array<OneD, NekDouble>(npoints,0.0);
        }
    }
    for (int i = 0; i < nvariables; ++i)
    {
        TraceDiff[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
    }

    // Diffusion term in physical rhs form
    // To notice, needs to firstly calculate volumeflux, traceflux uses it.
    m_diffusion->DiffuseCalculateDerivative(nvariables,m_fields,inarray,inarrayDiffderivative);
    m_diffusion->DiffuseVolumeFlux(nvariables, m_fields, inarray,inarrayDiffderivative, VolumeFlux);
    m_diffusion->DiffuseTraceFlux(nvariables, m_fields, inarray,inarrayDiffderivative,VolumeFlux,TraceFlux);

    for (int j = 0; j < nDim; ++j)
    {
        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Vadd(npoints, &VolumeDiff[j][i][0], 1, &VolumeFlux[j][i][0], 1,
                        &VolumeFlux[j][i][0], 1);
        }
    }
    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(nTracePts, &TraceDiff[i][0], 1, &TraceFlux[i][0], 1,
                   &TraceFlux[i][0], 1);
    }
}

/**
 *
 */
void ArtificialDiffusion::GetArtificialViscosity(
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble  >             &mu)
{
    v_GetArtificialViscosity(physfield, mu);
}

}
