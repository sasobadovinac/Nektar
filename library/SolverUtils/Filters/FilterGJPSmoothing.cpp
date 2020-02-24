///////////////////////////////////////////////////////////////////////////////
//
// File FilterGJPSmoothing.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Gradient Jump Penalty Smoothing. 
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Filters/FilterGJPSmoothing.h>

namespace Nektar
{
namespace SolverUtils
{
std::string FilterGJPSmoothing::className =
        GetFilterFactory().RegisterCreatorFunction(
                "GJPSmoothing", FilterGJPSmoothing::create);

FilterGJPSmoothing::FilterGJPSmoothing(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem>      &pEquation,
    const ParamMap &pParams)
    : Filter(pSession, pEquation),
      m_index(0)
{
    // Load sampling frequency
    auto it = pParams.find("SmoothingFrequency");
    if (it == pParams.end())
    {
        m_smoothingFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(
            m_session->GetInterpreter(), it->second);
        m_smoothingFrequency = round(equ.Evaluate());
    }

    m_traceNormals = pEquation->GetTraceNormals();
}
    

FilterGJPSmoothing::~FilterGJPSmoothing()
{
}

void FilterGJPSmoothing::v_Initialise
(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
 const NekDouble &time)
{
    boost::ignore_unused(pFields,time);


    // Check this is a Mixed CG-DG expansion
    ASSERTL0(boost::iequals(m_session->GetSolverInfo("Projection"),"Mixed_CG_Discontinuous"),
             "Currently the GJP is only set up for mixed CG/DG expansions");

    m_shapeDim = pFields[0]->GetShapeDimension();
}

void FilterGJPSmoothing::v_Update
(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
 const NekDouble &time)
{
    boost::ignore_unused(pFields,time);
    
    m_index++;
    if (m_index % m_smoothingFrequency > 0)
    {
        return;
    }

    int ncoeffs = pFields[f]->GetNcoeffs();
    int nphys   = pFields[f]->GetTotPoints();
    Array<OneD, NekDouble> FilterCoeffs(ncoeffs);

    Array<OneD, Array<OneD, NekDouble> > deriv(3,NullNekDouble1DArray);
    for(int i = 0; i < m_shapeDim; ++i)
    {
        deriv[i] = Array<OneD, NekDouble> (nphys); 
    }

    Array<OneD, NekDOuble> Fwd(nTracePts), Bwd(nTracePts); 
    Array<OneD, NekDOuble> GradJumpOnTrace(nTracePts,0.0); 

    for(int f = 0; f < pFields.num_elements(); ++i)
    {
        // calculate derivative 
        pFields[f]->PhysDeriv(pFields[f]->GetPhys(), deriv[0],deriv[1],deriv[2]);

        // Evaluate the summation of normal derivative on the trace
        for(int n = 0; n < m_shapeDim; ++n)
        {
            pFields[f]->GetFwdBwdTracePhys(deriv[n],Fwd,Bwd,UseFwdValsOnBndry);

            // Multiply by normal and add to trace evaluation
            Vmath::Vsub(nTracePts,Fwd,1,Bwd,1,Fwd,1);
            Vmath::Vvtvp(nTracePts,Fwd,1,m_normaTracel[n],1,GradJumpOnTrace,1,GradJumpOnTrace,1);
        }

        // Add filter to coeffs and bwdtrans 
        pFields[f]->FwdTrans_IterPerElmt(pFields[f]->GetPhys(),FilterCoeffs);
        pFields[f]->AddTraceNormWRTDerivBasis(GradJumpOnTrace,FilterCoeffs);
        pFields[f]->BwdTrans(FilterCoeffs,pFields[f]->UpdatePhys());
    }
}

void FilterGJPSmoothing::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(pFields,time);
}

bool FilterGJPSmoothing::v_IsTimeDependent()
{
    return true;
}
}
}
