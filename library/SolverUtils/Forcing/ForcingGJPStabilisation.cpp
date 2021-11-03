///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingGJPStabilisation.cpp
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
// Description: Body forcing
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Forcing/ForcingGJPStabilisation.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/AssemblyMap/LocTraceToTraceMap.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion2D.h>
#include <SolverUtils/EquationSystem.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
    std::string ForcingGJPStabilisation::classNameBody = GetForcingFactory().
        RegisterCreatorFunction("GJPStabilisation",
                                ForcingGJPStabilisation::create,
                                "Graient Jump Penalty Stablisation");
    std::string ForcingGJPStabilisation::classNameField = GetForcingFactory().
        RegisterCreatorFunction("GJPStabilization",
                                ForcingGJPStabilisation::create,
                                "Graient Jump Penalty Stablization");

    ForcingGJPStabilisation::ForcingGJPStabilisation(
                                                     const LibUtilities::SessionReaderSharedPtr &pSession,
                                                     const std::weak_ptr<EquationSystem>      &pEquation)
        : Forcing(pSession, pEquation)
    {
    }

    void ForcingGJPStabilisation::v_InitObject
    (const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
     const unsigned int& pNumForcingFields,
     const TiXmlElement* pForce)
    {

        const TiXmlElement *funcNameElmt = pForce->FirstChildElement("VELSCALING");

        if (funcNameElmt)
        {
            m_velScalingStr = funcNameElmt->GetText();
        }
        else
        {
            m_velScalingStr = std::string("DEFAULT");
        }

        
        if(m_session->GetComm()->GetRank() == 0)
        {
            NekDouble jumpScal; 
            m_session->LoadParameter("GJPJumpScale", jumpScal, 1.0);

            cout << "GJP Stabilisation:" << endl;
            cout << "\t Velocity-Scaling: " << m_velScalingStr << endl;
            cout << "\t jump-Scaling:    " << jumpScal << endl;
        }
        

        m_numForcingFields = pNumForcingFields;

        // set up GJPForcing 
        m_GJPData = MemoryManager<MultiRegions::GJPForcing>::
            AllocateSharedPtr(pFields[0]);

        m_traceNormals = m_GJPData->GetTraceNormals();
    }
    
    void ForcingGJPStabilisation::v_Apply
            (const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
             const Array<OneD, Array<OneD, NekDouble> > &inarray,
             Array<OneD, Array<OneD, NekDouble> > &outarray,
             const NekDouble &time)
    {
        boost::ignore_unused(time,pFields);

        int nTracePts = m_GJPData->GetNumTracePts(); 
        Array<OneD, NekDouble> unorm(nTracePts,1.0);
        
        // set up normal velocity if Specified
        if(boost::iequals(m_velScalingStr,"NormalVelocity"))
        {
            Array<OneD, NekDouble> Fwd(nTracePts), Bwd(nTracePts); 
            
            pFields[0]->GetFwdBwdTracePhys(inarray[0],Fwd,Bwd,true,true);
            Vmath::Vmul(nTracePts,Fwd,1,m_traceNormals[0],1,unorm,1);

            // Evaluate u.n on trace
            for(int f = 1; f < pFields[0]->GetCoordim(0); ++f)
            {
                pFields[0]->GetFwdBwdTracePhys(inarray[f],Fwd,Bwd,true,true);
                Vmath::Vvtvp(nTracePts,Fwd,1,m_traceNormals[f],1,
                             unorm,1,unorm,1);
            }
            Vmath::Vabs(nTracePts,unorm,1,unorm,1);
        }

        for(int f = 0; f < m_numForcingFields; ++f)
        {
            m_GJPData->Apply(inarray[f],outarray[f],false,unorm);
        }
    }
}
}
