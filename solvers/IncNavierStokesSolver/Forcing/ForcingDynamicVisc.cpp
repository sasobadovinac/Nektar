///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingDynamicVisc.cpp
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
// Description: Dynamics Forcing explicit forcing 
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Forcing/ForcingDynamicVisc.h>
#include <MultiRegions/ExpList.h>

using namespace std;

namespace Nektar
{
std::string ForcingDynamicVisc::className = SolverUtils::GetForcingFactory().
            RegisterCreatorFunction("DynamicVisc",
                                    ForcingDynamicVisc::create,
                                    "Moving Body Forcing");

ForcingDynamicVisc::ForcingDynamicVisc(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    : Forcing(pSession)
{
}

void ForcingDynamicVisc::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
        const unsigned int& pNumForcingFields,
        const TiXmlElement* pForce)
{
    int phystot = pFields[0]->GetTotPoints();
    m_kinvis = Array<OneD, NekDouble>(phystot,0.0);

    // assume number of fields is same as expansion dimension
    m_dim = pNumForcingFields;

    m_diff = Array<OneD, Array<OneD, NekDouble> > (3);
    for(int i = 0; i < m_dim; i++)
    {
        m_diff[i] = Array<OneD, NekDouble>(phystot);
    }
}


/*
  Laplacian d/dx kinvis d/dx + d/dy kinvis d/dy + d/dxz kinvis d/dz 
*/
void ForcingDynamicVisc::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>&  pFields,
        const Array<OneD, Array<OneD, NekDouble> >&         inarray,
              Array<OneD, Array<OneD, NekDouble> >&         outarray,
        const NekDouble&                                    time)
{

    int nquad = pFields[0]->GetTotPoints();
    Array<OneD, NekDouble> Lap(nquad);
    
    for(int i = 0; i < m_dim; ++i)
    {
        // Grad 
        if(m_dim == 2)
        {
            pFields[i]->PhysDeriv(inarray[i],m_diff[0],m_diff[1]);
        }
        else
        {
            pFields[i]->PhysDeriv(inarray[i],m_diff[0],m_diff[1],m_diff[2]);
        }

        Vmath::Zero(nquad,Lap,1);
        //Div kivis
        for(int j = 0; j < m_dim; ++j)
        {
            Vmath::Vmul(nquad,m_kinvis,1,m_diff[j],1,m_diff[j],1);
            pFields[i]->PhysDeriv(j,m_diff[j],m_diff[0]);
            Vmath::Vadd(nquad,m_diff[0],1,Lap,1,Lap,1);
        }
        
        Vmath::Vadd(nquad,Lap,1,outarray[i],1,outarray[i],1);
    }
}

void ForcingDynamicVisc::SetKinvis(Array<OneD, NekDouble> &kinvis)
{
    m_kinvis = kinvis;
}
    
Array<OneD, const NekDouble> &ForcingDynamicVisc::GetKinvis()
{
    return  m_kinvis;
}
    
}
