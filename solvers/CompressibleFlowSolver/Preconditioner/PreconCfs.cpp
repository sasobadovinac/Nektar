///////////////////////////////////////////////////////////////////////////////
//
// File:  PreconCfs.cpp
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
// Description:  PreconCfs definition
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/Preconditioner/PreconCfs.h>

using namespace std;

namespace Nektar
{
    PreconCfs::PreconCfs(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm)
    {
        m_comm = vComm;

        pSession->LoadParameter("PreconMatFreezNumb",     
            m_PreconMatFreezNumb    , 1);
    }
    void PreconCfs::v_InitObject()
    {

    }

    void PreconCfs::DoNullPrecon(
            const Array<OneD, NekDouble> &pInput,
            Array<OneD, NekDouble> &pOutput,
            const bool &flag)
    {
        boost::ignore_unused(flag);
        Vmath::Vcopy(pInput.size(), pInput, 1, pOutput, 1);
    }

    /**
     *
     */
    inline void PreconCfs::DoPreconCfs(
        const Array<OneD, NekDouble> &pInput,
        Array<OneD, NekDouble> &pOutput,
        const bool &flag)
    {
        ASSERTL0(pInput.size() == pOutput.size(), 
            "In and Out not the same size in DoNullPrecon");
        v_DoPreconCfs(pInput, pOutput, flag);
    }


    void PreconCfs::v_DoPreconCfs(
        const Array<OneD, NekDouble> &pInput,
        Array<OneD, NekDouble> &pOutput,
        const bool &flag)
    {
        NEKERROR(ErrorUtil::efatal, "v_DoPreconCfs not defined");
    }

    void PreconCfs::v_BuildPreconCfs()
    {
        NEKERROR(ErrorUtil::efatal, "v_BuildPreconCfs not defined");
    }

} // namespace Nektar