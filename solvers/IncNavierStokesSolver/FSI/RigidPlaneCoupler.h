///////////////////////////////////////////////////////////////////////////////
//
// File: RigidPlaneCoupler.h
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
// Description: FSICoupler with rigid 3DH1D planes deformation
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_FSI_RIGIDPLANECOUPLER
#define NEKTAR_SOLVERS_FSI_RIGIDPLANECOUPLER

#include <string>
#include <IncNavierStokesSolver/FSI/FSICoupler.h>

namespace Nektar
{

class RigidPlaneCoupler: public FSICoupler
{
public:

    friend class MemoryManager<RigidPlaneCoupler>;

    /// Creates an instance of this class
    static FSICouplerSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr        &pSession,
        const std::weak_ptr<SolverUtils::EquationSystem>   &pEquation,
        TiXmlElement                                *pFSI)
    {
        FSICouplerSharedPtr p =
                MemoryManager<RigidPlaneCoupler>::AllocateSharedPtr(pSession,
                                                                 pEquation);
        p->InitObject(pEquation.lock()->UpdateFields(), pFSI);
        return p;
    }

    ///Name of the class
    static std::string className;

protected:
    // Rank of process which will broadcast the displacement
    int     m_bcastRank;
    // Id of boundary from which we have to get the displacement
    int     m_bndId;
    // Constructor
    RigidPlaneCoupler(const LibUtilities::SessionReaderSharedPtr        &pSession,
              const std::weak_ptr<SolverUtils::EquationSystem>   &pEquation);


    // Virtual functions
    virtual void v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields,
              TiXmlElement* pFSI);

    virtual void v_CalculateDisplacement();

};

}

#endif
