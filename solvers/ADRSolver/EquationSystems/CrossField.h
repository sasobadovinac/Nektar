///////////////////////////////////////////////////////////////////////////////
//
// File CrossField.h
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
// Description: CrossField solve routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_CROSSFIELD_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_CROSSFIELD_H

#include <ADRSolver/EquationSystems/Laplace.h>

namespace Nektar
{
class CrossField : public Laplace
{
public:
    friend class MemoryManager<CrossField>;

    /// Creates an instance of this class
    static EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        EquationSystemSharedPtr p =
            MemoryManager<CrossField>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }
    /// Name of class
    static std::string className;

    virtual ~CrossField();

protected:
    CrossField(const LibUtilities::SessionReaderSharedPtr &pSession,
               const SpatialDomains::MeshGraphSharedPtr &pGraph);

    void CrossField::v_DoSolve()
};
}

#endif
