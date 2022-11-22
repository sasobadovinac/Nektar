///////////////////////////////////////////////////////////////////////////////
//
// File: PanditGilesDemir03.h
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
// Description: Pandit-Giles-Demir 2003 cell model
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_PANDIT_GILES_DEMIR_03_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_PANDIT_GILES_DEMIR_03_H

#include <CardiacEPSolver/CellModels/CellModel.h>
namespace Nektar
{
class PanditGilesDemir03 : public CellModel
{

public:
    /// Creates an instance of this class
    static CellModelSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession, const int nq)
    {
        return MemoryManager<PanditGilesDemir03>::AllocateSharedPtr(pSession,
                                                                    nq);
    }

    /// Name of class
    static std::string className;

    /// Constructor
    PanditGilesDemir03(const LibUtilities::SessionReaderSharedPtr &pSession,
                       const int nq);

    /// Desctructor
    virtual ~PanditGilesDemir03()
    {
    }

protected:
    /// Computes the reaction terms $f(u,v)$ and $g(u,v)$.
    virtual void v_Update(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const NekDouble time) override;

    /// Prints a summary of the model parameters.
    virtual void v_GenerateSummary(SummaryList &s) override;
};

} // namespace Nektar

#endif
