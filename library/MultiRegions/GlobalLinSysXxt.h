///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysXxt.h
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

/*
 * GlobalLinSysXxt.h
 *
 *  Created on: 19 Oct 2012
 *      Author: cc
 */

#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSXXT_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSXXT_H
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/GlobalLinSysKey.h>
#include <MultiRegions/MultiRegionsDeclspec.h>

namespace Xxt
{
struct crs_data;
}

namespace Nektar
{
namespace MultiRegions
{
// Forward declarations

// class AssemblyMapDG;
class ExpList;

class GlobalLinSysXxt : virtual public GlobalLinSys
{
public:
    /// Constructor for full direct matrix solve.
    MULTI_REGIONS_EXPORT GlobalLinSysXxt(
        const GlobalLinSysKey &pKey, const std::weak_ptr<ExpList> &pExp,
        const std::shared_ptr<AssemblyMap> &pLocToGloMap);

    MULTI_REGIONS_EXPORT virtual ~GlobalLinSysXxt();

protected:
    struct Xxt::crs_data *m_crsData;
    Array<OneD, unsigned int> m_Ai;
    Array<OneD, unsigned int> m_Aj;
    Array<OneD, double> m_Ar;

    Array<OneD, NekDouble> m_locToGloSignMult;

    Array<OneD, int> m_map;

    /// Solve the linear system for given input and output vectors.
    virtual void v_SolveLinearSystem(const int pNumRows,
                                     const Array<OneD, const NekDouble> &pInput,
                                     Array<OneD, NekDouble> &pOutput,
                                     const AssemblyMapSharedPtr &locToGloMap,
                                     const int pNumDir = 0);

    void GlobalToLocalNoSign(const Array<OneD, const NekDouble> &global,
                             Array<OneD, NekDouble> &local);

    void LocalToGlobalNoSign(const Array<OneD, const NekDouble> &local,
                             Array<OneD, NekDouble> &global);
};
} // namespace MultiRegions
} // namespace Nektar
#endif /* GLOBALLINSYSXXT_H_ */
