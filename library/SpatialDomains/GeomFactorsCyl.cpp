////////////////////////////////////////////////////////////////////////////////
//
//  File: GeomFactorsCyl.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Geometric factors Cylindrical coordinates.
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/GeomFactorsCyl.h>

namespace Nektar
{
namespace SpatialDomains
{

DerivStorage GeomFactorsCyl::ComputeDeriv(
    const LibUtilities::PointsKeyVector &keyTgt) const
{
    DerivStorage deriv = GeomFactors::ComputeDeriv(keyTgt);

    /*
    for (int i = 0; i < m_expDim; ++i)
    {
        for (int j = 0; j < m_coordDim; ++j)
        {
            MultByRadius(&deriv[i][j][0], keyTgt);
        }
    }
    */

    return deriv;
}

Array<OneD, NekDouble> GeomFactorsCyl::ComputeJac(
    const LibUtilities::PointsKeyVector &keyTgt) const
{
    Array<OneD, NekDouble> jac = GeomFactors::ComputeJac(keyTgt);
    MultByRadius(jac, keyTgt);
    return jac;
}

void GeomFactorsCyl::MultByRadius(
    Array<OneD, NekDouble> &work,
    const LibUtilities::PointsKeyVector &keyTgt) const
{
    LibUtilities::PointsKeyVector xmapKeys = m_xmap->GetPointsKeys();
    int j;

    int nTotPts = 1;

    for (j = 0; j < keyTgt.size(); ++j)
    {
        nTotPts *= keyTgt[j].GetNumPoints();
    }

    ASSERTL0(nTotPts == work.num_elements(), "should be equal");

    // Backward transform y-component to get radius.
    Array<OneD, NekDouble> tmp(m_xmap->GetTotPoints()), tmp2(nTotPts);
    m_xmap->BwdTrans(m_coords[1], tmp);

    // Interpolate onto target distribution
    Interp(xmapKeys, tmp, keyTgt, tmp2);
    Vmath::Vmul(nTotPts, &work[0], 1, &tmp2[0], 1, &work[0], 1);
}

}
}
