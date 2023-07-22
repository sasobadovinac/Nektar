///////////////////////////////////////////////////////////////////////////////
//
// File: NodalTetElec.h
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
// Description: Header file of 2D Nodal Triangle Fekete Points
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALTETELEC_H
#define NODALTETELEC_H

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/NodalUtil.h>
#include <memory>

namespace Nektar
{
namespace LibUtilities
{

class NodalTetElec : public Points<NekDouble>
{
public:
    virtual ~NodalTetElec()
    {
    }

    LIB_UTILITIES_EXPORT static std::shared_ptr<PointsBaseType> Create(
        const PointsKey &key);

    NodalTetElec(const PointsKey &key) : PointsBaseType(key)
    {
    }

protected:
    virtual const MatrixSharedPtrType v_GetI(const PointsKey &pkey) override
    {
        ASSERTL0(pkey.GetPointsDim() == 3,
                 "NodalTetElec Points can only interp to other 3d "
                 "point distributions");
        Array<OneD, const NekDouble> x, y, z;
        PointsManager()[pkey]->GetPoints(x, y, z);
        return GetI(x, y, z);
    }

    virtual const MatrixSharedPtrType v_GetI(
        const Array<OneD, const NekDouble> &x,
        const Array<OneD, const NekDouble> &y,
        const Array<OneD, const NekDouble> &z) override
    {
        size_t numpoints = x.size();
        size_t np        = GetTotNumPoints();

        Array<OneD, NekDouble> interp(np * numpoints);
        CalculateInterpMatrix(x, y, z, interp);

        NekDouble *d = interp.data();
        return MemoryManager<NekMatrix<NekDouble>>::AllocateSharedPtr(numpoints,
                                                                      np, d);
    }

private:
    static bool initPointsManager[];

    std::shared_ptr<NodalUtilTetrahedron> m_util;

    NodalTetElec()                           = delete;
    NodalTetElec(const NodalTetElec &points) = delete;

    void NodalPointReorder3d();

    virtual void v_CalculatePoints() override final;
    virtual void v_CalculateWeights() override final;
    virtual void v_CalculateDerivMatrix() override final;

    void CalculateInterpMatrix(const Array<OneD, const NekDouble> &xia,
                               const Array<OneD, const NekDouble> &yia,
                               const Array<OneD, const NekDouble> &zia,
                               Array<OneD, NekDouble> &interp);
};
} // namespace LibUtilities
} // namespace Nektar

#endif // NODALTETELEC_H
