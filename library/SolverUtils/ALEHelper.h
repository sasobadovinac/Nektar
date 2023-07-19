///////////////////////////////////////////////////////////////////////////////
//
// File: ALEHelper.h
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
// Description: Helper class for ALE process
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_ALEHELPER_H
#define NEKTAR_ALEHELPER_H

#include <MultiRegions/ExpList.h>

namespace Nektar
{

namespace SolverUtils
{

struct ALEBase;
typedef std::shared_ptr<ALEBase> ALEBaseShPtr;

class ALEHelper
{
public:
    void InitObject(int spaceDim,
                    Array<OneD, MultiRegions::ExpListSharedPtr> &fields);

    void UpdateGridVelocity(const NekDouble &time);

    void ALEPreMultiplyMass(Array<OneD, Array<OneD, NekDouble>> &fields);

    void ALEDoElmtInvMass(Array<OneD, Array<OneD, NekDouble>> &traceNormals,
                          Array<OneD, Array<OneD, NekDouble>> &fields,
                          NekDouble time);

    void ALEDoElmtInvMassBwdTrans(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    void MoveMesh(const NekDouble &time,
                  Array<OneD, Array<OneD, NekDouble>> &traceNormals);

    inline const Array<OneD, const Array<OneD, NekDouble>> &GetGridVelocity()
    {
        return m_gridVelocity;
    }

    const Array<OneD, const Array<OneD, NekDouble>> &GetGridVelocityTrace();

protected:
    Array<OneD, MultiRegions::ExpListSharedPtr> m_fieldsALE;
    Array<OneD, Array<OneD, NekDouble>> m_gridVelocity;
    Array<OneD, Array<OneD, NekDouble>> m_gridVelocityTrace;
    std::vector<ALEBaseShPtr> m_ALEs;
    bool m_ALESolver          = false;
    NekDouble m_prevStageTime = 0.0;
    int m_spaceDim;
};

struct ALEBase
{
    virtual ~ALEBase() = default;

    inline void UpdateGridVel(
        const NekDouble time,
        Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        Array<OneD, Array<OneD, NekDouble>> &gridVelocity)
    {
        v_UpdateGridVel(time, fields, gridVelocity);
    }

private:
    virtual void v_UpdateGridVel(
        const NekDouble time,
        Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        Array<OneD, Array<OneD, NekDouble>> &gridVelocity) = 0;
};

struct ALEFixed final : public ALEBase
{
    ALEFixed(SpatialDomains::ZoneBaseShPtr zone);

    virtual void v_UpdateGridVel(
        const NekDouble time,
        Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        Array<OneD, Array<OneD, NekDouble>> &gridVelocity) final;

private:
    SpatialDomains::ZoneFixedShPtr m_zone;
};

struct ALETranslate final : public ALEBase
{
    ALETranslate(SpatialDomains::ZoneBaseShPtr zone);

    virtual void v_UpdateGridVel(
        const NekDouble time,
        Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        Array<OneD, Array<OneD, NekDouble>> &gridVelocity) final;

private:
    SpatialDomains::ZoneTranslateShPtr m_zone;
};

struct ALERotate final : public ALEBase
{
    ALERotate(SpatialDomains::ZoneBaseShPtr zone);

    virtual void v_UpdateGridVel(
        const NekDouble time,
        Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        Array<OneD, Array<OneD, NekDouble>> &gridVelocity) final;

private:
    SpatialDomains::ZoneRotateShPtr m_zone;
};

struct ALEPrescribe final : public ALEBase
{
    ALEPrescribe(SpatialDomains::ZoneBaseShPtr zone);

    virtual void v_UpdateGridVel(
        const NekDouble time,
        Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        Array<OneD, Array<OneD, NekDouble>> &gridVelocity) final;

private:
    SpatialDomains::ZonePrescribeShPtr m_zone;
    std::map<int, std::map<int, Array<OneD, NekDouble>>> m_coords;
};

typedef std::shared_ptr<ALEFixed> ALEFixedShPtr;
typedef std::shared_ptr<ALETranslate> ALETranslateShPtr;
typedef std::shared_ptr<ALERotate> ALERotateShPtr;
typedef std::shared_ptr<ALEPrescribe> ALEPrescribeShPtr;

} // namespace SolverUtils

} // namespace Nektar
#endif // NEKTAR_ALEHELPER_H
