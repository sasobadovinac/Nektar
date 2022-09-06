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
};

typedef std::shared_ptr<ALEFixed> ALEFixedShPtr;
typedef std::shared_ptr<ALETranslate> ALETranslateShPtr;
typedef std::shared_ptr<ALERotate> ALERotateShPtr;
typedef std::shared_ptr<ALEPrescribe> ALEPrescribeShPtr;

} // namespace SolverUtils

} // namespace Nektar
#endif // NEKTAR_ALEHELPER_H
