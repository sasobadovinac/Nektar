#ifndef NEKTAR_ALEHELPER_H
#define NEKTAR_ALEHELPER_H

#include <MultiRegions/ExpList.h>

namespace Nektar
{

namespace SolverUtils
{

struct ALEBase;
typedef std::shared_ptr<ALEBase> ALEBaseShPtr;

class ALEHelper : public std::enable_shared_from_this<ALEHelper>
{
public:

    void InitObject(const SpatialDomains::MeshGraphSharedPtr &pGraph,
                    Array<OneD, MultiRegions::ExpListSharedPtr> &fields);

    void UpdateGridVelocity(NekDouble &time);

    inline Array<OneD, Array<OneD, NekDouble>> GetGridVelocity()
    {
        return m_gridVelocity;
    }

protected:
    Array<OneD, MultiRegions::ExpListSharedPtr> m_fieldsALE;
    std::map<int, int> m_elmtToExpId;
    Array<OneD, Array<OneD, NekDouble>> m_gridVelocity;
    std::vector<ALEBaseShPtr> m_ALEs;
};

struct ALEBase
{
    inline void UpdateGridVel(NekDouble time,
                              Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                              std::map<int, int> &elmtToExpId,
                              Array<OneD, Array<OneD, NekDouble>> &gridVelocity)
    {
        v_UpdateGridVel(time, fields, elmtToExpId, gridVelocity);
    }

private:
    inline virtual void v_UpdateGridVel(NekDouble time,
                                        Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                        std::map<int, int> &elmtToExpId,
                                        Array<OneD, Array<OneD, NekDouble>> &gridVelocity) = 0;
};

struct ALEFixed final : public ALEBase
{
    ALEFixed(SpatialDomains::ZoneBaseShPtr zone);

    virtual void v_UpdateGridVel(NekDouble time,
                                 Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                 std::map<int, int> &elmtToExpId,
                                 Array<OneD, Array<OneD, NekDouble>> &gridVelocity) final;

private:
    SpatialDomains::ZoneFixedShPtr m_zone;
};

struct ALETranslate final : public ALEBase
{
    ALETranslate(SpatialDomains::ZoneBaseShPtr zone);

    virtual void v_UpdateGridVel(NekDouble time,
                                 Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                 std::map<int, int> &elmtToExpId,
                                 Array<OneD, Array<OneD, NekDouble>> &gridVelocity) final;

private:
    SpatialDomains::ZoneTranslateShPtr m_zone;
};

struct ALERotate final : public ALEBase
{
    ALERotate(SpatialDomains::ZoneBaseShPtr zone);

    virtual void v_UpdateGridVel(NekDouble time,
                                 Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                 std::map<int, int> &elmtToExpId,
                                 Array<OneD, Array<OneD, NekDouble>> &gridVelocity) final;

private:
    SpatialDomains::ZoneRotateShPtr m_zone;
};

struct ALEPrescribe final : public ALEBase
{
    ALEPrescribe(SpatialDomains::ZoneBaseShPtr zone);

    virtual void v_UpdateGridVel(NekDouble time,
                                 Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                 std::map<int, int> &elmtToExpId,
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
