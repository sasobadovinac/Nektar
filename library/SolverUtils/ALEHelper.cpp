#include "ALEHelper.h"

namespace Nektar
{

namespace SolverUtils
{

void ALEHelper::InitObject(const SpatialDomains::MeshGraphSharedPtr &pGraph,
                     Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
{
    m_fieldsALE = fields;

    // Initialise grid velocity as 0s
    int spaceDim = pGraph->GetSpaceDimension();
    m_gridVelocity = Array<OneD, Array<OneD, NekDouble>>(spaceDim);
    for (int i = 0; i < spaceDim; ++i)
    {
        m_gridVelocity[i] = Array<OneD, NekDouble>(fields[0]->GetTotPoints(), 0.0);
    }

    // Create map from element ID to expansion ID
    auto exp = fields[0]->GetExp();
    for (int i = (*exp).size() - 1; i >= 0; --i)
    {
        m_elmtToExpId[(*exp)[i]->GetGeom()->GetGlobalID()] = i;
    }

    // Create ALE objects for each interface zone
    for (auto &zone : fields[0]->GetMovement()->GetZones())
    {
        switch(zone.second->GetMovementType())
        {
            case SpatialDomains::eFixed :
                m_ALEs.emplace_back(ALEFixedShPtr(MemoryManager<ALEFixed>::AllocateSharedPtr(zone.second)));
                break;
            case SpatialDomains::eTranslate :
                m_ALEs.emplace_back(ALETranslateShPtr(MemoryManager<ALETranslate>::AllocateSharedPtr(zone.second)));
                m_ALESolver = true;
                break;
            case SpatialDomains::eRotate :
                m_ALEs.emplace_back(ALERotateShPtr(MemoryManager<ALERotate>::AllocateSharedPtr(zone.second)));
                m_ALESolver = true;
                break;
            case SpatialDomains::ePrescribe :
                m_ALEs.emplace_back(ALEPrescribeShPtr(MemoryManager<ALEPrescribe>::AllocateSharedPtr(zone.second)));
                m_ALESolver = true;
                break;
        }
    }
}

void ALEHelper::UpdateGridVelocity(const NekDouble &time)
{
    // Reset grid velocity to 0
    int spaceDim = m_fieldsALE[0]->GetGraph()->GetSpaceDimension();
    for (int i = 0; i < spaceDim; ++i)
    {
        std::fill(m_gridVelocity[i].begin(), m_gridVelocity[i].end(), 0.0);
    }

    // Now update for each movement zone, adding the grid velocities
    for (auto &ALE : m_ALEs)
    {
        ALE->UpdateGridVel(time, m_fieldsALE, m_elmtToExpId, m_gridVelocity);
    }
}

void ALEHelper::ALEPreMultiplyMass(Array<OneD, Array<OneD, NekDouble> > &fields)
{
    const int nm = m_fieldsALE[0]->GetNcoeffs();
    MultiRegions::GlobalMatrixKey mkey(StdRegions::eMass);

    // Premultiply each field by the mass matrix
    for (int i = 0; i < m_fieldsALE.size(); ++i)
    {
        fields[i] = Array<OneD, NekDouble>(nm);
        m_fieldsALE[i]->GeneralMatrixOp_IterPerExp(
            mkey, m_fieldsALE[i]->GetCoeffs(), fields[i]);
    }
}

void ALEHelper::ALEDoElmtInvMass(Array<OneD, Array<OneD, NekDouble> > &traceNormals, Array<OneD, Array<OneD, NekDouble> > &fields, NekDouble time)
{
    m_fieldsALE[0]->Reset();
    m_fieldsALE[0]->GetTrace()->GetNormals(traceNormals);
    UpdateGridVelocity(time);

    // Update m_fields with u^n by multiplying by inverse mass
    // matrix. That's then used in e.g. checkpoint output and L^2 error
    // calculation.

    for (int i = 0; i < m_fieldsALE.size(); ++i)
    {
        m_fieldsALE[i]->MultiplyByElmtInvMass(
            fields[i], m_fieldsALE[i]->UpdateCoeffs());
        m_fieldsALE[i]->BwdTrans(
            m_fieldsALE[i]->GetCoeffs(), m_fieldsALE[i]->UpdatePhys());
    }
}

void ALEHelper::ALEDoOdeProjection(const NekDouble &time, Array<OneD, Array<OneD, NekDouble> > &traceNormals)
{
    if (time != m_prevStageTime)
    {
        for (auto & i : m_fieldsALE)
        {
            i->GetMovement()->PerformMovement(time);
            i->Reset();
        }

        m_fieldsALE[0]->SetUpPhysNormals();
        m_fieldsALE[0]->GetTrace()->GetNormals(traceNormals);

        // Recompute grid velocity.
        UpdateGridVelocity(time);

        m_prevStageTime = time;
    }
}

void ALEHelper::ALEDoOdeRhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                            Array<OneD,        Array<OneD, NekDouble> >&outarray,
                            const NekDouble time,
                            AdvectionSharedPtr advObject,
                            Array<OneD, Array<OneD, NekDouble> > &velocity)
{
    const int nc = m_fieldsALE[0]->GetNcoeffs();
    int nVariables = inarray.size();

    // General idea is that we are time-integrating the quantity (Mu), so we
    // need to multiply input by inverse mass matrix to get coefficients u,
    // and then backwards transform so we can apply the DG operator.
    Array<OneD, NekDouble> tmp(nc);
    Array<OneD, Array<OneD, NekDouble>> tmpin(nVariables);

    for (int i = 0; i < nVariables; ++i)
    {
        tmpin[i] = Array<OneD, NekDouble>(m_fieldsALE[0]->GetNpoints());
        m_fieldsALE[i]->MultiplyByElmtInvMass(inarray[i], tmp);
        m_fieldsALE[i]->BwdTrans(tmp, tmpin[i]);
    }

    // RHS computation using the new advection base class
    advObject->AdvectCoeffs(nVariables, m_fieldsALE, velocity, tmpin,
                                 outarray, time);
}

ALEFixed::ALEFixed(SpatialDomains::ZoneBaseShPtr zone) :
    m_zone(std::static_pointer_cast<SpatialDomains::ZoneFixed>(zone))
{
}

void ALEFixed::v_UpdateGridVel(NekDouble time,
                                   Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                   std::map<int, int> &elmtToExpId,
                                   Array<OneD, Array<OneD, NekDouble>> &gridVelocity)
{
    boost::ignore_unused(time, fields, elmtToExpId, gridVelocity);
}

ALETranslate::ALETranslate(SpatialDomains::ZoneBaseShPtr zone) :
      m_zone(std::static_pointer_cast<SpatialDomains::ZoneTranslate>(zone))
{
}

void ALETranslate::v_UpdateGridVel(NekDouble time,
                                   Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                   std::map<int, int> &elmtToExpId,
                                   Array<OneD, Array<OneD, NekDouble>> &gridVelocity)
{
    boost::ignore_unused(time);

    auto vel = m_zone->GetVel();
    auto exp = fields[0]->GetExp();

    auto ids = m_zone->GetElementIds();
    for (auto id : ids)
    {
        int indx       = elmtToExpId[id];
        int offset     = fields[0]->GetPhys_Offset(indx);
        auto expansion = (*exp)[indx];

        int nq = expansion->GetTotPoints();
        for (int i = 0; i < nq; ++i)
        {
            gridVelocity[0][offset + i] += vel[0];
            gridVelocity[1][offset + i] += vel[1];
        }
    }
}

ALERotate::ALERotate(SpatialDomains::ZoneBaseShPtr zone) :
    m_zone(std::static_pointer_cast<SpatialDomains::ZoneRotate>(zone))
{
}

void ALERotate::v_UpdateGridVel(NekDouble time,
                               Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                               std::map<int, int> &elmtToExpId,
                               Array<OneD, Array<OneD, NekDouble>> &gridVelocity)
{
    boost::ignore_unused(time);

    auto angVel = m_zone->GetAngularVel();
    auto exp = fields[0]->GetExp();

    auto ids = m_zone->GetElementIds();
    for (auto id : ids)
    {
        int indx       = elmtToExpId[id];
        int offset     = fields[0]->GetPhys_Offset(indx);
        auto expansion = (*exp)[indx];

        int nq = expansion->GetTotPoints();
        Array<OneD, NekDouble> xc(nq, 0.0), yc(nq, 0.0), zc(nq, 0.0);
        expansion->GetCoords(xc, yc, zc);

        for (int i = 0; i < nq; ++i)
        {
            gridVelocity[0][offset + i] += -angVel*yc[i]; // @TODO: Update to include axis and origin.
            gridVelocity[1][offset + i] += angVel*xc[i];
        }
    }
}

ALEPrescribe::ALEPrescribe(SpatialDomains::ZoneBaseShPtr zone) :
    m_zone(std::static_pointer_cast<SpatialDomains::ZonePrescribe>(zone))
{
}

void ALEPrescribe::v_UpdateGridVel(NekDouble time,
                                Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                std::map<int, int> &elmtToExpId,
                                Array<OneD, Array<OneD, NekDouble>> &gridVelocity)
{
    boost::ignore_unused(time, fields, elmtToExpId, gridVelocity);
}

}
}