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
                break;
            case SpatialDomains::eRotate :
                m_ALEs.emplace_back(ALERotateShPtr(MemoryManager<ALERotate>::AllocateSharedPtr(zone.second)));
                break;
            case SpatialDomains::ePrescribe :
                m_ALEs.emplace_back(ALEPrescribeShPtr(MemoryManager<ALEPrescribe>::AllocateSharedPtr(zone.second)));
                break;
        }
    }
}

void ALEHelper::UpdateGridVelocity(NekDouble &time)
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
    boost::ignore_unused(time, fields, elmtToExpId, gridVelocity);
    /*SpatialDomains::ZoneRotateShPtr interfaceRotate =
                    std::static_pointer_cast<
                        SpatialDomains::ZoneRotate>(interface);
                NekDouble angVel = interfaceRotate->GetAngularVel();

                auto ids = interface->GetElementIds();
                for (auto id : ids)
                {
                    int indx       = elmtToExpId[id];
                    int offset     = m_fields[0]->GetPhys_Offset(indx);
                    auto expansion = (*exp)[indx];

                    int nq = expansion->GetTotPoints();
                    Array<OneD, NekDouble> xc(nq, 0.0), yc(nq, 0.0),
                        zc(nq, 0.0);
                    expansion->GetCoords(xc, yc, zc);

                    for (int i = 0; i < nq; ++i)
                    {
                        m_gridVelocity[0][offset + i] = -angVel*yc[i];
                        m_gridVelocity[1][offset + i] = angVel*xc[i];
                        //std::cout << "Coordinate: (" << xc[i] << ", " << yc[i] << ") has velocity = "
                        //          << m_gridVelocity[0][offset + i] << ", " << m_gridVelocity[1][offset + i] << std::endl;
                    }
                }
     */
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