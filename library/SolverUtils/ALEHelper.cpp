#include "ALEHelper.h"

namespace Nektar
{

namespace SolverUtils
{

void ALEHelper::InitObject(const SpatialDomains::MeshGraphSharedPtr &pGraph,
                     Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
{
    m_fieldsALE = fields;

    // Initialise grid velocities as 0s
    int spaceDim = pGraph->GetSpaceDimension();
    m_gridVelocity = Array<OneD, Array<OneD, NekDouble>>(spaceDim);
    m_gridVelocityTrace = Array<OneD, Array<OneD, NekDouble>>(spaceDim);
    for (int i = 0; i < spaceDim; ++i)
    {
        m_gridVelocity[i] = Array<OneD, NekDouble>(fields[0]->GetTotPoints(), 0.0);
        m_gridVelocityTrace[i] = Array<OneD, NekDouble>(fields[0]->GetTrace()->GetTotPoints(), 0.0);
    }

    // Create map from element ID to expansion ID
    auto exp = fields[0]->GetExp();
    for (int i = (*exp).size() - 1; i >= 0; --i)
    {
        m_elmtToExpId[(*exp)[i]->GetGeom()->GetGlobalID()] = i;
    }

    // Create map from element ID to expansion ID for the trace
    auto expTrace = fields[0]->GetTrace()->GetExp();
    for (int i = (*expTrace).size() - 1; i >= 0; --i)
    {
        m_elmtToExpIdTrace[(*expTrace)[i]->GetGeom()->GetGlobalID()] = i;
    }

    // Create ALE objects for each interface zone
    for (auto &zone : fields[0]->GetMovement()->GetZones())
    {
        switch(zone.second->GetMovementType())
        {
            case SpatialDomains::MovementType::eFixed :
                m_ALEs.emplace_back(ALEFixedShPtr(MemoryManager<ALEFixed>::AllocateSharedPtr(zone.second)));
                break;
            case SpatialDomains::MovementType::eTranslate :
                m_ALEs.emplace_back(ALETranslateShPtr(MemoryManager<ALETranslate>::AllocateSharedPtr(zone.second)));
                m_ALESolver = true;
                break;
            case SpatialDomains::MovementType::eRotate :
                m_ALEs.emplace_back(ALERotateShPtr(MemoryManager<ALERotate>::AllocateSharedPtr(zone.second)));
                m_ALESolver = true;
                break;
            case SpatialDomains::MovementType::ePrescribe :
                m_ALEs.emplace_back(ALEPrescribeShPtr(MemoryManager<ALEPrescribe>::AllocateSharedPtr(zone.second)));
                m_ALESolver = true;
                break;
            case SpatialDomains::MovementType::eNone :
                WARNINGL0(false, "Zone cannot have movement type of 'None'.")
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
    boost::ignore_unused(time, traceNormals);
    // Update m_fields with u^n by multiplying by inverse mass
    // matrix. That's then used in e.g. checkpoint output and L^2 error
    // calculation.

    // @TODO: Look at geometric factor and junk only what is needed
    // @TODO: Look at collections and see if they offer a speed up
    for (int i = 0; i < m_fieldsALE.size(); ++i)
    {
        m_fieldsALE[i]->MultiplyByElmtInvMass(
            fields[i], m_fieldsALE[i]->UpdateCoeffs()); // @TODO: Potentially matrix free?
        m_fieldsALE[i]->BwdTrans(
            m_fieldsALE[i]->GetCoeffs(), m_fieldsALE[i]->UpdatePhys());
    }
}

void ALEHelper::ALEDoElmtInvMassBwdTrans(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    const int nc   = m_fieldsALE[0]->GetNcoeffs();
    int nVariables = inarray.size();

    // General idea is that we are time-integrating the quantity (Mu), so we
    // need to multiply input by inverse mass matrix to get coefficients u,
    // and then backwards transform to physical space so we can apply the DG
    // operator.
    Array<OneD, NekDouble> tmp(nc);
    for (int i = 0; i < nVariables; ++i)
    {
        outarray[i] = Array<OneD, NekDouble>(m_fieldsALE[0]->GetNpoints());
        m_fieldsALE[i]->MultiplyByElmtInvMass(inarray[i], tmp);
        m_fieldsALE[i]->BwdTrans(tmp, outarray[i]);
    }
}

void ALEHelper::MoveMesh(const NekDouble &time, Array<OneD, Array<OneD, NekDouble> > &traceNormals)
{

    auto curvedEdges = m_fieldsALE[0]->GetGraph()->GetCurvedEdges();
    auto curvedFaces = m_fieldsALE[0]->GetGraph()->GetCurvedFaces();

    if (time != m_prevStageTime)
    {
        // The order of the resets below is v important to avoid errors
        for (auto &field : m_fieldsALE)
        {
            field->GetMovement()->PerformMovement(time);
            field->ResetMatrices();
        }

        // Loop over all elements and faces and edges and reset geometry information.
        // Only need to do this on the first field as the geometry information is shared.
        for (auto &zone : m_fieldsALE[0]->GetMovement()->GetZones())
        {
            if (zone.second->GetMoved())
            {
                auto conEl = zone.second->GetConstituentElements();
                for (const auto &i : conEl)
                {
                    for (const auto &j : i)
                    {
                        j->ResetNonRecursive(curvedEdges, curvedFaces);
                    }
                }

                // We need to rebuild geometric factors on the trace elements
                for (const auto &i : conEl[zone.second->GetConstituentElements().size() -1]) // This only takes the trace elements
                {
                    m_fieldsALE[0]->GetTrace()->GetExp(m_elmtToExpIdTrace[i->GetGlobalID()])->Reset();
                }
            }
        }

        // @TODO: I don't know why I need to reset the elements expansions for every field but the trace expansions on just one.
        // @TODO: I guess we need the metric info on element expansions, but only the geomfactors for the trace geometry.
        for (auto &field : m_fieldsALE)
        {
            for (auto &zone : field->GetMovement()->GetZones())
            {
                if (zone.second->GetMoved())
                {
                    auto conEl = zone.second->GetConstituentElements();
                    // Loop over zone elements expansions and rebuild geometric factors
                    for (const auto &i : conEl[0]) // This only takes highest dimensioned elements
                    {
                        field->GetExp(m_elmtToExpId[i->GetGlobalID()])->Reset();
                    }
                }
            }
        }

        for (auto &zone : m_fieldsALE[0]->GetMovement()->GetZones())
        {
            if (zone.second->GetMoved())
            {
                auto conEl = zone.second->GetConstituentElements();
                // Loop over zone elements expansions and rebuild geometric factors and recalc trace normals
                for (const auto &i : conEl[0]) // This only takes highest dimensioned elements
                {
                    for (int j = 0; j < m_fieldsALE[0]->GetExp(m_elmtToExpId[i->GetGlobalID()])->GetNtraces();++j)
                    {
                        m_fieldsALE[0]->GetExp(m_elmtToExpId[i->GetGlobalID()])->ComputeTraceNormal(j);
                    }
                }
            }
        }

        for (auto &field : m_fieldsALE)
        {
            // Reset collections
            field->CreateCollections(Collections::eNoImpType);
        }

        // Reload new trace normals in to the solver cache
        m_fieldsALE[0]->GetTrace()->GetNormals(traceNormals);

        // Recompute grid velocity.
        UpdateGridVelocity(time);

        // Updates trace grid velocity
        for (int i = 0; i < m_gridVelocityTrace.size(); ++i)
        {
            m_fieldsALE[0]->ExtractTracePhys(m_gridVelocity[i], m_gridVelocityTrace[i]);
        }

        m_prevStageTime = time;
    }
}

const Array<OneD, const Array<OneD, NekDouble> > &ALEHelper::GetGridVelocityTrace()
{
    return m_gridVelocityTrace;
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
            for (int j = 0; j < gridVelocity.size(); ++j)
            {
                gridVelocity[j][offset + i] += vel[j];
            }
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

    auto angVel = m_zone->GetAngularVel(time);
    auto axis = m_zone->GetAxis();
    auto origin = m_zone->GetOrigin();

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
            // Vector from origin to point
            DNekVec pointMinOrigin = {xc[i] - origin(0),
                                      yc[i] - origin(1),
                                      zc[i] - origin(2)};

            // Vector orthogonal to plane formed by axis and point
            DNekVec norm = pointMinOrigin.Cross(axis);
            //std::cout << xc[i] << " " << yc[i] << " " << zc[i] << " " << norm[0] << " " << norm[1] << " " << norm[2] << std::endl;

            // Distance between point and axis
            //DNekVec dist = pointMinOrigin - (pointMinOrigin.Dot(axis) * axis);

            // @TODO: Not sure here? multiply normal vector by distance from axis and angular velocity?
            //norm.Normalize();
            norm = norm * -angVel; // We negate here as by convention a positive angular velocity is counter-clockwise

            for (int j = 0; j < gridVelocity.size(); ++j)
            {
                gridVelocity[j][offset + i] = norm[j];
            }
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