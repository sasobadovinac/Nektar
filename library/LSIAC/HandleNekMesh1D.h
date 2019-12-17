#pragma once

#include "HandleNekMesh.h"
/// Handles Nektar 1D meshes.
/**
 */
namespace Nektar
{
namespace LSIAC
{

class HandleNekMesh1D : public HandleNekMesh
{
public:
    HandleNekMesh1D(FieldUtils::FieldSharedPtr fldSharedPtr)
        : HandleNekMesh(fldSharedPtr)
    {
    }
    HandleNekMesh1D(LibUtilities::SessionReaderSharedPtr sessionPtr)
        : HandleNekMesh(sessionPtr)
    {
        //			cout << "into HandleNekMesh1D constructor" <<
        // endl;
        //			m_expansions.push_back(MemoryManager<MultiRegions::ContField1D>::
        //					AllocateSharedPtr(m_session,m_graph,
        // m_session->GetVariable(0)) ); 			cout << "Expansion
        // Type: " << m_expansions[0]->GetExpType() << endl;
    }

protected:
    virtual bool v_GetKnotVec(const int degree,
                              const Array<OneD, NekDouble> &coord,
                              const Array<OneD, NekDouble> &direction,
                              Array<OneD, NekDouble> &knotVec,
                              NekDouble &shift);
    virtual bool v_Get1DVec(vector<NekDouble> &coords);

    virtual bool v_GetBreakPts(
        const NekDouble xcen_offset, const NekDouble ycen_offset,
        const NekDouble zcen_offset, const Array<OneD, NekDouble> &direction,
        const NekDouble t_offset_min, const NekDouble t_offset_max,
        vector<NekDouble> &xPos, vector<NekDouble> &yPos,
        vector<NekDouble> &zPos, vector<NekDouble> &tPos);

    virtual bool v_GetBreakPts_Without_Tmin_Tmax(
        const NekDouble xcen_offset, const NekDouble ycen_offset,
        const NekDouble zcen_offset, const Array<OneD, NekDouble> &direction,
        const NekDouble t_offset_min, const NekDouble t_offset_max,
        vector<NekDouble> &xPos, vector<NekDouble> &yPos,
        vector<NekDouble> &zPos, vector<NekDouble> &tPos);

    virtual bool v_CanTRangebeApplied(
        const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
        const NekDouble scaling, const NekDouble tmin, const NekDouble tmax,
        NekDouble &tminUpdate, NekDouble &tmaxUpdate);
    virtual bool v_CanTRangebeApplied(
        const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
        const Array<OneD, NekDouble> &direction, const NekDouble tmin,
        const NekDouble tmax, NekDouble &meshShift);
    virtual bool v_CanTRangebeAppliedWOMeshShift(
        const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
        const Array<OneD, NekDouble> &direction, const NekDouble tmin,
        const NekDouble tmax);
    virtual bool v_LoadData(string Filename, vector<string> &variables);
    virtual bool v_LoadMesh(string var);

    virtual bool v_EvaluateAt(const Array<OneD, NekDouble> &xPos,
                              const Array<OneD, NekDouble> &yPos,
                              const Array<OneD, NekDouble> &zPos, const int gID,
                              const int eID, Array<OneD, NekDouble> &values,
                              int varNum = 0);

    virtual bool v_EvaluateAt(const NekDouble xPos, const NekDouble yPos,
                              const NekDouble zPos, int gID, int eID,
                              NekDouble &value, int varNum = 0);

    virtual bool v_GetListOfGIDs(const NekDouble xPos, const NekDouble yPos,
                                 const NekDouble zPos,
                                 const Array<OneD, NekDouble> &direction,
                                 const vector<NekDouble> t_breaks,
                                 vector<int> &t_GIDs,
                                 vector<int> &t_EIDs) const;
    virtual void v_LoadExpListIntoRTree()
    {
        assert(false && "Not implemented, should not use it.");
    }

    // virtual NekDouble v_GetJacobian(const int eID);
    virtual NekDouble v_GetLargestEdgeLength(const int eid);

    virtual NekDouble v_GetDynamicScaling(Array<OneD, NekDouble> glCoord,
                                          int eid = -1, NekDouble mu = 1.0);

    virtual bool v_InitializeMetricTensor();
    virtual bool v_GetMTScalingOfGIDs(vector<int> &t_GIDs,
                                      Array<OneD, NekDouble> &direction,
                                      vector<NekDouble> &scalings);
};
} // namespace LSIAC
} // namespace Nektar
