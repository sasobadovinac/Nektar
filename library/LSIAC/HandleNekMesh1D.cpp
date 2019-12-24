///////////////////////////////////////////////////////////////////////////////
//
// File HandleNekMesh1D.cpp
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
// Description: HandleNekMesh1D definition
//
///////////////////////////////////////////////////////////////////////////////
#include "HandleNekMesh1D.h"

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/DisContField1D.h>
#include <cmath>

namespace Nektar
{
namespace LSIAC
{

//! This function given range of tmin and tmax returns the element break points.
/*
        \param xcen_offset,ycen_offset,zcen_offset,direction,t_offsetmin
   t_offset_max \param [out] xPos,yPos,zPos,tPos This function makes few
   assumptions for simplicity.
        -> This function does not gaurentee if all of breakpoints tmin and tmax
   are returned.
        -> This function always includes tmin and tmax as break points while
   returning.
        -> There is no significance to bool in this function. It always return
   true.
        -> This function assumes all the seg element are non-curve elements.
*/
bool HandleNekMesh1D::v_LoadMesh(string var)
{
    SpatialDomains::ExpansionMap expansions = m_graph->GetExpansions();
    m_expansions.push_back(
        MemoryManager<MultiRegions::DisContField1D>::AllocateSharedPtr(
            m_session, m_graph, var));
    return true;
}

bool HandleNekMesh1D::v_LoadData(string filename, vector<string> &variables)
{
    SpatialDomains::ExpansionMap expansions = m_graph->GetExpansions();
    for (int i = 0; i < variables.size(); i++)
    {
        m_expansions.push_back(
            MemoryManager<MultiRegions::DisContField1D>::AllocateSharedPtr(
                m_session, m_graph, variables[i]));
    }

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> rFieldDef;
    std::vector<std::vector<NekDouble>> rFieldData;

    Array<OneD, int> ElementGIDs(expansions.size());
    SpatialDomains::ExpansionMap::const_iterator expIt;
    int i = 0;
    for (expIt = expansions.begin(); expIt != expansions.end(); ++expIt)
    {
        ElementGIDs[i++] = expIt->second->m_geomShPtr->GetGlobalID();
    }

    m_fld->Import(filename, rFieldDef, rFieldData,
                  LibUtilities::NullFieldMetaDataMap, ElementGIDs);
    for (int i = 0; i < rFieldDef.size(); i++)
    {
        for (int e = 0; e < variables.size(); e++)
        {
            m_expansions[e]->ExtractDataToCoeffs(
                rFieldDef[i], rFieldData[i], variables[e],
                m_expansions[e]->UpdateCoeffs());
        }
    }
    for (auto exp : m_expansions)
    {
        exp->BwdTrans(exp->GetCoeffs(), exp->UpdatePhys());
    }
    return true;
}

bool HandleNekMesh1D::v_Get1DVec(vector<NekDouble> &coords)
{
    coords.clear();
    SpatialDomains::PointGeomMap ptGM = m_graph->GetAllPointGeoms();
    vector<NekDouble> x_coords;
    NekDouble x, y, z;
    for (std::map<int, SpatialDomains::PointGeomSharedPtr>::iterator ptGM_it =
             ptGM.begin();
         ptGM_it != ptGM.end(); ++ptGM_it)
    {
        ptGM_it->second->GetCoords(x, y, z);
        coords.push_back(x);
    }
    std::sort(x_coords.begin(), x_coords.end());
    return true;
}

bool HandleNekMesh1D::v_GetKnotVec(const int degree,
                                   const Array<OneD, NekDouble> &coord,
                                   const Array<OneD, NekDouble> &direction,
                                   Array<OneD, NekDouble> &knotVec,
                                   NekDouble &shift)
{
    boost::ignore_unused(direction); // since it is in 1D
    // Assume all dimension in x-axis;
    SpatialDomains::PointGeomMap ptGM = m_graph->GetAllPointGeoms();
    vector<NekDouble> x_coords;
    NekDouble x, y, z;
    for (std::map<int, SpatialDomains::PointGeomSharedPtr>::iterator ptGM_it =
             ptGM.begin();
         ptGM_it != ptGM.end(); ++ptGM_it)
    {
        ptGM_it->second->GetCoords(x, y, z);
        x_coords.push_back(x);
    }
    std::sort(x_coords.begin(), x_coords.end());
    auto low = std::lower_bound(x_coords.begin(), x_coords.end(), coord[0]);
    auto up  = std::upper_bound(x_coords.begin(), x_coords.end(), coord[0]);
    if (low == up)
    {
        low--;
    }
    int distanceFromStart = std::distance(x_coords.begin(), low);
    int distanceToEnd     = std::distance(low, x_coords.end());

    int cknotnum, fknotnum;
    if (degree % 2 == 1)
    {
        cknotnum = ceil((3.0 * degree + 2.0) / 2.0);
        fknotnum = floor((3.0 * degree + 2.0) / 2.0);
    }
    else
    {
        cknotnum = ((3.0 * degree + 2.0) / 2.0);
        fknotnum = ((3.0 * degree + 2.0) / 2.0);
    }

    if (distanceFromStart < cknotnum - 1 || distanceToEnd < fknotnum + 2)
    {
        return false;
    }

    auto startIt = low;
    startIt      = low - cknotnum + 1;

    if (degree % 2 == 1)
    {
        NekDouble diff = *low;
        for (int i = 0; i < 3 * degree + 2; i++)
        {
            knotVec[i] = *(startIt + i) - diff;
        }
        shift = coord[0] - diff;
    }
    else
    {
        NekDouble diff = (*low + *(low + 1)) / 2.0;
        for (int i = 0; i < 3 * degree + 2; i++)
        {
            knotVec[i] = *(startIt + i) - diff;
        }
        shift = coord[0] - diff;
    }

    return true;
}

bool HandleNekMesh1D::v_GetBreakPts_Without_Tmin_Tmax(
    const NekDouble xcen_offset, const NekDouble ycen_offset,
    const NekDouble zcen_offset, const Array<OneD, NekDouble> &direction,
    const NekDouble tmin, const NekDouble tmax, vector<NekDouble> &xPos,
    vector<NekDouble> &yPos, vector<NekDouble> &zPos, vector<NekDouble> &tPos)
{
    // create a 1D geometry of the element.
    NekDouble x, y, z, xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = xcen_offset + tmin * direction[0] - TOLERENCE;
    xmax = xcen_offset + tmax * direction[0] + TOLERENCE;
    ymin = ycen_offset + tmin * direction[1] - TOLERENCE;
    ymax = ycen_offset + tmax * direction[1] + TOLERENCE;
    zmin = zcen_offset + tmin * direction[2] - TOLERENCE;
    zmax = zcen_offset + tmax * direction[2] + TOLERENCE;

    // Get all SegGeoms of the mesh.
    // loop through all segGeoms and match segment to segment overlap.
    xPos.clear();
    yPos.clear();
    zPos.clear();
    tPos.clear();
    SpatialDomains::PointGeomMap ptGM = m_graph->GetAllPointGeoms();
    for (std::map<int, SpatialDomains::PointGeomSharedPtr>::iterator ptGM_it =
             ptGM.begin();
         ptGM_it != ptGM.end(); ++ptGM_it)
    {
        ptGM_it->second->GetCoords(x, y, z);
        switch (m_spacedim)
        {
            case 1:
                if ((xmin < x && xmax > x))
                {
                    if ((xmin + 2 * TOLERENCE > x) ||
                        (xmax - 2 * TOLERENCE < x))
                    { // ignore it is on of the end points.
                    }
                    else
                    { // Add the point between the interval.
                        xPos.push_back(x);
                        tPos.push_back((x - xcen_offset) / direction[0]);
                    }
                }
                break;
            case 2:
                if ((xmin < x && xmax > x) && (ymin < y && ymax > y))
                {
                    if (((xmin + 2 * TOLERENCE > x) ||
                         (xmax - 2 * TOLERENCE < x)) &&
                        ((ymin + 2 * TOLERENCE > y) ||
                         (ymax - 2 * TOLERENCE < y)))
                    { // ignore it is on of the end points.
                    }
                    else
                    { // Add the point between the interval.
                        xPos.push_back(x);
                        yPos.push_back(y);
                        if (std::abs(direction[0]) > TOLERENCE)
                        {
                            tPos.push_back((x - xcen_offset) / direction[0]);
                        }
                        else
                        {
                            tPos.push_back((y - ycen_offset) / direction[0]);
                        }
                    }
                }
                break;
            case 3:
                if ((xmin < x && xmax > x) && (ymin < y && ymax > y) &&
                    (zmin < z && zmax > z))
                {
                    if (((xmin + 2 * TOLERENCE > x) ||
                         (xmax - 2 * TOLERENCE < x)) &&
                        ((ymin + 2 * TOLERENCE > y) ||
                         (ymax - 2 * TOLERENCE < y)) &&
                        ((zmin + 2 * TOLERENCE > z) ||
                         (zmax - 2 * TOLERENCE < z)))
                    { // ignore it is on of the end points.
                    }
                    else
                    { // Add the point between the interval.
                        xPos.push_back(x);
                        yPos.push_back(y);
                        zPos.push_back(z);
                        if (std::abs(direction[0]) > TOLERENCE)
                        {
                            tPos.push_back((x - xcen_offset) / direction[0]);
                        }
                        else if (std::abs(direction[1]) > TOLERENCE)
                        {
                            tPos.push_back((y - ycen_offset) / direction[1]);
                        }
                        else
                        {
                            tPos.push_back((z - zcen_offset) / direction[2]);
                        }
                    }
                }
                break;
            default:
                break;
        }
    }
    return true;
}

// This function does not verify if tmin to tmax is properly given.
bool HandleNekMesh1D::v_GetBreakPts(
    const NekDouble xcen_offset, const NekDouble ycen_offset,
    const NekDouble zcen_offset, const Array<OneD, NekDouble> &direction,
    const NekDouble tmin, const NekDouble tmax, vector<NekDouble> &xPos,
    vector<NekDouble> &yPos, vector<NekDouble> &zPos, vector<NekDouble> &tPos)
{
    // create a 1D geometry of the element.
    NekDouble x, y, z, xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = xcen_offset + tmin * direction[0] - TOLERENCE;
    xmax = xcen_offset + tmax * direction[0] + TOLERENCE;
    ymin = ycen_offset + tmin * direction[1] - TOLERENCE;
    ymax = ycen_offset + tmax * direction[1] + TOLERENCE;
    zmin = zcen_offset + tmin * direction[2] - TOLERENCE;
    zmax = zcen_offset + tmax * direction[2] + TOLERENCE;

    // Get all SegGeoms of the mesh.
    // loop through all segGeoms and match segment to segment overlap.
    xPos.clear();
    yPos.clear();
    zPos.clear();
    tPos.clear();
    xPos.push_back(xmin + TOLERENCE);
    yPos.push_back(ymin + TOLERENCE);
    zPos.push_back(zmin + TOLERENCE);
    tPos.push_back(tmin);
    SpatialDomains::PointGeomMap ptGM = m_graph->GetAllPointGeoms();
    for (std::map<int, SpatialDomains::PointGeomSharedPtr>::iterator ptGM_it =
             ptGM.begin();
         ptGM_it != ptGM.end(); ++ptGM_it)
    {
        ptGM_it->second->GetCoords(x, y, z);
        switch (m_spacedim)
        {
            case 1:
                if ((xmin < x && xmax > x))
                {
                    if ((xmin + 2 * TOLERENCE > x) ||
                        (xmax - 2 * TOLERENCE < x))
                    { // ignore it is on of the end points.
                    }
                    else
                    { // Add the point between the interval.
                        xPos.push_back(x);
                        tPos.push_back((x - xcen_offset) / direction[0]);
                    }
                }
                break;
            case 2:
                if ((xmin < x && xmax > x) && (ymin < y && ymax > y))
                {
                    if (((xmin + 2 * TOLERENCE > x) ||
                         (xmax - 2 * TOLERENCE < x)) &&
                        ((ymin + 2 * TOLERENCE > y) ||
                         (ymax - 2 * TOLERENCE < y)))
                    { // ignore it is on of the end points.
                    }
                    else
                    { // Add the point between the interval.
                        xPos.push_back(x);
                        yPos.push_back(y);
                        if (std::abs(direction[0]) > TOLERENCE)
                        {
                            tPos.push_back((x - xcen_offset) / direction[0]);
                        }
                        else
                        {
                            tPos.push_back((y - ycen_offset) / direction[0]);
                        }
                    }
                }
                break;
            case 3:
                if ((xmin < x && xmax > x) && (ymin < y && ymax > y) &&
                    (zmin < z && zmax > z))
                {
                    if (((xmin + 2 * TOLERENCE > x) ||
                         (xmax - 2 * TOLERENCE < x)) &&
                        ((ymin + 2 * TOLERENCE > y) ||
                         (ymax - 2 * TOLERENCE < y)) &&
                        ((zmin + 2 * TOLERENCE > z) ||
                         (zmax - 2 * TOLERENCE < z)))
                    { // ignore it is on of the end points.
                    }
                    else
                    { // Add the point between the interval.
                        xPos.push_back(x);
                        yPos.push_back(y);
                        zPos.push_back(z);
                        if (std::abs(direction[0]) > TOLERENCE)
                        {
                            tPos.push_back((x - xcen_offset) / direction[0]);
                        }
                        else if (std::abs(direction[1]) > TOLERENCE)
                        {
                            tPos.push_back((y - ycen_offset) / direction[1]);
                        }
                        else
                        {
                            tPos.push_back((z - zcen_offset) / direction[2]);
                        }
                    }
                }
                break;
            default:
                break;
        }
    }
    xPos.push_back(xmax - TOLERENCE);
    yPos.push_back(ymax - TOLERENCE);
    zPos.push_back(zmax - TOLERENCE);
    tPos.push_back(tmax);
    return true;
}

bool HandleNekMesh1D::v_CanTRangebeApplied(
    const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
    const NekDouble scaling, const NekDouble tmin, const NekDouble tmax,
    NekDouble &tminUpdate, NekDouble &tmaxUpdate)
{
    boost::ignore_unused(ptsY, ptsZ); // since it is in 1D
    int nq = m_expansions[0]->GetTotPoints();
    Array<OneD, NekDouble> xc0(nq), xc1(nq), xc2(nq);
    switch (m_expansions[0]->GetCoordim(0))
    {
        case 1:
            m_expansions[0]->GetCoords(xc0);
            break;
        case 2:
            m_expansions[0]->GetCoords(xc0, xc1);
            break;
        case 3:
            m_expansions[0]->GetCoords(xc0, xc1, xc2);
            break;
        default:
            break;
    }

    if ((xc0[0] <= ptsX + scaling * tmin) &&
        (xc0[nq - 1] >= ptsX + scaling * tmax))
    {
        return true;
    }
    else
    {
        if (xc0[0] > ptsX + scaling * tmin)
        {
            tminUpdate = xc0[0] / scaling;
        }
        if (xc0[nq - 1] < ptsX + scaling * tmax)
        {
            tmaxUpdate = xc0[nq - 1] / scaling;
        }
        return false;
    }
}

bool HandleNekMesh1D::v_CanTRangebeAppliedWOMeshShift(
    const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
    const Array<OneD, NekDouble> &direction, const NekDouble tmin,
    const NekDouble tmax)
{
    boost::ignore_unused(ptsX, ptsY, ptsZ, direction, tmin, tmax);
    NEKERROR(ErrorUtil::efatal, "Need to code this to use it");
    return false;
}

bool HandleNekMesh1D::v_CanTRangebeApplied(
    const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
    const Array<OneD, NekDouble> &direction, const NekDouble tmin,
    const NekDouble tmax, NekDouble &meshTShift)
{
    boost::ignore_unused(ptsY, ptsZ);
    int nq = m_expansions[0]->GetTotPoints();
    Array<OneD, NekDouble> xc0(nq), xc1(nq), xc2(nq);
    switch (m_expansions[0]->GetCoordim(0))
    {
        case 1:
            m_expansions[0]->GetCoords(xc0);
            break;
        case 2:
            m_expansions[0]->GetCoords(xc0, xc1);
            break;
        case 3:
            m_expansions[0]->GetCoords(xc0, xc1, xc2);
            break;
        default:
            break;
    }
    NekDouble xMeshMin = Vmath::Vmin(nq, xc0, 1);
    NekDouble xMeshMax = Vmath::Vmax(nq, xc0, 1);

    if ((std::abs(direction[0]) > TOLERENCE) &&
        (std::abs(direction[1]) < TOLERENCE) &&
        (std::abs(direction[2]) < TOLERENCE))
    {
        if ((xMeshMin <= ptsX + direction[0] * tmin) &&
            (xMeshMax >= ptsX + direction[0] * tmax))
        {
            return true;
        }
        else
        {
            if (xMeshMin > ptsX + direction[0] * tmin)
            { // This is mesh shift. kernel shift should be in the oppsite
              // direction.
                meshTShift =
                    (xMeshMin - (ptsX + direction[0] * tmin)) / direction[0];
            }
            if (xMeshMax < ptsX + direction[0] * tmax)
            {
                meshTShift =
                    (xMeshMax - (ptsX + direction[0] * tmax)) / direction[0];
            }
            if ((xMeshMin > ptsX + direction[0] * tmin) &&
                (xMeshMax < ptsX + direction[0] * tmax))
            {
                NEKERROR(ErrorUtil::efatal, "Kernel bigger than mesh");
            }
            return false;
        }
    }
    else
    {
        NEKERROR(ErrorUtil::efatal, "Only x-direction supported for 1D meshes");
        return false;
    }
}

bool HandleNekMesh1D::v_EvaluateAt(const NekDouble xPos, const NekDouble yPos,
                                   const NekDouble zPos, int gID, int eID,
                                   NekDouble &value, int varNum)
{
    boost::ignore_unused(xPos, yPos, zPos, gID, eID, value, varNum);
    NEKERROR(ErrorUtil::efatal, "Need to be implemented");
    return false;
}

bool HandleNekMesh1D::v_EvaluateAt(const Array<OneD, NekDouble> &xPos,
                                   const Array<OneD, NekDouble> &yPos,
                                   const Array<OneD, NekDouble> &zPos,
                                   const int gID, const int eID,
                                   Array<OneD, NekDouble> &values, int varNum)
{
    boost::ignore_unused(gID); // reserved for global id if implemented.
    if (!m_expansions[varNum]->GetPhysState())
    {
        m_expansions[varNum]->BwdTrans(m_expansions[varNum]->GetCoeffs(),
                                       m_expansions[varNum]->UpdatePhys());
    }

    // The reason for asking gID will be useful if we are using MPI.
    LocalRegions::ExpansionSharedPtr lexp = m_expansions[varNum]->GetExp(eID);
    const int phys_offset = m_expansions[varNum]->GetPhys_Offset(eID);
    const Array<OneD, NekDouble> g_Phys = m_expansions[varNum]->GetPhys();

    const Array<OneD, NekDouble> el_Phys =
        g_Phys.CreateWithOffset(g_Phys, phys_offset);
    Array<OneD, NekDouble> glCoord(3, 0.0), lCoord(3, 0.0);
    for (int i = 0; i < values.num_elements(); i++)
    {
        glCoord[0] = xPos[i];
        glCoord[1] = yPos[i];
        glCoord[2] = zPos[i];
        values[i]  = lexp->PhysEvaluate(glCoord, el_Phys);
    }
    return true;
}

bool HandleNekMesh1D::v_GetListOfGIDs(
    const NekDouble xPos, const NekDouble yPos, const NekDouble zPos,
    const Array<OneD, NekDouble> &direction, const vector<NekDouble> t_breaks,
    vector<int> &t_GIDs, vector<int> &t_EIDs) const
{
    t_GIDs.clear();
    t_EIDs.clear();
    t_GIDs.resize(t_breaks.size());
    t_EIDs.resize(t_breaks.size());
    SpatialDomains::ExpansionMap expansions = m_graph->GetExpansions();
    for (int i = 0; i < t_breaks.size() - 1; i++)
    {
        t_GIDs[i]         = -1;
        t_EIDs[i]         = -1;
        NekDouble t_break = (t_breaks[i] + t_breaks[i + 1]) / 2.0;
        Array<OneD, NekDouble> locCoord(3, 0.0);
        locCoord[0] = xPos + t_break * direction[0];
        locCoord[1] = yPos + t_break * direction[1];
        locCoord[2] = zPos + t_break * direction[2];
        t_GIDs[i]   = m_expansions[0]->GetExpIndex(locCoord);
        t_EIDs[i]   = t_GIDs[i];
    }

    return true;
}

NekDouble HandleNekMesh1D::v_GetLargestEdgeLength(const int eid)
{
    SpatialDomains::GeometrySharedPtr geomSPtr =
        m_expansions[0]->GetExp(eid)->GetGeom();
    // One Edge

    NekDouble edgeLength =
        geomSPtr->GetVertex(0)->dist(*(geomSPtr->GetVertex(1)));
    return edgeLength;
}

NekDouble HandleNekMesh1D::v_GetDynamicScaling(Array<OneD, NekDouble> glCoord,
                                               int eid, NekDouble mu)
{
    // if eid <0 find a elid.
    if (eid < 0)
    {
        eid = GetExpansionIndexUsingRTree(glCoord);
    }
    ASSERTL0(eid >= 0, "Point out of mesh");
    NekDouble result = -1.0;
    // Get local coordinates.
    // Depending on number of vertices triangle or quad.
    // use locCoordinates as barycentric coordinates
    SpatialDomains::GeometrySharedPtr geomSPtr =
        m_expansions[0]->GetExp(eid)->GetGeom();
    Array<OneD, NekDouble> lCoord(3, 0.0);
    geomSPtr->GetLocCoords(glCoord, lCoord);
    if (geomSPtr->GetShapeType() == Nektar::LibUtilities::eSegment)
    {
        int Vid0          = geomSPtr->GetVid(0);
        int Vid1          = geomSPtr->GetVid(1);
        NekDouble lambda1 = (lCoord[0] + 1.0) / 2.0;
        NekDouble lambda0 = 1.0 - lambda1;
        result =
            lambda0 * m_dynVertScaling[Vid0] + lambda1 * m_dynVertScaling[Vid1];
    }
    return mu * result;
}

bool HandleNekMesh1D::v_InitializeMetricTensor()
{
    m_metricTensor = new MetricTensor();
    m_metricTensor->LoadMetricTensor(this);
    m_MTDefined = true;
    return true;
}

bool HandleNekMesh1D::v_GetMTScalingOfGIDs(vector<int> &t_GIDs,
                                           Array<OneD, NekDouble> &direction,
                                           vector<NekDouble> &scalings)
{
    if (!m_MTDefined)
    {
        this->InitializeMetricTensor();
    }
    NekDouble lambda;
    for (auto gid : t_GIDs)
    {
        m_metricTensor->GetScaleForDirection(gid, direction, lambda);
        scalings.push_back(lambda);
        // scalings.push_back(0.1);
    }
    return true;
}
} // namespace LSIAC
} // namespace Nektar
