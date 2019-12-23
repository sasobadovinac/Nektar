///////////////////////////////////////////////////////////////////////////////
//
// File HandleNekMesh2D.cpp
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
// Description: HandleNekMesh2D definition
//
///////////////////////////////////////////////////////////////////////////////
#include "HandleNekMesh2D.h"

#include <MultiRegions/ContField2D.h>
#include <MultiRegions/DisContField2D.h>
//#include <SpatialDomains/MeshGraph2D.h>
#include <cmath>
#include <iomanip> // std::setprecision

namespace Nektar
{
namespace LSIAC
{
NekDouble HandleNekMesh2D::v_GetElLargestEdgeSize(const NekDouble ptsx,
                                                  const NekDouble ptsy,
                                                  const NekDouble ptsz,
                                                  int Elid)
{
    boost::ignore_unused(ptsz); // used in case of 3D data.
    if (Elid < 0)
    {
        // Find Element Id.
        Array<OneD, NekDouble> glCord(2, 0.0);
        glCord[0] = ptsx;
        glCord[1] = ptsy;
        // glCord[2] = ptsz;
        if (m_useRTree)
        {
            Elid = GetExpansionIndexUsingRTree(glCord);
        }
        else
        {
            Elid = m_expansions[0]->GetExpIndex(glCord);
        }
        NEKERROR(ErrorUtil : efatal, "Something wrong. Point outside boundary");
    }

    SpatialDomains::GeometrySharedPtr gEl =
        (m_expansions[0]->GetExp(Elid))->GetGeom();

    int numEdges        = gEl->GetNumEdges();
    NekDouble maxLength = -1.0;
    for (int i = 0; i < numEdges; i++)
    {
        int Eid = gEl->GetEid(i);
        SpatialDomains::PointGeomSharedPtr p0 =
            m_segMap.find(Eid)->second->GetVertex(0);
        SpatialDomains::PointGeomSharedPtr p1 =
            m_segMap.find(Eid)->second->GetVertex(1);
        SpatialDomains::PointGeom p2(*p1);
        NekDouble len = p0->dist(p2);
        if (len > maxLength)
        {
            maxLength = len;
        }
    }

    ASSERTL0(maxLength > 0 && "max Length > 0 ");
    return maxLength;
}

//! This function given range of tmin and tmax returns the element break points.
/*
    \param xcen_offset,ycen_offset,zcen_offset,direction,t_offsetmin
   t_offset_max
    \param [out] xPos,yPos,zPos,tPos This function makes few
   assumptions for simplicity.
        -> This function does not gaurentee if all of breakpoints tmin and tmax
   are returned.
        -> This function always includes tmin and tmax as break points while
   returning.
        -> There is no significance to bool in this function. It always return
   true.
        -> This function assumes all the seg element are non-curve elements.
*/
bool HandleNekMesh2D::v_LoadMesh(string var)
{
    SpatialDomains::ExpansionMap expansions = m_graph->GetExpansions();
    m_expansions.push_back(
        MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(
            m_session, m_graph, var));
    return true;
}

bool HandleNekMesh2D::v_LoadData(string filename, vector<string> &variables)
{
    SpatialDomains::ExpansionMap expansions = m_graph->GetExpansions();
    for (int i = 0; i < variables.size(); i++)
    {
        m_expansions.push_back(
            MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(
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
        m_Arrays.push_back(exp->GetPhys());
    }
    return true;
}

// This function does not verify if tmin to tmax is properly given.
// This function in 2D does not fill xpos,ypos and zpos yet.
bool HandleNekMesh2D::v_GetBreakPts(
    const NekDouble xcen_offset, const NekDouble ycen_offset,
    const NekDouble zcen_offset, const Array<OneD, NekDouble> &direction,
    const NekDouble tmin, const NekDouble tmax, vector<NekDouble> &xPos,
    vector<NekDouble> &yPos, vector<NekDouble> &zPos, vector<NekDouble> &tPos)
{
    boost::ignore_unused(xPos, yPos, zPos); // Only returning tPos for now.
    Array<OneD, NekDouble> point(3);
    point[0] = xcen_offset;
    point[1] = ycen_offset;
    point[2] = zcen_offset;

    if (m_useRTree)
    {
        IntersectWithEdgesUsingRTree(m_segMap, m_pointMap, direction, point,
                                     tmin, tmax, tPos);
    }
    else
    {
        IntersectWithEdges(m_segMap, m_pointMap, direction, point, tmin, tmax,
                           tPos);
    }
    if (0 != tPos.size())
    {
        if (!compare2NekDoublesH(tmin, tPos.front()))
        {
            tPos.insert(tPos.begin(), tmin);
        }
        if (!compare2NekDoublesH(tmax, tPos.back()))
        {
            tPos.insert(tPos.end(), tmax);
        }
    }
    else
    {
        tPos.push_back(tmin);
        tPos.push_back(tmax);
    }

    return true;
}

bool HandleNekMesh2D::v_GetBreakPts_Without_Tmin_Tmax(
    const NekDouble xcen_offset, const NekDouble ycen_offset,
    const NekDouble zcen_offset, const Array<OneD, NekDouble> &direction,
    const NekDouble tmin, const NekDouble tmax, vector<NekDouble> &xPos,
    vector<NekDouble> &yPos, vector<NekDouble> &zPos, vector<NekDouble> &tPos)
{
    boost::ignore_unused(xPos, yPos, zPos); // Only returning tPos for now.
    Array<OneD, NekDouble> point(3);
    point[0] = xcen_offset;
    point[1] = ycen_offset;
    point[2] = zcen_offset;

    if (m_useRTree)
    {
        IntersectWithEdgesUsingRTree(m_segMap, m_pointMap, direction, point,
                                     tmin, tmax, tPos);
    }
    else
    {
        IntersectWithEdges(m_segMap, m_pointMap, direction, point, tmin, tmax,
                           tPos);
    }

    return true;
}

bool HandleNekMesh2D::v_CanTRangebeApplied(
    const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
    const NekDouble scaling, const NekDouble tmin, const NekDouble tmax,
    NekDouble &tminUpdate, NekDouble &tmaxUpdate)
{
    boost::ignore_unused(ptsX, ptsY, ptsZ, scaling, tmin, tmax, tminUpdate,
                         tmaxUpdate);
    NEKERROR(ErrorUtil : efatal, "Not coded.");
    return false;
}

bool HandleNekMesh2D::v_CanTRangebeAppliedWOMeshShift(
    const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
    const Array<OneD, NekDouble> &direction, const NekDouble tmin,
    const NekDouble tmax)
{
    boost::ignore_unused(ptsZ); // used in case of 3D data.
    // pointer to expansion list.
    // Claculate left point.
    Array<OneD, NekDouble> pl(2), pr(2);
    vector<NekDouble> tPos;
    pl[0] = ptsX + direction[0] * tmin;
    pl[1] = ptsY + direction[1] * tmin;
    pr[0] = ptsX + direction[0] * tmax;
    pr[1] = ptsY + direction[1] * tmax;

    int pl_index, pr_index;
    if (m_useRTree)
    {
        pl_index = GetExpansionIndexUsingRTree(pl);
        pr_index = GetExpansionIndexUsingRTree(pr);
    }
    else
    {
        pl_index = m_expansions[0]->GetExpIndex(pl, TOLERENCE);
        pr_index = m_expansions[0]->GetExpIndex(pr, TOLERENCE);
    }

    if ((0 > pl_index) || (0 > pr_index))
    { // goes out of boundary
        return false;
    }
    return true;
}

bool HandleNekMesh2D::v_WhatIsTRange(const NekDouble PtsX, const NekDouble PtsY,
                                     const NekDouble PtsZ,
                                     const Array<OneD, NekDouble> &direction,
                                     NekDouble &tmin, NekDouble &tmax, int &num)
{
    Array<OneD, NekDouble> pl(2), pr(2);
    vector<NekDouble> tPos;
    pl[0] = PtsX + direction[0] * tmin;
    pl[1] = PtsY + direction[1] * tmin;
    pr[0] = PtsX + direction[0] * tmax;
    pr[1] = PtsY + direction[1] * tmax;

    int pl_index, pr_index;
    if (m_useRTree)
    {
        pl_index = GetExpansionIndexUsingRTree(pl);
        pr_index = GetExpansionIndexUsingRTree(pr);
    }
    else
    {
        pl_index = m_expansions[0]->GetExpIndex(pl, TOLERENCE);
        pr_index = m_expansions[0]->GetExpIndex(pr, TOLERENCE);
    }

    Array<OneD, NekDouble> point(3);
    point[0] = PtsX;
    point[1] = PtsY;
    point[2] = PtsZ;
    if (m_useRTree)
    {
        IntersectWithEdgesUsingRTree(m_segMap, m_pointMap, direction, point,
                                     tmin, tmax, tPos);
    }
    else
    {
        IntersectWithEdges(m_segMap, m_pointMap, direction, point, tmin, tmax,
                           tPos);
    }
    num = tPos.size();
    if (0 > pl_index || 0 > pr_index)
    {
        tmin = tPos[0];
        tmax = tPos[tPos.size() - 1];
        return false;
    }
    else
    {
        return true;
    }
}

bool HandleNekMesh2D::v_CanTRangebeApplied(
    const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
    const Array<OneD, NekDouble> &direction, const NekDouble tmin,
    const NekDouble tmax, NekDouble &meshTShift)
{
    // pointer to expansion list.
    // Claculate left point.
    Array<OneD, NekDouble> pl(2), pr(2);
    vector<NekDouble> tPos;
    pl[0] = PtsX + direction[0] * tmin;
    pl[1] = PtsY + direction[1] * tmin;
    pr[0] = PtsX + direction[0] * tmax;
    pr[1] = PtsY + direction[1] * tmax;
    int pl_index, pr_index;
    if (m_useRTree)
    {
        pl_index = GetExpansionIndexUsingRTree(pl);
        pr_index = GetExpansionIndexUsingRTree(pr);
    }
    else
    {
        pl_index = m_expansions[0]->GetExpIndex(pl, TOLERENCE);
        pr_index = m_expansions[0]->GetExpIndex(pr, TOLERENCE);
    }

    if ((0 > pl_index) || (0 > pr_index))
    { // goes out of boundary
        Array<OneD, NekDouble> point(3);
        point[0] = PtsX;
        point[1] = PtsY;
        point[2] = PtsZ;
        if (m_useRTree)
        {
            IntersectWithEdgesUsingRTree(m_segMap, m_pointMap, direction, point,
                                         tmin, tmax, tPos);
        }
        else
        {
            IntersectWithEdges(m_segMap, m_pointMap, direction, point, tmin,
                               tmax, tPos);
        }
        if ((0 > pl_index))
        {
            if (tPos.size() == 0)
            {
                return false;
            }
            meshTShift = tPos.front() - tmin;
        }
        if ((0 > pr_index))
        {
            meshTShift = tPos.back() - tmax;
        }
        if ((0 > pl_index) && (0 > pr_index))
        {
        }
        return false;
    }
    return true;
}

bool HandleNekMesh2D::v_EvaluateAt(const NekDouble xPos, const NekDouble yPos,
                                   const NekDouble zPos, int gID, int eID,
                                   NekDouble &value, int varNum)
{
    boost::ignore_unused(gID, zPos); // reserved for global id if implemented.
    Array<OneD, NekDouble> lcoord(2, 0.0);
    lcoord[0] = xPos;
    lcoord[1] = yPos;
    if (eID < 0)
    {
        if (m_useRTree)
        {
            eID = GetExpansionIndexUsingRTree(lcoord);
        }
        else
        {
            eID = m_expansions[0]->GetExpIndex(lcoord, TOLERENCE);
        }
        ASSERTL0(eID != -1 && "Input point is out of Mesh");
    }
    LocalRegions::ExpansionSharedPtr lexp = m_expansions[0]->GetExp(eID);
    const int phys_offset = m_expansions[0]->GetPhys_Offset(eID);
    const Array<OneD, NekDouble> el_Phys =
        m_Arrays[varNum].CreateWithOffset(m_Arrays[varNum], phys_offset);
    value = lexp->PhysEvaluate(lcoord, el_Phys);

    return true;
}

bool HandleNekMesh2D::v_EvaluateAt(const Array<OneD, NekDouble> &xPos,
                                   const Array<OneD, NekDouble> &yPos,
                                   const Array<OneD, NekDouble> &zPos,
                                   const int gID, const int eID,
                                   Array<OneD, NekDouble> &values, int varNum)
{
    boost::ignore_unused(gID); // reserved for global id if implemented.
    ASSERTL0(gID >= 0 && eID >= 0 && "Input paramerters are out of scope;");
    // The reason for asking gID will be useful if we are using MPI.
    LocalRegions::ExpansionSharedPtr lexp = m_expansions[0]->GetExp(eID);
    const int phys_offset = m_expansions[0]->GetPhys_Offset(eID);

    const Array<OneD, NekDouble> el_Phys =
        m_Arrays[varNum].CreateWithOffset(m_Arrays[varNum], phys_offset);
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

bool HandleNekMesh2D::v_GetListOfGIDs(
    const NekDouble xPos, const NekDouble yPos, const NekDouble zPos,
    const Array<OneD, NekDouble> &direction, const vector<NekDouble> t_breaks,
    vector<int> &t_GIDs, vector<int> &t_EIDs) const
{
    t_GIDs.clear();
    t_EIDs.clear();
    t_GIDs.resize(t_breaks.size());
    t_EIDs.resize(t_breaks.size());
    for (int i = 0; i < t_breaks.size() - 1; i++)
    {
        t_GIDs[i]         = -1;
        t_EIDs[i]         = -1;
        NekDouble t_break = (t_breaks[i] + t_breaks[i + 1]) / 2.0;
        Array<OneD, NekDouble> locCoord(2, 0.0);
        locCoord[0] = xPos + t_break * direction[0];
        locCoord[1] = yPos + t_break * direction[1];
        if (m_useRTree)
        {
            t_GIDs[i] = GetExpansionIndexUsingRTree(locCoord);
        }
        else
        {
            t_GIDs[i] = m_expansions[0]->GetExpIndex(locCoord, TOLERENCE);
        }
        t_EIDs[i] = t_GIDs[i];
        ASSERTL0(t_GIDs[i] >= 0 && "Cannot find element Id.");
    }
    return true;
}

// All private functions doing the bulk of calcualtions.
vector<NekDouble> HandleNekMesh2D::cross_Math(const vector<NekDouble> &r,
                                              const vector<NekDouble> &s)
{
    vector<NekDouble> ans(3);
    ans[0] = r[1] * s[2] - r[2] * s[1];
    ans[1] = r[2] * s[0] - r[0] * s[2];
    ans[2] = r[0] * s[1] - r[1] * s[0];
    return ans;
}

vector<NekDouble> HandleNekMesh2D::sub_Math(vector<NekDouble> &p2,
                                            vector<NekDouble> &p1)
{
    vector<NekDouble> r(3);
    r[0] = p2[0] - p1[0];
    r[1] = p2[1] - p1[1];
    r[2] = p2[2] - p1[2];
    return r;
}
NekDouble HandleNekMesh2D::dot_Math(vector<NekDouble> &p, vector<NekDouble> &q)
{
    return p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
}

NekDouble HandleNekMesh2D::norm2_Math(vector<NekDouble> p)
{
    return p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
}

bool HandleNekMesh2D::intersect(vector<NekDouble> &p1, vector<NekDouble> &p2,
                                vector<NekDouble> &q1, vector<NekDouble> &q2,
                                vector<NekDouble> &i1, vector<NekDouble> &i2)
{
    // Assuming all points are 3D.
    // 1. compute r and s ; r= P2-P1; s = Q2-Q1;
    i1.clear();
    i2.clear();
    vector<NekDouble> r, s, rCs, pMq, pMq_Cr, qMp, q2Mp, qMp_Cs, qMp_Cr;
    NekDouble t0, t1;
    r      = sub_Math(p2, p1);
    s      = sub_Math(q2, q1);
    rCs    = cross_Math(r, s);
    pMq    = sub_Math(p1, q1);
    pMq_Cr = cross_Math(pMq, r);
    if ((norm2_Math(rCs) < TOLERENCE) && (norm2_Math(pMq_Cr) < TOLERENCE))
    {
        // line segements are linear and coinside.
        // find i1 and i2;
        t0   = -1.0 * dot_Math(pMq, r) / dot_Math(r, r);
        q2Mp = sub_Math(q2, p1);
        t1   = dot_Math(q2Mp, r) / dot_Math(r, r);
        if (t0 > t1)
        {
            NekDouble temp = t0;
            t0             = t1;
            t1             = temp;
        }
        if (t0 < 0.0 && t1 < 0.0)
        { // ignore
            return false;
        }
        else if (t0 <= 0.0 && t1 <= 1.0)
        { // p1 and t1
            i1.push_back(p1[0]);
            i1.push_back(p1[1]);
            i1.push_back(p1[2]);
            i2.push_back(p1[0] + t1 * r[0]);
            i2.push_back(p1[1] + t1 * r[1]);
            i2.push_back(p1[2] + t1 * r[2]);
            return true;
        }
        else if (t0 <= 0.0 && t1 >= 1.0)
        {
            i1.push_back(p1[0]);
            i1.push_back(p1[1]);
            i1.push_back(p1[2]);
            i2.push_back(p2[0]);
            i2.push_back(p2[1]);
            return true;
        }
        else if (t0 <= 1.0 && t1 >= 1.0)
        { // t0 and P2
            i1.push_back(p1[0] + t0 * r[0]);
            i1.push_back(p1[1] + t0 * r[1]);
            i1.push_back(p1[2] + t0 * r[2]);
            i2.push_back(p2[0]);
            i2.push_back(p2[1]);
            i2.push_back(p2[2]);
            return true;
        }
        else // one case left t0>1 and t1>1
        {    // ignore
            return false;
        }
    }
    if ((norm2_Math(rCs) < TOLERENCE) && (norm2_Math(pMq_Cr) > TOLERENCE))
    {
        return false;
    }
    if (norm2_Math(rCs) > TOLERENCE)
    {
        qMp            = sub_Math(q1, p1);
        qMp_Cs         = cross_Math(qMp, s);
        qMp_Cr         = cross_Math(qMp, r);
        rCs            = cross_Math(r, s);
        NekDouble t    = std::sqrt(norm2_Math(qMp_Cs) / norm2_Math(rCs));
        NekDouble terr = (std::abs(qMp_Cs[0] - t * rCs[0]) +
                          std::abs(qMp_Cs[1] - t * rCs[1]) +
                          std::abs(qMp_Cs[2] - t * rCs[2]));
        NekDouble u    = std::sqrt(norm2_Math(qMp_Cr) / norm2_Math(rCs));
        NekDouble uerr = (std::abs(qMp_Cr[0] - u * rCs[0]) +
                          std::abs(qMp_Cr[1] - u * rCs[1]) +
                          std::abs(qMp_Cr[2] - u * rCs[2]));
        if ((t >= 0 && t <= 1) && (u >= 0 && u <= 1) && (terr < TOLERENCE) &&
            (uerr < TOLERENCE))
        {
            i1.push_back(p1[0] + t * r[0]);
            i1.push_back(p1[1] + t * r[1]);
            i1.push_back(p1[2] + t * r[2]);
            return true;
        }
        else
        {
            return false;
        }
    }
    return true;
}

void HandleNekMesh2D::IntersectWithEdges(
    const SpatialDomains::SegGeomMap &segMap,
    const SpatialDomains::PointGeomMap &pointMap,
    const Array<OneD, NekDouble> &dir, const Array<OneD, NekDouble> &point,
    const NekDouble t1, const NekDouble t2, vector<NekDouble> &tvalT)
{
    tvalT.clear();
    vector<NekDouble> p1(3), p2(3), i1(3), i2(3);
    int dirID = -1;
    p1[0]     = point[0] + t1 * dir[0];
    p1[1]     = point[1] + t1 * dir[1];
    p1[2]     = point[2] + t1 * dir[2];
    p2[0]     = point[0] + t2 * dir[0];
    p2[1]     = point[1] + t2 * dir[1];
    p2[2]     = point[2] + t2 * dir[2];
    for (int i = 0; i < 3; i++)
    {
        if (std::abs(dir[i]) > TOLERENCE)
        {
            dirID = i;
            break;
        }
    }
    ASSERTL0(dirID >= 0 && "Direction is not right.");
    for (int s = 0; s < segMap.size(); s++)
    {
        SpatialDomains::SegGeomSharedPtr segPtr =
            segMap.find(s)->second; //->second;
        vector<NekDouble> q1(3), q2(3);
        pointMap.find(segPtr->GetVid(0))
            ->second->GetCoords(q1[0], q1[1], q1[2]);
        pointMap.find(segPtr->GetVid(1))
            ->second->GetCoords(q2[0], q2[1], q2[2]);
        bool b = intersect(p1, p2, q1, q2, i1, i2);
        if (b)
        {
            NekDouble t = 0;
            t           = (i1[dirID] - point[dirID]) / dir[dirID];
            tvalT.push_back(t);
            if (i2.size() > 0)
            {
                NekDouble t = 0;
                t           = (i2[dirID] - point[dirID]) / dir[dirID];
                tvalT.push_back(t);
            }
        }
    }

    // sort and keep them unique
    std::sort(tvalT.begin(), tvalT.end());
    std::vector<NekDouble>::iterator ittt;
    ittt = std::unique(tvalT.begin(), tvalT.end(), this->compare2NekDoublesH);
    tvalT.resize(std::distance(tvalT.begin(), ittt));
}

void HandleNekMesh2D::IntersectWithEdgesUsingRTree(
    const SpatialDomains::SegGeomMap &segMap,
    const SpatialDomains::PointGeomMap &pointMap,
    const Array<OneD, NekDouble> &dir, const Array<OneD, NekDouble> &point,
    const NekDouble t1, const NekDouble t2, vector<NekDouble> &tvalT)
{
    boost::ignore_unused(segMap); // Moved away from using geometry;
    vector<NekDouble> p1(3), p2(3), i1(3), i2(3);
    p1[0] = point[0] + t1 * dir[0];
    p1[1] = point[1] + t1 * dir[1];
    p1[2] = point[2] + t1 * dir[2];
    p2[0] = point[0] + t2 * dir[0];
    p2[1] = point[1] + t2 * dir[1];
    p2[2] = point[2] + t2 * dir[2];

    vector<int> elIds, glIds;
    IntersectWithBoxUsingRTree(std::min(p1[0], p2[0]), std::min(p1[1], p2[1]),
                               std::min(p1[2], p2[2]) - TOLERENCE,
                               std::max(p1[0], p2[0]), std::max(p1[1], p2[1]),
                               std::max(p1[2], p2[2]) + TOLERENCE, elIds,
                               glIds);
    vector<int> EdgeIds;

    BOOST_FOREACH (int eId, elIds)
    {
        SpatialDomains::GeometrySharedPtr geomSPtr =
            m_expansions[0]->GetExp(eId)->GetGeom();
        for (int l = 0; l < geomSPtr->GetNumEdges(); l++)
        {
            int edgeId = geomSPtr->GetEid(l);
            EdgeIds.push_back(edgeId);
        }
    }

    IntersectWithFewEdges(EdgeIds, pointMap, dir, point, t1, t2, tvalT);
}

void HandleNekMesh2D::IntersectWithFewEdges(
    const vector<int> EdgeIds, const SpatialDomains::PointGeomMap &pointMap,
    const Array<OneD, NekDouble> &dir, const Array<OneD, NekDouble> &point,
    const NekDouble t1, const NekDouble t2, vector<NekDouble> &tvalT)
{
    tvalT.clear();
    vector<NekDouble> p1(3), p2(3), i1(3), i2(3);
    int dirID = -1;
    p1[0]     = point[0] + t1 * dir[0];
    p1[1]     = point[1] + t1 * dir[1];
    p1[2]     = point[2] + t1 * dir[2];
    p2[0]     = point[0] + t2 * dir[0];
    p2[1]     = point[1] + t2 * dir[1];
    p2[2]     = point[2] + t2 * dir[2];
    for (int i = 0; i < 3; i++)
    {
        if (std::abs(dir[i]) > TOLERENCE)
        {
            dirID = i;
            break;
        }
    }
    ASSERTL0(dirID >= 0 && "Direction is not right something is up ");
    for (int s = 0; s < EdgeIds.size(); s++)
    {
        SpatialDomains::SegGeomSharedPtr segPtr =
            m_segMap.find(EdgeIds[s])->second;
        vector<NekDouble> q1(3), q2(3);
        pointMap.find(segPtr->GetVid(0))
            ->second->GetCoords(q1[0], q1[1], q1[2]);
        pointMap.find(segPtr->GetVid(1))
            ->second->GetCoords(q2[0], q2[1], q2[2]);
        bool b = intersect(p1, p2, q1, q2, i1, i2);
        if (b)
        {
            NekDouble t = 0;
            t           = (i1[dirID] - point[dirID]) / dir[dirID];
            tvalT.push_back(t);
            if (i2.size() > 0)
            {
                NekDouble t = 0;
                t           = (i2[dirID] - point[dirID]) / dir[dirID];
                tvalT.push_back(t);
            }
        }
    }

    // sort and keep them unique
    std::sort(tvalT.begin(), tvalT.end());
    std::vector<NekDouble>::iterator ittt;
    ittt = std::unique(tvalT.begin(), tvalT.end(), this->compare2NekDoublesH);
    tvalT.resize(std::distance(tvalT.begin(), ittt));
}

void HandleNekMesh2D::FindElementIDForLineSegs(
    const vector<NekDouble> &tvalT, const Array<OneD, NekDouble> &point,
    const Array<OneD, NekDouble> &dir,
    const SpatialDomains::MeshGraphSharedPtr mesh_graph, vector<int> &EIDs)
{
    boost::ignore_unused(tvalT, point, dir, mesh_graph,
                         EIDs); // moved away from using geometry
    NEKERROR(ErrorUtil::efatal, "Not Yet implemented");
}

// Mostly funtions using RTree from here.
void HandleNekMesh2D::v_LoadExpListIntoRTree()
{
    // check if expansions are already loaded.
    ASSERTL0(m_expansions.size() > 0 &&
             "should have loaded atleast one expansion");
    MultiRegions::ExpListSharedPtr expList = m_expansions[0];
    int expSize                            = expList->GetExpSize();
    for (int e = 0; e < expSize; ++e)
    {
        SpatialDomains::GeometrySharedPtr geom = expList->GetExp(e)->GetGeom();
        int gid                                = geom->GetGlobalID();
        RBox b;
        GetBoundingOfElement(geom, b);
        m_rtree.insert(std::make_tuple(b, e, gid));
    }
    m_useRTree = true;
}

void HandleNekMesh2D::GetBoundingOfElement(
    SpatialDomains::GeometrySharedPtr elGeom, RBox &b)
{
    // Loop through all point in the vertices.
    // create min and max in all coordinates.
    // load the RBox data structure
    // return
    NekDouble elX, elY, elZ;
    NekDouble minX, minY, minZ;
    minX = std::numeric_limits<NekDouble>::max();
    minY = minX;
    minZ = minX;
    NekDouble maxX, maxY, maxZ;
    maxX       = std::numeric_limits<NekDouble>::lowest();
    maxY       = maxX;
    maxZ       = maxX;
    int numPts = elGeom->GetNumVerts();
    for (int v = 0; v < numPts; v++)
    {
        elGeom->GetVertex(v)->GetCoords(elX, elY, elZ);
        minX = getMin(minX, elX);
        minY = getMin(minY, elY);
        minZ = getMin(minZ, elZ);
        maxX = getMax(maxX, elX);
        maxY = getMax(maxY, elY);
        maxZ = getMax(maxZ, elZ);
    }

    Boostg::set<Boostg::min_corner, 0>(b, minX - TOLERENCE_MESH_COMP);
    Boostg::set<Boostg::min_corner, 1>(b, minY - TOLERENCE_MESH_COMP);
    Boostg::set<Boostg::min_corner, 2>(b, minZ - TOLERENCE_MESH_COMP);
    Boostg::set<Boostg::max_corner, 0>(b, maxX + TOLERENCE_MESH_COMP);
    Boostg::set<Boostg::max_corner, 1>(b, maxY + TOLERENCE_MESH_COMP);
    Boostg::set<Boostg::max_corner, 2>(b, maxZ + TOLERENCE_MESH_COMP);
    return;
}

void HandleNekMesh2D::IntersectWithBoxUsingRTree(
    const NekDouble minCornerX, const NekDouble minCornerY,
    const NekDouble minCornerZ, const NekDouble maxCornerX,
    const NekDouble maxCornerY, const NekDouble maxCornerZ, vector<int> &elIds,
    vector<int> &glIds)
{
    RBox b;
    Boostg::set<Boostg::min_corner, 0>(b, minCornerX);
    Boostg::set<Boostg::min_corner, 1>(b, minCornerY);
    Boostg::set<Boostg::min_corner, 2>(b, minCornerZ);
    Boostg::set<Boostg::max_corner, 0>(b, maxCornerX);
    Boostg::set<Boostg::max_corner, 1>(b, maxCornerY);
    Boostg::set<Boostg::max_corner, 2>(b, maxCornerZ);

    // rtree.query()
    std::vector<RValue> result_s;
    m_rtree.query(Boostgi::intersects(b), std::back_inserter(result_s));

    elIds.clear();
    glIds.clear();
    BOOST_FOREACH (RValue const &v, result_s)
    {
        int eId = std::get<1>(v);
        int gId = std::get<2>(v);
        elIds.push_back(eId);
        glIds.push_back(gId);
    }
}

int HandleNekMesh2D::v_GetExpansionIndexUsingRTree(
    const Array<OneD, NekDouble> &point) const
{
    int returnEid = -1;
    RBox b;
    Boostg::set<Boostg::min_corner, 0>(b, point[0] - TOLERENCE_MESH_COMP);
    Boostg::set<Boostg::min_corner, 1>(b, point[1] - TOLERENCE_MESH_COMP);
    Boostg::set<Boostg::min_corner, 2>(b, 0.0 - TOLERENCE_MESH_COMP);
    Boostg::set<Boostg::max_corner, 0>(b, point[0] + TOLERENCE_MESH_COMP);
    Boostg::set<Boostg::max_corner, 1>(b, point[1] + TOLERENCE_MESH_COMP);
    Boostg::set<Boostg::max_corner, 2>(b, 0.0 + TOLERENCE_MESH_COMP);

    vector<unsigned> res;
    MySearchCallback2 callback(res);
    m_rtree.query(Boostgi::intersects(b),
                  boost::make_function_output_iterator(callback));
    std::sort(res.begin(), res.end());
    BOOST_FOREACH (unsigned const &eId, res)
    {
        if (m_expansions[0]->GetExp(eId)->GetGeom()->ContainsPoint(point,
                                                                   TOLERENCE))
        {
            returnEid = eId;
            break;
        }
    }
    return returnEid;
}

bool HandleNekMesh2D::v_InitializeMetricTensor()
{
    m_metricTensor = new MetricTensor();
    m_metricTensor->LoadMetricTensor(this);
    m_MTDefined = true;
    return true;
}

NekDouble HandleNekMesh2D::v_GetDynamicScaling(Array<OneD, NekDouble> glCoord,
                                               int eid, NekDouble mu)
{

    // if eid <0 find a elid.
    if (eid < 0)
    {
        if (m_useRTree)
        {
            Array<OneD, NekDouble> dup_glCoord(2, 0.0); // glCoord--3D--Array
            dup_glCoord[0] = glCoord[0];
            dup_glCoord[1] = glCoord[1];
            eid            = GetExpansionIndexUsingRTree(dup_glCoord);
        }
        else
        {
            Array<OneD, NekDouble> dup_glCoord(2, 0.0); // glCoord--3D--Array
            dup_glCoord[0] = glCoord[0];
            dup_glCoord[1] = glCoord[1];
            eid = m_expansions[0]->GetExpIndex(dup_glCoord, TOLERENCE);
        }
    }
    ASSERTL0(eid >= 0 && "Point out of mesh");
    NekDouble result = -1.0;
    // Get local coordinates.
    // Depending on number of vertices triangle or quad.
    // use locCoordinates as barycentric coordinates
    SpatialDomains::GeometrySharedPtr geomSPtr =
        m_expansions[0]->GetExp(eid)->GetGeom();
    Array<OneD, NekDouble> lCoord(2, 0.0);
    geomSPtr->GetLocCoords(glCoord, lCoord);
    if (geomSPtr->GetShapeType() == Nektar::LibUtilities::eTriangle)
    {
        int Vid0          = geomSPtr->GetVid(0);
        int Vid1          = geomSPtr->GetVid(1);
        int Vid2          = geomSPtr->GetVid(2);
        NekDouble lambda1 = (lCoord[0] + 1.0) / 2.0;
        NekDouble lambda2 = (lCoord[1] + 1.0) / 2.0;
        NekDouble lambda0 = 1.0 - lambda1 - lambda2;
        result            = lambda0 * m_dynVertScaling[Vid0] +
                 lambda1 * m_dynVertScaling[Vid1] +
                 lambda2 * m_dynVertScaling[Vid2];
    }
    else if (geomSPtr->GetShapeType() == Nektar::LibUtilities::eQuadrilateral)
    {
        int Vid0     = geomSPtr->GetVid(0);
        NekDouble v0 = m_dynVertScaling[Vid0];
        int Vid1     = geomSPtr->GetVid(1);
        NekDouble v1 = m_dynVertScaling[Vid1];
        int Vid2     = geomSPtr->GetVid(2);
        NekDouble v2 = m_dynVertScaling[Vid2];
        int Vid3     = geomSPtr->GetVid(3);
        NekDouble v3 = m_dynVertScaling[Vid3];
        NekDouble a  = (lCoord[0] + 1.0) / 2.0;
        NekDouble b  = (lCoord[1] + 1.0) / 2.0;
        result =
            (a * v1 + (1 - a) * v0) * (1 - b) + (a * v2 + (1 - a) * v3) * b;
    }
    else
    {
        NEKERROR(ErrorUtil : efatal, "Shape not accounted for");
    }

    return mu * result;
}

bool HandleNekMesh2D::v_GetMTScalingOfGIDs(vector<int> &t_GIDs,
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
    }
    return true;
}

} // namespace LSIAC
} // namespace Nektar
