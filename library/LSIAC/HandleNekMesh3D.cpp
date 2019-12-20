///////////////////////////////////////////////////////////////////////////////
//
// File HandleNekMesh3D.cpp
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
// Description: HandleNekMesh3D definition
//
///////////////////////////////////////////////////////////////////////////////
#include "HandleNekMesh3D.h"

#include <MultiRegions/ContField3D.h>
#include <MultiRegions/DisContField3D.h>
//#include <SpatialDomains/MeshGraph3D.h>
#include <boost/timer.hpp>
#include <cmath>
//#include <omp.h>

namespace Nektar
{
namespace LSIAC
{
NekDouble HandleNekMesh3D::v_GetElLargestEdgeSize(const NekDouble ptsx,
                                                  const NekDouble ptsy,
                                                  const NekDouble ptsz,
                                                  int Elid)
{
    if (Elid < 0)
    {
        // Find Element Id.
        Array<OneD, NekDouble> glCord(3, 0.0);
        glCord[0] = ptsx;
        glCord[1] = ptsy;
        glCord[2] = ptsz;
        Elid      = m_expansions[0]->GetExpIndex(glCord, TOLERENCE);
        assert(Elid > 0 && "Something wrong. Point outside boundary");
    }

    SpatialDomains::GeometrySharedPtr gEl =
        (m_expansions[0]->GetExp(Elid))->GetGeom();

    int numEdges        = gEl->GetNumEdges();
    NekDouble maxLength = -1.0;
    // Array<OneD,NekDouble> p1loc(3,0.0),p1loc(3,0.0);
    for (int i = 0; i < numEdges; i++)
    {
        int Eid = gEl->GetEid(i);
        // int Vid0  = m_segMap.find(Eid)->second->GetVid(0);
        // int Vid1  = m_segMap.find(Eid)->second->GetVid(1);
        SpatialDomains::PointGeomSharedPtr p0 =
            m_segMap.find(Eid)->second->GetVertex(0);
        SpatialDomains::PointGeomSharedPtr p1 =
            m_segMap.find(Eid)->second->GetVertex(1);
        SpatialDomains::PointGeom p2(*p1);
        // PointGeomSharedPtr p0 = m_pointMap.find(Vid0);
        // PointGeomSharedPtr p1 = m_pointMap.find(Vid1);
        NekDouble len = p0->dist(p2);
        if (len > maxLength)
        {
            maxLength = len;
        }
    }
    assert(maxLength > 0 && "max Length > 0 ");
    return maxLength;
}

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
bool HandleNekMesh3D::v_LoadMesh(string var)
{
    SpatialDomains::ExpansionMap expansions = m_graph->GetExpansions();
    //	cout << "expansion size: " <<expansions.size() << endl;
    m_expansions.push_back(
        MemoryManager<MultiRegions::DisContField3D>::AllocateSharedPtr(
            m_session, m_graph, var));
    return true;
}

bool HandleNekMesh3D::v_LoadData(string filename, vector<string> &variables)
{
    SpatialDomains::ExpansionMap expansions = m_graph->GetExpansions();
    //	cout << "expansion size: " <<expansions.size() << endl;
    for (int i = 0; i < variables.size(); i++)
    {
        m_expansions.push_back(
            MemoryManager<MultiRegions::DisContField3D>::AllocateSharedPtr(
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

    /*
            // want to check if the file was loaded successfully.
            cout << "variables.size: " << variables.size()<< endl;
            for ( int i =0; i < variables.size(); i++)
            {
                    cout << variables[i] << endl;
            }
            Array<OneD,NekDouble> uExp_Coeffs =  m_expansions[0]->GetCoeffs();
            Array<OneD,NekDouble> uExp_Phys =  m_expansions[0]->GetPhys();

            cout << "coefficients" << endl;
            for ( auto c : uExp_Coeffs)
            {
                    cout << c <<"\t" ;
            }
            cout << endl;
            cout << "Phys" << endl;
            for (auto c: uExp_Phys)
            {
                    cout << c << "\t";
            }
            cout << endl;
    */
    return true;
}

// This function does not verify if tmin to tmax is properly given.
// This function in 3D does not fill xpos,ypos and zpos yet.
bool HandleNekMesh3D::v_GetBreakPts(
    const NekDouble xcen_offset, const NekDouble ycen_offset,
    const NekDouble zcen_offset, const Array<OneD, NekDouble> &direction,
    const NekDouble tmin, const NekDouble tmax, vector<NekDouble> &xPos,
    vector<NekDouble> &yPos, vector<NekDouble> &zPos, vector<NekDouble> &tPos)
{
    boost::ignore_unused(xPos, yPos, zPos); // only tPos is used.
    Array<OneD, NekDouble> point(3);
    point[0] = xcen_offset;
    point[1] = ycen_offset;
    point[2] = zcen_offset;
    // IntersectWithEdges ( m_segMap, m_pointMap, direction, point, tmin, tmax,
    // tPos);
    if (m_useRTree)
    {
        IntersectWithFacesUsingRTree(direction, point, tmin, tmax, tPos);
    }
    else
    {
        IntersectWithFaces(direction, point, tmin, tmax, tPos);
    }

    // cout << tPos.size()<< endl;
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

bool HandleNekMesh3D::v_CanTRangebeApplied(
    const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
    const NekDouble scaling, const NekDouble tmin, const NekDouble tmax,
    NekDouble &tminUpdate, NekDouble &tmaxUpdate)
{
    boost::ignore_unused(ptsX, ptsY, ptsZ, scaling, tmin, tmax, tminUpdate,
                         tmaxUpdate);
    NEKERROR(ErrorUtil::efatal, "Not yet coded");
    return false;
}

bool HandleNekMesh3D::v_CanTRangebeAppliedWOMeshShift(
    const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
    const Array<OneD, NekDouble> &direction, const NekDouble tmin,
    const NekDouble tmax)
{
    // pointer to expansion list.
    // Claculate left point.
    Array<OneD, NekDouble> pl(3), pr(3);
    vector<NekDouble> tPos;
    pl[0] = ptsX + direction[0] * tmin;
    pl[1] = ptsY + direction[1] * tmin;
    pl[2] = ptsZ + direction[2] * tmin;
    pr[0] = ptsX + direction[0] * tmax;
    pr[1] = ptsY + direction[1] * tmax;
    pr[2] = ptsZ + direction[2] * tmax;

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

    // cout << pl_index << endl;
    // cout << pr_index << endl;
    if ((0 > pl_index) || (0 > pr_index))
    { // goes out of boundary
        return false;
    }
    return true;
}

bool HandleNekMesh3D::v_CanTRangebeApplied(
    const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
    const Array<OneD, NekDouble> &direction, const NekDouble tmin,
    const NekDouble tmax, NekDouble &meshTShift)
{
    // pointer to expansion list.
    // Claculate left point.
    Array<OneD, NekDouble> pl(3), pr(3);
    vector<NekDouble> tPos;
    pl[0] = ptsX + direction[0] * tmin;
    pl[1] = ptsY + direction[1] * tmin;
    pl[2] = ptsZ + direction[2] * tmin;
    pr[0] = ptsX + direction[0] * tmax;
    pr[1] = ptsY + direction[1] * tmax;
    pr[2] = ptsZ + direction[2] * tmax;
    // int pl_index = m_expansions[0]->GetExpIndex(pl);
    // int pr_index = m_expansions[0]->GetExpIndex(pr);
    int pl_index, pr_index;
    if (m_useRTree)
    {
        pl_index = GetExpansionIndexUsingRTree(pl);
        pr_index = GetExpansionIndexUsingRTree(pr);
        // pl_index = m_expansions[0]->GetExpIndex(pl,TOLERENCE);
        // pr_index = m_expansions[0]->GetExpIndex(pr,TOLERENCE);
    }
    else
    {
        pl_index = m_expansions[0]->GetExpIndex(pl, TOLERENCE);
        pr_index = m_expansions[0]->GetExpIndex(pr, TOLERENCE);
    }

    // cout << pl_index << endl;
    // cout << pr_index << endl;
    if ((0 > pl_index) || (0 > pr_index))
    { // goes out of boundary
        Array<OneD, NekDouble> point(3);
        point[0] = ptsX;
        point[1] = ptsY;
        point[2] = ptsZ;
        // IntersectWithEdges ( m_segMap, m_pointMap, direction, point, tmin,
        // tmax, tPos);
        if (m_useRTree)
        {
            IntersectWithFacesUsingRTree(direction, point, tmin, tmax, tPos);
        }
        else
        {
            IntersectWithFaces(direction, point, tmin, tmax, tPos);
        }
        // A bug can be present here. After the shift. The filter might be out
        // of the mesh again.
        if ((0 > pl_index))
        {
            meshTShift = tPos.front() - tmin;
        }
        if ((0 > pr_index))
        {
            meshTShift = tPos.back() - tmax;
        }
        if ((0 > pl_index) && (0 > pr_index))
        {
            assert(false && "Not enough mesh size to apply.");
        }
        return false;
    }
    // assert( false && "Need more coding" );
    return true;
}

bool HandleNekMesh3D::v_EvaluateAt(const NekDouble xPos, const NekDouble yPos,
                                   const NekDouble zPos, int gID, int eID,
                                   NekDouble &value, int varNum)
{
    boost::ignore_unused(xPos, yPos, zPos, gID, eID, value, varNum);
    NEKERROR(ErrorUtil::efatal, "Not yet coded");
    return false;
}

bool HandleNekMesh3D::v_EvaluateAt(const Array<OneD, NekDouble> &xPos,
                                   const Array<OneD, NekDouble> &yPos,
                                   const Array<OneD, NekDouble> &zPos,
                                   const int gID, const int eID,
                                   Array<OneD, NekDouble> &values, int varNum)
{
    // The reason for asking gID will be useful if we are using MPI.
    boost::ignore_unused(gID); // reserved for global id if implemented.
    assert(gID >= 0 && eID >= 0 && "Input paramerters are out of scope;");
    LocalRegions::ExpansionSharedPtr lexp = m_expansions[0]->GetExp(eID);
    const int phys_offset = m_expansions[0]->GetPhys_Offset(eID);

    const Array<OneD, NekDouble> el_Phys =
        m_Arrays[varNum].CreateWithOffset(m_Arrays[varNum], phys_offset);

    Array<OneD, NekDouble> glCoord(3);
    for (int i = 0; i < values.num_elements(); i++)
    {
        glCoord[0] = xPos[i];
        glCoord[1] = yPos[i];
        glCoord[2] = zPos[i];
        values[i]  = lexp->PhysEvaluate(glCoord, el_Phys);
    }
    return true;
}

bool HandleNekMesh3D::v_GetListOfGIDs(
    const NekDouble xPos, const NekDouble yPos, const NekDouble zPos,
    const Array<OneD, NekDouble> &direction, const vector<NekDouble> t_breaks,
    vector<int> &t_GIDs, vector<int> &t_EIDs) const
{
    t_GIDs.clear();
    t_EIDs.clear();
    t_GIDs.resize(t_breaks.size());
    t_EIDs.resize(t_breaks.size());
    // SpatialDomains::ExpansionMap expansions = m_graph->GetExpansions();
    for (int i = 0; i < t_breaks.size() - 1; i++)
    {
        t_GIDs[i]         = -1;
        t_EIDs[i]         = -1;
        NekDouble t_break = (t_breaks[i] + t_breaks[i + 1]) / 2.0;
        Array<OneD, NekDouble> locCoord(3, 0.0);
        locCoord[0] = xPos + t_break * direction[0];
        locCoord[1] = yPos + t_break * direction[1];
        locCoord[2] = zPos + t_break * direction[2];
        if (m_useRTree)
        {
            t_GIDs[i] = GetExpansionIndexUsingRTree(locCoord);
            // t_GIDs[i] = m_expansions[0]->GetExpIndex(locCoord,TOLERENCE);
        }
        else
        {
            t_GIDs[i] = m_expansions[0]->GetExpIndex(locCoord, TOLERENCE);
        }
        if (t_GIDs[i] == -1)
        {
            t_GIDs[i] = m_expansions[0]->GetExpIndex(locCoord, TOLERENCE);
        }
        t_EIDs[i] = t_GIDs[i];
        if (t_GIDs[i] < 0)
        {
            cout << "Somehting is wrong" << endl;
            cout << "t_breaks" << endl;
            printNekArray(t_breaks);
            cout << std::setprecision(29) << "xPos: " << xPos
                 << " yPos: " << yPos << " zPos: " << zPos << endl;
            cout << std::setprecision(29) << "xPos: " << locCoord[0]
                 << " yPos: " << locCoord[1] << " zPos: " << locCoord[2]
                 << endl;
            cout << "t_GIDs" << endl;
            printNekArray(t_GIDs);
        }
        assert(t_GIDs[i] >= 0 && "Will fail down the line");
    }

    return true;
}

// All private functions doing the bulk of calcualtions.
vector<NekDouble> HandleNekMesh3D::cross_Math(const vector<NekDouble> &r,
                                              const vector<NekDouble> &s)
{
    vector<NekDouble> ans(3);
    ans[0] = r[1] * s[2] - r[2] * s[1];
    ans[1] = r[2] * s[0] - r[0] * s[2];
    ans[2] = r[0] * s[1] - r[1] * s[0];
    return ans;
}

vector<NekDouble> HandleNekMesh3D::sub_Math(vector<NekDouble> &p2,
                                            vector<NekDouble> &p1)
{
    vector<NekDouble> r(3);
    r[0] = p2[0] - p1[0];
    r[1] = p2[1] - p1[1];
    r[2] = p2[2] - p1[2];
    return r;
}
NekDouble HandleNekMesh3D::dot_Math(vector<NekDouble> &p, vector<NekDouble> &q)
{
    return p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
}

NekDouble HandleNekMesh3D::norm2_Math(vector<NekDouble> p)
{
    return p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
}

bool HandleNekMesh3D::intersect(vector<NekDouble> &p1, vector<NekDouble> &p2,
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
            // cout << "p7" << endl;
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
            // cout << "p8" << endl;
            return true;
        }
        else if (t0 <= 0.0 && t1 >= 1.0)
        {
            i1.push_back(p1[0]);
            i1.push_back(p1[1]);
            i1.push_back(p1[2]);
            i2.push_back(p2[0]);
            i2.push_back(p2[1]);
            // cout << "p11" << endl;
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
            // cout << "p9" << endl;
            return true;
        }
        else // one case left t0>1 and t1>1
        {    // ignore
            // cout << "p10" << endl;
            return false;
        }
    }
    if ((norm2_Math(rCs) < TOLERENCE) && (norm2_Math(pMq_Cr) > TOLERENCE))
    {
        // cout << "p3" << endl;
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
            //  cout << t << "\t" << u << "\t"<< r[0] << "\t" << r[1] << "\t" <<
            //  r[2] << endl; cout << "p4" << endl;
            return true;
        }
        else
        {
            //  cout << "p5" << endl;
            return false;
        }
    }
    // cout << "p6" << endl;
    return true;
}

void HandleNekMesh3D::IntersectWithEdges(
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
    // iterate through all the edges.
    // pick a direction not zero.
    for (int i = 0; i < 3; i++)
    {
        if (std::abs(dir[i]) > TOLERENCE)
        {
            dirID = i;
            break;
        }
    }
    assert(dirID >= 0 && "Direction is not right something is up ");
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
            // nbc.printNekArray(i1,0);
            NekDouble t = 0;
            t           = (i1[dirID] - point[dirID]) / dir[dirID];
            tvalT.push_back(t);
            if (i2.size() > 0)
            {
                //  nbc.printNekArray(i2,0);
                NekDouble t = 0;
                t           = (i2[dirID] - point[dirID]) / dir[dirID];
                tvalT.push_back(t);
            }
        }
    }

    // sort and keep them unique
    std::sort(tvalT.begin(), tvalT.end());
    std::vector<NekDouble>::iterator ittt;
    ittt = std::unique(tvalT.begin(), tvalT.end(), compare2NekDoublesH);
    tvalT.resize(std::distance(tvalT.begin(), ittt));
}

void HandleNekMesh3D::FindElementIDForLineSegs(
    const vector<NekDouble> &tvalT, const Array<OneD, NekDouble> &point,
    const Array<OneD, NekDouble> &dir,
    const SpatialDomains::MeshGraphSharedPtr mesh_graph, vector<int> &EIDs)
{
    boost::ignore_unused(mesh_graph, EIDs);
    Array<OneD, NekDouble> temp(3);
    // for every point iteration
    for (int i = 0; i < tvalT.size() - 1; i++)
    {
        temp[0] = point[0] + 0.5 * (tvalT[i] + tvalT[i + 1]) * dir[0];
        temp[1] = point[1] + 0.5 * (tvalT[i] + tvalT[i + 1]) * dir[1];
        temp[2] = point[2] + 0.5 * (tvalT[i] + tvalT[i + 1]) * dir[2];
    }
    NEKERROR(ErrorUtil::efatal, "Not yet coded");
}

void HandleNekMesh3D::IntersectWithFaces(const Array<OneD, NekDouble> &dir,
                                         const Array<OneD, NekDouble> &point,
                                         const NekDouble t1, const NekDouble t2,
                                         vector<NekDouble> &tvalT)
{
    tvalT.clear();
    // Go through all the Triangles;
    // Go through all the quad faces;
    for (int t = 0; t < m_triMap.size(); t++)
    {
        SpatialDomains::TriGeomSharedPtr triPtr =
            m_triMap.find(t)->second; //->second;
        IntersectLineSegWithFace(triPtr, dir, point, t1, t2, tvalT);
    }
    for (int q = 0; q < m_quadMap.size(); q++)
    {
        SpatialDomains::QuadGeomSharedPtr quadPtr = m_quadMap.find(q)->second;
        IntersectLineSegWithFace(quadPtr, dir, point, t1, t2, tvalT);
    }

    // sort and keep them unique
    std::sort(tvalT.begin(), tvalT.end());
    std::vector<NekDouble>::iterator ittt;
    ittt = std::unique(tvalT.begin(), tvalT.end(), compare2NekDoublesH);
    tvalT.resize(std::distance(tvalT.begin(), ittt));
    // End points are not added if not in mesh.
}

void HandleNekMesh3D::PlaneEquationOfFace(
    const SpatialDomains::GeometrySharedPtr geomEl, NekDouble &a, NekDouble &b,
    NekDouble &c, NekDouble &d)
{
    Array<OneD, NekDouble> P0(3, 0.0), P1(3, 0.0), P2(3, 0.0);
    SpatialDomains::PointGeomSharedPtr psh0 = geomEl->GetVertex(0);
    SpatialDomains::PointGeomSharedPtr psh1 = geomEl->GetVertex(1);
    SpatialDomains::PointGeomSharedPtr psh2 = geomEl->GetVertex(2);
    psh0->GetCoords(P0);
    psh1->GetCoords(P1);
    psh2->GetCoords(P2);
    // normal calculation a,b,c
    vector<NekDouble> v1(3, 0.0), v2(3, 0.0);
    v1[0] = P1[0] - P0[0];
    v1[1] = P1[1] - P0[1];
    v1[2] = P1[2] - P0[2];
    v2[0] = P2[0] - P0[0];
    v2[1] = P2[1] - P0[1];
    v2[2] = P2[2] - P0[2];
    a     = v1[1] * v2[2] - v1[2] * v2[1];
    b     = v1[2] * v2[0] - v1[0] * v2[2];
    c     = v1[0] * v2[1] - v1[1] * v2[0];
    d     = -1.0 * (a * P0[0] + b * P0[1] + c * P0[2]);
    return;
}

void HandleNekMesh3D::IntersectLineSegWithFace(
    const SpatialDomains::GeometrySharedPtr geomEl,
    const Array<OneD, NekDouble> &dir, const Array<OneD, NekDouble> &point,
    const NekDouble t1, const NekDouble t2, vector<NekDouble> &tvalT)
{
    // Input is assumed to be face.
    // Form a plane equation interms of abcd;
    // Get plane equation of face.
    // tvalT.clear();   // Assuming tvalT is collecting all the intersection.
    NekDouble a, b, c, d;
    PlaneEquationOfFace(geomEl, a, b, c, d);
    vector<NekDouble> p1(3), p2(3), i1(3), i2(3);
    p1[0] = point[0] + t1 * dir[0];
    p1[1] = point[1] + t1 * dir[1];
    p1[2] = point[2] + t1 * dir[2];
    p2[0] = point[0] + t2 * dir[0];
    p2[1] = point[1] + t2 * dir[1];
    p2[2] = point[2] + t2 * dir[2];
    // Point and direction are point an ddir.
    NekDouble Vd = a * dir[0] + b * dir[1] + c * dir[2];
    NekDouble Vn = a * point[0] + b * point[1] + c * point[2] + d;
    // case 1.
    if (((std::abs(Vd) > TOLERENCE) && (std::abs(Vn) > TOLERENCE)) ||
        ((std::abs(Vd) > TOLERENCE) && (std::abs(Vn) < TOLERENCE)))
    {
        // calcualte point.
        NekDouble t = -1.0 * Vn / Vd;
        // check if point is between t1 and t2.
        if (t >= t1 && t <= t2)
        {
            // Form physCoord and check using ContainsPoint if it is present.
            Array<OneD, NekDouble> gloC(3, 0.0);
            gloC[0] = point[0] + t * dir[0];
            gloC[1] = point[1] + t * dir[1];
            gloC[2] = point[2] + t * dir[2];
            // Then check if the point is on the face.
            if (geomEl->ContainsPoint(gloC, TOLERENCE))
            { // point on face. Add to tvalT and exit with true.
                tvalT.push_back(t);
                // cout << "Check point out 1" << endl;
                return;
            }
            else
            {
                // cout << "Check point out 2" << endl;
                return;
            }
        }
    }
    else if ((std::abs(Vd) < TOLERENCE) && (std::abs(Vn) > TOLERENCE))
    {
        // cout << "Check point out 3" << endl;
        return;
        // the line segment is parallel to plane but present at an offset.
    }
    else if ((std::abs(Vd) < TOLERENCE) && (std::abs(Vn) < TOLERENCE))
    {
        // go through all the edges of the element and check for intersections.
        // Accumulate all the intersections.
        // Do end point checks. Normally end points are added at end. So not
        // needed.
        int numE  = geomEl->GetNumEdges();
        int dirID = -1;
        for (int i = 0; i < 3; i++)
        {
            if (std::abs(dir[i]) > TOLERENCE)
            {
                dirID = i;
                break;
            }
        }
        if (dirID < 0)
        {
            cout << "something is wroing" << endl;
        }
        assert(dirID >= 0 && "Direction is not right something is up ");

        for (int e = 0; e < numE; e++)
        {
            int Eid = geomEl->GetEid(e);
            SpatialDomains::SegGeomSharedPtr segPtr =
                m_segMap.find(Eid)->second; //->second;
            vector<NekDouble> q1(3), q2(3);
            m_pointMap.find(segPtr->GetVid(0))
                ->second->GetCoords(q1[0], q1[1], q1[2]);
            m_pointMap.find(segPtr->GetVid(1))
                ->second->GetCoords(q2[0], q2[1], q2[2]);
            bool b = intersect(p1, p2, q1, q2, i1, i2);
            if (b)
            {
                // nbc.printNekArray(i1,0);
                NekDouble t = 0;
                t           = (i1[dirID] - point[dirID]) / dir[dirID];
                tvalT.push_back(t);
                if (i2.size() > 0)
                {
                    //  nbc.printNekArray(i2,0);
                    NekDouble t = 0;
                    t           = (i2[dirID] - point[dirID]) / dir[dirID];
                    tvalT.push_back(t);
                }
            }
        }
        // point is in plane of face. Need to work with all the edges seperatly.
        return;
    }
}

void HandleNekMesh3D::v_LoadExpListIntoRTree()
{
    // check if expansions are already loaded.
    assert(m_expansions.size() > 0 &&
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

void HandleNekMesh3D::GetBoundingOfElement(
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

void HandleNekMesh3D::BoundingBoxOfLineSeg(const Array<OneD, NekDouble> &dir,
                                           const Array<OneD, NekDouble> &pt,
                                           const NekDouble t1,
                                           const NekDouble t2, RBox &b)
{
    NekDouble pl_0 = pt[0] + dir[0] * t1;
    NekDouble pl_1 = pt[1] + dir[1] * t1;
    NekDouble pl_2 = pt[2] + dir[2] * t1;

    NekDouble pr_0 = pt[0] + dir[0] * t2;
    NekDouble pr_1 = pt[1] + dir[1] * t2;
    NekDouble pr_2 = pt[2] + dir[2] * t2;

    Boostg::set<Boostg::min_corner, 0>(b, getMin(pl_0, pr_0));
    Boostg::set<Boostg::min_corner, 1>(b, getMin(pl_1, pr_1));
    Boostg::set<Boostg::min_corner, 2>(b, getMin(pl_2, pr_2));
    Boostg::set<Boostg::max_corner, 0>(b, getMax(pl_0, pr_0));
    Boostg::set<Boostg::max_corner, 1>(b, getMax(pl_1, pr_1));
    Boostg::set<Boostg::max_corner, 2>(b, getMax(pl_2, pr_2));
}
/*
// used with Boost.Geometry R-tree
struct MySearchCallback2
{
    template <typename Value>
    void operator()(Value const& v)
    {
        res.push_back(v.second);
    }
};
*/

void HandleNekMesh3D::IntersectWithFacesUsingRTree(
    const Array<OneD, NekDouble> &dir, const Array<OneD, NekDouble> &pt,
    const NekDouble t1, const NekDouble t2, vector<NekDouble> &tvalT)
{
    RBox b;
    BoundingBoxOfLineSeg(dir, pt, t1, t2, b);

    // MySearchCallback2 callback;
    vector<unsigned> res;
    MySearchCallback2 callback(res);

    m_rtree.query(Boostgi::intersects(b),
                  boost::make_function_output_iterator(callback));

    BOOST_FOREACH (unsigned const &eId, res)
    {
        Nektar::SpatialDomains::Geometry3D *geomSPtr =
            dynamic_cast<Nektar::SpatialDomains::Geometry3D *>(
                (m_expansions[0]->GetExp(eId)->GetGeom()).get());
        for (int f = 0; f < geomSPtr->GetNumFaces(); ++f)
        {
            SpatialDomains::GeometrySharedPtr facePtr = geomSPtr->GetFace(f);
            IntersectLineSegWithFace(facePtr, dir, pt, t1, t2, tvalT);
        }
    }

    // sort and keep them unique
    std::sort(tvalT.begin(), tvalT.end());
    std::vector<NekDouble>::iterator ittt;
    ittt = std::unique(tvalT.begin(), tvalT.end(), compare2NekDoublesH);
    tvalT.resize(std::distance(tvalT.begin(), ittt));
}

int HandleNekMesh3D::v_GetExpansionIndexUsingRTree(
    const Array<OneD, NekDouble> &point) const
{
    int returnEid = -1;
    RBox b;
    // BoundingBoxOfLineSeg( dir, pt, t1, t2,b);
    Boostg::set<Boostg::min_corner, 0>(b, point[0] - TOLERENCE);
    Boostg::set<Boostg::min_corner, 1>(b, point[1] - TOLERENCE);
    Boostg::set<Boostg::min_corner, 2>(b, point[2] - TOLERENCE);
    Boostg::set<Boostg::max_corner, 0>(b, point[0] + TOLERENCE);
    Boostg::set<Boostg::max_corner, 1>(b, point[1] + TOLERENCE);
    Boostg::set<Boostg::max_corner, 2>(b, point[2] + TOLERENCE);

    vector<unsigned> res;
    MySearchCallback2 callback(res);
    // boost::timer t;
    m_rtree.query(Boostgi::intersects(b),
                  boost::make_function_output_iterator(callback));
    // double s = t.elapsed();
    //    std::cout <<"time elapsed \t "<< s << std::endl;
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

void HandleNekMesh3D::IntersectWithBoxUsingRTree(
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
    // assert( false && "Not done developing yet :( ");
    BOOST_FOREACH (RValue const &v, result_s)
    {
        int eId = std::get<1>(v);
        int gId = std::get<2>(v);
        elIds.push_back(eId);
        glIds.push_back(gId);
    }
}

NekDouble HandleNekMesh3D::v_GetDynamicScaling(Array<OneD, NekDouble> glCoord,
                                               int eid, NekDouble mu)
{
    // This will be different depending on the element type.
    if (eid < 0)
    {
        eid = GetExpansionIndexUsingRTree(glCoord);
    }

    assert(eid >= 0 && "Point out of mesh");
    NekDouble result = 0;

    SpatialDomains::GeometrySharedPtr geomSPtr =
        m_expansions[0]->GetExp(eid)->GetGeom();

    Array<OneD, NekDouble> lCoord(3, 0.0);
    geomSPtr->GetLocCoords(glCoord, lCoord);
    if (geomSPtr->GetShapeType() == Nektar::LibUtilities::eTetrahedron)
    {
        int Vid0          = geomSPtr->GetVid(0);
        int Vid1          = geomSPtr->GetVid(1);
        int Vid2          = geomSPtr->GetVid(2);
        int Vid3          = geomSPtr->GetVid(3);
        NekDouble lambda1 = (lCoord[0] + 1.0) / 2.0;
        NekDouble lambda2 = (lCoord[1] + 1.0) / 2.0;
        NekDouble lambda3 = (lCoord[2] + 1.0) / 2.0;
        NekDouble lambda0 = 1.0 - lambda1 - lambda2 - lambda3;
        result            = lambda0 * m_dynVertScaling[Vid0] +
                 lambda1 * m_dynVertScaling[Vid1] +
                 lambda2 * m_dynVertScaling[Vid2] +
                 lambda3 * m_dynVertScaling[Vid3];
    }
    else if (geomSPtr->GetShapeType() == Nektar::LibUtilities::eHexahedron)
    {
        int Vid0     = geomSPtr->GetVid(0);
        NekDouble v0 = m_dynVertScaling[Vid0];
        int Vid1     = geomSPtr->GetVid(1);
        NekDouble v1 = m_dynVertScaling[Vid1];
        int Vid2     = geomSPtr->GetVid(2);
        NekDouble v2 = m_dynVertScaling[Vid2];
        int Vid3     = geomSPtr->GetVid(3);
        NekDouble v3 = m_dynVertScaling[Vid3];
        int Vid4     = geomSPtr->GetVid(4);
        NekDouble v4 = m_dynVertScaling[Vid4];
        int Vid5     = geomSPtr->GetVid(5);
        NekDouble v5 = m_dynVertScaling[Vid5];
        int Vid6     = geomSPtr->GetVid(6);
        NekDouble v6 = m_dynVertScaling[Vid6];
        int Vid7     = geomSPtr->GetVid(7);
        NekDouble v7 = m_dynVertScaling[Vid7];
        NekDouble a  = (lCoord[0] + 1.0) / 2.0;
        NekDouble b  = (lCoord[1] + 1.0) / 2.0;
        NekDouble c  = (lCoord[2] + 1.0) / 2.0;
        result       = (((1.0 - a) * v0 + a * v1) * (1 - b) +
                  ((1.0 - a) * v3 + a * v2) * b) *
                     (1 - c) +
                 (((1.0 - a) * v4 + a * v5) * (1 - b) +
                  ((1.0 - a) * v7 + a * v6) * b) *
                     c;
    }
    else if (geomSPtr->GetShapeType() == Nektar::LibUtilities::ePrism)
    {

        // cout << "Into prism :"<< endl;
        // cout << "Get Num Verts" << geomSPtr->GetNumVerts()<<endl;
        int Vid0 = geomSPtr->GetVid(0);
        // NekDouble v0 = m_dynVertScaling[Vid0];
        int Vid1 = geomSPtr->GetVid(1);
        // NekDouble v1 = m_dynVertScaling[Vid1];
        int Vid2 = geomSPtr->GetVid(2);
        // NekDouble v2 = m_dynVertScaling[Vid2];
        int Vid3 = geomSPtr->GetVid(3);
        // NekDouble v3 = m_dynVertScaling[Vid3];
        int Vid4 = geomSPtr->GetVid(4);
        // NekDouble v4 = m_dynVertScaling[Vid4];
        int Vid5 = geomSPtr->GetVid(5);
        // NekDouble v5 = m_dynVertScaling[Vid5];

        NekDouble a   = (lCoord[0] + 1.0) / 2.0;
        NekDouble b   = (lCoord[2] + 1.0) / 2.0;
        NekDouble c   = 1.0 - a - b;
        NekDouble d   = (lCoord[1] + 1.0) / 2.0;
        NekDouble top = c * m_dynVertScaling[Vid0] +
                        a * m_dynVertScaling[Vid1] + b * m_dynVertScaling[Vid4];
        NekDouble bot = c * m_dynVertScaling[Vid3] +
                        a * m_dynVertScaling[Vid2] + b * m_dynVertScaling[Vid5];
        result = top * (1 - d) + bot * (d);

        // 6 verts are 000 010 001 100 110 101
    }

    // ****************************************************************** ///
    // ***This algorithm has few flaws does not work on hexaderal meshes.///
    // Algo: Using Spatial::GeomFactors. //
    // get coords of all vertices.
    // Create 2 double arrays to calcualte areas.
    // 	 1) For actual coordinates.
    // 	 2) For substituting each coord with given coord.
    // Loop through all the vertices to evaluate barycentric coordinates for all
    // vertices. Loop through all the vertices to evaluate final result.
    // multiply final result by mu.
    // ****************************************************************** ///

    return mu * result;
}

} // namespace LSIAC
} // namespace Nektar
