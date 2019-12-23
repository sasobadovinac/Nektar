///////////////////////////////////////////////////////////////////////////////
//
// File HandleNekMesh.h
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
// Description: HandleNekMesh definition
//
///////////////////////////////////////////////////////////////////////////////
#pragma once
#include "HandleMesh.h"
#include "MetricTensor.h"
#include <FieldUtils/Field.hpp>
#include <FieldUtils/FieldUtilsDeclspec.h>
#include <SpatialDomains/Geometry.h>
#include <memory>

namespace Nektar
{
namespace LSIAC
{

// Forward declartion
class MetricTensor;

/**
 * @brief This class deals with mesh queries of nektarMeshes
 *
 * This class, in particular, can load and save only meshes supported by Nektar.
 */
class HandleNekMesh : public HandleMesh
{
private:
    vector<std::shared_ptr<MultiRegions::ExpList>> Explist;

protected:
    int m_meshdim;
    int m_spacedim;
    bool m_MTDefined;

public:
    // used with Boost.Geometry R-tree
    struct MySearchCallback2
    {
        MySearchCallback2(vector<unsigned> &res) : m_res(&res)
        {
        }

        template <typename Value> void operator()(Value const &v)
        {
            m_res->push_back(std::get<1>(v));
        }

        vector<unsigned> *m_res;
    };

    LibUtilities::SessionReaderSharedPtr m_session;
    LibUtilities::FieldIOSharedPtr m_fld;
    SpatialDomains::MeshGraphSharedPtr m_graph;
    vector<MultiRegions::ExpListSharedPtr> m_expansions;
    SpatialDomains::SegGeomMap m_segMap;
    SpatialDomains::PointGeomMap m_pointMap;
    vector<Array<OneD, NekDouble>> m_Arrays;
    std::map<int, NekDouble> m_dynVertScaling;
    MetricTensor *m_metricTensor;
    FieldUtils::FieldSharedPtr m_f;

    /**
     * @brief Constructor created using FieldShared pointer.
     */
    HandleNekMesh(FieldUtils::FieldSharedPtr fldSharedPtr) : m_f(fldSharedPtr)
    {
        m_graph      = m_f->m_graph;
        m_meshdim    = m_graph->GetMeshDimension();
        m_spacedim   = m_graph->GetSpaceDimension();
        m_segMap     = m_graph->GetAllSegGeoms();
        m_pointMap   = m_graph->GetAllPointGeoms();
        m_MTDefined  = false;
        m_expansions = m_f->m_exp;
    };

    /**
     * @brief Constructor created using SessionReaderShared pointer.
     */
    HandleNekMesh(LibUtilities::SessionReaderSharedPtr sessionPtr)
        : m_session(sessionPtr)
    {
        // cout << "into HandleNekMesh constructol" << endl;
        m_fld       = LibUtilities::FieldIO::CreateDefault(m_session);
        m_graph     = SpatialDomains::MeshGraph::Read(m_session);
        m_meshdim   = m_graph->GetMeshDimension();
        m_spacedim  = m_graph->GetSpaceDimension();
        m_segMap    = m_graph->GetAllSegGeoms();
        m_pointMap  = m_graph->GetAllPointGeoms();
        m_MTDefined = false;
    };

    HandleNekMesh();
    /**
     * @brief Adds an expansion the fields using the field name.
     */
    bool LoadMesh(string var)
    {
        return v_LoadMesh(var);
    }

    /**
     * @brief Loads the data fields from the specified.
     */
    bool LoadData(string FileName, vector<string> &variables)
    {
        return v_LoadData(FileName, variables);
    }

    /**
     * @brief Load expansions into Rtree for speed up.
     */
    void LoadExpListIntoRTree()
    {
        v_LoadExpListIntoRTree();
    }

    /**
     * @brief Get element ID given a point in the mesh using RTree lookup.
     */
    int GetExpansionIndexUsingRTree(const Array<OneD, NekDouble> &point) const
    {
        return v_GetExpansionIndexUsingRTree(point);
    }

    /**
     * @brief Evaluates the value of the field at various points on the mesh.
     *
     * \param xPos Array of x-coordinates for the given points.
     * \param yPos Array of y-coordinates for the given points.
     * \param zPos Array of z-coordinates for the given points.
     * \param gID Global element id of the element.
     * \param eID Local element id of the element.
     * \param values Evaluations are written to this array.
     * \param varNum The index of the field to be evaluated.
     *
     * Note: All the points passed should be present in a single element given
     * by eID and gID.
     */
    bool EvaluateAt(const Array<OneD, NekDouble> &xPos,
                    const Array<OneD, NekDouble> &yPos,
                    const Array<OneD, NekDouble> &zPos, const int gID,
                    const int eID, Array<OneD, NekDouble> &values,
                    int varNum = 0)
    {
        return v_EvaluateAt(xPos, yPos, zPos, gID, eID, values, varNum);
    }

    /**
     * @brief Evaluates the value of the field at a given point in the mesh.
     *
     * \param xPos x-coordinate for the given points.
     * \param yPos y-coordinate for the given points.
     * \param zPos z-coordinate for the given points.
     * \param gID Global element id of the element.
     * \param eID Local element id of the element.
     * \param values Evaluation is written to this variable.
     * \param varNum The index of the field to be evaluated.
     */
    bool EvaluateAt(const NekDouble xPos, const NekDouble yPos,
                    const NekDouble zPos, int gID, int eID, NekDouble &value,
                    int varNum = 0)
    {
        return v_EvaluateAt(xPos, yPos, zPos, gID, eID, value, varNum);
    }

    /**
     * @brief Evaluate mesh intersection with the line segment.
     *
     * This function evaluates points at which the given line segment
     * intersects the element boundaries of the mesh.
     * The point (xcen_offset, ycen_offset, zcen_offset), direction, and bounds
     * (t_offset_min, t_offset_max) define the line segment. The intersections
     * are returned in the parametrized form (tPos) and global locations (xPos,
     * yPos, and zPos).
     *
     * Note: The endpoint tmin and tmax are not included in the set of
     * intersections unless they exist on the mesh boundaries.
     */
    bool GetBreakPts_Without_Tmin_Tmax(
        const NekDouble xcen_offset, const NekDouble ycen_offset,
        const NekDouble zcen_offset, const Array<OneD, NekDouble> &direction,
        const NekDouble t_offset_min, const NekDouble t_offset_max,
        vector<NekDouble> &xPos, vector<NekDouble> &yPos,
        vector<NekDouble> &zPos, vector<NekDouble> &tPos)
    {

        return v_GetBreakPts_Without_Tmin_Tmax(
            xcen_offset, ycen_offset, zcen_offset, direction, t_offset_min,
            t_offset_max, xPos, yPos, zPos, tPos);
    };

    /**
     * @brief Evaluate mesh intersection with the line segment.
     *
     * This function evaluates points at which the given line segment
     * intersects the element boundaries of the mesh.
     * The point (xcen_offset, ycen_offset, zcen_offset), direction, and bounds
     * (t_offset_min, t_offset_max) define the line segment. The intersections
     * are returned in the parametrized form (tPos) and global locations (xPos,
     * yPos, and zPos).
     *
     * Note: The endpoint tmin and tmax are included in the set of
     * intersections even if they are not located at the mesh boundaries.
     */
    bool GetBreakPts(const NekDouble xcen_offset, const NekDouble ycen_offset,
                     const NekDouble zcen_offset,
                     const Array<OneD, NekDouble> &direction,
                     const NekDouble t_offset_min, const NekDouble t_offset_max,
                     vector<NekDouble> &xPos, vector<NekDouble> &yPos,
                     vector<NekDouble> &zPos, vector<NekDouble> &tPos)
    {

        return v_GetBreakPts(xcen_offset, ycen_offset, zcen_offset, direction,
                             t_offset_min, t_offset_max, xPos, yPos, zPos,
                             tPos);
    };

    /**
     * @brief This function determines if the symmetric kernel can be applied.
     *
     * Returns true if the symmetric kernel can be applied. Returns false if
     * the symmetric kernel cannot be used. If returned false, The tminUpdate
     * and tmaxUpdate will contain the closest value to tmin and tmax,
     * respectively, for which the function can return true.
     */
    bool CanTRangebeApplied(const NekDouble ptsX, const NekDouble ptsY,
                            const NekDouble ptsZ, const NekDouble scaling,
                            const NekDouble tmin, const NekDouble tmax,
                            NekDouble &tminUpdate, NekDouble &tmaxUpdate)
    {
        return v_CanTRangebeApplied(ptsX, ptsY, ptsZ, scaling, tmin, tmax,
                                    tminUpdate, tmaxUpdate);
    };

    /**
     * @brief This function determines if the symmetric kernel can be applied.
     *
     * Returns true if the symmetric kernel can be applied. Returns false if
     * the symmetric kernel cannot be applied. If returned false, Meshshift
     * will contain shift required by Meshknots
     *
     * Note: Meshknots are Kernelknots*meshScaling.
     */
    bool CanTRangebeApplied(const NekDouble ptsX, const NekDouble ptsY,
                            const NekDouble ptsZ,
                            const Array<OneD, NekDouble> &direction,
                            const NekDouble tmin, const NekDouble tmax,
                            NekDouble &meshShift)
    {
        return v_CanTRangebeApplied(ptsX, ptsY, ptsZ, direction, tmin, tmax,
                                    meshShift);
    };

    /**
     * @brief This function determines if the symmetric kernel can be applied.
     *
     * Returns true if the symmetric kernel can be applied. Returns false if
     * the symmetric kernel cannot be used.
     * */
    bool CanTRangebeAppliedWOMeshShift(const NekDouble ptsX,
                                       const NekDouble ptsY,
                                       const NekDouble ptsZ,
                                       const Array<OneD, NekDouble> &direction,
                                       const NekDouble tmin,
                                       const NekDouble tmax)
    {
        return v_CanTRangebeAppliedWOMeshShift(ptsX, ptsY, ptsZ, direction,
                                               tmin, tmax);
    };

    /**
     * @brief Research Phase.
     *
     * In the reserach phase, behavior may change in time.
     *
     */
    bool GetMTScalingOfGIDs(vector<int> &t_GIDs,
                            Array<OneD, NekDouble> &direction,
                            vector<NekDouble> &scalings)
    {
        return v_GetMTScalingOfGIDs(t_GIDs, direction, scalings);
    };

    /**
     * @brief Get the list of element ids across a given line from the mesh.
     *
     * The line is defined by the point (xPos, yPos, and zPos) and direction.
     * The parameter t_breaks represent the line segments across the elements,
     * such that each line segment is located entirely within a single element.
     * The t_GIDs and t_EIDs return the global and local element ids of the
     * line segments.
     */
    bool GetListOfGIDs(const NekDouble xPos, const NekDouble yPos,
                       const NekDouble zPos,
                       const Array<OneD, NekDouble> &direction,
                       const vector<NekDouble> t_breaks, vector<int> &t_GIDs,
                       vector<int> &t_EIDs) const
    {
        return v_GetListOfGIDs(xPos, yPos, zPos, direction, t_breaks, t_GIDs,
                               t_EIDs);
    };

    /**
     * @brief This function returns the length of the largest edge within an
     * element.
     *
     * If the  Elid <0, then it finds the element which contains the given
     * point.
     */
    NekDouble GetElLargestEdgeSize(const NekDouble ptsx,
                                   const NekDouble ptsy = 0.0,
                                   const NekDouble ptsz = 0.0, int Elid = -1)
    {
        return v_GetElLargestEdgeSize(ptsx, ptsy, ptsz, Elid);
    };

    /**
     * @brief Returns the jacobian determinant for an element.
     *
     * This function returns jacobian determinant for an element. It also
     * assumes that the jacobian is constant across the element.
     */
    NekDouble GetJacobian(const int eID)
    {
        return v_GetJacobian(eID);
    };

    /**
     * @brief This function returns the length of the largest edge within an
     * element.
     *
     * If the  Elid <0, then it finds the element which contains the given
     * point.
     */
    NekDouble GetLargestEdgeLength(const int eID)
    {
        return v_GetLargestEdgeLength(eID);
    };

    /**
     * @brief This function returns the length of the largest edge in the mesh.
     */
    NekDouble GetMeshLargestEdgeLength()
    {
        return v_GetMeshLargestEdgeLength();
    };

    /**
     * @brief This function performs the preprocessing tasks required to
     * calculate adaptive scaling at any point in the mesh.
     */
    bool CalculateDynamicScaling()
    {
        return v_CalculateDynamicScaling();
    };

    /**
     * @brief This function returns the adaptive characteristic length at a
     * point.
     *
     * If eid <0, then the point (glCoord) is used to calculate the element id.
     */
    NekDouble GetDynamicScaling(Array<OneD, NekDouble> &glCoord, int eid = -1,
                                NekDouble mu = 1.0)
    {
        return v_GetDynamicScaling(glCoord, eid, mu);
    };

    /**
     * @brief Research Phase.
     *
     * In the reserach phase, behavior may change in time.
     *
     */
    bool GetKnotVec(const int degree, const Array<OneD, NekDouble> &coord,
                    const Array<OneD, NekDouble> &direction,
                    Array<OneD, NekDouble> &knotVec, NekDouble &shift)
    {
        return v_GetKnotVec(degree, coord, direction, knotVec, shift);
    }

    /**
     * @brief Research Phase.
     *
     * In the reserach phase, behavior may change in time.
     *
     */
    bool Get1DVec(vector<NekDouble> &coords)
    {
        return v_Get1DVec(coords);
    }

    /**
     * @brief Return the boundaries of the line segments with the mesh.
     *
     * The function returns the boundaries of the line segment in the
     * parametrized form (tmin, tmax), such that the line segment is within the
     * mesh boundaries.
     */
    bool WhatIsTRange(const NekDouble ptsX, const NekDouble ptsY,
                      const NekDouble ptsZ,
                      const Array<OneD, NekDouble> &direction, NekDouble &tmin,
                      NekDouble &tmax, int &num)
    {
        return v_WhatIsTRange(ptsX, ptsY, ptsZ, direction, tmin, tmax, num);
    }

    bool InitializeMetricTensor()
    {
        return v_InitializeMetricTensor();
    }

protected:
    virtual bool v_InitializeMetricTensor()
    {
        assert(false && "v_InitializeMetricTensor");
        return false;
    }

    virtual bool v_GetKnotVec(const int degree,
                              const Array<OneD, NekDouble> &coord,
                              const Array<OneD, NekDouble> &direction,
                              Array<OneD, NekDouble> &knotVec, NekDouble &shift)
    {
        boost::ignore_unused(degree, coord, direction, knotVec, shift);
        NEKERROR(ErrorUtil::efatal, "v_GetKnotVec");
        return false;
    }
    virtual bool v_Get1DVec(vector<NekDouble> &coords)
    {
        boost::ignore_unused(coords);
        NEKERROR(ErrorUtil::efatal, "v_Get1DVec");
        return false;
    }

    virtual int v_GetExpansionIndexUsingRTree(
        const Array<OneD, NekDouble> &point) const
    {
        boost::ignore_unused(point);
        NEKERROR(ErrorUtil::efatal, "v_GetExpansionIndexUsingRTree");
        return -1;
    }

    virtual bool v_GetBreakPts_Without_Tmin_Tmax(
        const NekDouble xcen_offset, const NekDouble ycen_offset,
        const NekDouble zcen_offset, const Array<OneD, NekDouble> &direction,
        const NekDouble t_offset_min, const NekDouble t_offset_max,
        vector<NekDouble> &xPos, vector<NekDouble> &yPos,
        vector<NekDouble> &zPos, vector<NekDouble> &tPos)
    {
        boost::ignore_unused(xcen_offset, ycen_offset, zcen_offset, direction,
                             t_offset_min, t_offset_max, xPos, yPos, zPos,
                             tPos);
        NEKERROR(ErrorUtil::efatal, "v_GetBreakPts_Without_Tmin_Tmax");
        return false;
    };

    virtual bool v_GetBreakPts(const NekDouble xcen_offset,
                               const NekDouble ycen_offset,
                               const NekDouble zcen_offset,
                               const Array<OneD, NekDouble> &direction,
                               const NekDouble t_offset_min,
                               const NekDouble t_offset_max,
                               vector<NekDouble> &xPos, vector<NekDouble> &yPos,
                               vector<NekDouble> &zPos, vector<NekDouble> &tPos)
    {
        boost::ignore_unused(xcen_offset, ycen_offset, zcen_offset, direction,
                             t_offset_min, t_offset_max, xPos, yPos, zPos,
                             tPos);
        NEKERROR(ErrorUtil::efatal, "v_GetBreakPts");
        return false;
    };
    virtual bool v_CanTRangebeApplied(
        const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
        const NekDouble scaling, const NekDouble tmin, const NekDouble tmax,
        NekDouble &tminUpdate, NekDouble &tmaxUpdate)
    {
        boost::ignore_unused(ptsX, ptsY, ptsZ, scaling, tmin, tmax, tminUpdate,
                             tmaxUpdate);
        NEKERROR(ErrorUtil::efatal, "v_CanTRangebeApplied");
        return false;
    };

    virtual bool v_CanTRangebeApplied(
        const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
        const Array<OneD, NekDouble> &direction, const NekDouble tmin,
        const NekDouble tmax, NekDouble &meshShift)
    {
        boost::ignore_unused(ptsX, ptsY, ptsZ, direction, tmin, tmax,
                             meshShift);
        NEKERROR(ErrorUtil::efatal, "v_CanTRangebeApplied");
        return false;
    };
    virtual bool v_CanTRangebeAppliedWOMeshShift(
        const NekDouble ptsX, const NekDouble ptsY, const NekDouble ptsZ,
        const Array<OneD, NekDouble> &direction, const NekDouble tmin,
        const NekDouble tmax)
    {
        boost::ignore_unused(ptsX, ptsY, ptsZ, direction, tmin, tmax);
        NEKERROR(ErrorUtil::efatal, "v_CanTRangebeAppliedWOMeshShift");
        return false;
    };
    virtual bool v_EvaluateAt(const Array<OneD, NekDouble> &xPos,
                              const Array<OneD, NekDouble> &yPos,
                              const Array<OneD, NekDouble> &zPos, const int gID,
                              const int eID, Array<OneD, NekDouble> &values,
                              int varNum = 0)
    {
        boost::ignore_unused(xPos, yPos, zPos, gID, eID, values, varNum);
        NEKERROR(ErrorUtil::efatal, "v_EvaluateAt");
        return false;
    }
    virtual bool v_EvaluateAt(const NekDouble xPos, const NekDouble yPos,
                              const NekDouble zPos, int gID, int eID,
                              NekDouble &value, int varNum = 0)
    {
        boost::ignore_unused(xPos, yPos, zPos, gID, eID, value, varNum);
        NEKERROR(ErrorUtil::efatal, "v_EvaluateAt");
        return false;
    }

    virtual bool v_LoadMesh(string var)
    {
        boost::ignore_unused(var);
        NEKERROR(ErrorUtil::efatal, "v_LoadMesh");
        return false;
    }
    virtual bool v_LoadData(string Filename, vector<string> &variables)
    {
        boost::ignore_unused(Filename, variables);
        NEKERROR(ErrorUtil::efatal, "v_LoadData");
        return false;
    }
    virtual bool v_GetListOfGIDs(const NekDouble xPos, const NekDouble yPos,
                                 const NekDouble zPos,
                                 const Array<OneD, NekDouble> &direction,
                                 const vector<NekDouble> t_breaks,
                                 vector<int> &t_GIDs, vector<int> &t_EIDs) const
    {
        boost::ignore_unused(xPos, yPos, zPos, direction, t_breaks, t_GIDs,
                             t_EIDs);
        NEKERROR(ErrorUtil::efatal, "v_GetListOfGIDs");
        return false;
    }
    virtual NekDouble v_GetElLargestEdgeSize(const NekDouble ptsx,
                                             const NekDouble ptsy,
                                             const NekDouble ptsz,
                                             int Elid = -1)
    {
        boost::ignore_unused(ptsx, ptsy, ptsz, Elid);
        NEKERROR(ErrorUtil::efatal, "v_GetElLargestEdgeSize");
        return false;
    }
    static bool compare2NekDoublesH(NekDouble x, NekDouble y)
    {
        return ((x - TOLERENCE < y) && (x + TOLERENCE > y));
    }

    virtual void v_LoadExpListIntoRTree()
    {
        NEKERROR(ErrorUtil::efatal, "v_LoadExpListIntoRTree");
        return false;
    }

    virtual NekDouble v_GetJacobian(const int eID)
    {
        SpatialDomains::GeometrySharedPtr geom =
            m_expansions[0]->GetExp(eID)->GetGeom();
        NekDouble area = 0.0;
        Nektar::SpatialDomains::GeomFactorsSharedPtr geomFactor =
            geom->GetGeomFactors();
        // Array<OneD,NekDouble> jacAtQuads =
        // geomFactor->GetJac(geom->GetPointsKeys());
        Array<OneD, NekDouble> jacAtQuads =
            geomFactor->GetJac(m_expansions[0]->GetExp(eID)->GetPointsKeys());
        if (jacAtQuads.num_elements() == 1)
        {
            area = jacAtQuads[0];
        }
        else
        {
            area = jacAtQuads[0];
            // assert( false && "Not a regular element");
        }
        return area;
    };

    virtual NekDouble v_GetMeshLargestEdgeLength()
    {

        //			assert(false && "The subclasses routine should
        // be called."); 			return 0.0;
        NekDouble max = 0.0;
        for (int eid = 0; eid < m_expansions[0]->GetExpSize(); eid++)
        {
            SpatialDomains::GeometrySharedPtr geomSPtr =
                m_expansions[0]->GetExp(eid)->GetGeom();
            for (int edgeid = 0; edgeid < geomSPtr->GetNumEdges(); edgeid++)
            {
                NekDouble edgeLength =
                    m_segMap[geomSPtr->GetEid(edgeid)]->GetVertex(0)->dist(
                        *(m_segMap[geomSPtr->GetEid(edgeid)]->GetVertex(1)));
                max = std::max(max, edgeLength);
            }
        }
        return max;
    };

    virtual NekDouble v_GetLargestEdgeLength(const int eid)
    {
        //			assert(false && "The subclasses routine should
        // be called."); 			return 0.0;
        NekDouble max = 0.0;
        SpatialDomains::GeometrySharedPtr geomSPtr =
            m_expansions[0]->GetExp(eid)->GetGeom();
        for (int edgeid = 0; edgeid < geomSPtr->GetNumEdges(); edgeid++)
        {
            NekDouble edgeLength =
                m_segMap[geomSPtr->GetEid(edgeid)]->GetVertex(0)->dist(
                    *(m_segMap[geomSPtr->GetEid(edgeid)]->GetVertex(1)));
            max = std::max(max, edgeLength);
        }
        return max;
    };

    virtual bool v_CalculateDynamicScaling()
    {
        // 1. Algorithm 2.
        std::map<int, NekDouble> Sigma_m_a;
        std::map<int, NekDouble> Sigma_a;
        for (int eid = 0; eid < m_expansions[0]->GetExpSize(); eid++)
        {
            SpatialDomains::GeometrySharedPtr geomSPtr =
                m_expansions[0]->GetExp(eid)->GetGeom();
            NekDouble m = GetLargestEdgeLength(eid);
            NekDouble a = GetJacobian(eid);
            for (int vid = 0; vid < geomSPtr->GetNumVerts(); vid++)
            {
                int Vid = geomSPtr->GetVid(vid);
                if (Sigma_m_a.find(Vid) == Sigma_m_a.end())
                {
                    Sigma_m_a.insert(std::make_pair(Vid, m * a));
                    Sigma_a.insert(std::make_pair(Vid, a));
                }
                else
                {
                    Sigma_m_a.find(Vid)->second += m * a;
                    Sigma_a.find(Vid)->second += a;
                }
            }
        }

        // Scaling at vertices
        for (std::map<int, NekDouble>::iterator it = Sigma_m_a.begin();
             it != Sigma_m_a.end(); it++)
        {
            NekDouble totalArea = Sigma_a.find(it->first)->second;
            m_dynVertScaling.insert(
                std::make_pair(it->first, it->second / totalArea));
            // Currently in Debug mode
            // m_dynVertScaling.insert(std::make_pair(it->first,it->first));
        }
        return true;
    };

    virtual NekDouble v_GetDynamicScaling(Array<OneD, NekDouble> glCoord,
                                          int eid = -1, NekDouble mu = 1.0)
    {
        boost::ignore_unused(glCoord, eid, mu);
        NEKERROR(ErrorUtil::efatal, "v_GetDynamicScaling");
        return -1.0;
    };
    virtual bool v_WhatIsTRange(const NekDouble ptsX, const NekDouble ptsY,
                                const NekDouble ptsZ,
                                const Array<OneD, NekDouble> &direction,
                                NekDouble &tmin, NekDouble &tmax, int &num)
    {
        boost::ignore_unused(ptsX, ptsY, ptsZ, direction, tmin, tmax, num);
        NEKERROR(ErrorUtil::efatal, "v_WhatIsTRange");
        return false;
    };

    virtual bool v_GetMTScalingOfGIDs(vector<int> &t_GIDs,
                                      Array<OneD, NekDouble> &direction,
                                      vector<NekDouble> &scalings)
    {
        boost::ignore_unused(t_GIDs, direction, scalings);
        NEKERROR(ErrorUtil::efatal, "v_GetMTScalingOfGIDs");
        return false;
    };
};

} // namespace LSIAC
} // namespace Nektar
