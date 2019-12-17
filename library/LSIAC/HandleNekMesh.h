#pragma once
#include "HandleMesh.h"
#include "MetricTensor.h"
#include <FieldUtils/Field.hpp>
#include <FieldUtils/FieldUtilsDeclspec.h>
#include <SpatialDomains/Geometry.h>
#include <memory>
/// Class deals mesh queries of nektarMeshes.
/** This class in particular can load an Save only NektarMeshes.
        This class also help in selecting elements in the given neighbourhood.
   etc.
*/

namespace Nektar
{
namespace LSIAC
{

// Forward declartion
class MetricTensor;

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
    bool LoadMesh(string var)
    {
        return v_LoadMesh(var);
    }
    bool LoadData(string FileName, vector<string> &variables)
    {
        return v_LoadData(FileName, variables);
    }

    void LoadExpListIntoRTree()
    {
        v_LoadExpListIntoRTree();
    }

    int GetExpansionIndexUsingRTree(const Array<OneD, NekDouble> &point) const
    {
        return v_GetExpansionIndexUsingRTree(point);
    }

    //		bool SaveData(string FileName,
    // vector<nektar::MultiRegions::ExpList> expList);
    //! This function gives break points for one particular point in the mesh.
    /*! This function is optimized to work faster when arbitrary points need to
       be evaluated. Please do not use this function in a for loop to process a
       bunch of points. It can be really slow.
     */
    bool GetFilterOverlapElemIds(const NekDouble xcen_offset,
                                 const NekDouble ycen_offset,
                                 const NekDouble zcen_offset,
                                 const NekDouble *direction,
                                 const NekDouble t_offset_min,
                                 const NekDouble t_offset_max,
                                 vector<int> &Elid_list);
    //! This function gives break points for list of points in the mesh, which
    //! are close to each other.
    /*! This function is optimized to work faster when large number of points
       need to be evaluated for finding break points. This function tries to
       take advantage that there are points close to each other which needs to
       be evaluated. It can help reducing the computation time.
     */
    bool GetFilterOverlapElemIds(const Array<OneD, NekDouble> xcen_offset,
                                 const Array<OneD, NekDouble> ycen_offset,
                                 const Array<OneD, NekDouble> zcen_offset,
                                 const NekDouble *direction,
                                 const NekDouble t_offset_min,
                                 const NekDouble t_offset_max,
                                 vector<vector<int>> &Elid_list);
    //! To evaluate values of variables at different points on the mesh.
    bool EvaluateAt(const Array<OneD, NekDouble> &xPos,
                    const Array<OneD, NekDouble> &yPos,
                    const Array<OneD, NekDouble> &zPos, const int gID,
                    const int eID, Array<OneD, NekDouble> &values,
                    int varNum = 0)
    {
        return v_EvaluateAt(xPos, yPos, zPos, gID, eID, values, varNum);
    }

    //! To evaluate values of variables at a single point in the mesh.
    bool EvaluateAt(const NekDouble xPos, const NekDouble yPos,
                    const NekDouble zPos, int gID, int eID, NekDouble &value,
                    int varNum = 0)
    {
        return v_EvaluateAt(xPos, yPos, zPos, gID, eID, value, varNum);
    }

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

    //! This function gives break points for one particular point in the mesh.
    /*! This function is optimized to work faster when arbitrary points need to
       be evaluated. Please do not use this function in a for loop to process a
       bunch of points. It can be really slow.
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
    //! This function gives break points for list of points in the mesh, which
    //! are close to each other.
    /*! This function is optimized to work faster when large number of points
       need to be evaluated for finding break points. This function tries to
       take advantage that there are points close to each other which needs to
       be evaluated. It can help reducing the computation time.
     */
    bool GetBreakPts(const Array<OneD, NekDouble> xcen_offset,
                     const Array<OneD, NekDouble> ycen_offset,
                     const Array<OneD, NekDouble> zcen_offset,
                     const NekDouble *direction, const NekDouble t_offset_min,
                     const NekDouble t_offset_max,
                     vector<Array<OneD, NekDouble>> &xPos,
                     vector<Array<OneD, NekDouble>> &yPos,
                     vector<Array<OneD, NekDouble>> &zPos);

    bool CanTRangebeApplied(const NekDouble ptsX, const NekDouble ptsY,
                            const NekDouble ptsZ, const NekDouble scaling,
                            const NekDouble tmin, const NekDouble tmax,
                            NekDouble &tminUpdate, NekDouble &tmaxUpdate)
    {
        return v_CanTRangebeApplied(ptsX, ptsY, ptsZ, scaling, tmin, tmax,
                                    tminUpdate, tmaxUpdate);
    };

    //! This function determines if Symmetric kernel can be applied.
    /*! Returns true if Symmetric Kernel can be applied.
     *  Returns false if Symmetric Kernek cannot be applied.
     *  If returned false. Meshshift will contain shift requried by Meshknots
     *  Note: Meshknots are Kernelknots*meshScaling.
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

    bool GetMTScalingOfGIDs(vector<int> &t_GIDs,
                            Array<OneD, NekDouble> &direction,
                            vector<NekDouble> &scalings)
    {
        return v_GetMTScalingOfGIDs(t_GIDs, direction, scalings);
    };

    bool GetListOfGIDs(const NekDouble xPos, const NekDouble yPos,
                       const NekDouble zPos,
                       const Array<OneD, NekDouble> &direction,
                       const vector<NekDouble> t_breaks, vector<int> &t_GIDs,
                       vector<int> &t_EIDs) const
    {
        return v_GetListOfGIDs(xPos, yPos, zPos, direction, t_breaks, t_GIDs,
                               t_EIDs);
    };

    NekDouble GetElLargestEdgeSize(const NekDouble ptsx,
                                   const NekDouble ptsy = 0.0,
                                   const NekDouble ptsz = 0.0, int Elid = -1)
    {
        return v_GetElLargestEdgeSize(ptsx, ptsy, ptsz, Elid);
    };

    NekDouble GetJacobian(const int eID)
    {
        return v_GetJacobian(eID);
    };

    NekDouble GetLargestEdgeLength(const int eID)
    {
        return v_GetLargestEdgeLength(eID);
    };

    NekDouble GetMeshLargestEdgeLength()
    {
        return v_GetMeshLargestEdgeLength();
    };
    bool CalculateDynamicScaling()
    {
        return v_CalculateDynamicScaling();
    };

    NekDouble GetDynamicScaling(Array<OneD, NekDouble> &glCoord, int eid = -1,
                                NekDouble mu = 1.0)
    {
        return v_GetDynamicScaling(glCoord, eid, mu);
    };

    bool GetKnotVec(const int degree, const Array<OneD, NekDouble> &coord,
                    const Array<OneD, NekDouble> &direction,
                    Array<OneD, NekDouble> &knotVec, NekDouble &shift)
    {
        return v_GetKnotVec(degree, coord, direction, knotVec, shift);
    }

    bool Get1DVec(vector<NekDouble> &coords)
    {
        return v_Get1DVec(coords);
    }

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
