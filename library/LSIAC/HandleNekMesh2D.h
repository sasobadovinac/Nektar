#pragma once
#include "HandleNekMesh.h"
#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace Boostg  = boost::geometry;
namespace Boostgi = boost::geometry::index;

typedef Boostg::model::point<NekDouble, 3, Boostg::cs::cartesian> RPoint;
typedef Boostg::model::box<RPoint> RBox;
typedef std::tuple<RBox, unsigned, unsigned> RValue;
/// typedef Boostgi::rtree<RValue, Boostgi::rstar<4> > RTree;
typedef Boostgi::rtree<RValue, Boostgi::quadratic<10, 3>> RTree;

namespace Nektar
{
namespace LSIAC
{

/// Handles Nektar 2D meshes.
class HandleNekMesh2D : public HandleNekMesh
{

private:
    bool m_useRTree;

protected:
    RTree m_rtree;

public:
    HandleNekMesh2D(FieldUtils::FieldSharedPtr fldSharedPtr)
        : HandleNekMesh(fldSharedPtr), m_useRTree(false)
    {
    }

    HandleNekMesh2D(LibUtilities::SessionReaderSharedPtr sessionPtr)
        : HandleNekMesh(sessionPtr), m_useRTree(false)
    {
    }

    void IntersectWithBoxUsingRTree(const NekDouble minCornerX,
                                    const NekDouble minCornerY,
                                    const NekDouble minCornerZ,
                                    const NekDouble maxCornerX,
                                    const NekDouble maxCornerY,
                                    const NekDouble maxCornerZ,
                                    vector<int> &elIds, vector<int> &glIds);

protected:
    virtual int v_GetExpansionIndexUsingRTree(
        const Array<OneD, NekDouble> &point) const;

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
    // functions used by RTree
    virtual void v_LoadExpListIntoRTree();

private:
    vector<NekDouble> cross_Math(const vector<NekDouble> &r,
                                 const vector<NekDouble> &s);

    vector<NekDouble> sub_Math(vector<NekDouble> &p2, vector<NekDouble> &p1);

    NekDouble dot_Math(vector<NekDouble> &p, vector<NekDouble> &q);

    NekDouble norm2_Math(vector<NekDouble> p);

    bool intersect(vector<NekDouble> &p1, vector<NekDouble> &p2,
                   vector<NekDouble> &q1, vector<NekDouble> &q2,
                   vector<NekDouble> &i1, vector<NekDouble> &i2);

    void IntersectWithEdges(const SpatialDomains::SegGeomMap &segMap,
                            const SpatialDomains::PointGeomMap &pointMap,
                            const Array<OneD, NekDouble> &dir,
                            const Array<OneD, NekDouble> &point,
                            const NekDouble t1, const NekDouble t2,
                            vector<NekDouble> &tvalT);

    void IntersectWithFewEdges(const vector<int> EdgeIds,
                               const SpatialDomains::PointGeomMap &pointMap,
                               const Array<OneD, NekDouble> &dir,
                               const Array<OneD, NekDouble> &point,
                               const NekDouble t1, const NekDouble t2,
                               vector<NekDouble> &tvalT);

    void IntersectWithEdgesUsingRTree(
        const SpatialDomains::SegGeomMap &segMap,
        const SpatialDomains::PointGeomMap &pointMap,
        const Array<OneD, NekDouble> &dir, const Array<OneD, NekDouble> &point,
        const NekDouble t1, const NekDouble t2, vector<NekDouble> &tvalT);

    void FindElementIDForLineSegs(
        const vector<NekDouble> &tvalT, const Array<OneD, NekDouble> &point,
        const Array<OneD, NekDouble> &dir,
        const SpatialDomains::MeshGraphSharedPtr mesh_graph, vector<int> &EIDs);

    virtual NekDouble v_GetElLargestEdgeSize(const NekDouble Ptsx = 0.0,
                                             const NekDouble Ptsy = 0.0,
                                             const NekDouble Ptsz = 0.0,
                                             int Elid             = -1);

    /*
                    // functions used by RTree
            void BoundingBoxOfLineSeg( const Array<OneD,NekDouble> &dir,
                const Array<OneD,NekDouble> &pt, const NekDouble t1, const
       NekDouble t2, RBox &b );


            void IntersectWithFacesUsingRTree( const Array<OneD,NekDouble> &dir,
       const Array<OneD,NekDouble> &point, const NekDouble t1, const NekDouble
       t2, vector<NekDouble> &tvalT );



            void GetBoundingOfElement(SpatialDomains::GeometrySharedPtr elGeom,
       RBox &b);

            NekDouble getMax( NekDouble a, NekDouble b)
            {
                if (a<b)
                    return b;
                else
                    return a;
            }
            NekDouble getMin( NekDouble a, NekDouble b)
            {
                if (a<b)
                    return a;
                else
                    return b;
            }
    */

    virtual void GetBoundingOfElement(SpatialDomains::GeometrySharedPtr elGeom,
                                      RBox &b);

    NekDouble getMax(NekDouble a, NekDouble b)
    {
        if (a < b)
            return b;
        else
            return a;
    }
    NekDouble getMin(NekDouble a, NekDouble b)
    {
        if (a < b)
            return a;
        else
            return b;
    }

    // virtual bool v_CalculateDynamicScaling();
    // virtual NekDouble v_GetJacobian(const int eID);
    virtual NekDouble v_GetDynamicScaling(Array<OneD, NekDouble> glCoord,
                                          int eid = -1, NekDouble mu = 1.0);
    virtual bool v_WhatIsTRange(const NekDouble PtsX, const NekDouble PtsY,
                                const NekDouble PtsZ,
                                const Array<OneD, NekDouble> &direction,
                                NekDouble &tmin, NekDouble &tmax, int &num);
    virtual bool v_InitializeMetricTensor();
    virtual bool v_GetMTScalingOfGIDs(vector<int> &t_GIDs,
                                      Array<OneD, NekDouble> &direction,
                                      vector<NekDouble> &scalings);
};
} // namespace LSIAC
} // namespace Nektar
