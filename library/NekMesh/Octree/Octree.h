////////////////////////////////////////////////////////////////////////////////
//
//  File: Octree.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: octree object header
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_MESHUTILS_OCTREE_OCTREE
#define NEKTAR_MESHUTILS_OCTREE_OCTREE

#include "SourcePoint.hpp"
#include "Octant.h"
#include <NekMesh/MeshElements/Mesh.h>
#include <NekMesh/Module/Log.hpp>

#include <string>

namespace Nektar
{
namespace NekMesh
{

/**
 * @brief this struct defines a CAD curve object to be used for defining spatial
 * refinment where a fixed element edge length is set to any element on or
 * within a specified distance from it.
 *
 * @param curve the CAD curve object - currently only supported using the OCE /
 * OpenCASCADE backend engine.
 * @param R the distance from the CAD curve where a fixed element edge length is
 * to be applied.
 * @param delta the fixed elementedge length.
 */
struct curvesource
{
    CADCurveSharedPtr curve;
    NekDouble R, delta;
    curvesource(CADCurveSharedPtr c,
               NekDouble r,
               NekDouble d)
        : curve(c), R(r), delta(d)
    {
    }

    /**
     * @brief tests if a point is within a specified range, R, from the CAD
     * curve.
     *
     * @param p array with the x,y,z position of the point to be tested.
     * @return true is within the specified distance and false if not.
     */
    bool withinRange(Array<OneD, NekDouble> p)
    {
        if(curve->GetMinDistance(p) <= R)
        {
            return true;
        }
        else
            return false;
    }
};

/**
 * @brief this struct defines two points that create a line to be used for defining
 * spatial refinment where a fixed element edge length is set to any element on
 * or within a specified distance from it.
 *
 * @param x1, x2 arrays containing the x,y,z positions of the two points used to
 * define the line.
 * @param R the distance from the line where a fixed element edge length is to
 * be applied.
 * @param delta the fixed elementedge length.
 */
struct linesource
{
    Array<OneD, NekDouble> x1, x2;
    NekDouble R, delta;
    linesource(Array<OneD, NekDouble> p1,
               Array<OneD, NekDouble> p2,
               NekDouble r,
               NekDouble d)
        : x1(p1), x2(p2), R(r), delta(d)
    {
    }

    /**
     * @brief tests if a point is within a specified range, R, from the line
     *
     * @param p array with the x,y,z position of the point to be tested.
     * @return true is within the specified distance and false if not.
     */
    bool withinRange(Array<OneD, NekDouble> p)
    {
        Array<OneD, NekDouble> Le(3), Re(3), s(3);
        for (int i = 0; i < 3; i++)
        {
            Le[i] = p[i] - x1[i];
            Re[i] = p[i] - x2[i];
            s[i]  = x2[i] - x1[i];
        }

        // check distances to endpoints
        if (Le[0]*Le[0] + Le[1]*Le[1] + Le[2]*Le[2] < R * R)
        {
            return true;
        }
        if (Re[0]*Re[0] + Re[1]*Re[1] + Re[2]*Re[2] < R * R)
        {
            return true;
        }

        // Do an orthogonal projection of p onto the line segment
        Array<OneD, NekDouble> dev(3);
        // (p-x1) \times (p-x2)
        dev[0] = Le[1] * Re[2] - Re[1] * Le[2];
        dev[1] = Le[2] * Re[0] - Re[2] * Le[0];
        dev[2] = Le[0] * Re[1] - Re[0] * Le[1];

        NekDouble dist =
            sqrt(dev[0] * dev[0] + dev[1] * dev[1] + dev[2] * dev[2]) / Length();

        NekDouble t = -1.0 * ((x1[0] - p[0]) * s[0] + (x1[1] - p[1]) * s[1] +
                              (x1[2] - p[2]) * s[2]) / Length() / Length();

        if (dist < R && !(t > 1) && !(t < 0))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    /**
     * @brief returns the length of the line defining the linesource
     */
    NekDouble Length()
    {
        return sqrt((x1[0] - x2[0]) * (x1[0] - x2[0]) +
                    (x1[1] - x2[1]) * (x1[1] - x2[1]) +
                    (x1[2] - x2[2]) * (x1[2] - x2[2]));
    }
};

/**
 * @brief class for octree
 *
 * This class contains the routines to generate and query a automatically
 * generated set of mesh spacing parameters based on the CAD
 */
class Octree
{
public:

    Octree(MeshSharedPtr m, Logger log) : m_mesh(m), m_log(log)
    {
        m_log.SetPrefix("Octree");
    }

    /**
     * @brief builds the octree based on curvature sampling and user defined
     * spacing
     */
    void Process();

    /**
     * @brief once constructed queryies the octree based on x,y,z location
     * to get a mesh spacing
     *
     * @param loc array of x,y,z
     * @return mesh spacing parameter
     */
    NekDouble Query(Array<OneD, NekDouble> loc);

    /**
     * @brief returns the miminum spacing in the octree (for meshing purposes)
     *
     * @return miminum delta in octree
     */
    NekDouble GetMinDelta();

    /**
     * @brief sets the parameters used for curvature sampling
     *
     * @param min minimum spacing to be found in the mesh
     * @param max maximum spacing to be found in the mesh
     * @param eps curvature sensivity relating radius of curvature to spacing
     */
    void SetParameters(NekDouble &min, NekDouble &max, NekDouble &ep)
    {
        m_minDelta = min;
        m_maxDelta = max;
        m_eps = ep;
    }

    /**
     * @brief populates the mesh m with a invalid hexahedral mesh based on the
     *        octree, used for visualisation
     * @param nm name of the mesh file to be made
     */
    void WriteOctree(std::string nm);

    /**
     * @brief informs the octree there is a user defined line spacing file
     *
     * @param nm name of the user defined spacing file
     */
    void Refinement(std::string nm)
    {
        m_refinement = nm;
    }

    /**
     * @brief informs the octree there is a user defined curve spacing file
     *
     * @param nm name of the user defined spacing file
     */
    void CurveRefinement(std::string nm)
    {
        m_curverefinement = nm;
    }

private:

    /**
     * @brief Smooths specification over all octants to a gradation criteria
     */
    void SmoothAllOctants();

    /**
     * @brief gets an optimum number of curvature sampling points and
     * calculates the curavture at these points
     */
    void CompileSourcePointList();

    /**
     * @brief Function which initiates and controls the subdivision process
     */
    void SubDivide();

    /**
     * @brief Smooths specification over the surface encompasing octants to a
     *        gradation criteria
     */
    void SmoothSurfaceOctants();

    /**
     * @brief takes the mesh specification from surface octants and
     *        progates that through the domain so all octants have a
     * specification
     *        using gradiation crieteria
     */
    void PropagateDomain();

    /**
     * @brief estimates the number of elements to be created in the mesh
     */
    int CountElemt();

    /**
     * @brief Calculates the difference in delta divided by the difference
     *        in location between two octants i and j
     */
    NekDouble ddx(OctantSharedPtr i, OctantSharedPtr j);

    /**
     * @brief Looks over all leaf octants and checks that their neigbour
     *        assigments are valid
     */
    bool VerifyNeigbours();

    /// minimum delta in the octree
    NekDouble m_minDelta;
    /// maximum delta in the octree
    NekDouble m_maxDelta;
    /// curavture sensivity paramter
    NekDouble m_eps;
    /// x,y,z location of the center of the octree
    Array<OneD, NekDouble> m_centroid;
    /// physical size of the octree
    NekDouble m_dim;
    /// list of source points
    std::vector<SPBaseSharedPtr> m_SPList;
    /// list of leaf octants
    std::vector<OctantSharedPtr> m_octants;
    /// master octant for searching
    OctantSharedPtr m_masteroct;
    /// number of octants made, used for id index
    int m_numoct;
    /// Mesh object
    MeshSharedPtr m_mesh;
    /// Logger
    Logger m_log;

    std::string m_refinement;
    std::string m_curverefinement;
    std::vector<linesource> m_lsources;
    std::vector<curvesource> m_csources;
};
typedef std::shared_ptr<Octree> OctreeSharedPtr;

}
}

#endif
