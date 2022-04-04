////////////////////////////////////////////////////////////////////////////////
//
//  File: Zones.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
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
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/Movement/Zones.h>
#include <tinyxml.h>

namespace Nektar
{
namespace SpatialDomains
{

ZoneBase::ZoneBase(MovementType type, int indx, CompositeMap domain,
                   int coordDim)
    : m_type(type), m_id(indx), m_domain(domain), m_coordDim(coordDim)
{
    // Fill elements and element IDs from domain
    int shapeDim = coordDim; // @TODO: We are using coordim when this should be shapedim, find easy way in MPI to determine shapedim
    m_constituentElements = Array<OneD,std::set<GeometrySharedPtr>>(shapeDim);

    for (auto &comp : domain)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            m_elementIds.emplace_back(geom->GetGlobalID());
            m_elements.emplace_back(geom);

            // Fill consituent elements (i.e. faces and edges)
            switch (shapeDim)
            {
                case 3:
                    for(int i = 0; i < geom->GetNumFaces(); ++i)
                    {
                        m_constituentElements[2].insert(geom->GetFace(i));
                    }
                    /* fall through */
                case 2:
                    for(int i = 0; i < geom->GetNumEdges(); ++i)
                    {
                        m_constituentElements[1].insert(geom->GetEdge(i));
                    }
                    /* fall through */
                case 1:
                    m_constituentElements[0].insert(geom);
                    /* fall through */
                default:
                    break;
            }
        }
    }
}

ZoneRotate::ZoneRotate(int id, const CompositeMap &domain, const int coordDim,
                       const NekPoint<NekDouble> &origin, const DNekVec &axis,
                       const LibUtilities::EquationSharedPtr &angularVelEqn)
    : ZoneBase(MovementType::eRotate, id, domain, coordDim), m_origin(origin),
      m_axis(axis), m_angularVelEqn(angularVelEqn)
{
    std::set<int> seenVerts, seenEdges, seenFaces;
    for (auto &comp : m_domain)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            for (int i = 0; i < geom->GetNumVerts(); ++i)
            {
                PointGeomSharedPtr vert = geom->GetVertex(i);

                if (seenVerts.find(vert->GetGlobalID()) != seenVerts.end())
                {
                    continue;
                }

                seenVerts.insert(vert->GetGlobalID());
                m_rotateVerts.emplace_back(vert);
            }

            for (int i = 0; i < geom->GetNumEdges(); ++i)
            {
                SegGeomSharedPtr edge =
                    std::static_pointer_cast<SegGeom>(geom->GetEdge(i));

                if (seenEdges.find(edge->GetGlobalID()) != seenEdges.end())
                {
                    continue;
                }

                seenEdges.insert(edge->GetGlobalID());

                CurveSharedPtr curve = edge->GetCurve();
                if (!curve)
                {
                    continue;
                }

                m_rotateCurves.emplace_back(curve);
            }

            for (int i = 0; i < geom->GetNumFaces(); ++i)
            {
                Geometry2DSharedPtr face =
                    std::static_pointer_cast<Geometry2D>(geom->GetFace(i));

                if (seenFaces.find(face->GetGlobalID()) != seenFaces.end())
                {
                    continue;
                }

                seenFaces.insert(face->GetGlobalID());

                CurveSharedPtr curve = face->GetCurve();
                if (!curve)
                {
                    continue;
                }

                m_rotateCurves.emplace_back(curve);
            }
        }
    }

    // Copy points so we know original positions.
    for (auto &pt : m_rotateVerts)
    {
        m_origPosition.emplace_back(*pt);
    }

    for (auto &curve : m_rotateCurves)
    {
        for (auto &pt : curve->m_points)
        {
            m_origPosition.emplace_back(*pt);
        }
    }

    // Construct rotation matrix
    m_W(0, 1) = -m_axis[2];
    m_W(0, 2) = m_axis[1];
    m_W(1, 0) = m_axis[2];
    m_W(1, 2) = -m_axis[0];
    m_W(2, 0) = -m_axis[1];
    m_W(2, 1) = m_axis[0];

    m_W2 = m_W * m_W;

    // std::cout << "W matrix: " << std::endl;
    // std::cout << m_W << std::endl << std::endl;
    // std::cout << "W^2 matrix: " << std::endl;
    // std::cout << m_W2 << std::endl << std::endl;
    // std::cout << "Origin: " << std::endl;
    // std::cout << m_origin << std::endl << std::endl;
}

ZoneTranslate::ZoneTranslate(int id, const CompositeMap &domain,
                             const int coordDim,
                             const std::vector<NekDouble> &velocity)
    : ZoneBase(MovementType::eTranslate, id, domain, coordDim),
      m_velocity(velocity)
{
    std::set<int> seenVerts, seenEdges;
    for (auto &comp : m_domain)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            for (int i = 0; i < geom->GetNumVerts(); ++i)
            {
                PointGeomSharedPtr vert = geom->GetVertex(i);

                if (seenVerts.find(vert->GetGlobalID()) != seenVerts.end())
                {
                    continue;
                }

                seenVerts.insert(vert->GetGlobalID());
                m_slideVerts.emplace_back(vert);
            }

            for (int i = 0; i < geom->GetNumEdges(); ++i)
            {
                SegGeomSharedPtr edge =
                    std::static_pointer_cast<SegGeom>(geom->GetEdge(i));

                if (seenEdges.find(edge->GetGlobalID()) != seenEdges.end())
                {
                    continue;
                }

                seenEdges.insert(edge->GetGlobalID());

                CurveSharedPtr curve = edge->GetCurve();
                if (!curve)
                {
                    continue;
                }

                m_slideCurves.emplace_back(curve);
            }
        }
    }

    // Copy points so we know original positions.
    for (auto &pt : m_slideVerts)
    {
        m_origPosition.emplace_back(*pt);
    }

    for (auto &curve : m_slideCurves)
    {
        for (auto &pt : curve->m_points)
        {
            m_origPosition.emplace_back(*pt);
        }
    }
}

ZonePrescribe::ZonePrescribe(int id, const CompositeMap &domain,
                             const int coordDim,
                             LibUtilities::EquationSharedPtr xDeform,
                             LibUtilities::EquationSharedPtr yDeform,
                             LibUtilities::EquationSharedPtr zDeform)
    : ZoneBase(MovementType::ePrescribe, id, domain, coordDim),
      m_xDeform(xDeform), m_yDeform(yDeform), m_zDeform(zDeform)
{
    std::set<int> seenVerts, seenEdges;
    for (auto &comp : m_domain)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            for (int i = 0; i < geom->GetNumVerts(); ++i)
            {
                PointGeomSharedPtr vert = geom->GetVertex(i);

                if (seenVerts.find(vert->GetGlobalID()) != seenVerts.end())
                {
                    continue;
                }

                seenVerts.insert(vert->GetGlobalID());
                m_interiorVerts.emplace_back(vert);
            }

            for (int i = 0; i < geom->GetNumEdges(); ++i)
            {
                SegGeomSharedPtr edge =
                    std::static_pointer_cast<SegGeom>(geom->GetEdge(i));

                if (seenEdges.find(edge->GetGlobalID()) != seenEdges.end())
                {
                    continue;
                }

                seenEdges.insert(edge->GetGlobalID());

                CurveSharedPtr curve = edge->GetCurve();
                if (!curve)
                {
                    continue;
                }

                for (auto &point : curve->m_points)
                {
                    if (seenVerts.find(point->GetGlobalID()) != seenVerts.end())
                    {
                        m_interiorVerts.emplace_back(point);
                    }
                }
            }
        }
    }

    // Copy points so we know original positions.
    for (auto &pt : m_interiorVerts)
    {
        m_origPosition.emplace_back(*pt);
    }
}

NekDouble ZoneRotate::GetAngularVel(NekDouble &time) const
{
    int spinUpTime = 10;
    if(time < spinUpTime)
    {
        return m_angularVelEqn->Evaluate(0,0,0,time) * spinUpTime / 10;
    }
    else
    {
        return m_angularVelEqn->Evaluate(0,0,0,time);
    }
}

// Calculate new location of points using Rodrigues formula
bool ZoneRotate::v_Move(NekDouble time)
{
    NekDouble angle = GetAngularVel(time) * time;
    // TODO: For none constant angular velocity this doesn't work ^^
    // @TODO: I need to take into account the total angle, summing timesteps works here
    // @TODO: but then it doesn't work for FieldConvert where only the checkpoint time is known
    // TODO: I want to integrate m_angularVelEqn up to current time
    //NekDouble period = 20.0;
    //NekDouble angle = ceil(1.0/period*time - 1.0) * time -  period / 2.0 * (ceil(1.0/period*time - 1.0) *  ceil(1.0/period*time));

    // Identity matrix
    DNekMat rot(3, 3, 0.0);
    rot(0, 0) = 1.0;
    rot(1, 1) = 1.0;
    rot(2, 2) = 1.0;

    rot = rot + sin(angle) * m_W + (1 - cos(angle)) * m_W2;

    //std::cout << "Rotation matrix: " << std::endl;
    //std::cout << rot << std::endl;
    int cnt = 0;
    for (auto &vert : m_rotateVerts)
    {
        NekPoint<NekDouble> pnt = m_origPosition[cnt] - m_origin;
        DNekVec pntVec          = {pnt[0], pnt[1], pnt[2]};

        DNekVec newLoc = rot * pntVec;

        vert->UpdatePosition(newLoc(0) + m_origin[0], newLoc(1) + m_origin[1],
                             newLoc(2) + m_origin[2]);
        cnt++;
    }

    for (auto &curve : m_rotateCurves)
    {
        for (auto &vert : curve->m_points)
        {
            NekPoint<NekDouble> pnt = m_origPosition[cnt] - m_origin;
            DNekVec pntVec          = {pnt[0], pnt[1], pnt[2]};

            DNekVec newLoc = rot * pntVec;

            vert->UpdatePosition(newLoc(0) + m_origin[0],
                                 newLoc(1) + m_origin[1],
                                 newLoc(2) + m_origin[2]);
            cnt++;
        }
    }

    // Clear bounding boxes (these will be regenerated next time GetBoundingBox is called)
    for (auto &el : m_elements)
    {
        el->ClearBoundingBox();

        int nfaces = el->GetNumFaces();
        for (int i = 0; i < nfaces; ++i)
        {
            el->GetFace(i)->ClearBoundingBox();
        }

        int nedges = el->GetNumEdges();
        for (int i = 0; i < nedges; ++i)
        {
            el->GetEdge(i)->ClearBoundingBox();
        }
    }

    return true;
}

bool ZoneTranslate::v_Move(NekDouble timeStep)
{
    Array<OneD, NekDouble> dist(3, 0.0);
    for (int i = 0; i < m_coordDim; ++i)
    {
        dist[i] = m_velocity[i] * timeStep;
    }

    int cnt = 0;
    for (auto &vert : m_slideVerts)
    {
        Array<OneD, NekDouble> newLoc(3, 0.0);
        auto pnt = m_origPosition[cnt];

        for (int i = 0; i < m_coordDim; ++i)
        {
            newLoc[i] = pnt(i) + dist[i];
        }

        vert->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
        cnt++;
    }

    for (auto &curve : m_slideCurves)
    {
        for (auto &vert : curve->m_points)
        {
            Array<OneD, NekDouble> newLoc(3, 0.0);
            auto pnt = m_origPosition[cnt];

            for (int i = 0; i < m_coordDim; ++i)
            {
                newLoc[i] = pnt(i) + dist[i];
            }

            vert->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
            cnt++;
        }
    }

    // Clear bounding boxes (these will be regenerated next time GetBoundingBox is called)
    for (auto &el : m_elements)
    {
        el->ClearBoundingBox();

        int nfaces = el->GetNumFaces();
        for (int i = 0; i < nfaces; ++i)
        {
            el->GetFace(i)->ClearBoundingBox();
        }

        int nedges = el->GetNumEdges();
        for (int i = 0; i < nedges; ++i)
        {
            el->GetEdge(i)->ClearBoundingBox();
        }
    }

    return true;
}

bool ZoneFixed::v_Move(NekDouble time)
{
    boost::ignore_unused(time);
    return false;
}

bool ZonePrescribe::v_Move(NekDouble time)
{
    int cnt = 0;
    for (auto &vert : m_interiorVerts)
    {
        auto pnt = m_origPosition[cnt++];

        Array<OneD, NekDouble> coords(3, 0.0);
        vert->GetCoords(coords);

        Array<OneD, NekDouble> newLoc(3, 0.0);
        newLoc[0] = m_xDeform->Evaluate(coords[0], coords[1], coords[2], time) + pnt(0);
        newLoc[1] = m_yDeform->Evaluate(coords[0], coords[1], coords[2], time) + pnt(1);
        newLoc[2] = m_zDeform->Evaluate(coords[0], coords[1], coords[2], time) + pnt(2);

        vert->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
    }

    // Clear bounding boxes (these will be regenerated next time GetBoundingBox is called)
    for (auto &el : m_elements)
    {
        el->ClearBoundingBox();

        int nfaces = el->GetNumFaces();
        for (int i = 0; i < nfaces; ++i)
        {
            el->GetFace(i)->ClearBoundingBox();
        }

        int nedges = el->GetNumEdges();
        for (int i = 0; i < nedges; ++i)
        {
            el->GetEdge(i)->ClearBoundingBox();
        }
    }

    return true;
}

}
}