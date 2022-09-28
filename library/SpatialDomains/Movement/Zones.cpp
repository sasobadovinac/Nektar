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
//  Description: Zones used in the non-conformal interfaces
//               and ALE implementations
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
    for (auto &comp : domain)
    {
        for (auto &geom : comp.second->m_geomVec)
        {
            m_elements.emplace_back(geom);

            // Fill constituent elements (i.e. faces and edges)
            switch (geom->GetShapeDim())
            {
                case 3:
                    for (int i = 0; i < geom->GetNumFaces(); ++i)
                    {
                        m_constituentElements[2].insert(geom->GetFace(i));
                    }
                    /* fall through */
                case 2:
                    for (int i = 0; i < geom->GetNumEdges(); ++i)
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

    // The seenVerts/Edges/Faces keeps track of geometry so duplicates aren't
    // emplaced in to the storage vector
    std::unordered_set<int> seenVerts, seenEdges, seenFaces;

    // Fill verts/edges/faces vector storage that are to be moved each timestep
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
                m_verts.emplace_back(vert);
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

                m_curves.emplace_back(curve);
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

                m_curves.emplace_back(curve);
            }
        }
    }

    // Copy points so we know original positions.
    for (auto &pt : m_verts)
    {
        m_origVerts.emplace_back(*pt);
    }

    for (auto &curve : m_curves)
    {
        for (auto &pt : curve->m_points)
        {
            m_origVerts.emplace_back(*pt);
        }
    }
}

ZoneRotate::ZoneRotate(int id, const CompositeMap &domain, const int coordDim,
                       const NekPoint<NekDouble> &origin, const DNekVec &axis,
                       const LibUtilities::EquationSharedPtr &angularVelEqn)
    : ZoneBase(MovementType::eRotate, id, domain, coordDim), m_origin(origin),
      m_axis(axis), m_angularVelEqn(angularVelEqn)
{
    // Construct rotation matrix
    m_W(0, 1) = -m_axis[2];
    m_W(0, 2) = m_axis[1];
    m_W(1, 0) = m_axis[2];
    m_W(1, 2) = -m_axis[0];
    m_W(2, 0) = -m_axis[1];
    m_W(2, 1) = m_axis[0];

    m_W2 = m_W * m_W;
}

void ZoneBase::ClearBoundingBoxes()
{
    // Clear bboxes (these will be regenerated next time GetBoundingBox is
    // called)
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
}

NekDouble ZoneRotate::GetAngularVel(NekDouble &time) const
{
    NekDouble rampTime = 1;
    if (time < rampTime)
    {
        return m_angularVelEqn->Evaluate(0, 0, 0, time) * time/rampTime;
    }
    else
    {
        return m_angularVelEqn->Evaluate(0, 0, 0, time);
    }
}

// Calculate new location of points using Rodrigues formula
bool ZoneRotate::v_Move(NekDouble time)
{
    // Currently only valid for constant angular velocity
    NekDouble angle = GetAngularVel(time) * time;
    // TODO: For none constant angular velocity this doesn't work ^^
    // @TODO: I need to take into account the total angle, summing timesteps
    // works here
    // @TODO: but then it doesn't work for FieldConvert where only the
    // checkpoint time is known
    // TODO: I want to integrate m_angularVelEqn up to current time

    // Identity matrix
    DNekMat rot(3, 3, 0.0);
    rot(0, 0) = 1.0;
    rot(1, 1) = 1.0;
    rot(2, 2) = 1.0;

    // Rodrigues' rotation formula in matrix form
    rot = rot + sin(angle) * m_W + (1 - cos(angle)) * m_W2;

    int cnt = 0;
    for (auto &vert : m_verts)
    {
        NekPoint<NekDouble> pnt = m_origVerts[cnt] - m_origin;
        DNekVec pntVec          = {pnt[0], pnt[1], pnt[2]};

        DNekVec newLoc = rot * pntVec;

        vert->UpdatePosition(newLoc(0) + m_origin[0], newLoc(1) + m_origin[1],
                             newLoc(2) + m_origin[2]);
        cnt++;
    }

    for (auto &curve : m_curves)
    {
        for (auto &vert : curve->m_points)
        {
            NekPoint<NekDouble> pnt = m_origVerts[cnt] - m_origin;
            DNekVec pntVec          = {pnt[0], pnt[1], pnt[2]};

            DNekVec newLoc = rot * pntVec;

            vert->UpdatePosition(newLoc(0) + m_origin[0],
                                 newLoc(1) + m_origin[1],
                                 newLoc(2) + m_origin[2]);
            cnt++;
        }
    }

    ClearBoundingBoxes();

    return true;
}

std::vector<NekDouble> ZoneTranslate::GetVel(NekDouble &time) const
{
    std::vector<NekDouble> vel(m_coordDim);
    for (int i = 0; i < m_coordDim; ++i)
    {
        vel[i] = m_velocityEqns[i]->Evaluate(0, 0, 0, time);
    }

    return vel;
}

std::vector<NekDouble> ZoneTranslate::GetDisp(NekDouble &time) const
{
    std::vector<NekDouble> disp(m_coordDim);
    for (int i = 0; i < m_coordDim; ++i)
    {
        disp[i] = m_displacementEqns[i]->Evaluate(0, 0, 0, time);
    }

    return disp;
}

bool ZoneTranslate::v_Move(NekDouble time)
{
    auto disp = GetDisp(time);

    int cnt = 0;
    for (auto &vert : m_verts)
    {
        Array<OneD, NekDouble> newLoc(3, 0.0);
        auto pnt = m_origVerts[cnt];

        for (int i = 0; i < m_coordDim; ++i)
        {
            newLoc[i] = pnt(i) + disp[i];
        }

        vert->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
        cnt++;
    }

    for (auto &curve : m_curves)
    {
        for (auto &vert : curve->m_points)
        {
            Array<OneD, NekDouble> newLoc(3, 0.0);
            auto pnt = m_origVerts[cnt];

            for (int i = 0; i < m_coordDim; ++i)
            {
                newLoc[i] = pnt(i) + disp[i];
            }

            vert->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
            cnt++;
        }
    }

    ClearBoundingBoxes();

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
    for (auto &vert : m_verts)
    {
        auto pnt = m_origVerts[cnt++];

        Array<OneD, NekDouble> coords(3, 0.0);
        vert->GetCoords(coords);

        Array<OneD, NekDouble> newLoc(3, 0.0);
        /*newLoc[0] =
            m_xDeform->Evaluate(coords[0], coords[1], coords[2], time) + pnt(0);
        newLoc[1] =
            m_yDeform->Evaluate(coords[0], coords[1], coords[2], time) + pnt(1);
        newLoc[2] =
            m_zDeform->Evaluate(coords[0], coords[1], coords[2], time) + pnt(2);*/

        newLoc[0] = pnt(0) + 0.5 * sin(2 * M_PI * time / sqrt(50)) *
                                   sin(2 * M_PI * coords[0] / 20) *
                                   sin(2 * M_PI * coords[1] / 20);
        newLoc[1] = pnt(1) + 0.5 * sin(2 * M_PI * time / sqrt(50)) *
                                   sin(2 * M_PI * coords[0] / 20) *
                                   sin(2 * M_PI * coords[1] / 20);
        newLoc[2] = 0.0;
        vert->UpdatePosition(newLoc[0], newLoc[1], newLoc[2]);
    }

    ClearBoundingBoxes();

    return true;
}

} // namespace SpatialDomains
} // namespace Nektar
