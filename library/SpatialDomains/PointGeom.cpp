////////////////////////////////////////////////////////////////////////////////
//
//  File:  PointGeom.cpp
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
//  Description: Point geometry information
//
////////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include <StdRegions/StdRegions.hpp>
#include <SpatialDomains/PointGeom.h>
#include <SpatialDomains/SegGeom.h>

namespace Nektar
{
namespace SpatialDomains
{
PointGeom::PointGeom() : NekPoint<NekDouble>(0.0, 0.0, 0.0)
{
    m_shapeType = LibUtilities::ePoint;
    m_coordim = 0;
    m_globalID = 0;
}

PointGeom::PointGeom(
    const int coordim, const int vid, NekDouble x, NekDouble y, NekDouble z)
    : NekPoint<NekDouble>(x, y, z)
{
    m_shapeType = LibUtilities::ePoint;
    m_coordim = coordim;
    m_globalID = vid;

    (*this)(0) = x;
    (*this)(1) = y;
    (*this)(2) = z;
}

// copy constructor
PointGeom::PointGeom(const PointGeom &T) : NekPoint<NekDouble>(T)
{
    m_shapeType = T.m_shapeType;
    m_globalID = T.m_globalID;
    m_coordim = T.m_coordim;
}

PointGeom::~PointGeom()
{
}

void PointGeom::GetCoords(NekDouble &x, NekDouble &y, NekDouble &z)
{
    switch (m_coordim)
    {
        case 3:
            z = (*this)(2);
        case 2:
            y = (*this)(1);
        case 1:
            x = (*this)(0);
            break;
    }
}

void PointGeom::GetCoords(Array<OneD, NekDouble> &coords)
{
    switch (m_coordim)
    {
        case 3:
            coords[2] = (*this)(2);
        case 2:
            coords[1] = (*this)(1);
        case 1:
            coords[0] = (*this)(0);
            break;
    }
}

void PointGeom::UpdatePosition(NekDouble x, NekDouble y, NekDouble z)
{
    (*this)(0) = x;
    (*this)(1) = y;
    (*this)(2) = z;
}

// _this = a + b
void PointGeom::Add(PointGeom &a, PointGeom &b)
{
    (*this)(0) = a[0] + b[0];
    (*this)(1) = a[1] + b[1];
    (*this)(2) = a[2] + b[2];
    m_coordim = std::max(a.GetCoordim(), b.GetCoordim());
}

// _this = a + b
void PointGeom::Sub(PointGeom &a, PointGeom &b)
{
    (*this)(0) = a[0] - b[0];
    (*this)(1) = a[1] - b[1];
    (*this)(2) = a[2] - b[2];
    m_coordim = std::max(a.GetCoordim(), b.GetCoordim());
}

// _this = a x b
void PointGeom::Mult(PointGeom &a, PointGeom &b)
{
    (*this)(0) = a[1] * b[2] - a[2] * b[1];
    (*this)(1) = a[2] * b[0] - a[0] * b[2];
    (*this)(2) = a[0] * b[1] - a[1] * b[0];
    m_coordim = 3;
}

// _output = this.a
NekDouble PointGeom::dist(PointGeom &a)
{
    return sqrt((x() - a.x()) * (x() - a.x()) + (y() - a.y()) * (y() - a.y()) +
                (z() - a.z()) * (z() - a.z()));
}

// _output = this.a
NekDouble PointGeom::dot(PointGeom &a)
{
    return (x() * a.x() + y() * a.y() + z() * a.z());
}

/// Determine equivalence by the ids.  No matter what the position,
/// if the ids are the same, then they are equivalent, and vice versa.
bool operator==(const PointGeom &x, const PointGeom &y)
{
    return (x.m_globalID == y.m_globalID);
}

bool operator==(const PointGeom &x, const PointGeom *y)
{
    return (x.m_globalID == y->m_globalID);
}

bool operator==(const PointGeom *x, const PointGeom &y)
{
    return (x->m_globalID == y.m_globalID);
}

bool operator!=(const PointGeom &x, const PointGeom &y)
{
    return (x.m_globalID != y.m_globalID);
}

bool operator!=(const PointGeom &x, const PointGeom *y)
{
    return (x.m_globalID != y->m_globalID);
}

bool operator!=(const PointGeom *x, const PointGeom &y)
{
    return (x->m_globalID != y.m_globalID);
}

PointGeomSharedPtr PointGeom::v_GetVertex(int i) const
{
    ASSERTL0(i == 0, "Index other than 0 is meaningless.");
    // shared_this_ptr() returns const PointGeom, which cannot be
    // returned.
    return PointGeomSharedPtr(new PointGeom(*this));
}

void PointGeom::v_GenGeomFactors()
{
}

}
}
