////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshGraphXmlCompressed.h
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
//  Description:
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_MGXMLCOM_H
#define NEKTAR_SPATIALDOMAINS_MGXMLCOM_H

#include "MeshGraphXml.h"

namespace Nektar
{
namespace SpatialDomains
{

class MeshGraphXmlCompressed : public MeshGraphXml
{
public:
    MeshGraphXmlCompressed()
    {
    }

    virtual ~MeshGraphXmlCompressed()
    {
    }

    static MeshGraphSharedPtr create()
    {
        return MemoryManager<MeshGraphXmlCompressed>::AllocateSharedPtr();
    }

    static std::string className;

protected:
    virtual void v_ReadVertices() override;
    virtual void v_ReadCurves() override;

    virtual void v_ReadEdges() override;
    virtual void v_ReadFaces() override;

    virtual void v_ReadElements1D() override;
    virtual void v_ReadElements2D() override;
    virtual void v_ReadElements3D() override;

    virtual void v_WriteVertices(TiXmlElement *geomTag,
                                 PointGeomMap &verts) override;
    virtual void v_WriteEdges(TiXmlElement *geomTag,
                              SegGeomMap &edges) override;
    virtual void v_WriteTris(TiXmlElement *faceTag, TriGeomMap &tris) override;
    virtual void v_WriteQuads(TiXmlElement *faceTag,
                              QuadGeomMap &quads) override;
    virtual void v_WriteHexs(TiXmlElement *elmtTag, HexGeomMap &hexs) override;
    virtual void v_WritePrisms(TiXmlElement *elmtTag,
                               PrismGeomMap &pris) override;
    virtual void v_WritePyrs(TiXmlElement *elmtTag, PyrGeomMap &pyrs) override;
    virtual void v_WriteTets(TiXmlElement *elmtTag, TetGeomMap &tets) override;
    virtual void v_WriteCurves(TiXmlElement *geomTag, CurveMap &edges,
                               CurveMap &faces) override;
};

} // namespace SpatialDomains
} // namespace Nektar

#endif
