////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshGraphXml.h
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
//  Description:  This file contains the base class specification for the
//                MeshGraphXml class.
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_MGXML_H
#define NEKTAR_SPATIALDOMAINS_MGXML_H

#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/MeshPartition.h>

namespace Nektar
{
namespace SpatialDomains
{

class MeshGraphXml : public MeshGraph
{
public:
    MeshGraphXml()
    {
    }

    virtual ~MeshGraphXml()
    {
    }

    SPATIAL_DOMAINS_EXPORT void WriteXMLGeometry(
        std::string outname, std::vector<std::set<unsigned int>> elements,
        std::vector<unsigned int> partitions);

    static MeshGraphSharedPtr create()
    {
        return MemoryManager<MeshGraphXml>::AllocateSharedPtr();
    }

    static std::string className;

protected:
    // some of these functions are going to be virtual because they will be
    // inherited by the XmlCompressed version

    SPATIAL_DOMAINS_EXPORT virtual void v_WriteGeometry(
        std::string &outfilename, bool defaultExp = false,
        const LibUtilities::FieldMetaDataMap &metadata =
            LibUtilities::NullFieldMetaDataMap) override;
    SPATIAL_DOMAINS_EXPORT virtual void v_ReadGeometry(
        LibUtilities::DomainRangeShPtr rng, bool fillGraph) override;
    SPATIAL_DOMAINS_EXPORT virtual void v_PartitionMesh(
        LibUtilities::SessionReaderSharedPtr session) override;

    virtual void v_ReadVertices();
    virtual void v_ReadCurves();
    void ReadDomain();

    virtual void v_ReadEdges();
    virtual void v_ReadFaces();

    void ReadElements();
    void ReadComposites();

    virtual void v_ReadElements1D();
    virtual void v_ReadElements2D();
    virtual void v_ReadElements3D();

    void ResolveGeomRef(const std::string &prevToken, const std::string &token,
                        CompositeSharedPtr &composite);
    void ResolveGeomRef1D(const std::string &prevToken,
                          const std::string &token,
                          CompositeSharedPtr &composite);
    void ResolveGeomRef2D(const std::string &prevToken,
                          const std::string &token,
                          CompositeSharedPtr &composite);
    void ResolveGeomRef3D(const std::string &prevToken,
                          const std::string &token,
                          CompositeSharedPtr &composite);

    virtual void v_WriteVertices(TiXmlElement *geomTag, PointGeomMap &verts);
    virtual void v_WriteEdges(TiXmlElement *geomTag, SegGeomMap &edges);
    virtual void v_WriteTris(TiXmlElement *faceTag, TriGeomMap &tris);
    virtual void v_WriteQuads(TiXmlElement *faceTag, QuadGeomMap &quads);
    virtual void v_WriteHexs(TiXmlElement *elmtTag, HexGeomMap &hexs);
    virtual void v_WritePrisms(TiXmlElement *elmtTag, PrismGeomMap &pris);
    virtual void v_WritePyrs(TiXmlElement *elmtTag, PyrGeomMap &pyrs);
    virtual void v_WriteTets(TiXmlElement *elmtTag, TetGeomMap &tets);
    virtual void v_WriteCurves(TiXmlElement *geomTag, CurveMap &edges,
                               CurveMap &faces);
    void WriteComposites(TiXmlElement *geomTag, CompositeMap &comps,
                         std::map<int, std::string> &compLabels);
    void WriteDomain(TiXmlElement *geomTag,
                     std::map<int, CompositeMap> &domain);
    void WriteDefaultExpansion(TiXmlElement *root);

    CompositeOrdering CreateCompositeOrdering();
};

} // namespace SpatialDomains
} // namespace Nektar

#endif
