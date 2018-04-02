////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry.h
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
//  Description:  This file contains the base class specification for the
//                Geometry class.
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_MGHDF5_H
#define NEKTAR_SPATIALDOMAINS_MGHDF5_H

#include "MeshGraph.h"

#include <LibUtilities/BasicUtils/H5.h>

namespace Nektar
{
namespace SpatialDomains
{

class MeshGraphHDF5 : public MeshGraph
{
public:
    MeshGraphHDF5()
    {
    }

    SPATIAL_DOMAINS_EXPORT virtual void WriteGeometry(
        std::string &outfilename,
        bool defaultExp = false,
        const LibUtilities::FieldMetaDataMap &metadata
                                         = LibUtilities::NullFieldMetaDataMap);

    SPATIAL_DOMAINS_EXPORT virtual void WriteGeometry(
        std::string outname,
        std::vector<std::set<unsigned int>> elements,
        std::vector<unsigned int> partitions);

    virtual ~MeshGraphHDF5()
    {
    }

    static MeshGraphSharedPtr create()
    {
        return MemoryManager<MeshGraphHDF5>::AllocateSharedPtr();
    }

    static std::string className;

protected:
    SPATIAL_DOMAINS_EXPORT virtual void ReadGeometry(
        DomainRangeShPtr rng,
        bool             fillGraph);
    SPATIAL_DOMAINS_EXPORT virtual void PartitionMesh(
        LibUtilities::SessionReaderSharedPtr session);

private:

    void ReadVertices();
    void ReadCurves();
    void ReadDomain();

    void ReadEdges();
    void ReadFaces();

    void ReadElements();
    void ReadComposites();

    void WriteVertices(PointGeomMap &verts);
    void WriteEdges(SegGeomMap &edges);
    void WriteTris(TriGeomMap &tris);
    void WriteQuads(QuadGeomMap &quads);
    void WriteHexs(HexGeomMap &hexs);
    void WritePrisms(PrismGeomMap &pris);
    void WritePyrs(PyrGeomMap &pyrs);
    void WriteTets(TetGeomMap &tets);
    void WriteCurves(CurveMap &edges, CurveMap &faces);
    void WriteComposites(CompositeMap &comps);
    void WriteDomain(vector<CompositeMap> &domain);

    string m_hdf5Name;
    LibUtilities::H5::FileSharedPtr m_file;
    LibUtilities::H5::GroupSharedPtr m_mesh;
    LibUtilities::H5::GroupSharedPtr m_maps;
};

} // end of namespace
} // end of namespace

#endif
