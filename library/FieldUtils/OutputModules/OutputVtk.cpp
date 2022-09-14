////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputVtk.cpp
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
//  Description: VTK file format output using VTK library.
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>
#include <boost/format.hpp>

#include "OutputVtk.h"
#include <FieldUtils/ProcessModules/ProcessEquiSpacedOutput.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>

#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey OutputVtk::m_className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eOutputModule, "vtu"), OutputVtk::create, "Writes a VTU file.");

OutputVtk::OutputVtk(FieldSharedPtr f) : OutputVtkBase(std::move(f))
{
    m_requireEquiSpaced = true;
    m_config["legacy"] =
        ConfigOption(true, "0", "Output using legacy manual file writers");
    m_config["highorder"] = ConfigOption(
        true, "0", "Output using new high-order Lagrange elements");
    m_config["multiblock"] = ConfigOption(
        true, "0", "Output using multi-blocks to separate composites");
    m_config["uncompress"] = ConfigOption(true, "0", "Uncompress xml sections");
}

// Anonymous namespace for spectral order -> VTK order mapping functions
// be aware that these orderings are very similar to gmsh but not the same!
namespace
{
// Hashing function for the ordering map taking two integers
inline size_t key2(int i, int j)
{
    return (size_t)i << 32 | (unsigned int)j;
}

// Hashing function for the ordering map taking three integers
inline size_t key3(int i, int j, int k)
{
    return (size_t)i << 10 ^ j << 5 ^ k;
}

// Map that takes (a,b,c) -> m
typedef std::tuple<int, int, int> Mode;

struct cmpop
{
    bool operator()(Mode const &a, Mode const &b) const
    {
        if (std::get<0>(a) < std::get<0>(b))
        {
            return true;
        }
        if (std::get<0>(a) > std::get<0>(b))
        {
            return false;
        }
        if (std::get<1>(a) < std::get<1>(b))
        {
            return true;
        }
        if (std::get<1>(a) > std::get<1>(b))
        {
            return false;
        }
        if (std::get<2>(a) < std::get<2>(b))
        {
            return true;
        }

        return false;
    }
};

void Rotate(int nrot, std::vector<long long> &surfVerts)
{
    int n, i, j, cnt;
    int np = static_cast<int>(
        (sqrt(8.0 * static_cast<int>(surfVerts.size()) + 1.0) - 1) / 2);
    std::vector<long long> tmp(np * np);

    for (n = 0; n < nrot; ++n)
    {
        for (cnt = i = 0; i < np; ++i)
        {
            for (j = 0; j < np - i; ++j, cnt++)
            {
                tmp[i * np + j] = surfVerts[cnt];
            }
        }
        for (cnt = i = 0; i < np; ++i)
        {
            for (j = 0; j < np - i; ++j, cnt++)
            {
                surfVerts[cnt] = tmp[(np - 1 - i - j) * np + i];
            }
        }
    }
}

void Reflect(std::vector<long long> &surfVerts)
{
    int i, j, cnt;
    int np = static_cast<int>(
        (sqrt(8.0 * static_cast<double>(surfVerts.size()) + 1.0) - 1) / 2);
    std::vector<long long> tmp(np * np);

    for (cnt = i = 0; i < np; ++i)
    {
        for (j = 0; j < np - i; ++j, cnt++)
        {
            tmp[i * np + np - i - 1 - j] = surfVerts[cnt];
        }
    }

    for (cnt = i = 0; i < np; ++i)
    {
        for (j = 0; j < np - i; ++j, cnt++)
        {
            surfVerts[cnt] = tmp[i * np + j];
        }
    }
}

void Align(std::vector<long long> thisVertId, std::vector<long long> vertId,
           std::vector<long long> &surfVerts)
{
    if (vertId[0] == thisVertId[0])
    {
        if (vertId[1] == thisVertId[1] || vertId[1] == thisVertId[2])
        {
            if (vertId[1] == thisVertId[2])
            {
                Rotate(1, surfVerts);
                Reflect(surfVerts);
            }
        }
    }
    else if (vertId[0] == thisVertId[1])
    {
        if (vertId[1] == thisVertId[0] || vertId[1] == thisVertId[2])
        {
            if (vertId[1] == thisVertId[0])
            {
                Reflect(surfVerts);
            }
            else
            {
                Rotate(2, surfVerts);
            }
        }
    }
    else if (vertId[0] == thisVertId[2])
    {
        if (vertId[1] == thisVertId[0] || vertId[1] == thisVertId[1])
        {
            if (vertId[1] == thisVertId[1])
            {
                Rotate(2, surfVerts);
                Reflect(surfVerts);
            }
            else
            {
                Rotate(1, surfVerts);
            }
        }
    }
}

std::vector<long long> quadTensorNodeOrdering(
    const std::vector<long long> &nodes)
{
    int nN = static_cast<int>(nodes.size());
    int n  = static_cast<int>(sqrt(nN));

    std::vector<long long> nodeList(nN);

    // Vertices
    nodeList[0] = nodes[0];
    if (n > 1)
    {
        nodeList[n - 1]       = nodes[1];
        nodeList[n * n - 1]   = nodes[2];
        nodeList[n * (n - 1)] = nodes[3];
    }

    if (n > 2)
    {
        // Edge 0 -> 1
        for (int i = 1; i < n - 1; ++i)
        {
            nodeList[i] = nodes[4 + i - 1];
        }

        // Edge 3 -> 2
        for (int i = 1; i < n - 1; ++i)
        {
            nodeList[n * (n - 1) + i] = nodes[4 + 2 * (n - 2) + i - 1];
        }

        for (int j = 1; j < n - 1; ++j)
        {
            // Edge 0 -> 3
            nodeList[n * (n - j - 1)] = nodes[4 + 3 * (n - 2) + n - 2 - j];

            // Interior
            for (int i = 1; i < n - 1; ++i)
            {
                nodeList[j * n + i] = nodes[(i - 3) - 2 * j + (3 + j) * n];
            }

            // Edge 1 -> 2
            nodeList[(j + 1) * n - 1] = nodes[4 + (n - 2) + j - 1];
        }
    }

    return nodeList;
}

std::vector<long long> triTensorNodeOrdering(
    const std::vector<long long> &nodes)
{
    int nN = static_cast<int>(nodes.size());
    int n  = static_cast<int>(std::sqrt(2 * nN));

    std::vector<long long> nodeList(nN);
    int cnt2;

    // Vertices
    nodeList[0] = nodes[0];
    if (n > 1)
    {
        nodeList[n - 1]               = nodes[1];
        nodeList[n * (n + 1) / 2 - 1] = nodes[2];
    }

    // Edges
    int cnt = n;
    for (int i = 1; i < n - 1; ++i)
    {
        nodeList[i]               = nodes[3 + i - 1];
        nodeList[cnt]             = nodes[3 + 3 * (n - 2) - i];
        nodeList[cnt + n - i - 1] = nodes[3 + (n - 2) + i - 1];
        cnt += n - i;
    }

    // Interior (recursion)
    if (n > 3)
    {
        // Reorder interior nodes
        std::vector<long long> interior((n - 3) * (n - 2) / 2);
        std::copy(nodes.begin() + 3 + 3 * (n - 2), nodes.end(),
                  interior.begin());
        interior = triTensorNodeOrdering(interior);

        // Copy into full node list
        cnt  = n;
        cnt2 = 0;
        for (int j = 1; j < n - 2; ++j)
        {
            for (int i = 0; i < n - j - 2; ++i)
            {
                nodeList[cnt + i + 1] = interior[cnt2 + i];
            }
            cnt += n - j;
            cnt2 += n - 2 - j;
        }
    }

    return nodeList;
}

std::vector<long long> tetTensorNodeOrdering(
    const std::vector<long long> &nodes)
{
    int nN = static_cast<int>(nodes.size());
    int n  = static_cast<int>(
        std::cbrt(3 * nN + std::sqrt(9 * nN * nN - 1 / 27)) +
        std::cbrt(3 * nN - std::sqrt(9 * nN * nN - 1 / 27)) - 0.5);

    std::vector<long long> nodeList(nN);
    int nTri = n * (n + 1) / 2;
    int nTet = n * (n + 1) * (n + 2) / 6;

    // Vertices
    nodeList[0] = nodes[0];
    if (n == 1)
    {
        return nodeList;
    }

    nodeList[n - 1]    = nodes[1];
    nodeList[nTri - 1] = nodes[2];
    nodeList[nTet - 1] = nodes[3];

    if (n == 2)
    {
        return nodeList;
    }

    // Set up a map that takes (a,b,c) -> m to help us figure out where things
    // are inside the tetrahedron.
    std::map<Mode, int, cmpop> tmp;
    for (int k = 0, cnt = 0; k < n; ++k)
    {
        for (int j = 0; j < n - k; ++j)
        {
            for (int i = 0; i < n - k - j; ++i)
            {
                tmp[Mode(i, j, k)] = cnt++;
            }
        }
    }

    // Edges first 3
    for (int i = 1; i < n - 1; ++i)
    {
        int eI                               = i - 1;
        nodeList[tmp[Mode(i, 0, 0)]]         = nodes[4 + eI];
        nodeList[tmp[Mode(n - 1 - i, i, 0)]] = nodes[4 + (n - 2) + eI];
        nodeList[tmp[Mode(0, n - 1 - i, 0)]] = nodes[4 + 2 * (n - 2) + eI];
    }

    // Edges last 3 reversed (compared with NekMesh)
    for (int i = 1; i < n - 1; ++i)
    {
        int eI                               = (n - 1 - i) - 1;
        nodeList[tmp[Mode(0, 0, n - 1 - i)]] = nodes[4 + 3 * (n - 2) + eI];
        nodeList[tmp[Mode(i, 0, n - 1 - i)]] = nodes[4 + 4 * (n - 2) + eI];
        nodeList[tmp[Mode(0, i, n - 1 - i)]] = nodes[4 + 5 * (n - 2) + eI];
    }

    if (n == 3)
    {
        return nodeList;
    }

    // For faces, we use the triTensorNodeOrdering routine to make our lives
    // slightly easier.
    int nFacePts = (n - 3) * (n - 2) / 2;

    // Grab face points and reorder into a tensor-product type format
    std::vector<std::vector<long long>> tmpNodes(4);
    int offset = 4 + 6 * (n - 2);

    for (int i = 0; i < 4; ++i)
    {
        tmpNodes[i].resize(nFacePts);
        for (int j = 0; j < nFacePts; ++j)
        {
            tmpNodes[i][j] = nodes[offset++];
        }
        tmpNodes[i] = triTensorNodeOrdering(tmpNodes[i]);
    }

    if (n > 4)
    {
        // Now align faces (different to NekMesh)
        std::vector<long long> triVertId(3), toAlign(3);
        triVertId[0] = 0;
        triVertId[1] = 1;
        triVertId[2] = 2;

        toAlign[0] = 0;
        toAlign[1] = 2;
        toAlign[2] = 1;
        Align(triVertId, toAlign, tmpNodes[2]);
        Align(triVertId, toAlign, tmpNodes[3]);

        toAlign[0] = 2;
        toAlign[1] = 0;
        toAlign[2] = 1;
        Align(triVertId, toAlign, tmpNodes[1]);
    }

    // Reordered from NekMesh to put base last
    for (int j = 1, cnt = 0; j < n - 2; ++j)
    {
        for (int i = 1; i < n - j - 1; ++i, ++cnt)
        {
            nodeList[tmp[Mode(i, j, 0)]]             = tmpNodes[3][cnt];
            nodeList[tmp[Mode(i, 0, j)]]             = tmpNodes[0][cnt];
            nodeList[tmp[Mode(n - 1 - i - j, i, j)]] = tmpNodes[1][cnt];
            nodeList[tmp[Mode(0, i, j)]]             = tmpNodes[2][cnt];
        }
    }

    if (n == 4)
    {
        return nodeList;
    }

    // Finally, recurse on interior volume
    std::vector<long long> intNodes, tmpInt;
    for (int i = offset; i < nTet; ++i)
    {
        intNodes.emplace_back(nodes[i]);
    }

    tmpInt = tetTensorNodeOrdering(intNodes);

    for (int k = 1, cnt = 0; k < n - 2; ++k)
    {
        for (int j = 1; j < n - k - 1; ++j)
        {
            for (int i = 1; i < n - k - j - 1; ++i)
            {
                nodeList[tmp[Mode(i, j, k)]] = tmpInt[cnt++];
            }
        }
    }

    return nodeList;
}

// This prism mapping is NOT using the standard spectral mapping! The point
// output from ProcessEquiSpacedOutput is not the same as the rest of the code.
std::vector<long long> prismTensorNodeOrdering(
    const std::vector<long long> &nodes)
{
    int nN = static_cast<int>(nodes.size());
    int n  = static_cast<int>(std::cbrt(2 * nN));

    std::vector<long long> nodeList(nN);
    int nQuad  = n * n;
    int nPrism = n * n * (n + 1) / 2;
    int edge   = n - 2;

    // Vertices
    nodeList[0] = nodes[0];
    if (n > 1)
    {
        nodeList[nPrism - n] = nodes[1];
        nodeList[n - 1]      = nodes[2];
        nodeList[nQuad - n]  = nodes[3];
        nodeList[nPrism - 1] = nodes[4];
        nodeList[nQuad - 1]  = nodes[5];
    }

    if (n == 2)
    {
        return nodeList;
    }

    int nPts = 6;
    int cnt  = 0;
    for (int i = 0; i < n - 2; ++i)
    {
        // Triangle 1 edges
        cnt += n * (n - i); // This is adding the quads as they reduce a side by
                            // one each time
        nodeList[cnt]            = nodes[nPts + i];                // Edge 1
        nodeList[cnt + edge - i] = nodes[nPts + 2 * edge - i - 1]; // Edge 2
        nodeList[i + 1]          = nodes[nPts + 3 * edge - i - 1]; // Edge 3
    }

    cnt = nQuad;
    for (int i = 1; i < n - 1; ++i)
    {
        // Triangle 2 edges
        cnt += n * (n - i); // This is adding the quads as they reduce a side by
                            // one each time
        nodeList[cnt - n + i]   = nodes[nPts + 3 * edge + i - 1]; // Edge 1 (4)
        nodeList[cnt - 1]       = nodes[nPts + 5 * edge - i];     // Edge 2 (5)
        nodeList[nQuad - i - 1] = nodes[nPts + 5 * edge + i - 1]; // Edge 3 (6)
    }

    // Vertical edges
    for (int i = 1; i < n - 1; ++i)
    {
        // Vertical prism edges
        nodeList[n * i]          = nodes[nPts + 6 * edge + i - 1]; // Edge 1 (7)
        nodeList[nPrism - n + i] = nodes[nPts + 7 * edge + i - 1]; // Edge 2 (8)
        nodeList[n * i + n - 1]  = nodes[nPts + 8 * edge + i - 1]; // Edge 3 (9)
    }

    // Do quad faces first as when n == 3 we have no tri faces
    int nSmallTri  = (n - 3) * (n - 2) / 2;
    int nSmallQuad = (n - 2) * (n - 2);
    nPts           = 6 + 9 * edge + 2 * nSmallTri;
    cnt            = nQuad;
    for (int i = 1; i < n - 1; ++i)
    {
        for (int j = 0; j < n - 2; ++j)
        {
            nodeList[cnt + (n - i) + (n - i) * j] =
                nodes[nPts + j * edge + (i - 1)];

            nodeList[cnt + 2 * (n - i) - 1 + (n - i) * j] =
                nodes[nPts + nSmallQuad + (j + 1) * edge - i];

            nodeList[i * n + j + 1] =
                nodes[nPts + 2 * nSmallQuad + i * edge - j - 1];
        }

        cnt += n * (n - i);
    }

    if (n == 3)
    {
        return nodeList;
    }

    // Triangular alignment
    std::vector<long long> tmpNodes(nSmallTri);
    std::iota(tmpNodes.begin(), tmpNodes.end(), 0);
    if (n > 4)
    {
        std::vector<long long> triVertId(3), toAlign(3);
        triVertId[0] = 0;
        triVertId[1] = 1;
        triVertId[2] = 2;

        toAlign[0] = 0;
        toAlign[1] = 2;
        toAlign[2] = 1;
        Align(triVertId, toAlign, tmpNodes);
    }

    // Triangle face 1
    nPts     = 6 + 9 * edge;
    cnt      = nQuad;
    int cnt2 = 0;
    for (int i = 1; i < n - 2; ++i)
    {
        for (int j = 1; j < n - i - 1; ++j)
        {
            nodeList[cnt + j] = nPts + tmpNodes[cnt2++];
        }

        cnt += n * (n - i);
    }

    // Triangle face 2
    cnt  = nQuad;
    cnt2 = 0;
    for (int i = 1; i < n - 2; ++i)
    {
        cnt += n * (n - i);
        for (int j = 1; j < n - i - 1; ++j)
        {
            nodeList[cnt - (n - i) + j] = nPts + nSmallTri + tmpNodes[cnt2++];
        }
    }

    // Triangular volume (uses same alignment mapping)
    nPts = 6 + 9 * edge + 2 * nSmallTri + 3 * nSmallQuad;
    for (int k = 1; k < n - 1; ++k)
    {
        cnt  = nQuad;
        cnt2 = 0;
        for (int i = 1; i < n - 2; ++i)
        {
            for (int j = 1; j < n - i - 1; ++j)
            {
                nodeList[cnt + k * (n - i) + j] = nPts + tmpNodes[cnt2++];
            }

            cnt += n * (n - i);
        }

        nPts += nSmallTri;
    }

    return nodeList;
}

std::vector<long long> hexTensorNodeOrdering(
    const std::vector<long long> &nodes)
{
    int nN = static_cast<int>(nodes.size());
    int n  = static_cast<int>(std::cbrt(nN));

    std::vector<long long> nodeList(nN);

    // Vertices
    nodeList[0] = nodes[0];
    if (n == 1)
    {
        return nodeList;
    }

    // Vertices: same order as Nektar++
    nodeList[n - 1]                         = nodes[1];
    nodeList[n * n - 1]                     = nodes[2];
    nodeList[n * (n - 1)]                   = nodes[3];
    nodeList[n * n * (n - 1)]               = nodes[4];
    nodeList[n - 1 + n * n * (n - 1)]       = nodes[5];
    nodeList[n * n - 1 + n * n * (n - 1)]   = nodes[6];
    nodeList[n * (n - 1) + n * n * (n - 1)] = nodes[7];

    if (n == 2)
    {
        return nodeList;
    }

    int hexEdges[12][2] = {{0, 1},
                           {n - 1, n},
                           {n * n - 1, -1},
                           {n * (n - 1), -n},
                           {0, n * n},
                           {n - 1, n * n},
                           {n * n - 1, n * n},
                           {n * (n - 1), n * n},
                           {n * n * (n - 1), 1},
                           {n * n * (n - 1) + n - 1, n},
                           {n * n * n - 1, -1},
                           {n * n * (n - 1) + n * (n - 1), -n}};
    int hexFaces[6][3]  = {{0, 1, n},         {0, 1, n * n},
                          {n - 1, n, n * n}, {n * (n - 1), 1, n * n},
                          {0, n, n * n},     {n * n * (n - 1), 1, n}};

    int gmshToNekEdge[12] = {0, 1, -2, -3, 8, 9, -10, -11, 4, 5, 6, 7};

    // Edges
    int offset = 8;
    for (int i : gmshToNekEdge)
    {
        int e = abs(i);

        if (i >= 0)
        {
            for (int j = 1; j < n - 1; ++j)
            {
                nodeList[hexEdges[e][0] + j * hexEdges[e][1]] = nodes[offset++];
            }
        }
        else
        {
            for (int j = 1; j < n - 1; ++j)
            {
                nodeList[hexEdges[e][0] + (n - j - 1) * hexEdges[e][1]] =
                    nodes[offset++];
            }
        }
    }

    // Faces
    int gmsh2NekFace[6] = {4, 2, 1, 3, 0, 5};

    // Map which defines orientation between Gmsh and Nektar++ faces.
    StdRegions::Orientation faceOrient[6] = {
        StdRegions::eDir1FwdDir2_Dir2FwdDir1,
        StdRegions::eDir1FwdDir1_Dir2FwdDir2,
        StdRegions::eDir1FwdDir2_Dir2FwdDir1,
        StdRegions::eDir1FwdDir1_Dir2FwdDir2,
        StdRegions::eDir1BwdDir1_Dir2FwdDir2,
        StdRegions::eDir1FwdDir1_Dir2FwdDir2};

    for (int i = 0; i < 6; ++i)
    {
        int n2   = (n - 2) * (n - 2);
        int face = gmsh2NekFace[i];
        offset   = 8 + 12 * (n - 2) + i * n2;

        // Create a list of interior face nodes for this face only.
        std::vector<long long> faceNodes(n2);
        for (int j = 0; j < n2; ++j)
        {
            faceNodes[j] = nodes[offset + j];
        }

        // Now get the reordering of this face, which puts Gmsh
        // recursive ordering into Nektar++ row-by-row order.
        std::vector<long long> tmp(n2);

        // Finally reorient the face according to the geometry
        // differences.
        if (faceOrient[i] == StdRegions::eDir1FwdDir1_Dir2FwdDir2)
        {
            // Orientation is the same, just copy.
            tmp = faceNodes;
        }
        else if (faceOrient[i] == StdRegions::eDir1FwdDir2_Dir2FwdDir1)
        {
            // Tranposed faces
            for (int j = 0; j < n - 2; ++j)
            {
                for (int k = 0; k < n - 2; ++k)
                {
                    tmp[j * (n - 2) + k] = faceNodes[k * (n - 2) + j];
                }
            }
        }
        else if (faceOrient[i] == StdRegions::eDir1BwdDir1_Dir2FwdDir2)
        {
            for (int j = 0; j < n - 2; ++j)
            {
                for (int k = 0; k < n - 2; ++k)
                {
                    tmp[j * (n - 2) + k] = faceNodes[j * (n - 2) + (n - k - 3)];
                }
            }
        }

        // Now put this into the right place in the output array
        for (int k = 1; k < n - 1; ++k)
        {
            for (int j = 1; j < n - 1; ++j)
            {
                nodeList[hexFaces[face][0] + j * hexFaces[face][1] +
                         k * hexFaces[face][2]] =
                    faceNodes[(k - 1) * (n - 2) + j - 1];
            }
        }
    }

    // Finally, recurse on interior volume
    std::vector<long long> intNodes, tmpInt;
    for (int i = 8 + 12 * (n - 2) + 6 * (n - 2) * (n - 2); i < n * n * n; ++i)
    {
        intNodes.push_back(nodes[i]);
    }

    if (!intNodes.empty())
    {
        tmpInt = hexTensorNodeOrdering(intNodes);
        for (int k = 1, cnt = 0; k < n - 1; ++k)
        {
            for (int j = 1; j < n - 1; ++j)
            {
                for (int i = 1; i < n - 1; ++i)
                {
                    nodeList[i + j * n + k * n * n] = tmpInt[cnt++];
                }
            }
        }
    }

    return nodeList;
}

std::vector<long long> lowOrderMapping(int nDim, Array<OneD, int> nquad)
{
    std::vector<long long> p;
    switch (nDim)
    {
        case 3:
            for (int j = 0; j < nquad[0] - 1; ++j)
            {
                for (int k = 0; k < nquad[1] - 1; ++k)
                {
                    for (int l = 0; l < nquad[2] - 1; ++l)
                    {
                        p.insert(
                            p.end(),
                            {l * nquad[0] * nquad[1] + k * nquad[0] + j,
                             l * nquad[0] * nquad[1] + k * nquad[0] + j + 1,
                             l * nquad[0] * nquad[1] + (k + 1) * nquad[0] + j +
                                 1,
                             l * nquad[0] * nquad[1] + (k + 1) * nquad[0] + j,
                             (l + 1) * nquad[0] * nquad[1] + k * nquad[0] + j,
                             (l + 1) * nquad[0] * nquad[1] + k * nquad[0] + j +
                                 1,
                             (l + 1) * nquad[0] * nquad[1] +
                                 (k + 1) * nquad[0] + j + 1,
                             (l + 1) * nquad[0] * nquad[1] +
                                 (k + 1) * nquad[0] + j});
                    }
                }
            }
            break;
        case 2:
            for (int j = 0; j < nquad[0] - 1; ++j)
            {
                for (int k = 0; k < nquad[1] - 1; ++k)
                {
                    p.insert(p.end(), {k * nquad[0] + j, k * nquad[0] + j + 1,
                                       (k + 1) * nquad[0] + j + 1,
                                       (k + 1) * nquad[0] + j});
                }
            }
            break;
        case 1:
            for (int j = 0; j < nquad[0] - 1; ++j)
            {
                p.insert(p.end(), {j, j + 1});
            }
            break;
        default:
            break;
    }

    return p;
}

int GetHighOrderVtkCellType(int sType)
{
#if VTK_MAJOR_VERSION >= 8
    // For deformed elements this is a map of shape type to VTK type
    static const std::map<int, int> vtkCellType = {
        {{LibUtilities::eSegment, VTK_LAGRANGE_CURVE},
         {LibUtilities::eTriangle, VTK_LAGRANGE_TRIANGLE},
         {LibUtilities::eQuadrilateral, VTK_LAGRANGE_QUADRILATERAL},
         {LibUtilities::eTetrahedron, VTK_LAGRANGE_TETRAHEDRON},
         {LibUtilities::ePyramid, VTK_LAGRANGE_PYRAMID},
         {LibUtilities::ePrism, VTK_LAGRANGE_WEDGE},
         {LibUtilities::eHexahedron, VTK_LAGRANGE_HEXAHEDRON}}};

    return vtkCellType.at(sType);
#else
    boost::ignore_unused(sType);
    NEKERROR(
        ErrorUtil::efatal,
        "High-order VTK output requires minimum VTK library version of 8.0")
    return 0;
#endif
}

} // namespace

vtkSmartPointer<vtkUnstructuredGrid> OutputVtk::OutputFromExpLowOrder()
{
    // Mesh file to output
    vtkSmartPointer<vtkUnstructuredGrid> vtkMesh =
        vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Cache ordering so we aren't recreating mappings
    std::unordered_map<size_t, std::vector<long long>> mappingCache;

    // Choose cell types and num quad points based on dimension
    int nDim     = m_f->m_exp[0]->GetShapeDimension();
    int nHomoDir = m_f->m_numHomogeneousDir;

    auto type = (nDim + nHomoDir == 3)   ? VTK_HEXAHEDRON
                : (nDim + nHomoDir == 2) ? VTK_QUAD
                                         : VTK_LINE;
    int nQpts = (nDim + nHomoDir == 3) ? 8 : (nDim + nHomoDir == 2) ? 4 : 2;

    // Fill homogeneous info
    if (nHomoDir != 0)
    {
        ASSERTL0(nHomoDir == 1 &&
                     m_f->m_exp[0]->GetExpType() == MultiRegions::e3DH1D,
                 "Only regular expansions and the 3DH1D homogeneous expansion "
                 "are supported in the new VTK writer. Please use the 'legacy' "
                 "option for all other expansion types.")

        m_numPlanes = m_f->m_exp[0]->GetHomogeneousBasis()->GetNumModes();
        // Extra plane if fourier
        m_extraPlane = (m_f->m_exp[0]->GetHomogeneousBasis()->GetBasisType() ==
                            LibUtilities::eFourier &&
                        m_f->m_exp[0]->GetHomogeneousBasis()->GetPointsType() ==
                            LibUtilities::eFourierEvenlySpaced);
    }

    // Insert points
    vtkNew<vtkPoints> vtkPoints;
    vtkPoints->SetDataType(VTK_DOUBLE);

    int offset = 0;
    int nElmts = static_cast<int>(m_f->m_exp[0]->GetNumElmts());
    for (int i = 0; i < nElmts; ++i)
    {
        Array<OneD, int> nquad(3, 0);
        for (int j = 0; j < nDim; ++j)
        {
            nquad[j] = m_f->m_exp[0]->GetExp(i)->GetNumPoints(j);
        }

        // Add homogeneous direction num points
        if (nHomoDir != 0)
        {
            nquad[nDim] = m_numPlanes;

            // Add extra plane if fourier
            if (m_extraPlane)
            {
                nquad[nDim]++;
            }
        }

        // Get total points by multiplying all nquads together (accounts for
        // homo)
        int nPts = nquad[0];
        for (int j = 1; j < nDim + nHomoDir; ++j)
        {
            nPts *= nquad[j];
        }

        // Get all coordinates at quadrature points in this element
        Array<OneD, NekDouble> x(nPts, 0.0), y(nPts, 0.0), z(nPts, 0.0);
        m_f->m_exp[0]->GetCoords(i, x, y, z);

        // If add extra plane for fourier need to fill last plane coordinates
        if (m_extraPlane)
        {
            // Copy x & y to extra plane
            Array<OneD, NekDouble> tmp;
            Vmath::Vcopy(nquad[0] * nquad[1], x, 1,
                         tmp = x + (nquad[2] - 1) * nquad[0] * nquad[1], 1);
            Vmath::Vcopy(nquad[0] * nquad[1], y, 1,
                         tmp = y + (nquad[2] - 1) * nquad[0] * nquad[1], 1);

            // Fill z on extra plane
            NekDouble zDouble = z[nquad[0] * nquad[1] * (nquad[2] - 1) - 1] +
                                (z[nquad[0] * nquad[1]] - z[0]);
            Vmath::Fill(nquad[0] * nquad[1], zDouble,
                        tmp = z + (nquad[2] - 1) * nquad[0] * nquad[1], 1);
        }

        // Insert points in to vtk mesh
        for (int j = 0; j < nPts; ++j)
        {
            vtkPoints->InsertNextPoint(x[j], y[j], z[j]);
        }

        // Insert elements
        auto oIt = mappingCache.find(key3(nquad[0], nquad[1], nquad[2]));
        if (oIt == mappingCache.end())
        {
            auto p = lowOrderMapping(nDim + nHomoDir, nquad);
            oIt    = mappingCache
                      .insert(
                          std::make_pair(key3(nquad[0], nquad[1], nquad[2]), p))
                      .first;
        }

        auto p = oIt->second;
        std::for_each(p.begin(), p.end(),
                      [offset](long long &d) { d += offset; });
        for (int j = 0; j < p.size(); j += nQpts)
        {
            vtkMesh->InsertNextCell(type, nQpts, &p[j]);
        }

        offset += nPts;
    }

    vtkMesh->SetPoints(vtkPoints.GetPointer());

    return vtkMesh;
}

void OutputVtk::AddFieldDataToVTKLowOrder(
    po::variables_map &vm, std::string &filename,
    vtkSmartPointer<vtkUnstructuredGrid> &vtkMesh)
{
    // Insert field information we iterate by element from first plane to last,
    // while the fld file iterates by plane from first element to last
    int nElmts         = static_cast<int>(m_f->m_exp[0]->GetNumElmts());
    int numPtsPerPlane = m_f->m_exp[0]->GetTotPoints() / m_numPlanes;
    for (int v = 0; v < m_f->m_variables.size(); ++v)
    {
        vtkNew<vtkDoubleArray> fieldData;
        for (int i = 0; i < nElmts; ++i)
        {
            int elmtOffset          = m_f->m_exp[v]->GetPhys_Offset(i);
            int nPtsPerElmtPerPlane = m_f->m_exp[v]->GetExp(i)->GetTotPoints();

            for (int j = 0; j < m_numPlanes; ++j)
            {
                int planeOffset = j * numPtsPerPlane;
                for (int k = 0; k < nPtsPerElmtPerPlane; ++k)
                {
                    fieldData->InsertNextValue(
                        m_f->m_exp[v]->GetPhys()[elmtOffset + planeOffset + k]);
                }
            }

            // if extra plane we copy the first plane values in to last plane
            if (m_extraPlane)
            {
                for (int k = 0; k < nPtsPerElmtPerPlane; ++k)
                {
                    fieldData->InsertNextValue(
                        m_f->m_exp[v]->GetPhys()[elmtOffset + k]);
                }
            }
        }

        fieldData->SetName(&m_f->m_variables[v][0]);
        vtkMesh->GetPointData()->AddArray(fieldData.GetPointer());
    }

    WriteVTK(vtkMesh, filename, vm);
}

void OutputVtk::OutputFromExpLowOrderMultiBlock(po::variables_map &vm,
                                                std::string &filename)
{
    ASSERTL0(
        m_f->m_numHomogeneousDir == 0,
        "Multi block VTK is not implemented for homogeneous expansion types.")

    ASSERTL0(m_f->m_comm->IsSerial(),
             "Multi block VTK is not implemented in parallel.")

    int dim = m_f->m_graph->GetMeshDimension();

    // Create mappings from geometry id to expansion ids
    std::array<std::map<int, std::pair<int, int>>, 4> geomIdToExpId;
    int nElmts = static_cast<int>(m_f->m_exp[0]->GetNumElmts());
    for (int i = 0; i < nElmts; ++i)
    {
        auto geom = m_f->m_exp[0]->GetExp(i)->GetGeom();
        geomIdToExpId[geom->GetShapeDim()][geom->GetGlobalID()] =
            std::make_pair(i, -1);

        for (int j = 0; j < geom->GetNumFaces(); ++j)
        {
            geomIdToExpId[2][geom->GetFid(j)] = std::make_pair(i, j);
        }

        for (int j = 0; j < geom->GetNumEdges(); ++j)
        {
            geomIdToExpId[1][geom->GetEid(j)] = std::make_pair(i, j);
        }
    }

    // Cache ordering so we aren't recreating mappings
    std::unordered_map<size_t, std::vector<long long>> mappingCache;

    std::map<int, vtkNew<vtkUnstructuredGrid>> vtkMesh;
    std::map<int, SpatialDomains::CompositeSharedPtr> composites =
        m_f->m_graph->GetComposites();
    std::map<int, std::string> compositeNames;
    for (auto &comp : composites)
    {
        // Vector of field data
        std::vector<vtkNew<vtkDoubleArray>> fieldData(m_f->m_variables.size());

        // Insert points
        vtkNew<vtkPoints> vtkPoints;
        vtkPoints->SetDataType(VTK_DOUBLE);

        int compId = comp.first;
        std::vector<std::shared_ptr<SpatialDomains::Geometry>> geomVec =
            comp.second->m_geomVec;

        unsigned int offset = 0;
        for (auto &geom : geomVec)
        {
            int geomId = geom->GetGlobalID();
            int sDim   = geom->GetShapeDim();
            auto type  = (sDim == 3)   ? VTK_HEXAHEDRON
                         : (sDim == 2) ? VTK_QUAD
                                       : VTK_LINE;
            int nQpts  = (sDim == 3) ? 8 : (sDim == 2) ? 4 : 2;

            LocalRegions::ExpansionSharedPtr exp =
                m_f->m_exp[0]->GetExp(geomIdToExpId[sDim][geomId].first);

            unsigned int nPts = exp->GetTotPoints();
            Array<OneD, NekDouble> x(nPts, 0.0), y(nPts, 0.0), z(nPts, 0.0);
            exp->GetCoords(x, y, z);

            int offsetPhys = m_f->m_exp[0]->GetPhys_Offset(
                geomIdToExpId[sDim][geomId].first);
            if (sDim == dim)
            {
                for (int j = 0; j < nPts; ++j)
                {
                    vtkPoints->InsertNextPoint(x[j], y[j], z[j]);

                    // Add field data
                    for (int k = 0; k < m_f->m_variables.size(); ++k)
                    {
                        fieldData[k]->InsertNextValue(
                            m_f->m_exp[k]->GetPhys()[j + offsetPhys]);
                    }
                }
            }
            else
            {
                Array<OneD, int> pointsMap;
                exp->GetTracePhysMap(geomIdToExpId[sDim][geomId].second,
                                     pointsMap);
                for (int j : pointsMap)
                {
                    vtkPoints->InsertNextPoint(x[j], y[j], z[j]);

                    // Add field data
                    for (int k = 0; k < m_f->m_variables.size(); ++k)
                    {
                        fieldData[k]->InsertNextValue(
                            m_f->m_exp[k]->GetPhys()[offsetPhys + j]);
                    }
                }

                exp  = exp->GetTraceExp(geomIdToExpId[sDim][geomId].second);
                nPts = pointsMap.size();
            }

            Array<OneD, int> nquad(3, 0);
            for (int j = 0; j < sDim; ++j)
            {
                nquad[j] = exp->GetNumPoints(j);
            }

            auto oIt = mappingCache.find(key3(nquad[0], nquad[1], nquad[2]));
            if (oIt == mappingCache.end())
            {
                auto p = lowOrderMapping(sDim, nquad);
                oIt    = mappingCache
                          .insert(std::make_pair(
                              key3(nquad[0], nquad[1], nquad[2]), p))
                          .first;
            }

            auto p = oIt->second;
            std::for_each(p.begin(), p.end(),
                          [offset](long long &d) { d += offset; });

            for (int j = 0; j < p.size(); j += nQpts)
            {
                vtkMesh[compId]->InsertNextCell(type, nQpts, &p[j]);
            }

            offset += nPts;
        }

        vtkMesh[compId]->SetPoints(vtkPoints.GetPointer());

        // Add all fields to vtkMesh
        for (int k = 0; k < m_f->m_variables.size(); ++k)
        {
            fieldData[k]->SetName(&m_f->m_variables[k][0]);
            vtkMesh[compId]->GetPointData()->AddArray(
                fieldData[k].GetPointer());
        }

        // Find composite names if exists and store in map
        // Set name as boundary label if it exists otherwise index ID
        compositeNames[compId] = "Composite ID " + std::to_string(compId);
        auto clabels           = m_f->m_graph->GetCompositesLabels();
        auto oIt               = clabels.find(compId);
        if (oIt != clabels.end())
        {
            compositeNames[compId] = oIt->second;
        }
    }

    // Main multi-block set
    vtkNew<vtkMultiBlockDataSet> mainBlock;

    std::set<int> compSet;

    // Fill domain multi blocks from composites
    vtkNew<vtkMultiBlockDataSet> mainDomainBlock;
    auto domains = m_f->m_graph->GetDomain();
    std::vector<vtkNew<vtkMultiBlockDataSet>> domainMultiBlocks(domains.size());
    for (int i = 0; i < domains.size(); ++i)
    {
        auto dom = domains[i];

        // Loop over composites and see if in domain
        for (auto &comp : composites)
        {
            int compId = comp.first;
            if (dom.find(compId) != dom.end())
            {
                unsigned int nBlock = domainMultiBlocks[i]->GetNumberOfBlocks();
                domainMultiBlocks[i]->SetBlock(nBlock,
                                               vtkMesh[compId].GetPointer());
                domainMultiBlocks[i]->GetMetaData(nBlock)->Set(
                    vtkCompositeDataSet::NAME(),
                    compositeNames[compId].c_str());
                compSet.insert(compId);
            }
        }

        unsigned int nBlock = mainDomainBlock->GetNumberOfBlocks();
        mainDomainBlock->SetBlock(nBlock, domainMultiBlocks[i].GetPointer());
        mainDomainBlock->GetMetaData(nBlock)->Set(
            vtkCompositeDataSet::NAME(),
            ("Domain ID " + std::to_string(i)).c_str());
    }

    if (mainDomainBlock->GetNumberOfBlocks() != 0)
    {
        auto nBlock = static_cast<unsigned int>(mainBlock->GetNumberOfBlocks());
        mainBlock->SetBlock(nBlock, mainDomainBlock.GetPointer());
        mainBlock->GetMetaData(nBlock)->Set(vtkCompositeDataSet::NAME(),
                                            "Domains");
    }

    // Fill boundary multi blocks from composites
    int cnt = 0;
    if(m_f->m_session)
    {
        SpatialDomains::BoundaryConditions bcs(m_f->m_session,
                                               m_f->m_exp[0]->GetGraph());
        const SpatialDomains::BoundaryRegionCollection &bregions =
            bcs.GetBoundaryRegions();

        vtkNew<vtkMultiBlockDataSet> mainBoundaryBlock;
        std::vector<vtkNew<vtkMultiBlockDataSet>> boundaryMultiBlocks(
            bregions.size());
        for (auto &boundary : bregions)
        {
            // Loop over composites and see if in boundary
            for (auto &comp : composites)
            {
                int compId = comp.first;
                if (boundary.second->find(compId) != boundary.second->end())
                {
                    unsigned int nBlock =
                        boundaryMultiBlocks[cnt]->GetNumberOfBlocks();
                    boundaryMultiBlocks[cnt]->SetBlock(
                        nBlock, vtkMesh[compId].GetPointer());
                    boundaryMultiBlocks[cnt]->GetMetaData(nBlock)->Set(
                        vtkCompositeDataSet::NAME(),
                        compositeNames[compId].c_str());
                    compSet.insert(compId);
                }
            }

            unsigned int nBlock = mainBoundaryBlock->GetNumberOfBlocks();
            mainBoundaryBlock->SetBlock(
                nBlock, boundaryMultiBlocks[cnt++].GetPointer());

            // Set name as boundary label if it exists otherwise index ID
            std::string name = "Boundary ID " + std::to_string(boundary.first);
            auto blabels     = bcs.GetBoundaryLabels();
            auto oIt         = blabels.find(boundary.first);
            if (oIt != blabels.end())
            {
                name = oIt->second;
            }

            mainBoundaryBlock->GetMetaData(nBlock)->Set(
                vtkCompositeDataSet::NAME(), name.c_str());
        }

        if (mainBoundaryBlock->GetNumberOfBlocks() != 0)
        {
            auto nBlock =
                static_cast<unsigned int>(mainBlock->GetNumberOfBlocks());
            mainBlock->SetBlock(nBlock, mainBoundaryBlock.GetPointer());
            mainBlock->GetMetaData(nBlock)->Set(vtkCompositeDataSet::NAME(),
                                                "Boundaries");
        }
    }

    // Include all other composites (not domains or boundaries)
    vtkNew<vtkMultiBlockDataSet> compsMultiBlocks;
    cnt = 0;
    for (auto &comp : composites)
    {
        int compId = comp.first;
        if (compSet.find(compId) == compSet.end())
        {
            unsigned int nBlock = compsMultiBlocks->GetNumberOfBlocks();
            compsMultiBlocks->SetBlock(nBlock, vtkMesh[compId].GetPointer());
            compsMultiBlocks->GetMetaData(nBlock)->Set(
                vtkCompositeDataSet::NAME(), compositeNames[compId].c_str());
        }
    }

    if (compsMultiBlocks->GetNumberOfBlocks() != 0)
    {
        auto nBlock = static_cast<unsigned int>(mainBlock->GetNumberOfBlocks());
        mainBlock->SetBlock(nBlock, compsMultiBlocks.GetPointer());
        mainBlock->GetMetaData(nBlock)->Set(vtkCompositeDataSet::NAME(),
                                            "Other composites");
    }

    WriteVTK(mainBlock.GetPointer(), filename, vm);
}

vtkSmartPointer<vtkUnstructuredGrid> OutputVtk::OutputFromExpHighOrder(
    po::variables_map &vm)
{
    ASSERTL0(
        m_f->m_numHomogeneousDir == 0,
        "High order VTK is not implemented for homogeneous expansion types.")

    // Save shapetypes before convert to equispaced because that nukes the
    // explist
    std::vector<LibUtilities::ShapeType> sType;
    for (auto &i : *m_f->m_exp[0]->GetExp())
    {
        sType.emplace_back(i->GetGeom()->GetShapeType());
    }

    // Convert expansion to an equispaced points field
    auto equispaced =
        MemoryManager<ProcessEquiSpacedOutput>::AllocateSharedPtr(m_f);
    equispaced->Process(vm);

    // Insert points
    vtkNew<vtkPoints> vtkPoints;
    LibUtilities::PtsFieldSharedPtr fPts = m_f->m_fieldPts;
    vtkPoints->SetDataType(VTK_DOUBLE);

    Array<OneD, Array<OneD, NekDouble>> pts;
    fPts->GetPts(pts);

    int dim  = static_cast<int>(fPts->GetDim());
    int nPts = static_cast<int>(fPts->GetNpoints());
    switch (dim)
    {
        case 3:
            for (int i = 0; i < nPts; ++i)
            {
                vtkPoints->InsertNextPoint(pts[0][i], pts[1][i], pts[2][i]);
            }
            break;
        case 2:
            for (int i = 0; i < nPts; ++i)
            {
                vtkPoints->InsertNextPoint(pts[0][i], pts[1][i], 0.0);
            }
            break;
        case 1:
            for (int i = 0; i < nPts; ++i)
            {
                vtkPoints->InsertNextPoint(pts[0][i], 0.0, 0.0);
            }
            break;
        default:
            break;
    }

    // Mesh file to output
    vtkSmartPointer<vtkUnstructuredGrid> vtkMesh =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkMesh->SetPoints(vtkPoints.GetPointer());

    // Cache ordering for shape type & npts so we aren't recreating mappings
    std::unordered_map<size_t, std::vector<long long>> mappingCache;

    // Get offset per element in a vector from ppe
    std::vector<int> ppe = m_f->m_fieldPts->GetPointsPerElement();
    Array<OneD, int> ppeOffset(ppe.size() + 1, 0.0);
    std::partial_sum(ppe.begin(), ppe.end(), &ppeOffset[1]);

    // Insert elements
    for (int i = 0; i < ppe.size(); ++i)
    {
        // Construct inverse of input reordering.
        // First try to find it in our mapping cache.
        auto oIt = mappingCache.find(key2(sType[i], ppe[i]));
        if (oIt == mappingCache.end())
        {
            std::vector<long long> p(ppe[i]);
            std::iota(p.begin(), p.end(), 0);
            switch (sType[i])
            {
                case LibUtilities::eQuadrilateral:
                    p = quadTensorNodeOrdering(p);
                    break;
                case LibUtilities::eTriangle:
                    p = triTensorNodeOrdering(p);
                    break;
                case LibUtilities::eTetrahedron:
                    p = tetTensorNodeOrdering(p);
                    break;
                case LibUtilities::eHexahedron:
                    p = hexTensorNodeOrdering(p);
                    break;
                case LibUtilities::ePrism:
                    p = prismTensorNodeOrdering(p);
                    break;
                default:
                    NEKERROR(ErrorUtil::efatal,
                             "High-order VTU output not set up for the " +
                                 static_cast<std::string>(
                                     LibUtilities::ShapeTypeMap[sType[i]]) +
                                 " shape type.");
                    break;
            }

            // Invert the ordering as this is for spectral -> recursive (VTU)
            std::vector<long long> inv(ppe[i]);
            for (int j = 0; j < ppe[i]; ++j)
            {
                inv[p[j]] = j;
            }

            oIt =
                mappingCache.insert(std::make_pair(key2(sType[i], ppe[i]), inv))
                    .first;
        }

        // Add offset to reordering
        std::vector<long long> p = oIt->second;
        for (long long &j : p)
        {
            j += ppeOffset[i];
        }

        vtkMesh->InsertNextCell(GetHighOrderVtkCellType(sType[i]), ppe[i],
                                &p[0]);
    }

    return vtkMesh;
}

void OutputVtk::AddFieldDataToVTKHighOrder(
    po::variables_map &vm, std::string &filename,
    vtkSmartPointer<vtkUnstructuredGrid> &vtkMesh)
{
    // Convert expansion to an equispaced points field
    if (m_cachedMesh)
    {
        auto equispaced =
            MemoryManager<ProcessEquiSpacedOutput>::AllocateSharedPtr(m_f);
        equispaced->Process(vm);
    }
    LibUtilities::PtsFieldSharedPtr fPts = m_f->m_fieldPts;

    int nPts = static_cast<int>(fPts->GetNpoints());
    int dim  = static_cast<int>(fPts->GetDim());

    Array<OneD, Array<OneD, NekDouble>> pts;
    fPts->GetPts(pts);

    // Insert field information
    for (int i = 0; i < fPts->GetNFields(); ++i)
    {
        vtkNew<vtkDoubleArray> fieldData;
        fieldData->SetArray(&pts[dim + i][0], nPts, 1);
        fieldData->SetName(&fPts->GetFieldName(i)[0]);
        vtkMesh->GetPointData()->AddArray(fieldData.GetPointer());
    }

    WriteVTK(vtkMesh, filename, vm);
}

void OutputVtk::WriteVTK(vtkDataObject *vtkMesh, std::string &filename,
                         po::variables_map &vm)
{
    // Initialise base writer class as default structured grid writer
    vtkSmartPointer<vtkXMLWriter> vtkMeshWriter =
        vtkNew<vtkXMLUnstructuredGridWriter>().GetPointer();

    // If multi block we need to switch to the MultiBlock writer
    if (m_config["multiblock"].m_beenSet)
    {
        vtkMeshWriter = vtkNew<vtkXMLMultiBlockDataWriter>().GetPointer();
    }

    // Choose the correct extension based on which mesh writer is used
    int dot          = static_cast<int>(filename.find_last_of('.'));
    std::string body = filename.substr(0, dot + 1);
    filename         = body + vtkMeshWriter->GetDefaultFileExtension();

    // Set time information if exists
    if (!m_f->m_fieldMetaDataMap["Time"].empty())
    {
        vtkMesh->GetInformation()->Set(
            vtkDataObject::DATA_TIME_STEP(),
            std::stod(m_f->m_fieldMetaDataMap["Time"]));
    }
    else if (!m_f->m_fieldMetaDataMap["FinalTime"].empty())
    {
        vtkMesh->GetInformation()->Set(
            vtkDataObject::DATA_TIME_STEP(),
            std::stod(m_f->m_fieldMetaDataMap["FinalTime"]));
    }

    // Write out the new mesh in XML format (don't support legacy
    // format here as we still have standard OutputVtk.cpp)
    vtkMeshWriter->SetFileName(filename.c_str());

#if VTK_MAJOR_VERSION <= 5
    vtkMeshWriter->SetInput(vtkMesh);
#else
    vtkMeshWriter->SetInputData(vtkMesh);
#endif

    if (m_config["uncompress"].m_beenSet)
    {
        vtkMeshWriter->SetDataModeToAscii();
    }

    vtkMeshWriter->Update();

    std::cout << "Written file: " << filename << std::endl;

    // Output parallel file information if needed - using a manual write.
    // We could use the VTK lib to do this, but that requires VTK with MPI
    // enabled & messing about with the parallel controller & changing
    // our file naming scheme as VTK forces _${proc-number} as a suffix...
    if (m_f->m_comm->TreatAsRankZero() && !m_f->m_comm->IsSerial())
    {
        WritePVtu(vm);
    }
}

void OutputVtk::WritePVtu(po::variables_map &vm)
{
    std::string filename = m_config["outfile"].as<std::string>();
    int dot              = static_cast<int>(filename.find_last_of('.'));
    std::string body     = filename.substr(0, dot);
    filename             = body + ".pvtu";

    std::ofstream outfile(filename.c_str());

    int nprocs = m_f->m_comm->GetSize();
    std::string path =
        LibUtilities::PortablePath(OutputVtkBase::GetPath(filename, vm));

    outfile << "<?xml version=\"1.0\"?>" << endl;
    outfile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" "
            << "byte_order=\"LittleEndian\">" << endl;
    outfile << "  <PUnstructuredGrid GhostLevel=\"0\">" << endl;

    // Add time if exists
    if (!m_f->m_fieldMetaDataMap["Time"].empty())
    {
        outfile << "    <FieldData> " << endl;
        outfile << "      <DataArray type=\"Float64\" Name=\"TimeValue\" "
                   "NumberOfTuples=\"1\" format=\"ascii\">"
                << endl;
        outfile << m_f->m_fieldMetaDataMap["Time"] << std::endl;
        outfile << "      </DataArray>" << std::endl;
        outfile << "    </FieldData> " << endl;
    }

    // Add point coordinates
    outfile << "    <PPoints> " << endl;
    outfile << "      <PDataArray type=\"Float64\" NumberOfComponents=\"" << 3
            << "\"/> " << endl;
    outfile << "    </PPoints>" << endl;

    // Add cell information
    outfile << "    <PCells>" << endl;
    outfile << "      <PDataArray type=\"Int64\" Name=\"connectivity\" "
               "NumberOfComponents=\"1\"/>"
            << endl;
    outfile << "      <PDataArray type=\"Int64\" Name=\"offsets\"      "
               "NumberOfComponents=\"1\"/>"
            << endl;
    outfile << "      <PDataArray type=\"UInt8\" Name=\"types\"        "
               "NumberOfComponents=\"1\"/>"
            << endl;
    outfile << "    </PCells>" << endl;

    // Add fields (point data)
    outfile << "    <PPointData>" << endl;
    for (auto &var : m_f->m_variables)
    {
        outfile << "      <PDataArray type=\"Float64\" Name=\"" << var << "\"/>"
                << endl;
    }
    outfile << "    </PPointData>" << endl;

    // Add parallel files
    for (int i = 0; i < nprocs; ++i)
    {
        boost::format pad("P%1$07d.vtu");
        pad % i;
        outfile << "    <Piece Source=\"" << path << "/" << pad.str() << "\"/>"
                << endl;
    }
    outfile << "  </PUnstructuredGrid>" << endl;
    outfile << "</VTKFile>" << endl;

    cout << "Written file: " << filename << endl;
}

void OutputVtk::OutputFromData(po::variables_map &vm)
{
    boost::ignore_unused(vm);
    NEKERROR(ErrorUtil::efatal, "OutputVtk can't write using only FieldData.");
}

void OutputVtk::OutputFromPts(po::variables_map &vm)
{
    OutputVtkBase::OutputFromPts(vm);
}

void OutputVtk::OutputFromExp(po::variables_map &vm)
{
    if (m_config["legacy"].m_beenSet)
    {
        ASSERTL0(!m_config["multiblock"].m_beenSet,
                 "Multi block VTK is not implemented for legacy output.")

        ASSERTL0(!m_config["highorder"].m_beenSet,
                 "High order VTK is not implemented for legacy output.")

        // No caching of mesh data in legacy output
        OutputVtkBase::OutputFromExp(vm);
        return;
    }

    // Extract the output filename and extension
    std::string filename = OutputVtkBase::PrepareOutput(vm);

    // Save mesh state (if using filter this allows us to only ProcessEquispaced
    // if needed)
    m_cachedMesh = m_vtkMesh;

    if (m_config["highorder"].m_beenSet)
    {
        ASSERTL0(!m_config["multiblock"].m_beenSet,
                 "Multi block VTK is not implemented for high-order output.")

        if (!m_cachedMesh)
        {
            m_vtkMesh = OutputFromExpHighOrder(vm);
        }

        AddFieldDataToVTKHighOrder(vm, filename, m_vtkMesh);
    }
    else if (m_config["multiblock"].m_beenSet)
    {
        // No mesh caching for multiblock due to complexity of field data
        OutputFromExpLowOrderMultiBlock(vm, filename);
    }
    else
    {
        if (!m_cachedMesh)
        {
            m_vtkMesh = OutputFromExpLowOrder();
        }

        AddFieldDataToVTKLowOrder(vm, filename, m_vtkMesh);
    }
}

} // namespace FieldUtils
} // namespace Nektar
