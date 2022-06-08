////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputVtkNew.cpp
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
//  Description: VTK file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include "OutputVtkNew.h"
#include <FieldUtils/ProcessModules/ProcessEquiSpacedOutput.h>

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey OutputVtkNew::m_className = GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eOutputModule, "vtuNew"), OutputVtkNew::create, "Writes a VTU file.");

OutputVtkNew::OutputVtkNew(FieldSharedPtr f) : OutputVtk(std::move(f))
{
    m_requireEquiSpaced = true;
    m_config["highorder"] = ConfigOption(true, "0", "Output using new high-order Lagrange elements");
    m_config["uncompress"] = ConfigOption(true, "0", "Uncompress xml sections");
    m_config["compressionlevel"] = ConfigOption(false, "5", "Compression level for the VTU output: 1-9");
}

// Anonymous namespace for spectral order -> VTK order mapping functions
namespace
{
// Hashing function for the ordering map
inline size_t key(int i,int j) {return (size_t) i << 32 | (unsigned int) j;}

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

void Align(std::vector<long long> thisVertId,
                 std::vector<long long> vertId,
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

std::vector<long long> quadTensorNodeOrdering(const std::vector<long long> &nodes)
{
    int nN = static_cast<int>(nodes.size());
    int n = static_cast<int>(sqrt(nN));

    std::vector<long long> nodeList(nN);

    // Vertices
    nodeList[0] = nodes[0];
    if (n > 1)
    {
        nodeList[n - 1]     = nodes[1];
        nodeList[n * n - 1] = nodes[2];
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

            //Edge 1 -> 2
            nodeList[(j + 1) * n - 1] = nodes[4 + (n - 2) + j - 1];
        }
    }

    return nodeList;
}

std::vector<long long> triTensorNodeOrdering(const std::vector<long long> &nodes)
{
    int nN = static_cast<int>(nodes.size());
    int n = static_cast<int>(std::sqrt(2 * nN));

    std::vector<long long> nodeList(nN);
    int cnt2;

    // Vertices
    nodeList[0] = nodes[0];
    if (n > 1)
    {
        nodeList[n - 1] = nodes[1];
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
        std::copy(
                nodes.begin() + 3 + 3 * (n - 2), nodes.end(), interior.begin());
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

std::vector<long long> tetTensorNodeOrdering(const std::vector<long long> &nodes)
{
    int nN = static_cast<int>(nodes.size());
    int n = static_cast<int>(std::cbrt(3 * nN + std::sqrt(9 * nN * nN - 1 / 27)) +
                             std::cbrt(3 * nN - std::sqrt(9 * nN * nN - 1 / 27)) - 0.5);

    std::vector<long long> nodeList(nN);
    int nTri = n*(n+1)/2;
    int nTet = n*(n+1)*(n+2)/6;

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
                tmp[Mode(i,j,k)] = cnt++;
            }
        }
    }

    // Edges first 3
    for (int i = 1; i < n-1; ++i)
    {
        int eI = i-1;
        nodeList[tmp[Mode(i,0,0)]]     = nodes[4 + eI];
        nodeList[tmp[Mode(n-1-i,i,0)]] = nodes[4 + (n-2) + eI];
        nodeList[tmp[Mode(0,n-1-i,0)]] = nodes[4 + 2*(n-2) + eI];
    }

    // Edges last 3 reversed (compared with NekMesh)
    for (int i = 1; i < n-1; ++i)
    {
        int eI = (n - 1 - i) - 1;
        nodeList[tmp[Mode(0,0,n-1-i)]] = nodes[4 + 3*(n-2) + eI];
        nodeList[tmp[Mode(i,0,n-1-i)]] = nodes[4 + 4*(n-2) + eI];
        nodeList[tmp[Mode(0,i,n-1-i)]] = nodes[4 + 5*(n-2) + eI];
    }

    if (n == 3)
    {
        return nodeList;
    }

    // For faces, we use the triTensorNodeOrdering routine to make our lives
    // slightly easier.
    int nFacePts = (n-3)*(n-2)/2;

    // Grab face points and reorder into a tensor-product type format
    std::vector<std::vector<long long>> tmpNodes(4);
    int offset = 4 + 6*(n-2);

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
    for (int j = 1, cnt = 0; j < n-2; ++j)
    {
        for (int i = 1; i < n-j-1; ++i, ++cnt)
        {
            nodeList[tmp[Mode(i,j,0)]]       = tmpNodes[3][cnt];
            nodeList[tmp[Mode(i,0,j)]]       = tmpNodes[0][cnt];
            nodeList[tmp[Mode(n-1-i-j,i,j)]] = tmpNodes[1][cnt];
            nodeList[tmp[Mode(0,i,j)]]       = tmpNodes[2][cnt];
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
                nodeList[tmp[Mode(i,j,k)]] = tmpInt[cnt++];
            }
        }
    }

    return nodeList;
}

}

int OutputVtkNew::GetVtkCellType(int sType, SpatialDomains::GeomType gType)
{
    static const std::map<int, std::map<int, int>> vtkCellType = {
        {SpatialDomains::eRegular,
            {
                {LibUtilities::eSegment, VTK_LINE},
                {LibUtilities::eTriangle, VTK_TRIANGLE},
                {LibUtilities::eQuadrilateral, VTK_QUAD},
                {LibUtilities::eTetrahedron, VTK_TETRA},
                {LibUtilities::ePyramid, VTK_PYRAMID},
                {LibUtilities::ePrism, VTK_WEDGE},
                {LibUtilities::eHexahedron, VTK_HEXAHEDRON}
            }
        },
        {SpatialDomains::eDeformed,
            {
                 {LibUtilities::eSegment, VTK_LAGRANGE_CURVE},
                 {LibUtilities::eTriangle, VTK_LAGRANGE_TRIANGLE},
                 {LibUtilities::eQuadrilateral, VTK_LAGRANGE_QUADRILATERAL},
                 {LibUtilities::eTetrahedron, VTK_LAGRANGE_TETRAHEDRON},
                 {LibUtilities::ePyramid, VTK_LAGRANGE_PYRAMID},
                 {LibUtilities::ePrism, VTK_LAGRANGE_WEDGE},
                 {LibUtilities::eHexahedron, VTK_LAGRANGE_HEXAHEDRON}
            }
        }
    };

   return vtkCellType.at(gType).at(sType);
}

void OutputVtkNew::OutputFromData(po::variables_map &vm)
{
    boost::ignore_unused(vm);
    NEKERROR(ErrorUtil::efatal, "OutputVtkNew can't write using only FieldData.");
}

void OutputVtkNew::OutputFromPts(po::variables_map &vm)
{
    boost::ignore_unused(vm);
    NEKERROR(ErrorUtil::efatal, "OutputVtkNew can't write using only PtsData.");
}

void OutputVtkNew::OutputFromExp(po::variables_map &vm)
{
    // Extract the output filename and extension
    std::string filename = OutputVtk::PrepareOutput(vm);

    if (m_config["highorder"].m_beenSet)
    {
        OutputFromExpHighOrder(vm, filename);
    }
    else
    {
        OutputFromExpLowOrder(vm, filename);
    }
}

void OutputVtkNew::OutputFromExpLowOrder(po::variables_map &vm, std::string &filename)
{
    boost::ignore_unused(vm);

    // Insert points
    vtkNew<vtkPoints> vtkPoints;
    vtkNew<vtkUnstructuredGrid> vtkMesh;

    int nElmts = static_cast<int>(m_f->m_exp[0]->GetNumElmts());
    for (int i = 0; i < nElmts; ++i)
    {
        int nPts = static_cast<int>(m_f->m_exp[0]->GetExp(i)->GetTotPoints());

        Array<OneD, NekDouble> x(nPts), y(nPts), z(nPts);
        m_f->m_exp[0]->GetExp(i)->GetCoords(x, y, z);

        for (int j = 0; j < nPts; ++j)
        {
          vtkPoints->InsertNextPoint(x[j], y[j], z[j]);
        }

        int nDim = m_f->m_exp[0]->GetShapeDimension();

        Array<OneD, NekDouble> nquad(nDim);
        for (int j = 0; j < nDim; ++j)
        {
            nquad[j] = m_f->m_exp[0]->GetExp(i)->GetBasisNumModes(0);
        }

        // Insert elements
        int offset = m_f->m_exp[0]->GetPhys_Offset(i);
        switch(nDim)
        {
            case 3:
                for (int j = 0; j < nquad[0] - 1; ++j)
                {
                    for (int k = 0; k < nquad[1] - 1; ++k)
                    {
                        for (int l = 0; l < nquad[2] - 1; ++l)
                        {
                            long long p[8] = {
                                    l * nquad[0] * nquad[1] + k * nquad[0] + j + offset,
                                    l * nquad[0] * nquad[1] + k * nquad[0] + j + 1 + offset,
                                    l * nquad[0] * nquad[1] + (k + 1) * nquad[0] + j + 1 + offset,
                                    l * nquad[0] * nquad[1] + (k + 1) * nquad[0] + j + offset,
                                    (l + 1) * nquad[0] * nquad[1] + k * nquad[0] + j + offset,
                                    (l + 1) * nquad[0] * nquad[1] + k * nquad[0] + j + 1 + offset,
                                    (l + 1) * nquad[0] * nquad[1] + (k + 1) * nquad[0] + j + 1 + offset,
                                    (l + 1) * nquad[0] * nquad[1] + (k + 1) * nquad[0] + j + offset
                            };

                            vtkMesh->InsertNextCell(GetVtkCellType(LibUtilities::eHexahedron,
                                                                   SpatialDomains::eRegular), 8, p);
                        }
                    }
                }
                break;
            case 2:
                for (int j = 0; j < nquad[0] - 1; ++j)
                {
                    for (int k = 0; k < nquad[1] - 1; ++k)
                    {
                        long long p[4] = {
                                k * nquad[0] + j + offset,
                                k * nquad[0] + j + 1 + offset,
                                (k + 1) * nquad[0] + j + 1 + offset,
                                (k + 1) * nquad[0] + j + offset
                        };

                        vtkMesh->InsertNextCell(GetVtkCellType(LibUtilities::eQuadrilateral,
                                                               SpatialDomains::eRegular), 4, p);
                    }
                }
                break;
            case 1:
                for (int j = 0; j < nquad[0] - 1; ++j)
                {
                    long long p[2] = {
                            j + offset,
                            j + 1 + offset
                    };

                    vtkMesh->InsertNextCell(GetVtkCellType(LibUtilities::eSegment,
                                                           SpatialDomains::eRegular), 2, p);
                }
                break;
            default:
                break;
        }
    }

    vtkMesh->SetPoints(vtkPoints);

    // Insert field information
    int totPts = m_f->m_exp[0]->GetTotPoints();
    for (int i = 0; i < m_f->m_variables.size(); ++i)
    {
        vtkNew<vtkDoubleArray> fieldData;
        for (int j = 0; j < totPts; ++j)
        {
            fieldData->InsertNextValue(m_f->m_exp[i]->GetPhys()[j]);
        }

        fieldData->SetName(&m_f->m_variables[i][0]);
        vtkMesh->GetPointData()->AddArray(fieldData);
    }

    WriteVTK(vtkMesh, filename);
}

void OutputVtkNew::OutputFromExpHighOrder(po::variables_map &vm, std::string &filename)
{
    // Save shapetypes before convert to equispaced because that nukes the explist
    std::vector<LibUtilities::ShapeType> sType;
    for (auto &i : *m_f->m_exp[0]->GetExp())
    {
        sType.emplace_back(i->GetGeom()->GetShapeType());
    }

    // Convert expansion to an equispaced points field
    auto equispaced = MemoryManager<ProcessEquiSpacedOutput>::AllocateSharedPtr(m_f);
    equispaced->Process(vm);

    // Insert points
    vtkNew<vtkPoints> vtkPoints;
    LibUtilities::PtsFieldSharedPtr fPts = m_f->m_fieldPts;

    Array<OneD, Array<OneD, NekDouble>> pts;
    fPts->GetPts(pts);

    int dim = static_cast<int>(fPts->GetDim());
    int nPts = static_cast<int>(fPts->GetNpoints());
    switch (dim)
    {
        case 3:
            for (int i = 0; i < nPts; ++i)
            {
                vtkPoints->InsertNextPoint(pts[0][i],pts[1][i],pts[2][i]);
            }
            break;
        case 2:
            for (int i = 0; i < nPts; ++i)
            {
                vtkPoints->InsertNextPoint(pts[0][i],pts[1][i], 0.0);
            }
            break;
        case 1:
            for (int i = 0; i < nPts; ++i)
            {
                vtkPoints->InsertNextPoint(pts[0][i],0.0,0.0);
            }
            break;
        default:
            break;
    }

    vtkNew<vtkUnstructuredGrid> vtkMesh;
    vtkMesh->SetPoints(vtkPoints);

    // Cache ordering for shape type & npts so we aren't recreating mappings
    std::unordered_map<size_t, std::vector<long long>> mappingCache;

    // Get points per element and offset per element in a vector
    std::vector<int> ppe = m_f->m_fieldPts->GetPointsPerElement();
    Array<OneD, int> ppeOffset(ppe.size() + 1, 0.0);
    std::partial_sum(ppe.begin(), ppe.end(), &ppeOffset[1]);

    // Insert elements
    for (int i = 0; i < ppe.size(); ++i)
    {
        // Construct inverse of input reordering.
        // First try to find it in our mapping cache.
        auto oIt = mappingCache.find(key(sType[i], ppe[i]));
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
                default:
                    NEKERROR(ErrorUtil::efatal,
                             "VTU output not set up for this shape type.");
                    break;
            }

            // Invert the ordering as this is for spectral -> recursive (VTU)
            std::vector<long long> inv(ppe[i]);
            for (int j = 0; j < ppe[i]; ++j)
            {
                inv[p[j]] = j;
            }

            oIt = mappingCache.insert(std::make_pair(key(sType[i], ppe[i]), inv)).first;
        }

        // Add offset to reordering
        std::vector<long long> p = oIt->second;
        std::for_each(p.begin(), p.end(),[j = ppeOffset[i]](long long &d) { d += j; });

        vtkMesh->InsertNextCell(
                GetVtkCellType(sType[i], SpatialDomains::eDeformed),
                ppe[i], &p[0]);
    }

    // Insert field information
    for (int i = 0; i < fPts->GetNFields(); ++i)
    {
        vtkNew<vtkDoubleArray> fieldData;
        fieldData->SetArray(&pts[dim + i][0], nPts, 1);
        fieldData->SetName(&fPts->GetFieldName(i)[0]);
        vtkMesh->GetPointData()->AddArray(fieldData);
    }

    WriteVTK(vtkMesh, filename);
}

void OutputVtkNew::WriteVTK(vtkUnstructuredGrid* vtkMesh, std::string &filename)
{
    // Write out the new mesh in XML format (don't support legacy
    // format here as we still have standard OutputVtk.cpp)

    vtkNew<vtkXMLUnstructuredGridWriter> vtkMeshWriter;
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

    if (m_config["compressionlevel"].m_beenSet)
    {
        vtkMeshWriter->SetCompressionLevel(
                std::stoi(m_config["compressionlevel"].m_value));
    }

    vtkMeshWriter->Update();

    cout << "Written file: " << filename << endl;
}

}
}
