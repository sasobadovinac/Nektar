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

#include <set>
#include <string>
#include <iomanip>

#include <boost/core/ignore_unused.hpp>
#include <boost/format.hpp>

#include <LibUtilities/BasicUtils/FileSystem.h>

#include "OutputVtkNew.h"

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellType.h>
#include <vtkNew.h>
#include <vtkDoubleArray.h>

// Hashing function for the ordering map
inline size_t key(int i,int j) {return (size_t) i << 32 | (unsigned int) j;}

namespace Nektar
{
namespace FieldUtils
{

ModuleKey OutputVtkNew::m_className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eOutputModule, "vtuNew"), OutputVtkNew::create, "Writes a VTU file.");

OutputVtkNew::OutputVtkNew(FieldSharedPtr f) : OutputVtk(f)
{
    m_requireEquiSpaced = true;
    m_config["uncompress"] = ConfigOption(true, "0", "Uncompress xml sections");
    m_config["compressionlevel"] = ConfigOption(false, "5", "Compression level for the VTU output: 1-9");
}

std::vector<long long> triTensorNodeOrdering(const std::vector<long long> &nodes, int n)
{
    std::vector<long long> nodeList(nodes.size());
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
        interior = triTensorNodeOrdering(interior, n - 3);

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

int OutputVtkNew::GetVtkCellType(int sType, SpatialDomains::GeomType gType)
{
    static std::map<int, std::map<int, int>> vtkCellType = {
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

   return vtkCellType[gType][sType];
}

std::vector<long long> OutputVtkNew::QuadrilateralNodes(int &ppe)
{
    int n = (int)sqrt(ppe);

    std::vector<long long> p(ppe, 0);
    // Write vertices
    p[1] += (n - 1);
    p[2] += (n * n - 1);
    p[3] += (n * (n - 1));

    // Write edge interior
    int cnt = 4;
    for (int j = 0; j < 4; ++j)
    {
        if (j == 0)
        {
            for (int k = 1; k < n - 1; ++k)
            {
                p[cnt++] += k;
            }
        }
        else if (j == 1)
        {
            for (int k = 2 * n - 1; k < n * (n - 1); k+=n)
            {
                p[cnt++] += k;
            }
        }
        else if (j == 2)
        {
            for (int k = n * (n - 1) + 1; k < n * n - 1; ++k)
            {
                p[cnt++] += k;
            }
        }
        else if (j == 3)
        {
            for (int k = n; k < n * (n - 1); k+=n)
            {
                p[cnt++] += k;
            }
        }
    }

    // Write surface interior
    for (int j = 1; j < n - 1; ++j)
    {
        for (int k = 1; k < n - 1; ++k)
        {
            // Fetch interior nodes from quad->curve
            p[cnt++] += (j * n + k);
        }
    }

    return p;
}

std::vector<long long> OutputVtkNew::TriangleNodes(int &ppe)
{
    // Calculate order from triangle number
    //sqrt(2 * ppe + 0.25) - 0.5 -> (int)sqrt(2 * ppe)
    int n = (int)sqrt(2 * ppe);

    std::vector<long long> p(ppe);
    std::iota(p.begin(), p.end(), 0);

    p = triTensorNodeOrdering(p, n);

    // Invert the ordering as this is for spectral -> recursive (VTU) and add offset
    std::vector<long long> inv(ppe);
    for (int j = 0; j < ppe; ++j)
    {
        inv[p[j]] = j;
    }

    return inv;
}

void OutputVtkNew::OutputFromExp(po::variables_map &vm)
{
    boost::ignore_unused(vm);
    NEKERROR(ErrorUtil::efatal, "OutputVtkNew can't write using only ExpData.");
}

void OutputVtkNew::OutputFromPts(po::variables_map &vm)
{
    // Extract the output filename and extension
    std::string filename = OutputVtk::PrepareOutput(vm);

    // Insert points
    vtkNew<vtkPoints> vtkPoints;
    LibUtilities::PtsFieldSharedPtr fPts = m_f->m_fieldPts;

    Array<OneD, Array<OneD, NekDouble>> pts;
    fPts->GetPts(pts);

    int nPts = fPts->GetNpoints();
    for (int i = 0; i < nPts; ++i) // @TODO: Change to allow 3D
    {
        vtkPoints->InsertNextPoint(pts[0][i],pts[1][i],0.0);
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
        auto geom = m_f->m_exp[0]->GetExp(i)->GetGeom();
        auto sType = geom->GetShapeType();

        // Construct inverse of input reordering.
        // First try to find it in our mapping cache.
        auto oIt = mappingCache.find(key(sType, ppe[i]));
        if (oIt == mappingCache.end())
        {
            std::vector<long long> p;
            switch (sType)
            {
                case LibUtilities::eQuadrilateral:
                    p = QuadrilateralNodes(ppe[i]);
                    break;
                case LibUtilities::eTriangle:
                    p = TriangleNodes(ppe[i]);
                    break;
                default:
                    NEKERROR(ErrorUtil::efatal,
                             "VTU output not set up for this shape type.");
                    break;
            }

            oIt = mappingCache.insert(std::make_pair(key(sType, ppe[i]), p)).first;
        }

        // Add offset to reordering
        std::vector<long long> p = oIt->second;
        std::for_each(p.begin(), p.end(),[j = ppeOffset[i]](long long &d) { d += j; });

        vtkMesh->InsertNextCell(
                GetVtkCellType(sType, geom->GetGeomFactors()->GetGtype()),
                ppe[i], &p[0]);
    }

    // Insert field information
    int dim = fPts->GetDim();
    for (int i = 0; i < fPts->GetNFields(); ++i)
    {
        vtkNew<vtkDoubleArray> fieldData;
        fieldData->SetArray(&pts[dim + i][0], nPts, 1);
        fieldData->SetName(&fPts->GetFieldName(i)[0]);
        vtkMesh->GetPointData()->AddArray(fieldData);
    }

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

void OutputVtkNew::OutputFromData(po::variables_map &vm)
{
    boost::ignore_unused(vm);
    NEKERROR(ErrorUtil::efatal, "OutputVtkNew can't write using only FieldData.");
}

}
}
