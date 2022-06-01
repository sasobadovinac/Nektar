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
using namespace std;

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
}

int OutputVtkNew::GetVtkCellType(LibUtilities::ShapeType sType, SpatialDomains::GeomType gType)
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

void OutputVtkNew::OutputFromExp(po::variables_map &vm)
{
    // Extract the output filename and extension
    string filename = OutputVtk::PrepareOutput(vm);
    
    vtkUnstructuredGrid *vtkMesh   = vtkUnstructuredGrid::New();
    vtkPoints           *vtkPoints = vtkPoints::New();

    // Save geometry information to VTU (assuming same expansion for each field!)
    int meshDim = m_f->m_graph->GetMeshDimension();
    for (int i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
    {
        auto exp = m_f->m_exp[0]->GetExp(i);
        int offset = m_f->m_exp[0]->GetPhys_Offset(i);

        int ntot = exp->GetTotPoints();
        Array<OneD, NekDouble> coords[3];
        coords[0] = Array<OneD, NekDouble>(ntot, 0.0);
        coords[1] = Array<OneD, NekDouble>(ntot, 0.0);
        coords[2] = Array<OneD, NekDouble>(ntot, 0.0);
        exp->GetCoords(coords[0], coords[1], coords[2]);

        for (int j = 0; j < ntot; ++j)
        {
            vtkPoints->InsertPoint(offset + j,coords[0][j],coords[1][j],coords[2][j]);
        }

        Array<OneD, int> nquad(meshDim);
        for (int j = 0; j < meshDim; ++j)
        {
            nquad[j] = exp->GetNumPoints(j);
        }

        std::vector<long long> p;
        switch (exp->GetGeom()->GetShapeType())
        {
            case LibUtilities::eQuadrilateral:
                p = QuadrilateralNodes(nquad);
                break;
            case LibUtilities::eTriangle:
                p = TriangleNodes(nquad);
                break;
            default:
                NEKERROR(ErrorUtil::efatal,
                         "VTU output not set up for this shape type.");
                break;
        }

        // Add offset to every value in node list
        std::for_each(p.begin(), p.end(), [&offset](long long& d) { d+=offset;});

        vtkMesh->InsertNextCell(
                GetVtkCellType(exp->GetGeom()->GetShapeType(),
                               exp->GetGeom()->GetGeomFactors()->GetGtype()),
                p.size(), &p[0]);
    }

    vtkMesh->SetPoints(vtkPoints);

    // Save field information to VTU
    for (int i = 0; i < m_f->m_variables.size(); ++i)
    {
        auto &expList = m_f->m_exp[i];

        int nPts = vtkPoints->GetNumberOfPoints();

        vtkNew<vtkDoubleArray> fieldData;
        fieldData->SetNumberOfComponents(1);

        for(int j = 0; j < nPts; ++j)
        {
            fieldData->InsertNextValue(expList->GetPhys()[j]);
        }

        fieldData->SetName(&m_f->m_variables[i][0]);
        vtkMesh->GetPointData()->AddArray(fieldData);
    }

    // Write out the new mesh in XML format (don't support legacy
    // format here as we still have standard OutputVtk.cpp)
    vtkXMLUnstructuredGridWriter *vtkMeshWriter = vtkXMLUnstructuredGridWriter::New();
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

    cout << "Written file: " << filename << endl;
}

std::vector<long long> OutputVtkNew::QuadrilateralNodes(Array<OneD, int> &nquad)
{
    std::vector<long long> p(4);
    // Write vertices
    p[0] = 0;
    p[1] = nquad[0] - 1;
    p[2] = nquad[0] * nquad[1] - 1;
    p[3] = nquad[0] * (nquad[1] - 1);

    // Write edge interior
    for (int j = 0; j < 4; ++j)
    {
        if (j == 0)
        {
            for (int k = 1; k < nquad[0] - 1; ++k)
            {
                p.emplace_back(k);
            }
        }
        else if (j == 1)
        {
            for (int k = 2 * nquad[0] - 1; k < nquad[0] * (nquad[1] - 1); k+=nquad[0])
            {
                p.emplace_back(k);
            }
        }
        else if (j == 2)
        {
            for (int k = nquad[0] * (nquad[1] - 1) + 1; k < nquad[0] * nquad[1] - 1; ++k)
            {
                p.emplace_back(k);
            }
        }
        else if (j == 3)
        {
            for (int k = nquad[0]; k < nquad[0] * (nquad[1] - 1); k+=nquad[0])
            {
                p.emplace_back(k);
            }
        }
    }

    // Write surface interior
    for (int j = 1; j < nquad[0] - 1; ++j)
    {
        for (int k = 1; k < nquad[1] - 1; ++k)
        {
            // Fetch interior nodes from quad->curve
            p.emplace_back(j * nquad[0] + k);
        }
    }

    return p;
};

std::vector<long long> OutputVtkNew::TriangleNodes(Array<OneD, int> &nquad)
{
    std::vector<long long> p(3);
    return p;
};

void OutputVtkNew::OutputFromPts(po::variables_map &vm)
{
    boost::ignore_unused(vm);
    NEKERROR(ErrorUtil::efatal, "OutputVtkNew can't write using only PointData.");
}

void OutputVtkNew::OutputFromData(po::variables_map &vm)
{
    boost::ignore_unused(vm);
    NEKERROR(ErrorUtil::efatal, "OutputVtkNew can't write using only FieldData.");
}

}
}
