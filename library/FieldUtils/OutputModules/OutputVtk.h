////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputVtk.h
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
//  Description: Vtk output module using the VTK library
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_OUTPUTVTK
#define FIELDUTILS_OUTPUTVTK

#include "OutputVtkBase.h"
#include <tinyxml.h>

#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

namespace Nektar
{
namespace FieldUtils
{

/// Converter from fld to vtk.
class OutputVtk final : public OutputVtkBase
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(const FieldSharedPtr &f)
    {
        return MemoryManager<OutputVtk>::AllocateSharedPtr(f);
    }

    static ModuleKey m_className;

    explicit OutputVtk(FieldSharedPtr f);

    ~OutputVtk() final = default;

protected:
    virtual std::string v_GetModuleName() override final
    {
        return "OutputVtk";
    }

    /// Write from pts to output file.
    virtual void v_OutputFromPts(po::variables_map &vm) override final;

    /// Write from m_exp to output file.
    virtual void v_OutputFromExp(po::variables_map &vm) override final;

    /// Write from data to output file.
    virtual void v_OutputFromData(po::variables_map &vm) override final;

    /// Cache file for unstructured grid VTK mesh data
    vtkSmartPointer<vtkUnstructuredGrid> m_vtkMesh;

    /// Number of planes if homogeneous
    int m_numPlanes = 1;

    /// Flag if extra plane in case of fourier expansion in homogeneous dir
    bool m_extraPlane = false;

    /// Flag if mesh has been cached
    bool m_cachedMesh = false;

private:
    /// Prepare high order Lagrange VTK output
    vtkSmartPointer<vtkUnstructuredGrid> OutputFromExpHighOrder(
        po::variables_map &vm);

    /// Add field data to high order Lagrange VTK output
    void AddFieldDataToVTKHighOrder(
        po::variables_map &vm, std::string &filename,
        vtkSmartPointer<vtkUnstructuredGrid> &vtkMesh);

    /// Prepare low order VTK output
    vtkSmartPointer<vtkUnstructuredGrid> OutputFromExpLowOrder();

    /// Add field data to low order VTK output
    void AddFieldDataToVTKLowOrder(
        po::variables_map &vm, std::string &filename,
        vtkSmartPointer<vtkUnstructuredGrid> &vtkMesh);

    /// Prepare low order multi-block VTK output & add field data
    void OutputFromExpLowOrderMultiBlock(po::variables_map &vm,
                                         std::string &filename);

    /// Write VTK file using @param vtkMesh
    void WriteVTK(vtkDataObject *vtkMesh, std::string &filename,
                  po::variables_map &vm);

    /// Write the parallel .pvtu file
    void WritePVtu(po::variables_map &vm);
};
} // namespace FieldUtils
} // namespace Nektar

#endif
