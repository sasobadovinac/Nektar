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
//  Description: Vtk output module
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_OUTPUTVTK
#define UTILITIES_NEKMESH_OUTPUTVTK

#include <NekMesh/Module/Module.h>

namespace Nektar
{
namespace NekMesh
{

/// Converter for Gmsh files.
class OutputVtk : public NekMesh::OutputModule
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(NekMesh::MeshSharedPtr m)
    {
        return MemoryManager<OutputVtk>::AllocateSharedPtr(m);
    }
    static NekMesh::ModuleKey className;

    OutputVtk(NekMesh::MeshSharedPtr m);
    virtual ~OutputVtk();

    /// Write mesh to output file.
    virtual void Process();
};
}
}

#endif
