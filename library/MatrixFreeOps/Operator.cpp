///////////////////////////////////////////////////////////////////////////////
//
// File: Operator.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include "Operator.hpp"

namespace Nektar
{
namespace MatrixFree
{

OperatorFactory &GetOperatorFactory()
{
    static OperatorFactory tmp;
    return tmp;
}

std::string GetOpstring(LibUtilities::ShapeType shape, bool deformed)
{

    std::string op_string = "_";

    if (shape == LibUtilities::eSegment)
    {
        op_string += "Seg";
    }
    else if (shape == LibUtilities::eTriangle)
    {
        op_string += "Tri";
    }
    else if (shape == LibUtilities::eQuadrilateral)
    {
        op_string += "Quad";
    }
    else if (shape == LibUtilities::eTetrahedron)
    {
        op_string += "Tet";
    }
    else if (shape == LibUtilities::ePyramid)
    {
        op_string += "Pyr";
    }
    else if (shape == LibUtilities::ePrism)
    {
        op_string += "Prism";
    }
    else if (shape == LibUtilities::eHexahedron)
    {
        op_string += "Hex";
    }

    if (deformed)
    {
        op_string += "_Deformed";
    }
    else
    {
        op_string += "_Regular";
    }

    return op_string;
}

} // namespace MatrixFree
} // namespace Nektar
