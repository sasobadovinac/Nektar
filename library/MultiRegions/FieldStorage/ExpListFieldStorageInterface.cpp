///////////////////////////////////////////////////////////////////////////////
//
// File ExpListFieldStorageInterface.cpp
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
// Description: Provide an decoupling interface to ExpList from FieldStorage.
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList.h>
#include <MultiRegions/FieldStorage/ExpListFieldStorageInterface.h>

namespace Nektar
{
namespace MultiRegions
{
namespace details
{

ExpListFieldStorageInterface::ExpListFieldStorageInterface(
    std::shared_ptr<ExpList> e)
{
    m_e = e;
}

ExpListFieldStorageInterface::ExpListFieldStorageInterface(
    const ExpListFieldStorageInterface &src)
    : m_e(src.m_e)
{
}

int ExpListFieldStorageInterface::GetNpoints()
{
    return m_e->GetNpoints();
}

int ExpListFieldStorageInterface::GetNcoeffs()
{
    return m_e->GetNcoeffs();
}

} // namespace details
} // namespace MultiRegions
} // namespace Nektar
