///////////////////////////////////////////////////////////////////////////////
//
// File FieldStorage.hpp
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
// Description: Container for solution fields.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_FIELDSTORAGE_H
#define NEKTAR_LIBS_MULTIREGIONS_FIELDSTORAGE_H

#include <iostream>
#include <vector>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/FieldStorage/ExpListFieldStorageInterface.h>

namespace Nektar
{
namespace MultiRegions
{

enum DataLayout
{
    eField, // data for variable 1 followed by variable 2 etc
    eElement,
    eVariable
};

enum StorageType
{
    ePhys,
    eCoeff
};

class ExpList;

template <typename TData, StorageType stype, DataLayout order = eField>
class FieldStorage
{
public:
    FieldStorage(std::shared_ptr<ExpList> exp) : m_expIF(exp), m_numVariables(1)
    {
        boost::ignore_unused(exp);
        if (stype == ePhys)
        {
            m_storage = Array<OneD, TData>(m_expIF->GetNpoints());
        }
        else if (stype == eCoeff)
        {
            m_storage = Array<OneD, TData>(m_expIF->GetNcoeffs());
        }
    }

    FieldStorage(const FieldStorage &F)
        : m_expIF(F.m_expIF), m_sType(F.m_sType),
          m_storage(Array<OneD, TData>(F.m_storage->size(), F.m_storage)),
          m_numVariables(F.m_numVariables)
    {
    }

    ~FieldStorage()
    {
        // nothing to do... yet...
    }

    Array<OneD, TData> &UpdateData()
    {
        return m_storage;
    }

    const Array<OneD, const TData> GetData() const
    {
        return m_storage;
    }

private:
    std::shared_ptr<ExpListFieldStorageInterface> m_expIF;
    enum StorageType m_sType;
    Array<OneD, TData> m_storage;
    //    std::vector<size_t>            m_offsets;      // Offset of element i
    //    in the array
    int m_numVariables;
    //    int num_elements;
    //    int std::array<int, num_variables> dofs;
};

} // namespace MultiRegions
} // namespace Nektar
#endif
