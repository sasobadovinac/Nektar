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

template <typename T_Data, StorageType T_Stype, DataLayout T_Order = eField>
class FieldStorage
{
public:
    FieldStorage(std::shared_ptr<MultiRegions::ExpList> exp, T_Data defval = 0) :
        m_numVariables(1)
    {
        m_expIF = std::make_shared<MultiRegions::details::
                                   ExpListFieldStorageInterface>(exp);

        if (T_Stype == ePhys)
        {
            m_storage = Array<OneD, T_Data>(m_expIF->GetNpoints(),defval);
        }
        else if (T_Stype == eCoeff)
        {
            m_storage = Array<OneD, T_Data>(m_expIF->GetNcoeffs(),defval);
        }
    }

    FieldStorage(const FieldStorage &F)
        : m_expIF(F.m_expIF),
          m_storage(Array<OneD, T_Data>(F.m_storage->size(), F.m_storage)),
          m_numVariables(F.m_numVariables)
    {
    }
    
    ~FieldStorage()
    {
        // nothing to do... yet...
    }

    const Array<OneD, const T_Data> GetData() const
    {
        return m_storage;
    }
    
    Array<OneD, T_Data> &UpdateData()
    {
        return m_storage;
    }
  
    std::shared_ptr<MultiRegions::ExpList> GetExpList() const
    {
        return m_expIF->GetExpList();
    }

protected: 
    /// interface to allow access to ExpList 
    std::shared_ptr<MultiRegions::details::ExpListFieldStorageInterface> m_expIF;
    
    /// Storage type 
        //enum StorageType            m_sType;
    
    /// native storage in field
    Array<OneD, T_Data>          m_storage;
    
    /// Offset of element i in the array
    //    std::vector<size_t>   m_offsets; 
    
    /// number of variables in storage 
    int m_numVariables;

    //    int num_elements;
    //    int std::array<int, num_variables> dofs;
};

} // namespace Nektar
#endif
