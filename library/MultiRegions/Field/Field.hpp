///////////////////////////////////////////////////////////////////////////////
//
// File Field.hpp
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
#include <MultiRegions/Field/ExpListFieldInterface.h>

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

template <typename TData, StorageType TStype>
class Field
{
public:

    // default contructor
    Field(): m_numVariables(1)
    {
    };

    Field(std::shared_ptr<MultiRegions::ExpList> exp,
                 TData defval = 0, int nvar = 1, DataLayout Order = eField) :
        m_numVariables(nvar),
        m_dataOrder(Order)
    {
        m_expIF = std::make_shared<MultiRegions::details::
                                   ExpListFieldInterface>(exp);

        m_storage = Array<OneD, Array<OneD, NekDouble> >(nvar); 
        int varsize; 
        if (TStype == ePhys)
        {
            varsize = m_expIF->GetNpoints(); 
        }
        else if (TStype == eCoeff)
        {
            varsize = m_expIF->GetNcoeffs(); 
        }

        m_storage[m_numVariables-1] =Array<OneD, TData>(varsize*m_numVariables,defval);

        for(int i = m_numVariables-1; i > 0; --i)
        {
            m_storage[i-1] = m_storage[i] + varsize; 
        }
    }

    Field(const Field &F)
        : m_expIF(F.m_expIF),
          m_storage(Array<OneD, Array<OneD, TData>>(F.m_numVariables)),
          m_numVariables(F.m_numVariables)
    {
        // Fecalre field in reverse order so that data could be
        // accessed in one block but m_storatge[0] will return correct
        // array length for one variable since this is often used as a
        // sizing variable in code. Probably there is a better
        // solution but hack for now.
        int varsize = F.m_storage[0].size();
        m_storage[m_numVariables-1] =Array<OneD, TData>(varsize*m_numVariables,F.m_storage[0]);

        for(int i = m_numVariables-1; i > 0; --i)
        {
            m_storage[i-1] = m_storage[i] + varsize; 
        }
    }
    
    ~Field()
    {
        // nothing to do... yet...
    }

    const Array<OneD, const TData> GetData(int varid=0) const
    {
        ASSERTL1(varid < m_numVariables, "variable id (varid) is out of range");
        return m_storage[varid];
    }

    Array<OneD, TData> &UpdateData(int varid = 0)
    {
        ASSERTL1(varid < m_numVariables, "variable id (varid) is out of range");
        return m_storage[varid];
    }

    const Array<OneD, const TData> GetArray1D(int varid=0) const
    {
        ASSERTL1(varid < m_numVariables, "variable id (varid) is out of range");
        ASSERTL1(m_dataOrder == eField,"Data must be in Array OneD format for"
                 " this method");
        return m_storage[varid];
    }
    
    Array<OneD, TData> &UpdateArray1D(int varid = 0)
    {
        ASSERTL1(varid < m_numVariables, "variable id (varid) is out of range");
        ASSERTL1(m_dataOrder == eField,"Data must be in Array OneD format for"
                 " this method");
        return m_storage[varid];
    }

    std::shared_ptr<MultiRegions::ExpList> GetExpList() const
    {
        return m_expIF->GetExpList();
    }

    inline int GetNumVariables() const 
    {
        return m_numVariables; 
    }
protected: 
    /// interface to allow access to ExpList 
    std::shared_ptr<MultiRegions::details::ExpListFieldInterface> m_expIF;
    /// native storage in field
    Array<OneD, Array<OneD, TData> >   m_storage;
    
    /// number of variables in storage 
    int m_numVariables;

    /// Data ordering format 
    DataLayout m_dataOrder; 
};

static Field<NekDouble,ePhys> ZeroFieldPhys; 
static Field<NekDouble,eCoeff> ZerolFieldCoef; 
    
} // namespace Nektar
#endif
