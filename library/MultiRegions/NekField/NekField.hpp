///////////////////////////////////////////////////////////////////////////////
//
// File NekField.hpp
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
#include <MultiRegions/NekField/ExpListNekFieldInterface.h>

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
class NekField
{
public:

    // default contructor
    NekField(): m_numVariables(1)
    {
    };

    NekField(std::shared_ptr<MultiRegions::ExpList> exp,
                 TData defval = 0, int nvar = 1, DataLayout Order = eField) :
        m_numVariables(nvar),
        m_dataOrder(Order)
    {
        m_expIF.push_back(std::make_shared<MultiRegions::details::
                          ExpListNekFieldInterface>(exp));
        // soft copy passed expansion into nvar 
        for(int i = 1; i < nvar; ++i)
        {
            m_expIF.push_back(m_expIF[0]);
        }

        m_storage = Array<OneD, Array<OneD, NekDouble> >(nvar); 
        int varsize; 
        if (TStype == ePhys)
        {
            varsize = m_expIF[0]->GetNpoints(); 
        }
        else if (TStype == eCoeff)
        {
            varsize = m_expIF[0]->GetNcoeffs(); 
        }

        m_storage[m_numVariables-1] =Array<OneD, TData>(varsize*m_numVariables,defval);

        for(int i = m_numVariables-1; i > 0; --i)
        {
            m_storage[i-1] = m_storage[i] + varsize; 
        }
    }

    NekField(Array<OneD, std::shared_ptr<MultiRegions::ExpList>> exp,
                 TData defval = 0, DataLayout Order = eField) :
        m_numVariables(exp.size()),
        m_dataOrder(Order)
    {
        for(int i = 0; i < m_numVariables; ++i)
        {
            m_expIF.push_back(std::make_shared<MultiRegions::details::
                          ExpListNekFieldInterface>(exp[i]));
        }

        m_storage = Array<OneD, Array<OneD, NekDouble> >(m_numVariables); 
        int varsize; 
        if (TStype == ePhys)
        {
            varsize = m_expIF[0]->GetNpoints(); 
        }
        else if (TStype == eCoeff)
        {
            varsize = m_expIF[0]->GetNcoeffs(); 
        }

        m_storage[m_numVariables-1] =Array<OneD, TData>(varsize*m_numVariables,defval);

        for(int i = m_numVariables-1; i > 0; --i)
        {
            m_storage[i-1] = m_storage[i] + varsize; 
        }
    }
    
    NekField(const NekField &F)
        : m_expIF(F.m_expIF),
          m_storage(Array<OneD, Array<OneD, TData>>(F.m_numVariables)),
          m_numVariables(F.m_numVariables),
          m_dataOrder(F.m_dataOrder)
    {
        // Fecalre field in reverse order so that data could be
        // accessed in one block but m_storatge[0] will return correct
        // array length for one variable since this is often used as a
        // sizing variable in code. Probably there is a better
        // solution but hack for now.
        int varsize = F.m_storage[0].size();
        m_storage[m_numVariables-1] =Array<OneD, TData>(varsize*m_numVariables,F.m_storage.GetData());

        for(int i = m_numVariables-1; i > 0; --i)
        {
            m_storage[i-1] = m_storage[i] + varsize; 
        }
    }

    /// \brief Creates a reference to rhs.
    NekField &operator=(const NekField &rhs)
    {
        m_expIF        = rhs.m_expIF; 
        m_storage      = rhs.m_storage; 
        m_numVariables = rhs.m_numVariables; 
        m_dataOrder    = rhs.m_dataOrder; 

    }
    
    ~NekField()
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

    std::vector<std::shared_ptr<MultiRegions::ExpList>> GetExpList() const
    {
        return m_expIF[0]->GetExpList();
    }

    std::shared_ptr<MultiRegions::ExpList> GetExpList(int varid) const
    {
        return m_expIF[varid]->GetExpList();
    }

    inline int GetNumVariables() const 
    {
        return m_numVariables; 
    }
protected: 
    /// interface to allow access to ExpList 
    std::vector<std::shared_ptr<MultiRegions::details::ExpListNekFieldInterface>> m_expIF;
    
    /// native storage in field
    Array<OneD, Array<OneD, TData> >   m_storage;
    
    /// number of variables in storage 
    int m_numVariables;

    /// Data ordering format 
    DataLayout m_dataOrder; 
};

static NekField<NekDouble,ePhys> ZeroNekFieldPhys; 
static NekField<NekDouble,eCoeff> ZerolNekField; 

typedef std::shared_ptr<NekField<NekDouble,eCoeff>> NekFieldCoeffSharedPtr;
typedef std::shared_ptr<NekField<NekDouble,ePhys >> NekFieldPhysSharedPtr; 
} // Namespace
#endif
