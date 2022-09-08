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
#include <LibUtilities/BasicUtils/Vmath.hpp>
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
    NekField(): m_numVariables(0), m_varSize(0)
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

        if (TStype == ePhys)
        {
            m_varSize = m_expIF[0]->GetNpoints(); 
        }
        else if (TStype == eCoeff)
        {
            m_varSize = m_expIF[0]->GetNcoeffs(); 
        }

        m_storage.resize(m_numVariables*m_varSize, defval);
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

        if (TStype == ePhys)
        {
            m_varSize = m_expIF[0]->GetNpoints(); 
        }
        else if (TStype == eCoeff)
        {
            m_varSize = m_expIF[0]->GetNcoeffs(); 
        }

        m_storage.resize(m_numVariables*m_varSize,defval);
    }

    NekField(std::vector<std::shared_ptr<MultiRegions::ExpList>> exp,
                 TData defval = 0, DataLayout Order = eField) :
        m_numVariables(exp.size()),
        m_dataOrder(Order)
    {
        for(int i = 0; i < m_numVariables; ++i)
        {
            m_expIF.push_back(std::make_shared<MultiRegions::details::
                          ExpListNekFieldInterface>(exp[i]));
        }

        if (TStype == ePhys)
        {
            m_varSize = m_expIF[0]->GetNpoints(); 
        }
        else if (TStype == eCoeff)
        {
            m_varSize = m_expIF[0]->GetNcoeffs(); 
        }

        m_storage.resize(m_numVariables*m_varSize, defval);
    }
    
    NekField(const NekField &F)
        : m_expIF(F.m_expIF),
          m_storage(F.m_storage),
          m_numVariables(F.m_numVariables),
          m_varSize(F.m_varSize),
          m_dataOrder(F.m_dataOrder)
    {
    }

    /// \brief Creates a reference to rhs.
    NekField &operator=(const NekField &rhs)
    {
        m_expIF        = rhs.m_expIF; 
        m_storage      = rhs.m_storage; 
        m_numVariables = rhs.m_numVariables; 
        m_varSize      = rhs.m_varSize; 
        m_dataOrder    = rhs.m_dataOrder; 
    }
    
    ~NekField()
    {
        // nothing to do... yet...
    }

    void AddVariable(std::shared_ptr<MultiRegions::ExpList> exp,
                           TData defval = 0)
    {

        m_expIF.push_back(std::make_shared<MultiRegions::details::
                          ExpListNekFieldInterface>(exp));

        m_numVariables +=1;

        // set up varSize in case initialised from default constructor. 
        if(m_varSize == 0)
        {
            if (TStype == ePhys)
            {
                m_varSize = m_expIF[0]->GetNpoints(); 
            }
            else if (TStype == eCoeff)
            {
                m_varSize = m_expIF[0]->GetNcoeffs(); 
            }
        }

        m_storage.resize(m_varSize*m_numVariables, defval);
    }


    void AddVariable(std::vector<std::shared_ptr<MultiRegions::ExpList>> exp,
                     TData defval = 0)
    {

        int nvar = exp.size();
        for(int i = 0; i < nvar; ++i)
        {
            m_expIF.push_back(std::make_shared<MultiRegions::details::
                              ExpListNekFieldInterface>(exp[i]));
        }

        m_numVariables += nvar;


        // set up varSize in case initialised from default constructor. 
        if(m_varSize == 0)
        {
            if (TStype == ePhys)
            {
                m_varSize = m_expIF[0]->GetNpoints(); 
            }
            else if (TStype == eCoeff)
            {
                m_varSize = m_expIF[0]->GetNcoeffs(); 
            }
        }
        
        m_storage.resize(m_varSize*m_numVariables, defval);
    }
    
    const Array<OneD, const TData> GetData(int varid=0) const
    {
        ASSERTL1(varid < m_numVariables, "variable id (varid) is out of range");
        m_array1D = Array<OneD, TData>(m_varSize,
                                       &m_storage[varid*m_varSize],true);
        return m_array1D; 
    }

    Array<OneD, TData> &UpdateData(int varid = 0)
    {
        ASSERTL1(varid < m_numVariables, "variable id (varid) is out of range");
        m_array1D = Array<OneD, TData>(m_varSize,
                                       &m_storage[varid*m_varSize],true);
        return m_array1D; 
    }

    const Array<OneD, const TData> GetArray1D(int varid=0) const
    {
        ASSERTL1(varid < m_numVariables, "variable id (varid) is out of range");
        ASSERTL1(m_dataOrder == eField,"Data must be in Array OneD format for"
                 " this method");
        
        return Array<OneD, TData> (m_varSize, &m_storage[varid*m_varSize]);
    }
    
    Array<OneD, TData> UpdateArray1D(int varid = 0)
    {
        ASSERTL1(varid < m_numVariables, "variable id (varid) is out of range");
        ASSERTL1(m_dataOrder == eField,"Data must be in Array OneD format for"
                 " this method");

        return Array<OneD, TData>(m_varSize,
                                       &m_storage[varid*m_varSize],true);
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

    inline int GetVarSize() const 
    {
        return m_varSize; 
    }
protected: 
    /// interface to allow access to ExpList 
    std::vector<std::shared_ptr<MultiRegions::details::ExpListNekFieldInterface>> m_expIF;
    
    /// native storage in field
    std::vector<NekDouble>   m_storage;

    Array<OneD, TData> m_array1D;
    
    /// number of variables in storage 
    int m_numVariables;

    /// size of storage  per variable
    int m_varSize;

    /// Data ordering format 
    DataLayout m_dataOrder; 
};

static NekField<NekDouble,ePhys> ZeroNekFieldPhys; 
static NekField<NekDouble,eCoeff> ZerolNekField; 

typedef std::shared_ptr<NekField<NekDouble,eCoeff>> NekFieldCoeffSharedPtr;
typedef std::shared_ptr<NekField<NekDouble,ePhys >> NekFieldPhysSharedPtr; 
} // Namespace
#endif
