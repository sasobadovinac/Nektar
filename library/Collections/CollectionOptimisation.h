///////////////////////////////////////////////////////////////////////////////
//
// File: CollectionOptimisation.h
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
// Description: Collection optimization definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBRARY_COLLECTIONS_COLLECTIONOPTIMISATION_H
#define NEKTAR_LIBRARY_COLLECTIONS_COLLECTIONOPTIMISATION_H

#include <Collections/Collection.h>
#include <Collections/CollectionsDeclspec.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <tinyxml.h>

namespace Nektar
{
namespace Collections
{

class OpImpTimingKey
{
public:
    /// Constructor
    OpImpTimingKey(StdRegions::StdExpansionSharedPtr pExp, int ngeoms,
                   int nbases)
        : m_exp(pExp), m_ngeoms(ngeoms), m_nbasis(nbases)
    {
    }

    /// Destructor
    ~OpImpTimingKey(void)
    {
    }

    LibUtilities::ShapeType GetShapeType(void) const
    {
        return m_exp->DetShapeType();
    }

    int GetExpOrder() const
    {
        return m_exp->GetBasis(0)->GetNumModes();
    }

    int GetNGeoms() const
    {
        return m_ngeoms;
    }

    bool operator<(const OpImpTimingKey &rhs) const
    {

        if (m_nbasis < rhs.m_nbasis)
        {
            return true;
        }

        if (m_nbasis > rhs.m_nbasis)
        {
            return false;
        }

        for (int i = 0; i < m_nbasis; ++i)
        {
            if (m_exp->GetBasis(i)->GetBasisKey() !=
                rhs.m_exp->GetBasis(i)->GetBasisKey())
            {
                return (m_exp->GetBasis(i)->GetBasisKey() <
                        rhs.m_exp->GetBasis(i)->GetBasisKey());
            }
        }

        if ((m_ngeoms < 100) && (rhs.m_ngeoms < 100))
        {
            if (m_ngeoms < rhs.m_ngeoms)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        return false;
    }

    StdRegions::StdExpansionSharedPtr m_exp;
    int m_ngeoms;
    int m_nbasis;

private:
};

class CollectionOptimisation
{
public:
    // Constuctor
    COLLECTIONS_EXPORT CollectionOptimisation(
        LibUtilities::SessionReaderSharedPtr pSession, const int shapedim,
        ImplementationType defaultType = eStdMat);

    ~CollectionOptimisation(){};

    ImplementationType GetDefaultImplementationType()
    {
        return m_defaultType;
    }

    unsigned int GetMaxCollectionSize()
    {
        return m_maxCollSize;
    }

    bool IsUsingAutotuning()
    {
        return m_autotune;
    }

    /// Get Operator Implementation Map from XMl or using default;
    COLLECTIONS_EXPORT OperatorImpMap
    GetOperatorImpMap(StdRegions::StdExpansionSharedPtr pExp);

    // Get Map by doing autotuning testing.
    COLLECTIONS_EXPORT OperatorImpMap
    SetWithTimings(std::vector<StdRegions::StdExpansionSharedPtr> pGeom,
                   OperatorImpMap &impTypes, bool verbose = true);

    // Wite out autotuning testing to file
    COLLECTIONS_EXPORT void UpdateOptFile(std::string sessName,
                                          LibUtilities::CommSharedPtr &comm);

private:
    typedef std::pair<LibUtilities::ShapeType, int> ElmtOrder;
    typedef std::map<OperatorType, std::map<ElmtOrder, ImplementationType>>
        GlobalOpMap;
    static std::map<OpImpTimingKey, OperatorImpMap> m_opImpMap;
    GlobalOpMap m_global;
    bool m_autotune;
    ImplementationType m_defaultType;
    unsigned int m_maxCollSize;
    unsigned int m_shapeDim;

    void ReadCollOps(TiXmlElement *xmlCol, GlobalOpMap &global, bool verbose);
};

} // namespace Collections
} // namespace Nektar
#endif
