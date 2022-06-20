#include <iostream>
#include <vector>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

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

template<typename TData,
         StorageType stype,
         DataLayout order = eField>
class FieldStorage{

    FieldStorage(MultiRegions::ExpListSharedPtr exp)
        : m_exp(exp), m_numVariables(1)
    {
        if (stype == ePhys)
        {
            m_storage = Array<OneD, TData>(m_exp->GetNpoints());
        }
        else if (stype == eCoeff)
        {
            m_storage = Array<OneD, TData>(m_exp->GetNcoeffs());
        }
    }

    FieldStorage(const FieldStorage& F)
        : m_exp(F.m_exp), m_sType(F.m_sType), 
          m_storage(Array<One, TData>(F.m_storage->size(), F.m_storage)), 
          m_numVariables(F.m_numVariables)
    {
    }

    FieldStorage(const FieldStorage&& F)
    {
    }

    ~FieldStorage()
    {
        // nothing to do... yet...
    }

    Array<OneD, TData>& UpdateData()
    {
        return m_storage;
    }

    const Array<OneD, const TData> GetData() const
    {
        return m_storage;
    }
    
private:
    ExpListSharedPtr m_exp;
    enum StorageType               m_sType;
    Array<OneD, TData>  m_storage;
//    std::vector<size_t>            m_offsets;      // Offset of element i in the array
    int m_numVariables;
//    int num_elements;
//    int std::array<int, num_variables> dofs;
};

}
}