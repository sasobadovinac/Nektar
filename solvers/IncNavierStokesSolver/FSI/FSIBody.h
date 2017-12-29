///////////////////////////////////////////////////////////////////////////////
//
// File: FSIBody.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Abstract base class for fluid-structure interaction bodies.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_FSI_FSIBODY
#define NEKTAR_SOLVERS_FSI_FSIBODY

#include <string>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <GlobalMapping/GlobalMappingDeclspec.h>

namespace Nektar
{

//  Forward declaration
class FSIBody;

/// A shared pointer to a FSIBody object
GLOBAL_MAPPING_EXPORT typedef std::shared_ptr<FSIBody> FSIBodySharedPtr;

/// Declaration of the FSIBody factory
typedef LibUtilities::NekFactory<std::string, FSIBody,
        const LibUtilities::SessionReaderSharedPtr&,
        const Array<OneD, MultiRegions::ExpListSharedPtr>&,
        const std::map<std::string, std::string>&> FSIBodyFactory;

/// Declaration of the FSIBody factory singleton
GLOBAL_MAPPING_EXPORT FSIBodyFactory& GetFSIBodyFactory();

/**
 * @class FSIBody
 * @brief Base class for updating the mapping in FSI problems
 */
class FSIBody
{
public:
    /// @brief Destructor
    GLOBAL_MAPPING_EXPORT virtual ~FSIBody() {}

    /// @brief Initialise the FSIBody object
    void InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const std::map<std::string, std::string>          &pParams);

    ///
    GLOBAL_MAPPING_EXPORT void Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pDisplFields,
        const NekDouble                                      &time);

protected:
    /// Session reader
    LibUtilities::SessionReaderSharedPtr m_session;
    /// Determines if a given Boundary Region is part of this body
    std::vector<int>                     m_boundaryRegionIsInList;

    /// @brief Constructor
    FSIBody(
        const LibUtilities::SessionReaderSharedPtr           &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pFields);

    // Virtual functions
    virtual void v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pFields,
        const std::map<std::string, std::string>             &pParams);

    virtual void v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pDisplFields,
        const NekDouble                                      &time) = 0;

};

inline void FSIBody::InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>     &pFields,
        const std::map<std::string, std::string>              &pParams)
{
    v_InitObject(pFields, pParams);
}

inline void FSIBody::Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pDisplFields,
        const NekDouble                                      &time)
{
    v_Apply(pFields, pDisplFields, time);
}

}

#endif
