///////////////////////////////////////////////////////////////////////////////
//
// File: FSICoupler.h
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
// Description: Abstract base class for fluid-structure interaction coupler.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_FSI_FSICOUPLER
#define NEKTAR_SOLVERS_FSI_FSICOUPLER

#include <string>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationWrapper.h>
#include <MultiRegions/ExpList.h>
#include <GlobalMapping/GlobalMappingDeclspec.h>
#include <GlobalMapping/Mapping.h>
#include <IncNavierStokesSolver/FSI/FSIBody.h>

namespace Nektar
{

//  Forward declaration
class FSICoupler;

/// A shared pointer to a FSICoupler object
GLOBAL_MAPPING_EXPORT typedef std::shared_ptr<FSICoupler> FSICouplerSharedPtr;

/// Declaration of the FSICoupler factory
typedef LibUtilities::NekFactory<std::string, FSICoupler,
        const LibUtilities::SessionReaderSharedPtr&,
        const Array<OneD, MultiRegions::ExpListSharedPtr>&,
              TiXmlElement*> FSICouplerFactory;

/// Declaration of the FSICoupler factory singleton
GLOBAL_MAPPING_EXPORT FSICouplerFactory& GetFSICouplerFactory();

/**
 * @class FSICoupler
 * @brief Base class for updating the mapping in FSI problems
 */
class FSICoupler
{
public:
    /// @brief Destructor
    GLOBAL_MAPPING_EXPORT virtual ~FSICoupler() {}

    /// @brief Initialise the FSICoupler object
    void InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
              TiXmlElement* pFSI);

    ///
    GLOBAL_MAPPING_EXPORT void Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        GlobalMapping::MappingSharedPtr                   &pMapping,
        const NekDouble                                   &time);

protected:
    /// Session reader
    LibUtilities::SessionReaderSharedPtr        m_session;
    /// Vector of moving bodies
    std::vector<FSIBodySharedPtr>               m_bodies;
    /// Explist for the displacement of the coordinates
    Array<OneD, MultiRegions::ExpListSharedPtr> m_displFields;

    /// Backward-difference coefficients for evaluating coordinates velocity
    static NekDouble                            BDF_Alpha_Coeffs[3][3];
    static NekDouble                            BDF_Gamma0_Coeffs[3];

    // Time step
    NekDouble                                   m_timestep;
    // Dimensions
    int                                         m_expDim;
    int                                         m_spaceDim;
    // Time integration order
    int                                         m_intSteps;

    // Coordinates of the transformed system
    Array<OneD, Array<OneD, NekDouble>>         m_coords;
    // Coordinates of the mesh
    Array<OneD, Array<OneD, NekDouble>>         m_meshCoords;
    // Coordinates velocity
    Array<OneD, Array<OneD, NekDouble>>         m_coordsVel;
    // Previous coordinates for evaluating m_coordsVel
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> m_oldCoords;

    /// @brief Constructor
    FSICoupler(
        const LibUtilities::SessionReaderSharedPtr&          pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields);

    ///
    void CreateDisplacementFields(
        const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields);

    /// 
    void CalculateDisplacement();

    ///
    void UpdateCoordinates();

    ///
    void CalculateCoordVel();

    ///
    void ReadBodies(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
              TiXmlElement* pFSI);

    // Virtual functions
    virtual void v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields,
              TiXmlElement* pFSI);

    virtual void v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        GlobalMapping::MappingSharedPtr                   &pMapping,
        const NekDouble                                   &time);

    virtual void v_CalculateDisplacement() = 0;
};

inline void FSICoupler::InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>&     pFields,
              TiXmlElement* pFSI)
{
    v_InitObject(pFields, pFSI);
}

inline void FSICoupler::Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        GlobalMapping::MappingSharedPtr                   &pMapping,
        const NekDouble                                   &time)
{
    v_Apply(pFields, pMapping, time);
}

inline void FSICoupler::CalculateDisplacement()
{
    v_CalculateDisplacement();
}

}

#endif
