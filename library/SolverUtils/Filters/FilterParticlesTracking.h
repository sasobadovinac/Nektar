///////////////////////////////////////////////////////////////////////////////
//
// File   FilterParticlesTracking.h
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
// Description: Particle tracking filter.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERPARTICLESTRACKING_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERPARTICLESTRACKING_H

#include <SolverUtils/Filters/Filter.h>

namespace Nektar
{
namespace SolverUtils
{

/// Struct holding information of a particle for particle tracking
struct Particle
{
    /// Constructor
    Particle(
        int                     dim,
        int                     numFields,
        int                     intOrder,
        Array<OneD, NekDouble>  newCoord)
        : m_dim(dim), m_eId(-1), m_used(true)
    {
        // Initialise arrays
        m_newCoord          = Array<OneD, NekDouble> (3, 0.0);
        m_oldCoord          = Array<OneD, NekDouble> (3, 0.0);
        m_locCoord          = Array<OneD, NekDouble> (3, 0.0);
        m_fluidVelocity     = Array<OneD, NekDouble> (m_dim, 0.0);
        m_fields            = Array<OneD, NekDouble> (numFields, 0.0);
        // Force, Velocity, Torque and Angular velocity with intOrder entries
        m_force            = Array<OneD, Array<OneD,NekDouble>> (intOrder);
        m_particleVelocity = Array<OneD, Array<OneD,NekDouble>> (intOrder);
        m_torque           = Array<OneD, Array<OneD,NekDouble>> (intOrder);
        m_angularVelocity  = Array<OneD, Array<OneD,NekDouble>> (intOrder);
        for (int i = 0; i < intOrder; ++i)
        {
            m_force[i]            = Array<OneD, NekDouble> (m_dim, 0.0);
            m_particleVelocity[i] = Array<OneD, NekDouble> (m_dim, 0.0);
            m_torque[i]           = Array<OneD, NekDouble> (3, 0.0);
            m_angularVelocity[i]  = Array<OneD, NekDouble> (3, 0.0);
        }
        // gradient array  
        m_grad = Array<OneD,Array<OneD,NekDouble> >(m_dim*m_dim);
        for (int i = 0; i < m_dim*m_dim; ++i)
        {
            m_grad[i] = Array<OneD, NekDouble> (m_dim, 0.0);
        }

        // Store coordinates
        SetCoord(newCoord);

        // Obtain a unique id
        SetNewId();
        // Times particles counter
        m_advanceCalls = 0;
    }

    /// Function to assign a new id to the particle
    void SetNewId()
    {
        m_id = m_nextId++;
    }

    /// Function to assign a new id to the particle
    void SetCoord(Array<OneD, NekDouble>  newCoord)
    {
        m_newCoord[0] = newCoord[0];
        m_newCoord[1] = newCoord[1];
        m_newCoord[2] = newCoord[2];
    }

    /// Static counter for numbering particles
    static int                      m_nextId;
    /// Id of the particle
    int                             m_id;
    /// Spatial dimension
    int                             m_dim;
    /// Global coordinate
    Array<OneD, NekDouble>          m_newCoord;
    /// Previous coordinates
    Array<OneD, NekDouble>          m_oldCoord;
    /// Coordinate in the standard element
    Array<OneD, NekDouble>          m_locCoord;
    /// Velocity of the particle
    Array<OneD, Array<OneD, NekDouble>> m_particleVelocity;
    /// Angular Velocity of the particle
    Array<OneD, Array<OneD, NekDouble>> m_angularVelocity;
    /// Fluid velocity
    Array<OneD, NekDouble>          m_fluidVelocity;
    /// Id of the element currently holding the particle
    int                             m_eId;
    /// Flag marking if particle is still in use (it may have left the domain)
    bool                            m_used;
    /// Simulation solution at the particle location
    Array<OneD, NekDouble>          m_fields;
    /// Force acting on the particle
    Array<OneD, Array<OneD, NekDouble>> m_force;
    /// Torque acting on the particle
    Array<OneD, Array<OneD, NekDouble>> m_torque;
    /// Gradiend acting on the particle
    Array<OneD, Array<OneD, NekDouble>> m_grad;
    /// Counter of times particles were advanced
    unsigned int                    m_advanceCalls;
};



class FilterParticlesTracking : public Filter
{
public:
    friend class MemoryManager<FilterParticlesTracking>;

    /// Creates an instance of this class
    static FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem>      &pEquation,
        const std::map<std::string, std::string> &pParams)
    {
        FilterSharedPtr p = MemoryManager<FilterParticlesTracking>::
                           AllocateSharedPtr(pSession, pEquation, pParams);
        return p;
    }

    ///Name of the class
    static std::string className;

    SOLVER_UTILS_EXPORT FilterParticlesTracking(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem>      &pEquation,
        const ParamMap                             &pParams);
    SOLVER_UTILS_EXPORT virtual ~FilterParticlesTracking();

protected:
    SOLVER_UTILS_EXPORT virtual void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
       const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual bool v_IsTimeDependent();

private:
    /// Counter of times filter was called
    unsigned int                            m_index;
    /// Location(s) where new points are created
    SpatialDomains::PointGeomVector         m_seedPoints;
    /// Location(s) of crossing points are created
    SpatialDomains::PointGeomVector         m_crossPoints;
    /// Stringstream for temporarily holding seed points coordinates
    std::stringstream                       m_seedPointStream;
    /// Frequency for adding new points
    unsigned int                            m_seedFrequency;
    /// Time integration order
    unsigned int                            m_infoSteps;

    /// Flag marking if tracking should be done during or after the simulation
    bool                                    m_postProc;
    /// Number of steps (for post-processing tracking)
    unsigned int                            m_numSteps;
    /// Frequency for updating when tracking during the simulation 
    unsigned int                            m_updateFrequency;

    /// Time-step
    NekDouble                               m_timestep;
    /// Time integration order
    unsigned int                            m_intOrder;

    /// Flag marking if tracking fluid or solid particles
    bool                                    m_fluidParticles;
    /// Flag marking for wear evaluation
    bool                                    m_wear;
    /// Particles specific gravity (for solid particles)
    NekDouble                               m_SG;
    /// Particles Surface roughness angle
    NekDouble                               m_SRA;
    /// Particles Gravity (for solid particles)
    NekDouble                               m_gravityX;
    /// Particles Gravity (for solid particles)
    NekDouble                               m_gravityY;
    /// Particles Gravity (for solid particles)
    NekDouble                               m_gravityZ;
    /// Particles Seeding Velocity  (for solid particles)
    NekDouble                               m_SV;
    /// Particles diameter (for solid particles)
    NekDouble                               m_diameter;
    /// Kinematic viscosity
    NekDouble                               m_kinvis;
    /// Flag marking for dragforce
    bool                                    m_dragforce;
    /// Flag marking for shear lift force
    bool                                    m_shearliftforce;
    /// Flag marking for rotational lift force
    bool                                    m_rotliftforce;
    /// Flag marking for virtual mass
    bool                                    m_virtualmass;
    /// Flag marking for outputvelocity 
    bool                                    m_outputvelocity;
    /// Flag marking for outputforce
    bool                                    m_outputforce;
    /// Flag marking for outputrank
    bool                                    m_outputrank;

    /// Variables for output file
    std::string                             m_outputFile;
    std::ofstream                           m_outputStream;
    unsigned int                            m_outputFrequency;
    std::string                             m_WearFile;
    std::ofstream                           m_WearStream;
    
    /// Variables for output file
    std::string                     m_BoundaryString;
    
    /// ID's of boundary regions where we want the forces
    std::vector<unsigned int>       m_boundaryRegionsIdList;
    /// Determines if a given Boundary Region is in
    /// m_boundaryRegionsIdList
    std::vector<bool>               m_boundaryRegionIsInList;

    /// Variables for bounding box determining domain of interest for particles
    bool                                    m_useBoundingBox;
    Array<OneD, NekDouble>                  m_boundingBox;

    /// Vector storing Particles
    std::vector<Particle>                   m_particles;

    static NekDouble AdamsBashforth_coeffs[4][4];
    static NekDouble AdamsMoulton_coeffs[4][4];

    /// Main function to advance the particles by one time-step
    void AdvanceParticles(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields);
    /// Add seed points to m_particles
    void AddSeedPoints(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields);
	/// Add cross points to m_particles
//    void AddCrossParticles(
//        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields);
    /// Update location of particle (eId and locCoords)
    void UpdateLocCoord(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        Particle &particle);
    /// Interpolate solution to particle location
    void InterpSolution(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        Particle &particle);
    /// Update particle position
    void UpdatePosition(Particle &particle);
    /// Set the particle velocity to match the fluid velocity
    void SetToFluidVelocity(Particle &particle);
    /// Set the particle velocity and force for crossing particles
    void SetToVelForce(const Array<OneD, NekDouble>  VelForce, Particle &particle);
    /// Update particle velocity
    void UpdateVelocity(Particle &particle);
    /// Calculate force (for solid particles)
    void CalculateForce(Particle &particle);
    /// Collision modelling
    void HandleCollision(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        Particle &particle);
    /// Check if particle left the domain of interest
    void CheckBoundingBox(Particle &particle);
    /// Write output information
    void OutputParticles(const NekDouble &time,
		 const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields);
    /// Rotate array for high order time integration
    void RollOver(Array<OneD, Array<OneD, NekDouble> > &input);
};

}
}

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERPARTICLESTRACKING */
