///////////////////////////////////////////////////////////////////////////////
//
// File FilterParticlesTracking.cpp
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

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <iomanip>
#include <SolverUtils/Filters/FilterParticlesTracking.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <boost/format.hpp>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
std::string FilterParticlesTracking::className =
    GetFilterFactory().RegisterCreatorFunction(
        "ParticlesTracking",FilterParticlesTracking::create);

int Particle::m_nextId = 0;

/**
 *
 */
FilterParticlesTracking::FilterParticlesTracking(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const ParamMap &pParams):
    Filter(pSession)
{
    // Read parameters
    ParamMap::const_iterator it;

    // Seed points
    it = pParams.find("SeedPoints");
    ASSERTL0(it != pParams.end(),
                "Missing parameter 'SeedPoints' for FilterParticleTracking.");
    m_seedPointStream.str(it->second);

    // Frequency for adding seed points
    it = pParams.find("SeedFrequency");
    if (it == pParams.end())
    {
        // By default only add points once at the beginning
        m_seedFrequency = 0;
    }
    else
    {
        LibUtilities::Equation equ(m_session, it->second);
        m_seedFrequency = round(equ.Evaluate());
    }

    // Determine if tracking is performed during or after the simulation
    it = pParams.find("PostProcessing");
    if (it == pParams.end())
    {
        // By default perform tracking as a post-processing step
        //      TO DO: we could change this to false later
        m_postProc = true;
    }
    else
    {
        std::string sOption =
                        it->second.c_str();
        m_postProc = ( boost::iequals(sOption,"true")) ||
                     ( boost::iequals(sOption,"yes"));
    }

    // Time integration parameters
    if( m_postProc)
    {
        it = pParams.find("NumSteps");
        ASSERTL0(it != pParams.end(),
                "Missing parameter 'NumSteps' for FilterParticleTracking.");
        LibUtilities::Equation equ(m_session, it->second);
        m_numSteps = round(equ.Evaluate());

        it = pParams.find("DeltaT");
        ASSERTL0(it != pParams.end(),
                "Missing parameter 'DeltaT' for FilterParticleTracking.");
        LibUtilities::Equation equ2(m_session, it->second);
        m_timestep = equ2.Evaluate();

        // Use m_updateFrequency = 0 to avoid running during the simulation
        m_updateFrequency = 0;
    }
    else
    {
        it = pParams.find("UpdateFrequency");
        if (it == pParams.end())
        {
            m_updateFrequency = 1;
        }
        else
        {
            LibUtilities::Equation equ(m_session, it->second);
            m_updateFrequency = round(equ.Evaluate());
        }

        m_session->LoadParameter("TimeStep", m_timestep);
        m_timestep *= m_updateFrequency;
    }

    // Determine if we are tracking fluid or solid particles
    it = pParams.find("FluidParticles");
    if (it == pParams.end())
    {
        // By default track fluid particles
        //      TO DO: we could change this to false later
        m_fluidParticles = true;
    }
    else
    {
        std::string sOption =
                        it->second.c_str();
        m_fluidParticles = ( boost::iequals(sOption,"true")) ||
                     ( boost::iequals(sOption,"yes"));
    }

    // Read parameters for solid particles
    if( !m_fluidParticles)
    {
        it = pParams.find("Density");
        ASSERTL0(it != pParams.end(),
                "Missing parameter 'Density' for FilterParticleTracking.");
        LibUtilities::Equation equ(m_session, it->second);
        m_density = equ.Evaluate();

        it = pParams.find("Diameter");
        ASSERTL0(it != pParams.end(),
                "Missing parameter 'Diameter' for FilterParticleTracking.");
        LibUtilities::Equation equ2(m_session, it->second);
        m_diameter = equ2.Evaluate();
    }

    // Read variables for output
    it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        m_outputFile = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        m_outputFile = it->second;
    }
    if (!(m_outputFile.length() >= 4
          && m_outputFile.substr(m_outputFile.length() - 4) == ".his"))
    {
        m_outputFile += ".his";
    }

    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = m_updateFrequency;
    }
    else
    {
        LibUtilities::Equation equ(m_session, it->second);
        m_outputFrequency = round(equ.Evaluate());
    }

    // Initialise m_index
    m_index = 0;
}


/**
 *
 */
FilterParticlesTracking::~FilterParticlesTracking()
{
}


/**
 *
 */
void FilterParticlesTracking::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{	
    int dim = pFields[0]->GetGraph()->GetSpaceDimension();
    // Open output stream
    m_outputStream.open(m_outputFile.c_str());
    m_outputStream << "Time, Id, x, y, z, particleU, particleV";
    if(dim == 3)
    {
        m_outputStream << ", particleW";
    }
    for(int n = 0; n < pFields.num_elements(); ++n)
    {
        m_outputStream << ", " << m_session->GetVariables()[n];
    }
    m_outputStream << endl;

    // Read seed points
    Array<OneD, NekDouble>  gloCoord(3,0.0);

    int i = 0;
    while (!m_seedPointStream.fail())
    {
        m_seedPointStream >> gloCoord[0]
                          >> gloCoord[1]
                          >> gloCoord[2];

        if (!m_seedPointStream.fail())
        {
            SpatialDomains::PointGeomSharedPtr vert
                = MemoryManager<SpatialDomains::PointGeom>
                    ::AllocateSharedPtr(dim, i, gloCoord[0],
                                    gloCoord[1], gloCoord[2]);

            m_seedPoints.push_back(vert);
            ++i;
        }
    }

    // Add seed points to m_particles
    AddSeedPoints(pFields);

    // Write initial particle results
    if( m_postProc)
    {
        OutputParticles(0.0);
    }
    else
    {
        OutputParticles(time);
    }

    // Advance particles
    v_Update(pFields, time);
}


/**
 *
 */
void FilterParticlesTracking::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    ++m_index;
    if (m_updateFrequency && (m_index) % m_updateFrequency)
    {
        AdvanceParticles(pFields);

        if (m_outputFrequency && (m_index) % m_outputFrequency)
        {
            OutputParticles(time);
        }

        if (m_seedFrequency && (m_index) % m_seedFrequency)
        {
            // Introduce new points in the domain;
            AddSeedPoints(pFields);
        }
    }	
}

/**
 *
 */
void FilterParticlesTracking::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    if (m_postProc)
    {
        NekDouble trackTime = 0.0;
        for(int n = 0; n < m_numSteps; ++n)
        {
            AdvanceParticles(pFields);
            trackTime += m_timestep;

            ++m_index;
            if (m_outputFrequency && (m_index) % m_outputFrequency)
            {
                OutputParticles(trackTime);
            }

            if (m_seedFrequency && (m_index) % m_seedFrequency)
            {
                // Introduce new points in the domain;
                AddSeedPoints(pFields);
            }
        }
    }

    if (pFields[0]->GetComm()->GetRank() == 0)
    {
        m_outputStream.close();
    }
}


/**
 *
 */
bool FilterParticlesTracking::v_IsTimeDependent()
{
    return true;
}

/**
 *
 */
void FilterParticlesTracking::AdvanceParticles(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    for (auto &particle : m_particles)
    {
        if(particle.m_used)
        {
            // Obtain solution (an fluid velocity) at the particle location
            InterpSolution(pFields, particle);
            // Set particle velocity
            //    - equal to fluid velocity for fluid particles
            //    - using dv/dt = F/m for solid particles
            UpdateVelocity(particle);
            // Integrate velocity to obtain position
            UpdatePosition(particle);
            // Update element containing particle and coordinate in std element
            UpdateLocCoord(pFields, particle);
            // Handle collisions if particle left the domain
            HandleCollision(particle);
        }
    }
}

/**
 *
 */
void FilterParticlesTracking::AddSeedPoints(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    int insertedPoints = 0;
    int totalPoints    = m_seedPoints.size();
    int dim            = pFields[0]->GetGraph()->GetSpaceDimension();
    int numFields      = pFields.num_elements();
    Array<OneD, NekDouble>  gloCoord(3,0.0);

    // Avoid problems with empty list of points...
    if( totalPoints == 0)
    {
        return;
    }

    // First include points in unused positions of m_particles
    for (auto &particle : m_particles)
    {
        if( particle.m_used == false)
        {
            m_seedPoints[insertedPoints++]->GetCoords(gloCoord[0],
                                                    gloCoord[1],
                                                    gloCoord[2]);
            // Change coordinates of particle
            particle.SetCoord(gloCoord);
            // Obtain new id for particle
            particle.SetNewId();
            // Change m_used flag
            particle.m_used = true;
            // Find location of new particle
            UpdateLocCoord(pFields, particle);
            // Initialise particle velocity to match fluid velocity
            InterpSolution(pFields, particle);
            SetToFluidVelocity(particle);
        }
        if ( insertedPoints == totalPoints)
        {
            return;
        }
    }

    // Add remaining particles to the end of m_particles
    int newSize = m_particles.size() + totalPoints - insertedPoints;
    m_particles.reserve(newSize);
    for( int i = insertedPoints; i < totalPoints; ++i)
    {
        m_seedPoints[i]->GetCoords(gloCoord[0],
                                   gloCoord[1],
                                   gloCoord[2]);
        m_particles.emplace_back(dim, numFields, gloCoord);
        // Find location of new particle
        UpdateLocCoord(pFields, m_particles.back());
        // Initialise particle velocity to match fluid velocity
        InterpSolution(pFields, m_particles.back());
        SetToFluidVelocity(m_particles.back());
    }
}

/**
 *
 */
void FilterParticlesTracking::UpdateLocCoord(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    Particle &particle)
{
    // TO DO: This would probably be more efficient if we first check
    //        if the particle is still inside element m_eId
    particle.m_eId = pFields[0]->GetExpIndex(particle.m_gloCoord,
                                             particle.m_locCoord,
                                             NekConstants::kNekZeroTol);
}

/**
 *
 */
void FilterParticlesTracking::InterpSolution(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    Particle &particle)
{
    Array<OneD, NekDouble> physvals;
    for (int i = 0; i < pFields.num_elements(); ++i)
    {
        physvals = pFields[i]->UpdatePhys()
                    + pFields[i]->GetPhys_Offset(particle.m_eId);

        particle.m_fields[i] =
                pFields[i]->GetExp(particle.m_eId)->StdPhysEvaluate(
                        particle.m_locCoord,physvals);
    }

    // TO DO: This could be changed later to allow solvers using different
    //        variables (e.g. compressible flow solver)
    for (int i = 0; i < particle.m_dim; ++i)
    {
        particle.m_fluidVelocity[i] = particle.m_fields[i];
    }
}

/**
 *
 */
void FilterParticlesTracking::UpdatePosition(Particle &particle)
{
    // TO DO: Generalise time integration to higher order
    for (int i = 0; i < particle.m_dim; ++i)
    {
        particle.m_gloCoord[i] = particle.m_gloCoord[i] +
                                 m_timestep * particle.m_particleVelocity[i];
    }
}

/**
 *
 */
void FilterParticlesTracking::SetToFluidVelocity(
    Particle &particle)
{
    for (int i = 0; i < particle.m_dim; ++i)
    {
        particle.m_particleVelocity[i] = particle.m_fluidVelocity[i];
    }
}

/**
 *
 */
void FilterParticlesTracking::UpdateVelocity(Particle &particle)
{
    if( m_fluidParticles)
    {
        SetToFluidVelocity(particle);
    }
    else
    {
        CalculateForce(particle);
        for (int i = 0; i < particle.m_dim; ++i)
        {
            // TO DO: v_{n+1} = v{n} + dT/m * F
        }
    }
}

/**
 *
 */
void FilterParticlesTracking::CalculateForce(Particle &particle)
{
    // TO DO
}

/**
 *
 */
void FilterParticlesTracking::HandleCollision(Particle &particle)
{
    // TO DO
    // For now, just mark particles outside the domain as unused
    if( particle.m_eId == -1)
    {
        cout << "Particle " << particle.m_id << " left the domain." << endl;
        particle.m_used = false;
    }
}

/**
 *
 */
void FilterParticlesTracking::OutputParticles(const NekDouble &time)
{
    for (auto &particle : m_particles)
    {
        if( particle.m_used)
        {
            m_outputStream << boost::format("%25.19e") % time;
            m_outputStream << ", " << particle.m_id;

            for(int n = 0; n < particle.m_dim; ++n)
            {
                m_outputStream << ", " << boost::format("%25.19e") %
                                            particle.m_particleVelocity[n];
            }

            for(int n = 0; n < particle.m_fields.num_elements(); ++n)
            {
                m_outputStream << ", " << boost::format("%25.19e") %
                                            particle.m_fields[n];
            }
        }
    }
}

}
}
