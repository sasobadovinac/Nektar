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

NekDouble FilterParticlesTracking::AdamsBashforth_coeffs[4][4] = {
        { 1.0       , 0.0       , 0.0       , 0.0       },
        { 3.0/2.0   ,-1.0/2.0   , 0.0       , 0.0       },
        { 23.0/12.0 ,-4.0/3.0   , 5.0/12.0  , 0.0       },
        { 55.0/24.0 ,-59.0/24.0 , 37.0/24.0 ,-3.0/8.0   }};
NekDouble FilterParticlesTracking::AdamsMoulton_coeffs[4][4] = {
        { 1.0       ,  0.0      , 0.0       , 0.0       },
        { 1.0/2.0   ,  1.0      , 0.0       , 0.0       },
        { 5.0/12.0  ,  2.0/3.0  ,-1.0/12.0  , 0.0       },
        { 3.0/8.0   ,  19.0/24.0,-5.0/24.0  , 1.0/24.0  }};

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

    m_session->LoadParameter("Kinvis", m_kinvis);

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

    // Determine if tracking is performed during or after the simulation
    it = pParams.find("TimeIntegrationOrder");
    if (it == pParams.end())
    {
        m_intOrder = 1;
    }
    else
    {
        LibUtilities::Equation equ(m_session, it->second);
        m_intOrder = round(equ.Evaluate());
        ASSERTL0(m_intOrder >= 1 && m_intOrder <= 4,
                "TimeIntegrationOrder must be between 1 and 4.");
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

        it = pParams.find("InfoSteps");
        if (it == pParams.end())
        {
            m_infoSteps = 0;
        }
        else
        {
            LibUtilities::Equation equ3(m_session, it->second);
            m_infoSteps = round(equ3.Evaluate());
        }

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
          && m_outputFile.substr(m_outputFile.length() - 4) == ".csv"))
    {
        m_outputFile += ".csv";
    }

    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        if(m_postProc)
        {
            m_outputFrequency = 1;
        }
        else
        {
            m_outputFrequency = m_updateFrequency;
        }
    }
    else
    {
        LibUtilities::Equation equ(m_session, it->second);
        m_outputFrequency = round(equ.Evaluate());
    }

    // Read bounding box, if defined
    it = pParams.find("BoundingBox");
    if (it == pParams.end())
    {
        m_useBoundingBox = false;
    }
    else
    {
        m_useBoundingBox = true;
        m_boundingBox    = Array<OneD, NekDouble> (6);
        std::stringstream bbStream(it->second);
        bbStream >> m_boundingBox[0] >> m_boundingBox[1]
                 >> m_boundingBox[2] >> m_boundingBox[3]
                 >> m_boundingBox[4] >> m_boundingBox[5];

        ASSERTL0( !bbStream.fail(), "BoundingBox needs six entries: "
                                    "xmin, xmax, ymin, ymax, zmin, zmax");
    }

    // Initialise m_index
    m_index = 0;
    m_advanceCalls = 0;
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
    m_outputStream << ", Fx, Fy";
    if(dim == 3)
    {
        m_outputStream << ", Fz";
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
    if (m_updateFrequency && !(m_index % m_updateFrequency))
    {
        AdvanceParticles(pFields);

        if (m_outputFrequency && !(m_index % m_outputFrequency))
        {
            OutputParticles(time + m_timestep);
        }

        if (m_seedFrequency && !(m_index % m_seedFrequency))
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
            if (m_outputFrequency && !( (n+1) % m_outputFrequency))
            {
                OutputParticles(trackTime);
            }

            if (m_seedFrequency && !( (n+1) % m_seedFrequency))
            {
                // Introduce new points in the domain;
                AddSeedPoints(pFields);
            }

            if (m_infoSteps && !( (n+1) % m_infoSteps))
            {
                cout << "Tracking steps: " << setw(8)  << left << n+1 << " "
                     << "Time: "  << setw(12) << left << trackTime << endl;
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
        m_advanceCalls++;
        if(particle.m_used)
        {
             //Store original position of the particle
            for (int i = 0; i < particle.m_dim; ++i)
            {
                particle.m_oldCoord[i] = particle.m_gloCoord[i];
            }
            
            // Rotate force array
            RollOver(particle.m_force);
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
            HandleCollision(pFields,particle);
            // Check if particle left the domain of interest
            CheckBoundingBox(particle);
            
            
            //Plot data for debuging 
            //for(int n = 0; n < 3; ++n)
            //{
                //cout << ", " << boost::format("%25.19e") %
                                        //particle.m_gloCoord[n];
            //}
            //cout <<endl;
            
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
        m_particles.emplace_back(dim, numFields, m_intOrder, gloCoord);
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
    // First check if still in the same element
    bool found = false;
    if(particle.m_eId >= 0)
    {
        found = pFields[0]->GetExp(particle.m_eId)->GetGeom()->ContainsPoint(
                                              particle.m_gloCoord,
                                              particle.m_locCoord,
                                              NekConstants::kNekZeroTol);
        
    }
    // If it changed elements just search again
    if (!found)
    {
        particle.m_eId = pFields[0]->GetExpIndex(particle.m_gloCoord,
                                             particle.m_locCoord,
                                             NekConstants::kNekZeroTol);
    }
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
    int order = min(m_advanceCalls, m_intOrder);
    if( m_fluidParticles)
    {
        // Velocity is at time n, so we use Adams-Bashforth
        for (int i = 0; i < particle.m_dim; ++i)
        {
            for(int j = 0; j < order; ++j)
            {
                particle.m_gloCoord[i] += m_timestep *
                                        AdamsBashforth_coeffs[order-1][j] *
                                        particle.m_particleVelocity[j][i];
            }
        }
    }
    else
    {
        // We already updated velocity to v_{n+1}, so we use Adams-Moulton
        for (int i = 0; i < particle.m_dim; ++i)
        {
            for(int j = 0; j < order; ++j)
            {
                particle.m_gloCoord[i] += m_timestep *
                                          AdamsMoulton_coeffs[order-1][j] *
                                          particle.m_particleVelocity[j][i];
            }
        }
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
        //particle.m_particleVelocity[0][i] = 0.0;
        particle.m_particleVelocity[0][i] = particle.m_fluidVelocity[i];
    }
}

/**
 *
 */
void FilterParticlesTracking::UpdateVelocity(Particle &particle)
{
    if( m_fluidParticles)
    {
        RollOver(particle.m_particleVelocity);
        SetToFluidVelocity(particle);
    }
    else
    {
        CalculateForce(particle);
        RollOver(particle.m_particleVelocity);
        int order = min(m_advanceCalls, m_intOrder);
        for (int i = 0; i < particle.m_dim; ++i)
        {
            for(int j = 0; j < order; ++j)
            {
                particle.m_particleVelocity[0][i] += m_timestep *
                                        AdamsBashforth_coeffs[order-1][j] *
                                        particle.m_force[j][i];
            }
        }
    }
}

/**
 *
 */
void FilterParticlesTracking::CalculateForce(Particle &particle)
{
    // TO DO: update particle.m_force[0][i] with force per unit mass
    
    NekDouble   Re = 0.0, Cd = 0.0, Fd = 0.0;     
    NekDouble   nu = 1.0e-6;     // m2/s
    
    for (int i = 0; i < particle.m_dim; ++i)
    {    
        Re += pow(particle.m_fluidVelocity[i]
                  -particle.m_particleVelocity[0][i],2.0);
    }
    //Re = sqrt(Re)*m_diameter/m_kinvis;
    Re = sqrt(Re)*m_diameter/nu;
    
    
    ////Calcule Drag coeficient
    
    //// Openfoam 
    ////if (Re < 0.1)
    ////{
        ////Cd = 1E4;
    ////}
    ////else if (Re >1000.0)
    ////{
        ////Cd = 0.424;
    ////}
    ////else
    ////{
        ////Cd = (24.0/Re)*(1.0+pow(Re,2.0/3.0)/6.0);
    ////}
    
    //// Crowe et al. (1998) 
    //if (Re == 0.0)
    //{
        //Cd = 1E4; //este valor es muy restrigindo!!!!
    //}
    //else if (Re < 0.5  )
    //{
        //Cd = 24.0 / Re;
    //}
    //else if (Re <1000.0)
    //{
        //Cd = (24.0 / Re) * ( 1.0 + 0.15 * pow(Re,0.687) );
    //}
    //else
    //{
        //Cd = 0.44;
    //}
    
    Cd = 0.47;
    
    ////////Fd = (18.0*m_kinvis/((m_density)
          ////*pow(m_diameter,2.0)))*(Cd*Re/24.0);
    
    Fd = (18.0 * nu / ((m_density*1000)
          * pow(m_diameter,2.0))) * (Cd * Re / 24.0);
    
    std::cout<<"Re "<<Re<<", "<<"Cd "<<Cd<<", "<<"Fd "<<Fd<<std::endl;
    
    for (int i = 0; i < particle.m_dim; ++i)
    {
        particle.m_force[0][i] = Fd * ( particle.m_fluidVelocity[i]
                               - particle.m_particleVelocity[0][i]);
    }
    
    cout<<"Particle Velocity [ "<<particle.m_particleVelocity[0][0]<<","
                             <<particle.m_particleVelocity[0][1]<<","
                             <<particle.m_particleVelocity[0][2]<<"]"<<endl;
    
    //Add gravity effects
    //particle.m_force[0][1] = 0.0;
    particle.m_force[0][1] -= 10 * (1.0 - 1.0/m_density);

    cout<<"Force = "<<particle.m_force[0][1]<<endl;

}

/**
 **********************************************************************
 */
void FilterParticlesTracking::HandleCollision(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    Particle &particle)
{
    cout<<endl<< "Punto 1: {"<<particle.m_oldCoord[0]<<", "
                             <<particle.m_oldCoord[1]<<", "
                             <<particle.m_oldCoord[2]<<"}"<<endl;

    cout<<endl<< "Punto 2: {"<<particle.m_gloCoord[0]<<", "
                             <<particle.m_gloCoord[1]<<", "
                             <<particle.m_gloCoord[2]<<"}"<<endl<<endl; 
    
    // TO DO
    // For now, just mark particles outside the domain as unused
while( particle.m_eId == -1 && particle.m_used == true)
{
    cout << "Particle " << particle.m_id << " collisioned." << endl; 
    //particle.m_used = false;    
 
    //Boundary Expansion
    Array<OneD,const MultiRegions::ExpListSharedPtr> bndExp;
    bndExp = pFields[0]->GetBndCondExpansions();
    LocalRegions::ExpansionSharedPtr elemBndExp;
    
    NekDouble   minDist = 0.0, dist = 0.0, distN = 0.0; 
    NekDouble   maxDotProd = 0.0, dotProd = 0.0, ScaleDP = 0.0;
                
    Array<OneD,double> minNormal(3);
    int minPnt = -1, minBnd = -1;
    
    //Loop over each boundary
    for (int nb = 0; nb < bndExp.num_elements(); ++nb)
    {
        ////Boundary normals - normals[dir][point]
        Array<OneD,Array<OneD, double>> normals;
        pFields[0]->GetBoundaryNormals(nb, normals);
        
        ////Coordinates 
        int npoints    = bndExp[nb]->GetTotPoints();
        Array<OneD,Array<OneD, NekDouble>> coords(3);
        for (int j = 0; j < 3; ++j)
        {
            coords[j] = Array<OneD,NekDouble>(npoints);
        }
        bndExp[nb]->GetCoords(coords[0],coords[1],coords[2]);
        
    //Loop for each integration point on  any boundary
    for (int j = 0; j < npoints; ++j)
    {
        //cout<<"Coordinates {"
        //<< x0[i]<<", "<< x1[i]<<", "<< x2[i]<<"}";
        
        //cout<<" Normals: [";
        //for (int i = 0; i < particle.m_dim; ++j)
        //{    
            //cout << normals[i][j]<<", ";//nq-1
        //}
        //cout<<']'<<endl;
        
        
        ////cout<<"Coordinates {"
                ////<< coords[0][j]<<", "
                ////<< coords[1][j]<<", "
                ////<< coords[2][j]<<"}";
                ////cout <<" Normals: [";

                
        dist = 0.0; distN = 0.0;
        for (int i = 0; i < particle.m_dim; ++i)
        {
            dist +=    (coords[i][j]
                       - particle.m_oldCoord[i])
                       * normals[i][j];
                             
            distN +=   (coords[i][j]
                       - particle.m_gloCoord[i])
                       * normals[i][j];
        }
        ////cout<<"]"<<endl;

        //if (distN*dist < 0 && (abs(dist) < minDist || minDist == 0.0) )
        
        if (dist*distN < 0.0 )
        {   
            dotProd = 0.0; ScaleDP = 0.0;
            for (int i = 0; i < particle.m_dim; ++i)
            {
                dotProd += ( particle.m_gloCoord[i]
                           -   particle.m_oldCoord[i] )
                           * ( coords[i][j]
                           -   particle.m_oldCoord[i] );
                
                
                ScaleDP += pow( particle.m_oldCoord[i]
                               -  coords[i][j] , 2 );
            }
            dotProd /= sqrt(ScaleDP);
            
            //cout<<endl<<"dotProd = "<<dotProd<<endl;
            
            if (dotProd > maxDotProd)
            {
                maxDotProd = dotProd;
                minDist    = abs(dist);
                minPnt     =    j;
                minBnd     =   nb;
                
                for (int i = 0; i < particle.m_dim; ++i)
                {
                    minNormal[i] = normals[i][j];
                }
            }
         
        }
    //cout<<endl;
    
    }
    }
    
    cout<<endl<<"Min Boundary: "<<minBnd<<" Min Point: "<<minPnt<<
    " Distance: "<<minDist<<" Min normal: { "<<minNormal[0]<<","<<
    minNormal[1]<<"}"<<endl<<endl;

    //// Collision point cordinates
    Array<OneD, NekDouble> collPnt(3,0.0);
    NekDouble absVel = 0.0;

    for (int i = 0; i < particle.m_dim; ++i)
    {
        absVel += pow((particle.m_gloCoord[i]
                   -  particle.m_oldCoord[i])*minNormal[i],2);
    }
    absVel = sqrt(absVel);
    
    cout<<"Punto collision: {";
    for (int i = 0; i < particle.m_dim; ++i)
    {
        collPnt[i] = particle.m_oldCoord[i]
                + (  particle.m_gloCoord[i]
                -    particle.m_oldCoord[i])
                *  minDist/absVel;
                
    cout<<collPnt[i]<<", ";
    }
    cout<<"}"<<endl<<endl;
    
    
    //New coordinates and Velocities
    
    cout<<"Vel antes"<<particle.m_particleVelocity[0][1]<<endl;
    cout<<"Force antes"<<particle.m_force[0][1]<<endl;

    // evaluation of the restitution coeficients
    int order = min(m_advanceCalls, m_intOrder);
    //cout<<"New positions: {";
    
     
    NekDouble dotProdCoord, dotProdVel, dotProdForce; 
     
    for(int j = 0; j < order; ++j)
    {   
        dotProdCoord = 0.0;
        dotProdVel = 0.0;
        dotProdForce = 0.0;
    
        for (int i = 0; i < particle.m_dim; ++i)
        {
            
            dotProdCoord += (particle.m_gloCoord[i]
                             - collPnt[i]) * minNormal[i];

            //particle.m_gloCoord[i] =     particle.m_gloCoord[i]
                              //- abs( particle.m_gloCoord[i]
                                   //- collPnt[i] )
                               //* 2 * minNormal[i];
                               
            particle.m_oldCoord[i] =  collPnt[i];
        
        dotProdVel += particle.m_particleVelocity[j][i] * minNormal[i];
        dotProdForce += particle.m_force[j][i] * minNormal[i];
        }
        
        // Velocities
        for (int i = 0; i < particle.m_dim; ++i)
        {
        particle.m_gloCoord[i]  = collPnt[i]
                                + (particle.m_gloCoord[i] - collPnt[i] )
                                - dotProdCoord * 2 * minNormal[i];
        
        particle.m_particleVelocity[j][i] -=dotProdVel * 2 * minNormal[i];
            
        particle.m_force[j][i] -= dotProdForce * 2 * minNormal[i];
        }

    
    //cout<<particle.m_gloCoord[i]<<", ";      
    
    
    
                         
    }
    //cout<<"}"<<endl<<endl;
    
    cout<<"Vel despues"<<particle.m_particleVelocity[0][1]<<endl;
    cout<<"Force despues"<<particle.m_force[0][1]<<endl;
    
    
    UpdateLocCoord(pFields, particle);
}
}

/**
 *
 */
void FilterParticlesTracking::CheckBoundingBox(Particle &particle)
{
    if( m_useBoundingBox)
    {
        for(int i = 0; i < particle.m_dim; ++i)
        {
            if( particle.m_gloCoord[i] < m_boundingBox[2*i] ||
                particle.m_gloCoord[i] > m_boundingBox[2*i+1])
            {
                cout << "Particle " << particle.m_id
                     << " left the bounding box." << endl;
                particle.m_used = false;
            }
        }
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

            for(int n = 0; n < 3; ++n)
            {
                m_outputStream << ", " << boost::format("%25.19e") %
                                            particle.m_gloCoord[n];
            }

            for(int n = 0; n < particle.m_dim; ++n)
            {
                m_outputStream << ", " << boost::format("%25.19e") %
                                            particle.m_particleVelocity[0][n];
            }

            for(int n = 0; n < particle.m_fields.num_elements(); ++n)
            {
                m_outputStream << ", " << boost::format("%25.19e") %
                                            particle.m_fields[n];
            }
            for(int n = 0; n < particle.m_dim; ++n)
            {
                m_outputStream << ", " << boost::format("%25.19e") %
                                            particle.m_force[0][n];
            }
            
            m_outputStream << endl;
        }
    }
}

/**
 *
 */
void FilterParticlesTracking::RollOver(
        Array<OneD, Array<OneD, NekDouble> > &input)
{
    int  nlevels = input.num_elements();
    Array<OneD, NekDouble> tmp = input[nlevels-1];

    for(int n = nlevels-1; n > 0; --n)
    {
        input[n] = input[n-1];
    }

    input[0] = tmp;
}

}
}
