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

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <SolverUtils/Filters/FilterParticlesTracking.h>
#include <boost/format.hpp>
#include <iomanip>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
std::string FilterParticlesTracking::className =
    GetFilterFactory().RegisterCreatorFunction("ParticlesTracking",
                                               FilterParticlesTracking::create);

int Particle::m_nextId = 0;

NekDouble FilterParticlesTracking::AdamsBashforth_coeffs[4][4] = {
    {1.0, 0.0, 0.0, 0.0},
    {3.0 / 2.0, -1.0 / 2.0, 0.0, 0.0},
    {23.0 / 12.0, -4.0 / 3.0, 5.0 / 12.0, 0.0},
    {55.0 / 24.0, -59.0 / 24.0, 37.0 / 24.0, -3.0 / 8.0}};
NekDouble FilterParticlesTracking::AdamsMoulton_coeffs[4][4] = {
    {1.0, 0.0, 0.0, 0.0},
    {1.0 / 2.0, 1.0 / 2.0, 0.0, 0.0},
    {5.0 / 12.0, 2.0 / 3.0, -1.0 / 12.0, 0.0},
    {3.0 / 8.0, 19.0 / 24.0, -5.0 / 24.0, 1.0 / 24.0}};

/**
 *
 */
FilterParticlesTracking::FilterParticlesTracking(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem>      &pEquation,
    const ParamMap &pParams) :
    Filter(pSession, pEquation)
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
        LibUtilities::Equation equ(
        m_session->GetExpressionEvaluator(), it->second);
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
        std::string sOption = it->second.c_str();
        m_postProc          = (boost::iequals(sOption, "true")) ||
                     (boost::iequals(sOption, "yes"));
    }

    // Determine if tracking is performed during or after the simulation
    it = pParams.find("TimeIntegrationOrder");
    if (it == pParams.end())
    {
        m_intOrder = 1;
    }
    else
    {
        LibUtilities::Equation equ1(
            m_session->GetExpressionEvaluator(), it->second);
        m_intOrder = round(equ1.Evaluate());
        ASSERTL0(m_intOrder >= 1 && m_intOrder <= 4,
                 "TimeIntegrationOrder must be between 1 and 4.");
    }

    // Time integration parameters
    if (m_postProc)
    {
        it = pParams.find("NumSteps");
        ASSERTL0(it != pParams.end(),
                 "Missing parameter 'NumSteps' for FilterParticleTracking.");
        LibUtilities::Equation equ2(
            m_session->GetExpressionEvaluator(), it->second);
        m_numSteps = round(equ2.Evaluate());

        it = pParams.find("DeltaT");
        ASSERTL0(it != pParams.end(),
                 "Missing parameter 'DeltaT' for FilterParticleTracking.");
        LibUtilities::Equation equ3(
            m_session->GetExpressionEvaluator(), it->second);
        m_timestep = equ3.Evaluate();

        it = pParams.find("InfoSteps");
        if (it == pParams.end())
        {
            m_infoSteps = 0;
        }
        else
        {
            LibUtilities::Equation equ4(
            m_session->GetExpressionEvaluator(), it->second);
            m_infoSteps = round(equ4.Evaluate());
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
            LibUtilities::Equation equ5(
            m_session->GetExpressionEvaluator(), it->second);
            m_updateFrequency = round(equ5.Evaluate());
        }

        m_session->LoadParameter("TimeStep", m_timestep);
        m_timestep *= m_updateFrequency;
    }

    // Determine if we are tracking fluid or solid particles
    it = pParams.find("Diameter");
    if (it == pParams.end())
    {
        // By default track fluid particles
        m_fluidParticles = true;
        //cout << " - Working with fluid particles" << endl;
    }
    else
    {
        LibUtilities::Equation equ6(
        m_session->GetExpressionEvaluator(), it->second);
        m_diameter       = equ6.Evaluate();
        m_fluidParticles = false;

        // Determine if SpecificGravity has effect
        it = pParams.find("SpecificGravity");
        if (it == pParams.end())
        {
            m_SG = 1.0; // By default SG = 1
            //cout << "- Particles with density of fluid" << endl;
        }
        else
        {
            LibUtilities::Equation equ7(
                m_session->GetExpressionEvaluator(), it->second);
            m_SG = equ7.Evaluate();
        }

        // Determine if gravity has effect
        it = pParams.find("Gravity");
        if (it == pParams.end())
        {
            m_gravity = 0.0; // By default the value of gravity is zero
            //cout << "- Working without Gravity" << endl;
        }
        else
        {
            LibUtilities::Equation equ8(
                m_session->GetExpressionEvaluator(), it->second);
            m_gravity = equ8.Evaluate();
        }
        // Determine if InitialVelocity 
        it = pParams.find("SeedingVelocity");
        if (it == pParams.end())
        {
            m_SV = 1.0; // By default the value of gravity is zero
            //cout << "- Working without Gravity" << endl;
        }
        else
        {
            LibUtilities::Equation equ9(
                m_session->GetExpressionEvaluator(), it->second);
            m_SV = equ9.Evaluate();
        }
        
         // Boundary (to evaluate colision)
         it = pParams.find("BoundaryID");
         ASSERTL0(it != pParams.end(), "Missing parameter 'BoundaryID'");
         ASSERTL0(it->second.length() > 0, "Empty parameter 'BoundaryID'.");
         m_BoundaryString = it->second;
    }
    
    // Set Wear evaluation
    it = pParams.find("Wear");
    if (it == pParams.end())
    {
        // By default track fluid particles
        m_wear = false;
    }
    else
    {
        std::string sOption = it->second.c_str();
        m_wear          = (boost::iequals(sOption, "true")) ||
                          (boost::iequals(sOption, "yes" ));
    }

    // Read variables for output
    it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        m_outputFile    = m_session->GetSessionName();
        m_WearFile = m_session->GetSessionName() + ".wear";
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        m_outputFile    = it->second;
        m_WearFile = it->second;
    }


    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        if (m_postProc)
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
        LibUtilities::Equation equ9(
           m_session->GetExpressionEvaluator(), it->second);
        m_outputFrequency = round(equ9.Evaluate());
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
        m_boundingBox    = Array<OneD, NekDouble>(6);
        std::stringstream bbStream(it->second);
        bbStream >> m_boundingBox[0] >> m_boundingBox[1] >> m_boundingBox[2] >>
            m_boundingBox[3] >> m_boundingBox[4] >> m_boundingBox[5];

        ASSERTL0(!bbStream.fail(), "BoundingBox needs six entries: "
                                   "xmin, xmax, ymin, ymax, zmin, zmax");
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
    // Parse the boundary regions into a list.
    std::string::size_type FirstInd = m_BoundaryString.find_first_of('[') + 1;
    std::string::size_type LastInd  = m_BoundaryString.find_last_of(']') - 1;

    ASSERTL0(FirstInd <= LastInd,
             (std::string("Error reading boundary region definition:") +
              m_BoundaryString)
                 .c_str());
	if (m_wear)
	{
		std::string IndString =
			m_BoundaryString.substr(FirstInd, LastInd - FirstInd + 1);
		bool parseGood =
			ParseUtils::GenerateSeqVector(IndString, m_boundaryRegionsIdList);
		ASSERTL0(parseGood && !m_boundaryRegionsIdList.empty(),
				 (std::string("Unable to read boundary regions index "
							  "range for FilterParticleTracking: ") +
				  IndString)
					 .c_str());
	}
	
    // determine what boundary regions need to be considered
    int cnt;
    unsigned int numBoundaryRegions =
        pFields[0]->GetBndConditions().num_elements();
    m_boundaryRegionIsInList.insert(m_boundaryRegionIsInList.end(),
                                    numBoundaryRegions, 0);

    SpatialDomains::BoundaryConditions bcs(m_session, pFields[0]->GetGraph());
    const SpatialDomains::BoundaryRegionCollection &bregions =
        bcs.GetBoundaryRegions();



   
    cnt = 0;
    for (auto &it : bregions)
    {
        if (std::find(m_boundaryRegionsIdList.begin(),
                      m_boundaryRegionsIdList.end(),
                      it.first) != m_boundaryRegionsIdList.end())
        {
            m_boundaryRegionIsInList[cnt] = 1;
        }
        cnt++;
    }
    // Create map for element and edge/face of each boundary expansion
    //~ pFields[0]->GetBoundaryToElmtMap(m_BCtoElmtID,m_BCtoTraceID);

    // Read seed points
    Array<OneD, NekDouble> newCoord(3, 0.0);
    int dim = pFields[0]->GetGraph()->GetSpaceDimension();
    
    int i = 0;
    while (!m_seedPointStream.fail())
    {
        m_seedPointStream >> newCoord[0] >> newCoord[1] >> newCoord[2];

        if (!m_seedPointStream.fail())
        {
            SpatialDomains::PointGeomSharedPtr vert =
                MemoryManager<SpatialDomains::PointGeom>::AllocateSharedPtr(
                    dim, i, newCoord[0], newCoord[1], newCoord[2]);

            m_seedPoints.push_back(vert);
            ++i;
        }
    }
    	
    //Parallel Comm
    int vRank  = pFields[0]->GetComm()->GetRowComm()->GetRank();

    if (!(m_outputFile.length() >= 4 &&
          m_outputFile.substr(m_outputFile.length() - 4) == ".csv"))
    {
        m_outputFile = m_outputFile +"."+to_string(vRank*0.000001).substr(2)+".csv";
    }
    if (!(m_WearFile.length() >= 4 &&
          m_WearFile.substr(m_WearFile.length() - 4) == ".csv"))
    {
        m_WearFile = m_WearFile+"."+to_string(vRank*0.000001).substr(2)+".csv";
    }

	m_outputStream.open(m_outputFile.c_str());
	if (m_wear)
	{
		m_WearStream.open(m_WearFile.c_str()); 
	}
	
    if (vRank == 0)
    {
		cout << "Tracking Particles..." <<endl;
		
		//Write headers in the csv file
		m_outputStream << "Time, Rank, Id, x, y, z, particleU, particleV";
	   
		if (dim == 3)
		{
			m_outputStream << ", particleW";
		}
		
		for (int n = 0; n < pFields.num_elements(); ++n)
		{
			m_outputStream << ", " << m_session->GetVariables()[n];
		}
		m_outputStream << ", Fx, Fy";
		if (dim == 3)
		{
			m_outputStream << ", Fz";
		}
		m_outputStream << endl;
		
		
		if (m_wear)
	    {
			if (dim == 2)
			{
				m_WearStream<<"x, y, Velocity, angle"<<endl;
			}
			if (dim == 3)
			{
				m_WearStream<<"x, y, z, Velocity, angle"<<endl;
			}
		}
	}


  // Add seed points to m_particles
    AddSeedPoints(pFields);
    OutputParticles(0.0,pFields);
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
            OutputParticles(time + m_timestep,pFields);
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
        for (int n = 0; n < m_numSteps; ++n)
        {
			
            AdvanceParticles(pFields);
            trackTime += m_timestep;

            ++m_index;
            if (m_outputFrequency && !((n + 1) % m_outputFrequency))
            {
                OutputParticles(trackTime,pFields);
            }

            if (m_seedFrequency && !((n + 1) % m_seedFrequency))
            {
                // Introduce new points in the domain;
                AddSeedPoints(pFields);
            }

            if (m_infoSteps && !((n + 1) % m_infoSteps))
            {
                cout << "Tracking steps: " << setw(8) << left << n + 1 << " "
                     << "Time: " << setw(12) << left << trackTime << endl;
            }
        }
    }

    if (pFields[0]->GetComm()->GetRank() == 0)
    {
        m_outputStream.close();
        if(m_wear)
        {
        m_WearStream.close();
        }
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
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
    int vRank     = vComm->GetRowComm()->GetRank();
    int vRankSize = vComm->GetSize();
    int dim       = pFields[0]->GetGraph()->GetSpaceDimension();
    int order     = m_intOrder;//min(particle.m_advanceCalls, m_intOrder);
	    
    Array<OneD, int> particleCounts(vRankSize, 0);
    Array<OneD, int> particleOffsets(vRankSize, 0);
	
    //Loop over particles
    for (auto &particle : m_particles)
    {
        if (particle.m_used)
        {
      	   particle.m_advanceCalls++;
            // Store original position of the particle
            for (int i = 0; i < particle.m_dim; ++i)
            {
                particle.m_oldCoord[i] = particle.m_newCoord[i];
            }
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
            //~ // Check if particle left the domain of interest
            //~ CheckBoundingBox(particle);
            //~ // Handle collisions if particle left the domain
            HandleCollision(pFields, particle);
            if (particle.m_eId ==-1)
	         {
      		    particleCounts[vRank]++; 
	         }
         }	
    }
    
   // Manage the particle leaving 
   vComm->AllReduce(particleCounts, LibUtilities::ReduceSum);
   int nTotparticles = Vmath::Vsum(vRankSize,particleCounts,1);
   if (nTotparticles>0)
   {
      for (int i = 1; i < vRankSize; ++i)
      {
         particleOffsets[i] = particleOffsets[i-1] + particleCounts[i-1];
      }
    
      Array<OneD, Array<OneD, NekDouble>> crossParticles(3);
      Array<OneD, int> crossRankSend(nTotparticles, -1);
      Array<OneD, int> crossParticlesId(nTotparticles, -1);
      for (int i = 0; i < 3; ++i)
      {
          crossParticles[i] = Array<OneD, NekDouble>(nTotparticles,0.0);
      }

      //Get cordinates from crossing particles and rank
      int j = 0;
      for (auto &particle : m_particles)
      {	
          if (particle.m_eId==-1 && particle.m_used)
          {
         crossRankSend[j+particleOffsets[vRank]] = vRank;
         crossParticlesId[j+particleOffsets[vRank]] = particle.m_id;
         for (int i = 0; i < 3; ++i)
         {
            crossParticles[i][j+particleOffsets[vRank]] =  particle.m_newCoord[i];
         }
         j++;
          }
      }
       
      vComm->AllReduce(crossRankSend, LibUtilities::ReduceMax);
      vComm->AllReduce(crossParticlesId, LibUtilities::ReduceMax);
      vComm->AllReduce(crossParticles[0], LibUtilities::ReduceSum);
      vComm->AllReduce(crossParticles[1], LibUtilities::ReduceSum);
      vComm->AllReduce(crossParticles[2], LibUtilities::ReduceSum);

      //Read the particle crossing array
      Array<OneD, NekDouble> newCoord(3, 0.0);
      Array<OneD, NekDouble> locCoord(3, 0.0);
      Array<OneD, int> crossRankRecv(nTotparticles, -1);
      Array<OneD, int> listId(nTotparticles, -1);

      for(int j=0;j<nTotparticles;j++)
      {
          if (crossRankSend[j]!=vRank)
          {
             newCoord[0] = crossParticles[0][j];
             newCoord[1] = crossParticles[1][j];
             newCoord[2] = crossParticles[2][j];
             listId[j] = pFields[0]->GetExpIndex(newCoord,locCoord,NekConstants::kNekZeroTol);
                         
             if(listId[j]>=0)
             {
                crossRankRecv[j] = vRank;
             }
          }
      }
          
      vComm->AllReduce(crossRankRecv, LibUtilities::ReduceMax);
      vComm->AllReduce(listId, LibUtilities::ReduceMax);

      //Comunication
      for(int j=0; j < nTotparticles; ++j)
      {
          if (listId[j]>=0)
          {
           if (crossRankSend[j]>=0 && crossRankRecv[j]>=0)
           {   
             Array<OneD, NekDouble> VelForce(2*dim*order, 0.0); 
             if(vRank == crossRankSend[j])
             {
                for (auto &particle : m_particles)
                {	
                   if(particle.m_id == crossParticlesId[j])
                   {   
                       int order = min(particle.m_advanceCalls, m_intOrder);
                       for (int k = 0; k < order; ++k)
                       {
                           for (int i = 0; i < dim; ++i)
                           {
                               VelForce[i+k*dim] = particle.m_particleVelocity[k][i];
                               VelForce[(dim*order)+i+k*dim] = particle.m_force[k][i];
                           }
                       }
                    particle.m_used = false;
                    vComm->GetRowComm()->Send(crossRankRecv[j], VelForce);
                    break;
                   }	
                }
             }
             if(vRank == crossRankRecv[j])
             {
               vComm->GetRowComm()->Recv(crossRankSend[j], VelForce);
            
               newCoord[0] = crossParticles[0][j];
               newCoord[1] = crossParticles[1][j];
               newCoord[2] = crossParticles[2][j];
             
               bool endInsert = true;
               // Found a space inside the particle list 
               //to  include points in unused positions of m_particles
            
              for (auto &particle : m_particles)
              {
                 if (particle.m_used == false)
                 {
                  // Change coordinates of particle
                  particle.SetCoord(newCoord);
                  // Obtain new id for particle
                  particle.SetNewId();
                  // Change m_used flag
                  particle.m_used = true;
                  // Find location of new particle
                  UpdateLocCoord(pFields, particle);
                  // Initialise particle velocity to match fluid velocity
                  InterpSolution(pFields, particle);
                  // Initialise particle velocity to match fluid velocity
                  SetToVelForce(VelForce, particle);
                  endInsert = false; 
                  break;
                }	
              }
              // Add particle to the end of m_particles
              if (endInsert)
              {
                m_particles.reserve(1);
                m_particles.emplace_back(dim,pFields.num_elements(),m_intOrder, newCoord);
                // Find location of new particle	
                UpdateLocCoord(pFields, m_particles.back());
                // Initialise particle velocity to match fluid velocity
                InterpSolution(pFields, m_particles.back());
                // Initialise particle velocity to match fluid velocity
                SetToVelForce(VelForce, m_particles.back());
               }
             }	
            }
          }
          else
          {
               // Just the rank with the particle evaluate collision
               if (crossRankSend[j]==vRank)
               {
                   for (auto &particle : m_particles)
                   {	
                        if(particle.m_id == crossParticlesId[j])
                        {   
                           // Handle collisions if particle left the domain
                           // HandleCollision(pFields, particle);

                            // Particle is unused because I wanna know why
                            cout<<"Particle unused  because I wanna know why ID: "
                                <<particle.m_id<<" Cross: "<<particle.m_eId<<endl;
                            particle.m_used         = false;
                            particle.m_advanceCalls = 0;

                        }	
                   }
               }
          }
       }	
   }					
}

/**
 *
 */
void FilterParticlesTracking::AddSeedPoints(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    int insertedPoints 	  = 0;
    int totalPoints    	  = m_seedPoints.size();
    int dim            	  = pFields[0]->GetGraph()->GetSpaceDimension();
    int numFields      	  = pFields.num_elements();
    Array<OneD, NekDouble> newCoord(3, 0.0);
    Array<OneD, NekDouble> locCoord(3, 0.0);
    
    
    // Avoid problems with empty list of points...
    if (totalPoints == 0)
    {	
		ASSERTL0(false,"Empty list of points");
        return;
    }
    
	// Parallel Comm
	 LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
	 int vRank     = vComm->GetRowComm()->GetRank();
	 int vRankSize = vComm->GetSize();
	 Array<OneD, int> listId(vRankSize, -2);	

	 // First include points in unused positions of m_particles
    for (auto &particle : m_particles)
    {
        if (particle.m_used == false)
        {
            m_seedPoints[insertedPoints++]-> GetCoords(newCoord[0],
 	                 					       newCoord[1], newCoord[2]);
                                                      
		      listId[vRank]  = pFields[0]->GetExpIndex(newCoord,
		       		     		locCoord,NekConstants::kNekZeroTol);
							
            if(listId[vRank]>=0)
            {
               vComm->AllReduce(listId, LibUtilities::ReduceMax);
               // Change coordinates of particle
               particle.SetCoord(newCoord);
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
            else
            {
               vComm->AllReduce(listId, LibUtilities::ReduceMax);
               ASSERTL0((Vmath::Vsum(vRankSize, listId, 1)!=(-vRankSize)),
                                          "Point is not in the domain ");
            }
	     }	
        if (insertedPoints == totalPoints)
        {
            return;
        }	 
	  }
	
    // Add remaining particles to the end of m_particles
    int newSize = m_particles.size() + totalPoints - insertedPoints;
    m_particles.reserve(newSize);
    for (int i = insertedPoints; i < totalPoints; ++i)
    {
         m_seedPoints[i]->GetCoords(newCoord[0], newCoord[1], newCoord[2]);	
         listId[vRank]  = pFields[0]->GetExpIndex(newCoord,
                      locCoord,NekConstants::kNekZeroTol);
         
         if(listId[vRank]>=0)
         {
            vComm->AllReduce(listId, LibUtilities::ReduceMax);
            m_particles.emplace_back(dim, numFields, m_intOrder, newCoord);
            // Find location of new particle	
            UpdateLocCoord(pFields, m_particles.back());
            // Initialise particle velocity to match fluid velocity
            InterpSolution(pFields, m_particles.back());
            SetToFluidVelocity(m_particles.back());
         }
         else
         {
            vComm->AllReduce(listId, LibUtilities::ReduceMax);
            ASSERTL0((Vmath::Vsum(vRankSize, listId, 1)!=(-vRankSize)),
                                         "Point is not in the domain ");
          }	
     }
}


/**
 *
 */
void FilterParticlesTracking::AddCrossParticles(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    //bool inserted = false;
    //Array<OneD, NekDouble> newCoord(3, 0.0);
    //Array<OneD, NekDouble> locCoord(3, 0.0);
    
    //// Parallel Comm
     //LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
     //int vRank     = vComm->GetRowComm()->GetRank();
     //int vRankSize = vComm->GetSize();
     //Array<OneD, int> listId(vRankSize, -2);	

     //~ // Found a space inside the particle list 
     //~ //to  include points in unused positions of m_particles
    //~ for (auto &particle : m_particles)
    //~ {
        //~ if (particle.m_used == false)
        //~ {
	    //~ // Change coordinates of particle
	    //~ particle.SetCoord(newCoord);
	    //~ // Obtain new id for particle
	    //~ particle.SetNewId();
	    //~ // Change m_used flag
	    //~ particle.m_used = true;
	    //~ // Find location of new particle
	    //~ UpdateLocCoord(pFields, particle);
	    //~ // Initialise particle velocity to match fluid velocity
	    //~ InterpSolution(pFields, particle);
	    //~ particle.m_particleVelocity[j][i]
	    
	    
	    //~ UpdateVelocity(particle);
	    // Integrate velocity to obtain position
	    //~ UpdatePosition(particle);
	    //~ // Update element containing particle and coordinate in std element
	    //~ UpdateLocCoord(pFields, particle);
	    //~ inserted = true; 
	    //~ break;
	//~ }	
    //~ }
	
    //~ // Add particle to the end of m_particles
    //~ if (inserted)
    //~ {
	//~ m_particles.reserve(1);
	//~ UpdateLocCoord(pFields, m_particles.back());
	//~ // Initialise particle velocity to match fluid velocity
	//~ InterpSolution(pFields, m_particles.back());
	//~ UpdateVelocity(particle);
	//~ // Integrate velocity to obtain position
	//~ UpdatePosition(particle);
	//~ // Update element containing particle and coordinate in std element
	//~ UpdateLocCoord(pFields, particle);
    //~ }
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
    if (particle.m_eId >= 0)
    {
        found = pFields[0]
                    ->GetExp(particle.m_eId)
                    ->GetGeom()
                    ->ContainsPoint(particle.m_newCoord, particle.m_locCoord,
                                    NekConstants::kNekZeroTol);
    }
    // If it changed elements just search again
    if (!found)
    {
        particle.m_eId =
            pFields[0]->GetExpIndex(particle.m_newCoord, particle.m_locCoord,
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
	if (particle.m_used == true)
     {	
    	Array<OneD, NekDouble> physvals;
		for (int i = 0; i < pFields.num_elements(); ++i)
		{
			physvals = pFields[i]->UpdatePhys() +
					   pFields[i]->GetPhys_Offset(particle.m_eId);

			particle.m_fields[i] =
				pFields[i]
					->GetExp(particle.m_eId)
					->StdPhysEvaluate(particle.m_locCoord, physvals);
		}

		// TO DO: This could be changed later to allow solvers using different
		//        variables (e.g. compressible flow solver)
		for (int i = 0; i < particle.m_dim; ++i)
		{
			particle.m_fluidVelocity[i] = particle.m_fields[i];
		}
	}	
}

/**
 *
 */
void FilterParticlesTracking::UpdatePosition(Particle &particle)
{
    int order = min(particle.m_advanceCalls, m_intOrder);
    if (m_fluidParticles)
    {
        // Velocity is at time n, so we use Adams-Bashforth
        for (int i = 0; i < particle.m_dim; ++i)
        {
            for (int j = 0; j < order; ++j)
            {
                particle.m_newCoord[i] += m_timestep *
                                          AdamsBashforth_coeffs[order - 1][j] *
                                          particle.m_particleVelocity[j][i];
            }
        }
    }
    else
    {
        // We already updated velocity to v_{n+1}, so we use Adams-Moulton
        for (int i = 0; i < particle.m_dim; ++i)
        {
            for (int j = 0; j < order; ++j)
            {
                particle.m_newCoord[i] += m_timestep *
                                          AdamsMoulton_coeffs[order - 1][j] *
                                          particle.m_particleVelocity[j][i];
            }
        }
    }
}

/**
 *
 */
void FilterParticlesTracking::SetToFluidVelocity(Particle &particle)
{
    for (int i = 0; i < particle.m_dim; ++i)
    {
        particle.m_particleVelocity[0][i] = m_SV *particle.m_fluidVelocity[i];
    }
}

void FilterParticlesTracking::SetToVelForce(Array<OneD, NekDouble> VelForce, Particle &particle)
{
    int dim = particle.m_dim;
    int order = m_intOrder;
    for (int k = 0; k < order; ++k)
    {
        for (int i = 0; i < dim; ++i)
        {
            particle.m_particleVelocity[k][i] = VelForce[i+k*dim];
            particle.m_force[k][i] = VelForce[(dim*order)+i+k*dim];
        }
    }
   //particle.m_advanceCalls = order;
}


/**
 *
 */
void FilterParticlesTracking::UpdateVelocity(Particle &particle)
{
    if (m_fluidParticles)
    {
        RollOver(particle.m_particleVelocity);
        SetToFluidVelocity(particle);
    }
    else
    {
        // Rotate force array
        RollOver(particle.m_force);
        CalculateForce(particle);

        RollOver(particle.m_particleVelocity);

        int order = min(particle.m_advanceCalls, m_intOrder);
        for (int i = 0; i < particle.m_dim; ++i)
        {
            if (order != 1)
            {
                particle.m_particleVelocity[0][i] =
                    particle.m_particleVelocity[1][i];
            }

            for (int j = 0; j < order; ++j)
            {
                particle.m_particleVelocity[0][i] +=
                    m_timestep * AdamsBashforth_coeffs[order - 1][j] *
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
    // Update particle.m_force[0][i] with force per unit mass

    // Particular Re evaluation
    NekDouble Re = 0.0, Cd = 0.0, Fd = 0.0;
    for (int i = 0; i < particle.m_dim; ++i)
    {
        Re +=
            pow(particle.m_fluidVelocity[i] - particle.m_particleVelocity[0][i],
                2.0);
    }
    Re = sqrt(Re) * m_diameter / m_kinvis;

    // Calcule Drag coeficient: Crowe et al. (1998), Lain et al. (2009)
    if (Re < 0.1)
    {
        Re = 0.1;
        Cd = 240.0;
    }
    else if (Re < 0.5)
    {
        Cd = 24.0 / Re;
    }
    else if (Re < 1000.0)
    {
        Cd = (24.0 / Re) * (1.0 + 0.15 * pow(Re, 0.687));
    }
    else
    {
        Cd = 0.44;
    }
    // Evaluate the drag force
    Fd = (3.0 * m_kinvis * Cd * Re) / (4.0 * m_SG * pow(m_diameter, 2.0));

    for (int i = 0; i < particle.m_dim; ++i)
    {
        particle.m_force[0][i] = Fd * (particle.m_fluidVelocity[i] -
                                       particle.m_particleVelocity[0][i]);
    }

    // Add gravity and buoyancy effects on y direction
    /* particle.m_force[0][1] += m_gravity * (1.0 - 1.0 / m_SG); */
    
    // Add gravity and buoyancy effects on x direction
    particle.m_force[0][0] += m_gravity * (1.0 - 1.0 / m_SG);
}

/**
 **********************************************************************
 */
void FilterParticlesTracking::HandleCollision(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    Particle &particle)
{
    NekDouble distN = 0.0;
//    while (particle.m_eId == -1 && particle.m_used == true)
    if (particle.m_eId == -1 && particle.m_used == true)
    {
        // Boundary Expansion:  Boundary list
        Array<OneD, const MultiRegions::ExpListSharedPtr> bndExp;
        bndExp = pFields[0]->GetBndCondExpansions();

        NekDouble minDist = 0.0, dist = 0.0;
        NekDouble maxDotProd = 0.0, dotProd = 0.0, ScaleDP = 0.0;
        Array<OneD, double> minNormal(3);
        int minBnd = -1;
        
        // Loop over each boundary finding cross point
        for (int nb = 0; nb < bndExp.num_elements(); ++nb)
        {
    	      // Boundary normals on each quadrature points in each nb 
    	      // normals[dir][point]
            Array<OneD, Array<OneD, double>> normals;
            pFields[0]->GetBoundaryNormals(nb, normals);

            // Get Coordinates of quadrature points
            int npoints = bndExp[nb]->GetTotPoints();
            Array<OneD, Array<OneD, NekDouble>> coords(3);
            for (int j = 0; j < 3; ++j)
            {
                coords[j] = Array<OneD, NekDouble>(npoints);
            }
            bndExp[nb]->GetCoords(coords[0], coords[1], coords[2]);
				
            // Loop for each integration point on  any boundary
            // Evaluating if the point cross the boundary at each
            // quadrature point
            for (int j = 0; j < npoints; ++j)
            {
              dist  = 0.0; distN = 0.0;
              for (int i = 0; i < particle.m_dim; ++i)
              {
                  dist += (coords[i][j] - particle.m_oldCoord[i])
                          * normals[i][j];

                  distN += (coords[i][j] - particle.m_newCoord[i])
                           * normals[i][j];
              }
              // Check if the wall is crossed
              if (dist * distN < 0.0 ||
                       abs(dist) < NekConstants::kNekZeroTol ||
                         abs(distN) < NekConstants::kNekZeroTol)
               {
                  dotProd = 0.0; ScaleDP = 0.0;
                  // Evaluate the dot Product
                  for (int i = 0; i < particle.m_dim; ++i)
                  {
                      dotProd +=
                          (particle.m_newCoord[i] - particle.m_oldCoord[i]) *
                          (coords[i][j] - particle.m_oldCoord[i]);

                      ScaleDP +=
                       pow(coords[i][j] - particle.m_oldCoord[i], 2);
                  }
                  dotProd /= sqrt(ScaleDP);
					
                  // Save the max dot Product data
                  if (abs(dotProd) > maxDotProd)
                  {
                      maxDotProd = dotProd;
                      minDist    = abs(dist);
                      minBnd     = nb;

                      for (int i = 0; i < particle.m_dim; ++i)
                      {
                          minNormal[i] = normals[i][j];
                      }
                  }              
                }
              }
        }
	
        // Check is the particle collision or leave the domain
        if (minBnd != -1 )
         {
              if (m_boundaryRegionIsInList[minBnd] == 1 )
              {
                  //Particle collisioned
                  // Magnitude and directions of incident vector
                  NekDouble VecMag = 0.0;
                  Array<OneD, NekDouble> VecDir(3, 0.0);          
                  for (int i = 0; i < particle.m_dim; ++i)
                  {
                      VecMag +=  pow(particle.m_newCoord[i]
                                 - particle.m_oldCoord[i],2);
                  }
                  VecMag = sqrt(VecMag);
                  for (int i = 0; i < particle.m_dim; ++i)
                  {
                      VecDir[i] =  (particle.m_newCoord[i]
                                    - particle.m_oldCoord[i])/VecMag;
                  } 
                          
                  // Collision point cordinates
                  Array<OneD, NekDouble> collPnt(3, 0.0);
                  for (int i = 0; i < particle.m_dim; ++i)
                  {
                    collPnt[i] = particle.m_oldCoord[i] + VecDir[i] * minDist; 
                  }
            
            
                  //~ int CPeId = pFields[0]->GetExpIndex(collPnt, 
                  //~ particle.m_locCoord, NekConstants::kNekZeroTol);
                  //~ cout<<"CP: ";
                  //~ for (int j = 0; j < 3; ++j)
                  //~ {
                  //~ 	cout<< collPnt[j]<< " ";   
                  //~ }
                  //~ cout<<  CPeId <<endl<<endl;

                  //~ xr2=[-2*u(1)*u(2)*x(2) - 2*u(1)*u(3)*x(3) + (-2*u(1)*u(1)+1)*x(1)
                  //~ -2*u(2)*u(1)*x(1) - 2*u(2)*u(3)*x(3) + (-2*u(2)*u(2)+1)*x(2)
                  //~ -2*u(3)*u(1)*x(1) - 2*u(3)*u(2)*x(2) + (-2*u(3)*u(3)+1)*x(3)]

                  // New coordinates and Velocities
                  // Evaluation of the restitution coeficients
                  int order              = min(particle.m_advanceCalls, m_intOrder);
                  NekDouble dotProdCoord = 0.0, dotProdVel = 0.0, dotProdForce = 0.0;
                  NekDouble angle = 0.0, Vel = 0.0;
       
                  for (int i = 0; i < particle.m_dim; ++i)
                  {
                      angle +=  abs(VecDir[i] * minNormal[i]);
                      Vel += pow(particle.m_particleVelocity[0][i],2);
                  }
                  angle = asinf(angle); Vel = sqrt(Vel); 
                  
                  // Evaluate dot products to make the reflection
                  for (int i = 0; i < particle.m_dim; ++i)
                  {
                      dotProdCoord +=
                          (particle.m_newCoord[i] - collPnt[i]) * minNormal[i];
                  }
                  
                  // Update coordinates
                  for (int i = 0; i < particle.m_dim; ++i)
                  {
                      //particle.m_oldCoord[i] = collPnt[i];

                      particle.m_newCoord[i] = collPnt[i] +
                                               (particle.m_newCoord[i] - collPnt[i]) -
                                               dotProdCoord * 2 * minNormal[i];
                  }
                  UpdateLocCoord(pFields, particle);
             
                  for (int j = 0; j < order; ++j)
                  {
                      // Evaluate dot products to make the reflection
                      for (int i = 0; i < particle.m_dim; ++i)
                      {
                          dotProdVel +=
                              particle.m_particleVelocity[j][i] * minNormal[i];
                              
                          dotProdForce += particle.m_force[j][i] * minNormal[i];
                      }

                      // Update velocites and forces
                      for (int i = 0; i < particle.m_dim; ++i)
                      {
                          particle.m_particleVelocity[j][i] -=
                                                      dotProdVel * 2 * minNormal[i];

                          particle.m_force[j][i] -= dotProdForce * 2 * minNormal[i];
                      }
                  }

                  // Evaluate the erosion
                  if (m_wear)
                  {
                      // Output the collision information
                      for (int n = 0; n < particle.m_dim; ++n)
                      {
                        m_WearStream << "   "
                           << boost::format("%25.19e") % collPnt[n];
                      }
                        m_WearStream <<"   "<< boost::format("%25.19e") % Vel
                            <<"   "<< boost::format("%25.19e") % angle  <<endl;
                  }
                  
                  // Evaluate the change in the collision position
                  distN = 1.0;
                  for (int i = 0; i < particle.m_dim; ++i)
                  {
                      //distN += pow(collPnt[i] - collPntOLD[i], 2);
                      distN += pow(collPnt[i] - particle.m_newCoord[i], 2);
                  }

                  if (distN < m_diameter / 2.0)
                  {
                      // Particle  is stalled
                     cout<<"Particle stalled after collision ID: "<<particle.m_id<<" distance: "<<distN <<endl;
                     particle.m_used         = false;
                     particle.m_advanceCalls = 0;
                  }
              }
         }
     }
     else
     {
         // Evaluate the change in the position when the parcticle no collisioned
         if (particle.m_advanceCalls == -1 )
            {
               distN = 0.0;
               for (int i = 0; i < particle.m_dim; ++i)
               {
                   //distN += pow(collPnt[i] - collPntOLD[i], 2);
                   distN += pow(particle.m_oldCoord[i] - particle.m_newCoord[i], 2);
                   //cout<<" ID: "<<particle.m_id<<" distance: "<<distN <<endl;
               }

               if (distN <m_diameter)
               {
                   // Particle  is stalled
                   cout<<" Particle stalled without collision ID: "<<particle.m_id<<" distance: "<<distN <<endl;
                   cout<<" iter "<<particle.m_advanceCalls<<endl;
                particle.m_used         = false;
                particle.m_advanceCalls = 0;
            }
         }
     } 
}

/**
 *
 */
void FilterParticlesTracking::CheckBoundingBox(Particle &particle)
{
    if (m_useBoundingBox)
    {
        for (int i = 0; i < particle.m_dim; ++i)
        {
            if (particle.m_newCoord[i] < m_boundingBox[2 * i] ||
                particle.m_newCoord[i] > m_boundingBox[2 * i + 1])
            {
               // "Particle left the bounding box." 
                particle.m_used         = false;
                particle.m_advanceCalls = 0;
            }
        }
    }
}

/**
 *
 */
void FilterParticlesTracking::OutputParticles(const NekDouble &time,
const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
	for (auto &particle : m_particles)
    {
        if (particle.m_used)
        {
            m_outputStream << boost::format("%25.19e") % time;
            m_outputStream << ", "<< pFields[0]->GetComm()
									->GetRowComm()->GetRank();
            m_outputStream << ", " << particle.m_id;

            for (int n = 0; n < 3; ++n)
            {
                m_outputStream
                    << ", "
                    << boost::format("%25.19e") % particle.m_newCoord[n];
            }

            for (int n = 0; n < particle.m_dim; ++n)
            {
                m_outputStream << ", "
                               << boost::format("%25.19e") %
                                      particle.m_particleVelocity[0][n];
            }

            for (int n = 0; n < particle.m_fields.num_elements(); ++n)
            {
                m_outputStream
                    << ", " << boost::format("%25.19e") % particle.m_fields[n];
            }
            for (int n = 0; n < particle.m_dim; ++n)
            {
                m_outputStream
                    << ", "
                    << boost::format("%25.19e") % particle.m_force[0][n];
            }
            m_outputStream << endl;
        }
    }
}

/**
   * Function to roll time-level storages to the next step layout.
   * The stored data associated with the oldest time-level
   * (not required anymore) are moved to the top, where they will
   * be overwritten as the solution process progresses.
   *
**/
void FilterParticlesTracking::RollOver(
    Array<OneD, Array<OneD, NekDouble>> &input)
{
    int nlevels = input.num_elements();
    Array<OneD, NekDouble> tmp = input[nlevels - 1];

    for (int n = nlevels - 1; n > 0; --n)
    {
        input[n] = input[n - 1];
    }

    input[0] = tmp;
}
}
}
