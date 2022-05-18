///////////////////////////////////////////////////////////////////////////////
//
// File FilterAeroForces.cpp
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
// Description: Output values of aerodynamic forces during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <iomanip>
#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string.hpp>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <MultiRegions/ExpList.h>     
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <SolverUtils/Filters/FilterAeroForces.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
std::string FilterAeroForces::className =
        GetFilterFactory().RegisterCreatorFunction(
                "AeroForces", FilterAeroForces::create);

/**
 *
 */
FilterAeroForces::FilterAeroForces(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem>      &pEquation,
    const ParamMap &pParams) :
    Filter(pSession, pEquation)
{
    // OutputFile
    auto it = pParams.find("OutputFile");
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
          && m_outputFile.substr(m_outputFile.length() - 4) == ".fce"))
    {
        m_outputFile += ".fce";
    }

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(
            m_session->GetInterpreter(), it->second);
        m_outputFrequency = round(equ.Evaluate());
    }

    // Time after which we need to calculate the forces
    it = pParams.find("StartTime");
    if (it == pParams.end())
    {
        m_startTime = 0;
    }
    else
    {
        LibUtilities::Equation equ(
            m_session->GetInterpreter(), it->second);
        m_startTime = equ.Evaluate();
    }

    // For 3DH1D, OutputAllPlanes or just average forces?
    m_session->MatchSolverInfo("Homogeneous", "1D",
                               m_isHomogeneous1D, false);
    if(m_isHomogeneous1D)
    {
        it = pParams.find("OutputAllPlanes");
        if (it == pParams.end())
        {
            m_outputAllPlanes = false;
        }
        else
        {
            std::string sOption =
                            it->second.c_str();
            m_outputAllPlanes = ( boost::iequals(sOption,"true")) ||
                                ( boost::iequals(sOption,"yes"));
        }
    }
    else
    {
        m_outputAllPlanes = false;
    }

    // Boundary (to calculate forces on)
    it = pParams.find("Boundary");
    ASSERTL0(it != pParams.end(),   "Missing parameter 'Boundary");
    ASSERTL0(it->second.length() > 0, "Empty parameter 'Boundary'.");
    m_BoundaryString = it->second;

    //
    // Directions (to project forces)
    //

    // Allocate m_directions and m_directions0
    // m_directions0 is to save the original direction defined by the user
    // and projected every time to the new orientation which is then be saved
    // in the m_direction, note if moving reference frame is not being used,
    // m_directions and m_directions0 will be the same
    m_directions = Array<OneD, Array<OneD, NekDouble> > (3);
    m_directions0 = Array<OneD, Array<OneD, NekDouble> > (3);
    //Initialise directions to default values (ex, ey, ez)
    for (int i = 0; i < 3; ++i)
    {
        m_directions[i] = Array<OneD, NekDouble>(3, 0.0);
        m_directions[i][i] = 1.0;

        m_directions0[i] = Array<OneD, NekDouble>(3, 0.0);
        m_directions0[i][i] = 1.0;

    }
    std::stringstream       directionStream;
    std::string             directionString;
    //Override with input from xml file (if defined)
    for (int i = 0; i < 3; ++i)
    {
        std::stringstream tmp;
        tmp << i+1;
        std::string dir = "Direction" + tmp.str();
        it = pParams.find(dir);
        if ( it != pParams.end() )
        {
            ASSERTL0(!(it->second.empty()),
                     "Missing parameter '"+dir+"'.");
            directionStream.str(it->second);
            // Guarantee the stream is in its start position
            //      before extracting
            directionStream.clear();
            // normalisation factor
            NekDouble norm = 0.0;
            for (int j = 0; j < 3; j++)
            {
                directionStream >> directionString;
                if (!directionString.empty())
                {
                    LibUtilities::Equation equ(
                        m_session->GetInterpreter(), directionString);
                    m_directions[i][j] = equ.Evaluate();
                    norm += m_directions[i][j]*m_directions[i][j];
                }
            }
            // Normalise direction
            for( int j = 0; j < 3; j++)
            {
                m_directions[i][j] /= sqrt(norm);
            }
        }
    }

    for(int i = 0; i < 3; ++i)
    {
        for(int j = 0; j < 3; ++j)
        {
            m_directions0[i][j] = m_directions[i][j];
        }
    }

    // Point around which we compute the moments (default to the origin)
    m_pivotPoint = Array<OneD, NekDouble> (3, 0.0);
    it = pParams.find("PivotPoint");
    if( it == pParams.end() )
    {
        it = pParams.find("MomentPoint");
    }

    if ( it != pParams.end() )
    {
        ASSERTL0(!(it->second.empty()), "Missing parameter 'PivotPoint'.");
        std::stringstream       pivotPointStream;
        std::string             pivotPointString;
        pivotPointStream.str(it->second);

        for (int j = 0; j < 3; ++j)
        {
            pivotPointStream >> pivotPointString;
            if (!pivotPointString.empty())
            {
                LibUtilities::Equation equ(m_session->GetInterpreter(),
                                           pivotPointString);
                m_pivotPoint[j] = equ.Evaluate();
            }
        }
    }

}


/**
 *
 */
FilterAeroForces::~FilterAeroForces()
{

}

/**
 *
 */
void FilterAeroForces::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Load mapping
    m_mapping = GlobalMapping::Mapping::Load(m_session, pFields);

    // Parse the boundary regions into a list.
    std::string::size_type FirstInd =
                            m_BoundaryString.find_first_of('[') + 1;
    std::string::size_type LastInd =
                            m_BoundaryString.find_last_of(']') - 1;

    ASSERTL0(FirstInd <= LastInd,
            (std::string("Error reading boundary region definition:") +
             m_BoundaryString).c_str());

    std::string IndString =
            m_BoundaryString.substr(FirstInd, LastInd - FirstInd + 1);
    bool parseGood = ParseUtils::GenerateSeqVector(IndString,
                                                   m_boundaryRegionsIdList);
    ASSERTL0(parseGood && !m_boundaryRegionsIdList.empty(),
             (std::string("Unable to read boundary regions index "
              "range for FilterAeroForces: ") + IndString).c_str());

    // determine what boundary regions need to be considered
    int cnt;
    unsigned int numBoundaryRegions =
                        pFields[0]->GetBndConditions().size();
    m_boundaryRegionIsInList.insert(m_boundaryRegionIsInList.end(),
                                    numBoundaryRegions, 0);

    SpatialDomains::BoundaryConditions bcs(m_session,
                                            pFields[0]->GetGraph());
    const SpatialDomains::BoundaryRegionCollection &bregions =
                                            bcs.GetBoundaryRegions();

    cnt = 0;
    for (auto &it : bregions)
    {
        if ( std::find(m_boundaryRegionsIdList.begin(),
                       m_boundaryRegionsIdList.end(), it.first) !=
                m_boundaryRegionsIdList.end() )
        {
            m_boundaryRegionIsInList[cnt] = 1;
        }
        cnt++;
    }

    // Create map for element and edge/face of each boundary expansion
    if(m_isHomogeneous1D)
    {
        pFields[0]->GetPlane(0)->GetBoundaryToElmtMap
                                        (m_BCtoElmtID,m_BCtoTraceID);
    }
    else
    {
        pFields[0]->GetBoundaryToElmtMap(m_BCtoElmtID,m_BCtoTraceID);
    }

    // Define number of planes  to calculate the forces
    //     in the Homogeneous direction ( if m_outputAllPlanes is false,
    //      consider only first plane in wave space)
    // If flow has no Homogeneous direction, use 1 to make code general
    if(m_isHomogeneous1D &&(m_outputAllPlanes || m_mapping->IsDefined()))
    {
        m_nPlanes = pFields[0]->GetHomogeneousBasis()->
                                            GetZ().size();
    }
    else
    {
        m_nPlanes = 1;
    }

    // Create map for Planes ID for Homogeneous direction
    //    If flow has no Homogeneous direction, create trivial map
    int j;
    m_planesID = Array<OneD, int> (m_nPlanes,-1);
    if(m_isHomogeneous1D)
    {
        Array<OneD, const unsigned int> IDs = pFields[0]->GetZIDs();
        //Loop through all planes
        for(cnt = 0; cnt < m_nPlanes; cnt++)
        {
            for(j = 0; j < IDs.size(); ++j)
            {
                //Look for current plane ID in this process
                if(IDs[j] == cnt)
                {
                    break;
                }
            }
            // Assign ID to planesID
            // If plane is not found in this process, leave it with -1
            if(j != IDs.size())
            {
                m_planesID[cnt] = j;
            }
        }
    }
    else
    {
        m_planesID[0] = 0;
    }

    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

    // Write header
    int expdim = pFields[0]->GetGraph()->GetMeshDimension();
    int momdim = (expdim == 2) ? 1 : 3;
    if (vComm->GetRank() == 0)
    {
        // Open output stream
        bool adaptive;
        m_session->MatchSolverInfo("Driver", "Adaptive",
                                    adaptive, false);
        if (adaptive)
        {
            m_outputStream.open(m_outputFile.c_str(), ofstream::app);
        }
        else
        {
            m_outputStream.open(m_outputFile.c_str());
        }
        m_outputStream << "# Forces and moments acting on bodies" << endl;
        for( int i = 0; i < expdim; i++ )
        {
            m_outputStream << "#" << " Direction" << i+1 << " = (";
            m_outputStream.width(8);
            m_outputStream << setprecision(4) << m_directions[i][0];
            m_outputStream.width(8);
            m_outputStream << setprecision(4) << m_directions[i][1];
            m_outputStream.width(8);
            m_outputStream << setprecision(4) << m_directions[i][2];
            m_outputStream << ")" << endl;
        }
        
        m_outputStream << "#" << " Moments around" << " (";
        m_outputStream.width(8);
        m_outputStream << setprecision(4) << m_pivotPoint[0];
        m_outputStream.width(8);
        m_outputStream << setprecision(4) << m_pivotPoint[1];
        m_outputStream.width(8);
        m_outputStream << setprecision(4) << m_pivotPoint[2];
        m_outputStream << ")" << endl;
        
        m_outputStream << "# Boundary regions: " << IndString.c_str() << endl;
        m_outputStream << "#";
        m_outputStream.width(7);
        m_outputStream << "Time";
        for( int i = 1; i <= expdim; i++ )
        {
            m_outputStream.width(8);
            m_outputStream <<  "F" << i << "-press";
            m_outputStream.width(9);
            m_outputStream <<  "F" << i << "-visc";
            m_outputStream.width(8);
            m_outputStream <<  "F" << i << "-total";
        }
        if (momdim == 1)
        {
            // 2D case: we only have M3 (z-moment)
            m_outputStream.width(8);
            m_outputStream <<  "M" << 3 << "-press";
            m_outputStream.width(9);
            m_outputStream <<  "M" << 3 << "-visc";
            m_outputStream.width(8);
            m_outputStream <<  "M" << 3 << "-total";
        }
        else
        {
            for( int i = 1; i <= momdim; i++ )
            {
                m_outputStream.width(8);
                m_outputStream <<  "M" << i << "-press";
                m_outputStream.width(9);
                m_outputStream <<  "M" << i << "-visc";
                m_outputStream.width(8);
                m_outputStream <<  "M" << i << "-total";
            }
        }
        if( m_outputAllPlanes )
        {
            m_outputStream.width(10);
            m_outputStream << "Plane";
        }
        if (m_session->DefinesSolverInfo("HomoStrip"))
        {
            ASSERTL0(m_outputAllPlanes==false,
                    "Output forces on all planes not compatible with HomoStrips");
            m_outputStream.width(10);
            m_outputStream << "Strip";
        }

        m_outputStream << endl;
    }

    m_lastTime = -1;
    m_index = 0;
    v_Update(pFields, time);

}

/**
 *
 */
void FilterAeroForces::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Only output every m_outputFrequency.
    if ( m_outputFrequency == 0)
    {
        return;
    }
    if ((m_index++) % m_outputFrequency  || (time < m_startTime))
    {
        return;
    }

    // Calculate the forces
    CalculateForces(pFields, time);


    // Calculate forces including all planes
    int expdim = pFields[0]->GetGraph()->GetMeshDimension();
    int momdim = (expdim == 2) ? 1 : 3;
    Array<OneD, NekDouble>  Fp(expdim,0.0);
    Array<OneD, NekDouble>  Fv(expdim,0.0);
    Array<OneD, NekDouble>  Ft(expdim,0.0);
    for( int i = 0; i < expdim; i++)
    {
        Fp[i] = Vmath::Vsum(m_nPlanes, m_Fpplane[i], 1) / m_nPlanes;
        Fv[i] = Vmath::Vsum(m_nPlanes, m_Fvplane[i], 1) / m_nPlanes;
        Ft[i] = Fp[i] + Fv[i];
    }
    
    Array<OneD, NekDouble>  Mp(momdim,0.0);
    Array<OneD, NekDouble>  Mv(momdim,0.0);
    Array<OneD, NekDouble>  Mt(momdim,0.0);
    for( int i = 0; i < momdim; i++)
    {
        Mp[i] = Vmath::Vsum(m_nPlanes, m_Mpplane[i], 1) / m_nPlanes;
        Mv[i] = Vmath::Vsum(m_nPlanes, m_Mvplane[i], 1) / m_nPlanes;
        Mt[i] = Mp[i] + Mv[i];
    }

    // Communicators to exchange results
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

    // get the total number of planes per strip for the case of homoStrip method
    int nZ{1};
    int nstrips{1};
    int colSize;
    int idOffset{0};

    if( m_session->DefinesSolverInfo("HomoStrip"))
    {
        m_session->LoadParameter("Strip_Z", nstrips);
        colSize = vComm->GetColumnComm()->GetSize();
        idOffset = colSize/nstrips;
    }


    if(m_isHomogeneous1D && m_session->DefinesSolverInfo("HomoStrip"))
    {
        m_session->LoadParameter("HomModesZ",nZ);
    }

    //Write Results
    if (vComm->GetRank() == 0)
    {
        // Write result in each plane
        if( m_outputAllPlanes)
        {
            for( int plane = 0; plane < m_nPlanes; plane++)
            {
                // Write time
                m_outputStream.width(8);
                m_outputStream << setprecision(6) << time;
                // Write forces
                for( int i = 0; i < expdim; i++ )
                {
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << m_Fpplane[i][plane];
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << m_Fvplane[i][plane];
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << m_Ftplane[i][plane];
                }
                // Write moments
                for( int i = 0; i < momdim; i++ )
                {
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << m_Mpplane[i][plane];
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << m_Mvplane[i][plane];
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << m_Mtplane[i][plane];
                }
                m_outputStream.width(10);
                m_outputStream << plane;
                m_outputStream << endl;
            }
        }
        // Output average (or total) force
        m_outputStream.width(8);
        m_outputStream << setprecision(6) << time;
        for( int i = 0; i < expdim; i++)
        {
            m_outputStream.width(15);
            m_outputStream << setprecision(8) << Fp[i];
            m_outputStream.width(15);
            m_outputStream << setprecision(8) << Fv[i];
            m_outputStream.width(15);
            m_outputStream << setprecision(8) << Ft[i];
        }
        for( int i = 0; i < momdim; i++)
        {
            m_outputStream.width(15);
            m_outputStream << setprecision(8) << Mp[i];
            m_outputStream.width(15);
            m_outputStream << setprecision(8) << Mv[i];
            m_outputStream.width(15);
            m_outputStream << setprecision(8) << Mt[i];
        }
        if( m_outputAllPlanes)
        {
            m_outputStream.width(10);
            m_outputStream << "average";
        }

        if( m_session->DefinesSolverInfo("HomoStrip"))
        {
            // The result we already wrote is for strip 0
            m_outputStream.width(10);
            m_outputStream << 0;
            m_outputStream << endl;

            // Now get result from other strips

            for(int i = 1; i<nstrips; i++)
            {
                int id = i * idOffset;
                vComm->GetColumnComm()->Recv(id, Fp);
                vComm->GetColumnComm()->Recv(id, Fv);
                vComm->GetColumnComm()->Recv(id, Ft);
                vComm->GetColumnComm()->Recv(id, Mp);
                vComm->GetColumnComm()->Recv(id, Mv);
                vComm->GetColumnComm()->Recv(id, Mt);

                m_outputStream.width(8);
                m_outputStream << setprecision(6) << time;
                for( int j = 0; j < expdim; j++)
                {
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << Fp[j];
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << Fv[j];
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << Ft[j];
                }
                for( int j = 0; j < momdim; j++)
                {
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << Mp[j];
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << Mv[j];
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << Mt[j];
                }
            m_outputStream.width(10);
            m_outputStream << i;
            m_outputStream << endl;
            }
        }
        else
        {
            m_outputStream << endl;
        }
    }
    else
    {
        // In the homostrips case, we need to send result to root
        if (m_session->DefinesSolverInfo("HomoStrip") && 
                (vComm->GetRowComm()->GetRank() == 0) )
        {
            // we start from strip 1 not 0, since we already have the data of strip 0
            for(int i=1; i<nstrips; ++i)
            {
                int rid = i * idOffset;
                int sid = vComm->GetColumnComm()->GetRank();
                if( rid == sid )
                {
                    vComm->GetColumnComm()->Send(0, Fp);
                    vComm->GetColumnComm()->Send(0, Fv);
                    vComm->GetColumnComm()->Send(0, Ft);
                    vComm->GetColumnComm()->Send(0, Mp);
                    vComm->GetColumnComm()->Send(0, Mv);
                    vComm->GetColumnComm()->Send(0, Mt);
                }
            }

        }
    }

}


/**
 *
 */
void FilterAeroForces::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(time);

    if (pFields[0]->GetComm()->GetRank() == 0)
    {
        m_outputStream.close();
    }
}


/**
 *
 */
bool FilterAeroForces::v_IsTimeDependent()
{
    return true;
}

/**
 *     This function outputs the force on all planes of the current
 *          process, in the format required by ForcingMovingBody
 */
void FilterAeroForces::GetForces(
                    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
                    Array<OneD, NekDouble> &Aeroforces,
                    const NekDouble &time)
{
    // Calculate forces if the result we have is not up-to-date
    if(time > m_lastTime)
    {
        CalculateForces(pFields, time);
    }
    // Get information to write result
    Array<OneD, unsigned int> ZIDs = pFields[0]->GetZIDs();
    int local_planes = ZIDs.size();
    int expdim = pFields[0]->GetGraph()->GetMeshDimension();

    // Copy results to Aeroforces
    if (m_outputAllPlanes)
    {
        for(int plane = 0 ; plane < local_planes; plane++)
        {
            for(int dir =0; dir < expdim; dir++)
            {
                Aeroforces[plane + dir*local_planes] =
                        m_Ftplane[dir][ZIDs[plane]];
            }
        }
    }
    else
    {
        for(int plane = 0 ; plane < local_planes; plane++)
        {
            for(int dir =0; dir < expdim; dir++)
            {
                Aeroforces[plane + dir*local_planes] = m_Ft[dir];
                       // m_Ftplane[dir][0];
            }
        }
    }
}

/**
 *     This function outputs the moments of force on all planes of the current
 *          process
 */
void FilterAeroForces::GetMoments(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    Array<OneD, NekDouble> &moments, const NekDouble &time)
{
    // Calculate forces if the result we have is not up-to-date
    if(time > m_lastTime)
    {
        CalculateForces(pFields, time);
    }

    // Get information to write result
    Array<OneD, unsigned int> ZIDs = pFields[0]->GetZIDs();
    int local_planes = ZIDs.size();
    int expdim = pFields[0]->GetGraph()->GetMeshDimension();
    int momdim = (expdim == 2) ? 1 : 3;

    // Copy results to moments
    if (m_outputAllPlanes)
    {
        for(int plane = 0 ; plane < local_planes; plane++)
        {
            for(int dir = 0; dir < momdim; ++dir)
            {
                moments[plane + dir*local_planes] =
                        m_Mtplane[dir][ZIDs[plane]];
            }
        }
    }
    else
    {
        for(int plane = 0 ; plane < local_planes; ++plane)
        {
            for(int dir = 0; dir < momdim; dir++)
            {
                moments[plane + dir*local_planes] = m_Mt[dir];
            }
        }
    }
}

/**
 *     This function calculates the forces
 */
void FilterAeroForces::CalculateForces(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    // Store time so we can check if result is up-to-date
    m_lastTime = time;

    // If a mapping is defined, call specific function
    //   Note: CalculateForcesMapping should work without a mapping,
    //         but since it is not very efficient the way it is now,
    //         it is only used when actually required
    if (m_mapping->IsDefined())
    {
        CalculateForcesMapping( pFields, time);
        return;
    }

    // Lock equation system weak pointer
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer expired");

    auto fluidEqu = std::dynamic_pointer_cast<FluidInterface>(equ);
    ASSERTL0(fluidEqu, "Aero forces filter is incompatible with this solver.");

    // update the direction vectors
    // only effective if we use moving reference frame
    {
        bnu::matrix<NekDouble> projMat = bnu::identity_matrix<NekDouble>(3,3);
        // get the projection matrix to transform between moving frame and
        // inertial stationary frame
        // note: it will be unity in case we are not using the moving frame
        fluidEqu->GetMovingFrameProjectionMat(projMat);
        // update the direction vectors
        // loop over the directions (ex, ey, ez)
        for(int idir=0; idir < 3; ++idir)
        {
            bnu::vector<NekDouble> v0 = bnu::zero_vector<NekDouble>(3);
            bnu::vector<NekDouble> v1 = bnu::zero_vector<NekDouble>(3);
            // copy the directions
            for (int j = 0; j < 3; ++j)
            {
                v0(j) = m_directions0[idir][j];
            }

            v1 = bnu::prec_prod(projMat, v0);
            // update the direction matrix
            for (int j = 0; j < 3; ++j)
            {
                m_directions[idir][j] = v1(j);
            }
        }
    }

    int cnt , elmtid, nq, offset, boundary;
    // Get number of quadrature points and dimensions
    int physTot = pFields[0]->GetNpoints();
    int expdim = pFields[0]->GetGraph()->GetMeshDimension();
    int momdim = (expdim == 2) ? 1 : 3;
    int nVel = expdim;
    if( m_isHomogeneous1D )
    {
        nVel = nVel + 1;
    }

    LocalRegions::ExpansionSharedPtr elmt;

    // Fields used to calculate forces (a single plane for 3DH1D)
    Array<OneD, MultiRegions::ExpListSharedPtr>
                                    fields( pFields.size() );

    // Arrays of variables in field
    Array<OneD, Array<OneD, NekDouble> > physfields(pFields.size());
    Array<OneD, Array<OneD, NekDouble> > velocity(nVel);
    Array<OneD, NekDouble>               pressure;

    // Arrays of variables in the element
    Array<OneD, Array<OneD, NekDouble> >       velElmt(expdim);
    Array<OneD, NekDouble>                     pElmt(physTot);

    // Velocity gradient
    Array<OneD, Array<OneD, NekDouble> >       grad( expdim*expdim);
    Array<OneD, NekDouble>                     div;

    Array<OneD, Array<OneD, NekDouble> >  coords(3);

    // Values at the boundary
    Array<OneD, NekDouble>                     Pb;
    Array<OneD, Array<OneD, NekDouble> >       gradb( expdim*expdim);
    Array<OneD, Array<OneD, NekDouble> >  coordsb(3);

    // Communicators to exchange results
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
    LibUtilities::CommSharedPtr rowComm = vComm->GetRowComm();
    LibUtilities::CommSharedPtr colComm =
                            m_session->DefinesSolverInfo("HomoStrip") ?
                                vComm->GetColumnComm()->GetColumnComm():
                                vComm->GetColumnComm();

    // Arrays with forces in each plane
    m_Fp = Array<OneD, NekDouble> (expdim,0.0);
    m_Fv = Array<OneD, NekDouble> (expdim,0.0);
    m_Ft = Array<OneD, NekDouble> (expdim,0.0);
    m_Mp = Array<OneD, NekDouble> (momdim,0.0);
    m_Mv = Array<OneD, NekDouble> (momdim,0.0);
    m_Mt = Array<OneD, NekDouble> (momdim,0.0);

    m_Fpplane = Array<OneD, Array<OneD, NekDouble> >  (expdim);
    m_Fvplane = Array<OneD, Array<OneD, NekDouble> >  (expdim);
    m_Ftplane = Array<OneD, Array<OneD, NekDouble> >  (expdim);
    for(int i = 0; i < expdim; i++)
    {
        m_Fpplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
        m_Fvplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
        m_Ftplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
    }
    
    // Arrays with moments in each plane
    m_Mpplane = Array<OneD, Array<OneD, NekDouble> >  (momdim);
    m_Mvplane = Array<OneD, Array<OneD, NekDouble> >  (momdim);
    m_Mtplane = Array<OneD, Array<OneD, NekDouble> >  (momdim);
    for( int i = 0; i < momdim; ++i)
    {
        m_Mpplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
        m_Mvplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
        m_Mtplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
    }

    // Forces per element length in a boundary
    Array<OneD, Array<OneD, NekDouble> >       fp( expdim );
    Array<OneD, Array<OneD, NekDouble> >       fv( expdim );

    // Moments per element length in a boundary
    Array<OneD, Array<OneD, NekDouble> >       mp( momdim );
    Array<OneD, Array<OneD, NekDouble> >       mv( momdim );

    // Get viscosity
    NekDouble mu{1}, rho{1.};
    if(m_session->DefinesParameter("Kinvis"))
    {
        rho = (m_session->DefinesParameter("rho"))
            ? (m_session->GetParameter("rho")) : 1;
        mu = rho * m_session->GetParameter("Kinvis");
    }
    else
    {
        mu = m_session->GetParameter("mu");
    }
    NekDouble lambda = -2.0/3.0;

    // Perform BwdTrans: when we only want the mean force in a 3DH1D
    //     we work in wavespace, otherwise we use physical space
    for(int i = 0; i < pFields.size(); ++i)
    {
        if (m_isHomogeneous1D && m_outputAllPlanes)
        {
            pFields[i]->SetWaveSpace(false);
        }
        pFields[i]->BwdTrans(pFields[i]->GetCoeffs(),
                             pFields[i]->UpdatePhys());
        pFields[i]->SetPhysState(true);
    }

    // Define boundary expansions
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
    Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
    if(m_isHomogeneous1D)
    {
        BndConds = pFields[0]->GetPlane(0)->GetBndConditions();
        BndExp   = pFields[0]->GetPlane(0)->GetBndCondExpansions();
    }
    else
    {
        BndConds = pFields[0]->GetBndConditions();
        BndExp   = pFields[0]->GetBndCondExpansions();
    }

    // For Homogeneous, calculate force on each 2D plane
    // Otherwise, m_nPlanes = 1, and loop only runs once
    for(int plane = 0; plane < m_nPlanes; plane++ )
    {
        // Check if plane is in this proc
        if( m_planesID[plane] != -1 )
        {
            // For Homogeneous, consider the 2D expansion
            //      on the current plane
            if(m_isHomogeneous1D)
            {
                for(int n = 0; n < pFields.size(); n++)
                {
                    fields[n] = pFields[n]->GetPlane(m_planesID[plane]);
                }
            }
            else
            {
                for(int n = 0; n < pFields.size(); n++)
                {
                    fields[n] = pFields[n];
                }
            }

            // Get velocity and pressure values
            for(int n = 0; n < physfields.size(); ++n)
            {
                physfields[n] = fields[n]->GetPhys();
            }
            for(int n = 0; n < nVel; ++n)
            {
                velocity[n] = Array<OneD, NekDouble>(fields[n]->GetTotPoints());
            }
            pressure = Array<OneD, NekDouble>(fields[0]->GetTotPoints());
            fluidEqu->GetVelocity(physfields, velocity);
            fluidEqu->GetPressure(physfields, pressure);

            //Loop all the Boundary Regions
            for(int n = cnt = 0; n < BndConds.size(); n++)
            {
                if(m_boundaryRegionIsInList[n] == 1)
                {
                    for (int i=0; i < BndExp[n]->GetExpSize(); ++i,cnt++)
                    {
                        elmtid = m_BCtoElmtID[cnt];
                        elmt   = fields[0]->GetExp(elmtid);
                        nq     = elmt->GetTotPoints();
                        offset = fields[0]->GetPhys_Offset(elmtid);

                        // Extract  fields on this element
                        for(int j=0; j<expdim; j++)
                        {
                            velElmt[j] = velocity[j] + offset;
                        }
                        pElmt = pressure + offset;

                        // Compute the velocity gradients
                        div = Array<OneD, NekDouble>(nq,0.0);
                        for (int j=0; j<expdim; j++)
                        {
                            for (int k=0; k<expdim; k++)
                            {
                                grad[j*expdim+k] =
                                        Array<OneD, NekDouble>(nq,0.0);
                                elmt->PhysDeriv(k,velElmt[j],
                                        grad[j*expdim+k]);

                                if( j == k)
                                {
                                    Vmath::Vadd(nq, grad[j*expdim+k], 1,
                                                    div, 1, div, 1);
                                }
                            }
                        }
                        // Scale div by lambda (for compressible flows)
                        Vmath::Smul(nq, lambda, div, 1, div, 1);

                        // Get Coordinates
                        for(int j =0; j<3; ++j)
                        {
                            coords[j] = Array<OneD, NekDouble>(nq, 0.0);
                        }
                        elmt->GetCoords(coords[0], coords[1], coords[2]);

                        // identify boundary of element
                        boundary = m_BCtoTraceID[cnt];

                        // Dimension specific part for obtaining values
                        //   at boundary and normal vector
                        Array<OneD, Array<OneD, NekDouble> > normals = 
                            elmt->GetTraceNormal(boundary);
                        
                        // Get expansion on boundary
                        LocalRegions::ExpansionSharedPtr bc = BndExp[n]->GetExp(i); 

                        // Get number of points on the boundary
                        int nbc = bc->GetTotPoints();

                        // Extract values at boundary
                        Pb = Array<OneD, NekDouble>(nbc,0.0);
                        elmt->GetTracePhysVals(boundary,bc,pElmt,Pb);

                        for(int j = 0; j < expdim*expdim; ++j)
                        {
                            gradb[j] = Array<OneD, NekDouble>(nbc,0.0);
                            elmt->GetTracePhysVals(boundary,bc,grad[j],gradb[j]);
                        }

                        for(int j=0; j<3; ++j)
                        {
                            coordsb[j] = Array<OneD, NekDouble> (nbc, 0.0);
                            elmt->GetTracePhysVals(boundary,bc,coords[j], coordsb[j]);
                        
                            // subtract m_pivotPoint
                            Vmath::Sadd(nbc, -1.0*m_pivotPoint[j],
                                    coordsb[j], 1, coordsb[j], 1);
                        }

                        // Calculate forces per unit length

                        // Pressure component: fp[j] = rho* p*n[j]
                        for (int j = 0; j < expdim; j++)
                        {
                            fp[j] = Array<OneD, NekDouble> (nbc,0.0);
                            Vmath::Vmul (nbc, Pb, 1,
                                              normals[j], 1,
                                              fp[j], 1);
                            Vmath::Smul(nbc, rho, fp[j], 1, fp[j], 1);
                        }

                        // Viscous component:
                        //     fv[j] = -mu*{(grad[k,j]+grad[j,k]) *n[k]}
                        for (int j = 0; j < expdim; j++ )
                        {
                            fv[j] = Array<OneD, NekDouble> (nbc,0.0);
                            for (int k = 0; k < expdim; k++ )
                            {
                                Vmath::Vvtvp (nbc, gradb[k*expdim+j], 1,
                                                   normals[k], 1,
                                                   fv[j], 1,
                                                   fv[j], 1);
                                Vmath::Vvtvp (nbc, gradb[j*expdim+k], 1,
                                                   normals[k], 1,
                                                   fv[j], 1,
                                                   fv[j], 1);
                            }
                            if(!fluidEqu->HasConstantDensity())
                            {
                                // Add gradient term
                                Vmath::Vvtvp (nbc, div, 1,
                                                   normals[j], 1,
                                                   fv[j], 1,
                                                   fv[j], 1);
                            }
                            Vmath::Smul(nbc, -mu, fv[j], 1, fv[j], 1);
                        }

                        // calculate moments per unit length
                        if( momdim == 1)
                        {
                        
                            mp[0] = Array<OneD, NekDouble> (nbc,0.0);
                            mv[0] = Array<OneD, NekDouble> (nbc,0.0);

                            Array<OneD, NekDouble> fptmp(nbc, 0.0);
                            Array<OneD, NekDouble> fvtmp(nbc, 0.0);
                            Vmath::Smul(nbc, -1.0, fp[0], 1, fptmp, 1);
                            Vmath::Smul(nbc, -1.0, fv[0], 1, fvtmp, 1);
                            
                            // Mz = Fy * x - Fx * y
                            Vmath::Vvtvvtp(nbc, fp[1], 1, coordsb[0], 1,
                                       fptmp, 1, coordsb[1], 1,
                                       mp[0], 1);

                            Vmath::Vvtvvtp(nbc, fv[1], 1, coordsb[0], 1,
                                       fvtmp, 1, coordsb[1], 1,
                                       mv[0], 1);
                        }
                        else
                        {
                            Array<OneD, NekDouble> fptmp(nbc, 0.0);
                            Array<OneD, NekDouble> fvtmp(nbc, 0.0);

                            // Mx = Fz * y - Fy * z
                            mp[0] = Array<OneD, NekDouble> (nbc,0.0);
                            mv[0] = Array<OneD, NekDouble> (nbc,0.0);


                            Vmath::Smul(nbc, -1.0, fp[1], 1, fptmp, 1);
                            Vmath::Smul(nbc, -1.0, fv[1], 1, fvtmp, 1);
                            Vmath::Vvtvvtp(nbc, fp[2], 1, coordsb[1], 1,
                                           fptmp, 1, coordsb[2], 1,
                                           mp[0], 1);
                            Vmath::Vvtvvtp(nbc, fv[2], 1, coordsb[1], 1,
                                           fvtmp, 1, coordsb[2], 1,
                                           mv[0], 1);
                            // My = Fx * z - Fz * x
                            mp[1] = Array<OneD, NekDouble> (nbc,0.0);
                            mv[1] = Array<OneD, NekDouble> (nbc,0.0);

                            Vmath::Smul(nbc, -1.0, fp[2], 1, fptmp, 1);
                            Vmath::Smul(nbc, -1.0, fv[2], 1, fvtmp, 1);
                            Vmath::Vvtvvtp(nbc, fp[0], 1, coordsb[2], 1,
                                           fptmp, 1, coordsb[0], 1,
                                           mp[1], 1);
                            Vmath::Vvtvvtp(nbc, fv[0], 1, coordsb[2], 1,
                                           fvtmp, 1, coordsb[0], 1,
                                           mv[1], 1);
                            // Mz = Fy * x - Fx * y
                            mp[2] = Array<OneD, NekDouble> (nbc,0.0);
                            mv[2] = Array<OneD, NekDouble> (nbc,0.0);

                            Vmath::Smul(nbc, -1.0, fp[0], 1, fptmp, 1);
                            Vmath::Smul(nbc, -1.0, fv[0], 1, fvtmp, 1);
                            Vmath::Vvtvvtp(nbc, fp[1], 1, coordsb[0], 1,
                                           fptmp, 1, coordsb[1], 1,
                                           mp[2], 1);
                            Vmath::Vvtvvtp(nbc, fv[1], 1, coordsb[0], 1,
                                           fvtmp, 1, coordsb[1], 1,
                                           mv[2], 1);
                        }



                        // Integrate to obtain force
                        for (int j = 0; j < expdim; j++)
                        {
                            m_Fpplane[j][plane] += BndExp[n]->GetExp(i)->
                                                    Integral(fp[j]);
                            m_Fvplane[j][plane] += BndExp[n]->GetExp(i)->
                                                    Integral(fv[j]);
                        }
                        for ( int j = 0; j < momdim; ++j)
                        {
                            m_Mpplane[j][plane] += BndExp[n]->GetExp(i)->
                                Integral(mp[j]);
                            m_Mvplane[j][plane] += BndExp[n]->GetExp(i)->
                                Integral(mv[j]);
                        }

                    }
                }
                else
                {
                    cnt += BndExp[n]->GetExpSize();
                }
            }
        }
    }

    // Combine contributions from different processes
    //    this is split between row and col comm because of
    //      homostrips case, which only keeps its own strip
    for(int i = 0; i < expdim; i++)
    {
        rowComm->AllReduce(m_Fpplane[i], LibUtilities::ReduceSum);
        colComm->AllReduce(m_Fpplane[i], LibUtilities::ReduceSum);

        rowComm->AllReduce(m_Fvplane[i], LibUtilities::ReduceSum);
        colComm->AllReduce(m_Fvplane[i], LibUtilities::ReduceSum);
    }
    for( int i = 0; i < momdim; ++i)
    {
        rowComm->AllReduce(m_Mpplane[i], LibUtilities::ReduceSum);
        colComm->AllReduce(m_Mpplane[i], LibUtilities::ReduceSum);

        rowComm->AllReduce(m_Mvplane[i], LibUtilities::ReduceSum);
        colComm->AllReduce(m_Mvplane[i], LibUtilities::ReduceSum);
    }
    

    // Project results to new directions
    for(int plane = 0; plane < m_nPlanes; plane++)
    {
        Array< OneD, NekDouble> tmpP(expdim, 0.0);
        Array< OneD, NekDouble> tmpV(expdim, 0.0);
        for(int i = 0; i < expdim; i++)
        {
            for(int j = 0; j < expdim; j++ )
            {
                tmpP[i] += m_Fpplane[j][plane]*m_directions[i][j];
                tmpV[i] += m_Fvplane[j][plane]*m_directions[i][j];
            }
        }
        // Copy result
        for(int i = 0; i < expdim; i++)
        {
            m_Fpplane[i][plane] = tmpP[i];
            m_Fvplane[i][plane] = tmpV[i];
        }
        
        // Project moments only in 3D, since 2D moment is always in z direction
        if (momdim == 3)
        {
            for( int i = 0; i < 3; ++i)
            {
                tmpP[i] = 0.0;
                tmpV[i] = 0.0;
                for( int j = 0; j < 3; ++j )
                {
                    tmpP[i] += m_Mpplane[j][plane]*m_directions[i][j];
                    tmpV[i] += m_Mvplane[j][plane]*m_directions[i][j];
                }
            }
            // Copy result
            for( int i = 0; i < 3; ++i)
            {
                m_Mpplane[i][plane] = tmpP[i];
                m_Mvplane[i][plane] = tmpV[i];
            }
        }

    }

    // Sum viscous and pressure components
    for(int plane = 0; plane < m_nPlanes; plane++)
    {
        for(int i = 0; i < expdim; i++)
        {
            m_Ftplane[i][plane] = m_Fpplane[i][plane] + m_Fvplane[i][plane];
        }

        if( momdim == 3 )
        {           
            for(int i = 0; i < 3; ++i)
            {
                m_Mtplane[i][plane] = m_Mpplane[i][plane] + m_Mvplane[i][plane];
            }
        }
        else
        {
            m_Mtplane[0][plane] = m_Mpplane[0][plane] + m_Mvplane[0][plane];

        }

    }

    // combine planes
     for( int i = 0; i < expdim; i++)
    {
        m_Fp[i] = Vmath::Vsum(m_nPlanes, m_Fpplane[i], 1) / m_nPlanes;
        m_Fv[i] = Vmath::Vsum(m_nPlanes, m_Fvplane[i], 1) / m_nPlanes;
        m_Ft[i] = m_Fp[i] + m_Fv[i];
    }
    for( int i = 0; i < momdim; i++)
    {
        m_Mp[i] = Vmath::Vsum(m_nPlanes, m_Mpplane[i], 1) / m_nPlanes;
        m_Mv[i] = Vmath::Vsum(m_nPlanes, m_Mvplane[i], 1) / m_nPlanes;
        m_Mt[i] = m_Mp[i] + m_Mv[i];
    }


    // Put results back to wavespace, if necessary
    if( m_isHomogeneous1D && m_outputAllPlanes )
    {
        for (int i = 0; i < pFields.size(); ++i)
        {
            pFields[i]->SetWaveSpace(true);
            pFields[i]->HomogeneousFwdTrans(pFields[i]->GetPhys(),
                                            pFields[i]->UpdatePhys());
        }
    }
}

/**
 *     This function calculates the forces when we have a mapping
 *         defining a coordinate system transformation
 */
void FilterAeroForces::CalculateForcesMapping(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    boost::ignore_unused(time);

    int cnt, elmtid, offset, boundary;
    // Get number of quadrature points and dimensions
    int physTot = pFields[0]->GetNpoints();
    int expdim = pFields[0]->GetGraph()->GetMeshDimension();
    int momdim = ( expdim == 2 ) ? 1 : 3;
    int nVel = expdim;
    if( m_isHomogeneous1D )
    {
        nVel = nVel + 1;
    }

    LocalRegions::ExpansionSharedPtr elmt;

    // Pressure stress tensor
    //    (global, in a plane, in element and boundary)
    Array<OneD, MultiRegions::ExpListSharedPtr>  P      ( nVel*nVel);
    Array<OneD, MultiRegions::ExpListSharedPtr>  PPlane ( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         PElmt  ( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         PBnd  ( nVel*nVel);
    // Velocity gradient
    Array<OneD, MultiRegions::ExpListSharedPtr>  grad     ( nVel*nVel);
    Array<OneD, MultiRegions::ExpListSharedPtr>  gradPlane( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         gradElmt ( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         gradBnd  ( nVel*nVel);

    // Coordinates
    Array<OneD, MultiRegions::ExpListSharedPtr>  coords      (3);
    Array<OneD, MultiRegions::ExpListSharedPtr>  coordsPlane (3);
    Array<OneD, Array<OneD, NekDouble> >         coordsElmt  (3);
    Array<OneD, Array<OneD, NekDouble> >         coordsBnd   (3);

    // Transformation to cartesian system
    Array<OneD, MultiRegions::ExpListSharedPtr>  C     ( nVel*nVel);
    Array<OneD, MultiRegions::ExpListSharedPtr>  CPlane( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         CElmt ( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         CBnd  ( nVel*nVel);

    // Jacobian
    MultiRegions::ExpListSharedPtr  Jac;
    MultiRegions::ExpListSharedPtr  JacPlane;
    Array<OneD, NekDouble>          JacElmt;
    Array<OneD, NekDouble>          JacBnd;

    // Communicators to exchange results
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
    LibUtilities::CommSharedPtr rowComm = vComm->GetRowComm();
    LibUtilities::CommSharedPtr colComm =
                            m_session->DefinesSolverInfo("HomoStrip") ?
                                vComm->GetColumnComm()->GetColumnComm():
                                vComm->GetColumnComm();

    // Arrays to store the plane average forces in each direction
    m_Fp = Array<OneD, NekDouble> (expdim,0.0);
    m_Fv = Array<OneD, NekDouble> (expdim,0.0);
    m_Ft = Array<OneD, NekDouble> (expdim,0.0);
    m_Mp = Array<OneD, NekDouble> (momdim,0.0);
    m_Mv = Array<OneD, NekDouble> (momdim,0.0);
    m_Mt = Array<OneD, NekDouble> (momdim,0.0);


    // Arrays with forces in each plane
    m_Fpplane = Array<OneD, Array<OneD, NekDouble> >  (expdim);
    m_Fvplane = Array<OneD, Array<OneD, NekDouble> >  (expdim);
    m_Ftplane = Array<OneD, Array<OneD, NekDouble> >  (expdim);
    for(int i = 0; i < expdim; i++)
    {
        m_Fpplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
        m_Fvplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
        m_Ftplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
    }

    // Arrays with moments in each plane
    m_Mpplane = Array<OneD, Array<OneD, NekDouble> >  (momdim);
    m_Mvplane = Array<OneD, Array<OneD, NekDouble> >  (momdim);
    m_Mtplane = Array<OneD, Array<OneD, NekDouble> >  (momdim);
    for( int i = 0; i < momdim; ++i)
    {
        m_Mpplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
        m_Mvplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
        m_Mtplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
    }


    // Forces per element length in a boundary
    Array<OneD, Array<OneD, NekDouble> >       fp( nVel );
    Array<OneD, Array<OneD, NekDouble> >       fv( nVel );

    // moments per element length in a boundary
    Array<OneD, Array<OneD, NekDouble> >       mp( nVel );
    Array<OneD, Array<OneD, NekDouble> >       mv( nVel );


    // Get viscosity
    NekDouble mu{1}, rho{1.};
    if(m_session->DefinesParameter("Kinvis"))
    {
        rho = (m_session->DefinesParameter("rho"))
            ? (m_session->GetParameter("rho")) : 1;
        mu = rho * m_session->GetParameter("Kinvis");
    }
    else
    {
        mu = m_session->GetParameter("mu");
    }
    //NekDouble lambda = -2.0/3.0;

    // Perform BwdTrans: for case with mapping, we can never work
    //                   in wavespace
    for(int i = 0; i < pFields.size(); ++i)
    {
        if (m_isHomogeneous1D)
        {
            pFields[i]->SetWaveSpace(false);
        }
        pFields[i]->BwdTrans(pFields[i]->GetCoeffs(),
                             pFields[i]->UpdatePhys());
        pFields[i]->SetPhysState(true);
    }

    // Define boundary expansions
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
    Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
    if(m_isHomogeneous1D)
    {
        BndConds = pFields[0]->GetPlane(0)->GetBndConditions();
        BndExp   = pFields[0]->GetPlane(0)->GetBndCondExpansions();
    }
    else
    {
        BndConds = pFields[0]->GetBndConditions();
        BndExp   = pFields[0]->GetBndCondExpansions();
    }

    //
    // Calculate pressure stress tensor, velocity gradient
    //      and get informations about the mapping

    // Initialise variables
    switch (expdim)
    {
        case 2:
        {
            if (m_isHomogeneous1D)
            {
                MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;
                Exp3DH1 = std::dynamic_pointer_cast
                                <MultiRegions::ExpList3DHomogeneous1D>
                                                    (pFields[0]);
                for(int i = 0; i < nVel*nVel; i++)
                {
                    grad[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*Exp3DH1);

                    P[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*Exp3DH1);

                    C[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*Exp3DH1);
                }

                for(int i = 0; i < 3; ++i)
                {
                    coords[i] = 
                        MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*Exp3DH1);
                }
                Jac = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*Exp3DH1);
            }
            else
            {
                MultiRegions::ExpListSharedPtr Exp2D;
                Exp2D = std::dynamic_pointer_cast
                                <MultiRegions::ExpList>
                                                    (pFields[0]);
                for(int i = 0; i < nVel*nVel; i++)
                {
                    grad[i] = MemoryManager<MultiRegions::ExpList>::
                                AllocateSharedPtr(*Exp2D);

                    P[i] = MemoryManager<MultiRegions::ExpList>::
                                AllocateSharedPtr(*Exp2D);

                    C[i] = MemoryManager<MultiRegions::ExpList>::
                                AllocateSharedPtr(*Exp2D);
                }
                for(int i = 0; i < 3; ++i)
                {
                    coords[i] = 
                        MemoryManager<MultiRegions::ExpList>::
                                AllocateSharedPtr(*Exp2D);
                }
                Jac = MemoryManager<MultiRegions::ExpList>::
                                AllocateSharedPtr(*Exp2D);
            }
            break;
        }
        case 3:
        {
            MultiRegions::ExpListSharedPtr Exp3D;
            Exp3D = std::dynamic_pointer_cast
                            <MultiRegions::ExpList>
                                                (pFields[0]);
            for(int i = 0; i < nVel*nVel; i++)
            {
                grad[i] = MemoryManager<MultiRegions::ExpList>::
                            AllocateSharedPtr(*Exp3D);

                P[i] = MemoryManager<MultiRegions::ExpList>::
                            AllocateSharedPtr(*Exp3D);

                C[i] = MemoryManager<MultiRegions::ExpList>::
                            AllocateSharedPtr(*Exp3D);
            }
            for(int i = 0; i < 3; ++i )
            {
                coords[i] = 
                    MemoryManager<MultiRegions::ExpList>::
                            AllocateSharedPtr(*Exp3D);
            }
            Jac = MemoryManager<MultiRegions::ExpList>::
                            AllocateSharedPtr(*Exp3D);

            break;
        }
        default:
            ASSERTL0(false,"Expansion dimension not supported by FilterAeroForces");
            break;
    }

    // Get g^ij
    Array<OneD, Array<OneD, NekDouble> >        tmp( nVel*nVel );
    m_mapping->GetInvMetricTensor(tmp);

    // Get Cartesian coordinates
    m_mapping->GetCartesianCoordinates(coords[0]->UpdatePhys(),
            coords[1]->UpdatePhys(), coords[2]->UpdatePhys());

    // Calculate P^ij = g^ij*p
    for (int i = 0; i < nVel*nVel; ++i)
    {
        Vmath::Vmul(physTot, tmp[i], 1,
                            pFields[nVel]->GetPhys(), 1,
                            P[i]->UpdatePhys(), 1);
    }

    // Calculate covariant derivatives of U = u^i_,k
    Array<OneD, Array<OneD, NekDouble> >        wk( nVel );
    for (int i=0; i<nVel; ++i)
    {
        wk[i] = Array<OneD, NekDouble>(physTot, 0.0);
        Vmath::Vcopy(physTot, pFields[i]->GetPhys(), 1,
                            wk[i], 1);
    }
    m_mapping->ApplyChristoffelContravar(wk, tmp);
    for (int i=0; i< nVel; ++i)
    {
        for (int k=0; k< nVel; ++k)
        {
            pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[k],
                                   wk[i], grad[i*nVel+k]->UpdatePhys());

            Vmath::Vadd(physTot,tmp[i*nVel+k],1,
                                grad[i*nVel+k]->UpdatePhys(), 1,
                                grad[i*nVel+k]->UpdatePhys(), 1);
        }
    }
    // Raise index to obtain Grad^ij = g^jk u^i_,k
    for (int i=0; i< nVel; ++i)
    {
        for (int k=0; k< nVel; ++k)
        {
            Vmath::Vcopy(physTot, grad[i*nVel+k]->GetPhys(), 1,
                                  wk[k], 1);
        }
        m_mapping->RaiseIndex(wk, wk);
        for (int j=0; j<nVel; ++j)
        {
            Vmath::Vcopy(physTot, wk[j], 1,
                                  grad[i*nVel+j]->UpdatePhys(), 1);
        }
    }

    // Get Jacobian
    m_mapping->GetJacobian( Jac->UpdatePhys());

    // Get transformation to Cartesian system
    for (int i=0; i< nVel; ++i)
    {
        // Zero wk storage
        for (int k=0; k< nVel; ++k)
        {
            wk[k] = Array<OneD, NekDouble>(physTot, 0.0);
        }
        // Make wk[i] = 1
        wk[i] = Array<OneD, NekDouble>(physTot, 1.0);
        // Transform wk to Cartesian
        m_mapping->ContravarToCartesian(wk,wk);
        // Copy result to a column in C
        for (int k=0; k< nVel; k++)
        {
            Vmath::Vcopy(physTot, wk[k], 1,
                                  C[k*nVel+i]->UpdatePhys(), 1);
        }
    }

    //
    // Calculate forces
    //

    // For Homogeneous, calculate force on each 2D plane
    // Otherwise, m_nPlanes = 1, and loop only runs once
    for(int plane = 0; plane < m_nPlanes; ++plane )
    {
        // Check if plane is in this proc
        if( m_planesID[plane] != -1 )
        {
            // For Homogeneous, consider the 2D expansion
            //      on the current plane
            if(m_isHomogeneous1D)
            {
                for(int n = 0; n < nVel*nVel; ++n)
                {
                    PPlane[n]    = P[n]->GetPlane(m_planesID[plane]);
                    gradPlane[n] = grad[n]->GetPlane(m_planesID[plane]);
                    CPlane[n]    = C[n]->GetPlane(m_planesID[plane]);
                }
                for(int n = 0; n < 3; ++n)
                {
                    coordsPlane[n] = coords[n]->GetPlane(m_planesID[plane]);
                }
                JacPlane = Jac->GetPlane(m_planesID[plane]);
            }
            else
            {
                for(int n = 0; n < nVel*nVel; ++n)
                {
                    PPlane[n]    = P[n];
                    gradPlane[n] = grad[n];
                    CPlane[n] = C[n];
                }
                for(int n = 0; n < 3; ++n)
                {
                    coordsPlane[n] = coords[n];
                }
                JacPlane = Jac;
            }

            //Loop all the Boundary Regions
            for( int n = cnt = 0; n < BndConds.size(); n++)
            {
                if(m_boundaryRegionIsInList[n] == 1)
                {
                    for (int i=0; i < BndExp[n]->GetExpSize(); ++i,cnt++)
                    {
                        elmtid = m_BCtoElmtID[cnt];
                        elmt   = PPlane[0]->GetExp(elmtid);
                        offset = PPlane[0]->GetPhys_Offset(elmtid);

                        // Extract  fields on this element
                        for(int j=0; j<nVel*nVel; j++)
                        {
                            PElmt[j]    = PPlane[j]->GetPhys()
                                        + offset;
                            gradElmt[j] = gradPlane[j]->GetPhys()
                                        + offset;
                            CElmt[j]    = CPlane[j]->GetPhys()
                                        + offset;
                        }
                        for(int j = 0; j < 3; ++j)
                        {
                            coordsElmt[j] = coordsPlane[j]->GetPhys() + offset;
                        }
                        JacElmt = JacPlane->GetPhys() + offset;

                        // identify boundary of element
                        boundary = m_BCtoTraceID[cnt];

                        // Dimension specific part for obtaining values
                        //   at boundary and normal vector
                        Array<OneD, Array<OneD, NekDouble> > normals;
                        // Get normals
                        normals = elmt->GetTraceNormal(boundary);

                        // Get expansion on boundary
                        LocalRegions::ExpansionSharedPtr bc = 
                            BndExp[n]->GetExp(i);
                        
                        // Get number of points on the boundary
                        int nbc = bc->GetTotPoints();

                        // Extract values at boundary
                        for(int j = 0; j < nVel*nVel; ++j)
                        {
                            gradBnd[j] = Array<OneD, NekDouble>(nbc,0.0);
                            elmt->GetTracePhysVals(boundary,bc,gradElmt[j],gradBnd[j]);
                            
                            PBnd[j] = Array<OneD, NekDouble>(nbc,0.0);
                            elmt->GetTracePhysVals(boundary,bc,PElmt[j],PBnd[j]);

                            CBnd[j] = Array<OneD, NekDouble>(nbc,0.0);
                            elmt->GetTracePhysVals(boundary,bc,CElmt[j],CBnd[j]);
                        }
                        for(int j = 0; j < 3; ++j)
                        {
                            coordsBnd[j] = Array<OneD, NekDouble>(nbc, 0.0);
                            elmt->GetTracePhysVals(boundary,
                                    bc, coordsElmt[j], coordsBnd[j]);

                            // substract the m_pivotPoint
                            Vmath::Sadd(nbc, -1.0*m_pivotPoint[j],
                                    coordsBnd[j], 1, coordsBnd[j], 1);
                        }
                        JacBnd = Array<OneD, NekDouble>(nbc,0.0);
                        elmt->GetTracePhysVals(boundary,bc,JacElmt,JacBnd);

                        // Calculate forces per unit length

                        // Pressure component: fp[j] = rho*P[j,k]*n[k]
                        for ( int j = 0; j < nVel; j++)
                        {
                            fp[j] = Array<OneD, NekDouble> (nbc,0.0);
                            // Normals only has expdim elements
                            for (int k = 0; k < expdim; k++)
                            {
                                Vmath::Vvtvp (nbc, PBnd[ j*nVel + k], 1,
                                                   normals[k], 1,
                                                   fp[j], 1,
                                                   fp[j], 1);
                            }
                            Vmath::Smul(nbc, rho, fp[j], 1, fp[j], 1);
                        }

                        // Viscous component:
                        //     fv[j] = -mu*{(grad[k,j]+grad[j,k]) *n[k]}
                        for (int j = 0; j < nVel; j++ )
                        {
                            fv[j] = Array<OneD, NekDouble> (nbc,0.0);
                            for (int k = 0; k < expdim; k++ )
                            {
                                Vmath::Vvtvp (nbc,gradBnd[k*nVel+j],1,
                                                   normals[k], 1,
                                                   fv[j], 1,
                                                   fv[j], 1);
                                Vmath::Vvtvp (nbc,gradBnd[j*nVel+k],1,
                                                   normals[k], 1,
                                                   fv[j], 1,
                                                   fv[j], 1);
                            }
                            Vmath::Smul(nbc, -mu, fv[j], 1, fv[j], 1);
                        }

                        // Multiply by Jacobian
                        for (int k = 0; k < nVel; k++ )
                        {
                            Vmath::Vmul(nbc, JacBnd, 1, fp[k], 1,
                                                        fp[k], 1);
                            Vmath::Vmul(nbc, JacBnd, 1, fv[k], 1,
                                                        fv[k], 1);
                        }

                        // Convert to cartesian system
                        for (int k = 0; k < nVel; k++ )
                        {
                            wk[k] = Array<OneD, NekDouble>(physTot,0.0);
                            for (int j = 0; j < nVel; j++ )
                            {
                                Vmath::Vvtvp(nbc, CBnd[k*nVel+j], 1,
                                                    fp[j], 1,
                                                    wk[k], 1,
                                                    wk[k], 1);
                            }
                        }
                        for (int k = 0; k < nVel; k++ )
                        {
                            Vmath::Vcopy(nbc, wk[k], 1, fp[k], 1);
                        }

                        for (int k = 0; k < nVel; k++ )
                        {
                            wk[k] = Array<OneD, NekDouble>(physTot,0.0);
                            for (int j = 0; j < nVel; j++ )
                            {
                                Vmath::Vvtvp(nbc, CBnd[k*nVel+j], 1,
                                                    fv[j], 1,
                                                    wk[k], 1,
                                                    wk[k], 1);
                            }
                        }
                        for (int k = 0; k < nVel; k++ )
                        {
                            Vmath::Vcopy(nbc, wk[k], 1, fv[k], 1);
                        }

                        // Calculate moments per unit length
                        if( momdim == 1 )
                        {
                            Array<OneD, NekDouble> fptmp(nbc, 0.0);
                            Array<OneD, NekDouble> fvtmp(nbc, 0.0);

                            mp[0] = Array<OneD, NekDouble> (nbc, 0.0);
                            mv[0] = Array<OneD, NekDouble> (nbc, 0.0);

                            Vmath::Smul(nbc, -1.0, fp[0], 1, fptmp, 1);
                            Vmath::Smul(nbc, -1.0, fv[0], 1, fvtmp, 1);
                            // Mz = Fy * x - Fx * y
                            Vmath::Vvtvvtp(nbc,
                                    fp[1], 1, coordsBnd[0], 1,
                                    fptmp, 1, coordsBnd[1], 1,
                                    mp[0], 1);
                            Vmath::Vvtvvtp(nbc,
                                    fv[1], 1, coordsBnd[0], 1,
                                    fvtmp, 1, coordsBnd[1], 1,
                                    mv[0], 1);
                        }
                        else
                        {
                            Array<OneD, NekDouble> fptmp(nbc, 0.0);
                            Array<OneD, NekDouble> fvtmp(nbc, 0.0);

                            // Mx = Fz * y - Fy * z
                            mp[0] = Array<OneD, NekDouble> (nbc, 0.0);
                            mv[0] = Array<OneD, NekDouble> (nbc, 0.0);

                            Vmath::Smul(nbc, -1.0, fp[1], 1, fptmp, 1);
                            Vmath::Smul(nbc, -1.0, fv[1], 1, fvtmp, 1);
                            Vmath::Vvtvvtp(nbc,
                                    fp[2], 1, coordsBnd[1], 1, 
                                    fptmp, 1, coordsBnd[2], 1,
                                    mp[0], 1);
                            Vmath::Vvtvvtp(nbc,
                                    fv[2], 1, coordsBnd[1], 1, 
                                    fvtmp, 1, coordsBnd[2], 1,
                                    mv[0], 1);
                            // My = Fx * z - Fz * x
                            mp[1] = Array<OneD, NekDouble> (nbc, 0.0);
                            mv[1] = Array<OneD, NekDouble> (nbc, 0.0);

                            Vmath::Smul(nbc, -1.0, fp[2], 1, fptmp, 1);
                            Vmath::Smul(nbc, -1.0, fv[2], 1, fvtmp, 1);
                            Vmath::Vvtvvtp(nbc,
                                    fp[0], 1, coordsBnd[2], 1, 
                                    fptmp, 1, coordsBnd[0], 1,
                                    mp[1], 1);
                            Vmath::Vvtvvtp(nbc,
                                    fv[0], 1, coordsBnd[2], 1, 
                                    fvtmp, 1, coordsBnd[0], 1,
                                    mv[1], 1);
                            // Mz = Fy * x - Fx * y
                            mp[2] = Array<OneD, NekDouble> (nbc, 0.0);
                            mv[2] = Array<OneD, NekDouble> (nbc, 0.0);

                            Vmath::Smul(nbc, -1.0, fp[0], 1, fptmp, 1);
                            Vmath::Smul(nbc, -1.0, fv[0], 1, fvtmp, 1);
                            Vmath::Vvtvvtp(nbc,
                                    fp[1], 1, coordsBnd[0], 1, 
                                    fptmp, 1, coordsBnd[1], 1,
                                    mp[2], 1);
                            Vmath::Vvtvvtp(nbc,
                                    fv[1], 1, coordsBnd[0], 1, 
                                    fvtmp, 1, coordsBnd[1], 1,
                                    mv[2], 1);
                        }
                            


                        // Integrate to obtain force
                        for (int j = 0; j < expdim; j++)
                        {
                            m_Fpplane[j][plane] += BndExp[n]->GetExp(i)->
                                                    Integral(fp[j]);
                            m_Fvplane[j][plane] += BndExp[n]->GetExp(i)->
                                                    Integral(fv[j]);
                        }

                        // Integrate to obtain moments
                        for(int j = 0; j < momdim; ++j)
                        {
                            m_Mpplane[j][plane] += BndExp[n]->GetExp(i)->
                                                     Integral(mp[j]);
                            m_Mvplane[j][plane] += BndExp[n]->GetExp(i)->
                                                     Integral(mv[j]);
                        }
                    }
                }
                else
                {
                    cnt += BndExp[n]->GetExpSize();
                }
            }
        }
    }

    // Combine contributions from different processes
    //    this is split between row and col comm because of
    //      homostrips case, which only keeps its own strip
    for(int i = 0; i < expdim; ++i)
    {
        rowComm->AllReduce(m_Fpplane[i], LibUtilities::ReduceSum);
        colComm->AllReduce(m_Fpplane[i], LibUtilities::ReduceSum);

        rowComm->AllReduce(m_Fvplane[i], LibUtilities::ReduceSum);
        colComm->AllReduce(m_Fvplane[i], LibUtilities::ReduceSum);
    }
    for(int i = 0; i < momdim; ++i)
    {
        rowComm->AllReduce(m_Mpplane[i], LibUtilities::ReduceSum);
        colComm->AllReduce(m_Mpplane[i], LibUtilities::ReduceSum);

        rowComm->AllReduce(m_Mvplane[i], LibUtilities::ReduceSum);
        colComm->AllReduce(m_Mvplane[i], LibUtilities::ReduceSum);
    }

    // Project results to new directions
    for(int plane = 0; plane < m_nPlanes; plane++)
    {
        Array< OneD, NekDouble> tmpP(expdim, 0.0);
        Array< OneD, NekDouble> tmpV(expdim, 0.0);
        for(int  i = 0; i < expdim; i++)
        {
            for(int  j = 0; j < expdim; j++ )
            {
                tmpP[i] += m_Fpplane[j][plane]*m_directions[i][j];
                tmpV[i] += m_Fvplane[j][plane]*m_directions[i][j];
            }
        }
        // Copy result
        for(int i = 0; i < expdim; i++)
        {
            m_Fpplane[i][plane] = tmpP[i];
            m_Fvplane[i][plane] = tmpV[i];
        }
        // 
        // For 3D moment, project the results, note in 2D case the moment is always in z direction
        if( momdim == 3)
        {
            for(int i = 0; i < 3; ++i)
            {
                tmpP[i] = 0.0;
                tmpV[i] = 0.0;
                for(int j = 0; j < 3; ++j)
                {
                    tmpP[i] += m_Mpplane[j][plane]*m_directions[i][j];
                    tmpV[i] += m_Mvplane[j][plane]*m_directions[i][j];
                }
            }

            // copy the result
            for(int i = 0; i < 3; ++i)
            {
                m_Mpplane[i][plane] = tmpP[i];
                m_Mvplane[i][plane] = tmpV[i];
            }
        }
    }

    // Sum viscous and pressure components
    for(int plane = 0; plane < m_nPlanes; plane++)
    {
        for(int  i = 0; i < expdim; i++)
        {
            m_Ftplane[i][plane] = m_Fpplane[i][plane] + m_Fvplane[i][plane];
        }

        if( momdim == 3 )
        {
            for(int i = 0; i < 3; ++i)
            {
                m_Mtplane[i][plane] = m_Mpplane[i][plane] + m_Mvplane[i][plane];
            }
        }
        else
        {
            m_Mtplane[0][plane] = m_Mpplane[0][plane] + m_Mvplane[0][plane];
        }

    }
    
    // combine planes
    for( int i = 0; i < expdim; i++)
    {
        m_Fp[i] = Vmath::Vsum(m_nPlanes, m_Fpplane[i], 1) / m_nPlanes;
        m_Fv[i] = Vmath::Vsum(m_nPlanes, m_Fvplane[i], 1) / m_nPlanes;
        m_Ft[i] = m_Fp[i] + m_Fv[i];
    }
    for( int i = 0; i < momdim; i++)
    {
        m_Mp[i] = Vmath::Vsum(m_nPlanes, m_Mpplane[i], 1) / m_nPlanes;
        m_Mv[i] = Vmath::Vsum(m_nPlanes, m_Mvplane[i], 1) / m_nPlanes;
        m_Mt[i] = m_Mp[i] + m_Mv[i];
    }

    // Put results back to wavespace, if necessary
    if( m_isHomogeneous1D)
    {
        for (int i = 0; i < pFields.size(); ++i)
        {
            pFields[i]->SetWaveSpace(true);
            pFields[i]->HomogeneousFwdTrans(pFields[i]->GetPhys(),
                                            pFields[i]->UpdatePhys());
        }
    }
}

}
}
