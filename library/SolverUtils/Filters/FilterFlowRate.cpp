///////////////////////////////////////////////////////////////////////////////
//
// File FilterFlowRate.cpp
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
// Description: Output values of flow rate on specific boundaries during time-stepping.
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
#include <SolverUtils/Filters/FilterFlowRate.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
std::string FilterFlowRate::className =
        GetFilterFactory().RegisterCreatorFunction(
                "FlowRate", FilterFlowRate::create);

/**
 *
 */
FilterFlowRate::FilterFlowRate(
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
          && m_outputFile.substr(m_outputFile.length() - 4) == ".rate"))
    {
        m_outputFile += ".rate";
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
}


/**
 *
 */
FilterFlowRate::~FilterFlowRate()
{

}

/**
 *
 */
void FilterFlowRate::v_Initialise(
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
              "range for FilterFlowRate: ") + IndString).c_str());

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
        m_outputStream << "# Flow rate across the bodies" << endl;
        
        m_outputStream << "# Boundary regions: " << IndString.c_str() << endl;
        m_outputStream << "#";
        m_outputStream.width(7);
        m_outputStream << "Time";
        m_outputStream.width(9);
        m_outputStream <<  "FlowRate";
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
void FilterFlowRate::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Only output every m_outputFrequency.
    if ((m_index++) % m_outputFrequency  || (time < m_startTime))
    {
        return;
    }
    // Calculate the forces
    CalculateForces(pFields, time);

    // Calculate forces including all planes
    NekDouble Fv = Vmath::Vsum(m_nPlanes, m_FratePlane, 1) / m_nPlanes;

    // Communicators to exchange results
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

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
                m_outputStream.width(15);
                m_outputStream << setprecision(8)
                                << m_FratePlane[plane];
                m_outputStream.width(10);
                m_outputStream << plane;
                m_outputStream << endl;
            }
        }
        // Output average (or total) force
        m_outputStream.width(8);
        m_outputStream << setprecision(6) << time;
        m_outputStream.width(15);
        m_outputStream << setprecision(8) << Fv;
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
            int nstrips;
            m_session->LoadParameter("Strip_Z", nstrips);
            for(int i = 1; i<nstrips; i++)
            {
                vComm->GetColumnComm()->Recv(i, Fv);

                m_outputStream.width(8);
                m_outputStream << setprecision(6) << time;
                m_outputStream.width(15);
                m_outputStream << setprecision(8) << Fv;
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
                vComm->GetColumnComm()->Send(0, Fv);
        }
    }

}


/**
 *
 */
void FilterFlowRate::v_Finalise(
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
bool FilterFlowRate::v_IsTimeDependent()
{
    return true;
}

/**
 *     This function outputs the force on all planes of the current
 *          process, in the format required by ForcingMovingBody
 */
void FilterFlowRate::GetForces(
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
                
            }
            Aeroforces[plane] =
                        m_FratePlane[ZIDs[plane]];
        }
    }
    else
    {
        for(int plane = 0 ; plane < local_planes; plane++)
        {
            Aeroforces[plane] =
                    m_FratePlane[0];
        }
    }
}

/**
 *     This function calculates the forces
 */
void FilterFlowRate::CalculateForces(
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

    int i, j, n, cnt, elmtid, offset, boundary, plane;
    // Get number of quadrature points and dimensions
    int expdim = pFields[0]->GetGraph()->GetMeshDimension();
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
    Array<OneD, NekDouble> density;

    // Arrays of variables in the element
    Array<OneD, Array<OneD, NekDouble> >       velElmt(expdim);

    // Velocity gradient
    Array<OneD, Array<OneD, NekDouble> >       grad( expdim*expdim);

    // Values at the boundary
    Array<OneD, Array<OneD, NekDouble> >       gradb( expdim);

    // Communicators to exchange results
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
    LibUtilities::CommSharedPtr rowComm = vComm->GetRowComm();
    LibUtilities::CommSharedPtr colComm =
                            m_session->DefinesSolverInfo("HomoStrip") ?
                                vComm->GetColumnComm()->GetColumnComm():
                                vComm->GetColumnComm();

    // Arrays with forces in each plane
    m_FratePlane = Array<OneD, NekDouble>(m_nPlanes,0.0);

    // Forces per element length in a boundary
    Array<OneD, NekDouble>       fv;

    // Perform BwdTrans: when we only want the mean force in a 3DH1D
    //     we work in wavespace, otherwise we use physical space
    for(i = 0; i < pFields.size(); ++i)
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
    for(plane = 0; plane < m_nPlanes; plane++ )
    {
        // Check if plane is in this proc
        if( m_planesID[plane] != -1 )
        {
            // For Homogeneous, consider the 2D expansion
            //      on the current plane
            if(m_isHomogeneous1D)
            {
                for(n = 0; n < pFields.size(); n++)
                {
                    fields[n] = pFields[n]->GetPlane(m_planesID[plane]);
                }
            }
            else
            {
                for(n = 0; n < pFields.size(); n++)
                {
                    fields[n] = pFields[n];
                }
            }

            // Get velocity and pressure values
            for(n = 0; n < physfields.size(); ++n)
            {
                physfields[n] = fields[n]->GetPhys();
            }
            for(n = 0; n < nVel; ++n)
            {
                velocity[n] = Array<OneD, NekDouble>(fields[n]->GetTotPoints());
            }
            density = Array<OneD, NekDouble>(fields[n]->GetTotPoints());
            fluidEqu->GetVelocity(physfields, velocity);
            fluidEqu->GetDensity(physfields, density);

            for(n = 0; n < nVel; ++n)
            {
                Vmath::Vmul(fields[n]->GetTotPoints(), density, 1, 
                            velocity[n], 1,
                            velocity[n], 1);
            }

            //Loop all the Boundary Regions
            for( cnt = n = 0; n < BndConds.size(); n++)
            {
                if(m_boundaryRegionIsInList[n] == 1)
                {
                    for (i=0; i < BndExp[n]->GetExpSize(); ++i,cnt++)
                    {
                        elmtid = m_BCtoElmtID[cnt];
                        elmt   = fields[0]->GetExp(elmtid);
                        offset = fields[0]->GetPhys_Offset(elmtid);

                        // Extract  fields on this element
                        for( j=0; j<expdim; j++)
                        {
                            velElmt[j] = velocity[j] + offset;
                        }

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

                        for(j = 0; j < expdim; ++j)
                        {
                            gradb[j] = Array<OneD, NekDouble>(nbc,0.0);
                            elmt->GetTracePhysVals(boundary,bc,velElmt[j],gradb[j]);
                        }

                        // Calculate forces per unit length
                        fv = Array<OneD, NekDouble> (nbc,0.0);
                        for ( j = 0; j < expdim; j++ )
                        {

                            Vmath::Vvtvp (nbc, gradb[j], 1,
                                                normals[j], 1,
                                                fv, 1,
                                                fv, 1);
                        }

                        // Integrate to obtain force
                        m_FratePlane[plane] += BndExp[n]->GetExp(i)->
                                                Integral(fv);
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
    rowComm->AllReduce(m_FratePlane, LibUtilities::ReduceSum);
    colComm->AllReduce(m_FratePlane, LibUtilities::ReduceSum);

    // Put results back to wavespace, if necessary
    if( m_isHomogeneous1D && m_outputAllPlanes )
    {
        for (i = 0; i < pFields.size(); ++i)
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
void FilterFlowRate::CalculateForcesMapping(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    boost::ignore_unused(pFields, time);

    NEKERROR(ErrorUtil::efatal, "CalculateForcesMapping Not implemented.");
}

}
}
