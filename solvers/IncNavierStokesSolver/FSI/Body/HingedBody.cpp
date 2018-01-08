///////////////////////////////////////////////////////////////////////////////
//
// File: HingedBody.cpp
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
// Description: Rigid body motion for FSI problems
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/FSI/Body/HingedBody.h>
#include <LibUtilities/BasicUtils/Equation.h>
#include <GlobalMapping/Mapping.h>

namespace Nektar
{

std::string HingedBody::className =
    GetFSIBodyFactory().RegisterCreatorFunction("Hinged",
    HingedBody::create, "Hinged body motion");

NekDouble HingedBody::AdamsBashforth_coeffs[3][3] = {
        { 1.0       , 0.0       , 0.0     },
        { 3.0/2.0   ,-1.0/2.0   , 0.0     },
        { 23.0/12.0 ,-4.0/3.0   , 5.0/12.0}};
NekDouble HingedBody::AdamsMoulton_coeffs[3][3] = {
        { 1.0       ,  0.0      , 0.0     },
        { 1.0/2.0   ,  1.0/2.0  , 0.0     },
        { 5.0/12.0  ,  2.0/3.0  ,-1.0/12.0}};


/**
 * @class HingedBody
 */
HingedBody::HingedBody(
        const LibUtilities::SessionReaderSharedPtr          &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields)
    : FSIBody(pSession, pFields)
{
}

/**
 *
 */
void HingedBody::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields,
        const std::map<std::string, std::string>            &pParams)
{
    FSIBody::v_InitObject(pFields, pParams);

    std::stringstream      sStream;
    std::string            sString;

    // Parameter list for creating FilterAeroForces
    std::map<std::string, std::string> vParams;

    // Get timestep
    m_session->LoadParameter("TimeStep", m_timestep);

    // OutputFrequency
    auto it = pParams.find("StartTime");
    if (it == pParams.end())
    {
        m_startTime = 0.0;
    }
    else
    {
        LibUtilities::Equation equ(m_session, it->second);
        m_startTime = equ.Evaluate();
    }

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = 1;
    }
    else
    {
        vParams[it->first] = it->second;
        LibUtilities::Equation equ1(m_session, it->second);
        m_outputFrequency = round(equ1.Evaluate());
    }

    // OutputFile
    it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        m_doOutput = false;
        vParams["OutputFrequency"] = "0";
    }
    else
    {
        m_doOutput = true;
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        vParams["OutputFile"] = it->second;
        m_outputFile = it->second;
    }
    if (!(m_outputFile.length() >= 4
          && m_outputFile.substr(m_outputFile.length() - 4) == ".mot"))
    {
        m_outputFile += ".mot";
    }

    // Moment of Inertia
    it = pParams.find("I");
    ASSERTL0(it != pParams.end(), "Missing parameter 'I'.");
    LibUtilities::Equation equ2(m_session, it->second);
    m_I = equ2.Evaluate();

    // Torsional spring coefficient
    it = pParams.find("K");
    ASSERTL0(it != pParams.end(), "Missing parameter 'K'.");
    LibUtilities::Equation equ3(m_session, it->second);
    m_K = equ3.Evaluate();

    // Torsional damping coefficient
    it = pParams.find("C");
    ASSERTL0(it != pParams.end(), "Missing parameter 'C'.");
    LibUtilities::Equation equ4(m_session, it->second);
    m_C = equ4.Evaluate();

    // Hinge point (default is origin)
    m_hingePoint = Array<OneD, NekDouble> (3, 0.0);
    it = pParams.find("HingePoint");
    if ( it != pParams.end() )
    {
        ASSERTL0(!(it->second.empty()), "Missing parameter 'HingePoint'.");
        // Add to vParams so that FilterAeroForces already return
        //    moment around this point
        vParams["MomentPoint"] = it->second;

        sStream.str(it->second);
        sStream.clear();
        for (int j = 0; j < 3; ++j)
        {
            sStream >> sString;
            if (!sString.empty())
            {
                LibUtilities::Equation equ6(m_session, sString);
                m_hingePoint[j] = equ6.Evaluate();
            }
        }
    }

    // Axis direction (only makes sense for 3D; default is z-axis)
    m_axis = Array<OneD, NekDouble> (3, 0.0);
    m_axis[2] = 1.0;
    it = pParams.find("HingeAxis");
    if ( it != pParams.end() )
    {
        ASSERTL0(!(it->second.empty()), "Missing parameter 'HingeAxis'.");

        sStream.str(it->second);
        sStream.clear();
        // Normalisation factor
        NekDouble norm = 0.0;
        for (int j = 0; j < 3; ++j)
        {
            sStream >> sString;
            if (!sString.empty())
            {
                LibUtilities::Equation equ7(m_session, sString);
                m_axis[j] = equ7.Evaluate();
                norm += m_axis[j]*m_axis[j];
            }
        }
        // Normalise direction
        for( int j = 0; j < 3; ++j)
        {
            m_axis[j] /= sqrt(norm);
        }
    }

    // Create FilterAeroForces object
    it = pParams.find("Boundary");
    ASSERTL0(it != pParams.end(), "Missing parameter 'Boundary'.");
    vParams[it->first] = it->second;
    m_filterForces = MemoryManager<SolverUtils::FilterAeroForces>::
                                AllocateSharedPtr(m_session, vParams);
    m_filterForces->Initialise(pFields, 0.0);

    // Determine time integration order
    std::string intMethod = m_session->GetSolverInfo("TIMEINTEGRATIONMETHOD");
    if(boost::iequals(intMethod, "IMEXOrder1"))
    {
        m_intSteps = 1;
    }
    else if (boost::iequals(intMethod, "IMEXOrder2"))
    {
        m_intSteps = 2;
    }
    else if (boost::iequals(intMethod, "IMEXOrder3"))
    {
        m_intSteps = 3;
    }
    else
    {
        ASSERTL0(false, "Time integration method not supported.");
    }

    // Initialise variables
    m_angle    = 0.0;
    m_velocity = Array<OneD,NekDouble> (m_intSteps, 0.0);
    m_moment   = Array<OneD,NekDouble> (m_intSteps, 0.0);

    // Get initial displacement and velocity
    GetInitialCondition(pFields);

    // Open outputstream and write header
    if( m_doOutput)
    {
        m_index = 0;
        if ( pFields[0]->GetComm()->GetRank() == 0)
        {
            // Open output stream
            bool adaptive;
            m_session->MatchSolverInfo("Driver", "Adaptive", adaptive, false);
            if (adaptive)
            {
                m_outputStream.open(m_outputFile.c_str(), ofstream::app);
            }
            else
            {
                m_outputStream.open(m_outputFile.c_str());
            }
            // Write header
            m_outputStream << "# Angular velocity and angle of hinged bodies"
                           << endl;
            m_outputStream << "#" << " Hinge Position = " << " (";
            m_outputStream.width(8);
            m_outputStream << setprecision(4) << m_hingePoint[0];
            m_outputStream.width(8);
            m_outputStream << setprecision(4) << m_hingePoint[1];
            m_outputStream.width(8);
            m_outputStream << setprecision(4) << m_hingePoint[2];
            m_outputStream << ")" << endl;
            m_outputStream << "#" << " Rotation  Axis = " << " (";
            m_outputStream.width(8);
            m_outputStream << setprecision(4) << m_axis[0];
            m_outputStream.width(8);
            m_outputStream << setprecision(4) << m_axis[1];
            m_outputStream.width(8);
            m_outputStream << setprecision(4) << m_axis[2];
            m_outputStream << ")" << endl;
            m_outputStream << "# Boundary regions: "
                           << m_bondaryString.c_str() << endl;
            m_outputStream << "#";
            m_outputStream.width(7);
            m_outputStream << "Time";
            m_outputStream.width(15);
            m_outputStream <<  "Angular_Vel";
            m_outputStream.width(15);
            m_outputStream <<  "Angle";
            m_outputStream << endl;
        }
    }
}

void HingedBody::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pDisplFields,
        const NekDouble                                      &time)
{
    // Do not move before m_startTime
    if(time < m_startTime)
    {
        return;
    }

    // Determine order for this time step
    static int nCalls = 0;
    ++nCalls;
    int order = min(nCalls, m_intSteps);

    int expdim = pFields[0]->GetGraph()->GetMeshDimension();

    // Get aerodynamic forces
    Array<OneD, NekDouble> moments(expdim, 0.0);
    m_filterForces->GetTotalMoments(pFields, moments, time);

    // Shift moment storage
    for(int n = m_intSteps-1; n > 0; --n)
    {
        m_moment[n] = m_moment[n-1];
    }

    // Calculate total moment
    if( expdim == 2)
    {
        // Only have moment in z-direction
        m_moment[0] = moments[0];
    }
    else
    {
        // Project to axis direction
        m_moment[0] = 0;
        for(int j = 0; j < 3; ++j)
        {
            m_moment[0] += moments[j] * m_axis[j];
        }
    }

    // Account for torsional spring and damping contributions
    m_moment[0] -= m_K * m_angle - m_C * m_velocity[0];

    // Shift velocity storage
    for(int n = m_intSteps-1; n > 0; --n)
    {
        m_velocity[n] = m_velocity[n-1];
    }

    // Update velocity
    for(int j = 0; j < order; ++j)
    {
        m_velocity[0] += m_timestep *
            AdamsBashforth_coeffs[order-1][j] * m_moment[j] / m_I;
    }
    // Update position
    for(int j = 0; j < order; ++j)
    {
        m_angle += m_timestep *
            AdamsMoulton_coeffs[order-1][j] * m_velocity[j];
    }

    // One minus cosine m_angle
    NekDouble oneMinusCosD = 1.0 - cos(m_angle);
    // Sine m_angle
    NekDouble sinD         = sin(m_angle);

    // Create matrix representing displacement from rotation: rotMat = (R - I)
    Array< OneD, Array< OneD, NekDouble>> rotMat(3);
    for(int j = 0; j < 3; ++j)
    {
        rotMat[j] = Array< OneD, NekDouble> (3, 0.0);
        // Diagonal terms
        rotMat[j][j] = oneMinusCosD * (m_axis[j]*m_axis[j] - 1.0);
    }
    // Off diagonal terms
    rotMat[0][1] = m_axis[0]*m_axis[1] * oneMinusCosD - m_axis[2] * sinD;
    rotMat[1][0] = m_axis[1]*m_axis[0] * oneMinusCosD + m_axis[2] * sinD;
    rotMat[0][2] = m_axis[0]*m_axis[2] * oneMinusCosD + m_axis[1] * sinD;
    rotMat[2][0] = m_axis[2]*m_axis[0] * oneMinusCosD - m_axis[1] * sinD;
    rotMat[1][2] = m_axis[1]*m_axis[2] * oneMinusCosD - m_axis[0] * sinD;
    rotMat[2][1] = m_axis[2]*m_axis[1] * oneMinusCosD + m_axis[0] * sinD;

    // Loop coordinates
    for( int i = 0; i < pDisplFields.num_elements(); ++i)
    {
        // Get boundary expansions
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr> bndConds =
                pDisplFields[i]->GetBndConditions();
        Array<OneD, MultiRegions::ExpListSharedPtr> bndExp =
                pDisplFields[i]->GetBndCondExpansions();

        // Loop on boundary regions
        for( int n = 0; n < bndConds.num_elements(); ++n)
        {
            // Only modify boundary regions corresponding to this body
            if(m_boundaryRegionIsInList[n] == 1)
            {
                // Number of points on this boundary
                int nPts = bndExp[n]->GetTotPoints();

                // Get position relative to hinge point (on mesh coordinates)
                Array<OneD, Array<OneD, NekDouble> >  coords(3);
                for( int j = 0; j < 3; ++j)
                {
                    coords[j] = Array<OneD, NekDouble> (nPts, 0.0);
                }
                bndExp[n]->GetCoords(coords[0], coords[1], coords[2]);
                for( int j = 0; j < 3; ++j)
                {
                    Vmath::Sadd (nPts, -1.0*m_hingePoint[j],
                                        coords[j], 1, coords[j], 1);
                }

                // Calculate displacement in this direction
                Vmath::Zero(nPts, bndExp[n]->UpdatePhys(), 1);
                for( int j = 0; j < 3; ++j)
                {
                    Vmath::Svtvp(nPts, rotMat[i][j], coords[j], 1,
                        bndExp[n]->UpdatePhys(), 1, bndExp[n]->UpdatePhys(), 1);
                }

                // Update coefficients
                bndExp[n]->FwdTrans_BndConstrained(bndExp[n]->GetPhys(),
                                            bndExp[n]->UpdateCoeffs());
                if (pDisplFields[i]->GetExpType() == MultiRegions::e3DH1D)
                {
                    bndExp[n]->HomogeneousFwdTrans(bndExp[n]->GetCoeffs(),
                                            bndExp[n]->UpdateCoeffs());
                }
            }
        }
    }

    if( m_doOutput)
    {
        m_filterForces->Update(pFields, time);
        ++m_index;
        if ( !(m_index % m_outputFrequency) )
        {
            if ( pFields[0]->GetComm()->GetRank() == 0)
            {
                m_outputStream.width(8);
                m_outputStream << setprecision(6) << time;
                m_outputStream.width(15);
                m_outputStream << setprecision(8) << m_velocity[0];
                m_outputStream.width(15);
                m_outputStream << setprecision(8) << m_angle;
                m_outputStream << endl;
            }
        }
    }
}

void HingedBody::GetInitialCondition(
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pFields)
{
    // Check if this partition contains the boundary, and get a bnd id
    int bndId = -1;
    for( int n = 0; n < pFields[0]->GetBndConditions().num_elements(); ++n)
    {
        if(m_boundaryRegionIsInList[n] == 1)
        {
            bndId = n;
        }
    }
    // Pick one process to broadcast the result
    LibUtilities::CommSharedPtr comm = pFields[0]->GetComm()->GetRowComm();
    int bcastRank = -1;
    if (bndId != -1)
    {
        bcastRank = comm->GetRank();
    }
    comm->AllReduce(bcastRank, LibUtilities::ReduceMax);
    ASSERTL0(bcastRank >= 0, "Boundary not found.");

    // Get mapping, which will have the initial coordinates and velocities
    GlobalMapping::MappingSharedPtr mapping =
        GlobalMapping::Mapping::Load(m_session, pFields);

    // Process bcastRank calculates the initial conditions
    if( comm->GetRank() == bcastRank)
    {
        // Allocate storage
        int totPts = pFields[0]->GetTotPoints();
        int bndPts = pFields[0]->GetBndCondExpansions()[bndId]->GetTotPoints();
        Array<OneD, Array<OneD,NekDouble>> tmp1(3);
        Array<OneD, Array<OneD,NekDouble>> tmp2(3);
        Array<OneD, Array<OneD,NekDouble>> tmp3(3);
        Array<OneD, Array<OneD,NekDouble>> rBnd(3);
        Array<OneD, Array<OneD,NekDouble>> rCartBnd(3);
        Array<OneD, Array<OneD,NekDouble>> velBnd(3);
        for (int i = 0; i < 3; ++i)
        {
            tmp1[i]     = Array<OneD,NekDouble> (totPts, 0.0);
            tmp2[i]     = Array<OneD,NekDouble> (totPts, 0.0);
            tmp3[i]     = Array<OneD,NekDouble> (totPts, 0.0);
            rBnd[i]     = Array<OneD,NekDouble> (bndPts, 0.0);
            rCartBnd[i] = Array<OneD,NekDouble> (bndPts, 0.0);
            velBnd[i]   = Array<OneD,NekDouble> (bndPts, 0.0);
        }

        // tmp1 = Position relative to hinge in Cartesian coordinates
        mapping->GetCartesianCoordinates(tmp1[0], tmp1[1], tmp1[2]);
        // tmp2 = Position relative to hinge in mesh coordinates
        pFields[0]->GetCoords(tmp2[0],tmp2[1],tmp2[2]);
        for (int i = 0; i < 3; ++i)
        {
            Vmath::Sadd(totPts, -1.0*m_hingePoint[i], tmp1[i], 1, tmp1[i], 1);
            Vmath::Sadd(totPts, -1.0*m_hingePoint[i], tmp2[i], 1, tmp2[i], 1);
        }
        // tmp3 = coordinate velocities
        mapping->GetCoordVelocity(tmp3);

        // Extract values at the boundary
        for (int i = 0; i < 3; ++i)
        {
            pFields[0]->ExtractPhysToBnd(bndId, tmp1[i], rCartBnd[i]);
            pFields[0]->ExtractPhysToBnd(bndId, tmp2[i], rBnd[i]);
            pFields[0]->ExtractPhysToBnd(bndId, tmp3[i], velBnd[i]);
        }
        // Pick a point not coincident with the hinge
        int pt;
        for (pt = 0; pt < bndPts; ++pt)
        {
            if( abs(rCartBnd[0][pt]) > 1e-6 || abs(rCartBnd[1][pt]) > 1e-6)
            {
                break;
            }
        }
        ASSERTL0(pt != bndPts, "Could not find valid point");

        // Calculate angle and velocity. TO DO: Extend to 3D
        NekDouble thetaCart = atan2(rCartBnd[1][pt], rCartBnd[0][pt]);
        NekDouble theta     = atan2(rBnd[1][pt]    , rBnd[0][pt]);
        m_angle             = thetaCart - theta;

        // The magnitude of the angular velocity is abs(V) / R
        NekDouble radius = sqrt( rCartBnd[0][pt]*rCartBnd[0][pt] +
                                 rCartBnd[1][pt]*rCartBnd[1][pt]);
        m_velocity[0]    = (velBnd[0][pt]*velBnd[0][pt] +
                            velBnd[1][pt]*velBnd[1][pt]) / radius;
        // Now check the sign
        //    (check both conditions in case either x or y is zero)
        if( velBnd[1][pt]*rCartBnd[0][pt] < 0 ||
            velBnd[0][pt]*rCartBnd[1][pt] > 0)
        {
            m_velocity[0] = -m_velocity[0];
        }
    }
    // Broadcast the initial conditions
    comm->Bcast(m_angle, bcastRank);
    comm->Bcast(m_velocity[0] , bcastRank);
}

}
