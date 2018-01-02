///////////////////////////////////////////////////////////////////////////////
//
// File: RigidBody.cpp
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

#include <IncNavierStokesSolver/FSI/RigidBody.h>
#include <LibUtilities/BasicUtils/Equation.h>

namespace Nektar
{

std::string RigidBody::className =
    GetFSIBodyFactory().RegisterCreatorFunction("Rigid",
    RigidBody::create, "Rigid body motion");

NekDouble RigidBody::AdamsBashforth_coeffs[3][3] = {
        { 1.0       , 0.0       , 0.0     },
        { 3.0/2.0   ,-1.0/2.0   , 0.0     },
        { 23.0/12.0 ,-4.0/3.0   , 5.0/12.0}};
NekDouble RigidBody::AdamsMoulton_coeffs[3][3] = {
        { 1.0       ,  0.0      , 0.0     },
        { 1.0/2.0   ,  1.0/2.0  , 0.0     },
        { 5.0/12.0  ,  2.0/3.0  ,-1.0/12.0}};


/**
 * @class RigidBody
 */
RigidBody::RigidBody(
        const LibUtilities::SessionReaderSharedPtr          &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields)
    : FSIBody(pSession, pFields)
{
}

/**
 *
 */
void RigidBody::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields,
        const std::map<std::string, std::string>            &pParams)
{
    std::stringstream      sStream;
    std::string            sString;

    // Parameter list for creating FilterAeroForces
    std::map<std::string, std::string> vParams;

    // Get timestep
    m_session->LoadParameter("TimeStep", m_timestep);

    // OutputFrequency
    auto it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = 1;
    }
    else
    {
        vParams[it->first] = it->second;
        LibUtilities::Equation equ(m_session, it->second);
        m_outputFrequency = round(equ.Evaluate());
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

    // Number of degrees of freedom
    it = pParams.find("TranslationDOFs");
    ASSERTL0(it != pParams.end(), "Missing parameter 'TranslationDOFs'.");
    LibUtilities::Equation equ(m_session, it->second);
    m_nDof = round(equ.Evaluate());

    // Mass
    it = pParams.find("M");
    ASSERTL0(it != pParams.end(), "Missing parameter 'M'.");
    equ = LibUtilities::Equation (m_session, it->second);
    m_M = equ.Evaluate();

    // Spring coefficient
    m_K = Array<OneD, NekDouble> (m_nDof);
    it = pParams.find("K");
    ASSERTL0(it != pParams.end(), "Missing parameter 'K'.");
    sStream.str(it->second);
    sStream.clear();
    for (int i = 0; i < m_nDof; ++i)
    {
        sStream >> sString;
        if (!sString.empty())
        {
            equ = LibUtilities::Equation (m_session, it->second);
            m_K[i] = equ.Evaluate();
        }
        else
        {
            ASSERTL0(i != 0, "Empty parameter 'K'.");
            // One value provided: assume K is the same in all directions
            m_K[i] = m_K[0];
        }
    }

    // Damping coefficient
    m_C = Array<OneD, NekDouble> (m_nDof);
    it = pParams.find("C");
    ASSERTL0(it != pParams.end(), "Missing parameter 'C'.");
    sStream.str(it->second);
    sStream.clear();
    for (int i = 0; i < m_nDof; ++i)
    {
        sStream >> sString;
        if (!sString.empty())
        {
            equ = LibUtilities::Equation (m_session, it->second);
            m_C[i] = equ.Evaluate();
        }
        else
        {
            ASSERTL0(i != 0, "Empty parameter 'C'.");
            // One value provided: assume C is the same in all directions
            m_C[i] = m_C[0];
        }
    }

    // Allocate m_directions
    m_directions = Array<OneD, Array<OneD, NekDouble> > (3);
    //Initialise directions to default values (ex, ey, ez)
    for (int i = 0; i < 3; ++i)
    {
        m_directions[i]    = Array<OneD, NekDouble>(3, 0.0);
        m_directions[i][i] = 1.0;
    }
    //Override with input from xml file if defined
    for (int i = 0; i < 3; ++i)
    {
        std::stringstream tmp;
        tmp << i+1;
        std::string dir = "Direction" + tmp.str();
        it = pParams.find(dir);
        if ( it != pParams.end() )
        {
            ASSERTL0(!(it->second.empty()), "Missing parameter '"+dir+"'.");
            // Add to vParams so that FilterAeroForces already return
            //    forces on these directions
            vParams[it->first] = it->second;

            sStream.str(it->second);
            sStream.clear();
            // Normalisation factor
            NekDouble norm = 0.0;
            for (int j = 0; j < 3; ++j)
            {
                sStream >> sString;
                if (!sString.empty())
                {
                    equ = LibUtilities::Equation (m_session, it->second);
                    m_directions[i][j] = equ.Evaluate();
                    norm += m_directions[i][j]*m_directions[i][j];
                }
            }
            // Normalise direction
            for( int j = 0; j < 3; ++j)
            {
                m_directions[i][j] /= sqrt(norm);
            }
        }
    }

    // Create FilterAeroForces object
    it = pParams.find("Boundary");
    ASSERTL0(it != pParams.end(), "Missing parameter 'Boundary'.");
    vParams[it->first] = it->second;
    m_filterForces = MemoryManager<SolverUtils::FilterAeroForces>::
                                AllocateSharedPtr(m_session, vParams);

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

    // Allocate variables
    m_displacement = Array<OneD,NekDouble> (m_nDof, 0.0);
    m_velocity     = Array<OneD, Array<OneD,NekDouble>> (m_intSteps);
    m_force        = Array<OneD, Array<OneD,NekDouble>> (m_intSteps);
    for (int i = 0; i < m_intSteps; ++i)
    {
        m_velocity[i] = Array<OneD, NekDouble> (m_nDof, 0.0);
        m_force[i]    = Array<OneD, NekDouble> (m_nDof, 0.0);
    }

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
            m_outputStream << "# Displacement of bodies" << endl;
            for( int i = 0; i < m_nDof; i++ )
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
            m_outputStream << "# Boundary regions: "
                           << m_bondaryString.c_str() << endl;
            m_outputStream << "#";
            m_outputStream.width(7);
            m_outputStream << "Time";
            for( int i = 1; i <= m_nDof; ++i )
            {
                m_outputStream.width(14);
                m_outputStream <<  "Displ_" << i;
            }
        }
    }
}

void RigidBody::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr>    &pDisplFields,
        const NekDouble                                      &time)
{
    // Determine order for this time step
    static int nCalls = 0;
    ++nCalls;
    int order = min(nCalls, m_intSteps);

    int expdim = pFields[0]->GetGraph()->GetMeshDimension();

    // Get aerodynamic forces
    Array<OneD, NekDouble> aeroForces(expdim, 0.0);
    m_filterForces->GetTotalForces(pFields, aeroForces, time);

    // Rotate force storage
    Array<OneD, NekDouble> tmp = m_force[m_intSteps-1];
    for(int n = m_intSteps-1; n > 0; --n)
    {
        m_force[n] = m_force[n-1];
    }
    m_force[0] = tmp;

    // Calculate total force
    for( int i = 0; i < m_nDof; ++i)
    {
        m_force[0][i] = aeroForces[i] -
                m_K[i] * m_displacement[i] - m_C[i] * m_velocity[0][i];
    }

    // Rotate velocity storage, keeping value of velocity[0]
    Vmath::Vcopy(m_nDof, m_velocity[0], 1, m_velocity[m_intSteps-1], 1);
    tmp = m_velocity[m_intSteps-1];
    for(int n = m_intSteps-1; n > 0; --n)
    {
        m_velocity[n] = m_force[n-1];
    }
    m_velocity[0] = tmp;

    for( int i = 0; i < m_nDof; ++i)
    {
        // Update velocity
        for(int j = 0; j < order; ++j)
        {
            m_velocity[0][i] += m_timestep *
                AdamsBashforth_coeffs[order-1][j] * m_force[j][i] / m_M;
        }
        // Update position
        for(int j = 0; j < order; ++j)
        {
            m_displacement[i] += m_timestep *
                AdamsMoulton_coeffs[order-1][j] * m_velocity[j][i];
        }
    }

    // Loop coordinates
    for( int i = 0; i < pDisplFields.num_elements(); ++i)
    {
        // Get boundary expansions
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr> bndConds =
                pDisplFields[i]->GetBndConditions();
        Array<OneD, MultiRegions::ExpListSharedPtr> bndExp =
                pDisplFields[i]->GetBndCondExpansions();

        // Calculate displacement in this direction
        NekDouble displ = 0;
        for (int j = 0; j < m_nDof; ++j)
        {
            displ += m_displacement[j] * m_directions[j][i];
        }
        // Loop on boundary regions
        for( int n = 0; n < bndConds.num_elements(); ++n)
        {
            // Only modify boundary regions corresponding to this body
            if(m_boundaryRegionIsInList[n] == 1)
            {
                // Number of points on this boundary
                int nPts = bndExp[n]->GetTotPoints();

                // Fill boundary conditions
                Vmath::Fill(nPts, displ, bndExp[n]->UpdatePhys(), 1);

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
        if ( m_index % m_outputFrequency )
        {
            if ( pFields[0]->GetComm()->GetRank() == 0)
            {
                m_outputStream.width(8);
                m_outputStream << setprecision(6) << time;
                for( int i = 1; i <= expdim; i++ )
                {
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << m_displacement[i];
                }
            }
        }
    }
}

}
