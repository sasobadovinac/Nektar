///////////////////////////////////////////////////////////////////////////////
//
// File FilterBodyFittedVelocity.cpp
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
// Description: Maximun/minimum solution fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Filters/FilterBodyFittedVelocity.h>
#include <CompressibleFlowSolver/Misc/VariableConverter.h>

#include "FieldUtils/ProcessModules/ProcessBodyFittedVelocity.h"

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ContField.h>

using std::cout;
using std::endl;

namespace Nektar
{
namespace SolverUtils
{
std::string FilterBodyFittedVelocity::className =
        GetFilterFactory().RegisterCreatorFunction(
                "BodyFittedVelocity", FilterBodyFittedVelocity::create);

FilterBodyFittedVelocity::FilterBodyFittedVelocity(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem>        &pEquation,
    const ParamMap                             &pParams)
    : FilterFieldConvert(pSession, pEquation, pParams)
{
    // Load sampling frequency
    auto it = pParams.find("SampleFrequency");
    if (it == pParams.end())
    {
        m_sampleFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        m_sampleFrequency = round(equ.Evaluate());
    }

    // Load flag for filter type
    it = pParams.find("OriginalOrMaxOrMin");
    if (it == pParams.end())
    {
        m_filterType = eOriginal;
    }
    else
    {
        std::string sOptionType = it->second.c_str();
        if (boost::iequals(sOptionType, "maximum") ||
            boost::iequals(sOptionType, "max"))
        {
            m_filterType = eMax;
        }
        else if (boost::iequals(sOptionType, "minumun") ||
                 boost::iequals(sOptionType, "min"))
        {
            m_filterType = eMin;
        }
        else if (boost::iequals(sOptionType, "original"))
        {
            m_filterType = eOriginal;
        }
        else
        {
            WARNINGL0(false, "Detailed filter type is not found, use original");
            m_filterType = eOriginal;
        }
    }

    // Load bodyFittedCoordinate from a file generated by FieldConvert -m bodyFittedVelocity 
    it = pParams.find("BodyFittedCoordinateFile");
    ASSERTL0(it->second.length() > 0, "Missing parameter 'BodyFittedCoordinateFile'.");

    if ( it->second.find_last_of('.') != std::string::npos)
    {

         m_bodyFittedCooriateFile = it->second;
    }
    else
    {
        std::stringstream outname;
        outname << it->second << ".fld";
        m_bodyFittedCooriateFile = outname.str();
    }

    m_initialized = m_restartFile != "";

    std::cout << "m_filterType = " << m_filterType << std::endl;
    std::cout << "m_bodyFittedCooriateFile = " << m_bodyFittedCooriateFile << std::endl;
}

FilterBodyFittedVelocity::~FilterBodyFittedVelocity()
{
}


void FilterBodyFittedVelocity::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Initialise output arrays
    FilterFieldConvert::v_Initialise(pFields, time);


    // Generate the body-fitted coordinate system and distance field
    // bfcsDir[i][j][k]
    //   i - dir vec:   0-main tangential, 1-normal, (2-minor tangential)
    //   j - component: 0-x, 1-y, (2-z)
    //   k - point id
    const int npoints = pFields[0]->GetTotPoints();
    //m_distance = Array<OneD, NekDouble>(npoints, 9999); // should be deleted later
    m_bfcsDir  = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(m_spaceDim);
    for (int i=0; i<m_spaceDim; ++i)
    {
        m_bfcsDir[i] = Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
        for (int j=0; j<m_spaceDim; ++j)
        {
            m_bfcsDir[i][j] = Array<OneD, NekDouble>(npoints);
        }
    }


    // Initialize the body fitted coordinate in the file
    // Ref: FilterFieldConvert.cpp (v_Initialise) & ProcessInterpField.cpp (Process)
    m_bfsVars.resize(m_spaceDim * m_spaceDim); // 4 or 9
    if (m_spaceDim == 2)
    {
        m_bfsVars[0] = "bfc_dir0_x";
        m_bfsVars[1] = "bfc_dir0_y";
        m_bfsVars[2] = "bfc_dir1_x";
        m_bfsVars[3] = "bfc_dir1_y";
    }
    else
    {
        m_bfsVars[0] = "bfc_dir0_x";
        m_bfsVars[1] = "bfc_dir0_y";
        m_bfsVars[2] = "bfc_dir0_z";
        m_bfsVars[3] = "bfc_dir1_x";
        m_bfsVars[4] = "bfc_dir1_y";
        m_bfsVars[5] = "bfc_dir1_z";
        m_bfsVars[6] = "bfc_dir2_x";
        m_bfsVars[7] = "bfc_dir2_y";
        m_bfsVars[8] = "bfc_dir2_z";
    }


    std::cout << "Loading the bfc file." << std::endl;
    // Load file
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> fieldDef;
    std::vector<std::vector<NekDouble> > fieldData;
    LibUtilities::FieldMetaDataMap fieldMetaData;
    LibUtilities::FieldIOSharedPtr fld =
    LibUtilities::FieldIO::CreateForFile(m_session, m_bodyFittedCooriateFile);
    fld->Import(m_bodyFittedCooriateFile, fieldDef, fieldData, fieldMetaData);

    // Get a tmp array to hold the coeffs from the nput file
    Array<OneD, NekDouble> tmp_coeffs;
    tmp_coeffs = Array<OneD, NekDouble>(pFields[0]->GetNcoeffs(), 0.0);

    // Extract data and convert to physical space
    for (int k=0; k<m_bfsVars.size(); ++k)
    {
        // In the fieldDef, fieldDef[0] is the quad and fieldDef[1] is the tri region
        // Each of fieldDef[i] contains all the fields we expected in an fld
        for (int i = 0; i < fieldData.size(); ++i)
        {
            //std::cout << "  i = " << i << std::endl;
            pFields[0]->ExtractDataToCoeffs(
                fieldDef[i],
                fieldData[i],
                m_bfsVars[k],
                tmp_coeffs);
        }

        // Bwd transorm the data and put them into m_bfcsDir[][]
        int ii, jj;
        jj = k % m_spaceDim;
        ii = k / m_spaceDim;

        std::cout << "  ii = " << ii << ", jj = " << jj << std::endl;

        pFields[0]->BwdTrans(tmp_coeffs, m_bfcsDir[ii][jj]);


    }




    /*
    // We change the flag for expansion to DisContField so that the copy
    // constructor for DisContField will be called in AppendExpList() when
    // CreateFields() is called. In the constructor of DisContField, the
    // m_bndCondExpansions will be generated and it is needed in
    // GenPntwiseBodyFittedCoordSys().
    // After the body-fitted coordinate system is generated, we set the flags
    // back so that the following modules will not be influenced. 
    // *Note: The flag m_declareExpansionAsContField is set true in the base
    // class (ProcessBoundaryExtract) constructor of ProcessBodyFittedVelocity
    // class. The ClearField will not change this flag, and it will lead to an
    // error in dynamic_pointer_cast in m_f->AppendExpList in CreateFields() in
    // FilterFieldConvert::OutputField. Thus, we explicitly set it false.
    bool isExpContField    = m_f->m_declareExpansionAsContField;    // save initial status
    bool isExpDisContField = m_f->m_declareExpansionAsDisContField;

    m_f->m_declareExpansionAsDisContField = true;
    FilterFieldConvert::CreateFields(pFields); // Generate m_bndCondExpansions
    
    ProcessBodyFittedVelocity bfv(m_f);
    //bfv.GenPntwiseBodyFittedCoordSys(m_wallRefId, m_assistVec, m_distance, m_bfcsDir,
    //    m_isAngleCheck, m_distTol, 1.0e-12, 1.0e-4, 1.0e-12);

    m_f->m_declareExpansionAsContField    = isExpContField;    // was set true in bfv
    m_f->m_declareExpansionAsDisContField = isExpDisContField; // was set true aobve

    m_f->ClearField();                         // clear the m_f for output useage
    */
    
    // Update distanceToWall, uvw_bfc in m_outFields
    // (The restart file does not include them so the values are useless)
    // With rst file    - default variables are available, other bfs and dist are not correct
    // Without rst file - everything is zero in m_outFields
    // Solution - discard the data in restart file, use current data 


    // Allocate storage
    m_curFieldsVels_Car.resize(m_spaceDim);
    m_curFieldsVels.resize(m_spaceDim); // m_spaceDim = m_nAddFields-1
    m_outFieldsVels.resize(m_spaceDim);

    for (int n = 0; n < m_spaceDim; ++n)
    {
        m_curFieldsVels_Car[n] = Array<OneD, NekDouble>(npoints, 0.0);
        m_curFieldsVels[n]     = Array<OneD, NekDouble>(npoints, 0.0);
        m_outFieldsVels[n]     = Array<OneD, NekDouble>(npoints, 0.0);

        // Re-initialize the coeff
        //pFields[0]->FwdTransLocalElmt(m_outFieldsVels[n], m_outFields[m_nVars+1+n]);
    }

    // m_outFields contains the initialization values
    if (m_initialized)
    {
        for (int n = 0; n < m_spaceDim; ++n)
        {
            pFields[0]->BwdTrans(m_outFields[m_nVars+1+n], m_outFieldsVels[n]);
            if (pFields[0]->GetWaveSpace())
            {
                pFields[0]->HomogeneousBwdTrans(m_outFieldsVels[n],
                                                m_outFieldsVels[n]);
            }
        }
    }

}



void FilterBodyFittedVelocity::v_FillVariablesName(
    const Array<OneD, const MultiRegions::ExpListSharedPtr>& pFields)
{
    // Fill in the var names using the routine in base class
    FilterFieldConvert::v_FillVariablesName(pFields);
    

    // Check type of problem
    std::vector<std::string> varNames = pFields[0]->GetSession()->GetVariables();
    if (boost::iequals(varNames[0], "u"))
    {
        m_problemType = eIncompressible;
        m_spaceDim    = pFields.size() - 1;
    }
    else if (boost::iequals(varNames[0], "rho") &&
             boost::iequals(varNames[1], "rhou"))
    {
        m_problemType = eCompressible;
        m_spaceDim    = pFields.size() - 2;
    }
    else
    {
        WARNINGL0(false, "Problem type is not claear. Please check");
        m_problemType = eOthers;
        m_spaceDim    = 2;
    }

    m_nVars      = m_variables.size(); // Include extra vars
    m_nFields    = pFields.size();     // Conservative vars for cfs
    m_nAddFields = 3;

    // Fill in the body-fitted velocity names
    // bfc is for body-fitted coordinate
    m_variables.push_back("distanceToWall");
    m_variables.push_back("u_bfc");
    m_variables.push_back("v_bfc");
    if (m_spaceDim == 3)
    {
        m_variables.push_back("w_bfc");
        ++m_nAddFields;
    }
    else if (m_spaceDim != 2)
    {
        ASSERTL0(false, "Unsupported dimension");
    }
    
}

// For compressible flows, also compute max/min for extra variabless, which can
// include (u, v, w, p, T. s, a, Mach, sensor, ArtificialVisc), but not sure
// for variables (divVelSquare, curlVelSquare, Ducros) when Ducros sensor is
// turned on. Now presume Ducros sensor is turned off.
// For other cases, including incompressible flows, only compute max/min for
// variables in pFields.
// 
// curFieldsCoeffs is the coeffs     at current time step, size=[*][m_ncoeffs]
// m_curFieldsPhys is the phys value at current time step, size=[*][m_npoints]
// m_outFields     is the coeffs     for output fld
// m_outFieldsPhys is the phys value for output fld
// PS. curFieldsCoeffs is not set to be a member variable since its size will
//     be increased by CompressibleFlowSystem::v_ExtraFldOutput to also include
//     coeffs for extra variables. fieldcoeffs is not used.

void FilterBodyFittedVelocity::v_ProcessSample(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
          std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
    const NekDouble &time)
{
    boost::ignore_unused(time);


    // Use shift to get the current u,v,w in Cartesian coordinate
    // cfs_2D: rho,rhou,rhov,E,u,v,p,T,s,a,Mach,Sensor,distanceToWall,u_bfc,v_bfc
    // inc_2D: u,v,p
    int shift = (m_problemType == eCompressible) ? (m_spaceDim + 2) : 0;
    
    // Get current u,v,w in Cartesian coordinate
    for (int n = 0; n < m_spaceDim; ++n)
    {
        pFields[0]->BwdTrans(fieldcoeffs[shift+n], m_curFieldsVels_Car[n]);
        if (pFields[0]->GetWaveSpace())
        {
            pFields[0]->HomogeneousBwdTrans(m_curFieldsVels_Car[n],
                                            m_curFieldsVels_Car[n]);
        }
    }


    /*if (m_problemType == eCompressible)
    {
        std::vector<Array<OneD, NekDouble> > curFieldsCoeffs(m_nFields);
        std::vector<std::string> tmpVariables;    // dummy vector 
        for (int n = 0; n < m_nFields; ++n)
        {
            curFieldsCoeffs[n] = pFields[n]->GetCoeffs();
        }

        // Get extra variables, then curFieldsCoeffs.size() == m_nVars
        auto equ = m_equ.lock();
        ASSERTL0(equ, "Weak pointer expired");
        equ->ExtraFldOutput(curFieldsCoeffs, tmpVariables);

        // curFieldsCoeffs and m_outFields has different size
        // curFieldsCoeffs has conservative and the default variables
        // m_outFields has (m_nVars+m_nAddFields) variables

        // Updata m_outFieldsPhys and m_curFieldsPhys (vel in Cartesian system) 
        for (int n = 0; n < m_spaceDim; ++n)
        {
            pFields[0]->BwdTrans(curFieldsCoeffs[m_spaceDim+2+n],
                                 m_curFieldsVels_Car[n]); // cur uvw_Car
            pFields[0]->BwdTrans(m_outFields[m_nVars+1+n],
                                 m_outFieldsVels[n]);     // out uvw_bfc
        }

        // Fill coeffs for other fields
        for (int n=0; n<m_nVars; ++n)
        {
            m_outFields[n] = curFieldsCoeffs[n];
        }
    }
    else
    {
        // Updata m_curFieldsPhys and m_outFieldsPhys
        for (int n = 0; n < m_spaceDim; ++n)
        {
            m_curFieldsVels_Car[n] = pFields[n]->GetPhys();

            pFields[0]->BwdTrans(m_outFields[n], m_outFieldsVels[n]);
            if (pFields[0]->GetWaveSpace())
            {
                pFields[0]->HomogeneousBwdTrans(m_outFieldsVels[n], m_outFieldsVels[n]);
            }          
        }

        // Fill coeffs for other fields
        for (int n = 0; n < m_nVars; ++n)
        {
            m_outFields[n] = pFields[n]->GetCoeffs();
        }
    }*/

    // Project the velocity into the body-fitted coordinate system
    int npoints = pFields[0]->GetTotPoints();
    Array<OneD, NekDouble> vel_tmp(npoints);
    
    for (int i=0; i<m_spaceDim; ++i)     // loop for bfc velocity
    {
        Vmath::Zero(npoints, m_curFieldsVels[i], 1);

        for (int j=0; j<m_spaceDim; ++j)
        {
            Vmath::Vmul(npoints, m_curFieldsVels_Car[j], 1, m_bfcsDir[i][j], 1, vel_tmp,  1);
            Vmath::Vadd(npoints, vel_tmp, 1,  m_curFieldsVels[i],  1, m_curFieldsVels[i], 1);
        }
    }


    // Get max/min/original for the u/v/w_bfc field
    for(int n = 0; n < m_spaceDim; ++n)
    {
        size_t length = m_outFieldsVels[n].size();
        
        if (m_filterType == eMax)
        {
            // Compute max
            for (int i = 0; i < length; ++i)
            {
                if (m_curFieldsVels[n][i] > m_outFieldsVels[n][i])
                {
                    m_outFieldsVels[n][i] = m_curFieldsVels[n][i];
                }
            }
        }
        else if (m_filterType == eMin)
        {
            // Compute min
            for (int i = 0; i < length; ++i)
            {
                if (m_curFieldsVels[n][i] < m_outFieldsVels[n][i])
                {
                    m_outFieldsVels[n][i] = m_curFieldsVels[n][i];
                }
            }
        }
        else
        {
            // Original field
            m_outFieldsVels[n] = m_curFieldsVels[n];
        }
    }
    
    // Forward transform and put into m_outFields
    //pFields[0]->FwdTransLocalElmt(m_distance, m_outFields[m_nVars]);
    for (int n = 0; n < m_spaceDim; ++n)
    {
        pFields[0]->FwdTransLocalElmt(m_outFieldsVels[n], m_outFields[m_nVars+1+n]);
    }

}

void FilterBodyFittedVelocity::v_PrepareOutput(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(pFields, time);

    m_fieldMetaData["NumberOfFieldDumps"] =
        boost::lexical_cast<std::string>(m_numSamples);
}

NekDouble FilterBodyFittedVelocity::v_GetScale()
{
    return 1.0;
}

}
}
