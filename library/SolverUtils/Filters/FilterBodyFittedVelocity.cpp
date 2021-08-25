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

    // Load wall reference id
    it = pParams.find("WallReferenceId");
    if (it == pParams.end())
    {
        WARNINGL0(false, "Wall ref id not found, use 0 by default");
        m_wallRefId = 0;
    }
    else
    {
        LibUtilities::Equation equId(m_session->GetInterpreter(), it->second);
        m_wallRefId = round(equId.Evaluate());
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

    // Get the assistant vector
    m_assistVec = Array<OneD, NekDouble>(3, 0.0);
    it = pParams.find("AssistDirection");
    if (it == pParams.end())
    {
        m_assistVec[0] = 0.0;
        m_assistVec[1] = 0.0;
        m_assistVec[2] = 1.0;
    }
    else
    {
        m_assistDir.str(it->second);
        m_assistDir >> m_assistVec[0] >> m_assistVec[1] >> m_assistVec[2];

        if (m_assistDir.fail())
        {
            std::string ss;
            ss = "The input assistant direction [" + m_assistDir.str()
               + "] has incorrect dimensions, or contains illegal comma or "
               + "quotation marks. Use [0,0,1] by default.";
            WARNINGL0(false, ss);

            m_assistVec[0] = 0.0;
            m_assistVec[1] = 0.0;
            m_assistVec[2] = 1.0;
        }
    }

    // Check if should use angle check
    it = pParams.find("AngleCheck");
    if (it == pParams.end())
    {
        m_isAngleCheck = false;
    }
    else
    {
        std::string sOptionChk = it->second.c_str();
        m_isAngleCheck = (boost::iequals(sOptionChk, "true")) ||
                         (boost::iequals(sOptionChk, "yes"));
    }

    // Get the distance tolerence
    it = pParams.find("DistTol");
    if (it == pParams.end())
    {
        m_distTol = 1.0e-12;
    }
    else
    {
        LibUtilities::Equation equDistTol(m_session->GetInterpreter(), it->second);
        m_distTol = equDistTol.Evaluate();
    }


    std::cout << "m_wallRefId    = " << m_wallRefId <<std::endl;
    std::cout << ", m_filterType = " << m_filterType << std::endl;
    std::cout << "assistVec      = " << m_assistVec[0] << ", " << m_assistVec[1] <<", "<< m_assistVec[2] << std::endl;
    std::cout << "angleChk       = " << m_isAngleCheck << std::endl;
    std::cout << "distTol        = " << m_distTol << std::endl;
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
    const int npoints = pFields[0]->GetTotPoints();
    std::cout << "npts 2 = "<< npoints  << std::endl;

    m_distance = Array<OneD, NekDouble>(npoints, 9999);
    m_bfcsDir  = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(m_spaceDim);
    for (int i=0; i<m_spaceDim; ++i)
    {
        m_bfcsDir[i] = Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
        for (int j=0; j<m_spaceDim; ++j)
        {
            m_bfcsDir[i][j] = Array<OneD, NekDouble>(npoints);
        }
    }

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

    std::cout << "Begin " << std::endl;
    std::cout << "Before: isCon = " << m_f->m_declareExpansionAsContField 
              << ", isDisCon = "    << m_f->m_declareExpansionAsDisContField << std::endl;

    
    bool isExpContField    = m_f->m_declareExpansionAsContField;    // save initial status
    bool isExpDisContField = m_f->m_declareExpansionAsDisContField;

    m_f->m_declareExpansionAsDisContField = true;
    FilterFieldConvert::CreateFields(pFields); // Generate m_bndCondExpansions
    
    ProcessBodyFittedVelocity bfv(m_f);
    bfv.GenPntwiseBodyFittedCoordSys(m_wallRefId, m_assistVec, m_distance, m_bfcsDir,
        m_isAngleCheck, m_distTol, 1.0e-12, 1.0e-4, 1.0e-12);

    m_f->m_declareExpansionAsContField    = isExpContField;    // was set true in bfv
    m_f->m_declareExpansionAsDisContField = isExpDisContField; // was set true aobve

    m_f->ClearField();                         // clear the m_f for output useage

    std::cout << "After: isCon = "      << m_f->m_declareExpansionAsContField 
              << ", isDisCon = " << m_f->m_declareExpansionAsDisContField << std::endl;
    std::cout << "End " << std::endl;

    
    //===============================================
    // Updated above
    //===============================================

    


    // Update distanceToWall, uvw_bfc in m_outFields
    // (The restart file does not include them so the values are useless)
    // With rst file    - default variables are available, other bfs and dist are not correct
    // Without rst file - everything is zero in m_outFields
    // Solution - discard the data in restart file, use current data 

    cout << "spaceDim   = " << m_spaceDim   << endl;
    cout << "nAddFields = " << m_nAddFields << endl;

    // Allocate storage
    NekDouble initVal;
    if (m_filterType == eMin)
    {
        initVal = 9999.0;
    }
    else if (m_filterType == eMax)
    {
        initVal = -9999.0;
    }
    else
    {
        initVal = 0.0;
    }

    m_curFieldsPhys.resize(m_nAddFields - 1);
    m_outFieldsPhys.resize(m_nAddFields - 1);

    for (int n = 0; n < (m_nAddFields-1); ++n)
    {
        m_curFieldsPhys[n] = Array<OneD, NekDouble>(
                                 pFields[0]->GetTotPoints(), 0.0);
        m_outFieldsPhys[n] = Array<OneD, NekDouble>(
                                 pFields[0]->GetTotPoints(), initVal);
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

    m_nVars      = m_variables.size();
    m_nFields    = pFields.size();
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
    boost::ignore_unused(fieldcoeffs, time, pFields);
    
    /*
    int nFields, nVars;
 
    if (m_problemType==eCompressible)
    {
        nFields = pFields.size();       // Conservative vars
        nVars   = m_variables.size();   // Include extra vars

        std::vector<Array<OneD, NekDouble> > curFieldsCoeffs(nFields);
        std::vector<std::string> tmpVariables;    // dummy vector 
        for (int n = 0; n < nFields; ++n)
        {
            curFieldsCoeffs[n] = pFields[n]->GetCoeffs();
        }

        // Get extra variables, then curFieldsCoeffs.size() == nVars
        auto equ = m_equ.lock();
        ASSERTL0(equ, "Weak pointer expired");
        equ->ExtraFldOutput(curFieldsCoeffs, tmpVariables);

        // curFieldsCoeffs and m_outFields has different size
        // curFieldsCoeffs has conservative and the default variables
        // m_outFields has nVar variables, include distanceToWall, ut,vn,w

        // Updata m_curFieldsPhys and m_outFieldsPhys
        for (int n = 0; n < nVars; ++n)
        {
            pFields[0]->BwdTrans(curFieldsCoeffs[n], m_curFieldsPhys[n]);
            pFields[0]->BwdTrans(m_outFields[n],     m_outFieldsPhys[n]);
        }
    }
    else
    {
        nFields = pFields.size();
        nVars   = nFields;

        // Updata m_curFieldsPhys and m_outFieldsPhys
        for (int n = 0; n < nVars; ++n)
        {
            m_curFieldsPhys[n] = pFields[n]->GetPhys();

            pFields[0]->BwdTrans(m_outFields[n], m_outFieldsPhys[n]);
            if (pFields[0]->GetWaveSpace())
            {
                pFields[0]->HomogeneousBwdTrans(m_outFieldsPhys[n], m_outFieldsPhys[n]);
            }          
        }
    }


    // Get max/min for each field
    for(int n = 0; n < nVars; ++n)
    {
        size_t length = m_outFieldsPhys[n].size();
        
        if (m_isMax)
        {
            // Compute max
            for (int i = 0; i < length; ++i)
            {
                if (m_curFieldsPhys[n][i] > m_outFieldsPhys[n][i])
                {
                    m_outFieldsPhys[n][i] = m_curFieldsPhys[n][i];
                }
            }
        }
        else
        {
            // Compute min
            for (int i = 0; i < length; ++i)
            {
                if (m_curFieldsPhys[n][i] < m_outFieldsPhys[n][i])
                {
                    m_outFieldsPhys[n][i] = m_curFieldsPhys[n][i];
                }
            }
        }
    }
    */
    
    // Forward transform and put into m_outFields
    pFields[0]->FwdTrans_IterPerExp(m_distance, m_outFields[m_nVars]);
    for (int n = 0; n < (m_nAddFields-1); ++n)
    {
        pFields[0]->FwdTrans_IterPerExp(m_outFieldsPhys[n], m_outFields[m_nVars+1+n]);
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
