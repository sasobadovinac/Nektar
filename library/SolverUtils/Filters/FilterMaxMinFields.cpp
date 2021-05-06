///////////////////////////////////////////////////////////////////////////////
//
// File FilterMaxMinFields.cpp
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

#include <SolverUtils/Filters/FilterMaxMinFields.h>
#include <CompressibleFlowSolver/Misc/VariableConverter.h>

namespace Nektar
{
namespace SolverUtils
{
std::string FilterMaxMinFields::className =
        GetFilterFactory().RegisterCreatorFunction(
                "MaxMinFields", FilterMaxMinFields::create);

FilterMaxMinFields::FilterMaxMinFields(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem>      &pEquation,
    const ParamMap &pParams)
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

    // Load flag for max or min
    it = pParams.find("MaxOrMin");
    std::string sOption = it->second.c_str();
    if ( boost::iequals(sOption, "maximun") ||
         boost::iequals(sOption, "max") )
    {
        m_isMax = true;
    }
    else if ( boost::iequals(sOption, "minimum") ||
              boost::iequals(sOption, "min") )
    {
        m_isMax = false;
    }
    else
    {
        ASSERTL1(false, "MaxOrMin needs to be max or min.");
    }


    // Initialize other member vars
    m_problemType = eCompressible;
}

FilterMaxMinFields::~FilterMaxMinFields()
{
}


void FilterMaxMinFields::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Initialise output arrays
    FilterFieldConvert::v_Initialise(pFields, time);

    // Allocate storage
    m_curFieldsPhys.resize(m_variables.size());
    m_outFieldsPhys.resize(m_variables.size());
    for (int n = 0; n < m_variables.size(); ++n)
    {
        m_curFieldsPhys[n] = Array<OneD, NekDouble>(
                                 pFields[0]->GetTotPoints(), 0.0);
        m_outFieldsPhys[n] = Array<OneD, NekDouble>(
                                 pFields[0]->GetTotPoints(), 0.0);
    }

    // Check type of problem
    std::string firstVarName = pFields[0]->GetSession()->GetVariable(0);
    if (boost::iequals(firstVarName, "u"))
    {
        m_problemType = eIncompressible;
    }
    else if (boost::iequals(firstVarName, "rho"))
    {
        m_problemType = eCompressible;
    }
    else
    {
        m_problemType = eOthers;
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

void FilterMaxMinFields::v_ProcessSample(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
          std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
    const NekDouble &time)
{
    boost::ignore_unused(fieldcoeffs, time);

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

    // Forward transform and put into m_outFields
    for (int n = 0; n < nVars; ++n)
    {
        pFields[0]->FwdTrans_IterPerExp(m_outFieldsPhys[n], m_outFields[n]);
    } 

}

void FilterMaxMinFields::v_PrepareOutput(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(pFields, time);

    m_fieldMetaData["NumberOfFieldDumps"] =
        boost::lexical_cast<std::string>(m_numSamples);
}

NekDouble FilterMaxMinFields::v_GetScale()
{
    return 1.0;
}

}
}
