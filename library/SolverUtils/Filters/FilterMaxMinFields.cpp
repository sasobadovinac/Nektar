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
    it = pParams.find("MaxMin");
    if (it == pParams.end())
    {
        m_isMax = true;
    }
    else
    {
        std::string sOption = it->second.c_str();
        m_isMax = !(boost::iequals(sOption, "minimum") || 
                    boost::iequals(sOption, "min"));
    }

    // Initialize other member vars
    m_problemType = eCompressible;
    m_spaceDim    = 2;
}

FilterMaxMinFields::~FilterMaxMinFields()
{
}


void FilterMaxMinFields::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{

    // Initialise output arrays
    // FilterMaxMinFields::v_FillVariablesName is called inside
    FilterFieldConvert::v_Initialise(pFields, time);

    // Allocate storage
    m_curFieldsPhys.resize(m_variables.size());
    m_outFieldsPhys.resize(m_variables.size());

    for (int n = 0; n < m_variables.size(); ++n)
    {
        m_curFieldsPhys[n] = Array<OneD, NekDouble>(pFields[0]->GetTotPoints(), 0.0);
        m_outFieldsPhys[n] = Array<OneD, NekDouble>(pFields[0]->GetTotPoints(), 0.0);
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

    // Check the space dim
    int nFields = pFields.size();
    if (m_problemType == eCompressible)
    {
        m_spaceDim = nFields==5 ? 3 : 2;
    }
    else if (m_problemType == eIncompressible)
    {
        m_spaceDim = nFields==4 ? 3 : 2; // will not be used
    }
    else
    {
        m_spaceDim = 2;                 // Will not be used
    }
        

}

/*
// Add fields for filtered variables
void FilterMaxMinFields::v_FillVariablesName(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    // General variable names 
    FilterFieldConvert::v_FillVariablesName(pFields);

    // Check type of problem
    std::string firstVarName = pFields[0]->GetSession()->GetVariable(0);
    if (boost::iequals(firstVarName, "u"))
    {
        m_problemType = eIncompressible;
    }
    else if (boost::iequals(firstVarName, "rho"))
    {
        m_problemType = eCompressible;
        // Remove sensor from [rho,rhou,rhov,E,u,v,p,T,s,a,Mach,Sensor]
        m_variables.pop_back();
    }
    else
    {
        m_problemType = eOthers;
    }

    // Check the space dim
    int nFields = pFields.size();
    if (m_problemType == eCompressible)
    {
        m_spaceDim = nFields==5 ? 3 : 2;
    }
    else if (m_problemType == eIncompressible)
    {
        m_spaceDim = nFields==4 ? 3 : 2; // will not be used
    }
    else
    {
        m_spaceDim = 2;                 // Will not be used
    }
    
    // Change the variable names
    
    for (int i=0; i<m_variables.size(); ++i)
    {
        m_variables[i] = m_variables[i] + v_GetFileSuffix();
        std::cout << m_variables[i] << std::endl;
    }
    


}
*/





void FilterMaxMinFields::v_ProcessSample(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
          std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
    const NekDouble &time)
{
    boost::ignore_unused(fieldcoeffs, time);

    // Ref 
    // FilterReynoldsStresses::v_FillVariablesName in solvers/IncNavierStokesSolver/Filters/FilterReynoldsStresses.cpp
    // FilterReynoldsStresses::v_Initialise


    //     fieldcoeffs is current coeffs at current time step, size=[n][m_ncoeffs]
    // m_curFieldsPhys is physical value at current time step, size=[n][m_npoints]
    // m_outFields is the coeffs for output fld, ref. FilterFieldConvert.cpp -> line 234,242
    // m_outFieldsPhys is the physical value for output fld

    // Generate array for physical variables
    const int nFields = pFields.size();       // field vars, rho, rhou, rhov, rhow, E
    const int nVars   = m_variables.size();   // extra vars, u, v, w, p, T. s, Mach, sensor

    const int nPhys   = pFields[0]->GetNpoints();
    const int nCoeffs = pFields[0]->GetNcoeffs();

    /*
    std::cout << "Size of filtered fields    = " << nFields << std::endl;
    std::cout << "Size of filtered variables = " << nVars   << std::endl;


    std::vector<Array<OneD, NekDouble> > fieldPhys(nFields);                                              
    for (int n = 0; n < nFields; ++n)
    {
        fieldPhys[n] = pFields[n]->GetPhys();
    }
    
    std::vector<Array<OneD, NekDouble> > outFieldsPhys(nFields);
    for (int n = 0; n < nFields; ++n)
    {
        outFieldsPhys[n] = Array<OneD, NekDouble>(pFields[0]->GetTotPoints(), 0.0);
    }
    */
    


    // Need to create a dummy coeffs vector to get extra variables names...
    std::vector<Array<OneD, NekDouble> > tmpCoeffs(nFields);
    std::string tmpStr;
    for (int n = 0; n < nFields; ++n)
    {
        tmpCoeffs[n] = pFields[n]->GetCoeffs();
    }
    // Get extra variables
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer expired");
    equ->ExtraFldOutput(tmpCoeffs, tmpStr);



    /*
    // tmp for field vars in type of Array of Array
    Array<OneD, Array<OneD, NekDouble> > tmp(nFields);
    for (int i=0; i < pFields.size(); ++i)
    {
        tmp[i] = pFields[i]->GetPhys();
    }


    Array<OneD, Array<OneD, NekDouble> > velocity(m_spaceDim);
    //Array<OneD, Array<OneD, NekDouble> > velFwd  (m_spaceDim);
    for (int i = 0; i < m_spaceDim; ++i)
    {
        velocity[i] = Array<OneD, NekDouble> (pFields[0]->GetTotPoints());
    }


    Array<OneD, NekDouble> pressure(nPhys), temperature(nPhys);
    Array<OneD, NekDouble> entropy(nPhys);
    Array<OneD, NekDouble> soundSpeed(nPhys), mach(nPhys);
    Array<OneD, NekDouble> sensor(nPhys), SensorKappa(nPhys);

    auto equ = m_equ.lock();
    VariableConverterSharedPtr varConv = equ->GetVarConv();
    //MemoryManager<VariableConverter>::AllocateSharedPtr(
    //                m_session, m_spacedim);
    varConv->GetVelocityVector(tmp, velocity);
    varConv->GetPressure      (tmp, pressure);
    varConv->GetTemperature   (tmp, temperature);
    varConv->GetEntropy       (tmp, entropy);
    varConv->GetSoundSpeed    (tmp, soundSpeed);
    varConv->GetMach          (tmp, soundSpeed, mach);
    */
    
    // Updata m_curFieldsPhys
    /*
    for (int n=0; n<nFields; ++n)
    {
        m_curFieldsPhys[n] = pFields[n]->GetPhys();
    }
    */





    // Update m_outFieldsPhys
    // What is the following if() for?
    // Ref. FilterReynoldsStresses.cpp
    for (int j = 0; j < nVars; ++j)
    {
        pFields[0]->BwdTrans(m_outFields[j], m_outFieldsPhys[j]);
        if (pFields[0]->GetWaveSpace())
        {
            pFields[0]->HomogeneousBwdTrans(m_outFieldsPhys[j], m_outFieldsPhys[j]);
        }
    }
    /*
    for (int j = 0; j < nFields; ++j)
    {
        pFields[0]->BwdTrans(m_outFields[j], outFieldsPhys[j]);
        if (pFields[0]->GetWaveSpace())
        {
            pFields[0]->HomogeneousBwdTrans(outFieldsPhys[j], outFieldsPhys[j]);
        }
    }
    */

    // Get max/min for each field
    // Save the result in outFieldsPhys
    for(int n = 0; n < nVars; ++n)
    {
        size_t length = m_outFieldsPhys[n].size();
        
        if (m_isMax)
        {
            // Compute max
            for (int i=0; i<length; ++i)
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
            for (int i=0; i<length; ++i)
            {
                if (m_curFieldsPhys[n][i] < m_outFieldsPhys[n][i])
                {
                    m_outFieldsPhys[n][i] = m_curFieldsPhys[n][i];
                }
            }
        }

    }
    /*
    for(int n = 0; n < nFields; ++n)
    {
        size_t length = outFieldsPhys[n].size();
        
        if (m_isMax)
        {
            // Compute max
            for (int i=0; i<length; ++i)
            {
                if (fieldPhys[n][i] > outFieldsPhys[n][i])
                {
                    outFieldsPhys[n][i] = fieldPhys[n][i];
                }
            }
        }
        else
        {
            // Compute min
            for (int i=0; i<length; ++i)
            {
                if (fieldPhys[n][i] < outFieldsPhys[n][i])
                {
                    outFieldsPhys[n][i] = fieldPhys[n][i];
                }
            }
        }

    }
    */

    // Forward transform and put into m_outFields (except pressure)
    for (int n = 0; n < nVars; ++n)
    {
        pFields[0]->FwdTrans_IterPerExp(m_outFieldsPhys[n], m_outFields[n]);
    }
    /*
    for (int n = 0; n < nFields; ++n)
    {
        pFields[0]->FwdTrans_IterPerExp(outFieldsPhys[n], m_outFields[n]);
    }
    */   

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
