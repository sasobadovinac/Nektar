///////////////////////////////////////////////////////////////////////////////
//
// File FilterAverageField.cpp
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
// Description: Average solution fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Filters/FilterCleanReynoldsStresses.h>
namespace Nektar
{
namespace SolverUtils
{
std::string FilterCleanReynoldsStresses::className =
        GetFilterFactory().RegisterCreatorFunction(
                "CleanReynoldsStresses", FilterCleanReynoldsStresses::create);

FilterCleanReynoldsStresses::FilterCleanReynoldsStresses(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const ParamMap &pParams)
    : FilterFieldConvert(pSession, pParams)
{
    // Averaged file
    auto it = pParams.find("AvgFile");
    ASSERTL0(it != pParams.end(), "Missing parameter 'AvgFile'.");
    if ( it->second.find_last_of('.') != string::npos)
    {   
        m_avgFile = it->second;
    }
    else
    {   
        std::stringstream outname;
        outname << it->second << ".fld";
        m_avgFile = outname.str();
    }
    // Restart file
    it = pParams.find("RestartFile");
    if (it == pParams.end())
    {   
        m_restartFile = "";
    }
    else
    {   
        ASSERTL0(it->second.length() > 0, "Missing parameter 'RestartFile'.");
        if ( it->second.find_last_of('.') != string::npos)
        {   
            m_restartFile = it->second;
        }
        else
        {   
            std::stringstream outname;
            outname << it->second << ".fld";
            m_restartFile = outname.str();
        }
    }


}

FilterCleanReynoldsStresses::~FilterCleanReynoldsStresses()
{
}

void FilterCleanReynoldsStresses::v_Initialise(
		    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
		        const NekDouble &time)
{
    v_FillVariablesName(pFields);
    // m_variables need to be filled by a derived class
    m_avgFields.resize(origFields);
    m_avgFields_phys.resize(origFields);
    m_outFields.resize(m_variables.size());
    m_outFields_phys.resize(m_variables.size());
    m_delta.resize(m_variables.size());
    int nfield;
    for (int n = 0; n < m_variables.size(); n++)
    {
        nfield              = (n < pFields.num_elements())? n: 0;
        m_outFields[n]      = Array<OneD, NekDouble>(pFields[nfield]->GetNcoeffs(), 0.0);
        m_outFields_phys[n] = Array<OneD, NekDouble>(pFields[nfield]->GetTotPoints(), 0.0);
	m_delta[n]          = Array<OneD, NekDouble>(pFields[nfield]->GetTotPoints(), 0.0);
    }
    for (int n = 0; n < origFields; n++ )
    {
        //nfield = (n < pFields.num_elements())? n: 0;
	m_avgFields[n]       = Array<OneD, NekDouble>(pFields[n]->GetNcoeffs(), 0.0);
	m_avgFields_phys[n]  = Array<OneD, NekDouble>(pFields[n]->GetTotPoints(), 0.0);
    }

    vel    = Array<OneD, NekDouble>(pFields[0]->GetTotPoints(), 0.0);

    m_fieldMetaData["InitialTime"] = boost::lexical_cast<std::string>(time);
    // read averaged file
    if ( m_avgFile != "" )
    {
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> fieldDef_avg;
        std::vector<std::vector<NekDouble> > fieldData_avg;
        LibUtilities::FieldMetaDataMap fieldMetaData_avg;
        LibUtilities::FieldIOSharedPtr fld =
        LibUtilities::FieldIO::CreateForFile(m_session, m_avgFile);
        fld->Import(m_avgFile, fieldDef_avg, fieldData_avg, fieldMetaData_avg);
	for ( int j = 0; j < origFields; j++ )
	{
        	for ( int n = 0; n < fieldData_avg.size(); n++  )
		{
	     		pFields[j]->ExtractDataToCoeffs(
                	    		fieldDef_avg[n],
                    			fieldData_avg[n],
                    			m_variables[j],
                    			m_avgFields[j]); 
		}
		pFields[j]->BwdTrans(m_avgFields[j], m_avgFields_phys[j]);
                if (pFields[j]->GetWaveSpace())
           	{   
               	 	pFields[j]->HomogeneousBwdTrans(m_avgFields_phys[j], m_avgFields_phys[j]);
           	}
	}
    }
    int  k;
    nfield = -1;
    // read the restart file
    if ( m_restartFile != "" )
    {
	// Load file
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> fieldDef;
        std::vector<std::vector<NekDouble> > fieldData;
        LibUtilities::FieldMetaDataMap fieldMetaData;
        LibUtilities::FieldIOSharedPtr fld =
            LibUtilities::FieldIO::CreateForFile(m_session, m_restartFile);
        fld->Import(m_restartFile, fieldDef, fieldData, fieldMetaData);
        // Extract fields to output
        for (int j = 0; j < m_variables.size(); ++j)
        {
            for(k = 0; k < pFields.num_elements(); ++k)
            {
                if(pFields[k]->GetSession()->GetVariable(k) == m_variables[j])
                {
                    nfield = k;
                    break;
                }
            }
            if(nfield == -1)
            {
                nfield = 0;
            }

            for (int i = 0; i < fieldData.size(); ++i)
            {
                pFields[nfield]->ExtractDataToCoeffs(
                    			fieldDef[i],
                    			fieldData[i],
                    			m_variables[j],
                    			m_outFields[j]);
            }
	    pFields[nfield]->BwdTrans(m_outFields[j], m_outFields_phys[j]);
            if (pFields[nfield]->GetWaveSpace())
	    {
		pFields[nfield]->HomogeneousBwdTrans(m_outFields_phys[j], m_outFields_phys[j]);
	    }
        }
        if (fieldMetaData.count("NumberOfFieldDumps"))
        {
            m_numSamples = atoi(fieldMetaData["NumberOfFieldDumps"].c_str());
        }
        else
        {
            m_numSamples = 1;
        }

        if(fieldMetaData.count("InitialTime"))
        {
            m_fieldMetaData["InitialTime"] = fieldMetaData["InitialTime"];
        }
        NekDouble scale = v_GetScale();
        for (int n = 0; n < m_outFields_phys.size(); ++n)
        {
            Vmath::Smul(m_outFields_phys[n].num_elements(),
                        1.0/scale, m_outFields_phys[n],
                        1, m_outFields_phys[n], 1);
        }
    }
}

void FilterCleanReynoldsStresses::v_FillVariablesName(
		    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    dim        = pFields.num_elements() - 1;
    origFields = pFields.num_elements();
    for (int n = 0; n < origFields; n++)
    {
        m_variables.push_back(pFields[n]->GetSession()->GetVariable(n));
    }
    if (dim == 2)
    {
        m_variables.push_back("uu");
        m_variables.push_back("uv");
        m_variables.push_back("vv");
    }
    else if (dim == 3)
    {
        m_variables.push_back("uu");
        m_variables.push_back("uv");
        m_variables.push_back("uw");
        m_variables.push_back("vv");
        m_variables.push_back("vw");
        m_variables.push_back("ww");
    }
    else
    {
        ASSERTL0(false, "Unsupported dimension");
    }
}
void FilterCleanReynoldsStresses::v_ProcessSample(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    int i, j, n;
    int nq             = pFields[0]->GetTotPoints();
    int dim            = pFields.num_elements() - 1;
    bool waveSpace     = pFields[0]->GetWaveSpace();

    // Sampling velocity and their differences
    for(n = 0; n < origFields; n++)
    {
	// sum up the velocity
        // keep primitive variable fixed
	Vmath::Vadd(nq,m_avgFields_phys[n],1,m_outFields_phys[n],1,m_outFields_phys[n], 1);
	// calculate differences between time-dependent vels and averaged vels
        if (waveSpace)
        {
            pFields[n]->HomogeneousBwdTrans(pFields[n]->GetPhys(), vel);
	    Vmath::Svtsvtp(nq, 1.0, vel, 1, -1.0,  m_avgFields_phys[n], 1, m_delta[n],1);
        }
        else
        {
	    Vmath::Svtsvtp(nq, 1.0, pFields[n]->GetPhys(), 1, -1.0,  m_avgFields_phys[n], 1, m_delta[n],1);
        }
    }
    // Sampling Reynolds stresses
    n = origFields;
    for (i = 0; i < dim; i++)
    {
        for (j = i; j < dim; j++)
        {
            Vmath::Vmul(nq, m_delta[i], 1, m_delta[j], 1, m_delta[n], 1);
	    n++;
        }
    }
    // Saving Reynolds stresses
    for ( i = origFields; i < m_variables.size(); i++ )
    {
	    Vmath::Vadd(nq,
			m_delta[i], 1,
			m_outFields_phys[i], 1,
			m_outFields_phys[i], 1);
    }
}

void FilterCleanReynoldsStresses::v_PrepareOutput(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    m_fieldMetaData["NumberOfFieldDumps"] =
        boost::lexical_cast<std::string>(m_numSamples);
    // Set wavespace to false, as calculations were performed in physical space
    bool waveSpace = pFields[0]->GetWaveSpace();
    pFields[0]->SetWaveSpace(false);
    // Forward transform and put into m_outFields (except pressure)
    for (int i = 0; i < m_outFields.size(); ++i)
    {
        if (i != dim)
        {
            pFields[0]->FwdTrans_IterPerExp(m_outFields_phys[i], m_outFields[i]);
        }
    }
    // Restore waveSpace
    pFields[0]->SetWaveSpace(waveSpace);
}

NekDouble FilterCleanReynoldsStresses::v_GetScale()
{
    return 1.0 / m_numSamples;
}

}
}
