///////////////////////////////////////////////////////////////////////////////
//
// File FilterReynoldsStresses.cpp
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
// Description: Append Reynolds stresses to the average fields
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/Filters/FilterReynoldsStresses.h>

namespace Nektar
{
namespace SolverUtils
{
std::string FilterReynoldsStresses::className =
    GetFilterFactory().RegisterCreatorFunction("ReynoldsStresses",
        FilterReynoldsStresses::create);

/**
 * @class FilterReynoldsStresses
 *
 * @brief Calculate Favre stresses and the average fields
 *
 * This class calculates the averaged conservative variables and the Favre averaged stresses of the form
 * \f$ \widetilde{u' v'} \f$ or \f$ sqrt(\widetilde{u' u'}) \f$ .
 *
 * This is achieved by calculating
 * \f$ \widetilde{u' v'} = \widetilde{u v}-\widetilde{u}\widetilde{v} \f$,
 * \f$ sqrt(\widetilde{u' u'}) = sqrt(\widetilde{u^2}-\widetilde{u}^2\f$,
 * in which \widetilde{u} and \widetilde{v} can be directly calculated by the 
 * averaged conservative variables. 
 * 
 * Take the 2D case as example. The sampling process calculates the 
 * \f$ \sum{\rho}, \sum{\rho u}, \sum{\rho v}, \sum{E}, \sum{\rho u^2}, 
 * \sum{\rho uv}, \sum{\rho v^2} \f$ and the number of sampling (in m_numSamples). 
 * 
 * The output variables would be \f$ \bar{\rho}, \bar{\rho u}, \bar{\rho v}, 
 * \bar{E}, sqrt(\widetilde{u' u'}), \widetilde{u' v'}, \widetilde{v' v'} \f$,
 * based on the variables above.
 * 
 * For restarted simulations sampling variables need to be recovered based on 
 * the output variables. For example, \f$ \sum{\rho u} = \widetilde{u} * 
 * \widetilde{\rho} * m_numSamples \f$ and \f$ \sum{\rho uv} = \sum{\rho} * 
 * (\widetilde{u' v'} + \sum{\rho u} * \sum{\rho v} / \sum{\rho} / \sum{\rho})
 * \f$
 * 
 * The FilterSampler base class then divides the result by n, leading
 * to the Reynolds stress.
 *
 * It is also possible to perform the averages using an exponential moving
 *  average, in which case either the moving average parameter \f$ \alpha \f$
 * or the time constant \f$ \tau \f$ must be prescribed.
 */
FilterReynoldsStresses::FilterReynoldsStresses(
    const LibUtilities::SessionReaderSharedPtr         &pSession,
    const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
    const std::map<std::string, std::string> &pParams)
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
}

FilterReynoldsStresses::~FilterReynoldsStresses()
{
}

void FilterReynoldsStresses::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    int dim          = pFields.size() - 2;
    int nExtraFields = dim == 2 ? 3 : 6;
    int origFields   = pFields.size();

    // Allocate storage
    m_fields.resize(origFields + nExtraFields);
    m_delta.resize(dim);

    for (int n = 0; n < m_fields.size(); ++n)
    {
        m_fields[n] = Array<OneD, NekDouble>(pFields[0]->GetTotPoints(), 0.0);
    }
    for (int n = 0; n < m_delta.size(); ++n)
    {
        m_delta[n] = Array<OneD, NekDouble>(pFields[0]->GetTotPoints(), 0.0);
    }

    // Initialise output arrays
    FilterFieldConvert::v_Initialise(pFields, time);

    // Update m_fields if using restart file
    if (m_numSamples)
    {
        for (int j = 0; j < m_fields.size(); ++j)
        {
            pFields[0]->BwdTrans(m_outFields[j], m_fields[j]);
            if (pFields[0]->GetWaveSpace())
            {
                pFields[0]->HomogeneousBwdTrans(m_fields[j], m_fields[j]);
            }
        }
        TransOutputVarsIntoSampleVars(origFields);
    }
}

void FilterReynoldsStresses::v_FillVariablesName(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    int dim          = pFields.size() - 2;
    int origFields   = pFields.size();

    // Fill name of variables
    for (int n = 0; n < origFields; ++n)
    {
        m_variables.push_back(pFields[n]->GetSession()->GetVariable(n));
    }
    int nstressTensor = 0;
    if (dim == 2)
    {
        m_variables.push_back("uu");
        m_variables.push_back("uv");
        m_variables.push_back("vv");

        nstressTensor = 3;
        m_stressTensor = Array<OneD, Array<OneD, int>> {nstressTensor};
        for (int i = 0; i < nstressTensor; ++i)
        {
            m_stressTensor[i] =  Array<OneD, int> {2, 0};
        }
        m_stressTensor[0][0] = 0;
        m_stressTensor[0][1] = 0;

        m_stressTensor[1][0] = 0;
        m_stressTensor[1][1] = 1;
        
        m_stressTensor[2][0] = 1;
        m_stressTensor[2][1] = 1;

    }
    else if (dim == 3)
    {
        m_variables.push_back("uu");
        m_variables.push_back("uv");
        m_variables.push_back("uw");
        m_variables.push_back("vv");
        m_variables.push_back("vw");
        m_variables.push_back("ww");

        nstressTensor = 6;
        m_stressTensor = Array<OneD, Array<OneD, int>> {nstressTensor};
        for (int i = 0; i < nstressTensor; ++i)
        {
            m_stressTensor[i] =  Array<OneD, int> {2, 0};
        }
        m_stressTensor[0][0] = 0;
        m_stressTensor[0][1] = 0;
        
        m_stressTensor[1][0] = 0;
        m_stressTensor[1][1] = 1;
        
        m_stressTensor[2][0] = 0;
        m_stressTensor[2][1] = 2;
        
        m_stressTensor[3][0] = 1;
        m_stressTensor[3][1] = 1;
        
        m_stressTensor[4][0] = 1;
        m_stressTensor[4][1] = 2;
        
        m_stressTensor[5][0] = 2;
        m_stressTensor[5][1] = 2;
    }
    else
    {
        ASSERTL0(false, "Unsupported dimension");
    }

    m_indexScaleStt = 0;
    m_indexScaleEnd = origFields;
}

void FilterReynoldsStresses::v_ProcessSample(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
          std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
    const NekDouble &time)
{
    boost::ignore_unused(fieldcoeffs, time);
    
    int nq             = pFields[0]->GetTotPoints();

    Array<OneD, NekDouble> tmp(nq);

    Array<OneD, NekDouble> invDens{nq};
    Vmath::Sdiv(nq, 1.0, pFields[0]->GetPhys(), 1, invDens, 1);

    //Calculate Favre averaged variabls: Divide all sums by the sum of densities.
    int origFields   = pFields.size();
    int totalFields = m_fields.size();
    for (int j = origFields; j < totalFields; ++j)
    {
        int nstress = j - origFields;
        int moment0 = m_stressTensor[nstress][0] + 1;
        int moment1 = m_stressTensor[nstress][1] + 1;

        Vmath::Vmul(nq, 
            pFields[moment0]->GetPhys(), 1, 
            pFields[moment1]->GetPhys(), 1, 
            tmp, 1);
        
        Vmath::Vvtvp(nq, tmp, 1, invDens, 1, m_fields[j], 1, m_fields[j], 1);

    }

    for (int j = 0; j < origFields; ++j)
    {
        Vmath::Vadd(nq, 
            pFields[j]->GetPhys(), 1, 
            m_fields[j], 1, 
            m_fields[j], 1);
    }
}

void FilterReynoldsStresses::v_PrepareOutput(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(time);

    m_fieldMetaData["NumberOfFieldDumps"] =
        boost::lexical_cast<std::string>(m_numSamples);

    // Set wavespace to false, as calculations were performed in physical space
    bool waveSpace = pFields[0]->GetWaveSpace();
    pFields[0]->SetWaveSpace(false);

    TransSampleVarsIntoOutputVars(pFields.size());

    // Forward transform and put into m_outFields (except pressure)
    for (int i = 0; i < m_fields.size(); ++i)
    {
        pFields[0]->FwdTrans_IterPerExp(m_fields[i], m_outFields[i]);
    }

    // Restore waveSpace
    pFields[0]->SetWaveSpace(waveSpace);

    TransOutputVarsIntoSampleVars(pFields.size());
}

NekDouble FilterReynoldsStresses::v_GetScale()
{
    
    return 1.0 / m_numSamples;
}

/**
 * Transform sampling variables into output variables.
 * Sampling variables are: the sum of conservative variables and the sum of 
 * the product two momentems divided by density.
 * 
 * The output variables are: the averaged conservative variables and the Reynolds
 * stresses calculated based on Favre averaging.
 * 
*/

void FilterReynoldsStresses::TransSampleVarsIntoOutputVars(
    int origFields)
{
    size_t totalFields = m_fields.size();
    size_t npnts = m_fields[totalFields - 1].size();
    size_t nVeloc = origFields - 2;

    Array<OneD, NekDouble> tmp{npnts};
    Vmath::Sdiv(npnts, 1.0, m_fields[0], 1, tmp, 1);

    Array<OneD, Array<OneD, NekDouble>> tmpVol{nVeloc};
    for (int i = 0; i < nVeloc; ++i)
    {
        tmpVol[i] = Array<OneD, NekDouble> {npnts};
        Vmath::Vmul(npnts,
                    tmp,
                    1,
                    m_fields[i + 1],
                    1,
                    tmpVol[i],
                    1);
    }

    //Calculate Favre averaged variabls: Divide all sums by the sum of densities.
    for (int j = origFields; j < totalFields; ++j)
    {
        int nstress = j - origFields;
        int moment0 = m_stressTensor[nstress][0];
        int moment1 = m_stressTensor[nstress][1];

        Vmath::Vmul(npnts,
                    tmp,
                    1,
                    m_fields[j],
                    1,
                    m_fields[j],
                    1);

        Vmath::Vvtvm(npnts,
                    tmpVol[moment0],
                    1,
                    tmpVol[moment1],
                    1,
                    m_fields[j],
                    1,
                    m_fields[j],
                    1);
        Vmath::Neg(npnts, m_fields[j], 1);
    }
}

void FilterReynoldsStresses::TransOutputVarsIntoSampleVars(
    int origFields)
{
    size_t totalFields = m_fields.size();
    size_t npnts = m_fields[totalFields - 1].size();
    size_t nVeloc = origFields - 2;

    Array<OneD, NekDouble> tmp{npnts};
    Vmath::Sdiv(npnts, 1.0, m_fields[0], 1, tmp, 1);

    Array<OneD, Array<OneD, NekDouble>> tmpVol{nVeloc};
    for (int i = 0; i < nVeloc; ++i)
    {
        tmpVol[i] = Array<OneD, NekDouble> {npnts};
        Vmath::Vmul(npnts,
                    tmp,
                    1,
                    m_fields[i + 1],
                    1,
                    tmpVol[i],
                    1);
    }

    //Calculate Favre averaged variabls: Divide all sums by the sum of densities.
    for (int j = origFields; j < totalFields; ++j)
    {
        int nstress = j - origFields;
        int moment0 = m_stressTensor[nstress][0];
        int moment1 = m_stressTensor[nstress][1];

        Vmath::Vvtvp(npnts,
                    tmpVol[moment0],
                    1,
                    tmpVol[moment1],
                    1,
                    m_fields[j],
                    1,
                    m_fields[j],
                    1);

        Vmath::Vmul(npnts,
                    m_fields[0],
                    1,
                    m_fields[j],
                    1,
                    m_fields[j],
                    1);
    }
}

}
}
