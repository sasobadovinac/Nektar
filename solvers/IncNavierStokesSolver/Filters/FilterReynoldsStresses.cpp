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

#include <IncNavierStokesSolver/Filters/FilterReynoldsStresses.h>

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
 * @brief Append Reynolds stresses to the average fields
 *
 * This class appends the average fields with the Reynolds stresses of the form
 * \f$ \overline{u' v'} \f$.
 *
 * For the default case, this is achieved by calculating
 * \f$ C_{n} = \Sigma_{i=1}^{n} (u_i - \bar{u}_n)(v_i - \bar{v}_n)\f$
 * using the recursive relation:
 *
 * \f[ C_{n} = C_{n-1} + \frac{n}{n-1} (u_n - \bar{u}_n)(v_n - \bar{v}_n) \f]
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

    // Check if should use moving average
    it = pParams.find("MovingAverage");
    if (it == pParams.end())
    {
        m_movAvg = false;
    }
    else
    {
        std::string sOption = it->second.c_str();
        m_movAvg = (boost::iequals(sOption, "true")) ||
                   (boost::iequals(sOption, "yes"));
    }

    if (!m_movAvg)
    {
        return;
    }

    // Load alpha parameter for moving average
    it = pParams.find("alpha");
    if (it == pParams.end())
    {
        it = pParams.find("tau");
        if (it == pParams.end())
        {
            ASSERTL0(false, "MovingAverage needs either alpha or tau.");
        }
        else
        {
            // Load time constant
            LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
            NekDouble tau = equ.Evaluate();
            // Load delta T between samples
            NekDouble dT;
            m_session->LoadParameter("TimeStep", dT);
            dT = dT * m_sampleFrequency;
            // Calculate alpha
            m_alpha = dT / (tau + dT);
        }
    }
    else
    {
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        m_alpha = equ.Evaluate();
        // Check if tau was also defined
        it = pParams.find("tau");
        if (it != pParams.end())
        {
            ASSERTL0(false,
                     "Cannot define both alpha and tau in MovingAverage.");
        }
    }
    // Check bounds of m_alpha
    ASSERTL0(m_alpha > 0 && m_alpha < 1, "Alpha out of bounds.");
    m_mu = pSession->GetParameter("Kinvis");
}

FilterReynoldsStresses::~FilterReynoldsStresses()
{
}

void FilterReynoldsStresses::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    MultiRegions::ExpansionType exptype = pFields[0]->GetExpType();
    if (exptype==MultiRegions::e3DH2D || exptype==MultiRegions::e3DH1D)
    {
        m_dim = 3;
    }
    else
    {
        m_dim = pFields[0]->GetGraph()->GetSpaceDimension();
    }
    int nFields      = pFields.size();
    int nExtraFields = nFields * (nFields + 1) / 2;
    nExtraFields += m_dim + 1;

    // Allocate storage
    m_fields.resize(nFields + nExtraFields);
    m_delta.resize(nFields);

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
    }
}

void FilterReynoldsStresses::v_FillVariablesName(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    int nFields      = pFields.size();

    // Fill name of variables
    for (int n = 0; n < nFields; ++n)
    {
        m_variables.push_back(pFields[n]->GetSession()->GetVariable(n));
    }
    for (int i=0; i < nFields; ++i)
    {
        for (int j=0; j<=i; ++j)
        {
            std::string name = m_variables[i] + m_variables[j];
            m_variables.push_back(name);
        }
    }
    for (int i=0; i < m_dim; ++i)
    {
        std::string name = "u_iu_i"+ m_variables[i];
        m_variables.push_back(name);
    }
    m_variables.push_back("dissipation");
}

void FilterReynoldsStresses::v_ProcessSample(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
          std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
    const NekDouble &time)
{
    int i, j, n;
    int nq             = pFields[0]->GetTotPoints();
    int nFields        = pFields.size();
    NekDouble nSamples = (NekDouble)m_numSamples;

    // For moving average, take first sample as initial vector
    NekDouble alpha = m_alpha;
    if (m_numSamples == 1)
    {
        alpha = 1.0;
    }

    // Define auxiliary constants for averages
    NekDouble facOld, facAvg, facStress, facDelta;
    NekDouble fac1Skew, fac2Skew;
    if (m_movAvg)
    {
        facOld    = 1.0 - alpha;
        facAvg    = alpha;
        facStress = alpha;
        facDelta  = 1.0;
        fac1Skew  = 0.;
        fac2Skew  = alpha;
    }
    else
    {
        facOld    = 1.0;
        facAvg    = 1.0;
        facStress = nSamples / (nSamples - 1);
        facDelta  = 1.0 / nSamples;
        fac1Skew  = -1.0 / (nSamples - 1);
        fac2Skew  = facStress *  (nSamples + 1.) / (nSamples - 1.);
    }

    Array<OneD, NekDouble> vel(nq);
    // Update original velocities in phys space and calculate (\bar{u} - u_n)
    for (n = 0; n < nFields; ++n)
    {
        if (pFields[n]->GetWaveSpace())
        {
            pFields[n]->HomogeneousBwdTrans(pFields[n]->GetPhys(), vel);
        }
        else
        {
            vel = pFields[n]->GetPhys();
        }
        Vmath::Svtsvtp(
            nq, facAvg, vel, 1, facOld, m_fields[n], 1, m_fields[n], 1);
        Vmath::Svtvm(nq, facDelta, m_fields[n], 1, vel, 1, m_delta[n], 1);
    }

    // Ignore Reynolds stress for first sample (its contribution is zero)
    if (m_numSamples == 1)
    {
        return;
    }

    // Calculate C_{n} = facOld * C_{n-1} + facStress * deltaI * deltaJ
    Array<OneD, NekDouble> tmp = vel;
    std::vector<std::vector<int>> cov(nFields); //vu, wu, wv
    for (i = 0, n = nFields; i < nFields; ++i)
    {
        cov[i].resize(nFields);
        for (j = 0; j <= i; ++j, ++n)
        {
            Vmath::Vmul(nq, m_delta[i], 1, m_delta[j], 1, tmp, 1);
            Vmath::Svtsvtp(
                nq, facStress, tmp, 1, facOld, m_fields[n], 1, m_fields[n], 1);
            cov[i][j] = cov[j][i] = n;
        }
    }
    // Calculate u_i u_i u_j, j=0,1,...,m_dim
    for (i=0; i < m_dim; ++i, ++n)
    {
        for (j=0; j < m_dim; ++j)
        {
            Vmath::Vmul(nq, m_delta[i], 1, m_delta[i], 1, tmp, 1);
            Vmath::Vmul(nq, m_delta[j], 1, tmp, 1, tmp, 1);
            Vmath::Svtsvtp(
                nq, fac2Skew, tmp, 1, facOld, m_fields[n], 1, m_fields[n], 1);
            if (fac1Skew==0.)
            {
                continue;
            }
            Vmath::Vmul(nq, m_fields[cov[i][j]], 1, m_delta[i], 1, tmp, 1);
            Vmath::Svtsvtp(
                nq, 2.*fac1Skew, tmp, 1, facOld, m_fields[n], 1, m_fields[n], 1);
            Vmath::Vmul(nq, m_fields[cov[i][i]], 1, m_delta[j], 1, tmp, 1);
            Vmath::Svtsvtp(
                nq, fac1Skew, tmp, 1, facOld, m_fields[n], 1, m_fields[n], 1);
        }

    }
    // Calculate dissipation term
    Array<OneD, Array<OneD, NekDouble>> dvel(3);
    for (i=0; i < m_dim; ++i)
    {
        dvel[i] = Array<OneD, NekDouble>(nq);
    }
    Vmath::Zero(nq, tmp, 1);
    for (i=0; i < m_dim; ++i)
    {
        if (pFields[i]->GetWaveSpace())
        {
            pFields[i]->HomogeneousFwdTrans(m_delta[i], m_delta[i]);
        }
        pFields[i]->PhysDeriv(m_delta[i], dvel[0], dvel[1], dvel[2]);
        for (j=0; j < m_dim; ++j)
        {
            if (pFields[i]->GetWaveSpace())
            {
                pFields[i]->HomogeneousBwdTrans(dvel[i], dvel[i]);
            }
            Vmath::Vvtvp(nq, dvel[i], 1, dvel[i], 1, tmp, 1, tmp, 1);
        }
    }
    int id = m_fields.size()-1;
    Vmath::Svtsvtp(nq, m_mu*facStress, tmp, 1, facOld,
        m_fields[id], 1, m_fields[id], 1);
}

void FilterReynoldsStresses::v_PrepareOutput(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    m_fieldMetaData["NumberOfFieldDumps"] =
        boost::lexical_cast<std::string>(m_numSamples);

    // Forward transform and put into m_outFields (except pressure)
    for (int i = 0; i < m_fields.size(); ++i)
    {
        int n = i >= pFields.size() ? 0 : i;
        pFields[n]->FwdTrans_IterPerExp(m_fields[i], m_outFields[i]);
        if (pFields[n]->GetWaveSpace())
        {
            pFields[n]->HomogeneousFwdTrans(m_outFields[i], m_outFields[i]);
        }
    }
}

NekDouble FilterReynoldsStresses::v_GetScale()
{
    if (m_movAvg)
    {
        return 1.0;
    }
    else
    {
        return 1.0 / m_numSamples;
    }
}

}
}
