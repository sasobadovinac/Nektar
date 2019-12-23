///////////////////////////////////////////////////////////////////////////////
//
// File Splines.cpp
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
// Description: Splines definition
//
///////////////////////////////////////////////////////////////////////////////
#include "Splines.h"

namespace Nektar
{
namespace LSIAC
{

Splines::Splines(int deg) : m_deg(deg)
{
}

Splines::Splines()
{
}

void Splines::Initialize(int deg, int n_Bspl,
                         const Array<OneD, NekDouble> &coeffs)
{
    boost::ignore_unused(n_Bspl); // used only for debug purposes
    m_deg = deg;
    m_cpts.clear();
    ASSERTL0(n_Bspl == coeffs.num_elements() &&
             "Calculation depends they being equal");
    for (int i = 0; i < coeffs.num_elements(); i++)
    {
        m_cpts.push_back(coeffs[i]);
    }
    std::reverse(m_cpts.begin(), m_cpts.end());
    m_knots.clear();
    // Number of knots  num_cpts + deg+1;
    int num_knots = coeffs.num_elements() + deg + 1;

    for (int i = 0; i < num_knots; i++)
    {
        m_knots.push_back(-1.0 * ((NekDouble)num_knots - 1.0) / 2.0 +
                          (NekDouble)i);
    }
    this->expandSupport();
}

void Splines::findks(const NekDouble u, int &k, int &s) const
{
    k = 0;
    s = 0;
    for (k = 0; k < m_knots.size(); k++)
    {
        if ((std::abs(m_knots[k] - u) < TOLERENCE))
        {
            s++;
        }
        else if (u < m_knots[k])
        {
            break;
        }
    }
    k--;
}

void Splines::expandSupport()
{
    for (int i = 0; i < m_deg; i++)
    {
        m_knots.insert(m_knots.begin(), m_knots[0] - 1.0);
        m_knots.insert(m_knots.end(), m_knots.back() + 1.0);
    }
    m_cpts.insert(m_cpts.begin(), m_deg, 0.0);
    m_cpts.insert(m_cpts.end(), m_deg, 0.0);
}

bool Splines::validk(int k, int s)
{
    if ((k > m_deg - 1) && (k < (m_knots.size() - m_deg - 1 + s)))
    {
        return true;
    }
    else
    {
        return false;
    }
}

/// This function currently has been only calibrated for Symmetric SIAC.
/// Usage of this function beyond Symmetric SIAC needs to OneSIdedSIAC needs to
/// be calibrated and tested. Will assert if meshShift !=0.0;
void Splines::EvaluateUArr(const vector<NekDouble> &uAr,
                           vector<NekDouble> &solAr,
                           vector<NekDouble> &pts_eval, NekDouble meshScaling,
                           NekDouble meshShift, int m_nthDer)
{
    boost::ignore_unused(meshShift); // For debug purpose.
    ASSERTL0(std::abs(meshShift) < TOLERENCE &&
             "Meshshift != 0.0 means used by Onesided SIAC. Not calibrated for "
             "OneSidedSIAC");

    int k, s;
    NekDouble u;
    NekDouble mulfactorForCpts = 1.0 / std::pow(meshScaling, 1 + m_nthDer);

    for (int i = 0; i < uAr.size(); i++)
    {
        u = uAr[i] / meshScaling;
        this->findks(u, k, s);

        // 2. Handle if  s = deg+1;
        if (m_deg + 1 == s)
        {
            if (k == m_knots.size() - 1)
            {
                solAr[i] = m_cpts[k - m_deg - 1];
            }
            continue;
        }

        // 3. If k is invalid return 0.
        if (!this->validk(k, s))
        {
            solAr[i] = 0.0;
            continue;
        }

        // 4. Evaualte sol.
        memcpy(&pts_eval[0], &m_cpts[k - m_deg],
               (m_deg - s + 1) * sizeof(NekDouble));
        NekDouble alpha = 0.0;
        int ptIndex     = 0;
        for (int r = 1; r <= m_deg - s; r++)
        {
            for (int i = k - m_deg + r; i <= k - s; i++)
            {
                alpha = (u - m_knots[i]) /
                        (m_knots[i + m_deg - r + 1] - m_knots[i]);
                ptIndex           = i - k + m_deg - r;
                pts_eval[ptIndex] = (1 - alpha) * pts_eval[ptIndex] +
                                    alpha * pts_eval[ptIndex + 1];
            }
        }
        solAr[i] = pts_eval[0] * mulfactorForCpts;
    }
    return;
}

void Splines::EvaluateUA(const NekDouble u, NekDouble &sol,
                         vector<NekDouble> &pts_eval)
{
    int k, s;
    this->findks(u, k, s);


    // 2. Handle if  s = deg+1;
    if (m_deg + 1 == s)
    {
        sol = m_cpts[k - m_deg];
        if (k == m_knots.size() - 1)
        {
            sol = m_cpts[k - m_deg - 1];
        }
        return;
    }

    // 3. If k is invalid return 0.
    if (!this->validk(k, s))
    {
        sol = 0.0;
        return;
    }

    // 4. Evaualte sol.
    memcpy(&pts_eval[0], &m_cpts[k - m_deg],
           (m_deg - s + 1) * sizeof(NekDouble));
    NekDouble alpha = 0.0;
    int ptIndex     = 0;
    for (int r = 1; r <= m_deg - s; r++)
    {
        for (int i = k - m_deg + r; i <= k - s; i++)
        {
            alpha =
                (u - m_knots[i]) / (m_knots[i + m_deg - r + 1] - m_knots[i]);
            ptIndex = i - k + m_deg - r;
            pts_eval[ptIndex] =
                (1 - alpha) * pts_eval[ptIndex] + alpha * pts_eval[ptIndex + 1];
        }
    }
    sol = pts_eval[0];
    return;
}

void Splines::EvaluateU(const NekDouble u, NekDouble &sol)
{
    int k, s;
    this->findks(u, k, s);

    // 2. Handle if  s = deg+1;
    if (m_deg + 1 == s)
    {
        sol = m_cpts[k - m_deg];
        if (k == m_knots.size() - 1)
        {
            sol = m_cpts[k - m_deg - 1];
        }
        return;
    }

    // 3. If k is invalid return 0.
    if (!this->validk(k, s))
    {
        sol = 0.0;
        return;
    }

    // 4. Evaualte sol.
    vector<NekDouble> pts_eval(m_cpts.begin() + k - m_deg,
                               m_cpts.begin() + k - s + 1);
    NekDouble alpha = 0.0;
    int ptIndex     = 0;
    for (int r = 1; r <= m_deg - s; r++)
    {
        for (int i = k - m_deg + r; i <= k - s; i++)
        {
            alpha =
                (u - m_knots[i]) / (m_knots[i + m_deg - r + 1] - m_knots[i]);
            ptIndex = i - k + m_deg - r;
            pts_eval[ptIndex] =
                (1 - alpha) * pts_eval[ptIndex] + alpha * pts_eval[ptIndex + 1];
        }
    }
    sol = pts_eval[0];
    return;
}
} // namespace LSIAC
} // namespace Nektar
