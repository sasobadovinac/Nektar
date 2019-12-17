#pragma once
#include "LSIACPostProcessor.h"

/// This class evaluates Splines.
/*** This class is a base class implemented here to take advantage of de-boorg
   algorithm to calcualte splines.
*/

namespace Nektar
{
namespace LSIAC
{
class Splines : LSIACPostProcessor
{
public:
    int m_deg;
    vector<NekDouble> m_cpts;
    vector<NekDouble> m_knots;

    Splines();
    Splines(int deg);
    void findks(const NekDouble u, int &k, int &s) const;
    void expandSupport();
    bool validk(int k, int s);

    void Initialize(int deg, int n_Bspl, const Array<OneD, NekDouble> &coeffs);
    void EvaluateUArr(const vector<NekDouble> &uAr, vector<NekDouble> &solAr,
                      vector<NekDouble> &pts_eval, NekDouble meshScale = 1.0,
                      NekDouble meshShift = 0.0, int m_nthDer = 0);
    void EvaluateUA(const NekDouble u, NekDouble &sol,
                    vector<NekDouble> &pts_eval);
    void EvaluateU(const NekDouble u, NekDouble &sol);
};

} // namespace LSIAC
} // namespace Nektar
