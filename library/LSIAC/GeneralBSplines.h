#pragma once
#include "BSplines.h"
#include <iostream>
#include <vector>

/// This class evaluates any general BSplines at given location when knots and
/// order are specified.

/** To evaluate General BSplines at any location, one needs knots and Order of
 * Bsplines.
 */

namespace Nektar
{
namespace LSIAC
{

class GeneralBSplines : public BSplines
{
    // data
public:
    Array<OneD, NekDouble> m_knotVector;
    int m_order;
    // functions
private:
    bool BSplinesBasis(const Array<OneD, NekDouble> &t_pos, const int k,
                       const int j, Array<OneD, NekDouble> &t_val,
                       const NekDouble scaling     = 0.0,
                       const NekDouble meshScaling = 1.0) const;

    bool BSplinesBasis(const Array<OneD, NekDouble> &t_pos,
                       const vector<NekDouble> &kvector, const int k,
                       const int j, Array<OneD, NekDouble> &t_val,
                       const NekDouble scaling     = 0.0,
                       const NekDouble meshScaling = 1.0) const;

protected:
public:
    GeneralBSplines(const int order);
    GeneralBSplines(const Array<OneD, NekDouble> &knots, const int order);

    bool SetKnotVector(const Array<OneD, NekDouble> &knots);

    bool GetKnotVector(Array<OneD, NekDouble> &knots);

    bool SetOrder(const int order);

    int GetOrder() const;

    // Not used yet!
    bool EvaluateBSplines(const Array<OneD, NekDouble> &t_pos,
                          const Array<OneD, NekDouble> &knots, const int j_th,
                          Array<OneD, NekDouble> &t_values);

    // Not used yet!
    bool EvaluateBSplines(const Array<OneD, NekDouble> &t_pos,
                          const Array<OneD, NekDouble> &knots, const int j_th,
                          const NekDouble shift,
                          Array<OneD, NekDouble> &t_values);

    bool EvaluateBSplines(const Array<OneD, NekDouble> &t_pos, const int j_th,
                          Array<OneD, NekDouble> &t_values,
                          const NekDouble shift       = 0.0,
                          const NekDouble meshScaling = 1.0) const;

    bool EvaluateBSplines(const Array<OneD, NekDouble> &t_pos,
                          const std::vector<NekDouble> &kvec, const int j_th,
                          Array<OneD, NekDouble> &t_values,
                          const NekDouble shift       = 0.0,
                          const NekDouble meshScaling = 1.0) const;
};
} // namespace LSIAC
} // namespace Nektar
