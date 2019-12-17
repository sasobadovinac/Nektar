#include "CentralBSplines.h"

namespace Nektar
{
namespace LSIAC
{

CentralBSplines::CentralBSplines(int Order, NekDouble Shift, NekDouble scale)
    : GeneralBSplines(Order), m_shift(Shift), m_scale(scale)
{
    Array<OneD, NekDouble> knots(Order + 1, 0.0);
    for (int i = 0; i < Order + 1; i++)
    {
        knots[i] = -Order / 2.0 + i;
    }
    this->SetKnotVector(knots);
    GeneralBSplines::SetKnotVector(knots);
}

CentralBSplines::CentralBSplines(int Order)
    : GeneralBSplines(Order), m_shift(0.0), m_scale(1.0)
{
    Array<OneD, NekDouble> knots(Order + 1, 0.0);
    for (int i = 0; i < Order + 1; i++)
    {
        knots[i] = -Order / 2.0 + i;
    }
    this->SetKnotVector(knots);
    GeneralBSplines::SetKnotVector(knots);
}

CentralBSplines::CentralBSplines(int Order, NekDouble shift)
    : GeneralBSplines(Order), m_shift(shift), m_scale(1.0)
{
    Array<OneD, NekDouble> knots(Order + 1, 0.0);
    for (int i = 0; i < Order + 1; i++)
    {
        knots[i] = -Order / 2.0 + i;
    }
    this->SetKnotVector(knots);
    GeneralBSplines::SetKnotVector(knots);
}

bool CentralBSplines::EvaluateBSplines(const Array<OneD, NekDouble> &t_pos,
                                       Array<OneD, NekDouble> &t_val,
                                       const NekDouble shift,
                                       const NekDouble meshScaling) const
{
    // For central BSplines. The j is always zero.
    // The spline needs to be shifted and evaluated at x_pos locations.
    GeneralBSplines::EvaluateBSplines(t_pos, 0, t_val, shift, meshScaling);
    return true;
}

} // namespace LSIAC
} // namespace Nektar
