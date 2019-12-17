#pragma once
#include "GeneralBSplines.h"

/// The class evaluates CentralB-Splines at given locations.
/**	This calss automatically calculates knot positions, when using the order
   specified. All the knot positions calcualted will range -(k+1)/2 to k+1)/2
*/
namespace Nektar
{
namespace LSIAC
{
class CentralBSplines : public GeneralBSplines
{

private:
protected:
public:
    NekDouble m_shift;
    NekDouble m_scale;

    CentralBSplines(int Order);
    CentralBSplines(int Order, NekDouble shift);
    CentralBSplines(int Order, NekDouble shift, NekDouble scale);

    bool GetShift(NekDouble &shift) const;
    bool SetShift(NekDouble shift);

    bool GetScale(NekDouble &scale) const;
    bool SetScale(NekDouble scale);

    bool SetOrder(int Order);
    int GetOrder() const;

    //	bool EvaluateBSplines( const Array<OneD,NekDouble> &t_pos,const
    // NekDouble shift, Array<OneD,NekDouble> &t_vals)const;
    bool EvaluateBSplines(const Array<OneD, NekDouble> &t_pos,
                          Array<OneD, NekDouble> &t_vals,
                          const NekDouble shift       = 0.0,
                          const NekDouble meshScaling = 1.0) const;
    //	bool EvaluateBSplines( Array<OneD,NekDouble> &t_pos,
    // Array<OneD,NekDouble> &t_vals);
};
} // namespace LSIAC
} // namespace Nektar
