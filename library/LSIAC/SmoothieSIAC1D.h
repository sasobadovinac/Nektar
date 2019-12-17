#pragma once
#include "HandleNekMesh.h"
#include "HandleNekMesh1D.h"
#include "SmoothieSIAC.h"
#include "SymmetricSIAC.h"

/// This class can postprocess 1D Meshes.

namespace Nektar
{
namespace LSIAC
{
class SmoothieSIAC1D : public SmoothieSIAC
{
private:
protected:
public:
    //	SmoothieSIAC1D( const FilterType filter,HandleNekMesh* meshHandle,
    //								const int
    // Order);

    SmoothieSIAC1D(const FilterType filter, HandleNekMesh *meshHandle,
                   const int Order, NekDouble meshSpacing = 1.0,
                   const int derivative = 0);

private:
    NekDouble m_meshSpacing;

protected:
    virtual bool v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                              const NekDouble PtsZ, NekDouble &valX,
                              NekDouble &valY, NekDouble &valZ);

    virtual bool v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                              const NekDouble PtsZ, NekDouble &valX,
                              NekDouble &valY, NekDouble &valZ,
                              Array<OneD, NekDouble> &direction,
                              NekDouble meshSpacing = -1.0, int varNum = 0);

    virtual bool v_EvaluateNonSymAt(const NekDouble PtsX, const NekDouble PtsY,
                                    const NekDouble PtsZ, NekDouble &valX,
                                    NekDouble &valY, NekDouble &valZ);

    virtual bool v_EvaluateAt_NUK_MetricTensor(
        const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        Array<OneD, NekDouble> &direction, NekDouble meshSpacing = -1.0,
        int varNum = 0);

    virtual bool v_Cal_NUK_ConstMetricTensor(const NekDouble PtsX,
                                             const NekDouble PtsY,
                                             const NekDouble PtsZ,
                                             const NekDouble meshSpacing,
                                             Array<OneD, NekDouble> &direction,
                                             Array<OneD, NekDouble> &knotVec);
};
} // namespace LSIAC
} // namespace Nektar
