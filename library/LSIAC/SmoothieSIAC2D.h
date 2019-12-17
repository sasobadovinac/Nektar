#pragma once
#include "SmoothieSIAC.h"

/// This class can postprocess 2D Meshes.
namespace Nektar
{
namespace LSIAC
{
class SmoothieSIAC2D : public SmoothieSIAC
{
private:
protected:
public:
    SmoothieSIAC2D(const FilterType filter, HandleNekMesh *meshHandle,
                   const int Order, NekDouble meshSpacing = 1.0,
                   const int derivative = 0);

private:
    NekDouble m_meshSpacing;

protected:
    virtual bool v_EvaluateAt(const NekDouble ptsX, const NekDouble PtsY,
                              const NekDouble PtsZ, NekDouble &valX,
                              NekDouble &valY, NekDouble &valZ);

    virtual bool v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                              const NekDouble PtsZ, NekDouble &valX,
                              NekDouble &valY, NekDouble &valZ,
                              Array<OneD, NekDouble> &direction,
                              NekDouble meshSpacing = -1.0, int varNum = 0);

    virtual bool v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                              const NekDouble PtsZ, NekDouble &valX,
                              NekDouble &valY, NekDouble &valZ,
                              vector<Array<OneD, NekDouble>> &directions,
                              NekDouble meshSpacing = -1.0, int varNum = 0);

    virtual bool v_EvaluateAt_NSK_FixedNumBSpl(
        const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        Array<OneD, NekDouble> &direction, NekDouble meshSpacing = -1.0,
        int varNum = 0);

    virtual bool v_EvaluateAt_NUK_MetricTensor(
        const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        Array<OneD, NekDouble> &direction, NekDouble meshSpacing = -1.0,
        int varNum = 0);

    virtual bool v_EvaluateAt_NSK_GivenFilterLength_v1(
        const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        Array<OneD, NekDouble> &direction, NekDouble meshSpacing = -1.0,
        int varNum = 0);

    virtual bool v_EvaluateAt_NSK_GivenFilterLength_v2(
        const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        Array<OneD, NekDouble> &direction, NekDouble meshSpacing = -1.0,
        int varNum = 0);

    virtual bool v_EvaluateAt_SymY(const NekDouble PtsX, const NekDouble PtsY,
                                   const NekDouble PtsZ, NekDouble &valX,
                                   NekDouble &valY, NekDouble &valZ,
                                   Array<OneD, NekDouble> &direction,
                                   NekDouble meshSpacing = -1.0,
                                   int varNum            = 0);

    virtual bool v_EvaluateRecursiveAt(
        const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
        NekDouble &valX, NekDouble &valY, NekDouble &valZ,
        const vector<std::shared_ptr<SmoothieSIAC>> &Sms,
        vector<Array<OneD, NekDouble>> directions,
        const vector<NekDouble> &meshSpacings = vector<NekDouble>(),
        const vector<int> &varNums = vector<int>(), const int curLevel = 0);

    virtual bool v_Cal_NUK_ConstMetricTensor(const NekDouble PtsX,
                                             const NekDouble PtsY,
                                             const NekDouble PtsZ,
                                             const NekDouble meshSpacing,
                                             Array<OneD, NekDouble> &direction,
                                             Array<OneD, NekDouble> &knotVec);
    /*
            virtual bool v_SetupLineForLSIAC( const Array<OneD,NekDouble>
       &direction, const vector<NekDouble> &stPoint, const NekDouble tmin, const
       NekDouble tmax, const int n_quadPts, vector<NekDouble> &HvalT,
       vector<int> &t_GIDs, vector<int> &t_EIDs, Array<OneD,NekDouble>
       &t_LineElm ); virtual bool v_SetupLineForLSIAC_ReSamp( const
       Array<OneD,NekDouble> &direction, const vector<NekDouble> &stPoint, const
       NekDouble tmin, const NekDouble tmax, const int n_quadPts, const int
       n_quadPts_Resample, vector<NekDouble> &HvalT, vector<int> &t_GIDs,
       vector<int> &t_EIDs, Array<OneD,NekDouble> &t_LineElm,
       Array<OneD,NekDouble> &t_LineElm_Resample );

            virtual bool v_GetVLineForLSIAC(	const int n_quadPts, const
       vector<NekDouble> &stPoint, const Array<OneD,NekDouble> &direction, const
       vector<NekDouble> &HvalT, const vector<int> &t_GIDs, const vector<int>
       &t_EIDs, const Array<OneD,NekDouble> &t_LineElm, Array<OneD,NekDouble>
       tv_LineElm, int varNum ); virtual bool v_GetVLineForLSIAC_resample(
       const int n_quadPts,const int n_quadPts_resample, const vector<NekDouble>
       &stPoint, const Array<OneD,NekDouble> &direction, const vector<NekDouble>
       &HvalT, const vector<int> &t_GIDs, const vector<int> &t_EIDs, const
       Array<OneD,NekDouble> &t_LineElm, const Array<OneD,NekDouble> tv_LineElm,
                                                            const
       Array<OneD,NekDouble> &t_LineElm_resample, Array<OneD,NekDouble>
       tv_LineElm_resample, int varNum ); virtual bool v_EvaluateUsingLineAt(
       const vector<NekDouble> &stPoint, const Array<OneD,NekDouble> &direction,
                                                                    const int
       n_quadPts, const NekDouble meshScaling, const vector<NekDouble> &tparams,
                                                                    vector<NekDouble>
       &tvals, int varNum); virtual bool v_EvaluateUsingLineAt_v2( const
       vector<NekDouble> &stPoint, const Array<OneD,NekDouble> &direction, const
       int n_quadPts,const int n_quadPts_resample, const NekDouble meshScaling,
       const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int varNum);
            virtual bool v_EvaluateUsingLineAt_v3( const vector<NekDouble>
       &stPoint, const Array<OneD,NekDouble> &direction, const int
       n_quadPts,const int n_quadPts_resample, const NekDouble meshScaling,
       const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int varNum);
    */
    virtual bool v_EvaluateL2UsingLineAt_v3(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const int n_quadPts_resample, const NekDouble meshScaling,
        const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int varNum);

    virtual bool v_EvaluateL2UsingLineAt(
        const vector<NekDouble> &stPoint,
        const Array<OneD, NekDouble> &direction, const int n_quadPts,
        const NekDouble meshScaling, const vector<NekDouble> &tparams,
        vector<NekDouble> &tvals, int varNum);
};

} // namespace LSIAC
} // namespace Nektar
