#include "SmoothieSIAC.h"
#include <algorithm>
/*
SmoothieSIAC::SmoothieSIAC( const SmoothieSIAC::FilterType filter,const
HandleMesh meshHandle, const int Order, Array<OneD,NekDouble> &direction):
                                                                        m_fitlerType(filter),
m_MeshHandle(meshHandle), m_order(order)
{
}
*/

namespace Nektar
{
namespace LSIAC
{
bool compare2NekDoubles(NekDouble x, NekDouble y)
{
    return ((x - TOLERENCE < y) && (x + TOLERENCE > y));
}

bool SmoothieSIAC::mergeBreakPts(const vector<NekDouble> &HvalT,
                                 const vector<NekDouble> &SvalT,
                                 vector<NekDouble> &TvalT)
{
    // insert two vector into TvalT
    TvalT = HvalT;
    TvalT.insert(TvalT.end(), SvalT.begin(), SvalT.end());
    // sort;
    sort(TvalT.begin(), TvalT.end());
    // unique;
    std::vector<NekDouble>::iterator it;
    it = std::unique(TvalT.begin(), TvalT.end(), compare2NekDoubles);
    TvalT.resize(std::distance(TvalT.begin(), it));

    return true;
}

bool SmoothieSIAC::mergeBreakPtsAdv(const vector<NekDouble> &HvalT,
                                    const vector<NekDouble> &SvalT,
                                    vector<NekDouble> &TvalT,
                                    const vector<int> &t_h_GIDs,
                                    const vector<int> &t_h_EIDs,
                                    vector<int> &t_GIDs, vector<int> &t_EIDs)
{
    mergeBreakPts(HvalT, SvalT, TvalT);
    t_GIDs.clear();
    t_EIDs.clear();
    t_GIDs.resize(TvalT.size(), -2);
    t_EIDs.resize(TvalT.size(), -2);
    /*
            cout << "HvalT "<< endl;
            printNekArray(HvalT,0);
            cout << "SvalT "<< endl;
            printNekArray(SvalT,0);
            cout << "TvalT "<< endl;
            printNekArray(TvalT,0);
            cout << "t_h_GIDs "<< endl;
            printNekArray(t_h_GIDs,0);
            cout << "t_h_EIDs "<< endl;
            printNekArray(t_h_EIDs,0);
    */
    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        NekDouble tval = (TvalT[i] + TvalT[i + 1]) / 2.0;
        // cout << tval << endl;
        for (int j = 0; j < HvalT.size() - 1; j++)
        {
            if (tval > HvalT[j] && tval < HvalT[j + 1])
            {
                t_GIDs[i] = t_h_GIDs[j];
                t_EIDs[i] = t_h_EIDs[j];
                continue;
            }
        }
        assert(true && "Something is wrong. Should never come till here.");
    }
    /*
            cout << "t_GIDs "<< endl;
            printNekArray(t_GIDs,0);
            cout << "t_EIDs "<< endl;
            printNekArray(t_EIDs,0);
    */
    return true;
}

bool SmoothieSIAC::v_EvaluateLineForLSIAC_v1_dynScaling(
    const int n_quadPts, const vector<NekDouble> &tparams,
    const vector<NekDouble> &t_dynScaling, const NekDouble t_mesh_min,
    const NekDouble t_mesh_max, const Array<OneD, NekDouble> &t_LineElm,
    const Array<OneD, NekDouble> &tv_LineElm, const vector<NekDouble> &HvalT,
    vector<NekDouble> &tvals)
{
    // Set up one more array for quadrature weights or jacobian values.
    NekDouble t_min, t_max, meshTShift = -1.0;
    LibUtilities::PointsKey quadPointsKey(
        n_quadPts, Nektar::LibUtilities::eGaussGaussLegendre);
    Array<OneD, NekDouble> t_quad_weights =
        LibUtilities::PointsManager()[quadPointsKey]->GetW();
    Array<OneD, NekDouble> t_quadPts(n_quadPts), t_quad_vals(n_quadPts);

    for (int i = 0; i < tparams.size(); i++)
    {
        NekDouble meshSpacing = t_dynScaling[i];
        // Figuring out which type of filter need to be used.
        bool b_symMesh = true;
        NekDouble t    = tparams[i];
        m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, t_min, t_max);
        if (t + t_min < t_mesh_min && t + t_max > t_mesh_max)
        {
            assert("No filter can be applied");
            return false;
        }
        if (t + t_min < t_mesh_min)
        { // one sided filter need to be applied.
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min,
                                                      t_max);
            meshTShift = t_mesh_min - t - t_min;
            b_symMesh  = false;
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min, t_max,
                                                      meshTShift);
        }
        if (t + t_max > t_mesh_max)
        { // one sided filter need to be applied.
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min,
                                                      t_max);
            meshTShift = t_mesh_max - t - t_max;
            b_symMesh  = false;
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min, t_max,
                                                      meshTShift);
        }

        // cout << "for parameter t " << t << "\t tmin: \t" << t_min << "\t
        // tmax: \t" << t_max << endl;

        assert(t + t_min + TOLERENCE > t_mesh_min &&
               "Above code should have fixed this issue");
        assert(t + t_max - TOLERENCE < t_mesh_max &&
               "Above code should have fixed this issue");
        // if symmetric mesh. Could be true for both.
        int startElmIndex, endElmIndex;
        // start ElementId and End Element Id;
        // const std::vector<NekDouble>::iterator
        const auto low =
            std::lower_bound(HvalT.begin(), HvalT.end(), t + t_min + TOLERENCE);
        const auto up =
            std::upper_bound(HvalT.begin(), HvalT.end(), t + t_max - TOLERENCE);
        // int stIndex = low - HvalT.begin()-1;
        // int edIndex = up -HvalT.begin()-1;
        int stIndex   = std::distance(HvalT.begin(), low);
        int edIndex   = std::distance(HvalT.begin(), up);
        startElmIndex = stIndex;
        endElmIndex   = edIndex;
        //		if (startElmIndex != stIndex)
        //		cout << "startElmIndex =\t" << startElmIndex <<
        //"\tstIndex
        //=\t" <<stIndex <<endl; 	if (endElmIndex != edIndex)
        // cout << "endElmIndex =\t" << endElmIndex << "\tedIndex =\t" <<edIndex
        //<<endl; 	cout << "for parameter t " << t << "\t sElIndex: \t" <<
        // startElmIndex<< "\t tmax: \t" << endElmIndex<< endl;

        assert(startElmIndex <= endElmIndex && "Wrong");
        NekDouble sum = 0.0;
        for (int el = startElmIndex; el <= endElmIndex; el++)
        {
            for (int qp = 0; qp < n_quadPts; qp++)
            {
                t_quadPts[qp] = t_LineElm[el * n_quadPts + qp] - t;
            }
            if (b_symMesh)
            {
                m_siacFilterPtrs[m_SymID]->EvaluateFilter(
                    t_quadPts, t_quad_vals, meshSpacing);
            }
            else
            {
                m_siacFilterPtrs[m_OneID]->EvaluateFilter(
                    t_quadPts, t_quad_vals, meshSpacing, meshTShift, true);
            }
            NekDouble integral = 0.0;
            for (int qp = 0; qp < n_quadPts; qp++)
            {
                integral +=
                    0.5 * t_quad_weights[qp] * tv_LineElm[el * n_quadPts + qp] *
                    std::abs(HvalT[el + 1] - HvalT[el]) * t_quad_vals[qp];
            }
            sum += integral;
        }
        tvals.push_back(sum);
    }
    return true;
}
bool SmoothieSIAC::v_EvaluateLineForLSIAC_v1(
    const int n_quadPts, const vector<NekDouble> &tparams,
    const NekDouble meshSpacing, const NekDouble t_mesh_min,
    const NekDouble t_mesh_max, const Array<OneD, NekDouble> &t_LineElm,
    const Array<OneD, NekDouble> &tv_LineElm, const vector<NekDouble> &HvalT,
    vector<NekDouble> &tvals)
{
    // Set up one more array for quadrature weights or jacobian values.
    // tW_LineElm
    // tJ_LineElm
    // NekDouble meshSpacing = 0.1; // NekDouble.
    NekDouble t_min, t_max, meshTShift;
    LibUtilities::PointsKey quadPointsKey(
        n_quadPts, Nektar::LibUtilities::eGaussGaussLegendre);
    // Array<OneD,NekDouble> quad_points =
    // LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    Array<OneD, NekDouble> t_quad_weights =
        LibUtilities::PointsManager()[quadPointsKey]->GetW();
    Array<OneD, NekDouble> t_quadPts(n_quadPts), t_quad_vals(n_quadPts);

    for (int i = 0; i < tparams.size(); i++)
    {
        // Figuring out which type of filter need to be used.
        bool b_symMesh = true;
        NekDouble t    = tparams[i];
        m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, t_min, t_max);
        if (t + t_min < t_mesh_min && t + t_max > t_mesh_max)
        {
            assert("No filter can be applied");
            return false;
        }
        if (t + t_min < t_mesh_min)
        { // one sided filter need to be applied.
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min,
                                                      t_max);
            meshTShift = t_mesh_min - t - t_min;
            b_symMesh  = false;
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min, t_max,
                                                      meshTShift);
        }
        if (t + t_max > t_mesh_max)
        { // one sided filter need to be applied.
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min,
                                                      t_max);
            meshTShift = t_mesh_max - t - t_max;
            b_symMesh  = false;
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min, t_max,
                                                      meshTShift);
        }

        // cout << "for parameter t " << t << "\t tmin: \t" << t_min << "\t
        // tmax: \t" << t_max << endl;

        assert(t + t_min + TOLERENCE > t_mesh_min &&
               "Above code should have fixed this issue");
        assert(t + t_max - TOLERENCE < t_mesh_max &&
               "Above code should have fixed this issue");
        // if symmetric mesh. Could be true for both.
        int startElmIndex, endElmIndex;
        // start ElementId and End Element Id;
        // const std::vector<NekDouble>::iterator
        const auto low =
            std::lower_bound(HvalT.begin(), HvalT.end(), t + t_min + TOLERENCE);
        const auto up =
            std::upper_bound(HvalT.begin(), HvalT.end(), t + t_max - TOLERENCE);
        // int stIndex = low - HvalT.begin()-1;
        // int edIndex = up -HvalT.begin()-1;
        int stIndex   = std::distance(HvalT.begin(), low);
        int edIndex   = std::distance(HvalT.begin(), up);
        startElmIndex = stIndex;
        endElmIndex   = edIndex;
        //		if (startElmIndex != stIndex)
        //		cout << "startElmIndex =\t" << startElmIndex <<
        //"\tstIndex
        //=\t" <<stIndex <<endl; 	if (endElmIndex != edIndex)
        // cout << "endElmIndex =\t" << endElmIndex << "\tedIndex =\t" <<edIndex
        //<<endl;
        // cout << "for parameter t " << t << "\t sElIndex: \t" <<
        // startElmIndex<< "\t tmax: \t" << endElmIndex<< endl;

        assert(startElmIndex <= endElmIndex && "Wrong");
        NekDouble sum = 0.0;
        for (int el = startElmIndex; el <= endElmIndex; el++)
        {
            for (int qp = 0; qp < n_quadPts; qp++)
            {
                t_quadPts[qp] = t_LineElm[el * n_quadPts + qp] - t;
            }
            if (b_symMesh)
            {
                m_siacFilterPtrs[m_SymID]->EvaluateFilter(
                    t_quadPts, t_quad_vals, meshSpacing);
            }
            else
            {
                m_siacFilterPtrs[m_OneID]->EvaluateFilter(
                    t_quadPts, t_quad_vals, meshSpacing, meshTShift, true);
            }
            NekDouble integral = 0.0;
            for (int qp = 0; qp < n_quadPts; qp++)
            {
                integral +=
                    0.5 * t_quad_weights[qp] * tv_LineElm[el * n_quadPts + qp] *
                    std::abs(HvalT[el + 1] - HvalT[el]) * t_quad_vals[qp];
            }
            sum += integral;
        }
        tvals.push_back(sum);
    }
    return false;
}

bool SmoothieSIAC::v_EvaluateLineForLSIAC_v3_dynScaling(
    const int n_quadPts, const int n_quadPts_resample,
    const vector<NekDouble> &tparams, const vector<NekDouble> &t_dynScaling,
    const NekDouble t_mesh_min, const NekDouble t_mesh_max,
    const Array<OneD, NekDouble> &t_LineElm,
    const Array<OneD, NekDouble> &tv_LineElm, const vector<NekDouble> &HvalT,
    vector<NekDouble> &tvals)
{
    boost::ignore_unused(t_LineElm); // not updated
    LibUtilities::PointsKey quadPointsKey(
        n_quadPts, Nektar::LibUtilities::eGaussGaussLegendre);
    LibUtilities::BasisKey bk = LibUtilities::BasisKey(
        Nektar::LibUtilities::eGauss_Lagrange, n_quadPts, quadPointsKey);
    Nektar::StdRegions::StdExpansionSharedPtr segExp_std =
        MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(bk);

    LibUtilities::PointsKey quadPointsKey_resample(
        n_quadPts_resample, Nektar::LibUtilities::eGaussGaussLegendre);
    Array<OneD, NekDouble> q_quadR_weights =
        LibUtilities::PointsManager()[quadPointsKey_resample]->GetW();
    Array<OneD, NekDouble> q_quadR_Pts =
        LibUtilities::PointsManager()[quadPointsKey_resample]->GetZ();
    Array<OneD, NekDouble> t_quadPts(n_quadPts_resample),
        t_quad_vals(n_quadPts_resample);
    //	Array<OneD,NekDouble> t_quad(n_quadPts_resample);

    Array<OneD, NekDouble> lcoord(3, 0.0);

    for (int ii = 0; ii < tparams.size(); ii++)
    {
        // Figuring out which type of filter need to be used.
        bool b_symMesh        = true;
        NekDouble meshSpacing = t_dynScaling[ii];
        NekDouble t           = tparams[ii];
        NekDouble t_min, t_max, meshTShift = -1.0;
        m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, t_min, t_max);
        if (t + t_min < t_mesh_min && t + t_max > t_mesh_max)
        {
            assert("No filter can be applied");
            return false;
        }
        if (t + t_min < t_mesh_min)
        { // one sided filter need to be applied.
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min,
                                                      t_max);
            meshTShift = t_mesh_min - t - t_min;
            b_symMesh  = false;
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min, t_max,
                                                      meshTShift);
        }
        if (t + t_max > t_mesh_max)
        { // one sided filter need to be applied.
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min,
                                                      t_max);
            meshTShift = t_mesh_max - t - t_max;
            b_symMesh  = false;
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min, t_max,
                                                      meshTShift);
        }

        // cout << "for parameter t " << t << "\t tmin: \t" << t_min << "\t
        // tmax: \t" << t_max << endl;

        assert(t + t_min + TOLERENCE > t_mesh_min &&
               "Above code should have fixed this issue");
        assert(t + t_max - TOLERENCE < t_mesh_max &&
               "Above code should have fixed this issue");

        // if symmetric mesh. Could be true for both.
        // start ElementId and End Element Id;
        int startElmIndex, endElmIndex;
        const auto low =
            std::lower_bound(HvalT.begin(), HvalT.end(), t + t_min + TOLERENCE);
        const auto up =
            std::upper_bound(HvalT.begin(), HvalT.end(), t + t_max - TOLERENCE);
        // int stIndex = low - HvalT.begin()-1;
        // int edIndex = up -HvalT.begin()-1;
        int stIndex   = std::distance(HvalT.begin(), low);
        int edIndex   = std::distance(HvalT.begin(), up);
        startElmIndex = stIndex;
        endElmIndex   = edIndex;

        // cout << "for parameter t " << t << "\t sElIndex: \t" <<
        // startElmIndex<< "\t tmax: \t" << endElmIndex<< endl;

        assert(startElmIndex <= endElmIndex && "Wrong");
        NekDouble sum = 0.0;
        vector<NekDouble> SvalT, LvalT, TvalT;
        if (b_symMesh)
        {
            m_siacFilterPtrs[m_SymID]->GetBreakPts(meshSpacing, SvalT);
        }
        else
        {
            m_siacFilterPtrs[m_OneID]->GetBreakPts(meshSpacing, SvalT,
                                                   meshTShift);
        }

        // Insert all mesh breakpoints including min and max using
        // HvalT into LvalT
        LvalT.clear();
        LvalT.push_back(t_min);
        for (int j = startElmIndex; j <= endElmIndex + 1; j++)
        {
            if ((HvalT[j] > t + t_min) && (HvalT[j] < t + t_max))
            {
                LvalT.push_back(HvalT[j] - t);
            }
        }
        LvalT.push_back(t_max);
        mergeBreakPts(SvalT, LvalT, TvalT);

        // m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, diretion, TvalT,
        // t_GIDs,t_EIDs); Need to be done.
        for (int i = 0; i < TvalT.size() - 1; i++)
        {
            NekDouble a = TvalT[i];
            NekDouble b = TvalT[i + 1];
            //	int gID = t_GIDs[i]; need elmIndex on line. Hence code below.
            //	int eID = t_EIDs[i];
            int elmIndex = -1;
            for (int j = startElmIndex; j <= endElmIndex + 1; j++)
            {
                if (HvalT[j] > t + (a + b) / 2.0)
                {
                    elmIndex = j - 1;
                    break;
                }
            }
            assert(elmIndex != -1 && "Elm Index is wrong");
            for (int j = 0; j < n_quadPts_resample; j++)
            {
                t_quadPts[j] = (b - a) * (q_quadR_Pts[j] + 1.0) / 2.0 + a;
            }
            if (b_symMesh)
            {
                m_siacFilterPtrs[m_SymID]->EvaluateFilter(
                    t_quadPts, t_quad_vals, meshSpacing);
            }
            else
            {
                m_siacFilterPtrs[m_OneID]->EvaluateFilter(
                    t_quadPts, t_quad_vals, meshSpacing, meshTShift, true);
            }
            //	m_meshHandlePtr->EvaluateAt(t_xyz_vals);

            NekDouble integral = 0.0;
            Array<OneD, NekDouble> elv =
                tv_LineElm.CreateWithOffset(tv_LineElm, elmIndex * n_quadPts);
            for (int k = 0; k < n_quadPts_resample; k++)
            {
                // need to calcuate value of t_xyz_vals.
                // scaling of HvalT[j-1] to -1 and HvalT[j] to +1
                //	lcoord[0] = (t+ t_quadPts[k] - HvalT[elmIndex-1])/
                //(HvalT[elmIndex]-HvalT[elmIndex-1])*2.0 -1.0;
                lcoord[0] = (t + t_quadPts[k] - HvalT[elmIndex]) /
                                (HvalT[elmIndex + 1] - HvalT[elmIndex]) * 2.0 -
                            1.0;
                NekDouble t_xyz_val = segExp_std->PhysEvaluate(lcoord, elv);
                // SegExpasionPhysEvaluate
                integral += 0.5 * q_quadR_weights[k] * t_quad_vals[k] *
                            std::abs(b - a) * t_xyz_val;
            }
            sum += integral;
        }
        tvals.push_back(sum);
    }
    return true;
}
bool SmoothieSIAC::v_EvaluateLineForLSIAC_v3(
    const int n_quadPts, const int n_quadPts_resample,
    const vector<NekDouble> &tparams, const NekDouble meshSpacing,
    const NekDouble t_mesh_min, const NekDouble t_mesh_max,
    const Array<OneD, NekDouble> &t_LineElm,
    const Array<OneD, NekDouble> &tv_LineElm, const vector<NekDouble> &HvalT,
    vector<NekDouble> &tvals)
{
    boost::ignore_unused(t_LineElm); // not updated
    LibUtilities::PointsKey quadPointsKey(
        n_quadPts, Nektar::LibUtilities::eGaussGaussLegendre);
    LibUtilities::BasisKey bk = LibUtilities::BasisKey(
        Nektar::LibUtilities::eGauss_Lagrange, n_quadPts, quadPointsKey);
    Nektar::StdRegions::StdExpansionSharedPtr segExp_std =
        MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(bk);

    LibUtilities::PointsKey quadPointsKey_resample(
        n_quadPts_resample, Nektar::LibUtilities::eGaussGaussLegendre);
    Array<OneD, NekDouble> q_quadR_weights =
        LibUtilities::PointsManager()[quadPointsKey_resample]->GetW();
    Array<OneD, NekDouble> q_quadR_Pts =
        LibUtilities::PointsManager()[quadPointsKey_resample]->GetZ();
    Array<OneD, NekDouble> t_quadPts(n_quadPts_resample),
        t_quad_vals(n_quadPts_resample);
    //	Array<OneD,NekDouble> t_quad(n_quadPts_resample);

    Array<OneD, NekDouble> lcoord(3, 0.0);

    for (int ii = 0; ii < tparams.size(); ii++)
    {
        // Figuring out which type of filter need to be used.
        bool b_symMesh = true;
        NekDouble t    = tparams[ii];
        NekDouble t_min, t_max, meshTShift = -1.0;
        m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, t_min, t_max);
        if (t + t_min < t_mesh_min && t + t_max > t_mesh_max)
        {
            assert("No filter can be applied");
            return false;
        }
        if (t + t_min < t_mesh_min)
        { // one sided filter need to be applied.
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min,
                                                      t_max);
            meshTShift = t_mesh_min - t - t_min;
            b_symMesh  = false;
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min, t_max,
                                                      meshTShift);
        }
        if (t + t_max > t_mesh_max)
        { // one sided filter need to be applied.
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min,
                                                      t_max);
            meshTShift = t_mesh_max - t - t_max;
            b_symMesh  = false;
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, t_min, t_max,
                                                      meshTShift);
        }

        // cout << "for parameter t " << t << "\t tmin: \t" << t_min << "\t
        // tmax: \t" << t_max << endl;

        assert(t + t_min + TOLERENCE > t_mesh_min &&
               "Above code should have fixed this issue");
        assert(t + t_max - TOLERENCE < t_mesh_max &&
               "Above code should have fixed this issue");

        // if symmetric mesh. Could be true for both.
        // start ElementId and End Element Id;
        int startElmIndex, endElmIndex;
        const auto low =
            std::lower_bound(HvalT.begin(), HvalT.end(), t + t_min + TOLERENCE);
        const auto up =
            std::upper_bound(HvalT.begin(), HvalT.end(), t + t_max - TOLERENCE);
        // int stIndex = low - HvalT.begin()-1;
        // int edIndex = up -HvalT.begin()-1;
        int stIndex   = std::distance(HvalT.begin(), low);
        int edIndex   = std::distance(HvalT.begin(), up);
        startElmIndex = stIndex;
        endElmIndex   = edIndex;

        // cout << "for parameter t " << t << "\t sElIndex: \t" <<
        // startElmIndex<< "\t tmax: \t" << endElmIndex<< endl;

        assert(startElmIndex <= endElmIndex && "Wrong");
        NekDouble sum = 0.0;
        vector<NekDouble> SvalT, LvalT, TvalT;
        if (b_symMesh)
        {
            m_siacFilterPtrs[m_SymID]->GetBreakPts(meshSpacing, SvalT);
        }
        else
        {
            m_siacFilterPtrs[m_OneID]->GetBreakPts(meshSpacing, SvalT,
                                                   meshTShift);
        }

        // Insert all mesh breakpoints including min and max using
        // HvalT into LvalT
        LvalT.clear();
        LvalT.push_back(t_min);
        for (int j = startElmIndex; j <= endElmIndex + 1; j++)
        {
            if ((HvalT[j] > t + t_min) && (HvalT[j] < t + t_max))
            {
                LvalT.push_back(HvalT[j] - t);
            }
        }
        LvalT.push_back(t_max);
        mergeBreakPts(SvalT, LvalT, TvalT);

        // m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, diretion, TvalT,
        // t_GIDs,t_EIDs); Need to be done.
        for (int i = 0; i < TvalT.size() - 1; i++)
        {
            NekDouble a = TvalT[i];
            NekDouble b = TvalT[i + 1];
            //	int gID = t_GIDs[i]; need elmIndex on line. Hence code below.
            //	int eID = t_EIDs[i];
            int elmIndex = -1;
            for (int j = startElmIndex; j <= endElmIndex + 1; j++)
            {
                if (HvalT[j] > t + (a + b) / 2.0)
                {
                    elmIndex = j - 1;
                    break;
                }
            }
            assert(elmIndex != -1 && "Elm Index is wrong");
            for (int j = 0; j < n_quadPts_resample; j++)
            {
                t_quadPts[j] = (b - a) * (q_quadR_Pts[j] + 1.0) / 2.0 + a;
            }
            if (b_symMesh)
            {
                m_siacFilterPtrs[m_SymID]->EvaluateFilter(
                    t_quadPts, t_quad_vals, meshSpacing);
            }
            else
            {
                m_siacFilterPtrs[m_OneID]->EvaluateFilter(
                    t_quadPts, t_quad_vals, meshSpacing, meshTShift, true);
            }
            //	m_meshHandlePtr->EvaluateAt(t_xyz_vals);

            NekDouble integral = 0.0;
            Array<OneD, NekDouble> elv =
                tv_LineElm.CreateWithOffset(tv_LineElm, elmIndex * n_quadPts);
            for (int k = 0; k < n_quadPts_resample; k++)
            {
                // need to calcuate value of t_xyz_vals.
                // scaling of HvalT[j-1] to -1 and HvalT[j] to +1
                //	lcoord[0] = (t+ t_quadPts[k] - HvalT[elmIndex-1])/
                //(HvalT[elmIndex]-HvalT[elmIndex-1])*2.0 -1.0;
                lcoord[0] = (t + t_quadPts[k] - HvalT[elmIndex]) /
                                (HvalT[elmIndex + 1] - HvalT[elmIndex]) * 2.0 -
                            1.0;
                NekDouble t_xyz_val = segExp_std->PhysEvaluate(lcoord, elv);
                // SegExpasionPhysEvaluate
                integral += 0.5 * q_quadR_weights[k] * t_quad_vals[k] *
                            std::abs(b - a) * t_xyz_val;
            }
            sum += integral;
        }
        tvals.push_back(sum);
    }
    return true;
}

bool SmoothieSIAC::v_EvaluateUsingLineAt_v1DynScaling(
    const vector<NekDouble> &stPoint, const Array<OneD, NekDouble> &direction,
    const int n_quadPts, const NekDouble meshSpacing,
    const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int varNum)
{
    tvals.clear();
    // set up line points.
    // evaluate funtion values at these points.
    NekDouble t_mesh_min = *(std::min_element(tparams.begin(), tparams.end()));
    NekDouble t_mesh_max = *(std::max_element(tparams.begin(), tparams.end()));
    vector<NekDouble> HvalT;
    vector<int> Gids, Eids;
    Array<OneD, NekDouble> t_LineElm;

    // (1) Set up the quadrature.
    SetupLineForLSIAC(direction, stPoint, t_mesh_min, t_mesh_max, n_quadPts,
                      HvalT, Gids, Eids, t_LineElm);

    // (2) Get the values of quadrature.
    Array<OneD, NekDouble> tv_LineElm(t_LineElm.num_elements());
    GetVLineForLSIAC(n_quadPts, stPoint, direction, HvalT, Gids, Eids,
                     t_LineElm, tv_LineElm, varNum);

    // (2a) Get dynscaling for quadrature.
    vector<NekDouble> t_dynScaling;
    GetDynScalingForLSIAC(stPoint, direction, tparams, t_dynScaling,
                          meshSpacing, 1.0);

    // (3) Evaluate at param locations.
    EvaluateLineForLSIAC_v1_dynScaling(n_quadPts, tparams, t_dynScaling,
                                       t_mesh_min, t_mesh_max, t_LineElm,
                                       tv_LineElm, HvalT, tvals);
    // EvaluateLineForLSIAC_v1(n_quadPts, tparams, meshSpacing, t_mesh_min,
    // t_mesh_max, 						t_LineElm, tv_LineElm,
    // HvalT, tvals
    // );

    return true;
}

bool SmoothieSIAC::v_GetDynScalingForLSIAC(
    const vector<NekDouble> &stPoint, const Array<OneD, NekDouble> &direction,
    const vector<NekDouble> &tparams, vector<NekDouble> &t_dynScaling,
    const NekDouble meshSpacing, const NekDouble mu)
{
    Array<OneD, NekDouble> sampleCoord(3, 0.0);
    for (size_t i = 0; i < tparams.size(); i++)
    {
        // Currently algorithm could be slowing down. Check if needed.
        sampleCoord[0]       = stPoint[0] + tparams[i] * direction[0];
        sampleCoord[1]       = stPoint[1] + tparams[i] * direction[1];
        sampleCoord[2]       = stPoint[2] + tparams[i] * direction[2];
        NekDouble dynScaling = m_meshHandlePtr->GetDynamicScaling(sampleCoord);
        if (meshSpacing > 0)
        {
            if (mu * dynScaling > meshSpacing)
            {
                t_dynScaling.push_back(meshSpacing);
            }
            else
            {
                t_dynScaling.push_back(mu * dynScaling);
            }
        }
        else
        {
            t_dynScaling.push_back(mu * dynScaling);
        }
    }
    return true;
}

bool SmoothieSIAC::v_EvaluateUsingLineAt(
    const vector<NekDouble> &stPoint, const Array<OneD, NekDouble> &direction,
    const int n_quadPts, const NekDouble meshSpacing,
    const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int varNum)
{
    tvals.clear();
    // set up line points.
    // evaluate funtion values at these points.
    NekDouble t_mesh_min = *(std::min_element(tparams.begin(), tparams.end()));
    NekDouble t_mesh_max = *(std::max_element(tparams.begin(), tparams.end()));
    vector<NekDouble> HvalT;
    vector<int> Gids, Eids;
    Array<OneD, NekDouble> t_LineElm;

    // (1) Set up the quadrature.
    SetupLineForLSIAC(direction, stPoint, t_mesh_min, t_mesh_max, n_quadPts,
                      HvalT, Gids, Eids, t_LineElm);

    // (2) Get the values of quadrature.
    Array<OneD, NekDouble> tv_LineElm(t_LineElm.num_elements());
    GetVLineForLSIAC(n_quadPts, stPoint, direction, HvalT, Gids, Eids,
                     t_LineElm, tv_LineElm, varNum);

    // (3) Evaluate at param locations.
    EvaluateLineForLSIAC_v1(n_quadPts, tparams, meshSpacing, t_mesh_min,
                            t_mesh_max, t_LineElm, tv_LineElm, HvalT, tvals);

    return true;
}

bool SmoothieSIAC::v_EvaluateUsingLineAt_v2DynScaling(
    const vector<NekDouble> &stPoint, const Array<OneD, NekDouble> &direction,
    const int n_quadPts, const int n_quadPts_resample,
    const NekDouble meshSpacing, const vector<NekDouble> &tparams,
    vector<NekDouble> &tvals, int varNum)
{
    tvals.clear();
    // set up line points.
    // evaluate funtion values at these points.
    NekDouble t_mesh_min = *(std::min_element(tparams.begin(), tparams.end()));
    NekDouble t_mesh_max = *(std::max_element(tparams.begin(), tparams.end()));
    vector<NekDouble> HvalT;
    vector<int> Gids, Eids;
    Array<OneD, NekDouble> t_LineElm;
    Array<OneD, NekDouble> t_LineElm_resample;
    SetupLineForLSIAC_ReSamp(direction, stPoint, t_mesh_min, t_mesh_max,
                             n_quadPts, n_quadPts_resample, HvalT, Gids, Eids,
                             t_LineElm, t_LineElm_resample);
    Array<OneD, NekDouble> tv_LineElm(t_LineElm.num_elements());
    Array<OneD, NekDouble> tv_LineElm_resample(
        t_LineElm_resample.num_elements());
    GetVLineForLSIAC(n_quadPts, stPoint, direction, HvalT, Gids, Eids,
                     t_LineElm, tv_LineElm, varNum);
    GetVLineForLSIAC_resample(n_quadPts, n_quadPts_resample, stPoint, direction,
                              HvalT, Gids, Eids, t_LineElm, tv_LineElm,
                              t_LineElm_resample, tv_LineElm_resample, varNum);

    // (2a) Get dynscaling for quadrature.
    vector<NekDouble> t_dynScaling;
    GetDynScalingForLSIAC(stPoint, direction, tparams, t_dynScaling,
                          meshSpacing, 1.0);

    // Set up the quadrature.
    // SetupLineForLSIAC( direction, stPoint, t_mesh_min, t_mesh_max,n_quadPts,
    // HvalT, Gids, Eids, t_LineElm);
    // Get the values of quadrature.
    EvaluateLineForLSIAC_v1_dynScaling(
        n_quadPts_resample, tparams, t_dynScaling, t_mesh_min, t_mesh_max,
        t_LineElm_resample, tv_LineElm_resample, HvalT, tvals);
    return true;
}
bool SmoothieSIAC::v_EvaluateUsingLineAt_v2(
    const vector<NekDouble> &stPoint, const Array<OneD, NekDouble> &direction,
    const int n_quadPts, const int n_quadPts_resample,
    const NekDouble meshSpacing, const vector<NekDouble> &tparams,
    vector<NekDouble> &tvals, int varNum)
{
    tvals.clear();
    // set up line points.
    // evaluate funtion values at these points.
    NekDouble t_mesh_min = *(std::min_element(tparams.begin(), tparams.end()));
    NekDouble t_mesh_max = *(std::max_element(tparams.begin(), tparams.end()));
    vector<NekDouble> HvalT;
    vector<int> Gids, Eids;
    Array<OneD, NekDouble> t_LineElm;
    Array<OneD, NekDouble> t_LineElm_resample;
    SetupLineForLSIAC_ReSamp(direction, stPoint, t_mesh_min, t_mesh_max,
                             n_quadPts, n_quadPts_resample, HvalT, Gids, Eids,
                             t_LineElm, t_LineElm_resample);
    Array<OneD, NekDouble> tv_LineElm(t_LineElm.num_elements());
    Array<OneD, NekDouble> tv_LineElm_resample(
        t_LineElm_resample.num_elements());
    GetVLineForLSIAC(n_quadPts, stPoint, direction, HvalT, Gids, Eids,
                     t_LineElm, tv_LineElm, varNum);
    GetVLineForLSIAC_resample(n_quadPts, n_quadPts_resample, stPoint, direction,
                              HvalT, Gids, Eids, t_LineElm, tv_LineElm,
                              t_LineElm_resample, tv_LineElm_resample, varNum);

    // Set up the quadrature.
    // SetupLineForLSIAC( direction, stPoint, t_mesh_min, t_mesh_max,n_quadPts,
    // HvalT, Gids, Eids, t_LineElm);

    // Get the values of quadrature.
    EvaluateLineForLSIAC_v1(n_quadPts_resample, tparams, meshSpacing,
                            t_mesh_min, t_mesh_max, t_LineElm_resample,
                            tv_LineElm_resample, HvalT, tvals);

    return true;
}

/// The aim of this funtion is take advantage of re-sampling and using it for
/// consistent integration.
bool SmoothieSIAC::v_EvaluateUsingLineAt_v3DynScaling(
    const vector<NekDouble> &stPoint, const Array<OneD, NekDouble> &direction,
    const int n_quadPts, const int n_quadPts_resample,
    const NekDouble meshSpacing, const vector<NekDouble> &tparams,
    vector<NekDouble> &tvals, int varNum)
{
    // set up line points.
    // evaluate funtion values at these points.
    NekDouble t_mesh_min = *(std::min_element(tparams.begin(), tparams.end()));
    NekDouble t_mesh_max = *(std::max_element(tparams.begin(), tparams.end()));
    vector<NekDouble> HvalT;
    vector<int> Gids, Eids;
    Array<OneD, NekDouble> t_LineElm;

    // Set up the quadrature.
    SetupLineForLSIAC(direction, stPoint, t_mesh_min, t_mesh_max, n_quadPts,
                      HvalT, Gids, Eids, t_LineElm);

    vector<NekDouble> t_dynScaling;
    GetDynScalingForLSIAC(stPoint, direction, tparams, t_dynScaling,
                          meshSpacing, 1.0);

    Array<OneD, NekDouble> tv_LineElm(t_LineElm.num_elements());
    // Get the values of quadrature.
    GetVLineForLSIAC(n_quadPts, stPoint, direction, HvalT, Gids, Eids,
                     t_LineElm, tv_LineElm, varNum);

    // (2a)

    EvaluateLineForLSIAC_v3_dynScaling(n_quadPts, n_quadPts_resample, tparams,
                                       t_dynScaling, t_mesh_min, t_mesh_max,
                                       t_LineElm, tv_LineElm, HvalT, tvals);
    // std::cout.precision(5);
    // cout << "part3\t" << el1 << "\t" << el2 << "\t" << el3 << "\t" << endl;
    return true;
}

bool SmoothieSIAC::v_EvaluateUsingLineAt_v3(
    const vector<NekDouble> &stPoint, const Array<OneD, NekDouble> &direction,
    const int n_quadPts, const int n_quadPts_resample,
    const NekDouble meshSpacing, const vector<NekDouble> &tparams,
    vector<NekDouble> &tvals, int varNum)
{
    // set up line points.
    // evaluate funtion values at these points.
    NekDouble t_mesh_min = *(std::min_element(tparams.begin(), tparams.end()));
    NekDouble t_mesh_max = *(std::max_element(tparams.begin(), tparams.end()));
    vector<NekDouble> HvalT;
    vector<int> Gids, Eids;
    Array<OneD, NekDouble> t_LineElm;

    // Set up the quadrature.
    SetupLineForLSIAC(direction, stPoint, t_mesh_min, t_mesh_max, n_quadPts,
                      HvalT, Gids, Eids, t_LineElm);
    Array<OneD, NekDouble> tv_LineElm(t_LineElm.num_elements());
    // Get the values of quadrature.
    GetVLineForLSIAC(n_quadPts, stPoint, direction, HvalT, Gids, Eids,
                     t_LineElm, tv_LineElm, varNum);

    EvaluateLineForLSIAC_v3(n_quadPts, n_quadPts_resample, tparams, meshSpacing,
                            t_mesh_min, t_mesh_max, t_LineElm, tv_LineElm,
                            HvalT, tvals);
    return true;
}

bool SmoothieSIAC::v_SetupLineForLSIAC(
    const Array<OneD, NekDouble> &direction, const vector<NekDouble> &stPoint,
    const NekDouble tmin, const NekDouble tmax, const int n_quadPts,
    vector<NekDouble> &HvalT, vector<int> &t_GIDs, vector<int> &t_EIDs,
    Array<OneD, NekDouble> &t_LineElm)
{
    // 1) Find all the elements on the line segment.
    // 2) Find all the break points for these elements.
    // 3) Get Elem Ids.
    // 4) Quadrature points.
    vector<NekDouble> HvalX, HvalY, HvalZ;
    m_meshHandlePtr->GetBreakPts(stPoint[0], stPoint[1], stPoint[2], direction,
                                 tmin, tmax, HvalX, HvalY, HvalZ, HvalT);
    // printNekArray(HvalT,0.0);
    m_meshHandlePtr->GetListOfGIDs(stPoint[0], stPoint[1], stPoint[2],
                                   direction, HvalT, t_GIDs, t_EIDs);
    // Number of quadrature points should be input.
    // int n_quadPts = 4;
    LibUtilities::PointsKey quadPointsKey(
        n_quadPts, Nektar::LibUtilities::eGaussGaussLegendre);
    Array<OneD, NekDouble> quad_points =
        LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    Array<OneD, NekDouble> quad_weights =
        LibUtilities::PointsManager()[quadPointsKey]->GetW();

    int n_meshT_quadPts = n_quadPts * (t_EIDs.size() - 1);
    t_LineElm           = Array<OneD, NekDouble>(n_meshT_quadPts);
    for (int i = 0; i < t_EIDs.size() - 1; i++)
    {
        for (int j = 0; j < n_quadPts; j++)
        {
            t_LineElm[i * n_quadPts + j] =
                HvalT[i] +
                (quad_points[j] + 1.0) * (HvalT[i + 1] - HvalT[i]) / 2.0;
        }
    }
    return true;
}
bool SmoothieSIAC::v_SetupLineForLSIAC_ReSamp(
    const Array<OneD, NekDouble> &direction, const vector<NekDouble> &stPoint,
    const NekDouble tmin, const NekDouble tmax, const int n_quadPts,
    const int n_quadPts_Resample, vector<NekDouble> &HvalT, vector<int> &t_GIDs,
    vector<int> &t_EIDs, Array<OneD, NekDouble> &t_LineElm,
    Array<OneD, NekDouble> &t_LineElm_Resample)
{
    // 1) Find all the elements on the line segment.
    // 2) Find all the break points for these elements.
    // 3) Get Elem Ids.
    // 4) Quadrature points.
    vector<NekDouble> HvalX, HvalY, HvalZ;
    m_meshHandlePtr->GetBreakPts(stPoint[0], stPoint[1], stPoint[2], direction,
                                 tmin, tmax, HvalX, HvalY, HvalZ, HvalT);
    // printNekArray(HvalT,0.0);
    m_meshHandlePtr->GetListOfGIDs(stPoint[0], stPoint[1], stPoint[2],
                                   direction, HvalT, t_GIDs, t_EIDs);
    // Number of quadrature points should be input.
    // int n_quadPts = 4;
    {
        LibUtilities::PointsKey quadPointsKey(
            n_quadPts, Nektar::LibUtilities::eGaussGaussLegendre);
        Array<OneD, NekDouble> quad_points =
            LibUtilities::PointsManager()[quadPointsKey]->GetZ();
        Array<OneD, NekDouble> quad_weights =
            LibUtilities::PointsManager()[quadPointsKey]->GetW();
        int n_meshT_quadPts = n_quadPts * (t_EIDs.size() - 1);
        t_LineElm           = Array<OneD, NekDouble>(n_meshT_quadPts);
        for (int i = 0; i < t_EIDs.size() - 1; i++)
        {
            for (int j = 0; j < n_quadPts; j++)
            {
                t_LineElm[i * n_quadPts + j] =
                    HvalT[i] +
                    (quad_points[j] + 1.0) * (HvalT[i + 1] - HvalT[i]) / 2.0;
            }
        }
    }
    {
        // int n_quadPts_Resample = 4;
        LibUtilities::PointsKey quadPointsKey(
            n_quadPts_Resample, Nektar::LibUtilities::eGaussGaussLegendre);
        Array<OneD, NekDouble> quad_points =
            LibUtilities::PointsManager()[quadPointsKey]->GetZ();
        Array<OneD, NekDouble> quad_weights =
            LibUtilities::PointsManager()[quadPointsKey]->GetW();
        int n_meshT_quadPts_Resample = n_quadPts_Resample * (t_EIDs.size() - 1);
        t_LineElm_Resample = Array<OneD, NekDouble>(n_meshT_quadPts_Resample);
        for (int i = 0; i < t_EIDs.size() - 1; i++)
        {
            for (int j = 0; j < n_quadPts_Resample; j++)
            {
                t_LineElm_Resample[i * n_quadPts_Resample + j] =
                    HvalT[i] +
                    (quad_points[j] + 1.0) * (HvalT[i + 1] - HvalT[i]) / 2.0;
            }
        }
    }

    return true;
}

bool SmoothieSIAC::v_GetVLineForLSIAC(
    const int n_quadPts, const vector<NekDouble> &stPoint,
    const Array<OneD, NekDouble> &direction, const vector<NekDouble> &HvalT,
    const vector<int> &t_GIDs, const vector<int> &t_EIDs,
    const Array<OneD, NekDouble> &t_LineElm, Array<OneD, NekDouble> tv_LineElm,
    int varNum)
{
    boost::ignore_unused(HvalT); // Unsed variable
    // loop through elements.
    // Find values for each quadrature point.
    // create tx_LineElm,ty_LineElm, tz_LineElm);
    Array<OneD, NekDouble> tx_LineElm(n_quadPts * (t_EIDs.size() - 1));
    Array<OneD, NekDouble> ty_LineElm(n_quadPts * (t_EIDs.size() - 1));
    Array<OneD, NekDouble> tz_LineElm(n_quadPts * (t_EIDs.size() - 1));

    for (int i = 0; i < t_LineElm.num_elements(); i++)
    {
        tx_LineElm[i] = stPoint[0] + direction[0] * t_LineElm[i];
        ty_LineElm[i] = stPoint[1] + direction[1] * t_LineElm[i];
        tz_LineElm[i] = stPoint[2] + direction[2] * t_LineElm[i];
    }

    Array<OneD, NekDouble> elx, ely, elz;
    Array<OneD, NekDouble> elv(n_quadPts, 0.0);
    Array<OneD, NekDouble> glCoord(3);
    for (int i = 0; i < t_EIDs.size() - 1; i++)
    {
        elx = tx_LineElm.CreateWithOffset(tx_LineElm, i * n_quadPts);
        ely = ty_LineElm.CreateWithOffset(ty_LineElm, i * n_quadPts);
        elz = tz_LineElm.CreateWithOffset(tz_LineElm, i * n_quadPts);
        //	elv = tz_LineElm.CreateWithOffset( tv_LineElm, i*n_quadPts);
        m_meshHandlePtr->EvaluateAt(elx, ely, elz, t_GIDs[i], t_EIDs[i], elv,
                                    varNum);
        memcpy(&tv_LineElm[i * n_quadPts], &elv[0],
               sizeof(NekDouble) * n_quadPts);
    }
    return true;
}

bool SmoothieSIAC::v_GetVLineForLSIAC_resample(
    const int n_quadPts, const int n_quadPts_resample,
    const vector<NekDouble> &stPoint, const Array<OneD, NekDouble> &direction,
    const vector<NekDouble> &HvalT, const vector<int> &t_GIDs,
    const vector<int> &t_EIDs, const Array<OneD, NekDouble> &t_LineElm,
    const Array<OneD, NekDouble> tv_LineElm,
    const Array<OneD, NekDouble> &t_LineElm_resample,
    Array<OneD, NekDouble> tv_LineElm_resample, int varNum)
{
    boost::ignore_unused(HvalT, t_GIDs, varNum); // Did not need it.
    // loop through elements.
    // Find values for each quadrature point.
    // create tx_LineElm,ty_LineElm, tz_LineElm);
    Array<OneD, NekDouble> tx_LineElm_resample(n_quadPts * (t_EIDs.size() - 1));
    Array<OneD, NekDouble> ty_LineElm_resample(n_quadPts * (t_EIDs.size() - 1));
    Array<OneD, NekDouble> tz_LineElm_resample(n_quadPts * (t_EIDs.size() - 1));

    for (int i = 0; i < t_LineElm.num_elements(); i++)
    {
        tx_LineElm_resample[i] =
            stPoint[0] + direction[0] * t_LineElm_resample[i];
        ty_LineElm_resample[i] =
            stPoint[1] + direction[1] * t_LineElm_resample[i];
        tz_LineElm_resample[i] =
            stPoint[2] + direction[2] * t_LineElm_resample[i];
    }

    Array<OneD, NekDouble> elx, ely, elz;
    Array<OneD, NekDouble> elv(n_quadPts, 0.0);
    Array<OneD, NekDouble> elv_resample(n_quadPts_resample, 0.0);
    Array<OneD, NekDouble> glCoord(3);

    // Create a STD element.
    // create Basiskey
    // m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetPointsType(0);
    // Array<OneD,NekDouble> quad_points =
    // LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    // Array<OneD,NekDouble> quad_weights  =
    // LibUtilities::PointsManager()[quadPointsKey]->GetW();
    // LibUtilities::PointsKey pk = LibUtilities::PointsKey( n_quadpts_resample,
    //				m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetPointsType(0)
    //)
    //;
    LibUtilities::PointsKey quadPointsKey(
        n_quadPts, Nektar::LibUtilities::eGaussGaussLegendre);
    LibUtilities::BasisKey bk = LibUtilities::BasisKey(
        Nektar::LibUtilities::eGauss_Lagrange, n_quadPts, quadPointsKey);
    Nektar::StdRegions::StdExpansionSharedPtr segExp_std =
        MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(bk);

    LibUtilities::PointsKey quadPointsKey_resample(
        n_quadPts_resample, Nektar::LibUtilities::eGaussGaussLegendre);
    Array<OneD, NekDouble> quad_points_resample =
        LibUtilities::PointsManager()[quadPointsKey_resample]->GetZ();
    Array<OneD, NekDouble> quad_weights_resample =
        LibUtilities::PointsManager()[quadPointsKey_resample]->GetW();

    // BasisKey
    // StdSegExp ();
    // Create new quandrature points.
    // 1. First for loop for each element.
    // 2. Second for loop to evaluate at new quadrature points.

    // cout << "test1" << segExp_std->GetNcoeffs() << "ptstest2" <<
    // segExp_std->GetTotPoints() << endl; return true;

    Array<OneD, NekDouble> Lcoord(3, 0.0);
    for (int i = 0; i < t_EIDs.size() - 1; i++)
    {
        // for each element get the values and evaluate at more points.

        //	elx = tx_LineElm.CreateWithOffset( tx_LineElm, i*n_quadPts);
        elv = tv_LineElm.CreateWithOffset(tv_LineElm, i * n_quadPts);
        // memcpy( &elv[0], &tv_LineElm[i*n_quadPts],
        // sizeof(NekDouble)*n_quadPts);
        for (int j = 0; j < n_quadPts_resample; j++)
        {
            Lcoord[0]       = quad_points_resample[j];
            elv_resample[j] = segExp_std->PhysEvaluate(Lcoord, elv);
        }
        memcpy(&tv_LineElm_resample[i * n_quadPts_resample], &elv_resample[0],
               sizeof(NekDouble) * n_quadPts_resample);
    }

    return true;
}

bool SmoothieSIAC::v_EvaluatePt_vNonSymKnots(
    const vector<NekDouble> &stPoint, const Array<OneD, NekDouble> &direction,
    const int n_quadPts, const NekDouble meshSpacing, vector<NekDouble> &tvals,
    int varNum)
{
    tvals.clear();

    // Parameter tparams does not exist.
    // meshSpacing = smallest element size;
    NekDouble tmin, tmax;
    m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, tmin, tmax);
    Array<OneD, NekDouble> knotVec(3 * (m_order - 1) + 2, 0.0);
    int NumHalfKnots = ceil((3.0 * (m_order - 1.0) + 2.0) / 2.0);

    NekDouble tTemp = 0.0;
    int numElOnL    = 0;
    bool loopCondL  = m_meshHandlePtr->WhatIsTRange(
        stPoint[0], stPoint[1], stPoint[2], direction, tmin, tTemp, numElOnL);
    while (NumHalfKnots > numElOnL && loopCondL)
    {
        tmin = tmin * 2.0;
        loopCondL =
            m_meshHandlePtr->WhatIsTRange(stPoint[0], stPoint[1], stPoint[2],
                                          direction, tmin, tTemp, numElOnL);
    }

    tTemp          = 0.0;
    int numElOnR   = 0;
    bool loopCondR = m_meshHandlePtr->WhatIsTRange(
        stPoint[0], stPoint[1], stPoint[2], direction, tTemp, tmax, numElOnR);
    while (NumHalfKnots > numElOnR && loopCondR)
    {
        tmax = tmax * 2.0;
        loopCondR =
            m_meshHandlePtr->WhatIsTRange(stPoint[0], stPoint[1], stPoint[2],
                                          direction, tTemp, tmax, numElOnR);
    }
    if (NumHalfKnots > numElOnR || NumHalfKnots > numElOnL)
    {
        tvals.push_back(-1);
        return false;
    }

    NekDouble t_mesh_min = tmin;
    NekDouble t_mesh_max = tmax;
    vector<NekDouble> tparams;
    tparams.push_back(0.0);

    // NekDouble t_mesh_min = *(std::min_element( tparams.begin(),
    // tparams.end())); NekDouble t_mesh_max = *(std::max_element(
    // tparams.begin(), tparams.end()));

    vector<NekDouble> HvalT;
    vector<int> Gids, Eids;
    Array<OneD, NekDouble> t_LineElm;

    SetupLineForLSIAC(direction, stPoint, t_mesh_min, t_mesh_max, n_quadPts,
                      HvalT, Gids, Eids, t_LineElm);

    Array<OneD, NekDouble> tv_LineElm(t_LineElm.num_elements());
    // Get the values of quadrature.
    GetVLineForLSIAC(n_quadPts, stPoint, direction, HvalT, Gids, Eids,
                     t_LineElm, tv_LineElm, varNum);

    const int n_quadPts_resample = ceil((n_quadPts + m_order) / 2.0);

    LibUtilities::PointsKey quadPointsKey(
        n_quadPts, Nektar::LibUtilities::eGaussGaussLegendre);
    LibUtilities::BasisKey bk = LibUtilities::BasisKey(
        Nektar::LibUtilities::eGauss_Lagrange, n_quadPts, quadPointsKey);
    Nektar::StdRegions::StdExpansionSharedPtr segExp_std =
        MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(bk);

    LibUtilities::PointsKey quadPointsKey_resample(
        n_quadPts_resample, Nektar::LibUtilities::eGaussGaussLegendre);
    Array<OneD, NekDouble> q_quadR_weights =
        LibUtilities::PointsManager()[quadPointsKey_resample]->GetW();
    Array<OneD, NekDouble> q_quadR_Pts =
        LibUtilities::PointsManager()[quadPointsKey_resample]->GetZ();
    Array<OneD, NekDouble> t_quadPts(n_quadPts_resample),
        t_quad_vals(n_quadPts_resample);
    //	Array<OneD,NekDouble> t_quad(n_quadPts_resample);

    Array<OneD, NekDouble> lcoord(3, 0.0);

    for (int ii = 0; ii < tparams.size(); ii++)
    {
        NekDouble t_min, t_max;
        NekDouble t = tparams[ii];
        Array<OneD, NekDouble> knotVec(3 * (m_order - 1) + 2, 0.0);
        NekDouble nonSymShift = 0.0;
        bool b_symMesh = CalculateKnotVec(t, HvalT, knotVec, nonSymShift);

        vector<NekDouble> SvalT, LvalT, TvalT;
        NekDouble sum = 0.0;

        if (b_symMesh)
        {
            SvalT.insert(SvalT.end(), &knotVec[0],
                         &knotVec[3 * (m_order - 1) + 2]);
            t_min = SvalT.front();
            t_max = SvalT.back();
        }
        else
        {
            tvals.push_back(-1);
            continue;
        }

        int startElmIndex, endElmIndex;
        const auto low =
            std::lower_bound(HvalT.begin(), HvalT.end(), t + t_min + TOLERENCE);
        const auto up =
            std::upper_bound(HvalT.begin(), HvalT.end(), t + t_max - TOLERENCE);
        // int stIndex = low - HvalT.begin()-1;
        // int edIndex = up -HvalT.begin()-1;
        int stIndex   = std::distance(HvalT.begin(), low);
        int edIndex   = std::distance(HvalT.begin(), up);
        startElmIndex = stIndex;
        endElmIndex   = edIndex;

        // Insert all mesh breakpoints including min and max using
        // HvalT into LvalT
        LvalT.clear();
        LvalT.push_back(t_min);
        // for ( int j= startElmIndex; j <= endElmIndex+1; j++)
        for (int j = startElmIndex; j < endElmIndex; j++)
        {
            if ((HvalT[j] > t + t_min) && (HvalT[j] < t + t_max))
            {
                LvalT.push_back(HvalT[j] - t);
            }
        }
        LvalT.push_back(t_max);
        mergeBreakPts(SvalT, LvalT, TvalT);

        // m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, diretion, TvalT,
        // t_GIDs,t_EIDs); Need to be done.
        for (int i = 0; i < TvalT.size() - 1; i++)
        {
            NekDouble a = TvalT[i];
            NekDouble b = TvalT[i + 1];
            //	int gID = t_GIDs[i]; need elmIndex on line. Hence code below.
            //	int eID = t_EIDs[i];
            int elmIndex = -1;
            for (int j = startElmIndex; j <= endElmIndex + 1; j++)
            {
                if (HvalT[j] > t + (a + b) / 2.0)
                {
                    elmIndex = j - 1;
                    break;
                }
            }
            assert(elmIndex != -1 && "Elm Index is wrong");
            for (int j = 0; j < n_quadPts_resample; j++)
            {
                t_quadPts[j] = (b - a) * (q_quadR_Pts[j] + 1.0) / 2.0 + a;
            }
            if (b_symMesh)
            {
                m_siacFilterPtrs[m_SymID]->EvaluateFilterWknots(
                    t_quadPts, t_quad_vals, knotVec, 1.0);
            }
            else
            {
                // m_siacFilterPtrs[m_OneID]->EvaluateFilter(t_quadPts,t_quad_vals,
                // meshSpacing, meshTShift, true);
                assert(false && "Should no be here.");
            }
            //	m_meshHandlePtr->EvaluateAt(t_xyz_vals);

            NekDouble integral = 0.0;
            Array<OneD, NekDouble> elv =
                tv_LineElm.CreateWithOffset(tv_LineElm, elmIndex * n_quadPts);
            for (int k = 0; k < n_quadPts_resample; k++)
            {
                // need to calcuate value of t_xyz_vals.
                // scaling of HvalT[j-1] to -1 and HvalT[j] to +1
                //	lcoord[0] = (t+ t_quadPts[k] - HvalT[elmIndex-1])/
                //(HvalT[elmIndex]-HvalT[elmIndex-1])*2.0 -1.0;
                lcoord[0] = (t + t_quadPts[k] - HvalT[elmIndex]) /
                                (HvalT[elmIndex + 1] - HvalT[elmIndex]) * 2.0 -
                            1.0;
                NekDouble t_xyz_val = segExp_std->PhysEvaluate(lcoord, elv);
                // SegExpasionPhysEvaluate
                integral += 0.5 * q_quadR_weights[k] * t_quad_vals[k] *
                            std::abs(b - a) * t_xyz_val;
            }
            sum += integral;
        }
        tvals.push_back(sum);
    }
    return true;
}

bool SmoothieSIAC::v_EvaluateUsingLineAt_vNonSymKnots(
    const vector<NekDouble> &stPoint, const Array<OneD, NekDouble> &direction,
    const int n_quadPts, const NekDouble meshSpacing,
    const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int varNum)
{
    boost::ignore_unused(meshSpacing);
    tvals.clear();
    NekDouble t_mesh_min = *(std::min_element(tparams.begin(), tparams.end()));
    NekDouble t_mesh_max = *(std::max_element(tparams.begin(), tparams.end()));
    vector<NekDouble> HvalT;
    vector<int> Gids, Eids;
    Array<OneD, NekDouble> t_LineElm;

    SetupLineForLSIAC(direction, stPoint, t_mesh_min, t_mesh_max, n_quadPts,
                      HvalT, Gids, Eids, t_LineElm);

    Array<OneD, NekDouble> tv_LineElm(t_LineElm.num_elements());
    // Get the values of quadrature.
    GetVLineForLSIAC(n_quadPts, stPoint, direction, HvalT, Gids, Eids,
                     t_LineElm, tv_LineElm, varNum);

    const int n_quadPts_resample = ceil((n_quadPts + m_order) / 2.0);

    LibUtilities::PointsKey quadPointsKey(
        n_quadPts, Nektar::LibUtilities::eGaussGaussLegendre);
    LibUtilities::BasisKey bk = LibUtilities::BasisKey(
        Nektar::LibUtilities::eGauss_Lagrange, n_quadPts, quadPointsKey);
    Nektar::StdRegions::StdExpansionSharedPtr segExp_std =
        MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(bk);

    LibUtilities::PointsKey quadPointsKey_resample(
        n_quadPts_resample, Nektar::LibUtilities::eGaussGaussLegendre);
    Array<OneD, NekDouble> q_quadR_weights =
        LibUtilities::PointsManager()[quadPointsKey_resample]->GetW();
    Array<OneD, NekDouble> q_quadR_Pts =
        LibUtilities::PointsManager()[quadPointsKey_resample]->GetZ();
    Array<OneD, NekDouble> t_quadPts(n_quadPts_resample),
        t_quad_vals(n_quadPts_resample);
    //	Array<OneD,NekDouble> t_quad(n_quadPts_resample);

    Array<OneD, NekDouble> lcoord(3, 0.0);

    for (int ii = 0; ii < tparams.size(); ii++)
    {
        NekDouble t_min, t_max;
        NekDouble t = tparams[ii];
        Array<OneD, NekDouble> knotVec(3 * (m_order - 1) + 2, 0.0);
        NekDouble nonSymShift = 0.0;
        bool b_symMesh = CalculateKnotVec(t, HvalT, knotVec, nonSymShift);

        vector<NekDouble> SvalT, LvalT, TvalT;
        NekDouble sum = 0.0;

        if (b_symMesh)
        {
            // SvalT.resize(knotVec.num_elements());
            SvalT.insert(SvalT.end(), &knotVec[0],
                         &knotVec[3 * (m_order - 1) + 2]);
            t_min = SvalT.front();
            t_max = SvalT.back();
        }
        else
        {
            tvals.push_back(-1);
            continue;
        }

        int startElmIndex, endElmIndex;
        const auto low =
            std::lower_bound(HvalT.begin(), HvalT.end(), t + t_min + TOLERENCE);
        const auto up =
            std::upper_bound(HvalT.begin(), HvalT.end(), t + t_max - TOLERENCE);
        // int stIndex = low - HvalT.begin()-1;
        // int edIndex = up -HvalT.begin()-1;
        int stIndex   = std::distance(HvalT.begin(), low);
        int edIndex   = std::distance(HvalT.begin(), up);
        startElmIndex = stIndex;
        endElmIndex   = edIndex;

        // Insert all mesh breakpoints including min and max using
        // HvalT into LvalT
        LvalT.clear();
        LvalT.push_back(t_min);
        // for ( int j= startElmIndex; j <= endElmIndex+1; j++)
        for (int j = startElmIndex; j < endElmIndex; j++)
        {
            if ((HvalT[j] > t + t_min) && (HvalT[j] < t + t_max))
            {
                LvalT.push_back(HvalT[j] - t);
            }
        }
        LvalT.push_back(t_max);
        mergeBreakPts(SvalT, LvalT, TvalT);

        // m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, diretion, TvalT,
        // t_GIDs,t_EIDs); Need to be done.
        for (int i = 0; i < TvalT.size() - 1; i++)
        {
            NekDouble a = TvalT[i];
            NekDouble b = TvalT[i + 1];
            //	int gID = t_GIDs[i]; need elmIndex on line. Hence code below.
            //	int eID = t_EIDs[i];
            int elmIndex = -1;
            for (int j = startElmIndex; j <= endElmIndex + 1; j++)
            {
                if (HvalT[j] > t + (a + b) / 2.0)
                {
                    elmIndex = j - 1;
                    break;
                }
            }
            assert(elmIndex != -1 && "Elm Index is wrong");
            for (int j = 0; j < n_quadPts_resample; j++)
            {
                t_quadPts[j] = (b - a) * (q_quadR_Pts[j] + 1.0) / 2.0 + a;
            }
            if (b_symMesh)
            {
                m_siacFilterPtrs[m_SymID]->EvaluateFilterWknots(
                    t_quadPts, t_quad_vals, knotVec, 1.0);
            }
            else
            {
                // m_siacFilterPtrs[m_OneID]->EvaluateFilter(t_quadPts,t_quad_vals,
                // meshSpacing, meshTShift, true);
                assert(false && "Should no be here.");
            }
            //	m_meshHandlePtr->EvaluateAt(t_xyz_vals);

            NekDouble integral = 0.0;
            Array<OneD, NekDouble> elv =
                tv_LineElm.CreateWithOffset(tv_LineElm, elmIndex * n_quadPts);
            for (int k = 0; k < n_quadPts_resample; k++)
            {
                // need to calcuate value of t_xyz_vals.
                // scaling of HvalT[j-1] to -1 and HvalT[j] to +1
                //	lcoord[0] = (t+ t_quadPts[k] - HvalT[elmIndex-1])/
                //(HvalT[elmIndex]-HvalT[elmIndex-1])*2.0 -1.0;
                lcoord[0] = (t + t_quadPts[k] - HvalT[elmIndex]) /
                                (HvalT[elmIndex + 1] - HvalT[elmIndex]) * 2.0 -
                            1.0;
                NekDouble t_xyz_val = segExp_std->PhysEvaluate(lcoord, elv);
                // SegExpasionPhysEvaluate
                integral += 0.5 * q_quadR_weights[k] * t_quad_vals[k] *
                            std::abs(b - a) * t_xyz_val;
            }
            sum += integral;
        }
        tvals.push_back(sum);
    }
    return true;
}

bool SmoothieSIAC::CalculateKnotVec(NekDouble t, vector<NekDouble> &HvalT,
                                    Array<OneD, NekDouble> &knotVec,
                                    NekDouble &shift)
{
    int degree = m_order - 1;
    auto low   = std::lower_bound(HvalT.begin(), HvalT.end(), t);
    auto up    = std::upper_bound(HvalT.begin(), HvalT.end(), t);
    if (low == up)
    {
        low--;
    }
    int distanceFromStart = std::distance(HvalT.begin(), low);
    int distanceToEnd     = std::distance(low, HvalT.end());

    int cknotnum, fknotnum;
    if (degree % 2 == 1)
    {
        cknotnum = ceil((3.0 * degree + 2.0) / 2.0);
        fknotnum = floor((3.0 * degree + 2.0) / 2.0);
    }
    else
    {
        cknotnum = ((3.0 * degree + 2.0) / 2.0);
        fknotnum = ((3.0 * degree + 2.0) / 2.0);
    }

    if (distanceFromStart < cknotnum - 1 || distanceToEnd < fknotnum + 2)
    {
        return false;
    }

    auto startIt = low;
    startIt      = low - cknotnum + 1;

    if (degree % 2 == 1)
    {
        // Check if with symmetric Boundary.
        if (distanceFromStart < cknotnum - 1 || distanceToEnd < fknotnum + 2)
        {
            return false;
        }

        /*
                        // Contant in each line segment.
                NekDouble diff = *low;
                for(int i =0; i< 3*degree+2; i++)
                {
                    //knotVec[i] = *(startIt+i)-diff;
                    knotVec[i] = *(startIt+i);
                }
                shift= t-diff;
        */
        // variable with each line segment.
        NekDouble ratio = (t - *low) / (*(low + 1) - *low);
        for (int i = 0; i < 3 * degree + 2; i++)
        {
            // knotVec[i] = *(startIt+i)-diff;
            NekDouble frontknot = *(startIt + i);
            NekDouble backknot  = *(startIt + i + 1);
            knotVec[i] = frontknot + ratio * (backknot - frontknot) - t;
        }
        shift = t - *low;
        // knotVec[2] = ratio; // Debug purpose
    }
    else
    {
        // Check if with symmetric Boundary.
        if (distanceFromStart < cknotnum || distanceToEnd < fknotnum + 2)
        {
            return false;
        }
        /*
                NekDouble diff = (*low+*(low+1))/2.0;
                for(int i =0; i< 3*degree+2; i++)
                        {
                    //knotVec[i] = *(startIt+i)-diff;
                    knotVec[i] = *(startIt+i);
                }
                shift= t-diff;
        */
        NekDouble ratio = (t - *low) / (*(low + 1) - *low);
        if (ratio > 0.5)
        {
            for (int i = 0; i < 3 * degree + 2; i++)
            {
                // knotVec[i] = *(startIt+i)-diff;
                NekDouble frontknot = *(startIt + i);
                NekDouble backknot  = *(startIt + i + 1);
                knotVec[i] =
                    frontknot + (ratio - 0.5) * (backknot - frontknot) - t;
            }
        }
        else
        {
            for (int i = 0; i < 3 * degree + 2; i++)
            {
                // knotVec[i] = *(startIt+i)-diff;
                NekDouble frontknot = *(startIt + i - 1);
                NekDouble backknot  = *(startIt + i);
                knotVec[i] =
                    frontknot + (0.5 + ratio) * (backknot - frontknot) - t;
            }
        }
        shift = t - *low;
        // knotVec[2] = ratio; // Debug purpose
    }

    return true;
}

} // namespace LSIAC
} // namespace Nektar
