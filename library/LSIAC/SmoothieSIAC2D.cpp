#include "SmoothieSIAC2D.h"
#include "OneSidedSIAC.h"
#include <LibUtilities/Foundations/ManagerAccess.h> // for Points Manager, etc
#include <algorithm>
#include <boost/timer.hpp>
#include <numeric>

namespace Nektar
{
namespace LSIAC
{
SmoothieSIAC2D::SmoothieSIAC2D(const FilterType filter,
                               HandleNekMesh *meshHandle, const int order,
                               NekDouble meshSpacing, const int derivative)
    : SmoothieSIAC(filter, order), m_meshSpacing(meshSpacing)
{
    // Need to Create SIAC filter depending on filterType.
    // Set the filterType flags to be used in the code.
    m_meshHandlePtr = meshHandle;
    switch (filter)
    {
        case eNONE:
            m_OneID = -1;
            m_SymID = -1;
            break;
        case eSYM_2kp1:
            // Set up symmetric siac filter and any other parameters needed for
            // case scenario. m_symSIACptr( new SymmetricSIAC(m_order));
            // m_siacFilterPtrs.push_back( SIACFilter() );
            //			cout << "Into SmoothieSIAC1D constructor working
            //:)
            //???
            //"<<m_order << endl;
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(order));
            m_OneID = -1;
            break;
        case eSYM_2kp1_1SIDED_2kp1:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(order));
            m_siacFilterPtrs.emplace_back(new OneSidedSIAC(
                order, OneSidedSIAC::OneSidedFilterType::BASIC_SIAC_2kp1));
            break;
        case eSYM_2kp1_1SIDED_4kp1:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(order));
            m_siacFilterPtrs.emplace_back(new OneSidedSIAC(
                order, OneSidedSIAC::OneSidedFilterType::VAN_SIAC_4kp1));
            break;
        case eSYM_4kp1:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(order, 4 * order));
            m_OneID = -1;
            assert(false && "symmetric 4k+1 filter is some how screwed up.");
            break;
        case eSYM_DER_2kp1_1SIDED_2kp1:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(
                order,
                SymmetricSIAC::SymFilterType::CUSTOM_SMOOTH_Derivative_SIAC,
                derivative));
            m_siacFilterPtrs.emplace_back(new OneSidedSIAC(
                order,
                OneSidedSIAC::OneSidedFilterType::Der_SMOOTH_BASIC_SIAC_2kp1,
                derivative));
            // m_OneID = -1; // Need to be replaced by derivative filter.
            break;
        case eSYM_DER_2kp1_1SIDED_4kp1:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(
                order,
                SymmetricSIAC::SymFilterType::CUSTOM_SMOOTH_Derivative_SIAC,
                derivative));
            m_siacFilterPtrs.emplace_back(new OneSidedSIAC(
                order, OneSidedSIAC::OneSidedFilterType::Der_BASIC_SIAC_4kp1,
                derivative));
            // m_OneID = -1; // Need to be replaced by derivative filter.
            break;
        case eSYM_2kp1_1SIDED_2kp2:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(order));
            m_siacFilterPtrs.emplace_back(new OneSidedSIAC(
                order, OneSidedSIAC::OneSidedFilterType::XLi_SIAC_2kp2));
            break;
        case eSYM_DER_2kp1_1SIDED_2kp2:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(
                order,
                SymmetricSIAC::SymFilterType::CUSTOM_SMOOTH_Derivative_SIAC,
                derivative));
            m_siacFilterPtrs.emplace_back(new OneSidedSIAC(
                order, OneSidedSIAC::OneSidedFilterType::Der_XLi_SIAC_2kp2,
                derivative));
            break;
        case eSYM_NDER_2kp1_1SIDED_2kp1:
            m_siacFilterPtrs.emplace_back(new SymmetricSIAC(
                order,
                SymmetricSIAC::SymFilterType::
                    CUSTOM_SMOOTH_Derivative_SIAC_WOUT_DIVDIFF,
                derivative));
            m_siacFilterPtrs.emplace_back(new OneSidedSIAC(
                order,
                OneSidedSIAC::OneSidedFilterType::N_Der_SMOOTH_BASIC_SIAC_2kp1,
                derivative));
            break;
        case eSYM_UNEVEN_2kp1:
            m_siacFilterPtrs.emplace_back(new NonSymmetricSIAC(order));
            m_OneID = -1;
            break;
        default:
            cout << "Not implemented yet " << endl;
    }
}

bool SmoothieSIAC2D::v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                                  const NekDouble PtsZ, NekDouble &valX,
                                  NekDouble &valY, NekDouble &valZ)
{
    Array<OneD, NekDouble> direction(3, 0.0);
    direction[0] = 1.0;
    return EvaluateAt(PtsX, PtsY, PtsZ, valX, valY, valZ, direction,
                      m_meshSpacing);
}
/*
bool SmoothieSIAC2D::v_EvaluateAt_SYM(const NekDouble PtsX, const NekDouble
PtsY, const NekDouble PtsZ, NekDouble &valX, NekDouble &valY, NekDouble &valZ,
Array<OneD,NekDouble> &direction, NekDouble meshSpacing, int varNum)
{
                // call HNM
                // TO do . Direction needs to picked up from meshptr. ???
        if (meshSpacing <0)
        { // No parameter was specified.
                meshSpacing = m_meshSpacing;
        }
        NekDouble meshTShift=0.0;
        NekDouble tmin,tmax;
        vector<NekDouble> HvalX,HvalY,HvalZ,HvalT,SvalT,TvalT;
                // HNM to get list of break points.
                // check if a symmetric fitler can be applied.
                // in position, symmetric range

                // Get filter range given scaling.
        m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing,tmin,tmax);

// break filter into two parts.
        // given a point direction tmin and tmax.
        //


        bool b_symMesh =
m_meshHandlePtr->CanTRangebeApplied(PtsX,PtsY,PtsZ,direction,tmin,tmax,meshTShift);
                // Made an assumtion that a symmetric filter can be applied ???
                // return List of mesh breakpoints in that range.
                // call SSS/SSO
        if (b_symMesh)
        {
                //	cout << "into symmetric condition loop" << endl;
                m_siacFilterPtrs[m_SymID]->GetBreakPts(meshSpacing , SvalT);
//"??? specify parameters" }else{
                        //cout << "PtsX: "<< PtsX << " before tmin : " << tmin;
                        //cout << " tmax : " << tmax << endl;
                        //cout << "after tmin : " << tmin;
                        //cout << " tmax : " << tmax << endl;
                SvalT.clear();
                if (m_OneID >= 0)
                { 	// OneSided Filter is defined and given by .
                        m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing,tmin,tmax);
                        m_meshHandlePtr->CanTRangebeApplied(PtsX,PtsY,PtsZ,direction,tmin,tmax,meshTShift);
                        // Next step could be equal to tmin = tmin+meshTShift;
tmax = tmax+meshTShift;
                        m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing,tmin,tmax,meshTShift);
                        m_siacFilterPtrs[m_OneID]->GetBreakPts(meshSpacing ,
SvalT,meshTShift); //"??? specify parameters"
                        //printNekArray( SvalT, 0);
                //	cout << "ptsX: " << PtsX ;
                //	cout << "\tmeshTShift: " << meshTShift << endl;
                }else
                { 	// OneSided Filter is not defined.
                        cout << "Symmetric filter does not fit and OneSided is
not defined here." << endl; cout << "Ideally keep the original value. -1 is
returned at end " << endl; valX = -1; return false;
                }

        }

        m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ,direction, tmin,tmax,
HvalX, HvalY, HvalZ, HvalT);
        // del vector<int> t_GIDs,t_EIDs;
        // del m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, direction, TvalT,
t_GIDs, t_EIDs); vector<int> t_GIDs,t_EIDs; vector<int> t_h_GIDs, t_h_EIDs;
        m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, HvalT,
t_h_GIDs, t_h_EIDs); mergeBreakPtsAdv(HvalT,SvalT, TvalT, t_h_GIDs, t_h_EIDs,
t_GIDs, t_EIDs);
        //mergeBreakPts( HvalT, SvalT, TvalT);

                // BPTS
        //printNekArray(HvalT,0);
        //printNekArray(SvalT,0);
        //printNekArray(TvalT,0);

        // Loop forsumation.

        NekDouble sum=0.0;
        // For getting the local gauss quadrature use the getCoords fo
StdRegions objects.

                                //cout << "////////Working till here////////"
<<endl; int quadratureOfMesh =
m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(0) +
                                m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(1);
                        //
m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(2);

        int quad_npoints = ceil( (m_order + quadratureOfMesh +1.0)/2.0);
        //int quad_npoints = ceil((m_order+3)/2) +
m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetTotPoints()+2;
        //cout << "quad_npoints" << quad_npoints << endl;
        LibUtilities::PointsKey quadPointsKey(quad_npoints,
                                                m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetPointsType(0));

        Array<OneD,NekDouble> quad_points
                                                        =
LibUtilities::PointsManager()[quadPointsKey]->GetZ(); Array<OneD,NekDouble>
quad_weights = LibUtilities::PointsManager()[quadPointsKey]->GetW();
        Array<OneD,NekDouble> t_quad(quad_npoints), t_x(quad_npoints),
                                                t_y(quad_npoints),t_z(quad_npoints),
                                                t_quad_vals(quad_npoints),
t_xyz_vals(quad_npoints);
//	vector<int> t_GIDs,t_EIDs;
//	m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, direction, TvalT, t_GIDs,
t_EIDs);

        for (int i =0; i < TvalT.size()-1; i++)
        {
                // Create a seg element
                        // Actual quadrature points needed are for order =
polyorder+Bspline order. NekDouble a = TvalT[i]; NekDouble b = TvalT[i+1]; int
gID = t_GIDs[i]; int eID = t_EIDs[i]; for (int j=0; j< quad_npoints; j++)
                {
                        t_quad[j] = (b-a)*(quad_points[j]+1.0)/2.0 +a;
                        t_x[j] = PtsX+t_quad[j]*direction[0];
                        t_y[j] = PtsY+t_quad[j]*direction[1];
                        t_z[j] = PtsZ+t_quad[j]*direction[2];
                }
                // Evaluate filter  at these points.
                if (b_symMesh)
                {
                        m_siacFilterPtrs[m_SymID]->EvaluateFilter(t_quad,t_quad_vals,
meshSpacing ); }else{
                        m_siacFilterPtrs[m_OneID]->EvaluateFilter(t_quad,t_quad_vals,
meshSpacing, meshTShift, true);
                }
                m_meshHandlePtr->EvaluateAt( t_x, t_y,t_z,gID,eID,
t_xyz_vals,varNum); NekDouble integral = 0;
                // integral
                for( int k=0; k < quad_npoints;k++)
                {
                        //integral += 0.5*quad_weights[k] * t_quad_vals[k]
*(b-a); // *  t_xyz_vals[k]; integral += 0.5*quad_weights[k] * t_quad_vals[k]
*std::abs((b-a)) *  t_xyz_vals[k];
                }
                sum += integral;
        }
        valX = sum;
        return true;
}
*/

bool SmoothieSIAC2D::v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                                  const NekDouble PtsZ, NekDouble &valX,
                                  NekDouble &valY, NekDouble &valZ,
                                  Array<OneD, NekDouble> &direction,
                                  NekDouble meshSpacing, int varNum)
{
    boost::ignore_unused(valY, valZ); // reserved for derivative application
    // call HNM
    // TO do . Direction needs to picked up from meshptr. ???
    if (meshSpacing < 0)
    { // No parameter was specified.
        meshSpacing = m_meshSpacing;
    }
    NekDouble meshTShift = 0.0;
    NekDouble tmin, tmax;
    vector<NekDouble> HvalX, HvalY, HvalZ, HvalT, SvalT, TvalT;
    // HNM to get list of break points.
    // check if a symmetric fitler can be applied.
    // in position, symmetric range

    // Get filter range given scaling.
    m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, tmin, tmax);

    bool b_symMesh = m_meshHandlePtr->CanTRangebeApplied(
        PtsX, PtsY, PtsZ, direction, tmin, tmax, meshTShift);
    if (b_symMesh)
    {
        m_siacFilterPtrs[m_SymID]->GetBreakPts(
            meshSpacing, SvalT); //"??? specify parameters"
    }
    else
    {
        SvalT.clear();
        if (m_OneID >= 0)
        { // OneSided Filter is defined and given by .
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, tmin, tmax);
            m_meshHandlePtr->CanTRangebeApplied(PtsX, PtsY, PtsZ, direction,
                                                tmin, tmax, meshTShift);
            // Next step could be equal to tmin = tmin+meshTShift; tmax =
            // tmax+meshTShift;
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, tmin, tmax,
                                                      meshTShift);
            if (!m_meshHandlePtr->CanTRangebeAppliedWOMeshShift(
                    PtsX, PtsY, PtsZ, direction, tmin, tmax))
            {
                // Filter does not fit. Hence returning the dG value.
                m_meshHandlePtr->EvaluateAt(PtsX, PtsY, PtsZ, -1, -1, valX,
                                            varNum);
                return false;
            }
            m_siacFilterPtrs[m_OneID]->GetBreakPts(
                meshSpacing, SvalT, meshTShift); //"??? specify parameters"
        }
        else
        { // OneSided Filter is not defined.
          // valX = -1;
            m_meshHandlePtr->EvaluateAt(PtsX, PtsY, PtsZ, -1, -1, valX, varNum);
            return false;
        }
    }

    m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ, direction, tmin, tmax, HvalX,
                                 HvalY, HvalZ, HvalT);
    // del vector<int> t_GIDs,t_EIDs;
    // del m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, direction, TvalT,
    // t_GIDs, t_EIDs);
    vector<int> t_GIDs, t_EIDs;
    vector<int> t_h_GIDs, t_h_EIDs;
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, HvalT, t_h_GIDs,
                                   t_h_EIDs);
    mergeBreakPtsAdv(HvalT, SvalT, TvalT, t_h_GIDs, t_h_EIDs, t_GIDs, t_EIDs);
    // mergeBreakPts( HvalT, SvalT, TvalT);

    // BPTS
    // printNekArray(HvalT,0);
    // printNekArray(SvalT,0);
    // printNekArray(TvalT,0);

    // Loop forsumation.

    NekDouble sum = 0.0;
    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.

    // cout << "////////Working till here////////" <<endl;
    int quadratureOfMesh =
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(0) +
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(1);
    //	m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(2);

    int quad_npoints = ceil((m_order + quadratureOfMesh + 1.0) / 2.0);
    // int quad_npoints = ceil((m_order+3)/2) +
    // m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetTotPoints()+2; cout <<
    // "quad_npoints" << quad_npoints << endl;
    // LibUtilities::PointsKey quadPointsKey(quad_npoints,
    //						m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetPointsType(0));
    LibUtilities::PointsKey quadPointsKey(
        quad_npoints, LibUtilities::PointsType::eGaussGaussLegendre);

    Array<OneD, NekDouble> quad_points =
        LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    Array<OneD, NekDouble> quad_weights =
        LibUtilities::PointsManager()[quadPointsKey]->GetW();
    Array<OneD, NekDouble> t_quad(quad_npoints), t_x(quad_npoints),
        t_y(quad_npoints), t_z(quad_npoints), t_quad_vals(quad_npoints),
        t_xyz_vals(quad_npoints);
    //	vector<int> t_GIDs,t_EIDs;
    //	m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, direction, TvalT, t_GIDs,
    // t_EIDs);

    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        // Create a seg element
        // Actual quadrature points needed are for order = polyorder+Bspline
        // order.
        NekDouble a = TvalT[i];
        NekDouble b = TvalT[i + 1];
        int gID     = t_GIDs[i];
        int eID     = t_EIDs[i];
        for (int j = 0; j < quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (quad_points[j] + 1.0) / 2.0 + a;
            t_x[j]    = PtsX + t_quad[j] * direction[0];
            t_y[j]    = PtsY + t_quad[j] * direction[1];
            t_z[j]    = PtsZ + t_quad[j] * direction[2];
        }
        // Evaluate filter  at these points.
        if (b_symMesh)
        {
            m_siacFilterPtrs[m_SymID]->EvaluateFilter(t_quad, t_quad_vals,
                                                      meshSpacing);
        }
        else
        {
            m_siacFilterPtrs[m_OneID]->EvaluateFilter(
                t_quad, t_quad_vals, meshSpacing, meshTShift, true);
        }
        m_meshHandlePtr->EvaluateAt(t_x, t_y, t_z, gID, eID, t_xyz_vals,
                                    varNum);
        NekDouble integral = 0;
        // integral
        for (int k = 0; k < quad_npoints; k++)
        {
            // integral += 0.5*quad_weights[k] * t_quad_vals[k] *(b-a); // *
            // t_xyz_vals[k];
            integral += 0.5 * quad_weights[k] * t_quad_vals[k] *
                        std::abs((b - a)) * t_xyz_vals[k];
        }
        sum += integral;
    }
    valX = sum;
    return true;
}

bool SmoothieSIAC2D::v_EvaluateAt_NSK_FixedNumBSpl(
    const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
    NekDouble &valX, NekDouble &valY, NekDouble &valZ,
    Array<OneD, NekDouble> &direction, NekDouble meshSpacing, int varNum)
{
    // call HNM
    // TO do . Direction needs to picked up from meshptr. ???
    if (meshSpacing < 0)
    { // No parameter was specified.
        meshSpacing = m_meshSpacing;
    }
    NekDouble tmin, tmax;
    vector<NekDouble> HvalX, HvalY, HvalZ, HvalT, SvalT, TvalT, LvalT;
    // HNM to get list of break points.
    // check if a symmetric fitler can be applied.
    // in position, symmetric range

    // Get filter range given scaling.
    m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, tmin, tmax);

    // Figuring out the num nodes problem.
    Array<OneD, NekDouble> knotVec(3 * (m_order - 1) + 2, 0.0);
    int NumHalfKnots = ceil((3.0 * (m_order - 1.0) + 2.0) / 2.0);

    NekDouble tTemp = 0.0;
    int numElOnL    = 0;
    bool loopCondL  = m_meshHandlePtr->WhatIsTRange(PtsX, PtsY, PtsZ, direction,
                                                   tmin, tTemp, numElOnL);
    while ((NumHalfKnots >= (numElOnL - 2)) && loopCondL)
    {
        tmin      = tmin * 2.0;
        loopCondL = m_meshHandlePtr->WhatIsTRange(PtsX, PtsY, PtsZ, direction,
                                                  tmin, tTemp, numElOnL);
    }

    tTemp          = 0.0;
    int numElOnR   = 0;
    bool loopCondR = m_meshHandlePtr->WhatIsTRange(PtsX, PtsY, PtsZ, direction,
                                                   tTemp, tmax, numElOnR);
    while ((NumHalfKnots >= (numElOnR - 2)) && loopCondR)
    {
        tmax      = tmax * 2.0;
        loopCondR = m_meshHandlePtr->WhatIsTRange(PtsX, PtsY, PtsZ, direction,
                                                  tTemp, tmax, numElOnR);
    }
    if ((NumHalfKnots >= (numElOnR - 2)) || (NumHalfKnots >= (numElOnL - 2)))
    {
        valX = -1;
        return false;
    }

    // This is function which returns HvalT. HvalT is the filter knots.
    // m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ,direction, tmin,tmax,
    // HvalX, HvalY, HvalZ, HvalT);
    m_meshHandlePtr->GetBreakPts_Without_Tmin_Tmax(
        PtsX, PtsY, PtsZ, direction, tmin, tmax, HvalX, HvalY, HvalZ, HvalT);

    // set Tmin and Tmax at exact number of nodes.
    // HOW TO DO IT.
    NekDouble nonSymShift = 0.0;
    bool b_symMesh        = CalculateKnotVec(0.0, HvalT, knotVec, nonSymShift);

    if (b_symMesh)
    {
        SvalT.insert(SvalT.end(), &knotVec[0], &knotVec[3 * (m_order - 1) + 2]);
        tmin = SvalT.front();
        tmax = SvalT.back();
    }
    else
    {
        assert("Should have been returned before coming here.");
    }

    m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ, direction, tmin, tmax, HvalX,
                                 HvalY, HvalZ, LvalT);
    valY = tmax - tmin;
    valZ = knotVec.num_elements();

    vector<int> t_GIDs, t_EIDs;
    vector<int> t_h_GIDs, t_h_EIDs;
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, LvalT, t_h_GIDs,
                                   t_h_EIDs);
    mergeBreakPtsAdv(LvalT, SvalT, TvalT, t_h_GIDs, t_h_EIDs, t_GIDs, t_EIDs);

    NekDouble sum = 0.0;
    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.
    int quadratureOfMesh =
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(0) +
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(1);
    //	m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(2);

    int quad_npoints = ceil((m_order + quadratureOfMesh + 1.0) / 2.0);
    LibUtilities::PointsKey quadPointsKey(
        quad_npoints, LibUtilities::PointsType::eGaussGaussLegendre);
    Array<OneD, NekDouble> quad_points =
        LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    Array<OneD, NekDouble> quad_weights =
        LibUtilities::PointsManager()[quadPointsKey]->GetW();
    Array<OneD, NekDouble> t_quad(quad_npoints), t_x(quad_npoints),
        t_y(quad_npoints), t_z(quad_npoints), t_quad_vals(quad_npoints),
        t_xyz_vals(quad_npoints);

    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        // Create a seg element
        // Actual quadrature points needed are for order = polyorder+Bspline
        // order.
        NekDouble a = TvalT[i];
        NekDouble b = TvalT[i + 1];
        int gID     = t_GIDs[i];
        int eID     = t_EIDs[i];
        for (int j = 0; j < quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (quad_points[j] + 1.0) / 2.0 + a;
            t_x[j]    = PtsX + t_quad[j] * direction[0];
            t_y[j]    = PtsY + t_quad[j] * direction[1];
            t_z[j]    = PtsZ + t_quad[j] * direction[2];
        }
        // Evaluate filter  at these points.
        if (b_symMesh)
        {
            //	m_siacFilterPtrs[m_SymID]->EvaluateFilter(t_quad,t_quad_vals,
            // meshSpacing );
            //		m_siacFilterPtrs[m_SymID]->EvaluateFilter_GivenNumSplines(t_quad,t_quad_vals,
            // meshSpacing );
            m_siacFilterPtrs[m_SymID]->EvaluateFilterWknots(t_quad, t_quad_vals,
                                                            knotVec, 1.0);
        }
        else
        {
            //	m_siacFilterPtrs[m_OneID]->EvaluateFilter(t_quad,t_quad_vals,
            // meshSpacing, meshTShift, true);
        }
        m_meshHandlePtr->EvaluateAt(t_x, t_y, t_z, gID, eID, t_xyz_vals,
                                    varNum);
        NekDouble integral = 0;
        // integral
        for (int k = 0; k < quad_npoints; k++)
        {
            // integral += 0.5*quad_weights[k] * t_quad_vals[k] *(b-a); // *
            // t_xyz_vals[k];
            integral += 0.5 * quad_weights[k] * t_quad_vals[k] *
                        std::abs((b - a)) * t_xyz_vals[k];
        }
        sum += integral;
    }
    valX = sum;
    return true;
}

/*
 * This function adds extra knots at the ends of the filter to have supply the
 * knots at the ends of the filter.
 */
bool SmoothieSIAC2D::v_EvaluateAt_NSK_GivenFilterLength_v2(
    const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
    NekDouble &valX, NekDouble &valY, NekDouble &valZ,
    Array<OneD, NekDouble> &direction, NekDouble meshSpacing, int varNum)
{
    boost::ignore_unused(valY); // reserved for derivative.
    // call HNM
    // TO do . Direction needs to picked up from meshptr. ???
    if (meshSpacing < 0)
    { // No parameter was specified.
        meshSpacing = m_meshSpacing;
    }
    NekDouble meshTShift = 0.0;
    NekDouble tmin, tmax;
    vector<NekDouble> HvalX, HvalY, HvalZ, HvalT, SvalT, TvalT;
    // HNM to get list of break points.
    // check if a symmetric fitler can be applied.
    // in position, symmetric range

    // Get filter range given scaling.
    m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, tmin, tmax);

    bool b_symMesh = m_meshHandlePtr->CanTRangebeApplied(
        PtsX, PtsY, PtsZ, direction, tmin, tmax, meshTShift);
    // Made an assumtion that a symmetric filter can be applied ???
    // return List of mesh breakpoints in that range.
    // call SSS/SSO
    if (!b_symMesh)
    {
        valX = -1;
        return false;
    }

    // This is function which returns HvalT. HvalT is the filter knots.
    // m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ,direction, tmin,tmax,
    // HvalX, HvalY, HvalZ, HvalT);
    m_meshHandlePtr->GetBreakPts_Without_Tmin_Tmax(
        PtsX, PtsY, PtsZ, direction, tmin, tmax, HvalX, HvalY, HvalZ, HvalT);
    valZ = HvalT.size();

    // del vector<int> t_GIDs,t_EIDs;
    // del m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, direction, TvalT,
    // t_GIDs, t_EIDs);
    vector<int> t_GIDs, t_EIDs;
    vector<int> t_h_GIDs, t_h_EIDs;

    // Assuming this might not be needed.
    // m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, HvalT,
    // t_h_GIDs, t_h_EIDs); mergeBreakPtsAdv(HvalT,SvalT, TvalT, t_h_GIDs,
    // t_h_EIDs, t_GIDs, t_EIDs); mergeBreakPts( HvalT, SvalT, TvalT);

    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, HvalT, t_GIDs,
                                   t_EIDs);
    TvalT = HvalT;
    // Loop forsumation.

    NekDouble sum = 0.0;
    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.

    int quadratureOfMesh =
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(0) +
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(1);

    int quad_npoints = ceil((m_order + quadratureOfMesh + 1.0) / 2.0);
    LibUtilities::PointsKey quadPointsKey(
        quad_npoints, LibUtilities::PointsType::eGaussGaussLegendre);

    Array<OneD, NekDouble> quad_points =
        LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    Array<OneD, NekDouble> quad_weights =
        LibUtilities::PointsManager()[quadPointsKey]->GetW();
    Array<OneD, NekDouble> t_quad(quad_npoints), t_x(quad_npoints),
        t_y(quad_npoints), t_z(quad_npoints), t_quad_vals(quad_npoints),
        t_xyz_vals(quad_npoints);
    //	vector<int> t_GIDs,t_EIDs;
    //	m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, direction, TvalT, t_GIDs,
    // t_EIDs);
    Array<OneD, NekDouble> TvalT_NekArray(TvalT.size());
    for (int i = 0; i < TvalT.size(); i++)
    {
        TvalT_NekArray[i] = TvalT[i];
    }
    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        // Create a seg element
        // Actual quadrature points needed are for order = polyorder+Bspline
        // order.
        NekDouble a = TvalT[i];
        NekDouble b = TvalT[i + 1];
        int gID     = t_GIDs[i];
        int eID     = t_EIDs[i];
        for (int j = 0; j < quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (quad_points[j] + 1.0) / 2.0 + a;
            t_x[j]    = PtsX + t_quad[j] * direction[0];
            t_y[j]    = PtsY + t_quad[j] * direction[1];
            t_z[j]    = PtsZ + t_quad[j] * direction[2];
        }
        // Evaluate filter  at these points.
        if (b_symMesh)
        {
            //	m_siacFilterPtrs[m_SymID]->EvaluateFilter(t_quad,t_quad_vals,
            // meshSpacing );
            //		m_siacFilterPtrs[m_SymID]->EvaluateFilter_GivenNumSplines(t_quad,t_quad_vals,
            // meshSpacing );
            m_siacFilterPtrs[m_SymID]->EvaluateFilterWknots_GivenNumSplines(
                t_quad, t_quad_vals, TvalT_NekArray, 1.0);
        }
        else
        {
            //	m_siacFilterPtrs[m_OneID]->EvaluateFilter(t_quad,t_quad_vals,
            // meshSpacing, meshTShift, true);
        }
        m_meshHandlePtr->EvaluateAt(t_x, t_y, t_z, gID, eID, t_xyz_vals,
                                    varNum);
        NekDouble integral = 0;
        // integral
        for (int k = 0; k < quad_npoints; k++)
        {
            // integral += 0.5*quad_weights[k] * t_quad_vals[k] *(b-a); // *
            // t_xyz_vals[k];
            integral += 0.5 * quad_weights[k] * t_quad_vals[k] *
                        std::abs((b - a)) * t_xyz_vals[k];
        }
        sum += integral;
    }
    valX = sum;
    return true;
}

bool SmoothieSIAC2D::v_EvaluateAt_NSK_GivenFilterLength_v1(
    const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
    NekDouble &valX, NekDouble &valY, NekDouble &valZ,
    Array<OneD, NekDouble> &direction, NekDouble meshSpacing, int varNum)
{
    // call HNM
    // TO do . Direction needs to picked up from meshptr. ???
    if (meshSpacing < 0)
    { // No parameter was specified.
        meshSpacing = m_meshSpacing;
    }
    NekDouble meshTShift = 0.0;
    NekDouble tmin, tmax;
    vector<NekDouble> HvalX, HvalY, HvalZ, HvalT, SvalT, TvalT;
    // HNM to get list of break points.
    // check if a symmetric fitler can be applied.
    // in position, symmetric range

    // Get filter range given scaling.
    m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, tmin, tmax);

    bool b_symMesh = m_meshHandlePtr->CanTRangebeApplied(
        PtsX, PtsY, PtsZ, direction, tmin, tmax, meshTShift);
    // Made an assumtion that a symmetric filter can be applied ???
    // return List of mesh breakpoints in that range.
    // call SSS/SSO
    if (!b_symMesh)
    {
        valX = -1;
        return false;
    }

    // This is function which returns HvalT. HvalT is the filter knots.
    m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ, direction, tmin, tmax, HvalX,
                                 HvalY, HvalZ, HvalT);
    valY = HvalT.back() - HvalT.front();
    valZ = HvalT.size();
    // m_meshHandlePtr->GetBreakPts_Without_Tmin_Tmax(PtsX, PtsY,
    // PtsZ,direction, tmin,tmax, HvalX, HvalY, HvalZ, HvalT);

    // del vector<int> t_GIDs,t_EIDs;
    // del m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, direction, TvalT,
    // t_GIDs, t_EIDs);
    vector<int> t_GIDs, t_EIDs;
    vector<int> t_h_GIDs, t_h_EIDs;

    // Assuming this might not be needed.
    // m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, HvalT,
    // t_h_GIDs, t_h_EIDs); mergeBreakPtsAdv(HvalT,SvalT, TvalT, t_h_GIDs,
    // t_h_EIDs, t_GIDs, t_EIDs); mergeBreakPts( HvalT, SvalT, TvalT);

    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, HvalT, t_GIDs,
                                   t_EIDs);
    TvalT = HvalT;
    // Loop forsumation.

    NekDouble sum = 0.0;
    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.

    // cout << "////////Working till here////////" <<endl;
    int quadratureOfMesh =
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(0) +
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(1);
    //	m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(2);

    int quad_npoints = ceil((m_order + quadratureOfMesh + 1.0) / 2.0);
    // int quad_npoints = ceil((m_order+3)/2) +
    // m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetTotPoints()+2; cout <<
    // "quad_npoints" << quad_npoints << endl;
    // LibUtilities::PointsKey quadPointsKey(quad_npoints,
    //						m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetPointsType(0));
    LibUtilities::PointsKey quadPointsKey(
        quad_npoints, LibUtilities::PointsType::eGaussGaussLegendre);

    Array<OneD, NekDouble> quad_points =
        LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    Array<OneD, NekDouble> quad_weights =
        LibUtilities::PointsManager()[quadPointsKey]->GetW();
    Array<OneD, NekDouble> t_quad(quad_npoints), t_x(quad_npoints),
        t_y(quad_npoints), t_z(quad_npoints), t_quad_vals(quad_npoints),
        t_xyz_vals(quad_npoints);
    //	vector<int> t_GIDs,t_EIDs;
    //	m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, direction, TvalT, t_GIDs,
    // t_EIDs);
    Array<OneD, NekDouble> TvalT_NekArray(TvalT.size());
    for (int i = 0; i < TvalT.size(); i++)
    {
        TvalT_NekArray[i] = TvalT[i];
    }
    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        // Create a seg element
        // Actual quadrature points needed are for order = polyorder+Bspline
        // order.
        NekDouble a = TvalT[i];
        NekDouble b = TvalT[i + 1];
        int gID     = t_GIDs[i];
        int eID     = t_EIDs[i];
        for (int j = 0; j < quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (quad_points[j] + 1.0) / 2.0 + a;
            t_x[j]    = PtsX + t_quad[j] * direction[0];
            t_y[j]    = PtsY + t_quad[j] * direction[1];
            t_z[j]    = PtsZ + t_quad[j] * direction[2];
        }
        // Evaluate filter  at these points.
        if (b_symMesh)
        {
            //	m_siacFilterPtrs[m_SymID]->EvaluateFilter(t_quad,t_quad_vals,
            // meshSpacing );
            //		m_siacFilterPtrs[m_SymID]->EvaluateFilter_GivenNumSplines(t_quad,t_quad_vals,
            // meshSpacing );
            m_siacFilterPtrs[m_SymID]->EvaluateFilterWknots_GivenNumSplines(
                t_quad, t_quad_vals, TvalT_NekArray, 1.0);
        }
        else
        {
            //	m_siacFilterPtrs[m_OneID]->EvaluateFilter(t_quad,t_quad_vals,
            // meshSpacing, meshTShift, true);
        }
        m_meshHandlePtr->EvaluateAt(t_x, t_y, t_z, gID, eID, t_xyz_vals,
                                    varNum);
        NekDouble integral = 0;
        // integral
        for (int k = 0; k < quad_npoints; k++)
        {
            // integral += 0.5*quad_weights[k] * t_quad_vals[k] *(b-a); // *
            // t_xyz_vals[k];
            integral += 0.5 * quad_weights[k] * t_quad_vals[k] *
                        std::abs((b - a)) * t_xyz_vals[k];
        }
        sum += integral;
    }
    valX = sum;
    return true;
}

bool SmoothieSIAC2D::v_EvaluateRecursiveAt(
    const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
    NekDouble &valX, NekDouble &valY, NekDouble &valZ,
    const vector<std::shared_ptr<SmoothieSIAC>> &Sms,
    vector<Array<OneD, NekDouble>> directions,
    const vector<NekDouble> &meshSpacings, const vector<int> &varNums,
    const int curLevel)
{
    int totLevels = directions.size();
    assert((totLevels > curLevel) && (curLevel >= 0) &&
           "Some parameters are not right.");
    NekDouble meshSpacing            = meshSpacings[curLevel];
    Array<OneD, NekDouble> direction = directions[curLevel];
    if (totLevels - 1 == curLevel)
    {
        return Sms[curLevel]->EvaluateAt(
            PtsX, PtsY, PtsZ, valX, valY, valZ, directions[curLevel],
            meshSpacings[curLevel], varNums[curLevel]);
    }

    // call HNM
    // TO do . Direction needs to picked up from meshptr. ???
    if (meshSpacing < 0)
    { // No parameter was specified.
        meshSpacing = m_meshSpacing;
    }
    NekDouble meshTShift = 0.0;
    NekDouble tmin, tmax;
    vector<NekDouble> HvalX, HvalY, HvalZ, HvalT, SvalT, TvalT;
    // HNM to get list of break points.
    // check if a symmetric fitler can be applied.
    // in position, symmetric range

    // Get filter range given scaling.
    m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, tmin, tmax);

    bool b_symMesh = m_meshHandlePtr->CanTRangebeApplied(
        PtsX, PtsY, PtsZ, direction, tmin, tmax, meshTShift);
    // Made an assumtion that a symmetric filter can be applied ???
    // return List of mesh breakpoints in that range.
    // call SSS/SSO
    if (b_symMesh)
    {
        //	cout << "into symmetric condition loop" << endl;
        m_siacFilterPtrs[m_SymID]->GetBreakPts(
            meshSpacing, SvalT); //"??? specify parameters"
    }
    else
    {
        // cout << "PtsX: "<< PtsX << " before tmin : " << tmin;
        // cout << " tmax : " << tmax << endl;
        // cout << "after tmin : " << tmin;
        // cout << " tmax : " << tmax << endl;
        SvalT.clear();
        if (m_OneID >= 0)
        { // OneSided Filter is defined and given by .
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, tmin, tmax);
            m_meshHandlePtr->CanTRangebeApplied(PtsX, PtsY, PtsZ, direction,
                                                tmin, tmax, meshTShift);
            // Next step could be equal to tmin = tmin+meshTShift; tmax =
            // tmax+meshTShift;
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, tmin, tmax,
                                                      meshTShift);
            m_siacFilterPtrs[m_OneID]->GetBreakPts(
                meshSpacing, SvalT, meshTShift); //"??? specify parameters"
                                                 // printNekArray( SvalT, 0);
            //	cout << "ptsX: " << PtsX ;
            //	cout << "\tmeshTShift: " << meshTShift << endl;
        }
        else
        { // OneSided Filter is not defined.
            // cout << "Symmetric filter does not fit and OneSided is not
            // defined here." << endl; cout << "Ideally keep the original value.
            // -1 is returned at end " << endl;
            valX = -1;
            return false;
        }
    }

    m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ, direction, tmin, tmax, HvalX,
                                 HvalY, HvalZ, HvalT);
    // del vector<int> t_GIDs,t_EIDs;
    // del m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, direction, TvalT,
    // t_GIDs, t_EIDs);
    vector<int> t_GIDs, t_EIDs;
    vector<int> t_h_GIDs, t_h_EIDs;
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, HvalT, t_h_GIDs,
                                   t_h_EIDs);
    mergeBreakPtsAdv(HvalT, SvalT, TvalT, t_h_GIDs, t_h_EIDs, t_GIDs, t_EIDs);
    // mergeBreakPts( HvalT, SvalT, TvalT);

    // BPTS
    // printNekArray(HvalT,0);
    // printNekArray(SvalT,0);
    // printNekArray(TvalT,0);

    // Loop forsumation.

    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.

    // cout << "////////Working till here////////" <<endl;
    int quadratureOfMesh =
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(0) +
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(1);
    //	m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(2);

    int quad_npoints = ceil((m_order + quadratureOfMesh + 1.0) / 2.0);
    // int quad_npoints = ceil((m_order+3)/2) +
    // m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetTotPoints()+2; cout <<
    // "quad_npoints" << quad_npoints << endl; LibUtilities::PointsKey
    // quadPointsKey(quad_npoints,
    //					m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetPointsType(0));
    LibUtilities::PointsKey quadPointsKey(
        quad_npoints, LibUtilities::PointsType::eGaussGaussLegendre);

    Array<OneD, NekDouble> quad_points =
        LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    Array<OneD, NekDouble> quad_weights =
        LibUtilities::PointsManager()[quadPointsKey]->GetW();
    Array<OneD, NekDouble> t_quad(quad_npoints), t_x(quad_npoints),
        t_y(quad_npoints), t_z(quad_npoints), t_quad_vals(quad_npoints),
        t_xyz_vals(quad_npoints);
    //	vector<int> t_GIDs,t_EIDs;
    //	m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, direction, TvalT, t_GIDs,
    // t_EIDs);

    Array<OneD, NekDouble> integralAry(TvalT.size(), 0.0);
    NekDouble sum = 0.0;

    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        // Create a seg element
        // Actual quadrature points needed are for order = polyorder+Bspline
        // order.
        NekDouble a = TvalT[i];
        NekDouble b = TvalT[i + 1];
        for (int j = 0; j < quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (quad_points[j] + 1.0) / 2.0 + a;
            t_x[j]    = PtsX + t_quad[j] * direction[0];
            t_y[j]    = PtsY + t_quad[j] * direction[1];
            t_z[j]    = PtsZ + t_quad[j] * direction[2];
        }
        // Evaluate filter  at these points.
        if (b_symMesh)
        {
            m_siacFilterPtrs[m_SymID]->EvaluateFilter(t_quad, t_quad_vals,
                                                      meshSpacing);
        }
        else
        {
            m_siacFilterPtrs[m_OneID]->EvaluateFilter(
                t_quad, t_quad_vals, meshSpacing, meshTShift, true);
        }
        // m_meshHandlePtr->EvaluateAt( t_x, t_y,t_z,gID,eID,
        // t_xyz_vals,varNum);
        for (int j = 0; j < quad_npoints; j++)
        {
            Sms[curLevel]->EvaluateRecursiveAt(
                t_x[j], t_y[j], t_z[j], t_xyz_vals[j], valY, valZ, Sms,
                directions, meshSpacings, varNums, curLevel + 1);
        }

        // NekDouble integral = 0;
        // integral
        for (int k = 0; k < quad_npoints; k++)
        {
            // integral += 0.5*quad_weights[k] * t_quad_vals[k] *(b-a); // *
            // t_xyz_vals[k]; integral += 0.5*quad_weights[k] * t_quad_vals[k]
            // *std::abs((b-a)) *  t_xyz_vals[k];
            integralAry[i] += 0.5 * quad_weights[k] * t_quad_vals[k] *
                              std::abs((b - a)) * t_xyz_vals[k];
        }
    }
    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        sum += integralAry[i];
    }
    valX = sum;
    return true;
}

bool SmoothieSIAC2D::v_EvaluateAt_SymY(const NekDouble PtsX,
                                       const NekDouble PtsY,
                                       const NekDouble PtsZ, NekDouble &valX,
                                       NekDouble &valY, NekDouble &valZ,
                                       Array<OneD, NekDouble> &direction,
                                       NekDouble meshSpacing, int varNum)
{
    // call HNM
    // TO do . Direction needs to picked up from meshptr. ???
    if (meshSpacing < 0)
    { // No parameter was specified.
        meshSpacing = m_meshSpacing;
    }
    NekDouble meshTShift = 0.0;
    NekDouble tmin, tmax;
    vector<NekDouble> HvalX, HvalY, HvalZ, HvalT, SvalT, TvalT;
    // HNM to get list of break points.
    // check if a symmetric fitler can be applied.
    // in position, symmetric range

    // Get filter range given scaling.
    m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, tmin, tmax);

    bool b_symMesh = m_meshHandlePtr->CanTRangebeApplied(
        PtsX, PtsY, PtsZ, direction, tmin, tmax, meshTShift);
    // Made an assumtion that a symmetric filter can be applied ???
    // return List of mesh breakpoints in that range.
    // call SSS/SSO
    if (b_symMesh)
    {
        //	cout << "into symmetric condition loop" << endl;
        // m_siacFilterPtrs[m_SymID]->GetBreakPts(meshSpacing , SvalT); //"???
        // specify parameters"
        return v_EvaluateAt(PtsX, PtsY, PtsZ, valX, valY, valZ, direction,
                            meshSpacing, varNum);
    }
    else
    {
        valX = -0.1;
        return false;
    }
}

bool SmoothieSIAC2D::v_EvaluateL2UsingLineAt(
    const vector<NekDouble> &stPoint, const Array<OneD, NekDouble> &direction,
    const int n_quadPts, const NekDouble meshSpacing,
    const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int varNum)
{
    NekDouble meshTShift = 0.0;
    NekDouble tmin, tmax;
    // vector<NekDouble> SvalT,TvalT;
    NekDouble PtsX = stPoint[0];
    NekDouble PtsY = stPoint[1];
    NekDouble PtsZ = stPoint[2];

    m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, tmin, tmax);

    bool b_symMesh = m_meshHandlePtr->CanTRangebeApplied(
        PtsX, PtsY, PtsZ, direction, tmin, tmax, meshTShift);

    if (b_symMesh)
    {
        // m_siacFilterPtrs[m_SymID]->GetBreakPts(meshSpacing , SvalT); //"???
        // specify parameters"
    }
    else
    {
        // SvalT.clear();
        if (m_OneID >= 0)
        {
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, tmin, tmax);
            m_meshHandlePtr->CanTRangebeApplied(PtsX, PtsY, PtsZ, direction,
                                                tmin, tmax, meshTShift);
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, tmin, tmax,
                                                      meshTShift);
            //      m_siacFilterPtrs[m_OneID]->GetBreakPts(meshSpacing ,
            //      SvalT,meshTShift); //"??? specify parameters"
        }
        else
        {
            assert(false && "Should not be using one sided filter. ");
            return false;
        }
    }

    // Now tmin and tmax are boundaries needed for mesh Handler.

    tvals.clear();
    // set up line points.
    // evaluate funtion values at these points.
    vector<NekDouble> HvalT;
    vector<int> Gids, Eids;
    Array<OneD, NekDouble> t_LineElm;
    NekDouble t_mesh_min = tmin;
    NekDouble t_mesh_max = tmax;
    // Set up the quadrature.
    SetupLineForLSIAC(direction, stPoint, t_mesh_min, t_mesh_max, n_quadPts,
                      HvalT, Gids, Eids, t_LineElm);
    Array<OneD, NekDouble> tv_LineElm(t_LineElm.num_elements());
    // Get the values of quadrature.
    GetVLineForLSIAC(n_quadPts, stPoint, direction, HvalT, Gids, Eids,
                     t_LineElm, tv_LineElm, varNum);

    // Set up one more array for quadrature weights or jacobian values.
    // tW_LineElm
    // tJ_LineElm
    // NekDouble meshSpacing = 0.1; // NekDouble.
    NekDouble t_min, t_max;
    // NekDouble meshTShift;
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
        for (int j = 0; j < HvalT.size(); j++)
        {
            if (HvalT[j] > t + t_min + TOLERENCE)
            {
                startElmIndex = j - 1;
                break;
            }
        }

        for (int j = 0; j < HvalT.size(); j++)
        {
            if (HvalT[j] > t + t_max - TOLERENCE)
            {
                endElmIndex = j - 1;
                break;
            }
        }

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
    // cout << "part1\t"<< el1 << "\t\t" << el2 << "\t\t" << el3 << "\tTotal\t"
    // << elT << endl;
    return true;
}

bool SmoothieSIAC2D::v_EvaluateL2UsingLineAt_v3(
    const vector<NekDouble> &stPoint, const Array<OneD, NekDouble> &direction,
    const int n_quadPts, const int n_quadPts_resample,
    const NekDouble meshSpacing, const vector<NekDouble> &tparams,
    vector<NekDouble> &tvals, int varNum)
{

    tvals.clear();
    NekDouble meshTShift = 0.0;
    NekDouble tmin, tmax;
    // vector<NekDouble> SvalT,TvalT;
    NekDouble PtsX = stPoint[0];
    NekDouble PtsY = stPoint[1];
    NekDouble PtsZ = stPoint[2];

    m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, tmin, tmax);

    bool b_symMesh = m_meshHandlePtr->CanTRangebeApplied(
        PtsX, PtsY, PtsZ, direction, tmin, tmax, meshTShift);

    if (b_symMesh)
    {
        // m_siacFilterPtrs[m_SymID]->GetBreakPts(meshSpacing , SvalT); //"???
        // specify parameters"
    }
    else
    {
        // SvalT.clear();
        if (m_OneID >= 0)
        {
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, tmin, tmax);
            m_meshHandlePtr->CanTRangebeApplied(PtsX, PtsY, PtsZ, direction,
                                                tmin, tmax, meshTShift);
            m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, tmin, tmax,
                                                      meshTShift);
            //      m_siacFilterPtrs[m_OneID]->GetBreakPts(meshSpacing ,
            //      SvalT,meshTShift); //"??? specify parameters"
        }
        else
        {
            assert(false && "Should not be using one sided filter. ");
            return false;
        }
    }

    // set up line points.
    // evaluate funtion values at these points.
    NekDouble t_mesh_min = tmin;
    NekDouble t_mesh_max =
        tmax; //*(std::max_element( tparams.begin(), tparams.end()));
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
        int startElmIndex = -1, endElmIndex = -1;
        // start ElementId and End Element Id;
        for (int j = 0; j < HvalT.size(); j++)
        {
            if (HvalT[j] > t + t_min + TOLERENCE)
            {
                startElmIndex = j - 1;
                break;
            }
        }

        for (int j = 0; j < HvalT.size(); j++)
        {
            if (HvalT[j] > t + t_max - TOLERENCE)
            {
                endElmIndex = j - 1;
                break;
            }
        }
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

        /*


        for (int el = startElmIndex; el<=endElmIndex; el++)
        {

                // Get break points for filter.


                for (int qp = 0; qp < n_quadPts; qp++)
                {
                        t_quadPts[qp] = t_LineElm[el*n_quadPts+qp] -t;
                }
                if(b_symMesh)
                {
                        m_siacFilterPtrs[m_SymID]->EvaluateFilter(t_quadPts,t_quad_vals,
        meshSpacing ); }else
                {
                        m_siacFilterPtrs[m_OneID]->EvaluateFilter(t_quadPts,t_quad_vals,
        meshSpacing, meshTShift, true);
                }
                NekDouble integral = 0.0;
                for (int qp = 0; qp < n_quadPts; qp++)
                {
                        integral += 0.5*t_quad_weights[qp]*
        tv_LineElm[el*n_quadPts+qp]* std::abs(HvalT[el+1]-HvalT[el])*
        t_quad_vals[qp];
                }
                sum += integral;
        }
        tvals.push_back(sum);
        */
    }

    return true;
}

bool SmoothieSIAC2D::v_Cal_NUK_ConstMetricTensor(
    const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
    const NekDouble meshSpacing, Array<OneD, NekDouble> &direction,
    Array<OneD, NekDouble> &knotVec)
{
    vector<NekDouble> HvalX, HvalY, HvalZ, HvalT, SvalT; //,TvalT,LvalT;
    // Array<OneD,NekDouble> coord(3,0.0);
    // coord[0] = PtsX; coord[1] = PtsY; coord[2] = PtsZ;
    NekDouble tmin, tmax;
    m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, tmin, tmax);
    m_siacFilterPtrs[m_SymID]->GetBreakPts(1.0,
                                           SvalT); //"??? specify parameters"

    // UK_pos
    vector<NekDouble>::iterator UK_pos_iter;
    for (std::vector<NekDouble>::iterator it = SvalT.begin(); it != SvalT.end();
         ++it)
    {
        if (*it > 0.0 + 1e-9)
        {
            UK_pos_iter = it;
            break;
        }
    }

    // Find NUK on the positive side
    m_meshHandlePtr->GetBreakPts_Without_Tmin_Tmax(
        PtsX, PtsY, PtsZ, direction, 0.0, tmax, HvalX, HvalY, HvalZ, HvalT);
    if (HvalT.size() == 0)
    {
        return false;
    }
    assert(HvalT.size() > 0 && "Make sure more neightbourhood is searched");
    vector<int> t_GIDs, t_EIDs;

    HvalT.insert(HvalT.begin(),
                 0.0); // Adding zero to make sure GetListOfGIDs work.
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, HvalT, t_GIDs,
                                   t_EIDs);
    t_GIDs.pop_back();
    t_EIDs.pop_back();
    HvalT.erase(
        HvalT.begin()); // Removing the first element to maintain consistency.

    // Design a method which return direction for a list of element ratios for
    // the line segments. Use Metric Tensor here.
    vector<NekDouble> scale; // initializing them to 1.0 for test case. (Can do
                             // mesh scaling also.)
    m_meshHandlePtr->GetMTScalingOfGIDs(t_EIDs, direction, scale);

    // Need ratios also.
    vector<NekDouble> tvalues = HvalT;
    vector<NekDouble> NUKpos, NUKneg;

    // Determine zero knot
    if ((m_order % 2) == 0)
    {
        NUKpos.push_back(0.0);
    }

    // Initialize loop
    NekDouble rem_KL, rem_EL, accum_KL, step;
    int elmIndex = 0; // location of the element needs to calculated.
    rem_KL       = *UK_pos_iter;
    accum_KL     = 0.0;
    rem_EL       = tvalues[elmIndex];
    bool b_inMesh;
    int multi = 0;

    while (true)
    {
        step = rem_EL; // Here not req to mul with scale[elmIndex];
        if (step > rem_KL * scale[elmIndex])
        {
            // cout << "Debug Info: Over elementBoundary" << endl;
            // cout << "rem_KL\trem_EL\taccum_KL\telmIndex" <<endl;
            // cout << rem_KL<< "\t" << rem_EL<<"\t"<<accum_KL<<"\t"<<elmIndex
            // <<endl;
            accum_KL += rem_KL * scale[elmIndex];
            // cout << "rem_KL\t\taccum_KL" <<endl;
            // cout << rem_KL - rem_KL<< "\t\t"<<accum_KL<<"\t"<<accum_KL
            // <<endl;
            rem_EL -= rem_KL * scale[elmIndex];
            NUKpos.push_back(accum_KL);
            UK_pos_iter++;
            if (UK_pos_iter != SvalT.end())
            {
                accum_KL = 0.0;
                rem_KL   = *UK_pos_iter - *(UK_pos_iter - 1);
            }
            else
            {
                break;
            }
        }
        else // step <= rem_KL*scale[elmIndex]
        {
            ////cout << "Debug Info: Same Element" << endl;
            ////cout << "rem_KL\trem_EL\taccum_KL\telmIndex" <<endl;
            ////cout << rem_KL<< "\t" << rem_EL<<"\t"<<accum_KL<<"\t"<<elmIndex
            ///<<endl;
            accum_KL += rem_EL;
            rem_KL = rem_KL - (rem_EL) / scale[elmIndex];
            // cout << "rem_KL\t\taccum_KL Next Elm" <<endl;
            // cout << rem_KL<< "\t\t"<<accum_KL<<"\t"<<
            // rem_KL*scale[elmIndex]+accum_KL <<"\t\t"<<rem_EL<<endl;
            // Check if the next element there or not.
            if (elmIndex + 1 >= tvalues.size())
            {
                // assert("Not coded yet");
                b_inMesh = true;
                HvalT.clear();
                while (HvalT.size() == 0 && b_inMesh == true)
                {
                    multi = multi + 1;
                    m_meshHandlePtr->GetBreakPts_Without_Tmin_Tmax(
                        PtsX, PtsY, PtsZ, direction, tmax * multi,
                        tmax * (multi + 1), HvalX, HvalY, HvalZ, HvalT);
                    b_inMesh = m_meshHandlePtr->CanTRangebeAppliedWOMeshShift(
                        PtsX, PtsY, PtsZ, direction, tmax * multi,
                        tmax * (multi + 1));
                }
                if (HvalT.size() == 0 && b_inMesh == false)
                {
                    // Debug info
                    // cout << "Exiting because at end of mesh.";
                    // cout << NUKpos.size() << endl;
                    // cout << multi << endl;
                    // for (auto it = NUKpos.begin(); it!=NUKpos.end(); ++it)
                    //{
                    //	cout << *it << endl;
                    //}
                    // cout << "OneSideFilter Not done yet " << endl;
                    return false;
                }
                else
                {
                    HvalT.insert(HvalT.begin(),
                                 tmax * multi); // Adding zero to make sure
                                                // GetListOfGIDs work.
                    t_GIDs.clear();
                    t_EIDs.clear();
                    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction,
                                                   HvalT, t_GIDs, t_EIDs);
                    t_GIDs.pop_back();
                    t_EIDs.pop_back();
                    HvalT.erase(HvalT.begin()); // Removing the first element to
                                                // maintain consistency.
                    vector<NekDouble>
                        scale_append; // initializing them to 1.0 for test case.
                                      // (Can do mesh scaling also.)
                    m_meshHandlePtr->GetMTScalingOfGIDs(t_EIDs, direction,
                                                        scale_append);

                    tvalues.insert(tvalues.end(), HvalT.begin(), HvalT.end());
                    scale.insert(scale.end(), scale_append.begin(),
                                 scale_append.end());
                }
            }
            elmIndex++;
            rem_EL = tvalues[elmIndex] - tvalues[elmIndex - 1];
        }
    }

    // Find NUK on the negative side
    m_meshHandlePtr->GetBreakPts_Without_Tmin_Tmax(
        PtsX, PtsY, PtsZ, direction, tmin, 0.0, HvalX, HvalY, HvalZ, HvalT);
    if (HvalT.size() == 0)
    {
        return false;
    }
    assert(HvalT.size() > 0 && "Make sure more neightbourhood is searched");
    t_GIDs.clear();
    t_EIDs.clear();

    HvalT.push_back(0.0); // Adding zero to make sure GetListOfGIDs work.
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, HvalT, t_GIDs,
                                   t_EIDs);
    t_GIDs.pop_back();
    t_EIDs.pop_back();
    HvalT.pop_back(); // Removing the last element to maintain consistency.

    // Design a method which return direction for a list of element ratios for
    // the line segments. Use Metric Tensor here.
    scale.clear(); // initializing them to 1.0 for test case. (Can do mesh
                   // scaling also.)
    m_meshHandlePtr->GetMTScalingOfGIDs(t_EIDs, direction, scale);

    // Need ratios also.
    tvalues = HvalT;
    std::reverse(std::begin(tvalues), std::end(tvalues));

    // Initialize loop
    rem_KL   = 0;
    rem_EL   = 0;
    accum_KL = 0;
    step     = 0;
    elmIndex = 0; // location of the element needs to calculated.
                  // UK_neg_iter
    vector<NekDouble>::reverse_iterator UK_neg_iter;
    for (std::vector<NekDouble>::reverse_iterator it = SvalT.rbegin();
         it != SvalT.rend(); ++it)
    {
        if (*it < 0.0 - 1e-9)
        {
            UK_neg_iter = it;
            break;
        }
    }
    rem_KL   = *(UK_neg_iter - 1) - *UK_neg_iter;
    accum_KL = 0.0;
    rem_EL   = -1 * tvalues[elmIndex]; // tvalues are negative.

    multi = 0;
    while (true)
    {
        step = rem_EL; // Here not req to mul with scale[elmIndex];
        if (step > rem_KL * scale[elmIndex])
        {
            // cout << "Debug Info: Over elementBoundary" << endl;
            // cout << "rem_KL\trem_EL\taccum_KL\telmIndex" <<endl;
            // cout << rem_KL<< "\t" << rem_EL<<"\t"<<accum_KL<<"\t"<<elmIndex
            // <<endl;
            accum_KL += rem_KL * scale[elmIndex];
            // cout << "rem_KL\t\taccum_KL" <<endl;
            // cout << rem_KL - rem_KL<< "\t\t"<<accum_KL<<"\t"<<accum_KL
            // <<endl;
            rem_EL -= rem_KL * scale[elmIndex];
            NUKneg.push_back(-1.0 * accum_KL);
            // NUKneg.insert(NUKneg.begin(),-1*accum_KL-1e-3); // accum_KL is
            // positive and knots should be negtive here.
            UK_neg_iter++;
            if (UK_neg_iter != SvalT.rend())
            {
                accum_KL = 0.0;
                rem_KL   = *(UK_neg_iter - 1) - *(UK_neg_iter);
            }
            else
            {
                break;
            }
        }
        else // step <= rem_KL*scale[elmIndex]
        {
            // cout << "Debug Info: Same Element" << endl;
            // cout << "rem_KL\trem_EL\taccum_KL\telmIndex" <<endl;
            // cout << rem_KL<< "\t" << rem_EL<<"\t"<<accum_KL<<"\t"<<elmIndex
            // <<endl;
            accum_KL += rem_EL;
            rem_KL = rem_KL - (rem_EL) / scale[elmIndex];
            // cout << "rem_KL\t\taccum_KL Next Elm" <<endl;
            // cout << rem_KL<< "\t\t"<<accum_KL<<"\t"<<
            // rem_KL*scale[elmIndex]+accum_KL <<"\t\t"<<rem_EL<<endl;
            // Check if the next element there or not.
            if (elmIndex + 1 >= tvalues.size())
            {
                // assert("Not coded yet");
                b_inMesh = true;
                HvalT.clear();
                while (HvalT.size() == 0 && b_inMesh == true)
                {
                    multi = multi + 1;
                    m_meshHandlePtr->GetBreakPts_Without_Tmin_Tmax(
                        PtsX, PtsY, PtsZ, direction, tmin * (multi + 1),
                        tmin * (multi), HvalX, HvalY, HvalZ, HvalT);
                    b_inMesh = m_meshHandlePtr->CanTRangebeAppliedWOMeshShift(
                        PtsX, PtsY, PtsZ, direction, tmin * (multi + 1),
                        tmin * (multi));
                }
                if (HvalT.size() == 0 && b_inMesh == false)
                {
                    // Debug info
                    // cout << "Exiting because at end of mesh.";
                    // cout << NUKneg.size() << endl;
                    // cout << multi << endl;
                    for (auto it = NUKneg.begin(); it != NUKneg.end(); ++it)
                    {
                        // cout << *it << endl;
                    }
                    return false;
                }
                else
                {
                    // HvalT.insert(HvalT.rbegin(),tmin*multi); // Adding zero
                    // to make sure GetListOfGIDs work.
                    HvalT.push_back(
                        tmin *
                        multi); // Adding zero to make sure GetListOfGIDs work.
                    t_GIDs.clear();
                    t_EIDs.clear();
                    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction,
                                                   HvalT, t_GIDs, t_EIDs);
                    // HvalT.erase(HvalT.rbegin()); //Removing the first element
                    // to maintain consistency.
                    HvalT.pop_back(); // erase(HvalT.rbegin()); //Removing the
                                      // first element to maintain consistency.
                    vector<NekDouble>
                        scale_append; // initializing them to 1.0 for test case.
                                      // (Can do mesh scaling also.)
                    m_meshHandlePtr->GetMTScalingOfGIDs(t_EIDs, direction,
                                                        scale_append);

                    tvalues.insert(tvalues.end(), HvalT.rbegin(), HvalT.rend());
                    scale.insert(scale.end(), scale_append.rbegin(),
                                 scale_append.rend());
                }
            }
            elmIndex++;
            rem_EL = tvalues[elmIndex - 1] - tvalues[elmIndex];
        }
    }

    vector<NekDouble> NUK(NUKpos.size());
    vector<NekDouble> NUKNEG(NUKneg.size());
    std::partial_sum(NUKpos.begin(), NUKpos.end(), NUK.begin());
    std::partial_sum(NUKneg.begin(), NUKneg.end(), NUKNEG.begin());
    for (auto it = NUKNEG.begin(); it != NUKNEG.end(); ++it)
    {
        NUK.insert(NUK.begin(), *it);
    }
    for (int i = 0; i < NUK.size(); i++)
    {
        knotVec[i] = NUK[i];
        // knotVec[i] = -1*NUK[NUK.size()-i-1];
        // cout << knotVec[i]<<endl;
    }
    // Need to rewrite NUK into knotVec;
    return true;
}

bool SmoothieSIAC2D::v_EvaluateAt_NUK_MetricTensor(
    const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
    NekDouble &valX, NekDouble &valY, NekDouble &valZ,
    Array<OneD, NekDouble> &direction, NekDouble meshSpacing, int varNum)
{
    boost::ignore_unused(valY, valZ); // reserved for derivative.
    // call HNM
    // TO do . Direction needs to picked up from meshptr. ???
    if (meshSpacing < 0)
    { // No parameter was specified.
        meshSpacing = m_meshSpacing;
    }
    NekDouble tmin, tmax;
    // m_meshHandlePtr->GetBreakPts_Without_Tmin_Tmax(PtsX, PtsY,
    // PtsZ,direction, tmin,tmax, HvalX, HvalY, HvalZ, HvalT);

    Array<OneD, NekDouble> knotVec(3 * (m_order - 1) + 2, 0.0);
    bool b_symMesh = this->Cal_NUK_ConstMetricTensor(
        PtsX, PtsY, PtsZ, m_meshSpacing, direction, knotVec);

    // Rewrite the whole code above.
    // Get information for kontVec, b_sysMesh, given Metric tensor for all the
    // elment;
    vector<NekDouble> HvalX, HvalY, HvalZ, SvalT, TvalT, LvalT;

    if (b_symMesh)
    {
        SvalT.insert(SvalT.end(), &knotVec[0], &knotVec[3 * (m_order - 1) + 2]);
        tmin = SvalT.front();
        tmax = SvalT.back();
    }
    else
    {
        assert("Should have been returned before coming here.");
        valX = -1;
        return false;
    }

    m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ, direction, tmin, tmax, HvalX,
                                 HvalY, HvalZ, LvalT);
    // valY = tmax-tmin;  // for debugging
    // valZ = knotVec.num_elements(); // for debugging

    vector<int> t_GIDs, t_EIDs;
    vector<int> t_h_GIDs, t_h_EIDs;
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, LvalT, t_h_GIDs,
                                   t_h_EIDs);
    mergeBreakPtsAdv(LvalT, SvalT, TvalT, t_h_GIDs, t_h_EIDs, t_GIDs, t_EIDs);

    NekDouble sum = 0.0;
    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.
    int quadratureOfMesh =
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(0) +
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(1);
    //	m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(2);

    int quad_npoints = ceil((m_order + quadratureOfMesh + 1.0) / 2.0);
    LibUtilities::PointsKey quadPointsKey(
        quad_npoints, LibUtilities::PointsType::eGaussGaussLegendre);
    Array<OneD, NekDouble> quad_points =
        LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    Array<OneD, NekDouble> quad_weights =
        LibUtilities::PointsManager()[quadPointsKey]->GetW();
    Array<OneD, NekDouble> t_quad(quad_npoints), t_x(quad_npoints),
        t_y(quad_npoints), t_z(quad_npoints), t_quad_vals(quad_npoints),
        t_xyz_vals(quad_npoints);

    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        // Create a seg element
        // Actual quadrature points needed are for order = polyorder+Bspline
        // order.
        NekDouble a = TvalT[i];
        NekDouble b = TvalT[i + 1];
        int gID     = t_GIDs[i];
        int eID     = t_EIDs[i];
        for (int j = 0; j < quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (quad_points[j] + 1.0) / 2.0 + a;
            t_x[j]    = PtsX + t_quad[j] * direction[0];
            t_y[j]    = PtsY + t_quad[j] * direction[1];
            t_z[j]    = PtsZ + t_quad[j] * direction[2];
        }
        // Evaluate filter  at these points.
        if (b_symMesh)
        {
            //	m_siacFilterPtrs[m_SymID]->EvaluateFilter(t_quad,t_quad_vals,
            // meshSpacing );
            //		m_siacFilterPtrs[m_SymID]->EvaluateFilter_GivenNumSplines(t_quad,t_quad_vals,
            // meshSpacing );
            m_siacFilterPtrs[m_SymID]->EvaluateFilterWknots(t_quad, t_quad_vals,
                                                            knotVec, 1.0);
        }
        else
        {
            //	m_siacFilterPtrs[m_OneID]->EvaluateFilter(t_quad,t_quad_vals,
            // meshSpacing, meshTShift, true);
        }
        m_meshHandlePtr->EvaluateAt(t_x, t_y, t_z, gID, eID, t_xyz_vals,
                                    varNum);
        NekDouble integral = 0;
        // integral
        for (int k = 0; k < quad_npoints; k++)
        {
            // integral += 0.5*quad_weights[k] * t_quad_vals[k] *(b-a); // *
            // t_xyz_vals[k];
            integral += 0.5 * quad_weights[k] * t_quad_vals[k] *
                        std::abs((b - a)) * t_xyz_vals[k];
        }
        sum += integral;
    }
    valX = sum;
    return true;
}

bool SmoothieSIAC2D::v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                                  const NekDouble PtsZ, NekDouble &valX,
                                  NekDouble &valY, NekDouble &valZ,
                                  vector<Array<OneD, NekDouble>> &directions,
                                  NekDouble meshSpacing, int varNum)
{
    boost::ignore_unused(valY, valZ); // reserved for derivative.
    // call HNM
    // TO do . Direction needs to picked up from meshptr. ???
    if (meshSpacing < 0)
    { // No parameter was specified.
        meshSpacing = m_meshSpacing;
    }
    NekDouble meshTShift = 0.0;
    NekDouble tmin, tmax;
    vector<NekDouble> HvalX, HvalY, HvalZ, HvalT, SvalT, TvalT;
    // HNM to get list of break points.
    // check if a symmetric fitler can be applied.
    // in position, symmetric range

    // Get filter range given scaling.
    m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, tmin, tmax);

    bool b_symMesh = false;
    Array<OneD, NekDouble> direction(3, 0.0);
    for (int i = 0; i < directions.size(); i++)
    {
        direction = directions[i];
        b_symMesh = m_meshHandlePtr->CanTRangebeApplied(
            PtsX, PtsY, PtsZ, direction, tmin, tmax, meshTShift);
        if (b_symMesh)
            break;
    }
    // bool b_symMesh =
    // m_meshHandlePtr->CanTRangebeApplied(PtsX,PtsY,PtsZ,direction,tmin,tmax,meshTShift);
    // Made an assumtion that a symmetric filter can be applied ???
    // return List of mesh breakpoints in that range.
    // call SSS/SSO
    if (b_symMesh)
    {
        //	cout << "into symmetric condition loop" << endl;
        m_siacFilterPtrs[m_SymID]->GetBreakPts(
            meshSpacing, SvalT); //"??? specify parameters"
    }
    else
    {
        // cout << "PtsX: "<< PtsX << " before tmin : " << tmin;
        // cout << " tmax : " << tmax << endl;
        // cout << "after tmin : " << tmin;
        // cout << " tmax : " << tmax << endl;
        SvalT.clear();
        if (m_OneID >= 0)
        { // OneSided Filter is defined and given by .
            bool foundOneDirection = false;
            for (int i = 0; i < directions.size(); i++)
            {
                direction = directions[i];
                m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, tmin,
                                                          tmax);
                m_meshHandlePtr->CanTRangebeApplied(PtsX, PtsY, PtsZ, direction,
                                                    tmin, tmax, meshTShift);
                m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing, tmin,
                                                          tmax, meshTShift);
                if (m_meshHandlePtr->CanTRangebeAppliedWOMeshShift(
                        PtsX, PtsY, PtsZ, direction, tmin, tmax))
                {
                    foundOneDirection = true;
                    break;
                }
            }
            // Next step could be equal to tmin = tmin+meshTShift; tmax =
            // tmax+meshTShift;
            if (foundOneDirection)
            {
                m_meshHandlePtr->EvaluateAt(PtsX, PtsY, PtsZ, -1, -1, valX,
                                            varNum);
                return false;
            }
            m_siacFilterPtrs[m_OneID]->GetBreakPts(
                meshSpacing, SvalT, meshTShift); //"??? specify parameters"
        }
        else
        { // OneSided Filter is not defined.
            // cout << "Symmetric filter does not fit and OneSided is not
            // defined here." << endl; cout << "Ideally keep the original value.
            // -1 is returned at end " << endl;
            // m_meshHandlePtr->EvaluateAt(PtsX,PtsY,PtsZ,-1,-1,valX,varNum);
            valX = -1;
            return false;
        }
    }

    m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ, direction, tmin, tmax, HvalX,
                                 HvalY, HvalZ, HvalT);
    // del vector<int> t_GIDs,t_EIDs;
    // del m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, direction, TvalT,
    // t_GIDs, t_EIDs);
    vector<int> t_GIDs, t_EIDs;
    vector<int> t_h_GIDs, t_h_EIDs;
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, HvalT, t_h_GIDs,
                                   t_h_EIDs);
    mergeBreakPtsAdv(HvalT, SvalT, TvalT, t_h_GIDs, t_h_EIDs, t_GIDs, t_EIDs);
    // mergeBreakPts( HvalT, SvalT, TvalT);

    // BPTS
    // printNekArray(HvalT,0);
    // printNekArray(SvalT,0);
    // printNekArray(TvalT,0);

    // Loop forsumation.

    NekDouble sum = 0.0;
    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.

    // cout << "////////Working till here////////" <<endl;
    int quadratureOfMesh =
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(0) +
        m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(1);
    //	m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(2);

    int quad_npoints = ceil((m_order + quadratureOfMesh + 1.0) / 2.0);
    // int quad_npoints = ceil((m_order+3)/2) +
    // m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetTotPoints()+2; cout <<
    // "quad_npoints" << quad_npoints << endl;
    // LibUtilities::PointsKey quadPointsKey(quad_npoints,
    //						m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetPointsType(0));
    LibUtilities::PointsKey quadPointsKey(
        quad_npoints, LibUtilities::PointsType::eGaussGaussLegendre);

    Array<OneD, NekDouble> quad_points =
        LibUtilities::PointsManager()[quadPointsKey]->GetZ();
    Array<OneD, NekDouble> quad_weights =
        LibUtilities::PointsManager()[quadPointsKey]->GetW();
    Array<OneD, NekDouble> t_quad(quad_npoints), t_x(quad_npoints),
        t_y(quad_npoints), t_z(quad_npoints), t_quad_vals(quad_npoints),
        t_xyz_vals(quad_npoints);
    //	vector<int> t_GIDs,t_EIDs;
    //	m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, direction, TvalT, t_GIDs,
    // t_EIDs);

    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        // Create a seg element
        // Actual quadrature points needed are for order = polyorder+Bspline
        // order.
        NekDouble a = TvalT[i];
        NekDouble b = TvalT[i + 1];
        int gID     = t_GIDs[i];
        int eID     = t_EIDs[i];
        for (int j = 0; j < quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (quad_points[j] + 1.0) / 2.0 + a;
            t_x[j]    = PtsX + t_quad[j] * direction[0];
            t_y[j]    = PtsY + t_quad[j] * direction[1];
            t_z[j]    = PtsZ + t_quad[j] * direction[2];
        }
        // Evaluate filter  at these points.
        if (b_symMesh)
        {
            m_siacFilterPtrs[m_SymID]->EvaluateFilter(t_quad, t_quad_vals,
                                                      meshSpacing);
        }
        else
        {
            m_siacFilterPtrs[m_OneID]->EvaluateFilter(
                t_quad, t_quad_vals, meshSpacing, meshTShift, true);
        }
        m_meshHandlePtr->EvaluateAt(t_x, t_y, t_z, gID, eID, t_xyz_vals,
                                    varNum);
        NekDouble integral = 0;
        // integral
        for (int k = 0; k < quad_npoints; k++)
        {
            // integral += 0.5*quad_weights[k] * t_quad_vals[k] *(b-a); // *
            // t_xyz_vals[k];
            integral += 0.5 * quad_weights[k] * t_quad_vals[k] *
                        std::abs((b - a)) * t_xyz_vals[k];
        }
        sum += integral;
    }
    valX = sum;
    return true;
}

} // namespace LSIAC
} // namespace Nektar
