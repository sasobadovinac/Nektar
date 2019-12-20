///////////////////////////////////////////////////////////////////////////////
//
// File SmoothieSIAC3D.cpp
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
// Description: SmoothieSIAC3D definition
//
///////////////////////////////////////////////////////////////////////////////
#include "SmoothieSIAC3D.h"
#include "OneSidedSIAC.h"
#include <LibUtilities/Foundations/ManagerAccess.h> // for Points Manager, etc
#include <boost/timer.hpp>

namespace Nektar
{
namespace LSIAC
{
SmoothieSIAC3D::SmoothieSIAC3D(const FilterType filter,
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
        case eSYM_UNEVEN_2kp1:
            m_siacFilterPtrs.emplace_back(new NonSymmetricSIAC(order));
            m_OneID = -1;
            break;
        default:
            cout << "Not implemented yet " << endl;
    }
}

bool SmoothieSIAC3D::v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                                  const NekDouble PtsZ, NekDouble &valX,
                                  NekDouble &valY, NekDouble &valZ)
{
    Array<OneD, NekDouble> direction(3, 0.0);
    direction[0] = 1.0;
    return EvaluateAt(PtsX, PtsY, PtsZ, valX, valY, valZ, direction,
                      m_meshSpacing);
}

bool SmoothieSIAC3D::v_EvaluateAt(const NekDouble PtsX, const NekDouble PtsY,
                                  const NekDouble PtsZ, NekDouble &valX,
                                  NekDouble &valY, NekDouble &valZ,
                                  Array<OneD, NekDouble> &direction,
                                  NekDouble meshSpacing, int varNum)
{
    boost::ignore_unused(valY, valZ);
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
            cout << "Symmetric filter does not fit and OneSided is not defined "
                    "here."
                 << endl;
            cout << "Ideally keep the original value. -1 is returned at end "
                 << endl;
            valX = -1;
            return false;
        }
    }

    m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ, direction, tmin, tmax, HvalX,
                                 HvalY, HvalZ, HvalT);
    mergeBreakPts(HvalT, SvalT, TvalT);

    // BPTS
    // printNekArray(HvalT,0);
    // printNekArray(SvalT,0);
    // printNekArray(TvalT,0);

    // Loop forsumation.

    NekDouble sum = 0.0;
    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.

    // cout << "////////Working till here////////" <<endl;
    if (m_calculateQuadrature)
    {
        // int quadratureOfMesh =
        // std::min(m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetTotPoints()+2,21);
        // int quadratureOfMesh =
        // std::min(m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetNcoeffs()+2,21);
        int quadratureOfMesh =
            m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(0) +
            m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(1) +
            m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(2);

        m_quad_npoints = ceil((m_order + quadratureOfMesh + 1) / 2);

        // cout << "m_quad_npoints" << m_quad_npoints << endl;
        // LibUtilities::PointsKey quadPointsKey(m_quad_npoints,
        //					m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetPointsType(0));
        LibUtilities::PointsKey quadPointsKey(
            m_quad_npoints, LibUtilities::PointsType::eGaussGaussLegendre);

        m_quad_points  = LibUtilities::PointsManager()[quadPointsKey]->GetZ();
        m_quad_weights = LibUtilities::PointsManager()[quadPointsKey]->GetW();
        cout << "came here into SmoSI3d:line 154 quadpoints are "
             << m_quad_npoints << "\tNumModes\t" << quadratureOfMesh << endl;
        m_calculateQuadrature = false;
    }

    Array<OneD, NekDouble> t_quad(m_quad_npoints), t_x(m_quad_npoints),
        t_y(m_quad_npoints), t_z(m_quad_npoints), t_quad_vals(m_quad_npoints),
        t_xyz_vals(m_quad_npoints);
    vector<int> t_GIDs, t_EIDs;
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, TvalT, t_GIDs,
                                   t_EIDs);

    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        // Create a seg element
        // Actual quadrature points needed are for order = polyorder+Bspline
        // order.
        NekDouble a = TvalT[i];
        NekDouble b = TvalT[i + 1];
        int gID     = t_GIDs[i];
        int eID     = t_EIDs[i];
        for (int j = 0; j < m_quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (m_quad_points[j] + 1.0) / 2.0 + a;
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
        for (int k = 0; k < m_quad_npoints; k++)
        {
            // integral += 0.5*quad_weights[k] * t_quad_vals[k] *(b-a); // *
            // t_xyz_vals[k];
            integral += 0.5 * m_quad_weights[k] * t_quad_vals[k] *
                        std::abs((b - a)) * t_xyz_vals[k];
        }
        sum += integral;
    }
    valX = sum;
    return true;
}
/*
bool SmoothieSIAC3D::v_EvaluateRecursiveAt( const NekDouble PtsX, const
NekDouble PtsY, const NekDouble PtsZ, NekDouble &valX, NekDouble &valY,
NekDouble &valZ, vector<Array<OneD,NekDouble>> directions, const
vector<NekDouble>& meshSpacings, const vector<int>& varNums)
{
        return true;
}
*/

bool SmoothieSIAC3D::v_EvaluateRecursiveAt(
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
    // cout <<"Recursive depthi:\t" << curLevel << endl;
    if (totLevels - 1 == curLevel)
    {
        //	cout << "Last Recursive Step :\t"<< curLevel << endl;
        // This is equivalent to normal evaluateAt. So Directly call Evaluate
        // At.
        return Sms[curLevel]->EvaluateAt(
            PtsX, PtsY, PtsZ, valX, valY, valZ, directions[curLevel],
            meshSpacings[curLevel], varNums[curLevel]);
    }

    // cout << "came out here successfully.Now get to work" << endl;
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
            cout << "Symmetric filter does not fit and OneSided is not defined "
                    "here."
                 << endl;
            cout << "Ideally keep the original value. -1 is returned at end "
                 << endl;
            valX = -1;
            return false;
        }
    }

    m_meshHandlePtr->GetBreakPts(PtsX, PtsY, PtsZ, direction, tmin, tmax, HvalX,
                                 HvalY, HvalZ, HvalT);
    mergeBreakPts(HvalT, SvalT, TvalT);

    // BPTS
    // printNekArray(HvalT,0);
    // printNekArray(SvalT,0);
    // printNekArray(TvalT,0);

    // Loop forsumation.

    NekDouble sum = 0.0;
    // For getting the local gauss quadrature use the getCoords fo StdRegions
    // objects.

    // cout << "////////Working till here////////" <<endl;
    if (m_calculateQuadrature)
    {
        // int quadratureOfMesh =
        // std::min(m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetTotPoints()+2,21);
        // int quadratureOfMesh =
        // std::min(m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetNcoeffs()+2,21);
        int quadratureOfMesh =
            m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(0) +
            m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(1) +
            m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetBasisNumModes(2);

        m_quad_npoints = ceil((m_order + quadratureOfMesh + 1.0) / 2.0);

        // cout << "m_quad_npoints" << m_quad_npoints << endl;
        // LibUtilities::PointsKey quadPointsKey(m_quad_npoints,
        //					m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetPointsType(0));
        LibUtilities::PointsKey quadPointsKey(
            m_quad_npoints, LibUtilities::PointsType::eGaussGaussLegendre);

        m_quad_points  = LibUtilities::PointsManager()[quadPointsKey]->GetZ();
        m_quad_weights = LibUtilities::PointsManager()[quadPointsKey]->GetW();
        cout << "came here into SmoSI3d:line 154 quadpoints are "
             << m_quad_npoints << "\tNumModes\t" << quadratureOfMesh << endl;
        m_calculateQuadrature = false;
    }

    Array<OneD, NekDouble> t_quad(m_quad_npoints), t_x(m_quad_npoints),
        t_y(m_quad_npoints), t_z(m_quad_npoints), t_quad_vals(m_quad_npoints),
        t_xyz_vals(m_quad_npoints);
    vector<int> t_GIDs, t_EIDs;
    m_meshHandlePtr->GetListOfGIDs(PtsX, PtsY, PtsZ, direction, TvalT, t_GIDs,
                                   t_EIDs);

    for (int i = 0; i < TvalT.size() - 1; i++)
    {
        // Create a seg element
        // Actual quadrature points needed are for order = polyorder+Bspline
        // order.
        NekDouble a = TvalT[i];
        NekDouble b = TvalT[i + 1];
        for (int j = 0; j < m_quad_npoints; j++)
        {
            t_quad[j] = (b - a) * (m_quad_points[j] + 1.0) / 2.0 + a;
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
        for (int j = 0; j < m_quad_npoints; j++)
        {
            Sms[curLevel]->EvaluateRecursiveAt(
                t_x[j], t_y[j], t_z[j], t_xyz_vals[j], valY, valZ, Sms,
                directions, meshSpacings, varNums, curLevel + 1);
        }

        NekDouble integral = 0;
        // integral
        for (int k = 0; k < m_quad_npoints; k++)
        {
            // integral += 0.5*quad_weights[k] * t_quad_vals[k] *(b-a); // *
            // t_xyz_vals[k];
            integral += 0.5 * m_quad_weights[k] * t_quad_vals[k] *
                        std::abs((b - a)) * t_xyz_vals[k];
        }
        sum += integral;
    }
    valX = sum;
    return true;
}

/*
bool SmoothieSIAC3D::v_SetupLineForLSIAC( const Array<OneD,NekDouble>
&direction, const vector<NekDouble> &stPoint, const NekDouble tmin, const
NekDouble tmax, const int n_quadPts, vector<NekDouble> &HvalT, vector<int>
&t_GIDs, vector<int> &t_EIDs, Array<OneD,NekDouble> &t_LineElm )
{
        //1) Find all the elements on the line segment.
        //2) Find all the break points for these elements.
        //3) Get Elem Ids.
        //4) Quadrature points.
        vector<NekDouble> HvalX, HvalY, HvalZ ;
        m_meshHandlePtr->GetBreakPts(stPoint[0], stPoint[1],
stPoint[2],direction, tmin,tmax, HvalX, HvalY, HvalZ, HvalT);
        //printNekArray(HvalT,0.0);
        m_meshHandlePtr->GetListOfGIDs(stPoint[0], stPoint[1], stPoint[2],
direction, HvalT, t_GIDs, t_EIDs);
        // Number of quadrature points should be input.
        //int n_quadPts = 4;
        LibUtilities::PointsKey
quadPointsKey(n_quadPts,Nektar::LibUtilities::eGaussGaussLegendre );
        Array<OneD,NekDouble> quad_points =
LibUtilities::PointsManager()[quadPointsKey]->GetZ(); Array<OneD,NekDouble>
quad_weights  = LibUtilities::PointsManager()[quadPointsKey]->GetW();

        int n_meshT_quadPts = n_quadPts* (t_EIDs.size()-1);
        t_LineElm = Array<OneD,NekDouble> (n_meshT_quadPts);
        for ( int i =0; i < t_EIDs.size()-1; i++)
        {
                for ( int j =0; j< n_quadPts; j++)
                {
                        t_LineElm[i*n_quadPts+j] = HvalT[i] +
(quad_points[j]+1.0) *(HvalT[i+1] - HvalT[i])/2.0;
                }
        }
        return true	;
}

bool SmoothieSIAC3D::v_SetupLineForLSIAC_ReSamp( const Array<OneD,NekDouble>
&direction, const vector<NekDouble> &stPoint, const NekDouble tmin, const
NekDouble tmax, const int n_quadPts, const int n_quadPts_Resample,
                                                                        vector<NekDouble>
&HvalT, vector<int> &t_GIDs, vector<int> &t_EIDs, Array<OneD,NekDouble>
&t_LineElm, Array<OneD,NekDouble> &t_LineElm_Resample )
{
        //1) Find all the elements on the line segment.
        //2) Find all the break points for these elements.
        //3) Get Elem Ids.
        //4) Quadrature points.
        vector<NekDouble> HvalX, HvalY, HvalZ ;
        m_meshHandlePtr->GetBreakPts(stPoint[0], stPoint[1],
stPoint[2],direction, tmin,tmax, HvalX, HvalY, HvalZ, HvalT);
        //printNekArray(HvalT,0.0);
        m_meshHandlePtr->GetListOfGIDs(stPoint[0], stPoint[1], stPoint[2],
direction, HvalT, t_GIDs, t_EIDs);
        // Number of quadrature points should be input.
        //int n_quadPts = 4;
        {
        LibUtilities::PointsKey
quadPointsKey(n_quadPts,Nektar::LibUtilities::eGaussGaussLegendre );
        Array<OneD,NekDouble> quad_points =
LibUtilities::PointsManager()[quadPointsKey]->GetZ(); Array<OneD,NekDouble>
quad_weights  = LibUtilities::PointsManager()[quadPointsKey]->GetW(); int
n_meshT_quadPts = n_quadPts* (t_EIDs.size()-1); t_LineElm =
Array<OneD,NekDouble> (n_meshT_quadPts); for ( int i =0; i < t_EIDs.size()-1;
i++)
        {
                for ( int j =0; j< n_quadPts; j++)
                {
                        t_LineElm[i*n_quadPts+j] = HvalT[i] +
(quad_points[j]+1.0) *(HvalT[i+1] - HvalT[i])/2.0;
                }
        }
        }
        {
        //int n_quadPts_Resample = 4;
        LibUtilities::PointsKey
quadPointsKey(n_quadPts_Resample,Nektar::LibUtilities::eGaussGaussLegendre );
        Array<OneD,NekDouble> quad_points =
LibUtilities::PointsManager()[quadPointsKey]->GetZ(); Array<OneD,NekDouble>
quad_weights  = LibUtilities::PointsManager()[quadPointsKey]->GetW(); int
n_meshT_quadPts_Resample = n_quadPts_Resample* (t_EIDs.size()-1);
        t_LineElm_Resample = Array<OneD,NekDouble> (n_meshT_quadPts_Resample);
        for ( int i =0; i < t_EIDs.size()-1; i++)
        {
                for ( int j =0; j< n_quadPts_Resample; j++)
                {
                        t_LineElm_Resample[i*n_quadPts_Resample+j] = HvalT[i] +
(quad_points[j]+1.0) *(HvalT[i+1] - HvalT[i])/2.0;
                }
        }
        }

        return true	;
}

bool SmoothieSIAC3D::v_GetVLineForLSIAC(	const int n_quadPts, const
vector<NekDouble> &stPoint, const Array<OneD,NekDouble> &direction, const
vector<NekDouble> &HvalT, const vector<int> &t_GIDs, const vector<int> &t_EIDs,
                                                        const
Array<OneD,NekDouble> &t_LineElm, Array<OneD,NekDouble> tv_LineElm, int varNum )
{
        // loop through elements.
                // Find values for each quadrature point.
        // create tx_LineElm,ty_LineElm, tz_LineElm);
        Array<OneD,NekDouble> tx_LineElm( n_quadPts*(t_EIDs.size()-1));
        Array<OneD,NekDouble> ty_LineElm( n_quadPts*(t_EIDs.size()-1));
        Array<OneD,NekDouble> tz_LineElm( n_quadPts*(t_EIDs.size()-1));

        for( int i =0; i< t_LineElm.num_elements(); i++)
        {
                tx_LineElm[i] = stPoint[0] + direction[0] * t_LineElm[i];
                ty_LineElm[i] = stPoint[1] + direction[1] * t_LineElm[i];
                tz_LineElm[i] = stPoint[2] + direction[2] * t_LineElm[i];
        }

        Array<OneD,NekDouble> elx,ely,elz;
        Array<OneD,NekDouble> elv(n_quadPts,0.0);
        Array<OneD,NekDouble> glCoord(3);
        for ( int i =0; i < t_EIDs.size()-1; i++)
        {
                elx = tx_LineElm.CreateWithOffset( tx_LineElm, i*n_quadPts);
                ely = ty_LineElm.CreateWithOffset( ty_LineElm, i*n_quadPts);
                elz = tz_LineElm.CreateWithOffset( tz_LineElm, i*n_quadPts);
        //	elv = tz_LineElm.CreateWithOffset( tv_LineElm, i*n_quadPts);
                m_meshHandlePtr->EvaluateAt( elx,ely,elz, t_GIDs[i], t_EIDs[i],
elv, varNum); memcpy( &tv_LineElm[i*n_quadPts], &elv[0],
sizeof(NekDouble)*n_quadPts);
        }
        return true;
}

bool SmoothieSIAC3D::v_EvaluateUsingLineAt( const vector<NekDouble> &stPoint,
const Array<OneD,NekDouble> &direction, const int n_quadPts, const NekDouble
meshSpacing, const vector<NekDouble> &tparams, vector<NekDouble> &tvals, int
varNum)
{
        //boost::timer tim1, timTotal;
        // set up line points.
        // evaluate funtion values at these points.
        NekDouble t_mesh_min = *(std::min_element( tparams.begin(),
tparams.end())); NekDouble t_mesh_max = *(std::max_element( tparams.begin(),
tparams.end())); vector<NekDouble> HvalT; vector<int> Gids, Eids;
        Array<OneD,NekDouble> t_LineElm;

        // Set up the quadrature.
        SetupLineForLSIAC( direction, stPoint, t_mesh_min, t_mesh_max,n_quadPts,
HvalT, Gids, Eids, t_LineElm);
        //double el1 = tim1.elapsed();

        //boost::timer tim2;
        Array<OneD,NekDouble> tv_LineElm(t_LineElm.num_elements());
        // Get the values of quadrature.
        GetVLineForLSIAC( n_quadPts, stPoint, direction, HvalT, Gids, Eids,
t_LineElm, tv_LineElm, varNum);
        //double el2 = tim2.elapsed();

        // Set up one more array for quadrature weights or jacobian values.
                //tW_LineElm
                //tJ_LineElm
        //NekDouble meshSpacing = 0.1; // NekDouble.
        //boost::timer tim3;
        NekDouble t_min, t_max,meshTShift;
        LibUtilities::PointsKey
quadPointsKey(n_quadPts,Nektar::LibUtilities::eGaussGaussLegendre );
        //Array<OneD,NekDouble> quad_points =
LibUtilities::PointsManager()[quadPointsKey]->GetZ(); Array<OneD,NekDouble>
t_quad_weights  = LibUtilities::PointsManager()[quadPointsKey]->GetW();
        Array<OneD,NekDouble> t_quadPts(n_quadPts), t_quad_vals(n_quadPts);

        for (int i =0; i < tparams.size(); i++)
        {
                // Figuring out which type of filter need to be used.
                bool b_symMesh = true;
                NekDouble t = tparams[i];
                m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, t_min,
t_max ); if (t + t_min < t_mesh_min && t+ t_max > t_mesh_max)
                {
                        assert("No filter can be applied");
                        return false;
                }
                if (t + t_min< t_mesh_min)
                { // one sided filter need to be applied.
                        m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing,t_min,t_max);
                        meshTShift = t_mesh_min - t - t_min;
                        b_symMesh= false;
                        m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing,t_min,t_max,meshTShift);
                }
                if (t + t_max > t_mesh_max)
                { // one sided filter need to be applied.
                        m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing,t_min,t_max);
                        meshTShift = t_mesh_max - t -t_max;
                        b_symMesh = false;
                        m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing,t_min,t_max,meshTShift);
                }

                //cout << "for parameter t " << t << "\t tmin: \t" << t_min <<
"\t tmax: \t" << t_max << endl;

                assert( t+ t_min +TOLERENCE > t_mesh_min && "Above code should
have fixed this issue" ); assert( t+ t_max -TOLERENCE < t_mesh_max && "Above
code should have fixed this issue" );
                // if symmetric mesh. Could be true for both.
                int startElmIndex, endElmIndex;
                // start ElementId and End Element Id;
                for ( int j= 0; j < HvalT.size(); j++)
                {
                        if( HvalT[j] > t + t_min+ TOLERENCE)
                        {
                                startElmIndex = j-1;
                                break;
                        }
                }

                for ( int j= 0; j < HvalT.size(); j++)
                {
                        if( HvalT[j] > t + t_max- TOLERENCE)
                        {
                                endElmIndex = j-1;
                                break;
                        }
                }

                //cout << "for parameter t " << t << "\t sElIndex: \t" <<
startElmIndex<< "\t tmax: \t" << endElmIndex<< endl;


                assert( startElmIndex <= endElmIndex && "Wrong");
                NekDouble sum= 0.0;
                for (int el = startElmIndex; el<=endElmIndex; el++)
                {
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
tv_LineElm[el*n_quadPts+qp]* std::abs(HvalT[el+1]-HvalT[el])* t_quad_vals[qp];
                        }
                        sum += integral;
                }
                tvals.push_back(sum);
        }
        //double el3 = tim3.elapsed();
        //double elT = timTotal.elapsed();
        //cout << "part1\t"<< el1 << "\t\t" << el2 << "\t\t" << el3 <<
"\tTotal\t" << elT << endl; return true;
}

bool SmoothieSIAC3D::v_EvaluateUsingLineAt_v2( const vector<NekDouble> &stPoint,
const Array<OneD,NekDouble> &direction, const int n_quadPts, const int
n_quadPts_resample, const NekDouble meshSpacing, const vector<NekDouble>
&tparams, vector<NekDouble> &tvals, int varNum)
{

        tvals.clear();
        boost::timer tim1, timTotal;
        // set up line points.
        // evaluate funtion values at these points.
        NekDouble t_mesh_min = *(std::min_element( tparams.begin(),
tparams.end())); NekDouble t_mesh_max = *(std::max_element( tparams.begin(),
tparams.end())); vector<NekDouble> HvalT; vector<int> Gids, Eids;
        Array<OneD,NekDouble> t_LineElm;
        Array<OneD,NekDouble> t_LineElm_resample;
        SetupLineForLSIAC_ReSamp( direction, stPoint, t_mesh_min,
t_mesh_max,n_quadPts, n_quadPts_resample, HvalT, Gids, Eids, t_LineElm,
t_LineElm_resample); Array<OneD,NekDouble> tv_LineElm(t_LineElm.num_elements());
        Array<OneD,NekDouble>
tv_LineElm_resample(t_LineElm_resample.num_elements()); GetVLineForLSIAC(
n_quadPts, stPoint, direction, HvalT, Gids, Eids, t_LineElm, tv_LineElm,
varNum); GetVLineForLSIAC_resample( n_quadPts, n_quadPts_resample, stPoint,
direction, HvalT, Gids, Eids, t_LineElm, tv_LineElm, t_LineElm_resample,
tv_LineElm_resample, varNum);

        // Set up the quadrature.
        //SetupLineForLSIAC( direction, stPoint, t_mesh_min,
t_mesh_max,n_quadPts, HvalT, Gids, Eids, t_LineElm); double el1 =
tim1.elapsed();


        boost::timer tim2;
        // Get the values of quadrature.




        // Get values for not extra points but new points.
        // Set up more points for each setup.
                // 1. Can use eariler code to set it up.
        //	SetupLineForLSIAC( direction, stPoint, t_mesh_min,
t_mesh_max,n_quadPts_resample, HvalT, Gids, Eids, t_LineElm_resample);
                // 2. Update values at resample points


        double el2 = tim2.elapsed();

        // Set up one more array for quadrature weights or jacobian values.
                //tW_LineElm
                //tJ_LineElm
        //NekDouble meshSpacing = 0.1; // NekDouble.
        boost::timer tim3;
        NekDouble t_min, t_max,meshTShift;
        LibUtilities::PointsKey
quadPointsKey(n_quadPts_resample,Nektar::LibUtilities::eGaussGaussLegendre );
        //Array<OneD,NekDouble> quad_points =
LibUtilities::PointsManager()[quadPointsKey]->GetZ(); Array<OneD,NekDouble>
t_quad_weights  = LibUtilities::PointsManager()[quadPointsKey]->GetW();
        Array<OneD,NekDouble> t_quadPts(n_quadPts_resample),
t_quad_vals(n_quadPts_resample);

        for (int i =0; i < tparams.size(); i++)
        {
                // Figuring out which type of filter need to be used.
                bool b_symMesh = true;
                NekDouble t = tparams[i];
                m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, t_min,
t_max ); if (t + t_min < t_mesh_min && t+ t_max > t_mesh_max)
                {
                        assert("No filter can be applied");
                        return false;
                }
                if (t + t_min< t_mesh_min)
                { // one sided filter need to be applied.
                        m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing,t_min,t_max);
                        meshTShift = t_mesh_min - t - t_min;
                        b_symMesh= false;
                        m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing,t_min,t_max,meshTShift);
                }
                if (t + t_max > t_mesh_max)
                { // one sided filter need to be applied.
                        m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing,t_min,t_max);
                        meshTShift = t_mesh_max - t -t_max;
                        b_symMesh = false;
                        m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing,t_min,t_max,meshTShift);
                }

                //cout << "for parameter t " << t << "\t tmin: \t" << t_min <<
"\t tmax: \t" << t_max << endl;

                assert( t+ t_min +TOLERENCE > t_mesh_min && "Above code should
have fixed this issue" ); assert( t+ t_max -TOLERENCE < t_mesh_max && "Above
code should have fixed this issue" );
                // if symmetric mesh. Could be true for both.
                int startElmIndex, endElmIndex;
                // start ElementId and End Element Id;
                for ( int j= 0; j < HvalT.size(); j++)
                {
                        if( HvalT[j] > t + t_min+ TOLERENCE)
                        {
                                startElmIndex = j-1;
                                break;
                        }
                }

                for ( int j= 0; j < HvalT.size(); j++)
                {
                        if( HvalT[j] > t + t_max- TOLERENCE)
                        {
                                endElmIndex = j-1;
                                break;
                        }
                }

                //cout << "for parameter t " << t << "\t sElIndex: \t" <<
startElmIndex<< "\t tmax: \t" << endElmIndex<< endl;

                assert( startElmIndex <= endElmIndex && "Wrong");
                NekDouble sum= 0.0;
                for (int el = startElmIndex; el<=endElmIndex; el++)
                {
                        for (int qp = 0; qp < n_quadPts_resample; qp++)
                        {
                                t_quadPts[qp] =
t_LineElm_resample[el*n_quadPts_resample+qp] -t;
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
                        for (int qp = 0; qp < n_quadPts_resample; qp++)
                        {
                                integral += 0.5*t_quad_weights[qp]*
tv_LineElm_resample[el*n_quadPts_resample+qp]* std::abs(HvalT[el+1]-HvalT[el])*
                                                                                        t_quad_vals[qp];
                        }
                        sum += integral;
                }
                tvals.push_back(sum);
        }
        double el3 = tim3.elapsed();
        double elT = timTotal.elapsed();
        //cout << "part1\t"<< el1 << "\t\t" << el2 << "\t\t" << el3 <<
"\tTotal\t" << elT << endl; return true;

}


bool SmoothieSIAC3D::v_EvaluateUsingLineAt_v3( const vector<NekDouble> &stPoint,
const Array<OneD,NekDouble> &direction, const int n_quadPts, const int
n_quadPts_resample, const NekDouble meshSpacing, const vector<NekDouble>
&tparams, vector<NekDouble> &tvals, int varNum)
{
        boost::timer tim1, timTotal;
        // set up line points.
        // evaluate funtion values at these points.
        NekDouble t_mesh_min = *(std::min_element( tparams.begin(),
tparams.end())); NekDouble t_mesh_max = *(std::max_element( tparams.begin(),
tparams.end())); vector<NekDouble> HvalT; vector<int> Gids, Eids;
        Array<OneD,NekDouble> t_LineElm;

        // Set up the quadrature.
        SetupLineForLSIAC( direction, stPoint, t_mesh_min, t_mesh_max,n_quadPts,
HvalT, Gids, Eids, t_LineElm); double el1 = tim1.elapsed();

        boost::timer tim2;
        Array<OneD,NekDouble> tv_LineElm(t_LineElm.num_elements());
        // Get the values of quadrature.
        GetVLineForLSIAC( n_quadPts, stPoint, direction, HvalT, Gids, Eids,
t_LineElm, tv_LineElm, varNum);

        LibUtilities::PointsKey
quadPointsKey(n_quadPts,Nektar::LibUtilities::eGaussGaussLegendre );
        LibUtilities::BasisKey bk =
LibUtilities::BasisKey(Nektar::LibUtilities::eGauss_Lagrange , n_quadPts,
quadPointsKey ) ; Nektar::StdRegions::StdExpansionSharedPtr segExp_std=
MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(bk);

        LibUtilities::PointsKey
quadPointsKey_resample(n_quadPts_resample,Nektar::LibUtilities::eGaussGaussLegendre
); Array<OneD,NekDouble> q_quadR_weights  =
LibUtilities::PointsManager()[quadPointsKey_resample]->GetW();
        Array<OneD,NekDouble> q_quadR_Pts=
LibUtilities::PointsManager()[quadPointsKey_resample]->GetZ();
        Array<OneD,NekDouble> t_quadPts(n_quadPts_resample),
t_quad_vals(n_quadPts_resample);
//	Array<OneD,NekDouble> t_quad(n_quadPts_resample);

        Array<OneD,NekDouble> lcoord(3,0.0);

        for (int ii =0; ii < tparams.size(); ii++)
        {
                // Figuring out which type of filter need to be used.
                bool b_symMesh = true;
                NekDouble t = tparams[ii];
                NekDouble t_min, t_max,meshTShift;
                m_siacFilterPtrs[m_SymID]->GetFilterRange(meshSpacing, t_min,
t_max ); if (t + t_min < t_mesh_min && t+ t_max > t_mesh_max)
                {
                        assert("No filter can be applied");
                        return false;
                }
                if (t + t_min< t_mesh_min)
                { // one sided filter need to be applied.
                        m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing,t_min,t_max);
                        meshTShift = t_mesh_min - t - t_min;
                        b_symMesh= false;
                        m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing,t_min,t_max,meshTShift);
                }
                if (t + t_max > t_mesh_max)
                { // one sided filter need to be applied.
                        m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing,t_min,t_max);
                        meshTShift = t_mesh_max - t -t_max;
                        b_symMesh = false;
                        m_siacFilterPtrs[m_OneID]->GetFilterRange(meshSpacing,t_min,t_max,meshTShift);
                }

                //cout << "for parameter t " << t << "\t tmin: \t" << t_min <<
"\t tmax: \t" << t_max << endl;

                assert( t+ t_min +TOLERENCE > t_mesh_min && "Above code should
have fixed this issue" ); assert( t+ t_max -TOLERENCE < t_mesh_max && "Above
code should have fixed this issue" );

                // if symmetric mesh. Could be true for both.
                int startElmIndex, endElmIndex;
                // start ElementId and End Element Id;
                for ( int j= 0; j < HvalT.size(); j++)
                {
                        if( HvalT[j] > t + t_min+ TOLERENCE)
                        {
                                startElmIndex = j-1;
                                break;
                        }
                }

                for ( int j= 0; j < HvalT.size(); j++)
                {
                        if( HvalT[j] > t + t_max- TOLERENCE)
                        {
                                endElmIndex = j-1;
                                break;
                        }
                }

                //cout << "for parameter t " << t << "\t sElIndex: \t" <<
startElmIndex<< "\t tmax: \t" << endElmIndex<< endl;

                assert( startElmIndex <= endElmIndex && "Wrong");
                NekDouble sum= 0.0;
                vector<NekDouble> SvalT, LvalT,TvalT;
                if ( b_symMesh )
                {
                        m_siacFilterPtrs[m_SymID]->GetBreakPts( meshSpacing,
SvalT);
                }
                else
                {
                        m_siacFilterPtrs[m_OneID]->GetBreakPts( meshSpacing,
SvalT, meshTShift );
                }

                // Insert all mesh breakpoints including min and max using
                // HvalT into LvalT
                LvalT.clear();
                LvalT.push_back(t_min);
                for ( int j= startElmIndex; j <= endElmIndex+1; j++)
                {
                        if ( ( HvalT[j] > t+t_min ) && (HvalT[j] < t+t_max) )
                        {
                                LvalT.push_back(HvalT[j]-t);
                        }
                }
                LvalT.push_back(t_max);
                mergeBreakPts( SvalT, LvalT, TvalT);

                // m_meshHandlePtr->GetListOfGIDs(PtsX,PtsY,PtsZ, diretion,
TvalT, t_GIDs,t_EIDs);
                // Need to be done.
                for (int i =0; i < TvalT.size()-1; i++)
                {
                        NekDouble a = TvalT[i];
                        NekDouble b = TvalT[i+1];
                        //	int gID = t_GIDs[i]; need elmIndex on line.
Hence code below.
                        //	int eID = t_EIDs[i];
                        int elmIndex = -1;
                        for ( int j= startElmIndex; j <= endElmIndex+1; j++)
                        {
                                if ( HvalT[j] > t+(a+b)/2.0 )
                                {
                                        elmIndex = j-1;
                                        break;
                                }
                        }
                        assert( elmIndex != -1 && "Elm Index is wrong" );
                        for ( int j =0; j < n_quadPts_resample; j++)
                        {
                                t_quadPts[j] = (b-a)*(q_quadR_Pts[j]+1.0)/2.0 +
a;
                        }
                        if (b_symMesh)
                        {
                                m_siacFilterPtrs[m_SymID]->EvaluateFilter(t_quadPts,t_quad_vals,
meshSpacing );
                        }
                        else{
                                m_siacFilterPtrs[m_OneID]->EvaluateFilter(t_quadPts,t_quad_vals,
meshSpacing, meshTShift, true);
                        }
                //	m_meshHandlePtr->EvaluateAt(t_xyz_vals);

                        NekDouble integral = 0.0;
                        Array<OneD,NekDouble> elv = tv_LineElm.CreateWithOffset(
tv_LineElm, elmIndex*n_quadPts); for ( int  k=0; k < n_quadPts_resample; k++)
                        {
                                // need to calcuate value of t_xyz_vals.
                                // scaling of HvalT[j-1] to -1 and HvalT[j] to
+1
                                //	lcoord[0] = (t+ t_quadPts[k] -
HvalT[elmIndex-1])/ (HvalT[elmIndex]-HvalT[elmIndex-1])*2.0 -1.0; lcoord[0] =
(t+ t_quadPts[k] - HvalT[elmIndex])/ (HvalT[elmIndex+1]-HvalT[elmIndex])*2.0
-1.0; NekDouble t_xyz_val = segExp_std->PhysEvaluate(lcoord, elv);
                        // SegExpasionPhysEvaluate
                                integral += 0.5*q_quadR_weights[k] *
t_quad_vals[k] *std::abs(b-a) *t_xyz_val;
                        }
                        sum += integral;
                }
                tvals.push_back(sum);


        }

        return true;
}


bool SmoothieSIAC3D::v_GetVLineForLSIAC_resample(	const int n_quadPts,
const int n_quadPts_resample, const vector<NekDouble> &stPoint, const
Array<OneD,NekDouble> &direction, const vector<NekDouble> &HvalT, const
vector<int> &t_GIDs, const vector<int> &t_EIDs, const Array<OneD,NekDouble>
&t_LineElm, const Array<OneD,NekDouble> tv_LineElm, const Array<OneD,NekDouble>
&t_LineElm_resample, Array<OneD,NekDouble> tv_LineElm_resample, int varNum )
{
        // loop through elements.
                // Find values for each quadrature point.
        // create tx_LineElm,ty_LineElm, tz_LineElm);
        Array<OneD,NekDouble> tx_LineElm_resample( n_quadPts*(t_EIDs.size()-1));
        Array<OneD,NekDouble> ty_LineElm_resample( n_quadPts*(t_EIDs.size()-1));
        Array<OneD,NekDouble> tz_LineElm_resample( n_quadPts*(t_EIDs.size()-1));

        for( int i =0; i< t_LineElm.num_elements(); i++)
        {
                tx_LineElm_resample[i] = stPoint[0] + direction[0] *
t_LineElm_resample[i]; ty_LineElm_resample[i] = stPoint[1] + direction[1] *
t_LineElm_resample[i]; tz_LineElm_resample[i] = stPoint[2] + direction[2] *
t_LineElm_resample[i];
        }

        Array<OneD,NekDouble> elx,ely,elz;
        Array<OneD,NekDouble> elv(n_quadPts,0.0);
        Array<OneD,NekDouble> elv_resample(n_quadPts_resample,0.0);
        Array<OneD,NekDouble> glCoord(3);

        // Create a STD element.
                //create Basiskey
        //m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetPointsType(0);
        //Array<OneD,NekDouble> quad_points =
LibUtilities::PointsManager()[quadPointsKey]->GetZ();
        //Array<OneD,NekDouble> quad_weights  =
LibUtilities::PointsManager()[quadPointsKey]->GetW();
        //LibUtilities::PointsKey pk = LibUtilities::PointsKey(
n_quadpts_resample,
        //
m_meshHandlePtr->m_expansions[0]->GetExp(0)->GetPointsType(0) ) ;
        LibUtilities::PointsKey
quadPointsKey(n_quadPts,Nektar::LibUtilities::eGaussGaussLegendre );
        LibUtilities::BasisKey bk =
LibUtilities::BasisKey(Nektar::LibUtilities::eGauss_Lagrange , n_quadPts,
quadPointsKey ) ; Nektar::StdRegions::StdExpansionSharedPtr segExp_std=
MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(bk);

        LibUtilities::PointsKey quadPointsKey_resample(n_quadPts_resample
,Nektar::LibUtilities::eGaussGaussLegendre ); Array<OneD,NekDouble>
quad_points_resample =
LibUtilities::PointsManager()[quadPointsKey_resample]->GetZ();
        Array<OneD,NekDouble> quad_weights_resample  =
LibUtilities::PointsManager()[quadPointsKey_resample]->GetW();

        //BasisKey
        //StdSegExp ();
        // Create new quandrature points.
        // 1. First for loop for each element.
        // 2. Second for loop to evaluate at new quadrature points.

        //cout << "test1" << segExp_std->GetNcoeffs() << "ptstest2" <<
segExp_std->GetTotPoints() << endl;
        //return true;

        Array<OneD,NekDouble> Lcoord(3,0.0);
        for ( int i =0; i < t_EIDs.size()-1; i++)
        {
                // for each element get the values and evaluate at more points.

        //	elx = tx_LineElm.CreateWithOffset( tx_LineElm, i*n_quadPts);
                elv = tv_LineElm.CreateWithOffset( tv_LineElm, i*n_quadPts);
                //memcpy( &elv[0], &tv_LineElm[i*n_quadPts],
sizeof(NekDouble)*n_quadPts); for (int j=0; j < n_quadPts_resample; j++)
                {
                        Lcoord[0] = quad_points_resample[j] ;
                        elv_resample[j] = segExp_std->PhysEvaluate(Lcoord, elv);
                }
                memcpy( &tv_LineElm_resample[i*n_quadPts_resample],
&elv_resample[0],  sizeof(NekDouble)*n_quadPts_resample);
        }


        return true;
}

*/

} // namespace LSIAC
} // namespace Nektar
