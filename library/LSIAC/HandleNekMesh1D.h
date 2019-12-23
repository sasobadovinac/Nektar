///////////////////////////////////////////////////////////////////////////////
//
// File HandleNekMesh1D.h
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
// Description: HandleNekMesh1D definition
//
///////////////////////////////////////////////////////////////////////////////
#pragma once

#include "HandleNekMesh.h"
/// Handles Nektar 1D meshes.
/**
 */
namespace Nektar
{
namespace LSIAC
{

class HandleNekMesh1D : public HandleNekMesh
{
public:
    HandleNekMesh1D(FieldUtils::FieldSharedPtr fldSharedPtr)
        : HandleNekMesh(fldSharedPtr)
    {
    }

    HandleNekMesh1D(LibUtilities::SessionReaderSharedPtr sessionPtr)
        : HandleNekMesh(sessionPtr)
    {
    }

protected:
    virtual bool v_GetKnotVec(const int degree,
                              const Array<OneD, NekDouble> &coord,
                              const Array<OneD, NekDouble> &direction,
                              Array<OneD, NekDouble> &knotVec,
                              NekDouble &shift);

    virtual bool v_Get1DVec(vector<NekDouble> &coords);

    virtual bool v_GetBreakPts(
        const NekDouble xcen_offset, const NekDouble ycen_offset,
        const NekDouble zcen_offset, const Array<OneD, NekDouble> &direction,
        const NekDouble t_offset_min, const NekDouble t_offset_max,
        vector<NekDouble> &xPos, vector<NekDouble> &yPos,
        vector<NekDouble> &zPos, vector<NekDouble> &tPos);

    virtual bool v_GetBreakPts_Without_Tmin_Tmax(
        const NekDouble xcen_offset, const NekDouble ycen_offset,
        const NekDouble zcen_offset, const Array<OneD, NekDouble> &direction,
        const NekDouble t_offset_min, const NekDouble t_offset_max,
        vector<NekDouble> &xPos, vector<NekDouble> &yPos,
        vector<NekDouble> &zPos, vector<NekDouble> &tPos);

    virtual bool v_CanTRangebeApplied(
        const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
        const NekDouble scaling, const NekDouble tmin, const NekDouble tmax,
        NekDouble &tminUpdate, NekDouble &tmaxUpdate);

    virtual bool v_CanTRangebeApplied(
        const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
        const Array<OneD, NekDouble> &direction, const NekDouble tmin,
        const NekDouble tmax, NekDouble &meshShift);

    virtual bool v_CanTRangebeAppliedWOMeshShift(
        const NekDouble PtsX, const NekDouble PtsY, const NekDouble PtsZ,
        const Array<OneD, NekDouble> &direction, const NekDouble tmin,
        const NekDouble tmax);

    virtual bool v_LoadData(string Filename, vector<string> &variables);

    virtual bool v_LoadMesh(string var);

    virtual bool v_EvaluateAt(const Array<OneD, NekDouble> &xPos,
                              const Array<OneD, NekDouble> &yPos,
                              const Array<OneD, NekDouble> &zPos, const int gID,
                              const int eID, Array<OneD, NekDouble> &values,
                              int varNum = 0);

    virtual bool v_EvaluateAt(const NekDouble xPos, const NekDouble yPos,
                              const NekDouble zPos, int gID, int eID,
                              NekDouble &value, int varNum = 0);

    virtual bool v_GetListOfGIDs(const NekDouble xPos, const NekDouble yPos,
                                 const NekDouble zPos,
                                 const Array<OneD, NekDouble> &direction,
                                 const vector<NekDouble> t_breaks,
                                 vector<int> &t_GIDs,
                                 vector<int> &t_EIDs) const;
    virtual void v_LoadExpListIntoRTree()
    {
        NEKERROR(ErrorUtil::efatal, "Not implemented for 1D meshes.");
    }

    virtual NekDouble v_GetLargestEdgeLength(const int eid);

    virtual NekDouble v_GetDynamicScaling(Array<OneD, NekDouble> glCoord,
                                          int eid = -1, NekDouble mu = 1.0);

    virtual bool v_InitializeMetricTensor();
    virtual bool v_GetMTScalingOfGIDs(vector<int> &t_GIDs,
                                      Array<OneD, NekDouble> &direction,
                                      vector<NekDouble> &scalings);
};
} // namespace LSIAC
} // namespace Nektar
