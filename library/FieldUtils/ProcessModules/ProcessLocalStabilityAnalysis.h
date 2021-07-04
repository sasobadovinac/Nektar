////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessLocalStabilityAnalysis.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Local linear stability analysis of compressible flow (for now).
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_PROCESSLOCALSTABILITYANALYSIS
#define FIELDUTILS_PROCESSLOCALSTABILITYANALYSIS

#include "ProcessBoundaryExtract.h"

//#include "/disk_two/Nek_Test/nektar++/build_f90/dist/include/nektar++/ThirdParty/CoPSE3d_LST_V1.h"


namespace Nektar
{
namespace FieldUtils
{

/**
 * @brief This processing module calculates the wall shear stress and adds it
 * as an extra-field to the output file, and writes it to a surface output file.
 */

// Define the 

extern "C" 
{
    void F77NAME(copse3d) ();
    void F77NAME(helloworld) ();
}


class ProcessLocalStabilityAnalysis : public ProcessBoundaryExtract
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<ProcessLocalStabilityAnalysis>::AllocateSharedPtr(f);
    }
    static ModuleKey className;

    ProcessLocalStabilityAnalysis(FieldSharedPtr f);
    virtual ~ProcessLocalStabilityAnalysis();

    /// Write mesh to output file.
    virtual void Process(po::variables_map &vm);

    virtual std::string GetModuleName()
    {
        return "ProcessLocalStabilityAnalysis";
    }

    virtual std::string GetModuleDescription()
    {
        return "Calculating wall shear stress";
    }

    // test
    inline void call_hello()
    {
        F77NAME(helloworld) ();
    }


    inline void call_lst()
    {
        F77NAME(copse3d) ();
    }

    // necessasy routines
    bool isInProjectedArea2D(
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array<OneD, const NekDouble > & gloCoord,
        const int projDir);
    
    bool isInProjectedArea3D(
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array<OneD, const NekDouble > & gloCoord,
        const int projDir);

    bool BisectionForLocCoordOnBndElmt(
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array<OneD, const NekDouble > & gloCoord,
        const Array<OneD, const Array<OneD, NekDouble> > & pts,
        const int projDir,
        Array< OneD, NekDouble > & locCoord,
        NekDouble & dist);

    bool NewtonIterationForLocCoordOnBndElmt(
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array<OneD, const NekDouble> &coords,
        const Array<OneD, const Array<OneD, NekDouble> > &pts,
        const int projDir,
        Array<OneD, NekDouble> &Lcoords,
        NekDouble &dist);

    bool BndElmtContainsPoint(
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array< OneD, const NekDouble > & gloCoord,
        Array< OneD, NekDouble > & locCoord,
        const bool isUseY, 
        const NekDouble geomTol,
        NekDouble & dist);

    void GetNormals(
        SpatialDomains::GeometrySharedPtr bndGeom,
        Array< OneD, NekDouble > & locCoord, 
        Array< OneD, NekDouble > & normals);

protected:

private:
    int m_spacedim;

};
}
}

#endif
