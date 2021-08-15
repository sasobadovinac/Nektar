////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessBodyFittedVelocity.h
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
//  Description: Convert velocity components into the body-fitted coordinate.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_PROCESSBODYFITTEDVELOCITY
#define FIELDUTILS_PROCESSBODYFITTEDVELOCITY

#include "ProcessBoundaryExtract.h"

namespace Nektar
{
namespace FieldUtils
{

/**
 * @brief This processing module calculates the wall shear stress and adds it
 * as an extra-field to the output file, and writes it to a surface output file.
 */
class ProcessBodyFittedVelocity : public ProcessBoundaryExtract
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<ProcessBodyFittedVelocity>::AllocateSharedPtr(f);
    }
    static ModuleKey className;

    ProcessBodyFittedVelocity(FieldSharedPtr f);
    virtual ~ProcessBodyFittedVelocity();

    /// Write mesh to output file.
    virtual void Process(po::variables_map &vm);

    virtual std::string GetModuleName()
    {
        return "ProcessBodyFittedVelocity";
    }

    virtual std::string GetModuleDescription()
    {
        return "Get the wall-normal data at a given origin.";
    }

protected:

private:
    int m_spacedim;

    NekDouble PntToBndElmtPntDistance(
        const Array<OneD, Array<OneD, NekDouble> > & pts,
        const int pId,
        const Array<OneD, Array<OneD, NekDouble> > & bndPts);

    bool LocCoordForNearestPntOnBndElmt_2D(
        const Array<OneD, const NekDouble > & inGloCoord,
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array<OneD, Array<OneD, NekDouble> > & pts,
        Array<OneD, NekDouble > & locCoord,
        Array<OneD, NekDouble > & gloCoord,
        NekDouble & dist,
        const NekDouble iterTol = 1.0e-6,
        const int iterMax = 51);

    bool LocCoordForNearestPntOnBndElmt(
        const Array<OneD, const NekDouble > & inGloCoord,
        SpatialDomains::GeometrySharedPtr bndGeom,
        Array<OneD, NekDouble > & locCoord,
        Array<OneD, NekDouble > & gloCoord,
        NekDouble & dist,
        const NekDouble iterTol = 1.0e-6,
        const int iterMax = 51);
    
    void ScaledCrosssProduct(
        const Array<OneD, NekDouble > & vec1,
        const Array<OneD, NekDouble > & vec2,
        Array<OneD, NekDouble > & vec3);

    //===========================================================
    

    /**
     * @brief Use iteration to get the locCoord. This routine should be used after
     *        we have checked the projected point is inside the projected element.
     * @param bndGeom      Geometry to get the xmap.
     * @param gloCoord     Global coordinate of the point. size=3.
     * @param pts          Global coordinate of the vertices of the elmt. size=2/3.
     * @param dieUse       The main direction(s) used to compute local coordinate
     * @param locCoord     Iteration results for local coordinate(s) 
     * @param dist         Returned distance in physical space if the collapsed 
     *                     locCoord is out of range [-1,1].
     * @param iterTol      Tolerence for iteration.
     * @param iterMax      Maximum iteration steps
     * @return             Converged (true) or not (false)
     */
    bool BisectionForLocCoordOnBndElmt(
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array<OneD, const NekDouble > & gloCoord,
        const Array<OneD, const Array<OneD, NekDouble> > & pts,
        const Array<OneD, const int > & dirUse,
        Array<OneD, NekDouble > & locCoord,
        const NekDouble iterTol = 1.0e-8,
        const int iterMax = 51);
    
    bool NewtonIterForLocCoordOnBndElmt(
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array<OneD, const NekDouble> & gloCoord,
        const Array<OneD, const Array<OneD, NekDouble> > & pts,
        const Array<OneD, const int > & dirUse,
        Array<OneD, NekDouble> & locCoord,
        NekDouble & dist,
        const NekDouble iterTol = 1.0e-8,
        const int iterMax = 51);
    
    /**
    * @brief Check if a point can be projected onto an oundary element in a given
    *        direction. If yes, give the local coordinates of the projected point.
    *        we have checked the projected point is inside the projected element.
    * @param bndGeom      Pointer to the geometry of the boundary element.
    * @param gloCoord     Global coordinate of the point. size=3.
    * @param projDir      Projection direction, which is used as the reference
    *                     direction in the 3D routine. size=3, norm=1. 
    * @param locCoord     Iteration results for local coordinates (if inside).
    * @param projDist     Projection distance betweem the point to the wall point.
    * @param maxDist      Disntance to check if the wall point is desired.
    * @param iterTol      Tolerence for iteration.
    * @return             Inside (true) or not (false)
    */
    /*
    bool BndElmtContainsPoint(
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array<OneD, const NekDouble > & gloCoord,
        const Array<OneD, const NekDouble > & projDir,
        Array< OneD, NekDouble > & locCoord,
        NekDouble & projDist,
        const NekDouble maxDist = 1.0,
        const NekDouble iterTol = 1.0e-8);
    */


};
}
}

#endif
