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


        /**
     * @brief At each quadrature point inside the domian, compute the body-fitted
     *        coordinate system with respect to the input boundary id.
     * @param targetBndId  Target boundary id.
     * @param assistVec    Unit assistant vector, the cross product of inward-
     *                     pointing wall normalId and wihch gives of the main
     *                     tangential direction of the body-fitted system.
     * @param bfcsDir      Pointwise body-fitted coordinate system.
     * @param isPerpendicularCondition Flag for using perpendicular check or not
     * @param distTol      Distance tolerence. Used to find the boundary elements
     *                     for local coordinate iteration.
     * @param iterTol      Iteration tolerence. Used to check iteration convergence.
     * @param dirTol       Direction tolerencce. Used to check if the inner product
     *                     of two unit vectors is cloes enough to 1.0.                    
     * @param geoTol       Geometry tolerence. Used as the relative tolerence for
     *                     local coord and distance absolute tolerence.
     */ 
    void GenPntwiseBodyFittedCoordSys(
        const int targetBndId,
        const Array<OneD, NekDouble> assistVec,
        Array<OneD, NekDouble> & distance,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > & bfcsDir,
        const bool isCheckAngle,
        const NekDouble distTol = 1.0e-12,
        const NekDouble iterTol = 1.0e-12,
        const NekDouble dirTol  = 1.0e-4,
        const NekDouble geoTol  = 1.0e-12);

protected:

private:

    /**
     * @brief Compute the local coordinate for the nearest point on the given
     *        2D boundary element to the input point.
     * @param inGloCoord  Global coordinate for the input point.
     * @param bndGeom     Geometry of the boundary element to search from.
     * @param pts         Global coordinate of the quadrature points in the boundary element.
     * @param locCoord    Local coordinate of the result.
     * @param gloCoord    Global coordinate of the result.
     * @param dist        Distance from the input point to the boundary element.
     * @param iterTol     Iteration tolerence.
     * @param iterMax     Max iteration steps.
     */ 
    NekDouble PntToBndElmtPntDistance(
        const Array<OneD, Array<OneD, NekDouble> > & pts,
        const int pId,
        const Array<OneD, Array<OneD, NekDouble> > & bndPts);

    
    /**
     * @brief Compute the local coordinate for the nearest point on the given
     *        2D boundary element to the input point.
     * @param inGloCoord  Global coordinate for the input point.
     * @param bndGeom     Geometry of the boundary element to search from.
     * @param pts         Global coordinate of the quadrature points in the boundary element.
     * @param locCoord    Local coordinate of the result.
     * @param gloCoord    Global coordinate of the result.
     * @param dist        Distance from the input point to the boundary element.
     * @param iterTol     Iteration tolerence.
     * @param iterMax     Max iteration steps.
     */ 
    bool LocCoordForNearestPntOnBndElmt_2D(
        const Array<OneD, const NekDouble > & inGloCoord,
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array<OneD, Array<OneD, NekDouble> > & pts,
        Array<OneD, NekDouble > & locCoord,
        Array<OneD, NekDouble > & gloCoord,
        NekDouble & dist,
        const NekDouble iterTol = 1.0e-12,
        const int iterMax = 51);


    /**
     * @brief Compute the local coordinate for the nearest point on the given
     *        boundary element to the input point. The locCoord is the position to
     *        set up the body-fitted coordinate. This function works as a driver.
     * @param inGloCoord  Global coordinate for the input point.
     * @param bndGeom     Geometry of the boundary element to search from.
     * @param locCoord    Local coordinate of the result.
     * @param gloCoord    Global coordinate of the result.
     * @param dist        Distance from the input point to the boundary element.
     * @param iterTol     Iteration tolerence.
     * @param iterMax     Max iteration steps.
     */ 
    bool LocCoordForNearestPntOnBndElmt(
        const Array<OneD, const NekDouble > & inGloCoord,
        SpatialDomains::GeometrySharedPtr bndGeom,
        Array<OneD, NekDouble > & locCoord,
        Array<OneD, NekDouble > & gloCoord,
        NekDouble & dist,
        const NekDouble iterTol = 1.0e-12,
        const int iterMax = 51);


    /**
     * @brief Compute the normalized cross product for two 2D vectors.
     *        vec3 = vec1 x vec2
     */ 
    void ScaledCrosssProduct(
        const Array<OneD, NekDouble > & vec1,
        const Array<OneD, NekDouble > & vec2,
        Array<OneD, NekDouble > & vec3);


    /**
     * @brief Get velocity and convert to Cartesian system, if it is still in
     *        transformed system. It is copied and modified from from 
     *        ProcessGrad.cpp
     */ 
    void GetVelAndConvertToCartSys(
        Array<OneD, Array<OneD, NekDouble> > & vel);

   
};
}
}

#endif
