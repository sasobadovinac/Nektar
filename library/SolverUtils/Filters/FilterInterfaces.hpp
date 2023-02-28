///////////////////////////////////////////////////////////////////////////////
//
// File: FilterInterfaces.hpp
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
// Description: Interface class for solvers that support fluid physics
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERINTERFACES_HPP
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERINTERFACES_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <SolverUtils/SolverUtilsDeclspec.h>

namespace Nektar
{
namespace SolverUtils
{

class FluidInterface
{
public:
    virtual ~FluidInterface() = default;

    /// Extract array with velocity from physfield
    SOLVER_UTILS_EXPORT void GetVelocity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &velocity);

    SOLVER_UTILS_EXPORT bool HasConstantDensity();

    /// Extract array with density from physfield
    SOLVER_UTILS_EXPORT void GetDensity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &density);

    /// Extract array with pressure from physfield
    SOLVER_UTILS_EXPORT void GetPressure(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &pressure);

    // gave access and set to the moving frame velocity
    // for Moving reference frame formulation
    SOLVER_UTILS_EXPORT void SetMovingFrameVelocities(
        const Array<OneD, NekDouble> &vFrameVels);

    SOLVER_UTILS_EXPORT void GetMovingFrameVelocities(
        Array<OneD, NekDouble> &vFrameVels);

    // gave access and set to the projection matrix that transfers
    // between stationary inertial frame and moving reference frame
    SOLVER_UTILS_EXPORT void SetMovingFrameProjectionMat(
        const boost::numeric::ublas::matrix<NekDouble> &vProjMat);

    SOLVER_UTILS_EXPORT void GetMovingFrameProjectionMat(
        boost::numeric::ublas::matrix<NekDouble> &vProjMat);

    // gave access and set the angle between moving frame and stationary
    // inertial frame
    SOLVER_UTILS_EXPORT void SetMovingFrameAngles(
        const Array<OneD, NekDouble> &vFrameTheta);

    SOLVER_UTILS_EXPORT void GetMovingFrameAngles(
        Array<OneD, NekDouble> &vFrameTheta);

protected:
    SOLVER_UTILS_EXPORT virtual void v_GetVelocity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &velocity)      = 0;
    SOLVER_UTILS_EXPORT virtual bool v_HasConstantDensity() = 0;
    SOLVER_UTILS_EXPORT virtual void v_GetDensity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &density) = 0;
    SOLVER_UTILS_EXPORT virtual void v_GetPressure(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &pressure) = 0;
    SOLVER_UTILS_EXPORT virtual void v_SetMovingFrameVelocities(
        const Array<OneD, NekDouble> &vFrameVels)
    {
        boost::ignore_unused(vFrameVels);
    }
    SOLVER_UTILS_EXPORT virtual void v_GetMovingFrameVelocities(
        Array<OneD, NekDouble> &vFrameVels)
    {
        boost::ignore_unused(vFrameVels);
    }
    SOLVER_UTILS_EXPORT virtual void v_SetMovingFrameProjectionMat(
        const boost::numeric::ublas::matrix<NekDouble> &vProjMat)
    {
        boost::ignore_unused(vProjMat);
    }
    SOLVER_UTILS_EXPORT virtual void v_GetMovingFrameProjectionMat(
        boost::numeric::ublas::matrix<NekDouble> &vProjMat)
    {
        boost::ignore_unused(vProjMat);
    }
    SOLVER_UTILS_EXPORT virtual void v_SetMovingFrameAngles(
        const Array<OneD, NekDouble> &vFrameTheta)
    {
        boost::ignore_unused(vFrameTheta);
    }
    SOLVER_UTILS_EXPORT virtual void v_GetMovingFrameAngles(
        Array<OneD, NekDouble> &vFrameTheta)
    {
        boost::ignore_unused(vFrameTheta);
    }
};

/**
 *
 */
inline void FluidInterface::GetVelocity(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &velocity)
{
    v_GetVelocity(physfield, velocity);
}

/**
 *
 */
inline bool FluidInterface::HasConstantDensity()
{
    return v_HasConstantDensity();
}

/**
 *
 */
inline void FluidInterface::GetDensity(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &density)
{
    v_GetDensity(physfield, density);
}

/**
 *
 */
inline void FluidInterface::GetPressure(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &pressure)
{
    v_GetPressure(physfield, pressure);
}

/**
 *
 */
inline void FluidInterface::SetMovingFrameVelocities(
    const Array<OneD, NekDouble> &vFrameVels)
{
    v_SetMovingFrameVelocities(vFrameVels);
}

/**
 *
 */
inline void FluidInterface::GetMovingFrameVelocities(
    Array<OneD, NekDouble> &vFrameVels)
{
    v_GetMovingFrameVelocities(vFrameVels);
}

/**
 *
 */
inline void FluidInterface::SetMovingFrameProjectionMat(
    const boost::numeric::ublas::matrix<NekDouble> &vProjMat)
{
    v_SetMovingFrameProjectionMat(vProjMat);
}

/**
 *
 */
inline void FluidInterface::GetMovingFrameProjectionMat(
    boost::numeric::ublas::matrix<NekDouble> &vProjMat)
{
    v_GetMovingFrameProjectionMat(vProjMat);
}

/**
 *
 */
inline void FluidInterface::SetMovingFrameAngles(
    const Array<OneD, NekDouble> &vFrameTheta)
{
    v_SetMovingFrameAngles(vFrameTheta);
}

inline void FluidInterface::GetMovingFrameAngles(
    Array<OneD, NekDouble> &vFrameTheta)
{
    v_GetMovingFrameAngles(vFrameTheta);
}

} // namespace SolverUtils
} // namespace Nektar

#endif
