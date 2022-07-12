///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingReferenceFrame.h
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
// Description: Allows for a moving frame of reference, through adding c * du/dx
// to the body force, where c is the frame velocity vector
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCINGMOVINGREFERENCEFRAME
#define NEKTAR_SOLVERUTILS_FORCINGMOVINGREFERENCEFRAME

#include <string>

#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cmath>

namespace Nektar
{
namespace SolverUtils
{

namespace bn = boost::numeric;

class ForcingMovingReferenceFrame : public Forcing
{

public:
    friend class MemoryManager<ForcingMovingReferenceFrame>;

    /// Creates an instance of this class
    SOLVER_UTILS_EXPORT static ForcingSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const unsigned int &pNumForcingFields, const TiXmlElement *pForce)
    {
        ForcingSharedPtr p =
            MemoryManager<ForcingMovingReferenceFrame>::AllocateSharedPtr(
                pSession, pEquation);
        p->InitObject(pFields, pNumForcingFields, pForce);
        return p;
    }

    /// Name of the class
    static std::string classNameBody;
    static std::string classNameField;

protected:
    SOLVER_UTILS_EXPORT virtual void v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const unsigned int &pNumForcingFields, const TiXmlElement *pForce);

    SOLVER_UTILS_EXPORT virtual void v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble &time);

    SOLVER_UTILS_EXPORT virtual void v_PreApply(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble &time);

private:
    // name of the function for linear and angular velocities in the session
    // file
    std::string m_funcName;
    std::string m_velFuncName;
    std::string m_omegaFuncName;

    // prescribed functions in the session file
    std::map<int, LibUtilities::EquationSharedPtr> m_frameFunction;
    std::map<int, LibUtilities::EquationSharedPtr> m_velFunction;
    std::map<int, LibUtilities::EquationSharedPtr> m_omegaFunction;

    // a boolean switch indicating for which direction the velocities are
    // available. The available velocites could be different from the
    // precscribed one because of the rotation which result in change of basis
    // vector of local frame to the inertial frame.
    std::map<int, bool> m_hasVel;
    std::map<int, bool> m_hasOmega;

    // frame linear velocities in inertial frame
    std::map<int, NekDouble> m_velXYZ;

    // frame linear velocities in local translating-rotating frame
    std::map<int, NekDouble> m_velxyz;

    // frame angular velocities in inertial frame
    std::map<int, NekDouble> m_omegaXYZ;

    // frame angular velocities in local translating-rotating frame
    std::map<int, NekDouble> m_omegaxyz;
    // coordinate vector
    Array<OneD, Array<OneD, NekDouble>> m_coords;

    // pivot point
    Array<OneD, NekDouble> m_pivotPoint;

    // rotation angel
    Array<OneD, NekDouble> m_theta;

    // Projection matrix for transformation of vectors between inertial and
    // moving reference frames
    bn::ublas::matrix<NekDouble> m_ProjMatX;
    bn::ublas::matrix<NekDouble> m_ProjMatY;
    bn::ublas::matrix<NekDouble> m_ProjMatZ;

    bool m_hasRotation;
    bool m_isH1d;
    bool m_hasPlane0;
    bool m_isH2d;
    NekInt m_spacedim;

    ForcingMovingReferenceFrame(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation);

    virtual ~ForcingMovingReferenceFrame(void){};

    void Update(const NekDouble &time);
    void UpdateTheta(const NekDouble &time);
    void CheckForRestartTheta(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        Array<OneD, NekDouble> &theta);

    void addRotation(int npoints,
                     const Array<OneD, Array<OneD, NekDouble>> &inarray0,
                     NekDouble angVelScale,
                     const Array<OneD, Array<OneD, NekDouble>> &inarray1,
                     Array<OneD, Array<OneD, NekDouble>> &outarray);
};

} // namespace SolverUtils
} // namespace Nektar

#endif
