///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingReferenceFrame.cpp
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
// Description: Solving the absolute flow in a moving body frame,
// by adding (U0 + Omega X (x - x0)) . grad u - Omega X u
// as the body force.
// U0 is the translational velocity of the body frame.
// Omega is the angular velocity.
// x0 is the rotation pivot in the body frame.
// All vectors use the basis of the body frame.
// Translational motion is allowed for all dimensions.
// Rotation is not allowed for 1D, 2DH1D, 3DH2D.
// Rotation in z direction is allowed for 2D and 3DH1D.
// Rotation in 3 directions are allowed for 3D.
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Forcing/ForcingMovingReferenceFrame.h>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <SolverUtils/Filters/FilterInterfaces.hpp>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

std::string ForcingMovingReferenceFrame::classNameBody = GetForcingFactory().
        RegisterCreatorFunction("MovingReferenceFrame",
                                ForcingMovingReferenceFrame::create,
                                "Moving Frame");
std::string ForcingMovingReferenceFrame::classNameField = GetForcingFactory().
        RegisterCreatorFunction("Field",
                                ForcingMovingReferenceFrame::create,
                                "Field Forcing");

/**
 * @brief
 * @param pSession
 * @param pEquation
 */
ForcingMovingReferenceFrame::ForcingMovingReferenceFrame(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation)
        : Forcing(pSession, pEquation)
{
}


/**
 * @brief Initialise the forcing module
 * @param pFields
 * @param pNumForcingFields
 * @param pForce
 */
void ForcingMovingReferenceFrame::v_InitObject(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const unsigned int &pNumForcingFields, const TiXmlElement *pForce)
{
    boost::ignore_unused(pNumForcingFields);
    m_session->MatchSolverInfo("Homogeneous", "1D", m_isH1d, false);
    m_session->MatchSolverInfo("Homogeneous", "2D", m_isH2d, false);
    bool singleMode, halfMode;
    m_session->MatchSolverInfo("ModeType", "SingleMode", singleMode, false);
    m_session->MatchSolverInfo("ModeType", "HalfMode", halfMode, false);
    if (singleMode || halfMode)
    {
        m_isH1d = false;
    }

    int npoints   = pFields[0]->GetNpoints();
    int expdim    = m_isH2d ? 1 : pFields[0]->GetGraph()->GetMeshDimension();
    m_spacedim    = expdim + (m_isH1d ? 1 : 0) + (m_isH2d ? 2 : 0);
    m_NumVariable = m_spacedim;

    m_hasPlane0 = true;
    if (m_isH1d)
    {
        m_hasPlane0 = pFields[0]->GetZIDs()[0] == 0;
    }

    // linear velocities
    const TiXmlElement *funcNameElmt;

    funcNameElmt = pForce->FirstChildElement("LinearVelocity");
    ASSERTL0(funcNameElmt, "Requires LinearVelocity tag specifying function "
                           "name which prescribes velocity of the moving "
                           "frame.");

    m_velFuncName = funcNameElmt->GetText();
    ASSERTL0(m_session->DefinesFunction(m_velFuncName),
             "Function '" + m_velFuncName + "' is not defined in the session.");

    // angular velocities
    // Initiate the rotation to false, will be updated later if rotation defined
    m_hasRotation = false;
    funcNameElmt  = pForce->FirstChildElement("AngularVelocity");
    if (funcNameElmt)
    {
        m_omegaFuncName = funcNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(m_omegaFuncName),
                 "Function '" + m_omegaFuncName +
                     "' is not defined in the session.");
        m_hasRotation = true;
    }

    // initiate the linear velocity values in both local and inertial frames
    for (int i = 0; i < 3; ++i)
    {
        m_velXYZ[i] = 0.;
        m_velxyz[i] = 0.;
    }

    // initiate the angular velocity values in both local and inertial frames
    for (int i = 0; i < 3; ++i)
    {
        m_omegaXYZ[i] = 0;
        m_omegaxyz[i] = 0;
    }

    // initiate the available frame velocities all to false
    // this will be updated after reading the frame velocities and checking
    // for the rotation
    for (int i = 0; i < 3; ++i)
    {
        m_hasVel[i]   = false;
        m_hasOmega[i] = false;
    }

    // initiate the rotation angle to zero
    // m_theta = {theta_X, theta_Y, theta_Z}
    m_theta = Array<OneD, NekDouble>(3, 0.0);

    // initialize theta from the metadata of restart file or zero
    CheckForRestartTheta(pFields, m_theta);

    // Initialize theta with the data from NS class
    // This ensure correct moving coordinate orientation with respect to the
    // stationary inertial frame when we restart the simulation
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer to the equation system is expired");
    auto FluidEq = std::dynamic_pointer_cast<FluidInterface>(equ);
    // Set the angles
    FluidEq->SetMovingFrameAngles(m_theta);

    // initiate the rotation matrix to zero,
    // Note that the rotation matrix is assumed for the rotation around z axis
    // TODO: Generalize this for the 3D case with possiblity of rotation around
    // each of axis. Probabley we only can support rotation around one axis. In
    // that case the generalization means the user can provide Omega for one of
    // x, y or z axis. Not sure how complicated to consider a general axis of
    // rotation
    //
    // Note that these rotation matrices should be extrinsic rotations
    m_ProjMatX = bn::ublas::zero_matrix<NekDouble>(3, 3);
    m_ProjMatY = bn::ublas::zero_matrix<NekDouble>(3, 3);
    m_ProjMatZ = bn::ublas::zero_matrix<NekDouble>(3, 3);
    // populate the rotation matrix R(z)
    {
        NekDouble sn, cs;
        sn = sin(m_theta[2]);
        cs = cos(m_theta[2]);

        m_ProjMatZ(0, 0) = cs;
        m_ProjMatZ(0, 1) = sn;
        m_ProjMatZ(1, 0) = -1.*sn;
        m_ProjMatZ(1, 1) = cs;
        m_ProjMatZ(2, 2) = 1.0;
    }

    // frame linear velocity
    for (int i = 0; i < m_spacedim; ++i)
    {
        std::string s_FieldStr = m_session->GetVariable(i);

        if (m_session->DefinesFunction(m_velFuncName, s_FieldStr))
        {
            m_velFunction[i] =
                m_session->GetFunction(m_velFuncName, s_FieldStr);
            m_hasVel[i] = true;
        }
    }
    if (expdim == 1)
    {
        return;
    }

    // initialize the poivot coordinate with a default at origin
    m_pivotPoint = Array<OneD, NekDouble>(3, 0.0);
    if (m_hasRotation)
    {
        // frame angular velocity
        std::vector<std::string> angularVar = {"Omega_x", "Omega_y", "Omega_z"};
        for (int i = 0; i < 3; ++i)
        {
            std::string s_FieldStr = angularVar[i];
            if (m_session->DefinesFunction(m_omegaFuncName, s_FieldStr))
            {
                // m_hasRotation = true;
                m_omegaFunction[i] =
                    m_session->GetFunction(m_omegaFuncName, s_FieldStr);
                m_hasOmega[i] = true;
                // m_hasVel[m_spacedim+i]=true;
            }
        }

        // for now only support Omega_z
        // TODO: add the support for general rotation
        for (int i = 0; i < 2; ++i)
        {
            ASSERTL0( !m_hasOmega[i] , "Currently only Omega_z is supported");
        }

        // Reading the pivot point coordinate for rotation
        const TiXmlElement *pivotElmt = pForce->FirstChildElement("PivotPoint");
        // if not defined, zero would be the default
        if (pivotElmt)
        {
            std::stringstream pivotPointStream;
            std::string pivotPointStr = pivotElmt->GetText();
            pivotPointStream.str(pivotPointStr);

            for (int j = 0; j < m_spacedim; ++j)
            {
                pivotPointStream >> pivotPointStr;
                if (!pivotPointStr.empty())
                {
                    LibUtilities::Equation equ(m_session->GetInterpreter(),
                                               pivotPointStr);
                    m_pivotPoint[j] = equ.Evaluate();
                }
            }
        }

        m_coords = Array<OneD, Array<OneD, NekDouble>>(3);
        for (int j = 0; j < m_spacedim; ++j)
        {
            m_coords[j] = Array<OneD, NekDouble>(npoints);
        }
        pFields[0]->GetCoords(m_coords[0], m_coords[1], m_coords[2]);
        // move the origin to the pivot point
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Sadd(npoints, -m_pivotPoint[i], m_coords[i], 1, m_coords[i], 1);
        }

        // account for the effect of rotation
        // Omega_X results in having v and w even if not defined by user
        // Omega_Y results in having u and w even if not defined by user
        // Omega_Z results in having u and v even if not defined by user
        //
        for (int i = 0; i < 3; ++i)
        {
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;

            if (m_hasOmega[i])
            {
                m_hasVel[j] = true;
                m_hasVel[k] = true;
            }
        }
    }
}

/**
 * @brief Updates the forcing array with the current required forcing.
 * @param pFields
 * @param time
 */
 void ForcingMovingReferenceFrame::Update(const NekDouble &time)
{
    // compute the velociites whos functions are provided in inertial frame
    for (auto it: m_velFunction)
    {
        m_velXYZ[it.first]=it.second->Evaluate(0.,0.,0.,time);
    }
    for (auto it: m_omegaFunction)
    {
        m_omegaXYZ[it.first]=it.second->Evaluate(0.,0.,0.,time);
    }
    // include the effect of rotation
    if(m_hasRotation)
    {
        UpdateTheta(time);

        // transform the translation velocities to the moving frame
        bn::ublas::vector<NekDouble> v0 = bn::ublas::zero_vector<NekDouble>(3);
        bn::ublas::vector<NekDouble> v1 = bn::ublas::zero_vector<NekDouble>(3);
        for (int i = 0; i < m_spacedim; ++i)
        {
            v0(i) = m_velXYZ[i];
        }
        v1 = bn::ublas::prec_prod(m_ProjMatZ, v0);
        for (int i = 0; i < 3; ++i)
        {
            m_velxyz[i] = v1(i);
        }

        // transform the angular velocities to moving frame
        v0 = bn::ublas::zero_vector<NekDouble>(3);
        v1 = bn::ublas::zero_vector<NekDouble>(3);
        for (auto it: m_omegaXYZ)
        {
            v0(it.first)=it.second;
        }
        v1 = bn::ublas::prec_prod(m_ProjMatZ, v0);
        for (int i = 0; i < 3; ++i)
        {
            m_omegaxyz[i] = v1(i);
        }
    }
    else
    {
        // for translation only,
        for(int i=0; i<m_spacedim; ++i)
        {
            m_velxyz[i] = m_velXYZ[i];
        }

    }

    // set the velocities and rotation matrix for enforcement of boundary
    // conditions in the Incompressible Naveri-Stokes
    Array<OneD, NekDouble> vel(6,0.0);
    // save the linear and angular velocities to be used in enforcing bc
    for (int i = 0; i < 3; ++i)
    {
        vel[i]     = m_velxyz[i];
        vel[i + 3] = m_omegaxyz[i];
    }

    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer to the equation system is expired");
    auto FluidEq = std::dynamic_pointer_cast<FluidInterface>(equ);
    // update the frame velocities
    FluidEq->SetMovingFrameVelocities(vel);
    // update the projection matrix
    FluidEq->SetMovingFrameProjectionMat(m_ProjMatZ);
    // update the frame angle (with respect to the inertial frame)
    // this angle is used to update the meta data,
    // on the other hand, for boundary conditions the projection matrix is used
    FluidEq->SetMovingFrameAngles(m_theta);

}

void ForcingMovingReferenceFrame::UpdateTheta(const NekDouble& time)
{

    NekDouble dt = m_session->GetParameter("TimeStep");
    NekDouble t1 = time + dt;
    std::map<int,NekDouble> Omega;

    // Calculate angular velocities at dt
    // TODO: Generalize to the case with a general 3D rotation support
    for (int i = 0; i < 3; ++i)
    {
        if (m_hasOmega[i])
        {
            Omega[i] = 0.5 * (m_omegaXYZ[i] +
                              m_omegaFunction[i]->Evaluate(0., 0., 0., t1));
        }
    }

    // Using the Omega_z
    m_theta[2] += (Omega[2]*dt);

    // update the rotation matrix
    {
        NekDouble sn, cs;
        sn = sin(m_theta[2]);
        cs = cos(m_theta[2]);

        m_ProjMatZ(0,0)=cs;
        m_ProjMatZ(0,1)=sn;
        m_ProjMatZ(1,0)=-sn;
        m_ProjMatZ(1,1)=cs;
        m_ProjMatZ(2,2)=1.0;
    }
}

/**
 * @brief Adds the body force, -Omega X u.
 * @param fields
 * @param inarray
 * @param outarray
 * @param time
 */
 void ForcingMovingReferenceFrame::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> > &inarray,
              Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble &time)
{
    boost::ignore_unused(time);
    // If there is no rotation, body force is zero,
    // nothing needs to be done here.
    if (!m_hasRotation)
    {
        return;
    }
    //frame velocities are already updated in pre_Apply
    int  npoints = fields[0]->GetNpoints();
    boost::ignore_unused(npoints);
    addRotation(npoints, outarray, -1., inarray, outarray);
}

/**
 * @brief outarray = inarray0 + angVelScale Omega x inarray1
 */
void ForcingMovingReferenceFrame::addRotation(
    int nPnts, // number of points
    const Array<OneD, Array<OneD, NekDouble>> &inarray0, NekDouble angVelScale,
    const Array<OneD, Array<OneD, NekDouble>> &inarray1,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    ASSERTL0(&inarray1 != &outarray, "inarray1 and outarray "
                                     "should not be the same.");

    // TODO: In case of having support for all three components of Omega,
    // they should be transformed into the rotating frame first!

    // In case that the inarray0 and outarry are different, to avoid using
    // un-initialized array, copy the array first
    if (&inarray0 != &outarray)
    {
        ASSERTL0(inarray0.size() == outarray.size(),
                 "inarray0 and outarray must have same dimentions");
        for (int i = 0; i < inarray0.size(); ++i)
        {
            Vmath::Vcopy(nPnts, inarray0[i], 1, outarray[i], 1);
        }
    }

    if (m_spacedim >= 2 && m_hasOmega[2])
    {
        NekDouble cp = m_omegaxyz[2] * angVelScale;
        NekDouble cm = -1. * cp;

        Vmath::Svtvp(nPnts, cm, inarray1[1], 1, outarray[0], 1, outarray[0], 1);
        Vmath::Svtvp(nPnts, cp, inarray1[0], 1, outarray[1], 1, outarray[1], 1);
    }

    if (m_spacedim == 3 && m_hasOmega[0])
    {
        NekDouble cp = m_omegaxyz[0] * angVelScale;
        NekDouble cm = -1. * cp;

        Vmath::Svtvp(nPnts, cp, inarray1[1], 1, outarray[2], 1, outarray[2], 1);
        Vmath::Svtvp(nPnts, cm, inarray1[2], 1, outarray[1], 1, outarray[1], 1);
    }

    if (m_spacedim == 3 && m_hasOmega[1])
    {
        NekDouble cp = m_omegaxyz[1] * angVelScale;
        NekDouble cm = -1. * cp;

        Vmath::Svtvp(nPnts, cp, inarray1[2], 1, outarray[0], 1, outarray[0], 1);
        Vmath::Svtvp(nPnts, cm, inarray1[0], 1, outarray[2], 1, outarray[2], 1);
    }
}

void ForcingMovingReferenceFrame::v_PreApply(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble &time)
{
    Update(time);
    int npoints = fields[0]->GetNpoints();
    if (m_isH2d && fields[0]->GetWaveSpace())
    {
        for (int i=0; i<m_spacedim; ++i)
        {
            if (m_hasVel[i])
            {
                Array<OneD, NekDouble> tmpphys(npoints, -m_velxyz[i]);
                Array<OneD, NekDouble> tmpcoef(npoints);
                fields[0]->HomogeneousFwdTrans(tmpphys, tmpcoef);
                Vmath::Vadd(npoints, tmpcoef, 1, inarray[i], 1,
                outarray[i], 1);
            }
            else if (&inarray != &outarray)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }
        }
    }
    else
    {
        int npoints0 = npoints;
        if (m_isH1d && fields[0]->GetWaveSpace())
        {
            npoints0 = m_hasPlane0 ? fields[0]->GetPlane(0)->GetNpoints()
            : 0;
        }
        for (int i=0; i<m_spacedim; ++i)
        {
            if (m_hasVel[i])
            {
                Vmath::Sadd(npoints0, -m_velxyz[i], inarray[i], 1,
                outarray[i], 1); if (&inarray != &outarray && npoints !=
                npoints0)
                {
                    Array<OneD, NekDouble> tmp = outarray[i]+npoints0;
                    Vmath::Vcopy(npoints - npoints0, inarray[i] +
                    npoints0, 1, tmp, 1);
                }
            }
            else if (&inarray != &outarray)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }
        }
        if (m_hasRotation)
        {
            addRotation(npoints0, outarray, -1., m_coords, outarray);
        }
    }
}

void ForcingMovingReferenceFrame::CheckForRestartTheta(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    Array<OneD, NekDouble> &theta)
{
    std::map<std::string, std::string> fieldMetaDataMap;

    // initialize theta to zero
    for(auto &t: theta)
    {
        t=0.0;
    }

    if (m_session->DefinesFunction("InitialConditions"))
    {
        for (int i = 0; i < pFields.size(); ++i)
        {
            LibUtilities::FunctionType vType;

            vType = m_session->GetFunctionType("InitialConditions",
                                               m_session->GetVariable(i));

            if (vType == LibUtilities::eFunctionTypeFile)
            {
                std::string filename = m_session->GetFunctionFilename(
                    "InitialConditions", m_session->GetVariable(i));

                fs::path pfilename(filename);

                // redefine path for parallel file which is in directory
                if (fs::is_directory(pfilename))
                {
                    fs::path metafile("Info.xml");
                    fs::path fullpath = pfilename / metafile;
                    filename          = LibUtilities::PortablePath(fullpath);
                }
                LibUtilities::FieldIOSharedPtr fld =
                    LibUtilities::FieldIO::CreateForFile(m_session, filename);
                fld->ImportFieldMetaData(filename, fieldMetaDataMap);

                // check to see if theta is defined
                if (fieldMetaDataMap != LibUtilities::NullFieldMetaDataMap)
                {
                    std::vector<std::string> vSuffix = {"_x", "_y", "_z"};

                    for (int j = 0; j < 3; ++j)
                    {
                        std::string sTheta = "Theta" + vSuffix[j];
                        auto iter = fieldMetaDataMap.find(sTheta);
                        if (iter != fieldMetaDataMap.end())
                        {
                            theta[j] =
                                boost::lexical_cast<NekDouble>(iter->second);
                        }
                    }

                    break;
                }

                break;
            }
        }
    }
}


}
}
