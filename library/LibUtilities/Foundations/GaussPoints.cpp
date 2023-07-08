///////////////////////////////////////////////////////////////////////////////
//
// File: GaussPoints.cpp
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
// Description: GaussPoints Definitions
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/GaussPoints.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Polylib/Polylib.h>

namespace Nektar
{
namespace LibUtilities
{
bool GaussPoints::initPointsManager[] = {
    PointsManager().RegisterCreator(PointsKey(0, eGaussGaussLegendre),
                                    GaussPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMLegendre),
                                    GaussPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eGaussRadauPLegendre),
                                    GaussPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eGaussLobattoLegendre),
                                    GaussPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eGaussGaussChebyshev),
                                    GaussPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMChebyshev),
                                    GaussPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eGaussRadauPChebyshev),
                                    GaussPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eGaussLobattoChebyshev),
                                    GaussPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta1),
                                    GaussPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta2),
                                    GaussPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMAlpha1Beta0),
                                    GaussPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMAlpha2Beta0),
                                    GaussPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eGaussKronrodLegendre),
                                    GaussPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eGaussRadauKronrodMLegendre),
                                    GaussPoints::Create),
    PointsManager().RegisterCreator(
        PointsKey(0, eGaussRadauKronrodMAlpha1Beta0), GaussPoints::Create),
    PointsManager().RegisterCreator(PointsKey(0, eGaussLobattoKronrodLegendre),
                                    GaussPoints::Create)};

void GaussPoints::v_CalculatePoints()
{
    // Allocate the storage for points and weights
    PointsBaseType::v_CalculatePoints();
    PointsBaseType::v_CalculateWeights();

    size_t numpoints = m_pointsKey.GetNumPoints();

    switch (m_pointsKey.GetPointsType())
    {
        case eGaussGaussLegendre:
            Polylib::zwgj(m_points[0].data(), m_weights.data(), numpoints, 0.0,
                          0.0);
            break;

        case eGaussRadauMLegendre:
            Polylib::zwgrjm(m_points[0].data(), m_weights.data(), numpoints,
                            0.0, 0.0);
            break;

        case eGaussRadauPLegendre:
            Polylib::zwgrjp(m_points[0].data(), m_weights.data(), numpoints,
                            0.0, 0.0);
            break;

        case eGaussLobattoLegendre:
            Polylib::zwglj(m_points[0].data(), m_weights.data(), numpoints, 0.0,
                           0.0);
            break;

        case eGaussGaussChebyshev:
            Polylib::zwgj(m_points[0].data(), m_weights.data(), numpoints, -0.5,
                          -0.5);
            break;

        case eGaussRadauMChebyshev:
            Polylib::zwgrjm(m_points[0].data(), m_weights.data(), numpoints,
                            -0.5, -0.5);
            break;

        case eGaussRadauPChebyshev:
            Polylib::zwgrjp(m_points[0].data(), m_weights.data(), numpoints,
                            -0.5, -0.5);
            break;

        case eGaussLobattoChebyshev:
            Polylib::zwglj(m_points[0].data(), m_weights.data(), numpoints,
                           -0.5, -0.5);
            break;

        case eGaussRadauMAlpha0Beta1:
            Polylib::zwgrjm(m_points[0].data(), m_weights.data(), numpoints,
                            0.0, 1.0);
            break;

        case eGaussRadauMAlpha0Beta2:
            Polylib::zwgrjm(m_points[0].data(), m_weights.data(), numpoints,
                            0.0, 2.0);
            break;

        case eGaussRadauMAlpha1Beta0:
            Polylib::zwgrjm(m_points[0].data(), m_weights.data(), numpoints,
                            1.0, 0.0);
            break;

        case eGaussRadauMAlpha2Beta0:
            Polylib::zwgrjm(m_points[0].data(), m_weights.data(), numpoints,
                            2.0, 0.0);
            break;

        case eGaussKronrodLegendre:
            Polylib::zwgk(m_points[0].data(), m_weights.data(), numpoints, 0.0,
                          0.0);
            break;

        case eGaussRadauKronrodMLegendre:
            Polylib::zwrk(m_points[0].data(), m_weights.data(), numpoints, 0.0,
                          0.0);
            break;

        case eGaussRadauKronrodMAlpha1Beta0:
            Polylib::zwrk(m_points[0].data(), m_weights.data(), numpoints, 1.0,
                          0.0);
            break;

        case eGaussLobattoKronrodLegendre:
            Polylib::zwlk(m_points[0].data(), m_weights.data(), numpoints, 0.0,
                          0.0);
            break;

        default:
            NEKERROR(ErrorUtil::efatal,
                     "Unknown Gauss quadrature point distribution requested");
    }
}

void GaussPoints::v_CalculateWeights()
{
    // For Gauss Quadrature, this is done as part of the points computation
}

void GaussPoints::v_CalculateDerivMatrix()
{
    // Allocate the derivative matrix
    PointsBaseType::v_CalculateDerivMatrix();

    size_t numpoints = m_pointsKey.GetNumPoints();
    size_t totpoints = m_pointsKey.GetTotNumPoints();
    Array<OneD, NekDouble> dmtempSharedArray =
        Array<OneD, NekDouble>(totpoints * totpoints);
    NekDouble *dmtemp = dmtempSharedArray.data();

    switch (m_pointsKey.GetPointsType())
    {
        case eGaussGaussLegendre:
            Polylib::Dgj(dmtemp, m_points[0].data(), numpoints, 0.0, 0.0);
            break;

        case eGaussRadauMLegendre:
            Polylib::Dgrjm(dmtemp, m_points[0].data(), numpoints, 0.0, 0.0);
            break;

        case eGaussRadauPLegendre:
            Polylib::Dgrjp(dmtemp, m_points[0].data(), numpoints, 0.0, 0.0);
            break;

        case eGaussLobattoLegendre:
            Polylib::Dglj(dmtemp, m_points[0].data(), numpoints, 0.0, 0.0);
            break;

        case eGaussGaussChebyshev:
            Polylib::Dgj(dmtemp, m_points[0].data(), numpoints, -0.5, -0.5);
            break;

        case eGaussRadauMChebyshev:
            Polylib::Dgrjm(dmtemp, m_points[0].data(), numpoints, -0.5, -0.5);
            break;

        case eGaussRadauPChebyshev:
            Polylib::Dgrjp(dmtemp, m_points[0].data(), numpoints, -0.5, -0.5);
            break;

        case eGaussLobattoChebyshev:
            Polylib::Dglj(dmtemp, m_points[0].data(), numpoints, -0.5, -0.5);
            break;

        case eGaussRadauMAlpha0Beta1:
            Polylib::Dgrjm(dmtemp, m_points[0].data(), numpoints, 0.0, 1.0);
            break;

        case eGaussRadauMAlpha0Beta2:
            Polylib::Dgrjm(dmtemp, m_points[0].data(), numpoints, 0.0, 2.0);
            break;

        case eGaussRadauMAlpha1Beta0:
            Polylib::Dgrjm(dmtemp, m_points[0].data(), numpoints, 1.0, 0.0);
            break;

        case eGaussRadauMAlpha2Beta0:
            Polylib::Dgrjm(dmtemp, m_points[0].data(), numpoints, 2.0, 0.0);
            break;

        case eGaussKronrodLegendre:
        case eGaussRadauKronrodMLegendre:
        case eGaussRadauKronrodMAlpha1Beta0:
        case eGaussLobattoKronrodLegendre:
        {
            for (size_t i = 0; i < m_pointsKey.GetNumPoints(); ++i)
            {
                for (size_t j = 0; j < m_pointsKey.GetNumPoints(); ++j)
                {
                    (*m_derivmatrix[0])(i, j) = Polylib::laginterpderiv(
                        m_points[0][i], j, &m_points[0][0],
                        m_pointsKey.GetNumPoints());
                }
            }
            return;
        }
        break;

        default:
            NEKERROR(ErrorUtil::efatal,
                     "Unknown Gauss quadrature point distribution requested");
    }

    std::copy(dmtemp, dmtemp + totpoints * totpoints,
              m_derivmatrix[0]->begin());
}

void GaussPoints::CalculateInterpMatrix(
    size_t npts, const Array<OneD, const NekDouble> &xpoints,
    Array<OneD, NekDouble> &interp)
{
    switch (m_pointsKey.GetPointsType())
    {
        case eGaussGaussLegendre:
            Polylib::Imgj(interp.data(), m_points[0].data(), xpoints.data(),
                          GetNumPoints(), npts, 0.0, 0.0);
            break;

        case eGaussRadauMLegendre:
            Polylib::Imgrjm(interp.data(), m_points[0].data(), xpoints.data(),
                            GetNumPoints(), npts, 0.0, 0.0);
            break;

        case eGaussRadauPLegendre:
            Polylib::Imgrjp(interp.data(), m_points[0].data(), xpoints.data(),
                            GetNumPoints(), npts, 0.0, 0.0);
            break;

        case eGaussLobattoLegendre:
            Polylib::Imglj(interp.data(), m_points[0].data(), xpoints.data(),
                           GetNumPoints(), npts, 0.0, 0.0);
            break;

        case eGaussGaussChebyshev:
            Polylib::Imgj(interp.data(), m_points[0].data(), xpoints.data(),
                          GetNumPoints(), npts, -0.5, -0.5);
            break;

        case eGaussRadauMChebyshev:
            Polylib::Imgrjm(interp.data(), m_points[0].data(), xpoints.data(),
                            GetNumPoints(), npts, -0.5, -0.5);
            break;

        case eGaussRadauPChebyshev:
            Polylib::Imgrjp(interp.data(), m_points[0].data(), xpoints.data(),
                            GetNumPoints(), npts, -0.5, -0.5);
            break;

        case eGaussLobattoChebyshev:
            Polylib::Imglj(interp.data(), m_points[0].data(), xpoints.data(),
                           GetNumPoints(), npts, -0.5, -0.5);
            break;

        case eGaussRadauMAlpha0Beta1:
            Polylib::Imgrjm(interp.data(), m_points[0].data(), xpoints.data(),
                            GetNumPoints(), npts, 0.0, 1.0);
            break;

        case eGaussRadauMAlpha0Beta2:
            Polylib::Imgrjm(interp.data(), m_points[0].data(), xpoints.data(),
                            GetNumPoints(), npts, 0.0, 2.0);
            break;

        case eGaussRadauMAlpha1Beta0:
            Polylib::Imgrjm(interp.data(), m_points[0].data(), xpoints.data(),
                            GetNumPoints(), npts, 1.0, 0.0);
            break;

        case eGaussRadauMAlpha2Beta0:
            Polylib::Imgrjm(interp.data(), m_points[0].data(), xpoints.data(),
                            GetNumPoints(), npts, 2.0, 0.0);
            break;

        case eGaussKronrodLegendre:
        case eGaussRadauKronrodMLegendre:
        case eGaussRadauKronrodMAlpha1Beta0:
        case eGaussLobattoKronrodLegendre:
        {
            for (size_t i = 0; i < npts; ++i)
            {
                for (size_t j = 0; j < m_pointsKey.GetNumPoints(); ++j)
                {
                    interp[i + j * npts] =
                        Polylib::laginterp(xpoints[i], j, &m_points[0][0],
                                           m_pointsKey.GetNumPoints());
                }
            }
        }
        break;

        default:
            NEKERROR(ErrorUtil::efatal,
                     "Unknown Gauss quadrature point distribution requested");
    }
}

std::shared_ptr<Points<NekDouble>> GaussPoints::Create(const PointsKey &pkey)
{
    std::shared_ptr<Points<NekDouble>> returnval(
        MemoryManager<GaussPoints>::AllocateSharedPtr(pkey));

    returnval->Initialize();

    return returnval;
}

std::shared_ptr<NekMatrix<NekDouble>> GaussPoints::CreateMatrix(
    const PointsKey &pkey)
{
    size_t numpoints = pkey.GetNumPoints();
    Array<OneD, const NekDouble> xpoints;

    PointsManager()[pkey]->GetPoints(xpoints);

    // Delegate to function below
    return GetI(numpoints, xpoints);
}

const std::shared_ptr<NekMatrix<NekDouble>> GaussPoints::v_GetI(
    const PointsKey &pkey)
{
    ASSERTL0(pkey.GetPointsDim() == 1,
             "Gauss Points can only interp to other 1d point distributions");

    return m_InterpManager[pkey];
}

const std::shared_ptr<NekMatrix<NekDouble>> GaussPoints::v_GetI(
    const Array<OneD, const NekDouble> &x)
{
    size_t numpoints = 1;

    // Delegate to function below
    return GetI(numpoints, x);
}

const std::shared_ptr<NekMatrix<NekDouble>> GaussPoints::v_GetI(
    size_t numpoints, const Array<OneD, const NekDouble> &x)
{
    Array<OneD, NekDouble> interp(GetNumPoints() * numpoints);

    CalculateInterpMatrix(numpoints, x, interp);

    NekDouble *t = interp.data();
    size_t np    = GetNumPoints();
    std::shared_ptr<NekMatrix<NekDouble>> returnval(
        MemoryManager<NekMatrix<NekDouble>>::AllocateSharedPtr(numpoints, np,
                                                               t));

    return returnval;
}

const std::shared_ptr<NekMatrix<NekDouble>> GaussPoints::
    v_GetGalerkinProjection(const PointsKey &pkey)
{
    return m_GalerkinProjectionManager[pkey];
}

std::shared_ptr<NekMatrix<NekDouble>> GaussPoints::CreateGPMatrix(
    const PointsKey &pkey)
{
    std::shared_ptr<NekMatrix<NekDouble>> returnval =
        CalculateGalerkinProjectionMatrix(pkey);

    // Delegate to function below
    return returnval;
}

std::shared_ptr<NekMatrix<NekDouble>> GaussPoints::
    CalculateGalerkinProjectionMatrix(const PointsKey &pkey)
{
    size_t numpointsfrom = pkey.GetNumPoints();
    size_t numpointsto   = GetNumPoints();

    Array<OneD, const NekDouble> weightsfrom;

    weightsfrom = PointsManager()[pkey]->GetW();

    std::shared_ptr<NekMatrix<NekDouble>> Interp = GetI(pkey);

    Array<OneD, NekDouble> GalProj(numpointsfrom * numpointsto);

    // set up inner product matrix and multiply by inverse of
    // diagaonal mass matrix
    for (size_t i = 0; i < numpointsto; ++i)
    {
        Vmath::Vmul(numpointsfrom, Interp->GetPtr().get() + i * numpointsfrom,
                    1, &weightsfrom[0], 1, &GalProj[0] + i, numpointsto);
        Vmath::Smul(numpointsfrom, 1.0 / m_weights[i], &GalProj[0] + i,
                    numpointsto, &GalProj[0] + i, numpointsto);
    }

    NekDouble *t = GalProj.data();
    std::shared_ptr<NekMatrix<NekDouble>> returnval(
        MemoryManager<NekMatrix<NekDouble>>::AllocateSharedPtr(
            numpointsto, numpointsfrom, t));

    return returnval;
}

} // end of namespace LibUtilities
} // end of namespace Nektar
