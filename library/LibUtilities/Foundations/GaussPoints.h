///////////////////////////////////////////////////////////////////////////////
//
// File: GaussPoints.h
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
// Description: Header file of GaussPoints Distributions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef GAUSSPOINTS_H
#define GAUSSPOINTS_H

#include <LibUtilities/Foundations/Points.h>

namespace Nektar
{
namespace LibUtilities
{
class GaussPoints : public Points<NekDouble>
{
public:
    virtual ~GaussPoints()
    {
    }

    LIB_UTILITIES_EXPORT static std::shared_ptr<Points<NekDouble>> Create(
        const PointsKey &pkey);

    LIB_UTILITIES_EXPORT std::shared_ptr<NekMatrix<NekDouble>> CreateMatrix(
        const PointsKey &pkey);

    LIB_UTILITIES_EXPORT std::shared_ptr<NekMatrix<NekDouble>> CreateGPMatrix(
        const PointsKey &pkey);

    GaussPoints(const PointsKey &pkey) : PointsBaseType(pkey)
    {
        namespace pl = std::placeholders;
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussGaussLegendre),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussRadauMLegendre),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussRadauPLegendre),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussLobattoLegendre),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussGaussChebyshev),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussRadauMChebyshev),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussRadauPChebyshev),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussLobattoChebyshev),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussRadauMAlpha0Beta1),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussRadauMAlpha0Beta2),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussRadauMAlpha1Beta0),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussRadauMAlpha2Beta0),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussKronrodLegendre),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussRadauKronrodMLegendre),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussRadauKronrodMAlpha1Beta0),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eGaussLobattoKronrodLegendre),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eFourierEvenlySpaced),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, ePolyEvenlySpaced),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eBoundaryLayerPoints),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_InterpManager.RegisterCreator(
            PointsKey(0, eBoundaryLayerPointsRev),
            std::bind(&GaussPoints::CreateMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussGaussLegendre),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussRadauMLegendre),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussRadauPLegendre),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussLobattoLegendre),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussGaussChebyshev),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussRadauMChebyshev),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussRadauPChebyshev),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussLobattoChebyshev),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussRadauMAlpha0Beta1),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussRadauMAlpha0Beta2),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussRadauMAlpha1Beta0),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussRadauMAlpha2Beta0),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussKronrodLegendre),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussRadauKronrodMLegendre),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussRadauKronrodMAlpha1Beta0),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eGaussLobattoKronrodLegendre),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eFourierEvenlySpaced),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, ePolyEvenlySpaced),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eBoundaryLayerPoints),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
        m_GalerkinProjectionManager.RegisterCreator(
            PointsKey(0, eBoundaryLayerPointsRev),
            std::bind(&GaussPoints::CreateGPMatrix, this, pl::_1));
    }

protected:
    LIB_UTILITIES_EXPORT virtual const std::shared_ptr<NekMatrix<NekDouble>> v_GetI(
        const PointsKey &pkey) override;
    LIB_UTILITIES_EXPORT virtual const std::shared_ptr<NekMatrix<NekDouble>> v_GetI(
        const Array<OneD, const NekDouble> &x) override;
    LIB_UTILITIES_EXPORT virtual const std::shared_ptr<NekMatrix<NekDouble>> v_GetI(
        size_t numpoints, const Array<OneD, const NekDouble> &x) override;
    LIB_UTILITIES_EXPORT virtual const std::shared_ptr<NekMatrix<NekDouble>>
    v_GetGalerkinProjection(const PointsKey &pkey) override;

private:
    static bool initPointsManager[];

    GaussPoints()                          = delete;
    GaussPoints(const GaussPoints &points) = delete;

    virtual void v_CalculatePoints() override final;
    virtual void v_CalculateWeights() override final;
    virtual void v_CalculateDerivMatrix() override final;

    void CalculateInterpMatrix(size_t npts,
                               const Array<OneD, const NekDouble> &xpoints,
                               Array<OneD, NekDouble> &interp);

    std::shared_ptr<NekMatrix<NekDouble>> CalculateGalerkinProjectionMatrix(
        const PointsKey &pkey);

}; // class GaussPoints
} // namespace LibUtilities
} // namespace Nektar

#endif // GAUSSPOINTS_H
