///////////////////////////////////////////////////////////////////////////////
//
// File MetricTensor.cpp
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
// Description: MetricTensor definition
//
///////////////////////////////////////////////////////////////////////////////
#include "MetricTensor.h"
#include <math.h>

#define PI 3.14159265

namespace Nektar
{
namespace LSIAC
{

MetricTensor::MetricTensor()
{
}

Eigen::Matrix2d MetricTensor::GetDynamicMetricTensor(
    Array<OneD, NekDouble> glCoord, int eid)
{
    // if eid <0 find a elid.
    if (eid < 0)
    {
        eid = m_meshHandlePtr->GetExpansionIndexUsingRTree(glCoord);
        ASSERTL0(eid >= 0 && "Point out of mesh");
    }
    ASSERTL0(eid >= 0 && "Point out of mesh");
    Eigen::Matrix2d result;
    // Get local coordinates.
    // Depending on number of vertices triangle or quad.
    // use locCoordinates as barycentric coordinates
    SpatialDomains::GeometrySharedPtr geomSPtr =
        m_meshHandlePtr->m_expansions[0]->GetExp(eid)->GetGeom();
    Array<OneD, NekDouble> lCoord(3, 0.0);
    geomSPtr->GetLocCoords(glCoord, lCoord);
    if (geomSPtr->GetShapeType() == Nektar::LibUtilities::eTriangle)
    {
        int Vid0          = geomSPtr->GetVid(0);
        int Vid1          = geomSPtr->GetVid(1);
        int Vid2          = geomSPtr->GetVid(2);
        NekDouble lambda1 = (lCoord[0] + 1.0) / 2.0;
        NekDouble lambda2 = (lCoord[1] + 1.0) / 2.0;
        NekDouble lambda0 = 1.0 - lambda1 - lambda2;
        Eigen::Matrix2d temp0 = m_metricTensorLogAtNode[Vid0];
        Eigen::Matrix2d temp1 = m_metricTensorLogAtNode[Vid1];
        Eigen::Matrix2d temp2 = m_metricTensorLogAtNode[Vid2];
        result = lambda0 * temp0 + lambda1 * temp1 + lambda2 * temp2;
    }
    else if (geomSPtr->GetShapeType() == Nektar::LibUtilities::eQuadrilateral)
    {
        ASSERTL0("Not designed for quadrilateral elements");
    }
    else
    {
        ASSERTL0("Shape not accounted for");
    }

    return result.exp();
}

bool MetricTensor::CalculateMetricTensorAtNodes()
{
    // Can calculate for only triangles as of now.
    std::map<int, Eigen::Matrix2d> MT_m_a;
    std::map<int, NekDouble> MT_a;
    for (size_t eid = 0; eid < m_meshHandlePtr->m_expansions[0]->GetExpSize();
         eid++)
    {
        SpatialDomains::GeometrySharedPtr geomSPtr =
            m_meshHandlePtr->m_expansions[0]->GetExp(eid)->GetGeom();
        Eigen::Matrix2d temp = m_metricTensorAtElm[eid];
        NekDouble a          = m_meshHandlePtr->GetJacobian(eid);
        for (size_t vid = 0; vid < geomSPtr->GetNumVerts(); vid++)
        {
            int Vid = geomSPtr->GetVid(vid);
            if (MT_m_a.find(Vid) == MT_m_a.end())
            {
                MT_m_a.insert(std::make_pair(Vid, a * (temp.log())));
                MT_a.insert(std::make_pair(Vid, a));
            }
            else
            {
                MT_m_a.find(Vid)->second += temp.log() * a;
                MT_a.find(Vid)->second += a;
            }
        }
    }

    // Area weighted at vertices.
    for (std::map<int, Eigen::Matrix2d>::iterator it = MT_m_a.begin();
         it != MT_m_a.end(); it++)
    {
        NekDouble totalArea  = MT_a.find(it->first)->second;
        Eigen::Matrix2d temp = (MT_m_a.find(it->first)->second) / totalArea;
        m_metricTensorLogAtNode.insert(std::make_pair(it->first, temp));
    }
    return true;
}
bool MetricTensor::CalculateDirection(SpatialDomains::PointGeomSharedPtr v0,
                                      SpatialDomains::PointGeomSharedPtr v1,
                                      Array<OneD, NekDouble> &eig)
{
    SpatialDomains::PointGeom vec;
    vec.Sub(*v0, *v1);
    NekDouble mag = v0->dist(*v1);
    vec.GetCoords(eig);
    eig[0] = eig[0] / mag;
    eig[1] = eig[1] / mag;
    eig[2] = eig[2] / mag;

    return true;
}
void setEigenVector(NekDouble &eig1, NekDouble &eig2)
{
    if (eig2 < 0)
    {
        eig1 = -1 * eig1;
        eig2 = -1 * eig2;
    }
}

bool MetricTensor::LoadMetricTensor(HandleNekMesh *HNM)
{
    // Loop through all the elements find the largest and smallest eigen values.
    // For quadrilateral mesh.
    // Assuming it is coordinate aligned or all angles are 90 degrees.
    // Find largest and smallest edges.
    // All the largest edge length are lambda 1 and direction is eigen 1.
    // All the smallest edge length are lambda2 and direction is eigen 2.
    m_nOfQPE        = 1;
    m_meshHandlePtr = HNM;
    int numElm      = HNM->m_expansions[0]->GetExpSize();
    // initialize the values.
    m_metricTensorAtElm.resize(numElm); // = Array<OneD,NekDouble>(numElm*4); //
                                        // since 4 elements for each tensor.
    m_eigenValue1 = Array<OneD, NekDouble>(m_nOfQPE * 3 * numElm);
    m_eigenValue2 = Array<OneD, NekDouble>(m_nOfQPE * 3 * numElm);
    m_lambda1     = Array<OneD, NekDouble>(m_nOfQPE * numElm);
    m_lambda2     = Array<OneD, NekDouble>(m_nOfQPE * numElm);
    Array<OneD, NekDouble> eig1(3, 0.0);
    Array<OneD, NekDouble> eig2(3, 0.0);

    // this is for only quad meshes with right angles.
    for (int eid = 0; eid < HNM->m_expansions[0]->GetExpSize(); eid++)
    {
        SpatialDomains::GeometrySharedPtr geomSPtr =
            HNM->m_expansions[0]->GetExp(eid)->GetGeom();
        switch (geomSPtr->GetShapeType())
        {
            case (LibUtilities::eQuadrilateral):
            {
                NekDouble maxV = 0.0, minV = 1e100;
                for (int edgeid = 0; edgeid < geomSPtr->GetNumEdges(); edgeid++)
                {
                    NekDouble edgeLength =
                        HNM->m_segMap[geomSPtr->GetEid(edgeid)]
                            ->GetVertex(0)
                            ->dist(*(HNM->m_segMap[geomSPtr->GetEid(edgeid)]
                                         ->GetVertex(1)));
                    if (maxV < edgeLength)
                    {
                        maxV = edgeLength;
                        CalculateDirection(
                            HNM->m_segMap[geomSPtr->GetEid(edgeid)]->GetVertex(
                                0),
                            HNM->m_segMap[geomSPtr->GetEid(edgeid)]->GetVertex(
                                1),
                            eig1);
                    }
                    if (minV > edgeLength)
                    {
                        minV = edgeLength;
                        CalculateDirection(
                            HNM->m_segMap[geomSPtr->GetEid(edgeid)]->GetVertex(
                                0),
                            HNM->m_segMap[geomSPtr->GetEid(edgeid)]->GetVertex(
                                1),
                            eig2);
                    }
                }
                if (((std::abs(eig1[0] - eig2[0]) < 1e-8) &&
                     (std::abs(eig1[1] - eig2[1]) < 1e-8)) ||
                    ((std::abs(eig1[0] + eig2[0]) < 1e-8) &&
                     (std::abs(eig1[1] + eig2[1]) < 1e-8)))
                { // if both vectors are same direction.
                    // Works only for 2D.
                    eig2[0] = -1.0 * eig1[1];
                    eig2[1] = eig1[0];
                }
                if (eig1[1] < 0)
                {
                    eig1[0] = -1.0 * eig1[0];
                    eig1[1] = -1.0 * eig1[1];
                }
                if (eig2[1] < 0)
                {
                    eig2[0] = -1.0 * eig2[0];
                    eig2[1] = -1.0 * eig2[1];
                }
                if (std::abs(eig1[0] + 1) < 1e-8)
                {
                    eig1[0] = 1;
                }
                if (std::abs(eig2[0] + 1) < 1e-8)
                {
                    eig2[0] = 1;
                }
                // Store the calculated values.
                m_lambda1[m_nOfQPE * eid]             = maxV;
                m_lambda2[m_nOfQPE * eid]             = minV;
                m_eigenValue1[m_nOfQPE * 3 * eid + 0] = eig1[0];
                m_eigenValue2[m_nOfQPE * 3 * eid + 0] = eig2[0];
                m_eigenValue1[m_nOfQPE * 3 * eid + 1] = eig1[1];
                m_eigenValue2[m_nOfQPE * 3 * eid + 1] = eig2[1];
                m_eigenValue1[m_nOfQPE * 3 * eid + 2] = eig1[2];
                m_eigenValue2[m_nOfQPE * 3 * eid + 2] = eig2[2];
            }
            break;
            case (LibUtilities::eTriangle):
            {
                SpatialDomains::PointGeomSharedPtr aANek =
                    geomSPtr->GetVertex(0);
                SpatialDomains::PointGeomSharedPtr aBNek =
                    geomSPtr->GetVertex(1);
                SpatialDomains::PointGeomSharedPtr aCNek =
                    geomSPtr->GetVertex(2);
                // Using eigen library to calculate metric;
                Eigen::Vector2d aA, aB, aC, rA, rB, rC, eA, eB, eC;
                rA << 1, 0;
                rB << 0, 0;
                rC << 0, 1;
                eA << 0, sqrt(3.0 / 4.0);
                eB << 0.5, 0.0;
                eC << -0.5, 0.0;
                NekDouble aA_x, aA_y, aA_z, aB_x, aB_y, aB_z, aC_x, aC_y, aC_z;
                aANek->GetCoords(aA_x, aA_y, aA_z);
                aBNek->GetCoords(aB_x, aB_y, aB_z);
                aCNek->GetCoords(aC_x, aC_y, aC_z);
                aA << aA_x, aA_y;
                aB << aB_x, aB_y;
                aC << aC_x, aC_y;
                Eigen::Matrix2d Tr_e, Tr_a;

                Tr_e(0, 0) = eA(0) - eC(0);
                Tr_e(1, 0) = eA(1) - eC(1);
                Tr_e(0, 1) = eB(0) - eC(0);
                Tr_e(1, 1) = eB(1) - eC(1);

                Tr_a(0, 0) = aA(0) - aC(0);
                Tr_a(1, 0) = aA(1) - aC(1);
                Tr_a(0, 1) = aB(0) - aC(0);
                Tr_a(1, 1) = aB(1) - aC(1);

                Eigen::Matrix2d Te_a =
                    Tr_e.transpose().inverse() * (Tr_a.transpose());
                Te_a.transposeInPlace();

                Eigen::Matrix2d Ta_e =
                    Tr_a.transpose().inverse() * (Tr_e.transpose());

                Eigen::Matrix2d Metric = Ta_e * Ta_e.transpose();

                Eigen::EigenSolver<Eigen::Matrix2d> es(Metric);
                m_metricTensorAtElm[eid] = Metric.pow(-0.5);

                if ((es.eigenvectors().col(0)(0).real()) <
                    (es.eigenvectors().col(0)(1).real()))
                {
                    m_lambda1[m_nOfQPE * eid] =
                        1 / std::sqrt(es.eigenvalues().col(0)(0).real());
                    m_lambda2[m_nOfQPE * eid] =
                        1 / std::sqrt(es.eigenvalues().col(0)(1).real());
                    m_eigenValue1[m_nOfQPE * 3 * eid + 0] =
                        es.eigenvectors().col(0)(0).real();
                    m_eigenValue1[m_nOfQPE * 3 * eid + 1] =
                        es.eigenvectors().col(0)(1).real();
                    m_eigenValue1[m_nOfQPE * 3 * eid + 2] = 0.0;

                    m_eigenValue2[m_nOfQPE * 3 * eid + 0] =
                        es.eigenvectors().col(1)(0).real();
                    m_eigenValue2[m_nOfQPE * 3 * eid + 1] =
                        es.eigenvectors().col(1)(1).real();
                    m_eigenValue2[m_nOfQPE * 3 * eid + 2] = 0.0;
                }
                else
                {
                    m_lambda1[m_nOfQPE * eid] =
                        1 / std::sqrt(es.eigenvalues().col(0)(1).real());
                    m_lambda2[m_nOfQPE * eid] =
                        1 / std::sqrt(es.eigenvalues().col(0)(0).real());

                    m_eigenValue1[m_nOfQPE * 3 * eid + 0] =
                        es.eigenvectors().col(1)(0).real();
                    m_eigenValue1[m_nOfQPE * 3 * eid + 1] =
                        es.eigenvectors().col(1)(1).real();
                    m_eigenValue1[m_nOfQPE * 3 * eid + 2] = 0.0;

                    m_eigenValue2[m_nOfQPE * 3 * eid + 0] =
                        es.eigenvectors().col(0)(0).real();
                    m_eigenValue2[m_nOfQPE * 3 * eid + 1] =
                        es.eigenvectors().col(0)(1).real();
                    m_eigenValue2[m_nOfQPE * 3 * eid + 2] = 0.0;
                }
            }
            break;
            case (LibUtilities::eSegment):
            {
                SpatialDomains::PointGeomSharedPtr aANek =
                    geomSPtr->GetVertex(0);
                SpatialDomains::PointGeomSharedPtr aBNek =
                    geomSPtr->GetVertex(1);
                NekDouble aA_x, aA_y, aA_z, aB_x, aB_y, aB_z;
                aANek->GetCoords(aA_x, aA_y, aA_z);
                aBNek->GetCoords(aB_x, aB_y, aB_z);
                // Assming all the elements are in x-direction here.
                m_eigenValue1[m_nOfQPE * 3 * eid + 0] = 1.0;
                m_eigenValue1[m_nOfQPE * 3 * eid + 1] = 0.0;
                m_eigenValue1[m_nOfQPE * 3 * eid + 2] = 0.0;
                m_lambda1[m_nOfQPE * eid]             = std::abs(aA_x - aB_x);
            }
            break; // some code in curly braces to help with cross-initializtion
                   // of variables such as aNEK.
            default:
            {
                NEKERROR(ErrorUtil : efatal, "Not accounted for");
            }
            break;
        }
    }
    return true;
}

bool MetricTensor::GetEigenPair(Array<OneD, NekDouble> coord, int eigNum,
                                NekDouble &lambda,
                                Array<OneD, NekDouble> &eigen) const
{
    int eid = -1;
    return GetEigenPair(coord, eid, eigNum, lambda, eigen);
}
bool MetricTensor::GetEigenPair(Array<OneD, NekDouble> coord, int eid,
                                int eigNum, NekDouble &lambda,
                                Array<OneD, NekDouble> &eigen) const
{
    if (eid == -1)
    {
        // Find eid.
        eid = m_meshHandlePtr->GetExpansionIndexUsingRTree(coord);
        ASSERTL0(eid >= 0 && "Did not initialize Rtrees");
    }
    // Logically this should have been a quadrature point.
    // For now assume 1.
    int quadId = 0;
    switch (eigNum)
    {
        case 1:
            lambda   = m_lambda1[eid * m_nOfQPE + quadId];
            eigen[0] = m_eigenValue1[eid * m_nOfQPE * 3 + quadId * 3 + 0];
            eigen[1] = m_eigenValue1[eid * m_nOfQPE * 3 + quadId * 3 + 1];
            eigen[2] = m_eigenValue1[eid * m_nOfQPE * 3 + quadId * 3 + 2];
            break;
        case 2:
            lambda   = m_lambda2[eid * m_nOfQPE + quadId];
            eigen[0] = m_eigenValue2[eid * m_nOfQPE * 3 + quadId * 3 + 0];
            eigen[1] = m_eigenValue2[eid * m_nOfQPE * 3 + quadId * 3 + 1];
            eigen[2] = m_eigenValue2[eid * m_nOfQPE * 3 + quadId * 3 + 2];
            break;
        case 3:
            lambda   = m_lambda3[eid * m_nOfQPE + quadId];
            eigen[0] = m_eigenValue3[eid * m_nOfQPE * 3 + quadId * 3 + 0];
            eigen[1] = m_eigenValue3[eid * m_nOfQPE * 3 + quadId * 3 + 1];
            eigen[2] = m_eigenValue3[eid * m_nOfQPE * 3 + quadId * 3 + 2];
            break;
        default:
            NEKERROR(ErrorUtil : efatal, "Wrong input value");
            return false;
    }
    return true;
}

bool MetricTensor::GetEigenPairAtTheta(Array<OneD, NekDouble> coord,
                                       NekDouble theta_degrees,
                                       NekDouble &lambda,
                                       Array<OneD, NekDouble> &eigen) const
{
    int eid = -1;
    return GetEigenPairAtTheta(coord, eid, theta_degrees, lambda, eigen);
}

bool MetricTensor::GetScaleForDirection(int eid,
                                        Array<OneD, NekDouble> direction,
                                        NekDouble &lambda) const
{
    ASSERTL0(eid >= 0 && "Not valid element id");
    ASSERTL0(std::abs(direction[0] * direction[0] +
                      direction[1] * direction[1] - 1) < 1e-9 &&
             "direction is not normallized");
    int quadId = 0;
    // Assuming only 2D for now.
    NekDouble lambda1, lambda2;
    Array<OneD, NekDouble> eigen1(3, 0.0), eigen2(3, 0.0);
    lambda1   = m_lambda1[eid * m_nOfQPE + quadId];
    eigen1[0] = m_eigenValue1[eid * m_nOfQPE * 3 + quadId * 3 + 0];
    eigen1[1] = m_eigenValue1[eid * m_nOfQPE * 3 + quadId * 3 + 1];
    eigen1[2] = m_eigenValue1[eid * m_nOfQPE * 3 + quadId * 3 + 2];

    lambda2   = m_lambda2[eid * m_nOfQPE + quadId];
    eigen2[0] = m_eigenValue2[eid * m_nOfQPE * 3 + quadId * 3 + 0];
    eigen2[1] = m_eigenValue2[eid * m_nOfQPE * 3 + quadId * 3 + 1];
    eigen2[2] = m_eigenValue2[eid * m_nOfQPE * 3 + quadId * 3 + 2];

    NekDouble vecX =
        lambda1 * (eigen1[0] * direction[0] + eigen1[1] * direction[1]);
    NekDouble vecY =
        lambda2 * (eigen2[0] * direction[0] + eigen2[1] * direction[1]);

    lambda = sqrt(vecX * vecX + vecY * vecY);
    return true;
}
bool MetricTensor::GetEigenPairAtTheta(int eid, NekDouble theta_degrees,
                                       NekDouble &lambda,
                                       Array<OneD, NekDouble> &eigen) const
{
    ASSERTL0(eid >= 0 && "Not valid element id");
    int quadId          = 0;
    NekDouble theta_rad = theta_degrees * PI / 180.0;
    // Assuming only 2D for now.
    NekDouble lambda1, lambda2;
    Array<OneD, NekDouble> eigen1(3, 0.0), eigen2(3, 0.0);
    lambda1   = m_lambda1[eid * m_nOfQPE + quadId];
    eigen1[0] = m_eigenValue1[eid * m_nOfQPE * 3 + quadId * 3 + 0];
    eigen1[1] = m_eigenValue1[eid * m_nOfQPE * 3 + quadId * 3 + 1];
    eigen1[2] = m_eigenValue1[eid * m_nOfQPE * 3 + quadId * 3 + 2];

    lambda2   = m_lambda2[eid * m_nOfQPE + quadId];
    eigen2[0] = m_eigenValue2[eid * m_nOfQPE * 3 + quadId * 3 + 0];
    eigen2[1] = m_eigenValue2[eid * m_nOfQPE * 3 + quadId * 3 + 1];
    eigen2[2] = m_eigenValue2[eid * m_nOfQPE * 3 + quadId * 3 + 2];

    NekDouble vecX =
        lambda1 * (eigen1[0] * cos(theta_rad) + eigen1[1] * sin(theta_rad));
    NekDouble vecY =
        lambda2 * (eigen2[0] * cos(theta_rad) + eigen2[1] * sin(theta_rad));

    lambda   = sqrt(vecX * vecX + vecY * vecY);
    eigen[0] = cos(theta_rad);
    eigen[1] = sin(theta_rad);
    return true;
}

bool MetricTensor::GetEigenPairAtTheta(Array<OneD, NekDouble> coord, int eid,
                                       NekDouble theta_degrees,
                                       NekDouble &lambda,
                                       Array<OneD, NekDouble> &eigen) const
{
    if (eid == -1)
    {
        // Find eid.
        eid = m_meshHandlePtr->GetExpansionIndexUsingRTree(coord);
        ASSERTL0(eid >= 0 && "Did not initialize Rtrees");
    }
    int quadId          = 0;
    NekDouble theta_rad = theta_degrees * PI / 180.0;
    // Assuming only 2D for now.
    NekDouble lambda1, lambda2;
    Array<OneD, NekDouble> eigen1(3, 0.0), eigen2(3, 0.0);
    lambda1   = m_lambda1[eid * m_nOfQPE + quadId];
    eigen1[0] = m_eigenValue1[eid * m_nOfQPE * 3 + quadId * 3 + 0];
    eigen1[1] = m_eigenValue1[eid * m_nOfQPE * 3 + quadId * 3 + 1];
    eigen1[2] = m_eigenValue1[eid * m_nOfQPE * 3 + quadId * 3 + 2];

    lambda2   = m_lambda2[eid * m_nOfQPE + quadId];
    eigen2[0] = m_eigenValue2[eid * m_nOfQPE * 3 + quadId * 3 + 0];
    eigen2[1] = m_eigenValue2[eid * m_nOfQPE * 3 + quadId * 3 + 1];
    eigen2[2] = m_eigenValue2[eid * m_nOfQPE * 3 + quadId * 3 + 2];

    NekDouble vecX =
        lambda1 * (eigen1[0] * cos(theta_rad) + eigen1[1] * sin(theta_rad));
    NekDouble vecY =
        lambda2 * (eigen2[0] * cos(theta_rad) + eigen2[1] * sin(theta_rad));

    lambda   = sqrt(vecX * vecX + vecY * vecY);
    eigen[0] = cos(theta_rad);
    eigen[1] = sin(theta_rad);
    return false;
}

bool MetricTensor::GetScaling(Array<OneD, NekDouble> coord, int eigNum,
                              NekDouble &scaling) const
{
    int eid  = -1;
    bool ret = GetScaling(coord, eid, eigNum, scaling);
    return ret;
}

bool MetricTensor::GetScaling(Array<OneD, NekDouble> coord, int eid, int eigNum,
                              NekDouble &scaling) const
{
    boost::ignore_unused(coord, eid, eigNum, scaling);
    return false;
}

bool MetricTensor::GetEigenPairUsingIP(Array<OneD, NekDouble> coord, int eigNum,
                                       NekDouble &lambda,
                                       Array<OneD, NekDouble> &eigen)
{
    int eid = -1;
    return GetEigenPairUsingIP(coord, eid, eigNum, lambda, eigen);
}
bool MetricTensor::GetEigenPairUsingIP(Array<OneD, NekDouble> coord, int eid,
                                       int eigNum, NekDouble &lambda,
                                       Array<OneD, NekDouble> &eigen)
{

    // First calculate interpolated MetricTensor.
    Eigen::Matrix2d Metric = GetDynamicMetricTensor(coord, eid);
    // Calculate the eigen values and vectors of the Tensor.
    Eigen::EigenSolver<Eigen::Matrix2d> es(Metric);

    NekDouble lambda1, lambda2;
    Array<OneD, NekDouble> eigen1(3, 0.0), eigen2(3, 0.0);
    if ((es.eigenvectors().col(0)(0).real()) >
        (es.eigenvectors().col(0)(1).real()))
    {
        lambda1   = es.eigenvalues().col(0)(0).real();
        lambda2   = es.eigenvalues().col(0)(1).real();
        eigen1[0] = es.eigenvectors().col(0)(0).real();
        eigen1[1] = es.eigenvectors().col(0)(1).real();
        eigen2[0] = es.eigenvectors().col(1)(0).real();
        eigen2[1] = es.eigenvectors().col(1)(1).real();
    }
    else
    {
        lambda2   = es.eigenvalues().col(0)(0).real();
        lambda1   = es.eigenvalues().col(0)(1).real();
        eigen2[0] = es.eigenvectors().col(0)(0).real();
        eigen2[1] = es.eigenvectors().col(0)(1).real();
        eigen1[0] = es.eigenvectors().col(1)(0).real();
        eigen1[1] = es.eigenvectors().col(1)(1).real();
    }

    switch (eigNum)
    {
        case 1:
            lambda = lambda1;
            eigen  = eigen1;
            break;
        case 2:
            lambda = lambda2;
            eigen  = eigen2;
            break;
        default:
            ASSERTL0("wrong input for eigNUM");
            break;
    }
    return true;
}

} // namespace LSIAC
} // namespace Nektar
