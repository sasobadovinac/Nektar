///////////////////////////////////////////////////////////////////////////////
//
// File MetricTensor.h
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
#pragma once

#include "HandleNekMesh.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

namespace Nektar
{
namespace LSIAC
{

// forward declaration
class HandleNekMesh;
/// This subclasses of the class can be used to postprocess any point on a given
/// mesh.
/** This class and its subclasses enable users to post process different points
   the mesh. These classes contain methods which are used by the end user of the
   tool.
*/
class MetricTensor : public LSIACPostProcessor
{
private:
    bool CalculateDirection(SpatialDomains::PointGeomSharedPtr v0,
                            SpatialDomains::PointGeomSharedPtr v1,
                            Array<OneD, NekDouble> &eig);

protected:
public:
    MetricTensor();

    // First iteration Everthing under public
    // 	1.	No of  Quadrature points.
    // 		2.	Location of Quadrature points.
    // 			3.	Bool for 1D,2D,3D
    // 				4.	Array<OneD,NekDouble>  Ne* Nq* Nm
    // 					5.	Array<OneD, NekDouble> 
    // Ne*Nq*Neig1 						6.	Array<OneD,
    // NekDouble> 
    // Ne*Nq*Neig2 							7.
    // Array<OneD, NekDouble>
    // Ne*Nq*Nlam 								8.
    // Array<OneD, NekDouble> 
    // Ne*Nq*Neig3 									9.
    // Array<OneD,NekDouble>  Metric tensor
    int m_nOfQPE; // Assumption that there only one type of element and all ...
                  //  of them have same number of quadrature points.
    Array<OneD, NekDouble> m_locOfquadPts; // in reference space.
    unsigned int dimension;                // 1,2,3 corresponds to dimension
    bool m_constTensor; // If true, We store one Metric tensor per element.
                        // Hence ... Discontinous .i.e m_noOfQpts=1;
    std::vector<Eigen::Matrix2d> m_metricTensorAtElm; // Ne*Nq*Nm;
    std::map<int, Eigen::Matrix2d> m_metricTensorLogAtNode;
    Array<OneD, NekDouble> m_eigenValue1; // Ne*Nq*3;
    Array<OneD, NekDouble> m_eigenValue2; // Ne*Nq*3;
    Array<OneD, NekDouble> m_eigenValue3; // Ne*Nq*3;
    Array<OneD, NekDouble> m_lambda1;     // Ne*Nq
    Array<OneD, NekDouble> m_lambda2;     // Ne*Nq
    Array<OneD, NekDouble> m_lambda3;     // Ne*Nq
    HandleNekMesh *m_meshHandlePtr;

    bool LoadMetricTensor(HandleNekMesh *HNM);

    Eigen::Matrix2d GetDynamicMetricTensor(Array<OneD, NekDouble> glCoord,
                                           int eid);

    bool CalculateMetricTensorAtNodes();
    // What functions do we need ?
    // Input: Given a point, eld ID?
    // Input: Gven a point only ?
    // Output: Get only Lambda1 and eigen1.
    // Output: Get only lambda2 and eigen2.
    // Output: Get only scaling given a direction ?
    bool GetEigenPair(Array<OneD, NekDouble> coord, int EigNum,
                      NekDouble &lambda, Array<OneD, NekDouble> &eigen) const;
    bool GetEigenPair(Array<OneD, NekDouble> coord, int eid, int EigNum,
                      NekDouble &lambda, Array<OneD, NekDouble> &eigen) const;

    bool GetEigenPairUsingIP(Array<OneD, NekDouble> coord, int EigNum,
                             NekDouble &lambda, Array<OneD, NekDouble> &eigen);
    bool GetEigenPairUsingIP(Array<OneD, NekDouble> coord, int eid, int EigNum,
                             NekDouble &lambda, Array<OneD, NekDouble> &eigen);

    bool GetEigenPairAtTheta(Array<OneD, NekDouble> coord,
                             NekDouble theta_degrees, NekDouble &lambda,
                             Array<OneD, NekDouble> &eigen) const;
    bool GetEigenPairAtTheta(Array<OneD, NekDouble> coord, int eid,
                             NekDouble theta_degrees, NekDouble &lambda,
                             Array<OneD, NekDouble> &eigen) const;

    bool GetEigenPairAtTheta(int eid, NekDouble theta_degrees,
                             NekDouble &lambda,
                             Array<OneD, NekDouble> &eigen) const;
    bool GetScaleForDirection(int eid, Array<OneD, NekDouble> direction,
                              NekDouble &lambda) const;

    bool GetScaling(Array<OneD, NekDouble> coord, int eigNum,
                    NekDouble &scaling) const;
    bool GetScaling(Array<OneD, NekDouble> coord, int eid, int eigNum,
                    NekDouble &scaling) const;
};

} // namespace LSIAC
} // namespace Nektar
