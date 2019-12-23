///////////////////////////////////////////////////////////////////////////////
//
// File SIACFilter.cpp
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
// Description: SIACFilter definition
//
///////////////////////////////////////////////////////////////////////////////
#include "SIACFilter.h"
#include "GeneralBSplines.h"
#include <Eigen/Dense>
#include <LibUtilities/Foundations/ManagerAccess.h> // for points Manager,etc

namespace Nektar
{
namespace LSIAC
{

void SIACFilter::CalCoeffForStandardKernel(
    int deg, Nektar::Array<OneD, NekDouble> &coeffs, const NekDouble shift)
{
    Nektar::Array<OneD, NekDouble> tknots(3 * deg + 2);
    for (int i = 0; i < tknots.num_elements(); i++)
    {
        tknots[i] = -(deg + 1.0) / 2.0 - deg + i + shift;
    }
    CalCoeffForWideSymKernel(deg, tknots, coeffs, shift);
}

void SIACFilter::CalCoeffForKnotMatrixVec(
    int deg, const std::vector<std::vector<NekDouble>> &kMatrix,
    Nektar::Array<OneD, NekDouble> &coeffs)
{
    // Note shift is dummy variable in this function.
    int nBSpl = kMatrix.size();
    Eigen::MatrixXd M0(nBSpl, nBSpl); // deg + jj+1, deg + jj+1);
    Eigen::VectorXd e1(nBSpl);        // (deg + jj+1);
    for (int r = 0; r < nBSpl; r++)
    {
        for (int l = 0; l < nBSpl; l++)
        {
            int delta = r;
            int gamma = l;
            M0(r, l) =
                dividedDiff(0, deg + 1, kMatrix[gamma], 0.0, deg + delta + 1);
        }
    }
    for (int i = 0; i < nBSpl; i++)
    {
        e1(i) = 0.0;
    }
    e1(0)              = 1.0;
    Eigen::VectorXd x  = M0.colPivHouseholderQr().solve(e1);
    Eigen::VectorXd x1 = M0.colPivHouseholderQr().solve(e1);
    Eigen::VectorXd x2 = M0.colPivHouseholderQr().solve(e1 - M0 * x1);
    x                  = x1 + x2;
    for (int i = 0; i < nBSpl; i++)
    {
        coeffs[i] = x(i);
    }
}

void SIACFilter::CalCoeffForWideSymKernel(
    const int deg, const int der, const int nSpl,
    Nektar::Array<OneD, NekDouble> &coeffs, const NekDouble shift)
{
    Nektar::Array<OneD, NekDouble> tknots(nSpl + deg - der + 1);
    for (int i = 0; i < tknots.num_elements(); i++)
    {
        tknots[i] = -(deg + 1.0 - der) / 2.0 - (nSpl - 1.0) / 2.0 + i + shift;
    }
    CalCoeffForWideSymKernel(deg, tknots, coeffs, shift);
}

void SIACFilter::CalCoeffForWideSymKernel(
    const int deg, const int nSpl, Nektar::Array<OneD, NekDouble> &coeffs,
    const NekDouble shift)
{
    Nektar::Array<OneD, NekDouble> tknots(nSpl + deg + 1);
    for (int i = 0; i < tknots.num_elements(); i++)
    {
        tknots[i] = -(deg + 1.0) / 2.0 - (nSpl - 1.0) / 2.0 + i + shift;
    }
    CalCoeffForWideSymKernel(deg, tknots, coeffs, shift);
}

void SIACFilter::CalCoeffForWideSymKernel(
    int deg, const Nektar::Array<OneD, NekDouble> &tknots,
    Nektar::Array<OneD, NekDouble> &coeffs, const NekDouble shift)
{
    boost::ignore_unused(shift); // does not not depend on shift
    // Note shift is dummy variable in this function.
    int jj = 0;
    jj     = tknots.num_elements() - 2 * deg - 2;
    Eigen::MatrixXd M0(deg + jj + 1, deg + jj + 1);
    Eigen::VectorXd e1(deg + jj + 1);
    for (int r = 0; r < deg + jj + 1; r++)
    {
        for (int l = 0; l < deg + jj + 1; l++)
        {
            int delta = r;
            int gamma = l;
            M0(r, l)  = dividedDiff(gamma, gamma + deg + 1, tknots, 0.0,
                                   deg + delta + 1);
        }
    }
    for (int i = 0; i < deg + jj + 1; i++)
    {
        e1(i) = 0.0;
    }
    e1(0)             = 1.0;
    Eigen::VectorXd x = M0.colPivHouseholderQr().solve(e1);

    // Trying to increase the precision
    Eigen::VectorXd x1 = M0.colPivHouseholderQr().solve(e1);
    Eigen::VectorXd x2 = M0.colPivHouseholderQr().solve(e1 - M0 * x1);
    x                  = x1 + x2;

    for (int i = 0; i < deg + jj + 1; i++)
    {
        coeffs[i] = x(i);
    }
}

NekDouble SIACFilter::dividedDiff(int ii, int jj,
                                  const std::vector<NekDouble> &tknots,
                                  NekDouble x, NekDouble D)
{
    NekDouble result = 0.0;
    if (ii == jj)
    {
        result = Monomial(tknots[ii], x, 0, D);
    }
    else if (abs(tknots[ii] - tknots[jj]) < TOLERENCE)
    {
        result = Monomial(tknots[ii], x, jj - ii, D) /
                 ((NekDouble)factorial(jj - ii));
    }
    else
    {
        result = (dividedDiff(ii + 1, jj, tknots, x, D) -
                  dividedDiff(ii, jj - 1, tknots, x, D)) /
                 (tknots[jj] - tknots[ii]);
    }
    return result;
}

NekDouble SIACFilter::dividedDiff(int ii, int jj,
                                  const Nektar::Array<OneD, NekDouble> &tknots,
                                  NekDouble x, NekDouble D)
{
    NekDouble result = 0.0;
    if (ii == jj)
    {
        result = Monomial(tknots[ii], x, 0, D);
    }
    else if (abs(tknots[ii] - tknots[jj]) < TOLERENCE)
    {
        result = Monomial(tknots[ii], x, jj - ii, D) /
                 ((NekDouble)factorial(jj - ii));
    }
    else
    {
        result = (dividedDiff(ii + 1, jj, tknots, x, D) -
                  dividedDiff(ii, jj - 1, tknots, x, D)) /
                 (tknots[jj] - tknots[ii]);
    }
    return result;
}

NekDouble SIACFilter::Monomial(NekDouble t, NekDouble x, int der, int degree)
{
    NekDouble y      = t - x;
    NekDouble result = pow(y, degree - der) * ((NekDouble)factorial(degree)) /
                       ((NekDouble)factorial(degree - der));
    return result;
}

int SIACFilter::factorial(int jj)
{
    int fact = 1;
    for (int i = 1; i <= jj; i++)
    {
        fact *= i;
    }
    return fact;
}

void SIACFilter::CalCoeffForCenBSplDerivatives(const int degree,
                                               const int alpha, int &nBSpl,
                                               Array<OneD, NekDouble> &coeffs)
{
    boost::ignore_unused(degree); // does not depend on degree;
    int tnBSpl = nBSpl - alpha;
    for (int a = 1; a <= alpha; a++)
    {
        NekDouble temp1 = coeffs[0];
        NekDouble temp2 = 0.0;
        for (int i = 1; i < tnBSpl; i++)
        {
            temp2     = coeffs[i];
            coeffs[i] = -1 * temp1 + temp2;
            temp1     = temp2;
        }
        coeffs[tnBSpl] = -1 * temp1;
        tnBSpl         = tnBSpl + 1;
    }
}

void SIACFilter::CalCoeffForKnotMatrixVec_Hanieh(
    int deg, const std::vector<std::vector<NekDouble>> &kMatrix,
    Nektar::Array<OneD, NekDouble> &coeffs)
{
    // Note shift is dummy variable in this function.
    int nBSpl = kMatrix.size();
    Eigen::MatrixXd M0(nBSpl, nBSpl); // deg + jj+1, deg + jj+1);
    Eigen::VectorXd e1(nBSpl);        // (deg + jj+1);
    for (int r = 0; r < nBSpl; r++)
    {
        for (int l = 0; l < nBSpl; l++)
        {
            int quad_npoints = std::ceil((deg + 1 + r + 3) / 2);
            LibUtilities::PointsKey quadPointsKey(
                quad_npoints, LibUtilities::eGaussGaussLegendre);
            Array<OneD, NekDouble> vals(quad_npoints, 0.0),
                xpos(quad_npoints, 0.0);
            Array<OneD, NekDouble> quad_points =
                LibUtilities::PointsManager()[quadPointsKey]->GetZ();
            Array<OneD, NekDouble> quad_weights =
                LibUtilities::PointsManager()[quadPointsKey]->GetW();
            M0(r, l) = 0;
            GeneralBSplines gSPl(deg + 1);
            for (int i = 0; i < kMatrix[l].size() - 1; i++)
            {
                NekDouble bgamma1 = kMatrix[l][i];
                NekDouble bgamma2 = kMatrix[l][i + 1];
                for (int i = 0; i < quad_npoints; i++)
                {
                    xpos[i] = bgamma1 +
                              (quad_points[i] + 1) * (bgamma2 - bgamma1) / 2.0;
                }
                gSPl.EvaluateBSplines(xpos, kMatrix[l], 0, vals);
                NekDouble sum = 0;
                for (int i = 0; i < quad_npoints; i++)
                {
                    sum += quad_weights[i] * vals[i] * std::pow(xpos[i], r) *
                           std::abs(bgamma2 - bgamma1) / 2.0;
                }
                M0(r, l) += sum;
            }
            // Calcualte the matrix given polynomial order and monomial order.
        }
    }
    for (int i = 0; i < nBSpl; i++)
    {
        e1(i) = 0.0;
    }
    e1(0)             = 1.0;
    Eigen::VectorXd x = M0.colPivHouseholderQr().solve(e1);
    Eigen::VectorXd x1 = M0.colPivHouseholderQr().solve(e1);
    Eigen::VectorXd x2 = M0.colPivHouseholderQr().solve(e1 - M0 * x1);
    x                  = x1 + x2;
    for (int i = 0; i < nBSpl; i++)
    {
        coeffs[i] = x(i);
    }
}

bool AreVectorsEqual(vector<double> t1, vector<double> t2)
{
    if (t1.size() != t2.size())
        return false;
    for (int i = 0; i < t1.size(); i++)
    {
        if (TOLERENCE < std::abs(t1[i] - t2[i]))
            return false;
    }
    return true;
}

void AddVectoMatrix(vector<double> &tv, vector<vector<double>> &temp)
{
    bool t = false;
    for (int i = 0; i < temp.size(); i++)
    {
        t = t || AreVectorsEqual(tv, temp[i]);
    }
    if (!t)
        temp.push_back(tv);
}

void calculateDerivativeKnotMatrix(const std::vector<std::vector<NekDouble>> &M,
                                   std::vector<std::vector<NekDouble>> &Result)
{
    vector<vector<NekDouble>> temp;
    temp.clear();
    for (int i = 0; i < M.size(); i++)
    {
        vector<NekDouble> tv1(M[i].begin(), M[i].end() - 1);
        vector<NekDouble> tv2(M[i].begin() + 1, M[i].end());
        AddVectoMatrix(tv1, temp);
        AddVectoMatrix(tv2, temp);
    }
    Result = temp;
}

bool AccumlateCoeffForMatrix(double coeff, vector<double> &tv,
                             vector<vector<double>> &mat,
                             vector<double> &coeffsArray)
{
    bool ret = false;
    for (int i = 0; i < mat.size(); i++)
    {
        if (AreVectorsEqual(tv, mat[i]))
        {
            coeffsArray[i] += coeff;
            ret = true;
            break;
        }
    }
    return ret;
}

void SIACFilter::CalDerivativesForKnotMatrixVec_Hanieh(
    const int deg, const int derivative,
    std::vector<std::vector<NekDouble>> &kMatrix,
    Nektar::Array<OneD, NekDouble> &coeffs)
{
    boost::ignore_unused(deg); // doesnot depend on degree;
    // loop through number of derivative times.
    // calcualte a new derivative matrix.
    // loop through knotmatrix
    // calculate knots and coefficients.
    // Add coefficients into matrix.
    // calculate a new coefficient matrix.
    // Replace original matrix by derivative matrix.
    // Replace old coefficient by new coffecients.
    std::vector<std::vector<NekDouble>> tempIn, tempOut;
    tempIn = kMatrix;
    std::vector<NekDouble> tempIn_coeffs(coeffs.data(),
                                         coeffs.data() + tempIn.size());
    for (int d = 1; d <= derivative; d++)
    {
        std::vector<NekDouble> tempOut_coeffs(tempIn.size() + 1, 0.0);
        calculateDerivativeKnotMatrix(tempIn, tempOut);
        for (int k = 0; k < tempIn.size(); k++)
        {
            int deg_v = tempIn[0].size() - 2;
            // vector under derivative is tempIn[k]
            // Cal coefficients.
            std::vector<NekDouble> tv1(tempIn[k].begin(), tempIn[k].end() - 1);
            std::vector<NekDouble> tv2(tempIn[k].begin() + 1, tempIn[k].end());
            if (TOLERENCE < std::abs(tv1.back() - tv1.front()))
            {
                NekDouble coeff1 = tempIn_coeffs[k] * (NekDouble)(deg_v) /
                                   (tv1.back() - tv1.front());
                AccumlateCoeffForMatrix(coeff1, tv1, tempOut, tempOut_coeffs);
            }
            if (TOLERENCE < std::abs(tv2.back() - tv2.front()))
            {
                NekDouble coeff2 = -1.0 * tempIn_coeffs[k] *
                                   (NekDouble)(deg_v) /
                                   (tv2.back() - tv2.front());
                AccumlateCoeffForMatrix(coeff2, tv2, tempOut, tempOut_coeffs);
            }
        }
        tempIn        = tempOut;
        tempIn_coeffs = tempOut_coeffs;
    }
    kMatrix = tempIn;
    coeffs  = Array<OneD, NekDouble>(tempIn_coeffs.size());
    for (int i = 0; i < tempIn_coeffs.size(); i++)
    {
        coeffs[i] = tempIn_coeffs[i];
    }
}

} // namespace LSIAC
} // namespace Nektar
