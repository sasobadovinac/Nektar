////////////////////////////////////////////////////////////////////////////////
//
// File: Interpolator.cpp
//
// For more information, please see: http://www.nektar.info/
//
// The MIT License
//
// Copyright (c) 2016 Kilian Lackhove
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
// Description: Interpolator
//
////////////////////////////////////////////////////////////////////////////////

#include <FieldUtils/Interpolator.h>
#include <boost/geometry.hpp>

using namespace std;

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

namespace Nektar
{
namespace FieldUtils
{

/**
 * @brief Interpolate from one expansion to an other
 *
 * @param expInField    input field
 * @param expOutField   output field
 *
 *
 * In and output fields must have the same dimension and number of fields.
 * Weights are currently not stored for later use.
 * The interpolation is performed by evaluating the expInField at the quadrature
 * points of expOutField, so only eNoMethod is supported.
 * If both expansions use the same mesh, use LibUtilities/Foundations/Interp.h
 * instead.
 */
template <typename T>
void Interpolator<T>::Interpolate(const T expInField, T &expOutField,
                                  NekDouble def_value)
{
    ASSERTL0(expInField.size() == expOutField.size(),
             "number of fields does not match");
    ASSERTL0(expInField[0]->GetCoordim(0) <= GetDim(),
             "too many dimensions in inField");
    ASSERTL0(expOutField[0]->GetCoordim(0) <= GetDim(),
             "too many dimensions in outField");
    ASSERTL0(GetInterpMethod() == LibUtilities::eNoMethod,
             "only direct evaluation supported for this interpolation");

    int nFields      = max((int)expInField.size(), (int)expOutField.size());
    int nOutPts      = expOutField[0]->GetTotPoints();
    int outNumHomDir = 0;
    if (expOutField[0]->GetExpType() == MultiRegions::e3DH1D ||
        expOutField[0]->GetExpType() == MultiRegions::e2DH1D)
    {
        outNumHomDir = 1;
    }
    else if (expOutField[0]->GetExpType() == MultiRegions::e3DH2D)
    {
        outNumHomDir = 2;
    }
    int outDim = expOutField[0]->GetCoordim(0) + outNumHomDir;

    // create intermediate Ptsfield that wraps the expOutField
    Array<OneD, Array<OneD, NekDouble>> pts(outDim);
    for (int i = 0; i < outDim; ++i)
    {
        pts[i] = Array<OneD, NekDouble>(nOutPts);
    }
    if (outDim == 1)
    {
        expOutField[0]->GetCoords(pts[0]);
    }
    else if (outDim == 2)
    {
        expOutField[0]->GetCoords(pts[0], pts[1]);
    }
    else if (outDim == 3)
    {
        expOutField[0]->GetCoords(pts[0], pts[1], pts[2]);
    }

    LibUtilities::PtsFieldSharedPtr tmpPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(outDim, pts);
    for (int f = 0; f < expOutField.size(); ++f)
    {
        tmpPts->AddField(expOutField[f]->GetPhys(), "DefaultVar");
    }
    // interpolate m_ptsInField to this intermediate field
    Interpolate(expInField, tmpPts, def_value);
    // write the intermediate fields data into our expOutField
    for (int i = 0; i < nFields; i++)
    {
        int ptsi = outDim + i;
        for (int j = 0; j < nOutPts; ++j)
        {
            expOutField[i]->UpdatePhys()[j] = tmpPts->GetPointVal(ptsi, j);
        }
    }
}

/**
 * @brief Interpolate from an expansion to a pts field
 *
 * @param expInField    input field
 * @param ptsOutField   output field
 *
 * In and output fields must have the same dimension and number of fields.
 * Weights are currently not stored for later use.
 * The interpolation is performed by evaluating the expInField at the points
 * of ptsOutField, so only eNoMethod is supported.
 */
template <typename T>
void Interpolator<T>::Interpolate(const T expInField,
                                  LibUtilities::PtsFieldSharedPtr &ptsOutField,
                                  NekDouble def_value)
{
    ASSERTL0(expInField.size() == ptsOutField->GetNFields(),
             "number of fields does not match");
    ASSERTL0(expInField[0]->GetCoordim(0) <= GetDim(),
             "too many dimesions in inField");
    ASSERTL0(ptsOutField->GetDim() <= GetDim(),
             "too many dimesions in outField");
    ASSERTL0(ptsOutField->GetDim() >= expInField[0]->GetCoordim(0),
             "too few dimesions in outField");
    ASSERTL0(GetInterpMethod() == LibUtilities::eNoMethod,
             "only direct evaluation supported for this interpolation");
    ASSERTL0(expInField[0]->GetExpType() != MultiRegions::e3DH2D,
             "interpolation from 3DH2D expansion unsupported");

    m_ptsOutField = ptsOutField;

    int nInDim   = expInField[0]->GetCoordim(0);
    int nOutPts  = m_ptsOutField->GetNpoints();
    int lastProg = 0;

    int elmtid = -1;
    for (int i = 0; i < nOutPts; ++i)
    {
        Array<OneD, NekDouble> Lcoords(nInDim, 0.0);
        Array<OneD, NekDouble> coords(3, 0.0);
        for (int j = 0; j < m_ptsOutField->GetDim(); ++j)
        {
            coords[j] = m_ptsOutField->GetPointVal(j, i);
        }

        // Obtain Element and LocalCoordinate to interpolate.
        elmtid = expInField[0]->GetExpIndex(
            coords, Lcoords, NekConstants::kGeomFactorsTol, true, elmtid,
            NekConstants::kGeomFactorsTol * 1e3);

        // we use kGeomFactorsTol as tolerance, while StdPhysEvaluate has
        // kNekZeroTol hardcoded, so we need to limit Lcoords to not produce
        // a ton of warnings
        for (int j = 0; j < nInDim; ++j)
        {
            Lcoords[j] = std::max(Lcoords[j], -1.0);
            Lcoords[j] = std::min(Lcoords[j], 1.0);
        }

        if (elmtid >= 0)
        {
            int offset = expInField[0]->GetPhys_Offset(elmtid);

            for (int f = 0; f < expInField.size(); ++f)
            {
                NekDouble value;

                if (expInField[f]->GetExpType() == MultiRegions::e3DH1D ||
                    expInField[f]->GetExpType() == MultiRegions::e2DH1D)
                {
                    ASSERTL0(expInField[f]->GetWaveSpace(),
                             "interpolation from 3DH1D/2DH1D requires field in "
                             "wavespace");
                    NekDouble lHom = expInField[f]->GetHomoLen();
                    NekDouble BetaT =
                        2. * M_PI * fmod(coords[nInDim], lHom) / lHom;
                    int nPlanes =
                        expInField[f]->GetHomogeneousBasis()->GetZ().size();
                    NekDouble coeff = 0.;
                    Array<OneD, const unsigned int> planes =
                        expInField[f]->GetZIDs();
                    value = 0.;

                    for (size_t n = 0; n < planes.size(); ++n)
                    {
                        auto planeExp = expInField[f]->GetPlane(planes[n]);
                        coeff = planeExp->GetExp(elmtid)->StdPhysEvaluate(
                            Lcoords, planeExp->GetPhys() + offset);

                        if (planes[n] == 0)
                        {
                            value += coeff;
                        }
                        else if (planes[n] == 1)
                        {
                            value += cos(0.5 * nPlanes * BetaT) * coeff;
                        }
                        else if (planes[n] % 2 == 0)
                        {
                            NekDouble phase = (planes[n] >> 1) * BetaT;
                            value += cos(phase) * coeff;
                        }
                        else
                        {
                            NekDouble phase = (planes[n] >> 1) * BetaT;
                            value += -sin(phase) * coeff;
                        }
                    }
                }
                else
                {
                    value = expInField[f]->GetExp(elmtid)->StdPhysEvaluate(
                        Lcoords, expInField[f]->GetPhys() + offset);
                }

                if ((boost::math::isnan)(value))
                {
                    ASSERTL0(false, "new value is not a number");
                }
                else
                {
                    m_ptsOutField->SetPointVal(m_ptsOutField->GetDim() + f, i,
                                               value);
                }
            }
        }
        else
        {
            for (int f = 0; f < expInField.size(); ++f)
            {
                m_ptsOutField->SetPointVal(m_ptsOutField->GetDim() + f, i,
                                           def_value);
            }
        }

        int progress = int(100 * i / nOutPts);
        if (m_progressCallback && progress > lastProg)
        {
            m_progressCallback(i, nOutPts);
            lastProg = progress;
        }
    }
}

/**
 * @brief Interpolate from a pts field to an expansion
 *
 * @param ptsInField    input field
 * @param expOutField   output field
 *
 * In and output fields must have the same dimension and number of fields.
 */
template <typename T>
void Interpolator<T>::Interpolate(
    const LibUtilities::PtsFieldSharedPtr ptsInField, T &expOutField)
{
    ASSERTL0(expOutField.size() == ptsInField->GetNFields(),
             "number of fields does not match");
    ASSERTL0(ptsInField->GetDim() <= GetDim(), "too many dimesions in inField");
    ASSERTL0(expOutField[0]->GetCoordim(0) <= GetDim(),
             "too many dimesions in outField");

    m_ptsInField = ptsInField;

    int nFields = max((int)ptsInField->GetNFields(), (int)expOutField.size());
    int nOutPts = expOutField[0]->GetTotPoints();
    int outNumHomDir = 0;
    if (expOutField[0]->GetExpType() == MultiRegions::e3DH1D ||
        expOutField[0]->GetExpType() == MultiRegions::e2DH1D)
    {
        outNumHomDir = 1;
    }
    else if (expOutField[0]->GetExpType() == MultiRegions::e3DH2D)
    {
        outNumHomDir = 2;
    }
    int outDim = expOutField[0]->GetCoordim(0) + outNumHomDir;

    // create intermediate Ptsfield that wraps the expOutField
    Array<OneD, Array<OneD, NekDouble>> pts(outDim);
    for (int i = 0; i < outDim; ++i)
    {
        pts[i] = Array<OneD, NekDouble>(nOutPts);
    }
    if (outDim == 1)
    {
        expOutField[0]->GetCoords(pts[0]);
    }
    else if (outDim == 2)
    {
        expOutField[0]->GetCoords(pts[0], pts[1]);
    }
    else if (outDim == 3)
    {
        expOutField[0]->GetCoords(pts[0], pts[1], pts[2]);
    }

    LibUtilities::PtsFieldSharedPtr tmpPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(outDim, pts);

    for (int f = 0; f < expOutField.size(); ++f)
    {
        tmpPts->AddField(expOutField[f]->GetPhys(),
                         m_ptsInField->GetFieldName(f));
    }
    // interpolate m_ptsInField to this intermediate field
    LibUtilities::Interpolator::Interpolate(m_ptsInField, tmpPts);
    // write the intermediate fields data into our expOutField
    for (int i = 0; i < nFields; i++)
    {
        int ptsi = outDim + i;
        for (int j = 0; j < nOutPts; ++j)
        {
            expOutField[i]->UpdatePhys()[j] = tmpPts->GetPointVal(ptsi, j);
        }
    }
}

template <typename T>
void Interpolator<T>::InterpExp1ToExp2(const T exp1, T &exp2)
{

    // Interpolation from exp1 -> exp2 assuming that exp1 and exp2 are the
    // same explists, but at potentially different polynomial orders.
    if (exp1.size() != exp2.size())
    {
        NEKERROR(ErrorUtil::efatal, "not the same mesh")
    }

    for (int n = 0; n < exp1.size(); ++n)
    {
        // Interpolation from exp1 -> exp2 assuming that exp1 and exp2 are the
        // same explists, but at potentially different polynomial orders.
        if (exp1[n]->GetExpSize() != exp2[n]->GetExpSize())
        {
            NEKERROR(ErrorUtil::efatal, "not the same mesh")
        }

        // If same polynomial orders, simply copy solution
        if (exp1[n]->GetTotPoints() == exp2[n]->GetTotPoints())
        {
            Vmath::Vcopy(exp1[n]->GetTotPoints(), exp1[n]->GetPhys(), 1,
                         exp2[n]->UpdatePhys(), 1);
        }
        // If different polynomial orders, interpolate solution
        else
        {
            // Transform solution from physical to coefficient space
            exp1[n]->FwdTransLocalElmt(exp1[n]->GetPhys(),
                                       exp1[n]->UpdateCoeffs());

            for (int i = 0; i < exp1[n]->GetExpSize(); ++i)
            {
                // Get the elements
                LocalRegions::ExpansionSharedPtr elmt1 = exp1[n]->GetExp(i),
                                                 elmt2 = exp2[n]->GetExp(i);

                // Get the offset of elements in the storage arrays.
                int offset1 = exp1[n]->GetCoeff_Offset(i);
                int offset2 = exp2[n]->GetCoeff_Offset(i);

                // Get number of modes
                Array<OneD, LibUtilities::BasisSharedPtr> base1 =
                    elmt1->GetBase();
                std::vector<unsigned int> nummodes1(base1.size());
                std::vector<LibUtilities::BasisType> btype1(base1.size());
                for (int j = 0; j < nummodes1.size(); ++j)
                {
                    nummodes1[j] = base1[j]->GetNumModes();
                    btype1[j]    = base1[j]->GetBasisType();
                }

                // Extract data from exp1 -> exp2.
                elmt2->ExtractDataToCoeffs(
                    &exp1[n]->GetCoeffs()[offset1], nummodes1, 0,
                    &exp2[n]->UpdateCoeffs()[offset2], btype1);
            }

            // Transform solution back to physical space
            exp2[n]->BwdTrans(exp2[n]->GetCoeffs(), exp2[n]->UpdatePhys());
        }
    }
}

template <typename T>
void Interpolator<T>::Interpolate(
    const LibUtilities::PtsFieldSharedPtr ptsInField,
    LibUtilities::PtsFieldSharedPtr &ptsOutField)
{
    LibUtilities::Interpolator::Interpolate(ptsInField, ptsOutField);
}

template class Interpolator<std::vector<MultiRegions::ExpListSharedPtr>>;
template class Interpolator<Array<OneD, MultiRegions::ExpListSharedPtr>>;
} // namespace FieldUtils
} // namespace Nektar
