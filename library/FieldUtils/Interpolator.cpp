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

#include <boost/geometry.hpp>
#include <FieldUtils/Interpolator.h>

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
template <typename ArrayExpListSharedPtr>
void Interpolator<ArrayExpListSharedPtr>::Interpolate(
    const ArrayExpListSharedPtr expInField,
    ArrayExpListSharedPtr &expOutField,
    NekDouble def_value)
{
    ASSERTL0(expInField.size() == expOutField.size(),
             "number of fields does not match");
    ASSERTL0(expInField[0]->GetCoordim(0) <= GetDim(),
             "too many dimesions in inField");
    ASSERTL0(expOutField[0]->GetCoordim(0) <= GetDim(),
             "too many dimesions in outField");
    ASSERTL0(GetInterpMethod() == LibUtilities::eNoMethod,
             "only direct evaluation supported for this interpolation");

    m_expInField  = expInField;
    m_expOutField = expOutField;

    int nFields = max((int)expInField.size(), (int)m_expOutField.size());
    int nOutPts = m_expOutField[0]->GetTotPoints();
    int outNumHomDir = 0;
    if (m_expOutField[0]->GetExpType() == MultiRegions::e3DH1D ||
        m_expOutField[0]->GetExpType() == MultiRegions::e2DH1D)
    {
        outNumHomDir = 1;
    }
    else if(m_expOutField[0]->GetExpType() == MultiRegions::e3DH2D)
    {
        outNumHomDir = 2;
    }
    int outDim  = m_expOutField[0]->GetCoordim(0) + outNumHomDir;

    // create intermediate Ptsfield that wraps the expOutField
    Array<OneD, Array<OneD, NekDouble> > pts(outDim);
    for (int i = 0; i < outDim; ++i)
    {
        pts[i] = Array<OneD, NekDouble>(nOutPts);
    }
    if (outDim == 1)
    {
        m_expOutField[0]->GetCoords(pts[0]);
    }
    else if (outDim == 2)
    {
        m_expOutField[0]->GetCoords(pts[0], pts[1]);
    }
    else if (outDim == 3)
    {
        m_expOutField[0]->GetCoords(pts[0], pts[1], pts[2]);
    }

    LibUtilities::PtsFieldSharedPtr tmpPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(outDim, pts);
    for (int f = 0; f < expOutField.size(); ++f)
    {
        tmpPts->AddField(m_expOutField[f]->GetPhys(), "DefaultVar");
    }

    // interpolate m_ptsInField to this intermediate field
    Interpolate(m_expInField, tmpPts, def_value);

    // write the intermediate fields data into our expOutField
    for (int i = 0; i < nFields; i++)
    {
        int ptsi = outDim + i;
        for (int j = 0; j < nOutPts; ++j)
        {
            m_expOutField[i]->UpdatePhys()[j] = tmpPts->GetPointVal(ptsi, j);
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
template <typename ArrayExpListSharedPtr>
void Interpolator<ArrayExpListSharedPtr>::Interpolate(
    const ArrayExpListSharedPtr expInField,
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

    m_expInField  = expInField;
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
        elmtid = m_expInField[0]->GetExpIndex(
            coords, Lcoords,
            NekConstants::kGeomFactorsTol, true, elmtid,
            NekConstants::kGeomFactorsTol*1e3);

        // we use kGeomFactorsTol as tolerance, while StdPhysEvaluate has
        // kNekZeroTol hardcoded, so we need to limit Lcoords to not produce
        // a ton of warnings
        for(int j = 0; j < nInDim; ++j)
        {
            Lcoords[j] = std::max(Lcoords[j], -1.0);
            Lcoords[j] = std::min(Lcoords[j], 1.0);
        }

        if (elmtid >= 0)
        {
            int offset = m_expInField[0]->GetPhys_Offset(elmtid);

            for (int f = 0; f < m_expInField.size(); ++f)
            {
                NekDouble value;
                if (m_expInField[f]->GetExpType() == MultiRegions::e3DH1D ||
                    m_expInField[f]->GetExpType() == MultiRegions::e2DH1D)
                {
                    ASSERTL0(m_expInField[f]->GetWaveSpace(),
                        "interpolation from 3DH1D/2DH1D requires field in wavespace");
                    NekDouble lHom = m_expInField[f]->GetHomoLen();
                    NekDouble BetaT = 2.*M_PI*fmod (coords[nInDim], lHom) / lHom;
                    int nPlanes = m_expInField[f]->GetHomogeneousBasis()->GetZ().size();
                    NekDouble coeff = 0.;
                    Array<OneD, const unsigned int> planes = m_expInField[f]->GetZIDs();
                    value = 0.;
                    for ( size_t n = 0; n < planes.size(); ++n)
                    {
                        auto planeExp = m_expInField[f]->GetPlane(planes[n]);
                        coeff = planeExp->GetExp(elmtid)->StdPhysEvaluate(
                            Lcoords, planeExp->GetPhys() + offset);
                        if (planes[n] == 0)
                        {
                            value += coeff;
                        }
                        else if (planes[n] == 1)
                        {
                            value += cos(0.5*nPlanes*BetaT)*coeff;
                        }
                        else if (planes[n]%2 == 0)
                        {
                            NekDouble phase = (planes[n]>>1) * BetaT;
                            value += cos(phase)*coeff;
                        }
                        else
                        {
                            NekDouble phase = (planes[n]>>1) * BetaT;
                            value += - sin(phase)*coeff;
                        }
                    }
                }
                else
                {
                    value = m_expInField[f]->GetExp(elmtid)->StdPhysEvaluate(
                        Lcoords, m_expInField[f]->GetPhys() + offset);
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
            for (int f = 0; f < m_expInField.size(); ++f)
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
template <typename ArrayExpListSharedPtr>
void Interpolator<ArrayExpListSharedPtr>::Interpolate(
    const LibUtilities::PtsFieldSharedPtr ptsInField,
    ArrayExpListSharedPtr &expOutField)
{
    ASSERTL0(expOutField.size() == ptsInField->GetNFields(),
             "number of fields does not match");
    ASSERTL0(ptsInField->GetDim() <= GetDim(), "too many dimesions in inField");
    ASSERTL0(expOutField[0]->GetCoordim(0) <= GetDim(),
             "too many dimesions in outField");

    m_ptsInField  = ptsInField;
    m_expOutField = expOutField;

    int nFields = max((int)ptsInField->GetNFields(), (int)m_expOutField.size());
    int nOutPts = m_expOutField[0]->GetTotPoints();
    int outNumHomDir = 0;
    if (m_expOutField[0]->GetExpType() == MultiRegions::e3DH1D ||
        m_expOutField[0]->GetExpType() == MultiRegions::e2DH1D)
    {
        outNumHomDir = 1;
    }
    else if(m_expOutField[0]->GetExpType() == MultiRegions::e3DH2D)
    {
        outNumHomDir = 2;
    }
    int outDim  = m_expOutField[0]->GetCoordim(0) + outNumHomDir;

    // create intermediate Ptsfield that wraps the expOutField
    Array<OneD, Array<OneD, NekDouble> > pts(outDim);
    for (int i = 0; i < outDim; ++i)
    {
        pts[i] = Array<OneD, NekDouble>(nOutPts);
    }
    if (outDim == 1)
    {
        m_expOutField[0]->GetCoords(pts[0]);
    }
    else if (outDim == 2)
    {
        m_expOutField[0]->GetCoords(pts[0], pts[1]);
    }
    else if (outDim == 3)
    {
        m_expOutField[0]->GetCoords(pts[0], pts[1], pts[2]);
    }

    LibUtilities::PtsFieldSharedPtr tmpPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(outDim, pts);
    for (int f = 0; f < expOutField.size(); ++f)
    {
        tmpPts->AddField(m_expOutField[f]->GetPhys(),
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
            m_expOutField[i]->UpdatePhys()[j] = tmpPts->GetPointVal(ptsi, j);
        }
    }
}

template <typename ArrayExpListSharedPtr>
void Interpolator<ArrayExpListSharedPtr>::Interpolate(
    const LibUtilities::PtsFieldSharedPtr ptsInField,
    LibUtilities::PtsFieldSharedPtr &ptsOutField)
{
    LibUtilities::Interpolator::Interpolate(ptsInField, ptsOutField);
}


}
}
