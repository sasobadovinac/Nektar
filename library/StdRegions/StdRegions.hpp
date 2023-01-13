///////////////////////////////////////////////////////////////////////////////
//
// File: StdRegions.hpp
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
// Description: Definition of enum lists and constants
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDREGIONS_HPP
#define STDREGIONS_HPP

#include <LibUtilities/BasicUtils/HashUtils.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <map>

namespace Nektar
{

/** \brief The namespace associated with the the StdRegions library
 * (\ref pageStdRegions "StdRegions introduction")
 */
namespace StdRegions
{
enum ElementType
{
    // eStdPointExp,
    eStdSegExp,
    eSegExp,
    eStdQuadExp,
    eStdTriExp,
    eStdNodalTriExp,
    eQuadExp,
    eTriExp,
    eNodalTriExp,
    eStdHexExp,
    eStdPrismExp,
    eStdPyrExp,
    eStdTetExp,
    eStdNodalTetExp,
    eHexExp,
    ePrismExp,
    ePyrExp,
    eTetExp,
    eNodalTetExp,
    SIZE_ElementType
};

const char *const ElementTypeMap[] = {
    //"StdPointExp",
    "StdSegExp", "SegExp",    "StdQuadExp",     "StdTriExp", "StdNodalTriExp",
    "QuadExp",   "TriExp",    "NodalTriExp",    "StdHexExp", "StdPrismExp",
    "StdPyrExp", "StdTetExp", "StdNodalTetExp", "HexExp",    "PrismExp",
    "PyrExp",    "TetExp",    "NodalTetExp",
};

/** @todo we need to tidy up matrix construction approach
 *  probably using a factory type approach
 */
enum MatrixType
{
    eNoMatrixType,
    eMass,
    eMassGJP,
    eInvMass,
    eLaplacian,
    eLaplacian00,
    eLaplacian01,
    eLaplacian02,
    eLaplacian10,
    eLaplacian11,
    eLaplacian12,
    eLaplacian20,
    eLaplacian21,
    eLaplacian22,
    eInvLaplacianWithUnityMean,
    eWeakDeriv0,
    eWeakDeriv1,
    eWeakDeriv2,
    eWeakDirectionalDeriv,
    eMassLevelCurvature,
    eLinearAdvectionReaction,
    eLinearAdvectionDiffusionReaction,
    eNBasisTrans,
    eInvNBasisTrans,
    eBwdTrans,
    eBwdMat,
    eIProductWRTBase,
    eIProductWRTDerivBase0,
    eIProductWRTDerivBase1,
    eIProductWRTDerivBase2,
    eDerivBase0,
    eDerivBase1,
    eDerivBase2,
    eHelmholtz,
    eHelmholtzGJP,
    eHybridDGHelmholtz,
    eInvHybridDGHelmholtz,
    eHybridDGHelmBndLam,
    eHybridDGLamToQ0,
    eHybridDGLamToQ1,
    eHybridDGLamToQ2,
    eHybridDGLamToU,
    eFwdTrans,
    ePreconR,
    ePreconRMass,
    ePreconLinearSpace,
    ePreconLinearSpaceMass,
    eInterpGauss,
    eGaussDG,
    ePhysInterpToEquiSpaced,
    eEquiSpacedToCoeffs,
    eNormDerivOnTrace,
    SIZE_MatrixType
};

const char *const MatrixTypeMap[] = {"NoMatrixType",
                                     "Mass",
                                     "Mass wiht Diagonal GJP",
                                     "InvMass",
                                     "Laplacian",
                                     "Laplacian00",
                                     "Laplacian01",
                                     "Laplacian02",
                                     "Laplacian10",
                                     "Laplacian11",
                                     "Laplacian12",
                                     "Laplacian20",
                                     "Laplacian21",
                                     "Laplacian22",
                                     "InvLaplacianWithUnityMean",
                                     "WeakDeriv0",
                                     "WeakDeriv1",
                                     "WeakDeriv2",
                                     "WeakDirectionalDeriv",
                                     "MassLevelCurvature",
                                     "LinearAdvectionReaction",
                                     "LinearAdvectionDiffusionReaction",
                                     "NBasisTrans",
                                     "InvNBasisTrans",
                                     "BwdTrans",
                                     "BwdMat",
                                     "IProductWRTBase",
                                     "IProductWRTDerivBase0",
                                     "IProductWRTDerivBase1",
                                     "IProductWRTDerivBase2",
                                     "DerivBase0",
                                     "DerivBase1",
                                     "DerivBase2",
                                     "Helmholtz",
                                     "Helmholtz with Diagonal GJP",
                                     "HybridDGHelmholz",
                                     "InvHybridDGHelmholtz",
                                     "HybridDGHelmBndLam",
                                     "HybridDGLamToQ0",
                                     "HybridDGLamToQ1",
                                     "HybridDGLamToQ2",
                                     "HybridDGLamToU",
                                     "FwdTrans",
                                     "PreconR",
                                     "PreconRMass",
                                     "PreconLinearSpace",
                                     "PreconLinearSpaceMass",
                                     "InterpGauss",
                                     "GaussDG",
                                     "PhysInterpToEquiSpaced",
                                     "EquiSpacedToCoeffs",
                                     "NormDerivOnTrace"};

enum VarCoeffType
{
    eVarCoeffMass,
    eVarCoeffLaplacian,
    eVarCoeffWeakDeriv,
    eVarCoeffD00,
    eVarCoeffD11,
    eVarCoeffD22,
    eVarCoeffD01,
    eVarCoeffD02,
    eVarCoeffD12,
    eVarCoeffVelX,
    eVarCoeffVelY,
    eVarCoeffVelZ,
    eVarCoeffMF1x,
    eVarCoeffMF1y,
    eVarCoeffMF1z,
    eVarCoeffMF1Div,
    eVarCoeffMF1Mag,
    eVarCoeffMF2x,
    eVarCoeffMF2y,
    eVarCoeffMF2z,
    eVarCoeffMF2Div,
    eVarCoeffMF2Mag,
    eVarCoeffMF3x,
    eVarCoeffMF3y,
    eVarCoeffMF3z,
    eVarCoeffMF3Div,
    eVarCoeffMF3Mag,
    eVarCoeffMF,
    eVarCoeffMFDiv,
    eVarCoeffGmat,
    eVarCoeffGJPNormVel,
    SIZE_VarCoeffType
};

const char *const VarCoeffTypeMap[] = {
    "VarCoeffMass",      "VarCoeffLaplacian", "VarCoeffWeakDeriv",
    "VarCoeffD00",       "VarCoeffD11",       "VarCoeffD22",
    "VarCoeffD01",       "VarCoeffD02",       "VarCoeffD12",
    "VarCoeffVelX",      "VarCoeffVelY",      "VarCoeffVelZ",
    "VarCoeffMF1x",      "VarCoeffMF1y",      "VarCoeffMF1z",
    "VarCoeffMF1Div",    "VarCoeffMF1Mag",    "VarCoeffMF2x",
    "VarCoeffMF2y",      "VarCoeffMF2z",      "VarCoeffMF2Div",
    "VarCoeffMF2Mag",    "VarCoeffMF3x",      "VarCoeffMF3y",
    "VarCoeffMF3z",      "VarCoeffMF3Div",    "VarCoeffMF3Mag",
    "VarCoeffMF",        "VarCoeffMFDiv",     "VarCoeffGmat",
    "VarCoeffGJPNormVel"};

/**
 * @brief Representation of a variable coefficient.
 *
 * Variable coefficients are entries stored inside a #VarCoeffMap which
 * are defined at each quadrature/solution point within an element. This
 * class wraps that concept, storing the values in #m_coeffs, but also
 * stores alongside this a hash of the data in #m_hash. This is then
 * used within MultiRegions::GlobalMatrixKey to efficiently distinguish
 * between matrix keys that have variable coefficients defined, but
 * whose entries are different.
 *
 * For that reason the entries here are deliberately protected by const
 * references; i.e. the entries inside of #m_coeffs should not be
 * modified in-place, but a new array copied in so that the hash can be
 * recalculated.
 */
struct VarCoeffEntry
{
    /// Default constructor.
    VarCoeffEntry() = default;

    /**
     * @brief Copy an array of values into this entry.
     *
     * Upon copy into #m_coeffs, compute the hash of this entry using
     * #ComputeHash.
     *
     * @param input  Variable coefficients to be defined at each
     *               solution point.
     */
    VarCoeffEntry(const Array<OneD, const NekDouble> &input) : m_coeffs(input)
    {
        ComputeHash();
    }

    /**
     * @brief Access an entry @p idx within #m_coeffs.
     *
     * @param idx  Index of the entry to access.
     */
    const NekDouble &operator[](std::size_t idx) const
    {
        return m_coeffs[idx];
    }

    /**
     * @brief Assignment operator given an array @p rhs.
     *
     * Upon copy into #m_coeffs, compute the hash of this entry using
     * #ComputeHash.
     *
     * @param rhs    Variable coefficients to be defined at each
     *               solution point.
     */
    void operator=(const Array<OneD, const NekDouble> &rhs)
    {
        m_coeffs = rhs;
        ComputeHash();
    }

    /**
     * @brief Returns a const reference to the coefficients.
     */
    const Array<OneD, const NekDouble> &GetValue() const
    {
        return m_coeffs;
    }

    /**
     * @brief Returns the hash of this entry.
     */
    std::size_t GetHash() const
    {
        return m_hash;
    }

protected:
    /**
     * @brief Computes the hash of this entry using #hash_range.
     */
    void ComputeHash()
    {
        if (m_coeffs.size() == 0)
        {
            m_hash = 0;
            return;
        }

        m_hash = hash_range(m_coeffs.begin(), m_coeffs.end());
    }

    /// Hash of the entries inside #m_coeffs.
    std::size_t m_hash = 0;

    /// Storage for the variable coefficient entries.
    Array<OneD, NekDouble> m_coeffs;
};

typedef std::map<StdRegions::VarCoeffType, VarCoeffEntry> VarCoeffMap;
static VarCoeffMap NullVarCoeffMap;

inline VarCoeffMap RestrictCoeffMap(const VarCoeffMap &m, size_t offset,
                                    size_t cnt)
{
    VarCoeffMap ret;

    for (auto &x : m)
    {
        ret[x.first] =
            Array<OneD, NekDouble>(cnt, x.second.GetValue() + offset);
    }

    return ret;
}

enum ConstFactorType
{
    eFactorLambda,
    eFactorCoeffD00,
    eFactorCoeffD11,
    eFactorCoeffD22,
    eFactorCoeffD01,
    eFactorCoeffD02,
    eFactorCoeffD12,
    eFactorTau,
    eFactorTime,
    eFactorSVVCutoffRatio,
    eFactorSVVDiffCoeff,
    eFactorSVVPowerKerDiffCoeff,
    eFactorSVVDGKerDiffCoeff,
    eFactorGaussVertex,
    eFactorGaussEdge,
    eFactorGJP,
    eFactorConst,
    SIZE_ConstFactorType
};

const char *const ConstFactorTypeMap[] = {"FactorLambda",
                                          "FactorCoeffD00",
                                          "FactorCoeffD11",
                                          "FactorCoeffD22",
                                          "FactorCoeffD01",
                                          "FactorCoeffD02",
                                          "FactorCoeffD12",
                                          "FactorTau",
                                          "FactorTime",
                                          "FactorSVVCutoffRatio",
                                          "FactorSVVDiffCoeff",
                                          "FactorSVVPowerKerDiffCoeff",
                                          "FactorSVVDGKerDiffCoeff",
                                          "FactorGaussVertex",
                                          "FactorGaussEdge",
                                          "FactorGJP",
                                          "FactorConstant"};
typedef std::map<ConstFactorType, NekDouble> ConstFactorMap;
static ConstFactorMap NullConstFactorMap;

// FactorMap
typedef ConstFactorMap FactorMap;
static FactorMap NullFactorMap;

enum Orientation
{
    eNoOrientation,
    eFwd,
    eBwd,
    eForwards,
    eBackwards,
    eDir1FwdDir1_Dir2FwdDir2, // These flags are interpreted as
    eDir1FwdDir1_Dir2BwdDir2, // taking the second direction to the
    eDir1BwdDir1_Dir2FwdDir2, // first direction. So Dir1FwdDir2 takes
    eDir1BwdDir1_Dir2BwdDir2, // direction 2 and makes it backward
    eDir1FwdDir2_Dir2FwdDir1, // to direction 1 in the mapped face.
    eDir1FwdDir2_Dir2BwdDir1, // Note be careful not to flip this
    eDir1BwdDir2_Dir2FwdDir1, // convention especially when using
    eDir1BwdDir2_Dir2BwdDir1, // transposed mappings.
    SIZE_Orientation
};

const char *const OrientationMap[] = {"NoOrientation",
                                      "Fwd",
                                      "Bwd",
                                      "Forwards",
                                      "Backwards",
                                      "Dir1FwdDir1_Dir2FwdDir2",
                                      "Dir1FwdDir1_Dir2BwdDir2",
                                      "Dir1BwdDir1_Dir2FwdDir2",
                                      "Dir1BwdDir1_Dir2BwdDir2",
                                      "Dir1FwdDir2_Dir2FwdDir1",
                                      "Dir1FwdDir2_Dir2BwdDir1",
                                      "Dir1BwdDir2_Dir2FwdDir1",
                                      "Dir1BwdDir2_Dir2BwdDir1"};

// Defines a "fast find"
// Assumes that first/last define the beginning/ending of
// a continuous range of classes, and that start is
// an iterator between first and last

template <class InputIterator, class EqualityComparable>
InputIterator find(InputIterator first, InputIterator last,
                   InputIterator startingpoint, const EqualityComparable &value)
{
    InputIterator val;

    if (startingpoint == first)
    {
        val = find(first, last, value);
    }
    else
    {
        val = find(startingpoint, last, value);
        if (val == last)
        {
            val = find(first, startingpoint, value);
            if (val == startingpoint)
            {
                val = last;
            }
        }
    }
    return val;
}

// Optimized Kernel Entries
const int kSVVDGFiltermodesmin = 3;
const int kSVVDGFiltermodesmax = 11;
// Optimized Kernel Entries for p = 2 - 10
const NekDouble kSVVDGFilter[9][11] = {
    {0, 0.36212, 1, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0.70546, 0.078836, 1, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0.49411, 0.072394, 1, 0, 0, 0, 0, 0, 0},
    {0, 0, 0.000073566, 0.40506, 0.094122, 1, 0, 0, 0, 0, 0},
    {0, 0, 0, 0.0001422, 0.36863, 0.11815, 1, 0, 0, 0, 0},
    {0, 0, 0, 0, 0.00019497, 0.41397, 0.16927, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0.0009762, 0.12747, 0.13763, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0.0023592, 0.23683, 0.17196, 1, 0},
    {0, 0, 0, 0, 0, 0, 0, 0.0026055, 0.28682, 0.22473, 1}};

} // namespace StdRegions
} // namespace Nektar

#endif // STDREGIONS_H
