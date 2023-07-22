///////////////////////////////////////////////////////////////////////////////
//
// File: StdPointExp.cpp
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
// Description: Definition of a Point expansion
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <StdRegions/StdPointExp.h>

using namespace std;

namespace Nektar
{
namespace StdRegions
{

StdPointExp::StdPointExp()
{
}

StdPointExp::StdPointExp(const LibUtilities::BasisKey &Ba)
    : StdExpansion(Ba.GetNumModes(), 1, Ba),
      StdExpansion0D(Ba.GetNumModes(), Ba)
{
}

/** \brief Copy Constructor */

StdPointExp::StdPointExp(const StdPointExp &T)
    : StdExpansion(T), StdExpansion0D(T)
{
}

StdPointExp::~StdPointExp()
{
}

LibUtilities::ShapeType StdPointExp::v_DetShapeType() const
{
    return LibUtilities::ePoint;
}

void StdPointExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                              Array<OneD, NekDouble> &coords_1,
                              Array<OneD, NekDouble> &coords_2)
{
    boost::ignore_unused(coords_1, coords_2);
    Blas::Dcopy(GetNumPoints(0), (m_base[0]->GetZ()).get(), 1, &coords_0[0], 1);
}

void StdPointExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                             Array<OneD, NekDouble> &outarray)
{
    int nquad = m_base[0]->GetNumPoints();

    if (m_base[0]->Collocation())
    {
        Vmath::Vcopy(nquad, inarray, 1, outarray, 1);
    }
    else
    {
        Blas::Dgemv('N', nquad, m_base[0]->GetNumModes(), 1.0,
                    (m_base[0]->GetBdata()).get(), nquad, &inarray[0], 1, 0.0,
                    &outarray[0], 1);
    }
}

void StdPointExp::v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                             Array<OneD, NekDouble> &outarray)
{
    if (m_base[0]->Collocation())
    {
        Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
    }
    else
    {
        v_IProductWRTBase(inarray, outarray);

        // get Mass matrix inverse
        StdMatrixKey masskey(eInvMass, v_DetShapeType(), *this);
        DNekMatSharedPtr matsys = GetStdMatrix(masskey);

        NekVector<NekDouble> in(m_ncoeffs, outarray, eCopy);
        NekVector<NekDouble> out(m_ncoeffs, outarray, eWrapper);

        out = (*matsys) * in;
    }
}

void StdPointExp::v_BwdTrans_SumFac(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
{
    v_BwdTrans(inarray, outarray);
}

// Inner product
void StdPointExp::v_IProductWRTBase(const Array<OneD, const NekDouble> &base,
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray,
                                    int coll_check)
{
    int nquad = m_base[0]->GetNumPoints();
    Array<OneD, NekDouble> tmp(nquad);
    Array<OneD, const NekDouble> z = m_base[0]->GetZ();
    Array<OneD, const NekDouble> w = m_base[0]->GetW();

    Vmath::Vmul(nquad, inarray, 1, w, 1, tmp, 1);

    if (coll_check && m_base[0]->Collocation())
    {
        Vmath::Vcopy(nquad, tmp, 1, outarray, 1);
    }
    else
    {
        Blas::Dgemv('T', nquad, m_ncoeffs, 1.0, base.get(), nquad, &tmp[0], 1,
                    0.0, outarray.get(), 1);
    }
}

void StdPointExp::v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
{
    v_IProductWRTBase(m_base[0]->GetBdata(), inarray, outarray, 1);
}

void StdPointExp::v_IProductWRTDerivBase(
    const int dir, const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray)
{
    boost::ignore_unused(dir);
    ASSERTL1(dir >= 0 && dir < 1, "input dir is out of range");
    v_IProductWRTBase(m_base[0]->GetDbdata(), inarray, outarray, 1);
}

void StdPointExp::v_IProductWRTBase_SumFac(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, bool multiplybyweights)
{
    boost::ignore_unused(multiplybyweights);
    v_IProductWRTBase(m_base[0]->GetBdata(), inarray, outarray, 1);
}

DNekMatSharedPtr StdPointExp::v_GenMatrix(const StdMatrixKey &mkey)
{
    DNekMatSharedPtr Mat;
    MatrixType mattype;

    switch (mattype = mkey.GetMatrixType())
    {
        case eFwdTrans:
        {
            Mat =
                MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs, m_ncoeffs);
            StdMatrixKey iprodkey(eIProductWRTBase, v_DetShapeType(), *this);
            DNekMat &Iprod = *GetStdMatrix(iprodkey);
            StdMatrixKey imasskey(eInvMass, v_DetShapeType(), *this);
            DNekMat &Imass = *GetStdMatrix(imasskey);

            (*Mat) = Imass * Iprod;
        }
        break;
        default:
        {
            Mat = StdExpansion::CreateGeneralMatrix(mkey);
        }
        break;
    }

    return Mat;
}

DNekMatSharedPtr StdPointExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
{
    return v_GenMatrix(mkey);
}

} // namespace StdRegions
} // namespace Nektar
