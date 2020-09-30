///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion1D.h
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
// Description: Daughter of StdExpansion. This class contains routine
// which are common to 1d expansion. Typically this inolves physical
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDEXP1D_H
#define STDEXP1D_H

#include <StdRegions/StdExpansion.h>
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar
{
namespace StdRegions
{

class StdExpansion1D : virtual public StdExpansion
{

public:
    STD_REGIONS_EXPORT StdExpansion1D();
    STD_REGIONS_EXPORT StdExpansion1D(int numcoeffs,
                                      const LibUtilities::BasisKey &Ba);
    STD_REGIONS_EXPORT StdExpansion1D(const StdExpansion1D &T);
    STD_REGIONS_EXPORT virtual ~StdExpansion1D();

    /** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
     *  physical quadrature points given by \a inarray and return in
     *  \a outarray.
     *
     *  \param inarray array of a function evaluated at the quadrature
     *  points
     *  \param outarray the resulting array of the derivative \f$
     *  du/d_{\xi_1}|_{\xi_{1i}} \f$ will be stored in the array
     *  \a outarray as output of the function
     */
    STD_REGIONS_EXPORT void PhysTensorDeriv(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray);

     // find derivative of u (inarray) at all coords points
     STD_REGIONS_EXPORT inline NekDouble BaryTensorDeriv(
        const Array<OneD, NekDouble> &coord,
        const Array<OneD, const NekDouble> &inarray, NekDouble &out_d0)
    {
        return StdExpansion::BaryEvaluate<0, true>(coord[0], &inarray[0],
                                                   out_d0);
    }

    STD_REGIONS_EXPORT virtual Array<OneD, NekDouble> v_PhysEvaluateBasis(
        const Array<OneD, const Array<OneD, NekDouble>> coords,
        const Array<OneD, Array<OneD, NekDouble>> storage, int mode) final;

protected:
    STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
        const Array<OneD, const NekDouble> &coords,
        const Array<OneD, const NekDouble> &physvals) override;

    STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluateOld(
        const Array<OneD, const NekDouble> &coords,
        const Array<OneD, const NekDouble> &physvals) override;

    STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
        const Array<OneD, NekDouble> coord,
        const Array<OneD, const NekDouble> &inarray, NekDouble &out_d0,
        NekDouble &out_d1, NekDouble &out_d2) override;

    STD_REGIONS_EXPORT virtual Array<OneD, Array<OneD, NekDouble>>
    v_GetPhysEvaluateStorage() final;

    STD_REGIONS_EXPORT virtual Array<OneD, NekDouble> v_PhysEvaluateBasis(
        const Array<OneD, const Array<OneD, NekDouble>> coords,
        Array<OneD, Array<OneD, NekDouble>> storage,
        Array<OneD, NekDouble> &out_d0,
        Array<OneD, NekDouble> &out_d1 = NullNekDouble1DArray,
        Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray) final;

private:
    // Virtual Functions ----------------------------------------

    virtual int v_GetCoordim(void) override
    {
        return 1;
    }

    virtual int v_GetShapeDimension() const override
    {
        return 1;
    }
};

typedef std::shared_ptr<StdExpansion1D> StdExpansion1DSharedPtr;

} // namespace StdRegions
} // namespace Nektar

#endif // STDEXP1D_H
