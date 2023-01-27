///////////////////////////////////////////////////////////////////////////////
//
// File: ExpList2DHomogeneous2D.h
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
// Description: A 1D field which is homogeneous in 2 directions and so
// uses much of the functionality from a ExpList2D and its daughters
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPLIST1DHOMO2D_H
#define EXPLIST1DHOMO2D_H

#include <MultiRegions/ExpListHomogeneous2D.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <vector>

namespace Nektar
{
namespace MultiRegions
{

// Forward declaration for typedefs
class ExpList2DHomogeneous2D;

/// Shared pointer to an ExpList2DHomogeneous2D object.
typedef std::shared_ptr<ExpList2DHomogeneous2D> ExpList2DHomogeneous2DSharedPtr;
/// Vector of pointers to ExpList2DHomogeneous2D objects.
typedef std::vector<ExpList2DHomogeneous2DSharedPtr>
    ExpList2DHomogeneous2DVector;

/// Abstraction of a one-dimensional multi-elemental expansion which
/// is merely a collection of local expansions.
class ExpList2DHomogeneous2D : public ExpListHomogeneous2D
{
public:
    /// Default constructor.
    MULTI_REGIONS_EXPORT ExpList2DHomogeneous2D();

    MULTI_REGIONS_EXPORT ExpList2DHomogeneous2D(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::BasisKey &HomoBasis_y,
        const LibUtilities::BasisKey &HomoBasis_z, const NekDouble lhom_y,
        const NekDouble lhom_z, const bool useFFT, const bool dealiasing,
        const Array<OneD, ExpListSharedPtr> &points);

    /// Copy constructor.
    MULTI_REGIONS_EXPORT ExpList2DHomogeneous2D(
        const ExpList2DHomogeneous2D &In);

    /// Destructor.
    MULTI_REGIONS_EXPORT virtual ~ExpList2DHomogeneous2D();

    // MULTI_REGIONS_EXPORT void HomoFwdTrans2D(const Array<OneD, const
    // NekDouble> &inarray, Array<OneD, NekDouble> &outarray);

protected:
    /// Definition of the total number of degrees of freedom and
    /// quadrature points. Sets up the storage for \a m_coeff and \a
    ///  m_phys.
    void SetCoeffPhys(void);

    //  virtual functions

    virtual void v_GetCoords(Array<OneD, NekDouble> &coord_0,
                             Array<OneD, NekDouble> &coord_1,
                             Array<OneD, NekDouble> &coord_2) override;

    virtual void v_GetCoords(const int eid, Array<OneD, NekDouble> &xc0,
                             Array<OneD, NekDouble> &xc1,
                             Array<OneD, NekDouble> &xc2) override;

    virtual void v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                            Array<OneD, NekDouble> &outarray) override;

    virtual void v_FwdTransLocalElmt(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    // This is same as fwdtrans for this expansion
    virtual void v_FwdTransBndConstrained(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray) override;

    virtual void v_WriteTecplotZone(std::ostream &outfile,
                                    int expansion) override;

    virtual void v_WriteVtkPieceHeader(std::ostream &outfile, int expansion,
                                       int istrip) override;

private:
};

} // namespace MultiRegions
} // namespace Nektar

#endif // EXPLIST3DHOMO1D_H
