///////////////////////////////////////////////////////////////////////////////
//
// File: ExpList2DHomogeneous1D.h
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
// Description: A 2D field which is homogeneous in 1 direction and so
// uses much of the functionality from a ExpList1D and its daughters
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPLIST2DHOMO1D_H
#define EXPLIST2DHOMO1D_H

#include <MultiRegions/ExpListHomogeneous1D.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <vector>

namespace Nektar
{
namespace MultiRegions
{

// Forward declaration for typedefs
class ExpList2DHomogeneous1D;

/// Shared pointer to an ExpList2DHomogeneous1D object.
typedef std::shared_ptr<ExpList2DHomogeneous1D> ExpList2DHomogeneous1DSharedPtr;
/// Vector of pointers to ExpList2DHomogeneous1D objects.
typedef std::vector<ExpList2DHomogeneous1DSharedPtr>
    ExpList2DHomogeneous1DVector;

/// Abstraction of a two-dimensional multi-elemental expansion which
/// is merely a collection of local expansions.
class ExpList2DHomogeneous1D : public ExpListHomogeneous1D
{
public:
    /// Default constructor.
    MULTI_REGIONS_EXPORT ExpList2DHomogeneous1D();

    /// Sets up a list of local expansions based on an input mesh.
    MULTI_REGIONS_EXPORT ExpList2DHomogeneous1D(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::BasisKey &HomoBasis, const NekDouble lz,
        const bool useFFT, const bool dealiasing,
        const SpatialDomains::MeshGraphSharedPtr &graph1D,
        const Collections::ImplementationType ImpType =
            Collections::eNoImpType);

    MULTI_REGIONS_EXPORT ExpList2DHomogeneous1D(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::BasisKey &HomoBasis, const NekDouble lhom,
        const bool useFFT, const bool dealiasing,
        const Array<OneD, ExpListSharedPtr> &planes,
        const LibUtilities::CommSharedPtr comm = LibUtilities::CommSharedPtr());

    /// Copy constructor.
    MULTI_REGIONS_EXPORT ExpList2DHomogeneous1D(
        const ExpList2DHomogeneous1D &In);

    /// Destructor.
    MULTI_REGIONS_EXPORT virtual ~ExpList2DHomogeneous1D();

protected:
    /// Definition of the total number of degrees of freedom and
    /// quadrature points. Sets up the storage for \a m_coeff and \a
    /// m_phys.
    void SetCoeffPhys(void);

    //  virtual functions

    virtual void v_GetCoords(Array<OneD, NekDouble> &coord_0,
                             Array<OneD, NekDouble> &coord_1,
                             Array<OneD, NekDouble> &coord_2) override;

    virtual void v_GetCoords(const int eid, Array<OneD, NekDouble> &xc0,
                             Array<OneD, NekDouble> &xc1,
                             Array<OneD, NekDouble> &xc2) override;

    virtual void v_WriteTecplotZone(std::ostream &outfile,
                                    int expansion) override;

    virtual void v_WriteVtkPieceHeader(std::ostream &outfile, int expansion,
                                       int istrip) override;

    virtual void v_GetNormals(
        Array<OneD, Array<OneD, NekDouble>> &normals) override;

private:
};

} // namespace MultiRegions
} // namespace Nektar

#endif // EXPLIST3DHOMO1D_H
