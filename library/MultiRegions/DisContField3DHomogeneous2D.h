///////////////////////////////////////////////////////////////////////////////
//
// File: DisContField3DHomogeneous2D.h
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
// Description: Field definition in three-dimensions for a discontinuous
// LDG-H expansion with 2 homogeneous directions in 2D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD3DHOMO2D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD3DHOMO2D_H

#include <MultiRegions/ExpList3DHomogeneous2D.h>
#include <MultiRegions/MultiRegionsDeclspec.h>

namespace Nektar
{
namespace MultiRegions
{
class DisContField3DHomogeneous2D : public ExpList3DHomogeneous2D
{
public:
    MULTI_REGIONS_EXPORT DisContField3DHomogeneous2D();

    MULTI_REGIONS_EXPORT DisContField3DHomogeneous2D(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::BasisKey &HomoBasis_y,
        const LibUtilities::BasisKey &HomoBasis_z, const NekDouble lhom_y,
        const NekDouble lhom_z, const bool useFFT, const bool dealiasing,
        const Collections::ImplementationType ImpType =
            Collections::eNoImpType);

    MULTI_REGIONS_EXPORT DisContField3DHomogeneous2D(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::BasisKey &HomoBasis_y,
        const LibUtilities::BasisKey &HomoBasis_z, const NekDouble lhom_y,
        const NekDouble lhom_z, const bool useFFT, const bool dealiasing,
        const SpatialDomains::MeshGraphSharedPtr &graph1D,
        const std::string &variable,
        const Collections::ImplementationType ImpType =
            Collections::eNoImpType);

    /// Copy constructor.
    MULTI_REGIONS_EXPORT DisContField3DHomogeneous2D(
        const DisContField3DHomogeneous2D &In,
        const bool DeclareLinesSetCoeffPhys = true);

    /// Destructor.
    MULTI_REGIONS_EXPORT virtual ~DisContField3DHomogeneous2D();

    MULTI_REGIONS_EXPORT void SetupBoundaryConditions(
        const LibUtilities::BasisKey &HomoBasis_y,
        const LibUtilities::BasisKey &HomoBasis_z, const NekDouble lhom_y,
        const NekDouble lhom_z, SpatialDomains::BoundaryConditions &bcs,
        const std::string variable);

    /// Storage space for the boundary to element and boundary to trace map.
    /// This member variable is really allocated just in case a boundary
    /// expansion recasting is required at the solver level. Otherwise is the 2
    /// vectors are not filled up. If is needed all the funcitons whihc require
    /// to use this map do not have to recalculate it anymore.
    Array<OneD, int> m_BCtoElmMap;
    Array<OneD, int> m_BCtoEdgMap;

protected:
    Array<OneD, MultiRegions::ExpListSharedPtr> m_bndCondExpansions;

    Array<OneD, NekDouble> m_bndCondBndWeight;

    Array<OneD, SpatialDomains::BoundaryConditionShPtr> m_bndConditions;

    virtual GlobalLinSysKey v_HelmSolve(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const StdRegions::ConstFactorMap &factors,
        const StdRegions::VarCoeffMap &varcoeff,
        const MultiRegions::VarFactorsMap &varfactors,
        const Array<OneD, const NekDouble> &dirForcing,
        const bool PhysSpaceForcing) override;

    virtual void v_GetBndElmtExpansion(
        int i, std::shared_ptr<ExpList> &result,
        const bool DeclareCoeffPhysArrays) override;

    virtual void v_GetBoundaryToElmtMap(Array<OneD, int> &ElmtID,
                                        Array<OneD, int> &EdgeID) override;

    /// @todo Fix Robin BCs for homogeneous case
    virtual std::map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo() override;

    virtual void v_EvaluateBoundaryConditions(
        const NekDouble time = 0.0, const std::string varName = "",
        const NekDouble x2_in = NekConstants::kNekUnsetDouble,
        const NekDouble x3_in = NekConstants::kNekUnsetDouble) override;

    virtual const Array<OneD, const std::shared_ptr<ExpList>>
        &v_GetBndCondExpansions(void) override;

    virtual const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
        &v_GetBndConditions() override;

    virtual std::shared_ptr<ExpList> &v_UpdateBndCondExpansion(int i) override;

    virtual Array<OneD, SpatialDomains::BoundaryConditionShPtr>
        &v_UpdateBndConditions() override;

    virtual void v_SetBndCondBwdWeight(const int index,
                                       const NekDouble value) override;
};

typedef std::shared_ptr<DisContField3DHomogeneous2D>
    DisContField3DHomogeneous2DSharedPtr;

} // namespace MultiRegions
} // namespace Nektar

#endif // MULTIERGIONS_DISCONTFIELD3DHOMO2D_H
