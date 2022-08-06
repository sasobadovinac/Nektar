//////////////////////////////////////////////////////////////////////////////
//
// File DisContField3DHomogeneous2D.cpp
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
// Description: Field definition for 3D domain with boundary
// conditions using LDG flux and a 2D homogeneous directions
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <MultiRegions/DisContField.h>
#include <MultiRegions/DisContField3DHomogeneous2D.h>
#include <MultiRegions/ExpList1DHomogeneous2D.h>

namespace Nektar
{
namespace MultiRegions
{

DisContField3DHomogeneous2D::DisContField3DHomogeneous2D(void)
    : ExpList3DHomogeneous2D(), m_bndCondExpansions(), m_bndCondBndWeight(),
      m_bndConditions()
{
}

DisContField3DHomogeneous2D::DisContField3DHomogeneous2D(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::BasisKey &HomoBasis_y,
    const LibUtilities::BasisKey &HomoBasis_z, const NekDouble lhom_y,
    const NekDouble lhom_z, const bool useFFT, const bool dealiasing,
    const Collections::ImplementationType ImpType)
    : ExpList3DHomogeneous2D(pSession, HomoBasis_y, HomoBasis_z, lhom_y, lhom_z,
                             useFFT, dealiasing, ImpType),
      m_bndCondExpansions(), m_bndCondBndWeight(), m_bndConditions()
{
}

DisContField3DHomogeneous2D::DisContField3DHomogeneous2D(
    const DisContField3DHomogeneous2D &In, const bool DeclareLinesSetCoeffPhys)
    : ExpList3DHomogeneous2D(In, false),
      m_bndCondExpansions(In.m_bndCondExpansions),
      m_bndCondBndWeight(In.m_bndCondBndWeight),
      m_bndConditions(In.m_bndConditions)
{
    if (DeclareLinesSetCoeffPhys)
    {
        DisContFieldSharedPtr zero_line =
            std::dynamic_pointer_cast<DisContField>(In.m_lines[0]);

        for (int n = 0; n < m_lines.size(); ++n)
        {
            m_lines[n] =
                MemoryManager<DisContField>::AllocateSharedPtr(*zero_line);
        }

        SetCoeffPhys();
    }
}

DisContField3DHomogeneous2D::DisContField3DHomogeneous2D(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::BasisKey &HomoBasis_y,
    const LibUtilities::BasisKey &HomoBasis_z, const NekDouble lhom_y,
    const NekDouble lhom_z, const bool useFFT, const bool dealiasing,
    const SpatialDomains::MeshGraphSharedPtr &graph1D,
    const std::string &variable, const Collections::ImplementationType ImpType)
    : ExpList3DHomogeneous2D(pSession, HomoBasis_y, HomoBasis_z, lhom_y, lhom_z,
                             useFFT, dealiasing, ImpType),
      m_bndCondExpansions(), m_bndCondBndWeight(), m_bndConditions()
{
    int i, n, nel;
    DisContFieldSharedPtr line_zero;
    SpatialDomains::BoundaryConditions bcs(pSession, graph1D);

    //
    m_lines[0] = line_zero = MemoryManager<DisContField>::AllocateSharedPtr(
        pSession, graph1D, variable, ImpType);

    m_exp = MemoryManager<LocalRegions::ExpansionVector>::AllocateSharedPtr();
    nel   = m_lines[0]->GetExpSize();

    for (i = 0; i < nel; ++i)
    {
        (*m_exp).push_back(m_lines[0]->GetExp(i));
    }

    int nylines = m_homogeneousBasis_y->GetNumPoints();
    int nzlines = m_homogeneousBasis_z->GetNumPoints();

    for (n = 1; n < nylines * nzlines; ++n)
    {
        m_lines[n] = MemoryManager<DisContField>::AllocateSharedPtr(
            pSession, graph1D, variable, ImpType);
        for (i = 0; i < nel; ++i)
        {
            (*m_exp).push_back((*m_exp)[i]);
        }
    }

    // Setup Default optimisation information.
    nel = GetExpSize();

    SetCoeffPhys();

    SetupBoundaryConditions(HomoBasis_y, HomoBasis_z, lhom_y, lhom_z, bcs,variable);
}

DisContField3DHomogeneous2D::~DisContField3DHomogeneous2D()
{
}

void DisContField3DHomogeneous2D::SetupBoundaryConditions(
    const LibUtilities::BasisKey &HomoBasis_y,
    const LibUtilities::BasisKey &HomoBasis_z, const NekDouble lhom_y,
    const NekDouble lhom_z, SpatialDomains::BoundaryConditions &bcs,
    const std::string variable)
{
    // Setup an ExpList1DHomogeneous2D expansion for boundary
    // conditions and link to class declared in m_lines.

    size_t nlines = m_lines.size();

    const SpatialDomains::BoundaryRegionCollection &bregions =
        bcs.GetBoundaryRegions();

    size_t nbnd = bregions.size();

    m_bndCondExpansions = Array<OneD, MultiRegions::ExpListSharedPtr>(nbnd);
    m_bndCondBndWeight  = Array<OneD, NekDouble>{nbnd, 0.0};

    Array<OneD, MultiRegions::ExpListSharedPtr> LinesBndCondExp(nlines);

    m_bndConditions = m_lines[0]->UpdateBndConditions();
#if EXPLISTDATA
#else            
    int bndsize = m_bndCondExpansions.size();
    m_bndCondFieldCoeff = Array<OneD, NekFieldCoeffSharedPtr>(bndsize);
    m_bndCondFieldPhys  = Array<OneD, NekFieldPhysSharedPtr>(bndsize);

    for (int n = 0; n < nlines; ++n)
    {
        std::dynamic_pointer_cast<DisContField>(m_lines[n])
            ->m_bndCondFieldCoeff =
            Array<OneD, NekFieldCoeffSharedPtr>(bndsize);
    }
#endif

    for (int i = 0; i < nbnd; ++i)
    {
        for (int n = 0; n < nlines; ++n)
        {
            LinesBndCondExp[n] = m_lines[n]->UpdateBndCondExpansion(i);
        }

        m_bndCondExpansions[i] =
            MemoryManager<ExpList1DHomogeneous2D>::AllocateSharedPtr(
                m_session, HomoBasis_y, HomoBasis_z, lhom_y, lhom_z, m_useFFT,
                false, LinesBndCondExp);

#if EXPLISTDATA
#else        
        m_bndCondFieldCoeff[i] = std::make_shared<
            NekField<NekDouble, eCoeff>>(m_bndCondExpansions[i]);
        m_bndCondFieldPhys[i] = std::make_shared<
            NekField<NekDouble, ePhys>>(m_bndCondExpansions[i]);

        // point m_bndCondFieldCoeff in m_planes to member of 1D array
        // attached to m_bndCondFieldCoeff - very hacky!
        int bnd_ncoeffs = LinesBndCondExp[0]->GetNcoeffs(); 
        for (int n = 0; n < nlines; ++n)
        {
            DisContFieldSharedPtr dgfield =
                std::dynamic_pointer_cast<DisContField>(m_lines[n]);

            dgfield->m_bndCondFieldCoeff[i] =
                std::make_shared<NekField<NekDouble, eCoeff>>
                (dgfield->m_bndCondExpansions[i]);
            dgfield->m_bndCondFieldCoeff[i]->UpdateArray1D() =
                m_bndCondFieldCoeff[i]->UpdateArray1D() + n*bnd_ncoeffs;
        }
#endif
    }

    v_EvaluateBoundaryConditions(0.0, variable);
}

void DisContField3DHomogeneous2D::v_HelmSolve(
    const Array<OneD, const NekDouble> &inarray,
    Array<OneD, NekDouble> &outarray, const StdRegions::ConstFactorMap &factors,
    const StdRegions::VarCoeffMap &varcoeff,
    const MultiRegions::VarFactorsMap &varfactors,
    const Array<OneD, const NekDouble> &dirForcing, const bool PhysSpaceForcing)
{
    int n, m;
    int cnt          = 0;
    int cnt1         = 0;
    int nhom_modes_y = m_homogeneousBasis_y->GetNumModes();
    int nhom_modes_z = m_homogeneousBasis_z->GetNumModes();
    NekDouble beta_y;
    NekDouble beta_z;
    StdRegions::ConstFactorMap new_factors;

    Array<OneD, NekDouble> e_out;
    Array<OneD, NekDouble> fce(inarray.size());
    Array<OneD, const NekDouble> wfce;

    // Fourier transform forcing function
    if (m_WaveSpace)
    {
        fce = inarray;
    }
    else
    {
        HomogeneousFwdTrans(inarray, fce);
    }

    for (n = 0; n < nhom_modes_z; ++n)
    {
        for (m = 0; m < nhom_modes_y; ++m)
        {
            beta_z      = 2 * M_PI * (n / 2) / m_lhom_z;
            beta_y      = 2 * M_PI * (m / 2) / m_lhom_y;
            new_factors = factors;
            new_factors[StdRegions::eFactorLambda] +=
                beta_y * beta_y + beta_z * beta_z;

            wfce = (PhysSpaceForcing) ? fce + cnt : fce + cnt1;
            m_lines[n]->HelmSolve(wfce, e_out = outarray + cnt1, new_factors,
                                  varcoeff, varfactors, dirForcing,
                                  PhysSpaceForcing);

            cnt += m_lines[n]->GetTotPoints();
            cnt1 += m_lines[n]->GetNcoeffs();
        }
    }
}

#if 1
void DisContField3DHomogeneous2D::v_EvaluateBoundaryConditions(
    const NekDouble time, const std::string varName, const NekDouble x2_in,
    const NekDouble x3_in)
{
    boost::ignore_unused(x2_in, x3_in);
    int i;
    int npoints;
    int nbnd = m_bndCondExpansions.size();
    MultiRegions::ExpListSharedPtr locExpList;

    for (i = 0; i < nbnd; ++i)
    {
        if (time == 0.0 || m_bndConditions[i]->IsTimeDependent())
        {
            locExpList = m_bndCondExpansions[i];
            npoints    = locExpList->GetNpoints();

            Array<OneD, NekDouble> x0(npoints, 0.0);
            Array<OneD, NekDouble> x1(npoints, 0.0);
            Array<OneD, NekDouble> x2(npoints, 0.0);
            Array<OneD, NekDouble> valuesFile(npoints, 1.0),
                valuesExp(npoints, 1.0);

            locExpList->GetCoords(x0, x1, x2);

            if (m_bndConditions[i]->GetBoundaryConditionType() ==
                SpatialDomains::eDirichlet)
            {
                SpatialDomains::DirichletBCShPtr bcPtr =
                    std::static_pointer_cast<
                        SpatialDomains::DirichletBoundaryCondition>(
                        m_bndConditions[i]);
                std::string bcfilename  = bcPtr->m_filename;
                std::string exprbcs     = bcPtr->m_expr;

                if (bcfilename != "")
                {
#if EXPLISTDATA
                    locExpList->ExtractCoeffsFromFile
                        (bcfilename, bcPtr->GetComm(), varName,
                         locExpList->UpdateCoeffs());
                    locExpList->BwdTrans(locExpList->GetCoeffs(),
                                         locExpList->UpdatePhys());
                    
                    valuesFile = locExpList->GetPhys();
#else
                    locExpList->ExtractCoeffsFromFile
                        (bcfilename, bcPtr->GetComm(), varName,
                         m_bndCondFieldCoeff[i]->UpdateArray1D());
                    
                    valuesFile = m_bndCondFieldPhys[i]->UpdateArray1D();
                    locExpList->BwdTrans(m_bndCondFieldCoeff[i]->GetArray1D(),
                                         valuesFile);
#endif
                }

                if (exprbcs != "")
                {
                    LibUtilities::Equation condition =
                        std::static_pointer_cast<
                            SpatialDomains::DirichletBoundaryCondition>(
                            m_bndConditions[i])
                            ->m_dirichletCondition;

                    condition.Evaluate(x0, x1, x2, time, valuesExp);
                }

#if EXPLISTDATA
                Vmath::Vmul(npoints, valuesExp, 1, valuesFile, 1,
                            locExpList->UpdatePhys(), 1);
                
                // set wave space to false since have set up phys values
                locExpList->SetWaveSpace(false);
                
                locExpList->FwdTransBndConstrained(
                        locExpList->GetPhys(), locExpList->UpdateCoeffs());
#else
                Vmath::Vmul(npoints, valuesExp, 1, valuesFile, 1,
                            m_bndCondFieldPhys[i]->UpdateArray1D(),1);
                
                // set wave space to false since have set up phys values
                locExpList->SetWaveSpace(false);
                
                locExpList->FwdTransBndConstrained(
                                    m_bndCondFieldPhys[i]->GetArray1D(),
                                    m_bndCondFieldCoeff[i]->UpdateArray1D());
#endif 

            }
            else if (m_bndConditions[i]->GetBoundaryConditionType() ==
                     SpatialDomains::eNeumann)
            {
                SpatialDomains::NeumannBCShPtr bcPtr = std::static_pointer_cast<
                    SpatialDomains::NeumannBoundaryCondition>(
                    m_bndConditions[i]);

                std::string bcfilename = bcPtr->m_filename;
                
                if (bcfilename != "")
                {
#if EXPLISTDATA
                        locExpList->ExtractCoeffsFromFile
                            (bcfilename, bcPtr->GetComm(), varName,
                             locExpList->UpdateCoeffs());
                        locExpList->BwdTrans(locExpList->GetCoeffs(),
                                             locExpList->UpdatePhys());
#else
                        locExpList->ExtractCoeffsFromFile
                            (bcfilename, bcPtr->GetComm(), varName,
                             m_bndCondFieldCoeff[i]->UpdateArray1D());

                        locExpList->BwdTrans(m_bndCondFieldCoeff[i]->GetArray1D(),
                                           m_bndCondFieldPhys[i]->UpdateArray1D());
#endif
                }
                else
                {
                    LibUtilities::Equation condition =
                        std::static_pointer_cast<
                            SpatialDomains::NeumannBoundaryCondition>(
                            m_bndConditions[i])
                            ->m_neumannCondition;

#if EXPLISTDATA
                    condition.Evaluate(x0, x1, x2, time,
                                       locExpList->UpdatePhys());
#else
                    condition.Evaluate(x0, x1, x2, time,
                                       m_bndCondFieldPhys[i]->UpdateArray1D());
#endif
                }

#if EXPLISTDATA
                    locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                locExpList->UpdateCoeffs());
#else
                    locExpList->IProductWRTBase(m_bndCondFieldPhys[i]->GetArray1D(),
                                           m_bndCondFieldCoeff[i]->UpdateArray1D());
#endif
            }
            else if (m_bndConditions[i]->GetBoundaryConditionType() ==
                     SpatialDomains::eRobin)
            {
                SpatialDomains::RobinBCShPtr bcPtr = std::static_pointer_cast<
                    SpatialDomains::RobinBoundaryCondition>(m_bndConditions[i]);

                std::string bcfilename = bcPtr->m_filename;
                
                if (bcfilename != "")
                {
#if EXPLISTDATA
                    locExpList->ExtractCoeffsFromFile
                        (bcfilename, bcPtr->GetComm(), varName,
                         locExpList->UpdateCoeffs());
                    locExpList->BwdTrans(locExpList->GetCoeffs(),
                                         locExpList->UpdatePhys());
#else
                    locExpList->ExtractCoeffsFromFile
                        (bcfilename, bcPtr->GetComm(), varName,
                         m_bndCondFieldCoeff[i]->UpdateArray1D());
                    
                    locExpList->BwdTrans(m_bndCondFieldCoeff[i]->GetArray1D(),
                                         m_bndCondFieldPhys[i]->UpdateArray1D());
#endif
                }
                else
                {
                    LibUtilities::Equation condition =
                        std::static_pointer_cast<
                            SpatialDomains::RobinBoundaryCondition>(
                            m_bndConditions[i])
                            ->m_robinFunction;

#if EXPLISTDATA
                    condition.Evaluate(x0, x1, x2, time,
                                       locExpList->UpdatePhys());
#else
                    condition.Evaluate(x0, x1, x2, time,
                                       m_bndCondFieldPhys[i]->UpdateArray1D());
#endif
                }

#if EXPLISTDATA
                locExpList->IProductWRTBase(locExpList->GetPhys(),
                                            locExpList->UpdateCoeffs());
#else
                locExpList->IProductWRTBase(m_bndCondFieldPhys[i]->GetArray1D(),
                                            m_bndCondFieldCoeff[i]->UpdateArray1D());
#endif
            }
            else if (m_bndConditions[i]->GetBoundaryConditionType() ==
                     SpatialDomains::ePeriodic)
            {
                continue;
            }
            else
            {
                ASSERTL0(false, "This type of BC not implemented yet");
            }
        }
    }
}
#else
    void DisContField3DHomogeneous2D::v_EvaluateBoundaryConditions(
    const NekDouble time, const std::string varName, const NekDouble x2_in,
    const NekDouble x3_in)
{
    boost::ignore_unused(x2_in, x3_in);
    int n, m;
    const Array<OneD, const NekDouble> y = m_homogeneousBasis_y->GetZ();
    const Array<OneD, const NekDouble> z = m_homogeneousBasis_z->GetZ();

    for (n = 0; n < m_nz; ++n)
    {
        for (m = 0; m < m_ny; ++m)
        {
            m_lines[m + (n * m_ny)]->EvaluateBoundaryConditions(
                time, varName, 0.5 * m_lhom_y * (1.0 + y[m]),
                0.5 * m_lhom_z * (1.0 + z[n]));
        }
    }

#if EXPLISTDATA
    // Fourier transform coefficient space boundary values
    for (n = 0; n < m_bndCondExpansions.size(); ++n)
    {
        if (time == 0.0 || m_bndConditions[n]->IsTimeDependent())
        {
            m_bndCondBndWeight[n] = 1.0;
            m_bndCondExpansions[n]->HomogeneousFwdTrans(
                m_bndCondExpansions[n]->GetCoeffs(),
                m_bndCondExpansions[n]->UpdateCoeffs());
        }
    }
#else
    // Fourier transform coefficient space boundary values
    for (n = 0; n < m_bndCondExpansions.size(); ++n)
    {
        if (time == 0.0 || m_bndConditions[n]->IsTimeDependent())
        {
            m_bndCondBndWeight[n] = 1.0;
            m_bndCondExpansions[n]->HomogeneousFwdTrans(
                m_bndCondFieldCoeff[n]->GetArray1D(),
                m_bndCondFieldPhys[n]->UpdateArray1D());
        }
    }
#endif
}
#endif

const Array<OneD, const std::shared_ptr<ExpList>>
    &DisContField3DHomogeneous2D::v_GetBndCondExpansions(void)
{
    return GetBndCondExpansions();
}

const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
    &DisContField3DHomogeneous2D::v_GetBndConditions()
{
    return GetBndConditions();
}

std::shared_ptr<ExpList> &DisContField3DHomogeneous2D::v_UpdateBndCondExpansion(
    int i)
{
    return UpdateBndCondExpansion(i);
}

Array<OneD, SpatialDomains::BoundaryConditionShPtr>
    &DisContField3DHomogeneous2D::v_UpdateBndConditions()
{
    return UpdateBndConditions();
}

void DisContField3DHomogeneous2D::GetBoundaryToElmtMap(Array<OneD, int> &ElmtID,
                                                       Array<OneD, int> &EdgeID)
{
    if (m_BCtoElmMap.size() == 0)
    {
        Array<OneD, int> ElmtID_tmp;
        Array<OneD, int> EdgeID_tmp;

        m_lines[0]->GetBoundaryToElmtMap(ElmtID_tmp, EdgeID_tmp);
        int nel_per_lines = m_lines[0]->GetExpSize();
        int nlines        = m_lines.size();

        int MapSize = ElmtID_tmp.size();

        m_BCtoElmMap = Array<OneD, int>(nlines * MapSize);
        m_BCtoEdgMap = Array<OneD, int>(nlines * MapSize);
        if (MapSize > 0)
        {
            int i, j, n, cnt;
            int cntLine = 0;
            for (cnt = n = 0; n < m_bndCondExpansions.size(); ++n)
            {
                int lineExpSize =
                    m_lines[0]->GetBndCondExpansions()[n]->GetExpSize();
                for (i = 0; i < lineExpSize; ++i, ++cntLine)
                {
                    for (j = 0; j < nlines; j++)
                    {
                        m_BCtoElmMap[cnt + i + j * lineExpSize] =
                            ElmtID_tmp[cntLine] + j * nel_per_lines;
                        m_BCtoEdgMap[cnt + i + j * lineExpSize] =
                            EdgeID_tmp[cntLine];
                    }
                }
                cnt += m_bndCondExpansions[n]->GetExpSize();
            }
        }
    }
    ElmtID = m_BCtoElmMap;
    EdgeID = m_BCtoEdgMap;
}

void DisContField3DHomogeneous2D::v_GetBndElmtExpansion(
    int i, std::shared_ptr<ExpList> &result, const bool DeclareCoeffPhysArrays)
{
    int n, cnt, nq;
    int offsetOld, offsetNew;

    std::vector<unsigned int> eIDs;
    Array<OneD, int> ElmtID, EdgeID;
    GetBoundaryToElmtMap(ElmtID, EdgeID);

    // Skip other boundary regions
    for (cnt = n = 0; n < i; ++n)
    {
        cnt += m_bndCondExpansions[n]->GetExpSize();
    }

    // Populate eIDs with information from BoundaryToElmtMap
    for (n = 0; n < m_bndCondExpansions[i]->GetExpSize(); ++n)
    {
        eIDs.push_back(ElmtID[cnt + n]);
    }

    // Create expansion list
    result =
        MemoryManager<ExpList3DHomogeneous2D>::AllocateSharedPtr(*this, eIDs);

    // Copy phys and coeffs to new explist
    if (DeclareCoeffPhysArrays)
    {
#if EXPLISTDATA
        Array<OneD, NekDouble> tmp1, tmp2;
        for (n = 0; n < result->GetExpSize(); ++n)
        {
            nq        = GetExp(ElmtID[cnt + n])->GetTotPoints();
            offsetOld = GetPhys_Offset(ElmtID[cnt + n]);
            offsetNew = result->GetPhys_Offset(n);
            Vmath::Vcopy(nq, tmp1 = GetPhys() + offsetOld, 1,
                         tmp2 = result->UpdatePhys() + offsetNew, 1);

            nq        = GetExp(ElmtID[cnt + n])->GetNcoeffs();
            offsetOld = GetCoeff_Offset(ElmtID[cnt + n]);
            offsetNew = result->GetCoeff_Offset(n);
            Vmath::Vcopy(nq, tmp1 = GetCoeffs() + offsetOld, 1,
                         tmp2 = result->UpdateCoeffs() + offsetNew, 1);
        }
#else
        boost::ignore_unused(offsetOld,offsetNew,nq);
        NEKERROR(ErrorUtil::efatal,
                 "This method needs updating for FieldStorage usage");
        
#endif
    }

    // Set wavespace value
    result->SetWaveSpace(GetWaveSpace());
}

} // namespace MultiRegions
} // namespace Nektar
