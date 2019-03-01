///////////////////////////////////////////////////////////////////////////////
//
// File DisContField2D.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Field definition in two-dimensions for a discontinuous LDG-H
// expansion.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD2D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD2D_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/DisContField.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/AssemblyMap/LocTraceToTraceMap.h>
#include <SpatialDomains/Conditions.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class DisContField2D : public DisContField
        {
        public:
            MULTI_REGIONS_EXPORT DisContField2D();

            MULTI_REGIONS_EXPORT DisContField2D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr   &graph2D,
                const std::string                          &variable,
                const bool SetUpJustDG            = true,
                const bool DeclareCoeffPhysArrays = true,
                const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);
            
            MULTI_REGIONS_EXPORT DisContField2D(
                const DisContField2D                     &In,
                const SpatialDomains::MeshGraphSharedPtr &graph2D,
                const std::string                        &variable,
                const bool SetUpJustDG            = false,
                const bool DeclareCoeffPhysArrays = true);
            
            MULTI_REGIONS_EXPORT DisContField2D(
                const DisContField2D &In,
                const bool DeclareCoeffPhysArrays = true);

            MULTI_REGIONS_EXPORT virtual ~DisContField2D();
            
            MULTI_REGIONS_EXPORT NekDouble L2_DGDeriv(
                const int                           dir,
                const Array<OneD, const NekDouble> &soln);

            MULTI_REGIONS_EXPORT void EvaluateHDGPostProcessing(
                Array<OneD, NekDouble> &outarray);
            
            virtual ExpListSharedPtr &v_GetTrace()
            {
                if(m_trace == NullExpListSharedPtr)
                {
                    SetUpDG();
                }
                
                return m_trace;
            }

            Array<OneD, int> m_BCtoElmMap;
            Array<OneD, int> m_BCtoEdgMap;

        protected:

            Array<OneD, LibUtilities::BasisSharedPtr> m_base; /**< Bases needed for the expansion */

            /** \brief This function gets the shared point to basis
             *
             *  \return returns the shared pointer to the bases
             */
            inline const Array<OneD, const LibUtilities::BasisSharedPtr>& GetBase() const
            {
                return(m_base);
            }

            /** \brief This function returns the type of basis used in
             *  the \a dir direction
             *
             *  The different types of bases implemented in the code
             *  are defined in the LibUtilities::BasisType enumeration
             *  list. As a result, the function will return one of the
             *  types of this enumeration list.
             *
             *  \param dir the direction \return returns the type of
             *  basis used in the \a dir direction
             */
            inline  LibUtilities::BasisType GetBasisType(const int dir) const
            {
                ASSERTL1(dir < m_base.num_elements(), "dir is larger than m_numbases");
                return(m_base[dir]->GetBasisType());
            }


            Array<OneD, Array<OneD, unsigned int> > m_mapEdgeToElmn;
            Array<OneD, Array<OneD, unsigned int> > m_signEdgeToElmn;
            Array<OneD,StdRegions::Orientation>     m_edgedir;

            virtual void v_AddFwdBwdTraceIntegral(
                const Array<OneD, const NekDouble> &Fwd, 
                const Array<OneD, const NekDouble> &Bwd, 
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_GeneralMatrixOp(
                const GlobalMatrixKey             &gkey,
                const Array<OneD,const NekDouble> &inarray,
                Array<OneD,      NekDouble> &outarray);
            virtual void v_GetBoundaryToElmtMap(
                Array<OneD, int> &ElmtID,
                Array<OneD, int> &EdgeID);
            virtual void v_GetBndElmtExpansion(int i,
                            std::shared_ptr<ExpList> &result,
                            const bool DeclareCoeffPhysArrays);
            virtual void v_Reset();

            /**
             * @brief Obtain a copy of the periodic edges and vertices for this
             * field.
             */
            virtual void v_GetPeriodicEntities(
                PeriodicMap &periodicVerts,
                PeriodicMap &periodicEdges,
                PeriodicMap &periodicFaces)
            {
                periodicVerts = m_periodicVerts;
                periodicEdges = m_periodicEdges;
            }

            
            virtual AssemblyMapDGSharedPtr &v_GetTraceMap()
            {
                return m_traceMap;
            }

            virtual const Array<OneD,const MultiRegions::ExpListSharedPtr>
                &v_GetBndCondExpansions()
            {
                return m_bndCondExpansions;
            }

            virtual const 
                Array<OneD,const SpatialDomains::BoundaryConditionShPtr>
                &v_GetBndConditions()
            {
                return m_bndConditions;
            }

            virtual MultiRegions::ExpListSharedPtr 
                &v_UpdateBndCondExpansion(int i)
            {
                return m_bndCondExpansions[i];
            }

            virtual Array<OneD, SpatialDomains::BoundaryConditionShPtr> 
                &v_UpdateBndConditions()
            {
                return m_bndConditions;
            }

            virtual std::map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo();
        };
        
        typedef std::shared_ptr<DisContField2D>   DisContField2DSharedPtr;
    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD2D_H
