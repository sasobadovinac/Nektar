///////////////////////////////////////////////////////////////////////////////
//
// File: AdvectionWeakDG.cpp
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
// Description: Weak DG advection class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Advection/AdvectionWeakDG.h>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string AdvectionWeakDG::type = GetAdvectionFactory().
            RegisterCreatorFunction("WeakDG", AdvectionWeakDG::create);

        AdvectionWeakDG::AdvectionWeakDG()
        {
        }

        /**
         * @brief Initialise AdvectionWeakDG objects and store them before
         * starting the time-stepping.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void AdvectionWeakDG::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            Advection::v_InitObject(pSession, pFields);
            v_SetupFluxLength      (pSession, pFields);
        }
        

        
        /**
         * @brief Setup the metric terms to compute the contravariant
         * fluxes. (i.e. this special metric terms transform the fluxes
         * at the interfaces of each element from the physical space to
         * the standard space).
         *
         * This routine calls the function #GetEdgeQFactors to compute and
         * store the metric factors following an anticlockwise conventions
         * along the edges/faces of each element. Note: for 1D problem
         * the transformation is not needed.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         *
         * \todo Add the metric terms for 3D Hexahedra.
         */
        void AdvectionWeakDG::v_SetupFluxLength(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            int nquad0, nquad1;
            int physOffset;
            int nLocalPts;
            int nElements   = pFields[0]->GetExpSize();
            int nDimensions = pFields[0]->GetCoordim(0);
            int nTotalPts   = pFields[0]->GetTotPoints();
            int nTracePts   = pFields[0]->GetTrace()->GetTotPoints();
            

            Array<TwoD, const NekDouble> gmat;
            Array<OneD, const NekDouble> jac;
            
            m_dx = Array<OneD, NekDouble>(nTotalPts);
            
            // Auxiliary array
            Array<OneD, NekDouble> auxArray1;
            
            // Base and point information
            Array<OneD, LibUtilities::BasisSharedPtr> base;
            LibUtilities::PointsKeyVector             ptsKeys;
            
            // Load local coords into coords
            Array<OneD, Array<OneD, NekDouble> > coords(nDimensions);
            coords[0] = Array<OneD, NekDouble> (nTotalPts);

            Array<OneD, NekDouble> tmpDX(nTotalPts);
            
            switch (nDimensions)
            {
            case 1:
            {
                for (int n = 0; n < nElements; ++n)
                {
                    pFields[0]->GetCoords(coords[0]);
                    
                    for (int i = 0; i < nTotalPts-1; ++i)
                    {
                        tmpDX[i] = (coords[0][i] - coords[0][i+1]);
                    }
                    tmpDX[nTotalPts] = (coords[0][nTotalPts-1] - coords[0][nTotalPts]);
                    
                    for (int i = 0; i < nTotalPts; ++i)
                    {
                        std::cout << "X     = " << coords[0][i] << std::endl;
                        std::cout << "tmpDX = " << tmpDX[i]     << std::endl;
                    }
                        
                    ptsKeys   = pFields[0]->GetExp(n)->GetPointsKeys();
                    nLocalPts = pFields[0]->GetExp(n)->GetTotPoints();
    
                    physOffset = pFields[0]->GetPhys_Offset(n);
                    base       = pFields[0]->GetExp(n)->GetBase();

                    jac = pFields[0]->GetExp(n)->
                            as<LocalRegions::Expansion1D>()->GetGeom1D()->
                                GetMetricInfo()->GetJac(ptsKeys);
                        
                    Array<OneD, const NekDouble> z0;
                    Array<OneD, const NekDouble> w0;
                    
                    base[0]->GetZW(z0, w0);

                    for (int i = 0; i < nLocalPts; ++i)
                    {
                        m_dx[i+physOffset] = 2.0 * (z0[1]-z0[0])*jac[n];
                    }
                }
                break;
            }
            case 2:
            {
                LibUtilities::PointsKeyVector ptsKeys;
                    
                for (int n = 0; n < nElements; ++n)
                {
                    base        = pFields[0]->GetExp(n)->GetBase();
                    nquad0      = base[0]->GetNumPoints();
                    nquad1      = base[1]->GetNumPoints();
                    
                    coords[1] = Array<OneD, NekDouble> (nTotalPts);
                    pFields[0]->GetCoords(coords[0], coords[1]);
                    
                    for (int i = 0; i < nTotalPts; ++i)
                    {
                        std::cout << "X = " << coords[0][i]
                        << "Y = " << coords[1][i] << std::endl;
                    }
                    
                    Array<OneD, const NekDouble> z0;
                    Array<OneD, const NekDouble> w0;   
                    Array<OneD, const NekDouble> z1;
                    Array<OneD, const NekDouble> w1;
                        
                    base[0]->GetZW(z0, w0);
                    base[1]->GetZW(z1, w1);

                    ptsKeys = pFields[0]->GetExp(n)->GetPointsKeys();
                    nLocalPts = pFields[0]->GetExp(n)->GetTotPoints();
                    physOffset = pFields[0]->GetPhys_Offset(n);
                        
                    jac  = pFields[0]->GetExp(n)
                        ->as<LocalRegions::Expansion2D>()->GetGeom2D()
                        ->GetMetricInfo()->GetJac(ptsKeys);

                        
                    if (pFields[0]->GetExp(n)->as<LocalRegions::Expansion2D>()->
                        GetGeom2D()->GetMetricInfo()->GetGtype()
                        == SpatialDomains::eDeformed)
                    {
                        for (int i = 0; i < nLocalPts; ++i)
                        {
                            m_dx[i+physOffset] = 2.0 * (z0[1]-z0[0])*jac[i];
                        }
                    }
                    else
                    {
                        for (int i = 0; i < nLocalPts; ++i)
                        {
                            m_dx[i+physOffset] = 2.0 * (z0[1]-z0[0])*jac[n];
                        }
                    }
                }
                break;
            }
            case 3:
            {
                coords[2] = Array<OneD, NekDouble> (nTotalPts);
                pFields[0]->GetCoords(coords[0], coords[1], coords[3]);
                ASSERTL0(false,"3D flux-lenght terms not implemented (yet)");
                break;
            }
            default:
            {
                ASSERTL0(false, "Expansion dimension not recognised");
                break;
            }
            }
        }

        /**
         * @brief Compute the advection term at each time-step using the
         * Discontinuous Glaerkin approach (DG).
         *
         * @param nConvectiveFields   Number of fields.
         * @param fields              Pointer to fields.
         * @param advVel              Advection velocities.
         * @param inarray             Solution at the previous time-step.
         * @param outarray            Advection term to be passed at the
         *                            time integration class.
         */
        void AdvectionWeakDG::v_Advect(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const NekDouble                                   &time,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            int nPointsTot      = fields[0]->GetTotPoints();
            int nCoeffs         = fields[0]->GetNcoeffs();
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();
            int nElements       = fields[0]->GetExpSize();

            int i, j;
            
            Array<OneD, Array<OneD, NekDouble> > tmp(nConvectiveFields);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > fluxvector(
                nConvectiveFields);


            std::cout << std::setprecision(16);
            
            // Allocate storage for flux vector F(u).
            for (i = 0; i < nConvectiveFields; ++i)
            {
                fluxvector[i] =
                    Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
                for (j = 0; j < m_spaceDim; ++j)
                {
                    fluxvector[i][j] = Array<OneD, NekDouble>(nPointsTot);
                }
            }

            ASSERTL1(m_riemann,
                     "Riemann solver must be provided for AdvectionWeakDG.");

            m_fluxVector(inarray, fluxvector);

            // Get the advection part (without numerical flux)
            for(i = 0; i < nConvectiveFields; ++i)
            {
                tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);

                fields[i]->IProductWRTDerivBase(fluxvector[i],tmp[i]);
            }

            // Store forwards/backwards space along trace space
            Array<OneD, Array<OneD, NekDouble> > Fwd    (nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > Bwd    (nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > numflux(nConvectiveFields);
            
            if (pFwd == NullNekDoubleArrayofArray ||
                pBwd == NullNekDoubleArrayofArray)
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    Bwd[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    numflux[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
                }
            }
            else
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd[i]     = pFwd[i];
                    Bwd[i]     = pBwd[i];
                    numflux[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                }
            }
            
            m_riemann->Solve(m_spaceDim, Fwd, Bwd, numflux, m_dx);

            // Evaulate <\phi, \hat{F}\cdot n> - OutField[i]
            for(i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Neg                      (nCoeffs, tmp[i], 1);
                fields[i]->AddTraceIntegral     (numflux[i], tmp[i]);
                fields[i]->MultiplyByElmtInvMass(tmp[i], tmp[i]);
                fields[i]->BwdTrans             (tmp[i], outarray[i]);
            }
        }
    }//end of namespace SolverUtils
}//end of namespace Nektar
