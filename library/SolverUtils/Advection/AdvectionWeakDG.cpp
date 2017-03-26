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
            int i, cntx, cnty;
            int nquad0, nquad1;
            int physOffset;
            int nLocalPts;
            int nElements   = pFields[0]->GetExpSize();
            int nDimensions = pFields[0]->GetCoordim(0);
            int nTotalPts   = pFields[0]->GetTotPoints();
            int nTracePts   = pFields[0]->GetTrace()->GetTotPoints();
            
            // Jacobian and its determinant
            Array<TwoD, const NekDouble> gmat;
            Array<OneD, const NekDouble> jac;
            
            // Store forwards/backwards space along trace space
            m_dxFwd = Array<OneD, NekDouble>(nTracePts, 0.0);
            m_dxBwd = Array<OneD, NekDouble>(nTracePts, 0.0);
            m_dx    = Array<OneD, NekDouble>(nTotalPts, 0.0);
            m_dyFwd = Array<OneD, NekDouble>(nTracePts, 0.0);
            m_dyBwd = Array<OneD, NekDouble>(nTracePts, 0.0);
            m_dy    = Array<OneD, NekDouble>(nTotalPts, 0.0);
            
            // Auxiliary array
            Array<OneD, NekDouble> auxArray1;
            
            // Base and point information
            Array<OneD, LibUtilities::BasisSharedPtr> base;
            LibUtilities::PointsKeyVector             ptsKeys;
            
            // Load local coords into coords
            Array<OneD, Array<OneD, NekDouble> > coords(nDimensions);
            coords[0] = Array<OneD, NekDouble> (nTotalPts);
            
            switch (nDimensions)
            {
            case 1:
            {
                // Approach on physical space
                pFields[0]->GetCoords(coords[0]);
                
                for (int n = 0; n < nElements; ++n)
                {
                    base       = pFields[0]->GetExp(n)->GetBase();
                    nquad0     = base[0]->GetNumPoints();
                    nLocalPts  = pFields[0]->GetExp(n)->GetTotPoints();
                    physOffset = pFields[0]->GetPhys_Offset(n);

                    for (i = 0; i < nLocalPts; ++i)
                    {
                        std::cout << "i  = "    << i
                        << ",    nLocalPts = "  << nLocalPts
                        << ",    physOffset = " << physOffset
                        << ",    nquad0 = "     << nquad0
                        << std::endl;
                        
                        if (i == nquad0-1)
                        {
                            m_dx[i+physOffset] =
                                (coords[0][i+physOffset] -
                                 coords[0][i-1+physOffset]);
                            std::cout << "IF = " << m_dx[i+physOffset] << std::endl;
                            std::cout << "i  = " << i
                                      << ",    physOffset = " << physOffset
                                      << ",    nquad0 = "     << nquad0
                                      << std::endl;
                        }
                        else
                        {
                            m_dx[i+physOffset] =
                                (coords[0][i+physOffset] -
                                 coords[0][i+1+physOffset]);
                            
                            std::cout << "ELSE = " << m_dx[i+physOffset] << std::endl;
                        }
                    }
                }
                Vmath::Vabs(nTotalPts, m_dx, 1, m_dx, 1);
                pFields[0]->GetFwdBwdTracePhys(m_dx, m_dxFwd, m_dxBwd);

                for (i = 0; i < nTotalPts; ++i)
                {
                    std::cout << "x    = " << coords[0][i] << std::endl;
                    std::cout << "m_dx = " << m_dx[i]      << std::endl;
                }
                
                // Approach on standard space
                /*
                for (int n = 0; n < nElements; ++n)
                {
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
                 */
                break;
            }
            case 2:
            {
                coords[1] = Array<OneD, NekDouble> (nTotalPts);

                // Approach on physical space
                pFields[0]->GetCoords(coords[0], coords[1]);
                
                for (int n = 0; n < nElements; ++n)
                {
                    cntx = 0;
                    cnty = 1;
                    base       = pFields[0]->GetExp(n)->GetBase();
                    nquad0     = base[0]->GetNumPoints();
                    nquad1     = base[1]->GetNumPoints();
                    nLocalPts  = pFields[0]->GetExp(n)->GetTotPoints();
                    physOffset = pFields[0]->GetPhys_Offset(n);
                    
                    for (i = 0; i < nLocalPts; ++i)
                    {
                        std::cout << "i  = "    << i
                        << ",    nLocalPts = "  << nLocalPts
                        << ",    physOffset = " << physOffset
                        << ",    nquad0 = "     << nquad0
                        << std::endl;
                        
                        if (i == (nquad0-1)+cntx*nquad0 && i != nLocalPts)
                        {
                            cnty++;
                        }
                        if (i == (nquad0-1)+cntx*nquad0)
                        {
                            cntx++;

                            m_dx[i+physOffset] =
                            (coords[0][i+physOffset] -
                             coords[0][i-1+physOffset]);
                            std::cout << "IF = " << m_dx[i+physOffset]
                                                 << std::endl;
                            std::cout << "i  = " << i
                            << ",    physOffset = " << physOffset
                            << ",    nquad0 = "     << nquad0
                                                    << std::endl;
                        }
                        else
                        {
                            m_dx[i+physOffset] =
                            (coords[0][i+physOffset] -
                             coords[0][i+1+physOffset]);
                            
                            std::cout << "ELSE = " << m_dx[i+physOffset] << std::endl;
                        }
                        
                        if (cnty == nquad1)
                        {
                            std::cout << "cnty = " << cnty << std::endl;
                            m_dy[i+physOffset] = (coords[1][i+physOffset] -
                                                  coords[1][i-nquad0+physOffset]);
                        }
                        else{
                            m_dy[i+physOffset] =
                                (coords[1][i+nquad0+physOffset] -
                                 coords[1][i+physOffset]);
                        }
                        std::cout << "ELSE DY = " << m_dy[i+physOffset] << std::endl;
                    }
                }
                Vmath::Vabs(nTotalPts, m_dx, 1, m_dx, 1);
                Vmath::Vabs(nTotalPts, m_dy, 1, m_dy, 1);

                for (i = 0; i < nTotalPts; ++i)
                {
                    m_dx[i] = std::min(m_dx[i], m_dy[i]);
                }
                pFields[0]->GetFwdBwdTracePhys(m_dx, m_dxFwd, m_dxBwd);
                
                for (i = 0; i < nTotalPts; ++i)
                {
                    std::cout << "x    = " << coords[0][i] << std::endl;
                    std::cout << "y    = " << coords[1][i] << std::endl;
                    std::cout << "m_dx - min = " << m_dx[i]      << std::endl;
                    std::cout << "m_dy = " << m_dx[i]      << std::endl;

                }
                
                // Approach on standard space
                /*
                LibUtilities::PointsKeyVector ptsKeys;
                for (int n = 0; n < nElements; ++n)
                {
                    base        = pFields[0]->GetExp(n)->GetBase();
                    nquad0      = base[0]->GetNumPoints();
                    nquad1      = base[1]->GetNumPoints();
                    
                    pFields[0]->GetCoords(coords[0], coords[1]);
                    
                    for (i = 0; i < nTotalPts; ++i)
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
                        for (i = 0; i < nLocalPts; ++i)
                        {
                            m_dx[i+physOffset] = 2.0 * (z0[1]-z0[0])*jac[i];
                        }
                    }
                    else
                    {
                        for (i = 0; i < nLocalPts; ++i)
                        {
                            m_dx[i+physOffset] = 2.0 * (z0[1]-z0[0])*jac[n];
                        }
                    }
                }
                 */
                break;
            }
            case 3:
            {
                coords[1] = Array<OneD, NekDouble> (nTotalPts);
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
            
            m_riemann->Solve(m_spaceDim, Fwd, Bwd, numflux, m_dxFwd, m_dxBwd);

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
