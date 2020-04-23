///////////////////////////////////////////////////////////////////////////////
//
// File TimeIntegrationScheme.cpp
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
// Description: implementation of time integration key class
//
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>
#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>  

namespace Nektar
{
    namespace LibUtilities 
    {  
        TimeIntegrationSchemeManagerT &TimeIntegrationSchemeManager(void)
        {
            static TimeIntegrationSchemeManagerT instance;
            instance.RegisterGlobalCreator(TimeIntegrationScheme::Create);
            return instance;
        }
        

        bool operator==(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs)
        {
            return (lhs.m_method == rhs.m_method);
        }
        
        bool operator<(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs)
        {
            return (lhs.m_method < rhs.m_method);
        }
        
        bool TimeIntegrationSchemeKey::opLess::operator()(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs) const
        {
            return (lhs.m_method < rhs.m_method);
        }

        std::ostream& operator<<(std::ostream& os, const TimeIntegrationSchemeKey& rhs)
        {
            os << "Time Integration Scheme: " << TimeIntegrationMethodMap[rhs.GetIntegrationMethod()] << endl;

            return os;
        }

        TimeIntegrationSolution::TimeIntegrationSolution(const TimeIntegrationSchemeKey &key, 
                                                         const DoubleArray& y, 
                                                         const NekDouble time, 
                                                         const NekDouble timestep):
            m_scheme(TimeIntegrationSchemeManager()[key]),
            m_solVector(m_scheme->GetNsteps()),
            m_t(m_scheme->GetNsteps())
        {
            m_solVector[0] = y;
            m_t[0] = time;

            int nsteps         = m_scheme->GetNsteps();
            int nvar           = y.num_elements();
            int npoints        = y[0].num_elements();
            int nMultiStepVals = m_scheme->GetNmultiStepValues();
            const Array<OneD, const unsigned int>& timeLevels = GetTimeLevelOffset(); 
            for(int i = 1; i < nsteps; i++)
            {
                m_solVector[i] = Array<OneD, Array<OneD, NekDouble> >(nvar);
                for(int j = 0; j < nvar; j++)
                {
                    m_solVector[i][j] = Array<OneD,NekDouble>(npoints,0.0);
                }
                if(i < nMultiStepVals)
                {
                    m_t[i] = time - i*timestep*timeLevels[i];
                }
                else
                {
                    m_t[i] = timestep;
                }
            }
        }

        TimeIntegrationSolution::TimeIntegrationSolution(const TimeIntegrationSchemeKey &key, 
                                                         const TripleArray& y, 
                                                         const Array<OneD, NekDouble>& t):
            m_scheme(TimeIntegrationSchemeManager()[key]),
            m_solVector(y),
            m_t(t)
        {
            ASSERTL1(y.num_elements()==m_scheme->GetNsteps(),"Amount of Entries does not match number of (multi-) steps");
        }

        TimeIntegrationSolution::TimeIntegrationSolution(const TimeIntegrationSchemeKey &key, 
                                                         unsigned int nvar,
                                                         unsigned int npoints):
            m_scheme(TimeIntegrationSchemeManager()[key]),
            m_solVector(m_scheme->GetNsteps()),
            m_t(m_scheme->GetNsteps())
        {
            for(int i = 0; i < m_scheme->GetNsteps(); i++)
            {
                m_solVector[i] = Array<OneD, Array<OneD, NekDouble> >(nvar);
                for(int j = 0; j < nvar; j++)
                {
                    m_solVector[i][j] = Array<OneD,NekDouble>(npoints);
                }
            }         
        }

        TimeIntegrationSolution::TimeIntegrationSolution(const TimeIntegrationSchemeKey &key):
            m_scheme(TimeIntegrationSchemeManager()[key]),
            m_solVector(m_scheme->GetNsteps()),
            m_t(m_scheme->GetNsteps())
        {      
        }
        
        TimeIntegrationSchemeSharedPtr TimeIntegrationScheme::Create(const TimeIntegrationSchemeKey &key)
        {
            TimeIntegrationSchemeSharedPtr returnval(
                MemoryManager<TimeIntegrationScheme>::AllocateSharedPtr(key));
            return returnval;
        }

        TimeIntegrationScheme::TimeIntegrationScheme(const TimeIntegrationSchemeKey &key):
            m_schemeKey(key),
            m_initialised(false)
        {
            switch(key.GetIntegrationMethod())
            {
            case eForwardEuler:
            case eAdamsBashforthOrder1:
                {
                    m_numsteps = 1;
                    m_numstages = 1;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,1.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps ,m_numsteps, 1.0);

                    m_schemeType = eExplicit;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            case eAdamsBashforthOrder2:
                {
                    m_numsteps = 2;
                    m_numstages = 1;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps,0.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps ,m_numsteps,0.0);

                    m_B[0][0][0] = 3.0/2.0;
                    m_B[0][1][0] = 1.0;

                    m_U[0][0] = 1.0;

                    m_V[0][0] = 1.0;
                    m_V[0][1] = -0.5;

                    m_schemeType = eExplicit;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 1;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                    m_timeLevelOffset[1] = 1;
                }
                break;
            case eAdamsBashforthOrder3:
                {
                    m_numsteps = 4;
                    m_numstages = 1;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps,0.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps ,m_numsteps,0.0);

                    m_B[0][1][0] = 1.0;

                    m_U[0][0] = 1.0;
                    m_U[0][1] = 23.0/12.0;
                    m_U[0][2] = -4.0/3.0;
                    m_U[0][3] = 5.0/12.0;

                    m_V[0][0] = 1.0;
                    m_V[0][1] = 23.0/12.0;
                    m_V[0][2] = -4.0/3.0;
                    m_V[0][3] = 5.0/12.0;
                    m_V[2][1] = 1.0;
                    m_V[3][2] = 1.0;

                    m_schemeType = eExplicit;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 3;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                    m_timeLevelOffset[1] = 1;
                    m_timeLevelOffset[2] = 2;
                    m_timeLevelOffset[3] = 3;
                }
                break;
            case eBackwardEuler:
            case eBDFImplicitOrder1:
            case eAdamsMoultonOrder1:
                {
                    m_numsteps = 1;
                    m_numstages = 1;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,1.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,1.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps ,m_numsteps, 1.0);

                    m_schemeType = eDiagonallyImplicit;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0; // SJS: Not sure whether this is correct
                }
                break;
            case eIMEXOrder1:
                {
                    m_numsteps  = 2;
                    m_numstages = 1;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,1.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,1.0);
                    m_A[1] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[1] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,1.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps ,m_numsteps, 0.0);
                    
                    m_B[0][1][0] = 0.0;
                    m_B[1][0][0] = 0.0;
                    m_V[0][0] = 1.0;
                    m_V[0][1] = 1.0;

                    m_schemeType = eIMEX;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 1;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                    m_timeLevelOffset[1] = 0;
                }
                break;
            case eIMEXOrder2:
                {
                    NekDouble third = 1.0/3.0;
                    m_numsteps  = 4;
                    m_numstages = 1;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,2*third);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,0.0);
                    m_A[1] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[1] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 4*third);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps ,m_numsteps, 0.0);
                    
                    m_B[0][0][0] = 2*third;
                    m_B[1][2][0] = 1.0;
                    m_U[0][1] = -third;
                    m_U[0][3] = -2*third;

                    m_V[0][0] =  4*third;
                    m_V[0][1] = -third;
                    m_V[0][2] =  4*third;
                    m_V[0][3] = -2*third;
                    m_V[1][0] =  1.0;
                    m_V[3][2] =  1.0;

                    m_schemeType = eIMEX;
                    m_numMultiStepValues = 2;
                    m_numMultiStepDerivs = 2;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                    m_timeLevelOffset[1] = 1;
                    m_timeLevelOffset[2] = 0;
                    m_timeLevelOffset[3] = 1;
                }
                break;
            case eIMEXOrder3:
                {
                    NekDouble eleventh = 1.0/11.0;
                    m_numsteps  = 6;
                    m_numstages = 1;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,6*eleventh);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,0.0);
                    m_A[1] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[1] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 18*eleventh);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps ,m_numsteps, 0.0);
                    
                    m_B[0][0][0] = 6*eleventh;
                    m_B[1][3][0] = 1.0;
                    m_U[0][1] = -9*eleventh;
                    m_U[0][2] =  2*eleventh;
                    m_U[0][4] = -18*eleventh;
                    m_U[0][5] =  6*eleventh;

                    m_V[0][0] =  18*eleventh;
                    m_V[0][1] = -9*eleventh;
                    m_V[0][2] =  2*eleventh;
                    m_V[0][3] =  18*eleventh;
                    m_V[0][4] = -18*eleventh;
                    m_V[0][5] =  6*eleventh;
                    m_V[1][0] =  1.0;
                    m_V[2][1] =  1.0;
                    m_V[4][3] =  1.0;
                    m_V[5][4] =  1.0;

                    m_schemeType = eIMEX;
                    m_numMultiStepValues = 3;
                    m_numMultiStepDerivs = 3;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                    m_timeLevelOffset[1] = 1;
                    m_timeLevelOffset[2] = 2;
                    m_timeLevelOffset[3] = 0;
                    m_timeLevelOffset[4] = 1;
                    m_timeLevelOffset[5] = 2;
                }
                break;
            case eAdamsMoultonOrder2:
                {
                    m_numsteps  = 2;
                    m_numstages = 1;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.5);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 0.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 0.0);

                    m_B[0][0][0] = 0.5;
                    m_B[0][1][0] = 1.0;

                    m_U[0][0] = 1.0;
                    m_U[0][1] = 0.5;

                    m_V[0][0] = 1.0;
                    m_V[0][1] = 0.5;


                    m_schemeType = eDiagonallyImplicit;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 1;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                    m_timeLevelOffset[1] = 0;
                }
                break;
            case eBDFImplicitOrder2:
                {
                    NekDouble third = 1.0/3.0;
                    m_numsteps  = 2;
                    m_numstages = 1;
                    
                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);
                    
                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,2*third);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 0.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 0.0);
                    
                    m_B[0][0][0] = 2*third;
                    m_B[0][1][0] = 0.0;
                    
                    m_U[0][0] = 4*third;
                    m_U[0][1] = -third;
                    
                    m_V[0][0] = 4*third;
                    m_V[0][1] = -third;
                    m_V[1][0] = 1.0;
                    
                    m_schemeType = eDiagonallyImplicit;
                    m_numMultiStepValues = 2;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                    m_timeLevelOffset[1] = 1;
                }
                break;
            case eMidpoint:
            case eRungeKutta2:
                {
                    m_numsteps = 1;
                    m_numstages = 2;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                    m_A[0][1][0] = 0.5;
                    m_B[0][0][1] = 1.0;

                    m_schemeType = eExplicit;
                    m_numMultiStepValues = 1; 
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            case eRungeKutta2_ImprovedEuler:
            case eRungeKutta2_SSP:
                {
                    m_numsteps = 1;
                    m_numstages = 2;
                    
                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);
                    
                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);
                    
                    m_A[0][1][0] = 1.0;
                    
                    m_B[0][0][0] = 0.5;
                    m_B[0][0][1] = 0.5;
                    
                    m_schemeType = eExplicit;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            case eRungeKutta3_SSP:
                {
                    m_numsteps = 1;
                    m_numstages = 3;
                    
                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);
                    
                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);
                    
                    m_A[0][1][0] = 1.0;
                    m_A[0][2][0] = 0.25;
                    m_A[0][2][1] = 0.25;
                    
                    m_B[0][0][0] = 1.0/6.0;
                    m_B[0][0][1] = 1.0/6.0;
                    m_B[0][0][2] = 2.0/3.0;
                    
                    m_schemeType = eExplicit;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            case eClassicalRungeKutta4:
            case eRungeKutta4:
                {
                    m_numsteps = 1;
                    m_numstages = 4;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                    m_A[0][1][0] = 0.5;
                    m_A[0][2][1] = 0.5;
                    m_A[0][3][2] = 1.0;

                    m_B[0][0][0] = 1.0/6.0;
                    m_B[0][0][1] = 1.0/3.0;
                    m_B[0][0][2] = 1.0/3.0;
                    m_B[0][0][3] = 1.0/6.0;

                    m_schemeType = eExplicit;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            case eDIRKOrder2:
                {
                    m_numsteps = 1;
                    m_numstages = 2;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                    NekDouble lambda = (2.0-sqrt(2.0))/2.0;

                    m_A[0][0][0] = lambda;
                    m_A[0][1][0] = 1.0 - lambda;
                    m_A[0][1][1] = lambda;

                    m_B[0][0][0] = 1.0 - lambda;
                    m_B[0][0][1] = lambda;

                    m_schemeType = eDiagonallyImplicit;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            case eDIRKOrder3:
                {
                    m_numsteps = 1;
                    m_numstages = 3;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                    NekDouble lambda = 0.4358665215;

                    m_A[0][0][0] = lambda;
                    m_A[0][1][0] = 0.5  * (1.0 - lambda);
                    m_A[0][2][0] = 0.25 * (-6.0*lambda*lambda + 16.0*lambda - 1.0);
                    m_A[0][1][1] = lambda;
                    m_A[0][2][1] = 0.25 * ( 6.0*lambda*lambda - 20.0*lambda + 5.0);
                    m_A[0][2][2] = lambda;

                    m_B[0][0][0] = 0.25 * (-6.0*lambda*lambda + 16.0*lambda - 1.0);
                    m_B[0][0][1] = 0.25 * ( 6.0*lambda*lambda - 20.0*lambda + 5.0);
                    m_B[0][0][2] = lambda;

                    m_schemeType = eDiagonallyImplicit;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
             case eDIRKOrder3Stage4Embedded5:
                {
                    // See "A FAMILY OF ESDIRK INTEGRATION METHODS" 
                    // ESDIRK 4/3 with 5 stages
                    //Or original paper "SINGLY DIAGONALLY IMPLICIT RUNGE–KUTTA METHODS WITH AN EXPLICIT FIRST STAGE"
                    m_numsteps = 1;
                    m_numstages = 4;
                     m_EmbeddedTemporalError=true;
                    //m_DirectTemporalError=true;
                    m_Order=3;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                    m_A[0][0][0] = 0.0;
                    m_A[0][1][0] = 0.43586652150846;
                    m_A[0][2][0] = 0.14073777472471;
                    m_A[0][3][0] = 0.10239940061991;

                    m_A[0][1][1] = 0.43586652150846;
                    m_A[0][2][1] = -0.10836555138132;
                    m_A[0][3][1] = -0.37687845225556;

                    m_A[0][2][2] = 0.43586652150846;
                    m_A[0][3][2] = 0.83861253012719;

                    m_A[0][3][3] = 0.43586652150846;

                    m_B[0][0][0] = m_A[0][3][0];
                    m_B[0][0][1] = m_A[0][3][1];
                    m_B[0][0][2] = m_A[0][3][2];
                    m_B[0][0][3] = m_A[0][3][3];

                    m_schemeType = eDiagonallyImplicit;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;

                    //Extra stage for estimating third order error
                    if(m_EmbeddedTemporalError)
                    {
                        m_A_Paired = Array<OneD, Array<TwoD,NekDouble> >(1);
                        m_B_Paired = Array<OneD, Array<TwoD,NekDouble> >(1);
                        //Extra one stage
                        m_A_Paired[0] = Array<TwoD,NekDouble>(1,m_numstages+1,0.0);
                        m_B_Paired[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages+1,0.0);
                        m_U_Paired    = Array<TwoD,NekDouble>(1,m_numsteps, 1.0);
                        //Unused 
                        m_V_Paired    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);
            
                        m_A_Paired[0][0][0] = 0.15702489786032;
                        m_A_Paired[0][0][1] = 0.11733044137044;
                        m_A_Paired[0][0][2] = 0.61667803039212;
                        m_A_Paired[0][0][3] = -0.32689989113134;
                        m_A_Paired[0][0][4] = 0.43586652150846;
                        m_B_Paired[0][0][0] = m_A_Paired[0][0][0];
                        m_B_Paired[0][0][1] = m_A_Paired[0][0][1];
                        m_B_Paired[0][0][2] = m_A_Paired[0][0][2];
                        m_B_Paired[0][0][3] = m_A_Paired[0][0][3];
                        m_B_Paired[0][0][4] = m_A_Paired[0][0][4];
                    }
                    

                    //To Do: Here just help debuging, remove
                    if(m_DirectTemporalError)
                    {
                        m_numsteps_Paired = 1;
                        m_numstages_Paired = 5;

                        m_A_Paired = Array<OneD, Array<TwoD,NekDouble> >(1);
                        m_B_Paired = Array<OneD, Array<TwoD,NekDouble> >(1);

                        m_A_Paired[0] = Array<TwoD,NekDouble>(m_numstages_Paired,m_numstages_Paired,0.0);
                        m_B_Paired[0] = Array<TwoD,NekDouble>(m_numsteps_Paired, m_numstages_Paired,0.0);
                        m_U_Paired    = Array<TwoD,NekDouble>(m_numstages_Paired,m_numsteps_Paired, 1.0);
                        m_V_Paired    = Array<TwoD,NekDouble>(m_numsteps_Paired, m_numsteps_Paired, 1.0);
                        
                        m_A_Paired[0][0][0] = 0.0;
                        m_A_Paired[0][1][0] = 0.43586652150846;
                        m_A_Paired[0][2][0] = 0.14073777472471;
                        m_A_Paired[0][3][0] = 0.10239940061991;
                        m_A_Paired[0][4][0] = 0.15702489786032;

                        m_A_Paired[0][1][1] = 0.43586652150846;
                        m_A_Paired[0][2][1] = -0.10836555138132;
                        m_A_Paired[0][3][1] = -0.37687845225556;
                        m_A_Paired[0][4][1] = 0.11733044137044;

                        m_A_Paired[0][2][2] = 0.43586652150846;
                        m_A_Paired[0][3][2] = 0.83861253012719;
                        m_A_Paired[0][4][2] = 0.61667803039212;

                        m_A_Paired[0][3][3] = 0.43586652150846;
                        m_A_Paired[0][4][3] = -0.32689989113134;

                        m_A_Paired[0][4][4] = 0.43586652150846;

                        m_B_Paired[0][0][0] = m_A_Paired[0][4][0];
                        m_B_Paired[0][0][1] = m_A_Paired[0][4][1];
                        m_B_Paired[0][0][2] = m_A_Paired[0][4][2];
                        m_B_Paired[0][0][3] = m_A_Paired[0][4][3];
                        m_B_Paired[0][0][4] = m_A_Paired[0][4][4];
                    }
                }
                break;
            case eDIRKOrder3Stage5:
                {
                    // See Kennedey&Carpenter 2016, Table 16
                    m_numsteps = 1;
                    m_numstages = 5;
                    m_DirectTemporalError=true;
                    m_Order=3;
                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                    NekDouble lambda = 9.0/40.0;
                    NekDouble ConstSqrt2  = sqrt(2.0);

                    m_A[0][0][0] = 0.0;
                    m_A[0][1][0] = lambda;
                    m_A[0][2][0] = 9.0*(1.0+ConstSqrt2)/80.0;
                    m_A[0][3][0] = (22.0+15.0*ConstSqrt2)/(80.0*(1+ConstSqrt2));
                    m_A[0][4][0] = (2398.0+1205.0*ConstSqrt2)/(2835.0*(4+3.0*ConstSqrt2));

                    m_A[0][1][1] = m_A[0][1][0];
                    m_A[0][2][1] = m_A[0][2][0];
                    m_A[0][3][1] = m_A[0][3][0];
                    m_A[0][4][1] = m_A[0][4][0];

                    m_A[0][2][2] = lambda;
                    m_A[0][3][2] = -7.0/(40.0*(1.0+ConstSqrt2));
                    m_A[0][4][2] = -2374*(2.0+ConstSqrt2)/(2835.0*(4.0+3.0*ConstSqrt2));

                    m_A[0][3][3] = lambda;
                    m_A[0][4][3] = 5827.0/7560.0;

                    m_A[0][4][4] = lambda;

                    m_B[0][0][0] = m_A[0][4][0];
                    m_B[0][0][1] = m_A[0][4][1];
                    m_B[0][0][2] = m_A[0][4][2];
                    m_B[0][0][3] = m_A[0][4][3];
                    m_B[0][0][4] = m_A[0][4][4];

                    m_schemeType = eDiagonallyImplicit;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;

                    //Paired Scheme move here DirkOrder4Stage6
                    if(m_DirectTemporalError)
                    {
                        m_numsteps_Paired = 1;
                        m_numstages_Paired = 6;

                        m_A_Paired = Array<OneD, Array<TwoD,NekDouble> >(1);
                        m_B_Paired = Array<OneD, Array<TwoD,NekDouble> >(1);

                        m_A_Paired[0] = Array<TwoD,NekDouble>(m_numstages_Paired,m_numstages_Paired,0.0);
                        m_B_Paired[0] = Array<TwoD,NekDouble>(m_numsteps_Paired, m_numstages_Paired,0.0);
                        m_U_Paired    = Array<TwoD,NekDouble>(m_numstages_Paired,m_numsteps_Paired, 1.0);
                        m_V_Paired    = Array<TwoD,NekDouble>(m_numsteps_Paired, m_numsteps_Paired, 1.0);

                        lambda = 1.0/4.0;
                        ConstSqrt2  = sqrt(2.0);

                        m_A_Paired[0][0][0] = 0.0;
                        m_A_Paired[0][1][0] = lambda;
                        m_A_Paired[0][2][0] = (1.0-ConstSqrt2)/8.0;
                        m_A_Paired[0][3][0] = (5.0-7.0*ConstSqrt2)/64.0;
                        m_A_Paired[0][4][0] = (-13796.0-54539*ConstSqrt2)/125000.0;
                        m_A_Paired[0][5][0] = (1181.0-987.0*ConstSqrt2)/ 13782.0;

                        m_A_Paired[0][1][1] = m_A_Paired[0][1][0];
                        m_A_Paired[0][2][1] = m_A_Paired[0][2][0];
                        m_A_Paired[0][3][1] = m_A_Paired[0][3][0];
                        m_A_Paired[0][4][1] = m_A_Paired[0][4][0];
                        m_A_Paired[0][5][1] = m_A_Paired[0][5][0];

                        m_A_Paired[0][2][2] = lambda;
                        m_A_Paired[0][3][2] = 7.0*(1.0+ConstSqrt2)/32.0;
                        m_A_Paired[0][4][2] = (506605.0+132109.0*ConstSqrt2)/437500.0;
                        m_A_Paired[0][5][2] = 47.0*(-267.0+1783.0*ConstSqrt2)/273343.0;

                        m_A_Paired[0][3][3] = lambda;
                        m_A_Paired[0][4][3] = 166.0*(-97.0+376.0*ConstSqrt2)/109375.0;
                        m_A_Paired[0][5][3] = -16.0*(-22922.0+3525.0*ConstSqrt2)/571953.0;

                        m_A_Paired[0][4][4] = lambda;
                        m_A_Paired[0][5][4] = -15625.0*(97.0+376.0*ConstSqrt2)/90749876.0;

                        m_A_Paired[0][5][5] = lambda;

                        m_B_Paired[0][0][0] = m_A_Paired[0][5][0];
                        m_B_Paired[0][0][1] = m_A_Paired[0][5][1];
                        m_B_Paired[0][0][2] = m_A_Paired[0][5][2];
                        m_B_Paired[0][0][3] = m_A_Paired[0][5][3];
                        m_B_Paired[0][0][4] = m_A_Paired[0][5][4];
                        m_B_Paired[0][0][5] = m_A_Paired[0][5][5];
                    }
                }
                break;
            case eDIRKOrder4Stage6:
                {
                    // See Kennedey&Carpenter 2016, Table 16
                    m_numsteps = 1;
                    m_numstages = 6;
                    m_Order=4;
                    m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                    NekDouble lambda = 1.0/4.0;
                    const NekDouble ConstSqrt2  = sqrt(2.0);

                    m_A[0][0][0] = 0.0;
                    m_A[0][1][0] = lambda;
                    m_A[0][2][0] = (1.0-ConstSqrt2)/8.0;
                    m_A[0][3][0] = (5.0-7.0*ConstSqrt2)/64.0;
                    m_A[0][4][0] = (-13796.0-54539*ConstSqrt2)/125000.0;
                    m_A[0][5][0] = (1181.0-987.0*ConstSqrt2)/ 13782.0;

                    m_A[0][1][1] = m_A[0][1][0];
                    m_A[0][2][1] = m_A[0][2][0];
                    m_A[0][3][1] = m_A[0][3][0];
                    m_A[0][4][1] = m_A[0][4][0];
                    m_A[0][5][1] = m_A[0][5][0];

                    m_A[0][2][2] = lambda;
                    m_A[0][3][2] = 7.0*(1.0+ConstSqrt2)/32.0;
                    m_A[0][4][2] = (506605.0+132109.0*ConstSqrt2)/437500.0;
                    m_A[0][5][2] = 47.0*(-267.0+1783.0*ConstSqrt2)/273343.0;

                    m_A[0][3][3] = lambda;
                    m_A[0][4][3] = 166.0*(-97.0+376.0*ConstSqrt2)/109375.0;
                    m_A[0][5][3] = -16.0*(-22922.0+3525.0*ConstSqrt2)/571953.0;

                    m_A[0][4][4] = lambda;
                    m_A[0][5][4] = -15625.0*(97.0+376.0*ConstSqrt2)/90749876.0;

                    m_A[0][5][5] = lambda;

                    m_B[0][0][0] = m_A[0][5][0];
                    m_B[0][0][1] = m_A[0][5][1];
                    m_B[0][0][2] = m_A[0][5][2];
                    m_B[0][0][3] = m_A[0][5][3];
                    m_B[0][0][4] = m_A[0][5][4];
                    m_B[0][0][5] = m_A[0][5][5];

                    m_schemeType = eDiagonallyImplicit;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            case eDIRKOrder4Stage6Embedded7:
            {
                // See "A FAMILY OF ESDIRK INTEGRATION METHODS" 
                // ESDIRK 5/4 with 7 stages
                //Or original paper "SINGLY DIAGONALLY IMPLICIT RUNGE–KUTTA METHODS WITH AN EXPLICIT FIRST STAGE"
                m_numsteps = 1;
                m_numstages = 6;
                m_EmbeddedTemporalError=true;
                //m_DirectTemporalError=true;
                m_Order=4;

                m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
                m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

                m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                m_A[0][0][0] = 0.0;
                m_A[0][1][0] =0.27000000000000000;
                m_A[0][2][0] = 0.13500000000000000;
                m_A[0][3][0] = 0.24814211234447322;
                m_A[0][4][0] = 0.25494479822150471;
                m_A[0][5][0] = 0.17549975523182941;

                m_A[0][1][1] = 0.27000000000000000;
                m_A[0][2][1] = 0.87265371804359686;
                m_A[0][3][1] = 0.13282088522859322;
                m_A[0][4][1] = 0.13106196422347200;
                m_A[0][5][1] = 0.0;

                m_A[0][2][2] = 0.27000000000000000;
                m_A[0][3][2] = -0.03886686658917771;
                m_A[0][4][2] = -0.04522093930235708;
                m_A[0][5][2] = -0.01641725931492383;

                m_A[0][3][3] = 0.27000000000000000;
                m_A[0][4][3] = 0.03389121682051642;
                m_A[0][5][3] = 3.59357175290010625;
                
                m_A[0][4][4] = 0.27000000000000000;
                m_A[0][5][4] = -3.02265424881701182;

                m_A[0][5][5] = 0.27000000000000000;


                m_B[0][0][0] = m_A[0][5][0];
                m_B[0][0][1] = m_A[0][5][1];
                m_B[0][0][2] = m_A[0][5][2];
                m_B[0][0][3] = m_A[0][5][3];
                m_B[0][0][4] = m_A[0][5][4];
                m_B[0][0][5] = m_A[0][5][5];

                m_schemeType = eDiagonallyImplicit;
                m_numMultiStepValues = 1;
                m_numMultiStepDerivs = 0;
                m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                m_timeLevelOffset[0] = 0;

                //Extra stage for estimating third order error
                if(m_EmbeddedTemporalError)
                {
                    m_A_Paired = Array<OneD, Array<TwoD,NekDouble> >(1);
                    m_B_Paired = Array<OneD, Array<TwoD,NekDouble> >(1);
                    //Extra one stage
                    m_A_Paired[0] = Array<TwoD,NekDouble>(1,m_numstages+1,0.0);
                    m_B_Paired[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages+1,0.0);
                    m_U_Paired    = Array<TwoD,NekDouble>(1,m_numsteps, 1.0);
                    //Unused 
                    m_V_Paired    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);
        
                    m_A_Paired[0][0][0] = 0.15847612643670410;
                    m_A_Paired[0][0][1] = 0.0;
                    m_A_Paired[0][0][2] = -0.07384703732094983;
                    m_A_Paired[0][0][3] =  5.26056776397634893;
                    m_A_Paired[0][0][4] = -4.83946947758407500;
                    m_A_Paired[0][0][5] = 0.22427262449197180;
                    m_A_Paired[0][0][6] = 0.27;
                    m_B_Paired[0][0][0] = m_A_Paired[0][0][0];
                    m_B_Paired[0][0][1] = m_A_Paired[0][0][1];
                    m_B_Paired[0][0][2] = m_A_Paired[0][0][2];
                    m_B_Paired[0][0][3] = m_A_Paired[0][0][3];
                    m_B_Paired[0][0][4] = m_A_Paired[0][0][4];
                    m_B_Paired[0][0][5] = m_A_Paired[0][0][5];
                    m_B_Paired[0][0][6] = m_A_Paired[0][0][6];
                }
            }
            break;
            case eIMEXdirk_2_3_2:
                {
                    m_numsteps  = 1;
                    m_numstages = 3;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_A[1] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[1] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                    NekDouble lambda = (2.0-sqrt(2.0))/2.0;
                    NekDouble delta = -2.0*sqrt(2.0)/3.0;

                    m_A[0][1][1] = lambda;
                    m_A[0][2][1] = 1.0 - lambda;
                    m_A[0][2][2] = lambda;

                    m_B[0][0][1] = 1.0 - lambda;
                    m_B[0][0][2] = lambda;

                    m_A[1][1][0] = lambda;
                    m_A[1][2][0] = delta;
                    m_A[1][2][1] = 1.0 - delta;

                    m_B[1][0][1] = 1.0 - lambda;
                    m_B[1][0][2] = lambda;

                    m_schemeType = eIMEX;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            case eIMEXdirk_3_4_3:
                {
                    m_numsteps  = 1;
                    m_numstages = 4;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_A[1] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[1] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                    NekDouble lambda = 0.4358665215;

                    m_A[0][1][1] = lambda;
                    m_A[0][2][1] = 0.5  * (1.0 - lambda);
                    m_A[0][3][1] = 0.25 * (-6.0*lambda*lambda + 16.0*lambda - 1.0);
                    m_A[0][2][2] = lambda;
                    m_A[0][3][2] = 0.25 * ( 6.0*lambda*lambda - 20.0*lambda + 5.0);
                    m_A[0][3][3] = lambda;

                    m_B[0][0][1] = 0.25 * (-6.0*lambda*lambda + 16.0*lambda - 1.0);
                    m_B[0][0][2] = 0.25 * ( 6.0*lambda*lambda - 20.0*lambda + 5.0);
                    m_B[0][0][3] = lambda;

                    m_A[1][1][0] = 0.4358665215;
                    m_A[1][2][0] = 0.3212788860;
                    m_A[1][2][1] = 0.3966543747;
                    m_A[1][3][0] =-0.105858296;
                    m_A[1][3][1] = 0.5529291479;
                    m_A[1][3][2] = 0.5529291479;

                    m_B[1][0][1] = 0.25 * (-6.0*lambda*lambda + 16.0*lambda - 1.0);
                    m_B[1][0][2] = 0.25 * ( 6.0*lambda*lambda - 20.0*lambda + 5.0);
                    m_B[1][0][3] = lambda;

                    m_schemeType = eIMEX;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            case eCNAB:
                {
                    NekDouble secondth = 1.0/2.0;
                    m_numsteps  = 4;
                    m_numstages = 1;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,secondth);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,0.0);
                    m_A[1] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[1] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps,secondth);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps ,m_numsteps, 0.0);

                    m_B[0][0][0] = secondth;
                    m_B[0][1][0] = 1.0;
                    m_B[1][2][0] = 1.0;
                    m_U[0][0] = 2*secondth;
                    m_U[0][2] = 3*secondth;
                    m_U[0][3] = -1*secondth;

                    m_V[0][0] = 2*secondth;
                    m_V[0][1] = secondth;
                    m_V[0][2] = 3*secondth;
                    m_V[0][3] = -1*secondth;
                    m_V[3][2] =  1.0;

                    m_schemeType = eIMEX;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 3;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                    m_timeLevelOffset[1] = 0;
                    m_timeLevelOffset[2] = 0;
                    m_timeLevelOffset[3] = 1;
                }
                break;
            case eIMEXGear:
                {
                    NekDouble twothirdth = 2.0/3.0;
		    		m_numsteps  = 2;
                    m_numstages = 1;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,twothirdth);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,0.0);
                    m_A[1] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[1] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps,twothirdth);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps ,m_numsteps, 0.0);

                    m_B[0][0][0] = twothirdth;
                    m_B[1][0][0] = twothirdth;
                    m_U[0][0] = 2*twothirdth;
                    m_U[0][1] = -0.5*twothirdth;

                    m_V[0][0] = 2*twothirdth;
                    m_V[0][1] = -0.5*twothirdth;
                    m_V[1][0] = 1.0;

                    m_schemeType = eIMEX;
                    m_numMultiStepValues = 2;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                    m_timeLevelOffset[1] = 1;
                }
                break;
            case eMCNAB:
                {
                    NekDouble sixthx = 9.0/16.0;
                    m_numsteps  = 5;
                    m_numstages = 1;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,sixthx);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,0.0);
                    m_A[1] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[1] = Array<TwoD,NekDouble>(m_numsteps ,m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 0.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps ,m_numsteps, 0.0);

                    m_B[0][0][0] = sixthx;
                    m_B[0][1][0] = 1.0;
                    m_B[1][3][0] = 1.0;
                    m_U[0][0] = 1.0;
                    m_U[0][1] = 6.0/16.0;
                    m_U[0][2] = 1.0/16.0;
                    m_U[0][3] = 1.5;
                    m_U[0][4] = -0.5;

                    m_V[0][0] = 1.0;
                    m_V[0][1] = 6.0/16.0;
                    m_V[0][2] = 1.0/16.0;
                    m_V[0][3] = 1.5;
                    m_V[0][4] = -0.5;
                    m_V[2][1] = 1.0;
                    m_V[4][3] = 1.0;

                    m_schemeType = eIMEX;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 4;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                    m_timeLevelOffset[1] = 0;
                    m_timeLevelOffset[2] = 1;
                    m_timeLevelOffset[3] = 0;
                    m_timeLevelOffset[4] = 1;
                }
                break;
            case eIMEXdirk_2_2_2:
                {
                    m_numsteps  = 1;
                    m_numstages = 3;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_A[1] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[1] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                    NekDouble glambda =  0.788675134594813;
                    NekDouble gdelta =  0.366025403784439;

                    m_A[0][1][1] = glambda;
                    m_A[0][2][1] = 1.0 - glambda;
                    m_A[0][2][2] = glambda;

                    m_B[0][0][1] = 1.0 - glambda;
                    m_B[0][0][2] = glambda;

                    m_A[1][1][0] = glambda;
                    m_A[1][2][0] = gdelta;
                    m_A[1][2][1] = 1.0 - gdelta;

                    m_B[1][0][0] = gdelta;
                    m_B[1][0][1] = 1.0 - gdelta;

                    m_schemeType = eIMEX;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            case eIMEXdirk_2_3_3:
                {
                    m_numsteps  = 1;
                    m_numstages = 3;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_A[1] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[1] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                    NekDouble glambda =  0.788675134594813;

                    m_A[0][1][1] = glambda;
                    m_A[0][2][1] = 1.0 - 2.0*glambda;
                    m_A[0][2][2] = glambda;

                    m_B[0][0][1] = 0.5;
                    m_B[0][0][2] = 0.5;

                    m_A[1][1][0] = glambda;
                    m_A[1][2][0] = glambda - 1.0;
                    m_A[1][2][1] = 2.0*(1-glambda);

                    m_B[1][0][1] = 0.5;
                    m_B[1][0][2] = 0.5;

                    m_schemeType = eIMEX;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            case eIMEXdirk_1_1_1:
                {
                    m_numsteps  = 1;
                    m_numstages = 2;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_A[1] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[1] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                    m_A[0][1][1] = 1.0;

                    m_B[0][0][1] = 1.0;

                    m_A[1][1][0] = 1.0;

                    m_B[1][0][0] = 1.0;

                    m_schemeType = eIMEX;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            case eIMEXdirk_1_2_1:
                {
                    m_numsteps  = 1;
                    m_numstages = 2;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_A[1] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[1] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);


                    m_A[0][1][1] = 1.0;

                    m_B[0][0][1] = 1.0;

                    m_A[1][1][0] = 1.0;

                    m_B[1][0][1] = 1.0;

                    m_schemeType = eIMEX;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            case eIMEXdirk_1_2_2:
                {
                    m_numsteps  = 1;
                    m_numstages = 2;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_A[1] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[1] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                    m_A[0][1][1] = 1.0/2.0;

                    m_B[0][0][1] = 1.0;

                    m_A[1][1][0] = 1.0/2.0;

                    m_B[1][0][1] = 1.0;

                    m_schemeType = eIMEX;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            case eIMEXdirk_4_4_3:
                {
                    m_numsteps  = 1;
                    m_numstages = 5;

                    m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
                    m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

                    m_A[0] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[0] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_A[1] = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B[1] = Array<TwoD,NekDouble>(m_numsteps, m_numstages,0.0);
                    m_U    = Array<TwoD,NekDouble>(m_numstages,m_numsteps, 1.0);
                    m_V    = Array<TwoD,NekDouble>(m_numsteps, m_numsteps, 1.0);

                    m_A[0][1][1] = 1.0/2.0;
                    m_A[0][2][1] = 1.0/6.0;
                    m_A[0][2][2] = 1.0/2.0;
                    m_A[0][3][1] = -1.0/2.0;
                    m_A[0][3][2] = 1.0/2.0;
                    m_A[0][3][3] = 1.0/2.0;
                    m_A[0][4][1] = 3.0/2.0;
                    m_A[0][4][2] = -3.0/2.0;
                    m_A[0][4][3] = 1.0/2.0;
                    m_A[0][4][4] = 1.0/2.0;

                    m_B[0][0][1] = 3.0/2.0;
                    m_B[0][0][2] = -3.0/2.0;
                    m_B[0][0][3] = 1.0/2.0;
                    m_B[0][0][4] = 1.0/2.0;

                    m_A[1][1][0] = 1.0/2.0;
                    m_A[1][2][0] = 11.0/18.0;
                    m_A[1][2][1] = 1.0/18.0;
                    m_A[1][3][0] = 5.0/6.0;
                    m_A[1][3][1] = -5.0/6.0;
                    m_A[1][3][2] = 1.0/2.0;
                    m_A[1][4][0] = 1.0/4.0;
                    m_A[1][4][1] = 7.0/4.0;
                    m_A[1][4][2] = 3.0/4.0;
                    m_A[1][4][3] = -7.0/4.0;

                    m_B[1][0][0] = 1.0/4.0;
                    m_B[1][0][1] = 7.0/4.0;
                    m_B[1][0][2] = 3.0/4.0;
                    m_B[1][0][3] = -7.0/4.0;

                    m_schemeType = eIMEX;
                    m_numMultiStepValues = 1;
                    m_numMultiStepDerivs = 0;
                    m_timeLevelOffset = Array<OneD,unsigned int>(m_numsteps);
                    m_timeLevelOffset[0] = 0;
                }
                break;
            default:
                {
                    NEKERROR(ErrorUtil::efatal,"Invalid Time Integration Scheme");
                }
            }

            m_firstStageEqualsOldSolution = CheckIfFirstStageEqualsOldSolution(m_A,m_B,m_U,m_V);
            m_lastStageEqualsNewSolution  = CheckIfLastStageEqualsNewSolution(m_A,m_B,m_U,m_V);

            ASSERTL1(VerifyIntegrationSchemeType(m_schemeType,m_A,m_B,m_U,m_V),
                     "Time integration scheme coefficients do not match its type");
        }


        bool TimeIntegrationScheme::
        VerifyIntegrationSchemeType(TimeIntegrationSchemeType type,
                                    const Array<OneD, const Array<TwoD, NekDouble> >& A,
                                    const Array<OneD, const Array<TwoD, NekDouble> >& B,
                                    const Array<TwoD, const NekDouble>& U,
                                    const Array<TwoD, const NekDouble>& V) const
        {
            int i;
            int j;
            int m;
            int  IMEXdim = A.num_elements();
            int  dim     = A[0].GetRows();

            Array<OneD, TimeIntegrationSchemeType> vertype(IMEXdim,eExplicit);

            for(m = 0; m < IMEXdim; m++)
            {
                for(i = 0; i < dim; i++)
                {
                    if( fabs(A[m][i][i]) > NekConstants::kNekZeroTol )
                    {
                        vertype[m] = eDiagonallyImplicit;
                    }
                }
                
                for(i = 0; i < dim; i++)
                {
                    for(j = i+1; j < dim; j++)
                    {
                        if( fabs(A[m][i][j]) > NekConstants::kNekZeroTol )
                        {
                            vertype[m] = eImplicit;
                            ASSERTL1(false,"Fully Implicit schemes cannnot be handled by the TimeIntegrationScheme class");
                        }
                    }
                }
            }

            if(IMEXdim == 2)
            {
                ASSERTL1(B.num_elements()==2,"Coefficient Matrix B should have an implicit and explicit part for IMEX schemes");
                if((vertype[0] == eDiagonallyImplicit) &&
                   (vertype[1] == eExplicit))
                {
                    vertype[0] = eIMEX;
                }
                else
                {
                    ASSERTL1(false,"This is not a proper IMEX scheme");
                }
            }

            return (vertype[0] == type);
        }

        TimeIntegrationSolutionSharedPtr 
        TimeIntegrationScheme::InitializeScheme(const NekDouble   timestep,
                                                ConstDoubleArray  &y_0    ,
                                                const NekDouble   time    ,
                                                const TimeIntegrationSchemeOperators &op)
        {
            // create a TimeIntegrationSolution object based upon the
            // initial value. Initialise all other multi-step values
            // and derivatives to zero
            TimeIntegrationSolutionSharedPtr y_out = 
                MemoryManager<TimeIntegrationSolution>::AllocateSharedPtr(m_schemeKey,y_0,time,timestep); 

            if( GetIntegrationSchemeType() == eExplicit)
            {
                // ensure initial solution is in correct space
                op.DoProjection(y_0,y_out->UpdateSolution(),time);
            }

            // calculate the initial derivative, if is part of the
            // solution vector of the current scheme
            if(m_numMultiStepDerivs)
            {
                if(m_timeLevelOffset[m_numMultiStepValues] == 0)
                {
                    int i;
                    int nvar    = y_0.num_elements();
                    int npoints = y_0[0].num_elements();
                    DoubleArray f_y_0(nvar);
                    for(i = 0; i < nvar; i++)
                    {
                        f_y_0[i] = Array<OneD,NekDouble>(npoints);
                    }
                    // calculate the derivative of the initial value
                    op.DoOdeRhs(y_0,f_y_0,time);
                    
                    // multiply by the step size
                    for(i = 0; i < nvar; i++)
                    {
                        Blas::Dscal(npoints,timestep,f_y_0[i].get(),1);
                    }
                    y_out->SetDerivative(0,f_y_0,timestep);
                }
            }
           
            return y_out;
        }
        
        TimeIntegrationScheme::ConstDoubleArray& 
        TimeIntegrationScheme::TimeIntegrate(const NekDouble    timestep, 
                                             TimeIntegrationSolutionSharedPtr &solvector,
                                             const TimeIntegrationSchemeOperators   &op)
        {
            ASSERTL1(!(GetIntegrationSchemeType() == eImplicit),
                     "Fully Implicit integration scheme cannot be handled by this routine.");
            
            int nvar    = solvector->GetFirstDim ();
            int npoints = solvector->GetSecondDim();

            if( (solvector->GetIntegrationScheme()).get() != this )
            {
                // This branch will be taken when the solution vector
                // (solvector) is set up for a different scheme than
                // the object this method is called from.  (typically
                // needed to calculate the first time-levels of a
                // multi-step scheme)

                // To do this kind of 'non-matching' integration, we
                // perform the following three steps:
                //
                // 1: copy the required input information from the
                //    solution vector of the master scheme to the
                //    input solution vector of the current scheme
                //
                // 2: time-integrate for one step using the current
                //    scheme
                //
                // 3: copy the information contained in the output
                //    vector of the current scheme to the solution
                //    vector of the master scheme

                // STEP 1: copy the required input information from
                //          the solution vector of the master scheme
                //          to the input solution vector of the
                //          current scheme

                // 1.1 Determine which information is required for the
                // current scheme
                int n;
                DoubleArray  y_n;
                NekDouble    t_n = 0;
                DoubleArray  dtFy_n;
                unsigned int nCurSchemeVals  = m_numMultiStepValues; // number of required values of the current scheme
                unsigned int nCurSchemeDers  = m_numMultiStepDerivs; // number of required derivs of the current scheme
                unsigned int nCurSchemeSteps = m_numsteps;  // number of steps in the current scheme
                unsigned int nMasterSchemeVals  = solvector->GetNvalues(); // number of values of the master scheme
                unsigned int nMasterSchemeDers  = solvector->GetNderivs(); // number of derivs of the master scheme
                // The arrays below contains information to which
                // time-level the values and derivatives of the
                // schemes belong
                const Array<OneD, const unsigned int>& curTimeLevels = m_timeLevelOffset; 
                const Array<OneD, const unsigned int>& masterTimeLevels = solvector->GetTimeLevelOffset(); 

                // 1.2 Copy the required information from the master
                //     solution vector to the input solution vector of
                //     the current scheme
                TimeIntegrationSolutionSharedPtr solvector_in = MemoryManager<TimeIntegrationSolution>::
                    AllocateSharedPtr(GetIntegrationSchemeKey()); // input solution vector of the current scheme

                for(n = 0; n < nCurSchemeVals; n++)
                {
                    // Get the required value out of the master solution vector
                    //DoubleArray& y_n = solvector->GetValue    ( curTimeLevels[n] );
                    //NekDouble    t_n = solvector->GetValueTime( curTimeLevels[n] );

                    y_n = solvector->GetValue    ( curTimeLevels[n] );
                    t_n = solvector->GetValueTime( curTimeLevels[n] );

                    // Set the required value in the input solution
                    // vector of the current scheme
                    solvector_in->SetValue(curTimeLevels[n],y_n,t_n);
                }
                for(n = nCurSchemeVals; n < nCurSchemeSteps; n++)
                {
                    // Get the required derivative out of the master
                    // solution vector
                    //DoubleArray& dtFy_n = solvector->GetDerivative    ( curTimeLevels[n] );
                    dtFy_n = solvector->GetDerivative    ( curTimeLevels[n] );

                    // Set the required derivative in the input
                    // solution vector of the current scheme
                    solvector_in->SetDerivative(curTimeLevels[n],dtFy_n,timestep);
                }

                // STEP 2: time-integrate for one step using the
                // current scheme
                TimeIntegrationSolutionSharedPtr solvector_out = MemoryManager<TimeIntegrationSolution>:: AllocateSharedPtr(GetIntegrationSchemeKey(),nvar,npoints);  // output solution vector of the current scheme

                // integrate
                TimeIntegrate(timestep, solvector_in->GetSolutionVector(),
                              solvector_in->GetTimeVector(),  
                              solvector_out->UpdateSolutionVector(),
                              solvector_out->UpdateTimeVector(),op);


                // STEP 3: copy the information contained in the
                //         output vector of the current scheme to the
                //         solution vector of the master scheme

                // 3.1 Check whether the current time scheme updates
                //     the most recent derivative that should be
                //     updated in the master scheme.  If not,
                //     calculate the derivative. This can be done
                //     based upon the corresponding value and the
                //     DoOdeRhs operator.
                int j;
                bool CalcNewDeriv = false; // flag inidicating whether the new derivative is availble in the output of
                                           // of the current scheme or whether it should be calculated
                if( nMasterSchemeDers > 0 )
                {
                    if(nCurSchemeDers == 0)
                    {
                        CalcNewDeriv = true;
                    }
                    else 
                    {
                        if( masterTimeLevels[nMasterSchemeVals] < curTimeLevels[nCurSchemeVals] )
                        {
                            CalcNewDeriv = true;
                        }
                    }
                }

                if(CalcNewDeriv)
                {
                    int newDerivTimeLevel = masterTimeLevels[nMasterSchemeVals]; // contains the time level at which
                                                                                 // we want to know the derivative of the
                                                                                 // master scheme
                    //DoubleArray  y_n;
                    //NekDouble    t_n;
                    // if the  time level correspond to 0, calculate the derivative based upon the solution value
                    // at the new time-level
                    if (newDerivTimeLevel == 0)
                    {
                        y_n = solvector_out->GetValue(0);
                        t_n = solvector_out->GetValueTime(0);
                    }
                    // if the  time level correspond to 1, calculate the derivative based upon the solution value
                    // at the new old-level
                    else if( newDerivTimeLevel == 1 ) 
                    {
                        y_n = solvector->GetValue(0);
                        t_n = solvector->GetValueTime(0);
                    }
                    else
                    {
                        ASSERTL1(false,"Problems with initialising scheme");
                    }
                    
                    DoubleArray  f_n(nvar);        
                    for(j = 0; j < nvar; j++)
                    {
                        f_n[j]   = Array<OneD, NekDouble>(npoints);
                    }
                    
                    // calculate the derivative
                    op.DoOdeRhs(y_n, f_n, t_n);
                    
                    // multiply by dt (as required by the General Linear Method framework)
                    for(j = 0; j < nvar; j++)
                    {
                        Vmath::Smul(npoints,timestep,f_n[j],1,
                                    f_n[j],1);
                    }
                    
                    // Rotate the solution vector 
                    // (i.e. updating without calculating/inserting new values)
                    solvector->RotateSolutionVector();
                    // Set the calculated derivative in the master solution vector
                    solvector->SetDerivative(newDerivTimeLevel,f_n,timestep);
                }
                else
                {
                    // Rotate the solution vector (i.e. updating
                    // without calculating/inserting new values)
                    solvector->RotateSolutionVector();
                }


                // 1.2 Copy the information calculated using the
                //     current scheme from the output solution vector
                //     to the master solution vector
                for(n = 0; n < nCurSchemeVals; n++)
                {
                    // Get the calculated value out of the output
                    // solution vector of the current scheme
                    //DoubleArray& y_n = solvector_out->GetValue    ( curTimeLevels[n] );
                    //NekDouble    t_n = solvector_out->GetValueTime( curTimeLevels[n] );
                    y_n = solvector_out->GetValue    ( curTimeLevels[n] );
                    t_n = solvector_out->GetValueTime( curTimeLevels[n] );

                    // Set the calculated value in the master solution vector
                    solvector->SetValue(curTimeLevels[n],y_n,t_n);
                }

                for(n = nCurSchemeVals; n < nCurSchemeSteps; n++)
                {
                    // Get the calculated derivative out of the output
                    // solution vector of the current scheme
                    // DoubleArray& dtFy_n =
                    // solvector_out->GetDerivative (curTimeLevels[n]);
                    dtFy_n = solvector_out->GetDerivative    ( curTimeLevels[n] );

                    // Set the calculated derivative in the master
                    // solution vector
                    solvector->SetDerivative(curTimeLevels[n],dtFy_n,timestep);
                }
            }
            else
            {
                const TimeIntegrationSchemeKey& key = solvector->GetIntegrationSchemeKey();
                
                TimeIntegrationSolutionSharedPtr solvector_new = MemoryManager<TimeIntegrationSolution>::AllocateSharedPtr(key,nvar,npoints); 
                /////////////////////////////////////////////////////////////////////////////////////
                //Direct Error calculating need back up initial solution
                if(m_TemporalErrorState && m_DirectTemporalError)
                {
                    //Back up initial solvector
                    //For temporary store original scheme
                    //FirstDim is variable number, you can check the codes of its definition
                    int nvar = GetFirstDim(solvector->GetSolutionVector());
                    //SecondDim is npoints, you can check the codes of its definition
                    int npoints= GetSecondDim(solvector->GetSolutionVector());
                    m_ttmp=Array<OneD,NekDouble> (m_numsteps,0.0);
                    m_solVectortmp=Array<OneD,Array<OneD, Array<OneD, NekDouble> > >(m_numsteps);
                    for(int i = 0; i < m_numsteps; i++)
                    {
                        //To Do: Because m_var have not defined yet
                        m_solVectortmp[i] = Array<OneD, Array<OneD, NekDouble> >(nvar);
                        for(int j = 0; j < nvar; j++)
                        {
                            m_solVectortmp[i][j] = Array<OneD,NekDouble>(npoints,0.0);
                        }
                    }
                    Vmath::Vcopy(solvector->GetTimeVector().num_elements(),solvector->GetTimeVector(),1,m_ttmp,1);
                    for(int i=0;i<m_solVectortmp.num_elements();i++)
                    {
                        for(int j=0;j<m_solVectortmp[i].num_elements();j++)
                        {
                            Vmath::Vcopy(npoints,solvector->GetSolutionVector()[i][j],1,m_solVectortmp[i][j],1);
                        }
                    }
                }
                ///////////////////////////////////////////////////////////////////////////////////////
               
                m_RealTimeStepFlag=true;
                TimeIntegrate(timestep,solvector->GetSolutionVector(),
                              solvector->GetTimeVector(),
                              solvector_new->UpdateSolutionVector(),
                              solvector_new->UpdateTimeVector(),op); 
                
                solvector = solvector_new;

                if(m_TemporalErrorState && m_DirectTemporalError)
                {
                    //Begin Calculating Error
                    //FirstDim is variable number, you can check the codes of its definition
                    m_nvar_Paired = GetFirstDim(solvector->GetSolutionVector());
                    //SecondDim is npoints, you can check the codes of its definition
                    m_npoints_Paired = GetSecondDim(solvector->GetSolutionVector());             
                    int nvar=m_nvar_Paired;
                    int npoints=m_npoints_Paired;

                    m_TemporalError = Array<OneD, Array<OneD, NekDouble> >(nvar);
                    for(int i = 0; i < nvar; i++)
                    {
                        m_TemporalError[i] = Array<OneD,NekDouble>(npoints,0.0);
                    }

                    m_RealTimeStepFlag=false;//This flag is used for tolerance 
                    PairedImplicitTimeIntegrate(timestep, 
                                        m_solVectortmp,
                                        m_ttmp,
                                        m_TemporalError,op);
                    
                    ASSERTL0(1==solvector->GetSolutionVector().num_elements(),"Only work for one-step DIRK scheme");
                    for(int i=0;i<m_TemporalError.num_elements();i++)
                    {
                        Vmath::Vsub(npoints,m_TemporalError[i],1,solvector->GetSolutionVector()[0][i],1,m_TemporalError[i],1);   
                        Vmath::Vabs(npoints,m_TemporalError[i],1,m_TemporalError[i],1);                    
                    }
                    m_TemporalErrorState=false;
                }               

            }
            return solvector->GetSolution();
        }
		
	    //To Do: Currently use Implicit DIRK schemes as Paired schemes. Cancel many flags from TimeIntegrate
       //Please check if the codes are general after testing
       void TimeIntegrationScheme::PairedImplicitTimeIntegrate(const NekDouble    timestep,
                                                  ConstTripleArray   &y_old  ,
                                                  ConstSingleArray   &t_old  ,
                                                  DoubleArray        &y_new  ,
                                                  const TimeIntegrationSchemeOperators &op)
        {
            unsigned int i,j,k;
            Array<OneD,NekDouble> t_new(t_old.num_elements(),0.0);

            ///////////////////////////////////////////////////////////////
            //Initialize tempory storage
            int nvar=m_nvar_Paired;
            int npoints= m_npoints_Paired;
            int nsteps=m_numsteps_Paired;
            int nstages=m_numstages_Paired;
            DoubleArray Y;      
            DoubleArray tmp;    
            TripleArray F;      
            NekDouble   T;  

            // Check if storage has already been initialised.
            // If so, we just zero the temporary storage.
            tmp = DoubleArray(nvar);
            for(j = 0; j < nvar; j++)
            {
                tmp[j]   = Array<OneD, NekDouble>(npoints,0.0);
            }

            Y = DoubleArray(nvar);
            for(j = 0; j < nvar; j++)
            {
                Y[j] =  Array<OneD, NekDouble>(npoints,0.0);
            }

            // Different storage for every stage derivative as the data
            // will be re-used to update the solution -> m_F is a TripleArray
            F = TripleArray(nstages);
            for(i = 0; i < nstages; ++i)
            {
                F[i]   = DoubleArray(nvar);
                for(j = 0; j < nvar; j++)
                {
                    F[i][j] = Array<OneD, NekDouble>(npoints,0.0);
                }
            }
         
			
            // The loop below calculates the stage values and derivatives
            for(i = 0; i < nstages; i++)
            {
                if( i==0 )
                {
                    for(k = 0; k < nvar; k++)
                    {
                        Vmath::Vcopy(npoints,y_old[0][k],1,Y[k],1);
                    }
                    T = t_old[0];
                }
                else
                {
                    // The stage values m_Y are a linear combination of:
                    // 1: the stage derivatives
					
                    if( i != 0 )
                    {
                        for(k = 0; k < nvar; k++)
                        {
                            Vmath::Smul(npoints,timestep*A_Paired(i,0),F[0][k],1,
                                        tmp[k],1);
                        }
                    }          
                    T = A_Paired(i,0)*timestep;
                        
                    for( j = 1; j < i; j++ )
                    {
                        for(k = 0; k < nvar; k++)
                        {
                            Vmath::Svtvp(npoints,timestep*A_Paired(i,j),F[j][k],1,
                                         tmp[k],1,tmp[k],1);
                        }          
                        
                        T += A_Paired(i,j)*timestep;
                    }
                    
                    // 2: the imported multi-step solution of the
                    // previous time level
                    for(k = 0; k < nvar; k++)
                    {
                        //U_Paired(i,0) is 1
                        Vmath::Svtvp(npoints,U_Paired(i,0),y_old[0][k],1,
                                        tmp[k],1,tmp[k],1);
                    }
                }
      
                // Calculate the stage derivative based upon the stage value
                if(i==0)
                {
                    // cout << "(i==0) && m_firstStageEqualsOldSolution "<<i<<endl;
                    op.DoProjection(Y,Y,T);
                    op.DoOdeRhs(Y, F[i], T);        
                }
                else
                {
                    if(nstages==1)
                    {
                        T= t_old[0]+timestep;
                    }
                    else 
                    {
                        T= t_old[0];
                        for(int j=0; j<=i; ++j)
                        {
                            T += A_Paired(i,j)*timestep;
                        }
                    }
                    
                    op.DoImplicitSolve(tmp, Y, T, A_Paired(i,i)*timestep);

                    for(k = 0; k < nvar; k++)
                    {
                        Vmath::Vsub(npoints,Y[k],1,tmp[k],1,F[i][k],1);
                        Vmath::Smul(npoints,1.0/(A_Paired(i,i)*timestep),F[i][k],1,F[i][k],1);
                    }
                }

            }

            // Next, the solution vector y at the new time level will
            // be calculated.
            //
            // For multi-step methods, this includes updating the
            // values of the auxiliary parameters
            //
            // The loop below calculates the solution at the new time
            // level
            //
            // If last stage equals the new solution, the new solution
            // needs not be calculated explicitly but can simply be
            // copied. This saves a solve. 	   	
            //Seems no need
            // for(k = 0; k < nvar; k++)
            // {
            //     Vmath::Vcopy(npoints,Y[k],1,y_new[k],1);
            // }
    

            // t_new[0] = B_Paired(0,0)*timestep;
            // for(j = 1; j < nstages; j++)
            // {
            //     t_new[0] += B_Paired(0,j)*timestep;
            // }

            		
            // The solution at the new time level is a linear
            // combination of: 
            // 1: the stage derivatives
            //Only one time step
            for(k = 0; k < nvar; k++)
            {
                Vmath::Smul(npoints,timestep*B_Paired(0,0),F[0][k],1,
                            y_new[k],1);
                
            }

            if(nstages != 1 )
            {
                t_new[0] = B_Paired(0,0)*timestep;
            }
            
    
            for(j = 1; j < nstages; j++)
            {
                for(k = 0; k < nvar; k++)
                {					
                    Vmath::Svtvp(npoints,timestep*B_Paired(0,j),F[j][k],1,
                                    y_new[k],1,y_new[k],1);
                }
                if(nstages != 1)
                {
                    t_new[0] += B_Paired(0,j)*timestep; 
                }
            }			
            
            // 2: the imported multi-step solution of the previous
            // time level
            for(k = 0; k < nvar; k++)
            {
                Vmath::Svtvp(npoints,V_Paired(0,0),y_old[0][k],1,
                                y_new[k],1,y_new[k],1);
            }
            if(nstages != 1)
            {
                t_new[0] += V_Paired(0,0)*t_old[0];
            }
        }

        void TimeIntegrationScheme::TimeIntegrate(const NekDouble    timestep,
                                                  ConstTripleArray   &y_old  ,
                                                  ConstSingleArray   &t_old  ,
                                                  TripleArray        &y_new  ,
                                                  SingleArray        &t_new  ,
                                                  const TimeIntegrationSchemeOperators &op)
        {
            ASSERTL1(CheckTimeIntegrateArguments(timestep,y_old,t_old,y_new,t_new,op), "Arguments not well defined");    
            
            unsigned int i,j,k;
            TimeIntegrationSchemeType type = GetIntegrationSchemeType();

            // Check if storage has already been initialised.
            // If so, we just zero the temporary storage.
            if (m_initialised && m_nvar == GetFirstDim(y_old)
                              && m_npoints == GetSecondDim(y_old))
            {
                for(j = 0; j < m_nvar; j++)
                {
                    Vmath::Zero(m_npoints, m_tmp[j], 1);
                }
            }
            else
            {
                m_nvar = GetFirstDim(y_old);
                m_npoints = GetSecondDim(y_old);
       
                // First, we are going to calculate the various stage
                // values and stage derivatives (this is the multi-stage
                // part of the method)
                // - m_Y   corresponds to the stage values
                // - m_F   corresponds to the stage derivatives
                // - m_T   corresponds to the time at the different stages
                // - m_tmp corresponds to the explicit right hand side of
                //   each stage equation
                //   (for explicit schemes, this correspond to m_Y)

                // Allocate memory for the arrays m_Y and m_F and m_tmp The same
                // storage will be used for every stage -> m_Y is a
                // DoubleArray
                m_tmp = DoubleArray(m_nvar);
                for(j = 0; j < m_nvar; j++)
                {
                    m_tmp[j]   = Array<OneD, NekDouble>(m_npoints,0.0);
                }

                // The same storage will be used for every stage -> m_tmp is
                // a DoubleArray
                if(type == eExplicit)
                {
                    m_Y = m_tmp;
                }
                else
                {
                    m_Y = DoubleArray(m_nvar);
                    for(j = 0; j < m_nvar; j++)
                    {
                        m_Y[j] =  Array<OneD, NekDouble>(m_npoints,0.0);
                    }
                }

                // Different storage for every stage derivative as the data
                // will be re-used to update the solution -> m_F is a TripleArray
                m_F = TripleArray(m_numstages);
                for(i = 0; i < m_numstages; ++i)
                {
                    m_F[i]   = DoubleArray(m_nvar);
                    for(j = 0; j < m_nvar; j++)
                    {
                        m_F[i][j] = Array<OneD, NekDouble>(m_npoints,0.0);
                    }
                }

                if(type == eIMEX)
                {
                    m_F_IMEX = TripleArray(m_numstages);
                    for(i = 0; i < m_numstages; ++i)
                    {
                        m_F_IMEX[i]   = DoubleArray(m_nvar);
                        for(j = 0; j < m_nvar; j++)
                        {
                            m_F_IMEX[i][j] =  Array<OneD, NekDouble>(m_npoints,0.0);
                        }
                    }
                }

                // Finally, flag that we have initialised the memory.
                m_initialised = true;
            }
			
            // The loop below calculates the stage values and derivatives
            for(i = 0; i < m_numstages; i++)
            {
                if( (i==0) && m_firstStageEqualsOldSolution )
                {
                    for(k = 0; k < m_nvar; k++)
                    {
                        Vmath::Vcopy(m_npoints,y_old[0][k],1,m_Y[k],1);
                    }
                    m_T = t_old[0];
                }
                else
                {
                    // The stage values m_Y are a linear combination of:
                    // 1: the stage derivatives
					
                    if( i != 0 )
                    {
                        for(k = 0; k < m_nvar; k++)
                        {
                            Vmath::Smul(m_npoints,timestep*A(i,0),m_F[0][k],1,
                                        m_tmp[k],1);
                            
                            if(type == eIMEX)       
                            {
                                Vmath::Svtvp(m_npoints,timestep*A_IMEX(i,0),
                                             m_F_IMEX[0][k],1,
                                             m_tmp[k],1,m_tmp[k],1);
                            }
                        }
                    }          
                    m_T = A(i,0)*timestep;
                        
                    for( j = 1; j < i; j++ )
                    {
                        for(k = 0; k < m_nvar; k++)
                        {
                            Vmath::Svtvp(m_npoints,timestep*A(i,j),m_F[j][k],1,
                                         m_tmp[k],1,m_tmp[k],1);
                            if(type == eIMEX)       
                            {
                                Vmath::Svtvp(m_npoints,timestep*A_IMEX(i,j),
                                             m_F_IMEX[j][k],1,
                                             m_tmp[k],1,m_tmp[k],1);
                            }
                        }          
                        
                        m_T += A(i,j)*timestep;
                    }
                    
                    // 2: the imported multi-step solution of the
                    // previous time level
                    for(j = 0; j < m_numsteps; j++)
                    {
                        for(k = 0; k < m_nvar; k++)
                        {
                            Vmath::Svtvp(m_npoints,U(i,j),y_old[j][k],1,
                                         m_tmp[k],1,m_tmp[k],1);
                        }
                        m_T += U(i,j)*t_old[j];
                    } 
                }
      
                // Calculate the stage derivative based upon the stage value
                if(type == eDiagonallyImplicit)
                {
                    if( (i==0) && m_firstStageEqualsOldSolution )
                    {
                        // cout << "(i==0) && m_firstStageEqualsOldSolution "<<i<<endl;
                        op.DoProjection(m_Y,m_Y,m_T);
                        op.DoOdeRhs(m_Y, m_F[i], m_T);      
                        // cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
                        // cout<<"Real Integrate"<<endl;  
                        // for(int m=0;m<m_nvar;m++)
                        // {    cout<<"Var["<<m<<"]"<<endl;
                        //     for(int n=0;n<m_npoints;n++)
                        //     {
                        //         cout<<std::setprecision(16)<<m_F[0][m][n]<<endl;
                        //     }
                        // }
                    }
                    else
                    {
                        if(m_numstages==1)
                        {
                            m_T= t_old[0]+timestep;
                        }
                        else 
                        {
                            m_T= t_old[0];
                            for(int j=0; j<=i; ++j)
                            {
                                m_T += A(i,j)*timestep;
                            }
                        }
                        
                        op.DoImplicitSolve(m_tmp, m_Y, m_T, A(i,i)*timestep);

                        for(k = 0; k < m_nvar; k++)
                        {
                            Vmath::Vsub(m_npoints,m_Y[k],1,m_tmp[k],1,m_F[i][k],1);
                            Vmath::Smul(m_npoints,1.0/(A(i,i)*timestep),m_F[i][k],1,m_F[i][k],1);
                        }
                    }
                }
                else if(type == eIMEX)
                { 
                    if(m_numstages==1)
                    {
                        m_T= t_old[0]+timestep;
                    }
                    else 
                    {
                        m_T= t_old[0];
                        for(int j=0; j<=i; ++j)
                        {
                            m_T += A(i,j)*timestep;
                        }
                    }	
                    
                    if(fabs(A(i,i)) > NekConstants::kNekZeroTol)
                    {
                        op.DoImplicitSolve(m_tmp, m_Y, m_T, A(i,i)*timestep);
                        
                        for(k = 0; k < m_nvar; k++)
                        {
                            Vmath::Vsub(m_npoints,m_Y[k],1,m_tmp[k],1,m_F[i][k],1);
                            Vmath::Smul(m_npoints,1.0/(A(i,i)*timestep),
                                        m_F[i][k],1,m_F[i][k],1);
                        }
                    }
                    op.DoOdeRhs(m_Y, m_F_IMEX[i], m_T);
                }
                else if( type == eExplicit)
                {
                    // Avoid projecting the same solution twice
                    if( ! ((i==0) && m_firstStageEqualsOldSolution) )
                    {
                        // ensure solution is in correct space
                        op.DoProjection(m_Y,m_Y,m_T);
                    }
                    op.DoOdeRhs(m_Y, m_F[i], m_T);        
                }
            }

            // Next, the solution vector y at the new time level will
            // be calculated.
            //
            // For multi-step methods, this includes updating the
            // values of the auxiliary parameters
            //
            // The loop below calculates the solution at the new time
            // level
            //
            // If last stage equals the new solution, the new solution
            // needs not be calculated explicitly but can simply be
            // copied. This saves a solve. 	   	
            int i_start = 0;
            if( (m_lastStageEqualsNewSolution) ) 
            {
                for(k = 0; k < m_nvar; k++)
                {
                    Vmath::Vcopy(m_npoints,m_Y[k],1,y_new[0][k],1);
                }
		
                if (m_numstages==1 && type == eIMEX)
                {
                    t_new[0] = t_old[0]+timestep;
                }
                else 
                {
                    t_new[0] = B(0,0)*timestep;
                    for(j = 1; j < m_numstages; j++)
                    {
                        t_new[0] += B(0,j)*timestep;
                    }
                    for(j = 0; j < m_numsteps; j++)
                    {
                        t_new[0] += V(0,j)*t_old[j];
                    }
                }
                i_start = 1;
            }
            
            for(i = i_start; i < m_numsteps; i++)
            {			
                // The solution at the new time level is a linear
                // combination of: 
                // 1: the stage derivatives
                for(k = 0; k < m_nvar; k++)
                {
                    Vmath::Smul(m_npoints,timestep*B(i,0),m_F[0][k],1,
                                y_new[i][k],1);
                    
                    if(type == eIMEX)
                    {
                        Vmath::Svtvp(m_npoints,timestep*B_IMEX(i,0),
                                     m_F_IMEX[0][k],1,y_new[i][k],1,
                                     y_new[i][k],1);
                    }
                }
                if(m_numstages != 1 || type != eIMEX)
                {
                    t_new[i] = B(i,0)*timestep;
                }
                
		
                for(j = 1; j < m_numstages; j++)
                {
                    for(k = 0; k < m_nvar; k++)
                    {					
                        Vmath::Svtvp(m_npoints,timestep*B(i,j),m_F[j][k],1,
                                     y_new[i][k],1,y_new[i][k],1);
                        
                        if(type == eIMEX)
                        {
                            Vmath::Svtvp(m_npoints,timestep*B_IMEX(i,j),
                                         m_F_IMEX[j][k],1,y_new[i][k],1,
                                         y_new[i][k],1);
                        }
                    }
                    if(m_numstages != 1 || type != eIMEX)
                    {
                        t_new[i] += B(i,j)*timestep; 
                    }
                }			
                
                // 2: the imported multi-step solution of the previous
                // time level
                for(j = 0; j < m_numsteps; j++)
                {
                    for(k = 0; k < m_nvar; k++)
                    {
                        Vmath::Svtvp(m_npoints,V(i,j),y_old[j][k],1,
                                     y_new[i][k],1,y_new[i][k],1);
                    }
                    if(m_numstages != 1 || type != eIMEX)
                    {
                        t_new[i] += V(i,j)*t_old[j];
                    }
                }
            }
            
            // Ensure that the new solution is projected if necessary
            if(type == eExplicit)
            {
                op.DoProjection(y_new[0],y_new[0],t_new[0]);
            }

            if(m_TemporalErrorState && m_EmbeddedTemporalError)
            {

                m_TemporalError = Array<OneD, Array<OneD, NekDouble> >(m_nvar);
                for(int i = 0; i < m_nvar; i++)
                {
                    m_TemporalError[i] = Array<OneD,NekDouble>(m_npoints,0.0);
                }


                Array<OneD,Array<OneD,NekDouble>> F_Paired(m_nvar);
                for(int k = 0; k < m_nvar; k++)
                {
                    F_Paired[k]= Array<OneD, NekDouble>(m_npoints,0.0);
                }
                
                //Calculate Extra stage
                for(int k = 0; k < m_nvar; k++)
                {
                    Vmath::Smul(m_npoints,timestep*A_Paired(0,0),m_F[0][k],1,
                                m_tmp[k],1);
                }
                m_T = A_Paired(0,0)*timestep;

                for( int j = 1; j < m_numstages; j++ )
                {
                    for(int k = 0; k < m_nvar; k++)
                    {
                        Vmath::Svtvp(m_npoints,timestep*A_Paired(0,j),m_F[j][k],1,
                                        m_tmp[k],1,m_tmp[k],1);
                    }          
                    
                    m_T += A_Paired(0,j)*timestep;
                }

                // 2: the imported multi-step solution of the
                // previous time level
                for(int k = 0; k < m_nvar; k++)
                {
                    //U_Paired(i,0) is 1
                    Vmath::Svtvp(m_npoints,U_Paired(0,0),y_old[0][k],1,
                                    m_tmp[k],1,m_tmp[k],1);
                }

                m_T= t_old[0];
                for(int j=0; j<=m_numstages; ++j)
                {
                    m_T += A_Paired(0,j)*timestep;
                }

                
                op.DoImplicitSolve(m_tmp, m_Y, m_T, A_Paired(0,m_numstages)*timestep);

                for(int k = 0; k < m_nvar; k++)
                {
                    Vmath::Vsub(m_npoints,m_Y[k],1,m_tmp[k],1,F_Paired[k],1);
                    Vmath::Smul(m_npoints,1.0/(A_Paired(0,m_numstages)*timestep),F_Paired[k],1,F_Paired[k],1);
                }

                for(int i= 0; i < m_numstages; i++)
                {
                    NekDouble diff=(m_B_Paired[0][0][i]-m_B[0][0][i])*timestep;
                    for(k = 0; k < m_nvar; k++)
                    {
                        Vmath::Svtvp(m_npoints,diff,m_F[i][k],1,m_TemporalError[k],1,m_TemporalError[k],1);
                    }
                }
                
                //Extra stage error
                for(k = 0; k < m_nvar; k++)
                {
                    NekDouble diff=m_B_Paired[0][0][m_numstages]*timestep;
                    Vmath::Svtvp(m_npoints,diff,F_Paired[k],1,m_TemporalError[k],1,m_TemporalError[k],1);
                    Vmath::Vabs(m_npoints,m_TemporalError[k],1,m_TemporalError[k],1);
                }
                
                m_TemporalErrorState=false;
            }
        }
        
        bool TimeIntegrationScheme::CheckIfFirstStageEqualsOldSolution(const Array<OneD, const Array<TwoD, NekDouble> >& A,
                                                                       const Array<OneD, const Array<TwoD, NekDouble> >& B,
                                                                       const Array<TwoD, const NekDouble>& U,
                                                                       const Array<TwoD, const NekDouble>& V) const
        {
            int i,m;
            // First stage equals old solution if:
            // 1. the first row of the coefficient matrix A consists of zeros
            // 2. U[0][0] is equal to one and all other first row entries of U are zero
 
            // 1. Check first condition
            for(m = 0; m < A.num_elements(); m++)
            {
                for(i = 0; i < m_numstages; i++)
                {
                    if( fabs(A[m][0][i]) > NekConstants::kNekZeroTol )
                    {
                        return false;
                    }
                }
            }

            // 2. Check second condition
            if( fabs(U[0][0] - 1.0) > NekConstants::kNekZeroTol )
            {
                return false;
            }
            for(i = 1; i < m_numsteps; i++)
            {
                if( fabs(U[0][i]) > NekConstants::kNekZeroTol )
                {
                    return false;
                }
            }

            return true;
        }
            
        bool TimeIntegrationScheme::CheckIfLastStageEqualsNewSolution(const Array<OneD, const Array<TwoD, NekDouble> >& A,
                                                                      const Array<OneD, const Array<TwoD, NekDouble> >& B,
                                                                      const Array<TwoD, const NekDouble>& U,
                                                                      const Array<TwoD, const NekDouble>& V) const
        {
            int i,m;
            // Last stage equals new solution if:
            // 1. the last row of the coefficient matrix A is equal to the first row of matrix B
            // 2. the last row of the coefficient matrix U is equal to the first row of matrix V
 
            // 1. Check first condition
            for(m = 0; m < A.num_elements(); m++)
            {
                for(i = 0; i < m_numstages; i++)
                {
                    if( fabs(A[m][m_numstages-1][i]-B[m][0][i]) > NekConstants::kNekZeroTol )
                    {
                        return false;
                    }
                }
            }

            // 2. Check second condition
            for(i = 0; i < m_numsteps; i++)
            {
                if( fabs(U[m_numstages-1][i]-V[0][i]) > NekConstants::kNekZeroTol )
                {
                    return false;
                }
            }

            return true;
        }
        
        bool TimeIntegrationScheme::CheckTimeIntegrateArguments(const NekDouble                      timestep,      
                                                                      ConstTripleArray               &y_old  ,
                                                                      ConstSingleArray               &t_old  ,
                                                                      TripleArray                    &y_new  ,
                                                                      SingleArray                    &t_new  ,
                                                                const TimeIntegrationSchemeOperators &op) const
        {
            // Check if arrays are all of consistent size
            ASSERTL1(y_old.num_elements()==m_numsteps,"Non-matching number of steps.");    
            ASSERTL1(y_new.num_elements()==m_numsteps,"Non-matching number of steps."); 
            
            ASSERTL1(y_old[0].   num_elements()==y_new[0].   num_elements(),"Non-matching number of variables.");  
            ASSERTL1(y_old[0][0].num_elements()==y_new[0][0].num_elements(),"Non-matching number of coefficients."); 
            
            ASSERTL1(t_old.num_elements()==m_numsteps,"Non-matching number of steps.");    
            ASSERTL1(t_new.num_elements()==m_numsteps,"Non-matching number of steps."); 

            return true;
        }
        
        std::ostream& operator<<(std::ostream& os, const TimeIntegrationSchemeSharedPtr& rhs)
        {
            return operator<<(os,*rhs);
        }

        std::ostream& operator<<(std::ostream& os, const TimeIntegrationScheme& rhs)
        {
            int i,j;
            int r = rhs.GetNsteps();
            int s = rhs.GetNstages();
            TimeIntegrationSchemeType type = rhs.GetIntegrationSchemeType();

            int oswidth = 9;
            int osprecision = 6;

            os << "Time Integration Scheme: " << TimeIntegrationMethodMap[rhs.GetIntegrationMethod()] << endl;
            os << "- number of steps:  " << r << endl;
            os << "- number of stages: " << s << endl;
            os << "- type of scheme:   " << TimeIntegrationSchemeTypeMap[rhs.GetIntegrationSchemeType()] << endl;
            os << "General linear method tableau: " << endl;

            for(i = 0; i < s; i++)
            {
                for(j = 0; j < s; j++)
                {
                    os.width(oswidth);
                    os.precision(osprecision);
                    os << right << rhs.A(i,j) << " ";
                }
                if(type == eIMEX)
                {
                    os << " '"; 
                    for(j = 0; j < s; j++)
                    {
                        os.width(oswidth);
                        os.precision(osprecision);
                        os << right << rhs.A_IMEX(i,j) << " ";
                    }
                }
                os << " |"; 

                for(j = 0; j < r; j++)
                {
                    os.width(oswidth);
                    os.precision(osprecision);
                    os << right << rhs.U(i,j);
                }
                os << endl;
            }
            int imexflag = (type == eIMEX)?2:1;
            for(int i = 0; i < (r+imexflag*s)*(oswidth+1)+imexflag*2-1; i++)
            {
                os << "-";
            }
            os << endl;
            for(i = 0; i < r; i++)
            {
                for(j = 0; j < s; j++)
                {
                    os.width(oswidth);
                    os.precision(osprecision);
                    os << right << rhs.B(i,j) << " ";
                }
                if(type == eIMEX)
                {
                    os << " '"; 
                    for(j = 0; j < s; j++)
                    {
                        os.width(oswidth);
                        os.precision(osprecision);
                        os << right << rhs.B_IMEX(i,j) << " ";
                    }
                }
                os << " |"; 

                for(j = 0; j < r; j++)
                {
                    os.width(oswidth);
                    os.precision(osprecision);
                    os << right << rhs.V(i,j);
                }
                os << endl;
            }
            return os;
        }
	}
}
