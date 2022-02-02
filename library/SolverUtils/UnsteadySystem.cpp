///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadySystem.cpp
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
// Description: Generic timestepping for Unsteady solvers
//
///////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <iomanip>
using namespace std;

#include <boost/core/ignore_unused.hpp>
#include <boost/format.hpp>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/Forcing/Forcing.h>

namespace Nektar
{
    namespace SolverUtils
    {
        using namespace Vmath;
        // Debuggin function; prints values at various positions
        void print_some_array_values(NekInt N, NekInt inc, const Array<OneD, NekDouble> &arr, const string &str)
        {
            cout << endl << str << endl;
            for( int i = 0; i<N; i+=inc)
            {
                cout << arr[i] << endl;
            }
        }

        void print_array_values(NekInt N, const Array<OneD, NekDouble> &arr)
        {
            for( int i = 0; i<N; i++ )
            {
                // cout << arr[i] << "\t";
                printf("%.1f \t", arr[i]);
            }
            cout << endl;
        }

        // Copied from Extrapolate::RollOver
        void RollOver(Array<OneD, Array<OneD, NekDouble> > &input)
        {
            int  nlevels = input.size();

            Array<OneD, NekDouble> tmp;

            tmp = input[nlevels-1];

            for(int n = nlevels-1; n > 0; --n)
            {
                input[n] = input[n-1];
            }

            input[0] = tmp;
        }

        // From (Standard-)Extrapolate
        NekDouble StifflyStable_Betaq_Coeffs[3][3] = {
            { 1.0,  0.0, 0.0},{ 2.0, -1.0, 0.0},{ 3.0, -3.0, 1.0}
            };

        // 1D rollover
        void RollOverOneD(Array<OneD, NekDouble> &input)
        {
            int  nlevels = input.size();
            NekDouble tmp;
            tmp = input[nlevels-1];

            for(int n = nlevels-1; n > 0; --n)
            {
                input[n] = input[n-1];
            }

            input[0] = tmp;
        }

        void ExtrapolateArray(Array<OneD, NekDouble> &array, int nint)
        {
            // int nint     = 2
            int nlevels  = array.size();
            // int nPts     = array[0].size();

            // Update array
            cout << "before rollover" << endl;
            print_array_values(nlevels, array);
            RollOverOneD(array); // Move oldest from array[nlevels-1] to array[0]
            cout << "after rollover" << endl;
            print_array_values(nlevels, array);
            cout << endl;

            // Extrapolate to outarray
            // Smul(nPts, StifflyStable_Betaq_Coeffs[nint-1][nint-1], // beta_q
            //                 array[nint-1],    1, 
            //                 array[nlevels-1], 1);
            cout << "betaq  * value = " << StifflyStable_Betaq_Coeffs[nint-1][nint-1] << " * " << array[nint-1] << endl;
            array[nlevels-1] = StifflyStable_Betaq_Coeffs[nint-1][nint-1] * array[nint-1];

            for(int n = 0; n < nint-1; ++n)
            {
                cout << "betaq  * value = " << StifflyStable_Betaq_Coeffs[nint-1][n] << " * " << array[n] << endl;
                array[nlevels-1] = StifflyStable_Betaq_Coeffs[nint-1][n] * array[n] + array[nlevels-1];
                // Svtvp(nPts, StifflyStable_Betaq_Coeffs[nint-1][n],
                //             array[n],1, array[nlevels-1],1,
                //             array[nlevels-1],1);
            }
        }

        void test_ExtrapolateArray()
        {
            // Test extrapolate
            int sz = 2; // Defines order of BDF scheme
            Array<OneD, NekDouble> testarray(sz);
            for( int i = 0; i < sz; i++)
            {
                testarray[i] = i+1;
            }
            cout << "Input testarray:" << endl;
            print_array_values(sz, testarray);
            cout << endl;
            
            ExtrapolateArray(testarray, 2);
            cout << "After 1st ExtrapolateArray:" << endl;
            print_array_values(sz, testarray);
            cout << endl;
            
            testarray[sz-1] = 10;
            cout << "New input testarray:" << endl;
            print_array_values(sz, testarray);
            cout << endl;

            ExtrapolateArray(testarray, 2);
            cout << "After 2nd ExtrapolateArray:" << endl;
            print_array_values(sz, testarray);
            cout << endl;
            // exit(EXIT_SUCCESS); // Optionally
        }

        // Define mask, if none is defined in session file
        void GetUnmaskFunction(std::vector<std::vector<LibUtilities::EquationSharedPtr> > & unmaskfun, 
                               LibUtilities::SessionReaderSharedPtr m_session)
        {
            string Unmask0("Unmask0");
            string C0("C0");
            for(size_t i=0; 1; ++i)
            {
                Unmask0[Unmask0.size()-1] = '0' + i;
                if(!m_session->DefinesFunction(Unmask0))
                {
                    break;
                }
                for(size_t j=0; 1; ++j)
                {
                    C0[C0.size()-1] = '0' + j;
                    if(!m_session->DefinesFunction(Unmask0, C0))
                    {
                        break;
                    }
                    if(j==0)
                    {
                        unmaskfun.push_back(std::vector<LibUtilities::EquationSharedPtr>());
                    }
                    unmaskfun[unmaskfun.size()-1].push_back(m_session->GetFunction(Unmask0, C0));
                }
            }
        }

        void MaskInit(MultiRegions::ExpListSharedPtr &field, NekInt nfields, 
                      Array<OneD, NekDouble> &limits, 
                      Array<OneD, NekDouble> &maskCoeffs, 
                      Array<OneD, NekDouble> &maskPhys)
        {
            // std::vector<std::vector<LibUtilities::EquationSharedPtr> > unmaskfun;
            // GetUnmaskFunction(unmaskfun, session);
            // if(unmaskfun.size()==0)
            // {
            //     return;
            // }

            // Init mask fields
            int ncoef   = field->GetNcoeffs();
            int nphys   = field->GetNpoints();
            maskCoeffs = Array<OneD, NekDouble>(ncoef*nfields, 1.);
            maskPhys = Array<OneD, NekDouble>(nphys*nfields, 1.); // Masking for all fields individually

            // Get Boundary conditions and expansions
            Array<OneD, const SpatialDomains::BoundaryConditionShPtr> bndConds;
            Array<OneD, MultiRegions::ExpListSharedPtr> bndExp;
            bndConds = field->GetBndConditions();
            bndExp = field->GetBndCondExpansions();

            // Get center coordinates of each expansion touching the boundary
            // NekInt nbnds = bndExp.size();
            // NekInt nbexps = 0;
            // // NekDouble tol = 1e-8;
            // for(int i = 0; i < nbnds; i++)
            // {
            //     nbexps += bndExp[i]->GetExpSize();
            // }
            // Array<OneD, Array<OneD, NekDouble> > bndCoords(3);
            // for (int i = 0; i < 3; i++)
            // {
            //     bndCoords[i] = Array<OneD, NekDouble>(nbexps);
            // }

            // Loop each boundary
            // int cnt = 0;
            // for(int i = 0; i < nbnds; i++)
            // {
            //     int nbexp = bndExp[i]->GetExpSize();
            //     // cout << "bndExp[" << i << "]->GetExpSize()\t" << n << endl;
            //     // Loop expansions within boundary
            //     for(int j = 0; j < nbexp; j++)
            //     {
            //         LocalRegions::ExpansionSharedPtr exp = bndExp[i]->GetExp(j);
            //         // cout << "exp->GetElmtId()\t" << exp->GetElmtId() << endl; // ID on Bnd
            //         SpatialDomains::GeometrySharedPtr geom = exp->GetGeom();
            //         int nv = geom->GetNumVerts();
            //         NekDouble gcb[3] = {0.,0.,0.};
            //         NekDouble gct[3] = {0.,0.,0.};
            //         // Compute element's center coordinates
            //         for(size_t k=0; k<nv; ++k)
            //         {
            //             SpatialDomains::PointGeomSharedPtr vertex = geom->GetVertex(k);
            //             vertex->GetCoords(gct[0],gct[1],gct[2]);
            //             gcb[0] += gct[0]/NekDouble(nv);
            //             gcb[1] += gct[1]/NekDouble(nv);
            //             gcb[2] += gct[2]/NekDouble(nv);
            //         }
            //         bndCoords[0][cnt] = gcb[0];
            //         bndCoords[1][cnt] = gcb[1];
            //         bndCoords[2][cnt] = gcb[2];
            //         // cout << "Element_Coords: " << bndCoords[0][cnt] << ", " << bndCoords[1][cnt] << ", " << bndCoords[2][cnt] << endl;
            //         cnt += 1;
            //     }
            //     cout << endl << endl;
            // }

            // Loop over each element = Exp(i)
            for(size_t i=0; i<field->GetExpSize(); ++i)
            {
                LocalRegions::ExpansionSharedPtr exp = field->GetExp(i);
                SpatialDomains::GeometrySharedPtr geom = exp->GetGeom();
                int nv = geom->GetNumVerts();
                NekDouble gc[3] = {0.,0.,0.};
                NekDouble gct[3] = {0.,0.,0.};

                // Compute element's center coordinates
                for(size_t j=0; j<nv; ++j)
                {
                    SpatialDomains::PointGeomSharedPtr vertex = geom->GetVertex(j);
                    vertex->GetCoords(gct[0],gct[1],gct[2]);
                    gc[0] += gct[0]/NekDouble(nv);
                    gc[1] += gct[1]/NekDouble(nv);
                    gc[2] += gct[2]/NekDouble(nv);
                }
                
                // Masking logic
                int unmask = 1;
                
                // Mask by comparison to each boundary expansion
                // for (int k = 0; k < cnt; k++)
                // {
                //     cout << "gc[0] - bndCoords[0]: " << gc[0] << " - " << bndCoords[0][k] << " = " << gc[0] - bndCoords[0][k] << endl;
                //     cout << "gc[1] - bndCoords[1]: " << gc[1] << " - " << bndCoords[1][k] << " = " << gc[1] - bndCoords[1][k] << endl;
                //     cout << "gc[2] - bndCoords[2]: " << gc[2] << " - " << bndCoords[2][k] << " = " << gc[2] - bndCoords[2][k] << endl;
                //     if( (bndCoords[0][k] + tol > gc[0] && gc[0] > bndCoords[0][k] - tol) && 
                //         (bndCoords[1][k] + tol > gc[1] && gc[1] > bndCoords[1][k] - tol) && 
                //         (bndCoords[2][k] + tol > gc[2] && gc[2] > bndCoords[2][k] - tol) )
                //     {
                //         unmask = 0;
                //         cout << "unmask = 0" << endl;
                //         break;
                //     }
                //     cout << endl << endl;
                // }
                // if(unmask == 0)
                // {
                //     exit(EXIT_SUCCESS);
                // }
                
                // Mask by coordinate
                if( gc[0] < limits[0] ) { unmask = 0;} // x-
                if( gc[0] > limits[1] ) { unmask = 0;} // x+
                if( gc[1] < limits[2] ) { unmask = 0;} // y-
                if( gc[1] > limits[3] ) { unmask = 0;} // y+

                // If set mask? array[i] = 0 : array[i] = 1 (So, multiply with other fields to apply masking?)
                if(unmask==0)
                {
                    for(int j=0; j<nfields; ++j)
                    {
                        Vmath::Fill(exp->GetNcoeffs()  , 0., &maskCoeffs[field->GetCoeff_Offset(i) + j*ncoef], 1);
                        Vmath::Fill(exp->GetTotPoints(), 0., &maskPhys[field->GetPhys_Offset(i)  + j*nphys], 1);
                    }
                }
            }
        }



        /**
         * @class UnsteadySystem
         *
         * Provides the underlying timestepping framework for unsteady solvers
         * including the general timestepping routines. This class is not
         * intended to be directly instantiated, but rather is a base class
         * on which to define unsteady solvers.
         *
         * For details on implementing unsteady solvers see
         * \ref sectionADRSolverModuleImplementation here
         */

        /**
         * Processes SolverInfo parameters from the session file and sets up
         * timestepping-specific code.
         * @param   pSession        Session object to read parameters from.
         */
        UnsteadySystem::UnsteadySystem(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr& pGraph)
            : EquationSystem(pSession, pGraph),
              m_infosteps(10)

        {
        }

        /**
         * Initialization object for UnsteadySystem class.
         */
        void UnsteadySystem::v_InitObject()
        {
            EquationSystem::v_InitObject();

            m_initialStep = 0;

            // Load SolverInfo parameters
            m_session->MatchSolverInfo("DIFFUSIONADVANCEMENT","Explicit",
                                       m_explicitDiffusion,true);
            m_session->MatchSolverInfo("ADVECTIONADVANCEMENT","Explicit",
                                       m_explicitAdvection,true);
            m_session->MatchSolverInfo("REACTIONADVANCEMENT", "Explicit",
                                       m_explicitReaction, true);
            
            m_session->MatchSolverInfo("FLAGIMPLICITITSSTATISTICS", "True",
                                       m_flagImplicitItsStatistics, false);

            m_session->LoadParameter("CheckAbortSteps", m_abortSteps, 1);
            // Steady state tolerance
            m_session->LoadParameter("SteadyStateTol", m_steadyStateTol, 0.0);
            // Frequency for checking steady state
            m_session->LoadParameter("SteadyStateSteps",
                                          m_steadyStateSteps, 1);

            // For steady problems, we do not initialise the time integration
            if (m_session->DefinesSolverInfo("TimeIntegrationMethod") ||
                m_session->DefinesTimeIntScheme())
            {
                LibUtilities::TimeIntScheme timeInt;
                if (m_session->DefinesTimeIntScheme())
                {
                    timeInt = m_session->GetTimeIntScheme();
                }
                else
                {
                    timeInt.method = m_session->GetSolverInfo(
                        "TimeIntegrationMethod");
                }

                m_intScheme = LibUtilities::GetTimeIntegrationSchemeFactory()
                    .CreateInstance(timeInt.method,
                                    timeInt.variant,
                                    timeInt.order,
                                    timeInt.freeParams);

                // Load generic input parameters
                m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
                m_session->LoadParameter("IO_FiltersInfoSteps",
                    m_filtersInfosteps, 10.0 * m_infosteps);
                m_session->LoadParameter("CFL", m_cflSafetyFactor, 0.0);
                m_session->LoadParameter("CFLEnd", m_CFLEnd, 0.0);
                m_session->LoadParameter("CFLGrowth", m_CFLGrowth, 1.0);
                m_session->LoadParameter("CFLGrowth", m_CFLGrowth, 1.0);


                // Time tolerance between filter update time and time integration
                m_session->LoadParameter("FilterTimeWarning",
                                         m_filterTimeWarning, 1);

                // Ensure that there is no conflict of parameters
                if(m_cflSafetyFactor > 0.0)
                {
                    // Check final condition
                    ASSERTL0(m_fintime == 0.0 || m_steps == 0,
                             "Final condition not unique: "
                             "fintime > 0.0 and Nsteps > 0");
                    // Check timestep condition
                    ASSERTL0(m_timestep == 0.0,
                             "Timestep not unique: timestep > 0.0 & CFL > 0.0");
                }
                else
                {
                    ASSERTL0(m_timestep != 0.0,
                             "Need to set either TimeStep or CFL");
                }

                // Ensure that there is no conflict of parameters
                if (m_CFLGrowth > 1.0)
                {
                    // Check final condition
                    ASSERTL0(m_CFLEnd >= m_cflSafetyFactor,
                             "m_CFLEnd >= m_cflSafetyFactor required");
                }

                // Set up time to be dumped in field information
                m_fieldMetaDataMap["Time"] =
                        boost::lexical_cast<std::string>(m_time);
            }

            // By default attempt to forward transform initial condition.
            m_homoInitialFwd = true;

            // Set up filters
            for (auto &x : m_session->GetFilters())
            {
                m_filters.push_back(make_pair(x.first, GetFilterFactory().CreateInstance(
                        x.first, m_session, shared_from_this(), x.second)));
            }
        }

        /**
         * Destructor for the class UnsteadyAdvection.
         */
        UnsteadySystem::~UnsteadySystem()
        {
        }

        /**
         * @brief Returns the maximum time estimator for CFL control.
         */
        NekDouble UnsteadySystem::MaxTimeStepEstimator()
        {
            return m_intScheme->GetTimeStability();
        }

        /**
         * @brief Initialises the time integration scheme (as specified in the
         * session file), and perform the time integration.
         */
        void UnsteadySystem::v_DoSolve()
        {
            ASSERTL0(m_intScheme != 0, "No time integration scheme.");

            int i = 1;
            int nvariables = 0;
            int nfields = m_fields.size();

            if (m_intVariables.empty())
            {
                for (i = 0; i < nfields; ++i)
                {
                    m_intVariables.push_back(i);
                }
                nvariables = nfields;
            }
            else
            {
                nvariables = m_intVariables.size();
            }

            // Integrate in wave-space if using homogeneous1D
            if(m_HomogeneousType != eNotHomogeneous && m_homoInitialFwd)
            {
                for(i = 0; i < nfields; ++i)
                {
                    m_fields[i]->HomogeneousFwdTrans(m_fields[i]->GetPhys(),
                                                     m_fields[i]->UpdatePhys());
                    m_fields[i]->SetWaveSpace(true);
                    m_fields[i]->SetPhysState(false);
                }
            }

            // Set up wrapper to fields data storage.
            Array<OneD, Array<OneD, NekDouble> > fields(nvariables);

            // Order storage to list time-integrated fields first.
            for(i = 0; i < nvariables; ++i)
            {
                fields[i] = m_fields[m_intVariables[i]]->GetPhys();
                m_fields[m_intVariables[i]]->SetPhysState(false);
            }

            // Initialise time integration scheme
            m_intScheme->InitializeScheme( m_timestep, fields, m_time, m_ode );

            // Initialise filters
            for( auto &x : m_filters )
            {
                x.second->Initialise(m_fields, m_time);
            }

            LibUtilities::Timer         timer;
            bool      doCheckTime       = false;
            int       step              = m_initialStep;
            int       stepCounter       = 0;
            int       restartStep      = -1;
            NekDouble intTime           = 0.0;
            NekDouble cpuTime           = 0.0;
            NekDouble cpuPrevious       = 0.0;
            NekDouble elapsed           = 0.0;
            NekDouble totFilterTime     = 0.0;

            m_lastCheckTime             = 0.0;

            m_TotNewtonIts  = 0;
            m_TotLinIts   = 0;
            m_TotImpStages  = 0;

           Array<OneD, int> abortFlags(2, 0);
            string    abortFile     = "abort";
            if (m_session->DefinesSolverInfo("CheckAbortFile"))
            {
                abortFile = m_session->GetSolverInfo("CheckAbortFile");
            }

            NekDouble tmp_cflSafetyFactor = m_cflSafetyFactor;

            // Initiate arrays for operator residuals, used only with -v arg, TODO create member m_variables in UnsteadySystem.h
            NekInt nout = 5*2 + 2 + 4; // No of output fields (momentum equation 2D + residuals + velocity gradients 2D)
            Array<OneD, NekDouble> globalResc, globalResm, times;
            Array<OneD, Array<OneD, NekDouble> > gradient, tmp, outfields;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > tfields;
            if(m_verbose)
            {
                initialise_operator(nfields, nvariables, nout, globalResc, globalResm, times, gradient, tmp, outfields, tfields);
            }

            m_timestepMax = m_timestep;
            while ((step   < m_steps ||
                   m_time < m_fintime - NekConstants::kNekZeroTol) &&
                   abortFlags[1] == 0)
            {
                restartStep++;
                
                if(m_CFLGrowth > 1.0&&m_cflSafetyFactor<m_CFLEnd)
                {
                    tmp_cflSafetyFactor = 
                        min(m_CFLEnd,m_CFLGrowth*tmp_cflSafetyFactor);
                }

                m_flagUpdatePreconMat = true;

                // Flag to update AV
                m_CalcPhysicalAV = true;
                // Frozen preconditioner checks
                if (UpdateTimeStepCheck())
                {
                    m_cflSafetyFactor = tmp_cflSafetyFactor;

                    if (m_cflSafetyFactor)
                    {
                        m_timestep = GetTimeStep(fields);
                    }

                    // Ensure that the final timestep finishes at the final
                    // time, or at a prescribed IO_CheckTime.
                    if (m_time + m_timestep > m_fintime && m_fintime > 0.0)
                    {
                        m_timestep = m_fintime - m_time;
                    }
                    else if (m_checktime &&
                        m_time + m_timestep - m_lastCheckTime >= m_checktime)
                    {
                        m_lastCheckTime += m_checktime;
                        m_timestep     = m_lastCheckTime - m_time;
                        doCheckTime    = true;
                    }
                }

                if (m_TimeIncrementFactor > 1.0)
                {
                    NekDouble timeincrementFactor = m_TimeIncrementFactor;
                    m_timestep  *=  timeincrementFactor;

                    if (m_time + m_timestep > m_fintime && m_fintime > 0.0)
                    {
                        m_timestep = m_fintime - m_time;
                    }
                }

                // Perform any solver-specific pre-integration steps
                timer.Start();
                if (v_PreIntegrate(step))
                {
                    break;
                }

                m_StagesPerStep = 0;
                m_TotLinItePerStep = 0;

                ASSERTL0(m_timestep > 0, "m_timestep < 0");

                fields =
                    m_intScheme->TimeIntegrate( stepCounter, m_timestep, m_ode);
                timer.Stop();
                m_time  += m_timestep;
                elapsed  = timer.TimePerTest(1);
                intTime += elapsed;
                cpuTime += elapsed;

                // Write out status information
                if (m_session->GetComm()->GetRank() == 0 &&
                    !((step+1) % m_infosteps))
                {
                    cout << "Steps: " << setw(8)  << left << step+1 << " "
                         << "Time: "  << setw(12) << left << m_time;

                    if (m_cflSafetyFactor)
                    {
                        cout << " Time-step: " << setw(12)
                             << left << m_timestep;
                    }

                    stringstream ss;
                    ss << cpuTime << "s";
                    cout << " CPU Time: " << setw(8) << left
                         << ss.str() << endl;
                    cpuPrevious = cpuTime;
                    cpuTime = 0.0;

                    if(m_flagImplicitItsStatistics && m_flagImplicitSolver)
                    {
                        cout 
                             << "       &&" 
                             << " TotImpStages= " << m_TotImpStages 
                             << " TotNewtonIts= " << m_TotNewtonIts
                             << " TotLinearIts = " << m_TotLinIts  
                             << endl;
                    }
                }

                // Transform data into coefficient space
                for (i = 0; i < nvariables; ++i)
                {
                    // copy fields into ExpList::m_phys and assign the new
                    // array to fields
                    m_fields[m_intVariables[i]]->SetPhys(fields[i]);
                    fields[i] = m_fields[m_intVariables[i]]->UpdatePhys();
                    if( v_RequireFwdTrans() )
                    {
                        m_fields[m_intVariables[i]]->FwdTrans_IterPerExp(
                            fields[i],
                            m_fields[m_intVariables[i]]->UpdateCoeffs());
                    }
                    m_fields[m_intVariables[i]]->SetPhysState(false);
                }
                
                // Evaluate operator
                if(m_verbose)
                {
                    evaluate_operator(nfields, nvariables, stepCounter, nout, globalResc, globalResm, times, gradient, tmp, outfields, tfields);
                }

                // Perform any solver-specific post-integration steps
                if (v_PostIntegrate(step))
                {
                    break;
                }

                // Check for steady-state
                if (m_steadyStateTol > 0.0 && (!((step+1)%m_steadyStateSteps)) )
                {
                    if (CheckSteadyState(step,intTime))
                    {
                        if (m_comm->GetRank() == 0)
                        {
                            cout << "Reached Steady State to tolerance "
                                 << m_steadyStateTol << endl;
                        }
                        break;
                    }
                }

                // test for abort conditions (nan, or abort file)
                if (m_abortSteps && !((step+1) % m_abortSteps) )
                {
                    abortFlags[0] = 0;
                    for (i = 0; i < nvariables; ++i)
                    {
                        if (Vmath::Nnan(fields[i].size(),
                                fields[i], 1) > 0)
                        {
                            abortFlags[0] = 1;
                        }
                    }

                    //rank zero looks for abort file and deltes it
                    //if it exists. The communicates the abort
                    if(m_session->GetComm()->GetRank() == 0)
                    {
                        if(boost::filesystem::exists(abortFile))
                        {
                            boost::filesystem::remove(abortFile);
                            abortFlags[1] = 1;
                        }
                    }

                    m_session->GetComm()->AllReduce(abortFlags,
                                LibUtilities::ReduceMax);

                    ASSERTL0 (!abortFlags[0],
                                "NaN found during time integration.");
                }

                // Update filters
                for (auto &x : m_filters)
                {
                    timer.Start();
                    x.second->Update(m_fields, m_time);
                    timer.Stop();
                    elapsed = timer.TimePerTest(1);
                    totFilterTime += elapsed;

                    // Write out individual filter status information
                    if(m_session->GetComm()->GetRank() == 0 &&
                    !((step+1) % m_filtersInfosteps) && !m_filters.empty() &&
                    m_session->DefinesCmdLineArgument("verbose"))
                    {
                        stringstream s0;
                        s0 << x.first << ":";
                        stringstream s1;
                        s1 << elapsed << "s";
                        stringstream s2;
                        s2 << elapsed / cpuPrevious * 100 << "%";
                        cout << "CPU time for filter " << setw(25) << left
                            << s0.str() << setw(12) << left << s1.str() <<
                            endl << "\t Percentage of time integration:     "
                             << setw(10) << left << s2.str() << endl;
                    }
                }

                // Write out overall filter status information
                if (m_session->GetComm()->GetRank() == 0 &&
                    !((step+1) % m_filtersInfosteps) && !m_filters.empty())
                 {
                    stringstream ss;
                    ss << totFilterTime << "s";
                    cout << "Total filters CPU Time:\t\t\t     " << setw(10)
                        << left << ss.str() << endl;
                 }
                totFilterTime = 0.0;

                // Write out checkpoint files
                if ((m_checksteps && !((step + 1) % m_checksteps)) ||
                     doCheckTime)
                {
                    if(m_HomogeneousType != eNotHomogeneous)
                    {
                        vector<bool> transformed(nfields, false);
                        for(i = 0; i < nfields; i++)
                        {
                            if (m_fields[i]->GetWaveSpace())
                            {
                                m_fields[i]->SetWaveSpace(false);
                                m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                                      m_fields[i]->UpdatePhys());
                                m_fields[i]->SetPhysState(true);
                                transformed[i] = true;
                            }
                        }
                        Checkpoint_Output(m_nchk);
                        m_nchk++;
                        for(i = 0; i < nfields; i++)
                        {
                            if (transformed[i])
                            {
                                m_fields[i]->SetWaveSpace(true);
                                m_fields[i]->HomogeneousFwdTrans(
                                    m_fields[i]->GetPhys(),
                                    m_fields[i]->UpdatePhys());
                                m_fields[i]->SetPhysState(false);
                            }
                        }
                    }
                    else
                    {
                        Checkpoint_Output(m_nchk);
                        m_nchk++;
                    }
                    doCheckTime = false;

                    if(m_verbose)
                    {
                        // Write field of residual data to file
                        int ncoe = GetNcoeffs();
                        vector<Array<OneD, NekDouble> > coeffs(nout);
                        vector<string> variables(coeffs.size());
                        variables[0] = "resc";
                        variables[1] = "resm";
                        variables[2] = "advx";
                        variables[3] = "px";
                        variables[4] = "difx";
                        variables[5] = "utx";
                        variables[6] = "resmx";
                        variables[7] = "advy";
                        variables[8] = "py";
                        variables[9] = "dify";
                        variables[10] = "uty";
                        variables[11] = "resmy";
                        variables[12] = "dudx";
                        variables[13] = "dudy";
                        variables[14] = "dvdx";
                        variables[15] = "dvdy";
                        for (int l = 0; l < nout; l++)
                        {
                            coeffs[l] = Array<OneD, NekDouble>(ncoe);
                            m_fields[0]->FwdTrans_IterPerExp(outfields[l], coeffs[l]); // Which FwdTrans to use?
                        }
                        string filename = "res_" + boost::lexical_cast<string>(m_nchk - 1) + ".fld";
                        WriteFld(filename, m_fields[0], coeffs, variables);
                    }
                }

                // // Exit for debugging
                // if(step == 1)
                // {
                //     exit(EXIT_SUCCESS);
                // }

                // Step advance
                ++step;
                ++stepCounter;
                if(m_verbose)
                {
                    cout << "Step " << step << endl;
                }
            }


            // Write operator residuals to file
            if(m_verbose)
            {
                write_operator(stepCounter, globalResc, globalResm, times);
            }
            

            // Print out summary statistics
            if (m_session->GetComm()->GetRank() == 0)
            {
                if (m_cflSafetyFactor > 0.0)
                {
                    cout << "CFL safety factor : " << m_cflSafetyFactor << endl
                         << "CFL time-step     : " << m_timestep        << endl;
                }

                if (m_session->GetSolverInfo("Driver") != "SteadyState")
                {
                    cout << "Time-integration  : " << intTime  << "s"   << endl;
                }

                if(m_flagImplicitItsStatistics && m_flagImplicitSolver)
                {
                    cout 
                    << "-------------------------------------------" << endl
                    << "Total Implicit Stages: " << m_TotImpStages << endl
                    << "Total Newton Its     : " << m_TotNewtonIts << endl
                    << "Total Linear Its     : " << m_TotLinIts  << endl 
                    << "-------------------------------------------" << endl;
                }
            }

            // If homogeneous, transform back into physical space if necessary.
            if(m_HomogeneousType != eNotHomogeneous)
            {
                for(i = 0; i < nfields; i++)
                {
                    if (m_fields[i]->GetWaveSpace())
                    {
                        m_fields[i]->SetWaveSpace(false);
                        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                              m_fields[i]->UpdatePhys());
                        m_fields[i]->SetPhysState(true);
                    }
                }
            }
            else
            {
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[m_intVariables[i]]->SetPhys(fields[i]);
                    m_fields[m_intVariables[i]]->SetPhysState(true);
                }
            }

            // Finalise filters
            for (auto &x : m_filters)
            {
                x.second->Finalise(m_fields, m_time);
            }

            // Print for 1D problems
            if(m_spacedim == 1)
            {
                v_AppendOutput1D(fields);
            }
        }

        /**
         * @brief Sets the initial conditions.
         */
        void UnsteadySystem::v_DoInitialise()
        {
            CheckForRestartTime(m_time, m_nchk);
            SetBoundaryConditions(m_time);
            SetInitialConditions(m_time);
            InitializeSteadyState();
        }

        /**
         * @brief Prints a summary with some information regards the
         * time-stepping.
         */
        void UnsteadySystem::v_GenerateSummary(SummaryList& s)
        {
            EquationSystem::v_GenerateSummary(s);
            AddSummaryItem(s, "Advection",
                           m_explicitAdvection ? "explicit" : "implicit");

            if(m_session->DefinesSolverInfo("AdvectionType"))
            {
                AddSummaryItem(s, "AdvectionType",
                               m_session->GetSolverInfo("AdvectionType"));
            }

            AddSummaryItem(s, "Diffusion",
                           m_explicitDiffusion ? "explicit" : "implicit");

            if (m_session->GetSolverInfo("EQTYPE")
                    == "SteadyAdvectionDiffusionReaction")
            {
                AddSummaryItem(s, "Reaction",
                               m_explicitReaction  ? "explicit" : "implicit");
            }

            AddSummaryItem( s, "Time Step", m_timestep );
            AddSummaryItem( s, "No. of Steps", m_steps );
            AddSummaryItem( s, "Checkpoints (steps)", m_checksteps );
            AddSummaryItem( s, "Integration Type", m_intScheme->GetName() );
        }

        /**
         * Stores the solution in a file for 1D problems only. This method has
         * been implemented to facilitate the post-processing for 1D problems.
         */
        void UnsteadySystem::v_AppendOutput1D(
            Array<OneD, Array<OneD, NekDouble> > &solution1D)
        {
            // Coordinates of the quadrature points in the real physical space
            Array<OneD,NekDouble> x(GetNpoints());
            Array<OneD,NekDouble> y(GetNpoints());
            Array<OneD,NekDouble> z(GetNpoints());
            m_fields[0]->GetCoords(x, y, z);

            // Print out the solution in a txt file
            ofstream outfile;
            outfile.open("solution1D.txt");
            for(int i = 0; i < GetNpoints(); i++)
            {
                outfile << scientific << setw (17) << setprecision(16) << x[i]
                        << "  " << solution1D[0][i] << endl;
            }
            outfile << endl << endl;
            outfile.close();
        }

        void UnsteadySystem::CheckForRestartTime(NekDouble &time, int &nchk)
        {
            if (m_session->DefinesFunction("InitialConditions"))
            {
                for (int i = 0; i < m_fields.size(); ++i)
                {
                    LibUtilities::FunctionType vType;

                    vType = m_session->GetFunctionType(
                        "InitialConditions", m_session->GetVariable(i));

                    if (vType == LibUtilities::eFunctionTypeFile)
                    {
                        std::string filename
                            = m_session->GetFunctionFilename(
                                "InitialConditions", m_session->GetVariable(i));

                        fs::path pfilename(filename);

                        // redefine path for parallel file which is in directory
                        if(fs::is_directory(pfilename))
                        {
                            fs::path metafile("Info.xml");
                            fs::path fullpath = pfilename / metafile;
                            filename = LibUtilities::PortablePath(fullpath);
                        }
                        LibUtilities::FieldIOSharedPtr fld =
                            LibUtilities::FieldIO::CreateForFile(
                                m_session, filename);
                        fld->ImportFieldMetaData(filename, m_fieldMetaDataMap);

                        // check to see if time defined
                        if (m_fieldMetaDataMap !=
                                LibUtilities::NullFieldMetaDataMap)
                        {
                            auto iter = m_fieldMetaDataMap.find("Time");
                            if (iter != m_fieldMetaDataMap.end())
                            {
                                time = boost::lexical_cast<NekDouble>(
                                    iter->second);
                            }

                            iter = m_fieldMetaDataMap.find("ChkFileNum");
                            if (iter != m_fieldMetaDataMap.end())
                            {
                                nchk = boost::lexical_cast<NekDouble>(
                                    iter->second);
                            }
                        }

                        break;
                    }
                }
            }
        }

        /**
         * @brief Return the timestep to be used for the next step in the
         * time-marching loop.
         *
         * This function can be overloaded to facilitate solver which utilise a
         * CFL (or other) parameter to determine a maximum timestep under which
         * the problem remains stable.
         */
        NekDouble UnsteadySystem::GetTimeStep(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray)
        {
            return v_GetTimeStep(inarray);
        }

        /**
         * @brief Return the timestep to be used for the next step in the
         * time-marching loop.
         *
         * @see UnsteadySystem::GetTimeStep
         */
        NekDouble UnsteadySystem::v_GetTimeStep(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray)
        {
            boost::ignore_unused(inarray);
            NEKERROR(ErrorUtil::efatal, "Not defined for this class");
            return 0.0;
        }

        bool UnsteadySystem::v_PreIntegrate(int step)
        {
            boost::ignore_unused(step);
            return false;
        }

        bool UnsteadySystem::v_PostIntegrate(int step)
        {
            boost::ignore_unused(step);
            return false;
        }

        void UnsteadySystem::SVVVarDiffCoeff(
            const Array<OneD, Array<OneD, NekDouble> >  vel,
                  StdRegions::VarCoeffMap              &varCoeffMap)
        {
            int phystot = m_fields[0]->GetTotPoints();
            int nvel = vel.size();

            Array<OneD, NekDouble> varcoeff(phystot),tmp;

            // calculate magnitude of v
            Vmath::Vmul(phystot,vel[0],1,vel[0],1,varcoeff,1);
            for(int n = 1; n < nvel; ++n)
            {
                Vmath::Vvtvp(phystot,vel[n],1,vel[n],1,varcoeff,1,varcoeff,1);
            }
            Vmath::Vsqrt(phystot,varcoeff,1,varcoeff,1);

            for(int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
            {
                int offset = m_fields[0]->GetPhys_Offset(i);
                int nq = m_fields[0]->GetExp(i)->GetTotPoints();
                Array<OneD, NekDouble> unit(nq,1.0);

                int nmodes = 0;

                for(int n = 0; n < m_fields[0]->GetExp(i)->GetNumBases(); ++n)
                {
                    nmodes = max(nmodes,
                                 m_fields[0]->GetExp(i)->GetBasisNumModes(n));
                }

                NekDouble h = m_fields[0]->GetExp(i)->Integral(unit);
                h = pow(h,(NekDouble) (1.0/nvel))/((NekDouble) nmodes);

                Vmath::Smul(nq,h,varcoeff+offset,1,tmp = varcoeff+offset,1);
            }

            // set up map with eVarCoffLaplacian key
            varCoeffMap[StdRegions::eVarCoeffLaplacian] = varcoeff;
        }

        void UnsteadySystem::InitializeSteadyState()
        {
            if (m_steadyStateTol > 0.0)
            {
                const int nPoints = m_fields[0]->GetTotPoints();
                m_previousSolution = Array<OneD, Array<OneD, NekDouble> > (
                            m_fields.size());

                for (int i = 0; i < m_fields.size(); ++i)
                {
                    m_previousSolution[i] = Array<OneD, NekDouble>(nPoints);
                    Vmath::Vcopy(nPoints, m_fields[i]->GetPhys(), 1,
                                          m_previousSolution[i], 1);
                }

                if (m_comm->GetRank() == 0)
                {
                    std::string fName = m_session->GetSessionName() +
                        std::string(".resd");
                    m_errFile.open(fName.c_str());
                    m_errFile << setw(26) << left << "# Time";

                    m_errFile << setw(26) << left << "CPU_Time";

                    m_errFile << setw(26) << left << "Step";

                    for (int i = 0; i < m_fields.size(); ++i)
                    {
                       m_errFile << setw(26) << m_session->GetVariables()[i];
                    }

                    m_errFile << endl;
                }
            }
        }

        // // Virtual function to return forcing from IncNavierStokes.h
        // std::vector<SolverUtils::ForcingSharedPtr>* v_GetForcing()
        // {
        //     std::vector<SolverUtils::ForcingSharedPtr>* forcing;
        //     return forcing;
        // }

        /**
        * @brief Calculate whether the system has reached a steady state by
        * observing residuals to a user-defined tolerance.
        */
        bool UnsteadySystem::CheckSteadyState(int step)
        {
            return CheckSteadyState(step,0.0);
        }

        bool UnsteadySystem::CheckSteadyState(int step, NekDouble totCPUTime)
        {
            const int nPoints = GetTotPoints();
            const int nFields = m_fields.size();

            // Holds L2 errors.
            Array<OneD, NekDouble> L2       (nFields);

            SteadyStateResidual(step, L2);

            if (m_comm->GetRank() == 0 && ( ((step+1) % m_infosteps == 0)||((step== m_initialStep)) ))
            {
                // Output time
                m_errFile << boost::format("%25.19e") % m_time;

                m_errFile << " "<< boost::format("%25.19e") % totCPUTime;

                int stepp = step +1;

                m_errFile << " "<< boost::format("%25.19e") % stepp;

                // Output residuals
                for (int i = 0; i < nFields; ++i)
                {
                    m_errFile << " " << boost::format("%25.19e") % L2[i];
                }

                m_errFile << endl;
            }

            // Calculate maximum L2 error
            NekDouble maxL2 = Vmath::Vmax(nFields, L2, 1);

            if (m_session->DefinesCmdLineArgument("verbose") &&
                m_comm->GetRank() == 0 && ((step+1) % m_infosteps == 0))
            {
                cout << "-- Maximum L^2 residual: " << maxL2 << endl;
            }

            if (maxL2 <= m_steadyStateTol)
            {
                return true;
            }

            for (int i = 0; i < m_fields.size(); ++i)
            {
                Vmath::Vcopy(nPoints, m_fields[i]->GetPhys(), 1,
                                      m_previousSolution[i], 1);
            }

            m_steadyStateRes0 = m_steadyStateRes;
            m_steadyStateRes = maxL2;

            return false;
        }

        void UnsteadySystem::v_SteadyStateResidual(
            int                         step, 
            Array<OneD, NekDouble>      &L2)
        {
            boost::ignore_unused(step);
            const int nPoints = GetTotPoints();
            const int nFields = m_fields.size();

            // Holds L2 errors.
            Array<OneD, NekDouble> RHSL2    (nFields);
            Array<OneD, NekDouble> residual (nFields);
            Array<OneD, NekDouble> reference(nFields);

            for (int i = 0; i < nFields; ++i)
            {
                Array<OneD, NekDouble> tmp(nPoints);

                Vmath::Vsub(nPoints, m_fields[i]->GetPhys(), 1,
                                     m_previousSolution[i], 1, tmp, 1);
                Vmath::Vmul(nPoints, tmp, 1, tmp, 1, tmp, 1);
                residual[i] = Vmath::Vsum(nPoints, tmp, 1);

                Vmath::Vmul(nPoints, m_previousSolution[i], 1,
                                     m_previousSolution[i], 1, tmp, 1);
                reference[i] = Vmath::Vsum(nPoints, tmp, 1);
            }

            m_comm->AllReduce(residual , LibUtilities::ReduceSum);
            m_comm->AllReduce(reference, LibUtilities::ReduceSum);

            // L2 error
            for (int i = 0; i < nFields; ++i)
            {
                reference[i] = (reference[i] == 0) ? 1 : reference[i];
                L2[i] = sqrt(residual[i] / reference[i]);
            }
        }

        void UnsteadySystem::initialise_operator(int nfields, int nvariables, int nout,
                                                 Array<OneD, NekDouble> &globalResc,
                                                 Array<OneD, NekDouble> &globalResm,
                                                 Array<OneD, NekDouble> &times,
                                                 Array<OneD, Array<OneD, NekDouble> > &gradient,
                                                 Array<OneD, Array<OneD, NekDouble> > &tmp,
                                                 Array<OneD, Array<OneD, NekDouble> > &outfields,
                                                 Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &tfields)
        {
            // Definition of array sizes (excludes nfields and navariables)
            NekInt npoints = m_fields[0]->GetNpoints();
            NekInt gradsize = nfields * nvariables;
            NekInt nts = 2; // number of previous time fields to save

            // Initialise required storage arrays
            globalResc = Array<OneD, NekDouble>(m_steps);
            globalResm = Array<OneD, NekDouble>(m_steps);
            times = Array<OneD, NekDouble>(m_steps);

            // Initialise gradient array_
            gradient = Array<OneD, Array<OneD, NekDouble> >(gradsize);
            for (int i = 0; i < gradsize; i++)
            {
                gradient[i] = Array<OneD, NekDouble>(npoints, 0.0);
            }

            // Initialise temporary storage arrays
            tmp = Array<OneD, Array<OneD, NekDouble> >(nvariables);
            for (int i = 0; i < nvariables; i++)
            {
                tmp[i] = Array<OneD, NekDouble>(npoints, 0.0);
            }

            // Storage array of velocity fields for time derivative approximation
            tfields = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(nfields);
            for (int i = 0; i < nfields; ++i)
            {
                tfields[i] = Array<OneD, Array<OneD, NekDouble> >(nts);
                for(int n = 0; n < nts; ++n)
                {
                    tfields[i][n] = Array<OneD, NekDouble>(npoints, 0.0);

                    // Store initial conditions (n=0)
                    if(n == 0)
                    {
                        tfields[i][0] = m_fields[i]->GetPhys();
                    }
                }
            }
            
            // // Masking setup
            // Array<OneD, NekDouble> maskC; // Coefficient mask
            // Array<OneD, NekDouble> maskP; // Physical mask
            // bool MASKING = false;
            // if(MASKING)
            // {
            //     // Define masking coordinate limits +-x, +-y, +-z (TODO: move to session file)
            //     Array<OneD, NekDouble> clims(nvariables * 2);
            //     clims[0] = -20.; // x-
            //     clims[1] = 40; // x+
            //     clims[2] = -20.; // y-
            //     clims[3] = 20; // y+
            //     if( nvariables == 3 )
            //     {
            //         clims[4] = 1e-8; // z+
            //         clims[5] = -clims[4]; // z-
            //     }
            //     // Initialise masks
            //     MaskInit(m_fields[0], nfields, clims, maskC, maskP);
            // }

            // IF output fields of operator terms
            // Initialise output arrays
            outfields = Array<OneD, Array<OneD, NekDouble> >(nout);
            for(int i=0; i<nout; i++)
            {
                outfields[i] = Array<OneD, NekDouble>(npoints, 0.0);
            }
        }

        void UnsteadySystem::evaluate_operator(int nfields, int nvariables, int stepCounter, int nout,
                                               Array<OneD, NekDouble> &globalResc,
                                               Array<OneD, NekDouble> &globalResm,
                                               Array<OneD, NekDouble> &times,
                                               Array<OneD, Array<OneD, NekDouble> > &gradient,
                                               Array<OneD, Array<OneD, NekDouble> > &tmp,
                                               Array<OneD, Array<OneD, NekDouble> > &outfields,
                                               Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &tfields)
        {
            using namespace Vmath;

            // Low-level array size
            int npoints = m_fields[0]->GetNpoints();
            NekDouble kinvis = 0.0;
            
            // Read kinvis from session file
            if( m_session->DefinesParameter("Kinvis") )
            {
                kinvis = m_session->GetParameter("Kinvis");
            }
            else
            {
                ASSERTL0(0, "kinvis undefined; must be defined for operator");
            }

            // Temporary storage array
            int tmp_dims = 2;
            Array<OneD, Array<OneD, NekDouble> > temporary_array(tmp_dims);
            for (int i = 0; i < tmp_dims; i++)
            {
                temporary_array[i] = Array<OneD, NekDouble>(npoints, 0.0);
            }
            
            // Get new fields
            for( int i = 0; i<nfields; i++)
            {
                RollOver(tfields[i]);
                tfields[i][0] = m_fields[i]->GetPhys();
            }

            // Evaluate new derivatives
            for( int i = 0; i<nfields; i++)
            {
                
                if( nvariables == 2)
                {
                    m_fields[i]->PhysDeriv(m_fields[i]->GetPhys(),
                                        gradient[i*nvariables],
                                        gradient[i*nvariables+1]);
                }
                else if( nvariables == 3)
                {
                    m_fields[i]->PhysDeriv(m_fields[i]->GetPhys(),
                                        gradient[i*nvariables],
                                        gradient[i*nvariables+1],
                                        gradient[i*nvariables+2]);
                }
                else
                {
                    cout << "Cannot handle nvariables = " << nvariables << endl;
                    quick_exit(EXIT_SUCCESS);
                }
            }
            // Write gradient field out
            for(int i = 0; i < nvariables; i++)
            {
                Vcopy(npoints, gradient[i*nvariables], 1, outfields[nout-(2-i)*nvariables], 1);
                Vcopy(npoints, gradient[i*nvariables+1], 1, outfields[nout-(2-i)*nvariables+1], 1);
            }
            cout << "du/dx sum = " << Vsum(npoints, gradient[0], 1) << endl;
            cout << "du/dy sum = " << Vsum(npoints, gradient[1], 1) << endl;
            cout << "dv/dx sum = " << Vsum(npoints, gradient[2], 1) << endl;
            cout << "dv/dy sum = " << Vsum(npoints, gradient[3], 1) << endl;
            
            cout << "u sum = " << Vsum(npoints, m_fields[0]->GetPhys(), 1) << endl;
            cout << "v sum = " << Vsum(npoints, m_fields[1]->GetPhys(), 1) << endl;
            cout << "p sum = " << Vsum(npoints, m_fields[2]->GetPhys(), 1) << endl;

            // Zero and compute continuity residual
            Zero(npoints, temporary_array[0], 1);
            for( int i = 0; i < nvariables; i++ )
            {
                Vadd(npoints, temporary_array[0], 1, gradient[i*(nvariables+1)], 1, temporary_array[0], 1); // /nabla /cdot u
            }
            Vabs(npoints, temporary_array[0], 1, temporary_array[0], 1);
            
            // if(MASKING)
            // {
            //     Vmul(npoints, resc, 1, maskP, 1, resc, 1); // Apply mask
            // }  
            
            // If output residual field at each time step
            Vcopy(npoints, temporary_array[0], 1, outfields[0], 1);
            
            globalResc[stepCounter] = Vsum(npoints, temporary_array[0], 1);

            // Setup momentum residual
            for(int i=0; i<nvariables; i++)
            {
                Zero(npoints, temporary_array[i], 1);
                Zero(npoints, tmp[i], 1);
            }

            // Compute momentum residual
            for( int i = 0; i < nvariables; i++)
            {
                // Advection
                for( int k = 0; k < nvariables; k++ )
                {
                    Vvtvp(npoints, tfields[k][0], 1, gradient[i*nvariables+k], 1, tmp[i], 1, tmp[i], 1);
                }
                // Save x and y advection fields
                if( i == 0 )
                {                    
                    Vcopy(npoints, tmp[i], 1, outfields[0+2], 1); // u dx + v dy + w dz
                }
                else if( i == 1 )
                {
                    Vcopy(npoints, tmp[i], 1, outfields[0+2+5], 1);
                }
                
                // Pressure
                Vadd(npoints, gradient[nvariables*nvariables+i], 1, tmp[i], 1, tmp[i], 1); // + dP
                // Save x and y pressure gradient fields
                if( i == 0 )
                {
                    Vcopy(npoints, gradient[nvariables*nvariables+i], 1, outfields[1+2], 1);
                }
                else if( i == 1 )
                {
                    Vcopy(npoints, gradient[nvariables*nvariables+i], 1, outfields[1+2+5], 1);
                }

                // Diffusion
                for( int k = 0; k < nvariables; k++ )
                {
                    m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[k], gradient[i*nvariables+k], temporary_array[1]); // 2nd derivatives
                    Vadd(npoints, temporary_array[0], 1, temporary_array[1], 1, temporary_array[0], 1); // Sum to Laplacian
                }
                // Save x and y diffusion fields
                if( i == 0 )
                {
                    Smul(npoints, -kinvis, temporary_array[0], 1, outfields[2+2], 1);
                }
                else if( i == 1 )
                {
                    Smul(npoints, -kinvis, temporary_array[0], 1, outfields[2+2+5], 1);
                }
                Svtvp(npoints, -kinvis, temporary_array[0], 1, tmp[i], 1, tmp[i], 1); // add to operator balance
                

                // Time derivative
                Vsub(npoints, tfields[i][0], 1, tfields[i][1], 1, temporary_array[0], 1); // Vn - Vn-1
                Svtvp(npoints, 1.0/m_timestep, temporary_array[0], 1, tmp[i], 1, tmp[i], 1);
                // Save x and y time derivatives
                if( i == 0 )
                {
                    Smul(npoints, 1.0/m_timestep, temporary_array[0], 1, outfields[3+2], 1);
                }
                else if( i == 1 )
                {
                    Smul(npoints, 1.0/m_timestep, temporary_array[0], 1, outfields[3+2+5], 1);
                }

                // Add forcing to operator
                

                // Zero temporary array
                for(int i=0; i<tmp_dims; i++)
                {
                    Zero(npoints, temporary_array[i], 1);
                }
            }

            // Add momentum terms for output residual
            for(int i = 2; i < nout-5-4-1; i++)
            {
                Vadd(npoints, outfields[6], 1, outfields[i], 1, outfields[6], 1); // resm-X
            }
            for(int i = 2+5; i < nout-4-1; i++)
            {
                Vadd(npoints, outfields[11], 1, outfields[i], 1, outfields[11], 1); // resm-Y
            }                

            // L2 norm of momentum residuals
            for( int i = 0; i<nvariables; i++ )
            {
                Vvtvp(npoints, tmp[i], 1, tmp[i], 1, temporary_array[0], 1, temporary_array[0], 1);
            }
            Vsqrt(npoints, temporary_array[0], 1, temporary_array[0], 1);
            
            // if(MASKING) 
            // { 
            //     Vmul(npoints, resm, 1, maskP, 1, resm, 1); // Apply mask
            // }
            
            // Store for field output
            Vcopy(npoints, temporary_array[0], 1, outfields[1], 1);
            
            globalResm[stepCounter] = Vsum(npoints, temporary_array[0], 1);

            // Save time (for plotting residuals)
            times[stepCounter] = m_time;

            // Print operator residuals
            cout << "resc: " << globalResc[stepCounter] << "\tresm: " << globalResm[stepCounter] << endl;
        }

        void UnsteadySystem::write_operator(int stepCounter,
                            Array<OneD,NekDouble> &globalResc,
                            Array<OneD,NekDouble> &globalResm,
                            Array<OneD,NekDouble> &times)
        {
            ofstream opfile;
            opfile.open("nsResiduals.dat");
            opfile << "time,resc,resm" << endl;
            for( int i = 0; i < stepCounter; i++ )
            {
                opfile << times[i] << "," << globalResc[i] << "," << globalResm[i] << endl;
            }
            opfile.close();
        }
    }
}