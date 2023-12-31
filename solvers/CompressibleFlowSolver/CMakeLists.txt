CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_COMPRESSIBLE_FLOW
    "Build the Compressible Flow Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_COMPRESSIBLE_FLOW OFF CACHE INTERNAL "")
ENDIF()

IF( NEKTAR_SOLVER_COMPRESSIBLE_FLOW )
    SET(CompressibleFlowSolverSource
        ArtificialDiffusion/ArtificialDiffusion.cpp
        ArtificialDiffusion/NonSmoothShockCapture.cpp
        BoundaryConditions/CFSBndCond.cpp
        BoundaryConditions/ExtrapOrder0BC.cpp
        BoundaryConditions/IsentropicVortexBC.cpp
        BoundaryConditions/PressureInflowFileBC.cpp
        BoundaryConditions/PressureMachTemperatureBC.cpp
        BoundaryConditions/PressureOutflowBC.cpp
        BoundaryConditions/PressureOutflowNonReflectiveBC.cpp
        BoundaryConditions/RiemannInvariantBC.cpp
        BoundaryConditions/RinglebFlowBC.cpp
        BoundaryConditions/StagnationInflowBC.cpp
        BoundaryConditions/SymmetryBC.cpp
        BoundaryConditions/TimeDependentBC.cpp
        BoundaryConditions/WallBC.cpp
        BoundaryConditions/WallViscousBC.cpp
        CompressibleFlowSolver.cpp
        Diffusion/DiffusionLDGNS.cpp
        EquationSystems/CompressibleFlowSystem.cpp
        EquationSystems/CompressibleFlowSystemImplicit.cpp
        EquationSystems/EulerCFE.cpp
        EquationSystems/EulerImplicitCFE.cpp
        EquationSystems/NavierStokesCFE.cpp
        EquationSystems/NavierStokesCFEAxisym.cpp
        EquationSystems/NavierStokesImplicitCFE.cpp
        Forcing/ForcingAxiSymmetric.cpp
        Forcing/ForcingQuasi1D.cpp
       )

    ADD_SOLVER_EXECUTABLE(CompressibleFlowSolver
        SOURCES ${CompressibleFlowSolverSource}
        LIBRARY_SOURCES
        Misc/EquationOfState.cpp
        Misc/IdealGasEoS.cpp
        Misc/PengRobinsonEoS.cpp
        Misc/RedlichKwongEoS.cpp
        Misc/VanDerWaalsEoS.cpp
        Misc/VariableConverter.cpp
        RiemannSolvers/AverageSolver.cpp
        RiemannSolvers/AUSM0Solver.cpp
        RiemannSolvers/AUSM1Solver.cpp
        RiemannSolvers/AUSM2Solver.cpp
        RiemannSolvers/AUSM3Solver.cpp
        RiemannSolvers/CompressibleSolver.cpp
        RiemannSolvers/ExactSolverToro.cpp
        RiemannSolvers/HLLSolver.cpp
        RiemannSolvers/HLLCSolver.cpp
        RiemannSolvers/LaxFriedrichsSolver.cpp
        RiemannSolvers/RoeSolver.cpp
        RiemannSolvers/RoeSolverSIMD.cpp
        Preconditioner/PreconCfs.cpp
        Preconditioner/PreconCfsBRJ.cpp
        )
    IF(NEKTAR_BUILD_UNIT_TESTS)
      ADD_SUBDIRECTORY(RiemannSolvers/UnitTests)
    ENDIF()

    SUBDIRS(Utilities)

    ADD_NEKTAR_TEST(CylinderSubsonic_FRDG_SEM)
    ADD_NEKTAR_TEST(CylinderSubsonic_FRDG_GAUSS)
    ADD_NEKTAR_TEST(CylinderSubsonic_FRHU_SEM LENGTHY)
    ADD_NEKTAR_TEST(CylinderSubsonicMix)
    ADD_NEKTAR_TEST(CylinderSubsonic_P3)
    ADD_NEKTAR_TEST(CylinderSubsonic_P8 LENGTHY)
    ADD_NEKTAR_TEST(Euler1D)
    ADD_NEKTAR_TEST(IsentropicVortex16_P1)
    ADD_NEKTAR_TEST(IsentropicVortex16_P3)
    ADD_NEKTAR_TEST(IsentropicVortex_FRDG_SEM)
    ADD_NEKTAR_TEST(IsentropicVortex_FRSD_SEM LENGTHY)
    ADD_NEKTAR_TEST(IsentropicVortex_FRHU_SEM LENGTHY)
    ADD_NEKTAR_TEST(IsentropicVortex16_P8 LENGTHY)
    IF (NOT APPLE OR CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
        ADD_NEKTAR_TEST(IsentropicVortex16Implicit_P8 LENGTHY)
    ENDIF()
    ADD_NEKTAR_TEST(IsentropicVortex16_P3_GAUSS LENGTHY)
    ADD_NEKTAR_TEST(IsentropicVortex16_P8_GAUSS LENGTHY)
    ADD_NEKTAR_TEST(IsentropicVortex_WeakDG_HexDeformed)
    ADD_NEKTAR_TEST(RinglebFlow_P3)
    ADD_NEKTAR_TEST(RinglebFlow_P8 LENGTHY)
    #ADD_NEKTAR_TEST(Couette_WeakDG_LDG_MODIFIED)
    ADD_NEKTAR_TEST(Couette_WeakDG_IP_MODIFIED_IM)
    ADD_NEKTAR_TEST(Couette_WeakDG_IP_MODIFIED_IM_PRECON)
    ADD_NEKTAR_TEST(Couette_WeakDG_IP_MODIFIED_IM_PRECON_NonDim_Pri)
    ADD_NEKTAR_TEST(Couette_WeakDG_IP_MODIFIED_IM_PRECON_NonDim_Tet)
    ADD_NEKTAR_TEST(Couette_WeakDG_IP_MODIFIED_Adiabatic)
    ADD_NEKTAR_TEST(Couette_WeakDG_LDG_SEM)
    ADD_NEKTAR_TEST(Couette_WeakDG_LDG_SEM_VariableMu)
    #ADD_NEKTAR_TEST(Couette_WeakDG_LDG_GAUSS)
    #ADD_NEKTAR_TEST(Couette_FRDG_LFRDG_GAUSS)
    ADD_NEKTAR_TEST(Couette_FRDG_LDG_GAUSS LENGTHY)
    ADD_NEKTAR_TEST(Couette_FRHU_LDG_SEM LENGTHY)
    ADD_NEKTAR_TEST(Couette_WeakDG_LFRHU_SEM)
    ADD_NEKTAR_TEST(Couette_FRHU_LFRHU_SEM LENGTHY)
    ADD_NEKTAR_TEST(Couette_WeakDG_LDG_adiabatic)
    # disabled as fails on 32-bit linux  -- cc
    #ADD_NEKTAR_TEST(Couette_FRDG_LFRDG_adiabatic)
    # disabled as fails on 32-bit linux  -- cc
    #ADD_NEKTAR_TEST(CylinderSubsonic_NS_WeakDG_LDG_SEM LENGTHY)
    ADD_NEKTAR_TEST(CylinderSubsonic_NS_WeakDG_LDG_GAUSS LENGTHY)
    ADD_NEKTAR_TEST(CylinderSubsonic_NS_WeakDG_IP_GAUSS LENGTHY)
    ADD_NEKTAR_TEST(CylinderSubsonic_NS_FRDG_LFRDG_GAUSS LENGTHY)
    ADD_NEKTAR_TEST(CylinderSubsonic_NS_WeakDG_LDG_SEM_VariableMu)
    ADD_NEKTAR_TEST(CylinderSubsonic_NS_WeakDG_IP_SEM_VariableMu)
    ADD_NEKTAR_TEST(Couette_WeakDG_LDG_SEM_3DHOMO1D_MVM)
    ADD_NEKTAR_TEST(CylinderSubsonic_NS_WeakDG_LDG_SEM_3DHomo1D_MVM)
    ADD_NEKTAR_TEST(Couette_FRDG_LFRDG_GLL_LAGRANGE_3DHOMO1D_MVM)
    ADD_NEKTAR_TEST(CylinderSubsonic_NS_FRDG_LFRDG_GLL_LAGRANGE_3DHOMO1D_MVM)
    ADD_NEKTAR_TEST(Couette_FRHU_LFRHU_SEM_3DHOMO1D_MVM)
    ADD_NEKTAR_TEST(CylinderSubsonic_NS_FRHU_LFRHU_SEM_3DHOMO1D_MVM)
    ADD_NEKTAR_TEST(Couette_FRSD_LFRSD_MODIFIED_3DHOMO1D_MVM)
    ADD_NEKTAR_TEST(CylinderSubsonic_NS_FRSD_LFRSD_MODIFIED_3DHOMO1D_MVM)
    ADD_NEKTAR_TEST(CylinderSubsonic_WeakDG_Implicit)
    # disabled since it uses the old implementation of the StagnationInflowBC
    #ADD_NEKTAR_TEST(Nozzle_AxiSym_NoSwirl)
    # disabled since it uses the old implementation of the StagnationInflowBC
    #ADD_NEKTAR_TEST(Nozzle_AxiSym_Swirl)
    ADD_NEKTAR_TEST(implicitSolverCallsSFD LENGTHY)
    ADD_NEKTAR_TEST(Nozzle_Quasi1D_P6)
    ADD_NEKTAR_TEST(PipeFlow_NSAxisym)
    ADD_NEKTAR_TEST(hump3D_GLL)
    ADD_NEKTAR_TEST(hump3D_SEM)
    ADD_NEKTAR_TEST(hump3D_SEM_IP) 
    ADD_NEKTAR_TEST(Rarefaction_vanderWaals LENGTHY)
    ADD_NEKTAR_TEST(MMS_Compressible_Poiseuille_testIP)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Lap)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Lap_NS)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_AUSM0)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_AUSM1)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_AUSM2)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_AUSM3)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_ExactToro)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_HLL)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_HLLC)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_LaxFriedrichs)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_Roe)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Phys_Dil_LDGNS)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Phys_Mod_IP)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Phys_Mod_LDGNS)
    ADD_NEKTAR_TEST(ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_Impl)
    ADD_NEKTAR_TEST(TimeDependentAbsorption)
    ADD_NEKTAR_TEST(SquareDomain_Euler_2D_AxialFlow)
    ADD_NEKTAR_TEST(SquareDomain_Euler_2D_DiagonalFlow)
    ADD_NEKTAR_TEST(Duct_Euler_2D_AxiSym_ForcedVortexSwirl)
    ADD_NEKTAR_TEST(EmbededCylinder2D)
    ADD_NEKTAR_TEST(Movement_fixed_couette)

    IF (NEKTAR_USE_MPI)
        #ADD_NEKTAR_TEST(Perturbation_M05_square_CBC_par LENGTHY)
        #ADD_NEKTAR_TEST(Perturbation_M05_square_CBC_back_par LENGTHY)
        #ADD_NEKTAR_TEST(Perturbation_M15_square_CBC_par LENGTHY)
        #ADD_NEKTAR_TEST(Perturbation_M15_square_CBC_back_par LENGTHY)
        ADD_NEKTAR_TEST(Perturbation_M05_circle_CBC_back_par LENGTHY)
        ADD_NEKTAR_TEST(Perturbation_M15_circle_CBC_back_par LENGTHY)
        #ADD_NEKTAR_TEST(Perturbation_NS_M05_square_CBC_par LENGTHY)
        #ADD_NEKTAR_TEST(Perturbation_NS_M15_square_CBC_par LENGTHY)
        ADD_NEKTAR_TEST(RAE5240_BSF_LE_bcsfromfile_par)
        ADD_NEKTAR_TEST(Couette_WeakDG_IP_MODIFIED_IM_PRECON_par)
        # disabled as intermittently fails
        #ADD_NEKTAR_TEST(Couette_WeakDG_LDG_Exact_par)
        ADD_NEKTAR_TEST(Couette_WeakDG_IP_MODIFIED_IM_par)
        ADD_NEKTAR_TEST(Movement_fixed_couette_par)
    ENDIF (NEKTAR_USE_MPI)

    IF (NEKTAR_USE_FFTW)
	ADD_NEKTAR_TEST(Couette_WeakDG_LDG_SEM_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(CylinderSubsonic_NS_WeakDG_LDG_SEM_3DHomo1D_FFT)
	ADD_NEKTAR_TEST(IsentropicVortex16_WeakDG_SEM_3DHomo1D_FFT)
        ADD_NEKTAR_TEST(Couette_FRDG_LFRDG_GLL_LAGRANGE_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(CylinderSubsonic_NS_FRDG_LFRDG_GLL_LAGRANGE_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(IsentropicVortex16_FRDG_GLL_LAGRANGE_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(Couette_FRHU_LFRHU_SEM_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(CylinderSubsonic_NS_FRHU_LFRHU_SEM_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(IsentropicVortex16_FRHU_SEM_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(Couette_FRSD_LFRSD_MODIFIED_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(CylinderSubsonic_NS_FRSD_LFRSD_MODIFIED_3DHOMO1D_FFT)
	ADD_NEKTAR_TEST(IsentropicVortex16_FRSD_MODIFIED_3DHOMO1D_FFT)
    ENDIF (NEKTAR_USE_FFTW)
ENDIF( NEKTAR_SOLVER_COMPRESSIBLE_FLOW )
