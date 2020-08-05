CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_VORTEXWAVE
    "Build the Vortex-Wave interaction solver." ON
    "NEKTAR_BUILD_SOLVERS;NEKTAR_SOLVER_INCNAVIERSTOKES" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS OR NOT NEKTAR_SOLVER_INCNAVIERSTOKES)
    SET(NEKTAR_SOLVER_VORTEXWAVE OFF CACHE INTERNAL "")
ENDIF()

IF( NEKTAR_SOLVER_VORTEXWAVE )
    SET(VortexWaveInteractionSolverSource
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusion.cpp
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusionReaction.cpp
       ./VortexWaveInteraction.cpp
       ./VortexWaveInteractionSolver.cpp
       $<TARGET_OBJECTS:IncNavierStokesSolver-obj>
       )

    SET(CalcL2ToLinfPressureSource
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusion.cpp
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusionReaction.cpp
       ./CalcL2ToLinfPressure.cpp
       ./VortexWaveInteraction.cpp
       $<TARGET_OBJECTS:IncNavierStokesSolver-obj>
       )

    SET(CalcVWIplusPressureSource
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusion.cpp
       ../ADRSolver/EquationSystems/SteadyAdvectionDiffusionReaction.cpp
       ./CalcVWIplusPressure.cpp
       ./VortexWaveInteraction.cpp
       $<TARGET_OBJECTS:IncNavierStokesSolver-obj>
       )

    SUBDIRS(utilities)

    ADD_SOLVER_EXECUTABLE(VortexWaveInteractionSolver SOURCES
	${VortexWaveInteractionSolverSource})
    TARGET_LINK_LIBRARIES(VortexWaveInteractionSolver PUBLIC SolverUtils)
    TARGET_LINK_LIBRARIES(VortexWaveInteractionSolver PUBLIC ADRSolver-plugin)
    ADD_SOLVER_UTILITY(CalcL2ToLinfPressure VortexWaveInteractionSolver SOURCES
	${CalcL2ToLinfPressureSource})
    TARGET_LINK_LIBRARIES(CalcL2ToLinfPressure PUBLIC SolverUtils)
    TARGET_LINK_LIBRARIES(CalcL2ToLinfPressure PUBLIC ADRSolver-plugin)
    ADD_SOLVER_UTILITY(CalcVWIplusPressure VortexWaveInteractionSolver SOURCES
	${CalcVWIplusPressureSource})
    TARGET_LINK_LIBRARIES(CalcVWIplusPressure PUBLIC SolverUtils)
    TARGET_LINK_LIBRARIES(CalcVWIplusPressure PUBLIC ADRSolver-plugin)
ENDIF( NEKTAR_SOLVER_VORTEXWAVE )
