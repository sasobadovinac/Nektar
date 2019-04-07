CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_GALERKINBOLTZMANN
    "Build the Galerkin Boltzmann solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_GALERKINBOLTZMANN OFF CACHE INTERNAL "")
ENDIF()

IF( NEKTAR_SOLVER_GALERKINBOLTZMANN )
    ADD_SOLVER_EXECUTABLE(GalerkinBoltzmannSolver
        SOURCES ./GalerkinBoltzmannSolver.cpp
        LIBRARY_SOURCES
        ./EquationSystems/GalerkinBoltzmann.cpp
        ./RiemannSolvers/BoltzmannSolver.cpp
        ./RiemannSolvers/UpwindSolver.cpp
        )       

ENDIF( NEKTAR_SOLVER_GALERKINBOLTZMANN )
