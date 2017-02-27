CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_OPTIONPRICING
    "Build the Option Pricing Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_OPTIONPRICING )
    SET(OptionPricingSolverSource 
       ./OptionPricingSolver.cpp
       ./EquationSystems/BS.cpp)

    ADD_SOLVER_EXECUTABLE(OptionPricingSolver solvers 
        ${OptionPricingSolverSource})

ENDIF( NEKTAR_SOLVER_OPTIONPRICING )
