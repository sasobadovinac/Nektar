SET(ADR_SOLVER_EQUATION_SYSTEMS_SOURCES
    Helmholtz.cpp
    Laplace.cpp
    Poisson.cpp
    Projection.cpp
    SteadyAdvectionDiffusion.cpp
    SteadyAdvectionDiffusionReaction.cpp
    UnsteadyAdvection.cpp
    MMFAdvection.cpp
    UnsteadyDiffusion.cpp
    UnsteadyAdvectionDiffusion.cpp
    UnsteadyInviscidBurger.cpp
    UnsteadyReactionDiffusion.cpp
    UnsteadyViscousBurgers.cpp
    EigenValuesAdvection.cpp
)

SET(ADR_SOLVER_EQUATION_SYSTEMS_HEADERS
    Helmholtz.h
    Laplace.h
    Poisson.h
    Projection.h
    SteadyAdvectionDiffusion.h
    SteadyAdvectionDiffusionReaction.h
    UnsteadyAdvection.h
    MMFAdvection.h
    UnsteadyDiffusion.h
    UnsteadyAdvectionDiffusion.h
    UnsteadyInviscidBurger.h
    UnsteadyReactionDiffusion.h
    UnsteadyViscousBurgers.h
    EigenValuesAdvection.h
)


ADD_LIBRARY(ADRSolverEquationSystems SHARED ${ADR_SOLVER_EQUATION_SYSTEMS_SOURCES})
ADD_DEPENDENCIES(ADRSolverEquationSystems SolverUtils SpatialDomains)

TARGET_LINK_LIBRARIES(ADRSolverEquationSystems
  PUBLIC SolverUtils
  PUBLIC SpatialDomains
)
