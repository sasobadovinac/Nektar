# Main library sub-directories, required by all of Nektar++.
SUBDIRS(LibUtilities LocalRegions SpatialDomains StdRegions Collections
    MultiRegions MatrixFreeOps SolverUtils GlobalMapping FieldUtils NekMesh)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

IF (NEKTAR_BUILD_UNIT_TESTS)
    SUBDIRS(UnitTests)
    SUBDIRS(UnitTests/SIMD)
ENDIF()

IF (NEKTAR_BUILD_DEMOS)
    SUBDIRS(Demos)
ENDIF()

IF (NEKTAR_BUILD_PYTHON)
    INCLUDE_DIRECTORIES(Python)
ENDIF()

