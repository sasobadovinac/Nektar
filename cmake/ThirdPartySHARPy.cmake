########################################################################
#
# ThirdParty configuration for Nektar++
#
# SHARPy
#
########################################################################

OPTION(NEKTAR_USE_SHARPY
    "Use SHARPy routines for evaluating structural solutions" OFF)

IF (NEKTAR_USE_SHARPY)
    FIND_LIBRARY(SHARPY_LIBRARY NAMES "sharpy" PATHS /opt/local/lib)

    IF (SHARPY_LIBRARY)
        MESSAGE(STATUS "Found SHARPy: ${SHARPY_LIBRARY}")
        MARK_AS_ADVANCED(SHARPY_LIBRARY)
    ELSE()
        MESSAGE(FATAL_ERROR "Could not find SHARPY")
    ENDIF()
ENDIF()

