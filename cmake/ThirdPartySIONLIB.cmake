########################################################################
#
# ThirdParty configuration for Nektar++
#
# SIONlib
#
########################################################################

OPTION(NEKTAR_USE_SIONLIB
    "Enable SIONlib I/O support." ON)

IF (NEKTAR_USE_SIONLIB)
    INCLUDE (CheckIncludeFiles)

    IF (NOT NEKTAR_USE_MPI)
        MESSAGE(FATAL_ERROR "SIONlib requires Nektar++ to be configured with NEKTAR_USE_MPI for MPI support.")
    ENDIF()

    # Try to find SIONlib package.
    FIND_PACKAGE(SIONLIB)

    IF (NOT SIONLIB_FOUND)
        MESSAGE(FATAL_ERROR "SIONlib not detected.")
    ENDIF()
    
    SET(BUILD_SIONLIB OFF)

    SET(SIONLIB_CONFIG_INCLUDE_DIR ${SIONLIB_INCLUDE_DIR})
    
    MARK_AS_ADVANCED(SIONLIB_LIBRARIES)
    MARK_AS_ADVANCED(SIONLIB_INCLUDE_DIR)
    INCLUDE_DIRECTORIES(SYSTEM ${SIONLIB_INCLUDE_DIR})
ENDIF()
