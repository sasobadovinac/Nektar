########################################################################
#
# ThirdParty configuration for Nektar++
#
# KOKKOS
#
########################################################################

OPTION(NEKTAR_USE_KOKKOS
    "Use KOKKOS routines for performing the Fast Fourier Transform." OFF)

IF (NEKTAR_USE_KOKKOS)
    #SET(KOKKOS_SEARCH_PATHS $ENV{LD_LIBRARY_PATH} $ENV{KOKKOS_HOME}/lib)
    #FIND_LIBRARY(KOKKOS_LIBRARY NAMES fftw3 fftw3f PATHS ${KOKKOS_SEARCH_PATHS})

    #IF (KOKKOS_LIBRARY)
    #    GET_FILENAME_COMPONENT(KOKKOS_PATH ${KOKKOS_LIBRARY} PATH)
    #    SET(KOKKOS_INCLUDE_DIR ${KOKKOS_PATH}/../include CACHE FILEPATH "KOKKOS include directory.")
    #    SET(BUILD_KOKKOS OFF)
    #ELSE()
    #    SET(BUILD_KOKKOS ON)
    #ENDIF ()

    #CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_KOKKOS
    #    "Build KOKKOS from ThirdParty" ${BUILD_KOKKOS}
    #    "NEKTAR_USE_KOKKOS" OFF)

    #IF (THIRDPARTY_BUILD_KOKKOS)
        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            kokkos-master
            GIT_REPOSITORY https://github.com/kokkos/kokkos.git
            GIT_TAG master
            STAMP_DIR ${TPBUILD}/stamp
            SOURCE_DIR ${TPSRC}/kokkos-master
            BINARY_DIR ${TPBUILD}/kokkos-master
            TMP_DIR ${TPBUILD}/kokkos-master-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER} ${TPSRC}/kokkos-master/generate_makefile.bash --prefix=${TPDIST}
            BUILD_COMMAND make lib
        )

        SET(KOKKOS_LIBRARY kokkos CACHE FILEPATH
            "KOKKOS library" FORCE)
        SET(KOKKOS_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "KOKKOS include" FORCE)

        LINK_DIRECTORIES(${TPDIST}/lib)

        MESSAGE(STATUS "Build KOKKOS: ${TPDIST}/lib/lib${KOKKOS_LIBRARY}.a")
        SET(KOKKOS_CONFIG_INCLUDE_DIR ${TPINC})
    #ELSE ()
    #    ADD_CUSTOM_TARGET(kokkos-master ALL)
    #    MESSAGE(STATUS "Found KOKKOS: ${KOKKOS_LIBRARY}")
    #    SET(KOKKOS_CONFIG_INCLUDE_DIR ${KOKKOS_INCLUDE_DIR})
    #ENDIF()
ENDIF( NEKTAR_USE_KOKKOS )

INCLUDE_DIRECTORIES(SYSTEM ${KOKKOS_INCLUDE_DIR})

MARK_AS_ADVANCED(KOKKOS_LIBRARY)
MARK_AS_ADVANCED(KOKKOS_INCLUDE_DIR)
