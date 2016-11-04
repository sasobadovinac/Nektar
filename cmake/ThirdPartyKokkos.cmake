########################################################################
#
# ThirdParty configuration for Nektar++
#
# KOKKOS
#
########################################################################

OPTION(NEKTAR_USE_KOKKOS
    "Use KOKKOS routines for solver acceleration." OFF)

IF (NEKTAR_USE_KOKKOS)
    #SET(KOKKOS_SEARCH_PATHS $ENV{LD_LIBRARY_PATH} $ENV{KOKKOS_HOME}/lib)
    FIND_LIBRARY(KOKKOS_LIBRARY NAMES kokkos)

    IF (KOKKOS_LIBRARY)
        GET_FILENAME_COMPONENT(KOKKOS_PATH ${KOKKOS_LIBRARY} PATH)
        SET(KOKKOS_INCLUDE_DIR ${KOKKOS_PATH}/../include CACHE FILEPATH "KOKKOS include directory.")
        SET(BUILD_KOKKOS OFF)
    ELSE()
        SET(BUILD_KOKKOS ON)
    ENDIF ()

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_KOKKOS
        "Build KOKKOS from ThirdParty" ${BUILD_KOKKOS}
        "NEKTAR_USE_KOKKOS" OFF)

    IF (THIRDPARTY_BUILD_KOKKOS)
        CMAKE_DEPENDENT_OPTION(NEKTAR_USE_KOKKOS_CUDA
            "Build KOKKOS with CUDA support" OFF
            "NEKTAR_USE_KOKKOS" OFF)

        SET(KOKKOS_OPTIONS "--prefix=${TPDIST}" --cxxflags=-fPIC --with-serial --with-pthread --with-openmp)
        #SET(KOKKOS_OPTIONS "--prefix=${TPDIST}" --cxxflags=-fPIC --with-serial --with-pthread)
        
        IF (NEKTAR_USE_KOKKOS_CUDA)
            FIND_PACKAGE(CUDA REQUIRED VERSION 8.0)
            #SET(KOKKOS_OPTIONS ${KOKKOS_OPTIONS} --with-cuda=${CUDA_TOOLKIT_ROOT_DIR})
            SET(KOKKOS_OPTIONS
                ${KOKKOS_OPTIONS} --with-cuda=${CUDA_TOOLKIT_ROOT_DIR} --arch=Maxwell52 --with-cuda-options=enable_lambda)
        ENDIF()

        #SET(KOKKOS_BRANCH_NAME master)
        SET(KOKKOS_BRANCH_NAME develop)

        #message("KOKKOS BRANCH NAME = ${KOKKOS_BRANCH_NAME}")

        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            kokkos-${KOKKOS_BRANCH_NAME}
            #GIT_REPOSITORY https://github.com/kokkos/kokkos.git
            #GIT_TAG ${KOKKOS_BRANCH_NAME}
            URL ${CMAKE_SOURCE_DIR}/../kokkos-develop_old.zip

            STAMP_DIR ${TPBUILD}/stamp
            SOURCE_DIR ${TPSRC}/kokkos-${KOKKOS_BRANCH_NAME}
            BINARY_DIR ${TPBUILD}/kokkos-${KOKKOS_BRANCH_NAME}
            TMP_DIR ${TPBUILD}/kokkos-${KOKKOS_BRANCH_NAME}-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER} ${TPSRC}/kokkos-${KOKKOS_BRANCH_NAME}/generate_makefile.bash ${KOKKOS_OPTIONS}
            BUILD_COMMAND make lib
        )

        SET(KOKKOS_LIBRARY kokkos CACHE FILEPATH
            "KOKKOS library" FORCE)
        SET(KOKKOS_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "KOKKOS include" FORCE)

        LINK_DIRECTORIES(${TPDIST}/lib)

        MESSAGE(STATUS "Build KOKKOS: ${TPDIST}/lib/lib${KOKKOS_LIBRARY}.a")
        SET(KOKKOS_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE ()
        ADD_CUSTOM_TARGET(kokkos-${KOKKOS_BRANCH_NAME} ALL)
        MESSAGE(STATUS "Found KOKKOS: ${KOKKOS_LIBRARY}")
        SET(KOKKOS_CONFIG_INCLUDE_DIR ${KOKKOS_INCLUDE_DIR})
    ENDIF()

    ADD_DEFINITIONS(-DNEKTAR_USE_KOKKOS)

ENDIF( NEKTAR_USE_KOKKOS )

INCLUDE_DIRECTORIES(SYSTEM ${KOKKOS_INCLUDE_DIR})

MARK_AS_ADVANCED(KOKKOS_LIBRARY)
MARK_AS_ADVANCED(KOKKOS_INCLUDE_DIR)
