########################################################################
#
# ThirdParty configuration for Nektar++
#
# Eigen 
#
########################################################################

# First search for system Eigen installs. Hint eigen3 for MacPorts.
FIND_PACKAGE(Eigen3)

# If we have our library then don't build Eigen.
IF (NOT EIGEN3_INCLUDE_DIR)
    SET(BUILD_EIGEN3 ON)
ELSE()
    SET(BUILD_EIGEN3 OFF)
ENDIF ()

OPTION(THIRDPARTY_BUILD_EIGEN3
    "Build Eigen library from ThirdParty." ${BUILD_EIGEN3})

IF (THIRDPARTY_BUILD_EIGEN3)
    INCLUDE(ExternalProject)

    EXTERNALPROJECT_ADD(
        eigen3-3.3.7
        PREFIX ${TPSRC}
        URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.bz2
        URL_MD5 b9e98a200d2455f06db9c661c5610496
        STAMP_DIR ${TPBUILD}/stamp
        DOWNLOAD_DIR ${TPSRC}
        SOURCE_DIR ${TPSRC}/eigen-3.3.7
        BINARY_DIR ${TPBUILD}/eigen-3.3.7
        TMP_DIR ${TPBUILD}/eigen-3.3.7-tmp
        INSTALL_DIR ${TPDIST}
        CONFIGURE_COMMAND ${CMAKE_COMMAND}
            -G ${CMAKE_GENERATOR}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
            ${TPSRC}/eigen-3.3.7
        )

    SET(EIGEN3_INCLUDE_DIR ${TPDIST}/include/eigen3 CACHE FILEPATH
        "eigen include" FORCE)
    SET(EIGEN3_CONFIG_INCLUDE_DIR ${TPINC})
    MESSAGE(STATUS "Build Eigen: ${EIGEN3_INCLUDE_DIR}")
ELSE()
    ADD_CUSTOM_TARGET(eigen-3.3.7 ALL)
    MESSAGE(STATUS "Found Eigen: ${EIGEN3_INCLUDE_DIR}")
    SET(EIGEN3_CONFIG_INCLUDE_DIR ${EIGEN3_INCLUDE_DIR})
ENDIF()

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

MARK_AS_ADVANCED(EIGEN3_INCLUDE_DIR)
MARK_AS_ADVANCED(EIGEN3_CONFIG_INCLUDE_DIR)

