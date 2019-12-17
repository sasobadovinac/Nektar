########################################################################
#
# ThirdParty configuration for Nektar++
#
# Eigen 
#
########################################################################

# First search for system Eigen installs. Hint /opt/local for MacPorts.
FIND_PATH   (EIGEN_INCLUDE_DIR Eigen/Dense)
FIND_LIBRARY(EIGEN_LIBRARY NAMES "Eigen")

# If we have our library then don't build Eigen.
IF (EIGEN_INCLUDE_DIR AND EIGEN_LIBRARY)
    SET(BUILD_EIGEN OFF)
ELSE()
    SET(BUILD_EIGEN ON)
ENDIF ()

OPTION(THIRDPARTY_BUILD_EIGEN
    "Build eigen library from ThirdParty." ${BUILD_EIGEN})

IF (THIRDPARTY_BUILD_EIGEN)
    INCLUDE(ExternalProject)

        EXTERNALPROJECT_ADD(
            eigen-3.3.7 
            PREFIX ${TPSRC}
            URL ${TPURL2}/eigen-3.3.7.tar.bz2
           # URL_MD5 240beaeb45f63b154c9801eef7561eac
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

  # THIRDPARTY_LIBRARY(EIGEN_LIBRARY STATIC eigen DESCRIPTION "Eigen library")
    SET(EIGEN_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
        "eigen include" FORCE)
  # MESSAGE(STATUS "Build Eigen: ${EIGEN_LIBRARY}")
    SET(EIGEN_CONFIG_INCLUDE_DIR ${TPINC})
ELSE()
    ADD_CUSTOM_TARGET(eigen-3.3.7 ALL)
 #   MESSAGE(STATUS "Found Eigen: ${EIGEN_LIBRARY}")
    SET(EIGEN_CONFIG_INCLUDE_DIR ${EIGEN_INCLUDE_DIR})
ENDIF (THIRDPARTY_BUILD_EIGEN)

INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIR})

MARK_AS_ADVANCED(EIGEN_INCLUDE_DIR)
#MARK_AS_ADVANCED(EIGEN_LIBRARY)
MARK_AS_ADVANCED(EIGEN_CONFIG_INCLUDE_DIR)
