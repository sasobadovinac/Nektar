########################################################################
#
# ThirdParty configuration for Nektar++
#
# TETGEN
#
########################################################################

IF(NEKTAR_USE_MESHGEN)
    # Search for system-installed Triangle installation
    FIND_LIBRARY(TETGEN_LIBRARY NAMES tet)
    FIND_PATH(TETGEN_INCLUDE_DIR tetgen.h)

    IF(TETGEN_LIBRARY AND TETGEN_INCLUDE_DIR)
        SET(BUILD_TETGEN OFF)
    ELSE()
        SET(BUILD_TETGEN ON)
    ENDIF()

    OPTION(THIRDPARTY_BUILD_TETGEN
        "Build TetGen library from ThirdParty." ${BUILD_TETGEN})

    IF (THIRDPARTY_BUILD_TETGEN)
        INCLUDE(ExternalProject)

        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
            MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build Tetgen.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        EXTERNALPROJECT_ADD(
            tetgen-1.5
            PREFIX ${TPSRC}
            URL ${TPURL}/tetgen-1.5.zip
            URL_MD5 6d62e63f9b1e7a8ce53d5bc87e6a0a09
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/tetgen-1.5
            BINARY_DIR ${TPBUILD}/tetgen-1.5
            TMP_DIR ${TPBUILD}/tetgen-1.5-tmp
            INSTALL_DIR ${TPDIST}
            PATCH_COMMAND ${PATCH} -p1 < ${PROJECT_SOURCE_DIR}/cmake/thirdparty-patches/tetgen-snprintf.patch
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
            -G ${CMAKE_GENERATOR}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
            ${TPSRC}/tetgen-1.5
            )
        THIRDPARTY_LIBRARY(TETGEN_LIBRARY STATIC tetgen
            DESCRIPTION "Tetgen library")
        SET(TETGEN_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "TetGen include" FORCE)
        ADD_DEFINITIONS(-DTETGEN_HAS_DEINITIALIZE)
        MESSAGE(STATUS "Build TetGen: ${TETGEN_LIBRARY}")
        SET(TETGEN_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE()
        ADD_CUSTOM_TARGET(tetgen-1.5 ALL)
        MESSAGE(STATUS "Found TetGen: ${TETGEN_LIBRARY}")
        SET(TETGEN_CONFIG_INCLUDE_DIR ${TETGEN_INCLUDE_DIR})

        # Test whether tetgen includes the deinitialize() routine (not include in 1.6.0
        # but is present in 1.5.0.
        INCLUDE(CheckCXXSourceCompiles)
        SET(CMAKE_REQUIRED_INCLUDES "${TETGEN_INCLUDE_DIR}")
        SET(CMAKE_REQUIRED_LIBRARIES "${TETGEN_LIBRARY}")
        CHECK_CXX_SOURCE_COMPILES("
            #include <tetgen.h>
            int main() {
                tetgenio tmp;
                tmp.deinitialize();
                return 0;
            }
            " TETGEN_HAS_DEINITIALIZE)
        IF (TETGEN_HAS_DEINITIALIZE)
            ADD_DEFINITIONS(-DTETGEN_HAS_DEINITIALIZE)
        ENDIF()
    ENDIF()

    MARK_AS_ADVANCED(TETGEN_LIBRARY)
    MARK_AS_ADVANCED(TETGEN_INCLUDE_DIR)
    INCLUDE_DIRECTORIES(${TETGEN_INCLUDE_DIR})
ENDIF(NEKTAR_USE_MESHGEN)
