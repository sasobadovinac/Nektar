########################################################################
#
# ThirdParty configuration for Nektar++
#
# PETSc
#
########################################################################

CMAKE_DEPENDENT_OPTION(NEKTAR_USE_SAENA
    "Enable Saena parallel matrix solver support." OFF
    "NEKTAR_USE_MKL;NEKTAR_USE_MPI" ON)

IF (NEKTAR_USE_SAENA)
    SET(BUILD_SAENA ON)

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_SAENA
        "Build Saena if needed" ${BUILD_SAENA}
        "NEKTAR_USE_SAENA" OFF)

    IF (THIRDPARTY_BUILD_SAENA)
        INCLUDE(ExternalProject)

        EXTERNALPROJECT_ADD(
            saena
            PREFIX ${TPSRC}
            STAMP_DIR ${TPBUILD}/stamp
            GIT_REPOSITORY https://github.com/mdave/Saena_Public.git
            GIT_TAG 298fc0249aa05fd03c27108615c8f2afe8c111bf
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPBUILD}/saena
            TMP_DIR ${TPBUILD}/saena-tmp
            INSTALL_DIR ${TPDIST}
            BINARY_DIR ${TPBUILD}/saena
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ./install.sh
            )

        THIRDPARTY_LIBRARY(SAENA_LIBRARY SHARED saena
            DESCRIPTION "Saena library")
        SET(SAENA_INCLUDE_DIR ${TPDIST}/include CACHE FILEPATH
            "Saenae includes" FORCE)
        MESSAGE(STATUS "Build Saena: ${SAENA_LIBRARY}")
        SET(SAENA_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE()
        MESSAGE(FATAL "Saena only available through third-party install")
    ENDIF()

    ADD_DEFINITIONS(-DNEKTAR_USING_SAENA)
    INCLUDE_DIRECTORIES(${SAENA_INCLUDE_DIR})
    MARK_AS_ADVANCED(SAENA_LIBRARY SAENA_INCLUDE_DIR)
ENDIF()
