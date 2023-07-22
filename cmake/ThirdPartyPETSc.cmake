########################################################################
#
# ThirdParty configuration for Nektar++
#
# PETSc
#
########################################################################

OPTION(NEKTAR_USE_PETSC
    "Enable PETSc parallel matrix solver support." OFF)

IF (NEKTAR_USE_PETSC)
    SET(PETSC_FIND_QUIETLY ON)
    FIND_PACKAGE(PkgConfig REQUIRED)
    PKG_SEARCH_MODULE(PETSC PETSc IMPORTED_TARGET)

    IF (PETSC_FOUND)
        # Prefer to use absolute paths to avoid separate link directories on
        # some platforms (e.g. macOS). However not available in all CMake
        # versions.
        IF (DEFINED PETSC_LINK_LIBRARIES)
            SET(PETSC_LIBRARIES ${PETSC_LINK_LIBRARIES} CACHE INTERNAL "")
        ENDIF()
        MESSAGE(STATUS "Found PETSc ${PETSC_VERSION}: ${PETSC_LIBRARIES}")
        SET(BUILD_PETSC OFF)
    ELSE()
        MESSAGE(STATUS "Could not find system-wide PETSc installation.")
        SET(BUILD_PETSC ON)
    ENDIF()

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_PETSC
        "Build PETSc if needed" ${BUILD_PETSC}
        "NEKTAR_USE_PETSC" OFF)

    IF (THIRDPARTY_BUILD_PETSC)
        INCLUDE(ExternalProject)

        FIND_PACKAGE(PythonInterp 2 REQUIRED)

        UNSET(PATCH CACHE)
        FIND_PROGRAM(PATCH patch)
        IF(NOT PATCH)
            MESSAGE(FATAL_ERROR
                "'patch' tool for modifying files not found. Cannot build PETSc.")
        ENDIF()
        MARK_AS_ADVANCED(PATCH)

        SET(PETSC_C_COMPILER   "${CMAKE_C_COMPILER}")
        SET(PETSC_CXX_COMPILER "${CMAKE_CXX_COMPILER}")
        SET(PETSC_Fortran_COMPILER "${CMAKE_Fortran_COMPILER}")

        IF (NEKTAR_USE_MPI)
            IF (NOT MPI_BUILTIN)
                SET(PETSC_C_COMPILER   "${MPI_C_COMPILER}")
                SET(PETSC_CXX_COMPILER "${MPI_CXX_COMPILER}")
                SET(PETSC_Fortran_COMPILER "${MPI_Fortran_COMPILER}")
            ENDIF (NOT MPI_BUILTIN)
        ELSE (NEKTAR_USE_MPI)
            SET(PETSC_NO_MPI "--with-mpi=0")
        ENDIF (NEKTAR_USE_MPI)

        IF(CMAKE_Fortran_COMPILER AND NEKTAR_USE_MPI)
            IF(NOT MPI_Fortran_COMPILER)
                MESSAGE(ERROR "MPI_Fortran_COMPILER not set")
            ENDIF()
            # we use a MUMPS build in ordering here, in the future it might make
            # sense to hook it up with metis/scotch since this MIGHT be faster
            SET(PETSC_MUMPS --download-scalapack --download-mumps)
            SET(PETSC_DEPS "")

            IF(NOT BLAS_LAPACK_BUILTIN )
                LIST(GET BLAS_LAPACK 0 LAPACK)
                LIST(GET BLAS_LAPACK 1 BLAS)
                GET_FILENAME_COMPONENT(BLAS_LAPACK_DIR ${BLAS} PATH)
                IF( NEKTAR_USE_MKL OR NEKTAR_USE_SYSTEM_BLAS_LAPACK)
                    SET(PETSC_MUMPS ${PETSC_MUMPS} --with-blas-lib=${BLAS}
                        --with-lapack-lib=${LAPACK})
                    IF(THIRDPARTY_BUILD_BLAS_LAPACK)
                        SET(PETSC_DEPS ${PETSC_DEPS} lapack-3.7.1)
                    ENDIF()
                ELSE()
                    MESSAGE(STATUS "No suitable blas/lapack found, downloading")
                    SET(PETSC_MUMPS ${PETSC_MUMPS} --download-fblaslapack)
                ENDIF()
            ELSE()
                MESSAGE(STATUS "No suitable blas/lapack found, downloading")
                SET(PETSC_MUMPS ${PETSC_MUMPS} --download-fblaslapack)
            ENDIF()

        ELSE()
            MESSAGE(WARNING "No MPI and/or Fortran support. Building PETSc without MUMPS support")
            SET(PETSC_Fortran_COMPILER "0")
        ENDIF()

        EXTERNALPROJECT_ADD(
            petsc-3.11.4
            DEPENDS ${PETSC_DEPS}
            PREFIX ${TPSRC}
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPBUILD}/petsc-3.11.4
            TMP_DIR ${TPBUILD}/petsc-3.11.4-tmp
            INSTALL_DIR ${TPDIST}
            BINARY_DIR ${TPBUILD}/petsc-3.11.4
            URL https://www.nektar.info/thirdparty/petsc-lite-3.11.4.tar.gz
            URL_MD5 "33da4b6a430d5d9e13b19871d707af0f"
            CONFIGURE_COMMAND
                OMPI_FC=${CMAKE_Fortran_COMPILER}
                OMPI_CC=${CMAKE_C_COMPILER}
                OMPI_CXX=${CMAKE_CXX_COMPILER}
                ${PYTHON_EXECUTABLE} ./configure
                MAKEFLAGS=$MAKEFLAGS
                CFLAGS="-w"
                CXXFLAGS="-w"
                --with-fc=${PETSC_Fortran_COMPILER}
                --with-cc=${PETSC_C_COMPILER}
                --with-cxx=${PETSC_CXX_COMPILER}
                --with-shared-libraries=1
                --with-make-np=1
                --with-pic=1
                --with-x=0
                --with-ssl=0
                --prefix=${TPDIST}
                --with-petsc-arch=c-opt
                ${PETSC_MUMPS}
                ${PETSC_NO_MPI}
            BUILD_COMMAND $(MAKE)
            TEST_COMMAND $(MAKE)
                PETSC_DIR=${TPDIST} PETSC_ARCH=c-opt test)

        THIRDPARTY_LIBRARY(PETSC_LIBRARIES SHARED petsc
            DESCRIPTION "PETSc library")
        SET(PETSC_INCLUDE_DIRS ${TPDIST}/include CACHE FILEPATH
            "PETSc includes" FORCE)
        MESSAGE(STATUS "Build PETSc: ${PETSC_LIBRARIES}")
        SET(PETSC_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE (THIRDPARTY_BUILD_PETSC)
        SET(PETSC_CONFIG_INCLUDE_DIR ${PETSC_INCLUDE_DIRS})
        ADD_CUSTOM_TARGET(petsc-3.11.4 ALL)
    ENDIF (THIRDPARTY_BUILD_PETSC)

    ADD_DEFINITIONS(-DNEKTAR_USING_PETSC)
    INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIRS})
    IF (NOT NEKTAR_USE_MPI)
        INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIRS}/petsc/mpiuni)
    ENDIF (NOT NEKTAR_USE_MPI)

    MARK_AS_ADVANCED(PETSC_CURRENT PETSC_DIR PETSC_LIBRARIES PETSC_INCLUDES)
ENDIF( NEKTAR_USE_PETSC )
