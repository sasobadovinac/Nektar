# FFTW
OPTION(NEKTAR_USE_FFTW
    "Use FFTW routines for performing the Fast Fourier Transform." OFF)

CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_FFTW
    "Build FFTW from ThirdParty" OFF
    "NEKTAR_USE_FFTW" OFF)

IF( NEKTAR_USE_FFTW )
    IF (THIRDPARTY_BUILD_FFTW)
        INCLUDE(ExternalProject)
        
        EXTERNALPROJECT_ADD(
            fftw-3.2.2
            URL ${TPURL}/fftw-3.2.2.tar.gz
            URL_MD5 "b616e5c91218cc778b5aa735fefb61ae"
            STAMP_DIR ${TPSRC}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/fftw-3.2.2
            BINARY_DIR ${TPBUILD}/fftw-3.2.2
            TMP_DIR ${TPBUILD}/fftw-3.2.2-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ${TPSRC}/fftw-3.2.2/configure --prefix=${TPDIST} --quiet --enable-shared --disable-dependency-tracking
        )
        SET(FFTW_LIB fftw3)
        MARK_AS_ADVANCED(FFTW_LIB)
        INCLUDE_DIRECTORIES(${TPDIST}/include)
        MESSAGE(STATUS "Build FFTW: ${TPDIST}/lib/lib${FFTW_LIB}.so")
    ELSE ()
        INCLUDE (FindFFTW)
    ENDIF()
ENDIF( NEKTAR_USE_FFTW )
