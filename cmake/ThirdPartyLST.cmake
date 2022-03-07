########################################################################
#
# ThirdParty confiuration for Nektar++
#
# Linear stability analysis (LST) module in FieldConvert
#
########################################################################

# TODO:
# 1. Set default off
# 2. hello world in C++ [done]
# 3. hello world in Fortran (f90) [done]
# 4. Shahid's code [done]

OPTION(NEKTAR_USE_LST "Use LST module in FieldConvert" ON)

IF(NEKTAR_USE_LST)

    INCLUDE(ExternalProject)
    
    EXTERNALPROJECT_ADD(
        lst-0.1
        URL                ${TPURL}/lst_v1.6.zip
        URL_MD5            5820c5c37016f036b0ecffe289b2f761
        PREFIX             ${TPSRC}
        STAMP_DIR          ${TPBUILD}/stamp
        DOWNLOAD_DIR       ${TPSRC}
        SOURCE_DIR         ${TPSRC}/lst-0.1
        BINARY_DIR         ${TPBUILD}/lst-0.1
        TMP_DIR            ${TPBUILD}/lst-0.1-tmp
        INSTALL_DIR        ${TPDIST}
        
        CONFIGURE_COMMAND ${CMAKE_COMMAND} 
        -G ${CMAKE_GENERATOR}
        -DCMAKE_INSTALL_LIBDIR:PATH=${TPDIST}/lib
        -DCMAKE_INSTALL_INCDIR:PATH=${TPDIST}/include
        ${TPSRC}/lst-0.1 )

    THIRDPARTY_LIBRARY(LST_LIBRARY SHARED lst DESCRIPTION "Linear stability analysis library")
    MARK_AS_ADVANCED(LST_LIBRARY)
    MESSAGE(STATUS "Build LST: ${LST_LIBRARY}")

ENDIF(NEKTAR_USE_LST)




