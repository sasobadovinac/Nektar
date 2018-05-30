# Module that checks whether SIONLIB is available.
#
# Variables used by this module which you may want to set:
# SIONLIB_ROOTDIR   Path list to search for SIONLIB
# SIONLIB_SUFFIX    suffix to the library name , e.g. gcc or something
# SIONLIB_INCDIR    directory with SIONLIB headers inside
# SIONLIB_LIBDIR    directory with SIONLIB libraries inside
#
# Sets the following variables
#
# SIONLIB_INCLUDE_DIR    Path to the SIONLIB include dir
# SIONLIB_MPI_LIBRARY
# SIONLIB_GEN_LIBRARY
# SIONLIB_SER_LIBRARY
# SIONLIB_COM_LIBRARY
# SIONLIB_COMLOCKNONE_LIBRARY
# SIONLIB_LIBRARY_DIR
# SIONLIB_FOUND          True if SIONLIB was found and usable
# HAVE_SIONLIB           True if SIONLIB was found and usable
# SIONLIB_LIBRARIES      Names of the SIONLIB libraries
#

set(SIONLIB_ROOTDIR $ENV{SIONLIB_DIR} CACHE PATH "Path list to search for SIONLIB" FORCE)
set(SIONLIB_SUFFIX "_lib64" CACHE STRING "suffix to the library name , e.g. gcc or something")
set(SIONLIB_INCDIR $ENV{SIONLIB_INC} CACHE PATH "directory with SIONLIB headers inside" FORCE)
set(SIONLIB_LIBDIR $ENV{SIONLIB_LIB} CACHE PATH "directory with SIONLIB libraries inside" FORCE)

MESSAGE(STATUS "SIONLIB_ROOTDIR=${SIONLIB_ROOTDIR}")
MESSAGE(STATUS "SIONLIB_INCDIR=${SIONLIB_INCDIR}")
MESSAGE(STATUS "SIONLIB_LIBDIR=${SIONLIB_LIBDIR}")

mark_as_advanced(SIONLIB_ROOTDIR SIONLIB_SUFFIX SIONLIB_INCDIR SIONLIB_LIBDIR)


#look for header files at positions given by the user
find_path(SIONLIB_INCLUDE_DIR
  NAMES "sion.h"
  PATHS ${SIONLIB_INCDIR}
  PATH_SUFFIXES "include" 
  NO_DEFAULT_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
)

MESSAGE(STATUS "SIONLIB_INCLUDE_DIR=${SIONLIB_INCLUDE_DIR}")

# check header usability
include(CMakePushCheckState)
cmake_push_check_state()
set(CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${MPI_INCLUDE_PATH} ${SIONLIB_INCLUDE_DIR})
set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${MPI_LIBRARIES})
check_include_files(sion.h SIONLIB_HEADER_USABLE)


find_library(SIONLIB_MPI_LIBRARY
  NAMES "sionmpi_64"
  PATHS ${SIONLIB_LIBDIR}
  PATH_SUFFIXES "lib"
  NO_DEFAULT_PATH
)

MESSAGE(STATUS "SIONLIB_MPI_LIBRARY=${SIONLIB_MPI_LIBRARY}")

find_library(SIONLIB_GEN_LIBRARY
  NAMES "siongen_64"
  PATHS ${SIONLIB_LIBDIR}
  PATH_SUFFIXES "lib"
  NO_DEFAULT_PATH
)

MESSAGE(STATUS "SIONLIB_GEN_LIBRARY=${SIONLIB_GEN_LIBRARY}")

find_library(SIONLIB_SER_LIBRARY
  NAMES "sionser_64"
  PATHS ${SIONLIB_LIBDIR}
  PATH_SUFFIXES "lib"
  NO_DEFAULT_PATH
)

MESSAGE(STATUS "SIONLIB_SER_LIBRARY=${SIONLIB_SER_LIBRARY}")

find_library(SIONLIB_COM_LIBRARY
  NAMES "sioncom_64"
  PATHS ${SIONLIB_LIBDIR}
  PATH_SUFFIXES "lib"
  NO_DEFAULT_PATH
)

MESSAGE(STATUS "SIONLIB_COM_LIBRARY=${SIONLIB_COM_LIBRARY}")

find_library(SIONLIB_COMLOCKNONE_LIBRARY
  NAMES "sioncom_64_lock_none"
  PATHS ${SIONLIB_LIBDIR}
  PATH_SUFFIXES "lib"
  NO_DEFAULT_PATH
)

MESSAGE(STATUS "SIONLIB_COMLOCKNONE_LIBRARY=${SIONLIB_COMLOCKNONE_LIBRARY}")


cmake_pop_check_state()

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "SIONlib"
  DEFAULT_MSG
  SIONLIB_INCLUDE_DIR
  SIONLIB_MPI_LIBRARY
  SIONLIB_GEN_LIBRARY
  SIONLIB_SER_LIBRARY
  SIONLIB_COM_LIBRARY
  SIONLIB_COMLOCKNONE_LIBRARY
)

mark_as_advanced(SIONLIB_INCLUDE_DIR SIONLIB_MPI_LIBRARY SIONLIB_GEN_LIBRARY SIONLIB_SER_LIBRARY SIONLIB_COM_LIBRARY SIONLIB_COMLOCKNONE_LIBRARY)

# if both headers and library are found, store results
if(SIONLIB_FOUND)
  set(SIONLIB_LIBRARY_DIR ${SIONLIB_LIBDIR})
  set(SIONLIB_LIBRARIES ${SIONLIB_MPI_LIBRARY} ${SIONLIB_GEN_LIBRARY} ${SIONLIB_SER_LIBRARY} ${SIONLIB_COM_LIBRARY} ${SIONLIB_COMLOCKNONE_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of SIONLIB succeded:\n"
    "Include directory: ${SIONLIB_INCLUDE_DIR}\n"
    "Library directory: ${SIONLIB_LIBRARY_DIR}\n\n")
  set(SIONLIB_COMPILE_FLAGS "-I${SIONLIB_INCLUDE_DIR} -D_SION_XT -DSION_MPI"
    CACHE STRING "Compile Flags used by Nektar++ when compiling with SIONLIB libraries")
  set(SIONLIB_LIBRARIES ${SIONLIB_LIBRARIES} 
    CACHE STRING "Libraries used by Nektar++ when linking with SIONLIB libraries")
else(SIONLIB_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of SIONLIB failed:\n"
    "Include directory: ${SIONLIB_INCLUDE_DIR}\n"
    "Library directory: ${SIONLIB_LIBRARY_DIR}\n\n")
endif(SIONLIB_FOUND)

MESSAGE(STATUS "SIONLIB_LIBRARIES=${SIONLIB_LIBRARIES}")

#set HAVE_SIONLIB for config.h
set(HAVE_SIONLIB ${SIONLIB_FOUND})

#add all sionlib related flags to ALL_PKG_FLAGS, this must happen regardless of a target using add_sionlib_flags
if(SIONLIB_FOUND)
  set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "${SIONLIB_COMPILE_FLAGS}")
  foreach(dir "${SIONLIB_INCLUDE_DIR}")
    set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "-I${dir}")
  endforeach()
endif()