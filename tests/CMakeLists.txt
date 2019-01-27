CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/tests/Tester.cpp.in
               ${CMAKE_BINARY_DIR}/tests/Tester.cpp
	       @ONLY)
SET(TESTER_SOURCES
  Metric.cpp
  MetricEigenvalue.cpp
  MetricFile.cpp
  MetricL2.cpp
  MetricLInf.cpp
  MetricPrecon.cpp
  MetricRegex.cpp
  TestData.cpp
  ${CMAKE_BINARY_DIR}/tests/Tester.cpp
  sha1.cpp
)

SET(TESTER_HEADERS
  Metric.h
  MetricEigenvalue.h
  MetricFile.h
  MetricL2.h
  MetricLInf.h
  MetricPrecon.h
  MetricRegex.h
  TestData.h
  Tester.h
  sha1.h
)

ADD_DEFINITIONS(-DBUILD_PATH="${CMAKE_BINARY_DIR}")

IF(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  ADD_DEFINITIONS(-DRELWITHDEBINFO)
ENDIF()


ADD_NEKTAR_EXECUTABLE(Tester tests TESTER_SOURCES)

TARGET_LINK_LIBRARIES(Tester
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_REGEX_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_THREAD_LIBRARY}
  optimized ${TINYXML_LIBRARY} debug ${TINYXML_LIBRARY}
)
ADD_DEPENDENCIES(Tester boost tinyxml-2.6.2)

IF( ${CMAKE_SYSTEM} MATCHES "Linux.*")
  # The boost thread library needs pthread on linux.
  GET_TARGET_PROPERTY(THE_COMPILE_FLAGS Tester COMPILE_FLAGS)
  GET_TARGET_PROPERTY(THE_LINK_FLAGS Tester LINK_FLAGS)

  # It is possible these flags have not been set yet.
  IF(NOT THE_COMPILE_FLAGS)
      SET(THE_COMPILE_FLAGS "")
  ENDIF(NOT THE_COMPILE_FLAGS)

  IF(NOT THE_LINK_FLAGS )
      SET(THE_LINK_FLAGS "")
  ENDIF(NOT THE_LINK_FLAGS)

  SET_TARGET_PROPERTIES(Tester
          PROPERTIES COMPILE_FLAGS "${THE_COMPILE_FLAGS} -pthread")
  SET_TARGET_PROPERTIES(Tester
          PROPERTIES LINK_FLAGS "${THE_LINK_FLAGS} -pthread")

  TARGET_LINK_LIBRARIES(Tester optimized rt debug rt)
ENDIF()

IF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")
  SET_TARGET_PROPERTIES(Tester
    PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup -Wl,-rpath,${CMAKE_INSTALL_PREFIX}/${LIB_DIR} -Wl,-rpath,${Boost_LIBRARY_DIRS}")
ENDIF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")

SET_TARGET_PROPERTIES(Tester
    PROPERTIES COMPILE_FLAGS "-UNEKTAR_USE_MPI")
