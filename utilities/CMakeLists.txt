MACRO(ADD_UTILITIES_EXECUTABLE name component source)
    ADD_EXECUTABLE(${name} ${source} ${ARGN})
    SET_COMMON_PROPERTIES(${name})

    IF( ${CMAKE_SYSTEM} MATCHES "Linux.*" )
        # The boost thread library needs pthread on linux.
        GET_TARGET_PROPERTY(THE_COMPILE_FLAGS ${name} COMPILE_FLAGS)
        GET_TARGET_PROPERTY(THE_LINK_FLAGS ${name} LINK_FLAGS)

        # It is possible these flags have not been set yet.
        IF(NOT THE_COMPILE_FLAGS)
            SET(THE_COMPILE_FLAGS "")
        ENDIF(NOT THE_COMPILE_FLAGS)

        IF(NOT THE_LINK_FLAGS )
            SET(THE_LINK_FLAGS "")
        ENDIF(NOT THE_LINK_FLAGS)

        SET_TARGET_PROPERTIES(${name}
                PROPERTIES COMPILE_FLAGS "${THE_COMPILE_FLAGS} -pthread")
        SET_TARGET_PROPERTIES(${name}
                PROPERTIES LINK_FLAGS "${THE_LINK_FLAGS} -pthread")

    ENDIF( ${CMAKE_SYSTEM} MATCHES "Linux.*" )

    TARGET_LINK_LIBRARIES(${name} MultiRegions)
    SET_PROPERTY(TARGET ${name} PROPERTY FOLDER ${component})
    INSTALL(TARGETS ${name}
	RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${component})
ENDMACRO()

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library ${CMAKE_SOURCE_DIR}/utilities)

ADD_UTILITIES_EXECUTABLE(FldAddFld         util FldAddFld.cpp)
ADD_UTILITIES_EXECUTABLE(FldAddScalGrad    util FldAddScalGrad.cpp)
ADD_UTILITIES_EXECUTABLE(FldAddScalGrad_elmt util FldAddScalGrad_elmt.cpp)
ADD_UTILITIES_EXECUTABLE(FldAddWSS         util FldAddWSS.cpp)
ADD_UTILITIES_EXECUTABLE(FldAddMultiShear  util FldAddMultiShear.cpp)
ADD_UTILITIES_EXECUTABLE(MeshConvert       util MeshConvert.cpp)
ADD_UTILITIES_EXECUTABLE(ProbeFld          util-extra ProbeFld.cpp)
ADD_UTILITIES_EXECUTABLE(XmlToTecplot      util XmlToTecplot.cpp)
ADD_UTILITIES_EXECUTABLE(XmlToTecplotWireFrame  util XmlToTecplotWireFrame.cpp)
ADD_UTILITIES_EXECUTABLE(XmlToVtk          util XmlToVtk.cpp)

#SUBDIRS(PostProcessing FieldConvert NekMesh SimpleMeshGen Extras)

# Generate list of available utilities
FILE(GLOB dir_list "*")
LIST(REMOVE_ITEM dir_list "${CMAKE_SOURCE_DIR}/.svn")
FOREACH(dir ${dir_list})
    IF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
        ADD_SUBDIRECTORY(${dir})
    ENDIF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
ENDFOREACH(dir ${dir_list})
