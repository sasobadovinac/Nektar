MACRO(ADD_SOLVER_EXECUTABLE name component source)

    ADD_EXECUTABLE(${name} ${source} ${ARGN})
    
    SET_TARGET_PROPERTIES(${name} PROPERTIES DEBUG_POSTFIX -g)
    SET_TARGET_PROPERTIES(${name} PROPERTIES MINSIZEREL_POSTFIX -ms)
    SET_TARGET_PROPERTIES(${name} PROPERTIES RELWITHDEBINFO_POSTFIX -rg)

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

    TARGET_LINK_LIBRARIES(${name} SolverUtils)

    SET_PROPERTY(TARGET ${name} PROPERTY FOLDER ${component})
    INSTALL(TARGETS ${name} 
            RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${component} OPTIONAL)
ENDMACRO(ADD_SOLVER_EXECUTABLE name component source)


# Generate list of available subdirectories
FILE(GLOB dir_list "*")
FOREACH(dir ${dir_list})
    IF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)    
        ADD_SUBDIRECTORY(${dir})
    ENDIF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
ENDFOREACH(dir ${dir_list})


SET(CPACK_COMPONENTS_ALL solvers)

SET(CPACK_COMPONENT_SOLVERS_NAME "nektar++-solvers")
SET(CPACK_COMPONENT_SOLVERS_DESCRIPTION "Solvers for Nektar++")
#SET(CPACK_COMPONENT_SOLVERS_DEPENDS lib)

INCLUDE(CPack)
