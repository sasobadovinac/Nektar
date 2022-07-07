SET(MATRIXFREEOPS_HEADERS
    Operator.hpp
)

SET(MATRIXFREEOPS_SOURCES
    Operator.cpp
)

# Boilerplate cpp implementation files for registering an operator
# with the factory. There are individual files for each operator shape
# so to reduce compile time timeout failures due to the large number
# of template impementations.

SET(OPERATORS BwdTrans PhysDeriv Helmholtz IProduct IProductWRTDerivBase )
SET(SHAPES Seg Tri Quad Hex Tet Prism Pyr )
SET(DIMENSIONS 1 2 2 3 3 3 3 )

# Note the commented out code is future simpler code for when cmake 3.17 can
# be used on the tester.

foreach(OPERATOR IN LISTS OPERATORS)
#  if(${CMAKE_VERSION} VERSION_LESS "3.17.0")
#    Remove the code in this part of the conditional once CMake 3.17
#    can be used on the CI tester.
     foreach(SHAPE IN LISTS SHAPES)
        IF( "${SHAPE}" STREQUAL "Seg" )
	    SET(DIMENSION 1)
        ELSEIF( "${SHAPE}" STREQUAL "Tri" OR "${SHAPE}" STREQUAL "Quad" )
	    SET(DIMENSION 2)
	ELSE()
	    SET(DIMENSION 3)
        ENDIF()
#  else()
#     Use this code once CMake 3.17 can be used on the CI tester.
#     foreach(SHAPE DIMENSION IN ZIP_LISTS SHAPES DIMENSIONS)
#  endif()

	# Uppercase is used for the #defines in the header files.
	string(TOUPPER ${SHAPE} TYPE)	

        # INST sets the macro used in NektarExpMacros.h
        # The BwdTrans only has regular elements.
        IF( "${OPERATOR}" STREQUAL "BwdTrans" )
            SET(INST 1)
        # All other operators have regular and deformed elements.
        ELSE()
            SET(INST 2)
        ENDIF()

        # The Helmholtz operator does not have a segment shape so skip.
        IF( "${OPERATOR}" STREQUAL "Helmholtz" AND "${SHAPE}" STREQUAL "Seg" )
             continue()
        ENDIF()

#        message(STATUS "operator=${OPERATOR} shape=${SHAPE} inst=${INST}")

        # Create the cpp file from the implementation file.
        configure_file(OperatorImp.cpp.in ${OPERATOR}${SHAPE}.cpp)

        # Add the respective headers and cpp file.
        SET(MATRIXFREEOPS_HEADERS ${MATRIXFREEOPS_HEADERS}
            ${OPERATOR}.h)
        SET(MATRIXFREEOPS_HEADERS ${MATRIXFREEOPS_HEADERS}
            ${OPERATOR}Kernels.hpp)
        SET(MATRIXFREEOPS_SOURCES ${MATRIXFREEOPS_SOURCES}
            ${CMAKE_CURRENT_BINARY_DIR}/${OPERATOR}${SHAPE}.cpp)

    endforeach()
endforeach()

# The current source dir is needed to get the header files.
ADD_DEFINITIONS(-DMATRIXFREE_EXPORTS -I${CMAKE_CURRENT_SOURCE_DIR})

ADD_NEKTAR_LIBRARY(MatrixFreeOps
    SOURCES ${MATRIXFREEOPS_SOURCES}
    HEADERS ${MATRIXFREEOPS_HEADERS}
    DEPENDS LibUtilities
    SUMMARY "Nektar++ matrix free operators library"
    DESCRIPTION "This library provides key matrix free operators that can use vector instructions.")

INSTALL(DIRECTORY ./
    DESTINATION ${NEKTAR_INCLUDE_DIR}/MatrixFreeOps
    COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
