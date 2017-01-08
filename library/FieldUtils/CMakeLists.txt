SET(FieldUtilsHeaders
    Module.h
    Field.hpp
    Interpolator.h
    InputModules/InputDat.h
    InputModules/InputFld.h
    InputModules/InputXml.h
    InputModules/InputPts.h
    InputModules/InputNek5000.h
    OutputModules/OutputInfo.h
    OutputModules/OutputTecplot.h
    OutputModules/OutputVtk.h
    OutputModules/OutputFld.h
    OutputModules/OutputStdOut.h
    OutputModules/OutputPts.h
    OutputModules/OutputXml.h
    ProcessModules/ProcessAddCompositeID.h
    ProcessModules/ProcessAddFld.h
    ProcessModules/ProcessBoundaryExtract.h
    ProcessModules/ProcessCombineAvg.h
    ProcessModules/ProcessConcatenateFld.h
    ProcessModules/ProcessDeform.h
    ProcessModules/ProcessDisplacement.h
    ProcessModules/ProcessEquiSpacedOutput.h
    ProcessModules/ProcessGrad.h
    ProcessModules/ProcessHomogeneousPlane.h
    ProcessModules/ProcessHomogeneousStretch.h
    ProcessModules/ProcessInnerProduct.h
    ProcessModules/ProcessInterpField.h
    ProcessModules/ProcessInterpPoints.h
    ProcessModules/ProcessInterpPointDataToFld.h
    ProcessModules/ProcessIsoContour.h
    ProcessModules/ProcessJacobianEnergy.h
    ProcessModules/ProcessMapping.h
    ProcessModules/ProcessNumModes.h
    ProcessModules/ProcessMeanMode.h
    ProcessModules/ProcessPointDataToFld.h
    ProcessModules/ProcessPrintFldNorms.h
    ProcessModules/ProcessScaleInFld.h
    ProcessModules/ProcessSurfDistance.h
    ProcessModules/ProcessVorticity.h
    ProcessModules/ProcessScalGrad.h
    ProcessModules/ProcessMultiShear.h
    ProcessModules/ProcessWSS.h
    ProcessModules/ProcessC0Projection.h
    ProcessModules/ProcessQualityMetric.h
)

SET(FieldUtilsSources
    Module.cpp
    Interpolator.cpp
    InputModules/InputDat.cpp
    InputModules/InputFld.cpp
    InputModules/InputXml.cpp
    InputModules/InputPts.cpp
    InputModules/InputNek5000.cpp
    OutputModules/OutputInfo.cpp
    OutputModules/OutputTecplot.cpp
    OutputModules/OutputVtk.cpp
    OutputModules/OutputFld.cpp
    OutputModules/OutputStdOut.cpp
    OutputModules/OutputPts.cpp
    OutputModules/OutputXml.cpp
    ProcessModules/ProcessAddCompositeID.cpp
    ProcessModules/ProcessAddFld.cpp
    ProcessModules/ProcessBoundaryExtract.cpp
    ProcessModules/ProcessCombineAvg.cpp
    ProcessModules/ProcessConcatenateFld.cpp
    ProcessModules/ProcessDeform.cpp
    ProcessModules/ProcessDisplacement.cpp
    ProcessModules/ProcessEquiSpacedOutput.cpp
    ProcessModules/ProcessGrad.cpp
    ProcessModules/ProcessHomogeneousPlane.cpp
    ProcessModules/ProcessHomogeneousStretch.cpp
    ProcessModules/ProcessInnerProduct.cpp
    ProcessModules/ProcessInterpField.cpp
    ProcessModules/ProcessInterpPoints.cpp
    ProcessModules/ProcessInterpPointDataToFld.cpp
    ProcessModules/ProcessIsoContour.cpp
    ProcessModules/ProcessJacobianEnergy.cpp
    ProcessModules/ProcessMapping.cpp
    ProcessModules/ProcessNumModes.cpp
    ProcessModules/ProcessMeanMode.cpp
    ProcessModules/ProcessPointDataToFld.cpp
    ProcessModules/ProcessPrintFldNorms.cpp
    ProcessModules/ProcessScaleInFld.cpp
    ProcessModules/ProcessVorticity.cpp
    ProcessModules/ProcessScalGrad.cpp
    ProcessModules/ProcessSurfDistance.cpp
    ProcessModules/ProcessMultiShear.cpp
    ProcessModules/ProcessWSS.cpp
    ProcessModules/ProcessC0Projection.cpp
    ProcessModules/ProcessQualityMetric.cpp
)

ADD_NEKTAR_LIBRARY(FieldUtils lib ${NEKTAR_LIBRARY_TYPE} ${FieldUtilsSources} ${FieldUtilsHeaders})
TARGET_LINK_LIBRARIES(FieldUtils LINK_PUBLIC GlobalMapping)

ADD_DEFINITIONS(-DFIELD_UTILS_EXPORTS)

SET_PROPERTY(SOURCE Module.cpp
    APPEND PROPERTY COMPILE_DEFINITIONS FIELDUTILS_EXTENSION_DIR="${CMAKE_INSTALL_PREFIX}/${NEKTAR_LIB_DIR}/FieldUtils")

MACRO(ADD_FIELDUTILS_PLUGIN name)
    CMAKE_PARSE_ARGUMENTS(FIELDMOD "" "" "SOURCES" ${ARGN})
    ADD_LIBRARY(${name} SHARED ${FIELDMOD_SOURCES})

    INSTALL(TARGETS ${name}
        ARCHIVE DESTINATION ${NEKTAR_LIB_DIR}/FieldUtils OPTIONAL
        LIBRARY DESTINATION ${NEKTAR_LIB_DIR}/FieldUtils OPTIONAL)

    TARGET_LINK_LIBRARIES(${name} LINK_PUBLIC FieldUtils)
ENDMACRO()

INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/FieldUtils COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
