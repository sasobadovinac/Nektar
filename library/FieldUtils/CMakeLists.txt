SET(FieldUtilsHeaders
    Module.h
    Field.hpp
    Interpolator.h
    Octree.h
    InputModules/InputDat.h
    InputModules/InputFld.h
    InputModules/InputXml.h
    InputModules/InputPts.h
    InputModules/InputNek5000.h
    InputModules/InputSemtex.h
    OutputModules/OutputFileBase.h
    OutputModules/OutputInfo.h
    OutputModules/OutputTecplot.h
    OutputModules/OutputVtkBase.h
    OutputModules/OutputFld.h
    OutputModules/OutputStdOut.h
    OutputModules/OutputPts.h
    OutputModules/OutputXml.h
    ProcessModules/ProcessAddCompositeID.h
    ProcessModules/ProcessFieldFromString.h
    ProcessModules/ProcessAddFld.h
    ProcessModules/ProcessBoundaryExtract.h
    ProcessModules/ProcessCFL.h
    ProcessModules/ProcessCombineAvg.h
    ProcessModules/ProcessConcatenateFld.h
    ProcessModules/ProcessCreateExp.h
    ProcessModules/ProcessDeform.h
    ProcessModules/ProcessDisplacement.h
    ProcessModules/ProcessDOF.h
    ProcessModules/ProcessEquiSpacedOutput.h
    ProcessModules/ProcessGrad.h
    ProcessModules/ProcessHalfModeToFourier.h
    ProcessModules/ProcessHomogeneousPlane.h
    ProcessModules/ProcessHomogeneousStretch.h
    ProcessModules/ProcessInnerProduct.h
    ProcessModules/ProcessInterpField.h
    ProcessModules/ProcessInterpPoints.h
    ProcessModules/ProcessInterpPointDataToFld.h
    ProcessModules/ProcessInterpPtsToPts.h
    ProcessModules/ProcessIsoContour.h
    ProcessModules/ProcessJacobianEnergy.h
    ProcessModules/ProcessL2Criterion.h
    ProcessModules/ProcessMapping.h
    ProcessModules/ProcessNumModes.h
    ProcessModules/ProcessMean.h
    ProcessModules/ProcessMeanMode.h
    ProcessModules/ProcessPhiFromFile.h
    ProcessModules/ProcessPointDataToFld.h
    ProcessModules/ProcessPrintFldNorms.h
    ProcessModules/ProcessRemoveField.h
    ProcessModules/ProcessScaleInFld.h
    ProcessModules/ProcessStreamFunction.h
    ProcessModules/ProcessSurfDistance.h
    ProcessModules/ProcessVorticity.h
    ProcessModules/ProcessScalGrad.h
    ProcessModules/ProcessMultiShear.h
    ProcessModules/ProcessWSS.h
    ProcessModules/ProcessWallNormalData.h
    ProcessModules/ProcessC0Projection.h
    ProcessModules/ProcessQCriterion.h
    ProcessModules/ProcessQualityMetric.h
)

SET(FieldUtilsSources
    Module.cpp
    Interpolator.cpp
    Octree.cpp
    InputModules/InputDat.cpp
    InputModules/InputFld.cpp
    InputModules/InputXml.cpp
    InputModules/InputPts.cpp
    InputModules/InputNek5000.cpp
    InputModules/InputSemtex.cpp
    OutputModules/OutputFileBase.cpp
    OutputModules/OutputInfo.cpp
    OutputModules/OutputTecplot.cpp
    OutputModules/OutputVtkBase.cpp
    OutputModules/OutputFld.cpp
    OutputModules/OutputStdOut.cpp
    OutputModules/OutputPts.cpp
    OutputModules/OutputXml.cpp
    ProcessModules/ProcessAddCompositeID.cpp
    ProcessModules/ProcessFieldFromString.cpp
    ProcessModules/ProcessAddFld.cpp
    ProcessModules/ProcessBoundaryExtract.cpp
    ProcessModules/ProcessCFL.cpp
    ProcessModules/ProcessCombineAvg.cpp
    ProcessModules/ProcessConcatenateFld.cpp
    ProcessModules/ProcessCreateExp.cpp
    ProcessModules/ProcessDeform.cpp
    ProcessModules/ProcessDisplacement.cpp
    ProcessModules/ProcessDOF.cpp
    ProcessModules/ProcessEquiSpacedOutput.cpp
    ProcessModules/ProcessGrad.cpp
    ProcessModules/ProcessHalfModeToFourier.cpp
    ProcessModules/ProcessHomogeneousPlane.cpp
    ProcessModules/ProcessHomogeneousStretch.cpp
    ProcessModules/ProcessInnerProduct.cpp
    ProcessModules/ProcessInterpField.cpp
    ProcessModules/ProcessInterpPoints.cpp
    ProcessModules/ProcessInterpPointDataToFld.cpp
    ProcessModules/ProcessInterpPtsToPts.cpp
    ProcessModules/ProcessIsoContour.cpp
    ProcessModules/ProcessJacobianEnergy.cpp
    ProcessModules/ProcessL2Criterion.cpp
    ProcessModules/ProcessMapping.cpp
    ProcessModules/ProcessNumModes.cpp
    ProcessModules/ProcessMean.cpp
    ProcessModules/ProcessMeanMode.cpp
    ProcessModules/ProcessPhiFromFile.cpp
    ProcessModules/ProcessPointDataToFld.cpp
    ProcessModules/ProcessPrintFldNorms.cpp
    ProcessModules/ProcessRemoveField.cpp    
    ProcessModules/ProcessScaleInFld.cpp
    ProcessModules/ProcessVorticity.cpp
    ProcessModules/ProcessScalGrad.cpp
    ProcessModules/ProcessStreamFunction.cpp
    ProcessModules/ProcessSurfDistance.cpp
    ProcessModules/ProcessMultiShear.cpp
    ProcessModules/ProcessWSS.cpp
    ProcessModules/ProcessWallNormalData.cpp
    ProcessModules/ProcessC0Projection.cpp
    ProcessModules/ProcessQCriterion.cpp
    ProcessModules/ProcessQualityMetric.cpp
)

IF (NEKTAR_USE_VTK)
    SET(FieldUtilsSources ${FieldUtilsSources}
            OutputModules/OutputVtk.cpp)
    SET(FieldUtilsHeaders ${FieldUtilsHeaders}
            OutputModules/OutputVtk.h)
ENDIF()

ADD_NEKTAR_LIBRARY(FieldUtils
    SUMMARY "Nektar++ FieldUtils library"
    DESCRIPTION "Collection of post-processing modules for the FieldConvert utility and filter"
    SOURCES ${FieldUtilsSources}
    HEADERS ${FieldUtilsHeaders}
    DEPENDS GlobalMapping)

ADD_DEFINITIONS(-DFIELD_UTILS_EXPORTS)

IF (NEKTAR_USE_VTK)
    IF (VTK_MAJOR_VERSION LESS 6)
        TARGET_LINK_LIBRARIES(FieldUtils LINK_PUBLIC vtkCommon vtkIO vtkGraphics vtkFiltersCore)
    ELSEIF (VTK_MAJOR_VERSION  GREATER_EQUAL 9)
        TARGET_LINK_LIBRARIES(FieldUtils LINK_PUBLIC VTK::CommonCore VTK::IOLegacy VTK::IOXML VTK::FiltersCore)
        vtk_module_autoinit(
                TARGETS FieldUtils
                MODULES VTK::CommonCore VTK::IOLegacy VTK::IOXML VTK::FiltersCore)
    ELSE ()
        TARGET_LINK_LIBRARIES(FieldUtils LINK_PUBLIC vtkCommonCore vtkIOLegacy vtkIOXML vtkFiltersCore)
    ENDIF ()
ENDIF (NEKTAR_USE_VTK)

IF (NEKTAR_BUILD_PYTHON)
    SUBDIRS(Python)
ENDIF()

INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/FieldUtils COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
