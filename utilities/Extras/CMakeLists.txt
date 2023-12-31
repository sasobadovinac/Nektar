ADD_UTILITIES_EXECUTABLE(CheckXmlFile COMPONENT utilities-extra SOURCES CheckXmlFile.cpp)

IF (NEKTAR_USE_VTK)
    ADD_UTILITIES_EXECUTABLE(VtkToFld COMPONENT utilities-extra SOURCES VtkToFld.cpp)
    ADD_UTILITIES_EXECUTABLE(VtkStripsToPolys COMPONENT utilities-extra SOURCES VtkStripsToPolys.cpp)

    IF (VTK_MAJOR_VERSION LESS 6)
        TARGET_LINK_LIBRARIES(VtkToFld LINK_PUBLIC vtkCommon vtkIO vtkGraphics)
        TARGET_LINK_LIBRARIES(VtkStripsToPolys LINK_PUBLIC vtkCommon vtkIO vtkGraphics)
    ELSEIF (VTK_MAJOR_VERSION  GREATER_EQUAL 9)
        TARGET_LINK_LIBRARIES(VtkToFld LINK_PUBLIC VTK::CommonCore VTK::IOLegacy VTK::IOXML
                              VTK::FiltersCore )
        vtk_module_autoinit(
            TARGETS VtkToFld
            MODULES VTK::CommonCore VTK::IOLegacy VTK::IOXML VTK::FiltersCore)
        TARGET_LINK_LIBRARIES(VtkStripsToPolys LINK_PUBLIC VTK::CommonCore VTK::IOLegacy VTK::IOXML)
        vtk_module_autoinit(
            TARGETS VtkStripsToPolys
            MODULES VTK::CommonCore VTK::IOLegacy VTK::IOXML)
    ELSE ()
        TARGET_LINK_LIBRARIES(VtkToFld LINK_PUBLIC vtkCommonCore vtkIOLegacy
                                              vtkFiltersGeometry vtkFiltersCore)
        TARGET_LINK_LIBRARIES(VtkStripsToPolys LINK_PUBLIC vtkCommonCore vtkIOLegacy)
    ENDIF ()

    ADD_UTILITIES_EXECUTABLE(VtkToPng COMPONENT utilities-extra SOURCES VtkToPng.cpp)
    IF (VTK_MAJOR_VERSION LESS 6)
        TARGET_LINK_LIBRARIES(VtkToPng LINK_PUBLIC vtkCommon vtksys vtkViews vtkWidgets
                                                   vtkRendering vtkIO)
    ELSEIF (VTK_MAJOR_VERSION  GREATER_EQUAL 9)
        TARGET_LINK_LIBRARIES(VtkToPng LINK_PUBLIC VTK::CommonCore VTK::RenderingCore
                              VTK::IOXML VTK::IOImage)
        vtk_module_autoinit(
            TARGETS VtkToPng
            MODULES VTK::CommonCore VTK::RenderingCore VTK::IOXML  VTK::IOImage) 
        ELSE ()
        TARGET_LINK_LIBRARIES(VtkToPng LINK_PUBLIC ${VTK_LIBRARIES})
    ENDIF ()
ENDIF (NEKTAR_USE_VTK)
