SET(LSIAC_SOURCES
LSIACPostProcessor.cpp
MetricTensor.cpp
HandleNekMesh1D.cpp
HandleNekMesh2D.cpp
HandleNekMesh3D.cpp
GeneralBSplines.cpp
CentralBSplines.cpp
SmoothieSIAC.cpp 
SmoothieSIAC1D.cpp 
SmoothieSIAC2D.cpp 
SmoothieSIAC3D.cpp 
SIACFilter.cpp
SymmetricSIAC.cpp
NonSymmetricSIAC.cpp
OneSidedSIAC.cpp
Splines.cpp
)

SET(LSIAC_HEADERS
LSIACPostProcessor.h
MetricTensor.h
HandleMesh.h
HandleNekMesh.h
HandleNekMesh1D.h
HandleNekMesh2D.h
HandleNekMesh3D.h
Splines.h
BSplines.h
GeneralBSplines.h
CentralBSplines.h
Smoothie.h 
SmoothieSIAC.h 
SmoothieSIAC1D.h 
SmoothieSIAC2D.h 
SmoothieSIAC3D.h 
FilterType.h
SIACFilter.h 
SymmetricSIAC.h 
NonSymmetricSIAC.h
OneSidedSIAC.h
)

ADD_DEFINITIONS(-DLSIAC_EXPORTS)

ADD_NEKTAR_LIBRARY(LSIAC
    SOURCES ${LSIAC_SOURCES}
    HEADERS ${LSIAC_HEADERS}
    DEPENDS Collections
    SUMMARY "Nektar++ LSIAC library designed for postprocessing"
    DESCRIPTION "This library provides tools to apply the LSIAC filter on simulation data.")

TARGET_LINK_LIBRARIES(LSIAC LINK_PUBLIC Collections)
    
ADD_DEPENDENCIES(LSIAC eigen-3.3.7)


INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/LSIAC
    COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
