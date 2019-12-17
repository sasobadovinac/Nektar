SUBDIRS(LibUtilities StdRegions SpatialDomains LocalRegions Collections MultiRegions LSIAC)

IF(NEKTAR_BUILD_PYTHON)
    SUBDIRS(Python)
ENDIF()
