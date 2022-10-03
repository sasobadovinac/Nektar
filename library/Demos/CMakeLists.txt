SUBDIRS(LibUtilities StdRegions SpatialDomains LocalRegions Collections MultiRegions MatrixFreeOps LSIAC)

IF(NEKTAR_BUILD_PYTHON)
    SUBDIRS(Python)
ENDIF()
