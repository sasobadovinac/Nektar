ADD_NEKTAR_EXECUTABLE(PartitionAnalyse
    COMPONENT demos DEPENDS SpatialDomains SOURCES PartitionAnalyse.cpp)
ADD_NEKTAR_EXECUTABLE(InterfaceTest
        COMPONENT demos DEPENDS SpatialDomains SOURCES InterfaceTest.cpp)

ADD_NEKTAR_EXECUTABLE(InterfaceDuplicateEdge
        COMPONENT demos DEPENDS MultiRegions SOURCES InterfaceDuplicateEdge.cpp)