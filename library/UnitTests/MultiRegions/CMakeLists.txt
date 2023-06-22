SET(MultiRegionsUnitTestsSources
    main.cpp
    TestVarcoeffHashing.cpp
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
ADD_NEKTAR_EXECUTABLE(MultiRegionsUnitTests
    COMPONENT unit-test DEPENDS MultiRegions SOURCES ${MultiRegionsUnitTestsSources})
ADD_TEST(NAME MultiRegionsUnitTests COMMAND MultiRegionsUnitTests --detect_memory_leaks=0)
