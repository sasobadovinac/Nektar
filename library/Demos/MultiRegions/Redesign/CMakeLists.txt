ADD_NEKTAR_EXECUTABLE(Helmholtz
    COMPONENT demos DEPENDS MultiRegions SOURCES Helmholtz.cpp)

ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBcs_sc)
