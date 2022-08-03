ADD_NEKTAR_EXECUTABLE(Helmholtz
    COMPONENT demos DEPENDS MultiRegions SOURCES Helmholtz.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz
    COMPONENT demos DEPENDS MultiRegions SOURCES HDGHelmholtz.cpp)

ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBcs_sc)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBcs_Vec_sc)
ADD_NEKTAR_TEST(Helmholtz2D_HDG_P7_Modes_AllBcs)
