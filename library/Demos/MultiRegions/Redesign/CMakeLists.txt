ADD_NEKTAR_EXECUTABLE(Helmholtz
    COMPONENT demos DEPENDS MultiRegions SOURCES Helmholtz.cpp)
ADD_NEKTAR_EXECUTABLE(HDGHelmholtz
    COMPONENT demos DEPENDS MultiRegions SOURCES HDGHelmholtz.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz3DHomo1D
    COMPONENT demos DEPENDS MultiRegions SOURCES Helmholtz3DHomo1D.cpp)
ADD_NEKTAR_EXECUTABLE(Helmholtz3DHomo2D
    COMPONENT demos DEPENDS MultiRegions SOURCES Helmholtz3DHomo2D.cpp)

ADD_NEKTAR_TEST(Helmholtz1D_CG_P8)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBcs_sc)
ADD_NEKTAR_TEST(Helmholtz2D_CG_P7_Modes_AllBcs_Vec_sc)
ADD_NEKTAR_TEST(Helmholtz2D_HDG_P7_Modes_AllBcs)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Homo1D)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Homo2D)
ADD_NEKTAR_TEST(Helmholtz3D_CG_Tet)
