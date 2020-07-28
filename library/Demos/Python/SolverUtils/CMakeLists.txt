IF (NEKTAR_BUILD_UNIT_TESTS)
    ADD_NEKPY_EXECUTABLE(test_solverutils_equationsystem.py test_solverutils_equationsystem.py)

    ADD_NEKPY_TEST(test_solverutils_equationsystem)
ENDIF()
