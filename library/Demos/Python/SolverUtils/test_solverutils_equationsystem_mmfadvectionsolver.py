from NekPy.SolverUtils    import EquationSystem
from NekPy.SolverUtils    import Plugin
from NekPy.SolverUtils    import LoadSolverPlugin
from NekPy.SpatialDomains import MeshGraph
from NekPy.LibUtilities   import SessionReader

import unittest

class TestEquationSystemMMFAdvectionSolver(unittest.TestCase):

    def getCN(self):
        return self.__class__.__name__

    @classmethod
    def setUpClass(cls):
        print("setUpClass() was called!!!")
        TestEquationSystemMMFAdvectionSolver.adrsolver_plugin = \
                                    LoadSolverPlugin('ADRSolver-plugin')

    def setUp(self):
        # load the solver
        # self.adrsolver_plugin = LoadSolverPlugin('ADRSolver-plugin')
        print("setUp() was called!!!")

    def testEquationSystemCreateMMFAdvectionSolver(self):
        self.session_name = ["test_solverutils_equationsystem_mmfadvectionsolver.py",\
                             "MMFAdvSphere.xml"]
        self.session      = SessionReader.CreateInstance(self.session_name)
        self.graph        = MeshGraph.Read(self.session)
        msg = self.getCN() + "::testEquationSystemCreateMMFAdvectionSolver: "
        try:
            eq = EquationSystem.Create("MMFAdvection",
                                        self.session, self.graph)
            eq.DoSolve()
            msg += "PASS"
            print(msg)
        except :
            msg += "FAIL"
            print(msg)
            raise

if __name__ == '__main__':
    unittest.main()
