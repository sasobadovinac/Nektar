from NekPy.SolverUtils    import EquationSystem
from NekPy.SolverUtils    import Plugin
from NekPy.SolverUtils    import LoadSolverPlugin
from NekPy.SpatialDomains import MeshGraph
from NekPy.LibUtilities   import SessionReader

import unittest

class TestEquationSystemUnsteadyViscousBurgersSolver(unittest.TestCase):

    def getCN(self):
        return self.__class__.__name__

    @classmethod
    def setUpClass(cls):
        print("setUpClass() was called!!!")
        TestEquationSystemUnsteadyViscousBurgersSolver.adrsolver_plugin = \
                                        LoadSolverPlugin('ADRSolver-plugin')

    def setUp(self):
        # load the solver
        # self.adrsolver_plugin = LoadSolverPlugin('ADRSolver-plugin')
        print("setUp() was called!!!")

    def testEquationSystemCreateUnsteadyViscousBurgersSolver(self):
        self.session_name = ["test_solverutils_equationsystem_unsteadyviscousburgerssolver.py",\
                             "UnsteadyAdvectionDiffusion_3DHomo1D_MVM.xml"]
        self.session      = SessionReader.CreateInstance(self.session_name)
        self.graph        = MeshGraph.Read(self.session)
        msg = self.getCN() + "::testEquationSystemCreateUnsteadyViscousBurgersSolver: "
        try:
            eq = EquationSystem.Create("UnsteadyViscousBurgers",
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
