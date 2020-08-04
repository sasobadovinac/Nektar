from NekPy.SolverUtils    import EquationSystem
from NekPy.SpatialDomains import MeshGraph
from NekPy.LibUtilities   import SessionReader

import unittest

class TestEquationSystem(unittest.TestCase):

    def getCN(self):
        return self.__class__.__name__

    def setUp(self):
        self.session_name = ["test_solverutils_equationsystem.py",\
                             "newsquare_2x2.xml"]
        self.session      = SessionReader.CreateInstance(self.session_name)
        self.graph        = MeshGraph.Read(self.session)

    def testEquationSystemConstructor(self):
        msg = self.getCN() + "::testEquationSystemConstructor: "
        try:
            eq1 = EquationSystem(self.session, self.graph)
            eq2 = EquationSystem.Create(self.session, self.graph)
            msg += "PASS"
            print(msg)
        except :
            msg += "FAIL"
            print(msg)
            raise


if __name__ == '__main__':
    unittest.main()
