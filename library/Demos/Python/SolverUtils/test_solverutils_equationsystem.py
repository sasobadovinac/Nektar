from NekPy.SolverUtils    import EquationSystem
from NekPy.SolverUtils    import Plugin
from NekPy.SolverUtils    import LoadSolverPlugin
from NekPy.SpatialDomains import MeshGraph
from NekPy.LibUtilities   import SessionReader

import unittest

class TestEquationSystem(unittest.TestCase):

    def getCN(self):
        return self.__class__.__name__

    @classmethod
    def setUpClass(cls):
        print("setUpClass() was called!!!")
        TestEquationSystem.adrsolver_plugin = LoadSolverPlugin('ADRSolver-plugin')

    def setUp(self):
        # load the solver
        # self.adrsolver_plugin = LoadSolverPlugin('ADRSolver-plugin')
        print("setUp() was called!!!")

    def testEquationSystemCreateLaplaceSolver(self):
        self.session_name = ["test_solverutils_equationsystem.py",\
                             "Helmholtz3D_nodal.xml"]
        self.session      = SessionReader.CreateInstance(self.session_name)
        self.graph        = MeshGraph.Read(self.session)
        msg = self.getCN() + "::testEquationSystemCreateLaplaceSolver: "
        try:
            eq = EquationSystem.Create("Laplace", self.session, self.graph)
            eq.DoSolve()
            msg += "PASS"
            print(msg)
        except :
            msg += "FAIL"
            print(msg)
            raise

    def testEquationSystemCreatePoissonSolver(self):
        self.session_name = ["test_solverutils_equationsystem.py",\
                             "Helmholtz3D_nodal.xml"]
        self.session      = SessionReader.CreateInstance(self.session_name)
        self.graph        = MeshGraph.Read(self.session)
        msg = self.getCN() + "::testEquationSystemCreatePoissonSolver: "
        try:
            eq = EquationSystem.Create("Poisson", self.session, self.graph)
            eq.DoSolve()
            msg += "PASS"
            print(msg)
        except :
            msg += "FAIL"
            print(msg)
            raise

    def testEquationSystemCreateHelmholtzSolver(self):
        self.session_name = ["test_solverutils_equationsystem.py",\
                             "Helmholtz3D_nodal.xml"]
        self.session      = SessionReader.CreateInstance(self.session_name)
        self.graph        = MeshGraph.Read(self.session)
        msg = self.getCN() + "::testEquationSystemCreateHelmholtzSolver: "
        try:
            eq = EquationSystem.Create("Helmholtz", self.session, self.graph)
            eq.DoSolve()
            msg += "PASS"
            print(msg)
        except :
            msg += "FAIL"
            print(msg)
            raise

    def testEquationSystemCreateProjectionSolver(self):
        self.session_name = ["test_solverutils_equationsystem.py",\
                             "Helmholtz3D_nodal.xml"]
        self.session      = SessionReader.CreateInstance(self.session_name)
        self.graph        = MeshGraph.Read(self.session)
        msg = self.getCN() + "::testEquationSystemCreateProjectionSolver: "
        try:
            eq = EquationSystem.Create("Projection", self.session, self.graph)
            eq.DoSolve()
            msg += "PASS"
            print(msg)
        except :
            msg += "FAIL"
            print(msg)
            raise

    def testEquationSystemCreateSteadyAdvectionDiffusionSolver(self):
        self.session_name = ["test_solverutils_equationsystem.py",\
                             "SteadyAdvDiffReact2D_modal.xml"]
        self.session      = SessionReader.CreateInstance(self.session_name)
        self.graph        = MeshGraph.Read(self.session)
        msg = self.getCN() + "::testEquationSystemCreateSteadyAdvectionDiffusionSolver: "
        try:
            eq = EquationSystem.Create("SteadyAdvectionDiffusion",
                                        self.session, self.graph)
            eq.DoSolve()
            msg += "PASS"
            print(msg)
        except :
            msg += "FAIL"
            print(msg)
            raise

    def testEquationSystemCreateSteadyAdvectionDiffusionReactionSolver(self):
        self.session_name = ["test_solverutils_equationsystem.py",\
                             "SteadyAdvDiffReact2D_modal.xml"]
        self.session      = SessionReader.CreateInstance(self.session_name)
        self.graph        = MeshGraph.Read(self.session)
        msg = self.getCN() + "::testEquationSystemCreateSteadyAdvectionDiffusionReactionSolver: "
        try:
            eq = EquationSystem.Create("SteadyAdvectionDiffusionReaction",
                                        self.session, self.graph)
            eq.DoSolve()
            msg += "PASS"
            print(msg)
        except :
            msg += "FAIL"
            print(msg)
            raise

    def testEquationSystemCreateUnsteadyAdvectionSolver(self):
        self.session_name = ["test_solverutils_equationsystem.py",\
                             "UnsteadyAdvection_FRDG_3DHomo1D_MVM.xml"]
        self.session      = SessionReader.CreateInstance(self.session_name)
        self.graph        = MeshGraph.Read(self.session)
        msg = self.getCN() + "::testEquationSystemCreateUnsteadyAdvectionSolver: "
        try:
            eq = EquationSystem.Create("UnsteadyAdvection",
                                        self.session, self.graph)
            eq.DoSolve()
            msg += "PASS"
            print(msg)
        except :
            msg += "FAIL"
            print(msg)
            raise

    def testEquationSystemCreateMMFAdvectionSolver(self):
        self.session_name = ["test_solverutils_equationsystem.py",\
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

    def testEquationSystemCreateUnsteadyDiffusionSolver(self):
        self.session_name = ["test_solverutils_equationsystem.py",\
                             "UnsteadyDiffusion_LDG_3DHomo1D_MVM.xml"]
        self.session      = SessionReader.CreateInstance(self.session_name)
        self.graph        = MeshGraph.Read(self.session)
        msg = self.getCN() + "::testEquationSystemCreateUnsteadyDiffusionSolver: "
        try:
            eq = EquationSystem.Create("UnsteadyDiffusion",
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
