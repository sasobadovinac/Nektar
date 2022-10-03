###############################################################################
##
## File: test_nekmesh_module.py
##
## For more information, please see: http://www.nektar.info
##
## The MIT License
##
## Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
## Department of Aeronautics, Imperial College London (UK), and Scientific
## Computing and Imaging Institute, University of Utah (USA).
##
## Permission is hereby granted, free of charge, to any person obtaining a
## copy of this software and associated documentation files (the "Software"),
## to deal in the Software without restriction, including without limitation
## the rights to use, copy, modify, merge, publish, distribute, sublicense,
## and/or sell copies of the Software, and to permit persons to whom the
## Software is furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included
## in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
## OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
## THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
## FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
## DEALINGS IN THE SOFTWARE.
##
## Description: Unit tests for the Module class.
##
###############################################################################

from NekPy.LibUtilities import ShapeType
from NekPy.NekMesh import Node, Element, ElmtConfig, NodeSet, Mesh, \
                          Module, ProcessModule, ModuleType,        \
                          InputModule, OutputModule, NekMeshError
import numpy as np
import unittest
import sys


class InheritFromInputModuleTest(InputModule):
    def __init__(self, mesh):
        super(InputModule, self).__init__(mesh)

    def Process(self):
        # Call the Module functions to create all of the edges, faces and
        # composites.
        self.ProcessVertices()
        self.ProcessEdges()
        self.ProcessFaces()
        self.ProcessElements()
        self.ProcessComposites()

class InheritFromProcessModuleTest(ProcessModule):
    def __init__(self, mesh):
        super(ProcessModule, self).__init__(mesh)

    def Process(self):
        # Call the Module functions to create all of the edges, faces and
        # composites.
        self.ProcessVertices()
        self.ProcessEdges()
        self.ProcessFaces()
        self.ProcessElements()
        self.ProcessComposites()

class InheritFromOutputModuleTest(OutputModule):
    def __init__(self, mesh):
        super(OutputModule, self).__init__(mesh)

    def Process(self):
        # Call the Module functions to create all of the edges, faces and
        # composites.
        self.ProcessVertices()
        self.ProcessEdges()
        self.ProcessFaces()
        self.ProcessElements()
        self.ProcessComposites()


class TestModule(unittest.TestCase):
    def _initialize_static_values(self):
        self.coord_1x   = 0.0
        self.coord_1y   = 1.0
        self.coord_2y   = 2.0
        self.coord_2x   = 3.0
        self.nx         = 11
        self.ny         = 7
        self.comp_ID    = 2
        self.shape_type = ShapeType.Triangle
        self.x_points   = np.linspace(self.coord_1x, self.coord_2x, self.nx)
        self.y_points   = np.linspace(self.coord_1y, self.coord_2y, self.ny)
        self.expDim     = 2
        self.spaceDim   = 2
        self.nummode    = 5
        self.verbose    = True

    def _create_mesh(self):
        self.mesh = Mesh()
        self.mesh.expDim   = self.expDim
        self.mesh.spaceDim = self.spaceDim
        self.mesh.nummode  = self.nummode
        self.mesh.verbose  = self.verbose

    def _initialize_nodes(self):
        self.nodes = []
        id_cnt = 0

        for y in range(self.ny):
            tmp = []
            for x in range(self.nx):
                tmp.append(Node(id_cnt, self.x_points[x], self.y_points[y], 0.0))
                id_cnt += 1
            self.nodes.append(tmp)

    def _create_triangular_elements(self):
        for y in range(self.ny - 1):
            for x in range(self.nx - 1):
                config = ElmtConfig(ShapeType.Triangle, 1, False, False)
                self.mesh.element[2].append(
                    Element.Create(
                    config,
                    [self.nodes[y][x], self.nodes[y+1][x+1], self.nodes[y+1][x]],
                    [self.comp_ID]))
                self.mesh.element[2].append(
                    Element.Create(
                    config,
                    [self.nodes[y][x], self.nodes[y][x+1], self.nodes[y+1][x+1]],
                    [self.comp_ID]))

    def _setUpMesh(self):
        self._initialize_static_values()
        self._create_mesh()
        self._initialize_nodes()
        self._create_triangular_elements()

    def setUp(self):
        self._setUpMesh()

    def testModuleProcessRuntimeError(self):
        mod = ProcessModule(self.mesh)
        try:
            mod.Process()
        except RuntimeError:
            pass

    def testInheritFromInputModuleTest(self):
        mod = InheritFromInputModuleTest(self.mesh)
        mod.Process()

    def testInheritFromProcessModuleTest(self):
        mod = InheritFromProcessModuleTest(self.mesh)
        mod.Process()

    def testInheritFromOutputModuleTest(self):
        mod = InheritFromOutputModuleTest(self.mesh)
        mod.Process()

    def testCreateModuleUnknownParameter(self):
        mod1 = ProcessModule.Create("jac", self.mesh, unknown_parameter=False).Process()

    def testCreateExceptionUnknownModule(self):
        try:
            mod1 = ProcessModule.Create("unknown_module", self.mesh)
        except NekMeshError:
            pass

    def testCreateExceptionWrongArgs(self):
        try:
            mod1 = InputModule.Create("xml", self.mesh)
        except NekMeshError:
            pass

    def testExceptionNoMesh(self):
        try:
            mod1 = ProcessModule.Create("jac", "wrong_argument")
        except NekMeshError:
            pass

if __name__ == '__main__':
    unittest.main()
