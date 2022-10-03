###############################################################################
##
## File: test_nekmesh_element.py
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
## Description: Unit tests for the Element class.
##
###############################################################################

from NekPy.NekMesh import ElmtConfig, Node, Element
from NekPy.LibUtilities import ShapeType, PointsType

import unittest

class TestElmtConfig(unittest.TestCase):
    def testElmtConfigConstructor(self):
        self.shapeType = ShapeType.Triangle
        self.order     = 2
        self.faceNodes = True
        self.volNodes  = True
        self.reorient  = 2
        self.edgeNodeType = PointsType.GaussLobattoLegendre
        self.faceNodeType = PointsType.GaussGaussLegendre
        self.config = ElmtConfig(self.shapeType, self.order,
                                 self.faceNodes, self.volNodes,
                                 self.reorient,  self.edgeNodeType,
                                 self.faceNodeType )

    def testNoDefaultConstructor(self):
        try:
            default_config = ElmtConfig()
        except:
            pass

class TestElement(unittest.TestCase):
    def setUp(self):
        self.x  = float(1.1)
        self.y  = float(2.2)
        self.z  = float(3.3)
        self.node_a = Node(1, self.x      , self.y      , self.z)
        self.node_b = Node(2, self.x + 1.0, self.y      , self.z)
        self.node_c = Node(3, self.x      , self.y + 1.0, self.z)
        self.config = ElmtConfig(ShapeType.Triangle, 1, False, False)
        self.comp_ID = 2
        self.element = Element.Create(self.config,
                                  [self.node_a, self.node_b, self.node_c],
                                  [self.comp_ID])

    def testElementGetId(self):
        # Possibly needs a better test?
        elmt_id = self.element.GetId()

    def testElementGetDim(self):
        self.assertEqual(self.element.GetDim(), 2)

    def testElementGetShapeType(self):
        self.assertEqual( self.element.GetShapeType(), ShapeType.Triangle)

    def testElementGetTag(self):
        self.assertEqual( self.element.GetTag(), "T")

if __name__ == '__main__':
    unittest.main()
