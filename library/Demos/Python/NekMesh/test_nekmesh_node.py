###############################################################################
##
## File: test_nekmesh_node.py
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
## Description: Unit tests for the Node class.
##
###############################################################################

from NekPy.NekMesh import Node, NodeSet
import unittest

class TestNode(unittest.TestCase):
    def setUp(self):
        self.id = int(1)
        self.x  = float(1.1)
        self.y  = float(2.2)
        self.z  = float(3.3)
        self.node = Node(self.id, self.x, self.y, self.z)

    def testNodeConstructor(self):
        self.assertEqual(self.node.id, self.id)
        self.assertEqual(self.node.x,  self.x)
        self.assertEqual(self.node.y,  self.y)
        self.assertEqual(self.node.z,  self.z)

    def testNodeGetID(self):
        self.assertEqual(self.node.GetID(), self.id)

    def testNodeSetID(self):
        self.id = 2
        self.node.SetID(self.id)
        self.assertEqual(self.node.GetID(), self.id)

    def testNodeDistance(self):
        node = Node(self.id + 1, self.x + 1, self.y + 2, self.z + 2)
        distance = self.node.Distance(node)
        expected_distance = 3.0
        self.assertAlmostEqual(distance, expected_distance)

    def testNodeGetLoc(self):
        expected_locations = [self.x, self.y, self.z]
        locations = self.node.GetLoc()
        for i in range(len(locations)):
            self.assertEqual(locations[i], expected_locations[i])

    def testNodeAbs2(self):
        expected_abs2 = self.x**2 + self.y**2 + self.z**2
        abs2 = self.node.abs2()
        self.assertAlmostEqual(abs2, expected_abs2)

    def testNodeFieldAccess(self):
        id = int(2)
        x  = float(11.1)
        y  = float(12.2)
        z  = float(13.3)

        self.node.id = id
        self.node.x  = x
        self.node.y  = y
        self.node.z  = z
        self.assertEqual(self.node.id, id)
        self.assertEqual(self.node.x,  x)
        self.assertEqual(self.node.y,  y)
        self.assertEqual(self.node.z,  z)


class TestNodeSet(unittest.TestCase):
    def setUp(self):
        self.nodeset = NodeSet()
        self.nodeset_def_len = 10
        for i in range(self.nodeset_def_len):
            id = i
            x  = float(i)
            y  = float(self.nodeset_def_len + i)
            z  = float(self.nodeset_def_len**2 + i)
            n  = Node(id, x, y, z)
            self.nodeset.add(n)

    def testNodeSet__len__(self):
        slen = len(self.nodeset)
        self.assertEqual(slen, self.nodeset_def_len)

    def testNodeSetClear(self):
        self.nodeset.clear()
        slen = len(self.nodeset)
        self.assertEqual(slen, 0)

    def testNodeSet__iter__(self):
        for node in self.nodeset:
            self.assertTrue(node.GetID() <= self.nodeset_def_len)

    def testNodeSet__contains__(self):
        for node in self.nodeset:
            self.assertTrue(node in self.nodeset)

    def testNodeSetAdd(self):
        new_node = Node(len(self.nodeset) + 1, 1.0, 2.0, 3.0)
        self.nodeset.add(new_node)
        self.assertTrue(new_node in self.nodeset)

if __name__ == '__main__':
    unittest.main()
