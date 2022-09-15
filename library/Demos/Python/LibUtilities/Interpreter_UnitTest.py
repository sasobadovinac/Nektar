###############################################################################
##
## File: Interpreter_UnitTest.py
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
## Description: Unit tests for the Interpreter class.
##
###############################################################################

from NekPy.LibUtilities import Interpreter
import numpy as np
import sys, unittest, time

class TestInterpreter(unittest.TestCase):
    def setUp(self):
        self.interpreter = Interpreter()

    # Test's the GetConstant and AddConstant functions
    def testInterpreterGetConstant(self):
        self.interpreter.AddConstant("abc", 356)
        self.assertEqual(self.interpreter.GetConstant("abc"), 356)

    # Test's the GetConstants and AddConstant functions
    def testInterpreterGetConstants(self):
        self.interpreter.AddConstants(a=23, b=34, c=41)
        self.assertEqual(self.interpreter.GetConstant("b"), 34)
        self.assertEqual(self.interpreter.GetConstant("c"), 41)

    # Test's the SetParameter and GetParameter functions
    def testInterpreterSetParameter(self):
        self.interpreter.SetParameter("A", 3.14)
        self.interpreter.SetParameter("B", 45)
        self.assertEqual(self.interpreter.GetParameter("A"), 3.14)
        self.assertEqual(self.interpreter.GetParameter("B"), 45)

    # Test's the SetParameters and Getparameter functions
    def testInterpreterSetParameters(self):
        self.interpreter.SetParameters(a1=2345, a2=2354, a3=2534, a4=5234)
        self.assertEqual(self.interpreter.GetParameter("a1"), 2345)

    # Test's the Evaluate and DefineFuncion functions
    def testInterpreterEvaluate(self):
        x = np.linspace(0, np.pi, 50)
        y = np.linspace(0, np.pi, 50)
        z = np.linspace(0, np.pi, 50)
        t = np.linspace(0, 1, 50)

        # Using numpy to create an evaluate function
        numpy_func = lambda x, y, z, t: np.cos(x*y*z)*t
        numpy_sol = numpy_func(x, y, z, t)

        # Using Nektar++ to create and evaluate the same
        nekpy_func = self.interpreter.DefineFunction("x y z t", "cos(x*y*x)*t")
        nekpy_sol = self.interpreter.Evaluate(nekpy_func, x, y, z, t)

        self.assertEqual(np.allclose(nekpy_sol, numpy_sol), True)

if __name__ == '__main__':
    unittest.main()
