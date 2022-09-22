###############################################################################
##
## File: Equation_UnitTest.py
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
## Description: Unit tests for the Equation class.
##
###############################################################################

from NekPy.LibUtilities import Interpreter, Equation
import sys, unittest, time
import numpy as np

class TestEquation(unittest.TestCase):
    def testEquationEvaluate(self):
        interp = Interpreter()

        interp.AddConstant('a', 100.0)
        eq1 = Equation(interp, 'x*y+a', 'x y')

        # Takes care of floating point errors
        self.assertAlmostEqual(eq1.Evaluate(1.0, 2.0, 0, 0), 102.0)

    def testEquationGetExpression(self):
        interp = Interpreter()

        eq = Equation(interp, 'x+y+sin(x)', 'x y')

        self.assertEqual(eq.GetExpression() == 'x+y+sin(x)', True)

    def testEquationGetVlist(self):
        interp = Interpreter()

        eq = Equation(interp, 'x+y+sin(x*y*z)-t^2', 'x y z t')

        self.assertEqual(eq.GetVlist() == 'x y z t', True)

    def testEquationSetParameter(self):
        interp = Interpreter()

        eq = Equation(interp, 'x', 'x y z')
        eq.SetParameter("a", 5)
        eq1 = Equation(interp, 'a', 'x y z')

        self.assertEqual(eq1.Evaluate(0.0, 0.0, 0.0, 0.0), 5.0)

    def testEquationSetConstants(self):
        interp = Interpreter()

        eq1 = Equation(interp, 'x', 'x y z')
        eq1.SetConstants(a=23, b=34,c=41)

        self.assertEqual(interp.GetConstant("a"), 23)
        self.assertEqual(interp.GetConstant("b"), 34)
        self.assertEqual(interp.GetConstant("c"), 41)

    def testEquationEvaluate(self):
        interp = Interpreter()
        eq = Equation(interp, 'x+y+z+t', 'x y z t')

        x = np.linspace(0, 1, 100)
        y = np.linspace(0, 1, 100)
        z = np.linspace(0, 1, 100)
        t = np.linspace(0, 1, 100)
        exact = x+y+z+t

        # Checks if the solution from the equations is the same as the exact
        self.assertEqual(np.allclose(exact, eq.Evaluate(x, y, z, t)), True)

if __name__ == '__main__':
    unittest.main()
