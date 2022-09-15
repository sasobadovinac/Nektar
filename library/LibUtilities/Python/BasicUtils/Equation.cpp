///////////////////////////////////////////////////////////////////////////////
//
// File: Equation.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Python wrapper for Equation.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/Interpreter/Interpreter.h>
#include <LibUtilities/Python/NekPyConfig.hpp>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

#include <boost/python/raw_function.hpp>

using namespace Nektar;
using namespace Nektar::LibUtilities;

/**
 * @brief Wrapper for Equation::SetConstants.
 *
 * boost.python does not know how (by default) to convert from any Python
 * datatype to a C++ map, so we add a pythonic way to set parameters through
 * keyword arguments.
 *
 * @param args    Python function positional arguments.
 * @param kwargs  Python function keyword arguments.
 *
 * @return None (null py::object).
 */
py::object Equation_SetConstants(py::tuple args, py::dict kwargs)
{
    // To extract the 'self' object
    Equation &equation = py::extract<Equation &>(args[0]);

    // Map that will be passed to C++.
    std::map<std::string, NekDouble> constants;

    // Loop over the keys inside the kwargs dictionary
    py::list keys = kwargs.keys();
    for (int i = 0; i < py::len(keys); ++i)
    {
        py::object arg = kwargs[keys[i]];
        if (arg)
        {
            constants[py::extract<std::string>(keys[i])] =
                py::extract<NekDouble>(arg);
        }
    }

    equation.SetConstants(constants);

    // Return nothing.
    return py::object();
}

/**
 * @brief Construct an equation object from an expression string @p expr and a
 * list of variables @p vlist.
 *
 * @param evaluator   Interpreter object
 * @param expr        String contaning the expression
 * @param vlist       String contining the list of variables in @p expr.
 *
 * @return An #Equation object.
 */
std::shared_ptr<Equation> ConstructEquation(
    std::shared_ptr<Interpreter> evaluator, std::string expr, std::string vlist)
{
    return std::make_shared<Equation>(evaluator, expr, vlist);
}

/**
 * @brief Wrapper for Equation::Evaluate (overloaded for no parameters).
 *
 * @param equation  Equation object
 *
 * @return Value of @p equation.
 */
NekDouble Equation_Evaluate1(std::shared_ptr<Equation> equation)
{
    return equation->Evaluate();
}

/**
 * @brief Wrapper for Equation::Evaluate (overloaded for constant parameters).
 *
 * @param equation     Equation object from Python
 * @param x            x-coordinate within the expression
 * @param y            y-coordinate within the expression
 * @param z            z-coordinate within the expression
 * @param t            value of time within the expression
 *
 * @return Value of evaluated expression.
 */
NekDouble Equation_Evaluate2(std::shared_ptr<Equation> equation,
                             const NekDouble x, const NekDouble y = 0,
                             const NekDouble z = 0, const NekDouble t = 0)
{
    return equation->Evaluate(x, y, z, t);
}

/**
 * @brief Wrapper for Equation::Evaluate (overloaded for Array parameters).
 *
 * @param equation     Equation object from Python
 * @param x            x-coordinates within the expression
 * @param y            y-coordinates within the expression
 * @param z            z-coordinates within the expression
 *
 * @return Array containing values of evaluated expression.
 */
Array<OneD, NekDouble> Equation_Evaluate3(std::shared_ptr<Equation> equation,
                                          const Array<OneD, const NekDouble> &x,
                                          const Array<OneD, const NekDouble> &y,
                                          const Array<OneD, const NekDouble> &z)
{
    Array<OneD, NekDouble> tmp(x.size());
    equation->Evaluate(x, y, z, tmp);
    return tmp;
}

/**
 * @brief Wrapper for Equation::Evaluate (overloaded for Array parameters +
 * constant time).
 *
 * @param equation     Equation object from Python
 * @param x            x-coordinates within the expression
 * @param y            y-coordinates within the expression
 * @param z            z-coordinates within the expression
 * @param t            Value of time
 *
 * @return Array containing values of evaluated expression.
 */
Array<OneD, NekDouble> Equation_Evaluate4(std::shared_ptr<Equation> equation,
                                          const Array<OneD, const NekDouble> &x,
                                          const Array<OneD, const NekDouble> &y,
                                          const Array<OneD, const NekDouble> &z,
                                          const NekDouble t)
{
    Array<OneD, NekDouble> tmp(x.size());
    equation->Evaluate(x, y, z, t, tmp);
    return tmp;
}

/**
 * @brief Wrapper for Equation::Evaluate (overloaded for Array parameters +
 * Array time).
 *
 * @param equation     Equation object from Python
 * @param x            x-coordinates within the expression
 * @param y            y-coordinates within the expression
 * @param z            z-coordinates within the expression
 * @param t            Time values within the expression.
 *
 * @return Array containing values of evaluated expression.
 */
Array<OneD, NekDouble> Equation_Evaluate5(std::shared_ptr<Equation> equation,
                                          const Array<OneD, const NekDouble> &x,
                                          const Array<OneD, const NekDouble> &y,
                                          const Array<OneD, const NekDouble> &z,
                                          const Array<OneD, const NekDouble> &t)
{
    Array<OneD, NekDouble> tmp(x.size());
    equation->Evaluate(x, y, z, t, tmp);
    return tmp;
}

void export_Equation()
{
    py::class_<Equation, std::shared_ptr<Equation>, boost::noncopyable>(
        "Equation", py::no_init)

        .def("__init__", py::make_constructor(ConstructEquation))

        .def("Evaluate", Equation_Evaluate1)
        .def("Evaluate", Equation_Evaluate2)
        .def("Evaluate", Equation_Evaluate3)
        .def("Evaluate", Equation_Evaluate4)
        .def("Evaluate", Equation_Evaluate5)

        .def("SetParameter", &Equation::SetParameter)
        .def("SetConstants", py::raw_function(Equation_SetConstants))
        .def("GetExpression", &Equation::GetExpression)
        .def("GetVlist", &Equation::GetVlist)

        .def("GetTime", &Equation::GetTime);
}
