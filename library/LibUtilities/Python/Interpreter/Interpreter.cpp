///////////////////////////////////////////////////////////////////////////////
//
// File: Interpreter.cpp
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
// Description: Python wrapper for Interpreter.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Interpreter/Interpreter.h>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <boost/python/raw_function.hpp>

using namespace Nektar;
using namespace Nektar::LibUtilities;

/**
 * @brief Wrapper for Interpreter::AddConstants.
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
py::object Interpreter_AddConstants(py::tuple args, py::dict kwargs)
{
    // To extract the 'self' object
    Interpreter &interpreter = py::extract<Interpreter &>(args[0]);

    // Map that will be passed to C++
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

    interpreter.AddConstants(constants);
    return py::object();
}

/**
 * @brief Wrapper for Interpreter::SetParameters.
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
py::object Interpreter_SetParameters(py::tuple args, py::dict kwargs)
{

    // To extract the 'self' object
    Interpreter &interpreter = py::extract<Interpreter &>(args[0]);

    // Map that will be passed to C++
    std::map<std::string, NekDouble> parameters;

    // Loop over the keys inside the kwargs dictionary
    py::list keys = kwargs.keys();
    for (int i = 0; i < py::len(keys); ++i)
    {
        py::object arg = kwargs[keys[i]];
        if (arg)
        {
            parameters[py::extract<std::string>(keys[i])] =
                py::extract<NekDouble>(arg);
        }
    }
    interpreter.SetParameters(parameters);
    return py::object();
}

/**
 * @brief Wrapper for Interpreter::GetParameter.
 *
 * @param interpreter  Interpreter object
 * @param paramName    Name of parameter
 *
 * @return Parameter defined by @a paramName.
 */
NekDouble Interpreter_GetParameter(std::shared_ptr<Interpreter> interpreter,
                                   std::string paramName)
{

    return interpreter->GetParameter(paramName);
}

/**
 * @brief Wrapper for Interpreter::GetConstant.
 *
 * @param interpreter  Interpreter object
 * @param constantName Name of constant
 *
 * @return Constant defined by @a constantName.
 */
NekDouble Interpreter_GetConstant(std::shared_ptr<Interpreter> interpreter,
                                  std::string constantName)
{

    return interpreter->GetConstant(constantName);
}

/**
 * @brief Wrapper for Interpreter::Evaluate (only constant or parameter).
 *
 * @param interpreter  Interpreter object
 * @param id           id of the expression
 *
 * @return Value corresponding to @a id
 */
NekDouble Interpreter_Evaluate(std::shared_ptr<Interpreter> interpreter,
                               const int id)
{
    return interpreter->Evaluate(id);
}

/**
 * @brief Wrapper for Interpreter::Evaluate (only constant parameters).
 *
 * @param interpreter  Interpreter object
 * @param id           id of the expression
 * @param x            x-coordinate within the expression
 * @param y            y-coordinate within the expression
 * @param z            z-coordinate within the expression
 * @param t            value of time within the expression
 *
 * @return Value of evaluated expression.
 */
NekDouble Interpreter_Evaluate2(std::shared_ptr<Interpreter> interpreter,
                                const int id, const NekDouble x,
                                const NekDouble y, const NekDouble z,
                                const NekDouble t)
{
    return interpreter->Evaluate(id, x, y, z, t);
}

/**
 * @brief Wrapper for Interpreter::Evaluate (vectorised version of the
 * evaluation method that will allow the same function to be evaluated many
 * times).
 *
 * @param interpreter  Interpreter object
 * @param id           id of the expression
 * @param x            x-coordinates within the expression
 * @param y            y-coordinates within the expression
 * @param z            z-coordinates within the expression
 * @param t            values of time within the expression
 *
 * @return Values of evaluated expression.
 */
Array<OneD, NekDouble> Interpreter_Evaluate3(
    std::shared_ptr<Interpreter> interpreter, const int id,
    const Array<OneD, const NekDouble> &x,
    const Array<OneD, const NekDouble> &y,
    const Array<OneD, const NekDouble> &z,
    const Array<OneD, const NekDouble> &t)
{
    Array<OneD, NekDouble> tmp(x.size());
    interpreter->Evaluate(id, x, y, z, t, tmp);
    return tmp;
}

/**
 * @brief Wrapper for Interpreter::Evaluate (zero or more arrays).
 *
 * @param expression_id  id of the expression
 * @param points         vector containing arrays of values required for the
 *                       expression.
 *
 * @return Values of evaluated expression.
 */
Array<OneD, NekDouble> Interpreter_Evaluate4(
    std::shared_ptr<Interpreter> interpreter, const int expression_id,
    const std::vector<Array<OneD, const NekDouble>> &points)
{
    Array<OneD, NekDouble> tmp(points.size());
    interpreter->Evaluate(expression_id, points, tmp);
    return tmp;
}

void export_Interpreter()
{
    py::class_<Interpreter, std::shared_ptr<Interpreter>, boost::noncopyable>(
        "Interpreter", py::init<>())

        .def("SetRandomSeed", &Interpreter::SetRandomSeed)

        .def("AddConstants", py::raw_function(Interpreter_AddConstants))
        .def("AddConstant", &Interpreter::AddConstant)
        .def("GetConstant", Interpreter_GetConstant)

        .def("SetParameters", py::raw_function(Interpreter_SetParameters))
        .def("SetParameter", &Interpreter::SetParameter)
        .def("GetParameter", Interpreter_GetParameter)

        .def("GetTime", &Interpreter::GetTime)
        .def("DefineFunction", &Interpreter::DefineFunction)

        .def("Evaluate", Interpreter_Evaluate)
        .def("Evaluate", Interpreter_Evaluate2)
        .def("Evaluate", Interpreter_Evaluate3)
        .def("Evaluate", Interpreter_Evaluate4)

        ;
}
