////////////////////////////////////////////////////////////////////////////////
//
//  File:  Equation.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:  Wrapper to Interpreter class.
//
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBUTILITIES_BASICUTILS_EQUATION_H
#define NEKTAR_LIBUTILITIES_BASICUTILS_EQUATION_H

#include <map>
#include <string>

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Interpreter/Interpreter.h>

namespace Nektar
{
namespace LibUtilities
{

/**
 * @class Equation
 *
 * This class stores a string form of a symbolic expression to be evaluated
 * e.g. for the boundary conditions, the unique numeric ID of that expression
 * and a reference to the unique static instance of the Interpreter.
 *
 * The scenario is that for multiple copies of Equation class holding their
 * symbolic expressions in the std::string form, there is a unique instance of
 * Interpreter, which holds a set of pre-processed symbolic expressions in the
 * form of sequential containers of execution functors, ready for fast
 * evaluation.
 *
 * The interpreter also keeps all constants and parameters specified in an XML
 * file. There should be only one copy of Equation class per each symbolic
 * expression specified in XML file, modulo possible bugs. Classes Equation and
 * Interpreter are complementary in a sense that the expression ID stored in the
 * Equation class is generated by the Interpreter which holds ordered container
 * of pre-processed expressions.
 */
class Equation
{
public:
    LIB_UTILITIES_EXPORT Equation(const Equation &) = default;
    LIB_UTILITIES_EXPORT Equation(InterpreterSharedPtr evaluator,
                                  const std::string &expr  = "",
                                  const std::string &vlist = "");

    LIB_UTILITIES_EXPORT Equation &operator=(const Equation &src);

    LIB_UTILITIES_EXPORT NekDouble Evaluate() const;

    LIB_UTILITIES_EXPORT NekDouble Evaluate(NekDouble x, NekDouble y = 0,
                                            NekDouble z = 0,
                                            NekDouble t = 0) const;

    LIB_UTILITIES_EXPORT void Evaluate(const Array<OneD, const NekDouble> &x,
                                       const Array<OneD, const NekDouble> &y,
                                       const Array<OneD, const NekDouble> &z,
                                       Array<OneD, NekDouble> &result) const;

    LIB_UTILITIES_EXPORT void Evaluate(const Array<OneD, const NekDouble> &x,
                                       const Array<OneD, const NekDouble> &y,
                                       const Array<OneD, const NekDouble> &z,
                                       const NekDouble t,
                                       Array<OneD, NekDouble> &result) const;

    LIB_UTILITIES_EXPORT void Evaluate(const Array<OneD, const NekDouble> &x,
                                       const Array<OneD, const NekDouble> &y,
                                       const Array<OneD, const NekDouble> &z,
                                       const Array<OneD, const NekDouble> &t,
                                       Array<OneD, NekDouble> &result) const;

    LIB_UTILITIES_EXPORT void Evaluate(
        const std::vector<Array<OneD, const NekDouble>> points,
        Array<OneD, NekDouble> &result) const;

    LIB_UTILITIES_EXPORT void SetParameter(const std::string &name,
                                           NekDouble value);

    LIB_UTILITIES_EXPORT void SetConstants(
        const std::map<std::string, NekDouble> &constants);

    LIB_UTILITIES_EXPORT std::string GetExpression(void) const;

    LIB_UTILITIES_EXPORT std::string GetVlist(void) const;

    /// Returns time spend on expression evaluation at points (it does not
    /// include parse/pre-processing time).
    LIB_UTILITIES_EXPORT NekDouble GetTime() const;

private:
    std::string m_vlist;
    std::string m_expr;
    int m_expr_id;
    InterpreterSharedPtr m_evaluator;
};

typedef std::shared_ptr<Equation> EquationSharedPtr;

} // namespace LibUtilities
} // namespace Nektar

#endif // NEKTAR_LIBUTILITIES_EQUATION_HPP
