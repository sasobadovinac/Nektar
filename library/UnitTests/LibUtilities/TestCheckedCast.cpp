///////////////////////////////////////////////////////////////////////////////
//
// File: TestCheckedCast.cpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/CheckedCast.hpp>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>

namespace Nektar
{
namespace LibUtilities
{
namespace CheckCastUnitTest
{

BOOST_AUTO_TEST_CASE(TestDoubleToInt)
{
    // expecting to convert
    {
        double adouble = std::numeric_limits<int>::max();
        int aint       = checked_cast<int>(adouble);
        BOOST_CHECK_EQUAL(aint, adouble);
    }
    {
        double adouble = std::numeric_limits<int>::min();
        int aint       = checked_cast<int>(adouble);
        BOOST_CHECK_EQUAL(aint, adouble);
    }

    // expecting to fail and throw
    try
    {
        double adouble = std::numeric_limits<int>::max() + 1.0;
        int aint       = checked_cast<int>(adouble);
        BOOST_CHECK_EQUAL(aint, adouble);
    }
    catch (std::runtime_error &e)
    {
        std::string errmss = e.what();
        BOOST_CHECK_EQUAL("Level 0 assertion violation", errmss.substr(0, 27));
    }

    try
    {
        double adouble = std::numeric_limits<int>::min() - 1.0;
        int aint       = checked_cast<int>(adouble);
        BOOST_CHECK_EQUAL(aint, adouble);
    }
    catch (std::runtime_error &e)
    {
        std::string errmss = e.what();
        BOOST_CHECK_EQUAL("Level 0 assertion violation", errmss.substr(0, 27));
    }
}

} // namespace CheckCastUnitTest
} // namespace LibUtilities
} // namespace Nektar