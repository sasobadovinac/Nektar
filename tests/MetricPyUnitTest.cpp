///////////////////////////////////////////////////////////////////////////////
//
// File: MetricPyUnitTest.cpp
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
// Description: Implementation of a metric for Python's unittest.
//
///////////////////////////////////////////////////////////////////////////////

#include <MetricPyUnitTest.h>

namespace Nektar
{
std::string MetricPyUnitTest::type = GetMetricFactory().RegisterCreatorFunction(
    "PYUNITTEST", MetricPyUnitTest::create);

MetricPyUnitTest::MetricPyUnitTest(TiXmlElement *metric, bool generate)
    : MetricRegex(metric, generate)
{
    // Set up the regular expression.
    m_regex = "^(test.*) \\(.*\\) ... (ERROR|FAIL|ok).*";

    // Python's unittest framework prints to stderr by default.
    m_useStderr = true;

    // Find the functions to match against.
    TiXmlElement *func = metric->FirstChildElement("function");
    ASSERTL0(func || m_generate,
             "Missing function tag for Python unittest metric!");

    while (func)
    {
        ASSERTL0(!EmptyString(func->GetText()),
                 "Missing function name in Python unittest metric.");

        if (!m_generate)
        {
            std::vector<MetricRegexFieldValue> tmp(2);
            tmp[0] = MetricRegexFieldValue(func->GetText());
            tmp[1] = MetricRegexFieldValue("ok");
            m_matches.push_back(tmp);
        }

        func = func->NextSiblingElement("function");
    }
}

void MetricPyUnitTest::v_Generate(std::istream &pStdout, std::istream &pStderr)
{
    // Run MetricRegex to generate matches.
    MetricRegex::v_Generate(pStdout, pStderr);

    // First remove all existing values.
    m_metric->Clear();

    // Now create new values.
    for (int i = 0; i < m_matches.size(); ++i)
    {
        ASSERTL0(m_matches[i].size() == 2,
                 "Wrong number of matches for regular expression.");
        ASSERTL0(m_matches[i][1].m_value == "ok",
                 "Test " + m_matches[i][0].m_value + " has failed");

        TiXmlElement *func = new TiXmlElement("function");
        func->LinkEndChild(new TiXmlText(m_matches[i][0].m_value));
        m_metric->LinkEndChild(func);
    }
}
} // namespace Nektar
