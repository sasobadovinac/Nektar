///////////////////////////////////////////////////////////////////////////////
//
// File: MetricFileExists.cpp
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
// Description: Implementation of the FileExists metric.
//
///////////////////////////////////////////////////////////////////////////////

#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/regex.hpp>

#include <MetricFileExists.h>

using namespace boost::filesystem;

namespace Nektar
{
std::string MetricFileExists::type = GetMetricFactory().RegisterCreatorFunction(
    "FILEEXISTS", MetricFileExists::create);

MetricFileExists::MetricFileExists(TiXmlElement *metric, bool generate)
    : Metric(metric, generate)
{
    TiXmlElement *file = metric->FirstChildElement("file");
    ASSERTL0(file, "Missing file tag for FileExists metric!");

    // Read metric and populate list of patterns to search for.
    while (file)
    {
        std::string pattern, count;

        // Check the pattern has been defined as an attribute and store.
        ASSERTL0(file->Attribute("pattern"), "Missing filename for file tag!");
        pattern = file->Attribute("pattern");

        // If we are testing, extract and store the expected file count
        // from the content portion of the tag.
        if (!m_generate)
        {
            count                 = file->GetText();
            m_fileCounts[pattern] = std::stoi(count);
        }
        // If we are generating, put a default value of zero so as to still
        // have the pattern in the map.
        else
        {
            m_fileCounts[pattern] = 0;
        }

        file = file->NextSiblingElement("file");
    }
}

bool MetricFileExists::v_Test(std::istream &pStdout, std::istream &pStderr)
{
    boost::ignore_unused(pStdout, pStderr);

    bool success = true;
    auto pwd     = boost::filesystem::current_path();

    // Check each pattern in turn
    for (auto it = m_fileCounts.begin(); it != m_fileCounts.end(); ++it)
    {
        int cnt = 0;
        boost::regex r(it->first.c_str());

        // Examine each file in the current path and check if it matches the
        // pattern provided. Count the number of files which match.
        for (auto &e : boost::make_iterator_range(directory_iterator(pwd), {}))
        {
            boost::cmatch matches;
            if (boost::regex_match(e.path().string().c_str(), matches, r))
            {
                if (matches.size() == 1)
                {
                    cnt++;
                }
            }
        }

        // Check if the count matches what we expect.
        if (it->second != cnt)
        {
            std::cerr << "Failed test." << std::endl;
            std::cerr << "  Expected file matches: " << it->second << std::endl;
            std::cerr << "  Found file matches:    " << cnt << std::endl;
            success = false;
        }
    }

    return success;
}

void MetricFileExists::v_Generate(std::istream &pStdout, std::istream &pStderr)
{
    boost::ignore_unused(pStdout, pStderr);

    // Update File counts.
    auto pwd = boost::filesystem::current_path();

    for (auto it = m_fileCounts.begin(); it != m_fileCounts.end(); ++it)
    {
        int cnt = 0;
        boost::regex r(it->first.c_str());
        for (auto &e : boost::make_iterator_range(directory_iterator(pwd), {}))
        {
            boost::cmatch matches;
            if (boost::regex_match(e.path().string().c_str(), matches, r))
            {
                if (matches.size() == 1)
                {
                    cnt++;
                }
            }
        }
        m_fileCounts[it->first] = cnt;
    }

    // Write new XML structure.
    TiXmlElement *file = m_metric->FirstChildElement("file");
    while (file)
    {
        std::string pattern = file->Attribute("pattern");
        file->Clear();

        ASSERTL0(m_fileCounts.count(pattern) != 0, "Couldn't find pattern " +
                                                       pattern +
                                                       " in list of calculated"
                                                       "hashes");

        file->LinkEndChild(
            new TiXmlText(std::to_string(m_fileCounts[pattern]).c_str()));
        file = file->NextSiblingElement("file");
    }
}
} // namespace Nektar
