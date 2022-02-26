///////////////////////////////////////////////////////////////////////////////
//
// File: Timer.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
// Description: Time getting class
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Communication/CommSerial.h>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <iomanip>
#include <tuple>

namespace Nektar
{
namespace LibUtilities
{

void Timer::Start()
{
    ASSERTL0(!m_isactive, "Call to Timer::Start() done when timer is active.");
    m_isactive = true;
    m_start = Clock::now();
}

void Timer::Stop()
{
    m_end = Clock::now();
    ASSERTL0(m_isactive, "Call to Timer::Stop() done when timer is inactive.");
    m_isactive = false;
}

Timer::Seconds Timer::Elapsed()
{
    ASSERTL0(!m_isactive, 
        "Call to Timer::Elapsed() done before Timer::Stop().");
    return std::chrono::duration_cast<Seconds>(m_end - m_start);
}

NekDouble Timer::TimePerTest(unsigned int n)
{
    return Elapsed().count() / static_cast<NekDouble>(n);
}

void Timer::AccumulateRegion(std::string region, int iolevel)
{
    // search for region
    auto search = m_elapsedRegion.find(region);
    if (search == m_elapsedRegion.end())
    {
        m_elapsedRegion.insert({region,
         std::make_tuple<Timer::Seconds, size_t>(this->Elapsed(),1, iolevel)});
    }
    else
    {
        std::get<0>(search->second) += this->Elapsed();
        std::get<1>(search->second) += 1;
    }
}

void Timer::PrintElapsedRegions()
{
    std::string  def("default");
    char *argv = new char [def.length()+1];
    std::strcpy(argv,def.c_str());
    LibUtilities::CommSharedPtr comm = 
        MemoryManager<LibUtilities::CommSerial>:: AllocateSharedPtr(1,&argv);

    PrintElapsedRegions(comm);
}
    
void Timer::PrintElapsedRegions(LibUtilities::CommSharedPtr comm,
                                std::ostream &o,
                                int iolevel)
{
    std::vector<std::string> labels{
        "Region",
        "Elapsed time Avg (s)",
        "Min (s)",
        "Max (s)",
        "Count"};

    // Set width of each column (minimum 10 characters)
    std::vector<size_t> widths;
    for(const auto& label : labels){
        widths.push_back(std::max<size_t>(label.size()+2, 10));
    }
    // Make sure that names for each "Region" fits
    for(const auto& entry : m_elapsedRegion){
        widths[0] = std::max<size_t>(entry.first.size()+2, widths[0]);
    }

    // Print header with labels
    if (comm->GetRank() == 0 &&
        m_elapsedRegion.begin() != m_elapsedRegion.end())
    {
        o << "-------------------------------------------\n";
        for(int i=0; i<labels.size(); ++i){
            o << std::setw(widths[i]) << labels[i];
        }
        o << '\n';
    }
    
    // first write out execute time
    auto item = m_elapsedRegion.find("Execute");
    if(item != m_elapsedRegion.end())
    {
        auto elapsedAve = std::get<0>(item->second).count();
        comm->AllReduce(elapsedAve, LibUtilities::ReduceSum);
        elapsedAve /= comm->GetSize();
        auto elapsedMin = std::get<0>(item->second).count();
        comm->AllReduce(elapsedMin, LibUtilities::ReduceMin);
        auto elapsedMax = std::get<0>(item->second).count();
        comm->AllReduce(elapsedMax, LibUtilities::ReduceMax);
        
        if (comm->GetRank() == 0)
        {
            o << std::setw(widths[0]) << item->first
              << std::setw(widths[1]) << elapsedAve
              << std::setw(widths[2]) << elapsedMin
              << std::setw(widths[3]) << elapsedMax
              << std::setw(widths[4]) << std::get<1>(item->second) << '\n';
        }
    }            

    // write out all other timings order alphabetically on string
    for (auto item = m_elapsedRegion.begin();
            item != m_elapsedRegion.end(); ++item)
    {
        if(std::get<2>(item->second) < iolevel)
        {
            if(boost::iequals(item->first,"Execute"))
            {
                continue;
            }

            auto elapsedAve = std::get<0>(item->second).count();
            comm->AllReduce(elapsedAve, LibUtilities::ReduceSum);
            elapsedAve /= comm->GetSize();
            auto elapsedMin = std::get<0>(item->second).count();
            comm->AllReduce(elapsedMin, LibUtilities::ReduceMin);
            auto elapsedMax = std::get<0>(item->second).count();
            comm->AllReduce(elapsedMax, LibUtilities::ReduceMax);

            if (comm->GetRank() == 0)
            {
                o << std::setw(widths[0]) << item->first
                  << std::setw(widths[1]) << elapsedAve
                  << std::setw(widths[2]) << elapsedMin
                  << std::setw(widths[3]) << elapsedMax
                  << std::setw(widths[4]) << std::get<1>(item->second) << '\n';
            }
        }
    }
}
// static members init
std::map<std::string, std::tuple<Timer::Seconds, size_t, int>>
    Timer::m_elapsedRegion{};

}
} // end Nektar namespace
