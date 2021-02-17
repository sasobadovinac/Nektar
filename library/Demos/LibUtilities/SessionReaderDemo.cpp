////////////////////////////////////////////////////////////////////////////////
//
//  File: SessionReader.cpp
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
//  Description: Check session reader functionality - currently GlobalSolve 
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>

#include <sstream>
#include <vector>

#include <boost/core/ignore_unused.hpp>

using namespace Nektar;

// Set up default values
LibUtilities::SolveConfigMap XxtMLSC = {{"TYPE","XXTMULTILEVELSTATICCOND"}};
static std::string XxtMLSCGSInfo = LibUtilities::SessionReader::
    RegisterSolveConfigInfo("XXTMULTILEVELSTATICCONDDEFAULT", XxtMLSC);

LibUtilities::SolveConfigMap XxtSC = {{"TYPE","XXTSTATICCOND"}};
static std::string XxtSCGSInfo = LibUtilities::SessionReader::
    RegisterSolveConfigInfo("XXTSTATICCONDDEFAULT", XxtSC);

LibUtilities::SolveConfigMap IterSC =
    {{"TYPE","ITERATIVESTATICCOND"},
     {"METHOD","ConjugateGradient"},
     {"PRECON","Diagonal"},
     {"ABSTOL","0.0"},
     {"RELTOL","1e-9"},
     {"MAXITER","5000"}
    };
static std::string IterSCGSInfo = LibUtilities::SessionReader::
    RegisterSolveConfigInfo("ITERATIVESTATICCONDDEFAULT", IterSC);

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session =
        LibUtilities::SessionReader::CreateInstance(argc,argv);
    session->InitSession();

    LibUtilities::SolveConfigMap
        GSmap = session->GetSolveConfigInfo("PConfig");
    
    std::cout << "PConfig: " << std::endl;
    for(auto& x: GSmap)
    {
        std::cout << "\t" << x.first << " : " << x.second <<  std::endl;
    }


    LibUtilities::SolveMap
        Smap = session->GetSolveInfo("VelocitySolve");
    std::cout << "VelocitySolve: " << std::endl;
    
    for(auto& x: Smap)
    {
        std::cout << "\t" << x.first << " : " << x.second.first << " " <<
            x.second.second << std::endl;
    }


    return 1;
}
