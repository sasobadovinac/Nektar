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
LibUtilities::GlobalSolveMap XxtMLSC = {{"TYPE","XXTMULTILEVELSTATICCOND"}};
static std::string XxtMLSCGSInfo = LibUtilities::SessionReader::
    RegisterGlobalSolveInfo("XXTMULTILEVELSTATICCONDDEFAULT", XxtMLSC);

LibUtilities::GlobalSolveMap XxtSC = {{"TYPE","XXTSTATICCOND"}};
static std::string XxtSCGSInfo = LibUtilities::SessionReader::
    RegisterGlobalSolveInfo("XXTSTATICCONDDEFAULT", XxtSC);

LibUtilities::GlobalSolveMap IterSC =
    {{"TYPE","ITERATIVESTATICCOND"},
     {"METHOD","ConjugateGradient"},
     {"PRECON","Diagonal"},
     {"ABSTOL","0.0"},
     {"RELTOL","1e-9"},
     {"MAXITER","5000"}
    };
static std::string IterSCGSInfo = LibUtilities::SessionReader::
    RegisterGlobalSolveInfo("ITERATIVESTATICCONDDEFAULT", IterSC);

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session =
        LibUtilities::SessionReader::CreateInstance(argc,argv);
    session->InitSession();

    LibUtilities::GlobalSolveMap
        GSmap = session->GetGlobalSolveInfo("PConfig");
    
    for(auto& x: GSmap)
    {
        std::cout << x.first << " : " << x.second <<  std::endl;
    }
    return 1;
}
