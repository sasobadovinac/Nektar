////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessLocalStabilityAnalysis.h
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
//  Description: Local linear stability analysis of compressible flow (for now).
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_PROCESSLOCALSTABILITYANALYSIS
#define FIELDUTILS_PROCESSLOCALSTABILITYANALYSIS

#include "ProcessBoundaryExtract.h"

//#include "/disk_two/Nek_Test/nektar++/build_f90/dist/include/nektar++/ThirdParty/CoPSE3d_LST_V1.h"


namespace Nektar
{
namespace FieldUtils
{

/**
 * @brief This processing module calculates the wall shear stress and adds it
 * as an extra-field to the output file, and writes it to a surface output file.
 */

// Define the 

extern "C" 
{
    void F77NAME(helloworld) ();
    void F77NAME(copse3d) (const long int  &option, 
                           const NekDouble* finalValue1,
                           const long int  &numStep1,
                           const NekDouble* finalValue2,
                           const long int  &numStep2);
}


class ProcessLocalStabilityAnalysis : public ProcessBoundaryExtract
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<ProcessLocalStabilityAnalysis>::AllocateSharedPtr(f);
    }
    static ModuleKey className;

    ProcessLocalStabilityAnalysis(FieldSharedPtr f);
    virtual ~ProcessLocalStabilityAnalysis();

    /// Write mesh to output file.
    virtual void Process(po::variables_map &vm);

    virtual std::string GetModuleName()
    {
        return "ProcessLocalStabilityAnalysis";
    }

    virtual std::string GetModuleDescription()
    {
        return "Calculating wall shear stress";
    }

    // test
    inline void call_hello()
    {
        F77NAME(helloworld) ();
    }


    inline void call_lst(const long int  &option, 
                         const NekDouble* finalValue1,
                         const long int  &numStep1,
                         const NekDouble* finalValue2,
                         const long int  &numStep2)
    {
        F77NAME(copse3d) (option, finalValue1, numStep1,
                                  finalValue2, numStep2);
    }

protected:

private:

};
}
}

#endif
