////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputFileBase.h
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
//  Description: Base class for outputting to a file
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_OUTPUTFILEBASE
#define FIELDUTILS_OUTPUTFILEBASE

#include "../Module.h"
#include <tinyxml.h>

namespace Nektar
{
namespace FieldUtils
{

/// Converter from fld to vtk.
class OutputFileBase : public OutputModule
{
public:
    OutputFileBase(FieldSharedPtr f);
    virtual ~OutputFileBase();

protected:
    /// Write fld to output file.
    virtual void v_Process(po::variables_map &vm) override;

    virtual std::string v_GetModuleName() override
    {
        return "OutputFileBase";
    }

    virtual std::string v_GetModuleDescription() override
    {
        return "Writing file";
    }

    virtual ModulePriority v_GetModulePriority() override
    {
        return eOutput;
    }

    /// Write from pts to output file.
    virtual void v_OutputFromPts(po::variables_map &vm) = 0;

    /// Write from m_exp to output file.
    virtual void v_OutputFromExp(po::variables_map &vm) = 0;

    /// Write from data to output file.
    virtual void v_OutputFromData(po::variables_map &vm) = 0;

    virtual fs::path v_GetPath(std::string &filename, po::variables_map &vm)
    {
        boost::ignore_unused(filename, vm);
        NEKERROR(ErrorUtil::efatal, "v_GetPath not coded");
        return fs::path();
    }
    fs::path GetPath(std::string &filename, po::variables_map &vm)
    {
        return v_GetPath(filename, vm);
    }

    virtual fs::path v_GetFullOutName(std::string &filename,
                                      po::variables_map &vm)
    {
        boost::ignore_unused(filename, vm);
        NEKERROR(ErrorUtil::efatal, "v_OutputFromExp not coded");
        return fs::path();
    }
    fs::path GetFullOutName(std::string &filename, po::variables_map &vm)
    {
        return v_GetFullOutName(filename, vm);
    }

    bool m_requireEquiSpaced;

private:
    bool WriteFile(std::string &filename, po::variables_map &vm);

    void ConvertExpToEquispaced(po::variables_map &vm);

    void PrintErrorFromPts();

    void PrintErrorFromExp();
};
} // namespace FieldUtils
} // namespace Nektar

#endif
