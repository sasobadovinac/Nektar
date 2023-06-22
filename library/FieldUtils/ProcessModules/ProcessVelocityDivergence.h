////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessVelocityDivergence.h
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
//  Description: Computes the divergence of the velocity field.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_PROCESSVELOCITYDIVERGENCE
#define FIELDUTILS_PROCESSVELOCITYDIVERGENCE

#include "../Module.h"

namespace Nektar
{
namespace FieldUtils
{
/**
 * @brief This processing module calculates the divergence of the velocity
 * field and adds it as an extra-field to the output file.
 */
class ProcessVelocityDivergence : public ProcessModule
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<ProcessVelocityDivergence>::AllocateSharedPtr(f);
    }
    static ModuleKey className;

    ProcessVelocityDivergence(FieldSharedPtr f);
    virtual ~ProcessVelocityDivergence();

protected:
    /// Write mesh to output file.
    virtual void v_Process(po::variables_map &vm) override;

    virtual std::string v_GetModuleName() override
    {
        return "ProcessVelocityDivergence";
    }

    virtual std::string v_GetModuleDescription() override
    {
        return "Calculating velocity divergence";
    }

    virtual ModulePriority v_GetModulePriority() override
    {
        return eModifyExp;
    }

    void GetVelocity(Array<OneD, Array<OneD, NekDouble>> &vel, int strip = 0);

private:
    int m_spacedim;
};
} // namespace FieldUtils
} // namespace Nektar

#endif
