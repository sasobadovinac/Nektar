////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessLSIAC.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
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
//  Description: Computes LSIAC field.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_PROCESSLSIAC
#define FIELDUTILS_PROCESSLSIAC

#include "../Module.h"
#include <tinyxml.h>

namespace Nektar
{
namespace FieldUtils
{
struct LSIACParams
{
    vector<string> inVar;
    vector<string> outVar;
    Array<OneD, NekDouble> direction;
    int derivative;
    int order;
    double maxChLen;
    int filterType;
};
/**
 * @brief This processing module calculates the LSIAC and adds it
 * as an extra-field to the output file.
 */
class ProcessLSIAC : public ProcessModule
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<ProcessLSIAC>::AllocateSharedPtr(f);
    }
    static ModuleKey className;

    ProcessLSIAC(FieldSharedPtr f);
    virtual ~ProcessLSIAC();

    void XmlType(TiXmlNode *pParent);
    /// Write mesh to output file.
    virtual void Process(po::variables_map &vm);

    virtual std::string GetModuleName()
    {
        return "ProcessLSIAC";
    }

    virtual std::string GetModuleDescription()
    {
        return "Calculating LSIAC";
    }

    virtual ModulePriority GetModulePriority()
    {
        return eModifyExp;
    }

    // protected:
    //  void GetVelocity( Array<OneD, Array<OneD, NekDouble> > &vel, int strip =
    //  0);

private:
    int m_spacedim;
    vector<string> getVecFromStr(string str);
    Array<OneD, NekDouble> getDirFromStr(string str);
    string getStrFromNode(TiXmlNode *pChild, string str);
    void ApplyLSIAC(LSIACParams set1);
    bool getIntVFromVals(vector<string> InVar, vector<int> &outV);
    void printHelpMessage();
};

} // namespace FieldUtils
} // namespace Nektar

#endif
