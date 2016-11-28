////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.h
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
//  Description: Calculate jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_PROCESSVAROPTI
#define UTILITIES_NEKMESH_PROCESSVAROPTI

#include <NekMeshUtils/Module/Module.h>

#include "ElUtil.h"

namespace Nektar
{
namespace Utilities
{

struct DerivUtil
{
    NekMatrix<NekDouble> VdmD[3];
    NekMatrix<NekDouble> VdmDStd[3]; //deriv matrix without interp
    NekVector<NekDouble> quadW;

    Array<OneD, Array<OneD, NekDouble> > basisDeriv;

    int pts;
    int ptsStd;
};
typedef boost::shared_ptr<DerivUtil> DerivUtilSharedPtr;

enum optimiser
{
    eLinEl,
    eWins,
    eRoca,
    eHypEl
};

struct Residual
{
    NekDouble val;
    int n;
    int nDoF;
    int startInv;
    int nReset[3];
    NekDouble worstJac;
    NekDouble func;
    NekDouble alphaAvg;
    int alphaI;
};

typedef boost::shared_ptr<Residual> ResidualSharedPtr;

class ProcessVarOpti : public ProcessModule
{
public:
    /// Creates an instance of this class
    static boost::shared_ptr<Module> create(MeshSharedPtr m)
    {
        return MemoryManager<ProcessVarOpti>::AllocateSharedPtr(m);
    }
    static ModuleKey className;

    ProcessVarOpti(MeshSharedPtr m);
    virtual ~ProcessVarOpti();

    virtual void Process();

private:
    typedef std::map<int, std::vector<ElUtilSharedPtr> > NodeElMap;

    void BuildDerivUtil();
    void GetElementMap();
    std::vector<ElementSharedPtr> GetLockedElements(NekDouble thres);
    std::vector<Array<OneD, NekDouble> > MappingIdealToRef(ElementSharedPtr el);
    std::vector<std::vector<NodeSharedPtr> > GetColouredNodes(std::vector<ElementSharedPtr> elLock);

    NodeElMap nodeElMap;
    std::vector<ElUtilSharedPtr> dataSet;

    ResidualSharedPtr res;
    std::map<LibUtilities::ShapeType,DerivUtilSharedPtr> derivUtil;
    optimiser opti;
};

}
}

#endif
