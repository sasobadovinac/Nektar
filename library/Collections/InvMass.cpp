///////////////////////////////////////////////////////////////////////////////
//
// File: BwdTrans.cpp
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
// Description: BwdTrans operator implementations
//
///////////////////////////////////////////////////////////////////////////////

#include <Collections/CoalescedGeomData.h>
#include <Collections/Operator.h>
#include <Collections/MatrixFreeBase.h>
#include <MatrixFreeOps/Operator.hpp>
#include <boost/core/ignore_unused.hpp>

using namespace std;

namespace Nektar
{
namespace Collections
{

using LibUtilities::eHexahedron;
using LibUtilities::ePrism;
using LibUtilities::ePyramid;
using LibUtilities::eQuadrilateral;
using LibUtilities::eSegment;
using LibUtilities::eTetrahedron;
using LibUtilities::eTriangle;

/**
 * @brief Backward transform operator using matrix free operators.
 */
class InvMass_MatrixFree final : public Operator, MatrixFreeOneInOneOut
{
public:
    OPERATOR_CREATE(InvMass_MatrixFree)

    ~InvMass_MatrixFree() final
    {
    }

    void operator()(const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output0,
                    Array<OneD, NekDouble> &output1,
                    Array<OneD, NekDouble> &output2,
                    Array<OneD, NekDouble> &wsp) final
    {
        boost::ignore_unused(output1, output2, wsp);

        if (m_isPadded)
        {
            // copy into padded vector
            Vmath::Vcopy(m_nIn, input, 1, m_input, 1);
            // call op
            (*m_oper)(m_input, m_output);
            // copy out of padded vector
            Vmath::Vcopy(m_nOut, m_output, 1, output0, 1);
        }
        else
        {
            (*m_oper)(input, output0);
        }
    }

    void operator()(int dir, const Array<OneD, const NekDouble> &input,
                    Array<OneD, NekDouble> &output,
                    Array<OneD, NekDouble> &wsp) final
    {
        boost::ignore_unused(dir, input, output, wsp);
        NEKERROR(ErrorUtil::efatal,
                 "InvMass_MatrixFree: Not valid for this operator.");
    }
    
    virtual void CheckFactors(StdRegions::FactorMap factors,
                              int coll_phys_offset)
    {
        boost::ignore_unused(factors, coll_phys_offset);
        ASSERTL0(false, "Not valid for this operator.");
    }


private:
    std::shared_ptr<MatrixFree::InvMass> m_oper;

    InvMass_MatrixFree(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                       CoalescedGeomDataSharedPtr pGeomData,
                       StdRegions::FactorMap factors)
        : Operator(pCollExp, pGeomData, factors),
          MatrixFreeOneInOneOut(pCollExp[0]->GetStdExp()->GetNcoeffs(),
                                pCollExp[0]->GetStdExp()->GetTotPoints(),
                                pCollExp.size())
    {
        // Basis vector.
        const auto dim = pCollExp[0]->GetStdExp()->GetShapeDimension();
        std::vector<LibUtilities::BasisSharedPtr> basis(dim);
        for (auto i = 0; i < dim; ++i)
        {
            basis[i] = pCollExp[0]->GetBasis(i);
        }

        // Get shape type
        auto shapeType = pCollExp[0]->GetStdExp()->DetShapeType();

        // Generate operator string and create operator.
        std::string op_string = "InvMass";
        op_string += MatrixFree::GetOpstring(shapeType, false);
        auto oper = MatrixFree::GetOperatorFactory().CreateInstance(
            op_string, basis, m_nElmtPad);

        m_oper = std::dynamic_pointer_cast<MatrixFree::InvMass>(oper);

        // Set Jacobian
        oper->SetJac(pGeomData->GetJacInterLeave(pCollExp,m_nElmtPad));

        ASSERTL0(m_oper, "Failed to cast pointer.");
    }
};

OperatorKey InvMass_MatrixFree::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eInvMass, eMatrixFree, false),
        InvMass_MatrixFree::create, "InvMass_MatrixFree_Quad"),
};


} // namespace Collections

} // namespace Nektar
