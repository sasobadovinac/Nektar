///////////////////////////////////////////////////////////////////////////////
//
// File: PreconCfsOp.h
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
// Description: PreconCfsOp header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_PRECONCFSOP
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_PRECONCFSOP

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>

namespace Nektar
{
// =====================================================================
// ==== Defines operators needed by operator based preconditioner in the
// ==== CFS solver
// =====================================================================
class NekPreconCfsOperators
{
public:
    typedef const Array<OneD, NekDouble> InArrayType;
    typedef Array<OneD, NekDouble> OutArrayType;

    typedef std::function<void(
        const Array<OneD, const Array<OneD, NekDouble>> &,
        Array<OneD, Array<OneD, SNekBlkMatSharedPtr>> &, SNekBlkMatSharedPtr &,
        Array<OneD, SNekBlkMatSharedPtr> &, Array<OneD, SNekBlkMatSharedPtr> &,
        Array<OneD, Array<OneD, NekSingle>> &, TensorOfArray4D<NekSingle> &,
        TensorOfArray4D<NekSingle> &, TensorOfArray5D<NekSingle> &)>
        Functor;

    typedef Array<OneD, Functor> FunctorArray;
    static const int nfunctor = 1;

    NekPreconCfsOperators(void) : m_functors(nfunctor)
    {
    }

    NekPreconCfsOperators(const NekPreconCfsOperators &in)
        : m_functors(nfunctor)
    {
        for (int i = 0; i < nfunctor; ++i)
        {
            m_functors[i] = in.m_functors[i];
        }
    }

    NekPreconCfsOperators &operator=(const NekPreconCfsOperators &in)
    {
        for (int i = 0; i < nfunctor; ++i)
        {
            m_functors[i] = in.m_functors[i];
        }
        return *this;
    }

    template <typename FuncPointerT, typename ObjectPointerT>
    void DefineCalcPreconMatBRJCoeff(FuncPointerT func, ObjectPointerT obj)
    {
        m_functors[0] = std::bind(func, obj, std::placeholders::_1,
                                  std::placeholders::_2, std::placeholders::_3,
                                  std::placeholders::_4, std::placeholders::_5,
                                  std::placeholders::_6, std::placeholders::_7,
                                  std::placeholders::_8, std::placeholders::_9);
    }

    inline void DoCalcPreconMatBRJCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, SNekBlkMatSharedPtr>> &gmtxarray,
        SNekBlkMatSharedPtr &gmtVar, Array<OneD, SNekBlkMatSharedPtr> &TraceJac,
        Array<OneD, SNekBlkMatSharedPtr> &TraceJacDeriv,
        Array<OneD, Array<OneD, NekSingle>> &TraceJacDerivSign,
        TensorOfArray4D<NekSingle> &TraceJacArray,
        TensorOfArray4D<NekSingle> &TraceJacDerivArray,
        TensorOfArray5D<NekSingle> &TraceIPSymJacArray)
    {
        ASSERTL1(m_functors[0], "DoNekSysResEval should be defined");
        m_functors[0](inarray, gmtxarray, gmtVar, TraceJac, TraceJacDeriv,
                      TraceJacDerivSign, TraceJacArray, TraceJacDerivArray,
                      TraceIPSymJacArray);
    }

protected:
    FunctorArray m_functors;
};

} // namespace Nektar

#endif
