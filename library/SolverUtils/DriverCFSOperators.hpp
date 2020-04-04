///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionTerm.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Driver class for the stability solver
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DRIVERCFSOPERATORS_HPP
#define NEKTAR_SOLVERUTILS_DRIVERCFSOPERATORS_HPP

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <SolverUtils/DriverCFSOperators.hpp>

namespace Nektar
{
    namespace SolverUtils
    {
        class DriverOperators
        {
        public:

            typedef const Array<OneD, const Array<OneD, NekDouble>> InArrayType;
            typedef       Array<OneD, Array<OneD,NekDouble>>       OutArrayType;
            typedef std::function< void (InArrayType&, OutArrayType&, const NekDouble)>  FunctorType;
            
            DriverOperators(void):
            m_functors(2)
            {
            }

            DriverOperators(DriverOperators &in):
            m_functors(2)
            {
                for (int i = 0; i < 2; i++)
                {
                    m_functors[i] = in.m_functors[i];
                }
            }

            //Set a functor that m_equ[0] can do Projection
            template<typename FuncPointerT, typename ObjectPointerT> 
            void DefineMultiOrderProjection(FuncPointerT func, ObjectPointerT obj)
            {
                 m_functors[0]=  std::bind(func, obj, std::placeholders::_1, std::placeholders::_2,std::placeholders::_3);
            }

            inline void DoMultiOrderProjection(InArrayType     &inarray, 
                                           OutArrayType    &outarray,
                                               const NekDouble time) const
            {
                ASSERTL1(m_functors[0],"DoHigherOrderOdeRhs should be defined for this time integration scheme");
                m_functors[0](inarray,outarray,time);
            }

            //Set a functor that m_equ[0] can extract Rhs from m_equ[1]
            template<typename FuncPointerT, typename ObjectPointerT> 
            void DefineMultiOrderOdeRhs(FuncPointerT func, ObjectPointerT obj)
            {
                 m_functors[1]=  std::bind(func, obj, std::placeholders::_1, std::placeholders::_2,std::placeholders::_3);
            }

            inline void DoMultiOrderOdeRhs(InArrayType     &inarray, 
                                           OutArrayType    &outarray,
                                               const NekDouble time) const
            {
                ASSERTL1(m_functors[1],"DoHigherOrderOdeRhs should be defined for this time integration scheme");
                m_functors[1](inarray,outarray,time);
            }

        protected:
            Array<OneD,FunctorType> m_functors;
        private:

        };
    }	
} //end of namespace

#endif //
