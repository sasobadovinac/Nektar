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

#ifndef NEKTAR_SOLVERUTILS_DRIVEROPERATORS_HPP
#define NEKTAR_SOLVERUTILS_DRIVEROPERATORS_HPP

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    namespace SolverUtils
    {
        class DriverOperators
        {
        public:

            typedef const Array<OneD, NekDouble> InArrayType1;
            typedef       Array<OneD,NekDouble>  OutArrayType1;
            typedef std::function< void (InArrayType1&, OutArrayType1&, const bool, const int)>  FunctorType1;

            typedef const Array<OneD, const Array<OneD, NekDouble>> InArrayType2;
            typedef       Array<OneD, Array<OneD,NekDouble>>       OutArrayType2;
            typedef std::function< void (InArrayType2&, OutArrayType2&, const NekDouble)>  FunctorType2;

            typedef std::function< void (InArrayType2&, const NekDouble, const NekDouble, const int)>  FunctorType3;
            
            DriverOperators(void):
            functors1(1),
            functors2(2),
            functors3(1)
            {
            }

            DriverOperators(DriverOperators &in):
            functors1(1),
            functors2(2),
            functors3(1)
            {
                for (int i = 0; i < 1; i++)
                {
                    functors1[i] = in.functors1[i];
                }

                for (int i = 0; i < 2; i++)
                {
                    functors2[i] = in.functors2[i];
                }


                for (int i = 0; i < 1; i++)
                {
                    functors3[i] = in.functors3[i];
                }
            }

            //Set a functor that m_equ[k] can access m_equ[k+1]'s MultiLevel 
            template<typename FuncPointerT, typename ObjectPointerT> 
            void DefineMultiLevel(FuncPointerT func, ObjectPointerT obj)
            {
                 functors1[0]=  std::bind(func, obj, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
            }

            inline void MultiLevel(InArrayType1     &inarray, 
                                   OutArrayType1    &outarray,
                                   const bool UpDateOperatorflag,
                                   const int  Level) const
            {
                ASSERTL1(functors1[0],"MultiLevel should be defined in InitObject");
                functors1[0](inarray, outarray, UpDateOperatorflag, Level);
            }

            //Set a functor that m_equ[0] can do Projection
            template<typename FuncPointerT, typename ObjectPointerT> 
            void DefineMultiOrderProjection(FuncPointerT func, ObjectPointerT obj)
            {
                 functors2[0]=  std::bind(func, obj, std::placeholders::_1, std::placeholders::_2,std::placeholders::_3);
            }

            inline void DoMultiOrderProjection(InArrayType2     &inarray, 
                                               OutArrayType2    &outarray,
                                               const NekDouble time) const
            {
                ASSERTL1(functors2[0],"DoMultiOrderProjection should be defined for this time integration scheme");
                functors2[0](inarray,outarray,time);
            }

            //Set a functor that m_equ[0] can extract Rhs from m_equ[1]
            template<typename FuncPointerT, typename ObjectPointerT> 
            void DefineMultiOrderOdeRhs(FuncPointerT func, ObjectPointerT obj)
            {
                 functors2[1]=  std::bind(func, obj, std::placeholders::_1, std::placeholders::_2,std::placeholders::_3);
            }

            inline void DoMultiOrderOdeRhs(InArrayType2     &inarray, 
                                           OutArrayType2    &outarray,
                                               const NekDouble time) const
            {
                ASSERTL1(functors2[1],"DoMultiOrderOdeRhs should be defined for this time integration scheme");
                functors2[1](inarray,outarray,time);
            }

            //Set a functor that calculate Next level's Stored matrices
            template<typename FuncPointerT, typename ObjectPointerT> 
            void DefineCalculateNextLevelPreconditioner(FuncPointerT func, ObjectPointerT obj)
            {
                 functors3[0]=  std::bind(func, obj, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
            }

            inline void  CalculateNextLevelPreconditioner(InArrayType2     &inarray, 
                                                          const NekDouble time,
                                                          const NekDouble lambda,
                                                          const int       NextLevel) const
            {
                ASSERTL1(functors3[0],"CalculateNextLevelPreconditioner has not been defined");
                functors3[0](inarray,time,lambda,NextLevel);
            }

        protected:
            Array<OneD,FunctorType1> functors1;
            Array<OneD,FunctorType2> functors2;
            Array<OneD,FunctorType3> functors3;
        private:

        };
    }	
} //end of namespace

#endif //
