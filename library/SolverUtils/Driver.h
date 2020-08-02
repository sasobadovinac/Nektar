///////////////////////////////////////////////////////////////////////////////
//
// File Driver.h
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
// Description: Base class for Drivers.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef SOLVERUTILS_DRIVER_H
#define SOLVERUTILS_DRIVER_H

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <SolverUtils/SolverUtils.hpp>
#include <SolverUtils/EquationSystem.h>
#include<SolverUtils/DriverOperators.hpp>

namespace Nektar
{
namespace SolverUtils
{

class Driver;

/// A shared pointer to a Driver object
typedef std::shared_ptr<Driver> DriverSharedPtr;

/// Datatype of the NekFactory used to instantiate classes derived from
/// the Driver class.
typedef LibUtilities::NekFactory<
    std::string, Driver,
    const LibUtilities::SessionReaderSharedPtr &,
    const SpatialDomains::MeshGraphSharedPtr &> DriverFactory;

SOLVER_UTILS_EXPORT DriverFactory& GetDriverFactory();

/// Base class for the development of solvers.
class Driver
{
public:
    /// Destructor
    virtual ~Driver();

    SOLVER_UTILS_EXPORT inline void DoMultiOrderProjection(
        const Array<OneD,const Array<OneD, NekDouble>>  &inarray,
        Array<OneD, Array<OneD,NekDouble>>              &outarray,
        const NekDouble time)
    {
        v_DoMultiOrderProjection(inarray,outarray,time);
    }

    SOLVER_UTILS_EXPORT inline void DoMultiOrderOdeRhs(
        const Array<OneD,const Array<OneD, NekDouble>>  &inarray,
        Array<OneD, Array<OneD,NekDouble>>              &outarray,
        const NekDouble time)
    {
        v_DoMultiOrderOdeRhs(inarray,outarray,time);
    }

    SOLVER_UTILS_EXPORT inline void MultiLevel(
        const TensorOfArray1D<NekDouble>    &inarray,
        TensorOfArray1D<NekDouble>          &outarray,
        const bool                          updateOperatorflag,
        const int                           level,
        const TensorOfArray2D<NekDouble>    &refSolution,
        const NekDouble                     time,
        const NekDouble                     dtlamda)
    {
        v_MultiLevel(inarray,outarray, updateOperatorflag, level, refSolution, 
            time, dtlamda);
    }

    SOLVER_UTILS_EXPORT inline void MultiLvlJacMultiplyMatFree(
        const int                         Level,
        const TensorOfArray1D<NekDouble>  &inarray, 
        TensorOfArray1D<NekDouble>        &out, 
        const NekDouble                   time, 
        const NekDouble                   dtlamda, 
        const TensorOfArray2D<NekDouble>  &refFields, 
        const bool                        flagUpdateJac,
        const bool                        flagDoAdv = true,
        const bool                        flagDoVis = true,
        const bool                        flagSourc = true)
    {
        v_MultiLvlJacMultiplyMatFree(
            Level, inarray, out, time, dtlamda, refFields, flagUpdateJac,
            flagDoAdv, flagDoVis, flagSourc);
    }

    /// Initialise Object
    SOLVER_UTILS_EXPORT inline void InitObject(std::ostream &out = std::cout);

    /// Execute driver
    SOLVER_UTILS_EXPORT inline void Execute(std::ostream &out = std::cout);

    SOLVER_UTILS_EXPORT inline Array<OneD, EquationSystemSharedPtr> GetEqu();

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> GetRealEvl(void);
    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> GetImagEvl(void);

protected:
    Array<OneD,Array<OneD,DNekMatSharedPtr>>    m_RestrictionResidualMatrix;

    Array<OneD,Array<OneD,DNekMatSharedPtr>>    m_RestrictionMatrix;

    Array<OneD,Array<OneD,DNekMatSharedPtr>>    m_ProlongationMatrix;

    Array<OneD,int>                             m_MultiLevelCoeffs;
    /// Communication object
    LibUtilities::CommSharedPtr                 m_comm;

    /// Session reader object
    LibUtilities::SessionReaderSharedPtr        m_session;

    /// I the Coupling between SFD and arnoldi
    LibUtilities::SessionReaderSharedPtr        session_LinNS;

    /// MeshGraph object
    SpatialDomains::MeshGraphSharedPtr          m_graph;

    /// Equation system to solve
    Array<OneD, EquationSystemSharedPtr>        m_equ;
    
    //For use of MultiLevel
    int                                         m_nLevels;

    ///number of equations
    int m_nequ;

    SolverUtils::DriverOperators                m_driverOperator;

    Array<OneD,int>                             m_multilvlncoeffOffset;
    Array<OneD,int>                             m_multilvlnphyscOffset;


    ///Evolution Operator
    enum EvolutionOperatorType m_EvolutionOperator;

    /// Initialises EquationSystem class members.
    Driver(const LibUtilities::SessionReaderSharedPtr pSession,
           const SpatialDomains::MeshGraphSharedPtr   pGraph);
    
    SOLVER_UTILS_EXPORT virtual void v_DoMultiOrderProjection(
        const Array<OneD,const Array<OneD, NekDouble>>  &inarray,
        Array<OneD, Array<OneD,NekDouble>>              &outarray,
        const NekDouble                                 time)
    {
         ASSERTL0(false,"This routine is not valid in this class");
    }

    SOLVER_UTILS_EXPORT virtual void v_DoMultiOrderOdeRhs(
        const Array<OneD,const Array<OneD, NekDouble>>  &inarray,
        Array<OneD, Array<OneD,NekDouble>>              &outarray,
        const NekDouble                                 time)
    {
         ASSERTL0(false,"This routine is not valid in this class");
    }

    SOLVER_UTILS_EXPORT virtual void v_MultiLevel(
        const TensorOfArray1D<NekDouble>    &inarray,
        TensorOfArray1D<NekDouble>          &outarray,
        const bool                          updateOperatorflag,
        const int                           level,
        const TensorOfArray2D<NekDouble>    &refSolution,
        const NekDouble                     time,
        const NekDouble                     dtlamda)
    {
         ASSERTL0(false,"This routine is not valid in this class");
    }

    SOLVER_UTILS_EXPORT virtual void v_MultiLvlJacMultiplyMatFree(
        const int                         Level,
        const TensorOfArray1D<NekDouble>  &inarray, 
        TensorOfArray1D<NekDouble>        &out, 
        const NekDouble                   time, 
        const NekDouble                   dtlamda, 
        const TensorOfArray2D<NekDouble>  &refFields, 
        const bool                        flagUpdateJac,
        const bool                        flagDoAdv,
        const bool                        flagDoVis,
        const bool                        flagSourc)
    {
         ASSERTL0(false,"MultiLvlJacMultiplyMatFree is not valid");
    }

    SOLVER_UTILS_EXPORT virtual void v_InitObject(std::ostream &out = std::cout);

    /// Virtual function for solve implementation.
    SOLVER_UTILS_EXPORT virtual void v_Execute(std::ostream &out = std::cout) = 0;


    SOLVER_UTILS_EXPORT virtual Array<OneD, NekDouble> v_GetRealEvl(void);
    SOLVER_UTILS_EXPORT virtual Array<OneD, NekDouble> v_GetImagEvl(void);


    static std::string evolutionOperatorLookupIds[];
    static std::string evolutionOperatorDef;
    static std::string driverDefault;

};

inline void Driver::InitObject(std::ostream &out)
{
    v_InitObject(out);
}

inline void Driver::Execute(std::ostream &out)
{
    v_Execute(out);
}

inline Array<OneD, EquationSystemSharedPtr> Driver::GetEqu()
{
    return m_equ;
}

inline Array<OneD, NekDouble> Driver::GetRealEvl()
{
    return v_GetRealEvl();
}

inline Array<OneD, NekDouble> Driver::GetImagEvl()
{
    return v_GetImagEvl();
}

}
} //end of namespace

#endif //NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H

