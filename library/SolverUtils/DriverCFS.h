///////////////////////////////////////////////////////////////////////////////
//
// File DriverCFS.h
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

#ifndef SOLVERUTILS_DRIVERCFS_H
#define SOLVERUTILS_DRIVERCFS_H

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <SolverUtils/CFSSolverUtils.hpp>
#include <SolverUtils/EquationSystem.h>
#include<SolverUtils/DriverCFSOperators.hpp>

namespace Nektar
{
namespace SolverUtils
{

class DriverCFS;

/// A shared pointer to a DriverCFS object
typedef std::shared_ptr<DriverCFS> DriverCFSSharedPtr;

/// Datatype of the NekFactory used to instantiate classes derived from
/// the DriverCFS class.
typedef LibUtilities::NekFactory<
    std::string, DriverCFS,
    const LibUtilities::SessionReaderSharedPtr &,
    const SpatialDomains::MeshGraphSharedPtr &,
    const SpatialDomains::MeshGraphSharedPtr &> DriverCFSFactory;

SOLVER_UTILS_EXPORT DriverCFSFactory& GetDriverCFSFactory();

/// Base class for the development of solvers.
class DriverCFS
{
public:
    /// Destructor
    virtual ~DriverCFS();

    SOLVER_UTILS_EXPORT inline void DoMultiOrderOdeRhs(const Array<OneD,const Array<OneD, NekDouble>> &inarray,
                                                                   Array<OneD, Array<OneD,NekDouble>> &outarray,
                                                                                           const NekDouble time)
    {
        v_DoMultiOrderOdeRhs(inarray,outarray,time);
    }
        

    /// Initialise Object
    SOLVER_UTILS_EXPORT inline void InitObject(std::ostream &out = std::cout);

    /// Execute driver
    SOLVER_UTILS_EXPORT inline void Execute(std::ostream &out = std::cout);

    SOLVER_UTILS_EXPORT inline Array<OneD, EquationSystemSharedPtr> GetEqu();

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> GetRealEvl(void);
    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> GetImagEvl(void);


protected:
    /// Communication object
    LibUtilities::CommSharedPtr                 m_comm;

    /// Session reader object
    LibUtilities::SessionReaderSharedPtr        m_session;

    /// I the Coupling between SFD and arnoldi
    LibUtilities::SessionReaderSharedPtr        session_LinNS;

    /// MeshGraph object
    SpatialDomains::MeshGraphSharedPtr          m_graph;
    SpatialDomains::MeshGraphSharedPtr          m_HigherOrdergraph;

    /// Equation system to solve
    Array<OneD, EquationSystemSharedPtr>        m_equ;

    ///number of equations
    int m_nequ;
    
    SolverUtils::DriverOperators                    m_driverOperator;

    ///Evolution Operator
    enum EvolutionOperatorType m_EvolutionOperator;

    /// Initialises EquationSystem class members.
    DriverCFS(const LibUtilities::SessionReaderSharedPtr pSession,
           const SpatialDomains::MeshGraphSharedPtr   pGraph,
           const SpatialDomains::MeshGraphSharedPtr   pHigherOrderGraph);

    SOLVER_UTILS_EXPORT virtual void v_DoMultiOrderOdeRhs(const Array<OneD,const Array<OneD, NekDouble>> &inarray,
                                                                   Array<OneD, Array<OneD,NekDouble>> &outarray,
                                                                              const NekDouble time)
    {
         ASSERTL0(false,"This routine is not valid in this class");
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

inline void DriverCFS::InitObject(std::ostream &out)
{
    v_InitObject(out);
}

inline void DriverCFS::Execute(std::ostream &out)
{
    v_Execute(out);
}

inline Array<OneD, EquationSystemSharedPtr> DriverCFS::GetEqu()
{
    return m_equ;
}

inline Array<OneD, NekDouble> DriverCFS::GetRealEvl()
{
    return v_GetRealEvl();
}

inline Array<OneD, NekDouble> DriverCFS::GetImagEvl()
{
    return v_GetImagEvl();
}

}
} //end of namespace

#endif //NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H

