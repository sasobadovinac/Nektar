///////////////////////////////////////////////////////////////////////////////
//
// File: Driver.h
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
// Description: Base class for Drivers.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef SOLVERUTILS_DRIVER_H
#define SOLVERUTILS_DRIVER_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/Comm.h>

#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/SolverUtils.hpp>

namespace Nektar
{
namespace SolverUtils
{

class Driver;

/// A shared pointer to a Driver object
typedef std::shared_ptr<Driver> DriverSharedPtr;

/// Datatype of the NekFactory used to instantiate classes derived from
/// the Driver class.
typedef LibUtilities::NekFactory<std::string, Driver,
                                 const LibUtilities::SessionReaderSharedPtr &,
                                 const SpatialDomains::MeshGraphSharedPtr &>
    DriverFactory;

SOLVER_UTILS_EXPORT DriverFactory &GetDriverFactory();

/// Base class for the development of solvers.
class Driver
{
public:
    /// Destructor
    virtual ~Driver();

    /// Initialise Object
    SOLVER_UTILS_EXPORT inline void InitObject(std::ostream &out = std::cout);

    /// Execute driver
    SOLVER_UTILS_EXPORT inline void Execute(std::ostream &out = std::cout);

    SOLVER_UTILS_EXPORT inline Array<OneD, EquationSystemSharedPtr> GetEqu();

protected:
    /// Communication object
    LibUtilities::CommSharedPtr m_comm;

    /// Session reader object
    LibUtilities::SessionReaderSharedPtr m_session;

    /// Coupling between SFD and arnoldi
    LibUtilities::SessionReaderSharedPtr session_LinNS;

    /// MeshGraph object
    SpatialDomains::MeshGraphSharedPtr m_graph;

    /// Equation system to solve
    Array<OneD, EquationSystemSharedPtr> m_equ;

    /// number of equations
    int m_nequ;

    /// Evolution Operator
    enum EvolutionOperatorType m_EvolutionOperator;

    /// Initialises EquationSystem class members.
    Driver(const LibUtilities::SessionReaderSharedPtr pSession,
           const SpatialDomains::MeshGraphSharedPtr pGraph);

    /// Virtual function for initialisation implementation.
    SOLVER_UTILS_EXPORT virtual void v_InitObject(
        std::ostream &out = std::cout);

    /// Virtual function for solve implementation.
    SOLVER_UTILS_EXPORT virtual void v_Execute(
        std::ostream &out = std::cout) = 0;

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

} // namespace SolverUtils
} // namespace Nektar

#endif // NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H
